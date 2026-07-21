# internal variable rescaling for badly scaled data (new in 0.7-2)
#
# When the observed variables have wildly different scales (variance
# ratios of 1e6 and beyond), the quasi-Newton optimizer often fails or
# needs thousands of iterations: several starting values are absolute
# (raw-metric) constants, the parameter vector mixes magnitudes the
# step-control heuristics cannot handle, and for the (D)WLS/GLS family
# the weight matrices become numerically ill-conditioned. The textbook
# advice is to rescale the offending variables before the analysis.
#
# This module automates that advice EXACTLY: each model variable is
# divided by a power of two close to its standard deviation (powers of
# two are exactly representable, so the transformation is lossless in
# floating point), the model is fitted in the well-scaled metric by a
# recursive lavaan() call, and the solution is mapped back to the
# original metric through the per-parameter monomials of the
# equivariance map (loadings scale by c_ov / c_lv, variances by
# c_lhs * c_rhs, regressions by c_lhs / c_rhs, intercepts by c_lhs,
# where a latent variable inherits the scale of its marker indicator).
# Fixed nonzero parameter values and user starting values are
# transformed on the way in and back on the way out, so the fitted
# model is exactly the user's model -- the rescaling is a pure change
# of units, invisible in the output. Everything downstream of the
# optimizer (implied moments, standard errors, test statistics,
# fit measures) is computed by the ordinary pipeline in the original
# metric, at the returned solution.
#
# The map is exact whenever every model 'constant' can be transformed
# along: this fails only for equality/inequality constraints that tie
# parameters with different monomials (e.g., equal loadings across
# indicators with different scales), which lav_rescale_prep() detects
# from the parameter table (the rescaling is then skipped).

# decide whether to rescale, and prepare the scale factors + monomials;
# returns NULL (do not rescale) or list(cj = , mono = )
lav_rescale_prep <- function(lavoptions = NULL, lavdata = NULL,
                             lavmodel = NULL, lavpartable = NULL) {
  opt <- lavoptions$rescale.data
  if (is.null(opt) || isFALSE(opt) || identical(opt, "none")) {
    return(NULL)
  }
  forced <- isTRUE(opt)

  # a helper to bail out (with a message if the user forced rescaling)
  skip <- function(reason) {
    if (forced) {
      lav_msg_warn(gettextf(
        "rescale.data = TRUE is not available for this model (%s);
        continuing without rescaling.", reason))
    }
    NULL
  }

  # hard requirements
  if (!isTRUE(lavoptions$do.fit) || lavmodel@nx.free == 0L) {
    return(NULL)
  }
  if (lavdata@data.type != "full") {
    return(skip(gettext("raw data required")))
  }
  if (lavdata@nlevels > 1L) {
    return(skip(gettext("multilevel data")))
  }
  if (lavmodel@categorical) {
    return(skip(gettext("categorical data")))
  }
  if (lavmodel@representation != "LISREL") {
    return(skip(gettext("representation is not LISREL")))
  }
  if (lavmodel@correlation || isTRUE(lavoptions$.correlation.ml)) {
    # correlation structures have their own scale handling
    return(NULL)
  }
  if (lavmodel@composites || any(lavpartable$op == "<~")) {
    return(skip(gettext("composites")))
  }
  if (!is.null(lavpartable$efa) && any(nchar(lavpartable$efa) > 0L)) {
    return(skip(gettext("efa blocks")))
  }
  if (any(grepl(":", c(lavpartable$lhs, lavpartable$rhs), fixed = TRUE))) {
    return(skip(gettext("product (interaction) terms")))
  }
  est_group <- lav_options_estimatorgroup(lavoptions$estimator)
  if (!est_group %in% c("ML", "GLS", "WLS", "DWLS", "NTRLS")) {
    # ULS is not scale-equivariant; noniterative estimators do not
    # need this
    return(skip(gettextf("estimator %s", est_group)))
  }
  if (!lavoptions$optim.method %in% c("nlminb", "nlminb0", "nlminb1",
                                      "gn")) {
    return(skip(gettextf("optim.method %s", lavoptions$optim.method)))
  }
  if (!lavoptions$missing %in% c("listwise", "ml", "ml.x")) {
    return(skip(gettextf("missing = %s", lavoptions$missing)))
  }
  if (any(lavpartable$op %in% c("<", ">"))) {
    # inequality constraints: the active-set information (con.lambda)
    # would live in the wrong metric
    return(skip(gettext("inequality constraints")))
  }

  # collect the model variables and their (pooled) standard deviations
  ov_names <- unique(unlist(c(
    lav_pt_vnames(lavpartable, "ov"),
    lav_pt_vnames(lavpartable, "ov.x")
  )))
  sds <- setNames(rep(NA_real_, length(ov_names)), ov_names)
  for (v in ov_names) {
    vals <- numeric(0L)
    for (g in seq_len(lavdata@ngroups)) {
      pos <- match(v, lavdata@ov.names[[g]])
      if (!is.na(pos) && length(lavdata@X) >= g &&
          !is.null(lavdata@X[[g]])) {
        vals <- c(vals, lavdata@X[[g]][, pos])
      } else {
        pos <- match(v, lavdata@ov.names.x[[g]])
        if (!is.na(pos) && length(lavdata@eXo) >= g &&
            !is.null(lavdata@eXo[[g]])) {
          vals <- c(vals, lavdata@eXo[[g]][, pos])
        }
      }
    }
    if (length(vals) > 1L) {
      sds[v] <- sd(vals, na.rm = TRUE)
    }
  }
  if (anyNA(sds) || any(!is.finite(sds)) || any(sds <= 0)) {
    return(skip(gettext("could not compute standard deviations")))
  }

  # trigger (auto): only for clearly badly scaled data
  if (!forced) {
    if (max(sds) / min(sds) < 1000 &&
        max(sds) < 1000 && min(sds) > 0.001) {
      return(NULL)
    }
  }

  # power-of-two scale factors (lossless in floating point)
  cj <- 2^round(log2(sds))
  if (all(cj == 1)) {
    return(NULL)
  }

  # per-row monomials of the equivariance map
  mono <- lav_rescale_mono(lavpartable = lavpartable, cj = cj)
  if (is.null(mono)) {
    return(skip(gettext(
      "constraints tie parameters with different scaling factors")))
  }

  list(cj = cj, mono = mono)
}

# per-row monomials: est_original = est_scaled * mono; returns NULL if
# the model has constraints that do not survive the rescaling
lav_rescale_mono <- function(lavpartable = NULL, cj = NULL) {
  n <- length(lavpartable$lhs)
  mono <- rep(1.0, n)
  group <- lavpartable$group
  if (is.null(group)) {
    group <- rep(1L, n)
  }
  group_values <- unique(group[group > 0L])

  for (g in group_values) {
    g_rows <- which(group == g)
    lv_names <- unique(lavpartable$lhs[g_rows][
      lavpartable$op[g_rows] == "=~"])
    # resolve the latent scale factors: a latent variable inherits the
    # scale of its marker (the first fixed nonzero loading); latent
    # variables without a marker (std.lv, or variance-anchored) keep
    # scale 1 -- any choice is valid, the map transforms all fixed
    # values consistently
    s_lv <- setNames(rep(NA_real_, length(lv_names)), lv_names)
    scale_of <- function(v) {
      if (v %in% names(cj)) {
        cj[[v]]
      } else if (v %in% names(s_lv)) {
        s_lv[[v]]
      } else {
        1.0
      }
    }
    if (length(lv_names) > 0L) {
      for (rep_i in 1:5) {
        for (l in lv_names) {
          if (!is.na(s_lv[[l]])) next
          m_idx <- which(group == g & lavpartable$op == "=~" &
            lavpartable$lhs == l & lavpartable$free == 0L &
            !is.na(lavpartable$ustart) & lavpartable$ustart != 0)
          if (length(m_idx) == 0L) {
            s_lv[[l]] <- 1.0
            next
          }
          mk <- lavpartable$rhs[m_idx[1]]
          s_lv[[l]] <- scale_of(mk) # NA if mk is a still-unresolved lv
        }
        if (!anyNA(s_lv)) break
      }
      s_lv[is.na(s_lv)] <- 1.0
    }

    for (i in g_rows) {
      op <- lavpartable$op[i]
      if (op == "=~") {
        mono[i] <- scale_of(lavpartable$rhs[i]) /
          scale_of(lavpartable$lhs[i])
      } else if (op == "~") {
        mono[i] <- scale_of(lavpartable$lhs[i]) /
          scale_of(lavpartable$rhs[i])
      } else if (op == "~~") {
        mono[i] <- scale_of(lavpartable$lhs[i]) *
          scale_of(lavpartable$rhs[i])
      } else if (op == "~1") {
        mono[i] <- scale_of(lavpartable$lhs[i])
      }
      # all other ops (==, :=, ...) keep mono 1
    }
  }

  # constraint compatibility checks -----------------------------------
  # (a) shared free indices (simple equality constraints, e.g.
  #     group.equal): all rows sharing a free index must share the
  #     same monomial
  free <- lavpartable$free
  dup_ids <- unique(free[free > 0L][duplicated(free[free > 0L])])
  for (id in dup_ids) {
    m_id <- mono[free == id]
    if (max(abs(m_id / m_id[1] - 1)) > 1e-12) {
      return(NULL)
    }
  }
  # (b) explicit equality constraints: only plain label == label (or
  #     plabel) pairs between parameters with identical monomials
  con_rows <- which(lavpartable$op == "==")
  if (length(con_rows) > 0L) {
    lab_mono <- function(lab) {
      # numeric constant?
      val <- suppressWarnings(as.numeric(lab))
      if (!is.na(val)) {
        # 'label == constant' pins a raw-metric value; only safe if
        # the labeled side has monomial 1 -- handled by the caller
        # returning the constant marker
        return("constant")
      }
      if (!grepl("^[A-Za-z_.][A-Za-z0-9_.]*$", lab)) {
        return(NULL) # an expression: not supported
      }
      idx <- which((!is.na(lavpartable$label) &
                    lavpartable$label == lab) |
                   (!is.null(lavpartable$plabel) &
                    !is.na(lavpartable$plabel) &
                    lavpartable$plabel == lab))
      # exclude the constraint/def rows themselves
      idx <- idx[!lavpartable$op[idx] %in% c("==", "<", ">", ":=")]
      if (length(idx) == 0L) {
        return(NULL) # references a := definition (or unknown): skip
      }
      m_l <- mono[idx]
      if (max(abs(m_l / m_l[1] - 1)) > 1e-12) {
        return(NULL)
      }
      m_l[1]
    }
    for (i in con_rows) {
      m_lhs <- lab_mono(lavpartable$lhs[i])
      m_rhs <- lab_mono(lavpartable$rhs[i])
      if (is.null(m_lhs) || is.null(m_rhs)) {
        return(NULL)
      }
      if (identical(m_lhs, "constant") && identical(m_rhs, "constant")) {
        next
      }
      if (identical(m_lhs, "constant") || identical(m_rhs, "constant")) {
        # label == constant: the raw-metric value is pinned; exact only
        # if the label's monomial is 1
        m_par <- if (identical(m_lhs, "constant")) m_rhs else m_lhs
        if (abs(m_par - 1) > 1e-12) {
          return(NULL)
        }
      } else if (abs(m_lhs / m_rhs - 1) > 1e-12) {
        return(NULL)
      }
    }
  }

  mono
}

# clone lavdata with the model variables divided by cj
lav_rescale_data <- function(lavdata = NULL, cj = NULL) {
  lavdata2 <- lavdata
  for (g in seq_len(lavdata@ngroups)) {
    if (length(lavdata2@X) >= g && !is.null(lavdata2@X[[g]]) &&
        ncol(lavdata2@X[[g]]) > 0L) {
      pos <- match(lavdata@ov.names[[g]], names(cj))
      hit <- which(!is.na(pos))
      if (length(hit) > 0L) {
        lavdata2@X[[g]][, hit] <-
          sweep(lavdata2@X[[g]][, hit, drop = FALSE], 2L,
                cj[pos[hit]], "/")
      }
    }
    if (length(lavdata2@eXo) >= g && !is.null(lavdata2@eXo[[g]]) &&
        ncol(lavdata2@eXo[[g]]) > 0L) {
      pos <- match(lavdata@ov.names.x[[g]], names(cj))
      hit <- which(!is.na(pos))
      if (length(hit) > 0L) {
        lavdata2@eXo[[g]][, hit] <-
          sweep(lavdata2@eXo[[g]][, hit, drop = FALSE], 2L,
                cj[pos[hit]], "/")
      }
    }
  }
  lavdata2
}

# the scaled-metric input parameter table: transform the fixed values
# and user starting values, drop the computed columns
lav_rescale_pt <- function(lavpartable = NULL, mono = NULL) {
  pt <- as.data.frame(lavpartable, stringsAsFactors = FALSE)
  if (!is.null(pt$ustart)) {
    idx <- which(!is.na(pt$ustart) & mono != 1)
    if (length(idx) > 0L) {
      pt$ustart[idx] <- pt$ustart[idx] / mono[idx]
    }
  }
  # drop everything that will be recomputed in the scaled metric
  pt$est <- NULL
  pt$se <- NULL
  pt$start <- NULL
  pt$lower <- NULL
  pt$upper <- NULL
  pt
}

# run the rescaled inner fit; returns the back-transformed 'x' vector
# (with the attributes lav_model_est() would set), or NULL on failure
lav_rescale_estimate <- function(lavdata = NULL, lavmodel = NULL,
                                 lavcache = NULL, lavsamplestats = NULL,
                                 lavoptions = NULL, lavpartable = NULL,
                                 prep = NULL) {
  cj <- prep$cj
  mono <- prep$mono

  lavdata2 <- lav_rescale_data(lavdata = lavdata, cj = cj)
  pt2 <- lav_rescale_pt(lavpartable = lavpartable, mono = mono)

  lavoptions2 <- lavoptions
  lavoptions2$rescale.data <- FALSE # recursion guard
  lavoptions2$baseline <- FALSE
  lavoptions2$.rescale.cj <- NULL
  lavoptions2$.rescale.vcov <- NULL
  lavoptions2$.rescale.test <- NULL
  # for a "plain" ML fit (standard se and test), the outer pipeline can
  # compute the se/test in the original metric (the information matrix
  # inversion is diagonally preconditioned); the inner fit then skips
  # them. For everything else (robust flavours, the (D)WLS/GLS family,
  # bootstrap) the se/test are taken FROM the inner fit: the vcov maps
  # back by the exact monomial congruence, and the test statistics are
  # invariant -- computing them in the well-scaled metric avoids the
  # raw-metric numerics of the sandwich/weight-matrix machinery
  # altogether. (Note: for the (D)WLS/GLS estimators the se/test
  # settings also co-determine which NACOV/weight-matrix quantities the
  # sample statistics compute, so the inner options must keep them to
  # reproduce the outer OBJECTIVE exactly.)
  plain_ml <- lav_options_estimatorgroup(lavoptions$estimator) == "ML" &&
    lavoptions$se %in% c("none", "standard") &&
    all(lavoptions$test %in% c("none", "standard")) &&
    lavoptions$missing == "listwise"
  if (plain_ml) {
    lavoptions2$se <- "none"
    lavoptions2$test <- "none"
  } else {
    lavoptions2$store.vcov <- TRUE
  }

  fit_in <- try(
    lavaan(
      model = pt2,
      slot_data = lavdata2,
      slot_options = lavoptions2
    ),
    silent = TRUE
  )
  if (inherits(fit_in, "try-error")) {
    return(NULL)
  }
  if (!fit_in@optim$converged) {
    return(NULL)
  }

  # back-transform the free parameters: one row per free index
  free <- lavpartable$free
  free_rows <- which(free > 0L)
  free_rows <- free_rows[!duplicated(free[free_rows])]
  free_rows <- free_rows[order(free[free_rows])]
  x_in <- fit_in@optim$x
  if (length(x_in) != length(free_rows)) {
    return(NULL)
  }
  x <- as.numeric(x_in) * mono[free_rows]

  # the objective value in the ORIGINAL metric (for ML this equals the
  # scaled value exactly; for the *LS family the weight matrices differ
  # by the corresponding congruence, so the value is also invariant)
  lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
  fx <- try(
    lav_model_objective(
      lavmodel = lavmodel2,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache
    ),
    silent = TRUE
  )
  if (inherits(fx, "try-error") || !is.finite(fx[1])) {
    return(NULL)
  }

  attr(x, "iterations") <- fit_in@optim$iterations
  attr(x, "converged") <- TRUE
  attr(x, "warn.txt") <- ""
  attr(x, "control") <- lavoptions$control
  attr(x, "dx") <- numeric(0L)
  attr(x, "fx") <- fx
  attr(x, "rescale.cj") <- cj

  # non-plain-ML: hand the (exactly back-transformed) vcov and the
  # (invariant) test statistics of the inner fit to the outer pipeline
  if (!plain_ml) {
    mono_x <- mono[free_rows]
    vcov_in <- fit_in@vcov$vcov
    if (!is.null(vcov_in) && lavoptions$se != "none") {
      vcov_orig <- vcov_in * tcrossprod(mono_x)
      boot_in <- try(fit_in@boot$coef, silent = TRUE)
      if (!inherits(boot_in, "try-error") && !is.null(boot_in) &&
          is.matrix(boot_in) && ncol(boot_in) == length(mono_x)) {
        attr(vcov_orig, "BOOT.COEF") <-
          sweep(boot_in, 2L, mono_x, "*")
      }
      attr(x, "rescale.vcov") <- vcov_orig
    }
    if (length(fit_in@test) > 0L &&
        !(length(lavoptions$test) == 1L && lavoptions$test == "none")) {
      attr(x, "rescale.test") <- fit_in@test
    }

    # missing data: the outer (raw-metric) saturated EM is unreliable
    # for badly scaled data; hand over the inner h1 as well, with the
    # moments mapped back and the log-likelihood corrected by the
    # (pattern-aware) Jacobian term of the rescaling
    if (lavoptions$missing != "listwise" && length(fit_in@h1) > 0L) {
      h1_in <- fit_in@h1
      log_c <- log(cj)
      for (g in seq_len(lavdata@ngroups)) {
        cvec <- cj[lavdata@ov.names[[g]]]
        if (!is.null(h1_in$implied$cov[[g]])) {
          h1_in$implied$cov[[g]] <-
            h1_in$implied$cov[[g]] * tcrossprod(cvec)
        }
        if (!is.null(h1_in$implied$mean[[g]])) {
          h1_in$implied$mean[[g]] <- h1_in$implied$mean[[g]] * cvec
        }
        # Jacobian: sum over cases of sum over OBSERVED variables of
        # log(c_j), via the missing-pattern table
        adj_g <- 0
        mp <- lavdata@Mp[[g]]
        if (!is.null(mp) && length(mp$pat) > 0L) {
          lc_g <- log_c[lavdata@ov.names[[g]]]
          for (pp in seq_len(nrow(mp$pat))) {
            adj_g <- adj_g + mp$freq[pp] * sum(lc_g[mp$pat[pp, ]])
          }
        } else {
          adj_g <- lavdata@nobs[[g]] * sum(log_c[lavdata@ov.names[[g]]])
        }
        if (length(h1_in$logl$loglik.group) >= g) {
          h1_in$logl$loglik.group[g] <-
            h1_in$logl$loglik.group[g] - adj_g
        }
      }
      h1_in$logl$loglik <- sum(h1_in$logl$loglik.group)
      attr(x, "rescale.h1") <- h1_in
    }
  }
  x
}
