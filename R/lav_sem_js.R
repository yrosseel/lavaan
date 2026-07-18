# James-Stein type (conditional-expectation) estimation
#
# Burghgraeve, E., De Neve, J., & Rosseel, Y. (2021). Estimating structural
# equation models using James-Stein type shrinkage estimators. Psychometrika,
# 86(1), 96-130.
#
# In every (measurement or structural) equation, the latent variables are
# replaced by an estimate of their conditional expectation given observed
# proxies: the scaling indicators (estimator = "JS"), or reliability-weighted
# aggregates of the indicators (estimator = "JSA"); the equation is then
# estimated by least squares. Everything is computed from the sample moments
# ('attenuation-corrected' least squares): regressing on the estimated
# conditional expectations is algebraically identical to correcting the
# normal equations for the measurement-error variance of the proxies.
#
# The equation decomposition (dependent variables, scaling-indicator
# substitution) is shared with the (MI)IV estimator (lav_model_find_iv);
# the two estimators differ only in how each equation is estimated (2SLS
# with instruments versus least squares on conditional expectations).
# The undirected (variance/covariance) parameters are estimated by the
# same second stage as the IV estimator (lav_sem_miiv_varcov).
#
# YR/EB 08 Feb 2023: first version (lav_cfa_jamesstein.R), cfa only
# YR    18 Jul 2026: rewritten in moment form: structural models, observed
#                    covariates, mean structure; stage 2 via the IV
#                    machinery; closed-form (alphamax) aggregation weights

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_js_internal <- function(lavmodel = NULL, lavh1 = NULL,
                                lavsamplestats = NULL,
                                lavpartable = NULL,
                                lavdata = NULL, lavoptions = NULL) {
  lavpta <- lav_pt_attributes(lavpartable)
  lavpartable <- lav_pt_set_cache(lavpartable, lavpta)

  # aggregated (JSA) or scaling (JS)?
  aggregated <- identical(lavoptions$estimator, "JSA")

  # model-class checks
  lav_sem_js_check(
    lavmodel = lavmodel, lavpartable = lavpartable, lavpta = lavpta,
    lavdata = lavdata, lavoptions = lavoptions, aggregated = aggregated
  )

  # equations: reuse the (MI)IV decomposition (each dependent variable, with
  # the latent variables replaced by their scaling indicators); the
  # instrument sets (eq$iv/eq$miiv) are ignored here
  eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)

  # directed versus undirected (free) parameters
  undirected_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  directed_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~", "~1"))
  free_directed_idx <- unique(lavpartable$free[directed_idx])
  free_undirected_idx <- unique(lavpartable$free[undirected_idx])

  # run the (complete) estimation map at the sample moments
  x <- lav_sem_js_estimate(
    vec = NULL, eqs = eqs,
    lavmodel = lavmodel, lavpartable = lavpartable,
    lavsamplestats = lavsamplestats, lavh1 = lavh1,
    lavoptions = lavoptions,
    free_directed_idx = free_directed_idx,
    free_undirected_idx = free_undirected_idx,
    aggregated = aggregated
  )
  eqs <- attr(x, "eqs")
  x <- as.numeric(x)

  # apply bounds (if any)
  if (!is.null(lavpartable$lower)) {
    lower_x <- lavpartable$lower[lavpartable$free > 0 &
      !duplicated(lavpartable$free)]
    too_small_idx <- which(x < lower_x)
    if (length(too_small_idx) > 0L) {
      x[too_small_idx] <- lower_x[too_small_idx]
    }
  }
  if (!is.null(lavpartable$upper)) {
    upper_x <- lavpartable$upper[lavpartable$free > 0 &
      !duplicated(lavpartable$free)]
    too_large_idx <- which(x > upper_x)
    if (length(too_large_idx) > 0L) {
      x[too_large_idx] <- upper_x[too_large_idx]
    }
  }

  attr(x, "eqs") <- eqs
  x
}


# the COMPLETE estimation map, as a function of the sample moments: stage 1
# (directed parameters), stage 2 (variances/covariances) and the joint mean
# solve. When 'vec' is NULL, the lavh1 moments are used (the point
# estimates); when 'vec' contains a stacked moment vector (as produced by
# lav_implied_to_vec), the same map is evaluated at those moments -- this
# is the function that is differentiated (numerically) for the standard
# errors, so that the uncertainty of ALL ingredients (the Spearman
# reliability estimates, the aggregation weights, the pooled solves, the
# stage-2 weights and the mean solve) is propagated automatically
lav_sem_js_estimate <- function(vec = NULL, eqs = NULL,
                                lavmodel = NULL, lavpartable = NULL,
                                lavsamplestats = NULL, lavh1 = NULL,
                                lavoptions = NULL,
                                free_directed_idx = NULL,
                                free_undirected_idx = NULL,
                                aggregated = FALSE) {
  samplestats <- !is.null(vec)
  if (samplestats) {
    implied <- lav_vec_to_implied(vec, lavmodel = lavmodel)
  } else {
    implied <- lavh1$implied
  }

  # initial parameter vector
  x <- lav_model_get_parameters(lavmodel)

  ########################################
  # first stage: directed parameters     #
  ########################################
  eqs_out <- eqs
  if (length(free_directed_idx) > 0L) {
    theta1 <- lav_sem_js_stage1_samp(
      x = vec, samplestats = samplestats, eqs = eqs,
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats, lavh1 = lavh1,
      lavoptions = lavoptions,
      free_directed_idx = free_directed_idx,
      aggregated = aggregated
    )
    eqs_out <- attr(theta1, "eqs")
    x[free_directed_idx] <- as.numeric(theta1)
  }

  #######################################
  # second stage: undirected parameters #
  #######################################
  js_varcov_method <- toupper(lavoptions$estimator.args[["js_varcov_method"]])
  if (length(js_varcov_method) == 0L) {
    js_varcov_method <- "RLS"
  }
  if (length(free_undirected_idx) > 0L && js_varcov_method != "NONE") {
    lavmodel_tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
    delta_list <- lav_sem_miiv_delta(lavmodel_tmp)
    if (lavmodel@nblocks == 1L) {
      theta2 <- lav_sem_miiv_varcov_block(
        b = 1L, fu = free_undirected_idx, delta_b = delta_list[[1L]],
        implied = implied, lavmodel = lavmodel, x = x,
        iv_varcov_method = js_varcov_method, return_h = FALSE
      )
      x[free_undirected_idx] <- as.numeric(theta2)
    } else {
      # multiple groups: each group's moments contribute its own block
      # (as for the IV estimator; cross-group constraints among the
      # variances, if any, are handled per block)
      pt_block <- if (!is.null(lavpartable$block)) {
        lavpartable$block
      } else {
        rep(1L, length(lavpartable$op))
      }
      for (b in seq_len(lavmodel@nblocks)) {
        fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
          lavpartable$free > 0L & pt_block == b])
        fu_b <- fu_b[fu_b %in% free_undirected_idx]
        if (length(fu_b) == 0L) {
          next
        }
        tb <- lav_sem_miiv_varcov_block(
          b = b, fu = fu_b, delta_b = delta_list[[b]],
          implied = implied, lavmodel = lavmodel, x = x,
          iv_varcov_method = js_varcov_method, return_h = FALSE
        )
        x[fu_b] <- as.numeric(tb)
      }
    }
  } else if (length(free_undirected_idx) > 0L) {
    x[free_undirected_idx] <- as.numeric(NA)
  }

  # mean structure: by default, re-estimate all free mean parameters
  # (observed intercepts and latent means) jointly by GLS, exactly as the IV
  # estimator does; this also provides the means of the exogenous latent
  # variables (which have no equation of their own)
  js_mean_structure <- tolower(lavoptions$estimator.args[["js_mean_structure"]])
  if (length(js_mean_structure) == 0L) {
    js_mean_structure <- "wls"
  }
  if (lavmodel@meanstructure && identical(js_mean_structure, "wls")) {
    free_mean_idx <- unique(lavpartable$free[lavpartable$op == "~1" &
      lavpartable$free > 0L & !duplicated(lavpartable$free)])
    if (length(free_mean_idx) > 0L) {
      x[free_mean_idx] <- lav_sem_miiv_mean_wls(
        lavmodel = lavmodel, lavsamplestats = lavsamplestats, x = x,
        free_mean_idx = free_mean_idx,
        sample_mean = implied$mean
      )
    }
  }

  attr(x, "eqs") <- eqs_out
  x
}


# VCOV for the free parameters: the delta method over the sample moments.
# The full estimation map theta(s) is differentiated numerically with
# respect to the stacked moment vector s = [mean, vech(cov)], and
# sandwiched with the asymptotic covariance of the sample moments:
#
#   vcov = J Gamma J' / n
#
# with Gamma the normal-theory moment covariance (js_gamma = "nt", the
# default; evaluated at the model-implied moments unless
# js_vcov_gamma_modelbased = FALSE) or the distribution-free (ADF) moment
# covariance (js_gamma = "adf", requires raw data). Because the Jacobian
# runs through the complete map, the extra variability caused by
# estimating the reliability matrix (Burghgraeve et al. 2021, Theorem 2)
# -- and, for JSA, the aggregation weights -- is included, and the
# directed/undirected cross-covariances are available (e.g. for defined
# parameters).
lav_sem_js_vcov <- function(lavmodel = NULL, lavsamplestats = NULL,
                            lavoptions = NULL, lavpartable = NULL,
                            lavimplied = NULL,
                            lavh1 = NULL, lavdata = NULL, eqs = NULL) {
  lavpta <- lav_pt_attributes(lavpartable)
  lavpartable <- lav_pt_set_cache(lavpartable, lavpta)
  aggregated <- identical(lavoptions$estimator, "JSA")

  # directed versus undirected (free) parameters
  undirected_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) &
    lavpartable$op == "~~")
  directed_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) &
    lavpartable$op %in% c("=~", "~", "~1"))
  free_directed_idx <- unique(lavpartable$free[directed_idx])
  free_undirected_idx <- unique(lavpartable$free[undirected_idx])

  if (is.null(eqs)) {
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
  }

  # Jacobian of the complete estimation map over the stacked moments
  vec0 <- lav_implied_to_vec(
    implied = lavh1$implied, lavmodel = lavmodel, drop_list = TRUE
  )
  # r = 2 Richardson extrapolation is accurate to ~1e-6 relative here and
  # halves the number of map evaluations (the map itself is smooth)
  jac <- numDeriv::jacobian(
    func = function(v) {
      as.numeric(lav_sem_js_estimate(
        vec = v, eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavsamplestats = lavsamplestats, lavh1 = lavh1,
        lavoptions = lavoptions,
        free_directed_idx = free_directed_idx,
        free_undirected_idx = free_undirected_idx,
        aggregated = aggregated
      ))
    },
    x = vec0, method.args = list(r = 2L)
  )

  # moment covariance: block-diagonal over the groups,
  # gamma_big = bdiag(Gamma_b / nobs_b)
  js_gamma <- tolower(lavoptions$estimator.args[["js_gamma"]])
  if (length(js_gamma) == 0L) {
    js_gamma <- "nt"
  }
  gamma_modelbased <-
    !isFALSE(lavoptions$estimator.args[["js_vcov_gamma_modelbased"]])
  two_stage <- any(lavoptions$missing == c("two.stage", "robust.two.stage"))
  if (js_gamma == "adf") {
    if (lavdata@data.type != "full" || two_stage) {
      lav_msg_warn(gettext(
        "[JS] js_gamma = \"adf\" requires complete raw data; using the
         default moment covariance instead."))
      js_gamma <- "nt"
    }
  }
  nblocks <- lavmodel@nblocks
  gamma_g <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    if (two_stage) {
      # two-stage missing data: the moments are the (saturated) EM
      # estimates; their asymptotic covariance is the inverse saturated
      # information (two.stage; Savalei & Bentler 2009) or its sandwich
      # (robust.two.stage; Savalei & Falk 2014)
      mu_b <- lavh1$implied$mean[[b]]
      sig_b <- lavh1$implied$cov[[b]]
      if (identical(tolower(lavoptions$missing), "robust.two.stage")) {
        gg <- lav_mvn_mi_h1_omega_sw(
          y = lavdata@X[[b]], mp = lavdata@Mp[[b]],
          yp = lavsamplestats@missing[[b]],
          mu = mu_b, sigma_1 = sig_b,
          x_idx = integer(0L), information = "observed"
        )
      } else {
        i1 <- lav_mvnorm_missing_information_observed_samplestats(
          yp = lavsamplestats@missing[[b]],
          mu = mu_b, sigma_1 = sig_b, x_idx = integer(0L)
        )
        gg <- lav_mat_sym_inverse(i1)
      }
      gamma_g[[b]] <- gg / lavsamplestats@nobs[[b]]
    } else if (js_gamma == "adf") {
      gamma_g[[b]] <- lav_samp_gamma(
        m_y = lavdata@X[[b]],
        meanstructure = lavmodel@meanstructure
      ) / lavsamplestats@nobs[[b]]
    } else {
      cov_b <- NULL
      if (gamma_modelbased) {
        cov_b <- lavimplied$cov[[b]]
        # the model-implied covariance may be non-PD (eg a Heywood case)
        if (inherits(try(chol(cov_b), silent = TRUE), "try-error")) {
          cov_b <- NULL
        }
      }
      if (is.null(cov_b)) {
        cov_b <- lavh1$implied$cov[[b]]
      }
      nd_b <- if (lavmodel@meanstructure) {
        nrow(cov_b) + nrow(cov_b) * (nrow(cov_b) + 1L) / 2L
      } else {
        nrow(cov_b) * (nrow(cov_b) + 1L) / 2L
      }
      gamma_g[[b]] <- lav_mat_k_gammant_kt(
        m_k = diag(nd_b), s = cov_b,
        meanstructure = lavmodel@meanstructure,
        x_idx = integer(0L)
      ) / lavsamplestats@nobs[[b]]
    }
  }
  gamma_big <- if (nblocks == 1L) {
    gamma_g[[1L]]
  } else {
    lav_mat_bdiag(gamma_g)
  }

  vcov <- jac %*% gamma_big %*% t(jac)
  vcov <- (vcov + t(vcov)) / 2

  # the rest of lavaan expects the vcov of a ceq.simple.only model in the
  # 'unco' space (one row/column per non-collapsed parameter); the vcov
  # above is built in the compact (nx.free) space, so expand it via the
  # ceq.simple.K mapping (as for the IV estimator)
  if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L) {
    vcov <- lavmodel@ceq.simple.K %*% vcov %*% t(lavmodel@ceq.simple.K)
  }

  vcov
}


# model-class checks for the JS/JSA estimators; stop with a clear message
# for anything the (phase-1) implementation does not cover
lav_sem_js_check <- function(lavmodel = NULL, lavpartable = NULL,
                             lavpta = NULL, lavdata = NULL,
                             lavoptions = NULL, aggregated = FALSE) {
  label <- if (aggregated) "JSA" else "JS"

  stopifnot(lavdata@nlevels == 1L)
  if (lavmodel@categorical) {
    lav_msg_stop(gettextf(
      "estimator %s does not support categorical data (yet).", label))
  }
  if (lavoptions$std.lv) {
    lav_msg_stop(gettextf(
      "estimator %s requires marker/scaling indicators;
       std.lv = TRUE is not supported.", label))
  }
  if (lavmodel@conditional.x) {
    lav_msg_stop(gettextf(
      "estimator %s does not support conditional.x = TRUE (yet).", label))
  }
  if (lavmodel@correlation) {
    lav_msg_stop(gettextf(
      "estimator %s does not support correlation structures (yet).", label))
  }
  if (isTRUE(lavmodel@composites)) {
    lav_msg_stop(gettextf(
      "estimator %s does not support composites (yet).", label))
  }
  if (lavmodel@nefa > 0L) {
    lav_msg_stop(gettextf(
      "estimator %s does not support efa blocks (yet).", label))
  }
  if (any(!vapply(lavdata@weights, is.null, logical(1L)))) {
    lav_msg_stop(gettextf(
      "estimator %s does not support sampling weights (yet).", label))
  }

  # the structural part must be recursive: the conditional expectations are
  # only valid regressors if the disturbance of each equation is
  # uncorrelated with (functions of) its right-hand-side variables
  if (any(vapply(lavmodel@modprop$acyclic, isFALSE, logical(1L)))) {
    lav_msg_stop(gettextf(
      "estimator %s requires a recursive structural model
       (no feedback loops).", label))
  }

  pt_block <- if (!is.null(lavpartable$block)) {
    lavpartable$block
  } else {
    rep(1L, length(lavpartable$op))
  }
  for (b in seq_len(lavmodel@nblocks)) {
    lav_sem_js_check_block(
      lavpartable = lavpartable, lavpta = lavpta, b = b,
      pt_block = pt_block, label = label, aggregated = aggregated
    )
  }

  invisible(TRUE)
}


# per-block model-structure checks
lav_sem_js_check_block <- function(lavpartable = NULL, lavpta = NULL,
                                   b = 1L, pt_block = NULL, label = "JS",
                                   aggregated = FALSE) {
  in_b <- pt_block == b
  lv_names <- lavpta$vnames$lv.regular[[b]]

  # higher-order factors: a latent variable measured by latent variables
  if (any(lavpartable$op == "=~" & in_b &
          lavpartable$rhs %in% lv_names)) {
    lav_msg_stop(gettextf(
      "estimator %s does not support higher-order factors (yet).", label))
  }

  # active covariance rows (free, or fixed to a nonzero value)
  cov_idx <- which(lavpartable$op == "~~" & in_b &
    lavpartable$lhs != lavpartable$rhs &
    (lavpartable$free > 0L |
      (!is.na(lavpartable$ustart) & lavpartable$ustart != 0)))

  # no disturbance of a dependent variable may be correlated with one of its
  # regressors, or with any ancestor of its regressors (the conditional
  # expectations condition on proxies of these variables)
  reg_idx <- which(lavpartable$op == "~" & in_b)
  if (length(reg_idx) > 0L && length(cov_idx) > 0L) {
    # ancestor sets via the regression graph (edge: rhs -> lhs)
    parents <- split(
      lavpartable$rhs[reg_idx],
      factor(lavpartable$lhs[reg_idx],
             levels = unique(lavpartable$lhs[reg_idx]))
    )
    ancestors <- function(vars) {
      anc <- character(0L)
      todo <- vars
      while (length(todo) > 0L) {
        v <- todo[1L]
        todo <- todo[-1L]
        if (v %in% anc) {
          next
        }
        anc <- c(anc, v)
        todo <- c(todo, parents[[v]])
      }
      anc
    }
    for (e in unique(lavpartable$lhs[reg_idx])) {
      rhs_e <- lavpartable$rhs[reg_idx][lavpartable$lhs[reg_idx] == e]
      anc_e <- ancestors(rhs_e)
      bad <- which((lavpartable$lhs[cov_idx] == e &
                      lavpartable$rhs[cov_idx] %in% anc_e) |
                   (lavpartable$rhs[cov_idx] == e &
                      lavpartable$lhs[cov_idx] %in% anc_e))
      if (length(bad) > 0L) {
        lav_msg_stop(gettextf(
          "estimator %1$s is not consistent when the disturbance of %2$s is
           correlated with (an ancestor of) one of its regressors: %3$s.",
          label, e,
          lav_msg_view(unique(c(lavpartable$lhs[cov_idx][bad],
                                lavpartable$rhs[cov_idx][bad])), "none")))
      }
    }
  }

  # residual covariances among the observed variables are handled by
  # per-equation clean conditioning sets (see lav_sem_js_eq_plans) and by
  # masking the contaminated tetrads in the reliability estimates (see
  # lav_sem_js_theta); unidentifiable configurations (no clean proxy, no
  # clean tetrad) are caught there with informative messages

  invisible(TRUE)
}


# residual-covariance ('contamination') map for one block: a symmetric
# logical nvar x nvar matrix flagging pairs of residual-bearing observed
# variables (indicators or endogenous observed variables) with an active
# (free, or fixed to a nonzero value) residual covariance. Such a pair
# invalidates any conditioning set that combines one member with the
# other member's equation, and any Spearman tetrad through the pair.
lav_sem_js_contam <- function(lavpartable = NULL, pt_block = NULL, b = 1L,
                              ov_names = NULL, lavpta = NULL) {
  nvar <- length(ov_names)
  contam <- matrix(FALSE, nrow = nvar, ncol = nvar)
  in_b <- pt_block == b
  ind_names <- unique(lavpartable$rhs[lavpartable$op == "=~" & in_b])
  res_names <- unique(c(ind_names, lavpta$vnames$eqs.y[[b]]))
  cov_idx <- which(lavpartable$op == "~~" & in_b &
    lavpartable$lhs != lavpartable$rhs &
    (lavpartable$free > 0L |
      (!is.na(lavpartable$ustart) & lavpartable$ustart != 0)) &
    lavpartable$lhs %in% res_names & lavpartable$rhs %in% res_names)
  if (length(cov_idx) > 0L) {
    ii <- match(lavpartable$lhs[cov_idx], ov_names)
    jj <- match(lavpartable$rhs[cov_idx], ov_names)
    ok <- !is.na(ii) & !is.na(jj)
    contam[cbind(ii[ok], jj[ok])] <- TRUE
    contam[cbind(jj[ok], ii[ok])] <- TRUE
  }
  contam
}


# Spearman/communality residual variances with tetrad masking: like
# lav_cfa_theta_spearman (h2_i = mean over tetrads r_ia * r_ib / r_ab),
# but tetrads that run through a contaminated pair (an active residual
# covariance) are excluded from the averaging. When a variable has no
# clean within-factor tetrad left, the tetrads are formed with EXTERNAL
# variables z instead (Kano-1990 style: any observed model variable
# outside the factor whose residual is uncorrelated with both members;
# the formula r_ia * r_iz / r_az identifies the communality through the
# factor correlations). Returns NA for a variable without any clean
# tetrad at all. With no contaminated pairs within the factor block,
# this reproduces lav_cfa_theta_spearman exactly.
#
# s_full/contam cover ALL observed variables of the block; idx is the
# factor's indicator set; zpool the external candidates.
lav_sem_js_theta_spearman_masked <- function(s_full = NULL, idx = NULL,
                                             contam = NULL,
                                             zpool = integer(0L),
                                             bounds = "wide") {
  if (!any(contam[idx, idx])) {
    return(lav_cfa_theta_spearman(s_full[idx, idx, drop = FALSE],
      bounds = bounds))
  }
  r <- cov2cor(s_full)
  p <- length(idx)
  out <- rep(as.numeric(NA), p)
  for (k in seq_len(p)) {
    i <- idx[k]
    others <- setdiff(idx, i)
    # clean within-factor tetrads
    vals <- numeric(0L)
    if (length(others) >= 2L) {
      pairs <- utils::combn(others, 2L)
      for (m in seq_len(ncol(pairs))) {
        a <- pairs[1L, m]
        b <- pairs[2L, m]
        if (contam[i, a] || contam[i, b] || contam[a, b]) {
          next
        }
        vals <- c(vals, r[i, a] * r[i, b] / r[a, b])
      }
    }
    # fall back to external tetrads when no clean within-factor tetrad
    if (length(vals) == 0L && length(zpool) > 0L) {
      for (a in others) {
        if (contam[i, a]) {
          next
        }
        for (z in zpool) {
          if (contam[i, z] || contam[a, z]) {
            next
          }
          vals <- c(vals, r[i, a] * r[i, z] / r[a, z])
        }
      }
    }
    if (length(vals) == 0L) {
      next # no clean tetrad: NA
    }
    h2 <- mean(vals)
    if (bounds == "standard") {
      h2[h2 < 0] <- 0
      h2[h2 > 1] <- 1
    } else if (bounds == "wide") {
      h2[h2 < -0.05] <- -0.05
      h2[h2 > +1.20] <- +1.20
    }
    out[k] <- (1 - h2) * s_full[i, i]
  }
  out
}


# residual (measurement-error) variances of the observed indicators, needed
# for the shrinkage (and, for JSA, the aggregation weights); by default the
# Spearman/communality estimator applied per factor (Burghgraeve et al.,
# 2021, Appendix B), with contaminated tetrads masked (see
# lav_sem_js_theta_spearman). Indicators without a clean tetrad get NA;
# whether that is a problem depends on whether they are needed as a
# conditioning proxy (checked in lav_sem_js_eq_plans).
lav_sem_js_theta <- function(s = NULL, ind_pat = NULL, ov_names = NULL,
                             lv_names = NULL, lavoptions = NULL,
                             aggregated = FALSE, contam = NULL) {
  label <- if (aggregated) "JSA" else "JS"
  nvar <- ncol(s)
  ea <- lavoptions$estimator.args
  js_theta <- tolower(ea[["js_theta"]])
  if (length(js_theta) == 0L) {
    js_theta <- "spearman"
  }
  js_theta_bounds <- tolower(ea[["js_theta_bounds"]])
  if (length(js_theta_bounds) == 0L) {
    js_theta_bounds <- "wide"
  }

  theta <- rep(as.numeric(NA), nvar)
  if (js_theta == "user") {
    values <- ea[["js_theta_values"]]
    ind_idx <- which(rowSums(ind_pat != 0) > 0L)
    missing_idx <- ind_idx[!ov_names[ind_idx] %in% names(values)]
    if (length(missing_idx) > 0L) {
      lav_msg_stop(gettextf(
        "js_theta_values should provide a residual variance for every
         indicator; missing: %s.",
        lav_msg_view(ov_names[missing_idx], "none")))
    }
    theta[ind_idx] <- as.numeric(values[ov_names[ind_idx]])
  } else {
    # spearman, per factor
    for (f in seq_len(ncol(ind_pat))) {
      idx <- which(ind_pat[, f] != 0)
      if (length(idx) < 3L) {
        lav_msg_stop(gettextf(
          "estimator %1$s needs at least 3 indicators per factor to estimate
           the reliability of the scaling indicators; factor %2$s has only
           %3$d. Alternatively, provide the residual variances directly via
           estimator.args = list(js_theta = \"user\", js_theta_values = ...).",
          label, lv_names[f], length(idx)))
      }
      if (is.null(contam) || !any(contam[idx, idx])) {
        theta[idx] <- lav_cfa_theta_spearman(s[idx, idx, drop = FALSE],
          bounds = js_theta_bounds
        )
      } else {
        # contaminated within-factor tetrads: mask them; when a variable
        # has none left, use external tetrads (all observed variables
        # outside this factor as candidate z's)
        theta[idx] <- lav_sem_js_theta_spearman_masked(
          s_full = s, idx = idx, contam = contam,
          zpool = setdiff(seq_len(nvar), idx),
          bounds = js_theta_bounds
        )
      }
    }
  }

  # non-indicator variables (observed covariates, endogenous observed
  # variables) are error-free proxies of themselves; indicators without a
  # clean tetrad keep NA
  nonind_idx <- which(rowSums(ind_pat != 0) == 0L)
  theta[nonind_idx] <- 0

  # keep the variances within [0, var(y)] (NA-safe)
  diag_s <- diag(s)
  theta[which(theta < 0)] <- 0
  too_large_idx <- which(theta > diag_s)
  if (length(too_large_idx) > 0L) {
    theta[too_large_idx] <- diag_s[too_large_idx]
  }

  theta
}


# reliability-maximizing (alphamax; Bentler, 1968) aggregation weights: the
# maximizer of (w's w - w'diag(theta) w) / (w's w) is the generalized
# eigenvector of (diag(theta), s) with the smallest eigenvalue; normalized
# so that the weights sum to one
lav_sem_js_weights <- function(s = NULL, theta = NULL) {
  p <- ncol(s)
  if (p == 1L) {
    return(1)
  }
  equal_w <- rep(1 / p, p)
  r_mat <- try(chol(s), silent = TRUE)
  if (inherits(r_mat, "try-error")) {
    return(equal_w)
  }
  r_inv <- backsolve(r_mat, diag(p))
  m_sym <- crossprod(r_inv, theta * r_inv) # t(Ri) %*% diag(theta) %*% Ri
  ev <- eigen((m_sym + t(m_sym)) / 2, symmetric = TRUE)
  w <- drop(r_inv %*% ev$vectors[, p]) # smallest eigenvalue
  sw <- sum(w)
  if (!is.finite(sw) || abs(sw) < .Machine$double.eps^0.5) {
    return(equal_w)
  }
  w / sw
}


# build the per-equation conditioning-set 'plans' for one block: for every
# latent regressor of every equation, the set of observed variables used as
# its conditioning proxy. The plans depend only on the MODEL STRUCTURE
# (which residual covariances are active), never on the data, so the
# estimation map remains a smooth function of the sample moments (needed
# for the delta-method standard errors).
#
# Default proxies: the scaling indicator (scaling mode) or the pure
# indicators of the factor excluding the dependent variable (aggregated
# mode). A proxy variable is invalid for an equation when it has an active
# residual covariance with the equation's dependent variable ('tainted'),
# when its own reliability cannot be estimated (theta NA), or when it has
# an active residual covariance with a proxy variable used elsewhere in
# the same conditioning set (which would contaminate the joint conditional
# expectation). Invalid variables are replaced by clean pure indicators of
# the same factor whose loadings are identified (fixed, or estimable from
# their own clean equation). An equation is 'deferred' when any of its
# proxies is not simply the (fixed-loading) scaling indicator: such
# equations need loading estimates from the clean equations and are solved
# in a second pass.
lav_sem_js_eq_plans <- function(eqs_b = NULL, aggregated = FALSE,
                                contam = NULL, ind_pat = NULL,
                                single_pat = NULL, markers = NULL,
                                ov_names = NULL, lv_names = NULL,
                                theta = NULL, lavpartable = NULL,
                                pt_block = NULL, b = 1L, label = "JS") {
  nvar <- length(ov_names)
  theta_na <- which(is.na(theta))
  any_contam <- any(contam)

  # loading availability: lambda(v, f) is identified before the deferred
  # pass when it is fixed, or when v's own (single-proxy) measurement
  # equation conditions on a clean scaling indicator
  mm_idx <- which(lavpartable$op == "=~" & pt_block == b)
  lam_free <- matrix(FALSE, nrow = nvar, ncol = length(lv_names))
  if (length(mm_idx) > 0L) {
    lam_free[cbind(
      match(lavpartable$rhs[mm_idx], ov_names),
      match(lavpartable$lhs[mm_idx], lv_names)
    )] <- lavpartable$free[mm_idx] > 0L
  }
  lam_available <- function(v, f) {
    if (!lam_free[v, f]) {
      return(TRUE) # fixed value
    }
    m_f <- match(markers[[lv_names[f]]], ov_names)
    !contam[m_f, v] && !v %in% theta_na
  }

  plans <- vector("list", length(eqs_b))
  deferred <- logical(length(eqs_b))
  for (j in seq_along(eqs_b)) {
    eq <- eqs_b[[j]]
    if (identical(eq$rhs_new, "1")) {
      next
    }
    y_idx <- match(eq$lhs_new, ov_names)
    is_lv <- eq$rhs %in% lv_names
    if (!any(is_lv)) {
      plans[[j]] <- vector("list", length(eq$rhs))
      next
    }
    tainted <- if (any_contam) which(contam[, y_idx]) else integer(0L)

    # candidate pool per latent regressor: clean pure indicators with an
    # estimable reliability and an identified loading
    fac <- match(eq$rhs, lv_names) # NA for observed regressors
    cand <- vector("list", length(eq$rhs))
    marker_idx <- rep(NA_integer_, length(eq$rhs))
    for (i in which(is_lv)) {
      f <- fac[i]
      marker_idx[i] <- match(markers[[eq$rhs[i]]], ov_names)
      pool <- which(ind_pat[, f] != 0 & single_pat)
      pool <- setdiff(pool, c(y_idx, tainted, theta_na))
      pool <- pool[vapply(pool, lam_available, logical(1L), f = f)]
      cand[[i]] <- pool
    }

    # initial sets
    sets <- vector("list", length(eq$rhs))
    for (i in which(is_lv)) {
      if (!aggregated) {
        sets[[i]] <- if (marker_idx[i] %in% cand[[i]]) {
          marker_idx[i]
        } else {
          cand[[i]]
        }
      } else {
        sets[[i]] <- cand[[i]]
      }
    }

    # fixpoint: drop (both members of) any contaminated pair within or
    # across the chosen sets; dropped variables never come back, so this
    # terminates
    if (any_contam) {
      dropped <- integer(0L)
      repeat {
        used <- unlist(sets[is_lv])
        pair_idx <- which(contam[used, used, drop = FALSE], arr.ind = TRUE)
        if (nrow(pair_idx) == 0L) {
          break
        }
        bad_vars <- unique(used[unique(as.vector(pair_idx))])
        dropped <- c(dropped, bad_vars)
        for (i in which(is_lv)) {
          sets[[i]] <- setdiff(sets[[i]], dropped)
          if (length(sets[[i]]) == 0L) {
            # refill from the candidate pool (minus everything dropped)
            sets[[i]] <- setdiff(cand[[i]], dropped)
          }
        }
        if (any(vapply(sets[is_lv], length, integer(1L)) == 0L)) {
          break
        }
      }
    }

    # identifiable?
    for (i in which(is_lv)) {
      if (length(sets[[i]]) == 0L) {
        lav_msg_stop(gettextf(
          "estimator %1$s cannot construct a clean conditioning proxy for
           %2$s in the equation for %3$s: too many residual covariances
           involve its indicators. Consider fixing some residual
           covariances, or provide reliabilities via js_theta_values.",
          label, eq$rhs[i], eq$lhs_new))
      }
    }

    plans[[j]] <- sets
    # deferred: any proxy that is not simply the fixed-loading scaling
    # indicator (the plain-marker solve needs no loading estimates)
    deferred[j] <- any(vapply(which(is_lv), function(i) {
      !(length(sets[[i]]) == 1L && sets[[i]] == marker_idx[i] &&
          !lam_free[marker_idx[i], fac[i]])
    }, logical(1L)))
  }

  list(plans = plans, deferred = deferred)
}


# stage 1: estimate all directed (loading/regression/intercept) parameters,
# equation by equation, by least squares on the estimated conditional
# expectations -- computed directly from the sample moments.
# if samplestats = TRUE, x should contain the stacked sample statistics
# (as produced by lav_implied_to_vec); this mode is used to compute the
# moment Jacobian for the standard errors
lav_sem_js_stage1_samp <- function(x = NULL, samplestats = FALSE,
                                   eqs = NULL, lavmodel = NULL,
                                   lavpartable = NULL,
                                   lavsamplestats = NULL, lavh1 = NULL,
                                   lavoptions = NULL,
                                   free_directed_idx = NULL,
                                   aggregated = FALSE) {
  # sample moments
  if (samplestats) {
    implied <- lav_vec_to_implied(x, lavmodel = lavmodel)
  } else {
    implied <- lavh1$implied
  }
  x <- lav_model_get_parameters(lavmodel)

  lavpta <- lav_pt_attributes(lavpartable)
  pt_block <- if (!is.null(lavpartable$block)) {
    lavpartable$block
  } else {
    rep(1L, length(lavpartable$op))
  }

  for (b in seq_len(lavmodel@nblocks)) {
    s_mat <- implied$cov[[b]]
    m_vec <- if (lavmodel@meanstructure) implied$mean[[b]] else NULL
    n <- lavsamplestats@nobs[[b]]
    ov_names <- lavpta$vnames$ov[[b]]
    lv_names <- lavpta$vnames$lv.regular[[b]]
    markers <- lavpta$vnames$lv.marker[[b]] # named by lv
    nvar <- length(ov_names)
    nfac <- length(lv_names)

    # indicator pattern (nvar x nfac)
    ind_pat <- matrix(0, nrow = nvar, ncol = nfac)
    mm_idx <- which(lavpartable$op == "=~" & pt_block == b)
    if (length(mm_idx) > 0L && nfac > 0L) {
      ind_pat[cbind(
        match(lavpartable$rhs[mm_idx], ov_names),
        match(lavpartable$lhs[mm_idx], lv_names)
      )] <- 1
    }

    single_pat <- rowSums(ind_pat != 0) == 1L
    label <- if (aggregated) "JSA" else "JS"

    # residual-covariance (contamination) map + residual variances of the
    # proxies (per block; contaminated tetrads are masked)
    contam <- lav_sem_js_contam(
      lavpartable = lavpartable, pt_block = pt_block, b = b,
      ov_names = ov_names, lavpta = lavpta
    )
    theta <- lav_sem_js_theta(
      s = s_mat, ind_pat = ind_pat, ov_names = ov_names,
      lv_names = lv_names,
      lavoptions = lavoptions, aggregated = aggregated, contam = contam
    )

    # small-sample (Efron-Morris) shrinkage factor
    csmall <- (n - 3) / (n - 1)
    if (isFALSE(lavoptions$estimator.args[["js_small_sample"]])) {
      csmall <- 1
    }

    # conditioning-set plans for the scaling passes
    pl <- lav_sem_js_eq_plans(
      eqs_b = eqs[[b]], aggregated = FALSE, contam = contam,
      ind_pat = ind_pat, single_pat = single_pat, markers = markers,
      ov_names = ov_names, lv_names = lv_names, theta = theta,
      lavpartable = lavpartable, pt_block = pt_block, b = b, label = label
    )

    # pass 1a: scaling; the equations that condition on plain (fixed-
    # loading) scaling indicators -- this also provides the loadings needed
    # to scale the other proxies
    lambda_mat <- lav_sem_js_lambda_mat(
      x = x, lavpartable = lavpartable, mm_idx = mm_idx,
      ov_names = ov_names, lv_names = lv_names
    )
    out <- lav_sem_js_eqs_pass(
      eqs_b = eqs[[b]], plans = pl$plans, only = which(!pl$deferred),
      s_mat = s_mat, m_vec = m_vec, n = n,
      ov_names = ov_names, lv_names = lv_names, markers = markers,
      contam = contam, theta = theta, csmall = csmall,
      lambda_mat = lambda_mat, lavpartable = lavpartable,
      lavmodel = lavmodel, x = x
    )
    x <- out$x
    eqs[[b]] <- out$eqs

    # pass 1b: the deferred equations (clean-proxy conditioning sets that
    # need the pass-1a loading estimates)
    if (any(pl$deferred)) {
      lambda_mat <- lav_sem_js_lambda_mat(
        x = x, lavpartable = lavpartable, mm_idx = mm_idx,
        ov_names = ov_names, lv_names = lv_names
      )
      out <- lav_sem_js_eqs_pass(
        eqs_b = eqs[[b]], plans = pl$plans, only = which(pl$deferred),
        s_mat = s_mat, m_vec = m_vec, n = n,
        ov_names = ov_names, lv_names = lv_names, markers = markers,
        contam = contam, theta = theta, csmall = csmall,
        lambda_mat = lambda_mat, lavpartable = lavpartable,
        lavmodel = lavmodel, x = x
      )
      x <- out$x
      eqs[[b]] <- out$eqs
    }

    # pass 2 (JSA): aggregated proxies, scaled by the scaling-pass loadings
    if (aggregated) {
      pla <- lav_sem_js_eq_plans(
        eqs_b = eqs[[b]], aggregated = TRUE, contam = contam,
        ind_pat = ind_pat, single_pat = single_pat, markers = markers,
        ov_names = ov_names, lv_names = lv_names, theta = theta,
        lavpartable = lavpartable, pt_block = pt_block, b = b, label = label
      )
      lambda_mat <- lav_sem_js_lambda_mat(
        x = x, lavpartable = lavpartable, mm_idx = mm_idx,
        ov_names = ov_names, lv_names = lv_names
      )
      out <- lav_sem_js_eqs_pass(
        eqs_b = eqs[[b]], plans = pla$plans,
        only = seq_along(eqs[[b]]),
        s_mat = s_mat, m_vec = m_vec, n = n,
        ov_names = ov_names, lv_names = lv_names, markers = markers,
        contam = contam, theta = theta, csmall = csmall,
        lambda_mat = lambda_mat, lavpartable = lavpartable,
        lavmodel = lavmodel, x = x
      )
      x <- out$x
      eqs[[b]] <- out$eqs
    }
  } # blocks

  # pooled (system) solve for shared free parameters among the directed
  # coefficients, across all blocks (a no-op when nothing is shared); this
  # handles simple (cross-group) equality constraints and general linear
  # equality constraints among the directed coefficients
  x <- lav_sem_miiv_apply_directed_pool(
    eqs_b = do.call(c, eqs), x = x, lavmodel = lavmodel,
    lavpartable = lavpartable
  )

  theta1 <- x[free_directed_idx]
  attr(theta1, "eqs") <- eqs
  theta1
}


# the current loading estimates as an nvar x nfac matrix (fixed values
# from the parameter table; free values from the parameter vector x)
lav_sem_js_lambda_mat <- function(x = NULL, lavpartable = NULL,
                                  mm_idx = NULL, ov_names = NULL,
                                  lv_names = NULL) {
  lambda_mat <- matrix(0, nrow = length(ov_names), ncol = length(lv_names))
  if (length(mm_idx) > 0L && length(lv_names) > 0L) {
    lambda_val <- lavpartable$ustart[mm_idx]
    lambda_val[is.na(lambda_val)] <- 0
    free_mm <- which(lavpartable$free[mm_idx] > 0L)
    lambda_val[free_mm] <- x[lavpartable$free[mm_idx][free_mm]]
    lambda_mat[cbind(
      match(lavpartable$rhs[mm_idx], ov_names),
      match(lavpartable$lhs[mm_idx], lv_names)
    )] <- lambda_val
  }
  lambda_mat
}


# one pass over (a subset of) the equations of a single block: for each
# equation, build the conditioning set Z from its plan (a weighted
# aggregate over the planned proxy variables per latent regressor; the
# observed regressors themselves), and solve the least-squares regression
# of the dependent variable on the estimated conditional expectations --
# all in moment form:
#
#   design  d = G S_zz^{-1} Z   with row l of G = Cov(eta_l, Z)
#   amat = Var(d)  = G S_zz^{-1} G'
#   bvec = Cov(d, y) = G S_zz^{-1} Cov(Z, y)
#   slopes = solve(amat, bvec)
#
# for an observed regressor, the corresponding row of G equals the row of
# S_zz, so its 'conditional expectation' is the variable itself
lav_sem_js_eqs_pass <- function(eqs_b = NULL, plans = NULL, only = NULL,
                                s_mat = NULL, m_vec = NULL, n = NULL,
                                ov_names = NULL, lv_names = NULL,
                                markers = NULL, contam = NULL,
                                theta = NULL, csmall = 1,
                                lambda_mat = NULL, lavpartable = NULL,
                                lavmodel = NULL, x = NULL) {
  meanstructure <- lavmodel@meanstructure
  nvar <- length(ov_names)
  if (is.null(only)) {
    only <- seq_along(eqs_b)
  }

  for (j in only) {
    eq <- eqs_b[[j]]
    y_idx <- match(eq$lhs_new, ov_names)

    # intercept-only equation
    if (identical(eq$rhs_new, "1")) {
      if (meanstructure && length(eq$ptint) == 1L) {
        free_int_idx <- lavpartable$free[eq$ptint]
        if (free_int_idx > 0L) {
          x[free_int_idx] <- m_vec[y_idx]
        }
      }
      eqs_b[[j]]$coef <- if (meanstructure) m_vec[y_idx] else numeric(0L)
      eqs_b[[j]]$nobs <- n
      next
    }

    nrhs <- length(eq$rhs)
    is_lv <- eq$rhs %in% lv_names

    # conditioning set: one row of u_mat (over the observed variables) per
    # regressor, plus its scale (kappa) and error variance (errvar); the
    # proxy variable sets come from the (structural) plan
    u_mat <- matrix(0, nrow = nrhs, ncol = nvar)
    kappa <- rep(1, nrhs)
    errvar <- rep(0, nrhs)
    for (i in seq_len(nrhs)) {
      if (!is_lv[i]) {
        u_mat[i, match(eq$rhs_new[i], ov_names)] <- 1
        next
      }
      f <- match(eq$rhs[i], lv_names)
      marker_idx <- match(markers[[eq$rhs[i]]], ov_names)
      set <- plans[[j]][[i]]
      if (length(set) == 1L) {
        w <- 1
      } else {
        w <- lav_sem_js_weights(
          s = s_mat[set, set, drop = FALSE],
          theta = theta[set]
        )
      }
      kap <- sum(w * lambda_mat[set, f])
      if (!is.finite(kap) || abs(kap) < .Machine$double.eps^0.5) {
        # degenerate proxy scale; fall back to the scaling indicator when
        # it is admissible for this equation
        if (!contam[marker_idx, y_idx] && !is.na(theta[marker_idx]) &&
            !identical(set, marker_idx)) {
          set <- marker_idx
          w <- 1
          kap <- lambda_mat[marker_idx, f]
        }
        if (!is.finite(kap) || abs(kap) < .Machine$double.eps^0.5) {
          lav_msg_stop(gettextf(
            "the conditioning proxy for %1$s in the equation for %2$s has a
             degenerate scale (aggregate loading near zero).",
            eq$rhs[i], eq$lhs_new))
        }
      }
      u_mat[i, set] <- w
      kappa[i] <- kap
      errvar[i] <- sum(w * w * theta[set])
    }

    # moments of the conditioning variables
    s_zz <- u_mat %*% s_mat %*% t(u_mat)
    s_zy <- drop(u_mat %*% s_mat[, y_idx])

    # G: covariance of the regressors (latent or observed) with Z
    g_mat <- s_zz
    for (i in which(is_lv)) {
      g_mat[i, i] <- g_mat[i, i] - csmall * errvar[i]
    }
    if (any(is_lv)) {
      g_mat[is_lv, ] <- g_mat[is_lv, , drop = FALSE] / kappa[is_lv]
    }

    # normal equations of the regression on the conditional expectations
    s_zz_inv <- try(solve((s_zz + t(s_zz)) / 2), silent = TRUE)
    if (inherits(s_zz_inv, "try-error")) {
      lav_msg_stop(gettextf(
        "the conditioning variables of the equation for %s are collinear.",
        eq$lhs_new))
    }
    gsi <- g_mat %*% s_zz_inv
    amat <- gsi %*% t(g_mat)
    amat <- (amat + t(amat)) / 2
    bvec <- drop(gsi %*% s_zy)

    # means of the design (for the intercepts)
    if (meanstructure) {
      z_mean <- drop(u_mat %*% m_vec)
      x_bar <- z_mean
      x_bar[is_lv] <- z_mean[is_lv] / kappa[is_lv]
      y_bar <- m_vec[y_idx]
    } else {
      x_bar <- numeric(nrhs)
      y_bar <- 0
    }

    # free/fixed partition and solve
    eq_free_idx <- lavpartable$free[eq$pt]
    b_full <- lavpartable$ustart[eq$pt]
    b_full[is.na(b_full)] <- 0
    fix_idx <- which(eq_free_idx == 0L)
    if (length(fix_idx) > 0L) {
      amat_free <- amat[-fix_idx, -fix_idx, drop = FALSE]
      bvec_free <- drop(bvec[-fix_idx] -
        amat[-fix_idx, fix_idx, drop = FALSE] %*% b_full[fix_idx])
      gcol_free <- eq_free_idx[-fix_idx]
    } else {
      amat_free <- amat
      bvec_free <- bvec
      gcol_free <- eq_free_idx
    }
    if (length(gcol_free) > 0L) {
      slopes_free <- drop(solve(amat_free, bvec_free))
      x[gcol_free] <- slopes_free
      b_full[eq_free_idx > 0L] <- slopes_free
      # slope system for the pooled solve (simple equality constraints)
      eqs_b[[j]]$slope_block <- list(
        amat = (n - 1) * amat_free, bvec = (n - 1) * bvec_free,
        gcol = gcol_free, x_idx = integer(0L),
        x_bar = x_bar, y_bar = y_bar
      )
    }

    # intercept
    beta0 <- y_bar - sum(b_full * x_bar)
    if (meanstructure && length(eq$ptint) == 1L) {
      free_int_idx <- lavpartable$free[eq$ptint]
      if (free_int_idx > 0L) {
        x[free_int_idx] <- beta0
      }
    }

    eqs_b[[j]]$coef <- c(if (meanstructure) beta0 else NULL, b_full)
    eqs_b[[j]]$nobs <- n
  }

  list(x = x, eqs = eqs_b)
}
