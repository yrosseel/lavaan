# Standardize parameter estimates of a fitted lavaan model.
#
# Three "modes" of standardization are supported, all implemented on top of a
# single observed-variable engine (lav_standardize_all):
#   - std.lv  : only the latent variables are standardized
#   - std.all : latent *and* observed variables are standardized
#   - std.nox : like std.all, but the (observed) exogenous 'x' variables are
#               left unstandardized
# std.nox is in fact a special case of a more general mechanism: standardize
# all observed variables *except* a given set (here: the exogenous 'x'). The
# 'ov_std' argument of lav_standardize_all() exposes this generality, so that a
# caller can request that only a user-specified subset of observed variables is
# standardized.

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

# Resolve the standard set of ingredients (lavmodel, lavpartable, partable,
# est, glist, cov_x) from either a fitted lavobject or the individual pieces.
lav_standardize_setup <- function(lavobject = NULL, lavmodel = NULL,
                                  lavpartable = NULL, partable = NULL,
                                  est = NULL, glist = NULL, cov_x = NULL,
                                  need_cov_x = FALSE) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel), !is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est
                # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_inspect_est(lavobject)
    }
    if (need_cov_x && lavmodel@conditional.x && is.null(cov_x)) {
      if (!is.null(lavobject@implied$cov.x[[1]])) {
        cov_x <- lavobject@implied$cov.x
                # if this changes, tag @TDJorgensen in commit message
      } else {
        # perhaps a lavaanList object: extract it from GLIST per block
        cov_x <- vector("list", length = lavmodel@nblocks)
        for (b in seq_len(lavmodel@nblocks)) {
          mm_in_block <- (seq_len(lavmodel@nmat[b]) +
            cumsum(c(0, lavmodel@nmat))[b])
          mlist <- lavmodel@GLIST[mm_in_block]
          cov_x[[b]] <- mlist[["cov.x"]]
        }
      }
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }

  list(lavmodel = lavmodel, lavpartable = lavpartable, partable = partable,
       est = est, glist = glist, cov_x = cov_x)
}

# Evaluate the standardized values of the defined (:=), equality (==) and
# inequality (< >) rows, given the (already standardized) 'out' vector.
lav_standardize_constraints <- function(out, partable, lavmodel) {
  x <- out[partable$free & !duplicated(partable$free)]

  idx <- which(partable$op == ":=")
  if (length(idx) > 0L) {
    out[idx] <- lavmodel@def.function(x)
  }
  idx <- which(partable$op == "==")
  if (length(idx) > 0L) {
    out[idx] <- lavmodel@ceq.function(x)
  }
  idx <- which(partable$op == "<" | partable$op == ">")
  if (length(idx) > 0L) {
    cin <- lavmodel@cin.function(x)
    # when the model also has box bounds, the inequality function holds
    # the bound rows (first) as well; only the explicit rows have a row
    # in the parameter table
    bound_idx <- attr(cin, "bound.idx", exact = TRUE)
    if (length(bound_idx) > 0L) {
      cin <- cin[-bound_idx]
    }
    if (length(cin) == length(idx)) {
      out[idx] <- cin
    }
  }

  out
}

# Observed-variable standard deviations (and variances) for block 'g'.
# Returns the (possibly conditional.x-extended) ov.names together with the
# vectors of SDs ('ov') and variances ('ov2').
lav_standardize_ov_sd <- function(lavpartable, g, vy = NULL, ov_var = NULL,
                                  conditional_x = FALSE, cov_x = NULL) {
  ov_names <- lav_pt_vnames(lavpartable, "ov", block = g) # not user

  if (is.null(ov_var)) {
    ov2 <- vy[[g]]
    # replace zero values by NA (but keep negative values)
    zero_idx <- which(abs(ov2) < .Machine$double.eps)
    if (length(zero_idx) > 0L) {
      ov2[zero_idx] <- as.numeric(NA)
    }
    # replace negative values by NA (for sqrt)
    tmp_ov2 <- ov2
    neg_idx <- which(tmp_ov2 < 0)
    if (length(neg_idx) > 0L) {
      tmp_ov2[neg_idx] <- as.numeric(NA)
    }
    ov <- sqrt(tmp_ov2)
  } else {
    ov2 <- ov_var[[g]]
    ov <- sqrt(ov2)
  }

  if (conditional_x) {
    # extend OV with ov.names.x
    ov_names_x <- lav_pt_vnames(lavpartable, "ov.x", block = g)
    ov_names_nox <- lav_pt_vnames(lavpartable, "ov.nox", block = g)
    ov_names <- c(ov_names_nox, ov_names_x)
    ov2 <- c(ov2, diag(cov_x[[g]]))
    ov <- c(ov, sqrt(diag(cov_x[[g]])))
  }

  list(ov_names = ov_names, ov = ov, ov2 = ov2)
}

# scaling corrections for two-way product terms (interactions X:Z, quadratic
# terms X:X) in block 'g', based on Brandt et al. (2015, Eq. 13) and Kelava &
# Brandt (2022). The standardized solution treats a product term as the
# product of the centered, standardized components, whereas the model
# parameters describe the raw product A:B. Three corrections are needed, all
# exact (no distributional assumptions):

# (1) scaling SD (FV): the SD of a standardized product term is
# SD(A)*SD(B), not SD(A:B). var_sd[A:B] is replaced accordingly, so that the
# generic scaling in the caller applies the right factor everywhere A:B
# appears.

# (2) simple main effects (KS, 08/05/26): with uncentered components, b1 in
# y ~ b1*x + b2*z + b3*x:z + b4*x:x is the slope at x = z = 0. The
# standardized solution reports the slope at the means, i.e. the partial
# derivative dy/dx = b1 + b3*z + 2*b4*x evaluated at E(x), E(z). Looping
# over the product terms, each term adds beta_ab * E(other component) to
# the main effect of each of its components.

# (3) product moments (FV, 05/07/26): the "~~"/"~1" rows of A:B describe the
# uncentered product. Without centering them first, the standardized
# moments would depend on how the mean structure is parameterized. With
# ac = A - E(A) and bc = B - E(B) we have the exact identity
# ac*bc = A:B - E(B)*A - E(A)*B + E(A)*E(B), so the centered product is
# a linear combination of the variables. Collect the weights in a
# matrix K (one row per variable, regular variables keep their identity
# row), e.g. for (X, Z, X:Z, X:X) with mx = E(X), mz = E(Z):
# row X:Z = (-mz, -mx, 1, 0), row X:X = (-2*mx, 0, 0, 1). Then
# K %*% C %*% t(K) holds all centered (co)variances at once (constants
# do not affect them), and the centered means follow from
# E(ac*bc) = E(A:B) - E(A)*E(B). C holds the "~~" parameters (PSI for
# latent variables): for a residual covariance with an endogenous
# variable Y, PSI[X, Y] = Cov(X, zeta_Y) is exactly what the centering
# needs. The corrected rows are also written into est_c, because the
# cov_std block in the caller derives its sqrt(variance) divisors from
# est_c (not out, which is scaled by then).

# The variables are described by parallel vectors: var_names, their scaling
# SDs (var_sd), their means (var_mean) and, optionally, var_scale: the
# factor that maps the model metric of a variable to its raw metric
# (the ~*~ scale under the Delta-augmented ML correlation structure, 1
# otherwise). C is expressed in the model metric; the centering identity
# holds in the raw metric, so with var_scale the (co)variances are
# converted to the raw metric, centered, and converted back.

# A component may be a latent or an observed variable; a product term is
# processed when it is two-way and both components are in var_names.
# Returns the corrected (out, est_c, var_sd)
lav_standardize_int_terms <- function(out, est_c, var_sd, est, partable,
                                      var_names, int_names, var_mean,
                                      cov_mat, g, var_scale = NULL) {
  # component indices of each usable product term (int = A:B, a = A, b = B)
  int_names_ok <- character(0L)
  int_comp <- list()
  for (int_name in int_names) {
    components <- strsplit(int_name, ":", fixed = TRUE)[[1L]]
    # only two-way terms: a three-way (or higher) product would need its own
    # decomposition, and the product of three SDs
    if (length(components) != 2L) next
    idx_int <- match(int_name, var_names)
    idx_a   <- match(components[1L], var_names)
    idx_b   <- match(components[2L], var_names)
    # need a, b and the interaction term before moving on
    if (is.na(idx_a) || is.na(idx_b) || is.na(idx_int)) next
    int_names_ok <- c(int_names_ok, int_name)
    int_comp <- c(int_comp,
      list(list(int = idx_int, a = idx_a, b = idx_b)))
  }
  if (length(int_comp) == 0L) {
    return(list(out = out, est_c = est_c, var_sd = var_sd))
  }

  # only used by the shift in (2), which needs nonzero means, but guarding
  # it with any(var_mean != 0) would cost as much as building it
  reg_row <- partable$op == "~" & partable$block == g
  for (i in seq_along(int_comp)) {
    p <- int_comp[[i]]

    # (1) replace the SD of the product term by SD(A)*SD(B)
    var_sd[p$int] <- var_sd[p$a] * var_sd[p$b]

    # (2) shift the simple main effects to the mean-evaluated slopes
    exp_a <- var_mean[p$a]
    exp_b <- var_mean[p$b]
    if (exp_a != 0 || exp_b != 0) {
      int_row <- reg_row & partable$rhs == int_names_ok[i]
      for (dep in unique(partable$lhs[int_row])) {
        # beta_ab is read from est (not partable$est): in the *_x wrappers
        # (delta method, bootstrap) est is the perturbed parameter vector, and
        # the gradient of the shifted main effects wrt beta_ab must not be lost
        beta_ab <- est[int_row & partable$lhs == dep]
        dep_row <- reg_row & partable$lhs == dep
        idx_beta_a <- which(dep_row & partable$rhs == var_names[p$a])
        idx_beta_b <- which(dep_row & partable$rhs == var_names[p$b])
        if (is.null(var_scale)) {
          out[idx_beta_a] <- out[idx_beta_a] + beta_ab * exp_b
          out[idx_beta_b] <- out[idx_beta_b] + beta_ab * exp_a
        } else {
          # the main effect of A is per unit of A in the model metric
          # (= raw A / var_scale[A]); the product is in the raw metric
          out[idx_beta_a] <- out[idx_beta_a] + beta_ab * exp_b * var_scale[p$a]
          out[idx_beta_b] <- out[idx_beta_b] + beta_ab * exp_a * var_scale[p$b]
        }
      }
    }
  }

  # (3) center the (co)variances and means of the product terms
  if (!is.null(cov_mat) && any(var_mean != 0)) {
    kmat <- diag(nrow(cov_mat))
    for (p in int_comp) {
      kmat[p$int, p$a] <- kmat[p$int, p$a] - var_mean[p$b]
      kmat[p$int, p$b] <- kmat[p$int, p$b] - var_mean[p$a]
    }
    # shift_mat holds, for every pair of variables, the mean-driven part
    # of their (co)variance (zero unless a product term is involved)
    if (is.null(var_scale)) {
      shift_mat <- cov_mat - kmat %*% cov_mat %*% t(kmat)
    } else {
      scale_mat <- tcrossprod(var_scale)
      cov_raw <- cov_mat * scale_mat
      shift_mat <- cov_mat - (kmat %*% cov_raw %*% t(kmat)) / scale_mat
    }

    # "~~" rows involving at least one product term: this covers the
    # variance of a product term (X:Z ~~ X:Z), its covariance with a
    # regular variable (X:Z ~~ V), and product-product covariances
    # (X:Z ~~ X:X). Rows whose other side is not in var_names are left
    # alone (keep)
    idx_cov <- which(partable$op == "~~" & partable$block == g &
      (partable$lhs %in% int_names_ok |
       partable$rhs %in% int_names_ok))
    l_idx <- match(partable$lhs[idx_cov], var_names)
    r_idx <- match(partable$rhs[idx_cov], var_names)
    keep <- which(!is.na(l_idx) & !is.na(r_idx))
    shift <- shift_mat[cbind(l_idx[keep], r_idx[keep])]
    out[idx_cov[keep]] <- out[idx_cov[keep]] - shift
    # of est_c only the variance entries are read back later, but centering
    # all touched rows keeps est_c fully consistent and a diagonal-only
    # filter would cost as much as it saves
    est_c[idx_cov[keep]] <- est_c[idx_cov[keep]] - shift

    # "~1" rows of the product terms: E(ac*bc) = E(A:B) - E(A)*E(B). After
    # the generic scaling this equals, e.g., cor(X, Z) for X:Z
    idx_mean <- which(partable$op == "~1" & partable$block == g &
      partable$lhs %in% int_names_ok)
    for (row_i in idx_mean) {
      p <- int_comp[[match(partable$lhs[row_i], int_names_ok)]]
      out[row_i] <- out[row_i] - var_mean[p$a] * var_mean[p$b]
    }
  }

  list(out = out, est_c = est_c, var_sd = var_sd)
}

# How are product terms with an observed component (x1:x2, f1:x1)
# standardized? "product": the product term is a variable in its own right
# (scaled by its own SD, no centering). "component": like latent product
# terms (see lav_standardize_int_terms). Controlled by the option
# std.ov.interaction; objects without it (older fits, no object) use
# "product".
lav_standardize_ov_int_mode <- function(lavobject = NULL) {
  mode <- NULL
  if (!is.null(lavobject) && isS4(lavobject) && .hasSlot(lavobject, "Options")) {
    mode <- lavobject@Options$std.ov.interaction
  }
  if (is.null(mode)) {
    mode <- "product"
  }
  mode
}

# TRUE if the model has at least one product term with an observed
# component (an observed product x1:x2, or a mixed product f1:x1)
lav_standardize_ov_int_any <- function(lavpartable) {
  int_names <- unique(unlist(c(
    lav_pt_vnames(lavpartable, "lv.interaction"),
    lav_pt_vnames(lavpartable, "ov.interaction")
  )))
  if (length(int_names) == 0L) {
    return(FALSE)
  }
  ov_names <- unique(unlist(lav_pt_vnames(lavpartable, "ov")))
  comps <- unique(unlist(strsplit(int_names, ":", fixed = TRUE)))
  any(comps %in% ov_names)
}

# the "~~" parameters of block g as a symmetric matrix over var_names
# (zero where the model has no such parameter)
lav_standardize_cov_from_partable <- function(partable, est, var_names, g) {
  n <- length(var_names)
  cov_mat <- matrix(0, n, n)
  idx <- which(partable$op == "~~" & partable$block == g &
    partable$lhs %in% var_names & partable$rhs %in% var_names)
  if (length(idx) > 0L) {
    l_idx <- match(partable$lhs[idx], var_names)
    r_idx <- match(partable$rhs[idx], var_names)
    cov_mat[cbind(l_idx, r_idx)] <- est[idx]
    cov_mat[cbind(r_idx, l_idx)] <- est[idx]
  }
  cov_mat
}

# per-block description of the observed variables for the product-term
# corrections: names, scaling SDs, means, and the model-to-raw scale
# (the ~*~ parameter, 1 if absent). This 'unit' version leaves the
# observed variables in their own metric (std.lv): SD 1, mean 0
lav_standardize_ov_info_unit <- function(lavpartable, est, nblocks) {
  ov_info <- vector("list", length = nblocks)
  for (g in seq_len(nblocks)) {
    ov_names <- unique(unlist(lav_pt_vnames(lavpartable, "ov", block = g)))
    scale_g <- rep(1, length(ov_names))
    idx <- which(lavpartable$op == "~*~" & lavpartable$block == g)
    if (length(idx) > 0L && length(est) == length(lavpartable$lhs)) {
      m_idx <- match(lavpartable$lhs[idx], ov_names)
      ok <- !is.na(m_idx)
      scale_g[m_idx[ok]] <- est[idx[ok]]
    }
    ov_info[[g]] <- list(names = ov_names, sd = scale_g,
      mean = rep(0, length(ov_names)), scale = scale_g)
  }
  ov_info
}

# model-implied means of the observed variables per block (named vectors),
# in the raw metric; NULL entries when they cannot be computed
lav_standardize_ov_mean <- function(lavobject = NULL, lavmodel = NULL,
                                    lavpartable = NULL, glist = NULL) {
  nblocks <- lavmodel@nblocks
  res <- vector("list", length = nblocks)
  if (is.null(lavobject) || !isS4(lavobject) ||
      !.hasSlot(lavobject, "SampleStats") ||
      lavmodel@representation != "LISREL") {
    return(res)
  }
  lavsamplestats <- lavobject@SampleStats
  if (lavmodel@meanstructure) {
    ey <- lav_model_ey(lavmodel = lavmodel, glist = glist,
      lavsamplestats = lavsamplestats, delta = TRUE)
  } else {
    # no mean structure: the (raw) sample means. lav_model_ey() would
    # reconstruct them too, but then rescale them by the ~*~ parameters
    # under the Delta-augmented ML correlation structure (they are
    # already in the raw metric); conditional.x always has a mean
    # structure, so the block/group layout is simple here
    ey <- lavsamplestats@mean
  }
  for (g in seq_len(nblocks)) {
    if (length(ey) < g || is.null(ey[[g]])) next
    if (lavmodel@conditional.x) {
      ov_names_nox <- lav_pt_vnames(lavpartable, "ov.nox", block = g)
      ov_names_x <- lav_pt_vnames(lavpartable, "ov.x", block = g)
      mean_x <- lavsamplestats@mean.x[[g]]
      if (length(ey[[g]]) != length(ov_names_nox) ||
          length(mean_x) != length(ov_names_x)) next
      mean_g <- c(ey[[g]], mean_x)
      names(mean_g) <- c(ov_names_nox, ov_names_x)
    } else {
      ov_names <- lav_pt_vnames(lavpartable, "ov", block = g)
      if (length(ey[[g]]) != length(ov_names)) next
      mean_g <- ey[[g]]
      names(mean_g) <- ov_names
    }
    res[[g]] <- mean_g
  }
  res
}

# Given a parameter vector 'x', set it as the model parameters and return the
# (glist, est) pair to feed into the standardization engines. Handles the
# (EFA) rotation case. Used by the bootstrap/Monte-Carlo machinery.
lav_standardize_prepare_x <- function(x, lavobject, rotation = FALSE) {
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est_rot <- lav_model_efa_rotate_x(
      x = x,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init_rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    glist <- attr(est_rot, "extra")$glist
    attributes(est_rot) <- NULL
    est <- est_rot
  } else {
    glist <- lavmodel@GLIST
                    # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  list(glist = glist, est = est)
}

# ---------------------------------------------------------------------------
# 'x' wrappers (parameter vector -> standardized values), used for bootstrap
# and Monte-Carlo standard errors.
# ---------------------------------------------------------------------------

lav_standardize_lv_x <- function(x, lavobject, partable = NULL, cov_std = TRUE,
                                 lv_var = NULL, rotation = FALSE) {
  ge <- lav_standardize_prepare_x(x, lavobject, rotation = rotation)
  lav_standardize_lv(
    lavobject = lavobject, partable = partable, est = ge$est,
    glist = ge$glist, cov_std = cov_std, lv_var = lv_var
  )
}

lav_standardize_all_x <- function(x, lavobject, partable = NULL, cov_std = TRUE,
                                  ov_std = NULL, nox = FALSE, rotation = FALSE) {
  ge <- lav_standardize_prepare_x(x, lavobject, rotation = rotation)
  lav_standardize_all(
    lavobject = lavobject, partable = partable, est = ge$est,
    est_std = NULL, glist = ge$glist, cov_std = cov_std,
    ov_std = ov_std, nox = nox
  )
}

lav_standardize_all_nox_x <- function(x, lavobject, partable = NULL,
                                      cov_std = TRUE, rotation = FALSE) {
  lav_standardize_all_x(
    x = x, lavobject = lavobject, partable = partable, cov_std = cov_std,
    nox = TRUE, rotation = rotation
  )
}

# ---------------------------------------------------------------------------
# std.lv : standardize the latent variables only
# ---------------------------------------------------------------------------

lav_standardize_lv <- function(lavobject = NULL,
                               partable = NULL, est = NULL, glist = NULL,
                               cov_std = TRUE, lv_var = NULL,
                               lavmodel = NULL, lavpartable = NULL,
                               ov_info = NULL) {
  s <- lav_standardize_setup(
    lavobject = lavobject, lavmodel = lavmodel, lavpartable = lavpartable,
    partable = partable, est = est, glist = glist
  )
  lavmodel <- s$lavmodel
  lavpartable <- s$lavpartable
  partable <- s$partable
  est <- s$est
  glist <- s$glist

  out <- est
  n <- length(est)
  stopifnot(n == length(partable$lhs))

  # working copy of est: for interaction terms (e.g., X:Z), the (co)variance
  # and mean rows are 'centered' below (the mean-driven part of the product
  # moments is removed) before any scaling takes place. The cov_std block
  # further down needs these centered -- but still unscaled -- variances to
  # compute its sqrt(variance) divisors, so we keep them in est_c.
  est_c <- est

  # product terms with an observed component (x1:x2, f1:x1): "product"
  # (default: the product term is a variable in its own right) or
  # "component" (treated like latent product terms, option
  # std.ov.interaction). In component mode, the corrections of
  # lav_standardize_int_terms() need the SDs and means of the observed
  # components: 'ov_info' (per block: names, sd, mean, scale) is supplied
  # by lav_standardize_all(); a standalone (std.lv) call leaves the
  # observed variables in their own metric (SD 1, mean 0)
  internal_call <- !is.null(ov_info)
  component_mode <- lav_standardize_ov_int_mode(lavobject) == "component" &&
    lav_standardize_ov_int_any(lavpartable)
  if (component_mode && is.null(ov_info)) {
    ov_info <- lav_standardize_ov_info_unit(lavpartable, est,
      lavmodel@nblocks)
  }
  # the "~~" parameters used for the centering must be the MODEL's: when
  # the caller hands in another table (e.g., the modification indices with
  # est = epc), fall back to the estimates stored in lavpartable
  if (component_mode) {
    if (identical(partable$lhs, lavpartable$lhs) &&
        identical(partable$op, lavpartable$op) &&
        identical(partable$rhs, lavpartable$rhs) &&
        length(est) == length(lavpartable$lhs)) {
      cov_est <- est
    } else if (!is.null(lavpartable$est)) {
      cov_est <- lavpartable$est
    } else {
      cov_est <- est
    }
    ov_int_sd <- vector("list", length = lavmodel@nblocks)
  }

  # any latent product terms (in any block)? only their corrections (see
  # lav_standardize_int_terms) need E(ETA), psi and the mm indices. Regular
  # models skip these, and E(ETA) is LISREL-only anyway (RAM would stop)
  lv_int_any <- length(lav_pt_vnames(lavpartable, "lv.interaction")) > 0L
  if (lv_int_any) {
    mm_idx <- lav_model_group_mm_indices(lavmodel@nmat)
  }

  # compute ETA
  if (is.null(lv_var)) {
    lv_eta <- lav_model_veta(
      lavmodel = lavmodel,
      glist = glist
    )

    if (lv_int_any) {
      lv_eeta <- lav_model_eeta(
        lavmodel = lavmodel,
        glist = glist,
        lavsamplestats = lavobject@SampleStats
      )
    }
  }

  for (g in 1:lavmodel@nblocks) {
    ov_names <- lav_pt_vnames(lavpartable, "ov", block = g) # not user,
    # which may be incomplete
    lv_names <- lav_pt_vnames(lavpartable, "lv", block = g)

    # observed product terms in this block (component mode only)
    ov_int_names <- character(0L)
    if (component_mode) {
      ov_int_names <- lav_pt_vnames(lavpartable, "ov.interaction", block = g)
    }

    # shortcut: no latents in this block, nothing to do
    if (length(lv_names) == 0L && length(ov_int_names) == 0L) {
      next
    }

    if (length(lv_names) == 0L) {
      # observed product terms only (component mode)
      veta_g <- NULL
      eta2 <- numeric(0L)
    } else if (is.null(lv_var)) {
      veta_g <- lv_eta[[g]] # full model-implied V(ETA), only kept as a
                            # psi fallback for the centering below (RAM)
      eta2 <- diag(veta_g)
    } else {
      veta_g <- NULL
      eta2 <- lv_var[[g]]
    }
    # change negative values to NA
    eta2[eta2 < 0] <- as.numeric(NA)
    eta <- sqrt(eta2)

    # two-way product terms (X:Z, X:X): shift the simple main effects,
    # scale by SD(A)*SD(B), and center the product moments (see
    # lav_standardize_int_terms() for the rationale)
    if (lv_int_any) {
      lv_int_names <- lav_pt_vnames(lavpartable, "lv.interaction",
        block = g
      )
    } else {
      lv_int_names <- character(0L)
    }
    if (length(lv_int_names) > 0L || length(ov_int_names) > 0L) {
      if (length(lv_int_names) > 0L && is.null(lv_var)) {
        eeta <- lv_eeta[[g]]
      } else {
        eeta <- numeric(length(eta2))
      }
    }
    if (component_mode &&
        (length(lv_int_names) > 0L || length(ov_int_names) > 0L)) {
      # latent, mixed and observed product terms together, over the
      # combined (latent + observed) variable set of this block
      nlv <- length(lv_names)
      ovi <- ov_info[[g]]
      var_names <- c(lv_names, ovi$names)
      var_sd <- c(eta[seq_len(nlv)], ovi$sd)
      var_mean <- c(eeta[seq_len(nlv)], ovi$mean)
      if (all(ovi$scale == 1)) {
        var_scale <- NULL
      } else {
        var_scale <- c(rep(1, nlv), ovi$scale)
      }
      cov_mat <- lav_standardize_cov_from_partable(
        partable = lavpartable, est = cov_est, var_names = var_names, g = g
      )
      tmp <- lav_standardize_int_terms(
        out = out, est_c = est_c, var_sd = var_sd, est = est,
        partable = partable, var_names = var_names,
        int_names = c(lv_int_names, ov_int_names), var_mean = var_mean,
        cov_mat = cov_mat, g = g, var_scale = var_scale
      )
      out <- tmp$out
      est_c <- tmp$est_c
      eta[seq_len(nlv)] <- tmp$var_sd[seq_len(nlv)]
      # the scaling SDs of the observed product terms, for the ov engine
      if (length(ov_int_names) > 0L) {
        int_sd <- tmp$var_sd[nlv + match(ov_int_names, ovi$names)]
        names(int_sd) <- ov_int_names
        ov_int_sd[[g]] <- int_sd
      }
    } else if (length(lv_int_names) > 0L) {
      # latent product terms with latent components only
      psi_g <- glist[mm_idx[[g]]]$psi
      if (is.null(psi_g)) { # e.g., RAM representation
        psi_g <- veta_g
      }
      tmp <- lav_standardize_int_terms(
        out = out, est_c = est_c, var_sd = eta, est = est,
        partable = partable, var_names = lv_names, int_names = lv_int_names,
        var_mean = eeta, cov_mat = psi_g, g = g
      )
      out <- tmp$out
      est_c <- tmp$est_c
      eta <- tmp$var_sd
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * eta[match(partable$lhs[idx], lv_names)]

    # 1b. "=~" regular higher-order lv indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% ov_names) &
      partable$block == g)
    out[idx] <- (out[idx] * eta[match(partable$lhs[idx], lv_names)]
      / eta[match(partable$rhs[idx], lv_names)])

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] / eta[match(partable$lhs[idx], lv_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] * eta[match(partable$rhs[idx], lv_names)]

    # 3b. "~~" lv variances
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of ETA
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'PSI' diagonal elements!!
    rv_idx <- which(partable$op == "~~" & partable$rhs %in% lv_names &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] / eta[match(partable$lhs[rv_idx], lv_names)]
      / eta[match(partable$rhs[rv_idx], lv_names)])

    # covariances lv
    # three types:
    # - only lhs is LV (and fixed.x = FALSE)
    # - only rhs is LV (and fixed.x = FALSE)
    # - both lhs and rhs are LV (regular case)
    if (cov_std) {
      # est_c: est with the interaction-term variances centered (see above)
      if (!is.complex(est_c[rv_idx])) {
        rv <- sqrt(abs(est_c[rv_idx])) # abs in case of heywood cases
      } else {
        rv <- sqrt(est_c[rv_idx])
      }
      rv_names <- partable$lhs[rv_idx]
    }

    # left
    idx_lhs <- which(partable$op == "~~" &
      partable$lhs %in% lv_names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_lhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_lhs] <-
          (out[idx_lhs] / eta[match(partable$lhs[idx_lhs], lv_names)])
      } else {
        out[idx_lhs] <-
          (out[idx_lhs] / rv[match(partable$lhs[idx_lhs], rv_names)])
      }
    }

    # right
    idx_rhs <- which(partable$op == "~~" &
      partable$rhs %in% lv_names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_rhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_rhs] <-
          (out[idx_rhs] / eta[match(partable$rhs[idx_rhs], lv_names)])
      } else {
        out[idx_rhs] <-
          (out[idx_rhs] / rv[match(partable$rhs[idx_rhs], rv_names)])
      }
    }

    # 4b. "~1" lv
    idx <- which(partable$op == "~1" & partable$lhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] / eta[match(partable$lhs[idx], lv_names)]
  }

  # 5. defined/(in)equality constraints
  out <- lav_standardize_constraints(out, partable, lavmodel)

  # internal call from lav_standardize_all(): hand back the centered est
  # and the scaling SDs of the observed product terms
  if (internal_call) {
    attr(out, "est_c") <- est_c
    attr(out, "ov_int_sd") <- if (component_mode) ov_int_sd else NULL
  }

  out
}

# ---------------------------------------------------------------------------
# std.all / std.nox / user subset : standardize observed (and latent)
# variables.
#
# 'ov_std' and 'nox' control *which* observed variables are standardized:
#   - nox = FALSE, ov_std = NULL      -> all observed variables (std.all)
#   - nox = TRUE                      -> all but the exogenous 'x' (std.nox)
#   - nox = FALSE, ov_std = c(...)    -> only the observed variables listed in
#                                        'ov_std' (a generalization of std.nox)
# In all cases the latent variables are fully standardized (via
# lav_standardize_lv()).
# ---------------------------------------------------------------------------

lav_standardize_all <- function(lavobject = NULL,
                                partable = NULL, est = NULL, est_std = NULL,
                                glist = NULL, cov_std = TRUE, ov_var = NULL,
                                lv_var = NULL,
                                lavmodel = NULL, lavpartable = NULL,
                                cov_x = NULL, ov_std = NULL, nox = FALSE) {
  s <- lav_standardize_setup(
    lavobject = lavobject, lavmodel = lavmodel, lavpartable = lavpartable,
    partable = partable, est = est, glist = glist, cov_x = cov_x,
    need_cov_x = TRUE
  )
  lavmodel <- s$lavmodel
  lavpartable <- s$lavpartable
  partable <- s$partable
  est <- s$est
  glist <- s$glist
  cov_x <- s$cov_x

  vy <- lav_model_vy(
    lavmodel = lavmodel, glist = glist,
    diagonal_only = TRUE
  )

  # which observed variables should NOT be standardized in block g?
  # (std.all: none; std.nox: the exogenous 'x'; user: all but 'ov_std')
  ov_nostd_block <- function(g, ov_names) {
    if (nox) {
      if (lavmodel@fixed.x) {
        lav_pt_vnames(lavpartable, "ov.x", block = g)
      } else {
        # if fixed.x = FALSE, the exogenous 'x' are treated as random
        # variables like any other, and std.nox reduces to std.all
        character(0L)
      }
    } else if (is.null(ov_std)) {
      character(0L)
    } else {
      setdiff(ov_names, ov_std)
    }
  }

  # product terms with an observed component in "component" mode (option
  # std.ov.interaction): the corrections are applied by
  # lav_standardize_lv(), which needs the SDs and means of the observed
  # components (in the metric requested here: a non-standardized variable
  # keeps SD 1 and mean 0). The std.lv result handed in by the caller
  # (est_std) was computed without that information, so it is recomputed
  component_mode <- lav_standardize_ov_int_mode(lavobject) == "component" &&
    lav_standardize_ov_int_any(lavpartable)
  est_c <- est
  ov_int_sd <- NULL
  if (component_mode) {
    ov_mean <- lav_standardize_ov_mean(
      lavobject = lavobject, lavmodel = lavmodel, lavpartable = lavpartable,
      glist = glist
    )
    ov_info <- vector("list", length = lavmodel@nblocks)
    for (g in 1:lavmodel@nblocks) {
      tmp <- lav_standardize_ov_sd(
        lavpartable = lavpartable, g = g, vy = vy, ov_var = ov_var,
        conditional_x = lavmodel@conditional.x, cov_x = cov_x
      )
      ov_names_g <- tmp$ov_names
      n_ov <- length(ov_names_g)
      # model-to-raw scale: the ~*~ parameter (Delta-augmented ML
      # correlation structure), 1 otherwise
      scale_g <- rep(1, n_ov)
      idx <- which(partable$op == "~*~" & partable$block == g)
      if (length(idx) > 0L) {
        m_idx <- match(partable$lhs[idx], ov_names_g)
        ok <- !is.na(m_idx)
        scale_g[m_idx[ok]] <- est[idx[ok]]
      }
      sd_g <- tmp$ov * scale_g
      mean_g <- rep(0, n_ov)
      if (!is.null(ov_mean[[g]])) {
        m_idx <- match(ov_names_g, names(ov_mean[[g]]))
        ok <- !is.na(m_idx)
        mean_g[ok] <- ov_mean[[g]][m_idx[ok]]
      }
      # variables that keep their metric: SD = scale (model -> raw),
      # mean 0 (no centering)
      nostd <- ov_names_g %in% ov_nostd_block(g, ov_names_g)
      sd_g[nostd] <- scale_g[nostd]
      mean_g[nostd] <- 0
      ov_info[[g]] <- list(names = ov_names_g, sd = sd_g, mean = mean_g,
        scale = scale_g)
    }
    est_std <- lav_standardize_lv(
      lavobject = lavobject,
      partable = partable, est = est, glist = glist,
      cov_std = cov_std, lv_var = lv_var, lavmodel = lavmodel,
      lavpartable = lavpartable, ov_info = ov_info
    )
    est_c <- attr(est_std, "est_c")
    ov_int_sd <- attr(est_std, "ov_int_sd")
    attributes(est_std) <- NULL
  } else if (is.null(est_std)) {
    est_std <- lav_standardize_lv(
      lavobject = lavobject,
      partable = partable, est = est, glist = glist,
      cov_std = cov_std, lv_var = lv_var, lavmodel = lavmodel,
      lavpartable = lavpartable
    )
  }

  out <- est_std
  n <- length(est_std)
  stopifnot(n == length(partable$lhs))

  for (g in 1:lavmodel@nblocks) {
    lv_names <- lav_pt_vnames(lavpartable, "lv", block = g)

    tmp <- lav_standardize_ov_sd(
      lavpartable = lavpartable, g = g, vy = vy, ov_var = ov_var,
      conditional_x = lavmodel@conditional.x, cov_x = cov_x
    )
    ov_names <- tmp$ov_names
    ov <- tmp$ov
    ov2 <- tmp$ov2

    # which observed variables should NOT be standardized in this block?
    ov_nostd <- ov_nostd_block(g, ov_names)

    # observed product terms in component mode: scaled by SD(A)*SD(B)
    # (as computed by lav_standardize_lv), and they follow their
    # components rather than the ov_nostd set
    if (component_mode && length(ov_int_sd[[g]]) > 0L) {
      int_sd <- ov_int_sd[[g]]
      m_idx <- match(names(int_sd), ov_names)
      ov[m_idx] <- int_sd
      ov2[m_idx] <- int_sd^2
      ov_nostd <- setdiff(ov_nostd, names(int_sd))
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      !(partable$rhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$rhs[idx], ov_names)]

    # 1b. "=~" regular higher-order lv indicators -> handled in std.lv
    # 1c. "=~" indicators that are both in ov and lv -> not standardized here

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov_names & !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov_names & !(partable$rhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$rhs[idx], ov_names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] /
      ov2[match(partable$lhs[rv_idx], ov_names)])

    # covariances ov
    # three types:
    # - only lhs is OV (and fixed.x = FALSE)
    # - only rhs is OV (and fixed.x = FALSE)
    # - both lhs and rhs are OV (regular case)
    if (cov_std) {
      # est_c: est with the product-term variances centered (component
      # mode); otherwise identical to est
      if (!is.complex(est_c[rv_idx])) {
        rv <- sqrt(abs(est_c[rv_idx]))
      } else {
        rv <- sqrt(est_c[rv_idx])
      }
      rv_names <- partable$lhs[rv_idx]
    }

    # left
    idx_lhs <- which(partable$op == "~~" &
      !(partable$lhs %in% lv_names) & !(partable$lhs %in% ov_nostd) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_lhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_lhs] <-
          (out[idx_lhs] / ov[match(partable$lhs[idx_lhs], ov_names)])
      } else {
        out[idx_lhs] <-
          (out[idx_lhs] / rv[match(partable$lhs[idx_lhs], rv_names)])
      }
    }

    # right
    idx_rhs <- which(partable$op == "~~" &
      !(partable$rhs %in% lv_names) & !(partable$rhs %in% ov_nostd) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_rhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_rhs] <-
          (out[idx_rhs] / ov[match(partable$rhs[idx_rhs], ov_names)])
      } else {
        out[idx_rhs] <-
          (out[idx_rhs] / rv[match(partable$rhs[idx_rhs], rv_names)])
      }
    }

    # 3b. "~~" lv -> handled in std.lv

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # Delta-augmented ML correlation mode (FREE ~*~ scales): the mean
    # parameters already live in the STANDARDIZED metric (the implied
    # mean is Mu = Delta mu*), and VY is in the unit-diagonal P metric
    # -- so the generic division above (by VY = 1) is a no-op and the
    # "~1" rows are already correct: nothing extra to do here.

    # 4b. "~1" lv -> handled in std.lv

    # 4c. "|" thresholds
    idx <- which(partable$op == "|" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # 4d. "~*~" scales
    idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_nostd) &
      partable$block == g)
    out[idx] <- 1.0
  }

  # 5. defined/(in)equality constraints
  lav_standardize_constraints(out, partable, lavmodel)
}

# std.nox: standardize all observed variables except the exogenous 'x'.
# Thin wrapper around lav_standardize_all().
lav_standardize_all_nox <- function(lavobject = NULL,
                                    partable = NULL, est = NULL, est_std = NULL,
                                    glist = NULL, cov_std = TRUE, ov_var = NULL,
                                    lv_var = NULL,
                                    lavmodel = NULL, lavpartable = NULL,
                                    cov_x = NULL) {
  lav_standardize_all(
    lavobject = lavobject, partable = partable, est = est, est_std = est_std,
    glist = glist, cov_std = cov_std, ov_var = ov_var, lv_var = lv_var,
    lavmodel = lavmodel, lavpartable = lavpartable, cov_x = cov_x,
    nox = TRUE
  )
}

# ---------------------------------------------------------------------------
# unstandardize observed variables (used by e.g. simulateData)
# ---------------------------------------------------------------------------

lav_unstandardize_ov <- function(partable, ov_var = NULL, cov_std = TRUE) {
  # check if ustart is missing; if so, look for est
  if (is.null(partable$ustart)) {
    partable$ustart <- partable$est
  }

  # check if block is missing
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$ustart))
  }

  stopifnot(!any(is.na(partable$ustart)))
  est <- out <- partable$ustart

  # nblocks
  nblocks <- lav_pt_nblocks(partable)

  # if ov.var is NOT a list, make a list
  if (!is.list(ov_var)) {
    tmp <- ov_var
    ov_var <- vector("list", length = nblocks)
    ov_var[1:nblocks] <- list(tmp)
  }

  for (g in 1:nblocks) {
    ov_names <- lav_pt_vnames(partable, "ov", block = g) # not user
    lv_names <- lav_pt_vnames(partable, "lv", block = g)

    ov <- sqrt(ov_var[[g]])

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$rhs[idx], ov_names)]

    # 1b. "=~" regular higher-order lv indicators
    # 1c. "=~" indicators that are both in ov and lv

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$lhs[idx], ov_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$rhs[idx], ov_names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] * ov[match(partable$lhs[rv_idx], ov_names)]
      * ov[match(partable$rhs[rv_idx], ov_names)])

    # covariances
    idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx) > 0L) {
      if (cov_std == FALSE) {
        out[idx] <- (out[idx] * ov[match(partable$lhs[idx], ov_names)]
          * ov[match(partable$rhs[idx], ov_names)])
      } else {
        if (!is.complex(out[rv_idx])) {
          rv <- sqrt(abs(out[rv_idx]))
        } else {
          rv <- sqrt(out[rv_idx])
        }
        rv_names <- partable$lhs[rv_idx]
        out[idx] <- (out[idx] * rv[match(partable$lhs[idx], rv_names)]
          * rv[match(partable$rhs[idx], rv_names)])
      }
    }

    # 3b. "~~" lv

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$lhs[idx], ov_names)]

    # 4b. "~1" lv
  }

  # 5a ":="
  # 5b "=="
  # 5c. "<" or ">"

  out
}
