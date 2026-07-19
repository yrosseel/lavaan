# the 'multiple group' method (MGM) as described in Guttman, 1952
#
# Guttman, L. (1952). Multiple group methods for common-factor analysis,
# their basis, computation, and interpretation. Psychometrika, 17(2) 209--222
#
# YR 02 Feb 2023: - first version in lavaan, using quadprog (no std.lv yet)
# YR 19 Jul 2026: - rewrite:
#                   * the constrained LS problem is separable under simple
#                     structure, so the quadprog dependency is replaced by an
#                     exact closed-form solve (the normalized factor structure
#                     coefficients, averaged within equality classes)
#                   * multiple groups, including cross-group equality
#                     constraints among the loadings (measurement invariance)
#                   * optional second stage (shared with the IV/JS
#                     estimators) for all variances/covariances, enabling
#                     correlated residuals
#                   * mean structure (via the joint WLS mean solve)

# per-block MGM ingredients: given the sample covariance matrix 's', the
# 0/1 pattern matrix 'm_b' (nvar x nfac) and the residual variances 'theta',
# compute the reduced covariance matrix (SminTheta), the normalized
# sum-score correlation matrix (phi) and the normalized factor structure
# matrix (ys_cor, nfac x nvar)
lav_cfa_mgm_block <- function(s, m_b, theta,
                              theta_bounds = FALSE,
                              force_pd = FALSE,
                              nobs = 20L,
                              force_pd_shrink = NULL) {
  # force_pd_shrink: externally supplied shrinkage factor for the pd
  # correction (used by the pooled multigroup solve, where a JOINT factor
  # -- computed from the smallest root and the smallest group size over
  # all groups -- keeps the tied loadings comparable across groups);
  # NA means "no shrinkage needed"
  nvar <- ncol(s)

  # do we first 'clip' the theta values so they are within standard bounds?
  # (Question: do we need the 0.01 and 0.99 multipliers?)
  diag_s <- diag(s)
  if (theta_bounds) {
    # lower bound
    lower_bound <- diag_s * 0 # * 0.01
    too_small_idx <- which(theta < lower_bound)
    if (length(too_small_idx) > 0L) {
      theta[too_small_idx] <- lower_bound[too_small_idx]
    }

    # upper bound
    upper_bound <- diag_s * 1 # * 0.99
    too_large_idx <- which(theta > upper_bound)
    if (length(too_large_idx) > 0L) {
      theta[too_large_idx] <- upper_bound[too_large_idx]
    }
  }

  # compute SminTheta: S where we replace diagonal with 'communalities'
  diag_theta <- diag(theta, nvar)
  s_min_theta <- s - diag_theta
  if (force_pd && !is.null(force_pd_shrink)) {
    if (!is.na(force_pd_shrink)) {
      s_min_theta <- s - force_pd_shrink * diag_theta
    } # else: no shrinkage needed
  } else if (force_pd) {
    lambda <- try(lav_mat_sym_diff_smallest_root(s, diag_theta),
      silent = TRUE
    )
    if (inherits(lambda, "try-error")) {
      lav_msg_warn(gettext("failed to compute lambda"))
      s_min_theta <- s - diag_theta # and hope for the best
    } else {
      cutoff <- 1 + 1 / (nobs - 1)
      if (lambda < cutoff) {
        lambda_star <- lambda - 1 / (nobs - 1)
        s_min_theta <- s - lambda_star * diag_theta
      } else {
        s_min_theta <- s - diag_theta
      }
    }
  } else {
    # at least we force the diagonal elements of SminTheta to be nonnegative
    lower_bound <- diag_s * 0.001
    too_small_idx <- which(diag(s_min_theta) < lower_bound)
    if (length(too_small_idx) > 0L) {
      diag(s_min_theta)[too_small_idx] <- lower_bound[too_small_idx]
    }
  }

  # compute covariances among 1) (corrected) variables, and
  #                           2) (corrected) sum-scores
  ys_cov <- s_min_theta %*% m_b

  # compute covariance matrix of corrected sum-scores
  # SS.COV <- t(B) %*% SminTheta %*% B
  ss_cov <- crossprod(m_b, ys_cov)

  # scaling factors
  # D.inv.sqrt <- diag(1/sqrt(diag(SS.COV)))
  d_inv_sqrt <- 1 / sqrt(diag(ss_cov))

  # factor correlation matrix
  # PHI <- D.inv.sqrt %*% SS.COV %*% D.inv.sqrt
  phi <- t(ss_cov * d_inv_sqrt) * d_inv_sqrt

  # factor *structure* matrix
  # (covariances corrected Y & corrected normalized sum-scores)
  # YS.COR <- YS.COV %*% D.inv.sqrt
  ys_cor <- t(ys_cov) * d_inv_sqrt # transposed!

  list(
    theta = theta, s_min_theta = s_min_theta,
    phi = phi, ys_cor = ys_cor, d_inv_sqrt = d_inv_sqrt
  )
}

# 'purge' the residual covariances from the sample covariance matrix:
# for each residual-covariance pair (i, j), the common-factor part
# lambda_i lambda_j psi is identified by the tetrad ratio
#
#    s_ia * s_jb / s_ab
#
# with (a, b) clean anchor cells: for a within-factor pair both anchors
# come from the same factor (minus i, j); for a between-factor pair the
# anchor a is taken from j's factor and b from i's factor (so that the
# ratio identifies lambda_i lambda_j psi_fg). The purged matrix S0 equals
# S with each contaminated cell replaced by (the average of) these ratios
# -- i.e., S minus the residual covariances -- so that the row sums, the
# sum-score covariances and the Spearman tetrads downstream are all clean.
# Fixed (nonzero) residual covariances are subtracted directly.
#
# fac_of: for each observed variable, the (single) factor it loads on
#         (NA if none); pairs: data.frame(i, j, value) with value = NA
#         for a free residual covariance; contam: symmetric logical
#         matrix flagging all contaminated cells
lav_cfa_mgm_purge <- function(s, fac_of, pairs, contam) {
  s0 <- s
  failed <- integer(0L) # rows of 'pairs' without any clean anchor pair
  for (r in seq_len(nrow(pairs))) {
    i <- pairs$i[r]
    j <- pairs$j[r]
    val <- pairs$value[r]
    if (!is.na(val)) {
      s0[i, j] <- s0[j, i] <- s[i, j] - val
      next
    }
    fi <- fac_of[i]
    fj <- fac_of[j]
    if (is.na(fi) || is.na(fj)) {
      # e.g., a cross-loading indicator: no single-factor tetrad ratio
      failed <- c(failed, r)
      next
    }
    # admissible ordered anchor pairs (a, b): writing the tetrad ratio in
    # terms of the model, s_ia s_jb / s_ab = lambda_i lambda_j
    # psi[fi, g_a] psi[fj, g_b] / psi[g_a, g_b], which reduces to
    # lambda_i lambda_j psi[fi, fj] iff
    #  - within-factor pair (fi == fj): at least one anchor loads on fi
    #  - between-factor pair: both anchors load on fi or fj, except the
    #    orientation (g_a, g_b) = (fi, fj)
    cand <- setdiff(seq_along(fac_of), c(i, j))
    vals <- numeric(0L)
    for (a in cand) {
      ga <- fac_of[a]
      if (is.na(ga)) {
        next
      }
      for (b in cand) {
        if (a == b) {
          next
        }
        gb <- fac_of[b]
        if (is.na(gb)) {
          next
        }
        if (fi == fj) {
          if (ga != fi && gb != fi) {
            next
          }
        } else {
          if (!(ga %in% c(fi, fj)) || !(gb %in% c(fi, fj)) ||
              (ga == fi && gb == fj)) {
            next
          }
        }
        if (contam[i, a] || contam[j, b] || contam[a, b]) {
          next
        }
        vals <- c(vals, s[i, a] * s[j, b] / s[a, b])
      }
    }
    if (length(vals) == 0L) {
      failed <- c(failed, r)
      next
    }
    s0[i, j] <- s0[j, i] <- mean(vals)
  }
  attr(s0, "failed") <- failed
  s0
}

# the ML/Bartlett mapping function for PSI: given LAMBDA (marker metric),
# the residual variances and the reduced covariance matrix
lav_cfa_mgm_psi_mapping <- function(lambda, theta, s_min_theta) {
  nvar <- nrow(lambda)
  ti <- 1 / theta
  zero_theta_idx <- which(abs(theta) < 0.01) # be conservative
  if (length(zero_theta_idx) > 0L) {
    ti[zero_theta_idx] <- 1
  }

  # ML mapping function
  m <- solve(t(lambda) %*% diag(ti, nvar) %*% lambda) %*%
    t(lambda) %*% diag(ti, nvar)
  m %*% s_min_theta %*% t(m)
}

# stage 1 for a block WITH cross-loadings.
#
# The sum-score composites lose none of the loading information: at the
# population the unconstrained MGM solution equals Lambda* = Lambda G for
# some nonsingular k x k mixing G (G = Psi A' (A Psi A')^-1 with
# A = X'Lambda; G is diagonal only under simple structure). The mixing is
# recovered from the model's own zero pattern:
#
# (1) rotation-to-pattern: column h_f of H is the least-squares null
#     vector of the CLEAN zero rows of Lambda* (indicators without
#     cross-loadings that have a specified zero on factor f), scaled so
#     that the marker row gives 1; Lambda = Lambda* H. A wrong
#     communality for a cross-loading indicator perturbs only its OWN
#     row of Lambda* (the perturbation of SminTheta is diagonal), so the
#     clean rows -- and hence H -- do not depend on it.
# (2) PSI from the Bartlett/ML mapping over the clean rows only.
# (3) each cross-loading row is re-estimated against sum scores that
#     EXCLUDE the cross-loading indicators:
#        cov(y_i, t_masked) = lambda_i' Psi (X_masked' Lambda)'
#     which contains no theta_i (y_i is not in its own composites) and
#     is linear in lambda_i: a small exact/LS solve per row.
# (4) theta_i of a cross-loading indicator from the model identity
#     s_ii - lambda_i' Psi lambda_i.
#
# All steps are closed-form; with exact input moments the true values
# are recovered exactly.
lav_cfa_mgm_cross_block <- function(blk, m_b, marker_idx, s_use, theta,
                                    cross_idx,
                                    failed_vars = integer(0L)) {
  nvar <- nrow(m_b)
  nfac <- ncol(m_b)
  clean <- setdiff(seq_len(nvar), cross_idx)
  # rows touching an unpurged residual covariance are not usable as
  # zero references (their Lambda* rows are contaminated)
  zclean <- setdiff(clean, failed_vars)

  # (1) unconstrained solution (any column scaling; absorbed by H)
  lstar <- t(solve(blk$phi, blk$ys_cor))

  hmat <- matrix(0, nfac, nfac)
  for (f in seq_len(nfac)) {
    zrows <- intersect(which(m_b[, f] == 0L), zclean)
    if (length(zrows) < nfac - 1L) {
      lav_msg_stop(gettext(
        "unable to identify the factors from the zero pattern: the
         rotation step needs, for every factor, at least (nfac - 1)
         indicators without cross-loadings that do not load on it"))
    }
    zf <- lstar[zrows, , drop = FALSE]
    mrow <- lstar[marker_idx[f], ]
    kkt <- rbind(cbind(crossprod(zf), mrow), c(mrow, 0))
    sol <- try(solve(kkt, c(numeric(nfac), 1)), silent = TRUE)
    if (inherits(sol, "try-error")) {
      lav_msg_stop(gettext(
        "the rotation-to-pattern step failed (singular system); the
         zero pattern may not identify the factors"))
    }
    hmat[, f] <- sol[seq_len(nfac)]
  }
  lambda <- (lstar %*% hmat) * m_b # zero the off-pattern LS noise

  # (2) PSI from the clean part
  psi <- lav_cfa_mgm_psi_mapping(
    lambda = lambda[clean, , drop = FALSE], theta = theta[clean],
    s_min_theta = blk$s_min_theta[clean, clean, drop = FALSE]
  )

  # (3) cross-loading rows from masked composites
  x_masked <- m_b
  x_masked[cross_idx, ] <- 0
  if (length(failed_vars) > 0L) {
    x_masked[failed_vars, ] <- 0 # unpurged cells stay out of the sums
  }
  bmat <- psi %*% t(crossprod(x_masked, lambda))
  for (i in cross_idx) {
    fidx <- which(m_b[i, ] == 1L)
    ci <- drop(s_use[i, , drop = FALSE] %*% x_masked)
    bi <- bmat[fidx, , drop = FALSE]
    lam_i <- try(solve(tcrossprod(bi), bi %*% ci), silent = TRUE)
    if (inherits(lam_i, "try-error")) {
      lav_msg_stop(gettext(
        "unable to estimate the loadings of a cross-loading indicator
         (singular masked composite system)"))
    }
    lambda[i, ] <- 0
    lambda[i, fidx] <- lam_i
  }

  # (4) residual variances of the cross-loading indicators
  for (i in cross_idx) {
    theta[i] <- s_use[i, i] -
      drop(lambda[i, , drop = FALSE] %*% psi %*% lambda[i, ])
  }

  list(lambda = lambda, psi = psi, theta = theta)
}

lav_cfa_guttman1952 <- function(s,
                                marker_idx = NULL,
                                lambda_nonzero_idx = NULL,
                                theta = NULL, # vector!
                                theta_bounds = FALSE,
                                force_pd = FALSE,
                                zero_after_efa = FALSE,
                                quadprog = FALSE,
                                psi_mapping = FALSE,
                                nobs = 20L) { # for cutoff
  # dimensions
  nvar <- ncol(s)
  nfac <- length(marker_idx)
  stopifnot(length(theta) == nvar)

  # overview of lambda structure
  m_b <- matrix(0, nvar, nfac)
  lambda_marker_idx <- (seq_len(nfac) - 1L) * nvar + marker_idx
  m_b[lambda_marker_idx] <- 1L
  m_b[lambda_nonzero_idx] <- 1L

  # this method does not support crossloadings!
  if (any(rowSums(m_b) > 1)) {
    lav_msg_stop(gettext("the guttman1952 procedure does not support ",
                       "crossloadings; consider fabin or bentler1982 instead"))
  }

  # if we wish to keep SminTheta PD, we must keep theta within bounds
  if (force_pd) {
    theta_bounds <- TRUE
  }
  if (psi_mapping) {
    theta_bounds <- TRUE
    force_pd <- TRUE
  }

  blk <- lav_cfa_mgm_block(
    s = s, m_b = m_b, theta = theta,
    theta_bounds = theta_bounds, force_pd = force_pd, nobs = nobs
  )
  theta <- blk$theta
  s_min_theta <- blk$s_min_theta
  phi <- blk$phi
  ys_cor <- blk$ys_cor

  if (zero_after_efa) {
    # we initially assume a saturated LAMBDA (like EFA)
    # then, we just fix the zero-elements to zero

    mm_lambda <- t(solve(phi, ys_cor)) # = unconstrained EFA version
    # force zeroes
    mm_lambda <- mm_lambda * m_b
  } else if (quadprog) {
    # constrained version: only the pattern (nonzero) cells are free, the
    # others are fixed to zero. Under simple structure this quadratic
    # program is separable per variable, and (because the normalized
    # sum-score correlation matrix has a unit diagonal) the solution for
    # each pattern cell is simply its normalized structure coefficient --
    # no quadratic-program solver is needed.
    # (historically this branch was solved with quadprog::solve.QP; the
    # closed form gives the identical solution)
    phi <- cov2cor(lav_mat_sym_force_pd(phi, tol = 1e-04)) # (kept: psi
    #                                     below is computed from this phi)
    mm_lambda <- t(ys_cor) * m_b
  } else {
    # default, if no (in)equality constraints
    mm_lambda <- t(solve(phi, ys_cor)) # = unconstrained EFA version
    #YS.COR0 <- YS.COR
    #YS.COR0[t(B) != 1] <- 0
    #LAMBDA <- t(YS.COR0)
  }

  # rescale LAMBDA, so that 'marker' indicator == 1
  marker_lambda <- mm_lambda[lambda_marker_idx]
  lambda_1 <- t(t(mm_lambda) * (1 / marker_lambda))

  # rescale PHI, covariance metric
  psi <- t(phi * marker_lambda) * marker_lambda

  # redo psi using ML mapping function?
  if (psi_mapping) {
    psi <- lav_cfa_mgm_psi_mapping(
      lambda = lambda_1, theta = theta,
      s_min_theta = s_min_theta
    )
  }

  list(lambda = lambda_1, theta = theta, psi = psi)
}

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_cfa_guttman1952_internal <- function(lavobject = NULL, # convenience
                                         # internal slot
                                         lavmodel = NULL,
                                         lavsamplestats = NULL,
                                         lavpartable = NULL,
                                         lavdata = NULL,
                                         lavh1 = NULL,
                                         lavoptions = NULL,
                                         theta_bounds = TRUE,
                                         force_pd = TRUE,
                                         zero_after_efa = FALSE,
                                         quadprog = FALSE,
                                         psi_mapping = TRUE,
                                         mgm_varcov = NULL) {
  lavpta <- NULL
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))

    # extract slots
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavpartable <- lav_pt_set_cache(lavobject@ParTable, lavobject@pta)
    lavpta <- lavobject@pta
    lavdata <- lavobject@Data
    lavh1 <- lavobject@h1
    lavoptions <- lavobject@Options
  }
  if (is.null(lavpta)) {
    lavpta <- lav_pt_attributes(lavpartable)
    lavpartable <- lav_pt_set_cache(lavpartable, lavpta)
  }

  # options; NOTE: read with [[ ]] -- $ would partially match
  ea <- lavoptions$estimator.args
  if (missing(zero_after_efa) && !is.null(ea[["zero.after.efa"]])) {
    zero_after_efa <- ea[["zero.after.efa"]]
  }
  if (missing(psi_mapping) && !is.null(ea[["psi.mapping"]])) {
    psi_mapping <- ea[["psi.mapping"]]
  }
  if (missing(quadprog) && !is.null(ea[["quadprog"]])) {
    quadprog <- ea[["quadprog"]]
  }
  if (is.null(mgm_varcov)) {
    mgm_varcov <- ea[["mgm.varcov"]]
    if (is.null(mgm_varcov)) {
      mgm_varcov <- ea[["mgm_varcov"]]
    }
    if (is.null(mgm_varcov)) {
      mgm_varcov <- "default"
    }
  }
  mgm_varcov <- toupper(mgm_varcov)

  # no structural part!
  if (any(lavpartable$op == "~")) {
    lav_msg_stop(gettext(
      "the MGM (guttman1952) estimator is only available for CFA models"))
  }
  # no BETA matrix! (i.e., no higher-order factors)
  if (!is.null(lavmodel@GLIST$beta)) {
    lav_msg_stop(gettext(
      "the MGM (guttman1952) estimator is not available for models that
       require a BETA matrix"))
  }
  # no std.lv = TRUE for now
  if (lavoptions$std.lv) {
    lav_msg_stop(gettext(
      "the MGM (guttman1952) estimator is not available if std.lv = TRUE"))
  }
  # only simple (a == b) equality constraints are supported (those are
  # absorbed by ceq.simple = TRUE; anything left over is more general)
  if (length(lavmodel@ceq.linear.idx) > 0L ||
      length(lavmodel@ceq.nonlinear.idx) > 0L) {
    lav_msg_stop(gettext(
      "the MGM (guttman1952) estimator only supports simple (a == b)
       equality constraints"))
  }

  nblocks <- lav_pt_nblocks(lavpartable)
  pt_block <- if (!is.null(lavpartable$block)) {
    lavpartable$block
  } else {
    rep(1L, length(lavpartable$op))
  }

  # ------------------------------------------------------------------
  # collect the per-block ingredients
  # ------------------------------------------------------------------
  block_info <- vector("list", nblocks)
  cells <- NULL # all pattern (nonzero) lambda cells, over all blocks
  res_cov_flag <- FALSE
  for (b in seq_len(nblocks)) {
    ov_names <- lavpta$vnames$ov[[b]]
    lv_names <- lavpta$vnames$lv.regular[[b]]
    nvar <- length(ov_names)
    nfac <- length(lv_names)
    marker_idx <- lavpta$vidx$lv.marker[[b]]
    sample_cov <- lavsamplestats@cov[[b]]
    nobs_b <- lavsamplestats@nobs[[b]]

    # free lambda cells (from the parameter table; with ceq.simple = TRUE
    # equality-constrained cells share the same 'free' index)
    lam_pt_idx <- which(lavpartable$op == "=~" & pt_block == b &
      lavpartable$free > 0L)
    lam_ov <- match(lavpartable$rhs[lam_pt_idx], ov_names)
    lam_fac <- match(lavpartable$lhs[lam_pt_idx], lv_names)
    if (anyNA(lam_ov) || anyNA(lam_fac)) {
      lav_msg_stop(gettext(
        "the MGM (guttman1952) estimator requires all loadings to connect
         an observed indicator to a (first-order) latent variable"))
    }

    # pattern matrix (markers + free cells)
    m_b <- matrix(0, nvar, nfac)
    lambda_marker_idx <- (seq_len(nfac) - 1L) * nvar + marker_idx
    m_b[lambda_marker_idx] <- 1L
    m_b[(lam_fac - 1L) * nvar + lam_ov] <- 1L

    # cross-loadings are handled by the rotation-to-pattern step below,
    # but the scaling (marker) indicators must be free of them (their row
    # anchors the metric of exactly one factor)
    cross_idx <- which(rowSums(m_b) > 1L)
    if (any(marker_idx %in% cross_idx)) {
      lav_msg_stop(gettext(
        "the MGM (guttman1952) estimator requires the scaling (marker)
         indicators to be free of cross-loadings; choose a different
         scaling indicator"))
    }

    # correlated residuals (free, or fixed to a nonzero value)? purge them
    # from the sample covariance matrix, so that the row sums, the
    # sum-score covariances and the Spearman tetrads downstream are clean
    contam <- lav_sem_js_contam(
      lavpartable = lavpartable, pt_block = pt_block, b = b,
      ov_names = ov_names, lavpta = lavpta
    )
    s_use <- sample_cov
    if (any(contam)) {
      res_cov_flag <- TRUE
      fac_of <- rep(NA_integer_, nvar)
      for (f in seq_len(nfac)) {
        fac_of[which(m_b[, f] == 1L)] <- f
      }
      if (length(cross_idx) > 0L) {
        # cross-loading indicators cannot serve as purge anchors, and a
        # pair involving one cannot be purged by the single-factor ratio
        fac_of[cross_idx] <- NA_integer_
      }
      cov_idx <- which(lavpartable$op == "~~" & pt_block == b &
        lavpartable$lhs != lavpartable$rhs &
        (lavpartable$free > 0L |
          (!is.na(lavpartable$ustart) & lavpartable$ustart != 0)) &
        lavpartable$lhs %in% ov_names & lavpartable$rhs %in% ov_names)
      pairs <- data.frame(
        i = match(lavpartable$lhs[cov_idx], ov_names),
        j = match(lavpartable$rhs[cov_idx], ov_names),
        value = ifelse(lavpartable$free[cov_idx] > 0L,
          NA_real_, lavpartable$ustart[cov_idx])
      )
      s_use <- lav_cfa_mgm_purge(
        s = s_use, fac_of = fac_of, pairs = pairs, contam = contam
      )
    }
    failed_idx <- attr(s_use, "failed") # pairs without clean anchors
    attr(s_use, "failed") <- NULL
    failed_vars <- integer(0L)
    contam_left <- matrix(FALSE, nvar, nvar)
    if (length(failed_idx) > 0L) {
      lav_msg_warn(gettext(
        "unable to purge some residual covariances (no clean anchor
         cells); the corresponding loadings may be distorted"))
      failed_vars <- unique(c(pairs$i[failed_idx], pairs$j[failed_idx]))
      contam_left[cbind(pairs$i[failed_idx], pairs$j[failed_idx])] <- TRUE
      contam_left[cbind(pairs$j[failed_idx], pairs$i[failed_idx])] <- TRUE
    }

    # initial estimate for (the diagonal elements of) THETA: Spearman per
    # factor on the purged matrix (all tetrads are clean after purging);
    # factors with fewer than 3 indicators -- or factors touching an
    # unpurged residual covariance or a cross-loading -- use
    # masked/external tetrads. A single cross-loading ANCHOR cancels out
    # of the tetrad ratio (like the purge anchors), but a pair of
    # cross-loading anchors does not: mask those pairs. The (biased)
    # tetrad value for a cross-loading indicator itself is only a
    # placeholder: a wrong communality there perturbs nothing but its own
    # row of the unconstrained solution, which the rotation step ignores,
    # and its final theta is recomputed from the model identity.
    both_cross <- matrix(FALSE, nvar, nvar)
    if (length(cross_idx) > 1L) {
      both_cross[cross_idx, cross_idx] <- TRUE
    }
    contam_theta <- contam_left | both_cross
    theta <- numeric(nvar)
    for (f in seq_len(nfac)) {
      ov_idx <- which(m_b[, f] == 1L)
      if (length(ov_idx) >= 3L && !any(ov_idx %in% failed_vars) &&
          !any(ov_idx %in% cross_idx)) {
        s_fac <- s_use[ov_idx, ov_idx, drop = FALSE]
        theta[ov_idx] <- lav_cfa_theta_spearman(s_fac, bounds = "wide")
      } else {
        theta[ov_idx] <- lav_sem_js_theta_spearman_masked(
          s_full = s_use, idx = ov_idx, contam = contam_theta,
          zpool = setdiff(seq_len(nvar), ov_idx), bounds = "wide"
        )
      }
    }
    if (anyNA(theta)) {
      lav_msg_stop(gettext(
        "unable to obtain an initial residual-variance estimate for some
         indicators (not enough clean tetrads); factors need at least two
         indicators, and a two-indicator factor must be correlated with
         the rest of the model"))
    }

    block_info[[b]] <- list(
      ov_names = ov_names, lv_names = lv_names, nvar = nvar, nfac = nfac,
      marker_idx = marker_idx, m_b = m_b, theta = theta,
      sample_cov = s_use, nobs = nobs_b,
      cross_idx = cross_idx, failed_vars = failed_vars,
      lambda_nonzero_idx = (lam_fac - 1L) * nvar + lam_ov
    )

    if (length(lam_pt_idx) > 0L) {
      cells <- rbind(cells, data.frame(
        block = b, ov = lam_ov, fac = lam_fac,
        facname = lavpartable$lhs[lam_pt_idx],
        free = lavpartable$free[lam_pt_idx]
      ))
    }
  }

  # equality constraints among the loadings? (shared 'free' indices)
  if (is.null(cells)) {
    cells <- data.frame(
      block = integer(0L), ov = integer(0L), fac = integer(0L),
      facname = character(0L), free = integer(0L)
    )
  }
  lam_free_ids <- unique(cells$free)
  tie_classes <- split(seq_len(nrow(cells)), cells$free)
  tie_classes <- tie_classes[vapply(tie_classes, length, 1L) > 1L]
  has_ties <- length(tie_classes) > 0L

  # a loading tied to a non-loading parameter cannot be honored here
  other_free <- lavpartable$free[lavpartable$op != "=~" &
    lavpartable$free > 0L]
  if (any(lam_free_ids %in% other_free)) {
    lav_msg_stop(gettext(
      "the MGM (guttman1952) estimator does not support equality
       constraints between a factor loading and another type of
       parameter"))
  }

  # tie classes must stay within a single factor (the constrained solve
  # operates in the normalized sum-score metric, which is factor-specific)
  if (has_ties) {
    for (class_idx in tie_classes) {
      if (length(unique(cells$facname[class_idx])) > 1L) {
        lav_msg_stop(gettext(
          "the MGM (guttman1952) estimator does not support equality
           constraints among loadings of different factors"))
      }
    }
  }

  # second stage for the variances/covariances?
  if (mgm_varcov == "DEFAULT") {
    mgm_varcov <- if (res_cov_flag) "RLS" else "NONE"
  } else if (mgm_varcov == "NONE" && res_cov_flag) {
    # residual covariances require the second stage
    mgm_varcov <- "RLS"
  }

  # ------------------------------------------------------------------
  # stage 1: the loadings (+ initial theta/psi)
  # ------------------------------------------------------------------
  gl_lambda_idx <- which(names(lavmodel@GLIST) == "lambda")
  gl_theta_idx <- which(names(lavmodel@GLIST) == "theta")
  gl_psi_idx <- which(names(lavmodel@GLIST) == "psi")

  any_cross <- any(vapply(block_info,
    function(bi) length(bi$cross_idx) > 0L, logical(1L)))

  # effective theta_bounds/force_pd (same implications as the worker)
  tb_eff <- theta_bounds
  fp_eff <- force_pd
  if (fp_eff) {
    tb_eff <- TRUE
  }
  if (psi_mapping) {
    tb_eff <- TRUE
    fp_eff <- TRUE
  }

  # JOINT pd-correction factor over all blocks (smallest root, smallest
  # group size): with cross-group ties, the same fraction of theta must
  # be subtracted in every group, so that the pooled (averaged)
  # estimates remain comparable across groups
  force_pd_shrink <- NULL
  if (fp_eff && nblocks > 1L && (has_ties || quadprog) && !any_cross) {
    roots <- numeric(nblocks)
    roots_ok <- TRUE
    nobs_min <- Inf
    for (b in seq_len(nblocks)) {
      bi <- block_info[[b]]
      th_b <- bi$theta
      if (tb_eff) {
        th_b <- pmin(pmax(th_b, 0), diag(bi$sample_cov))
      }
      r_b <- try(lav_mat_sym_diff_smallest_root(
        bi$sample_cov, diag(th_b, bi$nvar)), silent = TRUE)
      if (inherits(r_b, "try-error")) {
        roots_ok <- FALSE
        break
      }
      roots[b] <- r_b
      nobs_min <- min(nobs_min, bi$nobs)
    }
    if (roots_ok) {
      cutoff <- 1 + 1 / (nobs_min - 1)
      root_min <- min(roots)
      force_pd_shrink <- if (root_min < cutoff) {
        root_min - 1 / (nobs_min - 1)
      } else {
        NA_real_ # no shrinkage needed
      }
    }
  }

  if (any_cross) {
    # cross-loadings: rotation-to-pattern. The unconstrained MGM solution
    # recovers the loading space exactly (Lambda* = Lambda G for some
    # nonsingular k x k mixing G); the mixing is undone per factor from
    # the clean zero rows, and the cross-loading rows are re-estimated
    # against composites that exclude them (see lav_cfa_mgm_cross_block).
    # PSI comes from the Bartlett/ML mapping over the clean rows (the
    # sum-score correlation matrix is not usable here: with
    # cross-loadings it mixes the factors); the psi.mapping /
    # zero.after.efa / quadprog options are ignored for such blocks.
    #
    # NOTE: the pd correction is DISABLED here. Its root is driven by the
    # (deliberately rough) communality placeholders of the cross-loading
    # indicators, and shrinking ALL residual variances would perturb the
    # clean rows of the unconstrained solution -- the rotation step
    # requires that perturbations of SminTheta stay confined to the
    # cross-loading rows (and none of the steps below needs SminTheta to
    # be positive definite).
    for (b in seq_len(nblocks)) {
      bi <- block_info[[b]]
      blk <- lav_cfa_mgm_block(
        s = bi$sample_cov, m_b = bi$m_b, theta = bi$theta,
        theta_bounds = TRUE, force_pd = FALSE, nobs = bi$nobs
      )
      out <- lav_cfa_mgm_cross_block(
        blk = blk, m_b = bi$m_b, marker_idx = bi$marker_idx,
        s_use = bi$sample_cov, theta = blk$theta,
        cross_idx = bi$cross_idx, failed_vars = bi$failed_vars
      )
      block_info[[b]]$lambda <- out$lambda
      block_info[[b]]$psi <- out$psi
      block_info[[b]]$theta <- out$theta
    }

    # equality ties: average the tied cells (final, marker metric), then
    # refresh the cross rows' residual variances with the tied values
    if (has_ties) {
      cells$value <- numeric(nrow(cells))
      for (r in seq_len(nrow(cells))) {
        cells$value[r] <-
          block_info[[cells$block[r]]]$lambda[cells$ov[r], cells$fac[r]]
      }
      for (class_idx in tie_classes) {
        cells$value[class_idx] <- mean(cells$value[class_idx])
      }
      for (r in seq_len(nrow(cells))) {
        block_info[[cells$block[r]]]$lambda[cells$ov[r], cells$fac[r]] <-
          cells$value[r]
      }
      for (b in seq_len(nblocks)) {
        bi <- block_info[[b]]
        for (i in bi$cross_idx) {
          block_info[[b]]$theta[i] <- bi$sample_cov[i, i] -
            drop(bi$lambda[i, , drop = FALSE] %*% bi$psi %*%
                 bi$lambda[i, ])
        }
      }
    }

    for (b in seq_len(nblocks)) {
      bi <- block_info[[b]]
      lavmodel@GLIST[[gl_lambda_idx[b]]] <- bi$lambda
      lavmodel@GLIST[[gl_theta_idx[b]]] <- diag(bi$theta, bi$nvar)
      lavmodel@GLIST[[gl_psi_idx[b]]] <- bi$psi
    }
  } else if (!has_ties && !quadprog) {
    # no equality constraints: the classic per-block computation
    for (b in seq_len(nblocks)) {
      bi <- block_info[[b]]
      out <- lav_cfa_guttman1952(
        s = bi$sample_cov, marker_idx = bi$marker_idx,
        lambda_nonzero_idx = bi$lambda_nonzero_idx,
        theta = bi$theta,
        theta_bounds = theta_bounds,
        force_pd = force_pd,
        zero_after_efa = zero_after_efa,
        quadprog = FALSE,
        psi_mapping = psi_mapping,
        nobs = bi$nobs
      )
      lavmodel@GLIST[[gl_lambda_idx[b]]] <- out$lambda
      lavmodel@GLIST[[gl_theta_idx[b]]] <- diag(out$theta, bi$nvar)
      lavmodel@GLIST[[gl_psi_idx[b]]] <- out$psi
    }
  } else {
    # equality constraints among the loadings: the constrained LS problem
    # in the normalized sum-score metric is separable per cell, so tied
    # cells are estimated by AVERAGING their normalized structure
    # coefficients. The normalized marker loadings (the rescaling factors)
    # are pooled over all groups that are connected through at least one
    # tied loading of that factor, so that tied loadings remain tied after
    # the rescaling to the marker metric.
    blks <- vector("list", nblocks)
    lstar <- vector("list", nblocks) # normalized structure coefficients
    for (b in seq_len(nblocks)) {
      bi <- block_info[[b]]
      blks[[b]] <- lav_cfa_mgm_block(
        s = bi$sample_cov, m_b = bi$m_b, theta = bi$theta,
        theta_bounds = tb_eff, force_pd = fp_eff, nobs = bi$nobs,
        force_pd_shrink = force_pd_shrink
      )
      block_info[[b]]$theta <- blks[[b]]$theta # after clipping
      lstar[[b]] <- t(blks[[b]]$ys_cor) # nvar x nfac
    }

    # per factor name: partition the blocks into connected components
    # (blocks joined by at least one tied loading of that factor); the
    # normalized marker loading is pooled within each component
    fac_names <- unique(cells$facname)
    mhat <- vector("list", nblocks) # per block: named vector, per factor
    for (b in seq_len(nblocks)) {
      bi <- block_info[[b]]
      mh <- numeric(bi$nfac)
      for (f in seq_len(bi$nfac)) {
        mh[f] <- lstar[[b]][bi$marker_idx[f], f]
      }
      names(mh) <- bi$lv_names
      mhat[[b]] <- mh
    }
    for (fn in fac_names) {
      # components over blocks, connected by tie classes of this factor
      comps <- as.list(seq_len(nblocks))
      for (class_idx in tie_classes) {
        if (cells$facname[class_idx[1L]] != fn) {
          next
        }
        bset <- unique(cells$block[class_idx])
        # merge all components that intersect bset
        hit <- which(vapply(comps, function(cc) any(cc %in% bset), TRUE))
        if (length(hit) > 1L) {
          comps[[hit[1L]]] <- sort(unique(unlist(comps[hit])))
          comps[hit[-1L]] <- NULL
        }
      }
      for (cc in comps) {
        if (length(cc) < 2L) {
          next
        }
        # pool the normalized marker loadings over the component
        vals <- numeric(0L)
        for (b in cc) {
          if (fn %in% names(mhat[[b]])) {
            vals <- c(vals, mhat[[b]][[fn]])
          }
        }
        if (length(vals) > 1L) {
          for (b in cc) {
            if (fn %in% names(mhat[[b]])) {
              mhat[[b]][[fn]] <- mean(vals)
            }
          }
        }
      }
    }

    # pooled values for the tied cells (in the normalized metric)
    cells$value <- numeric(nrow(cells))
    for (r in seq_len(nrow(cells))) {
      cells$value[r] <- lstar[[cells$block[r]]][cells$ov[r], cells$fac[r]]
    }
    for (class_idx in tie_classes) {
      cells$value[class_idx] <- mean(cells$value[class_idx])
    }

    # assemble the per-block matrices
    for (b in seq_len(nblocks)) {
      bi <- block_info[[b]]
      lambda_b <- matrix(0, bi$nvar, bi$nfac)
      # markers: exactly 1 (pooled marker / own pooled marker)
      lambda_b[(seq_len(bi$nfac) - 1L) * bi$nvar + bi$marker_idx] <- 1
      rows_b <- which(cells$block == b)
      for (r in rows_b) {
        lambda_b[cells$ov[r], cells$fac[r]] <-
          cells$value[r] / mhat[[b]][[cells$facname[r]]]
      }

      # psi: from the normalized sum-score correlations, rescaled by the
      # (pooled) normalized marker loadings; or via the ML mapping
      if (psi_mapping) {
        psi_b <- lav_cfa_mgm_psi_mapping(
          lambda = lambda_b, theta = block_info[[b]]$theta,
          s_min_theta = blks[[b]]$s_min_theta
        )
      } else {
        phi_b <- cov2cor(lav_mat_sym_force_pd(blks[[b]]$phi, tol = 1e-04))
        mh <- as.numeric(mhat[[b]][bi$lv_names])
        psi_b <- t(phi_b * mh) * mh
      }

      lavmodel@GLIST[[gl_lambda_idx[b]]] <- lambda_b
      lavmodel@GLIST[[gl_theta_idx[b]]] <- diag(block_info[[b]]$theta,
                                                bi$nvar)
      lavmodel@GLIST[[gl_psi_idx[b]]] <- psi_b
    }
  }

  # extract free parameters only
  x <- lav_model_get_parameters(lavmodel)

  # ------------------------------------------------------------------
  # stage 2 (optional): re-estimate ALL variances/covariances (including
  # residual covariances and the factor (co)variances) holding the
  # loadings fixed -- the same second stage as the IV/JS estimators
  # ------------------------------------------------------------------
  implied <- NULL
  if (!is.null(lavh1$implied$cov)) {
    implied <- lavh1$implied
  } else {
    implied <- list(cov = lavsamplestats@cov, mean = lavsamplestats@mean)
  }
  if (mgm_varcov != "NONE") {
    free_undirected_idx <- unique(lavpartable$free[lavpartable$op == "~~" &
      lavpartable$free > 0L])
    if (length(free_undirected_idx) > 0L) {
      lavmodel_tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
      delta_list <- lav_sem_miiv_delta(lavmodel_tmp)
      for (b in seq_len(nblocks)) {
        fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
          lavpartable$free > 0L & pt_block == b])
        fu_b <- fu_b[fu_b %in% free_undirected_idx]
        if (length(fu_b) == 0L) {
          next
        }
        tb <- lav_sem_miiv_varcov_block(
          b = b, fu = fu_b, delta_b = delta_list[[b]],
          implied = implied, lavmodel = lavmodel, x = x,
          iv_varcov_method = mgm_varcov, return_h = FALSE
        )
        x[fu_b] <- as.numeric(tb)
      }
    }
  }

  # mean structure: re-estimate all free mean parameters (observed
  # intercepts and latent means) jointly by GLS, as the IV/JS estimators do
  if (lavmodel@meanstructure) {
    free_mean_idx <- unique(lavpartable$free[lavpartable$op == "~1" &
      lavpartable$free > 0L])
    if (length(free_mean_idx) > 0L) {
      x[free_mean_idx] <- lav_sem_miiv_mean_wls(
        lavmodel = lavmodel, lavsamplestats = lavsamplestats, x = x,
        free_mean_idx = free_mean_idx,
        sample_mean = implied$mean
      )
    }
  }

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

  x
}
