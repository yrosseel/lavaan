# the 'multiple group' method as described in Guttman, 1952
#
# Guttman, L. (1952). Multiple group methods for common-factor analysis,
# their basis, computation, and interpretation. Psychometrika, 17(2) 209--222
#
# YR 02 Feb 2023: - first version in lavaan, using quadprog (no std.lv yet)

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
  if (force_pd) {
    lambda <- try(lav_matrix_symmetric_diff_smallest_root(s, diag_theta),
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

  if (zero_after_efa) {
    # we initially assume a saturated LAMBDA (like EFA)
    # then, we just fix the zero-elements to zero

    mm_lambda <- t(solve(phi, ys_cor)) # = unconstrained EFA version
    # force zeroes
    mm_lambda <- mm_lambda * m_b
  } else if (quadprog) {
    # constained version using quadprog
    # only useful if (in)equality constraints are needed (TODo)

    # PHI MUST be positive-definite
    phi <- cov2cor(lav_matrix_symmetric_force_pd(phi,
      tol = 1e-04
    )) # option?
    dmat <- lav_matrix_bdiag(rep(list(phi), nvar))
    dvec <- as.vector(ys_cor)
    eq_idx <- which(t(m_b) != 1) # these must be zero (row-wise!)
    rmat <- diag(nrow(dmat))[eq_idx, , drop = FALSE]
    bvec <- rep(0, length(eq_idx)) # optional, 0=default
    out <- try(quadprog::solve.QP(
      Dmat = dmat, dvec = dvec, Amat = t(rmat),
      meq = length(eq_idx), bvec = bvec
    ), silent = TRUE)
    if (inherits(out, "try-error")) {
      lav_msg_warn(gettext("solve.QP failed to find a solution"))
      lambda_1 <- m_b
      lambda_1[lambda_nonzero_idx] <- as.numeric(NA)
      theta_1 <- diag(rep(as.numeric(NA), nvar), nvar)
      psi <- matrix(as.numeric(NA), nfac, nfac)
      return(list(lambda = lambda_1, theta = theta_1, psi = psi))
    } else {
      mm_lambda <- matrix(out$solution, nrow = nvar, ncol = nfac, byrow = TRUE)
      # zap almost zero elements
      mm_lambda[abs(mm_lambda) < sqrt(.Machine$double.eps)] <- 0
    }
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
    ti <- 1 / theta
    zero_theta_idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero_theta_idx) > 0L) {
      ti[zero_theta_idx] <- 1
    }

    # ML mapping function
    m <- solve(t(lambda_1) %*% diag(ti, nvar) %*% lambda_1) %*%
         t(lambda_1) %*% diag(ti, nvar)
    psi <- m %*% s_min_theta %*% t(m)
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
                                         lavoptions = NULL,
                                         theta_bounds = TRUE,
                                         force_pd = TRUE,
                                         zero_after_efa = FALSE,
                                         quadprog = FALSE,
                                         psi_mapping = TRUE) {
  lavpta <- NULL
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))

    # extract slots
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavpartable <- lav_partable_set_cache(lavobject@ParTable, lavobject@pta)
    lavpta <- lavobject@pta
    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
  }
  if (is.null(lavpta)) {
    lavpta <- lav_partable_attributes(lavpartable)
    lavpartable <- lav_partable_set_cache(lavpartable, lavpta)
  }

  if (missing(zero_after_efa) &&
    !is.null(lavoptions$estimator.args$zero.after.efa)) {
    zero_after_efa <- lavoptions$estimator.args$zero.after.efa
  }

  if (missing(psi_mapping) &&
    !is.null(lavoptions$estimator.args$psi.mapping)) {
    psi_mapping <- lavoptions$estimator.args$psi.mapping
  }

  if (missing(quadprog) &&
    !is.null(lavoptions$estimator.args$quadprog)) {
    quadprog <- lavoptions$estimator.args$quadprog
  }

  # no structural part!
  if (any(lavpartable$op == "~")) {
    lav_msg_stop(gettext("GUTTMAN1952 estimator only available for CFA models"))
  }
  # no BETA matrix! (i.e., no higher-order factors)
  if (!is.null(lavmodel@GLIST$beta)) {
    lav_msg_stop(gettext(
   "GUTTMAN1952 estimator not available for models that require a BETA matrix"))
  }
  # no std.lv = TRUE for now
  if (lavoptions$std.lv) {
    lav_msg_stop(gettext(
      "GUTTMAN1952 estimator not available if std.lv = TRUE"))
  }

  nblocks <- lav_partable_nblocks(lavpartable)
  stopifnot(nblocks == 1L) # for now
  b <- 1L
  sample_cov <- lavsamplestats@cov[[b]]
  nvar <- nrow(sample_cov)
  lv_names <- lavpta$vnames$lv.regular[[b]]
  nfac <- length(lv_names)
  marker_idx <- lavpta$vidx$lv.marker[[b]]
  lambda_idx <- which(names(lavmodel@GLIST) == "lambda")
  lambda_nonzero_idx <- lavmodel@m.free.idx[[lambda_idx]]
  # only diagonal THETA for now...
  # because if we have correlated residuals, we should remove the
  # corresponding variables as instruments before we estimate lambda...
  # (see MIIV)
  theta_idx <- which(names(lavmodel@GLIST) == "theta") # usually '2'
  m_theta <- lavmodel@m.free.idx[[theta_idx]]
  nondiag_idx <- m_theta[!m_theta %in% lav_matrix_diag_idx(nvar)]
  if (length(nondiag_idx) > 0L) {
    lav_msg_warn(gettext(
      "this implementation of FABIN does not handle correlated residuals yet!"))
  }

  # 1. obtain estimate for (diagonal elements of) THETA
  #    for now we use Spearman per factor
  m_b <- matrix(0, nvar, nfac)
  lambda_marker_idx <- (seq_len(nfac) - 1L) * nvar + marker_idx
  m_b[lambda_marker_idx] <- 1L
  m_b[lambda_nonzero_idx] <- 1L
  theta <- numeric(nvar)
  for (f in seq_len(nfac)) {
    ov_idx <- which(m_b[, f] == 1L)
    s_fac <- sample_cov[ov_idx, ov_idx, drop = FALSE]
    theta[ov_idx] <- lav_cfa_theta_spearman(s_fac, bounds = "wide")
  }

  # 2. run Guttman1952 'Multiple Groups' algorithm
  out <- lav_cfa_guttman1952(
    s = sample_cov, marker_idx = marker_idx,
    lambda_nonzero_idx = lambda_nonzero_idx,
    theta = theta,
    # experimental
    theta_bounds = theta_bounds,
    force_pd = force_pd,
    zero_after_efa = zero_after_efa,
    quadprog = quadprog,
    psi_mapping = psi_mapping,
    #
    nobs = lavsamplestats@ntotal
  )
  mm_lambda <- out$lambda
  mm_theta <- diag(out$theta, nvar)
  mm_psi <- out$psi

  # store matrices in lavmodel@GLIST
  lavmodel@GLIST$lambda <- mm_lambda
  lavmodel@GLIST$theta <- mm_theta
  lavmodel@GLIST$psi <- mm_psi

  # extract free parameters only
  x <- lav_model_get_parameters(lavmodel)

  # apply bounds (if any)
  if (!is.null(lavpartable$lower)) {
    lower_x <- lavpartable$lower[lavpartable$free > 0]
    too_small_idx <- which(x < lower_x)
    if (length(too_small_idx) > 0L) {
      x[too_small_idx] <- lower_x[too_small_idx]
    }
  }
  if (!is.null(lavpartable$upper)) {
    upper_x <- lavpartable$upper[lavpartable$free > 0]
    too_large_idx <- which(x > upper_x)
    if (length(too_large_idx) > 0L) {
      x[too_large_idx] <- upper_x[too_large_idx]
    }
  }

  x
}
