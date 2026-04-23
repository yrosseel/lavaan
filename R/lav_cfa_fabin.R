# FABIN = factor analysis by instrumental variables
# Hagglund 1982 (efa), 1986 (cfa)

lav_cfa_fabin2 <- function(s, marker_idx = NULL, lambda_nonzero_idx = NULL) {
  nvar <- ncol(s)
  nfac <- length(marker_idx)

  # overview of free/fixed
  mm_lambda <- matrix(0, nvar, nfac)
  mm_lambda[lambda_nonzero_idx] <- -1L

  lambda <- matrix(0, nvar, nfac)
  for (i in 1:nvar) {
    if (i %in% marker_idx) {
      lambda[i, marker_idx == i] <- 1.0
      next
    }
    free_idx <- mm_lambda[i, ] == -1L
    idx3 <- (1:nvar)[-c(i, marker_idx)]
    s23 <- s[i, idx3]
    fac_idx <- marker_idx[free_idx]

    if (length(fac_idx) == 1L) { # most common scenario in CFA
      s31 <- s13 <- s[idx3, fac_idx]
      lambda[i, free_idx] <- sum(s23 * s31) / sum(s13 * s13)
    } else {
      s31 <- s[idx3, fac_idx, drop = FALSE]
      s13 <- s[fac_idx, idx3, drop = FALSE]
      lambda[i, free_idx] <- solve(s13 %*% s31, drop(s23 %*% s31))
    }
  }

  lambda
}

lav_cfa_fabin3 <- function(s, marker_idx = NULL, lambda_nonzero_idx = NULL) {
  nvar <- ncol(s)
  nfac <- length(marker_idx)

  # overview of free/fixed
  mm_lambda <- matrix(0, nvar, nfac)
  mm_lambda[lambda_nonzero_idx] <- -1L

  s33_inv_1 <- try(solve(s[-marker_idx, -marker_idx, drop = FALSE]),
    silent = TRUE
  )
  if (inherits(s33_inv_1, "try-error")) {
    lav_msg_warn(gettext("fabin3 failed; switching to fabin2"))
    return(lav_cfa_fabin2(
      s = s, marker_idx = marker_idx,
      lambda_nonzero_idx = lambda_nonzero_idx
    ))
  }

  lambda <- matrix(0, nvar, nfac)
  rm3_idx <- 0L
  for (i in 1:nvar) {
    if (i %in% marker_idx) {
      lambda[i, marker_idx == i] <- 1.0
      next
    }
    free_idx <- mm_lambda[i, ] == -1L
    idx3 <- (1:nvar)[-c(i, marker_idx)]
    # s33 <- s[idx3, idx3, drop = FALSE]
    s23 <- s[i, idx3]
    fac_idx <- marker_idx[free_idx]
    rm3_idx <- rm3_idx + 1L
    # update inverse
    s33_inv <- lav_matrix_symmetric_inverse_update(
      S.inv = s33_inv_1,
      rm.idx = rm3_idx
    )

    if (length(fac_idx) == 1L) { # most common scenario in CFA
      s31 <- s13 <- s[idx3, fac_idx]
      tmp <- s33_inv %*% s31 # or colSums(s33.inv * S31)
      lambda[i, free_idx] <- sum(s23 * tmp) / sum(s13 * tmp)
    } else {
      s31 <- s[idx3, fac_idx, drop = FALSE]
      s13 <- s[fac_idx, idx3, drop = FALSE]
      tmp <- s33_inv %*% s31
      # lambda[i, free.idx] <- ( s23 %*% solve(S33) %*% S31 %*%
      #                          solve(S13 %*% solve(S33) %*% S31) )
      lambda[i, free_idx] <- solve(s13 %*% tmp, drop(s23 %*% tmp))
    }
  }

  lambda
}

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_cfa_fabin_internal <- function(lavmodel = NULL, lavsamplestats = NULL,
                                   lavpartable = NULL,
                                   lavdata = NULL, lavoptions = NULL) {
  lavpta <- lav_partable_attributes(lavpartable)
  lavpartable <- lav_partable_set_cache(lavpartable, lavpta)
  # no structural part!
  if (any(lavpartable$op == "~")) {
    lav_msg_stop(gettext("FABIN estimator only available for CFA models"))
  }
  # no BETA matrix! (i.e., no higher-order factors)
  if (!is.null(lavmodel@GLIST$beta)) {
    lav_msg_stop(gettext(
    "FABIN estimator not available for models that require a BETA matrix"))
  }
  # no std.lv = TRUE for now
  if (lavoptions$std.lv) {
    lav_msg_stop(
      gettext("FABIN estimator not available if std.lv = TRUE"))
  }

  nblocks <- lav_partable_nblocks(lavpartable)
  stopifnot(nblocks == 1L) # for now
  b <- 1L
  sample_cov <- lavsamplestats@cov[[b]]
  nvar <- nrow(sample_cov)
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
      "this implementation of FABIN does not handle correlated residuals yet!"
      ))
  }

  # 1. estimate LAMBDA
  if (lavoptions$estimator == "FABIN2") {
    mm_lambda <- lav_cfa_fabin2(
      s = sample_cov, marker_idx = marker_idx,
      lambda_nonzero_idx = lambda_nonzero_idx
    )
  } else {
    mm_lambda <- lav_cfa_fabin3(
      s = sample_cov, marker_idx = marker_idx,
      lambda_nonzero_idx = lambda_nonzero_idx
    )
  }

  # 2. simple ULS method to get THETA and PSI (for now)
  gls_flag <- FALSE
  psi_mapping_ml_flag <- FALSE
  if (!is.null(lavoptions$estimator.args$thetapsi.method) &&
    lavoptions$estimator.args$thetapsi.method %in% c("GLS", "GLS.ML")) {
    gls_flag <- TRUE
  }
  if (!is.null(lavoptions$estimator.args$thetapsi.method) &&
    lavoptions$estimator.args$thetapsi.method %in% c("ULS.ML", "GLS.ML")) {
    psi_mapping_ml_flag <- TRUE
  }
  out <- lav_cfa_lambda2thetapsi(
    lambda = mm_lambda, s = sample_cov,
    s_inv = lavsamplestats@icov[[b]],
    gls = gls_flag,
    psi_mapping_ml = psi_mapping_ml_flag,
    nobs = lavsamplestats@ntotal
  )
  mm_theta <- diag(out$theta)
  mm_psi <- out$psi

  # 3. correlated residuals (if any) are just the difference between
  #    Sigma and S
  # if(length(nondiag.idx) > 0L) {
  #    Sigma <- LAMBDA %*% PSI %*% t(LAMBDA) + THETA
  #    THETA[nondiag.idx] <- (sample.cov - Sigma)[nondiag.idx]
  # }

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
