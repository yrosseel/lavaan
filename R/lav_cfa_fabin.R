# FABIN = factor analysis by instrumental variables
# Hagglund 1982 (efa), 1986 (cfa)

lav_cfa_fabin2 <- function(S, marker.idx = NULL, lambda.nonzero.idx = NULL) {
  nvar <- ncol(S)
  nfac <- length(marker.idx)

  # overview of free/fixed
  LAMBDA <- matrix(0, nvar, nfac)
  LAMBDA[lambda.nonzero.idx] <- -1L

  lambda <- matrix(0, nvar, nfac)
  for (i in 1:nvar) {
    if (i %in% marker.idx) {
      lambda[i, marker.idx == i] <- 1.0
      next
    }
    free.idx <- LAMBDA[i, ] == -1L
    idx3 <- (1:nvar)[-c(i, marker.idx)]
    s23 <- S[i, idx3]
    fac.idx <- marker.idx[free.idx]

    if (length(fac.idx) == 1L) { # most common scenario in CFA
      S31 <- S13 <- S[idx3, fac.idx]
      lambda[i, free.idx] <- sum(s23 * S31) / sum(S13 * S13)
    } else {
      S31 <- S[idx3, fac.idx, drop = FALSE]
      S13 <- S[fac.idx, idx3, drop = FALSE]
      lambda[i, free.idx] <- solve(S13 %*% S31, drop(s23 %*% S31))
    }
  }

  lambda
}

lav_cfa_fabin3 <- function(S, marker.idx = NULL, lambda.nonzero.idx = NULL) {
  nvar <- ncol(S)
  nfac <- length(marker.idx)

  # overview of free/fixed
  LAMBDA <- matrix(0, nvar, nfac)
  LAMBDA[lambda.nonzero.idx] <- -1L

  S33.inv <- try(solve(S[-marker.idx, -marker.idx, drop = FALSE]),
    silent = TRUE
  )
  if (inherits(S33.inv, "try-error")) {
    lav_msg_warn(gettext("fabin3 failed; switching to fabin2"))
    return(lav_cfa_fabin2(
      S = S, marker.idx = marker.idx,
      lambda.nonzero.idx = lambda.nonzero.idx
    ))
  }

  lambda <- matrix(0, nvar, nfac)
  rm3.idx <- 0L
  for (i in 1:nvar) {
    if (i %in% marker.idx) {
      lambda[i, marker.idx == i] <- 1.0
      next
    }
    free.idx <- LAMBDA[i, ] == -1L
    idx3 <- (1:nvar)[-c(i, marker.idx)]
    S33 <- S[idx3, idx3, drop = FALSE]
    s23 <- S[i, idx3]
    fac.idx <- marker.idx[free.idx]
    rm3.idx <- rm3.idx + 1L
    # update inverse
    s33.inv <- lav_matrix_symmetric_inverse_update(
      S.inv = S33.inv,
      rm.idx = rm3.idx
    )

    if (length(fac.idx) == 1L) { # most common scenario in CFA
      S31 <- S13 <- S[idx3, fac.idx]
      tmp <- s33.inv %*% S31 # or colSums(s33.inv * S31)
      lambda[i, free.idx] <- sum(s23 * tmp) / sum(S13 * tmp)
    } else {
      S31 <- S[idx3, fac.idx, drop = FALSE]
      S13 <- S[fac.idx, idx3, drop = FALSE]
      tmp <- s33.inv %*% S31
      # lambda[i, free.idx] <- ( s23 %*% solve(S33) %*% S31 %*%
      #                          solve(S13 %*% solve(S33) %*% S31) )
      lambda[i, free.idx] <- solve(S13 %*% tmp, drop(s23 %*% tmp))
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
  sample.cov <- lavsamplestats@cov[[b]]
  nvar <- nrow(sample.cov)
  lv.names <- lavpta$vnames$lv.regular[[b]]
  nfac <- length(lv.names)
  marker.idx <- lavpta$vidx$lv.marker[[b]]
  lambda.idx <- which(names(lavmodel@GLIST) == "lambda")
  lambda.nonzero.idx <- lavmodel@m.free.idx[[lambda.idx]]
  # only diagonal THETA for now...
  # because if we have correlated residuals, we should remove the
  # corresponding variables as instruments before we estimate lambda...
  # (see MIIV)
  theta.idx <- which(names(lavmodel@GLIST) == "theta") # usually '2'
  m.theta <- lavmodel@m.free.idx[[theta.idx]]
  nondiag.idx <- m.theta[!m.theta %in% lav_matrix_diag_idx(nvar)]
  if (length(nondiag.idx) > 0L) {
    lav_msg_warn(gettext(
      "this implementation of FABIN does not handle correlated residuals yet!"
      ))
  }

  # 1. estimate LAMBDA
  if (lavoptions$estimator == "FABIN2") {
    LAMBDA <- lav_cfa_fabin2(
      S = sample.cov, marker.idx = marker.idx,
      lambda.nonzero.idx = lambda.nonzero.idx
    )
  } else {
    LAMBDA <- lav_cfa_fabin3(
      S = sample.cov, marker.idx = marker.idx,
      lambda.nonzero.idx = lambda.nonzero.idx
    )
  }

  # 2. simple ULS method to get THETA and PSI (for now)
  GLS.flag <- FALSE
  psi.mapping.ML.flag <- FALSE
  if (!is.null(lavoptions$estimator.args$thetapsi.method) &&
    lavoptions$estimator.args$thetapsi.method %in% c("GLS", "GLS.ML")) {
    GLS.flag <- TRUE
  }
  if (!is.null(lavoptions$estimator.args$thetapsi.method) &&
    lavoptions$estimator.args$thetapsi.method %in% c("ULS.ML", "GLS.ML")) {
    psi.mapping.ML.flag <- TRUE
  }
  out <- lav_cfa_lambda2thetapsi(
    lambda = LAMBDA, S = sample.cov,
    S.inv = lavsamplestats@icov[[b]],
    GLS = GLS.flag,
    psi.mapping.ML = psi.mapping.ML.flag,
    nobs = lavsamplestats@ntotal
  )
  THETA <- diag(out$theta)
  PSI <- out$psi

  # 3. correlated residuals (if any) are just the difference between
  #    Sigma and S
  # if(length(nondiag.idx) > 0L) {
  #    Sigma <- LAMBDA %*% PSI %*% t(LAMBDA) + THETA
  #    THETA[nondiag.idx] <- (sample.cov - Sigma)[nondiag.idx]
  # }

  # store matrices in lavmodel@GLIST
  lavmodel@GLIST$lambda <- LAMBDA
  lavmodel@GLIST$theta <- THETA
  lavmodel@GLIST$psi <- PSI

  # extract free parameters only
  x <- lav_model_get_parameters(lavmodel)

  # apply bounds (if any)
  if (!is.null(lavpartable$lower)) {
    lower.x <- lavpartable$lower[lavpartable$free > 0]
    too.small.idx <- which(x < lower.x)
    if (length(too.small.idx) > 0L) {
      x[too.small.idx] <- lower.x[too.small.idx]
    }
  }
  if (!is.null(lavpartable$upper)) {
    upper.x <- lavpartable$upper[lavpartable$free > 0]
    too.large.idx <- which(x > upper.x)
    if (length(too.large.idx) > 0L) {
      x[too.large.idx] <- upper.x[too.large.idx]
    }
  }

  x
}
