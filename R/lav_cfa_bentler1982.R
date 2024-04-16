# a partial implementation of the Bentler (1982) non-iterative method for CFA
#
# Bentler, P. M. (1982). Confirmatory factor-analysis via noniterative
# estimation - a fast, inexpensive method. Journal of Marketing Research,
# 19(4), 417-424. https://doi.org/10.1177/002224378201900403
#
#
# YR 03 Feb 2023: - first version in lavaan: simple setting only,
#                   no constraints, no 'fixed' (but nonzero) values,
#                   no correlated residuals (ie diagonal-theta only!)
# YR 23 Apr 2023: - quadprog is not needed if we have no (in)equality
#                   constraints

lav_cfa_bentler1982 <- function(S,
                                marker.idx = NULL,
                                lambda.nonzero.idx = NULL,
                                GLS = FALSE,
                                bounds = TRUE,
                                min.reliability.marker = 0.1,
                                quadprog = FALSE,
                                nobs = 20L) { # for cutoff
  # dimensions
  nvar <- ncol(S)
  nfac <- length(marker.idx)

  # lambda structure
  B <- matrix(0, nvar, nfac)
  lambda.marker.idx <- (seq_len(nfac) - 1L) * nvar + marker.idx
  B[lambda.marker.idx] <- 1L
  B[lambda.nonzero.idx] <- 1L

  # partition sample covariance matrix: marker vs non-marker
  S.xx <- S[marker.idx, marker.idx, drop = FALSE]
  S.yx <- S[-marker.idx, marker.idx, drop = FALSE]
  S.xy <- S[marker.idx, -marker.idx, drop = FALSE]
  S.yy <- S[-marker.idx, -marker.idx, drop = FALSE]
  p <- nvar - nfac
  B.y <- B[-marker.idx, , drop = FALSE]

  # check for p = 0?

  # phase 1: initial estimate for Sigma.yx
  Sigma.yx.hat <- S.yx

  # phase 2: using GLS/ULS to obtain PSI and Theta
  if (GLS) {
    W <- try(solve(S.yy), silent = TRUE)
    if (inherits(W, "try-error")) {
      lav_msg_warn(gettext("could not inverte S.yy; switching to ULS"))
      W <- diag(p)
    }
    WS.yx <- W %*% S.yx
    xy.SWS.yx <- crossprod(S.yx, WS.yx)
    G <- WS.yx %*% solve(xy.SWS.yx) %*% t(WS.yx)
  } else {
    Ip <- diag(p)
    xy.SS.yx <- crossprod(S.yx)
    G <- S.yx %*% solve(xy.SS.yx) %*% t(S.yx)
  }

  # only needed if theta.y is not diagonal:
  # q <- 6 # all free
  # # dimension P: q x p*p where q is the number of free elements theta.y
  # theta.fy <- function(x) {
  #     theta.y <- matrix(0, p, p)
  #     # insert 'free' parameters only
  #     diag(theta.y) <- x
  #     lav_matrix_vec(theta.y)
  # }
  # P <- t(numDeriv::jacobian(func = theta.fy, x = rep(1, q)))

  # tmp1 <- P %*% ((W %x% W) - (G %x% G)) %*% t(P)
  # NOTE:
  # if only the 'diagonal' element of Theta are free (as usual), then we
  # can write tmp1 as
  if (GLS) {
    tmp1 <- W * W - G * G
  } else {
    tmp1 <- Ip - G * G
  }

  # only needed if fixed values
  # Theta.F <- matrix(0, p, p) # all free
  # tmp2 <- W %*% (S.yy - Theta.F) %*% W - G %*% (S.yy - Theta.F) %*% G
  if (GLS) {
    tmp2 <- W %*% S.yy %*% W - G %*% S.yy %*% G
  } else {
    tmp2 <- S.yy - G %*% S.yy %*% G
  }

  # Theta.f    <- as.numeric(solve(tmp1) %*% P %*% lav_matrix_vec(tmp2))
  # Note:
  # if only the 'diagonal' element of Theta are free (as usual), then we
  # can write Theta.f as
  Theta.f <- solve(tmp1, diag(tmp2))
  Theta.f.nobounds <- Theta.f # store unbounded Theta.f values

  # ALWAYS apply standard bounds to proceed
  too.small.idx <- which(Theta.f < 0)
  if (length(too.small.idx) > 0L) {
    Theta.f[too.small.idx] <- 0
  }
  too.large.idx <- which(Theta.f > diag(S.yy))
  if (length(too.large.idx) > 0L) {
    Theta.f[too.large.idx] <- diag(S.yy)[too.large.idx] * 1
  }

  # create diagonal matrix with Theta.f elements on diagonal
  Theta.yhat <- diag(Theta.f, p)

  # force (S.yy - Theta.yhat) to be positive definite
  lambda <- try(lav_matrix_symmetric_diff_smallest_root(S.yy, Theta.yhat),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    SminTheta <- S.yy - Theta.yhat # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda.star <- lambda - 1 / (nobs - 1)
      SminTheta <- S.yy - lambda.star * Theta.yhat
    } else {
      SminTheta <- S.yy - Theta.yhat
    }
  }

  # estimate Phi
  if (GLS) {
    tmp1 <- xy.SWS.yx
    tmp2 <- t(WS.yx) %*% SminTheta %*% WS.yx
  } else {
    tmp1 <- xy.SS.yx
    tmp2 <- t(S.yx) %*% SminTheta %*% S.yx
  }
  PSI <- tmp1 %*% solve(tmp2, tmp1)
  PSI.nobounds <- PSI

  # ALWAYS apply bounds to proceed
  lower.bounds.psi <- diag(S.xx) - (1 - min.reliability.marker) * diag(S.xx)
  toolow.idx <- which(diag(PSI) < lower.bounds.psi)
  if (length(toolow.idx) > 0L) {
    diag(PSI)[toolow.idx] <- lower.bounds.psi[toolow.idx]
  }
  too.large.idx <- which(diag(PSI) > diag(S.xx))
  if (length(too.large.idx) > 0L) {
    diag(PSI)[too.large.idx] <- diag(S.xx)[too.large.idx] * 1
  }

  # in addition, force PSI to be PD
  PSI <- lav_matrix_symmetric_force_pd(PSI, tol = 1e-04)

  # residual variances markers
  Theta.x <- diag(S.xx - PSI)

  # create theta vector
  theta.nobounds <- numeric(nvar)
  theta.nobounds[marker.idx] <- Theta.x
  theta.nobounds[-marker.idx] <- Theta.f.nobounds

  # compute LAMBDA for non-marker items

  if (quadprog) {
    # only really needed if we need to impose (in)equality constraints
    # (TODO)
    Dmat <- lav_matrix_bdiag(rep(list(PSI), p))
    dvec <- as.vector(t(S.yx))
    eq.idx <- which(t(B.y) != 1) # these must be zero (row-wise!)
    Rmat <- diag(nrow(Dmat))[eq.idx, , drop = FALSE]
    bvec <- rep(0, length(eq.idx)) # optional, 0=default
    out <- try(quadprog::solve.QP(
      Dmat = Dmat, dvec = dvec, Amat = t(Rmat),
      meq = length(eq.idx), bvec = bvec
    ), silent = TRUE)
    if (inherits(out, "try-error")) {
      lav_msg_warn(gettext("solve.QP failed to find a solution"))
      Lambda <- matrix(0, nvar, nfac)
      Lambda[marker.idx, ] <- diag(nfac)
      Lambda[lambda.nonzero.idx] <- as.numeric(NA)
      Theta <- numeric(nvar)
      Theta[marker.idx] <- Theta.x
      Theta[-marker.idx] <- Theta.f
      Psi <- PSI
      return(list(
        lambda = Lambda, theta = theta.nobounds,
        psi = PSI.nobounds
      ))
    } else {
      LAMBDA.y <- matrix(out$solution,
        nrow = p, ncol = nfac,
        byrow = TRUE
      )
      # zap almost zero elements
      LAMBDA.y[abs(LAMBDA.y) < sqrt(.Machine$double.eps)] <- 0
    }
  } else {
    # simple version
    LAMBDA.y <- t(t(S.yx) / diag(PSI)) * B.y
  }


  # assemble matrices
  LAMBDA <- matrix(0, nvar, nfac)
  LAMBDA[marker.idx, ] <- diag(nfac)
  LAMBDA[-marker.idx, ] <- LAMBDA.y

  list(lambda = LAMBDA, theta = theta.nobounds, psi = PSI.nobounds)
}

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_cfa_bentler1982_internal <- function(lavobject = NULL, # convenience
                                         # internal slot
                                         lavmodel = NULL,
                                         lavsamplestats = NULL,
                                         lavpartable = NULL,
                                         lavdata = NULL,
                                         lavoptions = NULL,
                                         GLS = TRUE,
                                         min.reliability.marker = 0.1,
                                         quadprog = FALSE,
                                         nobs = 20L) {
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
  # no structural part!
  if (any(lavpartable$op == "~")) {
    lav_msg_stop(gettext("bentler1982 estimator only available for CFA models"))
  }
  # no BETA matrix! (i.e., no higher-order factors)
  if (!is.null(lavmodel@GLIST$beta)) {
    lav_msg_stop(gettext("bentler1982 estimator not available
                 for models that require a BETA matrix"))
  }
  # no std.lv = TRUE for now
  if (lavoptions$std.lv) {
    lav_msg_stop(gettext(
      "bentler1982 estimator not available if std.lv = TRUE"))
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
      "this implementation of FABIN does not handle correlated residuals yet!"))
  }

  if (!missing(GLS)) {
    GLS.flag <- GLS
  } else {
    GLS.flag <- FALSE
    if (!is.null(lavoptions$estimator.args$GLS) &&
      lavoptions$estimator.args$GLS) {
      GLS.flag <- TRUE
    }
  }

  if (missing(quadprog) &&
    !is.null(lavoptions$estimator.args$quadprog)) {
    quadprog <- lavoptions$estimator.args$quadprog
  }

  # run bentler1982 non-iterative CFA algorithm
  out <- lav_cfa_bentler1982(
    S = sample.cov, marker.idx = marker.idx,
    lambda.nonzero.idx = lambda.nonzero.idx,
    GLS = GLS.flag,
    min.reliability.marker = 0.1,
    quadprog = quadprog,
    nobs = lavsamplestats@ntotal
  )
  LAMBDA <- out$lambda
  THETA <- diag(out$theta, nvar)
  PSI <- out$psi

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
