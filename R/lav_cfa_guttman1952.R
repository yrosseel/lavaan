# the 'multiple group' method as described in Guttman, 1952
#
# Guttman, L. (1952). Multiple group methods for common-factor analysis,
# their basis, computation, and interpretation. Psychometrika, 17(2) 209--222
#
# YR 02 Feb 2023: - first version in lavaan, using quadprog (no std.lv yet)

lav_cfa_guttman1952 <- function(S,
                                marker.idx = NULL,
                                lambda.nonzero.idx = NULL,
                                theta = NULL, # vector!
                                theta.bounds = FALSE,
                                force.pd = FALSE,
                                zero.after.efa = FALSE,
                                quadprog = FALSE,
                                psi.mapping = FALSE,
                                nobs = 20L) { # for cutoff
  # dimensions
  nvar <- ncol(S)
  nfac <- length(marker.idx)
  stopifnot(length(theta) == nvar)

  # overview of lambda structure
  B <- matrix(0, nvar, nfac)
  lambda.marker.idx <- (seq_len(nfac) - 1L) * nvar + marker.idx
  B[lambda.marker.idx] <- 1L
  B[lambda.nonzero.idx] <- 1L

  # if we wish to keep SminTheta PD, we must keep theta within bounds
  if (force.pd) {
    theta.bounds <- TRUE
  }
  if (psi.mapping) {
    theta.bounds <- TRUE
    force.pd <- TRUE
  }

  # do we first 'clip' the theta values so they are within standard bounds?
  # (Question: do we need the 0.01 and 0.99 multipliers?)
  diagS <- diag(S)
  if (theta.bounds) {
    # lower bound
    lower.bound <- diagS * 0 # * 0.01
    too.small.idx <- which(theta < lower.bound)
    if (length(too.small.idx) > 0L) {
      theta[too.small.idx] <- lower.bound[too.small.idx]
    }

    # upper bound
    upper.bound <- diagS * 1 # * 0.99
    too.large.idx <- which(theta > upper.bound)
    if (length(too.large.idx) > 0L) {
      theta[too.large.idx] <- upper.bound[too.large.idx]
    }
  }

  # compute SminTheta: S where we replace diagonal with 'communalities'
  diag.theta <- diag(theta, nvar)
  SminTheta <- S - diag.theta
  if (force.pd) {
    lambda <- try(lav_matrix_symmetric_diff_smallest_root(S, diag.theta),
      silent = TRUE
    )
    if (inherits(lambda, "try-error")) {
      lav_msg_warn(gettext("failed to compute lambda"))
      SminTheta <- S - diag.theta # and hope for the best
    } else {
      cutoff <- 1 + 1 / (nobs - 1)
      if (lambda < cutoff) {
        lambda.star <- lambda - 1 / (nobs - 1)
        SminTheta <- S - lambda.star * diag.theta
      } else {
        SminTheta <- S - diag.theta
      }
    }
  } else {
    # at least we force the diagonal elements of SminTheta to be nonnegative
    lower.bound <- diagS * 0.001
    too.small.idx <- which(diag(SminTheta) < lower.bound)
    if (length(too.small.idx) > 0L) {
      diag(SminTheta)[too.small.idx] <- lower.bound[too.small.idx]
    }
  }

  # compute covariances among 1) (corrected) variables, and
  #                           2) (corrected) sum-scores
  YS.COV <- SminTheta %*% B

  # compute covariance matrix of corrected sum-scores
  # SS.COV <- t(B) %*% SminTheta %*% B
  SS.COV <- crossprod(B, YS.COV)

  # scaling factors
  # D.inv.sqrt <- diag(1/sqrt(diag(SS.COV)))
  d.inv.sqrt <- 1 / sqrt(diag(SS.COV))

  # factor correlation matrix
  # PHI <- D.inv.sqrt %*% SS.COV %*% D.inv.sqrt
  PHI <- t(SS.COV * d.inv.sqrt) * d.inv.sqrt

  # factor *structure* matrix
  # (covariances corrected Y & corrected normalized sum-scores)
  # YS.COR <- YS.COV %*% D.inv.sqrt
  YS.COR <- t(YS.COV) * d.inv.sqrt # transposed!

  if (zero.after.efa) {
    # we initially assume a saturated LAMBDA (like EFA)
    # then, we just fix the zero-elements to zero

    LAMBDA <- t(solve(PHI, YS.COR)) # = unconstrained EFA version
    # force zeroes
    LAMBDA <- LAMBDA * B
  } else if (quadprog) {
    # constained version using quadprog
    # only useful if (in)equality constraints are needed (TODo)

    # PHI MUST be positive-definite
    PHI <- cov2cor(lav_matrix_symmetric_force_pd(PHI,
      tol = 1e-04
    )) # option?
    Dmat <- lav_matrix_bdiag(rep(list(PHI), nvar))
    dvec <- as.vector(YS.COR)
    eq.idx <- which(t(B) != 1) # these must be zero (row-wise!)
    Rmat <- diag(nrow(Dmat))[eq.idx, , drop = FALSE]
    bvec <- rep(0, length(eq.idx)) # optional, 0=default
    out <- try(quadprog::solve.QP(
      Dmat = Dmat, dvec = dvec, Amat = t(Rmat),
      meq = length(eq.idx), bvec = bvec
    ), silent = TRUE)
    if (inherits(out, "try-error")) {
      lav_msg_warn(gettext("solve.QP failed to find a solution"))
      Lambda <- B
      Lambda[lambda.nonzero.idx] <- as.numeric(NA)
      Theta <- diag(rep(as.numeric(NA), nvar), nvar)
      Psi <- matrix(as.numeric(NA), nfac, nfac)
      return(list(lambda = Lambda, theta = Theta, psi = Psi))
    } else {
      LAMBDA <- matrix(out$solution, nrow = nvar, ncol = nfac, byrow = TRUE)
      # zap almost zero elements
      LAMBDA[abs(LAMBDA) < sqrt(.Machine$double.eps)] <- 0
    }
  } else {
    # default, if no (in)equality constraints

    YS.COR0 <- YS.COR
    YS.COR0[t(B) != 1] <- 0

    LAMBDA <- t(YS.COR0)
  }

  # rescale LAMBDA, so that 'marker' indicator == 1
  marker.lambda <- LAMBDA[lambda.marker.idx]
  Lambda <- t(t(LAMBDA) * (1 / marker.lambda))

  # rescale PHI, covariance metric
  Psi <- t(PHI * marker.lambda) * marker.lambda

  # redo psi using ML mapping function?
  if (psi.mapping) {
    Ti <- 1 / theta
    zero.theta.idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero.theta.idx) > 0L) {
      Ti[zero.theta.idx] <- 1
    }

    # ML mapping function
    M <- solve(t(Lambda) %*% diag(Ti, nvar) %*% Lambda) %*% t(Lambda) %*% diag(Ti, nvar)
    Psi <- M %*% SminTheta %*% t(M)
  }

  list(lambda = Lambda, theta = theta, psi = Psi)
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
                                         theta.bounds = TRUE,
                                         force.pd = TRUE,
                                         zero.after.efa = FALSE,
                                         quadprog = FALSE,
                                         psi.mapping = TRUE) {
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

  if (missing(psi.mapping) &&
    !is.null(lavoptions$estimator.args$psi.mapping)) {
    psi.mapping <- lavoptions$estimator.args$psi.mapping
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

  # 1. obtain estimate for (diagonal elements of) THETA
  #    for now we use Spearman per factor
  B <- matrix(0, nvar, nfac)
  lambda.marker.idx <- (seq_len(nfac) - 1L) * nvar + marker.idx
  B[lambda.marker.idx] <- 1L
  B[lambda.nonzero.idx] <- 1L
  theta <- numeric(nvar)
  for (f in seq_len(nfac)) {
    ov.idx <- which(B[, f] == 1L)
    S.fac <- sample.cov[ov.idx, ov.idx, drop = FALSE]
    theta[ov.idx] <- lav_cfa_theta_spearman(S.fac, bounds = "wide")
  }

  # 2. run Guttman1952 'Multiple Groups' algorithm
  out <- lav_cfa_guttman1952(
    S = sample.cov, marker.idx = marker.idx,
    lambda.nonzero.idx = lambda.nonzero.idx,
    theta = theta,
    # experimental
    theta.bounds = theta.bounds,
    force.pd = force.pd,
    zero.after.efa = zero.after.efa,
    quadprog = quadprog,
    psi.mapping = psi.mapping,
    #
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
