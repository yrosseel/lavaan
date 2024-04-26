# James-Stein estimator
#
# Burghgraeve, E., De Neve, J., & Rosseel, Y. (2021). Estimating structural
# equation models using James-Stein type shrinkage estimators. Psychometrika,
# 86(1), 96-130.
#
# YR 08 Feb 2023: - first version in lavaan, cfa only (for now)

lav_cfa_jamesstein <- function(S,
                               Y = NULL, # raw data
                               marker.idx = NULL,
                               lambda.nonzero.idx = NULL,
                               theta = NULL, # vector!
                               theta.bounds = TRUE,
                               aggregated = FALSE) { # aggregated?
  # dimensions
  nvar <- ncol(S)
  nfac <- length(marker.idx)
  stopifnot(length(theta) == nvar)
  N <- nrow(Y)
  stopifnot(ncol(Y) == nvar)

  # overview of lambda structure
  B <- LAMBDA <- B.nomarker <- matrix(0, nvar, nfac)
  lambda.marker.idx <- (seq_len(nfac) - 1L) * nvar + marker.idx
  B[lambda.marker.idx] <- LAMBDA[lambda.marker.idx] <- 1L
  B[lambda.nonzero.idx] <- B.nomarker[lambda.nonzero.idx] <- 1L

  # Nu
  NU <- numeric(nvar)

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

  # compute conditional expectation conditional on the scaling indicator
  E.JS1 <- lav_cfa_jamesstein_ce(
    Y = Y, marker.idx = marker.idx,
    resvars.markers = theta[marker.idx]
  )

  # compute LAMBDA
  for (f in seq_len(nfac)) {
    nomarker.idx <- which(B.nomarker[, f] == 1)
    Y.nomarker.f <- Y[, nomarker.idx, drop = FALSE]

    # regress no.marker.idx data on E(\eta|Y)
    fit <- lm(Y.nomarker.f ~ E.JS1[, f, drop = FALSE])

    # extract 'lambda' values
    LAMBDA[nomarker.idx, f] <- drop(coef(fit)[-1, ])

    # (optional) extract means
    # NU[nomarker.idx] <- drop(coef(fit)[1,])

    if (aggregated) {
      # local copy of 'scaling' LAMBDA
      LAMBDA.scaling <- LAMBDA

      J <- length(nomarker.idx)
      for (j in seq_len(J)) {
        # data without this indicator
        j.idx <- nomarker.idx[j]
        no.j.idx <- c(marker.idx[f], nomarker.idx[-j])
        Y.agg <- Y[, no.j.idx, drop = FALSE]
        Y.j <- Y[, j.idx, drop = FALSE]

        # retrieve estimated values scaling JS
        lambda.JS.scaling <- LAMBDA.scaling[no.j.idx, f, drop = FALSE]

        # optimize the weights
        starting.weights <- rep(1 / J, times = J)
        w <- optim(
          par = starting.weights,
          fn = lav_cfa_jamesstein_rel,
          data = Y.agg,
          resvars = theta[no.j.idx]
        )$par

        # make sure the weights sum up to 1
        w.optim <- w / sum(w)

        # compute aggregated indicator using the optimal weights
        y_agg <- t(t(w.optim) %*% t(Y.agg))

        # compute error variance of the aggregated indicator
        var_eps_agg <- drop(t(w.optim) %*%
          diag(theta[no.j.idx], nrow = length(no.j.idx)) %*% w.optim)

        # compute conditional expectation using aggregated indicator
        tmp <- lav_cfa_jamesstein_ce(
          Y = y_agg, marker.idx = 1L,
          resvars.markers = var_eps_agg
        )
        CE_agg <- tmp / drop(w.optim %*% lambda.JS.scaling)

        # compute factor loading
        fit <- lm(Y.j ~ CE_agg)
        LAMBDA[j.idx, f] <- drop(coef(fit)[-1])

        # (optional) extract means
        # NU[j.idx] <- drop(coef(fit)[1,])
      } # j
    } # aggregate
  } # f

  out <- list(lambda = LAMBDA, nu = NU)
}

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_cfa_jamesstein_internal <- function(lavobject = NULL, # convenience
                                        # internal slot
                                        lavmodel = NULL,
                                        lavsamplestats = NULL,
                                        lavpartable = NULL,
                                        lavdata = NULL,
                                        lavoptions = NULL,
                                        theta.bounds = TRUE) {
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
    lav_msg_stop(gettext(
      "JS(A) estimator only available for CFA models (for now)"))
  }
  # no BETA matrix! (i.e., no higher-order factors)
  if (!is.null(lavmodel@GLIST$beta)) {
    lav_msg_stop(gettext(
      "JS(A) estimator not available for models that require a BETA matrix"))
  }
  # no std.lv = TRUE for now
  if (lavoptions$std.lv) {
    lav_msg_stop(gettext("S(A) estimator not available if std.lv = TRUE"))
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
  theta.idx <- which(names(lavmodel@GLIST) == "theta") # usually '2'
  m.theta <- lavmodel@m.free.idx[[theta.idx]]
  nondiag.idx <- m.theta[!m.theta %in% lav_matrix_diag_idx(nvar)]
  if (length(nondiag.idx) > 0L) {
    lav_msg_warn(gettext(
      "this implementation of JS/JSA does not handle correlated residuals yet!"
      ))
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
  THETA <- diag(theta, nrow = nvar)

  # 2. run James-Stein algorithm
  Y <- lavdata@X[[1]] # raw data
  aggregated <- FALSE
  if (lavoptions$estimator == "JSA") {
    aggregated <- TRUE
  }
  out <- lav_cfa_jamesstein(
    S = sample.cov, Y = Y, marker.idx = marker.idx,
    lambda.nonzero.idx = lambda.nonzero.idx,
    theta = theta,
    # experimental
    theta.bounds = theta.bounds,
    #
    aggregated = aggregated
  )
  LAMBDA <- out$lambda

  # 3. PSI
  PSI <- lav_cfa_lambdatheta2psi(
    lambda = LAMBDA, theta = theta,
    S = sample.cov, mapping = "ML"
  )

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


# Conditional expectation (Section 2.1, eq. 10)
lav_cfa_jamesstein_ce <- function(Y = NULL,
                                  marker.idx = NULL,
                                  resvars.markers = NULL) {
  Y <- as.matrix(Y)

  # sample size
  N <- nrow(Y)
  N1 <- N - 1
  N3 <- N - 3

  # markers only
  Y.marker <- Y[, marker.idx, drop = FALSE]

  # means and variances
  MEAN <- colMeans(Y.marker, na.rm = TRUE)
  VAR <- apply(Y.marker, 2, var, na.rm = TRUE)

  # 1 - R per maker
  oneminR <- N3 * resvars.markers / (N1 * VAR)

  # R per marker
  R <- 1 - oneminR

  # create E(\eta | Y)
  E.eta.cond.Y <- t(t(Y.marker) * R + oneminR * MEAN)

  E.eta.cond.Y
}

# Reliability function used to obtain the weights (Section 4, Aggregation)
lav_cfa_jamesstein_rel <- function(w = NULL, data = NULL, resvars = NULL) {
  # construct weight vector
  w <- matrix(w, ncol = 1)

  # construct aggregated indicator: y_agg = t(w) %*% y_i
  y_agg <- t(t(w) %*% t(data))

  # calculate variance of aggregated indicator
  var_y_agg <- var(y_agg)

  # calculate error variance of the aggregated indicator
  var_eps_agg <- t(w) %*% diag(resvars) %*% w

  # reliability function to be maximized
  rel <- (var_y_agg - var_eps_agg) %*% solve(var_y_agg)

  # return value
  return(-rel)
}
