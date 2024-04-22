# the univariate (weighted) linear model

# - scores/gradient/hessian
# - including the residual variance!

# YR - 30 Dec 2019 (replacing the old lav_ols.R routines)

lav_uvreg_fit <- function(y = NULL,
                          X = NULL,
                          wt = NULL,
                          optim.method = "nlminb",
                          control = list(),
                          output = "list") {
  # check weights
  if (is.null(wt)) {
    wt <- rep(1, length(y))
  } else {
    if (length(y) != length(wt)) {
      lav_msg_stop(gettext("length y is not the same as length wt"))
    }
    if (any(wt < 0)) {
      lav_msg_stop(gettext("all weights should be positive"))
    }
  }

  # optim.method
  minObjective <- lav_uvreg_min_objective
  minGradient <- lav_uvreg_min_gradient
  minHessian <- lav_uvreg_min_hessian
  if (optim.method == "nlminb" || optim.method == "nlminb2") {
    # nothing to do
  } else if (optim.method == "nlminb0") {
    minGradient <- minHessian <- NULL
  } else if (optim.method == "nlminb1") {
    minHessian <- NULL
  }

  # create cache environment
  cache <- lav_uvreg_init_cache(y = y, X = X, wt = wt)

  # optimize -- only changes from defaults
  control.nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control.nlminb <- modifyList(control.nlminb, control)

  optim <- nlminb(
    start = cache$theta, objective = minObjective,
    gradient = minGradient, hessian = minHessian,
    control = control.nlminb, cache = cache
  )

  if (output == "cache") {
    return(cache)
  }

  # return results as a list (to be compatible with lav_polychor.R)
  out <- list(
    theta = optim$par,
    nexo = cache$nexo,
    int.idx = cache$int.idx,
    slope.idx = cache$slope.idx,
    beta.idx = cache$beta.idx,
    var.idx = cache$var.idx,
    y = cache$y,
    wt = cache$wt,
    X = cache$X1[, -1L, drop = FALSE],
    yhat = cache$yhat
  )
}

# prepare cache environment
lav_uvreg_init_cache <- function(y = NULL,
                                 X = NULL,
                                 wt = rep(1, length(y)),
                                 parent = parent.frame()) {
  # y
  y <- as.vector(y)

  # X
  if (is.null(X)) {
    nexo <- 0L
    X1 <- matrix(1, length(y), 1)
  } else {
    X <- unname(X)
    nexo <- ncol(X)
    X1 <- cbind(1, X, deparse.level = 0)

    # new in 0.6-17: check if X is full rank
    if (!anyNA(X)) {
      if (qr(X)$rank < ncol(X)) {
        lav_msg_stop(gettext("matrix of exogenous covariates is rank deficient!
                    (i.e., some x variables contain redundant information)"))
      }
    }
  }

  # nobs
  if (is.null(wt)) {
    N <- length(y)
  } else {
    N <- sum(wt)
  }


  # indices of free parameters
  int.idx <- 1L
  slope.idx <- seq_len(nexo) + 1L
  beta.idx <- c(int.idx, slope.idx)
  var.idx <- 1L + nexo + 1L

  # starting values + crossprod
  if (any(is.na(y)) || any(is.na(X1))) {
    missing.idx <- which(apply(cbind(y, X1), 1, function(x) any(is.na(x))))
    y.tmp <- y[-missing.idx]
    X1.tmp <- X1[-missing.idx, , drop = FALSE]
    wt.tmp <- wt[-missing.idx]
    fit.lm <- stats::lm.wfit(y = y.tmp, x = X1.tmp, w = wt.tmp)
    theta.evar <- sum(fit.lm$residuals * wt.tmp * fit.lm$residuals) / sum(wt.tmp)

    lav_crossprod <- lav_matrix_crossprod
  } else {
    fit.lm <- stats::lm.wfit(y = y, x = X1, w = wt)
    theta.evar <- sum(fit.lm$residuals * wt * fit.lm$residuals) / sum(wt)

    lav_crossprod <- base::crossprod
  }
  theta.beta <- unname(fit.lm$coefficients)
  theta <- c(theta.beta, theta.evar)

  out <- list2env(
    list(
      y = y, X1 = X1, wt = wt, N = N,
      int.idx = int.idx, beta.idx = beta.idx,
      var.idx = var.idx, slope.idx = slope.idx, nexo = nexo,
      lav_crossprod = lav_crossprod,
      theta = theta
    ),
    parent = parent
  )

  out
}

# compute total (log)likelihood
lav_uvreg_loglik <- function(y = NULL,
                             X = NULL,
                             wt = rep(1, length(y)),
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, X = X, wt = wt, output = "cache")
  }
  lav_uvreg_loglik_cache(cache = cache)
}

lav_uvreg_loglik_cache <- function(cache = NULL) {
  with(cache, {
    # free parameters
    beta <- theta[beta.idx]
    evar <- theta[var.idx]

    yhat <- drop(X1 %*% beta)
    logliki <- dnorm(y, mean = yhat, sd = sqrt(evar), log = TRUE)

    # total weighted log-likelihood
    loglik <- sum(wt * logliki, na.rm = TRUE)

    return(loglik)
  })
}

# casewise scores
lav_uvreg_scores <- function(y = NULL,
                             X = NULL,
                             wt = rep(1, length(y)),
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, X = X, wt = wt, output = "cache")
  }
  lav_uvreg_scores_cache(cache = cache)
}

lav_uvreg_scores_cache <- function(cache = NULL) {
  with(cache, {
    res <- y - yhat
    resw <- res * wt
    evar2 <- evar * evar

    scores.beta <- 1 / evar * X1 * resw
    scores.evar <- -wt / (2 * evar) + 1 / (2 * evar2) * res * resw

    return(cbind(scores.beta, scores.evar, deparse.level = 0))
  })
}

# gradient
lav_uvreg_gradient <- function(y = NULL,
                               X = NULL,
                               wt = rep(1, length(y)),
                               cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, X = X, wt = wt, output = "cache")
  }
  lav_uvreg_gradient_cache(cache = cache)
}

lav_uvreg_gradient_cache <- function(cache = NULL) {
  with(cache, {
    res <- y - yhat
    resw <- res * wt
    evar2 <- evar * evar

    dx.beta <- colSums(1 / evar * X1 * resw, na.rm = TRUE)
    dx.var <- sum(-wt / (2 * evar) + 1 / (2 * evar2) * res * resw, na.rm = TRUE)

    return(c(dx.beta, dx.var))
  })
}

# compute total Hessian
lav_uvreg_hessian <- function(y = NULL,
                              X = NULL,
                              wt = rep(1, length(y)),
                              cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, X = X, wt = wt, output = "cache")
  }
  lav_uvreg_hessian_cache(cache = cache)
}

lav_uvreg_hessian_cache <- function(cache = NULL) {
  with(cache, {
    dx2.beta <- -1 / evar * lav_crossprod(X1 * wt, X1)
    dx.beta.var <- -1 / (evar2) * lav_crossprod(X1, resw)

    sq.evar <- sqrt(evar)
    sq.evar6 <- sq.evar * sq.evar * sq.evar * sq.evar * sq.evar * sq.evar
    dx2.var <- (sum(wt, na.rm = TRUE) / (2 * evar2) -
      1 / sq.evar6 * sum(resw * res, na.rm = TRUE))

    Hessian <- rbind(cbind(dx2.beta, dx.beta.var, deparse.level = 0),
      cbind(t(dx.beta.var), dx2.var, deparse.level = 0),
      deparse.level = 0
    )
    return(Hessian)
  })
}

# compute total (log)likelihood, for specific 'x' (nlminb)
lav_uvreg_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  -1 * lav_uvreg_loglik_cache(cache = cache) / cache$N
}

# compute gradient, for specific 'x' (nlminb)
lav_uvreg_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_uvreg_loglik_cache(cache = cache)
  }
  -1 * lav_uvreg_gradient_cache(cache = cache) / cache$N
}

# compute hessian, for specific 'x' (nlminb)
lav_uvreg_min_hessian <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_uvreg_loglik_cache(cache = cache)
    tmp <- lav_uvreg_gradient_cache(cache = cache)
  }
  -1 * lav_uvreg_hessian_cache(cache = cache) / cache$N
}

# update fit object with new parameters
lav_uvreg_update_fit <- function(fit.y = NULL,
                                 evar.new = NULL, beta.new = NULL) {
  if (is.null(evar.new) && is.null(beta.new)) {
    return(fit.y)
  }

  if (!is.null(evar.new)) {
    fit.y$theta[fit.y$var.idx] <- evar.new
  }
  if (!is.null(beta.new)) {
    fit.y$theta[fit.y$beta.idx] <- beta.new
  }

  beta <- fit.y$theta[fit.y$beta.idx]
  X <- fit.y$X
  X1 <- cbind(1, X, deparse.level = 0)

  fit.y$yhat <- drop(X1 %*% beta)

  fit.y
}
