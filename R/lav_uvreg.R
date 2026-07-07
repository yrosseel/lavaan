# the univariate (weighted) linear model
#
# - scores/gradient/hessian
# - including the residual variance!
#
# YR 30 Dec 2019 (replacing the old lav_ols.R routines)
# YR/LDW 2026: refactored (see lav_uvbv_common.R for the shared cache/minfns
#              design); fixed a dormant bug in lav_uvreg_update_fit() (yhat
#              lost the slopes; never triggered, as all callers pass
#              intercept-only models)
#
# The parameter vector is theta = c(beta, evar), with beta = c(int, slopes)
# and evar the residual variance; yhat = X1 %*% beta with X1 = [1, X].
#
# fit object (list) contract -- see lav_uvreg_fit():
#   theta, nexo, int_idx, slope_idx, beta_idx, var_idx, y, wt, x, yhat
# The yhat values are those of the *last* evaluation inside nlminb;
# lav_uvreg_update_fit() recomputes them for new parameter values.

lav_uvreg_fit <- function(y = NULL,
                          x = NULL,
                          wt = NULL,
                          optim_method = "nlminb",
                          control = list(),
                          output = "list") {
  # check weights
  wt <- lav_uvbv_wt_check(y = y, wt = wt)

  # optim.method
  minfns <- lav_uvbv_optim_fns(optim_method, lav_uvreg_min_fns())

  # create cache environment
  cache <- lav_uvreg_init_cache(y = y, x = x, wt = wt)

  # optimize -- only changes from defaults
  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control_nlminb <- modifyList(control_nlminb, control)

  optim <- nlminb(
    start = cache$theta, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control_nlminb, cache = cache
  )

  if (output == "cache") {
    return(cache)
  }

  # return results as a list (to be compatible with lav_bvord.R)
  list(
    theta = optim$par,
    nexo = cache$nexo,
    int_idx = cache$int_idx,
    slope_idx = cache$slope_idx,
    beta_idx = cache$beta_idx,
    var_idx = cache$var_idx,
    y = cache$y,
    wt = cache$wt,
    x = cache$x1[, -1L, drop = FALSE],
    yhat = cache$yhat
  )
}

# prepare cache environment
lav_uvreg_init_cache <- function(y = NULL,
                                 x = NULL,
                                 wt = rep(1, length(y)),
                                 parent = parent.frame()) {
  # y
  y <- as.vector(y)

  # X
  if (is.null(x)) {
    nexo <- 0L
    x1 <- matrix(1, length(y), 1)
  } else {
    x <- unname(x)
    nexo <- ncol(x)
    x1 <- cbind(1, x, deparse.level = 0)

    # new in 0.6-17: check if X is full rank
    if (!anyNA(x)) {
      if (qr(x)$rank < ncol(x)) {
        lav_msg_stop(gettext("matrix of exogenous covariates is rank deficient!
                    (i.e., some x variables contain redundant information)"))
      }
    }
  }

  # nobs
  if (is.null(wt)) {
    n <- length(y)
  } else {
    n <- sum(wt)
  }


  # indices of free parameters
  int_idx <- 1L
  slope_idx <- seq_len(nexo) + 1L
  beta_idx <- c(int_idx, slope_idx)
  var_idx <- 1L + nexo + 1L

  # starting values + crossprod
  if (any(is.na(y)) || any(is.na(x1))) {
    missing_idx <- which(apply(cbind(y, x1), 1, function(x) any(is.na(x))))
    y_tmp <- y[-missing_idx]
    x1_tmp <- x1[-missing_idx, , drop = FALSE]
    wt_tmp <- wt[-missing_idx]
    fit_lm <- stats::lm.wfit(y = y_tmp, x = x1_tmp, w = wt_tmp)
    theta_evar <- sum(fit_lm$residuals * wt_tmp * fit_lm$residuals) /
                  sum(wt_tmp)

    lav_crossprod <- lav_mat_crossprod
  } else {
    fit_lm <- stats::lm.wfit(y = y, x = x1, w = wt)
    theta_evar <- sum(fit_lm$residuals * wt * fit_lm$residuals) / sum(wt)

    lav_crossprod <- base::crossprod
  }
  theta_beta <- unname(fit_lm$coefficients)
  theta <- c(theta_beta, theta_evar)

  list2env(
    list(
      y = y, x1 = x1, wt = wt, n = n,
      int_idx = int_idx, beta_idx = beta_idx,
      var_idx = var_idx, slope_idx = slope_idx, nexo = nexo,
      lav_crossprod = lav_crossprod,
      theta = theta
    ),
    parent = parent
  )
}

# compute total (log)likelihood
lav_uvreg_loglik <- function(y = NULL,
                             x = NULL,
                             wt = rep(1, length(y)),
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, x = x, wt = wt, output = "cache")
  }
  lav_uvreg_loglik_cache(cache = cache)
}

lav_uvreg_loglik_cache <- function(cache = NULL) {
  with(cache, {    # nolint start
    # free parameters
    beta <- theta[beta_idx]
    evar <- theta[var_idx]

    yhat <- drop(x1 %*% beta)
    logliki <- dnorm(y, mean = yhat, sd = sqrt(evar), log = TRUE)

    # total weighted log-likelihood
    loglik <- sum(wt * logliki, na.rm = TRUE)

    return(loglik)
  })               # nolint end
}

# casewise scores
lav_uvreg_sc <- function(y = NULL,
                             x = NULL,
                             wt = rep(1, length(y)),
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, x = x, wt = wt, output = "cache")
  }
  lav_uvreg_sc_cache(cache = cache)
}

lav_uvreg_sc_cache <- function(cache = NULL) {
  with(cache, {   # nolint start
    res <- y - yhat
    resw <- res * wt
    evar2 <- evar * evar

    scores_beta <- 1 / evar * x1 * resw
    scores_evar <- -wt / (2 * evar) + 1 / (2 * evar2) * res * resw

    return(cbind(scores_beta, scores_evar, deparse.level = 0))
  })              # nolint end
}

# gradient
lav_uvreg_grad <- function(y = NULL,
                               x = NULL,
                               wt = rep(1, length(y)),
                               cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, x = x, wt = wt, output = "cache")
  }
  lav_uvreg_grad_cache(cache = cache)
}

lav_uvreg_grad_cache <- function(cache = NULL) {
  with(cache, {    # nolint start
    res <- y - yhat
    resw <- res * wt
    evar2 <- evar * evar

    dx_beta <- colSums(1 / evar * x1 * resw, na.rm = TRUE)
    dx_var <- sum(-wt / (2 * evar) + 1 / (2 * evar2) * res * resw, na.rm = TRUE)

    return(c(dx_beta, dx_var))
  })               # nolint end
}

# compute total Hessian
# note: this evaluates at the *current* cache state (as left behind by the
# last nlminb evaluation) -- the beta/evar/yhat/res quantities must already
# be present in the cache
lav_uvreg_hessian <- function(y = NULL,
                              x = NULL,
                              wt = rep(1, length(y)),
                              cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvreg_fit(y = y, x = x, wt = wt, output = "cache")
  }
  lav_uvreg_hessian_cache(cache = cache)
}

lav_uvreg_hessian_cache <- function(cache = NULL) {
  with(cache, {      # nolint start
    dx2_beta <- -1 / evar * lav_crossprod(x1 * wt, x1)
    dx_beta_var <- -1 / (evar2) * lav_crossprod(x1, resw)

    sq_evar <- sqrt(evar)
    sq_evar6 <- sq_evar * sq_evar * sq_evar * sq_evar * sq_evar * sq_evar
    dx2_var <- (sum(wt, na.rm = TRUE) / (2 * evar2) -
      1 / sq_evar6 * sum(resw * res, na.rm = TRUE))

    hessian <- rbind(cbind(dx2_beta, dx_beta_var, deparse.level = 0),
      cbind(t(dx_beta_var), dx2_var, deparse.level = 0),
      deparse.level = 0
    )
    return(hessian)
  })                    # nolint end
}

# nlminb objective/gradient/hessian (see lav_uvbv_common.R)
lav_uvreg_min_fns <- function() {
  lav_uvbv_min_fns(
    logl_fun = lav_uvreg_loglik_cache,
    grad_fun = lav_uvreg_grad_cache,
    hessian_fun = lav_uvreg_hessian_cache
  )
}

# update fit object with new parameters
lav_uvreg_update_fit <- function(fit_y = NULL,
                                 evar_new = NULL, beta_new = NULL) {
  if (is.null(evar_new) && is.null(beta_new)) {
    return(fit_y)
  }

  if (!is.null(evar_new)) {
    fit_y$theta[fit_y$var_idx] <- evar_new
  }
  if (!is.null(beta_new)) {
    fit_y$theta[fit_y$beta_idx] <- beta_new
  }

  beta <- fit_y$theta[fit_y$beta_idx]
  x <- fit_y$x # was fit_y$X (always NULL): yhat lost the slopes
  x1 <- cbind(1, x, deparse.level = 0)

  fit_y$yhat <- drop(x1 %*% beta)

  fit_y
}
