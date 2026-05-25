# the weighted bivariate linear regression model
# YR 14 March 2020 ((replacing the old lav_pearson.R + lav_binorm.R routines)
#
# - bivariate standard normal
# - pearson correlation
# - bivariate linear regression
# - using sampling weights wt


# density of a bivariate __standard__ normal
lav_dbinorm <- function(u, v, rho, force_zero = FALSE) {
  # dirty hack to handle extreme large values for rho
  # note that u, v, and rho are vectorized!
  rho_limit <- 0.9999
  abs_rho <- abs(rho)
  idx <- which(abs_rho > rho_limit)
  if (length(idx) > 0L) {
    rho[idx] <- sign(rho[idx]) * rho_limit
  }

  r <- 1 - rho * rho
  out <- 1 / (2 * pi * sqrt(r)) *
       exp(-0.5 * (u * u - 2 * rho * u * v + v * v) / r)

  # if abs(u) or abs(v) are very large (say, >10), set result equal
  # to exactly zero
  idx <- which(abs(u) > 10 | abs(v) > 10)
  if (length(idx) > 0L && force_zero) {
    out[idx] <- 0
  }

  out
}


# switch between pbivnorm, mnormt, ...
pbinorm <- function(upper_x = NULL, upper_y = NULL, rho = 0.0,
                    lower_x = -Inf, lower_y = -Inf, check = FALSE) {
  pbinorm2(
    upper_x = upper_x, upper_y = upper_y, rho = rho,
    lower_x = lower_x, lower_y = lower_y, check = check
  )
}

# using vectorized version (a la pbivnorm)
pbinorm2 <- function(upper_x = NULL, upper_y = NULL, rho = 0.0,
                     lower_x = -Inf, lower_y = -Inf, check = FALSE) {
  n <- length(upper_x)
  stopifnot(length(upper_y) == n)
  if (n > 1L) {
    if (length(rho) == 1L) {
      rho <- rep(rho, n)
    }
    if (length(lower_x) == 1L) {
      lower_x <- rep(lower_x, n)
    }
    if (length(lower_y) == 1L) {
      lower_y <- rep(lower_y, n)
    }
  }

  upper_only <- all(lower_x == -Inf & lower_y == -Inf)
  if (upper_only) {
    upper_x[upper_x == +Inf] <- exp(10) # better pnorm?
    upper_y[upper_y == +Inf] <- exp(10)
    upper_x[upper_x == -Inf] <- -exp(10)
    upper_y[upper_y == -Inf] <- -exp(10)
    res <- pbivnorm(upper_x, upper_y, rho = rho)
  } else {
    # pbivnorm does not handle -Inf well...
    lower_x[lower_x == -Inf] <- -exp(10)
    lower_y[lower_y == -Inf] <- -exp(10)
    res <- pbivnorm(upper_x, upper_y, rho = rho) -
      pbivnorm(lower_x, upper_y, rho = rho) -
      pbivnorm(upper_x, lower_y, rho = rho) +
      pbivnorm(lower_x, lower_y, rho = rho)
  }

  res
}



# pearson correlation
# if no missing, solution is just cor(Y1,Y2) or cor(e1,e2)
# but if missing, two-step solution is NOT the same as cor(Y1,Y2) or cor(e1,e2)
lav_bvreg_cor_twostep_fit <- function(y1, y2, exo = NULL, wt = NULL,
                                      fit_y1 = NULL, fit_y2 = NULL,
                                      y1_name = NULL, y2_name = NULL,
                                      optim_method = "nlminb1",
                                      # optim.method = "none",
                                      optim_scale = 1,
                                      init_theta = NULL,
                                      control = list()) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvreg_fit(y = y2, x = exo, wt = wt)
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt)

  # the complete case is trivial
  if (!anyNA(fit_y1$y) && !anyNA(fit_y2$y)) {
    return(cache$theta[1L])
  }

  # optim.method
  min_objective <- lav_bvreg_min_objective
  min_gradient <- lav_bvreg_min_grad
  min_hessian <- lav_bvreg_min_hessian
  if (optim_method == "nlminb" || optim_method == "nlminb2") {
    # nothing to do
  } else if (optim_method == "nlminb0") {
    min_gradient <- min_hessian <- NULL
  } else if (optim_method == "nlminb1") {
    min_hessian <- NULL
  } else if (optim_method == "none") {
    return(cache$theta[1L])
  }

  # optimize
  if (is.null(control$trace)) {
    control$trace <- ifelse(lav_verbose(), 1, 0)
  }

  # init theta?
  if (!is.null(init_theta)) {
    start_x <- init_theta
  } else {
    start_x <- cache$theta
  }

  # try 1
  optim <- nlminb(
    start = start_x, objective = min_objective,
    gradient = min_gradient, hessian = min_hessian,
    control = control,
    scale = optim_scale, lower = -0.999, upper = +0.999,
    cache = cache
  )

  # try 2 (scale = 10)
  if (optim$convergence != 0L) {
    optim <- nlminb(
      start = start_x, objective = min_objective,
      gradient = min_gradient, hessian = min_hessian,
      control = control,
      scale = 10, lower = -0.999, upper = +0.999,
      cache = cache
    )
  }

  # try 3 (start = 0, step.min = 0.1)
  if (optim$convergence != 0L) {
    control$step.min <- 0.1
    min_gradient <- lav_bvreg_min_grad
    # try again, with different starting value
    optim <- nlminb(
      start = 0, objective = min_objective,
      gradient = min_gradient, hessian = NULL,
      control = control,
      scale = optim_scale, lower = -0.999, upper = +0.999,
      cache = cache
    )
  }

  # check convergence
  if (optim$convergence != 0L) {
    if (!is.null(y1_name) && !is.null(y2_name)) {
      lav_msg_warn(gettextf(
        "estimation pearson correlation did not converge for variables %1$s
        and %2$s.", y1_name, y2_name))
    } else {
      lav_msg_warn(gettext(
        "estimation pearson correlation(s) did not always converge"))
    }

    # use init (as we always did in < 0.6-6; this is also what Mplus does)
    rho <- start_x
  } else {
    # store result
    rho <- optim$par
  }

  rho
}

# Y1 = linear
# Y2 = linear
lav_bvreg_init_cache <- function(fit_y1 = NULL,
                                 fit_y2 = NULL,
                                 wt = NULL,
                                 scores = FALSE,
                                 parent = parent.frame()) {
  # data
  y1 <- fit_y1$y
  y2 <- fit_y2$y
  exo <- fit_y1$x

  # Y1
  y1c <- y1 - fit_y1$yhat
  evar_y1 <- fit_y1$theta[fit_y1$var_idx]
  sd_y1 <- sqrt(evar_y1)
  eta_y1 <- fit_y1$yhat

  # Y2
  y2c <- y2 - fit_y2$yhat
  evar_y2 <- fit_y2$theta[fit_y1$var_idx]
  sd_y2 <- sqrt(evar_y2)
  eta_y2 <- fit_y2$yhat


  # exo?
  if (is.null(exo)) {
    nexo <- 0L
  } else {
    nexo <- ncol(exo)
  }

  # nobs
  if (is.null(wt)) {
    n <- length(y1)
  } else {
    n <- sum(wt)
  }

  # starting value
  if (fit_y1$nexo > 0L) {
    e1 <- y1 - fit_y1$yhat
    e2 <- y2 - fit_y2$yhat
    if (is.null(wt)) {
      rho_init <- cor(e1, e2, use = "pairwise.complete.obs")
    } else {
      tmp <- na.omit(cbind(e1, e2, wt))
      rho_init <- cov.wt(tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
    }
  } else {
    if (is.null(wt)) {
      rho_init <- cor(y1, y2, use = "pairwise.complete.obs")
    } else {
      tmp <- na.omit(cbind(y1, y2, wt))
      rho_init <- cov.wt(tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
    }
  }

  # sanity check
  if (is.na(rho_init) || abs(rho_init) >= 1.0) {
    rho_init <- 0.0
  }

  # parameter vector
  theta <- rho_init # only

  # different cache if scores or not
  if (scores) {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, n = n,
        y1c = y1c, y2c = y2c, exo = exo,
        evar_y1 = evar_y1, sd_y1 = sd_y1, eta_y1 = eta_y1,
        evar_y2 = evar_y2, sd_y2 = sd_y2, eta_y2 = eta_y2
      ),
      parent = parent
    )
  } else {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, n = n,
        y1c = y1c, y2c = y2c,
        evar_y1 = evar_y1, sd_y1 = sd_y1, eta_y1 = eta_y1,
        evar_y2 = evar_y2, sd_y2 = sd_y2, eta_y2 = eta_y2
      ),
      parent = parent
    )
  }

  out
}


# casewise likelihoods, unweighted!
lav_bvreg_lik_cache <- function(cache = NULL) {
  with(cache, {   # nolint start
    rho <- theta[1L]

    cov_y12 <- rho * sqrt(evar_y1) * sqrt(evar_y2)
    sigma <- matrix(c(evar_y1, cov_y12, cov_y12, evar_y2), 2L, 2L)
    lik <- exp(lav_mvn_loglik_data(
      y = cbind(y1c, y2c), wt = NULL,
      mu = c(0, 0), sigma_1 = sigma,
      casewise = TRUE
    ))
    # catch very small values -- TODO: sqrt() too aggressive?
    lik_toosmall_idx <- which(lik < sqrt(.Machine$double.eps))
    lik[lik_toosmall_idx] <- as.numeric(NA)

    return(lik)
  })             # nolint end
}

lav_bvreg_logl_cache <- function(cache = NULL) {
  with(cache, {       # nolint start
    lik <- lav_bvreg_lik_cache(cache) # unweighted!

    if (!is.null(wt)) {
      logl <- sum(wt * log(lik), na.rm = TRUE)
    } else {
      logl <- sum(log(lik), na.rm = TRUE)
    }

    return(logl)
  })                # nolint end
}

lav_bvreg_grad_cache <- function(cache = NULL) {
  with(cache, {     # nolint start
    rho <- theta[1L]
    r <- (1 - rho * rho)

    sd_y1_y2 <- sd_y1 * sd_y2
    t1 <- (y1c * y2c) / sd_y1_y2
    t2 <- (y1c * y1c) / evar_y1 - (2 * rho * t1) + (y2c * y2c) / evar_y2
    dx <- (rho + t1 - t2 * rho / r) / r

    # to be consistent with (log)lik_cache
    if (length(lik_toosmall_idx) > 0L) {
      dx[lik_toosmall_idx] <- as.numeric(NA)
    }

    if (is.null(wt)) {
      dx_rho <- sum(dx, na.rm = TRUE)
    } else {
      dx_rho <- sum(wt * dx, na.rm = TRUE)
    }

    return(dx_rho)
  })                  # nolint end
}

lav_bvreg_hessian_cache <- function(cache = NULL) {
  with(cache, {       # nolint start
    rho <- theta[1L]

    rho2 <- rho * rho
    r2 <- r * r
    r3 <- r * r * r

    h <- 1 / r - (2 * rho2 * t2) / r3 +
      2 * rho2 * (1 - t2 / r) / r2 + 4 * rho * t1 / r2 - t2 / r2

    # to be consistent with (log)lik_cache
    if (length(lik_toosmall_idx) > 0L) {
      h[lik_toosmall_idx] <- as.numeric(NA)
    }

    if (is.null(wt)) {
      h_1 <- sum(h, na.rm = TRUE)
    } else {
      h_1 <- sum(wt * h, na.rm = TRUE)
    }
    dim(h_1) <- c(1L, 1L) # for nlminb

    return(h_1)
  })               # nolint end
}

# compute total (log)likelihood, for specific 'x' (nlminb)
lav_bvreg_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  -1 * lav_bvreg_logl_cache(cache = cache) / cache$n
}

# compute gradient, for specific 'x' (nlminb)
lav_bvreg_min_grad <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_bvreg_logl_cache(cache = cache)
    rm(tmp)
  }
  -1 * lav_bvreg_grad_cache(cache = cache) / cache$n
}

# compute hessian, for specific 'x' (nlminb)
lav_bvreg_min_hessian <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    tmp <- lav_bvreg_logl_cache(cache = cache)
    tmp <- lav_bvreg_grad_cache(cache = cache)
    rm(tmp)
  }
  -1 * lav_bvreg_hessian_cache(cache = cache) / cache$n
}

# casewise scores - cache
# FIXME: should we also set 'lik.toosmall.idx' cases to NA?
lav_bvreg_cor_sc_cache <- function(cache = NULL) {
  with(cache, {           # nolint start
    rho <- theta[1L]
    r <- (1 - rho * rho)

    # mu.y1
    dx_mu_y1 <- (2 * y1c / evar_y1 - 2 * rho * y2c / (sd_y1 * sd_y2)) / (2 * r)
    if (!is.null(wt)) {
      dx_mu_y1 <- wt * dx_mu_y1
    }

    # mu.y2
    dx_mu_y2 <- -(2 * rho * y1c / (sd_y1 * sd_y2) - 2 * y2c / evar_y2) / (2 * r)
    if (!is.null(wt)) {
      dx_mu_y2 <- wt * dx_mu_y2
    }

    # evar_y1
    dx_var_y1 <- -(0.5 / evar_y1 - ((y1c * y1c) / (evar_y1 * evar_y1) -
      rho * y1c * y2c / (evar_y1 * sd_y1 * sd_y2)) / (2 * r))
    if (!is.null(wt)) {
      dx_var_y1 <- wt * dx_var_y1
    }

    # var.y2
    dx_var_y2 <-
          -(0.5 / evar_y2 + (rho * y1c * y2c / (evar_y2 * sd_y1 * sd_y2) -
           (y2c * y2c) / (evar_y2 * evar_y2)) / (2 * r))
    if (!is.null(wt)) {
      dx_var_y2 <- wt * dx_var_y2
    }

    # sl.y1
    dx_sl_y1 <- NULL
    if (nexo > 0L) {
      dx_sl_y1 <- dx_mu_y1 * exo # weights already included in dx.mu.y1
    }

    # sl.y2
    dx_sl_y2 <- NULL
    if (nexo > 0L) {
      dx_sl_y2 <- dx_mu_y2 * exo # weights already included in dx.mu.y2
    }

    # rho
    z <- (y1c * y1c) / evar_y1 -
         2 * rho * y1c * y2c / (sd_y1 * sd_y2) + (y2c * y2c) / evar_y2
    dx_rho <- rho / r + (y1c * y2c / (sd_y1 * sd_y2 * r) - z * rho / (r * r))
    if (!is.null(wt)) {
      dx_rho <- wt * dx_rho
    }

    out <- list(
      dx_mu_y1 = dx_mu_y1, dx_var_y1 = dx_var_y1,
      dx_mu_y2 = dx_mu_y2, dx_var_y2 = dx_var_y2,
      dx_sl_y1 = dx_sl_y1, dx_sl_y2 = dx_sl_y2,
      dx_rho = dx_rho
    )
    return(out)
  })                       # nolint end
}

# casewise scores
#
# Y1 = linear
# Y2 = linear
lav_bvreg_cor_sc <- function(y1, y2, exo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit_y1 = NULL, fit_y2 = NULL,
                                 evar_y1 = NULL, beta_y1 = NULL,
                                 evar_y2 = NULL, beta_y2 = NULL) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvreg_fit(y = y2, x = exo, wt = wt)
  }

  # user specified parameters
  if (!is.null(evar_y1) || !is.null(beta_y1)) {
    fit_y1 <- lav_uvreg_update_fit(
      fit_y = fit_y1,
      evar_new = evar_y1, beta_new = beta_y1
    )
  }
  if (!is.null(evar_y2) || !is.null(beta_y2)) {
    fit_y2 <- lav_uvreg_update_fit(
      fit_y = fit_y2,
      evar_new = evar_y2, beta_new = beta_y2
    )
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  sc <- lav_bvreg_cor_sc_cache(cache = cache)

  sc
}

# logl - no cache
lav_bvreg_logl <- function(y1, y2, exo = NULL, wt = NULL,
                           rho = NULL,
                           fit_y1 = NULL, fit_y2 = NULL,
                           evar_y1 = NULL, beta_y1 = NULL,
                           evar_y2 = NULL, beta_y2 = NULL) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvreg_fit(y = y2, x = exo, wt = wt)
  }

  # user specified parameters
  if (!is.null(evar_y1) || !is.null(beta_y1)) {
    fit_y1 <- lav_uvreg_update_fit(
      fit_y = fit_y1,
      evar_new = evar_y1, beta_new = beta_y1
    )
  }
  if (!is.null(evar_y2) || !is.null(beta_y2)) {
    fit_y2 <- lav_uvreg_update_fit(
      fit_y = fit_y2,
      evar_new = evar_y2, beta_new = beta_y2
    )
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  lav_bvreg_logl_cache(cache = cache)
}

# lik - no cache
lav_bvreg_lik <- function(y1, y2, exo = NULL, wt = NULL,
                          rho = NULL,
                          fit_y1 = NULL, fit_y2 = NULL,
                          evar_y1 = NULL, beta_y1 = NULL,
                          evar_y2 = NULL, beta_y2 = NULL,
                          .log = FALSE) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvreg_fit(y = y2, x = exo, wt = wt)
  }

  # user specified parameters
  if (!is.null(evar_y1) || !is.null(beta_y1)) {
    fit_y1 <- lav_uvreg_update_fit(
      fit_y = fit_y1,
      evar_new = evar_y1, beta_new = beta_y1
    )
  }
  if (!is.null(evar_y2) || !is.null(beta_y2)) {
    fit_y2 <- lav_uvreg_update_fit(
      fit_y = fit_y2,
      evar_new = evar_y2, beta_new = beta_y2
    )
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  lik <- lav_bvreg_lik_cache(cache = cache)
  if (.log) {
    lik <- log(lik)
  }

  if (!is.null(wt)) {
    if (.log) {
      lik <- wt * lik
    } else {
      tmp <- wt * log(lik)
      lik <- exp(tmp)
    }
  }

  lik
}
