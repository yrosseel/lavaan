# the weighted bivariate linear regression model
# YR 14 March 2020 ((replacing the old lav_pearson.R + lav_binorm.R routines)
#
# - bivariate standard normal
# - pearson correlation
# - bivariate linear regression
# - using sampling weights wt


# density of a bivariate __standard__ normal
lav_dbinorm <- dbinorm <- function(u, v, rho, force.zero = FALSE) {
  # dirty hack to handle extreme large values for rho
  # note that u, v, and rho are vectorized!
  RHO.limit <- 0.9999
  abs.rho <- abs(rho)
  idx <- which(abs.rho > RHO.limit)
  if (length(idx) > 0L) {
    rho[idx] <- sign(rho[idx]) * RHO.limit
  }

  R <- 1 - rho * rho
  out <- 1 / (2 * pi * sqrt(R)) * exp(-0.5 * (u * u - 2 * rho * u * v + v * v) / R)

  # if abs(u) or abs(v) are very large (say, >10), set result equal
  # to exactly zero
  idx <- which(abs(u) > 10 | abs(v) > 10)
  if (length(idx) > 0L && force.zero) {
    out[idx] <- 0
  }

  out
}

# partial derivative - rho
lav_dbinorm_drho <- function(u, v, rho) {
  R <- 1 - rho * rho
  dbinorm(u, v, rho) * (u * v * R - rho * (u * u - 2 * rho * u * v + v * v) + rho * R) / (R * R)
}

# partial derivative - u
lav_dbinorm_du <- function(u, v, rho) {
  R <- 1 - rho * rho
  -dbinorm(u, v, rho) * (u - rho * v) / R
}

# partial derivative - v
lav_dbinorm_dv <- function(u, v, rho) {
  R <- 1 - rho * rho
  -dbinorm(u, v, rho) * (v - rho * u) / R
}


# CDF of bivariate standard normal
# function pbinorm(upper.x, upper.y, rho)

# partial derivative pbinorm - upper.x
lav_pbinorm_dupperx <- function(upper.x, upper.y, rho = 0.0) {
  R <- 1 - rho * rho
  dnorm(upper.x) * pnorm((upper.y - rho * upper.x) / sqrt(R))
}

lav_pbinorm_duppery <- function(upper.x, upper.y, rho = 0.0) {
  R <- 1 - rho * rho
  dnorm(upper.y) * pnorm((upper.x - rho * upper.y) / sqrt(R))
}

lav_pbinorm_drho <- function(upper.x, upper.y, rho = 0.0) {
  dbinorm(upper.x, upper.y, rho)
}


# switch between pbivnorm, mnormt, ...
pbinorm <- function(upper.x = NULL, upper.y = NULL, rho = 0.0,
                    lower.x = -Inf, lower.y = -Inf, check = FALSE) {
  pbinorm2(
    upper.x = upper.x, upper.y = upper.y, rho = rho,
    lower.x = lower.x, lower.y = lower.y, check = check
  )
}

# using vectorized version (a la pbivnorm)
pbinorm2 <- function(upper.x = NULL, upper.y = NULL, rho = 0.0,
                     lower.x = -Inf, lower.y = -Inf, check = FALSE) {
  N <- length(upper.x)
  stopifnot(length(upper.y) == N)
  if (N > 1L) {
    if (length(rho) == 1L) {
      rho <- rep(rho, N)
    }
    if (length(lower.x) == 1L) {
      lower.x <- rep(lower.x, N)
    }
    if (length(lower.y) == 1L) {
      lower.y <- rep(lower.y, N)
    }
  }

  upper.only <- all(lower.x == -Inf & lower.y == -Inf)
  if (upper.only) {
    upper.x[upper.x == +Inf] <- exp(10) # better pnorm?
    upper.y[upper.y == +Inf] <- exp(10)
    upper.x[upper.x == -Inf] <- -exp(10)
    upper.y[upper.y == -Inf] <- -exp(10)
    res <- pbivnorm(upper.x, upper.y, rho = rho)
  } else {
    # pbivnorm does not handle -Inf well...
    lower.x[lower.x == -Inf] <- -exp(10)
    lower.y[lower.y == -Inf] <- -exp(10)
    res <- pbivnorm(upper.x, upper.y, rho = rho) -
      pbivnorm(lower.x, upper.y, rho = rho) -
      pbivnorm(upper.x, lower.y, rho = rho) +
      pbivnorm(lower.x, lower.y, rho = rho)
  }

  res
}



# pearson correlation
# if no missing, solution is just cor(Y1,Y2) or cor(e1,e2)
# but if missing, two-step solution is NOT the same as cor(Y1,Y2) or cor(e1,e2)
lav_bvreg_cor_twostep_fit <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                      fit.y1 = NULL, fit.y2 = NULL,
                                      Y1.name = NULL, Y2.name = NULL,
                                      optim.method = "nlminb1",
                                      # optim.method = "none",
                                      optim.scale = 1,
                                      init.theta = NULL,
                                      control = list(),
                                      verbose = FALSE) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvreg_fit(y = Y2, X = eXo, wt = wt)
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt)

  # the complete case is trivial
  if (!anyNA(fit.y1$y) && !anyNA(fit.y2$y)) {
    return(cache$theta[1L])
  }

  # optim.method
  minObjective <- lav_bvreg_min_objective
  minGradient <- lav_bvreg_min_gradient
  minHessian <- lav_bvreg_min_hessian
  if (optim.method == "nlminb" || optim.method == "nlminb2") {
    # nothing to do
  } else if (optim.method == "nlminb0") {
    minGradient <- minHessian <- NULL
  } else if (optim.method == "nlminb1") {
    minHessian <- NULL
  } else if (optim.method == "none") {
    return(cache$theta[1L])
  }

  # optimize
  if (is.null(control$trace)) {
    control$trace <- ifelse(verbose, 1, 0)
  }

  # init theta?
  if (!is.null(init.theta)) {
    start.x <- init.theta
  } else {
    start.x <- cache$theta
  }

  # try 1
  optim <- nlminb(
    start = start.x, objective = minObjective,
    gradient = minGradient, hessian = minHessian,
    control = control,
    scale = optim.scale, lower = -0.999, upper = +0.999,
    cache = cache
  )

  # try 2 (scale = 10)
  if (optim$convergence != 0L) {
    optim <- nlminb(
      start = start.x, objective = minObjective,
      gradient = minGradient, hessian = minHessian,
      control = control,
      scale = 10, lower = -0.999, upper = +0.999,
      cache = cache
    )
  }

  # try 3 (start = 0, step.min = 0.1)
  if (optim$convergence != 0L) {
    control$step.min <- 0.1
    minGradient <- lav_bvreg_min_gradient
    # try again, with different starting value
    optim <- nlminb(
      start = 0, objective = minObjective,
      gradient = minGradient, hessian = NULL,
      control = control,
      scale = optim.scale, lower = -0.999, upper = +0.999,
      cache = cache
    )
  }

  # check convergence
  if (optim$convergence != 0L) {
    if (!is.null(Y1.name) && !is.null(Y2.name)) {
      lav_msg_warn(gettextf(
        "estimation pearson correlation did not converge for variables %1$s
        and %2$s.", Y1.name, Y2.name))
    } else {
      lav_msg_warn(gettext(
        "estimation pearson correlation(s) did not always converge"))
    }

    # use init (as we always did in < 0.6-6; this is also what Mplus does)
    rho <- start.x
  } else {
    # store result
    rho <- optim$par
  }

  rho
}

# Y1 = linear
# Y2 = linear
lav_bvreg_init_cache <- function(fit.y1 = NULL,
                                 fit.y2 = NULL,
                                 wt = NULL,
                                 scores = FALSE,
                                 parent = parent.frame()) {
  # data
  Y1 <- fit.y1$y
  Y2 <- fit.y2$y
  eXo <- fit.y1$X

  # Y1
  Y1c <- Y1 - fit.y1$yhat
  evar.y1 <- fit.y1$theta[fit.y1$var.idx]
  sd.y1 <- sqrt(evar.y1)
  eta.y1 <- fit.y1$yhat

  # Y2
  Y2c <- Y2 - fit.y2$yhat
  evar.y2 <- fit.y2$theta[fit.y1$var.idx]
  sd.y2 <- sqrt(evar.y2)
  eta.y2 <- fit.y2$yhat


  # exo?
  if (is.null(eXo)) {
    nexo <- 0L
  } else {
    nexo <- ncol(eXo)
  }

  # nobs
  if (is.null(wt)) {
    N <- length(Y1)
  } else {
    N <- sum(wt)
  }

  # starting value
  if (fit.y1$nexo > 0L) {
    E1 <- Y1 - fit.y1$yhat
    E2 <- Y2 - fit.y2$yhat
    if (is.null(wt)) {
      rho.init <- cor(E1, E2, use = "pairwise.complete.obs")
    } else {
      tmp <- na.omit(cbind(E1, E2, wt))
      rho.init <- cov.wt(tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
    }
  } else {
    if (is.null(wt)) {
      rho.init <- cor(Y1, Y2, use = "pairwise.complete.obs")
    } else {
      tmp <- na.omit(cbind(Y1, Y2, wt))
      rho.init <- cov.wt(tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
    }
  }

  # sanity check
  if (is.na(rho.init) || abs(rho.init) >= 1.0) {
    rho.init <- 0.0
  }

  # parameter vector
  theta <- rho.init # only

  # different cache if scores or not
  if (scores) {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, N = N,
        Y1c = Y1c, Y2c = Y2c, eXo = eXo,
        evar.y1 = evar.y1, sd.y1 = sd.y1, eta.y1 = eta.y1,
        evar.y2 = evar.y2, sd.y2 = sd.y2, eta.y2 = eta.y2
      ),
      parent = parent
    )
  } else {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, N = N,
        Y1c = Y1c, Y2c = Y2c,
        evar.y1 = evar.y1, sd.y1 = sd.y1, eta.y1 = eta.y1,
        evar.y2 = evar.y2, sd.y2 = sd.y2, eta.y2 = eta.y2
      ),
      parent = parent
    )
  }

  out
}


# casewise likelihoods, unweighted!
lav_bvreg_lik_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    cov.y12 <- rho * sqrt(evar.y1) * sqrt(evar.y2)
    sigma <- matrix(c(evar.y1, cov.y12, cov.y12, evar.y2), 2L, 2L)
    lik <- exp(lav_mvnorm_loglik_data(
      Y = cbind(Y1c, Y2c), wt = NULL,
      Mu = c(0, 0), Sigma = sigma,
      casewise = TRUE
    ))
    # catch very small values
    lik.toosmall.idx <- which(lik < sqrt(.Machine$double.eps))
    lik[lik.toosmall.idx] <- as.numeric(NA)

    return(lik)
  })
}

lav_bvreg_logl_cache <- function(cache = NULL) {
  with(cache, {
    lik <- lav_bvreg_lik_cache(cache) # unweighted!

    if (!is.null(wt)) {
      logl <- sum(wt * log(lik), na.rm = TRUE)
    } else {
      logl <- sum(log(lik), na.rm = TRUE)
    }

    return(logl)
  })
}

lav_bvreg_gradient_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]
    R <- (1 - rho * rho)

    sd.y1.y2 <- sd.y1 * sd.y2
    t1 <- (Y1c * Y2c) / sd.y1.y2
    t2 <- (Y1c * Y1c) / evar.y1 - (2 * rho * t1) + (Y2c * Y2c) / evar.y2
    dx <- (rho + t1 - t2 * rho / R) / R

    # to be consistent with (log)lik_cache
    if (length(lik.toosmall.idx) > 0L) {
      dx[lik.toosmall.idx] <- as.numeric(NA)
    }

    if (is.null(wt)) {
      dx.rho <- sum(dx, na.rm = TRUE)
    } else {
      dx.rho <- sum(wt * dx, na.rm = TRUE)
    }

    return(dx.rho)
  })
}

lav_bvreg_hessian_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    rho2 <- rho * rho
    R2 <- R * R
    R3 <- R * R * R

    h <- 1 / R - (2 * rho2 * t2) / R3 + 2 * rho2 * (1 - t2 / R) / R2 + 4 * rho * t1 / R2 - t2 / R2

    # to be consistent with (log)lik_cache
    if (length(lik.toosmall.idx) > 0L) {
      h[lik.toosmall.idx] <- as.numeric(NA)
    }

    if (is.null(wt)) {
      H <- sum(h, na.rm = TRUE)
    } else {
      H <- sum(wt * h, na.rm = TRUE)
    }
    dim(H) <- c(1L, 1L) # for nlminb

    return(H)
  })
}

# compute total (log)likelihood, for specific 'x' (nlminb)
lav_bvreg_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  -1 * lav_bvreg_logl_cache(cache = cache) / cache$N
}

# compute gradient, for specific 'x' (nlminb)
lav_bvreg_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_bvreg_logl_cache(cache = cache)
  }
  -1 * lav_bvreg_gradient_cache(cache = cache) / cache$N
}

# compute hessian, for specific 'x' (nlminb)
lav_bvreg_min_hessian <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    tmp <- lav_bvreg_logl_cache(cache = cache)
    tmp <- lav_bvreg_gradient_cache(cache = cache)
  }
  -1 * lav_bvreg_hessian_cache(cache = cache) / cache$N
}

# casewise scores - cache
# FIXME: should we also set 'lik.toosmall.idx' cases to NA?
lav_bvreg_cor_scores_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]
    R <- (1 - rho * rho)

    # mu.y1
    dx.mu.y1 <- (2 * Y1c / evar.y1 - 2 * rho * Y2c / (sd.y1 * sd.y2)) / (2 * R)
    if (!is.null(wt)) {
      dx.mu.y1 <- wt * dx.mu.y1
    }

    # mu.y2
    dx.mu.y2 <- -(2 * rho * Y1c / (sd.y1 * sd.y2) - 2 * Y2c / evar.y2) / (2 * R)
    if (!is.null(wt)) {
      dx.mu.y2 <- wt * dx.mu.y2
    }

    # evar.y1
    dx.var.y1 <- -(0.5 / evar.y1 - ((Y1c * Y1c) / (evar.y1 * evar.y1) -
      rho * Y1c * Y2c / (evar.y1 * sd.y1 * sd.y2)) / (2 * R))
    if (!is.null(wt)) {
      dx.var.y1 <- wt * dx.var.y1
    }

    # var.y2
    dx.var.y2 <- -(0.5 / evar.y2 + (rho * Y1c * Y2c / (evar.y2 * sd.y1 * sd.y2) -
      (Y2c * Y2c) / (evar.y2 * evar.y2)) / (2 * R))
    if (!is.null(wt)) {
      dx.var.y2 <- wt * dx.var.y2
    }

    # sl.y1
    dx.sl.y1 <- NULL
    if (nexo > 0L) {
      dx.sl.y1 <- dx.mu.y1 * eXo # weights already included in dx.mu.y1
    }

    # sl.y2
    dx.sl.y2 <- NULL
    if (nexo > 0L) {
      dx.sl.y2 <- dx.mu.y2 * eXo # weights already included in dx.mu.y2
    }

    # rho
    z <- (Y1c * Y1c) / evar.y1 - 2 * rho * Y1c * Y2c / (sd.y1 * sd.y2) + (Y2c * Y2c) / evar.y2
    dx.rho <- rho / R + (Y1c * Y2c / (sd.y1 * sd.y2 * R) - z * rho / (R * R))
    if (!is.null(wt)) {
      dx.rho <- wt * dx.rho
    }

    out <- list(
      dx.mu.y1 = dx.mu.y1, dx.var.y1 = dx.var.y1,
      dx.mu.y2 = dx.mu.y2, dx.var.y2 = dx.var.y2,
      dx.sl.y1 = dx.sl.y1, dx.sl.y2 = dx.sl.y2,
      dx.rho = dx.rho
    )
    return(out)
  })
}

# casewise scores
#
# Y1 = linear
# Y2 = linear
lav_bvreg_cor_scores <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit.y1 = NULL, fit.y2 = NULL,
                                 evar.y1 = NULL, beta.y1 = NULL,
                                 evar.y2 = NULL, beta.y2 = NULL) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvreg_fit(y = Y2, X = eXo, wt = wt)
  }

  # user specified parameters
  if (!is.null(evar.y1) || !is.null(beta.y1)) {
    fit.y1 <- lav_uvreg_update_fit(
      fit.y = fit.y1,
      evar.new = evar.y1, beta.new = beta.y1
    )
  }
  if (!is.null(evar.y2) || !is.null(beta.y2)) {
    fit.y2 <- lav_uvreg_update_fit(
      fit.y = fit.y2,
      evar.new = evar.y2, beta.new = beta.y2
    )
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(
    fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  SC <- lav_bvreg_cor_scores_cache(cache = cache)

  SC
}

# logl - no cache
lav_bvreg_logl <- function(Y1, Y2, eXo = NULL, wt = NULL,
                           rho = NULL,
                           fit.y1 = NULL, fit.y2 = NULL,
                           evar.y1 = NULL, beta.y1 = NULL,
                           evar.y2 = NULL, beta.y2 = NULL) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvreg_fit(y = Y2, X = eXo, wt = wt)
  }

  # user specified parameters
  if (!is.null(evar.y1) || !is.null(beta.y1)) {
    fit.y1 <- lav_uvreg_update_fit(
      fit.y = fit.y1,
      evar.new = evar.y1, beta.new = beta.y1
    )
  }
  if (!is.null(evar.y2) || !is.null(beta.y2)) {
    fit.y2 <- lav_uvreg_update_fit(
      fit.y = fit.y2,
      evar.new = evar.y2, beta.new = beta.y2
    )
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(
    fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  lav_bvreg_logl_cache(cache = cache)
}

# lik - no cache
lav_bvreg_lik <- function(Y1, Y2, eXo = NULL, wt = NULL,
                          rho = NULL,
                          fit.y1 = NULL, fit.y2 = NULL,
                          evar.y1 = NULL, beta.y1 = NULL,
                          evar.y2 = NULL, beta.y2 = NULL,
                          .log = FALSE) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvreg_fit(y = Y2, X = eXo, wt = wt)
  }

  # user specified parameters
  if (!is.null(evar.y1) || !is.null(beta.y1)) {
    fit.y1 <- lav_uvreg_update_fit(
      fit.y = fit.y1,
      evar.new = evar.y1, beta.new = beta.y1
    )
  }
  if (!is.null(evar.y2) || !is.null(beta.y2)) {
    fit.y2 <- lav_uvreg_update_fit(
      fit.y = fit.y2,
      evar.new = evar.y2, beta.new = beta.y2
    )
  }

  # create cache environment
  cache <- lav_bvreg_init_cache(
    fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
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
