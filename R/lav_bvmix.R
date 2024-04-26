# the weighted bivariate ordinal/linear model
# YR 08 March 2020 (replacing the old lav_polyserial.R routines)
#
# - polyserial (and biserial) correlations
# - bivariate ordinal/linear regression
# - using sampling weights wt


# polyserial correlation
#
# Y1 = linear
# Y2 = ordinal
lav_bvmix_cor_twostep_fit <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                      fit.y1 = NULL, fit.y2 = NULL,
                                      Y1.name = NULL, Y2.name = NULL,
                                      optim.method = "nlminb1", # 0.6-7
                                      optim.scale = 1.0,
                                      init.theta = NULL,
                                      control = list(),
                                      verbose = FALSE) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # create cache environment
  cache <- lav_bvmix_init_cache(fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt)

  # optim.method
  minObjective <- lav_bvmix_min_objective
  minGradient <- lav_bvmix_min_gradient
  minHessian <- lav_bvmix_min_hessian
  if (optim.method == "nlminb" || optim.method == "nlminb2") {
    # nothing to do
  } else if (optim.method == "nlminb0") {
    minGradient <- minHessian <- NULL
  } else if (optim.method == "nlminb1") {
    minHessian <- NULL
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
    scale = optim.scale, lower = -0.995, upper = +0.995,
    cache = cache
  )

  # try 2
  if (optim$convergence != 0L) {
    optim <- nlminb(
      start = start.x, objective = minObjective,
      gradient = NULL, hessian = NULL,
      control = control,
      scale = optim.scale, lower = -0.995, upper = +0.995,
      cache = cache
    )
  }

  # try 3
  if (optim$convergence != 0L) {
    optim <- nlminb(
      start = 0, objective = minObjective,
      gradient = NULL, hessian = NULL,
      control = control,
      scale = 10, lower = -0.995, upper = +0.995,
      cache = cache
    )
  }

  # try 4 -- new in 0.6-8
  if (optim$convergence != 0L) {
    optim <- optimize(
      f = minObjective, interval = c(-0.995, +0.995),
      cache = cache, tol = .Machine$double.eps
    )
    if (is.finite(optim$minimum)) {
      optim$convergence <- 0L
      optim$par <- optim$minimum
    }
  }

  # check convergence
  if (optim$convergence != 0L) {
    if (!is.null(Y1.name) && !is.null(Y2.name)) {
      lav_msg_warn(gettextf(
        "estimation polyserial correlation did not converge
        for variables %1$s and %2$s", Y1.name, Y2.name))
    } else {
      lav_msg_warn(gettext(
        "estimation polyserial correlation(s) did not always converge"))
    }
    rho <- cache$theta # starting value
  } else {
    rho <- optim$par
  }

  rho
}


# Y1 = linear
# Y2 = ordinal
lav_bvmix_init_cache <- function(fit.y1 = NULL,
                                 fit.y2 = NULL,
                                 wt = NULL,
                                 scores = FALSE,
                                 parent = parent.frame()) {
  # data
  Y1 <- fit.y1$y
  Y2 <- fit.y2$y
  eXo <- fit.y1$X

  # extract parameters

  # Y1
  y1.VAR <- fit.y1$theta[fit.y1$var.idx]
  y1.SD <- sqrt(y1.VAR)
  y1.ETA <- fit.y1$yhat
  Z <- (Y1 - y1.ETA) / y1.SD

  # Y2
  th.y2 <- fit.y2$theta[fit.y2$th.idx]

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

  # starting value -- Olsson 1982 eq 38
  if (nexo > 0L) {
    # exo
    if (is.null(wt)) {
      COR <- cor(Z, Y2, use = "pairwise.complete.obs")
      SD <- sd(Y2, na.rm = TRUE) * sqrt((N - 1) / N)
    } else {
      tmp <- na.omit(cbind(Z, Y2, wt))
      COR <- cov.wt(x = tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
      SD <- sqrt(lav_matrix_var_wt(tmp[, 2], wt = tmp[, 3]))
    }
    rho.init <- (COR * SD / sum(dnorm(th.y2)))
  } else {
    # no exo
    if (is.null(wt)) {
      COR <- cor(Y1, Y2, use = "pairwise.complete.obs")
      SD <- sd(Y2, na.rm = TRUE) * sqrt((N - 1) / N)
    } else {
      tmp <- na.omit(cbind(Y1, Y2, wt))
      COR <- cov.wt(x = tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
      SD <- sqrt(lav_matrix_var_wt(tmp[, 2], wt = tmp[, 3]))
    }
    rho.init <- (COR * SD / sum(dnorm(th.y2)))
  }

  # sanity check
  if (is.na(rho.init)) {
    rho.init <- 0.0
  } else if (abs(rho.init) > 0.9) {
    rho.init <- rho.init / 2
  }

  # parameter vector
  theta <- rho.init # only

  # different cache if scores or not
  if (scores) {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, N = N,
        y1.VAR = y1.VAR, eXo = eXo,
        y2.Y1 = fit.y2$Y1, y2.Y2 = fit.y2$Y2,
        Y1 = Y1, y1.SD = y1.SD, y1.ETA = y1.ETA, Z = Z,
        fit.y2.z1 = fit.y2$z1, fit.y2.z2 = fit.y2$z2
      ),
      parent = parent
    )
  } else {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, N = N,
        Y1 = Y1, y1.SD = y1.SD, y1.ETA = y1.ETA, Z = Z,
        fit.y2.z1 = fit.y2$z1, fit.y2.z2 = fit.y2$z2
      ),
      parent = parent
    )
  }

  out
}


# casewise likelihoods, unweighted!
lav_bvmix_lik_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]
    R <- sqrt(1 - rho * rho)

    # p(Y2|Y1)
    tauj.star <- (fit.y2.z1 - rho * Z) / R
    tauj1.star <- (fit.y2.z2 - rho * Z) / R
    py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
    # TODO, check when to use 1 - pnorm()
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps

    # p(Y1)
    py1 <- dnorm(Y1, mean = y1.ETA, sd = y1.SD)

    # lik
    lik <- py1 * py2y1

    # catch very small values
    lik.toosmall.idx <- which(lik < sqrt(.Machine$double.eps))
    lik[lik.toosmall.idx] <- as.numeric(NA)

    return(lik)
  })
}

lav_bvmix_logl_cache <- function(cache = NULL) {
  with(cache, {
    lik <- lav_bvmix_lik_cache(cache) # unweighted!

    if (!is.null(wt)) {
      logl <- sum(wt * log(lik), na.rm = TRUE)
    } else {
      logl <- sum(log(lik), na.rm = TRUE)
    }

    return(logl)
  })
}

lav_bvmix_gradient_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]

    y.Z1 <- dnorm(tauj.star)
    y.Z2 <- dnorm(tauj1.star)
    pyx.inv.R3 <- 1 / (py2y1 * R * R * R)

    # rho
    d1 <- fit.y2.z1 * rho - Z
    d2 <- fit.y2.z2 * rho - Z
    dx <- pyx.inv.R3 * (y.Z1 * d1 - y.Z2 * d2)

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

# YR 29 March 2020
# obtained by using 'Deriv' (from package Deriv) on the
# gradient function, and cleaning up
# correct, but not good enough
lav_bvmix_hessian_cache <- function(cache = NULL) {
  with(cache, {
    rho <- theta[1L]
    R2 <- R * R

    t1 <- Z - rho * tauj.star / R
    t2 <- Z - rho * tauj1.star / R

    tmp <- (y.Z1 * (d1 * ((3 * rho / R2) + tauj.star * t1 / R)
      + fit.y2.z1 + dx * R2 * t1)

    - y.Z2 * (d2 * ((3 * rho / R2) + tauj1.star * t2 / R)
        + fit.y2.z2 + dx * R2 * t2)
    )

    # to be consistent with (log)lik_cache
    if (length(lik.toosmall.idx) > 0L) {
      tmp[lik.toosmall.idx] <- as.numeric(NA)
    }

    if (is.null(wt)) {
      H <- sum(tmp * pyx.inv.R3, na.rm = TRUE)
    } else {
      H <- sum(wt * (tmp * pyx.inv.R3), na.rm = TRUE)
    }
    dim(H) <- c(1L, 1L) # for nlminb

    return(H)
  })
}

# compute total (log)likelihood, for specific 'x' (nlminb)
lav_bvmix_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  -1 * lav_bvmix_logl_cache(cache = cache) / cache$N
}

# compute gradient, for specific 'x' (nlminb)
lav_bvmix_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_bvmix_logl_cache(cache = cache)
  }
  -1 * lav_bvmix_gradient_cache(cache = cache) / cache$N
}

# compute hessian, for specific 'x' (nlminb)
lav_bvmix_min_hessian <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    tmp <- lav_bvmix_logl_cache(cache = cache)
    tmp <- lav_bvmix_gradient_cache(cache = cache)
  }
  -1 * lav_bvmix_hessian_cache(cache = cache) / cache$N
}


lav_bvmix_cor_scores_cache <- function(cache = NULL,
                                       sigma.correction = FALSE,
                                       na.zero = FALSE) {
  with(cache, {
    rho <- theta[1L]
    R <- sqrt(1 - rho * rho)

    tauj.star <- (fit.y2.z1 - rho * Z) / R
    tauj1.star <- (fit.y2.z2 - rho * Z) / R
    y.Z1 <- dnorm(tauj.star)
    y.Z2 <- dnorm(tauj1.star)

    # p(Y2|Y1)
    py2y1 <- pnorm(tauj.star) - pnorm(tauj1.star)
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps
    pyx.inv <- 1 / py2y1

    # mu.y1
    y.Z1.y.Z2 <- y.Z1 - y.Z2
    dx.mu.y1 <- 1 / y1.SD * (Z + (pyx.inv * (rho / R) * y.Z1.y.Z2))
    if (!is.null(wt)) {
      dx.mu.y1 <- wt * dx.mu.y1
    }

    # var.y1
    dx.var.y1 <- 1 / (2 * y1.VAR) * (((Z * Z) - 1) +
      (pyx.inv * rho * Z / R) * y.Z1.y.Z2)
    if (!is.null(wt)) {
      dx.var.y1 <- wt * dx.var.y1
    }

    # th.y2
    dx.th.y2 <- (y2.Y1 * y.Z1 - y2.Y2 * y.Z2) * 1 / R * pyx.inv
    if (!is.null(wt)) {
      dx.th.y2 <- wt * dx.th.y2
    }

    # sl.y1
    dx.sl.y1 <- NULL
    if (nexo > 0L) {
      dx.sl.y1 <- dx.mu.y1 * eXo
      # if(!is.null(wt)) {
      # dx.mu.y1 had already been weighted
      # }
    }

    # sl.y2
    dx.sl.y2 <- NULL
    if (nexo > 0L) {
      dx.sl.y2 <- (y.Z2 - y.Z1) * eXo * 1 / R * pyx.inv
      if (!is.null(wt)) {
        dx.sl.y2 <- wt * dx.sl.y2
      }
    }

    # rho
    TAUj <- y.Z1 * (fit.y2.z1 * rho - Z)
    TAUj1 <- y.Z2 * (fit.y2.z2 * rho - Z)
    dx.rho <- pyx.inv * 1 / (R * R * R) * (TAUj - TAUj1)
    if (!is.null(wt)) {
      dx.rho <- wt * dx.rho
    }

    # FIXME: only tested for non_exo!
    # used by pml_deriv1()
    if (sigma.correction) {
      dx.rho.orig <- dx.rho
      dx.var.y1.orig <- dx.var.y1

      # sigma
      dx.rho <- dx.rho.orig / y1.SD

      # var
      COV <- rho * y1.SD
      dx.var.y1 <- (dx.var.y1.orig -
        1 / 2 * COV / y1.VAR * 1 / y1.SD * dx.rho.orig)
    }

    out <- list(
      dx.mu.y1 = dx.mu.y1, dx.var.y1 = dx.var.y1,
      dx.th.y2 = dx.th.y2,
      dx.sl.y1 = dx.sl.y1, dx.sl.y2 = dx.sl.y2,
      dx.rho = dx.rho
    )

    return(out)
  })
}


# casewise scores
#
# Y1 = linear
# Y2 = ordinal
lav_bvmix_cor_scores <- function(Y1, Y2, eXo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit.y1 = NULL, fit.y2 = NULL,
                                 evar.y1 = NULL, beta.y1 = NULL,
                                 th.y2 = NULL, sl.y2 = NULL,
                                 sigma.correction = FALSE,
                                 na.zero = FALSE) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # update z1/z2 if needed (used in pml_deriv1() in lav_model_gradient_pml.R)
  fit.y1 <- lav_uvreg_update_fit(
    fit.y = fit.y1, evar.new = evar.y1,
    beta.new = beta.y1
  )
  fit.y2 <- lav_uvord_update_fit(
    fit.y = fit.y2, th.new = th.y2,
    sl.new = sl.y2
  )

  # create cache environment
  cache <- lav_bvmix_init_cache(
    fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  SC <- lav_bvmix_cor_scores_cache(
    cache = cache,
    sigma.correction = sigma.correction,
    na.zero = na.zero
  )

  SC
}

# logl - no cache
lav_bvmix_logl <- function(Y1, Y2, eXo = NULL, wt = NULL,
                           rho = NULL,
                           fit.y1 = NULL, fit.y2 = NULL,
                           evar.y1 = NULL, beta.y1 = NULL,
                           th.y2 = NULL, sl.y2 = NULL) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # update z1/z2 if needed (used in pml_deriv1() in lav_model_gradient_pml.R)
  fit.y1 <- lav_uvreg_update_fit(
    fit.y = fit.y1, evar.new = evar.y1,
    beta.new = beta.y1
  )
  fit.y2 <- lav_uvord_update_fit(
    fit.y = fit.y2, th.new = th.y2,
    sl.new = sl.y2
  )

  # create cache environment
  cache <- lav_bvmix_init_cache(
    fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  lav_bvmix_logl_cache(cache = cache)
}

# lik - no cache
lav_bvmix_lik <- function(Y1, Y2, eXo = NULL, wt = NULL,
                          rho = NULL,
                          fit.y1 = NULL, fit.y2 = NULL,
                          evar.y1 = NULL, beta.y1 = NULL,
                          th.y2 = NULL, sl.y2 = NULL,
                          .log = FALSE) {
  if (is.null(fit.y1)) {
    fit.y1 <- lav_uvreg_fit(y = Y1, X = eXo, wt = wt)
  }
  if (is.null(fit.y2)) {
    fit.y2 <- lav_uvord_fit(y = Y2, X = eXo, wt = wt)
  }

  # update z1/z2 if needed (used in pml_deriv1() in lav_model_gradient_pml.R)
  fit.y1 <- lav_uvreg_update_fit(
    fit.y = fit.y1, evar.new = evar.y1,
    beta.new = beta.y1
  )
  fit.y2 <- lav_uvord_update_fit(
    fit.y = fit.y2, th.new = th.y2,
    sl.new = sl.y2
  )

  # create cache environment
  cache <- lav_bvmix_init_cache(
    fit.y1 = fit.y1, fit.y2 = fit.y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  lik <- lav_bvmix_lik_cache(cache = cache) # unweighted
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
