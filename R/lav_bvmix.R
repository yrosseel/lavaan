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
#
#  info concerning lintr package, Luc DW May 7, 2026
#  there is a known problem in codetools, transferred to lintr, that
#  variables assigned in a with() are marked as 'may not be used';
#  to this end the code in with() statements are excluded from linting!
#
lav_bvmix_cor_twostep_fit <- function(y1, y2, exo = NULL, wt = NULL,
                                      fit_y1 = NULL, fit_y2 = NULL,
                                      y1_name = NULL, y2_name = NULL,
                                      optim_method = "nlminb1", # 0.6-7
                                      optim_scale = 1.0,
                                      init_theta = NULL,
                                      control = list()) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvord_fit(y = y2, x = exo, wt = wt)
  }

  # create cache environment
  cache <- lav_bvmix_init_cache(fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt)

  # optim.method
  min_objective <- lav_bvmix_min_objective
  min_gradient <- lav_bvmix_min_gradient
  min_hessian <- lav_bvmix_min_hessian
  if (optim_method == "nlminb" || optim_method == "nlminb2") {
    # nothing to do
  } else if (optim_method == "nlminb0") {
    min_gradient <- min_hessian <- NULL
  } else if (optim_method == "nlminb1") {
    min_hessian <- NULL
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
    scale = optim_scale, lower = -0.995, upper = +0.995,
    cache = cache
  )

  # try 2
  if (optim$convergence != 0L) {
    optim <- nlminb(
      start = start_x, objective = min_objective,
      gradient = NULL, hessian = NULL,
      control = control,
      scale = optim_scale, lower = -0.995, upper = +0.995,
      cache = cache
    )
  }

  # try 3
  if (optim$convergence != 0L) {
    optim <- nlminb(
      start = 0, objective = min_objective,
      gradient = NULL, hessian = NULL,
      control = control,
      scale = 10, lower = -0.995, upper = +0.995,
      cache = cache
    )
  }

  # try 4 -- new in 0.6-8
  if (optim$convergence != 0L) {
    optim <- optimize(
      f = min_objective, interval = c(-0.995, +0.995),
      cache = cache, tol = .Machine$double.eps
    )
    if (is.finite(optim$minimum)) {
      optim$convergence <- 0L
      optim$par <- optim$minimum
    }
  }

  # check convergence
  if (optim$convergence != 0L) {
    if (!is.null(y1_name) && !is.null(y2_name)) {
      lav_msg_warn(gettextf(
        "estimation polyserial correlation did not converge
        for variables %1$s and %2$s", y1_name, y2_name))
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
lav_bvmix_init_cache <- function(fit_y1 = NULL,
                                 fit_y2 = NULL,
                                 wt = NULL,
                                 scores = FALSE,
                                 parent = parent.frame()) {
  # data
  y1 <- fit_y1$y
  y2 <- fit_y2$y
  exo <- fit_y1$x

  # extract parameters

  # Y1
  y1_var <- fit_y1$theta[fit_y1$var_idx]
  y1_sd <- sqrt(y1_var)
  y1_eta <- fit_y1$yhat
  z <- (y1 - y1_eta) / y1_sd

  # Y2
  th_y2 <- fit_y2$theta[fit_y2$th_idx]

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

  # starting value -- Olsson 1982 eq 38
  if (nexo > 0L) {
    # exo
    if (is.null(wt)) {
      cor_1 <- cor(z, y2, use = "pairwise.complete.obs")
      sd_1 <- sd(y2, na.rm = TRUE) * sqrt((n - 1) / n)
    } else {
      tmp <- na.omit(cbind(z, y2, wt))
      cor_1 <- cov.wt(x = tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
      sd_1 <- sqrt(lav_matrix_var_wt(tmp[, 2], wt = tmp[, 3]))
    }
    rho_init <- (cor_1 * sd_1 / sum(dnorm(th_y2)))
  } else {
    # no exo
    if (is.null(wt)) {
      cor_1 <- cor(y1, y2, use = "pairwise.complete.obs")
      sd_1 <- sd(y2, na.rm = TRUE) * sqrt((n - 1) / n)
    } else {
      tmp <- na.omit(cbind(y1, y2, wt))
      cor_1 <- cov.wt(x = tmp[, 1:2], wt = tmp[, 3], cor = TRUE)$cor[2, 1]
      sd_1 <- sqrt(lav_matrix_var_wt(tmp[, 2], wt = tmp[, 3]))
    }
    rho_init <- (cor_1 * sd_1 / sum(dnorm(th_y2)))
  }

  # sanity check
  if (is.na(rho_init)) {
    rho_init <- 0.0
  } else if (abs(rho_init) > 0.9) {
    rho_init <- rho_init / 2
  }

  # parameter vector
  theta <- rho_init # only

  # different cache if scores or not
  if (scores) {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, n = n,
        y1_var = y1_var, exo = exo,
        y2_y1 = fit_y2$y1, y2_y2 = fit_y2$y2,
        y1 = y1, y1_sd = y1_sd, y1_eta = y1_eta, z = z,
        fit_y2_z1 = fit_y2$z1, fit_y2_z2 = fit_y2$z2
      ),
      parent = parent
    )
  } else {
    out <- list2env(
      list(
        nexo = nexo, theta = theta, n = n,
        y1 = y1, y1_sd = y1_sd, y1_eta = y1_eta, z = z,
        fit_y2_z1 = fit_y2$z1, fit_y2_z2 = fit_y2$z2
      ),
      parent = parent
    )
  }

  out
}


# casewise likelihoods, unweighted!
lav_bvmix_lik_cache <- function(cache = NULL) {
  with(cache, {                # nolint start
    rho <- theta[1L]
    r <- sqrt(1 - rho * rho)

    # p(Y2|Y1)
    tauj_star <- (fit_y2_z1 - rho * z) / r
    tauj1_star <- (fit_y2_z2 - rho * z) / r
    py2y1 <- pnorm(tauj_star) - pnorm(tauj1_star)
    # TODO, check when to use 1 - pnorm()
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps

    # p(Y1)
    py1 <- dnorm(y1, mean = y1_eta, sd = y1_sd)

    # lik
    lik <- py1 * py2y1

    # catch very small values -- TODO: sqrt() too aggressive?
    lik_toosmall_idx <- which(lik < sqrt(.Machine$double.eps))
    lik[lik_toosmall_idx] <- as.numeric(NA)

    return(lik)
  })                                 # nolint end
}

lav_bvmix_logl_cache <- function(cache = NULL) {
  with(cache, {                      # nolint start
    lik <- lav_bvmix_lik_cache(cache) # unweighted!

    if (!is.null(wt)) {
      logl <- sum(wt * log(lik), na.rm = TRUE)
    } else {
      logl <- sum(log(lik), na.rm = TRUE)
    }

    return(logl)
  })                                # nolint end
}

lav_bvmix_gradient_cache <- function(cache = NULL) {
  with(cache, {                     # nolint start
    rho <- theta[1L]

    y_z1 <- dnorm(tauj_star)
    y_z2 <- dnorm(tauj1_star)
    pyx_inv_r3 <- 1 / (py2y1 * r * r * r)

    # rho
    d1 <- fit_y2_z1 * rho - z
    d2 <- fit_y2_z2 * rho - z
    dx <- pyx_inv_r3 * (y_z1 * d1 - y_z2 * d2)

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
  })                                  # nolint end
}

# YR 29 March 2020
# obtained by using 'Deriv' (from package Deriv) on the
# gradient function, and cleaning up
# correct, but not good enough
lav_bvmix_hessian_cache <- function(cache = NULL) {
  with(cache, {                       # nolint start
    rho <- theta[1L]
    r2 <- r * r

    t1 <- z - rho * tauj_star / r
    t2 <- z - rho * tauj1_star / r

    tmp <- (y_z1 * (d1 * ((3 * rho / r2) + tauj_star * t1 / r)
      + fit_y2_z1 + dx * r2 * t1)

    - y_z2 * (d2 * ((3 * rho / r2) + tauj1_star * t2 / r)
        + fit_y2_z2 + dx * r2 * t2)
    )

    # to be consistent with (log)lik_cache
    if (length(lik_toosmall_idx) > 0L) {
      tmp[lik_toosmall_idx] <- as.numeric(NA)
    }

    if (is.null(wt)) {
      h <- sum(tmp * pyx_inv_r3, na.rm = TRUE)
    } else {
      h <- sum(wt * (tmp * pyx_inv_r3), na.rm = TRUE)
    }
    dim(h) <- c(1L, 1L) # for nlminb

    return(h)
  })                                    # nolint end
}

# compute total (log)likelihood, for specific 'x' (nlminb)
lav_bvmix_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  -1 * lav_bvmix_logl_cache(cache = cache) / cache$n
}

# compute gradient, for specific 'x' (nlminb)
lav_bvmix_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_bvmix_logl_cache(cache = cache)
    rm(tmp)
  }
  -1 * lav_bvmix_gradient_cache(cache = cache) / cache$n
}

# compute hessian, for specific 'x' (nlminb)
lav_bvmix_min_hessian <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    tmp <- lav_bvmix_logl_cache(cache = cache)
    tmp <- lav_bvmix_gradient_cache(cache = cache)
    rm(tmp)
  }
  -1 * lav_bvmix_hessian_cache(cache = cache) / cache$n
}


lav_bvmix_cor_scores_cache <- function(cache = NULL,
                                       sigma_correction = FALSE,
                                       na_zero = FALSE) {
  with(cache, {                           # nolint start
    rho <- theta[1L]
    r <- sqrt(1 - rho * rho)

    tauj_star <- (fit_y2_z1 - rho * z) / r
    tauj1_star <- (fit_y2_z2 - rho * z) / r
    y_z1 <- dnorm(tauj_star)
    y_z2 <- dnorm(tauj1_star)

    # p(Y2|Y1)
    py2y1 <- pnorm(tauj_star) - pnorm(tauj1_star)
    py2y1[py2y1 < .Machine$double.eps] <- .Machine$double.eps
    pyx_inv <- 1 / py2y1

    # mu.y1
    y_z1_y_z2 <- y_z1 - y_z2
    dx_mu_y1 <- 1 / y1_sd * (z + (pyx_inv * (rho / r) * y_z1_y_z2))
    if (!is.null(wt)) {
      dx_mu_y1 <- wt * dx_mu_y1
    }

    # var.y1
    dx_var_y1 <- 1 / (2 * y1_var) * (((z * z) - 1) +
      (pyx_inv * rho * z / r) * y_z1_y_z2)
    if (!is.null(wt)) {
      dx_var_y1 <- wt * dx_var_y1
    }

    # th.y2
    dx_th_y2 <- (y2_y1 * y_z1 - y2_y2 * y_z2) * 1 / r * pyx_inv
    if (!is.null(wt)) {
      dx_th_y2 <- wt * dx_th_y2
    }

    # sl.y1
    dx_sl_y1 <- NULL
    if (nexo > 0L) {
      dx_sl_y1 <- dx_mu_y1 * exo
      # if(!is.null(wt)) {
      # dx.mu.y1 had already been weighted
      # }
    }

    # sl.y2
    dx_sl_y2 <- NULL
    if (nexo > 0L) {
      dx_sl_y2 <- (y_z2 - y_z1) * exo * 1 / r * pyx_inv
      if (!is.null(wt)) {
        dx_sl_y2 <- wt * dx_sl_y2
      }
    }

    # rho
    tauj <- y_z1 * (fit_y2_z1 * rho - z)
    tauj1 <- y_z2 * (fit_y2_z2 * rho - z)
    dx_rho <- pyx_inv * 1 / (r * r * r) * (tauj - tauj1)
    if (!is.null(wt)) {
      dx_rho <- wt * dx_rho
    }

    # FIXME: only tested for non_exo!
    # used by lav_pml_dploglik_dimplied()
    if (sigma_correction) {
      dx_rho_orig <- dx_rho
      dx_var_y1_orig <- dx_var_y1

      # sigma
      dx_rho <- dx_rho_orig / y1_sd

      # var
      cov_1 <- rho * y1_sd
      dx_var_y1 <- (dx_var_y1_orig -
        1 / 2 * cov_1 / y1_var * 1 / y1_sd * dx_rho_orig)
    }

    out <- list(
      dx_mu_y1 = dx_mu_y1, dx_var_y1 = dx_var_y1,
      dx_th_y2 = dx_th_y2,
      dx_sl_y1 = dx_sl_y1, dx_sl_y2 = dx_sl_y2,
      dx_rho = dx_rho
    )

    return(out)
  })                                        # nolint end
}


# casewise scores
#
# Y1 = linear
# Y2 = ordinal
lav_bvmix_cor_scores <- function(y1, y2, exo = NULL, wt = NULL,
                                 rho = NULL,
                                 fit_y1 = NULL, fit_y2 = NULL,
                                 evar_y1 = NULL, beta_y1 = NULL,
                                 th_y2 = NULL, sl_y2 = NULL,
                                 sigma_correction = FALSE,
                                 na_zero = FALSE) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvord_fit(y = y2, x = exo, wt = wt)
  }

  # update z1/z2 if needed (used in lav_pml_dploglik_dimplied()
  #             in lav_model_gradient_pml.R)
  fit_y1 <- lav_uvreg_update_fit(
    fit_y = fit_y1, evar_new = evar_y1,
    beta_new = beta_y1
  )
  fit_y2 <- lav_uvord_update_fit(
    fit_y = fit_y2, th_new = th_y2,
    sl_new = sl_y2
  )

  # create cache environment
  cache <- lav_bvmix_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  sc <- lav_bvmix_cor_scores_cache(
    cache = cache,
    sigma_correction = sigma_correction,
    na_zero = na_zero
  )

  sc
}

# logl - no cache
lav_bvmix_logl <- function(y1, y2, exo = NULL, wt = NULL,
                           rho = NULL,
                           fit_y1 = NULL, fit_y2 = NULL,
                           evar_y1 = NULL, beta_y1 = NULL,
                           th_y2 = NULL, sl_y2 = NULL) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvord_fit(y = y2, x = exo, wt = wt)
  }

  # update z1/z2 if needed (used in lav_pml_dploglik_dimplied() in
  #                                     lav_model_gradient_pml.R)
  fit_y1 <- lav_uvreg_update_fit(
    fit_y = fit_y1, evar_new = evar_y1,
    beta_new = beta_y1
  )
  fit_y2 <- lav_uvord_update_fit(
    fit_y = fit_y2, th_new = th_y2,
    sl_new = sl_y2
  )

  # create cache environment
  cache <- lav_bvmix_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  lav_bvmix_logl_cache(cache = cache)
}

# lik - no cache
lav_bvmix_lik <- function(y1, y2, exo = NULL, wt = NULL,
                          rho = NULL,
                          fit_y1 = NULL, fit_y2 = NULL,
                          evar_y1 = NULL, beta_y1 = NULL,
                          th_y2 = NULL, sl_y2 = NULL,
                          take_log = FALSE) {
  if (is.null(fit_y1)) {
    fit_y1 <- lav_uvreg_fit(y = y1, x = exo, wt = wt)
  }
  if (is.null(fit_y2)) {
    fit_y2 <- lav_uvord_fit(y = y2, x = exo, wt = wt)
  }

  # update z1/z2 if needed (used in lav_pml_dploglik_dimplied() in
  #                                      lav_model_gradient_pml.R)
  fit_y1 <- lav_uvreg_update_fit(
    fit_y = fit_y1, evar_new = evar_y1,
    beta_new = beta_y1
  )
  fit_y2 <- lav_uvord_update_fit(
    fit_y = fit_y2, th_new = th_y2,
    sl_new = sl_y2
  )

  # create cache environment
  cache <- lav_bvmix_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, wt = wt,
    scores = TRUE
  )
  cache$theta <- rho

  lik <- lav_bvmix_lik_cache(cache = cache) # unweighted
  if (take_log) {
    lik <- log(lik)
  }

  if (!is.null(wt)) {
    if (take_log) {
      lik <- wt * lik
    } else {
      tmp <- wt * log(lik)
      lik <- exp(tmp)
    }
  }

  lik
}
