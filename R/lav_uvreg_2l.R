# the univariate two-level Gaussian (random-intercept) model
#
#     y_ij = mu + xw_ij' bw + xb_j' bb + b_j + e_ij
#     b_j ~ N(0, vb),  e_ij ~ N(0, vw)
#
# YR 2026 (two-level WLS, phase 2; covariates added in phase 7)
#
# parameter vector: theta = c(mu, vw, vb[, bw, bb])
# the slope parameters (bw = within-level covariates, bb = between-level
# covariates) are APPENDED, so that the positions of (mu, vw, vb) do not
# depend on the presence of covariates
#
# this is stage 1 of the stage-wise unrestricted two-level model for a
# continuous variable (cfr. lav_uvreg.R for the single-level case)
#
# fit object (list) contract -- see lav_uvreg_2l_fit():
#   theta, mu_idx (1), wvar_idx (2), bvar_idx (3), wslope_idx, bslope_idx,
#   nexo_w, nexo_b, nexo, between, y, x_w, x_b, z, resid,
#   cluster_idx, cs (no covariates only), logl, converged

lav_uvreg_2l_fit <- function(y = NULL,
                             x_w = NULL, # within-level covariates
                             x_b = NULL, # between-level covariates
                             cluster_idx = NULL,
                             between = TRUE, # vb free? (FALSE for
                                             # within-only variables)
                             optim_method = "nlminb1",
                             control = list(),
                             output = "list") {
  # create cache environment
  cache <- lav_uvreg_2l_init_cache(y = y, x_w = x_w, x_b = x_b,
                                   cluster_idx = cluster_idx,
                                   between = between)
  nexo <- cache$nexo

  # optim
  minfns <- lav_uvbv_optim_fns(
    optim_method, lav_uvreg_2l_min_fns(covariates = (nexo > 0L))
  )

  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control_nlminb <- modifyList(control_nlminb, control)

  lower <- c(-Inf, .Machine$double.eps^0.75, 0, rep(-Inf, nexo))
  upper <- c(+Inf, +Inf, +Inf, rep(+Inf, nexo))
  if (!between) {
    # pin vb at zero
    lower[3] <- upper[3] <- 0
    cache$theta[3] <- 0
  }

  optim <- nlminb(
    start = cache$theta, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control_nlminb,
    lower = lower, upper = upper,
    cache = cache
  )

  if (output == "cache") {
    return(cache)
  }

  # loglik (and residuals) at the optimum
  if (nexo > 0L) {
    cache$theta <- optim$par
    logl <- lav_uvreg_2l_x_loglik_cache(cache = cache)
    resid <- cache$e
  } else {
    logl <- -0.5 * sum(lav_2l_gauss_m2ll(
      cs = cache$cs, mu = optim$par[1],
      sigma_w = matrix(optim$par[2], 1, 1),
      sigma_b = matrix(optim$par[3], 1, 1)
    ))
    resid <- NULL
  }

  list(
    theta = optim$par,
    mu_idx = 1L, wvar_idx = 2L, bvar_idx = 3L,
    wslope_idx = if (cache$nexo_w > 0L) 3L + seq_len(cache$nexo_w)
                 else integer(0L),
    bslope_idx = if (cache$nexo_b > 0L) 3L + cache$nexo_w +
                   seq_len(cache$nexo_b)
                 else integer(0L),
    nexo_w = cache$nexo_w, nexo_b = cache$nexo_b, nexo = nexo,
    between = between,
    y = cache$y,
    x_w = cache$x_w, x_b = cache$x_b, z = cache$z,
    resid = resid,
    cluster_idx = cache$cluster_idx,
    cs = cache$cs,
    logl = logl,
    converged = (optim$convergence == 0L)
  )
}

# prepare cache environment
lav_uvreg_2l_init_cache <- function(y = NULL,
                                    x_w = NULL,
                                    x_b = NULL,
                                    cluster_idx = NULL,
                                    between = TRUE,
                                    parent = parent.frame()) {
  y <- as.numeric(y)
  if (anyNA(y)) {
    lav_msg_stop(gettext(
      "missing values are not supported (yet) in lav_uvreg_2l_fit()."))
  }
  nexo_w <- if (is.null(x_w)) 0L else NCOL(x_w)
  nexo_b <- if (is.null(x_b)) 0L else NCOL(x_b)
  nexo <- nexo_w + nexo_b

  if (nexo == 0L) {
    cs <- lav_2l_cluster_stats(y = y, cluster_idx = cluster_idx)

    # ANOVA-type starting values
    nobs <- cs$nobs
    nclusters <- cs$nclusters
    njs <- cs$njs
    mu0 <- mean(y)
    vw0 <- sum(cs$wtot) / (nobs - nclusters)
    msb <- sum(njs * (cs$ybar[, 1] - mu0)^2) / (nclusters - 1)
    n0 <- (nobs - sum(njs^2) / nobs) / (nclusters - 1)
    vb0 <- max((msb - vw0) / n0, 0.01 * vw0)
    if (!between) {
      vb0 <- 0
      vw0 <- sum((y - mu0)^2) / nobs
    }
    theta <- c(mu0, vw0, vb0)

    return(list2env(
      list(
        y = y, x_w = NULL, x_b = NULL, z = NULL,
        cluster_idx = cluster_idx, cs = cs,
        nexo_w = 0L, nexo_b = 0L, nexo = 0L,
        n = nobs, theta = theta
      ),
      parent = parent
    ))
  }

  # with covariates
  if (!between && nexo_b > 0L) {
    lav_msg_stop(gettext(
      "a within-only variable cannot be regressed on between-level
      covariates."))
  }
  z <- cbind(1, x_w, x_b)
  storage.mode(z) <- "double"
  if (anyNA(z)) {
    lav_msg_stop(gettext(
      "missing values are not supported (yet) in lav_uvreg_2l_fit()."))
  }
  nobs <- length(y)
  njs <- tabulate(cluster_idx)
  nclusters <- length(njs)

  # OLS starting values for the mean parameters; ANOVA-type starting
  # values for (vw, vb) from the OLS residuals
  theta_m0 <- qr.coef(qr(z), y)
  theta_m0[is.na(theta_m0)] <- 0
  r <- y - drop(z %*% theta_m0)
  rbar <- drop(rowsum.default(r, group = cluster_idx, reorder = TRUE)) / njs
  sr2 <- drop(rowsum.default(r * r, group = cluster_idx, reorder = TRUE))
  vw0 <- sum(sr2 - njs * rbar^2) / (nobs - nclusters)
  msb <- sum(njs * rbar^2) / (nclusters - 1)
  n0 <- (nobs - sum(njs^2) / nobs) / (nclusters - 1)
  vb0 <- max((msb - vw0) / n0, 0.01 * vw0)
  if (!between) {
    vb0 <- 0
    vw0 <- sum(r * r) / nobs
  }
  theta <- c(theta_m0[1], vw0, vb0, theta_m0[-1])

  list2env(
    list(
      y = y, x_w = x_w, x_b = x_b, z = z,
      cluster_idx = cluster_idx, njs = njs, nclusters = nclusters,
      nexo_w = nexo_w, nexo_b = nexo_b, nexo = nexo,
      n = nobs, theta = theta
    ),
    parent = parent
  )
}

lav_uvreg_2l_loglik_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    mu <- theta[1]
    vw <- theta[2]
    vb <- theta[3]
    m2ll_j <- lav_2l_gauss_m2ll(
      cs = cs, mu = mu,
      sigma_w = matrix(vw, 1, 1), sigma_b = matrix(vb, 1, 1)
    )
    loglik <- -0.5 * sum(m2ll_j)
    return(loglik)
  })                     # nolint end
}

# the (dmu, dSw, dSb) patterns of theta = c(mu, vw, vb)
lav_uvreg_2l_params <- function() {
  list(
    list(dmu = 1, dsw = NULL, dsb = NULL),
    list(dmu = NULL, dsw = matrix(1, 1, 1), dsb = NULL),
    list(dmu = NULL, dsw = NULL, dsb = matrix(1, 1, 1))
  )
}

lav_uvreg_2l_grad_cache <- function(cache = NULL) {
  with(cache, {           # nolint start
    mu <- theta[1]
    vw <- theta[2]
    vb <- theta[3]
    dm2ll <- lav_2l_gauss_dm2ll(
      cs = cs, mu = mu,
      sigma_w = matrix(vw, 1, 1), sigma_b = matrix(vb, 1, 1),
      params = lav_uvreg_2l_params()
    )
    return(-0.5 * colSums(dm2ll))
  })                      # nolint end
}

# --- with covariates ---
#
# per-cluster -2 loglik (univariate, mean structure mu_ij = z_ij' theta_m):
# with e_ij = y_ij - mu_ij, ebar_j the cluster residual means, sww_j the
# within-cluster residual scatter and Vj = vw + nj vb,
#
#   -2 logl_j = nj log(2 pi) + (nj - 1) log(vw) + log(Vj)
#               + sww_j / vw + nj ebar_j^2 / Vj
#
# scores: the GLS-weighted residuals m_ij = [Sigma_j^{-1} e_j]_ij have the
# closed form m_ij = (e_ij - h_j ebar_j) / vw with h_j = nj vb / Vj, and
# the score of the loglik w.r.t. any mean parameter with observation-level
# design column z is sum_ij z_ij m_ij

lav_uvreg_2l_x_loglik_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    vw <- theta[2]
    vb <- theta[3]
    if (vw <= 0 || vb < 0) {
      return(-Inf)
    }
    theta_m <- theta[c(1L, 3L + seq_len(nexo))]
    e <- y - drop(z %*% theta_m)
    ebar <- drop(rowsum.default(e, group = cluster_idx,
                                reorder = TRUE)) / njs
    se2 <- drop(rowsum.default(e * e, group = cluster_idx,
                               reorder = TRUE))
    sww <- se2 - njs * ebar^2
    vj <- vw + njs * vb
    m2ll_j <- njs * log(2 * pi) + (njs - 1) * log(vw) + log(vj) +
      sww / vw + njs * ebar^2 / vj
    loglik <- -0.5 * sum(m2ll_j)
    return(loglik)
  })                     # nolint end
}

# assumes the loglik has been evaluated at cache$theta (e, ebar, sww, vj)
lav_uvreg_2l_x_grad_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    vw <- theta[2]
    vb <- theta[3]
    h <- njs * vb / vj
    m <- (e - (h * ebar)[cluster_idx]) / vw
    g_mean <- drop(crossprod(z, m))
    g_vw <- -0.5 * sum((njs - 1) / vw + 1 / vj - sww / vw^2 -
                         njs * ebar^2 / vj^2)
    g_vb <- -0.5 * sum(njs / vj - njs^2 * ebar^2 / vj^2)
    return(c(g_mean[1], g_vw, g_vb, g_mean[-1]))
  })                     # nolint end
}

# cluster-wise scores, covariate case; columns (mu, vw, vb, slopes)
lav_uvreg_2l_x_sc <- function(y = NULL, z = NULL, cluster_idx = NULL,
                              theta = NULL, nexo = NULL) {
  njs <- tabulate(cluster_idx)
  vw <- theta[2]
  vb <- theta[3]
  theta_m <- theta[c(1L, 3L + seq_len(nexo))]
  e <- y - drop(z %*% theta_m)
  ebar <- drop(rowsum.default(e, group = cluster_idx, reorder = TRUE)) / njs
  se2 <- drop(rowsum.default(e * e, group = cluster_idx, reorder = TRUE))
  sww <- se2 - njs * ebar^2
  vj <- vw + njs * vb
  h <- njs * vb / vj
  m <- (e - (h * ebar)[cluster_idx]) / vw
  sc_mean <- rowsum.default(z * m, group = cluster_idx, reorder = TRUE)
  sc_vw <- -0.5 * ((njs - 1) / vw + 1 / vj - sww / vw^2 -
                     njs * ebar^2 / vj^2)
  sc_vb <- -0.5 * (njs / vj - njs^2 * ebar^2 / vj^2)
  cbind(sc_mean[, 1], sc_vw, sc_vb, sc_mean[, -1, drop = FALSE])
}

# slope column names (within covariates first, then between covariates)
lav_uvreg_2l_slope_names <- function(nexo_w = 0L, nexo_b = 0L) {
  c(
    if (nexo_w > 0L) paste0("wsl", seq_len(nexo_w)) else character(0L),
    if (nexo_b > 0L) paste0("bsl", seq_len(nexo_b)) else character(0L)
  )
}

# cluster-wise scores (of the loglik) at given (or fitted) parameter values
lav_uvreg_2l_sc <- function(fit_y = NULL, theta = NULL) {
  if (is.null(theta)) {
    theta <- fit_y$theta
  }
  if (fit_y$nexo > 0L) {
    sc <- lav_uvreg_2l_x_sc(
      y = fit_y$y, z = fit_y$z, cluster_idx = fit_y$cluster_idx,
      theta = theta, nexo = fit_y$nexo
    )
    colnames(sc) <- c(
      "mu", "wvar", "bvar",
      lav_uvreg_2l_slope_names(fit_y$nexo_w, fit_y$nexo_b)
    )
    return(sc)
  }
  dm2ll <- lav_2l_gauss_dm2ll(
    cs = fit_y$cs, mu = theta[1],
    sigma_w = matrix(theta[2], 1, 1), sigma_b = matrix(theta[3], 1, 1),
    params = lav_uvreg_2l_params()
  )
  sc <- -0.5 * dm2ll
  colnames(sc) <- c("mu", "wvar", "bvar")
  sc
}

# nlminb objective/gradient (see lav_uvbv_common.R)
lav_uvreg_2l_min_fns <- function(covariates = FALSE) {
  if (covariates) {
    lav_uvbv_min_fns(
      logl_fun = lav_uvreg_2l_x_loglik_cache,
      grad_fun = lav_uvreg_2l_x_grad_cache
    )
  } else {
    lav_uvbv_min_fns(
      logl_fun = lav_uvreg_2l_loglik_cache,
      grad_fun = lav_uvreg_2l_grad_cache
    )
  }
}
