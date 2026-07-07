# the univariate two-level Gaussian (random-intercept) model
#
#     y_ij = mu + b_j + e_ij,   b_j ~ N(0, vb),  e_ij ~ N(0, vw)
#
# YR 2026 (two-level WLS, phase 2)
#
# parameter vector: theta = c(mu, vw, vb)
# this is stage 1 of the stage-wise unrestricted two-level model for a
# continuous variable (cfr. lav_uvreg.R for the single-level case)
#
# fit object (list) contract -- see lav_uvreg_2l_fit():
#   theta, mu_idx (1), wvar_idx (2), bvar_idx (3),
#   y, cluster_idx, cs (cluster sufficient stats), logl, converged

lav_uvreg_2l_fit <- function(y = NULL,
                             cluster_idx = NULL,
                             optim_method = "nlminb1",
                             control = list(),
                             output = "list") {
  # create cache environment
  cache <- lav_uvreg_2l_init_cache(y = y, cluster_idx = cluster_idx)

  # optim
  minfns <- lav_uvbv_optim_fns(optim_method, lav_uvreg_2l_min_fns())

  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control_nlminb <- modifyList(control_nlminb, control)

  optim <- nlminb(
    start = cache$theta, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control_nlminb,
    lower = c(-Inf, .Machine$double.eps^0.75, 0),
    cache = cache
  )

  if (output == "cache") {
    return(cache)
  }

  list(
    theta = optim$par,
    mu_idx = 1L, wvar_idx = 2L, bvar_idx = 3L,
    y = cache$y,
    cluster_idx = cache$cluster_idx,
    cs = cache$cs,
    logl = -0.5 * sum(lav_2l_gauss_m2ll(
      cs = cache$cs, mu = optim$par[1],
      sigma_w = matrix(optim$par[2], 1, 1),
      sigma_b = matrix(optim$par[3], 1, 1)
    )),
    converged = (optim$convergence == 0L)
  )
}

# prepare cache environment
lav_uvreg_2l_init_cache <- function(y = NULL,
                                    cluster_idx = NULL,
                                    parent = parent.frame()) {
  y <- as.numeric(y)
  if (anyNA(y)) {
    lav_msg_stop(gettext(
      "missing values are not supported (yet) in lav_uvreg_2l_fit()."))
  }
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
  theta <- c(mu0, vw0, vb0)

  list2env(
    list(
      y = y, cluster_idx = cluster_idx, cs = cs,
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

# cluster-wise scores (of the loglik) at given (or fitted) parameter values
lav_uvreg_2l_sc <- function(fit_y = NULL, theta = NULL) {
  if (is.null(theta)) {
    theta <- fit_y$theta
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
lav_uvreg_2l_min_fns <- function() {
  lav_uvbv_min_fns(
    logl_fun = lav_uvreg_2l_loglik_cache,
    grad_fun = lav_uvreg_2l_grad_cache
  )
}
