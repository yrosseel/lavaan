# the bivariate two-level Gaussian model, univariate parameters fixed
#
#     (y_a, y_b)_ij = mu + b_j + e_ij
#     b_j ~ N(0, Sigma_b),  e_ij ~ N(0, Sigma_w)
#
#     Sigma_w = | vw_a  cw   |     Sigma_b = | vb_a  cb   |
#               | cw    vw_b |               | cb    vb_b |
#
# YR 2026 (two-level WLS, phase 2)
#
# This is stage 2 of the stage-wise unrestricted two-level model for a
# pair of continuous variables: the univariate parameters (mu, vw, vb of
# both variables) are held fixed at their stage-1 values, and only the
# two covariances theta = c(cw, cb) are estimated ('twostep', cfr.
# lav_bvreg.R for the single-level case; note that here a pair has TWO
# free parameters).
#
# lav_bvreg_2l_sc() returns the cluster-wise scores with respect to ALL
# parameters (both fixed univariate and the two covariances), as needed
# for the stage-wise sandwich (A21 blocks).

lav_bvreg_2l_cov_twostep_fit <- function(fit_y1 = NULL,
                                         fit_y2 = NULL,
                                         optim_method = "nlminb1",
                                         control = list(),
                                         y1_name = NULL, y2_name = NULL,
                                         output = "theta") {
  # create cache environment
  cache <- lav_bvreg_2l_init_cache(fit_y1 = fit_y1, fit_y2 = fit_y2)

  # degenerate between level? (one of the between variances ~ 0)
  cb_bound <- 0.999 * sqrt(cache$vb1 * cache$vb2)
  cw_bound <- 0.999 * sqrt(cache$vw1 * cache$vw2)
  cb_fixed <- (cb_bound < 1e-08)

  # optim
  minfns <- lav_uvbv_optim_fns(optim_method, lav_bvreg_2l_min_fns())

  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control_nlminb <- modifyList(control_nlminb, control)

  if (cb_fixed) {
    cache$theta[2] <- 0
    lower <- c(-cw_bound, 0)
    upper <- c(+cw_bound, 0)
  } else {
    lower <- c(-cw_bound, -cb_bound)
    upper <- c(+cw_bound, +cb_bound)
  }

  # try 1
  optim <- nlminb(
    start = cache$theta, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control_nlminb, lower = lower, upper = upper,
    cache = cache
  )

  # try 2: start from zero
  if (optim$convergence != 0L) {
    optim <- nlminb(
      start = c(0, 0), objective = minfns$objective,
      gradient = NULL, hessian = NULL,
      control = control_nlminb, lower = lower, upper = upper,
      cache = cache
    )
  }

  if (optim$convergence != 0L) {
    if (!is.null(y1_name) && !is.null(y2_name)) {
      lav_msg_warn(gettextf(
        "estimation of the two-level covariances did not converge for
        variables %1$s and %2$s", y1_name, y2_name))
    } else {
      lav_msg_warn(gettext(
        "estimation of the two-level covariances did not always converge"))
    }
  }

  if (output == "cache") {
    cache$theta <- optim$par
    return(cache)
  }

  # theta = c(cw, cb)
  optim$par
}

# prepare cache environment
# fit_y1/fit_y2: lav_uvreg_2l_fit() objects of the two variables
#
# with covariates (fit$nexo > 0) the pair likelihood is evaluated on the
# RESIDUALS e = y - z theta_m (all mean parameters fixed at their stage-1
# values): the pair model is then the no-covariate model with mu = 0
lav_bvreg_2l_init_cache <- function(fit_y1 = NULL,
                                    fit_y2 = NULL,
                                    parent = parent.frame()) {
  stopifnot(length(fit_y1$y) == length(fit_y2$y))
  cluster_idx <- fit_y1$cluster_idx

  # fixed univariate parameters
  mu1 <- fit_y1$theta[1]; vw1 <- fit_y1$theta[2]; vb1 <- fit_y1$theta[3]
  mu2 <- fit_y2$theta[1]; vw2 <- fit_y2$theta[2]; vb2 <- fit_y2$theta[3]

  # residualize (covariates)
  e1 <- fit_y1$y
  if (fit_y1$nexo > 0L) {
    theta_m1 <- fit_y1$theta[c(1L, 3L + seq_len(fit_y1$nexo))]
    e1 <- fit_y1$y - drop(fit_y1$z %*% theta_m1)
    mu1 <- 0
  }
  e2 <- fit_y2$y
  if (fit_y2$nexo > 0L) {
    theta_m2 <- fit_y2$theta[c(1L, 3L + seq_len(fit_y2$nexo))]
    e2 <- fit_y2$y - drop(fit_y2$z %*% theta_m2)
    mu2 <- 0
  }

  # pair sufficient statistics
  cs <- lav_2l_cluster_stats(
    y = cbind(e1, e2),
    cluster_idx = cluster_idx
  )

  # ANOVA-type starting values for (cw, cb)
  nobs <- cs$nobs
  nclusters <- cs$nclusters
  njs <- cs$njs
  # within cross-products: wvech columns are (11, 21, 22)
  cw0 <- cs$wtot[2] / (nobs - nclusters)
  csb <- sum(njs * (cs$ybar[, 1] - mean(e1)) *
               (cs$ybar[, 2] - mean(e2))) / (nclusters - 1)
  n0 <- (nobs - sum(njs^2) / nobs) / (nclusters - 1)
  cb0 <- (csb - cw0) / n0
  # keep the starting values admissible
  cw_max <- 0.95 * sqrt(vw1 * vw2)
  cb_max <- 0.95 * sqrt(vb1 * vb2)
  cw0 <- max(min(cw0, cw_max), -cw_max)
  cb0 <- max(min(cb0, cb_max), -cb_max)
  if (!is.finite(cb0)) {
    cb0 <- 0
  }

  theta <- c(cw0, cb0)

  list2env(
    list(
      cs = cs, n = nobs, cluster_idx = cluster_idx,
      e = cbind(e1, e2),
      z1 = fit_y1$z, z2 = fit_y2$z,
      nexo1 = fit_y1$nexo, nexo2 = fit_y2$nexo,
      nexo_w1 = fit_y1$nexo_w, nexo_b1 = fit_y1$nexo_b,
      nexo_w2 = fit_y2$nexo_w, nexo_b2 = fit_y2$nexo_b,
      mu1 = mu1, vw1 = vw1, vb1 = vb1,
      mu2 = mu2, vw2 = vw2, vb2 = vb2,
      theta = theta
    ),
    parent = parent
  )
}

lav_bvreg_2l_loglik_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    cw <- theta[1]
    cb <- theta[2]
    sigma_w <- matrix(c(vw1, cw, cw, vw2), 2, 2)
    sigma_b <- matrix(c(vb1, cb, cb, vb2), 2, 2)
    m2ll_j <- lav_2l_gauss_m2ll(
      cs = cs, mu = c(mu1, mu2),
      sigma_w = sigma_w, sigma_b = sigma_b
    )
    loglik <- -0.5 * sum(m2ll_j)
    return(loglik)
  })                     # nolint end
}

# the (dmu, dSw, dSb) patterns; order:
#   1 mu_y1, 2 wvar_y1, 3 bvar_y1, 4 mu_y2, 5 wvar_y2, 6 bvar_y2,
#   7 wcov, 8 bcov
lav_bvreg_2l_params <- function() {
  e11 <- matrix(c(1, 0, 0, 0), 2, 2)
  e22 <- matrix(c(0, 0, 0, 1), 2, 2)
  e12 <- matrix(c(0, 1, 1, 0), 2, 2)
  list(
    mu_y1 = list(dmu = c(1, 0), dsw = NULL, dsb = NULL),
    wvar_y1 = list(dmu = NULL, dsw = e11, dsb = NULL),
    bvar_y1 = list(dmu = NULL, dsw = NULL, dsb = e11),
    mu_y2 = list(dmu = c(0, 1), dsw = NULL, dsb = NULL),
    wvar_y2 = list(dmu = NULL, dsw = e22, dsb = NULL),
    bvar_y2 = list(dmu = NULL, dsw = NULL, dsb = e22),
    wcov = list(dmu = NULL, dsw = e12, dsb = NULL),
    bcov = list(dmu = NULL, dsw = NULL, dsb = e12)
  )
}

lav_bvreg_2l_grad_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    cw <- theta[1]
    cb <- theta[2]
    sigma_w <- matrix(c(vw1, cw, cw, vw2), 2, 2)
    sigma_b <- matrix(c(vb1, cb, cb, vb2), 2, 2)
    dm2ll <- lav_2l_gauss_dm2ll(
      cs = cs, mu = c(mu1, mu2),
      sigma_w = sigma_w, sigma_b = sigma_b,
      params = lav_bvreg_2l_params()[c("wcov", "bcov")]
    )
    return(-0.5 * colSums(dm2ll))
  })                     # nolint end
}

# nlminb objective/gradient (see lav_uvbv_common.R)
lav_bvreg_2l_min_fns <- function() {
  lav_uvbv_min_fns(
    logl_fun = lav_bvreg_2l_loglik_cache,
    grad_fun = lav_bvreg_2l_grad_cache
  )
}

# cluster-wise scores (of the loglik) w.r.t. ALL parameters
# (mu_y1, wvar_y1, bvar_y1, mu_y2, wvar_y2, bvar_y2, wcov, bcov
#  [, y1 slopes, y2 slopes])
lav_bvreg_2l_sc <- function(fit_y1 = NULL, fit_y2 = NULL,
                            cw = NULL, cb = NULL) {
  cache <- lav_bvreg_2l_init_cache(fit_y1 = fit_y1, fit_y2 = fit_y2)
  sigma_w <- matrix(c(cache$vw1, cw, cw, cache$vw2), 2, 2)
  sigma_b <- matrix(c(cache$vb1, cb, cb, cache$vb2), 2, 2)
  dm2ll <- lav_2l_gauss_dm2ll(
    cs = cache$cs, mu = c(cache$mu1, cache$mu2),
    sigma_w = sigma_w, sigma_b = sigma_b,
    params = lav_bvreg_2l_params()
  )
  sc <- -0.5 * dm2ll
  colnames(sc) <- names(lav_bvreg_2l_params())

  # slope columns (covariates): rowsum(z * m) with m the GLS-weighted
  # residuals
  if (cache$nexo1 > 0L || cache$nexo2 > 0L) {
    m <- lav_2l_gauss_mresid(
      e = cache$e, cluster_idx = cache$cluster_idx,
      sigma_w = sigma_w, sigma_b = sigma_b
    )
    if (cache$nexo1 > 0L) {
      sc_sl1 <- rowsum.default(
        cache$z1[, -1, drop = FALSE] * m[, 1],
        group = cache$cluster_idx, reorder = TRUE
      )
      colnames(sc_sl1) <- paste0(
        "y1_", lav_uvreg_2l_slope_names(cache$nexo_w1, cache$nexo_b1)
      )
      sc <- cbind(sc, sc_sl1)
    }
    if (cache$nexo2 > 0L) {
      sc_sl2 <- rowsum.default(
        cache$z2[, -1, drop = FALSE] * m[, 2],
        group = cache$cluster_idx, reorder = TRUE
      )
      colnames(sc_sl2) <- paste0(
        "y2_", lav_uvreg_2l_slope_names(cache$nexo_w2, cache$nexo_b2)
      )
      sc <- cbind(sc, sc_sl2)
    }
  }
  sc
}

# total loglik at given values - no cache
lav_bvreg_2l_logl <- function(fit_y1 = NULL, fit_y2 = NULL,
                              cw = NULL, cb = NULL) {
  cache <- lav_bvreg_2l_init_cache(fit_y1 = fit_y1, fit_y2 = fit_y2)
  cache$theta <- c(cw, cb)
  lav_bvreg_2l_loglik_cache(cache = cache)
}
