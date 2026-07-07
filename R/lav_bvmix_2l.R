# the bivariate two-level linear/ordinal model, univariate parameters fixed
#
#     y_a,ij  = mu_a + b_a,j + e_a,ij           (continuous)
#     y*_c,ij = w_c,ij + b_c,j                  (ordinal, probit)
#
#     Cov(e_a, w_c) = cw (within),  Cov(b_a, b_c) = cb (between)
#     Var(e_a) = vw_a, Var(w_c) = 1, Var(b_a) = vb_a, Var(b_c) = vb_c
#
# YR 2026 (two-level WLS, phase 3)
#
# This is stage 2 of the stage-wise unrestricted two-level model for a
# mixed continuous/ordinal pair: the univariate parameters (mu_a, vw_a,
# vb_a from lav_uvreg_2l_fit; tau_c, vb_c from lav_uvord_2l_fit) are held
# fixed, and theta = c(cw, cb) is estimated by ML with ONE-dimensional
# Gauss-Hermite quadrature:
#
# the probit means depend on the cluster effects only through the scalar
#
#     gamma_j = b_c,j - (cw/vw_a) b_a,j
#
# because, conditional on (b_a, b_c) and y_a,
#
#     y*_c,ij | . ~ N( gamma_j + (cw/vw_a)(y_a,ij - mu_a),  1 - cw^2/vw_a )
#
# with independence over i. Hence
#
#     L_j = f(y_a,j) * int prod_i P(y_c,ij | gamma) phi(gamma; m_j, v_j)
#
# where f(y_a,j) is the (fixed-parameter) univariate two-level Gaussian
# density of y_a in cluster j, and gamma | y_a,j is normal with
#
#     m_j = c_g nj r_aj / v_aj,  v_j = Var(gamma) - c_g^2 nj / v_aj
#     c_g = cb - (cw/vw_a) vb_a,  Var(gamma) = vb_c + (cw/vw_a)^2 vb_a
#                                              - 2 (cw/vw_a) cb
#     v_aj = vw_a + nj vb_a,  r_aj = ybar_a,j - mu_a
#
# The gradient is posterior-weighted, with both a direct term (cw enters
# the conditional probit) and a chain term through the per-cluster node
# positions gamma_jq = m_j + sqrt(v_j) x_q.

lav_bvmix_2l_cov_twostep_fit <- function(fit_y1 = NULL, # continuous
                                         fit_y2 = NULL, # ordinal
                                         ngh = 21L,
                                         optim_method = "nlminb1",
                                         control = list(),
                                         y1_name = NULL, y2_name = NULL,
                                         output = "theta") {
  # create cache environment
  cache <- lav_bvmix_2l_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, ngh = ngh
  )

  # bounds
  cw_bound <- 0.999 * sqrt(cache$vw_a)
  cb_bound <- 0.999 * sqrt(cache$vb_a * cache$vb_c)
  cb_fixed <- (cb_bound < 1e-08)
  if (cb_fixed) {
    cache$theta[2] <- 0
    lower <- c(-cw_bound, 0)
    upper <- c(+cw_bound, 0)
  } else {
    lower <- c(-cw_bound, -cb_bound)
    upper <- c(+cw_bound, +cb_bound)
  }

  # optim
  minfns <- lav_uvbv_optim_fns(optim_method, lav_bvmix_2l_min_fns())

  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control_nlminb <- modifyList(control_nlminb, control)

  # try 1
  optim <- nlminb(
    start = cache$theta, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control_nlminb, lower = lower, upper = upper,
    cache = cache
  )

  # try 2: start from independence
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
        "estimation of the two-level polyserial parameters did not converge
        for variables %1$s and %2$s", y1_name, y2_name))
    } else {
      lav_msg_warn(gettext(
        "estimation of the two-level polyserial parameters did not always
        converge"))
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
# fit_y1: lav_uvreg_2l_fit() object (continuous)
# fit_y2: lav_uvord_2l_fit() object (ordinal)
lav_bvmix_2l_init_cache <- function(fit_y1 = NULL,
                                    fit_y2 = NULL,
                                    ngh = 21L,
                                    parent = parent.frame()) {
  stopifnot(length(fit_y1$y) == length(fit_y2$y))
  y_a <- fit_y1$y
  y_c <- fit_y2$y
  cluster_idx <- fit_y1$cluster_idx
  nobs <- length(y_a)

  # fixed univariate parameters
  mu_a <- fit_y1$theta[fit_y1$mu_idx]
  vw_a <- fit_y1$theta[fit_y1$wvar_idx]
  vb_a <- fit_y1$theta[fit_y1$bvar_idx]
  tau_c <- fit_y2$theta[fit_y2$th_idx]
  vb_c <- fit_y2$theta[fit_y2$bvar_idx]
  nth <- length(tau_c)

  # cluster structure (reuse the stage-1 sufficient statistics of y_a)
  cs_a <- fit_y1$cs
  njs <- cs_a$njs
  nclusters <- cs_a$nclusters
  cl_sorted <- sort(unique(cluster_idx))
  cl_pos <- match(cluster_idx, cl_sorted)

  # fixed per-cluster quantities
  v_aj <- vw_a + njs * vb_a
  r_aj <- cs_a$ybar[, 1] - mu_a

  # fixed constant: log f(y_a,j) (does not depend on cw/cb)
  logf_a <- -0.5 * lav_2l_gauss_m2ll(
    cs = cs_a, mu = mu_a,
    sigma_w = matrix(vw_a, 1, 1), sigma_b = matrix(vb_a, 1, 1)
  )

  # per-item pieces
  dev_a <- y_a - mu_a
  o12 <- lav_uvord_o12(y = y_c, nth = nth)
  y1 <- matrix(seq_len(nth), nobs, nth, byrow = TRUE) == y_c
  y2 <- matrix(seq_len(nth), nobs, nth, byrow = TRUE) == (y_c - 1L)

  # standardized GH nodes
  gh <- lav_2l_gh_xw(ngh = ngh)

  # starting values: crude split of the total correlation
  r_tot <- suppressWarnings(cor(y_a, y_c))
  if (!is.finite(r_tot)) {
    r_tot <- 0
  }
  cw0 <- max(min(r_tot * sqrt(vw_a), 0.9 * sqrt(vw_a)),
             -0.9 * sqrt(vw_a))
  cb_max <- 0.9 * sqrt(vb_a * vb_c)
  cb0 <- max(min(r_tot * sqrt(vb_a * vb_c), cb_max), -cb_max)
  theta <- c(cw0, cb0)

  list2env(
    list(
      y_c = y_c, cluster_idx = cluster_idx, cl_pos = cl_pos,
      nobs = nobs, n = nobs, nclusters = nclusters, njs = njs,
      mu_a = mu_a, vw_a = vw_a, vb_a = vb_a,
      tau_c = tau_c, vb_c = vb_c, nth = nth,
      v_aj = v_aj, r_aj = r_aj, logf_a = logf_a,
      dev_a = dev_a, o1 = o12$o1, o2 = o12$o2, y1 = y1, y2 = y2,
      gh = gh, ngh = ngh,
      theta = theta
    ),
    parent = parent
  )
}

# the (cw, cb)-dependent per-cluster quantities
lav_bvmix_2l_gamma_pars <- function(cache = NULL) {
  with(cache, {          # nolint start
    cw <- theta[1]
    cb <- theta[2]
    a_sl <- cw / vw_a
    sd_c2 <- 1 - cw * cw / vw_a
    if (sd_c2 <= 0) {
      return(NULL)
    }
    sd_c <- sqrt(sd_c2)
    c_g <- cb - a_sl * vb_a
    var_g <- vb_c + a_sl * a_sl * vb_a - 2 * a_sl * cb
    m_j <- c_g * njs * r_aj / v_aj
    v_gj <- var_g - c_g * c_g * njs / v_aj
    if (any(v_gj <= 0)) {
      return(NULL)
    }
    out <- list(
      a_sl = a_sl, sd_c = sd_c, c_g = c_g, var_g = var_g,
      m_j = m_j, v_gj = v_gj, s_gj = sqrt(v_gj)
    )
    return(out)
  })                     # nolint end
}

lav_bvmix_2l_loglik_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    gp <- lav_bvmix_2l_gamma_pars(cache)
    if (is.null(gp)) {
      loglik <- -Inf
      return(loglik)
    }

    th_s <- tau_c / gp$sd_c
    base_i <- gp$a_sl * dev_a / gp$sd_c # per item, fixed over nodes

    ll <- matrix(0, nclusters, ngh)
    for (q in seq_len(ngh)) {
      gamma_j <- gp$m_j + gp$s_gj * gh$x[q]
      offset_i <- gamma_j[cl_pos] / gp$sd_c + base_i
      z12 <- lav_uvord_z12(y = y_c, th = th_s, o1 = o1, o2 = o2,
                           offset = offset_i)
      pi_i <- pnorm(z12$z1) - pnorm(z12$z2)
      large_idx <- which(z12$z2 > 1)
      if (length(large_idx) > 0L) {
        pi_i[large_idx] <- (pnorm(z12$z2[large_idx], lower.tail = FALSE) -
          pnorm(z12$z1[large_idx], lower.tail = FALSE))
      }
      pi_i[pi_i < .Machine$double.eps] <- .Machine$double.eps
      ll[, q] <- rowsum.default(log(pi_i), group = cluster_idx,
                                reorder = TRUE)
    }

    lwll <- sweep(ll, 2L, gh$logw, "+")
    mmax <- apply(lwll, 1L, max)
    logl_j <- mmax + log(rowSums(exp(lwll - mmax)))
    pw <- exp(lwll - logl_j)

    loglik <- sum(logl_j + logf_a)
    return(loglik)
  })                     # nolint end
}

lav_bvmix_2l_grad_cache <- function(cache = NULL) {
  sc <- lav_bvmix_2l_sc_cache(cache = cache)
  colSums(sc[, c("cw", "cb"), drop = FALSE])
}

# cluster-wise scores of the loglik w.r.t. (tau_c, cw, cb)
# (the remaining fixed-parameter columns are added in the sandwich stage)
# assumes the loglik has been evaluated at cache$theta (pw available)
lav_bvmix_2l_sc_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    cw <- theta[1]
    cb <- theta[2]
    gp <- lav_bvmix_2l_gamma_pars(cache)

    a_sl <- gp$a_sl
    sd_c <- gp$sd_c
    c_g <- gp$c_g

    # derivatives of the ingredients
    da_dcw <- 1 / vw_a
    dsd_dcw <- -cw / (vw_a * sd_c)
    dcg_dcw <- -vb_a / vw_a
    dcg_dcb <- 1
    dvarg_dcw <- (2 * a_sl * vb_a - 2 * cb) / vw_a
    dvarg_dcb <- -2 * a_sl

    # per-cluster chain pieces: gamma_jq = m_j + s_gj x_q
    dm_dcw <- dcg_dcw * njs * r_aj / v_aj
    dm_dcb <- dcg_dcb * njs * r_aj / v_aj
    dvg_dcw <- dvarg_dcw - 2 * c_g * dcg_dcw * njs / v_aj
    dvg_dcb <- dvarg_dcb - 2 * c_g * dcg_dcb * njs / v_aj
    dsg_dcw <- dvg_dcw / (2 * gp$s_gj)
    dsg_dcb <- dvg_dcb / (2 * gp$s_gj)

    th_s <- tau_c / gp$sd_c
    base_i <- a_sl * dev_a / sd_c

    sc <- matrix(0, nclusters, nth + 2L)
    for (q in seq_len(ngh)) {
      gamma_j <- gp$m_j + gp$s_gj * gh$x[q]
      offset_i <- gamma_j[cl_pos] / sd_c + base_i
      z12 <- lav_uvord_z12(y = y_c, th = th_s, o1 = o1, o2 = o2,
                           offset = offset_i)
      z1 <- z12$z1
      z2 <- z12$z2
      pi_i <- pnorm(z1) - pnorm(z2)
      large_idx <- which(z2 > 1)
      if (length(large_idx) > 0L) {
        pi_i[large_idx] <- (pnorm(z2[large_idx], lower.tail = FALSE) -
          pnorm(z1[large_idx], lower.tail = FALSE))
      }
      pi_i[pi_i < .Machine$double.eps] <- .Machine$double.eps
      inv_pi <- 1 / pi_i
      p1 <- dnorm(z1)
      p2 <- dnorm(z2)

      # thresholds (direct): dz/dtau_k = y_k / sd_c
      g_tau <- rowsum.default(
        (y1 * p1 - y2 * p2) * (inv_pi / sd_c),
        group = cluster_idx, reorder = TRUE
      )

      # gamma: dz/dgamma = -1/sd_c
      g_gam <- rowsum.default((p2 - p1) * (inv_pi / sd_c),
                              group = cluster_idx, reorder = TRUE)

      # direct cw term:
      # dz1/dcw = -(dev_i da_dcw)/sd_c - z1 * dsd_dcw/sd_c
      # (z1 includes the +-100 offsets of the extreme categories, but
      #  there p1 = dnorm(z1) is exactly zero, so the product vanishes)
      dz1_dcw <- -(dev_a * da_dcw) / sd_c - z1 * dsd_dcw / sd_c
      dz2_dcw <- -(dev_a * da_dcw) / sd_c - z2 * dsd_dcw / sd_c
      g_cw_dir <- rowsum.default((p1 * dz1_dcw - p2 * dz2_dcw) * inv_pi,
                                 group = cluster_idx, reorder = TRUE)

      # chain terms through the node positions
      dgam_dcw <- dm_dcw + gh$x[q] * dsg_dcw
      dgam_dcb <- dm_dcb + gh$x[q] * dsg_dcb

      pwq <- pw[, q]
      sc[, seq_len(nth)] <- sc[, seq_len(nth)] + pwq * g_tau
      sc[, nth + 1L] <- sc[, nth + 1L] +
        pwq * (g_cw_dir + g_gam * dgam_dcw)
      sc[, nth + 2L] <- sc[, nth + 2L] + pwq * (g_gam * dgam_dcb)
    }

    colnames(sc) <- c(paste0("tau_c", seq_len(nth)), "cw", "cb")
    return(sc)
  })                     # nolint end
}

# cluster-wise scores at given values - no cache
lav_bvmix_2l_sc <- function(fit_y1 = NULL, fit_y2 = NULL,
                            cw = NULL, cb = NULL, ngh = 21L) {
  cache <- lav_bvmix_2l_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, ngh = ngh
  )
  cache$theta <- c(cw, cb)
  tmp <- lav_bvmix_2l_loglik_cache(cache = cache) # populates pw
  lav_bvmix_2l_sc_cache(cache = cache)
}

# total loglik at given values - no cache
lav_bvmix_2l_logl <- function(fit_y1 = NULL, fit_y2 = NULL,
                              cw = NULL, cb = NULL, ngh = 21L) {
  cache <- lav_bvmix_2l_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, ngh = ngh
  )
  cache$theta <- c(cw, cb)
  lav_bvmix_2l_loglik_cache(cache = cache)
}

# nlminb objective/gradient (see lav_uvbv_common.R)
lav_bvmix_2l_min_fns <- function() {
  lav_uvbv_min_fns(
    logl_fun = lav_bvmix_2l_loglik_cache,
    grad_fun = lav_bvmix_2l_grad_cache
  )
}
