# the bivariate two-level ordinal probit model, univariate parameters fixed
#
#     y*_a,ij = w_a,ij + b_a,j     y*_c,ij = w_c,ij + b_c,j
#
#     (w_a, w_c) ~ N(0, | 1     rho_w |)    (b_a, b_c) ~ N(0, | vb_a  cb   |)
#                       | rho_w 1     |                       | cb    vb_c |
#
# YR 2026 (two-level WLS, phase 3)
#
# This is stage 2 of the stage-wise unrestricted two-level model for a
# pair of ordinal variables: the univariate parameters (thresholds tau_a,
# tau_c and between variances vb_a, vb_c) are held fixed at their stage-1
# values (lav_uvord_2l_fit), and theta = c(rho_w, cb) -- the WITHIN
# polychoric correlation and the BETWEEN covariance -- is estimated by ML
# with two-dimensional Gauss-Hermite quadrature over the cluster effects:
#
#     b_a = s_a u1,   b_c = (cb/s_a) u1 + w u2,  w = sqrt(vb_c - cb^2/vb_a)
#
# (the Cholesky of Sigma_b; u1, u2 iid standard normal nodes). The
# per-node conditional kernel is the shifted-threshold bivariate probit
# rectangle probability (pbinorm), exactly as in the exo branch of the
# single-level lav_bvord.R.
#
# Computations are organized per response PATTERN (category pair): the
# conditional kernel takes only ncat_a * ncat_c distinct values per node,
# and the cluster log-likelihoods follow from the (nclusters x npatterns)
# count matrix -- the cost is (nearly) independent of the number of
# observations.
#
# lav_bvord_2l_sc() returns cluster-wise scores with respect to ALL
# parameters (tau_a, tau_c, vb_a, vb_c, rho_w, cb): the fixed-parameter
# columns are the cross-derivatives needed for the stage-wise sandwich.

lav_bvord_2l_cor_twostep_fit <- function(fit_y1 = NULL,
                                         fit_y2 = NULL,
                                         ngh = 13L,
                                         optim_method = "nlminb1",
                                         control = list(),
                                         y1_name = NULL, y2_name = NULL,
                                         output = "theta") {
  # create cache environment
  cache <- lav_bvord_2l_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, ngh = ngh
  )

  # bounds
  cb_bound <- 0.999 * sqrt(cache$vb1 * cache$vb2)
  cb_fixed <- (cb_bound < 1e-08)
  if (cb_fixed) {
    cache$theta[2] <- 0
    lower <- c(-0.999, 0)
    upper <- c(+0.999, 0)
  } else {
    lower <- c(-0.999, -cb_bound)
    upper <- c(+0.999, +cb_bound)
  }

  # optim
  minfns <- lav_uvbv_optim_fns(optim_method, lav_bvord_2l_min_fns())

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
        "estimation of the two-level polychoric parameters did not converge
        for variables %1$s and %2$s", y1_name, y2_name))
    } else {
      lav_msg_warn(gettext(
        "estimation of the two-level polychoric parameters did not always
        converge"))
    }
  }

  if (output == "cache") {
    cache$theta <- optim$par
    return(cache)
  }

  # theta = c(rho_w, cb)
  optim$par
}

# prepare cache environment
# fit_y1/fit_y2: lav_uvord_2l_fit() objects of the two variables
lav_bvord_2l_init_cache <- function(fit_y1 = NULL,
                                    fit_y2 = NULL,
                                    ngh = 13L,
                                    parent = parent.frame()) {
  stopifnot(length(fit_y1$y) == length(fit_y2$y))
  y1 <- fit_y1$y
  y2 <- fit_y2$y
  cluster_idx <- fit_y1$cluster_idx
  nclusters <- length(unique(cluster_idx))
  nobs <- length(y1)

  # fixed univariate parameters
  th1 <- fit_y1$theta[fit_y1$th_idx]
  vb1 <- fit_y1$theta[fit_y1$bvar_idx]
  th2 <- fit_y2$theta[fit_y2$th_idx]
  vb2 <- fit_y2$theta[fit_y2$bvar_idx]

  nth1 <- length(th1)
  nth2 <- length(th2)
  ncat1 <- nth1 + 1L
  ncat2 <- nth2 + 1L
  npat <- ncat1 * ncat2

  # response patterns: pattern p = (k, l) at index k + (l-1)*ncat1
  yp1 <- rep(seq_len(ncat1), times = ncat2)
  yp2 <- rep(seq_len(ncat2), each = ncat1)

  # cluster x pattern count matrix
  pid <- y1 + (y2 - 1L) * ncat1
  cl_sorted <- sort(unique(cluster_idx))
  cl_pos <- match(cluster_idx, cl_sorted)
  cmat <- unclass(table(
    factor(cl_pos, levels = seq_len(nclusters)),
    factor(pid, levels = seq_len(npat))
  ))
  dimnames(cmat) <- NULL
  storage.mode(cmat) <- "double"

  # pattern-level unshifted rectangle bounds
  o12_1 <- lav_uvord_o12(y = yp1, nth = nth1)
  o12_2 <- lav_uvord_o12(y = yp2, nth = nth2)
  z12_1 <- lav_uvord_z12(y = yp1, th = th1, o1 = o12_1$o1, o2 = o12_1$o2)
  z12_2 <- lav_uvord_z12(y = yp2, th = th2, o1 = o12_2$o1, o2 = o12_2$o2)

  # pattern-level threshold indicator matrices
  y1_1 <- matrix(seq_len(nth1), npat, nth1, byrow = TRUE) == yp1
  y1_2 <- matrix(seq_len(nth1), npat, nth1, byrow = TRUE) == (yp1 - 1L)
  y2_1 <- matrix(seq_len(nth2), npat, nth2, byrow = TRUE) == yp2
  y2_2 <- matrix(seq_len(nth2), npat, nth2, byrow = TRUE) == (yp2 - 1L)

  # 2-D tensor grid of standardized GH nodes
  gh <- lav_2l_gh_xw(ngh = ngh)
  u1 <- rep(gh$x, times = ngh)
  u2 <- rep(gh$x, each = ngh)
  logw <- rep(gh$logw, times = ngh) + rep(gh$logw, each = ngh)

  # starting values: crude total correlation, split within/between
  r_tot <- suppressWarnings(cor(y1, y2))
  if (!is.finite(r_tot)) {
    r_tot <- 0
  }
  rho0 <- max(min(r_tot, 0.80), -0.80)
  cb0 <- r_tot * sqrt(vb1 * vb2)
  cb0 <- max(min(cb0, 0.95 * sqrt(vb1 * vb2)), -0.95 * sqrt(vb1 * vb2))
  theta <- c(rho0, cb0)

  list2env(
    list(
      cluster_idx = cluster_idx,
      nobs = nobs, n = nobs, nclusters = nclusters, npat = npat,
      th1 = th1, th2 = th2, nth1 = nth1, nth2 = nth2,
      vb1 = vb1, vb2 = vb2, s1 = sqrt(vb1), s2 = sqrt(vb2),
      cmat = cmat,
      z1a0 = z12_1$z1, z2a0 = z12_1$z2,
      z1c0 = z12_2$z1, z2c0 = z12_2$z2,
      y1_1 = y1_1, y1_2 = y1_2, y2_1 = y2_1, y2_2 = y2_2,
      u1 = u1, u2 = u2, logw = logw, nq = ngh * ngh, ngh = ngh,
      theta = theta
    ),
    parent = parent
  )
}

# full-grid (pattern x node) rectangle bounds and probabilities; all
# quantities are stored in the cache (pattern index moves fastest)
lav_bvord_2l_grid_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    rho_w <- theta[1]
    cb <- theta[2]

    # Cholesky of Sigma_b
    l21 <- cb / s1
    l22 <- sqrt(max(vb2 - cb * cb / vb1, 1e-12))
    b_a <- s1 * u1
    b_c <- l21 * u1 + l22 * u2

    ba_g <- rep(b_a, each = npat)
    bc_g <- rep(b_c, each = npat)
    z1a_g <- rep(z1a0, times = nq) - ba_g
    z2a_g <- rep(z2a0, times = nq) - ba_g
    z1c_g <- rep(z1c0, times = nq) - bc_g
    z2c_g <- rep(z2c0, times = nq) - bc_g

    p_g <- pbinorm(
      upper_x = z1a_g, upper_y = z1c_g,
      lower_x = z2a_g, lower_y = z2c_g, rho = rho_w
    )
    p_g[p_g < .Machine$double.eps] <- .Machine$double.eps
    inv_p_g <- 1 / p_g

    return(invisible(NULL))
  })                     # nolint end
}

lav_bvord_2l_loglik_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    tmp <- lav_bvord_2l_grid_cache(cache)

    # cluster x node log-likelihood contributions
    ll <- cmat %*% matrix(log(p_g), npat, nq)

    lwll <- sweep(ll, 2L, logw, "+")
    mmax <- apply(lwll, 1L, max)
    logl_j <- mmax + log(rowSums(exp(lwll - mmax)))
    pw <- exp(lwll - logl_j)

    loglik <- sum(logl_j)
    return(loglik)
  })                     # nolint end
}

lav_bvord_2l_grad_cache <- function(cache = NULL) {
  sc <- lav_bvord_2l_sc_cache(cache = cache, free_only = TRUE)
  colSums(sc)
}

# cluster-wise scores of the loglik; assumes the loglik has been
# evaluated at cache$theta (grid + pw available)
# free_only = FALSE: columns tau_a (nth1), tau_c (nth2), vb_a, vb_c,
#                    rho_w, cb
# free_only = TRUE:  columns rho_w, cb (for the optimizer)
lav_bvord_2l_sc_cache <- function(cache = NULL, free_only = FALSE) {
  cache$free_only <- free_only
  with(cache, {          # nolint start
    rho_w <- theta[1]
    cb <- theta[2]
    r <- sqrt(1 - rho_w * rho_w)

    l21 <- cb / s1
    l22 <- sqrt(max(vb2 - cb * cb / vb1, 1e-12))

    # derivatives of the node positions w.r.t. (vb_a, vb_c, cb)
    # b_a = s1 u1; b_c = (cb/s1) u1 + l22 u2
    dbc_dcb <- u1 / s1 - u2 * cb / (vb1 * l22)

    # rectangle derivatives w.r.t. the c-bounds (needed for cb)
    d_z1c_g <- dnorm(z1c_g) * (pnorm((z1a_g - rho_w * z1c_g) / r) -
                               pnorm((z2a_g - rho_w * z1c_g) / r))
    d_z2c_g <- dnorm(z2c_g) * (pnorm((z1a_g - rho_w * z2c_g) / r) -
                               pnorm((z2a_g - rho_w * z2c_g) / r))
    g_c <- cmat %*% matrix(-(d_z1c_g - d_z2c_g) * inv_p_g, npat, nq)

    # within correlation: four-corner densities
    phi_rho_g <- (lav_dbinorm(z1a_g, z1c_g, rho_w) -
                  lav_dbinorm(z2a_g, z1c_g, rho_w) -
                  lav_dbinorm(z1a_g, z2c_g, rho_w) +
                  lav_dbinorm(z2a_g, z2c_g, rho_w))
    g_rho <- cmat %*% matrix(phi_rho_g * inv_p_g, npat, nq)

    sc_rho <- rowSums(pw * g_rho)
    sc_cb <- rowSums(pw * sweep(g_c, 2L, dbc_dcb, "*"))

    if (free_only) {
      sc <- cbind(sc_rho, sc_cb)
      colnames(sc) <- c("rho_w", "cb")
      return(sc)
    }

    # remaining node-position derivatives
    dba_dvba <- u1 / (2 * s1)
    dbc_dvba <- -u1 * cb / (2 * vb1 * s1) +
      u2 * (cb * cb / (vb1 * vb1)) / (2 * l22)
    dbc_dvbc <- u2 / (2 * l22)

    # rectangle derivatives w.r.t. the a-bounds
    d_z1a_g <- dnorm(z1a_g) * (pnorm((z1c_g - rho_w * z1a_g) / r) -
                               pnorm((z2c_g - rho_w * z1a_g) / r))
    d_z2a_g <- dnorm(z2a_g) * (pnorm((z1c_g - rho_w * z2a_g) / r) -
                               pnorm((z2c_g - rho_w * z2a_g) / r))
    g_a <- cmat %*% matrix(-(d_z1a_g - d_z2a_g) * inv_p_g, npat, nq)

    npar <- nth1 + nth2 + 4L
    sc <- matrix(0, nclusters, npar)

    # thresholds: dP/dtau_ak = y1a_k d_z1a - y2a_k d_z2a
    # (the pattern-level indicator columns recycle down the grid columns)
    m_z1a <- matrix(d_z1a_g * inv_p_g, npat, nq)
    m_z2a <- matrix(d_z2a_g * inv_p_g, npat, nq)
    m_z1c <- matrix(d_z1c_g * inv_p_g, npat, nq)
    m_z2c <- matrix(d_z2c_g * inv_p_g, npat, nq)
    for (k in seq_len(nth1)) {
      g_k <- cmat %*% (y1_1[, k] * m_z1a - y1_2[, k] * m_z2a)
      sc[, k] <- rowSums(pw * g_k)
    }
    for (k in seq_len(nth2)) {
      g_k <- cmat %*% (y2_1[, k] * m_z1c - y2_2[, k] * m_z2c)
      sc[, nth1 + k] <- rowSums(pw * g_k)
    }

    sc[, nth1 + nth2 + 1L] <-
      rowSums(pw * (sweep(g_a, 2L, dba_dvba, "*") +
                    sweep(g_c, 2L, dbc_dvba, "*")))
    sc[, nth1 + nth2 + 2L] <- rowSums(pw * sweep(g_c, 2L, dbc_dvbc, "*"))
    sc[, nth1 + nth2 + 3L] <- sc_rho
    sc[, nth1 + nth2 + 4L] <- sc_cb

    colnames(sc) <- c(
      paste0("tau_a", seq_len(nth1)), paste0("tau_c", seq_len(nth2)),
      "vb_a", "vb_c", "rho_w", "cb"
    )
    return(sc)
  })                     # nolint end
}

# cluster-wise scores at given values - no cache
lav_bvord_2l_sc <- function(fit_y1 = NULL, fit_y2 = NULL,
                            rho_w = NULL, cb = NULL, ngh = 13L) {
  cache <- lav_bvord_2l_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, ngh = ngh
  )
  cache$theta <- c(rho_w, cb)
  tmp <- lav_bvord_2l_loglik_cache(cache = cache) # populates pw
  lav_bvord_2l_sc_cache(cache = cache)
}

# total loglik at given values - no cache
lav_bvord_2l_logl <- function(fit_y1 = NULL, fit_y2 = NULL,
                              rho_w = NULL, cb = NULL, ngh = 13L) {
  cache <- lav_bvord_2l_init_cache(
    fit_y1 = fit_y1, fit_y2 = fit_y2, ngh = ngh
  )
  cache$theta <- c(rho_w, cb)
  lav_bvord_2l_loglik_cache(cache = cache)
}

# nlminb objective/gradient (see lav_uvbv_common.R)
lav_bvord_2l_min_fns <- function() {
  lav_uvbv_min_fns(
    logl_fun = lav_bvord_2l_loglik_cache,
    grad_fun = lav_bvord_2l_grad_cache
  )
}
