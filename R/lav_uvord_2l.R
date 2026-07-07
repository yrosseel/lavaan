# the univariate two-level ordinal probit (random-intercept) model
#
#     y*_ij = xw_ij' pw + xb_j' pb + y_w,ij + b_j,
#     y_w,ij ~ N(0, 1),  b_j ~ N(0, vb)
#     y_ij = k  <=>  tau_k < y*_ij <= tau_{k+1}
#
# i.e. conditional on the cluster effect b_j:
#
#     P(y_ij = k | b_j) = Phi(tau_k+1 - eta_ij - b_j)
#                         - Phi(tau_k - eta_ij - b_j)
#
# with eta_ij = xw_ij' pw + xb_j' pb (cfr. the single-level convention in
# lav_uvord.R: the linear predictor enters with a minus sign relative to
# the thresholds). The within-level latent residual variance is fixed to
# 1 (theta parameterization); there is no mean parameter (absorbed by the
# thresholds). This is stage 1 of the stage-wise unrestricted two-level
# model for a categorical variable.
#
# YR 2026 (two-level WLS, phase 2; covariates added in phase 7)
#
# parameter vector: theta = c(tau_1, ..., tau_nth, vb[, pw, pb])
# the slope parameters are APPENDED, so that the positions of the
# thresholds and vb do not depend on the presence of covariates
#
# The cluster likelihood integrates over b_j with (standardized)
# Gauss-Hermite quadrature; the per-node conditional kernel reuses
# lav_uvord_z12() with offset = eta + b_q:
#
#     L_j = sum_q w_q prod_i P(y_ij | b_q),   b_q = sqrt(vb) x_q
#
# computed per cluster on the log scale (log-sum-exp). The analytic
# gradient is posterior-weighted: with pw_jq = w_q exp(ll_jq) / L_j,
#
#   d logL_j / d tau = sum_q pw_jq * sum_{i in j} (y1 p1 - y2 p2)/pi |_q
#   d logL_j / d vb  = sum_q pw_jq * [sum_{i in j} (p2 - p1)/pi |_q]
#                                  * b_q / (2 vb)
#   d logL_j / d pw  = sum_q pw_jq * sum_{i in j} x_ij (p2 - p1)/pi |_q
#
# (the standardized nodes x_q and weights w_q do not depend on vb).
#
# WITHOUT covariates the conditional kernel only depends on the data
# through the response CATEGORY, so all per-node quantities are computed
# per category and aggregated per cluster via the (nclusters x ncat)
# count matrix -- the cost is (nearly) independent of the number of
# observations. With covariates every observation has its own linear
# predictor, so the kernels are computed per observation.
#
# fit object (list) contract -- see lav_uvord_2l_fit():
#   theta, nth, th_idx, bvar_idx, wslope_idx, bslope_idx,
#   nexo_w, nexo_b, nexo, y, x_w, x_b, z, cluster_idx, ngh,
#   logl, converged

lav_uvord_2l_fit <- function(y = NULL,
                             x_w = NULL, # within-level covariates
                             x_b = NULL, # between-level covariates
                             cluster_idx = NULL,
                             ngh = 21L,
                             optim_method = "nlminb1",
                             control = list(),
                             output = "list") {
  # create cache environment
  cache <- lav_uvord_2l_init_cache(
    y = y, x_w = x_w, x_b = x_b, cluster_idx = cluster_idx, ngh = ngh
  )
  nexo <- cache$nexo

  # optim
  minfns <- lav_uvbv_optim_fns(
    optim_method, lav_uvord_2l_min_fns(covariates = (nexo > 0L))
  )

  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control_nlminb <- modifyList(control_nlminb, control)

  lower <- c(rep(-Inf, cache$nth), 1e-6, rep(-Inf, nexo))

  # try 1
  optim <- nlminb(
    start = cache$theta, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control_nlminb, lower = lower,
    cache = cache
  )

  # try 2: different starting value for vb
  if (optim$convergence != 0L) {
    start_2 <- cache$theta_start
    start_2[cache$nth + 1L] <- 0.5
    optim <- nlminb(
      start = start_2, objective = minfns$objective,
      gradient = NULL, hessian = NULL,
      control = control_nlminb, lower = lower,
      cache = cache
    )
  }

  if (optim$convergence != 0L) {
    lav_msg_warn(gettext(
      "estimation of the two-level ordinal probit model did not converge"))
  }

  if (output == "cache") {
    return(cache)
  }

  list(
    theta = optim$par,
    nth = cache$nth,
    th_idx = seq_len(cache$nth),
    bvar_idx = cache$nth + 1L,
    wslope_idx = if (cache$nexo_w > 0L) cache$nth + 1L +
                   seq_len(cache$nexo_w)
                 else integer(0L),
    bslope_idx = if (cache$nexo_b > 0L) cache$nth + 1L + cache$nexo_w +
                   seq_len(cache$nexo_b)
                 else integer(0L),
    nexo_w = cache$nexo_w, nexo_b = cache$nexo_b, nexo = nexo,
    y = cache$y,
    x_w = cache$x_w, x_b = cache$x_b, z = cache$z,
    cluster_idx = cache$cluster_idx,
    ngh = ngh,
    logl = {
      cache$theta <- optim$par
      if (nexo > 0L) {
        lav_uvord_2l_x_loglik_cache(cache = cache)
      } else {
        lav_uvord_2l_loglik_cache(cache = cache)
      }
    },
    converged = (optim$convergence == 0L)
  )
}

# prepare cache environment
lav_uvord_2l_init_cache <- function(y = NULL,
                                    x_w = NULL,
                                    x_b = NULL,
                                    cluster_idx = NULL,
                                    ngh = 21L,
                                    parent = parent.frame()) {
  # y
  if (!is.integer(y)) {
    y <- as.integer(y)
  }
  if (anyNA(y)) {
    lav_msg_stop(gettext(
      "missing values are not supported (yet) in lav_uvord_2l_fit()."))
  }
  if (!min(y) == 1L) {
    y <- as.integer(ordered(y))
  }

  nobs <- length(y)
  y_freq <- tabulate(y)
  ncat <- length(y_freq)
  nth <- ncat - 1L
  y_prop <- y_freq / sum(y_freq)

  nexo_w <- if (is.null(x_w)) 0L else NCOL(x_w)
  nexo_b <- if (is.null(x_b)) 0L else NCOL(x_b)
  nexo <- nexo_w + nexo_b

  # cluster x category count matrix
  cl_sorted <- sort(unique(cluster_idx))
  cl_pos <- match(cluster_idx, cl_sorted)
  nclusters <- length(cl_sorted)
  cmat <- unclass(table(
    factor(cl_pos, levels = seq_len(nclusters)),
    factor(y, levels = seq_len(ncat))
  ))
  dimnames(cmat) <- NULL
  storage.mode(cmat) <- "double"

  # category-level offsets and indicator matrices
  ycat <- seq_len(ncat)
  o12 <- lav_uvord_o12(y = ycat, nth = nth)
  y1 <- matrix(seq_len(nth), ncat, nth, byrow = TRUE) == ycat
  y2 <- matrix(seq_len(nth), ncat, nth, byrow = TRUE) == (ycat - 1L)

  # standardized GH nodes
  gh <- lav_2l_gh_xw(ngh = ngh)

  # starting values: marginal thresholds (and slopes) live on the
  # sd = sqrt(1 + vb) scale; take vb0 = 0.2
  vb0 <- 0.2
  s0 <- sqrt(1 + vb0)
  if (nexo == 0L) {
    th_start <- qnorm(cumsum(y_prop[-length(y_prop)])) * s0
    theta <- c(th_start, vb0)
    z <- NULL
    o1_i <- o2_i <- NULL
  } else {
    z <- cbind(x_w, x_b)
    storage.mode(z) <- "double"
    if (anyNA(z)) {
      lav_msg_stop(gettext(
        "missing values are not supported (yet) in lav_uvord_2l_fit()."))
    }
    # single-level probit fit for the starting values
    fit0 <- lav_uvord_fit(y = y, x = z)
    theta <- c(fit0$theta[fit0$th_idx] * s0, vb0,
               fit0$theta[fit0$slope_idx] * s0)
    # observation-level offsets (for the per-observation kernels)
    o12_i <- lav_uvord_o12(y = y, nth = nth)
    o1_i <- o12_i$o1
    o2_i <- o12_i$o2
  }

  list2env(
    list(
      y = y, x_w = x_w, x_b = x_b, z = z,
      cluster_idx = cluster_idx, cmat = cmat,
      ycat = ycat, ncat = ncat,
      o1 = o12$o1, o2 = o12$o2, y1 = y1, y2 = y2,
      o1_i = o1_i, o2_i = o2_i,
      nexo_w = nexo_w, nexo_b = nexo_b, nexo = nexo,
      nobs = nobs, n = nobs, nth = nth, nclusters = nclusters,
      gh = gh, ngh = ngh,
      theta = theta, theta_start = theta
    ),
    parent = parent
  )
}

# per-node conditional cell probabilities (with the same numerical guards
# as the single-level lav_uvord_loglik_cache)
lav_uvord_2l_pi_q <- function(y = NULL, th = NULL, o1 = NULL, o2 = NULL,
                              b_q = NULL) {
  z12 <- lav_uvord_z12(y = y, th = th, o1 = o1, o2 = o2, offset = b_q)
  z1 <- z12$z1
  z2 <- z12$z2
  pi_i <- pnorm(z1) - pnorm(z2)
  large_idx <- which(z2 > 1)
  if (length(large_idx) > 0L) {
    pi_i[large_idx] <- (pnorm(z2[large_idx], lower.tail = FALSE) -
      pnorm(z1[large_idx], lower.tail = FALSE))
  }
  # guard against zero cells
  pi_i[pi_i < .Machine$double.eps] <- .Machine$double.eps
  list(pi_i = pi_i, z1 = z1, z2 = z2)
}

lav_uvord_2l_loglik_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    th <- theta[1:nth]
    vb <- theta[nth + 1L]
    s_b <- sqrt(vb)

    # category x node conditional log-probabilities
    logp <- matrix(0, ncat, ngh)
    for (q in seq_len(ngh)) {
      b_q <- s_b * gh$x[q]
      pq <- lav_uvord_2l_pi_q(y = ycat, th = th, o1 = o1, o2 = o2,
                              b_q = b_q)
      logp[, q] <- log(pq$pi_i)
    }

    # cluster x node log-likelihood contributions
    ll <- cmat %*% logp

    # log-sum-exp per cluster
    lwll <- sweep(ll, 2L, gh$logw, "+")
    mmax <- apply(lwll, 1L, max)
    logl_j <- mmax + log(rowSums(exp(lwll - mmax)))

    # posterior node weights (for the gradient)
    pw <- exp(lwll - logl_j)

    loglik <- sum(logl_j)
    return(loglik)
  })                     # nolint end
}

lav_uvord_2l_grad_cache <- function(cache = NULL) {
  colSums(lav_uvord_2l_sc_cache(cache = cache))
}

# cluster-wise scores (of the loglik); assumes the loglik has been
# evaluated at cache$theta (pw available)
lav_uvord_2l_sc_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    th <- theta[1:nth]
    vb <- theta[nth + 1L]
    s_b <- sqrt(vb)

    sc <- matrix(0, nclusters, nth + 1L)
    for (q in seq_len(ngh)) {
      b_q <- s_b * gh$x[q]
      pq <- lav_uvord_2l_pi_q(y = ycat, th = th, o1 = o1, o2 = o2,
                              b_q = b_q)
      p1 <- dnorm(pq$z1)
      p2 <- dnorm(pq$z2)
      inv_pi <- 1 / pq$pi_i

      # category-level score pieces -> cluster level via the count matrix
      # thresholds: (y1 p1 - y2 p2)/pi
      sc_th_q <- cmat %*% ((y1 * p1 - y2 * p2) * inv_pi)
      # cluster intercept: (p2 - p1)/pi
      g_b_q <- cmat %*% ((p2 - p1) * inv_pi)

      sc[, 1:nth] <- sc[, 1:nth] + pw[, q] * sc_th_q
      sc[, nth + 1L] <- sc[, nth + 1L] +
        pw[, q] * g_b_q * (b_q / (2 * vb))
    }

    return(sc)
  })                     # nolint end
}

# --- with covariates (per-observation kernels) ---

lav_uvord_2l_x_loglik_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    th <- theta[1:nth]
    vb <- theta[nth + 1L]
    s_b <- sqrt(vb)
    beta <- theta[nth + 1L + seq_len(nexo)]
    eta <- drop(z %*% beta)

    # cluster x node log-likelihood contributions
    ll <- matrix(0, nclusters, ngh)
    for (q in seq_len(ngh)) {
      pq <- lav_uvord_2l_pi_q(y = y, th = th, o1 = o1_i, o2 = o2_i,
                              b_q = eta + s_b * gh$x[q])
      ll[, q] <- drop(rowsum.default(log(pq$pi_i), group = cluster_idx,
                                     reorder = TRUE))
    }

    # log-sum-exp per cluster
    lwll <- sweep(ll, 2L, gh$logw, "+")
    mmax <- apply(lwll, 1L, max)
    logl_j <- mmax + log(rowSums(exp(lwll - mmax)))

    # posterior node weights (for the gradient)
    pw <- exp(lwll - logl_j)

    loglik <- sum(logl_j)
    return(loglik)
  })                     # nolint end
}

lav_uvord_2l_x_grad_cache <- function(cache = NULL) {
  colSums(lav_uvord_2l_x_sc_cache(cache = cache))
}

# cluster-wise scores, covariate case; assumes the loglik has been
# evaluated at cache$theta (pw available);
# columns: (tau_1..nth, vb, slopes)
lav_uvord_2l_x_sc_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    th <- theta[1:nth]
    vb <- theta[nth + 1L]
    s_b <- sqrt(vb)
    beta <- theta[nth + 1L + seq_len(nexo)]
    eta <- drop(z %*% beta)

    # z1 involves tau_y (y <= nth), z2 involves tau_{y-1} (y >= 2)
    idx1 <- which(y <= nth)
    idx2 <- which(y >= 2L)

    sc <- matrix(0, nclusters, nth + 1L + nexo)
    for (q in seq_len(ngh)) {
      pq <- lav_uvord_2l_pi_q(y = y, th = th, o1 = o1_i, o2 = o2_i,
                              b_q = eta + s_b * gh$x[q])
      p1 <- dnorm(pq$z1)
      p2 <- dnorm(pq$z2)
      inv_pi <- 1 / pq$pi_i

      # observation-level score pieces
      # thresholds: d logp_i / d tau_k = (p1 1{y=k} - p2 1{y=k+1}) / pi
      kern_th <- matrix(0, nobs, nth)
      kern_th[cbind(idx1, y[idx1])] <- (p1 * inv_pi)[idx1]
      kern_th[cbind(idx2, y[idx2] - 1L)] <-
        kern_th[cbind(idx2, y[idx2] - 1L)] - (p2 * inv_pi)[idx2]
      # offset (cluster intercept/slopes): (p2 - p1)/pi
      g_i <- (p2 - p1) * inv_pi

      cl <- rowsum.default(cbind(kern_th, g_i, g_i * z),
                           group = cluster_idx, reorder = TRUE)
      sc[, 1:nth] <- sc[, 1:nth] + pw[, q] * cl[, 1:nth, drop = FALSE]
      sc[, nth + 1L] <- sc[, nth + 1L] +
        pw[, q] * cl[, nth + 1L] * (s_b * gh$x[q] / (2 * vb))
      sc[, nth + 1L + seq_len(nexo)] <- sc[, nth + 1L + seq_len(nexo)] +
        pw[, q] * cl[, nth + 1L + seq_len(nexo), drop = FALSE]
    }

    return(sc)
  })                     # nolint end
}

# cluster-wise scores at given (or fitted) parameter values - no cache
lav_uvord_2l_sc <- function(fit_y = NULL, theta = NULL) {
  cache <- lav_uvord_2l_init_cache(
    y = fit_y$y, x_w = fit_y$x_w, x_b = fit_y$x_b,
    cluster_idx = fit_y$cluster_idx, ngh = fit_y$ngh
  )
  if (is.null(theta)) {
    theta <- fit_y$theta
  }
  cache$theta <- theta
  if (cache$nexo > 0L) {
    tmp <- lav_uvord_2l_x_loglik_cache(cache = cache) # populates pw
    sc <- lav_uvord_2l_x_sc_cache(cache = cache)
    colnames(sc) <- c(
      paste0("th", seq_len(cache$nth)), "bvar",
      lav_uvreg_2l_slope_names(cache$nexo_w, cache$nexo_b)
    )
  } else {
    tmp <- lav_uvord_2l_loglik_cache(cache = cache) # populates pw
    sc <- lav_uvord_2l_sc_cache(cache = cache)
    colnames(sc) <- c(paste0("th", seq_len(cache$nth)), "bvar")
  }
  sc
}

# total loglik - no cache
lav_uvord_2l_logl <- function(fit_y = NULL, theta = NULL) {
  cache <- lav_uvord_2l_init_cache(
    y = fit_y$y, x_w = fit_y$x_w, x_b = fit_y$x_b,
    cluster_idx = fit_y$cluster_idx, ngh = fit_y$ngh
  )
  if (is.null(theta)) {
    theta <- fit_y$theta
  }
  cache$theta <- theta
  if (cache$nexo > 0L) {
    lav_uvord_2l_x_loglik_cache(cache = cache)
  } else {
    lav_uvord_2l_loglik_cache(cache = cache)
  }
}

# reject descending thresholds (cfr. lav_uvord_pre_objective)
lav_uvord_2l_pre_objective <- function(x, cache = NULL) {
  if (cache$nth > 1L && x[1] > x[2]) {
    return(+Inf)
  }
  if (cache$nth > 2L && x[2] > x[3]) {
    return(+Inf)
  }
  if (cache$nth > 3L && x[3] > x[4]) {
    return(+Inf)
  }
  NULL
}

# nlminb objective/gradient (see lav_uvbv_common.R)
lav_uvord_2l_min_fns <- function(covariates = FALSE) {
  if (covariates) {
    lav_uvbv_min_fns(
      logl_fun = lav_uvord_2l_x_loglik_cache,
      grad_fun = lav_uvord_2l_x_grad_cache,
      pre_objective = lav_uvord_2l_pre_objective
    )
  } else {
    lav_uvbv_min_fns(
      logl_fun = lav_uvord_2l_loglik_cache,
      grad_fun = lav_uvord_2l_grad_cache,
      pre_objective = lav_uvord_2l_pre_objective
    )
  }
}
