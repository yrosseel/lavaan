# the univariate ordinal (probit/logit) regression model
#
# - probit regression / ordinal probit regression
# - logit regression / ordinal logit regression (logistic = TRUE)
# - using case weights wt
#
# YR 25 Nov 2019 (replacing the old lav_probit.R routines)
# YR/LDW 2026: refactored (see lav_uvbv_common.R for the shared cache/minfns
#              design); no change in behavior
#
# Note: the idea of using 'o1' and 'o2' when computing z1/z2 comes from
#       the dissertation of Christensen, 2012 (see also his `ordinal' package)
#
# The model for case i with response category y_i = k (k = 1, .., nth + 1):
#
#     P(y_i = k) = pfun(z1_i) - pfun(z2_i)
#
# with the latent-response 'rectangle' bounds
#
#     z1_i = th[k]     - eta_i  (+100 if k == nth + 1, i.e. th[nth+1] = +Inf)
#     z2_i = th[k - 1] - eta_i  (-100 if k == 1,       i.e. th[0]     = -Inf)
#
# and linear predictor eta_i = X_i %*% beta (0 if no covariates).
# The parameter vector is theta = c(th, beta).
#
# fit object (list) contract -- see lav_uvord_fit():
#   theta, nexo, nth, th_idx, slope_idx, missing_idx,
#   y (integer, recoded to 1..ncat), wt, y1, y2 (case x th indicator
#   matrices), z1, z2 (bounds at theta), x
# The z1/z2 values are those of the *last* evaluation inside nlminb;
# lav_uvord_update_fit() recomputes them for new parameter values.

lav_uvord_fit <- function(y = NULL,
                          x = NULL,
                          wt = rep(1, length(y)),
                          lower = -Inf,
                          upper = +Inf,
                          optim_method = "nlminb",
                          logistic = FALSE, # probit is the default
                          control = list(),
                          output = "list") {
  # y
  if (!is.integer(y)) {
    # brute force, no checking! (this is a lower-level function)
    y <- as.integer(y)
  }
  if (!min(y, na.rm = TRUE) == 1L) {
    y <- as.integer(ordered(y))
  }

  # check weights
  wt <- lav_uvbv_wt_check(y = y, wt = wt)

  # optim.method
  minfns <- lav_uvbv_optim_fns(optim_method, lav_uvord_min_fns())
  y_ncat <- length(tabulate(y)) # number of response categories
  if (y_ncat <= 1L) {
    minfns$gradient <- minfns$hessian <- NULL
  }

  # create cache environment
  cache <- lav_uvord_init_cache(y = y, x = x, wt = wt, logistic = logistic)

  # optimize -- only changes from defaults
  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control_nlminb <- modifyList(control_nlminb, control)

  optim <- nlminb(
    start = cache$theta, objective = minfns$objective,
    gradient = minfns$gradient, hessian = minfns$hessian,
    control = control_nlminb, lower = lower, upper = upper,
    cache = cache
  )

  if (output == "cache") {
    return(cache)
  }

  # return results as a list (to be compatible with lav_bvord.R)
  list(
    theta = optim$par,
    nexo = cache$nexo,
    nth = cache$nth,
    th_idx = seq_len(cache$nth),
    slope_idx = seq_along(optim$par)[-seq_len(cache$nth)],
    missing_idx = cache$missing_idx,
    y = cache$y,
    wt = cache$wt,
    y1 = cache$y1,
    y2 = cache$y2,
    z1 = cache$z1,
    z2 = cache$z2,
    x = cache$x
  )
}

# shortcut to get (possibly weighted) thresholds only, if no eXo
lav_uvord_th <- function(y = NULL, wt = NULL) {
  y_freq <- tabulate(y) # unweighted
  y_ncat <- length(y_freq) # number of response categories
  if (is.null(wt)) {
    y_prop <- y_freq / sum(y_freq)
  } else {
    y_freq <- numeric(y_ncat) # numeric! weights...
    for (cat in seq_len(y_ncat)) {
      y_freq[cat] <- sum(wt[y == cat], na.rm = TRUE)
    }
    y_prop <- y_freq / sum(y_freq)
  }

  qnorm(cumsum(y_prop[-length(y_prop)]))
}

# the -Inf/+Inf 'offsets': add +100 (-100) to z1 (z2) for the extreme
# categories, so that pfun() evaluates to 1 (0) there
lav_uvord_o12 <- function(y = NULL, nth = NULL) {
  list(
    o1 = ifelse(y == nth + 1, 100, 0),
    o2 = ifelse(y == 1, -100, 0)
  )
}

# the latent-response rectangle bounds z1/z2, given thresholds th, (optional)
# slopes beta for covariates x, and an (optional) additional per-case offset
# (the offset is not used by the single-level code, but allows a cluster
#  random effect to shift eta in a two-level extension)
lav_uvord_z12 <- function(y = NULL, th = NULL, o1 = NULL, o2 = NULL,
                          x = NULL, beta = NULL, offset = NULL) {
  th_1 <- c(0, th, 0) # th[0] and th[nth+1]; o1/o2 take over there
  eta <- NULL
  if (!is.null(x) && length(beta) > 0L) {
    eta <- drop(x %*% beta)
  }
  if (!is.null(offset)) {
    eta <- if (is.null(eta)) offset else eta + offset
  }
  if (!is.null(eta)) {
    z1 <- th_1[y + 1L] - eta + o1
    z2 <- th_1[y] - eta + o2
  } else {
    z1 <- th_1[y + 1L] + o1
    z2 <- th_1[y] + o2
  }
  list(z1 = z1, z2 = z2)
}

# prepare cache environment
lav_uvord_init_cache <- function(y = NULL,
                                 x = NULL,
                                 wt = rep(1, length(y)),
                                 logistic = FALSE,
                                 parent = parent.frame()) {
  nobs <- length(y)

  # number of response categories
  y_ncat <- length(tabulate(y)) # unweighted
  # number of thresholds
  nth <- y_ncat - 1L

  # X
  if (is.null(x)) {
    nexo <- 0L
  } else {
    x <- unname(x)
    nexo <- ncol(x)
    # new in 0.6-17: check if X is full rank
    if (!anyNA(x)) {
      if (qr(x)$rank < ncol(x)) {
        lav_msg_stop(gettext(
          "matrix of exogenous covariates is rank deficient! (i.e., some x
          variables contain redundant information)"))
      }
    }
  }

  # nobs
  if (is.null(wt)) {
    n <- nobs
  } else {
    n <- sum(wt)
  }

  # frequencies (possibly weighted by wt)
  y_freq <- numeric(y_ncat) # numeric! weights...
  for (cat in seq_len(y_ncat)) {
    y_freq[cat] <- sum(wt[y == cat], na.rm = TRUE)
  }
  y_prop <- y_freq / sum(y_freq)

  # missing values
  missing_idx <- which(is.na(y))

  # missing values
  if (any(is.na(y)) || (!is.null(x) && any(is.na(x)))) {
    lav_crossprod <- lav_mat_crossprod
  } else {
    lav_crossprod <- base::crossprod
  }

  # distribution
  if (logistic) {
    pfun <- plogis
    dfun <- dlogis
    gfun <- function(x) {
      # FIXMe: is it worth making this work for abs(x) > 200?
      out <- numeric(length(x))
      out[is.na(x)] <- NA
      x_ok <- which(abs(x) < 200)
      e <- exp(-x[x_ok])
      e1 <- 1 + e
      e2 <- e1 * e1
      e4 <- e2 * e2
      out[x_ok] <- -e / e2 + e * (2 * (e * e1)) / e4
      out
    }
  } else {
    pfun <- pnorm
    dfun <- dnorm
    gfun <- function(x) {
      -x * dnorm(x)
    }
  }

  # offsets -Inf/+Inf
  o12 <- lav_uvord_o12(y = y, nth = nth)
  o1 <- o12$o1
  o2 <- o12$o2

  # TH matrices (Matrix logical?)
  if (nth > 0L) {
    y1 <- matrix(1:nth, nobs, nth, byrow = TRUE) == y
    y2 <- matrix(1:nth, nobs, nth, byrow = TRUE) == (y - 1L)
  } else {
    y1 <- y2 <- matrix(nrow = nobs, ncol = 0)
  }

  # starting values
  if (nexo == 0L && nth > 0L) {
    if (logistic) {
      th_start <- qlogis(cumsum(y_prop[-length(y_prop)]))
    } else {
      th_start <- qnorm(cumsum(y_prop[-length(y_prop)]))
    }
  } else if ((nth == 1L && nexo > 0L) || nth == 0L) {
    th_start <- 0
  } else {
    if (logistic) {
      th_start <- qlogis((1:nth) / (nth + 1))
    } else {
      th_start <- qnorm((1:nth) / (nth + 1))
    }
  }
  beta_start <- rep(0, nexo)
  theta <- c(th_start, beta_start)

  list2env(
    list(
      y = y, x = x, wt = wt, o1 = o1, o2 = o2,
      missing_idx = missing_idx, n = n,
      pfun = pfun, dfun = dfun, gfun = gfun,
      lav_crossprod = lav_crossprod,
      nth = nth, nobs = nobs, y.ncat = y_ncat, nexo = nexo,
      y1 = y1, y2 = y2,
      theta = theta
    ),
    parent = parent
  )
}

# compute total (log)likelihood
lav_uvord_loglik <- function(y = NULL,
                             x = NULL,
                             wt = rep(1, length(y)),
                             logistic = FALSE,
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvord_fit(
      y = y, x = x, wt = wt,
      logistic = logistic, output = "cache"
    )
  }
  lav_uvord_loglik_cache(cache = cache)
}

lav_uvord_loglik_cache <- function(cache = NULL) {
  with(cache, {      # nolint start
    # Note: we could treat the binary case separately,
    # avoiding calling pfun() twice

    # free parameters
    th <- theta[1:nth]
    beta <- theta[-c(1:nth)]
    z12 <- lav_uvord_z12(y = y, th = th, o1 = o1, o2 = o2, x = x, beta = beta)
    z1 <- z12$z1
    z2 <- z12$z2
    pi_i <- pfun(z1) - pfun(z2)

    # avoid numerical degradation if z2 (and therefore z1) are both 'large'
    # and the pfuns are close to 1.0
    large_idx <- which(z2 > 1)
    if (length(large_idx) > 0L) {
      pi_i[large_idx] <- (pfun(z2[large_idx], lower.tail = FALSE) -
        pfun(z1[large_idx], lower.tail = FALSE))
    }

    loglik <- sum(wt * log(pi_i), na.rm = TRUE)

    return(loglik)
  })                      # nolint end
}

# casewise scores
lav_uvord_sc <- function(y = NULL,
                             x = NULL,
                             wt = rep(1, length(y)),
                             use_weights = TRUE,
                             logistic = FALSE,
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvord_fit(
      y = y, x = x, wt = wt,
      logistic = logistic, output = "cache"
    )
  }
  sc <- lav_uvord_sc_cache(cache = cache)

  if (!is.null(wt) && use_weights) {
    sc <- sc * wt
  }

  sc
}

lav_uvord_sc_cache <- function(cache = NULL) {
  with(cache, {
    # d logl / d pi
    dldpi <- 1 / pi_i # unweighted!

    # we assume z1/z2 are available
    p1 <- dfun(z1)
    p2 <- dfun(z2)

    # th
    scores_th <- dldpi * (y1 * p1 - y2 * p2)

    # beta
    if (nexo > 0L) {
      scores_beta <- dldpi * (-x) * (p1 - p2)
      return(cbind(scores_th, scores_beta, deparse.level = 0))
    } else {
      return(scores_th)
    }
  })
}

lav_uvord_grad_cache <- function(cache = NULL) {
  with(cache, {
    # d logl / d pi
    wtp <- wt / pi_i
    p1 <- dfun(z1)
    p2 <- dfun(z2)

    # th
    dxa <- y1 * p1 - y2 * p2
    scores_th <- wtp * dxa

    # beta
    if (nexo > 0L) {
      dxb <- x * (p1 - p2) # == x*p1 - x*p2
      scores_beta <- wtp * (-dxb)
      return(colSums(cbind(scores_th, scores_beta, deparse.level = 0),
        na.rm = TRUE
      ))
    } else {
      return(colSums(scores_th, na.rm = TRUE))
    }
  })
}


# compute total Hessian
lav_uvord_hessian <- function(y = NULL,
                              x = NULL,
                              wt = rep(1, length(y)),
                              logistic = FALSE,
                              cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvord_fit(
      y = y, x = x, wt = wt,
      logistic = logistic, output = "cache"
    )
  }
  tmp <- lav_uvord_loglik_cache(cache = cache)
  tmp <- lav_uvord_grad_cache(cache = cache)
  lav_uvord_hessian_cache(cache = cache)
}

lav_uvord_hessian_cache <- function(cache = NULL) {
  with(cache, {          # nolint start
    wtp2 <- wt / (pi_i * pi_i)

    g1w <- gfun(z1) * wtp
    g2w <- gfun(z2) * wtp

    y1gw <- y1 * g1w
    y2gw <- y2 * g2w

    dx2_tau <- (lav_crossprod(y1gw, y1) - lav_crossprod(y2gw, y2) -
      lav_crossprod(dxa, dxa * wtp2))

    if (nexo == 0L) {
      return(dx2_tau)
    }

    dxb2 <- dxb * wtp2
    dx2_beta <- (lav_crossprod(x * g1w, x) - lav_crossprod(x * g2w, x) -
      lav_crossprod(dxb, dxb2))

    dx_taubeta <- (-lav_crossprod(y1gw, x) + lav_crossprod(y2gw, x) +
      lav_crossprod(dxa, dxb2))

    hessian <- rbind(cbind(dx2_tau, dx_taubeta, deparse.level = 0),
      cbind(t(dx_taubeta), dx2_beta, deparse.level = 0),
      deparse.level = 0
    )

    return(hessian)
  })                             # nolint end
}

# reject descending values for the first few thresholds (new in 0.6-8)
lav_uvord_pre_objective <- function(x, cache = NULL) {
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

# nlminb objective/gradient/hessian (see lav_uvbv_common.R)
lav_uvord_min_fns <- function() {
  lav_uvbv_min_fns(
    logl_fun = lav_uvord_loglik_cache,
    grad_fun = lav_uvord_grad_cache,
    hessian_fun = lav_uvord_hessian_cache,
    pre_objective = lav_uvord_pre_objective
  )
}

# get 'z1' and 'z2' values, given (new) values for the parameters
# only needed for lav_bvord_cor_sc(), which is called from
# lav_pml_dploglik_dimplied() in lav_model_gradient_pml.R
lav_uvord_update_fit <- function(fit_y = NULL, th_new = NULL, sl_new = NULL) {
  # return fit.y with 'update' z1/z2 values
  if (is.null(th_new) && is.null(sl_new)) {
    return(fit_y)
  }

  if (!is.null(th_new)) {
    fit_y$theta[fit_y$th_idx] <- th_new
  }
  if (!is.null(sl_new)) {
    fit_y$theta[fit_y$slope_idx] <- sl_new
  }

  nth <- length(fit_y$th_idx)
  o12 <- lav_uvord_o12(y = fit_y$y, nth = nth)

  theta <- fit_y$theta
  th <- theta[1:nth]
  beta <- theta[-c(1:nth)]
  z12 <- lav_uvord_z12(
    y = fit_y$y, th = th, o1 = o12$o1, o2 = o12$o2,
    x = fit_y$x, beta = beta
  )
  fit_y$z1 <- z12$z1
  fit_y$z2 <- z12$z2

  fit_y
}
