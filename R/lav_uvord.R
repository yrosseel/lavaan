# functions to deal with binary/ordinal univariate data

# - probit regression
# - ordinal probit regression
# - logit regression
# - ordinal logit regression

# Note: the idea of using 'o1' and 'o2' when computing z1/z2 comes from
#       the dissertation of Christensen, 2012 (see also his `ordinal' package)

# YR - 25 Nov 2019 (replacing the old lav_probit.R routines)

lav_uvord_fit <- function(y = NULL,
                          X = NULL,
                          wt = rep(1, length(y)),
                          lower = -Inf,
                          upper = +Inf,
                          optim.method = "nlminb",
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
  if (is.null(wt)) {
    wt <- rep(1, length(y))
  } else {
    if (length(y) != length(wt)) {
      lav_msg_stop(gettext("length y is not the same as length wt"))
    }
    if (any(wt < 0)) {
      lav_msg_stop(gettext("all weights should be positive"))
    }
  }

  # check lower/upper
  # TODO

  # optim.method
  minObjective <- lav_uvord_min_objective
  minGradient <- lav_uvord_min_gradient
  minHessian <- lav_uvord_min_hessian
  if (optim.method == "nlminb" || optim.method == "nlminb2") {
    # nothing to do
  } else if (optim.method == "nlminb0") {
    minGradient <- minHessian <- NULL
  } else if (optim.method == "nlminb1") {
    minHessian <- NULL
  }

  # create cache environment
  cache <- lav_uvord_init_cache(y = y, X = X, wt = wt, logistic = logistic)

  # optimize -- only changes from defaults
  control.nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = 0L, abs.tol = (.Machine$double.eps * 10)
  )
  control.nlminb <- modifyList(control.nlminb, control)

  optim <- nlminb(
    start = cache$theta, objective = minObjective,
    gradient = minGradient, hessian = minHessian,
    control = control.nlminb, lower = lower, upper = upper,
    cache = cache
  )

  if (output == "cache") {
    return(cache)
  }

  # return results as a list (to be compatible with lav_polychor.R)
  out <- list(
    theta = optim$par,
    nexo = cache$nexo,
    nth = cache$nth,
    th.idx = seq_len(cache$nth),
    slope.idx = seq_len(length(optim$par))[-seq_len(cache$nth)],
    missing.idx = cache$missing.idx,
    y = cache$y,
    wt = cache$wt,
    Y1 = cache$Y1,
    Y2 = cache$Y2,
    z1 = cache$z1,
    z2 = cache$z2,
    X = cache$X
  )
}

# shortcut to get (possibly weighted) thresholds only, if no eXo
lav_uvord_th <- function(y = NULL, wt = NULL) {
  y.freq <- tabulate(y) # unweighted
  y.ncat <- length(y.freq) # number of response categories
  if (is.null(wt)) {
    y.prop <- y.freq / sum(y.freq)
  } else {
    y.freq <- numeric(y.ncat) # numeric! weights...
    for (cat in seq_len(y.ncat)) {
      y.freq[cat] <- sum(wt[y == cat], na.rm = TRUE)
    }
    y.prop <- y.freq / sum(y.freq)
  }

  qnorm(cumsum(y.prop[-length(y.prop)]))
}


# prepare cache environment
lav_uvord_init_cache <- function(y = NULL,
                                 X = NULL,
                                 wt = rep(1, length(y)),
                                 logistic = FALSE,
                                 parent = parent.frame()) {
  nobs <- length(y)

  # number of response categories
  y.ncat <- length(tabulate(y)) # unweighted
  # number of thresholds
  nth <- y.ncat - 1L

  # X
  if (is.null(X)) {
    nexo <- 0L
  } else {
    X <- unname(X)
    nexo <- ncol(X)
    # new in 0.6-17: check if X is full rank
    if (!anyNA(X)) {
      if (qr(X)$rank < ncol(X)) {
        lav_msg_stop(gettext(
          "matrix of exogenous covariates is rank deficient!(i.e., some x
          variables contain redundant information)"))
      }
    }
  }

  # nobs
  if (is.null(wt)) {
    N <- nobs
  } else {
    N <- sum(wt)
  }

  # frequencies (possibly weighted by wt)
  y.freq <- numeric(y.ncat) # numeric! weights...
  for (cat in seq_len(y.ncat)) {
    y.freq[cat] <- sum(wt[y == cat], na.rm = TRUE)
  }
  y.prop <- y.freq / sum(y.freq)

  # missing values
  missing.idx <- which(is.na(y))

  # missing values
  if (any(is.na(y)) || (!is.null(X) && any(is.na(X)))) {
    lav_crossprod <- lav_matrix_crossprod
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
      x.ok <- which(abs(x) < 200)
      e <- exp(-x[x.ok])
      e1 <- 1 + e
      e2 <- e1 * e1
      e4 <- e2 * e2
      out[x.ok] <- -e / e2 + e * (2 * (e * e1)) / e4
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
  o1 <- ifelse(y == nth + 1, 100, 0)
  o2 <- ifelse(y == 1, -100, 0)

  # TH matrices (Matrix logical?)
  Y1 <- matrix(1:nth, nobs, nth, byrow = TRUE) == y
  Y2 <- matrix(1:nth, nobs, nth, byrow = TRUE) == (y - 1L)

  # starting values
  if (nexo == 0L) {
    if (logistic) {
      th.start <- qlogis(cumsum(y.prop[-length(y.prop)]))
    } else {
      th.start <- qnorm(cumsum(y.prop[-length(y.prop)]))
    }
  } else if (nth == 1L && nexo > 0L) {
    th.start <- 0
  } else {
    if (logistic) {
      # th.start <- seq(-1, 1, length = nth) / 2
      th.start <- qlogis((1:nth) / (nth + 1))
    } else {
      # th.start <- seq(-1, 1, length = nth) / 2
      th.start <- qnorm((1:nth) / (nth + 1))
    }
  }
  beta.start <- rep(0, nexo)
  theta <- c(th.start, beta.start)

  # parameter labels (for pretty output only)
  # th.lab <- paste("th", seq_len(nth), sep = "")
  # sl.lab <- character(0L)
  # if(nexo > 0L) {
  #    sl.lab <- paste("beta", seq_len(nexo), sep = "")
  # }
  # theta.labels <- c(th.lab, sl.lab)

  out <- list2env(
    list(
      y = y, X = X, wt = wt, o1 = o1, o2 = o2,
      missing.idx = missing.idx, N = N,
      pfun = pfun, dfun = dfun, gfun = gfun,
      lav_crossprod = lav_crossprod,
      nth = nth, nobs = nobs, y.ncat = y.ncat, nexo = nexo,
      Y1 = Y1, Y2 = Y2,
      theta = theta
    ),
    parent = parent
  )

  out
}

# compute total (log)likelihood
lav_uvord_loglik <- function(y = NULL,
                             X = NULL,
                             wt = rep(1, length(y)),
                             logistic = FALSE,
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvord_fit(
      y = y, X = X, wt = wt,
      logistic = logistic, output = "cache"
    )
  }
  lav_uvord_loglik_cache(cache = cache)
}

lav_uvord_loglik_cache <- function(cache = NULL) {
  with(cache, {
    # Note: we could treat the binary case separately,
    # avoiding calling pfun() twice

    # free parameters
    th <- theta[1:nth]
    TH <- c(0, th, 0)
    beta <- theta[-c(1:nth)]
    if (nexo > 0L) {
      eta <- drop(X %*% beta)
      z1 <- TH[y + 1L] - eta + o1
      z2 <- TH[y] - eta + o2
    } else {
      z1 <- TH[y + 1L] + o1
      z2 <- TH[y] + o2
    }
    pi.i <- pfun(z1) - pfun(z2)

    # avoid numerical degradation if z2 (and therefore z1) are both 'large'
    # and the pfuns are close to 1.0
    large.idx <- which(z2 > 1)
    if (length(large.idx) > 0L) {
      pi.i[large.idx] <- (pfun(z2[large.idx], lower.tail = FALSE) -
        pfun(z1[large.idx], lower.tail = FALSE))
    }

    loglik <- sum(wt * log(pi.i), na.rm = TRUE)

    return(loglik)
  })
}

# casewise scores
lav_uvord_scores <- function(y = NULL,
                             X = NULL,
                             wt = rep(1, length(y)),
                             use.weights = TRUE,
                             logistic = FALSE,
                             cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvord_fit(
      y = y, X = X, wt = wt,
      logistic = logistic, output = "cache"
    )
  }
  SC <- lav_uvord_scores_cache(cache = cache)

  if (!is.null(wt) && use.weights) {
    SC <- SC * wt
  }

  SC
}

lav_uvord_scores_cache <- function(cache = NULL) {
  with(cache, {
    # d logl / d pi
    dldpi <- 1 / pi.i # unweighted!

    # we assume z1/z2 are available
    p1 <- dfun(z1)
    p2 <- dfun(z2)

    # th
    scores.th <- dldpi * (Y1 * p1 - Y2 * p2)

    # beta
    if (nexo > 0L) {
      scores.beta <- dldpi * (-X) * (p1 - p2)
      return(cbind(scores.th, scores.beta, deparse.level = 0))
    } else {
      return(scores.th)
    }
  })
}

lav_uvord_gradient_cache <- function(cache = NULL) {
  with(cache, {
    # d logl / d pi
    wtp <- wt / pi.i
    p1 <- dfun(z1)
    p2 <- dfun(z2)

    # th
    dxa <- Y1 * p1 - Y2 * p2
    scores.th <- wtp * dxa

    # beta
    if (nexo > 0L) {
      dxb <- X * (p1 - p2) # == X*p1 - X*p2
      scores.beta <- wtp * (-dxb)
      return(colSums(cbind(scores.th, scores.beta, deparse.level = 0),
        na.rm = TRUE
      ))
    } else {
      return(colSums(scores.th, na.rm = TRUE))
    }
  })
}


# compute total Hessian
lav_uvord_hessian <- function(y = NULL,
                              X = NULL,
                              wt = rep(1, length(y)),
                              logistic = FALSE,
                              cache = NULL) {
  if (is.null(cache)) {
    cache <- lav_uvord_fit(
      y = y, X = X, wt = wt,
      logistic = logistic, output = "cache"
    )
  }
  tmp <- lav_uvord_loglik_cache(cache = cache)
  tmp <- lav_uvord_gradient_cache(cache = cache)
  lav_uvord_hessian_cache(cache = cache)
}

lav_uvord_hessian_cache <- function(cache = NULL) {
  with(cache, {
    wtp2 <- wt / (pi.i * pi.i)

    g1w <- gfun(z1) * wtp
    g2w <- gfun(z2) * wtp

    Y1gw <- Y1 * g1w
    Y2gw <- Y2 * g2w

    dx2.tau <- (lav_crossprod(Y1gw, Y1) - lav_crossprod(Y2gw, Y2) -
      lav_crossprod(dxa, dxa * wtp2))

    if (nexo == 0L) {
      return(dx2.tau)
    }

    dxb2 <- dxb * wtp2
    dx2.beta <- (lav_crossprod(X * g1w, X) - lav_crossprod(X * g2w, X) -
      lav_crossprod(dxb, dxb2))

    dx.taubeta <- (-lav_crossprod(Y1gw, X) + lav_crossprod(Y2gw, X) +
      lav_crossprod(dxa, dxb2))

    Hessian <- rbind(cbind(dx2.tau, dx.taubeta, deparse.level = 0),
      cbind(t(dx.taubeta), dx2.beta, deparse.level = 0),
      deparse.level = 0
    )

    return(Hessian)
  })
}


# compute total (log)likelihood, for specific 'x' (nlminb)
lav_uvord_min_objective <- function(x, cache = NULL) {
  # check order of first 2 thresholds; if x[1] > x[2], return Inf
  # new in 0.6-8
  if (cache$nth > 1L && x[1] > x[2]) {
    return(+Inf)
  }
  if (cache$nth > 2L && x[2] > x[3]) {
    return(+Inf)
  }
  if (cache$nth > 3L && x[3] > x[4]) {
    return(+Inf)
  }
  cache$theta <- x
  -1 * lav_uvord_loglik_cache(cache = cache) / cache$N
}

# compute gradient, for specific 'x' (nlminb)
lav_uvord_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_uvord_loglik_cache(cache = cache)
  }
  -1 * lav_uvord_gradient_cache(cache = cache) / cache$N
}

# compute hessian, for specific 'x' (nlminb)
lav_uvord_min_hessian <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    tmp <- lav_uvord_loglik_cache(cache = cache)
    tmp <- lav_uvord_gradient_cache(cache = cache)
  }
  -1 * lav_uvord_hessian_cache(cache = cache) / cache$N
}

# get 'z1' and 'z2' values, given (new) values for the parameters
# only needed for lav_bvord_cor_scores(), which is called from
# pml_deriv1() in lav_model_gradient_pml.R
lav_uvord_update_fit <- function(fit.y = NULL, th.new = NULL, sl.new = NULL) {
  # return fit.y with 'update' z1/z2 values
  if (is.null(th.new) && is.null(sl.new)) {
    return(fit.y)
  }

  if (!is.null(th.new)) {
    fit.y$theta[fit.y$th.idx] <- th.new
  }
  if (!is.null(sl.new)) {
    fit.y$theta[fit.y$slope.idx] <- sl.new
  }

  nth <- length(fit.y$th.idx)
  o1 <- ifelse(fit.y$y == nth + 1, 100, 0)
  o2 <- ifelse(fit.y$y == 1, -100, 0)

  theta <- fit.y$theta
  th <- theta[1:nth]
  TH <- c(0, th, 0)
  beta <- theta[-c(1:nth)]
  y <- fit.y$y
  X <- fit.y$X
  if (length(fit.y$slope.idx) > 0L) {
    eta <- drop(X %*% beta)
    fit.y$z1 <- TH[y + 1L] - eta + o1
    fit.y$z2 <- TH[y] - eta + o2
  } else {
    fit.y$z1 <- TH[y + 1L] + o1
    fit.y$z2 <- TH[y] + o2
  }

  fit.y
}
