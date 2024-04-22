# Factor extraction method(s)
# YR Feb 2020
#
# - ULS_corner only (for now)
# - just to get better starting values for ESEM

# YR July 2020
# - adding generic function lav_efa_extraction, using eigenvalue based
#   approach; ML and ULS
# - 'corner' is an option

lav_efa_extraction <- function(S, nfactors = 1L,
                               method = "ULS", # or ML
                               corner = FALSE,
                               reflect = FALSE, order.lv.by = "none",
                               verbose = FALSE,
                               min.var = 0.0001) {
  stopifnot(is.matrix(S))
  S <- unname(S)
  method <- tolower(method)

  # extract variances
  S.var <- diag(S)

  # force S to be pd (eg if we have polychoric correlations)
  S <- lav_matrix_symmetric_force_pd(S, tol = 1e-08)

  # convert to correlation matrix (ULS is not scale invariant!)
  R <- cov2cor(S)

  # optim.method
  if (method == "uls") {
    minObjective <- efa_extraction_uls_min_objective
    minGradient <- efa_extraction_uls_min_gradient
    cache <- efa_extraction_uls_init_cache(R = R, nfactors = nfactors)
  } else if (method == "ml") {
    minObjective <- efa_extraction_ml_min_objective
    minGradient <- efa_extraction_ml_min_gradient
    cache <- efa_extraction_ml_init_cache(R = R, nfactors = nfactors)
  } else {
    lav_msg_stop(gettext("method must be uls or ml (for now)"))
  }
  minHessian <- NULL

  # optimize
  control.nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = if (verbose) {
      1L
    } else {
      0L
    },
    abs.tol = (.Machine$double.eps * 10)
  )
  out <- nlminb(
    start = cache$theta, objective = minObjective,
    gradient = minGradient, hessian = minHessian,
    control = control.nlminb, lower = min.var, upper = +1,
    cache = cache
  )

  # extract LAMBDA/THETA
  if (method == "uls") {
    THETA <- diag(out$par * out$par)

    # compute LAMBDA
    A <- R
    diag(A) <- diag(A) - (out$par * out$par)
    EV <- eigen(A, symmetric = TRUE)
    Omega.1 <- EV$vectors[, 1:nfactors]
    gamma.1 <- EV$values[1:nfactors]

    # LAMBDA <- Omega.1 %*% diag(sqrt(gamma.1))
    LAMBDA <- t(t(Omega.1) * sqrt(gamma.1))

    # rescale if the input matrix was not a correlation matrix
    LAMBDA <- sqrt(S.var) * LAMBDA
    diag(THETA) <- S.var * diag(THETA)
  } else if (method == "ml") {
    THETA <- diag(out$par * out$par)

    # compute LAMBDA
    psi <- out$par
    A <- t(psi * cache$R.inv) * psi
    EV <- eigen(A, symmetric = TRUE)
    Omega.1 <- EV$vectors[, 1L + cache$nvar - seq_len(cache$nfactors),
      drop = FALSE
    ]
    gamma.1 <- EV$values[1L + cache$nvar - seq_len(cache$nfactors)]

    # LAMBDA <- diag(psi) %*% Omega.1 %*%sqrt(solve(Gamma.1)-diag(nfactors))
    tmp1 <- psi * Omega.1
    LAMBDA <- t(t(tmp1) * sqrt((1 / gamma.1) - 1))

    # rescale if the input matrix was not a correlation matrix
    LAMBDA <- sqrt(S.var) * LAMBDA
    diag(THETA) <- S.var * diag(THETA)
  }

  # corner?
  if (corner) {
    # rotate to echelon pattern (see echelon() in GPArotation package)
    HEAD <- LAMBDA[seq_len(nfactors), , drop = FALSE]
    POST <- try(solve(HEAD, t(chol(tcrossprod(HEAD)))), silent = TRUE)
    okflag <- FALSE
    if (inherits(POST, "try-error")) { # new in 0.6-18
      # this will happen if we have identical elements in the columns
      # of HEAD (perhaps the data is artificial?)
      # -> add some fuzz and try again
      SD <- sqrt(mean(abs(HEAD))) * 1e-04
      fuzz <- matrix(rnorm(nfactors * nfactors, 0, SD), nfactors, nfactors)
      HEAD2 <- HEAD + fuzz
      POST <- try(solve(HEAD2, t(chol(tcrossprod(HEAD2)))), silent = TRUE)
      if (!inherits(POST, "try-error")) {
        okflag <- TRUE
      }
    } else {
      okflag <- TRUE
    }
    if (okflag) {
      LAMBDA <- LAMBDA %*% POST
    } else {
      lav_msg_warn(gettext(
        "rotation of initial factor solution to echelon pattern failed."))
    }
  }

  # ALWAYS change the sign so that largest element in the column is positive
  # neg.max <- apply(LAMBDA, 2, function(x) { sign(x[which.max(abs(x))]) })
  # neg.idx <- which(neg.max < 0)
  # if(length(neg.idx) > 0L) {
  #    LAMBDA[, neg.idx] <- -1 * LAMBDA[, neg.idx, drop = FALSE]
  # }

  # ALWAYS change the sign so that diag(LAMBDA) is positive
  neg.idx <- which(diag(LAMBDA) < 0)
  if (length(neg.idx) > 0L) {
    LAMBDA[, neg.idx] <- -1 * LAMBDA[, neg.idx, drop = FALSE]
  }

  # reflect so that column sum is always positive
  if (reflect) {
    SUM <- colSums(LAMBDA)
    neg.idx <- which(SUM < 0)
    if (length(neg.idx) > 0L) {
      LAMBDA[, neg.idx] <- -1 * LAMBDA[, neg.idx, drop = FALSE]
    }
  }

  # reorder the columns
  if (order.lv.by == "sumofsquares") {
    L2 <- LAMBDA * LAMBDA
    order.idx <- base::order(colSums(L2), decreasing = TRUE)
  } else if (order.lv.by == "index") {
    # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
    max.loading <- apply(abs(LAMBDA), 2, max)
    # 1: per factor, number of the loadings that are at least 0.8 of the
    #    highest loading of the factor
    # 2: mean of the index numbers
    average.index <- sapply(seq_len(ncol(LAMBDA)), function(i) {
      mean(which(abs(LAMBDA[, i]) >= 0.8 * max.loading[i]))
    })
    # order of the factors
    order.idx <- base::order(average.index)
  } else if (order.lv.by == "none") {
    order.idx <- seq_len(ncol(LAMBDA))
  } else {
    lav_msg_stop(gettext("order must be index, sumofsquares or none"))
  }
  LAMBDA <- LAMBDA[, order.idx, drop = FALSE]

  list(LAMBDA = LAMBDA, THETA = THETA)
}



efa_extraction_uls_init_cache <- function(R = NULL,
                                          nfactors = 1L,
                                          parent = parent.frame()) {
  R.inv <- solve(R)
  nvar <- ncol(R)

  # starting values for diagonal elements of THETA
  # using Joreskog (1966) suggestion:
  theta.init <- (1 - nfactors / (2 * nvar)) * 1 / diag(R.inv)
  theta <- sqrt(theta.init)

  out <- list2env(
    list(
      R = R, nfactors = nfactors,
      theta = theta
    ),
    parent = parent
  )
  out
}

# x is here the sqrt() of theta!
efa_extraction_uls_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  with(cache, {
    A <- R
    diag(A) <- diag(A) - (theta * theta)
    EV <- eigen(A, symmetric = TRUE, only.values = TRUE)
    gamma.2 <- EV$values[-seq_len(nfactors)]

    res <- 0.5 * sum(gamma.2 * gamma.2)
    return(res)
  })
}

efa_extraction_uls_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    # nothing to do
  }
  with(cache, {
    A <- R
    diag(A) <- diag(A) - (theta * theta)
    EV <- eigen(A, symmetric = TRUE)
    Omega.2 <- EV$vectors[, -seq_len(nfactors)]
    gamma.2 <- EV$values[-seq_len(nfactors)]

    res <- -2 * theta * colSums(t(Omega.2 * Omega.2) * gamma.2)
    return(res)
  })
}


# ML
efa_extraction_ml_init_cache <- function(R = NULL,
                                         nfactors = 1L,
                                         parent = parent.frame()) {
  R.inv <- solve(R)
  nvar <- ncol(R)

  # starting values for diagonal elements of THETA
  # using Joreskog (1966) suggestion:
  theta.init <- (1 - nfactors / (2 * nvar)) * 1 / diag(R.inv)
  theta <- sqrt(theta.init)

  out <- list2env(
    list(
      R = R, nfactors = nfactors, R.inv = R.inv,
      nvar = nvar, # for ML only
      theta = theta
    ),
    parent = parent
  )
  out
}


# x is here the sqrt of theta
efa_extraction_ml_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  with(cache, {
    psi <- theta
    # A <- diag(psi) %*% R.inv %*% diag(psi)
    A <- t(R.inv * psi) * psi
    EV <- eigen(A, symmetric = TRUE, only.values = TRUE)
    gamma.2 <- EV$values[(nvar - nfactors):1L]

    res <- sum(log(gamma.2) + 1 / gamma.2 - 1)
    return(res)
  })
}

efa_extraction_ml_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    # nothing to do
  }
  with(cache, {
    psi <- theta
    # A <- diag(psi) %*% solve(S) %*% diag(psi)
    A <- t(R.inv * psi) * psi
    EV <- eigen(A, symmetric = TRUE)

    omega.2 <- EV$vectors[, (nvar - nfactors):1L, drop = FALSE]
    gamma.2 <- EV$values[(nvar - nfactors):1L]

    res <- colSums(t(omega.2 * omega.2) * (1 - 1 / gamma.2))
    return(res)
  })
}




# ULS estimation
#
# - but resulting in a upper-corner all zeroes LAMBDA matrix
# - not using eigenvalues/vectors, but minimizing the residuals
#   directly

# - should give the same results as MINRES (after an orthogonal transformation)
# - unless there are heywood cases; this function allows for negative variances!

lav_efa_extraction_uls_corner <- function(S, nfactors = 1L, reflect = TRUE,
                                          order.lv.by = "none",
                                          verbose = TRUE) {
  stopifnot(is.matrix(S))
  S <- unname(S)
  nvar <- nrow(S)

  # extract variances
  S.var <- diag(S)

  # convert to correlation matrix (ULS is not scale invariant!)
  R <- cov2cor(S)
  # R.inv <- solve(R)

  # eigenvalue decomposition (to get starting values for LAMBDA)
  EV <- eigen(R, symmetric = TRUE)

  # extract first nfac components (assuming no measurement error)
  PC <- (EV$vectors[, seq_len(nfactors), drop = FALSE] %*%
    diag(sqrt(EV$values[seq_len(nfactors)])))

  # rotate to echelon pattern (see echelon() in GPArotation package)
  HEAD <- PC[seq_len(nfactors), , drop = FALSE]
  LAMBDA <- PC %*% solve(HEAD, t(chol(tcrossprod(HEAD))))

  THETA <- diag(nvar)
  if (nfactors > 1L) {
    corner.idx <- which(row(LAMBDA) < nfactors & col(LAMBDA) > row(LAMBDA))
    lambda.idx <- seq_len(nvar * nfactors)[-corner.idx]
    LAMBDA[corner.idx] <- 0 # to make them exactly zero
  } else {
    corner.idx <- integer(0L)
    lambda.idx <- seq_len(nvar)
  }

  # optim.method
  minObjective <- efa_extraction_uls_corner_min_objective
  minGradient <- efa_extraction_uls_corner_min_gradient
  minHessian <- NULL

  # create cache environment
  cache <- efa_extraction_uls_corner_init_cache(
    LAMBDA = LAMBDA,
    lambda.idx = lambda.idx, R = R
  )

  control.nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = if (verbose) {
      1L
    } else {
      0L
    },
    abs.tol = (.Machine$double.eps * 10)
  )

  # optimize
  out <- nlminb(
    start = cache$theta, objective = minObjective,
    gradient = minGradient, hessian = minHessian,
    control = control.nlminb, lower = -1, upper = +1,
    cache = cache
  )

  LAMBDA[lambda.idx] <- out$par
  diag(THETA) <- 1 - diag(tcrossprod(LAMBDA))

  # rescale if the input matrix was not a correlation matrix
  LAMBDA <- sqrt(S.var) * LAMBDA
  diag(THETA) <- S.var * diag(THETA)

  # reflect so that column sum is always positive
  if (reflect) {
    SUM <- colSums(LAMBDA)
    neg.idx <- which(SUM < 0)
    if (length(neg.idx) > 0L) {
      LAMBDA[, neg.idx] <- -1 * LAMBDA[, neg.idx, drop = FALSE]
    }
  }

  # reorder the columns
  if (order.lv.by == "sumofsquares") {
    L2 <- LAMBDA * LAMBDA
    order.idx <- base::order(colSums(L2), decreasing = TRUE)
  } else if (order.lv.by == "index") {
    # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
    max.loading <- apply(abs(LAMBDA), 2, max)
    # 1: per factor, number of the loadings that are at least 0.8 of the
    #    highest loading of the factor
    # 2: mean of the index numbers
    average.index <- sapply(seq_len(ncol(LAMBDA)), function(i) {
      mean(which(abs(LAMBDA[, i]) >= 0.8 * max.loading[i]))
    })
    # order of the factors
    order.idx <- base::order(average.index)
  } else if (order.lv.by == "none") {
    order.idx <- seq_len(ncol(LAMBDA))
  } else {
    lav_msg_stop(gettext("order must be index, sumofsquares or none"))
  }
  LAMBDA <- LAMBDA[, order.idx, drop = FALSE]

  list(LAMBDA = LAMBDA, THETA = THETA)
}


efa_extraction_uls_corner_init_cache <- function(LAMBDA = NULL,
                                                 lambda.idx = NULL,
                                                 R = NULL,
                                                 parent = parent.frame()) {
  theta <- LAMBDA[lambda.idx]
  out <- list2env(
    list(
      LAMBDA = LAMBDA,
      lambda.idx = lambda.idx,
      R = R,
      theta = theta
    ),
    parent = parent
  )
  out
}

efa_extraction_uls_corner_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  with(cache, {
    LAMBDA[lambda.idx] <- theta
    res1 <- lav_matrix_vech(R - tcrossprod(LAMBDA), diagonal = FALSE)
    res2 <- res1 * res1
    return(sum(res2))
  })
}

efa_extraction_uls_corner_min_gradient <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    # nothing to do
  }
  with(cache, {
    LAMBDA[lambda.idx] <- theta
    Sigma <- tcrossprod(LAMBDA)
    diag(Sigma) <- 1 # diagonal is ignored
    tmp <- -2 * (R - Sigma) %*% LAMBDA
    return(tmp[lambda.idx])
  })
}
