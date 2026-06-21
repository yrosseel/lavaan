# Factor extraction method(s)
# YR Feb 2020
#
# - ULS_corner only (for now)
# - just to get better starting values for ESEM

# YR July 2020
# - adding generic function lav_efa_extraction, using eigenvalue based
#   approach; ML and ULS
# - 'corner' is an option

lav_efa_extraction <- function(s, nfactors = 1L,
                               method = "ULS", # or ML
                               corner = FALSE,
                               reflect = FALSE, order_lv_by = "none",
                               min_var = 0.0001) {
  stopifnot(is.matrix(s))
  s <- unname(s)
  method <- tolower(method)

  # extract variances
  s_var <- diag(s)

  # force S to be pd (eg if we have polychoric correlations)
  s <- lav_mat_sym_force_pd(s, tol = 1e-08)

  # convert to correlation matrix (ULS is not scale invariant!)
  r <- cov2cor(s)

  # optim.method
  if (method == "uls") {
    min_objective <- lav_efa_uls_min_objective
    min_gradient <- lav_efa_uls_min_grad
    cache <- lav_efa_uls_init_cache(r = r, nfactors = nfactors)
  } else if (method == "ml") {
    min_objective <- lav_efa_ml_min_objective
    min_gradient <- lav_efa_ml_min_grad
    cache <- lav_efa_ml_init_cache(r = r, nfactors = nfactors)
  } else {
    lav_msg_stop(gettext("method must be uls or ml (for now)"))
  }
  min_hessian <- NULL

  # optimize
  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = if (lav_verbose()) {
      1L
    } else {
      0L
    },
    abs.tol = (.Machine$double.eps * 10)
  )
  if (lav_verbose()) {
    cat("\n")
    cat("factor extraction iterations using ", method, ":\n")
  }
  out <- nlminb(
    start = cache$theta, objective = min_objective,
    gradient = min_gradient, hessian = min_hessian,
    control = control_nlminb, lower = min_var, upper = +1,
    cache = cache
  )

  # extract LAMBDA/THETA
  if (method == "uls") {
    mm_theta <- diag(out$par * out$par)

    # compute LAMBDA
    a <- r
    diag(a) <- diag(a) - (out$par * out$par)
    ev <- eigen(a, symmetric = TRUE)
    omega_1 <- ev$vectors[, 1:nfactors]
    gamma_1 <- ev$values[1:nfactors]

    # LAMBDA <- Omega.1 %*% diag(sqrt(gamma.1))
    mm_lambda <- t(t(omega_1) * sqrt(gamma_1))

    # rescale if the input matrix was not a correlation matrix
    mm_lambda <- sqrt(s_var) * mm_lambda
    diag(mm_theta) <- s_var * diag(mm_theta)
  } else if (method == "ml") {
    mm_theta <- diag(out$par * out$par)

    # compute LAMBDA
    psi <- out$par
    a <- t(psi * cache$r_inv) * psi
    ev <- eigen(a, symmetric = TRUE)
    omega_1 <- ev$vectors[, 1L + cache$nvar - seq_len(cache$nfactors),
      drop = FALSE
    ]
    gamma_1 <- ev$values[1L + cache$nvar - seq_len(cache$nfactors)]

    # LAMBDA <- diag(psi) %*% Omega.1 %*%sqrt(solve(Gamma.1)-diag(nfactors))
    tmp1 <- psi * omega_1
    mm_lambda <- t(t(tmp1) * sqrt((1 / gamma_1) - 1))

    # rescale if the input matrix was not a correlation matrix
    mm_lambda <- sqrt(s_var) * mm_lambda
    diag(mm_theta) <- s_var * diag(mm_theta)
  }

  # corner?
  if (corner) {
    # rotate to echelon pattern (see echelon() in GPArotation package)
    head_1 <- mm_lambda[seq_len(nfactors), , drop = FALSE]
    post <- try(solve(head_1, t(chol(tcrossprod(head_1)))), silent = TRUE)
    okflag <- FALSE
    if (inherits(post, "try-error")) { # new in 0.6-18
      # this will happen if we have identical elements in the columns
      # of HEAD (perhaps the data is artificial?)
      # -> add some fuzz and try again
      sd_1 <- sqrt(mean(abs(head_1))) * 1e-04
      fuzz <- matrix(rnorm(nfactors * nfactors, 0, sd_1), nfactors, nfactors)
      head2 <- head_1 + fuzz
      post <- try(solve(head2, t(chol(tcrossprod(head2)))), silent = TRUE)
      if (!inherits(post, "try-error")) {
        okflag <- TRUE
      }
    } else {
      okflag <- TRUE
    }
    if (okflag) {
      mm_lambda <- mm_lambda %*% post
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
  neg_idx <- which(diag(mm_lambda) < 0)
  if (length(neg_idx) > 0L) {
    mm_lambda[, neg_idx] <- -1 * mm_lambda[, neg_idx, drop = FALSE]
  }

  # reflect so that column sum is always positive
  if (reflect) {
    sum_1 <- colSums(mm_lambda)
    neg_idx <- which(sum_1 < 0)
    if (length(neg_idx) > 0L) {
      mm_lambda[, neg_idx] <- -1 * mm_lambda[, neg_idx, drop = FALSE]
    }
  }

  # reorder the columns
  if (order_lv_by == "sumofsquares") {
    l2 <- mm_lambda * mm_lambda
    order_idx <- base::order(colSums(l2), decreasing = TRUE)
  } else if (order_lv_by == "index") {
    # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
    max_loading <- apply(abs(mm_lambda), 2, max)
    # 1: per factor, number of the loadings that are at least 0.8 of the
    #    highest loading of the factor
    # 2: mean of the index numbers
    average_index <- sapply(seq_len(ncol(mm_lambda)), function(i) {
      mean(which(abs(mm_lambda[, i]) >= 0.8 * max_loading[i]))
    })
    # order of the factors
    order_idx <- base::order(average_index)
  } else if (order_lv_by == "none") {
    order_idx <- seq_len(ncol(mm_lambda))
  } else {
    lav_msg_stop(gettext("order must be index, sumofsquares or none"))
  }
  mm_lambda <- mm_lambda[, order_idx, drop = FALSE]

  list(LAMBDA = mm_lambda, THETA = mm_theta)
}



lav_efa_uls_init_cache <- function(r = NULL,
                                          nfactors = 1L,
                                          parent = parent.frame()) {
  r_inv <- solve(r)
  nvar <- ncol(r)

  # starting values for diagonal elements of THETA
  # using Joreskog (1966) suggestion:
  theta_init <- (1 - nfactors / (2 * nvar)) * 1 / diag(r_inv)
  theta <- sqrt(theta_init)

  out <- list2env(
    list(
      r = r, nfactors = nfactors,
      theta = theta
    ),
    parent = parent
  )
  out
}

# x is here the sqrt() of theta!
lav_efa_uls_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  with(cache, {    # nolint start
    a <- r
    diag(a) <- diag(a) - (theta * theta)
    ev <- eigen(a, symmetric = TRUE, only.values = TRUE)
    gamma_2 <- ev$values[-seq_len(nfactors)]

    res <- 0.5 * sum(gamma_2 * gamma_2)
    return(res)
  })               # nolint end
}

lav_efa_uls_min_grad <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    # nothing to do
  }
  with(cache, {       # nolint start
    a <- r
    diag(a) <- diag(a) - (theta * theta)
    ev <- eigen(a, symmetric = TRUE)
    omega_2 <- ev$vectors[, -seq_len(nfactors)]
    gamma_2 <- ev$values[-seq_len(nfactors)]

    res <- -2 * theta * colSums(t(omega_2 * omega_2) * gamma_2)
    return(res)
  })                  # nolint end
}


# ML
lav_efa_ml_init_cache <- function(r = NULL,
                                  nfactors = 1L,
                                  parent = parent.frame()) {
  r_inv <- solve(r)
  nvar <- ncol(r)

  # starting values for diagonal elements of THETA
  # using Joreskog (1966) suggestion:
  theta_init <- (1 - nfactors / (2 * nvar)) * 1 / diag(r_inv)
  theta <- sqrt(theta_init)

  out <- list2env(
    list(
      r = r, nfactors = nfactors, r_inv = r_inv,
      nvar = nvar, # for ML only
      theta = theta
    ),
    parent = parent
  )
  out
}


# x is here the sqrt of theta
lav_efa_ml_min_objective <- function(x, cache = NULL) {
  cache$theta <- x
  with(cache, {       # nolint start
    psi <- theta
    # A <- diag(psi) %*% R.inv %*% diag(psi)
    a <- t(r_inv * psi) * psi
    ev <- eigen(a, symmetric = TRUE, only.values = TRUE)
    gamma_2 <- ev$values[(nvar - nfactors):1L]

    res <- sum(log(gamma_2) + 1 / gamma_2 - 1)
    return(res)
  })                  # nolint end
}

lav_efa_ml_min_grad <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    # nothing to do
  }
  with(cache, {          # nolint start
    psi <- theta
    # A <- diag(psi) %*% solve(S) %*% diag(psi)
    a <- t(r_inv * psi) * psi
    ev <- eigen(a, symmetric = TRUE)

    omega_2 <- ev$vectors[, (nvar - nfactors):1L, drop = FALSE]
    gamma_2 <- ev$values[(nvar - nfactors):1L]

    res <- colSums(t(omega_2 * omega_2) * (1 - 1 / gamma_2))
    return(res)
  })                     # nolint end
}




# ULS estimation
#
# - but resulting in a upper-corner all zeroes LAMBDA matrix
# - not using eigenvalues/vectors, but minimizing the residuals
#   directly

# - should give the same results as MINRES (after an orthogonal transformation)
# - unless there are heywood cases; this function allows for negative variances!

lav_efa_extraction_uls_corner <- function(s, nfactors = 1L, reflect = TRUE,
                                          order_lv_by = "none") {
  stopifnot(is.matrix(s))
  s <- unname(s)
  nvar <- nrow(s)

  # extract variances
  s_var <- diag(s)

  # convert to correlation matrix (ULS is not scale invariant!)
  r <- cov2cor(s)
  # R.inv <- solve(R)

  # eigenvalue decomposition (to get starting values for LAMBDA)
  ev <- eigen(r, symmetric = TRUE)

  # extract first nfac components (assuming no measurement error)
  pc <- (ev$vectors[, seq_len(nfactors), drop = FALSE] %*%
         diag(sqrt(ev$values[seq_len(nfactors)])))

  # rotate to echelon pattern (see echelon() in GPArotation package)
  head_1 <- pc[seq_len(nfactors), , drop = FALSE]
  mm_lambda <- pc %*% solve(head_1, t(chol(tcrossprod(head_1))))

  mm_theta <- diag(nvar)
  if (nfactors > 1L) {
    corner_idx <- which(row(mm_lambda) < nfactors &
                        col(mm_lambda) > row(mm_lambda))
    lambda_idx <- seq_len(nvar * nfactors)[-corner_idx]
    mm_lambda[corner_idx] <- 0 # to make them exactly zero
  } else {
    corner_idx <- integer(0L)
    lambda_idx <- seq_len(nvar)
  }

  # optim.method
  min_objective <- lav_efa_uls_corner_min_objective
  min_gradient <- lav_efa_uls_corner_min_grad
  min_hessian <- NULL

  # create cache environment
  cache <- lav_efa_uls_corner_init_cache(
    mm_lambda = mm_lambda,
    lambda_idx = lambda_idx, r = r
  )

  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L,
    trace = if (lav_verbose()) {
      1L
    } else {
      0L
    },
    abs.tol = (.Machine$double.eps * 10)
  )

  # optimize
  out <- nlminb(
    start = cache$theta, objective = min_objective,
    gradient = min_gradient, hessian = min_hessian,
    control = control_nlminb, lower = -1, upper = +1,
    cache = cache
  )

  mm_lambda[lambda_idx] <- out$par
  diag(mm_theta) <- 1 - diag(tcrossprod(mm_lambda))

  # rescale if the input matrix was not a correlation matrix
  mm_lambda <- sqrt(s_var) * mm_lambda
  diag(mm_theta) <- s_var * diag(mm_theta)

  # reflect so that column sum is always positive
  if (reflect) {
    sum_1 <- colSums(mm_lambda)
    neg_idx <- which(sum_1 < 0)
    if (length(neg_idx) > 0L) {
      mm_lambda[, neg_idx] <- -1 * mm_lambda[, neg_idx, drop = FALSE]
    }
  }

  # reorder the columns
  if (order_lv_by == "sumofsquares") {
    l2 <- mm_lambda * mm_lambda
    order_idx <- base::order(colSums(l2), decreasing = TRUE)
  } else if (order_lv_by == "index") {
    # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
    max_loading <- apply(abs(mm_lambda), 2, max)
    # 1: per factor, number of the loadings that are at least 0.8 of the
    #    highest loading of the factor
    # 2: mean of the index numbers
    average_index <- sapply(seq_len(ncol(mm_lambda)), function(i) {
      mean(which(abs(mm_lambda[, i]) >= 0.8 * max_loading[i]))
    })
    # order of the factors
    order_idx <- base::order(average_index)
  } else if (order_lv_by == "none") {
    order_idx <- seq_len(ncol(mm_lambda))
  } else {
    lav_msg_stop(gettext("order must be index, sumofsquares or none"))
  }
  mm_lambda <- mm_lambda[, order_idx, drop = FALSE]

  list(LAMBDA = mm_lambda, THETA = mm_theta)
}


lav_efa_uls_corner_init_cache <- function(mm_lambda = NULL,
                                          lambda_idx = NULL,
                                          r = NULL,
                                          parent = parent.frame()) {
  theta <- mm_lambda[lambda_idx]
  out <- list2env(
    list(
      mm_lambda = mm_lambda,
      lambda_idx = lambda_idx,
      r = r,
      theta = theta
    ),
    parent = parent
  )
  out
}

lav_efa_uls_corner_min_objective <- function(x, cache = NULL) { # nolint
  cache$theta <- x
  with(cache, {               # nolint start
    mm_lambda[lambda_idx] <- theta
    res1 <- lav_mat_vech(r - tcrossprod(mm_lambda), diagonal = FALSE)
    res2 <- res1 * res1
    return(sum(res2))
  })                          # nolint end
}

lav_efa_uls_corner_min_grad <- function(x, cache = NULL) {
  # check if x has changed
  if (!all(x == cache$theta)) {
    cache$theta <- x
    # nothing to do
  }
  with(cache, {                            # nolint start
    mm_lambda[lambda_idx] <- theta
    sigma_1 <- tcrossprod(mm_lambda)
    diag(sigma_1) <- 1 # diagonal is ignored
    tmp <- -2 * (r - sigma_1) %*% mm_lambda
    return(tmp[lambda_idx])
  })                                       # nolint end
}
