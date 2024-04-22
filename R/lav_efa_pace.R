# Albert (1944a/b) & Ihara & Kano 1986 method to estimate residual variances
# of indicators using a PArtitioned Covariance matrix Estimator (PACE)
#
# The implementation is based on Cudeck 1991:

# Cudeck, R. (1991). Noniterative factor analysis estimators, with algorithms
# for subset and instrumental variable selection. Journal of Educational
# Statistics, 16(1), 35-52.

# YR -- 14 FEB 2020

# - 'fast' version; only (2*nfactors + 1) iterations are needed
# - scale-invariant (by default)
# - always assuming unit variances for the factors
lav_efa_pace <- function(S, nfactors = 1L, p.idx = seq_len(ncol(S)),
                         reflect = TRUE, order.lv.by = "none",
                         use.R = TRUE, theta.only = TRUE) {
  S <- unname(S)
  nvar <- ncol(S)
  theta <- numeric(nvar)
  stopifnot(nfactors < nvar / 2)

  # because subset selection is not scale-invariant, we transform
  # S to R, compute theta based on R, and then rescale again
  if (use.R) {
    s.var <- diag(S)
    R <- stats::cov2cor(S)
  } else {
    R <- S
  }

  # find principal variables ('largest' sub-block)
  A <- R

  # row indices
  v.r <- integer(0L)
  # column indices
  v.c <- integer(0L)

  for (h in seq_len(nfactors)) {
    # mask
    mask.idx <- c(v.r, v.c)
    tmp <- abs(A)
    if (length(mask.idx) > 0L) {
      tmp[mask.idx, ] <- 0
      tmp[, mask.idx] <- 0
    }
    diag(tmp) <- 0

    # find maximum off-diagonal element
    idx <- which(tmp == max(tmp), arr.ind = TRUE, useNames = FALSE)[1, ]
    k <- idx[1]
    l <- idx[2]
    v.r <- c(v.r, k)
    v.c <- c(v.c, l)

    # non-symmetric sweep operator
    a.kl <- A[k, l]
    if (abs(a.kl) < sqrt(.Machine$double.eps)) {
      out <- A
      out[k, ] <- 0
      out[, l] <- 0
    } else {
      out <- A - tcrossprod(A[, l], A[k, ]) / a.kl
      out[k, ] <- A[k, ] / a.kl
      out[, l] <- -A[, l] / a.kl
      out[k, l] <- 1 / a.kl
    }
    A <- out
  }
  # diagonal elements are estimates of theta
  # for all variables not in (v.r, v.c)
  all.idx <- seq_len(nvar)
  v.r.init <- v.r
  v.c.init <- v.c
  other.idx <- all.idx[-c(v.r, v.c)]
  theta[other.idx] <- diag(A)[other.idx]


  # now fill in theta for the 2*m remaining variables in c(v.r.init, v.c.init)
  for (i in p.idx) {
    if (i %in% other.idx) {
      next
    }

    # row indices
    v.r <- integer(0L)
    # column indices
    v.c <- integer(0L)

    A <- R
    for (h in seq_len(nfactors)) {
      # mask
      mask.idx <- c(i, v.r, v.c)
      tmp <- abs(A)
      tmp[mask.idx, ] <- 0
      tmp[, mask.idx] <- 0
      diag(tmp) <- 0

      # find maximum off-diagonal element
      idx <- which(tmp == max(tmp), arr.ind = TRUE, useNames = FALSE)[1, ]
      k <- idx[1]
      l <- idx[2]
      v.r <- c(v.r, k)
      v.c <- c(v.c, l)

      # non-symmetric sweep operator
      a.kl <- A[k, l]
      if (abs(a.kl) < sqrt(.Machine$double.eps)) {
        out <- A
        out[k, ] <- 0
        out[, l] <- 0
      } else {
        out <- A - tcrossprod(A[, l], A[k, ]) / a.kl
        out[k, ] <- A[k, ] / a.kl
        out[, l] <- -A[, l] / a.kl
        out[k, l] <- 1 / a.kl
      }
      A <- out
    }

    # diagonal element is estimate of theta
    theta[i] <- A[i, i]
  }

  # return theta elements only
  if (theta.only) {
    # rescale back to S metric
    if (use.R) {
      theta <- theta * s.var
    }
    return(theta[p.idx])
  }

  # compute LAMBDA using the 'eigenvalue' method
  EV <- eigen(R, symmetric = TRUE)
  S.sqrt <- EV$vectors %*% sqrt(diag(EV$values)) %*% t(EV$vectors)
  S.inv.sqrt <- EV$vectors %*% sqrt(diag(1 / EV$values)) %*% t(EV$vectors)
  RTR <- S.inv.sqrt %*% diag(theta) %*% S.inv.sqrt

  EV <- eigen(RTR, symmetric = TRUE)
  Omega.m <- EV$vectors[, 1L + nvar - seq_len(nfactors), drop = FALSE]
  gamma.m <- EV$values[1L + nvar - seq_len(nfactors)]
  Gamma.m <- diag(gamma.m, nrow = nfactors, ncol = nfactors)

  # Cuceck 1991 page 37 bottom of the page:
  LAMBDA.dot <- S.sqrt %*% Omega.m %*% sqrt(diag(nfactors) - Gamma.m)

  if (use.R) {
    # IF (and only if) the input is a correlation matrix,
    # we must rescale so that the diag(R.implied) == 1

    # R.unscaled <- tcrossprod(LAMBDA.dot) + diag(theta)
    # r.var.inv <- 1/diag(R.unscaled)

    # LAMBDA/THETA in correlation metric
    # LAMBDA.R <- sqrt(r.var.inv) * LAMBDA.dot
    # THETA.R <- diag(r.var.inv * theta)

    # convert to 'S' metric
    LAMBDA <- sqrt(s.var) * LAMBDA.dot
    THETA <- diag(s.var * theta)
  } else {
    LAMBDA <- LAMBDA.dot
    THETA <- diag(theta)
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
