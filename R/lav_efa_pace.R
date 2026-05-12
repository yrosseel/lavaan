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
lav_efa_pace <- function(s, nfactors = 1L, p_idx = seq_len(ncol(s)),
                         reflect = TRUE, order_lv_by = "none",
                         use_r = TRUE, theta_only = TRUE) {
  s <- unname(s)
  nvar <- ncol(s)
  theta <- numeric(nvar)
  stopifnot(nfactors < nvar / 2)

  # because subset selection is not scale-invariant, we transform
  # S to R, compute theta based on R, and then rescale again
  if (use_r) {
    s_var <- diag(s)
    r <- stats::cov2cor(s)
  } else {
    r <- s
  }

  # find principal variables ('largest' sub-block)
  a <- r

  # row indices
  v_r <- integer(0L)
  # column indices
  v_c <- integer(0L)

  for (h in seq_len(nfactors)) {
    # mask
    mask_idx <- c(v_r, v_c)
    tmp <- abs(a)
    if (length(mask_idx) > 0L) {
      tmp[mask_idx, ] <- 0
      tmp[, mask_idx] <- 0
    }
    diag(tmp) <- 0

    # find maximum off-diagonal element
    idx <- which(tmp == max(tmp), arr.ind = TRUE, useNames = FALSE)[1, ]
    k <- idx[1]
    l <- idx[2]
    v_r <- c(v_r, k)
    v_c <- c(v_c, l)

    # non-symmetric sweep operator
    a_kl <- a[k, l]
    if (abs(a_kl) < sqrt(.Machine$double.eps)) {
      out <- a
      out[k, ] <- 0
      out[, l] <- 0
    } else {
      out <- a - tcrossprod(a[, l], a[k, ]) / a_kl
      out[k, ] <- a[k, ] / a_kl
      out[, l] <- -a[, l] / a_kl
      out[k, l] <- 1 / a_kl
    }
    a <- out
  }
  # diagonal elements are estimates of theta
  # for all variables not in (v.r, v.c)
  all_idx <- seq_len(nvar)
  # v_r_init <- v_r
  # v_c_init <- v_c
  other_idx <- all_idx[-c(v_r, v_c)]
  theta[other_idx] <- diag(a)[other_idx]


  # now fill in theta for the 2*m remaining variables in c(v.r.init, v.c.init)
  for (i in p_idx) {
    if (i %in% other_idx) {
      next
    }

    # row indices
    v_r <- integer(0L)
    # column indices
    v_c <- integer(0L)

    a <- r
    for (h in seq_len(nfactors)) {
      # mask
      mask_idx <- c(i, v_r, v_c)
      tmp <- abs(a)
      tmp[mask_idx, ] <- 0
      tmp[, mask_idx] <- 0
      diag(tmp) <- 0

      # find maximum off-diagonal element
      idx <- which(tmp == max(tmp), arr.ind = TRUE, useNames = FALSE)[1, ]
      k <- idx[1]
      l <- idx[2]
      v_r <- c(v_r, k)
      v_c <- c(v_c, l)

      # non-symmetric sweep operator
      a_kl <- a[k, l]
      if (abs(a_kl) < sqrt(.Machine$double.eps)) {
        out <- a
        out[k, ] <- 0
        out[, l] <- 0
      } else {
        out <- a - tcrossprod(a[, l], a[k, ]) / a_kl
        out[k, ] <- a[k, ] / a_kl
        out[, l] <- -a[, l] / a_kl
        out[k, l] <- 1 / a_kl
      }
      a <- out
    }

    # diagonal element is estimate of theta
    theta[i] <- a[i, i]
  }

  # return theta elements only
  if (theta_only) {
    # rescale back to S metric
    if (use_r) {
      theta <- theta * s_var
    }
    return(theta[p_idx])
  }

  # compute LAMBDA using the 'eigenvalue' method
  ev <- eigen(r, symmetric = TRUE)
  s_sqrt <- ev$vectors %*% sqrt(diag(ev$values)) %*% t(ev$vectors)
  s_inv_sqrt <- ev$vectors %*% sqrt(diag(1 / ev$values)) %*% t(ev$vectors)
  rtr <- s_inv_sqrt %*% diag(theta) %*% s_inv_sqrt

  ev <- eigen(rtr, symmetric = TRUE)
  omega_m <- ev$vectors[, 1L + nvar - seq_len(nfactors), drop = FALSE]
  gamma_m <- ev$values[1L + nvar - seq_len(nfactors)]
  gamma_m_1 <- diag(gamma_m, nrow = nfactors, ncol = nfactors)

  # Cudeck 1991 page 37 bottom of the page:
  lambda_dot <- s_sqrt %*% omega_m %*% sqrt(diag(nfactors) - gamma_m_1)

  if (use_r) {
    # IF (and only if) the input is a correlation matrix,
    # we must rescale so that the diag(R.implied) == 1

    # R.unscaled <- tcrossprod(LAMBDA.dot) + diag(theta)
    # r.var.inv <- 1/diag(R.unscaled)

    # LAMBDA/THETA in correlation metric
    # LAMBDA.R <- sqrt(r.var.inv) * LAMBDA.dot
    # THETA.R <- diag(r.var.inv * theta)

    # convert to 'S' metric
    mm_lambda <- sqrt(s_var) * lambda_dot
    mm_theta <- diag(s_var * theta)
  } else {
    mm_lambda <- lambda_dot
    mm_theta <- diag(theta)
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
