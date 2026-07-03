# the multivariate normal distribution

# 1) loglikelihood (from raw data, or sample statistics)
# 2) derivatives with respect to mu, Sigma, vech(Sigma)
# 3) casewise scores with respect to mu, vech(Sigma), mu + vech(Sigma)
# 4) hessian mu + vech(Sigma)
# 5) information h0 mu + vech(Sigma)
#    5a: (unit)    expected information
#    5b: (unit)    observed information
#    5c: (unit) first.order information
# 6) inverted information h0 mu + vech(Sigma)
#    6a: (unit) inverted expected information
#    6b: /
#    6c: /
# 7) ACOV h0 mu + vech(Sigma)
#    7a: 1/N * inverted expected    information
#    7b: 1/N * inverted observed    information
#    7c: 1/N * inverted first-order information
#    7d: sandwich acov

# YR 07 Feb 2016: first version
# YR 24 Mar 2016: added firstorder information, hessian logl
# YR 19 Jan 2017: added lav_mvnorm_inverted_information_expected
# YR 27 Jun 2018: adding cluster.idx= argument for information_firstorder
# YR 04 Okt 2018: adding wt= argument, and missing meanstructure=
# YR 24 Jul 2022: adding correlation= argument for information_expected
#                 (only for catml; not for correlation = TRUE!)
# YR 03 Jul 2026: use shared kernels from lav_mvnorm_kernels.R

# 0. densities
lav_mvn_dmvnorm <- function(y = NULL,
                               wt = NULL,
                               mu = NULL,
                               sigma_1 = NULL,
                               sigma_inv = NULL,
                               sinv_method = "eigen",
                               x_idx = integer(0L),
                               x_mean = NULL,
                               x_cov = NULL,
                               log = TRUE) {
  if (is.matrix(y)) {
    if (is.null(mu) && is.null(sigma_1) && is.null(sigma_inv)) {
      out <- lav_mvn_loglik_data_z(y = y, casewise = TRUE)
    } else {
      out <- lav_mvn_loglik_data(
        y = y, mu = mu, sigma_1 = sigma_1,
        casewise = TRUE,
        sinv_method = sinv_method
      )
    }
  } else {
    # just one
    p <- length(y)
    log_2pi <- log(2 * pi)

    if (is.null(mu) && is.null(sigma_1) && is.null(sigma_inv)) {
      # mahalanobis distance
      dist_1 <- sum(y * y)
      out <- -(p * log_2pi + dist_1) / 2
    } else {
      sigma_inv <- lav_mvn_sigma_inv(
        sigma_1 = sigma_1, sigma_inv = sigma_inv,
        sinv_method = sinv_method, logdet = TRUE
      )
      logdet <- attr(sigma_inv, "logdet")

      # mahalanobis distance
      yc <- y - mu
      dist_1 <- sum(yc %*% sigma_inv * yc)
      out <- -(p * log_2pi + logdet + dist_1) / 2
    }
  }

  if (!is.null(wt)) {
    out <- out * wt
  }

  # x.idx?
  if (length(x_idx) > 0L) {
    if (is.null(sigma_1) && is.null(x_cov)) {
      lav_msg_stop(gettext("when x.idx is not empty, we need Sigma or x.cov"))
    }
    if (is.matrix(y)) {
      x <- y[, x_idx, drop = FALSE]
    } else {
      x <- y[x_idx]
    }

    mu_x <- x_mean
    sigma_x <- x_cov
    if (is.null(x_mean)) {
      mu_x <- as.numeric(mu)[x_idx]
    }
    if (is.null(x_cov)) {
      sigma_x <- sigma_1[x_idx, x_idx, drop = FALSE]
    }

    logl_x <- lav_mvn_dmvnorm(
      y = x, wt = wt, mu = mu_x, sigma_1 = sigma_x,
      sigma_inv = NULL,
      sinv_method = sinv_method,
      x_idx = integer(0L), log = TRUE
    )

    # subtract logl.X
    out <- out - logl_x
  }

  if (!log) {
    out <- exp(out)
  }

  out
}

# 1. likelihood

# 1a: input is raw data
# (note casewise = TRUE same as: dmvnorm(Y, mean, sigma, log = TRUE))
lav_mvn_loglik_data <- function(y = NULL,
                                   wt = NULL,
                                   mu = NULL,
                                   sigma_1 = NULL,
                                   x_idx = integer(0L),
                                   x_mean = NULL,
                                   x_cov = NULL,
                                   casewise = FALSE,
                                   sinv_method = "eigen") {
  # Y must be a matrix (use lav_mvn_dmvnorm() for non-matrix input)
  stopifnot(is.matrix(y))

  n <- lav_mvn_nobs(y = y, wt = wt)
  p <- NCOL(y)
  mu <- as.numeric(mu)

  if (casewise) {
    log_2pi <- log(2 * pi)

    # invert Sigma
    if (sinv_method == "chol") {
      c_s <- chol(sigma_1)
      ic_s <- backsolve(c_s, diag(p))
      yc <- t(t(y) - mu)
      dist_1 <- rowSums((yc %*% ic_s)^2)
      logdet <- -2 * sum(log(diag(ic_s)))
    } else {
      sigma_inv <- lav_mat_sym_inverse(
        s = sigma_1, logdet = TRUE,
        sinv_method = sinv_method
      )
      logdet <- attr(sigma_inv, "logdet")
      # mahalanobis distance
      yc <- t(t(y) - mu)
      dist_1 <- rowSums(yc %*% sigma_inv * yc)
    }

    loglik <- -(p * log_2pi + logdet + dist_1) / 2

    # weights
    if (!is.null(wt)) {
      loglik <- loglik * wt
    }
  } else {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = TRUE,
      sinv_method = sinv_method
    )
    samp <- lav_mvn_samp_stats(y = y, wt = wt)
    loglik <- lav_mvn_loglik_samp(
      sample_mean = samp$mean,
      sample_cov = samp$cov,
      sample_nobs = n,
      mu = mu,
      sigma_inv = sigma_inv
    )
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    mu_x <- x_mean
    sigma_x <- x_cov
    if (is.null(x_mean)) {
      mu_x <- as.numeric(mu)[x_idx]
    }
    if (is.null(x_cov)) {
      sigma_x <- sigma_1[x_idx, x_idx, drop = FALSE]
    }
    loglik_x <- lav_mvn_loglik_data(
      y = y[, x_idx, drop = FALSE],
      wt = wt, mu = mu_x, sigma_1 = sigma_x,
      x_idx = integer(0L), casewise = casewise,
      sinv_method = sinv_method
    )
    # subtract logl.X
    loglik <- loglik - loglik_x
  }

  loglik
}



# 1b: input are sample statistics (mean, cov, N) only
lav_mvn_loglik_samp <- function(sample_mean = NULL,
                                          sample_cov = NULL,
                                          sample_nobs = NULL,
                                          mu = NULL,
                                          sigma_1 = NULL,
                                          x_idx = integer(0L),
                                          x_mean = NULL,
                                          x_cov = NULL,
                                          sinv_method = "eigen",
                                          sigma_inv = NULL) {
  p <- length(sample_mean)
  n <- sample_nobs
  mu <- as.numeric(mu)
  sample_mean <- as.numeric(sample_mean)
  log_2pi <- log(2 * pi)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method, logdet = TRUE
  )
  logdet <- attr(sigma_inv, "logdet")

  # tr(Sigma^{-1} %*% S)
  dist1 <- sum(sigma_inv * sample_cov)
  # (ybar - mu)^T %*% Sigma.inv %*% (ybar - mu)
  diff_1 <- as.numeric(sample_mean - mu)
  dist2 <- sum(as.numeric(crossprod(diff_1, sigma_inv)) * diff_1)

  loglik <- -n / 2 * (p * log_2pi + logdet + dist1 + dist2)

  # fixed.x?
  if (length(x_idx) > 0L) {
    mu_x <- x_mean
    sigma_x <- x_cov
    if (is.null(x_mean)) {
      mu_x <- mu[x_idx]
    }
    if (is.null(x_cov)) {
      sigma_x <- sigma_1[x_idx, x_idx, drop = FALSE]
    }
    sample_mean_x <- sample_mean[x_idx]
    sample_cov_x <- sample_cov[x_idx, x_idx, drop = FALSE]
    loglik_x <-
      lav_mvn_loglik_samp(
        sample_mean = sample_mean_x,
        sample_cov = sample_cov_x,
        sample_nobs = sample_nobs,
        mu = mu_x, sigma_1 = sigma_x,
        x_idx = integer(0L),
        sinv_method = sinv_method
      )
    # subtract logl.X
    loglik <- loglik - loglik_x
  }

  loglik
}

# 1c special case: Mu = 0, Sigma = I
lav_mvn_loglik_data_z <- function(y = NULL,
                                     wt = NULL,
                                     casewise = FALSE) {
  n <- lav_mvn_nobs(y = y, wt = wt)
  p <- NCOL(y)
  log_2pi <- log(2 * pi)

  if (casewise) {
    dist_1 <- rowSums(y * y)
    loglik <- -(p * log_2pi + dist_1) / 2
    if (!is.null(wt)) {
      loglik <- loglik * wt
    }
  } else {
    samp <- lav_mvn_samp_stats(y = y, wt = wt)

    dist1 <- sum(diag(samp$cov))
    dist2 <- sum(samp$mean * samp$mean)

    loglik <- -n / 2 * (p * log_2pi + dist1 + dist2)
  }

  loglik
}





# 2. Derivatives

# 2a: derivative logl with respect to mu
lav_mvn_dlogl_dmu <- function(y = NULL,
                                 wt = NULL,
                                 mu = NULL,
                                 sigma_1 = NULL,
                                 x_idx = integer(0L),
                                 sinv_method = "eigen",
                                 sigma_inv = NULL) {
  mu <- as.numeric(mu)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  # subtract 'Mu' from Y
  yc <- t(t(y) - mu)

  # weights
  if (!is.null(wt)) {
    yc <- yc * wt
  }

  # derivative
  dmu <- as.numeric(sigma_inv %*% colSums(yc))

  # fixed.x?
  if (length(x_idx) > 0L) {
    dmu[x_idx] <- 0
  }

  dmu
}

# 2b: derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvn_dlogl_dsigma <- function(y = NULL,
                                    wt = NULL,
                                    mu = NULL,
                                    sigma_1 = NULL,
                                    x_idx = integer(0L),
                                    sinv_method = "eigen",
                                    sigma_inv = NULL) {
  n <- lav_mvn_nobs(y = y, wt = wt)
  mu <- as.numeric(mu)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  # W.tilde
  w_tilde <- lav_mvn_w_tilde(y = y, wt = wt, mu = mu)

  # derivative
  d_sigma <- -(n / 2) * (sigma_inv - (sigma_inv %*% w_tilde %*% sigma_inv))

  # fixed.x?
  if (length(x_idx) > 0L) {
    d_sigma[x_idx, x_idx] <- 0
  }

  d_sigma
}

# 2c: derivative logl with respect to vech(Sigma)
lav_mvn_dlogl_dvechsigma <- function(y = NULL,
                                        wt = NULL,
                                        mu = NULL,
                                        sigma_1 = NULL,
                                        x_idx = integer(0L),
                                        sinv_method = "eigen",
                                        sigma_inv = NULL) {
  d_sigma <- lav_mvn_dlogl_dsigma(
    y = y, wt = wt, mu = mu, sigma_1 = sigma_1,
    x_idx = x_idx, sinv_method = sinv_method, sigma_inv = sigma_inv
  )

  # vech
  as.numeric(lav_mat_dup_pre(as.matrix(lav_mat_vec(d_sigma))))
}

# 2d: : derivative logl with respect to Mu and vech(Sigma)
lav_mvn_dlogl_dmu_dvechsigma <- function(m_y = NULL,
                                            wt = NULL,
                                            m_mu = NULL,
                                            m_sigma = NULL,
                                            x_idx = integer(0L),
                                            sinv_method = "eigen",
                                            sigma_inv = NULL) {
  n <- lav_mvn_nobs(y = m_y, wt = wt)
  m_mu <- as.numeric(m_mu)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = m_sigma, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  # subtract Mu
  m_yc <- t(t(m_y) - m_mu)

  # m_w_tilde
  m_w_tilde <- lav_mvn_w_tilde(y = m_y, wt = wt, mu = m_mu)
  if (!is.null(wt)) {
    dmu <- as.numeric(sigma_inv %*% colSums(m_yc * wt))
  } else {
    dmu <- as.numeric(sigma_inv %*% colSums(m_yc))
  }

  # derivative (avoiding kronecker product)
  m_dsigma <- -(n / 2) * (sigma_inv - (sigma_inv %*% m_w_tilde %*% sigma_inv))

  # fixed.x?
  if (length(x_idx) > 0L) {
    m_dsigma[x_idx, x_idx] <- 0
    dmu[x_idx] <- 0
  }

  # vech
  dvech_sigma <- as.numeric(lav_mat_dup_pre(
    as.matrix(lav_mat_vec(m_dsigma))
  ))

  c(dmu, dvech_sigma)
}

# 3. Casewise scores

# 3a: casewise scores with respect to mu
lav_mvn_sc_mu <- function(y = NULL,
                                 wt = NULL,
                                 mu = NULL,
                                 x_idx = integer(0L),
                                 sigma_1 = NULL,
                                 sinv_method = "eigen",
                                 sigma_inv = NULL) {
  mu <- as.numeric(mu)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  # subtract Mu
  yc <- t(t(y) - mu)

  # postmultiply with Sigma.inv
  sc <- yc %*% sigma_inv

  # weights
  if (!is.null(wt)) {
    sc <- sc * wt
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    sc[, x_idx] <- 0
  }

  sc
}

# 3b: casewise scores with respect to vech(Sigma)
lav_mvn_sc_sigma <- function(y = NULL,
                                         wt = NULL,
                                         mu = NULL,
                                         sigma_1 = NULL,
                                         x_idx = integer(0L),
                                         sinv_method = "eigen",
                                         sigma_inv = NULL) {
  p <- NCOL(y)
  mu <- as.numeric(mu)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  # subtract Mu
  yc <- t(t(y) - mu)

  # score kernel
  sc <- lav_mvn_sc_kernel(yc = yc, sigma_inv = sigma_inv)$sigma

  # fixed.x?
  if (length(x_idx) > 0L) {
    sc[, lav_mat_vech_which_idx(n = p, idx = x_idx)] <- 0
  }

  # weights
  if (!is.null(wt)) {
    sc <- sc * wt
  }

  sc
}

# 3c: casewise scores with respect to mu + vech(Sigma)
lav_mvn_sc_mu_sigma <- function(y = NULL,
                                            wt = NULL,
                                            mu = NULL,
                                            sigma_1 = NULL,
                                            x_idx = integer(0L),
                                            sinv_method = "eigen",
                                            sigma_inv = NULL) {
  p <- NCOL(y)
  mu <- as.numeric(mu)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  # subtract Mu
  yc <- t(t(y) - mu)

  # score kernel
  sc2 <- lav_mvn_sc_kernel(yc = yc, sigma_inv = sigma_inv)
  yc <- sc2$mu
  sc <- sc2$sigma

  # fixed.x?
  if (length(x_idx) > 0L) {
    yc[, x_idx] <- 0
    sc[, lav_mat_vech_which_idx(n = p, idx = x_idx)] <- 0
  }

  out <- cbind(yc, sc)

  # weights
  if (!is.null(wt)) {
    out <- out * wt
  }

  out
}


# 4. hessian of logl

# 4a: hessian logl Mu and vech(Sigma) from raw data
lav_mvn_logl_hessian_data <- function(y = NULL,
                                         wt = NULL,
                                         mu = NULL,
                                         sigma_1 = NULL,
                                         x_idx = integer(0L),
                                         sinv_method = "eigen",
                                         sigma_inv = NULL,
                                         meanstructure = TRUE) {
  n <- lav_mvn_nobs(y = y, wt = wt)

  # observed information
  observed <- lav_mvn_info_observed_data(
    y = y, wt = wt, mu = mu,
    sigma_1 = sigma_1, x_idx = x_idx,
    sinv_method = sinv_method, sigma_inv = sigma_inv,
    meanstructure = meanstructure
  )

  -n * observed
}

# 4b: hessian Mu and vech(Sigma) from samplestats
lav_mvnorm_logl_hessian_samplestats <- # nolint
  function(sample_mean = NULL,
           sample_cov = NULL,
           sample_nobs = NULL,
           mu = NULL,
           sigma_1 = NULL,
           x_idx = integer(0L),
           sinv_method = "eigen",
           sigma_inv = NULL,
           meanstructure = TRUE) {
    n <- sample_nobs

    # observed information
    observed <- lav_mvn_info_observed_samp(
      sample_mean = sample_mean, sample_cov = sample_cov,
      mu = mu, sigma_1 = sigma_1,
      x_idx = x_idx, sinv_method = sinv_method, sigma_inv = sigma_inv,
      meanstructure = meanstructure
    )

    -n * observed
  }

# 5) Information h0

# 5a: unit expected information h0 Mu and vech(Sigma)
lav_mvn_info_expected <- function(y = NULL, # unused!
                                            wt = NULL, # unused!
                                            mu = NULL, # unused!
                                            sigma_1 = NULL,
                                            x_idx = integer(0L),
                                            sinv_method = "eigen",
                                            sigma_inv = NULL,
                                            meanstructure = TRUE,
                      correlation = FALSE) {
  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  i22 <- lav_mvn_kron_dup_half(sigma_inv, correlation = correlation)

  if (meanstructure) {
    out <- lav_mat_bdiag(sigma_inv, i22)
  } else {
    out <- i22
  }

  # fixed.x?
  out <- lav_mvn_zero_x_idx(out,
    n = NCOL(sigma_inv), x_idx = x_idx,
    meanstructure = meanstructure
  )

  out
}

# 5b: unit observed information h0
lav_mvn_info_observed_data <- function(y = NULL,
                                                 wt = NULL,
                                                 mu = NULL,
                                                 sigma_1 = NULL,
                                                 x_idx = integer(0L),
                                                 sinv_method = "eigen",
                                                 sigma_inv = NULL,
                                                 meanstructure = TRUE) {
  # sample statistics
  samp <- lav_mvn_samp_stats(y = y, wt = wt)

  lav_mvn_info_observed_samp(
    sample_mean = samp$mean,
    sample_cov = samp$cov, mu = mu, sigma_1 = sigma_1,
    x_idx = x_idx, sinv_method = sinv_method, sigma_inv = sigma_inv,
    meanstructure = meanstructure
  )
}

# 5b-bis: observed information h0 from sample statistics
lav_mvn_info_observed_samp <- function(
    sample_mean = NULL,
    sample_cov = NULL,
    mu = NULL,
    sigma_1 = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen",
    sigma_inv = NULL,
    meanstructure = TRUE) {

  sample_mean <- as.numeric(sample_mean)
  mu <- as.numeric(mu)

  sigma_inv <- lav_mvn_sigma_inv(
    sigma_1 = sigma_1, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  w_tilde <- lav_mvn_w_tilde_samp(
    sample_mean = sample_mean,
    sample_cov = sample_cov, mu = mu
  )

  info <- lav_mvn_info_obs_kernel(
    sigma_inv = sigma_inv, sigma_1 = sigma_1, w_tilde = w_tilde,
    mean_diff = if (meanstructure) sample_mean - mu else NULL
  )

  if (meanstructure) {
    out <- rbind(
      cbind(info$i11, t(info$i21)),
      cbind(info$i21, info$i22)
    )
  } else {
    out <- info$i22
  }

  # fixed.x?
  out <- lav_mvn_zero_x_idx(out,
    n = NCOL(sigma_inv), x_idx = x_idx,
    meanstructure = meanstructure
  )

  out
}

# 5c: unit first-order information h0
lav_mvn_info_firstorder <- function(y = NULL,
                                              wt = NULL,
                                              cluster_idx = NULL,
                                              mu = NULL,
                                              sigma_1 = NULL,
                                              x_idx = integer(0L),
                                              sinv_method = "eigen",
                                              sigma_inv = NULL,
                                              meanstructure = TRUE) {
  n <- lav_mvn_nobs(y = y, wt = wt)

  if (meanstructure) {
    sc <- lav_mvn_sc_mu_sigma(
      y = y, wt = wt,
      mu = mu, sigma_1 = sigma_1, x_idx = x_idx,
      sinv_method = sinv_method, sigma_inv = sigma_inv
    )
  } else {
    # the caller should use Mu = sample_mean
    sc <- lav_mvn_sc_sigma(
      y = y, wt = wt,
      mu = mu, sigma_1 = sigma_1,
      sinv_method = sinv_method, sigma_inv = sigma_inv
    )
  }

  # handle clustering
  if (!is.null(cluster_idx)) {
    # take the sum within each cluster
    sc <- rowsum(sc, group = cluster_idx, reorder = FALSE, na.rm = TRUE)

    # lower bias if number of clusters is not very high
    # FIXME: reference?
    n_c <- nrow(sc)
    correction_factor <- n_c / (n_c - 1)
    sc <- sc * sqrt(correction_factor)
  }

  # unit information
  out <- crossprod(sc) / n

  out
}


# 6: inverted information h0

# 6a: inverted unit expected information h0 Mu and vech(Sigma)
#
#     Note: this is the same as lav_samp_gamma_nt()
#           but where COV=Sigma and MEAN=Mu
#
lav_mvn_inv_info_expected <- function(y = NULL, # unused!
                                                     wt = NULL, # unused!
                                                     mu = NULL, # unused!
                                                     sigma_1 = NULL,
                                                     x_idx = integer(0L),
                                                     meanstructure = TRUE) {
  if (length(x_idx) > 0L) {
    # cov(Y|X) = A - B m_c^{-1} B'
    # where A = cov(Y), B = cov(Y,X), m_c = cov(X)
    a <- sigma_1[-x_idx, -x_idx, drop = FALSE]
    m_b <- sigma_1[-x_idx,  x_idx, drop = FALSE]
    m_c <- sigma_1[x_idx, x_idx, drop = FALSE]
    ybar_x <- a - m_b %*% solve(m_c, t(m_b))

    # reinsert YbarX in Y+X (residual) covariance matrix
    ybar_x_aug <- matrix(0, nrow = NROW(sigma_1), ncol = NCOL(sigma_1))
    ybar_x_aug[-x_idx, -x_idx] <- ybar_x

    # take difference
    r <- sigma_1 - ybar_x_aug

    ss <- lav_mvn_kron_dup_ginv2(sigma_1)
    rr <- lav_mvn_kron_dup_ginv2(r)
    i22 <- ss - rr

    if (meanstructure) {
      i11 <- ybar_x_aug
      out <- lav_mat_bdiag(i11, i22)
    } else {
      out <- i22
    }
  } else {
    i22 <- lav_mvn_kron_dup_ginv2(sigma_1)
    if (meanstructure) {
      i11 <- sigma_1
      out <- lav_mat_bdiag(i11, i22)
    } else {
      out <- i22
    }
  }

  out
}

# 6b:  inverted unit observed information h0

# one could use the inverse of a partitioned matrix, but that does not
# seem to help much... unless we can find an expression for solve(I22)

# 6c: inverted unit first-order information h0
# /


# 7) ACOV h0 mu + vech(Sigma)
# not implemented, as too trivial

#    7a: 1/N * inverted expected    information

#    7b: 1/N * inverted observed    information

#    7c: 1/N * inverted first-order information

#    7d: sandwich acov
