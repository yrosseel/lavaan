# the multivariate normal distribution, unrestricted (h1)
# - everything is evaluated under the MLEs: Mu = ybar, Sigma = S

# 1) loglikelihood h1 (from raw data, or sample statistics)
# 4) hessian h1 around MLEs
# 5) information h1 (restricted Sigma/mu)
#    5a: (unit)    expected information h1 (A1 = gamma_nt^{-1})
#    5b: (unit)    observed information h1 (A1 = gamma_nt^{-1})
#    5c: (unit) first.order information h1 (B1 = A1 %*% Gamma %*% A1)
# 6) inverted information h1 mu + vech(Sigma)
#    6a: (unit) inverted expected information (A1.inv = gamma_nt)
#    6b: (unit) inverted observed information (A1.inv = gamma_nt)
#    6c: (unit) inverted first-order information (B1.inv)
# 7) ACOV h1 mu + vech(Sigma)
#    7a: 1/N * gamma_nt
#    7b: 1/N * gamma_nt
#    7c: 1/N * (gamma_nt * Gamma^{-1} * gamma_nt)
#    7d: 1/N * Gamma (sandwich)


# YR 25 Mar 2016: first version
# YR 19 Jan 2017: added 6) + 7)
# YR 04 Jan 2020: adjust for sum(wt) != N
# YR 22 Jul 2022: adding correlation= argument for information_expected
#                 (only for catml; not used if correlation = TRUE!)

# 1. log-likelihood h1

# 1a: input is raw data
lav_mvn_h1_loglik_data <- function(
    y = NULL,
    x_idx = integer(0L),
    casewise = FALSE,
    wt = NULL,
    sinv_method = "eigen") {

  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }
  p <- NCOL(y)

  # sample statistics
  if (!is.null(wt)) {
    out <- stats::cov.wt(y, wt = wt, method = "ML")
    sample_mean <- out$center
    sample_cov <- out$cov
  } else {
    sample_mean <- base::.colMeans(y, m = n, n = p)
    sample_cov <- lav_mat_cov(y)
  }

  if (casewise) {
    log_2pi <- log(2 * pi)

    # invert sample_cov
    if (sinv_method == "chol") {
      c_s <- chol(sample_cov)
      ic_s <- backsolve(c_s, diag(p))
      yc <- t(t(y) - sample_mean)
      dist_1 <- rowSums((yc %*% ic_s)^2)
      logdet <- -2 * sum(log(diag(ic_s)))
    } else {
      sample_cov_inv <- lav_mat_sym_inverse(
        s = sample_cov,
        logdet = TRUE, sinv_method = sinv_method
      )
      logdet <- attr(sample_cov_inv, "logdet")
      # mahalanobis distance
      yc <- t(t(y) - sample_mean)
      dist_1 <- rowSums(yc %*% sample_cov_inv * yc)
    }

    loglik <- -(p * log_2pi + logdet + dist_1) / 2

    # weights
    if (!is.null(wt)) {
      loglik <- loglik * wt
    }
  } else {
    # invert sample_cov
    sample_cov_inv <- lav_mat_sym_inverse(
      s = sample_cov,
      logdet = TRUE, sinv_method = sinv_method
    )
    logdet <- attr(sample_cov_inv, "logdet")

    loglik <-
      lav_mvn_h1_loglik_samp(
        sample_cov_logdet = logdet,
        sample_nvar = p,
        sample_nobs = n
      )
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    loglik_x <- lav_mvn_h1_loglik_data(
      y = y[, x_idx, drop = FALSE],
      wt = wt, x_idx = integer(0L),
      casewise = casewise,
      sinv_method = sinv_method
    )
    # subtract logl.X
    loglik <- loglik - loglik_x
  }

  loglik
}



# 1b: input are sample statistics only (logdet, N and P)
lav_mvn_h1_loglik_samp <- function(
    sample_cov_logdet = NULL,
    sample_nvar = NULL,
    sample_nobs = NULL,
    # or
    sample_cov = NULL,
    x_idx = integer(0L),
    x_cov = NULL,
    sinv_method = "eigen") {

  if (is.null(sample_nvar)) {
    p <- NCOL(sample_cov)
  } else {
    p <- sample_nvar # number of variables
  }

  n <- sample_nobs
  stopifnot(!is.null(p), !is.null(n))

  log_2pi <- log(2 * pi)

  # all we need is the logdet
  if (is.null(sample_cov_logdet)) {
    sample_cov_inv <- lav_mat_sym_inverse(
      s = sample_cov,
      logdet = TRUE, sinv_method = sinv_method
    )
    logdet <- attr(sample_cov_inv, "logdet")
  } else {
    logdet <- sample_cov_logdet
  }

  loglik <- -n / 2 * (p * log_2pi + logdet + p)

  # fixed.x?
  if (length(x_idx) > 0L) {
    if (is.null(sample_cov)) {
      if (is.null(x_cov)) {
        lav_msg_stop(gettext(
          "when x_idx is not empty, we need sample_cov or x.cov"
        ))
      } else {
        sample_cov_x <- x_cov
      }
    } else {
      sample_cov_x <- sample_cov[x_idx, x_idx, drop = FALSE]
    }

    loglik_x <-
      lav_mvn_h1_loglik_samp(
        sample_cov = sample_cov_x,
        sample_nobs = sample_nobs,
        x_idx = integer(0L),
        sinv_method = sinv_method
      )
    # subtract logl.X
    loglik <- loglik - loglik_x
  }

  loglik
}


# 4. hessian of logl (around MLEs of Mu and Sigma)

# 4a: hessian logl Mu and vech(Sigma) from raw data
lav_mvn_h1_logl_hessian_data <- function(
    y = NULL,
    wt = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen",
    sample_cov_inv = NULL,
    meanstructure = TRUE) {

  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }

  # observed information
  observed <- lav_mvn_h1_info_observed_data(
    y = y, wt = wt,
    x_idx = x_idx,
    sinv_method = sinv_method,
    sample_cov_inv = sample_cov_inv,
    meanstructure = meanstructure
  )

  -n * observed
}

# 4b: hessian Mu and vech(Sigma) from samplestats
lav_mvn_h1_logl_hessian_samp <- function(
    sample_mean = NULL, # unused!
    sample_cov = NULL,
    sample_nobs = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen",
    sample_cov_inv = NULL,
    meanstructure = TRUE) {

  n <- sample_nobs

  # observed information
  observed <- lav_mvn_h1_info_observed_samp(
    sample_mean = sample_mean, sample_cov = sample_cov,
    x_idx = x_idx,
    sinv_method = sinv_method, sample_cov_inv = sample_cov_inv,
    meanstructure = meanstructure
  )

  -n * observed
}



# 5) Information h1 (note: expected == observed if data is complete!)

# 5a: unit expected information h1
lav_mvn_h1_info_expected <- function(
    y = NULL,
    wt = NULL,
    sample_cov = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen",
    sample_cov_inv = NULL,
    meanstructure = TRUE,
  correlation = FALSE) {

  if (is.null(sample_cov_inv)) {
    if (is.null(sample_cov)) {
      if (is.null(wt)) {
        sample_mean <- base::.colMeans(y, m = NROW(y), n = NCOL(y))
        sample_cov <- lav_mat_cov(y)
      } else {
        out <- stats::cov.wt(y, wt = wt, method = "ML")
        sample_cov <- out$cov
      }
    }

    # invert sample_cov
    sample_cov_inv <- lav_mat_sym_inverse(
      s = sample_cov,
      logdet = FALSE, sinv_method = sinv_method
    )
  }

  i11 <- sample_cov_inv
  if (correlation) {
    # if (lav_use_lavaanC()) {
    #   i22 <- lavaanC::m_kronecker_dup_cor_pre_post(sample.cov.inv,
    #                                                multiplicator = 0.5)
    # } else {
      i22 <- 0.5 * lav_mat_dup_cor_pre_post(sample_cov_inv %x%
                                                     sample_cov_inv)
    # }
  } else {
    # if (lav_use_lavaanC()) {
    #   i22 <- lavaanC::m_kronecker_dup_pre_post(sample.cov.inv,
    #                                            multiplicator = 0.5)
    # } else {
      i22 <- 0.5 * lav_mat_dup_pre_post(sample_cov_inv %x%
                                                   sample_cov_inv)
    # }
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    pstar_x <- lav_mat_vech_which_idx(
      n = NCOL(sample_cov_inv), idx = x_idx
    )
    i22[pstar_x, ] <- 0
    i22[, pstar_x] <- 0
  }

  if (meanstructure) {
    # fixed.x?
    if (length(x_idx) > 0L) {
      i11[x_idx, ] <- 0
      i11[, x_idx] <- 0
    }
    out <- lav_mat_bdiag(i11, i22)
  } else {
    out <- i22
  }

  out
}

# 5b: unit observed information h1
lav_mvn_h1_info_observed_data <- function(
    y = NULL,
    wt = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen",
    sample_cov_inv = NULL,
    meanstructure = TRUE) {

  lav_mvn_h1_info_expected(
    y = y, sinv_method = sinv_method,
    wt = wt, x_idx = x_idx,
    sample_cov_inv = sample_cov_inv,
    meanstructure = meanstructure
  )
}

# 5b-bis: observed information h1 from sample statistics
lav_mvn_h1_info_observed_samp <- function(
    sample_mean = NULL, # unused!
    sample_cov = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen",
    sample_cov_inv = NULL,
    meanstructure = TRUE) {

  if (is.null(sample_cov_inv)) {
    # invert sample_cov
    sample_cov_inv <- lav_mat_sym_inverse(
      s = sample_cov,
      logdet = FALSE, sinv_method = sinv_method
    )
  }

  i11 <- sample_cov_inv
  # fixed.x?
  if (length(x_idx) > 0L) {
    i11[x_idx, ] <- 0
    i11[, x_idx] <- 0
  }

  # if (lav_use_lavaanC()) {
  #   i22 <- lavaanC::m_kronecker_dup_pre_post(sample.cov.inv,
  #                                            multiplicator = 0.5)
  # } else {
    i22 <- 0.5 * lav_mat_dup_pre_post(sample_cov_inv %x%
    sample_cov_inv)
  # }

  # fixed.x?
  if (length(x_idx) > 0L) {
    pstar_x <- lav_mat_vech_which_idx(
      n = NCOL(sample_cov_inv), idx = x_idx
    )
    i22[pstar_x, ] <- 0
    i22[, pstar_x] <- 0
  }

  if (meanstructure) {
    out <- lav_mat_bdiag(i11, i22)
  } else {
    out <- i22
  }

  out
}

# 5c: unit first-order information h1
#     note: first order information h1 == A1 %*% Gamma %*% A1
#           (where A1 = obs/exp information h1)
lav_mvn_h1_info_firstorder <- function(
    y = NULL,
    wt = NULL,
    sample_cov = NULL,
    x_idx = integer(0L),
    cluster_idx = NULL,
    sinv_method = "eigen",
    sample_cov_inv = NULL,
    gamma_1 = NULL,
    meanstructure = TRUE) {

  if (!is.null(wt)) {
    out <- stats::cov.wt(y, wt = wt, method = "ML")
    res <- lav_mvn_info_firstorder(
      y = y, wt = wt,
      cluster_idx = cluster_idx,
      mu = out$center, sigma_1 = out$cov, x_idx = x_idx,
      meanstructure = meanstructure
    )
    return(res)
  }

   # sample.cov.inv
  if (is.null(sample_cov_inv)) {
    # invert sample_cov
    if (is.null(sample_cov)) {
      sample_mean <- base::.colMeans(y, m = NROW(y), n = NCOL(y))
      sample_cov <- lav_mat_cov(y)
    }
    sample_cov_inv <- lav_mat_sym_inverse(
      s = sample_cov,
      logdet = FALSE, sinv_method = sinv_method
    )
  }

  # question: is there any benefit computing Gamma/A1 instead of just
  # calling lav_mvn_info_firstorder()?
  # answer (2014): probably not; it is just reassuring that the expression
  # J = A1 %*% gamma_1 %*% A1 seems to hold

  # Gamma
  # FIXME: what about the 'unbiased = TRUE' option?
  if (is.null(gamma_1)) {
    if (length(x_idx) > 0L) {
      gamma_1 <- lav_samp_gamma(y,
        x_idx = x_idx, fixed_x = TRUE,
        cluster_idx = cluster_idx,
        meanstructure = meanstructure
      )
    } else {
      gamma_1 <- lav_samp_gamma(y,
        meanstructure = meanstructure,
        cluster_idx = cluster_idx
      )
    }
  }

  # sample.cov.inv
  if (is.null(sample_cov_inv)) {
    # invert sample_cov
    if (is.null(sample_cov)) {
      sample_mean <- base::.colMeans(y, m = NROW(y), n = NCOL(y))
      sample_cov <- lav_mat_cov(y)
    }
    sample_cov_inv <- lav_mat_sym_inverse(
      s = sample_cov,
      logdet = FALSE, sinv_method = sinv_method
    )
  }

  # A1
  a1 <- lav_mvn_h1_info_expected(
    y = y, sinv_method = sinv_method,
    sample_cov_inv = sample_cov_inv,
    x_idx = x_idx,
    meanstructure = meanstructure
  )

  a1 %*% gamma_1 %*% a1
}

# 6) inverted information h1 mu + vech(Sigma) (not used?)

#    6a: (unit) inverted expected information (A1.inv = gamma_nt)
#    6b: (unit) inverted observed information (A1.inv = gamma_nt)

lav_mvnorm_h1_inverted_information_expected <-              # nolint
lav_mvn_h1_inv_info_observed <- function(
  y = NULL,
  wt = NULL,
  sample_cov = NULL,
  x_idx = integer(0L)) {

  # sample_cov
  if (is.null(sample_cov)) {
    if (is.null(wt)) {
      sample_mean <- base::.colMeans(y, m = NROW(y), n = NCOL(y))
      sample_cov <- lav_mat_cov(y)
    } else {
      out <- stats::cov.wt(y, wt = wt, method = "ML")
      sample_cov <- out$cov
    }
  }

  if (length(x_idx) > 0L) {
    gamma_nt <- lav_samp_gamma_nt(
      m_y = y, wt = wt, x_idx = x_idx,
      m_cov = sample_cov,
      meanstructure = TRUE,
      fixed_x = TRUE
    )
  } else {
    i11 <- sample_cov
    # if (lav_use_lavaanC()) {
    #   i22 <- lavaanC::m_kronecker_dup_ginv_pre_post(sample_cov,
    #                                                 multiplicator = 2.0)
    # } else {
      i22 <- 2 * lav_mat_dup_ginv_pre_post(sample_cov %x% sample_cov)
    # }
    gamma_nt <- lav_mat_bdiag(i11, i22)
  }

  gamma_nt
}

#    6c: (unit) inverted first-order information (B1.inv) (not used?)
#        J1.inv = gamma_nt %*% solve(Gamma) %*%  gamma_nt
#
lav_mvn_h1_inv_info_firstorder <- function(
    y = NULL,
    wt = NULL,
    sample_cov = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen",
    sample_cov_inv = NULL,
    gamma_1 = NULL) {

  # lav_samp_gamma() has no wt argument (yet)
  if (!is.null(wt)) {
    lav_msg_stop(gettext("function not supported if wt is not NULL"))
  }

  # Gamma
  # what about the 'unbiased = TRUE' option?
  if (is.null(gamma_1)) {
    if (length(x_idx) > 0L) {
      gamma_1 <- lav_samp_gamma(y,
        x_idx = x_idx, fixed_x = TRUE,
        meanstructure = TRUE
      )
    } else {
      gamma_1 <- lav_samp_gamma(y, meanstructure = TRUE)
    }
  }

  # gamma_nt
  gamma_nt <-
    lav_mvnorm_h1_inverted_information_expected(
      y = y,
      sample_cov = sample_cov,
      x_idx = x_idx
    )
  if (length(x_idx) > 0L) {
    # FIXME: surely there is better way
    out <- gamma_nt %*% MASS::ginv(gamma_1) %*% gamma_nt
  } else {
    out <- gamma_nt %*% solve(gamma_1, gamma_nt)
  }

  out
}


# 7) ACOV h1 mu + vech(Sigma) (not used?)

#    7a: 1/N * gamma_nt
#    7b: 1/N * gamma_nt
lav_mvn_h1_acov_expected <-
lav_mvn_h1_acov_observed <- function(
    y = NULL,
    wt = NULL,
    sample_cov = NULL,
    x_idx = integer(0L)) {

  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }

  gamma_nt <-
    lav_mvnorm_h1_inverted_information_expected(
      y = y,
      wt = wt,
      sample_cov = sample_cov,
      x_idx = x_idx
    )

  (1 / n) * gamma_nt
}

#    7c: 1/N * (gamma_nt * Gamma^{-1} * gamma_nt)
lav_mvn_h1_acov_firstorder <- function(
    y = NULL,
    wt = NULL,
    sample_cov = NULL,
    sinv_method = "eigen",
    x_idx = integer(0L),
    sample_cov_inv = NULL,
    gamma_1 = NULL) {

  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }

  j1_inv <- lav_mvn_h1_inv_info_firstorder(
    y = y, wt = wt,
    sample_cov = sample_cov,
    x_idx = x_idx, sinv_method = sinv_method,
    sample_cov_inv = sample_cov_inv, gamma_1 = gamma_1
  )

  (1 / n) * j1_inv
}

#    7d: 1/N * Gamma (sandwich)
lav_mvn_h1_acov_sandwich <- function(
    y = NULL,
    wt = NULL,
    sample_cov = NULL,
    x_idx = integer(0L),
    gamma_1 = NULL) {

  # lav_samp_gamma() has no wt argument (yet)
  if (!is.null(wt)) {
    lav_msg_stop(gettext("function not supported if wt is not NULL"))
  }

  # if(!is.null(wt)) {
  #    N <- sum(wt)
  # } else {
  n <- NROW(y)
  # }

  # Gamma
  if (is.null(gamma_1)) {
    if (length(x_idx) > 0L) {
      gamma_1 <- lav_samp_gamma(y,
        x_idx = x_idx, fixed_x = TRUE,
        meanstructure = TRUE
      )
    } else {
      gamma_1 <- lav_samp_gamma(y, meanstructure = TRUE)
    }
  }

  (1 / n) * gamma_1
}
