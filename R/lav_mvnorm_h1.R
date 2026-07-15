# the multivariate normal distribution, unrestricted (h1)
# - everything is evaluated under the MLEs: Mu = ybar, Sigma = S

# 1) loglikelihood h1 (from raw data, or sample statistics)
# 4) hessian h1 around MLEs
# 5) information h1 (restricted Sigma/mu)
#    5a: (unit)    expected information h1 (A1 = gamma_nt^{-1})
#    5b: (unit)    observed information h1 (A1 = gamma_nt^{-1})
#    5c: (unit) first.order information h1 (B1 = A1 %*% Gamma %*% A1)


# YR 25 Mar 2016: first version
# YR 19 Jan 2017: added 6) + 7)
# YR 04 Jan 2020: adjust for sum(wt) != N
# YR 22 Jul 2022: adding correlation= argument for information_expected
#                 (only for catml; not used if correlation = TRUE!)
# YR 03 Jul 2026: use shared kernels from lav_mvnorm_kernels.R
# YR 08 Jul 2026: removed the unused 6) inverted-information and 7) ACOV
#                 families (see note at the end of this file)

# 1. log-likelihood h1

# 1a: input is raw data
lav_mvn_h1_loglik_data <- function(
    y = NULL,
    x_idx = integer(0L),
    casewise = FALSE,
    wt = NULL,
    sinv_method = "eigen") {

  n <- lav_mvn_nobs(y = y, wt = wt)
  p <- NCOL(y)

  # sample statistics
  samp <- lav_mvn_samp_stats(y = y, wt = wt)
  sample_mean <- samp$mean
  sample_cov <- samp$cov

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

  n <- lav_mvn_nobs(y = y, wt = wt)

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
      sample_cov <- lav_mvn_samp_stats(y = y, wt = wt)$cov
    }

    # invert sample_cov
    sample_cov_inv <- lav_mat_sym_inverse(
      s = sample_cov,
      logdet = FALSE, sinv_method = sinv_method
    )
  }

  i22 <- lav_mvn_kron_dup_half(sample_cov_inv, correlation = correlation)

  if (meanstructure) {
    out <- lav_mat_bdiag(sample_cov_inv, i22)
  } else {
    out <- i22
  }

  # fixed.x?
  out <- lav_mvn_zero_x_idx(out,
    n = NCOL(sample_cov_inv), x_idx = x_idx,
    meanstructure = meanstructure
  )

  out
}

# 5b: unit observed information h1
#     (identical to the expected information, as the data is complete
#      and we evaluate at the MLEs Mu = ybar, Sigma = S;
#      note: the correlation= argument is NOT forwarded here)
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

  sample_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = sample_cov, sigma_inv = sample_cov_inv,
    sinv_method = sinv_method
  )

  i22 <- lav_mvn_kron_dup_half(sample_cov_inv)

  if (meanstructure) {
    out <- lav_mat_bdiag(sample_cov_inv, i22)
  } else {
    out <- i22
  }

  # fixed.x?
  out <- lav_mvn_zero_x_idx(out,
    n = NCOL(sample_cov_inv), x_idx = x_idx,
    meanstructure = meanstructure
  )

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

  # A1
  a1 <- lav_mvn_h1_info_expected(
    y = y, sinv_method = sinv_method,
    sample_cov_inv = sample_cov_inv,
    x_idx = x_idx,
    meanstructure = meanstructure
  )

  a1 %*% gamma_1 %*% a1
}

# NOTE: the former 6) inverted-information and 7) ACOV families (gamma_nt,
# B1.inv, 1/N * Gamma, ...) were unused reference implementations of the
# identity Gamma_ADF = A1^{-1} J1 A1^{-1}; they were removed in 0.7-1.
# The saturated-model Gamma itself lives in lav_samplestats_gamma.R
# (lav_samp_gamma / lav_samp_gamma_nt), and Gamma/N for a fitted model in
# lav_model_h1_acov() (lav_model_h1_information.R).
