# the multivariate linear model using maximum likelihood

# 1) loglikelihood (from raw data, or sample statistics)
# 2) derivatives with respect to Beta, res.cov, vech(res.cov)
# 3) casewise scores with respect to Beta, vech(res.cov), Beta + vech(res.cov)
# 4) hessian Beta + vech(res.cov)
# 5) information h0 Beta + vech(res.cov)
#    5a: (unit)    expected information
#    5b: (unit)    observed information
#    5c: (unit) first.order information

# YR 24 Mar 2016: first version
# YR 20 Jan 2017: removed added 'N' in many equations, to be consistent with
#                 lav_mvnorm_*
# YR 18 Okt 2018: add 'information' functions, change arguments
#                 (X -> eXo, Sigma -> res.cov, Beta -> res.int + res.slopes)

# 1. loglikelihood

# 1a. input is raw data
lav_mvreg_loglik_data <- function(y = NULL,
                                  exo = NULL, # no intercept
                                  beta_1 = NULL,
                                  res_int = NULL,
                                  res_slopes = NULL,
                                  res_cov = NULL,
                                  casewise = FALSE,
                                  sinv_method = "eigen") {
  y <- unname(y)
  q_1 <- NCOL(y)
  n <- NROW(y)
  x <- cbind(1, unname(exo))

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  if (casewise) {
    log_2pi <- log(2 * pi)

    # invert res.cov
    if (sinv_method == "chol") {
      c_s <- chol(res_cov)
      ic_s <- backsolve(c_s, diag(q_1))
      logdet <- -2 * sum(log(diag(ic_s)))

      res <- y - x %*% beta_1
      dist_1 <- rowSums((res %*% ic_s)^2)
    } else {
      res_cov_inv <- lav_matrix_symmetric_inverse(
        s = res_cov, logdet = TRUE,
        sinv_method = sinv_method
      )
      logdet <- attr(res_cov_inv, "logdet")

      res <- y - x %*% beta_1
      dist_1 <- rowSums(res %*% res_cov_inv * res)
    }

    loglik <- -(q_1 * log_2pi + logdet + dist_1) / 2
  } else {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = TRUE,
      sinv_method = sinv_method
    )
    logdet <- attr(res_cov_inv, "logdet")

    res <- y - x %*% beta_1
    # TOTAL <- TR( (Y - X%*%Beta) %*% res.cov.inv %*% t(Y - X%*%Beta) )
    total <- sum(rowSums(res %*% res_cov_inv * res))
    loglik <- -(n * q_1 / 2) * log(2 * pi) - (n / 2) * logdet - (1 / 2) * total
  }

  loglik
}


# 2b. input are sample statistics (res.int, res.slopes, res.cov, N) only
lav_mvreg_loglik_samplestats <- function(sample_res_int = NULL,
                                         sample_res_slopes = NULL,
                                         sample_res_cov = NULL,
                                         sample_mean_x = NULL,
                                         sample_cov_x = NULL,
                                         sample_nobs = NULL,
                                         beta_1 = NULL, # optional
                                         res_int = NULL,
                                         res_slopes = NULL,
                                         res_cov = NULL,
                                         sinv_method = "eigen",
                                         res_cov_inv = NULL) {
  q_1 <- NCOL(sample_res_cov)
  n <- sample_nobs
  # log_2pi <- log(2 * pi)

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  # construct 'saturated' (sample-based) B
  sample_b <- rbind(matrix(sample_res_int, nrow = 1), t(sample_res_slopes))

  # construct sample.xx = 1/N*crossprod(X1) (including intercept)
  sample_xx <- rbind(
    cbind(1, matrix(sample_mean_x, nrow = 1, )),
    cbind(
      matrix(sample_mean_x, ncol = 1),
      sample_cov_x + tcrossprod(sample_mean_x)
    )
  )

  # res.cov.inv
  if (is.null(res_cov_inv)) {
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = TRUE,
      sinv_method = sinv_method
    )
    logdet <- attr(res_cov_inv, "logdet")
  } else {
    logdet <- attr(res_cov_inv, "logdet")
    if (is.null(logdet)) {
      # compute - ln|res.cov.inv|
      ev <- eigen(res_cov_inv, symmetric = TRUE, only.values = TRUE)
      logdet <- -1 * sum(log(ev$values))
    }
  }

  # tr(res.cov^{-1} %*% S)
  dist1 <- sum(res_cov_inv * sample_res_cov)

  # tr( res.cov^{-1} (B-beta)' X'X (B-beta)
  diff_1 <- sample_b - beta_1
  dist2 <- sum(res_cov_inv * crossprod(diff_1, sample_xx) %*% diff_1)

  loglik <- -(n / 2) * (q_1 * log(2 * pi) + logdet + dist1 + dist2)

  loglik
}




# 2. Derivatives

# 2a. derivative logl with respect to Beta (=intercepts and slopes)
lav_mvreg_dlogl_dbeta <- function(y = NULL,
                                  exo = NULL,
                                  beta_1 = NULL,
                                  res_int = NULL,
                                  res_slopes = NULL,
                                  res_cov = NULL,
                                  sinv_method = "eigen",
                                  res_cov_inv = NULL) {
  y <- unname(y)
  x <- cbind(1, unname(exo))

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  # res.cov.inv
  if (is.null(res_cov_inv)) {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # subtract 'X %*% Beta' from Y
  res <- y - x %*% beta_1

  # derivative
  dbeta <- as.numeric(t(x) %*% res %*% res_cov_inv)

  dbeta
}

# 2b: derivative logl with respect to res.cov (full matrix, ignoring symmetry)
lav_mvreg_dlogl_drescov <- function(y = NULL,
                                    exo = NULL,
                                    beta_1 = NULL,
                                    res_cov = NULL,
                                    res_int = NULL,
                                    res_slopes = NULL,
                                    sinv_method = "eigen",
                                    res_cov_inv = NULL) {
  y <- unname(y)
  n <- NROW(y)
  x <- cbind(1, unname(exo))

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  # res.cov.in
  if (is.null(res_cov_inv)) {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # subtract 'X %*% Beta' from Y
  res <- y - x %*% beta_1

  # W.tilde
  w_tilde <- crossprod(res) / n

  # derivative
  dres_cov <- -(n / 2) * (res_cov_inv -
                         (res_cov_inv %*% w_tilde %*% res_cov_inv))

  dres_cov
}

# 2c: derivative logl with respect to vech(res.cov)
lav_mvreg_dlogl_dvechrescov <- function(y = NULL,
                                        exo = NULL,
                                        beta_1 = NULL,
                                        res_int = NULL,
                                        res_slopes = NULL,
                                        res_cov = NULL,
                                        sinv_method = "eigen",
                                        res_cov_inv = NULL) {
  y <- unname(y)
  n <- NROW(y)
  x <- cbind(1, unname(exo))

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  # res.cov.inv
  if (is.null(res_cov_inv)) {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # subtract 'X %*% Beta' from Y
  res <- y - x %*% beta_1

  # W.tilde
  w_tilde <- crossprod(res) / n

  # derivative
  dres_cov <- -(n / 2) * (res_cov_inv -
                         (res_cov_inv %*% w_tilde %*% res_cov_inv))
  dvechres_cov <- as.numeric(lav_matrix_duplication_pre(
    as.matrix(lav_matrix_vec(dres_cov))
  ))

  dvechres_cov
}


# 3. Casewise scores

# 3a: casewise scores with respect to Beta (=intercepts and slopes)
#     column order: Y1_int, Y1_x1, Y1_x2, ...| Y2_int, Y2_x1, Y2_x2, ... |
lav_mvreg_scores_beta <- function(y = NULL,
                                  exo = NULL,
                                  beta_1 = NULL,
                                  res_int = NULL,
                                  res_slopes = NULL,
                                  res_cov = NULL,
                                  sinv_method = "eigen",
                                  res_cov_inv = NULL) {
  y <- unname(y)
  q_1 <- NCOL(y)
  x <- cbind(1, unname(exo))
  p <- NCOL(x)

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  # res.cov.inv
  if (is.null(res_cov_inv)) {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # subtract Mu
  res <- y - x %*% beta_1

  # post-multiply with res.cov.inv
  res <- res %*% res_cov_inv

  sc_beta <- x[, rep(1:p, times = q_1), drop = FALSE] *
    res[, rep(1:q_1, each = p), drop = FALSE]

  sc_beta
}


# 3b: casewise scores with respect to vech(res.cov)
lav_mvreg_scores_vech_sigma <- function(y = NULL,
                                        exo = NULL,
                                        beta_1 = NULL,
                                        res_int = NULL,
                                        res_slopes = NULL,
                                        res_cov = NULL,
                                        sinv_method = "eigen",
                                        res_cov_inv = NULL) {
  y <- unname(y)
  q_1 <- NCOL(y)
  x <- cbind(1, unname(exo))

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  # res.cov.inv
  if (is.null(res_cov_inv)) {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # vech(res.cov.inv)
  isigma <- lav_matrix_vech(res_cov_inv)

  # subtract X %*% Beta
  res <- y - x %*% beta_1

  # postmultiply with res.cov.inv
  res <- res %*% res_cov_inv

  # tcrossprod
  idx1 <- lav_matrix_vech_col_idx(q_1)
  idx2 <- lav_matrix_vech_row_idx(q_1)
  z <- res[, idx1] * res[, idx2]

  # subtract isigma from each row
  sc <- t(t(z) - isigma)

  # adjust for vech (and avoiding the 1/2 factor)
  sc[, lav_matrix_diagh_idx(q_1)] <- sc[, lav_matrix_diagh_idx(q_1)] / 2

  sc
}


# 3c: casewise scores with respect to beta + vech(res.cov)
lav_mvreg_scores_beta_vech_sigma <- function(y = NULL,      # nolint
                                             exo = NULL,
                                             beta_1 = NULL,
                                             res_int = NULL,
                                             res_slopes = NULL,
                                             res_cov = NULL,
                                             sinv_method = "eigen",
                                             res_cov_inv = NULL) {
  y <- unname(y)
  q_1 <- NCOL(y)
  x <- cbind(1, unname(exo))
  p <- NCOL(x)

  # construct model-implied Beta
  if (is.null(beta_1)) {
    beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
  }

  # res.cov.inv
  if (is.null(res_cov_inv)) {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # vech(res.cov.inv)
  isigma <- lav_matrix_vech(res_cov_inv)

  # subtract X %*% Beta
  res <- y - x %*% beta_1

  # postmultiply with res.cov.inv
  res <- res %*% res_cov_inv

  sc_beta <- x[, rep(1:p, times = q_1), drop = FALSE] *
    res[, rep(1:q_1, each = p), drop = FALSE]

  # tcrossprod
  idx1 <- lav_matrix_vech_col_idx(q_1)
  idx2 <- lav_matrix_vech_row_idx(q_1)
  z <- res[, idx1] * res[, idx2]

  # subtract isigma from each row
  sc <- t(t(z) - isigma)

  # adjust for vech (and avoiding the 1/2 factor)
  sc[, lav_matrix_diagh_idx(q_1)] <- sc[, lav_matrix_diagh_idx(q_1)] / 2

  cbind(sc_beta, sc)
}

# 4. hessian of logl

# 4a. hessian logl Beta and vech(res.cov) from raw data
lav_mvreg_logl_hessian_data <- function(y = NULL,
                                        exo = NULL, # no int
                                        beta_1 = NULL, # int+slopes
                                        res_int = NULL,
                                        res_slopes = NULL,
                                        res_cov = NULL,
                                        res_cov_inv = NULL,
                                        sinv_method = "eigen") {
  # sample size
  n <- NROW(y)

  # observed information
  observed <- lav_mvreg_information_observed_data(
    y = y, exo = exo,
    beta_1 = beta_1, res_int = res_int, res_slopes = res_slopes,
    res_cov = res_cov, res_cov_inv = res_cov_inv,
    sinv_method = sinv_method
  )

  # hessian
  -n * observed
}

# 4b. hessian logl Beta and vech(res.cov) from samplestats
lav_mvreg_logl_hessian_samplestats <- function(  # nolint
      sample_res_int = NULL,
      sample_res_slopes = NULL,
      sample_res_cov = NULL,
      sample_mean_x = NULL,
      sample_cov_x = NULL,
      sample_nobs = NULL,
      beta_1 = NULL, # int + slopes
      res_int = NULL, # intercepts only
      res_slopes = NULL, # slopes only (y x x)
      res_cov = NULL, # res.cov
      sinv_method = "eigen",
      res_cov_inv = NULL) {
  # sample size
  n <- sample_nobs

  # information
  observed <- lav_mvreg_information_observed_samplestats(
    sample_res_int = sample_res_int, sample_res_slopes = sample_res_slopes,
    sample_res_cov = sample_res_cov, sample_mean_x = sample_mean_x,
    sample_cov_x = sample_cov_x, beta_1 = beta_1, res_int = res_int,
    res_slopes = res_slopes, res_cov = res_cov, sinv_method = sinv_method,
    res_cov_inv = res_cov_inv
  )

  # hessian
  -n * observed
}


# Information h0

# 5a: unit expected information h0 Beta and vech(res.cov)
lav_mvreg_information_expected <- function(y = NULL, # not used
                                           exo = NULL, # not used
                                           sample_mean_x = NULL,
                                           sample_cov_x = NULL,
                                           sample_nobs = NULL,
                                           beta_1 = NULL, # not used
                                           res_int = NULL, # not used
                                           res_slopes = NULL, # not used
                                           res_cov = NULL,
                                           res_cov_inv = NULL,
                                           sinv_method = "eigen") {
  exo <- unname(exo)

  # res.cov.inv
  if (is.null(res_cov_inv)) {
    # invert res.cov
    res_cov_inv <- lav_matrix_symmetric_inverse(
      s = res_cov, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # N
  if (is.null(sample_nobs)) {
    sample_nobs <- nrow(exo) # hopefully not NULL either
  } else {
    # n <- sample_nobs
  }

  # sample.mean.x + sample.cov.x
  if (is.null(sample_mean_x)) {
    sample_mean_x <- base::.colMeans(exo, m = NROW(exo), n = NCOL(exo))
  }
  if (is.null(sample_cov_x)) {
    sample_cov_x <- lav_matrix_cov(exo)
  }

  # construct sample.xx = 1/N*crossprod(X1) (including intercept)
  sample_xx <- rbind(
    cbind(1, matrix(sample_mean_x, nrow = 1, )),
    cbind(
      matrix(sample_mean_x, ncol = 1),
      sample_cov_x + tcrossprod(sample_mean_x)
    )
  )

  # expected information
  i11 <- res_cov_inv %x% sample_xx
  # if (lav_use_lavaanC()) {
  #   I22 <- lavaanC::m_kronecker_dup_pre_post(res.cov.inv, multiplicator = 0.5)
  # } else {
    i22 <- 0.5 * lav_matrix_duplication_pre_post(res_cov_inv %x% res_cov_inv)
  # }

  lav_matrix_bdiag(i11, i22)
}

# 5b: unit observed information h0
lav_mvreg_information_observed_data <- function(y = NULL,          # nolint
                                                exo = NULL, # no int
                                                beta_1 = NULL, # int+slopes
                                                res_int = NULL,
                                                res_slopes = NULL,
                                                res_cov = NULL,
                                                res_cov_inv = NULL,
                                                sinv_method = "eigen") {
  # create sample statistics
  y <- unname(y)
  x1 <- cbind(1, unname(exo))
  n <- NROW(y)

  # find 'B'
  qr_1 <- qr(x1)
  sample_b <- qr.coef(qr_1, y)

  sample_res_int <- as.numeric(sample_b[1, ])
  sample_res_slopes <- t(sample_b[-1, , drop = FALSE]) # transpose!
  sample_res_cov <- cov(qr.resid(qr_1, y)) * (n - 1) / n
  sample_mean_x <- base::.colMeans(exo, m = NROW(exo), n = NCOL(exo))
  sample_cov_x <- lav_matrix_cov(exo)

  lav_mvreg_information_observed_samplestats(
    sample_res_int = sample_res_int,
    sample_res_slopes = sample_res_slopes, sample_res_cov = sample_res_cov,
    sample_mean_x = sample_mean_x, sample_cov_x = sample_cov_x,
    beta_1 = beta_1, res_int = res_int, res_slopes = res_slopes,
    res_cov = res_cov, sinv_method = sinv_method, res_cov_inv = res_cov_inv
  )
}


# 5b-bis: observed information h0 from sample statistics
lav_mvreg_information_observed_samplestats <-               # nolint
  function(sample_res_int = NULL,
           sample_res_slopes = NULL,
           sample_res_cov = NULL,
           sample_mean_x = NULL,
           sample_cov_x = NULL,
           beta_1 = NULL, # int + slopes
           res_int = NULL, # intercepts only
           res_slopes = NULL, # slopes only (y x x)
           res_cov = NULL, # res.cov
           sinv_method = "eigen",
           res_cov_inv = NULL) {
    # construct model-implied Beta
    if (is.null(beta_1)) {
      beta_1 <- rbind(matrix(res_int, nrow = 1), t(res_slopes))
    }

    # construct 'saturated' (sample-based) B
    sample_b <- rbind(matrix(sample_res_int, nrow = 1), t(sample_res_slopes))

    # construct sample.xx = 1/N*crossprod(X1) (including intercept)
    sample_xx <- rbind(
      cbind(1, matrix(sample_mean_x, nrow = 1, )),
      cbind(
        matrix(sample_mean_x, ncol = 1),
        sample_cov_x + tcrossprod(sample_mean_x)
      )
    )

    # W.tilde = S + t(B - Beta) %*% (1/N)*X'X %*% (B - Beta)
    w_tilde <- (sample_res_cov +
      t(sample_b - beta_1) %*% sample_xx %*% (sample_b - beta_1))

    # res.cov.inv
    if (is.null(res_cov_inv)) {
      # invert res.cov
      res_cov_inv <- lav_matrix_symmetric_inverse(
        s = res_cov, logdet = FALSE,
        sinv_method = sinv_method
      )
    }

    h11 <- res_cov_inv %x% sample_xx
    h21 <- lav_matrix_duplication_pre(res_cov_inv %x%
      (res_cov_inv %*% (crossprod(sample_b - beta_1, sample_xx))))
    h12 <- t(h21)

    aaa <- res_cov_inv %*% (2 * w_tilde - res_cov) %*% res_cov_inv
    # if (lav_use_lavaanC()) {
    #   H22 <- lavaanC::m_kronecker_dup_pre_post(res.cov.inv, AAA, 0.5)
    # } else {
      h22 <- (1 / 2) * lav_matrix_duplication_pre_post(res_cov_inv %x% aaa)
    # }

    out <- rbind(
      cbind(h11, h12),
      cbind(h21, h22)
    )
    out
  }


# 5c: unit first-order information h0
lav_mvreg_information_firstorder <- function(y = NULL,        # nolint
                                             exo = NULL, # no int
                                             beta_1 = NULL, # int+slopes
                                             res_int = NULL,
                                             res_slopes = NULL,
                                             res_cov = NULL,
                                             res_cov_inv = NULL,
                                             sinv_method = "eigen") {
  n <- NROW(y)

  # scores
  sc <- lav_mvreg_scores_beta_vech_sigma(
    y = y, exo = exo, beta_1 = beta_1,
    res_int = res_int, res_slopes = res_slopes, res_cov = res_cov,
    sinv_method = sinv_method, res_cov_inv = res_cov_inv
  )

  crossprod(sc) / n
}


# 6: inverted information h0

# 6a: inverted unit expected information h0 Beta and vech(res.cov)
#
# lav_mvreg_inverted_information_expected <- function(Y       = NULL, # unused!
# }
