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
lav_mvreg_loglik_data <- function(Y = NULL,
                                  eXo = NULL, # no intercept
                                  Beta = NULL,
                                  res.int = NULL,
                                  res.slopes = NULL,
                                  res.cov = NULL,
                                  casewise = FALSE,
                                  Sinv.method = "eigen") {
  Y <- unname(Y)
  Q <- NCOL(Y)
  N <- NROW(Y)
  X <- cbind(1, unname(eXo))

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  if (casewise) {
    LOG.2PI <- log(2 * pi)

    # invert res.cov
    if (Sinv.method == "chol") {
      cS <- chol(res.cov)
      icS <- backsolve(cS, diag(Q))
      logdet <- -2 * sum(log(diag(icS)))

      RES <- Y - X %*% Beta
      DIST <- rowSums((RES %*% icS)^2)
    } else {
      res.cov.inv <- lav_matrix_symmetric_inverse(
        S = res.cov, logdet = TRUE,
        Sinv.method = Sinv.method
      )
      logdet <- attr(res.cov.inv, "logdet")

      RES <- Y - X %*% Beta
      DIST <- rowSums(RES %*% res.cov.inv * RES)
    }

    loglik <- -(Q * LOG.2PI + logdet + DIST) / 2
  } else {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = TRUE,
      Sinv.method = Sinv.method
    )
    logdet <- attr(res.cov.inv, "logdet")

    RES <- Y - X %*% Beta
    # TOTAL <- TR( (Y - X%*%Beta) %*% res.cov.inv %*% t(Y - X%*%Beta) )
    TOTAL <- sum(rowSums(RES %*% res.cov.inv * RES))
    loglik <- -(N * Q / 2) * log(2 * pi) - (N / 2) * logdet - (1 / 2) * TOTAL
  }

  loglik
}


# 2b. input are sample statistics (res.int, res.slopes, res.cov, N) only
lav_mvreg_loglik_samplestats <- function(sample.res.int = NULL,
                                         sample.res.slopes = NULL,
                                         sample.res.cov = NULL,
                                         sample.mean.x = NULL,
                                         sample.cov.x = NULL,
                                         sample.nobs = NULL,
                                         Beta = NULL, # optional
                                         res.int = NULL,
                                         res.slopes = NULL,
                                         res.cov = NULL,
                                         Sinv.method = "eigen",
                                         res.cov.inv = NULL) {
  Q <- NCOL(sample.res.cov)
  N <- sample.nobs
  LOG.2PI <- log(2 * pi)

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  # construct 'saturated' (sample-based) B
  sample.B <- rbind(matrix(sample.res.int, nrow = 1), t(sample.res.slopes))

  # construct sample.xx = 1/N*crossprod(X1) (including intercept)
  sample.xx <- rbind(
    cbind(1, matrix(sample.mean.x, nrow = 1, )),
    cbind(
      matrix(sample.mean.x, ncol = 1),
      sample.cov.x + tcrossprod(sample.mean.x)
    )
  )

  # res.cov.inv
  if (is.null(res.cov.inv)) {
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = TRUE,
      Sinv.method = Sinv.method
    )
    logdet <- attr(res.cov.inv, "logdet")
  } else {
    logdet <- attr(res.cov.inv, "logdet")
    if (is.null(logdet)) {
      # compute - ln|res.cov.inv|
      ev <- eigen(res.cov.inv, symmetric = TRUE, only.values = TRUE)
      logdet <- -1 * sum(log(ev$values))
    }
  }

  # tr(res.cov^{-1} %*% S)
  DIST1 <- sum(res.cov.inv * sample.res.cov)

  # tr( res.cov^{-1} (B-beta)' X'X (B-beta)
  Diff <- sample.B - Beta
  DIST2 <- sum(res.cov.inv * crossprod(Diff, sample.xx) %*% Diff)

  loglik <- -(N / 2) * (Q * log(2 * pi) + logdet + DIST1 + DIST2)

  loglik
}




# 2. Derivatives

# 2a. derivative logl with respect to Beta (=intercepts and slopes)
lav_mvreg_dlogl_dbeta <- function(Y = NULL,
                                  eXo = NULL,
                                  Beta = NULL,
                                  res.int = NULL,
                                  res.slopes = NULL,
                                  res.cov = NULL,
                                  Sinv.method = "eigen",
                                  res.cov.inv = NULL) {
  Y <- unname(Y)
  X <- cbind(1, unname(eXo))

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  # res.cov.inv
  if (is.null(res.cov.inv)) {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # substract 'X %*% Beta' from Y
  RES <- Y - X %*% Beta

  # derivative
  dbeta <- as.numeric(t(X) %*% RES %*% res.cov.inv)

  dbeta
}

# 2b: derivative logl with respect to res.cov (full matrix, ignoring symmetry)
lav_mvreg_dlogl_drescov <- function(Y = NULL,
                                    eXo = NULL,
                                    Beta = NULL,
                                    res.cov = NULL,
                                    res.int = NULL,
                                    res.slopes = NULL,
                                    Sinv.method = "eigen",
                                    res.cov.inv = NULL) {
  Y <- unname(Y)
  N <- NROW(Y)
  X <- cbind(1, unname(eXo))

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  # res.cov.in
  if (is.null(res.cov.inv)) {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # substract 'X %*% Beta' from Y
  RES <- Y - X %*% Beta

  # W.tilde
  W.tilde <- crossprod(RES) / N

  # derivative
  dres.cov <- -(N / 2) * (res.cov.inv - (res.cov.inv %*% W.tilde %*% res.cov.inv))

  dres.cov
}

# 2c: derivative logl with respect to vech(res.cov)
lav_mvreg_dlogl_dvechrescov <- function(Y = NULL,
                                        eXo = NULL,
                                        Beta = NULL,
                                        res.int = NULL,
                                        res.slopes = NULL,
                                        res.cov = NULL,
                                        Sinv.method = "eigen",
                                        res.cov.inv = NULL) {
  Y <- unname(Y)
  N <- NROW(Y)
  X <- cbind(1, unname(eXo))

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  # res.cov.inv
  if (is.null(res.cov.inv)) {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # substract 'X %*% Beta' from Y
  RES <- Y - X %*% Beta

  # W.tilde
  W.tilde <- crossprod(RES) / N

  # derivative
  dres.cov <- -(N / 2) * (res.cov.inv - (res.cov.inv %*% W.tilde %*% res.cov.inv))
  dvechres.cov <- as.numeric(lav_matrix_duplication_pre(
    as.matrix(lav_matrix_vec(dres.cov))
  ))

  dvechres.cov
}


# 3. Casewise scores

# 3a: casewise scores with respect to Beta (=intercepts and slopes)
#     column order: Y1_int, Y1_x1, Y1_x2, ...| Y2_int, Y2_x1, Y2_x2, ... |
lav_mvreg_scores_beta <- function(Y = NULL,
                                  eXo = NULL,
                                  Beta = NULL,
                                  res.int = NULL,
                                  res.slopes = NULL,
                                  res.cov = NULL,
                                  Sinv.method = "eigen",
                                  res.cov.inv = NULL) {
  Y <- unname(Y)
  Q <- NCOL(Y)
  X <- cbind(1, unname(eXo))
  P <- NCOL(X)

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  # res.cov.inv
  if (is.null(res.cov.inv)) {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # substract Mu
  RES <- Y - X %*% Beta

  # post-multiply with res.cov.inv
  RES <- RES %*% res.cov.inv

  SC.Beta <- X[, rep(1:P, times = Q), drop = FALSE] *
    RES[, rep(1:Q, each = P), drop = FALSE]

  SC.Beta
}


# 3b: casewise scores with respect to vech(res.cov)
lav_mvreg_scores_vech_sigma <- function(Y = NULL,
                                        eXo = NULL,
                                        Beta = NULL,
                                        res.int = NULL,
                                        res.slopes = NULL,
                                        res.cov = NULL,
                                        Sinv.method = "eigen",
                                        res.cov.inv = NULL) {
  Y <- unname(Y)
  Q <- NCOL(Y)
  X <- cbind(1, unname(eXo))

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  # res.cov.inv
  if (is.null(res.cov.inv)) {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # vech(res.cov.inv)
  isigma <- lav_matrix_vech(res.cov.inv)

  # substract X %*% Beta
  RES <- Y - X %*% Beta

  # postmultiply with res.cov.inv
  RES <- RES %*% res.cov.inv

  # tcrossprod
  idx1 <- lav_matrix_vech_col_idx(Q)
  idx2 <- lav_matrix_vech_row_idx(Q)
  Z <- RES[, idx1] * RES[, idx2]

  # substract isigma from each row
  SC <- t(t(Z) - isigma)

  # adjust for vech (and avoiding the 1/2 factor)
  SC[, lav_matrix_diagh_idx(Q)] <- SC[, lav_matrix_diagh_idx(Q)] / 2

  SC
}


# 3c: casewise scores with respect to beta + vech(res.cov)
lav_mvreg_scores_beta_vech_sigma <- function(Y = NULL,
                                             eXo = NULL,
                                             Beta = NULL,
                                             res.int = NULL,
                                             res.slopes = NULL,
                                             res.cov = NULL,
                                             Sinv.method = "eigen",
                                             res.cov.inv = NULL) {
  Y <- unname(Y)
  Q <- NCOL(Y)
  X <- cbind(1, unname(eXo))
  P <- NCOL(X)

  # construct model-implied Beta
  if (is.null(Beta)) {
    Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
  }

  # res.cov.inv
  if (is.null(res.cov.inv)) {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # vech(res.cov.inv)
  isigma <- lav_matrix_vech(res.cov.inv)

  # substract X %*% Beta
  RES <- Y - X %*% Beta

  # postmultiply with res.cov.inv
  RES <- RES %*% res.cov.inv

  SC.Beta <- X[, rep(1:P, times = Q), drop = FALSE] *
    RES[, rep(1:Q, each = P), drop = FALSE]

  # tcrossprod
  idx1 <- lav_matrix_vech_col_idx(Q)
  idx2 <- lav_matrix_vech_row_idx(Q)
  Z <- RES[, idx1] * RES[, idx2]

  # substract isigma from each row
  SC <- t(t(Z) - isigma)

  # adjust for vech (and avoiding the 1/2 factor)
  SC[, lav_matrix_diagh_idx(Q)] <- SC[, lav_matrix_diagh_idx(Q)] / 2

  cbind(SC.Beta, SC)
}

# 4. hessian of logl

# 4a. hessian logl Beta and vech(res.cov) from raw data
lav_mvreg_logl_hessian_data <- function(Y = NULL,
                                        eXo = NULL, # no int
                                        Beta = NULL, # int+slopes
                                        res.int = NULL,
                                        res.slopes = NULL,
                                        res.cov = NULL,
                                        res.cov.inv = NULL,
                                        Sinv.method = "eigen") {
  # sample size
  N <- NROW(Y)

  # observed information
  observed <- lav_mvreg_information_observed_data(
    Y = Y, eXo = eXo,
    Beta = Beta, res.int = res.int, res.slopes = res.slopes,
    res.cov = res.cov, res.cov.inv = res.cov.inv,
    Sinv.method = Sinv.method
  )

  # hessian
  -N * observed
}

# 4b. hessian logl Beta and vech(res.cov) from samplestats
lav_mvreg_logl_hessian_samplestats <- function(sample.res.int = NULL,
                                               sample.res.slopes = NULL,
                                               sample.res.cov = NULL,
                                               sample.mean.x = NULL,
                                               sample.cov.x = NULL,
                                               sample.nobs = NULL,
                                               Beta = NULL, # int + slopes
                                               res.int = NULL, # intercepts only
                                               res.slopes = NULL, # slopes only (y x x)
                                               res.cov = NULL, # res.cov
                                               Sinv.method = "eigen",
                                               res.cov.inv = NULL) {
  # sample size
  N <- sample.nobs

  # information
  observed <- lav_mvreg_information_observed_samplestats(
    sample.res.int = sample.res.int, sample.res.slopes = sample.res.slopes,
    sample.res.cov = sample.res.cov, sample.mean.x = sample.mean.x,
    sample.cov.x = sample.cov.x, Beta = Beta, res.int = res.int,
    res.slopes = res.slopes, res.cov = res.cov, Sinv.method = Sinv.method,
    res.cov.inv = res.cov.inv
  )

  # hessian
  -N * observed
}


# Information h0

# 5a: unit expected information h0 Beta and vech(res.cov)
lav_mvreg_information_expected <- function(Y = NULL, # not used
                                           eXo = NULL, # not used
                                           sample.mean.x = NULL,
                                           sample.cov.x = NULL,
                                           sample.nobs = NULL,
                                           Beta = NULL, # not used
                                           res.int = NULL, # not used
                                           res.slopes = NULL, # not used
                                           res.cov = NULL,
                                           res.cov.inv = NULL,
                                           Sinv.method = "eigen") {
  eXo <- unname(eXo)

  # res.cov.inv
  if (is.null(res.cov.inv)) {
    # invert res.cov
    res.cov.inv <- lav_matrix_symmetric_inverse(
      S = res.cov, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # N
  if (is.null(sample.nobs)) {
    sample.nobs <- nrow(eXo) # hopefully not NULL either
  } else {
    N <- sample.nobs
  }

  # sample.mean.x + sample.cov.x
  if (is.null(sample.mean.x)) {
    sample.mean.x <- base::.colMeans(eXo, m = NROW(eXo), n = NCOL(eXo))
  }
  if (is.null(sample.cov.x)) {
    sample.cov.x <- lav_matrix_cov(eXo)
  }

  # construct sample.xx = 1/N*crossprod(X1) (including intercept)
  sample.xx <- rbind(
    cbind(1, matrix(sample.mean.x, nrow = 1, )),
    cbind(
      matrix(sample.mean.x, ncol = 1),
      sample.cov.x + tcrossprod(sample.mean.x)
    )
  )

  # expected information
  I11 <- res.cov.inv %x% sample.xx
  I22 <- 0.5 * lav_matrix_duplication_pre_post(res.cov.inv %x% res.cov.inv)

  lav_matrix_bdiag(I11, I22)
}

# 5b: unit observed information h0
lav_mvreg_information_observed_data <- function(Y = NULL,
                                                eXo = NULL, # no int
                                                Beta = NULL, # int+slopes
                                                res.int = NULL,
                                                res.slopes = NULL,
                                                res.cov = NULL,
                                                res.cov.inv = NULL,
                                                Sinv.method = "eigen") {
  # create sample statistics
  Y <- unname(Y)
  X1 <- cbind(1, unname(eXo))
  N <- NROW(Y)

  # find 'B'
  QR <- qr(X1)
  sample.B <- qr.coef(QR, Y)

  sample.res.int <- as.numeric(sample.B[1, ])
  sample.res.slopes <- t(sample.B[-1, , drop = FALSE]) # transpose!
  sample.res.cov <- cov(qr.resid(QR, Y)) * (N - 1) / N
  sample.mean.x <- base::.colMeans(eXo, m = NROW(eXo), n = NCOL(eXo))
  sample.cov.x <- lav_matrix_cov(eXo)

  lav_mvreg_information_observed_samplestats(
    sample.res.int = sample.res.int,
    sample.res.slopes = sample.res.slopes, sample.res.cov = sample.res.cov,
    sample.mean.x = sample.mean.x, sample.cov.x = sample.cov.x,
    Beta = Beta, res.int = res.int, res.slopes = res.slopes,
    res.cov = res.cov, Sinv.method = Sinv.method, res.cov.inv = res.cov.inv
  )
}


# 5b-bis: observed information h0 from sample statistics
lav_mvreg_information_observed_samplestats <-
  function(sample.res.int = NULL,
           sample.res.slopes = NULL,
           sample.res.cov = NULL,
           sample.mean.x = NULL,
           sample.cov.x = NULL,
           Beta = NULL, # int + slopes
           res.int = NULL, # intercepts only
           res.slopes = NULL, # slopes only (y x x)
           res.cov = NULL, # res.cov
           Sinv.method = "eigen",
           res.cov.inv = NULL) {
    # construct model-implied Beta
    if (is.null(Beta)) {
      Beta <- rbind(matrix(res.int, nrow = 1), t(res.slopes))
    }

    # construct 'saturated' (sample-based) B
    sample.B <- rbind(matrix(sample.res.int, nrow = 1), t(sample.res.slopes))

    # construct sample.xx = 1/N*crossprod(X1) (including intercept)
    sample.xx <- rbind(
      cbind(1, matrix(sample.mean.x, nrow = 1, )),
      cbind(
        matrix(sample.mean.x, ncol = 1),
        sample.cov.x + tcrossprod(sample.mean.x)
      )
    )

    # W.tilde = S + t(B - Beta) %*% (1/N)*X'X %*% (B - Beta)
    W.tilde <- (sample.res.cov +
      t(sample.B - Beta) %*% sample.xx %*% (sample.B - Beta))

    # res.cov.inv
    if (is.null(res.cov.inv)) {
      # invert res.cov
      res.cov.inv <- lav_matrix_symmetric_inverse(
        S = res.cov, logdet = FALSE,
        Sinv.method = Sinv.method
      )
    }

    H11 <- res.cov.inv %x% sample.xx
    H21 <- lav_matrix_duplication_pre(res.cov.inv %x%
      (res.cov.inv %*% (crossprod(sample.B - Beta, sample.xx))))
    H12 <- t(H21)

    AAA <- res.cov.inv %*% (2 * W.tilde - res.cov) %*% res.cov.inv
    H22 <- (1 / 2) * lav_matrix_duplication_pre_post(res.cov.inv %x% AAA)

    out <- rbind(
      cbind(H11, H12),
      cbind(H21, H22)
    )
    out
  }


# 5c: unit first-order information h0
lav_mvreg_information_firstorder <- function(Y = NULL,
                                             eXo = NULL, # no int
                                             Beta = NULL, # int+slopes
                                             res.int = NULL,
                                             res.slopes = NULL,
                                             res.cov = NULL,
                                             res.cov.inv = NULL,
                                             Sinv.method = "eigen") {
  N <- NROW(Y)

  # scores
  SC <- lav_mvreg_scores_beta_vech_sigma(
    Y = Y, eXo = eXo, Beta = Beta,
    res.int = res.int, res.slopes = res.slopes, res.cov = res.cov,
    Sinv.method = Sinv.method, res.cov.inv = res.cov.inv
  )

  crossprod(SC) / N
}


# 6: inverted information h0

# 6a: inverted unit expected information h0 Beta and vech(res.cov)
#
# lav_mvreg_inverted_information_expected <- function(Y       = NULL, # unused!
# }
