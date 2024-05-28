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
# YR 04 Okt 2018: adding wt= argument, and missing meanstructure=
# YR 27 Jun 2018: adding cluster.idx= argument for information_firstorder
# YR 24 Jul 2022: adding correlation= argument for information_expected
# YR 27 May 2024: adding correlation= argument to most functions

# 0. densities
lav_mvnorm_dmvnorm <- function(Y = NULL,
                               wt = NULL,
                               Mu = NULL,
                               Sigma = NULL,
                               Sigma.inv = NULL,
                               Sinv.method = "eigen",
                               x.idx = integer(0L),
                               x.mean = NULL,
                               x.cov = NULL,
                               log = TRUE) {
  if (is.matrix(Y)) {
    if (is.null(Mu) && is.null(Sigma) && is.null(Sigma.inv)) {
      out <- lav_mvnorm_loglik_data_z(Y = Y, casewise = TRUE)
    } else {
      out <- lav_mvnorm_loglik_data(
        Y = Y, Mu = Mu, Sigma = Sigma,
        casewise = TRUE,
        Sinv.method = Sinv.method
      )
    }
  } else {
    # just one
    P <- length(Y)
    LOG.2PI <- log(2 * pi)

    if (is.null(Mu) && is.null(Sigma) && is.null(Sigma.inv)) {
      # mahalanobis distance
      DIST <- sum(Y * Y)
      out <- -(P * LOG.2PI + DIST) / 2
    } else {
      if (is.null(Sigma.inv)) {
        Sigma.inv <- lav_matrix_symmetric_inverse(
          S = Sigma,
          logdet = TRUE, Sinv.method = Sinv.method
        )
        logdet <- attr(Sigma.inv, "logdet")
      } else {
        logdet <- attr(Sigma.inv, "logdet")
        if (is.null(logdet)) {
          # compute - ln|Sigma.inv|
          ev <- eigen(Sigma.inv, symmetric = TRUE, only.values = TRUE)
          logdet <- -1 * sum(log(ev$values))
        }
      }

      # mahalanobis distance
      Yc <- Y - Mu
      DIST <- sum(Yc %*% Sigma.inv * Yc)
      out <- -(P * LOG.2PI + logdet + DIST) / 2
    }
  }

  if (!is.null(wt)) {
    out <- out * wt
  }

  # x.idx?
  if (length(x.idx) > 0L) {
    if (is.null(Sigma) && is.null(x.cov)) {
      lav_msg_stop(gettext("when x.idx is not empty, we need Sigma or x.cov"))
    }
    if (is.matrix(Y)) {
      X <- Y[, x.idx, drop = FALSE]
    } else {
      X <- Y[x.idx]
    }

    Mu.X <- x.mean
    Sigma.X <- x.cov
    if (is.null(x.mean)) {
      Mu.X <- as.numeric(Mu)[x.idx]
    }
    if (is.null(x.cov)) {
      Sigma.X <- Sigma[x.idx, x.idx, drop = FALSE]
    }

    logl.X <- lav_mvnorm_dmvnorm(
      Y = X, wt = wt, Mu = Mu.X, Sigma = Sigma.X,
      Sigma.inv = NULL,
      Sinv.method = Sinv.method,
      x.idx = integer(0L), log = TRUE
    )

    # subtract logl.X
    out <- out - logl.X
  }

  if (!log) {
    out <- exp(out)
  }

  out
}

# 1. likelihood

# 1a: input is raw data
# (note casewise = TRUE same as: dmvnorm(Y, mean, sigma, log = TRUE))
lav_mvnorm_loglik_data <- function(Y = NULL,
                                   wt = NULL,
                                   Mu = NULL,
                                   Sigma = NULL,
                                   x.idx = integer(0L),
                                   x.mean = NULL,
                                   x.cov = NULL,
                                   casewise = FALSE,
                                   Sinv.method = "eigen") {
  # Y must be a matrix (use lav_mvnorm_dmvnorm() for non-matrix input)
  stopifnot(is.matrix(Y))

  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  P <- NCOL(Y)
  Mu <- as.numeric(Mu)

  if (casewise) {
    LOG.2PI <- log(2 * pi)

    # invert Sigma
    if (Sinv.method == "chol") {
      cS <- chol(Sigma)
      icS <- backsolve(cS, diag(P))
      Yc <- t(t(Y) - Mu)
      DIST <- rowSums((Yc %*% icS)^2)
      logdet <- -2 * sum(log(diag(icS)))
    } else {
      Sigma.inv <- lav_matrix_symmetric_inverse(
        S = Sigma, logdet = TRUE,
        Sinv.method = Sinv.method
      )
      logdet <- attr(Sigma.inv, "logdet")
      # mahalanobis distance
      Yc <- t(t(Y) - Mu)
      DIST <- rowSums(Yc %*% Sigma.inv * Yc)
    }

    loglik <- -(P * LOG.2PI + logdet + DIST) / 2

    # weights
    if (!is.null(wt)) {
      loglik <- loglik * wt
    }
  } else {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = TRUE,
      Sinv.method = Sinv.method
    )
    if (!is.null(wt)) {
      out <- stats::cov.wt(Y, wt = wt, method = "ML")
      sample.mean <- out$center
      sample.cov <- out$cov
    } else {
      sample.mean <- base::.colMeans(Y, m = N, n = P)
      sample.cov <- lav_matrix_cov(Y)
    }
    loglik <- lav_mvnorm_loglik_samplestats(
      sample.mean = sample.mean,
      sample.cov = sample.cov,
      sample.nobs = N,
      Mu = Mu,
      Sigma.inv = Sigma.inv
    )
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    Mu.X <- x.mean
    Sigma.X <- x.cov
    if (is.null(x.mean)) {
      Mu.X <- as.numeric(Mu)[x.idx]
    }
    if (is.null(x.cov)) {
      Sigma.X <- Sigma[x.idx, x.idx, drop = FALSE]
    }
    loglik.x <- lav_mvnorm_loglik_data(
      Y = Y[, x.idx, drop = FALSE],
      wt = wt, Mu = Mu.X, Sigma = Sigma.X,
      x.idx = integer(0L), casewise = casewise,
      Sinv.method = Sinv.method
    )
    # subtract logl.X
    loglik <- loglik - loglik.x
  }

  loglik
}



# 1b: input are sample statistics (mean, cov, N) only
lav_mvnorm_loglik_samplestats <- function(sample.mean = NULL,
                                          sample.cov = NULL,
                                          sample.nobs = NULL,
                                          Mu = NULL,
                                          Sigma = NULL,
                                          x.idx = integer(0L),
                                          x.mean = NULL,
                                          x.cov = NULL,
                                          Sinv.method = "eigen",
                                          Sigma.inv = NULL) {
  P <- length(sample.mean)
  N <- sample.nobs
  Mu <- as.numeric(Mu)
  sample.mean <- as.numeric(sample.mean)
  LOG.2PI <- log(2 * pi)

  if (is.null(Sigma.inv)) {
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = TRUE,
      Sinv.method = Sinv.method
    )
    logdet <- attr(Sigma.inv, "logdet")
  } else {
    logdet <- attr(Sigma.inv, "logdet")
    if (is.null(logdet)) {
      # compute - ln|Sigma.inv|
      ev <- eigen(Sigma.inv, symmetric = TRUE, only.values = TRUE)
      logdet <- -1 * sum(log(ev$values))
    }
  }

  # tr(Sigma^{-1} %*% S)
  DIST1 <- sum(Sigma.inv * sample.cov)
  # (ybar - mu)^T %*% Sigma.inv %*% (ybar - mu)
  Diff <- as.numeric(sample.mean - Mu)
  DIST2 <- sum(as.numeric(crossprod(Diff, Sigma.inv)) * Diff)

  loglik <- -N / 2 * (P * LOG.2PI + logdet + DIST1 + DIST2)

  # fixed.x?
  if (length(x.idx) > 0L) {
    Mu.X <- x.mean
    Sigma.X <- x.cov
    if (is.null(x.mean)) {
      Mu.X <- Mu[x.idx]
    }
    if (is.null(x.cov)) {
      Sigma.X <- Sigma[x.idx, x.idx, drop = FALSE]
    }
    sample.mean.x <- sample.mean[x.idx]
    sample.cov.x <- sample.cov[x.idx, x.idx, drop = FALSE]
    loglik.x <-
      lav_mvnorm_loglik_samplestats(
        sample.mean = sample.mean.x,
        sample.cov = sample.cov.x,
        sample.nobs = sample.nobs,
        Mu = Mu.X, Sigma = Sigma.X,
        x.idx = integer(0L),
        Sinv.method = Sinv.method
      )
    # subtract logl.X
    loglik <- loglik - loglik.x
  }

  loglik
}

# 1c special case: Mu = 0, Sigma = I
lav_mvnorm_loglik_data_z <- function(Y = NULL,
                                     wt = NULL,
                                     casewise = FALSE) {
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  P <- NCOL(Y)
  LOG.2PI <- log(2 * pi)

  if (casewise) {
    DIST <- rowSums(Y * Y)
    loglik <- -(P * LOG.2PI + DIST) / 2
    if (!is.null(wt)) {
      loglik <- loglik * wt
    }
  } else {
    if (!is.null(wt)) {
      out <- stats::cov.wt(Y, wt = wt, method = "ML")
      sample.mean <- out$center
      sample.cov <- out$cov
    } else {
      sample.mean <- base::.colMeans(Y, m = N, n = P)
      sample.cov <- lav_matrix_cov(Y)
    }

    DIST1 <- sum(diag(sample.cov))
    DIST2 <- sum(sample.mean * sample.mean)

    loglik <- -N / 2 * (P * LOG.2PI + DIST1 + DIST2)
  }

  loglik
}





# 2. Derivatives

# 2a: derivative logl with respect to mu
lav_mvnorm_dlogl_dmu <- function(Y = NULL,
                                 wt = NULL,
                                 Mu = NULL,
                                 Sigma = NULL,
                                 x.idx = integer(0L),
                                 Sinv.method = "eigen",
                                 Sigma.inv = NULL) {
  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # substract 'Mu' from Y
  Yc <- t(t(Y) - Mu)

  # weights
  if (!is.null(wt)) {
    Yc <- Yc * wt
  }

  # derivative
  dmu <- as.numeric(Sigma.inv %*% colSums(Yc))

  # fixed.x?
  if (length(x.idx) > 0L) {
    dmu[x.idx] <- 0
  }

  dmu
}

# 2b: derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvnorm_dlogl_dSigma <- function(Y = NULL,
                                    wt = NULL,
                                    Mu = NULL,
                                    Sigma = NULL,
                                    x.idx = integer(0L),
                                    correlation = FALSE,
                                    Sinv.method = "eigen",
                                    Sigma.inv = NULL) {
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # W.tilde
  if (!is.null(wt)) {
    out <- stats::cov.wt(Y, wt = wt, method = "ML")
    SY <- out$cov
    MY <- out$center
    W.tilde <- SY + tcrossprod(MY - Mu)
  } else {
    # substract 'Mu' from Y
    # Yc <- t( t(Y) - Mu )
    # W.tilde <- crossprod(Yc) / N
    W.tilde <- lav_matrix_cov(Y, Mu = Mu)
  }

  # derivative
  dSigma <- -(N / 2) * (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))

  # correlation?
  if (correlation) {
    diag(dSigma) <- 0
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    dSigma[x.idx, x.idx] <- 0
  }

  dSigma
}

# 2c: derivative logl with respect to vech(Sigma)
lav_mvnorm_dlogl_dvechSigma <- function(Y = NULL,
                                        wt = NULL,
                                        Mu = NULL,
                                        Sigma = NULL,
                                        x.idx = integer(0L),
                                        correlation = FALSE,
                                        Sinv.method = "eigen",
                                        Sigma.inv = NULL) {
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # W.tilde
  if (!is.null(wt)) {
    out <- stats::cov.wt(Y, wt = wt, method = "ML")
    SY <- out$cov
    MY <- out$center
    W.tilde <- SY + tcrossprod(MY - Mu)
  } else {
    W.tilde <- lav_matrix_cov(Y, Mu = Mu)
  }

  # derivative (avoiding kronecker product)
  dSigma <- -(N / 2) * (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))

  # fixed.x?
  if (length(x.idx) > 0L) {
    dSigma[x.idx, x.idx] <- 0
  }

  # vech
  if (correlation) {
    dvechSigma <- as.numeric(lav_matrix_duplication_cor_pre(
      as.matrix(lav_matrix_vec(dSigma))
    ))
  } else {
    dvechSigma <- as.numeric(lav_matrix_duplication_pre(
      as.matrix(lav_matrix_vec(dSigma))
    ))
  }

  dvechSigma
}

# 2d: : derivative logl with respect to Mu and vech(Sigma)
lav_mvnorm_dlogl_dmu_dvechSigma <- function(Y = NULL,
                                            wt = NULL,
                                            Mu = NULL,
                                            Sigma = NULL,
                                            x.idx = integer(0L),
                                            correlation = FALSE,
                                            Sinv.method = "eigen",
                                            Sigma.inv = NULL) {
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # substract Mu
  Yc <- t(t(Y) - Mu)

  # W.tilde
  if (!is.null(wt)) {
    out <- stats::cov.wt(Y, wt = wt, method = "ML")
    SY <- out$cov
    MY <- out$center
    W.tilde <- SY + tcrossprod(MY - Mu)
    dmu <- as.numeric(Sigma.inv %*% colSums(Yc * wt))
  } else {
    W.tilde <- lav_matrix_cov(Y, Mu = Mu)
    dmu <- as.numeric(Sigma.inv %*% colSums(Yc))
  }

  # derivative (avoiding kronecker product)
  dSigma <- -(N / 2) * (Sigma.inv - (Sigma.inv %*% W.tilde %*% Sigma.inv))

  # fixed.x?
  if (length(x.idx) > 0L) {
    dSigma[x.idx, x.idx] <- 0
    dmu[x.idx] <- 0
  }

  # vech
  if (correlation) {
    dvechSigma <- as.numeric(lav_matrix_duplication_cor_pre(
      as.matrix(lav_matrix_vec(dSigma))
    ))
  } else {
    dvechSigma <- as.numeric(lav_matrix_duplication_pre(
      as.matrix(lav_matrix_vec(dSigma))
    ))
  }

  c(dmu, dvechSigma)
}

# 3. Casewise scores

# 3a: casewise scores with respect to mu
lav_mvnorm_scores_mu <- function(Y = NULL,
                                 wt = NULL,
                                 Mu = NULL,
                                 x.idx = integer(0L),
                                 Sigma = NULL,
                                 Sinv.method = "eigen",
                                 Sigma.inv = NULL) {
  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # substract Mu
  Yc <- t(t(Y) - Mu)

  # postmultiply with Sigma.inv
  SC <- Yc %*% Sigma.inv

  # weights
  if (!is.null(wt)) {
    SC <- SC * wt
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    SC[, x.idx] <- 0
  }

  SC
}

# 3b: casewise scores with respect to vech(Sigma)
lav_mvnorm_scores_vech_sigma <- function(Y = NULL,
                                         wt = NULL,
                                         Mu = NULL,
                                         Sigma = NULL,
                                         x.idx = integer(0L),
                                         correlation = FALSE,
                                         Sinv.method = "eigen",
                                         Sigma.inv = NULL) {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # vech(Sigma.inv)
  isigma <- lav_matrix_vech(Sigma.inv)

  # substract Mu
  Yc <- t(t(Y) - Mu)

  # postmultiply with Sigma.inv
  Yc <- Yc %*% Sigma.inv

  # tcrossprod
  idx1 <- lav_matrix_vech_col_idx(P)
  idx2 <- lav_matrix_vech_row_idx(P)
  Z <- Yc[, idx1] * Yc[, idx2]

  # substract isigma from each row
  SC <- t(t(Z) - isigma)

  # adjust for vech
  SC[, lav_matrix_diagh_idx(P)] <- SC[, lav_matrix_diagh_idx(P)] / 2

  # fixed.x?
  if (length(x.idx) > 0L) {
    SC[, lav_matrix_vech_which_idx(n = P, idx = x.idx)] <- 0
  }

  # correlation?
  if (correlation) {
    SC <- SC[, -lav_matrix_diagh_idx(n = P), drop = FALSE]
  }

  # weights
  if (!is.null(wt)) {
    SC <- SC * wt
  }

  SC
}

# 3c: casewise scores with respect to mu + vech(Sigma)
lav_mvnorm_scores_mu_vech_sigma <- function(Y = NULL,
                                            wt = NULL,
                                            Mu = NULL,
                                            Sigma = NULL,
                                            x.idx = integer(0L),
                                            correlation = FALSE,
                                            Sinv.method = "eigen",
                                            Sigma.inv = NULL) {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # vech(Sigma.inv)
  isigma <- lav_matrix_vech(Sigma.inv)

  # substract Mu
  Yc <- t(t(Y) - Mu)

  # postmultiply with Sigma.inv
  Yc <- Yc %*% Sigma.inv

  # tcrossprod
  idx1 <- lav_matrix_vech_col_idx(P)
  idx2 <- lav_matrix_vech_row_idx(P)
  Z <- Yc[, idx1] * Yc[, idx2]

  # substract isigma from each row
  SC <- t(t(Z) - isigma)

  # adjust for lav_matrix_duplication_pre (not vech!)
  SC[, lav_matrix_diagh_idx(P)] <- SC[, lav_matrix_diagh_idx(P)] / 2

  # fixed.x?
  if (length(x.idx) > 0L) {
    Yc[, x.idx] <- 0
    SC[, lav_matrix_vech_which_idx(n = P, idx = x.idx)] <- 0
  }

  # correlation?
  if (correlation) {
    SC <- SC[, -lav_matrix_diagh_idx(n = P), drop = FALSE]
  }

  out <- cbind(Yc, SC)

  # weights
  if (!is.null(wt)) {
    out <- out * wt
  }

  out
}


# 4. hessian of logl

# 4a: hessian logl Mu and vech(Sigma) from raw data
lav_mvnorm_logl_hessian_data <- function(Y = NULL,
                                         wt = NULL,
                                         Mu = NULL,
                                         Sigma = NULL,
                                         x.idx = integer(0L),
                                         correlation = FALSE,
                                         Sinv.method = "eigen",
                                         Sigma.inv = NULL,
                                         meanstructure = TRUE) {
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  # observed information
  observed <- lav_mvnorm_information_observed_data(
    Y = Y, wt = wt, Mu = Mu,
    Sigma = Sigma, x.idx = x.idx, correlation = correlation,
    Sinv.method = Sinv.method, Sigma.inv = Sigma.inv,
    meanstructure = meanstructure
  )

  -N * observed
}

# 4b: hessian Mu and vech(Sigma) from samplestats
lav_mvnorm_logl_hessian_samplestats <-
  function(sample.mean = NULL,
           sample.cov = NULL,
           sample.nobs = NULL,
           Mu = NULL,
           Sigma = NULL,
           x.idx = integer(0L),
           correlation = FALSE,
           Sinv.method = "eigen",
           Sigma.inv = NULL,
           meanstructure = TRUE) {
    N <- sample.nobs

    # observed information
    observed <- lav_mvnorm_information_observed_samplestats(
      sample.mean = sample.mean, sample.cov = sample.cov,
      Mu = Mu, Sigma = Sigma,
      x.idx = x.idx, correlation = correlation,
      Sinv.method = Sinv.method, Sigma.inv = Sigma.inv,
      meanstructure = meanstructure
    )

    -N * observed
  }

# 5) Information h0

# 5a: unit expected information h0 Mu and vech(Sigma)
lav_mvnorm_information_expected <- function(Y = NULL, # unused!
                                            wt = NULL, # unused!
                                            Mu = NULL, # unused!
                                            Sigma = NULL,
                                            x.idx = integer(0L),
                                            correlation = FALSE,
                                            Sinv.method = "eigen",
                                            Sigma.inv = NULL,
                                            meanstructure = TRUE) {
  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  if (correlation) {
    I22 <- 0.5 * lav_matrix_duplication_cor_pre_post(Sigma.inv %x%
      Sigma.inv)
  } else {
    I22 <- 0.5 * lav_matrix_duplication_pre_post(Sigma.inv %x% Sigma.inv)
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    pstar.x <- lav_matrix_vech_which_idx(
      n = NCOL(Sigma.inv), idx = x.idx,
      diagonal = !correlation
    )
    I22[pstar.x, ] <- 0
    I22[, pstar.x] <- 0
  }

  if (meanstructure) {
    I11 <- Sigma.inv
    # fixed.x?
    if (length(x.idx) > 0L) {
      I11[x.idx, ] <- 0
      I11[, x.idx] <- 0
    }
    out <- lav_matrix_bdiag(I11, I22)
  } else {
    out <- I22
  }

  out
}

# 5b: unit observed information h0
lav_mvnorm_information_observed_data <- function(Y = NULL,
                                                 wt = NULL,
                                                 Mu = NULL,
                                                 Sigma = NULL,
                                                 x.idx = integer(0L),
                                                 correlation = FALSE,
                                                 Sinv.method = "eigen",
                                                 Sigma.inv = NULL,
                                                 meanstructure = TRUE) {
  if (!is.null(wt)) {
    N <- sum(wt)
    out <- stats::cov.wt(Y, wt = wt, method = "ML")
    sample.cov <- out$cov
    sample.mean <- out$center
  } else {
    N <- NROW(Y)
    # sample statistics
    sample.mean <- colMeans(Y)
    sample.cov <- lav_matrix_cov(Y)
  }

  lav_mvnorm_information_observed_samplestats(
    sample.mean = sample.mean,
    sample.cov = sample.cov, Mu = Mu, Sigma = Sigma,
    x.idx = x.idx, correlation = correlation,
    Sinv.method = Sinv.method, Sigma.inv = Sigma.inv,
    meanstructure = meanstructure
  )
}

# 5b-bis: observed information h0 from sample statistics
lav_mvnorm_information_observed_samplestats <- function(
    sample.mean = NULL,
    sample.cov = NULL,
    Mu = NULL,
    Sigma = NULL,
    x.idx = integer(0L),
    correlation = FALSE,
    Sinv.method = "eigen",
    Sigma.inv = NULL,
    meanstructure = TRUE) {

  sample.mean <- as.numeric(sample.mean)
  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  W.tilde <- sample.cov + tcrossprod(sample.mean - Mu)

  if (meanstructure) {
    I11 <- Sigma.inv
	if (correlation) {
      I21 <- lav_matrix_duplication_cor_pre((Sigma.inv %*%
        (sample.mean - Mu)) %x% Sigma.inv)
	} else {
      I21 <- lav_matrix_duplication_pre((Sigma.inv %*%
        (sample.mean - Mu)) %x% Sigma.inv)
    }
    I12 <- t(I21)
  }

  AAA <- Sigma.inv %*% (2 * W.tilde - Sigma) %*% Sigma.inv
  if (correlation) {
    I22 <- (1 / 2) * lav_matrix_duplication_cor_pre_post(Sigma.inv %x% AAA)
  } else {
    I22 <- (1 / 2) * lav_matrix_duplication_pre_post(Sigma.inv %x% AAA)
  }

  if (meanstructure) {
    out <- rbind(
      cbind(I11, I12),
      cbind(I21, I22)
    )
  } else {
    out <- I22
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    not.x <- lav_matrix_vech_which_idx(
      n = NCOL(Sigma.inv), idx = x.idx,
      diagonal = !correlation,
      add.idx.at.start = meanstructure
    )
    out[, not.x] <- 0
    out[not.x, ] <- 0
  }

  out
}

# 5c: unit first-order information h0
lav_mvnorm_information_firstorder <- function(Y = NULL,
                                              wt = NULL,
                                              cluster.idx = NULL,
                                              Mu = NULL,
                                              Sigma = NULL,
                                              x.idx = integer(0L),
                                              correlation = FALSE,
                                              Sinv.method = "eigen",
                                              Sigma.inv = NULL,
                                              meanstructure = TRUE) {
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  if (meanstructure) {
    SC <- lav_mvnorm_scores_mu_vech_sigma(
      Y = Y, wt = wt,
      Mu = Mu, Sigma = Sigma, x.idx = x.idx, correlation = correlation,
      Sinv.method = Sinv.method, Sigma.inv = Sigma.inv
    )
  } else {
    # the caller should use Mu = sample.mean
    SC <- lav_mvnorm_scores_vech_sigma(
      Y = Y, wt = wt,
      Mu = Mu, Sigma = Sigma, correlation = correlation,
      Sinv.method = Sinv.method, Sigma.inv = Sigma.inv
    )
  }

  # handle clustering
  if (!is.null(cluster.idx)) {
    # take the sum within each cluster
    SC <- rowsum(SC, group = cluster.idx, reorder = FALSE, na.rm = TRUE)

    # lower bias if number of clusters is not very high
    # FIXME: reference?
    nC <- nrow(SC)
    correction.factor <- nC / (nC - 1)
    SC <- SC * sqrt(correction.factor)
  }

  # unit information
  out <- crossprod(SC) / N

  out
}


# 6: inverted information h0

# 6a: inverted unit expected information h0 Mu and vech(Sigma)
#
#     Note: this is the same as lav_samplestats_Gamma_NT()
#           but where COV=Sigma and MEAN=Mu
#
lav_mvnorm_inverted_information_expected <- function(Y = NULL, # unused!
                                                     wt = NULL, # unused!
                                                     Mu = NULL, # unused!
                                                     Sigma = NULL,
                                                     x.idx = integer(0L),
                                                     correlation = FALSE,
                                                     meanstructure = TRUE) {

  if (correlation) {
    Info <- lav_mvnorm_information_expected(Y = Y, wt = wt, Mu = Mu,
      Sigma = Sigma, x.idx = x.idx, correlation = correlation,
      meanstructure = meanstructure)
    out <- solve(Info)
    return(out)
  }

  if (length(x.idx) > 0L) {
    # cov(Y|X) = A - B C^{-1} B'
    # where A = cov(Y), B = cov(Y,X), C = cov(X)
    A <- Sigma[-x.idx, -x.idx, drop = FALSE]
    B <- Sigma[-x.idx,  x.idx, drop = FALSE]
    C <- Sigma[ x.idx,  x.idx, drop = FALSE]
    YbarX <- A - B %*% solve(C, t(B))

    # reinsert YbarX in Y+X (residual) covariance matrix
    YbarX.aug <- matrix(0, nrow = NROW(Sigma), ncol = NCOL(Sigma))
    YbarX.aug[-x.idx, -x.idx] <- YbarX

    # take difference
    R <- Sigma - YbarX.aug

    SS <- 2 * lav_matrix_duplication_ginv_pre_post(Sigma %x% Sigma)
    RR <- 2 * lav_matrix_duplication_ginv_pre_post(R %x% R)
    I22 <- SS - RR

    if (meanstructure) {
      I11 <- YbarX.aug
      out <- lav_matrix_bdiag(I11, I22)
    } else {
      out <- I22
    }
  } else {
    I22 <- 2 * lav_matrix_duplication_ginv_pre_post(Sigma %x% Sigma)
    if (meanstructure) {
      I11 <- Sigma
      out <- lav_matrix_bdiag(I11, I22)
    } else {
      out <- I22
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
