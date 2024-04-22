# the multivariate normal distribution + missing values
#  (so-called 'FIML')

# 1) loglikelihood (from raw data, or sample statitics)
# 2) derivatives with respect to mu, Sigma, vech(Sigma)
# 3) casewise scores with respect to mu, vech(Sigma), mu + vech(Sigma)
# 4) hessian of mu + vech(Sigma)
# 5) (unit) information of mu + vech(Sigma)
#    5a: (unit)    expected information
#    5b: (unit)    observed information
#    5c: (unit) first.order information
#    5d: lav_mvnorm_missing_information_both (both observed + first.order)

# 6) inverted information h0 mu + vech(Sigma)
#    6a: /
#    6b: /
#    6c: /
# 7) ACOV h0 mu + vech(Sigma)
#    7a: 1/N * inverted expected    information
#    7b: 1/N * inverted observed    information
#    7c: 1/N * inverted first-order information
#    7d: sandwich acov

# 10) additional functions
#      - lav_mvnorm_missing_impute_pattern
#      - lav_mvnorm_missing_estep


# YR 09 Feb 2016: first version
# YR 19 Mar 2017: 10)
# YR 03 Okt 2018: a few functions gain a wt= argument
# YR 01 Jul 2018: first_order functions gain cluster.idx= argument


# 1) likelihood

# 1a: input is raw data
#  - two strategies: 1) using missing patterns (pattern = TRUE)
#                    2) truly case per case    (pattern = FALSE)
#    depending on the sample size, missing patterns, etc... one can be
#    (much) faster than the other
lav_mvnorm_missing_loglik_data <- function(Y = NULL,
                                           Mu = NULL,
                                           wt = NULL,
                                           Sigma = NULL,
                                           x.idx = NULL,
                                           casewise = FALSE,
                                           pattern = TRUE,
                                           Sinv.method = "eigen",
                                           log2pi = TRUE,
                                           minus.two = FALSE) {
  if (pattern) {
    llik <- lav_mvnorm_missing_llik_pattern(
      Y = Y, wt = wt, Mu = Mu,
      Sigma = Sigma, x.idx = x.idx, Sinv.method = Sinv.method,
      log2pi = log2pi, minus.two = minus.two
    )
  } else {
    llik <- lav_mvnorm_missing_llik_casewise(
      Y = Y, wt = wt, Mu = Mu,
      Sigma = Sigma, x.idx = x.idx, Sinv.method = Sinv.method,
      log2pi = log2pi, minus.two = minus.two
    )
  }

  if (casewise) {
    loglik <- llik
  } else {
    loglik <- sum(llik, na.rm = TRUE)
  }

  loglik
}

# 1b: input are sample statistics (mean, cov, N) per pattern
lav_mvnorm_missing_loglik_samplestats <- function(Yp = NULL,
                                                  Mu = NULL,
                                                  Sigma = NULL,
                                                  x.idx = NULL,
                                                  x.mean = NULL,
                                                  x.cov = NULL,
                                                  Sinv.method = "eigen",
                                                  log2pi = TRUE,
                                                  minus.two = FALSE) {
  # if(!is.null(x.idx) && length(x.idx) > 0L) {
  #    #warning("lavaan WARNING: x.idx not supported yet (ignored)")
  # }

  LOG.2PI <- log(2 * pi)
  pat.N <- length(Yp)
  P <- length(Yp[[1]]$var.idx)

  # global inverse + logdet
  Sigma.inv <- lav_matrix_symmetric_inverse(
    S = Sigma, logdet = TRUE,
    Sinv.method = Sinv.method
  )
  Sigma.logdet <- attr(Sigma.inv, "logdet")

  # DIST/logdet per pattern
  DIST <- logdet <- P.LOG.2PI <- numeric(pat.N)

  # for each pattern, compute sigma.inv/logdet; compute DIST for all
  # observations of this pattern
  for (p in seq_len(pat.N)) {
    # observed variables for this pattern
    var.idx <- Yp[[p]]$var.idx

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # constant
    P.LOG.2PI[p] <- sum(var.idx) * LOG.2PI * Yp[[p]]$freq

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = TRUE, S.logdet = Sigma.logdet
      )
      logdet[p] <- attr(sigma.inv, "logdet") * Yp[[p]]$freq
    } else {
      sigma.inv <- Sigma.inv
      logdet[p] <- Sigma.logdet * Yp[[p]]$freq
    }

    TT <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])
    DIST[p] <- sum(sigma.inv * TT) * Yp[[p]]$freq
  }

  # loglikelihood all data
  if (log2pi) {
    loglik <- sum(-(P.LOG.2PI + logdet + DIST) / 2)
  } else {
    loglik <- sum(-(logdet + DIST) / 2)
  }

  if (minus.two) {
    loglik <- -2 * loglik
  }

  # x.idx
  if (length(x.idx) > 0L) {
    stopifnot(!is.null(x.cov))
    # Note: x.cov should be identical to Sigma[x.idx, x.idx]
    #       so we don't really need x.cov
    N <- sum(sapply(Yp, "[[", "freq"))
    loglik.x <- lav_mvnorm_h1_loglik_samplestats(
      sample.cov = x.cov,
      sample.nobs = N
    )

    loglik <- loglik - loglik.x
  }

  loglik
}

## casewise loglikelihoods

# casewise Sinv.method
lav_mvnorm_missing_llik_casewise <- function(Y = NULL,
                                             wt = NULL,
                                             Mu = NULL,
                                             Sigma = NULL,
                                             x.idx = NULL,
                                             Sinv.method = "eigen",
                                             log2pi = TRUE,
                                             minus.two = FALSE) {
  P <- NCOL(Y)
  LOG.2PI <- log(2 * pi)
  Mu <- as.numeric(Mu)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  NY <- NROW(Y)

  # global inverse + logdet
  Sigma.inv <- lav_matrix_symmetric_inverse(
    S = Sigma, logdet = TRUE,
    Sinv.method = Sinv.method
  )
  Sigma.logdet <- attr(Sigma.inv, "logdet")

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # DIST/logdet per case
  DIST <- logdet <- P.LOG.2PI <- rep(as.numeric(NA), NY)

  # missing pattern per case
  OBS <- !is.na(Y)
  P.i <- rowSums(OBS)

  # constant
  P.LOG.2PI <- P.i * LOG.2PI

  # complete cases first (only an advantage if we have mostly complete
  # observations)
  other.idx <- seq_len(NY)
  complete.idx <- which(P.i == P)
  if (length(complete.idx) > 0L) {
    other.idx <- other.idx[-complete.idx]
    DIST[complete.idx] <-
      rowSums(Yc[complete.idx, , drop = FALSE] %*% Sigma.inv *
        Yc[complete.idx, , drop = FALSE])
    logdet[complete.idx] <- Sigma.logdet
  }

  # non-complete cases
  for (i in other.idx) {
    na.idx <- which(!OBS[i, ])

    # catch empty cases
    if (length(na.idx) == P) next

    # invert Sigma for this pattern
    sigma.inv <- lav_matrix_symmetric_inverse_update(
      S.inv = Sigma.inv,
      rm.idx = na.idx, logdet = TRUE, S.logdet = Sigma.logdet
    )
    logdet[i] <- attr(sigma.inv, "logdet")

    # distance for this case
    DIST[i] <- sum(sigma.inv * crossprod(Yc[i, OBS[i, ], drop = FALSE]))
  }

  # compute casewise loglikelihoods
  if (log2pi) {
    llik <- -(P.LOG.2PI + logdet + DIST) / 2
  } else {
    llik <- -(logdet + DIST) / 2
  }

  # minus.two
  if (minus.two) {
    llik <- -2 * llik
  }

  # weights?
  if (!is.null(wt)) {
    llik <- llik * wt
  }

  # x.idx
  if (length(x.idx) > 0L) {
    llik.x <- lav_mvnorm_missing_llik_casewise(
      Y = Y[, x.idx, drop = FALSE],
      wt = wt, Mu = Mu[x.idx],
      Sigma = Sigma[x.idx, x.idx, drop = FALSE],
      x.idx = NULL, Sinv.method = Sinv.method,
      log2pi = log2pi, minus.two = minus.two
    )
    llik <- llik - llik.x
  }

  llik
}

# pattern-based, but casewise loglikelihoods
lav_mvnorm_missing_llik_pattern <- function(Y = NULL,
                                            Mp = NULL,
                                            wt = NULL,
                                            Mu = NULL,
                                            Sigma = NULL,
                                            x.idx = NULL,
                                            Sinv.method = "eigen",
                                            log2pi = TRUE,
                                            minus.two = FALSE) {
  P <- NCOL(Y)
  LOG.2PI <- log(2 * pi)
  Mu <- as.numeric(Mu)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  NY <- NROW(Y)

  # global inverse + logdet
  Sigma.inv <- lav_matrix_symmetric_inverse(
    S = Sigma, logdet = TRUE,
    Sinv.method = Sinv.method
  )
  Sigma.logdet <- attr(Sigma.inv, "logdet")

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # DIST/logdet per case
  DIST <- logdet <- P.LOG.2PI <- rep(as.numeric(NA), NY)

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # for each pattern, compute sigma.inv/logdet; compute DIST for all
  # observations of this pattern
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # identify cases with this pattern
    case.idx <- Mp$case.idx[[p]]

    # constant
    P.LOG.2PI[case.idx] <- sum(var.idx) * LOG.2PI

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = TRUE, S.logdet = Sigma.logdet
      )
      logdet[case.idx] <- attr(sigma.inv, "logdet")
    } else {
      sigma.inv <- Sigma.inv
      logdet[case.idx] <- Sigma.logdet
    }

    if (Mp$freq[p] == 1L) {
      DIST[case.idx] <- sum(sigma.inv *
        crossprod(Yc[case.idx, var.idx, drop = FALSE]))
    } else {
      DIST[case.idx] <-
        rowSums(Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv *
          Yc[case.idx, var.idx, drop = FALSE])
    }
  }

  # compute casewise loglikelihoods
  if (log2pi) {
    llik <- -(P.LOG.2PI + logdet + DIST) / 2
  } else {
    llik <- -(logdet + DIST) / 2
  }

  # minus.two
  if (minus.two) {
    llik <- -2 * llik
  }

  # weights?
  if (!is.null(wt)) {
    llik <- llik * wt
  }

  # x.idx -- using casewise (as patterns for Y may not be the same as
  #                          patterns for Y[,-x.idx])
  if (length(x.idx) > 0L) {
    llik.x <- lav_mvnorm_missing_llik_casewise(
      Y = Y[, x.idx, drop = FALSE],
      wt = wt, Mu = Mu[x.idx],
      Sigma = Sigma[x.idx, x.idx, drop = FALSE],
      x.idx = NULL, Sinv.method = Sinv.method,
      log2pi = log2pi, minus.two = minus.two
    )
    llik <- llik - llik.x
  }

  llik
}




# 2. Derivatives

# 2a: derivative logl with respect to mu
lav_mvnorm_missing_dlogl_dmu <- function(Y = NULL,
                                         wt = NULL,
                                         Mu = NULL,
                                         Sigma = NULL,
                                         x.idx = NULL,
                                         Sigma.inv = NULL,
                                         Sinv.method = "eigen") {
  SC <- lav_mvnorm_missing_scores_mu(
    Y = Y, wt = wt, Mu = Mu, Sigma = Sigma,
    x.idx = x.idx, Sigma.inv = Sigma.inv, Sinv.method = Sinv.method
  )

  colSums(SC, na.rm = TRUE)
}

# 2abis: using samplestats
lav_mvnorm_missing_dlogl_dmu_samplestats <- function(Yp = NULL,
                                                     Mu = NULL,
                                                     Sigma = NULL,
                                                     x.idx = NULL,
                                                     Sigma.inv = NULL,
                                                     Sinv.method = "eigen") {
  pat.N <- length(Yp)
  P <- length(Yp[[1]]$var.idx)

  if (is.null(Sigma.inv)) {
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # dmu
  dmu <- numeric(P)

  # for each pattern, compute sigma.inv
  for (p in seq_len(pat.N)) {
    # observed variables for this pattern
    var.idx <- Yp[[p]]$var.idx

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
    } else {
      sigma.inv <- Sigma.inv
    }

    # dmu for this pattern
    dmu.pattern <- as.numeric(sigma.inv %*% (Yp[[p]]$MY - Mu[var.idx]))

    # update mu
    dmu[var.idx] <- dmu[var.idx] + (dmu.pattern * Yp[[p]]$freq)
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    dmu[x.idx] <- 0
  }

  dmu
}



# 2b: derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvnorm_missing_dlogl_dSigma <- function(Y = NULL,
                                            Mp = NULL,
                                            wt = NULL,
                                            Mu = NULL,
                                            Sigma = NULL,
                                            x.idx = NULL,
                                            Sigma.inv = NULL,
                                            Sinv.method = "eigen") {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  NY <- NROW(Y)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # dvechSigma
  dSigma <- matrix(0, P, P)

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # for each pattern
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # cases with this pattern
    case.idx <- Mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
    } else {
      sigma.inv <- Sigma.inv
    }

    if (!is.null(wt)) {
      FREQ <- sum(wt[case.idx])
    } else {
      FREQ <- Mp$freq[p]
    }

    if (length(case.idx) > 1L) {
      if (!is.null(wt)) {
        out <- stats::cov.wt(Y[case.idx, var.idx, drop = FALSE],
          wt = wt[Mp$case.idx[[p]]], method = "ML"
        )
        SY <- out$cov
        MY <- out$center
        W.tilde <- SY + tcrossprod(MY - Mu[var.idx])
      } else {
        W.tilde <- crossprod(Yc[case.idx, var.idx, drop = FALSE]) / FREQ
      }
    } else {
      W.tilde <- tcrossprod(Yc[case.idx, var.idx])
    }

    # dSigma for this pattern
    dSigma.pattern <- matrix(0, P, P)
    dSigma.pattern[var.idx, var.idx] <- -(1 / 2) * (sigma.inv -
      (sigma.inv %*% W.tilde %*% sigma.inv))

    # update dSigma
    dSigma <- dSigma + (dSigma.pattern * FREQ)
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    dSigma[x.idx, x.idx] <- 0
  }

  dSigma
}

# 2bbis: using samplestats
lav_mvnorm_missing_dlogl_dSigma_samplestats <- function(Yp = NULL,
                                                        Mu = NULL,
                                                        Sigma = NULL,
                                                        x.idx = NULL,
                                                        Sigma.inv = NULL,
                                                        Sinv.method = "eigen") {
  pat.N <- length(Yp)
  P <- length(Yp[[1]]$var.idx)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # dvechSigma
  dSigma <- matrix(0, P, P)

  # for each pattern
  for (p in seq_len(pat.N)) {
    # observed variables for this pattern
    var.idx <- Yp[[p]]$var.idx

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
    } else {
      sigma.inv <- Sigma.inv
    }

    W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])

    # dSigma for this pattern
    dSigma.pattern <- matrix(0, P, P)
    dSigma.pattern[var.idx, var.idx] <- -(1 / 2) * (sigma.inv -
      (sigma.inv %*% W.tilde %*% sigma.inv))

    # update dSigma
    dSigma <- dSigma + (dSigma.pattern * Yp[[p]]$freq)
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    dSigma[x.idx, x.idx] <- 0
  }

  dSigma
}


# 2c: derivative logl with respect to vech(Sigma)
lav_mvnorm_missing_dlogl_dvechSigma <- function(Y = NULL,
                                                wt = NULL,
                                                Mu = NULL,
                                                x.idx = NULL,
                                                Sigma = NULL,
                                                Sigma.inv = NULL,
                                                Sinv.method = "eigen") {
  dSigma <- lav_mvnorm_missing_dlogl_dSigma(
    Y = Y, wt = wt, Mu = Mu,
    Sigma = Sigma, x.idx = x.idx, Sigma.inv = Sigma.inv,
    Sinv.method = Sinv.method
  )

  dvechSigma <- as.numeric(lav_matrix_duplication_pre(
    as.matrix(lav_matrix_vec(dSigma))
  ))

  dvechSigma
}

# 2cbis: using samplestats
lav_mvnorm_missing_dlogl_dvechSigma_samplestats <-
  function(Yp = NULL,
           Mu = NULL,
           Sigma = NULL,
           x.idx = NULL,
           Sigma.inv = NULL,
           Sinv.method = "eigen") {
    pat.N <- length(Yp)
    P <- length(Yp[[1]]$var.idx)

    if (is.null(Sigma.inv)) {
      # invert Sigma
      Sigma.inv <- lav_matrix_symmetric_inverse(
        S = Sigma, logdet = FALSE,
        Sinv.method = Sinv.method
      )
    }

    # dvechSigma
    dvechSigma <- numeric(P * (P + 1) / 2)

    # for each pattern
    for (p in seq_len(pat.N)) {
      # observed variables for this pattern
      var.idx <- Yp[[p]]$var.idx

      # missing values for this pattern
      na.idx <- which(!var.idx)

      # invert Sigma for this pattern
      if (length(na.idx) > 0L) {
        sigma.inv <- lav_matrix_symmetric_inverse_update(
          S.inv = Sigma.inv,
          rm.idx = na.idx, logdet = FALSE
        )
      } else {
        sigma.inv <- Sigma.inv
      }

      W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])

      # dSigma for this pattern
      dSigma.pattern <- matrix(0, P, P)
      dSigma.pattern[var.idx, var.idx] <- -(1 / 2) * (sigma.inv -
        (sigma.inv %*% W.tilde %*% sigma.inv))

      # fixed.x?
      if (length(x.idx) > 0L) {
        dSigma.pattern[x.idx, x.idx] <- 0
      }

      # convert to vechSigma
      dvechSigma.pattern <- as.numeric(lav_matrix_duplication_pre(
        as.matrix(lav_matrix_vec(dSigma.pattern))
      ))

      # update dvechSigma
      dvechSigma <- dvechSigma + (dvechSigma.pattern * Yp[[p]]$freq)
    }

    dvechSigma
  }




# 3. Casewise scores

# 3a: casewise scores with respect to mu
lav_mvnorm_missing_scores_mu <- function(Y = NULL,
                                         wt = NULL,
                                         Mp = NULL,
                                         Mu = NULL,
                                         Sigma = NULL,
                                         x.idx = NULL,
                                         Sigma.inv = NULL,
                                         Sinv.method = "eigen") {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  NY <- NROW(Y)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # dmu per case
  dmu <- matrix(as.numeric(NA), NY, P)

  # for each pattern, compute sigma.inv
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # missing values for this pattern
    na.idx <- which(!var.idx)

    case.idx <- Mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
    } else {
      sigma.inv <- Sigma.inv
    }

    # compute dMu for all observations of this pattern
    dmu[case.idx, var.idx] <-
      Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv
  }

  # weights
  if (!is.null(wt)) {
    dmu <- dmu * wt
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    dmu[, x.idx] <- 0
  }

  dmu
}

# 3b: casewise scores with respect to vech(Sigma)
lav_mvnorm_missing_scores_vech_sigma <- function(Y = NULL,
                                                 wt = NULL,
                                                 Mp = NULL,
                                                 Mu = NULL,
                                                 Sigma = NULL,
                                                 x.idx = NULL,
                                                 Sigma.inv = NULL,
                                                 Sinv.method = "eigen") {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  NY <- NROW(Y)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # for the tcrossprod
  idx1 <- lav_matrix_vech_col_idx(P)
  idx2 <- lav_matrix_vech_row_idx(P)

  # vech(Sigma.inv)
  iSigma <- lav_matrix_vech(Sigma.inv)

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # SC
  SC <- matrix(as.numeric(NA), nrow = NY, ncol = length(iSigma))

  # for each pattern
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # cases with this pattern
    case.idx <- Mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
      tmp <- matrix(0, P, P)
      tmp[var.idx, var.idx] <- sigma.inv
      isigma <- lav_matrix_vech(tmp)
    } else {
      sigma.inv <- Sigma.inv
      isigma <- iSigma
    }

    # postmultiply these cases with sigma.inv
    Yc[case.idx, var.idx] <- Yc[case.idx, var.idx] %*% sigma.inv

    # tcrossprod
    SC[case.idx, ] <- Yc[case.idx, idx1] * Yc[case.idx, idx2]

    # substract isigma from each row
    SC[case.idx, ] <- t(t(SC[case.idx, , drop = FALSE]) - isigma)
  }

  # adjust for vech
  SC[, lav_matrix_diagh_idx(P)] <- SC[, lav_matrix_diagh_idx(P)] / 2

  # weights
  if (!is.null(wt)) {
    SC <- SC * wt
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    not.x <- eliminate.pstar.idx(P, el.idx = x.idx)
    SC[, !not.x] <- 0
  }

  SC
}

# 3c: casewise scores with respect to mu + vech(Sigma)
lav_mvnorm_missing_scores_mu_vech_sigma <- function(Y = NULL,
                                                    Mp = NULL,
                                                    wt = NULL,
                                                    Mu = NULL,
                                                    Sigma = NULL,
                                                    x.idx = NULL,
                                                    Sigma.inv = NULL,
                                                    Sinv.method = "eigen") {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  NY <- NROW(Y)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # for the tcrossprod
  idx1 <- lav_matrix_vech_col_idx(P)
  idx2 <- lav_matrix_vech_row_idx(P)

  # vech(Sigma.inv)
  iSigma <- lav_matrix_vech(Sigma.inv)

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # dmu per case
  dmu <- matrix(as.numeric(NA), NY, P)

  # SC
  SC <- matrix(as.numeric(NA), nrow = NY, ncol = length(iSigma))

  # for each pattern, compute Yc %*% sigma.inv
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # cases with this pattern
    case.idx <- Mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
      tmp <- matrix(0, P, P)
      tmp[var.idx, var.idx] <- sigma.inv
      isigma <- lav_matrix_vech(tmp)
    } else {
      sigma.inv <- Sigma.inv
      isigma <- iSigma
    }

    # compute dMu for all observations of this pattern
    dmu[case.idx, var.idx] <-
      Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv

    # postmultiply these cases with sigma.inv
    Yc[case.idx, var.idx] <- Yc[case.idx, var.idx] %*% sigma.inv

    # tcrossprod
    SC[case.idx, ] <- Yc[case.idx, idx1] * Yc[case.idx, idx2]

    # substract isigma from each row
    SC[case.idx, ] <- t(t(SC[case.idx, , drop = FALSE]) - isigma)
  }

  # adjust for vech
  SC[, lav_matrix_diagh_idx(P)] <- SC[, lav_matrix_diagh_idx(P)] / 2

  out <- cbind(dmu, SC)

  # weights
  if (!is.null(wt)) {
    out <- out * wt
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    not.x <- eliminate.pstar.idx(P, el.idx = x.idx, meanstructure = TRUE)
    out[, !not.x] <- 0
  }

  out
}


# 4) Hessian of logl
lav_mvnorm_missing_logl_hessian_data <- function(Y = NULL,
                                                 Mp = NULL,
                                                 wt = NULL,
                                                 Mu = NULL,
                                                 Sigma = NULL,
                                                 x.idx = NULL,
                                                 Sinv.method = "eigen",
                                                 Sigma.inv = NULL) {
  # missing patterns
  Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp, wt = wt)

  lav_mvnorm_missing_logl_hessian_samplestats(
    Yp = Yp, Mu = Mu,
    Sigma = Sigma, x.idx = x.idx, Sinv.method = Sinv.method,
    Sigma.inv = Sigma.inv
  )
}

lav_mvnorm_missing_logl_hessian_samplestats <-
  function(Yp = NULL,
           # wt not needed
           Mu = NULL,
           Sigma = NULL,
           x.idx = NULL,
           Sinv.method = "eigen",
           Sigma.inv = NULL) {
    pat.N <- length(Yp)
    P <- length(Yp[[1]]$var.idx)

    if (is.null(Sigma.inv)) {
      Sigma.inv <- lav_matrix_symmetric_inverse(
        S = Sigma, logdet = FALSE,
        Sinv.method = Sinv.method
      )
    }

    H11 <- matrix(0, P, P)
    H21 <- matrix(0, P * (P + 1) / 2, P)
    H22 <- matrix(0, P * (P + 1) / 2, P * (P + 1) / 2)

    # for each pattern, compute sigma.inv
    for (p in seq_len(pat.N)) {
      # observed variables
      var.idx <- Yp[[p]]$var.idx
      pat.freq <- Yp[[p]]$freq

      # missing values for this pattern
      na.idx <- which(!var.idx)

      # invert Sigma for this pattern
      if (length(na.idx) > 0L) {
        sigma.inv <- lav_matrix_symmetric_inverse_update(
          S.inv = Sigma.inv,
          rm.idx = na.idx, logdet = FALSE
        )
      } else {
        sigma.inv <- Sigma.inv
      }

      S.inv <- matrix(0, P, P)
      S.inv[var.idx, var.idx] <- sigma.inv

      tmp21 <- matrix(0, P, 1)
      tmp21[var.idx, 1] <- sigma.inv %*% (Yp[[p]]$MY - Mu[var.idx])

      W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])
      AAA <- (sigma.inv %*%
        (2 * W.tilde - Sigma[var.idx, var.idx, drop = FALSE]) %*%
        sigma.inv)
      tmp22 <- matrix(0, P, P)
      tmp22[var.idx, var.idx] <- AAA

      i11 <- S.inv
      i21 <- lav_matrix_duplication_pre(tmp21 %x% S.inv)
      i22 <- (1 / 2) * lav_matrix_duplication_pre_post(S.inv %x% tmp22)

      H11 <- H11 + pat.freq * i11
      H21 <- H21 + pat.freq * i21
      H22 <- H22 + pat.freq * i22
    }

    H12 <- t(H21)

    out <- -1 * rbind(
      cbind(H11, H12),
      cbind(H21, H22)
    )

    # fixed.x?
    if (length(x.idx) > 0L) {
      not.x <- eliminate.pstar.idx(
        nvar = P, el.idx = x.idx,
        meanstructure = TRUE
      )
      out[, !not.x] <- 0
      out[!not.x, ] <- 0
    }

    out
  }




# 5) Information

# 5a: expected unit information Mu and vech(Sigma)
#     (only useful under MCAR)
# (old term: Abeta, expected)
lav_mvnorm_missing_information_expected <- function(Y = NULL,
                                                    Mp = NULL,
                                                    wt = NULL,
                                                    Mu = NULL, # unused
                                                    Sigma = NULL,
                                                    x.idx = NULL,
                                                    Sigma.inv = NULL,
                                                    Sinv.method = "eigen") {
  P <- NCOL(Y)

  if (is.null(Sigma.inv)) {
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # N
  if (!is.null(wt)) {
    if (length(Mp$empty.idx) > 0L) {
      N <- sum(wt) - sum(wt[Mp$empty.idx])
    } else {
      N <- sum(wt)
    }
  } else {
    N <- sum(Mp$freq) # removed empty cases!
  }

  I11 <- matrix(0, P, P)
  I22 <- matrix(0, P * (P + 1) / 2, P * (P + 1) / 2)

  # for each pattern, compute sigma.inv
  for (p in seq_len(Mp$npatterns)) {
    # observed variables
    var.idx <- Mp$pat[p, ]

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
    } else {
      sigma.inv <- Sigma.inv
    }

    S.inv <- matrix(0, P, P)
    S.inv[var.idx, var.idx] <- sigma.inv

    S2.inv <- 0.5 * lav_matrix_duplication_pre_post(S.inv %x% S.inv)

    if (!is.null(wt)) {
      FREQ <- sum(wt[Mp$case.idx[[p]]])
    } else {
      FREQ <- Mp$freq[p]
    }

    I11 <- I11 + FREQ * S.inv
    I22 <- I22 + FREQ * S2.inv
  }

  out <- lav_matrix_bdiag(I11, I22) / N

  # fixed.x?
  if (length(x.idx) > 0L) {
    not.x <- eliminate.pstar.idx(
      nvar = P, el.idx = x.idx,
      meanstructure = TRUE
    )
    out[!not.x, ] <- 0
    out[, !not.x] <- 0
  }

  out
}

# 5b: unit observed information Mu and vech(Sigma) from raw data
# (old term: Abeta, observed)
lav_mvnorm_missing_information_observed_data <- function(Y = NULL,
                                                         Mp = NULL,
                                                         wt = NULL,
                                                         Mu = NULL,
                                                         Sigma = NULL,
                                                         x.idx = NULL,
                                                         Sinv.method = "eigen",
                                                         Sigma.inv = NULL) {
  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # N
  if (!is.null(wt)) {
    if (length(Mp$empty.idx) > 0L) {
      N <- sum(wt) - sum(wt[Mp$empty.idx])
    } else {
      N <- sum(wt)
    }
  } else {
    N <- sum(Mp$freq) # removed empty cases!
  }

  # observed information
  observed <- lav_mvnorm_missing_logl_hessian_data(
    Y = Y, Mp = Mp, wt = wt,
    Mu = Mu, Sigma = Sigma, x.idx = x.idx,
    Sinv.method = Sinv.method, Sigma.inv = Sigma.inv
  )

  -observed / N
}

# 5b-bis: unit observed information Mu and vech(Sigma) from samplestats
lav_mvnorm_missing_information_observed_samplestats <-
  function(Yp = NULL,
           # wt not needed
           Mu = NULL,
           Sigma = NULL,
           x.idx = NULL,
           Sinv.method = "eigen",
           Sigma.inv = NULL) {
    N <- sum(sapply(Yp, "[[", "freq")) # implicitly: removed empty cases!

    # observed information
    observed <- lav_mvnorm_missing_logl_hessian_samplestats(
      Yp = Yp, Mu = Mu,
      Sigma = Sigma, x.idx = x.idx, Sinv.method = Sinv.method,
      Sigma.inv = Sigma.inv
    )

    -observed / N
  }

# 5c: unit first-order information Mu and vech(Sigma) from raw data
# (old term: Bbeta)
lav_mvnorm_missing_information_firstorder <- function(Y = NULL,
                                                      Mp = NULL,
                                                      wt = NULL,
                                                      cluster.idx = NULL,
                                                      Mu = NULL,
                                                      Sigma = NULL,
                                                      x.idx = NULL,
                                                      Sinv.method = "eigen",
                                                      Sigma.inv = NULL) {
  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # N
  if (!is.null(wt)) {
    if (length(Mp$empty.idx) > 0L) {
      N <- sum(wt) - sum(wt[Mp$empty.idx])
    } else {
      N <- sum(wt)
    }
  } else {
    N <- sum(Mp$freq) # removed empty cases!
  }

  SC <- lav_mvnorm_missing_scores_mu_vech_sigma(
    Y = Y, Mp = Mp, wt = wt,
    Mu = Mu, Sigma = Sigma, x.idx = x.idx, Sinv.method = Sinv.method,
    Sigma.inv = Sigma.inv
  )

  # handle clustering
  if (!is.null(cluster.idx)) {
    # take the sum within each cluster
    SC <- rowsum(SC, group = cluster.idx, reorder = FALSE, na.rm = TRUE)

    # lower bias is number of clusters is not very high
    nC <- nrow(SC)
    correction.factor <- nC / (nC - 1)
    SC <- SC * sqrt(correction.factor)
  }

  lav_matrix_crossprod(SC) / N
}

# 5d: both unit first-order information and expected/observed information
#     from raw data, in one go for efficiency
lav_mvnorm_missing_information_both <- function(Y = NULL,
                                                Mp = NULL,
                                                wt = NULL,
                                                cluster.idx = NULL,
                                                Mu = NULL,
                                                Sigma = NULL,
                                                x.idx = NULL,
                                                Sinv.method = "eigen",
                                                Sigma.inv = NULL,
                                                information = "observed") {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # for the tcrossprod
  idx1 <- lav_matrix_vech_col_idx(P)
  idx2 <- lav_matrix_vech_row_idx(P)

  # vech(Sigma.inv)
  iSigma <- lav_matrix_vech(Sigma.inv)

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  if (information == "observed") {
    Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp, wt = wt)
  }

  # N
  if (!is.null(wt)) {
    if (length(Mp$empty.idx) > 0L) {
      N <- sum(wt) - sum(wt[Mp$empty.idx])
    } else {
      N <- sum(wt)
    }
  } else {
    N <- sum(Mp$freq) # removed empty cases!
  }

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # dmu per case
  dmu <- matrix(as.numeric(NA), nrow = NROW(Y), ncol = P)

  # SC
  SC <- matrix(as.numeric(NA), nrow = NROW(Y), ncol = length(iSigma))

  # expected/observed information
  I11 <- matrix(0, P, P)
  I22 <- matrix(0, P * (P + 1) / 2, P * (P + 1) / 2)
  if (information == "observed") {
    I21 <- matrix(0, P * (P + 1) / 2, P)
  }

  # for each pattern, compute Yc %*% sigma.inv
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # cases with this pattern
    case.idx <- Mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na.idx) > 0L) {
      sigma.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = Sigma.inv,
        rm.idx = na.idx, logdet = FALSE
      )
      tmp <- matrix(0, P, P)
      tmp[var.idx, var.idx] <- sigma.inv
      isigma <- lav_matrix_vech(tmp)
    } else {
      sigma.inv <- Sigma.inv
      isigma <- iSigma
    }

    # information
    S.inv <- matrix(0, P, P)
    S.inv[var.idx, var.idx] <- sigma.inv

    if (!is.null(wt)) {
      FREQ <- sum(wt[case.idx])
    } else {
      FREQ <- Mp$freq[p]
    }

    if (information == "expected") {
      S2.inv <- 0.5 * lav_matrix_duplication_pre_post(S.inv %x% S.inv)

      I11 <- I11 + FREQ * S.inv
      I22 <- I22 + FREQ * S2.inv
    } else {
      pat.freq <- Yp[[p]]$freq

      tmp21 <- matrix(0, P, 1)
      tmp21[var.idx, 1] <- sigma.inv %*% (Yp[[p]]$MY - Mu[var.idx])

      W.tilde <- Yp[[p]]$SY + tcrossprod(Yp[[p]]$MY - Mu[var.idx])
      AAA <- (sigma.inv %*%
        (2 * W.tilde - Sigma[var.idx, var.idx, drop = FALSE]) %*%
        sigma.inv)
      tmp22 <- matrix(0, P, P)
      tmp22[var.idx, var.idx] <- AAA

      i11 <- S.inv
      i21 <- lav_matrix_duplication_pre(tmp21 %x% S.inv)
      i22 <- (1 / 2) * lav_matrix_duplication_pre_post(S.inv %x% tmp22)

      I11 <- I11 + pat.freq * i11
      I21 <- I21 + pat.freq * i21
      I22 <- I22 + pat.freq * i22
    }

    # compute dMu for all observations of this pattern
    dmu[case.idx, var.idx] <-
      Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv

    # postmultiply these cases with sigma.inv
    Yc[case.idx, var.idx] <- Yc[case.idx, var.idx] %*% sigma.inv

    # tcrossprod
    SC[case.idx, ] <- Yc[case.idx, idx1] * Yc[case.idx, idx2]

    # substract isigma from each row
    SC[case.idx, ] <- t(t(SC[case.idx, , drop = FALSE]) - isigma)
  }

  # adjust for vech
  SC[, lav_matrix_diagh_idx(P)] <- SC[, lav_matrix_diagh_idx(P)] / 2

  # add dmu
  SC <- cbind(dmu, SC)

  # weights
  if (!is.null(wt)) {
    SC <- SC * wt
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    not.x <- eliminate.pstar.idx(P, el.idx = x.idx, meanstructure = TRUE)
    SC[, !not.x] <- 0
  }

  # handle clustering
  if (!is.null(cluster.idx)) {
    # take the sum within each cluster
    SC <- rowsum(SC, group = cluster.idx, reorder = FALSE, na.rm = TRUE)

    # lower bias is number of clusters is not very high
    nC <- nrow(SC)
    correction.factor <- nC / (nC - 1)
    SC <- SC * sqrt(correction.factor)
  }

  # first order information
  Bbeta <- lav_matrix_crossprod(SC) / N

  # expected/observed information
  if (information == "expected") {
    Abeta <- lav_matrix_bdiag(I11, I22) / N
  } else {
    Abeta <- rbind(
      cbind(I11, t(I21)),
      cbind(I21, I22)
    ) / N
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    not.x <- eliminate.pstar.idx(
      nvar = P, el.idx = x.idx,
      meanstructure = TRUE
    )
    Abeta[!not.x, ] <- 0
    Abeta[, !not.x] <- 0
  }

  list(Abeta = Abeta, Bbeta = Bbeta)
}


# 6) inverted information h0 mu + vech(Sigma)

#    6a: (unit) inverted expected information
#    NOT USED: is not equal to solve(expected)
#    (although it does converge to the same solution eventually)
# lav_mvnorm_missing_inverted_information_expected <- function(Y  = NULL,
#                                                    Mp          = NULL,
#                                                    Mu          = NULL,# unused
#                                                    Sigma       = NULL) {
#    P <- NCOL(Y)
#
#    # missing patterns
#    if(is.null(Mp)) {
#        Mp <- lav_data_missing_patterns(Y)
#    }
#
#    # N
#    N <- sum(Mp$freq) # removed empty cases!
#
#    I11 <- matrix(0, P, P)
#    I22 <- matrix(0, P*(P+1)/2, P*(P+1)/2)
#
#    # for each pattern
#    for(p in seq_len(Mp$npatterns)) {
#
#        # observed variables
#        var.idx <- Mp$pat[p,]
#
#        sigma <- matrix(0, P, P)
#        sigma[var.idx, var.idx] <- Sigma[var.idx, var.idx]
#        sigma2 <- 2 * lav_matrix_duplication_ginv_pre_post(sigma %x% sigma)
#
#        I11 <- I11 + Mp$freq[p] * sigma
#        I22 <- I22 + Mp$freq[p] * sigma2
#    }
#
#    lav_matrix_bdiag(I11, I22)/N
# }

#    6b: /

#    6c: /


# 7) ACOV h0 mu + vech(Sigma)

#    7a: 1/N * inverted expected    information

#    7b: 1/N * inverted observed    information

#    7c: 1/N * inverted first-order information

#    7d: sandwich acov



# 10) other stuff

# single imputation missing cells, under the normal model, pattern-based
# FIXME: add wt
lav_mvnorm_missing_impute_pattern <- function(Y = NULL,
                                              Mp = NULL,
                                              Mu = NULL,
                                              Sigma = NULL,
                                              Sigma.inv = NULL,
                                              Sinv.method = "eigen") {
  Mu <- as.numeric(Mu)

  # complete data
  Y.complete <- Y

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # fill in data per pattern
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # if complete, nothing to do
    if (all(var.idx)) {
      next
    }

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # extract observed data for these (centered) cases
    Oc <- Yc[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

    # invert Sigma (Sigma_22, observed part only) for this pattern
    Sigma_22.inv <- try(
      lav_matrix_symmetric_inverse_update(
        S.inv =
          Sigma.inv, rm.idx = na.idx, logdet = FALSE
      ),
      silent = TRUE
    )
    if (inherits(Sigma_22.inv, "try-error")) {
      lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
    }

    # estimate missing values in this pattern
    Sigma_12 <- Sigma[!var.idx, var.idx, drop = FALSE]
    Y.missing <- t(Sigma_12 %*% Sigma_22.inv %*% t(Oc) + Mu[!var.idx])

    # complete data for this pattern
    Y.complete[Mp$case.idx[[p]], !var.idx] <- Y.missing
  }

  Y.complete
}


# E-step: expectations of sum, sum of squares, sum of crossproducts
# plus correction
lav_mvnorm_missing_estep <- function(Y = NULL,
                                     Mp = NULL,
                                     wt = NULL,
                                     Mu = NULL,
                                     Sigma = NULL,
                                     Sigma.inv = NULL,
                                     Sinv.method = "eigen") {
  P <- NCOL(Y)
  Mu <- as.numeric(Mu)

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  if (is.null(Sigma.inv)) {
    # invert Sigma
    Sigma.inv <- lav_matrix_symmetric_inverse(
      S = Sigma, logdet = FALSE,
      Sinv.method = Sinv.method
    )
  }

  # T1, T2
  T1 <- numeric(P)
  T2 <- matrix(0, P, P)

  # update T1 and T2 per pattern
  for (p in seq_len(Mp$npatterns)) {
    # observed values for this pattern
    var.idx <- Mp$pat[p, ]

    # extract observed data
    O <- Y[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

    # if complete, just compute first and second moments
    if (all(var.idx)) {
      if (!is.null(wt)) {
        WT <- wt[Mp$case.idx[[p]]]
        T1 <- T1 + colSums(WT * O)
        T2 <- T2 + crossprod(sqrt(WT) * O)
      } else {
        # complete pattern
        T1 <- T1 + colSums(O)
        T2 <- T2 + crossprod(O)
      }
      next
    }

    # missing values for this pattern
    na.idx <- which(!var.idx)

    # partition Sigma (1=missing, 2=complete)
    Sigma_11 <- Sigma[!var.idx, !var.idx, drop = FALSE]
    Sigma_12 <- Sigma[!var.idx, var.idx, drop = FALSE]
    Sigma_21 <- Sigma[var.idx, !var.idx, drop = FALSE]

    # invert Sigma (Sigma_22, observed part only) for this pattern
    Sigma_22.inv <- try(
      lav_matrix_symmetric_inverse_update(
        S.inv =
          Sigma.inv, rm.idx = na.idx, logdet = FALSE
      ),
      silent = TRUE
    )
    if (inherits(Sigma_22.inv, "try-error")) {
      lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
    }

    # estimate missing values in this pattern
    Oc <- t(t(O) - Mu[var.idx])
    Y.missing <- t(Sigma_12 %*% Sigma_22.inv %*% t(Oc) + Mu[!var.idx])

    # complete data for this pattern
    Y.complete <- matrix(0, Mp$freq[[p]], P)
    Y.complete[, var.idx] <- O
    Y.complete[, !var.idx] <- Y.missing

    if (!is.null(wt)) {
      WT <- wt[Mp$case.idx[[p]]]
      T1.pat <- colSums(WT * Y.complete)
      T2.pat <- crossprod(sqrt(WT) * Y.complete)
    } else {
      # 1. SUM `completed' pattern
      T1.pat <- colSums(Y.complete)

      # 2. CROSSPROD `completed' pattern
      T2.pat <- crossprod(Y.complete)
    }

    # correction for missing cells: conditional covariances
    T2.p11 <- Sigma_11 - (Sigma_12 %*% Sigma_22.inv %*% Sigma_21)
    if (!is.null(wt)) {
      T2.pat[!var.idx, !var.idx] <-
        T2.pat[!var.idx, !var.idx] + (T2.p11 * sum(WT))
    } else {
      T2.pat[!var.idx, !var.idx] <-
        T2.pat[!var.idx, !var.idx] + (T2.p11 * Mp$freq[[p]])
    }

    # accumulate
    T1 <- T1 + T1.pat
    T2 <- T2 + T2.pat
  }

  list(T1 = T1, T2 = T2)
}
