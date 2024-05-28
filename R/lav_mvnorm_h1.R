# the multivariate normal distribution, unrestricted (h1)
# - everything is evalued under the MLEs: Mu = ybar, Sigma = S

# 1) loglikelihood h1 (from raw data, or sample statistics)
# 4) hessian h1 around MLEs
# 5) information h1 (restricted Sigma/mu)
#    5a: (unit)    expected information h1 (A1 = Gamma.NT^{-1})
#    5b: (unit)    observed information h1 (A1 = Gamma.NT^{-1})
#    5c: (unit) first.order information h1 (B1 = A1 %*% Gamma %*% A1)
# 6) inverted information h1 mu + vech(Sigma)
#    6a: (unit) inverted expected information (A1.inv = Gamma.NT)
#    6b: (unit) inverted observed information (A1.inv = Gamma.NT)
#    6c: (unit) inverted first-order information (B1.inv)
# 7) ACOV h1 mu + vech(Sigma)
#    7a: 1/N * Gamma.NT
#    7b: 1/N * Gamma.NT
#    7c: 1/N * (Gamma.NT * Gamma^{-1} * Gamma.NT)
#    7d: 1/N * Gamma (sandwich)


# YR 25 Mar 2016: first version
# YR 19 Jan 2017: added 6) + 7)
# YR 04 Jan 2020: adjust for sum(wt) != N
# YR 22 Jul 2022: adding correlation= argument for information_expected

# 1. log-likelihood h1

# 1a: input is raw data
lav_mvnorm_h1_loglik_data <- function(
    Y = NULL,
    x.idx = integer(0L),
    casewise = FALSE,
    wt = NULL,
    Sinv.method = "eigen") {

  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  P <- NCOL(Y)

  # sample statistics
  if (!is.null(wt)) {
    out <- stats::cov.wt(Y, wt = wt, method = "ML")
    sample.mean <- out$center
    sample.cov <- out$cov
  } else {
    sample.mean <- base::.colMeans(Y, m = N, n = P)
    sample.cov <- lav_matrix_cov(Y)
  }

  if (casewise) {
    LOG.2PI <- log(2 * pi)

    # invert sample.cov
    if (Sinv.method == "chol") {
      cS <- chol(sample.cov)
      icS <- backsolve(cS, diag(P))
      Yc <- t(t(Y) - sample.mean)
      DIST <- rowSums((Yc %*% icS)^2)
      logdet <- -2 * sum(log(diag(icS)))
    } else {
      sample.cov.inv <- lav_matrix_symmetric_inverse(
        S = sample.cov,
        logdet = TRUE, Sinv.method = Sinv.method
      )
      logdet <- attr(sample.cov.inv, "logdet")
      # mahalanobis distance
      Yc <- t(t(Y) - sample.mean)
      DIST <- rowSums(Yc %*% sample.cov.inv * Yc)
    }

    loglik <- -(P * LOG.2PI + logdet + DIST) / 2

    # weights
    if (!is.null(wt)) {
      loglik <- loglik * wt
    }
  } else {
    # invert sample.cov
    sample.cov.inv <- lav_matrix_symmetric_inverse(
      S = sample.cov,
      logdet = TRUE, Sinv.method = Sinv.method
    )
    logdet <- attr(sample.cov.inv, "logdet")

    loglik <-
      lav_mvnorm_h1_loglik_samplestats(
        sample.cov.logdet = logdet,
        sample.nvar = P,
        sample.nobs = N
      )
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    loglik.x <- lav_mvnorm_h1_loglik_data(
      Y = Y[, x.idx, drop = FALSE],
      wt = wt, x.idx = integer(0L),
      casewise = casewise,
      Sinv.method = Sinv.method
    )
    # subtract logl.X
    loglik <- loglik - loglik.x
  }

  loglik
}



# 1b: input are sample statistics only (logdet, N and P)
lav_mvnorm_h1_loglik_samplestats <- function(
    sample.cov.logdet = NULL,
    sample.nvar = NULL,
    sample.nobs = NULL,
    # or
    sample.cov = NULL,
    x.idx = integer(0L),
    x.cov = NULL,
    Sinv.method = "eigen") {

  if (is.null(sample.nvar)) {
    P <- NCOL(sample.cov)
  } else {
    P <- sample.nvar # number of variables
  }

  N <- sample.nobs
  stopifnot(!is.null(P), !is.null(N))

  LOG.2PI <- log(2 * pi)

  # all we need is the logdet
  if (is.null(sample.cov.logdet)) {
    sample.cov.inv <- lav_matrix_symmetric_inverse(
      S = sample.cov,
      logdet = TRUE, Sinv.method = Sinv.method
    )
    logdet <- attr(sample.cov.inv, "logdet")
  } else {
    logdet <- sample.cov.logdet
  }

  loglik <- -N / 2 * (P * LOG.2PI + logdet + P)

  # fixed.x?
  if (length(x.idx) > 0L) {
    if (is.null(sample.cov)) {
      if (is.null(x.cov)) {
        lav_msg_stop(gettext(
          "when x.idx is not empty, we need sample.cov or x.cov"
        ))
      } else {
        sample.cov.x <- x.cov
      }
    } else {
      sample.cov.x <- sample.cov[x.idx, x.idx, drop = FALSE]
    }

    loglik.x <-
      lav_mvnorm_h1_loglik_samplestats(
        sample.cov = sample.cov.x,
        sample.nobs = sample.nobs,
        x.idx = integer(0L),
        Sinv.method = Sinv.method
      )
    # subtract logl.X
    loglik <- loglik - loglik.x
  }

  loglik
}


# 4. hessian of logl (around MLEs of Mu and Sigma)

# 4a: hessian logl Mu and vech(Sigma) from raw data
lav_mvnorm_h1_logl_hessian_data <- function(
    Y = NULL,
    wt = NULL,
    x.idx = integer(0L),
    Sinv.method = "eigen",
    sample.cov.inv = NULL,
    meanstructure = TRUE) {

  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  # observed information
  observed <- lav_mvnorm_h1_information_observed_data(
    Y = Y, wt = wt,
    x.idx = x.idx, Sinv.method = Sinv.method,
    sample.cov.inv = sample.cov.inv,
    meanstructure = meanstructure
  )

  -N * observed
}

# 4b: hessian Mu and vech(Sigma) from samplestats
lav_mvnorm_h1_logl_hessian_samplestats <- function(
    sample.mean = NULL, # unused!
    sample.cov = NULL,
    sample.nobs = NULL,
    x.idx = integer(0L),
    Sinv.method = "eigen",
    sample.cov.inv = NULL,
    meanstructure = TRUE) {

  N <- sample.nobs

  # observed information
  observed <- lav_mvnorm_h1_information_observed_samplestats(
    sample.mean =
      sample.mean, sample.cov = sample.cov, x.idx = x.idx,
    Sinv.method = Sinv.method, sample.cov.inv = sample.cov.inv,
    meanstructure = meanstructure
  )

  -N * observed
}



# 5) Information h1 (note: expected == observed if data is complete!)

# 5a: unit expected information h1
lav_mvnorm_h1_information_expected <- function(
    Y = NULL,
    wt = NULL,
    sample.cov = NULL,
    x.idx = integer(0L),
    Sinv.method = "eigen",
    sample.cov.inv = NULL,
    meanstructure = TRUE,
    correlation = FALSE) {

  if (is.null(sample.cov.inv)) {
    if (is.null(sample.cov)) {
      if (is.null(wt)) {
        sample.mean <- base::.colMeans(Y, m = NROW(Y), n = NCOL(Y))
        sample.cov <- lav_matrix_cov(Y)
      } else {
        out <- stats::cov.wt(Y, wt = wt, method = "ML")
        sample.cov <- out$cov
      }
    }

    # invert sample.cov
    sample.cov.inv <- lav_matrix_symmetric_inverse(
      S = sample.cov,
      logdet = FALSE, Sinv.method = Sinv.method
    )
  }

  I11 <- sample.cov.inv
  if (correlation) {
    I22 <- 0.5 * lav_matrix_duplication_cor_pre_post(sample.cov.inv %x%
      sample.cov.inv)
  } else {
    I22 <- 0.5 * lav_matrix_duplication_pre_post(sample.cov.inv %x%
      sample.cov.inv)
  }

  # fixed.x?
  if (length(x.idx) > 0L) {
    pstar.x <- lav_matrix_vech_which_idx(
      n = NCOL(sample.cov.inv), idx = x.idx,
      diagonal = !correlation
    )
    I22[pstar.x, ] <- 0
    I22[, pstar.x] <- 0
  }

  if (meanstructure) {
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

# 5b: unit observed information h1
lav_mvnorm_h1_information_observed_data <- function(
    Y = NULL,
    wt = NULL,
    x.idx = integer(0L),
    Sinv.method = "eigen",
    sample.cov.inv = NULL,
    meanstructure = TRUE,
    correlation = FALSE) { # todo

  lav_mvnorm_h1_information_expected(
    Y = Y, Sinv.method = Sinv.method,
    wt = wt, x.idx = x.idx,
    sample.cov.inv = sample.cov.inv,
    meanstructure = meanstructure
  )
}

# 5b-bis: observed information h1 from sample statistics
lav_mvnorm_h1_information_observed_samplestats <- function(
    sample.mean = NULL, # unused!
    sample.cov = NULL,
    x.idx = integer(0L),
    Sinv.method = "eigen",
    sample.cov.inv = NULL,
    meanstructure = TRUE,
    correlation = FALSE) { # todo

  if (is.null(sample.cov.inv)) {
    # invert sample.cov
    sample.cov.inv <- lav_matrix_symmetric_inverse(
      S = sample.cov,
      logdet = FALSE, Sinv.method = Sinv.method
    )
  }

  I11 <- sample.cov.inv
  # fixed.x?
  if (length(x.idx) > 0L) {
    I11[x.idx, ] <- 0
    I11[, x.idx] <- 0
  }

  I22 <- 0.5 * lav_matrix_duplication_pre_post(sample.cov.inv %x%
    sample.cov.inv)
  # fixed.x?
  if (length(x.idx) > 0L) {
    pstar.x <- lav_matrix_vech_which_idx(
      n = NCOL(sample.cov.inv),
      idx = x.idx, diagonal = !correlation
    )
    I22[pstar.x, ] <- 0
    I22[, pstar.x] <- 0
  }

  if (meanstructure) {
    out <- lav_matrix_bdiag(I11, I22)
  } else {
    out <- I22
  }

  out
}

# 5c: unit first-order information h1
#     note: first order information h1 == A1 %*% Gamma %*% A1
#           (where A1 = obs/exp information h1)
lav_mvnorm_h1_information_firstorder <- function(
    Y = NULL,
    wt = NULL,
    sample.cov = NULL,
    x.idx = integer(0L),
    cluster.idx = NULL,
    Sinv.method = "eigen",
    sample.cov.inv = NULL,
    Gamma = NULL,
    meanstructure = TRUE) {

  if (!is.null(wt)) {
    out <- stats::cov.wt(Y, wt = wt, method = "ML")
    res <- lav_mvnorm_information_firstorder(
      Y = Y, wt = wt,
      cluster.idx = cluster.idx,
      Mu = out$center, Sigma = out$cov, x.idx = x.idx,
      meanstructure = meanstructure
    )
    return(res)
  }

  # question: is there any benefit computing Gamma/A1 instead of just
  # calling lav_mvnorm_information_firstorder()?

  # Gamma
  # FIXME: what about the 'unbiased = TRUE' option?
  if (is.null(Gamma)) {
    if (length(x.idx) > 0L) {
      Gamma <- lav_samplestats_Gamma(Y,
        x.idx = x.idx, fixed.x = TRUE,
        cluster.idx = cluster.idx,
        meanstructure = meanstructure
      )
    } else {
      Gamma <- lav_samplestats_Gamma(Y,
        meanstructure = meanstructure,
        cluster.idx = cluster.idx
      )
    }
  }

  # sample.cov.inv
  if (is.null(sample.cov.inv)) {
    # invert sample.cov
    if (is.null(sample.cov)) {
      sample.mean <- base::.colMeans(Y, m = NROW(Y), n = NCOL(Y))
      sample.cov <- lav_matrix_cov(Y)
    }
    sample.cov.inv <- lav_matrix_symmetric_inverse(
      S = sample.cov,
      logdet = FALSE, Sinv.method = Sinv.method
    )
  }

  # A1
  A1 <- lav_mvnorm_h1_information_expected(
    Y = Y, Sinv.method = Sinv.method,
    sample.cov.inv = sample.cov.inv,
    x.idx = x.idx,
    meanstructure = meanstructure
  )

  A1 %*% Gamma %*% A1
}

# 6) inverted information h1 mu + vech(Sigma)

#    6a: (unit) inverted expected information (A1.inv = Gamma.NT)
#    6b: (unit) inverted observed information (A1.inv = Gamma.NT)

lav_mvnorm_h1_inverted_information_expected <-
lav_mvnorm_h1_inverted_information_observed <- function(
    Y = NULL,
    wt = NULL,
    sample.cov = NULL,
    x.idx = integer(0L)) {

  # sample.cov
  if (is.null(sample.cov)) {
    if (is.null(wt)) {
      sample.mean <- base::.colMeans(Y, m = NROW(Y), n = NCOL(Y))
      sample.cov <- lav_matrix_cov(Y)
    } else {
      out <- stats::cov.wt(Y, wt = wt, method = "ML")
      sample.cov <- out$cov
    }
  }

  if (length(x.idx) > 0L) {
    Gamma.NT <- lav_samplestats_Gamma_NT(
      Y = Y, wt = wt, x.idx = x.idx,
      COV = sample.cov,
      meanstructure = TRUE,
      fixed.x = TRUE
    )
  } else {
    I11 <- sample.cov
    I22 <- 2 * lav_matrix_duplication_ginv_pre_post(sample.cov %x% sample.cov)

    Gamma.NT <- lav_matrix_bdiag(I11, I22)
  }

  Gamma.NT
}

#    6c: (unit) inverted first-order information (B1.inv)
#        J1.inv = Gamma.NT %*% solve(Gamma) %*%  Gamma.NT
#
lav_mvnorm_h1_inverted_information_firstorder <- function(
    Y = NULL,
    wt = NULL,
    sample.cov = NULL,
    x.idx = integer(0L),
    Sinv.method = "eigen",
    sample.cov.inv = NULL,
    Gamma = NULL) {

  # lav_samplestats_Gamma() has no wt argument (yet)
  if (!is.null(wt)) {
    lav_msg_stop(gettext("function not supported if wt is not NULL"))
  }

  # Gamma
  # what about the 'unbiased = TRUE' option?
  if (is.null(Gamma)) {
    if (length(x.idx) > 0L) {
      Gamma <- lav_samplestats_Gamma(Y,
        x.idx = x.idx, fixed.x = TRUE,
        meanstructure = TRUE
      )
    } else {
      Gamma <- lav_samplestats_Gamma(Y, meanstructure = TRUE)
    }
  }

  # Gamma.NT
  Gamma.NT <-
    lav_mvnorm_h1_inverted_information_expected(
      Y = Y,
      sample.cov = sample.cov,
      x.idx = x.idx
    )
  if (length(x.idx) > 0L) {
    # FIXME: surely there is better way
    out <- Gamma.NT %*% MASS::ginv(Gamma) %*% Gamma.NT
  } else {
    out <- Gamma.NT %*% solve(Gamma, Gamma.NT)
  }

  out
}


# 7) ACOV h1 mu + vech(Sigma)

#    7a: 1/N * Gamma.NT
#    7b: 1/N * Gamma.NT
lav_mvnorm_h1_acov_expected <-
lav_mvnorm_h1_acov_observed <- function(
    Y = NULL,
    wt = NULL,
    sample.cov = NULL,
    x.idx = integer(0L)) {

  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  Gamma.NT <-
    lav_mvnorm_h1_inverted_information_expected(
      Y = Y,
      wt = wt,
      sample.cov = sample.cov,
      x.idx = x.idx
    )

  (1 / N) * Gamma.NT
}

#    7c: 1/N * (Gamma.NT * Gamma^{-1} * Gamma.NT)
lav_mvnorm_h1_acov_firstorder <- function(
    Y = NULL,
    wt = NULL,
    sample.cov = NULL,
    Sinv.method = "eigen",
    x.idx = integer(0L),
    sample.cov.inv = NULL,
    Gamma = NULL) {

  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  J1.inv <- lav_mvnorm_h1_inverted_information_firstorder(
    Y = Y, wt = wt,
    sample.cov = sample.cov,
    x.idx = x.idx, Sinv.method = Sinv.method,
    sample.cov.inv = sample.cov.inv, Gamma = Gamma
  )

  (1 / N) * J1.inv
}

#    7d: 1/N * Gamma (sandwich)
lav_mvnorm_h1_acov_sandwich <- function(
    Y = NULL,
    wt = NULL,
    sample.cov = NULL,
    x.idx = integer(0L),
    Gamma = NULL) {

  # lav_samplestats_Gamma() has no wt argument (yet)
  if (!is.null(wt)) {
    lav_msg_stop(gettext("function not supported if wt is not NULL"))
  }

  # if(!is.null(wt)) {
  #    N <- sum(wt)
  # } else {
  N <- NROW(Y)
  # }

  # Gamma
  if (is.null(Gamma)) {
    if (length(x.idx) > 0L) {
      Gamma <- lav_samplestats_Gamma(Y,
        x.idx = x.idx, fixed.x = TRUE,
        meanstructure = TRUE
      )
    } else {
      Gamma <- lav_samplestats_Gamma(Y, meanstructure = TRUE)
    }
  }

  (1 / N) * Gamma
}
