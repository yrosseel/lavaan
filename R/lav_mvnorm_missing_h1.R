# the Multivariate normal distribution, unrestricted (h1), missing values

# 1) loglikelihood --> same as h0 but where Mu and Sigma are unrestricted
# 2) 3) 4) 5)      --> (idem)

# YR 26 Mar 2016: first version
# YR 20 Jan 2017: added _h1_omega_sw()

# here, we estimate Mu and Sigma from Y with missing values, assuming normality
# this is a rewrite of the 'estimate.moments.EM' function in <= 0.5-22
lav_mvnorm_missing_h1_estimate_moments <- function(Y = NULL,
                                                   Mp = NULL,
                                                   Yp = NULL,
                                                   wt = NULL,
                                                   Sinv.method = "eigen",
                                                   verbose = FALSE,
                                                   max.iter = 500L,
                                                   tol = 1e-05,
                                                   warn = FALSE) {
  # check input
  Y <- as.matrix(Y)
  P <- NCOL(Y)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }

  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }
  if (is.null(Yp)) {
    Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp, wt = wt)
  }

  if (is.null(max.iter)) {
    max.iter <- 500L
  }
  if (is.null(tol)) {
    tol <- 1e-05
  }
  if (is.null(warn)) {
    warn <- FALSE
  }

  # remove empty cases
  N.full <- N
  if (length(Mp$empty.idx) > 0L) {
    if (!is.null(wt)) {
      N <- N - sum(wt[Mp$empty.idx])
    } else {
      N <- N - length(Mp$empty.idx)
    }
  }

  # verbose?
  if (verbose) {
    cat("\n")
    cat("lav_mvnorm_missing_h1_estimate_moments: start EM steps\n")
  }

  # starting values; zero covariances to guarantee a pd matrix
  if (!is.null(wt)) {
    tmp <- na.omit(cbind(wt, Y))
    if (nrow(tmp) > 2L) {
      Y.tmp <- tmp[, -1, drop = FALSE]
      wt.tmp <- tmp[, 1]
      out <- stats::cov.wt(Y.tmp, wt = wt.tmp, method = "ML")
      Mu0 <- out$center
      var0 <- diag(out$cov)
    } else {
      Mu0 <- base::.colMeans(Y, m = N.full, n = P, na.rm = TRUE)
      Yc <- t(t(Y) - Mu0)
      var0 <- base::.colMeans(Yc * Yc, m = N.full, n = P, na.rm = TRUE)
    }
  } else {
    Mu0 <- base::.colMeans(Y, m = N.full, n = P, na.rm = TRUE)
    Yc <- t(t(Y) - Mu0)
    var0 <- base::.colMeans(Yc * Yc, m = N.full, n = P, na.rm = TRUE)
  }
  # sanity check
  bad.idx <- which(!is.finite(var0) | var0 == 0)
  if (length(bad.idx) > 0L) {
    var0[bad.idx] <- 1
  }
  bad.idx <- which(!is.finite(Mu0))
  if (length(bad.idx) > 0L) {
    Mu0[bad.idx] <- 0
  }
  Sigma0 <- diag(x = var0, nrow = P)
  Mu <- Mu0
  Sigma <- Sigma0

  # report
  if (verbose) {
    # fx0 <- estimator.FIML(Sigma.hat=Sigma, Mu.hat=Mu, M=Yp)
    fx0 <- lav_mvnorm_missing_loglik_samplestats(
      Yp = Yp,
      Mu = Mu, Sigma = Sigma,
      log2pi = FALSE,
      minus.two = TRUE
    ) / N
    cat(
      "  EM iteration:", sprintf("%4d", 0),
      " fx = ", sprintf("%15.10f", fx0),
      "\n"
    )
  }

  # EM steps
  for (i in 1:max.iter) {
    # E-step
    Estep <- lav_mvnorm_missing_estep(
      Y = Y, Mp = Mp, wt = wt,
      Mu = Mu, Sigma = Sigma,
      Sinv.method = Sinv.method
    )
    T1 <- Estep$T1
    T2 <- Estep$T2

    # M-step
    Mu <- T1 / N
    Sigma <- T2 / N - tcrossprod(Mu)

    # check if Sigma is near-pd (+ poor fix)
    ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)
    evtol <- 1e-6 # FIXME!
    if (any(ev$values < evtol)) {
      # too.small <- which( ev$values < tol )
      # ev$values[too.small] <- tol
      # ev$values <- ev$values + tol
      # Sigma <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)

      # ridge
      diag(Sigma) <- diag(Sigma) + max(diag(Sigma)) * 1e-08
    }

    # max absolute difference in parameter values
    DELTA <- max(abs(c(Mu, lav_matrix_vech(Sigma)) -
      c(Mu0, lav_matrix_vech(Sigma0))))

    # report fx
    if (verbose) {
      # fx <- estimator.FIML(Sigma.hat=Sigma, Mu.hat=Mu, M=Yp)
      fx <- lav_mvnorm_missing_loglik_samplestats(
        Yp = Yp,
        Mu = Mu, Sigma = Sigma,
        log2pi = FALSE,
        minus.two = TRUE
      ) / N
      cat(
        "  EM iteration:", sprintf("%4d", i),
        " fx = ", sprintf("%15.10f", fx),
        " delta par = ", sprintf("%9.8f", DELTA),
        "\n"
      )
    }

    # convergence check: using parameter values:
    if (DELTA < tol) {
      break
    }

    # again
    Mu0 <- Mu
    Sigma0 <- Sigma
  } # EM iterations

  if (verbose) {
    cat("\nSigma:\n")
    print(Sigma)
    cat("\nMu:\n")
    print(Mu)
    cat("\n")
  }

  # compute fx if we haven't already
  if (!verbose) {
    # fx <- estimator.FIML(Sigma.hat = Sigma, Mu.hat = Mu, M = Yp)
    fx <- lav_mvnorm_missing_loglik_samplestats(
      Yp = Yp,
      Mu = Mu, Sigma = Sigma,
      log2pi = FALSE,
      minus.two = TRUE
    ) / N
  }

  # warning?
  if (warn && i == max.iter) {
    lav_msg_warn(
      gettext("Maximum number of iterations reached when computing the sample
              moments using EM; use the em.h1.iter.max= argument to increase
              the number of iterations")
    )
  }

  if (warn) {
    ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < 1e-05)) { # make an option?
      lav_msg_warn(
        gettext("The smallest eigenvalue of the EM estimated variance-covariance
                matrix (Sigma) is smaller than 1e-05; this may cause numerical
                instabilities; interpret the results with caution.")
      )
    }
  }

  list(Sigma = Sigma, Mu = Mu, fx = fx)
}

# compute N times ACOV(Mu, vech(Sigma))
# in the literature: - `Omega_{SW}'
#                    - `Gamma for incomplete data'
#                    - (N times the) sandwich estimator for acov(mu,vech(Sigma))
lav_mvnorm_missing_h1_omega_sw <- function(Y = NULL,
                                           Mp = NULL,
                                           wt = NULL,
                                           cluster.idx = NULL,
                                           Yp = NULL,
                                           Sinv.method = "eigen",
                                           Mu = NULL,
                                           Sigma = NULL,
                                           x.idx = integer(0L),
                                           Sigma.inv = NULL,
                                           information = "observed") {
  # missing patterns
  if (is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # sample stats per pattern
  if (is.null(Yp) && (information == "observed" || is.null(Sigma))) {
    Yp <- lav_samplestats_missing_patterns(Y = Y, Mp = Mp, wt = wt)
  }

  # Sigma and Mu
  if (is.null(Sigma) || is.null(Mu)) {
    out <- lav_mvnorm_missing_h1_estimate_moments(Y = Y, Mp = Mp, Yp = Yp)
    Mu <- out$Mu
    Sigma <- out$Sigma
  }

  # information matrices
  info <- lav_mvnorm_missing_information_both(
    Y = Y, Mp = Mp, Mu = Mu,
    wt = wt, cluster.idx = cluster.idx,
    Sigma = Sigma, x.idx = x.idx, Sinv.method = Sinv.method,
    Sigma.inv = Sigma.inv, information = information
  )

  A <- info$Abeta
  A.inv <- lav_matrix_symmetric_inverse(
    S = A, logdet = FALSE,
    Sinv.method = Sinv.method
  )
  B <- info$Bbeta

  # sandwich
  SW <- A.inv %*% B %*% A.inv

  SW
}
