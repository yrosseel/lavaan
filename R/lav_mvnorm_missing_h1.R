# the Multivariate normal distribution, unrestricted (h1), missing values

# 1) loglikelihood --> same as h0 but where Mu and Sigma are unrestricted
# 2) 3) 4) 5)      --> (idem)

# YR 26 Mar 2016: first version
# YR 20 Jan 2017: added _h1_omega_sw()

# here, we estimate Mu and Sigma from Y with missing values, assuming normality
# this is a rewrite of the 'estimate.moments.EM' function in <= 0.5-22
lav_mvn_mi_h1_est_moments <- function(y = NULL,
                                                   mp = NULL,
                                                   yp = NULL,
                                                   wt = NULL,
                                                   sinv_method = "eigen",
                                                   max_iter = 500L,
                                                   tol = 1e-05) {
  # check input
  y <- as.matrix(y)
  p <- NCOL(y)
  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }
  if (is.null(yp)) {
    yp <- lav_samp_mi_patterns(y = y, mp = mp, wt = wt)
  }

  # covariances with zero coverage (perhaps planned?)
  zero_coverage_flag <- FALSE
  zero_coverage_idx <- which(mp$coverage == 0) # as a vector
  if (length(zero_coverage_idx) > 0L) {
    zero_coverage_flag <- TRUE
    # this is a problem; the current implementation of the E-step does not
    # take this into account; for now, we need to switch to nlminb() instead
    lav_msg_warn(
      gettext("some covariances have zero coverage; the current implementation
               of the EM algorithm will ignore this!"))
  }

  if (is.null(max_iter)) {
    max_iter <- 500L
  }
  if (is.null(tol)) {
    tol <- 1e-05
  }

  # remove empty cases
  n_full <- n
  if (length(mp$empty.idx) > 0L) {
    if (!is.null(wt)) {
      n <- n - sum(wt[mp$empty.idx])
    } else {
      n <- n - length(mp$empty.idx)
    }
  }

  # verbose?
  if (lav_verbose()) {
    cat("\n")
    cat("lav_mvn_mi_h1_est_moments: start EM steps\n")
  }

  # starting values; zero covariances to guarantee a pd matrix
  if (!is.null(wt)) {
    tmp <- na.omit(cbind(wt, y))
    if (nrow(tmp) > 2L) {
      y_tmp <- tmp[, -1, drop = FALSE]
      wt_tmp <- tmp[, 1]
      out <- stats::cov.wt(y_tmp, wt = wt_tmp, method = "ML")
      mu0 <- out$center
      var0 <- diag(out$cov)
    } else {
      mu0 <- base::.colMeans(y, m = n_full, n = p, na.rm = TRUE)
      yc <- t(t(y) - mu0)
      var0 <- base::.colMeans(yc * yc, m = n_full, n = p, na.rm = TRUE)
    }
  } else {
    mu0 <- base::.colMeans(y, m = n_full, n = p, na.rm = TRUE)
    yc <- t(t(y) - mu0)
    var0 <- base::.colMeans(yc * yc, m = n_full, n = p, na.rm = TRUE)
  }
  # sanity check
  bad_idx <- which(!is.finite(var0) | var0 == 0)
  if (length(bad_idx) > 0L) {
    var0[bad_idx] <- 1
  }
  bad_idx <- which(!is.finite(mu0))
  if (length(bad_idx) > 0L) {
    mu0[bad_idx] <- 0
  }
  sigma0 <- diag(x = var0, nrow = p)
  mu <- mu0
  sigma_1 <- sigma0

  # report
  if (lav_verbose()) {
    # fx0 <- lav_model_objective_fiml(Sigma.hat=Sigma, Mu.hat=Mu, M=Yp)
    fx0 <- lav_mvn_mi_loglik_samp(
      yp = yp,
      mu = mu, sigma_1 = sigma_1,
      log2pi = FALSE,
      minus_two = TRUE
    ) / n
    cat(
      "  EM iteration:", sprintf("%4d", 0),
      " fx = ", sprintf("%15.10f", fx0),
      "\n"
    )
  }

  # EM steps
  for (i in 1:max_iter) {
    # E-step
    estep <- lav_mvn_mi_estep(
      y = y, mp = mp, wt = wt,
      mu = mu, sigma_1 = sigma_1,
      sinv_method = sinv_method
    )
    t1 <- estep$T1
    t2 <- estep$T2

    # M-step
    mu <- t1 / n
    sigma_1 <- t2 / n - tcrossprod(mu)

    # check if Sigma is near-pd (+ poor fix)
    ev <- eigen(sigma_1, symmetric = TRUE, only.values = TRUE)
    evtol <- 1e-6 # FIXME!
    if (any(ev$values < evtol)) {
      # too.small <- which( ev$values < tol )
      # ev$values[too.small] <- tol
      # ev$values <- ev$values + tol
      # Sigma <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)

      # ridge
      diag(sigma_1) <- diag(sigma_1) + max(diag(sigma_1)) * 1e-08
    }

    # max absolute difference in parameter values
    mm_delta <- max(abs(c(mu, lav_mat_vech(sigma_1)) -
      c(mu0, lav_mat_vech(sigma0))))

    # report fx
    if (lav_verbose()) {
      # fx <- lav_model_objective_fiml(Sigma.hat=Sigma, Mu.hat=Mu, M=Yp)
      fx <- lav_mvn_mi_loglik_samp(
        yp = yp,
        mu = mu, sigma_1 = sigma_1,
        log2pi = FALSE,
        minus_two = TRUE
      ) / n
      cat(
        "  EM iteration:", sprintf("%4d", i),
        " fx = ", sprintf("%15.10f", fx),
        " delta par = ", sprintf("%9.8f", mm_delta),
        "\n"
      )
    }

    # convergence check: using parameter values:
    if (mm_delta < tol) {
      break
    }

    # again
    mu0 <- mu
    sigma0 <- sigma_1
  } # EM iterations

  if (lav_verbose()) {
    cat("\nSigma:\n")
    print(sigma_1)
    cat("\nMu:\n")
    print(mu)
    cat("\n")
  }

  # compute fx if we haven't already
  if (!lav_verbose()) {
    # fx <- lav_model_objective_fiml(Sigma.hat = Sigma, Mu.hat = Mu, M = Yp)
    fx <- lav_mvn_mi_loglik_samp(
      yp = yp,
      mu = mu, sigma_1 = sigma_1,
      log2pi = FALSE,
      minus_two = TRUE
    ) / n
  }

  # warning?
  if (i == max_iter) {
    lav_msg_warn(
      gettext("Maximum number of iterations reached when computing the sample
              moments using EM; use the em.h1.iter.max= argument to increase
              the number of iterations")
    )
  }

  ev <- eigen(sigma_1, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < 1e-05)) { # make an option?
    lav_msg_warn(
      gettext("The smallest eigenvalue of the EM estimated variance-covariance
              matrix (Sigma) is smaller than 1e-05; this may cause numerical
              instabilities; interpret the results with caution.")
    )
  }

  list(Sigma = sigma_1, Mu = mu, fx = fx)
}

# if we cannot use the EM algorithm (zero coverage?), we can always use
# plain FIML and nlminb() instead
#
# single level only
lav_mvn_mi_h1_est_moments_chol <- function(lavdata = NULL,
                                                        lavsamplestats = NULL,
                                                        lavoptions = NULL,
                                                        group = 1L) {
  # not for multilevel
  stopifnot(lavdata@nlevels == 1L)

  # no ov.names.x in lavdata
  lavdata@ov.names.x <- vector("list", length = lavdata@ngroups)

  # construct unrestricted partable (using chol parameterization)
  # for this group only
  lavpartable <- lav_pt_unrestricted_chol(
                             lavdata = lavdata, lavoptions = lavoptions,
                             lavpta = NULL, group = group)

  lavoptions2 <- lavoptions
  lavoptions2$estimator <- "ML"
  lavoptions2$missing <- "ml"
  lavoptions2$se <- "none"
  lavoptions2$test <- "none"
  lavoptions2$do.fit <- TRUE
  lavoptions2$optim.method <- "nlminb"
  lavoptions2$h1 <- FALSE
  lavoptions2$implied <- TRUE
  lavoptions2$loglik <- TRUE
  lavoptions2$baseline <- FALSE
  lavoptions2$fixed.x <- FALSE # even if model uses fixed.x=TRUE
  lavoptions2$model.type <- "unrestricted"
  lavoptions2$optim.attempts <- 4L
  lavoptions2$check.gradient <- FALSE
  lavoptions2$optim.force.converged <- TRUE # for now...
  lavoptions2$control <- list(rel.tol = 1e-7)
  lavoptions2$start <- "simple" # add this point, we have no lavh1 yet!
  fit <- lavaan(lavpartable,
    slot_options = lavoptions2,
    slot_sample_stats = lavsamplestats,
    slot_data = lavdata,
    warn = FALSE
  )

  out <- list(Sigma = fit@implied$cov[[1]],
              Mu = fit@implied$mean[[1]],
              fx = fit@optim$fx)

  out
}

# compute N times ACOV(Mu, vech(Sigma))
# in the literature: - `Omega_{SW}'
#                    - `Gamma for incomplete data'
#                    - (N times the) sandwich estimator for acov(mu,vech(Sigma))
lav_mvn_mi_h1_omega_sw <- function(y = NULL,
                                           mp = NULL,
                                           wt = NULL,
                                           cluster_idx = NULL,
                                           yp = NULL,
                                           sinv_method = "eigen",
                                           mu = NULL,
                                           sigma_1 = NULL,
                                           x_idx = integer(0L),
                                           sigma_inv = NULL,
                                           information = "observed") {
  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # sample stats per pattern
  if (is.null(yp) && (information == "observed" || is.null(sigma_1))) {
    yp <- lav_samp_mi_patterns(y = y, mp = mp, wt = wt)
  }

  # Sigma and Mu
  if (is.null(sigma_1) || is.null(mu)) {
    out <- lav_mvn_mi_h1_est_moments(y = y, mp = mp, yp = yp)
    mu <- out$Mu
    sigma_1 <- out$Sigma
  }

  # information matrices
  info <- lav_mvn_mi_info_both(
    y = y, mp = mp, mu = mu,
    wt = wt, cluster_idx = cluster_idx,
    sigma_1 = sigma_1, x_idx = x_idx, sinv_method = sinv_method,
    sigma_inv = sigma_inv, information = information
  )

  a <- info$Abeta
  a_inv <- lav_mat_sym_inverse(
    s = a, logdet = FALSE,
    sinv_method = sinv_method
  )
  m_b <- info$Bbeta

  # sandwich
  sw <- a_inv %*% m_b %*% a_inv

  sw
}
