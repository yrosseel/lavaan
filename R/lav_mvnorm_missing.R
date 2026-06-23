# the multivariate normal distribution + missing values
#  (so-called 'FIML')

# 1) loglikelihood (from raw data, or sample statistics)
# 2) derivatives with respect to mu, Sigma, vech(Sigma)
# 3) casewise scores with respect to mu, vech(Sigma), mu + vech(Sigma)
# 4) hessian of mu + vech(Sigma)
# 5) (unit) information of mu + vech(Sigma)
#    5a: (unit)    expected information
#    5b: (unit)    observed information
#    5c: (unit) first.order information
#    5d: lav_mvn_mi_info_both (both observed + first.order)

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
lav_mvn_mi_loglik_data <- function(y = NULL,
                                           mu = NULL,
                                           wt = NULL,
                                           sigma_1 = NULL,
                                           x_idx = integer(0L),
                                           casewise = FALSE,
                                           pattern = TRUE,
                                           sinv_method = "eigen",
                                           log2pi = TRUE,
                                           minus_two = FALSE) {
  if (pattern) {
    llik <- lav_mvn_mi_llik_pattern(
      y = y, wt = wt, mu = mu,
      sigma_1 = sigma_1, x_idx = x_idx, sinv_method = sinv_method,
      log2pi = log2pi, minus_two = minus_two
    )
  } else {
    llik <- lav_mvn_mi_llik_casewise(
      y = y, wt = wt, mu = mu,
      sigma_1 = sigma_1, x_idx = x_idx, sinv_method = sinv_method,
      log2pi = log2pi, minus_two = minus_two
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
lav_mvn_mi_loglik_samp <- function(yp = NULL,
                                                  mu = NULL,
                                                  sigma_1 = NULL,
                                                  x_idx = integer(0L),
                                                  x_mean = NULL,
                                                  x_cov = NULL,
                                                  sinv_method = "eigen",
                                                  log2pi = TRUE,
                                                  minus_two = FALSE) {
  log_2pi <- log(2 * pi)
  pat_n <- length(yp)
  p_1 <- length(yp[[1]]$var.idx)

  # global inverse + logdet
  sigma_inv_1 <- lav_mat_sym_inverse(
    s = sigma_1, logdet = TRUE,
    sinv_method = sinv_method
  )
  sigma_logdet <- attr(sigma_inv_1, "logdet")

  # DIST/logdet per pattern
  dist_1 <- logdet <- p_log_2pi <- numeric(pat_n)

  # for each pattern, compute sigma.inv/logdet; compute DIST for all
  # observations of this pattern
  for (p in seq_len(pat_n)) {
    # observed variables for this pattern
    var_idx <- yp[[p]]$var.idx

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # constant
    p_log_2pi[p] <- sum(var_idx) * log_2pi * yp[[p]]$freq

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv_1,
        rm_idx = na_idx, logdet = TRUE, s_logdet = sigma_logdet
      )
      logdet[p] <- attr(sigma_inv, "logdet") * yp[[p]]$freq
    } else {
      sigma_inv <- sigma_inv_1
      logdet[p] <- sigma_logdet * yp[[p]]$freq
    }

    tt <- yp[[p]]$SY + tcrossprod(yp[[p]]$MY - mu[var_idx])
    dist_1[p] <- sum(sigma_inv * tt) * yp[[p]]$freq
  }

  # loglikelihood all data
  if (log2pi) {
    loglik <- sum(-(p_log_2pi + logdet + dist_1) / 2)
  } else {
    loglik <- sum(-(logdet + dist_1) / 2)
  }

  if (minus_two) {
    loglik <- -2 * loglik
  }

  # x.idx
  if (length(x_idx) > 0L) {
    # conditional loglikelihood (fixed.x = TRUE): subtract the marginal
    # loglikelihood of the x (exogenous) variables. This marginal must be
    # computed using FIML over the missing-data patterns of the x variables
    # themselves; using the complete-data formula (as we did before) is wrong
    # whenever the x variables contain missing values (eg missing = "fiml.x"),
    # as it ignores their missingness (see GitHub issue #226).
    # Note: mu[x_idx]/Sigma[x_idx, x_idx] equal the (model-implied) moments of
    # the x variables, so we do not need x_mean/x_cov.
    xp <- vector("list", length = length(yp))
    for (p in seq_len(length(yp))) {
      var_idx_full <- yp[[p]]$var.idx
      # which x variables are observed in this pattern?
      x_obs <- var_idx_full[x_idx]
      # skip patterns where no x variable is observed (no contribution)
      if (!any(x_obs)) {
        next
      }
      # positions of the observed x variables within the observed block
      obs_full <- which(var_idx_full)
      pos_in_obs <- match(x_idx[x_obs], obs_full)
      if (yp[[p]]$freq > 1L) {
        sy_x <- yp[[p]]$SY[pos_in_obs, pos_in_obs, drop = FALSE]
      } else {
        sy_x <- 0
      }
      xp[[p]] <- list(
        SY = sy_x, MY = yp[[p]]$MY[pos_in_obs],
        var.idx = x_obs, freq = yp[[p]]$freq
      )
    }
    # drop the (skipped) empty patterns
    xp <- xp[!vapply(xp, is.null, logical(1L))]

    loglik_x <- lav_mvn_mi_loglik_samp(
      yp = xp,
      mu = mu[x_idx],
      sigma_1 = sigma_1[x_idx, x_idx, drop = FALSE],
      x_idx = integer(0L),
      sinv_method = sinv_method,
      log2pi = log2pi,
      minus_two = minus_two
    )

    loglik <- loglik - loglik_x
  }

  loglik
}

## casewise loglikelihoods

# casewise Sinv.method
lav_mvn_mi_llik_casewise <- function(y = NULL,
                                             wt = NULL,
                                             mu = NULL,
                                             sigma_1 = NULL,
                                             x_idx = integer(0L),
                                             sinv_method = "eigen",
                                             log2pi = TRUE,
                                             minus_two = FALSE) {
  p <- NCOL(y)
  log_2pi <- log(2 * pi)
  mu <- as.numeric(mu)
  # if (!is.null(wt)) {
  #   n <- sum(wt)
  # } else {
  #   n <- NROW(y)
  # }
  ny <- NROW(y)

  # global inverse + logdet
  sigma_inv_1 <- lav_mat_sym_inverse(
    s = sigma_1, logdet = TRUE,
    sinv_method = sinv_method
  )
  sigma_logdet <- attr(sigma_inv_1, "logdet")

  # subtract Mu
  yc <- t(t(y) - mu)

  # DIST/logdet per case
  dist_1 <- logdet <- p_log_2pi <- rep(as.numeric(NA), ny)

  # missing pattern per case
  obs <- !is.na(y)
  p_i <- rowSums(obs)

  # constant
  p_log_2pi <- p_i * log_2pi

  # complete cases first (only an advantage if we have mostly complete
  # observations)
  other_idx <- seq_len(ny)
  complete_idx <- which(p_i == p)
  if (length(complete_idx) > 0L) {
    other_idx <- other_idx[-complete_idx]
    dist_1[complete_idx] <-
      rowSums(yc[complete_idx, , drop = FALSE] %*% sigma_inv_1 *
        yc[complete_idx, , drop = FALSE])
    logdet[complete_idx] <- sigma_logdet
  }

  # non-complete cases
  for (i in other_idx) {
    na_idx <- which(!obs[i, ])

    # catch empty cases
    if (length(na_idx) == p) next

    # invert Sigma for this pattern
    sigma_inv <- lav_mat_sym_inverse_update(
      s_inv = sigma_inv_1,
      rm_idx = na_idx, logdet = TRUE, s_logdet = sigma_logdet
    )
    logdet[i] <- attr(sigma_inv, "logdet")

    # distance for this case
    dist_1[i] <- sum(sigma_inv * crossprod(yc[i, obs[i, ], drop = FALSE]))
  }

  # compute casewise loglikelihoods
  if (log2pi) {
    llik <- -(p_log_2pi + logdet + dist_1) / 2
  } else {
    llik <- -(logdet + dist_1) / 2
  }

  # minus.two
  if (minus_two) {
    llik <- -2 * llik
  }

  # weights?
  if (!is.null(wt)) {
    llik <- llik * wt
  }

  # x.idx
  if (length(x_idx) > 0L) {
    llik_x <- lav_mvn_mi_llik_casewise(
      y = y[, x_idx, drop = FALSE],
      wt = wt, mu = mu[x_idx],
      sigma_1 = sigma_1[x_idx, x_idx, drop = FALSE],
      x_idx = integer(0L), sinv_method = sinv_method,
      log2pi = log2pi, minus_two = minus_two
    )
    llik <- llik - llik_x
  }

  llik
}

# pattern-based, but casewise loglikelihoods
lav_mvn_mi_llik_pattern <- function(y = NULL,
                                            mp = NULL,
                                            wt = NULL,
                                            mu = NULL,
                                            sigma_1 = NULL,
                                            x_idx = integer(0L),
                                            sinv_method = "eigen",
                                            log2pi = TRUE,
                                            minus_two = FALSE) {
  p_1 <- NCOL(y)
  log_2pi <- log(2 * pi)
  mu <- as.numeric(mu)
  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }
  ny <- NROW(y)

  # global inverse + logdet
  sigma_inv_1 <- lav_mat_sym_inverse(
    s = sigma_1, logdet = TRUE,
    sinv_method = sinv_method
  )
  sigma_logdet <- attr(sigma_inv_1, "logdet")

  # subtract Mu
  yc <- t(t(y) - mu)

  # DIST/logdet per case
  dist_1 <- logdet <- p_log_2pi <- rep(as.numeric(NA), ny)

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # for each pattern, compute sigma.inv/logdet; compute DIST for all
  # observations of this pattern
  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # identify cases with this pattern
    case_idx <- mp$case.idx[[p]]

    # constant
    p_log_2pi[case_idx] <- sum(var_idx) * log_2pi

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv_1,
        rm_idx = na_idx, logdet = TRUE, s_logdet = sigma_logdet
      )
      logdet[case_idx] <- attr(sigma_inv, "logdet")
    } else {
      sigma_inv <- sigma_inv_1
      logdet[case_idx] <- sigma_logdet
    }

    if (mp$freq[p] == 1L) {
      dist_1[case_idx] <- sum(sigma_inv *
        crossprod(yc[case_idx, var_idx, drop = FALSE]))
    } else {
      dist_1[case_idx] <-
        rowSums(yc[case_idx, var_idx, drop = FALSE] %*% sigma_inv *
          yc[case_idx, var_idx, drop = FALSE])
    }
  }

  # compute casewise loglikelihoods
  if (log2pi) {
    llik <- -(p_log_2pi + logdet + dist_1) / 2
  } else {
    llik <- -(logdet + dist_1) / 2
  }

  # minus.two
  if (minus_two) {
    llik <- -2 * llik
  }

  # weights?
  if (!is.null(wt)) {
    llik <- llik * wt
  }

  # x.idx -- using casewise (as patterns for Y may not be the same as
  #                          patterns for Y[,-x.idx])
  if (length(x_idx) > 0L) {
    llik_x <- lav_mvn_mi_llik_casewise(
      y = y[, x_idx, drop = FALSE],
      wt = wt, mu = mu[x_idx],
      sigma_1 = sigma_1[x_idx, x_idx, drop = FALSE],
      x_idx = integer(0L), sinv_method = sinv_method,
      log2pi = log2pi, minus_two = minus_two
    )
    llik <- llik - llik_x
  }

  llik
}




# 2. Derivatives

# 2a: derivative logl with respect to mu
lav_mvn_mi_dlogl_dmu <- function(y = NULL,
                                         wt = NULL,
                                         mu = NULL,
                                         sigma_1 = NULL,
                                         x_idx = integer(0L),
                                         sigma_inv = NULL,
                                         sinv_method = "eigen") {
  sc <- lav_mvn_mi_sc_mu(
    y = y, wt = wt, mu = mu, sigma_1 = sigma_1,
    x_idx = x_idx, sigma_inv = sigma_inv, sinv_method = sinv_method
  )

  colSums(sc, na.rm = TRUE)
}

# 2abis: using samplestats
lav_mvn_mi_dlogl_dmu_samp <- function(yp = NULL,
                                                     mu = NULL,
                                                     sigma_1 = NULL,
                                                     x_idx = integer(0L),
                                                     sigma_inv = NULL,
                                                     sinv_method = "eigen") {
  pat_n <- length(yp)
  p_1 <- length(yp[[1]]$var.idx)

  if (is.null(sigma_inv)) {
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # dmu
  dmu <- numeric(p_1)

  # for each pattern, compute sigma.inv
  for (p in seq_len(pat_n)) {
    # observed variables for this pattern
    var_idx <- yp[[p]]$var.idx

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv_1 <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
    } else {
      sigma_inv_1 <- sigma_inv
    }

    # dmu for this pattern
    dmu_pattern <- as.numeric(sigma_inv_1 %*% (yp[[p]]$MY - mu[var_idx]))

    # update mu
    dmu[var_idx] <- dmu[var_idx] + (dmu_pattern * yp[[p]]$freq)
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    dmu[x_idx] <- 0
  }

  dmu
}



# 2b: derivative logl with respect to Sigma (full matrix, ignoring symmetry)
lav_mvn_mi_dlogl_dsigma <- function(m_y = NULL,
                                            l_mp = NULL,
                                            wt = NULL,
                                            m_mu = NULL,
                                            m_sigma = NULL,
                                            x_idx = integer(0L),
                                            sigma_inv = NULL,
                                            sinv_method = "eigen") {
  p <- NCOL(m_y)
  m_mu <- as.numeric(m_mu)
  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(m_y)
  }
  ny <- NROW(m_y)

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = m_sigma, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # subtract Mu
  m_yc <- t(t(m_y) - m_mu)

  # dvechSigma
  dsigma <- matrix(0, p, p)

  # missing patterns
  if (is.null(l_mp)) {
    l_mp <- lav_data_mi_patterns(m_y)
  }

  # for each pattern
  for (p in seq_len(l_mp$npatterns)) {
    # observed values for this pattern
    var_idx <- l_mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # cases with this pattern
    case_idx <- l_mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      local_sigma_inv <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
    } else {
      local_sigma_inv <- sigma_inv
    }

    if (!is.null(wt)) {
      freq <- sum(wt[case_idx])
    } else {
      freq <- l_mp$freq[p]
    }

    if (length(case_idx) > 1L) {
      if (!is.null(wt)) {
        out <- stats::cov.wt(m_y[case_idx, var_idx, drop = FALSE],
          wt = wt[l_mp$case.idx[[p]]], method = "ML"
        )
        m_sy <- out$cov
        m_my <- out$center
        w_tilde <- m_sy + tcrossprod(m_my - m_mu[var_idx])
      } else {
        w_tilde <- crossprod(m_yc[case_idx, var_idx, drop = FALSE]) / freq
      }
    } else {
      w_tilde <- tcrossprod(m_yc[case_idx, var_idx])
    }

    # dSigma for this pattern
    dsigma_pattern <- matrix(0, p, p)
    dsigma_pattern[var_idx, var_idx] <- -(1 / 2) * (local_sigma_inv -
      (local_sigma_inv %*% w_tilde %*% local_sigma_inv))

    # update dSigma
    dsigma <- dsigma + (dsigma_pattern * freq)
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    dsigma[x_idx, x_idx] <- 0
  }

  dsigma
}

# 2bbis: using samplestats
lav_mvn_mi_dlogl_dsigma_samp <- function(yp = NULL,
                                                        mu = NULL,
                                                        sigma_1 = NULL,
                                                        x_idx = integer(0L),
                                                        sigma_inv = NULL,
                                                        sinv_method = "eigen") {
  pat_n <- length(yp)
  p_1 <- length(yp[[1]]$var.idx)

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # dvechSigma
  d_sigma <- matrix(0, p_1, p_1)

  # for each pattern
  for (p in seq_len(pat_n)) {
    # observed variables for this pattern
    var_idx <- yp[[p]]$var.idx

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv_1 <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
    } else {
      sigma_inv_1 <- sigma_inv
    }

    w_tilde <- yp[[p]]$SY + tcrossprod(yp[[p]]$MY - mu[var_idx])

    # dSigma for this pattern
    d_sigma_pattern <- matrix(0, p_1, p_1)
    d_sigma_pattern[var_idx, var_idx] <- -(1 / 2) * (sigma_inv_1 -
      (sigma_inv_1 %*% w_tilde %*% sigma_inv_1))

    # update dSigma
    d_sigma <- d_sigma + (d_sigma_pattern * yp[[p]]$freq)
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    d_sigma[x_idx, x_idx] <- 0
  }

  d_sigma
}


# 2c: derivative logl with respect to vech(Sigma)
lav_mvn_mi_dlogl_dvechsigma <- function(y = NULL,
                                                wt = NULL,
                                                mu = NULL,
                                                x_idx = integer(0L),
                                                sigma_1 = NULL,
                                                sigma_inv = NULL,
                                                sinv_method = "eigen") {
  d_sigma <- lav_mvn_mi_dlogl_dsigma(
    m_y = y, wt = wt, m_mu = mu,
    m_sigma = sigma_1, x_idx = x_idx, sigma_inv = sigma_inv,
    sinv_method = sinv_method
  )

  dvech_sigma <- as.numeric(lav_mat_dup_pre(
    as.matrix(lav_mat_vec(d_sigma))
  ))

  dvech_sigma
}

# 2cbis: using samplestats
lav_mvnorm_missing_dlogl_dvechsigma_samplestats <- # nolint
  function(yp = NULL,
           mu = NULL,
           sigma_1 = NULL,
           x_idx = integer(0L),
           sigma_inv = NULL,
           sinv_method = "eigen") {
    pat_n <- length(yp)
    p_1 <- length(yp[[1]]$var.idx)

    if (is.null(sigma_inv)) {
      # invert Sigma
      sigma_inv <- lav_mat_sym_inverse(
        s = sigma_1, logdet = FALSE,
        sinv_method = sinv_method
      )
    }

    # dvechSigma
    dvech_sigma <- numeric(p_1 * (p_1 + 1) / 2)

    # for each pattern
    for (p in seq_len(pat_n)) {
      # observed variables for this pattern
      var_idx <- yp[[p]]$var.idx

      # missing values for this pattern
      na_idx <- which(!var_idx)

      # invert Sigma for this pattern
      if (length(na_idx) > 0L) {
        sigma_inv_1 <- lav_mat_sym_inverse_update(
          s_inv = sigma_inv,
          rm_idx = na_idx, logdet = FALSE
        )
      } else {
        sigma_inv_1 <- sigma_inv
      }

      w_tilde <- yp[[p]]$SY + tcrossprod(yp[[p]]$MY - mu[var_idx])

      # dSigma for this pattern
      d_sigma_pattern <- matrix(0, p_1, p_1)
      d_sigma_pattern[var_idx, var_idx] <- -(1 / 2) * (sigma_inv_1 -
        (sigma_inv_1 %*% w_tilde %*% sigma_inv_1))

      # fixed.x?
      if (length(x_idx) > 0L) {
        d_sigma_pattern[x_idx, x_idx] <- 0
      }

      # convert to vechSigma
      dvech_sigma_pattern <- as.numeric(lav_mat_dup_pre(
        as.matrix(lav_mat_vec(d_sigma_pattern))
      ))

      # update dvechSigma
      dvech_sigma <- dvech_sigma + (dvech_sigma_pattern * yp[[p]]$freq)
    }

    dvech_sigma
  }




# 3. Casewise scores

# 3a: casewise scores with respect to mu
lav_mvn_mi_sc_mu <- function(y = NULL,
                                         wt = NULL,
                                         mp = NULL,
                                         mu = NULL,
                                         sigma_1 = NULL,
                                         x_idx = integer(0L),
                                         sigma_inv = NULL,
                                         sinv_method = "eigen") {
  p_1 <- NCOL(y)
  mu <- as.numeric(mu)
  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }
  ny <- NROW(y)

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # subtract Mu
  yc <- t(t(y) - mu)

  # dmu per case
  dmu <- matrix(as.numeric(NA), ny, p_1)

  # for each pattern, compute sigma.inv
  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    case_idx <- mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv_1 <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
    } else {
      sigma_inv_1 <- sigma_inv
    }

    # compute dMu for all observations of this pattern
    dmu[case_idx, var_idx] <-
      yc[case_idx, var_idx, drop = FALSE] %*% sigma_inv_1
  }

  # weights
  if (!is.null(wt)) {
    dmu <- dmu * wt
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    dmu[, x_idx] <- 0
  }

  dmu
}

# 3b: casewise scores with respect to vech(Sigma)
lav_mvn_mi_sc_sigma <- function(y = NULL,
                                                 wt = NULL,
                                                 mp = NULL,
                                                 mu = NULL,
                                                 sigma_1 = NULL,
                                                 x_idx = integer(0L),
                                                 sigma_inv = NULL,
                                                 sinv_method = "eigen") {
  p_1 <- NCOL(y)
  mu <- as.numeric(mu)
  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }
  ny <- NROW(y)

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # for the tcrossprod
  idx1 <- lav_mat_vech_col_idx(p_1)
  idx2 <- lav_mat_vech_row_idx(p_1)

  # vech(Sigma.inv)
  i_sigma <- lav_mat_vech(sigma_inv)

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # subtract Mu
  yc <- t(t(y) - mu)

  # SC
  sc <- matrix(as.numeric(NA), nrow = ny, ncol = length(i_sigma))

  # for each pattern
  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # cases with this pattern
    case_idx <- mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv_1 <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
      tmp <- matrix(0, p_1, p_1)
      tmp[var_idx, var_idx] <- sigma_inv_1
      isigma <- lav_mat_vech(tmp)
    } else {
      sigma_inv_1 <- sigma_inv
      isigma <- i_sigma
    }

    # postmultiply these cases with sigma.inv
    yc[case_idx, var_idx] <- yc[case_idx, var_idx] %*% sigma_inv_1

    # tcrossprod
    sc[case_idx, ] <- yc[case_idx, idx1] * yc[case_idx, idx2]

    # subtract isigma from each row
    sc[case_idx, ] <- t(t(sc[case_idx, , drop = FALSE]) - isigma)
  }

  # adjust for vech
  sc[, lav_mat_diagh_idx(p_1)] <- sc[, lav_mat_diagh_idx(p_1)] / 2

  # weights
  if (!is.null(wt)) {
    sc <- sc * wt
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    sc[, lav_mat_vech_which_idx(n = p_1, idx = x_idx)] <- 0
  }

  sc
}

# 3c: casewise scores with respect to mu + vech(Sigma)
lav_mvn_mi_sc_mu_sigma <- function(y = NULL,
                                                    mp = NULL,
                                                    wt = NULL,
                                                    mu = NULL,
                                                    sigma_1 = NULL,
                                                    x_idx = integer(0L),
                                                    sigma_inv = NULL,
                                                    sinv_method = "eigen") {
  p_1 <- NCOL(y)
  mu <- as.numeric(mu)
  if (!is.null(wt)) {
    n <- sum(wt)
  } else {
    n <- NROW(y)
  }
  ny <- NROW(y)

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # for the tcrossprod
  idx1 <- lav_mat_vech_col_idx(p_1)
  idx2 <- lav_mat_vech_row_idx(p_1)

  # vech(Sigma.inv)
  i_sigma <- lav_mat_vech(sigma_inv)

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # subtract Mu
  yc <- t(t(y) - mu)

  # dmu per case
  dmu <- matrix(as.numeric(NA), ny, p_1)

  # SC
  sc <- matrix(as.numeric(NA), nrow = ny, ncol = length(i_sigma))

  # for each pattern, compute Yc %*% sigma.inv
  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # cases with this pattern
    case_idx <- mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv_1 <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
      tmp <- matrix(0, p_1, p_1)
      tmp[var_idx, var_idx] <- sigma_inv_1
      isigma <- lav_mat_vech(tmp)
    } else {
      sigma_inv_1 <- sigma_inv
      isigma <- i_sigma
    }

    # compute dMu for all observations of this pattern
    dmu[case_idx, var_idx] <-
      yc[case_idx, var_idx, drop = FALSE] %*% sigma_inv_1

    # postmultiply these cases with sigma.inv
    yc[case_idx, var_idx] <- yc[case_idx, var_idx] %*% sigma_inv_1

    # tcrossprod
    sc[case_idx, ] <- yc[case_idx, idx1] * yc[case_idx, idx2]

    # subtract isigma from each row
    sc[case_idx, ] <- t(t(sc[case_idx, , drop = FALSE]) - isigma)
  }

  # adjust for vech
  sc[, lav_mat_diagh_idx(p_1)] <- sc[, lav_mat_diagh_idx(p_1)] / 2

  out <- cbind(dmu, sc)

  # weights
  if (!is.null(wt)) {
    out <- out * wt
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    not_x <- lav_mat_vech_which_idx(n = p_1, idx = x_idx,
                                       #diagonal = !correlation,
                                       add_idx_at_start = TRUE)
    out[, not_x] <- 0
  }

  out
}


# 4) Hessian of logl
lav_mvn_mi_logl_hessian_data <- function(y = NULL,
                                                 mp = NULL,
                                                 wt = NULL,
                                                 mu = NULL,
                                                 sigma_1 = NULL,
                                                 x_idx = integer(0L),
                                                 sinv_method = "eigen",
                                                 sigma_inv = NULL) {
  # missing patterns
  yp <- lav_samp_mi_patterns(y = y, mp = mp, wt = wt)

  lav_mvnorm_missing_logl_hessian_samplestats(
    yp = yp, mu = mu,
    sigma_1 = sigma_1, x_idx = x_idx, sinv_method = sinv_method,
    sigma_inv = sigma_inv
  )
}

lav_mvnorm_missing_logl_hessian_samplestats <- # nolint
  function(yp = NULL,
           # wt not needed
           mu = NULL,
           sigma_1 = NULL,
           x_idx = integer(0L),
           sinv_method = "eigen",
           sigma_inv = NULL) {
    pat_n <- length(yp)
    p_1 <- length(yp[[1]]$var.idx)

    if (is.null(sigma_inv)) {
      sigma_inv <- lav_mat_sym_inverse(
        s = sigma_1, logdet = FALSE,
        sinv_method = sinv_method
      )
    }

    h11 <- matrix(0, p_1, p_1)
    h21 <- matrix(0, p_1 * (p_1 + 1) / 2, p_1)
    h22 <- matrix(0, p_1 * (p_1 + 1) / 2, p_1 * (p_1 + 1) / 2)

    # for each pattern, compute sigma.inv
    for (p in seq_len(pat_n)) {
      # observed variables
      var_idx <- yp[[p]]$var.idx
      pat_freq <- yp[[p]]$freq

      # missing values for this pattern
      na_idx <- which(!var_idx)

      # invert Sigma for this pattern
      if (length(na_idx) > 0L) {
        sigma_inv_1 <- lav_mat_sym_inverse_update(
          s_inv = sigma_inv,
          rm_idx = na_idx, logdet = FALSE
        )
      } else {
        sigma_inv_1 <- sigma_inv
      }

      s_inv <- matrix(0, p_1, p_1)
      s_inv[var_idx, var_idx] <- sigma_inv_1

      tmp21 <- matrix(0, p_1, 1)
      tmp21[var_idx, 1] <- sigma_inv_1 %*% (yp[[p]]$MY - mu[var_idx])

      w_tilde <- yp[[p]]$SY + tcrossprod(yp[[p]]$MY - mu[var_idx])
      aaa <- (sigma_inv_1 %*%
        (2 * w_tilde - sigma_1[var_idx, var_idx, drop = FALSE]) %*%
        sigma_inv_1)
      tmp22 <- matrix(0, p_1, p_1)
      tmp22[var_idx, var_idx] <- aaa

      i11 <- s_inv
      i21 <- lav_mat_dup_pre(tmp21 %x% s_inv)
      # if (lav_use_lavaanC()) {
      #   i22 <- lavaanC::m_kronecker_dup_pre_post(S.inv, tmp22, 0.5)
      # } else {
        i22 <- (1 / 2) * lav_mat_dup_pre_post(s_inv %x% tmp22)
      # }
      h11 <- h11 + pat_freq * i11
      h21 <- h21 + pat_freq * i21
      h22 <- h22 + pat_freq * i22
    }

    h12 <- t(h21)

    out <- -1 * rbind(
      cbind(h11, h12),
      cbind(h21, h22)
    )

    # fixed.x?
    if (length(x_idx) > 0L) {
      not_x <- lav_mat_vech_which_idx(n = p_1, idx = x_idx,
                                         #diagonal = !correlation,
                                         add_idx_at_start = TRUE)
      out[, not_x] <- 0
      out[not_x, ] <- 0
    }

    out
  }




# 5) Information

# 5a: expected unit information Mu and vech(Sigma)
#     (only useful under MCAR)
# (old term: Abeta, expected)
lav_mvn_mi_info_expected <- function(y = NULL,
                                                    mp = NULL,
                                                    wt = NULL,
                                                    mu = NULL, # unused
                                                    sigma_1 = NULL,
                                                    x_idx = integer(0L),
                                                    sigma_inv = NULL,
                                                    sinv_method = "eigen") {
  p_1 <- NCOL(y)

  if (is.null(sigma_inv)) {
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # N
  if (!is.null(wt)) {
    if (length(mp$empty.idx) > 0L) {
      n <- sum(wt) - sum(wt[mp$empty.idx])
    } else {
      n <- sum(wt)
    }
  } else {
    n <- sum(mp$freq) # removed empty cases!
  }

  i11 <- matrix(0, p_1, p_1)
  i22 <- matrix(0, p_1 * (p_1 + 1) / 2, p_1 * (p_1 + 1) / 2)

  # for each pattern, compute sigma.inv
  for (p in seq_len(mp$npatterns)) {
    # observed variables
    var_idx <- mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv_1 <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
    } else {
      sigma_inv_1 <- sigma_inv
    }

    s_inv <- matrix(0, p_1, p_1)
    s_inv[var_idx, var_idx] <- sigma_inv_1

    # if (lav_use_lavaanC()) {
    #   S2.inv <- lavaanC::m_kronecker_dup_pre_post(S.inv, multiplicator = 0.5)
    # } else {
      s2_inv <- 0.5 * lav_mat_dup_pre_post(s_inv %x% s_inv)
    # }

    if (!is.null(wt)) {
      freq <- sum(wt[mp$case.idx[[p]]])
    } else {
      freq <- mp$freq[p]
    }

    i11 <- i11 + freq * s_inv
    i22 <- i22 + freq * s2_inv
  }

  out <- lav_mat_bdiag(i11, i22) / n

  # fixed.x?
  if (length(x_idx) > 0L) {
    not_x <- lav_mat_vech_which_idx(n = p_1, idx = x_idx,
                                       # diagonal = !correlation,
                                       add_idx_at_start = TRUE)
    out[not_x, ] <- 0
    out[, not_x] <- 0
  }

  out
}

# 5b: unit observed information Mu and vech(Sigma) from raw data
# (old term: Abeta, observed)
lav_mvn_mi_info_observed_data <- function(y = NULL,
                                                         mp = NULL,
                                                         wt = NULL,
                                                         mu = NULL,
                                                         sigma_1 = NULL,
                                                         x_idx = integer(0L),
                                                         sinv_method = "eigen",
                                                         sigma_inv = NULL) {
  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # N
  if (!is.null(wt)) {
    if (length(mp$empty.idx) > 0L) {
      n <- sum(wt) - sum(wt[mp$empty.idx])
    } else {
      n <- sum(wt)
    }
  } else {
    n <- sum(mp$freq) # removed empty cases!
  }

  # observed information
  observed <- lav_mvn_mi_logl_hessian_data(
    y = y, mp = mp, wt = wt,
    mu = mu, sigma_1 = sigma_1, x_idx = x_idx,
    sinv_method = sinv_method, sigma_inv = sigma_inv
  )

  -observed / n
}

# 5b-bis: unit observed information Mu and vech(Sigma) from samplestats
lav_mvnorm_missing_information_observed_samplestats <- # nolint
  function(yp = NULL,
           # wt not needed
           mu = NULL,
           sigma_1 = NULL,
           x_idx = integer(0L),
           sinv_method = "eigen",
           sigma_inv = NULL) {
    n <- sum(sapply(yp, "[[", "freq")) # implicitly: removed empty cases!

    # observed information
    observed <- lav_mvnorm_missing_logl_hessian_samplestats(
      yp = yp, mu = mu,
      sigma_1 = sigma_1, x_idx = x_idx, sinv_method = sinv_method,
      sigma_inv = sigma_inv
    )

    -observed / n
  }

# 5c: unit first-order information Mu and vech(Sigma) from raw data
# (old term: Bbeta)
lav_mvn_mi_info_firstorder <- function(y = NULL,
                                                      mp = NULL,
                                                      wt = NULL,
                                                      cluster_idx = NULL,
                                                      mu = NULL,
                                                      sigma_1 = NULL,
                                                      x_idx = integer(0L),
                                                      sinv_method = "eigen",
                                                      sigma_inv = NULL) {
  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # N
  if (!is.null(wt)) {
    if (length(mp$empty.idx) > 0L) {
      n <- sum(wt) - sum(wt[mp$empty.idx])
    } else {
      n <- sum(wt)
    }
  } else {
    n <- sum(mp$freq) # removed empty cases!
  }

  sc <- lav_mvn_mi_sc_mu_sigma(
    y = y, mp = mp, wt = wt,
    mu = mu, sigma_1 = sigma_1, x_idx = x_idx, sinv_method = sinv_method,
    sigma_inv = sigma_inv
  )

  # handle clustering
  if (!is.null(cluster_idx)) {
    # take the sum within each cluster
    sc <- rowsum(sc, group = cluster_idx, reorder = FALSE, na.rm = TRUE)

    # lower bias is number of clusters is not very high
    n_c <- nrow(sc)
    correction_factor <- n_c / (n_c - 1)
    sc <- sc * sqrt(correction_factor)
  }

  lav_mat_crossprod(sc) / n
}

# 5d: both unit first-order information and expected/observed information
#     from raw data, in one go for efficiency
lav_mvn_mi_info_both <- function(y = NULL,
                                                mp = NULL,
                                                wt = NULL,
                                                cluster_idx = NULL,
                                                mu = NULL,
                                                sigma_1 = NULL,
                                                x_idx = integer(0L),
                                                sinv_method = "eigen",
                                                sigma_inv = NULL,
                                                information = "observed") {
  p_1 <- NCOL(y)
  mu <- as.numeric(mu)

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # for the tcrossprod
  idx1 <- lav_mat_vech_col_idx(p_1)
  idx2 <- lav_mat_vech_row_idx(p_1)

  # vech(Sigma.inv)
  i_sigma <- lav_mat_vech(sigma_inv)

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  if (information == "observed") {
    yp <- lav_samp_mi_patterns(y = y, mp = mp, wt = wt)
  }

  # N
  if (!is.null(wt)) {
    if (length(mp$empty.idx) > 0L) {
      n <- sum(wt) - sum(wt[mp$empty.idx])
    } else {
      n <- sum(wt)
    }
  } else {
    n <- sum(mp$freq) # removed empty cases!
  }

  # subtract Mu
  yc <- t(t(y) - mu)

  # dmu per case
  dmu <- matrix(as.numeric(NA), nrow = NROW(y), ncol = p_1)

  # SC
  sc <- matrix(as.numeric(NA), nrow = NROW(y), ncol = length(i_sigma))

  # expected/observed information
  i11_1 <- matrix(0, p_1, p_1)
  i22_1 <- matrix(0, p_1 * (p_1 + 1) / 2, p_1 * (p_1 + 1) / 2)
  if (information == "observed") {
    i21_1 <- matrix(0, p_1 * (p_1 + 1) / 2, p_1)
  }

  # for each pattern, compute Yc %*% sigma.inv
  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # cases with this pattern
    case_idx <- mp$case.idx[[p]]

    # invert Sigma for this pattern
    if (length(na_idx) > 0L) {
      sigma_inv_1 <- lav_mat_sym_inverse_update(
        s_inv = sigma_inv,
        rm_idx = na_idx, logdet = FALSE
      )
      tmp <- matrix(0, p_1, p_1)
      tmp[var_idx, var_idx] <- sigma_inv_1
      isigma <- lav_mat_vech(tmp)
    } else {
      sigma_inv_1 <- sigma_inv
      isigma <- i_sigma
    }

    # information
    s_inv <- matrix(0, p_1, p_1)
    s_inv[var_idx, var_idx] <- sigma_inv_1

    if (!is.null(wt)) {
      freq <- sum(wt[case_idx])
    } else {
      freq <- mp$freq[p]
    }

    if (information == "expected") {
      # if (lav_use_lavaanC()) {
      #   S2.inv <-
      #         lavaanC::m_kronecker_dup_pre_post(S.inv, multiplicator = 0.5)
      # } else {
        s2_inv <- 0.5 * lav_mat_dup_pre_post(s_inv %x% s_inv)
      # }

      i11_1 <- i11_1 + freq * s_inv
      i22_1 <- i22_1 + freq * s2_inv
    } else {
      pat_freq <- yp[[p]]$freq

      tmp21 <- matrix(0, p_1, 1)
      tmp21[var_idx, 1] <- sigma_inv_1 %*% (yp[[p]]$MY - mu[var_idx])

      w_tilde <- yp[[p]]$SY + tcrossprod(yp[[p]]$MY - mu[var_idx])
      aaa <- (sigma_inv_1 %*%
        (2 * w_tilde - sigma_1[var_idx, var_idx, drop = FALSE]) %*%
        sigma_inv_1)
      tmp22 <- matrix(0, p_1, p_1)
      tmp22[var_idx, var_idx] <- aaa

      i11 <- s_inv
      i21 <- lav_mat_dup_pre(tmp21 %x% s_inv)
      # if (lav_use_lavaanC()) {
      #   i22 <- lavaanC::m_kronecker_dup_pre_post(S.inv, tmp22, 0.5)
      # } else {
        i22 <- (1 / 2) * lav_mat_dup_pre_post(s_inv %x% tmp22)
      # }

      i11_1 <- i11_1 + pat_freq * i11
      i21_1 <- i21_1 + pat_freq * i21
      i22_1 <- i22_1 + pat_freq * i22
    }

    # compute dMu for all observations of this pattern
    dmu[case_idx, var_idx] <-
      yc[case_idx, var_idx, drop = FALSE] %*% sigma_inv_1

    # postmultiply these cases with sigma.inv
    yc[case_idx, var_idx] <- yc[case_idx, var_idx] %*% sigma_inv_1

    # tcrossprod
    sc[case_idx, ] <- yc[case_idx, idx1] * yc[case_idx, idx2]

    # subtract isigma from each row
    sc[case_idx, ] <- t(t(sc[case_idx, , drop = FALSE]) - isigma)
  }

  # adjust for vech
  sc[, lav_mat_diagh_idx(p_1)] <- sc[, lav_mat_diagh_idx(p_1)] / 2

  # add dmu
  sc <- cbind(dmu, sc)

  # weights
  if (!is.null(wt)) {
    sc <- sc * wt
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    not_x <- lav_mat_vech_which_idx(n = p_1, idx = x_idx,
                                       #diagonal = !correlation,
                                       add_idx_at_start = TRUE)
    sc[, not_x] <- 0
  }

  # handle clustering
  if (!is.null(cluster_idx)) {
    # take the sum within each cluster
    sc <- rowsum(sc, group = cluster_idx, reorder = FALSE, na.rm = TRUE)

    # lower bias is number of clusters is not very high
    n_c <- nrow(sc)
    correction_factor <- n_c / (n_c - 1)
    sc <- sc * sqrt(correction_factor)
  }

  # first order information
  bbeta <- lav_mat_crossprod(sc) / n

  # expected/observed information
  if (information == "expected") {
    abeta <- lav_mat_bdiag(i11_1, i22_1) / n
  } else {
    abeta <- rbind(
      cbind(i11_1, t(i21_1)),
      cbind(i21_1, i22_1)
    ) / n
  }

  # fixed.x?
  if (length(x_idx) > 0L) {
    not_x <- lav_mat_vech_which_idx(n = p_1, idx = x_idx,
                                       # diagonal = !correlation,
                                       add_idx_at_start = TRUE)
    abeta[not_x, ] <- 0
    abeta[, not_x] <- 0
  }

  list(Abeta = abeta, Bbeta = bbeta)
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
#        Mp <- lav_data_mi_patterns(Y)
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
#        sigma2 <- 2 * lav_mat_dup_ginv_pre_post(sigma %x% sigma)
#
#        I11 <- I11 + Mp$freq[p] * sigma
#        I22 <- I22 + Mp$freq[p] * sigma2
#    }
#
#    lav_mat_bdiag(I11, I22)/N
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
lav_mvn_mi_impute_pattern <- function(y = NULL,
                                              mp = NULL,
                                              mu = NULL,
                                              sigma_1 = NULL,
                                              sigma_inv = NULL,
                                              sinv_method = "eigen") {
  mu <- as.numeric(mu)

  # complete data
  y_complete <- y

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # subtract Mu
  yc <- t(t(y) - mu)

  # fill in data per pattern
  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # if complete, nothing to do
    if (all(var_idx)) {
      next
    }

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # extract observed data for these (centered) cases
    oc <- yc[mp$case.idx[[p]], mp$pat[p, ], drop = FALSE]

    # invert Sigma (Sigma_22, observed part only) for this pattern
    sigma_22_inv <- try(
      lav_mat_sym_inverse_update(
        s_inv = sigma_inv, rm_idx = na_idx, logdet = FALSE
      ),
      silent = TRUE
    )
    if (inherits(sigma_22_inv, "try-error")) {
      lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
    }

    # estimate missing values in this pattern
    sigma_12 <- sigma_1[!var_idx, var_idx, drop = FALSE]
    y_missing <- t(sigma_12 %*% sigma_22_inv %*% t(oc) + mu[!var_idx])

    # complete data for this pattern
    y_complete[mp$case.idx[[p]], !var_idx] <- y_missing
  }

  y_complete
}


# E-step: expectations of sum, sum of squares, sum of crossproducts
# plus correction
lav_mvn_mi_estep <- function(y = NULL,
                                     mp = NULL,
                                     wt = NULL,
                                     mu = NULL,
                                     sigma_1 = NULL,
                                     sigma_inv = NULL,
                                     sinv_method = "eigen") {
  p_1 <- NCOL(y)
  mu <- as.numeric(mu)

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  if (is.null(sigma_inv)) {
    # invert Sigma
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = FALSE,
      sinv_method = sinv_method
    )
  }

  # T1, T2
  t1 <- numeric(p_1)
  t2 <- matrix(0, p_1, p_1)

  # update T1 and T2 per pattern
  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # extract observed data
    o <- y[mp$case.idx[[p]], mp$pat[p, ], drop = FALSE]

    # if complete, just compute first and second moments
    if (all(var_idx)) {
      if (!is.null(wt)) {
        wt_1 <- wt[mp$case.idx[[p]]]
        t1 <- t1 + colSums(wt_1 * o)
        t2 <- t2 + crossprod(sqrt(wt_1) * o)
      } else {
        # complete pattern
        t1 <- t1 + colSums(o)
        t2 <- t2 + crossprod(o)
      }
      next
    }

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # partition Sigma (1=missing, 2=complete)
    sigma_11 <- sigma_1[!var_idx, !var_idx, drop = FALSE]
    sigma_12 <- sigma_1[!var_idx, var_idx, drop = FALSE]
    sigma_21 <- sigma_1[var_idx, !var_idx, drop = FALSE]

    # invert Sigma (Sigma_22, observed part only) for this pattern
    sigma_22_inv <- try(
      lav_mat_sym_inverse_update(
        s_inv = sigma_inv, rm_idx = na_idx, logdet = FALSE
      ),
      silent = TRUE
    )
    if (inherits(sigma_22_inv, "try-error")) {
      lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
    }

    # estimate missing values in this pattern
    oc <- t(t(o) - mu[var_idx])
    y_missing <- t(sigma_12 %*% sigma_22_inv %*% t(oc) + mu[!var_idx])

    # complete data for this pattern
    y_complete <- matrix(0, mp$freq[[p]], p_1)
    y_complete[, var_idx] <- o
    y_complete[, !var_idx] <- y_missing

    if (!is.null(wt)) {
      wt_1 <- wt[mp$case.idx[[p]]]
      t1_pat <- colSums(wt_1 * y_complete)
      t2_pat <- crossprod(sqrt(wt_1) * y_complete)
    } else {
      # 1. SUM `completed' pattern
      t1_pat <- colSums(y_complete)

      # 2. CROSSPROD `completed' pattern
      t2_pat <- crossprod(y_complete)
    }

    # correction for missing cells: conditional covariances
    t2_p11 <- sigma_11 - (sigma_12 %*% sigma_22_inv %*% sigma_21)
    if (!is.null(wt)) {
      t2_pat[!var_idx, !var_idx] <-
        t2_pat[!var_idx, !var_idx] + (t2_p11 * sum(wt_1))
    } else {
      t2_pat[!var_idx, !var_idx] <-
        t2_pat[!var_idx, !var_idx] + (t2_p11 * mp$freq[[p]])
    }

    # accumulate
    t1 <- t1 + t1_pat
    t2 <- t2 + t2_pat
  }

  list(T1 = t1, T2 = t2)
}
