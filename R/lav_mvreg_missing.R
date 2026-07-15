# the multivariate linear model + missing (y) values
#  ('FIML' for the conditional.x case)
#
# the model is y_i | x_i ~ N(Beta' X1_i, res.cov), where X1_i = c(1, x_i);
# the x (exogenous) values are assumed to be complete; the y values may
# contain missing values, organized in missing patterns; each case
# contributes the loglikelihood of its observed y values only:
#   log f(y_obs,i | x_i) with mean (Beta' X1_i)[obs] and
#   covariance res.cov[obs, obs]

# 0) extended pattern statistics (sufficient statistics per pattern)
# 1) loglikelihood (from raw data, or pattern statistics)
# 2) derivatives with respect to Beta, res.cov, vech(res.cov)
# 3) casewise scores with respect to Beta + vech(res.cov)
# 4) hessian of Beta + vech(res.cov)
# 5) (unit) information of Beta + vech(res.cov)
#    5a: (unit)    expected information
#    5b: (unit)    observed information
#    5c: (unit) first.order information

# the parameter order is always: vec(Beta) | vech(res.cov), where
# vec(Beta) is the column-major vec of the (1+nexo) x nvar matrix Beta,
# i.e. per y variable: intercept, slope_1, slope_2, ...; this matches
# both lav_mvreg_sc_beta() and the conditional.x DELTA rows

# YR/CC 13 Jul 2026: first version

# 0. extended per-pattern sufficient statistics
#
#    because the conditional mean is linear in x, the crossproducts
#    WXX = sum_i X1_i X1_i', WXY = sum_i X1_i y_obs,i' and
#    WYY = sum_i y_obs,i y_obs,i' (sums over the cases within a pattern,
#    weighted if wt is given) are sufficient for the loglikelihood, the
#    gradient, the hessian and the expected information; only the
#    casewise scores need the raw data
#
#    the result extends the lav_samp_mi_patterns() output: each element
#    also contains the usual SY, MY, var.idx and freq fields
lav_samp_mi_patterns_res <- function(y = NULL,
                                     exo = NULL, # no intercept
                                     mp = NULL,
                                     wt = NULL) {
  y <- as.matrix(unname(y))
  x1 <- cbind(1, as.matrix(unname(exo)))

  # missing patterns (y side only; x is complete)
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y, sort_freq = FALSE, coverage = FALSE)
  }

  # the usual pattern statistics (SY, MY, var.idx, freq)
  yp <- lav_samp_mi_patterns(y = y, mp = mp, wt = wt)

  # add the crossproduct statistics
  for (p in seq_len(mp$npatterns)) {
    case_idx <- mp$case.idx[[p]]
    var_idx <- mp$pat[p, ]

    x1_p <- x1[case_idx, , drop = FALSE]
    y_p <- y[case_idx, var_idx, drop = FALSE]

    if (!is.null(wt)) {
      wx1_p <- x1_p * wt[case_idx]
      yp[[p]]$WXX <- crossprod(wx1_p, x1_p)
      yp[[p]]$WXY <- crossprod(wx1_p, y_p)
      yp[[p]]$WYY <- crossprod(y_p * wt[case_idx], y_p)
    } else {
      yp[[p]]$WXX <- crossprod(x1_p)
      yp[[p]]$WXY <- crossprod(x1_p, y_p)
      yp[[p]]$WYY <- crossprod(y_p)
    }
  }

  yp
}

# small helper: per-pattern beta kernel and total residual crossprod
#
#   C.p      = WXY - WXX %*% Beta[, obs]         (the beta score kernel)
#   W.tilde  = sum_i r_i r_i'  (total, NOT divided by freq), where
#              r_i = y_obs,i - Beta[, obs]' X1_i
lav_mvreg_mi_pat_kernel <- function(yp_1 = NULL, beta_1 = NULL) {
  b_o <- beta_1[, yp_1$var.idx, drop = FALSE]
  xyb <- yp_1$WXY - yp_1$WXX %*% b_o
  w_tilde <- (yp_1$WYY - crossprod(yp_1$WXY, b_o) - crossprod(b_o, yp_1$WXY)
              + crossprod(b_o, yp_1$WXX %*% b_o))
  list(c_p = xyb, w_tilde = w_tilde)
}

# 1. loglikelihood

# 1a: input is raw data (casewise = TRUE returns casewise loglikelihoods;
#     rows of empty cases remain NA)
lav_mvreg_mi_loglik_data <- function(y = NULL,
                                     exo = NULL, # no intercept
                                     mp = NULL,
                                     wt = NULL,
                                     beta_1 = NULL,
                                     res_int = NULL,
                                     res_slopes = NULL,
                                     res_cov = NULL,
                                     casewise = FALSE,
                                     sinv_method = "eigen",
                                     res_cov_inv = NULL,
                                     log2pi = TRUE,
                                     minus_two = FALSE) {
  y <- as.matrix(unname(y))
  x1 <- cbind(1, as.matrix(unname(exo)))
  ny <- NROW(y)
  log_2pi <- log(2 * pi)

  # construct model-implied Beta
  beta_1 <- lav_mvreg_beta1(beta_1, res_int, res_slopes)

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y, sort_freq = FALSE, coverage = FALSE)
  }

  # global inverse + logdet
  res_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = res_cov, sigma_inv = res_cov_inv,
    sinv_method = sinv_method, logdet = TRUE
  )
  res_cov_logdet <- attr(res_cov_inv, "logdet")

  # residuals (NAs propagate)
  res <- y - x1 %*% beta_1

  # DIST/logdet per case
  dist_1 <- logdet <- p_log_2pi <- rep(as.numeric(NA), ny)

  for (p in seq_len(mp$npatterns)) {
    # observed values for this pattern
    var_idx <- mp$pat[p, ]

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # cases with this pattern
    case_idx <- mp$case.idx[[p]]

    # constant
    p_log_2pi[case_idx] <- sum(var_idx) * log_2pi

    # invert res.cov for this pattern
    sigma_inv <- lav_mvn_mi_pattern_inv(
      sigma_inv = res_cov_inv, na_idx = na_idx,
      logdet = TRUE, sigma_logdet = res_cov_logdet
    )
    logdet[case_idx] <- attr(sigma_inv, "logdet")

    if (mp$freq[p] == 1L) {
      dist_1[case_idx] <- sum(sigma_inv *
        crossprod(res[case_idx, var_idx, drop = FALSE]))
    } else {
      dist_1[case_idx] <-
        rowSums(res[case_idx, var_idx, drop = FALSE] %*% sigma_inv *
          res[case_idx, var_idx, drop = FALSE])
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

  if (casewise) {
    loglik <- llik
  } else {
    loglik <- sum(llik, na.rm = TRUE)
  }

  loglik
}

# 1b: input are the extended pattern statistics (see
#     lav_samp_mi_patterns_res)
lav_mvreg_mi_loglik_samp <- function(yp = NULL,
                                     beta_1 = NULL,
                                     res_int = NULL,
                                     res_slopes = NULL,
                                     res_cov = NULL,
                                     sinv_method = "eigen",
                                     res_cov_inv = NULL,
                                     log2pi = TRUE,
                                     minus_two = FALSE) {
  log_2pi <- log(2 * pi)
  pat_n <- length(yp)

  # construct model-implied Beta
  beta_1 <- lav_mvreg_beta1(beta_1, res_int, res_slopes)

  # global inverse + logdet
  res_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = res_cov, sigma_inv = res_cov_inv,
    sinv_method = sinv_method, logdet = TRUE
  )
  res_cov_logdet <- attr(res_cov_inv, "logdet")

  # DIST/logdet per pattern
  dist_1 <- logdet <- p_log_2pi <- numeric(pat_n)

  for (p in seq_len(pat_n)) {
    # observed variables for this pattern
    var_idx <- yp[[p]]$var.idx

    # missing values for this pattern
    na_idx <- which(!var_idx)

    # constant
    p_log_2pi[p] <- sum(var_idx) * log_2pi * yp[[p]]$freq

    # invert res.cov for this pattern
    sigma_inv <- lav_mvn_mi_pattern_inv(
      sigma_inv = res_cov_inv, na_idx = na_idx,
      logdet = TRUE, sigma_logdet = res_cov_logdet
    )
    logdet[p] <- attr(sigma_inv, "logdet") * yp[[p]]$freq

    # total residual crossprod for this pattern
    kern <- lav_mvreg_mi_pat_kernel(yp_1 = yp[[p]], beta_1 = beta_1)
    dist_1[p] <- sum(sigma_inv * kern$w_tilde)
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

  loglik
}


# 2. Derivatives

# 2a: derivative logl with respect to Beta (=intercepts and slopes),
#     returned as vec(Beta) (column-major of the (1+nexo) x nvar matrix)
lav_mvreg_mi_dlogl_dbeta_samp <- function(yp = NULL,
                                          beta_1 = NULL,
                                          res_int = NULL,
                                          res_slopes = NULL,
                                          res_cov = NULL,
                                          sinv_method = "eigen",
                                          res_cov_inv = NULL) {
  # construct model-implied Beta
  beta_1 <- lav_mvreg_beta1(beta_1, res_int, res_slopes)
  p_1 <- NCOL(beta_1)

  # global inverse
  res_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = res_cov, sigma_inv = res_cov_inv,
    sinv_method = sinv_method
  )

  dbeta <- matrix(0, NROW(beta_1), p_1)

  for (p in seq_len(length(yp))) {
    var_idx <- yp[[p]]$var.idx

    # invert res.cov for this pattern
    sigma_inv <- lav_mvn_mi_pattern_inv(
      sigma_inv = res_cov_inv, na_idx = which(!var_idx)
    )

    kern <- lav_mvreg_mi_pat_kernel(yp_1 = yp[[p]], beta_1 = beta_1)
    dbeta[, var_idx] <- dbeta[, var_idx] + kern$c_p %*% sigma_inv
  }

  as.numeric(dbeta)
}

# 2b: derivative logl with respect to res.cov (full matrix, ignoring
#     symmetry)
lav_mvreg_mi_dlogl_drescov_samp <- function(yp = NULL,
                                            beta_1 = NULL,
                                            res_int = NULL,
                                            res_slopes = NULL,
                                            res_cov = NULL,
                                            sinv_method = "eigen",
                                            res_cov_inv = NULL) {
  # construct model-implied Beta
  beta_1 <- lav_mvreg_beta1(beta_1, res_int, res_slopes)
  p_1 <- NCOL(beta_1)

  # global inverse
  res_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = res_cov, sigma_inv = res_cov_inv,
    sinv_method = sinv_method
  )

  dres_cov <- matrix(0, p_1, p_1)

  for (p in seq_len(length(yp))) {
    var_idx <- yp[[p]]$var.idx
    pat_freq <- yp[[p]]$freq

    # invert res.cov for this pattern
    sigma_inv <- lav_mvn_mi_pattern_inv(
      sigma_inv = res_cov_inv, na_idx = which(!var_idx)
    )

    kern <- lav_mvreg_mi_pat_kernel(yp_1 = yp[[p]], beta_1 = beta_1)
    dres_cov[var_idx, var_idx] <- (dres_cov[var_idx, var_idx]
      - (pat_freq / 2) * (sigma_inv -
        (sigma_inv %*% (kern$w_tilde / pat_freq) %*% sigma_inv)))
  }

  dres_cov
}

# 2c: derivative logl with respect to vech(res.cov)
lav_mvreg_mi_dlogl_dvechrescov_samp <- function(yp = NULL,
                                                beta_1 = NULL,
                                                res_int = NULL,
                                                res_slopes = NULL,
                                                res_cov = NULL,
                                                sinv_method = "eigen",
                                                res_cov_inv = NULL) {
  dres_cov <- lav_mvreg_mi_dlogl_drescov_samp(
    yp = yp, beta_1 = beta_1,
    res_int = res_int, res_slopes = res_slopes, res_cov = res_cov,
    sinv_method = sinv_method, res_cov_inv = res_cov_inv
  )

  as.numeric(lav_mat_dup_pre(as.matrix(lav_mat_vec(dres_cov))))
}


# 3. Casewise scores

# 3a: casewise scores with respect to Beta + vech(res.cov)
#     column order: Y1_int, Y1_x1, ... | Y2_int, Y2_x1, ... | vech(res.cov)
#     NAs are left in place (empty cases, and the vech(res.cov) columns
#     that involve an unobserved variable); callers should zero them out
lav_mvreg_mi_sc_beta_sigma <- function(y = NULL,
                                       exo = NULL, # no intercept
                                       mp = NULL,
                                       wt = NULL,
                                       beta_1 = NULL,
                                       res_int = NULL,
                                       res_slopes = NULL,
                                       res_cov = NULL,
                                       sinv_method = "eigen",
                                       res_cov_inv = NULL) {
  y <- as.matrix(unname(y))
  q_1 <- NCOL(y)
  x1 <- cbind(1, as.matrix(unname(exo)))
  p <- NCOL(x1)

  # construct model-implied Beta
  beta_1 <- lav_mvreg_beta1(beta_1, res_int, res_slopes)

  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y, sort_freq = FALSE, coverage = FALSE)
  }

  # global inverse
  res_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = res_cov, sigma_inv = res_cov_inv,
    sinv_method = sinv_method
  )

  # residuals (NAs propagate)
  res <- y - x1 %*% beta_1

  # score engine (res %*% res.cov.inv per pattern + vech(res.cov) scores)
  eng <- lav_mvn_mi_sc_engine(
    y = res, mp = mp, mu = numeric(q_1),
    sigma_inv = res_cov_inv
  )

  sc_beta <- x1[, rep(1:p, times = q_1), drop = FALSE] *
    eng$dmu[, rep(1:q_1, each = p), drop = FALSE]

  out <- cbind(sc_beta, eng$sc)

  # weights
  if (!is.null(wt)) {
    out <- out * wt
  }

  out
}


# 4. hessian of logl

# 4a. hessian logl Beta and vech(res.cov) from the extended pattern
#     statistics
lav_mvreg_mi_logl_hessian_samp <- function(yp = NULL,
                                           beta_1 = NULL,
                                           res_int = NULL,
                                           res_slopes = NULL,
                                           res_cov = NULL,
                                           sinv_method = "eigen",
                                           res_cov_inv = NULL) {
  # construct model-implied Beta
  beta_1 <- lav_mvreg_beta1(beta_1, res_int, res_slopes)
  p_1 <- NCOL(beta_1)
  nx1 <- NROW(beta_1) # 1 + nexo
  pstar <- p_1 * (p_1 + 1) / 2

  # global inverse
  res_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = res_cov, sigma_inv = res_cov_inv,
    sinv_method = sinv_method
  )

  h11 <- matrix(0, nx1 * p_1, nx1 * p_1)
  h21 <- matrix(0, pstar, nx1 * p_1)
  h22 <- matrix(0, pstar, pstar)

  for (p in seq_len(length(yp))) {
    var_idx <- yp[[p]]$var.idx
    pat_freq <- yp[[p]]$freq

    # invert res.cov for this pattern
    sigma_inv_1 <- lav_mvn_mi_pattern_inv(
      sigma_inv = res_cov_inv, na_idx = which(!var_idx)
    )

    # zero-padded pattern inverse
    s_inv <- matrix(0, p_1, p_1)
    s_inv[var_idx, var_idx] <- sigma_inv_1

    kern <- lav_mvreg_mi_pat_kernel(yp_1 = yp[[p]], beta_1 = beta_1)

    # beta-beta block: vec(Beta) positions of the observed y variables
    obs_idx <- which(var_idx)
    idx_o <- as.vector(outer(seq_len(nx1), (obs_idx - 1L) * nx1, "+"))
    h11[idx_o, idx_o] <- (h11[idx_o, idx_o]
      + (sigma_inv_1 %x% yp[[p]]$WXX))

    # sigma-beta block: zero-padded C.p (per-pattern beta kernel)
    c_bar <- matrix(0, nx1, p_1)
    c_bar[, var_idx] <- kern$c_p
    h21 <- h21 + lav_mat_dup_pre(s_inv %x% (s_inv %*% t(c_bar)))

    # sigma-sigma block
    aaa <- (sigma_inv_1 %*%
      (2 * (kern$w_tilde / pat_freq) -
        res_cov[var_idx, var_idx, drop = FALSE]) %*%
      sigma_inv_1)
    tmp22 <- matrix(0, p_1, p_1)
    tmp22[var_idx, var_idx] <- aaa
    h22 <- h22 + pat_freq * lav_mvn_kron_dup_half(s_inv, tmp22)
  }

  h12 <- t(h21)

  -1 * rbind(
    cbind(h11, h12),
    cbind(h21, h22)
  )
}

# 4b. hessian logl Beta and vech(res.cov) from raw data
lav_mvreg_mi_logl_hessian_data <- function(y = NULL,
                                           exo = NULL, # no intercept
                                           mp = NULL,
                                           wt = NULL,
                                           beta_1 = NULL,
                                           res_int = NULL,
                                           res_slopes = NULL,
                                           res_cov = NULL,
                                           sinv_method = "eigen",
                                           res_cov_inv = NULL) {
  yp <- lav_samp_mi_patterns_res(y = y, exo = exo, mp = mp, wt = wt)

  lav_mvreg_mi_logl_hessian_samp(
    yp = yp, beta_1 = beta_1,
    res_int = res_int, res_slopes = res_slopes, res_cov = res_cov,
    sinv_method = sinv_method, res_cov_inv = res_cov_inv
  )
}


# 5. Information

# 5a: expected unit information Beta and vech(res.cov)
#     (expected given the observed missing patterns; as in the
#     unconditional case, this is only strictly valid under MCAR)
lav_mvreg_mi_info_expected <- function(yp = NULL,
                                       res_cov = NULL,
                                       res_cov_inv = NULL,
                                       sinv_method = "eigen") {
  p_1 <- NCOL(res_cov)
  nx1 <- NROW(yp[[1]]$WXX) # 1 + nexo
  pstar <- p_1 * (p_1 + 1) / 2

  # global inverse
  res_cov_inv <- lav_mvn_sigma_inv(
    sigma_1 = res_cov, sigma_inv = res_cov_inv,
    sinv_method = sinv_method
  )

  # N
  n <- sum(sapply(yp, "[[", "freq")) # implicitly: removed empty cases!

  i11 <- matrix(0, nx1 * p_1, nx1 * p_1)
  i22 <- matrix(0, pstar, pstar)

  for (p in seq_len(length(yp))) {
    var_idx <- yp[[p]]$var.idx
    pat_freq <- yp[[p]]$freq

    # invert res.cov for this pattern
    sigma_inv_1 <- lav_mvn_mi_pattern_inv(
      sigma_inv = res_cov_inv, na_idx = which(!var_idx)
    )

    # beta block
    obs_idx <- which(var_idx)
    idx_o <- as.vector(outer(seq_len(nx1), (obs_idx - 1L) * nx1, "+"))
    i11[idx_o, idx_o] <- (i11[idx_o, idx_o]
      + (sigma_inv_1 %x% yp[[p]]$WXX))

    # sigma block
    s_inv <- matrix(0, p_1, p_1)
    s_inv[var_idx, var_idx] <- sigma_inv_1
    i22 <- i22 + pat_freq * lav_mvn_kron_dup_half(s_inv)
  }

  lav_mat_bdiag(i11, i22) / n
}

# 5b: unit observed information Beta and vech(res.cov) from the extended
#     pattern statistics
lav_mvreg_mi_information_observed_samplestats <- # nolint
  function(yp = NULL,
           beta_1 = NULL,
           res_int = NULL,
           res_slopes = NULL,
           res_cov = NULL,
           sinv_method = "eigen",
           res_cov_inv = NULL) {
    n <- sum(sapply(yp, "[[", "freq")) # implicitly: removed empty cases!

    # observed information
    observed <- lav_mvreg_mi_logl_hessian_samp(
      yp = yp, beta_1 = beta_1,
      res_int = res_int, res_slopes = res_slopes, res_cov = res_cov,
      sinv_method = sinv_method, res_cov_inv = res_cov_inv
    )

    -observed / n
  }

# 5c: unit first-order information Beta and vech(res.cov) from raw data
lav_mvreg_mi_info_firstorder <- function(y = NULL,
                                         exo = NULL, # no intercept
                                         mp = NULL,
                                         wt = NULL,
                                         cluster_idx = NULL,
                                         beta_1 = NULL,
                                         res_int = NULL,
                                         res_slopes = NULL,
                                         res_cov = NULL,
                                         sinv_method = "eigen",
                                         res_cov_inv = NULL) {
  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y, sort_freq = FALSE, coverage = FALSE)
  }

  # N
  n <- lav_mvn_mi_nobs(mp = mp, wt = wt)

  sc <- lav_mvreg_mi_sc_beta_sigma(
    y = y, exo = exo, mp = mp, wt = wt,
    beta_1 = beta_1, res_int = res_int, res_slopes = res_slopes,
    res_cov = res_cov, sinv_method = sinv_method,
    res_cov_inv = res_cov_inv
  )

  # handle clustering
  if (!is.null(cluster_idx)) {
    # take the sum within each cluster
    sc <- rowsum(sc, group = cluster_idx, reorder = FALSE, na.rm = TRUE)

    # lower bias if number of clusters is not very high
    n_c <- nrow(sc)
    correction_factor <- n_c / (n_c - 1)
    sc <- sc * sqrt(correction_factor)
  } else {
    sc[is.na(sc)] <- 0
  }

  lav_mat_crossprod(sc) / n
}


# 6. sandwich ACOV of the (saturated) Beta + vech(res.cov) estimates:
#    A^{-1} B A^{-1} (Savalei & Falk, 2014, in the conditional metric),
#    used for se = "robust.two.stage" + conditional.x
lav_mvreg_mi_h1_omega_sw <- function(y = NULL,
                                     exo = NULL, # no intercept
                                     mp = NULL,
                                     wt = NULL,
                                     cluster_idx = NULL,
                                     yp = NULL,
                                     sinv_method = "eigen",
                                     res_int = NULL,
                                     res_slopes = NULL,
                                     res_cov = NULL,
                                     res_cov_inv = NULL,
                                     information = "observed") {
  # missing patterns
  if (is.null(mp)) {
    mp <- lav_data_mi_patterns(y)
  }

  # extended pattern statistics
  if (is.null(yp)) {
    yp <- lav_samp_mi_patterns_res(y = y, exo = exo, mp = mp, wt = wt)
  }

  # bread
  if (information == "expected") {
    a <- lav_mvreg_mi_info_expected(
      yp = yp, res_cov = res_cov,
      res_cov_inv = res_cov_inv, sinv_method = sinv_method
    )
  } else {
    a <- lav_mvreg_mi_information_observed_samplestats(
      yp = yp, res_int = res_int, res_slopes = res_slopes,
      res_cov = res_cov, res_cov_inv = res_cov_inv,
      sinv_method = sinv_method
    )
  }
  a_inv <- lav_mat_sym_inverse(
    s = a, logdet = FALSE,
    sinv_method = sinv_method
  )

  # meat
  m_b <- lav_mvreg_mi_info_firstorder(
    y = y, exo = exo, mp = mp, wt = wt,
    cluster_idx = cluster_idx,
    res_int = res_int, res_slopes = res_slopes, res_cov = res_cov,
    res_cov_inv = res_cov_inv, sinv_method = sinv_method
  )

  # sandwich
  a_inv %*% m_b %*% a_inv
}
