# shared computational kernels for the lav_mvnorm_* / lav_mvreg_* family
#
# these helpers centralize small blocks of code that used to be copy-pasted
# across lav_mvnorm.R, lav_mvnorm_h1.R, lav_mvnorm_missing.R,
# lav_mvnorm_missing_h1.R, lav_mvreg.R (and, for some, the cluster files)
#
# YR/CC 03 July 2026: first version (refactoring; no change in behavior)

# 1) ensure we have Sigma^{-1} (optionally with a "logdet" attribute)
#
#    - if sigma_inv is NULL, invert sigma_1 using lav_mat_sym_inverse()
#    - if sigma_inv is given but the "logdet" attribute is missing (and
#      logdet = TRUE), compute -log|Sigma^{-1}| from the eigenvalues
lav_mvn_sigma_inv <- function(sigma_1 = NULL,
                              sigma_inv = NULL,
                              sinv_method = "eigen",
                              logdet = FALSE) {
  if (is.null(sigma_inv)) {
    sigma_inv <- lav_mat_sym_inverse(
      s = sigma_1, logdet = logdet,
      sinv_method = sinv_method
    )
  } else if (logdet && is.null(attr(sigma_inv, "logdet"))) {
    # compute - ln|Sigma.inv|
    ev <- eigen(sigma_inv, symmetric = TRUE, only.values = TRUE)
    attr(sigma_inv, "logdet") <- -1 * sum(log(ev$values))
  }

  sigma_inv
}

# 2) total number of observations: sum(wt), or the number of rows of y
lav_mvn_nobs <- function(y = NULL, wt = NULL) {
  if (!is.null(wt)) {
    sum(wt)
  } else {
    NROW(y)
  }
}

# 2b) ML sample statistics (mean + cov) from raw data, weighted or not
lav_mvn_samp_stats <- function(y = NULL, wt = NULL) {
  if (!is.null(wt)) {
    out <- stats::cov.wt(y, wt = wt, method = "ML")
    list(mean = out$center, cov = out$cov)
  } else {
    list(
      mean = base::.colMeans(y, m = NROW(y), n = NCOL(y)),
      cov = lav_mat_cov(y)
    )
  }
}

# 3) W.tilde = S + (ybar - mu)(ybar - mu)' from raw data
lav_mvn_w_tilde <- function(y = NULL, wt = NULL, mu = NULL) {
  if (!is.null(wt)) {
    out <- stats::cov.wt(y, wt = wt, method = "ML")
    out$cov + tcrossprod(out$center - mu)
  } else {
    lav_mat_cov(y, mu = mu)
  }
}

# 3b) W.tilde from sample statistics
lav_mvn_w_tilde_samp <- function(sample_mean = NULL,
                                 sample_cov = NULL,
                                 mu = NULL) {
  sample_cov + tcrossprod(as.numeric(sample_mean) - as.numeric(mu))
}

# 4) casewise score kernel
#
#    input:  yc = centered data (y - mu), sigma_inv
#    output: list with
#            $mu    = casewise scores with respect to mu          (= yc Sinv)
#            $sigma = casewise scores with respect to vech(Sigma)
lav_mvn_sc_kernel <- function(yc = NULL, sigma_inv = NULL) {
  p <- NCOL(sigma_inv)

  # vech(Sigma.inv)
  isigma <- lav_mat_vech(sigma_inv)

  # postmultiply with Sigma.inv
  ysi <- yc %*% sigma_inv

  # tcrossprod
  idx1 <- lav_mat_vech_col_idx(p)
  idx2 <- lav_mat_vech_row_idx(p)
  z <- ysi[, idx1, drop = FALSE] * ysi[, idx2, drop = FALSE]

  # subtract isigma from each row
  sc <- t(t(z) - isigma)

  # adjust for lav_mat_dup_pre (not vech!)
  sc[, lav_mat_diagh_idx(p)] <- sc[, lav_mat_diagh_idx(p)] / 2

  list(mu = ysi, sigma = sc)
}

# 5) observed-information kernel (the i11/i21/i22 algebra)
#
#    given Sigma^{-1}, Sigma, W.tilde and (ybar - mu), return the blocks of
#    the unit observed information for mu + vech(Sigma):
#      i11 = Sigma^{-1}
#      i21 = D' ((Sigma^{-1} (ybar - mu)) %x% Sigma^{-1})
#      i22 = 1/2 D' (Sigma^{-1} %x% AAA) D,
#            AAA = Sigma^{-1} (2 W.tilde - Sigma) Sigma^{-1}
#    if mean_diff is NULL, only i22 is computed
lav_mvn_info_obs_kernel <- function(sigma_inv = NULL,
                                    sigma_1 = NULL,
                                    w_tilde = NULL,
                                    mean_diff = NULL) {
  aaa <- sigma_inv %*% (2 * w_tilde - sigma_1) %*% sigma_inv
  i22 <- lav_mvn_kron_dup_half(sigma_inv, aaa)

  if (is.null(mean_diff)) {
    list(i22 = i22)
  } else {
    i21 <- lav_mat_dup_pre((sigma_inv %*% mean_diff) %x% sigma_inv)
    list(i11 = sigma_inv, i21 = i21, i22 = i22)
  }
}

# 6) 1/2 * D' (x %x% y) D  (y defaults to x)
#
#    this is the ONE home for a (future) lavaanC fast path:
#    lavaanC::m_kronecker_dup_pre_post(x, y, 0.5) or
#    lavaanC::m_kronecker_dup_cor_pre_post(x, multiplicator = 0.5)
lav_mvn_kron_dup_half <- function(x, y = NULL, correlation = FALSE) {
  if (is.null(y)) {
    y <- x
  }
  if (correlation) {
    0.5 * lav_mat_dup_cor_pre_post(x %x% y)
  } else {
    0.5 * lav_mat_dup_pre_post(x %x% y)
  }
}

# 6b) 2 * D^+ (x %x% y) D^{+T}  (y defaults to x)
#
#     (single home for a future lavaanC fast path:
#      lavaanC::m_kronecker_dup_ginv_pre_post(x, y, multiplicator = 2.0))
lav_mvn_kron_dup_ginv2 <- function(x, y = NULL) {
  if (is.null(y)) {
    y <- x
  }
  2 * lav_mat_dup_ginv_pre_post(x %x% y)
}

# 7) zero out the rows/cols of an information matrix (for mu + vech(Sigma),
#    or vech(Sigma) only) that correspond to fixed exogenous variables
lav_mvn_zero_x_idx <- function(out,
                               n = NULL,
                               x_idx = integer(0L),
                               meanstructure = TRUE) {
  if (length(x_idx) == 0L) {
    return(out)
  }
  not_x <- lav_mat_vech_which_idx(
    n = n, idx = x_idx,
    add_idx_at_start = meanstructure
  )
  out[not_x, ] <- 0
  out[, not_x] <- 0

  out
}
