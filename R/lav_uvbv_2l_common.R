# shared infrastructure for the univariate/bivariate TWO-LEVEL modules:
#
#   lav_uvreg_2l.R  univariate two-level Gaussian (random intercept)
#   lav_uvord_2l.R  univariate two-level ordinal probit (random intercept)
#   lav_bvreg_2l.R  bivariate two-level Gaussian (two-step)
#
# YR 2026 (two-level WLS, phase 2)
#
# These are the stage-wise ('univariate + bivariate') fitters for the
# unrestricted two-level model. The Gaussian kernel below works with the
# per-cluster sufficient statistics
#
#     nj      cluster size
#     ybar_j  cluster mean vector (d)
#     W_j     within-cluster scatter sum (y_ij - ybar_j)(y_ij - ybar_j)'
#
# and evaluates, for Sigma_j = I_nj (x) Sigma_w + J_nj (x) Sigma_b,
#
#  -2 logl_j = nj*d*log(2pi) + (nj-1) log|Sw| + log|Vj|
#              + tr(Sw^{-1} W_j) + nj r_j' Vj^{-1} r_j
#
# with Vj = Sw + nj Sb and r_j = ybar_j - mu. The derivative engine takes
# the pattern matrices (dmu, dSw, dSb) of any parameter and returns the
# per-cluster derivative of -2 logl_j:
#
#  d(-2 logl_j) = (nj-1) tr(Sw^{-1} dSw) - tr(Sw^{-1} dSw Sw^{-1} W_j)
#               + tr(Vj^{-1} dVj) - nj r_j' Vj^{-1} dVj Vj^{-1} r_j
#               - 2 nj r_j' Vj^{-1} dmu
#
# with dVj = dSw + nj dSb.

# per-cluster sufficient statistics
# y: numeric matrix (n x d); cluster_idx: integer vector (1..J)
lav_2l_cluster_stats <- function(y = NULL, cluster_idx = NULL) {
  y <- as.matrix(y)
  d <- NCOL(y)
  njs <- tabulate(cluster_idx)
  nclusters <- length(njs)

  ysum <- rowsum.default(y, group = cluster_idx, reorder = TRUE)
  ybar <- ysum / njs

  # within scatter, stored as vech per cluster (J x d*(d+1)/2)
  yc <- y - ybar[cluster_idx, , drop = FALSE]
  cp_idx <- which(lower.tri(diag(d), diag = TRUE), arr.ind = TRUE)
  wmat <- matrix(0, NROW(y), NROW(cp_idx))
  for (k in seq_len(NROW(cp_idx))) {
    wmat[, k] <- yc[, cp_idx[k, 1]] * yc[, cp_idx[k, 2]]
  }
  wvech <- rowsum.default(wmat, group = cluster_idx, reorder = TRUE)

  list(
    d = d, njs = njs, nclusters = nclusters, nobs = NROW(y),
    ybar = ybar, wvech = wvech,
    wtot = colSums(wmat) # total within scatter (vech)
  )
}

# make a d x d symmetric matrix from its vech
lav_2l_vech2sym <- function(x, d) {
  s <- matrix(0, d, d)
  s[lower.tri(s, diag = TRUE)] <- x
  s[upper.tri(s)] <- t(s)[upper.tri(s)]
  s
}

# per-cluster -2 loglik of the two-level Gaussian model
# cs: output of lav_2l_cluster_stats(); mu (d), sigma_w, sigma_b (d x d)
# returns a vector of length nclusters (or +Inf if not admissible)
lav_2l_gauss_m2ll <- function(cs = NULL, mu = NULL,
                              sigma_w = NULL, sigma_b = NULL) {
  d <- cs$d
  njs <- cs$njs
  nsizes <- sort(unique(njs))

  sw_inv <- try(solve(sigma_w), silent = TRUE)
  det_sw <- det(sigma_w)
  if (inherits(sw_inv, "try-error") || det_sw <= 0) {
    return(rep(+Inf, cs$nclusters))
  }

  # within part: (nj-1) log|Sw| + tr(Sw^{-1} W_j)
  # tr(Sw^{-1} W_j) from the vech: diagonal entries once, off-diag twice
  swi_vech <- sw_inv[lower.tri(sw_inv, diag = TRUE)]
  mult <- ifelse(
    which(lower.tri(diag(d), diag = TRUE)) %in%
      which(lower.tri(diag(d), diag = FALSE)), 2, 1
  )
  tr_w <- drop(cs$wvech %*% (swi_vech * mult))
  out <- njs * d * log(2 * pi) + (njs - 1) * log(det_sw) + tr_w

  # between part per unique cluster size
  r <- cs$ybar - matrix(mu, cs$nclusters, d, byrow = TRUE)
  for (nj in nsizes) {
    idx <- which(njs == nj)
    vj <- sigma_w + nj * sigma_b
    vj_inv <- try(solve(vj), silent = TRUE)
    det_vj <- det(vj)
    if (inherits(vj_inv, "try-error") || det_vj <= 0) {
      return(rep(+Inf, cs$nclusters))
    }
    rj <- r[idx, , drop = FALSE]
    quad <- rowSums((rj %*% vj_inv) * rj)
    out[idx] <- out[idx] + log(det_vj) + nj * quad
  }

  out
}

# per-cluster derivative of -2 logl_j for a list of parameters
# params: list of list(dmu = (d), dsw = (d x d), dsb = (d x d)); any of the
#         three may be NULL (treated as zero)
# returns a (nclusters x length(params)) matrix
lav_2l_gauss_dm2ll <- function(cs = NULL, mu = NULL,
                               sigma_w = NULL, sigma_b = NULL,
                               params = NULL) {
  d <- cs$d
  njs <- cs$njs
  nsizes <- sort(unique(njs))
  npar <- length(params)

  sw_inv <- solve(sigma_w)
  r <- cs$ybar - matrix(mu, cs$nclusters, d, byrow = TRUE)

  out <- matrix(0, cs$nclusters, npar)

  # precompute Sw^{-1} dSw and Sw^{-1} dSw Sw^{-1} per parameter
  for (p in seq_len(npar)) {
    dmu <- params[[p]]$dmu
    dsw <- params[[p]]$dsw
    dsb <- params[[p]]$dsb

    # within part
    if (!is.null(dsw)) {
      a <- sw_inv %*% dsw
      tr_a <- sum(diag(a))
      b <- a %*% sw_inv # Sw^{-1} dSw Sw^{-1}
      b_vech <- b[lower.tri(b, diag = TRUE)]
      mult <- ifelse(
        which(lower.tri(diag(d), diag = TRUE)) %in%
          which(lower.tri(diag(d), diag = FALSE)), 2, 1
      )
      # note: tr(Sw^{-1} dSw Sw^{-1} W_j); B is symmetric iff dSw is
      tr_bw <- drop(cs$wvech %*% (b_vech * mult))
      out[, p] <- out[, p] + (njs - 1) * tr_a - tr_bw
    }

    # between part per cluster size
    for (nj in nsizes) {
      idx <- which(njs == nj)
      vj <- sigma_w + nj * sigma_b
      vj_inv <- solve(vj)
      dvj <- matrix(0, d, d)
      if (!is.null(dsw)) dvj <- dvj + dsw
      if (!is.null(dsb)) dvj <- dvj + nj * dsb
      rj <- r[idx, , drop = FALSE]
      rv <- rj %*% vj_inv
      if (!is.null(dsw) || !is.null(dsb)) {
        out[idx, p] <- out[idx, p] + sum(diag(vj_inv %*% dvj)) -
          njs[idx] * rowSums((rv %*% dvj) * rv)
      }
      if (!is.null(dmu)) {
        out[idx, p] <- out[idx, p] - 2 * njs[idx] * drop(rv %*% dmu)
      }
    }
  }

  out
}

# standardized Gauss-Hermite nodes for integrating over N(0, 1):
#   integral f(u) phi(u) du ~= sum_q w_q f(x_q)
# (the caller rescales: b_q = sqrt(vb) * x_q)
lav_2l_gh_xw <- function(ngh = 21L) {
  gh <- lav_integration_gauss_hermite(n = ngh, dnorm = TRUE,
                                      mean = 0, sd = 1)
  list(x = gh$x, w = gh$w, logw = log(gh$w))
}
