# model-implied statistics + Jacobian for TWO-LEVEL WLS with categorical
# data (the model-side counterpart of lav_samp_wls_2l_cat)
#
# YR 2026 (two-level WLS, phase 5)
#
# The unrestricted (stage-wise) statistics live on the scale where the
# WITHIN-level latent responses have variance 1. The structured model
# (theta parameterization) fixes the within RESIDUAL variances to 1, so
# its implied within variances sigma**_w,pp = 1 + (common part) exceed 1:
# a different scale convention. The model-implied statistics are
# therefore STANDARDIZED by the within-level scale factors
#
#     delta_p = 1 / sqrt(sigma**_w,pp)      (ordinal variables)
#     delta_p = 1                           (continuous variables)
#
# applied to the thresholds, the within covariance matrix, the between
# means AND the between covariance matrix (note: the between matrix is
# scaled by the WITHIN scale factors -- a cross-block standardization).
#
# statistic vector layout (per group; matches lav_samp_wls_2l_cat):
#   1. th: thresholds (ordinal, standardized); 0 for continuous variables
#   2. within variances of the continuous variables
#   3. within covariances/'correlations' (vech, no diagonal), standardized
#   4. between means (continuous; standardized 0 for ordinal)
#   5. vech(Sigma_b, diagonal = TRUE), standardized
#
# The Jacobian is assembled as H %*% Delta_raw, where Delta_raw contains
# the raw (unstandardized) derivative rows
#   [ th ; mean_w ; vech(Sigma_w) ; mean_b ; vech(Sigma_b) ]
# and H is the (implied-value dependent) sparse standardization map.
# Delta_raw is obtained from the standard continuous two-level Delta
# machinery (a 'continuous view' of the model) plus the threshold rows
# (exact finite differences: th is affine in the free parameters).

# a 'continuous view' of the model: same GLIST/parameters, but the
# categorical machinery switched off, so that lav_model_delta() and
# lav_model_implied() return the raw means/covariances per block
lav_m07_model_nocat <- function(lavmodel = NULL) {
  lavmodel2 <- lavmodel
  lavmodel2@categorical <- FALSE
  lavmodel2@correlation <- FALSE
  lavmodel2@th.idx <- lapply(lavmodel@th.idx, function(x) integer(0L))
  # no implied standardization:
  # - under parameterization = "theta", lav_model_set_parameters() stores
  #   the per-block scaling values 1/sqrt(diag(Sigma*)) in the GLIST
  #   'delta' matrices; switch to "delta" and reset those matrices to 1,
  #   so that all scaling is the identity and the implied moments (and
  #   their Jacobian) come out unstandardized (raw)
  lavmodel2@parameterization <- "delta"
  delta_mm <- which(names(lavmodel2@GLIST) == "delta")
  for (m in delta_mm) {
    lavmodel2@GLIST[[m]][] <- 1
  }
  lavmodel2
}

# raw thresholds per block: the tau parameters themselves, without any
# scaling (lav_model_th applies the delta/theta standardization);
# entries for numeric variables (th.idx == 0) are set to zero
lav_m07_th_raw <- function(lavmodel = NULL, glist = NULL) {
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }
  mm_idx <- lav_model_group_mm_indices(lavmodel@nmat)
  out <- vector("list", lavmodel@nblocks)
  for (b in seq_len(lavmodel@nblocks)) {
    th_idx_b <- lavmodel@th.idx[[b]]
    th_b <- numeric(length(th_idx_b))
    mlist <- glist[mm_idx[[b]]]
    if (!is.null(mlist$tau) && any(th_idx_b > 0L)) {
      th_b[th_idx_b > 0L] <- as.numeric(mlist$tau)
    }
    out[[b]] <- th_b
  }
  out
}

# raw model-implied moments per block (no delta/correlation scaling:
# delta = FALSE switches off the implied standardization)
lav_m07_implied_raw <- function(lavmodel = NULL, glist = NULL) {
  imp <- lav_model_implied(lavmodel, glist = glist, delta = FALSE)
  imp$th.raw <- lav_m07_th_raw(lavmodel = lavmodel, glist = glist)
  imp
}

# model-implied (standardized) statistic vector; one vector per GROUP
lav_m07_wls_est <- function(lavmodel = NULL, glist = NULL) {
  nlevels <- 2L
  ngroups <- lavmodel@nblocks %/% nlevels
  imp <- lav_m07_implied_raw(lavmodel, glist = glist)

  out <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    b1 <- (g - 1L) * nlevels + 1L
    b2 <- (g - 1L) * nlevels + 2L
    sigma_w <- imp$cov[[b1]]
    sigma_b <- imp$cov[[b2]]
    mu_b <- imp$mean[[b2]]
    # note: the thresholds are between-level parameters (the level-1
    # thresholds are fixed to zero); numeric entries of th are ignored
    # (the between means enter separately, unnegated)
    th_raw <- imp$th.raw[[b2]]
    th_idx <- lavmodel@th.idx[[b2]]

    nvar <- NROW(sigma_w)
    ord <- (seq_len(nvar) %in% th_idx)
    num_idx <- which(!ord)

    # within scale factors
    delta_1 <- ifelse(ord, 1 / sqrt(diag(sigma_w)), 1.0)

    # 1. thresholds (standardized); 0 entries for continuous variables
    th_std <- numeric(length(th_idx))
    ord_th <- which(th_idx > 0L)
    th_std[ord_th] <- th_raw[ord_th] * delta_1[th_idx[ord_th]]

    # 2. within variances of the continuous variables
    vw <- diag(sigma_w)[num_idx]

    # 3. standardized within covariances (vech, no diagonal)
    sw_std <- diag(delta_1) %*% sigma_w %*% diag(delta_1)
    cw <- lav_mat_vech(sw_std, diagonal = FALSE)

    # 4./5. standardized between means + covariance matrix
    mub_std <- delta_1 * mu_b
    sb_std <- diag(delta_1) %*% sigma_b %*% diag(delta_1)

    out[[g]] <- c(
      th_std, vw, cw, mub_std,
      lav_mat_vech(sb_std, diagonal = TRUE)
    )
  }

  out
}

# the standardization map H (per group): rows = statistics, columns = raw
# rows [ th (nth) ; mean_w (nvar) ; vech(Sigma_w) ; mean_b ; vech(Sigma_b) ]
lav_m07_hmat <- function(sigma_w = NULL, sigma_b = NULL, mu_b = NULL,
                         th_raw = NULL, th_idx = NULL) {
  nvar <- NROW(sigma_w)
  ord <- (seq_len(nvar) %in% th_idx)
  num_idx <- which(!ord)
  nth <- length(th_idx)
  pstar_d <- nvar * (nvar + 1L) / 2L # vech with diagonal
  pstar <- nvar * (nvar - 1L) / 2L # vech without diagonal

  delta_1 <- ifelse(ord, 1 / sqrt(diag(sigma_w)), 1.0)
  # d delta_p / d sigma_w,pp = -0.5 * delta_p^3 (ordinal only)
  ddelta <- ifelse(ord, -0.5 * delta_1^3, 0.0)

  # raw row indices
  r_th <- seq_len(nth)
  r_muw <- nth + seq_len(nvar)
  r_sw <- nth + nvar + seq_len(pstar_d)
  r_mub <- nth + nvar + pstar_d + seq_len(nvar)
  r_sb <- nth + nvar + pstar_d + nvar + seq_len(pstar_d)
  nraw <- nth + 2L * (nvar + pstar_d)

  # position of sigma[i,j] (i >= j) within vech(, diagonal = TRUE)
  vech_pos <- matrix(0L, nvar, nvar)
  vech_pos[lower.tri(vech_pos, diag = TRUE)] <- seq_len(pstar_d)
  vech_pos[upper.tri(vech_pos)] <- t(vech_pos)[upper.tri(vech_pos)]
  diag_pos <- diag(vech_pos)

  nstat <- nth + length(num_idx) + pstar + nvar + pstar_d
  h <- matrix(0, nstat, nraw)
  row <- 0L

  # 1. thresholds: th*_k = delta_p th_k
  for (k in seq_len(nth)) {
    row <- row + 1L
    p <- th_idx[k]
    if (p > 0L) {
      h[row, r_th[k]] <- delta_1[p]
      h[row, r_sw[diag_pos[p]]] <- th_raw[k] * ddelta[p]
    }
    # continuous entries: fixed at zero (zero row)
  }

  # 2. within variances of continuous variables
  for (p in num_idx) {
    row <- row + 1L
    h[row, r_sw[diag_pos[p]]] <- 1
  }

  # 3. standardized within covariances (i > j)
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      row <- row + 1L
      s_ij <- sigma_w[i, j]
      h[row, r_sw[vech_pos[i, j]]] <- delta_1[i] * delta_1[j]
      h[row, r_sw[diag_pos[i]]] <- h[row, r_sw[diag_pos[i]]] +
        s_ij * delta_1[j] * ddelta[i]
      h[row, r_sw[diag_pos[j]]] <- h[row, r_sw[diag_pos[j]]] +
        s_ij * delta_1[i] * ddelta[j]
    }
  }

  # 4. standardized between means
  for (p in seq_len(nvar)) {
    row <- row + 1L
    h[row, r_mub[p]] <- delta_1[p]
    h[row, r_sw[diag_pos[p]]] <- mu_b[p] * ddelta[p]
  }

  # 5. standardized between covariances (vech with diagonal, i >= j)
  for (j in seq_len(nvar)) {
    for (i in j:nvar) {
      row <- row + 1L
      s_ij <- sigma_b[i, j]
      h[row, r_sb[vech_pos[i, j]]] <- delta_1[i] * delta_1[j]
      h[row, r_sw[diag_pos[i]]] <- h[row, r_sw[diag_pos[i]]] +
        s_ij * delta_1[j] * ddelta[i]
      h[row, r_sw[diag_pos[j]]] <- h[row, r_sw[diag_pos[j]]] +
        s_ij * delta_1[i] * ddelta[j]
    }
  }

  h
}

# Jacobian of the standardized statistic vector; one matrix per GROUP
# (rows match lav_m07_wls_est, columns match lav_model_delta)
lav_m07_delta <- function(lavmodel = NULL, glist = NULL) {
  nlevels <- 2L
  ngroups <- lavmodel@nblocks %/% nlevels

  # raw continuous-style Delta: per group rbind of
  # [mean_w; vech(Sigma_w); mean_b; vech(Sigma_b)]
  lavmodel2 <- lav_m07_model_nocat(lavmodel)
  glist2 <- glist
  if (!is.null(glist2)) {
    # reset the delta matrices here as well (see lav_m07_model_nocat)
    delta_mm <- which(names(glist2) == "delta")
    for (m in delta_mm) {
      glist2[[m]][] <- 1
    }
  }
  delta_raw <- lav_model_delta(lavmodel = lavmodel2, glist = glist2)

  # raw implied (for the H map)
  imp <- lav_m07_implied_raw(lavmodel, glist = glist)

  # threshold Jacobian: the raw thresholds are affine in the free
  # parameters, so one-sided finite differences are exact
  x0 <- lav_model_get_parameters(lavmodel)
  th0 <- lav_m07_th_raw(lavmodel = lavmodel, glist = glist)
  nx <- length(x0)
  th_jac <- vector("list", lavmodel@nblocks)
  nth_b <- sapply(th0, length)
  for (b in seq_len(lavmodel@nblocks)) {
    th_jac[[b]] <- matrix(0, nth_b[b], nx)
  }
  for (k in seq_len(nx)) {
    x1 <- x0
    x1[k] <- x1[k] + 1
    lavmodel_k <- lav_model_set_parameters(lavmodel, x = x1)
    th1 <- lav_m07_th_raw(lavmodel = lavmodel_k,
                          glist = lavmodel_k@GLIST)
    for (b in seq_len(lavmodel@nblocks)) {
      if (nth_b[b] > 0L) {
        th_jac[[b]][, k] <- th1[[b]] - th0[[b]]
      }
    }
  }

  delta_out <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    b1 <- (g - 1L) * nlevels + 1L
    b2 <- (g - 1L) * nlevels + 2L
    sigma_w <- imp$cov[[b1]]
    sigma_b <- imp$cov[[b2]]
    mu_b <- imp$mean[[b2]]
    th_raw <- imp$th.raw[[b2]]
    th_idx <- lavmodel@th.idx[[b2]]

    # stack the raw rows: [th; mean_w; vech Sw; mean_b; vech Sb]
    draw <- rbind(th_jac[[b2]], delta_raw[[g]])

    h <- lav_m07_hmat(
      sigma_w = sigma_w, sigma_b = sigma_b, mu_b = mu_b,
      th_raw = th_raw, th_idx = th_idx
    )
    delta_out[[g]] <- h %*% draw
  }

  delta_out
}
