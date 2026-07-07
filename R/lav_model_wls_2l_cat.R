# model-implied statistics + Jacobian for TWO-LEVEL WLS with categorical
# data (the model-side counterpart of lav_samp_wls_2l_cat)
#
# YR 2026 (two-level WLS, phase 5/6)
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
#   2. within means of the within-only (continuous) variables
#   3. within variances of the continuous variables
#   4. within covariances/'correlations' (vech, no diagonal), standardized
#   5. between means (continuous; standardized 0 for ordinal) -- between
#      variables only
#   6. vech(Sigma_b, diagonal = TRUE), standardized -- between variables
#      only
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
  if (lavmodel@conditional.x) {
    # lav_model_implied() does not forward delta = FALSE to the slopes
    imp$res.slopes <- lav_model_pi(lavmodel = lavmodel, glist = glist,
                                   delta = FALSE)
  }
  imp$th.raw <- lav_m07_th_raw(lavmodel = lavmodel, glist = glist)
  imp
}

# per-group block info: the within (b1) variable layout, the ordinal
# positions, and the between (b2) -> within (b1) variable mapping
lav_m07_block_info <- function(lavmodel = NULL, g = NULL,
                               nw = NULL, nb = NULL) {
  b1 <- (g - 1L) * 2L + 1L
  b2 <- (g - 1L) * 2L + 2L

  # ordinal variables at the within level (level-1 threshold rows exist,
  # fixed to zero, with the same variable-major structure as level 2)
  th_idx_1 <- lavmodel@th.idx[[b1]]
  ord_1 <- sort(unique(th_idx_1[th_idx_1 > 0L]))

  # between -> within variable mapping
  if (nb == nw) {
    w2b <- seq_len(nw)
  } else {
    # from the dimNames of the (ov x ov) theta matrices per block
    mm_idx <- lav_model_group_mm_indices(lavmodel@nmat)
    nm <- names(lavmodel@GLIST)
    t1 <- mm_idx[[b1]][which(nm[mm_idx[[b1]]] == "theta")]
    t2 <- mm_idx[[b2]][which(nm[mm_idx[[b2]]] == "theta")]
    ov1 <- lavmodel@dimNames[[t1]][[1]]
    ov2 <- lavmodel@dimNames[[t2]][[1]]
    if (is.null(ov1) || is.null(ov2)) {
      lav_msg_stop(gettext(
        "cannot map the between-level variables to the within level
        (no dimNames); this should not happen."))
    }
    w2b <- match(ov2, ov1)
    if (anyNA(w2b)) {
      lav_msg_stop(gettext(
        "between-only variables are not supported (yet) for two-level
        (D)WLS estimation with categorical data."))
    }
  }

  between_flag <- seq_len(nw) %in% w2b
  list(
    b1 = b1, b2 = b2, th_idx_1 = th_idx_1, ord_1 = ord_1,
    w2b = w2b, between_flag = between_flag
  )
}

# model-implied (standardized) statistic vector; one vector per GROUP
lav_m07_wls_est <- function(lavmodel = NULL, glist = NULL) {
  nlevels <- 2L
  ngroups <- lavmodel@nblocks %/% nlevels
  conditional_x <- lavmodel@conditional.x
  imp <- lav_m07_implied_raw(lavmodel, glist = glist)

  out <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    if (conditional_x) {
      sigma_w <- imp$res.cov[[(g - 1L) * 2L + 1L]]
      sigma_b <- imp$res.cov[[(g - 1L) * 2L + 2L]]
    } else {
      sigma_w <- imp$cov[[(g - 1L) * 2L + 1L]]
      sigma_b <- imp$cov[[(g - 1L) * 2L + 2L]]
    }
    nvar <- NROW(sigma_w)
    nb <- NROW(sigma_b)
    bi <- lav_m07_block_info(lavmodel, g = g, nw = nvar, nb = nb)

    if (conditional_x) {
      mean_w <- imp$res.int[[bi$b1]]
      mu_b <- imp$res.int[[bi$b2]]
      pi_w <- imp$res.slopes[[bi$b1]]
      pi_b <- imp$res.slopes[[bi$b2]]
    } else {
      mean_w <- imp$mean[[bi$b1]]
      mu_b <- imp$mean[[bi$b2]]
      pi_w <- matrix(0, nvar, 0L)
      pi_b <- matrix(0, nb, 0L)
    }
    th_raw2 <- imp$th.raw[[bi$b2]]
    th_idx_2 <- lavmodel@th.idx[[bi$b2]]

    ord <- (seq_len(nvar) %in% bi$ord_1)

    # within scale factors
    delta_1 <- ifelse(ord, 1 / sqrt(diag(sigma_w)), 1.0)

    # 1. thresholds (standardized); 0 entries for continuous variables
    # (the k-th ordinal threshold at level 1 corresponds to the k-th
    #  threshold parameter at level 2)
    th_std <- numeric(length(bi$th_idx_1))
    ord_pos_1 <- which(bi$th_idx_1 > 0L)
    ord_pos_2 <- which(th_idx_2 > 0L)
    th_std[ord_pos_1] <- th_raw2[ord_pos_2] *
      delta_1[bi$th_idx_1[ord_pos_1]]

    # 2. within means/intercepts of the within-only (continuous) variables
    mu_w <- mean_w[!bi$between_flag]

    # 3. standardized within slopes (vec, column-major over (var, exo))
    piw_std <- as.numeric(delta_1 * pi_w)

    # 4. within variances of the continuous variables
    vw <- diag(sigma_w)[!ord]

    # 5. standardized within covariances (vech, no diagonal)
    sw_std <- diag(delta_1) %*% sigma_w %*% diag(delta_1)
    cw <- lav_mat_vech(sw_std, diagonal = FALSE)

    # 6./7./8. standardized between means/intercepts + slopes +
    # covariance matrix (all scaled by the WITHIN scale factors)
    delta_b <- delta_1[bi$w2b]
    mub_std <- delta_b * mu_b
    pib_std <- as.numeric(delta_b * pi_b)
    sb_std <- diag(delta_b, nrow = nb) %*% sigma_b %*%
      diag(delta_b, nrow = nb)

    out[[g]] <- c(
      th_std, mu_w, piw_std, vw, cw, mub_std, pib_std,
      lav_mat_vech(sb_std, diagonal = TRUE)
    )
  }

  out
}

# the standardization map H (per group): rows = statistics, columns = raw
# rows. The raw row layout follows the continuous-view lav_model_delta():
#
#   conditional.x = FALSE:
#     [ th (nth2) ; mean_w (nvar) ; vech(Sigma_w) ;
#       mean_b (nb) ; vech(Sigma_b) ]
#   conditional.x = TRUE (per block: interleaved intercept + slope rows,
#   one set per variable):
#     [ th (nth2) ; {int_p, slopes_p (nq1)} x nvar ; vech(Sigma_w) ;
#       {int_k, slopes_k (nq2)} x nb ; vech(Sigma_b) ]
#
# (the unconditional case is the nq1 = nq2 = 0 special case)
lav_m07_hmat <- function(sigma_w = NULL, sigma_b = NULL, mu_b = NULL,
                         pi_w = NULL, pi_b = NULL,
                         th_raw2 = NULL, th_idx_1 = NULL,
                         th_idx_2 = NULL, w2b = NULL) {
  nvar <- NROW(sigma_w)
  nb <- NROW(sigma_b)
  if (is.null(pi_w)) {
    pi_w <- matrix(0, nvar, 0L)
  }
  if (is.null(pi_b)) {
    pi_b <- matrix(0, nb, 0L)
  }
  nq1 <- NCOL(pi_w)
  nq2 <- NCOL(pi_b)
  between_flag <- seq_len(nvar) %in% w2b
  ord <- (seq_len(nvar) %in% th_idx_1[th_idx_1 > 0L])
  num_idx <- which(!ord)
  nth1 <- length(th_idx_1) # output th rows (within-block layout)
  nth2 <- length(th_idx_2) # raw th rows (between-block tau vector)
  pstar_d <- nvar * (nvar + 1L) / 2L # vech with diagonal (within)
  pstar <- nvar * (nvar - 1L) / 2L # vech without diagonal
  pstar_b <- nb * (nb + 1L) / 2L # vech with diagonal (between)

  delta_1 <- ifelse(ord, 1 / sqrt(diag(sigma_w)), 1.0)
  # d delta_p / d sigma_w,pp = -0.5 * delta_p^3 (ordinal only)
  ddelta <- ifelse(ord, -0.5 * delta_1^3, 0.0)

  # raw row indices
  r_th <- seq_len(nth2)
  r_int1 <- nth2 + (seq_len(nvar) - 1L) * (1L + nq1) + 1L
  # slope rows of variable p: r_int1[p] + 1..nq1
  r_sw <- nth2 + nvar * (1L + nq1) + seq_len(pstar_d)
  r_int2 <- nth2 + nvar * (1L + nq1) + pstar_d +
    (seq_len(nb) - 1L) * (1L + nq2) + 1L
  r_sb <- nth2 + nvar * (1L + nq1) + pstar_d + nb * (1L + nq2) +
    seq_len(pstar_b)
  nraw <- nth2 + nvar * (1L + nq1) + pstar_d + nb * (1L + nq2) + pstar_b

  # position of sigma[i,j] (i >= j) within vech(, diagonal = TRUE)
  vech_pos <- matrix(0L, nvar, nvar)
  vech_pos[lower.tri(vech_pos, diag = TRUE)] <- seq_len(pstar_d)
  vech_pos[upper.tri(vech_pos)] <- t(vech_pos)[upper.tri(vech_pos)]
  diag_pos <- diag(vech_pos)
  # ... and for the between block
  vech_pos_b <- matrix(0L, nb, nb)
  vech_pos_b[lower.tri(vech_pos_b, diag = TRUE)] <- seq_len(pstar_b)
  vech_pos_b[upper.tri(vech_pos_b)] <- t(vech_pos_b)[upper.tri(vech_pos_b)]

  nstat <- nth1 + sum(!between_flag) + nvar * nq1 + length(num_idx) +
    pstar + nb + nb * nq2 + pstar_b
  h <- matrix(0, nstat, nraw)
  row <- 0L

  # 1. thresholds: th*_k = delta_p th_k
  ord_pos_2 <- which(th_idx_2 > 0L)
  kk <- 0L
  for (k in seq_len(nth1)) {
    row <- row + 1L
    p <- th_idx_1[k]
    if (p > 0L) {
      kk <- kk + 1L
      h[row, r_th[ord_pos_2[kk]]] <- delta_1[p]
      h[row, r_sw[diag_pos[p]]] <- th_raw2[ord_pos_2[kk]] * ddelta[p]
    }
    # continuous entries: fixed at zero (zero row)
  }

  # 2. within means/intercepts of the within-only (continuous) variables
  for (p in which(!between_flag)) {
    row <- row + 1L
    h[row, r_int1[p]] <- 1
  }

  # 3. standardized within slopes: pi*_pq = delta_p pi_pq
  # (vec, column-major: all variables for exo 1, then exo 2, ...)
  for (q in seq_len(nq1)) {
    for (p in seq_len(nvar)) {
      row <- row + 1L
      h[row, r_int1[p] + q] <- delta_1[p]
      h[row, r_sw[diag_pos[p]]] <- pi_w[p, q] * ddelta[p]
    }
  }

  # 4. within variances of continuous variables
  for (p in num_idx) {
    row <- row + 1L
    h[row, r_sw[diag_pos[p]]] <- 1
  }

  # 5. standardized within covariances (i > j)
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

  # 6. standardized between means/intercepts (between variables only)
  for (k in seq_len(nb)) {
    row <- row + 1L
    p <- w2b[k]
    h[row, r_int2[k]] <- delta_1[p]
    h[row, r_sw[diag_pos[p]]] <- mu_b[k] * ddelta[p]
  }

  # 7. standardized between slopes: pi*_kq = delta_{w2b(k)} pi_kq
  for (q in seq_len(nq2)) {
    for (k in seq_len(nb)) {
      row <- row + 1L
      p <- w2b[k]
      h[row, r_int2[k] + q] <- delta_1[p]
      h[row, r_sw[diag_pos[p]]] <- pi_b[k, q] * ddelta[p]
    }
  }

  # 8. standardized between covariances (vech with diagonal, i >= j)
  for (jj in seq_len(nb)) {
    for (ii in jj:nb) {
      row <- row + 1L
      pi_1 <- w2b[ii]; pj_1 <- w2b[jj]
      s_ij <- sigma_b[ii, jj]
      h[row, r_sb[vech_pos_b[ii, jj]]] <- delta_1[pi_1] * delta_1[pj_1]
      h[row, r_sw[diag_pos[pi_1]]] <- h[row, r_sw[diag_pos[pi_1]]] +
        s_ij * delta_1[pj_1] * ddelta[pi_1]
      h[row, r_sw[diag_pos[pj_1]]] <- h[row, r_sw[diag_pos[pj_1]]] +
        s_ij * delta_1[pi_1] * ddelta[pj_1]
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
  # parameters, so one-sided finite differences are exact.
  # IMPORTANT: the base thresholds must be evaluated at the SAME
  # parameter vector as the perturbed ones (x0 from the lavmodel object,
  # which may differ from the current glist state!)
  x0 <- lav_model_get_parameters(lavmodel)
  lavmodel_0 <- lav_model_set_parameters(lavmodel, x = x0)
  th_base <- lav_m07_th_raw(lavmodel = lavmodel_0,
                            glist = lavmodel_0@GLIST)
  nx <- length(x0)
  th_jac <- vector("list", lavmodel@nblocks)
  nth_b <- sapply(th_base, length)
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
        th_jac[[b]][, k] <- th1[[b]] - th_base[[b]]
      }
    }
  }

  conditional_x <- lavmodel@conditional.x
  delta_out <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    if (conditional_x) {
      sigma_w <- imp$res.cov[[(g - 1L) * 2L + 1L]]
      sigma_b <- imp$res.cov[[(g - 1L) * 2L + 2L]]
    } else {
      sigma_w <- imp$cov[[(g - 1L) * 2L + 1L]]
      sigma_b <- imp$cov[[(g - 1L) * 2L + 2L]]
    }
    nvar <- NROW(sigma_w)
    nb <- NROW(sigma_b)
    bi <- lav_m07_block_info(lavmodel, g = g, nw = nvar, nb = nb)
    if (conditional_x) {
      mu_b <- imp$res.int[[bi$b2]]
      pi_w <- imp$res.slopes[[bi$b1]]
      pi_b <- imp$res.slopes[[bi$b2]]
    } else {
      mu_b <- imp$mean[[bi$b2]]
      pi_w <- pi_b <- NULL
    }
    th_raw2 <- imp$th.raw[[bi$b2]]
    th_idx_2 <- lavmodel@th.idx[[bi$b2]]

    # stack the raw rows: [th (b2); int/slope rows W; vech Sw;
    #                      int/slope rows B; vech Sb]
    draw <- rbind(th_jac[[bi$b2]], delta_raw[[g]])

    h <- lav_m07_hmat(
      sigma_w = sigma_w, sigma_b = sigma_b, mu_b = mu_b,
      pi_w = pi_w, pi_b = pi_b,
      th_raw2 = th_raw2, th_idx_1 = bi$th_idx_1,
      th_idx_2 = th_idx_2, w2b = bi$w2b
    )
    delta_out[[g]] <- h %*% draw
  }

  delta_out
}