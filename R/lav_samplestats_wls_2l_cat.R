# two-level WLS/DWLS/ULS with categorical data: sample 'statistics' and
# their asymptotic covariance matrix, via the stage-wise (univariate +
# bivariate) estimation of the unrestricted two-level model
# (lav_muthen2007)
#
# YR 2026 (two-level WLS, phase 5)
#
# The per-group statistic vector follows the lavaan two-level WLS layout
# (the within and between blocks stacked per group, matching the rows of
# the per-group Delta matrix):
#
#   within block (categorical layout, theta parameterization):
#     1. th: thresholds (ordinal) and -mu entries (continuous, all zero
#        here: the means of shared variables live at the between level)
#     2. variances of the continuous variables (vw)
#     3. within 'correlations' (vech, no diagonal): cw in the theta
#        metric (ordinal latent responses have within variance 1)
#   between block (continuous layout):
#     4. between means: mu (continuous), 0 (ordinal, fixed)
#     5. vech(Sigma_b, diagonal = TRUE): vb / cb
#
# note: the exact layout is aligned with lav_model_wls_est() /
# lav_model_delta() for the two-level categorical model.

lav_samp_wls_2l_cat <- function(lavsamplestats = NULL,
                                lavdata = NULL,
                                lavoptions = NULL) {
  stopifnot(lavdata@nlevels == 2L)
  estimator <- lavoptions$estimator
  ngroups <- lavsamplestats@ngroups
  nlevels <- 2L

  wls_obs <- vector("list", length = ngroups)
  wls_v <- vector("list", length = ngroups)
  wls_vd <- vector("list", length = ngroups)
  nacov <- vector("list", length = ngroups)

  # 'h1' (unrestricted model) implied statistics, in the standard
  # multilevel layout (within/between block per group); consumed by
  # fit measures (srmr), the baseline model (starting values), ...
  h1_implied <- list(
    cov = vector("list", ngroups * nlevels),
    mean = vector("list", ngroups * nlevels)
  )

  for (g in seq_len(ngroups)) {
    lp <- lavdata@Lp[[g]]
    nobs <- lavsamplestats@nobs[[g]]

    # phase 5 (v1) restrictions
    if (length(lavsamplestats@x.idx[[g]]) > 0L) {
      lav_msg_stop(gettext(
        "fixed.x = TRUE is not supported (yet) for two-level (D)WLS
        estimation; use fixed.x = FALSE."))
    }
    ov_names_1 <- lavdata@ov.names.l[[g]][[1]]
    ov_names_2 <- lavdata@ov.names.l[[g]][[2]]
    if (length(setdiff(ov_names_2, ov_names_1)) > 0L) {
      lav_msg_stop(gettext(
        "two-level (D)WLS estimation with categorical data does not
        support between-only variables (yet)."))
    }
    if (!identical(ov_names_2,
                   ov_names_1[ov_names_1 %in% ov_names_2])) {
      lav_msg_stop(gettext(
        "two-level (D)WLS estimation with categorical data requires the
        between-level variables to appear in the same relative order as
        at the within level."))
    }

    # variable info (in the within-block variable order)
    ov_types <- lavdata@ov$type[match(ov_names_1, lavdata@ov$name)]
    between_flag <- ov_names_1 %in% ov_names_2
    if (any(ov_types == "ordered" & !between_flag)) {
      lav_msg_stop(gettext(
        "two-level (D)WLS estimation does not support within-only ORDINAL
        variables (yet)."))
    }

    # stage-wise estimation of the unrestricted two-level model
    m07 <- lav_muthen2007(
      y = lavdata@X[[g]],
      ov_names = ov_names_1,
      ov_types = ov_types,
      cluster_idx = lp$cluster.idx[[2]],
      between_flag = between_flag,
      ngh = lavoptions$integration.ngh,
      ngh2 = 13L,
      wls_w = TRUE
    )

    # map the stage-wise statistic vector + ACOV to the lavaan layout
    out <- lav_m07_repackage(m07 = m07, ov_types = ov_types,
                             between_flag = between_flag)

    wls_obs[[g]] <- out$wls_obs
    nacov[[g]] <- nobs * out$acov

    if (estimator == "WLS") {
      keep <- out$free_idx
      if (lp$nclusters[[2]] < length(keep)) {
        lav_msg_stop(gettextf(
          "estimator WLS for two-level data needs more clusters (%1$s) than
          sample statistics (%2$s) in group %3$s; use estimator WLSMV
          (or ULSMV) instead.", lp$nclusters[[2]], length(keep), g))
      }
      w_full <- matrix(0, nrow(nacov[[g]]), ncol(nacov[[g]]))
      w_keep <- try(
        lav_mat_sym_inverse(nacov[[g]][keep, keep, drop = FALSE]),
        silent = TRUE
      )
      if (inherits(w_keep, "try-error")) {
        lav_msg_stop(gettextf(
          "could not invert Gamma (acov of the two-level sample statistics)
          in group %s; use estimator WLSMV (or ULSMV) instead.", g))
      }
      w_full[keep, keep] <- w_keep
      wls_v[[g]] <- w_full
    } else if (estimator == "DWLS") {
      dacov <- diag(nacov[[g]])
      idacov <- ifelse(dacov > 0, 1 / dacov, 0)
      wls_v[[g]] <- diag(idacov, nrow = length(idacov))
      wls_vd[[g]] <- idacov
    } else if (estimator == "ULS") {
      wls_v[[g]] <- diag(length(wls_obs[[g]]))
      wls_vd[[g]] <- rep(1, length(wls_obs[[g]]))
    }

    # store the stage-wise 'sample' statistics for inspection/start values
    lavsamplestats@th[[g]] <- out$th
    lavsamplestats@th.idx[[g]] <- out$th_idx

    # unrestricted implied statistics (theta metric: the diagonal of
    # SIGMA.W is 1 for the ordinal variables); the between block spans
    # the between-level variables only, and the means of within-only
    # variables live at the within level
    b_idx <- which(between_flag)
    mu_w1 <- numeric(NROW(m07$SIGMA.W))
    mu_w1[!between_flag] <- m07$MU[!between_flag]
    h1_implied$cov[[(g - 1) * nlevels + 1L]] <- m07$SIGMA.W
    h1_implied$cov[[(g - 1) * nlevels + 2L]] <-
      m07$SIGMA.B[b_idx, b_idx, drop = FALSE]
    h1_implied$mean[[(g - 1) * nlevels + 1L]] <- mu_w1
    h1_implied$mean[[(g - 1) * nlevels + 2L]] <- m07$MU[b_idx]
  }

  lavsamplestats@WLS.obs <- wls_obs
  lavsamplestats@WLS.V <- wls_v
  lavsamplestats@WLS.VD <- wls_vd
  if (!lavsamplestats@NACOV.user) {
    lavsamplestats@NACOV <- nacov
  }

  list(
    lavsamplestats = lavsamplestats,
    lavh1 = list(
      implied = h1_implied,
      logl = list(loglik = as.numeric(NA), loglik.group = as.numeric(NA))
    )
  )
}

# repackage the stage-wise parameter vector (per variable, then per pair)
# and its ACOV into the lavaan two-level WLS layout:
#
#   1. th (per within variable: thresholds for ordinal, fixed 0 for
#      continuous)
#   2. within means of the within-only (continuous) variables
#   3. within variances of the continuous variables
#   4. within covariances (vech, no diagonal)
#   5. between means (continuous: mu; ordinal: fixed 0) -- between
#      variables only
#   6. vech(Sigma_b, diagonal = TRUE) -- between variables only
lav_m07_repackage <- function(m07 = NULL, ov_types = NULL,
                              between_flag = NULL) {
  nvar <- length(ov_types)
  ord <- (ov_types == "ordered")
  if (is.null(between_flag)) {
    between_flag <- rep(TRUE, nvar)
  }
  npar_uni <- m07$npar.uni
  uni_end <- cumsum(npar_uni)
  uni_start <- uni_end - npar_uni + 1L
  ntheta <- sum(npar_uni)
  pstar <- nvar * (nvar - 1L) / 2L
  pair_neta <- m07$pair.neta
  if (is.null(pair_neta)) {
    pair_neta <- rep(2L, pstar)
  }
  eta_end <- cumsum(pair_neta)
  eta_start <- eta_end - pair_neta + 1L

  # source positions (in the stage-wise vector) of each target entry;
  # 0 = structurally-fixed entry (value stored in fixed_val)
  src <- integer(0L)
  fixed_val <- numeric(0L)
  th_idx <- integer(0L)

  # --- within block ---
  # 1. th: thresholds (ordinal); fixed 0 slots (continuous)
  for (p in seq_len(nvar)) {
    if (ord[p]) {
      nth_p <- npar_uni[p] - 1L
      src <- c(src, uni_start[p] - 1L + seq_len(nth_p))
      fixed_val <- c(fixed_val, rep(NA_real_, nth_p))
      th_idx <- c(th_idx, rep(p, nth_p))
    } else {
      src <- c(src, 0L)
      fixed_val <- c(fixed_val, 0)
      th_idx <- c(th_idx, 0L)
    }
  }
  # 2. within means of the within-only (continuous) variables
  for (p in seq_len(nvar)) {
    if (!between_flag[p]) {
      src <- c(src, uni_start[p]) # mu
      fixed_val <- c(fixed_val, NA_real_)
    }
  }
  # 3. variances of the continuous variables (vw)
  for (p in seq_len(nvar)) {
    if (!ord[p]) {
      src <- c(src, uni_start[p] + 1L) # (mu, vw[, vb])
      fixed_val <- c(fixed_val, NA_real_)
    }
  }
  # 4. within covariances/'correlations' (vech, no diagonal): cw
  for (pp in seq_len(pstar)) {
    src <- c(src, ntheta + eta_start[pp]) # cw is the first free param
    fixed_val <- c(fixed_val, NA_real_)
  }

  # --- between block (between variables only) ---
  b_idx <- which(between_flag)
  # 5. between means: mu (continuous); 0 (ordinal, fixed)
  for (p in b_idx) {
    if (ord[p]) {
      src <- c(src, 0L)
      fixed_val <- c(fixed_val, 0)
    } else {
      src <- c(src, uni_start[p]) # mu
      fixed_val <- c(fixed_val, NA_real_)
    }
  }
  # 6. vech(Sigma_b, diagonal = TRUE): vb (diag) and cb (off-diagonal)
  # in lav_mat_vech (column-major) order over the between variables
  pair_pos <- matrix(0L, nvar, nvar)
  pp <- 0L
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      pp <- pp + 1L
      pair_pos[i, j] <- pp
    }
  }
  nb <- length(b_idx)
  for (jj in seq_len(nb)) {
    for (ii in jj:nb) {
      i <- b_idx[ii]; j <- b_idx[jj]
      if (i == j) {
        src <- c(src, uni_end[j]) # vb is the last univariate parameter
        fixed_val <- c(fixed_val, NA_real_)
      } else {
        pp <- pair_pos[i, j]
        stopifnot(pair_neta[pp] == 2L)
        src <- c(src, ntheta + eta_start[pp] + 1L) # cb
        fixed_val <- c(fixed_val, NA_real_)
      }
    }
  }

  # assemble
  nstat <- length(src)
  wls_obs <- numeric(nstat)
  free <- (src > 0L)
  wls_obs[free] <- m07$theta[src[free]]
  wls_obs[!free] <- fixed_val[!free]

  acov <- matrix(0, nstat, nstat)
  acov[free, free] <- m07$WLS.W[src[free], src[free]]

  # th vector (thresholds + zero means, as in the single-level layout)
  nth_all <- length(th_idx)
  th <- wls_obs[seq_len(nth_all)]

  list(
    wls_obs = wls_obs, acov = acov,
    free_idx = which(free),
    th = th, th_idx = th_idx
  )
}
