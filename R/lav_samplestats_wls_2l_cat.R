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
  if (isTRUE(lavoptions$conditional.x)) {
    # note: we ALSO fill in the (residual, y-only) cov/mean entries:
    # several consumers (starting values, ...) look for
    # lavh1$implied$cov[[b]] regardless of conditional.x
    h1_implied <- list(
      cov = vector("list", ngroups * nlevels),
      mean = vector("list", ngroups * nlevels),
      res.cov = vector("list", ngroups * nlevels),
      res.int = vector("list", ngroups * nlevels),
      res.slopes = vector("list", ngroups * nlevels),
      cov.x = vector("list", ngroups * nlevels),
      mean.x = vector("list", ngroups * nlevels)
    )
  } else {
    h1_implied <- list(
      cov = vector("list", ngroups * nlevels),
      mean = vector("list", ngroups * nlevels)
    )
  }

  conditional_x <- isTRUE(lavoptions$conditional.x)

  for (g in seq_len(ngroups)) {
    lp <- lavdata@Lp[[g]]
    nobs <- lavsamplestats@nobs[[g]]

    if (!conditional_x && length(lavsamplestats@x.idx[[g]]) > 0L) {
      lav_msg_stop(gettext(
        "fixed.x = TRUE with conditional.x = FALSE is not supported (yet)
        for two-level (D)WLS estimation."))
    }

    # y variables per level (excluding exogenous covariates), in the
    # block variable order; the covariates (if any) come from the Lp
    # within.x/between.x index sets ('split' covariates -- appearing at
    # both levels -- have already been moved to the y set).
    # NOTE: all Lp index sets refer to the 'tilde' universe (lp$ov.names
    # = the y variables followed by the x variables), which is also the
    # column order of lavdata@X[[g]] for multilevel data
    ov_names_all <- lp$ov.names
    x_w <- x_b <- NULL
    xw_names <- xb_names <- NULL
    if (conditional_x) {
      y_idx1 <- lp$ov.y.idx[[1]]
      y_idx2 <- lp$ov.y.idx[[2]]
      xw_idx <- lp$within.x.idx[[1]]
      xb_idx <- lp$between.x.idx[[2]]
      if (length(xw_idx) > 0L) {
        x_w <- lavdata@X[[g]][, xw_idx, drop = FALSE]
        xw_names <- ov_names_all[xw_idx]
      }
      if (length(xb_idx) > 0L) {
        x_b <- lavdata@X[[g]][, xb_idx, drop = FALSE]
        xb_names <- ov_names_all[xb_idx]
      }
      # 'split' covariates (exogenous, but appearing at both levels) are
      # not supported (yet): every covariate must be within-only or
      # between-only
      split_x <- setdiff(lavdata@ov.names.x[[g]], c(xw_names, xb_names))
      if (length(split_x) > 0L) {
        lav_msg_stop(gettextf(
          "two-level (D)WLS estimation does not support exogenous
          covariates that appear at both levels (yet): %s. Each covariate
          should be used at a single level only (cfr. within/between
          covariates); alternatively, use fixed.x = FALSE.",
          lav_msg_view(split_x, "none")))
      }
    } else {
      y_idx1 <- lp$ov.idx[[1]]
      y_idx2 <- lp$ov.idx[[2]]
    }
    ov_names_1 <- ov_names_all[y_idx1]
    ov_names_2 <- ov_names_all[y_idx2]
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
      y = lavdata@X[[g]][, y_idx1, drop = FALSE],
      ov_names = ov_names_1,
      ov_types = ov_types,
      cluster_idx = lp$cluster.idx[[2]],
      between_flag = between_flag,
      x_w = x_w, x_b = x_b,
      xw_names = xw_names, xb_names = xb_names,
      ngh = lavoptions$integration.ngh,
      ngh2 = 13L,
      wls_w = TRUE
    )

    # map the stage-wise statistic vector + ACOV to the lavaan layout
    out <- lav_m07_repackage(m07 = m07, ov_types = ov_types,
                             between_flag = between_flag)

    wls_obs[[g]] <- out$wls_obs
    nacov[[g]] <- nobs * out$acov

    w <- lav_samp_wls_v_2l(
      gamma_1 = nacov[[g]], keep = out$free_idx,
      estimator = estimator, nclusters = lp$nclusters[[2]], group = g
    )
    wls_v[[g]] <- w$wls_v
    if (!is.null(w$wls_vd)) {
      wls_vd[[g]] <- w$wls_vd
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
    w_pos <- (g - 1) * nlevels + 1L
    b_pos <- (g - 1) * nlevels + 2L
    if (conditional_x) {
      h1_implied$cov[[w_pos]] <- m07$SIGMA.W
      h1_implied$cov[[b_pos]] <-
        m07$SIGMA.B[b_idx, b_idx, drop = FALSE]
      h1_implied$mean[[w_pos]] <- mu_w1
      h1_implied$mean[[b_pos]] <- m07$MU[b_idx]
      h1_implied$res.cov[[w_pos]] <- m07$SIGMA.W
      h1_implied$res.cov[[b_pos]] <-
        m07$SIGMA.B[b_idx, b_idx, drop = FALSE]
      h1_implied$res.int[[w_pos]] <- mu_w1
      h1_implied$res.int[[b_pos]] <- m07$MU[b_idx]
      h1_implied$res.slopes[[w_pos]] <- m07$SLOPES.W
      h1_implied$res.slopes[[b_pos]] <-
        m07$SLOPES.B[b_idx, , drop = FALSE]
      # observed covariate moments (within: observation level; between:
      # cluster level), with the covariate names as dimnames
      if (!is.null(x_w)) {
        n1 <- NROW(x_w)
        cx1 <- cov(x_w) * (n1 - 1) / n1
        dimnames(cx1) <- list(xw_names, xw_names)
        h1_implied$cov.x[[w_pos]] <- cx1
        h1_implied$mean.x[[w_pos]] <- setNames(colMeans(x_w), xw_names)
      } else {
        h1_implied$cov.x[[w_pos]] <- matrix(0, 0L, 0L)
        h1_implied$mean.x[[w_pos]] <- numeric(0L)
      }
      if (!is.null(x_b)) {
        xbj <- x_b[!duplicated(lp$cluster.idx[[2]]), , drop = FALSE]
        nj <- NROW(xbj)
        cx2 <- cov(xbj) * (nj - 1) / nj
        dimnames(cx2) <- list(xb_names, xb_names)
        h1_implied$cov.x[[b_pos]] <- cx2
        h1_implied$mean.x[[b_pos]] <- setNames(colMeans(xbj), xb_names)
      } else {
        h1_implied$cov.x[[b_pos]] <- matrix(0, 0L, 0L)
        h1_implied$mean.x[[b_pos]] <- numeric(0L)
      }
    } else {
      h1_implied$cov[[w_pos]] <- m07$SIGMA.W
      h1_implied$cov[[b_pos]] <-
        m07$SIGMA.B[b_idx, b_idx, drop = FALSE]
      h1_implied$mean[[w_pos]] <- mu_w1
      h1_implied$mean[[b_pos]] <- m07$MU[b_idx]
    }
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
#   2. within means/intercepts of the within-only (continuous) variables
#   3. vec(Pi_w): within slopes, column-major over (variable, covariate)
#   4. within variances of the continuous variables
#   5. within covariances (vech, no diagonal)
#   6. between means/intercepts (continuous: mu; ordinal: fixed 0) --
#      between variables only
#   7. vec(Pi_b): between slopes, column-major -- between variables only
#   8. vech(Sigma_b, diagonal = TRUE) -- between variables only
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
  nexo_w <- if (is.null(m07$nexo.w)) 0L else m07$nexo.w
  nexo_b <- if (is.null(m07$nexo.b)) 0L else m07$nexo.b

  # per-variable offsets (within the univariate block) of the 'base'
  # parameters and the slopes; the univariate block layout is
  #   ordinal:              (tau_1..K, vb, wsl..., bsl...)
  #   continuous (between): (mu, vw, vb, wsl..., bsl...)
  #   continuous (within):  (mu, vw, wsl...)
  nth_p <- integer(nvar) # thresholds per variable (0 = continuous)
  vb_off <- integer(nvar) # offset of vb (NA-like 0 if absent)
  wsl_off <- integer(nvar) # offset of the last param before wsl
  for (p in seq_len(nvar)) {
    if (ord[p]) {
      nth_p[p] <- npar_uni[p] - 1L - nexo_w -
        (if (between_flag[p]) nexo_b else 0L)
      vb_off[p] <- nth_p[p] + 1L
      wsl_off[p] <- nth_p[p] + 1L
    } else if (between_flag[p]) {
      vb_off[p] <- 3L
      wsl_off[p] <- 3L
    } else {
      vb_off[p] <- 0L
      wsl_off[p] <- 2L
    }
  }

  # source positions (in the stage-wise vector) of each target entry;
  # 0 = structurally-fixed entry (value stored in fixed_val)
  src <- integer(0L)
  fixed_val <- numeric(0L)
  th_idx <- integer(0L)

  # --- within block ---
  # 1. th: thresholds (ordinal); fixed 0 slots (continuous)
  for (p in seq_len(nvar)) {
    if (ord[p]) {
      src <- c(src, uni_start[p] - 1L + seq_len(nth_p[p]))
      fixed_val <- c(fixed_val, rep(NA_real_, nth_p[p]))
      th_idx <- c(th_idx, rep(p, nth_p[p]))
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
  # 3. vec(Pi_w): all variables for covariate 1, then covariate 2, ...
  for (q in seq_len(nexo_w)) {
    for (p in seq_len(nvar)) {
      src <- c(src, uni_start[p] - 1L + wsl_off[p] + q)
      fixed_val <- c(fixed_val, NA_real_)
    }
  }
  # 4. variances of the continuous variables (vw)
  for (p in seq_len(nvar)) {
    if (!ord[p]) {
      src <- c(src, uni_start[p] + 1L) # (mu, vw, ...)
      fixed_val <- c(fixed_val, NA_real_)
    }
  }
  # 5. within covariances/'correlations' (vech, no diagonal): cw
  for (pp in seq_len(pstar)) {
    src <- c(src, ntheta + eta_start[pp]) # cw is the first free param
    fixed_val <- c(fixed_val, NA_real_)
  }

  # --- between block (between variables only) ---
  b_idx <- which(between_flag)
  # 6. between means: mu (continuous); 0 (ordinal, fixed)
  for (p in b_idx) {
    if (ord[p]) {
      src <- c(src, 0L)
      fixed_val <- c(fixed_val, 0)
    } else {
      src <- c(src, uni_start[p]) # mu
      fixed_val <- c(fixed_val, NA_real_)
    }
  }
  # 7. vec(Pi_b) over the between variables
  for (q in seq_len(nexo_b)) {
    for (p in b_idx) {
      src <- c(src, uni_start[p] - 1L + wsl_off[p] + nexo_w + q)
      fixed_val <- c(fixed_val, NA_real_)
    }
  }
  # 8. vech(Sigma_b, diagonal = TRUE): vb (diag) and cb (off-diagonal)
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
        src <- c(src, uni_start[j] - 1L + vb_off[j])
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
