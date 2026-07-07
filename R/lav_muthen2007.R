# stage-wise estimation of the unrestricted TWO-LEVEL model for mixed
# continuous/ordinal data, and the asymptotic covariance matrix of the
# estimates (the two-level analogue of muthen1984())
#
# YR 2026 (two-level WLS, phase 4; covariates added in phase 7)
#
# stage 1: univariate two-level models per variable
#            continuous: (mu, vw, vb[, slopes]) -- lav_uvreg_2l_fit()
#            ordinal:    (tau_1..K, vb[, slopes]) -- lav_uvord_2l_fit()
#          (ordinal: within latent variance fixed to 1, no mean)
#          exogenous covariates (conditional.x): within-level covariates
#          x_w enter every regression; between-level covariates x_b only
#          the regressions of variables with a between part
# stage 2: bivariate two-level models per pair, univariate parameters
#          (incl. the slopes) fixed, estimating theta = (cw, cb) -- the
#          within and between RESIDUAL covariance -- jointly:
#            cont/cont: lav_bvreg_2l_cov_twostep_fit()
#            cont/ord:  lav_bvmix_2l_cov_twostep_fit()
#            ord/ord:   lav_bvord_2l_cor_twostep_fit()
# stage 3: sandwich ACOV of all estimates, with the CLUSTERS as the
#          independent units:
#
#              WLS.W = B^{-1} INNER B^{-T}
#
#          B = | A11   0  |  A11: block-diagonal univariate information
#              | A21  A22 |  A21: d(pair scores)/d(univariate params)
#                            A22: block-diagonal (2 x 2) pair information
#
#          all blocks estimated from the cluster-wise scores via the
#          information identities (cfr. lav_m84_* for the single-level
#          case); INNER = crossprod of the stacked cluster-wise scores.
#
# The parameter (statistic) vector s is ordered per variable first
# (mu/vw/vb or tau_1..K/vb, followed by the within + between slopes),
# then per pair (cw, cb), pairs in column-major lower-triangular order
# ((2,1), (3,1), ..., (3,2), ...). The mapping to the lavaan WLS.obs
# layout is done by the caller.

lav_muthen2007 <- function(y = NULL,
                           ov_names = NULL,
                           ov_types = NULL,
                           cluster_idx = NULL,
                           between_flag = NULL, # per variable: does it
                                                # have a between part?
                           x_w = NULL, # within-level covariates
                           x_b = NULL, # between-level covariates
                           xw_names = NULL,
                           xb_names = NULL,
                           ngh = 21L, # stage 1 + mixed pairs
                           ngh2 = 13L, # ordinal pairs (per dimension)
                           wls_w = TRUE) {
  y <- as.matrix(y)
  nvar <- NCOL(y)
  if (is.null(ov_names)) {
    ov_names <- colnames(y)
    if (is.null(ov_names)) {
      ov_names <- paste0("y", seq_len(nvar))
    }
  }
  ord_flag <- (ov_types == "ordered")
  if (is.null(between_flag)) {
    between_flag <- rep(TRUE, nvar)
  }
  if (any(ord_flag & !between_flag)) {
    lav_msg_stop(gettext(
      "within-only ORDINAL variables are not supported (yet) in the
      two-level stage-wise estimation."))
  }

  # covariates
  nexo_w <- if (is.null(x_w)) 0L else NCOL(x_w)
  nexo_b <- if (is.null(x_b)) 0L else NCOL(x_b)
  if (nexo_w > 0L && is.null(xw_names)) {
    xw_names <- colnames(x_w)
    if (is.null(xw_names)) {
      xw_names <- paste0("xw", seq_len(nexo_w))
    }
  }
  if (nexo_b > 0L && is.null(xb_names)) {
    xb_names <- colnames(x_b)
    if (is.null(xb_names)) {
      xb_names <- paste0("xb", seq_len(nexo_b))
    }
  }

  # ---- stage 1: univariate fits ----
  fit <- vector("list", nvar)
  for (p in seq_len(nvar)) {
    x_w_p <- if (nexo_w > 0L) x_w else NULL
    x_b_p <- if (nexo_b > 0L && between_flag[p]) x_b else NULL
    if (ord_flag[p]) {
      fit[[p]] <- lav_uvord_2l_fit(
        y = y[, p], x_w = x_w_p, x_b = x_b_p,
        cluster_idx = cluster_idx, ngh = ngh
      )
    } else {
      fit[[p]] <- lav_uvreg_2l_fit(
        y = y[, p], x_w = x_w_p, x_b = x_b_p,
        cluster_idx = cluster_idx,
        between = between_flag[p]
      )
    }
  }

  # univariate parameter blocks
  # free_idx_uni: positions of the FREE univariate parameters within
  # fit$theta (excludes the pinned vb of within-only variables);
  # sc_names_uni: the matching cluster-score column names
  npar_uni <- integer(nvar)
  free_idx_uni <- vector("list", nvar)
  sc_names_uni <- vector("list", nvar)
  th <- vector("list", nvar)
  th_idx <- vector("list", nvar)
  mu <- numeric(nvar)
  vw <- numeric(nvar)
  vb <- numeric(nvar)
  slopes_w <- matrix(0, nvar, nexo_w)
  slopes_b <- matrix(0, nvar, nexo_b)
  for (p in seq_len(nvar)) {
    f <- fit[[p]]
    sln <- lav_uvreg_2l_slope_names(f$nexo_w, f$nexo_b)
    if (ord_flag[p]) {
      th[[p]] <- f$theta[f$th_idx]
      th_idx[[p]] <- rep(p, length(th[[p]]))
      mu[p] <- 0
      vw[p] <- 1.0
      vb[p] <- f$theta[f$bvar_idx]
      free_idx_uni[[p]] <- c(f$th_idx, f$bvar_idx,
                             f$wslope_idx, f$bslope_idx)
      sc_names_uni[[p]] <- c(paste0("th", seq_along(f$th_idx)),
                             "bvar", sln)
    } else {
      th[[p]] <- f$theta[f$mu_idx]
      th_idx[[p]] <- 0L
      mu[p] <- f$theta[f$mu_idx]
      vw[p] <- f$theta[f$wvar_idx]
      vb[p] <- f$theta[f$bvar_idx]
      if (between_flag[p]) {
        free_idx_uni[[p]] <- c(1L, 2L, 3L, f$wslope_idx, f$bslope_idx)
        sc_names_uni[[p]] <- c("mu", "wvar", "bvar", sln)
      } else {
        # within-only variables: vb is fixed at 0 (not a parameter)
        free_idx_uni[[p]] <- c(1L, 2L, f$wslope_idx)
        sc_names_uni[[p]] <- c("mu", "wvar", sln)
      }
    }
    if (f$nexo_w > 0L) {
      slopes_w[p, ] <- f$theta[f$wslope_idx]
    }
    if (f$nexo_b > 0L) {
      slopes_b[p, ] <- f$theta[f$bslope_idx]
    }
    npar_uni[p] <- length(free_idx_uni[[p]])
  }

  # ---- stage 2: pairwise fits ----
  pstar <- nvar * (nvar - 1L) / 2L
  sigma_w <- diag(vw, nrow = nvar)
  sigma_b <- diag(vb, nrow = nvar)
  pair_theta <- matrix(0, pstar, 2L) # (cw, cb) per pair
  pair_vars <- matrix(0L, pstar, 2L) # (j, i) with j < i
  pair_type <- character(pstar)
  # pairs involving a within-only variable have cb fixed at 0: only cw
  # is a free parameter
  pair_neta <- integer(pstar)
  pp <- 0L
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      pp <- pp + 1L
      pair_vars[pp, ] <- c(j, i)
      if (!ord_flag[j] && !ord_flag[i]) {
        pair_type[pp] <- "bvreg"
        est <- lav_bvreg_2l_cov_twostep_fit(
          fit_y1 = fit[[j]], fit_y2 = fit[[i]],
          y1_name = ov_names[j], y2_name = ov_names[i]
        )
      } else if (ord_flag[j] && ord_flag[i]) {
        pair_type[pp] <- "bvord"
        est <- lav_bvord_2l_cor_twostep_fit(
          fit_y1 = fit[[j]], fit_y2 = fit[[i]], ngh = ngh2,
          y1_name = ov_names[j], y2_name = ov_names[i]
        )
      } else if (!ord_flag[j] && ord_flag[i]) {
        pair_type[pp] <- "bvmix"
        est <- lav_bvmix_2l_cov_twostep_fit(
          fit_y1 = fit[[j]], fit_y2 = fit[[i]], ngh = ngh,
          y1_name = ov_names[j], y2_name = ov_names[i]
        )
      } else {
        # ordinal (j) / continuous (i): continuous variable first
        pair_type[pp] <- "bvmix.swap"
        est <- lav_bvmix_2l_cov_twostep_fit(
          fit_y1 = fit[[i]], fit_y2 = fit[[j]], ngh = ngh,
          y1_name = ov_names[i], y2_name = ov_names[j]
        )
      }
      pair_theta[pp, ] <- est
      pair_neta[pp] <- if (between_flag[j] && between_flag[i]) 2L else 1L
      sigma_w[i, j] <- sigma_w[j, i] <- est[1]
      sigma_b[i, j] <- sigma_b[j, i] <- est[2]
    }
  }

  # the stage-wise parameter (statistic) vector
  s_uni <- unlist(lapply(seq_len(nvar), function(p) {
    fit[[p]]$theta[free_idx_uni[[p]]]
  }))
  s_pair <- unlist(lapply(seq_len(pstar), function(pp) {
    pair_theta[pp, seq_len(pair_neta[pp])]
  }))
  s_all <- c(s_uni, s_pair)
  par_names <- c(
    unlist(lapply(seq_len(nvar), function(p) {
      f <- fit[[p]]
      sl_nm <- character(0L)
      if (f$nexo_w > 0L) {
        sl_nm <- c(sl_nm, paste0(ov_names[p], "~", xw_names))
      }
      if (f$nexo_b > 0L) {
        sl_nm <- c(sl_nm, paste0(ov_names[p], "~", xb_names))
      }
      if (ord_flag[p]) {
        c(paste0(ov_names[p], "|t", seq_along(th[[p]])),
          paste0(ov_names[p], "~b~", ov_names[p]), sl_nm)
      } else if (between_flag[p]) {
        c(paste0(ov_names[p], "~1"),
          paste0(ov_names[p], "~w~", ov_names[p]),
          paste0(ov_names[p], "~b~", ov_names[p]), sl_nm)
      } else {
        c(paste0(ov_names[p], "~1"),
          paste0(ov_names[p], "~w~", ov_names[p]), sl_nm)
      }
    })),
    unlist(lapply(seq_len(pstar), function(pp) {
      j <- pair_vars[pp, 1]; i <- pair_vars[pp, 2]
      nm <- paste0(ov_names[i], "~w~", ov_names[j])
      if (pair_neta[pp] == 2L) {
        nm <- c(nm, paste0(ov_names[i], "~b~", ov_names[j]))
      }
      nm
    }))
  )
  names(s_all) <- par_names

  out <- list(
    FIT = fit,
    TH = th, TH.IDX = th_idx,
    MU = mu, SIGMA.W = sigma_w, SIGMA.B = sigma_b,
    SLOPES.W = slopes_w, SLOPES.B = slopes_b,
    pair.theta = pair_theta, pair.vars = pair_vars,
    pair.type = pair_type, pair.neta = pair_neta,
    between.flag = between_flag,
    nexo.w = nexo_w, nexo.b = nexo_b,
    theta = s_all, npar.uni = npar_uni,
    ord.flag = ord_flag,
    WLS.W = NULL, SC = NULL
  )

  if (!wls_w) {
    return(out)
  }

  # ---- stage 3: sandwich ----
  # global column layout of the univariate blocks
  uni_end <- cumsum(npar_uni)
  uni_start <- uni_end - npar_uni + 1L
  ntheta <- sum(npar_uni)
  eta_end <- cumsum(pair_neta)
  eta_start <- eta_end - pair_neta + 1L
  neta <- sum(pair_neta)
  npar_all <- ntheta + neta
  nclusters <- length(unique(cluster_idx))

  # univariate cluster-wise scores
  sc_all <- matrix(0, nclusters, npar_all)
  for (p in seq_len(nvar)) {
    if (ord_flag[p]) {
      sc_p <- lav_uvord_2l_sc(fit[[p]])
    } else {
      sc_p <- lav_uvreg_2l_sc(fit[[p]])
    }
    sc_p <- sc_p[, sc_names_uni[[p]], drop = FALSE]
    sc_all[, uni_start[p]:uni_end[p]] <- sc_p
  }

  # A11: block-diagonal univariate information (information identity)
  a11 <- matrix(0, ntheta, ntheta)
  a11_inv <- matrix(0, ntheta, ntheta)
  for (p in seq_len(nvar)) {
    idx <- uni_start[p]:uni_end[p]
    a11_p <- crossprod(sc_all[, idx, drop = FALSE])
    a11[idx, idx] <- a11_p
    a11_inv[idx, idx] <- lav_mat_sym_inverse(a11_p)
  }

  # pair-score column names for a member variable, in the same order as
  # its univariate block (sc_names_uni): base names + slope names
  m07_pair_var_names <- function(f, base_names, prefix_sl) {
    sln <- lav_uvreg_2l_slope_names(f$nexo_w, f$nexo_b)
    if (length(sln) == 0L) {
      return(base_names)
    }
    c(base_names, paste0(prefix_sl, sln))
  }

  # pairwise scores; A21 and A22
  a21 <- matrix(0, neta, ntheta)
  a22_inv <- matrix(0, neta, neta)
  for (pp in seq_len(pstar)) {
    j <- pair_vars[pp, 1]; i <- pair_vars[pp, 2]
    cw <- pair_theta[pp, 1]; cb <- pair_theta[pp, 2]
    type <- pair_type[pp]

    if (type == "bvreg") {
      sc_pair <- lav_bvreg_2l_sc(
        fit_y1 = fit[[j]], fit_y2 = fit[[i]], cw = cw, cb = cb
      )
      free_names <- c("wcov", "bcov")
      var1_names <- m07_pair_var_names(
        fit[[j]], c("mu_y1", "wvar_y1",
                    if (between_flag[j]) "bvar_y1"), "y1_")
      var2_names <- m07_pair_var_names(
        fit[[i]], c("mu_y2", "wvar_y2",
                    if (between_flag[i]) "bvar_y2"), "y2_")
    } else if (type == "bvord") {
      sc_pair <- lav_bvord_2l_sc(
        fit_y1 = fit[[j]], fit_y2 = fit[[i]], rho_w = cw, cb = cb,
        ngh = ngh2
      )
      free_names <- c("rho_w", "cb")
      var1_names <- m07_pair_var_names(
        fit[[j]], c(paste0("tau_a", seq_along(th[[j]])), "vb_a"), "a_")
      var2_names <- m07_pair_var_names(
        fit[[i]], c(paste0("tau_c", seq_along(th[[i]])), "vb_c"), "c_")
    } else if (type == "bvmix") {
      # y1 = continuous (j), y2 = ordinal (i)
      sc_pair <- lav_bvmix_2l_sc(
        fit_y1 = fit[[j]], fit_y2 = fit[[i]], cw = cw, cb = cb, ngh = ngh
      )
      free_names <- c("cw", "cb")
      var1_names <- m07_pair_var_names(
        fit[[j]], c("mu_a", "vw_a",
                    if (between_flag[j]) "vb_a"), "a_")
      var2_names <- m07_pair_var_names(
        fit[[i]], c(paste0("tau_c", seq_along(th[[i]])), "vb_c"), "c_")
    } else { # bvmix.swap: y1 = continuous (i), y2 = ordinal (j)
      sc_pair <- lav_bvmix_2l_sc(
        fit_y1 = fit[[i]], fit_y2 = fit[[j]], cw = cw, cb = cb, ngh = ngh
      )
      free_names <- c("cw", "cb")
      # variable j is the ordinal one here (the 'c' side of the module)
      var1_names <- m07_pair_var_names(
        fit[[j]], c(paste0("tau_c", seq_along(th[[j]])), "vb_c"), "c_")
      var2_names <- m07_pair_var_names(
        fit[[i]], c("mu_a", "vw_a",
                    if (between_flag[i]) "vb_a"), "a_")
    }

    free_idx <- match(free_names, colnames(sc_pair))
    var1_idx <- match(var1_names, colnames(sc_pair))
    var2_idx <- match(var2_names, colnames(sc_pair))
    if (anyNA(c(free_idx, var1_idx, var2_idx))) {
      wanted <- c(free_names, var1_names, var2_names)
      missing_names <- wanted[!wanted %in% colnames(sc_pair)]
      lav_msg_stop(gettextf(
        "internal error: pair score columns do not match the univariate
        parameter blocks (pair %1$s/%2$s, missing: %3$s; available: %4$s).",
        ov_names[j], ov_names[i],
        paste(missing_names, collapse = ", "),
        paste(colnames(sc_pair), collapse = ", ")))
    }

    # pairs with a within-only member: cb is fixed at 0, cw only
    free_idx <- free_idx[seq_len(pair_neta[pp])]

    eidx <- eta_start[pp]:eta_end[pp] # position within the eta block
    sc_all[, ntheta + eidx] <- sc_pair[, free_idx, drop = FALSE]

    # A22 block and its inverse
    a22_pp <- crossprod(sc_pair[, free_idx, drop = FALSE])
    a22_inv[eidx, eidx] <- solve(a22_pp)

    # A21 blocks: -E[d v / d theta] = E[v u_pair'] (information identity)
    g1 <- uni_start[j]:uni_end[j]
    g2 <- uni_start[i]:uni_end[i]
    a21[eidx, g1] <- crossprod(
      sc_pair[, free_idx, drop = FALSE],
      sc_pair[, var1_idx, drop = FALSE]
    )
    a21[eidx, g2] <- crossprod(
      sc_pair[, free_idx, drop = FALSE],
      sc_pair[, var2_idx, drop = FALSE]
    )
  }

  # block-triangular inverse of the bread
  a21_inv <- -a22_inv %*% a21 %*% a11_inv
  b_inv <- rbind(
    cbind(a11_inv, matrix(0, ntheta, neta)),
    cbind(a21_inv, a22_inv)
  )

  # meat: crossprod of the stacked cluster-wise scores
  inner <- crossprod(sc_all)

  # sandwich = ACOV of the stage-wise estimates
  wls_w_1 <- b_inv %*% inner %*% t(b_inv)
  wls_w_1 <- (wls_w_1 + t(wls_w_1)) / 2
  rownames(wls_w_1) <- colnames(wls_w_1) <- par_names

  out$WLS.W <- wls_w_1
  out$SC <- sc_all
  out$A11 <- a11
  out$A21 <- a21
  out
}
