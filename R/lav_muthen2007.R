# stage-wise estimation of the unrestricted TWO-LEVEL model for mixed
# continuous/ordinal data, and the asymptotic covariance matrix of the
# estimates (the two-level analogue of muthen1984())
#
# YR 2026 (two-level WLS, phase 4)
#
# stage 1: univariate two-level models per variable
#            continuous: (mu, vw, vb)        -- lav_uvreg_2l_fit()
#            ordinal:    (tau_1..K, vb)      -- lav_uvord_2l_fit()
#          (ordinal: within latent variance fixed to 1, no mean)
# stage 2: bivariate two-level models per pair, univariate parameters
#          fixed, estimating theta = (cw, cb) -- the within and between
#          covariance -- jointly:
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
# (mu/vw/vb or tau_1..K/vb), then per pair (cw, cb), pairs in
# column-major lower-triangular order ((2,1), (3,1), ..., (3,2), ...).
# The mapping to the lavaan WLS.obs layout is done by the caller.

lav_muthen2007 <- function(y = NULL,
                           ov_names = NULL,
                           ov_types = NULL,
                           cluster_idx = NULL,
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

  # ---- stage 1: univariate fits ----
  fit <- vector("list", nvar)
  for (p in seq_len(nvar)) {
    if (ord_flag[p]) {
      fit[[p]] <- lav_uvord_2l_fit(
        y = y[, p], cluster_idx = cluster_idx, ngh = ngh
      )
    } else {
      fit[[p]] <- lav_uvreg_2l_fit(
        y = y[, p], cluster_idx = cluster_idx
      )
    }
  }

  # univariate parameter blocks
  npar_uni <- integer(nvar)
  th <- vector("list", nvar)
  th_idx <- vector("list", nvar)
  mu <- numeric(nvar)
  vw <- numeric(nvar)
  vb <- numeric(nvar)
  for (p in seq_len(nvar)) {
    if (ord_flag[p]) {
      th[[p]] <- fit[[p]]$theta[fit[[p]]$th_idx]
      th_idx[[p]] <- rep(p, length(th[[p]]))
      mu[p] <- 0
      vw[p] <- 1.0
      vb[p] <- fit[[p]]$theta[fit[[p]]$bvar_idx]
      npar_uni[p] <- length(th[[p]]) + 1L
    } else {
      th[[p]] <- fit[[p]]$theta[fit[[p]]$mu_idx]
      th_idx[[p]] <- 0L
      mu[p] <- fit[[p]]$theta[fit[[p]]$mu_idx]
      vw[p] <- fit[[p]]$theta[fit[[p]]$wvar_idx]
      vb[p] <- fit[[p]]$theta[fit[[p]]$bvar_idx]
      npar_uni[p] <- 3L
    }
  }

  # ---- stage 2: pairwise fits ----
  pstar <- nvar * (nvar - 1L) / 2L
  sigma_w <- diag(vw, nrow = nvar)
  sigma_b <- diag(vb, nrow = nvar)
  pair_theta <- matrix(0, pstar, 2L) # (cw, cb) per pair
  pair_vars <- matrix(0L, pstar, 2L) # (j, i) with j < i
  pair_type <- character(pstar)
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
      sigma_w[i, j] <- sigma_w[j, i] <- est[1]
      sigma_b[i, j] <- sigma_b[j, i] <- est[2]
    }
  }

  # the stage-wise parameter (statistic) vector
  s_uni <- unlist(lapply(seq_len(nvar), function(p) {
    if (ord_flag[p]) {
      c(th[[p]], vb[p])
    } else {
      c(mu[p], vw[p], vb[p])
    }
  }))
  s_all <- c(s_uni, t(pair_theta)) # pairs appended as (cw, cb) each
  par_names <- c(
    unlist(lapply(seq_len(nvar), function(p) {
      if (ord_flag[p]) {
        c(paste0(ov_names[p], "|t", seq_along(th[[p]])),
          paste0(ov_names[p], "~b~", ov_names[p]))
      } else {
        c(paste0(ov_names[p], "~1"),
          paste0(ov_names[p], "~w~", ov_names[p]),
          paste0(ov_names[p], "~b~", ov_names[p]))
      }
    })),
    unlist(lapply(seq_len(pstar), function(pp) {
      j <- pair_vars[pp, 1]; i <- pair_vars[pp, 2]
      c(paste0(ov_names[i], "~w~", ov_names[j]),
        paste0(ov_names[i], "~b~", ov_names[j]))
    }))
  )
  names(s_all) <- par_names

  out <- list(
    FIT = fit,
    TH = th, TH.IDX = th_idx,
    MU = mu, SIGMA.W = sigma_w, SIGMA.B = sigma_b,
    pair.theta = pair_theta, pair.vars = pair_vars,
    pair.type = pair_type,
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
  neta <- 2L * pstar
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
      free_idx <- match(c("wcov", "bcov"), colnames(sc_pair))
      var1_idx <- match(c("mu_y1", "wvar_y1", "bvar_y1"), colnames(sc_pair))
      var2_idx <- match(c("mu_y2", "wvar_y2", "bvar_y2"), colnames(sc_pair))
    } else if (type == "bvord") {
      sc_pair <- lav_bvord_2l_sc(
        fit_y1 = fit[[j]], fit_y2 = fit[[i]], rho_w = cw, cb = cb,
        ngh = ngh2
      )
      nth_j <- length(th[[j]]); nth_i <- length(th[[i]])
      free_idx <- match(c("rho_w", "cb"), colnames(sc_pair))
      var1_idx <- c(seq_len(nth_j), # tau_a
                    match("vb_a", colnames(sc_pair)))
      var2_idx <- c(nth_j + seq_len(nth_i), # tau_c
                    match("vb_c", colnames(sc_pair)))
    } else if (type == "bvmix") {
      # y1 = continuous (j), y2 = ordinal (i)
      sc_pair <- lav_bvmix_2l_sc(
        fit_y1 = fit[[j]], fit_y2 = fit[[i]], cw = cw, cb = cb, ngh = ngh
      )
      nth_i <- length(th[[i]])
      free_idx <- match(c("cw", "cb"), colnames(sc_pair))
      var1_idx <- match(c("mu_a", "vw_a", "vb_a"), colnames(sc_pair))
      var2_idx <- c(3L + seq_len(nth_i), # tau_c
                    match("vb_c", colnames(sc_pair)))
    } else { # bvmix.swap: y1 = continuous (i), y2 = ordinal (j)
      sc_pair <- lav_bvmix_2l_sc(
        fit_y1 = fit[[i]], fit_y2 = fit[[j]], cw = cw, cb = cb, ngh = ngh
      )
      nth_j <- length(th[[j]])
      free_idx <- match(c("cw", "cb"), colnames(sc_pair))
      # variable j is the ordinal one here (the 'c' side of the module)
      var1_idx <- c(3L + seq_len(nth_j), match("vb_c", colnames(sc_pair)))
      var2_idx <- match(c("mu_a", "vw_a", "vb_a"), colnames(sc_pair))
    }

    eidx <- (pp - 1L) * 2L + 1:2 # position within the eta block
    sc_all[, ntheta + eidx] <- sc_pair[, free_idx, drop = FALSE]

    # A22 block (2 x 2) and its inverse
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
