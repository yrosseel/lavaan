# delta-method standard errors for the noniterative estimators
#
# YR 20 Jul 2026: first version (wired for estimator = "MGM")
#
# Every noniterative estimator is a smooth map of the per-group sample
# moments, x = g(mean_b, vech(S_b)), so the covariance matrix of the
# estimates is the block-diagonal moment sandwich
#
#    vcov = sum_b J_b Gamma_b J_b' / nobs_b
#
# with J_b the Jacobian of the complete estimation map over block b's
# moments, and Gamma_b the asymptotic covariance of those moments. Two
# Gamma flavors are available:
#  - "nt" (default): the normal-theory Gamma;
#  - "adf": the distribution-free (empirical) Gamma. Because the casewise
#    influence of a moment-based estimator is J times the casewise moment
#    contribution, this flavor equals the infinitesimal-jackknife
#    covariance of the estimator.
#
# The Jacobian is computed numerically (numDeriv over the stacked moment
# vector, as in the JS fallback path), so ALL branches of the estimation
# map -- purging, pooled equality constraints, cross-loadings, the std.lv
# rescaling, the restricted PSI fit, the second-stage RLS weights -- are
# propagated automatically. The inner solvers converge to tolerances far
# below the numDeriv step size, so the finite differences stay clean; the
# regime switches (pd correction, clipping, bounds) are measure-zero
# kinks -- when a bound is ACTIVE at the estimate, the corresponding
# standard error degenerates towards zero.
#
# Other noniterative estimators (FABIN, BENTLER1982) can reuse the engine
# by adding a branch to lav_noniter_estimate_from_vec().

# re-run the complete estimation map on a (perturbed) moment vector;
# vec = NULL evaluates the map at the sample moments
lav_noniter_estimate_from_vec <- function(vec = NULL,
                                          lavmodel = NULL,
                                          lavsamplestats = NULL,
                                          lavpartable = NULL,
                                          lavdata = NULL,
                                          lavh1 = NULL,
                                          lavoptions = NULL) {
  if (!is.null(vec)) {
    implied <- lav_vec_to_implied(vec, lavmodel = lavmodel)
    # stage 1 reads lavsamplestats@cov/@mean; the second stage (and the
    # mean solve) read lavh1$implied when available: both must move
    # together
    for (b in seq_len(lavmodel@nblocks)) {
      lavsamplestats@cov[[b]] <- implied$cov[[b]]
      if (lavmodel@meanstructure) {
        lavsamplestats@mean[[b]] <- implied$mean[[b]]
      }
    }
    if (!is.null(lavh1$implied$cov)) {
      lavh1$implied$cov <- implied$cov
      if (lavmodel@meanstructure) {
        lavh1$implied$mean <- implied$mean
      }
    }
  }

  if (lavoptions$estimator == "MGM") {
    x <- lav_cfa_guttman1952_internal(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats,
      lavpartable = lavpartable, lavdata = lavdata,
      lavh1 = lavh1, lavoptions = lavoptions
    )
  } else {
    lav_msg_stop(gettextf(
      "delta-method standard errors are not available for estimator %s.",
      lavoptions$estimator))
  }

  as.numeric(x)
}

# ---------------------------------------------------------------------
# ANALYTIC Jacobian of the MGM estimation map over the sample moments
# ---------------------------------------------------------------------
#
# Covered branches (everything else returns NULL, and the caller falls
# back to the numerical Jacobian):
#  - stage 1, classic (simple structure, marker = 1, saturated PSI):
#      theta_p = clip( s_pp - mean_{a<b} s_pa s_pb / s_ab )   (Spearman)
#      W       = S* - rho_eff * diag(theta)                   (pd corr.)
#      T       = X' W X,   U = W X T^{-1},   u_f = U[marker_f, f]
#      lambda_pf = U_pf / u_f,   psi_fg = T_fg u_f u_g
#      psi       = M W M', M = (L'KL)^{-1} L'K  (psi.mapping = TRUE)
#    (the sum-score normalization cancels exactly)
#  - stage 1, pooled (equality ties among the loadings, or quadprog):
#      lstar_pf = (WX)_pf / sqrt(T_ff)  (normalized structure coeff.)
#      tied cells averaged; markers pooled over connected components;
#      lambda = value / mhat;  psi_fg = phi_fg mhat_f mhat_g
#    including the JOINT pd-correction shrink (one factor over all
#    groups -- its derivative has cross-block support)
#  - residual covariances: the purge s0_ij = mean_(a,b) s_ia s_jb / s_ab
#    is a smooth map of S; its chain P is composed with the stage-1
#    derivatives (which are evaluated on the purged matrix)
#  - stage 2 (mgm.varcov = ULS/GLS/RLS/2RLS): reuses the JS stage-2
#    derivative (lav_sem_js_deriv_varcov: direct moment part, chain
#    through the loadings, RLS implicit-function step at the fixed
#    point), with the MGM loading derivatives as dx
#  - mean structure: reuses the JS joint-GLS mean-solve derivative
#    (lav_sem_js_deriv_mean_wls), so free intercepts, tied intercepts
#    (scalar invariance) and free latent means are all covered
#
# The pd-correction root rho is a generalized eigenvalue of the pencil
# (S, Theta): d rho = v' (dS - rho dTheta) v / (v' Theta v). Regime
# switches (Spearman clips, theta bounds, pd correction active or not,
# the zero-theta guard of the mapping) are frozen at the point estimate
# -- exactly what the numerical Jacobian sees locally. A final safety
# net compares the predicted estimates against the actual ones; any
# mismatch (an active bound clip, or a branch the coverage checks
# missed) also triggers the numerical fallback.

# Spearman theta chain for one block: the (clipped) initial residual
# variances and their derivative with respect to the local vech(S)
lav_cfa_mgm_jac_theta <- function(s_mat = NULL, m_b = NULL) {
  nvar <- nrow(m_b)
  nfac <- ncol(m_b)
  pstar <- nvar * (nvar + 1L) / 2L
  vpos <- function(i, j) {
    if (i < j) {
      tmp <- i
      i <- j
      j <- tmp
    }
    (j - 1L) * nvar - ((j - 1L) * j) %/% 2L + i
  }
  theta0 <- numeric(nvar)
  dtheta <- matrix(0, nvar, pstar)
  for (f in seq_len(nfac)) {
    f_idx <- which(m_b[, f] == 1L)
    for (p in f_idx) {
      others <- setdiff(f_idx, p)
      prs <- utils::combn(others, 2L)
      nk <- ncol(prs)
      asum <- 0
      drow <- numeric(pstar)
      for (r in seq_len(nk)) {
        a <- prs[1L, r]
        bb <- prs[2L, r]
        rat <- s_mat[p, a] * s_mat[p, bb] / s_mat[a, bb]
        asum <- asum + rat
        drow[vpos(p, a)] <- drow[vpos(p, a)] + s_mat[p, bb] / s_mat[a, bb]
        drow[vpos(p, bb)] <- drow[vpos(p, bb)] + s_mat[p, a] / s_mat[a, bb]
        drow[vpos(a, bb)] <- drow[vpos(a, bb)] - rat / s_mat[a, bb]
      }
      h2 <- asum / (nk * s_mat[p, p])
      if (h2 < -0.05 || h2 > 1.20) {
        # Spearman clip active: theta = (1 - c) s_pp
        h2c <- if (h2 < -0.05) -0.05 else 1.20
        theta0[p] <- (1 - h2c) * s_mat[p, p]
        dtheta[p, vpos(p, p)] <- 1 - h2c
      } else {
        # theta = s_pp - asum / nk
        theta0[p] <- s_mat[p, p] - asum / nk
        dtheta[p, ] <- -drow / nk
        dtheta[p, vpos(p, p)] <- dtheta[p, vpos(p, p)] + 1
      }
      # worker bounds: clip to [0, s_pp]
      if (theta0[p] < 0) {
        theta0[p] <- 0
        dtheta[p, ] <- 0
      } else if (theta0[p] > s_mat[p, p]) {
        theta0[p] <- s_mat[p, p]
        dtheta[p, ] <- 0
        dtheta[p, vpos(p, p)] <- 1
      }
    }
  }
  list(theta0 = theta0, dtheta = dtheta)
}

# smallest root of the pencil |S - rho Theta| = 0, plus the full-space
# eigenvector needed for its derivative; NULL when the root is
# (near-)multiple or the zero-theta Schur step fails
lav_cfa_mgm_jac_root <- function(s_mat = NULL, theta0 = NULL) {
  nvar <- nrow(s_mat)
  zero_idx <- which(abs(theta0) < sqrt(.Machine$double.eps))
  if (length(zero_idx) == nvar) {
    return(NULL) # no finite roots (theta ~ 0): let the caller fall back
  }
  pos <- setdiff(seq_len(nvar), zero_idx)
  a_mat <- s_mat[pos, pos, drop = FALSE]
  s22i_s21 <- NULL
  if (length(zero_idx) > 0L) {
    s22i_s21 <- try(solve(
      s_mat[zero_idx, zero_idx, drop = FALSE],
      s_mat[zero_idx, pos, drop = FALSE]
    ), silent = TRUE)
    if (inherits(s22i_s21, "try-error")) {
      return(NULL)
    }
    a_mat <- a_mat - s_mat[pos, zero_idx, drop = FALSE] %*% s22i_s21
  }
  ld <- 1 / sqrt(theta0[pos])
  ee <- eigen(t(ld * a_mat) * ld, symmetric = TRUE)
  np <- length(pos)
  rho <- ee$values[np]
  # a (near-)multiple smallest root is not differentiable
  if (np > 1L && (ee$values[np - 1L] - rho) < 1e-8 * max(1, abs(rho))) {
    return(NULL)
  }
  # full-space pencil eigenvector: S v = rho Theta v
  v <- numeric(nvar)
  v[pos] <- ld * ee$vectors[, np]
  if (length(zero_idx) > 0L) {
    v[zero_idx] <- -drop(s22i_s21 %*% v[pos])
  }
  list(rho = rho, v = v, denom = sum(theta0 * v * v))
}

# derivative of the pencil root over the local vech(S):
# d rho = v'(dS - rho dTheta)v / (v' Theta v)
lav_cfa_mgm_jac_droot <- function(root = NULL, dtheta = NULL) {
  v <- root$v
  nvar <- length(v)
  vv <- tcrossprod(v)
  dquad <- 2 * lav_mat_vech(vv)
  dquad[lav_mat_diagh_idx(nvar)] <- diag(vv)
  (dquad - drop((v * v) %*% dtheta) * root$rho) / root$denom
}

lav_cfa_mgm_jacobian <- function(lavmodel = NULL, lavsamplestats = NULL,
                                 lavoptions = NULL, lavpartable = NULL,
                                 implied0 = NULL) {
  ea <- lavoptions$estimator.args
  psi_mapping <- isTRUE(ea[["psi.mapping"]])
  quadprog <- isTRUE(ea[["quadprog"]])
  zae <- ea[["zero.after.efa"]]
  zae <- is.null(zae) || isTRUE(zae)
  if (lavmodel@group.w.free) {
    return(NULL)
  }

  free <- lavpartable$free
  op <- lavpartable$op
  lam_all <- which(op == "=~" & free > 0L)
  has_ties <- any(duplicated(free[lam_all]))
  pooled <- has_ties || quadprog
  # not covered: psi.mapping in the pooled branch, or over an unzeroed
  # (EFA-style) lambda
  if (psi_mapping && (pooled || !zae)) {
    return(NULL)
  }

  lavpta <- lav_pt_attributes(lavpartable)
  nblocks <- lav_pt_nblocks(lavpartable)
  pt_block <- if (!is.null(lavpartable$block)) {
    lavpartable$block
  } else {
    rep(1L, length(op))
  }

  co <- lav_sem_js_deriv_coords(lavmodel)
  m_total <- co$m_total
  npar <- lavmodel@nx.free
  dx <- matrix(0, nrow = npar, ncol = m_total)
  # the values these formulas predict; compared against the actual
  # estimates at the end
  x_pred <- rep(NA_real_, npar)

  # --------------------------------------------------------------------
  # per-block structures and coverage checks
  # --------------------------------------------------------------------
  bi_list <- vector("list", nblocks)
  res_cov_flag <- FALSE
  pin_flag <- FALSE # any ~~ cell fixed to a nonzero value
  psi_ids_all <- integer(0L)
  th_ids_all <- integer(0L)
  for (b in seq_len(nblocks)) {
    ov_names <- lavpta$vnames$ov[[b]]
    lv_names <- lavpta$vnames$lv.regular[[b]]
    nvar <- length(ov_names)
    nfac <- length(lv_names)
    pstar <- nvar * (nvar + 1L) / 2L

    # markers fixed to 1 for every factor
    marker_idx <- integer(nfac)
    for (f in seq_len(nfac)) {
      fx <- which(op == "=~" & pt_block == b &
        lavpartable$lhs == lv_names[f] & free == 0L &
        !is.na(lavpartable$ustart) & lavpartable$ustart != 0)
      if (length(fx) != 1L || lavpartable$ustart[fx] != 1) {
        return(NULL) # std.lv / fixed value c != 1 / multiple markers
      }
      marker_idx[f] <- match(lavpartable$rhs[fx], ov_names)
    }
    lam_pt_idx <- which(op == "=~" & pt_block == b & free > 0L)
    lam_ov <- match(lavpartable$rhs[lam_pt_idx], ov_names)
    lam_fac <- match(lavpartable$lhs[lam_pt_idx], lv_names)
    if (anyNA(lam_ov) || anyNA(lam_fac)) {
      return(NULL)
    }
    m_b <- matrix(0L, nvar, nfac)
    m_b[(seq_len(nfac) - 1L) * nvar + marker_idx] <- 1L
    m_b[(lam_fac - 1L) * nvar + lam_ov] <- 1L
    if (any(rowSums(m_b) > 1L)) {
      return(NULL) # cross-loadings
    }
    if (any(colSums(m_b) < 3L)) {
      return(NULL) # masked/external tetrads
    }
    # saturated PSI, all cells free
    psi_rows_pt <- which(op == "~~" & pt_block == b &
      lavpartable$lhs %in% lv_names & lavpartable$rhs %in% lv_names)
    if (length(psi_rows_pt) != nfac * (nfac + 1L) / 2L ||
        any(free[psi_rows_pt] == 0L)) {
      return(NULL) # fixed variance/covariance -> rescale or psi fit
    }
    psi_ids_all <- c(psi_ids_all, free[psi_rows_pt])
    # theta diagonal
    th_pt_idx <- which(op == "~~" & pt_block == b &
      lavpartable$lhs %in% ov_names &
      lavpartable$lhs == lavpartable$rhs & free > 0L)
    th_ids_all <- c(th_ids_all, free[th_pt_idx])
    if (any(op == "~~" & pt_block == b &
            lavpartable$lhs %in% ov_names &
            lavpartable$lhs == lavpartable$rhs & free == 0L &
            !is.na(lavpartable$ustart) & lavpartable$ustart != 0)) {
      pin_flag <- TRUE # fixed nonzero residual variance
    }

    # residual covariances -> purge chain
    s_raw <- implied0$cov[[b]]
    s_use <- s_raw
    p_chain <- NULL
    cov_idx <- which(op == "~~" & pt_block == b &
      lavpartable$lhs != lavpartable$rhs &
      (free > 0L | (!is.na(lavpartable$ustart) & lavpartable$ustart != 0)) &
      lavpartable$lhs %in% ov_names & lavpartable$rhs %in% ov_names)
    if (length(cov_idx) > 0L) {
      res_cov_flag <- TRUE
      if (any(free[cov_idx] == 0L)) {
        pin_flag <- TRUE # fixed nonzero residual covariance
      }
      contam <- lav_sem_js_contam(
        lavpartable = lavpartable, pt_block = pt_block, b = b,
        ov_names = ov_names, lavpta = lavpta
      )
      fac_of <- rep(NA_integer_, nvar)
      for (f in seq_len(nfac)) {
        fac_of[which(m_b[, f] == 1L)] <- f
      }
      pairs <- data.frame(
        i = match(lavpartable$lhs[cov_idx], ov_names),
        j = match(lavpartable$rhs[cov_idx], ov_names),
        value = ifelse(free[cov_idx] > 0L, NA_real_,
                       lavpartable$ustart[cov_idx])
      )
      s_use <- lav_cfa_mgm_purge(
        s = s_raw, fac_of = fac_of, pairs = pairs, contam = contam,
        return_anchors = TRUE
      )
      if (length(attr(s_use, "failed")) > 0L) {
        return(NULL) # unpurgeable pair (masked-spearman path)
      }
      anchors <- attr(s_use, "anchors")
      attr(s_use, "failed") <- NULL
      attr(s_use, "anchors") <- NULL
      # the purge chain: replaced cells are means of tetrad ratios of
      # the RAW matrix; everything else is the identity
      vpos <- function(i, j) {
        if (i < j) {
          tmp <- i
          i <- j
          j <- tmp
        }
        (j - 1L) * nvar - ((j - 1L) * j) %/% 2L + i
      }
      p_chain <- diag(pstar)
      for (r in seq_len(nrow(pairs))) {
        if (!is.na(pairs$value[r])) {
          next # fixed pair: s0 = s - value, identity row
        }
        i <- pairs$i[r]
        j <- pairs$j[r]
        an <- anchors[[r]]
        nk <- ncol(an)
        row <- numeric(pstar)
        for (k in seq_len(nk)) {
          a <- an[1L, k]
          bb <- an[2L, k]
          rat <- s_raw[i, a] * s_raw[j, bb] / s_raw[a, bb]
          row[vpos(i, a)] <- row[vpos(i, a)] +
            s_raw[j, bb] / s_raw[a, bb] / nk
          row[vpos(j, bb)] <- row[vpos(j, bb)] +
            s_raw[i, a] / s_raw[a, bb] / nk
          row[vpos(a, bb)] <- row[vpos(a, bb)] - rat / s_raw[a, bb] / nk
        }
        p_chain[vpos(i, j), ] <- row
      }
    }

    # Spearman theta chain (on the purged matrix)
    th_chain <- lav_cfa_mgm_jac_theta(s_mat = s_use, m_b = m_b)

    bi_list[[b]] <- list(
      ov_names = ov_names, lv_names = lv_names, nvar = nvar, nfac = nfac,
      pstar = pstar, m_b = m_b, marker_idx = marker_idx,
      s_use = s_use, p_chain = p_chain,
      theta0 = th_chain$theta0, dtheta = th_chain$dtheta,
      lam_pt_idx = lam_pt_idx, lam_ov = lam_ov, lam_fac = lam_fac,
      psi_rows_pt = psi_rows_pt, th_pt_idx = th_pt_idx,
      nobs = lavsamplestats@nobs[[b]]
    )
  }
  # ties among the factor (co)variances trigger the restricted PSI fit
  if (any(duplicated(psi_ids_all))) {
    return(NULL)
  }

  # resolved second stage (as in lav_cfa_guttman1952_internal)
  mgm_varcov <- ea[["mgm.varcov"]]
  if (is.null(mgm_varcov)) {
    mgm_varcov <- "default"
  }
  mgm_varcov <- toupper(mgm_varcov)
  if (mgm_varcov == "DEFAULT") {
    mgm_varcov <- if (res_cov_flag) "RLS" else "NONE"
  } else if (mgm_varcov == "NONE" && res_cov_flag) {
    mgm_varcov <- "RLS"
  }
  stage2 <- mgm_varcov != "NONE"
  # not covered: pinned nonzero ~~ cells make the stage-2 system affine
  # (the offset chain is not in lav_sem_js_deriv_varcov)
  if (stage2 && pin_flag) {
    return(NULL)
  }
  # not covered: tied residual variances without a second stage
  # (last-block-wins in the packed extraction)
  if (!stage2 && any(duplicated(th_ids_all))) {
    return(NULL)
  }

  # --------------------------------------------------------------------
  # stage 1
  # --------------------------------------------------------------------
  if (!pooled) {
    # classic per-block chain: per-block pd correction, T/U solve
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      nvar <- bi$nvar
      nfac <- bi$nfac
      pstar <- bi$pstar
      cc <- co$cov_vech[[b]]
      s_mat <- bi$s_use
      theta0 <- bi$theta0
      dtheta <- bi$dtheta
      xt <- bi$m_b

      rho_eff <- 1
      drho <- NULL
      root <- lav_cfa_mgm_jac_root(s_mat = s_mat, theta0 = theta0)
      if (is.null(root)) {
        return(NULL)
      }
      cutoff <- 1 + 1 / (bi$nobs - 1)
      if (root$rho < cutoff) {
        rho_eff <- root$rho - 1 / (bi$nobs - 1)
        drho <- lav_cfa_mgm_jac_droot(root = root, dtheta = dtheta)
      }

      w_mat <- s_mat - rho_eff * diag(theta0, nvar)
      t_mat <- crossprod(xt, w_mat %*% xt)
      t_inv <- try(solve(t_mat), silent = TRUE)
      if (inherits(t_inv, "try-error")) {
        return(NULL)
      }
      u_mat <- w_mat %*% xt %*% t_inv
      u_vec <- u_mat[cbind(bi$marker_idx, seq_len(nfac))]
      if (any(abs(u_vec) < .Machine$double.eps)) {
        return(NULL)
      }
      th0x <- theta0 * xt

      if (psi_mapping) {
        lam_full <- t(t(u_mat) / u_vec) * bi$m_b
        ti <- 1 / theta0
        ti_guard <- abs(theta0) < 0.01
        ti[ti_guard] <- 1
        b_map <- t(lam_full * ti) # L'K (nfac x nvar)
        a_map <- b_map %*% lam_full
        a_inv <- try(solve(a_map), silent = TRUE)
        if (inherits(a_inv, "try-error")) {
          return(NULL)
        }
        m_map <- a_inv %*% b_map
        wm <- w_mat %*% t(m_map) # W M' (nvar x nfac)
        mtm <- m_map %*% (theta0 * t(m_map)) # M Theta M'
      }

      th_ov <- match(lavpartable$lhs[bi$th_pt_idx], bi$ov_names)
      psi_f <- match(lavpartable$lhs[bi$psi_rows_pt], bi$lv_names)
      psi_g <- match(lavpartable$rhs[bi$psi_rows_pt], bi$lv_names)
      lam_free <- free[bi$lam_pt_idx]
      psi_free <- free[bi$psi_rows_pt]
      lam_ov <- bi$lam_ov
      lam_fac <- bi$lam_fac
      nlam <- length(lam_free)
      lam_val <- u_mat[cbind(lam_ov, lam_fac)] / u_vec[lam_fac]

      # theta rows (overwritten by stage 2 when it runs)
      if (length(bi$th_pt_idx) > 0L) {
        dx[free[bi$th_pt_idx], cc] <- dtheta[th_ov, , drop = FALSE]
        x_pred[free[bi$th_pt_idx]] <- theta0[th_ov]
      }
      x_pred[lam_free] <- lam_val
      if (psi_mapping) {
        psi_map_full <- m_map %*% wm
        x_pred[psi_free] <- psi_map_full[cbind(psi_f, psi_g)]
      } else {
        x_pred[psi_free] <- t_mat[cbind(psi_f, psi_g)] *
          u_vec[psi_f] * u_vec[psi_g]
      }

      # one pass over the local vech directions
      ii_of <- unlist(lapply(seq_len(nvar), function(j) seq(j, nvar)))
      jj_of <- unlist(lapply(seq_len(nvar), function(j) {
        rep(j, nvar - j + 1L)
      }))
      for (k in seq_len(pstar)) {
        i <- ii_of[k]
        j <- jj_of[k]
        dwx <- matrix(0, nvar, nfac)
        dwx[i, ] <- xt[j, ]
        if (i != j) {
          dwx[j, ] <- dwx[j, ] + xt[i, ]
        }
        if (!is.null(drho)) {
          dwx <- dwx - drho[k] * th0x
        }
        dthk <- dtheta[, k]
        nz <- which(dthk != 0)
        if (length(nz) > 0L) {
          dwx[nz, ] <- dwx[nz, ] -
            rho_eff * dthk[nz] * xt[nz, , drop = FALSE]
        }
        dt_mat <- crossprod(xt, dwx)
        du_mat <- (dwx - u_mat %*% dt_mat) %*% t_inv
        du_vec <- du_mat[cbind(bi$marker_idx, seq_len(nfac))]
        if (nlam > 0L) {
          dx[lam_free, cc[k]] <-
            (du_mat[cbind(lam_ov, lam_fac)] -
             lam_val * du_vec[lam_fac]) / u_vec[lam_fac]
        }
        if (psi_mapping) {
          dlam <- (du_mat - t(t(u_mat) * (du_vec / u_vec))) *
            bi$m_b / rep(u_vec, each = nvar)
          dti <- -dthk * ti * ti
          dti[ti_guard] <- 0
          db_map <- t(dlam * ti) + t(lam_full * dti)
          da_map <- db_map %*% lam_full + b_map %*% dlam
          dm_map <- a_inv %*% (db_map - da_map %*% m_map)
          if (i != j) {
            mdwm <- tcrossprod(m_map[, i], m_map[, j])
            mdwm <- mdwm + t(mdwm)
          } else {
            mdwm <- tcrossprod(m_map[, i])
          }
          if (!is.null(drho)) {
            mdwm <- mdwm - drho[k] * mtm
          }
          if (length(nz) > 0L) {
            mdwm <- mdwm - rho_eff *
              (m_map[, nz, drop = FALSE] %*%
               (dthk[nz] * t(m_map[, nz, drop = FALSE])))
          }
          dpsi_mat <- dm_map %*% wm + mdwm + t(dm_map %*% wm)
          dx[psi_free, cc[k]] <- dpsi_mat[cbind(psi_f, psi_g)]
        } else {
          du_f <- du_vec[psi_f]
          du_g <- du_vec[psi_g]
          dx[psi_free, cc[k]] <-
            dt_mat[cbind(psi_f, psi_g)] * u_vec[psi_f] * u_vec[psi_g] +
            t_mat[cbind(psi_f, psi_g)] *
              (du_f * u_vec[psi_g] + u_vec[psi_f] * du_g)
        }
      }
    }
  } else {
    # ---- pooled branch (loading ties / quadprog) -----------------------
    # pd correction: JOINT shrink over the blocks (nblocks > 1), or the
    # per-block root (single block); the joint-shrink derivative has
    # cross-block support
    shrink <- rep(1, nblocks) # effective multiplier of Theta, per block
    dshrink <- NULL # 1 x m_total (shared over the blocks), or NULL
    if (nblocks > 1L) {
      roots <- numeric(nblocks)
      root_l <- vector("list", nblocks)
      for (b in seq_len(nblocks)) {
        root_l[[b]] <- lav_cfa_mgm_jac_root(
          s_mat = bi_list[[b]]$s_use, theta0 = bi_list[[b]]$theta0
        )
        if (is.null(root_l[[b]])) {
          return(NULL)
        }
        roots[b] <- root_l[[b]]$rho
      }
      bstar <- which.min(roots)
      # the argmin block must be well separated
      if (sort(roots)[2L] - roots[bstar] < 1e-8 * max(1, abs(roots[bstar]))) {
        return(NULL)
      }
      nobs_min <- min(vapply(bi_list, "[[", numeric(1L), "nobs"))
      cutoff <- 1 + 1 / (nobs_min - 1)
      if (roots[bstar] < cutoff) {
        sh <- roots[bstar] - 1 / (nobs_min - 1)
        shrink <- rep(sh, nblocks)
        dshrink <- numeric(m_total)
        dshrink[co$cov_vech[[bstar]]] <- lav_cfa_mgm_jac_droot(
          root = root_l[[bstar]], dtheta = bi_list[[bstar]]$dtheta
        )
      }
    } else {
      root <- lav_cfa_mgm_jac_root(
        s_mat = bi_list[[1L]]$s_use, theta0 = bi_list[[1L]]$theta0
      )
      if (is.null(root)) {
        return(NULL)
      }
      cutoff <- 1 + 1 / (bi_list[[1L]]$nobs - 1)
      if (root$rho < cutoff) {
        shrink[1L] <- root$rho - 1 / (bi_list[[1L]]$nobs - 1)
        dshrink <- numeric(m_total)
        dshrink[co$cov_vech[[1L]]] <- lav_cfa_mgm_jac_droot(
          root = root, dtheta = bi_list[[1L]]$dtheta
        )
      }
    }

    # per block: T, WX, and full-width derivative rows for the needed
    # primitives (pattern cells of WX, and the T cells)
    cells <- NULL
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      nvar <- bi$nvar
      nfac <- bi$nfac
      cc <- co$cov_vech[[b]]
      xt <- bi$m_b
      w_mat <- bi$s_use - shrink[b] * diag(bi$theta0, nvar)
      wx <- w_mat %*% xt
      t_mat <- crossprod(xt, wx)
      # the historic cov2cor(force_pd(phi)) step must be an identity
      dvec <- 1 / sqrt(diag(t_mat))
      phi <- t(t_mat * dvec) * dvec
      evp <- eigen(phi, symmetric = TRUE, only.values = TRUE)$values
      if (min(evp) / max(evp) < 2e-4) {
        return(NULL) # force_pd about to clip
      }
      # dtheta placed at the global columns
      dth_g <- matrix(0, nrow = nvar, ncol = m_total)
      dth_g[, cc] <- bi$dtheta
      vpos <- function(i, j) {
        if (i < j) {
          tmp <- i
          i <- j
          j <- tmp
        }
        (j - 1L) * nvar - ((j - 1L) * j) %/% 2L + i
      }
      # dT rows (nfac x nfac, symmetric)
      dt_rows <- vector("list", nfac * nfac)
      for (f in seq_len(nfac)) {
        f_set <- which(xt[, f] == 1L)
        for (g in seq_len(f)) {
          g_set <- which(xt[, g] == 1L)
          row <- numeric(m_total)
          for (p in f_set) {
            for (j in g_set) {
              row[cc[vpos(p, j)]] <- row[cc[vpos(p, j)]] + 1
            }
          }
          if (f == g) {
            row <- row - shrink[b] * colSums(dth_g[f_set, , drop = FALSE])
            if (!is.null(dshrink)) {
              row <- row - sum(bi$theta0[f_set]) * dshrink
            }
          }
          dt_rows[[(f - 1L) * nfac + g]] <- row
          dt_rows[[(g - 1L) * nfac + f]] <- row
        }
      }
      # d(WX) rows for the pattern cells (markers included)
      pat <- which(xt == 1L)
      pat_ov <- ((pat - 1L) %% nvar) + 1L
      pat_fac <- ((pat - 1L) %/% nvar) + 1L
      dwx_rows <- matrix(0, nrow = length(pat), ncol = m_total)
      for (e in seq_along(pat)) {
        p <- pat_ov[e]
        f <- pat_fac[e]
        f_set <- which(xt[, f] == 1L)
        row <- numeric(m_total)
        for (j in f_set) {
          row[cc[vpos(p, j)]] <- row[cc[vpos(p, j)]] + 1
        }
        # p loads on f, so the diagonal cell (p, p) is always in the sum
        row <- row - shrink[b] * dth_g[p, ]
        if (!is.null(dshrink)) {
          row <- row - bi$theta0[p] * dshrink
        }
        dwx_rows[e, ] <- row
      }
      # normalized structure coefficients lstar = (WX)_pf / sqrt(T_ff)
      lstar_val <- numeric(length(pat))
      dlstar <- matrix(0, nrow = length(pat), ncol = m_total)
      for (e in seq_along(pat)) {
        f <- pat_fac[e]
        tf <- t_mat[f, f]
        lstar_val[e] <- wx[pat_ov[e], f] / sqrt(tf)
        dlstar[e, ] <- dwx_rows[e, ] / sqrt(tf) -
          0.5 * lstar_val[e] * dt_rows[[(f - 1L) * nfac + f]] / tf
      }
      cell_of <- matrix(0L, nrow = nvar, ncol = nfac)
      cell_of[cbind(pat_ov, pat_fac)] <- seq_along(pat)

      bi_list[[b]]$w_mat <- w_mat
      bi_list[[b]]$wx <- wx
      bi_list[[b]]$t_mat <- t_mat
      bi_list[[b]]$dt_rows <- dt_rows
      bi_list[[b]]$lstar_val <- lstar_val
      bi_list[[b]]$dlstar <- dlstar
      bi_list[[b]]$cell_of <- cell_of

      # theta rows
      if (length(bi$th_pt_idx) > 0L) {
        th_ov <- match(lavpartable$lhs[bi$th_pt_idx], bi$ov_names)
        dx[free[bi$th_pt_idx], cc] <- bi$dtheta[th_ov, , drop = FALSE]
        x_pred[free[bi$th_pt_idx]] <- bi$theta0[th_ov]
      }

      if (length(bi$lam_pt_idx) > 0L) {
        cells <- rbind(cells, data.frame(
          block = b, ov = bi$lam_ov, fac = bi$lam_fac,
          facname = lavpartable$lhs[bi$lam_pt_idx],
          free = free[bi$lam_pt_idx]
        ))
      }
    }

    # tie classes and connected components (as in the internal function)
    if (is.null(cells)) {
      cells <- data.frame(
        block = integer(0L), ov = integer(0L), fac = integer(0L),
        facname = character(0L), free = integer(0L)
      )
    }
    tie_classes <- split(seq_len(nrow(cells)), cells$free)
    tie_classes <- tie_classes[vapply(tie_classes, length, 1L) > 1L]
    fac_components <- function(fn) {
      comps <- as.list(seq_len(nblocks))
      for (class_idx in tie_classes) {
        if (cells$facname[class_idx[1L]] != fn) {
          next
        }
        bset <- unique(cells$block[class_idx])
        hit <- which(vapply(comps, function(x) any(x %in% bset), TRUE))
        if (length(hit) > 1L) {
          comps[[hit[1L]]] <- sort(unique(unlist(comps[hit])))
          comps[hit[-1L]] <- NULL
        }
      }
      comps
    }

    # (pooled) normalized marker loadings, value + derivative row
    mh_val <- vector("list", nblocks)
    mh_row <- vector("list", nblocks)
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      mh_val[[b]] <- numeric(bi$nfac)
      mh_row[[b]] <- vector("list", bi$nfac)
      for (f in seq_len(bi$nfac)) {
        e <- bi$cell_of[bi$marker_idx[f], f]
        mh_val[[b]][f] <- bi$lstar_val[e]
        mh_row[[b]][[f]] <- bi$dlstar[e, ]
      }
      names(mh_val[[b]]) <- bi$lv_names
      names(mh_row[[b]]) <- bi$lv_names
    }
    all_fn <- unique(cells$facname)
    for (fn in all_fn) {
      comps <- fac_components(fn)
      for (ccomp in comps) {
        bset <- ccomp[vapply(ccomp, function(b) {
          fn %in% bi_list[[b]]$lv_names
        }, logical(1L))]
        if (length(bset) < 2L) {
          next
        }
        vals <- vapply(bset, function(b) mh_val[[b]][[fn]], numeric(1L))
        rows <- lapply(bset, function(b) mh_row[[b]][[fn]])
        pooled_val <- mean(vals)
        pooled_row <- Reduce(`+`, rows) / length(bset)
        for (b in bset) {
          mh_val[[b]][[fn]] <- pooled_val
          mh_row[[b]][[fn]] <- pooled_row
        }
      }
    }

    # cell values (tied classes averaged) and the loading rows
    val <- numeric(nrow(cells))
    vrow <- matrix(0, nrow = nrow(cells), ncol = m_total)
    for (r in seq_len(nrow(cells))) {
      bi <- bi_list[[cells$block[r]]]
      e <- bi$cell_of[cells$ov[r], cells$fac[r]]
      val[r] <- bi$lstar_val[e]
      vrow[r, ] <- bi$dlstar[e, ]
    }
    for (class_idx in tie_classes) {
      val[class_idx] <- mean(val[class_idx])
      avg <- colMeans(vrow[class_idx, , drop = FALSE])
      vrow[class_idx, ] <- rep(avg, each = length(class_idx))
    }
    for (r in seq_len(nrow(cells))) {
      mh <- mh_val[[cells$block[r]]][[cells$facname[r]]]
      dmh <- mh_row[[cells$block[r]]][[cells$facname[r]]]
      lam_val_r <- val[r] / mh
      dx[cells$free[r], ] <- (vrow[r, ] - lam_val_r * dmh) / mh
      x_pred[cells$free[r]] <- lam_val_r
    }

    # psi rows: psi_fg = T_fg (WX)_{mf, f} (WX)_{mg, g} / (T_ff T_gg)
    # (= phi_fg mhat_own_f mhat_own_g, with the OWN -- unpooled --
    # normalized marker loadings)
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      nfac <- bi$nfac
      psi_f <- match(lavpartable$lhs[bi$psi_rows_pt], bi$lv_names)
      psi_g <- match(lavpartable$rhs[bi$psi_rows_pt], bi$lv_names)
      for (r in seq_along(bi$psi_rows_pt)) {
        f <- psi_f[r]
        g <- psi_g[r]
        e_f <- bi$cell_of[bi$marker_idx[f], f]
        e_g <- bi$cell_of[bi$marker_idx[g], g]
        wxf <- bi$wx[bi$marker_idx[f], f]
        wxg <- bi$wx[bi$marker_idx[g], g]
        tff <- bi$t_mat[f, f]
        tgg <- bi$t_mat[g, g]
        tfg <- bi$t_mat[f, g]
        dwxf <- bi$dlstar[e_f, ] * sqrt(tff) +
          0.5 * (wxf / tff) * bi$dt_rows[[(f - 1L) * nfac + f]]
        dwxg <- bi$dlstar[e_g, ] * sqrt(tgg) +
          0.5 * (wxg / tgg) * bi$dt_rows[[(g - 1L) * nfac + g]]
        aa <- wxf * wxg / (tff * tgg)
        psi_val <- tfg * aa
        drow <- bi$dt_rows[[(f - 1L) * nfac + g]] * aa +
          tfg * ((dwxf * wxg + wxf * dwxg) / (tff * tgg) -
                 aa * (bi$dt_rows[[(f - 1L) * nfac + f]] / tff +
                       bi$dt_rows[[(g - 1L) * nfac + g]] / tgg))
        dx[free[bi$psi_rows_pt[r]], ] <- drow
        x_pred[free[bi$psi_rows_pt[r]]] <- psi_val
      }
    }
  }

  # --------------------------------------------------------------------
  # compose the purge chains (all stage-1 derivatives were taken with
  # respect to the PURGED matrices)
  # --------------------------------------------------------------------
  for (b in seq_len(nblocks)) {
    if (!is.null(bi_list[[b]]$p_chain)) {
      cc <- co$cov_vech[[b]]
      dx[, cc] <- dx[, cc, drop = FALSE] %*% bi_list[[b]]$p_chain
    }
  }

  # --------------------------------------------------------------------
  # stage 2 (shared with the JS estimator) and the mean structure
  # --------------------------------------------------------------------
  x0 <- as.numeric(lav_model_get_parameters(lavmodel))
  lavmodel0 <- lav_model_set_parameters(lavmodel = lavmodel, x = x0)
  lavh1_use <- list(implied = implied0)
  if (stage2) {
    free_directed_idx <- unique(free[lam_all])
    delta_list0 <- lav_sem_miiv_delta(lavmodel0)
    for (b in seq_len(nblocks)) {
      fu_b <- unique(free[op == "~~" & free > 0L & pt_block == b])
      if (length(fu_b) == 0L) {
        next
      }
      # this block's own solution (shared undirected parameters are
      # overwritten by later blocks in the point estimation)
      theta2_b <- as.numeric(lav_sem_miiv_varcov_block(
        b = b, fu = fu_b, delta_b = delta_list0[[b]],
        implied = implied0, lavmodel = lavmodel, x = x0,
        iv_varcov_method = mgm_varcov, return_h = FALSE
      ))
      dx[fu_b, ] <- lav_sem_js_deriv_varcov(
        b = b, fu = fu_b, delta_b = delta_list0[[b]],
        lavmodel = lavmodel, lavmodel0 = lavmodel0,
        lavh1 = lavh1_use, x0 = x0,
        js_varcov_method = mgm_varcov,
        free_directed_idx = free_directed_idx,
        dx = dx, co_b = co$block[[b]], theta2_use = theta2_b
      )
      x_pred[fu_b] <- theta2_b
    }
  }

  if (lavmodel@meanstructure) {
    free_mean_idx <- unique(free[op == "~1" & free > 0L])
    if (length(free_mean_idx) > 0L) {
      dx[free_mean_idx, ] <- lav_sem_js_deriv_mean_wls(
        lavmodel = lavmodel, lavmodel0 = lavmodel0,
        lavsamplestats = lavsamplestats, lavh1 = lavh1_use,
        x0 = x0, free_mean_idx = free_mean_idx, dx = dx, co = co
      )
      x_pred[free_mean_idx] <- as.numeric(lav_sem_miiv_mean_wls(
        lavmodel = lavmodel, lavsamplestats = lavsamplestats, x = x0,
        free_mean_idx = free_mean_idx, sample_mean = implied0$mean
      ))
    }
  }

  # verify: every free parameter must be covered, and the predicted
  # estimates must match the actual ones (they differ when a bound clip
  # or a branch the checks above missed fired at estimation time)
  if (anyNA(x_pred)) {
    return(NULL)
  }
  if (max(abs(x_pred - x0)) > 1e-6 * (1 + max(abs(x0)))) {
    return(NULL)
  }

  dx
}
lav_noniter_vcov <- function(lavmodel = NULL, lavsamplestats = NULL,
                             lavoptions = NULL, lavpartable = NULL,
                             lavimplied = NULL, lavh1 = NULL,
                             lavdata = NULL) {
  # the engine handles continuous, unconditional moments only
  if (lavmodel@categorical || lavmodel@conditional.x ||
      lavmodel@correlation || lavmodel@group.w.free) {
    lav_msg_stop(gettextf(
      "delta-method standard errors for estimator %s require continuous,
       unconditional sample moments (no categorical, conditional.x,
       correlation or group.w.free structures).", lavoptions$estimator))
  }

  # Gamma flavor; NOTE: read with [[ ]] -- $ would partially match
  gamma_flavor <- lavoptions$estimator.args[["mgm.gamma"]]
  if (is.null(gamma_flavor)) {
    gamma_flavor <- "nt"
  }
  gamma_flavor <- tolower(gamma_flavor)
  if (gamma_flavor == "adf" && lavdata@data.type != "full") {
    lav_msg_warn(gettext(
      "the \"adf\" moment covariance requires complete raw data; using
       the normal-theory moment covariance instead."))
    gamma_flavor <- "nt"
  }

  # the sample moments (per block: mean if meanstructure, then vech(S)),
  # in the same layout that lav_vec_to_implied() unpacks
  implied0 <- lavh1$implied
  if (is.null(implied0$cov)) {
    implied0 <- list(cov = lavsamplestats@cov, mean = lavsamplestats@mean)
  }
  vec0 <- numeric(0L)
  for (b in seq_len(lavmodel@nblocks)) {
    if (lavmodel@meanstructure) {
      vec0 <- c(vec0, implied0$mean[[b]])
    }
    vec0 <- c(vec0, lav_mat_vech(implied0$cov[[b]]))
  }

  # Jacobian of the complete estimation map over the stacked moments:
  # analytic for the classic MGM branch (mgm.jacobian = "analytic", the
  # default; lav_cfa_mgm_jacobian returns NULL for the branches it does
  # not cover), numerical otherwise
  jac <- NULL
  mgm_jacobian <- lavoptions$estimator.args[["mgm.jacobian"]]
  if (is.null(mgm_jacobian)) {
    mgm_jacobian <- "analytic"
  }
  if (lavoptions$estimator == "MGM" &&
      tolower(mgm_jacobian) == "analytic") {
    jac <- try(lav_cfa_mgm_jacobian(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats,
      lavoptions = lavoptions, lavpartable = lavpartable,
      implied0 = implied0
    ), silent = TRUE)
    if (inherits(jac, "try-error")) {
      jac <- NULL
    }
  }
  # numerical fallback; r = 2 Richardson extrapolation is accurate to
  # ~1e-6 relative and halves the number of map evaluations (the map
  # itself is smooth); warnings raised on perturbed evaluations (purge
  # admissibility, RLS iteration counts) merely repeat what the point
  # estimation already reported, so they are suppressed here
  if (is.null(jac)) {
    jac <- suppressWarnings(numDeriv::jacobian(
      func = function(v) {
        lav_noniter_estimate_from_vec(
          vec = v, lavmodel = lavmodel, lavsamplestats = lavsamplestats,
          lavpartable = lavpartable, lavdata = lavdata,
          lavh1 = lavh1, lavoptions = lavoptions
        )
      },
      x = vec0, method.args = list(r = 2L)
    ))
  }

  # vcov = sum over blocks of J_b Gamma_b J_b' / nobs_b; for the (default)
  # normal-theory Gamma the sandwich is computed directly from the sample
  # covariance matrix (lav_mat_k_gammant_kt with K = J_b), without ever
  # forming the pstar x pstar Gamma matrix
  npar <- nrow(jac)
  vcov <- matrix(0, nrow = npar, ncol = npar)
  offset <- 0L
  for (b in seq_len(lavmodel@nblocks)) {
    nvar_b <- lavmodel@nvar[b]
    nd_b <- if (lavmodel@meanstructure) {
      nvar_b + nvar_b * (nvar_b + 1L) / 2L
    } else {
      nvar_b * (nvar_b + 1L) / 2L
    }
    jac_b <- jac[, offset + seq_len(nd_b), drop = FALSE]
    offset <- offset + nd_b
    if (gamma_flavor == "adf") {
      gg <- lav_samp_gamma(
        m_y = lavdata@X[[b]],
        meanstructure = lavmodel@meanstructure
      )
      vcov <- vcov + (jac_b %*% gg %*% t(jac_b)) / lavsamplestats@nobs[[b]]
    } else {
      vcov <- vcov + lav_mat_k_gammant_kt(
        m_k = jac_b, s = implied0$cov[[b]],
        meanstructure = lavmodel@meanstructure,
        x_idx = integer(0L)
      ) / lavsamplestats@nobs[[b]]
    }
  }
  vcov <- (vcov + t(vcov)) / 2

  # the rest of lavaan expects the vcov of a ceq.simple.only model in the
  # 'unco' space (one row/column per non-collapsed parameter); the vcov
  # above is built in the compact (nx.free) space, so expand it via the
  # ceq.simple.K mapping (as for the IV and JS estimators)
  if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L) {
    vcov <- lavmodel@ceq.simple.K %*% vcov %*% t(lavmodel@ceq.simple.K)
  }

  vcov
}
