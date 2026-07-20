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
#  - stage 1, classic (simple structure, saturated PSI):
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
#  - stage 1, cross-loadings (rotation-to-pattern; no pd correction in
#    this branch by design): the bordered KKT rotation solve, the
#    Bartlett mapping over the clean rows, the masked-composite solves
#    for the cross rows, and the model-identity thetas
#  - factors with fewer than 3 indicators, and factors touching a
#    cross-loading: the masked/external (Kano) tetrad thetas, via the
#    embedded derivative of lav_sem_js_theta_spearman_masked()
#  - the metric rescale (std.lv, fixed loading c != 1, tied factors):
#    a closed-form post-composition lambda -> lambda s,
#    psi -> psi / (s s'), with s = c or s = sqrt(psi_ff / v), pooled
#    per connected component
#  - the restricted PSI fit (fixed factor covariances, psi ties, loose
#    variance pins): implicit function theorem at the ML optimum,
#    d x = -H^{-1} (dg/dV) dV, with V the (rescaled) stage-1 psi
#  - residual covariances: the purge s0_ij = mean_(a,b) s_ia s_jb /
#    s_ab is composed with the stage-1 derivatives (which are evaluated
#    on the purged matrix)
#  - stage 2 (mgm.varcov = ULS/GLS/RLS/2RLS): reuses the JS stage-2
#    derivative (lav_sem_js_deriv_varcov), with the MGM loading
#    derivatives as dx
#  - mean structure: reuses the JS joint-GLS mean-solve derivative
#    (free/tied intercepts, free latent means)
#
# Remaining fallbacks: cross-loadings combined with equality ties,
# psi.mapping in the pooled branch, variance/covariance cells fixed to
# a NONZERO value combined with the second stage (the affine-offset
# chain), tied residual variances without a second stage, unpurgeable
# residual covariances, and (near-)non-smooth points (multiple pencil
# roots, a clipping force-pd step). Regime switches (Spearman clips,
# theta bounds, pd correction active or not, the zero-theta guard, the
# diagonal floor of the cross branch) are frozen at the point estimate.
# A final safety net compares the predicted estimates against the
# actual ones; any mismatch also triggers the numerical fallback.

# Spearman theta chain for one block: the (clipped) initial residual
# variances and their derivative with respect to the local vech(S).
# Factors with >= 3 indicators and no cross-loading member use the
# plain per-factor tetrads; the others use the masked/external tetrads
# (lav_sem_js_theta_spearman_masked, which carries its own derivative)
lav_cfa_mgm_jac_theta <- function(s_mat = NULL, m_b = NULL,
                                  cross_idx = integer(0L)) {
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
  both_cross <- matrix(FALSE, nvar, nvar)
  if (length(cross_idx) > 1L) {
    both_cross[cross_idx, cross_idx] <- TRUE
  }
  theta0 <- numeric(nvar)
  dtheta <- matrix(0, nvar, pstar)
  for (f in seq_len(nfac)) {
    f_idx <- which(m_b[, f] == 1L)
    if (length(f_idx) >= 3L && !any(f_idx %in% cross_idx)) {
      # plain Spearman on this factor
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
          drow[vpos(p, a)] <- drow[vpos(p, a)] +
            s_mat[p, bb] / s_mat[a, bb]
          drow[vpos(p, bb)] <- drow[vpos(p, bb)] +
            s_mat[p, a] / s_mat[a, bb]
          drow[vpos(a, bb)] <- drow[vpos(a, bb)] - rat / s_mat[a, bb]
        }
        h2 <- asum / (nk * s_mat[p, p])
        if (h2 < -0.05 || h2 > 1.20) {
          h2c <- if (h2 < -0.05) -0.05 else 1.20
          theta0[p] <- (1 - h2c) * s_mat[p, p]
          dtheta[p, vpos(p, p)] <- 1 - h2c
        } else {
          theta0[p] <- s_mat[p, p] - asum / nk
          dtheta[p, ] <- -drow / nk
          dtheta[p, vpos(p, p)] <- dtheta[p, vpos(p, p)] + 1
        }
      }
    } else {
      # masked/external tetrads (the JS helper carries the derivative)
      out <- lav_sem_js_theta_spearman_masked(
        s_full = s_mat, idx = f_idx, contam = both_cross,
        zpool = setdiff(seq_len(nvar), f_idx), bounds = "wide",
        deriv = TRUE
      )
      if (anyNA(out)) {
        return(NULL) # the estimation would have stopped
      }
      theta0[f_idx] <- as.numeric(out)
      dtheta[f_idx, ] <- attr(out, "dtheta")
    }
    # worker bounds: clip to [0, s_pp]
    for (p in f_idx) {
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
  pin_flag <- FALSE # any ~~ cell fixed to a NONZERO value
  th_ids_all <- integer(0L)
  any_cross <- FALSE
  for (b in seq_len(nblocks)) {
    ov_names <- lavpta$vnames$ov[[b]]
    lv_names <- lavpta$vnames$lv.regular[[b]]
    nvar <- length(ov_names)
    nfac <- length(lv_names)
    pstar <- nvar * (nvar + 1L) / 2L

    # scale identification per factor (as in the internal function):
    # a fixed nonzero loading ("marker", value c), a fixed positive
    # factor variance ("variance", std.lv), or ties ("tied")
    scale_type <- rep("tied", nfac)
    marker_value <- rep(1, nfac)
    scale_target <- rep(NA_real_, nfac)
    marker_idx <- rep(0L, nfac)
    lam_pt_idx <- which(op == "=~" & pt_block == b & free > 0L)
    lam_ov <- match(lavpartable$rhs[lam_pt_idx], ov_names)
    lam_fac <- match(lavpartable$lhs[lam_pt_idx], lv_names)
    if (anyNA(lam_ov) || anyNA(lam_fac)) {
      return(NULL)
    }
    m_b <- matrix(0L, nvar, nfac)
    m_b[(lam_fac - 1L) * nvar + lam_ov] <- 1L
    for (f in seq_len(nfac)) {
      fx <- which(op == "=~" & pt_block == b &
        lavpartable$lhs == lv_names[f] & free == 0L &
        !is.na(lavpartable$ustart) & lavpartable$ustart != 0)
      if (length(fx) > 1L) {
        return(NULL) # the estimation stops on this
      }
      if (length(fx) == 1L) {
        scale_type[f] <- "marker"
        marker_value[f] <- lavpartable$ustart[fx]
        marker_idx[f] <- match(lavpartable$rhs[fx], ov_names)
        m_b[(f - 1L) * nvar + marker_idx[f]] <- 1L
      } else {
        vrow <- which(op == "~~" & pt_block == b &
          lavpartable$lhs == lv_names[f] &
          lavpartable$rhs == lv_names[f] & free == 0L)
        if (length(vrow) > 0L) {
          scale_type[f] <- "variance"
          scale_target[f] <- lavpartable$ustart[vrow[1L]]
        }
      }
    }
    cross_idx <- which(rowSums(m_b) > 1L)
    if (length(cross_idx) > 0L) {
      any_cross <- TRUE
    }
    # provisional scaling indicators (cross-free)
    for (f in which(marker_idx == 0L)) {
      cand <- which(m_b[, f] == 1L)
      cand <- cand[!cand %in% cross_idx]
      if (length(cand) == 0L) {
        return(NULL) # the estimation stops on this
      }
      marker_idx[f] <- cand[1L]
    }
    if (any(marker_idx %in% cross_idx)) {
      return(NULL) # the estimation stops on this
    }

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

    # psi structure (restrictions are handled by the fit chain below)
    psi_rows_pt <- which(op == "~~" & pt_block == b &
      lavpartable$lhs %in% lv_names & lavpartable$rhs %in% lv_names)
    psi_id <- matrix(0L, nfac, nfac)
    psi_pin <- matrix(0, nfac, nfac)
    if (length(psi_rows_pt) > 0L) {
      ii <- match(lavpartable$lhs[psi_rows_pt], lv_names)
      jj <- match(lavpartable$rhs[psi_rows_pt], lv_names)
      for (r in seq_along(psi_rows_pt)) {
        fr <- free[psi_rows_pt[r]]
        if (fr > 0L) {
          psi_id[ii[r], jj[r]] <- psi_id[jj[r], ii[r]] <- fr
        } else {
          val <- lavpartable$ustart[psi_rows_pt[r]]
          if (is.na(val)) {
            val <- 0
          }
          psi_pin[ii[r], jj[r]] <- psi_pin[jj[r], ii[r]] <- val
        }
      }
    }
    if (any(psi_pin != 0)) {
      # nonzero psi pins beyond the (absorbable) variance targets
      offd <- psi_pin
      diag(offd) <- 0
      if (any(offd != 0)) {
        pin_flag <- TRUE # fixed nonzero factor covariance
      }
      dg <- diag(psi_pin)
      if (any(dg != 0 & scale_type != "variance")) {
        pin_flag <- TRUE # fixed variance on a marker-scaled factor
      }
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
      if (length(cross_idx) > 0L) {
        fac_of[cross_idx] <- NA_integer_
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
        return(NULL) # unpurgeable pair
      }
      anchors <- attr(s_use, "anchors")
      attr(s_use, "failed") <- NULL
      attr(s_use, "anchors") <- NULL
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
    th_chain <- lav_cfa_mgm_jac_theta(
      s_mat = s_use, m_b = m_b, cross_idx = cross_idx
    )
    if (is.null(th_chain)) {
      return(NULL)
    }

    bi_list[[b]] <- list(
      ov_names = ov_names, lv_names = lv_names, nvar = nvar, nfac = nfac,
      pstar = pstar, m_b = m_b, marker_idx = marker_idx,
      scale_type = scale_type, marker_value = marker_value,
      scale_target = scale_target, cross_idx = cross_idx,
      s_use = s_use, p_chain = p_chain,
      theta0 = th_chain$theta0, dtheta = th_chain$dtheta,
      lam_pt_idx = lam_pt_idx, lam_ov = lam_ov, lam_fac = lam_fac,
      psi_rows_pt = psi_rows_pt, psi_id = psi_id, psi_pin = psi_pin,
      th_pt_idx = th_pt_idx,
      nobs = lavsamplestats@nobs[[b]]
    )
  }
  pooled <- (has_ties || quadprog) && !any_cross
  # not covered (yet): cross-loadings combined with ties, psi.mapping in
  # the pooled branch or over an unzeroed lambda
  if (any_cross && has_ties) {
    return(NULL)
  }
  if (psi_mapping && (pooled || !zae)) {
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
  # not covered: pinned nonzero ~~ cells (incl. std.lv variance pins)
  # make the stage-2 system affine (offset chain not in the JS deriv)
  if (stage2 && (pin_flag ||
      any(vapply(bi_list, function(bi) any(diag(bi$psi_pin) != 0),
                 logical(1L))))) {
    return(NULL)
  }
  # not covered: tied residual variances without a second stage
  # (last-block-wins in the packed extraction)
  if (!stage2 && any(duplicated(th_ids_all))) {
    return(NULL)
  }
  if (pin_flag && !stage2) {
    # fixed nonzero theta diagonals without stage 2 are fine (they do
    # not enter the stage-1 chain); fixed nonzero rescov pairs are
    # purge-subtracted with zero derivative -- also fine
    pin_flag <- FALSE
  }

  # --------------------------------------------------------------------
  # stage 1: fills, per block, the full value/derivative structures
  #   lam  (nvar x nfac), dlam ((f-1)*nvar+p rows x m_total)
  #   psi  (nfac x nfac), dpsi (lower-tri rows x m_total)
  #   theta, dth (nvar x m_total)
  # in the provisional marker = 1 metric
  # --------------------------------------------------------------------
  pidx <- function(f, g, nfac) {
    if (f < g) {
      tmp <- f
      f <- g
      g <- tmp
    }
    (g - 1L) * nfac - ((g - 1L) * g) %/% 2L + f
  }

  if (any_cross) {
    # ---- cross-loadings: rotation-to-pattern (no pd correction) -------
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      nvar <- bi$nvar
      nfac <- bi$nfac
      pstar <- bi$pstar
      cc <- co$cov_vech[[b]]
      xt <- bi$m_b
      theta0 <- bi$theta0
      dtheta <- bi$dtheta
      s_mat <- bi$s_use
      cross_idx <- bi$cross_idx
      clean <- setdiff(seq_len(nvar), cross_idx)

      # W = S - diag(theta), with the (rare) diagonal floor of the
      # !force_pd path: diag(W) >= 0.001 * diag(S)
      w_mat <- s_mat - diag(theta0, nvar)
      floor_idx <- which(diag(w_mat) < 0.001 * diag(s_mat))
      if (length(floor_idx) > 0L) {
        diag(w_mat)[floor_idx] <- 0.001 * diag(s_mat)[floor_idx]
      }
      t_mat <- crossprod(xt, w_mat %*% xt)
      t_inv <- try(solve(t_mat), silent = TRUE)
      if (inherits(t_inv, "try-error")) {
        return(NULL)
      }
      u_mat <- w_mat %*% xt %*% t_inv
      dsq <- sqrt(diag(t_mat)) # D^-1 of the worker
      mm <- t(t(u_mat) * dsq) # the unconstrained solve t(solve(phi, ys_cor))

      # rotation: per factor, the bordered KKT solve on the clean zero
      # rows of mm
      hmat <- matrix(0, nfac, nfac)
      kkt_inv <- vector("list", nfac)
      zrows_l <- vector("list", nfac)
      for (f in seq_len(nfac)) {
        zrows <- intersect(which(xt[, f] == 0L), clean)
        if (length(zrows) < nfac - 1L) {
          return(NULL) # the estimation stops on this
        }
        zf <- mm[zrows, , drop = FALSE]
        mrow <- mm[bi$marker_idx[f], ]
        kkt <- rbind(cbind(crossprod(zf), mrow), c(mrow, 0))
        kkt_i <- try(solve(kkt), silent = TRUE)
        if (inherits(kkt_i, "try-error")) {
          return(NULL)
        }
        sol <- kkt_i %*% c(numeric(nfac), 1)
        hmat[, f] <- sol[seq_len(nfac)]
        kkt_inv[[f]] <- kkt_i
        zrows_l[[f]] <- zrows
      }
      lam <- (mm %*% hmat) * bi$m_b

      # PSI: Bartlett mapping over the clean rows
      lam_c <- lam[clean, , drop = FALSE]
      th_c <- theta0[clean]
      ti_c <- 1 / th_c
      ti_guard <- abs(th_c) < 0.01
      ti_c[ti_guard] <- 1
      b_map <- t(lam_c * ti_c) # L'K (nfac x nclean)
      a_map <- b_map %*% lam_c
      a_inv <- try(solve(a_map), silent = TRUE)
      if (inherits(a_inv, "try-error")) {
        return(NULL)
      }
      m_map <- a_inv %*% b_map # nfac x nclean
      w_c <- w_mat[clean, clean, drop = FALSE]
      wm <- w_c %*% t(m_map) # nclean x nfac
      psi <- m_map %*% wm

      # cross rows: masked-composite solves
      x_masked <- xt
      if (length(cross_idx) > 0L) {
        x_masked[cross_idx, ] <- 0
      }
      am_x <- crossprod(x_masked, lam) # X_m' Lambda (nfac x nfac)
      bmat <- psi %*% t(am_x)
      cross_info <- vector("list", length(cross_idx))
      for (e in seq_along(cross_idx)) {
        i <- cross_idx[e]
        fidx <- which(xt[i, ] == 1L)
        ci <- drop(s_mat[i, , drop = FALSE] %*% x_masked)
        b_i <- bmat[fidx, , drop = FALSE]
        g_i <- tcrossprod(b_i)
        g_inv <- try(solve(g_i), silent = TRUE)
        if (inherits(g_inv, "try-error")) {
          return(NULL)
        }
        lam_i <- drop(g_inv %*% (b_i %*% ci))
        lam[i, ] <- 0
        lam[i, fidx] <- lam_i
        theta0[i] <- s_mat[i, i] -
          drop(lam[i, , drop = FALSE] %*% psi %*% lam[i, ])
        cross_info[[e]] <- list(i = i, fidx = fidx, ci = ci,
                                b_i = b_i, g_inv = g_inv, lam_i = lam_i)
      }

      # derivative pass over the local vech directions
      dlam <- matrix(0, nrow = nvar * nfac, ncol = m_total)
      npsic <- nfac * (nfac + 1L) / 2L
      dpsi <- matrix(0, nrow = npsic, ncol = m_total)
      dth <- matrix(0, nrow = nvar, ncol = m_total)
      dth[, cc] <- dtheta
      ii_of <- unlist(lapply(seq_len(nvar), function(j) seq(j, nvar)))
      jj_of <- unlist(lapply(seq_len(nvar), function(j) {
        rep(j, nvar - j + 1L)
      }))
      pat_cells <- which(bi$m_b == 1L)
      for (k in seq_len(pstar)) {
        i <- ii_of[k]
        j <- jj_of[k]
        dthk <- dtheta[, k]
        # dW columns action: dW = dS_k - dTheta_k, with the floored
        # diagonal cells replaced by 0.001 dS_pp
        dwx <- matrix(0, nvar, nfac)
        dwx[i, ] <- xt[j, ]
        if (i != j) {
          dwx[j, ] <- dwx[j, ] + xt[i, ]
        }
        nz <- which(dthk != 0)
        nz <- setdiff(nz, floor_idx) # floored diag cells: no theta term
        if (length(nz) > 0L) {
          dwx[nz, ] <- dwx[nz, ] - dthk[nz] * xt[nz, , drop = FALSE]
        }
        if (i == j && i %in% floor_idx) {
          # floored diagonal: W_pp = 0.001 s_pp
          dwx[i, ] <- dwx[i, ] - xt[i, ] + 0.001 * xt[i, ]
        }
        dt_mat <- crossprod(xt, dwx)
        du_mat <- (dwx - u_mat %*% dt_mat) %*% t_inv
        dmm <- t(t(du_mat) * dsq) +
          t(t(u_mat) * (0.5 * diag(dt_mat) / dsq))
        # rotation chain
        dh <- matrix(0, nfac, nfac)
        for (f in seq_len(nfac)) {
          zrows <- zrows_l[[f]]
          zf <- mm[zrows, , drop = FALSE]
          dzf <- dmm[zrows, , drop = FALSE]
          dmrow <- dmm[bi$marker_idx[f], ]
          dkkt <- rbind(
            cbind(crossprod(dzf, zf) + crossprod(zf, dzf), dmrow),
            c(dmrow, 0)
          )
          solf <- kkt_inv[[f]] %*% c(numeric(nfac), 1)
          dsol <- -kkt_inv[[f]] %*% (dkkt %*% solf)
          dh[, f] <- dsol[seq_len(nfac)]
        }
        dlam_k <- (dmm %*% hmat + mm %*% dh) * bi$m_b
        # mapping chain over the clean rows
        dlam_ck <- dlam_k[clean, , drop = FALSE]
        dth_ck <- dthk[clean]
        dti_c <- -dth_ck * ti_c * ti_c
        dti_c[ti_guard] <- 0
        db_map <- t(dlam_ck * ti_c) + t(lam_c * dti_c)
        da_map <- db_map %*% lam_c + b_map %*% dlam_ck
        dm_map <- a_inv %*% (db_map - da_map %*% m_map)
        # M dW_c M'
        ci_pos <- match(i, clean)
        cj_pos <- match(j, clean)
        if (!is.na(ci_pos) && !is.na(cj_pos)) {
          if (i != j) {
            mdwm <- tcrossprod(m_map[, ci_pos], m_map[, cj_pos])
            mdwm <- mdwm + t(mdwm)
          } else {
            mdwm <- tcrossprod(m_map[, ci_pos])
          }
        } else {
          mdwm <- matrix(0, nfac, nfac)
        }
        nz_c <- which(dth_ck != 0)
        # floored cells: no theta term; diag dS scaled by 0.001
        if (length(floor_idx) > 0L) {
          fl_c <- match(intersect(floor_idx, clean), clean)
          nz_c <- setdiff(nz_c, fl_c)
          if (!is.na(ci_pos) && i == j && (i %in% floor_idx)) {
            mdwm <- mdwm * 0.001
          }
        }
        if (length(nz_c) > 0L) {
          mdwm <- mdwm - (m_map[, nz_c, drop = FALSE] %*%
            (dth_ck[nz_c] * t(m_map[, nz_c, drop = FALSE])))
        }
        dpsi_k <- dm_map %*% wm + mdwm + t(dm_map %*% wm)
        # cross rows
        dam_x <- crossprod(x_masked, dlam_k)
        dbmat <- dpsi_k %*% t(am_x) + psi %*% t(dam_x)
        for (e in seq_along(cross_idx)) {
          cinf <- cross_info[[e]]
          ic <- cinf$i
          fidx <- cinf$fidx
          dci <- numeric(nfac)
          if (ic == i) {
            dci <- dci + x_masked[j, ]
          }
          if (ic == j && i != j) {
            dci <- dci + x_masked[i, ]
          }
          db_i <- dbmat[fidx, , drop = FALSE]
          dg_i <- tcrossprod(db_i, cinf$b_i)
          dg_i <- dg_i + t(dg_i)
          dlam_i <- drop(cinf$g_inv %*%
            (db_i %*% cinf$ci + cinf$b_i %*% dci - dg_i %*% cinf$lam_i))
          dlam_k[ic, ] <- 0
          dlam_k[ic, fidx] <- dlam_i
          # theta identity: theta_i = s_ii - lam_i psi lam_i'
          dth_i <- (i == j && i == ic) * 1 -
            (2 * sum((psi %*% lam[ic, ]) * dlam_k[ic, ]) +
             drop(lam[ic, , drop = FALSE] %*% dpsi_k %*% lam[ic, ]))
          dth[ic, cc[k]] <- dth_i
        }
        # store
        dlam[pat_cells, cc[k]] <- dlam_k[pat_cells]
        for (f in seq_len(nfac)) {
          for (g in seq_len(f)) {
            dpsi[pidx(f, g, nfac), cc[k]] <- dpsi_k[f, g]
          }
        }
      }
      bi_list[[b]]$st <- list(lam = lam, dlam = dlam, psi = psi,
                              dpsi = dpsi, theta = theta0, dth = dth)
    }
  } else if (!pooled) {
    # ---- classic per-block chain --------------------------------------
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
        b_map <- t(lam_full * ti)
        a_map <- b_map %*% lam_full
        a_inv <- try(solve(a_map), silent = TRUE)
        if (inherits(a_inv, "try-error")) {
          return(NULL)
        }
        m_map <- a_inv %*% b_map
        wm <- w_mat %*% t(m_map)
        mtm <- m_map %*% (theta0 * t(m_map))
      }

      lam <- t(t(u_mat) / u_vec) * bi$m_b
      psi <- if (psi_mapping) {
        m_map %*% wm
      } else {
        t(t_mat * u_vec) * u_vec # T_fg u_f u_g
      }

      dlam <- matrix(0, nrow = nvar * nfac, ncol = m_total)
      npsic <- nfac * (nfac + 1L) / 2L
      dpsi <- matrix(0, nrow = npsic, ncol = m_total)
      dth <- matrix(0, nrow = nvar, ncol = m_total)
      dth[, cc] <- dtheta
      pat_cells <- which(bi$m_b == 1L)
      pat_ov <- ((pat_cells - 1L) %% nvar) + 1L
      pat_fac <- ((pat_cells - 1L) %/% nvar) + 1L
      lam_cell_val <- lam[pat_cells]
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
        dlam[pat_cells, cc[k]] <-
          (du_mat[pat_cells] - lam_cell_val * du_vec[pat_fac]) /
          u_vec[pat_fac]
        if (psi_mapping) {
          dlamk <- (du_mat - t(t(u_mat) * (du_vec / u_vec))) *
            bi$m_b / rep(u_vec, each = nvar)
          dti <- -dthk * ti * ti
          dti[ti_guard] <- 0
          db_map <- t(dlamk * ti) + t(lam_full * dti)
          da_map <- db_map %*% lam_full + b_map %*% dlamk
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
          dpsi_k <- dm_map %*% wm + mdwm + t(dm_map %*% wm)
          for (f in seq_len(nfac)) {
            for (g in seq_len(f)) {
              dpsi[pidx(f, g, nfac), cc[k]] <- dpsi_k[f, g]
            }
          }
        } else {
          for (f in seq_len(nfac)) {
            for (g in seq_len(f)) {
              dpsi[pidx(f, g, nfac), cc[k]] <-
                dt_mat[f, g] * u_vec[f] * u_vec[g] +
                t_mat[f, g] *
                  (du_vec[f] * u_vec[g] + u_vec[f] * du_vec[g])
            }
          }
        }
      }
      bi_list[[b]]$st <- list(lam = lam, dlam = dlam, psi = psi,
                              dpsi = dpsi, theta = theta0, dth = dth)
    }
  } else {
    # ---- pooled branch (loading ties / quadprog) -----------------------
    shrink <- rep(1, nblocks)
    dshrink <- NULL
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
      if (sort(roots)[2L] - roots[bstar] <
          1e-8 * max(1, abs(roots[bstar]))) {
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
      dvec <- 1 / sqrt(diag(t_mat))
      phi <- t(t_mat * dvec) * dvec
      evp <- eigen(phi, symmetric = TRUE, only.values = TRUE)$values
      if (min(evp) / max(evp) < 2e-4) {
        return(NULL) # force_pd about to clip
      }
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
        row <- row - shrink[b] * dth_g[p, ]
        if (!is.null(dshrink)) {
          row <- row - bi$theta0[p] * dshrink
        }
        dwx_rows[e, ] <- row
      }
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

      if (length(bi$lam_pt_idx) > 0L) {
        cells <- rbind(cells, data.frame(
          block = b, ov = bi$lam_ov, fac = bi$lam_fac,
          facname = lavpartable$lhs[bi$lam_pt_idx],
          free = free[bi$lam_pt_idx]
        ))
      }
    }
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
    for (fn in unique(cells$facname)) {
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

    # cell values (tied classes averaged) -> the assembled lambda
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
    # assemble the stage-1 structures per block
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      nvar <- bi$nvar
      nfac <- bi$nfac
      lam <- matrix(0, nvar, nfac)
      lam[(seq_len(nfac) - 1L) * nvar + bi$marker_idx] <- 1
      dlam <- matrix(0, nrow = nvar * nfac, ncol = m_total)
      rows_b <- which(cells$block == b)
      for (r in rows_b) {
        mh <- mh_val[[b]][[cells$facname[r]]]
        dmh <- mh_row[[b]][[cells$facname[r]]]
        lv <- val[r] / mh
        cell <- (cells$fac[r] - 1L) * nvar + cells$ov[r]
        lam[cell] <- lv
        dlam[cell, ] <- (vrow[r, ] - lv * dmh) / mh
      }
      # psi: T_fg (WX)_{mf, f} (WX)_{mg, g} / (T_ff T_gg) -- own mhat
      psi <- matrix(0, nfac, nfac)
      npsic <- nfac * (nfac + 1L) / 2L
      dpsi <- matrix(0, nrow = npsic, ncol = m_total)
      for (f in seq_len(nfac)) {
        for (g in seq_len(f)) {
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
          psi[f, g] <- psi[g, f] <- tfg * aa
          dpsi[pidx(f, g, nfac), ] <-
            bi$dt_rows[[(f - 1L) * nfac + g]] * aa +
            tfg * ((dwxf * wxg + wxf * dwxg) / (tff * tgg) -
                   aa * (bi$dt_rows[[(f - 1L) * nfac + f]] / tff +
                         bi$dt_rows[[(g - 1L) * nfac + g]] / tgg))
        }
      }
      dth <- matrix(0, nrow = nvar, ncol = m_total)
      dth[, co$cov_vech[[b]]] <- bi$dtheta
      bi_list[[b]]$st <- list(lam = lam, dlam = dlam, psi = psi,
                              dpsi = dpsi, theta = bi$theta0, dth = dth)
    }
  }

  # --------------------------------------------------------------------
  # metric rescale (std.lv, fixed loading c != 1, tied factors):
  # lambda[, f] * s_f, psi / (s_f s_g); s from the anchor blocks
  # --------------------------------------------------------------------
  need_rescale <- any(vapply(bi_list, function(bi) {
    any(bi$scale_type != "marker") || any(bi$marker_value != 1)
  }, logical(1L)))
  psi_pin_loose <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    psi_pin_loose[[b]] <- integer(0L)
  }
  if (need_rescale) {
    # component structure over the loading ties (empty tie set for the
    # classic/cross branches gives singleton components)
    cells_rs <- NULL
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      if (length(bi$lam_pt_idx) > 0L) {
        cells_rs <- rbind(cells_rs, data.frame(
          block = b, facname = lavpartable$lhs[bi$lam_pt_idx],
          free = free[bi$lam_pt_idx]
        ))
      }
    }
    tie_rs <- split(seq_len(nrow(cells_rs)), cells_rs$free)
    tie_rs <- tie_rs[vapply(tie_rs, length, 1L) > 1L]
    comp_rs <- function(fn) {
      comps <- as.list(seq_len(nblocks))
      for (class_idx in tie_rs) {
        if (cells_rs$facname[class_idx[1L]] != fn) {
          next
        }
        bset <- unique(cells_rs$block[class_idx])
        hit <- which(vapply(comps, function(x) any(x %in% bset), TRUE))
        if (length(hit) > 1L) {
          comps[[hit[1L]]] <- sort(unique(unlist(comps[hit])))
          comps[hit[-1L]] <- NULL
        }
      }
      comps
    }
    all_fn <- unique(unlist(lapply(bi_list, "[[", "lv_names")))
    for (fn in all_fn) {
      comps <- comp_rs(fn)
      for (ccomp in comps) {
        bset <- ccomp[vapply(ccomp, function(b) {
          fn %in% bi_list[[b]]$lv_names
        }, logical(1L))]
        if (length(bset) == 0L) {
          next
        }
        f_of <- vapply(bset, function(b) {
          match(fn, bi_list[[b]]$lv_names)
        }, integer(1L))
        types <- vapply(seq_along(bset), function(k) {
          bi_list[[bset[k]]]$scale_type[f_of[k]]
        }, character(1L))
        var_k <- which(types == "variance")
        ds_row <- NULL # NULL = constant s
        if (any(types == "marker")) {
          s <- mean(vapply(which(types == "marker"), function(k) {
            bi_list[[bset[k]]]$marker_value[f_of[k]]
          }, numeric(1L)))
          loose_k <- var_k
        } else {
          if (length(var_k) == 0L) {
            return(NULL) # unanchored: the estimation stops on this
          }
          ratios <- vapply(var_k, function(k) {
            b <- bset[k]
            f <- f_of[k]
            bi_list[[b]]$st$psi[f, f] / bi_list[[b]]$scale_target[f]
          }, numeric(1L))
          if (any(ratios <= 0)) {
            return(NULL)
          }
          s <- sqrt(mean(ratios))
          ds_row <- numeric(m_total)
          for (k in var_k) {
            b <- bset[k]
            f <- f_of[k]
            ds_row <- ds_row +
              bi_list[[b]]$st$dpsi[
                pidx(f, f, bi_list[[b]]$nfac), ] /
              bi_list[[b]]$scale_target[f]
          }
          ds_row <- ds_row / (length(var_k) * 2 * s)
          loose_k <- if (length(var_k) > 1L) var_k else integer(0L)
        }
        for (k in loose_k) {
          psi_pin_loose[[bset[k]]] <-
            c(psi_pin_loose[[bset[k]]], f_of[k])
        }
        if (s == 1 && is.null(ds_row)) {
          next
        }
        # apply to every block in the component
        for (k in seq_along(bset)) {
          b <- bset[k]
          f <- f_of[k]
          bi <- bi_list[[b]]
          nvar <- bi$nvar
          nfac <- bi$nfac
          st <- bi$st
          col_cells <- (f - 1L) * nvar + seq_len(nvar)
          old_col <- st$lam[, f]
          st$lam[, f] <- old_col * s
          if (is.null(ds_row)) {
            st$dlam[col_cells, ] <- st$dlam[col_cells, , drop = FALSE] * s
          } else {
            st$dlam[col_cells, ] <-
              st$dlam[col_cells, , drop = FALSE] * s +
              outer(old_col, ds_row)
          }
          for (g in seq_len(nfac)) {
            pp <- pidx(f, g, nfac)
            sfg <- if (g == f) s * s else s
            old_val <- st$psi[f, g]
            new_val <- old_val / sfg
            st$psi[f, g] <- st$psi[g, f] <- new_val
            if (is.null(ds_row)) {
              st$dpsi[pp, ] <- st$dpsi[pp, ] / sfg
            } else {
              mult <- if (g == f) 2 else 1
              st$dpsi[pp, ] <- st$dpsi[pp, ] / sfg -
                (new_val * mult / s) * ds_row
            }
          }
          bi_list[[b]]$st <- st
        }
      }
    }
  }

  # --------------------------------------------------------------------
  # restricted PSI fit: implicit function theorem at the ML optimum
  # --------------------------------------------------------------------
  psi_ids_vec <- unlist(lapply(bi_list, function(bi) {
    bi$psi_id[upper.tri(bi$psi_id, diag = TRUE) & bi$psi_id > 0L]
  }))
  psi_restricted <- any(duplicated(psi_ids_vec)) ||
    any(vapply(seq_len(nblocks), function(b) {
      length(psi_pin_loose[[b]]) > 0L
    }, logical(1L))) ||
    any(vapply(bi_list, function(bi) {
      bi$nfac > 1L && any(bi$psi_id[upper.tri(bi$psi_id)] == 0L)
    }, logical(1L)))
  if (psi_restricted) {
    ids <- sort(unique(psi_ids_vec))
    id_list <- lapply(bi_list, function(bi) {
      matrix(match(bi$psi_id, ids, nomatch = 0L),
             bi$nfac, bi$nfac)
    })
    pin_list <- lapply(bi_list, "[[", "psi_pin")
    for (b in seq_len(nblocks)) {
      for (f in psi_pin_loose[[b]]) {
        id_list[[b]][f, f] <- 0L
        pin_list[[b]][f, f] <- bi_list[[b]]$scale_target[f]
      }
    }
    nid <- length(ids)
    w_fit <- vapply(bi_list, "[[", numeric(1L), "nobs")
    w_fit <- w_fit / sum(w_fit)
    if (nid > 0L) {
      # replicate the fit (also feeds the safety net through x_pred)
      x0fit <- vapply(seq_len(nid), function(k) {
        vals <- numeric(0L)
        for (b in seq_len(nblocks)) {
          hit <- which(id_list[[b]] == k)
          if (length(hit) > 0L) {
            vals <- c(vals, bi_list[[b]]$st$psi[hit])
          }
        }
        mean(vals)
      }, numeric(1L))
      xfit <- lav_cfa_mgm_psi_fit(
        veta = lapply(bi_list, function(bi) bi$st$psi),
        id = id_list, pin = pin_list, w = w_fit, x0 = x0fit
      )
      # IFT: H dx = -(dg/dV) dV, with
      # g_k = sum_b w_b tr[(Psi^-1 - Psi^-1 V Psi^-1) E_k]
      psi_hat <- vector("list", nblocks)
      psi_hat_inv <- vector("list", nblocks)
      for (b in seq_len(nblocks)) {
        ph <- pin_list[[b]]
        idx <- which(id_list[[b]] > 0L)
        ph[idx] <- xfit[id_list[[b]][idx]]
        ch <- try(chol(ph), silent = TRUE)
        if (inherits(ch, "try-error")) {
          return(NULL)
        }
        psi_hat[[b]] <- ph
        psi_hat_inv[[b]] <- chol2inv(ch)
      }
      h_mat <- lav_cfa_mgm_psi_fit_hessian(
        x = xfit, veta = lapply(bi_list, function(bi) bi$st$psi),
        id = id_list, pin = pin_list, w = w_fit
      )
      dg <- matrix(0, nrow = nid, ncol = m_total)
      for (b in seq_len(nblocks)) {
        nfac <- bi_list[[b]]$nfac
        pinv <- psi_hat_inv[[b]]
        for (k in seq_len(nid)) {
          ek <- matrix(0, nfac, nfac)
          ek[id_list[[b]] == k] <- 1
          if (all(ek == 0)) {
            next
          }
          a_k <- pinv %*% ek %*% pinv
          # dg_k over the moments: -w_b tr[A_k dV_b]
          for (f in seq_len(nfac)) {
            for (g in seq_len(f)) {
              mult <- if (f == g) 1 else 2
              dg[k, ] <- dg[k, ] - w_fit[b] * mult * a_k[f, g] *
                bi_list[[b]]$st$dpsi[pidx(f, g, nfac), ]
            }
          }
        }
      }
      dxfit <- try(-solve(h_mat, dg), silent = TRUE)
      if (inherits(dxfit, "try-error")) {
        return(NULL)
      }
      # write the fitted psi (values + rows) back into the blocks
      for (b in seq_len(nblocks)) {
        nfac <- bi_list[[b]]$nfac
        st <- bi_list[[b]]$st
        st$psi <- psi_hat[[b]]
        for (f in seq_len(nfac)) {
          for (g in seq_len(f)) {
            k <- id_list[[b]][f, g]
            st$dpsi[pidx(f, g, nfac), ] <- if (k > 0L) {
              dxfit[k, ]
            } else {
              0
            }
          }
        }
        bi_list[[b]]$st <- st
      }
      # keep track of the fitted values for the safety net
      x_pred[ids] <- xfit
    } else {
      # PSI completely fixed
      for (b in seq_len(nblocks)) {
        st <- bi_list[[b]]$st
        st$psi <- pin_list[[b]]
        st$dpsi[] <- 0
        bi_list[[b]]$st <- st
      }
    }
    # refresh the cross-loading thetas (they depend on PSI)
    for (b in seq_len(nblocks)) {
      bi <- bi_list[[b]]
      if (length(bi$cross_idx) == 0L) {
        next
      }
      st <- bi$st
      nfac <- bi$nfac
      cc <- co$cov_vech[[b]]
      vpos_b <- function(i, j) {
        (j - 1L) * bi$nvar - ((j - 1L) * j) %/% 2L + i
      }
      for (i in bi$cross_idx) {
        lam_i <- st$lam[i, ]
        pl <- drop(st$psi %*% lam_i)
        st$theta[i] <- bi$s_use[i, i] - sum(lam_i * pl)
        row <- numeric(m_total)
        row[cc[vpos_b(i, i)]] <- 1
        for (f in seq_len(nfac)) {
          row <- row - 2 * pl[f] * st$dlam[(f - 1L) * bi$nvar + i, ]
          for (g in seq_len(f)) {
            mult <- if (f == g) 1 else 2
            row <- row - mult * lam_i[f] * lam_i[g] *
              st$dpsi[pidx(f, g, nfac), ]
          }
        }
        st$dth[i, ] <- row
      }
      bi_list[[b]]$st <- st
    }
  }

  # --------------------------------------------------------------------
  # write the stage-1 rows into dx (+ predicted values)
  # --------------------------------------------------------------------
  for (b in seq_len(nblocks)) {
    bi <- bi_list[[b]]
    st <- bi$st
    nvar <- bi$nvar
    nfac <- bi$nfac
    for (r in seq_along(bi$lam_pt_idx)) {
      cell <- (bi$lam_fac[r] - 1L) * nvar + bi$lam_ov[r]
      dx[free[bi$lam_pt_idx[r]], ] <- st$dlam[cell, ]
      x_pred[free[bi$lam_pt_idx[r]]] <- st$lam[cell]
    }
    if (length(bi$th_pt_idx) > 0L) {
      th_ov <- match(lavpartable$lhs[bi$th_pt_idx], bi$ov_names)
      dx[free[bi$th_pt_idx], ] <- st$dth[th_ov, , drop = FALSE]
      x_pred[free[bi$th_pt_idx]] <- st$theta[th_ov]
    }
    if (length(bi$psi_rows_pt) > 0L) {
      pf <- match(lavpartable$lhs[bi$psi_rows_pt], bi$lv_names)
      pg <- match(lavpartable$rhs[bi$psi_rows_pt], bi$lv_names)
      for (r in seq_along(bi$psi_rows_pt)) {
        if (free[bi$psi_rows_pt[r]] == 0L) {
          next
        }
        pp <- pidx(pf[r], pg[r], nfac)
        dx[free[bi$psi_rows_pt[r]], ] <- st$dpsi[pp, ]
        x_pred[free[bi$psi_rows_pt[r]]] <- st$psi[pf[r], pg[r]]
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

# summary() header rows describing the standard errors of the noniterative
# estimators (MGM, JS/JSA, IV). Every analytic standard-error flavor in
# this family is a first-order delta-method (sandwich) expression over the
# sample moments, J Gamma J' / n, so the common label is "Delta"; the
# remaining rows show the flavor of the moment ACOV (Gamma) and, for the
# IV estimator, the methods used for the two parts of the model. Returns
# a named character vector (names = the left-hand labels), or NULL when
# the generic header applies (bootstrap/external, or another estimator).
lav_noniter_se_rows <- function(object = NULL) {
  estimator <- object@Options$estimator
  se <- object@Options$se
  if (!estimator %in% c("MGM", "JS", "JSA", "IV") ||
      !se %in% c("standard", "robust")) {
    return(NULL)
  }
  ea <- object@Options$estimator.args
  two_stage <- any(object@Options$missing ==
                   c("two.stage", "robust.two.stage"))

  out <- c("Standard errors" = "Delta")
  if (estimator == "IV") {
    out <- c(out,
      "Regression part (stage 1)" = ea[["iv_vcov_stage1"]],
      "Variance part (stage 2)"   = ea[["iv_vcov_stage2"]])
    # type of Gamma (moment ACOV) used in the stage-2 standard errors;
    # not applicable when no stage-2 covariance is computed
    if (!identical(tolower(ea[["iv_vcov_stage2"]]), "none")) {
      out <- c(out, "Gamma matrix" = if (object@Model@categorical) {
        "ADF"
      } else if (two_stage) {
        "TS"
      } else {
        "NT"
      })
    }
  } else if (estimator %in% c("JS", "JSA")) {
    # under two-stage missing data the moment ACOV is the two-stage one,
    # whatever js_gamma says (see lav_sem_js_vcov)
    gamma_flavor <- toupper(ea[["js_gamma"]])
    if (two_stage) {
      gamma_flavor <- "TS"
    } else if (gamma_flavor == "ADF" && object@Data@data.type != "full") {
      gamma_flavor <- "NT" # the adf request fell back at fit time
    }
    out <- c(out, "Gamma matrix" = gamma_flavor)
  } else { # MGM
    gamma_flavor <- toupper(ea[["mgm.gamma"]])
    if (gamma_flavor == "ADF" && object@Data@data.type != "full") {
      gamma_flavor <- "NT" # the adf request fell back at fit time
    }
    out <- c(out, "Gamma matrix" = gamma_flavor)
  }

  out
}
