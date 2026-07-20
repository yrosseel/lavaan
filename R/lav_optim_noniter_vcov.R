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

# ANALYTIC Jacobian of the classic MGM branch (simple structure, no
# equality ties, no residual covariances, marker-scaled factors,
# saturated PSI, no second stage), per block:
#
#   theta_p = clip( s_pp - mean_{a<b} s_pa s_pb / s_ab )   (Spearman)
#   W       = S - rho_eff * diag(theta)                    (pd correction)
#   T       = X' W X,   U = W X T^{-1},   u_f = U[marker_f, f]
#   lambda_pf = U_pf / u_f
#   psi_fg    = T_fg u_f u_g            (sum-score psi; the normalized
#                                        metric cancels exactly)
#   psi       = M W M', M = (L'KL)^{-1} L'K   (psi.mapping = TRUE)
#
# all smooth rational functions of the sample moments; the pd-correction
# root rho is a simple generalized eigenvalue of the pencil (S, Theta)
# with d rho = v' (dS - rho dTheta) v / (v' Theta v). Regime switches
# (Spearman clips, theta bounds, pd correction active or not, the
# zero-theta guard of the mapping) are frozen at the point estimate --
# exactly what the numerical Jacobian sees locally.
#
# returns the npar x (stacked moments) Jacobian, or NULL when the model
# uses a branch this derivation does not cover (the caller then falls
# back to the numerical Jacobian)
lav_cfa_mgm_jacobian <- function(lavmodel = NULL, lavsamplestats = NULL,
                                 lavoptions = NULL, lavpartable = NULL,
                                 implied0 = NULL) {
  ea <- lavoptions$estimator.args

  # not covered: the quadprog (constrained-pooled) stage-1 branch
  if (isTRUE(ea[["quadprog"]])) {
    return(NULL)
  }
  psi_mapping <- isTRUE(ea[["psi.mapping"]])
  # not covered: the second stage for the variances/covariances
  mgm_varcov <- ea[["mgm.varcov"]]
  if (is.null(mgm_varcov)) {
    mgm_varcov <- "default"
  }
  if (!toupper(mgm_varcov) %in% c("DEFAULT", "NONE")) {
    return(NULL)
  }
  # not covered: equality ties (the pooled stage-1 branch)
  free <- lavpartable$free
  if (any(duplicated(free[free > 0L]))) {
    return(NULL)
  }
  # not covered: psi.mapping over an unzeroed (EFA-style) lambda
  zae <- ea[["zero.after.efa"]]
  if (psi_mapping && !is.null(zae) && !isTRUE(zae)) {
    return(NULL)
  }

  lavpta <- lav_pt_attributes(lavpartable)
  nblocks <- lav_pt_nblocks(lavpartable)
  pt_block <- if (!is.null(lavpartable$block)) {
    lavpartable$block
  } else {
    rep(1L, length(lavpartable$op))
  }

  npar <- lavmodel@nx.free
  ndtotal <- 0L
  for (b in seq_len(nblocks)) {
    nvar_b <- lavmodel@nvar[b]
    ndtotal <- ndtotal + nvar_b * (nvar_b + 1L) / 2L +
      if (lavmodel@meanstructure) nvar_b else 0L
  }
  jac <- matrix(0, nrow = npar, ncol = ndtotal)
  # the values these formulas predict; compared against the actual
  # estimates at the end -- any mismatch (an active bound clip, or a
  # branch the coverage checks above missed) triggers the numerical
  # fallback
  x_pred <- rep(NA_real_, npar)
  col0 <- 0L

  for (b in seq_len(nblocks)) {
    ov_names <- lavpta$vnames$ov[[b]]
    lv_names <- lavpta$vnames$lv.regular[[b]]
    nvar <- length(ov_names)
    nfac <- length(lv_names)
    pstar <- nvar * (nvar + 1L) / 2L
    mean_off <- if (lavmodel@meanstructure) nvar else 0L
    nobs_b <- lavsamplestats@nobs[[b]]
    s_mat <- implied0$cov[[b]]

    # ---- coverage checks for this block --------------------------------
    # markers fixed to 1 for every factor; no other fixed nonzero loadings
    marker_idx <- integer(nfac)
    for (f in seq_len(nfac)) {
      fx <- which(lavpartable$op == "=~" & pt_block == b &
        lavpartable$lhs == lv_names[f] & lavpartable$free == 0L &
        !is.na(lavpartable$ustart) & lavpartable$ustart != 0)
      if (length(fx) != 1L || lavpartable$ustart[fx] != 1) {
        return(NULL) # std.lv / fixed value c != 1 / multiple markers
      }
      marker_idx[f] <- match(lavpartable$rhs[fx], ov_names)
    }
    # pattern; no cross-loadings
    lam_pt_idx <- which(lavpartable$op == "=~" & pt_block == b &
      lavpartable$free > 0L)
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
    # every factor needs >= 3 indicators (plain Spearman path)
    if (any(colSums(m_b) < 3L)) {
      return(NULL)
    }
    # no residual covariances (free, or fixed nonzero)
    rescov_idx <- which(lavpartable$op == "~~" & pt_block == b &
      lavpartable$lhs %in% ov_names & lavpartable$rhs %in% ov_names &
      lavpartable$lhs != lavpartable$rhs &
      (lavpartable$free > 0L |
        (!is.na(lavpartable$ustart) & lavpartable$ustart != 0)))
    if (length(rescov_idx) > 0L) {
      return(NULL)
    }
    # saturated PSI: every lv-lv cell present and free
    psi_rows_pt <- which(lavpartable$op == "~~" & pt_block == b &
      lavpartable$lhs %in% lv_names & lavpartable$rhs %in% lv_names)
    if (length(psi_rows_pt) != nfac * (nfac + 1L) / 2L ||
        any(lavpartable$free[psi_rows_pt] == 0L)) {
      return(NULL) # fixed variance/covariance -> rescale or psi fit
    }
    # mean structure: only the saturated case (every observed intercept
    # free, no free latent means) reduces to nu = ybar exactly
    nu_rows_pt <- integer(0L)
    if (lavmodel@meanstructure) {
      m1_idx <- which(lavpartable$op == "~1" & pt_block == b &
        lavpartable$free > 0L)
      if (any(!lavpartable$lhs[m1_idx] %in% ov_names)) {
        return(NULL) # free latent means
      }
      if (length(m1_idx) != nvar) {
        return(NULL) # some intercepts fixed
      }
      nu_rows_pt <- m1_idx
    }

    # ---- theta: Spearman tetrads + clips -------------------------------
    # vech position of cell (i, j), lower triangle, column-major
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
          dtheta[p, ] <- 0
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

    # ---- pd correction: smallest root of |S - rho Theta| = 0 -----------
    zero_idx <- which(abs(theta0) < sqrt(.Machine$double.eps))
    rho_eff <- 1
    drho <- NULL # NULL = inactive (plain S - Theta)
    if (length(zero_idx) < nvar) {
      pos <- setdiff(seq_len(nvar), zero_idx)
      a_mat <- s_mat[pos, pos, drop = FALSE]
      if (length(zero_idx) > 0L) {
        s22i_s21 <- try(solve(
          s_mat[zero_idx, zero_idx, drop = FALSE],
          s_mat[zero_idx, pos, drop = FALSE]
        ), silent = TRUE)
        if (inherits(s22i_s21, "try-error")) {
          return(NULL)
        }
        a_mat <- a_mat -
          s_mat[pos, zero_idx, drop = FALSE] %*% s22i_s21
      }
      ld <- 1 / sqrt(theta0[pos])
      ee <- eigen(t(ld * a_mat) * ld, symmetric = TRUE)
      np <- length(pos)
      rho <- ee$values[np]
      # a (near-)multiple smallest root is not differentiable
      if (np > 1L &&
          (ee$values[np - 1L] - rho) < 1e-8 * max(1, abs(rho))) {
        return(NULL)
      }
      cutoff <- 1 + 1 / (nobs_b - 1)
      if (rho < cutoff) {
        rho_eff <- rho - 1 / (nobs_b - 1)
        # full-space pencil eigenvector: S v = rho Theta v
        v <- numeric(nvar)
        v[pos] <- ld * ee$vectors[, np]
        if (length(zero_idx) > 0L) {
          v[zero_idx] <- -drop(s22i_s21 %*% v[pos])
        }
        denom <- sum(theta0 * v * v)
        # d rho = v' (dS - rho dTheta) v / (v' Theta v)
        vv <- tcrossprod(v)
        dquad <- 2 * lav_mat_vech(vv)
        dquad[lav_mat_diagh_idx(nvar)] <- diag(vv)
        drho <- (dquad - drop((v * v) %*% dtheta) * rho) / denom
      }
    }

    # ---- the classic solve: T, U, and the parameter maps ---------------
    w_mat <- s_mat - rho_eff * diag(theta0, nvar)
    xt <- m_b # numeric use
    t_mat <- crossprod(xt, w_mat %*% xt)
    t_inv <- try(solve(t_mat), silent = TRUE)
    if (inherits(t_inv, "try-error")) {
      return(NULL)
    }
    u_mat <- w_mat %*% xt %*% t_inv
    u_vec <- u_mat[cbind(marker_idx, seq_len(nfac))]
    if (any(abs(u_vec) < .Machine$double.eps)) {
      return(NULL)
    }
    th0x <- theta0 * xt

    # psi.mapping ingredients (marker-metric lambda, mapping matrix)
    if (psi_mapping) {
      lam_full <- t(t(u_mat) / u_vec) * m_b
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

    # free-parameter rows of this block
    th_pt_idx <- which(lavpartable$op == "~~" & pt_block == b &
      lavpartable$lhs %in% ov_names &
      lavpartable$lhs == lavpartable$rhs & lavpartable$free > 0L)
    th_ov <- match(lavpartable$lhs[th_pt_idx], ov_names)
    psi_f <- match(lavpartable$lhs[psi_rows_pt], lv_names)
    psi_g <- match(lavpartable$rhs[psi_rows_pt], lv_names)

    # theta rows (directly from the Spearman chain)
    if (length(th_pt_idx) > 0L) {
      jac[lavpartable$free[th_pt_idx],
          col0 + mean_off + seq_len(pstar)] <- dtheta[th_ov, , drop = FALSE]
      x_pred[lavpartable$free[th_pt_idx]] <- theta0[th_ov]
    }
    # nu rows: nu = ybar exactly
    if (length(nu_rows_pt) > 0L) {
      nu_ov <- match(lavpartable$lhs[nu_rows_pt], ov_names)
      jac[cbind(lavpartable$free[nu_rows_pt], col0 + nu_ov)] <- 1
      x_pred[lavpartable$free[nu_rows_pt]] <- implied0$mean[[b]][nu_ov]
    }

    # loading and psi rows: one pass over the vech directions
    nlam <- length(lam_pt_idx)
    lam_free <- lavpartable$free[lam_pt_idx]
    lam_val <- u_mat[cbind(lam_ov, lam_fac)] / u_vec[lam_fac]
    psi_free <- lavpartable$free[psi_rows_pt]
    x_pred[lam_free] <- lam_val
    if (psi_mapping) {
      psi_map_full <- m_map %*% wm
      x_pred[psi_free] <- psi_map_full[cbind(psi_f, psi_g)]
    } else {
      x_pred[psi_free] <- t_mat[cbind(psi_f, psi_g)] *
        u_vec[psi_f] * u_vec[psi_g]
    }
    # (i, j) cell of vech direction k: column-major lower triangle
    ii_of <- unlist(lapply(seq_len(nvar), function(j) seq(j, nvar)))
    jj_of <- unlist(lapply(seq_len(nvar), function(j) rep(j, nvar - j + 1L)))
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
      du_vec <- du_mat[cbind(marker_idx, seq_len(nfac))]
      # loadings
      if (nlam > 0L) {
        jac[lam_free, col0 + mean_off + k] <-
          (du_mat[cbind(lam_ov, lam_fac)] -
           lam_val * du_vec[lam_fac]) / u_vec[lam_fac]
      }
      # psi
      if (psi_mapping) {
        dlam <- (du_mat - t(t(u_mat) * (du_vec / u_vec))) *
          m_b / rep(u_vec, each = nvar)
        dti <- -dthk * ti * ti
        dti[ti_guard] <- 0
        db_map <- t(dlam * ti) + t(lam_full * dti)
        da_map <- db_map %*% lam_full + b_map %*% dlam
        dm_map <- a_inv %*% (db_map - da_map %*% m_map)
        # M dW M' (dW = dS_k - drho_k Theta - rho_eff dTheta_k)
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
        jac[psi_free, col0 + mean_off + k] <-
          dpsi_mat[cbind(psi_f, psi_g)]
      } else {
        du_f <- du_vec[psi_f]
        du_g <- du_vec[psi_g]
        jac[psi_free, col0 + mean_off + k] <-
          dt_mat[cbind(psi_f, psi_g)] * u_vec[psi_f] * u_vec[psi_g] +
          t_mat[cbind(psi_f, psi_g)] *
            (du_f * u_vec[psi_g] + u_vec[psi_f] * du_g)
      }
    }

    col0 <- col0 + mean_off + pstar
  }

  # verify: every free parameter must be covered, and the predicted
  # estimates must match the actual ones (they differ when a bound clip
  # or a branch the checks above missed fired at estimation time)
  if (anyNA(x_pred)) {
    return(NULL)
  }
  x_act <- as.numeric(lav_model_get_parameters(lavmodel))
  if (max(abs(x_pred - x_act)) > 1e-6 * (1 + max(abs(x_act)))) {
    return(NULL)
  }

  jac
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
