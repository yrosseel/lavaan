# analytic Jacobian of the James-Stein (JS/JSA) estimation map with respect
# to the sample moments -- used by the delta-method standard errors
# (lav_sem_js_vcov). This replaces the numerical (numDeriv) Jacobian of
# lav_sem_js_estimate(): one forward-mode chain-rule sweep instead of
# thousands of map evaluations, so the standard errors keep the
# 'non-iterative' feel of the estimator.
#
# The chain runs through every ingredient of the map:
#   (1) the Spearman/communality reliability estimates (dtheta; computed
#       inside lav_sem_js_theta, see lav_sem_js.R),
#   (2) the alphamax aggregation weights (generalized-eigenvector
#       derivative via a bordered system),
#   (3) the per-equation least-squares solves on the conditional
#       expectations (all in moment form),
#   (4) the pooled (system) solve for equality-constrained directed
#       coefficients (plain or KKT),
#   (5) the second-stage variance/covariance solve (ULS/GLS/2RLS/RLS); for
#       RLS the derivative is taken at the fixed point via the implicit
#       function theorem,
#   (6) the joint GLS mean solve.
#
# Branch decisions that depend on the data (proxy-set fallbacks,
# equal-weights fallbacks, bound clippings) are *recorded* during one
# re-run of the point-estimation map (record = TRUE) and replayed here, so
# the derivative always matches the code path actually taken.
#
# YR 19 Jul 2026: first version

# main entry: the Jacobian of the free-parameter vector with respect to the
# stacked moment vector (as produced by lav_implied_to_vec), evaluated at
# the sample (lavh1) moments. Returns NULL when this implementation does
# not cover the model configuration (the caller then falls back to the
# numerical Jacobian).
lav_sem_js_jacobian <- function(lavmodel = NULL, lavsamplestats = NULL,
                                lavoptions = NULL, lavpartable = NULL,
                                lavh1 = NULL, eqs = NULL) {
  lavpta <- lav_pt_attributes(lavpartable)
  lavpartable <- lav_pt_set_cache(lavpartable, lavpta)
  aggregated <- identical(lavoptions$estimator, "JSA")

  # directed versus undirected (free) parameters
  undirected_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) &
    lavpartable$op == "~~")
  directed_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) &
    lavpartable$op %in% c("=~", "~", "~1"))
  free_directed_idx <- unique(lavpartable$free[directed_idx])
  free_undirected_idx <- unique(lavpartable$free[undirected_idx])

  js_varcov_method <- toupper(lavoptions$estimator.args[["js_varcov_method"]])
  if (length(js_varcov_method) == 0L) {
    js_varcov_method <- "RLS"
  }
  js_mean_structure <- tolower(lavoptions$estimator.args[["js_mean_structure"]])
  if (length(js_mean_structure) == 0L) {
    js_mean_structure <- "wls"
  }

  pt_block <- if (!is.null(lavpartable$block)) {
    lavpartable$block
  } else {
    rep(1L, length(lavpartable$op))
  }
  nblocks <- lavmodel@nblocks

  # not covered (yet): free group weights (extra rows in Delta / the moment
  # vector); fall back to the numerical Jacobian
  if (lavmodel@group.w.free) {
    return(NULL)
  }

  # not covered (yet): general linear equality constraints among the
  # undirected parameters (the stage-2 projection step); fall back to the
  # numerical Jacobian
  if (length(free_undirected_idx) > 0L && js_varcov_method != "NONE") {
    for (b in seq_len(nblocks)) {
      fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
        lavpartable$free > 0L & pt_block == b])
      fu_b <- fu_b[fu_b %in% free_undirected_idx]
      if (length(fu_b) > 0L &&
          !is.null(lav_sem_miiv_linear_con(lavmodel, fu_b))) {
        return(NULL)
      }
    }
  }

  if (is.null(eqs)) {
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
  }

  # 1. re-run the point-estimation map once, recording the resolved
  #    structure (proxy sets, weights, branch decisions, theta derivative)
  x0 <- lav_sem_js_estimate(
    vec = NULL, eqs = eqs,
    lavmodel = lavmodel, lavpartable = lavpartable,
    lavsamplestats = lavsamplestats, lavh1 = lavh1,
    lavoptions = lavoptions,
    free_directed_idx = free_directed_idx,
    free_undirected_idx = free_undirected_idx,
    aggregated = aggregated, record = TRUE
  )
  eqs_rec <- attr(x0, "eqs")
  js_trace <- attr(x0, "js_trace")
  x0 <- as.numeric(x0)
  nx <- length(x0)


  # 2. moment coordinates (matching lav_implied_to_vec)
  co <- lav_sem_js_deriv_coords(lavmodel)
  m_total <- co$m_total

  # 3. stage 1: directed parameters, pass by pass
  dx <- matrix(0, nrow = nx, ncol = m_total)
  d_slope_blocks <- vector("list", nblocks)
  if (length(free_directed_idx) > 0L) {
    for (b in seq_len(nblocks)) {
      bt <- js_trace[[b]]
      if (is.null(bt)) {
        next
      }
      s_mat <- lavh1$implied$cov[[b]]
      m_vec <- if (lavmodel@meanstructure) lavh1$implied$mean[[b]] else NULL
      nvar <- length(bt$ov_names)
      # dtheta: nvar x m_total (recorded local vech coordinates)
      dtheta <- matrix(0, nrow = nvar, ncol = m_total)
      if (!is.null(bt$dtheta)) {
        dtheta[, co$cov_vech[[b]]] <- bt$dtheta
      }
      # free-loading map (nvar x nfac): the free-parameter index of each
      # loading (0 = fixed)
      lam_free_map <- matrix(0L, nrow = nvar, ncol = length(bt$lv_names))
      if (length(bt$mm_idx) > 0L && length(bt$lv_names) > 0L) {
        lam_free_map[cbind(
          match(lavpartable$rhs[bt$mm_idx], bt$ov_names),
          match(lavpartable$lhs[bt$mm_idx], bt$lv_names)
        )] <- lavpartable$free[bt$mm_idx]
      }
      d_slope_blocks[[b]] <- vector("list", length(eqs_rec[[b]]))
      for (pass_id in names(bt$passes)) {
        # the map rebuilds lambda_mat once per pass: the loading values
        # (and hence their derivatives) are frozen at the pass start,
        # even when an equation earlier in the same pass re-estimates a
        # loading that a later equation uses as proxy scale
        dx_lam <- dx
        for (j in bt$passes[[pass_id]]) {
          tr <- eqs_rec[[b]][[j]]$js_trace[[pass_id]]
          if (is.null(tr)) {
            next
          }
          out <- lav_sem_js_deriv_eq(
            tr = tr, eq = eqs_rec[[b]][[j]], s_mat = s_mat, m_vec = m_vec,
            theta = bt$theta, dtheta = dtheta, csmall = bt$csmall,
            lam_free_map = lam_free_map, dx = dx, dx_lam = dx_lam,
            co_b = co$block[[b]],
            n = lavsamplestats@nobs[[b]],
            meanstructure = lavmodel@meanstructure,
            lavpartable = lavpartable
          )
          dx <- out$dx
          if (!is.null(out$d_slope_block)) {
            d_slope_blocks[[b]][[j]] <- out$d_slope_block
          }
        }
      }
    }

    # pooled (system) solve for shared directed coefficients
    dx <- lav_sem_js_deriv_pool(
      eqs_rec = eqs_rec, d_slope_blocks = d_slope_blocks, dx = dx,
      x0 = x0, lavmodel = lavmodel, lavpartable = lavpartable
    )
  }

  # 4. stage 2: undirected parameters
  lavmodel0 <- lav_model_set_parameters(lavmodel = lavmodel, x = x0)
  if (length(free_undirected_idx) > 0L && js_varcov_method != "NONE") {
    delta_list0 <- lav_sem_miiv_delta(lavmodel0)
    for (b in seq_len(nblocks)) {
      fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
        lavpartable$free > 0L & pt_block == b])
      fu_b <- fu_b[fu_b %in% free_undirected_idx]
      if (length(fu_b) == 0L) {
        next
      }
      # replicate the point solve of *this* block: when undirected
      # parameters are shared across blocks (cross-group equality), later
      # blocks overwrite the shared entries, so x0 does not contain this
      # block's own solution (needed for the residual and the weight
      # matrix of the derivative)
      theta2_b <- as.numeric(lav_sem_miiv_varcov_block(
        b = b, fu = fu_b, delta_b = delta_list0[[b]],
        implied = lavh1$implied, lavmodel = lavmodel, x = x0,
        iv_varcov_method = js_varcov_method, return_h = FALSE
      ))
      dx[fu_b, ] <- lav_sem_js_deriv_varcov(
        b = b, fu = fu_b, delta_b = delta_list0[[b]],
        lavmodel = lavmodel, lavmodel0 = lavmodel0,
        lavh1 = lavh1, x0 = x0,
        js_varcov_method = js_varcov_method,
        free_directed_idx = free_directed_idx,
        dx = dx, co_b = co$block[[b]],
        theta2_use = theta2_b
      )
    }
  } else if (length(free_undirected_idx) > 0L) {
    dx[free_undirected_idx, ] <- as.numeric(NA)
  }

  # 5. joint GLS mean solve
  if (lavmodel@meanstructure && identical(js_mean_structure, "wls")) {
    free_mean_idx <- unique(lavpartable$free[lavpartable$op == "~1" &
      lavpartable$free > 0L & !duplicated(lavpartable$free)])
    if (length(free_mean_idx) > 0L) {
      dx[free_mean_idx, ] <- lav_sem_js_deriv_mean_wls(
        lavmodel = lavmodel, lavmodel0 = lavmodel0,
        lavsamplestats = lavsamplestats, lavh1 = lavh1,
        x0 = x0, free_mean_idx = free_mean_idx, dx = dx, co = co
      )
    }
  }

  dx
}


# moment coordinates per block, matching the stacking order of
# lav_implied_to_vec / lav_samp_wls_obs for the JS model class (continuous,
# single level, not conditional.x): per group, first the group weight (if
# group.w.free), then the means (if meanstructure), then vech(cov)
lav_sem_js_deriv_coords <- function(lavmodel) {
  nblocks <- lavmodel@nblocks
  block <- vector("list", nblocks)
  cov_vech <- vector("list", nblocks)
  offset <- 0L
  for (b in seq_len(nblocks)) {
    nvar <- lavmodel@nvar[b]
    pstar <- nvar * (nvar + 1L) / 2L
    if (lavmodel@group.w.free) {
      offset <- offset + 1L
    }
    mean_idx <- integer(0L)
    if (lavmodel@meanstructure) {
      mean_idx <- offset + seq_len(nvar)
      offset <- offset + nvar
    }
    cv <- offset + seq_len(pstar)
    offset <- offset + pstar
    # (i, j) covariance entry -> global coordinate
    cm <- matrix(0L, nrow = nvar, ncol = nvar)
    cm[lower.tri(cm, diag = TRUE)] <- cv
    cm <- cm + t(cm) - diag(diag(cm))
    block[[b]] <- list(mean = mean_idx, cov = cm, cov_vech = cv,
                       nvar = nvar, pstar = pstar)
    cov_vech[[b]] <- cv
  }
  list(block = block, cov_vech = cov_vech, m_total = offset)
}


# derivative of the (sum-to-one normalized) alphamax aggregation weights
# with respect to the sample moments, given the derivative columns of the
# proxy submatrix (via the coordinate map) and of the theta values. Solves
# the bordered system of the generalized eigenproblem
#   diag(theta) w0 = lambda S w0,  w0' S w0 = 1
# for all moment directions at once, then chains the normalization
# w = w0 / sum(w0).
lav_sem_js_deriv_weights <- function(s_set = NULL, theta_set = NULL,
                                     set = NULL, co_b = NULL,
                                     dtheta_set = NULL, m_total = NULL) {
  k <- ncol(s_set)
  # recompute the eigen solution (same code path as lav_sem_js_weights)
  r_mat <- chol(s_set)
  r_inv <- backsolve(r_mat, diag(k))
  m_sym <- crossprod(r_inv, theta_set * r_inv)
  ev <- eigen((m_sym + t(m_sym)) / 2, symmetric = TRUE)
  lambda <- ev$values[k]
  w0 <- drop(r_inv %*% ev$vectors[, k]) # B-normalized: w0' S w0 = 1
  s0 <- sum(w0)

  # dS w0 (k x m): only the coordinates of the set x set entries
  d_sw0 <- matrix(0, nrow = k, ncol = m_total)
  for (a in seq_len(k)) {
    for (b in seq_len(a)) {
      cc <- co_b$cov[set[a], set[b]]
      d_sw0[a, cc] <- d_sw0[a, cc] + w0[b]
      if (a != b) {
        d_sw0[b, cc] <- d_sw0[b, cc] + w0[a]
      }
    }
  }
  # dTheta w0 (k x m)
  d_tw0 <- dtheta_set * w0

  # dlambda = w0' dTheta w0 - lambda * w0' dS w0
  dlambda <- colSums(w0 * d_tw0) - lambda * colSums(w0 * d_sw0)

  # bordered solve for dw0
  sw0 <- drop(s_set %*% w0)
  rhs <- rbind(
    -d_tw0 + outer(sw0, dlambda) + lambda * d_sw0,
    -0.5 * colSums(w0 * d_sw0)
  )
  bmat <- rbind(
    cbind(diag(theta_set, nrow = k) - lambda * s_set, sw0),
    c(sw0, 0)
  )
  sol <- solve(bmat, rhs)
  dw0 <- sol[seq_len(k), , drop = FALSE]

  # normalization w = w0 / sum(w0)
  dw0 / s0 - outer(w0, colSums(dw0)) / s0^2
}


# derivative of one recorded equation solve (one pass): updates the rows of
# dx for the free slopes and (if any) the free intercept, and returns the
# derivative version of the slope block for the pooled solve
lav_sem_js_deriv_eq <- function(tr = NULL, eq = NULL, s_mat = NULL,
                                m_vec = NULL, theta = NULL, dtheta = NULL,
                                csmall = 1, lam_free_map = NULL, dx = NULL,
                                dx_lam = NULL, co_b = NULL, n = NULL,
                                meanstructure = FALSE,
                                lavpartable = NULL) {
  m_total <- ncol(dx)

  # intercept-only equation: x[free_int] = mean(y)
  if (isTRUE(tr$int_only)) {
    if (meanstructure && length(eq$ptint) == 1L) {
      free_int_idx <- lavpartable$free[eq$ptint]
      if (free_int_idx > 0L) {
        dx[free_int_idx, ] <- 0
        dx[free_int_idx, co_b$mean[tr$y_idx]] <- 1
      }
    }
    return(list(dx = dx, d_slope_block = NULL))
  }

  u_mat <- tr$u_mat
  is_lv <- tr$is_lv
  kappa <- tr$kappa
  nrhs <- nrow(u_mat)
  y_idx <- tr$y_idx

  # derivative of the proxy weights (rows of u_mat)
  dw_list <- vector("list", nrhs)
  dkap <- matrix(0, nrow = nrhs, ncol = m_total)
  derr <- matrix(0, nrow = nrhs, ncol = m_total)
  for (i in which(is_lv)) {
    set <- tr$sets[[i]]
    w <- tr$w[[i]]
    f <- tr$fac[i]
    k <- length(set)
    if (k > 1L && !tr$w_fb[i]) {
      dw <- lav_sem_js_deriv_weights(
        s_set = s_mat[set, set, drop = FALSE], theta_set = theta[set],
        set = set, co_b = co_b,
        dtheta_set = dtheta[set, , drop = FALSE], m_total = m_total
      )
    } else {
      dw <- matrix(0, nrow = k, ncol = m_total)
    }
    dw_list[[i]] <- dw
    # kappa_i = sum(w * lambda[set, f])
    lam_vals <- tr$lam[[i]]
    dkap[i, ] <- colSums(dw * lam_vals)
    for (a in seq_len(k)) {
      fidx <- lam_free_map[set[a], f]
      if (fidx > 0L) {
        dkap[i, ] <- dkap[i, ] + w[a] * dx_lam[fidx, ]
      }
    }
    # errvar_i = sum(w^2 * theta[set])
    derr[i, ] <- colSums((2 * w * theta[set]) * dw) +
      colSums(w^2 * dtheta[set, , drop = FALSE])
  }

  # supports (nonzero columns of u_mat)
  supp <- which(colSums(abs(u_mat)) > 0)
  nvar <- ncol(u_mat)

  # d(s_zz) = dU S U' + U dS U' + U S dU'    (nrhs^2 x m, vec order)
  # d(s_zy) = dU S[, y] + U dS[, y]          (nrhs x m)
  d_szz <- matrix(0, nrow = nrhs * nrhs, ncol = m_total)
  d_szy <- matrix(0, nrow = nrhs, ncol = m_total)
  # -- U dS U' and U dS[, y]: loop over the coordinates that touch supp/y
  cols_ks <- sort(unique(c(supp, y_idx)))
  seen <- matrix(FALSE, nrow = nvar, ncol = nvar)
  for (k1 in cols_ks) {
    for (l1 in seq_len(nvar)) {
      kk <- max(k1, l1)
      ll <- min(k1, l1)
      if (seen[kk, ll]) {
        next
      }
      seen[kk, ll] <- TRUE
      cc <- co_b$cov[kk, ll]
      u_k <- u_mat[, kk]
      u_l <- u_mat[, ll]
      # U E_kl U' + U E_lk U' (single term when kk == ll)
      m1 <- outer(u_k, u_l)
      if (kk != ll) {
        m1 <- m1 + outer(u_l, u_k)
      }
      if (any(m1 != 0)) {
        d_szz[, cc] <- d_szz[, cc] + as.vector(m1)
      }
      # U dS[, y]
      if (ll == y_idx || kk == y_idx) {
        if (kk == ll) {
          d_szy[, cc] <- d_szy[, cc] + u_k
        } else if (ll == y_idx) {
          d_szy[, cc] <- d_szy[, cc] + u_k
        } else {
          d_szy[, cc] <- d_szy[, cc] + u_l
        }
      }
    }
  }
  # -- dU parts
  s_u <- s_mat %*% t(u_mat) # nvar x nrhs
  s_y <- s_mat[, y_idx]
  for (i in which(is_lv)) {
    dw <- dw_list[[i]]
    if (is.null(dw) || all(dw == 0)) {
      next
    }
    set <- tr$sets[[i]]
    # row i of dU S U': t(dw) %*% s_u[set, ] -> (m x nrhs)
    t1 <- crossprod(dw, s_u[set, , drop = FALSE]) # m x nrhs
    for (j2 in seq_len(nrhs)) {
      pos_ij <- (j2 - 1L) * nrhs + i
      pos_ji <- (i - 1L) * nrhs + j2
      d_szz[pos_ij, ] <- d_szz[pos_ij, ] + t1[, j2]
      if (j2 != i) {
        d_szz[pos_ji, ] <- d_szz[pos_ji, ] + t1[, j2]
      }
    }
    # note: the diagonal (i, i) gets dU S U' + U S dU' = 2 * t1[, i]; the
    # loop above adds t1[, i] once for pos_ij == pos_ji; add the transpose
    # contribution for the diagonal explicitly
    pos_ii <- (i - 1L) * nrhs + i
    d_szz[pos_ii, ] <- d_szz[pos_ii, ] + t1[, i]
    # dU S[, y]
    d_szy[i, ] <- d_szy[i, ] + crossprod(dw, s_y[set])
  }

  # d(G): g0 = s_zz with (lv, lv) diagonal minus csmall * errvar; then rows
  # scaled by 1 / kappa
  d_g <- d_szz
  for (i in which(is_lv)) {
    pos_ii <- (i - 1L) * nrhs + i
    d_g[pos_ii, ] <- d_g[pos_ii, ] - csmall * derr[i, ]
  }
  # reconstruct the unscaled g0 from the recorded (scaled) g_mat
  g_scaled <- tr$g_mat
  g0 <- g_scaled * kappa # rows scaled back
  for (i in which(is_lv)) {
    rows_i <- (seq_len(nrhs) - 1L) * nrhs + i
    d_g[rows_i, ] <- d_g[rows_i, ] / kappa[i] -
      outer(g0[i, ] / kappa[i]^2, dkap[i, ])
  }

  # d(amat), d(bvec)
  s_zz <- tr$s_zz
  szz_inv <- solve((s_zz + t(s_zz)) / 2)
  p_mat <- g_scaled %*% szz_inv
  u2 <- drop(szz_inv %*% (u_mat %*% s_y))
  tperm <- as.vector(t(matrix(seq_len(nrhs * nrhs), nrow = nrhs)))
  qq <- kronecker(p_mat, diag(nrhs)) %*% d_g
  damat <- qq + qq[tperm, , drop = FALSE] -
    kronecker(p_mat, p_mat) %*% d_szz
  dbvec <- kronecker(matrix(u2, nrow = 1L), diag(nrhs)) %*% d_g -
    kronecker(matrix(u2, nrow = 1L), p_mat) %*% d_szz +
    p_mat %*% d_szy

  # free/fixed partition
  fix_idx <- tr$fix_idx
  b_full <- tr$b_full
  if (length(fix_idx) > 0L) {
    rows_free <- setdiff(seq_len(nrhs), fix_idx)
  } else {
    rows_free <- seq_len(nrhs)
  }
  nf <- length(rows_free)
  vpos <- function(ii, jj) (jj - 1L) * nrhs + ii
  if (nf > 0L) {
    free_pairs <- as.matrix(expand.grid(rows_free, rows_free))
    damat_free <- damat[vpos(free_pairs[, 1L], free_pairs[, 2L]), ,
                        drop = FALSE]
    dbvec_free <- dbvec[rows_free, , drop = FALSE]
    if (length(fix_idx) > 0L) {
      for (jf in fix_idx) {
        if (b_full[jf] != 0) {
          dbvec_free <- dbvec_free -
            b_full[jf] * damat[vpos(rows_free, rep(jf, nf)), , drop = FALSE]
        }
      }
    }
    # d(slopes) = A_free^{-1} (dbvec_free - damat_free %*% slopes)
    amat <- tr$amat
    amat_free <- amat[rows_free, rows_free, drop = FALSE]
    slopes <- b_full[rows_free]
    contr <- matrix(0, nrow = nf, ncol = m_total)
    for (j2 in seq_len(nf)) {
      rows_j <- seq_len(nf) + (j2 - 1L) * nf
      contr <- contr + slopes[j2] * damat_free[rows_j, , drop = FALSE]
    }
    dslopes <- solve(amat_free, dbvec_free - contr)
    gcol_free <- tr$gcol_free
    dx[gcol_free, ] <- dslopes
  } else {
    dslopes <- matrix(0, nrow = 0L, ncol = m_total)
    damat_free <- matrix(0, nrow = 0L, ncol = m_total)
    dbvec_free <- matrix(0, nrow = 0L, ncol = m_total)
    gcol_free <- integer(0L)
  }

  # intercept and mean-part derivatives
  d_x_bar <- NULL
  d_y_bar <- NULL
  if (meanstructure) {
    # z_mean = U m; x_bar = z_mean (obs) or z_mean / kappa (lv)
    d_z <- matrix(0, nrow = nrhs, ncol = m_total)
    for (v in supp) {
      d_z[, co_b$mean[v]] <- d_z[, co_b$mean[v]] + u_mat[, v]
    }
    for (i in which(is_lv)) {
      dw <- dw_list[[i]]
      if (!is.null(dw) && !all(dw == 0)) {
        d_z[i, ] <- d_z[i, ] + crossprod(dw, m_vec[tr$sets[[i]]])
      }
    }
    x_bar <- tr$x_bar
    d_x_bar <- d_z
    for (i in which(is_lv)) {
      d_x_bar[i, ] <- d_z[i, ] / kappa[i] - x_bar[i] * dkap[i, ] / kappa[i]
    }
    d_y_bar <- matrix(0, nrow = 1L, ncol = m_total)
    d_y_bar[1L, co_b$mean[y_idx]] <- 1

    # db_full: free rows from dslopes, fixed rows zero
    db_full <- matrix(0, nrow = nrhs, ncol = m_total)
    if (nf > 0L) {
      db_full[rows_free, ] <- dslopes
    }
    dbeta0 <- d_y_bar -
      matrix(colSums(db_full * x_bar) + colSums(d_x_bar * b_full),
             nrow = 1L)
    if (length(eq$ptint) == 1L) {
      free_int_idx <- lavpartable$free[eq$ptint]
      if (free_int_idx > 0L) {
        dx[free_int_idx, ] <- dbeta0
      }
    }
  }

  d_slope_block <- NULL
  if (length(gcol_free) > 0L) {
    d_slope_block <- list(
      damat = (n - 1) * damat_free, dbvec = (n - 1) * dbvec_free,
      gcol = gcol_free, nf = nf,
      d_x_bar = d_x_bar, d_y_bar = d_y_bar
    )
  }

  list(dx = dx, d_slope_block = d_slope_block)
}


# derivative of the pooled directed solve (shared slope columns and/or
# linear equality constraints among the directed slopes) plus the pooled
# intercept recomputation -- mirrors lav_sem_miiv_apply_directed_pool()
lav_sem_js_deriv_pool <- function(eqs_rec = NULL, d_slope_blocks = NULL,
                                  dx = NULL, x0 = NULL, lavmodel = NULL,
                                  lavpartable = NULL) {
  eqs_all <- do.call(c, eqs_rec)
  d_all <- do.call(c, d_slope_blocks)
  slope_blocks <- lapply(eqs_all, `[[`, "slope_block")
  keep <- !vapply(slope_blocks, is.null, logical(1L))
  slope_blocks <- slope_blocks[keep]
  d_sb <- d_all[keep]
  eqs_kept <- eqs_all[keep]
  if (length(slope_blocks) == 0L) {
    return(dx)
  }
  all_gcol <- unlist(lapply(slope_blocks, `[[`, "gcol"))
  free_slope_idx <- sort(unique(all_gcol))
  con <- lav_sem_miiv_linear_con(lavmodel, free_slope_idx)
  all_int <- unlist(lapply(eqs_kept, function(e) {
    if (!is.null(e$ptint) && length(e$ptint) == 1L) {
      lavpartable$free[e$ptint]
    } else {
      integer(0L)
    }
  }))
  shared_int <- anyDuplicated(all_int[all_int > 0L]) > 0L
  if (anyDuplicated(all_gcol) == 0L && is.null(con) && !shared_int) {
    return(dx) # nothing pooled; the per-equation rows are final
  }

  m_total <- ncol(dx)
  npool <- length(free_slope_idx)

  # accumulate the pooled system and its derivative
  dmat <- matrix(0, nrow = npool, ncol = npool)
  dvec <- numeric(npool)
  d_dmat <- matrix(0, nrow = npool * npool, ncol = m_total)
  d_dvec <- matrix(0, nrow = npool, ncol = m_total)
  for (ib in seq_along(slope_blocks)) {
    bl <- slope_blocks[[ib]]
    dbl <- d_sb[[ib]]
    p <- match(bl$gcol, free_slope_idx)
    nf <- length(p)
    if (anyDuplicated(p) > 0L) {
      # a slope label repeated WITHIN one equation: collapse the
      # duplicated rows/columns (C'AC / C'b), exactly as in
      # lav_sem_miiv_pool_directed(); the derivative rows collapse the
      # same way, grouped by the global (row, col) cell
      a_c <- rowsum(bl$amat, group = p)
      a_c <- t(rowsum(t(a_c), group = p))
      b_c <- drop(rowsum(matrix(bl$bvec, ncol = 1L), group = p))
      pu <- sort(unique(p))
      dmat[pu, pu] <- dmat[pu, pu] + a_c
      dvec[pu] <- dvec[pu] + b_c
      grp <- as.vector(outer(p, p, function(i, j) (j - 1L) * npool + i))
      dd <- rowsum(dbl$damat, group = grp)
      gu <- sort(unique(grp))
      d_dmat[gu, ] <- d_dmat[gu, ] + dd
      d_dvec[pu, ] <- d_dvec[pu, ] + rowsum(dbl$dbvec, group = p)
    } else {
      dmat[p, p] <- dmat[p, p] + bl$amat
      dvec[p] <- dvec[p] + bl$bvec
      for (j2 in seq_len(nf)) {
        rows_j <- seq_len(nf) + (j2 - 1L) * nf
        d_dmat[(p[j2] - 1L) * npool + p, ] <-
          d_dmat[(p[j2] - 1L) * npool + p, ] +
          dbl$damat[rows_j, , drop = FALSE]
      }
      d_dvec[p, ] <- d_dvec[p, ] + dbl$dbvec
    }
  }

  # pooled point solution: the map already wrote it into x0 (and nothing
  # after the pooled solve touches the directed slopes)
  sol <- x0[free_slope_idx]

  # d(sol): contract d_dmat with sol
  contr <- matrix(0, nrow = npool, ncol = m_total)
  for (j2 in seq_len(npool)) {
    rows_j <- seq_len(npool) + (j2 - 1L) * npool
    contr <- contr + sol[j2] * d_dmat[rows_j, , drop = FALSE]
  }
  rhs <- d_dvec - contr
  if (!is.null(con) && nrow(con$jac) > 0L) {
    ncon <- nrow(con$jac)
    kkt <- rbind(
      cbind(dmat, t(con$jac)),
      cbind(con$jac, matrix(0, ncon, ncon))
    )
    dsol <- solve(kkt, rbind(rhs, matrix(0, ncon, m_total)))
    dsol <- dsol[seq_len(npool), , drop = FALSE]
  } else {
    dsol <- solve(dmat, rhs)
  }

  # write the pooled slope derivatives
  for (ib in seq_along(slope_blocks)) {
    bl <- slope_blocks[[ib]]
    p <- match(bl$gcol, free_slope_idx)
    dx[bl$gcol, ] <- dsol[p, , drop = FALSE]
  }

  # recompute the intercepts from the pooled slopes (nobs-weighted average
  # for shared intercepts), as in lav_sem_miiv_apply_directed_pool()
  if (lavmodel@meanstructure) {
    int_idx <- integer(0L)
    int_wgt <- numeric(0L)
    int_dval <- vector("list", 0L)
    for (ib in seq_along(eqs_kept)) {
      eq <- eqs_kept[[ib]]
      sb <- eq$slope_block
      dbl <- d_sb[[ib]]
      free_int_idx <- lavpartable$free[eq$ptint]
      if (length(free_int_idx) != 1L || free_int_idx <= 0L) {
        next
      }
      eq_free_idx <- lavpartable$free[eq$pt]
      b_full <- lavpartable$ustart[eq$pt]
      b_full[is.na(b_full)] <- 0
      free_pos <- which(eq_free_idx > 0L)
      b_full[free_pos] <- x0[eq_free_idx[free_pos]]
      db_full <- matrix(0, nrow = length(b_full), ncol = m_total)
      db_full[free_pos, ] <- dx[eq_free_idx[free_pos], , drop = FALSE]
      dbeta0 <- dbl$d_y_bar -
        matrix(colSums(db_full * sb$x_bar) +
               colSums(dbl$d_x_bar * b_full), nrow = 1L)
      int_idx <- c(int_idx, free_int_idx)
      int_wgt <- c(int_wgt, if (!is.null(eq$nobs)) eq$nobs else 1)
      int_dval <- c(int_dval, list(dbeta0))
    }
    if (length(int_idx) > 0L) {
      for (ii in unique(int_idx)) {
        sel <- which(int_idx == ii)
        wsum <- sum(int_wgt[sel])
        acc <- matrix(0, nrow = 1L, ncol = m_total)
        for (s2 in sel) {
          acc <- acc + int_wgt[s2] * int_dval[[s2]]
        }
        dx[ii, ] <- acc / wsum
      }
    }
  }

  dx
}


# GLIST entry map for one block: for every free parameter, its entries
# (row, col) in the lambda / beta / psi / theta / alpha / nu matrices of
# that block (symmetric matrices carry both triangles, which makes the
# accumulated E-matrix derivative automatically symmetric)
lav_sem_js_deriv_glist_map <- function(lavmodel, b) {
  mm_in_group <- seq_len(lavmodel@nmat[b]) + cumsum(c(0L, lavmodel@nmat))[b]
  unco2free <- NULL
  if (lavmodel@ceq.simple.only) {
    unco2free <- max.col(lavmodel@ceq.simple.K)
  }
  out <- list()
  for (mm in mm_in_group) {
    mname <- names(lavmodel@GLIST)[mm]
    midx <- lavmodel@m.free.idx[[mm]]
    if (length(midx) == 0L) {
      next
    }
    if (lavmodel@ceq.simple.only) {
      fidx <- unco2free[lavmodel@x.unco.idx[[mm]]]
    } else {
      fidx <- lavmodel@x.free.idx[[mm]]
    }
    rc <- arrayInd(midx, dim(lavmodel@GLIST[[mm]]))
    out[[mname]] <- cbind(row = rc[, 1L], col = rc[, 2L], free = fidx)
  }
  out
}


# dA_j = d(Lambda (I - B)^{-1}) / d x_j for the directed free parameters of
# one block: a list (indexed by free parameter number) of p x q matrices;
# only parameters with lambda/beta entries in this block get an entry
lav_sem_js_deriv_da <- function(gmap = NULL, a_mat = NULL, ib_inv = NULL,
                                free_idx = NULL) {
  p <- nrow(a_mat)
  q <- ncol(a_mat)
  da <- vector("list", 0L)
  add_da <- function(da, j, contrib) {
    key <- as.character(j)
    if (is.null(da[[key]])) {
      da[[key]] <- contrib
    } else {
      da[[key]] <- da[[key]] + contrib
    }
    da
  }
  if (!is.null(gmap$lambda)) {
    for (e in seq_len(nrow(gmap$lambda))) {
      j <- gmap$lambda[e, "free"]
      if (!j %in% free_idx) {
        next
      }
      r <- gmap$lambda[e, "row"]
      cc <- gmap$lambda[e, "col"]
      contrib <- matrix(0, nrow = p, ncol = q)
      contrib[r, ] <- ib_inv[cc, ]
      da <- add_da(da, j, contrib)
    }
  }
  if (!is.null(gmap$beta)) {
    for (e in seq_len(nrow(gmap$beta))) {
      j <- gmap$beta[e, "free"]
      if (!j %in% free_idx) {
        next
      }
      r <- gmap$beta[e, "row"]
      cc <- gmap$beta[e, "col"]
      da <- add_da(da, j, outer(a_mat[, r], ib_inv[cc, ]))
    }
  }
  da
}


# derivative of the second-stage (variance/covariance) solve for one block
# with respect to the sample moments: chains the direct svec dependence,
# the dependence on the directed estimates (through Delta2 and through the
# weight matrix), and -- for RLS -- solves the fixed-point sensitivity via
# the implicit function theorem
lav_sem_js_deriv_varcov <- function(b = 1L, fu = NULL, delta_b = NULL,
                                    lavmodel = NULL, lavmodel0 = NULL,
                                    lavh1 = NULL, x0 = NULL,
                                    js_varcov_method = "RLS",
                                    free_directed_idx = NULL,
                                    dx = NULL, co_b = NULL,
                                    theta2_use = NULL) {
  m_total <- ncol(dx)
  meanstructure <- lavmodel@meanstructure
  nvar <- co_b$nvar
  pstar <- co_b$pstar
  cov_rows <- if (meanstructure) nvar + seq_len(pstar) else seq_len(pstar)

  d2c <- delta_b[cov_rows, fu, drop = FALSE]
  n_u <- length(fu)
  # the coefficient vector of *this* solve step: the final estimates,
  # except for the first (ULS) step of 2RLS (theta2_use passed in)
  theta2 <- if (is.null(theta2_use)) x0[fu] else theta2_use
  svec_c <- lav_mat_vech(lavh1$implied$cov[[b]])
  r_vec <- drop(svec_c - d2c %*% theta2)
  r_mat <- lav_mat_vech_rev(r_vec)

  # vech weights (0.5 on the diagonal)
  wvech <- rep(1, pstar)
  wvech[lav_mat_diagh_idx(nvar)] <- 0.5

  # the weight matrix S of the final linearization, and the Delta whose
  # directed columns give dSigma/dx1 for the weight-matrix chain (for 2RLS
  # the weight Sigma is evaluated at the first-step -- ULS -- estimates)
  uls_like <- js_varcov_method == "ULS"
  theta2_uls <- NULL
  delta_w_b <- delta_b
  if (js_varcov_method == "GLS") {
    s_w <- lavh1$implied$cov[[b]]
  } else if (js_varcov_method == "RLS") {
    # at the fixed point, the weight Sigma is the model-implied Sigma at
    # *this block's* solution (which may differ from x0 for undirected
    # parameters shared across blocks)
    x_b <- x0
    x_b[fu] <- theta2
    lavmodel_w <- lav_model_set_parameters(lavmodel = lavmodel0, x = x_b)
    s_w <- lav_model_sigma(lavmodel = lavmodel_w, extra = FALSE)[[b]]
    # a non-PD weight matrix means the point estimation fell back;
    # signal the caller to use the numerical Jacobian
    if (inherits(try(chol(s_w), silent = TRUE), "try-error")) {
      lav_msg_stop(gettext("[JS] non-PD weight matrix in stage 2."))
    }
    # the implicit-function-theorem derivative below is only valid at the
    # fixed point; verify with one more map application (the point RLS
    # loop may have stopped -- with a warning -- before converging)
    svec_full <- if (meanstructure) {
      c(lavh1$implied$mean[[b]], svec_c)
    } else {
      svec_c
    }
    theta2_chk <- as.numeric(lav_utils_wls_linearization(
      delta = delta_b[, fu, drop = FALSE], s = s_w,
      meanstructure = meanstructure, categorical = FALSE,
      svec = svec_full, return_h = FALSE
    ))
    if (max(abs(theta2_chk - theta2)) > 1e-5 * (1 + max(abs(theta2)))) {
      lav_msg_stop(gettext(
        "[JS] stage-2 RLS did not reach its fixed point."))
    }
    delta_w_b <- lav_sem_miiv_delta(lavmodel_w)[[b]]
  } else if (js_varcov_method == "2RLS") {
    svec_full <- if (meanstructure) {
      c(lavh1$implied$mean[[b]], svec_c)
    } else {
      svec_c
    }
    theta2_uls <- as.numeric(lav_utils_wls_linearization(
      delta = delta_b[, fu, drop = FALSE],
      s = diag(1, nrow = nvar),
      meanstructure = meanstructure, categorical = FALSE,
      svec = svec_full, return_h = FALSE
    ))
    x_uls <- x0
    x_uls[fu] <- theta2_uls
    lavmodel_uls <- lav_model_set_parameters(lavmodel = lavmodel0,
                                             x = x_uls)
    s_w <- lav_model_sigma(lavmodel = lavmodel_uls, extra = FALSE)[[b]]
    if (inherits(try(chol(s_w), silent = TRUE), "try-error")) {
      lav_msg_stop(gettext("[JS] non-PD weight matrix in stage 2."))
    }
    delta_w_b <- lav_sem_miiv_delta(lavmodel_uls)[[b]]
  }

  # W-applied Delta2 (q_cov), A, Ainv
  if (uls_like) {
    q_cov <- wvech * d2c
    si <- NULL
  } else {
    si <- lav_mat_sym_inverse(s_w)
    q_cov <- matrix(0, nrow = pstar, ncol = n_u)
    for (k in seq_len(n_u)) {
      vmat <- lav_mat_vech_rev(d2c[, k])
      q_cov[, k] <- wvech * lav_mat_vech(si %*% vmat %*% si)
    }
  }
  a_small <- crossprod(d2c, q_cov)
  a_inv <- solve((a_small + t(a_small)) / 2)

  # (a) direct svec part: H = Ainv D2' W, placed at the cov coordinates
  j2 <- matrix(0, nrow = n_u, ncol = m_total)
  j2[, co_b$cov_vech] <- a_inv %*% t(q_cov)

  # (b) chain through the directed estimates (Delta2(x1) and, for
  #     RLS/2RLS, the weight matrix Sigma(x1, theta2))
  gmap <- lav_sem_js_deriv_glist_map(lavmodel0, b)
  glist_b <- lavmodel0@GLIST[
    seq_len(lavmodel0@nmat[b]) + cumsum(c(0L, lavmodel0@nmat))[b]
  ]
  lambda0 <- glist_b[["lambda"]]
  q_lv <- ncol(lambda0)
  ib_inv <- diag(q_lv)
  if (!is.null(glist_b[["beta"]])) {
    ib_inv <- solve(diag(q_lv) - glist_b[["beta"]])
  }
  a_mat <- lambda0 %*% ib_inv
  da_list <- lav_sem_js_deriv_da(
    gmap = gmap, a_mat = a_mat, ib_inv = ib_inv,
    free_idx = free_directed_idx
  )

  # psi entries of the undirected free parameters (for dDelta2 theta2 and
  # the per-parameter E structure)
  psi_entries <- gmap$psi
  psi_by_free <- NULL
  psi_acc <- matrix(0, nrow = q_lv, ncol = q_lv)
  if (!is.null(psi_entries)) {
    keep <- psi_entries[, "free"] %in% fu
    psi_entries <- psi_entries[keep, , drop = FALSE]
    if (nrow(psi_entries) > 0L) {
      psi_by_free <- split(
        seq_len(nrow(psi_entries)),
        factor(psi_entries[, "free"], levels = fu)
      )
      # weighted by the coefficient vector of *this* solve step
      psi_acc[cbind(psi_entries[, "row"], psi_entries[, "col"])] <-
        psi_acc[cbind(psi_entries[, "row"], psi_entries[, "col"])] +
        theta2[match(psi_entries[, "free"], fu)]
    }
  }

  # helper: apply the derivative of the W-weighting to a direction V
  # (symmetric p x p): returns w * vech(-(Si V Q + Q V Si)) with
  # Q = Si R Si; for ULS the weight is constant (returns NULL)
  qr_mat <- NULL
  if (!uls_like) {
    qr_mat <- si %*% r_mat %*% si
  }
  dw_apply <- function(v_mat) {
    z <- si %*% v_mat %*% qr_mat
    wvech * lav_mat_vech(-(z + t(z)))
  }

  # M: derivative of the stage-2 map with respect to the directed free
  # parameters (n_u x nx), nonzero only for the columns with dA entries
  fd_keys <- names(da_list)
  m_dir <- matrix(0, nrow = n_u, ncol = length(x0))
  if (length(fd_keys) > 0L) {
    aq <- if (uls_like) {
      crossprod(a_mat, r_mat)
    } else {
      crossprod(a_mat, qr_mat)
    } # q x p : A' Q
    psi_a <- psi_acc %*% t(a_mat) # q x p
    for (key in fd_keys) {
      j <- as.integer(key)
      da_j <- da_list[[key]]
      # (i) dD2' W r, for every undirected parameter k: tr(S_ab C_j) with
      #     C_j = A' Q dA_j (only psi columns; theta columns are constant)
      c_j <- aq %*% da_j # q x q
      t1 <- numeric(n_u)
      if (!is.null(psi_by_free)) {
        for (k in seq_len(n_u)) {
          ee <- psi_by_free[[k]]
          if (length(ee) == 0L) {
            next
          }
          t1[k] <- sum(c_j[cbind(psi_entries[ee, "col"],
                                 psi_entries[ee, "row"])])
        }
      }
      # (ii) D2' W (dD2 theta2), with dD2 theta2 = vech(dA_j PsiF A' +
      #      A PsiF dA_j')
      v_th <- da_j %*% psi_a + t(da_j %*% psi_a)
      if (uls_like) {
        wv <- wvech * lav_mat_vech(v_th)
      } else {
        z2 <- si %*% v_th %*% si
        wv <- wvech * lav_mat_vech(z2)
      }
      t2 <- drop(crossprod(d2c, wv))
      m_dir[, j] <- drop(a_inv %*% (t1 - t2))
      # (iii) the weight matrix depends on Sigma(x1, .) for RLS/2RLS
      if (js_varcov_method %in% c("RLS", "2RLS")) {
        # dSigma / dx_j = full-Delta directed column j (cov rows),
        # evaluated at the parameter values behind the weight matrix
        v_j <- lav_mat_vech_rev(delta_w_b[cov_rows, j])
        m_dir[, j] <- m_dir[, j] +
          drop(a_inv %*% crossprod(d2c, dw_apply(v_j)))
      }
    }
  }
  j2 <- j2 + m_dir %*% dx

  # (c) weight-matrix dependence on the *data* moments (GLS only)
  if (js_varcov_method == "GLS") {
    t_s <- matrix(0, nrow = n_u, ncol = pstar)
    for (cc in seq_len(pstar)) {
      e_v <- lav_mat_vech_rev(replace(numeric(pstar), cc, 1))
      t_s[, cc] <- drop(a_inv %*% crossprod(d2c, dw_apply(e_v)))
    }
    j2[, co_b$cov_vech] <- j2[, co_b$cov_vech] + t_s
  }

  # (d) fixed-point / two-step corrections
  if (js_varcov_method %in% c("RLS", "2RLS")) {
    # K: sensitivity of the map to theta2 through the weight matrix
    k_mat <- matrix(0, nrow = n_u, ncol = n_u)
    for (k in seq_len(n_u)) {
      v_k <- lav_mat_vech_rev(d2c[, k])
      k_mat[, k] <- drop(a_inv %*% crossprod(d2c, dw_apply(v_k)))
    }
    if (js_varcov_method == "RLS") {
      # implicit function theorem at the fixed point
      j2 <- solve(diag(n_u) - k_mat, j2)
    } else {
      # 2RLS: one ULS step, then one weighted step; chain the ULS
      # sensitivity through K (the ULS step has its own coefficient
      # vector theta2_uls)
      j2_uls <- lav_sem_js_deriv_varcov(
        b = b, fu = fu, delta_b = delta_b,
        lavmodel = lavmodel, lavmodel0 = lavmodel0,
        lavh1 = lavh1, x0 = x0,
        js_varcov_method = "ULS",
        free_directed_idx = free_directed_idx,
        dx = dx, co_b = co_b, theta2_use = theta2_uls
      )
      j2 <- j2 + k_mat %*% j2_uls
    }
  }

  j2
}


# derivative of the joint GLS mean solve (lav_sem_miiv_mean_wls) with
# respect to the sample moments, chaining through the directed estimates
# (design and residual), the undirected estimates (weight matrix) and the
# data means
lav_sem_js_deriv_mean_wls <- function(lavmodel = NULL, lavmodel0 = NULL,
                                      lavsamplestats = NULL, lavh1 = NULL,
                                      x0 = NULL, free_mean_idx = NULL,
                                      dx = NULL, co = NULL) {
  m_total <- ncol(dx)
  nfm <- length(free_mean_idx)
  nblocks <- lavmodel@nblocks

  # point quantities (as in lav_sem_miiv_mean_wls, at x0)
  implied0_x <- x0
  implied0_x[free_mean_idx] <- 0
  lavmodel00 <- lav_model_set_parameters(lavmodel0, x = implied0_x)
  implied <- lav_model_implied(lavmodel0)
  implied0 <- lav_model_implied(lavmodel00)
  delta_list <- lav_sem_miiv_delta(lavmodel0)

  amat <- matrix(0, nrow = nfm, ncol = nfm)
  bvec <- numeric(nfm)
  d_amat <- matrix(0, nrow = nfm * nfm, ncol = m_total)
  d_bvec <- matrix(0, nrow = nfm, ncol = m_total)

  for (g in seq_len(nblocks)) {
    nvar <- lavmodel@nvar[[g]]
    pstar <- nvar * (nvar + 1L) / 2L
    cov_rows <- if (lavmodel@meanstructure) {
      nvar + seq_len(pstar)
    } else {
      seq_len(pstar)
    }
    co_b <- co$block[[g]]
    n_g <- lavsamplestats@nobs[[g]]
    m_g <- delta_list[[g]][seq_len(nvar), free_mean_idx, drop = FALSE]
    mean_g <- lavh1$implied$mean[[g]]
    resid <- drop(mean_g - implied0$mean[[g]])
    w_g <- try(lav_mat_sym_inverse(implied$cov[[g]]), silent = TRUE)
    if (inherits(w_g, "try-error")) {
      w_g <- diag(nvar)
      si_flag <- FALSE
    } else {
      si_flag <- TRUE
    }
    si <- w_g
    w_g <- w_g * n_g
    aw <- w_g %*% m_g # nvar x nfm
    rw <- drop(w_g %*% resid) # nvar

    amat <- amat + crossprod(m_g, aw)
    bvec <- bvec + crossprod(m_g, rw)

    # dSigma_g (vech) as a function of the moments: all free columns
    # (only needed for the weight-matrix chain below)
    d_sig <- NULL
    if (si_flag) {
      d_sig <- delta_list[[g]][cov_rows, , drop = FALSE] %*% dx # pstar x m
    }

    # dA-matrices for the directed parameters of this block
    gmap <- lav_sem_js_deriv_glist_map(lavmodel0, g)
    glist_b <- lavmodel0@GLIST[
      seq_len(lavmodel0@nmat[g]) + cumsum(c(0L, lavmodel0@nmat))[g]
    ]
    lambda0 <- glist_b[["lambda"]]
    q_lv <- ncol(lambda0)
    ib_inv <- diag(q_lv)
    if (!is.null(glist_b[["beta"]])) {
      ib_inv <- solve(diag(q_lv) - glist_b[["beta"]])
    }
    a_mat <- lambda0 %*% ib_inv
    da_list <- lav_sem_js_deriv_da(
      gmap = gmap, a_mat = a_mat, ib_inv = ib_inv,
      free_idx = seq_along(x0) # all free parameters with lambda/beta entries
    )

    # alpha structure: entries of the free mean parameters and the
    # zeroed-free alpha vector of implied0
    alpha_entries <- gmap$alpha
    glist00_b <- lavmodel00@GLIST[
      seq_len(lavmodel00@nmat[g]) + cumsum(c(0L, lavmodel00@nmat))[g]
    ]
    alpha0 <- if (!is.null(glist00_b[["alpha"]])) {
      drop(glist00_b[["alpha"]])
    } else {
      numeric(q_lv)
    }

    # d(resid) = d(mean data) - d(mu0); d(mu0) = sum_j (dA_j alpha0) dx_j
    d_resid <- matrix(0, nrow = nvar, ncol = m_total)
    d_resid[, co_b$mean] <- diag(nvar)
    # d(m_g): only the alpha columns move: dm_g[, c] = sum_j dA_j[, a] dx_j
    d_mg <- vector("list", nfm) # each nvar x m_total (NULL = zero)
    if (length(da_list) > 0L) {
      for (key in names(da_list)) {
        j <- as.integer(key)
        da_j <- da_list[[key]]
        if (any(alpha0 != 0)) {
          d_resid <- d_resid - outer(drop(da_j %*% alpha0), dx[j, ])
        }
        if (!is.null(alpha_entries)) {
          for (e in seq_len(nrow(alpha_entries))) {
            cfree <- alpha_entries[e, "free"]
            pos <- match(cfree, free_mean_idx)
            if (is.na(pos)) {
              next
            }
            contrib <- outer(da_j[, alpha_entries[e, "row"]], dx[j, ])
            if (is.null(d_mg[[pos]])) {
              d_mg[[pos]] <- contrib
            } else {
              d_mg[[pos]] <- d_mg[[pos]] + contrib
            }
          }
        }
      }
    }

    # accumulate d_bvec and d_amat
    # (1) m' W d(resid)
    d_bvec <- d_bvec + crossprod(aw, d_resid)
    # (2) dm' W resid and dm' W m (+ transpose)
    for (cpos in seq_len(nfm)) {
      if (is.null(d_mg[[cpos]])) {
        next
      }
      d_bvec[cpos, ] <- d_bvec[cpos, ] + crossprod(d_mg[[cpos]], rw)
      for (c2 in seq_len(nfm)) {
        contrib <- crossprod(d_mg[[cpos]], aw[, c2])
        d_amat[(c2 - 1L) * nfm + cpos, ] <-
          d_amat[(c2 - 1L) * nfm + cpos, ] + t(contrib)
        d_amat[(cpos - 1L) * nfm + c2, ] <-
          d_amat[(cpos - 1L) * nfm + c2, ] + t(contrib)
      }
    }
    # (3) m' dW resid and m' dW m, with dW = -n W0 dSigma W0 (W0 = Si)
    if (si_flag) {
      # coefficient matrices over the vech coordinates
      # m' dW r: rows = nfm; -(n) * (aw0' dSig rw0 + symmetric)
      aw0 <- si %*% m_g # nvar x nfm
      rw0 <- drop(si %*% resid)
      coef_br <- matrix(0, nrow = nfm, ncol = pstar)
      coef_aa <- matrix(0, nrow = nfm * nfm, ncol = pstar)
      cc <- 0L
      for (jj in seq_len(nvar)) {
        for (ii in jj:nvar) {
          cc <- cc + 1L
          if (ii == jj) {
            coef_br[, cc] <- aw0[ii, ] * rw0[ii]
            coef_aa[, cc] <- as.vector(outer(aw0[ii, ], aw0[ii, ]))
          } else {
            coef_br[, cc] <- aw0[ii, ] * rw0[jj] + aw0[jj, ] * rw0[ii]
            coef_aa[, cc] <- as.vector(outer(aw0[ii, ], aw0[jj, ]) +
                                       outer(aw0[jj, ], aw0[ii, ]))
          }
        }
      }
      d_bvec <- d_bvec - n_g * (coef_br %*% d_sig)
      d_amat <- d_amat - n_g * (coef_aa %*% d_sig)
    }
  }

  a_inv <- lav_mat_sym_inverse(amat)
  out <- drop(a_inv %*% bvec)

  # d(out) = Ainv (d_bvec - d_amat out)
  contr <- matrix(0, nrow = nfm, ncol = m_total)
  for (j2 in seq_len(nfm)) {
    rows_j <- seq_len(nfm) + (j2 - 1L) * nfm
    contr <- contr + out[j2] * d_amat[rows_j, , drop = FALSE]
  }
  a_inv %*% (d_bvec - contr)
}
