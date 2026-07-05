# two-level ML estimation with random slopes (the rv() modifier)
#
# YR July 2026 (with help from Claude)
#
# References:
# - Asparouhov, T., & Muthen, B. (2003). Full-information maximum-likelihood
#   estimation of general two-level latent variable models with missing data.
#   Unpublished manuscript.
# - Rockwood, N. J. (2020). Maximum likelihood estimation of multilevel
#   structural equation models with random slopes for latent covariates.
#   Psychometrika, 85(2), 275-300.
#
# This file implements the observed-covariates-only setting: all random
# slopes involve *observed* level-1 covariates, so all random effects are
# 'linear' (Rockwood, 2020, section 3.2) and the marginal (observed-data)
# loglikelihood is available in closed form (no numerical integration).
#
# NOTE: the likelihood is *conditional* on the covariates that carry random
# slopes (and on the between-level exogenous covariates): their joint
# distribution with the outcomes is no longer normal.
#
# The model for cluster j (complete data, single group):
#
#   v_b = (y_b, z) ~ N(mu_v(w_j), Sigma_v)     [between; z = random slopes]
#   y_ij | v_b ~ N(mu_w + P x_ij + Q_ij v_b, Sigma_w),   independent over i
#
# where
#   - y_b     = between components of the 'both-level' variables (pb)
#   - z       = the random slopes, latent variables at level 2 (q)
#   - w_j     = between-level exogenous covariates (conditioned upon)
#   - x_ij    = level-1 covariates with random slopes (conditioned upon)
#   - Sigma_w = within-implied covariance of the level-1 y-vector (the
#               random-slope cells in GLIST are zero, so the standard
#               LISREL formula gives the *conditional* covariance)
#   - P       = reduced-form *fixed* slopes of y on x (from the within GLIST)
#   - Q_ij    = [ E_b , L_1 x_(1)ij , ..., L_q x_(q)ij ] with
#               L_r = (Lambda_w (I - B_w)^{-1})[, lhs(r)]
#
# Because Q_ij is linear in x_ij, all cluster-wise sums that appear in the
# loglikelihood (and in the EM E-step) reduce to per-cluster crossproducts
# of (y, x) -- computed once (lav_mvn_cl_rs_stats), so that each evaluation
# of the loglikelihood costs O(nclusters * (p + q)^3), without per-case
# loops. (cfr. the 'partial sufficient statistics' of Mplus.)


# collect and validate all static information needed for the random-slope
# kernels; returns an index map ('rs_info') that is computed only once
# (and cached in lavcache[[g]]$rs)
lav_mvn_cl_rs_info <- function(lavmodel = NULL, lavpartable = NULL,
                               lavdata = NULL) {
  # single group only (for now)
  if (lavdata@ngroups > 1L) {
    lav_msg_stop(gettext(
      "random slopes (rv() modifier) are not supported (yet) for
       multiple-group models."))
  }
  if (lavdata@nlevels != 2L) {
    lav_msg_stop(gettext(
      "random slopes (rv() modifier) require a two-level model."))
  }

  # random slopes for latent covariates? not yet
  # (this is where the Rockwood (2020) quadrature will plug in later)
  if (length(lavmodel@rv.lv) > 0L) {
    lav_msg_stop(gettext(
      "random slopes for latent covariates are not supported yet; only
       observed (within-level) covariates can have a random slope."))
  }

  # path list: from the parameter table (with validation), or from
  # lavmodel@rv.ov (no validation; used post-estimation)
  if (!is.null(lavpartable)) {
    pt <- lavpartable
    if (is.null(pt$level)) {
      lav_msg_stop(gettext("parameter table has no level column."))
    }

    # the random-slope paths at level 1
    rs_idx <- which(nchar(pt$rv) > 0L & pt$op == "~" & pt$level == 1L)
    if (length(rs_idx) == 0L) {
      lav_msg_stop(gettext("no random-slope (rv) regressions found at
                            level 1."))
    }
    # random slopes at level 2 are not allowed
    bad_idx <- which(nchar(pt$rv) > 0L & pt$op == "~" & pt$level > 1L)
    if (length(bad_idx) > 0L) {
      lav_msg_stop(gettext(
        "random slopes (rv() modifier) can only be used in the level-1
         (within) part of the model."))
    }

    path_rv <- pt$rv[rs_idx]
    path_lhs <- pt$lhs[rs_idx]
    path_rhs <- pt$rhs[rs_idx]
  } else {
    pt <- NULL
    path_rv <- names(lavmodel@rv.ov)
    path_lhs <- sapply(lavmodel@rv.ov, "[[", 1L)
    path_rhs <- sapply(lavmodel@rv.ov, "[[", 2L)
  }

  z_names <- unique(path_rv) # canonical order of the random slopes
  q <- length(z_names)

  # variable names per level
  ov_names_1 <- lavdata@ov.names.l[[1]][[1]] # within
  ov_names_2 <- lavdata@ov.names.l[[1]][[2]] # between
  both_names <- ov_names_1[ov_names_1 %in% ov_names_2]
  within_only <- ov_names_1[!ov_names_1 %in% ov_names_2]
  between_only <- ov_names_2[!ov_names_2 %in% ov_names_1]

  # the covariates carrying random slopes
  x_names <- unique(path_rhs)

  # validation (only if we have a parameter table)
  if (!is.null(pt)) {
    # exogenous variables (fixed.x)
    ov_x_1 <- lav_pt_vnames(pt, "ov.x", block = 1L)
    ov_x_2 <- lav_pt_vnames(pt, "ov.x", block = 2L)

    # latent variables at level 1
    lv_names_1 <- lav_pt_vnames(pt, "lv", block = 1L)

    # the rv covariates must be observed, within-only, and exogenous
    bad <- x_names[!x_names %in% within_only]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "covariate(s) with a random slope must be observed within-only
         variables (they should not appear in the level-2 part of the
         model): %s.", paste(bad, collapse = " ")))
    }
    bad <- x_names[!x_names %in% ov_x_1]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "covariate(s) with a random slope must be exogenous (fixed.x)
         variables: %s.", paste(bad, collapse = " ")))
    }

    # the lhs of each random-slope path must be a level-1 dependent
    # variable
    bad <- path_lhs[!(path_lhs %in% lv_names_1 |
                      path_lhs %in% ov_names_1)]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "left-hand side of random-slope regression(s) not found among the
         level-1 (latent or observed) variables: %s.",
        paste(bad, collapse = " ")))
    }
    bad <- path_lhs[path_lhs %in% x_names]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "left-hand side of random-slope regression(s) cannot be an
         exogenous covariate: %s.", paste(bad, collapse = " ")))
    }

    # the random slopes must be 'pure' latent variables at level 2:
    # no (real) indicators
    bad_idx <- which(pt$op == "=~" & pt$lhs %in% z_names &
                     pt$rhs != pt$lhs)
    if (length(bad_idx) > 0L) {
      lav_msg_stop(gettextf(
        "random slope(s) cannot have indicators at level 2: %s.",
        paste(unique(pt$lhs[bad_idx]), collapse = " ")))
    }

    # between-only *endogenous* observed variables (level-2 'z' outcomes
    # in Rockwood's notation) are not supported yet
    bad <- between_only[!between_only %in% ov_x_2]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "between-only endogenous (non-exogenous) observed variables are
         not supported (yet) in combination with random slopes: %s.",
        paste(bad, collapse = " ")))
    }
  }

  # the within y-vector: all within variables, except the rv covariates
  y_names <- ov_names_1[!ov_names_1 %in% x_names]
  p1 <- length(y_names)

  # path table
  path_tab <- data.frame(
    rv = path_rv,
    lhs = path_lhs,
    rhs = path_rhs,
    z.idx = match(path_rv, z_names),
    x.idx = match(path_rhs, x_names),
    stringsAsFactors = FALSE
  )

  # data column indices (in lavdata@X[[1]], columns = lavdata@ov.names[[1]])
  ov_names_data <- lavdata@ov.names[[1]]
  y_data_idx <- match(y_names, ov_names_data)
  x_data_idx <- match(x_names, ov_names_data)
  exo_b_names <- between_only # all between-only vars are exogenous here
  exo_b_data_idx <- match(exo_b_names, ov_names_data)

  # the both-level variables
  yb_names <- both_names
  pb <- length(yb_names)
  # position of the y_b variables within the y-vector
  yb_y_idx <- match(yb_names, y_names)

  # the between random vector v_b = (eta_b, eps_b) lives in 'factor'
  # space: the between latent variables (excluding the dummy lv's of the
  # between-level exogenous covariates, which are conditioned out),
  # followed by the between residuals (eps_b) of the both-level
  # variables -- but only those with a free (or fixed to a nonzero
  # value) residual (co)variance. This keeps v_b minimal, and Sigma_v
  # positive definite even when (some) between residual variances are
  # fixed to zero (as in Mplus ex9.8).
  mm_b <- 1:lavmodel@nmat[2] + lavmodel@nmat[1]
  glist_names_b <- names(lavmodel@GLIST)[mm_b]
  theta_mm <- mm_b[which(glist_names_b == "theta")]
  theta_b <- lavmodel@GLIST[[theta_mm]]
  free_mat <- matrix(FALSE, nrow(theta_b), ncol(theta_b))
  free_mat[lavmodel@m.free.idx[[theta_mm]]] <- TRUE
  nonzero_mat <- (theta_b != 0) | free_mat
  ov_names_b_theta <- lavmodel@dimNames[[theta_mm]][[1]]
  yb_b_idx <- match(yb_names, ov_names_b_theta)
  eps_flag <- apply(nonzero_mat[yb_b_idx, , drop = FALSE], 1L, any)
  eps_names <- yb_names[eps_flag]

  list(
    z.names = z_names, q = q,
    x.names = x_names, nx = length(x_names),
    y.names = y_names, p1 = p1,
    yb.names = yb_names, pb = pb, yb.y.idx = yb_y_idx,
    eps.names = eps_names, neps = length(eps_names),
    exo.b.names = exo_b_names, nexo.b = length(exo_b_names),
    path.tab = path_tab,
    y.data.idx = y_data_idx,
    x.data.idx = x_data_idx,
    exo.b.data.idx = exo_b_data_idx
  )
}

# per-cluster crossproduct statistics; computed once
lav_mvn_cl_rs_stats <- function(y1 = NULL, lp = NULL, rs_info = NULL) {
  cluster_idx <- lp$cluster.idx[[2]]
  nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]

  p1 <- rs_info$p1
  nx <- rs_info$nx

  yy <- y1[, rs_info$y.data.idx, drop = FALSE]
  xx <- y1[, rs_info$x.data.idx, drop = FALSE]

  # per-cluster sums and crossproducts
  sy <- rowsum.default(yy, cluster_idx, reorder = TRUE) # G x p1
  sx <- rowsum.default(xx, cluster_idx, reorder = TRUE) # G x nx

  sxx <- vector("list", nclusters)
  sxy <- vector("list", nclusters)
  syy <- vector("list", nclusters)

  # split row indices per cluster (loop over clusters; matrices are small)
  row_idx <- split.default(seq_len(nrow(y1)), cluster_idx)
  for (j in seq_len(nclusters)) {
    idx <- row_idx[[j]]
    yj <- yy[idx, , drop = FALSE]
    xj <- xx[idx, , drop = FALSE]
    syy[[j]] <- crossprod(yj)
    sxx[[j]] <- crossprod(xj)
    sxy[[j]] <- crossprod(xj, yj)
  }

  # between-level exogenous covariates: one row per cluster
  if (rs_info$nexo.b > 0L) {
    first_idx <- which(!duplicated(cluster_idx))
    # order by cluster index (1..G)
    first_idx <- first_idx[order(cluster_idx[first_idx])]
    exo_b <- y1[first_idx, rs_info$exo.b.data.idx, drop = FALSE]
  } else {
    exo_b <- matrix(0, nclusters, 0L)
  }

  list(
    nclusters = nclusters, cluster.size = cluster_size,
    sy = sy, sx = sx, sxx = sxx, sxy = sxy, syy = syy,
    exo.b = exo_b, ntotal = sum(cluster_size)
  )
}

# model-implied ingredients for the random-slope loglikelihood,
# computed from GLIST
#
# returns: sigma.w (p1 x p1), mu.w (p1), P (p1 x nx),
#          lmat (p1 x npaths), sigma.v (pv x pv), mu.v (pv),
#          cc (pv x nexo.b), mu.exo (nexo.b)
lav_mvn_cl_rs_implied <- function(lavmodel = NULL, glist = NULL,
                                  rs_info = NULL) {
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }
  nmat <- lavmodel@nmat

  # within block
  mm_w <- 1:nmat[1]
  mlist_w <- glist[mm_w]
  dimnames_w <- lavmodel@dimNames[mm_w]
  names(dimnames_w) <- names(glist)[mm_w]

  ov_names_w <- dimnames_w[["lambda"]][[1]]
  lv_names_w <- dimnames_w[["lambda"]][[2]]

  y_idx <- match(rs_info$y.names, ov_names_w)
  x_lv_idx <- match(rs_info$x.names, lv_names_w)

  lambda_w <- mlist_w$lambda
  psi_w <- mlist_w$psi
  theta_w <- mlist_w$theta
  beta_w <- mlist_w$beta
  nu_w <- mlist_w$nu
  alpha_w <- mlist_w$alpha

  nlv_w <- ncol(lambda_w)
  if (is.null(beta_w)) {
    ib_inv_w <- diag(nlv_w)
  } else {
    ib_inv_w <- try(solve(diag(nlv_w) - beta_w), silent = TRUE)
    if (inherits(ib_inv_w, "try-error")) {
      return(NULL)
    }
  }

  # reduced form: ov on lv
  rf_w <- lambda_w %*% ib_inv_w

  # condition on the x-covariates: zero out their psi rows/cols and alpha
  psi_w0 <- psi_w
  alpha_w0 <- if (is.null(alpha_w)) {
    matrix(0, nlv_w, 1L)
  } else {
    alpha_w
  }
  if (length(x_lv_idx) > 0L) {
    psi_w0[x_lv_idx, ] <- 0
    psi_w0[, x_lv_idx] <- 0
    alpha_w0[x_lv_idx, 1L] <- 0
  }

  # Sigma_w: conditional (on x and z) covariance of the y-vector
  rf_y <- rf_w[y_idx, , drop = FALSE]
  sigma_w <- rf_y %*% psi_w0 %*% t(rf_y) +
    theta_w[y_idx, y_idx, drop = FALSE]

  # mu_w: within intercept part of the y-vector
  nu_w0 <- if (is.null(nu_w)) {
    matrix(0, nrow(lambda_w), 1L)
  } else {
    nu_w
  }
  mu_w <- as.numeric(nu_w0[y_idx, 1L] + rf_y %*% alpha_w0)

  # P: reduced-form *fixed* slopes of y on the x-covariates
  # (the random-slope cells in beta are zero, so only the fixed paths
  #  contribute)
  pmat <- rf_y[, x_lv_idx, drop = FALSE]

  # L_r columns, one per random-slope path
  path_tab <- rs_info$path.tab
  npaths <- nrow(path_tab)
  lmat <- matrix(0, rs_info$p1, npaths)
  for (i in seq_len(npaths)) {
    lhs_lv_idx <- match(path_tab$lhs[i], lv_names_w)
    if (is.na(lhs_lv_idx)) {
      # lhs is an observed variable without a dummy lv: direct row
      lhs_ov_idx <- match(path_tab$lhs[i], rs_info$y.names)
      lmat[lhs_ov_idx, i] <- 1
    } else {
      lmat[, i] <- rf_y[, lhs_lv_idx]
    }
  }

  # between block
  mm_b <- 1:nmat[2] + nmat[1]
  mlist_b <- glist[mm_b]
  dimnames_b <- lavmodel@dimNames[mm_b]
  names(dimnames_b) <- names(glist)[mm_b]

  ov_names_b <- dimnames_b[["lambda"]][[1]]
  lv_names_b <- dimnames_b[["lambda"]][[2]]

  yb_idx <- match(rs_info$yb.names, ov_names_b)
  z_lv_idx <- match(rs_info$z.names, lv_names_b)
  exo_ov_idx <- match(rs_info$exo.b.names, ov_names_b)

  lambda_b <- mlist_b$lambda
  psi_b <- mlist_b$psi
  theta_b <- mlist_b$theta
  beta_b <- mlist_b$beta
  nu_b <- mlist_b$nu
  alpha_b <- mlist_b$alpha

  nlv_b <- ncol(lambda_b)
  if (is.null(beta_b)) {
    ib_inv_b <- diag(nlv_b)
  } else {
    ib_inv_b <- try(solve(diag(nlv_b) - beta_b), silent = TRUE)
    if (inherits(ib_inv_b, "try-error")) {
      return(NULL)
    }
  }

  e_eta <- if (is.null(alpha_b)) {
    matrix(0, nlv_b, 1L)
  } else {
    ib_inv_b %*% alpha_b
  }
  v_eta <- ib_inv_b %*% psi_b %*% t(ib_inv_b)

  nu_b0 <- if (is.null(nu_b)) {
    matrix(0, nrow(lambda_b), 1L)
  } else {
    nu_b
  }

  # the between random vector v_b = (eta_b, eps_b) in 'factor' space:
  # - eta_b: all between latent variables (including the random slopes,
  #   and excluding the dummy lv's of the between-level exogenous
  #   covariates, which are conditioned out)
  # - eps_b: the between residuals of the both-level variables with a
  #   non-zero residual (co)variance (see rs_info$eps.names)
  pb <- rs_info$pb
  q <- rs_info$q
  nexo <- rs_info$nexo.b
  neps <- rs_info$neps

  exo_lv_idx <- match(rs_info$exo.b.names, lv_names_b)
  if (nexo > 0L) {
    eta_idx <- seq_len(nlv_b)[-exo_lv_idx]
  } else {
    eta_idx <- seq_len(nlv_b)
  }
  meta <- length(eta_idx)
  pv <- meta + neps

  # position of the random slopes within v_b
  z_v_idx <- match(z_lv_idx, eta_idx)

  # condition eta on the between-level exogenous covariates
  if (nexo > 0L) {
    v_ee <- v_eta[exo_lv_idx, exo_lv_idx, drop = FALSE]
    v_ee_inv <- lav_mat_sym_inverse(v_ee)
    cc_eta <- v_eta[eta_idx, exo_lv_idx, drop = FALSE] %*% v_ee_inv
    v_eta_c <- v_eta[eta_idx, eta_idx, drop = FALSE] -
      cc_eta %*% v_eta[exo_lv_idx, eta_idx, drop = FALSE]
    v_eta_c <- (v_eta_c + t(v_eta_c)) / 2
    mu_exo <- as.numeric(nu_b0[exo_ov_idx, 1L] + e_eta[exo_lv_idx, 1L])
  } else {
    cc_eta <- matrix(0, meta, 0L)
    v_eta_c <- v_eta
    mu_exo <- numeric(0L)
  }

  # Sigma_v and (unconditional part of) mu_v
  sigma_v <- matrix(0, pv, pv)
  sigma_v[seq_len(meta), seq_len(meta)] <- v_eta_c
  mu_v <- numeric(pv)
  mu_v[seq_len(meta)] <- e_eta[eta_idx, 1L]
  cc <- matrix(0, pv, nexo)
  cc[seq_len(meta), ] <- cc_eta
  if (neps > 0L) {
    eps_b_idx <- match(rs_info$eps.names, ov_names_b)
    eps_pos <- meta + seq_len(neps)
    sigma_v[eps_pos, eps_pos] <-
      theta_b[eps_b_idx, eps_b_idx, drop = FALSE]
  }

  # constant loading matrix of v_b into the y-vector
  q0 <- matrix(0, rs_info$p1, pv)
  yb_y_idx <- rs_info$yb.y.idx
  q0[yb_y_idx, seq_len(meta)] <- lambda_b[yb_idx, eta_idx, drop = FALSE]
  if (neps > 0L) {
    eps_y_idx <- match(rs_info$eps.names, rs_info$y.names)
    for (i in seq_len(neps)) {
      q0[eps_y_idx[i], meta + i] <- 1
    }
  }

  # constant part of the mean of the y-vector: within intercepts +
  # between (fixed) intercepts of the both-level variables
  mu_y <- mu_w
  mu_y[yb_y_idx] <- mu_y[yb_y_idx] + nu_b0[yb_idx, 1L]

  list(
    sigma.w = sigma_w, mu.y = mu_y, P = pmat, lmat = lmat,
    q0 = q0, z.v.idx = z_v_idx, pv = pv,
    sigma.v = sigma_v, mu.v = mu_v, cc = cc, mu.exo = mu_exo
  )
}

# -2 * (observed-data, conditional-on-x) loglikelihood
# computed from the per-cluster crossproduct statistics
lav_mvn_cl_rs_loglik <- function(rs_stats = NULL, imp = NULL,
                                 rs_info = NULL,
                                 sinv_method = "eigen",
                                 log2pi = TRUE,
                                 minus_two = TRUE,
                                 per_cluster = FALSE) {
  # non-PD or failed implied?
  if (is.null(imp)) {
    return(+Inf)
  }

  p1 <- rs_info$p1
  nx <- rs_info$nx
  pv <- imp$pv
  path_tab <- rs_info$path.tab
  npaths <- nrow(path_tab)

  nclusters <- rs_stats$nclusters
  nj_all <- rs_stats$cluster.size

  # invert sigma.w
  sigma_w_inv <- lav_mat_sym_inverse(imp$sigma.w,
    logdet = TRUE, sinv_method = sinv_method
  )
  logdet_w <- attr(sigma_w_inv, "logdet")
  if (!is.finite(logdet_w)) {
    return(+Inf)
  }

  # NOTE: we never invert sigma.v; we use the identities
  #   logdet(Sigma_v) + logdet(Sigma_v^{-1} + A) = logdet(I + A Sigma_v)
  #   (Sigma_v^{-1} + A)^{-1} = Sigma_v (I + A Sigma_v)^{-1}
  # so that the loglikelihood remains well-defined when Sigma_v is
  # singular (eg between residual variances fixed to zero) or even
  # (mildly) indefinite (Heywood territory), as long as every
  # cluster-level marginal covariance matrix remains positive definite
  sigma_v <- imp$sigma.v
  if (anyNA(sigma_v)) {
    return(+Inf)
  }
  # is sigma.v PSD? (fast path: no per-cluster eigenvalue checks needed)
  ev_v <- eigen(sigma_v, symmetric = TRUE, only.values = TRUE)$values
  v_psd <- ev_v[length(ev_v)] > -1e-12

  w_mat <- sigma_w_inv
  lmat <- imp$lmat
  pmat <- imp$P
  mu_y <- imp$mu.y
  q0 <- imp$q0
  z_v_idx <- imp$z.v.idx

  # precompute
  wl <- w_mat %*% lmat # p1 x npaths
  ltwl <- crossprod(lmat, wl) # npaths x npaths
  wq0 <- w_mat %*% q0 # p1 x pv
  q0twq0 <- crossprod(q0, wq0) # pv x pv
  q0twl <- crossprod(q0, wl) # pv x npaths
  wp <- w_mat %*% pmat # p1 x nx
  ptwp <- crossprod(pmat, wp) # nx x nx
  mu_y_w <- as.numeric(w_mat %*% mu_y) # p1

  path_z <- path_tab$z.idx # 1..q
  path_x <- path_tab$x.idx # 1..nx
  path_zcol <- z_v_idx[path_z] # column in v_b, per path

  # more constants (hoisted out of the cluster loop)
  pmat_t <- t(pmat)
  wp_t <- t(wp)
  mu_y_w_p <- as.numeric(mu_y_w %*% pmat) # nx
  mu_y_quad <- sum(mu_y_w * mu_y)
  dpv <- diag(pv)
  sigma_v_local <- sigma_v
  mu_v <- imp$mu.v
  cc_mat <- imp$cc
  mu_exo <- imp$mu.exo
  nexo <- rs_info$nexo.b
  exo_b <- rs_stats$exo.b

  sy_all <- rs_stats$sy
  sx_all <- rs_stats$sx
  sxx_all <- rs_stats$sxx
  sxy_all <- rs_stats$sxy
  syy_all <- rs_stats$syy

  loglik_j <- numeric(nclusters)
  ok <- TRUE

  ok <- tryCatch({
  for (j in seq_len(nclusters)) {
    nj <- nj_all[j]
    sy_j <- sy_all[j, ]
    sx_j <- sx_all[j, ]
    sxx_j <- sxx_all[[j]]
    sxy_j <- sxy_all[[j]]
    syy_j <- syy_all[[j]]

    # centered sums: e~_i = y_i - mu_y - P x_i
    # sum_i e~_i
    se_j <- sy_j - nj * mu_y - as.numeric(pmat %*% sx_j)
    # sum_i x_ai e~_i' (nx x p1)
    sxe_j <- sxy_j - outer(sx_j, mu_y) - sxx_j %*% pmat_t

    # sum_i e~_i' W e~_i
    quad0 <- sum(w_mat * syy_j) -
      2 * sum(mu_y_w * sy_j) +
      nj * mu_y_quad -
      2 * sum(wp_t * sxy_j) +
      2 * sum(mu_y_w_p * sx_j) +
      sum(ptwp * sxx_j)

    # Q_ij = Q0 + sum_p L_p x_pij e_{z(p)}'
    # A_j = sum_i Q_ij' W Q_ij; p_j = sum_i Q_ij' W e~_i

    # constant part
    a_j <- nj * q0twq0
    p_j <- as.numeric(crossprod(wq0, se_j))

    # cross part: Q0' W (sum_i Lx_i) (only the z columns are involved)
    cz <- matrix(0, pv, pv)
    for (i in seq_len(npaths)) {
      zcol <- path_zcol[i]
      cz[, zcol] <- cz[, zcol] + q0twl[, i] * sx_j[path_x[i]]
    }
    a_j <- a_j + cz + t(cz)

    # z x z part: sum_i (Lx_i)' W (Lx_i)
    for (i in seq_len(npaths)) {
      for (k in seq_len(npaths)) {
        a_j[path_zcol[i], path_zcol[k]] <-
          a_j[path_zcol[i], path_zcol[k]] +
          ltwl[i, k] * sxx_j[path_x[i], path_x[k]]
      }
    }

    # p_j z part: sum_p L_p' W (sum_i x_pi e~_i)
    for (i in seq_len(npaths)) {
      zi <- path_zcol[i]
      p_j[zi] <- p_j[zi] + sum(wl[, i] * sxe_j[path_x[i], ])
    }

    # mean of v_b for this cluster (given the between covariates)
    if (nexo > 0L) {
      d_j <- mu_v + as.numeric(cc_mat %*% (exo_b[j, ] - mu_exo))
    } else {
      d_j <- mu_v
    }

    # center v_b at its mean
    p_tilde <- p_j - as.numeric(a_j %*% d_j)
    quad_j <- quad0 - 2 * sum(d_j * p_j) + sum(d_j * (a_j %*% d_j))

    # Z_j = I + A_j Sigma_v; logdet(Z_j) and Sigma_v Z_j^{-1} p~
    z_j <- dpv + a_j %*% sigma_v_local
    if (v_psd) {
      dt <- determinant(z_j, logarithm = TRUE)
      if (dt$sign <= 0 || !is.finite(dt$modulus)) {
        ok <- FALSE
        break
      }
      logdet_z <- as.numeric(dt$modulus)
    } else {
      # sigma.v indefinite: check that all eigenvalues of Z_j are
      # (real and) positive, so that the cluster-level marginal
      # covariance matrix is positive definite
      lambda <- Re(eigen(z_j, only.values = TRUE)$values)
      if (any(lambda < sqrt(.Machine$double.eps))) {
        ok <- FALSE
        break
      }
      logdet_z <- sum(log(lambda))
    }
    zp_j <- solve(z_j, p_tilde)

    # -2 loglik for this cluster
    loglik_j[j] <- nj * logdet_w + logdet_z +
      quad_j - sum(p_tilde * (sigma_v_local %*% zp_j))
    if (log2pi) {
      loglik_j[j] <- loglik_j[j] + nj * p1 * log(2 * pi)
    }
  }
  ok
  }, error = function(e) FALSE) # tryCatch around the cluster loop

  if (!ok) {
    return(+Inf)
  }

  out <- sum(loglik_j)
  if (!minus_two) {
    out <- out / (-2)
    loglik_j <- loglik_j / (-2)
  }
  if (per_cluster) {
    attr(out, "loglik.cluster") <- loglik_j
  }

  out
}
