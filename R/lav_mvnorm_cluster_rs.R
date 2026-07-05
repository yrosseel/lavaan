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

  # structural fixed paths from the rv covariates? (a nonzero 'P'
  # matrix); the EM version cannot handle these (the per-case mean
  # P x_{ji} does not fit in the two-group complete-data model);
  # nlminb handles the general case
  mm_w <- 1:lavmodel@nmat[1]
  glist_names_w <- names(lavmodel@GLIST)[mm_w]
  beta_mm <- mm_w[which(glist_names_w == "beta")]
  p_flag <- FALSE
  if (length(beta_mm) > 0L) {
    beta_w <- lavmodel@GLIST[[beta_mm]]
    free_mat_w <- matrix(FALSE, nrow(beta_w), ncol(beta_w))
    free_mat_w[lavmodel@m.free.idx[[beta_mm]]] <- TRUE
    lv_names_w <- lavmodel@dimNames[[beta_mm]][[2]]
    x_lv_idx_w <- match(x_names, lv_names_w)
    x_lv_idx_w <- x_lv_idx_w[!is.na(x_lv_idx_w)]
    if (length(x_lv_idx_w) > 0L) {
      p_flag <- any(free_mat_w[, x_lv_idx_w]) ||
        any(beta_w[, x_lv_idx_w] != 0)
    }
  }

  # shared rv labels (the same slope attached to more than one path)?
  # the EM version requires a one-to-one path <-> slope mapping
  shared_flag <- (nrow(path_tab) != q) ||
    !identical(path_tab$z.idx, seq_len(q))

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
    exo.b.data.idx = exo_b_data_idx,
    p.flag = p_flag, shared.flag = shared_flag
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
  mu_y_b <- numeric(rs_info$p1)
  mu_y_b[yb_y_idx] <- nu_b0[yb_idx, 1L]
  mu_y <- mu_w + mu_y_b

  list(
    sigma.w = sigma_w, mu.y = mu_y, mu.y.b = mu_y_b, P = pmat,
    lmat = lmat, q0 = q0, z.v.idx = z_v_idx, pv = pv,
    sigma.v = sigma_v, mu.v = mu_v, cc = cc, mu.exo = mu_exo,
    eta.names = lv_names_b[eta_idx], meta = meta
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


# ---------------------------------------------------------------------
# the EM algorithm (Asparouhov & Muthen, 2003), complete data
# ---------------------------------------------------------------------
#
# complete data = observed data + the between random vector v_b (in
# factor space, see above) + the 'within slope' variables
# z_w[r] = v_b[z_r] * x_r (one per random-slope path; completed with an
# arbitrary standard-normal marginal, cfr. A&M section 2).
#
# - E-step: the posterior of v_b given the cluster data is
#       Sigma0_j = Sigma_v (I + A_j Sigma_v)^{-1},
#       mu0_j    = d_j + Sigma0_j (p_j - A_j d_j)
#   (the same A_j/p_j workhorses as the marginal loglikelihood); the
#   expected complete-data moments then collapse onto the same
#   per-cluster crossproducts (sy, sx, sxx, sxy, syy)
# - M-step: a standard TWO-GROUP SEM, fitted to the expected moments:
#   group 1 ('within', N obs) contains the y-vector plus the phantom
#   z_w variables (unit variance, zero mean, fixed loading 1 on the
#   lhs of their path); group 2 ('between', G obs) contains the
#   both-level variables, the random slopes (now 'observed') and the
#   between covariates
#
# limitations (fall back to nlminb): fixed paths from the rv
# covariates (P != 0), shared rv labels; between-level variances that
# are fixed to (near) zero are floored at 1e-4 in the M-step model
# (cfr. the Mplus 'minimum variance'); the nlminb hand-off (enabled by
# default) polishes the solution of the *exact* model afterwards

# check if the EM version can handle this model; returns character(0)
# if ok, otherwise the reason(s)
lav_mvn_cl_rs_em_ok <- function(rs_info = NULL) {
  reason <- character(0L)
  if (rs_info$p.flag) {
    reason <- c(reason, gettext(
      "the model contains fixed regression paths from the random-slope
       covariates"))
  }
  if (rs_info$shared.flag) {
    reason <- c(reason, gettext(
      "the same rv() label is attached to more than one regression
       path"))
  }
  reason
}

# E-step: expected complete-data moments
#
# returns:
#   r1 (p1 + q), t1 (p1 + q square): mean and covariance (about r1)
#     of the within complete-data vector (y_w, z_w) -- N observations
#   r2 (pb + q + nexo), t2 (square): mean and covariance of the
#     between complete-data vector (y_b, z, w) -- G observations
lav_mvn_cl_rs_estep <- function(rs_stats = NULL, imp = NULL,
                                rs_info = NULL,
                                sinv_method = "eigen") {
  p1 <- rs_info$p1
  nx <- rs_info$nx
  q <- rs_info$q
  pb <- rs_info$pb
  nexo <- rs_info$nexo.b
  pv <- imp$pv
  path_tab <- rs_info$path.tab
  npaths <- nrow(path_tab)

  nclusters <- rs_stats$nclusters
  nj_all <- rs_stats$cluster.size
  ntotal <- rs_stats$ntotal

  sigma_w_inv <- lav_mat_sym_inverse(imp$sigma.w,
    logdet = TRUE, sinv_method = sinv_method
  )
  if (!is.finite(attr(sigma_w_inv, "logdet"))) {
    lav_msg_stop(gettext(
      "within covariance matrix is not positive definite in the
       E-step."))
  }
  w_mat <- sigma_w_inv
  sigma_v <- imp$sigma.v
  mu_y <- imp$mu.y
  nu_map <- imp$mu.y.b
  q0 <- imp$q0
  z_v_idx <- imp$z.v.idx
  path_zcol <- z_v_idx[path_tab$z.idx]
  path_x <- path_tab$x.idx

  # constants
  wq0 <- w_mat %*% q0
  q0twq0 <- crossprod(q0, wq0)
  wl <- w_mat %*% imp$lmat
  q0twl <- crossprod(q0, wl)
  ltwl <- crossprod(imp$lmat, wl)
  dpv <- diag(pv)

  # between mapping: (y_b, z) = c0 + G v
  gmat <- rbind(
    q0[rs_info$yb.y.idx, , drop = FALSE],
    dpv[z_v_idx, , drop = FALSE]
  )
  c0 <- c(nu_map[rs_info$yb.y.idx], rep(0, q))

  # accumulators
  sum1_w <- numeric(p1 + q) # within first moments
  sum2_w <- matrix(0, p1 + q, p1 + q) # within second moments
  sum1_b <- numeric(pb + q + nexo)
  sum2_b <- matrix(0, pb + q + nexo, pb + q + nexo)

  y_idx1 <- seq_len(p1)
  z_idx1 <- p1 + seq_len(q)
  b_idx2 <- seq_len(pb + q)
  w_idx2 <- pb + q + seq_len(nexo)

  for (j in seq_len(nclusters)) {
    nj <- nj_all[j]
    sy_j <- rs_stats$sy[j, ]
    sx_j <- rs_stats$sx[j, ]
    sxx_j <- rs_stats$sxx[[j]]
    sxy_j <- rs_stats$sxy[[j]]
    syy_j <- rs_stats$syy[[j]]

    # A_j and p_j (centered at mu.y; no P term: the EM requires P = 0)
    se_j <- sy_j - nj * mu_y
    a_j <- nj * q0twq0
    p_j <- as.numeric(crossprod(wq0, se_j))
    cz <- matrix(0, pv, pv)
    for (i in seq_len(npaths)) {
      zcol <- path_zcol[i]
      cz[, zcol] <- cz[, zcol] + q0twl[, i] * sx_j[path_x[i]]
    }
    a_j <- a_j + cz + t(cz)
    for (i in seq_len(npaths)) {
      for (k in seq_len(npaths)) {
        a_j[path_zcol[i], path_zcol[k]] <-
          a_j[path_zcol[i], path_zcol[k]] +
          ltwl[i, k] * sxx_j[path_x[i], path_x[k]]
      }
    }
    sxe_j <- sxy_j - outer(sx_j, mu_y)
    for (i in seq_len(npaths)) {
      zi <- path_zcol[i]
      p_j[zi] <- p_j[zi] + sum(wl[, i] * sxe_j[path_x[i], ])
    }

    # prior mean of v_b for this cluster
    if (nexo > 0L) {
      w_j <- rs_stats$exo.b[j, ]
      d_j <- imp$mu.v + as.numeric(imp$cc %*% (w_j - imp$mu.exo))
    } else {
      w_j <- numeric(0L)
      d_j <- imp$mu.v
    }

    # posterior of v_b
    z_j <- dpv + a_j %*% sigma_v
    sigma0_j <- sigma_v %*% solve(z_j)
    sigma0_j <- (sigma0_j + t(sigma0_j)) / 2
    p_tilde <- p_j - as.numeric(a_j %*% d_j)
    mu0_j <- d_j + as.numeric(sigma0_j %*% p_tilde)
    m0_j <- sigma0_j + tcrossprod(mu0_j) # E[v v']

    # ---------- within complete-data moments: (y_w, z_w) ----------
    # y_w = a - Q0 v, with a = y - nu_map (observed)
    sa_j <- sy_j - nj * nu_map # sum_i a_i
    saa_j <- syy_j - outer(sy_j, nu_map) - outer(nu_map, sy_j) +
      nj * tcrossprod(nu_map) # sum_i a_i a_i'
    q0mu <- as.numeric(q0 %*% mu0_j)
    q0m0q0 <- q0 %*% m0_j %*% t(q0)

    sum1_w[y_idx1] <- sum1_w[y_idx1] + sa_j - nj * q0mu
    sum2_w[y_idx1, y_idx1] <- sum2_w[y_idx1, y_idx1] +
      saa_j - outer(sa_j, q0mu) - outer(q0mu, sa_j) + nj * q0m0q0

    # z_w[r] = v[z_r] * (the covariate of path r)
    # (one-to-one path <-> slope mapping, so path r defines slope r;
    #  but note that two slopes may share the same covariate)
    sxa_j <- sxy_j - outer(sx_j, nu_map) # sum_i x_ai a_i'
    for (r in seq_len(q)) {
      zc_r <- z_v_idx[r]
      xr <- path_x[r] # covariate (column in the x-block) of path r
      sum1_w[z_idx1[r]] <- sum1_w[z_idx1[r]] + mu0_j[zc_r] * sx_j[xr]
      for (s in seq_len(q)) {
        zc_s <- z_v_idx[s]
        xs <- path_x[s]
        sum2_w[z_idx1[r], z_idx1[s]] <- sum2_w[z_idx1[r], z_idx1[s]] +
          m0_j[zc_r, zc_s] * sxx_j[xr, xs]
      }
      # cross y_w x z_w[r]:
      # sum_i x_(xr)i ( a_i mu0[z_r] - Q0 M0[, z_r] )
      cross_r <- sxa_j[xr, ] * mu0_j[zc_r] -
        sx_j[xr] * as.numeric(q0 %*% m0_j[, zc_r])
      sum2_w[y_idx1, z_idx1[r]] <- sum2_w[y_idx1, z_idx1[r]] + cross_r
      sum2_w[z_idx1[r], y_idx1] <- sum2_w[z_idx1[r], y_idx1] + cross_r
    }

    # ---------- between complete-data moments: (y_b, z, w) ----------
    b_j <- c0 + as.numeric(gmat %*% mu0_j)
    bb_j <- gmat %*% m0_j %*% t(gmat) +
      outer(c0, b_j) + outer(b_j, c0) - tcrossprod(c0)

    sum1_b[b_idx2] <- sum1_b[b_idx2] + b_j
    sum2_b[b_idx2, b_idx2] <- sum2_b[b_idx2, b_idx2] + bb_j
    if (nexo > 0L) {
      sum1_b[w_idx2] <- sum1_b[w_idx2] + w_j
      sum2_b[b_idx2, w_idx2] <- sum2_b[b_idx2, w_idx2] +
        outer(b_j, w_j)
      sum2_b[w_idx2, b_idx2] <- sum2_b[w_idx2, b_idx2] +
        outer(w_j, b_j)
      sum2_b[w_idx2, w_idx2] <- sum2_b[w_idx2, w_idx2] +
        tcrossprod(w_j)
    }
  } # j

  r1 <- sum1_w / ntotal
  t1 <- sum2_w / ntotal - tcrossprod(r1)
  r2 <- sum1_b / nclusters
  t2 <- sum2_b / nclusters - tcrossprod(r2)

  list(r1 = r1, t1 = t1, r2 = r2, t2 = t2)
}

# construct the two-group M-step parameter table (computed once)
#
# returns: list(pt = the two-group partable (data.frame),
#               ov1/ov2 = the observed-variable order per group,
#               perm1/perm2 = permutations from the canonical E-step
#                             order to ov1/ov2,
#               free.idx = which(pt$free > 0),
#               fixed.x = any exo)
lav_mvn_cl_rs_mstep_pt <- function(lavpartable = NULL, rs_info = NULL,
                                   floor_value = 1e-4) {
  pt <- as.data.frame(lavpartable, stringsAsFactors = FALSE)
  # drop columns we do not need
  pt$est <- NULL
  pt$se <- NULL
  pt$start <- NULL
  pt$rv <- NULL
  pt$mod.idx <- NULL
  # two-group version: level -> group
  pt$group <- NULL
  names(pt)[names(pt) == "level"] <- "group"

  # remove all level-1 rows involving the rv covariates (their moments,
  # and the frozen random-slope paths)
  x_names <- rs_info$x.names
  rm_idx <- which(pt$group == 1L &
    (pt$lhs %in% x_names | pt$rhs %in% x_names))
  # remove the phantom marker rows (s =~ s) at the between level
  rm_idx <- c(rm_idx, which(pt$op == "=~" & pt$lhs == pt$rhs &
    pt$lhs %in% rs_info$z.names))
  if (length(rm_idx) > 0L) {
    pt <- pt[-rm_idx, , drop = FALSE]
  }

  # Mplus-style variance floor: between-level variances fixed to (near)
  # zero would make the complete-data model singular; floor them in the
  # M-step model only (the nlminb hand-off afterwards uses the exact
  # model)
  floor_idx <- which(pt$group == 2L & pt$op == "~~" &
    pt$lhs == pt$rhs & pt$free == 0L &
    !is.na(pt$ustart) & abs(pt$ustart) < floor_value)
  floored <- length(floor_idx) > 0L
  if (floored) {
    pt$ustart[floor_idx] <- floor_value
  }

  # phantom z_w variables at the within level; make sure the names do
  # not clash with existing variables
  znames <- paste0(rs_info$z.names, ".zw")
  taken <- unique(c(pt$lhs, pt$rhs))
  while (any(znames %in% taken)) {
    znames <- paste0(".", znames)
  }

  q <- rs_info$q
  path_tab <- rs_info$path.tab
  new_lhs <- character(0L)
  new_op <- character(0L)
  new_rhs <- character(0L)
  # the fixed unit loading of the path lhs on its z_w variable
  for (r in seq_len(q)) {
    new_lhs <- c(new_lhs, path_tab$lhs[r])
    new_op <- c(new_op, "~")
    new_rhs <- c(new_rhs, znames[r])
  }
  # unit variances, zero covariances, zero means (the standard-normal
  # completion of A&M)
  for (r in seq_len(q)) {
    for (s in r:q) {
      new_lhs <- c(new_lhs, znames[r])
      new_op <- c(new_op, "~~")
      new_rhs <- c(new_rhs, znames[s])
    }
  }
  new_lhs <- c(new_lhs, znames)
  new_op <- c(new_op, rep("~1", q))
  new_rhs <- c(new_rhs, rep("", q))
  nnew <- length(new_lhs)
  new_ustart <- numeric(nnew)
  new_ustart[seq_len(q)] <- 1 # the fixed loadings
  var_idx <- q + which(new_op[-seq_len(q)] == "~~" &
    new_lhs[-seq_len(q)] == new_rhs[-seq_len(q)])
  new_ustart[var_idx] <- 1 # the unit variances

  new_rows <- pt[rep(1L, nnew), , drop = FALSE]
  new_rows$id <- max(pt$id) + seq_len(nnew)
  new_rows$lhs <- new_lhs
  new_rows$op <- new_op
  new_rows$rhs <- new_rhs
  new_rows$user <- 2L
  new_rows$block <- 1L
  new_rows$group <- 1L
  new_rows$free <- 0L
  new_rows$ustart <- new_ustart
  new_rows$exo <- 0L
  if (!is.null(new_rows$label)) new_rows$label <- ""
  if (!is.null(new_rows$plabel)) new_rows$plabel <- ""
  if (!is.null(new_rows$efa)) new_rows$efa <- ""
  if (!is.null(new_rows$lower)) new_rows$lower <- -Inf
  if (!is.null(new_rows$upper)) new_rows$upper <- +Inf
  pt <- rbind(pt, new_rows)
  rownames(pt) <- NULL

  # observed-variable order per group
  ov1 <- lav_pt_vnames(pt, "ov", block = 1L)
  ov2 <- lav_pt_vnames(pt, "ov", block = 2L)

  # canonical E-step order: (y.names, znames) and (yb, z, exo)
  can1 <- c(rs_info$y.names, znames)
  can2 <- c(rs_info$yb.names, rs_info$z.names, rs_info$exo.b.names)
  if (length(ov1) != length(can1) || !all(ov1 %in% can1)) {
    lav_msg_fixme("within M-step variables do not match the E-step")
  }
  if (length(ov2) != length(can2) || !all(ov2 %in% can2)) {
    lav_msg_fixme("between M-step variables do not match the E-step")
  }

  list(
    pt = pt, ov1 = ov1, ov2 = ov2,
    perm1 = match(ov1, can1), perm2 = match(ov2, can2),
    free.idx = which(pt$free > 0L),
    fixed.x = any(pt$exo == 1L),
    znames = znames, floored = floored
  )
}

# EM estimation of the h0 model with random slopes
# (mirrors lav_mvn_cl_em_h0; single group, complete data)
lav_mvn_cl_rs_em_h0 <- function(lavsamplestats = NULL, lavdata = NULL,
                                lavpartable = NULL, lavmodel = NULL,
                                lavoptions = NULL, lavcache = NULL,
                                fx_tol = 1e-08, dx_tol = 1e-04,
                                max_iter = 5000L,
                                mstep_max_iter = 10000L,
                                mstep_rel_tol = 1e-10,
                                mstep_verbose = FALSE,
                                acceleration = "none") {
  stopifnot(lavdata@ngroups == 1L)

  # the rs cache
  rs <- lavcache[[1]]$rs
  if (is.null(rs)) {
    rs_info <- lav_mvn_cl_rs_info(lavmodel = lavmodel,
                                  lavdata = lavdata)
    rs_stats <- lav_mvn_cl_rs_stats(
      y1 = lavdata@X[[1]], lp = lavdata@Lp[[1]], rs_info = rs_info
    )
    rs <- list(info = rs_info, stats = rs_stats)
  }

  # feasibility
  reason <- lav_mvn_cl_rs_em_ok(rs$info)
  if (length(reason) > 0L) {
    lav_msg_stop(gettextf(
      "the EM algorithm is not available for this random-slope model
       (%s); use optim.method = \"nlminb\" instead.",
      paste(reason, collapse = "; ")))
  }

  # the two-group M-step model (built once)
  mstep <- lav_mvn_cl_rs_mstep_pt(
    lavpartable = lavpartable, rs_info = rs$info
  )
  nobs_list <- list(rs$stats$ntotal, rs$stats$nclusters)

  # M-step: fit the two-group model to the expected moments
  em_mstep <- function(estep, x_start) {
    local_pt <- mstep$pt
    local_pt$ustart[mstep$free.idx] <- x_start
    t1 <- estep$t1[mstep$perm1, mstep$perm1, drop = FALSE]
    r1 <- estep$r1[mstep$perm1]
    t2 <- estep$t2[mstep$perm2, mstep$perm2, drop = FALSE]
    r2 <- estep$r2[mstep$perm2]
    dimnames(t1) <- list(mstep$ov1, mstep$ov1)
    dimnames(t2) <- list(mstep$ov2, mstep$ov2)
    names(r1) <- mstep$ov1
    names(r2) <- mstep$ov2
    if (!mstep_verbose) {
      current_verbose <- lav_verbose()
      if (lav_verbose(FALSE)) {
        on.exit(lav_verbose(current_verbose), TRUE)
      }
    }
    lavaan(local_pt,
      sample_cov = list(within = t1, between = t2),
      sample_mean = list(within = r1, between = r2),
      sample_nobs = nobs_list,
      sample.cov.rescale = FALSE,
      control = list(
        iter.max = mstep_max_iter,
        rel.tol = mstep_rel_tol
      ),
      fixed.x = mstep$fixed.x,
      estimator = "ML",
      warn = FALSE,
      check.start = FALSE,
      check.post = FALSE,
      check.gradient = FALSE,
      check.vcov = FALSE,
      baseline = FALSE,
      h1 = FALSE,
      se = "none",
      test = "none"
    )
  }

  # one EM step as a map in the model parameter space: x -> x'
  em_step <- function(x, logl = FALSE) {
    lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
    imp <- lav_mvn_cl_rs_implied(lavmodel2, rs_info = rs$info)
    if (is.null(imp)) {
      lav_msg_stop(gettext(
        "model-implied matrices could not be computed in the E-step."))
    }
    estep <- lav_mvn_cl_rs_estep(
      rs_stats = rs$stats, imp = imp, rs_info = rs$info
    )
    local_fit <- em_mstep(estep, x_start = x)
    local_fit@optim$x
  }

  # observed-data loglikelihood at x
  em_logl <- function(x) {
    lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
    imp <- lav_mvn_cl_rs_implied(lavmodel2, rs_info = rs$info)
    out <- lav_mvn_cl_rs_loglik(
      rs_stats = rs$stats, imp = imp, rs_info = rs$info,
      log2pi = TRUE, minus_two = FALSE
    )
    as.numeric(out)
  }

  # initial values
  x_current <- lav_model_get_parameters(lavmodel)
  fx <- em_logl(x_current)

  if (lav_verbose()) {
    cat(
      "EM iter:", sprintf("%3d", 0),
      " fx =", sprintf("%17.10f", fx),
      "\n"
    )
  }

  if (acceleration %in% c("squarem", "qn")) {
    em_conv <- function(theta_old, theta, fx_old, fx) {
      if ((fx - fx_old) < fx_tol) {
        if (lav_verbose()) {
          cat("EM stopping rule reached: fx.delta < ", fx_tol, "\n")
        }
        return(TRUE)
      }
      # the gradient is numerical (2 x npar objective evaluations):
      # only check it when the loglikelihood is no longer moving much
      if ((fx - fx_old) > 0.01) {
        return(FALSE)
      }
      lavmodel2 <- lav_model_set_parameters(lavmodel, x = theta)
      dx <- lav_model_grad(lavmodel2,
        lavdata = lavdata,
        lavsamplestats = lavsamplestats,
        lavcache = lavcache
      )
      if (max(abs(dx)) < dx_tol) {
        if (lav_verbose()) {
          cat("EM stopping rule reached: max.dx < ", dx_tol, "\n")
        }
        return(TRUE)
      }
      FALSE
    }
    accel_fn <- if (acceleration == "qn") lav_em_qn else lav_em_squarem
    out_em <- accel_fn(
      theta = x_current, step_fn = em_step, logl_fn = em_logl,
      fx0 = fx, conv_fn = em_conv, max_iter = max_iter,
      logl_dec = 0 # inexact M-step: only accept monotone extrapolations
    )
    x <- out_em$theta
    fx <- out_em$fx
    em_converged <- out_em$converged || out_em$stalled
    em_iterations <- out_em$fpeval
  } else {
    # plain EM iterations
    fx_old <- fx
    x <- x_current
    for (i in 1:max_iter) {
      x <- em_step(x)
      fx <- em_logl(x)
      fx_delta <- fx - fx_old

      # the gradient is numerical (2 x npar objective evaluations):
      # only check it when the loglikelihood is no longer moving much
      if (fx_delta < 0.01) {
        lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
        dx <- lav_model_grad(lavmodel2,
          lavdata = lavdata,
          lavsamplestats = lavsamplestats,
          lavcache = lavcache
        )
        max_dx <- max(abs(dx))
      } else {
        max_dx <- as.numeric(NA)
      }

      if (lav_verbose()) {
        cat(
          "EM iter:", sprintf("%3d", i),
          " fx =", sprintf("%17.10f", fx),
          " fx.delta =", sprintf("%9.8f", fx_delta),
          " max.dx = ", sprintf("%9.8f", max_dx),
          "\n"
        )
      }

      if (fx_delta < fx_tol) {
        if (lav_verbose()) {
          cat("EM stopping rule reached: fx.delta < ", fx_tol, "\n")
        }
        break
      } else {
        fx_old <- fx
      }
      if (is.finite(max_dx) && max_dx < dx_tol) {
        if (lav_verbose()) {
          cat("EM stopping rule reached: max.dx < ", dx_tol, "\n")
        }
        break
      }
    }
    em_converged <- i < max_iter
    em_iterations <- i
  }

  # add attributes
  if (em_converged) {
    attr(x, "converged") <- TRUE
    attr(x, "warn.txt") <- ""
  } else {
    attr(x, "converged") <- FALSE
    attr(x, "warn.txt") <- paste("maximum number of iterations (",
      max_iter, ") ",
      "was reached without convergence.\n",
      sep = ""
    )
  }
  attr(x, "iterations") <- em_iterations
  attr(x, "control") <- list(
    em.args = list(
      max_iter = max_iter,
      fx_tol = fx_tol,
      dx_tol = dx_tol,
      acceleration = acceleration
    )
  )
  attr(fx, "fx.group") <- fx # single group for now
  attr(x, "fx") <- fx
  # if the M-step model used the variance floor, the EM fixed point is
  # NOT the optimum of the exact model: the nlminb hand-off must always
  # run (even if the gradient looks small)
  attr(x, "rs.floored") <- mstep$floored

  x
}


# ---------------------------------------------------------------------
# per-cluster scores and the first-order (`meat') information
# ---------------------------------------------------------------------

# per-cluster scores: d loglik_j / d theta (natural loglikelihood
# scale), computed numerically (central differences) from the
# per-cluster loglikelihood vector of lav_mvn_cl_rs_loglik()
#
# returns a (nclusters x npar) matrix
lav_mvn_cl_rs_scores <- function(lavmodel = NULL, rs = NULL,
                                 h = 1e-05) {
  x0 <- lav_model_get_parameters(lavmodel)
  npar <- length(x0)
  nclusters <- rs$stats$nclusters

  ll_cluster <- function(x) {
    lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
    imp <- lav_mvn_cl_rs_implied(lavmodel2, rs_info = rs$info)
    out <- lav_mvn_cl_rs_loglik(
      rs_stats = rs$stats, imp = imp, rs_info = rs$info,
      log2pi = FALSE, minus_two = FALSE, per_cluster = TRUE
    )
    attr(out, "loglik.cluster")
  }

  sc <- matrix(0, nclusters, npar)
  for (k in seq_len(npar)) {
    h_k <- h * max(1, abs(x0[k]))
    x_right <- x_left <- x0
    x_right[k] <- x0[k] + h_k
    x_left[k] <- x0[k] - h_k
    ll_right <- ll_cluster(x_right)
    ll_left <- ll_cluster(x_left)
    if (is.null(ll_right) || is.null(ll_left)) {
      lav_msg_stop(gettext(
        "per-cluster scores could not be computed (non-finite
         loglikelihood near the parameter estimates)."))
    }
    sc[, k] <- (ll_right - ll_left) / (2 * h_k)
  }

  sc
}


# ---------------------------------------------------------------------
# empirical Bayes (posterior mean) predictions
# ---------------------------------------------------------------------
#
# - level 2: the posterior means mu0_j of the between latent variables
#   (the between factors AND the random slopes) -- the same posterior
#   as the EM E-step
# - level 1: the within factor scores, conditional on the posterior
#   means of the between random vector; because everything is linear
#   and Gaussian, plugging in mu0_j gives the *exact* posterior mean
#   E(eta_wij | all data of cluster j)
lav_mvn_cl_rs_eb <- function(lavmodel = NULL, lavdata = NULL,
                             lavcache = NULL, se = FALSE) {
  # the rs cache
  rs <- lavcache[[1]]$rs
  if (is.null(rs)) {
    rs_info <- lav_mvn_cl_rs_info(lavmodel = lavmodel,
                                  lavdata = lavdata)
    rs_stats <- lav_mvn_cl_rs_stats(
      y1 = lavdata@X[[1]], lp = lavdata@Lp[[1]], rs_info = rs_info
    )
    rs <- list(info = rs_info, stats = rs_stats)
  }
  rs_info <- rs$info
  rs_stats <- rs$stats

  imp <- lav_mvn_cl_rs_implied(lavmodel, rs_info = rs_info)
  if (is.null(imp)) {
    lav_msg_stop(gettext(
      "model-implied matrices could not be computed."))
  }

  p1 <- rs_info$p1
  nx <- rs_info$nx
  pv <- imp$pv
  meta <- imp$meta
  path_tab <- rs_info$path.tab
  npaths <- nrow(path_tab)
  z_v_idx <- imp$z.v.idx
  path_zcol <- z_v_idx[path_tab$z.idx]
  path_x <- path_tab$x.idx
  nclusters <- rs_stats$nclusters
  nj_all <- rs_stats$cluster.size

  sigma_w_inv <- lav_mat_sym_inverse(imp$sigma.w, logdet = TRUE)
  if (!is.finite(attr(sigma_w_inv, "logdet"))) {
    lav_msg_stop(gettext(
      "within covariance matrix is not positive definite."))
  }
  w_mat <- sigma_w_inv
  sigma_v <- imp$sigma.v
  mu_y <- imp$mu.y

  # constants (as in the loglikelihood/E-step)
  wl <- w_mat %*% imp$lmat
  wq0 <- w_mat %*% imp$q0
  q0twq0 <- crossprod(imp$q0, wq0)
  q0twl <- crossprod(imp$q0, wl)
  ltwl <- crossprod(imp$lmat, wl)
  dpv <- diag(pv)

  # ---- per-cluster posterior of v_b: mu0_j (and posterior sd) ----
  mu0 <- matrix(0, nclusters, pv)
  sd0 <- if (se) matrix(0, nclusters, pv) else NULL
  for (j in seq_len(nclusters)) {
    nj <- nj_all[j]
    sy_j <- rs_stats$sy[j, ]
    sx_j <- rs_stats$sx[j, ]
    sxx_j <- rs_stats$sxx[[j]]
    sxy_j <- rs_stats$sxy[[j]]

    se_j <- sy_j - nj * mu_y
    a_j <- nj * q0twq0
    p_j <- as.numeric(crossprod(wq0, se_j))
    cz <- matrix(0, pv, pv)
    for (i in seq_len(npaths)) {
      zcol <- path_zcol[i]
      cz[, zcol] <- cz[, zcol] + q0twl[, i] * sx_j[path_x[i]]
    }
    a_j <- a_j + cz + t(cz)
    for (i in seq_len(npaths)) {
      for (k in seq_len(npaths)) {
        a_j[path_zcol[i], path_zcol[k]] <-
          a_j[path_zcol[i], path_zcol[k]] +
          ltwl[i, k] * sxx_j[path_x[i], path_x[k]]
      }
    }
    sxe_j <- sxy_j - outer(sx_j, mu_y)
    for (i in seq_len(npaths)) {
      zi <- path_zcol[i]
      p_j[zi] <- p_j[zi] + sum(wl[, i] * sxe_j[path_x[i], ])
    }

    if (rs_info$nexo.b > 0L) {
      d_j <- imp$mu.v +
        as.numeric(imp$cc %*% (rs_stats$exo.b[j, ] - imp$mu.exo))
    } else {
      d_j <- imp$mu.v
    }
    z_j <- dpv + a_j %*% sigma_v
    sigma0_j <- sigma_v %*% solve(z_j)
    sigma0_j <- (sigma0_j + t(sigma0_j)) / 2
    p_tilde <- p_j - as.numeric(a_j %*% d_j)
    mu0[j, ] <- d_j + as.numeric(sigma0_j %*% p_tilde)
    if (se) {
      sd0[j, ] <- sqrt(pmax(diag(sigma0_j), 0))
    }
  }

  # level-2 scores: the between latent variables (eta part of v_b)
  l2 <- mu0[, seq_len(meta), drop = FALSE]
  colnames(l2) <- imp$eta.names
  if (se) {
    se2 <- sd0[, seq_len(meta), drop = FALSE]
    colnames(se2) <- imp$eta.names
  }

  # ---- level-1 scores: within factor scores given (mu0_j, x_ij) ----
  # within model pieces (same construction as lav_mvn_cl_rs_implied)
  nmat <- lavmodel@nmat
  mm_w <- 1:nmat[1]
  mlist_w <- lavmodel@GLIST[mm_w]
  dimnames_w <- lavmodel@dimNames[mm_w]
  names(dimnames_w) <- names(lavmodel@GLIST)[mm_w]
  ov_names_w <- dimnames_w[["lambda"]][[1]]
  lv_names_w <- dimnames_w[["lambda"]][[2]]
  y_idx <- match(rs_info$y.names, ov_names_w)
  x_lv_idx <- match(rs_info$x.names, lv_names_w)
  lambda_w <- mlist_w$lambda
  psi_w <- mlist_w$psi
  beta_w <- mlist_w$beta
  alpha_w <- mlist_w$alpha
  nlv_w <- ncol(lambda_w)
  if (is.null(beta_w)) {
    ib_inv_w <- diag(nlv_w)
  } else {
    ib_inv_w <- solve(diag(nlv_w) - beta_w)
  }
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

  # 'real' latent variables: not an observed variable (dummy lv), and
  # not an rv covariate
  ov_names_data <- lavdata@ov.names[[1]]
  real_idx <- which(!lv_names_w %in% ov_names_data)

  y1 <- lavdata@X[[1]]
  yy <- y1[, rs_info$y.data.idx, drop = FALSE]
  xx <- y1[, rs_info$x.data.idx, drop = FALSE]
  cluster_idx <- lavdata@Lp[[1]]$cluster.idx[[2]]
  n1 <- nrow(y1)

  if (length(real_idx) > 0L) {
    # residuals: y_ij - E(y_ij | mu0_j, x_ij)
    vhat <- mu0[cluster_idx, , drop = FALSE]
    yhat <- matrix(mu_y, n1, p1, byrow = TRUE) +
      xx %*% t(imp$P) + vhat %*% t(imp$q0)
    # per-case slope inputs: z-hat * x
    zx <- matrix(0, n1, npaths)
    for (i in seq_len(npaths)) {
      zx[, i] <- vhat[, path_zcol[i]] * xx[, path_x[i]]
      yhat <- yhat + outer(zx[, i], imp$lmat[, i])
    }
    e_mat <- yy - yhat

    # conditional mean of the (full) within lv vector:
    # (I - B)^{-1} [ alpha0 + x (in the x-dummy slots)
    #                       + zhat*x (in the lhs equations) ]
    t_inp <- matrix(alpha_w0[, 1L], n1, nlv_w, byrow = TRUE)
    if (length(x_lv_idx) > 0L) {
      t_inp[, x_lv_idx] <- t_inp[, x_lv_idx, drop = FALSE] + xx
    }
    for (i in seq_len(npaths)) {
      lhs_col <- match(path_tab$lhs[i], lv_names_w)
      if (!is.na(lhs_col)) {
        t_inp[, lhs_col] <- t_inp[, lhs_col] + zx[, i]
      }
    }
    m_full <- t_inp %*% t(ib_inv_w)

    # regression update: + Cov(eta, y) Sigma_w^{-1} e
    veta0 <- ib_inv_w %*% psi_w0 %*% t(ib_inv_w)
    ceta_y <- (veta0 %*% t(lambda_w[y_idx, , drop = FALSE]))[real_idx, ,
      drop = FALSE
    ]
    l1 <- m_full[, real_idx, drop = FALSE] +
      e_mat %*% (w_mat %*% t(ceta_y))
    colnames(l1) <- lv_names_w[real_idx]
  } else {
    l1 <- matrix(0, n1, 0L)
  }

  out <- list(l1 = l1, l2 = l2, cluster.idx = cluster_idx)
  if (se) {
    out$se2 <- se2
  }
  out
}
