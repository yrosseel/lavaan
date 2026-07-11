# two-level ML estimation with random slopes (the rv() modifier)
#
# YR July 2026 (with help from Claude)
#
# Reference:
# - Rockwood, N. J. (2020). Maximum likelihood estimation of multilevel
#   structural equation models with random slopes for latent covariates.
#   Psychometrika, 85(2), 275-300.
#
# Two settings are implemented:
#
# 1) all random slopes involve *observed* level-1 covariates: all random
#    effects are 'linear' (Rockwood, 2020, section 3.2) and the marginal
#    (observed-data) loglikelihood is available in closed form (no
#    numerical integration);
#
# 2) one or more random slopes involve a *latent* level-1 covariate:
#    these slopes are 'nonlinear' random effects; conditional on their
#    values, the model reduces to setting (1), so the loglikelihood is
#    obtained by (Gauss-Hermite) quadrature over the nonlinear slopes
#    only, while all linear random effects are still integrated out
#    analytically (Rockwood, 2020, sections 3.3-3.6); see the
#    'quadrature' section below (lav_mvn_cl_rs_loglik_nl).
#
# NOTE: the likelihood is *conditional* on the covariates that carry random
# slopes (and on the between-level exogenous covariates): their joint
# distribution with the outcomes is no longer normal.
#
# Missing data (missing = "ml") is handled throughout: the per-case
# weight matrix becomes pattern-specific (see lav_mvn_cl_rs_loglik_m),
# the E-step additionally completes the missing y-values (see
# lav_mvn_cl_rs_estep_m), and lavPredict() conditions on the observed
# coordinates only.
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
                               lavdata = NULL, lavoptions = NULL) {
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

  # number of quadrature nodes per dimension (latent covariates only)
  ngh <- 21L
  if (!is.null(lavoptions) && !is.null(lavoptions$integration.ngh)) {
    ngh <- as.integer(lavoptions$integration.ngh)
  }

  # path list: from the parameter table (with validation), or from
  # lavmodel@rv.ov/@rv.lv (no validation; used post-estimation)
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

    # covariate type per path: latent ("lv") or observed ("ov")
    # (this mirrors the rv.lv/rv.ov classification in lav_model();
    #  'split' both-level observed covariates count as latent: the
    #  random slope multiplies their latent within-cluster component)
    lv_names_1 <- lav_pt_vnames(pt, "lv", block = 1L)
    ov_names_pt_1 <- lav_pt_vnames(pt, "ov", block = 1L)
    ov_names_pt_2 <- lav_pt_vnames(pt, "ov", block = 2L)
    ov_12 <- ov_names_pt_2[ov_names_pt_2 %in% ov_names_pt_1]
    path_type <- ifelse(path_rhs %in% c(lv_names_1, ov_12),
                        "lv", "ov")
  } else {
    pt <- NULL
    rv_all <- c(lavmodel@rv.ov, lavmodel@rv.lv)
    path_rv <- names(rv_all)
    path_lhs <- sapply(rv_all, "[[", 1L)
    path_rhs <- sapply(rv_all, "[[", 2L)
    path_type <- rep(c("ov", "lv"),
                     c(length(lavmodel@rv.ov), length(lavmodel@rv.lv)))
  }

  z_names <- unique(path_rv) # canonical order of the random slopes
  q <- length(z_names)

  # slope type: a slope multiplying a latent covariate is a 'nonlinear'
  # random effect (numerical integration); a single label may not mix
  # covariate types
  z_type <- character(q)
  for (r in seq_len(q)) {
    tps <- unique(path_type[path_rv == z_names[r]])
    if (length(tps) > 1L) {
      lav_msg_stop(gettextf(
        "random-slope label %s is attached to both an observed and a
         latent covariate; this is not supported (yet).", z_names[r]))
    }
    z_type[r] <- tps
  }
  nl_z_idx <- which(z_type == "lv")
  nl_flag <- length(nl_z_idx) > 0L

  # variable names per level
  ov_names_1 <- lavdata@ov.names.l[[1]][[1]] # within
  ov_names_2 <- lavdata@ov.names.l[[1]][[2]] # between
  both_names <- ov_names_1[ov_names_1 %in% ov_names_2]
  within_only <- ov_names_1[!ov_names_1 %in% ov_names_2]
  between_only <- ov_names_2[!ov_names_2 %in% ov_names_1]

  # the observed covariates that are conditioned upon: the (observed)
  # rv covariates PLUS all other within-only exogenous (fixed.x)
  # covariates -- the loglikelihood is conditional on ALL of them
  # (cfr. Mplus TYPE = TWOLEVEL RANDOM, which conditions on all
  # observed x-covariates); the fixed regression paths of the
  # non-rv covariates enter through the reduced-form P matrix
  x_rv <- unique(path_rhs[path_type == "ov"])
  x_extra <- within_only[within_only %in% lavdata@ov.names.x[[1]] &
                         !within_only %in% x_rv]
  x_names <- c(x_rv, x_extra)

  # validation (only if we have a parameter table)
  if (!is.null(pt)) {
    # note on 'split' both-level covariates (rhs in ov_12): the random
    # slope multiplies the *latent* within-cluster component of the
    # covariate (latent decomposition/centering, cfr. Rockwood's
    # simulation model); the covariate is then jointly modeled: it is
    # part of the y-vector, its within dummy lv is the beta_w
    # injection target, and its between dummy lv is an ordinary linear
    # random effect inside v_b -- no special handling needed below

    # exogenous variables (fixed.x)
    ov_x_1 <- lav_pt_vnames(pt, "ov.x", block = 1L)
    ov_x_2 <- lav_pt_vnames(pt, "ov.x", block = 2L)

    # the observed rv covariates must be observed, within-only, and
    # exogenous
    bad <- x_rv[!x_rv %in% within_only]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "covariate(s) with a random slope must be observed within-only
         variables (they should not appear in the level-2 part of the
         model): %s.", paste(bad, collapse = " ")))
    }
    bad <- x_rv[!x_rv %in% ov_x_1]
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

  }

  # the within y-vector: all within variables, except the rv covariates
  y_names <- ov_names_1[!ov_names_1 %in% x_names]
  p1 <- length(y_names)

  # path table
  path_tab <- data.frame(
    rv = path_rv,
    lhs = path_lhs,
    rhs = path_rhs,
    type = path_type,
    z.idx = match(path_rv, z_names),
    x.idx = match(path_rhs, x_names), # NA for latent covariates
    stringsAsFactors = FALSE
  )

  # data column indices (in lavdata@X[[1]], columns = lavdata@ov.names[[1]])
  ov_names_data <- lavdata@ov.names[[1]]
  y_data_idx <- match(y_names, ov_names_data)
  x_data_idx <- match(x_names, ov_names_data)

  # between-only variables: exogenous covariates (w; conditioned
  # upon) vs *endogenous* z-outcomes (Rockwood's z-variables: observed
  # once per cluster, and modeled -- e.g., indicators of a
  # between-level factor)
  mm_b0 <- 1:lavmodel@nmat[2] + lavmodel@nmat[1]
  glist_names_b0 <- names(lavmodel@GLIST)[mm_b0]
  lambda_mm_b <- mm_b0[which(glist_names_b0 == "lambda")]
  lv_names_b_all <- lavmodel@dimNames[[lambda_mm_b]][[2]]
  if (!is.null(pt)) {
    exo_b_names <- between_only[between_only %in% ov_x_2]
  } else {
    # without a parameter table: the exogenous covariates appear as
    # dummy latent variables in the between block; the z-outcomes
    # (pure indicators) do not
    exo_b_names <- between_only[between_only %in% lv_names_b_all]
  }
  zb_names <- between_only[!between_only %in% exo_b_names]
  kz <- length(zb_names)
  exo_b_data_idx <- match(exo_b_names, ov_names_data)
  zb_data_idx <- match(zb_names, ov_names_data)

  # z-outcomes must be pure indicators: a z-variable that is itself
  # involved in a regression gets a dummy latent variable (and a
  # zero residual variance), which does not fit the 'extra
  # observation block' treatment below
  bad <- zb_names[zb_names %in% lv_names_b_all]
  if (length(bad) > 0L) {
    lav_msg_stop(gettextf(
      "between-only endogenous observed variables can only appear as
       indicators of a between-level latent variable (for now); not
       as outcome or predictor in a regression: %s.",
      paste(bad, collapse = " ")))
  }

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

  # z-outcomes: the residual block Sigma_z must be positive definite
  # (all diagonal elements free or fixed to a nonzero value), and may
  # not covary with the between residuals of the both-level variables
  if (kz > 0L) {
    zb_b_idx <- match(zb_names, ov_names_b_theta)
    bad <- zb_names[!nonzero_mat[cbind(zb_b_idx, zb_b_idx)]]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "between-only endogenous observed variables must have a free
         (or fixed to a nonzero value) residual variance: %s.",
        paste(bad, collapse = " ")))
    }
    if (pb > 0L &&
        any(nonzero_mat[zb_b_idx, yb_b_idx, drop = FALSE])) {
      lav_msg_stop(gettext(
        "residual covariances between between-only endogenous
         observed variables and the both-level variables are not
         supported (yet)."))
    }
  }

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

  # latent covariates: GLIST cell positions of the (fixed-to-zero)
  # within regression paths; the quadrature driver injects the node
  # values for the nonlinear slopes there
  nl_beta_row <- nl_beta_col <- integer(0L)
  nl_beta_mm <- integer(0L)
  if (nl_flag) {
    if (length(beta_mm) == 0L) {
      lav_msg_stop(gettext(
        "no beta matrix found in the within part of the model
         (random slopes for latent covariates)."))
    }
    lv_names_w_all <- lavmodel@dimNames[[beta_mm]][[2]]
    nl_path_idx <- which(path_type == "lv")
    nl_beta_row <- match(path_lhs[nl_path_idx], lv_names_w_all)
    nl_beta_col <- match(path_rhs[nl_path_idx], lv_names_w_all)
    if (anyNA(nl_beta_row) || anyNA(nl_beta_col)) {
      lav_msg_stop(gettext(
        "unable to locate the latent-covariate random-slope path(s) in
         the within beta matrix."))
    }
    nl_beta_mm <- beta_mm
  }

  list(
    z.names = z_names, q = q,
    x.names = x_names, nx = length(x_names),
    y.names = y_names, p1 = p1,
    yb.names = yb_names, pb = pb, yb.y.idx = yb_y_idx,
    eps.names = eps_names, neps = length(eps_names),
    exo.b.names = exo_b_names, nexo.b = length(exo_b_names),
    zb.names = zb_names, kz = kz, zb.data.idx = zb_data_idx,
    path.tab = path_tab,
    y.data.idx = y_data_idx,
    x.data.idx = x_data_idx,
    exo.b.data.idx = exo_b_data_idx,
    p.flag = p_flag, shared.flag = shared_flag,
    # latent covariates (nonlinear random slopes; quadrature)
    z.type = z_type, nl.flag = nl_flag,
    nl.z.idx = nl_z_idx, q.nl = length(nl_z_idx),
    nl.beta.mm = nl_beta_mm,
    nl.beta.row = nl_beta_row, nl.beta.col = nl_beta_col,
    ngh = ngh
  )
}

# per-cluster crossproduct statistics; computed once
#
# with missing data (missing = "ml"), the y-block may contain NAs; the
# crossproducts are then computed per (cluster, missing-pattern) cell
# and stored in the 'mp' component; the covariates (both the rv
# covariates and the between-level exogenous covariates) must be
# complete: cases with missing values on exogenous variables have
# already been removed by lav_data (fixed.x = TRUE)
lav_mvn_cl_rs_stats <- function(y1 = NULL, lp = NULL, rs_info = NULL) {
  cluster_idx <- lp$cluster.idx[[2]]
  nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]

  p1 <- rs_info$p1
  nx <- rs_info$nx

  yy <- y1[, rs_info$y.data.idx, drop = FALSE]
  xx <- y1[, rs_info$x.data.idx, drop = FALSE]

  if (anyNA(xx)) {
    lav_msg_stop(gettext(
      "unexpected missing values in the covariate(s) with a random
       slope; such cases should have been removed (fixed.x = TRUE)."))
  }

  # between-level exogenous covariates: one row per cluster
  if (rs_info$nexo.b > 0L) {
    first_idx <- which(!duplicated(cluster_idx))
    # order by cluster index (1..G)
    first_idx <- first_idx[order(cluster_idx[first_idx])]
    exo_b <- y1[first_idx, rs_info$exo.b.data.idx, drop = FALSE]
    if (anyNA(exo_b)) {
      lav_msg_stop(gettext(
        "unexpected missing values in the between-level exogenous
         covariate(s); such cases should have been removed
         (fixed.x = TRUE)."))
    }
  } else {
    exo_b <- matrix(0, nclusters, 0L)
  }

  # between-only endogenous z-outcomes: one row per cluster
  if (rs_info$kz > 0L) {
    first_idx <- which(!duplicated(cluster_idx))
    first_idx <- first_idx[order(cluster_idx[first_idx])]
    zb <- y1[first_idx, rs_info$zb.data.idx, drop = FALSE]
    if (anyNA(zb)) {
      lav_msg_stop(gettext(
        "missing values in between-only endogenous observed variables
         are not supported (yet) in combination with random slopes."))
    }
  } else {
    zb <- matrix(0, nclusters, 0L)
  }

  # missing data in the y-block?
  mp <- NULL
  if (anyNA(yy)) {
    obs <- !is.na(yy)
    nobs_case <- rowSums(obs)

    # cases with no observed outcome values carry no information (the
    # likelihood is conditional on the covariates)
    empty_idx <- which(nobs_case == 0L)
    if (length(empty_idx) > 0L) {
      lav_msg_warn(gettextf(
        "%s case(s) have missing values on all within-level outcome
         variables; they do not contribute to the loglikelihood.",
        length(empty_idx)))
      keep_idx <- seq_len(nrow(yy))[-empty_idx]
    } else {
      keep_idx <- seq_len(nrow(yy))
    }

    # NA -> 0: the padded crossproduct entries are annihilated by the
    # (embedded) per-pattern weight matrices later on
    yy0 <- yy
    yy0[is.na(yy0)] <- 0

    # missing patterns (over the y-block only)
    case_id <- apply(1L * obs, 1L, paste, collapse = "")
    pat_ids <- unique(case_id[keep_idx])
    npatterns <- length(pat_ids)

    pattern <- vector("list", npatterns)
    for (p in seq_len(npatterns)) {
      case_p <- keep_idx[case_id[keep_idx] == pat_ids[p]]
      o <- which(obs[case_p[1L], ])
      clus_p <- cluster_idx[case_p]

      yp <- yy0[case_p, , drop = FALSE]
      xp <- xx[case_p, , drop = FALSE]

      # per-cluster sums for this pattern (rows ordered by cluster id)
      sy_p <- rowsum.default(yp, clus_p, reorder = TRUE)
      sx_p <- rowsum.default(xp, clus_p, reorder = TRUE)
      n_p <- as.numeric(rowsum.default(
        rep(1, length(case_p)), clus_p, reorder = TRUE))
      j1 <- sort(unique.default(clus_p))

      idx_split <- split.default(seq_along(case_p), clus_p)
      nclus_p <- length(j1)
      sxx_p <- vector("list", nclus_p)
      sxy_p <- vector("list", nclus_p)
      syy_p <- vector("list", nclus_p)
      for (jj in seq_len(nclus_p)) {
        ii <- idx_split[[jj]]
        yjj <- yp[ii, , drop = FALSE]
        xjj <- xp[ii, , drop = FALSE]
        syy_p[[jj]] <- crossprod(yjj)
        sxx_p[[jj]] <- crossprod(xjj)
        sxy_p[[jj]] <- crossprod(xjj, yjj)
      }

      pattern[[p]] <- list(
        o = o, j1 = j1, n = n_p, sy = sy_p, sx = sx_p,
        sxx = sxx_p, sxy = sxy_p, syy = syy_p
      )
    }

    mp <- list(
      npatterns = npatterns, pattern = pattern,
      nobs = sum(nobs_case), empty.idx = empty_idx
    )

    return(list(
      nclusters = nclusters, cluster.size = cluster_size,
      sy = NULL, sx = NULL, sxx = NULL, sxy = NULL, syy = NULL,
      exo.b = exo_b, zb = zb, ntotal = sum(cluster_size),
      ntotal.used = sum(cluster_size) - length(empty_idx), mp = mp
    ))
  }

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

  list(
    nclusters = nclusters, cluster.size = cluster_size,
    sy = sy, sx = sx, sxx = sxx, sxy = sxy, syy = syy,
    exo.b = exo_b, zb = zb, ntotal = sum(cluster_size),
    ntotal.used = sum(cluster_size), mp = NULL
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

  # between-only endogenous z-outcomes: conditional on v_b,
  #   z_j ~ N(mu_z + G v_b, sigma_z)
  # with G the between loadings of the z-block on the eta part of v_b
  # (the z-outcomes are pure indicators: no loadings on the eps part,
  # no dummy latent variables)
  kz <- rs_info$kz
  gmat <- matrix(0, kz, pv)
  mu_z <- numeric(kz)
  sigma_z <- matrix(0, kz, kz)
  if (kz > 0L) {
    zb_idx <- match(rs_info$zb.names, ov_names_b)
    gmat[, seq_len(meta)] <- lambda_b[zb_idx, eta_idx, drop = FALSE]
    mu_z <- as.numeric(nu_b0[zb_idx, 1L])
    sigma_z <- theta_b[zb_idx, zb_idx, drop = FALSE]
  }

  list(
    sigma.w = sigma_w, mu.y = mu_y, mu.y.b = mu_y_b, P = pmat,
    lmat = lmat, q0 = q0, z.v.idx = z_v_idx, pv = pv,
    sigma.v = sigma_v, mu.v = mu_v, cc = cc, mu.exo = mu_exo,
    eta.names = lv_names_b[eta_idx], meta = meta,
    mu.z = mu_z, gmat = gmat, sigma.z = sigma_z
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
  # missing data? use the pattern-based version
  if (!is.null(rs_stats$mp)) {
    return(lav_mvn_cl_rs_loglik_m(
      rs_stats = rs_stats, imp = imp, rs_info = rs_info,
      sinv_method = sinv_method, log2pi = log2pi,
      minus_two = minus_two, per_cluster = per_cluster
    ))
  }

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

  # between-only endogenous z-outcomes: one extra (per-cluster)
  # Gaussian observation block, z_j | v_b ~ N(mu_z + G v_b, Sigma_z)
  kz <- rs_info$kz
  if (kz > 0L) {
    wzb <- lav_mat_sym_inverse(imp$sigma.z,
      logdet = TRUE, sinv_method = sinv_method
    )
    logdet_zb <- attr(wzb, "logdet")
    if (!is.finite(logdet_zb)) {
      return(+Inf)
    }
    wzg <- wzb %*% imp$gmat # kz x pv
    a_zb <- crossprod(imp$gmat, wzg) # pv x pv
    zb_all <- rs_stats$zb
    mu_zb <- imp$mu.z
  }

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

    # z-outcome contributions (the extra observation block)
    if (kz > 0L) {
      e_zb <- zb_all[j, ] - mu_zb
      a_j <- a_j + a_zb
      p_j <- p_j + as.numeric(crossprod(wzg, e_zb))
      quad0 <- quad0 + sum(e_zb * (wzb %*% e_zb))
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
    if (kz > 0L) {
      loglik_j[j] <- loglik_j[j] + logdet_zb
    }
    if (log2pi) {
      loglik_j[j] <- loglik_j[j] + (nj * p1 + kz) * log(2 * pi)
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


# -2 * (observed-data, conditional-on-x) loglikelihood with missing
# data in the y-block (missing = "ml")
#
# the only change wrt the complete-data version: the per-case weight
# matrix is now pattern-specific, W_p = Sigma_w[o,o]^{-1} embedded in
# a p1 x p1 zero matrix; because the zero rows/columns of W_p
# annihilate the zero-padded entries of the per-(cluster, pattern)
# crossproducts, all accumulation formulas carry over unchanged --
# they are simply applied per pattern, and summed per cluster before
# the marginalization over v_b
lav_mvn_cl_rs_loglik_m <- function(rs_stats = NULL, imp = NULL,
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
  pv <- imp$pv
  path_tab <- rs_info$path.tab
  npaths <- nrow(path_tab)

  nclusters <- rs_stats$nclusters
  mp <- rs_stats$mp

  sigma_v <- imp$sigma.v
  if (anyNA(sigma_v)) {
    return(+Inf)
  }
  # is sigma.v PSD? (fast path: no per-cluster eigenvalue checks needed)
  ev_v <- eigen(sigma_v, symmetric = TRUE, only.values = TRUE)$values
  v_psd <- ev_v[length(ev_v)] > -1e-12

  lmat <- imp$lmat
  pmat <- imp$P
  mu_y <- imp$mu.y
  q0 <- imp$q0
  z_v_idx <- imp$z.v.idx

  path_z <- path_tab$z.idx # 1..q
  path_x <- path_tab$x.idx # 1..nx
  path_zcol <- z_v_idx[path_z] # column in v_b, per path

  pmat_t <- t(pmat)
  dpv <- diag(pv)
  mu_v <- imp$mu.v
  cc_mat <- imp$cc
  mu_exo <- imp$mu.exo
  nexo <- rs_info$nexo.b
  exo_b <- rs_stats$exo.b
  log2pi_c <- if (log2pi) log(2 * pi) else 0

  # between-only endogenous z-outcomes (always complete)
  kz <- rs_info$kz
  if (kz > 0L) {
    wzb <- lav_mat_sym_inverse(imp$sigma.z,
      logdet = TRUE, sinv_method = sinv_method
    )
    logdet_zb <- attr(wzb, "logdet")
    if (!is.finite(logdet_zb)) {
      return(+Inf)
    }
    wzg <- wzb %*% imp$gmat # kz x pv
    a_zb <- crossprod(imp$gmat, wzg) # pv x pv
    zb_all <- rs_stats$zb
    mu_zb <- imp$mu.z
  }

  # per-cluster accumulators
  ld_j <- numeric(nclusters) # sum_p n_jp (logdet_wp + p_o log(2 pi))
  quad_j <- numeric(nclusters)
  a_arr <- array(0, dim = c(pv, pv, nclusters))
  p_all <- matrix(0, nclusters, pv)
  loglik_j <- numeric(nclusters)

  ok <- tryCatch({
    # ---- accumulation: loop over patterns, then clusters within ----
    for (p in seq_len(mp$npatterns)) {
      cell <- mp$pattern[[p]]
      o <- cell$o
      po <- length(o)

      w_o <- lav_mat_sym_inverse(imp$sigma.w[o, o, drop = FALSE],
        logdet = TRUE, sinv_method = sinv_method
      )
      logdet_o <- attr(w_o, "logdet")
      if (!is.finite(logdet_o)) {
        # caught by the enclosing tryCatch() -> objective returns +Inf
        lav_msg_stop(gettext("sigma.w is not positive definite."))
      }
      # embed in the full p1 x p1 space
      w_mat <- matrix(0, p1, p1)
      w_mat[o, o] <- w_o

      # per-pattern constants
      wl <- w_mat %*% lmat # p1 x npaths
      ltwl <- crossprod(lmat, wl) # npaths x npaths
      wq0 <- w_mat %*% q0 # p1 x pv
      q0twq0 <- crossprod(q0, wq0) # pv x pv
      q0twl <- crossprod(q0, wl) # pv x npaths
      wp <- w_mat %*% pmat # p1 x nx
      ptwp <- crossprod(pmat, wp) # nx x nx
      wp_t <- t(wp)
      mu_y_w <- as.numeric(w_mat %*% mu_y) # p1
      mu_y_w_p <- as.numeric(mu_y_w %*% pmat) # nx
      mu_y_quad <- sum(mu_y_w * mu_y)
      ld_p <- logdet_o + po * log2pi_c

      j1 <- cell$j1
      for (jj in seq_along(j1)) {
        j <- j1[jj]
        n_jp <- cell$n[jj]
        sy_jp <- cell$sy[jj, ]
        sx_jp <- cell$sx[jj, ]
        sxx_jp <- cell$sxx[[jj]]
        sxy_jp <- cell$sxy[[jj]]
        syy_jp <- cell$syy[[jj]]

        ld_j[j] <- ld_j[j] + n_jp * ld_p

        # centered sums: e~_i = y_i - mu_y - P x_i
        se_jp <- sy_jp - n_jp * mu_y - as.numeric(pmat %*% sx_jp)
        sxe_jp <- sxy_jp - outer(sx_jp, mu_y) - sxx_jp %*% pmat_t

        # sum_i e~_i' W_p e~_i
        quad_j[j] <- quad_j[j] + sum(w_mat * syy_jp) -
          2 * sum(mu_y_w * sy_jp) +
          n_jp * mu_y_quad -
          2 * sum(wp_t * sxy_jp) +
          2 * sum(mu_y_w_p * sx_jp) +
          sum(ptwp * sxx_jp)

        # A contribution
        a_jp <- n_jp * q0twq0
        p_jp <- as.numeric(crossprod(wq0, se_jp))
        cz <- matrix(0, pv, pv)
        for (i in seq_len(npaths)) {
          zcol <- path_zcol[i]
          cz[, zcol] <- cz[, zcol] + q0twl[, i] * sx_jp[path_x[i]]
        }
        a_jp <- a_jp + cz + t(cz)
        for (i in seq_len(npaths)) {
          for (k in seq_len(npaths)) {
            a_jp[path_zcol[i], path_zcol[k]] <-
              a_jp[path_zcol[i], path_zcol[k]] +
              ltwl[i, k] * sxx_jp[path_x[i], path_x[k]]
          }
        }
        for (i in seq_len(npaths)) {
          zi <- path_zcol[i]
          p_jp[zi] <- p_jp[zi] + sum(wl[, i] * sxe_jp[path_x[i], ])
        }

        a_arr[, , j] <- a_arr[, , j] + a_jp
        p_all[j, ] <- p_all[j, ] + p_jp
      }
    }

    # ---- marginalization over v_b, per cluster ----
    for (j in seq_len(nclusters)) {
      a_j <- matrix(a_arr[, , j], pv, pv)
      p_j <- p_all[j, ]

      # z-outcome contributions (the extra observation block)
      if (kz > 0L) {
        e_zb <- zb_all[j, ] - mu_zb
        a_j <- a_j + a_zb
        p_j <- p_j + as.numeric(crossprod(wzg, e_zb))
        quad_j[j] <- quad_j[j] + sum(e_zb * (wzb %*% e_zb))
        ld_j[j] <- ld_j[j] + logdet_zb + kz * log2pi_c
      }

      # mean of v_b for this cluster (given the between covariates)
      if (nexo > 0L) {
        d_j <- mu_v + as.numeric(cc_mat %*% (exo_b[j, ] - mu_exo))
      } else {
        d_j <- mu_v
      }

      # center v_b at its mean
      p_tilde <- p_j - as.numeric(a_j %*% d_j)
      quad_full <- quad_j[j] - 2 * sum(d_j * p_j) +
        sum(d_j * (a_j %*% d_j))

      # Z_j = I + A_j Sigma_v; logdet(Z_j) and Sigma_v Z_j^{-1} p~
      z_j <- dpv + a_j %*% sigma_v
      if (v_psd) {
        dt <- determinant(z_j, logarithm = TRUE)
        if (dt$sign <= 0 || !is.finite(dt$modulus)) {
          # caught by the enclosing tryCatch() -> objective returns +Inf
          lav_msg_stop(gettext("Z matrix is not positive definite."))
        }
        logdet_z <- as.numeric(dt$modulus)
      } else {
        # sigma.v indefinite: check that all eigenvalues of Z_j are
        # (real and) positive
        lambda <- Re(eigen(z_j, only.values = TRUE)$values)
        if (any(lambda < sqrt(.Machine$double.eps))) {
          # caught by the enclosing tryCatch() -> objective returns +Inf
          lav_msg_stop(gettext("Z matrix is not positive definite."))
        }
        logdet_z <- sum(log(lambda))
      }
      zp_j <- solve(z_j, p_tilde)

      # -2 loglik for this cluster (log 2 pi already in ld_j)
      loglik_j[j] <- ld_j[j] + logdet_z +
        quad_full - sum(p_tilde * (sigma_v %*% zp_j))
    }
    TRUE
  }, error = function(e) FALSE) # tryCatch around both loops

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
# latent covariates: quadrature over the nonlinear random slopes
# ---------------------------------------------------------------------
#
# random slopes that multiply a *latent* level-1 covariate enter the
# within model nonlinearly (a product of two random variables); the
# marginal loglikelihood is no longer available in closed form.
# Following Rockwood (2020), we split the between random vector v_b
# into the nonlinear slopes z_nl and the remaining ('linear') part
# v_lin. Conditional on z_nl = b:
#
#   - the within model is the *linear* random-slope model with the
#     slope values b injected into the (fixed-to-zero) beta_w cells:
#     the SAME lav_mvn_cl_rs_implied() gives sigma_w(b), lmat(b), P(b)
#   - v_lin | z_nl = b is normal with moments given by standard
#     Gaussian conditioning (a common covariance matrix, and a mean
#     that is linear in b)
#   - the conditional (on b) loglikelihood is computed by the SAME
#     crossproduct-based lav_mvn_cl_rs_loglik() -- including the
#     missing-data version
#
# the marginal loglikelihood is then obtained by Gauss-Hermite
# quadrature over z_nl only (q.nl dimensions, typically 1):
#
#   loglik_j = logsumexp_q [ loglik_j(b_q) + log w_q ]
#
# with nodes b_q = mu_nl + chol(Sigma_nl)' t_q placed at the prior of
# the nonlinear slopes (Rockwood, 2020, eq. 36-37). Per-node cost is
# one closed-form evaluation (the crossproduct statistics are shared
# across the nodes).
#
# when between-level covariates are related to the nonlinear slopes
# ('Case B', e.g. a regression s ~ w; Rockwood's PISA model), the
# prior means mu_nl_j become cluster-specific. We then keep ONE
# shared grid -- centered at the average prior mean, with the grid
# covariance inflated by the between-cluster spread of the prior
# means (so the grid covers every cluster's prior) -- and repair the
# weights per cluster with the exact density ratio:
#
#   loglik_j = logsumexp_q [ loglik_j(b_q) + log w_q + log r_jq ],
#   log r_jq = log phi(b_q; mu_nl_j, Sigma_nl)
#              + log|A| - log phi(t_q; 0, I),
#
# with b_q = c + A' t_q the shared grid. This keeps the per-node cost
# structure of Case A (the within matrices are rebuilt per node, not
# per cluster x node); without between-level covariates (Case A) the
# corrections are skipped entirely (bit-identical to the plain
# scheme), and when covariates are present but unrelated to the
# slopes they vanish up to roundoff. Note that Sigma_nl (the
# conditional-on-w covariance) is always common across clusters:
# conditioning on the between covariates shifts means only.

# single entry point for the -2 loglikelihood: dispatches between the
# closed-form kernel (observed covariates only) and the quadrature
# version (one or more latent covariates)
lav_mvn_cl_rs_m2ll <- function(lavmodel = NULL, glist = NULL, rs = NULL,
                               sinv_method = "eigen", log2pi = TRUE,
                               minus_two = TRUE, per_cluster = FALSE) {
  if (isTRUE(rs$info$nl.flag)) {
    return(lav_mvn_cl_rs_loglik_nl(
      lavmodel = lavmodel, glist = glist,
      rs_stats = rs$stats, rs_info = rs$info,
      sinv_method = sinv_method, log2pi = log2pi,
      minus_two = minus_two, per_cluster = per_cluster
    ))
  }
  imp <- lav_mvn_cl_rs_implied(
    lavmodel = lavmodel, glist = glist, rs_info = rs$info
  )
  lav_mvn_cl_rs_loglik(
    rs_stats = rs$stats, imp = imp, rs_info = rs$info,
    sinv_method = sinv_method, log2pi = log2pi,
    minus_two = minus_two, per_cluster = per_cluster
  )
}

# build the conditional-model ingredients at all quadrature nodes:
# the reduced index map (rs_info_c), and one 'imp' list per node with
# the node values injected into beta_w and v_b conditioned on the
# nonlinear slopes; when the prior means of the nonlinear slopes are
# cluster-specific (Case B), exo_b (the per-cluster between-covariate
# rows) is needed, and the returned 'logr' matrix (nclusters x
# nnodes) holds the per-cluster log density-ratio corrections;
# returns NULL if the base implied matrices fail or Sigma_nl is not
# positive definite (individual nodes may still fail: those imp
# entries are NULL)
lav_mvn_cl_rs_cond <- function(lavmodel = NULL, glist = NULL,
                               rs_info = NULL, exo_b = NULL) {
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }

  # base implied ingredients; only the between pieces (sigma.v, mu.v,
  # cc, z.v.idx) are used here -- they do not depend on the injected
  # beta_w cells
  imp0 <- lav_mvn_cl_rs_implied(
    lavmodel = lavmodel, glist = glist, rs_info = rs_info
  )
  if (is.null(imp0)) {
    return(NULL)
  }

  path_tab <- rs_info$path.tab
  nl_z <- rs_info$nl.z.idx # indices in z.names
  nl_v <- imp0$z.v.idx[nl_z] # positions in v_b
  lin_v <- setdiff(seq_len(imp0$pv), nl_v)

  # prior (between) moments of the nonlinear slopes; the quadrature
  # requires a positive-definite Sigma_nl
  s_nl <- imp0$sigma.v[nl_v, nl_v, drop = FALSE]
  mu_nl <- imp0$mu.v[nl_v]
  if (anyNA(s_nl)) {
    return(NULL)
  }
  ch <- try(chol(s_nl), silent = TRUE) # s_nl = t(ch) %*% ch
  if (inherits(ch, "try-error")) {
    return(NULL)
  }
  s_nl_inv <- chol2inv(ch)

  # cluster-specific shifts of the nonlinear prior means (Case B):
  # mu_nl_j = mu_nl + delta_j, delta_j = cc[nl, ] (w_j - mu_exo).
  # NOTE: the branch is taken whenever between-level covariates are
  # present (a *structural* decision), not only when the current
  # cc[nl, ] values are nonzero: a free slope-on-covariate parameter
  # that happens to be zero (e.g., at the starting values) still has
  # a nonzero *derivative* of the correction terms, which the
  # gradient must see. When cc[nl, ] = 0, the corrections vanish (up
  # to roundoff) and the general path reduces to the plain scheme.
  q_nl <- length(nl_v)
  delta <- NULL
  case_b <- rs_info$nexo.b > 0L
  if (case_b) {
    if (is.null(exo_b)) {
      lav_msg_stop(gettext(
        "internal error: exo_b is needed when the model contains both
         between-level covariates and nonlinear random slopes."))
    }
    wc <- exo_b - rep(imp0$mu.exo, each = nrow(exo_b))
    delta <- wc %*% t(imp0$cc[nl_v, , drop = FALSE]) # G x q_nl
  }

  # conditional distribution of the linear part given the slopes:
  # v_lin | z_nl = b ~ N(mu_lin + R (b - mu_nl) + cc_c (w_j - mu_exo),
  #                      sigma_v_c)
  # with cc_c the *conditional* regression on the between covariates
  # (equal to cc[lin, ] when the slopes do not depend on them)
  r_mat <- imp0$sigma.v[lin_v, nl_v, drop = FALSE] %*% s_nl_inv
  sigma_v_c <- imp0$sigma.v[lin_v, lin_v, drop = FALSE] -
    r_mat %*% imp0$sigma.v[nl_v, lin_v, drop = FALSE]
  sigma_v_c <- (sigma_v_c + t(sigma_v_c)) / 2
  cc_c <- imp0$cc[lin_v, , drop = FALSE]
  if (case_b) {
    cc_c <- cc_c - r_mat %*% imp0$cc[nl_v, , drop = FALSE]
  }

  # conditional index map: only the observed-covariate paths remain
  ov_rows <- which(path_tab$type == "ov")
  z_keep <- which(rs_info$z.type == "ov")
  path_tab_c <- path_tab[ov_rows, , drop = FALSE]
  path_tab_c$z.idx <- match(path_tab$rv[ov_rows],
                            rs_info$z.names[z_keep])
  rs_info_c <- rs_info
  rs_info_c$path.tab <- path_tab_c
  rs_info_c$z.names <- rs_info$z.names[z_keep]
  rs_info_c$q <- length(z_keep)
  rs_info_c$nl.flag <- FALSE
  # positions of the remaining (observed-type) slopes in reduced v_b
  z_v_idx_c <- match(imp0$z.v.idx[z_keep], lin_v)

  # Gauss-Hermite grid (dnorm kernel; tensor product over q.nl dims)
  ngh <- rs_info$ngh
  if (is.null(ngh)) {
    ngh <- 21L
  }
  gh <- lav_integration_gauss_hermite(n = ngh, dnorm = TRUE)
  if (q_nl == 1L) {
    tmat <- matrix(gh$x, nrow = 1L)
    logw <- log(gh$w)
  } else {
    tmat <- t(as.matrix(expand.grid(rep(list(gh$x), q_nl))))
    logw <- rowSums(as.matrix(
      expand.grid(rep(list(log(gh$w)), q_nl))
    ))
  }
  nnodes <- ncol(tmat)

  # grid placement: b_q = c + A' t_q
  # - Case A: c = mu_nl, A = chol(Sigma_nl) (no correction needed)
  # - Case B: c = average prior mean; grid covariance = Sigma_nl +
  #   between-cluster covariance of the prior means (so the grid
  #   covers every cluster's prior); exact per-cluster density-ratio
  #   corrections below
  if (case_b) {
    grid_c <- mu_nl + colMeans(delta)
    dc <- delta - rep(colMeans(delta), each = nrow(delta))
    s_grid <- s_nl + crossprod(dc) / nrow(delta)
    a_up <- try(chol(s_grid), silent = TRUE)
    if (inherits(a_up, "try-error")) {
      return(NULL)
    }
  } else {
    grid_c <- mu_nl
    a_up <- ch
  }
  bmat <- grid_c + t(a_up) %*% tmat # q_nl x nnodes

  # per-cluster log density-ratio corrections (Case B only):
  # log r_jq = log phi(b_q; mu_nl + delta_j, Sigma_nl)
  #            + log|A| - log phi(t_q; 0, I)
  logr <- NULL
  if (case_b) {
    mu_mat <- rep(mu_nl, each = nrow(delta)) + delta # G x q_nl
    mb <- s_nl_inv %*% bmat # q_nl x nnodes
    bmb <- colSums(bmat * mb) # nnodes
    mm <- rowSums((mu_mat %*% s_nl_inv) * mu_mat) # G
    cross <- mu_mat %*% mb # G x nnodes
    logdet_s <- 2 * sum(log(diag(ch)))
    logdet_a <- sum(log(diag(a_up)))
    tt <- colSums(tmat * tmat) # nnodes
    logr <- -0.5 * (outer(mm, bmb, "+") - 2 * cross) -
      0.5 * logdet_s + logdet_a +
      rep(0.5 * tt, each = nrow(delta))
  }

  # node-value injection map: nonlinear path -> beta_w cell
  nl_path_rows <- which(path_tab$type == "lv")
  path_nl_z <- match(path_tab$z.idx[nl_path_rows], nl_z)
  beta_mm <- rs_info$nl.beta.mm
  beta_cells <- cbind(rs_info$nl.beta.row, rs_info$nl.beta.col)

  imp_list <- vector("list", nnodes)
  glist_i <- glist
  for (i in seq_len(nnodes)) {
    b_i <- bmat[, i]
    glist_i[[beta_mm]][beta_cells] <- b_i[path_nl_z]
    imp_i <- lav_mvn_cl_rs_implied(
      lavmodel = lavmodel, glist = glist_i, rs_info = rs_info
    )
    if (is.null(imp_i)) {
      next
    }
    imp_list[[i]] <- list(
      sigma.w = imp_i$sigma.w, mu.y = imp_i$mu.y,
      mu.y.b = imp_i$mu.y.b, P = imp_i$P,
      lmat = imp_i$lmat[, ov_rows, drop = FALSE],
      q0 = imp_i$q0[, lin_v, drop = FALSE],
      z.v.idx = z_v_idx_c, pv = length(lin_v),
      sigma.v = sigma_v_c,
      mu.v = imp0$mu.v[lin_v] + as.numeric(r_mat %*% (b_i - mu_nl)),
      cc = cc_c,
      mu.exo = imp0$mu.exo,
      # z-outcomes: the slopes have no indicators, so the nl columns
      # of G are structurally zero and simply drop out
      mu.z = imp_i$mu.z,
      gmat = imp_i$gmat[, lin_v, drop = FALSE],
      sigma.z = imp_i$sigma.z
    )
  }

  list(rs_info_c = rs_info_c, imp = imp_list, logw = logw,
       nnodes = nnodes, logr = logr,
       # extras for the EAP (lavPredict) machinery
       bmat = bmat, lin.v = lin_v, nl.v = nl_v,
       eta.names = imp0$eta.names, meta = imp0$meta, pv = imp0$pv,
       z.v.idx = imp0$z.v.idx)
}

# -2 * loglikelihood by Gauss-Hermite quadrature over the nonlinear
# (latent-covariate) random slopes; all failure modes return +Inf
# (mirroring lav_mvn_cl_rs_loglik)
lav_mvn_cl_rs_loglik_nl <- function(lavmodel = NULL, glist = NULL,
                                    rs_stats = NULL, rs_info = NULL,
                                    sinv_method = "eigen",
                                    log2pi = TRUE,
                                    minus_two = TRUE,
                                    per_cluster = FALSE) {
  cond <- lav_mvn_cl_rs_cond(
    lavmodel = lavmodel, glist = glist, rs_info = rs_info,
    exo_b = rs_stats$exo.b
  )
  if (is.null(cond)) {
    return(+Inf)
  }

  nclusters <- rs_stats$nclusters
  nnodes <- cond$nnodes
  hmat <- matrix(-Inf, nclusters, nnodes)

  for (i in seq_len(nnodes)) {
    if (is.null(cond$imp[[i]])) {
      next
    }
    m2_i <- lav_mvn_cl_rs_loglik(
      rs_stats = rs_stats, imp = cond$imp[[i]],
      rs_info = cond$rs_info_c,
      sinv_method = sinv_method, log2pi = log2pi,
      minus_two = TRUE, per_cluster = TRUE
    )
    if (!is.finite(m2_i)) {
      next
    }
    hmat[, i] <- -0.5 * attr(m2_i, "loglik.cluster") + cond$logw[i]
  }

  # cluster-specific prior means (Case B): density-ratio corrections
  if (!is.null(cond$logr)) {
    hmat <- hmat + cond$logr
  }

  # per-cluster log-sum-exp over the nodes
  hmax <- apply(hmat, 1L, max)
  if (any(!is.finite(hmax))) {
    return(+Inf)
  }
  loglik_j <- hmax + log(rowSums(exp(hmat - hmax)))

  out <- sum(loglik_j)
  if (minus_two) {
    out <- -2 * out
    loglik_j <- -2 * loglik_j
  }
  if (per_cluster) {
    attr(out, "loglik.cluster") <- loglik_j
  }

  out
}

# analytic (adjoint-based) gradient of the quadrature loglikelihood:
#
#   d loglik_j / d theta = sum_i pi_ji d loglik_j(b_i) / d theta
#
# with pi_ji the posterior node weights (the derivative of the
# log-sum-exp); the conditional (per-node) derivatives reuse the
# closed-form adjoint pass (lav_mvn_cl_rs_dfphi) on the conditional
# model, chained through a small numerical Jacobian of the per-node
# ingredient map theta -> phi_c(theta; node i) -- which includes the
# dependence of the node values b_i(theta) = mu_nl + chol(Sigma_nl)'t_i
# on theta (no cluster loops); returns d(-2 loglik)/d theta, or the
# per-cluster score matrix (natural scale) when per_cluster = TRUE;
# NULL on failure (callers fall back to numerical differentiation)
lav_mvn_cl_rs_grad_nl <- function(lavmodel = NULL, rs = NULL,
                                  sinv_method = "eigen",
                                  per_cluster = FALSE,
                                  h = 1e-06) {
  rs_info <- rs$info
  rs_stats <- rs$stats

  cond0 <- lav_mvn_cl_rs_cond(lavmodel = lavmodel, rs_info = rs_info,
                              exo_b = rs_stats$exo.b)
  if (is.null(cond0)) {
    return(NULL)
  }
  nnodes <- cond0$nnodes
  nclusters <- rs_stats$nclusters
  case_b <- !is.null(cond0$logr)

  # per-node conditional logliks + adjoints at the current parameters
  hmat <- matrix(-Inf, nclusters, nnodes)
  gphi_list <- vector("list", nnodes)
  for (i in seq_len(nnodes)) {
    imp_i <- cond0$imp[[i]]
    if (is.null(imp_i)) {
      return(NULL)
    }
    m2_i <- lav_mvn_cl_rs_loglik(
      rs_stats = rs_stats, imp = imp_i, rs_info = cond0$rs_info_c,
      sinv_method = sinv_method, log2pi = FALSE,
      minus_two = TRUE, per_cluster = TRUE
    )
    if (!is.finite(m2_i)) {
      return(NULL)
    }
    hmat[, i] <- -0.5 * attr(m2_i, "loglik.cluster") + cond0$logw[i]
    gphi_list[[i]] <- lav_mvn_cl_rs_dfphi(
      rs_stats = rs_stats, imp = imp_i, rs_info = cond0$rs_info_c,
      sinv_method = sinv_method
    )
    if (is.null(gphi_list[[i]])) {
      return(NULL)
    }
  }
  if (case_b) {
    hmat <- hmat + cond0$logr
  }
  hmax <- apply(hmat, 1L, max)
  if (any(!is.finite(hmax))) {
    return(NULL)
  }
  loglik_j <- hmax + log(rowSums(exp(hmat - hmax)))
  pw <- exp(hmat - loglik_j) # posterior node weights (rows sum to 1)

  # numerical Jacobian of the per-node ingredient vectors phi_c
  # (central differences; no cluster loops); in Case B, the log
  # density-ratio corrections logr also depend on theta -- their
  # derivatives ride along in the same loop
  x0 <- lav_model_get_parameters(lavmodel)
  npar <- length(x0)
  cond_of <- function(x) {
    lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
    cond <- lav_mvn_cl_rs_cond(lavmodel = lavmodel2, rs_info = rs_info,
                               exo_b = rs_stats$exo.b)
    if (is.null(cond)) {
      return(NULL)
    }
    list(
      phi = lapply(cond$imp, function(m) {
        if (is.null(m)) NULL else lav_mvn_cl_rs_phi_vec(m)
      }),
      logr = cond$logr
    )
  }
  nphi <- length(lav_mvn_cl_rs_phi_vec(cond0$imp[[1]]))
  jac <- lapply(seq_len(nnodes), function(i) matrix(0, nphi, npar))
  dlr_sc <- if (case_b) matrix(0, nclusters, npar) else NULL
  for (k in seq_len(npar)) {
    h_k <- h * max(1, abs(x0[k]))
    x_right <- x_left <- x0
    x_right[k] <- x0[k] + h_k
    x_left[k] <- x0[k] - h_k
    cond_r <- cond_of(x_right)
    cond_l <- cond_of(x_left)
    if (is.null(cond_r) || is.null(cond_l)) {
      return(NULL)
    }
    for (i in seq_len(nnodes)) {
      if (is.null(cond_r$phi[[i]]) || is.null(cond_l$phi[[i]])) {
        return(NULL)
      }
      jac[[i]][, k] <- (cond_r$phi[[i]] - cond_l$phi[[i]]) / (2 * h_k)
    }
    if (case_b) {
      if (is.null(cond_r$logr) || is.null(cond_l$logr)) {
        return(NULL)
      }
      dlr_k <- (cond_r$logr - cond_l$logr) / (2 * h_k)
      dlr_sc[, k] <- rowSums(pw * dlr_k)
    }
  }

  if (per_cluster) {
    # per-cluster scores, natural loglikelihood scale
    sc <- matrix(0, nclusters, npar)
    for (i in seq_len(nnodes)) {
      sc <- sc + (pw[, i] * (-0.5) * gphi_list[[i]]) %*% jac[[i]]
    }
    if (case_b) {
      sc <- sc + dlr_sc
    }
    sc
  } else {
    # d(-2 loglik)/d theta
    dx <- numeric(npar)
    for (i in seq_len(nnodes)) {
      dx <- dx +
        as.numeric(colSums(pw[, i] * gphi_list[[i]]) %*% jac[[i]])
    }
    if (case_b) {
      dx <- dx - 2 * colSums(dlr_sc)
    }
    dx
  }
}


# ---------------------------------------------------------------------
# the EM algorithm
# ---------------------------------------------------------------------
#
# complete data = observed data + the between random vector v_b (in
# factor space, see above) + the 'within slope' variables
# z_w[r] = v_b[z_r] * x_r (one per random-slope path; completed with an
# arbitrary standard-normal marginal).
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
lav_mvn_cl_rs_em_ok <- function(rs = NULL) {
  reason <- character(0L)
  if (isTRUE(rs$info$nl.flag)) {
    reason <- c(reason, gettext(
      "the model contains random slopes for latent covariates"))
  }
  if (isTRUE(rs$info$kz > 0L)) {
    reason <- c(reason, gettext(
      "the model contains between-only endogenous observed variables"))
  }
  if (rs$info$p.flag) {
    reason <- c(reason, gettext(
      "the model contains fixed regression paths from the random-slope
       covariates"))
  }
  if (rs$info$shared.flag) {
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
  # missing data? use the pattern-based version
  if (!is.null(rs_stats$mp)) {
    return(lav_mvn_cl_rs_estep_m(
      rs_stats = rs_stats, imp = imp, rs_info = rs_info,
      sinv_method = sinv_method
    ))
  }
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

# E-step with missing data in the y-block (missing = "ml")
#
# the complete data now include the missing y-values. Per pattern
# (observed set o, missing set m), define
#   B_mo = Sigma_w[m,o] Sigma_w[o,o]^{-1}
#   K_p  = (embedded m-selector) - (embedded B_mo)   [p1 x p1]
#   R_p  = I - K_p
#   C_p  = Schur complement Sigma_w[m,m] - B_mo Sigma_w[o,m], embedded
# then, conditional on (v_b, y_obs), the completed y-vector is
#   y_i = R_p y~_i + K_p (mu_y + Q_i v_b) + eps*,  Var(eps*) = C_p
# (y~_i = the zero-padded observed vector). Because K_p Q_i remains
# affine in x_i, all expected-moment sums again collapse onto the
# per-(cluster, pattern) crossproducts; for a complete pattern
# (K_p = 0, R_p = I, C_p = 0) every term reduces to the complete-data
# E-step. The M-step (the two-group SEM) is unchanged.
lav_mvn_cl_rs_estep_m <- function(rs_stats = NULL, imp = NULL,
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
  mp <- rs_stats$mp

  sigma_w <- imp$sigma.w
  sigma_v <- imp$sigma.v
  mu_y <- imp$mu.y
  nu_map <- imp$mu.y.b
  q0 <- imp$q0
  lmat <- imp$lmat
  z_v_idx <- imp$z.v.idx
  path_zcol <- z_v_idx[path_tab$z.idx]
  path_x <- path_tab$x.idx
  dpv <- diag(pv)

  # between mapping: (y_b, z) = c0 + G v
  gmat <- rbind(
    q0[rs_info$yb.y.idx, , drop = FALSE],
    dpv[z_v_idx, , drop = FALSE]
  )
  c0 <- c(nu_map[rs_info$yb.y.idx], rep(0, q))

  # ---- per-pattern constants ----
  npat <- mp$npatterns
  pc <- vector("list", npat)
  for (p in seq_len(npat)) {
    o <- mp$pattern[[p]]$o
    m <- seq_len(p1)[-o]
    if (length(m) == 0L) {
      m <- integer(0L)
    }
    w_o <- lav_mat_sym_inverse(sigma_w[o, o, drop = FALSE],
      logdet = TRUE, sinv_method = sinv_method
    )
    if (!is.finite(attr(w_o, "logdet"))) {
      lav_msg_stop(gettext(
        "within covariance matrix is not positive definite in the
         E-step."))
    }
    w_mat <- matrix(0, p1, p1)
    w_mat[o, o] <- w_o

    k_mat <- matrix(0, p1, p1)
    ct_mat <- matrix(0, p1, p1)
    if (length(m) > 0L) {
      b_mo <- sigma_w[m, o, drop = FALSE] %*% w_o
      k_mat[cbind(m, m)] <- 1
      k_mat[m, o] <- -b_mo
      ct_mat[m, m] <- sigma_w[m, m, drop = FALSE] -
        b_mo %*% sigma_w[o, m, drop = FALSE]
    }
    r_mat <- diag(p1) - k_mat

    pc[[p]] <- list(
      w = w_mat,
      k = k_mat, r = r_mat, ct = ct_mat,
      g = as.numeric(k_mat %*% mu_y) - nu_map,
      h0 = -r_mat %*% q0, # = K Q0 - Q0 (p1 x pv)
      kl = k_mat %*% lmat, # p1 x npaths
      # constants for the A_j/p_j accumulation
      wl = w_mat %*% lmat,
      wq0 = w_mat %*% q0,
      q0twq0 = crossprod(q0, w_mat %*% q0),
      q0twl = crossprod(q0, w_mat %*% lmat),
      ltwl = crossprod(lmat, w_mat %*% lmat)
    )
  }

  # ---- phase 1: pattern-aware A_j and p_j (P = 0) ----
  a_arr <- array(0, dim = c(pv, pv, nclusters))
  p_all <- matrix(0, nclusters, pv)
  sx_tot <- matrix(0, nclusters, nx)
  sxx_tot <- array(0, dim = c(nx, nx, nclusters))
  ntot_used <- 0

  for (p in seq_len(npat)) {
    cell <- mp$pattern[[p]]
    cst <- pc[[p]]
    j1 <- cell$j1
    for (jj in seq_along(j1)) {
      j <- j1[jj]
      n_jp <- cell$n[jj]
      sy_jp <- cell$sy[jj, ]
      sx_jp <- cell$sx[jj, ]
      sxx_jp <- cell$sxx[[jj]]
      sxy_jp <- cell$sxy[[jj]]

      ntot_used <- ntot_used + n_jp
      sx_tot[j, ] <- sx_tot[j, ] + sx_jp
      sxx_tot[, , j] <- sxx_tot[, , j] + sxx_jp

      se_jp <- sy_jp - n_jp * mu_y
      a_jp <- n_jp * cst$q0twq0
      p_jp <- as.numeric(crossprod(cst$wq0, se_jp))
      cz <- matrix(0, pv, pv)
      for (i in seq_len(npaths)) {
        zcol <- path_zcol[i]
        cz[, zcol] <- cz[, zcol] + cst$q0twl[, i] * sx_jp[path_x[i]]
      }
      a_jp <- a_jp + cz + t(cz)
      for (i in seq_len(npaths)) {
        for (k in seq_len(npaths)) {
          a_jp[path_zcol[i], path_zcol[k]] <-
            a_jp[path_zcol[i], path_zcol[k]] +
            cst$ltwl[i, k] * sxx_jp[path_x[i], path_x[k]]
        }
      }
      sxe_jp <- sxy_jp - outer(sx_jp, mu_y)
      for (i in seq_len(npaths)) {
        zi <- path_zcol[i]
        p_jp[zi] <- p_jp[zi] + sum(cst$wl[, i] * sxe_jp[path_x[i], ])
      }

      a_arr[, , j] <- a_arr[, , j] + a_jp
      p_all[j, ] <- p_all[j, ] + p_jp
    }
  }

  # ---- posterior of v_b, per cluster ----
  mu0 <- matrix(0, nclusters, pv)
  m0_list <- vector("list", nclusters)
  for (j in seq_len(nclusters)) {
    a_j <- matrix(a_arr[, , j], pv, pv)
    if (nexo > 0L) {
      d_j <- imp$mu.v +
        as.numeric(imp$cc %*% (rs_stats$exo.b[j, ] - imp$mu.exo))
    } else {
      d_j <- imp$mu.v
    }
    z_j <- dpv + a_j %*% sigma_v
    sigma0_j <- sigma_v %*% solve(z_j)
    sigma0_j <- (sigma0_j + t(sigma0_j)) / 2
    p_tilde <- p_all[j, ] - as.numeric(a_j %*% d_j)
    mu0_j <- d_j + as.numeric(sigma0_j %*% p_tilde)
    mu0[j, ] <- mu0_j
    m0_list[[j]] <- sigma0_j + tcrossprod(mu0_j) # E[v v']
  }

  # ---- accumulators ----
  sum1_w <- numeric(p1 + q)
  sum2_w <- matrix(0, p1 + q, p1 + q)
  sum1_b <- numeric(pb + q + nexo)
  sum2_b <- matrix(0, pb + q + nexo, pb + q + nexo)

  y_idx1 <- seq_len(p1)
  z_idx1 <- p1 + seq_len(q)
  b_idx2 <- seq_len(pb + q)
  w_idx2 <- pb + q + seq_len(nexo)

  # ---- phase 2: within y-moments, per (pattern, cluster) cell ----
  # completed within part: y_w = c_i + H_i v + eps*, with
  #   c_i = R y~_i + g,  H_i = H0 + sum_r x_ir (K L_r) e_{z_r}'
  for (p in seq_len(npat)) {
    cell <- mp$pattern[[p]]
    cst <- pc[[p]]
    j1 <- cell$j1
    for (jj in seq_along(j1)) {
      j <- j1[jj]
      n_jp <- cell$n[jj]
      sy_jp <- cell$sy[jj, ]
      sx_jp <- cell$sx[jj, ]
      sxx_jp <- cell$sxx[[jj]]
      sxy_jp <- cell$sxy[[jj]]
      syy_jp <- cell$syy[[jj]]
      mu0_j <- mu0[j, ]
      m0_j <- m0_list[[j]]

      # sum_i c_i and sum_i c_i c_i'
      rsy <- as.numeric(cst$r %*% sy_jp)
      sc <- rsy + n_jp * cst$g
      cc_mat <- cst$r %*% syy_jp %*% t(cst$r) +
        outer(rsy, cst$g) + outer(cst$g, rsy) +
        n_jp * tcrossprod(cst$g)

      # sum_i H_i (p1 x pv)
      sh <- n_jp * cst$h0
      for (r in seq_len(q)) {
        zc <- path_zcol[r]
        sh[, zc] <- sh[, zc] + sx_jp[path_x[r]] * cst$kl[, r]
      }

      # sum_i x_ir c_i, per slope r (p1 vectors)
      cx <- vector("list", q)
      for (r in seq_len(q)) {
        xr <- path_x[r]
        cx[[r]] <- as.numeric(cst$r %*% sxy_jp[xr, ]) +
          sx_jp[xr] * cst$g
      }

      # sum_i c_i (H_i mu0)'
      h0mu <- as.numeric(cst$h0 %*% mu0_j)
      ch <- outer(sc, h0mu)
      for (r in seq_len(q)) {
        ch <- ch + mu0_j[path_zcol[r]] * outer(cx[[r]], cst$kl[, r])
      }

      # sum_i H_i M0 H_i'
      h0m0 <- cst$h0 %*% m0_j # p1 x pv
      hmh <- n_jp * h0m0 %*% t(cst$h0)
      for (r in seq_len(q)) {
        zc_r <- path_zcol[r]
        h0m_r <- h0m0[, zc_r]
        hmh <- hmh + sx_jp[path_x[r]] *
          (outer(h0m_r, cst$kl[, r]) + outer(cst$kl[, r], h0m_r))
        for (s in seq_len(q)) {
          hmh <- hmh + sxx_jp[path_x[r], path_x[s]] *
            m0_j[zc_r, path_zcol[s]] *
            outer(cst$kl[, r], cst$kl[, s])
        }
      }

      sum1_w[y_idx1] <- sum1_w[y_idx1] + sc + as.numeric(sh %*% mu0_j)
      sum2_w[y_idx1, y_idx1] <- sum2_w[y_idx1, y_idx1] +
        cc_mat + ch + t(ch) + hmh + n_jp * cst$ct

      # cross y_w x z_w[r]
      for (r in seq_len(q)) {
        zc_r <- path_zcol[r]
        xr <- path_x[r]
        cross_r <- mu0_j[zc_r] * cx[[r]] +
          sx_jp[xr] * h0m0[, zc_r]
        for (s in seq_len(q)) {
          cross_r <- cross_r + sxx_jp[xr, path_x[s]] *
            m0_j[path_zcol[s], zc_r] * cst$kl[, s]
        }
        sum2_w[y_idx1, z_idx1[r]] <- sum2_w[y_idx1, z_idx1[r]] + cross_r
        sum2_w[z_idx1[r], y_idx1] <- sum2_w[z_idx1[r], y_idx1] + cross_r
      }
    }
  }

  # ---- z_w moments and between moments, per cluster ----
  for (j in seq_len(nclusters)) {
    mu0_j <- mu0[j, ]
    m0_j <- m0_list[[j]]

    for (r in seq_len(q)) {
      zc_r <- z_v_idx[r]
      xr <- path_x[r]
      sum1_w[z_idx1[r]] <- sum1_w[z_idx1[r]] +
        mu0_j[zc_r] * sx_tot[j, xr]
      for (s in seq_len(q)) {
        sum2_w[z_idx1[r], z_idx1[s]] <- sum2_w[z_idx1[r], z_idx1[s]] +
          m0_j[zc_r, z_v_idx[s]] * sxx_tot[xr, path_x[s], j]
      }
    }

    # between complete-data moments: (y_b, z, w)
    b_j <- c0 + as.numeric(gmat %*% mu0_j)
    bb_j <- gmat %*% m0_j %*% t(gmat) +
      outer(c0, b_j) + outer(b_j, c0) - tcrossprod(c0)
    sum1_b[b_idx2] <- sum1_b[b_idx2] + b_j
    sum2_b[b_idx2, b_idx2] <- sum2_b[b_idx2, b_idx2] + bb_j
    if (nexo > 0L) {
      w_j <- rs_stats$exo.b[j, ]
      sum1_b[w_idx2] <- sum1_b[w_idx2] + w_j
      sum2_b[b_idx2, w_idx2] <- sum2_b[b_idx2, w_idx2] +
        outer(b_j, w_j)
      sum2_b[w_idx2, b_idx2] <- sum2_b[w_idx2, b_idx2] +
        outer(w_j, b_j)
      sum2_b[w_idx2, w_idx2] <- sum2_b[w_idx2, w_idx2] +
        tcrossprod(w_j)
    }
  }

  r1 <- sum1_w / ntot_used
  t1 <- sum2_w / ntot_used - tcrossprod(r1)
  r2 <- sum1_b / nclusters
  t2 <- sum2_b / nclusters - tcrossprod(r2)

  list(r1 = r1, t1 = t1, r2 = r2, t2 = t2, ntotal.used = ntot_used)
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
  # completion)
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
# (mirrors lav_mvn_cl_em_h0; single group)
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
  reason <- lav_mvn_cl_rs_em_ok(rs)
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
  ntot_used <- rs$stats$ntotal.used
  if (is.null(ntot_used)) {
    ntot_used <- rs$stats$ntotal
  }
  nobs_list <- list(ntot_used, rs$stats$nclusters)

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
      fit.by.level = FALSE,
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
# analytic gradient (and per-cluster scores)
# ---------------------------------------------------------------------
#
# the gradient of the -2 loglikelihood is computed in two (chained)
# layers:
#
# 1) an ADJOINT pass (lav_mvn_cl_rs_dfphi): analytic derivatives of
#    F_j = sum_P n_jP ln|Sigma_w[o,o]| + ln|Z_j| + q~_j
#          - p~_j' Sigma_v Z_j^{-1} p~_j
#    with respect to the model-implied ingredients
#    phi = (Sigma_w, mu_y, P, lmat, q0, Sigma_v, mu_v, cc, mu_exo).
#    With u = Sigma_v Z^{-1} p~ (so mu0 = d + u, Sigma0 = Sigma_v
#    Z^{-1}), the derivatives wrt the per-cluster workhorses are
#    remarkably clean (Fisher's identity at work):
#       dF/dA      = Sigma0 + mu0 mu0' = M0
#       dF/dp      = -2 mu0
#       dF/dq      = 1
#       dF/dd      = -2 (p - A mu0)
#       dF/dSigma_v = sym(Z^{-1} A) - sym( (p~ - A u) (Z^{-1} p~)' )
#    (all inverse-free: valid for singular/indefinite Sigma_v);
#    chaining through A_j, p_j, q_j collapses -- once again -- onto
#    the per-(cluster, pattern) crossproducts: ONE pass over the data
#    cells, cost comparable to a single loglikelihood evaluation
#
# 2) the (small) Jacobian d phi / d theta of the implied map,
#    computed by central differences of lav_mvn_cl_rs_implied():
#    2 x npar evaluations of a few small matrix products -- no
#    cluster loops, negligible cost
#
# the same adjoint pass, kept per cluster, yields the per-cluster
# scores (for the robust/sandwich standard errors) at no extra cost

# flatten the theta-dependent implied ingredients into one vector
# (the adjoint pass uses the same layout)
lav_mvn_cl_rs_phi_vec <- function(imp) {
  c(
    as.numeric(imp$sigma.w), imp$mu.y, as.numeric(imp$P),
    as.numeric(imp$lmat), as.numeric(imp$q0), as.numeric(imp$sigma.v),
    imp$mu.v, as.numeric(imp$cc), imp$mu.exo,
    imp$mu.z, as.numeric(imp$gmat), as.numeric(imp$sigma.z)
  )
}

# adjoint pass: returns a (nclusters x nphi) matrix with the
# per-cluster derivatives of F_j wrt phi (NULL if something fails)
lav_mvn_cl_rs_dfphi <- function(rs_stats = NULL, imp = NULL,
                                rs_info = NULL,
                                sinv_method = "eigen") {
  p1 <- rs_info$p1
  nx <- rs_info$nx
  pv <- imp$pv
  path_tab <- rs_info$path.tab
  npaths <- nrow(path_tab)
  z_v_idx <- imp$z.v.idx
  path_zcol <- z_v_idx[path_tab$z.idx]
  path_x <- path_tab$x.idx
  nclusters <- rs_stats$nclusters
  nexo <- rs_info$nexo.b
  exo_b <- rs_stats$exo.b

  sigma_v <- imp$sigma.v
  mu_v <- imp$mu.v
  cc_mat <- imp$cc
  mu_exo <- imp$mu.exo
  mu_y <- imp$mu.y
  pmat <- imp$P
  pmat_t <- t(pmat)
  lmat <- imp$lmat
  q0 <- imp$q0
  dpv <- diag(pv)
  mu_y_op <- tcrossprod(mu_y)

  # flattened phi layout (must match lav_mvn_cl_rs_phi_vec)
  kz <- rs_info$kz
  sizes <- c(p1 * p1, p1, p1 * nx, p1 * npaths, p1 * pv,
             pv * pv, pv, pv * nexo, nexo,
             kz, kz * pv, kz * kz)
  off <- cumsum(c(0, sizes))
  i_sw <- off[1] + seq_len(sizes[1])
  i_muy <- off[2] + seq_len(sizes[2])
  i_p <- off[3] + seq_len(sizes[3])
  i_lm <- off[4] + seq_len(sizes[4])
  i_q0 <- off[5] + seq_len(sizes[5])
  i_sv <- off[6] + seq_len(sizes[6])
  i_mv <- off[7] + seq_len(sizes[7])
  i_cc <- off[8] + seq_len(sizes[8])
  i_me <- off[9] + seq_len(sizes[9])
  i_mz <- off[10] + seq_len(sizes[10])
  i_g <- off[11] + seq_len(sizes[11])
  i_sz <- off[12] + seq_len(sizes[12])
  nphi <- off[13]

  gphi <- matrix(0, nclusters, nphi)

  # z-outcome constants (the extra observation block)
  if (kz > 0L) {
    gmat_z <- imp$gmat
    mu_zb <- imp$mu.z
    zb_all <- rs_stats$zb
    wzb <- lav_mat_sym_inverse(imp$sigma.z,
      logdet = TRUE, sinv_method = sinv_method
    )
    if (!is.finite(attr(wzb, "logdet"))) {
      return(NULL)
    }
    wzg <- wzb %*% gmat_z # kz x pv
    a_zb <- crossprod(gmat_z, wzg) # pv x pv
  }

  # ---- z-block adjoints (mu_z, G, Sigma_z); the z-block is one
  #      extra 'cell' with constant loading G and weight Wz ----
  add_zcell <- function(j, g_a, g_p) {
    e_z <- zb_all[j, ] - mu_zb
    wze <- as.numeric(wzb %*% e_z)
    gga <- gmat_z %*% g_a # kz x pv
    ggp <- as.numeric(gmat_z %*% g_p) # kz
    # G adjoint (cfr. the q0 adjoint with n = 1, e~ = e_z)
    gg <- 2 * (wzb %*% gga) + outer(wze, g_p)
    # mu_z adjoint (cfr. the mu_y adjoint)
    gmz <- -as.numeric(wzb %*% (ggp + 2 * e_z))
    # Sigma_z adjoint (via Wz, plus the ln|Sigma_z| term)
    gw_z <- gga %*% t(gmat_z) + outer(ggp, e_z) + tcrossprod(e_z)
    gw_z <- (gw_z + t(gw_z)) / 2
    gsz <- wzb - wzb %*% gw_z %*% wzb
    gphi[j, i_mz] <<- gphi[j, i_mz] + gmz
    gphi[j, i_g] <<- gphi[j, i_g] + as.numeric(gg)
    gphi[j, i_sz] <<- gphi[j, i_sz] + as.numeric(gsz)
    invisible(NULL)
  }

  # ---- posterior pieces + (Sigma_v, mu_v, cc, mu_exo) adjoints ----
  # (returns the g_a/g_p adjoints needed for the cell collapse)
  post_adj <- function(j, a_j, p_j) {
    if (nexo > 0L) {
      d_j <- mu_v + as.numeric(cc_mat %*% (exo_b[j, ] - mu_exo))
    } else {
      d_j <- mu_v
    }
    p_t <- p_j - as.numeric(a_j %*% d_j)
    z_j <- dpv + a_j %*% sigma_v
    z_inv <- solve(z_j)
    zp <- as.numeric(z_inv %*% p_t)
    u <- as.numeric(sigma_v %*% zp)
    mu0 <- d_j + u

    s0 <- sigma_v %*% z_inv
    s0 <- (s0 + t(s0)) / 2
    g_a <- s0 + tcrossprod(mu0)
    g_p <- -2 * mu0
    g_d <- -2 * (p_j - as.numeric(a_j %*% mu0))
    za <- z_inv %*% a_j
    m1 <- outer(p_t - as.numeric(a_j %*% u), zp)
    g_sv <- (za + t(za)) / 2 - (m1 + t(m1)) / 2

    gphi[j, i_sv] <<- gphi[j, i_sv] + as.numeric(g_sv)
    gphi[j, i_mv] <<- gphi[j, i_mv] + g_d
    if (nexo > 0L) {
      gphi[j, i_cc] <<- gphi[j, i_cc] +
        as.numeric(outer(g_d, exo_b[j, ] - mu_exo))
      gphi[j, i_me] <<- gphi[j, i_me] -
        as.numeric(crossprod(cc_mat, g_d))
    }
    list(g_a = g_a, g_p = g_p)
  }

  # ---- collapse of one (cluster, pattern) cell onto the
  #      (Sigma_w, mu_y, P, lmat, q0) adjoints ----
  add_cell <- function(j, w_mat, n_c, sx_c, sxx_c, sy_c, sxy_c, syy_c,
                       se_c, sxe_c, g_a, g_p) {
    # SQ = sum_i Q_i
    sq <- n_c * q0
    for (r in seq_len(npaths)) {
      sq[, path_zcol[r]] <- sq[, path_zcol[r]] +
        sx_c[path_x[r]] * lmat[, r]
    }
    q0ga <- q0 %*% g_a

    # G_W (raw): <G_A, dA/dW> + <g_p, dp/dW> + dq/dW
    gw <- n_c * (q0ga %*% t(q0))
    for (r in seq_len(npaths)) {
      arv <- as.numeric(q0 %*% g_a[, path_zcol[r]])
      gw <- gw + sx_c[path_x[r]] *
        (outer(lmat[, r], arv) + outer(arv, lmat[, r]))
    }
    for (r in seq_len(npaths)) {
      for (s in seq_len(npaths)) {
        gw <- gw + sxx_c[path_x[r], path_x[s]] *
          g_a[path_zcol[r], path_zcol[s]] *
          outer(lmat[, r], lmat[, s])
      }
    }
    q0gp <- as.numeric(q0 %*% g_p)
    gw <- gw + outer(q0gp, se_c)
    for (r in seq_len(npaths)) {
      gw <- gw + g_p[path_zcol[r]] * outer(lmat[, r], sxe_c[path_x[r], ])
    }
    # sum_i e~_i e~_i'
    psx <- as.numeric(pmat %*% sx_c)
    see <- syy_c - outer(sy_c, mu_y) - outer(mu_y, sy_c) +
      n_c * mu_y_op -
      t(sxy_c) %*% pmat_t - pmat %*% sxy_c +
      outer(mu_y, psx) + outer(psx, mu_y) +
      pmat %*% sxx_c %*% pmat_t
    gw <- gw + see
    gw <- (gw + t(gw)) / 2

    # chain W = Sigma_w[o,o]^{-1} (embedded) -> Sigma_w,
    # including the n_c ln|Sigma_w[o,o]| term
    gsw <- n_c * w_mat - w_mat %*% gw %*% w_mat

    # q0 adjoint: sum_i (2 W Q_i G_A + W e~_i g_p')
    gq0 <- 2 * (w_mat %*% (sq %*% g_a)) +
      outer(as.numeric(w_mat %*% se_c), g_p)

    # lmat adjoint
    glm <- matrix(0, p1, npaths)
    for (r in seq_len(npaths)) {
      xq <- sx_c[path_x[r]] * q0
      for (s in seq_len(npaths)) {
        xq[, path_zcol[s]] <- xq[, path_zcol[s]] +
          sxx_c[path_x[r], path_x[s]] * lmat[, s]
      }
      glm[, r] <- 2 * as.numeric(w_mat %*% (xq %*% g_a[, path_zcol[r]])) +
        g_p[path_zcol[r]] * as.numeric(w_mat %*% sxe_c[path_x[r], ])
    }

    # mu_y adjoint: -sum_i (W Q_i g_p + 2 W e~_i)
    gmy <- -as.numeric(w_mat %*% (as.numeric(sq %*% g_p) + 2 * se_c))

    # P adjoint: -[ W sum_i (Q_i g_p) x_i' + 2 W sum_i e~_i x_i' ]
    mm <- outer(q0gp, sx_c)
    for (r in seq_len(npaths)) {
      mm <- mm + g_p[path_zcol[r]] * outer(lmat[, r], sxx_c[path_x[r], ])
    }
    gp_mat <- -(w_mat %*% (mm + 2 * t(sxe_c)))

    gphi[j, i_sw] <<- gphi[j, i_sw] + as.numeric(gsw)
    gphi[j, i_muy] <<- gphi[j, i_muy] + gmy
    if (nx > 0L) {
      gphi[j, i_p] <<- gphi[j, i_p] + as.numeric(gp_mat)
    }
    gphi[j, i_lm] <<- gphi[j, i_lm] + as.numeric(glm)
    gphi[j, i_q0] <<- gphi[j, i_q0] + as.numeric(gq0)
    invisible(NULL)
  }

  if (is.null(rs_stats$mp)) {
    # -------------------- complete data --------------------
    w_mat <- lav_mat_sym_inverse(imp$sigma.w,
      logdet = TRUE, sinv_method = sinv_method
    )
    if (!is.finite(attr(w_mat, "logdet"))) {
      return(NULL)
    }
    wl <- w_mat %*% lmat
    ltwl <- crossprod(lmat, wl)
    wq0 <- w_mat %*% q0
    q0twq0 <- crossprod(q0, wq0)
    q0twl <- crossprod(q0, wl)
    nj_all <- rs_stats$cluster.size

    for (j in seq_len(nclusters)) {
      nj <- nj_all[j]
      sy_j <- rs_stats$sy[j, ]
      sx_j <- rs_stats$sx[j, ]
      sxx_j <- rs_stats$sxx[[j]]
      sxy_j <- rs_stats$sxy[[j]]
      syy_j <- rs_stats$syy[[j]]

      se_j <- sy_j - nj * mu_y - as.numeric(pmat %*% sx_j)
      sxe_j <- sxy_j - outer(sx_j, mu_y) - sxx_j %*% pmat_t

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
      for (i in seq_len(npaths)) {
        zi <- path_zcol[i]
        p_j[zi] <- p_j[zi] + sum(wl[, i] * sxe_j[path_x[i], ])
      }

      if (kz > 0L) {
        a_j <- a_j + a_zb
        p_j <- p_j + as.numeric(crossprod(wzg, zb_all[j, ] - mu_zb))
      }

      adj <- post_adj(j, a_j, p_j)
      if (kz > 0L) {
        add_zcell(j, adj$g_a, adj$g_p)
      }
      add_cell(
        j, w_mat, nj, sx_j, sxx_j, sy_j, sxy_j, syy_j,
        se_j, sxe_j, adj$g_a, adj$g_p
      )
    }
  } else {
    # -------------------- missing data --------------------
    mp <- rs_stats$mp
    npat <- mp$npatterns

    # per-pattern (embedded) weight matrices and constants
    w_list <- vector("list", npat)
    for (p in seq_len(npat)) {
      o <- mp$pattern[[p]]$o
      w_o <- lav_mat_sym_inverse(imp$sigma.w[o, o, drop = FALSE],
        logdet = TRUE, sinv_method = sinv_method
      )
      if (!is.finite(attr(w_o, "logdet"))) {
        return(NULL)
      }
      w_p <- matrix(0, p1, p1)
      w_p[o, o] <- w_o
      w_list[[p]] <- w_p
    }

    # pass 1: per-cluster A_j and p_j
    a_arr <- array(0, dim = c(pv, pv, nclusters))
    p_all <- matrix(0, nclusters, pv)
    for (p in seq_len(npat)) {
      cell <- mp$pattern[[p]]
      w_p <- w_list[[p]]
      wl <- w_p %*% lmat
      ltwl <- crossprod(lmat, wl)
      wq0 <- w_p %*% q0
      q0twq0 <- crossprod(q0, wq0)
      q0twl <- crossprod(q0, wl)
      for (jj in seq_along(cell$j1)) {
        j <- cell$j1[jj]
        n_jp <- cell$n[jj]
        sy_jp <- cell$sy[jj, ]
        sx_jp <- cell$sx[jj, ]
        sxx_jp <- cell$sxx[[jj]]
        sxy_jp <- cell$sxy[[jj]]

        se_jp <- sy_jp - n_jp * mu_y - as.numeric(pmat %*% sx_jp)
        sxe_jp <- sxy_jp - outer(sx_jp, mu_y) - sxx_jp %*% pmat_t

        a_jp <- n_jp * q0twq0
        p_jp <- as.numeric(crossprod(wq0, se_jp))
        cz <- matrix(0, pv, pv)
        for (i in seq_len(npaths)) {
          zcol <- path_zcol[i]
          cz[, zcol] <- cz[, zcol] + q0twl[, i] * sx_jp[path_x[i]]
        }
        a_jp <- a_jp + cz + t(cz)
        for (i in seq_len(npaths)) {
          for (k in seq_len(npaths)) {
            a_jp[path_zcol[i], path_zcol[k]] <-
              a_jp[path_zcol[i], path_zcol[k]] +
              ltwl[i, k] * sxx_jp[path_x[i], path_x[k]]
          }
        }
        for (i in seq_len(npaths)) {
          zi <- path_zcol[i]
          p_jp[zi] <- p_jp[zi] + sum(wl[, i] * sxe_jp[path_x[i], ])
        }
        a_arr[, , j] <- a_arr[, , j] + a_jp
        p_all[j, ] <- p_all[j, ] + p_jp
      }
    }

    # z-outcome contributions (once per cluster)
    if (kz > 0L) {
      for (j in seq_len(nclusters)) {
        a_arr[, , j] <- a_arr[, , j] + a_zb
        p_all[j, ] <- p_all[j, ] +
          as.numeric(crossprod(wzg, zb_all[j, ] - mu_zb))
      }
    }

    # posterior adjoints per cluster
    ga_list <- vector("list", nclusters)
    gp_list <- vector("list", nclusters)
    for (j in seq_len(nclusters)) {
      adj <- post_adj(j, matrix(a_arr[, , j], pv, pv), p_all[j, ])
      ga_list[[j]] <- adj$g_a
      gp_list[[j]] <- adj$g_p
      if (kz > 0L) {
        add_zcell(j, adj$g_a, adj$g_p)
      }
    }

    # pass 2: cell collapse
    for (p in seq_len(npat)) {
      cell <- mp$pattern[[p]]
      w_p <- w_list[[p]]
      for (jj in seq_along(cell$j1)) {
        j <- cell$j1[jj]
        n_jp <- cell$n[jj]
        sy_jp <- cell$sy[jj, ]
        sx_jp <- cell$sx[jj, ]
        sxx_jp <- cell$sxx[[jj]]
        sxy_jp <- cell$sxy[[jj]]
        syy_jp <- cell$syy[[jj]]
        se_jp <- sy_jp - n_jp * mu_y - as.numeric(pmat %*% sx_jp)
        sxe_jp <- sxy_jp - outer(sx_jp, mu_y) - sxx_jp %*% pmat_t
        add_cell(
          j, w_p, n_jp, sx_jp, sxx_jp, sy_jp, sxy_jp, syy_jp,
          se_jp, sxe_jp, ga_list[[j]], gp_list[[j]]
        )
      }
    }
  }

  gphi
}

# gradient of the -2 loglikelihood wrt the free parameters
# (per_cluster = TRUE: the per-cluster scores d loglik_j / d theta,
#  on the natural loglikelihood scale)
lav_mvn_cl_rs_grad <- function(lavmodel = NULL, rs = NULL,
                               sinv_method = "eigen",
                               per_cluster = FALSE,
                               h = 1e-06) {
  # latent covariates: posterior-weighted adjoint pass per quadrature
  # node (numerical fallback on failure, via the callers)
  if (isTRUE(rs$info$nl.flag)) {
    out <- try(lav_mvn_cl_rs_grad_nl(
      lavmodel = lavmodel, rs = rs, sinv_method = sinv_method,
      per_cluster = per_cluster, h = h
    ), silent = TRUE)
    if (inherits(out, "try-error")) {
      return(NULL)
    }
    return(out)
  }
  imp0 <- lav_mvn_cl_rs_implied(lavmodel, rs_info = rs$info)
  if (is.null(imp0)) {
    return(NULL)
  }
  gphi <- lav_mvn_cl_rs_dfphi(
    rs_stats = rs$stats, imp = imp0, rs_info = rs$info,
    sinv_method = sinv_method
  )
  if (is.null(gphi)) {
    return(NULL)
  }

  # the (small) Jacobian of the implied map: central differences,
  # no cluster loops
  x0 <- lav_model_get_parameters(lavmodel)
  npar <- length(x0)
  phi0 <- lav_mvn_cl_rs_phi_vec(imp0)
  jac <- matrix(0, length(phi0), npar)
  for (k in seq_len(npar)) {
    h_k <- h * max(1, abs(x0[k]))
    x_right <- x_left <- x0
    x_right[k] <- x0[k] + h_k
    x_left[k] <- x0[k] - h_k
    imp_r <- lav_mvn_cl_rs_implied(
      lav_model_set_parameters(lavmodel, x = x_right),
      rs_info = rs$info
    )
    imp_l <- lav_mvn_cl_rs_implied(
      lav_model_set_parameters(lavmodel, x = x_left),
      rs_info = rs$info
    )
    if (is.null(imp_r) || is.null(imp_l)) {
      return(NULL)
    }
    jac[, k] <- (lav_mvn_cl_rs_phi_vec(imp_r) -
      lav_mvn_cl_rs_phi_vec(imp_l)) / (2 * h_k)
  }

  if (per_cluster) {
    -0.5 * (gphi %*% jac)
  } else {
    as.numeric(colSums(gphi) %*% jac)
  }
}

# H0 scaling correction factor (cfr. the Mplus MLR output):
#
#   c_H0 = tr( H^{-1} J ) / npar
#
# with H the total *observed* information (minus the second derivative
# of the loglikelihood; via lav_model_hessian() with method =
# "central": central differences of the analytic gradient, relative
# step size) and J = crossprod of the per-cluster scores (the
# first-order information); this is the ingredient for scaled
# (loglikelihood) difference tests
lav_mvn_cl_rs_scaling_h0 <- function(lavobject = NULL, rs = NULL,
                                     h = 1e-04) {
  lavmodel <- lavobject@Model
  sc <- try(lav_mvn_cl_rs_scores(lavmodel = lavmodel, rs = rs),
            silent = TRUE)
  if (inherits(sc, "try-error")) {
    return(as.numeric(NA))
  }
  j_mat <- crossprod(sc)

  # hessian of the objective ( = -2 loglik / (2 N) ): the observed
  # information is hessian * ntotal
  h_obj <- lav_model_hessian(
    lavmodel = lavmodel,
    lavsamplestats = lavobject@SampleStats,
    lavdata = lavobject@Data,
    lavoptions = lavobject@Options,
    lavcache = lavobject@Cache,
    h = h, method = "central"
  )
  if (anyNA(h_obj)) {
    return(as.numeric(NA))
  }
  h_mat <- h_obj * lavobject@SampleStats@ntotal

  npar <- lavmodel@nx.free
  out <- try(sum(diag(solve(h_mat, j_mat))) / npar, silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) {
    return(as.numeric(NA))
  }
  out
}


# ---------------------------------------------------------------------
# per-cluster scores and the first-order (`meat') information
# ---------------------------------------------------------------------

# per-cluster scores: d loglik_j / d theta (natural loglikelihood
# scale); analytic (adjoint pass), with a numerical fallback
# (central differences of the per-cluster loglikelihood vector)
#
# returns a (nclusters x npar) matrix
lav_mvn_cl_rs_scores <- function(lavmodel = NULL, rs = NULL,
                                 h = 1e-05) {
  # analytic scores
  sc <- try(
    lav_mvn_cl_rs_grad(lavmodel = lavmodel, rs = rs,
                       per_cluster = TRUE),
    silent = TRUE
  )
  if (!inherits(sc, "try-error") && !is.null(sc)) {
    return(sc)
  }

  # numerical fallback
  x0 <- lav_model_get_parameters(lavmodel)
  npar <- length(x0)
  nclusters <- rs$stats$nclusters

  ll_cluster <- function(x) {
    lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
    out <- lav_mvn_cl_rs_m2ll(
      lavmodel = lavmodel2, rs = rs,
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

  # latent covariates: EAP predictions (posterior-weighted mixture
  # over the quadrature nodes)
  if (isTRUE(rs_info$nl.flag)) {
    return(lav_mvn_cl_rs_eb_nl(
      lavmodel = lavmodel, lavdata = lavdata, rs = rs, se = se
    ))
  }

  imp <- lav_mvn_cl_rs_implied(lavmodel, rs_info = rs_info)
  if (is.null(imp)) {
    lav_msg_stop(gettext(
      "model-implied matrices could not be computed."))
  }

  core <- lav_mvn_cl_rs_eb_core(
    lavmodel = lavmodel, glist = lavmodel@GLIST, lavdata = lavdata,
    rs_stats = rs_stats, rs_info = rs_info, imp = imp, se = se
  )

  # level-2 scores: the between latent variables (eta part of v_b)
  meta <- imp$meta
  l2 <- core$mu0[, seq_len(meta), drop = FALSE]
  colnames(l2) <- imp$eta.names
  out <- list(l1 = core$l1, l2 = l2, cluster.idx = core$cluster.idx)
  if (se) {
    se2 <- sqrt(pmax(core$var0[, seq_len(meta), drop = FALSE], 0))
    colnames(se2) <- imp$eta.names
    out$se2 <- se2
  }
  out
}

# the (linear-model) EB workhorse: per-cluster posterior of v_b
# (mu0, and the posterior variances var0 = diag(Sigma0)) plus the
# level-1 (within) factor scores, for the model described by
# (glist, imp, rs_info) -- which may be the *conditional* model of
# one quadrature node (latent covariates); the caller assembles
lav_mvn_cl_rs_eb_core <- function(lavmodel = NULL, glist = NULL,
                                  lavdata = NULL,
                                  rs_stats = NULL, rs_info = NULL,
                                  imp = NULL, se = FALSE) {
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

  # between-only endogenous z-outcomes
  kz <- rs_info$kz
  if (kz > 0L) {
    wzb <- lav_mat_sym_inverse(imp$sigma.z)
    wzg <- wzb %*% imp$gmat
    a_zb <- crossprod(imp$gmat, wzg)
    zb_all <- rs_stats$zb
    mu_zb <- imp$mu.z
  }

  # ---- per-cluster posterior of v_b: mu0_j (and posterior sd) ----
  pmat <- imp$P
  pmat_t <- t(pmat)
  mu0 <- matrix(0, nclusters, pv)
  var0 <- if (se) matrix(0, nclusters, pv) else NULL

  if (!is.null(rs_stats$mp)) {
    # missing data: pattern-based accumulation of A_j and p_j
    # (cfr. lav_mvn_cl_rs_loglik_m)
    mp <- rs_stats$mp
    a_arr <- array(0, dim = c(pv, pv, nclusters))
    p_all <- matrix(0, nclusters, pv)
    for (p in seq_len(mp$npatterns)) {
      cell <- mp$pattern[[p]]
      o <- cell$o
      w_o <- lav_mat_sym_inverse(imp$sigma.w[o, o, drop = FALSE],
        logdet = TRUE)
      if (!is.finite(attr(w_o, "logdet"))) {
        lav_msg_stop(gettext(
          "within covariance matrix is not positive definite."))
      }
      w_p <- matrix(0, p1, p1)
      w_p[o, o] <- w_o
      wl_p <- w_p %*% imp$lmat
      wq0_p <- w_p %*% imp$q0
      q0twq0_p <- crossprod(imp$q0, wq0_p)
      q0twl_p <- crossprod(imp$q0, wl_p)
      ltwl_p <- crossprod(imp$lmat, wl_p)
      for (jj in seq_along(cell$j1)) {
        j <- cell$j1[jj]
        n_jp <- cell$n[jj]
        sy_jp <- cell$sy[jj, ]
        sx_jp <- cell$sx[jj, ]
        sxx_jp <- cell$sxx[[jj]]
        sxy_jp <- cell$sxy[[jj]]
        se_jp <- sy_jp - n_jp * mu_y - as.numeric(pmat %*% sx_jp)
        sxe_jp <- sxy_jp - outer(sx_jp, mu_y) - sxx_jp %*% pmat_t
        a_jp <- n_jp * q0twq0_p
        p_jp <- as.numeric(crossprod(wq0_p, se_jp))
        cz <- matrix(0, pv, pv)
        for (i in seq_len(npaths)) {
          zcol <- path_zcol[i]
          cz[, zcol] <- cz[, zcol] + q0twl_p[, i] * sx_jp[path_x[i]]
        }
        a_jp <- a_jp + cz + t(cz)
        for (i in seq_len(npaths)) {
          for (k in seq_len(npaths)) {
            a_jp[path_zcol[i], path_zcol[k]] <-
              a_jp[path_zcol[i], path_zcol[k]] +
              ltwl_p[i, k] * sxx_jp[path_x[i], path_x[k]]
          }
        }
        for (i in seq_len(npaths)) {
          zi <- path_zcol[i]
          p_jp[zi] <- p_jp[zi] + sum(wl_p[, i] * sxe_jp[path_x[i], ])
        }
        a_arr[, , j] <- a_arr[, , j] + a_jp
        p_all[j, ] <- p_all[j, ] + p_jp
      }
    }
    for (j in seq_len(nclusters)) {
      a_j <- matrix(a_arr[, , j], pv, pv)
      if (kz > 0L) {
        a_j <- a_j + a_zb
        p_all[j, ] <- p_all[j, ] +
          as.numeric(crossprod(wzg, zb_all[j, ] - mu_zb))
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
      p_tilde <- p_all[j, ] - as.numeric(a_j %*% d_j)
      mu0[j, ] <- d_j + as.numeric(sigma0_j %*% p_tilde)
      if (se) {
        var0[j, ] <- diag(sigma0_j)
      }
    }
  } else {
  for (j in seq_len(nclusters)) {
    nj <- nj_all[j]
    sy_j <- rs_stats$sy[j, ]
    sx_j <- rs_stats$sx[j, ]
    sxx_j <- rs_stats$sxx[[j]]
    sxy_j <- rs_stats$sxy[[j]]

    # centered sums: e~_i = y_i - mu_y - P x_i (P = the reduced-form
    # *fixed* slopes of y on the rv covariates)
    se_j <- sy_j - nj * mu_y - as.numeric(pmat %*% sx_j)
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
    sxe_j <- sxy_j - outer(sx_j, mu_y) - sxx_j %*% pmat_t
    for (i in seq_len(npaths)) {
      zi <- path_zcol[i]
      p_j[zi] <- p_j[zi] + sum(wl[, i] * sxe_j[path_x[i], ])
    }

    if (kz > 0L) {
      a_j <- a_j + a_zb
      p_j <- p_j + as.numeric(crossprod(wzg, zb_all[j, ] - mu_zb))
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
      var0[j, ] <- diag(sigma0_j)
    }
  }
  } # complete data

  # ---- level-1 scores: within factor scores given (mu0_j, x_ij) ----
  # within model pieces (same construction as lav_mvn_cl_rs_implied;
  # note: from 'glist', which may hold node-injected beta_w cells)
  nmat <- lavmodel@nmat
  mm_w <- 1:nmat[1]
  mlist_w <- glist[mm_w]
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
    if (is.null(rs_stats$mp)) {
      upd <- e_mat %*% (w_mat %*% t(ceta_y))
    } else {
      # missing data: the update uses the observed coordinates only
      # (cases with no observed outcomes get the plain conditional
      # mean, without a residual update)
      upd <- matrix(0, n1, length(real_idx))
      obs <- !is.na(yy)
      case_id <- apply(1L * obs, 1L, paste, collapse = "")
      for (pid in unique(case_id)) {
        ridx <- which(case_id == pid)
        o <- which(obs[ridx[1L], ])
        if (length(o) == 0L) {
          next
        }
        w_o <- lav_mat_sym_inverse(imp$sigma.w[o, o, drop = FALSE])
        upd[ridx, ] <- e_mat[ridx, o, drop = FALSE] %*%
          (w_o %*% t(ceta_y[, o, drop = FALSE]))
      }
    }
    l1 <- m_full[, real_idx, drop = FALSE] + upd
    colnames(l1) <- lv_names_w[real_idx]
  } else {
    l1 <- matrix(0, n1, 0L)
  }

  list(mu0 = mu0, var0 = var0, l1 = l1, cluster.idx = cluster_idx)
}

# EAP predictions for latent-covariate (quadrature) models: the
# posterior of the cluster-level random effects is a *mixture* over
# the quadrature nodes, with the posterior node weights
# pi_jq = exp(h_jq - loglik_j) (the same weights that drive the
# analytic gradient). Per node, the conditional model is
# linear-Gaussian, so the ordinary EB workhorse applies; the mixture
# is then assembled:
#
#   E[z_nl | y_j]   = sum_q pi_jq b_q
#   E[v_lin | y_j]  = sum_q pi_jq mu0_jq
#   Var[.. | y_j]   = sum_q pi_jq (Sigma0 + mu0 mu0') - mean mean'
#
# and analogously for the level-1 (within) factor scores, which are
# linear in the cluster-level effects per node
lav_mvn_cl_rs_eb_nl <- function(lavmodel = NULL, lavdata = NULL,
                                rs = NULL, se = FALSE) {
  rs_info <- rs$info
  rs_stats <- rs$stats

  cond <- lav_mvn_cl_rs_cond(
    lavmodel = lavmodel, rs_info = rs_info, exo_b = rs_stats$exo.b
  )
  if (is.null(cond)) {
    lav_msg_stop(gettext(
      "model-implied matrices could not be computed."))
  }
  nnodes <- cond$nnodes
  nclusters <- rs_stats$nclusters

  # posterior node weights (cfr. lav_mvn_cl_rs_grad_nl)
  hmat <- matrix(-Inf, nclusters, nnodes)
  for (i in seq_len(nnodes)) {
    if (is.null(cond$imp[[i]])) {
      next
    }
    m2_i <- lav_mvn_cl_rs_loglik(
      rs_stats = rs_stats, imp = cond$imp[[i]],
      rs_info = cond$rs_info_c,
      log2pi = FALSE, minus_two = TRUE, per_cluster = TRUE
    )
    if (!is.finite(m2_i)) {
      next
    }
    hmat[, i] <- -0.5 * attr(m2_i, "loglik.cluster") + cond$logw[i]
  }
  if (!is.null(cond$logr)) {
    hmat <- hmat + cond$logr
  }
  hmax <- apply(hmat, 1L, max)
  if (any(!is.finite(hmax))) {
    lav_msg_stop(gettext(
      "the (quadrature) loglikelihood could not be computed for one
       or more clusters."))
  }
  loglik_j <- hmax + log(rowSums(exp(hmat - hmax)))
  pw <- exp(hmat - loglik_j) # nclusters x nnodes; rows sum to 1

  # per-node conditional EB, mixed with the posterior node weights
  beta_mm <- rs_info$nl.beta.mm
  beta_cells <- cbind(rs_info$nl.beta.row, rs_info$nl.beta.col)
  path_tab <- rs_info$path.tab
  nl_path_rows <- which(path_tab$type == "lv")
  path_nl_z <- match(path_tab$z.idx[nl_path_rows], rs_info$nl.z.idx)

  pv_c <- length(cond$lin.v)
  m_lin <- matrix(0, nclusters, pv_c) # E[v_lin | y]
  v_lin <- if (se) matrix(0, nclusters, pv_c) else NULL # E[v^2] part
  l1_mix <- NULL
  glist_i <- lavmodel@GLIST
  for (i in seq_len(nnodes)) {
    if (is.null(cond$imp[[i]])) {
      # a failed node: its posterior weights are (numerically) zero
      next
    }
    b_i <- cond$bmat[, i]
    glist_i[[beta_mm]][beta_cells] <- b_i[path_nl_z]
    core <- lav_mvn_cl_rs_eb_core(
      lavmodel = lavmodel, glist = glist_i, lavdata = lavdata,
      rs_stats = rs_stats, rs_info = cond$rs_info_c,
      imp = cond$imp[[i]], se = se
    )
    m_lin <- m_lin + pw[, i] * core$mu0
    if (se) {
      v_lin <- v_lin + pw[, i] * (core$var0 + core$mu0^2)
    }
    if (is.null(l1_mix)) {
      l1_mix <- matrix(0, nrow(core$l1), ncol(core$l1))
      colnames(l1_mix) <- colnames(core$l1)
      cluster_idx <- core$cluster.idx
    }
    if (ncol(l1_mix) > 0L) {
      l1_mix <- l1_mix + pw[cluster_idx, i] * core$l1
    }
  }
  if (se) {
    v_lin <- v_lin - m_lin^2
  }

  # EAP of the nonlinear slopes: moments of the discrete posterior
  q_nl <- length(cond$nl.v)
  bt <- t(cond$bmat) # nnodes x q_nl
  m_nl <- pw %*% bt
  v_nl <- if (se) pw %*% bt^2 - m_nl^2 else NULL

  # assemble the level-2 output in the full eta ordering
  meta <- cond$meta
  l2 <- matrix(0, nclusters, meta)
  colnames(l2) <- cond$eta.names
  se2 <- if (se) l2 else NULL
  for (m in seq_len(meta)) {
    lin_pos <- match(m, cond$lin.v)
    if (!is.na(lin_pos)) {
      l2[, m] <- m_lin[, lin_pos]
      if (se) {
        se2[, m] <- sqrt(pmax(v_lin[, lin_pos], 0))
      }
    } else {
      nl_pos <- match(m, cond$nl.v)
      l2[, m] <- m_nl[, nl_pos]
      if (se) {
        se2[, m] <- sqrt(pmax(v_nl[, nl_pos], 0))
      }
    }
  }

  out <- list(l1 = l1_mix, l2 = l2, cluster.idx = cluster_idx)
  if (se) {
    out$se2 <- se2
  }
  out
}
