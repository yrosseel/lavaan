# James-Stein type (conditional-expectation) estimation
#
# Burghgraeve, E., De Neve, J., & Rosseel, Y. (2021). Estimating structural
# equation models using James-Stein type shrinkage estimators. Psychometrika,
# 86(1), 96-130.
#
# In every (measurement or structural) equation, the latent variables are
# replaced by an estimate of their conditional expectation given observed
# proxies: the scaling indicators (estimator = "JS"), or reliability-weighted
# aggregates of the indicators (estimator = "JSA"); the equation is then
# estimated by least squares. Everything is computed from the sample moments
# ('attenuation-corrected' least squares): regressing on the estimated
# conditional expectations is algebraically identical to correcting the
# normal equations for the measurement-error variance of the proxies.
#
# The equation decomposition (dependent variables, scaling-indicator
# substitution) is shared with the (MI)IV estimator (lav_model_find_iv);
# the two estimators differ only in how each equation is estimated (2SLS
# with instruments versus least squares on conditional expectations).
# The undirected (variance/covariance) parameters are estimated by the
# same second stage as the IV estimator (lav_sem_miiv_varcov).
#
# YR/EB 08 Feb 2023: first version (lav_cfa_jamesstein.R), cfa only
# YR    18 Jul 2026: rewritten in moment form: structural models, observed
#                    covariates, mean structure; stage 2 via the IV
#                    machinery; closed-form (alphamax) aggregation weights

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_js_internal <- function(lavmodel = NULL, lavh1 = NULL,
                                lavsamplestats = NULL,
                                lavpartable = NULL,
                                lavdata = NULL, lavoptions = NULL) {
  lavpta <- lav_pt_attributes(lavpartable)
  lavpartable <- lav_pt_set_cache(lavpartable, lavpta)

  # aggregated (JSA) or scaling (JS)?
  aggregated <- identical(lavoptions$estimator, "JSA")

  # model-class checks
  lav_sem_js_check(
    lavmodel = lavmodel, lavpartable = lavpartable, lavpta = lavpta,
    lavdata = lavdata, lavoptions = lavoptions, aggregated = aggregated
  )

  # equations: reuse the (MI)IV decomposition (each dependent variable, with
  # the latent variables replaced by their scaling indicators); the
  # instrument sets (eq$iv/eq$miiv) are ignored here
  eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)

  # directed versus undirected (free) parameters
  undirected_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  directed_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~", "~1"))
  free_directed_idx <- unique(lavpartable$free[directed_idx])
  free_undirected_idx <- unique(lavpartable$free[undirected_idx])

  # initial parameter vector
  x <- lav_model_get_parameters(lavmodel)

  ########################################
  # first stage: directed parameters     #
  ########################################
  theta1 <- numeric(0L)
  if (length(free_directed_idx) > 0L) {
    theta1 <- lav_sem_js_stage1_samp(
      x = NULL, samplestats = FALSE, eqs = eqs,
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats, lavh1 = lavh1,
      lavoptions = lavoptions,
      free_directed_idx = free_directed_idx,
      aggregated = aggregated
    )
    # update equations
    eqs <- attr(theta1, "eqs")
    theta1 <- as.numeric(theta1) # drop attributes
    # store theta1 elements in x
    x[free_directed_idx] <- theta1
  }

  #######################################
  # second stage: undirected parameters #
  #######################################
  js_varcov_method <- toupper(lavoptions$estimator.args[["js_varcov_method"]])
  if (length(js_varcov_method) == 0L) {
    js_varcov_method <- "RLS"
  }
  if (length(free_undirected_idx) > 0L && js_varcov_method != "NONE") {
    theta2 <- lav_sem_miiv_varcov(
      x = theta1,
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats,
      lavh1 = lavh1, free_directed_idx = free_directed_idx,
      free_undirected_idx = free_undirected_idx,
      iv_varcov_method = js_varcov_method
    )
  } else {
    theta2 <- rep(as.numeric(NA), length(free_undirected_idx))
  }
  # store theta2 elements in x
  x[free_undirected_idx] <- theta2

  # mean structure: by default, re-estimate all free mean parameters
  # (observed intercepts and latent means) jointly by GLS, exactly as the IV
  # estimator does; this also provides the means of the exogenous latent
  # variables (which have no equation of their own)
  js_mean_structure <- tolower(lavoptions$estimator.args[["js_mean_structure"]])
  if (length(js_mean_structure) == 0L) {
    js_mean_structure <- "wls"
  }
  if (lavmodel@meanstructure && identical(js_mean_structure, "wls")) {
    free_mean_idx <- unique(lavpartable$free[lavpartable$op == "~1" &
      lavpartable$free > 0L & !duplicated(lavpartable$free)])
    if (length(free_mean_idx) > 0L) {
      x[free_mean_idx] <- lav_sem_miiv_mean_wls(
        lavmodel = lavmodel, lavsamplestats = lavsamplestats, x = x,
        free_mean_idx = free_mean_idx
      )
    }
  }

  # apply bounds (if any)
  if (!is.null(lavpartable$lower)) {
    lower_x <- lavpartable$lower[lavpartable$free > 0 &
      !duplicated(lavpartable$free)]
    too_small_idx <- which(x < lower_x)
    if (length(too_small_idx) > 0L) {
      x[too_small_idx] <- lower_x[too_small_idx]
    }
  }
  if (!is.null(lavpartable$upper)) {
    upper_x <- lavpartable$upper[lavpartable$free > 0 &
      !duplicated(lavpartable$free)]
    too_large_idx <- which(x > upper_x)
    if (length(too_large_idx) > 0L) {
      x[too_large_idx] <- upper_x[too_large_idx]
    }
  }

  attr(x, "eqs") <- eqs
  x
}


# model-class checks for the JS/JSA estimators; stop with a clear message
# for anything the (phase-1) implementation does not cover
lav_sem_js_check <- function(lavmodel = NULL, lavpartable = NULL,
                             lavpta = NULL, lavdata = NULL,
                             lavoptions = NULL, aggregated = FALSE) {
  label <- if (aggregated) "JSA" else "JS"

  stopifnot(lavdata@nlevels == 1L)
  if (lavmodel@nblocks > 1L) {
    lav_msg_stop(gettextf(
      "estimator %s does not support multiple groups (yet).", label))
  }
  if (lavmodel@categorical) {
    lav_msg_stop(gettextf(
      "estimator %s does not support categorical data (yet).", label))
  }
  if (lavoptions$std.lv) {
    lav_msg_stop(gettextf(
      "estimator %s requires marker/scaling indicators;
       std.lv = TRUE is not supported.", label))
  }
  if (lavmodel@conditional.x) {
    lav_msg_stop(gettextf(
      "estimator %s does not support conditional.x = TRUE (yet).", label))
  }
  if (lavmodel@correlation) {
    lav_msg_stop(gettextf(
      "estimator %s does not support correlation structures (yet).", label))
  }
  if (isTRUE(lavmodel@composites)) {
    lav_msg_stop(gettextf(
      "estimator %s does not support composites (yet).", label))
  }
  if (lavmodel@nefa > 0L) {
    lav_msg_stop(gettextf(
      "estimator %s does not support efa blocks (yet).", label))
  }
  if (!is.null(lavdata@weights[[1L]])) {
    lav_msg_stop(gettextf(
      "estimator %s does not support sampling weights (yet).", label))
  }

  b <- 1L
  ov_names <- lavpta$vnames$ov[[b]]
  lv_names <- lavpta$vnames$lv.regular[[b]]

  # higher-order factors: a latent variable measured by latent variables
  if (any(lavpartable$op == "=~" & lavpartable$rhs %in% lv_names)) {
    lav_msg_stop(gettextf(
      "estimator %s does not support higher-order factors (yet).", label))
  }

  # the structural part must be recursive: the conditional expectations are
  # only valid regressors if the disturbance of each equation is
  # uncorrelated with (functions of) its right-hand-side variables
  if (isFALSE(lavmodel@modprop$acyclic[b])) {
    lav_msg_stop(gettextf(
      "estimator %s requires a recursive structural model
       (no feedback loops).", label))
  }

  # active covariance rows (free, or fixed to a nonzero value)
  cov_idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs != lavpartable$rhs &
    (lavpartable$free > 0L |
      (!is.na(lavpartable$ustart) & lavpartable$ustart != 0)))

  # no disturbance of a dependent variable may be correlated with one of its
  # regressors, or with any ancestor of its regressors (the conditional
  # expectations condition on proxies of these variables)
  reg_idx <- which(lavpartable$op == "~")
  if (length(reg_idx) > 0L && length(cov_idx) > 0L) {
    # ancestor sets via the regression graph (edge: rhs -> lhs)
    parents <- split(
      lavpartable$rhs[reg_idx],
      factor(lavpartable$lhs[reg_idx],
             levels = unique(lavpartable$lhs[reg_idx]))
    )
    ancestors <- function(vars) {
      anc <- character(0L)
      todo <- vars
      while (length(todo) > 0L) {
        v <- todo[1L]
        todo <- todo[-1L]
        if (v %in% anc) {
          next
        }
        anc <- c(anc, v)
        todo <- c(todo, parents[[v]])
      }
      anc
    }
    for (e in unique(lavpartable$lhs[reg_idx])) {
      rhs_e <- lavpartable$rhs[reg_idx][lavpartable$lhs[reg_idx] == e]
      anc_e <- ancestors(rhs_e)
      bad <- which((lavpartable$lhs[cov_idx] == e &
                      lavpartable$rhs[cov_idx] %in% anc_e) |
                   (lavpartable$rhs[cov_idx] == e &
                      lavpartable$lhs[cov_idx] %in% anc_e))
      if (length(bad) > 0L) {
        lav_msg_stop(gettextf(
          "estimator %1$s is not consistent when the disturbance of %2$s is
           correlated with (an ancestor of) one of its regressors: %3$s.",
          label, e,
          lav_msg_view(unique(c(lavpartable$lhs[cov_idx][bad],
                                lavpartable$rhs[cov_idx][bad])), "none")))
      }
    }
  }

  # residual covariances among the observed variables
  markers <- lavpta$vnames$lv.marker[[b]]
  ind_names <- unique(lavpartable$rhs[lavpartable$op == "=~"])
  res_names <- unique(c(ind_names, lavpta$vnames$eqs.y[[b]]))
  ov_cov_idx <- cov_idx[lavpartable$lhs[cov_idx] %in% res_names &
    lavpartable$rhs[cov_idx] %in% res_names]
  if (length(ov_cov_idx) > 0L) {
    # a residual covariance involving a scaling indicator invalidates the
    # conditional expectations (and the reliability estimates)
    marker_cov_idx <- ov_cov_idx[lavpartable$lhs[ov_cov_idx] %in% markers |
      lavpartable$rhs[ov_cov_idx] %in% markers]
    if (length(marker_cov_idx) > 0L) {
      lav_msg_stop(gettextf(
        "estimator %s does not support residual covariances that involve a
         marker/scaling indicator.", label))
    }
    if (aggregated) {
      # any correlated indicator residual contaminates the aggregates
      lav_msg_stop(gettext(
        "estimator JSA does not support residual covariances among the
         indicators (yet); use estimator JS instead."))
    }
    # remaining case (JS): the equations remain valid, but the (Spearman)
    # reliability estimates may be biased if the correlated residuals belong
    # to indicators of the same factor
    lav_msg_warn(gettext(
      "[JS] residual covariances among indicators may bias the (Spearman)
       reliability estimates of the scaling indicators."))
  }

  invisible(TRUE)
}


# residual (measurement-error) variances of the observed indicators, needed
# for the shrinkage (and, for JSA, the aggregation weights); by default the
# Spearman/communality estimator applied per factor (Burghgraeve et al.,
# 2021, Appendix B)
lav_sem_js_theta <- function(s = NULL, ind_pat = NULL, ov_names = NULL,
                             lv_names = NULL, lavoptions = NULL,
                             aggregated = FALSE) {
  label <- if (aggregated) "JSA" else "JS"
  nvar <- ncol(s)
  ea <- lavoptions$estimator.args
  js_theta <- tolower(ea[["js_theta"]])
  if (length(js_theta) == 0L) {
    js_theta <- "spearman"
  }
  js_theta_bounds <- tolower(ea[["js_theta_bounds"]])
  if (length(js_theta_bounds) == 0L) {
    js_theta_bounds <- "wide"
  }

  theta <- rep(as.numeric(NA), nvar)
  if (js_theta == "user") {
    values <- ea[["js_theta_values"]]
    ind_idx <- which(rowSums(ind_pat != 0) > 0L)
    missing_idx <- ind_idx[!ov_names[ind_idx] %in% names(values)]
    if (length(missing_idx) > 0L) {
      lav_msg_stop(gettextf(
        "js_theta_values should provide a residual variance for every
         indicator; missing: %s.",
        lav_msg_view(ov_names[missing_idx], "none")))
    }
    theta[ind_idx] <- as.numeric(values[ov_names[ind_idx]])
  } else {
    # spearman, per factor
    for (f in seq_len(ncol(ind_pat))) {
      idx <- which(ind_pat[, f] != 0)
      if (length(idx) < 3L) {
        lav_msg_stop(gettextf(
          "estimator %1$s needs at least 3 indicators per factor to estimate
           the reliability of the scaling indicators; factor %2$s has only
           %3$d. Alternatively, provide the residual variances directly via
           estimator.args = list(js_theta = \"user\", js_theta_values = ...).",
          label, lv_names[f], length(idx)))
      }
      theta[idx] <- lav_cfa_theta_spearman(s[idx, idx, drop = FALSE],
        bounds = js_theta_bounds
      )
    }
  }

  # non-indicator variables (observed covariates, endogenous observed
  # variables) are error-free proxies of themselves
  theta[is.na(theta)] <- 0

  # keep the variances within [0, var(y)]
  diag_s <- diag(s)
  theta[theta < 0] <- 0
  too_large_idx <- which(theta > diag_s)
  if (length(too_large_idx) > 0L) {
    theta[too_large_idx] <- diag_s[too_large_idx]
  }

  theta
}


# reliability-maximizing (alphamax; Bentler, 1968) aggregation weights: the
# maximizer of (w's w - w'diag(theta) w) / (w's w) is the generalized
# eigenvector of (diag(theta), s) with the smallest eigenvalue; normalized
# so that the weights sum to one
lav_sem_js_weights <- function(s = NULL, theta = NULL) {
  p <- ncol(s)
  if (p == 1L) {
    return(1)
  }
  equal_w <- rep(1 / p, p)
  r_mat <- try(chol(s), silent = TRUE)
  if (inherits(r_mat, "try-error")) {
    return(equal_w)
  }
  r_inv <- backsolve(r_mat, diag(p))
  m_sym <- crossprod(r_inv, theta * r_inv) # t(Ri) %*% diag(theta) %*% Ri
  ev <- eigen((m_sym + t(m_sym)) / 2, symmetric = TRUE)
  w <- drop(r_inv %*% ev$vectors[, p]) # smallest eigenvalue
  sw <- sum(w)
  if (!is.finite(sw) || abs(sw) < .Machine$double.eps^0.5) {
    return(equal_w)
  }
  w / sw
}


# stage 1: estimate all directed (loading/regression/intercept) parameters,
# equation by equation, by least squares on the estimated conditional
# expectations -- computed directly from the sample moments.
# if samplestats = TRUE, x should contain the stacked sample statistics
# (as produced by lav_implied_to_vec); this mode is used to compute the
# moment Jacobian for the standard errors
lav_sem_js_stage1_samp <- function(x = NULL, samplestats = FALSE,
                                   eqs = NULL, lavmodel = NULL,
                                   lavpartable = NULL,
                                   lavsamplestats = NULL, lavh1 = NULL,
                                   lavoptions = NULL,
                                   free_directed_idx = NULL,
                                   aggregated = FALSE) {
  # sample moments
  if (samplestats) {
    implied <- lav_vec_to_implied(x, lavmodel = lavmodel)
  } else {
    implied <- lavh1$implied
  }
  x <- lav_model_get_parameters(lavmodel)

  lavpta <- lav_pt_attributes(lavpartable)
  b <- 1L
  s_mat <- implied$cov[[b]]
  m_vec <- if (lavmodel@meanstructure) implied$mean[[b]] else NULL
  n <- lavsamplestats@nobs[[b]]
  ov_names <- lavpta$vnames$ov[[b]]
  lv_names <- lavpta$vnames$lv.regular[[b]]
  markers <- lavpta$vnames$lv.marker[[b]] # named by lv
  nvar <- length(ov_names)
  nfac <- length(lv_names)

  # indicator pattern (nvar x nfac)
  ind_pat <- matrix(0, nrow = nvar, ncol = nfac)
  mm_idx <- which(lavpartable$op == "=~")
  if (length(mm_idx) > 0L && nfac > 0L) {
    ind_pat[cbind(
      match(lavpartable$rhs[mm_idx], ov_names),
      match(lavpartable$lhs[mm_idx], lv_names)
    )] <- 1
  }

  # residual variances of the proxies
  theta <- lav_sem_js_theta(
    s = s_mat, ind_pat = ind_pat, ov_names = ov_names, lv_names = lv_names,
    lavoptions = lavoptions, aggregated = aggregated
  )

  # small-sample (Efron-Morris) shrinkage factor
  csmall <- (n - 3) / (n - 1)
  if (isFALSE(lavoptions$estimator.args[["js_small_sample"]])) {
    csmall <- 1
  }

  # pass 1: scaling (JS); this also provides the loadings needed to scale
  # the aggregated proxies
  out <- lav_sem_js_eqs_pass(
    eqs_b = eqs[[b]], aggregated = FALSE,
    s_mat = s_mat, m_vec = m_vec, n = n,
    ov_names = ov_names, lv_names = lv_names, markers = markers,
    ind_pat = ind_pat, theta = theta, csmall = csmall,
    lambda_mat = NULL, lavpartable = lavpartable, lavmodel = lavmodel,
    x = x
  )
  x <- out$x
  eqs[[b]] <- out$eqs

  # pass 2 (JSA): aggregated proxies, scaled by the pass-1 loadings
  if (aggregated) {
    lambda_mat <- matrix(0, nrow = nvar, ncol = nfac)
    if (length(mm_idx) > 0L) {
      lambda_val <- lavpartable$ustart[mm_idx]
      lambda_val[is.na(lambda_val)] <- 0
      free_mm <- which(lavpartable$free[mm_idx] > 0L)
      lambda_val[free_mm] <- x[lavpartable$free[mm_idx][free_mm]]
      lambda_mat[cbind(
        match(lavpartable$rhs[mm_idx], ov_names),
        match(lavpartable$lhs[mm_idx], lv_names)
      )] <- lambda_val
    }
    out <- lav_sem_js_eqs_pass(
      eqs_b = eqs[[b]], aggregated = TRUE,
      s_mat = s_mat, m_vec = m_vec, n = n,
      ov_names = ov_names, lv_names = lv_names, markers = markers,
      ind_pat = ind_pat, theta = theta, csmall = csmall,
      lambda_mat = lambda_mat, lavpartable = lavpartable,
      lavmodel = lavmodel, x = x
    )
    x <- out$x
    eqs[[b]] <- out$eqs
  }

  # pooled (system) solve for shared free parameters among the directed
  # coefficients (a no-op when nothing is shared)
  x <- lav_sem_miiv_apply_directed_pool(
    eqs_b = do.call(c, eqs), x = x, lavmodel = lavmodel,
    lavpartable = lavpartable
  )

  theta1 <- x[free_directed_idx]
  attr(theta1, "eqs") <- eqs
  theta1
}


# one pass over the equations of a single block: for each equation, build
# the conditioning set Z (scaling indicators or weighted aggregates for the
# latent regressors; the observed regressors themselves), and solve the
# least-squares regression of the dependent variable on the estimated
# conditional expectations -- all in moment form:
#
#   design  d = G S_zz^{-1} Z   with row l of G = Cov(eta_l, Z)
#   amat = Var(d)  = G S_zz^{-1} G'
#   bvec = Cov(d, y) = G S_zz^{-1} Cov(Z, y)
#   slopes = solve(amat, bvec)
#
# for an observed regressor, the corresponding row of G equals the row of
# S_zz, so its 'conditional expectation' is the variable itself
lav_sem_js_eqs_pass <- function(eqs_b = NULL, aggregated = FALSE,
                                s_mat = NULL, m_vec = NULL, n = NULL,
                                ov_names = NULL, lv_names = NULL,
                                markers = NULL, ind_pat = NULL,
                                theta = NULL, csmall = 1,
                                lambda_mat = NULL, lavpartable = NULL,
                                lavmodel = NULL, x = NULL) {
  meanstructure <- lavmodel@meanstructure
  nvar <- length(ov_names)
  single_pat <- rowSums(ind_pat != 0) == 1L

  for (j in seq_along(eqs_b)) {
    eq <- eqs_b[[j]]
    y_idx <- match(eq$lhs_new, ov_names)

    # intercept-only equation
    if (identical(eq$rhs_new, "1")) {
      if (meanstructure && length(eq$ptint) == 1L) {
        free_int_idx <- lavpartable$free[eq$ptint]
        if (free_int_idx > 0L) {
          x[free_int_idx] <- m_vec[y_idx]
        }
      }
      eqs_b[[j]]$coef <- if (meanstructure) m_vec[y_idx] else numeric(0L)
      eqs_b[[j]]$nobs <- n
      next
    }

    nrhs <- length(eq$rhs)
    is_lv <- eq$rhs %in% lv_names

    # conditioning set: one row of u_mat (over the observed variables) per
    # regressor, plus its scale (kappa) and error variance (errvar)
    u_mat <- matrix(0, nrow = nrhs, ncol = nvar)
    kappa <- rep(1, nrhs)
    errvar <- rep(0, nrhs)
    for (i in seq_len(nrhs)) {
      if (!is_lv[i]) {
        u_mat[i, match(eq$rhs_new[i], ov_names)] <- 1
        next
      }
      f <- match(eq$rhs[i], lv_names)
      marker_idx <- match(markers[[eq$rhs[i]]], ov_names)
      agg_idx <- integer(0L)
      if (aggregated) {
        # aggregate over the 'pure' (single-loading) indicators of this
        # factor, excluding the dependent variable itself (leave-one-out
        # for the measurement equations)
        agg_idx <- setdiff(which(ind_pat[, f] != 0 & single_pat), y_idx)
      }
      if (length(agg_idx) > 1L) {
        w <- lav_sem_js_weights(
          s = s_mat[agg_idx, agg_idx, drop = FALSE],
          theta = theta[agg_idx]
        )
        kap <- sum(w * lambda_mat[agg_idx, f])
        if (is.finite(kap) && abs(kap) > .Machine$double.eps^0.5) {
          u_mat[i, agg_idx] <- w
          kappa[i] <- kap
          errvar[i] <- sum(w * w * theta[agg_idx])
          next
        }
      }
      # scaling indicator (also the JSA fallback)
      u_mat[i, marker_idx] <- 1
      errvar[i] <- theta[marker_idx]
    }

    # moments of the conditioning variables
    s_zz <- u_mat %*% s_mat %*% t(u_mat)
    s_zy <- drop(u_mat %*% s_mat[, y_idx])

    # G: covariance of the regressors (latent or observed) with Z
    g_mat <- s_zz
    for (i in which(is_lv)) {
      g_mat[i, i] <- g_mat[i, i] - csmall * errvar[i]
    }
    if (any(is_lv)) {
      g_mat[is_lv, ] <- g_mat[is_lv, , drop = FALSE] / kappa[is_lv]
    }

    # normal equations of the regression on the conditional expectations
    s_zz_inv <- try(solve((s_zz + t(s_zz)) / 2), silent = TRUE)
    if (inherits(s_zz_inv, "try-error")) {
      lav_msg_stop(gettextf(
        "the conditioning variables of the equation for %s are collinear.",
        eq$lhs_new))
    }
    gsi <- g_mat %*% s_zz_inv
    amat <- gsi %*% t(g_mat)
    amat <- (amat + t(amat)) / 2
    bvec <- drop(gsi %*% s_zy)

    # means of the design (for the intercepts)
    if (meanstructure) {
      z_mean <- drop(u_mat %*% m_vec)
      x_bar <- z_mean
      x_bar[is_lv] <- z_mean[is_lv] / kappa[is_lv]
      y_bar <- m_vec[y_idx]
    } else {
      x_bar <- numeric(nrhs)
      y_bar <- 0
    }

    # free/fixed partition and solve
    eq_free_idx <- lavpartable$free[eq$pt]
    b_full <- lavpartable$ustart[eq$pt]
    b_full[is.na(b_full)] <- 0
    fix_idx <- which(eq_free_idx == 0L)
    if (length(fix_idx) > 0L) {
      amat_free <- amat[-fix_idx, -fix_idx, drop = FALSE]
      bvec_free <- drop(bvec[-fix_idx] -
        amat[-fix_idx, fix_idx, drop = FALSE] %*% b_full[fix_idx])
      gcol_free <- eq_free_idx[-fix_idx]
    } else {
      amat_free <- amat
      bvec_free <- bvec
      gcol_free <- eq_free_idx
    }
    if (length(gcol_free) > 0L) {
      slopes_free <- drop(solve(amat_free, bvec_free))
      x[gcol_free] <- slopes_free
      b_full[eq_free_idx > 0L] <- slopes_free
      # slope system for the pooled solve (simple equality constraints)
      eqs_b[[j]]$slope_block <- list(
        amat = (n - 1) * amat_free, bvec = (n - 1) * bvec_free,
        gcol = gcol_free, x_idx = integer(0L),
        x_bar = x_bar, y_bar = y_bar
      )
    }

    # intercept
    beta0 <- y_bar - sum(b_full * x_bar)
    if (meanstructure && length(eq$ptint) == 1L) {
      free_int_idx <- lavpartable$free[eq$ptint]
      if (free_int_idx > 0L) {
        x[free_int_idx] <- beta0
      }
    }

    eqs_b[[j]]$coef <- c(if (meanstructure) beta0 else NULL, b_full)
    eqs_b[[j]]$nobs <- n
  }

  list(x = x, eqs = eqs_b)
}
