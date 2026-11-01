# utility functions needed for lav_cfa_*

# compute THETA and PSI, given lambda using either ULS or GLS
# this function assumes:
# - THETA is diagonal
# - PSI is unrestricted
# - we assume W = S^{-1}
#
# YR 17 oct 2022: - add lower/upper bounds for theta (only to compute PSI)
#                 - use 'lambda' correction to ensure PSI is positive definite
# YR 02 feb 2023: - add psi.mapping.ML argument
lav_cfa_lambda2thetapsi <- function(lambda = NULL, s = NULL, s_inv = NULL,
                                    gls = FALSE, psi_mapping_ml = FALSE,
                                    nobs = 20L) {
  mm_lambda <- as.matrix(lambda)
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)

  if (gls) {
    # see Browne, 1974 section 4 case II
    if (is.null(s_inv)) {
      m_w <- solve(s)
    } else {
      m_w <- s_inv
    }
    t_lw <- crossprod(mm_lambda, m_w)
    m <- solve(t_lw %*% mm_lambda, t_lw) # GLS mapping
    # d <- W %*% mm_lambda %*% M # symmetric
    d <- crossprod(m, t_lw)
    # theta <- solve(W*W - d*d, diag(W %*% S %*% W - d %*% S %*% d))
    theta <- try(solve(m_w * m_w - d * d, diag(m_w - d)), # because W == S^{-1}
      silent = TRUE
    )
    if (inherits(theta, "try-error")) {
      # what to do?
      lav_msg_warn(gettext(
        "problem computing THETA values; trying pace algorithm"))
      theta <- lav_efa_pace(s = s, nfactors = nfac, theta_only = TRUE)
    }
  } else {
    # see Hagglund 1982, section 4
    m <- solve(crossprod(mm_lambda), t(mm_lambda)) # ULS mapping function
    d <- mm_lambda %*% m
    theta <- try(solve(diag(nvar) - d * d, diag(s - (d %*% s %*% d))),
      silent = TRUE
    )
    if (inherits(theta, "try-error")) {
      # what to do?
      lav_msg_warn(gettext(
        "problem computing THETA values; trying pace algorithm"))
      theta <- lav_efa_pace(s = s, nfactors = nfac, theta_only = TRUE)
    }
  }
  theta_nobounds <- theta

  # ALWAYS check bounds for theta (only to compute PSI)!
  theta_bounds <- TRUE
  if (theta_bounds) {
    diag_s <- diag(s)
    # lower bound
    lower_bound <- diag_s * 0 # * 0.01
    too_small_idx <- which(theta < lower_bound)
    if (length(too_small_idx) > 0L) {
      theta[too_small_idx] <- lower_bound[too_small_idx]
    }

    # upper bound
    upper_bound <- diag_s * 1 # * 0.99
    too_large_idx <- which(theta > upper_bound)
    if (length(too_large_idx) > 0L) {
      theta[too_large_idx] <- upper_bound[too_large_idx]
    }
  }

  # psi
  diag_theta <- diag(theta, nvar)
  lambda <- try(lav_mat_sym_diff_smallest_root(s, diag_theta),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    s_min_theta <- s - diag_theta # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda_star <- lambda - 1 / (nobs - 1)
      s_min_theta <- s - lambda_star * diag_theta
    } else {
      s_min_theta <- s - diag_theta
    }
  }

  # just like local SAM
  if (psi_mapping_ml) {
    ti <- 1 / theta
    zero_theta_idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero_theta_idx) > 0L) {
      ti[zero_theta_idx] <- 1
    }
    m <- solve(t(mm_lambda) %*% diag(ti, nvar) %*% mm_lambda) %*%
      t(mm_lambda) %*% diag(ti, nvar)
    mm_psi <- m %*% s_min_theta %*% t(m) # ML
  } else {
    mm_psi <- m %*% s_min_theta %*% t(m) # ULS/GLS
  }

  # we take care of the bounds later!
  list(lambda = mm_lambda, theta = theta_nobounds, psi = mm_psi)
}

# compute PSI, given lambda and theta using either ULS, GLS, ML
# this function assumes:
# - THETA is diagonal
# - PSI is unrestricted
#
# YR 08 Mar 2023: - first version
lav_cfa_lambdatheta2psi <- function(lambda = NULL, theta = NULL, # vector!
                                    s = NULL, s_inv = NULL,
                                    mapping = "ML", nobs = 20L) {
  mm_lambda <- as.matrix(lambda)
  nvar <- nrow(mm_lambda)

  # ALWAYS check bounds for theta to compute PSI
  diag_s <- diag(s)
  # lower bound
  lower_bound <- diag_s * 0 # * 0.01
  too_small_idx <- which(theta < lower_bound)
  if (length(too_small_idx) > 0L) {
    theta[too_small_idx] <- lower_bound[too_small_idx]
  }

  # upper bound
  upper_bound <- diag_s * 1 # * 0.99
  too_large_idx <- which(theta > upper_bound)
  if (length(too_large_idx) > 0L) {
    theta[too_large_idx] <- upper_bound[too_large_idx]
  }

  # psi
  diag_theta <- diag(theta, nvar)
  lambda <- try(lav_mat_sym_diff_smallest_root(s, diag_theta),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    s_min_theta <- s - diag_theta # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda_star <- lambda - 1 / (nobs - 1)
      s_min_theta <- s - lambda_star * diag_theta
    } else {
      s_min_theta <- s - diag_theta
    }
  }

  # mapping matrix
  if (mapping == "ML") {
    ti <- 1 / theta
    zero_theta_idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero_theta_idx) > 0L) {
      ti[zero_theta_idx] <- 1
    }
    m <- solve(t(mm_lambda) %*% diag(ti, nvar) %*% mm_lambda) %*%
         t(mm_lambda) %*% diag(ti, nvar)
  } else if (mapping == "GLS") {
    if (is.null(s_inv)) {
      s_inv <- try(solve(s), silent = TRUE)
    }
    if (inherits(s_inv, "try-error")) {
      m <- tcrossprod(solve(crossprod(mm_lambda)), mm_lambda)
    } else {
      m <- solve(t(mm_lambda) %*% s_inv %*% mm_lambda) %*%
           t(mm_lambda) %*% s_inv
    }
  } else if (mapping == "ULS") {
    m <- tcrossprod(solve(crossprod(mm_lambda)), mm_lambda)
  }

  # compute PSI
  mm_psi <- m %*% s_min_theta %*% t(m)

  mm_psi
}

# compute theta elements for a 1-factor model
#
# Spearman/communality estimator: h2_i is the mean over all tetrads
# r_ia r_ib / r_ab (a, b != i). This is the m = 1 case of the
# Albert/Ihara-Kano partitioned estimator implemented (for m = nfactors,
# with a single greedily chosen partition) by lav_efa_pace() in
# lav_efa_pace.R; see the comments there. The masked version below
# (lav_cfa_theta_spearman_masked, used by estimator = "JS"/"JSA"/"MGM"
# through lav_cfa_theta_pattern) averages only the tetrads that the model
# structure declares admissible, and reproduces this function exactly when
# all tetrads are clean.
lav_cfa_theta_spearman <- function(s, bounds = "wide") {
  p <- ncol(s)
  out <- numeric(p)
  r <- cov2cor(s)
  for (p_idx in seq_len(p)) {
    x <- r[, p_idx][-p_idx]
    aa <- lav_mat_vech(tcrossprod(x), diagonal = FALSE)
    ss <- lav_mat_vech(r[-p_idx, -p_idx, drop = FALSE], diagonal = FALSE)
    h2 <- mean(aa / ss) # communaliteit
    if (bounds == "standard") {
      h2[h2 < 0] <- 0
      h2[h2 > 1] <- 1
    } else if (bounds == "wide") {
      h2[h2 < -0.05] <- -0.05 # corresponds to lower bound ov.var "wide"
      h2[h2 > +1.20] <- +1.20 # corresponds to upper bound ov.var "wide"
    }
    out[p_idx] <- (1 - h2) * s[p_idx, p_idx]
  }
  out
}


# ----------------------------------------------------------------------
# Spearman/communality residual variances with tetrad masking, and the
# per-pattern dispatcher shared by the JS/JSA and MGM estimators
# (moved here from lav_sem_js.R; the model-structure helpers that produce
# the admissibility information -- lav_sem_js_contam, lav_sem_js_graph,
# lav_sem_js_sep -- remain in lav_sem_js.R)
# ----------------------------------------------------------------------
# Spearman/communality residual variances with tetrad masking: like
# lav_cfa_theta_spearman (h2_i = mean over tetrads r_ia * r_ib / r_ab),
# but only the admissible tetrads are averaged: a tetrad (i; a, b) is
# admissible when the three pairs (i, a), (i, b) and (a, b) are all
# admissible, where a pair is admissible when the model implies no
# association between its members other than through the factor
# (sep[a, b]; see lav_sem_js_sep), and when i is a pure indicator of the
# factor (pure[i]). The legacy form (contam) only masks pairs with an
# active residual covariance (sep = !contam), and treats every indicator
# as pure. When a variable has no admissible within-factor tetrad left,
# the tetrads are formed with EXTERNAL variables z instead (Kano-1990
# style: any observed model variable outside the factor; the formula
# r_ia * r_iz / r_az identifies the communality through the factor
# correlations), first pairing an indicator with an external variable
# and then (external_pairs = TRUE) two external variables -- the latter
# also covers the reliability of a single indicator identified through
# the structural model. Returns NA for a variable without any admissible
# tetrad at all. When all within-factor tetrads are admissible, this
# reproduces lav_cfa_theta_spearman exactly.
#
# The tetrad is the m = 1 case of the Albert/Ihara-Kano partitioned
# estimator s_ii - s_{i,C} S_{R,C}^{-1} s_{R,i} that lav_efa_pace()
# (lav_efa_pace.R) uses with m = nfactors and one greedily chosen
# partition; the admissibility conditions here (sep[i, a], sep[i, b],
# sep[a, b]) are the m = 1 form of requiring that no cell of that
# partition carries an association outside the factor. The two are kept
# separate because PACE needs no model structure and may select
# partitions from the data, whereas this function must be a smooth
# function of the sample moments (analytic Jacobian below).
#
# s_full/sep/contam/pure cover ALL observed variables of the block; idx is
# the factor's indicator set; zpool the external candidates.
#
# deriv = TRUE also returns the Jacobian of the theta values with respect
# to vech(s_full) as the "dtheta" attribute (length(idx) x pstar), using
# the moment form of the tetrad: r_ia r_ib / r_ab = S_ia S_ib / (S_ii S_ab)
# (used by the analytic delta-method standard errors; see
# lav_sem_js_deriv.R)
lav_cfa_theta_spearman_masked <- function(s_full = NULL, idx = NULL,
                                             contam = NULL,
                                             sep = NULL, pure = NULL,
                                             zpool = integer(0L),
                                             external_pairs = FALSE,
                                             bounds = "wide",
                                             deriv = FALSE) {
  r <- cov2cor(s_full)
  p <- length(idx)
  nvar <- nrow(s_full)
  if (is.null(sep)) {
    sep <- if (is.null(contam)) {
      matrix(TRUE, nrow = nvar, ncol = nvar)
    } else {
      !contam
    }
  }
  if (is.null(pure)) {
    pure <- rep(TRUE, nvar)
  }
  out <- rep(as.numeric(NA), p)
  dout <- NULL
  cm <- NULL
  if (deriv) {
    pstar <- nvar * (nvar + 1L) / 2L
    dout <- matrix(0, nrow = p, ncol = pstar)
    # vech coordinate of the (i, j) entry (lav_mat_vech order)
    cm <- matrix(0L, nrow = nvar, ncol = nvar)
    cm[lower.tri(cm, diag = TRUE)] <- seq_len(pstar)
    cm <- cm + t(cm) - diag(diag(cm))
  }
  for (k in seq_len(p)) {
    i <- idx[k]
    if (!pure[i]) {
      next # no tetrad identifies the communality of an impure indicator
    }
    others <- setdiff(idx, i)
    # admissible within-factor tetrads; tet collects the (a, b) pairs
    # actually used (for the derivative)
    vals <- numeric(0L)
    tet <- matrix(0L, nrow = 2L, ncol = 0L)
    if (length(others) >= 2L) {
      pairs <- utils::combn(others, 2L)
      for (m in seq_len(ncol(pairs))) {
        a <- pairs[1L, m]
        b <- pairs[2L, m]
        if (!sep[i, a] || !sep[i, b] || !sep[a, b]) {
          next
        }
        vals <- c(vals, r[i, a] * r[i, b] / r[a, b])
        if (deriv) {
          tet <- cbind(tet, c(a, b))
        }
      }
    }
    # fall back to external tetrads when no admissible within-factor
    # tetrad: an indicator paired with an external variable
    if (length(vals) == 0L && length(zpool) > 0L) {
      for (a in others) {
        if (!sep[i, a]) {
          next
        }
        for (z in zpool) {
          if (!sep[i, z] || !sep[a, z]) {
            next
          }
          vals <- c(vals, r[i, a] * r[i, z] / r[a, z])
          if (deriv) {
            tet <- cbind(tet, c(a, z))
          }
        }
      }
    }
    # last resort: two external variables
    if (length(vals) == 0L && external_pairs && length(zpool) >= 2L) {
      zpairs <- utils::combn(zpool, 2L)
      for (m in seq_len(ncol(zpairs))) {
        a <- zpairs[1L, m]
        b <- zpairs[2L, m]
        if (!sep[i, a] || !sep[i, b] || !sep[a, b]) {
          next
        }
        vals <- c(vals, r[i, a] * r[i, b] / r[a, b])
        if (deriv) {
          tet <- cbind(tet, c(a, b))
        }
      }
    }
    if (length(vals) == 0L) {
      next # no admissible tetrad: NA
    }
    h2 <- mean(vals)
    clipped <- FALSE
    if (bounds == "standard") {
      clipped <- h2 < 0 || h2 > 1
      h2[h2 < 0] <- 0
      h2[h2 > 1] <- 1
    } else if (bounds == "wide") {
      clipped <- h2 < -0.05 || h2 > +1.20
      h2[h2 < -0.05] <- -0.05
      h2[h2 > +1.20] <- +1.20
    }
    out[k] <- (1 - h2) * s_full[i, i]
    if (deriv) {
      # theta_k = (1 - h2) * S_ii with h2 the mean over the tetrads
      # v = S_ia S_ib / (S_ii S_ab); when h2 was clipped, dh2 = 0
      if (!clipped) {
        dh2 <- numeric(ncol(dout))
        s_ii <- s_full[i, i]
        for (m in seq_len(ncol(tet))) {
          a <- tet[1L, m]
          b <- tet[2L, m]
          s_ia <- s_full[i, a]
          s_ib <- s_full[i, b]
          s_ab <- s_full[a, b]
          dh2[cm[i, a]] <- dh2[cm[i, a]] + s_ib / (s_ii * s_ab)
          dh2[cm[i, b]] <- dh2[cm[i, b]] + s_ia / (s_ii * s_ab)
          dh2[cm[a, b]] <- dh2[cm[a, b]] - s_ia * s_ib / (s_ii * s_ab^2)
          dh2[cm[i, i]] <- dh2[cm[i, i]] - s_ia * s_ib / (s_ii^2 * s_ab)
        }
        dh2 <- dh2 / ncol(tet)
        dout[k, ] <- -s_full[i, i] * dh2
      }
      dout[k, cm[i, i]] <- dout[k, cm[i, i]] + (1 - h2)
    }
  }
  if (deriv) {
    attr(out, "dtheta") <- dout
  }
  out
}


# residual (measurement-error) variances of the observed indicators of a
# block, given the indicator pattern (ind_pat: nvar x nfac, nonzero for a
# loading); shared by the non-iterative estimators that need these
# variances as an ingredient: estimator = "JS"/"JSA" (shrinkage and, for
# JSA, the aggregation weights; lav_sem_js.R) and estimator = "MGM"
# (initial communalities for the sum-score composites;
# lav_cfa_guttman1952.R). By default the Spearman/communality estimator
# applied per factor (Burghgraeve et al., 2021, Appendix B), averaging only
# the tetrads that the model structure declares admissible (sep_list: per
# factor, the pairs of observed variables that are d-separated by the
# factor; pure_mat: the pure indicators of each factor; see lav_sem_js_sep
# and lav_cfa_theta_spearman_masked). Without sep_list, only the pairs
# with an active residual covariance (contam) are masked. Indicators
# without an admissible tetrad -- including indicators with
# cross-loadings or regressors, whose communality no tetrad identifies --
# get NA; whether that is a problem is up to the caller (JS: only when
# they are needed as a conditioning proxy, see lav_sem_js_eq_plans; MGM:
# always an error).
#
# Options are read from lavoptions$estimator.args, with the same names for
# all estimators: js_theta ("spearman" or "user"), js_theta_bounds
# ("wide", "standard" or "none"), js_theta_values (named vector, for
# js_theta = "user"). theta_fixed (per variable, NA when free) supplies
# the residual variances that are fixed in the model, used for factors
# with fewer than 3 indicators. deriv = TRUE adds the Jacobian with
# respect to vech(s) as the "dtheta" attribute (nvar x pstar).
lav_cfa_theta_pattern <- function(s = NULL, ind_pat = NULL, ov_names = NULL,
                                  lv_names = NULL, lavoptions = NULL,
                                  contam = NULL,
                                  sep_list = NULL, pure_mat = NULL,
                                  theta_fixed = NULL, deriv = FALSE) {
  nvar <- ncol(s)
  pstar <- nvar * (nvar + 1L) / 2L
  ea <- lavoptions$estimator.args
  js_theta <- tolower(ea[["js_theta"]])
  if (length(js_theta) == 0L) {
    js_theta <- "spearman"
  }
  js_theta_bounds <- tolower(ea[["js_theta_bounds"]])
  if (length(js_theta_bounds) == 0L) {
    js_theta_bounds <- "wide"
  }

  # derivative of theta with respect to vech(s) (analytic standard errors)
  dtheta <- if (deriv) matrix(0, nrow = nvar, ncol = pstar) else NULL

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
      # admissible pairs and pure indicators for this factor
      if (!is.null(sep_list)) {
        sep_f <- sep_list[[f]]
      } else if (!is.null(contam)) {
        sep_f <- !contam
      } else {
        sep_f <- matrix(TRUE, nrow = nvar, ncol = nvar)
      }
      pure_f <- if (!is.null(pure_mat)) {
        pure_mat[, f]
      } else {
        rep(TRUE, nvar)
      }
      # every within-factor tetrad admissible?
      all_clean <- length(idx) >= 3L && all(pure_f[idx]) &&
        all(sep_f[idx, idx][upper.tri(sep_f[idx, idx])])
      if (all_clean) {
        if (!deriv) {
          theta[idx] <- lav_cfa_theta_spearman(s[idx, idx, drop = FALSE],
            bounds = js_theta_bounds
          )
        } else {
          # same values as lav_cfa_theta_spearman (no tetrad is masked
          # when the factor block is clean), but with the derivative
          tf <- lav_cfa_theta_spearman_masked(
            s_full = s, idx = idx,
            sep = matrix(TRUE, nrow = nvar, ncol = nvar),
            zpool = integer(0L),
            bounds = js_theta_bounds, deriv = TRUE
          )
          theta[idx] <- as.numeric(tf)
          dtheta[idx, ] <- attr(tf, "dtheta")
        }
      } else if (length(idx) >= 1L) {
        # inadmissible within-factor tetrads are masked; when a variable
        # has none left -- including the 2-indicator case, which has no
        # within-factor tetrads at all, and the single-indicator case --
        # external tetrads take over (Kano-style, with the observed
        # variables outside the factor in the z role; this requires the
        # factor to be correlated with other model variables, and the
        # pairs to be d-separated by the factor)
        tf <- lav_cfa_theta_spearman_masked(
          s_full = s, idx = idx, sep = sep_f, pure = pure_f,
          zpool = setdiff(seq_len(nvar), idx),
          external_pairs = !is.null(sep_list),
          bounds = js_theta_bounds, deriv = deriv
        )
        theta[idx] <- as.numeric(tf)
        if (deriv) {
          dtheta[idx, ] <- attr(tf, "dtheta")
        }
      }
      # a model-fixed residual variance (below) overrides the tetrad
      # value for factors with fewer than 3 indicators; when nothing
      # provides the reliability, theta stays NA and an informative error
      # is raised when the proxy is actually needed
    }

    # residual variances that are FIXED in the model provide the
    # reliability directly for factors where the tetrad estimator is not
    # available (fewer than 3 indicators); e.g. a single-indicator factor
    # with its residual variance fixed to zero, or a 2-indicator factor
    # with one known reliability
    small_fac <- which(colSums(ind_pat != 0) < 3L)
    if (length(small_fac) > 0L && !is.null(theta_fixed)) {
      for (f in small_fac) {
        idx <- which(ind_pat[, f] != 0)
        fix_idx <- idx[!is.na(theta_fixed[idx])]
        theta[fix_idx] <- theta_fixed[fix_idx]
        if (deriv && length(fix_idx) > 0L) {
          dtheta[fix_idx, ] <- 0 # model constants
        }
      }
    }
  }

  # non-indicator variables (observed covariates, endogenous observed
  # variables) are error-free proxies of themselves; indicators without a
  # clean tetrad keep NA
  nonind_idx <- which(rowSums(ind_pat != 0) == 0L)
  theta[nonind_idx] <- 0
  if (deriv && length(nonind_idx) > 0L) {
    dtheta[nonind_idx, ] <- 0
  }

  # keep the variances within [0, var(y)] (NA-safe)
  diag_s <- diag(s)
  too_small_idx <- which(theta < 0)
  theta[too_small_idx] <- 0
  too_large_idx <- which(theta > diag_s)
  if (length(too_large_idx) > 0L) {
    theta[too_large_idx] <- diag_s[too_large_idx]
  }
  if (deriv) {
    if (js_theta == "user") {
      # user-given constants; a value clipped below may still become
      # data-dependent, so zero the derivative BEFORE the clip handling
      dtheta[] <- 0
    }
    if (length(too_small_idx) > 0L) {
      dtheta[too_small_idx, ] <- 0
    }
    if (length(too_large_idx) > 0L) {
      # theta_i = S_ii: derivative is the (i, i) vech coordinate
      cm_diag <- lav_mat_diagh_idx(nvar)
      dtheta[too_large_idx, ] <- 0
      dtheta[cbind(too_large_idx, cm_diag[too_large_idx])] <- 1
    }
    attr(theta, "dtheta") <- dtheta
  }

  theta
}
