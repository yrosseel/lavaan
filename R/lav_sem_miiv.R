# IV/MIIV estimation

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_miiv_internal <- function(lavmodel = NULL, lavh1 = NULL,
                                  lavsamplestats = NULL,
                                  lavpartable = NULL,
                                  lavdata = NULL, lavoptions = NULL) {
  # IV options
  iv_method <- toupper(lavoptions$estimator.args$iv_method)
  stopifnot(iv_method %in% "2SLS")
  iv_varcov_method <- toupper(lavoptions$estimator.args$iv_varcov_method)
  iv_samplestats <- lavoptions$estimator.args$iv_samplestats
  if (lavdata@data.type == "moment") {
    # force iv_samplestats = TRUE
    iv_samplestats <- TRUE
  }
  iv_vcov_stage1 <- tolower(lavoptions$estimator.args$iv_vcov_stage1)
  iv_vcov_stage2 <- tolower(lavoptions$estimator.args$iv_vcov_stage2)
  # [[ ]] exact match: $ would partial-match iv_sargan to iv_sargan_adjust
  iv_sargan <- lavoptions$estimator.args[["iv_sargan"]]
  iv_mean_structure <- tolower(lavoptions$estimator.args$iv_mean_structure)
  if (length(iv_mean_structure) == 0L) {
    iv_mean_structure <- "moments"
  }
  # ML divisor (N) instead of the unbiased N-P divisor for the residual
  # variances (see lav_options_est_iv())
  iv_mimic_ml <- isTRUE(lavoptions$estimator.args$iv_mimic_ml)
  # just in case
  if (lavmodel@categorical) {
    iv_samplestats <- TRUE
  }
  stopifnot(iv_varcov_method %in% c("ULS", "GLS", "2RLS", "RLS", "NONE"))

  # get lavpta
  # lavpta <- lav_pt_attributes(lavpartable)

  # we assume the blocks are independent groups for now
  stopifnot(lavdata@nlevels == 1L)

  # directed versus undirected (free) parameters
  undirected_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  directed_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~", "~1"))
  free_directed_idx <- unique(lavpartable$free[directed_idx])
  free_undirected_idx <- unique(lavpartable$free[undirected_idx])

  # external instruments with categorical data: equality constraints among the
  # directed coefficients (eg measurement invariance) are not supported, because
  # the constraint-aware standard errors would require a numerical moment
  # Jacobian, which is not available for categorical data (the analytic
  # per-equation Jacobian cannot represent the pooling). Fail cleanly rather
  # than report unreliable standard errors.
  if (lavmodel@categorical && length(unlist(lavdata@ov.names.aux)) > 0L) {
    dir_slope_free <- lavpartable$free[lavpartable$free > 0L &
      lavpartable$op %in% c("=~", "~")]
    if (anyDuplicated(dir_slope_free) > 0L ||
        length(lavmodel@ceq.linear.idx) > 0L) {
      lav_msg_stop(gettext(
        "external instruments with categorical data do not support equality
         constraints among the directed coefficients (eg measurement
         invariance) yet."))
    }
  }

  # general linear equality constraints are handled stage-wise: among the
  # directed slopes (stage 1) and among the undirected parameters (stage 2).
  # Warn once about constraints we cannot honor this way (constraints on
  # intercepts, or spanning both stages).
  if (length(lavmodel@ceq.linear.idx) > 0L) {
    noint_idx <- which(lavpartable$free > 0L & !duplicated(lavpartable$free) &
      lavpartable$op %in% c("=~", "~"))
    free_noint_idx <- unique(lavpartable$free[noint_idx])
    con_dir <- lav_sem_miiv_linear_con(lavmodel, free_noint_idx)
    con_und <- lav_sem_miiv_linear_con(lavmodel, free_undirected_idx)
    n_handled <- (if (is.null(con_dir)) 0L else nrow(con_dir$jac)) +
      (if (is.null(con_und)) 0L else nrow(con_und$jac))
    if (length(lavmodel@ceq.linear.idx) > n_handled) {
      lav_msg_warn(gettext(
        "[IV] some linear equality constraints involve intercepts or span both
         the directed and undirected parameters; these are not honored by the
         (sequential) IV estimator and are ignored."))
    }
  }

  # find ALL model-implied instrumental variables (miivs) per equation
  eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)

  # weak-instrument diagnostics (warn) or pruning, based on the first-stage
  # F-statistic of each overidentified equation
  # use [[ ]] (exact match): $ would partial-match iv_weak to iv_weak_threshold
  iv_weak <- tolower(lavoptions$estimator.args[["iv_weak"]])
  if (length(iv_weak) > 0L && iv_weak != "none") {
    eqs <- lav_sem_miiv_weak(
      eqs = eqs, lavpta = lav_pt_attributes(lavpartable),
      lavh1 = lavh1, lavsamplestats = lavsamplestats, lavdata = lavdata,
      threshold = lavoptions$estimator.args[["iv_weak_threshold"]],
      action = iv_weak
    )
  }

  # initial parameter vector
  x <- lav_model_get_parameters(lavmodel)

  ########################################
  # first stage: directed parameter only #
  ########################################
  theta1 <- numeric(0L)
  if (length(free_directed_idx) > 0L) {
    if (iv_samplestats) {
      theta1 <- lav_sem_miiv_2sls_samp(
        x = NULL, samplestats = FALSE, eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavdata = lavdata, lavsamplestats = lavsamplestats,
        lavh1 = lavh1, free_directed_idx = free_directed_idx,
        iv_vcov_stage1 = iv_vcov_stage1, iv_sargan = iv_sargan,
        iv_vcov_stage2 = iv_vcov_stage2, iv_mimic_ml = iv_mimic_ml
      )
    } else {
      theta1 <- lav_sem_miiv_2sls(
        eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavdata = lavdata, free_directed_idx = free_directed_idx,
        iv_vcov_stage1 = iv_vcov_stage1, iv_sargan = iv_sargan,
        iv_vcov_stage2 = iv_vcov_stage2, iv_mimic_ml = iv_mimic_ml
      )
    }
    # update equations
    eqs <- attr(theta1, "eqs")
    theta1 <- as.numeric(theta1) # drop attributes
    # store theta1 elements in x
    x[free_directed_idx] <- theta1
  }

  #######################################
  # second stage: undirected parameters #
  #######################################

  # compute theta2 using ULS/GLS/RLS/2RLS
  if (length(free_undirected_idx) > 0L && iv_varcov_method != "NONE") {
    theta2 <- lav_sem_miiv_varcov(
      x = theta1,
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats,
      lavh1 = lavh1, free_directed_idx = free_directed_idx,
      free_undirected_idx = free_undirected_idx,
      iv_varcov_method = iv_varcov_method
    )
  } else {
    theta2 <- rep(as.numeric(NA), length(free_undirected_idx))
  }
  # store theta2 elements in x
  x[free_undirected_idx] <- theta2

  # mean structure: re-estimate the free mean parameters (observed intercepts
  # and latent means) jointly by GLS (iv_mean_structure = "wls"), which matches
  # the ML mean estimate and pools shared intercepts/latent means across groups
  # (e.g. for scalar invariance). The default "moments" path keeps the
  # per-equation, nobs-pooled estimates from stage 1.
  if (lavmodel@meanstructure && identical(iv_mean_structure, "wls")) {
    free_mean_idx <- unique(lavpartable$free[lavpartable$op == "~1" &
      lavpartable$free > 0L & !duplicated(lavpartable$free)])
    if (length(free_mean_idx) > 0L) {
      x[free_mean_idx] <- lav_sem_miiv_mean_wls(
        lavmodel = lavmodel, lavsamplestats = lavsamplestats, x = x,
        free_mean_idx = free_mean_idx)
    }
  }

  attr(x, "eqs") <- eqs
  x
}

# joint mean-structure GLS solve. Given the estimated loadings/regressions and
# variances, re-estimate all free mean parameters (observed intercepts nu and
# latent means alpha) by a single GLS solve that fits the model-implied means
# to the sample means across all groups, weighting by nobs * Sigma^{-1}. The
# model-implied means are linear in the mean parameters, so this is a linear
# solve; shared mean parameters (e.g. equal intercepts for scalar invariance)
# are pooled jointly, reproducing the ML mean estimate (given Sigma). The "vcov"
# attribute holds the GLS information inverse A^{-1} (= the SE covariance for
# the continuous, complete-data case).
lav_sem_miiv_mean_wls <- function(lavmodel, lavsamplestats, x, free_mean_idx,
                                  sample_mean = NULL) {
  implied <- lav_model_implied(lav_model_set_parameters(lavmodel, x = x))
  delta_list <- lav_sem_miiv_delta(lav_model_set_parameters(lavmodel, x = x))
  # fixed mean contribution (free mean parameters set to 0)
  x0 <- x
  x0[free_mean_idx] <- 0
  implied0 <- lav_model_implied(lav_model_set_parameters(lavmodel, x = x0))

  nfree <- length(free_mean_idx)
  amat <- matrix(0, nfree, nfree)
  bvec <- numeric(nfree)
  for (g in seq_len(lavmodel@nblocks)) {
    nvar <- lavmodel@nvar[[g]]
    m_g <- delta_list[[g]][seq_len(nvar), free_mean_idx, drop = FALSE]
    mean_g <- if (is.null(sample_mean)) {
      lavsamplestats@mean[[g]]
    } else {
      sample_mean[[g]]
    }
    resid <- mean_g - implied0$mean[[g]]
    w_g <- try(lav_mat_sym_inverse(implied$cov[[g]]), silent = TRUE)
    if (inherits(w_g, "try-error")) {
      w_g <- diag(nvar)
    }
    w_g <- w_g * lavsamplestats@nobs[[g]]
    amat <- amat + crossprod(m_g, w_g %*% m_g)
    bvec <- bvec + crossprod(m_g, w_g %*% resid)
  }
  ainv <- lav_mat_sym_inverse(amat)
  out <- drop(ainv %*% bvec)
  attr(out, "vcov") <- ainv
  out
}


# external instruments (estimator = "IV"): build the saturated moments over
# [model variables, external instruments] for block 'b' by augmenting the
# model (saturated) covariance/mean with the extra observed variables
# (lavData@aux). The model-variable block is kept exactly equal to the supplied
# model moments; the instrument blocks (and means) come from the raw data.
# Model variables come first, instruments last (so the model sub-block can be
# recovered later). Returns NULL when there are no external instruments for
# this block (so callers keep their model-only moments unchanged).
lav_sem_miiv_aug_moments <- function(lavdata = NULL, b = 1L, ov_names = NULL,
                                     sample_cov = NULL, sample_mean = NULL) {
  aux_names_b <- if (length(lavdata@ov.names.aux) >= b) {
    lavdata@ov.names.aux[[b]]
  } else {
    character(0L)
  }
  if (length(aux_names_b) == 0L || length(lavdata@X) < b ||
    is.null(lavdata@X[[b]]) || is.null(lavdata@aux[[b]])) {
    return(NULL)
  }
  p_model <- length(ov_names)
  full_b <- cbind(lavdata@X[[b]], lavdata@aux[[b]])
  aug_cov <- stats::cov(full_b, use = "pairwise.complete.obs")
  aug_cov[seq_len(p_model), seq_len(p_model)] <- sample_cov
  aug_mean <- NULL
  if (!is.null(sample_mean)) {
    aug_mean <- colMeans(full_b, na.rm = TRUE)
    aug_mean[seq_len(p_model)] <- sample_mean
  }
  list(cov = aug_cov, mean = aug_mean, ov.names = c(ov_names, aux_names_b))
}


# external instruments, CATEGORICAL data: build the augmented sample statistics
# over [model variables, instruments] for block 'b'. The instruments must be
# ordinal (the IV/PIV estimator works on the all-ordinal threshold/correlation
# structure). Runs muthen1984() over the augmented data to obtain the augmented
# polychoric correlation matrix R, the thresholds, and the asymptotic
# covariance of the sample statistics (NACOV = WLS.W * nobs). The model
# correlation block (and the model thresholds) are overwritten with the model's
# own saturated values, so that the model-variable estimates are unchanged.
# Returns NULL when there are no external instruments for this block.
lav_sem_miiv_aug_cat_moments <- function(lavdata = NULL, b = 1L,
                                         model_cor = NULL, model_th = NULL,
                                         nobs = NULL) {
  aux_names_b <- if (length(lavdata@ov.names.aux) >= b) {
    lavdata@ov.names.aux[[b]]
  } else {
    character(0L)
  }
  if (length(aux_names_b) == 0L || length(lavdata@X) < b ||
    is.null(lavdata@X[[b]]) || is.null(lavdata@aux[[b]])) {
    return(NULL)
  }
  model_ov <- lavdata@ov.names[[b]]
  p_model <- length(model_ov)

  # variable types/levels for the augmented set
  m_type <- lavdata@ov$type[match(model_ov, lavdata@ov$name)]
  m_lev <- lavdata@ov$nlev[match(model_ov, lavdata@ov$name)]
  a_type <- lavdata@ov.aux$type[match(aux_names_b, lavdata@ov.aux$name)]
  a_lev <- lavdata@ov.aux$nlev[match(aux_names_b, lavdata@ov.aux$name)]

  # ordinal instruments only (for now)
  if (any(a_type != "ordered")) {
    lav_msg_stop(gettextf(
      "external instrument(s) must be ordered (categorical) when the model is
       categorical; continuous instrument(s): %s",
      paste(aux_names_b[a_type != "ordered"], collapse = " ")))
  }

  full_b <- cbind(lavdata@X[[b]], lavdata@aux[[b]])
  cat1 <- muthen1984(
    data_1 = full_b, ov_names = c(model_ov, aux_names_b),
    ov_types = c(m_type, a_type), ov_levels = c(m_lev, a_lev),
    wls_w = TRUE, allow_empty_cell = FALSE
  )

  R <- unname(cat1$COR)
  th <- unlist(cat1$TH)
  th_idx <- unlist(cat1$TH.IDX)
  # keep the model block exactly equal to the model's saturated statistics
  if (!is.null(model_cor)) {
    R[seq_len(p_model), seq_len(p_model)] <- model_cor
  }
  if (!is.null(model_th)) {
    th[th_idx <= p_model] <- model_th
  }
  nacov <- if (!is.null(cat1$WLS.W) && !is.null(nobs)) {
    cat1$WLS.W * nobs
  } else {
    NULL
  }

  list(R = R, th = th, th.idx = th_idx, nacov = nacov,
       ov.names = c(model_ov, aux_names_b), p.model = p_model)
}


# external instruments, CATEGORICAL data: the column indices, within the
# augmented categorical moment vector [thresholds, off-diagonal correlations],
# that correspond to the model moments (model thresholds + correlations among
# model variables). 'th.idx' is the per-threshold variable index of the
# augmented threshold block.
lav_sem_miiv_aug_cat_keep_idx <- function(p_model = 0L, p_aug = 0L,
                                          th_idx = integer(0L)) {
  nth <- length(th_idx)
  th_keep <- which(th_idx <= p_model)
  ri <- lav_mat_vech_row_idx(p_aug)
  ci <- lav_mat_vech_col_idx(p_aug)
  offdiag <- ri > ci
  ri_off <- ri[offdiag]
  ci_off <- ci[offdiag]
  cor_keep <- which(ri_off <= p_model & ci_off <= p_model)
  c(th_keep, nth + cor_keep)
}


# external instruments: assemble the full (augmented-moment) directed Jacobian
# jac_k for block 'b' from the per-equation k_mat_full matrices stored by
# lav_sem_miiv_2sls_samp(). This mirrors lav_sem_miiv_utils_jack_eqs(), but the
# columns span the AUGMENTED moments [model variables, instruments] (so the
# directed coefficients' dependence on the instrument moments is retained).
# Returns NULL when no augmented Jacobian is available for this block.
lav_sem_miiv_jack_eqs_aug <- function(eqs = NULL, block = 1L, lavmodel = NULL,
                                      lavpartable = NULL,
                                      free_directed_idx = integer(0L)) {
  b <- block
  ntheta1 <- length(free_directed_idx)
  # augmented moment dimension, taken from the first available k_mat_full
  m_aug <- NULL
  for (j in seq_along(eqs[[b]])) {
    if (!is.null(eqs[[b]][[j]]$k_mat_full)) {
      m_aug <- ncol(eqs[[b]][[j]]$k_mat_full)
      break
    }
  }
  if (is.null(m_aug)) {
    return(NULL)
  }
  k_mat <- matrix(0.0, nrow = ntheta1, ncol = m_aug)
  for (j in seq_along(eqs[[b]])) {
    eq <- eqs[[b]][[j]]
    if (is.null(eq$k_mat_full)) {
      next
    }
    tmp <- lavpartable$free[eq$pt]
    tmp <- tmp[tmp != 0]
    free_idx <- match(tmp, free_directed_idx)
    nx <- nrow(eq$k_mat_full)
    if (length(free_idx) > 0L && nx > 0L) {
      if (all(free_idx > 0L)) {
        k_mat[free_idx, ] <- eq$k_mat_full
      } else {
        zero_idx <- which(free_idx == 0L)
        free_idx <- free_idx[-zero_idx]
        if (length(free_idx) > 0L) {
          k_mat[free_idx, ] <- eq$k_mat_full[-zero_idx, , drop = FALSE]
        }
      }
    }
    if (lavmodel@meanstructure) {
      tmp <- lavpartable$free[eq$ptint]
      tmp <- tmp[tmp != 0]
      free_int_idx <- match(tmp, free_directed_idx)
      if (length(free_int_idx) > 0L && free_int_idx > 0L) {
        k_mat[free_int_idx, ] <- eq$k_mat_int_full
      }
    }
  }
  k_mat
}


# external instruments: map a model-moment quantity into the augmented moment
# layout [model variables, instruments]. Returns the column indices, within the
# augmented (mean, vech(cov)) vector, that correspond to the model moments.
lav_sem_miiv_aug_keep_idx <- function(p_model = 0L, p_aug = 0L,
                                      meanstructure = FALSE) {
  if (meanstructure) {
    lav_aux_moment_idx(p_model, p_aug)
  } else {
    ri <- lav_mat_vech_row_idx(p_aug)
    ci <- lav_mat_vech_col_idx(p_aug)
    which(ri <= p_model & ci <= p_model)
  }
}


# external instruments: build the augmented saturated moment vector at the
# estimate, concatenated over blocks (each block: [mean, vech(cov)] of the
# augmented [model variables, instruments] covariance). Returns the vector and
# the per-block augmented variable count (p_aug), as consumed by the aug_vec
# mode of lav_sem_miiv_2sls_samp().
lav_sem_miiv_aug_vec0 <- function(lavdata = NULL, lavmodel = NULL,
                                  lavh1 = NULL) {
  ms <- lavmodel@meanstructure
  nblocks <- lavmodel@nblocks
  vec <- numeric(0L)
  pdim <- integer(nblocks)
  for (b in seq_len(nblocks)) {
    augb <- lav_sem_miiv_aug_moments(lavdata = lavdata, b = b,
      ov_names = lavdata@ov.names[[b]],
      sample_cov = lavh1$implied$cov[[b]],
      sample_mean = if (ms) lavh1$implied$mean[[b]] else NULL)
    pdim[b] <- ncol(augb$cov)
    seg <- if (ms) {
      c(augb$mean, lav_matrix_vech(augb$cov))
    } else {
      lav_matrix_vech(augb$cov)
    }
    vec <- c(vec, seg)
  }
  list(vec = vec, pdim = pdim)
}


# external instruments: the full (augmented-moment) directed Jacobian computed
# NUMERICALLY, by differentiating the (constraint-aware) 2SLS directed solve
# with respect to the augmented moment vector. This is used when there are
# equality constraints among the directed coefficients (which the per-equation
# analytic Jacobian, lav_sem_miiv_jack_eqs_aug(), cannot represent because it
# ignores the cross-equation pooling). Returns an n_directed x ntot_aug matrix,
# with the columns in the same per-block [mean, vech(cov)] order used elsewhere.
lav_sem_miiv_jac_full_numeric <- function(eqs = NULL, lavmodel = NULL,
                                          lavpartable = NULL, lavdata = NULL,
                                          lavsamplestats = NULL, lavh1 = NULL,
                                          free_directed_idx = NULL) {
  av <- lav_sem_miiv_aug_vec0(lavdata = lavdata, lavmodel = lavmodel,
                              lavh1 = lavh1)
  numDeriv::jacobian(func = function(v) {
    lav_sem_miiv_2sls_samp(x = NULL, samplestats = TRUE, eqs = eqs,
      lavmodel = lavmodel, lavpartable = lavpartable, lavdata = lavdata,
      lavsamplestats = lavsamplestats, lavh1 = lavh1,
      free_directed_idx = free_directed_idx,
      aug_vec = v, aug_pdim = av$pdim)
  }, x = av$vec)
}


# assess instrument strength for each IV equation via the first-stage
# F-statistic (Staiger & Stock, 1997): for each endogenous regressor in an
# equation, regress it on the equation's instruments and compute the F for the
# instruments. A minimum F below 'threshold' (rule of thumb: 10) flags weak
# instruments. 'action' controls the response:
#  - "warn":  warn about the affected equations (default); the user can supply
#             stronger instruments manually via the |~ operator
#  - "prune": greedily drop the weakest excess instrument from weak,
#             overidentified equations (never below just-identified)
#  - "none":  do nothing (caller skips this function)
# Returns the (possibly pruned) eqs list.
lav_sem_miiv_weak <- function(eqs = NULL, lavpta = NULL, lavh1 = NULL,
                              lavsamplestats = NULL, lavdata = NULL,
                              threshold = 10,
                              action = "warn") {
  # first-stage R^2 of regressor 'xi' on instruments 'zi' (indices)
  r2_fun <- function(s, xi, zi) {
    if (length(zi) == 0L) {
      return(0)
    }
    beta <- try(solve(s[zi, zi, drop = FALSE], s[zi, xi]), silent = TRUE)
    if (inherits(beta, "try-error")) {
      return(NA_real_)
    }
    min(max(drop(s[xi, zi, drop = FALSE] %*% beta) / s[xi, xi], 0), 1 - 1e-12)
  }
  # first-stage F from an R^2 with nz instruments and N observations
  f_fun <- function(r2, nz, n) {
    df_res <- n - nz - 1L
    if (is.na(r2) || df_res <= 0L || nz <= 0L) {
      return(NA_real_)
    }
    (r2 / nz) / ((1 - r2) / df_res)
  }

  weak_msgs <- character(0L)
  pruned_msgs <- character(0L)

  for (b in seq_along(eqs)) {
    ov_names <- lavpta$vnames$ov[[b]]
    s <- lavh1$implied$cov[[b]]
    n <- lavsamplestats@nobs[[b]]
    # external instruments: use the covariance augmented with the extra
    # observed variables (lavData@aux), so instruments not in the model can be
    # assessed for strength
    if (!is.null(lavdata)) {
      aug <- lav_sem_miiv_aug_moments(
        lavdata = lavdata, b = b, ov_names = ov_names, sample_cov = s
      )
      if (!is.null(aug)) {
        s <- aug$cov
        ov_names <- aug$ov.names
      }
    }
    for (j in seq_along(eqs[[b]])) {
      eq <- eqs[[b]][[j]]
      iv_flag <- if (!is.null(eq$iv_type)) {
        eq$iv_type != "ols"
      } else {
        !identical(eq$rhs_new, eq$miiv)
      }
      if (!iv_flag || identical(eq$rhs_new, "1")) {
        next
      }
      endog <- setdiff(eq$rhs_new, eq$iv) # regressors needing instruments
      if (length(endog) == 0L) {
        next
      }
      iv <- eq$iv
      min_f <- function(iv_set) {
        zi <- match(iv_set, ov_names)
        min(vapply(endog, function(xv) {
          f_fun(r2_fun(s, match(xv, ov_names), zi), length(iv_set), n)
        }, numeric(1)), na.rm = TRUE)
      }
      cur_f <- min_f(iv)

      # prune the weakest excess instruments
      if (action == "prune" && is.finite(cur_f) && cur_f < threshold) {
        repeat {
          if (length(iv) <= length(eq$rhs_new) || min_f(iv) >= threshold) {
            break
          }
          # binding (weakest) regressor; drop the excluded instrument whose
          # removal best preserves that regressor's R^2
          f_each <- vapply(endog, function(xv) {
            f_fun(r2_fun(s, match(xv, ov_names), match(iv, ov_names)),
                  length(iv), n)
          }, numeric(1))
          xi <- match(endog[which.min(f_each)], ov_names)
          cand <- setdiff(iv, eq$rhs_new)
          if (length(cand) == 0L) {
            break
          }
          keep_r2 <- vapply(cand, function(z) {
            r2_fun(s, xi, match(setdiff(iv, z), ov_names))
          }, numeric(1))
          iv <- setdiff(iv, cand[which.max(keep_r2)])
        }
        if (!identical(iv, eq$iv)) {
          new_f <- min_f(iv)
          pruned_msgs <- c(pruned_msgs, sprintf(
            "%s ~ %s: dropped %s (first-stage F %.1f -> %.1f)",
            eq$lhs_new, paste(eq$rhs_new, collapse = " + "),
            paste(setdiff(eq$iv, iv), collapse = ", "), cur_f, new_f))
          eqs[[b]][[j]]$iv <- iv
          eqs[[b]][[j]]$iv_type <- "user" # fixed set downstream
          cur_f <- new_f
        }
      }

      if (is.finite(cur_f) && cur_f < threshold) {
        weak_msgs <- c(weak_msgs, sprintf(
          "%s ~ %s (first-stage F = %.1f)",
          eq$lhs_new, paste(eq$rhs_new, collapse = " + "), cur_f))
      }
    }
  }

  if (length(pruned_msgs) > 0L) {
    lav_msg_warn(gettextf(
      "[IV] weak instruments pruned (iv_weak = \"prune\"):\n%s",
      paste(" -", pruned_msgs, collapse = "\n")))
  }
  if (length(weak_msgs) > 0L) {
    lav_msg_warn(gettextf(
      "[IV] weak instruments (first-stage F < %g) for:\n%s\nConsider supplying
       stronger instruments via the |~ operator, or set
       estimator.args = list(..., iv_weak = \"prune\") to drop weak ones.",
      threshold, paste(" -", weak_msgs, collapse = "\n")))
  }

  eqs
}


# pool the directed (slope) coefficients across equations subject to simple
# equality constraints (ie parameters that share the same 'free' column, such
# as equal factor loadings in different equations).
#
# each equation j contributes a quadratic block to a stacked least-squares
# problem: amat_j = Xhat_j' Xhat_j (centered) and bvec_j = Xhat_j' y_j. Without
# constraints, the stacked system is block-diagonal and the joint solution
# reproduces the per-equation 2SLS estimates exactly. Simple equality
# constraints become linear equality restrictions across the blocks, which we
# enforce with quadprog::solve.QP() (meq = number of equality rows).
#
# 'blocks' is a list with, per equation, $amat (nx x nx), $bvec (nx) and
# $gcol (length nx: the global 'free' column for each slope). The per-equation
# blocks are accumulated into a single system D (over the distinct free slope
# columns) and d; simple equality constraints (parameters that share a free
# column) are handled automatically by this accumulation. General linear
# equality constraints are passed in 'con_jac' (n_con x n, columns aligned to
# sort(unique(gcol))) and 'con_rhs' (length n_con) and enforced with
# quadprog::solve.QP(). Returns a named numeric vector: the pooled estimate per
# distinct free slope column.
lav_sem_miiv_pool_directed <- function(blocks = NULL, con_jac = NULL,
                                       con_rhs = NULL) {
  # drop equations without slopes
  blocks <- blocks[vapply(blocks, function(b) length(b$gcol) > 0L, logical(1))]
  if (length(blocks) == 0L) {
    return(numeric(0L))
  }

  # distinct free slope columns
  free_slope_idx <- sort(unique(unlist(lapply(blocks, `[[`, "gcol"))))
  n <- length(free_slope_idx)

  # accumulate the per-equation blocks into one system (shared columns sum)
  dmat <- matrix(0, n, n)
  dvec <- numeric(n)
  for (bl in blocks) {
    p <- match(bl$gcol, free_slope_idx)
    dmat[p, p] <- dmat[p, p] + bl$amat
    dvec[p] <- dvec[p] + bl$bvec
  }

  if (!is.null(con_jac) && nrow(con_jac) > 0L) {
    sol <- quadprog::solve.QP(
      Dmat = dmat, dvec = dvec, Amat = t(con_jac),
      bvec = con_rhs, meq = nrow(con_jac)
    )$solution
  } else {
    sol <- solve(dmat, dvec)
  }

  names(sol) <- as.character(free_slope_idx)
  sol
}


# extract the linear equality constraints (rows of ceq.JAC) that involve ONLY
# the free columns in 'target_idx'. Returns NULL if there are none, else a list
# with $jac (n_con x length(target_idx), columns aligned to 'target_idx' in the
# given order) and $rhs (length n_con). Constraints whose support is not fully
# inside 'target_idx' (eg cross-stage or intercept constraints) are skipped.
lav_sem_miiv_linear_con <- function(lavmodel = NULL, target_idx = NULL) {
  lin_idx <- lavmodel@ceq.linear.idx
  if (length(lin_idx) == 0L || is.null(lavmodel@ceq.JAC) ||
      nrow(lavmodel@ceq.JAC) == 0L || length(target_idx) == 0L) {
    return(NULL)
  }
  jac <- lavmodel@ceq.JAC[lin_idx, , drop = FALSE]
  rhs <- lavmodel@ceq.rhs[lin_idx]
  keep <- apply(jac, 1L, function(r) {
    nz <- which(abs(r) > 0)
    length(nz) > 0L && all(nz %in% target_idx)
  })
  if (!any(keep)) {
    return(NULL)
  }
  list(jac = jac[keep, target_idx, drop = FALSE], rhs = rhs[keep])
}


# apply the pooled directed solve to the parameter vector 'x'. This is engaged
# when a free (directed slope) column is shared by more than one equation (a
# simple equality constraint, handled via column accumulation) OR when there is
# a general linear equality constraint among the directed slopes (handled via
# quadprog::solve.QP). When neither applies, the per-equation fills are already
# correct and 'x' is returned unchanged. Intercepts (if free) are recomputed
# from the pooled slopes.
lav_sem_miiv_apply_directed_pool <- function(eqs_b = NULL, x = NULL,
                                             lavmodel = NULL,
                                             lavpartable = NULL) {
  slope_blocks <- lapply(eqs_b, `[[`, "slope_block")
  slope_blocks <- slope_blocks[!vapply(slope_blocks, is.null, logical(1))]
  if (length(slope_blocks) == 0L) {
    return(x)
  }
  all_gcol <- unlist(lapply(slope_blocks, `[[`, "gcol"))
  free_slope_idx <- sort(unique(all_gcol))

  # general linear equality constraints among the directed slopes (if any)
  con <- lav_sem_miiv_linear_con(lavmodel, free_slope_idx)

  # intercepts that are shared across equations/groups (e.g. equal intercepts
  # for scalar measurement invariance) must also be pooled, not left as
  # last-write-wins
  all_int <- unlist(lapply(eqs_b, function(e) {
    if (!is.null(e$slope_block) && !is.null(e$ptint) &&
        length(e$ptint) == 1L) {
      lavpartable$free[e$ptint]
    } else {
      integer(0L)
    }
  }))
  shared_int <- anyDuplicated(all_int[all_int > 0L]) > 0L

  # nothing to pool?
  if (anyDuplicated(all_gcol) == 0L && is.null(con) && !shared_int) {
    return(x)
  }

  theta_pool <- lav_sem_miiv_pool_directed(
    slope_blocks,
    con_jac = if (is.null(con)) NULL else con$jac,
    con_rhs = if (is.null(con)) NULL else con$rhs
  )
  # 1. write the pooled slopes
  for (j in seq_along(eqs_b)) {
    sb <- eqs_b[[j]]$slope_block
    if (is.null(sb)) {
      next
    }
    x[sb$gcol] <- theta_pool[as.character(sb$gcol)]
  }
  # 2. recompute the intercepts from the pooled slopes. When an intercept free
  #    column is shared across equations/groups, pool the per-equation
  #    intercepts by an nobs-weighted average (consistent with the nobs-scaled
  #    cross-products used for the slopes); a non-shared intercept reduces to
  #    its single per-equation value.
  if (lavmodel@meanstructure) {
    int_idx <- integer(0L)
    int_val <- numeric(0L)
    int_wgt <- numeric(0L)
    for (j in seq_along(eqs_b)) {
      sb <- eqs_b[[j]]$slope_block
      if (is.null(sb)) {
        next
      }
      eq <- eqs_b[[j]]
      free_int_idx <- lavpartable$free[eq$ptint]
      if (length(free_int_idx) != 1L || free_int_idx <= 0L) {
        next
      }
      # full slope vector in rhs order (free pooled values + fixed ustart)
      eq_free_idx <- lavpartable$free[eq$pt]
      b_full <- lavpartable$ustart[eq$pt]
      b_full[is.na(b_full)] <- 0
      free_pos <- which(eq_free_idx > 0L)
      b_full[free_pos] <- x[eq_free_idx[free_pos]]
      beta0 <- as.vector(sb$y_bar - sb$x_bar %*% b_full)
      int_idx <- c(int_idx, free_int_idx)
      int_val <- c(int_val, beta0)
      int_wgt <- c(int_wgt, if (!is.null(eq$nobs)) eq$nobs else 1)
    }
    if (length(int_idx) > 0L) {
      num <- tapply(int_wgt * int_val, int_idx, sum)
      den <- tapply(int_wgt, int_idx, sum)
      pooled_int <- num / den
      x[as.integer(names(pooled_int))] <- pooled_int
    }
  }
  x
}


# stage 1: use 2SLS to find the regression coefficients of all
#          directed effects in the model -- continuous/raw-data version
lav_sem_miiv_2sls <- function(eqs = NULL, lavmodel = NULL, lavpartable = NULL,
                              lavdata = NULL, free_directed_idx = NULL,
                              iv_vcov_stage1 = "lm.vcov.dfres",
                              iv_vcov_stage2 = "delta",
                              iv_sargan = TRUE, iv_mimic_ml = FALSE) {
  # this function is for continuous/raw-data only
  stopifnot(lavdata@data.type == "full")
  stopifnot(!lavmodel@categorical)
  iv_vcov_stage1 <- tolower(iv_vcov_stage1)
  stopifnot(iv_vcov_stage1 %in% c("lm.vcov.dfres", "lm.vcov", "none"))

  # get lavpta
  lavpta <- lav_pt_attributes(lavpartable)

  if (is.null(eqs)) {
    # find ALL model-implied instrumental variables (miivs) per equation
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
  }

  # number of blocks
  nblocks <- lavmodel@nblocks

  # parameter vector
  x <- lav_model_get_parameters(lavmodel)

  # 2SLS per equation
  # for now: - no equality constraints (yet)
  #          - only OLS (not robust, no lasso, ...)
  for (b in seq_len(nblocks)) {
    # ov.names for this block
    ov_names <- lavpta$vnames$ov[[b]]

    # raw data for this block
    xy <- lavdata@X[[b]]

    # external instruments (estimator = "IV"): append the extra observed
    # variables (lavData@aux) as additional columns so that instruments which
    # are not part of the model can be used by the per-equation 2SLS below. The
    # model variables come first, the instruments last.
    aux_names_b <- if (length(lavdata@ov.names.aux) >= b) {
      lavdata@ov.names.aux[[b]]
    } else {
      character(0L)
    }
    if (length(aux_names_b) > 0L && !is.null(lavdata@aux[[b]])) {
      xy <- cbind(xy, lavdata@aux[[b]])
      ov_names <- c(ov_names, aux_names_b)
    }

    # ingredients for the robust (Browne residual-based) overidentification
    # test: the distribution-free (ADF) ACOV of the sample (co)variances and
    # the sample covariance matrix of all variables in this block. This test
    # is valid under non-normality, where the classical Sargan test is not.
    # The moment block of the ADF Gamma is vech(S) (with the diagonal), so the
    # index map runs over the full vech (including variances).
    gamma_block <- NULL
    cov_block <- NULL
    moment_idx_mat <- NULL
    weights_b <- lavdata@weights[[b]]
    if (iv_sargan && nrow(xy) > ncol(xy)) {
      gamma_block <- try(lav_samp_gamma(
        m_y = xy, meanstructure = FALSE, wt = weights_b
      ), silent = TRUE)
      if (inherits(gamma_block, "try-error")) {
        gamma_block <- NULL
      } else {
        nvar_b <- ncol(xy)
        # ML sample covariance (divisor N), matching the ADF Gamma scale
        if (is.null(weights_b)) {
          cov_block <- crossprod(scale(xy, center = TRUE, scale = FALSE)) /
            nrow(xy)
        } else {
          cov_block <- stats::cov.wt(xy, wt = weights_b, method = "ML")$cov
        }
        # vech index map (column-major lower triangle, WITH the diagonal)
        moment_idx_mat <- matrix(0L, nvar_b, nvar_b)
        moment_idx_mat[lower.tri(moment_idx_mat, diag = TRUE)] <-
          seq_len(nvar_b * (nvar_b + 1L) / 2L)
        moment_idx_mat[upper.tri(moment_idx_mat)] <-
          t(moment_idx_mat)[upper.tri(moment_idx_mat)]
      }
    }

    # estimation per equation
    for (j in seq_along(eqs[[b]])) {
      # this equation
      eq <- eqs[[b]][[j]]

      # iv_flag?
      iv_flag <- TRUE
      if (!is.null(eq$iv_type)) {
        iv_flag <- eq$iv_type != "ols"
      } else {
        # if rhs_new matches miiv, iv_flag is FALSE
        if (identical(eq$rhs_new, eq$miiv)) {
          iv_flag <- FALSE
        }
      }

      # Y: there is always an y variable
      y_idx <- match(eq$lhs_new, ov_names)
      yvec <- xy[, y_idx, drop = TRUE] # y-variable (always scalar)

      # X: usually, there are x variables (apart from the "1")
      if (identical(eq$rhs_new, "1")) {
        x_idx <- integer(0L)
        xmat <- matrix(0, nrow = nrow(xy), ncol = 0L)
      } else {
        x_idx <- match(eq$rhs_new, ov_names)
        xmat <- xy[, x_idx, drop = FALSE]
      }

      # Z: instruments
      if (iv_flag) {
        i_idx <- match(eq$iv, ov_names)
        imat <- cbind(1, xy[, i_idx, drop = FALSE]) # instruments
      }

      # weights
      weights <- lavdata@weights[[b]]

      # sargan vector (classic test) and browne vector (robust test)
      sargan <- rep(as.numeric(NA), 3L)
      names(sargan) <- c("stat", "df", "pvalue")
      browne <- sargan

      # 0. check
      if (iv_flag && length(eq$iv) < length(eq$rhs)) {
        # what to do? skip, or proceed anyway?
        eqs[[b]][[j]]$coef <- numeric(0L)
        eqs[[b]][[j]]$nobs <- nrow(xmat)
        eqs[[b]][[j]]$df_res <- as.integer(NA)
        eqs[[b]][[j]]$resvar <- as.numeric(NA)
        eqs[[b]][[j]]$resvar_df_res <- as.numeric(NA)
        eqs[[b]][[j]]$XX <- crossprod(cbind(1, xmat))
        eqs[[b]][[j]]$vcov <- matrix(0, 0L, 0L)
        eqs[[b]][[j]]$sargan <- sargan
        eqs[[b]][[j]]$browne <- browne
        next
      }

      # 1. regress x on instruments
      if (iv_flag) {
        if (is.null(weights)) {
          fit_x_on_z <- lm.fit(x = imat, y = xmat)
        } else {
          fit_x_on_z <- lm.wfit(x = imat, y = xmat, weights = weights)
        }
        # check for NA's
        if (anyNA(fit_x_on_z$coefficients)) {
          lav_msg_warn(gettextf(
            "regression %s on %s failed (NAs);
                                 redundant instruments?",
            paste(ov_names[x_idx], collapse = " + "),
            paste(ov_names[i_idx], collapse = " + ")
          ))
        }
        xhat <- as.matrix(fit_x_on_z$fitted.values)
      } else {
        # just plain regression: x == xhat
        xhat <- xmat
      }

      # 2. regress y on xhat
      if (is.null(weights)) {
        fit_y_on_xhat <- lm.fit(x = cbind(1, xhat), y = yvec)
      } else {
        fit_y_on_xhat <- lm.wfit(
          x = cbind(1, xhat), y = yvec,
          weights = weights
        )
      }

      # store the (centered) slope system for the pooled solve that handles
      # simple equality constraints across equations (see below)
      if (length(x_idx) > 0L) {
        if (is.null(weights)) {
          x_bar <- colMeans(xhat)
          y_bar <- mean(yvec)
          xc <- sweep(xhat, 2L, x_bar)
          yc <- yvec - y_bar
          amat <- crossprod(xc)
          bvec <- crossprod(xc, yc)
        } else {
          sw <- sum(weights)
          x_bar <- colSums(weights * xhat) / sw
          y_bar <- sum(weights * yvec) / sw
          xc <- sweep(xhat, 2L, x_bar)
          yc <- yvec - y_bar
          amat <- crossprod(xc * weights, xc)
          bvec <- crossprod(xc * weights, yc)
        }
        eq_free_idx <- lavpartable$free[eq$pt]
        fix_idx <- which(eq_free_idx == 0L)
        if (length(fix_idx) > 0L) {
          b_fix <- lavpartable$ustart[eq$pt][fix_idx]
          b_fix[is.na(b_fix)] <- 0
          amat_free <- amat[-fix_idx, -fix_idx, drop = FALSE]
          bvec_free <- bvec[-fix_idx, , drop = FALSE] -
            amat[-fix_idx, fix_idx, drop = FALSE] %*% b_fix
          gcol_free <- eq_free_idx[-fix_idx]
        } else {
          amat_free <- amat
          bvec_free <- bvec
          gcol_free <- eq_free_idx
        }
        if (length(gcol_free) > 0L) {
          eqs[[b]][[j]]$slope_block <- list(
            amat = amat_free, bvec = drop(bvec_free), gcol = gcol_free,
            x_idx = x_idx, x_bar = x_bar, y_bar = y_bar
          )
        }
      }

      # 3. fill estimates in x
      # - eq$pt contains partable/user rows
      # - lavpartable$free[eq$pt] should give the free idx
      free_idx <- lavpartable$free[eq$pt]
      if (length(x_idx) > 0L) {
        if (all(free_idx > 0L)) {
          x[free_idx] <- fit_y_on_xhat$coefficients[-1]
        } else {
          # remove non-free elements
          zero_idx <- which(free_idx == 0L)
          free_idx <- free_idx[-zero_idx]
          if (length(free_idx) > 0) {
            x[free_idx] <- fit_y_on_xhat$coefficients[-1][-zero_idx]
          }
        }
      }
      if (lavmodel@meanstructure) {
        free_int_idx <- lavpartable$free[eq$ptint]
        if (free_int_idx > 0L) {
          x[free_int_idx] <- fit_y_on_xhat$coefficients[1]
        }
      }

      # 4. resvar
      resvar_df_res <- resvar <- df_res <- NULL
      if (iv_vcov_stage1 %in% c("lm.vcov.dfres", "lm.vcov")) {
        df_res <- fit_y_on_xhat$df.residual
        notna_idx <- unname(which(!is.na(fit_y_on_xhat$coefficients)))
        ycoef <- fit_y_on_xhat$coefficients[notna_idx]
        res <- yvec - drop(cbind(1, xmat)[, notna_idx, drop = FALSE] %*% ycoef)
        if (is.null(weights)) {
          sse <- sum(res * res)
        } else {
          sse <- sum(weights * res * res)
        }
        resvar_df_res <- sse / df_res # what we should do...
        resvar <- sse / length(res)
        # iv_mimic_ml: use the ML divisor (N) for the residual variance
        # used in the standard errors / Sargan test as well
        if (iv_mimic_ml) {
          resvar_df_res <- resvar
        }
      }

      # 5. naive cov (for standard errors) (see summary.lm in base R)
      vcov <- NULL
      if (iv_vcov_stage1 %in% c("lm.vcov.dfres", "lm.vcov")) {
        p1 <- 1L:fit_y_on_xhat$rank
        r <- chol2inv(fit_y_on_xhat$qr$qr[p1, p1, drop = FALSE])
        if (iv_vcov_stage1 == "lm.vcov") {
          vcov <- r * resvar
        } else {
          vcov <- r * resvar_df_res
        }
      }

      # 6. Sargan test (see summary.ivreg.R 363--371)
      if (iv_flag && iv_sargan && iv_vcov_stage1 != "none") {
        sargan["df"] <- length(i_idx) - length(x_idx)
        if (sargan["df"] > 0L) {
          if (is.null(weights)) {
            fit_yres_on_z <- lm.fit(x = imat, y = res)
            rssr <- sum((res - mean(res))^2)
            sse2 <- sum(fit_yres_on_z$residuals * fit_yres_on_z$residuals)
          } else {
            fit_yres_on_z <- lm.wfit(x = imat, y = res, weights = weights)
            rssr <- sum(weights * (res - weighted.mean(res, weights))^2)
            sse2 <- sum(weights * fit_yres_on_z$residuals *
              fit_yres_on_z$residuals)
          }
          sargan["stat"] <- length(res) * (1 - sse2 / rssr)
          sargan["pvalue"] <- pchisq(sargan["stat"], sargan["df"],
            lower.tail = FALSE
          )
        }
      }

      # 6b. robust residual-based (Browne) overidentification test, using the
      # distribution-free (ADF) ACOV. The residual g = cov(z, y - x beta) is
      # linear in the sample covariances; the helper extracts the relevant
      # sub-block of the ADF Gamma and forms the test. Skipped when an
      # instrument is redundant (NA slope) or the ADF Gamma is unavailable.
      if (iv_flag && iv_sargan && !is.null(gamma_block) &&
        length(i_idx) - length(x_idx) > 0L) {
        beta_x <- fit_y_on_xhat$coefficients[-1]
        if (!anyNA(beta_x)) {
          browne <- lav_sem_miiv_browne_test(
            y_idx = y_idx, x_idx = x_idx, i_idx = i_idx,
            beta = beta_x, sample_cov = cov_block,
            nacov = gamma_block, nth = 0L,
            moment_idx_mat = moment_idx_mat, nobs = nrow(xy)
          )
        }
      }

      # add info to eqs list
      eqs[[b]][[j]]$coef <- unname(fit_y_on_xhat$coefficients)
      eqs[[b]][[j]]$nobs <- nrow(xmat)
      eqs[[b]][[j]]$df_res <- df_res
      eqs[[b]][[j]]$resvar <- resvar
      eqs[[b]][[j]]$resvar_df_res <- resvar_df_res
      # eqs[[b]][[j]]$XX <- crossprod(cbind(1, xmat)) # needed? or only vcov?
      eqs[[b]][[j]]$vcov <- vcov
      eqs[[b]][[j]]$sargan <- sargan
      eqs[[b]][[j]]$browne <- browne
    } # eqs
  } # nblocks

  # pooled (system) solve for simple equality constraints among the directed
  # slopes, across ALL blocks/groups so that cross-group constraints (e.g.
  # equal loadings for measurement invariance) are pooled jointly rather than
  # left as last-write-wins (no-op when nothing is shared)
  x <- lav_sem_miiv_apply_directed_pool(
    eqs_b = do.call(c, eqs), x = x, lavmodel = lavmodel,
    lavpartable = lavpartable
  )

  # return theta1
  theta1 <- x[free_directed_idx]

  # add equations as an attribute
  attr(theta1, "eqs") <- eqs

  theta1
}

# stage 1: use 2SLS to find the regression coefficients of all
#          directed effects in the model -- samplestats version
# if samplestats = TRUE, x should contain the sample statistics
# (taken from lavh1$implied)
lav_sem_miiv_2sls_samp <- function(x = NULL, samplestats = FALSE,
                                          eqs = NULL, lavmodel = NULL,
                                          lavpartable = NULL, lavdata = NULL,
                                          lavsamplestats = NULL, lavh1 = NULL,
                                          iv_vcov_stage1 = "lm.vcov.dfres",
                                          iv_vcov_stage2 = "delta",
                                          iv_sargan = TRUE,
                                          free_directed_idx = NULL,
                                          aug_vec = NULL, aug_pdim = NULL,
                                          iv_mimic_ml = FALSE) {
  # no conditional.x for now!
  stopifnot(!lavmodel@conditional.x)
  iv_vcov_stage1 <- tolower(iv_vcov_stage1)
  stopifnot(iv_vcov_stage1 %in%
    c("lm.vcov.dfres", "lm.vcov", "gamma", "none"))

  # get lavpta
  lavpta <- lav_pt_attributes(lavpartable)

  if (is.null(eqs)) {
    # find ALL model-implied instrumental variables (miivs) per equation
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
  }

  # external instruments: 'aug_vec' mode. The directed coefficients are
  # computed from an AUGMENTED moment vector [model variables, instruments]
  # supplied directly (one segment per block, concatenated), instead of from
  # the model moments plus the raw-data instrument blocks. This lets a
  # numerical Jacobian differentiate the (constraint-aware) directed solve with
  # respect to the instrument moments as well (see lav_sem_miiv_vcov()).
  use_aug_vec <- !is.null(aug_vec)
  if (use_aug_vec) {
    aug_ms <- lavmodel@meanstructure
    aug_md <- if (aug_ms) {
      aug_pdim + aug_pdim * (aug_pdim + 1L) / 2L
    } else {
      aug_pdim * (aug_pdim + 1L) / 2L
    }
    aug_off <- cumsum(c(0L, aug_md))
  } else if (samplestats) {
    implied <- lav_vec_to_implied(x, lavmodel = lavmodel)
  } else {
    implied <- lavh1$implied
  }

  # number of blocks
  nblocks <- lavmodel@nblocks

  # parameter vector -> all NA (for the free parameters)
  x <- lav_model_get_parameters(lavmodel) * as.numeric(NA)

  # 2SLS per equation
  # for now: - no equality constraints (yet)
  #          - only OLS (not robust, no lasso, ...)
  for (b in seq_len(nblocks)) {
    # ov.names for this block
    ov_names <- lavpta$vnames$ov[[b]]

    # nobs for this 'block'
    nobs <- lavsamplestats@nobs[[b]]

    p_model <- length(ov_names)
    aug_nth <- NULL      # external instruments, categorical: augmented #th
    aug_th_idx <- NULL   # external instruments, categorical: augmented th.idx
    if (use_aug_vec) {
      # reconstruct the augmented covariance/mean directly from aug_vec
      seg <- aug_vec[(aug_off[b] + 1L):aug_off[b + 1L]]
      p_aug <- aug_pdim[b]
      aux_names_b <- if (length(lavdata@ov.names.aux) >= b) {
        lavdata@ov.names.aux[[b]]
      } else {
        character(0L)
      }
      ov_names <- c(ov_names, aux_names_b)
      has_ext_iv <- length(aux_names_b) > 0L
      if (aug_ms) {
        sample_mean <- seg[seq_len(p_aug)]
        sample_cov <- lav_matrix_vech_reverse(seg[-seq_len(p_aug)])
      } else {
        sample_cov <- lav_matrix_vech_reverse(seg)
        # no mean structure: the directed slopes do not depend on the means,
        # but the engine still needs a (held-fixed) mean vector for centering
        mm <- lavh1$implied$mean[[b]]
        if (is.null(mm)) {
          mm <- numeric(p_model)
        }
        instr_mean <- if (has_ext_iv && !is.null(lavdata@aux[[b]])) {
          colMeans(lavdata@aux[[b]], na.rm = TRUE)
        } else {
          numeric(0L)
        }
        sample_mean <- c(mm, instr_mean)
      }
    } else {
      # sample statistics for this block
      sample_cov <- implied$cov[[b]]
      sample_mean <- implied$mean[[b]]

      # external instruments (estimator = "IV"): augment the saturated sample
      # statistics with the extra observed variables (lavData@aux) so the 2SLS
      # estimator can use instruments that are not part of the model. The
      # model-variable block is kept exactly equal to the saturated statistics;
      # the instrument blocks come from the data. Model variables come first
      # and instruments last, so that the model-moment sub-block can be
      # recovered for the (model-dimensioned) Jacobians (k_mat) below.
      if (lavmodel@categorical) {
        # categorical: augmented polychoric correlation matrix (the instruments
        # must be ordinal); the slopes use the correlations only
        augc <- lav_sem_miiv_aug_cat_moments(
          lavdata = lavdata, b = b, model_cor = sample_cov, nobs = nobs
        )
        has_ext_iv <- !is.null(augc)
        if (has_ext_iv) {
          sample_cov <- augc$R
          ov_names <- augc$ov.names
          sample_mean <- c(sample_mean,
            rep(0, length(ov_names) - length(sample_mean)))
          aug_nth <- length(augc$th)
          aug_th_idx <- augc$th.idx
        }
      } else {
        aug <- lav_sem_miiv_aug_moments(
          lavdata = lavdata, b = b, ov_names = ov_names,
          sample_cov = sample_cov, sample_mean = sample_mean
        )
        has_ext_iv <- !is.null(aug)
        if (has_ext_iv) {
          sample_cov <- aug$cov
          sample_mean <- aug$mean
          ov_names <- aug$ov.names
        }
      }
    }

    # ingredients for the robust (Browne residual-based) overidentification
    # test (see lav_sem_miiv_browne_test()). It needs the ACOV (Gamma) of the
    # sample statistics plus a map from a variable pair to its position in the
    # moment block of that ACOV.
    #  - categorical: the moments are polychoric correlations and the classic
    #    Sargan test is not valid; the ACOV is the polychoric ACOV over
    #    [thresholds, correlations] (off-diagonal vech, eq_nth thresholds).
    #  - continuous: the robust test uses the distribution-free (ADF) ACOV of
    #    the sample covariances (vech WITH diagonal, eq_nth = 0), computed from
    #    the (complete) raw data; this is valid under non-normality.
    # 'eq_cov' is the sample covariance/correlation matrix on the SAME scale as
    # the ACOV, used to form the residual moments. Not computed in the 'aug_vec'
    # (numerical-differentiation) calls, where the results are discarded anyway.
    eq_nacov <- NULL
    eq_nth <- 0L
    moment_idx_mat <- NULL
    eq_cov <- sample_cov
    if (!use_aug_vec && iv_sargan) {
      nvar_b <- nrow(sample_cov)
      if (lavmodel@categorical) {
        moment_idx_mat <- matrix(0L, nvar_b, nvar_b)
        moment_idx_mat[lower.tri(moment_idx_mat)] <-
          seq_len(nvar_b * (nvar_b - 1L) / 2L)
        moment_idx_mat <- moment_idx_mat + t(moment_idx_mat)
        if (has_ext_iv) {
          eq_nacov <- augc$nacov
          eq_nth <- aug_nth
        } else {
          eq_nacov <- lavsamplestats@NACOV[[b]]
          eq_nth <- nrow(eq_nacov) - nvar_b * (nvar_b - 1L) / 2L
        }
      } else {
        # continuous: ADF Gamma from the (complete) raw data, ordered to match
        # ov_names (model variables first, external instruments appended)
        xy_b <- if (length(lavdata@X) >= b) lavdata@X[[b]] else NULL
        if (has_ext_iv && !is.null(lavdata@aux[[b]])) {
          xy_b <- cbind(xy_b, lavdata@aux[[b]])
        }
        if (!is.null(xy_b) && !anyNA(xy_b) && nrow(xy_b) > ncol(xy_b)) {
          wt_b <- lavdata@weights[[b]]
          g_try <- try(lav_samp_gamma(
            m_y = xy_b, meanstructure = FALSE, wt = wt_b
          ), silent = TRUE)
          if (!inherits(g_try, "try-error")) {
            eq_nacov <- g_try
            eq_nth <- 0L
            pvar <- ncol(xy_b)
            moment_idx_mat <- matrix(0L, pvar, pvar)
            moment_idx_mat[lower.tri(moment_idx_mat, diag = TRUE)] <-
              seq_len(pvar * (pvar + 1L) / 2L)
            moment_idx_mat[upper.tri(moment_idx_mat)] <-
              t(moment_idx_mat)[upper.tri(moment_idx_mat)]
            # ML sample covariance (divisor N), matching the ADF Gamma scale
            eq_cov <- if (is.null(wt_b)) {
              crossprod(scale(xy_b, center = TRUE, scale = FALSE)) / nrow(xy_b)
            } else {
              stats::cov.wt(xy_b, wt = wt_b, method = "ML")$cov
            }
          }
        }
      }
    }

    # estimation per equation
    for (j in seq_along(eqs[[b]])) {
      # this equation
      eq <- eqs[[b]][[j]]

      # iv_flag?
      iv_flag <- TRUE
      if (!is.null(eq$iv_type)) {
        iv_flag <- eq$iv_type != "ols"
      } else {
        # if rhs_new matches miiv, iv_flag is FALSE
        if (identical(eq$rhs_new, eq$miiv)) {
          iv_flag <- FALSE
        }
      }

      # Y: there is always an y variable
      y_idx <- match(eq$lhs_new, ov_names)
      y_bar <- sample_mean[y_idx]

      # X: usually, there are x variables (apart from the "1")
      if (identical(eq$rhs_new, "1")) {
        x_idx <- integer(0L)
      } else {
        x_idx <- match(eq$rhs_new, ov_names)
        x_bar <- sample_mean[x_idx]
      }
      nx <- length(x_idx)

      # Z: instruments
      s_xy <- sample_cov[x_idx, y_idx, drop = FALSE]
      s_xx <- sample_cov[x_idx, x_idx, drop = FALSE]
      s_yy <- sample_cov[y_idx, y_idx, drop = FALSE]
      nz <- 0L
      if (iv_flag) {
        i_idx <- match(eq$iv, ov_names)
        nz <- length(i_idx)
        s_xz <- sample_cov[x_idx, i_idx, drop = FALSE]
        s_zx <- sample_cov[i_idx, x_idx, drop = FALSE]
        s_zz <- sample_cov[i_idx, i_idx, drop = FALSE]
        s_zy <- sample_cov[i_idx, y_idx, drop = FALSE]
      }

      # sargan vector (classic test) and browne vector (robust test)
      sargan <- rep(as.numeric(NA), 3L)
      names(sargan) <- c("stat", "df", "pvalue")
      browne <- sargan

      # 0. check
      if (iv_flag && length(eq$iv) < length(eq$rhs)) {
        # what to do? skip, or proceed anyway?
        eqs[[b]][[j]]$coef <- numeric(0L)
        eqs[[b]][[j]]$nobs <- nobs
        eqs[[b]][[j]]$df_res <- as.integer(NA)
        eqs[[b]][[j]]$resvar <- as.numeric(NA)
        eqs[[b]][[j]]$resvar_df_res <- as.numeric(NA)
        eqs[[b]][[j]]$XX <- matrix(0, 0L, 0L)
        eqs[[b]][[j]]$vcov <- matrix(0, 0L, 0L)
        eqs[[b]][[j]]$sargan <- sargan
        eqs[[b]][[j]]$browne <- browne
        next
      }

      # 1. + 2. OLS + 2SLS
      fit_y_on_xhat <- list()
      if (iv_flag) {
        # Step 1: compute S_ZZ^{-1} S_ZX and S_ZZ^{-1} S_Zy
        m_w <- lav_mat_sym_solve_spd(s_zz, s_zx)
        v <- lav_mat_sym_solve_spd(s_zz, s_zy)
        # Step 2: build reduced system
        amat <- s_xz %*% m_w
        bvec <- s_xz %*% v
        beta_slopes <- drop(lav_mat_sym_solve_spd(amat, bvec))
        beta0 <- as.vector(y_bar - t(x_bar) %*% beta_slopes)
      } else {
        if (nx > 0L) {
          amat <- s_xx
          bvec <- s_xy
          beta_slopes <- drop(lav_mat_sym_solve_spd(amat, bvec))
          beta0 <- as.vector(y_bar - t(x_bar) %*% beta_slopes)
        } else {
          beta_slopes <- numeric(0L)
          beta0 <- y_bar
        }
      }
      fit_y_on_xhat$coefficients <- c(beta0, beta_slopes)

      # store the (centered) slope system for the pooled solve that handles
      # simple equality constraints across equations (see below). fixed slopes
      # (free == 0) are moved to the right-hand side.
      if (nx > 0L) {
        eq_free_idx <- lavpartable$free[eq$pt]
        fix_idx <- which(eq_free_idx == 0L)
        if (length(fix_idx) > 0L) {
          b_fix <- lavpartable$ustart[eq$pt][fix_idx]
          b_fix[is.na(b_fix)] <- 0
          amat_free <- amat[-fix_idx, -fix_idx, drop = FALSE]
          bvec_free <- bvec[-fix_idx, , drop = FALSE] -
            amat[-fix_idx, fix_idx, drop = FALSE] %*% b_fix
          gcol_free <- eq_free_idx[-fix_idx]
        } else {
          amat_free <- amat
          bvec_free <- bvec
          gcol_free <- eq_free_idx
        }
        if (length(gcol_free) > 0L) {
          eqs[[b]][[j]]$slope_block <- list(
            amat = amat_free, bvec = drop(bvec_free), gcol = gcol_free,
            x_idx = x_idx, x_bar = x_bar, y_bar = y_bar
          )
        }
      }

      # 2b. jac_eqs_k
      k_mat <- NULL
      k_mat_int <- NULL
      # full (augmented-moment) directed Jacobians, kept only when external
      # instruments are present (see the reduction below and the stage-2
      # standard errors in lav_sem_miiv_vcov())
      k_mat_full <- NULL
      k_mat_int_full <- NULL
      if (iv_vcov_stage1 == "gamma" || iv_vcov_stage2 == "h2") {
        nvar <- nrow(sample_cov)
        pstar <- nvar * (nvar + 1) / 2
        # Ex[k, j] = 1 iff x.idx[k] == j
        ex <- matrix(0.0, nx, nvar)
        if (nx > 0L) {
          ex[cbind(seq_len(nx), x_idx)] <- 1.0
        }
        # ---- Jacobian of beta_slopes w.r.t. vech(S) --------------
        a_vec <- lav_mat_vech_row_idx(nvar)
        b_vec <- lav_mat_vech_col_idx(nvar)
        offdiag <- (a_vec > b_vec) # logical: strictly off-diagonal elements

        if (!iv_flag && nx == 0L) {
          j_slopes_s <- matrix(0.0, nrow = 0L, ncol = pstar)
        } else if (!iv_flag) {
          g <- -as.vector(crossprod(ex, beta_slopes))
          g[y_idx] <- g[y_idx] + 1.0
          # build M (p x m) column-by-column in vectorised form
          m <- sweep(ex[, b_vec, drop = FALSE], 2L, g[a_vec], `*`) +
            sweep(ex[, a_vec, drop = FALSE], 2L, g[b_vec] * offdiag, `*`)
          j_slopes_s <- lav_mat_sym_solve_spd(s_xx, m)
        } else {
          r_z <- drop(s_zy) - s_zx %*% beta_slopes
          u <- drop(lav_mat_sym_solve_spd(s_zz, r_z))
          ez <- matrix(0.0, nz, nvar)
          ez[cbind(seq_len(nz), i_idx)] <- 1.0
          ew <- crossprod(m_w, ez)
          cu <- as.vector(crossprod(ez, u))
          cb_x <- as.vector(crossprod(ex, beta_slopes))
          h <- -cu
          h[y_idx] <- h[y_idx] + 1.0
          # build M (p x m) vectorised
          m <-
            sweep(ex[, b_vec, drop = FALSE], 2L, cu[a_vec], `*`) +
            sweep(ew[, a_vec, drop = FALSE], 2L, h[b_vec], `*`) -
            sweep(ew[, b_vec, drop = FALSE], 2L, cb_x[a_vec], `*`) +
            sweep(ex[, a_vec, drop = FALSE], 2L, cu[b_vec] * offdiag, `*`) +
            sweep(ew[, b_vec, drop = FALSE], 2L, h[a_vec] * offdiag, `*`) -
            sweep(ew[, a_vec, drop = FALSE], 2L, cb_x[b_vec] * offdiag, `*`)
          j_slopes_s <- solve(amat, m)
        }

        # fill in k_mat_int and k_mat
        k_mat_int <- numeric(0L)
        if (!lavmodel@categorical && lavmodel@meanstructure) {
          k_mat_int <- numeric(nvar + pstar)
          k_mat_int[y_idx] <- 1.0
          if (nx > 0L) {
            k_mat_int[x_idx] <- -beta_slopes
            k_mat_int[nvar + seq_len(pstar)] <- -x_bar %*% j_slopes_s
          }
          k_mat <- matrix(0.0, nx, nvar + pstar)
          if (nx > 0L) {
            k_mat[, nvar + seq_len(pstar)] <- j_slopes_s
          }
        } else if (lavmodel@categorical) {
          # split: var_num, off
          # FIXME: all ordinal for now (ie no variances)
          tmp <- j_slopes_s[, -lav_mat_diagh_idx(nvar), drop = FALSE]
          # number of thresholds: for external instruments the moments span the
          # augmented [thresholds, correlations] set (aug_nth thresholds), else
          # the model NACOV size minus the correlation block
          nth <- if (has_ext_iv) {
            aug_nth
          } else {
            nrow(lavsamplestats@NACOV[[b]]) - ncol(tmp)
          }
          k_mat_int <- numeric(nth + ncol(tmp))
          k_mat <- matrix(0.0, nx, nth + ncol(tmp))
          if (nx > 0L) {
            k_mat[, nth + seq_len(ncol(tmp))] <- tmp
          }
        } else {
          k_mat <- matrix(0.0, nrow = nx, ncol = pstar)
          if (nx > 0L) {
            k_mat <- j_slopes_s
          }
        }

        # external instruments: the Jacobians here are over the AUGMENTED
        # moments [model variables, instruments]. Keep the full (augmented)
        # versions for the stage-2 standard errors (so the instrument sampling
        # variance is folded into the variance/covariance parameter standard
        # errors, see lav_sem_miiv_vcov()), and also produce model-moment-only
        # versions (the model variables come first) for the model-dimensioned
        # consumers that do not augment (e.g. lav_sem_miiv_utils_jack_eqs()).
        if (has_ext_iv && lavmodel@categorical) {
          k_mat_full <- k_mat
          k_mat_int_full <- k_mat_int
          keep <- lav_sem_miiv_aug_cat_keep_idx(p_model, nvar, aug_th_idx)
          k_mat <- k_mat[, keep, drop = FALSE]
          k_mat_int <- k_mat_int[keep]
        } else if (has_ext_iv && !lavmodel@categorical) {
          k_mat_full <- k_mat
          k_mat_int_full <- k_mat_int
          if (lavmodel@meanstructure) {
            keep <- lav_aux_moment_idx(p_model, nvar)
            k_mat <- k_mat[, keep, drop = FALSE]
            k_mat_int <- k_mat_int[keep]
          } else {
            ri <- lav_mat_vech_row_idx(nvar)
            ci <- lav_mat_vech_col_idx(nvar)
            vech_keep <- which(ri <= p_model & ci <= p_model)
            k_mat <- k_mat[, vech_keep, drop = FALSE]
          }
        }
        k_mat
      }

      # 3. fill estimates in x
      # - eq$pt contains partable/user rows
      # - lavpartable$free[eq$pt] should give the free idx
      free_idx <- lavpartable$free[eq$pt]
      if (nx > 0L) {
        if (all(free_idx > 0L)) {
          x[free_idx] <- fit_y_on_xhat$coefficients[-1]
        } else {
          # remove non-free elements
          zero_idx <- which(free_idx == 0L)
          free_idx <- free_idx[-zero_idx]
          if (length(free_idx) > 0) {
            x[free_idx] <- fit_y_on_xhat$coefficients[-1][-zero_idx]
          }
        }
      }
      if (lavmodel@meanstructure) {
        free_int_idx <- lavpartable$free[eq$ptint]
        if (free_int_idx > 0L) {
          x[free_int_idx] <- fit_y_on_xhat$coefficients[1]
        }
      }

      # 4. resvar
      resvar_df_res <- resvar <- df_res <- NULL
      if (iv_vcov_stage1 %in% c("lm.vcov.dfres", "lm.vcov") || iv_sargan) {
        if (nx > 0L) {
          tmp <- s_xx %*% beta_slopes
          resvar <- as.numeric(s_yy - 2 * t(beta_slopes) %*% s_xy +
            t(beta_slopes) %*% tmp)
        } else {
          resvar <- drop(s_yy)
        }
        df_res <- nobs - (nx + 1L)
        resvar_df_res <- resvar * nobs / df_res
        # iv_mimic_ml: use the ML divisor (N) instead of N-P, so the residual
        # variance used in the standard errors / Sargan test stays on the same
        # (ML) scale as the rescaled sample covariance (resvar = RSS / N)
        if (iv_mimic_ml) {
          resvar_df_res <- resvar
        }
      }

      # 5. naive cov (for standard errors) (see summary.lm in base R)
      vcov <- NULL
      if (iv_vcov_stage1 %in% c("lm.vcov.dfres", "lm.vcov")) {
        this_resvar <- resvar_df_res
        if (iv_vcov_stage1 == "lm.vcov") {
          this_resvar <- resvar
        }
        if (nx > 0L) {
          ainv <-
            lav_mat_sym_solve_spd(amat, diag(x = 1, nrow = nrow(amat)))
          vcov_slopes <- (this_resvar / nobs) * ainv
          vcov_beta0 <-
            (this_resvar / nobs) + t(x_bar) %*% vcov_slopes %*% x_bar
          cov_beta0_slopes <- -vcov_slopes %*% x_bar
          vcov <- matrix(0, nx + 1L, nx + 1L)
          vcov[1, 1] <- vcov_beta0
          vcov[1, -1] <- t(cov_beta0_slopes)
          vcov[-1, 1] <- cov_beta0_slopes
          vcov[-1, -1] <- vcov_slopes
        } else {
          # only intercept
          vcov <- matrix(this_resvar / nobs, 1L, 1L)
        }
      }

      # 6. overidentification test(s)
      if (iv_flag && iv_sargan && iv_vcov_stage1 != "none") {
        # classic Sargan test (see summary.ivreg.R 363--371); for categorical
        # data the p-value is left NA because the classical Sargan test is not
        # valid for polychoric correlations (use the Browne test below instead)
        sargan["df"] <- nz - nx
        if (sargan["df"] > 0L) {
          g <- s_zy - s_zx %*% beta_slopes
          w <- lav_mat_sym_solve_spd(s_zz, g)
          sargan["stat"] <- as.numeric(nobs * t(g) %*% w / resvar_df_res)
          if (!lavmodel@categorical) {
            sargan["pvalue"] <- pchisq(sargan["stat"], sargan["df"],
              lower.tail = FALSE
            )
          }
        }
        # robust residual-based (Browne) test using the ACOV of the sample
        # statistics (polychoric ACOV for categorical data, distribution-free
        # ADF ACOV for continuous data); only when that ACOV is available (it is
        # not in the 'aug_vec' numerical-differentiation calls, nor for missing
        # data, where eq_nacov stays NULL)
        if (!use_aug_vec && !is.null(eq_nacov)) {
          browne <- lav_sem_miiv_browne_test(
            y_idx = y_idx, x_idx = x_idx, i_idx = i_idx,
            beta = beta_slopes, sample_cov = eq_cov,
            nacov = eq_nacov, nth = eq_nth,
            moment_idx_mat = moment_idx_mat, nobs = nobs
          )
        }
      }

      # XX (needed?)
      # XX <- NULL
      #       if (!lavmodel@categorical) {
      #         if (nx > 0L) {
      #           XX <- tcrossprod(c(1, x.bar)) * nobs
      #           XX[-1, -1] <- (S_XX + tcrossprod(x.bar)) * nobs
      #         } else {
      #           XX <- matrix(nobs, 1L, 1L)
      #         }
      #       }

      # add info to eqs list
      eqs[[b]][[j]]$coef <- unname(fit_y_on_xhat$coefficients)
      eqs[[b]][[j]]$nobs <- nobs
      eqs[[b]][[j]]$df_res <- df_res
      eqs[[b]][[j]]$resvar <- resvar
      eqs[[b]][[j]]$resvar_df_res <- resvar_df_res
      # eqs[[b]][[j]]$XX <- XX
      eqs[[b]][[j]]$vcov <- vcov
      eqs[[b]][[j]]$sargan <- sargan
      eqs[[b]][[j]]$browne <- browne
      eqs[[b]][[j]]$k_mat_int <- k_mat_int
      eqs[[b]][[j]]$k_mat <- k_mat
      eqs[[b]][[j]]$k_mat_full <- k_mat_full
      eqs[[b]][[j]]$k_mat_int_full <- k_mat_int_full
    } # eqs
  } # nblocks

  # pooled (system) solve for simple equality constraints among the directed
  # slopes: if a free column is shared by >1 equation (within a group, or
  # across groups for measurement invariance), the per-equation fills above are
  # last-write-wins and must be replaced by a joint solution over ALL
  # blocks/groups (no-op when nothing is shared).
  x <- lav_sem_miiv_apply_directed_pool(
    eqs_b = do.call(c, eqs), x = x, lavmodel = lavmodel,
    lavpartable = lavpartable
  )

  # return theta1
  theta1 <- x[free_directed_idx]

  # add equations as an attribute
  if (!samplestats) {
    attr(theta1, "eqs") <- eqs
  }

  theta1
}


# second stage: compute theta_2 as follows:
# theta_2 = solve(t(Delta2) %*% W2 %*% Delta2) %*% t(Delta2) %*% W2 %*% svec
#
# this is the fast but single-block, no-equality-constraints version
# (for now)
#
# see lav_sem_miiv_varcov_old() for an alternative
# model derivative matrix Delta (per block), with the ceq.simple shared
# columns collapsed so that the columns are indexed by the *free* parameter
# number (matching free_directed_idx / free_undirected_idx). Under
# ceq.simple.only, lav_model_delta() returns nx.unco columns (one per
# unconstrained parameter); indexing those by free numbers would select the
# wrong columns and corrupt the second-stage (co)variance estimates (e.g. the
# residual-variance blow-up seen with equal factor loadings).
lav_sem_miiv_delta <- function(lavmodel, glist = NULL) {
  delta <- lav_model_delta(lavmodel = lavmodel, glist = glist)
  if (lavmodel@ceq.simple.only) {
    delta <- lapply(delta, function(d) d %*% lavmodel@ceq.simple.K)
  }
  delta
}

# solve the undirected (variance/covariance) parameters for a single block
# 'b', given that block's collapsed Delta and the free-parameter indices 'fu'
# of its undirected parameters. ULS/GLS use a single WLS linearization; 2RLS
# and RLS reweight with the model-implied covariance of block 'b'. Returns
# theta2 for 'fu' (with the H attribute when return_h = TRUE).
lav_sem_miiv_varcov_block <- function(b, fu, delta_b, implied, lavmodel, x,
                                      iv_varcov_method, return_h) {
  delta2 <- delta_b[, fu, drop = FALSE]

  # general linear equality constraints among the undirected parameters (if
  # any); aligned to the 'fu' (ie delta2 column) order
  con_und <- lav_sem_miiv_linear_con(lavmodel, fu)
  con_jac <- if (is.null(con_und)) NULL else con_und$jac
  con_rhs <- if (is.null(con_und)) NULL else con_und$rhs

  # svec
  sample_cov <- implied$cov[[b]]
  sample_mean <- implied$mean[[b]]
  if (lavmodel@categorical) {
    svec <- c(implied$th[[b]], lav_mat_vech(sample_cov, diagonal = FALSE))
  } else if (lavmodel@meanstructure) {
    svec <- c(sample_mean, lav_mat_vech(sample_cov))
  } else {
    svec <- lav_mat_vech(sample_cov)
  }

  # initial estimate for theta.2
  if (iv_varcov_method == "GLS") {
    s <- sample_cov
  } else { # ULS, or starting point for 2RLS/RLS
    s <- diag(1, nrow = nrow(sample_cov))
  }
  theta2 <- lav_utils_wls_linearization(delta = delta2, s = s,
                                        meanstructure = lavmodel@meanstructure,
                                        categorical = lavmodel@categorical,
                                        svec = svec, return_h = return_h,
                                        con_jac = con_jac, con_rhs = con_rhs)

  # 2RLS
  if (iv_varcov_method == "2RLS") {
    x[fu] <- theta2
    lavmodel_tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
    sigma <- lav_model_sigma(lavmodel = lavmodel_tmp, extra = FALSE)[[b]]
    # fall back to the ULS estimate if the model-implied Sigma (from the ULS
    # step) is not positive definite (eg a Heywood case)
    new_theta2 <- try(
      lav_utils_wls_linearization(delta = delta2, s = sigma,
                                  meanstructure = lavmodel@meanstructure,
                                  svec = svec, return_h = return_h,
                                  con_jac = con_jac, con_rhs = con_rhs),
      silent = TRUE
    )
    if (inherits(new_theta2, "try-error")) {
      lav_msg_warn(gettext("2RLS for variances/covariances fell back to the
                            ULS estimate: the model-implied covariance matrix
                            is not positive definite."))
    } else {
      theta2 <- new_theta2
    }
  # RLS
  } else if (iv_varcov_method == "RLS") {
    # typically, we need 5-15 iterations, so max = 200 should be enough
    for (i in seq_len(200L)) {
      old_x <- theta2
      x[fu] <- theta2
      lavmodel_tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
      sigma <- lav_model_sigma(lavmodel = lavmodel_tmp, extra = FALSE)[[b]]
      # the RLS weight needs a positive-definite model-implied Sigma; if the
      # current iterate yields a non-PD Sigma (eg a Heywood case), keep the
      # last valid estimate rather than aborting the whole estimator
      new_theta2 <- try(
        lav_utils_wls_linearization(delta = delta2, s = sigma,
                                    meanstructure = lavmodel@meanstructure,
                                    svec = svec, return_h = return_h,
                                    con_jac = con_jac, con_rhs = con_rhs),
        silent = TRUE
      )
      if (inherits(new_theta2, "try-error")) {
        lav_msg_warn(gettext("RLS for variances/covariances stopped early:
                              the model-implied covariance matrix is not
                              positive definite; returning the last valid
                              estimate."))
        theta2 <- old_x
        break
      }
      theta2 <- new_theta2
      sse <- Re(sum((old_x - theta2)^2))
      if (sse < 1e-12 * Re(1 + sum(theta2^2))) {
        break
      } else if (i == 200L) {
        lav_msg_warn(gettext("RLS for variances/covariances did not converge
                              after 200 iterations."))
      }
    }
  }

  theta2
}

lav_sem_miiv_varcov <- function(x = NULL, samplestats = FALSE,
                                lavmodel = NULL, lavpartable = NULL,
                                lavsamplestats = NULL,
                                lavh1 = NULL, free_directed_idx = NULL,
                                free_undirected_idx = NULL,
                                add_h2 = FALSE,
                                iv_varcov_method = "RLS") {

  # x = NULL?
  if (is.null(x)) {
    # calling from vcov?
    x <- lav_model_get_parameters(lavmodel)
    lavmodel_tmp <- lavmodel
    implied <- lavh1$implied
  } else if (samplestats) {
    # only used to compute jac_b numerically
    implied <- lav_vec_to_implied(x, lavmodel = lavmodel)
    x <- lav_model_get_parameters(lavmodel)
    lavmodel_tmp <- lavmodel
  } else {
    # first time, or to compute jac_a numerically
    theta1 <- x
    x <- lav_model_get_parameters(lavmodel)
    x[free_directed_idx] <- theta1
    lavmodel_tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
    implied <- lavh1$implied
  }

  # number of blocks
  nblocks <- lavmodel@nblocks
  delta_list <- lav_sem_miiv_delta(lavmodel_tmp)

  # single group: solve directly (also keeps the H attribute for the SE path)
  if (nblocks == 1L) {
    return(lav_sem_miiv_varcov_block(
      b = 1L, fu = free_undirected_idx, delta_b = delta_list[[1L]],
      implied = implied, lavmodel = lavmodel, x = x,
      iv_varcov_method = iv_varcov_method, return_h = add_h2))
  }

  # multiple groups: each group's moments contribute its own block; solve the
  # undirected parameters per block and assemble in free_undirected_idx order.
  # (Cross-group equality constraints among loadings/intercepts are imposed in
  # stage 1; cross-group constraints among variances, if any, are handled per
  # block here.)
  theta2 <- numeric(length(free_undirected_idx))
  for (b in seq_len(nblocks)) {
    fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
      lavpartable$free > 0L & lavpartable$block == b])
    fu_b <- fu_b[fu_b %in% free_undirected_idx]
    if (length(fu_b) == 0L) {
      next
    }
    tb <- lav_sem_miiv_varcov_block(
      b = b, fu = fu_b, delta_b = delta_list[[b]], implied = implied,
      lavmodel = lavmodel, x = x, iv_varcov_method = iv_varcov_method,
      return_h = FALSE)
    theta2[match(fu_b, free_undirected_idx)] <- as.numeric(tb)
  }

  theta2
}


# restricted-2SLS covariance for the directed parameters under (simple)
# equality constraints. Each equation contributes its per-equation 2SLS
# coefficient information (the inverse of eq$vcov) to the shared free-parameter
# positions; summing the information across equations that share a free
# parameter yields the variance-weighted restricted-2SLS covariance used by
# MIIVsem. This is algebraically the restricted-GLS covariance K (K' A K)^{-1}
# K' (with A the block-diagonal per-equation information and K the sharing
# map), and reproduces the top-left block of MIIVsem's KKT system
# [[A, R'], [R, 0]]^{-1}. Returns the covariance block aligned to
# free_directed_idx (in the compact / nx.free space).
lav_sem_miiv_directed_vcov_restricted <- function(eqs, lavpartable,
                                                  free_directed_idx, nfree) {
  info <- matrix(0, nfree, nfree)
  for (b in seq_along(eqs)) {
    for (eq in eqs[[b]]) {
      if (is.null(eq$vcov) || nrow(eq$vcov) == 0L) {
        next
      }
      # coefficient -> free parameter index; eq$vcov is ordered as
      # (intercept, slopes), matching c(eq$ptint, eq$pt)
      int_free <- if (!is.null(eq$ptint) && length(eq$ptint) > 0L) {
        lavpartable$free[eq$ptint]
      } else {
        0L
      }
      slope_free <- if (length(eq$rhs) > 0L) {
        lavpartable$free[eq$pt]
      } else {
        integer(0L)
      }
      coef_free <- c(int_free, slope_free)
      keep <- which(coef_free > 0L)
      if (length(keep) == 0L) {
        next
      }
      a_j <- solve(eq$vcov[keep, keep, drop = FALSE])
      fidx <- coef_free[keep]
      info[fidx, fidx] <- info[fidx, fidx] + a_j
    }
  }
  solve(info[free_directed_idx, free_directed_idx, drop = FALSE])
}


# K Gamma K' over multiple blocks with the block-diagonal NORMAL-THEORY
# moment covariance, without ever materializing the (pstar x pstar) Gamma
# blocks: per block the K-sandwich form of lav_mat_k_gammant_kt is used
# (Gamma_b = NT gamma of block b, divided by its nobs). For larger models
# building the explicit Gamma dominates the whole vcov computation.
lav_sem_miiv_k_gammant_kt_blocks <- function(k_mat = NULL, cov_list = NULL,
                                             md = NULL, moff = NULL,
                                             meanstructure = FALSE,
                                             x_idx_list = NULL,
                                             nobs = NULL) {
  out <- matrix(0, nrow = nrow(k_mat), ncol = nrow(k_mat))
  for (b in seq_along(md)) {
    k_b <- k_mat[, moff[b] + seq_len(md[b]), drop = FALSE]
    out <- out + lav_mat_k_gammant_kt(
      m_k = k_b, s = cov_list[[b]],
      meanstructure = meanstructure, x_idx = x_idx_list[[b]]
    ) / nobs[[b]]
  }
  out
}


# VCOV for free parameters
lav_sem_miiv_vcov <- function(lavmodel = NULL, lavsamplestats = NULL,
                              lavoptions = NULL, lavpartable = NULL,
                              lavimplied = NULL,
                              lavh1 = NULL, lavdata = NULL, eqs = NULL) {
  # iv options
  iv_vcov_stage1 <- tolower(lavoptions$estimator.args$iv_vcov_stage1)
  iv_vcov_stage2 <- tolower(lavoptions$estimator.args$iv_vcov_stage2)
  # iv_samplestats <- lavoptions$estimator.args$iv_samplestats
  iv_vcov_gamma_modelbased <-
    lavoptions$estimator.args$iv_vcov_gamma_modelbased
  iv_varcov_method <- toupper(lavoptions$estimator.args$iv_varcov_method)
  iv_vcov_jack_numerical <- lavoptions$estimator.args$iv_vcov_jack_numerical
  iv_vcov_jaca_numerical <- lavoptions$estimator.args$iv_vcov_jaca_numerical
  # iv_vcov_jacb_numerical <- lavoptions$estimator.args$iv_vcov_jacb_numerical

  # empty vcov
  vcov <- matrix(0, lavmodel@nx.free, lavmodel@nx.free)

  # nblocks
  nblocks <- lavmodel@nblocks

  # eqs?
  if (is.null(eqs)) {
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
  }

  # the raw-data engine (iv_samplestats = FALSE) produces the per-equation
  # covariances (eq$vcov) used for the directed standard errors, but not the
  # analytic moment Jacobian (eq$k_mat / eq$k_mat_full). When those are
  # missing, the directed Jacobian needed for the stage-2 (and gamma) standard
  # errors is computed numerically instead (from the sample statistics, which
  # does not require k_mat).
  eqs_have_kmat <- FALSE
  eqs_have_kmat_full <- FALSE
  for (eqs_b in eqs) {
    for (e in eqs_b) {
      if (!is.null(e$k_mat)) {
        eqs_have_kmat <- TRUE
      }
      if (!is.null(e$k_mat_full)) {
        eqs_have_kmat_full <- TRUE
      }
    }
  }
  if (!eqs_have_kmat) {
    iv_vcov_jack_numerical <- TRUE
  }

  # directed versus undirected (free) parameters
  undirected_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  directed_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~", "~1")) # WITH intercepts!!
  directed_noint_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~")) # NO intercepts!!
  if (lavmodel@categorical) {
    th_idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "|")
    if (length(unlist(lavmodel@num.idx)) > 0L) {
      lav_msg_warn(gettext(
        "[IV] no standard errors (yet) when variables are both ordinal
         and continuous"
      ))
      iv_vcov_stage1 <- "none"
      iv_vcov_stage2 <- "none"
    }
  }

  free_directed_idx <- unique(lavpartable$free[directed_idx])
  free_directed_noint_idx <- unique(lavpartable$free[directed_noint_idx])
  free_undirected_idx <- unique(lavpartable$free[undirected_idx])
  if (lavmodel@categorical) {
    free_th_idx <- unique(lavpartable$free[th_idx])
  }

  # directed free parameters only
  x <- lav_model_get_parameters(lavmodel, type = "free")
  # theta1 <- x[free_directed_idx]
  theta1_noint <- x[free_directed_noint_idx]
  theta2 <- x[free_undirected_idx] # final estimate

  # equality constraints? Two flavors:
  #  - simple (parameters sharing a free column, ceq.simple.only path)
  #  - general linear (rows of ceq.JAC)
  # The analytic Jacobians assume independent equations / no constraints, so
  # when a constraint touches the directed (undirected) parameters we fall back
  # to the numerical Jacobian of the constrained point-estimation function,
  # which differentiates the actual pooled / projected solve.
  dir_slope_free <- lavpartable$free[lavpartable$free > 0L &
    lavpartable$op %in% c("=~", "~")]
  has_shared_dir <- anyDuplicated(dir_slope_free) > 0L
  con_dir <- lav_sem_miiv_linear_con(lavmodel, free_directed_noint_idx)
  con_und <- lav_sem_miiv_linear_con(lavmodel, free_undirected_idx)
  has_dir_con <- has_shared_dir || !is.null(con_dir)
  has_und_con <- !is.null(con_und)

  # how to handle equality constraints among the *directed* parameters:
  #  - lm.vcov(.dfres) (default): variance-weighted restricted-2SLS covariance
  #    for the directed block (matches MIIVsem) -- see
  #    lav_sem_miiv_directed_vcov_restricted()
  #  - gamma: the original delta-method (moment-Jacobian) standard errors, with
  #    a numerical Jacobian of the constrained pooled solve
  # Either way stage 2 still needs the constraint-aware (numerical) directed
  # Jacobian to propagate the directed uncertainty into the undirected block.
  use_restricted_directed <- FALSE
  if (has_dir_con && length(free_directed_idx) > 0L &&
      iv_vcov_stage1 != "none") {
    if (iv_vcov_stage1 %in% c("lm.vcov", "lm.vcov.dfres")) {
      if (isTRUE(lavoptions$estimator.args$iv_samplestats)) {
        # the restricted-2SLS covariance only handles *simple* equality
        # constraints (parameters that share a free column, ceq.simple.only);
        # for general linear directed constraints (rows of ceq.JAC) fall back
        # to the delta method
        if (lavmodel@ceq.simple.only && has_shared_dir && is.null(con_dir)) {
          use_restricted_directed <- TRUE
        } else {
          iv_vcov_stage1 <- "gamma"
        }
      } else {
        lav_msg_warn(gettext(
          "[IV] equality constraints are ignored by the raw-data standard
           errors; use iv_samplestats = TRUE for constraint-aware
           standard errors."))
      }
    }
    iv_vcov_jack_numerical <- TRUE
  }

  # general linear constraints among the undirected parameters: the analytic
  # stage-2 Jacobians (and the h2 weight) ignore the constraint, so switch to
  # the delta method with numerical jac_a / jac_b (which differentiate the
  # constrained projection in lav_sem_miiv_varcov())
  # (iv_vcov_jaca_numerical gates both the numerical jac_a and jac_b below)
  if (has_und_con && iv_vcov_stage2 != "none") {
    iv_vcov_stage2 <- "delta"
    iv_vcov_jaca_numerical <- TRUE
  }

  # switch off iv_vcov_stage2?
  if (length(free_directed_idx) > 0L && iv_vcov_stage1 == "none") {
    iv_vcov_stage2 <- "none"
  }

  # ---- multiple groups (continuous data) ----------------------------------
  # The sample moments are independent across groups, so the moment covariance
  # is block-diagonal. We build that block-diagonal 'gamma_big' (each block its
  # NT moment covariance divided by its nobs), a block-diagonal stage-2 map h2,
  # the (already multiblock) numerical directed Jacobian jac_k, and the
  # restricted-2SLS directed covariance (which loops over all blocks/equations
  # and so handles cross-group constraints). The ceq.simple expansion at the
  # end maps the compact vcov to the parameter table.
  if (nblocks > 1L && iv_vcov_stage1 != "none") {
    two_stage <- any(lavoptions$missing == c("two.stage", "robust.two.stage"))
    # two-stage missing data and categorical data need the explicit moment
    # covariance for the directed block; complete continuous data uses the
    # (exact) restricted-2SLS directed covariance
    use_explicit_gamma <- two_stage || lavmodel@categorical
    vec_list <- lav_implied_to_vec(implied = lavh1$implied,
                                   lavmodel = lavmodel, drop_list = FALSE)
    md <- vapply(vec_list, length, integer(1L)) # moment dim per block
    moff <- cumsum(c(0L, md))                    # column offsets
    ntot_moment <- sum(md)

    # block-diagonal moment covariance gamma_big = bdiag(Gamma_b / nobs_b),
    # where Gamma_b is the (per-observation) asymptotic covariance of the
    # sample moments of group b: the categorical ADF gamma (NACOV) for ordered
    # data, or the two-stage gamma (Savalei & Bentler / Savalei & Falk) under
    # two-stage missing data. For complete continuous data (the normal-theory
    # gamma) gamma_big is never materialized: the K Gamma K' sandwiches below
    # use the per-block K-form of lav_mat_k_gammant_kt directly (see
    # lav_sem_miiv_k_gammant_kt_blocks), which only needs the per-block
    # covariance matrices (nt_cov_b)
    x_idx_list <- vector("list", nblocks)
    nt_cov_b <- vector("list", nblocks)
    gamma_big <- NULL
    if (use_explicit_gamma) {
      gamma_g <- vector("list", nblocks)
    }
    for (b in seq_len(nblocks)) {
      x_idx_b <-
        if (lavmodel@fixed.x) lavsamplestats@x.idx[[b]] else integer(0L)
      x_idx_list[[b]] <- x_idx_b
      if (lavmodel@categorical) {
        gamma_g[[b]] <- lavsamplestats@NACOV[[b]] / lavsamplestats@nobs[[b]]
      } else if (two_stage) {
        mu_b <- lavh1$implied$mean[[b]]
        sig_b <- lavh1$implied$cov[[b]]
        if (identical(tolower(lavoptions$missing), "robust.two.stage")) {
          gg <- lav_mvn_mi_h1_omega_sw(y = lavdata@X[[b]], mp = lavdata@Mp[[b]],
            yp = lavsamplestats@missing[[b]], mu = mu_b, sigma_1 = sig_b,
            x_idx = x_idx_b, information = "observed")
        } else {
          i1 <- lav_mvnorm_missing_information_observed_samplestats(
            yp = lavsamplestats@missing[[b]], mu = mu_b, sigma_1 = sig_b,
            x_idx = x_idx_b)
          gg <- lav_mat_sym_inverse(i1)
        }
        gamma_g[[b]] <- gg / lavsamplestats@nobs[[b]]
      } else {
        cov_b <- if (iv_vcov_gamma_modelbased) lavimplied$cov[[b]] else
          lavh1$implied$cov[[b]]
        if (inherits(try(chol(cov_b), silent = TRUE), "try-error")) {
          cov_b <- lavh1$implied$cov[[b]]
        }
        nt_cov_b[[b]] <- cov_b
      }
    }
    if (use_explicit_gamma) {
      gamma_big <- lav_mat_bdiag(gamma_g)
    }

    # external instruments with CATEGORICAL data (multiple groups): replace the
    # per-block model gamma (NACOV) and the directed Jacobian with their
    # augmented [thresholds, correlations] versions, so the instrument sampling
    # variability is folded into all standard errors. The augmented directed
    # Jacobian is assembled analytically from k_mat_full (the numerical route is
    # not valid for categorical data).
    ext_cat_mg <- lavmodel@categorical && length(free_directed_idx) > 0L &&
      length(unlist(lavdata@ov.names.aux)) > 0L && eqs_have_kmat_full
    md_aug <- NULL
    moff_aug <- NULL
    cat_keep_b <- NULL
    if (ext_cat_mg) {
      md_aug <- integer(nblocks)
      cat_keep_b <- vector("list", nblocks)
      gamma_aug_b <- vector("list", nblocks)
      for (b in seq_len(nblocks)) {
        augb <- lav_sem_miiv_aug_cat_moments(lavdata = lavdata, b = b,
          model_cor = lavh1$implied$cov[[b]], model_th = lavsamplestats@th[[b]],
          nobs = lavsamplestats@nobs[[b]])
        gamma_aug_b[[b]] <- augb$nacov / lavsamplestats@nobs[[b]]
        cat_keep_b[[b]] <- lav_sem_miiv_aug_cat_keep_idx(
          lavmodel@nvar[[b]], ncol(augb$R), augb$th.idx)
        md_aug[b] <- nrow(augb$nacov)
      }
      moff_aug <- cumsum(c(0L, md_aug))
      gamma_big <- lav_mat_bdiag(gamma_aug_b)
    }

    # numerical directed Jacobian over the stacked moments (jac_k differentiates
    # the per-group 2SLS solve, including the cross-group pooled constraints)
    jac_k <- NULL
    if (length(free_directed_idx) > 0L && ext_cat_mg) {
      # categorical external instruments: analytic augmented Jacobian
      ntot_aug <- sum(md_aug)
      jac_k <- matrix(0, length(free_directed_idx), ntot_aug)
      for (b in seq_len(nblocks)) {
        jb <- lav_sem_miiv_jack_eqs_aug(eqs = eqs, block = b,
          lavmodel = lavmodel, lavpartable = lavpartable,
          free_directed_idx = free_directed_idx)
        jac_k[, moff_aug[b] + seq_len(md_aug[b])] <- jb
      }
    } else if (length(free_directed_idx) > 0L) {
      vec <- unlist(vec_list)
      jac_k <- numDeriv::jacobian(lav_sem_miiv_2sls_samp, x = vec,
        samplestats = TRUE, eqs = eqs, lavmodel = lavmodel,
        lavpartable = lavpartable, lavsamplestats = lavsamplestats,
        lavh1 = lavh1, lavdata = lavdata,
        free_directed_idx = free_directed_idx)
    }

    # directed covariance. Complete continuous data: restricted-2SLS
    # (block-diagonal for unconstrained groups, pooled across groups for
    # cross-group constraints). Two-stage missing data: the explicit
    # jac_k %*% gamma %*% t(jac_k), since the per-equation covariances do not
    # carry the EM-uncertainty correction.
    if (length(free_directed_idx) > 0L) {
      v_dir <- NULL
      if (!use_explicit_gamma) {
        v_dir <- try(lav_sem_miiv_directed_vcov_restricted(eqs = eqs,
          lavpartable = lavpartable, free_directed_idx = free_directed_idx,
          nfree = lavmodel@nx.free), silent = TRUE)
      }
      if (!is.null(v_dir) && !inherits(v_dir, "try-error")) {
        vcov[free_directed_idx, free_directed_idx] <- v_dir
      } else if (use_explicit_gamma) {
        vcov[free_directed_idx, free_directed_idx] <-
          jac_k %*% gamma_big %*% t(jac_k)
      } else {
        # continuous NT gamma: per-block K-sandwich, no explicit gamma
        vcov[free_directed_idx, free_directed_idx] <-
          lav_sem_miiv_k_gammant_kt_blocks(
            k_mat = jac_k, cov_list = nt_cov_b, md = md, moff = moff,
            meanstructure = lavmodel@meanstructure,
            x_idx_list = x_idx_list, nobs = lavsamplestats@nobs
          )
      }
    }

    # undirected covariance: block-diagonal stage-2 map, propagating the
    # directed uncertainty via jac_k
    if (iv_vcov_stage2 != "none" && length(free_undirected_idx) > 0L &&
        iv_varcov_method != "NONE") {
      delta_list <- lav_sem_miiv_delta(lavmodel)
      h2_full <- matrix(0, length(free_undirected_idx), ntot_moment)
      delta1_full <- matrix(0, ntot_moment, length(free_directed_idx))
      for (b in seq_len(nblocks)) {
        cols <- moff[b] + seq_len(md[b])
        if (length(free_directed_idx) > 0L) {
          delta1_full[cols, ] <-
            delta_list[[b]][, free_directed_idx, drop = FALSE]
        }
        fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
          lavpartable$free > 0L & lavpartable$block == b])
        fu_b <- fu_b[fu_b %in% free_undirected_idx]
        if (length(fu_b) == 0L) {
          next
        }
        delta2_b <- delta_list[[b]][, fu_b, drop = FALSE]
        if (lavmodel@categorical || iv_varcov_method == "ULS") {
          s_mat <- diag(lavmodel@nvar[[b]])
        } else if (iv_varcov_method == "GLS") {
          s_mat <- lavh1$implied$cov[[b]]
        } else {
          s_mat <- lavimplied$cov[[b]]
          if (inherits(try(chol(s_mat), silent = TRUE), "try-error")) {
            s_mat <- lavh1$implied$cov[[b]]
          }
        }
        tmp_b <- lav_utils_wls_linearization(delta = delta2_b, s = s_mat,
          meanstructure = lavmodel@meanstructure,
          categorical = lavmodel@categorical,
          svec = vec_list[[b]], return_h = TRUE)
        h2_b <- attr(tmp_b, "H")
        # the stage-2 map covers the (co)variance/correlation moments; for
        # categorical data the threshold moments come first in each block and
        # are left at zero (the (co)variance parameters do not depend on them)
        nth_b <- md[b] - ncol(h2_b)
        h2_cols <- moff[b] + nth_b + seq_len(ncol(h2_b))
        h2_full[match(fu_b, free_undirected_idx), h2_cols] <- h2_b
      }
      # external instruments (continuous data): fold the instrument sampling
      # variance into the variance/covariance parameter standard errors, by
      # rebuilding the stage-2 map over the per-block AUGMENTED moments
      # [model variables, instruments] and sandwiching with the block-diagonal
      # normal-theory gamma of the augmented covariances (see the single-group
      # path for the rationale). The directed coefficients' dependence on the
      # instrument moments enters through the full (augmented) directed
      # Jacobian assembled analytically from k_mat_full.
      ext_aux_mg <- !use_explicit_gamma && length(free_directed_idx) > 0L &&
        length(unlist(lavdata@ov.names.aux)) > 0L
      if (ext_aux_mg) {
        ms <- lavmodel@meanstructure
        md_aug <- integer(nblocks)
        cov_aug_b <- vector("list", nblocks)
        keep_b <- vector("list", nblocks)
        for (b in seq_len(nblocks)) {
          augb <- lav_sem_miiv_aug_moments(lavdata = lavdata, b = b,
            ov_names = lavdata@ov.names[[b]],
            sample_cov = lavh1$implied$cov[[b]],
            sample_mean = if (ms) lavh1$implied$mean[[b]] else NULL)
          cov_aug_b[[b]] <- augb$cov
          p_model_b <- lavmodel@nvar[[b]]
          p_aug_b <- ncol(augb$cov)
          keep_b[[b]] <- lav_sem_miiv_aug_keep_idx(p_model_b, p_aug_b, ms)
          md_aug[b] <- if (ms) {
            p_aug_b + p_aug_b * (p_aug_b + 1L) / 2L
          } else {
            p_aug_b * (p_aug_b + 1L) / 2L
          }
        }
        moff_aug <- cumsum(c(0L, md_aug))
        ntot_aug <- sum(md_aug)
        jac_k_aug <- matrix(0, length(free_directed_idx), ntot_aug)
        delta1_aug <- matrix(0, ntot_aug, length(free_directed_idx))
        h2_aug <- matrix(0, length(free_undirected_idx), ntot_aug)
        for (b in seq_len(nblocks)) {
          cols_b <- moff_aug[b] + seq_len(md_aug[b])
          model_cols_b <- moff_aug[b] + keep_b[[b]]
          if (!has_dir_con && eqs_have_kmat_full) {
            jb <- lav_sem_miiv_jack_eqs_aug(eqs = eqs, block = b,
              lavmodel = lavmodel, lavpartable = lavpartable,
              free_directed_idx = free_directed_idx)
            jac_k_aug[, cols_b] <- jb
          }
          delta1_aug[model_cols_b, ] <-
            delta_list[[b]][, free_directed_idx, drop = FALSE]
          fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
            lavpartable$free > 0L & lavpartable$block == b])
          fu_b <- fu_b[fu_b %in% free_undirected_idx]
          if (length(fu_b) == 0L) {
            next
          }
          delta2_b <- delta_list[[b]][, fu_b, drop = FALSE]
          if (iv_varcov_method == "ULS") {
            s_mat <- diag(lavmodel@nvar[[b]])
          } else if (iv_varcov_method == "GLS") {
            s_mat <- lavh1$implied$cov[[b]]
          } else {
            s_mat <- lavimplied$cov[[b]]
            if (inherits(try(chol(s_mat), silent = TRUE), "try-error")) {
              s_mat <- lavh1$implied$cov[[b]]
            }
          }
          tmp_b <- lav_utils_wls_linearization(delta = delta2_b, s = s_mat,
            meanstructure = ms, categorical = FALSE,
            svec = vec_list[[b]], return_h = TRUE)
          h2_b <- attr(tmp_b, "H")
          nth_b <- md[b] - ncol(h2_b) # mean offset (continuous)
          h2_b <- cbind(matrix(0, nrow(h2_b), nth_b), h2_b)
          h2_aug[match(fu_b, free_undirected_idx), model_cols_b] <- h2_b
        }
        if (has_dir_con || !eqs_have_kmat_full) {
          # use the numerical augmented Jacobian when the analytic k_mat_full
          # is unavailable (raw-data engine) or when there are equality
          # constraints among the directed coefficients (the per-equation
          # analytic Jacobian cannot represent the cross-equation/cross-group
          # pooling)
          jac_k_aug <- lav_sem_miiv_jac_full_numeric(
            eqs = eqs, lavmodel = lavmodel, lavpartable = lavpartable,
            lavdata = lavdata, lavsamplestats = lavsamplestats,
            lavh1 = lavh1, free_directed_idx = free_directed_idx)
        }
        tmp <- h2_aug %*% (diag(ntot_aug) - delta1_aug %*% jac_k_aug)
        # continuous NT gamma over the augmented moments: per-block
        # K-sandwich, no explicit gamma
        vcov[free_undirected_idx, free_undirected_idx] <-
          lav_sem_miiv_k_gammant_kt_blocks(
            k_mat = tmp, cov_list = cov_aug_b, md = md_aug, moff = moff_aug,
            meanstructure = ms, x_idx_list = x_idx_list,
            nobs = lavsamplestats@nobs
          )
      } else if (ext_cat_mg) {
        # categorical external instruments: build the stage-2 map and directed
        # derivative over the per-block augmented [thresholds, correlations]
        # layout; gamma_big and jac_k are already augmented (see above)
        ntot_aug <- sum(md_aug)
        delta1_aug <- matrix(0, ntot_aug, length(free_directed_idx))
        h2_aug <- matrix(0, length(free_undirected_idx), ntot_aug)
        for (b in seq_len(nblocks)) {
          model_cols_b <- moff_aug[b] + cat_keep_b[[b]]
          delta1_aug[model_cols_b, ] <-
            delta_list[[b]][, free_directed_idx, drop = FALSE]
          fu_b <- unique(lavpartable$free[lavpartable$op == "~~" &
            lavpartable$free > 0L & lavpartable$block == b])
          fu_b <- fu_b[fu_b %in% free_undirected_idx]
          if (length(fu_b) == 0L) {
            next
          }
          delta2_b <- delta_list[[b]][, fu_b, drop = FALSE]
          s_mat <- diag(lavmodel@nvar[[b]])
          tmp_b <- lav_utils_wls_linearization(delta = delta2_b, s = s_mat,
            meanstructure = lavmodel@meanstructure, categorical = TRUE,
            svec = vec_list[[b]], return_h = TRUE)
          h2_b <- attr(tmp_b, "H")
          nth_b <- md[b] - ncol(h2_b)
          h2_b <- cbind(matrix(0, nrow(h2_b), nth_b), h2_b)
          h2_aug[match(fu_b, free_undirected_idx), model_cols_b] <- h2_b
        }
        tmp <- h2_aug %*% (diag(ntot_aug) - delta1_aug %*% jac_k)
        vcov[free_undirected_idx, free_undirected_idx] <-
          tmp %*% gamma_big %*% t(tmp)
      } else {
        if (length(free_directed_idx) > 0L) {
          tmp <- h2_full %*% (diag(ntot_moment) - delta1_full %*% jac_k)
        } else {
          tmp <- h2_full
        }
        if (use_explicit_gamma) {
          vcov[free_undirected_idx, free_undirected_idx] <-
            tmp %*% gamma_big %*% t(tmp)
        } else {
          # continuous NT gamma: per-block K-sandwich, no explicit gamma
          vcov[free_undirected_idx, free_undirected_idx] <-
            lav_sem_miiv_k_gammant_kt_blocks(
              k_mat = tmp, cov_list = nt_cov_b, md = md, moff = moff,
              meanstructure = lavmodel@meanstructure,
              x_idx_list = x_idx_list, nobs = lavsamplestats@nobs
            )
        }
      }
    }

    # thresholds (categorical): the threshold parameters equal the sample
    # thresholds, so their covariance is the corresponding block of the moment
    # covariance, delta_th' gamma delta_th
    if (lavmodel@categorical && length(free_th_idx) > 0L) {
      delta_list_th <- lav_sem_miiv_delta(lavmodel)
      delta_full <- do.call("rbind", delta_list_th)
      delta_th <- delta_full[, free_th_idx, drop = FALSE]
      gamma_th <- gamma_big
      if (ext_cat_mg) {
        # map the per-block model threshold derivative into the augmented layout
        delta_th <- matrix(0, sum(md_aug), length(free_th_idx))
        for (b in seq_len(nblocks)) {
          cols <- moff[b] + seq_len(md[b])
          delta_th[moff_aug[b] + cat_keep_b[[b]], ] <-
            delta_list_th[[b]][, free_th_idx, drop = FALSE]
        }
        gamma_th <- gamma_big
      }
      zero_idx <- which(rowSums(abs(delta_th)) == 0)
      if (length(zero_idx) > 0L) {
        delta_th <- delta_th[-zero_idx, , drop = FALSE]
        gamma_th <- gamma_th[-zero_idx, -zero_idx, drop = FALSE]
      }
      vcov[free_th_idx, free_th_idx] <-
        t(delta_th) %*% gamma_th %*% delta_th
    }

    # ceq.simple: expand the compact vcov to the 'unco' space
    if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L) {
      vcov <- lavmodel@ceq.simple.K %*% vcov %*% t(lavmodel@ceq.simple.K)
    }
    return(vcov)
  }

  # missing data: two-stage standard errors. Under missing = "two.stage" /
  # "robust.two.stage" the sample moments are the (saturated) EM estimates;
  # their asymptotic covariance ('gamma') is no longer the complete-data NT
  # gamma but the inverse saturated information (two.stage) or its sandwich
  # (robust.two.stage). We build that 'gamma_big' explicitly (in the same
  # (mean, vech(cov)) ordering as jac_k) and route the SEs through the same
  # explicit jac %*% gamma %*% t(jac) path used for categorical data.
  two_stage <- any(lavoptions$missing == c("two.stage", "robust.two.stage"))
  use_explicit_gamma <- lavmodel@categorical || two_stage

  # do we need gamma_big?
  if (iv_vcov_stage1 == "gamma" ||
     (iv_vcov_stage2 != "none" && length(free_undirected_idx) > 0L)) {
      gamma_g <- vector("list", lavmodel@ngroups)
    for (g in seq_len(lavmodel@ngroups)) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      # (model-implied or h1) covariance, used by the continuous NT-gamma
      # path; needed whether or not NACOV is available
      if (iv_vcov_gamma_modelbased) {
        cov_g <- lavimplied$cov[[g]]
        # the model-implied covariance may be non-PD (eg a Heywood case); the
        # continuous NT-gamma path needs a PD matrix, so fall back to the
        # (PD) sample covariance in that case
        if (!lavmodel@categorical &&
            inherits(try(chol(cov_g), silent = TRUE), "try-error")) {
          cov_g <- lavh1$implied$cov[[g]]
        }
      } else {
        cov_g <- lavh1$implied$cov[[g]]
      }
      if (two_stage) {
        # two-stage moment covariance, in (mean, vech(cov)) order
        mu_g <- lavh1$implied$mean[[g]]
        sig_g <- lavh1$implied$cov[[g]]
        x_idx_g <-
          if (lavmodel@fixed.x) lavsamplestats@x.idx[[g]] else integer(0L)
        if (identical(tolower(lavoptions$missing), "robust.two.stage")) {
          # sandwich I_1^{-1} J_1 I_1^{-1} (Savalei & Falk, 2014)
          gg <- lav_mvn_mi_h1_omega_sw(
            y = lavdata@X[[g]], mp = lavdata@Mp[[g]],
            yp = lavsamplestats@missing[[g]],
            mu = mu_g, sigma_1 = sig_g, x_idx = x_idx_g,
            information = "observed"
          )
        } else {
          # I_1^{-1} (Savalei & Bentler, 2009)
          i1 <- lav_mvnorm_missing_information_observed_samplestats(
            yp = lavsamplestats@missing[[g]],
            mu = mu_g, sigma_1 = sig_g, x_idx = x_idx_g
          )
          gg <- lav_mat_sym_inverse(i1)
        }
        gamma_g[[g]] <- fg * gg
      } else if (!is.null(lavsamplestats@NACOV[[g]])) {
        gamma_g[[g]] <- fg * lavsamplestats@NACOV[[g]]
      }
    }
    if (two_stage || !is.null(lavsamplestats@NACOV[[g]])) {
      gamma_big <- lav_mat_bdiag(gamma_g)
    }
  }

  # compute jac_k if needed
  if (length(free_directed_idx) > 0L && iv_vcov_stage1 != "none"
      && (iv_vcov_stage1 == "gamma" || iv_vcov_stage2 == "h2")) {
    if (iv_vcov_jack_numerical) {
      # compute K matrix
      vec <- lav_implied_to_vec(
        implied = lavh1$implied, lavmodel = lavmodel,
        drop_list = TRUE
      )
      jac_k <- numDeriv::jacobian(
        lav_sem_miiv_2sls_samp,
        x = vec, samplestats = TRUE, eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavsamplestats = lavsamplestats,
        lavh1 = lavh1, lavdata = lavdata,
        free_directed_idx = free_directed_idx
      )
    } else {
      # analytic version -- # FIXME: multiple blocks!!!
      stopifnot(lavmodel@nblocks == 1L)
      jac_k <- lav_sem_miiv_utils_jack_eqs(
        eqs = eqs, block = 1L,
        lavmodel = lavmodel, lavpartable = lavpartable,
        free_directed_idx = free_directed_idx
      )
    }
  }

  # external instruments with CATEGORICAL data (single group): replace the
  # model gamma (NACOV) and the directed Jacobian with their augmented
  # [thresholds, correlations] versions, so the instrument sampling variability
  # is folded into all standard errors (directed coefficients and the stage-2
  # parameters). The augmented directed Jacobian is assembled analytically from
  # k_mat_full (the numerical route reconstructs a Pearson covariance and is
  # not valid for categorical data; equality constraints among directed
  # coefficients are therefore not supported here).
  ext_cat <- lavmodel@categorical && lavmodel@nblocks == 1L &&
    length(unlist(lavdata@ov.names.aux)) > 0L &&
    length(free_directed_idx) > 0L && eqs_have_kmat_full
  cat_keep <- NULL
  if (ext_cat) {
    augc <- lav_sem_miiv_aug_cat_moments(
      lavdata = lavdata, b = 1L, model_cor = lavh1$implied$cov[[1]],
      model_th = lavsamplestats@th[[1]], nobs = lavsamplestats@nobs[[1]])
    # single-group explicit-gamma path uses fg * NACOV (fg = nobs/ntotal = 1),
    # i.e. the augmented NACOV on the same scale as lavsamplestats@NACOV
    gamma_big <- (lavsamplestats@nobs[[1]] / lavsamplestats@ntotal) * augc$nacov
    jac_k <- lav_sem_miiv_jack_eqs_aug(
      eqs = eqs, block = 1L, lavmodel = lavmodel,
      lavpartable = lavpartable, free_directed_idx = free_directed_idx)
    cat_keep <- lav_sem_miiv_aug_cat_keep_idx(
      lavmodel@nvar[[1]], ncol(augc$R), augc$th.idx)
  }

  # stage 1: directed effects
  if (length(free_directed_idx) > 0L && iv_vcov_stage1 != "none") {
    if (use_restricted_directed) {
      # variance-weighted restricted-2SLS covariance (matches MIIVsem)
      v_dir <- try(lav_sem_miiv_directed_vcov_restricted(
        eqs = eqs, lavpartable = lavpartable,
        free_directed_idx = free_directed_idx, nfree = lavmodel@nx.free),
        silent = TRUE)
      if (!inherits(v_dir, "try-error")) {
        vcov[free_directed_idx, free_directed_idx] <- v_dir
      }
    } else if (iv_vcov_stage1 == "gamma") {
      if (use_explicit_gamma) {
        # explicit moment covariance: categorical (ADF) or two-stage missing
        k_gamma_adf_kt <- jac_k %*% gamma_big %*% t(jac_k)
        vcov[free_directed_idx, free_directed_idx] <-
          k_gamma_adf_kt / lavsamplestats@ntotal
      } else {
        # vcov = K %*% Gamma_NT %*% t(K)
        # k_gammant_kt <- jac_k %*% gamma_big %*% t(jac_k)
        x_idx <-
          if (lavmodel@fixed.x) lavsamplestats@x.idx[[1]] else integer(0L)
        if (length(unlist(lavdata@ov.names.aux)) > 0L) {
          # external instruments: the 2SLS slopes depend on the instrument
          # moments (and not on the model moments), so the model-moment
          # Jacobian would be ~0. Use the full (augmented) directed Jacobian
          # and the augmented covariance instead.
          jac_full <- lav_sem_miiv_jac_full_numeric(
            eqs = eqs, lavmodel = lavmodel, lavpartable = lavpartable,
            lavdata = lavdata, lavsamplestats = lavsamplestats,
            lavh1 = lavh1, free_directed_idx = free_directed_idx)
          cov_aug <- lav_sem_miiv_aug_moments(
            lavdata = lavdata, b = 1L, ov_names = lavdata@ov.names[[1]],
            sample_cov = lavh1$implied$cov[[1]],
            sample_mean = if (lavmodel@meanstructure) {
              lavh1$implied$mean[[1]]
            } else {
              NULL
            })$cov
          k_gammant_kt <- lav_mat_k_gammant_kt(m_k = jac_full, s = cov_aug,
            meanstructure = lavmodel@meanstructure, x_idx = x_idx)
        } else {
          k_gammant_kt <- lav_mat_k_gammant_kt(m_k = jac_k, s = cov_g,
            meanstructure = lavmodel@meanstructure, x_idx = x_idx)
        }
        vcov[free_directed_idx, free_directed_idx] <-
          k_gammant_kt / lavsamplestats@ntotal
      }

    # iv_vcov_stage1 = lm.vcov or lm.vcov.dfres
    } else {
      for (b in seq_len(nblocks)) {
        neqs <- length(eqs[[b]])
        for (j in seq_len(neqs)) {
          eq <- eqs[[b]][[j]]

          eq_vcov <- NULL
          if (!is.null(eq$vcov) && nrow(eq$vcov) > 0L) {
            eq_vcov <- eq$vcov # always includes intercept
          } else if (!is.null(eq$XX) && !is.null(eq$resvar)) {
            # reconstruct using XX and resvar (less stable!)
            eq_vcov <- solve(crossprod(eq$XX)) * eq$resvar # or resvar_df_res
          }
          if (is.null(eq_vcov)) {
            next
          }

          # only keep free parameters
          free_idx <- lavpartable$free[eq$pt]
          if (length(eq$rhs) > 0L) {
            if (all(free_idx > 0L)) {
              vcov[free_idx, free_idx] <- eq_vcov[-1, -1]
            } else {
              # remove non-free elements
              zero_idx_1 <- which(free_idx == 0L)
              free_idx <- free_idx[-zero_idx_1]
              if (length(free_idx) > 0L) {
                vcov[free_idx, free_idx] <-
                  eq_vcov[-1, -1][-zero_idx_1, -zero_idx_1, drop = FALSE]
              }
            }
          }
          if (lavmodel@meanstructure) {
            free_int_idx <- lavpartable$free[eq$ptint]
            if (free_int_idx > 0L) {
              vcov[free_int_idx, free_int_idx] <- eq_vcov[1, 1]
            }
          }
        } # neqs
      } # nblocks
    } # continuous/vcov version
  } # stage 1

  # stage 2: undirected effects (note: intercepts have no effect!)

  if (iv_vcov_stage2 != "none" && length(free_undirected_idx) > 0L &&
    iv_varcov_method != "NONE") {
    if (iv_vcov_stage2 == "h2") {
      # delta_list <- lav_model_delta(lavmodel) # per block/group

      # # block-wise
      # delta1_block <- vector("list", length = nblocks)
      # delta2_block <- vector("list", length = nblocks)
      # w2_block <- vector("list", length = nblocks)
      # for (b in seq_len(nblocks)) {
      #   nvar <- lavmodel@nvar[b]
      #   # delta1 + delta2
      #   delta1_block[[b]] <- delta_list[[b]][, free_directed_idx,
      #                                                  drop = FALSE]
      #   delta2_block[[b]] <-
      #     delta_list[[b]][, free_undirected_idx, drop = FALSE]
      #   # w2
      #   if (lavmodel@categorical || iv_varcov_method == "ULS") {
      #     s.inv <- diag(1, nrow = nvar)
      #   } else if (iv_varcov_method == "GLS") {
      #     s.inv <- lavh1$implied$cov[[b]]
      #   } else {
      #     s.inv <- lavimplied$cov[[b]]
      #   }
      #   w2_22 <- 0.5 * lav_mat_dup_pre_post(s.inv %x% s.inv)
      #   if (lavmodel@meanstructure) {
      #     w2_block[[b]] <- lav_mat_bdiag(s.inv, w2_22)
      #   } else {
      #     w2_block[[b]] <- w2_22
      #   }
      # } # blocks
      # delta1 <- do.call("rbind", delta1_block)
      # delta2 <- do.call("rbind", delta2_block)
      # w2 <- lav_mat_bdiag(w2_block)
      # h2 <- solve(t(delta2) %*% w2 %*% delta2, t(delta2) %*% w2)

      # tmp <- lav_sem_miiv_varcov(
      #   x = NULL, lavmodel = lavmodel,
      #   lavpartable = lavpartable, lavsamplestats = lavsamplestats,
      #   lavh1 = lavh1, free_directed_idx = free_directed_idx,
      #   free_undirected_idx = free_undirected_idx, add_h2 = TRUE,
      #   iv_varcov_method = iv_varcov_method
      # )
      # h2 <- attr(tmp, "H")
      stopifnot(lavmodel@nblocks == 1L)
      delta <- lav_sem_miiv_delta(lavmodel)[[1L]]
      delta2 <- delta[, free_undirected_idx, drop = FALSE]
      svec <- lav_implied_to_vec(
          implied = lavh1$implied, lavmodel = lavmodel,
          drop_list = TRUE
      )
      if (lavmodel@categorical || iv_varcov_method == "ULS") {
        s_mat <- diag(lavmodel@nvar[[1]])
      } else if (iv_varcov_method == "GLS") {
        s_mat <- lavh1$implied$cov[[1]]
      } else {
        # 2RLS/RLS use the model-implied covariance, which may be non-PD (eg a
        # Heywood case); fall back to the (PD) sample covariance if so
        s_mat <- lavimplied$cov[[1]]
        if (!lavmodel@categorical &&
            inherits(try(chol(s_mat), silent = TRUE), "try-error")) {
          s_mat <- lavh1$implied$cov[[1]]
        }
      }
      tmp <- lav_utils_wls_linearization(delta = delta2, s = s_mat,
        meanstructure = lavmodel@meanstructure,
        categorical = lavmodel@categorical, svec = svec, return_h = TRUE)
      h2 <- attr(tmp, "H")
      if (length(free_directed_idx) > 0L) {
        # assuming a SINGLE block for now
        delta1 <- delta[, free_directed_idx, drop = FALSE]
        nth <- nrow(delta1) - ncol(h2)
        h2 <- cbind(matrix(0, nrow(h2), ncol = nth), h2)
        if (ext_cat) {
          # place the (model-moment) stage-2 map and directed derivative into
          # the augmented [thresholds, correlations] layout, then propagate the
          # directed uncertainty through the augmented Jacobian
          m_aug <- ncol(jac_k)
          h2_aug <- matrix(0, nrow(h2), m_aug)
          h2_aug[, cat_keep] <- h2
          d1_aug <- matrix(0, m_aug, ncol(delta1))
          d1_aug[cat_keep, ] <- delta1
          tmp <- h2_aug %*% (diag(m_aug) - d1_aug %*% jac_k)
        } else {
          idiag <- diag(1, nrow = nrow(delta1))
          tmp <- h2 %*% (idiag - delta1 %*% jac_k)
        }
      } else {
        tmp <- h2
      }
      if (use_explicit_gamma) {
        tmp_gamma_tmpt <- tmp %*% gamma_big %*% t(tmp)
        vcov[free_undirected_idx, free_undirected_idx] <-
          tmp_gamma_tmpt / lavsamplestats@ntotal
        # add thresholds (categorical only)
        nth <- if (lavmodel@categorical) length(free_th_idx) else 0L
        if (length(nth) > 0L && nth > 0L) {
          delta_th <- delta[, free_th_idx, drop = FALSE]
          if (ext_cat) {
            # map the model threshold derivative into the augmented layout
            dta <- matrix(0, ncol(jac_k), ncol(delta_th))
            dta[cat_keep, ] <- delta_th
            delta_th <- dta
            gamma_th <- gamma_big
          } else {
            gamma_th <- gamma_big
          }
          # FIXME: remove zero rows from delta
          zero_idx <- which(rowSums(abs(delta_th)) == 0)
          if (length(zero_idx) > 0L) {
            delta_th <- delta_th[-zero_idx, , drop = FALSE]
            gamma_th <- gamma_th[-zero_idx, -zero_idx, drop = FALSE]
          }
          th_gamma_tht <- t(delta_th) %*% gamma_th %*% delta_th
          vcov[free_th_idx, free_th_idx] <-
            th_gamma_tht / lavsamplestats@ntotal
        }
      } else {
        x_idx <-
          if (lavmodel@fixed.x) lavsamplestats@x.idx[[1]] else integer(0L)
        # external instruments: fold the instrument sampling variance into the
        # variance/covariance parameter standard errors. The directed
        # coefficients depend on the instrument moments, so we propagate the
        # uncertainty through the FULL (augmented-moment) directed Jacobian and
        # sandwich with the normal-theory gamma of the augmented covariance,
        # rather than treating the instrument moments as fixed.
        jac_full <- NULL
        if (length(unlist(lavdata@ov.names.aux)) > 0L &&
            length(free_directed_idx) > 0L) {
          if (has_dir_con || !eqs_have_kmat_full) {
            # use the numerical augmented Jacobian when the analytic k_mat_full
            # is unavailable (raw-data engine) or when there are equality
            # constraints among the directed coefficients (the per-equation
            # analytic Jacobian cannot represent the cross-equation pooling)
            jac_full <- lav_sem_miiv_jac_full_numeric(
              eqs = eqs, lavmodel = lavmodel, lavpartable = lavpartable,
              lavdata = lavdata, lavsamplestats = lavsamplestats,
              lavh1 = lavh1, free_directed_idx = free_directed_idx)
          } else {
            jac_full <- lav_sem_miiv_jack_eqs_aug(
              eqs = eqs, block = 1L, lavmodel = lavmodel,
              lavpartable = lavpartable, free_directed_idx = free_directed_idx)
          }
        }
        if (!is.null(jac_full)) {
          # the instrument moments have no model-implied counterpart, so for a
          # coherent augmented covariance (needed by the normal-theory gamma)
          # we use the saturated (H1) sample covariance for the whole matrix
          aug <- lav_sem_miiv_aug_moments(
            lavdata = lavdata, b = 1L, ov_names = lavdata@ov.names[[1]],
            sample_cov = lavh1$implied$cov[[1]],
            sample_mean = if (lavmodel@meanstructure) {
              lavh1$implied$mean[[1]]
            } else {
              NULL
            })
          cov_aug <- aug$cov
          p_model <- lavmodel@nvar[[1]]
          p_aug <- ncol(cov_aug)
          keep <- lav_sem_miiv_aug_keep_idx(p_model, p_aug,
            lavmodel@meanstructure)
          m_aug <- ncol(jac_full)
          h2_aug <- matrix(0, nrow(h2), m_aug)
          h2_aug[, keep] <- h2
          d1_aug <- matrix(0, m_aug, ncol(delta1))
          d1_aug[keep, ] <- delta1
          tmp_aug <- h2_aug %*% (diag(m_aug) - d1_aug %*% jac_full)
          tmp_gammant_tmpt <- lav_mat_k_gammant_kt(m_k = tmp_aug, s = cov_aug,
            meanstructure = lavmodel@meanstructure, x_idx = x_idx)
        } else {
          tmp_gammant_tmpt <- lav_mat_k_gammant_kt(m_k = tmp, s = cov_g,
            meanstructure = lavmodel@meanstructure, x_idx = x_idx)
        }
        #tmp_gammant_tmpt_bis <- tmp %*% gamma_big %*% t(tmp)
        vcov[free_undirected_idx, free_undirected_idx] <-
          tmp_gammant_tmpt / lavsamplestats@ntotal
      }

      # iv_vcov_stage2 == "delta"
    } else {
      # part a: effect of theta1
      if (!iv_vcov_jaca_numerical) {
        # use analytic expressions
        if (iv_varcov_method %in% c("ULS", "GLS")) {
          jac_a <- lav_sem_miiv_utils_jaca_uls_gls(
            lavmodel = lavmodel,
            lavpartable = lavpartable, lavh1 = lavh1,
            free_directed_idx = free_directed_noint_idx,
            free_undirected_idx = free_undirected_idx,
            iv_varcov_method = iv_varcov_method
          )
        } else if (iv_varcov_method == "2RLS") {
          jac_a <- lav_sem_miiv_utils_jaca_2rls(
            lavmodel = lavmodel,
            lavpartable = lavpartable, lavh1 = lavh1,
            free_directed_idx = free_directed_noint_idx,
            free_undirected_idx = free_undirected_idx
          )
        } else if (iv_varcov_method == "RLS") {
          jac_a <- lav_sem_miiv_utils_jaca_rls(
            lavmodel = lavmodel,
            lavpartable = lavpartable, lavh1 = lavh1,
            free_directed_idx = free_directed_noint_idx,
            free_undirected_idx = free_undirected_idx
          )
        }
        # for checking only: numerical approximation
      } else {
        jac_a <- lav_func_jacobian_complex(
          func = lav_sem_miiv_varcov,
          x = theta1_noint, samplestats = FALSE, lavmodel = lavmodel,
          lavpartable = lavpartable, lavsamplestats = lavsamplestats,
          lavh1 = lavh1, free_directed_idx = free_directed_noint_idx,
          free_undirected_idx = free_undirected_idx,
          iv_varcov_method = iv_varcov_method
        )

        # jac_a <- lav_func_jacobian_complex(
        #   func = lav_sem_miiv_varcov_old,
        #   x = theta1_noint, samplestats = FALSE, lavmodel = lavmodel,
        #   lavpartable = lavpartable, lavsamplestats = lavsamplestats,
        #   lavh1 = lavh1, free_directed_idx = free_directed_noint_idx,
        #   free_undirected_idx = free_undirected_idx,
        #   iv_varcov_method = iv_varcov_method
        # )
      }

      # compute vcov_a (due to theta1)
      vcov_directed <-
        vcov[free_directed_noint_idx, free_directed_noint_idx, drop = FALSE]
      vcov_a <- jac_a %*% vcov_directed %*% t(jac_a)

      # part b: effect of sample statistics
      # - get jacobian lav_sem_miiv_varcov wrt samplestats_vec
      # - get Gamma (NT or ADF)
      # - vcov_undirected_b <- jac_b %*% (1/N * gamma_mat) %*% t(jac_b)

      if (!iv_vcov_jaca_numerical) {
        delta_block <- lav_sem_miiv_delta(lavmodel = lavmodel)
        jac_b_block <- vector("list", length = nblocks)
        for (b in seq_len(nblocks)) {
          sample_cov <- lavh1$implied$cov[[b]]
          this_delta2 <- delta_block[[b]][, free_undirected_idx, drop = FALSE]
          if (iv_varcov_method == "ULS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_uls(sample_cov, delta2 = this_delta2)
          } else if (iv_varcov_method == "GLS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_gls(sample_cov, delta2 = this_delta2,
                                          theta2 = theta2)
          } else if (iv_varcov_method == "2RLS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_2rls(sample_cov, delta2 = this_delta2,
                                           theta2 = theta2)
          } else if (iv_varcov_method == "RLS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_rls(sample_cov, delta2 = this_delta2,
                                          theta2 = theta2)
          }
        } # blocks
        jac_b <- lav_mat_bdiag(jac_b_block)

      # for checking only: numerical approximation (very slow!)
      } else {
        vec <- lav_implied_to_vec(
          implied = lavh1$implied, lavmodel = lavmodel,
          drop_list = TRUE
        )

        jac_b <- lav_func_jacobian_complex(
          func = lav_sem_miiv_varcov,
          x = vec, samplestats = TRUE, lavmodel = lavmodel,
          lavpartable = lavpartable, lavsamplestats = lavsamplestats,
          lavh1 = lavh1, free_directed_idx = free_directed_noint_idx,
          free_undirected_idx = free_undirected_idx,
          iv_varcov_method = iv_varcov_method
        )

        # jac_b <- lav_func_jacobian_complex(
        #   func = lav_sem_miiv_varcov_old,
        #   x = vec, samplestats = TRUE, lavmodel = lavmodel,
        #   lavpartable = lavpartable, lavsamplestats = lavsamplestats,
        #   lavh1 = lavh1, free_directed_idx = free_directed_noint_idx,
        #   free_undirected_idx = free_undirected_idx,
        #   iv_varcov_method = iv_varcov_method
        # )
      }
      if (use_explicit_gamma) {
        jacb_gamma_jacbt <- jac_b %*% gamma_big %*% t(jac_b)
      } else {
        x_idx <-
          if (lavmodel@fixed.x) lavsamplestats@x.idx[[1]] else integer(0L)
        jacb_gamma_jacbt <- lav_mat_k_gammant_kt(m_k = jac_b, s = cov_g,
          meanstructure = lavmodel@meanstructure, x_idx = x_idx)
      }
      vcov_b <- jacb_gamma_jacbt / lavsamplestats@ntotal
      vcov_ab <- vcov_a + vcov_b
      vcov[free_undirected_idx, free_undirected_idx] <- vcov_ab
    } # continuous
  } # iv_vcov_stage2 != "none"

  # the rest of lavaan expects the vcov of a ceq.simple.only model to be in the
  # 'unco' space (one row/column per non-collapsed parameter, nx.unco), eg in
  # lav_model_vcov_se(). The vcov above is built in the compact (nx.free)
  # space, so expand it via the ceq.simple.K mapping: constrained parameters
  # then share their (co)variance, and the SEs map to the parameter table
  # correctly (without this, the constrained parameters get inconsistent SEs
  # and the trailing parameters get NA).
  if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L) {
    vcov <- lavmodel@ceq.simple.K %*% vcov %*% t(lavmodel@ceq.simple.K)
  }

  vcov
}
