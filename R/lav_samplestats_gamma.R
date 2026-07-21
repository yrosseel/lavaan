# Gamma: N times the asymptotic variance matrix of the sample statistics
# (means/thresholds/slopes first, then vech of the (co)variances).
#
# Two computational routes exist, tied together by the identity
#
#     Gamma_ADF = A1^{-1} %*% J1 %*% A1^{-1}
#
# with A1 the expected h1 information (A1^{-1} = Gamma_NT) and J1 the
# first-order h1 information (crossprod of the casewise h1 scores / n),
# both at the saturated sample moments:
#
# - the DIRECT route (lav_samp_gamma below): crossprod of the casewise
#   moment contributions. Since sc_i = A1 %*% zc_i, this equals the
#   sandwich exactly (machine precision) and is 1.5-3x faster (the
#   sandwich adds two pstar x pstar multiplications); it is the canonical
#   engine for complete single-level data.
# - the SANDWICH route: needed whenever the casewise contributions are not
#   iid draws of a simple moment vector -- sampling weights (the wt branch
#   of lav_samp_gamma below), missing data (lav_mvn_mi_h1_omega_sw), and
#   two-level data (lav_samp_gamma_2l_g), where it is the only route.
#
# (evidence + benchmarks: aja_identity.R in the gamma-refactor test dir)
#
# SCALING CONVENTIONS
# - @SampleStats@NACOV[[g]] stores Gamma_g = n_g * acov(s_g)
#   = Cov(sqrt(n_g) * s_g), per group. (For two-level data n_g is the
#   number of OBSERVATIONS, even though the clusters are the independent
#   units: Gamma_g = n_g * acov always holds.)
# - the h1 information matrices (lav_model_h1_information.R) are UNIT
#   (per-observation single-level, per-cluster two-level); acov versions
#   (lav_model_h1_acov) carry the 1/n_g; the categorical WLS.W is stored
#   as Gamma/n and rescaled at consumption.
# - multigroup consumers stack per-group blocks into one 'global' Gamma.
#   Two equivalent conventions are in use:
#     (a) Gamma_g / fg  (= Cov(sqrt(ntotal) * s_g)) with fg-weighted
#         information/weight blocks appearing TWICE in U
#         (satorra.bentler, yuan.bentler, catml_dwls fit indices);
#     (b) fg * Gamma_g with UNweighted U blocks (lav_test_diff).
#   These give identical tr(U Gamma) and tr((U Gamma)^2): every (g,h)
#   cross-term picks up exactly one fg and one fh under either placement
#   (verified numerically in fg_convention_check.R, gamma-refactor dir).

# rescale a per-group Gamma list from the Cov(sqrt(n_g) s_g) convention
# (as stored in @NACOV) to the Cov(sqrt(ntotal) s_g) convention
# (Gamma_g / fg); convention (a) above, used when stacking the per-group
# blocks into one 'global' Gamma
lav_gamma_rescale_ntotal <- function(gamma, nobs, ntotal = NULL) {
  nobs <- unlist(nobs)
  if (is.null(ntotal)) {
    ntotal <- sum(nobs)
  }
  fg <- nobs / ntotal
  for (g in seq_along(gamma)) {
    gamma[[g]] <- gamma[[g]] / fg[g]
  }
  gamma
}

# 'the' Gamma (NACOV) of a fitted object: the fit-time stored NACOV when
# available, otherwise recomputed with the fit-time settings (the recipe).
# The single entry point for consumers that need "the Gamma that was (or
# would have been) used by this object" (satorra.bentler/fmg tests, sam,
# lavInspect "gamma", ...).
lav_gamma_used <- function(lavobject = NULL,
                           # or individual slots
                           lavdata = NULL,
                           lavoptions = NULL,
                           lavsamplestats = NULL,
                           lavh1 = NULL,
                           lavimplied = NULL) {
  if (!is.null(lavobject)) {
    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
    lavsamplestats <- lavobject@SampleStats
    lavh1 <- lavobject@h1
    lavimplied <- lavobject@implied
  }

  # fit-time NACOV?
  if (length(lavsamplestats@NACOV) > 0L &&
      !is.null(lavsamplestats@NACOV[[1]])) {
    return(lavsamplestats@NACOV)
  }

  # recompute with the fit-time settings
  lav_object_gamma(
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavsamplestats = lavsamplestats,
    lavh1 = lavh1,
    lavimplied = lavimplied,
    adf = TRUE,
    model_based = FALSE
  )
}

# the ADF 'meat': crossproduct of the (centered) casewise moment
# contributions zc, divided by n; when cluster_idx is given, the
# contributions are aggregated within clusters first (the G/(G-1)
# finite-sample correction is applied by the caller). Shared by the three
# (plain / fixed.x / conditional.x) branches of lav_samp_gamma().
lav_samp_gamma_meat <- function(zc, cluster_idx = NULL, n = nrow(zc)) {
  if (length(cluster_idx) > 0L) {
    zc <- rowsum(zc, cluster_idx)
  }
  if (anyNA(zc)) {
    lav_mat_crossprod(zc) / n
  } else {
    base::crossprod(zc) / n
  }
}

# weighted crossproduct of casewise scores (the 'meat' of a sandwich) under
# sampling weights; "design" weights the crossproduct by wt^2, "frequency"
# by wt (as if each case were physically replicated wt times). Shared by
# the ADF Gamma and the (PML) first-order h1 information.
lav_samp_wt_crossprod <- function(sc, wt, sampling_weights_type = "design") {
  if (identical(sampling_weights_type, "frequency")) {
    crossprod(sqrt(wt) * sc)
  } else {
    crossprod(wt * sc)
  }
}

# DLS with a model-based GammaNT: the (inverted) weight matrix for one
# group, W_DLS^{-1}, where
#     W_DLS = (1 - a) * Gamma_ADF (= NACOV) + a * Gamma_NT(m_cov, m_mean)
# shared by the objective, gradient, h1-information and object-generate code
lav_dls_wls_v_g <- function(m_cov = NULL, m_mean = NULL,
                            nacov_g = NULL, dls_a = 1.0,
                            x_idx = integer(0L), fixed_x = FALSE,
                            conditional_x = FALSE, meanstructure = FALSE) {
  gamma_nt <- lav_samp_gamma_nt(
    m_cov          = m_cov,
    m_mean         = m_mean,
    x_idx          = x_idx,
    fixed_x        = fixed_x,
    conditional_x  = conditional_x,
    meanstructure  = meanstructure,
    slopestructure = conditional_x
  )
  w_dls <- (1 - dls_a) * nacov_g + dls_a * gamma_nt
  lav_mat_sym_inverse(w_dls)
}

# for internal use -- lavobject or internal slots
lav_object_gamma <- function(lavobject = NULL,
                             # or individual slots
                             lavdata = NULL,
                             lavoptions = NULL,
                             lavsamplestats = NULL,
                             lavh1 = NULL,
                             lavimplied = NULL,
                             # other options
                             adf = TRUE, model_based = FALSE,
                             mplus_wls = NULL) {
  # extract slots
  if (!is.null(lavobject)) {
    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
    lavsamplestats <- lavobject@SampleStats
    lavh1 <- lavobject@h1
    lavimplied <- lavobject@implied
  }

  # the Gamma recipe: the same settings as used at fit time
  recipe <- lav_gamma_recipe(lavoptions = lavoptions, lavdata = lavdata,
                             nacov_compute = TRUE)

  missing <- recipe$missing
  fixed_x <- recipe$fixed.x
  conditional_x <- recipe$conditional.x
  meanstructure <- recipe$meanstructure
  gamma_n_minus_one <- recipe$n.minus.one
  gamma_unbiased <- recipe$unbiased
  if (is.null(mplus_wls)) {
    # gamma.wls.mplus, but only where fit time applies it (plain ADF Gamma
    # for the (D)WLS-family estimators)
    mplus_wls <- recipe$mplus.wls
  }

  # output container
  out <- vector("list", length = lavdata@ngroups)

  # ---- two-level data: cluster sandwich at the h1 estimates --------------
  if (recipe$multilevel) {
    if (recipe$categorical) {
      lav_msg_stop(gettext(
        "Gamma can not be (re)computed for two-level categorical data; the
         (stage-wise) Gamma used at fit time (if any) is stored in the
         NACOV slot."))
    }
    if (model_based || !adf) {
      lav_msg_stop(gettext(
        "only the sandwich-type Gamma (at the saturated h1 estimates) is
         available for two-level data."))
    }
    if (is.null(lavh1) || length(lavh1$implied) == 0L) {
      lav_msg_stop(gettext(
        "the h1 (saturated) estimates are needed to compute the two-level
         Gamma, but they are not available in this object."))
    }
    for (g in seq_len(lavdata@ngroups)) {
      gamma_g <- lav_samp_gamma_2l_g(
        lavsamplestats = lavsamplestats,
        lavh1 = lavh1,
        lavdata = lavdata,
        g = g
      )
      attr(gamma_g, "keep") <- NULL
      out[[g]] <- gamma_g
    }
    return(out)
  }

  # ---- categorical data (single-level): muthen1984 three-stage -----------
  if (recipe$categorical) {
    if (!adf || model_based) {
      lav_msg_stop(gettext(
        "only the sample-based (ADF-type) Gamma is available for
         categorical data."))
    }
    data_ov <- lavdata@ov
    for (g in seq_len(lavdata@ngroups)) {
      ov_names_g <- lavdata@ov.names[[g]]
      ov_types_g <- data_ov$type[match(ov_names_g, data_ov$name)]
      ov_levels_g <- data_ov$nlev[match(ov_names_g, data_ov$name)]
      if (length(lavdata@cluster) > 0L) {
        cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster_idx <- NULL
      }
      cat_1 <- muthen1984(
        data_1 = lavdata@X[[g]],
        wt = lavdata@weights[[g]],
        sampling_weights_type = recipe$swt.type,
        ov_names = ov_names_g,
        ov_types = ov_types_g,
        ov_levels = ov_levels_g,
        ov_names_x = if (conditional_x) lavdata@ov.names.x[[g]] else NULL,
        exo = if (conditional_x) lavdata@eXo[[g]] else NULL,
        group = g, # for error messages only
        wls_w = TRUE,
        zero_add = lavoptions$zero.add,
        zero_keep_margins = lavoptions$zero.keep.margins,
        zero_cell_warn = FALSE,
        zero_cell_tables = TRUE,
        allow_empty_cell = lavoptions$allow.empty.cell,
        cluster_idx = cluster_idx
      )
      nobs_g <- lavsamplestats@nobs[[g]]
      if (!is.null(cat_1$WLS.W.cluster)) {
        # cluster-robust meat + G/(G-1) correction, as at fit time
        nc <- lavdata@Lp[[g]]$nclusters[[2]]
        gamma_g <- cat_1$WLS.W.cluster * nobs_g * (nc / (nc - 1))
      } else {
        gamma_g <- cat_1$WLS.W * nobs_g
        if (gamma_n_minus_one) {
          gamma_g <- gamma_g * (nobs_g / (nobs_g - 1L))
        }
      }
      if (lavoptions$estimator == "catML") {
        # remove all but the correlation part, as at fit time
        ntotal <- nrow(gamma_g)
        pstar <- nrow(cat_1$A22)
        nocor <- ntotal - pstar
        if (length(nocor) > 0L) {
          gamma_g <- gamma_g[-seq_len(nocor), -seq_len(nocor)]
        }
      }
      if (recipe$group.w.free) {
        gamma_g <- lav_mat_bdiag(matrix(1, 1, 1), gamma_g)
      }
      out[[g]] <- gamma_g
    }
    return(out)
  }

  # ---- incomplete data (ml / two.stage / robust.two.stage) ---------------
  # Gamma = N x acov of the saturated estimates under missingness, from the
  # missing-data h1 information at the EM estimates (or at the model-implied
  # moments if model_based): the sandwich I^{-1} J I^{-1} (adf) or the
  # inverted observed information I^{-1} (not adf). This mirrors both the
  # fit-time NACOV of the two.stage least-squares estimators and the
  # two-stage robust vcov machinery.
  if (!missing %in% c("listwise", "pairwise")) {
    sandwich <- adf
    if (missing == "two.stage") {
      sandwich <- FALSE
    } else if (missing == "robust.two.stage") {
      sandwich <- TRUE
    }
    for (g in seq_len(lavdata@ngroups)) {
      if (conditional_x) {
        # conditional metric: Gamma is the N x acov of the saturated
        # (vec(Beta), vech(res.cov)) estimates (x is complete)
        if (model_based) {
          res_int_g <- lavimplied$res.int[[g]]
          res_slopes_g <- lavimplied$res.slopes[[g]]
          res_cov_g <- lavimplied$res.cov[[g]]
        } else {
          res_int_g <- lavsamplestats@missing.h1[[g]]$res.int
          res_slopes_g <- lavsamplestats@missing.h1[[g]]$res.slopes
          res_cov_g <- lavsamplestats@missing.h1[[g]]$res.cov
        }
        if (is.null(res_cov_g)) {
          lav_msg_stop(gettext(
            "Gamma is not available: no EM estimates of the saturated
            conditional moments (zero coverage?)."))
        }
        if (sandwich) {
          gamma_g <- lav_mvreg_mi_h1_omega_sw(
            y = lavdata@X[[g]], exo = lavdata@eXo[[g]],
            mp = lavdata@Mp[[g]],
            yp = lavsamplestats@missing[[g]], wt = lavdata@weights[[g]],
            res_int = res_int_g, res_slopes = res_slopes_g,
            res_cov = res_cov_g,
            information = "observed"
          )
        } else {
          i1 <- lav_mvreg_mi_information_observed_samplestats(
            yp = lavsamplestats@missing[[g]],
            res_int = res_int_g, res_slopes = res_slopes_g,
            res_cov = res_cov_g
          )
          gamma_g <- lav_mat_sym_inverse(i1)
        }
        if (recipe$group.w.free) {
          gamma_g <- lav_mat_bdiag(matrix(1, 1, 1), gamma_g)
        }
        out[[g]] <- gamma_g
        next
      }
      if (model_based) {
        mu_g <- lavimplied$mean[[g]]
        sigma_g <- lavimplied$cov[[g]]
      } else {
        mu_g <- lavsamplestats@missing.h1[[g]]$mu
        sigma_g <- lavsamplestats@missing.h1[[g]]$sigma
      }
      x_idx_g <- if (fixed_x) lavsamplestats@x.idx[[g]] else integer(0L)
      if (sandwich) {
        gamma_g <- lav_mvn_mi_h1_omega_sw(
          y = lavdata@X[[g]], mp = lavdata@Mp[[g]],
          yp = lavsamplestats@missing[[g]], wt = lavdata@weights[[g]],
          mu = mu_g, sigma_1 = sigma_g, x_idx = x_idx_g,
          information = "observed"
        )
      } else {
        i1 <- lav_mvnorm_missing_information_observed_samplestats(
          yp = lavsamplestats@missing[[g]],
          mu = mu_g, sigma_1 = sigma_g, x_idx = x_idx_g
        )
        gamma_g <- lav_mat_sym_inverse(i1)
      }
      if (recipe$group.w.free) {
        gamma_g <- lav_mat_bdiag(matrix(1, 1, 1), gamma_g)
      }
      out[[g]] <- gamma_g
    }
    return(out)
  }

  # ---- continuous, complete-data, single-level ----------------------------
  if (adf && lavdata@data.type != "full") {
    lav_msg_stop(gettext(
      "the (ADF) Gamma matrix can not be computed without full data; please
       provide the NACOV= argument."))
  }
  if (adf && model_based && conditional_x) {
    lav_msg_stop(gettext(
      "ADF + model.based + conditional.x is not supported yet."))
  }
  if (adf && model_based && recipe$correlation) {
    lav_msg_stop(gettext(
      "ADF + model.based is not supported yet for correlation structures."))
  }

  # compute Gamma matrix for each group
  for (g in seq_len(lavdata@ngroups)) {
    x_idx <- lavsamplestats@x.idx[[g]]
    cov_1 <- mean_1 <- NULL
    if (!adf || model_based) {
      implied <- lavh1$implied # saturated/unstructured
      if (model_based) {
        implied <- lavimplied # model-based/structured
      }
      if (conditional_x) {
        # convert to joint COV/MEAN
        res_s <- implied$res.cov[[g]]
        res_slopes <- implied$res.slopes[[g]]
        res_int <- implied$res.int[[g]]
        s_xx <- implied$cov.x[[g]]
        m_x <- implied$mean.x[[g]]

        s_yy <- res_s + res_slopes %*% s_xx %*% t(res_slopes)
        s_yx <- res_slopes %*% s_xx
        s_xy <- s_xx %*% t(res_slopes)
        m_y <- res_int + res_slopes %*% m_x

        cov_1 <- rbind(cbind(s_yy, s_yx), cbind(s_xy, s_xx))
        mean_1 <- c(m_y, m_x)
      } else {
        # not conditional.x
        cov_1 <- implied$cov[[g]]
        mean_1 <- implied$mean[[g]]
      }
    } # COV/MEAN

    if (adf) {
      if (conditional_x) {
        y <- cbind(lavdata@X[[g]], lavdata@eXo[[g]])
      } else {
        y <- lavdata@X[[g]]
      }
      if (length(lavdata@cluster) > 0L) {
        cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster_idx <- NULL
      }
      if (recipe$correlation) {
        # (partial) correlation ADF Gamma via the delta method, as at fit
        # time (no cluster/weights support in the leaf function)
        cor_idx_g <- lav_gamma_recipe_cor_idx(recipe, lavdata@ov.names[[g]])
        out[[g]] <- lav_samp_partial_cor_gamma(
          m_y = y, cor_idx = cor_idx_g,
          meanstructure = meanstructure,
          fixed_x = fixed_x, x_idx = x_idx,
          conditional_x = conditional_x
        )
      } else {
        out[[g]] <- lav_samp_gamma(
          m_y = y,
          m_mu = mean_1,
          m_sigma = cov_1,
          x_idx = x_idx,
          cluster_idx = cluster_idx,
          fixed_x = fixed_x,
          conditional_x = conditional_x,
          meanstructure = meanstructure,
          slopestructure = conditional_x,
          gamma_n_minus_one = gamma_n_minus_one,
          unbiased = gamma_unbiased,
          mplus_wls = mplus_wls,
          wt = lavdata@weights[[g]],
          sampling_weights_type = recipe$swt.type
        )
      }
    } else {
      if (recipe$correlation) {
        cor_idx_g <- lav_gamma_recipe_cor_idx(recipe, lavdata@ov.names[[g]])
        out[[g]] <- lav_samp_partial_cor_gamma_nt(
          m_cov = cov_1, cor_idx = cor_idx_g,
          meanstructure = meanstructure,
          fixed_x = fixed_x, x_idx = x_idx,
          conditional_x = conditional_x,
          m_mean = mean_1
        )
      } else {
        out[[g]] <- lav_samp_gamma_nt(
          m_cov = cov_1, # joint!
          m_mean = mean_1, # joint!
          x_idx = x_idx,
          fixed_x = fixed_x,
          conditional_x = conditional_x,
          meanstructure = meanstructure,
          slopestructure = conditional_x
        )
      }
    }

    # group.w.free: prepend the group-weight element; the fit-time weight
    # a = group_w * ntotal / nobs is identically 1, so this matches the
    # NACOV stored at fit time
    if (recipe$group.w.free) {
      out[[g]] <- lav_mat_bdiag(matrix(1, 1, 1), out[[g]])
    }

  } # g

  out
}





# NOTE:
#  - three types:
#       1) plain         (conditional.x = FALSE, fixed.x = FALSE)
#       2) fixed.x       (conditional.x = FALSE, fixed.x = TRUE)
#       3) conditional.x (conditional.x = TRUE)
#  - if conditional.x = TRUE, we ignore fixed.x (can be TRUE or FALSE)

# the E[(1,x)(1,x)'] moment matrix of the (intercept, exo) regressors; the
# building block of the mean/slope block of the conditional.x NT Gamma
# (lav_samp_gamma_nt) and of its direct inverse (lav_samp_gamma_inverse_nt)
lav_samp_gamma_c3 <- function(m_mx, m_c) {
  rbind(
    c(1, m_mx),
    cbind(m_mx, m_c + tcrossprod(m_mx))
  )
}

# NORMAL-THEORY
lav_samp_gamma_nt <- function(m_y = NULL, # should include
                                     # eXo if
                                     # conditional.x=TRUE
                                     wt = NULL,
                                     m_cov = NULL, # joint!
                                     m_mean = NULL, # joint!
                                     rescale = FALSE,
                                     x_idx = integer(0L),
                                     fixed_x = FALSE,
                                     conditional_x = FALSE,
                                     meanstructure = FALSE,
                                     slopestructure = FALSE) {
  # check arguments
  if (length(x_idx) == 0L) {
    conditional_x <- FALSE
    fixed_x <- FALSE
  }

  # compute m_cov from m_y
  if (is.null(m_cov)) {
    stopifnot(!is.null(m_y))

    # coerce to matrix
    m_y <- unname(as.matrix(m_y))
    n <- nrow(m_y)
    if (is.null(wt)) {
      m_cov <- cov(m_y)
      if (rescale) {
        m_cov <- m_cov * (n - 1) / n # (normal) ML version
      }
    } else {
      out <- stats::cov.wt(m_y, wt = wt, method = "ML")
      m_cov <- out$cov
    }
  } else {
    if (!missing(rescale)) {
      lav_msg_warn(gettext("rescale= argument has no effect if m_cov is given"))
    }
    if (!missing(wt)) {
      lav_msg_warn(gettext("wt= argument has no effect if m_cov is given"))
    }
  }

  # if needed, compute m_mean from m_y
  if (conditional_x && length(x_idx) > 0L && is.null(m_mean) &&
    (meanstructure || slopestructure)) {
    stopifnot(!is.null(m_y))
    if (is.null(wt)) {
      m_mean <- colMeans(m_y, na.rm = TRUE)
    } else {
      m_mean <- out$center
    }
  }

  # rename
  m_s <- m_cov
  m_m <- m_mean

  # unconditional
  if (!conditional_x) {
    # unconditional - stochastic x
    if (!fixed_x) {
      # if (lav_use_lavaanC()) {
      #   m_gamma <- lavaanC::m_kronecker_dup_ginv_pre_post(m_s,
      #                                              multiplicator = 2.0)
      # } else {
        m_gamma <- 2 * lav_mat_dup_ginv_pre_post(m_s %x% m_s)
      # }
      if (meanstructure) {
        m_gamma <- lav_mat_bdiag(m_s, m_gamma)
      }

      # unconditional - fixed x
    } else {
      # handle fixed_x = TRUE

      # cov(m_y|X) = m_a - m_b m_c^{-1} m_b'
      # where m_a = cov(m_y), m_b = cov(m_y,X), m_c = cov(X)
      m_a <- m_s[-x_idx, -x_idx, drop = FALSE]
      m_b <- m_s[-x_idx, x_idx, drop = FALSE]
      m_c <- m_s[x_idx, x_idx, drop = FALSE]
      ybarx <- m_a - m_b %*% solve(m_c, t(m_b))

      # reinsert ybarx in m_y+X (residual) covariance matrix
      ybarx_aug <- matrix(0, nrow = NROW(m_s), ncol = NCOL(m_s))
      ybarx_aug[-x_idx, -x_idx] <- ybarx

      # take difference
      m_r <- m_s - ybarx_aug

      # if (lav_use_lavaanC()) {
      #   gamma_s <- lavaanC::m_kronecker_dup_ginv_pre_post(m_s,
      #                                               multiplicator = 2.0)
      #   gamma_r <- lavaanC::m_kronecker_dup_ginv_pre_post(m_r,
      #                                               multiplicator = 2.0)
      # } else {
        gamma_s <- 2 * lav_mat_dup_ginv_pre_post(m_s %x% m_s)
        gamma_r <- 2 * lav_mat_dup_ginv_pre_post(m_r %x% m_r)
      # }
      m_gamma <- gamma_s - gamma_r

      if (meanstructure) {
        m_gamma <- lav_mat_bdiag(ybarx_aug, m_gamma)
      }
    }
  } else {
    # conditional_x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    # regress m_y on X, and compute covariance of residuals 'm_r'
    m_a <- m_s[-x_idx, -x_idx, drop = FALSE]
    m_b <- m_s[-x_idx, x_idx, drop = FALSE]
    m_c <- m_s[x_idx, x_idx, drop = FALSE]
    cov_ybarx <- m_a - m_b %*% solve(m_c) %*% t(m_b)
    # if (lav_use_lavaanC()) {
    #   m_gamma <- lavaanC::m_kronecker_dup_ginv_pre_post(cov_ybarx,
    #                                                  multiplicator = 2.0)
    # } else {
      m_gamma <- 2 *
               lav_mat_dup_ginv_pre_post(cov_ybarx %x% cov_ybarx)
    # }

    if (meanstructure || slopestructure) {
      m_mx <- m_m[x_idx]
      c3 <- lav_samp_gamma_c3(m_mx, m_c)
      # B3 <- cbind(m_my, m_b + tcrossprod(m_my,m_mx))
    }

    if (meanstructure) {
      if (slopestructure) {
        # a11 <- solve(c3) %x% cov_ybarx
        a11 <- cov_ybarx %x% solve(c3)
      } else {
        # a11 <- solve(c3)[1, 1, drop=FALSE] %x% cov_ybarx
        a11 <- cov_ybarx %x% solve(c3)[1, 1, drop = FALSE]
      }
    } else {
      if (slopestructure) {
        # a11 <- solve(c3)[-1, -1, drop=FALSE] %x% cov_ybarx
        a11 <- cov_ybarx %x% solve(c3)[-1, -1, drop = FALSE]
      } else {
        a11 <- matrix(0, 0, 0)
      }
    }

    m_gamma <- lav_mat_bdiag(a11, m_gamma)
  }

  m_gamma
}

# NOTE:
#  - three types:
#       1) plain         (conditional_x = FALSE, fixed_x = FALSE)
#       2) fixed_x       (conditional_x = FALSE, fixed_x = TRUE)
#       3) conditional_x (conditional_x = TRUE)
#  - if conditional_x = TRUE, we ignore fixed_x (can be TRUE or FALSE)
#
#  - new in 0.6-1: if Mu/Sigma is provided, compute 'model-based' Gamma
#                  (only if conditional_x = FALSE, for now)
#  - new in 0.6-2: if cluster.idx is not NULL, correct for clustering
#  - new in 0.6-13: add unbiased = TRUE (for the 'plain' setting only)

# ADF THEORY
lav_samp_gamma <- function(m_y, # Y+X if cond!
                                  m_mu = NULL,
                                  m_sigma = NULL,
                                  x_idx = integer(0L),
                                  cluster_idx = NULL,
                                  fixed_x = FALSE,
                                  conditional_x = FALSE,
                                  meanstructure = FALSE,
                                  slopestructure = FALSE,
                                  gamma_n_minus_one = FALSE,
                                  unbiased = FALSE,
                                  mplus_wls = FALSE,
                                  wt = NULL,
                                  sampling_weights_type = "design") {
  # coerce to matrix
  m_y <- unname(as.matrix(m_y))
  n <- nrow(m_y)
  p <- ncol(m_y)

  # sampling weights?
  # Gamma is the asymptotic (co)variance of the sample moments, which for
  # the saturated (h1) model equals the sandwich I^{-1} J I^{-1}, where I is
  # the expected and J the first-order h1 information matrix. We reuse the
  # (weight-aware) expected information I, and construct the meat J according
  # to how the sampling weights are interpreted (sampling.weights.type):
  #  - "design"    : J = crossprod(wt * sc) / sum(wt), i.e. sum(wt^2)-weighted,
  #                  the design-based sandwich; the parameter vcov is then
  #                  invariant to the overall scale of the weights, matching
  #                  the continuous ML robust (huber.white) SEs.
  #  - "frequency" : J = crossprod(sqrt(wt) * sc) / sum(wt), i.e. sum(wt)-
  #                  weighted; treats the weights as frequencies, so the result
  #                  equals the Gamma from physically replicating each case wt
  #                  times. (NB: this also feeds the WLS/DWLS/DLS weight matrix,
  #                  so it shifts those point estimates relative to "design".)
  # Under normality J = I and this reduces to the normal-theory Gamma.
  if (!is.null(wt)) {
    if (conditional_x) {
      lav_msg_stop(gettext(
        "Gamma with sampling weights is not available (yet) if
         conditional.x = TRUE."))
    }
    if (length(cluster_idx) > 0L) {
      lav_msg_stop(gettext(
        "Gamma with sampling weights is not available (yet) for
         clustered data."))
    }
    out <- stats::cov.wt(m_y, wt = wt, method = "ML")
    i_mat <- lav_mvn_h1_info_expected(
      y = m_y, wt = wt, sample_cov = out$cov, x_idx = x_idx,
      meanstructure = meanstructure
    )
    # unit (unweighted) casewise scores of the saturated model, evaluated
    # at the weighted sample moments
    sc <- lav_mvn_sc_mu_sigma(
      y = m_y, mu = out$center, sigma_1 = out$cov, x_idx = x_idx
    )
    if (!meanstructure) {
      sc <- sc[, -seq_len(p), drop = FALSE]
    }
    j_mat <- lav_samp_wt_crossprod(sc, wt, sampling_weights_type) / sum(wt)
    # fixed.x zeroes the x-blocks of I, so use a pseudo-inverse there
    if (length(x_idx) > 0L) {
      i_inv <- MASS::ginv(i_mat)
    } else {
      i_inv <- solve(i_mat)
    }
    return(i_inv %*% j_mat %*% i_inv)
  }

  # unbiased?
  if (unbiased) {
    if (conditional_x || fixed_x || !is.null(m_sigma) ||
      !is.null(cluster_idx)) {
      lav_msg_stop(
        gettext("unbiased Gamma only available for the simple
        (not conditional_x or fixed_x or model-based or clustered) setting."))
    } else {
    # data really should be complete
      m_cov <- stats::cov(m_y, use = "pairwise.complete.obs")
      m_cov <- m_cov * (n - 1) / n
      cov_vech <- lav_mat_vech(m_cov)
    }
  }

  # model-based?
  if (!is.null(m_sigma)) {
    stopifnot(!conditional_x)
    model_based <- TRUE
    if (meanstructure) {
      stopifnot(!is.null(m_mu))
      sigma <- c(as.numeric(m_mu), lav_mat_vech(m_sigma))
    } else {
      m_mu <- colMeans(m_y, na.rm = TRUE) # for centering!
      sigma <- lav_mat_vech(m_sigma)
    }
  } else {
    model_based <- FALSE
  }

  # denominator
  if (gamma_n_minus_one) {
    n <- n - 1
  }

  # check arguments
  if (length(x_idx) == 0L) {
    conditional_x <- FALSE
    fixed_x <- FALSE
  }
  if (mplus_wls) {
    stopifnot(!conditional_x, !fixed_x)
  }

  if (!conditional_x && !fixed_x) {
    # center, so we can use crossprod instead of cov
    if (model_based) {
      m_yc <- t(t(m_y) - as.numeric(m_mu))
    } else {
      m_yc <- t(t(m_y) - colMeans(m_y, na.rm = TRUE))
    }

    # create z where the rows_i contain the following elements:
    #  - Y_i (if meanstructure is TRUE)
    #  - vech(Yc_i' %*% Yc_i) where Yc_i are the residuals
    idx1 <- lav_mat_vech_col_idx(p)
    idx2 <- lav_mat_vech_row_idx(p)
    if (meanstructure) {
      z <- cbind(m_y, m_yc[, idx1, drop = FALSE] *
        m_yc[, idx2, drop = FALSE])
    } else {
      z <- (m_yc[, idx1, drop = FALSE] *
        m_yc[, idx2, drop = FALSE])
    }

    if (model_based) {
      if (meanstructure) {
        stopifnot(!is.null(m_mu))
        sigma <- c(as.numeric(m_mu), lav_mat_vech(m_sigma))
      } else {
        sigma <- lav_mat_vech(m_sigma)
      }
      zc <- t(t(z) - sigma)
    } else {
      zc <- t(t(z) - colMeans(z, na.rm = TRUE))
    }

    # meat (aggregated within clusters first, if clustered)
    m_gamma <- lav_samp_gamma_meat(zc, cluster_idx = cluster_idx, n = n)
  } else if (!conditional_x && fixed_x) {
    if (model_based) {
      m_yc <- t(t(m_y) - as.numeric(m_mu))
      y_bar <- colMeans(m_y, na.rm = TRUE)
      res_cov <- (m_sigma[-x_idx, -x_idx, drop = FALSE] -
        m_sigma[-x_idx, x_idx, drop = FALSE] %*%
        solve(m_sigma[x_idx, x_idx, drop = FALSE]) %*%
        m_sigma[x_idx, -x_idx, drop = FALSE])
      res_slopes <- (solve(m_sigma[x_idx, x_idx, drop = FALSE]) %*%
        m_sigma[x_idx, -x_idx, drop = FALSE])
      res_int <- (y_bar[-x_idx] -
        as.numeric(colMeans(m_y[, x_idx, drop = FALSE],
          na.rm = TRUE
        ) %*% res_slopes))
      x_bar <- y_bar[x_idx]
      yhat__bar <- as.numeric(res_int + as.numeric(x_bar) %*% res_slopes)
      yhat_bar <- numeric(p)
      yhat_bar[-x_idx] <- yhat__bar
      yhat_bar[x_idx] <- x_bar
      yhat_cov <- m_sigma
      yhat_cov[-x_idx, -x_idx] <- m_sigma[-x_idx, -x_idx] - res_cov


      yhat <- cbind(1, m_y[, x_idx]) %*% rbind(res_int, res_slopes)
      y_hat <- m_y
      y_hat[, -x_idx] <- yhat
      # y_hat <- cbind(yhat, m_y[,x_idx])
      y_hatc <- t(t(y_hat) - yhat_bar)
      idx1 <- lav_mat_vech_col_idx(p)
      idx2 <- lav_mat_vech_row_idx(p)
      if (meanstructure) {
        z <- (cbind(m_y, m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE]) -
          cbind(y_hat, y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE]))
        sigma1 <- c(m_mu, lav_mat_vech(m_sigma))
        sigma2 <- c(yhat_bar, lav_mat_vech(yhat_cov))
      } else {
        z <- (m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE] -
          y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE])
        sigma1 <- lav_mat_vech(m_sigma)
        sigma2 <- lav_mat_vech(yhat_cov)
      }
      zc <- t(t(z) - (sigma1 - sigma2))
    } else {
      m_qr <- qr(cbind(1, m_y[, x_idx, drop = FALSE]))
      yhat <- qr.fitted(m_qr, m_y[, -x_idx, drop = FALSE])
      # y_hat <- cbind(yhat, m_y[,x_idx])
      y_hat <- m_y
      y_hat[, -x_idx] <- yhat

      m_yc <- t(t(m_y) - colMeans(m_y, na.rm = TRUE))
      y_hatc <- t(t(y_hat) - colMeans(y_hat, na.rm = TRUE))
      idx1 <- lav_mat_vech_col_idx(p)
      idx2 <- lav_mat_vech_row_idx(p)
      if (meanstructure) {
        z <- (cbind(m_y, m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE]) -
          cbind(y_hat, y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE]))
      } else {
        z <- (m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE] -
          y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE])
      }
      zc <- t(t(z) - colMeans(z, na.rm = TRUE))
    }

    # meat (aggregated within clusters first, if clustered)
    m_gamma <- lav_samp_gamma_meat(zc, cluster_idx = cluster_idx, n = n)
  } else {
    # conditional_x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    # regress m_y on m_x, and compute residuals
    m_x <- cbind(1, m_y[, x_idx, drop = FALSE])
    m_qr <- qr(m_x)
    m_res <- qr.resid(m_qr, m_y[, -x_idx, drop = FALSE])
    p <- ncol(m_res)

    idx1 <- lav_mat_vech_col_idx(p)
    idx2 <- lav_mat_vech_row_idx(p)

    if (meanstructure || slopestructure) {
      xtx_inv <- unname(solve(crossprod(m_x)))
      xi <- (m_x %*% xtx_inv) * n ## FIXME, shorter way?
      nc_x <- NCOL(m_x)
      nc_y <- NCOL(m_res)
    }

    if (meanstructure) {
      if (slopestructure) {
        # xi_idx <- rep(seq_len(nc_x), each  = nc_y)
        # res_idx <- rep(seq_len(nc_y), times = nc_x)
        xi_idx <- rep(seq_len(nc_x), times = nc_y)
        res_idx <- rep(seq_len(nc_y), each = nc_x)
        z <- cbind(
          xi[, xi_idx, drop = FALSE] *
            m_res[, res_idx, drop = FALSE],
          m_res[, idx1, drop = FALSE] *
            m_res[, idx2, drop = FALSE]
        )
      } else {
        xi_idx <- rep(1L, each = nc_y)
        z <- cbind(
          xi[, xi_idx, drop = FALSE] *
            m_res,
          m_res[, idx1, drop = FALSE] *
            m_res[, idx2, drop = FALSE]
        )
      }
    } else {
      if (slopestructure) {
        # xi_idx <- rep(seq_len(nc_x), each  = nc_y)
        # xi_idx <- xi_idx[ -seq_len(nc_y) ]
        xi_idx <- rep(seq(2, nc_x), times = nc_y)
        # res_idx <- rep(seq_len(nc_y), times = (nc_x - 1L))
        res_idx <- rep(seq_len(nc_y), each = (nc_x - 1L))
        z <- cbind(
          xi[, xi_idx, drop = FALSE] *
            m_res[, res_idx, drop = FALSE],
          m_res[, idx1, drop = FALSE] *
            m_res[, idx2, drop = FALSE]
        )
      } else {
        z <- m_res[, idx1, drop = FALSE] * m_res[, idx2, drop = FALSE]
      }
    }

    if (model_based) {
      zc <- t(t(z) - sigma)
    } else {
      zc <- t(t(z) - colMeans(z, na.rm = TRUE))
    }

    # meat (aggregated within clusters first, if clustered)
    m_gamma <- lav_samp_gamma_meat(zc, cluster_idx = cluster_idx, n = n)
  }


  # only to mimic Mplus when estimator = "WLS"
  if (mplus_wls && !fixed_x && !conditional_x) {
    # adjust G_22 (the varcov part)
    m_s <- cov(m_y, use = "pairwise")
    w <- lav_mat_vech(m_s)
    w_biased <- (n - 1) / n * w
    diff <- outer(w, w) - outer(w_biased, w_biased)
    if (meanstructure) {
      m_gamma[-seq_len(p), -seq_len(p)] <-
        m_gamma[-seq_len(p), -seq_len(p), drop = FALSE] - diff
    } else {
      m_gamma <- m_gamma - diff
    }

    if (meanstructure) {
      # adjust G_12/G_21 (third-order)
      # strange rescaling?
      n1 <- (n - 1) / n
      m_gamma[seq_len(p), -seq_len(p)] <- m_gamma[seq_len(p), -seq_len(p)] * n1
      m_gamma[-seq_len(p), seq_len(p)] <- m_gamma[-seq_len(p), seq_len(p)] * n1
    }
  }

  # clustered? G/(G-1) finite-sample correction, G = number of clusters
  if (length(cluster_idx) > 0L) {
    nc <- length(unique(cluster_idx))
    m_gamma <- m_gamma * nc / (nc - 1)
  }

  # unbiased?
  if (unbiased) {
    # normal-theory m_gamma (cov only)
    # if (lav_use_lavaanC()) {
    #   gammant_cov <- lavaanC::m_kronecker_dup_ginv_pre_post(m_cov,
    #                                               multiplicator = 2.0)
    # } else {
      gammant_cov <- 2 * lav_mat_dup_ginv_pre_post(m_cov %x% m_cov)
    # }

    if (meanstructure) {
      gamma_cov <- m_gamma[-(1:p), -(1:p), drop = FALSE]
      gamma_mean_cov <- m_gamma[1:p, -(1:p), drop = FALSE]
    } else {
      gamma_cov <- m_gamma
    }

    # Browne's unbiased DF estimator (m_cov part)
    gamma_u <- (n * (n - 1) / (n - 2) / (n - 3) * gamma_cov -
      n / (n - 2) / (n - 3) * (gammant_cov -
        2 / (n - 1) * tcrossprod(cov_vech)))
    if (meanstructure) {
      m_gamma <- lav_mat_bdiag(m_cov, gamma_u)

      # 3-rd order:
      m_gamma[1:p, (p + 1):ncol(m_gamma)] <- gamma_mean_cov * n / (n - 2)
      m_gamma[(p + 1):ncol(m_gamma), 1:p] <- t(gamma_mean_cov * n / (n - 2))
    } else {
      m_gamma <- gamma_u
    }
  } # unbiased

  m_gamma
}

# ADF Gamma for correlations
#
# 30 May 2024: basic version: fixed_x=FALSE, conditional_x=FALSE, ...
# Jacobian of the (partial) correlation moment vector
#   [ mean?, var(num_idx), off-diagonal (co)variances ]
# with respect to the raw moments [ mean?, vech(S) ], evaluated at 'm_s'.
# 'cor_idx' lists the variables that are standardized to unit variance;
# num_idx (= complement) keep their original metric. This is the building
# block of the delta-method correlation Gamma (both NT and ADF), and works
# for full (cor_idx = all) as well as partial correlation structures.
lav_partial_cor_jacobian <- function(m_s, cor_idx = NULL,
                                      meanstructure = FALSE) {
  p <- nrow(m_s)
  is_cor <- logical(p)
  is_cor[cor_idx] <- TRUE
  num_idx <- which(!is_cor)
  # scaling: sd for correlation variables, 1 for the others
  d <- rep(1.0, p)
  d[is_cor] <- sqrt(diag(m_s)[is_cor])

  # vech bookkeeping (full, with diagonal)
  vr <- lav_mat_vech_row_idx(p)
  vc <- lav_mat_vech_col_idx(p)
  pstar <- length(vr)
  diag_pos <- which(vr == vc) # variances, in variable order
  off_pos <- which(vr != vc) # covariances
  var_pos <- diag_pos[num_idx] # variances we keep (non-correlation vars)

  # Jacobian (covariance part): rows = [var(num_idx), offdiag], cols = vech(S)
  nkeep <- length(var_pos) + length(off_pos)
  j_cov <- matrix(0, nkeep, pstar)
  if (length(var_pos) > 0L) {
    j_cov[cbind(seq_along(var_pos), var_pos)] <- 1
  }
  off0 <- length(var_pos)
  for (i in seq_along(off_pos)) {
    pos <- off_pos[i]
    k <- vr[pos]
    l <- vc[pos]
    m_kl <- m_s[k, l] / (d[k] * d[l])
    row <- off0 + i
    j_cov[row, pos] <- 1 / (d[k] * d[l])
    if (is_cor[k]) {
      j_cov[row, diag_pos[k]] <- j_cov[row, diag_pos[k]] -
        0.5 * m_kl / m_s[k, k]
    }
    if (is_cor[l]) {
      j_cov[row, diag_pos[l]] <- j_cov[row, diag_pos[l]] -
        0.5 * m_kl / m_s[l, l]
    }
  }

  if (meanstructure) {
    # partial moments: [mean(p), var(num_idx), offdiag]; the mean of a
    # standardized variable is scaled by 1/sd (first order), others by 1
    j_mat <- matrix(0, p + nkeep, p + pstar)
    j_mat[cbind(seq_len(p), seq_len(p))] <- 1 / d
    j_mat[(p + 1L):(p + nkeep), (p + 1L):(p + pstar)] <- j_cov
  } else {
    j_mat <- j_cov
  }

  j_mat
}

# conditional.x + correlation structure: the standardized CONDITIONAL moment
# vector
#   [ vecr(cbind(res.int, res.slopes))  (if meanstructure)
#     or vecr(res.slopes)               (slopes only),
#     vech(res.cov, diagonal = FALSE) ]
# as a smooth function of the JOINT raw moments [ mean (if meanstructure),
# vech(S) ]: first scale all variables to unit variance, then partition into
# (y, x) and compute the residual moments. This is the moment layout of
# lav_samp_wls_obs() for the continuous conditional.x + correlation case.
# Every operation is complex-step safe (no crossprod/tcrossprod: those
# conjugate complex input as of R 4.3).
lav_cor_cond_stats <- function(stats, p = NULL,
                               x_idx = integer(0L),
                               meanstructure = FALSE) {
  if (meanstructure) {
    m_m <- stats[seq_len(p)]
    s_vech <- stats[-seq_len(p)]
  } else {
    s_vech <- stats
  }
  m_s <- lav_mat_vech_rev(s_vech)

  # standardize all variables (full correlation structure)
  d <- sqrt(diag(m_s))
  m_r <- m_s / outer(d, d)

  # partition into (y, x) and compute the residual moments
  r_yx <- m_r[-x_idx, x_idx, drop = FALSE]
  r_xx <- m_r[x_idx, x_idx, drop = FALSE]
  pi_1 <- r_yx %*% solve(r_xx)
  res_cov <- m_r[-x_idx, -x_idx, drop = FALSE] - pi_1 %*% t(r_yx)

  if (meanstructure) {
    m_m <- m_m / d
    res_int <- m_m[-x_idx] - as.vector(pi_1 %*% m_m[x_idx])
    out1 <- lav_mat_vecr(cbind(res_int, pi_1))
  } else {
    out1 <- lav_mat_vecr(pi_1)
  }

  c(out1, lav_mat_vech(res_cov, diagonal = FALSE))
}

# Jacobian of the conditional correlation moment vector with respect to the
# joint raw moments [ mean (if meanstructure), vech(S) ], evaluated at
# (m_mean, m_cov); complex-step differentiation of the map above
lav_cor_cond_jacobian <- function(m_cov, m_mean = NULL,
                                  x_idx = integer(0L),
                                  meanstructure = FALSE) {
  p <- nrow(m_cov)
  if (meanstructure) {
    if (is.null(m_mean)) {
      m_mean <- numeric(p)
    }
    x0 <- c(m_mean, lav_mat_vech(m_cov))
  } else {
    x0 <- lav_mat_vech(m_cov)
  }
  lav_func_jacobian_complex(
    func = function(x) {
      lav_cor_cond_stats(x, p = p, x_idx = x_idx,
                         meanstructure = meanstructure)
    },
    x = x0, fallback_simple = FALSE
  )
}

# ADF (sandwich) Gamma for a (partial) correlation structure, via the delta
# method: Gamma_cor = J %*% Gamma_full %*% t(J), where Gamma_full is the raw
# ADF Gamma (lav_samp_gamma, which also handles fixed.x via x_idx) and J is
# the Jacobian above. If conditional_x = TRUE, m_y must be the joint (y, x)
# data and the result is the Gamma of the standardized conditional moments
# (x treated as fixed, as everywhere in the conditional.x machinery).
lav_samp_partial_cor_gamma <- function(m_y, cor_idx = NULL,
                                       meanstructure = FALSE,
                                       fixed_x = FALSE,
                                       x_idx = integer(0L),
                                       conditional_x = FALSE) {
  m_y <- unname(as.matrix(m_y))
  n <- nrow(m_y)
  m_s <- cov(m_y) * (n - 1) / n

  if (conditional_x) {
    # conditional.x implies fixed.x: the x-moment rows of the joint Gamma
    # are zeroed, so only the sampling variability of the y-moments (given
    # x) propagates through the Jacobian
    j_mat <- lav_cor_cond_jacobian(
      m_cov = m_s, m_mean = colMeans(m_y),
      x_idx = x_idx, meanstructure = meanstructure
    )
    gamma_full <- lav_samp_gamma(m_y, meanstructure = meanstructure,
                                 x_idx = x_idx, fixed_x = TRUE)
  } else {
    j_mat <- lav_partial_cor_jacobian(m_s, cor_idx = cor_idx,
                                      meanstructure = meanstructure)
    gamma_full <- lav_samp_gamma(m_y, meanstructure = meanstructure,
                                 x_idx = x_idx, fixed_x = fixed_x)
  }

  j_mat %*% gamma_full %*% t(j_mat)
}

# Normal-theory Gamma for a (partial) correlation structure, via the delta
# method: Gamma_cor = J %*% Gamma_full_NT %*% t(J), where Gamma_full_NT is the
# raw NT Gamma (lav_samp_gamma_nt, which also handles fixed.x). If
# conditional_x = TRUE, m_cov/m_mean are the JOINT (y, x) moments and the
# result is the NT Gamma of the standardized conditional moments.
lav_samp_partial_cor_gamma_nt <- function(m_cov, cor_idx = NULL,
                                          meanstructure = FALSE,
                                          fixed_x = FALSE,
                                          x_idx = integer(0L),
                                          conditional_x = FALSE,
                                          m_mean = NULL) {
  if (conditional_x) {
    j_mat <- lav_cor_cond_jacobian(
      m_cov = m_cov, m_mean = m_mean,
      x_idx = x_idx, meanstructure = meanstructure
    )
    gamma_full <- lav_samp_gamma_nt(m_cov = m_cov,
                                    meanstructure = meanstructure,
                                    x_idx = x_idx, fixed_x = TRUE)
  } else {
    j_mat <- lav_partial_cor_jacobian(m_cov, cor_idx = cor_idx,
                                      meanstructure = meanstructure)
    gamma_full <- lav_samp_gamma_nt(m_cov = m_cov,
                                    meanstructure = meanstructure,
                                    x_idx = x_idx, fixed_x = fixed_x)
  }

  j_mat %*% gamma_full %*% t(j_mat)
}
