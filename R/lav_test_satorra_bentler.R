# - 0.6-13: fix multiple-group UG^2 bug (reported by Gronneberg, Foldnes and
#           Moss) when Satterthwaite = TRUE, ngroups > 1, and eq constraints.
#           Use ug2.old.approach = TRUE to get the old result
# - 0.7-1:  the traces are computed by the streamed worker
#           lav_test_ug_trace_stream() (lav_test_utils.R) unless the U or
#           UGamma matrices themselves are needed; shared helpers for the
#           E.inv-recompute check, failure entries and labels

lav_test_sb <- function(lavobject = NULL,
                                     lavsamplestats = NULL,
                                     lavmodel = NULL,
                                     lavimplied = NULL,
                                     lavoptions = NULL,
                                     lavdata = NULL,
                                     test_unscaled = NULL,
                                     e_inv = NULL,
                                     delta = NULL,
                                     wls_v = NULL,
                                     m_gamma = NULL,
                                     test = "satorra.bentler",
                                     method = "original",
                                     ug2_old_approach = FALSE,
                                     return_u = FALSE,
                                     return_ugamma = FALSE,
                                     gamma_full = NULL) {
  test_1 <- list()

  # if a full (cross-group) Gamma is supplied, we must use the 'original'
  # method (the only path that consumes a single assembled Gamma matrix)
  if (!is.null(gamma_full)) {
    method <- "original"
  }

  if (!is.null(lavobject)) {
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavoptions <- lavobject@Options
    lavimplied <- lavobject@implied
    lavdata <- lavobject@Data
    test_1$standard <- lavobject@test[[1]]
  } else {
    test_1$standard <- test_unscaled
  }
  npar <- lavmodel@nx.free

  # ug2.old.approach
  if (missing(ug2_old_approach)) {
    if (!is.null(lavoptions$ug2.old.approach)) {
      ug2_old_approach <- lavoptions$ug2.old.approach
    } else {
      ug2_old_approach <- FALSE
    }
  }

  # E.inv ok?
  lavoptions <- lav_test_change_options(lavoptions)
  e_inv_recompute <- attr(lavoptions, "recompute")
  attr(lavoptions, "recompute") <- NULL
  if (!is.null(e_inv) && !is.null(wls_v) && !is.null(delta)) {
    e_inv_recompute <- FALSE # user-provided
  }

  # check test
  if (!all(test %in% c(
    "satorra.bentler",
    "scaled.shifted",
    "mean.var.adjusted",
    "mean.var.adjusted.corrected",
    "scaled.shifted.corrected"
  ))) {
    lav_msg_warn(gettext(
      "test must be one of `satorra.bentler', `scaled.shifted',
      `mean.var.adjusted', `mean.var.adjusted.corrected' or
      `scaled.shifted.corrected'; will use `satorra.bentler' only"))
    test <- "satorra.bentler"
  }

  # corrected trace (Hayakawa 2018) tests need the U matrix and casewise data
  corrected_trace <- any(test %in% c(
    "mean.var.adjusted.corrected",
    "scaled.shifted.corrected"
  ))
  if (corrected_trace) {
    lav_test_hayakawa_check(
      lavoptions = lavoptions, lavdata = lavdata, lavmodel = lavmodel
    )
  }

  if (return_u || corrected_trace) {
    method <- "original"
  }

  # check method
  if (!all(method %in% c("original", "orthogonal.complement", "ABA"))) {
    lav_msg_warn(gettext("method must be one of `original', `ABA',
                         `orthogonal.complement'; will use `ABA'"))
    method <- "original"
  }

  # do we have E.inv, Delta, WLS.V?
  if (npar > 0L &&
    (is.null(e_inv) || is.null(delta) || is.null(wls_v) || e_inv_recompute)) {
    if (lavoptions$information.expected.mplus && lavoptions$estimator == "ML") {
      m_e <- lav_model_info_expected_mlm(
        lavmodel = lavmodel,
        augmented = FALSE, inverted = FALSE,
        lavsamplestats = lavsamplestats, extra = TRUE
      )
    } else {
      m_e <- lav_model_info(
        lavmodel = lavmodel,
        lavimplied = lavimplied,
        lavsamplestats = lavsamplestats, lavdata = lavdata,
        lavoptions = lavoptions, extra = TRUE
      )
    }
    e_inv <- try(lav_model_info_augment_invert(lavmodel,
      information = m_e, inverted = TRUE
    ), silent = TRUE)
    if (inherits(e_inv, "try-error")) {
      if (return_ugamma) {
        lav_msg_warn(gettext(
          "could not invert information matrix needed for UGamma"))
        return(NULL)
      } else if (return_u) {
        lav_msg_warn(gettext(
          "could not invert information matrix needed for UfromUGamma"))
        return(NULL)
      } else {
        na_out <- lav_test_scaled_na(test_1$standard,
          test1 = test[1],
          ngroups = lavdata@ngroups
        )
        test_1$standard <- na_out$standard
        test_1[[test[1]]] <- na_out$failed
        lav_msg_warn(gettext("could not invert information matrix needed for
                             robust test statistic"))
        return(test_1)
      }
    }
    delta <- attr(m_e, "Delta")
    wls_v <- attr(m_e, "WLS.V")
  }

  # catch df == 0
  if ((test_1$standard$df == 0L || test_1$standard$df < 0) &&
    !return_u && !return_ugamma) {
    na_out <- lav_test_scaled_na(test_1$standard,
      test1 = test[1],
      na_standard = FALSE, shift = FALSE
    )
    test_1[[test[1]]] <- na_out$failed
    return(test_1)
  }

  # Gamma
  if (is.null(m_gamma)) {
    # stored NACOV, else recompute with the fit-time settings
    if (!is.null(lavobject)) {
      m_gamma <- lav_gamma_used(lavobject)
    } else {
      m_gamma <- lav_gamma_used(
        lavdata = lavdata,
        lavoptions = lavoptions,
        lavsamplestats = lavsamplestats,
        lavh1 = NULL,
        lavimplied = NULL
      )
    }
  }

  if (lavoptions$information.expected.mplus && lavmodel@categorical) {
    for (g in 1:lavsamplestats@ngroups) {
      ng <- lavsamplestats@nobs[[g]]
      m_gamma[[g]] <- m_gamma[[g]] / ng * (ng - 1L)
    }
  }

  # ngroups
  ngroups <- lavsamplestats@ngroups

  # mean and variance adjusted?
  satterthwaite <- FALSE
  if (any(test %in% c(
    "mean.var.adjusted", "scaled.shifted",
    "mean.var.adjusted.corrected", "scaled.shifted.corrected"
  ))) {
    satterthwaite <- TRUE
  }

  if (npar == 0) {
    # catch npar == 0 (eg baseline model if correlation structure)
    u_all <- as.numeric(NA)
    # Gamma_g / fg: the Cov(sqrt(ntotal) s_g) stacking convention (see the
    # SCALING CONVENTIONS note in lav_samplestats_gamma.R)
    gamma_f <- lav_gamma_rescale_ntotal(m_gamma,
      nobs = lavsamplestats@nobs, ntotal = lavsamplestats@ntotal
    )
    gamma_all <- lav_mat_bdiag(gamma_f)
    ug <- gamma_all
    trace_ugamma <- sum(diag(gamma_all))
    trace_ugamma2 <- sum(ug * t(ug))
    out <- list(
      trace.UGamma = trace_ugamma, trace.UGamma2 = trace_ugamma2,
      UGamma = ug, UfromUGamma = u_all
    )
  } else if (method == "original") {
    out <- lav_test_sb_trace_original(
      m_gamma = m_gamma,
      delta = delta, wls_v = wls_v, e_inv = e_inv,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal, return_u = return_u || corrected_trace,
      return_ugamma = return_ugamma,
      ug2_old_approach = ug2_old_approach,
      satterthwaite = satterthwaite,
      gamma_full = gamma_full
    )
  } else if (method == "orthogonal.complement") {
    out <- lav_test_sb_trace_complement(
      m_gamma = m_gamma,
      delta = delta, wls_v = wls_v, lavmodel = lavmodel,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal,
      return_ugamma = return_ugamma,
      ug2_old_approach = ug2_old_approach,
      satterthwaite = satterthwaite
    )
  } else if (method == "ABA") {
    out <- lav_test_sb_trace_aba(
      m_gamma = m_gamma,
      m_delta = delta, wls_v = wls_v, e_inv = e_inv,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal,
      return_ugamma = return_ugamma,
      ug2_old_approach = ug2_old_approach,
      satterthwaite = satterthwaite
    )
  } else {
    lav_msg_stop(gettextf("method `%s' not supported", method))
  }
  trace_ugamma <- out$trace.UGamma
  trace_ugamma2 <- out$trace.UGamma2

  # corrected traces (Hayakawa 2018) from the casewise moment vectors
  if (corrected_trace) {
    u_mat <- out$UfromUGamma
    if (is.list(u_mat)) {
      u_mat <- u_mat[[1]] # single group only (checked above)
    }
    if (is.matrix(u_mat)) {
      trace_c <- lav_test_hayakawa_trace2(
        u = u_mat, m_y = lavdata@X[[1]],
        meanstructure = lavmodel@meanstructure
      )
      trace_ugamma_c <- trace_c$trace.UGamma
      trace_ugamma2_c <- trace_c$trace.UGamma2
    } else {
      trace_ugamma_c <- trace_ugamma2_c <- as.numeric(NA)
    }
    if (!is.finite(trace_ugamma2_c) || trace_ugamma2_c <= 0 ||
      !is.finite(trace_ugamma_c) || trace_ugamma_c <= 0) {
      lav_msg_warn(gettext(
        "corrected trace estimates are not positive; the corrected
        adjusted test statistics will be NA."))
      trace_ugamma_c <- trace_ugamma2_c <- as.numeric(NA)
    }
  }

  # Mplus-flavored labels?
  mplus_flag <- lavoptions$information.expected.mplus

  if ("satorra.bentler" %in% test) {
    # same df
    df_scaled <- test_1$standard$df

    # scaling factor
    scaling_factor <- trace_ugamma / df_scaled
    if (scaling_factor < 0) scaling_factor <- as.numeric(NA)

    # scaled test statistic per group
    stat_group <- test_1$standard$stat.group / scaling_factor

    # scaled test statistic global
    stat <- sum(stat_group)

    test_1$satorra.bentler <-
      list(
        test = "satorra.bentler",
        stat = stat,
        stat.group = stat_group,
        df = df_scaled,
        pvalue = 1 - pchisq(stat, df_scaled),
        trace.UGamma = trace_ugamma,
        scaling.factor = scaling_factor,
        scaled.test.stat = test_1$standard$stat,
        scaled.test = test_1$standard$test,
        label = lav_test_scaled_label(
          "satorra.bentler",
          lavoptions$estimator, mplus_flag
        )
      )
  }

  if ("mean.var.adjusted" %in% test) {
    if (mplus_flag) {
      df_scaled <- floor(trace_ugamma^2 / trace_ugamma2 + 0.5)
    } else {
      # more precise, fractional df
      df_scaled <- trace_ugamma^2 / trace_ugamma2
    }

    # scaling factor
    scaling_factor <- trace_ugamma / df_scaled
    if (scaling_factor < 0) scaling_factor <- as.numeric(NA)

    # scaled test statistic per group
    stat_group <- test_1$standard$stat.group / scaling_factor

    # scaled test statistic global
    if (ug2_old_approach) {
      stat <- sum(stat_group)
    } else {
      stat <- test_1$standard$stat / scaling_factor
    }

    test_1$mean.var.adjusted <-
      list(
        test = "mean.var.adjusted",
        stat = stat,
        stat.group = stat_group,
        df = df_scaled,
        pvalue = 1 - pchisq(stat, df_scaled),
        trace.UGamma = trace_ugamma,
        trace.UGamma2 = trace_ugamma2,
        scaling.factor = scaling_factor,
        scaled.test.stat = test_1$standard$stat,
        scaled.test = test_1$standard$test,
        label = lav_test_scaled_label(
          "mean.var.adjusted",
          lavoptions$estimator, mplus_flag
        )
      )
  }

  if ("scaled.shifted" %in% test) {
    # this is the T3 statistic as used by Mplus 6 and higher
    # see 'Simple Second Order Chi-Square Correction' 2010
    # www.statmodel.com

    # same df
    df_scaled <- test_1$standard$df

    # scaling factor
    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    a <- sqrt(df_scaled / trace_ugamma2)
    if (isTRUE(a < 0) || is.nan(a)) a <- as.numeric(NA)
    scaling_factor <- 1 / a
    if (isTRUE(scaling_factor < 0)) scaling_factor <- as.numeric(NA)

    if (ug2_old_approach) {
      # scaling factor
      shift_parameter <- fg * (df_scaled - a * trace_ugamma)
      # scaled test statistic per group
      stat_group <- (test_1$standard$stat.group * a + shift_parameter)
      # scaled test statistic global
      stat <- sum(stat_group)
    } else {
      shift_parameter <- df_scaled - a * trace_ugamma
      stat <- test_1$standard$stat * a + shift_parameter
      stat_group <- test_1$standard$stat.group * a + fg * shift_parameter
    }

    test_1$scaled.shifted <-
      list(
        test = "scaled.shifted",
        stat = stat,
        stat.group = stat_group,
        df = df_scaled,
        pvalue = 1 - pchisq(stat, df_scaled),
        trace.UGamma = trace_ugamma,
        trace.UGamma2 = trace_ugamma2,
        scaling.factor = scaling_factor,
        shift.parameter = shift_parameter,
        scaled.test.stat = test_1$standard$stat,
        scaled.test = test_1$standard$test,
        label = lav_test_scaled_label(
          "scaled.shifted",
          lavoptions$estimator, mplus_flag
        )
      )
  }

  if ("mean.var.adjusted.corrected" %in% test) {
    # Hayakawa (2018): mean and variance adjusted test with the unbiased
    # (Srivastava 2005; Himeno & Yamada 2014) estimator of tr(UGamma^2)
    df_scaled <- trace_ugamma_c^2 / trace_ugamma2_c

    # scaling factor
    scaling_factor <- trace_ugamma2_c / trace_ugamma_c
    if (isTRUE(scaling_factor < 0)) scaling_factor <- as.numeric(NA)

    # scaled test statistic per group and global
    stat_group <- test_1$standard$stat.group / scaling_factor
    stat <- test_1$standard$stat / scaling_factor

    test_1$mean.var.adjusted.corrected <-
      list(
        test = "mean.var.adjusted.corrected",
        stat = stat,
        stat.group = stat_group,
        df = df_scaled,
        pvalue = 1 - pchisq(stat, df_scaled),
        trace.UGamma = trace_ugamma_c,
        trace.UGamma2 = trace_ugamma2_c,
        trace.UGamma2.naive = trace_ugamma2,
        scaling.factor = scaling_factor,
        scaled.test.stat = test_1$standard$stat,
        scaled.test = test_1$standard$test,
        label =
          "mean and variance adjusted correction (corrected trace)"
      )
  }

  if ("scaled.shifted.corrected" %in% test) {
    # scaled and shifted test with the corrected estimator of tr(UGamma^2)
    # (not in Hayakawa 2018, but the same substitution)
    df_scaled <- test_1$standard$df

    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    a <- sqrt(df_scaled / trace_ugamma2_c)
    if (isTRUE(a < 0) || is.nan(a)) a <- as.numeric(NA)
    scaling_factor <- 1 / a
    if (isTRUE(scaling_factor < 0)) scaling_factor <- as.numeric(NA)

    shift_parameter <- df_scaled - a * trace_ugamma_c
    stat <- test_1$standard$stat * a + shift_parameter
    stat_group <- test_1$standard$stat.group * a + fg * shift_parameter

    test_1$scaled.shifted.corrected <-
      list(
        test = "scaled.shifted.corrected",
        stat = stat,
        stat.group = stat_group,
        df = df_scaled,
        pvalue = 1 - pchisq(stat, df_scaled),
        trace.UGamma = trace_ugamma_c,
        trace.UGamma2 = trace_ugamma2_c,
        trace.UGamma2.naive = trace_ugamma2,
        scaling.factor = scaling_factor,
        shift.parameter = shift_parameter,
        scaled.test.stat = test_1$standard$stat,
        scaled.test = test_1$standard$test,
        label = "simple second-order correction (corrected trace)"
      )
  }

  if (return_ugamma) {
    test_1$UGamma <- out$UGamma
  }

  if (return_u) {
    test_1$UfromUGamma <- out$UfromUGamma
  }

  test_1
}

# using the `classical' formula
# UG = Gamma * [V - V Delta E.inv Delta' V']
#
# if only the traces are needed (no return_u/return_ugamma, no cross-group
# gamma_full), they are computed by the streamed worker
# lav_test_ug_trace_stream() without ever forming U or U %*% Gamma
# (this matters in the categorical setting -- and whenever the number of
# variables grows -- where pstar is large)
lav_test_sb_trace_original <- function(m_gamma = NULL,
                                                    delta = NULL,
                                                    wls_v = NULL,
                                                    e_inv = NULL,
                                                    ngroups = NULL,
                                                    nobs = NULL,
                                                    ntotal = NULL,
                                                    return_u = FALSE,
                                                    return_ugamma = FALSE,
                                                    ug2_old_approach = FALSE,
                                                    satterthwaite = FALSE,
                                                    gamma_full = NULL) {
  # gamma_full: an already-assembled (full) Gamma matrix for the *stacked*
  # statistics of all groups. When supplied (e.g. by sam() with across-group
  # constraints, where Gamma has nonzero cross-group blocks), it is used as-is
  # instead of block-diagonalizing the per-group m_gamma list. This forces the
  # full (non per-group) computation path.
  fg <- unlist(nobs) / ntotal

  # fast path: only the traces are needed
  if (!return_u && !return_ugamma && is.null(gamma_full)) {
    out <- lav_test_ug_trace_stream(
      m_gamma = m_gamma, delta = delta, wls_v = wls_v, e_inv = e_inv,
      fg = fg, satterthwaite = satterthwaite,
      old_approach = ug2_old_approach
    )
    return(list(
      trace.UGamma = out$trace.UGamma, trace.UGamma2 = out$trace.UGamma2,
      UGamma = as.numeric(NA), UfromUGamma = as.numeric(NA)
    ))
  }

  # explicit path: the U and/or UGamma matrices themselves are needed
  # (lavInspect, Hayakawa corrected traces, FMG eigenvalues, gamma_full)
  if (ug2_old_approach) {
    # this is what we did <0.6-13: everything per group
    ufrom_ugamma <- ug <- vector("list", ngroups)
    trace_ugamma <- trace_ugamma2 <- rep(as.numeric(NA), ngroups)
    for (g in 1:ngroups) {
      gamma_g <- m_gamma[[g]] / fg[g]
      if (is.matrix(wls_v[[g]])) {
        wls_vg <- wls_v[[g]] * fg[g]
      } else {
        wls_vg <- diag(wls_v[[g]]) * fg[g]
      }

      u <- (wls_vg - wls_vg %*% delta[[g]] %*% e_inv %*%
        t(delta[[g]]) %*% wls_vg)
      trace_ugamma[g] <- sum(u * gamma_g)

      if (return_u) {
        ufrom_ugamma[[g]] <- u
      }

      ug <- NULL
      if (satterthwaite || return_ugamma) {
        ug_group <- u %*% gamma_g
        trace_ugamma2[g] <- sum(ug_group * t(ug_group))
        ug[[g]] <- ug_group
      }
    } # g
    # sum over groups
    trace_ugamma <- sum(trace_ugamma)
    trace_ugamma2 <- sum(trace_ugamma2)
    u_all <- ufrom_ugamma # group-specific
  } else {
    stacked <- lav_test_ug_stack(
      m_gamma = m_gamma, delta = delta, wls_v = wls_v,
      nobs = nobs, ntotal = ntotal, gamma_full = gamma_full
    )
    u_all <- stacked$v_all - stacked$v_all %*% stacked$delta_all %*%
      e_inv %*% t(stacked$delta_all) %*% stacked$v_all
    ug <- u_all %*% stacked$gamma_all

    trace_ugamma <- sum(u_all * stacked$gamma_all)
    trace_ugamma2 <- sum(ug * t(ug))
  }

  list(
    trace.UGamma = trace_ugamma, trace.UGamma2 = trace_ugamma2,
    UGamma = ug, UfromUGamma = u_all
  )
}

# using the orthogonal complement of Delta: Delta.c
# UG = [ (Delta.c' W Delta.c)^{-1} (Delta.c' Gamma Delta.c)
lav_test_sb_trace_complement <- function(m_gamma = NULL,
                                                      delta = NULL,
                                                      wls_v = NULL,
                                                      lavmodel = NULL,
                                                      ngroups = NULL,
                                                      nobs = NULL,
                                                      ntotal = NULL,
                                                      return_ugamma = FALSE,
                                                      ug2_old_approach = FALSE,
                                                      satterthwaite = FALSE) {
  # handle equality constraints
  # FIXME: inequality constraints are ignored!
  eq_basis <- lav_con_eq_basis(lavmodel)

  # this is what we did <0.6-13: everything per group
  # does not work when ngroups > 1 + equality constraints
  if (ug2_old_approach) {
    ug <- vector("list", ngroups)
    trace_ugamma <- trace_ugamma2 <- rep(as.numeric(NA), ngroups)
    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      gamma_g <- m_gamma[[g]] / fg
      delta_g <- delta[[g]]
      if (is.matrix(wls_v[[g]])) {
        wls_vg <- wls_v[[g]] * fg
      } else {
        wls_vg <- diag(wls_v[[g]]) * fg
      }

      if (!is.null(eq_basis)) {
        delta_g <- delta_g %*% eq_basis
      }

      # orthogonal complement of Delta.g
      delta_c <- lav_mat_ortho_complement(delta_g)

      ### FIXME: compute WLS.W directly, instead of using solve(WLS.V)

      tmp1 <- solve(t(delta_c) %*% solve(wls_vg) %*% delta_c)
      tmp2 <- t(delta_c) %*% gamma_g %*% delta_c

      trace_ugamma[g] <- sum(tmp1 * tmp2)
      ug <- NULL
      if (satterthwaite || return_ugamma) {
        ug_group <- tmp1 %*% tmp2
        trace_ugamma2[g] <- sum(ug_group * t(ug_group))
        ug[[g]] <- ug_group
      }
    }
    # sum over groups
    trace_ugamma <- sum(trace_ugamma)
    trace_ugamma2 <- sum(trace_ugamma2)
  } else {
    trace_ugamma2 <- ug <- as.numeric(NA)
    stacked <- lav_test_ug_stack(
      m_gamma = m_gamma, delta = delta, wls_v = wls_v,
      nobs = nobs, ntotal = ntotal
    )
    delta_all <- stacked$delta_all
    if (!is.null(eq_basis)) {
      delta_all <- delta_all %*% eq_basis
    }

    # orthogonal complement of Delta
    delta_c <- lav_mat_ortho_complement(delta_all)

    tmp1 <- solve(t(delta_c) %*% solve(stacked$v_all) %*% delta_c)
    tmp2 <- t(delta_c) %*% stacked$gamma_all %*% delta_c

    ug <- tmp1 %*% tmp2

    trace_ugamma <- sum(tmp1 * tmp2)
    trace_ugamma2 <- sum(ug * t(ug))
  }

  list(
    trace.UGamma = trace_ugamma, trace.UGamma2 = trace_ugamma2,
    UGamma = ug
  )
}

# using the ABA form
# UG = Gamma %*% [V - V %*% Delta %*% E.inv %*% tDelta %*% V]
#    = Gamma %*% V  - Gamma %*% V %*% Delta %*% E.inv %*% tDelta %*% V
#    = Gamma %*% A1 - Gamma %*% A1 %*% Delta %*% E.inv %*% tDelta %*% A1
# (define AGA1 := A1 %*% Gamma %*% A1)
# Note this is not identical to 'B1', (model-based) first-order information
#
#    = A1.inv %*% A1 %*% Gamma %*% A1 -
#      A1.inv %*% A1 %*% Gamma %*% A1 %*% Delta %*% E.inv %*% tDelta %*% A1
#
#    = A1.inv %*% AGA1 -
#      A1.inv %*% AGA1 %*% Delta %*% E.inv %*% tDelta %*% A1
#
# if only the trace is needed, we can use reduce the rhs (after the minus)
# to AGA1 %*% Delta %*% E.inv %*% tDelta (eliminating A1 and A1.inv)

# we write it like this to highlight the connection with MLR
#
lav_test_sb_trace_aba <- function(m_gamma = NULL,
                                               m_delta = NULL,
                                               wls_v = NULL,
                                               e_inv = NULL,
                                               ngroups = NULL,
                                               nobs = NULL,
                                               ntotal = NULL,
                                               return_ugamma = FALSE,
                                               ug2_old_approach = FALSE,
                                               satterthwaite = FALSE) {
  # per-group computation: either the <0.6-13 approach (also for the
  # UGamma^2 trace), or (default) for the first trace only
  if (ug2_old_approach || !(satterthwaite || return_ugamma)) {
    ug <- vector("list", ngroups)
    trace_ugamma <- trace_ugamma2 <- rep(as.numeric(NA), ngroups)

    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      gamma_g <- m_gamma[[g]] / fg
      delta_g <- m_delta[[g]]

      # diagonal wls_v? we check for this since 0.5-17
      diagonal <- FALSE
      if (is.matrix(wls_v[[g]])) {
        aa1 <- wls_v[[g]] * fg
        aga1 <- aa1 %*% gamma_g %*% aa1
      } else {
        diagonal <- TRUE
        a1 <- wls_v[[g]] * fg # numeric vector!
        aga1 <- gamma_g * tcrossprod(a1)
      }

      # note: we have aga1 at the end, to avoid ending up with
      # a transposed matrix (both parts are non-symmetric)
      if (diagonal) {
        ug_group <- t(gamma_g * a1) -
          (delta_g %*% tcrossprod(e_inv, delta_g) %*% aga1)
      } else {
        ug_group <- (gamma_g %*% aa1) -
          (delta_g %*% tcrossprod(e_inv, delta_g) %*% aga1)
      }

      trace_ugamma[g] <- sum(diag(ug_group))
      if (ug2_old_approach && satterthwaite) {
        trace_ugamma2[g] <- sum(ug_group * t(ug_group))
      }
    } # g
    # sum over groups
    trace_ugamma <- sum(trace_ugamma)
    if (ug2_old_approach && satterthwaite) {
      trace_ugamma2 <- sum(trace_ugamma2)
    } else {
      trace_ugamma2 <- as.numeric(NA)
    }
    ug <- as.numeric(NA)
  } else {
    # for trace_ugamma2, we can no longer compute the trace per group
    stacked <- lav_test_ug_stack(
      m_gamma = m_gamma, delta = m_delta, wls_v = wls_v,
      nobs = nobs, ntotal = ntotal
    )
    aga1 <- stacked$v_all %*% stacked$gamma_all %*% stacked$v_all

    ug <- (stacked$gamma_all %*% stacked$v_all) -
      (stacked$delta_all %*% tcrossprod(e_inv, stacked$delta_all) %*% aga1)

    trace_ugamma <- sum(diag(ug))
    trace_ugamma2 <- sum(ug * t(ug))
  }

  if (!return_ugamma) {
    ug <- NULL
  }

  list(
    trace.UGamma = trace_ugamma, trace.UGamma2 = trace_ugamma2,
    UGamma = ug
  )
}
