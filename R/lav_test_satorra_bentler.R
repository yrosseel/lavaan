# - 0.6-13: fix multiple-group UG^2 bug (reported by Gronneberg, Foldnes and
#           Moss) when Satterthwaite = TRUE, ngroups > 1, and eq constraints.
#           Use ug2.old.approach = TRUE to get the old result

lav_test_satorra_bentler <- function(lavobject = NULL,
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
                                     return_ugamma = FALSE) {
  test_1 <- list()

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
  if (length(lavoptions$information) == 1L &&
    length(lavoptions$h1.information) == 1L &&
    length(lavoptions$observed.information) == 1L) {
    e_inv_recompute <- FALSE
  } else if (
    (lavoptions$information[1] == lavoptions$information[2]) &&
      (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
      (lavoptions$information[2] == "expected" ||
        (lavoptions$observed.information[1] ==
          lavoptions$observed.information[2]))) {
    e_inv_recompute <- FALSE
  } else {
    e_inv_recompute <- TRUE
    # change information options
    lavoptions$information[1] <- lavoptions$information[2]
    lavoptions$h1.information[1] <- lavoptions$h1.information[2]
    lavoptions$observed.information[1] <- lavoptions$observed.information[2]
  }
  if (!is.null(e_inv) && !is.null(wls_v) && !is.null(delta)) {
    e_inv_recompute <- FALSE # user-provided
  }

  # check test
  if (!all(test %in% c(
    "satorra.bentler",
    "scaled.shifted",
    "mean.var.adjusted"
  ))) {
    lav_msg_warn(gettext(
      "test must be one of `satorra.bentler', `scaled.shifted' or
      `mean.var.adjusted'; will use `satorra.bentler' only"))
    test <- "satorra.bentler"
  }

  if (return_u) {
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
      m_e <- lav_model_information_expected_mlm(
        lavmodel = lavmodel,
        augmented = FALSE, inverted = FALSE,
        lavsamplestats = lavsamplestats, extra = TRUE
      )
    } else {
      m_e <- lav_model_information(
        lavmodel = lavmodel,
        lavimplied = lavimplied,
        lavsamplestats = lavsamplestats, lavdata = lavdata,
        lavoptions = lavoptions, extra = TRUE
      )
    }
    e_inv <- try(lav_model_information_augment_invert(lavmodel,
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
        test_1$standard$stat <- as.numeric(NA)
        test_1$standard$stat.group <- rep(as.numeric(NA), lavdata@ngroups)
        test_1$standard$pvalue <- as.numeric(NA)
        test_1[[test[1]]] <- c(test_1$standard,
          scaling.factor = as.numeric(NA),
          shift.parameter = as.numeric(NA),
          label = character(0)
        )
        lav_msg_warn(gettext("could not invert information matrix needed for
                             robust test statistic"))
        test_1[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when
                       # robust test is detected for some but not all models
        return(test_1)
      }
    }
    delta <- attr(m_e, "Delta")
    wls_v <- attr(m_e, "WLS.V")
  }

  # catch df == 0
  if ((test_1$standard$df == 0L || test_1$standard$df < 0) &&
    !return_u && !return_ugamma) {
    test_1[[test[1]]] <- c(test_1$standard,
      scaling.factor = as.numeric(NA),
      label = character(0)
    )
    test_1[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when
                   # robust test is detected for some but not all models
    return(test_1)
  }

  # Gamma
  if (is.null(m_gamma)) {
    m_gamma <- lavsamplestats@NACOV
    # still NULL? (perhaps estimator = ML)
    if (is.null(m_gamma[[1]])) {
      if (!is.null(lavobject)) {
        m_gamma <- lav_object_gamma(lavobject, model_based = FALSE)
      } else {
        m_gamma <- lav_object_gamma(
          lavobject = NULL,
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavsamplestats = lavsamplestats,
          lavh1 = NULL,
          lavimplied = NULL,
          adf = TRUE,
          model_based = FALSE
        )
      }
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
  if (any(test %in% c("mean.var.adjusted", "scaled.shifted"))) {
    satterthwaite <- TRUE
  }

  if (npar == 0) {
    # catch npar == 0 (eg baseline model if correlation structure)
    trace_ugamma <- trace_ugamma2 <- u_all <- ug <- as.numeric(NA)
    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    gamma_f <- m_gamma
    for (g in 1:ngroups) {
      gamma_f[[g]] <- 1 / fg[g] * m_gamma[[g]]
    }
    gamma_all <- lav_matrix_bdiag(gamma_f)
    ug <- gamma_all
    trace_ugamma <- sum(diag(gamma_all))
    trace_ugamma2 <- sum(ug * t(ug))
    out <- list(
      trace.UGamma = trace_ugamma, trace.UGamma2 = trace_ugamma2,
      UGamma = ug, UfromUGamma = u_all
    )
  } else if (method == "original") {
    out <- lav_test_satorra_bentler_trace_original(
      m_gamma = m_gamma,
      delta = delta, wls_v = wls_v, e_inv = e_inv,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal, return_u = return_u,
      return_ugamma = return_ugamma,
      ug2_old_approach = ug2_old_approach,
      satterthwaite = satterthwaite
    )
  } else if (method == "orthogonal.complement") {
    out <- lav_test_satorra_bentler_trace_complement(
      m_gamma = m_gamma,
      delta = delta, wls_v = wls_v, lavmodel = lavmodel,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal,
      return_ugamma = return_ugamma,
      ug2_old_approach = ug2_old_approach,
      satterthwaite = satterthwaite
    )
  } else if (method == "ABA") {
    out <- lav_test_satorra_bentler_trace_aba(
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

    # label
    if (lavoptions$information.expected.mplus) {
      if (lavoptions$estimator == "ML") {
        label <-
          "Satorra-Bentler correction (Mplus variant)"
      } else if (lavoptions$estimator == "DWLS") {
        label <-
          "Satorra-Bentler correction (WLSM)"
      } else if (lavoptions$estimator == "ULS") {
        label <-
          "Satorra-Bentler correction (ULSM)"
      }
    } else {
      label <- "Satorra-Bentler correction"
    }

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
        label = label
      )
  }

  if ("mean.var.adjusted" %in% test) {
    if (lavoptions$information.expected.mplus) {
      df_scaled <- floor(trace_ugamma^2 / trace_ugamma2 + 0.5)
    } else {
      # more precise, fractional df
      df_scaled <- trace_ugamma^2 / trace_ugamma2
    }

    # scaling factor
    scaling_factor <- trace_ugamma / df_scaled
    if (scaling_factor < 0) scaling_factor <- as.numeric(NA)

    if (ug2_old_approach) {
      # scaled test statistic per group
      stat_group <- test_1$standard$stat.group / scaling_factor

      # scaled test statistic global
      stat <- sum(stat_group)
    } else {
      # scaled test statistic per group
      stat_group <- test_1$standard$stat.group / scaling_factor

      # scaled test statistic global
      stat <- test_1$standard$stat / scaling_factor
    }

    # label
    if (lavoptions$information.expected.mplus) {
      if (lavoptions$estimator == "ML") {
        label <-
          "mean and variance adjusted correction (MLMV)"
      } else if (lavoptions$estimator == "DWLS") {
        label <-
          "mean and variance adjusted correction (WLSMV)"
      } else if (lavoptions$estimator == "ULS") {
        label <-
          "mean and variance adjusted correction (ULSMV)"
      }
    } else {
      label <- "mean and variance adjusted correction"
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
        label = label
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

    # label
    if (lavoptions$information.expected.mplus) {
      if (lavoptions$estimator == "ML") {
        label <-
          "simple second-order correction (MLMV)"
      } else if (lavoptions$estimator == "DWLS") {
        label <-
          "simple second-order correction (WLSMV)"
      } else if (lavoptions$estimator == "ULS") {
        label <-
          "simple second-order correction (ULSMV)"
      }
    } else {
      label <- "simple second-order correction"
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
        label = label
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
lav_test_satorra_bentler_trace_original <- function(m_gamma = NULL,  # nolint
                                                    delta = NULL,
                                                    wls_v = NULL,
                                                    e_inv = NULL,
                                                    ngroups = NULL,
                                                    nobs = NULL,
                                                    ntotal = NULL,
                                                    return_u = FALSE,
                                                    return_ugamma = FALSE,
                                                    ug2_old_approach = FALSE,
                                                    satterthwaite = FALSE) {
  # this is what we did <0.6-13: everything per group
  if (ug2_old_approach) {
    ufrom_ugamma <- ug <- vector("list", ngroups)
    trace_ugamma <- trace_ugamma2 <- rep(as.numeric(NA), ngroups)
    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      gamma_g <- m_gamma[[g]] / fg ## ?? check this
      # delta_g <- delta[[g]]
      if (is.matrix(wls_v[[g]])) {
        wls_vg <- wls_v[[g]] * fg
      } else {
        wls_vg <- diag(wls_v[[g]]) * fg
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
    trace_ugamma <- trace_ugamma2 <- u_all <- ug <- as.numeric(NA)
    fg <- unlist(nobs) / ntotal
    if (satterthwaite || return_ugamma || return_u) {
      # for trace.UGamma2, we can no longer compute the trace per group
      v_g <- wls_v
      for (g in 1:ngroups) {
        if (is.matrix(wls_v[[g]])) {
          v_g[[g]] <- fg[g] * wls_v[[g]]
        } else {
          v_g[[g]] <- fg[g] * diag(wls_v[[g]])
        }
      }
      v_all <- lav_matrix_bdiag(v_g)
      gamma_f <- m_gamma
      for (g in 1:ngroups) {
        gamma_f[[g]] <- 1 / fg[g] * m_gamma[[g]]
      }
      gamma_all <- lav_matrix_bdiag(gamma_f)
      delta_all <- do.call("rbind", delta)
      u_all <- v_all - v_all %*% delta_all %*% e_inv %*% t(delta_all) %*% v_all
      ug <- u_all %*% gamma_all

      trace_ugamma <- sum(u_all * gamma_all)
      trace_ugamma2 <- sum(ug * t(ug))
    } else {
      # we only need trace.UGamma - this can be done group-specific
      trace_ugamma_group <- numeric(ngroups)
      for (g in 1:ngroups) {
        gamma_g <- m_gamma[[g]] / fg[g]
        # delta_g <- delta[[g]]
        if (is.matrix(wls_v[[g]])) {
          wls_vg <- wls_v[[g]] * fg[g]
        } else {
          wls_vg <- diag(wls_v[[g]]) * fg[g]
        }

        u <- (wls_vg - wls_vg %*% delta[[g]] %*% e_inv %*%
          t(delta[[g]]) %*% wls_vg)
        trace_ugamma_group[g] <- sum(u * gamma_g)
      }
      trace_ugamma <- sum(trace_ugamma_group)
    }
  }

  list(
    trace.UGamma = trace_ugamma, trace.UGamma2 = trace_ugamma2,
    UGamma = ug, UfromUGamma = u_all
  )
}

# using the orthogonal complement of Delta: Delta.c
# UG = [ (Delta.c' W Delta.c)^{-1} (Delta.c' Gamma Delta.c)
lav_test_satorra_bentler_trace_complement <- function(m_gamma = NULL,  # nolint
                                                      delta = NULL,
                                                      wls_v = NULL,
                                                      lavmodel = NULL,
                                                      ngroups = NULL,
                                                      nobs = NULL,
                                                      ntotal = NULL,
                                                      return_ugamma = FALSE,
                                                      ug2_old_approach = FALSE,
                                                      satterthwaite = FALSE) {
  # this is what we did <0.6-13: everything per group
  # does not work when ngroups > 1 + equality constraints
  if (ug2_old_approach) {
    ug <- vector("list", ngroups)
    trace_ugamma <- trace_ugamma2 <- rep(as.numeric(NA), ngroups)
    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      gamma_g <- m_gamma[[g]] / fg ## ?? check this
      delta_g <- delta[[g]]
      if (is.matrix(wls_v[[g]])) {
        wls_vg <- wls_v[[g]] * fg
      } else {
        wls_vg <- diag(wls_v[[g]]) * fg
      }

      # handle equality constraints
      # FIXME: inequality constraints are ignored!
      if (lavmodel@eq.constraints) {
        delta_g <- delta_g %*% lavmodel@eq.constraints.K
      } else if (lavmodel@ceq.simple.only) {
        delta_g <- delta_g %*% lavmodel@ceq.simple.K
      }

      # orthogonal complement of Delta.g
      delta_c <- lav_matrix_orthogonal_complement(delta_g)

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
    trace_ugamma <- trace_ugamma2 <- ug <- as.numeric(NA)
    fg <- unlist(nobs) / ntotal

    v_g <- wls_v
    for (g in 1:ngroups) {
      if (is.matrix(wls_v[[g]])) {
        v_g[[g]] <- fg[g] * wls_v[[g]]
      } else {
        v_g[[g]] <- fg[g] * diag(wls_v[[g]])
      }
    }
    v_all <- lav_matrix_bdiag(v_g)
    gamma_f <- m_gamma
    for (g in 1:ngroups) {
      gamma_f[[g]] <- 1 / fg[g] * m_gamma[[g]]
    }
    gamma_all <- lav_matrix_bdiag(gamma_f)
    delta_all <- do.call("rbind", delta)

    # handle equality constraints
    # FIXME: inequality constraints are ignored!
    if (lavmodel@eq.constraints) {
      delta_all <- delta_all %*% lavmodel@eq.constraints.K
    } else if (lavmodel@ceq.simple.only) {
      delta_all <- delta_all %*% lavmodel@ceq.simple.K
    }

    # orthogonal complement of Delta.g
    delta_c <- lav_matrix_orthogonal_complement(delta_all)

    tmp1 <- solve(t(delta_c) %*% solve(v_all) %*% delta_c)
    tmp2 <- t(delta_c) %*% gamma_all %*% delta_c

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
lav_test_satorra_bentler_trace_aba <- function(m_gamma = NULL,       # nolint
                                               m_delta = NULL,
                                               wls_v = NULL,
                                               e_inv = NULL,
                                               ngroups = NULL,
                                               nobs = NULL,
                                               ntotal = NULL,
                                               return_ugamma = FALSE,
                                               ug2_old_approach = FALSE,
                                               satterthwaite = FALSE) {
  # this is what we did <0.6-13: everything per group
  if (ug2_old_approach) {
    # ufromugamma <-
                     ug <- vector("list", ngroups)
    trace_ugamma <- trace_ugamma2 <- rep(as.numeric(NA), ngroups)

    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      gamma_g <- m_gamma[[g]] / fg ## ?? check this
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
        ug <- t(gamma_g * a1) -
          (delta_g %*% tcrossprod(e_inv, delta_g) %*% aga1)
      } else {
        ug <- (gamma_g %*% aa1) -
          (delta_g %*% tcrossprod(e_inv, delta_g) %*% aga1)
      }

      trace_ugamma[g] <- sum(diag(ug))
      if (satterthwaite) {
        trace_ugamma2[g] <- sum(ug * t(ug))
      }
    }
    # sum over groups
    trace_ugamma <- sum(trace_ugamma)
    trace_ugamma2 <- sum(trace_ugamma2)
  } else {
    trace_ugamma <- trace_ugamma2 <- ug <- as.numeric(NA)
    fg <- unlist(nobs) / ntotal
    if (satterthwaite || return_ugamma) {
      # for trace_ugamma2, we can no longer compute the trace per group
      v_g <- wls_v
      for (g in 1:ngroups) {
        if (is.matrix(wls_v[[g]])) {
          v_g[[g]] <- fg[g] * wls_v[[g]]
        } else {
          v_g[[g]] <- fg[g] * diag(wls_v[[g]])
        }
      }
      v_all <- lav_matrix_bdiag(v_g)
      gamma_f <- m_gamma
      for (g in 1:ngroups) {
        gamma_f[[g]] <- 1 / fg[g] * m_gamma[[g]]
      }
      gamma_all <- lav_matrix_bdiag(gamma_f)
      delta_all <- do.call("rbind", m_delta)

      aga1 <- v_all %*% gamma_all %*% v_all

      ug <- (gamma_all %*% v_all) -
        (delta_all %*% tcrossprod(e_inv, delta_all) %*% aga1)

      trace_ugamma <- sum(diag(ug))
      trace_ugamma2 <- sum(ug * t(ug))
    } else {
      trace_ugamma_group <- numeric(ngroups)
      for (g in 1:ngroups) {
        fg <- nobs[[g]] / ntotal
        gamma_g <- m_gamma[[g]] / fg ## ?? check this
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
          ug <- t(gamma_g * a1) -
            (delta_g %*% tcrossprod(e_inv, delta_g) %*% aga1)
        } else {
          ug <- (gamma_g %*% aa1) -
            (delta_g %*% tcrossprod(e_inv, delta_g) %*% aga1)
        }

        trace_ugamma_group[g] <- sum(diag(ug))
      } # g
      trace_ugamma <- sum(trace_ugamma_group)
    }
  }

  if (!return_ugamma) {
    ug <- NULL
  }

  list(
    trace.UGamma = trace_ugamma, trace.UGamma2 = trace_ugamma2,
    UGamma = ug
  )
}
