# - 0.6-13: fix multiple-group UG^2 bug (reported by Gronneberg, Foldnes and
#           Moss) when Satterthwaite = TRUE, ngroups > 1, and eq constraints.
#           Use ug2.old.approach = TRUE to get the old result

lav_test_satorra_bentler <- function(lavobject = NULL,
                                     lavsamplestats = NULL,
                                     lavmodel = NULL,
                                     lavimplied = NULL,
                                     lavoptions = NULL,
                                     lavdata = NULL,
                                     TEST.unscaled = NULL,
                                     E.inv = NULL,
                                     Delta = NULL,
                                     WLS.V = NULL,
                                     Gamma = NULL,
                                     test = "satorra.bentler",
                                     mimic = "lavaan",
                                     method = "original",
                                     ug2.old.approach = FALSE,
                                     return.u = FALSE,
                                     return.ugamma = FALSE) {
  TEST <- list()

  if (!is.null(lavobject)) {
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavoptions <- lavobject@Options
    lavimplied <- lavobject@implied
    lavdata <- lavobject@Data
    TEST$standard <- lavobject@test[[1]]
  } else {
    TEST$standard <- TEST.unscaled
  }
  npar <- lavmodel@nx.free

  # ug2.old.approach
  if (missing(ug2.old.approach)) {
    if (!is.null(lavoptions$ug2.old.approach)) {
      ug2.old.approach <- lavoptions$ug2.old.approach
    } else {
      ug2.old.approach <- FALSE
    }
  }

  # E.inv ok?
  if (length(lavoptions$information) == 1L &&
    length(lavoptions$h1.information) == 1L &&
    length(lavoptions$observed.information) == 1L) {
    E.inv.recompute <- FALSE
  } else if (
    (lavoptions$information[1] == lavoptions$information[2]) &&
      (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
      (lavoptions$information[2] == "expected" ||
        (lavoptions$observed.information[1] ==
          lavoptions$observed.information[2]))) {
    E.inv.recompute <- FALSE
  } else {
    E.inv.recompute <- TRUE
    # change information options
    lavoptions$information[1] <- lavoptions$information[2]
    lavoptions$h1.information[1] <- lavoptions$h1.information[2]
    lavoptions$observed.information[1] <- lavoptions$observed.information[2]
  }
  if (!is.null(E.inv) && !is.null(WLS.V) && !is.null(Delta)) {
    E.inv.recompute <- FALSE # user-provided
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

  if (return.u) {
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
    (is.null(E.inv) || is.null(Delta) || is.null(WLS.V) || E.inv.recompute)) {
    if (mimic == "Mplus" && lavoptions$estimator == "ML") {
      E <- lav_model_information_expected_MLM(
        lavmodel = lavmodel,
        augmented = FALSE, inverted = FALSE,
        lavsamplestats = lavsamplestats, extra = TRUE
      )
    } else {
      E <- lav_model_information(
        lavmodel = lavmodel,
        lavimplied = lavimplied,
        lavsamplestats = lavsamplestats, lavdata = lavdata,
        lavoptions = lavoptions, extra = TRUE
      )
    }
    E.inv <- try(lav_model_information_augment_invert(lavmodel,
      information = E, inverted = TRUE
    ), silent = TRUE)
    if (inherits(E.inv, "try-error")) {
      if (return.ugamma) {
        lav_msg_warn(gettext(
          "could not invert information matrix needed for UGamma"))
        return(NULL)
      } else if (return.u) {
        lav_msg_warn(gettext(
          "could not invert information matrix needed for UfromUGamma"))
        return(NULL)
      } else {
        TEST$standard$stat <- as.numeric(NA)
        TEST$standard$stat.group <- rep(as.numeric(NA), lavdata@ngroups)
        TEST$standard$pvalue <- as.numeric(NA)
        TEST[[test[1]]] <- c(TEST$standard,
          scaling.factor = as.numeric(NA),
          shift.parameter = as.numeric(NA),
          label = character(0)
        )
        lav_msg_warn(gettext("could not invert information matrix needed for
                             robust test statistic"))
        TEST[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when robust test is detected for some but not all models
        return(TEST)
      }
    }
    Delta <- attr(E, "Delta")
    WLS.V <- attr(E, "WLS.V")
  }

  # catch df == 0
  if ((TEST$standard$df == 0L || TEST$standard$df < 0) &&
    !return.u && !return.ugamma) {
    TEST[[test[1]]] <- c(TEST$standard,
      scaling.factor = as.numeric(NA),
      label = character(0)
    )
    TEST[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when robust test is detected for some but not all models
    return(TEST)
  }

  # Gamma
  if (is.null(Gamma)) {
    Gamma <- lavsamplestats@NACOV
    # still NULL? (perhaps estimator = ML)
    if (is.null(Gamma[[1]])) {
      if (!is.null(lavobject)) {
        Gamma <- lav_object_gamma(lavobject, model.based = FALSE)
      } else {
        Gamma <- lav_object_gamma(
          lavobject = NULL,
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavsamplestats = lavsamplestats,
          lavh1 = NULL,
          lavimplied = NULL,
          ADF = TRUE,
          model.based = FALSE
        )
      }
    }
  }

  if (mimic == "Mplus" && lavmodel@categorical) {
    for (g in 1:lavsamplestats@ngroups) {
      Ng <- lavsamplestats@nobs[[g]]
      Gamma[[g]] <- Gamma[[g]] / Ng * (Ng - 1L)
    }
  }

  # ngroups
  ngroups <- lavsamplestats@ngroups

  # mean and variance adjusted?
  Satterthwaite <- FALSE
  if (any(test %in% c("mean.var.adjusted", "scaled.shifted"))) {
    Satterthwaite <- TRUE
  }

  if (npar == 0) {
    # catch npar == 0 (eg baseline model if correlation structure)
    trace.UGamma <- trace.UGamma2 <- U.all <- UG <- as.numeric(NA)
    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    Gamma.f <- Gamma
    for (g in 1:ngroups) {
      Gamma.f[[g]] <- 1 / fg[g] * Gamma[[g]]
    }
    Gamma.all <- lav_matrix_bdiag(Gamma.f)
    UG <- Gamma.all
    trace.UGamma <- sum(diag(Gamma.all))
    trace.UGamma2 <- sum(UG * t(UG))
    out <- list(
      trace.UGamma = trace.UGamma, trace.UGamma2 = trace.UGamma2,
      UGamma = UG, UfromUGamma = U.all
    )
  } else if (method == "original") {
    out <- lav_test_satorra_bentler_trace_original(
      Gamma = Gamma,
      Delta = Delta, WLS.V = WLS.V, E.inv = E.inv,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal, return.u = return.u,
      return.ugamma = return.ugamma,
      ug2.old.approach = ug2.old.approach,
      Satterthwaite = Satterthwaite
    )
  } else if (method == "orthogonal.complement") {
    out <- lav_test_satorra_bentler_trace_complement(
      Gamma = Gamma,
      Delta = Delta, WLS.V = WLS.V, lavmodel = lavmodel,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal,
      return.ugamma = return.ugamma,
      ug2.old.approach = ug2.old.approach,
      Satterthwaite = Satterthwaite
    )
  } else if (method == "ABA") {
    out <- lav_test_satorra_bentler_trace_ABA(
      Gamma = Gamma,
      Delta = Delta, WLS.V = WLS.V, E.inv = E.inv,
      ngroups = ngroups, nobs = lavsamplestats@nobs,
      ntotal = lavsamplestats@ntotal,
      return.ugamma = return.ugamma,
      ug2.old.approach = ug2.old.approach,
      Satterthwaite = Satterthwaite
    )
  } else {
    lav_msg_stop(gettextf("method `%s' not supported", method))
  }
  trace.UGamma <- out$trace.UGamma
  trace.UGamma2 <- out$trace.UGamma2

  if ("satorra.bentler" %in% test) {
    # same df
    df.scaled <- TEST$standard$df

    # scaling factor
    scaling.factor <- trace.UGamma / df.scaled
    if (scaling.factor < 0) scaling.factor <- as.numeric(NA)

    # scaled test statistic per group
    stat.group <- TEST$standard$stat.group / scaling.factor

    # scaled test statistic global
    stat <- sum(stat.group)

    # label
    if (mimic == "Mplus") {
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

    TEST$satorra.bentler <-
      list(
        test = "satorra.bentler",
        stat = stat,
        stat.group = stat.group,
        df = df.scaled,
        pvalue = 1 - pchisq(stat, df.scaled),
        trace.UGamma = trace.UGamma,
        scaling.factor = scaling.factor,
        scaled.test.stat = TEST$standard$stat,
        scaled.test = TEST$standard$test,
        label = label
      )
  }

  if ("mean.var.adjusted" %in% test) {
    if (mimic == "Mplus") {
      df.scaled <- floor(trace.UGamma^2 / trace.UGamma2 + 0.5)
    } else {
      # more precise, fractional df
      df.scaled <- trace.UGamma^2 / trace.UGamma2
    }

    # scaling factor
    scaling.factor <- trace.UGamma / df.scaled
    if (scaling.factor < 0) scaling.factor <- as.numeric(NA)

    if (ug2.old.approach) {
      # scaled test statistic per group
      stat.group <- TEST$standard$stat.group / scaling.factor

      # scaled test statistic global
      stat <- sum(stat.group)
    } else {
      # scaled test statistic per group
      stat.group <- TEST$standard$stat.group / scaling.factor

      # scaled test statistic global
      stat <- TEST$standard$stat / scaling.factor
    }

    # label
    if (mimic == "Mplus") {
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

    TEST$mean.var.adjusted <-
      list(
        test = "mean.var.adjusted",
        stat = stat,
        stat.group = stat.group,
        df = df.scaled,
        pvalue = 1 - pchisq(stat, df.scaled),
        trace.UGamma = trace.UGamma,
        trace.UGamma2 = trace.UGamma2,
        scaling.factor = scaling.factor,
        scaled.test.stat = TEST$standard$stat,
        scaled.test = TEST$standard$test,
        label = label
      )
  }

  if ("scaled.shifted" %in% test) {
    # this is the T3 statistic as used by Mplus 6 and higher
    # see 'Simple Second Order Chi-Square Correction' 2010
    # www.statmodel.com

    # same df
    df.scaled <- TEST$standard$df

    # scaling factor
    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    a <- sqrt(df.scaled / trace.UGamma2)
    scaling.factor <- 1 / a
    if (scaling.factor < 0) scaling.factor <- as.numeric(NA)

    if (ug2.old.approach) {
      # scaling factor
      shift.parameter <- fg * (df.scaled - a * trace.UGamma)
      # scaled test statistic per group
      stat.group <- (TEST$standard$stat.group * a + shift.parameter)
      # scaled test statistic global
      stat <- sum(stat.group)
    } else {
      shift.parameter <- df.scaled - a * trace.UGamma
      stat <- TEST$standard$stat * a + shift.parameter
      stat.group <- TEST$standard$stat.group * a + fg * shift.parameter
    }

    # label
    if (mimic == "Mplus") {
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

    TEST$scaled.shifted <-
      list(
        test = "scaled.shifted",
        stat = stat,
        stat.group = stat.group,
        df = df.scaled,
        pvalue = 1 - pchisq(stat, df.scaled),
        trace.UGamma = trace.UGamma,
        trace.UGamma2 = trace.UGamma2,
        scaling.factor = scaling.factor,
        shift.parameter = shift.parameter,
        scaled.test.stat = TEST$standard$stat,
        scaled.test = TEST$standard$test,
        label = label
      )
  }

  if (return.ugamma) {
    TEST$UGamma <- out$UGamma
  }

  if (return.u) {
    TEST$UfromUGamma <- out$UfromUGamma
  }

  TEST
}

# using the `classical' formula
# UG = Gamma * [V - V Delta E.inv Delta' V']
lav_test_satorra_bentler_trace_original <- function(Gamma = NULL,
                                                    Delta = NULL,
                                                    WLS.V = NULL,
                                                    E.inv = NULL,
                                                    ngroups = NULL,
                                                    nobs = NULL,
                                                    ntotal = NULL,
                                                    return.u = FALSE,
                                                    return.ugamma = FALSE,
                                                    ug2.old.approach = FALSE,
                                                    Satterthwaite = FALSE) {
  # this is what we did <0.6-13: everything per group
  if (ug2.old.approach) {
    UfromUGamma <- UG <- vector("list", ngroups)
    trace.UGamma <- trace.UGamma2 <- rep(as.numeric(NA), ngroups)
    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      Gamma.g <- Gamma[[g]] / fg ## ?? check this
      Delta.g <- Delta[[g]]
      if (is.matrix(WLS.V[[g]])) {
        WLS.Vg <- WLS.V[[g]] * fg
      } else {
        WLS.Vg <- diag(WLS.V[[g]]) * fg
      }

      U <- (WLS.Vg - WLS.Vg %*% Delta[[g]] %*% E.inv %*%
        t(Delta[[g]]) %*% WLS.Vg)
      trace.UGamma[g] <- sum(U * Gamma.g)

      if (return.u) {
        UfromUGamma[[g]] <- U
      }

      UG <- NULL
      if (Satterthwaite || return.ugamma) {
        UG.group <- U %*% Gamma.g
        trace.UGamma2[g] <- sum(UG.group * t(UG.group))
        UG[[g]] <- UG.group
      }
    } # g
    # sum over groups
    trace.UGamma <- sum(trace.UGamma)
    trace.UGamma2 <- sum(trace.UGamma2)
    U.all <- UfromUGamma # group-specific
  } else {
    trace.UGamma <- trace.UGamma2 <- U.all <- UG <- as.numeric(NA)
    fg <- unlist(nobs) / ntotal
    if (Satterthwaite || return.ugamma || return.u) {
      # for trace.UGamma2, we can no longer compute the trace per group
      V.g <- WLS.V
      for (g in 1:ngroups) {
        if (is.matrix(WLS.V[[g]])) {
          V.g[[g]] <- fg[g] * WLS.V[[g]]
        } else {
          V.g[[g]] <- fg[g] * diag(WLS.V[[g]])
        }
      }
      V.all <- lav_matrix_bdiag(V.g)
      Gamma.f <- Gamma
      for (g in 1:ngroups) {
        Gamma.f[[g]] <- 1 / fg[g] * Gamma[[g]]
      }
      Gamma.all <- lav_matrix_bdiag(Gamma.f)
      Delta.all <- do.call("rbind", Delta)
      U.all <- V.all - V.all %*% Delta.all %*% E.inv %*% t(Delta.all) %*% V.all
      UG <- U.all %*% Gamma.all

      trace.UGamma <- sum(U.all * Gamma.all)
      trace.UGamma2 <- sum(UG * t(UG))
    } else {
      # we only need trace.UGamma - this can be done group-specific
      trace.UGamma.group <- numeric(ngroups)
      for (g in 1:ngroups) {
        Gamma.g <- Gamma[[g]] / fg[g]
        Delta.g <- Delta[[g]]
        if (is.matrix(WLS.V[[g]])) {
          WLS.Vg <- WLS.V[[g]] * fg[g]
        } else {
          WLS.Vg <- diag(WLS.V[[g]]) * fg[g]
        }

        U <- (WLS.Vg - WLS.Vg %*% Delta[[g]] %*% E.inv %*%
          t(Delta[[g]]) %*% WLS.Vg)
        trace.UGamma.group[g] <- sum(U * Gamma.g)
      }
      trace.UGamma <- sum(trace.UGamma.group)
    }
  }

  list(
    trace.UGamma = trace.UGamma, trace.UGamma2 = trace.UGamma2,
    UGamma = UG, UfromUGamma = U.all
  )
}

# using the orthogonal complement of Delta: Delta.c
# UG = [ (Delta.c' W Delta.c)^{-1} (Delta.c' Gamma Delta.c)
lav_test_satorra_bentler_trace_complement <- function(Gamma = NULL,
                                                      Delta = NULL,
                                                      WLS.V = NULL,
                                                      lavmodel = NULL,
                                                      ngroups = NULL,
                                                      nobs = NULL,
                                                      ntotal = NULL,
                                                      return.ugamma = FALSE,
                                                      ug2.old.approach = FALSE,
                                                      Satterthwaite = FALSE) {
  # this is what we did <0.6-13: everything per group
  # does not work when ngroups > 1 + equality constraints
  if (ug2.old.approach) {
    UG <- vector("list", ngroups)
    trace.UGamma <- trace.UGamma2 <- rep(as.numeric(NA), ngroups)
    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      Gamma.g <- Gamma[[g]] / fg ## ?? check this
      Delta.g <- Delta[[g]]
      if (is.matrix(WLS.V[[g]])) {
        WLS.Vg <- WLS.V[[g]] * fg
      } else {
        WLS.Vg <- diag(WLS.V[[g]]) * fg
      }

      # handle equality constraints
      # FIXME: inequality constraints are ignored!
      if (lavmodel@eq.constraints) {
        Delta.g <- Delta.g %*% lavmodel@eq.constraints.K
      } else if (.hasSlot(lavmodel, "ceq.simple.only") &&
        lavmodel@ceq.simple.only) {
        Delta.g <- Delta.g %*% lavmodel@ceq.simple.K
      }

      # orthogonal complement of Delta.g
      Delta.c <- lav_matrix_orthogonal_complement(Delta.g)

      ### FIXME: compute WLS.W directly, instead of using solve(WLS.V)

      tmp1 <- solve(t(Delta.c) %*% solve(WLS.Vg) %*% Delta.c)
      tmp2 <- t(Delta.c) %*% Gamma.g %*% Delta.c

      trace.UGamma[g] <- sum(tmp1 * tmp2)
      UG <- NULL
      if (Satterthwaite || return.ugamma) {
        UG.group <- tmp1 %*% tmp2
        trace.UGamma2[g] <- sum(UG.group * t(UG.group))
        UG[[g]] <- UG.group
      }
    }
    # sum over groups
    trace.UGamma <- sum(trace.UGamma)
    trace.UGamma2 <- sum(trace.UGamma2)
  } else {
    trace.UGamma <- trace.UGamma2 <- UG <- as.numeric(NA)
    fg <- unlist(nobs) / ntotal

    V.g <- WLS.V
    for (g in 1:ngroups) {
      if (is.matrix(WLS.V[[g]])) {
        V.g[[g]] <- fg[g] * WLS.V[[g]]
      } else {
        V.g[[g]] <- fg[g] * diag(WLS.V[[g]])
      }
    }
    V.all <- lav_matrix_bdiag(V.g)
    Gamma.f <- Gamma
    for (g in 1:ngroups) {
      Gamma.f[[g]] <- 1 / fg[g] * Gamma[[g]]
    }
    Gamma.all <- lav_matrix_bdiag(Gamma.f)
    Delta.all <- do.call("rbind", Delta)

    # handle equality constraints
    # FIXME: inequality constraints are ignored!
    if (lavmodel@eq.constraints) {
      Delta.all <- Delta.all %*% lavmodel@eq.constraints.K
    } else if (.hasSlot(lavmodel, "ceq.simple.only") &&
      lavmodel@ceq.simple.only) {
      Delta.all <- Delta.all %*% lavmodel@ceq.simple.K
    }

    # orthogonal complement of Delta.g
    Delta.c <- lav_matrix_orthogonal_complement(Delta.all)

    tmp1 <- solve(t(Delta.c) %*% solve(V.all) %*% Delta.c)
    tmp2 <- t(Delta.c) %*% Gamma.all %*% Delta.c

    UG <- tmp1 %*% tmp2

    trace.UGamma <- sum(tmp1 * tmp2)
    trace.UGamma2 <- sum(UG * t(UG))
  }

  list(
    trace.UGamma = trace.UGamma, trace.UGamma2 = trace.UGamma2,
    UGamma = UG
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
lav_test_satorra_bentler_trace_ABA <- function(Gamma = NULL,
                                               Delta = NULL,
                                               WLS.V = NULL,
                                               E.inv = NULL,
                                               ngroups = NULL,
                                               nobs = NULL,
                                               ntotal = NULL,
                                               return.ugamma = FALSE,
                                               ug2.old.approach = FALSE,
                                               Satterthwaite = FALSE) {
  # this is what we did <0.6-13: everything per group
  if (ug2.old.approach) {
    UfromUGamma <- UG <- vector("list", ngroups)
    trace.UGamma <- trace.UGamma2 <- rep(as.numeric(NA), ngroups)

    for (g in 1:ngroups) {
      fg <- nobs[[g]] / ntotal
      Gamma.g <- Gamma[[g]] / fg ## ?? check this
      Delta.g <- Delta[[g]]

      # diagonal WLS.V? we check for this since 0.5-17
      diagonal <- FALSE
      if (is.matrix(WLS.V[[g]])) {
        A1 <- WLS.V[[g]] * fg
        AGA1 <- A1 %*% Gamma.g %*% A1
      } else {
        diagonal <- TRUE
        a1 <- WLS.V[[g]] * fg # numeric vector!
        AGA1 <- Gamma.g * tcrossprod(a1)
      }

      # note: we have AGA1 at the end, to avoid ending up with
      # a transposed matrix (both parts are non-symmetric)
      if (diagonal) {
        UG <- t(Gamma.g * a1) -
          (Delta.g %*% tcrossprod(E.inv, Delta.g) %*% AGA1)
      } else {
        UG <- (Gamma.g %*% A1) -
          (Delta.g %*% tcrossprod(E.inv, Delta.g) %*% AGA1)
      }

      trace.UGamma[g] <- sum(diag(UG))
      if (Satterthwaite) {
        trace.UGamma2[g] <- sum(UG * t(UG))
      }
    }
    # sum over groups
    trace.UGamma <- sum(trace.UGamma)
    trace.UGamma2 <- sum(trace.UGamma2)
  } else {
    trace.UGamma <- trace.UGamma2 <- UG <- as.numeric(NA)
    fg <- unlist(nobs) / ntotal
    if (Satterthwaite || return.ugamma) {
      # for trace.UGamma2, we can no longer compute the trace per group
      V.g <- WLS.V
      for (g in 1:ngroups) {
        if (is.matrix(WLS.V[[g]])) {
          V.g[[g]] <- fg[g] * WLS.V[[g]]
        } else {
          V.g[[g]] <- fg[g] * diag(WLS.V[[g]])
        }
      }
      V.all <- lav_matrix_bdiag(V.g)
      Gamma.f <- Gamma
      for (g in 1:ngroups) {
        Gamma.f[[g]] <- 1 / fg[g] * Gamma[[g]]
      }
      Gamma.all <- lav_matrix_bdiag(Gamma.f)
      Delta.all <- do.call("rbind", Delta)

      AGA1 <- V.all %*% Gamma.all %*% V.all

      UG <- (Gamma.all %*% V.all) -
        (Delta.all %*% tcrossprod(E.inv, Delta.all) %*% AGA1)

      trace.UGamma <- sum(diag(UG))
      trace.UGamma2 <- sum(UG * t(UG))
    } else {
      trace.UGamma.group <- numeric(ngroups)
      for (g in 1:ngroups) {
        fg <- nobs[[g]] / ntotal
        Gamma.g <- Gamma[[g]] / fg ## ?? check this
        Delta.g <- Delta[[g]]

        # diagonal WLS.V? we check for this since 0.5-17
        diagonal <- FALSE
        if (is.matrix(WLS.V[[g]])) {
          A1 <- WLS.V[[g]] * fg
          AGA1 <- A1 %*% Gamma.g %*% A1
        } else {
          diagonal <- TRUE
          a1 <- WLS.V[[g]] * fg # numeric vector!
          AGA1 <- Gamma.g * tcrossprod(a1)
        }

        # note: we have AGA1 at the end, to avoid ending up with
        # a transposed matrix (both parts are non-symmetric)
        if (diagonal) {
          UG <- t(Gamma.g * a1) -
            (Delta.g %*% tcrossprod(E.inv, Delta.g) %*% AGA1)
        } else {
          UG <- (Gamma.g %*% A1) -
            (Delta.g %*% tcrossprod(E.inv, Delta.g) %*% AGA1)
        }

        trace.UGamma.group[g] <- sum(diag(UG))
      } # g
      trace.UGamma <- sum(trace.UGamma.group)
    }
  }

  if (!return.ugamma) {
    UG <- NULL
  }

  list(
    trace.UGamma = trace.UGamma, trace.UGamma2 = trace.UGamma2,
    UGamma = UG
  )
}
