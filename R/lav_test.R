# chi-square test statistic:
# comparing the current model versus the saturated/unrestricted model
# TDJ 9 April 2024: Add a (hidden) function to update the @test slot
#                   when the user provides a custom h1 model.  Called by
#                   lav_object_summary and lav_fit(_check_baseline)

lavTest <- function(lavobject, test = "standard",               # nolint
                    scaled_test = "standard",
                    output = "list", drop_list_single = TRUE,
                  ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  # check object
  lavobject <- lav_object_check_version(lavobject)

  # check output
  output_valid <- c("list", "text")
  if (!any(output == output_valid)) {
      lav_msg_stop(gettextf(
        "%1$s argument must be either %2$s",
        "output", lav_msg_view(output_valid, "or")
      ))
  }

  # extract 'test' slot
  test_1 <- lavobject@test

  # backwards compatibility:
  # if (length(TEST) > 0L && is.null(names(TEST))) {
  #   names(TEST) <- sapply(TEST, "[[", "test")
  # }

  # which test?
  if (!missing(test)) {
    # check 'test'
    if (!is.character(test)) {
      lav_msg_stop(
        gettextf("%s should be a character string.", "test"))
    } else {
      test <- lav_test_rename(test, check = TRUE)
      # always include "standard"
      test <- unique(c("standard", test))
    }

    # check scaled_test
    if (!missing(scaled_test)) {
      if (!is.character(scaled_test)) {
        lav_msg_stop(
          gettextf("%s should be a character string.", "scaled_test"))
      } else {
        scaled_test <- lav_test_rename(scaled_test, check = TRUE)
      }

      # merge
      test <- unique(c(test, scaled_test))
    }

    # "standard" must always be first
    standard_idx <- which(test == "standard")
    if (length(standard_idx) > 0L && standard_idx != 1L) {
      test <- c("standard", test[-standard_idx])
    }

    if (test[1] == "none") {
      return(list())
    } else if (any(test %in% c("bootstrap", "bollen.stine"))) {
      lav_msg_stop(gettext(
      "please use lavBootstrap() to obtain a bootstrap based test statistic."
      ))
    }

    # check if we already have it:
    if (all(test %in% names(test_1))) {
      info_attr <- attr(test_1, "info")
      test_idx <- which(names(test_1) %in% test)
      test_1 <- test_1[test_idx]
      attr(test_1, "info") <- info_attr
    } else {
      # redo ALL of them, even if already have some in TEST
      # later, we will allow to also change the options (like information)
      # and this should be reflected in the 'info' attribute

      # fill-in test in Options slot
      lavobject@Options$test <- test

      # fill-in scaled_test in Options slot
      lavobject@Options$scaled.test <- scaled_test

      # get requested test statistics
      test_1 <- lav_model_test(lavobject = lavobject)
    }
  }

  if (output == "list") {
    # remove 'info' attribute
    attr(test_1, "info") <- NULL

    # select only those that were requested (eg remove standard)
    test_idx <- which(names(test_1) %in% test)
    test_1 <- test_1[test_idx]

    # if only 1 test, drop outer list
    if (length(test_1) == 1L && drop_list_single) {
      test_1 <- test_1[[1]]
    }

    return(test_1)
  } else {
    lav_test_print(test_1)
  }

  invisible(test_1)
}

# allow for 'flexible' names for the test statistics
# 0.6-13: if multiple names, order them in such a way
#         that the 'scaled' variants appear after the others
lav_test_rename <- function(test, check = FALSE) {
  test <- tolower(test)

  if (length(target_idx <- which(test %in%
    c("standard", "chisq", "chi", "chi-square", "chi.square"))) > 0L) {
    test[target_idx] <- "standard"
  }
  if (length(target_idx <- which(test %in%
    c(
      "satorra", "sb", "satorra.bentler", "satorra-bentler",
      "m.adjusted", "m", "mean.adjusted", "mean-adjusted"
    ))) > 0L) {
    test[target_idx] <- "satorra.bentler"
  }
  if (length(target_idx <- which(test %in%
    c("yuan", "yb", "yuan.bentler", "yuan-bentler"))) > 0L) {
    test[target_idx] <- "yuan.bentler"
  }
  if (length(target_idx <- which(test %in%
    c(
      "yuan.bentler.mplus", "yuan-bentler.mplus",
      "yuan-bentler-mplus"
    ))) > 0L) {
    test[target_idx] <- "yuan.bentler.mplus"
  }
  if (length(target_idx <- which(test %in%
    c("yuan.chan", "yuan-chan", "yuan_chan", "yuanchan", "yc"))) > 0L) {
    test[target_idx] <- "yuan.chan"
  }
  if (length(target_idx <- which(test %in%
    c(
      "mean.var.adjusted", "mean-var-adjusted", "mv", "second.order",
      "satterthwaite", "mv.adjusted"
    ))) > 0L) {
    test[target_idx] <- "mean.var.adjusted"
  }
  if (length(target_idx <- which(test %in%
    c(
      "mplus6", "scale.shift", "scaled.shifted",
      "scaled-shifted"
    ))) > 0L) {
    test[target_idx] <- "scaled.shifted"
  }
  if (length(target_idx <- which(test %in%
    c("bootstrap", "boot", "bollen.stine", "bollen-stine"))) > 0L) {
    test[target_idx] <- "bollen.stine"
  }
  if (length(target_idx <- which(test %in%
    c(
      "browne", "residual", "residuals", "browne.residual",
      "browne.residuals", "residual-based", "residual.based",
      "browne.residuals.adf", "browne.residual.adf"
    ))) > 0L) {
    test[target_idx] <- "browne.residual.adf"
  }
  if (length(target_idx <- which(test %in%
    c("browne.residuals.nt", "browne.residual.nt",
      "browne.nt.residuals", "browne.nt.residual"))) > 0L) {
    test[target_idx] <- "browne.residual.nt"
  }
  if (length(target_idx <- which(test %in%
    c("browne.residual.adf.model", "browne.residuals.adf.model",
      "browne.residual.model.adf", "browne.residuals.model.adf"))) > 0L) {
    test[target_idx] <- "browne.residual.adf.model"
  }
  if (length(target_idx <- which(test %in%
    c(
      "browne.residuals.nt.model", "browne.residual.nt.model",
      "rls", "browne.rls", "nt.rls", "nt-rls", "ntrls"
    ))) > 0L) {
    test[target_idx] <- "browne.residual.nt.model"
  }

  if (length(target_idx <- which(vapply(
    test, lav_test_fmg_is_preset, logical(1L)
  ))) > 0L) {
    test[target_idx] <- vapply(
      test[target_idx],
      lav_test_fmg_resolve_preset,
      character(1L)
    )
  }


  # check?
  if (check) {
    fmg_idx <- which(vapply(test, lav_test_fmg_is_fmg, logical(1L)))

    # report unknown values
    bad_idx <- which(!test %in% c(
      "standard", "none", "default",
      "satorra.bentler",
      "yuan.bentler",
      "yuan.bentler.mplus",
      "yuan.chan",
      "mean.adjusted",
      "mean.var.adjusted",
      "scaled.shifted",
      "bollen.stine",
      "browne.residual.nt",
      "browne.residual.nt.model",
      "browne.residual.adf",
      "browne.residual.adf.model"
    ))
    bad_idx <- setdiff(bad_idx, fmg_idx)
    if (length(bad_idx) > 0L) {
      lav_msg_stop(sprintf(
        ngettext(
          length(test[bad_idx]),
          "invalid value in %1$s argument: %2$s.",
          "invalid values in %1$s argument: %2$s."
        ),
        "test", lav_msg_view(test[bad_idx], log_sep = "none")
      ))
    }

    # if 'default' is included, length(test) must be 1
    if (length(test) > 1L && any("default" == test)) {
      lav_msg_stop(
        gettextf("if test= argument contains \"%s\", it cannot contain
                 additional elements", "default"))
    }

    # if 'none' is included, length(test) must be 1
    if (length(test) > 1L && any("none" == test)) {
      lav_msg_stop(
        gettextf("if test= argument contains \"%s\" it cannot contain
                 additional elements", "none"))
    }
  }

  # reorder: first nonscaled, then scaled
  nonscaled_idx <- which(test %in% c(
    "standard", "none", "default",
    "bollen.stine",
    "browne.residual.nt",
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  scaled_idx <- which(test %in% c(
    "satorra.bentler",
    "yuan.bentler",
    "yuan.bentler.mplus",
    "yuan.chan",
    "mean.adjusted",
    "mean.var.adjusted",
    "scaled.shifted"
  ))
  fmg_idx <- which(vapply(test, lav_test_fmg_is_fmg, logical(1L)))
  test <- c(test[nonscaled_idx], test[scaled_idx], test[fmg_idx])

  test
}

lav_model_test <- function(lavobject = NULL,
                           lavmodel = NULL,
                           lavpartable = NULL,
                           lavsamplestats = NULL,
                           lavimplied = NULL,
                           lavh1 = list(),
                           lavoptions = NULL,
                           x = NULL,
                           vcov_1 = NULL,
                           lavcache = NULL,
                           lavdata = NULL,
                           lavloglik = NULL,
                           test_ugamma_eigvals = FALSE) {
  # lavobject?
  if (!is.null(lavobject)) {
    lavmodel <- lavobject@Model
    lavpartable <- lav_pt_set_cache(lavobject@ParTable, lavobject@pta)
    lavsamplestats <- lavobject@SampleStats
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavoptions <- lavobject@Options
    x <- lavobject@optim$x
    fx <- lavobject@optim[["fx"]]
    fx_group <- lavobject@optim[["fx.group"]]
    attr(fx, "fx.group") <- fx_group
    attr(x, "fx") <- fx
    vcov_1 <- lavobject@vcov$vcov
    lavcache <- lavobject@Cache
    lavdata <- lavobject@Data
    lavloglik <- lavobject@loglik
  }

  # backwards compatibility
  if (is.null(lavoptions$scaled.test)) {
    lavoptions$scaled.test <- "standard"
  }

  test <- lavoptions$test

  # "yuan.chan" is a SAM-only (sam.method = "global") rescaled test statistic;
  # it is computed in lav_sam_global_test(), not here. If it ever reaches the
  # core test machinery (eg sem(test = "yuan.chan")), drop it with a note.
  if (any(test == "yuan.chan")) {
    lav_msg_warn(gettext(
      "test = \"yuan.chan\" is only available for sam(sam.method = \"global\");
       it will be ignored here."))
    test <- test[test != "yuan.chan"]
    if (length(test) == 0L) {
      test <- "standard"
    }
    lavoptions$test <- test
  }

  test_1 <- list()

  # degrees of freedom (ignoring constraints)
  df <- lav_pt_df(lavpartable)

  # handle equality constraints (note: we ignore inequality constraints,
  # active or not!)
  # we use the rank of con.jac (even if the constraints are nonlinear)
  if (!lavmodel@cin.simple.only && nrow(lavmodel@con.jac) > 0L) {
    ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
    if (length(ceq_idx) > 0L) {
      neq <- qr(lavmodel@con.jac[ceq_idx, , drop = FALSE])$rank
      df <- df + neq
    }
  } else if (lavmodel@ceq.simple.only) {
    # needed??
    ndat <- lav_pt_ndat(lavpartable)
    npar <- max(lavpartable$free)
    df <- ndat - npar
  }

  # shortcut: return empty list if one of the conditions below is true:
  # - test == "none"
  # - df < 0
  # - estimator == "MML"
  if (test[1] == "none" || df < 0L || lavoptions$estimator == "MML") {
    test_1[[1]] <- list(
      test = test[1],
      stat = as.numeric(NA),
      stat.group = as.numeric(NA),
      df = df,
      refdistr = "unknown",
      pvalue = as.numeric(NA)
    )

    if (length(test) > 1L) {
      test_1[[2]] <- list(
        test = test[2],
        stat = as.numeric(NA),
        stat.group = as.numeric(NA),
        df = df,
        refdistr = "unknown",
        pvalue = as.numeric(NA)
      )
    }

    attr(test_1, "info") <-
      list(
        ngroups = lavdata@ngroups, group.label = lavdata@group.label,
        information = lavoptions$information,
        h1.information = lavoptions$h1.information,
        observed.information = lavoptions$observed.information
      )

    return(test_1)
  }


  ######################
  ## TEST == STANDARD ##
  ######################

  # get chisq value, per group

  # PML
  if (lavoptions$estimator == "PML" && test[1] != "none") {
    # attention!
    # if the thresholds are saturated (ie, nuisance parameters)
    # we should use the lav_pml_plrt() function.
    #
    # BUT, if the thresholds are structured (eg equality constraints)
    # then we MUST use the lav_pml_plrt2() function.
    #
    # This was not done automatically < 0.6-6
    #


    thresholds_structured <- FALSE
    # check
    th_idx <- which(lavpartable$op == "|")
    if (any(lavpartable$free[th_idx] == 0L)) {
      thresholds_structured <- TRUE
    }

    eq_idx <- which(lavpartable$op == "==")
    if (length(eq_idx) > 0L) {
      th_labels <- lavpartable$plabel[th_idx]
      eq_labels <- unique(c(
        lavpartable$lhs[eq_idx],
        lavpartable$rhs[eq_idx]
      ))
      if (any(th_labels %in% eq_labels)) {
        thresholds_structured <- TRUE
      }
    }

    # switch between lav_pml_plrt() and lav_pml_plrt2()
    if (thresholds_structured) {
      pml_plrt <- lav_pml_plrt2
    } else {
      pml_plrt <- lav_pml_plrt
    }

    pml <- pml_plrt(
      lavobject = NULL,
      lavmodel = lavmodel,
      lavdata = lavdata,
      lavoptions = lavoptions,
      x = x,
      vcov_1 = vcov_1,
      lavcache = lavcache,
      lavsamplestats = lavsamplestats,
      lavpartable = lavpartable
    )
    # get chi.group from PML, since we compare to `unrestricted' model,
    # NOT observed data
    chisq_group <- pml$PLRTH0Sat.group

    # twolevel (ML: LRT versus the h1 model; the least-squares estimators
    # fall through to the fx-based statistic below)
  } else if (lavdata@nlevels > 1L &&
    !(lavoptions$estimator %in% c("WLS", "DWLS", "ULS"))) {
    if (length(lavh1) > 0L) {
      # LRT
      chisq_group <- -2 * (lavloglik$loglik.group - lavh1$logl$loglik.group)
    } else {
      chisq_group <- rep(as.numeric(NA), lavdata@ngroups)
    }
  } else if (lavdata@missing == "ml" && any(lavdata@Mp[[1]]$coverage == 0)) {
    if (length(lavh1) > 0L) {
      chisq_group <- -2 * (lavloglik$loglik.group - lavh1$logl$loglik.group)
    } else {
      chisq_group <- rep(as.numeric(NA), lavdata@ngroups)
    }
  } else if (lavoptions$estimator %in% c("IV")) {
    # no 'standard' chi-square statistic
    chisq_group <- rep(as.numeric(NA), lavdata@ngroups)
  } else {
    # get fx.group
    fx <- attr(x, "fx")
    fx_group <- attr(fx, "fx.group")

    # always compute `standard' test statistic
    ## FIXME: the NFAC is now implicit in the computation of fx...
    nfac <- 2 * unlist(lavsamplestats@nobs)
    if (lavoptions$estimator == "ML" && lavoptions$likelihood == "wishart") {
      # first divide by two
      nfac <- nfac / 2
      nfac <- nfac - 1
      nfac <- nfac * 2
    } else if (lavoptions$estimator == "DLS") {
      nfac <- nfac / 2
      nfac <- nfac - 1
      nfac <- nfac * 2
    }

    chisq_group <- fx_group * nfac
  }

  # check for negative values
  chisq_group[chisq_group < 0] <- 0.0

  # global test statistic
  chisq <- sum(chisq_group)

  # reference distribution: always chi-square, except for the
  # non-robust version of ULS and PML
  if (lavoptions$estimator %in% c("ULS", "DWLS", "PML")) {
    refdistr <- "unknown"
    pvalue <- as.numeric(NA)
  } else {
    refdistr <- "chisq"

    # pvalue  ### FIXME: what if df=0? NA? or 1? or 0?
    # this is not trivial, since
    # 1 - pchisq(0, df=0) = 1
    # but
    # 1 - pchisq(0.00000000001, df=0) = 0
    # and
    # 1 - pchisq(0, df=0, ncp=0) = 0
    #
    # This is due to different definitions of limits (from the left,
    # or from the right)
    #
    # From 0.5-17 onwards, we will use NA if df=0, to be consistent
    if (df == 0) {
      pvalue <- as.numeric(NA)
    } else {
      pvalue <- 1 - pchisq(chisq, df)
    }
  }

  test_1[["standard"]] <- list(
    test = "standard",
    stat = chisq,
    stat.group = chisq_group,
    df = df,
    refdistr = refdistr,
    pvalue = pvalue
  )

  if (length(test) == 1L && test == "standard") {
    # we are done
    attr(test_1, "info") <-
      list(
        ngroups = lavdata@ngroups, group.label = lavdata@group.label,
        information = lavoptions$information,
        h1.information = lavoptions$h1.information,
        observed.information = lavoptions$observed.information
      )
    return(test_1)
  } else {
    # strip 'standard' from test list
    if (length(test) > 1L) {
      standard_idx <- which(test == "standard")
      if (length(standard_idx) > 0L) {
        test <- test[-standard_idx]
      }
    }
  }

  ######################
  ## additional tests ## # new in 0.6-5
  ######################

  # first: sort tests, so that the unscaled tests come first, and then the
  # scaled (otherwise, we don't find the base test)
  unscaled_idx <- which(test %in% c(
      "browne.residual.adf",
      "browne.residual.adf.model",
      "browne.residual.nt",
      "browne.residual.nt.model"
    ))
  if (length(unscaled_idx) > 0L) {
    test <- c(test[unscaled_idx], test[-unscaled_idx])
  }

  for (this_test in test) {
    if (lav_test_fmg_is_fmg(this_test)) {
      parsed <- lav_test_fmg_parse(this_test)
      chisq <- lav_test_fmg_resolve_chisq(parsed, lavoptions = lavoptions)
      unscaled_test <- test_1[[1]]
      if (chisq == "rls") {
        unscaled_test <- lav_test_fmg_browne_nt_model(
          lavobject = lavobject,
          lavdata = lavdata,
          lavsamplestats = lavsamplestats,
          lavmodel = lavmodel,
          lavpartable = lavpartable,
          lavoptions = lavoptions,
          lavh1 = lavh1,
          lavimplied = lavimplied
        )
        test_1[[unscaled_test$test]] <- unscaled_test
      }

      test_1[[this_test]] <- lav_test_fmg(
        lavobject = lavobject,
        lavsamplestats = lavsamplestats,
        lavmodel = lavmodel,
        lavdata = lavdata,
        lavoptions = lavoptions,
        lavimplied = lavimplied,
        test_unscaled = test_1[[1]],
        test_chisq = unscaled_test,
        e_inv = attr(vcov_1, "E.inv"),
        delta = attr(vcov_1, "Delta"),
        wls_v = attr(vcov_1, "WLS.V"),
        gamma = attr(vcov_1, "Gamma"),
        test = this_test
      )
    } else if (lavoptions$estimator == "PML") {
      if (this_test == "mean.var.adjusted") {
        label <- "mean+var adjusted correction (PML)"
        test_1[[this_test]] <-
          list(
            test = this_test,
            stat = pml$stat,
            stat.group = test_1[[1]]$stat.group * pml$scaling.factor,
            df = pml$df,
            pvalue = pml$p.value,
            scaling.factor = 1 / pml$scaling.factor,
            label = label,
            shift.parameter = as.numeric(NA),
            trace.UGamma = as.numeric(NA),
            trace.UGamma4 = as.numeric(NA),
            trace.UGamma2 = as.numeric(NA),
            UGamma.eigenvalues = as.numeric(NA)
          )
      } else {
        lav_msg_warn(gettextf("test option %s not available for estimator PML",
                              this_test))
      }
    } else if (this_test %in% c(
      "browne.residual.adf",
      "browne.residual.adf.model",
      "browne.residual.nt",
      "browne.residual.nt.model"
    )) {
      adf <- TRUE
      if (this_test %in% c(
        "browne.residual.nt",
        "browne.residual.nt.model"
      )) {
        adf <- FALSE
      }
      model_based <- FALSE
      if (this_test %in% c(
        "browne.residual.adf.model",
        "browne.residual.nt.model"
      )) {
        model_based <- TRUE
      }

      out <- lav_test_browne(
        lavobject = NULL,
        lavdata = lavdata,
        lavsamplestats = lavsamplestats,
        lavmodel = lavmodel,
        lavpartable = lavpartable,
        lavoptions = lavoptions,
        lavh1 = lavh1,
        lavimplied = lavimplied,
        adf = adf,
        model_based = model_based
      )
      test_1[[this_test]] <- out
    } else if (this_test %in% c(
      "satorra.bentler",
      "mean.var.adjusted",
      "scaled.shifted"
    )) {
      # which test statistic shall we scale?
      unscaled_test <- test_1[[1]]
      if (lavoptions$scaled.test != "standard") {
        idx <- which(names(test_1) == lavoptions$scaled.test)
        if (length(idx) > 0L) {
          unscaled_test <- test_1[[idx[1]]]
        } else {
          lav_msg_warn(gettextf(
            "scaled.test [%1$s] not found among available (non scaled) tests:
            %2$s. Using standard test instead.",
            lavoptions$scaled.test, lav_msg_view(test)))
        }
      }

      out <- lav_test_sb(
        lavobject = NULL,
        lavsamplestats = lavsamplestats,
        lavmodel = lavmodel,
        lavimplied = lavimplied,
        lavdata = lavdata,
        lavoptions = lavoptions,
        test_unscaled = unscaled_test,
        e_inv = attr(vcov_1, "E.inv"),
        delta = attr(vcov_1, "Delta"),
        wls_v = attr(vcov_1, "WLS.V"),
        m_gamma = attr(vcov_1, "Gamma"),
        test = this_test,
        method = "original", # since 0.6-13
        return_ugamma = FALSE
      )
      test_1[[this_test]] <- out[[this_test]]
    } else if (this_test %in% c(
      "yuan.bentler",
      "yuan.bentler.mplus"
    )) {
      # which test statistic shall we scale?
      unscaled_test <- test_1[[1]]
      if (lavoptions$scaled.test != "standard") {
        idx <- which(names(test_1) == lavoptions$scaled.test)
        if (length(idx) > 0L) {
          unscaled_test <- test_1[[idx[1]]]
        } else {
          lav_msg_warn(gettextf(
            "scaled.test [%1$s] not found among available (non scaled) tests:
            %2$s. Using standard test instead.",
            lavoptions$scaled.test, lav_msg_view(test)))
          }
      }

      out <- lav_test_yb(
        lavobject = NULL,
        lavsamplestats = lavsamplestats,
        lavmodel = lavmodel,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        test_unscaled = unscaled_test,
        e_inv = attr(vcov_1, "E.inv"),
        b0_group = attr(vcov_1, "B0.group"),
        test = this_test,
        mimic = lavoptions$mimic,
        # method         = "default",
        return_ugamma = FALSE
      )
      test_1[[this_test]] <- out[[this_test]]
    } else if (this_test == "bollen.stine") {
      # check if we have bootstrap lavdata
      boot_test <- attr(vcov_1, "BOOT.TEST")
      if (is.null(boot_test)) {
        if (!is.null(lavoptions$bootstrap)) {
          r <- as.integer(lavoptions$bootstrap$R)
        } else {
          r <- 1000L
        }
        boot_type <- "bollen.stine"
        boot_test <-
          lav_bootstrap_internal(
            object = NULL,
            lavmodel = lavmodel,
            lavsamplestats = lavsamplestats,
            lavpartable = lavpartable,
            lavoptions = lavoptions,
            lavdata = lavdata,
            r = r,
            type = boot_type,
            fun = "test"
          )

        # new in 0.6-12: always warn for failed and nonadmissible
        error_idx <- attr(boot_test, "error.idx")
        nfailed <- length(attr(boot_test, "error.idx")) # zero if NULL
        if (nfailed > 0L) {
          lav_msg_warn(gettextf(
            "%d bootstrap runs failed or did not converge.", nfailed
          ))
        }

        notok <- length(attr(boot_test, "nonadmissible")) # zero if NULL
        if (notok > 0L) {
          lav_msg_warn(gettextf(
            "%d bootstrap runs resulted in nonadmissible solutions.", notok
          ))
        }

        if (length(error_idx) > 0L) {
          # new in 0.6-13: we must still remove them!
          boot_test <- boot_test[-error_idx, , drop = FALSE]
          # this also drops the attributes
        }

        boot_test <- drop(boot_test)
      }

      # bootstrap p-value
      boot_larger <- sum(boot_test > chisq)
      boot_length <- length(boot_test)
      pvalue_boot <- boot_larger / boot_length

      test_1[[this_test]] <- list(
        test = this_test,
        stat = chisq,
        stat.group = chisq_group,
        df = df,
        pvalue = pvalue_boot,
        refdistr = "bootstrap",
        boot.T = boot_test,
        boot.larger = boot_larger,
        boot.length = boot_length
      )
    }
  } # additional tests

  # add additional information as an attribute, needed for independent
  # printing
  attr(test_1, "info") <-
    list(
      ngroups = lavdata@ngroups, group.label = lavdata@group.label,
      information = lavoptions$information,
      h1.information = lavoptions$h1.information,
      observed.information = lavoptions$observed.information
    )

  test_1
}

lav_update_test_custom_h1 <- function(lav_obj_h0, lav_obj_h1) {
  stopifnot(inherits(lav_obj_h0, "lavaan"))
  stopifnot(inherits(lav_obj_h1, "lavaan"))

  ## this breaks if object not nested in (df >=) h1, so check df
  stopifnot(lav_obj_h0@test[[1]]$df >= lav_obj_h1@test[[1]]$df)

  ## remove any other (potentially hidden) h1 model from BOTH objects
  lav_obj_h0@external$h1.model <- NULL
  lav_obj_h1@external$h1.model <- NULL

  ## save old @test slot as template
  ## (so the @test[[1]]$df don't change while looping over tests to update)
  new_test <- lav_obj_h0@test

  ## assemble a call to lavTestLRT()
  lrt_call_template <- list(quote(lavTestLRT), object = quote(lav_obj_h0),
                          quote(lav_obj_h1)) # in ...

  ## can only update tests available in both objects
  test_names0 <- names(lav_obj_h0@test)
  test_names1 <- names(lav_obj_h1@test)
  test_names <- intersect(test_names0, test_names1)

  ## loop over those tests
  for (tn in test_names) {
    lrt_call <- lrt_call_template
    ## conditional arguments:
    if (tn == "standard") {
      lrt_call$method <- "standard"
    } else if (tn %in% c("scaled.shifted", "mean.var.adjusted")) {
      if (lav_obj_h0@Options$estimator == "PML") {
        lrt_call$method <- "mean.var.adjusted.PLRT"
      } else {
        lrt_call$method <- "satorra.2000"
      }
      lrt_call$scaled.shifted <- tn == "scaled.shifted"
    } else if (tn %in% c("satorra.bentler",
                         "yuan.bentler", "yuan.bentler.mplus")) {
      lrt_call$test <- tn
    } else if (grepl(pattern = "browne", x = tn)) {
      lrt_call$type <- tn
    } else {
      #TODO?
      #   - if (tn %in% c("bootstrap", "bollen.stine")) next
      #   - any other possibilities in @test?
    }

    ## get new test
    if (lav_obj_h0@test[[1]]$df == lav_obj_h1@test[[1]]$df) {
      ## suppress warning about == df
      anova_1 <- suppressWarnings(eval(as.call(lrt_call)))
    } else {
      ## maybe some other informative warning would be important to see
      anova_1 <- eval(as.call(lrt_call))
    }


    ## replace old @test[[tn]] values
    new_test[[tn]]$stat.group <- NULL # avoid wrong stats in summary() header?
    new_test[[tn]]$stat   <- anova_1["lav_obj_h0", "Chisq diff"]
    new_test[[tn]]$df     <- anova_1["lav_obj_h0", "Df diff"]
    new_test[[tn]]$pvalue <- anova_1["lav_obj_h0", "Pr(>Chisq)"]
    if (!is.null(new_test[[tn]]$scaling.factor)) {
      new_test[[tn]]$scaling.factor <- attr(anova_1, "scale")[2] #1st row is NA
    }
    if (!is.null(new_test[[tn]]$shift.parameter)) {
      new_test[[tn]]$shift.parameter <- attr(anova_1, "shift")[2] #1st row is NA
    } else {
      ## unless scaled.shifted, RMSEA is calculated from $standard$stat and
      ## df == sum($trace.UGamma).  Reverse-engineer from $scaling.factor:
      new_test[[tn]]$trace.UGamma <- new_test[[tn]]$df *
                                                   new_test[[tn]]$scaling.factor
    }
    ## should not be necessary to replace $trace.UGamma2
    ## nor to replace $scaling.factor.h0/h1
  } # end loop over tests

  ## assign updated @test slot and return
  lav_obj_h0@test <- new_test
  lav_obj_h0
}
