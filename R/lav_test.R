# chi-square test statistic:
# comparing the current model versus the saturated/unrestricted model
# TDJ 9 April 2024: Add a (hidden) function to update the @test slot
#                   when the user provides a custom h1 model.  Called by
#                   lav_object_summary and lav_fit_measures(_check_baseline)

lavTest <- function(lavobject, test = "standard",
                    scaled.test = "standard",
                    output = "list", drop.list.single = TRUE) {
  # check output
  output.valid <- c("list", "text")
  if (!any(output == output.valid)) {
      lav_msg_stop(gettextf(
        "%1$s argument must be either %2$s",
        "output", lav_msg_view(output.valid, "or")
      ))
  }
  # extract 'test' slot
  TEST <- lavobject@test

  # which test?
  if (!missing(test)) {
    # check 'test'
    if (!is.character(test)) {
      lav_msg_stop(
        gettextf("%s should be a character string.", "test"))
    } else {
      test <- lav_test_rename(test, check = TRUE)
    }

    # check scaled.test
    if (!missing(scaled.test)) {
      if (!is.character(scaled.test)) {
        lav_msg_stop(
          gettextf("%s should be a character string.", "scaled.test"))
      } else {
        scaled.test <- lav_test_rename(scaled.test, check = TRUE)
      }

      # merge
      test <- unique(c(test, scaled.test))

      # but "standard" must always be first
      standard.idx <- which(test == "standard")
      if (length(standard.idx) > 0L && standard.idx != 1L) {
        test <- c("standard", test[-standard.idx])
      }
    }

    if (test[1] == "none") {
      return(list())
    } else if (any(test %in% c("bootstrap", "bollen.stine"))) {
      lav_msg_stop(gettext(
      "please use bootstrapLavaan() to obtain a bootstrap based test statistic."
      ))
    }

    # check if we already have it:
    if (all(test %in% names(TEST))) {
      info.attr <- attr(TEST, "info")
      test.idx <- which(names(TEST) %in% test)
      TEST <- TEST[test.idx]
      attr(TEST, "info") <- info.attr
    } else {
      # redo ALL of them, even if already have some in TEST
      # later, we will allow to also change the options (like information)
      # and this should be reflected in the 'info' attribute

      # fill-in test in Options slot
      lavobject@Options$test <- test

      # fill-in scaled.test in Options slot
      lavobject@Options$scaled.test <- scaled.test

      # get requested test statistics
      TEST <- lav_model_test(lavobject = lavobject)
    }
  }

  if (output == "list") {
    # remove 'info' attribute
    attr(TEST, "info") <- NULL

    # select only those that were requested (eg remove standard)
    test.idx <- which(names(TEST) %in% test)
    TEST <- TEST[test.idx]

    # if only 1 test, drop outer list
    if (length(TEST) == 1L && drop.list.single) {
      TEST <- TEST[[1]]
    }

    return(TEST)
  } else {
    lav_test_print(TEST)
  }

  invisible(TEST)
}

# allow for 'flexible' names for the test statistics
# 0.6-13: if multiple names, order them in such a way
#         that the 'scaled' variants appear after the others
lav_test_rename <- function(test, check = FALSE) {
  test <- tolower(test)

  if (length(target.idx <- which(test %in%
    c("standard", "chisq", "chi", "chi-square", "chi.square"))) > 0L) {
    test[target.idx] <- "standard"
  }
  if (length(target.idx <- which(test %in%
    c(
      "satorra", "sb", "satorra.bentler", "satorra-bentler",
      "m.adjusted", "m", "mean.adjusted", "mean-adjusted"
    ))) > 0L) {
    test[target.idx] <- "satorra.bentler"
  }
  if (length(target.idx <- which(test %in%
    c("yuan", "yb", "yuan.bentler", "yuan-bentler"))) > 0L) {
    test[target.idx] <- "yuan.bentler"
  }
  if (length(target.idx <- which(test %in%
    c(
      "yuan.bentler.mplus", "yuan-bentler.mplus",
      "yuan-bentler-mplus"
    ))) > 0L) {
    test[target.idx] <- "yuan.bentler.mplus"
  }
  if (length(target.idx <- which(test %in%
    c(
      "mean.var.adjusted", "mean-var-adjusted", "mv", "second.order",
      "satterthwaite", "mv.adjusted"
    ))) > 0L) {
    test[target.idx] <- "mean.var.adjusted"
  }
  if (length(target.idx <- which(test %in%
    c(
      "mplus6", "scale.shift", "scaled.shifted",
      "scaled-shifted"
    ))) > 0L) {
    test[target.idx] <- "scaled.shifted"
  }
  if (length(target.idx <- which(test %in%
    c("bootstrap", "boot", "bollen.stine", "bollen-stine"))) > 0L) {
    test[target.idx] <- "bollen.stine"
  }
  if (length(target.idx <- which(test %in%
    c(
      "browne", "residual", "residuals", "browne.residual",
      "browne.residuals", "residual-based", "residual.based",
      "browne.residuals.adf", "browne.residual.adf"
    ))) > 0L) {
    test[target.idx] <- "browne.residual.adf"
  }
  if (length(target.idx <- which(test %in%
    c("browne.residuals.nt", "browne.residual.nt"))) > 0L) {
    test[target.idx] <- "browne.residual.nt"
  }
  if (length(target.idx <- which(test %in%
    c("browne.residual.adf.model"))) > 0L) {
    test[target.idx] <- "browne.residual.adf.model"
  }
  if (length(target.idx <- which(test %in%
    c(
      "browne.residuals.nt.model", "browne.residual.nt.model",
      "rls", "browne.rls", "nt.rls", "nt-rls", "ntrls"
    ))) > 0L) {
    test[target.idx] <- "browne.residual.nt.model"
  }


  # check?
  if (check) {
    # report unknown values
    bad.idx <- which(!test %in% c(
      "standard", "none", "default",
      "satorra.bentler",
      "yuan.bentler",
      "yuan.bentler.mplus",
      "mean.adjusted",
      "mean.var.adjusted",
      "scaled.shifted",
      "bollen.stine",
      "browne.residual.nt",
      "browne.residual.nt.model",
      "browne.residual.adf",
      "browne.residual.adf.model"
    ))
    if (length(bad.idx) > 0L) {
      lav_msg_stop(sprintf(
        ngettext(
          length(test[bad.idx]),
          "invalid value in %1$s argument: %2$s.",
          "invalid values in %1$s argument: %2$s."
        ),
        "test", lav_msg_view(test[bad.idx], log.sep = "none")
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
  nonscaled.idx <- which(test %in% c(
    "standard", "none", "default",
    "bollen.stine",
    "browne.residual.nt",
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  scaled.idx <- which(test %in% c(
    "satorra.bentler",
    "yuan.bentler",
    "yuan.bentler.mplus",
    "mean.adjusted",
    "mean.var.adjusted",
    "scaled.shifted"
  ))
  test <- c(test[nonscaled.idx], test[scaled.idx])

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
                           VCOV = NULL,
                           lavcache = NULL,
                           lavdata = NULL,
                           lavloglik = NULL,
                           test.UGamma.eigvals = FALSE) {
  # lavobject?
  if (!is.null(lavobject)) {
    lavmodel <- lavobject@Model
    lavpartable <- lav_partable_set_cache(lavobject@ParTable, lavobject@pta)
    lavsamplestats <- lavobject@SampleStats
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavoptions <- lavobject@Options
    x <- lavobject@optim$x
    fx <- lavobject@optim[["fx"]]
    fx.group <- lavobject@optim[["fx.group"]]
    attr(fx, "fx.group") <- fx.group
    attr(x, "fx") <- fx
    VCOV <- lavobject@vcov$vcov
    lavcache <- lavobject@Cache
    lavdata <- lavobject@Data
    lavloglik <- lavobject@loglik
  }

  # backwards compatibility
  if (is.null(lavoptions$scaled.test)) {
    lavoptions$scaled.test <- "standard"
  }

  test <- test.orig <- lavoptions$test

  TEST <- list()

  # degrees of freedom (ignoring constraints)
  df <- lav_partable_df(lavpartable)

  # handle equality constraints (note: we ignore inequality constraints,
  # active or not!)
  # we use the rank of con.jac (even if the constraints are nonlinear)
  if (nrow(lavmodel@con.jac) > 0L) {
    ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
    if (length(ceq.idx) > 0L) {
      neq <- qr(lavmodel@con.jac[ceq.idx, , drop = FALSE])$rank
      df <- df + neq
    }
  } else if (lavmodel@ceq.simple.only) {
    # needed??
    ndat <- lav_partable_ndat(lavpartable)
    npar <- max(lavpartable$free)
    df <- ndat - npar
  }

  # shortcut: return empty list if one of the conditions below is true:
  # - test == "none"
  # - df < 0
  # - estimator == "MML"
  if (test[1] == "none" || df < 0L || lavoptions$estimator == "MML") {
    TEST[[1]] <- list(
      test = test[1],
      stat = as.numeric(NA),
      stat.group = as.numeric(NA),
      df = df,
      refdistr = "unknown",
      pvalue = as.numeric(NA)
    )

    if (length(test) > 1L) {
      TEST[[2]] <- list(
        test = test[2],
        stat = as.numeric(NA),
        stat.group = as.numeric(NA),
        df = df,
        refdistr = "unknown",
        pvalue = as.numeric(NA)
      )
    }

    attr(TEST, "info") <-
      list(
        ngroups = lavdata@ngroups, group.label = lavdata@group.label,
        information = lavoptions$information,
        h1.information = lavoptions$h1.information,
        observed.information = lavoptions$observed.information
      )

    return(TEST)
  }


  ######################
  ## TEST == STANDARD ##
  ######################

  # get chisq value, per group

  # PML
  if (lavoptions$estimator == "PML" && test[1] != "none") {
    # attention!
    # if the thresholds are saturated (ie, nuisance parameters)
    # we should use the ctr_pml_plrt() function.
    #
    # BUT, if the thresholds are structured (eg equality constraints)
    # then we MUST use the ctr_pml_plrt2() function.
    #
    # This was not done automatically < 0.6-6
    #


    thresholds.structured <- FALSE
    # check
    th.idx <- which(lavpartable$op == "|")
    if (any(lavpartable$free[th.idx] == 0L)) {
      thresholds.structured <- TRUE
    }

    eq.idx <- which(lavpartable$op == "==")
    if (length(eq.idx) > 0L) {
      th.labels <- lavpartable$plabel[th.idx]
      eq.labels <- unique(c(
        lavpartable$lhs[eq.idx],
        lavpartable$rhs[eq.idx]
      ))
      if (any(th.labels %in% eq.labels)) {
        thresholds.structured <- TRUE
      }
    }

    # switch between ctr_pml_plrt() and ctr_pml_plrt2()
    if (thresholds.structured) {
      pml_plrt <- ctr_pml_plrt2
    } else {
      pml_plrt <- ctr_pml_plrt
    }

    PML <- pml_plrt(
      lavobject = NULL,
      lavmodel = lavmodel,
      lavdata = lavdata,
      lavoptions = lavoptions,
      x = x,
      VCOV = VCOV,
      lavcache = lavcache,
      lavsamplestats = lavsamplestats,
      lavpartable = lavpartable
    )
    # get chi.group from PML, since we compare to `unrestricted' model,
    # NOT observed data
    chisq.group <- PML$PLRTH0Sat.group

    # twolevel
  } else if (lavdata@nlevels > 1L) {
    if (length(lavh1) > 0L) {
      # LRT
      chisq.group <- -2 * (lavloglik$loglik.group - lavh1$logl$loglik.group)
    } else {
      chisq.group <- rep(as.numeric(NA), lavdata@ngroups)
    }
  } else {
    # get fx.group
    fx <- attr(x, "fx")
    fx.group <- attr(fx, "fx.group")

    # always compute `standard' test statistic
    ## FIXME: the NFAC is now implicit in the computation of fx...
    NFAC <- 2 * unlist(lavsamplestats@nobs)
    if (lavoptions$estimator == "ML" && lavoptions$likelihood == "wishart") {
      # first divide by two
      NFAC <- NFAC / 2
      NFAC <- NFAC - 1
      NFAC <- NFAC * 2
    } else if (lavoptions$estimator == "DLS") {
      NFAC <- NFAC / 2
      NFAC <- NFAC - 1
      NFAC <- NFAC * 2
    }

    chisq.group <- fx.group * NFAC
  }

  # check for negative values
  chisq.group[chisq.group < 0] <- 0.0

  # global test statistic
  chisq <- sum(chisq.group)

  # reference distribution: always chi-square, except for the
  # non-robust version of ULS and PML
  if (lavoptions$estimator == "ULS" || lavoptions$estimator == "PML") {
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

  TEST[["standard"]] <- list(
    test = "standard",
    stat = chisq,
    stat.group = chisq.group,
    df = df,
    refdistr = refdistr,
    pvalue = pvalue
  )

  if (length(test) == 1L && test == "standard") {
    # we are done
    attr(TEST, "info") <-
      list(
        ngroups = lavdata@ngroups, group.label = lavdata@group.label,
        information = lavoptions$information,
        h1.information = lavoptions$h1.information,
        observed.information = lavoptions$observed.information
      )
    return(TEST)
  } else {
    # strip 'standard' from test list
    if (length(test) > 1L) {
      standard.idx <- which(test == "standard")
      if (length(standard.idx) > 0L) {
        test <- test[-standard.idx]
      }
    }
  }




  ######################
  ## additional tests ## # new in 0.6-5
  ######################

  for (this.test in test) {
    if (lavoptions$estimator == "PML") {
      if (this.test == "mean.var.adjusted") {
        LABEL <- "mean+var adjusted correction (PML)"
        TEST[[this.test]] <-
          list(
            test = this.test,
            stat = PML$stat,
            stat.group = TEST[[1]]$stat.group * PML$scaling.factor,
            df = PML$df,
            pvalue = PML$p.value,
            scaling.factor = 1 / PML$scaling.factor,
            label = LABEL,
            shift.parameter = as.numeric(NA),
            trace.UGamma = as.numeric(NA),
            trace.UGamma4 = as.numeric(NA),
            trace.UGamma2 = as.numeric(NA),
            UGamma.eigenvalues = as.numeric(NA)
          )
      } else {
        lav_msg_warn(gettextf("test option %s not available for estimator PML",
                              this.test))
      }
    } else if (this.test %in% c(
      "browne.residual.adf",
      "browne.residual.adf.model",
      "browne.residual.nt",
      "browne.residual.nt.model"
    )) {
      ADF <- TRUE
      if (this.test %in% c(
        "browne.residual.nt",
        "browne.residual.nt.model"
      )) {
        ADF <- FALSE
      }
      model.based <- FALSE
      if (this.test %in% c(
        "browne.residual.adf.model",
        "browne.residual.nt.model"
      )) {
        model.based <- TRUE
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
        ADF = ADF,
        model.based = model.based
      )
      TEST[[this.test]] <- out
    } else if (this.test %in% c(
      "satorra.bentler",
      "mean.var.adjusted",
      "scaled.shifted"
    )) {
      # which test statistic shall we scale?
      unscaled.TEST <- TEST[[1]]
      if (lavoptions$scaled.test != "standard") {
        idx <- which(test.orig == lavoptions$scaled.test)
        if (length(idx) > 0L) {
          unscaled.TEST <- TEST[[idx[1]]]
        } else {
          lav_msg_warn(gettextf(
            "scaled.test [%1$s] not found among available (non scaled) tests:
            %2$s. Using standard test instead.",
            lavoptions$scaled.test, lav_msg_view(test)))
        }
      }

      out <- lav_test_satorra_bentler(
        lavobject = NULL,
        lavsamplestats = lavsamplestats,
        lavmodel = lavmodel,
        lavimplied = lavimplied,
        lavdata = lavdata,
        lavoptions = lavoptions,
        TEST.unscaled = unscaled.TEST,
        E.inv = attr(VCOV, "E.inv"),
        Delta = attr(VCOV, "Delta"),
        WLS.V = attr(VCOV, "WLS.V"),
        Gamma = attr(VCOV, "Gamma"),
        test = this.test,
        mimic = lavoptions$mimic,
        method = "original", # since 0.6-13
        return.ugamma = FALSE
      )
      TEST[[this.test]] <- out[[this.test]]
    } else if (this.test %in% c(
      "yuan.bentler",
      "yuan.bentler.mplus"
    )) {
      # which test statistic shall we scale?
      unscaled.TEST <- TEST[[1]]
      if (lavoptions$scaled.test != "standard") {
        idx <- which(test.orig == lavoptions$scaled.test)
        if (length(idx) > 0L) {
          unscaled.TEST <- TEST[[idx[1]]]
        } else {
          lav_msg_warn(gettextf(
            "scaled.test [%1$s] not found among available (non scaled) tests:
            %2$s. Using standard test instead.",
            lavoptions$scaled.test, lav_msg_view(test)))
          }
      }

      out <- lav_test_yuan_bentler(
        lavobject = NULL,
        lavsamplestats = lavsamplestats,
        lavmodel = lavmodel,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        TEST.unscaled = unscaled.TEST,
        E.inv = attr(VCOV, "E.inv"),
        B0.group = attr(VCOV, "B0.group"),
        test = this.test,
        mimic = lavoptions$mimic,
        # method         = "default",
        return.ugamma = FALSE
      )
      TEST[[this.test]] <- out[[this.test]]
    } else if (this.test == "bollen.stine") {
      # check if we have bootstrap lavdata
      BOOT.TEST <- attr(VCOV, "BOOT.TEST")
      if (is.null(BOOT.TEST)) {
        if (!is.null(lavoptions$bootstrap)) {
          R <- lavoptions$bootstrap
        } else {
          R <- 1000L
        }
        boot.type <- "bollen.stine"
        BOOT.TEST <-
          lav_bootstrap_internal(
            object = NULL,
            lavmodel. = lavmodel,
            lavsamplestats. = lavsamplestats,
            lavpartable. = lavpartable,
            lavoptions. = lavoptions,
            lavdata. = lavdata,
            R = R,
            verbose = lavoptions$verbose,
            type = boot.type,
            FUN = "test"
          )

        # new in 0.6-12: always warn for failed and nonadmissible
        error.idx <- attr(BOOT.TEST, "error.idx")
        nfailed <- length(attr(BOOT.TEST, "error.idx")) # zero if NULL
        if (nfailed > 0L && lavoptions$warn) {
          lav_msg_warn(gettextf(
            "%d bootstrap runs failed or did not converge.", nfailed
          ))
        }

        notok <- length(attr(BOOT.TEST, "nonadmissible")) # zero if NULL
        if (notok > 0L && lavoptions$warn) {
          lav_msg_warn(gettextf(
            "%d bootstrap runs resulted in nonadmissible solutions.", notok
          ))
        }

        if (length(error.idx) > 0L) {
          # new in 0.6-13: we must still remove them!
          BOOT.TEST <- BOOT.TEST[-error.idx, , drop = FALSE]
          # this also drops the attributes
        }

        BOOT.TEST <- drop(BOOT.TEST)
      }

      # bootstrap p-value
      boot.larger <- sum(BOOT.TEST > chisq)
      boot.length <- length(BOOT.TEST)
      pvalue.boot <- boot.larger / boot.length

      TEST[[this.test]] <- list(
        test = this.test,
        stat = chisq,
        stat.group = chisq.group,
        df = df,
        pvalue = pvalue.boot,
        refdistr = "bootstrap",
        boot.T = BOOT.TEST,
        boot.larger = boot.larger,
        boot.length = boot.length
      )
    }
  } # additional tests

  # add additional information as an attribute, needed for independent
  # printing
  attr(TEST, "info") <-
    list(
      ngroups = lavdata@ngroups, group.label = lavdata@group.label,
      information = lavoptions$information,
      h1.information = lavoptions$h1.information,
      observed.information = lavoptions$observed.information
    )

  TEST
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
  newTEST <- lav_obj_h0@test

  ## assemble a call to lavTestLRT()
  lrtCallTemplate <- list(quote(lavTestLRT), object = quote(lav_obj_h1),
                          quote(lav_obj_h0)) # in ...

  ## can only update tests available in both objects
  testNames0 <- names(lav_obj_h0@test)
  testNames1 <- names(lav_obj_h1@test)
  testNames <- intersect(testNames0, testNames1)

  ## loop over those tests
  for (tn in testNames) {
    lrtCall <- lrtCallTemplate
    ## conditional arguments:
    if (tn == "standard") {
      lrtCall$method <- "standard"
    } else if (tn %in% c("scaled.shifted","mean.var.adjusted")) {
      if (lav_obj_h0@Options$estimator == "PML") {
        lrtCall$method <- "mean.var.adjusted.PLRT"
      } else {
        lrtCall$method <- "satorra.2000"
      }
      lrtCall$scaled.shifted <- tn == "scaled.shifted"
    } else if (tn %in% c("satorra.bentler",
                         "yuan.bentler","yuan.bentler.mplus")) {
      lrtCall$test <- tn
    } else if (grepl(pattern = "browne", x = tn)) {
      lrtCall$type <- tn
    } else {
      #TODO?
      #   - if (tn %in% c("bootstrap", "bollen.stine")) next
      #   - any other possibilities in @test?
    }

    ## get new test
    if (lav_obj_h0@test[[1]]$df == lav_obj_h1@test[[1]]$df) {
      ## suppress warning about == df
      ANOVA <- suppressWarnings(eval(as.call(lrtCall)))
    } else {
      ## maybe some other informative warning would be important to see
      ANOVA <- eval(as.call(lrtCall))
    }


    ## replace old @test[[tn]] values
    newTEST[[tn]]$stat.group <- NULL # avoid wrong stats in summary() header?
    newTEST[[tn]]$stat   <- ANOVA["lav_obj_h0" , "Chisq diff"]
    newTEST[[tn]]$df     <- ANOVA["lav_obj_h0" , "Df diff"   ]
    newTEST[[tn]]$pvalue <- ANOVA["lav_obj_h0" , "Pr(>Chisq)"]
    if (!is.null(newTEST[[tn]]$scaling.factor)) {
      newTEST[[tn]]$scaling.factor <- attr(ANOVA, "scale")[2] # first row is NA
    }
    if (!is.null(newTEST[[tn]]$shift.parameter)) {
      newTEST[[tn]]$shift.parameter <- attr(ANOVA, "shift")[2] # first row is NA
    } else {
      ## unless scaled.shifted, RMSEA is calculated from $standard$stat and
      ## df == sum($trace.UGamma).  Reverse-engineer from $scaling factor:
      newTEST[[tn]]$trace.UGamma <- newTEST[[tn]]$df * newTEST[[tn]]$scaling.factor
    }
    ## should not be necessary to replace $trace.UGamma2
    ## nor to replace $scaling.factor.h0/h1
  } # end loop over tests

  ## assign updated @test slot and return
  lav_obj_h0@test <- newTEST
  lav_obj_h0
}


