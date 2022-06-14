lav_model_test <- function(lavmodel       = NULL,
                           lavpartable    = NULL,
                           lavpta         = NULL,
                           lavsamplestats = NULL,
                           lavimplied     = NULL,
                           lavh1          = list(),
                           lavoptions     = NULL,
                           x              = NULL,
                           VCOV           = NULL,
                           lavcache       = NULL,
                           lavdata        = NULL,
                           lavloglik      = NULL,
                           test.UGamma.eigvals = FALSE) {


    test <- lavoptions$test

    TEST <- list()

    # degrees of freedom (ignoring constraints)
    df <- lav_partable_df(lavpartable)

    # handle equality constraints (note: we ignore inequality constraints,
    # active or not!)
    # we use the rank of con.jac (even if the constraints are nonlinear)
    if(nrow(lavmodel@con.jac) > 0L) {
        ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
        if(length(ceq.idx) > 0L) {
            neq <- qr(lavmodel@con.jac[ceq.idx,,drop=FALSE])$rank
            df <- df + neq
        }
    } else if(lavmodel@ceq.simple.only) {
        ndat <- lav_partable_ndat(lavpartable)
        npar <- max(lavpartable$free)
        df <- ndat - npar
    }

    # shortcut: return empty list if one of the conditions below is true:
    # - test == "none"
    # - df < 0
    # - estimator == "MML"
    if(test[1] == "none" || df < 0L || lavoptions$estimator == "MML") {

        TEST[[1]] <- list(test       = test[1],
                          stat       = as.numeric(NA),
                          stat.group = as.numeric(NA),
                          df         = df,
                          refdistr   = "unknown",
                          pvalue     = as.numeric(NA))

        if(length(test) > 1L) {
            TEST[[2]] <- list(test       = test[2],
                              stat       = as.numeric(NA),
                              stat.group = as.numeric(NA),
                              df         = df,
                              refdistr   = "unknown",
                              pvalue     = as.numeric(NA))
        }

        attr(TEST, "info") <-
        list(ngroups = lavdata@ngroups, group.label = lavdata@group.label,
             information = lavoptions$information,
             h1.information = lavoptions$h1.information,
             observed.information = lavoptions$observed.information)

        return(TEST)
    }


    ######################
    ## TEST == STANDARD ##
    ######################

    # get chisq value, per group

    # PML
    if(lavoptions$estimator == "PML" && test[1] != "none") {

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
        if(any(lavpartable$free[th.idx] == 0L)) {
            thresholds.structured <- TRUE
        }

        eq.idx <- which(lavpartable$op == "==")
        if(length(eq.idx) > 0L) {
            th.labels <- lavpartable$plabel[th.idx]
            eq.labels <- unique(c(lavpartable$lhs[eq.idx],
                                  lavpartable$rhs[eq.idx]))
            if(any(th.labels %in% eq.labels)) {
                thresholds.structured <- TRUE
            }
        }

        # switch between ctr_pml_plrt() and ctr_pml_plrt2()
        if(thresholds.structured) {
            pml_plrt <- ctr_pml_plrt2
        } else {
            pml_plrt <- ctr_pml_plrt
        }

        PML <- pml_plrt(lavobject      = NULL,
                        lavmodel       = lavmodel,
                        lavdata        = lavdata,
                        lavoptions     = lavoptions,
                        lavpta         = lavpta,
                        x              = x,
                        VCOV           = VCOV,
                        lavcache       = lavcache,
                        lavsamplestats = lavsamplestats,
                        lavpartable    = lavpartable)
        # get chi.group from PML, since we compare to `unrestricted' model,
        # NOT observed data
        chisq.group <- PML$PLRTH0Sat.group

    # twolevel
    } else if(lavdata@nlevels > 1L) {

        if(length(lavh1) > 0L) {
            # LRT
            chisq.group <- -2 * (lavloglik$loglik.group - lavh1$loglik.group)
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
        if(lavoptions$estimator == "ML" && lavoptions$likelihood == "wishart") {
            # first divide by two
            NFAC <- NFAC / 2; NFAC <- NFAC - 1; NFAC <- NFAC * 2
        } else if(lavoptions$estimator == "DLS") {
            NFAC <- NFAC / 2; NFAC <- NFAC - 1; NFAC <- NFAC * 2
        }

        chisq.group <- fx.group * NFAC
    }

    # check for negative values
    chisq.group[ chisq.group < 0 ] <- 0.0

    # global test statistic
    chisq <- sum(chisq.group)

    # reference distribution: always chi-square, except for the
    # non-robust version of ULS and PML
    if(lavoptions$estimator == "ULS" || lavoptions$estimator == "PML") {
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
        if(df == 0) {
            pvalue <- as.numeric(NA)
        } else {
            pvalue <- 1 - pchisq(chisq, df)
        }
    }

    TEST[["standard"]] <- list(test       = "standard",
                               stat       = chisq,
                               stat.group = chisq.group,
                               df         = df,
                               refdistr   = refdistr,
                               pvalue     = pvalue)

    if(length(test) == 1L && test == "standard") {
        # we are done
        attr(TEST, "info") <-
        list(ngroups = lavdata@ngroups, group.label = lavdata@group.label,
             information = lavoptions$information,
             h1.information = lavoptions$h1.information,
             observed.information = lavoptions$observed.information)
        return(TEST)
    } else {
        # strip 'standard' from test list
        if(length(test) > 1L) {
            standard.idx <- which(test == "standard")
            if(length(standard.idx) > 0L) {
                test <- test[-standard.idx]
            }
        }
    }




    ######################
    ## additional tests ## # new in 0.6-5
    ######################

    for(this.test in test) {

        if(lavoptions$estimator == "PML") {
            if(this.test == "mean.var.adjusted") {
                LABEL <- "mean+var adjusted correction (PML)"
                TEST[[this.test]] <-
                    list(test                 = this.test,
                         stat                 = PML$stat,
                         stat.group           = TEST[[1]]$stat.group*PML$scaling.factor,
                         df                   = PML$df,
                         pvalue               = PML$p.value,
                         scaling.factor       = 1/PML$scaling.factor,
                         label                = LABEL,
                         shift.parameter      = as.numeric(NA),
                         trace.UGamma         = as.numeric(NA),
                         trace.UGamma4        = as.numeric(NA),
                         trace.UGamma2        = as.numeric(NA),
                         UGamma.eigenvalues   = as.numeric(NA))
            } else {
                warning("test option ", this.test,
                        " not available for estimator PML")
            }



        } else if(this.test %in% c("satorra.bentler",
                                   "mean.var.adjusted",
                                   "scaled.shifted")) {

            out <- lav_test_satorra_bentler(lavobject = NULL,
                             lavsamplestats = lavsamplestats,
                             lavmodel       = lavmodel,
                             lavimplied     = lavimplied,
                             lavdata        = lavdata,
                             lavoptions     = lavoptions,
                             TEST.unscaled  = TEST[[1]],
                             E.inv          = attr(VCOV, "E.inv"),
                             Delta          = attr(VCOV, "Delta"),
                             WLS.V          = attr(VCOV, "WLS.V"),
                             Gamma          = attr(VCOV, "Gamma"),
                             test           = this.test,
                             mimic          = lavoptions$mimic,
                             method         = "ABA",
                             return.ugamma  = FALSE)
            TEST[[this.test]] <- out[[this.test]]

        } else if(this.test %in% c("yuan.bentler",
                                   "yuan.bentler.mplus")) {

            out <- lav_test_yuan_bentler(lavobject = NULL,
                             lavsamplestats = lavsamplestats,
                             lavmodel       = lavmodel,
                             lavdata        = lavdata,
                             lavimplied     = lavimplied,
                             lavh1          = lavh1,
                             lavoptions     = lavoptions,
                             TEST.unscaled  = TEST[[1]],
                             E.inv          = attr(VCOV, "E.inv"),
                             B0.group       = attr(VCOV, "B0.group"),
                             test           = this.test,
                             mimic          = lavoptions$mimic,
                             #method         = "default",
                             return.ugamma  = FALSE)
            TEST[[this.test]] <- out[[this.test]]

        } else if(this.test == "bollen.stine") {

            # check if we have bootstrap lavdata
            BOOT.TEST <- attr(VCOV, "BOOT.TEST")
            if(is.null(BOOT.TEST)) {
                if(!is.null(lavoptions$bootstrap)) {
                    R <- lavoptions$bootstrap
                } else {
                    R <- 1000L
                }
                boot.type <- "bollen.stine"
                BOOT.TEST <-
                    lav_bootstrap_internal(object          = NULL,
                                           lavmodel.       = lavmodel,
                                           lavsamplestats. = lavsamplestats,
                                           lavpartable.    = lavpartable,
                                           lavoptions.     = lavoptions,
                                           lavdata.        = lavdata,
                                           R               = R,
                                           verbose         = lavoptions$verbose,
                                           type            = boot.type,
                                           FUN             = "test",
                                           warn            = -1L)
                BOOT.TEST <- drop(BOOT.TEST)
            }

            # bootstrap p-value
            boot.larger <- sum(BOOT.TEST > chisq)
            boot.length <- length(BOOT.TEST)
            pvalue.boot <- boot.larger/boot.length

            TEST[[this.test]] <- list(test        = this.test,
                                      stat        = chisq,
                                      stat.group  = chisq.group,
                                      df          = df,
                                      pvalue      = pvalue.boot,
                                      refdistr    = "bootstrap",
                                      boot.T      = BOOT.TEST,
                                      boot.larger = boot.larger,
                                      boot.length = boot.length)
        }

    } # additional tests

    # add additional information as an attribute, needed for independent
    # printing
    attr(TEST, "info") <-
        list(ngroups = lavdata@ngroups, group.label = lavdata@group.label,
             information = lavoptions$information,
             h1.information = lavoptions$h1.information,
             observed.information = lavoptions$observed.information)

    TEST
}


