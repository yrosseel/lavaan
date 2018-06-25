lav_model_test <- function(lavmodel       = NULL,
                           lavpartable    = NULL,
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


    mimic         <- lavoptions$mimic
    test          <- lavoptions$test
    information   <- lavoptions$information

    estimator     <- lavmodel@estimator
    meanstructure <- lavmodel@meanstructure


    TEST <- list()

    # degrees of freedom
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
    }


    if(test == "none" || df < 0L || estimator == "MML") {
        TEST[[1]] <- list(test=test,
                          stat=as.numeric(NA),
                          stat.group=as.numeric(NA),
                          df=df,
                          refdistr="unknown",
                          pvalue=as.numeric(NA))

        # just in case
        TEST[[2]] <- list(test=test,
                          stat=as.numeric(NA),
                          stat.group=as.numeric(NA),
                          df=df,
                          refdistr="unknown",
                          pvalue=as.numeric(NA))
        return(TEST)
    }

    if(lavoptions$estimator == "PML" && test != "none") {
        PML <- ctr_pml_plrt(lavobject = NULL,
                            lavmodel = lavmodel,
                            lavdata = lavdata,
                            lavoptions = lavoptions,
                            x = x, VCOV = VCOV,
                            lavcache = lavcache,
                            lavsamplestats = lavsamplestats,
                            lavpartable = lavpartable)
        # get chi.group from PML, since we compare to `unrestricted' model,
        # NOT observed data
        chisq.group <- PML$PLRTH0Sat.group
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
            NFAC <- NFAC / 2
            NFAC <- NFAC - 1
            NFAC <- NFAC * 2
        }
        chisq.group <- fx.group * NFAC
    }

    # check for negative values
    chisq.group[which(chisq.group < 0)] <- 0.0

    # global test statistic
    chisq <- sum(chisq.group)

    # reference distribution: always chi-square, except for the
    # non-robust version of ULS and PML
    if(estimator == "ULS" || estimator == "PML") {
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

    TEST[[1]] <- list(test="standard",
                      stat=chisq,
                      stat.group=chisq.group,
                      df=df,
                      refdistr=refdistr,
                      pvalue=pvalue)

    if(df == 0 && test %in% c("satorra.bentler", "yuan.bentler",
                              "yuan.bentler.mplus",
                              "mean.var.adjusted", "scaled.shifted")) {
        TEST[[2]] <- list(test=test, stat=chisq, stat.group=chisq.group,
                          df=df, refdistr=refdistr, pvalue=pvalue,
                          scaling.factor=as.numeric(NA))
        return(TEST)
    }

    # do we already know E.inv is singular?
    E.inv <- attr(VCOV, "E.inv")
    if(!is.null(E.inv) && inherits(E.inv, "try-error")) {
        TEST[[2]] <- list(test=test, stat=chisq, stat.group=chisq.group,
                          df=df, refdistr=refdistr, pvalue=pvalue,
                          scaling.factor=as.numeric(NA))
        return(TEST)
    }

    # some require meanstructure (for now)
    #if(test %in% c("satorra.bentler", "yuan.bentler") &&
    #   !lavoptions$meanstructure) {
    #    stop("test (", test, ") requires meanstructure (for now)")
    #}

    # fixed.x idx
    if(lavmodel@fixed.x && estimator == "ML" && !lavmodel@conditional.x) {
        x.idx <- lavsamplestats@x.idx
    } else {
        x.idx <- vector("list", length=lavsamplestats@ngroups)
        for(g in 1:lavsamplestats@ngroups) {
            x.idx[[g]] <- integer(0L)
        }
    }

    if(lavoptions$estimator == "PML") {
        if(test == "standard") {
            # nothing to do
        } else if(test == "mean.var.adjusted") {
            TEST[[2]] <- list(test               = test,
                              stat               = PML$stat,
                              stat.group         = TEST[[1]]$stat.group*PML$scaling.factor,
                              df                 = PML$df,
                              pvalue             = PML$p.value,
                              scaling.factor     = 1/PML$scaling.factor,
                              shift.parameter    = as.numeric(NA),
                              trace.UGamma       = as.numeric(NA),
                              trace.UGamma4      = as.numeric(NA),
                              trace.UGamma2      = as.numeric(NA),
                              UGamma.eigenvalues = as.numeric(NA))
        } else {
            warning("test option ", test, " not available for estimator PML")
        }
    } else if(test %in%
          c("satorra.bentler", "mean.var.adjusted", "scaled.shifted") &&
          df > 0 && lavoptions$estimator != "PML") {

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
                         test           = test,
                         mimic          = lavoptions$mimic,
                         method         = "ABA",
                         return.ugamma  = FALSE)
        TEST[[2]] <- out[[test]]

    } else if(test %in% c("yuan.bentler", "yuan.bentler.mplus") && df > 0 &&
              lavoptions$estimator != "PML") {

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
                         test           = test,
                         mimic          = lavoptions$mimic,
                         #method         = "default",
                         return.ugamma  = FALSE)
        TEST[[2]] <- out[[test]]

    } else if(test == "bootstrap" || test == "bollen.stine") {
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
                bootstrap.internal(object          = NULL,
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

        TEST[[2]] <- list(test="bootstrap",
                          stat=chisq,
                          stat.group=chisq.group,
                          df=df,
                          pvalue=pvalue.boot,
                          boot.T=BOOT.TEST,
                          boot.larger=boot.larger,
                          boot.length=boot.length)
    }

    TEST
}


