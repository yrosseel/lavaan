
testStatisticYuanBentler <- function(lavsamplestats =lavsamplestats,
                                     meanstructure  = TRUE,
                                     A1.group       = NULL,
                                     B1.group       = NULL,
                                     Delta          = NULL,
                                     E.inv          = NULL,
                                     x.idx          = list(integer(0)),
                                     Satterthwaite  = FALSE) {

    # we always assume a meanstructure (nope, not any longer, since 0.6)
    #meanstructure <- TRUE

    trace.UGamma  <- numeric( lavsamplestats@ngroups )
    trace.UGamma2 <- numeric( lavsamplestats@ngroups )
    trace.h1      <- numeric( lavsamplestats@ngroups )
    trace.h0      <- numeric( lavsamplestats@ngroups )

    for(g in 1:lavsamplestats@ngroups) {
        A1 <- A1.group[[g]]
        B1 <- B1.group[[g]]

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx[[g]]) > 0L) {
            nvar <- ncol(lavsamplestats@cov[[g]])
            idx <- eliminate.pstar.idx(nvar=nvar, el.idx=x.idx[[g]],
                                       meanstructure=meanstructure, type="all")
            A1 <- A1[idx,idx]
            B1 <- B1[idx,idx]
        }
        A1.inv <- solve(A1)

        trace.h1[g] <- sum( B1 * t( A1.inv ) )
        trace.h0[g] <- sum( B1 * Delta[[g]] %*% E.inv %*% t(Delta[[g]]) )
        trace.UGamma[g] <- trace.h1[g] - trace.h0[g]

        if(Satterthwaite) {
            UG <- (A1.inv %*% B1) - (A1.inv %*% B1 %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]]) %*% A1)
            trace.UGamma2[g] <- sum(UG * t(UG))
        }
    }

    # traces
    trace.UGamma <- sum(trace.UGamma)
    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0

    if(Satterthwaite) {
        attr(trace.UGamma, "trace.UGamma2") <- sum(trace.UGamma2)
    }

    trace.UGamma
}

testStatisticYuanBentler.Mplus <- function(lavsamplestats=lavsamplestats, 
                                           lavdata=lavdata,
                                           information="observed",
                                           B0.group=NULL,
                                           E.inv=NULL,
                                           x.idx=list(integer(0)),
                                           meanstructure = TRUE) {
    # typical for Mplus:
    # - do NOT use the YB formula, but use an approximation
    #   relying  on A0 ~= Delta'*A1*Delta and the same for B0

    # we always assume a meanstructure (no, not any longer since 0.6)
    #meanstructure <- TRUE
    ngroups <- lavsamplestats@ngroups

    trace.UGamma <- numeric( lavsamplestats@ngroups )
    trace.h1     <- numeric( lavsamplestats@ngroups )
    trace.h0     <- numeric( lavsamplestats@ngroups )

    for(g in 1:lavsamplestats@ngroups) {
        # if lavdata is complete, A1.22 is simply 0.5*D'(S.inv x S.inv)D
        if(lavsamplestats@missing.flag) {
            if(information == "expected") {
                A1 <- lav_mvnorm_missing_information_expected(
                          Y     = lavdata@X[[g]],
                          Mp    = lavdata@Mp[[g]],
                          wt    = lavdata@weights[[g]],
                          Mu    = lavsamplestats@missing.h1[[g]]$mu,
                          Sigma = lavsamplestats@missing.h1[[g]]$sigma)
            } else {
                A1 <- lav_mvnorm_missing_information_observed_samplestats(
                          Yp    = lavsamplestats@missing[[g]],
                          Mu    = lavsamplestats@missing.h1[[g]]$mu,
                          Sigma = lavsamplestats@missing.h1[[g]]$sigma)
            }
        } else {
            # data complete, under h1, expected == observed
            A1 <- lav_mvnorm_h1_information_observed_samplestats(
                      sample.cov     = lavsamplestats@cov[[g]],
                      sample.cov.inv = lavsamplestats@icov[[g]],
                      meanstructure  = meanstructure)
        }     

        if(lavsamplestats@missing.flag) {
            B1 <- lav_mvnorm_missing_information_firstorder(
                          Y     = lavdata@X[[g]],
                          Mp    = lavdata@Mp[[g]],
                          wt    = lavdata@weights[[g]],
                          Mu    = lavsamplestats@missing.h1[[g]]$mu,
                          Sigma = lavsamplestats@missing.h1[[g]]$sigma)
        } else {
            # B1 =  A1 %*% Gamma %*% A1
            B1 <- lav_mvnorm_h1_information_firstorder(
                          Y     = lavdata@X[[g]],
                          wt    = lavdata@weights[[g]],
                          Gamma = lavsamplestats@NACOV[[g]],
                          sample.cov = lavsamplestats@cov[[g]],
                          sample.cov.inv = lavsamplestats@icov[[g]],
                          meanstructure  = meanstructure)
        }

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx[[g]]) > 0L) {
            nvar <- ncol(lavsamplestats@cov[[g]])
            idx <- eliminate.pstar.idx(nvar=nvar, el.idx=x.idx[[g]],
                                       meanstructure=meanstructure, type="all")
            A1 <- A1[idx,idx]
            B1 <- B1[idx,idx]
        }
        A1.inv <- solve(A1)

        # if data is complete, why not just A1 %*% Gamma?
        trace.h1[g]     <- sum( B1 * t( A1.inv ) )
        trace.h0[g]     <- ( lavsamplestats@nobs[[g]]/lavsamplestats@ntotal *
                             sum( B0.group[[g]] * t(E.inv) ) )
        trace.UGamma[g] <- (trace.h1[g] - trace.h0[g])
    }

    # we take the sum here
    trace.UGamma <- sum(trace.UGamma)

    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0

    trace.UGamma
}


           

lav_model_test <- function(lavmodel       = NULL, 
                           lavpartable    = NULL, 
                           lavsamplestats = NULL, 
                           lavoptions     = NULL, 
                           x              = NULL, 
                           VCOV           = NULL, 
                           lavcache       = NULL,
                           lavdata        = NULL,
                           h1             = list(),
                           loglik         = NULL,
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
        if(length(h1) > 0L) {
            # LRT
            chisq.group <- -2 * (loglik$loglik.group - h1$loglik.group)
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
    # non-robust version of ULS
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
                           
    } else if(test == "yuan.bentler" && df > 0) {
        # try to extract attr from VCOV (if present)
        E.inv    <- attr(VCOV, "E.inv")
        B0.group <- attr(VCOV, "B0.group")

        if(is.null(E.inv)) {
            # if se="standard", information is probably expected
            # change it to observed
            # if(lavoptions$se != "robust.mlr") information <- "observed"
            # NO longer, since 0.6-1 
            E.inv <- lav_model_information(lavmodel       = lavmodel,
                                           lavsamplestats = lavsamplestats,
                                           lavdata        = lavdata,
                                           lavcache       = lavcache,
                                           lavoptions     = lavoptions,
                                           extra          = FALSE,
                                           augmented      = TRUE,
                                           inverted       = TRUE)
        }

        if(inherits(E.inv, "try-error")) {
            TEST[[2]] <- list(test=test, stat=as.numeric(NA),
                stat.group=rep(as.numeric(NA), lavsamplestats@ngroups),
                df=df, refdistr=refdistr, pvalue=as.numeric(NA),
                scaling.factor=as.numeric(NA))
            warning("lavaan WARNING: could not compute scaled test statistic\n")
            return(TEST)
        }

        if(mimic == "Mplus") { # since 0.6-1
            if(is.null(B0.group)) {
                B0 <- 
                    lav_model_information_firstorder(lavmodel = lavmodel,
                                             lavsamplestats = lavsamplestats,
                                             lavdata        = lavdata,
                                             lavcache       = lavcache,
                                             lavoptions     = lavoptions,
                                             extra          = TRUE,
                                             check.pd       = FALSE,
                                             augmented      = FALSE,
                                             inverted       = FALSE)
                B0.group <- attr(B0, "B0.group")
            }
            trace.UGamma <- 
                testStatisticYuanBentler.Mplus(lavsamplestats = lavsamplestats,
                                               lavdata        = lavdata,
                                               information    = information,
                                               B0.group       = B0.group,
                                               E.inv          = E.inv,
                                               x.idx          = x.idx,
                                               meanstructure  = meanstructure)
        } else {
            if(!is.null(lavdata@weights[[1]])) {
                stop("lavaan ERROR: weights not supported yet using mimic!")
            }
            Delta <- computeDelta(lavmodel = lavmodel)
            Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
            if(lavmodel@meanstructure) {
                Mu.hat <- computeMuHat(lavmodel = lavmodel)
            }
            A1.group <- vector("list", length=lavsamplestats@ngroups)
            B1.group <- vector("list", length=lavsamplestats@ngroups)
            for(g in 1:lavsamplestats@ngroups) {
                if(lavsamplestats@missing.flag) {
                    if(lavoptions$h1.information == "structured") {
                        SIGMA <- Sigma.hat[[g]]
                        MEAN  <- Mu.hat[[g]]
                    } else {
                        SIGMA <- lavsamplestats@missing.h1[[g]]$sigma
                        MEAN  <- lavsamplestats@missing.h1[[g]]$mu
                    }
                    out <- lav_mvnorm_missing_information_both(
                               Y  = lavdata@X[[g]],
                               Mp = lavdata@Mp[[g]],
                               wt = lavdata@weights[[g]],
                               Mu = MEAN,
                               Sigma = SIGMA,
                               information = information)
                    A1.group[[g]] <- out$Abeta
                    B1.group[[g]] <- out$Bbeta
                } else {
                    # different scenarios:
                    # - meanstructure + structured 
                    # - meanstructure + unstructured
                    # - no meanstructure + structured 
                    # - no meanstructure + unstructured
                    if(meanstructure &&
                       lavoptions$h1.information == "structured") {
                        if(information == "expected") {
                             A1.group[[g]] <- lav_mvnorm_information_expected(
                                Y = lavdata@X[[g]],
                                Mu = Mu.hat[[g]], Sigma = Sigma.hat[[g]])
                        } else {
                            A1.group[[g]] <- 
                                lav_mvnorm_information_observed_samplestats(
                                    sample.mean = lavsamplestats@mean[[g]],
                                    sample.cov = lavsamplestats@cov[[g]],
                                    Mu = Mu.hat[[g]], Sigma = Sigma.hat[[g]])
                        }
                        B1.group[[g]] <- lav_mvnorm_information_firstorder(
                                Y = lavdata@X[[g]],
                                Mu = Mu.hat[[g]], Sigma = Sigma.hat[[g]],
                                wt = lavdata@weights[[g]])
                    } else if(!meanstructure &&
                              lavoptions$h1.information == "structured") {
                        # no meanstructure
                        if(information == "expected") {
                             A1.group[[g]] <- lav_mvnorm_information_expected(
                                Y = lavdata@X[[g]],
                                Mu = lavsamplestats@mean[[g]],
                                Sigma = Sigma.hat[[g]],
                                meanstructure = FALSE)
                        } else {
                            A1.group[[g]] <-
                                lav_mvnorm_information_observed_samplestats(
                                    sample.mean = lavsamplestats@mean[[g]],
                                    sample.cov = lavsamplestats@cov[[g]],
                                    Mu = lavsamplestats@mean[[g]],
                                    Sigma = Sigma.hat[[g]],
                                    meanstructure = FALSE)
                        }
                        B1.group[[g]] <- lav_mvnorm_information_firstorder(
                                Y             = lavdata@X[[g]],
                                Mu            = lavsamplestats@mean[[g]],
                                Sigma         = Sigma.hat[[g]],
                                wt            = lavdata@weights[[g]],
                                meanstructure = FALSE)
                    } else if(meanstructure &&
                              !lavoptions$h1.information == "structured") {
                        # information expected == observed if h1!!
                        A1.group[[g]] <- lav_mvnorm_h1_information_expected(
                                Y = lavdata@X[[g]], # for wt
                                sample.cov.inv = lavsamplestats@icov[[g]],
                                wt = lavdata@weights[[g]],
                                meanstructure = TRUE)
                        B1.group[[g]] <- lav_mvnorm_h1_information_firstorder(
                                Y = lavdata@X[[g]], # for wt
                                Gamma = lavsamplestats@NACOV[[g]],
                                sample.cov.inv = lavsamplestats@icov[[g]],
                                wt = lavdata@weights[[g]],
                                meanstructure = TRUE)
                    } else if(!meanstructure &&
                              !lavoptions$h1.information == "structured") {
                        # information expected == observed if h1!!
                        A1.group[[g]] <- lav_mvnorm_h1_information_expected(
                                Y = lavdata@X[[g]], # for wt
                                sample.cov.inv = lavsamplestats@icov[[g]],
                                wt = lavdata@weights[[g]],
                                meanstructure = FALSE)
                        B1.group[[g]] <- lav_mvnorm_h1_information_firstorder(
                                Y = lavdata@X[[g]], # for wt
                                Gamma = lavsamplestats@NACOV[[g]],
                                sample.cov.inv = lavsamplestats@icov[[g]],
                                wt = lavdata@weights[[g]],
                                meanstructure = FALSE)
                    }
                }
            }
            trace.UGamma <-
                testStatisticYuanBentler(lavsamplestats = lavsamplestats,
                                         meanstructure  = meanstructure,
                                         A1.group       = A1.group,
                                         B1.group       = B1.group,
                                         Delta          = Delta,
                                         E.inv          = E.inv,
                                         x.idx          = x.idx,
                                         Satterthwaite  = TRUE) # for now
        }

        scaling.factor       <- sum(trace.UGamma) / df
        if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
        chisq.scaled         <- sum(chisq.group / scaling.factor)
        pvalue.scaled        <- 1 - pchisq(chisq.scaled, df)

        ndat <- lav_partable_ndat(lavpartable)
        npar <- lav_partable_npar(lavpartable)

        scaling.factor.h1    <- sum( attr(trace.UGamma, "h1") ) / ndat
        scaling.factor.h0    <- sum( attr(trace.UGamma, "h0") ) / npar
        trace.UGamma2         <- attr(trace.UGamma, "trace.UGamma2")
        attributes(trace.UGamma) <- NULL

        TEST[[2]] <- list(test=test,
                          stat=chisq.scaled,
                          stat.group=(chisq.group / scaling.factor),
                          df=df,
                          pvalue=pvalue.scaled,
                          scaling.factor=scaling.factor,
                          scaling.factor.h1=scaling.factor.h1,
                          scaling.factor.h0=scaling.factor.h0,
                          trace.UGamma=trace.UGamma,
                          trace.UGamma2=trace.UGamma2)

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


