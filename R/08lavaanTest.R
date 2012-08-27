testStatisticSatorraBentler <- function(samplestats=samplestats, 
                                        E.inv, Delta, WLS.V,
                                        x.idx=list(integer(0))) {

    # warn if fixed.x!
    #if(length(x.idx[[1L]] > 0L)) {
    #    warning("lavaan WARNING: SB scaling factor may not be correct in the presence of exogenous fixed.x covariates; either use fixed.x=FALSE or mimic=Mplus to get better results")
    #}

    Gamma <- samplestats@NACOV

    trace.UGamma  <- numeric( samplestats@ngroups )
    trace.UGamma2 <- numeric( samplestats@ngroups )

    for(g in 1:samplestats@ngroups) {
        tmp <- WLS.V[[g]] %*% Delta[[g]] %*% E.inv
        U <- WLS.V[[g]] - (tmp %*% t(Delta[[g]]) %*% WLS.V[[g]])
        UG <- U %*% Gamma[[g]]; tUG <- t(UG)
        trace.UGamma[g] <-
                samplestats@ntotal/samplestats@nobs[[g]] * sum( U * Gamma[[g]] )
        trace.UGamma2[g] <-
                samplestats@ntotal/samplestats@nobs[[g]] * sum( UG * tUG )
    }

    attr(trace.UGamma, "trace.UGamma2") <- trace.UGamma2

    trace.UGamma
}


testStatisticSatorraBentler.Mplus <- function(samplestats=samplestats, 
                                              E.inv, Delta, WLS.V,
                                              x.idx=list(integer(0))) {
    Gamma <- samplestats@NACOV
    
    trace.UGamma <- numeric( samplestats@ngroups )
    trace.UGamma2 <- numeric( samplestats@ngroups )

    for(g in 1:samplestats@ngroups) {
        A1 <- WLS.V[[g]]
        B1 <- A1 %*% Gamma[[g]] %*% A1

        # A0 <- t(Delta) %*% A1 %*% Delta  ## A0 = E
        B0 <- t(Delta[[g]]) %*% B1 %*% Delta[[g]]

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx[[g]]) > 0L) {
            nvar <- ncol(samplestats@cov[[g]])
            idx <- eliminate.pstar.idx(nvar=nvar, el.idx=x.idx[[g]],
                                       meanstructure=TRUE, type="all")
            A1 <- A1[idx,idx]
            B1 <- B1[idx,idx]
        }
        A1.inv <- solve(A1)

        trace.h1     <- sum( B1 * t(A1.inv) )
        trace.h0     <- sum( B0 * t(E.inv)  ) 
        trace.UGamma[g] <- 
            samplestats@ntotal/samplestats@nobs[[g]] * (trace.h1-trace.h0)
        UG <- (B1 %*% A1.inv) - (B1 %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]]))
        tUG <- t(UG)
        trace.UGamma2[g] <- 
            samplestats@ntotal/samplestats@nobs[[g]] * sum( UG * tUG )
    }        

    attr(trace.UGamma, "trace.UGamma2") <- trace.UGamma2

    trace.UGamma
}

testStatisticYuanBentler <- function(samplestats=samplestats,
                                     A1.group=NULL,
                                     B1.group=NULL,
                                     Delta=NULL,
                                     E.inv=NULL,
                                     x.idx=list(integer(0))) {

    # we always assume a meanstructure
    meanstructure <- TRUE

    trace.UGamma <- numeric( samplestats@ngroups )
    trace.h1     <- numeric( samplestats@ngroups )
    trace.h0     <- numeric( samplestats@ngroups )

    for(g in 1:samplestats@ngroups) {
        A1 <- A1.group[[g]]
        B1 <- B1.group[[g]]

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx[[g]]) > 0L) {
            nvar <- ncol(samplestats@cov[[g]])
            idx <- eliminate.pstar.idx(nvar=nvar, el.idx=x.idx[[g]],
                                       meanstructure=meanstructure, type="all")
            A1 <- A1[idx,idx]
            B1 <- B1[idx,idx]
        }

        trace.h1[g]     <- sum( B1 * t( solve(A1) ) )
        trace.h0[g]     <- sum( (B1 %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]]) %*% A1) * t( solve(A1) ) )
                           
        trace.UGamma[g] <- (trace.h1[g] - trace.h0[g])
    }

    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0

    trace.UGamma
}

testStatisticYuanBentler.Mplus <- function(samplestats=samplestats, 
                                           data=data,
                                           information="observed",
                                           B0.group=NULL,
                                           E.inv=NULL,
                                           x.idx=list(integer(0))) {
    # typical for Mplus:
    # - do NOT use the YB formula, but use an approximation
    #   relying  on A0 ~= Delta'*A1*Delta and the same for B0

    # we always assume a meanstructure
    meanstructure <- TRUE

    trace.UGamma <- numeric( samplestats@ngroups )
    trace.h1     <- numeric( samplestats@ngroups )
    trace.h0     <- numeric( samplestats@ngroups )

    for(g in 1:samplestats@ngroups) {
        # if data is complete, A1.22 is simply 0.5*D'(S.inv x S.inv)D
        A1 <- compute.A1.sample(samplestats=samplestats, group=g,
                               meanstructure=meanstructure,
                               information=information)
        B1 <- compute.B1.sample(samplestats=samplestats, data=data, group=g)

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx[[g]]) > 0L) {
            nvar <- ncol(samplestats@cov[[g]])
            idx <- eliminate.pstar.idx(nvar=nvar, el.idx=x.idx[[g]],
                                       meanstructure=meanstructure, type="all")
            A1 <- A1[idx,idx]
            B1 <- B1[idx,idx]
        }

        trace.h1[g]     <- sum( B1 * t( solve(A1) ) )
        trace.h0[g]     <- ( samplestats@nobs[[g]]/samplestats@ntotal *
                             sum( B0.group[[g]] * t(E.inv) ) )
        trace.UGamma[g] <- (trace.h1[g] - trace.h0[g])
    }

    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0

    trace.UGamma
}


           

computeTestStatistic <- function(object, partable=NULL, samplestats=NULL, 
                                 options=NULL, x=NULL, VCOV=NULL,
                                 data=NULL, control=list()) {

    mimic       <- options$mimic
    test        <- options$test
    information <- options$information
    estimator   <- options$estimator


    TEST <- list()

    # degrees of freedom
    df <- getDF(partable)

    # handle constraints
    if(nrow(object@con.jac) > 0L) {
        # Mplus seems to ignore inequality constraints
        if(options$mimic == "Mplus") {
            df <- ( df + length(attr(object@con.jac, "ceq.idx")) )
        } else {
            df <- ( df + length(attr(object@con.jac, "ceq.idx")) 
                       + length(attr(object@con.jac, "cin.idx"))
                       # FIXME:
                       # we should not count inactive constraints?
                       - length(attr(object@con.jac, "inactive.idx")) )
        }
    }

    if(test == "none") {
        TEST[[1]] <- list(test=test,
                          stat=as.numeric(NA),
                          stat.group=as.numeric(NA),
                          df=df,
                          pvalue=as.numeric(NA))
        return(TEST)
    }    

    # get fx.group
    fx <- attr(x, "fx")
    fx.group <- attr(fx, "fx.group")

    # always compute `standard' test statistic
    ## FIXME: the NFAC is now implicit in the computation of fx...
    NFAC <- 2 * unlist(samplestats@nobs)
    if(options$estimator == "ML" && options$likelihood == "wishart") {
        # first divide by two
        NFAC <- NFAC / 2
        NFAC <- NFAC - 1
        NFAC <- NFAC * 2
    }

    chisq.group <- fx.group * NFAC

    # check for negative values
    chisq.group[which(chisq.group < 0)] <- 0.0

    # global test statistic
    chisq <- sum(chisq.group)

    # pvalue  ### FIXME: what if df=0? NA? or 1? or 0?
    pvalue <- 1 - pchisq(chisq, df)

    TEST[[1]] <- list(test="standard",
                      stat=chisq, 
                      stat.group=chisq.group, 
                      df=df, 
                      pvalue=pvalue) 

    if(df == 0 && test %in% c("satorra.bentler", "yuan.bentler",
                              "mean.var.adjusted", "scaled.shifted")) {
        TEST[[2]] <- list(test=test,
                          stat=chisq, 
                          stat.group=chisq.group,
                          df=df, 
                          pvalue=pvalue,
                          scaling.factor=as.numeric(NA))
    }

    # some require meanstructure (for now)
    if(test %in% c("satorra.bentler", "yuan.bentler") &&
       !options$meanstructure) {
        stop("test (", test, ") requires meanstructure (for now)")
    }

    # fixed.x idx
    x.idx <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        if(options$fixed.x) {
            x.idx[[g]] <- match(vnames(partable, "ov.x", group=g), 
                                data@ov.names[[g]])
        } else {
            x.idx[[g]] <- integer(0L)
        }
    }

    if(test %in% c("satorra.bentler", "mean.var.adjusted", "scaled.shifted") &&
       df > 0) {
        # try to extract attr from VCOV (if present)
        E.inv <- attr(VCOV, "E.inv")
        Delta <- attr(VCOV, "Delta")
        WLS.V <- attr(VCOV, "WLS.V")

        # if not present (perhaps se.type="standard" or se.type="none")
        #  we need to compute these again
        if(is.null(E.inv) || is.null(Delta) || is.null(WLS.V)) {
            if(mimic == "Mplus" && estimator == "ML") {
                # special treatment for Mplus
                E <- computeExpectedInformationMLM(object, 
                                                   samplestats=samplestats)
            } else {
                E <- computeExpectedInformation(object, samplestats=samplestats,
                                                data=data,
                                                estimator=estimator,
                                                extra=TRUE)
            }
            E.inv <- solve(E)
            Delta <- attr(E, "Delta")
            WLS.V <- attr(E, "WLS.V")
        }

        if(mimic == "Mplus" && estimator == "ML") {
            ### TESTING ONLY FOR FIXED.X !!! ###
#            if(length(x.idx) > 0L) {
#                cat("\n\nDEBUG FIXED.X\n\n\n")
#                augUser <- user
#                idx <- which(augUser$exo > 0L)
#                augUser$exo[       idx ] <- 0L
#                augUser$free[      idx ] <- max(augUser$free) + 1:length(idx)
#                augUser$unco[idx ] <- max(augUser$unco) + 1:length(idx)
#                augModel <- Model(user           = augUser,
#                                  start          = getModelParameters(object, type="user"),
#                                  representation = object@representation,
#                                  debug          = FALSE)
#
#                Delta <- computeDelta(augModel)
#                E <- computeExpectedInformationMLM(object, samplestats=samplestats,
#                                                   Delta=Delta)
#                fixed.x.idx <- max(partable$free) + 1:length(idx)
#                free.idx    <- 1:max(partable$free)
#                E[free.idx, fixed.x.idx] <- 0.0
#                E[fixed.x.idx, free.idx] <- 0.0
#                E.inv <- solve(E)
#                x.idx <- integer(0)
#            }
            trace.UGamma <-
                testStatisticSatorraBentler.Mplus(samplestats = samplestats,
                                                  E.inv       = E.inv, 
                                                  Delta       = Delta, 
                                                  WLS.V       = WLS.V, 
                                                  x.idx       = x.idx)
        } else {
            trace.UGamma <- 
                testStatisticSatorraBentler(samplestats = samplestats,
                                            E.inv       = E.inv, 
                                            Delta       = Delta, 
                                            WLS.V       = WLS.V, 
                                            x.idx       = x.idx)
        }

        
        trace.UGamma2 <- attr(trace.UGamma, "trace.UGamma2")
        attributes(trace.UGamma) <- NULL

        # adjust df?
        if(test == "mean.var.adjusted") {
            if(mimic == "Mplus") {
                df <- floor((sum(trace.UGamma)^2 / sum(trace.UGamma2)) + 0.5)
            } else {
                # more precise, fractional df
                df <- sum(trace.UGamma)^2 / sum(trace.UGamma2)
            }
        } else if(test == "satorra.bentler") {
            trace.UGamma2 <- as.numeric(NA)
        }

        if(test == "scaled.shifted") {
            # this is the T3 statistic as used by Mplus 6 and higher
            # see 'Simple Second Order Chi-Square Correction' 2010
            # www.statmodel.com
            scaling.factor <- sqrt(sum(trace.UGamma2) /df)
            if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
            shift.parameter <- max(0, df - sqrt( df*sum(trace.UGamma)^2 /
                                      sum(trace.UGamma2) ) )
            stat.group      <- (chisq.group / scaling.factor + shift.parameter)
            chisq.scaled    <- sum(stat.group)
            pvalue.scaled   <- 1 - pchisq(chisq.scaled, df)
        } else {
            scaling.factor <- sum(trace.UGamma) / df
            if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
            stat.group     <- chisq.group / scaling.factor
            chisq.scaled   <- sum(stat.group)
            pvalue.scaled  <- 1 - pchisq(chisq.scaled, df)
            shift.parameter <- as.numeric(NA)
        }

        TEST[[2]] <- list(test=test,
                          stat=chisq.scaled,
                          stat.group=stat.group,
                          df=df,
                          pvalue=pvalue.scaled,
                          scaling.factor=scaling.factor,
                          shift.parameter=shift.parameter,
                          trace.UGamma=trace.UGamma,
                          trace.UGamma2=trace.UGamma2)

    } else if(test == "yuan.bentler" && df > 0) {
        # try to extract attr from VCOV (if present)
        E.inv    <- attr(VCOV, "E.inv")
        B0.group <- attr(VCOV, "B0.group")

        if(is.null(E.inv)) {
            # if se="standard", information is probably expected
            # change it to observed
            if(options$se != "robust.mlr") information <- "observed"
            E.inv <- Nvcov.standard(object      = object,
                                    samplestats = samplestats,
                                    data        = data,
                                    estimator   = "ML",
                                    information = information)
        }

        if(mimic == "Mplus" || mimic == "lavaan") {
            if(is.null(B0.group)) {
                Nvcov <- Nvcov.first.order(object      = object, 
                                           samplestats = samplestats,
                                           data        = data)
                B0.group <- attr(Nvcov, "B0.group")
            } 
            trace.UGamma <- 
                testStatisticYuanBentler.Mplus(samplestats = samplestats,
                                               data        = data,
                                               information = information,
                                               B0.group    = B0.group,
                                               E.inv       = E.inv,
                                               x.idx       = x.idx)
        } else {
            Delta <- computeDelta(object)
            Sigma.hat <- computeSigmaHat(object)
            if(object@meanstructure) Mu.hat <- computeMuHat(object)
            A1.group <- vector("list", length=samplestats@ngroups)
            B1.group <- vector("list", length=samplestats@ngroups)
            for(g in 1:samplestats@ngroups) {
                if(samplestats@missing.flag) {
                    X <- NULL
                    M <- samplestats@missing[[g]]
                } else {
                    X <- data@X[[g]]
                    M <- NULL
                }
                out <- compute.Abeta.Bbeta(Sigma.hat=Sigma.hat[[g]], 
                                           Mu.hat=Mu.hat[[g]], 
                                           X=X,
                                           M=M,
                                           Abeta=TRUE, Bbeta=TRUE, 
                                           information=information)
                A1.group[[g]] <- out$Abeta
                B1.group[[g]] <- out$Bbeta
            }
            trace.UGamma <-
                testStatisticYuanBentler(samplestats=samplestats,
                                         A1.group=A1.group,
                                         B1.group=B1.group,
                                         Delta=Delta,
                                         E.inv=E.inv,
                                         x.idx=x.idx)
        }

        scaling.factor       <- sum(trace.UGamma) / df
        if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
        chisq.scaled         <- sum(chisq.group / scaling.factor)
        pvalue.scaled        <- 1 - pchisq(chisq.scaled, df)

        ndat <- getNDAT(partable)
        npar <- getNPAR(partable)

        scaling.factor.h1    <- sum( attr(trace.UGamma, "h1") ) / ndat
        scaling.factor.h0    <- sum( attr(trace.UGamma, "h0") ) / npar

        TEST[[2]] <- list(test=test,
                          stat=chisq.scaled,
                          stat.group=(chisq.group / scaling.factor),
                          df=df,
                          pvalue=pvalue.scaled,
                          scaling.factor=scaling.factor,
                          scaling.factor.h1=scaling.factor.h1,
                          scaling.factor.h0=scaling.factor.h0,
                          trace.UGamma=trace.UGamma)

    } else if(test == "bootstrap" || test == "bollen.stine") {
        # check if we have bootstrap data
        BOOT.TEST <- attr(VCOV, "BOOT.TEST")
        if(is.null(BOOT.TEST)) {
            if(!is.null(options$bootstrap)) {
                R <- options$bootstrap
            } else {
                R <- 1000L
            }
            boot.type <- "bollen.stine"
            BOOT.TEST <- 
                bootstrap.internal(object=NULL,
                                   model.=object, samplestats.=samplestats, 
                                   partable.=partable,
                                   options.=options, data.=data,
                                   R=R, verbose=options$verbose,
                                   type=boot.type,
                                   FUN ="test",
                                   warn=-1L,
                                   parallel=control$parallel,
                                   ncpus=control$ncpus,
                                   cl=control$cl)
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


