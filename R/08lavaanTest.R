testStatisticSatorraBentler <- function(samplestats=samplestats, 
                                        E.inv, Delta, WLS.V,
                                        x.idx=list(integer(0))) {

    # UG = Gamma %*% [V - V %*% Delta %*% E.inv %*% tDelta %*% V]
    #    = Gamma %*% V  - Gamma %*% V %*% Delta %*% E.inv %*% tDelta %*% V
    #    = Gamma %*% A1 - Gamma %*% A1 %*% Delta %*% E.inv %*% tDelta %*% A1
    # (B1 = A1 %*% Gamma %*% A1)
    #    = B1 %*% A1.inv - B1 %*% A1.inv %*% Delta %*% E.inv %*% tDelta %*% A1
    #    
    # if only the trace is needed, we can use reduce the rhs (after the minus)
    # to B1 %*% Delta %*% E.inv %*% tDelta (eliminating A1 and A1.inv)

    # we write it like this to allow for fixed.x covariates which affect A1
    # and B1

    Gamma <- samplestats@NACOV
    nss <- ncol(Gamma[[1]])
    ngroups <- samplestats@ngroups

    UG.group  <- vector("list", length=ngroups)
    trace.UGamma2 <- numeric(ngroups)
    for(g in 1:ngroups) {
        fg  <-  samplestats@nobs[[g]]   /samplestats@ntotal
        fg1 <- (samplestats@nobs[[g]]-1)/samplestats@ntotal
        WLS.Vg  <- WLS.V[[g]] * fg
        Gamma.g <- Gamma[[g]] / fg  ## ?? check this
        Delta.g <- Delta[[g]]

        # just for testing: regular UG:
        UG1 <- Gamma.g %*% (WLS.Vg - WLS.Vg %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]]) %*% WLS.Vg)

        A1      <- WLS.V[[g]] * fg
        B1 <- A1 %*% Gamma.g %*% A1
        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx[[g]]) > 0L) {
            nvar <- ncol(samplestats@cov[[g]])
            idx <- eliminate.pstar.idx(nvar=nvar, el.idx=x.idx[[g]],
                                       meanstructure=TRUE, type="all")
            A1      <- A1[idx,idx]
            B1      <- B1[idx,idx]
            Delta.g <- Delta.g[idx,]
        } 
        A1.inv <- solve(A1)

        #B0 <- t(Delta[[g]]) %*% B1 %*% Delta[[g]]
        #trace.h1     <- sum( B1 * t(A1.inv) )
        #trace.h0     <- sum( B0 * t(E.inv)  ) 
        #trace.UGamma[g] <- (trace.h1-trace.h0)

        tmp <- (B1 %*% A1.inv) - (B1 %*% Delta.g %*% E.inv %*% t(Delta.g))
        # sanity check 1: sum(diag(UG1)) - sum(diag(tmp))
        # sanity check 2: sum(diag(UG1 %*% UG1)) - sum(diag(tmp %*% tmp))
        trace.UGamma2[g] <- sum(tmp * t(tmp))

        UG.group[[g]] <- tmp
    }

    # NOTE: if A, B, C are matrices
    # tr(A+B+C) = tr(A) + tr(B) + tr(C)
    #
    # BUT:
    # tr( (A+B+C)^2 ) != tr(A^2) + tr(B^2) + tr(C^2) 
    # it would seem that we need the latter... (trace.UGamma3) for MLMV and
    # friends

    UG <- UG.group[[1]]
    if(ngroups > 1L) {
        for(g in 2:ngroups) {
            UG <- UG + UG.group[[g]]
        }
    }

    # trace
    trace.UGamma <- sum(diag(UG))

    # for mean and variance adjusted tr UG^2 (per group)
    trace.UGamma4 <- trace.UGamma2
    trace.UGamma2 <- sum(trace.UGamma2) # at least for MLMV in Mplus
                                        # this is what is needed when multiple
                                        # groups are used?? 
    attr(trace.UGamma, "trace.UGamma2") <- trace.UGamma2
    attr(trace.UGamma, "trace.UGamma4") <- trace.UGamma4

    # testing only -- alternative interpretation of tr UG^2
    tUG <- t(UG); trace.UGamma3 <- sum(UG * tUG) # seems wrong?
    attr(trace.UGamma, "trace.UGamma3") <- trace.UGamma3

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
        A1.inv <- solve(A1)

        trace.h1[g] <- sum( B1 * t( A1.inv ) )
        trace.h0[g] <- sum( (B1 %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]])) )
        #trace.UGamma[g] <- (trace.h1[g] - trace.h0[g])
        UG <- (B1 %*% A1.inv) - (B1 %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]]))
    }

    # traces
    trace.UGamma <- sum(diag(UG))
    #tUG <- t(UG); trace.UGamma2 <- sum(UG * tUG)
    #attr(trace.UGamma, "trace.UGamma2") <- trace.UGamma2

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

    # we take the sum here
    trace.UGamma <- sum(trace.UGamma)

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

    # handle equality constraints (note: we ignore inequality constraints, 
    # active or not!)
    if(nrow(object@con.jac) > 0L) {
        df <- ( df + length(attr(object@con.jac, "ceq.idx")) )
    }

    if(test == "none") {
        TEST[[1]] <- list(test=test,
                          stat=as.numeric(NA),
                          stat.group=as.numeric(NA),
                          df=df,
                          refdist="unknown",
                          pvalue=as.numeric(NA))
        return(TEST)
    }    

    # get fx.group
    fx <- attr(x, "fx")
    fx.group <- attr(fx, "fx.group")

    if(estimator != "PML") {
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
    } else {
        # for estimator PML, we compute the loglikelihood for the 
        # `unrestricted' model, and then compute the LRT (-2 logl - logl_un)
        group.fx <- numeric( samplestats@ngroups )
        for(g in 1:samplestats@ngroups) {
            group.fx <- estimator.PML(Sigma.hat = samplestats@cov[[g]],
                                      TH        = samplestats@th[[g]],
                                      th.idx    = object@th.idx[[g]],
                                      num.idx   = object@num.idx[[g]],
                                      X         = data@X[[g]])
        }
        chisq.group <- 2 * (fx.group - group.fx) # LRT per group

        #cat("model fx = \n"); print(fx.group)
        #cat("unres fx = \n"); print(group.fx)
    }

    # check for negative values
    chisq.group[which(chisq.group < 0)] <- 0.0

    # global test statistic
    chisq <- sum(chisq.group)

    # reference distribution: always chi-square, except for the
    # non-robust version of ULS
    if(estimator == "ULS") {
        refdistr <- "unknown"
        pvalue <- as.numeric(NA)
    } else {
        refdistr <- "chisq"
        # pvalue  ### FIXME: what if df=0? NA? or 1? or 0?
        pvalue <- 1 - pchisq(chisq, df)
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

    # some require meanstructure (for now)
    if(test %in% c("satorra.bentler", "yuan.bentler") &&
       !options$meanstructure) {
        stop("test (", test, ") requires meanstructure (for now)")
    }

    # fixed.x idx
    x.idx <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        if(options$fixed.x && estimator == "ML") {
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
            E.inv <- try(solve(E))
            if(inherits(E.inv, "try-error")) {
                TEST[[2]] <- list(test=test, stat=as.numeric(NA), 
                    stat.group=rep(as.numeric(NA), samplestats@ngroups),
                    df=df, refdistr=refdistr, pvalue=as.numeric(NA),
                    scaling.factor=as.numeric(NA))
                warning("lavaan WARNING: could not compute scaled test statistic\n")
                return(TEST)                
            }
            Delta <- attr(E, "Delta")
            WLS.V <- attr(E, "WLS.V")
        }

#        if(mimic == "Mplus" && estimator == "ML") {
            ### TESTING ONLY FOR FIXED.X !!! ###
#            if(length(x.idx) > 0L) {
#                cat("\n\nDEBUG FIXED.X\n\n\n")
#                augUser <- user
#                idx <- which(augUser$exo > 0L)
#                augUser$exo[       idx ] <- 0L
#                augUser$free[      idx ] <- max(augUser$free) + 1:length(idx)
#                augUser$unco[idx ] <- max(augUser$unco) + 1:length(idx)
#                augModel <- Model(partable       = augUser,
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
#        } 

        trace.UGamma <- 
            testStatisticSatorraBentler(samplestats = samplestats,
                                        E.inv       = E.inv, 
                                        Delta       = Delta, 
                                        WLS.V       = WLS.V, 
                                        x.idx       = x.idx)
        trace.UGamma2 <- attr(trace.UGamma, "trace.UGamma2")
        trace.UGamma3 <- attr(trace.UGamma, "trace.UGamma3")
        trace.UGamma4 <- attr(trace.UGamma, "trace.UGamma4")
        attributes(trace.UGamma) <- NULL

        # adjust df?
        if(test == "mean.var.adjusted") {
            if(mimic == "Mplus") {
                df <- floor(trace.UGamma^2/trace.UGamma2 + 0.5)
            } else {
                # more precise, fractional df
                df <- trace.UGamma^2 / trace.UGamma2
            }
        } else if(test == "satorra.bentler") {
            trace.UGamma2 <- as.numeric(NA)
        }

        if(test == "scaled.shifted") {
            # this is the T3 statistic as used by Mplus 6 and higher
            # see 'Simple Second Order Chi-Square Correction' 2010
            # www.statmodel.com

            # however, for multiple groups, Mplus reports something else
            # YR. 30 Aug 2012 -- after much trial and error, it turns out
            # that the shift-parameter (b) is weighted (while a is not)??
            # however, the chisq.square per group are different; only
            # the sum seems ok??
            fg <- unlist(samplestats@nobs)/samplestats@ntotal
           
            a <- sqrt(df/trace.UGamma2)
            #dprime <- trace.UGamma^2 / trace.UGamma2
            #shift.parameter <- fg * (df - sqrt(df * dprime))
            shift.parameter <- fg * (df - a*trace.UGamma)
            scaling.factor  <- 1/a
            stat.group      <- (chisq.group * a + shift.parameter)
            chisq.scaled    <- sum(stat.group)
            pvalue.scaled   <- 1 - pchisq(chisq.scaled, df)
        } else {
            scaling.factor <- trace.UGamma/df
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
                          trace.UGamma4=trace.UGamma4,
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


