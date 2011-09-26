testStatisticSatorraBentler <- function(sample=sample, 
                                        E.inv, Delta, WLS.V, Gamma, 
                                        x.idx=integer(0)) {
    
    trace.UGamma <- numeric( sample@ngroups )

    for(g in 1:sample@ngroups) {
        tmp <- WLS.V[[g]] %*% Delta[[g]] %*% E.inv
        U <- WLS.V[[g]] - (tmp %*% t(Delta[[g]]) %*% WLS.V[[g]])
        trace.UGamma[g] <-
                sample@ntotal/sample@nobs[[g]] * sum( U * Gamma[[g]] )
    }

    trace.UGamma
}


testStatisticSatorraBentler.Mplus <- function(sample=sample, 
                                              E.inv, Delta, WLS.V, Gamma,
                                              x.idx=integer(0)) {
    
    trace.UGamma <- numeric( sample@ngroups )

    for(g in 1:sample@ngroups) {
        A1 <- WLS.V[[g]]
        B1 <- A1 %*% Gamma[[g]] %*% A1

        # A0 <- t(Delta) %*% A1 %*% Delta  ## A0 = E
        B0 <- t(Delta[[g]]) %*% B1 %*% Delta[[g]]

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx) > 0L) {
            idx <- eliminate.pstar.idx(nvar=sample@nvar, el.idx=x.idx,
                                       meanstructure=TRUE, type="all")
            A1 <- A1[idx,idx]
            B1 <- B1[idx,idx]
        }

        trace.h1     <- sum( B1 * t( solve(A1) ) )
        trace.h0     <- sum( B0 * t(E.inv)       ) 
        trace.UGamma[g] <- sample@ntotal/sample@nobs[[g]] * (trace.h1-trace.h0)
    }        

    trace.UGamma
}

testStatisticYuanBentler <- function(sample=sample,
                                     A1.group=NULL,
                                     B1.group=NULL,
                                     Delta=NULL,
                                     E.inv=NULL,
                                     x.idx=integer(0)) {

    # we always assume a meanstructure
    meanstructure <- TRUE

    trace.UGamma <- numeric( sample@ngroups )
    trace.h1     <- numeric( sample@ngroups )
    trace.h0     <- numeric( sample@ngroups )

    for(g in 1:sample@ngroups) {
        A1 <- A1.group[[g]]
        B1 <- B1.group[[g]]

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx) > 0L) {
            idx <- eliminate.pstar.idx(nvar=sample@nvar, el.idx=x.idx,
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

testStatisticYuanBentler.Mplus <- function(sample=sample, 
                                           information="observed",
                                           B0.group=NULL,
                                           E.inv=NULL,
                                           x.idx=integer(0)) {
    # typical for Mplus:
    # - do NOT use the YB formula, but use an approximation
    #   relying  on A0 ~= Delta'*A1*Delta and the same for B0

    # we always assume a meanstructure
    meanstructure <- TRUE

    trace.UGamma <- numeric( sample@ngroups )
    trace.h1     <- numeric( sample@ngroups )
    trace.h0     <- numeric( sample@ngroups )

    for(g in 1:sample@ngroups) {
        # if data is complete, A1.22 is simply 0.5*D'(S.inv x S.inv)D
        A1 <- compute.A1.sample(sample=sample, group=g,
                               meanstructure=meanstructure,
                               information=information)
        B1 <- compute.B1.sample(sample=sample, group=g)

        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx) > 0L) {
            idx <- eliminate.pstar.idx(nvar=sample@nvar, el.idx=x.idx,
                                       meanstructure=meanstructure, type="all")
            A1 <- A1[idx,idx]
            B1 <- B1[idx,idx]
        }

        trace.h1[g]     <- sum( B1 * t( solve(A1) ) )
        trace.h0[g]     <- ( sample@nobs[[g]]/sample@ntotal *
                             sum( B0.group[[g]] * t(E.inv) ) )
        trace.UGamma[g] <- (trace.h1[g] - trace.h0[g])
    }

    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0

    trace.UGamma
}


           

computeTestStatistic <- function(object, user=NULL, sample=NULL, 
                                 options=NULL, x=NULL, VCOV=NULL,
                                 data=NULL) {

    estimator   <- options$estimator
    mimic       <- options$mimic
    test        <- options$test
    information <- options$information

    TEST <- list()

    # degrees of freedom
    df <- getDF(user)

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
    NFAC <- 2 * unlist(sample@nobs)
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

    # pvalue
    pvalue <- 1 - pchisq(chisq, df)

    TEST[[1]] <- list(test="standard",
                      stat=chisq, 
                      stat.group=chisq.group, 
                      df=df, 
                      pvalue=pvalue) 

    if(df == 0 && test %in% c("satorra.bentler", "yuan.bentler")) {
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
    x.idx <- integer(0)
    if(options$fixed.x) {
         x.idx <- match(vnames(user, "ov.x"), sample@ov.names)
    }

    if(test == "satorra.bentler" && df > 0) {
        # try to extract attr from VCOV (if present)
        E.inv <- attr(VCOV, "E.inv")
        Delta <- attr(VCOV, "Delta")
        WLS.V <- attr(VCOV, "WLS.V")
        Gamma <- attr(VCOV, "Gamma")

        # if not present (perhaps se.type="standard" or se.type="none")
        #  we need to compute these again
        if(is.null(E.inv) || is.null(Delta) || is.null(WLS.V)) {
            if(!mimic == "Mplus") {
                E <- computeExpectedInformation(object, sample=sample,
                                                estimator="ML", Delta=NULL,
                                                extra=TRUE)
                E.inv <- solve(E)
                Delta <- attr(E, "Delta")
                WLS.V <- attr(E, "WLS.V")
            } else {
                # special treatment for Mplus
                E <- computeExpectedInformationMLM(object, sample=sample,
                                                   Delta=NULL)
                E.inv <- solve(E)
                Delta <- attr(E, "Delta")
                WLS.V <- attr(E, "WLS.V")
            }
        }

        if(is.null(Gamma)) {
            Gamma <- vector("list", length=sample@ngroups)
            for(g in 1:sample@ngroups) {
                Gamma[[g]] <- compute.Gamma(sample@data.obs[[g]],
                                            meanstructure=TRUE)
            }
        }
        
        if(mimic == "Mplus") {
            ### TESTING ONLY FOR FIXED.X !!! ###
#            if(length(x.idx) > 0L) {
#                cat("\n\nDEBUG FIXED.X\n\n\n")
#                augUser <- user
#                idx <- which(augUser$fixed.x > 0L)
#                augUser$fixed.x[   idx ] <- 0L
#                augUser$free[      idx ] <- max(augUser$free) + 1:length(idx)
#                augUser$free.uncon[idx ] <- max(augUser$free.uncon) + 1:length(idx)
#                augModel <- Model(user           = augUser,
#                                  start          = getModelParameters(object, type="user"),
#                                  representation = object@representation,
#                                  debug          = FALSE)
#
#                Delta <- computeDelta(augModel)
#                E <- computeExpectedInformationMLM(object, sample=sample,
#                                                   Delta=Delta)
#                fixed.x.idx <- max(user$free) + 1:length(idx)
#                free.idx    <- 1:max(user$free)
#                E[free.idx, fixed.x.idx] <- 0.0
#                E[fixed.x.idx, free.idx] <- 0.0
#                E.inv <- solve(E)
#                x.idx <- integer(0)
#            }
            trace.UGamma <-
                testStatisticSatorraBentler.Mplus(sample = sample,
                                                  E.inv  = E.inv, 
                                                  Delta  = Delta, 
                                                  WLS.V  = WLS.V, 
                                                  Gamma  = Gamma,
                                                  x.idx  = x.idx)
        } else if(test == "satorra.bentler" && mimic != "Mplus") {
            trace.UGamma <- 
                testStatisticSatorraBentler(sample = sample,
                                            E.inv  = E.inv, 
                                            Delta  = Delta, 
                                            WLS.V  = WLS.V, 
                                            Gamma  = Gamma,
                                            x.idx  = x.idx)
        }

        scaling.factor <- sum(trace.UGamma) / df
        if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
        chisq.scaled         <- sum(chisq.group / scaling.factor)
        pvalue.scaled        <- 1 - pchisq(chisq.scaled, df)

        TEST[[2]] <- list(test=test,
                          stat=chisq.scaled,
                          stat.group=(chisq.group / scaling.factor),
                          df=df,
                          pvalue=pvalue.scaled,
                          scaling.factor=scaling.factor,
                          trace.UGamma=trace.UGamma)

    } else if(test == "yuan.bentler" && df > 0) {
        # try to extract attr from VCOV (if present)
        E.inv    <- attr(VCOV, "E.inv")
        B0.group <- attr(VCOV, "B0.group")

        if(is.null(E.inv)) {
            # if se="standard", information is probably expected
            # change it to observed
            if(options$se != "robust.mlr") information <- "observed"
            E.inv <- Nvcov.standard(object      = object,
                                    sample      = sample,
                                    estimator   = "ML",
                                    information = information)
        }

        if(mimic == "Mplus" || mimic == "lavaan") {
            if(is.null(B0.group)) {
                Nvcov <- Nvcov.first.order(object = object, sample = sample)
                B0.group <- attr(Nvcov, "B0.group")
            } 
            trace.UGamma <- 
                testStatisticYuanBentler.Mplus(sample=sample,
                                               information=information,
                                               B0.group=B0.group,
                                               E.inv=E.inv,
                                               x.idx=x.idx)
        } else {
            Delta <- computeDelta(object)
            Sigma.hat <- computeSigmaHat(object)
            if(object@meanstructure) Mu.hat <- computeMuHat(object)
            A1.group <- vector("list", length=sample@ngroups)
            B1.group <- vector("list", length=sample@ngroups)
            for(g in 1:sample@ngroups) {
                if(sample@missing.flag[g]) {
                    X <- NULL
                    M <- sample@missing[[g]]
                } else {
                    X <- sample@data.obs[[g]]
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
                testStatisticYuanBentler(sample=sample,
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

        ndat <- getNDAT(user)
        npar <- getNPAR(user)

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
        BOOT <- attr(VCOV, "BOOT")
        if(is.null(BOOT)) {
            if(!is.null(options$bootstrap)) {
                R <- options$bootstrap
            } else {
                R <- 1000L
            }
            BOOT <- bootstrapParameters.internal(model=object, sample=sample,
                        options=options, data=data, R=R, 
                        verbose=options$verbose)
        }    

        fx.group <- attr(BOOT, "fx.group")
        if(sample@ngroups > 1L) {
            boot.group <- t(apply(fx.group, 1, '*', NFAC))
            boot.group[which(boot.group < 0)] <- 0.0
            T.boot <- apply(boot.group, 1, sum)
        } else {
            boot.group <- fx.group * NFAC
            boot.group[which(boot.group < 0)] <- 0.0
            T.boot <- boot.group
        }

        # p-value
        boot.larger <- sum(T.boot > chisq)
        boot.length <- length(T.boot)
        pvalue.boot <- boot.larger/boot.length

        TEST[[2]] <- list(test="bootstrap",
                          stat=chisq,
                          stat.group=chisq.group,
                          df=df,
                          pvalue=pvalue.boot,
                          boot.T=T.boot,
                          boot.larger=boot.larger,
                          boot.length=boot.length)
    }

    TEST
}


