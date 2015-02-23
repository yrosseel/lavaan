testStatisticSatorraBentler <- function(lavsamplestats=lavsamplestats, 
                                        E.inv, Delta, WLS.V,
                                        x.idx=list(integer(0)),
                                        test.UGamma.eigvals = TRUE) {

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

    Gamma <- lavsamplestats@NACOV
    nss <- ncol(Gamma[[1]])
    ngroups <- lavsamplestats@ngroups

    UG.group  <- vector("list", length=ngroups)
    trace.UGamma2 <- numeric(ngroups)
    for(g in 1:ngroups) {
        fg  <-  lavsamplestats@nobs[[g]]   /lavsamplestats@ntotal
        fg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@ntotal
        Gamma.g <- Gamma[[g]] / fg  ## ?? check this
        Delta.g <- Delta[[g]]

        # just for testing: regular UG:
        # WLS.Vg  <- WLS.V[[g]] * fg
        # UG1 <- Gamma.g %*% (WLS.Vg - WLS.Vg %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]]) %*% WLS.Vg)

        # diagonal WLS.V? we check for this since 0.5-17
        diagonal <- FALSE
        if(is.matrix(WLS.V[[g]])) {
            A1 <- WLS.V[[g]] * fg
            B1 <- A1 %*% Gamma.g %*% A1
        } else {
            diagonal <- TRUE
            a1 <- WLS.V[[g]] * fg # numeric vector!
            B1 <- Gamma.g * tcrossprod(a1)
        }
        # mask independent 'fixed-x' variables
        # note: this only affects the saturated H1 model
        if(length(x.idx[[g]]) > 0L) {
            nvar <- ncol(lavsamplestats@cov[[g]])
            idx <- eliminate.pstar.idx(nvar=nvar, el.idx=x.idx[[g]],
                                       meanstructure=TRUE, type="all")
            if(diagonal) {
                a1 <- a1[idx]
            } else {
                A1 <- A1[idx,idx]
            }
            B1      <- B1[idx,idx]
            Delta.g <- Delta.g[idx,]
        } 


        if(diagonal) {
            tmp <- t((1/a1) * B1) - (B1 %*% Delta.g %*% tcrossprod(E.inv, Delta.g))
        } else {
            A1.inv <- solve(A1)
            tmp <- (B1 %*% A1.inv) - (B1 %*% Delta.g %*% tcrossprod(E.inv, Delta.g))
        }
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
    # tUG <- t(UG); trace.UGamma3 <- sum(UG * tUG) # seems wrong?
    # attr(trace.UGamma, "trace.UGamma3") <- trace.UGamma3

    # eigen values
    # this is for the lavaan.survey pval.pFsum() function
    # but for large problems, this can take a loooong time; perhaps
    # we should make this optional?
    if(test.UGamma.eigvals) {
        attr(trace.UGamma, "eigenvalues") <-  
            Re(eigen(UG, only.values=TRUE)$values)
    }

    trace.UGamma
}

testStatisticYuanBentler <- function(lavsamplestats=lavsamplestats,
                                     A1.group=NULL,
                                     B1.group=NULL,
                                     Delta=NULL,
                                     E.inv=NULL,
                                     x.idx=list(integer(0))) {

    # we always assume a meanstructure
    meanstructure <- TRUE

    trace.UGamma <- numeric( lavsamplestats@ngroups )
    trace.h1     <- numeric( lavsamplestats@ngroups )
    trace.h0     <- numeric( lavsamplestats@ngroups )

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
        trace.h0[g] <- sum( (B1 %*% Delta[[g]] %*% E.inv %*% t(Delta[[g]])) )
        trace.UGamma[g] <- abs(trace.h1[g] - trace.h0[g])
    }

    # traces
    trace.UGamma <- sum(trace.UGamma)
    #tUG <- t(UG); trace.UGamma2 <- sum(UG * tUG)
    #attr(trace.UGamma, "trace.UGamma2") <- trace.UGamma2

    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0

    trace.UGamma
}

testStatisticYuanBentler.Mplus <- function(lavsamplestats=lavsamplestats, 
                                           lavdata=lavdata,
                                           information="observed",
                                           B0.group=NULL,
                                           E.inv=NULL,
                                           x.idx=list(integer(0))) {
    # typical for Mplus:
    # - do NOT use the YB formula, but use an approximation
    #   relying  on A0 ~= Delta'*A1*Delta and the same for B0

    # we always assume a meanstructure
    meanstructure <- TRUE
    ngroups <- lavsamplestats@ngroups

    trace.UGamma <- numeric( lavsamplestats@ngroups )
    trace.h1     <- numeric( lavsamplestats@ngroups )
    trace.h0     <- numeric( lavsamplestats@ngroups )

    for(g in 1:lavsamplestats@ngroups) {
        # if lavdata is complete, A1.22 is simply 0.5*D'(S.inv x S.inv)D
        A1 <- compute.A1.sample(lavsamplestats=lavsamplestats, group=g,
                                meanstructure=meanstructure,
                                information=information)
        B1 <- compute.B1.sample(lavsamplestats=lavsamplestats, 
                                lavdata=lavdata, group=g)

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
                           test.UGamma.eigvals = TRUE,
                           control        = list()) {


    mimic       <- lavoptions$mimic
    test        <- lavoptions$test
    information <- lavoptions$information
    estimator   <- lavoptions$estimator


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
        #### FIXME!!!! does not work for multiple groups yet!!!
        chisq.group <- PML$PLRTH0Sat
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

    # some require meanstructure (for now)
    if(test %in% c("satorra.bentler", "yuan.bentler") &&
       !lavoptions$meanstructure) {
        stop("test (", test, ") requires meanstructure (for now)")
    }

    # fixed.x idx
    x.idx <- vector("list", length=lavsamplestats@ngroups)
    for(g in 1:lavsamplestats@ngroups) {
        if(lavoptions$fixed.x && estimator == "ML") {
            x.idx[[g]] <- match(vnames(lavpartable, "ov.x", group=g), 
                                lavdata@ov.names[[g]])
        } else {
            x.idx[[g]] <- integer(0L)
        }
    }

    if(lavoptions$estimator == "PML") {
        if(test == "standard") {
            # nothing to do
        } else if(test == "mean.var.adjusted") {
            TEST[[2]] <- list(test               = test,
                              stat               = PML$stat,
                              stat.group         = as.numeric(NA), # for now
                              df                 = PML$df,
                              pvalue             = PML$p.value,
                              scaling.factor     = PML$scaling.factor,
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
        # try to extract attr from VCOV (if present)
        E.inv <- attr(VCOV, "E.inv")
        Delta <- attr(VCOV, "Delta")
        WLS.V <- attr(VCOV, "WLS.V")

        # if not present (perhaps se.type="standard" or se.type="none")
        #  we need to compute these again
        if(is.null(E.inv) || is.null(Delta) || is.null(WLS.V)) {
            if(mimic == "Mplus" && estimator == "ML") {
                # special treatment for Mplus
                E <- lav_model_information_expected_MLM(lavmodel = lavmodel,
                         lavsamplestats=lavsamplestats)
            } else {
                E <- lav_model_information_expected(lavmodel = lavmodel, 
                         lavsamplestats = lavsamplestats, lavdata = lavdata,
                         estimator = estimator, extra = TRUE)
            }
            E.inv <- try(solve(E), silent=TRUE)
            if(inherits(E.inv, "try-error")) {
                TEST[[2]] <- list(test=test, stat=as.numeric(NA), 
                    stat.group=rep(as.numeric(NA), lavsamplestats@ngroups),
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
#                augModel <- lav_model(lavpartable       = augUser,
#                                  start          = lav_model_get_parameters(object, type="user"),
#                                  representation = lavmodel@representation,
#                                  link           = lavmodel@link,
#                                  debug          = FALSE)
#
#                Delta <- computeDelta(lavmodel = augModel)
#                E <- lav_model_information_expected_MLM(object, lavsamplestats=lavsamplestats,
#                                                   Delta=Delta)
#                fixed.x.idx <- max(lavpartable$free) + 1:length(idx)
#                free.idx    <- 1:max(lavpartable$free)
#                E[free.idx, fixed.x.idx] <- 0.0
#                E[fixed.x.idx, free.idx] <- 0.0
#                E.inv <- solve(E)
#                x.idx <- integer(0)
#            }
#        } 

        trace.UGamma <- 
            testStatisticSatorraBentler(lavsamplestats = lavsamplestats,
                                        E.inv       = E.inv, 
                                        Delta       = Delta, 
                                        WLS.V       = WLS.V, 
                                        x.idx       = x.idx)
        trace.UGamma2 <- attr(trace.UGamma, "trace.UGamma2")
        # trace.UGamma3 <- attr(trace.UGamma, "trace.UGamma3")
        trace.UGamma4 <- attr(trace.UGamma, "trace.UGamma4")
        if(test.UGamma.eigvals) {
            UGamma.eigenvalues <- attr(trace.UGamma, "eigenvalues")
        } else {
            UGamma.eigenvalues <- numeric(0L)
        }
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
            fg <- unlist(lavsamplestats@nobs)/lavsamplestats@ntotal
           
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

        TEST[[2]] <- list(test               = test,
                          stat               = chisq.scaled,
                          stat.group         = stat.group,
                          df                 = df,
                          pvalue             = pvalue.scaled,
                          scaling.factor     = scaling.factor,
                          shift.parameter    = shift.parameter,
                          trace.UGamma       = trace.UGamma,
                          trace.UGamma4      = trace.UGamma4,
                          trace.UGamma2      = trace.UGamma2,
                          UGamma.eigenvalues = UGamma.eigenvalues)

    } else if(test == "yuan.bentler" && df > 0) {
        # try to extract attr from VCOV (if present)
        E.inv    <- attr(VCOV, "E.inv")
        B0.group <- attr(VCOV, "B0.group")

        if(is.null(E.inv)) {
            # if se="standard", information is probably expected
            # change it to observed
            if(lavoptions$se != "robust.mlr") information <- "observed"
            E.inv <- lav_model_information(lavmodel       = lavmodel,
                                           lavsamplestats = lavsamplestats,
                                           lavdata        = lavdata,
                                           estimator      = "ML",
                                           lavcache       = lavcache,
                                           information    = information,
                                           extra          = FALSE,
                                           augmented      = TRUE,
                                           inverted       = TRUE)
        }

        if(mimic == "Mplus" || mimic == "lavaan") {
            if(is.null(B0.group)) {
                B0 <- 
                    lav_model_information_firstorder(lavmodel = lavmodel,
                                             lavsamplestats = lavsamplestats,
                                             lavdata        = lavdata,
                                             estimator      = estimator,
                                             lavcache       = lavcache,
                                             check.pd       = FALSE,
                                             augmented      = FALSE,
                                             inverted       = FALSE)
                B0.group <- attr(B0, "B0.group")
            }
            trace.UGamma <- 
                testStatisticYuanBentler.Mplus(lavsamplestats = lavsamplestats,
                                               lavdata        = lavdata,
                                               information = information,
                                               B0.group    = B0.group,
                                               E.inv       = E.inv,
                                               x.idx       = x.idx)
        } else {
            Delta <- computeDelta(lavmodel = lavmodel)
            Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
            if(lavmodel@meanstructure) {
                Mu.hat <- computeMuHat(lavmodel = lavmodel)
            }
            A1.group <- vector("list", length=lavsamplestats@ngroups)
            B1.group <- vector("list", length=lavsamplestats@ngroups)
            for(g in 1:lavsamplestats@ngroups) {
                if(lavsamplestats@missing.flag) {
                    X <- NULL
                    M <- lavsamplestats@missing[[g]]
                } else {
                    X <- lavdata@X[[g]]
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
                testStatisticYuanBentler(lavsamplestats=lavsamplestats,
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

        ndat <- lav_partable_ndat(lavpartable)
        npar <- lav_partable_npar(lavpartable)

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
                                   warn            = -1L,
                                   parallel        = control$parallel,
                                   ncpus           = control$ncpus,
                                   cl              = control$cl)
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


