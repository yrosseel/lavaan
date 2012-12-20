Nvcov.standard <- function(object, samplestats=NULL, data=NULL, estimator="ML", 
                           cache=NULL, information="observed") {

    # compute information matrix
    if(information == "observed") {
        if(samplestats@missing.flag) {
            group.weight <- FALSE
        } else {
            group.weight <- TRUE
        }
        E <- computeObservedInformation(object, samplestats=samplestats, 
                                        X=data@X, cache=cache,
                                        estimator=estimator, type="free",
                                        group.weight=group.weight)
    } else {
        E <- computeExpectedInformation(object, samplestats=samplestats, 
                                        data=data, cache=cache,
                                        estimator=estimator)
    }

    # handle constraints
    if(nrow(object@con.jac) > 0L) {
        H <- object@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- object@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(E), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            E3 <- rbind( cbind(     E,  H10, t(H)),
                         cbind(t(H10),   DL,  H0),
                         cbind(     H,   H0,  H0)  )
            NVarCov <- MASS::ginv(E3)[1:ncol(E), 1:ncol(E)]
            # FIXME: better include inactive + slacks??
        } else {
            NVarCov <- solve(E)
        }
    } else {
        NVarCov <- solve(E)
    }

    if(estimator == "PML") {
        NVarCov <- NVarCov * samplestats@ntotal
    }

    NVarCov
}

Nvcov.bootstrap <- function(object, samplestats=NULL, options=NULL, data=NULL,
                            cache=NULL, partable=NULL, control=list()) {

    # number of bootstrap draws
    if(!is.null(options$bootstrap)) {
        R <- options$bootstrap
    } else {
        R <- 1000L
    }
  
    boot.type <- "ordinary"
    if(options$test == "bollen.stine") boot.type <- "bollen.stine"

    TEST <- NULL
    COEF <- bootstrap.internal(object=NULL,
                               model.=object, samplestats.=samplestats, 
                               partable.=partable, 
                               options.=options, data.=data,
                               R=R, verbose=options$verbose,
                               type=boot.type,
                               FUN=ifelse(boot.type == "bollen.stine",
                                          "coeftest", "coef"),
                               warn=-1L,
                               parallel=control$parallel,
                               ncpus=control$ncpus,
                               cl=control$cl)
    if(boot.type == "bollen.stine") {
        nc <- ncol(COEF)
        TEST <- COEF[,nc]
        COEF <- COEF[,-nc]
    }

    # FIXME: cov rescale? Yes for now
    nboot <- nrow(COEF)
    NVarCov <- samplestats@ntotal * (cov(COEF) * (nboot-1)/nboot )

    # save COEF and TEST (if any)
    attr(NVarCov, "BOOT.COEF") <- COEF
    attr(NVarCov, "BOOT.TEST") <- TEST
 
    NVarCov
}

Nvcov.first.order <- function(object, samplestats=NULL, data=NULL,
                              cache=NULL, estimator="ML") {

    B0.group <- vector("list", samplestats@ngroups)

    if(estimator == "PML") {
        Delta <- computeDelta(object)
        Sigma.hat <- computeSigmaHat(object)
        TH <- computeTH(object)
    } else {
        Sigma.hat <- computeSigmaHat(object); Mu.hat <- computeMuHat(object)
        Delta <- computeDelta(object)
    }

    for(g in 1:samplestats@ngroups) {
        if(estimator == "PML") {
            # slow approach: compute outer product of case-wise scores
            SC <- pml_deriv1(Sigma.hat = Sigma.hat[[g]],
                             TH        = TH[[g]],
                             th.idx    = object@th.idx[[g]],
                             num.idx   = object@num.idx[[g]],
                             X         = data@X[[g]],
                             cache     = cache,
                             scores    = TRUE,
                             negative  = FALSE)

            # chain rule
            group.SC <- SC %*% Delta[[g]]

            # outer product
            B0.group[[g]] <- crossprod(group.SC)

            # to get NACOV instead of ACOV
            B0.group[[g]] <- B0.group[[g]] / samplestats@ntotal

            # group weights (if any)
            #group.dx <- group.w[g] * group.dx 
        } else {
            B1 <- compute.Bbeta(Sigma.hat=Sigma.hat[[g]], 
                                Mu.hat=Mu.hat[[g]],
                                samplestats=samplestats, data=data, group=g)

            B0.group[[g]] <- t(Delta[[g]]) %*% B1 %*% Delta[[g]] 

            #if(type=="free.fixed.x" && object@fixed.x && length(object@x.idx)>0) {
            #    # for resid(., type="standardized")
            #    idx.all <- which(object$free > 0 | (object$fixed.x > 0 &
            #                     object$row >= object$col))
            #    idx.free <- which(object$free > 0)
            #    idx.x <- which(object$fixed.x > 0 & 
            #                   object$row >= object$col)
            #        
            #    info.free.idx  <- idx.all %in% idx.free
            #    info.fixed.idx <- idx.all %in% idx.x
            #    B0.group[[g]][info.fixed.idx, info.free.idx] <- 0
            #    B0.group[[g]][info.free.idx, info.fixed.idx] <- 0
            #}
        }
    } # g

    if(samplestats@ngroups > 1L) {
        # groups weights
        B0 <- (samplestats@nobs[[1]]/samplestats@ntotal) * B0.group[[1]]
        for(g in 2:samplestats@ngroups) {
            B0 <- B0 + (samplestats@nobs[[g]]/samplestats@ntotal) * B0.group[[g]]
        }
    } else {
        B0 <- B0.group[[1]]
    }

    # handle constraints
    E <- B0
    if(nrow(object@con.jac) > 0L) {
        H <- object@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- object@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(E), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            E3 <- rbind( cbind(     E,  H10, t(H)),
                         cbind(t(H10),   DL,  H0),
                         cbind(     H,   H0,  H0)  )
            NVarCov <- MASS::ginv(E3)[1:ncol(E), 1:ncol(E)]
            # FIXME: better include inactive + slacks??
        } else {
            NVarCov <- solve(E)
        }
    } else {
        NVarCov <- solve(E)
    }

    attr(NVarCov, "B0") <- B0
    attr(NVarCov, "B0.group") <- B0.group

    NVarCov
}

Nvcov.robust.sem <- function(object, samplestats=NULL, data=NULL, cache=NULL,
                             estimator = "ML", mimic = "lavaan") {

    # compute information matrix
    if(estimator == "ML" && mimic == "Mplus") {
        # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
        # - WLS.V is not based on Sigma.hat and Mu.hat (as it
        #   should be?), but on samplestats@cov and samplestats@mean...
        # - Gamma is not identical to what is used for WLS; closer to EQS
        # - N/N-1 bug in G11 for NVarCov (but not test statistic)
        # - we divide by N-1! (just like EQS)
        E <- computeExpectedInformationMLM(object, samplestats=samplestats)
        Gamma <- samplestats@NACOV
        # 'fix' G11 part of Gamma (NOTE: this is NOT needed for SB test 
        # statistic
        for(g in 1:samplestats@ngroups) {
            gg1 <- (samplestats@nobs[[g]]-1)/samplestats@nobs[[g]]
            G11 <- gg1 * samplestats@cov[[g]]
            Gamma[[g]][1:nrow(G11), 1:nrow(G11)] <- G11
        } # g
    } else {
        E <- computeExpectedInformation(object, samplestats=samplestats, 
                                        data=data, estimator=estimator, 
                                        extra=TRUE) 
        Gamma <- samplestats@NACOV
    }

    # handle constraints
    if(nrow(object@con.jac) > 0L) {
        H <- object@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- object@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(E), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            E3 <- rbind( cbind(     E,  H10, t(H)),
                         cbind(t(H10),   DL,  H0),
                         cbind(     H,   H0,  H0)  )
            E.inv <- MASS::ginv(E3)[1:ncol(E), 1:ncol(E)]
            # FIXME: better include inactive + slacks??
        } else {
            E.inv <- solve(E)
        }
    } else {
        E.inv <- solve(E)
    }

    Delta <- attr(E, "Delta")
    WLS.V <- attr(E, "WLS.V")

    tDVGVD <- matrix(0, ncol=ncol(E), nrow=nrow(E))
    for(g in 1:samplestats@ngroups) {
        fg  <-  samplestats@nobs[[g]]   /samplestats@ntotal
        fg1 <- (samplestats@nobs[[g]]-1)/samplestats@ntotal
        # fg twice for WLS.V, 1/fg1 once for GaMMA
        # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
        # t(Delta) * WLS.V %*% Gamma %*% WLS.V %*% Delta
        WD <- WLS.V[[g]] %*% Delta[[g]]
        tDVGVD <- tDVGVD + fg*fg/fg1 * crossprod(WD, Gamma[[g]] %*% WD)
    } # g
    NVarCov <- (E.inv %*% tDVGVD %*% E.inv)

    # to be reused by lavaanTest
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V

    NVarCov
}

Nvcov.robust.huber.white <- function(object, samplestats=NULL, data=NULL,
                                     information="observed", cache=NULL, 
                                     estimator="ML") {

    # compute standard Nvcov
    E.inv <- Nvcov.standard(object      = object,
                            samplestats = samplestats,
                            data        = data,
                            estimator   = estimator,
                            cache       = cache,
                            information = information)

    # compute first.order Nvcov
    Nvcov <- Nvcov.first.order(object      = object, 
                               samplestats = samplestats,
                               data        = data,
                               cache       = cache,
                               estimator   = estimator)
    B0 <- attr(Nvcov, "B0")
    B0.group <- attr(Nvcov, "B0.group")

    # compute sandwich
    NVarCov <- E.inv %*% B0 %*% E.inv

    attr(NVarCov, "B0.group") <- B0.group
    attr(NVarCov, "E.inv") <- E.inv

    NVarCov
}


estimateVCOV <- function(object, samplestats, options=NULL, data=NULL, 
                         partable=NULL, cache=NULL, control=list()) {

    estimator   <- options$estimator
    likelihood  <- options$likelihood
    information <- options$information
    se     <- options$se
    verbose     <- options$verbose
    mimic       <- options$mimic
  
    if(se == "none") return(NULL)

    # some require meanstructure (for now)
    if(se %in% c("first.order", "robust.sem", "robust.huber.white") && 
       !options$meanstructure) {
        stop("se (", se, ") requires meanstructure (for now)")
    }

    if(se == "standard") {
        NVarCov <- try( Nvcov.standard(object      = object, 
                                       samplestats = samplestats,
                                       data        = data,
                                       estimator   = estimator, 
                                       cache       = cache,
                                       information = information),
                         silent=TRUE )

    } else if(se == "first.order") {
        NVarCov <- try( Nvcov.first.order(object      = object,
                                          samplestats = samplestats,
                                          data        = data,
                                          cache       = cache,
                                          estimator   = estimator),
                         silent=TRUE )

    } else if(se == "robust.sem") {
        NVarCov <- try( Nvcov.robust.sem(object      = object,
                                         samplestats = samplestats,
                                         estimator   = estimator,
                                         mimic       = mimic,
                                         cache       = cache,
                                         data        = data),
                        silent=TRUE )

    } else if(se == "robust.huber.white") {
        NVarCov <- try( Nvcov.robust.huber.white(object      = object,
                                         samplestats = samplestats,
                                         data        = data,
                                         information = information,
                                         cache       = cache,
                                         estimator   = estimator),
                        silent=TRUE )

    } else if(se == "bootstrap") {
        NVarCov <- try( Nvcov.bootstrap(object      = object,
                                        samplestats = samplestats,
                                        options     = options,
                                        data        = data,
                                        cache       = cache,
                                        partable    = partable,
                                        control     = control),
                        silent=TRUE )
    }

    if(! inherits(NVarCov, "try-error") ) {

        # denominator!
        if(estimator == "ML" && likelihood == "normal") {
            N <- samplestats@ntotal
        } else {
            N <- samplestats@ntotal - samplestats@ngroups
        }

        VarCov <- 1/N * NVarCov

    } else {
        warning("lavaan WARNING: could not compute standard errors!\n")
        VarCov <- NULL
    }

    VarCov
}


