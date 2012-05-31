Nvcov.standard <- function(object, samplestats=NULL, data=NULL, estimator="ML", 
                           information="observed") {

    # compute information matrix
    if(information == "observed") {
        if(samplestats@missing.flag) {
            group.weight <- FALSE
        } else {
            group.weight <- TRUE
        }
        E <- computeObservedInformation(object, samplestats=samplestats, 
                                        data=data,
                                        estimator=estimator, type="free",
                                        group.weight=group.weight)
    } else {
        E <- computeExpectedInformation(object, samplestats=samplestats, 
                                        data=data,
                                        estimator=estimator)
    }

    E.inv <- solve(E)
    #E.inv <- MASS.ginv(E)
    if(nrow(object@con.jac) > 0L) {
        H <- object@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
        }
        if(nrow(H) > 0L) {
            NVarCov <- ( E.inv - E.inv %*% t(H) %*% 
                                 solve(H %*% E.inv %*% t(H)) %*%
                                 H %*% E.inv )
        } else {
            NVarCov <- E.inv
        }
    } else {
        NVarCov <- E.inv
    }

    NVarCov
}

Nvcov.bootstrap <- function(object, samplestats=NULL, options=NULL, data=NULL,
                            partable=NULL, control=list()) {

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

Nvcov.first.order <- function(object, samplestats=NULL, data=NULL) {

    B0.group <- vector("list", samplestats@ngroups)

    Sigma.hat <- computeSigmaHat(object); Mu.hat <- computeMuHat(object)
    Delta <- computeDelta(object)

    for(g in 1:samplestats@ngroups) {
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

    NVarCov <- solve( B0 )
    attr(NVarCov, "B0") <- B0
    attr(NVarCov, "B0.group") <- B0.group

    NVarCov
}

Nvcov.robust.mlm <- function(object, samplestats=NULL, data=NULL) {

    # compute information matrix
    E <- computeExpectedInformation(object, samplestats=samplestats, data=data,
                                    estimator="ML", extra=TRUE) 
    E.inv <- solve(E)
    Delta <- attr(E, "Delta")
    WLS.V <- attr(E, "WLS.V")

    # compute Gamma matrix (per group)
    Gamma <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        Gamma[[g]] <- compute.Gamma(data@X[[g]], meanstructure=TRUE)
    }

    NVarCov <- matrix(0, ncol=ncol(E), nrow=nrow(E))

    for(g in 1:samplestats@ngroups) {
        tmp <- WLS.V[[g]] %*% Delta[[g]] %*% E.inv
        NVarCov <- ( NVarCov +  (1/(samplestats@nobs[[g]]-1)) *
                              (t(tmp) %*% Gamma[[g]] %*% tmp) )
    } # g

    NVarCov <- NVarCov * (samplestats@ntotal-1)

    # to be reused by lavaanTest
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V
    attr(NVarCov, "Gamma") <- Gamma

    NVarCov
}


Nvcov.robust.mlm.mplus <- function(object, samplestats=NULL, data=NULL) {

    # compute information matrix with custom WLS.V
    E <- computeExpectedInformationMLM(object, samplestats=samplestats)
    E.inv <- solve(E)
    Delta <- attr(E, "Delta")
    WLS.V <- attr(E, "WLS.V")

    # compute Gamma matrix (per group)
    Gamma <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        Gamma[[g]] <- compute.Gamma(data@X[[g]], meanstructure=TRUE)
    }

    NVarCov <- matrix(0, ncol=ncol(E), nrow=nrow(E))

    for(g in 1:samplestats@ngroups) {
        # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
        # - WLS.V is not based on Sigma.hat and Mu.hat (as it
        #   should be?), but on samplestats@cov and samplestats@mean...
        # - Gamma is not identical to what is used for WLS; closer to EQS
        # - N/N-1 bug in G11 for NVarCov (but not test statistic)
        # - we divide by N-1! (just like EQS)

        tmp <- WLS.V[[g]] %*% Delta[[g]] %*% E.inv

        # should not be necessary, since already done!!
        # but otherwise, the std.errs for the intercepts are off...
        GAMMA <- Gamma[[g]]
        G11 <- ( samplestats@cov[[g]] * (samplestats@nobs[[g]]-1)/samplestats@nobs[[g]] ) 
        GAMMA[1:nrow(G11), 1:nrow(G11)] <- G11

        NVarCov <- ( NVarCov + (1/(samplestats@nobs[[g]]-1)) * 
                               (t(tmp) %*% GAMMA %*% tmp) )
    } # g

    NVarCov <- NVarCov * samplestats@ntotal

    # to be reused by lavaanTest
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V
    attr(NVarCov, "Gamma") <- Gamma

    NVarCov
}

Nvcov.robust.mlr <- function(object, samplestats=NULL, data=NULL,
                             information="observed") {

    # compute standard Nvcov
    E.inv <- Nvcov.standard(object      = object,
                            samplestats = samplestats,
                            data        = data,
                            estimator   = "ML",
                            information = information)

    # compute first.order Nvcov
    Nvcov <- Nvcov.first.order(object      = object, 
                               samplestats = samplestats,
                               data        = data)
    B0 <- attr(Nvcov, "B0")
    B0.group <- attr(Nvcov, "B0.group")

    # compute sandwich
    NVarCov <- E.inv %*% B0 %*% E.inv

    attr(NVarCov, "B0.group") <- B0.group
    attr(NVarCov, "E.inv") <- E.inv

    NVarCov
}

Nvcov.robust.wls <- function(object, samplestats=NULL) {

    # compute delta
    Delta <- computeDelta(object)

    NVarCov <- matrix(0, ncol=ncol(Delta[[1]]), nrow=ncol(Delta[[1]]))

    for(g in 1:samplestats@ngroups) {
        tmp <- Delta[[g]]
        acov <- samplestats@ACOV[[g]]
        NVarCov <- ( NVarCov +  (1/(samplestats@nobs[[g]]-1)) *
                              (t(tmp) %*% acov %*% tmp) )
    } # g

    NVarCov <- NVarCov * (samplestats@ntotal-1)

    # to be reused by lavaanTest
    attr(NVarCov, "Delta") <- Delta

    NVarCov
}
           

estimateVCOV <- function(object, samplestats, options=NULL, data=NULL, 
                         partable=NULL, control=list()) {

    estimator   <- options$estimator
    likelihood  <- options$likelihood
    information <- options$information
    se     <- options$se
    verbose     <- options$verbose
    mimic       <- options$mimic
  
    if(se == "none") return(NULL)

    # some require meanstructure (for now)
    if(se %in% c("first.order", "robust.mlm", "robust.mlr") && 
       !options$meanstructure) {
        stop("se (", se, ") requires meanstructure (for now)")
    }

    if(se == "standard") {
        NVarCov <- try( Nvcov.standard(object      = object, 
                                       samplestats = samplestats,
                                       data        = data,
                                       estimator   = estimator, 
                                       information = information) )

    } else if(se == "first.order") {
        NVarCov <- try( Nvcov.first.order(object      = object,
                                          samplestats = samplestats,
                                          data        = data) )

    } else if(se == "robust.mlm" && mimic == "Mplus") {
        NVarCov <- try( Nvcov.robust.mlm.mplus(object      = object,
                                               samplestats = samplestats,
                                               data        = data) )

    } else if(se == "robust.mlm" && mimic != "Mplus") {
        NVarCov <- try( Nvcov.robust.mlm(object      = object,
                                         samplestats = samplestats,
                                         data        = data) )

    } else if(se == "robust.mlr") {
        NVarCov <- try( Nvcov.robust.mlr(object      = object,
                                         samplestats = samplestats,
                                         data        = data,
                                         information = information) )

    } else if(se == "robust.wls") {
        NVarCov <- try( Nvcov.robust.wls(object      = object,
                                         samplestats = samplestats) )

    } else if(se == "bootstrap") {
        NVarCov <- try( Nvcov.bootstrap(object      = object,
                                        samplestats = samplestats,
                                        options     = options,
                                        data        = data,
                                        partable    = partable,
                                        control     = control) )
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


