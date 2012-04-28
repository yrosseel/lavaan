Nvcov.standard <- function(object, sample=NULL, data=NULL, estimator="ML", 
                           information="observed") {

    # compute information matrix
    if(information == "observed") {
        if(sample@missing.flag) {
            group.weight <- FALSE
        } else {
            group.weight <- TRUE
        }
        E <- computeObservedInformation(object, sample=sample, data=data,
                                        estimator=estimator, type="free",
                                        group.weight=group.weight)
    } else {
        E <- computeExpectedInformation(object, sample=sample, data=data,
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

Nvcov.bootstrap <- function(object, sample=NULL, options=NULL, data=NULL,
                            partable=NULL, control=list()) {

    # number of bootstrap draws
    if(!is.null(options$bootstrap.R)) {
        R <- options$bootstrap.R
    } else {
        R <- 1000L
    }
  
    boot.type <- "ordinary"
    if(options$test == "bollen.stine") boot.type <- "bollen.stine"

    TEST <- NULL
    COEF <- bootstrap.internal(object=NULL,
                               model.=object, sample.=sample, 
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
    NVarCov <- sample@ntotal * (cov(COEF) * (nboot-1)/nboot )

    # save COEF and TEST (if any)
    attr(NVarCov, "BOOT.COEF") <- COEF
    attr(NVarCov, "BOOT.TEST") <- TEST
 
    NVarCov
}

Nvcov.first.order <- function(object, sample=NULL, data=NULL) {

    B0.group <- vector("list", sample@ngroups)

    Sigma.hat <- computeSigmaHat(object); Mu.hat <- computeMuHat(object)
    Delta <- computeDelta(object)

    for(g in 1:sample@ngroups) {
        B1 <- compute.Bbeta(Sigma.hat=Sigma.hat[[g]], 
                            Mu.hat=Mu.hat[[g]],
                            sample=sample, data=data, group=g)

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

    if(sample@ngroups > 1L) {
        # groups weights
        B0 <- (sample@nobs[[1]]/sample@ntotal) * B0.group[[1]]
        for(g in 2:sample@ngroups) {
            B0 <- B0 + (sample@nobs[[g]]/sample@ntotal) * B0.group[[g]]
        }
    } else {
        B0 <- B0.group[[1]]
    }

    NVarCov <- solve( B0 )
    attr(NVarCov, "B0") <- B0
    attr(NVarCov, "B0.group") <- B0.group

    NVarCov
}

Nvcov.robust.mlm <- function(object, sample=NULL, data=NULL) {

    # compute information matrix
    E <- computeExpectedInformation(object, sample=sample, data=data,
                                    estimator="ML", extra=TRUE) 
    E.inv <- solve(E)
    Delta <- attr(E, "Delta")
    WLS.V <- attr(E, "WLS.V")

    # compute Gamma matrix (per group)
    Gamma <- vector("list", length=sample@ngroups)
    for(g in 1:sample@ngroups) {
        Gamma[[g]] <- compute.Gamma(data@X[[g]], meanstructure=TRUE)
    }

    NVarCov <- matrix(0, ncol=ncol(E), nrow=nrow(E))

    for(g in 1:sample@ngroups) {
        tmp <- WLS.V[[g]] %*% Delta[[g]] %*% E.inv
        NVarCov <- ( NVarCov +  (1/(sample@nobs[[g]]-1)) *
                              (t(tmp) %*% Gamma[[g]] %*% tmp) )
    } # g

    NVarCov <- NVarCov * (sample@ntotal-1)

    # to be reused by lavaanTest
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V
    attr(NVarCov, "Gamma") <- Gamma

    NVarCov
}


Nvcov.robust.mlm.mplus <- function(object, sample=NULL, data=NULL) {

    # compute information matrix with custom WLS.V
    E <- computeExpectedInformationMLM(object, sample=sample)
    E.inv <- solve(E)
    Delta <- attr(E, "Delta")
    WLS.V <- attr(E, "WLS.V")

    # compute Gamma matrix (per group)
    Gamma <- vector("list", length=sample@ngroups)
    for(g in 1:sample@ngroups) {
        Gamma[[g]] <- compute.Gamma(data@X[[g]], meanstructure=TRUE)
    }

    NVarCov <- matrix(0, ncol=ncol(E), nrow=nrow(E))

    for(g in 1:sample@ngroups) {
        # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
        # - WLS.V is not based on Sigma.hat and Mu.hat (as it
        #   should be?), but on sample@cov and sample@mean...
        # - Gamma is not identical to what is used for WLS; closer to EQS
        # - N/N-1 bug in G11 for NVarCov (but not test statistic)
        # - we divide by N-1! (just like EQS)

        tmp <- WLS.V[[g]] %*% Delta[[g]] %*% E.inv

        # should not be necessary, since already done!!
        # but otherwise, the std.errs for the intercepts are off...
        GAMMA <- Gamma[[g]]
        G11 <- ( sample@cov[[g]] * (sample@nobs[[g]]-1)/sample@nobs[[g]] ) 
        GAMMA[1:nrow(G11), 1:nrow(G11)] <- G11

        NVarCov <- ( NVarCov + (1/(sample@nobs[[g]]-1)) * 
                               (t(tmp) %*% GAMMA %*% tmp) )
    } # g

    NVarCov <- NVarCov * sample@ntotal

    # to be reused by lavaanTest
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "Delta") <- Delta
    attr(NVarCov, "WLS.V") <- WLS.V
    attr(NVarCov, "Gamma") <- Gamma

    NVarCov
}

Nvcov.robust.mlr <- function(object, sample=NULL, data=NULL,
                             information="observed") {

    # compute standard Nvcov
    E.inv <- Nvcov.standard(object      = object,
                            sample      = sample,
                            data        = data,
                            estimator   = "ML",
                            information = information)

    # compute first.order Nvcov
    Nvcov <- Nvcov.first.order(object = object, 
                               sample = sample,
                               data   = data)
    B0 <- attr(Nvcov, "B0")
    B0.group <- attr(Nvcov, "B0.group")

    # compute sandwich
    NVarCov <- E.inv %*% B0 %*% E.inv

    attr(NVarCov, "B0.group") <- B0.group
    attr(NVarCov, "E.inv") <- E.inv

    NVarCov
}
           

estimateVCOV <- function(object, sample, options=NULL, data=NULL, 
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
                                       sample      = sample,
                                       data        = data,
                                       estimator   = estimator, 
                                       information = information) )

    } else if(se == "first.order") {
        NVarCov <- try( Nvcov.first.order(object = object,
                                          sample = sample,
                                          data   = data) )

    } else if(se == "robust.mlm" && mimic == "Mplus") {
        NVarCov <- try( Nvcov.robust.mlm.mplus(object      = object,
                                               sample      = sample,
                                               data        = data) )

    } else if(se == "robust.mlm" && mimic != "Mplus") {
        NVarCov <- try( Nvcov.robust.mlm(object      = object,
                                         sample      = sample,
                                         data        = data) )

    } else if(se == "robust.mlr") {
        NVarCov <- try( Nvcov.robust.mlr(object      = object,
                                         sample      = sample,
                                         data        = data,
                                         information = information) )

    } else if(se == "bootstrap") {
        NVarCov <- try( Nvcov.bootstrap(object      = object,
                                        sample      = sample,
                                        options     = options,
                                        data        = data,
                                        partable    = partable,
                                        control     = control) )
    }

    if(! inherits(NVarCov, "try-error") ) {

        # denominator!
        if(estimator == "ML" && likelihood == "normal") {
            N <- sample@ntotal
        } else {
            N <- sample@ntotal - sample@ngroups
        }

        VarCov <- 1/N * NVarCov

    } else {
        warning("lavaan WARNING: could not compute standard errors!\n")
        VarCov <- NULL
    }

    VarCov
}


