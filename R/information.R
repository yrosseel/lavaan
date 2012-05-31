
computeExpectedInformation <- function(object, samplestats=NULL, data=NULL,
                                       estimator="ML",
                                       # is no Delta is provided, we compute 
                                       # Delta for the free parameters only
                                       Delta=computeDelta(object), 
                                       extra=FALSE) {

    # compute WLS.V
    WLS.V       <- vector("list", length=samplestats@ngroups)
    if(estimator == "GLS"  || 
       estimator == "WLS"  || 
       estimator == "DWLS" ||
       estimator == "ULS") {
        # for GLS, the WLS.V22 part is: 0.5 * t(D) %*% [S.inv %x% S.inv] %*% D
        # for WLS, the WLS.V22 part is: Gamma
        WLS.V <- samplestats@WLS.V
    } else if(estimator == "ML") {
        Sigma.hat <- computeSigmaHat(object)
        if(object@meanstructure) Mu.hat <- computeMuHat(object)
        for(g in 1:samplestats@ngroups) {
            if(samplestats@missing.flag) {
                WLS.V[[g]] <- compute.Abeta(Sigma.hat=Sigma.hat[[g]],
                                            Mu.hat=Mu.hat[[g]],
                                            samplestats=samplestats, 
                                            data=data, group=g,
                                            information="expected")
            } else {
                # WLS.V22 = 0.5*t(D) %*% [Sigma.hat.inv %x% Sigma.hat.inv]%*% D
                WLS.V[[g]] <- 
                    compute.Abeta.complete(Sigma.hat=Sigma.hat[[g]],
                                           meanstructure=object@meanstructure)
            }
        }
    } # ML


    # compute Information per group
    Info.group  <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        # take care of multiple groups
        WLS.V[[g]] <- samplestats@nobs[[g]]/samplestats@ntotal * WLS.V[[g]]

        # compute information for this group
        Info.group[[g]] <- t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]] 
    }

    # assemble over groups
    Information <- Info.group[[1]]
    if(samplestats@ngroups > 1) {
        for(g in 2:samplestats@ngroups) {
            Information <- Information + Info.group[[g]]
        }
    }

    if(extra) {
        attr(Information, "Delta") <- Delta
        attr(Information, "WLS.V") <- WLS.V
    }

    Information
}

# only for Mplus MLM
computeExpectedInformationMLM <- function(object, samplestats = NULL, 
                                          Delta = computeDelta(object)) {

    # compute WLS.V 
    WLS.V <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        WLS.V[[g]] <- compute.A1.sample(samplestats=samplestats, group=g,
                                        meanstructure=TRUE)
        # the same as GLS... (except for the N/N-1 scaling)
    }

    # compute Information per group
    Info.group  <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        # take care of multiple groups
        WLS.V[[g]] <- samplestats@nobs[[g]]/samplestats@ntotal * WLS.V[[g]]

        # compute information for this group
        Info.group[[g]] <- t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]]

    }

    # assemble over groups
    Information <- Info.group[[1]]
    if(samplestats@ngroups > 1) {
        for(g in 2:samplestats@ngroups) {
            Information <- Information + Info.group[[g]]
        }
    }

    # always
    attr(Information, "Delta") <- Delta
    attr(Information, "WLS.V") <- WLS.V

    Information
}





computeObservedInformation <- function(object, samplestats=NULL, data=NULL,
                                       type="free", estimator="ML", 
                                       group.weight=TRUE) {

    # computing the Richardson extrapolation
    # (note that this matrix is not fully symmetric --> do not use chol)
    Hessian <- matrix(0, object@nx.free, object@nx.free)
    x <- getModelParameters(object)
    for(j in 1:object@nx.free) {
        h.j <- 10e-6
        x.left <- x.left2 <- x.right <- x.right2 <- x
        x.left[j]  <- x[j] - h.j; x.left2[j]  <- x[j] - 2*h.j
        x.right[j] <- x[j] + h.j; x.right2[j] <- x[j] + 2*h.j

        g.left <- 
            computeGradient(object=object, GLIST=x2GLIST(object, x.left), 
                            samplestats=samplestats, type="free", estimator=estimator, 
                            group.weight=group.weight)
        g.left2 <-    
            computeGradient(object=object, GLIST=x2GLIST(object, x.left2),
                            samplestats=samplestats, type="free", estimator=estimator, 
                            group.weight=group.weight)

        g.right <- 
            computeGradient(object=object, GLIST=x2GLIST(object, x.right),
                            samplestats=samplestats, type="free", estimator=estimator,
                            group.weight=group.weight)

        g.right2 <- 
            computeGradient(object=object, GLIST=x2GLIST(object, x.right2),
                            samplestats=samplestats, type="free", estimator=estimator,
                            group.weight=group.weight)
    
        Hessian[,j] <- ( -1 * (g.left2 - 8*g.left + 
                                         8*g.right - 
                               g.right2)/(12*h.j) )
    }

    #cat("Hessian 1:\n")
    #print(Hessian)

    #x <- getModelParameters(object, type="free")
    #compute.fx <- function(x) {
    #    GLIST <- x2GLIST(object, x=x)
    #    fx <- computeObjective(object, GLIST=GLIST, samplestats=samplestats,
    #                           estimator=estimator)
    #    fx
    #}
    #Hessian <- - numDeriv::hessian(func=compute.fx, x=x)
    #cat("Hessian 2:\n")
    #print(Hessian)

    #stop("for now")


    # make symmetric (NEEDED? probably not)
    Hessian <- ( Hessian + t(Hessian) )/2.0
    Information <- ( -1 * Hessian )

    # augmented Information matrix (only for fixed.x)
    if(type == "free.fixed.x" && object@fixed.x && length(object@x.idx) > 0) {
        idx.all <- which(object$free > 0 | (object$fixed.x > 0 &
                         object$row >= object$col))
        idx.free <- which(object$free > 0)
        idx.x <- which(object$fixed.x > 0 & object$row >= object$col)

        info.free.idx  <- idx.all %in% idx.free
        info.fixed.idx <- idx.all %in% idx.x

        Information.big <- matrix(0, nrow=length(idx.all), ncol=length(idx.all))
        Information.big[info.free.idx, info.free.idx] <- Information

        #### FIXME FOR MULTIPLE GROUPS !!!!!
        x.idx <- object@x.idx
        A1.group <- vector("list", samplestats@ngroups)
        for(g in 1:samplestats@ngroups) {
            A1.group[[g]] <- 
                compute.A1.sample(samplestats=samplestats, group=g,
                                  meanstructure=object@meanstructure,
                                  idx=x.idx)
        }
        if(object@multigroup) {
            # groups weights
            A1 <- (samplestats@nobs[[1L]]/samplestats@ntotal) * A1.group[[1L]]
            for(g in 2:samplestats@ngroups) {
                A1 <- A1 + (samplestats@nobs[[g]]/samplestats@ntotal) * A1.group[[g]]
            }
        } else {
            A1 <- A1.group[[1L]]
        }

        Information.big[info.fixed.idx, info.fixed.idx] <- A1
        
        Information <- Information.big
    }

    Information
}
