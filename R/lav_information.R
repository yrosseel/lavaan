# model information
computeExpectedInformation <- function(lavmodel       = NULL, 
                                       lavsamplestats = NULL, 
                                       lavdata        = NULL,
                                       estimator      = "ML",
                                       # is no Delta is provided, we compute 
                                       # Delta for the free parameters only
                                       Delta = computeDelta(lavmodel=lavmodel), 
                                       lavcache       = NULL,
                                       extra          = FALSE) {

    # compute WLS.V
    WLS.V       <- vector("list", length=lavsamplestats@ngroups)
    if(estimator == "GLS"  || 
       estimator == "WLS"  || 
       estimator == "DWLS" ||
       estimator == "ULS") {
        # for GLS, the WLS.V22 part is: 0.5 * t(D) %*% [S.inv %x% S.inv] %*% D
        # for WLS, the WLS.V22 part is: Gamma
        WLS.V <- lavsamplestats@WLS.V
    } else if(estimator == "ML") {
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
        if(lavmodel@group.w.free) {
            GW <- unlist(computeGW(lavmodel = lavmodel))
        }
        if(lavsamplestats@missing.flag) {
            Mu.hat <- computeMuHat(lavmodel = lavmodel)
        }
        for(g in 1:lavsamplestats@ngroups) {
            if(lavsamplestats@missing.flag) {
                WLS.V[[g]] <- compute.Abeta(Sigma.hat=Sigma.hat[[g]],
                                            Mu.hat=Mu.hat[[g]],
                                            lavsamplestats=lavsamplestats, 
                                            lavdata=lavdata, group=g,
                                            information="expected")
            } else {
                # WLS.V22 = 0.5*t(D) %*% [Sigma.hat.inv %x% Sigma.hat.inv]%*% D
                WLS.V[[g]] <- 
                    compute.Abeta.complete(Sigma.hat=Sigma.hat[[g]],
                                           meanstructure=lavmodel@meanstructure)
            }
            if(lavmodel@group.w.free) {
                # unweight!!
                a <- exp(GW[g]) / lavsamplestats@nobs[[g]]
                # a <- exp(GW[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
                WLS.V[[g]] <- bdiag( matrix(a,1,1), WLS.V[[g]])
            }
        }
    } # ML


    # compute Information per group
    Info.group  <- vector("list", length=lavsamplestats@ngroups)
    for(g in 1:lavsamplestats@ngroups) {
        # note LISREL documentation suggest (Ng - 1) instead of Ng...
        fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
        # compute information for this group
        Info.group[[g]] <- fg * (t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]])
    }

    # assemble over groups
    Information <- Info.group[[1]]
    if(lavsamplestats@ngroups > 1) {
        for(g in 2:lavsamplestats@ngroups) {
            Information <- Information + Info.group[[g]]
        }
    }

    if(extra) {
        attr(Information, "Delta") <- Delta
        attr(Information, "WLS.V") <- WLS.V # unweighted
    }

    Information
}

# only for Mplus MLM
computeExpectedInformationMLM <- function(lavmodel = NULL, 
                                          lavsamplestats = NULL, 
                                          Delta = computeDelta(lavmodel = 
                                                               lavmodel)) {

    # compute WLS.V 
    WLS.V <- vector("list", length=lavsamplestats@ngroups)
    if(lavmodel@group.w.free) {
        GW <- unlist(computeGW(lavmodel = lavmodel))
    }
    for(g in 1:lavsamplestats@ngroups) {
        WLS.V[[g]] <- compute.A1.sample(lavsamplestats=lavsamplestats, group=g,
                                        meanstructure=TRUE, 
                                        information="expected")
        # the same as GLS... (except for the N/N-1 scaling)
        if(lavmodel@group.w.free) {
            # unweight!!
            a <- exp(GW[g]) / lavsamplestats@nobs[[g]]
            # a <- exp(GW[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
            WLS.V[[g]] <- bdiag( matrix(a,1,1), WLS.V[[g]])
        }
    }

    # compute Information per group
    Info.group  <- vector("list", length=lavsamplestats@ngroups)
    for(g in 1:lavsamplestats@ngroups) {
        fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
        # compute information for this group
        Info.group[[g]] <- fg * (t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]])
    }

    # assemble over groups
    Information <- Info.group[[1]]
    if(lavsamplestats@ngroups > 1) {
        for(g in 2:lavsamplestats@ngroups) {
            Information <- Information + Info.group[[g]]
        }
    }

    # always
    attr(Information, "Delta") <- Delta
    attr(Information, "WLS.V") <- WLS.V # unweighted

    Information
}





computeObservedInformation <- function(lavmodel       = NULL, 
                                       lavsamplestats = NULL, 
                                       lavdata        = NULL,
                                       type           = "free", 
                                       estimator      = "ML", 
                                       lavcache       = NULL,
                                       group.weight   = TRUE) {

    # computing the Richardson extrapolation
    # (note that this matrix is not fully symmetric --> do not use chol)
    Hessian <- matrix(0, lavmodel@nx.free, lavmodel@nx.free)
    x <- lav_model_get_parameters(lavmodel = lavmodel)
    for(j in 1:lavmodel@nx.free) {
        h.j <- 10e-6
        x.left <- x.left2 <- x.right <- x.right2 <- x
        x.left[j]  <- x[j] - h.j; x.left2[j]  <- x[j] - 2*h.j
        x.right[j] <- x[j] + h.j; x.right2[j] <- x[j] + 2*h.j

        g.left <- 
            lav_model_gradient(lavmodel       = lavmodel, 
                               GLIST          = lav_model_x2GLIST(lavmodel = 
                                                            lavmodel, x.left), 
                               lavsamplestats = lavsamplestats, 
                               lavdata        = lavdata, 
                               lavcache       = lavcache,
                               type           = "free", 
                               estimator      = estimator, 
                               group.weight   = group.weight)
        g.left2 <-    
            lav_model_gradient(lavmodel       = lavmodel,
                               GLIST          = lav_model_x2GLIST(lavmodel =
                                                            lavmodel, x.left2),
                               lavsamplestats = lavsamplestats, 
                               lavdata        = lavdata, 
                               lavcache       = lavcache,
                               type           = "free", 
                               estimator      = estimator,
                               group.weight   = group.weight)

        g.right <- 
            lav_model_gradient(lavmodel       = lavmodel,
                               GLIST          = lav_model_x2GLIST(lavmodel =
                                                            lavmodel, x.right),
                               lavsamplestats = lavsamplestats, 
                               lavdata        = lavdata, 
                               lavcache       = lavcache,
                               type           = "free", 
                               estimator      = estimator,
                               group.weight   = group.weight)

        g.right2 <- 
            lav_model_gradient(lavmodel       = lavmodel,
                               GLIST          = lav_model_x2GLIST(lavmodel =
                                                            lavmodel, x.right2),
                               lavsamplestats = lavsamplestats, 
                               lavdata        = lavdata, 
                               lavcache       = lavcache,
                               type           = "free", 
                               estimator      = estimator,
                               group.weight   = group.weight)
    
        Hessian[,j] <- ( -1 * (g.left2 - 8*g.left + 
                                         8*g.right - 
                               g.right2)/(12*h.j) )
    }

    #cat("Hessian 1:\n")
    #print(Hessian)

    #x <- lav_model_get_parameters(object, type="free")
    #compute.fx <- function(x) {
    #    GLIST <- lav_model_x2GLIST(object, x=x)
    #    fx <- lav_model_objective(lavmodel, GLIST=GLIST, 
    #                              lavlavsamplestats=lavlavsamplestats,
    #                              estimator=estimator,
    #                              link=link)
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
    # FIXME: what to do here if estimator = PML???
    if(type == "free.fixed.x" && lavmodel@fixed.x && 
       length(lavmodel@x.idx) > 0 && estimator != "PML") {
        idx.all <- which(object$free > 0 | (object$fixed.x > 0 &
                         object$row >= object$col))
        idx.free <- which(object$free > 0)
        idx.x <- which(object$fixed.x > 0 & object$row >= object$col)

        info.free.idx  <- idx.all %in% idx.free
        info.fixed.idx <- idx.all %in% idx.x

        Information.big <- matrix(0, nrow=length(idx.all), ncol=length(idx.all))
        Information.big[info.free.idx, info.free.idx] <- Information

        #### FIXME FOR MULTIPLE GROUPS !!!!!
        x.idx <- lavmodel@x.idx
        A1.group <- vector("list", lavsamplestats@ngroups)
        for(g in 1:lavsamplestats@ngroups) {
            A1.group[[g]] <- 
                compute.A1.sample(lavsamplestats=lavsamplestats, group=g,
                                  meanstructure=lavmodel@meanstructure,
                                  idx=x.idx)
        }
        if(lavmodel@multigroup) {
            # groups weights
            A1 <- (lavsamplestats@nobs[[1L]]/lavsamplestats@ntotal) * A1.group[[1L]]
            for(g in 2:lavsamplestats@ngroups) {
                A1 <- A1 + (lavsamplestats@nobs[[g]]/lavsamplestats@ntotal) * A1.group[[g]]
            }
        } else {
            A1 <- A1.group[[1L]]
        }

        Information.big[info.fixed.idx, info.fixed.idx] <- A1
        
        Information <- Information.big
    }

    Information
}
