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

    # compute/get WLS.V
    # if DWLS or ULS, this is the diagonal only! (since 0.5-17)
    WLS.V <- lav_model_wls_v(lavmodel       = lavmodel, 
                             lavsamplestats = lavsamplestats,
                             estimator      = estimator,
                             lavdata        = lavdata)

    # compute Information per group
    Info.group  <- vector("list", length=lavsamplestats@ngroups)
    for(g in 1:lavsamplestats@ngroups) {
        # note LISREL documentation suggest (Ng - 1) instead of Ng...
        fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
        # compute information for this group
        if(estimator %in% c("DWLS", "ULS")) {
            # diagonal weight matrix
            Delta2 <- sqrt(WLS.V[[g]]) * Delta[[g]]
            Info.group[[g]] <- fg * crossprod(Delta2)
        } else {
            # full weight matrix
            # Info.group[[g]] <- 
                # fg * (t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]])
            Info.group[[g]] <- 
                fg * ( crossprod(Delta[[g]], WLS.V[[g]]) %*% Delta[[g]] )
        }
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
                                       estimator      = "ML",
                                       lavcache       = NULL,
                                       group.weight   = TRUE) {

    Hessian <- lav_model_hessian(lavmodel       = lavmodel,
                                 lavsamplestats = lavsamplestats,
                                 lavdata        = lavdata,
                                 lavcache       = lavcache,
                                 estimator      = estimator,
                                 group.weight   = TRUE)

    # in lavaan, we ALWAYS minimize, so the Hessian is already pos def
    Information <- Hessian

    Information
}
