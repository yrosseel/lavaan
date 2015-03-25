# compute WLS.est (as a list per group)
lav_model_wls_est <- function(lavmodel = NULL, GLIST = NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    ngroups       <- lavmodel@ngroups
    meanstructure <- lavmodel@meanstructure
    categorical   <- lavmodel@categorical
    group.w.free  <- lavmodel@group.w.free
    fixed.x       <- lavmodel@fixed.x
    num.idx       <- lavmodel@num.idx

    # compute moments for all groups
    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, GLIST = GLIST,
                                 extra = FALSE)
    if(meanstructure && !categorical) {
        Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
    } else if(categorical) {
        TH <- computeTH(lavmodel = lavmodel, GLIST = GLIST)
        if(fixed.x)
            PI <- computePI(lavmodel = lavmodel, GLIST = GLIST)
    }
    if(group.w.free) {
        GW <- computeGW(lavmodel = lavmodel, GLIST = GLIST)
    }

    WLS.est <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        if(categorical) {
            # order of elements is important here:
            # 1. thresholds + means (interleaved)
            # 2. slopes (if any, columnwise per exo)
            # 3. variances (if any)
            # 4. correlations (no diagonal!)
            if(fixed.x) {
                wls.est <- c(TH[[g]],lav_matrix_vec(PI[[g]]),
                             diag(Sigma.hat[[g]])[num.idx[[g]]],
                             lav_matrix_vech(Sigma.hat[[g]], diagonal=FALSE))
            } else {
                wls.est <- c(TH[[g]],diag(Sigma.hat[[g]])[num.idx[[g]]],
                             lav_matrix_vech(Sigma.hat[[g]], diagonal=FALSE))
            }
        } else if(meanstructure) {
            wls.est <- c(Mu.hat[[g]], lav_matrix_vech(Sigma.hat[[g]]))
        } else {
            wls.est <- lav_matrix_vech(Sigma.hat[[g]])
        }
        if(group.w.free) {
            #wls.est <- c(log(GW[[g]]/GW[[samplestats@ngroups]]), wls.est)
             wls.est <- c(GW[[g]], wls.est)
        }

        WLS.est[[g]] <- wls.est
    }

    WLS.est
}

# compute WLS.V (as a list per group)
lav_model_wls_v <- function(lavmodel       = NULL, 
                            lavsamplestats = NULL,
                            estimator      = "ML",
                            lavdata        = NULL) {

    WLS.V <- vector("list", length=lavsamplestats@ngroups)

    # if we are using *LS, we already have WLS.V
    if(estimator == "GLS"  || estimator == "WLS") {
        # for GLS, the WLS.V22 part is: 0.5 * t(D) %*% [S.inv %x% S.inv] %*% D
        # for WLS, the WLS.V22 part is: Gamma
        WLS.V <- lavsamplestats@WLS.V
    } else if(estimator == "DWLS" || estimator == "ULS") {
        # diagonal only!!
        WLS.V <- lavsamplestats@WLS.VD


    # for ML, we need to recompute this, as it is function of Sigma (and Mu)
    } else if(estimator == "ML") {
        WLS.V <- vector("list", length=lavsamplestats@ngroups)

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

    WLS.V
}

