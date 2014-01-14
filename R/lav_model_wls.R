# compute WLS.est (as a list per group)
lav_model_wls_est <- function(object, GLIST = NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    ngroups       <- object@ngroups
    meanstructure <- object@meanstructure
    categorical   <- object@categorical
    group.w.free  <- object@group.w.free
    fixed.x       <- object@fixed.x
    num.idx       <- object@num.idx

    # compute moments for all groups
    Sigma.hat <- computeSigmaHat(object, GLIST=GLIST, extra=FALSE)
    if(meanstructure && !categorical) {
        Mu.hat <- computeMuHat(object, GLIST=GLIST)
    } else if(categorical) {
        TH <- computeTH(object, GLIST=GLIST)
        if(fixed.x)
            PI <- computePI(object, GLIST=GLIST)
    }
    if(group.w.free) {
        GW <- computeGW(object, GLIST=GLIST)
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
                wls.est <- c(TH[[g]],vec(PI[[g]]),
                             diag(Sigma.hat[[g]])[num.idx[[g]]],
                             vech(Sigma.hat[[g]], diagonal=FALSE))
            } else {
                wls.est <- c(TH[[g]],diag(Sigma.hat[[g]])[num.idx[[g]]],
                             vech(Sigma.hat[[g]], diagonal=FALSE))
            }
        } else if(meanstructure) {
            wls.est <- c(Mu.hat[[g]], vech(Sigma.hat[[g]]))
        } else {
            wls.est <- vech(Sigma.hat[[g]])
        }
        if(group.w.free) {
            #wls.est <- c(log(GW[[g]]/GW[[samplestats@ngroups]]), wls.est)
             wls.est <- c(GW[[g]], wls.est)
        }

        WLS.est[[g]] <- wls.est
    }

    WLS.est
}

