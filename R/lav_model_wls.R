# compute WLS.est (as a list per group)
lav_model_wls_est <- function(lavmodel = NULL, GLIST = NULL) { #,
                              #cov.x = NULL) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    nblocks       <- lavmodel@nblocks
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
    }
    if(lavmodel@conditional.x && lavmodel@nexo > 0L) {
        PI <- computePI(lavmodel = lavmodel, GLIST = GLIST)
    } 
    if(group.w.free) {
        GW <- computeGW(lavmodel = lavmodel, GLIST = GLIST)
    }

    WLS.est <- vector("list", length=nblocks)
    for(g in 1:nblocks) {

        # PI?
        if(lavmodel@conditional.x && lavmodel@nexo > 0L) {
            PI.g <- PI[[g]]
            # Sigma_yy = Sigma_yy|x + PI %*% cov.x %*% t(PI)
            #if(!categorical) {
            #    cov.xg <- cov.x[[g]]
            #    Sigma.hat[[g]] <- Sigma.hat[[g]] + PI.g %*% cov.xg %*% t(PI.g)
            #}
        } else {
            PI.g <- numeric(0L)
        }

        if(categorical) {
            # order of elements is important here:
            # 1. thresholds + means (interleaved)
            # 2. slopes (if any, columnwise per exo)
            # 3. variances (if any)
            # 4. correlations (no diagonal!)
            wls.est <- c(TH[[g]],
                         lav_matrix_vec(PI.g),
                         diag(Sigma.hat[[g]])[num.idx[[g]]],
                         lav_matrix_vech(Sigma.hat[[g]], diagonal=FALSE))
        } else if(!categorical && meanstructure) {
            wls.est <- c(Mu.hat[[g]],
                         lav_matrix_vec(PI.g),
                         lav_matrix_vech(Sigma.hat[[g]]))
        } else {
            wls.est <- c(lav_matrix_vec(PI.g),
                         lav_matrix_vech(Sigma.hat[[g]]))
        }
        if(group.w.free) {
            #wls.est <- c(log(GW[[g]]/GW[[samplestats@ngroups]]), wls.est)
             wls.est <- c(GW[[g]], wls.est)
        }

        WLS.est[[g]] <- wls.est
    }

    WLS.est
}

# Note: lav_model_wls_v() is replaced by lav_model_h1_information() in 0.6-1
