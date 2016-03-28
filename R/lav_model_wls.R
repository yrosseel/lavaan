# compute WLS.est (as a list per group)
lav_model_wls_est <- function(lavmodel = NULL, GLIST = NULL) { #,
                              #cov.x = NULL) {

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
    }
    if(lavmodel@conditional.x && lavmodel@nexo > 0L) {
        PI <- computePI(lavmodel = lavmodel, GLIST = GLIST)
    } 
    if(group.w.free) {
        GW <- computeGW(lavmodel = lavmodel, GLIST = GLIST)
    }

    WLS.est <- vector("list", length=ngroups)
    for(g in 1:ngroups) {

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
    } else if(estimator == "NTRLS") {
        stopifnot(!lavmodel@conditional.x)
        # compute moments for all groups
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, extra = TRUE)
        Mu.hat <- computeMuHat(lavmodel = lavmodel)

        for(g in 1:lavsamplestats@ngroups) {
            WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(
                              ICOV = attr(Sigma.hat[[g]],"inv")[,,drop=FALSE],
                              COV            = Sigma.hat[[g]][,,drop=FALSE],
                              MEAN           = Mu.hat[[g]],
                              x.idx          = lavsamplestats@x.idx[[g]],
                              fixed.x        = lavmodel@fixed.x,
                              conditional.x  = lavmodel@conditional.x,
                              meanstructure  = lavmodel@meanstructure,
                              slopestructure = lavmodel@conditional.x)
        }
    } else if(estimator == "DWLS" || estimator == "ULS") {
        # diagonal only!!
        WLS.V <- lavsamplestats@WLS.VD


    # for ML, we need to recompute this, as it is function of Sigma (and Mu)
    } else if(estimator == "ML") {
        WLS.V <- vector("list", length=lavsamplestats@ngroups)
        
        if(lavmodel@conditional.x) {
             Sigma.hat <- computeSigmaHatJoint(lavmodel = lavmodel,
                              lavsamplestats = lavsamplestats, extra = TRUE)
        } else {
             Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, extra = TRUE)
        }
        if(lavmodel@group.w.free) {
            GW <- unlist(computeGW(lavmodel = lavmodel))
        }
        if(lavsamplestats@missing.flag || lavmodel@conditional.x) {
            if(lavmodel@conditional.x) {
                Mu.hat <- computeMuHatJoint(lavmodel = lavmodel, 
                              lavsamplestats = lavsamplestats)
            } else {
                Mu.hat <- computeMuHat(lavmodel = lavmodel)
            }
        } else {
            Mu.hat <- NULL
        }
        for(g in 1:lavsamplestats@ngroups) {
            if(lavsamplestats@missing.flag) {
                stopifnot(!lavmodel@conditional.x)
                WLS.V[[g]] <- lav_mvnorm_missing_information_expected(
                                  Y     = lavdata@X[[g]],
                                  Mp    = lavdata@Mp[[g]],
                                  Mu    = Mu.hat[[g]],
                                  Sigma = Sigma.hat[[g]])
            } else {
                # WLS.V22 = 0.5*t(D) %*% [Sigma.hat.inv %x% Sigma.hat.inv]%*% D

                # NOTE: when fixed.x=TRUE, this will give slightly different
                # results for the SEs compared to <= 0.5-20 (and Mplus), 
                # because we set the rows/columns of exo variables to zero
                #
                # but this should be ok (although SEs for ~1 are somewhat
                # larger)

                WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(
                    ICOV           = attr(Sigma.hat[[g]],"inv")[,,drop=FALSE],
                    COV            = Sigma.hat[[g]][,,drop=FALSE],
                    MEAN           = Mu.hat[[g]],
                    x.idx          = lavsamplestats@x.idx[[g]],
                    fixed.x        = lavmodel@fixed.x,
                    conditional.x  = lavmodel@conditional.x,
                    meanstructure  = lavmodel@meanstructure,
                    slopestructure = lavmodel@conditional.x)
            }

            if(lavmodel@group.w.free) {
                # unweight!!
                a <- exp(GW[g]) / lavsamplestats@nobs[[g]]
                # a <- exp(GW[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
                WLS.V[[g]] <- lav_matrix_bdiag( matrix(a,1,1), WLS.V[[g]])
            }
        }
    } # ML

    WLS.V
}

