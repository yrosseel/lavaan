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

# compute WLS.V (as a list per group)
#
# three options:
#  1) *LS: WLS.V is already in lavsamplestats
#  2) NTRLS: WLS.V needs to recomputed after every iteration, using
#            the structured estimates of Sigma/Mu
#  3) ML: 3a: complete, structured (default)
#         3b: complete, unstructured
#         3c: incomplete, FIML, structured
#         3d: incomplete, FIML, unstructured
#         3e: incomplete, two.stage, structured
#         ef: incomplete, two.stage, unstructured (EM estimates)
lav_model_wls_v <- function(lavmodel       = NULL, 
                            lavsamplestats = NULL,
                            structured     = TRUE,
                            lavimplied     = NULL,
                            lavdata        = NULL) {

    WLS.V <- vector("list", length=lavsamplestats@ngroups)

    # 1) *LS: WLS.V is already in lavsamplestats
    if(lavmodel@estimator == "GLS"  || lavmodel@estimator == "WLS") {
        # for GLS, the WLS.V22 part is: 0.5 * t(D) %*% [S.inv %x% S.inv] %*% D
        # for WLS, the WLS.V22 part is: Gamma
        WLS.V <- lavsamplestats@WLS.V
    } else if(lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
        # diagonal only!!
        WLS.V <- lavsamplestats@WLS.VD


    # 2) NTRLS: based on structured estimates of Sigma/Mu
    } else if(lavmodel@estimator == "NTRLS") {
        stopifnot(!lavmodel@conditional.x)

        # by definition, we always use the 'structured' moments
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
        Mu.hat <- computeMuHat(lavmodel = lavmodel)

        for(g in 1:lavsamplestats@ngroups) {
            WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(
                              COV            = Sigma.hat[[g]][,,drop=FALSE],
                              MEAN           = Mu.hat[[g]],
                              x.idx          = lavsamplestats@x.idx[[g]],
                              fixed.x        = lavmodel@fixed.x,
                              conditional.x  = lavmodel@conditional.x,
                              meanstructure  = lavmodel@meanstructure,
                              slopestructure = lavmodel@conditional.x)
        }

    #  3) ML: 3a: complete, structured (default)
    #         3b: complete, unstructured
    #         3c: incomplete, FIML, structured
    #         3d: incomplete, FIML, unstructured
    #         3e: incomplete, two.stage, structured
    #         ef: incomplete, two.stage, unstructured (EM estimates)
    } else if(lavmodel@estimator == "ML") {
        WLS.V <- vector("list", length=lavsamplestats@ngroups)
        
        if(structured) {
            if(lavmodel@conditional.x) {
                 Sigma.hat <- computeSigmaHatJoint(lavmodel = lavmodel,
                                 lavsamplestats = lavsamplestats)
            } else {
                 Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
            }
            if(lavmodel@meanstructure) {
                if(lavmodel@conditional.x) {
                    Mu.hat <- computeMuHatJoint(lavmodel = lavmodel, 
                                  lavsamplestats = lavsamplestats)
                } else {
                    Mu.hat <- computeMuHat(lavmodel = lavmodel)
                }
            } else {
                Mu.hat <- NULL
            }
        } else {
            if(lavmodel@conditional.x) {
                 # FIXME: wahat to do here?
                 stop("lavaan ERROR: conditional.x = TRUE, but structured = FALSE?")
            } else {
                 # complete data: observed var/cov matrix 'S'
                 # two.stage + incomplete data: EM estimate of 'S'
                 Sigma.hat <- lavsamplestats@cov
            }
            if(lavmodel@meanstructure) {
                if(lavmodel@conditional.x) {
                    # FIXME!
                } else {
                    # complete data: observed mean vector 'ybar'
                    # two.stage + incomplete data: EM estimate of 'ybar'
                    Mu.hat <- lavsamplestats@mean
                }
            } else {
                Mu.hat <- NULL
            }
            
        }

        # GW?
        if(lavmodel@group.w.free) {
            GW <- unlist(computeGW(lavmodel = lavmodel))
        }


        # three options
        #  - complete data
        #  - incomplete data + FIML
        #  - incomplete data + two.stage
        # two variants:
        #  - using unstructured moments
        #  - using structured moments
        for(g in 1:lavsamplestats@ngroups) {
            if(lavsamplestats@missing.flag) {
                stopifnot(!lavmodel@conditional.x)
                WLS.V[[g]] <- lav_mvnorm_missing_information_expected(
                                  Y     = lavdata@X[[g]],
                                  Mp    = lavdata@Mp[[g]],
                                  Mu    = Mu.hat[[g]],
                                  Sigma = Sigma.hat[[g]])
            } else {
                WLS.V[[g]] <- lav_samplestats_Gamma_inverse_NT(
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

