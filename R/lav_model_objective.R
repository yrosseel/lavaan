# model objective

lav_model_objective <- function(lavmodel       = NULL,
                                GLIST          = NULL,
                                lavsamplestats = NULL,
                                lavdata        = NULL,
                                lavcache       = NULL,
                                estimator      = "ML",
                                verbose        = FALSE,
                                forcePD        = TRUE,
                                debug          = FALSE) {

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    # shortcut for data.type == "none" or estimator == "none"
    if(estimator == "none" || length(lavsamplestats@cov) == 0L) {
        fx <- as.numeric(NA)
        attr(fx, "fx.group") <- rep(as.numeric(NA), lavsamplestats@ngroups)
        return(fx)
    }

    meanstructure <- lavmodel@meanstructure
    categorical   <- lavmodel@categorical
    group.w.free  <- lavmodel@group.w.free
    fixed.x       <- lavmodel@fixed.x
    num.idx       <- lavmodel@num.idx
    th.idx        <- lavmodel@th.idx

    # do we need WLS.est?
    if(estimator == "GLS"  || estimator == "WLS"  || 
       estimator == "DWLS" || estimator == "ULS") {
        WLS.est <- lav_model_wls_est(lavmodel = lavmodel, GLIST = GLIST)
        if(debug) print(WLS.est)
    } else if(estimator == "ML" || estimator == "PML" || 
              estimator == "FML" || estimator == "REML") {
        # compute moments for all groups
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, GLIST = GLIST, 
                                     extra = (estimator %in% c("ML", "REML")))

        if(estimator == "REML") {
            LAMBDA <- computeLAMBDA(lavmodel = lavmodel, GLIST = GLIST)
        }

        # ridge?
        if( lavsamplestats@ridge > 0.0 ) {
            for(g in 1:lavsamplestats@ngroups) {
                diag(Sigma.hat[[g]]) <- diag(Sigma.hat[[g]]) + 
                                            lavsamplestats@ridge
            }
        }
        if(debug) print(Sigma.hat)

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
    } else if(estimator == "MML") {
        TH    <- computeTH(   lavmodel = lavmodel, GLIST = GLIST)
        THETA <- computeTHETA(lavmodel = lavmodel, GLIST = GLIST)
        GW    <- computeGW(   lavmodel = lavmodel, GLIST = GLIST)
    }
 
    fx <- 0.0
    fx.group <- numeric( lavsamplestats@ngroups )
    logl.group <- rep(as.numeric(NA), lavsamplestats@ngroups)
    for(g in 1:lavsamplestats@ngroups) {

        # incomplete data and fiml?
        if(lavsamplestats@missing.flag) {
            if(estimator == "ML") {
                # FIML
                if(!attr(Sigma.hat[[g]], "po")) return(Inf)
                group.fx <- estimator.FIML(Sigma.hat=Sigma.hat[[g]],
                                           Mu.hat=Mu.hat[[g]],
                                           M=lavsamplestats@missing[[g]],
                                           h1=lavsamplestats@missing.h1[[g]]$h1)
            } else {
                stop("this estimator: `", estimator, 
                     "' can not be used with incomplete data and the missing=\"ml\" option")
            }
        } else if(estimator == "ML") {
        # complete data
            # ML and friends
            group.fx <- estimator.ML(Sigma.hat=Sigma.hat[[g]], 
                                     Mu.hat=Mu.hat[[g]],
                                     data.cov=lavsamplestats@cov[[g]], 
                                     data.mean=lavsamplestats@mean[[g]], 
                                     data.cov.log.det=lavsamplestats@cov.log.det[[g]],
                                     meanstructure=meanstructure)
        } else if(estimator == "GLS"  || estimator == "WLS") {
            # full weight matrix
            group.fx <- estimator.WLS(WLS.est = WLS.est[[g]],
                                      WLS.obs = lavsamplestats@WLS.obs[[g]],
                                      WLS.V=lavsamplestats@WLS.V[[g]])
            attr(group.fx, "WLS.est") <- WLS.est[[g]]

        } else if(estimator == "DWLS" || estimator == "ULS") {
            # diagonal weight matrix
            group.fx <- estimator.DWLS(WLS.est = WLS.est[[g]],
                                      WLS.obs = lavsamplestats@WLS.obs[[g]], 
                                      WLS.VD = lavsamplestats@WLS.VD[[g]])
            attr(group.fx, "WLS.est") <- WLS.est[[g]]

        } else if(estimator == "PML") {
            # Pairwise maximum likelihood
            group.fx <- estimator.PML(Sigma.hat = Sigma.hat[[g]],
                                      TH        = TH[[g]],
                                      th.idx    = th.idx[[g]],
                                      num.idx   = num.idx[[g]],
                                      X         = lavdata@X[[g]],
                                      lavcache  = lavcache[[g]])
            logl.group[g] <- attr(group.fx, "logl")

        } else if(estimator == "FML") { 
            # Full maximum likelihood (underlying multivariate normal)
            group.fx <- estimator.FML(Sigma.hat = Sigma.hat[[g]],
                                      TH        = TH[[g]],
                                      th.idx    = th.idx[[g]],
                                      num.idx   = num.idx[[g]],
                                      X         = lavdata@X[[g]],
                                      lavcache  = lavcache[[g]])

        } else if(estimator == "MML") { 
            # marginal maximum likelihood
            group.fx <- estimator.MML(lavmodel    = lavmodel,
                                      GLIST       = GLIST,
                                      THETA       = THETA[[g]],
                                      TH          = TH[[g]],
                                      group       = g,
                                      lavdata     = lavdata,
                                      sample.mean = lavsamplestats@mean[[g]],
                                      lavcache    = lavcache)
        } else if(estimator == "REML") {
            # restricted/residual maximum likelihood
            group.fx <- estimator.REML(Sigma.hat = Sigma.hat[[g]],
                                       Mu.hat=Mu.hat[[g]],
                                       data.cov=lavsamplestats@cov[[g]],
                                       data.mean=lavsamplestats@mean[[g]],
                                       data.cov.log.det=lavsamplestats@cov.log.det[[g]],
                                       meanstructure=meanstructure,
                                       group = g,
                                       lavmodel = lavmodel,
                                       lavsamplestats = lavsamplestats,
                                       lavdata = lavdata)
        } else {
            stop("unsupported estimator: ", estimator)
        }

        if(estimator == "ML" || estimator == "REML") {
            group.fx <- 0.5 * group.fx ## FIXME
        } else if(estimator == "PML" || estimator == "FML" || 
                  estimator == "MML") {
            # do nothing
        } else {
            group.fx <- 0.5 * (lavsamplestats@nobs[[g]]-1)/lavsamplestats@nobs[[g]] * group.fx
        }

        fx.group[g] <- group.fx
    } # g

    if(lavsamplestats@ngroups > 1) {
        ## FIXME: if group.w.free, should we use group.w or nobs???
        ##  - if we use estimated group.w, gradient changes!!!!
        ##  - but, if group models are misspecified, the group weights
        ##    will be affected too... which is unwanted (I think)
        #if(group.w.free) {
           # nobs <- unlist(GW) * lavsamplestats@ntotal
        #   nobs <- exp(unlist(GW))
        #} else {
        if(estimator == "PML") {
            # no weighting needed! (since N_g is part of the logl per group)
            fx <- sum(fx.group)
        } else {
            nobs <- unlist(lavsamplestats@nobs)
            #}
            fx <- weighted.mean(fx.group, w=nobs)
        }
    } else { # single group
        fx <- fx.group[1]
    }

    # penalty for group.w + ML
    if(group.w.free && estimator %in% c("ML","MML","FML","PML", "REML")) {
        #obs.prop <- unlist(lavsamplestats@group.w)
        #est.prop <- unlist(GW)
        # if(estimator %in% c("WLS", "GLS", ...) {
        #    # X2 style discrepancy measures (aka GLS/WLS!!)
        #    fx.w <- sum ( (obs.prop-est.prop)^2/est.prop )
        # } else {
        #    # G2 style discrepancy measures (aka ML)
        #    # deriv is here -2 * (obs.prop - est.prop)
        #fx.w <- sum(obs.prop * log(obs.prop/est.prop) )
        # }
        
        # poisson kernel
        obs.freq <- unlist(lavsamplestats@group.w) * lavsamplestats@ntotal
        est.freq <- exp(unlist(GW))
        fx.w <- -1 * sum( obs.freq * log(est.freq) - est.freq )
        # divide by N (to be consistent with the rest of lavaan)
        fx.w <- fx.w / lavsamplestats@ntotal

        fx.sat <- sum( obs.freq * log(obs.freq) - obs.freq )
        fx.sat <- fx.sat / lavsamplestats@ntotal

        # saturated - poisson
        #fx.w <- sum(obs.freq * log(obs.freq/est.freq))
        # does not work without constraints?

        fx <- fx + (fx.w + fx.sat)
    }

    fx.value <- as.numeric(fx)

    attr(fx, "fx.group") <- fx.group
    if(estimator == "PML") {
        attr(fx, "logl.group") <- logl.group
        attr(fx, "fx.pml") <- fx.value
    }

    fx
}

