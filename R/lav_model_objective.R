# model objective

lav_model_objective <- function(object, GLIST=NULL, 
                                samplestats=NULL, X = NULL,
                                cache=NULL,
                                estimator="ML", link="logit",
                                verbose=FALSE, forcePD=TRUE,
                                debug=FALSE) {

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    # shortcut for data.type == "none" or estimator == "none"
    if(estimator == "none" || length(samplestats@cov) == 0L) {
        fx <- as.numeric(NA)
        attr(fx, "fx.group") <- rep(as.numeric(NA), samplestats@ngroups)
        return(fx)
    }

    meanstructure <- object@meanstructure
    categorical   <- object@categorical
    group.w.free  <- object@group.w.free
    fixed.x       <- object@fixed.x
    num.idx       <- object@num.idx
    th.idx        <- object@th.idx

    # do we need WLS.est?
    if(estimator == "GLS"  || estimator == "WLS"  || 
       estimator == "DWLS" || estimator == "ULS") {
        WLS.est <- lav_model_wls_est(object = object, GLIST = GLIST)
        if(debug) print(WLS.est)
    } else {
        # compute moments for all groups
        Sigma.hat <- computeSigmaHat(object, GLIST=GLIST, 
                                     extra=(estimator=="ML"))
        if(debug) print(Sigma.hat)
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
    }
 
    fx <- 0.0
    fx.group <- numeric( samplestats@ngroups )
    for(g in 1:samplestats@ngroups) {

        # ridge?
        if( samplestats@ridge > 0.0 ) {
            diag(Sigma.hat[[g]]) <- diag(Sigma.hat[[g]]) + samplestats@ridge
        }

        # incomplete data and fiml?
        if(samplestats@missing.flag) {
            if(estimator == "ML") {
                # FIML
                if(!attr(Sigma.hat[[g]], "po")) return(Inf)
                group.fx <- estimator.FIML(Sigma.hat=Sigma.hat[[g]],
                                           Mu.hat=Mu.hat[[g]],
                                           M=samplestats@missing[[g]],
                                           h1=samplestats@missing.h1[[g]]$h1)
            } else {
                stop("this estimator: `", estimator, 
                     "' can not be used with incomplete data and the missing=\"ml\" option")
            }
        } else if(estimator == "ML") {
        # complete data
            # ML and friends
            group.fx <- estimator.ML(Sigma.hat=Sigma.hat[[g]], 
                                     Mu.hat=Mu.hat[[g]],
                                     data.cov=samplestats@cov[[g]], 
                                     data.mean=samplestats@mean[[g]], 
                                     data.cov.log.det=samplestats@cov.log.det[[g]],
                                     meanstructure=meanstructure)
        } else if(estimator == "GLS"  || estimator == "WLS"  || 
                  estimator == "DWLS" || estimator == "ULS") {
            group.fx <- estimator.WLS(WLS.est = WLS.est[[g]],
                                      WLS.obs = samplestats@WLS.obs[[g]], 
                                      WLS.V=samplestats@WLS.V[[g]])  
            attr(group.fx, "WLS.est") <- WLS.est[[g]]
        } else if(estimator == "PML") {
            # Pairwise maximum likelihood
            group.fx <- estimator.PML(Sigma.hat = Sigma.hat[[g]],
                                      TH        = TH[[g]],
                                      th.idx    = th.idx[[g]],
                                      num.idx   = num.idx[[g]],
                                      X         = X[[g]],
                                      cache     = cache[[g]])
        } else if(estimator == "FML") { 
            # Full maximum likelihood (underlying multivariate normal)
            group.fx <- estimator.FML(Sigma.hat = Sigma.hat[[g]],
                                      TH        = TH[[g]],
                                      th.idx    = th.idx[[g]],
                                      num.idx   = num.idx[[g]],
                                      X         = X[[g]],
                                      cache     = cache[[g]])
        } else if(estimator == "MML") { 
            # marginal maximum likelihood
            group.fx <- estimator.MML(object    = object,
                                      g         = g,
                                      GLIST     = GLIST,
                                      sample.mean = samplestats@mean[[g]],
                                      link      = link,
                                      X         = X[[g]],
                                      cache     = cache[[g]])
        } else {
            stop("unsupported estimator: ", estimator)
        }

        if(estimator == "ML") {
            group.fx <- 0.5 * group.fx ## FIXME
        } else if(estimator == "PML" || estimator == "FML" || 
                  estimator == "MML") {
            # do nothing
        } else {
            group.fx <- 0.5 * (samplestats@nobs[[g]]-1)/samplestats@nobs[[g]] * group.fx
        }

        fx.group[g] <- group.fx
    } # g

    if(samplestats@ngroups > 1) {
        ## FIXME: if group.w.free, should we use group.w or nobs???
        ##  - if we use estimated group.w, gradient changes!!!!
        ##  - but, if group models are misspecified, the group weights
        ##    will be affected too... which is unwanted (I think)
        #if(group.w.free) {
           # nobs <- unlist(GW) * samplestats@ntotal
        #   nobs <- exp(unlist(GW))
        #} else {
           nobs <- unlist(samplestats@nobs)
        #}
        fx <- weighted.mean(fx.group, w=nobs)
    } else { # single group
        fx <- fx.group[1]
    }

    # penalty for group.w + ML
    if(group.w.free && estimator %in% c("ML","MML","FML","PML")) {
        #obs.prop <- unlist(samplestats@group.w)
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
        obs.freq <- unlist(samplestats@group.w) * samplestats@ntotal
        est.freq <- exp(unlist(GW))
        fx.w <- -1 * sum( obs.freq * log(est.freq) - est.freq )
        # divide by N (to be consistent with the rest of lavaan)
        fx.w <- fx.w / samplestats@ntotal

        fx.sat <- sum( obs.freq * log(obs.freq) - obs.freq )
        fx.sat <- fx.sat / samplestats@ntotal

        # saturated - poisson
        #fx.w <- sum(obs.freq * log(obs.freq/est.freq))
        # does not work without constraints?

        fx <- fx + (fx.w + fx.sat)
    }

    attr(fx, "fx.group") <- fx.group

    fx
}

