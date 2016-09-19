# compute model implied statistics
lav_model_implied <- function(lavmodel = NULL) {

    stopifnot(inherits(lavmodel, "Model"))
    
    # model-implied variance/covariance matrix ('sigma hat')
    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)

    # model-implied mean structure ('mu hat')
    Mu.hat <-    computeMuHat(lavmodel = lavmodel)

    # if conditional.x, slopes
    if(lavmodel@conditional.x) {
        SLOPES <- computePI(lavmodel = lavmodel)
    } else {
        SLOPES <- vector("list", length = lavmodel@ngroups)
    }

    # if categorical, model-implied thresholds
    if(lavmodel@categorical) {
        TH <- computeTH(lavmodel = lavmodel)
    } else {
        TH <- vector("list", length = lavmodel@ngroups)
    }
 
    if(lavmodel@group.w.free) {
        w.idx <- which(names(lavmodel@GLIST) == "gw")
        GW <- lavmodel@GLIST[ w.idx ]
    } else {
        GW <- vector("list", length = lavmodel@ngroups)
    }

    # FIXME: should we use 'res.cov', 'res.int', 'res.th' if conditionl.x??
    # Yes, since 0.5-22
    if(lavmodel@conditional.x) {
        implied <- list(res.cov = Sigma.hat, res.int = Mu.hat, res.slopes = SLOPES, res.th = TH, GW = GW)
    } else {
        implied <- list(cov = Sigma.hat, mean = Mu.hat, slopes = SLOPES, th = TH, GW = GW)
    }

    implied
}
