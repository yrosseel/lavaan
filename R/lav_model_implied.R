# compute model implied statistics
lav_model_implied <- function(lavmodel = NULL) {

    stopifnot(inherits(lavmodel, "Model"))
    
    # model-implied variance/covariance matrix ('sigma hat')
    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)

    # model-implied mean structure ('mu hat')
    Mu.hat <-    computeMuHat(lavmodel = lavmodel)

    # if categorical, model-implied thresholds
    if(lavmodel@categorical) {
        TH <- computeTH(lavmodel = lavmodel)
    } else {
        TH <- vector("list", length = lavmodel@ngroups)
    }

    list(cov = Sigma.hat, mean = Mu.hat, th = TH)
}
