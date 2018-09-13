# compute model implied statistics
# per block
lav_model_implied <- function(lavmodel = NULL, GLIST = NULL) {

    stopifnot(inherits(lavmodel, "lavModel"))

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    # model-implied variance/covariance matrix ('sigma hat')
    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, GLIST = GLIST)

    # model-implied mean structure ('mu hat')
    if(lavmodel@meanstructure) {
        Mu.hat <- computeMuHat(lavmodel = lavmodel,  GLIST = GLIST)
    } else {
        Mu.hat <- vector("list", length = lavmodel@nblocks)
    }

    # if conditional.x, slopes
    if(lavmodel@conditional.x) {
        SLOPES <- computePI(lavmodel = lavmodel,  GLIST = GLIST)
    } else {
        SLOPES <- vector("list", length = lavmodel@nblocks)
    }

    # if categorical, model-implied thresholds
    if(lavmodel@categorical) {
        TH <- computeTH(lavmodel = lavmodel,  GLIST = GLIST)
    } else {
        TH <- vector("list", length = lavmodel@nblocks)
    }

    if(lavmodel@group.w.free) {
        w.idx <- which(names(lavmodel@GLIST) == "gw")
        GW <- unname(GLIST[ w.idx ])
        GW <- lapply(GW, as.numeric)
    } else {
        GW <- vector("list", length = lavmodel@nblocks)
    }

    # FIXME: should we use 'res.cov', 'res.int', 'res.th' if conditionl.x??
    # Yes, since 0.5-22
    if(lavmodel@conditional.x) {
        implied <- list(res.cov = Sigma.hat, res.int = Mu.hat, res.slopes = SLOPES, res.th = TH, group.w = GW)
    } else {
        implied <- list(cov = Sigma.hat, mean = Mu.hat, th = TH, group.w = GW)
    }

    implied
}
