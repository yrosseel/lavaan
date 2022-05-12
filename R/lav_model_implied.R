# compute model implied statistics
# per block

# YR 7 May 2022: add cov.x and mean.x if conditional.x (so that we do
#                no longer depend on SampleStats)

lav_model_implied <- function(lavmodel = NULL, GLIST = NULL, delta = TRUE) {

    stopifnot(inherits(lavmodel, "lavModel"))

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    # model-implied variance/covariance matrix ('sigma hat')
    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, GLIST = GLIST,
                                 delta = delta)

    # model-implied mean structure ('mu hat')
    if(lavmodel@meanstructure) {
        Mu.hat <- computeMuHat(lavmodel = lavmodel,  GLIST = GLIST)
    } else {
        Mu.hat <- vector("list", length = lavmodel@nblocks)
    }

    # if conditional.x, slopes, cov.x, mean.x
    if(lavmodel@conditional.x) {
        SLOPES <- computePI(lavmodel = lavmodel,  GLIST = GLIST)

        # per block, because for some blocks, cov.x may not exist
        COV.X  <- vector("list", lavmodel@nblocks)
        MEAN.X <- vector("list", lavmodel@nblocks)
        for(b in seq_len(lavmodel@nblocks)) {
            mm.in.block <- ( seq_len(lavmodel@nmat[b]) +
                             cumsum(c(0, lavmodel@nmat))[b] )
            MLIST <- lavmodel@GLIST[ mm.in.block ]
            cov.x.idx <- which( names(MLIST) == "cov.x" )
            if(length(cov.x.idx) > 0L) {
                COV.X[[b]] <- MLIST[[cov.x.idx]]
            } else {
                COV.X[[b]] <- matrix(0, 0L, 0L)
            }
            mean.x.idx <- which( names(MLIST) == "mean.x" )
            if(length(mean.x.idx) > 0L) {
                MEAN.X[[b]] <- MLIST[[mean.x.idx]]
            } else {
                MEAN.X[[b]] <- matrix(0, 0L, 0L)
            }
        }
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

    if(lavmodel@conditional.x) {
        implied <- list(res.cov = Sigma.hat, res.int = Mu.hat,
                        res.slopes = SLOPES, cov.x = COV.X, mean.x = MEAN.X,
                        res.th = TH, group.w = GW)
    } else {
        implied <- list(cov = Sigma.hat, mean = Mu.hat, th = TH, group.w = GW)
    }

    implied
}
