# casewise likelihoods

lav_model_lik_mml <- function(lavmodel    = NULL,
                              THETA       = NULL,
                              TH          = NULL,
                              GLIST       = NULL,
                              group       = 1L,
                              lavdata     = NULL,
                              sample.mean = NULL,
                              control     = list(),
                              lavcache    = NULL) {

    # data for this group
    X <- lavdata@X[[group]]; nobs <- nrow(X); nvar <- ncol(X)
    eXo <- lavdata@eXo[[group]]

    # MLIST (for veta and yhat)
    mm.in.group <- 1:lavmodel@nmat[group] + cumsum(c(0,lavmodel@nmat))[group]
    MLIST <- GLIST[ mm.in.group ]

    # quadrature points
    GH <- lavcache[[group]]$GH; nGH <- nrow(GH$x)
    nfac <- ncol(GH$x)

    # compute VETAx (latent lv only)
    lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[group]],
                      lavmodel@ov.x.dummy.lv.idx[[group]])
    VETAx <- computeVETAx.LISREL(MLIST = MLIST,
                                 lv.dummy.idx = lv.dummy.idx)
    # check for negative values?
    if(any(diag(VETAx) < 0)) {
        warning("lavaan WARNING: --- VETAx contains negative values")
        print(VETAx)
        return(0)
    }

    # cholesky?
    if(is.null(control$cholesky)) {
        CHOLESKY <- TRUE    
    } else {
        CHOLESKY <- control$cholesky
        if(nfac > 1L && !CHOLESKY) {
            warning("lavaan WARNING: CHOLESKY is OFF but nfac > 1L")
        }
    }

    if(!CHOLESKY) {
        # we should still 'scale' the factors, if std.lv=FALSE
        ETA.sd <- sqrt( diag(VETAx) )
    } else {
        # cholesky takes care of scaling
        ETA.sd <- rep(1, nfac)
        chol.VETA <- try(chol(VETAx), silent = TRUE)
        if(inherits(chol.VETA, "try-error")) {
            warning("lavaan WARNING: --- VETAx not positive definite")
            print(VETAx)
            return(0)
        }
    }

    # compute (log)lik for each node, for each observation
    SUM.LOG.FY <- matrix(0, nrow=nGH, ncol=nobs)
    for(q in 1:nGH) {

        # current value(s) for ETA
        eta <- GH$x[q,,drop=FALSE]

        # rescale/unwhiten
        if(CHOLESKY) {
            eta <- eta %*% chol.VETA
        } else {
            # no unit scale? (un-standardize)
            eta <- sweep(eta, MARGIN=2, STATS=ETA.sd, FUN="*")
        }

        # eta_i = alpha + BETA eta_i + GAMMA eta_i + error
        eta <- computeETAx.LISREL(MLIST = MLIST, eXo = eXo, ETA = eta,
                   remove.dummy.lv = TRUE,
                   ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
                   ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]],
                   Nobs = 1L)

        # compute yhat for this node (eta)
        yhat <- computeYHATetax.LISREL(MLIST = MLIST, eXo = eXo,
                    ETA = eta, sample.mean = sample.mean,
                    ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
                    ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
                    ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
                    ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]])

        # compute fy.var, for this node (eta): P(Y_i =  y_i | eta_i, x_i)
        log.fy.var <- lav_predict_fy_internal(X = X, yhat = yhat,
                          TH = TH, THETA = THETA,
                          num.idx = lavmodel@num.idx[[group]],
                          th.idx  = lavmodel@th.idx[[group]],
                          link    = lavmodel@link, log. = TRUE)

        # if log, fy is just the sum of log.fy.var
        log.fy <- apply(log.fy.var, 1L, sum)

        # store log likelihoods for this node
        SUM.LOG.FY[q,] <- log.fy
    }

    # integration
    lik <- as.numeric( t(GH$w) %*% exp(SUM.LOG.FY) )

    # avoid underflow
    idx <- which(lik < exp(-600))
    if(length(idx) > 0L) {
        lik[idx] <- exp(-600)
    }

    lik
}
