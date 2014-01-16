# casewise likelihoods

lav_model_lik_mml <- function(lavmodel    = NULL,
                              THETA       = NULL,
                              TH          = NULL,
                              GLIST       = NULL,
                              group       = 1L,
                              lavdata     = NULL,
                              sample.mean = NULL,
                              link        = "logit",
                              lavcache    = NULL) {

    # data for this group
    X <- lavdata@X[[group]]; nobs <- nrow(X); nvar <- ncol(X)
    eXo <- lavdata@eXo[[group]]

    # MLIST (for veta and yhat)
    mm.in.group <- 1:lavmodel@nmat[group] + cumsum(c(0,lavmodel@nmat))[group]
    MLIST <- GLIST[ mm.in.group ]

    # quadrature points
    GH <- lavcache[[group]]$GH; nGH <- nrow(GH$x)

    # whitening?
    nfac <- ncol(GH$x)
    if(nfac > 1L) {
        # whitening
        VETA <- computeVETA.LISREL(MLIST = MLIST)
        chol.VETA <- chol(VETA)
        GHx <- GH$x %*% chol.VETA
    } else {
        GHx <- GH$x
    }

    # compute (log)lik for each node, for each observation
    SUM.LOG.FY <- matrix(0, nrow=nGH, ncol=nobs)
    for(q in 1:nGH) {

        # first, compute yhat for this node (eta)
        yhat <- computeYHATx.LISREL(MLIST  = MLIST,
                    eXo    = eXo,
                    ETA    = GHx[q,,drop=FALSE],
                    sample.mean = sample.mean,
                    ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
                    ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
                    ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
                    ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]])

        # compute fy.var, for this node (eta): P(Y_i =  y_i | eta_i, x_i)
        log.fy.var <- lav_predict_fy_internal(X = X, yhat = yhat,
                          TH = TH, THETA = THETA,
                          num.idx = lavmodel@num.idx[[group]],
                          th.idx  = lavmodel@th.idx[[group]],
                          link = link, log. = TRUE)

        # if log, fy is just the sum of log.fy.var
        log.fy <- apply(log.fy.var, 1L, sum)

        # store log likelihoods for this node
        SUM.LOG.FY[q,] <- log.fy
        SUM.LOG.FY[q,] <- log.fy
    }

    # integration
    lik <- as.numeric( t(GH$w) %*% exp(SUM.LOG.FY) )

    lik
}
