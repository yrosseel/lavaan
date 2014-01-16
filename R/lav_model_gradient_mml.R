lav_model_gradient_mml <- function(lavmodel    = NULL,
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

    # compute case-wise likelihoods 
    lik <- lav_model_lik_mml(lavmodel = lavmodel, THETA = THETA, TH = TH,
               GLIST = GLIST, group = group, lavdata = lavdata,
               sample.mean = sample.mean, link = link, lavcache = lavcache)

    # compute dL/dx for each node
    dLdx <-  matrix(0, nGH, lavmodel@nx.free)
    for(q in 1:nGH) {

        # again, compute yhat for this node (eta)
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

        dFYp <- dFYp_x(X = X, yhat = yhat, THETA = THETA, TH = TH,
                    lavmodel = lavmodel, g = group, # for computeDeltaDx
                    ETA = GHx[q,,drop=FALSE],
                    FY = exp(log.fy.var), ## FIXME log/exp/log/...
                    num.idx = lavmodel@num.idx[[group]],
                    th.idx  = lavmodel@th.idx[[group]],
                    link = link)

        dLdx[q,] <- apply(1/lik * dFYp, 2, sum)
    }
  
    # integration
    dx <- apply(as.numeric(GH$w) * dLdx, 2, sum)

    # minimize
    dx <- -1*dx

    dx
}

dFYp_x <- function(X        = NULL,
                   yhat     = NULL,
                   THETA    = NULL, 
                   TH       = NULL,
                   lavmodel = NULL,
                   g        = NULL,
                   ETA      = NULL,
                    FY      = NULL,
                   num.idx  = NULL,
                   th.idx   = NULL,
                   link     = NULL) {

    nobs <- nrow(X); nvar <- ncol(X)

    Delta.lambda <- computeDeltaDx(lavmodel, target="lambda")[[g]]
    Delta.th     <- computeDeltaDx(lavmodel, target="th"    )[[g]]
    Delta.mu     <- computeDeltaDx(lavmodel, target="mu"    )[[g]]
    #Delta.sigma  <- computeDeltaDx(lavmodel, target="sigma" )[[g]]
    Delta.theta  <- computeDeltaDx(lavmodel, target="theta" )[[g]]

    # ordered?
    ord.idx <- unique( th.idx[th.idx > 0L] )

    dFYp <- matrix(0, nobs, lavmodel@nx.free)
    for(p in 1:nvar) {

        nfac <- nrow(Delta.lambda)/nvar
        lambda.idx <- nvar*((1:nfac) - 1L) + p
        theta.idx <- lavaan:::diag.idx(nvar)[p]
        nu.idx <- p

        # prod minus p itself
        FYp <- as.numeric(apply(FY[,-p], 1L, prod))

        # numeric
        if(p %in% num.idx) {
            y <- X[,p]
            sd.v <- sqrt(THETA[p,p])
            dy <- dnorm(y, mean=yhat[,p], sd=sd.v)
            # lambda
            dlambda <- 1/sd.v^2 * (y - yhat[,p]) * ETA * dy
            # theta
            dsigma2 <- (1/(2*sd.v^4)*(y - yhat[,p])^2 - 1/(2*sd.v^2)) * dy
            # nu
            dnu <-  1/sd.v^2 * (y - yhat[,p]) * dy

            dFYp <- dFYp +
                ( (dlambda %*% Delta.lambda[lambda.idx,,drop=FALSE]) * FYp)
            # FIXME: sigma[j,j] or THETA[j,j] ????
            #sigma.idx <- lavaan:::diagh.idx(nvar)[p]
            #dFYp <- dFYp +
            #    ( (dsigma2 %*% Delta.sigma[sigma.idx,,drop=FALSE]) * FYp)
            dFYp <- dFYp +
                ( (dsigma2 %*% Delta.theta[theta.idx,,drop=FALSE]) * FYp)
            dFYp <- dFYp +
               ( (dnu %*% Delta.mu[nu.idx,,drop=FALSE]) * FYp)
        }

        # ordinal
        if(p %in% ord.idx) {

            # just in case we need theta[v,v] after all...
            sd.v <- sqrt(THETA[p,p])
            # lav_probit
            y <- X[,p]
            th.y <- TH[ th.idx == p]; TH.Y <- c(-Inf, th.y, Inf)
            ncat <- length(th.y) + 1L; nth <- ncat - 1L
            Y1 <- matrix(1:nth, nobs, nth, byrow=TRUE) == y
            Y2 <- matrix(1:nth, nobs, nth, byrow=TRUE) == (y - 1L)
            z1 <- pmin( 100, TH.Y[y+1L   ] - yhat[,p])
            z2 <- pmax(-100, TH.Y[y+1L-1L] - yhat[,p])
            p1 <- dnorm(z1)
            p2 <- dnorm(z2)
            # d tau
            dth <- -1 * (Y2*p2 - Y1*p1)/sd.v
            # d lambda
            dlambda <- -1*(p1-p2)*ETA

            dFYp <- dFYp +
                ( (dth %*% Delta.th[th.idx==p,,drop=FALSE]) * FYp)
            dFYp <- dFYp +
                ( (dlambda %*% Delta.lambda[lambda.idx,,drop=FALSE]) * FYp)
        }
    }

    dFYp
}

