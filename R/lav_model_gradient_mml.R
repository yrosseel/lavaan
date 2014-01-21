lav_model_gradient_mml <- function(lavmodel    = NULL,
                                   THETA       = NULL,
                                   TH          = NULL,
                                   GLIST       = NULL,
                                   group       = 1L,
                                   lavdata     = NULL,
                                   sample.mean = NULL,
                                   lavcache    = NULL) {

    if(lavmodel@link == "logit") 
        stop("logit link not implemented yet; use probit")

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
        #GHx <- GH$x %*% chol.VETA
    } else {
        chol.VETA <- NULL
    }

    # compute case-wise likelihoods 
    lik <- lav_model_lik_mml(lavmodel = lavmodel, THETA = THETA, TH = TH,
               GLIST = GLIST, group = group, lavdata = lavdata,
               sample.mean = sample.mean, lavcache = lavcache)

    # chol.VETA? we need EETAx (only if eXo!!)
    if(!is.null(chol.VETA) && !is.null(eXo)) {
        EETAx <- computeEETAx.LISREL(MLIST = MLIST, eXo = eXo,
            sample.mean = sample.mean,
            ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
            ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
            ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
            ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]])
    } else {
        EETAx <- NULL
    }

    # compute dL/dx for each node
    dLdx <-  matrix(0, nGH, lavmodel@nx.free)
    for(q in 1:nGH) {

        # again, compute yhat for this node (eta)
        yhat <- computeYHATx.LISREL(MLIST  = MLIST,
                    eXo               = eXo,
                    ETA               = GH$x[q,,drop=FALSE],
                    chol.VETA         = chol.VETA,
                    EETAx             = EETAx,
                    sample.mean       = sample.mean,
                    ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
                    ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
                    ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
                    ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]])

        # compute fy.var, for this node (eta): P(Y_i =  y_i | eta_i, x_i)
        log.fy.var <- lav_predict_fy_internal(X = X, yhat = yhat,
                          TH = TH, THETA = THETA,
                          num.idx = lavmodel@num.idx[[group]],
                          th.idx  = lavmodel@th.idx[[group]],
                          link = lavmodel@link, log. = TRUE)

        dFYp <- dFYp_x(X = X, yhat = yhat, MLIST = MLIST,
                    THETA = THETA, TH = TH,
                    lavmodel = lavmodel, g = group, # for computeDeltaDx
                    KSI = GH$x[q,,drop=FALSE], chol.VETA = chol.VETA,
                    EETAx = EETAx,
                    FY = exp(log.fy.var), ## FIXME log/exp/log/...
                    num.idx = lavmodel@num.idx[[group]],
                    th.idx  = lavmodel@th.idx[[group]])

        dLdx[q,] <- apply(1/lik * dFYp, 2, sum)
    }
  
    # integration
    dx <- apply(as.numeric(GH$w) * dLdx, 2, sum)

    # minimize
    dx <- -1*dx

    dx
}

dFYp_x <- function(X         = NULL,
                   yhat      = NULL,
                   MLIST     = NULL,
                   THETA     = NULL, 
                   TH        = NULL,
                   lavmodel  = NULL,
                   g         = NULL,
                   KSI       = NULL,
                   chol.VETA = NULL,
                   EETAx     = NULL,
                   FY        = NULL,
                   num.idx   = NULL,
                   th.idx    = NULL) {

    nobs <- nrow(X); nvar <- ncol(X)

    # chol.VETA?
    if(is.null(chol.VETA)) {
        ETA <- KSI
    } else {
        ETA <- (KSI %*% chol.VETA)
        if(!is.null(EETAx)) {
            #ETA <- ETA + EETAx  ### FIXME? dim!!
            ETA <- sweep(EETAx, MARGIN=1, STATS=ETA, FUN="+")
        }
    }

    # fix Lambda?
    LAMBDA <- computeLAMBDA.LISREL(MLIST = MLIST, 
                    ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                    ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                    ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                    ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]])

    Delta.lambda <- computeDeltaDx(lavmodel, target="lambda")[[g]]
    Delta.th     <- computeDeltaDx(lavmodel, target="th"    )[[g]]
    Delta.mu     <- computeDeltaDx(lavmodel, target="mu"    )[[g]]
    #Delta.sigma  <- computeDeltaDx(lavmodel, target="sigma" )[[g]]
    Delta.theta  <- computeDeltaDx(lavmodel, target="theta" )[[g]]

    x <- lav_model_get_parameters(lavmodel = lavmodel, GLIST = MLIST)
    dVetadx <- function(x, lavmodel = fit@Model, g = 1L) {
        GLIST <- lavaan:::lav_model_x2GLIST(lavmodel, x=x, type="free")
        VETA <- lavaan:::computeVETA(lavmodel, GLIST = GLIST)[[g]]
        S <- chol(VETA)
        S
    }
    Delta.S <- lavaan:::lavJacobianD(func=dVetadx, x=x, lavmodel = lavmodel, 
                                       g= g)

    

    # ordered?
    ord.idx <- unique( th.idx[th.idx > 0L] )

    dFYp <- matrix(0, nobs, lavmodel@nx.free)
    for(p in 1:nvar) {

        nfac <- nrow(Delta.lambda)/nvar
        lambda.idx <- nvar*((1:nfac) - 1L) + p
        theta.idx <- diagh.idx(nvar)[p]
        nu.idx <- p

        # prod minus p itself
        FYp <- as.numeric(apply(FY[,-p], 1L, prod))

        # numeric
        if(p %in% num.idx) {
            y <- X[,p]
            sd.v <- sqrt(THETA[p,p])
            dy <- dnorm(y, mean=yhat[,p], sd=sd.v)

            pre <- 

            # lambda
            dlambda <- dy * 1/sd.v^2 * ((y - yhat[,p]) %*% ETA)
            dFYp <- dFYp +
                ( (dlambda %*% Delta.lambda[lambda.idx,,drop=FALSE]) * FYp)

            # theta
            dsigma2 <- dy * (1/(2*sd.v^4)*(y - yhat[,p])^2 - 1/(2*sd.v^2))
            dFYp <- dFYp +
                ( (dsigma2 %*% Delta.theta[theta.idx,,drop=FALSE]) * FYp)

            # nu
            dnu <-  dy * 1/sd.v^2 * (y - yhat[,p])
            dFYp <- dFYp +
               ( (dnu %*% Delta.mu[nu.idx,,drop=FALSE]) * FYp)

            # psi
            dpsi <- 0 # FIXME

        }

        # ordinal
        if(p %in% ord.idx) {

            # just in case we need theta[v,v] after all...
            sd.v.inv <- 1/sqrt(THETA[p,p])

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
            dth <- -1 * (Y2*p2 - Y1*p1) * sd.v.inv
            dFYp <- dFYp +
                ( (dth %*% Delta.th[th.idx==p,,drop=FALSE]) * FYp)

            # d lambda
            dlambda <- (-1 * (p1 - p2) * sd.v.inv) %*% ETA
            dFYp <- dFYp +
                ( (dlambda %*% Delta.lambda[lambda.idx,,drop=FALSE]) * FYp)

            # d S (=chol psi)
            dpsi <- (-1*(p1-p2)*sd.v.inv) %*% 
                        kronecker(LAMBDA[p,,drop=FALSE], KSI)
            dFYp <- dFYp +
                ( (dpsi %*% Delta.S) * FYp)
        }
    }

    dFYp
}

