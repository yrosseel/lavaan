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
    if(is.null(lavmodel@control$cholesky)) {
        CHOLESKY <- TRUE
    } else {
        CHOLESKY <- as.logical(lavmodel@control$cholesky)
        #if(nfac > 1L && !CHOLESKY) {
        #    warning("lavaan WARNING: CHOLESKY is OFF but nfac > 1L")
        #}
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
        if(!is.null(MLIST$alpha) || !is.null(MLIST$gamma)) {
            EETAx <- computeEETAx.LISREL(MLIST = MLIST, eXo = eXo, N = nobs,
                        sample.mean = sample.mean,
                        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
                        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
                        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
                        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]])
            if(length(lv.dummy.idx) > 0L) {
                EETAx <- EETAx[,-lv.dummy.idx,drop=FALSE]
            }
        }
    }

    # compute case-wise likelihoods 
    lik <- lav_model_lik_mml(lavmodel = lavmodel, THETA = THETA, TH = TH,
               GLIST = GLIST, group = group, lavdata = lavdata,
               sample.mean = sample.mean, lavcache = lavcache)

    # prepare common stuff
    # fix Lambda?
    LAMBDA <- computeLAMBDA.LISREL(MLIST = MLIST,
                    ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
                    ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
                    ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
                    ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]])

    #### FIX th + mu!!!!!
    Delta.lambda <- computeDeltaDx(lavmodel, target="lambda")[[group]]
    if(!is.null(lavmodel@th.idx[[group]])) {
        Delta.th     <- computeDeltaDx(lavmodel, target="th"    )[[group]] ## FIXME? tau?
    }
    Delta.mu     <- computeDeltaDx(lavmodel, target="mu"    )[[group]] ## FIXME? nu?
    Delta.theta  <- computeDeltaDx(lavmodel, target="theta" )[[group]]
    Delta.beta   <- computeDeltaDx(lavmodel, target="beta"  )[[group]]
    ### FIXME -1?
    #Delta.beta <- -1*Delta.beta

    Delta.psi    <- computeDeltaDx(lavmodel, target="psi"   )[[group]]
    Delta.alpha  <- computeDeltaDx(lavmodel, target="alpha" )[[group]]
    Delta.gamma  <- computeDeltaDx(lavmodel, target="gamma" )[[group]]

    ## FIXME!!! do this analytically...
    #if(!is.null(chol.VETA)) {
        x <- lav_model_get_parameters(lavmodel = lavmodel, GLIST = MLIST)
        dVetadx <- function(x, lavmodel = fit@Model, g = 1L) {
            GLIST <- lavaan:::lav_model_x2GLIST(lavmodel, x=x, type="free")
            VETAx <- lavaan:::computeVETAx(lavmodel, GLIST = GLIST)[[g]]
            if(CHOLESKY) {
                S <- chol(VETAx)
            } else {
                S <- diag( sqrt(diag(VETAx)) )
            }
            S
        }
        Delta.S <- lavaan:::lavJacobianD(func=dVetadx, x=x, lavmodel = lavmodel,
                                         g = group)
    #} else {


    # remove dummy lv's from Delta.psi
    #lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[group]],
    #                  lavmodel@ov.x.dummy.lv.idx[[group]])
    #if(length(lv.dummy.idx) > 0L) {
    #    pstar <- nfac*(nfac+1)/2
    #    idx <- 1:pstar
    #    Delta.psi <- Delta.psi[idx,,drop=FALSE]
    #}

    #Delta.S <- Delta.psi
    # remove dummy lv's from Delta.psi
    #lv.dummy.idx <- c(lavmodel@ov.y.dummy.lv.idx[[group]],
    #                  lavmodel@ov.x.dummy.lv.idx[[group]])
    #if(length(lv.dummy.idx) > 0L) {
    #    pstar <- nfac*(nfac+1)/2
    #    idx <- 1:pstar
    #    Delta.S <- Delta.S[idx,,drop=FALSE]
    #}
    #}

    # compute dL/dx for each node
    dLdx <-  matrix(0, nGH, lavmodel@nx.free)
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
        #
        # - direct effect of BETA is already in VETAx, and hence chol.VETA
        # - need to add alpha, and GAMMA eta_i
        if(!is.null(MLIST$alpha) || !is.null(MLIST$gamma)) {
            eta <- sweep(EETAx, MARGIN=2, STATS=eta, FUN="+")
        }

        # again, compute yhat for this node (eta)
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
                          link = lavmodel@link, log. = TRUE)

        # FYp: FY prod minus p itself
        FY <- exp(log.fy.var) ### FIXME log/exp/log/...
        FYp <- matrix(0, nrow(FY), ncol(FY))
        for(p in 1:ncol(log.fy.var)) {
            FYp[,p] <- apply(FY[,-p], 1L, prod)
        }

        dFYp <- dFYp_x(X = X, eXo = eXo, yhat = yhat, MLIST = MLIST,
                       THETA = THETA, TH = TH, LAMBDA = LAMBDA,
                       Delta.lambda = Delta.lambda, 
                       Delta.th = Delta.th, Delta.mu = Delta.mu,
                       Delta.theta = Delta.theta, Delta.S = Delta.S,
                       Delta.psi = Delta.psi, Delta.beta = Delta.beta,
                       Delta.alpha = Delta.alpha, Delta.gamma = Delta.gamma,
                    ETA = eta, 
                    KSI = GH$x[q,,drop=FALSE],
                    chol.VETA = chol.VETA,
                    EETAx = EETAx,
                    FYp = FYp,
                    num.idx = lavmodel@num.idx[[group]],
                    th.idx  = lavmodel@th.idx[[group]],
                    ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
                    ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
                    ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
                    ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]])

        dLdx[q,] <- apply(1/lik * dFYp, 2, sum)
    }
  
    # integration
    dx <- apply(as.numeric(GH$w) * dLdx, 2, sum)

    # minimize
    dx <- -1*dx

    dx
}

dFYp_x <- function(X         = NULL,
                   eXo       = NULL,
                   yhat      = NULL,
                   MLIST     = NULL,
                   THETA     = NULL, 
                   TH        = NULL,
                   LAMBDA    = NULL,
                   Delta.lambda = NULL,
                   Delta.th     = NULL,
                   Delta.mu     = NULL,
                   Delta.theta  = NULL,
                   Delta.S      = NULL,
                   Delta.psi    = NULL,
                   Delta.beta   = NULL,
                   Delta.alpha  = NULL,
                   Delta.gamma  = NULL,
                   ETA       = NULL,
                   KSI       = NULL,
                   chol.VETA = NULL,
                   EETAx     = NULL,
                   FYp       = NULL,
                   num.idx   = NULL,
                   th.idx    = NULL,
                   ov.y.dummy.ov.idx = NULL,
                   ov.x.dummy.ov.idx = NULL,
                   ov.y.dummy.lv.idx = NULL,
                   ov.x.dummy.lv.idx = NULL) {

    ov.dummy.idx <- c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
    lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)
    nobs <- nrow(X); nvar <- ncol(X)

    # 'true' latent variables
    nfac <- ncol(MLIST$psi)
    if(length(lv.dummy.idx) > 0L) {
        nfac <- nfac - length(lv.dummy.idx)
    }

    # chol.VETA?
    #if(is.null(chol.VETA)) {
    #    ETA <- KSI
    #} else {
    #    ETA <- (KSI %*% chol.VETA)
    #    if(!is.null(EETAx)) {
    #        #ETA <- ETA + EETAx  ### FIXME? dim!!
    #        ETA <- sweep(EETAx, MARGIN=1, STATS=ETA, FUN="+")
    #    }
    #}

    # fix ALPHA
    ALPHA <- MLIST$alpha
    if(is.null(ALPHA)) {
        ALPHA <- numeric( nfac )
    } else if(length(lv.dummy.idx)) {
        ALPHA <- ALPHA[-lv.dummy.idx,,drop=FALSE]
    }

    # Beta?
    BETA <- MLIST$beta
    if(is.null(BETA)) {
        LAMBDA..IB.inv <- LAMBDA
    } else {
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
        LAMBDA..IB.inv <- MLIST$lambda %*% IB.inv ## no need to FIX???
        if(length(lv.dummy.idx) > 0L) {
            LAMBDA..IB.inv <- LAMBDA..IB.inv[,-lv.dummy.idx,drop=FALSE]
        }

        # fix BETA
        if(length(lv.dummy.idx)) {
            BETA <- MLIST$beta[-lv.dummy.idx, -lv.dummy.idx, drop=FALSE]
        }
        tmp <- -BETA; nr <- nrow(BETA); i <- seq_len(nr);
        tmp[cbind(i, i)] <- 1
        IB.inv <- solve(tmp)
    }

    # fix GAMMA
    GAMMA <- MLIST$gamma
    if(is.null(GAMMA)) {
        ALPHA.GAMMA.eXo <- matrix(as.numeric(ALPHA), nobs, nfac, byrow=TRUE)
    } else if(length(lv.dummy.idx)) {
        GAMMA <- GAMMA[-lv.dummy.idx,,drop=FALSE]
        ALPHA.GAMMA.eXo <- sweep(eXo %*% t(GAMMA),
                               MARGIN=2 ,STATS=as.numeric(ALPHA), FUN="+")
    }

    # ordered?
    ord.idx <- unique( th.idx[th.idx > 0L] )

    dFYp <- matrix(0, nobs, ncol(Delta.lambda))
    for(p in 1:nvar) {

        
        theta.idx <- diagh.idx(nvar)[p]

        # nu.idx
        nu.idx <- p
        D.mu <- Delta.mu[nu.idx,,drop=FALSE]

        # kappa idx
        if(p %in% ov.y.dummy.ov.idx) {
            nr <- nrow(MLIST$gamma); nc <- ncol(MLIST$gamma)
            ov.idx <- match(p, ov.y.dummy.ov.idx)
            pp <- ov.y.dummy.lv.idx[ov.idx]
            kappa.idx <- nr*((1:nc) - 1L) + pp
            D.kappa <- Delta.gamma[kappa.idx,,drop=FALSE]
        }

        # gamma.idx -- FIXME!!!!
        if(!is.null(MLIST$gamma)) {
            nr <- nrow(MLIST$gamma); nc <- ncol(MLIST$gamma)
            lv.idx <- 1:nfac
            # MUST BE ROWWISE!
            gamma.idx <- rep(nr*((1:nc) - 1L), times=length(lv.idx)) + rep(lv.idx, each=nc)
            D.gamma <- Delta.gamma[gamma.idx,,drop=FALSE]
        }

        # beta.idx
        if(!is.null(MLIST$beta)) {
            nr <- nc <- nrow(MLIST$beta)
            lv.idx <- 1:nfac
            # MUST BE ROWWISE!
            beta.idx <- rep(nr*((1:nfac) - 1L), times=nfac) + rep(lv.idx, each=nfac) 
            #beta.idx <- rep(nr*((1:nfac) - 1L), each=nfac) + lv.idx
            D.beta <- Delta.beta[beta.idx,,drop=FALSE]
        }

        # lambda idx
        if(p %in% ov.y.dummy.ov.idx) {
            nr <- nrow(MLIST$beta); nc <- nfac # only the first 1:nfac columns
            ov.idx <- match(p, ov.y.dummy.ov.idx)
            pp <- ov.y.dummy.lv.idx[ov.idx]
            beta.idx <- nr*((1:nc) - 1L) + pp
            D.lambda <- Delta.beta[beta.idx,,drop=FALSE]
        } else {
            lambda.idx <- nvar*((1:nfac) - 1L) + p
            D.lambda <- Delta.lambda[lambda.idx,,drop=FALSE]
        }

        # theta idx
        if(p %in% num.idx) {
            if(p %in% ov.y.dummy.ov.idx) {
                ov.idx <- match(p, ov.y.dummy.ov.idx)
                pp <- ov.y.dummy.lv.idx[ov.idx]
                psi.idx <- diagh.idx( ncol(MLIST$psi) )[pp]
                D.theta <- Delta.psi[psi.idx,,drop=FALSE]
            } else {
                theta.idx <- diagh.idx(nvar)[p]
                D.theta <-  Delta.theta[theta.idx,,drop=FALSE]
            }
        }

      

        # FYp
        fyp <- FYp[,p]

        # numeric
        if(p %in% num.idx) {
            y <- X[,p]
            sd.v <- sqrt(THETA[p,p])
            dy <- dnorm(y, mean=yhat[,p], sd=sd.v)

            # [nobs x 1]
            pre <- dy * 1/sd.v^2 * (y - yhat[,p])

            # lambda
            if(nrow(ETA) == 1L) {
                dlambda <- pre %*% ETA
            } else {
                dlambda <- pre * ETA
            }
            dFYp <- dFYp + ( (dlambda %*% D.lambda) * fyp)

            # theta
            dsigma2 <- dy * (1/(2*sd.v^4)*(y - yhat[,p])^2 - 1/(2*sd.v^2))
            dFYp <- dFYp + ( (dsigma2 %*% D.theta) * fyp)

            # nu
            dnu <- pre
            dFYp <- dFYp + ( (dnu %*% D.mu) * fyp)

            # psi
            if(nrow(KSI) == 1L) {
                dpsi <- pre %*% kronecker(LAMBDA[p,,drop=FALSE], KSI)
            } else {
                dpsi <- pre * kronecker(LAMBDA[p,,drop=FALSE], KSI)
            }
            dFYp <- dFYp + ( (dpsi %*% Delta.S) * fyp)
            
            # beta
            if(!is.null(BETA)) {
                tmp <- kronecker(LAMBDA[p,,drop=FALSE], ALPHA.GAMMA.eXo) %*%
                         t( kronecker(t(IB.inv), IB.inv) )
                         # kronecker(t(IB.inv), IB.inv)
                dbeta <- pre * tmp

                #dFYp <- dFYp + ((dpsi %*% Delta.S - dbeta %*% D.beta) * fyp)
                #dFYp <- dFYp - ((dpsi %*% Delta.S) * fyp)
                dFYp <- dFYp + ( (dbeta %*% D.beta) * fyp)
            }
            # kappa
            if(p %in% ov.y.dummy.ov.idx) {
                dkappa <- pre * eXo # [nobs x nexo]
                dFYp <- dFYp + ( (dkappa %*% D.kappa) * fyp)
            }

            # gamma
            if(!is.null(eXo)) {
                dgamma <- pre * kronecker(LAMBDA..IB.inv[p,,drop=FALSE], eXo)
                dFYp <- dFYp + ( (dgamma %*% D.gamma) * fyp)
            }
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

            # numeric(nobs)
            pre <- -1 * (p1 - p2) * sd.v.inv 

            # d tau
            # [nobx * n.th]
            dth <- -1 * (Y2*p2 - Y1*p1) * sd.v.inv
            dFYp <- dFYp +
                ( (dth %*% Delta.th[th.idx==p,,drop=FALSE]) * fyp)

            # d lambda
            # lambda
            if(nrow(ETA) == 1L) {
                dlambda <- pre %*% ETA
            } else {
                dlambda <- pre * ETA
            }
            dFYp <- dFYp + ( (dlambda %*% D.lambda) * fyp)

            # d S (=chol psi)
            if(nrow(KSI) == 1L) {
                dpsi <- pre %*% kronecker(LAMBDA[p,,drop=FALSE], KSI)
            } else {
                dpsi <- pre * kronecker(LAMBDA[p,,drop=FALSE], KSI)
            }
            dFYp <- dFYp + ( (dpsi %*% Delta.S) * fyp)

            # beta
            #if(!is.null(BETA)) {
            #    if(nrow(KSI) == 1L) {
            #        tmp.S <- kronecker(LAMBDA[p,,drop=FALSE], KSI) %*% Delta.S
            #        tmp1 <- kronecker(IB.inv, t(IB.inv))
            #        tmp2 <- kronecker(LAMBDA[p,,drop=FALSE], ALPHA.GAMMA.eXo)
            #        tmp.B <- tmp2 %*% t(tmp1)
            #        #tmp.B <- tmp2 %*% tmp1
            #    } else {
            #        stop("fix me!")
            #    }
            #    BB <- sweep(-(tmp.B %*% D.beta), MARGIN=2, 
            #                STATS=as.numeric(tmp.S %*% Delta.S), FUN="+")
            #    dFYp <- dFYp + ( pre * BB * fyp)
            #}

            # beta2
            if(!is.null(BETA)) {
                tmp <- kronecker(LAMBDA[p,,drop=FALSE], ALPHA.GAMMA.eXo) %*%
                        t( kronecker(t(IB.inv), IB.inv) )
                        # kronecker(t(IB.inv), IB.inv)
                dbeta <- pre * tmp

                #dFYp <- dFYp + ((dpsi %*% Delta.S - dbeta %*% D.beta) * fyp)
                #dFYp <- dFYp - ((dpsi %*% Delta.S) * fyp)
                dFYp <- dFYp + ( (dbeta %*% D.beta) * fyp)
            }



            # kappa
            if(p %in% ov.y.dummy.ov.idx) {
                dkappa <- pre * eXo
                dFYp <- dFYp + ( (dkappa %*% D.kappa) * fyp)
            }

            # gamma
            if(!is.null(eXo)) {
                dgamma <- pre * kronecker(LAMBDA..IB.inv[p,,drop=FALSE], eXo)
                dFYp <- dFYp + ( (dgamma %*% D.gamma) * fyp)
            }
        }
    }

    dFYp
}

