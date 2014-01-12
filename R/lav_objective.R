# fitting function for standard ML
estimator.ML <- function(Sigma.hat=NULL, Mu.hat=NULL, 
                         data.cov=NULL, data.mean=NULL,
                         data.cov.log.det=NULL,
                         meanstructure=FALSE) {

    # FIXME: WHAT IS THE BEST THING TO DO HERE??
    # CURRENTLY: return Inf  (at least for nlminb, this works well)
    # used to be: Sigma.hat <- force.pd(Sigma.hat) etc...
    if(!attr(Sigma.hat, "po")) return(Inf)


    Sigma.hat.inv     <- attr(Sigma.hat, "inv")
    Sigma.hat.log.det <- attr(Sigma.hat, "log.det")
    nvar <- ncol(Sigma.hat)

    if(!meanstructure) {
        fx <- (Sigma.hat.log.det + sum(data.cov * Sigma.hat.inv) -
               data.cov.log.det - nvar)
    } else {
        W.tilde <- data.cov + tcrossprod(data.mean - Mu.hat)
        fx <- (Sigma.hat.log.det + sum(W.tilde * Sigma.hat.inv) -
               data.cov.log.det - nvar)
    }

    # no negative values
    if(fx < 0.0) fx <- 0.0

    fx
}

# 'classic' fitting function for GLS, not used for now
estimator.GLS <- function(Sigma.hat=NULL, Mu.hat=NULL,
                          data.cov=NULL, data.mean=NULL, data.nobs=NULL,
                          meanstructure=FALSE) {

    W <- data.cov
    W.inv <- solve(data.cov)
    if(!meanstructure) {
        tmp <- ( W.inv %*% (W - Sigma.hat) )
        fx <- 0.5 * (data.nobs-1)/data.nobs * sum( tmp * t(tmp))
    } else {
        tmp <- W.inv %*% (W - Sigma.hat)
        tmp1 <- 0.5 * (data.nobs-1)/data.nobs * sum( tmp * t(tmp))
        tmp2 <- sum(diag( W.inv %*% tcrossprod(data.mean - Mu.hat) ))
        fx <- tmp1 + tmp2
    }

    # no negative values
    if(fx < 0.0) fx <- 0.0

    fx
}

# general WLS estimator (Muthen, Appendix 4, eq 99 single group)
estimator.WLS <- function(WLS.est=NULL, WLS.obs=NULL, WLS.V=NULL) {

    diff <- as.matrix(WLS.obs - WLS.est)
    fx <- as.numeric( t(diff) %*% WLS.V %*% diff )

    # no negative values
    if(fx < 0.0) fx <- 0.0

    fx
}

# Full Information ML estimator (FIML) handling the missing values 
estimator.FIML <- function(Sigma.hat=NULL, Mu.hat=NULL, M=NULL, h1=NULL) {

    npatterns <- length(M)

    fx.p <- numeric(npatterns)
     w.p <- numeric(npatterns)

    # for each missing pattern, combine cases and compute raw loglikelihood
    for(p in 1:npatterns) {

        SX <- M[[p]][["SX"]]
        MX <- M[[p]][["MX"]]
        w.p[p] <- nobs <- M[[p]][["nobs"]]
        var.idx <- M[[p]][["var.idx"]]

        # note: if a decent 'sweep operator' was available (in fortran)
        # we might win some time by 'updating' the inverse by sweeping
        # out the changed patterns... (but to get the logdet, we need
        # to do it one column at a time?)

        #cat("FIML: pattern ", p, "\n")
        #print(Sigma.hat[var.idx, var.idx])

        Sigma.inv <- inv.chol(Sigma.hat[var.idx, var.idx], logdet=TRUE)
        Sigma.log.det <- attr(Sigma.inv, "logdet")
        Mu <- Mu.hat[var.idx]

        TT <- SX + tcrossprod(MX - Mu)
        trace <- sum(Sigma.inv * TT)

        fx.p[p] <- Sigma.log.det + trace
    }

    fx <- weighted.mean(fx.p, w=w.p)

    # ajust for h1
    if(!is.null(h1)) {
        fx <- fx - h1

        # no negative values
        if(fx < 0.0) fx <- 0.0
    }

    fx
}

# pairwise maximum likelihood
# this is adapted from code written by Myrsini Katsikatsou
#
# some changes:
# - no distinction between x/y (ksi/eta)
# - loglikelihoods are computed case-wise
estimator.PML <- function(Sigma.hat = NULL,    # model-based var/cov/cor
                          TH        = NULL,    # model-based thresholds + means
                          th.idx    = NULL,    # threshold idx per variable
                          num.idx   = NULL,    # which variables are numeric
                          X         = NULL,    # raw data
                          cache     = NULL) {  # housekeeping stuff

    # YR 3 okt 2012
    # the idea is to compute for each pair of variables, the model-based 
    # probability (or likelihood in mixed case) (that we observe the data 
    # for this pair under the model) for *each case* 
    # after taking logs, the sum over the cases gives the probablity/likelihood 
    # for this pair
    # the sum over all pairs gives the final PML based logl

    # first of all: check if all correlations are within [-1,1]
    # if not, return Inf; (at least with nlminb, this works well)
    cors <- Sigma.hat[lower.tri(Sigma.hat)]

    #cat("[DEBUG objective\n]"); print(range(cors)); print(range(TH)); cat("\n")
    if(any(abs(cors) > 1)) {
        # question: what is the best approach here??
        return(+Inf) 
        #idx <- which( abs(cors) > 0.99 )
        #cors[idx] <- 0.99 # clip
        #cat("CLIPPING!\n")
    }

    nvar <- nrow(Sigma.hat)
    pstar <- nvar*(nvar-1)/2
    ov.types <- rep("ordered", nvar)
    if(length(num.idx) > 0L) ov.types[num.idx] <- "numeric"
    #print(Sigma.hat); print(TH); print(th.idx); print(num.idx); print(str(X))

    LIK <- matrix(0, nrow(X), pstar) # likelihood per case, per pair
    PSTAR <- matrix(0, nvar, nvar)   # utility matrix, to get indices
    PSTAR[vech.idx(nvar, diagonal=FALSE)] <- 1:pstar
    PROW <- row(PSTAR)
    PCOL <- col(PSTAR)

    # shortcut for all ordered - tablewise
    if(all(ov.types == "ordered")) {
        # prepare for Myrsini's vectorization scheme
        LONG2 <- LongVecTH.Rho(no.x               = nvar,
                               all.thres          = TH,
                               index.var.of.thres = th.idx, 
                               rho.xixj           = cors)
        # get expected probability per table, per pair
        PI <- pairwiseExpProbVec(ind.vec = cache$LONG, th.rho.vec=LONG2)
        # get frequency per table, per pair
        #LogLik <- sum(cache$bifreq * log(PI))
    
        # more convenient fit function
        prop <- cache$bifreq / cache$nobs
        # remove zero props # FIXME!!! or add 0.5???
        zero.idx <- which(prop == 0.0)
        if(length(zero.idx) > 0L) {
            prop <- prop[-zero.idx]
            PI   <- PI[-zero.idx]
        }
        Fmin <- sum( prop*log(prop/PI) )

    } else {
        # # order! first i, then j, vec(table)!
        for(i in seq_len(nvar-1L)) {
            for(j in (i+1L):nvar) {
                pstar.idx <- PSTAR[i,j]
                # cat("pstar.idx =", pstar.idx, "i = ", i, " j = ", j, "\n")
                if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                    # ordinary pearson correlation
                    stop("not done yet")
                } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                    # polyserial correlation
                    stop("not done yet")
                } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                    # polyserial correlation
                    stop("not done yet")
                } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                    # polychoric correlation
                    PI <- pc_PI(rho   = Sigma.hat[i,j], 
                                th.y1 = TH[ th.idx == i ],
                                th.y2 = TH[ th.idx == j ])
                    LIK[,pstar.idx] <- PI[ cbind(X[,i], X[,j]) ]
                }
                #cat("Done\n")
            }
        }

        # check for zero likelihoods/probabilities
        # FIXME: or should we replace them with a tiny number?
        if(any(LIK == 0.0)) return(Inf) # we minimize
 
        # loglikelihood
        LogLIK.cases <- log(LIK)
 
        # sum over cases
        LogLIK.pairs <- colSums(LogLIK.cases)

        # sum over pairs
        LogLik <- sum(LogLIK.pairs)
    }

    # function value as returned to the minimizer
    #fx <- -1 * LogLik
    fx <- Fmin

    fx
}

# full information maximum likelihood
# underlying multivariate normal approach (see Joreskog & Moustaki, 2001)
#
estimator.FML <- function(Sigma.hat = NULL,    # model-based var/cov/cor
                          TH        = NULL,    # model-based thresholds + means
                          th.idx    = NULL,    # threshold idx per variable
                          num.idx   = NULL,    # which variables are numeric
                          X         = NULL,    # raw data
                          cache     = NULL) {  # patterns

    # YR 27 aug 2013
    # just for fun, and to compare with PML for small models

    # first of all: check if all correlations are within [-1,1]
    # if not, return Inf; (at least with nlminb, this works well)
    cors <- Sigma.hat[lower.tri(Sigma.hat)]

    if(any(abs(cors) > 1)) {
        return(+Inf) 
    }

    nvar <- nrow(Sigma.hat)
    pstar <- nvar*(nvar-1)/2
    ov.types <- rep("ordered", nvar)
    if(length(num.idx) > 0L) ov.types[num.idx] <- "numeric"
    MEAN <- rep(0, nvar)

    # shortcut for all ordered - per pattern
    if(all(ov.types == "ordered")) {
        PAT <- cache$pat; npatterns <- nrow(PAT)
        freq <- as.numeric( rownames(PAT) )
        PI <- numeric(npatterns)
        TH.VAR <- lapply(1:nvar, function(x) c(-Inf, TH[th.idx==x], +Inf))
        # FIXME!!! ok to set diagonal to 1.0?
        diag(Sigma.hat) <- 1.0
        for(r in 1:npatterns) {
            # compute probability for each pattern
            lower <- sapply(1:nvar, function(x) TH.VAR[[x]][ PAT[r,x]      ])
            upper <- sapply(1:nvar, function(x) TH.VAR[[x]][ PAT[r,x] + 1L ])


            # how accurate must we be here???
            PI[r] <- sadmvn(lower, upper, mean=MEAN, varcov=Sigma.hat,
                            maxpts=10000*nvar, abseps = 1e-07)
        }
        # sum (log)likelihood over all patterns
        #LogLik <- sum(log(PI) * freq)

        # more convenient fit function
        prop <- freq/sum(freq)
        # remove zero props # FIXME!!! or add 0.5???
        zero.idx <- which(prop == 0.0)
        if(length(zero.idx) > 0L) {
            prop <- prop[-zero.idx]
            PI   <- PI[-zero.idx]
        }
        Fmin <- sum( prop*log(prop/PI) )

    } else { # case-wise
        PI <- numeric(nobs)
        for(i in 1:nobs) {
            # compute probability for each case
            PI[i] <- stop("not implemented")
        }
        # sum (log)likelihood over all observations
        LogLik <- sum(log(PI))
        stop("not implemented")
    }

    # function value as returned to the minimizer
    #fx <- -1 * LogLik
    fx <- Fmin

    fx
}

estimator.MML <- function(object    = NULL,    # object
                          GLIST     = NULL,
                          g         = 1L,      # group
                          link      = "logit",
                          sample.mean = NULL,
                          X         = NULL,    # raw data
                          cache     = NULL) {  # patterns


    # compute likelihood for every observation
    lik <- lav_est_MML(object = object, GLIST = GLIST, g = g, link = link, 
                       sample.mean = sample.mean, X = X, cache = cache)

    # log + sum over observations
    logl <- sum( log(lik) )

    # function value as returned to the minimizer
    fx <- -logl

    fx
}

