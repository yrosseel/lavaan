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
estimator.WLS <- function(Sigma.hat=NULL, Mu.hat=NULL,
                          w.vecs=NULL, data.mean,
                          WLS.V=NULL,
                          meanstructure=FALSE) {

    if(!meanstructure) {
        obs <- w.vecs
        est <- vech(Sigma.hat)
    } else {
        obs <- c(data.mean, w.vecs)
        est <- c(Mu.hat, vech(Sigma.hat))
    }

    diff <- as.matrix(obs - est)
    fx <- as.numeric( t(diff) %*% WLS.V %*% diff )

    # no negative values
    if(fx < 0.0) fx <- 0.0

    fx
}

# Full Information ML estimator (FIML) handling the missing values 
estimator.FIML <- function(Sigma.hat=NULL, Mu.hat=NULL, M=NULL) {

    npatterns <- M$npatterns

    fx.p <- numeric(npatterns)
     w.p <- numeric(npatterns)

    # for each missing pattern, combine cases and compute raw loglikelihood
    for(p in 1:npatterns) {

        SX <- M$data[[p]][["SX"]]
        MX <- M$data[[p]][["MX"]]
        w.p[p] <- nobs <- M$data[[p]][["nobs"]]
        var.idx <- M$data[[p]][["var.idx"]]

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

    # no negative values
    if(fx < 0.0) fx <- 0.0

    fx
}

