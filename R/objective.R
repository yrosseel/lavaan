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
                          X  = NULL) {  # data

    # YR 3 okt 2012
    # the idea is to compute for each pair of variables, the model-based 
    # probability (or likelihood in mixed case) (that we observe the data 
    # for this pair under the model) for *each case* 
    # after taking logs, the sum over the cases gives the probablity/likelihood 
    # for this pair
    # the sum over all pairs gives the final PML based logl

    nvar <- nrow(Sigma.hat)
    pstar <- nvar*(nvar-1)/2
    ov.types <- rep("ordered", nvar)
    if(length(num.idx) > 0L) ov.types[num.idx] <- "numeric"

    #print(Sigma.hat); print(TH); print(th.idx); print(num.idx); print(str(X))

    LIK <- matrix(0, nrow(X), pstar) # likelihood per case, per pair
    PSTAR <- matrix(0, nvar, nvar)   # utility matrix, to get indices
    PSTAR[lavaan:::vech.idx(nvar, diag=FALSE)] <- 1:pstar
    PROW <- row(PSTAR)
    PCOL <- col(PSTAR)

    # shortcut for all ordered
    if(all(ov.types == "ordered")) {
        PROW <- row(PSTAR)
        PCOL <- col(PSTAR)
        NTH <- tabulate(th.idx)
        NCAT <- NTH + 1
        PTH <- pnorm(TH)
 
        #cat("SHORTCUT for ordinal\n")
        V12 <- as.matrix(expand.grid(th.idx, th.idx))
        V34 <- as.matrix(expand.grid(TH, TH))
        # remove vars and upper
        idx <- which(V12[,1] > V12[,2])
        V12 <- V12[idx,]; V34 <- V34[idx,]
        V5  <- Sigma.hat[ as.matrix(V12) ]
        B2 <- pbivnorm:::pbivnorm(V34[,1], V34[,2], rho=V5)
        B1a <- pnorm(V34[,1])
        B1b <- pnorm(V34[,2])

        # for all pairs, get table, get PI, get freq, get logl
        res <- sapply(seq_len(pstar), function(x) {
            m.idx <- which(PSTAR == x)
            i <- PROW[m.idx]; j <- PCOL[m.idx]
            idx <- which(V12[,1] == i & V12[,2] == j)
            BI <- B2[idx]; dim(BI) <- c(NTH[i], NTH[j])
            pth.y1 <- PTH[th.idx == i]
            pth.y2 <- PTH[th.idx == j]
            BI <- rbind(0, BI, pth.y2, deparse.level = 0)
            BI <- cbind(0, BI, c(0, pth.y1, 1), deparse.level = 0)
            # get probabilities
            nr <- nrow(BI); nc <- ncol(BI)
            PI <- BI[-1L,-1L] - BI[-1L,-nc] - BI[-nr,-1L] + BI[-nr,-nc]
            PI[PI < .Machine$double.eps] <- .Machine$double.eps
            LIK <- PI[ cbind(X[,i], X[,j]) ]
            if(any(LIK == 0.0)) {
                logl <- Inf
            } else {
                logl <- sum(log(LIK))
            }
        })
        # sum over all pairs
        LogLik <- sum(res) 
    
    } else {
        for(j in seq_len(nvar-1L)) {
            for(i in (j+1L):nvar) {
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
    fx <- -1 * LogLik

    fx
}

