## This file contains various routines that are used when
## the data are incomplete
##
## initial verions Y.R. -- july 2010


# construct summary information of missing patterns
missing.patterns <- function (X, warn=FALSE) {

    ntotal <- nrow(X); nvar <- ncol(X)

    # create a 1/0 matrix where the 1's denote the missing values
    MISSING <- 1L * is.na(X)

    # compute pairwise coverage (for reporting)
    coverage <- crossprod(1 - MISSING)/ntotal

    # identify, label and sort missing patterns 
    id <- apply(MISSING, MARGIN = 1, 
                function(x) {
                    if(sum(x) == length(x)) {
                        out <- "empty"
                    } else {
                        paste(x, collapse = "")
                    }
                }
               )
     
    # identify and remove empty row
    empty.idx <- which(id == "empty")
    if(length(empty.idx) > 0) {
        if(warn) {
            warning("lavaan WARNING: some cases are empty and will be removed")
            cat("lavaan WARNING: empty cases: "); print( empty.idx )
        }
        MISSING <- MISSING[-empty.idx,]
              X <-       X[-empty.idx,]
             id <-      id[-empty.idx]
        # adjust ntotal
        ntotal <- ntotal - length(empty.idx)
    }

    # sort patterns (from high occurence to low occurence)
    TABLE <- sort(table(id), decreasing = TRUE)
    order <- names(TABLE)
    npatterns <- length(TABLE)
    pat <- 1L - MISSING[match(order, id), , drop = FALSE]
    storage.mode(pat) <- "logical"
    row.names(pat) <- as.character(TABLE)
    data <- vector("list", length = npatterns)

    # fill in pattern information
    for (p in 1:npatterns) {
        row.idx <- which(id == order[p])
        nobs <- length(row.idx)
        Xp <- X[row.idx, pat[p, ], drop = FALSE]
        if (nobs > 1) {
            M <- colMeans(Xp)
            S <- crossprod(Xp)/nobs - tcrossprod(M)
        }
        else {
            S <- 0
            M <- as.numeric(Xp)
        }
        data[[p]] <- list(X = Xp, SX = S, MX = M, nobs = nobs, 
            var.idx = pat[p, ])
    }

    # return a list
    out <- list(nobs=ntotal, nvar=nvar, data=data, 
                coverage=coverage, id=id, npatterns=npatterns,
                order=order, pat=pat)
    out
}


# estimate the `saturated' mu + sigma 
# using a quasi-newton method
# FIXME: this is a stub to get us going
#        we may prefer an EM-type algorithm here?
estimate.moments.fiml <- function (X = NULL, M = NULL, verbose = FALSE) {

    if(verbose) {
        cat("\n")
        cat("lavaan Sample: estimation saturated H1 model -- start\n")
    }

    nvar <- ncol(X); pstar <- nvar * (nvar + 1)/2

    # starting values
    start.cov <- force.pd(cov(X, use = "p")); dimnames(start.cov) <- NULL
    start.mean <- apply(X, 2, mean, na.rm = TRUE); names(start.mean) <- NULL

    # x2param
    lower.idx <- which(lower.tri(start.cov, diag = TRUE))
    upper.idx <- which(upper.tri(t(start.cov), diag = TRUE))
    x2param <- function(x) {
        mu <- x[1:nvar]
        sigma.el <- x[-(1:nvar)]
        sigma <- matrix(0, nvar, nvar)
        sigma[lower.idx] <- sigma.el
        sigma[upper.idx] <- t(sigma)[upper.idx]
        list(mu = mu, sigma = sigma)
    }

    # objective function
    minimize.this.function <- function(x, verbose = TRUE) {
        out <- x2param(x)
        ev <- eigen(out$sigma)$values
        if(any(ev < 0.0)) { # too strict?
            return(Inf)
            #cat("lavaan WARNING: small or negative eigenvalues in objective function for H1 estimation\n")
            #print(eigen(out$sigma)$values)
            #out$sigma <- force.pd(out$sigma)
        }
        fx <- estimator.FIML(Sigma.hat=out$sigma, Mu.hat=out$mu, M=M)
        if (verbose) {
            cat("Objective function  = ", sprintf("%10.8f", fx), 
                "\n", sep = "")
        }
        fx
    }

    # gradient - analytical version
    first.derivative.param <- function(x, verbose = FALSE) {
        out <- x2param(x)
        dx.out <- derivative.FIML(Sigma.hat=out$sigma, Mu.hat=out$mu, M=M)
        dx <- c(dx.out$dx.mu, vecs(dx.out$dx.Sigma))
        dx
    }

    # gradient - numerical version (to check; not used)
    first.derivative.param.numerical <- function(x, verbose = FALSE) {
        npar <- length(x)
        h <- 1e-04
        dx <- numeric(npar)
        fx <- minimize.this.function(x)
        for (i in 1:npar) {
            x.right <- x
            x.right[i] <- x.right[i] + h
            fx.right <- minimize.this.function(x.right)
            dx[i] <- (fx.right - fx)/h
        }
        print(dx)
        stop("for now")
        dx
    }

    # get staring values
    start.x <- c(start.mean, vecs(start.cov))

    # start iterations
    iter.max <- 100
    optim.out <- nlminb(start=start.x, objective=minimize.this.function, 
                        gradient=first.derivative.param, 
                        control=list(iter.max=iter.max, eval.max=iter.max*2, 
                                     trace=0), 
                        verbose = verbose)
    if(verbose) {
        cat("convergence status (0=ok): ", optim.out$convergence, "\n")
        cat("optim/nlminb message says: ", optim.out$message, "\n")
        cat("number of iterations: ", optim.out$iterations, "\n")
        cat("number of function evaluations [objective, gradient]: ", 
            optim.out$evaluations, "\n")
    }

    # collect information and return results
    x <- optim.out$par
    fx <- optim.out$objective
    out <- x2param(x); sigma <- out$sigma; mu <- out$mu

    if(verbose) {
        cat("estimated Sigma and Mu (H1):\n")
        cat("\nSigma:\n"); print(sigma)
        cat("\nMu:\n"); print(mu)
        cat("\n")
        cat("lavaan Sample: estimation saturated H1 model -- end\n\n")
    }

    list(sigma = sigma, mu = mu, fx = fx)
}

