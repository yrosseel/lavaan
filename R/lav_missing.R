## This file contains various routines that are used when
## the data are incomplete
##
## initial verions Y.R. -- july 2010
##
## -- added EM algorithm: Y.R. aug 2011

# mle using EM
estimate.moments.EM <- function (X = NULL, M = NULL, verbose = FALSE,
                                 max.iter = 500L, tol = 1e-05) {

    if(verbose) {
        cat("\n")
        cat("estimation saturated H1 model -- start EM steps\n")
    }

    nvar <- ncol(X); pstar <- nvar * (nvar + 1)/2
    npatterns <- length(M)
    N <- nrow(X)

    # starting values as used by Mplus
    mu0  <- apply(X, 2, base::mean, na.rm = TRUE); names(mu0) <- NULL
    var0 <- apply(X, 2,  stats::var, na.rm = TRUE); names(var0) <- NULL
    sigma0 <- diag(x=var0, nrow=length(var0))
    mu <- mu0; sigma <- sigma0

    # report
    if(verbose) {
        fx0 <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=M)
        cat("  EM iteration:", sprintf("%4d", 0),
            " fx = ", sprintf("%15.10f", fx0),
            "\n")
    }

    # EM steps
    for(i in 1:max.iter) {
        T1 <- numeric(nvar)
        T2 <- matrix(0, nvar, nvar)
        for(p in 1:npatterns) {
            X       <- M[[p]][["X"]]
            MX      <- M[[p]][["MX"]]
            nobs    <- M[[p]][["nobs"]]
            var.idx <- M[[p]][["var.idx"]]

            if(all(var.idx)) {
                # complete pattern
                T1 <- T1 + colSums(X)
                T2 <- T2 + crossprod(X)
                next
            }

            # partition Mu (1=missing, 2=complete)
            Mu_1 <- mu[!var.idx]
            Mu_2 <- mu[ var.idx]

            # partition Sigma (1=missing, 2=complete)
            Sigma_11 <- sigma[!var.idx, !var.idx, drop=FALSE]
            Sigma_12 <- sigma[!var.idx,  var.idx, drop=FALSE]
            Sigma_21 <- sigma[ var.idx, !var.idx, drop=FALSE]
            Sigma_22 <- sigma[ var.idx,  var.idx, drop=FALSE]
            Sigma_22.inv <- inv.chol(Sigma_22, logdet=FALSE)
            #Sigma_22.inv <- solve(Sigma_22)

            # estimate missing values in this pattern
            Diff <- apply(X, 1, '-', Mu_2)
            X_missing2 <- t(Sigma_12 %*% Sigma_22.inv %*% Diff)
            X_missing <- t(apply(X_missing2, 1, '+', Mu_1))

            # complete data for this pattern
            X_complete <- matrix(0, nobs, nvar)
            X_complete[, var.idx] <- X
            X_complete[,!var.idx] <- X_missing

            # 1. SUM `completed' pattern
            T1_p <- colSums(X_complete)
            T1 <- T1 + T1_p

            # 2. CROSSPROD `completed' pattern
            T2_p <- crossprod(X_complete)

            # correction for missing cells: conditional covariances
            T2_p11 <- Sigma_11 - (Sigma_12 %*% Sigma_22.inv %*% Sigma_21)
            T2_p[!var.idx, !var.idx] <- T2_p[!var.idx, !var.idx] + (T2_p11*nobs)
            T2 <- T2 + T2_p
        }

        # M-step -- Little & Rubin (2000) page 225: eq. 11.6
        # recompute mu and sigma
        mu    <- T1/N
        sigma <- T2/N - tcrossprod(mu)

        # max absolute difference in parameter values
        DELTA <- max(abs(c(mu,vech(sigma)) - c(mu0,vech(sigma0))))

        # report fx
        if(verbose) {
            fx <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=M)
            cat("  EM iteration:", sprintf("%4d", i),
                " fx = ", sprintf("%15.10f", fx), 
                " delta par = ", sprintf("%9.8f", DELTA),
                "\n")
        }

        # convergence check: using parameter values:
        if(DELTA < tol)
            break

        # again
        mu0 <- mu; sigma0 <- sigma
    }

    # compute fx if we haven't already
    if(!verbose) 
        fx <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=M)

    if(verbose) {
        cat("estimated Sigma and Mu (H1):\n")
        cat("\nSigma:\n"); print(sigma)
        cat("\nMu:\n"); print(mu)
        cat("\n")
        cat("estimation saturated H1 model -- end\n\n")
    }

    # fx <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=M)
    list(sigma = sigma, mu = mu, fx = fx)
}

# get missing patterns for a single group (X is a matrix)
getMissingPatterns <- function(X) {

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
        MISSING <- MISSING[-empty.idx,,drop=FALSE]
              X <-       X[-empty.idx,,drop=FALSE]
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

    # return a list
    out <- list(nobs=ntotal, nvar=nvar,
                coverage=coverage, id=id, npatterns=npatterns,
                order=order, pat=pat, empty.idx=empty.idx)
    out
}

# construct summary statistics (cov/mean) per missing pattern
getMissingPatternStats <- function (X = NULL, Mp = NULL) {

    npatterns <- Mp$npatterns
    id        <- Mp$id   
    order     <- Mp$order
    pat       <- Mp$pat

    # prepare
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

    data
}


# estimate the `saturated' mu + sigma 
# using a quasi-newton method
# FIXME: this is a stub to get us going
#        we may prefer an EM-type algorithm here?
estimate.moments.fiml <- function (X = NULL, M = NULL, verbose = FALSE) {

    if(verbose) {
        cat("\n")
        cat("estimation saturated H1 model -- start\n")
    }

    nvar <- ncol(X); pstar <- nvar * (nvar + 1)/2

    # starting values
    start.cov <- cov(X, use = "p"); dimnames(start.cov) <- NULL
    start.mean <- apply(X, 2, base::mean, na.rm = TRUE)
    names(start.mean) <- NULL

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
        dx <- c(dx.out$dx.mu, vech(dx.out$dx.Sigma))
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
    start.x <- c(start.mean, vech(start.cov))

    # start iterations
    iter.max <- 500
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
        cat("estimation saturated H1 model -- end\n\n")
    }

    list(sigma = sigma, mu = mu, fx = fx)
}

