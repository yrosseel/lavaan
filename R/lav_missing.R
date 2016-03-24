## This file contains various routines that are used when
## the data are incomplete
##
## initial verions Y.R. -- july 2010
##
## - added EM algorithm: Y.R. aug 2011


# mle using EM
estimate.moments.EM <- function (Y = NULL, Mp = NULL, Yp = NULL,
                                 verbose = FALSE, max.iter = 500L, tol = 1e-05) {

    if(verbose) {
        cat("\n")
        cat("estimation saturated H1 model -- start EM steps\n")
    }

    nvar <- ncol(Y); N <- nrow(Y)
    if(length(Mp$empty.idx) > 0L) {
        N <- N - length(Mp$empty.idx)
    }

    # starting values as used by Mplus
    mu0  <- apply(Y, 2, base::mean, na.rm = TRUE); names(mu0) <- NULL
    var0 <- apply(Y, 2,  stats::var, na.rm = TRUE); names(var0) <- NULL
    sigma0 <- diag(x=var0, nrow=length(var0))
    mu <- mu0; sigma <- sigma0

    # report
    if(verbose) {
        fx0 <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=Yp)
        cat("  EM iteration:", sprintf("%4d", 0),
            " fx = ", sprintf("%15.10f", fx0),
            "\n")
    }

    # EM steps
    for(i in 1:max.iter) {
        T1 <- numeric(nvar)
        T2 <- matrix(0, nvar, nvar)
        for(p in seq_len(Mp$npatterns)) {

            nobs    <- Yp[[p]]$freq
            var.idx <- Mp$pat[p,]

            # extract raw data for these cases
            X <- Y[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

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
            Sigma_22.inv <- try(inv.chol(Sigma_22, logdet=FALSE), silent = TRUE)
            if(inherits(Sigma_22.inv, "try-error")) {
                stop("lavaan ERROR: Sigma_22.inv cannot be inverted")
            }
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

        # check if sigma is near-pd (+ poor fix)
        ev <- eigen(sigma, symmetric = TRUE, only.values = TRUE)
        tol <- 1e-6 # FIXME!
        if(any(ev$values < tol)) {
            #too.small <- which( ev$values < tol )
            #ev$values[too.small] <- tol 
            #ev$values <- ev$values + tol
            #sigma <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)

            # ridge
            diag(sigma) <- diag(sigma) + max(diag(sigma))*1e-08
        }

        # max absolute difference in parameter values
        DELTA <- max(abs(c(mu, lav_matrix_vech(sigma)) - c(mu0, lav_matrix_vech(sigma0))))

        # report fx
        if(verbose) {
            fx <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=Yp)
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
        fx <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=Yp)

    if(verbose) {
        cat("estimated Sigma and Mu (H1):\n")
        cat("\nSigma:\n"); print(sigma)
        cat("\nMu:\n"); print(mu)
        cat("\n")
        cat("estimation saturated H1 model -- end\n\n")
    }

    # fx <- estimator.FIML(Sigma.hat=sigma, Mu.hat=mu, M=Yp)
    list(sigma = sigma, mu = mu, fx = fx)
}

