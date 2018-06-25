
# derivatives FIML unrestricted model (h1)
# pattern-based method
derivative.FIML <- function(Sigma.hat, Mu.hat, M) {

    ntotal    <- sum(sapply(M, "[[", "freq"))
    nvar      <- length(Mu.hat)
    npatterns <- length(M)
    dx.Sigma <- matrix(0, nvar, nvar)
    dx.Mu    <- matrix(0, nvar, 1)

    for(p in 1:npatterns) {
        SX <- M[[p]][["SY"]]
        MX <- M[[p]][["MY"]]
        nobs <- M[[p]][["freq"]]
        var.idx <- M[[p]][["var.idx"]]

        Sigma.inv <- inv.chol(Sigma.hat[var.idx, var.idx],
                              logdet=FALSE)
        Mu <- Mu.hat[var.idx]
        TT <- SX + tcrossprod(MX - Mu)

        dx.Mu[var.idx, 1] <-  ( dx.Mu[var.idx, 1] + nobs/ntotal *
                                  -2 * t(t(MX - Mu) %*% Sigma.inv) )

        dx.Sigma[var.idx, var.idx] <-
            ( dx.Sigma[var.idx, var.idx] - nobs/ntotal * 2 *
              # in the 'textbook' formula's, the Sigma.inv below is often
              # replaced by [0.5 * D'(Sigma.inv %x% Sigma.inv) D]
              # but we do not use the 'vecs' notation here, and
              # we 'compensate' for the symmetry later on
              (Sigma.inv %*%
              (TT - Sigma.hat[var.idx,var.idx]) %*% Sigma.inv ) )
    }

    # compensate for symmetry
    diag(dx.Sigma) <- diag(dx.Sigma)/2

    out <- list(dx.mu=dx.Mu, dx.Sigma=dx.Sigma)
    out

}


# X  <- matrix(rnorm(200*3), ncol=3)
# X2 <- matrix(rnorm(200*3), ncol=3)
# Sigma <- cov(X2); Mu <- colMeans(X2)

# logl of the MVM, no constants, factor 0.5, factor -1
# summary statistics version
logl.MVN.complete <- function(Sigma, Mu,
                              X=NULL, data.cov=NULL, data.mean=NULL) {

    if(is.null(data.cov)) {
        stopifnot(!is.null(X))
        nobs <- nrow(X)
        data.cov <- cov(X) * (nobs-1)/nobs
        data.mean <- colMeans(X)
    }

    Sigma.inv <- inv.chol(Sigma, logdet=TRUE)
    Sigma.log.det <- attr(Sigma.inv, "logdet")

    diff <- as.matrix(data.mean - Mu)
    TT <- data.cov + tcrossprod(diff)

    logl <- Sigma.log.det + sum(TT * Sigma.inv) # - S.log.det - nvar
    logl <- 0.5 * logl

    logl
}


# logl of the MVM, no constant, factor 0.5, factor -1, factor 1/nobs
# case-wise version
logl.MVN.casewise <- function(Sigma, Mu, X) {

    Sigma.inv <- inv.chol(Sigma, logdet=TRUE)
    Sigma.log.det <- attr(Sigma.inv, "logdet")

    tmp1 <- tmp2 <- logl <- 0.0
    for(i in 1:nrow(X)) {
        diff <- as.matrix(X[i,] - Mu)
        Tmp1 <- Sigma.log.det
        Tmp2 <- (t(diff) %*% Sigma.inv %*% diff)
        tmp1 <- tmp1 + Tmp1
        tmp2 <- tmp2 + Tmp2
    }

    logl <- -1 * as.numeric( -0.5 * tmp1  -0.5* tmp2 ) * 1/nrow(X)

    logl
}

# numerical/analytical hessian saturated model under MVN
hessian.MVN.saturated <- function(Sigma=NULL, Mu=NULL,
                                  X=NULL, data.cov=NULL, data.mean=NULL,
                                  meanstructure=TRUE,
                                  analytical=TRUE) {

    if(is.null(data.cov)) {
        stopifnot(!is.null(X))
        nobs <- nrow(X)
        data.cov <- cov(X) * (nobs-1)/nobs
        data.mean <- colMeans(X)
    }
    nvar <- ncol(data.cov)

    #if(!analytical) {
    #    lower.idx <- which(lower.tri(Sigma, diag = TRUE))
    #    upper.idx <- which(upper.tri(t(Sigma), diag = TRUE))
    #
    #    x2param <- function(x) {
    #        if(meanstructure) {
    #            mu <- x[1:nvar]
    #            sigma.el <- x[-(1:nvar)]
    #        } else {
    #            mu <- data.mean
    #            sigma.el <- x
    #        }
    #        sigma <- matrix(0, nvar, nvar)
    #        sigma[lower.idx] <- sigma.el
    #        sigma[upper.idx] <- t(sigma)[upper.idx]
    #        list(mu = mu, sigma = sigma)
    #    }
    #
    #    param2x <- function(Sigma, Mu) {
    #        if(meanstructure) {
    #            x1 <- as.numeric(Mu)
    #            x2 <- lav_matrix_vech(Sigma)
    #            x <- c(x1, x2)
    #        } else {
    #            x <- lav_matrix_vech(Sigma)
    #        }
    #        x
    #    }
    #
    #    objective.function <- function(x) {
    #        out <- x2param(x)
    #        mu.local <- out$mu
    #        sigma.local <- out$sigma
    #        fx <- logl.MVN.complete(Sigma=sigma.local, Mu=mu.local,
    #                                data.cov=data.cov, data.mean=data.mean)
    #        fx
    #    }
    #
    #    # compute numerical approximation of the Hessian
    #    H <- numDeriv::hessian(func=objective.function, x=param2x(Sigma,Mu))
    #
    #} else {

        Sigma.inv <- inv.chol(Sigma, logdet=FALSE)

        if(meanstructure) {
            diff <- as.matrix(data.mean - Mu)
            TT <- data.cov + tcrossprod(diff)

            H11 <- inv.chol(Sigma, logdet=FALSE)

            tmp <- t(diff) %*% Sigma.inv
            H12 <- lav_matrix_duplication_post(Sigma.inv %x% tmp)
            H21 <- t(H12)

            tmp <- (Sigma.inv %*% TT %*% Sigma.inv) - 0.5*Sigma.inv
            H22 <- lav_matrix_duplication_pre_post(Sigma.inv %x% tmp)

            H <- rbind( cbind(H11, H12),
                        cbind(H21, H22) )
        } else {
           TT <- data.cov
           tmp <- (Sigma.inv %*% TT %*% Sigma.inv) - 0.5*Sigma.inv
           H <- lav_matrix_duplication_pre_post(Sigma.inv %x% tmp)
        }
    # }

    H
}


