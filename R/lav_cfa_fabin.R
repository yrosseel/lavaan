# FABIN = factor analysis by instrumental variables
# Hagglund 1982 (efa), 1986 (cfa)

lav_cfa_fabin2 <- function(S, marker.idx = NULL, lambda.nonzero.idx = NULL) {
    nvar <- ncol(S); nfac <- length(marker.idx)

    # overview of free/fixed
    LAMBDA <- matrix(0, nvar, nfac)
    LAMBDA[lambda.nonzero.idx] <- -1L

    lambda <- matrix(0, nvar, nfac)
    for(i in 1:nvar) {
        if(i %in% marker.idx) {
            lambda[i, marker.idx == i] <- 1.0
            next
        }
        free.idx <- LAMBDA[i,] == -1L
        idx3 <- (1:nvar)[-c(i, marker.idx)]
        s23 <- S[i, idx3]
        fac.idx <- marker.idx[free.idx]

        if(length(fac.idx) == 1L) { # most common scenario in CFA
            S31 <- S13 <- S[idx3, fac.idx]
            lambda[i, free.idx] <- sum(s23 * S31) / sum(S13 * S13)
        } else {
            S31 <- S[idx3, fac.idx, drop = FALSE]
            S13 <- S[fac.idx, idx3, drop = FALSE]
            lambda[i, free.idx] <-  solve(S13 %*% S31, drop(s23 %*% S31))
        }
    }

    lambda
}

lav_cfa_fabin3 <- function(S, marker.idx = NULL, lambda.nonzero.idx = NULL) {
    nvar <- ncol(S); nfac <- length(marker.idx)

    # overview of free/fixed
    LAMBDA <- matrix(0, nvar, nfac)
    LAMBDA[lambda.nonzero.idx] <- -1L

    S33.inv <- try(solve(S[-marker.idx, -marker.idx, drop = FALSE]), 
                   silent = TRUE)
    if(inherits(S33.inv, "try-error")) {
        warning("lavaan WARNING: fabin3 failed; switching to fabin2")
        return(lav_cfa_fabin2(S = S, marker.idx = marker.idx,
                              lambda.nonzero.idx= lambda.nonzero.idx))
    }

    lambda <- matrix(0, nvar, nfac)
    rm3.idx <- 0L
    for(i in 1:nvar) {
        if(i %in% marker.idx) {
            lambda[i, marker.idx == i] <- 1.0
            next
        }
        free.idx <- LAMBDA[i,] == -1L
        idx3 <- (1:nvar)[-c(i, marker.idx)]
        S33 <- S[idx3, idx3, drop = FALSE]
        s23 <- S[i, idx3]
        fac.idx <- marker.idx[free.idx]
        rm3.idx <- rm3.idx + 1L
        # update inverse
        s33.inv <- lav_matrix_symmetric_inverse_update(S.inv = S33.inv,
                                                       rm.idx = rm3.idx)

        if(length(fac.idx) == 1L) { # most common scenario in CFA
            S31 <- S13 <- S[idx3, fac.idx]
            tmp <- s33.inv %*% S31 # or colSums(s33.inv * S31)
            lambda[i, free.idx] <- sum(s23 * tmp) / sum(S13 * tmp)
        } else {
            S31 <- S[idx3, fac.idx, drop = FALSE]
            S13 <- S[fac.idx, idx3, drop = FALSE]
            tmp <- s33.inv %*% S31
            # lambda[i, free.idx] <- ( s23 %*% solve(S33) %*% S31 %*%
            #                          solve(S13 %*% solve(S33) %*% S31) )
            lambda[i, free.idx] <-  solve(S13 %*% tmp, drop(s23 %*% tmp))
        }
    }

    lambda
}

# compute THETA and PSI, given lambda using either ULS or GLS
# this function assumes:
# - THETA is diagonal
# - PSI is unrestricted
# - we assume W = S^{-1}
lav_cfa_lambda2thetapsi <- function(lambda = NULL, S = NULL, S.inv = NULL, 
                                    GLS = FALSE) {
    LAMBDA <- as.matrix(lambda)
    nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)

    if(GLS) {
        # see Browne, 1974 section 4 case II
        if(is.null(S.inv)) {
            W <- solve(S)
        } else {   
            W <- S.inv
        }
        tLW <- crossprod(LAMBDA, W)
        M <- solve(tLW %*% LAMBDA, tLW) # GLS mapping
        #D <- W %*% LAMBDA %*% M # symmmetric
        D <- crossprod(M, tLW)
        #theta <- solve(W*W - D*D, diag(W %*% S %*% W - D %*% S %*% D))
        theta <- solve(W*W - D*D, diag(W - D)) # because W == S^{-1}
    } else {
        # see Hagglund 1982, section 4
        M <- solve(crossprod(LAMBDA), t(LAMBDA)) # ULS mapping function
        D <- LAMBDA %*% M
        theta <-  solve(diag(nvar) - D*D, diag(S - (D %*% S %*% D)))
    }

    # psi
    SminTheta <- S
    diag.idx <- lav_matrix_diag_idx(nvar)
    SminTheta[diag.idx] <- SminTheta[diag.idx] - theta
    PSI <- M %*% SminTheta %*% t(M)

    list(lambda = LAMBDA, theta = theta, psi = PSI)
}
