# Factor extraction method(s)
# YR Feb 2020
#
# - ULS_corner only (for now)
# - just to get better starting values for ESEM 


# ULS estimation
#
# - but resulting in a upper-corner all zeroes LAMBDA matrix
# - not using eigenvalues/vectors, but minimizing the residuals
#   directly

# - should give the same results as MINRES (after an orthogonal transformation)
# - unless there are heywood cases; this function allows for negative variances!

lav_efa_extraction_uls_corner <- function(S, nfactors = 1L) {

    stopifnot(is.matrix(S))
    S <- unname(S)
    nvar <- nrow(S)

    # extract variances
    S.var <- diag(S)

    # convert to correlation matrix (ULS is not scale invariant!)
    R <- cov2cor(S)
    R.inv <- solve(R)

    LAMBDA <- matrix(0.1, nvar, nfactors) # cannot be zero
    LAMBDA[,1] <- R[,1] * 0.7
    THETA <- diag(nvar)
    if(nfactors > 1L) {
        corner.idx <- which(row(LAMBDA) < nfactors & col(LAMBDA) > row(LAMBDA))
        lambda.idx <- seq_len(nvar*nfactors)[-corner.idx]
        LAMBDA[corner.idx] <- 0
    } else {
        corner.idx <- integer(0L)
        lambda.idx <- seq_len(nvar)
    }

    # x are here the elements of LAMBDA (except the upper-right corner)
    obj.f.uls <- function(x) {
        LAMBDA[lambda.idx] <- x
        res <- lav_matrix_vech(R - tcrossprod(LAMBDA), diagonal = FALSE)^2
        sum(res)
    }

    # x are here the elements of LAMBDA (except the upper-right corner)
    grad.f.uls <- function(x) {
        LAMBDA[lambda.idx] <- x
        Sigma <- tcrossprod(LAMBDA)
        diag(Sigma) <- 1 # diagonal is ignored
        tmp <- -2 * (R - Sigma) %*% LAMBDA
        tmp[lambda.idx]
    }

    # optimize
    out <- nlminb(start = LAMBDA[lambda.idx], objective = obj.f.uls,
                  gradient = grad.f.uls, lower = -1, upper = +1)
    LAMBDA[lambda.idx] <- out$par
    diag(THETA) <- 1 - diag(tcrossprod(LAMBDA))

    # rescale if the input matrix was not a correlation matrix
    LAMBDA <- sqrt(S.var) * LAMBDA
    diag(THETA) <- S.var * diag(THETA)

    list(LAMBDA = LAMBDA, THETA = THETA)
}

                                                                        

