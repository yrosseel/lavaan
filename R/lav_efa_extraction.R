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
    #R.inv <- solve(R)

    # eigenvalue decomposition (to get starting values for LAMBDA)
    EV <- eigen(R, symmetric = TRUE)

    # extract first nfac components (assuming no measurement error)
    PC <- ( EV$vectors[, seq_len(nfactors), drop = FALSE] %*%
                diag(sqrt(EV$values[seq_len(nfactors)])) )

    # rotate to echelon pattern (see echelon() in GPArotation package)
    HEAD <- PC[seq_len(nfactors), , drop = FALSE]
    LAMBDA <- PC %*% solve(HEAD, t(chol(tcrossprod(HEAD))))

    THETA <- diag(nvar)
    if(nfactors > 1L) {
        corner.idx <- which(row(LAMBDA) < nfactors & col(LAMBDA) > row(LAMBDA))
        lambda.idx <- seq_len(nvar*nfactors)[-corner.idx]
        LAMBDA[corner.idx] <- 0
    } else {
        corner.idx <- integer(0L)
        lambda.idx <- seq_len(nvar)
    }

    # optim.method
    minObjective <- efa_extraction_uls_corner_min_objective
    minGradient  <- efa_extraction_uls_corner_min_gradient
    minHessian   <- NULL

    # create cache environment
    cache <- efa_extraction_uls_corner_init_cache(LAMBDA = LAMBDA,
                                         lambda.idx = lambda.idx, R = R)

    control.nlminb <- list(eval.max = 20000L, iter.max = 10000L,
                           trace = 0L, abs.tol=(.Machine$double.eps * 10))

    # optimize
    out <- nlminb(start = cache$theta, objective = minObjective,
                  gradient = minGradient, hessian = minHessian,
                  control = control.nlminb, lower = -1, upper = +1,
                  cache = cache)

    LAMBDA[lambda.idx] <- out$par
    diag(THETA) <- 1 - diag(tcrossprod(LAMBDA))

    # rescale if the input matrix was not a correlation matrix
    LAMBDA <- sqrt(S.var) * LAMBDA
    diag(THETA) <- S.var * diag(THETA)

    list(LAMBDA = LAMBDA, THETA = THETA)
}


efa_extraction_uls_corner_init_cache <- function(LAMBDA = NULL,
                                                 lambda.idx = NULL,
                                                 R = NULL,
                                                 parent = parent.frame()) {
    theta <- LAMBDA[lambda.idx]
    out <- list2env(list(LAMBDA = LAMBDA,
                         lambda.idx = lambda.idx,
                         R = R,
                         theta = theta),
                    parent = parent)
    out
}

efa_extraction_uls_corner_min_objective <- function(x, cache = NULL) {
    cache$theta <- x
    with(cache, {
        LAMBDA[lambda.idx] <- theta
        res1 <- lav_matrix_vech(R - tcrossprod(LAMBDA), diagonal = FALSE)
        res2 <- res1 * res1
        return(sum(res2))
    })
}

efa_extraction_uls_corner_min_gradient <- function(x, cache = NULL) {
    # check if x has changed
    if(!all(x == cache$theta)) {
        cache$theta <- x
        # nothing to do
    }
    with(cache, {
        LAMBDA[lambda.idx] <- theta
        Sigma <- tcrossprod(LAMBDA)
        diag(Sigma) <- 1 # diagonal is ignored
        tmp <- -2 * (R - Sigma) %*% LAMBDA
        return(tmp[lambda.idx])
    })
}

