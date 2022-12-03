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
#
# YR 17 oct 2022: - add lower/upper bounds for theta
#                 - use 'lambda' correction to ensure PSI is positive definite
#
lav_cfa_lambda2thetapsi <- function(lambda = NULL, S = NULL, S.inv = NULL,
                                    GLS = FALSE, bounds = TRUE) {
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

    # check bounds for theta
    if(bounds) {
        diagS <- diag(S)

        # nonnegative
        too.small.idx <- which(theta < 0)
        theta[too.small.idx] <- 0

        # not larger than diag(S)
        too.large.idx <- which(theta > diagS)
        if(length(too.large.idx) > 0L) {
            theta[too.large.idx] <- diagS[too.large.idx] * 0.99
        }
    }

    # psi
    diag.theta <- diag(theta)
    lambda <- try(lav_matrix_symmetric_diff_smallest_root(S, diag.theta),
                  silent = TRUE)
    if(inherits(lambda, "try-error")) {
        warning("lavaan WARNING: failed to compute lambda")
        SminTheta <- S - diag.theta # and hope for the best
    } else {
        N <- 20L # conservative lower bound, no need to change
        cutoff <- 1 + 1/(N-1) # 1.052632
        if(lambda < cutoff) {
            lambda.star <- lambda - 1/(N - 1)
            SminTheta <- S - lambda.star * diag.theta
        } else {
            SminTheta <- S - diag.theta
        }
    }
    PSI <- M %*% SminTheta %*% t(M) # Just like local SAM

    list(lambda = LAMBDA, theta = theta, psi = PSI)
}

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_cfa_fabin_internal <- function(lavmodel = NULL, lavsamplestats = NULL,
                                   lavpartable = NULL, lavpta = NULL,
                                   lavdata = NULL, lavoptions = NULL) {
    # no structural part!
    if(any(lavpartable$op == "~")) {
        stop("lavaan ERROR: FABIN estimator only available for CFA models")
    }
    # no BETA matrix! (i.e., no higher-order factors)
    if(!is.null(lavmodel@GLIST$beta)) {
        stop("lavaan ERROR: FABIN estimator not available for models the require a BETA matrix")
    }
    # no std.lv = TRUE for now
    if(lavoptions$std.lv) {
        stop("lavaan ERROR: FABIN estimator not available if std.lv = TRUE")
    }

    nblocks <- lav_partable_nblocks(lavpartable)
    stopifnot(nblocks == 1L) # for now
    b <- 1L
    sample.cov <- lavsamplestats@cov[[b]]
    nvar <- nrow(sample.cov)
    lv.names <- lavpta$vnames$lv.regular[[b]]
    nfac <- length(lv.names)
    marker.idx <- lavpta$vidx$lv.marker[[b]]
    lambda.idx <- which(names(lavmodel@GLIST) == "lambda")
    lambda.nonzero.idx <- lavmodel@m.free.idx[[lambda.idx]]
    # only diagonal THETA for now...
    # because if we have correlated residuals, we should remove the
    # corresponding variables as instruments before we estimate lambda...
    # (see MIIV)
    theta.idx <- which(names(lavmodel@GLIST) == "theta") # usually '2'
    m.theta <- lavmodel@m.free.idx[[theta.idx]]
    nondiag.idx <- m.theta[!m.theta %in% lav_matrix_diag_idx(nvar)]
    if(length(nondiag.idx) > 0L) {
        warning("lavaan WARNING: this implementation of FABIN does not handle correlated residuals yet!")
    }

    # 1. estimate LAMBDA
    if(lavoptions$estimator == "FABIN2") {
        LAMBDA <- lav_cfa_fabin2(S = sample.cov, marker.idx = marker.idx,
                                 lambda.nonzero.idx = lambda.nonzero.idx)
    } else {
        LAMBDA <- lav_cfa_fabin3(S = sample.cov, marker.idx = marker.idx,
                                 lambda.nonzero.idx = lambda.nonzero.idx)
    }

    # 2. simple ULS method to get THETA and PSI (for now)
    out <- lav_cfa_lambda2thetapsi(lambda = LAMBDA, S = sample.cov,
                                   S.inv = lavsamplestats@icov[[b]],
                                   GLS = FALSE)
    THETA <- diag(out$theta)
    PSI <- out$psi

    # 3. correlated residuals (if any) are just the difference between
    #    Sigma and S
    #if(length(nondiag.idx) > 0L) {
    #    Sigma <- LAMBDA %*% PSI %*% t(LAMBDA) + THETA
    #    THETA[nondiag.idx] <- (sample.cov - Sigma)[nondiag.idx]
    #}

    # store matrices in lavmodel@GLIST
    lavmodel@GLIST$lambda <- LAMBDA
    lavmodel@GLIST$theta  <- THETA
    lavmodel@GLIST$psi    <- PSI

    # extract free parameters only
    x <- lav_model_get_parameters(lavmodel)
    x
}

