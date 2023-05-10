# utility functions for the sam() function
# YR 4 April 2023

# construct 'mapping matrix' M using either "ML", "GLS" or "ULS" method
# optionally return MTM (for ML)
#
# by construction, M %*% LAMBDA = I (the identity matrix)
lav_sam_mapping_matrix <- function(LAMBDA = NULL, THETA = NULL,
                                   S = NULL, S.inv = NULL,
                                   method = "ML", warn = TRUE) {
    # ULS
    # M == solve( t(LAMBDA) %*% LAMBDA ) %*% t(LAMBDA)
    #   == MASS:::ginv(LAMBDA)
    if(method == "ULS") {

        # M == solve( t(LAMBDA) %*% LAMBDA ) %*% t(LAMBDA)
        #   == MASS:::ginv(LAMBDA)
        M <- try(tcrossprod(solve(crossprod(LAMBDA)), LAMBDA),
                 silent = TRUE)
        if(inherits(M, "try-error")) {
            if(warn) {
                warning("lavaan WARNING: cannot invert crossprod(LAMBDA); using generalized inverse")
            }
            M <- MASS::ginv(LAMBDA)
        }


    # GLS
    # M == solve( t(LAMBDA) %*% S.inv %*% LAMBDA ) %*% t(LAMBDA) %*% S.inv
    } else if(method == "GLS") {
        if(is.null(S.inv)) {
            S.inv <- try(solve(S), silent = TRUE)
        }
        if(inherits(S.inv, "try-error")) {
            if(warn) {
                warning("lavaan WARNING: S is not invertible; switching to ULS metho")
            }
            M <- lav_sam_mapping_matrix(LAMBDA = LAMBDA, method = "ULS")
        } else {
            tLSinv  <- t(LAMBDA) %*% S.inv
            tLSinvL <- tLSinv %*% LAMBDA
            M <- try(solve(tLSinvL, tLSinv), silent = TRUE)
            if(inherits(M, "try-error")) {
                if(warn) {
                    warning("lavaan WARNING: problem contructing mapping matrix; switching to generalized inverse")
                }
                M <- MASS::ginv(tLSinvL) %*% tLSinv
            }
        }

    # ML
    # M == solve(t(LAMBDA) %*% THETA.inv %*% LAMBDA) %*% t(LAMBDA) %*% THETA.inv
    } else if(method == "ML") {

        # Problem: if THETA has zero elements on the diagonal, we cannot invert

        # As we do not have access to Sigma(.inv), we cannot
        # use the trick as in lavPredict() (where we replace THETA.inv by
        # Sigma.inv)

        # old method (<0.6-16): remove rows/cols with zero values
        # on the diagonal of THETA, invert the submatrix and
        # set the '0' diagonal elements to one. This resulted in (somewhat)
        # distorted results.

        # new in 0.6-16: we use the Wall & Amemiya (2000) method using
        # the so-called 'T' transformation

        zero.theta.idx <- which(abs(diag(THETA)) < 1e-4) # be conservative
        if(length(zero.theta.idx) == 0L) {
            # ok, no zero diagonal elements: try to invert THETA
            if(lav_matrix_is_diagonal(THETA)) {
                THETA.inv <- diag( 1/diag(THETA), nrow = nrow(THETA) )
            } else {
                THETA.inv <- try(solve(THETA), silent = TRUE)
                if(inherits(THETA, "try-error")) {
                    THETA.inv <- NULL
                }
            }
        } else {
            THETA.inv <- NULL
        }

        # could we invert THETA?
        if(!is.null(THETA.inv)) {
            # ha, all is good; compute M the usual way
            tLTi  <- t(LAMBDA) %*% THETA.inv
            tLTiL <- tLTi %*% LAMBDA
            M <- try(solve(tLTiL, tLTi), silent = TRUE)
            if(inherits(M, "try-error")) {
                if(warn) {
                    warning("lavaan WARNING: problem contructing ML mapping matrix; switching to ULS")
                }
                M <- lav_sam_mapping_matrix(LAMBDA = LAMBDA, method = "ULS")
            }
        } else {
            # use W&A2000's method using the 'T' transformation
            M <- try(lav_sam_mapping_matrix_tmat(LAMBDA = LAMBDA,
                                               THETA = THETA), silent = TRUE)
            if(inherits(M, "try-error")) {
                if(warn) {
                    warning("lavaan WARNING: problem contructing ML mapping matrix; switching to ULS")
                }
                M <- lav_sam_mapping_matrix(LAMBDA = LAMBDA, method = "ULS")
            }
        }
    } # ML

    M
}

# use 'T' transformation to create the 'Bartlett/ML' mapping matrix
# see Wall & Amemiya (2000) eq (7)
# see also Fuller 1987 page 357 (where T is called H), and page 364

# although Fuller/W&A always assumed that THETA is diagonal,
# their method seems to work equally well for non-diagonal THETA
#
# in our implementation:
# - we do NOT reorder the rows of LAMBDA
# - if std.lv = TRUE, we first rescale to 'create' marker indicators
#   and then rescale back at the end
#
lav_sam_mapping_matrix_tmat <- function(LAMBDA     = NULL,
                                        THETA      = NULL,
                                        marker.idx = NULL,
                                        std.lv     = NULL) {

    LAMBDA <- as.matrix.default(LAMBDA)
    nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)

    # do we have marker.idx?
    if(is.null(marker.idx)) {
        # 'marker' indicator has a single non-zero element in a row
        marker.idx <- lav_utils_get_marker(LAMBDA = LAMBDA, std.lv = TRUE)
        if(any(is.na(marker.idx))) {
            stop("lavaan ERROR: no clear markers in LAMBDA matrix")
        }
    }

    # std.lv TRUE or FALSE?
    if(is.null(std.lv)) {
        std.lv <- FALSE
        if(any(diag(LAMBDA[marker.idx,,drop = FALSE]) != 1)) {
            std.lv <- TRUE
        }
    }

    # if std.lv = TRUE, rescale
    if(std.lv) {
        MARKER <- LAMBDA[marker.idx,,drop = FALSE]
        marker.inv <- 1/diag(MARKER)
        LAMBDA <- t(t(LAMBDA) * marker.inv)
    }

    # compute 'T' matrix
    TMAT <- lav_sam_tmat(LAMBDA = LAMBDA, THETA = THETA,
                         marker.idx = marker.idx)

    # ML mapping matrix
    M <- TMAT[marker.idx,,drop = FALSE]

    if(std.lv) {
        M <- M * marker.inv
    }

    M
}

# create 'T' matrix (tmat) for T-transformation
#
# Notes: - here we assume that LAMBDA has unity markers (no std.lv = TRUE)
#        - TMAT is NOT symmetric!
#        - Yc %*% t(TMAT) transforms the data in such a way that we get:
#           1)  Bartlett factor scores in the marker columns
#           2) 'V' values in the non-marker columns, where:
#               V = Yc - Yc[,marker.idx] %*% t(LAMBDA)
#
lav_sam_tmat <- function(LAMBDA     = NULL,
                         THETA      = NULL,
                         marker.idx = NULL) {

    LAMBDA <- as.matrix.default(LAMBDA)
    nvar <- nrow(LAMBDA); nfac <- ncol(LAMBDA)

    # do we have marker.idx?
    if(is.null(marker.idx)) {
        # 'marker' indicator has a single 1 element in a row
        marker.idx <- lav_utils_get_marker(LAMBDA = LAMBDA, std.lv = FALSE)
        if(any(is.na(marker.idx))) {
            stop("lavaan ERROR: no clear markers in LAMBDA matrix")
        }
    }

    # construct 'C' matrix
    C2 <- diag(nvar); C2[, marker.idx] <- -1 * LAMBDA
    C <- C2[-marker.idx, , drop = FALSE]

    # compute Sigma.ve and Sigma.vv
    Sigma.ve <- C %*% THETA
    # Sigma.vv <- C %*% THETA %*% t(C)
    Sigma.vv <- Sigma.ve %*% t(C)

    # construct 'Gamma' (and Gamma2) matrix
    #Gamma <- (t(Sigma.ve) %*% solve(Sigma.vv))[marker.idx,, drop = FALSE]
    Gamma <- try(t(solve(Sigma.vv, Sigma.ve)[, marker.idx, drop = FALSE]),
                 silent = TRUE)
    if(inherits(Gamma, "try-error")) {
        tmp <- t(Sigma.ve) %*% MASS::ginv(Sigma.vv)
        Gamma <- tmp[marker.idx, , drop = FALSE]
    }
    Gamma2 <- matrix(0, nfac, nvar)
    Gamma2[,-marker.idx] <- Gamma
    Gamma2[, marker.idx] <- diag(nfac)

    # transformation matrix 'T' (we call it here 'Tmat')
    Tmat <- matrix(0, nvar, nvar)
    Tmat[-marker.idx, ] <- C
    Tmat[ marker.idx, ] <- -Gamma2 %*% C2

    Tmat
}


# compute VETA
# - if alpha.correction == 0     -> same as local SAM (or MOC)
# - if alpha.correction == (N-1) -> same as FSR+Bartlett
lav_sam_veta <- function(M = NULL, S = NULL, THETA = NULL,
                         alpha.correction = 0L,
                         lambda.correction = TRUE,
                         N = 20L, extra = FALSE) {

    # MSM
    MSM <- M %*% S %*% t(M)

    # MTM
    MTM <- M %*% THETA %*% t(M)

    # new in 0.6-16: make sure MTM is pd
    # (otherwise lav_matrix_symmetric_diff_smallest_root will fail)
    MTM <- zapsmall(lav_matrix_symmetric_force_pd(MTM, tol = 1e-04))

    # apply small sample correction (if requested)
    if(alpha.correction > 0) {
        alpha.N1 <- alpha.correction / (N - 1)
        if(alpha.N1 > 1.0) {
            alpha.N1 <- 1.0
        } else if(alpha.N1 < 0.0) {
            alpha.N1 <- 0.0
        }
        MTM <- (1 - alpha.N1) * MTM
        alpha <- alpha.correction
    } else {
        alpha <- alpha.correction
    }

    if(lambda.correction) {
        # use Fuller (1987) approach to ensure VETA is positive
        lambda <- try(lav_matrix_symmetric_diff_smallest_root(MSM, MTM),
                      silent = TRUE)
        if(inherits(lambda, "try-error")) {
            warning("lavaan WARNING: failed to compute lambda")
            VETA <- MSM - MTM # and hope for the best
        } else {
            cutoff <- 1 + 1/(N-1)
            if(lambda < cutoff) {
                lambda.star <- lambda - 1/(N - 1)
                VETA <- MSM - lambda.star * MTM
            } else {
                VETA <- MSM - MTM
            }
        }
    } else {
        VETA <- MSM - MTM
    }

    # extra attributes?
    if(extra) {
        attr(VETA, "lambda") <- lambda
        attr(VETA, "alpha")  <- alpha
    }

    VETA
}

# compute EETA = E(Eta) = M %*% [YBAR - NU]
lav_sam_eeta <- function(M = NULL, YBAR = NULL, NU = NULL) {
    EETA <- M %*% (YBAR - NU)
    EETA
}

