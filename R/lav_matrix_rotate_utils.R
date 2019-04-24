# collection of functions that deal with rotation matrices
# YR  3 April 2019 -- initial version

# generate random orthogonal rotation matrix
lav_matrix_rotate_gen <- function(M = 10L, orthogonal = TRUE) {

    # catch M=1
    if(M == 1L) {
        return(matrix(1, 1, 1))
    }

    # create random normal matrix
    tmp <- matrix(rnorm(M*M), nrow = M, ncol = M)

    if(orthogonal) {
        # use QR decomposition
        qr.out <- qr(tmp)

        # extra 'Q' part
        out <- qr.Q(qr.out)
    } else {
        # just normalize *columns* of tmp -> crossprod(out) has 1 on diagonal
        out <- t( t(tmp) / sqrt(diag(crossprod(tmp))) )
    }

    out
}

# check if ROT is an orthogonal matrix if orthogonal = TRUE, or normal if
# orthogonal = FALSE
lav_matrix_rotate_check <- function(ROT = NULL, orthogonal = TRUE,
                                    tolerance = sqrt(.Machine$double.eps)) {

    # we assume ROT is a matrix
    M <- nrow(ROT)

    # crossprod
    RR <- crossprod(ROT)

    # target
    if(orthogonal) {
        # ROT^T %*% ROT = I
        target <- diag(M)
    } else {
        # diagonal should be 1
        target <- RR
        diag(target) <- 1
    }

    # compare for near-equality
    res <- all.equal(target = target, current = RR, tolerance = tolerance)

    # return TRUE or FALSE
    if(is.logical(res) && res) {
        out <- TRUE
    } else {
        out <- FALSE
    }

    out
}

# get weights vector needed to weight the rows using Kaiser normalization
lav_matrix_rotate_kaiser_weights <- function(A = NULL) {
    1 / sqrt( rowSums(A * A) )
}

# get weights vector needed to weight the rows using Cureton & Mulaik (1975)
# standardization
# see also Browne (2001) page 128-129
#
# Note: the 'final' weights are mutliplied by the Kaiser weights (see CEFA)
#
lav_matrix_rotate_cm_weights <- function(A = NULL) {

    P <- nrow(A); M <- ncol(A)

    # first principal component of AA'
    A.eigen <- eigen(tcrossprod(A), symmetric = TRUE)
    a <- A.eigen$vectors[,1] * sqrt(A.eigen$values[1])

    Kaiser.weights <- 1/sqrt( rowSums(A * A) )
    a.star <- abs(a * Kaiser.weights) # always between 0 and 1

    m.sqrt.inv <- 1/sqrt(M)
    acos.m.sqrt.inv <- acos(m.sqrt.inv)

    delta <- numeric(P)
    delta[a.star < m.sqrt.inv] <- pi/2

    tmp <- (acos.m.sqrt.inv - acos(a.star))/(acos.m.sqrt.inv - delta) * (pi/2)

    # add constant (see Cureton & Mulaik, 1975, page 187)
    cm <- cos(tmp) * cos(tmp) + 0.001

    # final weights = weighted by Kaiser weights
    cm * Kaiser.weights
}


