# full inverse of a square matrix using Gauss-Jordan elimination
# without pivoting (interchanging rows)
#
# not for real use; only to check the GaussJordanPivot function
# this function will fail if a pivot is (near)zero!
GaussJordanInverseNoPivoting <- function(A.) {

    A <- A. # we will change it
    p <- (d <- dim(A))[1L]
    if (!is.numeric(A) || length(d) != 2L || p != d[2L])
        stop("'A' is not a square numeric matrix")

    # the code below is based on Matlab code presented in Figure 1
    # of the following paper: 
    #
    # "Efficient Matrix Inversion via Gauss-Jordan Elimination and its
    #  Parralelization" (1998). 
    # Authors: Quintana, Quintana, Sun & van de Geijn
    # retrieved from www.cs.utexas.edu/users/plapack/papers/inverse-tr.ps
    for(k in 1:p) {
        pivot <- A[k,k]
        # column scaling
        A[-k,k] <- -A[-k,k]/pivot; A[k,k] <- 0.0
        # rank-1 update
        A <- A + tcrossprod(A[,k], A[k,])
        # row scaling
        A[k,-k] <- A[k,-k]/pivot; A[k,k] <- 1.0/pivot
    }

    A
}

# just for fun, *with* pivoting of rows
GaussJordanInverse <- function(A.) {

    A <- A.
    p <- (d <- dim(A))[1L]
    if (!is.numeric(A) || length(d) != 2L || p != d[2L])
        stop("'A' is not a square numeric matrix")

    ipiv <- 1:p
    for(k in 1:p) {
        # look for largest element in the column
        tmp <- A[,k]; if(k > 1L) tmp[1:(k-1)] <- 0.0; k.max <- which.max(tmp)
        pivot <- A[k.max,k]
        # pivoting rows?
        if(k != k.max) {
            tmp <- A[k,]; A[k,] <- A[k.max,]; A[k.max,] <- tmp
            tmp <- ipiv[k]; ipiv[k] <- ipiv[k.max]; ipiv[k.max] <- tmp
        }

        # column scaling
        A[-k,k] <- -A[-k,k]/pivot; A[k,k] <- 0.0
        # rank-1 update
        A <- A + tcrossprod(A[,k], A[k,])
        # row scaling
        A[k,-k] <- A[k,-k]/pivot; A[k,k] <- 1.0/pivot
    }

    # backward permutation
    A[,ipiv] <- A
    A
}

# perform a single 'Gauss-Jordan' pivot on a 
# square matrix 'A' using a diagonal element A[k,k]
# as a pivot element
GaussJordanPivot <- function(A., k=1L) {

    A <- A.
    p <- (d <- dim(A))[1L]
    if (!is.numeric(A) || length(d) != 2L || p != d[2L])
        stop("'A' is not a square numeric matrix")

    pivot <- A[k,k]
    # column scaling
    A[-k,k] <- -A[-k,k]/pivot; A[k,k] <- 0.0
    # rank-1 update
    A <- A + tcrossprod(A[,k], A[k,])
    # row scaling
    A[k,-k] <- A[k,-k]/pivot; A[k,k] <- 1.0/pivot

    A
}


# fabin3 (2sls) for a single factor (that is what Mplus is doing)
# S is a selection of ov.names corresponding with the indicators
# of the (single) factor

# 'Hagglund' version (slow)
fabin3.uni2 <- function(S) {
    # Gosta Hagglund (Psychometrika, 1982, 47)
    nvar <- ncol(S)
    if(nvar < 3) {
        return( rep(1, nvar) )
    }

    out <- numeric( nvar ); out[1L] <- 1.0
    for(i in 2:nvar) {
        idx3 <- (1:nvar)[-c(i, 1L)]
        s23 <- S[i, idx3]
        S31 <- S[idx3, 1L]; S13 <- S[1L, idx3]
        S33 <- S[idx3,idx3]
        out[i] <- ( s23 %*% solve(S33) %*% S31 %*%
                    solve(S13 %*% solve(S33) %*% S31) )
    }

    out
}

# 'Jennrich' version (fast)
fabin3.uni <- function(S) {
    nvar <- ncol(S)
    if(nvar < 3) {
        return( rep(1, nvar) )
    }

    # zero out references
    S[1L, 1L] <- 0.0
    S.tilde <- try(solve(S), silent = TRUE)
    if(inherits(S.tilde, "try-error")) {
        return( rep(1, nvar) )
    }

    out <- numeric( nvar ); out[1L] <- 1.0
    for(i in 2:nvar) {
        S.bar <- GaussJordanPivot(S.tilde[c(i,1L), c(i,1L)], k=1L)
        out[i] <- -S.bar[1,-1] - (S[i,1L] %*% S.bar[-1,-1])
    }

    out
}



fabin2 <- function(S, ref.idx=NULL) {
    # Gosta Hagglund (Psychometrika, 1982, 47)
    nvar <- ncol(S); nfac <- length(ref.idx)
    stopifnot(nvar >= 3*nfac)

    out <- matrix(0, nvar, nfac)
    for(i in 1:nvar) {
        if(i %in% ref.idx) {
            out[i, ref.idx == i] <- 1.0
            next
        }
        idx3 <- (1:nvar)[-c(i, ref.idx)]
        s23 <- S[i, idx3]
        S31 <- S[idx3, ref.idx]; S13 <- S[ref.idx, idx3]
        out[i,] <- s23 %*% S31 %*% solve(S13 %*% S31)
    }

    out
}

fabin3 <- function(S, ref.idx=NULL) {
    # Gosta Hagglund (Psychometrika, 1982, 47)
    nvar <- ncol(S); nfac <- length(ref.idx)
    stopifnot(nvar >= 3*nfac)

    out <- matrix(0, nvar, nfac)
    for(i in 1:nvar) {
        if(i %in% ref.idx) {
            out[i, ref.idx == i] <- 1.0
            next
        }
        idx3 <- (1:nvar)[-c(i, ref.idx)]
        s23 <- S[i, idx3]
        S31 <- S[idx3, ref.idx]; S13 <- S[ref.idx, idx3]
        S33 <- S[idx3,idx3]
        out[i,] <- ( s23 %*% solve(S33) %*% S31 %*% 
                     solve(S13 %*% solve(S33) %*% S31) )
    }

    out
}

# using the 'jennrich' computations (much faster!)
fabin3.jennrich <- function(S, ref.idx=NULL) {
    # Robert I. Jennrich (Psychometrika, 1987, 52)
    nvar <- ncol(S); nfac <- length(ref.idx)
    stopifnot(nvar >= 3*nfac)

    # zero out references
    S[ref.idx, ref.idx] <- 0.0
    S.tilde <- solve(S)

    out <- matrix(0, nvar, nfac)
    for(i in 1:nvar) {
        if(i %in% ref.idx) {
            out[i, ref.idx == i] <- 1.0
            next
        }
        S.bar <- GaussJordanPivot(S.tilde[c(i,ref.idx), c(i,ref.idx)], k=1L)
        out[i,] <- -S.bar[1,-1] - (S[i,ref.idx] %*% S.bar[-1,-1])
    }
      
    out
}





# Hagglund table 1
JoreskogLawley1968.COR <- matrix(c(
1.000,  0.411,  0.479, 0.401, 0.370, 0.393,  0.078, 0.389, 0.411,
0.411,  1.000,  0.463, 0.223, 0.198, 0.244, -0.042, 0.169, 0.324,
0.479,  0.463,  1.000, 0.231, 0.272, 0.357, -0.126, 0.153, 0.307,
0.401,  0.223,  0.231, 1.000, 0.659, 0.688,  0.215, 0.221, 0.256,
0.370,  0.198,  0.272, 0.659, 1.000, 0.649,  0.293, 0.279, 0.324,
0.393,  0.244,  0.357, 0.688, 0.649, 1.000,  0.226, 0.298, 0.294,
0.078, -0.042, -0.126, 0.215, 0.293, 0.226,  1.000, 0.602, 0.446,
0.389,  0.169,  0.153, 0.221, 0.279, 0.298,  0.602, 1.000, 0.630,
0.411,  0.324,  0.307, 0.256, 0.324, 0.294,  0.446, 0.630, 1.000),
9, 9, byrow=TRUE)

# round(fabin2(JoreskogLawley1968.COR, ref.idx=c(1,4,7)), 3)
# round(fabin3(JoreskogLawley1968.COR, ref.idx=c(1,4,7)), 3)



# Jennrich 
Emmett.1949 <- matrix(c(
  1.0, 0.523, 0.395, 0.471, 0.346, 0.426, 0.576, 0.434, 0.639,
0.523,   1.0, 0.479, 0.506, 0.418, 0.462, 0.547, 0.283, 0.645,
0.395, 0.479,   1.0, 0.355, 0.270, 0.254, 0.452, 0.219, 0.504,
0.471, 0.506, 0.355,   1.0, 0.691, 0.791, 0.443, 0.285, 0.505,
0.346, 0.418, 0.270, 0.691,   1.0, 0.679, 0.383, 0.149, 0.409,
0.426, 0.462, 0.254, 0.791, 0.679,   1.0, 0.372, 0.314, 0.472,
0.576, 0.547, 0.452, 0.443, 0.383, 0.372,   1.0, 0.385, 0.680,
0.434, 0.283, 0.219, 0.285, 0.149, 0.314, 0.385,   1.0, 0.470,
0.639, 0.645, 0.504, 0.505, 0.409, 0.472, 0.680, 0.470,   1.0),
9, 9, byrow=TRUE)
colnames(Emmett.1949) <- rownames(Emmett.1949) <- paste("x",1:9,sep="")

# round( fabin3(Emmett.1949, ref.idx=c(1,4,7)), 3) 

