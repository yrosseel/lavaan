# Magnus & Neudecker (1999) style matrix operations
# YR - 11 may 2011: initial version
# YR - 19 okt 2014: rename functions using lav_matrix_ prefix

# vec operator
#
# the vec operator (for 'vectorization') transforms a matrix into 
# a vector by stacking the *columns* of the matrix one underneath the other
#
# M&N book: page 30
#
# note: we do not coerce to 'double/numeric' storage-mode (like as.numeric)
lav_matrix_vec <- function(A) {
    as.vector(A)
}


# vecr operator
#
# the vecr operator ransforms a matrix into 
# a vector by stacking the *rows* of the matrix one underneath the other
lav_matrix_vecr <- function(A) {
    lav_matrix_vec(t(A))
}


# vech 
# 
# the vech operator (for 'half vectorization') transforms a *symmetric* matrix 
# into a vector by stacking the *columns* of the matrix one underneath the 
# other, but eliminating all supradiagonal elements
#
# see Henderson & Searle, 1979
#
# M&N book: page 48-49
#
lav_matrix_vech <- function(S, diagonal = TRUE) {
    ROW <- row(S); COL <- col(S)
    if(diagonal) S[ROW >= COL] else S[ROW > COL]
}


# the vechr operator transforms a *symmetric* matrix 
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all supradiagonal elements
lav_matrix_vechr <- function(S, diagonal = TRUE) {
    S[lav_matrix_vechr_idx(n = ncol(S), diagonal = diagonal)]
}


# the vechu operator transforms a *symmetric* matrix 
# into a vector by stacking the *columns* of the matrix one after the
# other, but eliminating all infradiagonal elements
lav_matrix_vechu <- function(S, diagonal = TRUE) {
    S[lav_matrix_vechu_idx(n = ncol(S), diagonal = diagonal)]
}


# the vechru operator transforms a *symmetric* matrix 
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all infradiagonal elements
#
# same as vech (but using upper-diagonal elements)
lav_matrix_vechru <- function(S, diagonal = TRUE) {
    S[lav_matrix_vechru_idx(n = ncol(S), diagonal = diagonal)]
}




# return the *vector* indices of the lower triangular elements of a 
# symmetric matrix of size 'n'
lav_matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
    # FIXME: is there a way to avoid creating ROW/COL matrices?
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if(diagonal) which(ROW >= COL) else which(ROW > COL)
}

# return the *row* indices of the lower triangular elements of a 
# symmetric matrix of size 'n'
lav_matrix_vech_row_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    if(diagonal ) {
        unlist(lapply(seq_len(n), seq.int, n))
    } else {
        1 + unlist(lapply(seq_len(n-1), seq.int, n-1))
    }
}

# return the *col* indices of the lower triangular elements of a 
# symmetric matrix of size 'n'
lav_matrix_vech_col_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    if(!diagonal) {
        n <- n - 1L
    }
    rep.int(seq_len(n), times = rev(seq_len(n)))
}


 

# return the *vector* indices of the lower triangular elements of a 
# symmetric matrix of size 'n' -- ROW-WISE
lav_matrix_vechr_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n),   n, n, byrow = TRUE)
    tmp <- matrix(seq_len(n*n), n, n, byrow = TRUE)
    if(diagonal) tmp[ROW <= COL] else tmp[ROW < COL]
}

# return the *vector* indices of the upper triangular elements of a 
# symmetric matrix of size 'n' -- COLUMN-WISE
lav_matrix_vechu_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if(diagonal) which(ROW <= COL) else which(ROW < COL)
}

# return the *vector* indices of the upper triangular elements of a 
# symmetric matrix of size 'n' -- ROW-WISE
#
# FIXME!! make this more efficient (without creating 3 n*n matrices!)
#
lav_matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n),   n, n, byrow = TRUE)
    tmp <- matrix(seq_len(n*n), n, n, byrow = TRUE)
    if(diagonal) tmp[ROW >= COL] else tmp[ROW > COL]
}


# vech.reverse and vechru.reverse (aka `upper2full')
#
# given the output of vech(S) --or vechru(S) which is identical-- 
# reconstruct S
lav_matrix_vech_reverse <- lav_matrix_vechru_reverse <- 
lav_matrix_upper2full <- 
function(x, diagonal = TRUE) {
    # guess dimensions
    if(diagonal) {
        p <- (sqrt(1 + 8*length(x))-1)/2
    } else {
        p <- (sqrt(1 + 8*length(x))+1)/2
    }

    S <- numeric(p * p)
    S[lav_matrix_vech_idx(  p, diagonal = diagonal)] <- x
    S[lav_matrix_vechru_idx(p, diagonal = diagonal)] <- x

    attr(S, "dim") <- c(p, p)   
    S
}


# vechr.reverse vechu.reversie (aka `lower2full)
#
# given the output of vechr(S) --or vechu(S) which is identical--
# reconstruct S
lav_matrix_vechr_reverse <- lav_matrix_vechu_reverse <- 
lav_matrix_lower2full <- function(x, diagonal = TRUE) {
    # guess dimensions
    if(diagonal) {
        p <- (sqrt(1 + 8*length(x))-1)/2
    } else {
        p <- (sqrt(1 + 8*length(x))+1)/2
    }
    stopifnot(p == round(p,0))

    S <- numeric(p * p)
    S[lav_matrix_vechr_idx(p, diagonal = diagonal)] <- x
    S[lav_matrix_vechu_idx(p, diagonal = diagonal)] <- x

    attr(S, "dim") <- c(p, p)
    S
}


# return the *vector* indices of the diagonal elements of a symmetric
# matrix of size 'n'
lav_matrix_diag_idx <- function(n = 1L) {
    if(n < 1L) return(integer(0L))
    1L + 0L:(n-1L)*(n+1L)
}


# return the *vector* indices of the diagonal elements of the LOWER part
# of a symmatrix matrix of size 'n'
lav_matrix_diagh_idx <- function(n = 1L) {
    if(n < 1L) return(integer(0L))
    if(n == 1L) return(1L)
    c(1L, cumsum(n:2L) + 1L)
}


# return the *vector* indices of the ANTI diagonal elements of a symmetric
# matrix of size 'n'
lav_matrix_antidiag_idx <- function(n = 1L) {
    if(n < 1L) return(integer(0L))
    1L + seq_len(n)*(n-1L)
}



# create the duplication matrix (D_n): it 'duplicates' the elements
# in vech(S) to create vec(S) (where S is symmetric)
#
# D %*% vech(S) == vec(S)
#
# M&N book: pages 48-50
#
# note: several flavors: dup1, dup2, dup3, ...
# currently used: dup3

# first attempt
# dup1: working on the vector indices only
.dup1 <- function(n = 1L) {

    if ((n < 1L) | (round(n) != n)) {
        stop("n must be a positive integer")
    }
    if (n > 255L) {
        stop("n is too large")
    }

    # dimensions
    n2 <- n * n; nstar <- n * (n + 1)/2
    # THIS is the real bottleneck: allocating an ocean of zeroes...
    x <- numeric(n2 * nstar)

    # delta patterns
    r1 <- seq.int(from = n*n+1, by = -(n-1), length.out = n-1)
    r2 <- seq.int(from = n-1,   by = n-1,    length.out = n-1)
    r3 <- seq.int(from = 2*n+1, by = n,      length.out = n-1)

    # is there a more elegant way to do this?
    rr <- unlist(lapply((n-1):1, 
                function(x) { c(rbind(r1[1:x], r2[1:x]), r3[n-x]) }))

    idx <- c(1L, cumsum(rr) + 1L)

    # create matrix
    x[idx] <- 1.0
    attr(x, "dim") <- c(n2, nstar)

    x
}

# second attempt
# dup2: working on the row/col matrix indices
# (but only create the matrix at the very end)
.dup2 <- function (n = 1L) {

    if ((n < 1L) | (round(n) != n)) {
        stop("n must be a positive integer")
    }

    if(n > 255L) {
        stop("n is too large")
    }
    
    nstar <- n * (n+1)/2
    n2    <- n * n
    # THIS is the real bottleneck: allocating an ocean of zeroes...
    x <- numeric(n2 * nstar)

    idx1 <- lav_matrix_vech_idx(n)  + ((1L:nstar)-1L) * n2 # vector indices
    idx2 <- lav_matrix_vechru_idx(n) + ((1L:nstar)-1L) * n2 # vector indices

    x[idx1] <- 1.0
    x[idx2] <- 1.0

    attr(x, "dim") <- c(n2, nstar)
    x
}

# dup3: using col idx only
# D7 <- dup(7L); x<- apply(D7, 1, function(x) which(x > 0)); matrix(x,7,7)
.dup3 <- function(n = 1L) {
    if ((n < 1L) | (round(n) != n)) {
        stop("n must be a positive integer")
    }

    if(n > 255L) {
        stop("n is too large")
    }

    nstar <- n * (n+1)/2
    n2    <- n * n
    # THIS is the real bottleneck: allocating an ocean of zeroes...
    x <- numeric(n2 * nstar)

    tmp <- matrix(0L, n, n)
    tmp[lav_matrix_vech_idx(n)] <- 1:nstar
    tmp[lav_matrix_vechru_idx(n)] <- 1:nstar

    idx <- (1:n2) + (lav_matrix_vec(tmp)-1L) * n2

    x[idx] <- 1.0   

    attr(x, "dim") <- c(n2, nstar)
    x
}


# dup4: using Matrix package, returning a sparse matrix
#.dup4 <- function(n = 1L) {
#    if ((n < 1L) | (round(n) != n)) {
#        stop("n must be a positive integer")
#    }
#
#    if(n > 255L) {
#        stop("n is too large")
#    }
#
#    nstar <- n * (n+1)/2
#    #n2    <- n * n
#
#    tmp <- matrix(0L, n, n)
#    tmp[lav_matrix_vech_idx(n)] <- 1:nstar
#    tmp[lav_matrix_vechru_idx(n)] <- 1:nstar
#
#    x <- Matrix::sparseMatrix(i = 1:(n*n), j = vec(tmp), x = 1.0)
#
#    x
#}

# default dup:
lav_matrix_duplication <- .dup3


# compute t(D) %*% A (without explicitly computing D)
# sqrt(nrow(A)) is an integer
# A is not symmetric, and not even square, only n^2 ROWS
lav_matrix_duplication_pre <- function(A = matrix(0,0,0)) {

    # number of rows
    n2 <- nrow(A)

    # square nrow(A) only, n2 = n^2
    stopifnot(sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- lav_matrix_vech_idx(n); idx2 <- lav_matrix_vechru_idx(n)

    OUT <- A[idx1,,drop=FALSE] + A[idx2,,drop=FALSE]
    u <- which(idx1 %in% idx2); OUT[u,] <- OUT[u,] / 2.0
    
    OUT
}

# dupr_pre is faster...
lav_matrix_duplication_dup_pre2 <- function(A = matrix(0,0,0)) {

    # number of rows
    n2 <- nrow(A)

    # square nrow(A) only, n2 = n^2
    stopifnot(sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- lav_matrix_vech_idx(n); idx2 <- lav_matrix_vechru_idx(n)

    OUT <- A[idx1,,drop=FALSE]
    u <- which(!idx1 %in% idx2); OUT[u,] <- OUT[u,] + A[idx2[u],]

    OUT
}


# compute A %*% D (without explicitly computing D)
# sqrt(ncol(A)) must be an integer
# A is not symmetric, and not even square, only n^2 COLUMNS
lav_matrix_duplication_post <- function(A = matrix(0,0,0)) {

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- lav_matrix_vech_idx(n); idx2 <- lav_matrix_vechru_idx(n)

    OUT <- A[,idx1] + A[,idx2]
    u <- which(idx1 %in% idx2); OUT[,u] <- OUT[,u] / 2.0

    OUT
}

# compute t(D) %*% A %*% D (without explicitly computing D)
# A must be a square matrix and sqrt(ncol) an integer
lav_matrix_duplication_pre_post <- function(A = matrix(0,0,0)) {

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(nrow(A) == n2, sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- lav_matrix_vech_idx(n); idx2 <- lav_matrix_vechru_idx(n)

    OUT <- A[idx1,,drop=FALSE] + A[idx2,,drop=FALSE]
    u <- which(idx1 %in% idx2);     OUT[u,] <- OUT[u,] / 2.0
    OUT <- OUT[,idx1,drop=FALSE] + OUT[,idx2,drop=FALSE]
    OUT[,u] <- OUT[,u] / 2.0

    OUT
}

# create the generalized inverse of the duplication matrix (D^+_n): 
# it removes the duplicated elements in vec(S) to create vech(S)
#
# D^+ %*% vec(S) == vech(S)
#
# M&N book: page 49
#
# D^+ == solve(t(D_n %*% D_n) %*% t(D_n)

# create first t(DUP.ginv)
.dup_ginv1 <- function(n = 1L) {
    if ((n < 1L) | (round(n) != n)) {
        stop("n must be a positive integer")
    }

    if(n > 255L) {
        stop("n is too large")
    }

    nstar <- n * (n+1)/2
    n2    <- n * n
    # THIS is the real bottleneck: allocating an ocean of zeroes...
    x <- numeric(nstar * n2)

    tmp <- matrix(1:(n*n), n, n)
    idx1 <- lav_matrix_vech(tmp) + (0:(nstar-1L))*n2
    x[idx1] <- 0.5
    idx2 <- lav_matrix_vechru(tmp) + (0:(nstar-1L))*n2
    x[idx2] <- 0.5
    idx3 <- lav_matrix_diag_idx(n) + (lav_matrix_diagh_idx(n)-1L)*n2
    x[idx3] <- 1.0

    attr(x, "dim") <- c(n2, nstar)
    x <- t(x)
    x
}

# create DUP.ginv without transpose
.dup_ginv2 <- function(n = 1L) {
    if ((n < 1L) | (round(n) != n)) {
        stop("n must be a positive integer")
    }

    if(n > 255L) {
        stop("n is too large")
    }

    nstar <- n * (n+1)/2
    n2    <- n * n
    # THIS is the real bottleneck: allocating an ocean of zeroes...
    x <- numeric(nstar * n2)

    x[(lav_matrix_vech_idx(n)   - 1L)*nstar + 1:nstar] <- 0.5
    x[(lav_matrix_vechru_idx(n) - 1L)*nstar + 1:nstar] <- 0.5
    x[(lav_matrix_diag_idx(n) - 1L)*nstar + lav_matrix_diagh_idx(n)] <- 1.0

    attr(x, "dim") <- c(nstar, n2)
    x
}

lav_matrix_duplication_ginv <- .dup_ginv2

# pre-multiply with D^+
# number of rows in A must be 'square' (n*n)
lav_matrix_duplication_ginv_pre <- function(A = matrix(0,0,0)) {

    A <- as.matrix(A)

    # number of rows
    n2 <- nrow(A)

    # square nrow(A) only, n2 = n^2
    stopifnot(sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)
    nstar <- n * (n+1)/2

    idx1 <- lav_matrix_vech_idx(n); idx2 <- lav_matrix_vechru_idx(n)
    OUT <- (A[idx1,,drop=FALSE] + A[idx2,,drop=FALSE]) / 2
    OUT
}

# post-multiply with t(D^+)
# number of columns in A must be 'square' (n*n)
lav_matrix_duplication_ginv_post <- function(A = matrix(0,0,0)) {

    A <- as.matrix(A)

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    idx1 <- lav_matrix_vech_idx(n); idx2 <- lav_matrix_vechru_idx(n)
    OUT <- (A[,idx1,drop=FALSE] + A[,idx2,drop=FALSE]) / 2
    OUT
}

# pre AND post-multiply with D^+: D^+ %*% A %*% t(D^+)
# for square matrices only, with ncol = nrow = n^2
lav_matrix_duplication_ginv_pre_post <- function(A = matrix(0,0,0)) {

    A <- as.matrix(A)

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(nrow(A) == n2, sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)
   
    idx1 <- lav_matrix_vech_idx(n); idx2 <- lav_matrix_vechru_idx(n)
    OUT <- (A[idx1,,drop=FALSE] + A[idx2,,drop=FALSE]) / 2
    OUT <- (OUT[,idx1,drop=FALSE] + OUT[,idx2,drop=FALSE]) / 2
    OUT
}




# create the commutation matrix (K_mn) 
# the mn x mx commutation matrix is a permutation matrix which
# transforms vec(A) into vec(A')
#
# K_mn %*% vec(A) == vec(A')
#
# (in Henderson & Searle 1979, it is called the vec-permutation matrix)
# M&N book: pages 46-48
#
# note: K_mn is a permutation matrix, so it is orthogonal: t(K_mn) = K_mn^-1
#       K_nm %*% K_mn == I_mn
# 
# it is called the 'commutation' matrix because it enables us to interchange
# ('commute') the two matrices of a Kronecker product, eg
#   K_pm (A %x% B) K_nq == (B %x% A)
#
# important property: it allows us to transform a vec of a Kronecker product
# into the Kronecker product of the vecs (if A is m x n and B is p x q):
#   vec(A %x% B) == (I_n %x% K_qm %x% I_p)(vec A %x% vec B)

# first attempt
.com1 <- function(m = 1L, n = 1L) {

    if ((m < 1L) | (round(m) != m)) {
        stop("n must be a positive integer")
    }

    if ((n < 1L) | (round(n) != n)) {
        stop("n must be a positive integer")
    }
    
    p <- m*n
    x <- numeric( p*p )

    pattern <- rep(c(rep((m+1L)*n, (m-1L)), n+1L), n)
    idx <- c(1L, 1L + cumsum(pattern)[-p])

    x[idx] <- 1.0
    attr(x, "dim") <- c(p,p)

    x
}

lav_matrix_commutation <- .com1

# compute K_n %*% A without explicitly computing K
# K_n = K_nn, so sqrt(nrow(A)) must be an integer!
# = permuting the rows of A
lav_matrix_commutation_pre <- function(A = matrix(0,0,0)) {

    # number of rows of A
    n2 <- nrow(A)

    # K_nn only (n2 = m * n)
    stopifnot(sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # compute row indices
    #row.idx <- as.integer(t(matrix(1:n2, n, n)))
    row.idx <- rep(1:n, each=n) + (0:(n-1L))*n

    OUT <- A[row.idx,,drop=FALSE]
    OUT   
}

# compute K_mn %*% A without explicitly computing K
# = permuting the rows of A
lav_matrix_commutation_mn_pre <- function(A, m = 1L, n = 1L) {

    # number of rows of A
    mn <- nrow(A)
    stopifnot(mn == m * n)

    # compute row indices
    # row.idx <- as.integer(t(matrix(1:mn, m, n)))
    row.idx <- rep(1:m, each=n) + (0:(n-1L))*m

    OUT <- A[row.idx,,drop=FALSE]
    OUT
}

# N_n == 1/2 (I_n^2 + K_nn)
# see MN page 48
#
# N_n == D_n %*% D^+_n
#
lav_matrix_commutation_Nn <- function(n = 1L) {
     stop("not implemented yet")
}

# (simplified) kronecker product for square matrices
lav_matrix_kronecker_square <- function(A, check = TRUE) {

    dimA <- dim(A); n <- dimA[1L]; n2 <- n*n
    if(check) {
        stopifnot(dimA[2L] == n)
    }

    # all possible combinations
    out <- tcrossprod(as.vector(A))

    # break up in n*n pieces, and rearrange
    dim(out) <- c(n,n,n,n)
    out <- aperm(out, perm = c(3,1,4,2))
    
    # reshape again, to form n2 x n2 matrix
    dim(out) <- c(n2, n2)

    out
}

# (simplified) faster kronecker product for symmetric matrices
# note: not faster, but the logic extends to vech versions
lav_matrix_kronecker_symmetric <- function(S, check = TRUE) {
    
    dimS <- dim(S); n <- dimS[1L]; n2 <- n*n
    if(check) {
        stopifnot(dimS[2L] == n)
    }
    
    # all possible combinations
    out <- tcrossprod(as.vector(S))
    
    # break up in n*(n*n) pieces, and rearrange
    dim(out) <- c(n,n*n,n) 
    out <- aperm(out, perm = c(3L,2L,1L))
    
    # reshape again, to form n2 x n2 matrix
    dim(out) <- c(n2, n2)
    
    out
}

# shortcut for the idiom 't(S2) %*% (S %x% S) %*% S2'
# where S is symmetric, and the rows of S2 correspond to
# the elements of S
# eg - S2 = DELTA (the jacobian dS/dtheta)
lav_matrix_tS2_SxS_S2 <- function(S2, S, check = TRUE) {

    # size of S
    n <- nrow(S)

    if(check) {
        stopifnot(nrow(S2) == n*n)
    }

    A <- matrix(S %*% matrix(S2, n, ), n*n,)
    A2 <- A[rep(1:n, each=n) + (0:(n-1L))*n,,drop = FALSE]
    crossprod(A, A2)
}

# shortcut for the idiom 't(D) %*% (S %x% S) %*% S'
# where S is symmetric, and D is the duplication matrix
lav_matrix_tD_SxS_D <- function(S) {
}

# square root of a positive definite symmetric matrix
lav_matrix_symmetric_sqrt <- function(S = matrix(0,0,0)) {

    n <- nrow(S)

    # eigen decomposition
    S.eigen <- eigen(S)
    V <- S.eigen$vectors; d <- S.eigen$values

    # 'fix' slightly negative tiny numbers
    d[d < 0] <- 0.0

    # sqrt the eigenvalues and reconstruct
    S.sqrt <- V %*% diag(sqrt(d), n, n) %*% t(V)

    S.sqrt
}

# orthogonal complement of a matrix A
# see Satorra (1992). Sociological Methodology, 22, 249-278, footnote 3:
#
# To compute such an orthogonal matrix, consider the p* x p* matrix P = I -
# A(A'A)^-1A', which is idempotent of rank p* - q. Consider the singular value
# decomposition P = HVH', where H is a p* x (p* - q) matrix of full column rank,
# and V is a (p* - q) x (p* - q) diagonal matrix. It is obvious that H'A = 0;
# hence, H is the desired orthogonal complement. This method of constructing an
# orthogonal complement was proposed by Heinz Neudecker (1990, pers. comm.). 
#
# update YR 21 okt 2014:
# - note that A %*% solve(t(A) %*% A) %*% t(A) == tcrossprod(qr.Q(qr(A)))
# - if we are using qr, we can as well use qr.Q to get the complement
#
lav_matrix_orthogonal_complement <- function(A = matrix(0,0,0)) {

     QR <- qr(A)
     ranK <- QR$rank

     # following Heinz Neudecker:
     #n <- nrow(A)
     #P <- diag(n) - tcrossprod(qr.Q(QR))
     #OUT <- svd(P)$u[, seq_len(n - ranK), drop = FALSE]

     Q <- qr.Q(QR, complete = TRUE)
     # get rid of the first ranK columns
     OUT <- Q[, -seq_len(ranK), drop = FALSE]

     OUT
}

# construct block diagonal matrix
lav_matrix_bdiag <- function(...) {

    if(nargs() == 0L) return(matrix(0,0,0))
    dots <- list(...)

    # create list of matrices
    if(is.list(dots[[1]])) {
        mlist <- dots[[1]]
    } else {
        mlist <- dots
    }
    if(length(mlist) == 1L) return(mlist[[1]])

    # more than 1 matrix
    nmat <- length(mlist)
    nrows <- sapply(mlist, nrow); crows <- cumsum(nrows)
    ncols <- sapply(mlist, ncol); ccols <- cumsum(ncols)
    trows <- sum(nrows)
    tcols <- sum(ncols)
    x <- numeric(trows * tcols)
   
    for(m in seq_len(nmat)) {
        if(m > 1L) {
            rcoffset <- trows*ccols[m-1] + crows[m-1]
        } else {
            rcoffset <- 0L
        }
        m.idx <- ( rep((0:(ncols[m] - 1L))*trows, each=nrows[m]) +
                   rep(1:nrows[m], ncols[m]) + rcoffset )
        x[m.idx] <- mlist[[m]]        
    }

    attr(x, "dim") <- c(trows, tcols)
    x
}

# crossproduct, but handling NAs pairwise
lav_matrix_crossprod <- function(A, B) {
    if(missing(B)) {
        B <- A
    }
    apply(A, 2L, function(x) colSums(B * x, na.rm=TRUE))
}



