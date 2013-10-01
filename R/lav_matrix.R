# Magnus & Neudecker (1999) style matrix operations
# initial version: YR - 11 may 2011


# vec operator
#
# the vec operator (for 'vectorization') transforms a matrix into 
# a vector by stacking the *columns* of the matrix one underneath the other
#
# M&N book: page 30
#
# note: we do not coerce to 'double/numeric' storage-mode (like as.numeric)
vec <- function(A) {
    as.vector(A)
}

# vecr operator
#
# the vecr operator ransforms a matrix into 
# a vector by stacking the *rows* of the matrix one underneath the other
vecr <- function(A) {
    vec(t(A))
}

# vech 
# 
# the vech operator (for 'half vectorization') transforms a *symmetric* matrix 
# into a vector by stacking the *columns* of the matrix one underneath the 
# other, but eliminating all supradiagonal elements
#
# M&N book: page 48-49
#
vech <- function(S, diagonal = TRUE) {
    ROW <- row(S); COL <- col(S)
    if(diagonal) S[ROW >= COL] else S[ROW > COL]
}

# the vechru operator transforms a *symmetric* matrix 
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all infradiagonal elements
#
# same as vech (but using upper-diagonal elements)
vechru <- function(S, diagonal = TRUE) {
    S[vechru.idx(n = ncol(S), diagonal = diagonal)]
}

# return the the *vector* indices of the lower triangular elements of a 
# symmetric matrix of size 'n'
vech.idx <- function(n = 1L, diagonal = TRUE) {
    # FIXME: is there a way to avoid creating ROW/COL matrices?
    n <- as.integer(n)
    ROW <- matrix(1:n, n, n); COL <- matrix(1:n, n, n, byrow=TRUE)
    if(diagonal) which(ROW >= COL) else which(ROW > COL)
}
 

# return the the *vector* indices of the upper triangular elements of a 
# symmetric matrix of size 'n' -- ROW-WISE
#
# FIXME!! make this more efficient (without creating 3 n*n matrices!)
#
vechru.idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(1:n, n, n); COL <- matrix(1:n, n, n, byrow=TRUE)
    tmp <- matrix(1:(n*n), n, n, byrow=TRUE)
    if(diagonal) tmp[ROW >= COL] else tmp[ROW > COL]
}


# vech.reverse and vechru.reverse (aka `upper2full')
#
# given the output of vech(S) --or vechru(S) which is identical-- 
# reconstruct S
vech.reverse <- vechru.reverse <- upper2full <- function(x, diagonal = TRUE) {
    # guess dimensions
    if(diagonal) {
        p <- (sqrt(1 + 8*length(x))-1)/2
    } else {
        p <- (sqrt(1 + 8*length(x))+1)/2
    }

    S <- numeric(p*p)
    S[vech.idx( p, diagonal=diagonal)] <- x
    S[vechru.idx(p, diagonal=diagonal)] <- x

    attr(S, "dim") <- c(p, p)   
    S
}



# the vechr operator transforms a *symmetric* matrix 
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all supradiagonal elements
vechr <- function(S, diagonal = TRUE) {
    S[vechr.idx(n = ncol(S), diagonal = diagonal)]
}

# the vechu operator transforms a *symmetric* matrix 
# into a vector by stacking the *columns* of the matrix one after the
# other, but eliminating all infradiagonal elements
vechu <- function(S, diagonal = TRUE) {
    S[vechu.idx(n = ncol(S), diagonal = diagonal)]
}

# return the the *vector* indices of the lower triangular elements of a 
# symmetric matrix of size 'n' -- ROW-WISE
vechr.idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(1:n, n, n); COL <- matrix(1:n, n, n, byrow=TRUE)
    tmp <- matrix(1:(n*n), n, n, byrow=TRUE)
    if(diagonal) tmp[ROW <= COL] else tmp[ROW < COL]
}

# return the the *vector* indices of the upper triangular elements of a 
# symmetric matrix of size 'n' -- COLUMN-WISE
vechu.idx <- function(n = 1L, diagonal = TRUE) {
    n <- as.integer(n)
    ROW <- matrix(1:n, n, n); COL <- matrix(1:n, n, n, byrow=TRUE)
    if(diagonal) which(ROW <= COL) else which(ROW < COL)
}

# vechr.reverse vechu.reversie (aka `lower2full)
#
# given the output of vechr(S) --or vechu(S) which is identical--
# reconstruct S
vechr.reverse <- vechu.reverse <- lower2full <- function(x, diagonal = TRUE) {
    # guess dimensions
    if(diagonal) {
        p <- (sqrt(1 + 8*length(x))-1)/2
    } else {
        p <- (sqrt(1 + 8*length(x))+1)/2
    }
    stopifnot(p == round(p,0))

    S <- numeric(p*p)
    S[vechr.idx(p, diagonal=diagonal)] <- x
    S[vechu.idx(p, diagonal=diagonal)] <- x

    attr(S, "dim") <- c(p, p)
    S
}

# return the *vector* indices of the diagonal elements of a symmetric
# matrix of size 'n'
diag.idx <- function(n = 1L) {
    if(n < 1L) return(integer(0L))
    1L + 0L:(n-1L)*(n+1L)
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
dup1 <- function(n = 1L) {

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
dup2 <- function (n = 1L) {

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

    idx1 <- vech.idx(n)  + ((1L:nstar)-1L) * n2 # vector indices
    idx2 <- vechru.idx(n) + ((1L:nstar)-1L) * n2 # vector indices

    x[idx1] <- 1.0
    x[idx2] <- 1.0

    attr(x, "dim") <- c(n2, nstar)
    x
}

# dup3: using col idx only
# D7 <- dup(7L); x<- apply(D7, 1, function(x) which(x > 0)); matrix(x,7,7)
dup3 <- function(n = 1L) {
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
    tmp[vech.idx(n)] <- 1:nstar
    tmp[vechru.idx(n)] <- 1:nstar

    idx <- (1:n2) + (vec(tmp)-1L) * n2

    x[idx] <- 1.0   

    attr(x, "dim") <- c(n2, nstar)
    x
}


# dup4: using Matrix package, returning a sparse matrix
#dup4 <- function(n = 1L) {
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
#    tmp[vech.idx(n)] <- 1:nstar
#    tmp[vechru.idx(n)] <- 1:nstar
#
#    x <- Matrix::sparseMatrix(i = 1:(n*n), j = vec(tmp), x = 1.0)
#
#    x
#}

# default dup:
duplicationMatrix <- dup3


# compute t(D) %*% A (without explicitly computing D)
# sqrt(ncol(A)) is an integer
# A is not symmetric, and not even square
D.pre <- function(A = matrix(0,0,0)) {

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(nrow(A) == n2, sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- vech.idx(n); idx2 <- vechru.idx(n)

    OUT <- A[idx1,,drop=FALSE] + A[idx2,,drop=FALSE]
    u <- which(idx1 %in% idx2); OUT[u,] <- OUT[u,] / 2.0
    
    OUT
}

# dup.pre is faster...
dup.pre2 <- function(A = matrix(0,0,0)) {

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(nrow(A) == n2, sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- vech.idx(n); idx2 <- vechru.idx(n)

    OUT <- A[idx1,,drop=FALSE]
    u <- which(!idx1 %in% idx2); OUT[u,] <- OUT[u,] + A[idx2[u],]

    OUT
}


# compute A %*% D (without explicitly computing D)
# sqrt(ncol(A)) must be an integer
D.post <- function(A = matrix(0,0,0)) {

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(ncol(A) == n2, sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- vech.idx(n); idx2 <- vechru.idx(n)

    OUT <- A[,idx1] + A[,idx2]
    u <- which(idx1 %in% idx2); OUT[,u] <- OUT[,u] / 2.0

    OUT
}

# compute t(D) %*% A %*% D (without explicitly computing D)
# A must be a square matrix and sqrt(ncol) an integer
D.pre.post <- function(A = matrix(0,0,0)) {

    # number of columns
    n2 <- ncol(A)

    # square A only, n2 = n^2
    stopifnot(nrow(A) == n2, sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # dup idx
    idx1 <- vech.idx(n); idx2 <- vechru.idx(n)

    OUT <- A[idx1,,drop=FALSE] + A[idx2,,drop=FALSE]
    u <- which(idx1 %in% idx2);     OUT[u,] <- OUT[u,] / 2.0
    OUT <- OUT[,idx1,drop=FALSE] + OUT[,idx2,drop=FALSE]
    OUT[,u] <- OUT[,u] / 2.0

    OUT
}

# create the commutation matrix (K_mn) 
# the mn x mx commutation matrix is a permutation matrix which
# transforms vec(A) into vec(A')
#
# K %*% vec(A) == vec(A')
#
# M&N book: pages 46-48
#
# note: several flavors: com1, com2, ...
# currently used: com1

# first attempt
com1 <- function(m = 1L, n = 1L) {

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

commutationMatrix <- com1

# compute K_n %*% A without explicitly computing K
# K_n = K_nn, so sqrt(nrow(A)) must be an integer!
# = permuting the rows of A
K.n.pre <- function(A) {

    # number of rows of A
    n2 <- nrow(A)

    # K_nn only (n2 = m * n)
    stopifnot(sqrt(n2) == round(sqrt(n2)))

    # dimension
    n <- sqrt(n2)

    # compute row indices
    row.idx <- as.integer(t(matrix(1:n2, n, n)))

    OUT <- A[row.idx,,drop=FALSE]
    OUT   
}

# compute K_mn %*% A without explicitly computing K
# = permuting the rows of A
K.mn.pre <- function(A, m = 1L, n = 1L) {

    # number of rows of A
    mn <- nrow(A)
    stopifnot(mn == m * n)

    # compute row indices
    row.idx <- as.integer(t(matrix(1:mn, m, n)))

    OUT <- A[row.idx,,drop=FALSE]
    OUT
}

# square root of a positive definite symmetric matrix
sqrtSymmetricMatrix <- function(S = matrix(0,0,0)) {

    n <- nrow(S)

    # eigen decomposition
    S.eigen <- eigen(S)
    V <- S.eigen$vectors; d <- S.eigen$values

    # sqrt the eigenvalues and reconstruct
    S.sqrt <- V %*% diag(sqrt(d), n, n) %*% t(V)

    S.sqrt
}

# orthogonal complement of a matrix A
# see Satorra (1992). Sociological Methodology, 22, 249-278, footnote 3:
#
# To compute such an orthogonal matrix, consider the p* x p* matrix P = I -
# A(A'A)^1A', which is idempotent of rank p* - q. Consider the singular value
# decomposition P = HVH', where H is a p* x (p* - q) matrix of full column rank,
# and V is a (p* - q) x (p* - q) diagonal matrix. It is obvious that H'A = 0;
# hence, H is the desired orthogonal complement. This method of constructing an
# orthogonal complement was proposed by Heinz Neudecker (1990, pers. comm.). 
orthogonalComplement <- function(A = matrix(0,0,0)) {
     n <- nrow(A)
     ranK <- qr(A)$rank

     P <- diag(n) - A %*% solve(crossprod(A)) %*% t(A)
     H <- svd(P)$u[,1:(n-ranK),drop=FALSE]

     H
}



