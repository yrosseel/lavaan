# routines for numerical intregration

# return Gauss-Hermite quadrature rules for given order (n)
# return list: x = nodes, w = quadrature weights
#

# This is the Wilf 1962 method: nodes are given by the eigenvalues
# of the Jacobi matrix; weights are given by the squares of the
# first components of the (normalized) eigenvectors
#
# (This is NOT Golub & Welsch, 1968: they use a specific method
# tailored for tridiagonal symmetric matrices)
#
# approximation of the integral of 'f(x) * exp(-x*x)' from -inf to +inf
# by sum( f(x_i) * w_i )
#
# CHECK: sum(w_i) should be always sqrt(pi) = 1.772454
#
# Example: define g <- function(x) exp(-x*x) * sin(x^2)
# brute force solution: integrate(g, -100, 100) = 0.5703706
# using gauss_hermite, define f <- function(x) sin(x^2)
# XW <- lav_gauss_hermite_xw(10) 
# sum( f(XW$x) * XW$w ) = 0.5703706
#
lav_gauss_hermite_xw <- function(n = 100L, revert = FALSE) {

    # construct symmetric, tridiagonal Jacobi matrix
    # diagonal = 0, -1/+1 diagonal is sqrt(1:(n-1)/2)
    u <-  sqrt(1:(n-1)/2) # upper diagonal of J
    Jn <- matrix(0, n, n); didx <- diag.idx(n)
    Jn[(didx+1)[-n]] <- u
    #Jn[(didx-1)[-1]] <- u # only lower matrix is used anyway

    # eigen decomposition
    # FIXME: use specialized function for tridiagonal symmetrix matrix
    ev <- eigen(Jn, symmetric = TRUE)
    x <- ev$values
    tmp <- ev$vector[1L,]
    w <- sqrt(pi)*tmp*tmp

    # revert? (minus to plus)
    if(revert) {
        x <- -x
    }
    list(x=x, w=w)
}

# 'dnorm' version: weighting function is not exp(-x*x) but
# w(x) = 1/(sqrt(2*pi)) * exp(-0.5 * x*x)
#
# Example: define g <- function(x) sin(x^2) * dnorm(x, 0.3, 0.4)
# brute force solution: integrate(g, -100, 100) = 0.2227559
# using gauss_hermite_dnorm, define f <- function(x) sin(x^2)
# XW <- lav_gauss_hermite_xw_dnorm(n=10, mean=0.3, sd=0.4)
# sum( f(XW$x) * XW$w ) = 0.2227559
#
lav_gauss_hermite_xw_dnorm <- function(n = 100L, revert = FALSE, 
                                       mean = 0, sd = 1, ndim = 1L) {
    XW <- lav_gauss_hermite_xw(n = n, revert = revert)

    # scale/shift x
    x <- XW$x * sqrt(2) * sd + mean

    # scale w
    w <- XW$w / sqrt(pi)

    if(ndim > 1L) {
        # cartesian product
        x <- as.matrix(expand.grid(rep(list(x), ndim), KEEP.OUT.ATTRS = FALSE))
        w <- as.matrix(expand.grid(rep(list(w), ndim), KEEP.OUT.ATTRS = FALSE))
        w <- apply(w, 1, prod)
    } else {
        x <- as.matrix(x)
        w <- as.matrix(w)
    }

    list(x=x, w=w)
}
