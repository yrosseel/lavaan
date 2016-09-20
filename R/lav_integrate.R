# routines for numerical intregration

# integrate (-infty to +infty) a product of univariate Gaussian densities
# with givens means (mus) and standard deviations (sds) (or variances, vars)
lav_integration_gaussian_product <- function(mus = NULL, sds = NULL, vars = NULL) {

    n <- length(mus)
    if(is.null(vars)) {
        vars <- sds^2
    }

    # variance product
    var.prod <- 1/sum(1/vars)
     
    # mean product
    mu.prod <- sum(mus/vars)*var.prod

    # normalization constant
    const <- 1/sqrt((2*pi)^(n-1)) * sqrt(var.prod) * sqrt(1/prod(vars)) * exp(-0.5 * (sum(mus^2/vars) - mu.prod^2/var.prod))

    const
}



# return Gauss-Hermite quadrature rules for given order (n)
# return list: x = nodes, w = quadrature weights
#

# As noted by Wilf (1962, chapter 2, ex 9), the nodes are given by 
# the eigenvalues of the Jacobi matrix; weights are given by the squares of the
# first components of the (normalized) eigenvectors, multiplied by sqrt(pi)
#
# (This is NOT Golub & Welsch, 1968: as they used a specific method
# tailored for tridiagonal symmetric matrices)
#
# TODO: look at https://github.com/ajt60gaibb/FastGaussQuadrature.jl/blob/master/src/gausshermite.jl
# featuring the work of Ignace Bogaert (UGent)
#
# approximation of the integral of 'f(x) * exp(-x*x)' from -inf to +inf
# by sum( f(x_i) * w_i )
#
# CHECK: sum(w_i) should be always sqrt(pi) = 1.772454
#
# Example: define g <- function(x) dnorm(x,0,0.2) * exp(-x*x)
# plot:
#   x <- seq(-5,5,0.01)
#   plot(x, g(x))
# brute force solution: integrate(g, -100, 100) = 0.9622504
# using gauss_hermite, define f <- function(x) dnorm(x,0,0.2)
# XW <- lavaan:::lav_integration_gauss_hermite(200)
# sum( f(XW$x) * XW$w ) = 0.9622504
#
# note that we need N=200 to get decent accuracy here, because we have
# a rather spiky function

lav_integration_gauss_hermite <- function(n = 21L, revert = FALSE) {

    # force n to be an integer
    n <- as.integer(n); stopifnot(n > 0L)

    if(n == 1L) {
        x <- 0
        w <- sqrt(pi)
    } else {
        # construct symmetric, tridiagonal Jacobi matrix
        # diagonal = 0, -1/+1 diagonal is sqrt(1:(n-1)/2)
        u <-  sqrt(seq.int(n-1L)/2) # upper diagonal of J
        Jn <- matrix(0, n, n); didx <- lav_matrix_diag_idx(n)
        Jn[(didx+1)[-n]] <- u
        #Jn[(didx-1)[-1]] <- u # only lower matrix is used anyway

        # eigen decomposition
        # FIXME: use specialized function for tridiagonal symmetrix matrix
        ev <- eigen(Jn, symmetric = TRUE)
        x <- ev$values
        tmp <- ev$vector[1L,]
        w <- sqrt(pi)*tmp*tmp
    }

    # revert? (minus to plus)
    if(revert) {
        x <- -x
    }

    list(x = x, w = w)
}

# 'dnorm' version: weighting function is not exp(-x^2) but
# w(x) = 1/(sqrt(2*pi)) * exp(-0.5 * x^2)
#
# Example: define g <- function(x) dnorm(x,0,0.2) * dnorm(x, 0.3, 0.4)
# brute force solution: integrate(g, -100, 100) = 0.712326
# using gauss_hermite_dnorm, define f <- function(x) dnorm(x,0,0.2)
# XW <- lav_integration_gauss_hermite_dnorm(n=100, mean=0.3, sd=0.4)
# sum( f(XW$x) * XW$w ) = 0.712326
lav_integration_gauss_hermite_dnorm <- function(n = 21L, mean = 0, sd = 1, 
                                                ndim = 1L,
                                                revert = FALSE,
                                                prune = 0) {
    XW <- lav_integration_gauss_hermite(n = n, revert = revert)

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

    # prune?
    if(prune > 0) {
        lower.limit <- quantile(w, probs = prune)
        keep.idx <- which(w > lower.limit)
        w <- w[keep.idx]
        x <- x[keep.idx,, drop = FALSE]
    }

    list(x=x, w=w)
}

# plot 2-dim
# out <- lavaan:::lav_integration_gauss_hermite_dnorm(n = 20, ndim = 2)
# plot(out$x, cex = -10/log(out$w), col = "darkgrey", pch=19)
