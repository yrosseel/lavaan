


# return Gauss-Hermite quadrature rule for given order (n)
# return list: x = nodes, w = quadrature weights
#

# As noted by Wilf (1962, chapter 2, ex 9), the nodes are given by
# the eigenvalues of the Jacobi matrix; weights are given by the squares of the
# first components of the (normalized) eigenvectors, multiplied by sqrt(pi)
#
# (This is NOT identical to Golub & Welsch, 1968: as they used a specific
#  method tailored for tridiagonal symmetric matrices)
#
# TODO: look at internet site
#  github.com/ajt60gaibb/FastGaussQuadrature.jl/blob/master/src/gausshermite.jl
# featuring the work of Ignace Bogaert (UGent)
#
# approximation of the integral of 'f(x) * exp(-x*x)' from -inf to +inf
# by sum( f(x_i) * w_i )
#
# CHECK: sum(w_i) should be always sqrt(pi) = 1.772454
lav_integration_gauss_hermite_xw <- function(n = 21L, revert = FALSE) {  # nolint
  # force n to be an integer
  n <- as.integer(n)
  stopifnot(n > 0L)

  if (n == 1L) {
    x <- 0
    w <- sqrt(pi)
  } else {
    # construct symmetric, tridiagonal Jacobi matrix
    # diagonal = 0, -1/+1 diagonal is sqrt(1:(n-1)/2)
    u <- sqrt(seq.int(n - 1L) / 2) # upper diagonal of J
    jn <- matrix(0, n, n)
    didx <- lav_matrix_diag_idx(n)
    jn[(didx + 1)[-n]] <- u
    # Jn[(didx-1)[-1]] <- u # only lower matrix is used anyway

    # eigen decomposition
    # FIXME: use specialized function for tridiagonal symmetric matrix
    ev <- eigen(jn, symmetric = TRUE)
    x <- ev$values
    tmp <- ev$vectors[1L, ]
    w <- sqrt(pi) * tmp * tmp
  }

  # revert? (minus to plus)
  if (revert) {
    x <- -x
  }

  list(x = x, w = w)
}

# generate GH points + weights
lav_integration_gauss_hermite <- function(n = 21L,
                                          dnorm = FALSE,
                                          mean = 0, sd = 1,
                                          ndim = 1L,
                                          revert = TRUE,
                                          prune = FALSE) {
  xw <- lav_integration_gauss_hermite_xw(n = n, revert = revert)

  # dnorm kernel?
  if (dnorm) {
    # scale/shift x
    x <- xw$x * sqrt(2) * sd + mean

    # scale w
    w <- xw$w / sqrt(pi)
  } else {
    x <- xw$x
    w <- xw$w
  }

  if (ndim > 1L) {
    # cartesian product
    x <- as.matrix(expand.grid(rep(list(x), ndim), KEEP.OUT.ATTRS = FALSE))
    w <- as.matrix(expand.grid(rep(list(w), ndim), KEEP.OUT.ATTRS = FALSE))
    w <- apply(w, 1, prod)
  } else {
    x <- as.matrix(x)
    w <- as.matrix(w)
  }

  # prune?
  if (is.logical(prune) && prune) {
    # always divide by N=21
    lower_limit <- xw$w[1] * xw$w[floor((n + 1) / 2)] / 21
    keep_idx <- which(w > lower_limit)
    w <- w[keep_idx]
    x <- x[keep_idx, , drop = FALSE]
  } else if (is.numeric(prune) && prune > 0) {
    lower_limit <- quantile(w, probs = prune)
    keep_idx <- which(w > lower_limit)
    w <- w[keep_idx]
    x <- x[keep_idx, , drop = FALSE]
  }

  list(x = x, w = w)
}
