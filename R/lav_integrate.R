# routines for numerical intregration

# integrate (-infty to +infty) a product of univariate Gaussian densities
# with givens means (mus) and standard deviations (sds) (or variances, vars)
lav_integration_gaussian_product <- function(mus = NULL, sds = NULL, vars = NULL) {
  n <- length(mus)
  if (is.null(vars)) {
    vars <- sds^2
  }

  # variance product
  var.prod <- 1 / sum(1 / vars)

  # mean product
  mu.prod <- sum(mus / vars) * var.prod

  # normalization constant
  const <- 1 / sqrt((2 * pi)^(n - 1)) * sqrt(var.prod) * sqrt(1 / prod(vars)) * exp(-0.5 * (sum(mus^2 / vars) - mu.prod^2 / var.prod))

  const
}



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
# TODO: look at https://github.com/ajt60gaibb/FastGaussQuadrature.jl/blob/master/src/gausshermite.jl
# featuring the work of Ignace Bogaert (UGent)
#
# approximation of the integral of 'f(x) * exp(-x*x)' from -inf to +inf
# by sum( f(x_i) * w_i )
#
# CHECK: sum(w_i) should be always sqrt(pi) = 1.772454
lav_integration_gauss_hermite_xw <- function(n = 21L, revert = FALSE) {
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
    Jn <- matrix(0, n, n)
    didx <- lav_matrix_diag_idx(n)
    Jn[(didx + 1)[-n]] <- u
    # Jn[(didx-1)[-1]] <- u # only lower matrix is used anyway

    # eigen decomposition
    # FIXME: use specialized function for tridiagonal symmetrix matrix
    ev <- eigen(Jn, symmetric = TRUE)
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
  XW <- lav_integration_gauss_hermite_xw(n = n, revert = revert)

  # dnorm kernel?
  if (dnorm) {
    # scale/shift x
    x <- XW$x * sqrt(2) * sd + mean

    # scale w
    w <- XW$w / sqrt(pi)
  } else {
    x <- XW$x
    w <- XW$w
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
    lower.limit <- XW$w[1] * XW$w[floor((n + 1) / 2)] / 21
    keep.idx <- which(w > lower.limit)
    w <- w[keep.idx]
    x <- x[keep.idx, , drop = FALSE]
  } else if (is.numeric(prune) && prune > 0) {
    lower.limit <- quantile(w, probs = prune)
    keep.idx <- which(w > lower.limit)
    w <- w[keep.idx]
    x <- x[keep.idx, , drop = FALSE]
  }

  list(x = x, w = w)
}

# backwards compatibility
lav_integration_gauss_hermite_dnorm <- function(n = 21L, mean = 0, sd = 1,
                                                ndim = 1L,
                                                revert = TRUE,
                                                prune = FALSE) {
  lav_integration_gauss_hermite(
    n = n, dnorm = TRUE, mean = mean, sd = sd,
    ndim = ndim, revert = revert, prune = prune
  )
}

# plot 2-dim
# out <- lavaan:::lav_integration_gauss_hermite_dnorm(n = 20, ndim = 2)
# plot(out$x, cex = -10/log(out$w), col = "darkgrey", pch=19)

# integrand g(x) has the form g(x) = f(x) dnorm(x, m, s^2)
lav_integration_f_dnorm <- function(func = NULL, # often ly.prod
                                    dnorm.mean = 0, # dnorm mean
                                    dnorm.sd = 1, # dnorm sd
                                    XW = NULL, # GH points
                                    n = 21L, # number of nodes
                                    adaptive = FALSE, # adaptive?
                                    iterative = FALSE, # iterative?
                                    max.iter = 20L, # max iterations
                                    verbose = FALSE, # verbose?
                                    ...) { # optional args for 'f'

  # create GH rule
  if (is.null(XW)) {
    XW <- lav_integration_gauss_hermite_xw(n = n, revert = TRUE)
  }

  if (!adaptive) {
    w.star <- XW$w / sqrt(pi)
    x.star <- dnorm.sd * (sqrt(2) * XW$x) + dnorm.mean
    out <- sum(func(x.star, ...) * w.star)
  } else {
    # Naylor & Smith (1982, 1988)
    if (iterative) {
      mu.est <- 0
      sd.est <- 1

      for (i in 1:max.iter) {
        w.star <- sqrt(2) * sd.est * dnorm(sqrt(2) * sd.est * XW$x + mu.est, dnorm.mean, dnorm.sd) * exp(XW$x^2) * XW$w
        x.star <- sqrt(2) * sd.est * XW$x + mu.est
        LIK <- sum(func(x.star, ...) * w.star)

        # update mu
        mu.est <- sum(x.star * (func(x.star, ...) * w.star) / LIK)

        # update sd
        var.est <- sum(x.star^2 * (func(x.star, ...) * w.star) / LIK) - mu.est^2
        sd.est <- sqrt(var.est)

        if (verbose) {
          cat(
            "i = ", i, "LIK = ", LIK, "mu.est = ", mu.est,
            "sd.est = ", sd.est, "\n"
          )
        }
      }
      out <- LIK

      # Liu and Pierce (1994)
    } else {
      # integrand g(x) = func(x) * dnorm(x; m, s^2)
      log.g <- function(x, ...) {
        ## FIXME: should we take the log right away?
        log(func(x, ...) * dnorm(x, mean = dnorm.mean, sd = dnorm.sd))
      }
      # find mu hat and sd hat
      mu.est <- optimize(
        f = log.g, interval = c(-10, 10),
        maximum = TRUE, tol = .Machine$double.eps, ...
      )$maximum
      H <- as.numeric(numDeriv::hessian(func = log.g, x = mu.est, ...))
      sd.est <- sqrt(1 / -H)

      w.star <- sqrt(2) * sd.est * dnorm(sd.est * (sqrt(2) * XW$x) + mu.est, dnorm.mean, dnorm.sd) * exp(XW$x^2) * XW$w
      x.star <- sd.est * (sqrt(2) * XW$x) + mu.est

      out <- sum(func(x.star, ...) * w.star)
    }
  }

  out
}

# integrand g(z) has the form g(z) = f(sz+m) dnorm(z, 0, 1)
lav_integration_f_dnorm_z <- function(func = NULL, # often ly.prod
                                      f.mean = 0, # f mean
                                      f.sd = 1, # f sd
                                      XW = NULL, # GH points
                                      n = 21L, # number of nodes
                                      adaptive = FALSE, # adaptive?
                                      iterative = FALSE, # iterative?
                                      max.iter = 20L, # max iterations
                                      verbose = FALSE, # verbose?
                                      ...) { # optional args for 'f'

  # create GH rule
  if (is.null(XW)) {
    XW <- lav_integration_gauss_hermite_xw(n = n, revert = TRUE)
  }

  if (!adaptive) {
    w.star <- XW$w / sqrt(pi)
    x.star <- sqrt(2) * XW$x
    out <- sum(func(f.sd * x.star + f.mean, ...) * w.star)
  } else {
    # Naylor & Smith (1982, 1988)
    if (iterative) {
      mu.est <- 0
      sd.est <- 1

      for (i in 1:max.iter) {
        w.star <- sqrt(2) * sd.est * dnorm(sd.est * sqrt(2) * XW$x + mu.est, 0, 1) * exp(XW$x^2) * XW$w
        x.star <- sd.est * (sqrt(2) * XW$x) + mu.est
        LIK <- sum(func(f.sd * x.star + f.mean, ...) * w.star)

        # update mu
        mu.est <- sum(x.star * (func(f.sd * x.star + f.mean, ...) * w.star) / LIK)

        # update sd
        var.est <- sum(x.star^2 * (func(f.sd * x.star + f.mean, ...) * w.star) / LIK) - mu.est^2
        sd.est <- sqrt(var.est)

        if (verbose) {
          cat(
            "i = ", i, "LIK = ", LIK, "mu.est = ", mu.est,
            "sd.est = ", sd.est, "\n"
          )
        }
      }
      out <- LIK

      # Liu and Pierce (1994)
    } else {
      # integrand g(x) = func(x) * dnorm(x; m, s^2)
      log.gz <- function(x, ...) {
        ## FIXME: should we take the log right away?
        log(func(f.sd * x + f.mean, ...) * dnorm(x, mean = 0, sd = 1))
      }
      # find mu hat and sd hat
      mu.est <- optimize(
        f = log.gz, interval = c(-10, 10),
        maximum = TRUE, tol = .Machine$double.eps, ...
      )$maximum
      H <- as.numeric(numDeriv::hessian(func = log.gz, x = mu.est, ...))
      sd.est <- sqrt(1 / -H)

      w.star <- sqrt(2) * sd.est * dnorm(sd.est * (sqrt(2) * XW$x) + mu.est, 0, 1) * exp(XW$x^2) * XW$w
      x.star <- sd.est * (sqrt(2) * XW$x) + mu.est

      out <- sum(func(f.sd * x.star + f.mean, ...) * w.star)
    }
  }

  out
}
