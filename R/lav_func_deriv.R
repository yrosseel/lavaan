# numerical derivatives using complex numbers
# see Squire & Trapp 1998, siam rev 40(1) 110-112
# or Ridout, MS (2009), the american statistician 63(1) 66-74

# it would seem that you can choose h to be fairly small, without
# sacrifycing accuracy due to rounding errors

# YR 17 July 2012

lav_func_gradient_complex <- function(func, x,
                                      h = .Machine$double.eps, ...,
                                      fallback.simple = TRUE) {
  f0 <- try(func(x * (0 + 1i), ...), silent = TRUE)
  if (!is.complex(f0)) {
    if (fallback.simple) {
      dx <- lav_func_gradient_simple(func = func, x = x, h = sqrt(h), ...)
      return(dx)
    } else {
      lav_msg_stop(gettext(
        "function does not return a complex value")) # eg abs()
    }
  }
  if (inherits(f0, "try-error")) {
    if (fallback.simple) {
      dx <- lav_func_gradient_simple(func = func, x = x, h = sqrt(h), ...)
      return(dx)
    } else {
      lav_msg_stop(gettext(
        "function does not support non-numeric (complex) argument"))
    }
  }
  if (length(f0) != 1L) {
    lav_msg_stop(gettext(
      "function is not scalar and returns more than one element"))
  }

  nvar <- length(x)

  # determine 'h' per element of x
  h <- pmax(h, abs(h * x))

  # get exact h, per x
  tmp <- x + h
  h <- (tmp - x)

  # simple 'forward' method
  dx <- rep(as.numeric(NA), nvar)
  for (p in seq_len(nvar)) {
    dx[p] <- Im(func(x + h * 1i * (seq.int(nvar) == p), ...)) / h[p]
  }

  dx
}

# as a backup, if func() is not happy about non-numeric arguments
lav_func_gradient_simple <- function(func, x,
                                     h = sqrt(.Machine$double.eps), ...) {
  # check current point, see if it is a scalar function
  f0 <- func(x, ...)
  if (length(f0) != 1L) {
    lav_msg_stop(gettext(
      "function is not scalar and returns more than one element"))
  }

  nvar <- length(x)

  # determine 'h' per element of x
  h <- pmax(h, abs(h * x))

  # get exact h, per x
  tmp <- x + h
  h <- (tmp - x)

  # simple 'forward' method
  dx <- rep(as.numeric(NA), nvar)
  for (p in seq_len(nvar)) {
    dx[p] <- (func(x + h * (seq.int(nvar) == p), ...) - func(x, ...)) / h[p]
  }

  dx
}

lav_func_jacobian_complex <- function(func, x,
                                      h = .Machine$double.eps, ...,
                                      fallback.simple = TRUE) {
  f0 <- try(func(x * (0 + 1i), ...), silent = TRUE)
  if (!is.complex(f0)) {
    if (fallback.simple) {
      dx <- lav_func_jacobian_simple(func = func, x = x, h = sqrt(h), ...)
      return(dx)
    } else {
      lav_msg_stop(gettext(
        "function does not return a complex value")) # eg abs()
    }
  }
  if (inherits(f0, "try-error")) {
    if (fallback.simple) {
      dx <- lav_func_jacobian_simple(func = func, x = x, h = sqrt(h), ...)
      return(dx)
    } else {
      lav_msg_stop(gettext(
        "function does not support non-numeric (complex) argument"))
    }
  }
  nres <- length(f0)
  nvar <- length(x)

  # determine 'h' per element of x
  h <- pmax(h, abs(h * x))

  # get exact h, per x
  tmp <- x + h
  h <- (tmp - x)

  # simple 'forward' method
  dx <- matrix(as.numeric(NA), nres, nvar)
  for (p in seq_len(nvar)) {
    dx[, p] <- Im(func(x + h * 1i * (seq.int(nvar) == p), ...)) / h[p]
  }

  dx
}

lav_func_jacobian_simple <- function(func, x,
                                     h = sqrt(.Machine$double.eps), ...) {
  f0 <- func(x, ...)
  nres <- length(f0)
  nvar <- length(x)

  # determine 'h' per element of x
  h <- pmax(h, abs(h * x))

  # get exact h, per x
  tmp <- x + h
  h <- (tmp - x)

  # simple 'forward' method
  dx <- matrix(as.numeric(NA), nres, nvar)
  for (p in seq_len(nvar)) {
    dx[, p] <- (func(x + h * (seq.int(nvar) == p), ...) - func(x, ...)) / h[p]
  }

  dx
}

# this is based on the Ridout (2009) paper, and the code snippet for 'h4'
lav_func_hessian_complex <- function(func, x,
                                     h = .Machine$double.eps, ...) {
  f0 <- try(func(x * (0 + 1i), ...), silent = TRUE)
  if (!is.complex(f0)) {
    lav_msg_stop(gettext(
      "function does not return a complex value")) # eg abs()
  }
  if (inherits(f0, "try-error")) {
    lav_msg_stop(gettext(
      "function does not support non-numeric (complex) argument"))
  }
  if (length(f0) != 1L) {
    lav_msg_stop(gettext(
      "function is not scalar and returns more than one element"))
  }

  nvar <- length(x)

  # determine 'h' per element of x
  # delta1 <- pmax(h^(1/3), abs(h^(1/3)*x))
  # delta2 <- pmax(h^(1/5), abs(h^(1/5)*x))
  delta1 <- h^(1 / 3)
  delta2 <- h^(1 / 5)

  H <- matrix(as.numeric(NA), nvar, nvar)
  for (i in seq_len(nvar)) {
    for (j in 1:i) {
      if (i == j) {
        delta <- delta2
      } else {
        delta <- delta1
      }
      H[i, j] <- H[j, i] <-
        Im(func(x + delta * 1i * (seq.int(nvar) == i) * x +
          delta * (seq.int(nvar) == j) * x, ...) -
          func(x + delta * 1i * (seq.int(nvar) == i) * x -
            delta * (seq.int(nvar) == j) * x, ...)) /
          (2 * delta * delta * x[i] * x[j])
    }
  }

  H
}

lav_deriv_cov2corB <- function(COV = NULL) {
  nvar <- nrow(COV)
  dS.inv <- 1 / diag(COV)
  R <- cov2cor(COV)
  A <- -R %x% (0.5 * diag(dS.inv))
  B <- (0.5 * diag(dS.inv)) %x% -R
  DD <- diag(lav_matrix_vec(diag(nvar)))
  A2 <- A %*% DD
  B2 <- B %*% DD
  out <- A2 + B2 + diag(lav_matrix_vec(tcrossprod(sqrt(dS.inv))))
  D <- lav_matrix_duplication(nvar)
  out.vech <- 0.5 * (t(D) %*% out %*% D)
  out.vech
}

# quick and dirty (FIXME!!!) way to get
# surely there must be a more elegant way?
# see lav_deriv_cov2corB, if no num.idx...
# dCor/dCov
lav_deriv_cov2cor <- function(COV = NULL, num.idx = NULL) {
  # dCor/dvar1 = - cov / (2*var1 * sqrt(var1) * sqrt(var2))
  # dCor/dvar2 = - cov / (2*var2 * sqrt(var1) * sqrt(var2))
  # dCor/dcov  =  1/(sqrt(var1) * sqrt(var2))

  # diagonal: diag(lav_matrix_vech(tcrossprod(1/delta)))

  nvar <- ncol(COV)
  pstar <- nvar * (nvar + 1) / 2
  delta <- sqrt(diag(COV))
  if (length(num.idx) > 0L) {
    delta[num.idx] <- 1.0
  }

  A <- COV * -1 / (2 * delta * delta * tcrossprod(delta))
  if (length(num.idx) > 0L) {
    A[num.idx, ] <- 0
    A[cbind(num.idx, num.idx)] <- 1
  }
  A2 <- diag(nvar) %x% t(A)

  OUT <- diag(pstar)
  diag(OUT) <- lav_matrix_vech(tcrossprod(1 / delta))
  var.idx <- lav_matrix_diagh_idx(nvar)
  DUP <- lav_matrix_duplication(nvar)
  OUT[, var.idx] <- t(DUP) %*% A2[, lav_matrix_diag_idx(nvar)]

  if (length(num.idx) > 0L) {
    var.idx <- var.idx[-num.idx]
  }
  OUT[var.idx, var.idx] <- 0

  OUT
}


lav_deriv_cov2cor_numerical <- function(COV, num.idx = integer(0)) {
  compute.R <- function(x) {
    S <- lav_matrix_vech_reverse(x)
    diagS <- diag(S)
    delta <- 1 / sqrt(diagS)
    if (length(num.idx) > 0L) {
      delta[num.idx] <- 1.0
    }
    R <- diag(delta) %*% S %*% diag(delta)
    # R <- cov2cor(S)
    R.vec <- lav_matrix_vech(R, diagonal = TRUE)
    R.vec
  }

  x <- lav_matrix_vech(COV, diagonal = TRUE)
  dx <- lav_func_jacobian_complex(func = compute.R, x = x)

  dx
}
