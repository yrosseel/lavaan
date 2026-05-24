# numerical derivatives using complex numbers
# see Squire & Trapp 1998, siam rev 40(1) 110-112
# or Ridout, MS (2009), the american statistician 63(1) 66-74

# it would seem that you can choose h to be fairly small, without
# sacrificing accuracy due to rounding errors

# YR 17 July 2012

lav_func_gradient_complex <- function(func, x,                        # nolint start
                                      h = .Machine$double.eps, ...,
                                      fallback.simple = TRUE) {       #nolint end
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

lav_func_jacobian_complex <- function(func, x,                       # nolint start
                                      h = .Machine$double.eps, ...,
                                      fallback.simple = TRUE) {      # nolint end
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

# Hessian computation for a (possibly vector-valued) function func: R^k -> R^m
# Returns a list of length m, where the i-th element is the (k x k)
# symmetric Hessian matrix of the i-th component of func at x.
#
# The inner Jacobian is computed using the complex-step method (high
# accuracy), while the outer step uses forward finite differences.
# A small symmetrization step compensates for finite-difference asymmetry.
#
# Components whose Jacobian row is exactly zero at x (i.e. the function
# is locally constant for that component) are returned with a zero
# Hessian. Without this, forward-FD noise in the Jacobian can be
# amplified to O(1) when divided by the outer step.
lav_func_hessian_complex <- function(func, x,                          # nolint start
                                     h = sqrt(.Machine$double.eps), ...,
                                     fallback.simple = TRUE) {         # nolint end
  # Jacobian at x
  j0 <- try(lav_func_jacobian_complex(func = func, x = x, ...,
                                      fallback.simple = fallback.simple),
            silent = TRUE)
  if (inherits(j0, "try-error")) {
    return(j0)
  }
  nres <- nrow(j0)
  nvar <- length(x)

  # outer step per element of x
  h_vec <- pmax(h, abs(h * x))
  # exact h per element (avoid floating-point asymmetry)
  tmp <- x + h_vec
  h_vec <- tmp - x

  # detect components whose function is (locally) constant: their
  # Jacobian row is essentially zero. We will zero their Hessians too.
  # The threshold accommodates numerical noise from the fallback simple
  # forward FD, which has O(sqrt(eps)) accuracy per element.
  j0_row_max <- apply(abs(j0), 1, max)
  constant_idx <- which(j0_row_max < .Machine$double.eps^(1 / 3))

  # build Hessians by finite-differencing the Jacobian
  h_list <- vector("list", nres)
  for (i in seq_len(nres)) {
    h_list[[i]] <- matrix(0, nrow = nvar, ncol = nvar)
  }
  for (j in seq_len(nvar)) {
    x_plus <- x
    x_plus[j] <- x[j] + h_vec[j]
    jp <- try(lav_func_jacobian_complex(func = func, x = x_plus, ...,
                                        fallback.simple = fallback.simple),
              silent = TRUE)
    if (inherits(jp, "try-error")) {
      return(jp)
    }
    for (i in seq_len(nres)) {
      h_list[[i]][, j] <- (jp[i, ] - j0[i, ]) / h_vec[j]
    }
  }

  # symmetrize (small numerical asymmetry from finite differences)
  for (i in seq_len(nres)) {
    h_list[[i]] <- 0.5 * (h_list[[i]] + t(h_list[[i]]))
  }

  # zero out Hessians for locally-constant components
  for (i in constant_idx) {
    h_list[[i]][] <- 0
  }

  h_list
}

# 'simple' version (forward finite differences only), used as a fallback
# when the complex-step Jacobian cannot be used (e.g. functions that
# do not accept complex arguments). See lav_func_hessian_complex() for
# the rationale behind zeroing Hessians of locally-constant components.
lav_func_hessian_simple <- function(func, x,
                                    h = .Machine$double.eps^(1 / 4), ...) {
  j0 <- lav_func_jacobian_simple(func = func, x = x, ...)
  nres <- nrow(j0)
  nvar <- length(x)

  h_vec <- pmax(h, abs(h * x))
  tmp <- x + h_vec
  h_vec <- tmp - x

  # detect components whose function is (locally) constant (see complex
  # variant for rationale)
  j0_row_max <- apply(abs(j0), 1, max)
  constant_idx <- which(j0_row_max < .Machine$double.eps^(1 / 3))

  h_list <- vector("list", nres)
  for (i in seq_len(nres)) {
    h_list[[i]] <- matrix(0, nrow = nvar, ncol = nvar)
  }
  for (j in seq_len(nvar)) {
    x_plus <- x
    x_plus[j] <- x[j] + h_vec[j]
    jp <- lav_func_jacobian_simple(func = func, x = x_plus, ...)
    for (i in seq_len(nres)) {
      h_list[[i]][, j] <- (jp[i, ] - j0[i, ]) / h_vec[j]
    }
  }

  for (i in seq_len(nres)) {
    h_list[[i]] <- 0.5 * (h_list[[i]] + t(h_list[[i]]))
  }

  for (i in constant_idx) {
    h_list[[i]][] <- 0
  }

  h_list
}

lav_deriv_cov2cor_b <- function(m_cov = NULL) {
  nvar <- nrow(m_cov)
  ds_inv <- 1 / diag(m_cov)
  m_r <- cov2cor(m_cov)
  m_a <- -m_r %x% (0.5 * diag(ds_inv))
  m_b <- (0.5 * diag(ds_inv)) %x% -m_r
  m_dd <- diag(lav_matrix_vec(diag(nvar)))
  a2 <- m_a %*% m_dd
  b2 <- m_b %*% m_dd
  out <- a2 + b2 + diag(lav_matrix_vec(tcrossprod(sqrt(ds_inv))))
  m_d <- lav_matrix_duplication(nvar)
  out_vech <- 0.5 * (t(m_d) %*% out %*% m_d)
  out_vech
}

# quick and dirty (FIXME!!!) way to get
# surely there must be a more elegant way?
# see lav_deriv_cov2cor_b, if no num.idx...
# dCor/dCov
lav_deriv_cov2cor <- function(cov_1 = NULL, num_idx = NULL) {
  # dCor/dvar1 = - cov / (2*var1 * sqrt(var1) * sqrt(var2))
  # dCor/dvar2 = - cov / (2*var2 * sqrt(var1) * sqrt(var2))
  # dCor/dcov  =  1/(sqrt(var1) * sqrt(var2))

  # diagonal: diag(lav_matrix_vech(tcrossprod(1/delta)))

  nvar <- ncol(cov_1)
  pstar <- nvar * (nvar + 1) / 2
  delta <- sqrt(diag(cov_1))
  if (length(num_idx) > 0L) {
    delta[num_idx] <- 1.0
  }

  a <- cov_1 * -1 / (2 * delta * delta * tcrossprod(delta))
  if (length(num_idx) > 0L) {
    a[num_idx, ] <- 0
    a[cbind(num_idx, num_idx)] <- 1
  }
  a2 <- diag(nvar) %x% t(a)

  out <- diag(pstar)
  diag(out) <- lav_matrix_vech(tcrossprod(1 / delta))
  var_idx <- lav_matrix_diagh_idx(nvar)
  dup <- lav_matrix_duplication(nvar)
  out[, var_idx] <- t(dup) %*% a2[, lav_matrix_diag_idx(nvar)]

  if (length(num_idx) > 0L) {
    var_idx <- var_idx[-num_idx]
  }
  out[var_idx, var_idx] <- 0

  out
}
