# utility functions
#
# initial version: YR 25/03/2009

# multivariate normal random number generation
# replacement for MASS::mvrnorm for better cross-machine reproducibility
# see: https://blog.djnavarro.net/posts/2025-05-18_multivariate-normal-sampling-floating-point/
#
# the issue with MASS::mvrnorm is that it uses eigendecomposition with a
# transformation matrix (sqrt(Lambda) Q') that is NOT invariant to eigenvector
# sign flips. this leads to irreproducible results across machines even with
# the same random seed. this implementation follows mvtnorm::rmvnorm which uses:
# - eigen: Q sqrt(Lambda) Q' which IS invariant to sign flips (default)
# - svd: similar to eigen, also invariant
# - chol: unique factorization
# eigen is the default to match MASS::mvrnorm's approach while fixing the
# sign-flip issue with the invariant formula
lav_mvrnorm <- function(n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE,
                        method = "eigen", checkSymmetry = TRUE, byrow = FALSE) {
  p <- length(mu)

  # check symmetry (like mvtnorm)
  if (checkSymmetry && !isSymmetric(Sigma,
    tol = sqrt(.Machine$double.eps),
    check.attributes = FALSE
  )) {
    lav_msg_stop(gettext("'Sigma' must be a symmetric matrix in lav_mvrnorm"))
  }

  # check dimensions (like mvtnorm)
  if (p != nrow(Sigma)) {
    lav_msg_stop(gettext("incompatible arguments in lav_mvrnorm"))
  }

  # compute transformation matrix R based on method (following mvtnorm exactly)
  method <- match.arg(method, c("chol", "eigen", "svd"))

  R <- if (method == "eigen") {
    ev <- eigen(Sigma, symmetric = TRUE)
    if (!all(ev$values >= -tol * abs(ev$values[1L]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    # Q sqrt(Lambda) Q' - invariant to sign flips
    t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  } else if (method == "svd") {
    s. <- svd(Sigma)
    if (!all(s.$d >= -tol * abs(s.$d[1L]))) {
      warning("sigma is numerically not positive semidefinite")
    }
    t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
  } else if (method == "chol") {
    R <- chol(Sigma, pivot = TRUE)
    R[, order(attr(R, "pivot"))]
  }

  # for names (fallback to dimnames of Sigma if mu has no names)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) {
    nm <- dn[[1L]]
  }

  if (empirical) {
    # generate standard normal, then apply empirical transformation
    X <- matrix(stats::rnorm(p * n), n, p, byrow = byrow)
    X <- scale(X, center = TRUE, scale = FALSE) # center
    X <- X %*% svd(X, nu = 0)$v # orthogonalize
    X <- scale(X, center = FALSE, scale = TRUE) # unit variance

    # transform by R
    X <- sweep(X %*% R, 2, mu, "+")
    colnames(X) <- nm

    if (n == 1) drop(X) else X
  } else {
    # generate samples (following mvtnorm exactly)
    X <- matrix(stats::rnorm(n * p), nrow = n, byrow = byrow) %*% R
    X <- sweep(X, 2, mu, "+")
    colnames(X) <- nm

    if (n == 1) drop(X) else X
  }
}

# outlier detection based on inter-quartile range
# same as boxplot.stats, but returning the indices (not the values)
lav_sample_outlier_idx <- function(x, coef = 1.5) {
  if (coef < 0) {
    lav_msg_stop(gettext("'coef' must not be negative"))
  }
  stats <- stats::fivenum(x, na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  if (coef == 0) {
    return(seq_len(length(x)))
  } else {
    out <- if (!is.na(iqr)) {
      which(x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef * iqr))
    } else {
      which(!is.finite(x))
    }
  }
  out
}

# sd with trimming
lav_sample_trimmed_sd <- function(x, na.rm = TRUE, trim = 0) {
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (trim > 0 && n) {
    if (is.complex(x)) {
      lav_msg_stop(gettext("trimmed means are not defined for complex data"))
    }
    if (anyNA(x)) {
      return(NA_real_)
    }
    if (trim >= 0.5) {
      return(stats::median(x, na.rm = FALSE))
    }
    lo <- floor(n * trim) + 1
    hi <- n + 1 - lo
    x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}

# mdist = Mahalanobis distance
lav_sample_mdist <- function(Y, Mp = NULL, wt = NULL,
                             Mu = NULL, Sigma = NULL,
                             Sinv.method = "eigen", ginv = TRUE,
                             rescale = FALSE) {
  # check input
  Y <- as.matrix(Y)
  P <- NCOL(Y)
  if (!is.null(wt)) {
    N <- sum(wt)
  } else {
    N <- NROW(Y)
  }
  NY <- NROW(Y)

  # missing data?
  missing.flag <- anyNA(Y)

  # missing patterns?
  if (missing.flag && is.null(Mp)) {
    Mp <- lav_data_missing_patterns(Y)
  }

  # no Mu? compute sample mean
  if (is.null(Mu)) {
    Mu <- colMeans(Y, na.rm = TRUE)
  }

  # no Sigma?
  if (is.null(Sigma)) {
    if (missing.flag) {
      out <- lav_mvnorm_missing_h1_estimate_moments(
        Y = Y, Mp = Mp,
        wt = wt
      )
      Mu <- out$Mu
      Sigma <- out$Sigma
    } else {
      if (!is.null(wt)) {
        out <- stats::cov.wt(Y, wt = wt, method = "ML")
        Sigma <- out$cov
        Mu <- out$center
      } else {
        Sigma <- stats::cov(Y, use = "pairwise")
        # rescale?
        if (rescale) {
          Sigma <- ((N - 1) / N) * Sigma
        }
      }
    }
  }

  # subtract Mu
  Yc <- t(t(Y) - Mu)

  # DIST per case
  DIST <- rep(as.numeric(NA), NY)

  # invert Sigma
  if (ginv) {
    Sigma.inv <- MASS::ginv(Sigma)
  } else {
    Sigma.inv <-
      try(
        lav_matrix_symmetric_inverse(
          S = Sigma, logdet = FALSE,
          Sinv.method = Sinv.method
        ),
        silent = TRUE
      )
    if (inherits(Sigma.inv, "try-error")) {
      lav_msg_warn(gettext(
        "problem computing distances: could not invert Sigma"
      ))
      return(DIST)
    }
  }

  # complete data?
  if (!missing.flag) {
    # center factor scores
    Y.c <- t(t(Y) - Mu)
    # Mahalobis distance
    DIST <- rowSums((Y.c %*% Sigma.inv) * Y.c)

    # missing data?
  } else {
    # for each pattern, compute sigma.inv; compute DIST for all
    # observations of this pattern
    for (p in seq_len(Mp$npatterns)) {
      # observed values for this pattern
      var.idx <- Mp$pat[p, ]

      # missing values for this pattern
      na.idx <- which(!var.idx)

      # identify cases with this pattern
      case.idx <- Mp$case.idx[[p]]

      # invert Sigma for this pattern
      if (length(na.idx) > 0L) {
        if (ginv) {
          sigma.inv <- MASS::ginv(Sigma[-na.idx, -na.idx, drop = FALSE])
        } else {
          sigma.inv <-
            lav_matrix_symmetric_inverse_update(
              S.inv = Sigma.inv,
              rm.idx = na.idx, logdet = FALSE
            )
        }
      } else {
        sigma.inv <- Sigma.inv
      }

      if (Mp$freq[p] == 1L) {
        DIST[case.idx] <- sum(sigma.inv *
          crossprod(Yc[case.idx, var.idx, drop = FALSE]))
      } else {
        DIST[case.idx] <-
          rowSums(Yc[case.idx, var.idx, drop = FALSE] %*% sigma.inv *
            Yc[case.idx, var.idx, drop = FALSE])
      }
    } # patterns
  } # missing data

  # use weights? (no for now)
  # DIST <- DIST * wt

  DIST
}

# convert correlation matrix + standard deviations to covariance matrix
# based on cov2cor in package:stats
lav_cor2cov <- function(R, sds, names = NULL) {
  p <- (d <- dim(R))[1L]
  if (!is.numeric(R) || length(d) != 2L || p != d[2L]) {
    lav_msg_stop(gettext("'V' is not a square numeric matrix"))
  }

  if (any(!is.finite(sds))) {
    lav_msg_warn(gettext(
      "sds had 0 or NA entries; non-finite result is doubtful"
    ))
  }

  # if(sum(diag(R)) != p)
  #    stop("The diagonal of a correlation matrix should be all ones.")

  if (p != length(sds)) {
    lav_msg_stop(gettext("The standard deviation vector and correlation matrix
                         have a different number of variables"))
  }

  S <- R
  S[] <- sds * R * rep(sds, each = p)

  # optionally, add names
  if (!is.null(names)) {
    stopifnot(length(names) == p)
    rownames(S) <- colnames(S) <- names
  }

  S
}

# convert characters within single quotes to numeric vector
# eg. s <- '3 4.3 8e-3 2.0'
#     x <- lav_char2num(s)
lav_char2num <- function(s = "") {
  # first, strip all ',' or ';'
  s. <- gsub(",", " ", s)
  s. <- gsub(";", " ", s.)
  tc <- textConnection(s.)
  x <- scan(tc, quiet = TRUE)
  close(tc)
  x
}

# create full matrix based on lower.tri or upper.tri elements; add names
# always ROW-WISE!!
lav_getcov <- function(x, lower = TRUE, diagonal = TRUE, sds = NULL,
                       names = paste("V", 1:nvar, sep = "")) {
  # check x and sds
  if (is.character(x)) x <- lav_char2num(x)
  if (is.character(sds)) sds <- lav_char2num(sds)

  nels <- length(x)
  if (lower) {
    COV <- lav_matrix_lower2full(x, diagonal = diagonal)
  } else {
    COV <- lav_matrix_upper2full(x, diagonal = diagonal)
  }
  nvar <- ncol(COV)

  # if diagonal is false, assume unit diagonal
  if (!diagonal) diag(COV) <- 1

  # check if we have a sds argument
  if (!is.null(sds)) {
    stopifnot(length(sds) == nvar)
    COV <- lav_cor2cov(COV, sds)
  }

  # names
  stopifnot(length(names) == nvar)
  rownames(COV) <- colnames(COV) <- names

  COV
}

lav_char2hash <- function(s = "") {
  stopifnot(is.character(s))
  nums <- utf8ToInt(paste(s, collapse = "\n"))
  rval <- 0x7EDCBA98L
  for (i in nums) {
    positions <- 1L + (bitwAnd(i, 15L))
    bitsR <- bitwShiftL(1L, positions) - 1L
    wrap1 <- bitwShiftL(bitwAnd(bitsR, rval), 31L - positions)
    wrap2 <- bitwShiftR(bitwAnd(bitwNot(bitsR), rval), positions)
    saw <- bitwOr(wrap1, wrap2)
    rval <- bitwXor(saw, i)
  }
  as.hexmode(rval)
}

# vectorize all (h0 or h1) sample statistics, in the same order
# as Gamma
lav_implied_to_vec <- function(implied = NULL, lavmodel = NULL,
                               drop.list = TRUE) {
  ngroups <- lavmodel@ngroups

  wls_obs <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {

    var <- NULL
    res.var <- NULL
    if (lavmodel@conditional.x) {
      res.var <- diag(implied$res.cov[[g]])
      res.var <- res.var[res.var != 1]
    } else {
      var <- diag(implied$cov[[g]])
      var <- var[var != 1]
    }

    wls_obs[[g]] <- lav_samplestats_wls_obs(
      # plain
      mean.g = implied$mean[[g]],
      cov.g = implied$cov[[g]],
      var.g = var,
      th.g = implied$th[[g]],
      th.idx.g = lavmodel@th.idx[[g]],

      # conditional.x
      res.int.g = implied$res.int[[g]],
      res.cov.g = implied$res.cov[[g]],
      res.var.g = res.var,
      res.th.g = implied$res.th[[g]],
      res.slopes.g = implied$res.slopes[[g]],
      group.w.g = implied$group.w[[g]],

      # flags
      categorical = lavmodel@categorical,
      conditional.x = lavmodel@conditional.x,
      meanstructure = lavmodel@meanstructure,
      correlation = lavmodel@correlation,
      slopestructure = lavmodel@conditional.x,
      group.w.free = lavmodel@group.w.free
    )
  }

  if (drop.list) {
    out <- unlist(wls_obs)
  } else {
    out <- wls_obs
  }

  out
}

# given a vector of sample statistics, reconstruct implied
# entries (cov, mean, th, res.cov, ...)
lav_vec_to_implied <- function(x = NULL, lavmodel) {

  ngroups <- lavmodel@ngroups
  implied <- list()

  for (g in seq_len(ngroups)) {

    # number of variables
    nvar <- lavmodel@nvar[g]

    # if group.w.free, always comes first
    if (lavmodel@group.w.free) {
      idx <- 1L
      group.w <- x[idx]
      x <- x[-idx]
    } else {
      group.w <- 1
    }

    if (lavmodel@categorical) {
      if (lavmodel@conditional.x) {
        # TODO
        cat("\n NOT READY YET! \n")
      }

      # th
      nth <- length(lavmodel@th.idx[[g]])
      idx <- seq_len(nth)
      th_g <- x[idx]
      x <- x[-idx]

      # mean/var ov.num
      mean_g <- rep(0, nvar)
      var_g <- rep(1, nvar)
      if (any(lavmodel@th.idx[[g]] == 0)) {
        num.idx <- seq_len(nvar)
        num.idx <- num.idx[!num.idx %in% lavmodel@th.idx[[g]]]
        idx <- seq_len(length(num.idx))
        var_g[num.idx] <- x[idx]
        x <- x[-idx]
        # FIXME: change sign?
        mean_g[num.idx] <- th_g[lavmodel@th.idx[[g]] == 0]
      }

      # cov - lower only
      idx <- seq_len( nvar * (nvar - 1) / 2 )
      cov_g <- lav_matrix_vech_reverse(x[idx], diagonal = FALSE)
      x <- x[-idx]
      diag(cov_g) <- var_g

      # fill in
      implied$cov[[g]] <- cov_g
      implied$mean[[g]] <- mean_g
      implied$th[[g]] <- th_g

    } else {
      # diag flag
      diag_flag <- TRUE
      if (lavmodel@correlation) {
        diag_flag <- FALSE
      }
      if (lavmodel@conditional.x) {
        # TODO
        cat("\n NOT READY YET! \n")
      }
      if (lavmodel@meanstructure) {
        idx <- seq_len(nvar)
        mean_g <- x[idx]
        x <- x[-idx]
      } else {
        mean_g <- numeric(nvar)
      }
      if (diag_flag) {
        idx <- seq_len( nvar * (nvar + 1) / 2 )
      } else {
        idx <- seq_len( nvar * (nvar - 1) / 2 )
      }
      cov_g <- lav_matrix_vech_reverse(x[idx], diagonal = diag_flag)
      x <- x[-idx]

      # fill in
      implied$cov[[g]] <- cov_g
      implied$mean[[g]] <- mean_g
      implied$th[g] <- list(NULL)
    }

    implied$group.w[[g]] <- group.w
  }

  implied
}

# compute the 'information' matrix I = t(Delta) %*% W %*% Delta
#
# where W is a function of S.inv only, without constructing the kronecker
# product of S.inv %x% S.inv
#
# strategy:
# - use cholesky to write S.inv = t(U) %*% U
# - for the cov block: we use the fact that the ij element equals
#   1/2 * tr(S.inv X_i S.inv X_j) = 1/2 * tr(Z_i Z_j)
#   where Z_i = U X_i t(U), and X_i = unvech(jac_cov[,i])
#   since Z_i is symmetric, this equals t(h_i) %*% h_j, where
#   h_i = sqrt(w) * vech(Z_i), so we can write  I_cov = t(H) %*% H
# - for the mean block, we use I_mean = t(U J_mu) %*% (U J_mu)
#
# FIXME: what if fixed.x = TRUE?
#
lav_utils_deltat_w_delta <- function(Delta = NULL, S = NULL,
                                     meanstructure = FALSE) {
  nvar <- nrow(S)
  pstar <- nvar * (nvar + 1L) / 2L

  jac_cov <- Delta

  # split Delta
  if (meanstructure) {
    mean_idx <- seq_len(nvar)
    jac_mean <- Delta[ mean_idx, , drop = FALSE]
    jac_cov <- Delta[-mean_idx, , drop = FALSE]
  } else if (nrow(Delta) != pstar) {
    lav_msg_stop(gettext("nrow(Delta) != pstar"))
  }

  n_free  <- ncol(jac_cov)

  # vech
  r_s <- lav_matrix_vech_row_idx(nvar)
  c_s <- lav_matrix_vech_col_idx(nvar)
  sqrt_w <- ifelse(r_s == c_s, sqrt(0.5), 1.0)


  # special case: S = I (ULS)
  uls.flag <- isTRUE(all.equal(S, diag(nvar), check.attributes = FALSE))
  if (uls.flag) {
    H <- sqrt_w * jac_cov
    A <- crossprod(H)
    if (meanstructure) {
      A <- A + crossprod(jac_mean)
    }
  } else {
    cS <- chol(S); Si <- chol2inv(cS)
    U  <- chol(Si)

    # build mat_X (nvar^2 x n_free): vec of symmetrized X_i
    vec_lo  <- r_s + (c_s - 1L) * nvar     # lower-tri vec index
    vec_up  <- c_s + (r_s - 1L) * nvar     # upper-tri vec index
    mat_X   <- matrix(0, nvar * nvar, n_free)
    mat_X[vec_lo, ] <- jac_cov
    mat_X[vec_up, ] <- jac_cov

    # left-multiply each X_i by U
    mat_W <- U %*% matrix(mat_X, nvar, nvar * n_free)

    # right-multiply each block by t(U) using array reshape
    W_13 <- matrix(aperm(array(mat_W, c(nvar, nvar, n_free)), c(1L, 3L, 2L)),
                   nrow = nvar * n_free, ncol = nvar)
    mat_Z <- W_13 %*% t(U)

    # extract vech of each Z_i and scale by sqrt_w -> H
    i_off <- rep((seq_len(n_free) - 1L) * nvar, each = pstar)
    H <- matrix(mat_Z[cbind(rep(r_s, n_free) + i_off,
                            rep(c_s, n_free))], pstar, n_free) * sqrt_w
    A <- crossprod(H)

    if (meanstructure) {
      # t(jac_mean) %*% Si %*% jac_mean = crossprod(U %*% jac_mean)
      V <- U %*% jac_mean
      A <- A + crossprod(V)
    }
  }

  A
}


# find the solution of the linear equation:
#
# theta = solve(t(Delta) %*% W %*% Delta) %*% t(Delta) %*% W %*% svec
#
# where Delta is the Jacobian of Sigma wrt the free parameters,
# W is lav_matrix_bdiag(S.inv, S22_inv)
# where S22_inv = 0.5 * lav_matrix_duplication_pre_post(S.inv %x% S.inv)
#  - S can be I (ULS) (then W = 0.5 * t(D) %*% D)
#  - S can be sample.cov (GLS)
#  - S can be Sigma (2RLS or RLS)
#
# svec is the vector with sample statistics (first means, then vech(S))
#
# strategy:
# - first we compute Q = W %*% Delta, using the fact that W only contains S.inv
#   elements (without explicitly constructing the kronecker product)
#   note that for the cov block:
#     W %*% Delta = w * vech(S.inv * unvech(v) * S.inv)
#   where w is 1/2 for diagonal elements, and 1 elsewhere
# - compute A = t(Delta) %*% Q, and b = Q %*% svec
# - compute solve(A,b) using cholesky back/forward, as A is symmetric and pd
#
# optionally, return H = solve(t(Delta) %*% W %*% Delta) %*% t(Delta) %*% W
# as an attribute
lav_utils_wls_linearization <- function(Delta = NULL, S = NULL,
                                        meanstructure = FALSE,
                                        #fixed.x = FALSE, # FIXME: needed?
                                        #x.idx = integer(0L),
                                        svec = NULL,
                                        return_H = FALSE) {
  nvar <- nrow(S)
  pstar <- nvar * (nvar + 1L) / 2L

  jac_cov <- Delta; nc <- ncol(jac_cov)
  s_cov <- svec

  # vech
  w <- rep(1.0, pstar)
  w[lav_matrix_diagh_idx(nvar)] <- 0.5

  # split Delta/svec
  if (meanstructure) {
    mean_idx <- seq_len(nvar)
    jac_mean <- Delta[ mean_idx, , drop = FALSE]
    jac_cov <- Delta[-mean_idx, , drop = FALSE]
    s_mean <- svec[mean_idx]
    s_cov <- svec[-mean_idx]
  } else if (nrow(Delta) != pstar) {
    lav_msg_stop(gettext("nrow(Delta) != pstar"))
  }

  # special case: S = I (ULS)
  uls.flag <- isTRUE(all.equal(S, diag(nvar), check.attributes = FALSE))
  if (uls.flag) {
    # cov part
    q_cov <- w * jac_cov
    A <- crossprod(jac_cov, q_cov)
    b <- drop(crossprod(q_cov, s_cov))

    # mean part
    if (meanstructure) {
      A <- A + crossprod(jac_mean)
      b <- b + drop(crossprod(jac_mean, s_mean))
    }

    # tQ only needed if H is requested
    if (return_H) {
      tQ <- t(q_cov)
      if (meanstructure) tQ <- cbind(t(jac_mean), tQ)
    }
  } else {
    # S inverse
    cS <- chol(S); Si <- chol2inv(cS)

    # cov part
    q_cov <- matrix(0, pstar, nc)
    for (j in seq_len(nc)) {
      vmat <- lav_matrix_vech_reverse(jac_cov[, j])
      W <- Si %*% vmat %*% Si
      q_cov[, j] <- w * lav_matrix_vech(W)
    }

    A <- crossprod(jac_cov, q_cov)
    b <- drop(crossprod(q_cov, s_cov))

    # mean part
    if (meanstructure) {
      q_mean <- backsolve(cS, forwardsolve(t(cS), jac_mean))
      A <- A + crossprod(jac_mean, q_mean)
      b <- b + drop(crossprod(q_mean, s_mean))
    }

    if (return_H) {
      tQ <- t(q_cov)
      if (meanstructure) tQ <- cbind(t(q_mean), tQ)
    }
  }

  R <- chol(A)
  if (!return_H) {
    out <- backsolve(R, forwardsolve(t(R), b))
  } else {
    Ainv  <- chol2inv(R)
    out <- drop(Ainv %*% b)
    H <- Ainv %*% tQ
    attr(out, "H") <- H
  }

  out
}

