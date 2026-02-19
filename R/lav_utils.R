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
