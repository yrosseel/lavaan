# utility functions
#
# initial version: YR 25/03/2009

# outlier detection based on inter-quartile range
# same as boxplot.stats, but returning the indices (not the values)
lav_sample_outlier_idx <- function (x, coef = 1.5) {
  if (coef < 0)
   lav_msg_stop(gettext("'coef' must not be negative"))
  stats <- stats::fivenum(x, na.rm = TRUE)
  iqr <- diff(stats[c(2, 4)])
  if (coef == 0)
    return(seq_len(length(x)))
  else {
    out <- if (!is.na(iqr)) {
      which(x < (stats[2L] - coef * iqr) | x > (stats[4L] + coef * iqr))
    }
    else which(!is.finite(x))
  }
  out
}

# sd with trimming
lav_sample_trimmed_sd <- function(x, na.rm = TRUE, trim = 0) {
  if (isTRUE(na.rm))
    x <- x[!is.na(x)]
  n <- length(x)
  if (trim > 0 && n) {
      if (is.complex(x))
          lav_msg_stop(gettext("trimmed means are not defined for complex data"))
      if (anyNA(x))
          return(NA_real_)
      if (trim >= 0.5)
          return(stats::median(x, na.rm = FALSE))
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
        "problem computing distances: could not invert Sigma"))
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
      "sds had 0 or NA entries; non-finite result is doubtful"))
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
