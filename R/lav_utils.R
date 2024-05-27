# utility functions
#
# initial version: YR 25/03/2009

# get 'test'
# make sure we return a single element
lav_utils_get_test <- function(lavobject) {
  test <- lavobject@Options$test
  # 0.6.5: for now, we make sure that 'test' is a single element
  if (length(test) > 1L) {
    standard.idx <- which(test == "standard")
    if (length(standard.idx) > 0L) {
      test <- test[-standard.idx]
    }
    if (length(test) > 1L) {
      # only retain the first one
      test <- test[1]
    }
  }

  test
}

# check if we use a robust/scaled test statistic
lav_utils_get_scaled <- function(lavobject) {
  test.names <- unname(sapply(lavobject@test, "[[", "test"))
  scaled <- FALSE
  if (any(test.names %in% c(
    "satorra.bentler",
    "yuan.bentler", "yuan.bentler.mplus",
    "mean.var.adjusted", "scaled.shifted"
  ))) {
    scaled <- TRUE
  }

  scaled
}

# check for marker indicators:
#   - if std.lv = FALSE: a single '1' per factor, everything else zero
#   - if std.lv = TRUE: a single non-zero value per factor, everything else zero
lav_utils_get_marker <- function(LAMBDA = NULL, std.lv = FALSE) {
  LAMBDA <- as.matrix(LAMBDA)
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

  # round values
  LAMBDA <- round(LAMBDA, 3L)

  marker.idx <- numeric(nfac)
  for (f in seq_len(nfac)) {
    if (std.lv) {
      marker.idx[f] <- which(rowSums(cbind(
        LAMBDA[, f] != 0,
        LAMBDA[, -f] == 0
      )) == nfac)[1]
    } else {
      marker.idx[f] <- which(rowSums(cbind(
        LAMBDA[, f] == 1,
        LAMBDA[, -f] == 0
      )) == nfac)[1]
    }
  }

  marker.idx
}


# get npar (taking into account explicit equality constraints)
# (changed in 0.5-13)
lav_utils_get_npar <- function(lavobject) {
  npar <- lav_partable_npar(lavobject@ParTable)
  if (nrow(lavobject@Model@con.jac) > 0L) {
    ceq.idx <- attr(lavobject@Model@con.jac, "ceq.idx")
    if (length(ceq.idx) > 0L) {
      neq <- qr(lavobject@Model@con.jac[ceq.idx, , drop = FALSE])$rank
      npar <- npar - neq
    }
  } else if (.hasSlot(lavobject@Model, "ceq.simple.only") &&
    lavobject@Model@ceq.simple.only) {
    npar <- lavobject@Model@nx.free
  }

  npar
}

# N versus N-1 (or N versus N-G in the multiple group setting)
# Changed 0.5-15: suggestion by Mark Seeto
lav_utils_get_ntotal <- function(lavobject) {
  if (lavobject@Options$estimator %in% c("ML", "PML", "FML", "catML") &&
    lavobject@Options$likelihood %in% c("default", "normal")) {
    N <- lavobject@SampleStats@ntotal
  } else {
    N <- lavobject@SampleStats@ntotal - lavobject@SampleStats@ngroups
  }

  N
}

# compute log(sum(exp(x))) avoiding under/overflow
# using the identity: log(sum(exp(x)) = a + log(sum(exp(x - a)))
lav_utils_logsumexp <- function(x) {
  a <- max(x)
  a + log(sum(exp(x - a)))
}

# mdist = Mahalanobis distance
lav_mdist <- function(Y, Mp = NULL, wt = NULL,
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


# create matrix with indices to reconstruct the bootstrap samples
# per group
# (originally needed for BCa confidence intervals)
#
# rows are the (R) bootstrap runs
# columns are the (N) observations
#
# simple version: no strata, no weights
#
lav_utils_bootstrap_indices <- function(R = 0L,
                                        nobs = list(0L), # per group
                                        parallel = "no",
                                        ncpus = 1L,
                                        cl = NULL,
                                        iseed = NULL,
                                        merge.groups = FALSE,
                                        return.freq = FALSE) {
  # iseed must be set!
  stopifnot(!is.null(iseed))

  if (return.freq && !merge.groups) {
    lav_msg_stop(gettext("return.freq only available if merge.groups = TRUE"))
  }

  if (is.integer(nobs)) {
    nobs <- list(nobs)
  }

  # number of groups
  ngroups <- length(nobs)

  # mimic 'random' sampling from lav_bootstrap_internal:

  # the next 7 lines are borrowed from the boot package
  have_mc <- have_snow <- FALSE
  parallel <- parallel[1]
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) ncpus <- 1L
    loadNamespace("parallel") # before recording seed!
  }
  temp.seed <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    temp.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  if (!(ncpus > 1L && (have_mc || have_snow))) { # Only for serial
    set.seed(iseed)
  }

  # fn() returns indices per group
  fn <- function(b) {
    BOOT.idx <- vector("list", length = ngroups)
    OFFSet <- cumsum(c(0, unlist(nobs)))
    for (g in 1:ngroups) {
      stopifnot(nobs[[g]] > 1L)
      boot.idx <- sample.int(nobs[[g]], replace = TRUE)
      if (merge.groups) {
        BOOT.idx[[g]] <- boot.idx + OFFSet[g]
      } else {
        BOOT.idx[[g]] <- boot.idx
      }
    }
    BOOT.idx
  }

  RR <- R
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      RNGkind_old <- RNGkind() # store current kind
      RNGkind("L'Ecuyer-CMRG") # to allow for reproducible results
      set.seed(iseed)
      parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
    } else if (have_snow) {
      # list(...) # evaluate any promises
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        parallel::clusterSetRNGStream(cl, iseed = iseed)
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      } else {
        parallel::parLapply(cl, seq_len(RR), fn)
      }
    }
  } else {
    lapply(seq_len(RR), fn)
  }

  # restore old RNGkind()
  if (ncpus > 1L && have_mc) {
    RNGkind(RNGkind_old[1], RNGkind_old[2], RNGkind_old[3])
  }

  # handle temp.seed
  if (!is.null(temp.seed) && !identical(temp.seed, NA)) {
    assign(".Random.seed", temp.seed, envir = .GlobalEnv)
  } else if (is.null(temp.seed) && !(ncpus > 1L && (have_mc || have_snow))) {
    # serial
    rm(.Random.seed, pos = 1)
  } else if (is.null(temp.seed) && (ncpus > 1L && have_mc)) {
    # parallel/multicore only
    rm(.Random.seed, pos = 1) # because set used set.seed()
  }


  # assemble IDX
  BOOT.idx <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    # FIXME: handle failed runs
    BOOT.idx[[g]] <- do.call("rbind", lapply(res, "[[", g))
  }

  # merge groups
  if (merge.groups) {
    out <- do.call("cbind", BOOT.idx)
  } else {
    out <- BOOT.idx
  }

  # NOTE: the order of the indices is different from the boot package!
  # we fill in the matrix 'row-wise' (1 row = sample(N, replace = TRUE)),
  # while boot fills in the matrix 'column-wise'
  # this also explains why we get different results with return.boot = TRUE
  # despite using the same iseed

  # return frequencies instead?
  if (return.freq && merge.groups) {
    out <- t(apply(out, 1L, tabulate, ncol(out)))
  }

  out
}


# invert positive definite symmetric matrix (eg cov matrix)
# using choleski decomposition
# return log determinant as an attribute
inv.chol <- function(S, logdet = FALSE) {
  cS <- chol(S)
  # if( inherits(cS, "try-error") ) {
  #    print(S)
  #    warning("lavaan WARNING: symmetric matrix is not positive symmetric!")
  # }
  S.inv <- chol2inv(cS)
  if (logdet) {
    diag.cS <- diag(cS)
    attr(S.inv, "logdet") <- sum(log(diag.cS * diag.cS))
  }
  S.inv
}

# convert correlation matrix + standard deviations to covariance matrix
# based on cov2cor in package:stats
cor2cov <- function(R, sds, names = NULL) {
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
#     x <- char2num(s)
char2num <- function(s = "") {
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
getCov <- function(x, lower = TRUE, diagonal = TRUE, sds = NULL,
                   names = paste("V", 1:nvar, sep = "")) {
  # check x and sds
  if (is.character(x)) x <- char2num(x)
  if (is.character(sds)) sds <- char2num(sds)

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
    COV <- cor2cov(COV, sds)
  }

  # names
  stopifnot(length(names) == nvar)
  rownames(COV) <- colnames(COV) <- names

  COV
}


# translate row+col matrix indices to vec idx
rowcol2vec <- function(row.idx, col.idx, nrow, symmetric = FALSE) {
  idx <- row.idx + (col.idx - 1) * nrow
  if (symmetric) {
    idx2 <- col.idx + (row.idx - 1) * nrow
    idx <- unique(sort(c(idx, idx2)))
  }
  idx
}

# dummy function to 'pretty' print a vector with fixed width
pprint.vector <- function(x,
                          digits.after.period = 3,
                          ncols = NULL, max.col.width = 11,
                          newline = TRUE) {
  n <- length(x)
  var.names <- names(x)

  total.width <- getOption("width")

  max.width <- max(nchar(var.names))
  if (max.width < max.col.width) { # shrink
    max.col.width <- max(max.width, digits.after.period + 2)
  }

  # automatic number of columns
  if (is.null(ncols)) {
    ncols <- floor((total.width - 2) / (max.col.width + 2))
  }
  nrows <- ceiling(n / ncols)

  if (digits.after.period >= (max.col.width - 3)) {
    max.col.width <- digits.after.period + 3
  }
  string.format <- paste(" %", max.col.width, "s", sep = "")
  number.format <- paste(" %", max.col.width, ".", digits.after.period, "f", sep = "")

  for (nr in 1:nrows) {
    rest <- min(ncols, n)
    if (newline) cat("\n")
    # labels
    for (nc in 1:rest) {
      vname <- substr(var.names[(nr - 1) * ncols + nc], 1, max.col.width)
      cat(sprintf(string.format, vname))
    }
    cat("\n")
    for (nc in 1:rest) {
      cat(sprintf(number.format, x[(nr - 1) * ncols + nc]))
    }
    cat("\n")
    n <- n - ncols
  }
  if (newline) cat("\n")
}

# print only lower half of symmetric matrix
pprint.matrix.symm <- function(x,
                               digits.after.period = 3,
                               ncols = NULL, max.col.width = 11,
                               newline = TRUE) {
  n <- ncol <- ncol(x)
  nrow <- nrow(x)
  stopifnot(ncol == nrow)
  var.names <- rownames(x)

  total.width <- getOption("width")

  max.width <- max(nchar(var.names))
  if (max.width < max.col.width) { # shrink
    max.col.width <- max(max.width, digits.after.period + 2)
  }

  # automatic number of columns
  if (is.null(ncols)) {
    ncols <- floor((total.width - 2) / (max.col.width + 2))
  }

  nblocks <- ceiling(n / ncols)

  if (digits.after.period >= (max.col.width - 3)) {
    max.col.width <- digits.after.period + 3
  }
  fc.format <- paste(" %", min(max.width, max.col.width), "s", sep = "")
  string.format <- paste(" %", max.col.width, "s", sep = "")
  number.format <- paste(" %", max.col.width, ".", digits.after.period, "f", sep = "")

  for (nb in 1:nblocks) {
    rest <- min(ncols, n)
    if (newline) cat("\n")
    # empty column
    cat(sprintf(fc.format, ""))
    # labels
    for (nc in 1:rest) {
      vname <- substr(var.names[(nb - 1) * ncols + nc], 1, max.col.width)
      cat(sprintf(string.format, vname))
    }
    cat("\n")
    row.start <- (nb - 1) * ncols + 1
    for (nr in row.start:nrow) {
      # label
      vname <- substr(var.names[nr], 1, max.col.width)
      cat(sprintf(fc.format, vname))
      col.rest <- min(rest, (nr - row.start + 1))
      for (nc in 1:col.rest) {
        value <- x[nr, (nb - 1) * ncols + nc]
        cat(sprintf(number.format, value))
      }
      cat("\n")
    }
    n <- n - ncols
  }
  if (newline) cat("\n")
}

# elimination of rows/cols symmetric matrix
eliminate.rowcols <- function(x, el.idx = integer(0)) {
  if (length(el.idx) == 0) {
    return(x)
  }
  stopifnot(ncol(x) == nrow(x))
  stopifnot(min(el.idx) > 0 && max(el.idx) <= ncol(x))

  x[-el.idx, -el.idx]
}

# locate pstar idx of a given set of elements
#
# for example, if nvar = 4, pstar = 10 elements
# which elements correponds to the second and third element (in 1:nvar)?
#
# deprecated! replaced by lav_matrix_vech_which_idx()
#
eliminate.pstar.idx2 <- function(nvar = 1L, el.idx = integer(0),
                                   meanstructure = FALSE,
                                   correlation = FALSE,
                                   return.idx = FALSE) {
  if (length(el.idx) > 0) {
    stopifnot(min(el.idx) > 0 && max(el.idx) <= nvar)
  }

  # create col/row indices
  XX <- rbind(
    lav_matrix_vech_col_idx(nvar),
    lav_matrix_vech_row_idx(nvar)
  )

  # if correlation matrix, remove col/row corresponding to the variances
  if (correlation && nvar > 1L) {
    var.idx <- lav_matrix_diagh_idx(n = nvar)
    XX <- XX[, -var.idx, drop = FALSE]
  }

  # locate pstar indices (as logicals)
  idx <- XX[1, ] %in% el.idx & XX[2, ] %in% el.idx

  # if meanstructure, add location in mean vector
  if (meanstructure) {
    idx <- c((1:nvar %in% el.idx), idx)
  }

  # return indices instead of logicals
  if (return.idx) {
    idx <- which(idx)
  }

  idx
}

# elimination of rows/cols pstar symmetric matrix
#
# type = "all" -> only remove var(el.idx) and cov(el.idx)
# type = "any" -> remove all rows/cols of el.idx
eliminate.pstar.idx <- function(nvar = 1, el.idx = integer(0),
                                meanstructure = FALSE, type = "all") {
  if (length(el.idx) > 0) {
    stopifnot(min(el.idx) > 0 && max(el.idx) <= nvar)
  }

  XX <- utils::combn(1:(nvar + 1), 2)
  XX[2, ] <- XX[2, ] - 1

  if (type == "all") {
    idx <- !(apply(apply(XX, 2, function(x) {
      x %in% el.idx
    }), 2, all))
  } else {
    idx <- !(apply(apply(XX, 2, function(x) {
      x %in% el.idx
    }), 2, any))
  }

  if (meanstructure) {
    idx <- c(!(1:nvar %in% el.idx), idx)
    # idx <- c(rep(TRUE, nvar), idx)
  }

  idx
}

# construct 'augmented' covariance matrix
# based on the covariance matrix and the mean vector
augmented.covariance <- function(S., mean) {
  S <- as.matrix(S.)
  m <- as.matrix(mean)
  p <- ncol(S)

  if (nrow(m) != p) {
    lav_msg_stop(gettext("incompatible dimension of mean vector"))
  }

  out <- matrix(0, ncol = (p + 1), nrow = (p + 1))

  out[1:p, 1:p] <- S + m %*% t(m)
  out[p + 1, 1:p] <- t(m)
  out[1:p, p + 1] <- m
  out[p + 1, p + 1] <- 1

  out
}



# linesearch using 'armijo' backtracking
# to find a suitable `stepsize' (alpha)
linesearch.backtracking.armijo <- function(f.alpha, s.alpha, alpha = 10) {
  tau <- 0.5
  ftol <- 0.001

  f.old <- f.alpha(0)
  s.old <- s.alpha(0)

  armijo.condition <- function(alpha) {
    f.new <- f.alpha(alpha)

    # condition
    f.new > f.old + ftol * alpha * s.old
  }

  i <- 1
  while (armijo.condition(alpha)) {
    alpha <- alpha * tau
    f.new <- f.alpha(alpha)
    cat("... backtracking: ", i, "alpha = ", alpha, "f.new = ", f.new, "\n")
    i <- i + 1
  }

  alpha
}

steepest.descent <- function(start, objective, gradient, iter.max, verbose) {
  x <- start

  if (verbose) {
    cat("Steepest descent iterations\n")
    cat("iter        function  abs.change  rel.change     step.size       norm.gx\n")
    gx <- gradient(x)
    norm.gx <- sqrt(gx %*% gx)
    fx <- objective(x)
    cat(sprintf(
      "%4d   %11.7E                         %11.5E %11.5E",
      0, fx, 0, norm.gx
    ), "\n")
  }

  for (iter in 1:iter.max) {
    fx.old <- objective(x)

    # normalized gradient
    gx <- gradient(x)
    old.gx <- gx
    norm.gx <- sqrt(gx %*% gx)
    gradient.old <- gx / norm.gx
    direction.vector <- (-1) * gradient.old

    f.alpha <- function(alpha) {
      new.x <- x + alpha * direction.vector
      fx <- objective(new.x)
      # cat("  [stepsize]  iter ", iter, " step size = ", alpha,
      #    " fx = ", fx, "\n", sep="")
      # for optimize only
      if (is.infinite(fx)) {
        fx <- .Machine$double.xmax
      }
      fx
    }

    # s.alpha <- function(alpha) {
    #    new.x <- x + alpha * direction.vector
    #    gradient.new <- gradient(new.x)
    #    norm.gx <- sqrt( gradient.new %*% gradient.new)
    #    gradient.new <- gradient.new/norm.gx
    #    as.numeric(gradient.new %*% direction.vector)
    # }

    # find step size
    # alpha <- linesearch.backtracking.armijo(f.alpha, s.alpha, alpha=1)

    if (iter == 1) {
      alpha <- 0.1
    } else {
      alpha <- optimize(f.alpha, lower = 0.0, upper = 1)$minimum
      if (f.alpha(alpha) > fx.old) {
        alpha <- optimize(f.alpha, lower = -1, upper = 0.0)$minimum
      }
    }


    # steepest descent step
    old.x <- x
    x <- x + alpha * direction.vector
    gx.old <- gx
    gx <- gradient(x)
    dx.max <- max(abs(gx))

    # verbose
    if (verbose) {
      fx <- fx.old
      fx.new <- objective(x)
      abs.change <- fx.new - fx.old
      rel.change <- abs.change / fx.old
      norm.gx <- sqrt(gx %*% gx)
      if (verbose) {
        cat(
          sprintf(
            "%4d   %11.7E %10.7f %10.7f %11.5E %11.5E",
            iter, fx.new, abs.change, rel.change, alpha, norm.gx
          ),
          "\n"
        )
      }
    }

    # convergence check
    if (dx.max < 1e-05) {
      break
    }
  }

  x
}
