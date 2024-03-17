# various rotation criteria and their gradients
# YR 05 April 2019: initial version
# YR 14 June  2019: add more rotation criteria

# references:
#
# Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms
# and software for arbitrary rotation criteria in factor analysis. Educational
# and Psychological Measurement, 65(5), 676-696.
# old website: http://web.archive.org/web/20180708170331/http://www.stat.ucla.edu/research/gpa/splusfunctions.net
#
# Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
# analysis. Multivariate behavioral research, 36(1), 111-150.
#
# Mulaik, S. A. (2010). Foundations of factor analysis (Second Edition).
# Boca Raton: Chapman and Hall/CRC.

# Note: this is YR's implementation, not a copy of the GPArotation
#       package
#
# Why did I write my own functions (and not use the GPArotation):
# - to better understand what is going on
# - to have direct access to the gradient functions
# - to avoid yet another dependency
# - to simplify further experiments


# Orthomax family (Harman, 1960)
#
# gamma = 0   -> quartimax
# gamma = 1/2 -> biquartimax
# gamma = 1/P -> equamax
# gamma = 1   -> varimax
#
lav_matrix_rotate_orthomax <- function(LAMBDA = NULL, orthomax.gamma = 1,
                                       ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA
  # center L2 column-wise
  cL2 <- t(t(L2) - orthomax.gamma * colMeans(L2))
  out <- -1 * sum(L2 * cL2) / 4

  if (grad) {
    attr(out, "grad") <- -1 * LAMBDA * cL2
  }

  out
}

# Crawford-Ferguson (1970) family
#
# combine penalization for 1) row complexity, and 2) column complexity
# if combined with orthogonal rotation, this is equivalent to the
# orthomax family:
#
#     quartimax -> gamma = 0 (only row complexity)
#     varimax   -> gamma = 1/nrow
#     equamax   -> gamma = ncol/(2*nrow)
#     parsimax  -> gamma = (ncol - 1)/(nrow + ncol - 2)
#     factor parsimony -> gamma = 1 (only column complexity)
#
# the Crawford-Ferguson family is also equivalent to the oblimin family
# if the latter is restricted to orthogonal rotation
#
lav_matrix_rotate_cf <- function(LAMBDA = NULL, cf.gamma = 0, ...,
                                 grad = FALSE) {
  # check if gamma is between 0 and 1?
  nRow <- nrow(LAMBDA)
  nCol <- ncol(LAMBDA)
  ROW1 <- matrix(1.0, nCol, nCol)
  diag(ROW1) <- 0.0
  COL1 <- matrix(1.0, nRow, nRow)
  diag(COL1) <- 0.0

  L2 <- LAMBDA * LAMBDA
  LR <- L2 %*% ROW1
  LC <- COL1 %*% L2

  f1 <- sum(L2 * LR) / 4
  f2 <- sum(L2 * LC) / 4

  out <- (1 - cf.gamma) * f1 + cf.gamma * f2

  if (grad) {
    attr(out, "grad") <- ((1 - cf.gamma) * LAMBDA * LR) +
      (cf.gamma * LAMBDA * LC)
  }

  out
}

# Oblimin family (Carroll, 1960; Harman, 1976)
#
# quartimin   -> gamma = 0
# biquartimin -> gamma = 1/2
# covarimin   -> gamma = 1
#
# if combined with orthogonal rotation, this is equivalent to the
# orthomax family (they have the same optimizers):
#
# gamma = 0   -> quartimax
# gamma = 1/2 -> biquartimax
# gamma = 1   -> varimax
# gamma = P/2 -> equamax
#
lav_matrix_rotate_oblimin <- function(LAMBDA = NULL, oblimin.gamma = 0, ...,
                                      grad = FALSE) {
  nRow <- nrow(LAMBDA)
  nCol <- ncol(LAMBDA)
  ROW1 <- matrix(1.0, nCol, nCol)
  diag(ROW1) <- 0.0

  L2 <- LAMBDA * LAMBDA
  LR <- L2 %*% ROW1
  Jp <- matrix(1, nRow, nRow) / nRow

  # see Jennrich (2002, p. 11)
  tmp <- (diag(nRow) - oblimin.gamma * Jp) %*% LR

  # same as t( t(L2) - gamma * colMeans(L2) ) %*% ROW1

  out <- sum(L2 * tmp) / 4

  if (grad) {
    attr(out, "grad") <- LAMBDA * tmp
  }

  out
}


# quartimax criterion
# Carroll (1953); Saunders (1953) Neuhaus & Wrigley (1954); Ferguson (1954)
# we use here the equivalent 'Ferguson, 1954' variant
# (See Mulaik 2010, p. 303)
lav_matrix_rotate_quartimax <- function(LAMBDA = NULL, ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA
  out <- -1 * sum(L2 * L2) / 4

  if (grad) {
    attr(out, "grad") <- -1 * LAMBDA * L2
  }

  out
}


# varimax criterion
# Kaiser (1958, 1959)
#
# special case of the Orthomax family (Harman, 1960), where gamma = 1
# see Jennrich (2001, p. 296)
lav_matrix_rotate_varimax <- function(LAMBDA = NULL, ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA
  # center L2 column-wise
  cL2 <- t(t(L2) - colMeans(L2))
  out <- -1 * abs(sum(L2 * cL2)) / 4 # abs needed?

  if (grad) {
    attr(out, "grad") <- -1 * LAMBDA * cL2
  }

  out
}

# quartimin criterion (part of Carroll's oblimin family
lav_matrix_rotate_quartimin <- function(LAMBDA = NULL, ..., grad = FALSE) {
  nCol <- ncol(LAMBDA)
  ROW1 <- matrix(1.0, nCol, nCol)
  diag(ROW1) <- 0.0

  L2 <- LAMBDA * LAMBDA
  LR <- L2 %*% ROW1

  out <- sum(L2 * LR) / 4

  if (grad) {
    attr(out, "grad") <- LAMBDA * LR
  }

  out
}

# Browne's (2001) version of Yates (1984) geomin criterion
#
# we use the exp/log trick as in Bernaard & Jennrich (2005, p. 687)
lav_matrix_rotate_geomin <- function(LAMBDA = NULL, geomin.epsilon = 0.01,
                                     ..., grad = FALSE) {
  nCol <- ncol(LAMBDA)

  L2 <- LAMBDA * LAMBDA
  L2 <- L2 + geomin.epsilon

  if (geomin.epsilon < sqrt(.Machine$double.eps)) {
    # Yates's original formula
    tmp <- apply(L2, 1, prod)^(1 / nCol)
  } else {
    tmp <- exp(rowSums(log(L2)) / nCol)
  }

  out <- sum(tmp)

  if (grad) {
    attr(out, "grad") <- (2 / nCol) * LAMBDA / L2 * tmp
  }

  out
}

# simple entropy
# seems to only work for orthogonal rotation
lav_matrix_rotate_entropy <- function(LAMBDA = NULL, ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA

  # handle zero elements -> replace by '1', so log(1) == 0
  L2[L2 == 0] <- 1

  out <- -1 * sum(L2 * log(L2)) / 2

  if (grad) {
    attr(out, "grad") <- -LAMBDA * log(L2) - LAMBDA
  }

  out
}

# McCammon's (1966) Minimum Entropy Criterion
#
# for p-vector x, where x > 0 and sum(x) = 1, we have
# - entropy(x) == 0, if there is only one 1, and all zeroes
# - entropy(x) == max == log(p) if all elements are 1/p
# - entropy(x) is similar as complexity(x), but also measure of equality
#   of elements of x
#
# works only ok with orthogonal rotation!

lav_matrix_rotate_mccammon <- function(LAMBDA = NULL, ..., grad = FALSE) {
  nCol <- ncol(LAMBDA)
  nRow <- nrow(LAMBDA)
  L2 <- LAMBDA * LAMBDA

  # entropy function (Browne, 2001, eq 9)
  f_entropy <- function(x) {
    -1 * sum(ifelse(x > 0, x * log(x), 0))
  }

  # sums of rows/columns/all
  sumi. <- rowSums(L2)
  sum.j <- colSums(L2)
  sum.. <- sum(L2)


  Q1 <- f_entropy(t(L2) / sum.j) # encouraging columns with few large,
  # and many small elements
  Q2 <- f_entropy(sum.j / sum..) # encouraging equal column sums

  # minimize
  out <- log(Q1) - log(Q2)

  if (grad) { # See Bernaards and Jennrich 2005 page 685+686
    H <- -(log(t(t(L2) / sum.j)) + 1)
    G1 <- t(t(H) / sum.j - rowSums(t(L2 * H) / (sum.j * sum.j)))

    h <- -(log(sum.j / sum..) + 1)
    alpha <- as.numeric(h %*% sum.j) / (sum.. * sum..) # paper divides by
    # sum.., not sum..^2??
    G2 <- matrix(h / sum.. - alpha, nRow, nCol, byrow = TRUE)

    attr(out, "grad") <- 2 * LAMBDA * (G1 / Q1 - G2 / Q2)
  }

  out
}


# Infomax
# McKeon (1968, unpublished) and Browne (2001)
# Treat LAMBDA^2 as a contingency table, and use simplicity function based
# on tests for association; most effective was LRT for association
# (see Agresti, 1990, eq 3.13) which is maximized for max simplicity
#
# McKeon: criterion may be regarded as a measure of information about row
# categories conveyed by column categories (and vice versa); hence infomax

# - favors perfect cluster
# - discourages general factor
# - both for orthogonal and oblique rotation

#
# Note: typo in Browne (2001), see last  paragraph of Bernaards and
# Jennrich (2005) page 684
lav_matrix_rotate_infomax <- function(LAMBDA = NULL, ..., grad = FALSE) {
  nCol <- ncol(LAMBDA)
  nRow <- nrow(LAMBDA)
  L2 <- LAMBDA * LAMBDA

  # entropy function (Browne, 2001, eq 9)
  f_entropy <- function(x) {
    -1 * sum(ifelse(x > 0, x * log(x), 0))
  }

  # sums of rows/columns/all
  sumi. <- rowSums(L2)
  sum.j <- colSums(L2)
  sum.. <- sum(L2)

  Q1 <- f_entropy(L2 / sum..) # Bernaards & Jennrich version!! (Browne
  # divides by sum.j, like in McCammon)
  Q2 <- f_entropy(sum.j / sum..)
  Q3 <- f_entropy(sumi. / sum..)

  # minimize
  out <- log(nCol) + Q1 - Q2 - Q3

  if (grad) {
    H <- -(log(L2 / sum..) + 1)
    alpha <- sum(L2 * H) / (sum.. * sum..)
    G1 <- H / sum.. - alpha

    hj <- -(log(sum.j / sum..) + 1)
    alphaj <- as.numeric(hj %*% sum.j) / (sum.. * sum..)
    G2 <- matrix(hj, nRow, nCol, byrow = TRUE) / sum.. - alphaj

    hi <- -(log(sumi. / sum..) + 1)
    alphai <- as.numeric(sumi. %*% hi) / (sum.. * sum..)
    G3 <- matrix(hi, nRow, nCol) / sum.. - alphai

    attr(out, "grad") <- 2 * LAMBDA * (G1 - G2 - G3)
  }

  out
}

# oblimax
# Harman, 1976; Saunders, 1961
#
# for orthogonal rotation, oblimax is equivalent to quartimax
lav_matrix_rotate_oblimax <- function(LAMBDA = NULL, ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA

  # minimize version
  out <- -log(sum(L2 * L2)) + 2 * log(sum(L2))

  if (grad) {
    attr(out, "grad") <- (-4 * L2 * LAMBDA / (sum(L2 * L2))
      + 4 * LAMBDA / (sum(L2)))
  }

  out
}

# Bentler's Invariant Pattern Simplicity
# Bentler (1977)
#
#
lav_matrix_rotate_bentler <- function(LAMBDA = NULL, ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA

  L2tL2 <- crossprod(L2)
  L2tL2.inv <- lav_matrix_symmetric_inverse(S = L2tL2, logdet = TRUE)
  L2tL2.logdet <- attr(L2tL2.inv, "logdet")

  DIag <- diag(L2tL2)
  DIag.inv <- diag(1 / DIag)
  DIag.logdet <- sum(log(DIag)) # add small constant?

  # minimize version
  out <- -(L2tL2.logdet - DIag.logdet) / 4

  if (grad) {
    attr(out, "grad") <- -LAMBDA * (L2 %*% (L2tL2.inv - DIag.inv))
  }

  out
}


# The Tandem criteria
# Comrey (1967)
#
# only for sequential use:
# - tandem1 is used to determine the number of factors
#   (it removes the minor factors)
# - tandomII is used for final rotation
#
lav_matrix_rotate_tandem1 <- function(LAMBDA, ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA
  LL <- tcrossprod(LAMBDA)
  LL2 <- LL * LL

  # minimize version
  out <- -1 * sum(L2 * (LL2 %*% L2))

  if (grad) {
    tmp1 <- 4 * LAMBDA * (LL2 %*% L2)
    tmp2 <- 4 * (LL * (L2 %*% t(L2))) %*% LAMBDA
    attr(out, "grad") <- -tmp1 - tmp2
  }

  out
}

lav_matrix_rotate_tandem2 <- function(LAMBDA, ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA
  LL <- tcrossprod(LAMBDA)
  LL2 <- LL * LL

  # minimize version
  out <- sum(L2 * ((1 - LL2) %*% L2))

  if (grad) {
    tmp1 <- 4 * LAMBDA * ((1 - LL2) %*% L2)
    tmp2 <- 4 * (LL * tcrossprod(L2, L2)) %*% LAMBDA
    attr(out, "grad") <- tmp1 - tmp2
  }

  out
}





# simplimax
# Kiers (1994)
#
# oblique rotation method
# designed to rotate so that a given number 'k' of small loadings are
# as close to zero as possible
#
# may be viewed as partially specified target rotation with
# dynamically chosen weights
#
lav_matrix_rotate_simplimax <- function(LAMBDA = NULL, k = nrow(LAMBDA),
                                        ..., grad = FALSE) {
  L2 <- LAMBDA * LAMBDA

  # 'k' smallest element of L2
  small.element <- sort(L2)[k]

  # which elements are smaller than (or equal than) 'small.element'?
  ID <- sign(L2 <= small.element)

  # minimize version
  out <- sum(L2 * ID)

  if (grad) {
    attr(out, "grad") <- 2 * ID * LAMBDA
  }

  out
}


# target rotation
# Harman, 1976
#
# LAMBDA is rotated toward a specified target matrix 'target'
#
# Note: 'target' must be fully specified; if there are any NAs
#        use lav_matrix_rotate_pst() instead
#
lav_matrix_rotate_target <- function(LAMBDA = NULL, target = NULL,
                                     ..., grad = FALSE) {
  # squared difference
  DIFF <- LAMBDA - target
  DIFF2 <- DIFF * DIFF

  out <- sum(DIFF2, na.rm = TRUE)

  if (grad) {
    tmp <- 2 * DIFF
    # change NAs to zero
    tmp[is.na(tmp)] <- 0
    attr(out, "grad") <- tmp
  }

  out
}

# partially specified target rotation
#
# Browne 1972a, 1972b
#
# a pre-specified weight matrix W with ones/zeroes determines
# which elements of (LAMBDA - target) are used by the rotation criterion
#
# if 'target' contains NAs, they should correspond to '0' values in the
# target.mask matrix
#
lav_matrix_rotate_pst <- function(LAMBDA = NULL, target = NULL,
                                  target.mask = NULL, ..., grad = FALSE) {
  # mask target+LAMBDA
  target <- target.mask * target
  LAMBDA <- target.mask * LAMBDA

  # squared difference
  DIFF <- LAMBDA - target
  DIFF2 <- DIFF * DIFF

  # minimize
  out <- sum(DIFF2, na.rm = TRUE)

  if (grad) {
    tmp <- 2 * DIFF
    # change NAs to zero
    tmp[is.na(tmp)] <- 0
    attr(out, "grad") <- tmp
  }

  out
}


# bi-quartimin
#
# Jennrich & Bentler 2011
#
lav_matrix_rotate_biquartimin <- function(LAMBDA, ..., grad = FALSE) {
  # see Matlab code page 549
  stopifnot(ncol(LAMBDA) > 1L)

  # remove first column
  LAMBDA.group <- LAMBDA[, -1, drop = FALSE]

  # apply quartimin on the 'group' part
  out <- lav_matrix_rotate_quartimin(LAMBDA.group, ..., grad = grad)

  if (grad) {
    tmp <- attr(out, "grad")
    attr(out, "grad") <- cbind(0, tmp)
  }

  out
}


# bi-geomin
#
# Jennrich & Bentler 2012
#
lav_matrix_rotate_bigeomin <- function(LAMBDA, geomin.epsilon = 0.01, ...,
                                       grad = FALSE) {
  stopifnot(ncol(LAMBDA) > 1L)

  # remove first column
  LAMBDA.group <- LAMBDA[, -1, drop = FALSE]

  # apply geomin on the 'group' part
  out <- lav_matrix_rotate_geomin(LAMBDA.group,
    geomin.epsilon = geomin.epsilon, ...,
    grad = grad
  )

  if (grad) {
    tmp <- attr(out, "grad")
    attr(out, "grad") <- cbind(0, tmp)
  }

  out
}



# gradient check
ilav_matrix_rotate_grad_test <- function(crit = NULL, ...,
                                         LAMBDA = NULL,
                                         nRow = 20L, nCol = 5L,
                                         verbose = FALSE) {
  # test matrix
  if (is.null(LAMBDA)) {
    LAMBDA <- matrix(rnorm(nRow * nCol), nRow, nCol)
  }

  ff <- function(x, ...) {
    Lambda <- matrix(x, nRow, nCol)
    crit(Lambda, ..., grad = FALSE)
  }

  GQ1 <- matrix(
    numDeriv::grad(func = ff, x = as.vector(LAMBDA), ...),
    nRow, nCol
  )
  GQ2 <- attr(crit(LAMBDA, ..., grad = TRUE), "grad")

  if (verbose) {
    print(list(LAMBDA = LAMBDA, GQ1 = GQ1, GQ2 = GQ2))
  }

  all.equal(GQ1, GQ2, tolerance = 1e-07)
}

ilav_matrix_rotate_grad_test_all <- function() {
  # Orthomax family with various values for gamma
  for (gamma in seq(0, 1, 0.2)) {
    check <- ilav_matrix_rotate_grad_test(
      crit = lav_matrix_rotate_orthomax,
      gamma = gamma
    )
    if (is.logical(check) && check) {
      cat(
        "orthomax + gamma = ", sprintf("%3.1f", gamma),
        ": OK\n"
      )
    } else {
      cat(
        "orthomax + gamma = ", sprintf("%3.1f", gamma),
        ": FAILED\n"
      )
    }
  }

  # Crawford-Ferguson with various values for gamma
  for (gamma in seq(0, 1, 0.2)) {
    check <- ilav_matrix_rotate_grad_test(
      crit = lav_matrix_rotate_cf,
      gamma = gamma
    )
    if (is.logical(check) && check) {
      cat(
        "Crawford-Ferguson + gamma = ", sprintf("%3.1f", gamma),
        ": OK\n"
      )
    } else {
      cat(
        "Crawford-Ferguson + gamma = ", sprintf("%3.1f", gamma),
        ": FAILED\n"
      )
    }
  }

  # Oblimin family with various values for gamma
  for (gamma in seq(0, 1, 0.2)) {
    check <- ilav_matrix_rotate_grad_test(
      crit = lav_matrix_rotate_oblimin,
      gamma = gamma
    )
    if (is.logical(check) && check) {
      cat(
        "Oblimin + gamma = ", sprintf("%3.1f", gamma),
        ": OK\n"
      )
    } else {
      cat(
        "Oblimin + gamma = ", sprintf("%3.1f", gamma),
        ": FAILED\n"
      )
    }
  }

  # quartimax
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_quartimax)
  if (is.logical(check) && check) {
    cat("quartimax: OK\n")
  } else {
    cat("quartimax: FAILED\n")
  }

  # varimax
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_varimax)
  if (is.logical(check) && check) {
    cat("varimax: OK\n")
  } else {
    cat("varimax: FAILED\n")
  }

  # quartimin
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_quartimin)
  if (is.logical(check) && check) {
    cat("quartimin: OK\n")
  } else {
    cat("quartimin: FAILED\n")
  }

  # geomin
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_geomin)
  if (is.logical(check) && check) {
    cat("geomin: OK\n")
  } else {
    cat("geomin: FAILED\n")
  }

  # simple entropy
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_entropy)
  if (is.logical(check) && check) {
    cat("entropy: OK\n")
  } else {
    cat("entropy: FAILED\n")
  }

  # McCammon entropy criterion
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_mccammon)
  if (is.logical(check) && check) {
    cat("McCammon: OK\n")
  } else {
    cat("McCammon: FAILED\n")
  }

  # infomax
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_infomax)
  if (is.logical(check) && check) {
    cat("infomax: OK\n")
  } else {
    cat("infomax: FAILED\n")
  }

  # oblimax
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_oblimax)
  if (is.logical(check) && check) {
    cat("oblimax: OK\n")
  } else {
    cat("oblimax: FAILED\n")
  }

  # bentler
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_bentler)
  if (is.logical(check) && check) {
    cat("bentler: OK\n")
  } else {
    cat("bentler: FAILED\n")
  }

  # simplimax
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_simplimax)
  if (is.logical(check) && check) {
    cat("simplimax: OK\n")
  } else {
    cat("simplimax: FAILED\n")
  }

  # tandem1
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_tandem1)
  if (is.logical(check) && check) {
    cat("tandem1: OK\n")
  } else {
    cat("tandem1: FAILED\n")
  }

  # tandem2
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_tandem2)
  if (is.logical(check) && check) {
    cat("tandem2: OK\n")
  } else {
    cat("tandem2: FAILED\n")
  }

  # bi-quartimin
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_biquartimin)
  if (is.logical(check) && check) {
    cat("biquartimin: OK\n")
  } else {
    cat("biquartimin: FAILED\n")
  }

  # bi-quartimin
  check <- ilav_matrix_rotate_grad_test(crit = lav_matrix_rotate_bigeomin)
  if (is.logical(check) && check) {
    cat("bigeomin: OK\n")
  } else {
    cat("bigeomin: FAILED\n")
  }
}
