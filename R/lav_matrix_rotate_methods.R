# various rotation criteria and their gradients
# YR 05 April 2019: initial version
# YR 14 June  2019: add more rotation criteria
# YR 22 Dec   2025: add multigroup rotation criteria

# references:
#
# Bernaards, C. A., & Jennrich, R. I. (2005). Gradient projection algorithms
# and software for arbitrary rotation criteria in factor analysis. Educational
# and Psychological Measurement, 65(5), 676-696.
# old website: http://web.archive.org/web/20180708170331/http://www.stat.ucla.edu/research/gpa/splusfunctions.net # nolint
#
# Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
# analysis. Multivariate behavioral research, 36(1), 111-150.
#
# Mulaik, S. A. (2010). Foundations of factor analysis (Second Edition).
# Boca Raton: Chapman and Hall/CRC.
#
# De Roover, K., & Vermunt, J. K. (2019). On the exploratory road to unraveling
# factor loading non-invariance: A new multigroup rotation approach. Structural
# Equation Modeling: A Multidisciplinary Journal, 26(6), 905-923.

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
lav_mat_rotate_orthomax <- function(mm_lambda = NULL, orthomax_gamma = 1,
                                       ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda
  # center L2 column-wise
  c_l2 <- t(t(l2) - orthomax_gamma * colMeans(l2))
  out <- -1 * sum(l2 * c_l2) / 4

  if (grad) {
    attr(out, "grad") <- -1 * mm_lambda * c_l2
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
lav_mat_rotate_cf <- function(mm_lambda = NULL, cf_gamma = 0, ...,
                                 grad = FALSE) {
  # check if gamma is between 0 and 1?
  n_row <- nrow(mm_lambda)
  n_col <- ncol(mm_lambda)
  row1 <- matrix(1.0, n_col, n_col)
  diag(row1) <- 0.0
  col1 <- matrix(1.0, n_row, n_row)
  diag(col1) <- 0.0

  l2 <- mm_lambda * mm_lambda
  lr <- l2 %*% row1
  lc <- col1 %*% l2

  f1 <- sum(l2 * lr) / 4
  f2 <- sum(l2 * lc) / 4

  out <- (1 - cf_gamma) * f1 + cf_gamma * f2

  if (grad) {
    attr(out, "grad") <- ((1 - cf_gamma) * mm_lambda * lr) +
      (cf_gamma * mm_lambda * lc)
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
lav_mat_rotate_oblimin <- function(mm_lambda = NULL, oblimin_gamma = 0, ...,
                                      grad = FALSE) {
  n_row <- nrow(mm_lambda)
  n_col <- ncol(mm_lambda)
  row1 <- matrix(1.0, n_col, n_col)
  diag(row1) <- 0.0

  l2 <- mm_lambda * mm_lambda
  lr <- l2 %*% row1
  jp <- matrix(1, n_row, n_row) / n_row

  # see Jennrich (2002, p. 11)
  tmp <- (diag(n_row) - oblimin_gamma * jp) %*% lr

  # same as t( t(L2) - gamma * colMeans(L2) ) %*% ROW1

  out <- sum(l2 * tmp) / 4

  if (grad) {
    attr(out, "grad") <- mm_lambda * tmp
  }

  out
}


# quartimax criterion
# Carroll (1953); Saunders (1953) Neuhaus & Wrigley (1954); Ferguson (1954)
# we use here the equivalent 'Ferguson, 1954' variant
# (See Mulaik 2010, p. 303)
lav_mat_rotate_quartimax <- function(mm_lambda = NULL, ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda
  out <- -1 * sum(l2 * l2) / 4

  if (grad) {
    attr(out, "grad") <- -1 * mm_lambda * l2
  }

  out
}


# varimax criterion
# Kaiser (1958, 1959)
#
# special case of the Orthomax family (Harman, 1960), where gamma = 1
# see Jennrich (2001, p. 296)
lav_mat_rotate_varimax <- function(mm_lambda = NULL, ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda
  # center L2 column-wise
  c_l2 <- t(t(l2) - colMeans(l2))
  out <- -1 * abs(sum(l2 * c_l2)) / 4 # abs needed?

  if (grad) {
    attr(out, "grad") <- -1 * mm_lambda * c_l2
  }

  out
}

# quartimin criterion (part of Carroll's oblimin family
lav_mat_rotate_quartimin <- function(mm_lambda = NULL, ..., grad = FALSE) {
  n_col <- ncol(mm_lambda)
  row1 <- matrix(1.0, n_col, n_col)
  diag(row1) <- 0.0

  l2 <- mm_lambda * mm_lambda
  lr <- l2 %*% row1

  out <- sum(l2 * lr) / 4

  if (grad) {
    attr(out, "grad") <- mm_lambda * lr
  }

  out
}

# Browne's (2001) version of Yates (1984) geomin criterion
#
# we use the exp/log trick as in Bernaard & Jennrich (2005, p. 687)
lav_mat_rotate_geomin <- function(mm_lambda = NULL, geomin_epsilon = 0.01,
                                     ..., grad = FALSE) {
  n_col <- ncol(mm_lambda)

  l2 <- mm_lambda * mm_lambda
  l2 <- l2 + geomin_epsilon

  if (geomin_epsilon < sqrt(.Machine$double.eps)) {
    # Yates's original formula
    tmp <- apply(l2, 1, prod)^(1 / n_col)
  } else {
    tmp <- exp(rowSums(log(l2)) / n_col)
  }

  out <- sum(tmp)

  if (grad) {
    attr(out, "grad") <- (2 / n_col) * mm_lambda / l2 * tmp
  }

  out
}

# simple entropy
# seems to only work for orthogonal rotation
lav_mat_rotate_entropy <- function(mm_lambda = NULL, ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda

  # handle zero elements -> replace by '1', so log(1) == 0
  l2[l2 == 0] <- 1

  out <- -1 * sum(l2 * log(l2)) / 2

  if (grad) {
    attr(out, "grad") <- -mm_lambda * log(l2) - mm_lambda
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

lav_mat_rotate_mccammon <- function(mm_lambda = NULL, ..., grad = FALSE) {
  n_col <- ncol(mm_lambda)
  n_row <- nrow(mm_lambda)
  l2 <- mm_lambda * mm_lambda

  # entropy function (Browne, 2001, eq 9)
  f_entropy <- function(x) {
    -1 * sum(ifelse(x > 0, x * log(x), 0))
  }

  # sums of rows/columns/all
  sum_j <- colSums(l2)
  sum_1 <- sum(l2)


  q1 <- f_entropy(t(l2) / sum_j) # encouraging columns with few large,
  # and many small elements
  q2 <- f_entropy(sum_j / sum_1) # encouraging equal column sums

  # minimize
  out <- log(q1) - log(q2)

  if (grad) { # See Bernaards and Jennrich 2005 page 685+686
    h_1 <- -(log(t(t(l2) / sum_j)) + 1)
    g1 <- t(t(h_1) / sum_j - rowSums(t(l2 * h_1) / (sum_j * sum_j)))

    h <- -(log(sum_j / sum_1) + 1)
    alpha <- as.numeric(h %*% sum_j) / (sum_1 * sum_1) # paper divides by
    # sum.., not sum..^2??
    g2 <- matrix(h / sum_1 - alpha, n_row, n_col, byrow = TRUE)

    attr(out, "grad") <- 2 * mm_lambda * (g1 / q1 - g2 / q2)
  }

  out
}


# Infomax
# McKeon (1968, unpublished) and Browne (2001)
# Treat mm_lambda^2 as a contingency table, and use simplicity function based
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
lav_mat_rotate_infomax <- function(mm_lambda = NULL, ..., grad = FALSE) {
  n_col <- ncol(mm_lambda)
  n_row <- nrow(mm_lambda)
  l2 <- mm_lambda * mm_lambda

  # entropy function (Browne, 2001, eq 9)
  f_entropy <- function(x) {
    -1 * sum(ifelse(x > 0, x * log(x), 0))
  }

  # sums of rows/columns/all
  sumi <- rowSums(l2)
  sum_j <- colSums(l2)
  sum_1 <- sum(l2)

  q1 <- f_entropy(l2 / sum_1) # Bernaards & Jennrich version!! (Browne
  # divides by sum.j, like in McCammon)
  q2 <- f_entropy(sum_j / sum_1)
  q3 <- f_entropy(sumi / sum_1)

  # minimize
  out <- log(n_col) + q1 - q2 - q3

  if (grad) {
    h <- -(log(l2 / sum_1) + 1)
    alpha <- sum(l2 * h) / (sum_1 * sum_1)
    g1 <- h / sum_1 - alpha

    hj <- -(log(sum_j / sum_1) + 1)
    alphaj <- as.numeric(hj %*% sum_j) / (sum_1 * sum_1)
    g2 <- matrix(hj, n_row, n_col, byrow = TRUE) / sum_1 - alphaj

    hi <- -(log(sumi / sum_1) + 1)
    alphai <- as.numeric(sumi %*% hi) / (sum_1 * sum_1)
    g3 <- matrix(hi, n_row, n_col) / sum_1 - alphai

    attr(out, "grad") <- 2 * mm_lambda * (g1 - g2 - g3)
  }

  out
}

# oblimax
# Harman, 1976; Saunders, 1961
#
# for orthogonal rotation, oblimax is equivalent to quartimax
lav_mat_rotate_oblimax <- function(mm_lambda = NULL, ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda

  # minimize version
  out <- -log(sum(l2 * l2)) + 2 * log(sum(l2))

  if (grad) {
    attr(out, "grad") <- (-4 * l2 * mm_lambda / (sum(l2 * l2))
      + 4 * mm_lambda / (sum(l2)))
  }

  out
}

# Bentler's Invariant Pattern Simplicity
# Bentler (1977)
#
#
lav_mat_rotate_bentler <- function(mm_lambda = NULL, ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda

  l2t_l2 <- crossprod(l2)
  l2t_l2_inv <- lav_mat_sym_inverse(s = l2t_l2, logdet = TRUE)
  l2t_l2_logdet <- attr(l2t_l2_inv, "logdet")

  diag_1 <- diag(l2t_l2)
  diag_inv <- diag(1 / diag_1)
  diag_logdet <- sum(log(diag_1)) # add small constant?

  # minimize version
  out <- -(l2t_l2_logdet - diag_logdet) / 4

  if (grad) {
    attr(out, "grad") <- -mm_lambda * (l2 %*% (l2t_l2_inv - diag_inv))
  }

  out
}


# The Tandem criteria
# Comrey (1967)
#
# only for sequential use:
# - tandem1 is used to determine the number of factors
#   (it removes the minor factors)
# - tandem2 is used for final rotation
#
lav_mat_rotate_tandem1 <- function(mm_lambda, ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda
  ll <- tcrossprod(mm_lambda)
  ll2 <- ll * ll

  # minimize version
  out <- -1 * sum(l2 * (ll2 %*% l2))

  if (grad) {
    tmp1 <- 4 * mm_lambda * (ll2 %*% l2)
    tmp2 <- 4 * (ll * (l2 %*% t(l2))) %*% mm_lambda
    attr(out, "grad") <- -tmp1 - tmp2
  }

  out
}

lav_mat_rotate_tandem2 <- function(mm_lambda, ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda
  ll <- tcrossprod(mm_lambda)
  ll2 <- ll * ll

  # minimize version
  out <- sum(l2 * ((1 - ll2) %*% l2))

  if (grad) {
    tmp1 <- 4 * mm_lambda * ((1 - ll2) %*% l2)
    tmp2 <- 4 * (ll * tcrossprod(l2, l2)) %*% mm_lambda
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
lav_mat_rotate_simplimax <- function(mm_lambda = NULL, k = nrow(mm_lambda),
                                        ..., grad = FALSE) {
  l2 <- mm_lambda * mm_lambda

  # 'k' smallest element of L2
  small_element <- sort(l2)[k]

  # which elements are smaller than (or equal than) 'small.element'?
  id <- sign(l2 <= small_element)

  # minimize version
  out <- sum(l2 * id)

  if (grad) {
    attr(out, "grad") <- 2 * id * mm_lambda
  }

  out
}


# target rotation
# Harman, 1976
#
# mm_lambda is rotated toward a specified target matrix 'target'
#
# Note: 'target' must be fully specified; if there are any NAs
#        use lav_mat_rotate_pst() instead
#
lav_mat_rotate_target <- function(mm_lambda = NULL, target = NULL,
                                     ..., grad = FALSE) {
  # squared difference
  diff_1 <- mm_lambda - target
  diff2 <- diff_1 * diff_1

  out <- sum(diff2, na.rm = TRUE)

  if (grad) {
    tmp <- 2 * diff_1
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
# which elements of (mm_lambda - target) are used by the rotation criterion
#
# if 'target' contains NAs, they should correspond to '0' values in the
# target_mask matrix
#
lav_mat_rotate_pst <- function(mm_lambda = NULL, target = NULL,
                                  target_mask = NULL, ..., grad = FALSE) {
  # mask target+mm_lambda
  target <- target_mask * target
  mm_lambda <- target_mask * mm_lambda

  # squared difference
  diff_1 <- mm_lambda - target
  diff2 <- diff_1 * diff_1

  # minimize
  out <- sum(diff2, na.rm = TRUE)

  if (grad) {
    tmp <- 2 * diff_1
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
lav_mat_rotate_biquartimin <- function(mm_lambda, ..., grad = FALSE) {
  # see Matlab code page 549
  stopifnot(ncol(mm_lambda) > 1L)

  # remove first column
  lambda_group <- mm_lambda[, -1, drop = FALSE]

  # apply quartimin on the 'group' part
  out <- lav_mat_rotate_quartimin(lambda_group, ..., grad = grad)

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
lav_mat_rotate_bigeomin <- function(mm_lambda, geomin_epsilon = 0.01, ...,
                                       grad = FALSE) {
  stopifnot(ncol(mm_lambda) > 1L)

  # remove first column
  lambda_group <- mm_lambda[, -1, drop = FALSE]

  # apply geomin on the 'group' part
  out <- lav_mat_rotate_geomin(lambda_group,
    geomin_epsilon = geomin_epsilon, ...,
    grad = grad
  )

  if (grad) {
    tmp <- attr(out, "grad")
    attr(out, "grad") <- cbind(0, tmp)
  }

  out
}

# multiple group rotation + agreement
lav_mat_rotate_mg_agreement <- function(lambda_list, method_fname = "geomin",
                                   method_args = list(),
                                   mg_agreement_method = "procrustes",
                                   w = 0.5, scale = TRUE) {
  ngroups <- length(lambda_list)

  # 0. rescale, so all columns have the same (average) sum-of-squares
  if (scale) {
    ss_list <- lapply(lambda_list, function(x) colSums(x * x))
    ss_ave <- Reduce("+", ss_list) / length(ss_list)
    for (g in seq_len(ngroups)) {
      lambda_list[[g]] <- t(t(lambda_list[[g]]) /
        sqrt(ss_list[[g]]) * sqrt(ss_ave))
    }
  }

  # 1. simple structure crit for each group
  q_group <- numeric(ngroups)
  for (g in seq_len(ngroups)) {
    q_group[g] <- do.call(method_fname, c(
      list(mm_lambda = lambda_list[[g]]),
      method_args, list(grad = FALSE)
    ))
  }
  # FIXME: should we 'weight' the groups, based on sample size?
  q_mg <- sum(q_group)

  # 2. agreement (generalized procrustes) across all groups
  if (mg_agreement_method == "procrustes") {
    a_mg <- 0
    pairs <- utils::combn(seq_len(ngroups), 2L)
    for (p in seq_len(ncol(pairs))) {
      g1 <- pairs[1L, p]
      g2 <- pairs[2L, p]
      diff <- (lambda_list[[g1]] - lambda_list[[g2]])
      diff2 <- diff * diff
      # FIXME: should we 'weight' the groups, based on sample size?
      a_mg <- a_mg + sum(diff2)
    }
  } else {
    lav_msg_stop(gettext(
      "only mg.agreement.method procrustes is supported for now"))
  }

  out <- (w * a_mg) + ((1 - w) * q_mg)

  attr(out, "Q.mg") <- q_mg
  attr(out, "A.mg") <- a_mg

  out
}

# gradient check
ilav_mat_rotate_grad_test <- function(crit = NULL, ...,
                                         mm_lambda = NULL,
                                         n_row = 20L, n_col = 5L) {
  # test matrix
  if (is.null(mm_lambda)) {
    mm_lambda <- matrix(rnorm(n_row * n_col), n_row, n_col)
  }

  ff <- function(x, ...) {
    lambda <- matrix(x, n_row, n_col)
    crit(lambda, ..., grad = FALSE)
  }

  gq1 <- matrix(
    numDeriv::grad(func = ff, x = as.vector(mm_lambda), ...),
    n_row, n_col
  )
  gq2 <- attr(crit(mm_lambda, ..., grad = TRUE), "grad")

  if (lav_verbose()) {
    print(list(mm_lambda = mm_lambda, GQ1 = gq1, GQ2 = gq2))
  }

  all.equal(gq1, gq2, tolerance = 1e-07)
}

ilav_mat_rotate_grad_test_all <- function() {
  # Orthomax family with various values for gamma
  for (gamma in seq(0, 1, 0.2)) {
    check <- ilav_mat_rotate_grad_test(
      crit = lav_mat_rotate_orthomax,
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
    check <- ilav_mat_rotate_grad_test(
      crit = lav_mat_rotate_cf,
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
    check <- ilav_mat_rotate_grad_test(
      crit = lav_mat_rotate_oblimin,
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
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_quartimax)
  if (is.logical(check) && check) {
    cat("quartimax: OK\n")
  } else {
    cat("quartimax: FAILED\n")
  }

  # varimax
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_varimax)
  if (is.logical(check) && check) {
    cat("varimax: OK\n")
  } else {
    cat("varimax: FAILED\n")
  }

  # quartimin
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_quartimin)
  if (is.logical(check) && check) {
    cat("quartimin: OK\n")
  } else {
    cat("quartimin: FAILED\n")
  }

  # geomin
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_geomin)
  if (is.logical(check) && check) {
    cat("geomin: OK\n")
  } else {
    cat("geomin: FAILED\n")
  }

  # simple entropy
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_entropy)
  if (is.logical(check) && check) {
    cat("entropy: OK\n")
  } else {
    cat("entropy: FAILED\n")
  }

  # McCammon entropy criterion
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_mccammon)
  if (is.logical(check) && check) {
    cat("McCammon: OK\n")
  } else {
    cat("McCammon: FAILED\n")
  }

  # infomax
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_infomax)
  if (is.logical(check) && check) {
    cat("infomax: OK\n")
  } else {
    cat("infomax: FAILED\n")
  }

  # oblimax
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_oblimax)
  if (is.logical(check) && check) {
    cat("oblimax: OK\n")
  } else {
    cat("oblimax: FAILED\n")
  }

  # bentler
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_bentler)
  if (is.logical(check) && check) {
    cat("bentler: OK\n")
  } else {
    cat("bentler: FAILED\n")
  }

  # simplimax
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_simplimax)
  if (is.logical(check) && check) {
    cat("simplimax: OK\n")
  } else {
    cat("simplimax: FAILED\n")
  }

  # tandem1
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_tandem1)
  if (is.logical(check) && check) {
    cat("tandem1: OK\n")
  } else {
    cat("tandem1: FAILED\n")
  }

  # tandem2
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_tandem2)
  if (is.logical(check) && check) {
    cat("tandem2: OK\n")
  } else {
    cat("tandem2: FAILED\n")
  }

  # bi-quartimin
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_biquartimin)
  if (is.logical(check) && check) {
    cat("biquartimin: OK\n")
  } else {
    cat("biquartimin: FAILED\n")
  }

  # bi-quartimin
  check <- ilav_mat_rotate_grad_test(crit = lav_mat_rotate_bigeomin)
  if (is.logical(check) && check) {
    cat("bigeomin: OK\n")
  } else {
    cat("bigeomin: FAILED\n")
  }
}
