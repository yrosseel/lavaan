# Magnus & Neudecker (1999) style matrix operations
# YR - 11 may 2011: initial version
# YR - 19 okt 2014: rename functions using lav_matrix_ prefix

# vec operator
#
# the vec operator (for 'vectorization') transforms a matrix into
# a vector by stacking the *columns* of the matrix one underneath the other
#
# M&N book: page 30
#
# note: we do not coerce to 'double/numeric' storage-mode (like as.numeric)
lav_matrix_vec <- function(A) {
  as.vector(A)
}


# vecr operator
#
# the vecr operator ransforms a matrix into
# a vector by stacking the *rows* of the matrix one underneath the other
lav_matrix_vecr <- function(A) {
  # faster way??
  # nRow <- NROW(A); nCol <- NCOL(A)
  # idx <- (seq_len(nCol) - 1L) * nRow + rep(seq_len(nRow), each = nCol)

  lav_matrix_vec(t(A))
}


# vech
#
# the vech operator (for 'half vectorization') transforms a *symmetric* matrix
# into a vector by stacking the *columns* of the matrix one underneath the
# other, but eliminating all supradiagonal elements
#
# see Henderson & Searle, 1979
#
# M&N book: page 48-49
#
lav_matrix_vech <- function(S, diagonal = TRUE) {
  ROW <- row(S)
  COL <- col(S)
  if (diagonal) S[ROW >= COL] else S[ROW > COL]
}


# the vechr operator transforms a *symmetric* matrix
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all supradiagonal elements
lav_matrix_vechr <- function(S, diagonal = TRUE) {
  S[lav_matrix_vechr_idx(n = NCOL(S), diagonal = diagonal)]
}


# the vechu operator transforms a *symmetric* matrix
# into a vector by stacking the *columns* of the matrix one after the
# other, but eliminating all infradiagonal elements
lav_matrix_vechu <- function(S, diagonal = TRUE) {
  S[lav_matrix_vechu_idx(n = NCOL(S), diagonal = diagonal)]
}


# the vechru operator transforms a *symmetric* matrix
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all infradiagonal elements
#
# same as vech (but using upper-diagonal elements)
lav_matrix_vechru <- function(S, diagonal = TRUE) {
  S[lav_matrix_vechru_idx(n = NCOL(S), diagonal = diagonal)]
}




# return the *vector* indices of the lower triangular elements of a
# symmetric matrix of size 'n'
lav_matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (n < 100L) {
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if (diagonal) which(ROW >= COL) else which(ROW > COL)
  } else {
    # ldw version
    if (diagonal) {
      unlist(lapply(
        seq_len(n),
        function(j) (j - 1L) * n + seq.int(j, n)
      ))
    } else {
      unlist(lapply(
        seq_len(n - 1L),
        function(j) (j - 1L) * n + seq.int(j + 1L, n)
      ))
    }
  }
}

# return the *row* indices of the lower triangular elements of a
# symmetric matrix of size 'n'
lav_matrix_vech_row_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (diagonal) {
    unlist(lapply(seq_len(n), seq.int, n))
  } else {
    1 + unlist(lapply(seq_len(n - 1), seq.int, n - 1))
  }
}

# return the *col* indices of the lower triangular elements of a
# symmetric matrix of size 'n'
lav_matrix_vech_col_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (!diagonal) {
    n <- n - 1L
  }
  rep.int(seq_len(n), times = rev(seq_len(n)))
}




# return the *vector* indices of the lower triangular elements of a
# symmetric matrix of size 'n' -- ROW-WISE
lav_matrix_vechr_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (n < 100L) {
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
    if (diagonal) tmp[ROW <= COL] else tmp[ROW < COL]
  } else {
    if (diagonal) {
      unlist(lapply(
        seq_len(n),
        function(j) seq.int(1, j) * n - (n - j)
      ))
    } else {
      unlist(lapply(
        seq_len(n - 1L),
        function(j) seq.int(1, j) * n - (n - j) + 1
      ))
    }
  }
}

# return the *vector* indices of the upper triangular elements of a
# symmetric matrix of size 'n' -- COLUMN-WISE
lav_matrix_vechu_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (n < 100L) {
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    if (diagonal) which(ROW <= COL) else which(ROW < COL)
  } else {
    if (diagonal) {
      unlist(lapply(seq_len(n), function(j) seq.int(j) + (j - 1) * n))
    } else {
      unlist(lapply(seq_len(n - 1L), function(j) seq.int(j) + j * n))
    }
  }
}

# return the *vector* indices of the upper triangular elements of a
# symmetric matrix of size 'n' -- ROW-WISE
lav_matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  if (n < 100L) {
    # FIXME!! make this more efficient (without creating 3 n*n matrices!)
    ROW <- matrix(seq_len(n), n, n)
    COL <- matrix(seq_len(n), n, n, byrow = TRUE)
    tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
    if (diagonal) tmp[ROW >= COL] else tmp[ROW > COL]
  } else {
    # ldw version
    if (diagonal) {
      unlist(lapply(
        seq_len(n),
        function(j) seq.int(j - 1, n - 1) * n + j
      ))
    } else {
      unlist(lapply(
        seq_len(n - 1L),
        function(j) seq.int(j, n - 1) * n + j
      ))
    }
  }
}


# vech.reverse and vechru.reverse (aka `upper2full')
#
# given the output of vech(S) --or vechru(S) which is identical--
# reconstruct S
lav_matrix_vech_reverse <- lav_matrix_vechru_reverse <-
  lav_matrix_upper2full <-
  function(x, diagonal = TRUE) {
    # guess dimensions
    if (diagonal) {
      p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    } else {
      p <- (sqrt(1 + 8 * length(x)) + 1) / 2
    }

    S <- numeric(p * p)
    S[lav_matrix_vech_idx(p, diagonal = diagonal)] <- x
    S[lav_matrix_vechru_idx(p, diagonal = diagonal)] <- x

    attr(S, "dim") <- c(p, p)
    S
  }


# vechr.reverse vechu.reversie (aka `lower2full')
#
# given the output of vechr(S) --or vechu(S) which is identical--
# reconstruct S
lav_matrix_vechr_reverse <- lav_matrix_vechu_reverse <-
  lav_matrix_lower2full <- function(x, diagonal = TRUE) {
    # guess dimensions
    if (diagonal) {
      p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    } else {
      p <- (sqrt(1 + 8 * length(x)) + 1) / 2
    }
    stopifnot(p == round(p, 0))

    S <- numeric(p * p)
    S[lav_matrix_vechr_idx(p, diagonal = diagonal)] <- x
    S[lav_matrix_vechu_idx(p, diagonal = diagonal)] <- x

    attr(S, "dim") <- c(p, p)
    S
  }


# return the *vector* indices of the diagonal elements of a symmetric
# matrix of size 'n'
lav_matrix_diag_idx <- function(n = 1L) {
  # if(n < 1L) return(integer(0L))
  1L + (seq_len(n) - 1L) * (n + 1L)
}


# return the *vector* indices of the diagonal elements of the LOWER part
# of a symmatrix matrix of size 'n'
lav_matrix_diagh_idx <- function(n = 1L) {
  if (n < 1L) {
    return(integer(0L))
  }
  if (n == 1L) {
    return(1L)
  }
  c(1L, cumsum(n:2L) + 1L)
}


# return the *vector* indices of the ANTI diagonal elements of a symmetric
# matrix of size 'n'
lav_matrix_antidiag_idx <- function(n = 1L) {
  if (n < 1L) {
    return(integer(0L))
  }
  1L + seq_len(n) * (n - 1L)
}

# return the *vector* indices of 'idx' elements in a vech() matrix
#
# eg if n = 4 and type == "and" and idx = c(2,4)
#    we create matrix A =
#       [,1]  [,2]  [,3]  [,4]
# [1,] FALSE FALSE FALSE FALSE
# [2,] FALSE  TRUE FALSE  TRUE
# [3,] FALSE FALSE FALSE FALSE
# [4,] FALSE  TRUE FALSE  TRUE
#
# and the result is c(5,7,10)
#
# eg if n = 4 and type == "or" and idx = c(2,4)
#    we create matrix A =
#       [,1] [,2]  [,3] [,4]
# [1,] FALSE TRUE FALSE TRUE
# [2,]  TRUE TRUE  TRUE TRUE
# [3,] FALSE TRUE FALSE TRUE
# [4,]  TRUE TRUE  TRUE TRUE
#
# and the result is c(2, 4, 5, 6, 7, 9, 10)
#
lav_matrix_vech_which_idx <- function(n = 1L, diagonal = TRUE,
                                      idx = integer(0L), type = "and",
                                      add.idx.at.start = FALSE) {
  if (length(idx) == 0L) {
    return(integer(0L))
  }
  n <- as.integer(n)
  A <- matrix(FALSE, n, n)
  if (type == "and") {
    A[idx, idx] <- TRUE
  } else if (type == "or") {
    A[idx, ] <- TRUE
    A[, idx] <- TRUE
  }
  pstar.idx <- which(lav_matrix_vech(A, diagonal = diagonal))

  if (add.idx.at.start) {
    pstar.idx <- c(idx, pstar.idx + n)
  }

  pstar.idx
}

# similar to lav_matrix_vech_which_idx(), but
# - only 'type = and'
# - order of idx matters!
lav_matrix_vech_match_idx <- function(n = 1L, diagonal = TRUE,
                                      idx = integer(0L)) {
  if (length(idx) == 0L) {
    return(integer(0L))
  }
  n <- as.integer(n)
  pstar <- n * (n + 1) / 2
  A <- lav_matrix_vech_reverse(seq_len(pstar))
  B <- A[idx, idx, drop = FALSE]
  lav_matrix_vech(B, diagonal = diagonal)
}

# check if square matrix is diagonal (no tolerance!)
lav_matrix_is_diagonal <- function(A = NULL) {
  A <- as.matrix.default(A)
  stopifnot(nrow(A) == ncol(A))

  diag(A) <- 0
  all(A == 0)
}


# create the duplication matrix (D_n): it 'duplicates' the elements
# in vech(S) to create vec(S) (where S is symmetric)
#
# D %*% vech(S) == vec(S)
#
# M&N book: pages 48-50
#
# note: several flavors: dup1, dup2, dup3, ...
# currently used: dup3

# first attempt
# dup1: working on the vector indices only
.dup1 <- function(n = 1L) {
  if ((n < 1L) | (round(n) != n)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }
  if (n > 255L) {
    lav_msg_stop(gettext("n is too large"))
  }

  # dimensions
  n2 <- n * n
  nstar <- n * (n + 1) / 2
  # THIS is the real bottleneck: allocating an ocean of zeroes...
  x <- numeric(n2 * nstar)

  # delta patterns
  r1 <- seq.int(from = n * n + 1, by = -(n - 1), length.out = n - 1)
  r2 <- seq.int(from = n - 1, by = n - 1, length.out = n - 1)
  r3 <- seq.int(from = 2 * n + 1, by = n, length.out = n - 1)

  # is there a more elegant way to do this?
  rr <- unlist(lapply(
    (n - 1):1,
    function(x) {
      c(rbind(r1[1:x], r2[1:x]), r3[n - x])
    }
  ))

  idx <- c(1L, cumsum(rr) + 1L)

  # create matrix
  x[idx] <- 1.0
  attr(x, "dim") <- c(n2, nstar)

  x
}

# second attempt
# dup2: working on the row/col matrix indices
# (but only create the matrix at the very end)
.dup2 <- function(n = 1L) {
  if ((n < 1L) | (round(n) != n)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }

  if (n > 255L) {
    lav_msg_stop(gettext("n is too large"))
  }

  nstar <- n * (n + 1) / 2
  n2 <- n * n
  # THIS is the real bottleneck: allocating an ocean of zeroes...
  x <- numeric(n2 * nstar)

  idx1 <- lav_matrix_vech_idx(n) + ((1L:nstar) - 1L) * n2 # vector indices
  idx2 <- lav_matrix_vechru_idx(n) + ((1L:nstar) - 1L) * n2 # vector indices

  x[idx1] <- 1.0
  x[idx2] <- 1.0

  attr(x, "dim") <- c(n2, nstar)
  x
}

# dup3: using col idx only
# D7 <- dup(7L); x<- apply(D7, 1, function(x) which(x > 0)); matrix(x,7,7)
.dup3 <- function(n = 1L) {
  if ((n < 1L) | (round(n) != n)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }

  if (n > 255L) {
    lav_msg_stop(gettext("n is too large"))
  }

  nstar <- n * (n + 1) / 2
  n2 <- n * n
  # THIS is the real bottleneck: allocating an ocean of zeroes...
  x <- numeric(n2 * nstar)

  tmp <- matrix(0L, n, n)
  tmp[lav_matrix_vech_idx(n)] <- 1:nstar
  tmp[lav_matrix_vechru_idx(n)] <- 1:nstar

  idx <- (1:n2) + (lav_matrix_vec(tmp) - 1L) * n2

  x[idx] <- 1.0

  attr(x, "dim") <- c(n2, nstar)
  x
}


# dup4: using Matrix package, returning a sparse matrix
# .dup4 <- function(n = 1L) {
#    if ((n < 1L) | (round(n) != n)) {
#        stop("n must be a positive integer")
#    }
#
#    if(n > 255L) {
#        stop("n is too large")
#    }
#
#    nstar <- n * (n+1)/2
#    #n2    <- n * n
#
#    tmp <- matrix(0L, n, n)
#    tmp[lav_matrix_vech_idx(n)] <- 1:nstar
#    tmp[lav_matrix_vechru_idx(n)] <- 1:nstar
#
#    x <- Matrix::sparseMatrix(i = 1:(n*n), j = vec(tmp), x = 1.0)
#
#    x
# }

# default dup:
lav_matrix_duplication <- .dup3

# duplication matrix for correlation matrices:
# - it returns a matrix of size p^2 * (p*(p-1))/2
# - the columns corresponding to the diagonal elements have been removed
lav_matrix_duplication_cor <- function(n = 1L) {
  out <- lav_matrix_duplication(n = n)
  diag.idx <- lav_matrix_diagh_idx(n = n)
  out[, -diag.idx, drop = FALSE]
}

# compute t(D) %*% A (without explicitly computing D)
# sqrt(nrow(A)) is an integer
# A is not symmetric, and not even square, only n^2 ROWS
lav_matrix_duplication_pre <- function(A = matrix(0, 0, 0)) {
  # number of rows
  n2 <- NROW(A)

  # square nrow(A) only, n2 = n^2
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)

  OUT <- A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]
  u <- which(idx1 %in% idx2)
  OUT[u, ] <- OUT[u, ] / 2.0

  OUT
}

# dupr_pre is faster...
lav_matrix_duplication_dup_pre2 <- function(A = matrix(0, 0, 0)) {
  # number of rows
  n2 <- NROW(A)

  # square nrow(A) only, n2 = n^2
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)

  OUT <- A[idx1, , drop = FALSE]
  u <- which(!idx1 %in% idx2)
  OUT[u, ] <- OUT[u, ] + A[idx2[u], ]

  OUT
}

# compute t(D) %*% A (without explicitly computing D)
# sqrt(nrow(A)) is an integer
# A is not symmetric, and not even square, only n^2 ROWS
# correlation version: ignoring diagonal elements
lav_matrix_duplication_cor_pre <- function(A = matrix(0, 0, 0)) {
  # number of rows
  n2 <- NROW(A)

  # square nrow(A) only, n2 = n^2
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n, diagonal = FALSE)
  idx2 <- lav_matrix_vechru_idx(n, diagonal = FALSE)

  OUT <- A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]
  u <- which(idx1 %in% idx2)
  OUT[u, ] <- OUT[u, ] / 2.0

  OUT
}

# compute A %*% D (without explicitly computing D)
# sqrt(ncol(A)) must be an integer
# A is not symmetric, and not even square, only n^2 COLUMNS
lav_matrix_duplication_post <- function(A = matrix(0, 0, 0)) {
  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)

  OUT <- A[, idx1, drop = FALSE] + A[, idx2, drop = FALSE]
  u <- which(idx1 %in% idx2)
  OUT[, u] <- OUT[, u] / 2.0

  OUT
}

# compute A %*% D (without explicitly computing D)
# sqrt(ncol(A)) must be an integer
# A is not symmetric, and not even square, only n^2 COLUMNS
# correlation version: ignoring the diagonal elements
lav_matrix_duplication_cor_post <- function(A = matrix(0, 0, 0)) {
  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n, diagonal = FALSE)
  idx2 <- lav_matrix_vechru_idx(n, diagonal = FALSE)

  OUT <- A[, idx1, drop = FALSE] + A[, idx2, drop = FALSE]
  u <- which(idx1 %in% idx2)
  OUT[, u] <- OUT[, u] / 2.0

  OUT
}


# compute t(D) %*% A %*% D (without explicitly computing D)
# A must be a square matrix and sqrt(ncol) an integer
lav_matrix_duplication_pre_post <- function(A = matrix(0, 0, 0)) {
  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(NROW(A) == n2, sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)

  OUT <- A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]
  u <- which(idx1 %in% idx2)
  OUT[u, ] <- OUT[u, ] / 2.0
  OUT <- OUT[, idx1, drop = FALSE] + OUT[, idx2, drop = FALSE]
  OUT[, u] <- OUT[, u] / 2.0

  OUT
}

# compute t(D) %*% A %*% D (without explicitly computing D)
# A must be a square matrix and sqrt(ncol) an integer
# correlation version: ignoring diagonal elements
lav_matrix_duplication_cor_pre_post <- function(A = matrix(0, 0, 0)) {
  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(NROW(A) == n2, sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n, diagonal = FALSE)
  idx2 <- lav_matrix_vechru_idx(n, diagonal = FALSE)

  OUT <- A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]
  u <- which(idx1 %in% idx2)
  OUT[u, ] <- OUT[u, ] / 2.0
  OUT <- OUT[, idx1, drop = FALSE] + OUT[, idx2, drop = FALSE]
  OUT[, u] <- OUT[, u] / 2.0

  OUT
}

# create the generalized inverse of the duplication matrix (D^+_n):
# it removes the duplicated elements in vec(S) to create vech(S)
#
# D^+ %*% vec(S) == vech(S)
#
# M&N book: page 49
#
# D^+ == solve(t(D_n %*% D_n) %*% t(D_n)

# create first t(DUP.ginv)
.dup_ginv1 <- function(n = 1L) {
  if ((n < 1L) | (round(n) != n)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }

  if (n > 255L) {
    lav_msg_stop(gettext("n is too large"))
  }

  nstar <- n * (n + 1) / 2
  n2 <- n * n
  # THIS is the real bottleneck: allocating an ocean of zeroes...
  x <- numeric(nstar * n2)

  tmp <- matrix(1:(n * n), n, n)
  idx1 <- lav_matrix_vech(tmp) + (0:(nstar - 1L)) * n2
  x[idx1] <- 0.5
  idx2 <- lav_matrix_vechru(tmp) + (0:(nstar - 1L)) * n2
  x[idx2] <- 0.5
  idx3 <- lav_matrix_diag_idx(n) + (lav_matrix_diagh_idx(n) - 1L) * n2
  x[idx3] <- 1.0

  attr(x, "dim") <- c(n2, nstar)
  x <- t(x)
  x
}

# create DUP.ginv without transpose
.dup_ginv2 <- function(n = 1L) {
  if ((n < 1L) | (round(n) != n)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }

  if (n > 255L) {
    lav_msg_stop(gettext("n is too large"))
  }

  nstar <- n * (n + 1) / 2
  n2 <- n * n
  # THIS is the real bottleneck: allocating an ocean of zeroes...
  x <- numeric(nstar * n2)

  x[(lav_matrix_vech_idx(n) - 1L) * nstar + 1:nstar] <- 0.5
  x[(lav_matrix_vechru_idx(n) - 1L) * nstar + 1:nstar] <- 0.5
  x[(lav_matrix_diag_idx(n) - 1L) * nstar + lav_matrix_diagh_idx(n)] <- 1.0

  attr(x, "dim") <- c(nstar, n2)
  x
}

lav_matrix_duplication_ginv <- .dup_ginv2

# pre-multiply with D^+
# number of rows in A must be 'square' (n*n)
lav_matrix_duplication_ginv_pre <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix.default(A)

  # number of rows
  n2 <- NROW(A)

  # square nrow(A) only, n2 = n^2
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)
  nstar <- n * (n + 1) / 2

  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)
  OUT <- (A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]) / 2
  OUT
}

# post-multiply with t(D^+)
# number of columns in A must be 'square' (n*n)
lav_matrix_duplication_ginv_post <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix.default(A)

  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)
  OUT <- (A[, idx1, drop = FALSE] + A[, idx2, drop = FALSE]) / 2
  OUT
}

# pre AND post-multiply with D^+: D^+ %*% A %*% t(D^+)
# for square matrices only, with ncol = nrow = n^2
lav_matrix_duplication_ginv_pre_post <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix.default(A)

  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(NROW(A) == n2, sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)
  OUT <- (A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]) / 2
  OUT <- (OUT[, idx1, drop = FALSE] + OUT[, idx2, drop = FALSE]) / 2
  OUT
}

# pre AND post-multiply with D^+: D^+ %*% A %*% t(D^+)
# for square matrices only, with ncol = nrow = n^2
# - ignoring diagonal elements
lav_matrix_duplication_ginv_cor_pre_post <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix.default(A)

  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  stopifnot(NROW(A) == n2, sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  idx1 <- lav_matrix_vech_idx(n, diagonal = FALSE)
  idx2 <- lav_matrix_vechru_idx(n, diagonal = FALSE)
  OUT <- (A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]) / 2
  OUT <- (OUT[, idx1, drop = FALSE] + OUT[, idx2, drop = FALSE]) / 2
  OUT
}


# create the commutation matrix (K_mn)
# the mn x mx commutation matrix is a permutation matrix which
# transforms vec(A) into vec(A')
#
# K_mn %*% vec(A) == vec(A')
#
# (in Henderson & Searle 1979, it is called the vec-permutation matrix)
# M&N book: pages 46-48
#
# note: K_mn is a permutation matrix, so it is orthogonal: t(K_mn) = K_mn^-1
#       K_nm %*% K_mn == I_mn
#
# it is called the 'commutation' matrix because it enables us to interchange
# ('commute') the two matrices of a Kronecker product, eg
#   K_pm (A %x% B) K_nq == (B %x% A)
#
# important property: it allows us to transform a vec of a Kronecker product
# into the Kronecker product of the vecs (if A is m x n and B is p x q):
#   vec(A %x% B) == (I_n %x% K_qm %x% I_p)(vec A %x% vec B)

# first attempt
.com1 <- function(m = 1L, n = 1L) {
  if ((m < 1L) | (round(m) != m)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }

  if ((n < 1L) | (round(n) != n)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }

  p <- m * n
  x <- numeric(p * p)

  pattern <- rep(c(rep((m + 1L) * n, (m - 1L)), n + 1L), n)
  idx <- c(1L, 1L + cumsum(pattern)[-p])

  x[idx] <- 1.0
  attr(x, "dim") <- c(p, p)

  x
}

lav_matrix_commutation <- .com1

# compute K_n %*% A without explicitly computing K
# K_n = K_nn, so sqrt(nrow(A)) must be an integer!
# = permuting the rows of A
lav_matrix_commutation_pre <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix(A)

  # number of rows of A
  n2 <- nrow(A)

  # K_nn only (n2 = m * n)
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # compute row indices
  # row.idx <- as.integer(t(matrix(1:n2, n, n)))
  row.idx <- rep(1:n, each = n) + (0:(n - 1L)) * n

  OUT <- A[row.idx, , drop = FALSE]
  OUT
}

# compute A %*% K_n without explicitly computing K
# K_n = K_nn, so sqrt(ncol(A)) must be an integer!
# = permuting the columns of A
lav_matrix_commutation_post <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix(A)

  # number of columns of A
  n2 <- ncol(A)

  # K_nn only (n2 = m * n)
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # compute col indices
  # row.idx <- as.integer(t(matrix(1:n2, n, n)))
  col.idx <- rep(1:n, each = n) + (0:(n - 1L)) * n

  OUT <- A[, col.idx, drop = FALSE]
  OUT
}

# compute K_n %*% A %*% K_n without explicitly computing K
# K_n = K_nn, so sqrt(ncol(A)) must be an integer!
# = permuting both the rows AND columns of A
lav_matrix_commutation_pre_post <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix(A)

  # number of columns of A
  n2 <- NCOL(A)

  # K_nn only (n2 = m * n)
  stopifnot(sqrt(n2) == round(sqrt(n2)))

  # dimension
  n <- sqrt(n2)

  # compute col indices
  row.idx <- rep(1:n, each = n) + (0:(n - 1L)) * n
  col.idx <- row.idx

  OUT <- A[row.idx, col.idx, drop = FALSE]
  OUT
}



# compute K_mn %*% A without explicitly computing K
# = permuting the rows of A
lav_matrix_commutation_mn_pre <- function(A, m = 1L, n = 1L) {
  # number of rows of A
  mn <- NROW(A)
  stopifnot(mn == m * n)

  # compute row indices
  # row.idx <- as.integer(t(matrix(1:mn, m, n)))
  row.idx <- rep(1:m, each = n) + (0:(n - 1L)) * m

  OUT <- A[row.idx, , drop = FALSE]
  OUT
}

# N_n == 1/2 (I_n^2 + K_nn)
# see MN page 48
#
# N_n == D_n %*% D^+_n
#
lav_matrix_commutation_Nn <- function(n = 1L) {
  lav_msg_stop(gettext("not implemented yet"))
}

# (simplified) kronecker product for square matrices
lav_matrix_kronecker_square <- function(A, check = TRUE) {
  dimA <- dim(A)
  n <- dimA[1L]
  n2 <- n * n
  if (check) {
    stopifnot(dimA[2L] == n)
  }

  # all possible combinations
  out <- tcrossprod(as.vector(A))

  # break up in n*n pieces, and rearrange
  dim(out) <- c(n, n, n, n)
  out <- aperm(out, perm = c(3, 1, 4, 2))

  # reshape again, to form n2 x n2 matrix
  dim(out) <- c(n2, n2)

  out
}

# (simplified) faster kronecker product for symmetric matrices
# note: not faster, but the logic extends to vech versions
lav_matrix_kronecker_symmetric <- function(S, check = TRUE) {
  dimS <- dim(S)
  n <- dimS[1L]
  n2 <- n * n
  if (check) {
    stopifnot(dimS[2L] == n)
  }

  # all possible combinations
  out <- tcrossprod(as.vector(S))

  # break up in n*(n*n) pieces, and rearrange
  dim(out) <- c(n, n * n, n)
  out <- aperm(out, perm = c(3L, 2L, 1L))

  # reshape again, to form n2 x n2 matrix
  dim(out) <- c(n2, n2)

  out
}

# shortcut for the idiom 't(S2) %*% (S %x% S) %*% S2'
# where S is symmetric, and the rows of S2 correspond to
# the elements of S
# eg - S2 = DELTA (the jacobian dS/dtheta)
lav_matrix_tS2_SxS_S2 <- function(S2, S, check = TRUE) {
  # size of S
  n <- NROW(S)

  if (check) {
    stopifnot(NROW(S2) == n * n)
  }

  A <- matrix(S %*% matrix(S2, n, ), n * n, )
  A2 <- A[rep(1:n, each = n) + (0:(n - 1L)) * n, , drop = FALSE]
  crossprod(A, A2)
}

# shortcut for the idiom 't(D) %*% (S %x% S) %*% D'
# where S is symmetric, and D is the duplication matrix
# lav_matrix_tD_SxS_D <- function(S) {

# TODO!!

# }

# square root of a positive definite symmetric matrix
lav_matrix_symmetric_sqrt <- function(S = matrix(0, 0, 0)) {
  n <- NROW(S)

  # eigen decomposition, assume symmetric matrix
  S.eigen <- eigen(S, symmetric = TRUE)
  V <- S.eigen$vectors
  d <- S.eigen$values

  # 'fix' slightly negative tiny numbers
  d[d < 0] <- 0.0

  # sqrt the eigenvalues and reconstruct
  S.sqrt <- V %*% diag(sqrt(d), n, n) %*% t(V)

  S.sqrt
}

# orthogonal complement of a matrix A
# see Satorra (1992). Sociological Methodology, 22, 249-278, footnote 3:
#
# To compute such an orthogonal matrix, consider the p* x p* matrix P = I -
# A(A'A)^-1A', which is idempotent of rank p* - q. Consider the singular value
# decomposition P = HVH', where H is a p* x (p* - q) matrix of full column rank,
# and V is a (p* - q) x (p* - q) diagonal matrix. It is obvious that H'A = 0;
# hence, H is the desired orthogonal complement. This method of constructing an
# orthogonal complement was proposed by Heinz Neudecker (1990, pers. comm.).
#
# update YR 21 okt 2014:
# - note that A %*% solve(t(A) %*% A) %*% t(A) == tcrossprod(qr.Q(qr(A)))
# - if we are using qr, we can as well use qr.Q to get the complement
#
lav_matrix_orthogonal_complement <- function(A = matrix(0, 0, 0)) {
  QR <- qr(A)
  ranK <- QR$rank

  # following Heinz Neudecker:
  # n <- nrow(A)
  # P <- diag(n) - tcrossprod(qr.Q(QR))
  # OUT <- svd(P)$u[, seq_len(n - ranK), drop = FALSE]

  Q <- qr.Q(QR, complete = TRUE)
  # get rid of the first ranK columns
  OUT <- Q[, -seq_len(ranK), drop = FALSE]

  OUT
}

# construct block diagonal matrix from a list of matrices
# ... can contain multiple arguments, which will be coerced to a list
#     or it can be a single list (of matrices)
lav_matrix_bdiag <- function(...) {
  if (nargs() == 0L) {
    return(matrix(0, 0, 0))
  }
  dots <- list(...)

  # create list of matrices
  if (is.list(dots[[1]])) {
    mlist <- dots[[1]]
  } else {
    mlist <- dots
  }
  if (length(mlist) == 1L) {
    return(mlist[[1]])
  }

  # more than 1 matrix
  nmat <- length(mlist)
  nrows <- sapply(mlist, NROW)
  crows <- cumsum(nrows)
  ncols <- sapply(mlist, NCOL)
  ccols <- cumsum(ncols)
  trows <- sum(nrows)
  tcols <- sum(ncols)
  x <- numeric(trows * tcols)

  for (m in seq_len(nmat)) {
    if (m > 1L) {
      rcoffset <- trows * ccols[m - 1] + crows[m - 1]
    } else {
      rcoffset <- 0L
    }
    m.idx <- (rep((0:(ncols[m] - 1L)) * trows, each = nrows[m]) +
      rep(1:nrows[m], ncols[m]) + rcoffset)
    x[m.idx] <- mlist[[m]]
  }

  attr(x, "dim") <- c(trows, tcols)
  x
}

# trace of a single square matrix, or the trace of a product of (compatible)
# matrices resulting in a single square matrix
lav_matrix_trace <- function(..., check = TRUE) {
  if (nargs() == 0L) {
    return(as.numeric(NA))
  }
  dots <- list(...)

  # create list of matrices
  if (is.list(dots[[1]])) {
    mlist <- dots[[1]]
  } else {
    mlist <- dots
  }

  # number of matrices
  nMat <- length(mlist)

  # single matrix
  if (nMat == 1L) {
    S <- mlist[[1]]
    if (check) {
      # check if square
      stopifnot(NROW(S) == NCOL(S))
    }
    out <- sum(S[lav_matrix_diag_idx(n = NROW(S))])
  } else if (nMat == 2L) {
    # dimension check is done by '*'
    out <- sum(mlist[[1]] * t(mlist[[2]]))
  } else if (nMat == 3L) {
    A <- mlist[[1]]
    B <- mlist[[2]]
    C <- mlist[[3]]

    # A, B, C
    # below is the logic; to be coded inline
    # DIAG <- numeric( NROW(A) )
    # for(i in seq_len(NROW(A))) {
    #     DIAG[i] <- sum( rep(A[i,], times = NCOL(B)) *
    #                  as.vector(B) *
    #                  rep(C[,i], each=NROW(B)) )
    # }
    # out <- sum(DIAG)

    # FIXME:

    # dimension check is automatic
    B2 <- B %*% C
    out <- sum(A * t(B2))
  } else {
    # nRows <- sapply(mlist, NROW)
    # nCols <- sapply(mlist, NCOL)

    # check if product is ok
    # stopifnot(all(nCols[seq_len(nMat-1L)] == nRows[2:nMat]))

    # check if product is square
    # stopifnot(nRows[1] == nCols[nMat])

    M1 <- mlist[[1]]
    M2 <- mlist[[2]]
    for (m in 3L:nMat) {
      M2 <- M2 %*% mlist[[m]]
    }
    out <- sum(M1 * t(M2))
  }

  out
}

# crossproduct, but handling NAs pairwise, if needed
# otherwise, just call base::crossprod
lav_matrix_crossprod <- function(A, B) {
  # single argument?
  if (missing(B)) {
    if (!anyNA(A)) {
      return(base::crossprod(A))
    }
    B <- A
    # no missings?
  } else if (!anyNA(A) && !anyNA(B)) {
    return(base::crossprod(A, B))
  }

  # A and B must be matrices
  if (!inherits(A, "matrix")) {
    A <- matrix(A)
  }
  if (!inherits(B, "matrix")) {
    B <- matrix(B)
  }

  out <- apply(B, 2L, function(x) colSums(A * x, na.rm = TRUE))

  # only when A is a vector, and B is a matrix, we get back a vector
  # while the result should be a matrix with 1-row
  if (!is.matrix(out)) {
    out <- t(matrix(out))
  }

  out
}


# reduced row echelon form of A
lav_matrix_rref <- function(A, tol = sqrt(.Machine$double.eps)) {
  # MATLAB documentation says rref uses: tol = (max(size(A))*eps *norm(A,inf)
  if (missing(tol)) {
    A.norm <- max(abs(apply(A, 1, sum)))
    tol <- max(dim(A)) * A.norm * .Machine$double.eps
  }

  # check if A is a matrix
  stopifnot(is.matrix(A))

  # dimensions
  nRow <- NROW(A)
  nCol <- NCOL(A)
  pivot <- integer(0L)

  # catch empty matrix
  if (nRow == 0 && nCol == 0) {
    return(matrix(0, 0, 0))
  }

  rowIndex <- colIndex <- 1
  while (rowIndex <= nRow &&
    colIndex <= nCol) {
    # look for largest (in absolute value) element in this column:
    i.below <- which.max(abs(A[rowIndex:nRow, colIndex]))
    i <- i.below + rowIndex - 1L
    p <- A[i, colIndex]

    # check if column is empty
    if (abs(p) <= tol) {
      A[rowIndex:nRow, colIndex] <- 0L # clean up
      colIndex <- colIndex + 1
    } else {
      # store pivot column
      pivot <- c(pivot, colIndex)

      # do we need to swap column?
      if (rowIndex != i) {
        A[c(rowIndex, i), colIndex:nCol] <-
          A[c(i, rowIndex), colIndex:nCol]
      }

      # scale pivot to be 1.0
      A[rowIndex, colIndex:nCol] <- A[rowIndex, colIndex:nCol] / p

      # create zeroes below and above pivot
      other <- seq_len(nRow)[-rowIndex]
      A[other, colIndex:nCol] <-
        A[other, colIndex:nCol] - tcrossprod(
          A[other, colIndex],
          A[rowIndex, colIndex:nCol]
        )
      # next row/col
      rowIndex <- rowIndex + 1
      colIndex <- colIndex + 1
    }
  }

  # rounding?

  list(R = A, pivot = pivot)
}

# non-orthonoramal (left) null space basis, using rref
lav_matrix_orthogonal_complement2 <- function(
    A,
    tol = sqrt(.Machine$double.eps)) {
  # left
  A <- t(A)

  # compute rref
  out <- lav_matrix_rref(A = A, tol = tol)

  # number of free columns in R (if any)
  nfree <- NCOL(A) - length(out$pivot)

  if (nfree) {
    R <- out$R

    # remove all-zero rows
    zero.idx <- which(apply(R, 1, function(x) {
      all(abs(x) < tol)
    }))
    if (length(zero.idx) > 0) {
      R <- R[-zero.idx, , drop = FALSE]
    }

    FREE <- R[, -out$pivot, drop = FALSE]
    I <- diag(nfree)
    N <- rbind(-FREE, I)
  } else {
    N <- matrix(0, nrow = NCOL(A), ncol = 0L)
  }

  N
}


# inverse of a non-singular (not necessarily positive-definite) symmetric matrix
# FIXME: error handling?
lav_matrix_symmetric_inverse <- function(S, logdet = FALSE,
                                         Sinv.method = "eigen",
                                         zero.warn = FALSE) {
  # catch zero cols/rows
  zero.idx <- which(colSums(S) == 0 & diag(S) == 0 & rowSums(S) == 0)
  S.orig <- S
  if (length(zero.idx) > 0L) {
    if (zero.warn) {
      lav_msg_warn(gettext("matrix to be inverted contains zero cols/rows"))
    }
    S <- S[-zero.idx, -zero.idx, drop = FALSE]
  }

  P <- NCOL(S)

  if (P == 0L) {
    S.inv <- matrix(0, 0, 0)
    if (logdet) {
      attr(S.inv, "logdet") <- 0
    }
    return(S.inv)
  } else if (P == 1L) {
    tmp <- S[1, 1]
    S.inv <- matrix(1 / tmp, 1, 1)
    if (logdet) {
      if (tmp > 0) {
        attr(S.inv, "logdet") <- log(tmp)
      } else {
        attr(S.inv, "logdet") <- -Inf
      }
    }
  } else if (P == 2L) {
    a11 <- S[1, 1]
    a12 <- S[1, 2]
    a21 <- S[2, 1]
    a22 <- S[2, 2]
    tmp <- a11 * a22 - a12 * a21
    if (tmp == 0) {
      S.inv <- matrix(c(Inf, Inf, Inf, Inf), 2, 2)
      if (logdet) {
        attr(S.inv, "logdet") <- -Inf
      }
    } else {
      S.inv <- matrix(c(a22 / tmp, -a21 / tmp, -a12 / tmp, a11 / tmp), 2, 2)
      if (logdet) {
        if (tmp > 0) {
          attr(S.inv, "logdet") <- log(tmp)
        } else {
          attr(S.inv, "logdet") <- -Inf
        }
      }
    }
  } else if (Sinv.method == "eigen") {
    EV <- eigen(S, symmetric = TRUE)
    # V %*% diag(1/d) %*% V^{-1}, where V^{-1} = V^T
    S.inv <-
      tcrossprod(
        EV$vectors / rep(EV$values, each = length(EV$values)),
        EV$vectors
      )

    # 0.5 version
    # S.inv <- tcrossprod(sweep(EV$vectors, 2L,
    #                          STATS = (1/EV$values), FUN="*"), EV$vectors)

    if (logdet) {
      if (all(EV$values >= 0)) {
        attr(S.inv, "logdet") <- sum(log(EV$values))
      } else {
        attr(S.inv, "logdet") <- as.numeric(NA)
      }
    }
  } else if (Sinv.method == "solve") {
    S.inv <- solve.default(S)
    if (logdet) {
      ev <- eigen(S, symmetric = TRUE, only.values = TRUE)
      if (all(ev$values >= 0)) {
        attr(S.inv, "logdet") <- sum(log(ev$values))
      } else {
        attr(S.inv, "logdet") <- as.numeric(NA)
      }
    }
  } else if (Sinv.method == "chol") {
    # this will break if S is not positive definite
    cS <- chol.default(S)
    S.inv <- chol2inv(cS)
    if (logdet) {
      diag.cS <- diag(cS)
      attr(S.inv, "logdet") <- sum(log(diag.cS * diag.cS))
    }
  } else {
    lav_msg_stop(gettext("method must be either `eigen', `solve' or `chol'"))
  }

  if (length(zero.idx) > 0L) {
    logdet <- attr(S.inv, "logdet")
    tmp <- S.orig
    tmp[-zero.idx, -zero.idx] <- S.inv
    S.inv <- tmp
    attr(S.inv, "logdet") <- logdet
    attr(S.inv, "zero.idx") <- zero.idx
  }

  S.inv
}

# update inverse of A, after removing 1 or more rows (and corresponding
# colums) from A
#
# - this is just an application of the inverse of partitioned matrices
# - only removal for now
#
lav_matrix_inverse_update <- function(A.inv, rm.idx = integer(0L)) {
  ndel <- length(rm.idx)

  # rank-1 update
  if (ndel == 1L) {
    a <- A.inv[-rm.idx, rm.idx, drop = FALSE]
    b <- A.inv[rm.idx, -rm.idx, drop = FALSE]
    h <- A.inv[rm.idx, rm.idx]
    out <- A.inv[-rm.idx, -rm.idx, drop = FALSE] - (a %*% b) / h
  }

  # rank-n update
  else if (ndel < NCOL(A.inv)) {
    A <- A.inv[-rm.idx, rm.idx, drop = FALSE]
    B <- A.inv[rm.idx, -rm.idx, drop = FALSE]
    H <- A.inv[rm.idx, rm.idx, drop = FALSE]
    out <- A.inv[-rm.idx, -rm.idx, drop = FALSE] - A %*% solve.default(H, B)

    # erase all col/rows...
  } else if (ndel == NCOL(A.inv)) {
    out <- matrix(0, 0, 0)
  }

  out
}

# update inverse of S, after removing 1 or more rows (and corresponding
# colums) from S, a symmetric matrix
#
# - only removal for now!
#
lav_matrix_symmetric_inverse_update <- function(S.inv, rm.idx = integer(0L),
                                                logdet = FALSE,
                                                S.logdet = NULL) {
  ndel <- length(rm.idx)

  if (ndel == 0L) {
    out <- S.inv
    if (logdet) {
      attr(out, "logdet") <- S.logdet
    }
  }

  # rank-1 update
  else if (ndel == 1L) {
    h <- S.inv[rm.idx, rm.idx]
    a <- S.inv[-rm.idx, rm.idx, drop = FALSE] / sqrt(h)
    out <- S.inv[-rm.idx, -rm.idx, drop = FALSE] - tcrossprod(a)
    if (logdet) {
      attr(out, "logdet") <- S.logdet + log(h)
    }
  }

  # rank-n update
  else if (ndel < NCOL(S.inv)) {
    A <- S.inv[rm.idx, -rm.idx, drop = FALSE]
    H <- S.inv[rm.idx, rm.idx, drop = FALSE]
    out <- (S.inv[-rm.idx, -rm.idx, drop = FALSE] -
      crossprod(A, solve.default(H, A)))
    if (logdet) {
      # cH <- chol.default(Re(H)); diag.cH <- diag(cH)
      # H.logdet <- sum(log(diag.cH * diag.cH))
      H.logdet <- log(det(H))
      attr(out, "logdet") <- S.logdet + H.logdet
    }

    # erase all col/rows...
  } else if (ndel == NCOL(S.inv)) {
    out <- matrix(0, 0, 0)
  } else {
    lav_msg_stop(gettext("column indices exceed number of columns in S.inv"))
  }

  out
}


# update determinant of A, after removing 1 or more rows (and corresponding
# colums) from A
#
lav_matrix_det_update <- function(det.A, A.inv, rm.idx = integer(0L)) {
  ndel <- length(rm.idx)

  # rank-1 update
  if (ndel == 1L) {
    h <- A.inv[rm.idx, rm.idx]
    out <- det.A * h
  }

  # rank-n update
  else if (ndel < NCOL(A.inv)) {
    H <- A.inv[rm.idx, rm.idx, drop = FALSE]
    det.H <- det(H)
    out <- det.A * det.H

    # erase all col/rows...
  } else if (ndel == NCOL(A.inv)) {
    out <- matrix(0, 0, 0)
  }

  out
}

# update determinant of S, after removing 1 or more rows (and corresponding
# colums) from S, a symmetric matrix
#
lav_matrix_symmetric_det_update <- function(det.S, S.inv, rm.idx = integer(0L)) {
  ndel <- length(rm.idx)

  # rank-1 update
  if (ndel == 1L) {
    h <- S.inv[rm.idx, rm.idx]
    out <- det.S * h
  }

  # rank-n update
  else if (ndel < NCOL(S.inv)) {
    H <- S.inv[rm.idx, rm.idx, drop = FALSE]
    cH <- chol.default(H)
    diag.cH <- diag(cH)
    det.H <- prod(diag.cH * diag.cH)
    out <- det.S * det.H

    # erase all col/rows...
  } else if (ndel == NCOL(S.inv)) {
    out <- numeric(0L)
  }

  out
}

# update log determinant of S, after removing 1 or more rows (and corresponding
# colums) from S, a symmetric matrix
#
lav_matrix_symmetric_logdet_update <- function(S.logdet, S.inv,
                                               rm.idx = integer(0L)) {
  ndel <- length(rm.idx)

  # rank-1 update
  if (ndel == 1L) {
    h <- S.inv[rm.idx, rm.idx]
    out <- S.logdet + log(h)
  }

  # rank-n update
  else if (ndel < NCOL(S.inv)) {
    H <- S.inv[rm.idx, rm.idx, drop = FALSE]
    cH <- chol.default(H)
    diag.cH <- diag(cH)
    H.logdet <- sum(log(diag.cH * diag.cH))
    out <- S.logdet + H.logdet

    # erase all col/rows...
  } else if (ndel == NCOL(S.inv)) {
    out <- numeric(0L)
  }

  out
}

# compute `lambda': the smallest root of the determinantal equation
# |M - lambda*P| = 0 (see Fuller 1987, p.125 or p.172
#
# the function allows for zero rows/columns in P, by regressing them out
# this approach was suggested to me by Wayne A. Fuller, personal communication,
# 12 Nov 2020
#
lav_matrix_symmetric_diff_smallest_root <- function(M = NULL, P = NULL,
                                                    warn = FALSE) {
  # check input (we will 'assume' they are square and symmetric)
  stopifnot(is.matrix(M), is.matrix(P))

  # check if P is diagonal or not
  PdiagFlag <- FALSE
  tmp <- P
  diag(tmp) <- 0
  if (all(abs(tmp) < sqrt(.Machine$double.eps))) {
    PdiagFlag <- TRUE
  }

  # diagonal elements of P
  nP <- nrow(P)
  diagP <- P[lav_matrix_diag_idx(nP)]

  # force diagonal elements of P to be nonnegative (warn?)
  neg.idx <- which(diagP < 0)
  if (length(neg.idx) > 0L) {
    if (warn) {
      lav_msg_warn(gettext(
        "some diagonal elements of P are negative (and set to zero)"))
    }
    diag(P)[neg.idx] <- diagP[neg.idx] <- 0
  }

  # check for (near)zero diagonal elements
  zero.idx <- which(abs(diagP) < sqrt(.Machine$double.eps))

  # three cases:
  #    1. all elements are zero (P=0) -> lambda = 0
  #    2. no elements are zero
  #    3. some elements are zero -> regress out

  # 1. all elements are zero
  if (length(zero.idx) == nP) {
    return(0.0)
  }

  # 2. no elements are zero
  else if (length(zero.idx) == 0L) {
    if (PdiagFlag) {
      Ldiag <- 1 / sqrt(diagP)
      LML <- t(Ldiag * M) * Ldiag
    } else {
      L <- solve(lav_matrix_symmetric_sqrt(P))
      LML <- L %*% M %*% t(L)
    }

    # compute lambda
    lambda <- eigen(LML, symmetric = TRUE, only.values = TRUE)$values[nP]

    # 3. some elements are zero
  } else {
    # regress out M-block corresponding to zero diagonal elements in P

    # partition M accordingly: p = positive, n = negative
    M.pp <- M[-zero.idx, -zero.idx, drop = FALSE]
    M.pn <- M[-zero.idx, zero.idx, drop = FALSE]
    M.np <- M[zero.idx, -zero.idx, drop = FALSE]
    M.nn <- M[zero.idx, zero.idx, drop = FALSE]

    # create Mp.n
    Mp.n <- M.pp - M.pn %*% solve(M.nn) %*% M.np

    # extract positive part of P
    P.p <- P[-zero.idx, -zero.idx, drop = FALSE]

    # compute smallest root
    if (PdiagFlag) {
      diagPp <- diag(P.p)
      Ldiag <- 1 / sqrt(diagPp)
      LML <- t(Ldiag * Mp.n) * Ldiag
    } else {
      L <- solve(lav_matrix_symmetric_sqrt(P.p))
      LML <- L %*% Mp.n %*% t(L)
    }
    lambda <- eigen(LML,
      symmetric = TRUE,
      only.values = TRUE
    )$values[nrow(P.p)]
  }

  lambda
}

# force a symmetric matrix to be positive definite
# simple textbook version (see Matrix::nearPD for a more sophisticated version)
#
lav_matrix_symmetric_force_pd <- function(S, tol = 1e-06) {
  if (ncol(S) == 1L) {
    return(matrix(max(S[1, 1], tol), 1L, 1L))
  }

  # eigen decomposition
  S.eigen <- eigen(S, symmetric = TRUE)

  # eigen values
  ev <- S.eigen$values

  # replace small/negative eigen values
  ev[ev / abs(ev[1]) < tol] <- tol * abs(ev[1])

  # reconstruct
  out <- S.eigen$vectors %*% diag(ev) %*% t(S.eigen$vectors)

  out
}

# compute sample covariance matrix, divided by 'N' (not N-1, as in cov)
#
# Mu is not supposed to be ybar, but close
# if provided, we compute S as 1/N*crossprod(Y - Mu) instead of
#                              1/N*crossprod(Y - ybar)
lav_matrix_cov <- function(Y, Mu = NULL) {
  N <- NROW(Y)
  S1 <- stats::cov(Y) # uses a corrected two-pass algorithm
  S <- S1 * (N - 1) / N

  # Mu?
  if (!is.null(Mu)) {
    P <- NCOL(Y)
    ybar <- base::.colMeans(Y, m = N, n = P)
    S <- S + tcrossprod(ybar - Mu)
  }

  S
}

# transform a matrix to match a given target mean/covariance
lav_matrix_transform_mean_cov <- function(Y,
                                          target.mean = numeric(NCOL(Y)),
                                          target.cov = diag(NCOL(Y))) {
  # coerce to matrix
  Y <- as.matrix.default(Y)

  # convert to vector
  target.mean <- as.vector(target.mean)

  S <- lav_matrix_cov(Y)
  S.inv <- solve.default(S)
  S.inv.sqrt <- lav_matrix_symmetric_sqrt(S.inv)
  target.cov.sqrt <- lav_matrix_symmetric_sqrt(target.cov)

  # transform cov
  X <- Y %*% S.inv.sqrt %*% target.cov.sqrt

  # shift mean
  xbar <- colMeans(X)
  X <- t(t(X) - xbar + target.mean)

  X
}

# weighted column means
#
# for each column in Y: mean = sum(wt * Y)/sum(wt)
#
# if we have missing values, we use only the observations and weights
# that are NOT missing
#
lav_matrix_mean_wt <- function(Y, wt = NULL) {
  Y <- unname(as.matrix.default(Y))
  DIM <- dim(Y)

  if (is.null(wt)) {
    return(colMeans(Y, na.rm = TRUE))
  }

  if (anyNA(Y)) {
    WT <- wt * !is.na(Y)
    wN <- .colSums(WT, m = DIM[1], n = DIM[2])
    out <- .colSums(wt * Y, m = DIM[1], n = DIM[2], na.rm = TRUE) / wN
  } else {
    out <- .colSums(wt * Y, m = DIM[1], n = DIM[2]) / sum(wt)
  }

  out
}

# weighted column variances
#
# for each column in Y: var = sum(wt * (Y - w.mean(Y))^2) / N
#
# where N = sum(wt) - 1 (method = "unbiased") assuming wt are frequency weights
#   or  N = sum(wt)     (method = "ML")
#
# Note: another approach (when the weights are 'reliability weights' is to
#       use N = sum(wt) - sum(wt^2)/sum(wt) (not implemented here)
#
# if we have missing values, we use only the observations and weights
# that are NOT missing
#
lav_matrix_var_wt <- function(Y, wt = NULL, method = c("unbiased", "ML")) {
  Y <- unname(as.matrix.default(Y))
  DIM <- dim(Y)

  if (is.null(wt)) {
    wt <- rep(1, nrow(Y))
  }

  if (anyNA(Y)) {
    WT <- wt * !is.na(Y)
    wN <- .colSums(WT, m = DIM[1], n = DIM[2])
    w.mean <- .colSums(wt * Y, m = DIM[1], n = DIM[2], na.rm = TRUE) / wN
    Ytc <- t(t(Y) - w.mean)
    tmp <- .colSums(wt * Ytc * Ytc, m = DIM[1], n = DIM[2], na.rm = TRUE)
    out <- switch(match.arg(method),
      unbiased = tmp / (wN - 1),
      ML = tmp / wN
    )
  } else {
    w.mean <- .colSums(wt * Y, m = DIM[1], n = DIM[2]) / sum(wt)
    Ytc <- t(t(Y) - w.mean)
    tmp <- .colSums(wt * Ytc * Ytc, m = DIM[1], n = DIM[2])
    out <- switch(match.arg(method),
      unbiased = tmp / (sum(wt) - 1),
      ML = tmp / sum(wt)
    )
  }

  out
}

# weighted variance-covariance matrix
#
# always dividing by sum(wt) (for now) (=ML version)
#
# if we have missing values, we use only the observations and weights
# that are NOT missing
#
# same as cov.wt(Y, wt, method = "ML")
#
lav_matrix_cov_wt <- function(Y, wt = NULL) {
  Y <- unname(as.matrix.default(Y))
  DIM <- dim(Y)

  if (is.null(wt)) {
    wt <- rep(1, nrow(Y))
  }

  if (anyNA(Y)) {
    tmp <- na.omit(cbind(Y, wt))
    Y <- tmp[, seq_len(DIM[2]), drop = FALSE]
    wt <- tmp[, DIM[2] + 1L]
    DIM[1] <- nrow(Y)
    w.mean <- .colSums(wt * Y, m = DIM[1], n = DIM[2]) / sum(wt)
    Ytc <- t(t(Y) - w.mean)
    tmp <- crossprod(sqrt(wt) * Ytc)
    out <- tmp / sum(wt)
  } else {
    w.mean <- .colSums(wt * Y, m = DIM[1], n = DIM[2]) / sum(wt)
    Ytc <- t(t(Y) - w.mean)
    tmp <- crossprod(sqrt(wt) * Ytc)
    out <- tmp / sum(wt)
  }

  out
}

# compute (I-A)^{-1} where A is square
# using a (truncated) Neumann series:  (I-A)^{-1} = \sum_k=0^{\infty} A^k
#                                                 = I + A + A^2 + A^3 + ...
#
# note: this only works if the largest eigenvalue for A is < 1; but if A
#       represents regressions, the diagonal will be zero, and all eigenvalues
#       are zero
#
# as A is typically sparse, we can stop if all elements in A^k are zero for,
# say, k<=6
lav_matrix_inverse_iminus <- function(A = NULL) {
  nr <- nrow(A)
  nc <- ncol(A)
  stopifnot(nr == nc)

  # create I + A
  IA <- A
  diag.idx <- lav_matrix_diag_idx(nr)
  IA[diag.idx] <- IA[diag.idx] + 1

  # initial approximation
  IA.inv <- IA

  # first order
  A2 <- A %*% A
  if (all(A2 == 0)) {
    # we are done
    return(IA.inv)
  } else {
    IA.inv <- IA.inv + A2
  }

  # second order
  A3 <- A2 %*% A
  if (all(A3 == 0)) {
    # we are done
    return(IA.inv)
  } else {
    IA.inv <- IA.inv + A3
  }

  # third order
  A4 <- A3 %*% A
  if (all(A4 == 0)) {
    # we are done
    return(IA.inv)
  } else {
    IA.inv <- IA.inv + A4
  }

  # fourth order
  A5 <- A4 %*% A
  if (all(A5 == 0)) {
    # we are done
    return(IA.inv)
  } else {
    IA.inv <- IA.inv + A5
  }

  # fifth order
  A6 <- A5 %*% A
  if (all(A6 == 0)) {
    # we are done
    return(IA.inv)
  } else {
    # naive version (for now)
    tmp <- -A
    tmp[diag.idx] <- tmp[diag.idx] + 1
    IA.inv <- solve(tmp)
    return(IA.inv)
  }
}
