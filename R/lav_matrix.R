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
lav_matrix_vec <- function(A) {          # nolint
  as.vector(A)
}


# vecr operator
#
# the vecr operator transforms a matrix into
# a vector by stacking the *rows* of the matrix one underneath the other
lav_matrix_vecr <- function(A) {         # nolint
  # faster way??
  # nRow <- NROW(A); nCol <- NCOL(A)
  # idx <- (seq_len(nCol) - 1L) * nRow + rep(seq_len(nRow), each = nCol)
  # if (lav_use_lavaanC() && is.numeric(A)) {
  #   return(lavaanC::m_vecr(A))
  # }
  lav_matrix_vec(t(A))
}


# Delta' A Delta product used by information-matrix calculations.
lav_matrix_delta_A_delta <- function(delta, a1) {
  crossprod(delta, a1) %*% delta
}


# Return matrix vector indices and values from row/column/value triples.
lav_matrix_rowcol_idx <- function(row, col, value, nrow, ncol, symmetric = FALSE) {
  if (length(row) == 0L) {
    return(list(m.idx = integer(0L), x.idx = value[integer(0L)]))
  }

  row <- as.integer(row)
  col <- as.integer(col)
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)

  if (symmetric) {
    upper_idx <- row <= col
    row <- row[upper_idx]
    col <- col[upper_idx]
    value <- value[upper_idx]
  }

  m_idx <- row + (col - 1L) * nrow
  keep_idx <- !duplicated(m_idx, fromLast = TRUE)
  m_idx <- m_idx[keep_idx]
  row <- row[keep_idx]
  col <- col[keep_idx]
  value <- value[keep_idx]

  if (symmetric) {
    offdiag_idx <- row != col
    m_idx <- c(m_idx, col[offdiag_idx] + (row[offdiag_idx] - 1L) * nrow)
    value <- c(value, value[offdiag_idx])
  }

  keep_idx <- value > 0
  m_idx <- m_idx[keep_idx]
  value <- value[keep_idx]

  order_idx <- order(m_idx)
  list(m.idx = m_idx[order_idx], x.idx = value[order_idx])
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
lav_matrix_vech <- function(S, diagonal = TRUE) {     # nolint
  # if (lav_use_lavaanC() && is.numeric(S)) {
  #   return(lavaanC::m_vech(S, diagonal))
  # }
  row_1 <- row(S)
  col_1 <- col(S)
  if (diagonal) S[row_1 >= col_1] else S[row_1 > col_1]
}


# the vechr operator transforms a *symmetric* matrix
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all supradiagonal elements
lav_matrix_vechr <- function(S, diagonal = TRUE) {      # nolint
  # if (lav_use_lavaanC() && is.numeric(S)) {
  #   return(lavaanC::m_vechr(S, diagonal))
  # }
  S[lav_matrix_vechr_idx(n = NCOL(S), diagonal = diagonal)]
}


# the vechu operator transforms a *symmetric* matrix
# into a vector by stacking the *columns* of the matrix one after the
# other, but eliminating all infradiagonal elements
lav_matrix_vechu <- function(S, diagonal = TRUE) {     # nolint
  # if (lav_use_lavaanC() && is.numeric(S)) {
  #   return(lavaanC::m_vechu(S, diagonal))
  # }
  S[lav_matrix_vechu_idx(n = NCOL(S), diagonal = diagonal)]
}

# the vechru operator transforms a *symmetric* matrix
# into a vector by stacking the *rows* of the matrix one after the
# other, but eliminating all infradiagonal elements
#
# same as vech (but using upper-diagonal elements)
lav_matrix_vechru <- function(S, diagonal = TRUE) {     # nolint
  # if (lav_use_lavaanC() && is.numeric(S)) {
  #   return(lavaanC::m_vechru(S, diagonal))
  # }
  S[lav_matrix_vechru_idx(n = NCOL(S), diagonal = diagonal)]
}

# return the *vector* indices of the lower triangular elements of a
# symmetric matrix of size 'n'
lav_matrix_vech_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  # if (lav_use_lavaanC() && n > 1L) {
  #   return(lavaanC::m_vech_idx(n, diagonal))
  # }
  if (n < 1L) {
    return(integer(0L))
  }
  if (diagonal) {
    sequence(n:1L) + rep((seq_len(n) - 1L) * (n + 1L), n:1L)
  } else {
    if (n == 1L) {
      return(integer(0L))
    }
    (sequence((n - 1L):1L) +
      rep((seq_len(n - 1L) - 1L) * (n + 1L) + 1L, (n - 1L):1L))
  }
}

# return the *row* indices of the lower triangular elements of a
# symmetric matrix of size 'n'
lav_matrix_vech_row_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  # if (lav_use_lavaanC() && n > 1L) {
  #   return(lavaanC::m_vech_row_idx(n, diagonal))
  # }
  if (diagonal) {
    sequence(n:1L) + rep(seq_len(n) - 1L, n:1L)
  } else {
    if (n <= 1L) {
      return(integer(0L))
    }
    sequence((n - 1L):1L) + rep(seq_len(n - 1L), (n - 1L):1L)
  }
}

# return the *col* indices of the lower triangular elements of a
# symmetric matrix of size 'n'
lav_matrix_vech_col_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  # if (lav_use_lavaanC() && n > 1L) {
  #   return(lavaanC::m_vech_col_idx(n, diagonal))
  # }
  if (!diagonal) {
    n <- n - 1L
  }
  rep.int(seq_len(n), times = rev(seq_len(n)))
}


# return the *vector* indices of the lower triangular elements of a
# symmetric matrix of size 'n' -- ROW-WISE
lav_matrix_vechr_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  # if (lav_use_lavaanC() && n > 1L) {
  #   return(lavaanC::m_vechr_idx(n, diagonal))
  # }
  if (n < 1L) {
    return(integer(0L))
  }
  if (diagonal) {
    row_idx <- rep(seq_len(n), seq_len(n))
    col_idx <- sequence(seq_len(n))
  } else {
    if (n == 1L) {
      return(integer(0L))
    }
    row_idx <- rep(seq_len(n)[-1L], seq_len(n - 1L))
    col_idx <- sequence(seq_len(n - 1L))
  }
  (col_idx - 1L) * n + row_idx
}

# return the *vector* indices of the upper triangular elements of a
# symmetric matrix of size 'n' -- COLUMN-WISE
lav_matrix_vechu_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  # if (lav_use_lavaanC() && n > 1L) {
  #   return(lavaanC::m_vechu_idx(n, diagonal))
  # }
  if (n < 1L) {
    return(integer(0L))
  }
  if (diagonal) {
    sequence(seq_len(n)) + rep((seq_len(n) - 1L) * n, seq_len(n))
  } else {
    if (n == 1L) {
      return(integer(0L))
    }
    sequence(seq_len(n - 1L)) + rep(seq_len(n - 1L) * n, seq_len(n - 1L))
  }
}

# return the *vector* indices of the upper triangular elements of a
# symmetric matrix of size 'n' -- ROW-WISE
lav_matrix_vechru_idx <- function(n = 1L, diagonal = TRUE) {
  n <- as.integer(n)
  # if (lav_use_lavaanC() && n > 1L) {
  #   return(lavaanC::m_vechru_idx(n, diagonal))
  # }
  if (n < 1L) {
    return(integer(0L))
  }
  if (diagonal) {
    row_idx <- rep(seq_len(n), n:1L)
    col_idx <- sequence(n:1L) + rep(seq_len(n) - 1L, n:1L)
  } else {
    if (n == 1L) {
      return(integer(0L))
    }
    row_idx <- rep(seq_len(n - 1L), (n - 1L):1L)
    col_idx <- sequence((n - 1L):1L) + rep(seq_len(n - 1L), (n - 1L):1L)
  }
  (col_idx - 1L) * n + row_idx
}


# vech.reverse and vechru.reverse (aka `upper2full')
#
# given the output of vech(S) --or vechru(S) which is identical--
# reconstruct S
lav_matrix_vech_reverse <- lav_matrix_vechru_reverse <-
  lav_matrix_upper2full <-
  function(x, diagonal = TRUE) {
    # if (lav_use_lavaanC()) {
    #   return(lavaanC::m_vech_reverse(x, diagonal))
    # }
    # guess dimensions
    if (diagonal) {
      p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    } else {
      p <- (sqrt(1 + 8 * length(x)) + 1) / 2
    }
    stopifnot(p == round(p, 0))

    s <- numeric(p * p)
    s[lav_matrix_vech_idx(p, diagonal = diagonal)] <- x
    s[lav_matrix_vechru_idx(p, diagonal = diagonal)] <- x

    attr(s, "dim") <- c(p, p)
    s
  }


# vechr.reverse vechu.reverse (aka `lower2full')
#
# given the output of vechr(S) --or vechu(S) which is identical--
# reconstruct S
lav_matrix_vechr_reverse <- lav_matrix_vechu_reverse <-
  lav_matrix_lower2full <- function(x, diagonal = TRUE) {
    # if (lav_use_lavaanC()) {
    #   return(lavaanC::m_vechr_reverse(x, diagonal))
    # }
    # guess dimensions
    if (diagonal) {
      p <- (sqrt(1 + 8 * length(x)) - 1) / 2
    } else {
      p <- (sqrt(1 + 8 * length(x)) + 1) / 2
    }
    stopifnot(p == round(p, 0))

    s <- numeric(p * p)
    s[lav_matrix_vechr_idx(p, diagonal = diagonal)] <- x
    s[lav_matrix_vechu_idx(p, diagonal = diagonal)] <- x

    attr(s, "dim") <- c(p, p)
    s
  }


# return the *vector* indices of the diagonal elements of a symmetric
# matrix of size 'n'
lav_matrix_diag_idx <- function(n = 1L) {
  # if(n < 1L) return(integer(0L))
  n <- as.integer(n)
  if (n < 1L) {
    return(integer(0L))
  }
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_diag_idx(n))
  # }
  1L + (seq_len(n) - 1L) * (n + 1L)
}


# return the *vector* indices of the diagonal elements of the LOWER part
# of a symmetric matrix of size 'n'
lav_matrix_diagh_idx <- function(n = 1L) {
  n <- as.integer(n)
  if (n < 1L) {
    return(integer(0L))
  }
  if (n == 1L) {
    return(1L)
  }
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_diagh_idx(n))
  # }
  c(1L, cumsum(n:2L) + 1L)
}


# return the *vector* indices of the ANTI diagonal elements of a symmetric
# matrix of size 'n'
lav_matrix_antidiag_idx <- function(n = 1L) {
  if (n < 1L) {
    return(integer(0L))
  }
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_antidiag_idx(n))
  # }
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
                                      add_idx_at_start = FALSE) {
  if (length(idx) == 0L) {
    return(integer(0L))
  }
  n <- as.integer(n)
  idx <- as.integer(idx)
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_vech_which_idx(n, diagonal,
  #                         idx, type, add.idx.at.start))
  # }
  a <- matrix(FALSE, n, n)
  if (type == "and") {
    a[idx, idx] <- TRUE
  } else if (type == "or") {
    a[idx, ] <- TRUE
    a[, idx] <- TRUE
  }
  pstar_idx <- which(lav_matrix_vech(a, diagonal = diagonal))

  if (add_idx_at_start) {
    pstar_idx <- c(idx, pstar_idx + n)
  }

  pstar_idx
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
  idx <- as.integer(idx)
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_vech_match_idx(n, diagonal, idx))
  # }
  pstar <- n * (n + 1) / 2
  a <- lav_matrix_vech_reverse(seq_len(pstar))
  m_b <- a[idx, idx, drop = FALSE]
  lav_matrix_vech(m_b, diagonal = diagonal)
}

# check if square matrix is diagonal
lav_matrix_is_diagonal <- function(a = NULL, tol = .Machine$double.eps) {
  a <- as.matrix.default(a)
  stopifnot(nrow(a) == ncol(a))
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_is_diagonal(A))
  # }
  diag(a) <- 0
  all(abs(a) <= tol)
}


# create the duplication matrix (D_n): it 'duplicates' the elements
# in vech(S) to create vec(S) (where S is symmetric)
#
# D %*% vech(S) == vec(S)
#
# M&N book: pages 48-50
#

# dup3: using col idx only
# D7 <- dup(7L); x<- apply(D7, 1, function(x) which(x > 0)); matrix(x,7,7)
lav_matrix_duplication <- function(n = 1L) {
  n <- as.integer(n)
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_duplication(n))
  # }
  if ((n < 1L) || (round(n) != n)) {
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
  tmp[lav_matrix_vech_idx(n)] <- seq_len(nstar)
  tmp[lav_matrix_vechru_idx(n)] <- seq_len(nstar)

  idx <- seq_len(n2) + (lav_matrix_vec(tmp) - 1L) * n2

  x[idx] <- 1.0

  attr(x, "dim") <- c(n2, nstar)
  x
}

# compute t(D) %*% A (without explicitly computing D)
# sqrt(nrow(A)) is an integer
# A is not symmetric, and not even square, only n^2 ROWS
lav_matrix_duplication_pre <- function(A = matrix(0, 0, 0)) { # nolint
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_duplication_pre(A))
  # }
  # number of rows
  n2 <- NROW(A)

  # square nrow(A) only, n2 = n^2
  n <- as.integer(round(sqrt(n2)))
  stopifnot(n * n == n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)

  out <- A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]
  u <- lav_matrix_diagh_idx(n)
  out[u, ] <- out[u, ] / 2.0

  out
}

# compute A %*% D (without explicitly computing D)
# sqrt(ncol(A)) must be an integer
# A is not symmetric, and not even square, only n^2 COLUMNS
lav_matrix_duplication_post <- function(A = matrix(0, 0, 0)) { # nolint
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_duplication_post(A))
  # }
  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  n <- as.integer(round(sqrt(n2)))
  stopifnot(n * n == n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)

  out <- A[, idx1, drop = FALSE] + A[, idx2, drop = FALSE]
  u <- lav_matrix_diagh_idx(n)
  out[, u] <- out[, u] / 2.0

  out
}


# compute t(D) %*% A %*% D (without explicitly computing D)
# A must be a square matrix and sqrt(ncol) an integer
lav_matrix_duplication_pre_post <- function(A = matrix(0, 0, 0)) { # nolint
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_duplication_pre_post(A))
  # }
  # number of columns
  n2 <- NCOL(A)

  # square A only, n2 = n^2
  n <- as.integer(round(sqrt(n2)))
  stopifnot(NROW(A) == n2, n * n == n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)

  out <- A[idx1, , drop = FALSE] + A[idx2, , drop = FALSE]
  u <- lav_matrix_diagh_idx(n)
  out[u, ] <- out[u, ] / 2.0
  out <- out[, idx1, drop = FALSE] + out[, idx2, drop = FALSE]
  out[, u] <- out[, u] / 2.0

  out
}

# compute t(D) %*% A %*% D (without explicitly computing D)
# A must be a square matrix and sqrt(ncol) an integer
# correlation version: ignoring diagonal elements
lav_matrix_duplication_cor_pre_post <- function(a = matrix(0, 0, 0)) { # nolint
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_duplication_cor_pre_post(A))
  # }
  # number of columns
  n2 <- NCOL(a)

  # square A only, n2 = n^2
  n <- as.integer(round(sqrt(n2)))
  stopifnot(NROW(a) == n2, n * n == n2)

  # dup idx
  idx1 <- lav_matrix_vech_idx(n, diagonal = FALSE)
  idx2 <- lav_matrix_vechru_idx(n, diagonal = FALSE)

  out <- a[idx1, , drop = FALSE] + a[idx2, , drop = FALSE]
  out <- out[, idx1, drop = FALSE] + out[, idx2, drop = FALSE]

  out
}


# create the generalized inverse of the duplication matrix (D^+_n):
# it removes the duplicated elements in vec(S) to create vech(S)
# if S is symmetric
#
# D^+ %*% vec(S) == vech(S)
#
# M&N book: page 49
#
# D^+ == solve(t(D_n %*% D_n) %*% t(D_n)

# create DUP.ginv without transpose
lav_matrix_duplication_ginv <- function(n = 1L) {
  # if (lav_use_lavaanC()) {
  #   n <- as.integer(n)
  #   return(lavaanC::m_duplication_ginv(n))
  # }
  if ((n < 1L) || (round(n) != n)) {
    lav_msg_stop(gettext("n must be a positive integer"))
  }

  if (n > 255L) {
    lav_msg_stop(gettext("n is too large"))
  }

  nstar <- n * (n + 1) / 2
  n2 <- n * n
  # THIS is the real bottleneck: allocating an ocean of zeroes...
  x <- numeric(nstar * n2)

  x[(lav_matrix_vech_idx(n) - 1L) * nstar + seq_len(nstar)] <- 0.5
  x[(lav_matrix_vechru_idx(n) - 1L) * nstar + seq_len(nstar)] <- 0.5
  x[(lav_matrix_diag_idx(n) - 1L) * nstar + lav_matrix_diagh_idx(n)] <- 1.0

  attr(x, "dim") <- c(nstar, n2)
  x
}

# pre-multiply with D^+
# number of rows in A must be 'square' (n*n)
lav_matrix_duplication_ginv_pre <- function(A = matrix(0, 0, 0)) { # nolint
  m_a <- as.matrix.default(A)

  # number of rows
  n2 <- NROW(m_a)

  # square nrow(A) only, n2 = n^2
  n <- as.integer(round(sqrt(n2)))
  stopifnot(n * n == n2)
  nstar <- n * (n + 1) / 2

  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)
  out <- (m_a[idx1, , drop = FALSE] + m_a[idx2, , drop = FALSE]) / 2
  out
}

# post-multiply with t(D^+)
# number of columns in A must be 'square' (n*n)
lav_matrix_duplication_ginv_post <- function(A = matrix(0, 0, 0)) { #  nolint
  m_a <- as.matrix.default(A)

  # number of columns
  n2 <- NCOL(m_a)

  # square A only, n2 = n^2
  n <- as.integer(round(sqrt(n2)))
  stopifnot(n * n == n2)

  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)
  out <- (m_a[, idx1, drop = FALSE] + m_a[, idx2, drop = FALSE]) / 2
  out
}

# pre AND post-multiply with D^+: D^+ %*% A %*% t(D^+)
# for square matrices only, with ncol = nrow = n^2
lav_matrix_duplication_ginv_pre_post <- function(A = matrix(0, 0, 0)) { # nolint
  m_a <- as.matrix.default(A)
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_duplication_ginv_pre_post(A))
  # }
  # number of columns
  n2 <- NCOL(m_a)

  # square A only, n2 = n^2
  n <- as.integer(round(sqrt(n2)))
  stopifnot(NROW(m_a) == n2, n * n == n2)

  idx1 <- lav_matrix_vech_idx(n)
  idx2 <- lav_matrix_vechru_idx(n)
  out <- (m_a[idx1, , drop = FALSE] + m_a[idx2, , drop = FALSE]) / 2
  out <- (out[, idx1, drop = FALSE] + out[, idx2, drop = FALSE]) / 2
  out
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
lav_matrix_commutation <- function(m = 1L, n = 1L) {
  # if (lav_use_lavaanC()) {
  #   m <- as.integer(m)
  #   n <- as.integer(n)
  #   return(lavaanC::m_commutation(m, n))
  # }

  if ((m < 1L) || (round(m) != m)) {
    lav_msg_stop(gettext("m must be a positive integer"))
  }

  if ((n < 1L) || (round(n) != n)) {
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

# compute K_n %*% A without explicitly computing K
# K_n = K_nn, so sqrt(nrow(A)) must be an integer!
# = permuting the rows of A
lav_matrix_commutation_pre <- function(A = matrix(0, 0, 0)) {
  A <- as.matrix(A)

  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_commutation_pre(A))
  # }

  # number of rows of A
  n2 <- nrow(A)

  # K_nn only (n2 = m * n)
  n <- as.integer(round(sqrt(n2)))
  stopifnot(n * n == n2)

  # compute row indices
  # row.idx <- as.integer(t(matrix(1:n2, n, n)))
  row_idx <- rep(seq_len(n), each = n) + (seq_len(n) - 1L) * n

  out <- A[row_idx, , drop = FALSE]
  out
}

# compute A %*% K_n without explicitly computing K
# K_n = K_nn, so sqrt(ncol(A)) must be an integer!
# = permuting the columns of A
lav_matrix_commutation_post <- function(A = matrix(0, 0, 0)) { # nolint
  m_a <- as.matrix(A)

  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_commutation_post(A))
  # }

  # number of columns of A
  n2 <- ncol(m_a)

  # K_nn only (n2 = m * n)
  n <- as.integer(round(sqrt(n2)))
  stopifnot(n * n == n2)

  # compute col indices
  # row.idx <- as.integer(t(matrix(1:n2, n, n)))
  col_idx <- rep(seq_len(n), each = n) + (seq_len(n) - 1L) * n

  out <- m_a[, col_idx, drop = FALSE]
  out
}

# compute K_n %*% A %*% K_n without explicitly computing K
# K_n = K_nn, so sqrt(ncol(A)) must be an integer!
# = permuting both the rows AND columns of A
lav_matrix_commutation_pre_post <- function(A = matrix(0, 0, 0)) { # nolint
  m_a <- as.matrix(A)
  # if (lav_use_lavaanC()) {
  #   return(lavaanC::m_commutation_pre_post(A))
  # }
  # number of columns of A
  n2 <- NCOL(m_a)

  # K_nn only (n2 = m * n)
  n <- as.integer(round(sqrt(n2)))
  stopifnot(n * n == n2)

  # compute col indices
  row_idx <- rep(seq_len(n), each = n) + (seq_len(n) - 1L) * n
  col_idx <- row_idx

  out <- m_a[row_idx, col_idx, drop = FALSE]
  out
}


# compute K_mn %*% A without explicitly computing K
# = permuting the rows of A
lav_matrix_commutation_mn_pre <- function(A, m = 1L, n = 1L) { # nolint
  # number of rows of A
  mn <- NROW(A)
  stopifnot(mn == m * n)

  # compute row indices
  # row.idx <- as.integer(t(matrix(1:mn, m, n)))
  row_idx <- rep(seq_len(m), each = n) + (seq_len(n) - 1L) * m

  out <- A[row_idx, , drop = FALSE]
  out
}

# shortcut for the idiom 't(D) %*% (S %x% S) %*% D'
# where S is symmetric, and D is the duplication matrix
# lav_matrix_tD_SxS_D <- function(S) {

#TODO: continue coding here

# }

# square root of a positive definite symmetric matrix
lav_matrix_symmetric_sqrt <- function(S = matrix(0, 0, 0)) { # nolint
  n <- NROW(S)

  # eigen decomposition, assume symmetric matrix
  s_eigen <- eigen(S, symmetric = TRUE)
  v <- s_eigen$vectors
  d <- s_eigen$values

  # 'fix' slightly negative tiny numbers
  d[d < 0] <- 0.0

  # sqrt the eigenvalues and reconstruct (scale columns of V,
  #                                 avoid n*n diagonal alloc)
  s_sqrt <- tcrossprod(v * rep(sqrt(d), each = n), v)

  s_sqrt
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
lav_matrix_orthogonal_complement <- function(A = matrix(0, 0, 0)) { # nolint
  qr_1 <- qr(A)
  ran_k <- qr_1$rank

  # following Heinz Neudecker:
  # n <- nrow(A)
  # P <- diag(n) - tcrossprod(qr.Q(QR))
  # OUT <- svd(P)$u[, seq_len(n - ranK), drop = FALSE]

  q_1 <- qr.Q(qr_1, complete = TRUE)
  # get rid of the first ranK columns
  out <- q_1[, -seq_len(ran_k), drop = FALSE]

  out
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
  # empty names will shorten the c/r-names vector
  cnames <- unlist(lapply(mlist, colnames))
  rnames <- unlist(lapply(mlist, rownames))

  x <- numeric(trows * tcols)

  for (m in seq_len(nmat)) {
    if (m > 1L) {
      rcoffset <- trows * ccols[m - 1] + crows[m - 1]
    } else {
      rcoffset <- 0L
    }
    m_idx <- (rep((seq_len(ncols[m]) - 1L) * trows, each = nrows[m]) +
      rep(seq_len(nrows[m]), ncols[m]) + rcoffset)
    x[m_idx] <- mlist[[m]]
  }

  attr(x, "dim") <- c(trows, tcols)
  # col/row names (only if complete)
  if (nrow(x) == length(rnames) && ncol(x) == length(cnames)) {
    attr(x, "dimnames") <- list(rnames, cnames)
  }

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
  n_mat <- length(mlist)

  # single matrix
  if (n_mat == 1L) {
    s <- mlist[[1]]
    if (check) {
      # check if square
      stopifnot(NROW(s) == NCOL(s))
    }
    out <- sum(s[lav_matrix_diag_idx(n = NROW(s))])
  } else if (n_mat == 2L) {
    # dimension check is done by '*'
    out <- sum(mlist[[1]] * t(mlist[[2]]))
  } else if (n_mat == 3L) {
    a <- mlist[[1]]
    m_b <- mlist[[2]]
    m_c <- mlist[[3]]

    # A, B, m_c
    # below is the logic; to be coded inline
    # DIAG <- numeric( NROW(A) )
    # for(i in seq_len(NROW(A))) {
    #     DIAG[i] <- sum( rep(A[i,], times = NCOL(B)) *
    #                  as.vector(B) *
    #                  rep(m_c[,i], each=NROW(B)) )
    # }
    # out <- sum(DIAG)

    # FIXME:

    # dimension check is automatic
    b2 <- m_b %*% m_c
    out <- sum(a * t(b2))
  } else {
    # nRows <- sapply(mlist, NROW)
    # nCols <- sapply(mlist, NCOL)

    # check if product is ok
    # stopifnot(all(nCols[seq_len(nMat-1L)] == nRows[2:nMat]))

    # check if product is square
    # stopifnot(nRows[1] == nCols[nMat])

    m1 <- mlist[[1]]
    m2 <- mlist[[2]]
    for (m in 3L:n_mat) {
      m2 <- m2 %*% mlist[[m]]
    }
    out <- sum(m1 * t(m2))
  }

  out
}

# crossproduct, but handling NAs pairwise, if needed
# otherwise, just call base::crossprod
lav_matrix_crossprod <- function(a, m_b) {
  # single argument?
  if (missing(m_b)) {
    if (!anyNA(a)) {
      return(base::crossprod(a))
    }
    m_b <- a
    # no missings?
  } else if (!anyNA(a) && !anyNA(m_b)) {
    return(base::crossprod(a, m_b))
  }

  # A and B must be matrices
  if (!inherits(a, "matrix")) {
    a <- matrix(a)
  }
  if (!inherits(m_b, "matrix")) {
    m_b <- matrix(m_b)
  }

  out <- apply(m_b, 2L, function(x) colSums(a * x, na.rm = TRUE))

  # only when A is a vector, and B is a matrix, we get back a vector
  # while the result should be a matrix with 1-row
  if (!is.matrix(out)) {
    out <- t(matrix(out))
  }

  out
}


# reduced row echelon form of A
lav_matrix_rref <- function(a, tol = sqrt(.Machine$double.eps)) {
  # MATLAB documentation says rref uses: tol = (max(size(A))*eps *norm(A,inf)
  if (missing(tol)) {
    a_norm <- max(rowSums(abs(a)))
    tol <- max(dim(a)) * a_norm * .Machine$double.eps
  }

  # check if A is a matrix
  stopifnot(is.matrix(a))

  # dimensions
  n_row <- NROW(a)
  n_col <- NCOL(a)
  pivot <- integer(min(n_row, n_col))
  pivot_n <- 0L

  # catch empty matrix
  if (n_row == 0 && n_col == 0) {
    return(matrix(0, 0, 0))
  }

  row_index <- col_index <- 1
  while (row_index <= n_row &&
    col_index <= n_col) {
    # look for largest (in absolute value) element in this column:
    i_below <- which.max(abs(a[row_index:n_row, col_index]))
    i <- i_below + row_index - 1L
    p <- a[i, col_index]

    # check if column is empty
    if (abs(p) <= tol) {
      a[row_index:n_row, col_index] <- 0L # clean up
      col_index <- col_index + 1
    } else {
      # store pivot column
      pivot_n <- pivot_n + 1L
      pivot[pivot_n] <- col_index

      # do we need to swap column?
      if (row_index != i) {
        a[c(row_index, i), col_index:n_col] <-
          a[c(i, row_index), col_index:n_col]
      }

      # scale pivot to be 1.0
      a[row_index, col_index:n_col] <- a[row_index, col_index:n_col] / p

      # create zeroes below and above pivot
      other <- seq_len(n_row)[-row_index]
      a[other, col_index:n_col] <-
        a[other, col_index:n_col] - tcrossprod(
          a[other, col_index],
          a[row_index, col_index:n_col]
        )
      # next row/col
      row_index <- row_index + 1
      col_index <- col_index + 1
    }
  }

  # rounding?

  list(R = a, pivot = pivot[seq_len(pivot_n)])
}


# inverse of a non-singular (not necessarily positive-definite) symmetric matrix
# FIXME: error handling?
lav_matrix_symmetric_inverse <- function(s, logdet = FALSE,
                                         sinv_method = "eigen",
                                         zero_warn = FALSE) {
  # catch zero cols/rows (S is symmetric, so zero row iff zero col)
  zero_idx <- which(rowSums(abs(s)) == 0)
  s_orig <- s
  if (length(zero_idx) > 0L) {
    if (zero_warn) {
      lav_msg_warn(gettext("matrix to be inverted contains zero cols/rows"))
    }
    s <- s[-zero_idx, -zero_idx, drop = FALSE]
  }

  p <- NCOL(s)

  if (p == 0L) {
    s_inv <- matrix(0, 0, 0)
    if (logdet) {
      attr(s_inv, "logdet") <- 0
    }
    return(s_inv)
  } else if (p == 1L) {
    tmp <- s[1, 1]
    s_inv <- matrix(1 / tmp, 1, 1)
    if (logdet) {
      if (tmp > 0) {
        attr(s_inv, "logdet") <- log(tmp)
      } else {
        attr(s_inv, "logdet") <- -Inf
      }
    }
  } else if (p == 2L) {
    a11 <- s[1, 1]
    a12 <- s[1, 2]
    a21 <- s[2, 1]
    a22 <- s[2, 2]
    tmp <- a11 * a22 - a12 * a21
    if (tmp == 0) {
      s_inv <- matrix(c(Inf, Inf, Inf, Inf), 2, 2)
      if (logdet) {
        attr(s_inv, "logdet") <- -Inf
      }
    } else {
      s_inv <- matrix(c(a22 / tmp, -a21 / tmp, -a12 / tmp, a11 / tmp), 2, 2)
      if (logdet) {
        if (tmp > 0) {
          attr(s_inv, "logdet") <- log(tmp)
        } else {
          attr(s_inv, "logdet") <- -Inf
        }
      }
    }
  } else if (sinv_method == "eigen") {
    ev_1 <- eigen(s, symmetric = TRUE)
    # V %*% diag(1/d) %*% V^{-1}, where V^{-1} = V^T
    s_inv <-
      tcrossprod(
        ev_1$vectors / rep(ev_1$values, each = length(ev_1$values)),
        ev_1$vectors
      )

    # 0.5 version
    # S.inv <- tcrossprod(sweep(EV$vectors, 2L,
    #                          STATS = (1/EV$values), FUN="*"), EV$vectors)

    if (logdet) {
      if (all(ev_1$values >= 0)) {
        attr(s_inv, "logdet") <- sum(log(ev_1$values))
      } else {
        attr(s_inv, "logdet") <- as.numeric(NA)
      }
    }
  } else if (sinv_method == "solve") {
    s_inv <- solve.default(s)
    if (logdet) {
      ev <- eigen(s, symmetric = TRUE, only.values = TRUE)
      if (all(ev$values >= 0)) {
        attr(s_inv, "logdet") <- sum(log(ev$values))
      } else {
        attr(s_inv, "logdet") <- as.numeric(NA)
      }
    }
  } else if (sinv_method == "chol") {
    # this will break if S is not positive definite
    c_s <- chol.default(s)
    s_inv <- chol2inv(c_s)
    if (logdet) {
      diag_c_s <- diag(c_s)
      attr(s_inv, "logdet") <- sum(log(diag_c_s * diag_c_s))
    }
  } else {
    lav_msg_stop(gettext("method must be either `eigen', `solve' or `chol'"))
  }

  if (length(zero_idx) > 0L) {
    logdet <- attr(s_inv, "logdet")
    tmp <- s_orig
    tmp[-zero_idx, -zero_idx] <- s_inv
    s_inv <- tmp
    attr(s_inv, "logdet") <- logdet
    attr(s_inv, "zero.idx") <- zero_idx
  }

  s_inv
}

# solve(A) %*% B, where A is symmetric, and B is vector of matrix
lav_matrix_symmetric_solve_spd <- function(a, m_b, tol = 1e-10) {
  chol_a <- try(chol(a), silent = TRUE)
  if (!inherits(chol_a, "try-error")) {
    # Cholesky solve
    backsolve(chol_a, forwardsolve(t(chol_a), m_b))
  } else {
    # Eigen fallback (pseudo-inverse)
    eig <- eigen(a, symmetric = TRUE)
    keep <- eig$values > tol * max(eig$values)
    ainv <- eig$vectors[, keep, drop = FALSE] %*%
      diag(1 / eig$values[keep], length(eig$values[keep])) %*%
      t(eig$vectors[, keep, drop = FALSE])
    ainv %*% m_b
  }
}

# update inverse of S, after removing 1 or more rows (and corresponding
# colums) from S, a symmetric matrix
#
# - only removal for now!
#
lav_matrix_symmetric_inverse_update <- function(s_inv, rm_idx = integer(0L), # nolint
                                                logdet = FALSE,
                                                s_logdet = NULL) {
  ndel <- length(rm_idx)

  if (ndel == 0L) {
    out <- s_inv
    if (logdet) {
      attr(out, "logdet") <- s_logdet
    }

  # rank-1 update
  } else if (ndel == 1L) {
    h <- s_inv[rm_idx, rm_idx]
    if (h > 0) {
      a <- s_inv[-rm_idx, rm_idx, drop = FALSE] / sqrt(h)
      out <- s_inv[-rm_idx, -rm_idx, drop = FALSE] - tcrossprod(a)
      if (logdet) {
        attr(out, "logdet") <- s_logdet + log(h)
      }
    } else {
      out <- s_inv[-rm_idx, -rm_idx, drop = FALSE]
      if (logdet) {
        attr(out, "logdet") <- as.numeric(NA)
      }
    }
  }

  # rank-n update
  else if (ndel < NCOL(s_inv)) {
    a_1 <- s_inv[rm_idx, -rm_idx, drop = FALSE]
    h_1 <- s_inv[rm_idx, rm_idx, drop = FALSE]
    out <- (s_inv[-rm_idx, -rm_idx, drop = FALSE] -
      crossprod(a_1, solve.default(h_1, a_1)))
    if (logdet) {
      c_h <- try(chol.default(h_1), silent = TRUE)
      if (!inherits(c_h, "try-error")) {
        diag_c_h <- diag(c_h)
        h_logdet <- sum(log(diag_c_h * diag_c_h))
      } else {
        ev <- eigen(h_1, symmetric = TRUE, only.values = TRUE)$values
        h_logdet <- sum(log(abs(ev)))
      }
      attr(out, "logdet") <- s_logdet + h_logdet
    }

    # erase all col/rows...
  } else if (ndel == NCOL(s_inv)) {
    out <- matrix(0, 0, 0)
  } else {
    lav_msg_stop(gettext("column indices exceed number of columns in S.inv"))
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
lav_matrix_symmetric_diff_smallest_root <- function(m = NULL, p = NULL) { # nolint
  # check input (we will 'assume' they are square and symmetric)
  stopifnot(is.matrix(m), is.matrix(p))

  # check if P is diagonal or not
  pdiag_flag <- FALSE
  tmp <- p
  diag(tmp) <- 0
  if (all(abs(tmp) < sqrt(.Machine$double.eps))) {
    pdiag_flag <- TRUE
  }

  # diagonal elements of P
  n_p <- nrow(p)
  diag_p <- p[lav_matrix_diag_idx(n_p)]

  # force diagonal elements of P to be nonnegative (warn?)
  neg_idx <- which(diag_p < 0)
  if (length(neg_idx) > 0L) {
    lav_msg_warn(gettext(
      "some diagonal elements of P are negative (and set to zero)"
    ))
    diag(p)[neg_idx] <- diag_p[neg_idx] <- 0
  }

  # check for (near)zero diagonal elements
  zero_idx <- which(abs(diag_p) < sqrt(.Machine$double.eps))

  # three cases:
  #    1. all elements are zero (P=0) -> lambda = 0
  #    2. no elements are zero
  #    3. some elements are zero -> regress out

  # 1. all elements are zero
  if (length(zero_idx) == n_p) {
    return(0.0)
  # 2. no elements are zero
  } else if (length(zero_idx) == 0L) {
    if (pdiag_flag) {
      ldiag <- 1 / sqrt(diag_p)
      lml <- t(ldiag * m) * ldiag
    } else {
      eig_p <- eigen(p, symmetric = TRUE)
      d_p <- pmax(eig_p$values, 0)
      vt_mv <- crossprod(eig_p$vectors, m %*% eig_p$vectors)
      lml <- t(vt_mv / sqrt(d_p)) / sqrt(d_p)
    }

    # compute lambda
    lambda <- eigen(lml, symmetric = TRUE, only.values = TRUE)$values[n_p]

    # 3. some elements are zero
  } else {
    # regress out M-block corresponding to zero diagonal elements in P

    # partition M accordingly: p = positive, n = negative
    m_pp <- m[-zero_idx, -zero_idx, drop = FALSE]
    m_pn <- m[-zero_idx, zero_idx, drop = FALSE]
    m_np <- m[zero_idx, -zero_idx, drop = FALSE]
    m_nn <- m[zero_idx, zero_idx, drop = FALSE]

    # create Mp.n
    mp_n <- m_pp - m_pn %*% lav_matrix_symmetric_solve_spd(m_nn, m_np)

    # extract positive part of P
    p_p <- p[-zero_idx, -zero_idx, drop = FALSE]

    # compute smallest root
    if (pdiag_flag) {
      diag_pp <- diag(p_p)
      ldiag <- 1 / sqrt(diag_pp)
      lml <- t(ldiag * mp_n) * ldiag
    } else {
      eig_pp <- eigen(p_p, symmetric = TRUE)
      d_pp <- pmax(eig_pp$values, 0)
      vt_mv <- crossprod(eig_pp$vectors, mp_n %*% eig_pp$vectors)
      lml <- t(vt_mv / sqrt(d_pp)) / sqrt(d_pp)
    }
    lambda <- eigen(lml,
      symmetric = TRUE,
      only.values = TRUE
    )$values[nrow(p_p)]
  }

  lambda
}

# force a symmetric matrix to be positive definite
# simple textbook version (see Matrix::nearPD for a more sophisticated version)
#
lav_matrix_symmetric_force_pd <- function(s, tol = 1e-06) {
  if (ncol(s) == 1L) {
    return(matrix(max(s[1, 1], tol), 1L, 1L))
  }

  # eigen decomposition
  s_eigen <- eigen(s, symmetric = TRUE)

  # eigen values
  ev <- s_eigen$values

  # replace small/negative eigen values
  ev[ev / abs(ev[1]) < tol] <- tol * abs(ev[1])

  # reconstruct (scale columns of V, avoid n*n diagonal alloc)
  n <- ncol(s_eigen$vectors)
  out <- tcrossprod(s_eigen$vectors * rep(ev, each = n), s_eigen$vectors)

  out
}

# compute sample covariance matrix, divided by 'N' (not N-1, as in cov)
#
# Mu is not supposed to be ybar, but close
# if provided, we compute S as 1/N*crossprod(Y - Mu) instead of
#                              1/N*crossprod(Y - ybar)
lav_matrix_cov <- function(Y, Mu = NULL) {                  # nolint
  n <- NROW(Y)
  s1 <- stats::cov(Y) # uses a corrected two-pass algorithm
  s <- s1 * (n - 1) / n

  # Mu?
  if (!is.null(Mu)) {
    p <- NCOL(Y)
    ybar <- base::.colMeans(Y, m = n, n = p)
    s <- s + tcrossprod(ybar - Mu)
  }

  s
}

# transform a matrix to match a given target mean/covariance
lav_matrix_transform_mean_cov <- function(y,
                                          target_mean = numeric(NCOL(y)),
                                          target_cov = diag(NCOL(y))) {
  # coerce to matrix
  y <- as.matrix.default(y)

  # convert to vector
  target_mean <- as.vector(target_mean)

  s <- lav_matrix_cov(y)
  s_inv <- lav_matrix_symmetric_inverse(s)
  s_inv_sqrt <- lav_matrix_symmetric_sqrt(s_inv)
  target_cov_sqrt <- lav_matrix_symmetric_sqrt(target_cov)

  # transform cov
  x <- y %*% s_inv_sqrt %*% target_cov_sqrt

  # shift mean
  xbar <- colMeans(x)
  x <- t(t(x) - xbar + target_mean)

  x
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
lav_matrix_var_wt <- function(y, wt = NULL, method = c("unbiased", "ML")) {
  y <- unname(as.matrix.default(y))
  dim_1 <- dim(y)

  if (is.null(wt)) {
    wt <- rep(1, nrow(y))
  }

  if (anyNA(y)) {
    wt_1 <- wt * !is.na(y)
    w_n <- .colSums(wt_1, m = dim_1[1], n = dim_1[2])
    w_mean <- .colSums(wt * y, m = dim_1[1], n = dim_1[2], na.rm = TRUE) / w_n
    ytc <- t(t(y) - w_mean)
    tmp <- .colSums(wt * ytc * ytc, m = dim_1[1], n = dim_1[2], na.rm = TRUE)
    out <- switch(match.arg(method),
      unbiased = tmp / (w_n - 1),
      ML = tmp / w_n
    )
  } else {
    w_mean <- .colSums(wt * y, m = dim_1[1], n = dim_1[2]) / sum(wt)
    ytc <- t(t(y) - w_mean)
    tmp <- .colSums(wt * ytc * ytc, m = dim_1[1], n = dim_1[2])
    out <- switch(match.arg(method),
      unbiased = tmp / (sum(wt) - 1),
      ML = tmp / sum(wt)
    )
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
lav_matrix_inverse_iminus <- function(a = NULL) {
  nr <- nrow(a)
  nc <- ncol(a)
  stopifnot(nr == nc)

  # create I + A
  ia <- a
  diag_idx <- lav_matrix_diag_idx(nr)
  ia[diag_idx] <- ia[diag_idx] + 1

  # initial approximation
  ia_inv <- ia

  # first order
  a2 <- a %*% a
  if (all(a2 == 0)) {
    # we are done
    return(ia_inv)
  } else {
    ia_inv <- ia_inv + a2
  }

  # second order
  a3 <- a2 %*% a
  if (all(a3 == 0)) {
    # we are done
    return(ia_inv)
  } else {
    ia_inv <- ia_inv + a3
  }

  # third order
  a4 <- a3 %*% a
  if (all(a4 == 0)) {
    # we are done
    return(ia_inv)
  } else {
    ia_inv <- ia_inv + a4
  }

  # fourth order
  a5 <- a4 %*% a
  if (all(a5 == 0)) {
    # we are done
    return(ia_inv)
  } else {
    ia_inv <- ia_inv + a5
  }

  # fifth order
  a6 <- a5 %*% a
  if (all(a6 == 0)) {
    # we are done
    ia_inv
  } else {
    # naive version (for now)
    tmp <- -a
    tmp[diag_idx] <- tmp[diag_idx] + 1
    ia_inv <- try(solve(tmp), silent = TRUE)
    if (inherits(ia_inv, "try-error")) {
      lav_msg_warn(gettext("(I - A) is singular; using generalized inverse."))
      ia_inv <- MASS::ginv(tmp)
    }
    ia_inv
  }
}

# let Gamma_NT = 2 * D.plus %*% (S %x% S) %*% t(D.plus), where
# D.plus is the elimination matrix: D.plus = ginv(D)
# we wish to compute R = K %*% Gamma_NT %*% t(K)
# but without computing the kronecker product explicitly
# - we allow for meanstructure TRUE/FALSE
# - we allow for fixed.x = TRUE, if x.idx is not empty
# - but NOT for conditional.x = TRUE (for now)
lav_matrix_k_gammant_kt <- function(m_k, s, meanstructure = FALSE,
                                    x_idx = integer(0L)) {
  # dimensions of S
  nvar <- nrow(s)
  pstar <- nvar * (nvar + 1L) / 2L

  # dimension output
  q <- nrow(m_k)

  # sanity checks
  expected_q <- if (meanstructure) nvar + pstar else pstar
  if (ncol(m_k) != expected_q) {
    lav_msg_stop(gettextf(
      "ncol(K) must be %d when meanstructure = %s",
      expected_q, toupper(meanstructure)
    ))
  }
  x_idx <- as.integer(x_idx)
  if (length(x_idx) > 0L && any(x_idx < 1L | x_idx > nvar)) {
    lav_msg_stop(gettext("x.idx contains indices outside 1:nvar."))
  }

  y_idx <- setdiff(seq_len(nvar), x_idx)
  nvar_y <- length(y_idx)
  nvar_x <- length(x_idx)
  if (nvar_y == 0L) {
    lav_msg_stop(gettextf("No random variables remain after removing x.idx."))
  }

  # Cholesky of S (no check: we assume S is pd)
  l_s <- t(chol(s))

  # fixed.x: conditional covariance YbarX and low-rank factor Q for R = Q Q'
  if (nvar_x > 0L) {
    a_yy <- s[y_idx, y_idx, drop = FALSE]
    b_yx <- s[y_idx, x_idx, drop = FALSE]
    c_xx <- s[x_idx, x_idx, drop = FALSE]
    l_x <- t(chol(c_xx))

    bl <- forwardsolve(l_x, t(b_yx))
    ybar_x <- a_yy - t(bl) %*% bl

    f_full <- matrix(0, nvar, nvar_x)
    f_full[y_idx, ] <- b_yx
    f_full[x_idx, ] <- c_xx
    q_1 <- t(forwardsolve(l_x, t(f_full)))
  } else {
    ybar_x <- s
    q_1 <- matrix(0, nvar, 0L)
  }
  l_ybar_x <- t(chol(ybar_x))

  # build W (d^2 x p) with col j = vec(L' N_j L)
  # B2 is q x pstar (rows are vech vectors); L is nvar x d
  build_w <- function(b2, l) {
    d <- ncol(l)
    if (d == 0L) {
      return(matrix(0, 0L, q))
    }

    # expand vech row j of B2 into nvar x nvar_q block matrix V_mat
    v_mat <- matrix(0, nvar, nvar * q)
    for (j in seq_len(q)) {
      m <- lav_matrix_vech_reverse(b2[j, ])
      diag(m) <- diag(m) * 2
      v_mat[, (j - 1L) * nvar + seq_len(nvar)] <- m / 2
    }
    # batch compute t(L) N_j L via two (d x nvar_q) multiplies + aperm trick
    lt_v <- t(l) %*% v_mat
    lt_v3 <- array(lt_v, c(d, nvar, q))
    # aperm transposes each block: since M_i symmetric, t(L' M_j) = M_j L
    lt_vt <- matrix(aperm(lt_v3, c(2L, 1L, 3L)), nvar, d * q)
    matrix(t(l) %*% lt_vt, d^2, q)
  }

  # assemble result
  k2 <- if (meanstructure) {
    m_k[, nvar + seq_len(pstar), drop = FALSE]
  } else {
    m_k
  }
  w_s <- build_w(k2, l_s)
  w_r <- build_w(k2, q_1)
  out <- 2 * (crossprod(w_s) - crossprod(w_r))

  if (meanstructure) {
    k1 <- m_k[, seq_len(nvar), drop = FALSE]
    z <- k1[, y_idx, drop = FALSE] %*% l_ybar_x
    t1 <- tcrossprod(z)
    out <- out + t1
  }

  out
}
