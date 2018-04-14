\name{lav_matrix}
\alias{lav_matrix_vec}
\alias{lav_matrix_vecr}
\alias{lav_matrix_vech}
\alias{lav_matrix_vechr}
\alias{lav_matrix_vechu}
\alias{lav_matrix_vechru}
\alias{lav_matrix_vech_idx}
\alias{lav_matrix_vech_row_idx}
\alias{lav_matrix_vech_col_idx}
\alias{lav_matrix_vechr_idx}
\alias{lav_matrix_vechu_idx}
\alias{lav_matrix_vechru_idx}
\alias{lav_matrix_vech_reverse}
\alias{lav_matrix_vechru_reverse}
\alias{lav_matrix_upper2full}
\alias{lav_matrix_vechr_reverse}
\alias{lav_matrix_vechu_reverse}
\alias{lav_matrix_lower2full}
\alias{lav_matrix_diag_idx}
\alias{lav_matrix_diagh_idx}
\alias{lav_matrix_antidiag_idx}
\alias{lav_matrix_duplication}
\alias{lav_matrix_duplication_pre}
\alias{lav_matrix_duplication_post}
\alias{lav_matrix_duplication_pre_post}
\alias{lav_matrix_duplication_ginv}
\alias{lav_matrix_duplication_ginv_pre}
\alias{lav_matrix_duplication_ginv_post}
\alias{lav_matrix_duplication_ginv_pre_post}
\alias{lav_matrix_commutation}
\alias{lav_matrix_commutation_pre}
\alias{lav_matrix_commutation_mn_pre}
\alias{lav_matrix_symmetric_sqrt}
\alias{lav_matrix_orthogonal_complement}
\alias{lav_matrix_bdiag}
\alias{lav_matrix_trace}
\alias{lav_matrix_cov}

\title{Utility Functions: Matrices and Vectors}
\description{Utility functions for Matrix and Vector operations.}
\usage{
# matrix to vector
lav_matrix_vec(A)
lav_matrix_vecr(A)
lav_matrix_vech(S, diagonal = TRUE)
lav_matrix_vechr(S, diagonal = TRUE)

# matrix/vector indices
lav_matrix_vech_idx(n = 1L, diagonal = TRUE)
lav_matrix_vech_row_idx(n = 1L, diagonal = TRUE)
lav_matrix_vech_col_idx(n = 1L, diagonal = TRUE)
lav_matrix_vechr_idx(n = 1L, diagonal = TRUE)
lav_matrix_vechru_idx(n = 1L, diagonal = TRUE)
lav_matrix_diag_idx(n = 1L)
lav_matrix_diagh_idx(n = 1L)
lav_matrix_antidiag_idx(n = 1L)

# vector to matrix
lav_matrix_vech_reverse(x, diagonal = TRUE)
lav_matrix_vechru_reverse(x, diagonal = TRUE) 
lav_matrix_upper2full(x, diagonal = TRUE)
lav_matrix_vechr_reverse(x, diagonal = TRUE)
lav_matrix_vechu_reverse(x, diagonal = TRUE)
lav_matrix_lower2full(x, diagonal = TRUE)

# the duplication matrix
lav_matrix_duplication(n = 1L)
lav_matrix_duplication_pre(A = matrix(0,0,0))
lav_matrix_duplication_post(A = matrix(0,0,0))
lav_matrix_duplication_pre_post(A = matrix(0,0,0))
lav_matrix_duplication_ginv(n = 1L)
lav_matrix_duplication_ginv_pre(A = matrix(0,0,0))
lav_matrix_duplication_ginv_post(A = matrix(0,0,0))
lav_matrix_duplication_ginv_pre_post(A = matrix(0,0,0))

# the commutation matrix
lav_matrix_commutation(m = 1L, n = 1L)
lav_matrix_commutation_pre(A = matrix(0,0,0))
lav_matrix_commutation_mn_pre(A, m = 1L, n = 1L)

# sample statistics
lav_matrix_cov(Y, ybar = NULL)

# other matrix operations
lav_matrix_symmetric_sqrt(S = matrix(0,0,0))
lav_matrix_orthogonal_complement(A = matrix(0,0,0))
lav_matrix_bdiag(...)
lav_matrix_trace(..., check = TRUE)
}
\arguments{
\item{A}{A general matrix.}
\item{S}{A symmetric matrix.}
\item{Y}{A matrix representing a (numeric) dataset.}
\item{diagonal}{Logical. If TRUE, include the diagonal.}
\item{n}{Integer. When it is the only argument, the dimension of a square 
matrix. If m is also provided, the number of column of the matrix.}
\item{m}{Integer. The number of rows of a matrix.}
\item{x}{Numeric. A vector.}
\item{ybar}{Numeric. The mean vector of the matrix.}
\item{...}{One or more matrices, or a list of matrices.}
\item{check}{Logical. If \code{check = TRUE}, we check if the (final) matrix 
is square.}
}
\details{
These are a collection of lower-level matrix/vector related functions
that are used throughout the lavaan code. They are made public per
request of package developers. Below is a brief description of what
they do:

The \code{lav_matrix_vec} function implements the vec operator (for
'vectorization') and transforms a matrix into a vector by stacking the
columns of the matrix one underneath the other.

The \code{lav_matrix_vecr} function is similar to the \code{lav_matrix_vec}
function but transforms a matrix into a vector by stacking the
rows of the matrix one underneath the other.

The \code{lav_matrix_vech} function implements the vech operator
(for 'half vectorization') and transforms a symmetric matrix
into a vector by stacking the columns of the matrix one underneath the
other, but eliminating all supradiagonal elements. If diagonal = FALSE,
the diagonal elements are also eliminated.

The \code{lav_matrix_vechr} function is similar to the \code{lav_matrix_vech} 
function but transforms a matrix into a vector by stacking the
rows of the matrix one underneath the other, eliminating all 
supradiagonal elements.

The \code{lav_matrix_vech_idx} function returns the vector indices of the lower
triangular elements of a symmetric matrix of size n, column by column.

The \code{lav_matrix_vech_row_idx} function returns the row indices of the
lower triangular elements of a symmetric matrix of size n.

The \code{lav_matrix_vech_col_idx} function returns the column indices of the
lower triangular elements of a symmetric matrix of size n.

The \code{lav_matrix_vechr_idx} function returns the vector indices of the
lower triangular elements of a symmetric matrix of size n, row by row.

The \code{lav_matrix_vechu_idx} function returns the vector indices of the
upper triangular elements of a symmetric matrix of size n, column by column.

The \code{lav_matrix_vechru_idx} function returns the vector indices
of the upper triangular elements of a symmetric matrix of size n, row by row.

The \code{lav_matrix_diag_idx} function returns the vector indices of the
diagonal elements of a symmetric matrix of size n.

The \code{lav_matrix_diagh_idx} function returns the vector indices of
the lower part of a symmetric matrix of size n.

The \code{lav_matrix_antidiag_idx} function returns the vector indices of
the anti diagonal elements a symmetric matrix of size n.

The \code{lav_matrix_vech_reverse} function (alias:
\code{lav_matrix_vechru_reverse} and \code{lav_matrix_upper2full}) creates a
symmetric matrix, given only upper triangular elements, row by row. If
diagonal = FALSE, an diagonal with zero elements is added.

The \code{lav_matrix_vechr_reverse} (alias: \code{lav_matrix_vechu_reverse} and
\code{lav_matrix_lower2full}) creates a symmetric matrix, given only the lower
triangular elements, row by row. If diagonal = FALSE, an diagonal with zero
elements is added.

The \code{lav_matrix_duplication} function generates the duplication matrix
for a symmetric matrix of size n. This matrix duplicates the elements in
vech(S) to create vec(S) (where S is symmetric). This matrix is very
sparse, and should probably never be explicitly created. Use one of 
the functions below.

The \code{lav_matrix_duplication_pre} function computes the product of the
transpose of the duplication matrix and a matrix A. The A matrix should have
n*n rows, where n is an integer. The duplication matrix is not explicitly
created.

The \code{lav_matrix_duplication_post} function computes the product of a
matrix A with the duplication matrix. The A matrix should have n*n columns,
where n is an integer. The duplication matrix is not explicitly created.

The \code{lav_matrix_duplication_pre_post} function first pre-multiplies a
matrix A with the transpose of the duplication matrix, and then post multiplies
the result again with the duplication matrix. A must be square matrix with n*n
rows and columns, where n is an integer. The duplication matrix is not
explicitly created.
multiplies a matrix A with the

The \code{lav_matrix_duplication_ginv} function computes the generalized
inverse of the duplication matrix. The matrix removes the duplicated elements
in vec(S) to create vech(S). This matrix is very sparse, and should probably
never be explicitly created. Use one of the functions below.

The \code{lav_matrix_duplication_ginv_pre} function computes the product of the
generalized inverse of the duplication matrix and a matrix A with n*n rows,
where n is an integer. The generalized inverse of the duplication matrix
is not explicitly created.

The \code{lav_matrix_duplication_ginv_post} function computes the product of a
matrix A (with n*n columns, where n is an integer) and the transpose of the
generalized inverse of the duplication matrix. The generalized inverse of the
duplication matrix is not explicitly created.

The \code{lav_matrix_duplication_ginv_pre_post} function first pre-multiplies
a matrix A with the transpose of the generalized inverse of the duplication
matrix, and the post multiplies the result again with the transpose of the
generalized inverse matrix. The matrix A must be square with n*n rows and
columns, where n is an integer. The generalized inverse of the duplication
matrix is not explicitly created.

The \code{lav_matrix_commutation} function computes the commutation matrix which
is a permutation matrix which transforms vec(A) (with m rows and n columns) 
into vec(t(A)).

The \code{lav_matrix_commutation_pre} function computes the product of the
commutation matrix with a matrix A, without explicitly creating the commutation
matrix. The matrix A must have n*n rows, where n is an integer.

The \code{lav_matrix_commutation_mn_pre} function computes the product of the
commutation matrix with a matrix A, without explicitly creating the commutation
matrix. The matrix A must have m*n rows, where m and n are integers.

The \code{lav_matrix_cov} function computes the sample covariance matrix of
its input matrix, where the elements are divided by N (the number of rows).

The \code{lav_matrix_symmetric_sqrt} function computes the square root of a
positive definite symmetric matrix (using an eigen decomposition).  If some of
the eigenvalues are negative, they are silently fixed to zero.

The \code{lav_matrix_orthogonal_complement} function computes an orthogonal
complement of the matrix A, using a qr decomposition.

The \code{lav_matrix_bdiag} function constructs a block diagonal matrix from
its arguments.

The \code{lav_matrix_trace} function computes the trace (the sum of the
diagonal elements) of a single (square) matrix, or if multiple matrices are
provided (either as a list, or as multiple arguments), we first compute their
product (which must result in a square matrix), and then we compute the trace; 
if \code{check = TRUE}, we check if the (final) matrix is square.
}
\references{
Magnus, J. R. and H. Neudecker (1999). Matrix Differential Calculus
with Applications in Statistics and Econometrics, Second Edition, John Wiley.
}
\examples{
# upper elements of a 3 by 3 symmetric matrix (row by row)
x <- c(30, 16, 5, 10, 3, 1)
# construct full symmetric matrix
S <- lav_matrix_upper2full(x)

# compute the normal theory `Gamma' matrix given a covariance
# matrix (S), using the formula: Gamma = 2 * D^{+} (S %x% S) t(D^{+})
Gamma.NT <- 2 * lav_matrix_duplication_ginv_pre_post(S \%x\% S)
Gamma.NT
}
