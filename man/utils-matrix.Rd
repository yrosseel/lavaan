\name{utils-matrix}
\alias{vech}
\alias{vechru}
\alias{vech.reverse}
\alias{vechru.reverse}
\alias{upper2full}
\alias{vechr}
\alias{vechu}
\alias{vechr.reverse}
\alias{vechu.reverse}
\alias{lower2full}
\alias{duplicationMatrix}
\alias{commutationMatrix}
\alias{sqrtSymmetricMatrix}
\title{Utility Functions For Matrices.}
\description{Utility functions to deal with (mostly symmetric) matrices.}
\usage{
vech(S, diagonal = TRUE)
vechr(S, diagonal = TRUE)
vechu(S, diagonal = TRUE)
vechru(S, diagonal = TRUE)
vech.reverse(x, diagonal = TRUE)
vechru.reverse(x, diagonal = TRUE)
vechr.reverse(x, diagonal = TRUE)
vechu.reverse(x, diagonal = TRUE)
lower2full(x, diagonal = TRUE)
upper2full(x, diagonal = TRUE)
duplicationMatrix(n = 1L)
commutationMatrix(m = 1L, n = 1L)
sqrtSymmetricMatrix(S)
}
\arguments{
\item{S}{A symmetric matrix.}
\item{x}{A numeric vector containing the lower triangular or upper 
triangular elements of a symmetric
matrix, possibly including the diagonal elements.}
\item{diagonal}{Logical. If \code{TRUE}, the diagonal is included. If
\code{FALSE}, the diagonal is not included.}
\item{n}{Integer. Dimension of the symmetric matrix, or column dimension
of a non-square matrix.}
\item{m}{Integer. Row dimension of a matrix.}
}
\details{
The \code{vech} function implements the vech operator 
(for 'half vectorization') and transforms a symmetric matrix 
into a vector by stacking the columns of the matrix one underneath the 
other, but eliminating all supradiagonal elements.
 
The \code{vech.reverse} function does the reverse: given the output of
the \code{vech} function, it reconstructs the symmetric matrix. 

The \code{lower2full} function takes the lower
  
The \code{duplicationMatrix} function creates a duplication matrix D:
it duplicates the elements in vech(S) to create vec(S) 
(where S is symmetric), such that \code{D \%*\% vech(S) == vec(S)}.
 
The \code{commutationMatrix} function creates a commutation matrix (K):
this  mn x mx matrix is a permutation matrix which transforms vec(A) 
into vec(t(A)), such that \code{K \%*\% vec(A) == vec(t(A))}.
 
The \code{sqrtSymmetricMatrix} function computes the square root
of a (positive definite) symmetric matrix.
}
\references{
Magnus, J. R. and H. Neudecker (1999). Matrix Differential Calculus 
with Applications in Statistics and Econometrics, Second Edition, John Wiley.
}
\examples{
# lower.tri elements (including diagonal) of a symmetric matrix
x <- c(4,1,5,2,3,6)

# reconstruct full symmetric matrix (row-wise!)
S <- lower2full(x)

# extract the same lower.tri elements again in the same order
vechr(S)

# without diagonal elements
vechr(S, diagonal=FALSE)

# duplication matrix
nvar <- ncol(S)
vec <- as.vector
Dup <- duplicationMatrix(nvar)
Dup \%*\% vech(S) == vec(S) # should all be true

# commutation matrix
K <- commutationMatrix(nvar, nvar)
K \%*\% vec(S) == vec(t(S)) # should all be true

# take sqrt root of a symmetric matrix
S.sqrt <- sqrtSymmetricMatrix(S)
S.sqrt \%*\% S.sqrt
# should be equal to S again (ignoring some rounding-off errors)
}
