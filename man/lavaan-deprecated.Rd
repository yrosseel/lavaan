\name{lavaan-deprecated}
\alias{lavaan-deprecated}
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
\title{Deprecated Functions in the lavaan package}
\description{These functions have been renamed, or have been removed. They
are still included in this version for convenience, but they may be 
eventually removed.}
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
The \code{vech} function has been renamed \code{lav_matrix_vech}.
 
The \code{vech.reverse} function has been renamed \code{lav_matrix_vech_reverse}.

The \code{lower2full} function has been renamed \code{lav_matrix_lower2full}.
  
The \code{duplicationMatrix} function has been renamed \code{lav_matrix_duplication}.
 
The \code{commutationMatrix} function has been renamed \code{lav_matrix_commutation}.
 
The \code{sqrtSymmetricMatrix} function has been renamed \code{lav_matrix_symmetric_sqrt}.
}
