\name{estfun}
\alias{estfun.lavaan}
\alias{lavScores}
\title{Extract Empirical Estimating Functions}
\description{
A function for extracting the empirical estimating functions of a fitted 
lavaan model. This is the derivative of the objective function 
with respect to the parameter vector, evaluated at the observed (case-wise)
data. In other words, this function returns the case-wise scores, evaluated
at the fitted model parameters.
}
\usage{
estfun.lavaan(object, scaling = FALSE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{scaling}{If \code{TRUE}, the scores are scaled to reflect the specific
objective function used by lavaan. If \code{FALSE} (the default), the 
objective function is the loglikelihood function assuming multivariate 
normality.}
}
\value{A n x k matrix corresponding to n observations and k parameters.}
\author{Ed Merkle}
