\name{getCov}
\alias{getCov}
\alias{cor2cov}
\alias{char2num}
\title{Utility Functions For Covariance Matrices}
\description{Convenience functions to deal with covariance and correlation
matrices.}
\usage{
getCov(x, lower = TRUE, diagonal = TRUE, sds = NULL,
       names = paste("V", 1:nvar, sep=""))
char2num(s)
cor2cov(R, sds, names = NULL)
}
\arguments{
\item{x}{The elements of the covariance matrix. Either inside a character 
  string or as a numeric vector. In the former case, the function
  \code{char2num} is used to convert the numbers (inside the character
  string) to numeric values.}
\item{lower}{Logical. If \code{TRUE}, the numeric values in \code{x} are
the lower-triangular elements of the (symmetric) covariance matrix only.
If \code{FALSE}, \code{x} contains the upper triangular elements only.
Note we always assumed the elements are provided row-wise!}
\item{diagonal}{Logical. If \code{TRUE}, the numeric values in \code{x}
include the diagonal elements. If \code{FALSE}, a unit diagonal is assumed.}
\item{sds}{A numeric vector containing the standard deviations to be 
used to scale the elements in \code{x} or the correlation matrix \code{R}
into a covariance matrix.}
\item{names}{The variable names of the observed variables.}
\item{s}{Character string containing numeric values; comma's and semi-colons
are ignored.}
\item{R}{A correlation matrix, to be scaled into a covariance matrix.}
}
\details{
The \code{getCov} function is typically used to input the lower 
(or upper) triangular elements of a (symmetric) covariance matrix. In many 
examples found in handbooks, only those elements are shown. However, lavaan 
needs a full matrix to proceed. 

The \code{cor2cov} function is the inverse of the \code{\link[stats]{cov2cor}} 
function, and scales a correlation matrix into a covariance matrix given
the standard deviations of the variables. Optionally, variable names can
be given.
}
\examples{
# The classic Wheaton et. al. (1977) model 
# panel data on he stability of alienation
lower <- '
 11.834,
  6.947,    9.364,
  6.819,    5.091,   12.532,
  4.783,    5.028,    7.495,    9.986,
 -3.839,   -3.889,   -3.841,   -3.625,   9.610,
-21.899,  -18.831,  -21.748,  -18.775,  35.522,  450.288 '

# convert to a full symmetric covariance matrix with names
wheaton.cov <- getCov(lower, names=c("anomia67","powerless67", "anomia71",
                                     "powerless71","education","sei"))

# the model
wheaton.model <- '
  # measurement model
    ses     =~ education + sei
    alien67 =~ anomia67 + powerless67
    alien71 =~ anomia71 + powerless71

  # equations
    alien71 ~ alien67 + ses
    alien67 ~ ses

  # correlated residuals
    anomia67 ~~ anomia71
    powerless67 ~~ powerless71
'

# fitting the model
fit <- sem(wheaton.model, sample.cov=wheaton.cov, sample.nobs=932)

# showing the results
summary(fit, standardized=TRUE)
}
