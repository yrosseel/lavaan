\name{lavTestWald}
\alias{lavTestWald}
\alias{lavtestwald}
\alias{wald}
\alias{Wald}
\alias{lavWaldTest}
\title{Wald test}
\description{
Wald test for testing a linear hypothesis about the parameters of 
fitted lavaan object.}
\usage{
lavTestWald(object, constraints = NULL, verbose = FALSE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{constraints}{A character string (typically between single quotes) 
  containing one or more equality constraints. See examples for more details.}
\item{verbose}{Logical. If \code{TRUE}, print out the restriction matrix
  and the estimated restricted values.}
}
\details{
  The constraints are specified using the \code{"=="} operator. Both
  the left-hand side and the right-hand side of the equality can contain
  a linear combination of model parameters, or a constant (like zero).
  The model parameters must be specified by their user-specified labels. 
  Names of defined parameters (using the \code{":="} operator) can be 
  included too.
}
\value{
  A list containing three elements: the Wald test statistic (stat), the
  degrees of freedom (df), and a p-value under the chi-square distribution
  (p.value).
}
\examples{
HS.model <- '
    visual  =~ x1 + b1*x2 + x3
    textual =~ x4 + b2*x5 + x6
    speed   =~ x7 + b3*x8 + x9
'

fit <- cfa(HS.model, data=HolzingerSwineford1939)

# test 1: test about a single parameter
# this is the 'chi-square' version of the 
# z-test from the summary() output
lavTestWald(fit, constraints = "b1 == 0")

# test 2: several constraints
con = '
   2*b1 == b3
   b2 - b3 == 0
'
lavTestWald(fit, constraints = con)
}
