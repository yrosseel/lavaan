\name{lavTestWald}
\alias{lavTestWald}
\alias{lavtestwald}
\alias{wald}
\alias{Wald}
\alias{lavWaldTest}
\title{Wald test}
\description{
Wald test for testing a linear hypothesis about the parameters of
a fitted lavaan object.}
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

  The test statistic is based on the variance-covariance matrix of the
  estimated parameters (\code{vcov}) of the fitted object. Therefore, if
  robust standard errors were requested for the fitted object (e.g.,
  \code{se = "robust.sem"} or \code{se = "robust.huber.white"}), the main
  test statistic is already the generalized (\sQuote{robust}) Wald test
  described in Satorra (2000). In that case (or when a scaled test
  statistic was requested, so that the ingredients for a sandwich-type
  covariance matrix are available), additional versions of the Wald
  test statistic are computed as well: the standard (normal-theory)
  statistic (\code{stat.standard}), a mean-scaled statistic
  (\code{stat.scaled}), a mean-and-variance adjusted statistic with
  fractional degrees of freedom (\code{stat.adjusted}), and the
  generalized statistic (\code{stat.robust}); see Satorra (2000).
}
\value{
  A list containing at least three elements: the Wald test statistic
  (\code{stat}), the degrees of freedom (\code{df}), and a p-value under
  the chi-square distribution (\code{p.value}). If robust versions could
  be computed, the list also contains the standard, scaled, adjusted and
  generalized versions of the test statistic, together with their
  p-values (see Details).
}
\references{
Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample
analysis of moment structures. In Heijmans, R.D.H., Pollock, D.S.G. &
Satorra, A. (Eds.), Innovations in multivariate statistical analysis:
A festschrift for Heinz Neudecker (pp. 233-247). London: Kluwer Academic
Publishers.
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
