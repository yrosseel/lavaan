\name{lavTestScore}
\alias{lavTestScore}
\alias{lavtestscore}
\alias{score}
\alias{Score}
\alias{lavScoreTest}
\title{Score test}
\description{
Score test (or Lagrange Multiplier test) for releasing one or more
fixed or constrained parameters in model.}
\usage{
lavTestScore(object, add = NULL, release = NULL,
             univariate = TRUE, cumulative = FALSE,
             verbose = FALSE, warn = TRUE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{add}{Either a character string (typically between single quotes) or
  a parameter table containing additional (currently fixed-to-zero) parameters
  for which the score test must be computed.}
\item{release}{Vector of Integers. The indices of the constraints that
should be released. The indices correspond to the order of the equality
constraints as they appear in the parameter table.}
\item{univariate}{Logical. If \code{TRUE}, compute the univariate score
statistics, one for each constraints.}
\item{cumulative}{Logical. If \code{TRUE}, order the univariate score
statistics from large to small, and compute a series of 
multivariate score statistics, each time adding an additional constraint.}
\item{verbose}{Logical. Not used for now.} 
\item{warn}{Logical. If \code{TRUE}, print out warnings if they occur.}
}
\details{
    This function can be used to compute both multivariate and univariate
    score tests. There are two modes: 1) releasing fixed-to-zero parameters
    (using the \code{add} argument), and 2) releasing existing equality 
    constraints (using the \code{release} argument). The two modes can not
    be used simultaneously.

    When adding new parameters, they should not already be part of the model
    (i.e. not listed in the parameter table). If you want to test for 
    a parameter that was explicitly fixed to a constant (say to zero), 
    it is better to label the parameter, and use an explicit equality
    constraint.
}
\value{
  A list containing at least three elements: the Score test statistic 
  (stat), the degrees of freedom (df), and a p-value under the 
  chi-square distribution (p.value). If univariate tests were requested,
  an additional element (TS.univariate) containing a numeric vector of
  univariate score statistics. If cumulative tests were requested,
  an additional element (TS.order) showing the order of the univariate
  test statistics, and an element (TS.cumulative) containing a numeric
  vector of cumulative multivariate score statistics.
}
\examples{
HS.model <- '
    visual  =~ x1 + b1*x2 + x3
    textual =~ x4 + b2*x5 + x6
    speed   =~ x7 + b3*x8 + x9

    b1 == b2
    b2 == b3
'
fit <- cfa(HS.model, data=HolzingerSwineford1939)

# test 1: release both two equality constraints
lavTestScore(fit, cumulative = TRUE)

# test 2: the score test for adding two (currently fixed
# to zero) cross-loadings
newpar = '
    visual =~ x9
    textual =~ x3
'
lavTestScore(fit, add = newpar)
}
