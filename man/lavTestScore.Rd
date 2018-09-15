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
             univariate = TRUE, cumulative = FALSE, epc = FALSE,
             verbose = FALSE, warn = TRUE, information = "expected")
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
\item{epc}{Logical. If \code{TRUE}, and we are releasing existing constraints,
compute the expected parameter changes for the existing (free) parameters,
for each released constraint.}
\item{verbose}{Logical. Not used for now.}
\item{warn}{Logical. If \code{TRUE}, print out warnings if they occur.}
\item{information}{\code{character} indicating the type of information matrix
to use (check \code{\link{lavInspect}} for available options). \code{"expected"}
information is the default, which provides better control of Type I errors.}
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
  A list containing at least one \code{data.frame}:
  \itemize{
    \item{\code{$test}: The total score test, with columns for the score
      test statistic (\code{X2}), the degrees of freedom (\code{df}), and
      a \emph{p} value under the \eqn{\chi^2} distribution (\code{p.value}).}
    \item{\code{$uni}: Optional (if \code{univariate=TRUE}).
      Each 1-\emph{df} score test, equivalent to modification indices.}
    \item{\code{$cumulative}: Optional (if \code{cumulative=TRUE}).
      Cumulative score tests.}
    \item{\code{$epc}: Optional (if \code{epc=TRUE}). Parameter estimates,
      expected parameter changes, and expected parameter values if all
      the tested constraints were freed.}
  }
}
\references{
Bentler, P. M., & Chou, C. P. (1993). Some new covariance structure model
improvement statistics. Sage Focus Editions, 154, 235-255.
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

# equivalently, "add" can be a parameter table specifying parameters to free,
# but must include some additional information:
PT.add <- data.frame(lhs = c("visual","textual"),
                     op = c("=~","=~"),
                     rhs = c("x9","x3"),
                     user = 10L, # needed to identify new parameters
                     free = 1, # arbitrary numbers > 0
                     start = 0) # null-hypothesized value
PT.add
lavTestScore(fit, add = PT.add) # same result as above

}
