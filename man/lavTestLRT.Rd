\name{lavTestLRT}
\alias{lavTestLRT}
\alias{lavtestLRT}
\alias{LRT}
\alias{lavLRTTest}
\alias{lavLRT}
\alias{anova}
\title{LRT test}
\description{
LRT test for comparing (nested) lavaan models.}
\usage{
lavTestLRT(object, ..., method = "default", A.method = "delta",
           scaled.shifted = TRUE,
           H1 = TRUE, type = "Chisq", model.names = NULL)
anova(object, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{...}{additional objects of class \code{\linkS4class{lavaan}}.}
\item{method}{Character string. The possible options are 
\code{"satorra.bentler.2001"}, \code{"satorra.bentler.2010"},
\code{"satorra.2000"}, and \code{"standard"}. See details.}
\item{H1}{Not used yet}
\item{A.method}{Character string. The possible options are \code{"exact"}
and \code{"delta"}. This is only used when method = \code{"satorra.2000"}.
It determines how the Jacobian of the constraint function (the matrix A)
will be computed. Note that if \code{A.method = "exact"}, the models must
    be nested in the parameter sense, while if \code{A.method = "delta"}, they
    only need to be nested in the covariance matrix sense.}
\item{scaled.shifted}{Logical. Only used when method = \code{"satorra.2000"}.
If \code{TRUE}, we use a scaled and shifted test statistic; if \code{FALSE},
we use a mean and variance adjusted (Satterthwaite style) test statistic.}
\item{type}{Character. If \code{"Chisq"}, the test statistic for each
  model is the (scaled or unscaled) model fit test statistic. If \code{"Cf"}, 
  the test statistic for each model is computed by the 
  \code{\link{lavTablesFitCf}} function. If \code{"browne.residual.adf"} (alias
  \code{"browne"}) or \code{"browne.residual.nt"}, the standard chi-squared
  difference is calculated from each model's residual-based statistic.}
\item{model.names}{Character vector. If provided, use these model names in
  the first column of the anova table.}
}
\value{
    An object of class anova. When given a single argument, it simply 
    returns the test statistic of this model. When given a sequence of 
    objects, this function tests the models against one another, after
    reordering the models according to their degrees of freedom.
}
\details{
    The \code{anova} function for lavaan objects simply calls the
    \code{lavTestLRT} function, which has a few additional arguments.

    If \code{type = "Chisq"} and the test statistics are scaled, a
    special scaled difference test statistic is computed. If \code{method} is
    \code{"satorra.bentler.2001"}, a simple approximation is used
    described in Satorra & Bentler (2001). In some settings,
    this can lead to a negative test statistic. To ensure a positive
    test statistic, we can use the method proposed by 
    Satorra & Bentler (2010). Alternatively, when \code{method="satorra.2000"},
    the original formulas of Satorra (2000) are used.  The latter is used for
    model comparison, when \code{...} contain additional (nested) models.
    Even when test statistics are scaled in \code{object} or \code{...},
    users may request the  \code{method="standard"} test statistic,
    without a robust adjustment.
}
\references{
Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample
analysis of moment structures. In Heijmans, R.D.H., Pollock, D.S.G. & Satorra,
A. (eds.), \emph{Innovations in multivariate statistical analysis}:
\emph{A Festschrift for Heinz Neudecker} (pp.233-247). 
London, UK: Kluwer Academic Publishers.     

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
statistic for moment structure analysis. \emph{Psychometrika, 66}(4), 507-514.
\doi{10.1007/BF02296192}

Satorra, A., & Bentler, P. M. (2010). Ensuring postiveness of the scaled
difference chi-square test statistic. \emph{Psychometrika, 75}(2), 243-248.
\doi{10.1007/s11336-009-9135-y}
}
\examples{
HS.model <- '
    visual  =~ x1 + b1*x2 + x3
    textual =~ x4 + b2*x5 + x6
    speed   =~ x7 + b3*x8 + x9
'
fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
fit0 <- cfa(HS.model, data = HolzingerSwineford1939, 
            orthogonal = TRUE)
lavTestLRT(fit1, fit0)
}
