\name{lavTestLRT}
\alias{lavTestLRT}
\alias{lavtestLRT}
\alias{LRT}
\alias{lavLRTTest}
\alias{lavLRT}
\alias{anova}
\title{LRT test}
\description{
LRT test for comparing two (nested) lavaan models.}
\usage{
lavTestLRT(object, ..., SB.classic = TRUE, SB.H0 = FALSE, 
           type = "Chisq", model.names = NULL)
anova(object, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{...}{additional objects of class \code{\linkS4class{lavaan}}.}
\item{SB.classic}{Logical. Only used if the test statistics are scaled (eg. 
  Satorra-Bentler corrected). If \code{TRUE}, the scaled difference test
  statistics are computed by simple (approximate) formulas as described
  in Satorra & Bentler (2001). If \code{FALSE} (the default), 
  the original formulas are used as described in Satorra (2000).}
\item{SB.H0}{Not used yet}
\item{type}{Character. If \code{"Chisq"}, the test statistic for each
  model is the (scaled or unscaled) model fit test statistic. If \code{"Cf"}, 
  the test statistic for each model is computed by the 
  \code{\link{lavTablesFitCf}} function.} 
\item{model.names}{Character vector. If provided, use these model names in
  the first column of the anova table.}
}
\value{
    An object of class anova. When given a single argument, it simply 
    returns the test statistic of this model. When given a sequence of 
    objects, this function tests the models against one another in the order
    specified. 
}
\details{
    The \code{anova} function for lavaan objects simply calls the
    \code{lavTestLRT} function, which has a few additional arguments.

    If \code{type = "Chisq"} and the test statistics are scaled, a
    special scaled difference test is computed as described in 
    Satorra (2000).
}
\references{
Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample
analysis of moment structures. In Heijmans, R.D.H., Pollock, D.S.G. & Satorra,
A. (eds.), Innovations in multivariate statistical analysis. A Festschrift for
Heinz Neudecker (pp.233-247). London: Kluwer Academic Publishers.     

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
statistic for moment structure analysis. Psychometrika, 66(4), 507-514.
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
