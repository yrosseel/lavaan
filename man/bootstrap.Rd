\name{bootstrapLavaan}
\alias{bootstrapLavaan}
\alias{bootstrapLRT}
\title{Bootstrapping a Lavaan Model}
\description{Bootstrap functions to get bootstrap standard errors and
bootstrap test statistics}
\usage{
bootstrapLavaan(object, R = 1000, type = "ordinary", verbose = FALSE, 
                coef = TRUE, test = FALSE, FUN = NULL, ...)

bootstrapLRT(h0, h1, R = 1000, type="bollen.stine", verbose=FALSE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{h0}{An object of class \code{\linkS4class{lavaan}}. The restricted
model.}
\item{h1}{An object of class \code{\linkS4class{lavaan}}. The unrestricted
model.}
\item{R}{Integer. The number of bootstrap draws.}
\item{type}{If \code{"ordinary"} or \code{"nonparametric"}, the usual
(naive) bootstrap method is used. If \code{"bollen.stine"}, the
data is first transformed such that the null hypothesis holds exactly
in the resampling space. If \code{"parametric"}, the parametric bootstrap
approach is used. Currently, this is only valid for continuous data
following a multivariate normal distribution.
}
\item{verbose}{If \code{TRUE}, show information for each bootstrap draw.}
\item{coef}{If \code{TRUE}, return the estimated parameters for each
bootstrap draw as a matrix.}
\item{test}{If \code{TRUE}, return the test statistic for each bootstrap
draw as a numeric vector.}
\item{FUN}{Not used for now.}
\item{...}{Not used for now.}
}

\examples{
# fit the Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)

# get the test statistic for the original sample
T.orig <- fitMeasures(fit, "chisq")

# bootstrap to get bootstrap test statistics
# we only generate 10 bootstrap sample in this example; in practice
# you may wish to use a much higher number
T.boot <- bootstrapLavaan(fit, R=10, coef=FALSE, test=TRUE, 
                          type="bollen.stine")

# compute a bootstrap based p-value
pvalue.boot <- length(which(T.boot > T.orig))/length(T.boot)
}
