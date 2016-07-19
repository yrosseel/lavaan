\name{bootstrapLavaan}
\alias{bootstrapLavaan}
\alias{bootstrapLRT}
\title{Bootstrapping a Lavaan Model}
\description{Bootstrap the LRT, or any other statistic (or vector
of statistics) you can extract from a fitted lavaan object.}
\usage{
bootstrapLavaan(object, R = 1000L, type = "ordinary", verbose = FALSE, 
                FUN = "coef", warn = -1L, return.boot = FALSE,
                parallel = c("no", "multicore", "snow"),
                ncpus = 1L, cl = NULL, iseed = NULL, h0.rmsea = NULL, ...)

bootstrapLRT(h0 = NULL, h1 = NULL, R = 1000L, type="bollen.stine", 
             verbose = FALSE, return.LRT = FALSE, double.bootstrap = "no", 
             double.bootstrap.R = 500L, double.bootstrap.alpha = 0.05, 
             warn = -1L, parallel = c("no", "multicore", "snow"),
             ncpus = 1L, cl = NULL, iseed = NULL)
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
in the resampling space. If \code{"yuan"}, the data is first transformed 
by combining data and theory (model), such that the resampling space
is closer to the population space.
If \code{"parametric"}, the parametric bootstrap
approach is used; currently, this is only valid for continuous data
following a multivariate normal distribution. See references for more
details.}
\item{FUN}{A function which when applied to the \code{\linkS4class{lavaan}} 
object returns a vector containing the statistic(s) of interest.
The default is \code{FUN="coef"}, returning the estimated values of the
free parameters in the model.}
\item{\dots}{Other named arguments for \code{FUN} which are passed
    unchanged each time it is called.}
\item{verbose}{If \code{TRUE}, show information for each bootstrap draw.}
\item{warn}{Sets the handling of warning messages. See \code{\link{options}}.}
\item{return.boot}{Not used for now.}
\item{return.LRT}{If \code{TRUE}, return the LRT values as an attribute to the pvalue.}
\item{parallel}{The type of parallel operation to be used (if any).  If
missing, the default is \code{"no"}.  }
\item{ncpus}{Integer: number of processes to be used in parallel operation:
typically one would chose this to the number of available CPUs.  }
\item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
\code{parallel = "snow"}.  If not supplied, a cluster on the local machine is
created for the duration of the \code{bootstrapLavaan} or \code{bootstrapLRT}
call.}
\item{iseed}{An integer to set the seed. Or NULL if no reproducible seeds are
  needed. To make this work, make sure the first
  RNGkind() element is \code{"L'Ecuyer-CMRG"}. You can check this by
  typing \code{RNGkind()} in the console. You can set it by typing
  \code{RNGkind("L'Ecuyer-CMRG")}, before the boostrap functions are called.}
\item{h0.rmsea}{Only used if \code{type="yuan"}. Allows one to do the Yuan
bootstrap under the hypothesis that the population RMSEA equals a specified
value.}
\item{double.bootstrap}{If \code{"standard"} the genuine double bootstrap is
used to compute an additional set of plug-in p-values for each boostrap sample.
If \code{"FDB"}, the fast double bootstrap is used to compute   second level
LRT-values for each bootstrap sample. If \code{"no"}, no double bootstrap is
used. The default is set to \code{"FDB"}.}
\item{double.bootstrap.R}{Integer. The number of bootstrap draws to be use for
the double bootstrap.}
\item{double.bootstrap.alpha}{The significance level to compute the adjusted
alpha based on the plugin p-values.}
}
\author{Yves Rosseel and Leonard Vanbrabant. Ed Merkle contributed
Yuan's bootstrap. Improvements to Yuan's bootstrap were contributed
by Hao Wu and Chuchu Cheng.}
\details{The FUN function can return either a scalar or a numeric vector.
This function can be an existing function (for example \code{coef}) or
can be a custom defined function. For example:
\preformatted{
myFUN <- function(x) {
    # require(lavaan)
    modelImpliedCov <- fitted(x)$cov
    vech(modelImpliedCov)
}
}
If \code{parallel="snow"}, it is imperative that the \code{require(lavaan)}
is included in the custom function.}
\references{
Bollen, K. and Stine, R. (1992) Bootstrapping Goodness of Fit Measures in
        Structural Equation Models. Sociological Methods and Research, 21,
        205--229.

Yuan, K.-H., Hayashi, K., & Yanagihara, H. (2007). A class of population
covariance matrices in the bootstrap approach to covariance structure analysis.
Multivariate Behavioral Research, 42, 261--281.
}

\examples{
# fit the Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939, se="none")

# get the test statistic for the original sample
T.orig <- fitMeasures(fit, "chisq")

# bootstrap to get bootstrap test statistics
# we only generate 10 bootstrap sample in this example; in practice
# you may wish to use a much higher number
T.boot <- bootstrapLavaan(fit, R=10, type="bollen.stine",
                          FUN=fitMeasures, fit.measures="chisq")

# compute a bootstrap based p-value
pvalue.boot <- length(which(T.boot > T.orig))/length(T.boot)
}
