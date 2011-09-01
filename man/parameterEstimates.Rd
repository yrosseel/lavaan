\name{parameterEstimates}
\alias{parameterEstimates}
\alias{parameterestimates}
\title{Parameter Estimates}
\description{
Parameter estimates of a latent variable model.}
\usage{
parameterEstimates(object, level = 0.95, boot.ci.type = "bca.simple")
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{level}{the confidence level required.}
\item{boot.ci.type}{If bootstrapping was used, the type of interval required.
  The value should be one of \code{"norm"}, \code{"basic"}, \code{"perc"}, 
  or \code{"bca.simple"}. For the first three options, see the help page of
  the \code{\link[boot]{boot.ci}} function in the boot package. The 
  \code{"bca.simple"} option produces intervals using the adjusted bootstrap
  percentile (BCa) method, but with no correction for acceleration (only for 
  bias).}
}
\value{
  A data.frame containing both unstandardized and standardized model 
  parameters, standard errors (for the unstandardized estimates), and
  z-values.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
parameterEstimates(fit)
}
