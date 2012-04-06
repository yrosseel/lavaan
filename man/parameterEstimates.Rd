\name{parameterEstimates}
\alias{parameterEstimates}
\alias{parameterestimates}
\title{Parameter Estimates}
\description{
Parameter estimates of a latent variable model.}
\usage{
parameterEstimates(object, ci = TRUE, level = 0.95, 
                   boot.ci.type = "perc", standardized = FALSE,
                   fmi = "default")
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{ci}{If TRUE, confidence intervals are added to the output}
\item{level}{The confidence level required.}
\item{boot.ci.type}{If bootstrapping was used, the type of interval required.
  The value should be one of \code{"norm"}, \code{"basic"}, \code{"perc"}, 
  or \code{"bca.simple"}. For the first three options, see the help page of
  the \code{\link[boot]{boot.ci}} function in the boot package. The 
  \code{"bca.simple"} option produces intervals using the adjusted bootstrap
  percentile (BCa) method, but with no correction for acceleration (only for 
  bias).}
\item{standardized}{If TRUE, standardized estimates are added to the output}
\item{fmi}{Logical. If \code{TRUE}, an extra column is added containing the
fraction of missing information for each estimated parameter. If
\code{"default"}, the value is set to \code{TRUE} only if
\code{estimator="ML"}, \code{missing="(fi)ml"}, and \code{se="standard"}.}
}
\value{
  A data.frame containing the estimated parameters,
  parameters, standard errors, z-values, and (by default) 
  the lower and upper values of the confidence intervals. If requested,
  extra columns are added with standardized versions of the parameter
  estimates.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
parameterEstimates(fit)
}
