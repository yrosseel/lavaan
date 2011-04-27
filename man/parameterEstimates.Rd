\name{parameterEstimates}
\alias{parameterEstimates}
\alias{parameterestimates}
\title{Parameter Estimates}
\description{
Parameter estimates of a latent variable model.}
\usage{
parameterEstimates(object)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
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
