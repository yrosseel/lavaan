\name{standardizedSolution}
\alias{standardizedSolution}
\alias{standardizedsolution}
\title{Standardized Solution}
\description{
Standardized solution of a latent variable model.}
\usage{
standardizedSolution(object, std.YX = TRUE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{std.YX}{If \code{FALSE}, a single column (est.std) will contain the
standardized estimates based on the variances of the (continuous) latent
variables. If \code{TRUE}, an additional column (est.std.all) will contain 
the standardized estimates based
on both the variances of both (continuous) observed and latent variables.
(Residual) covariances are standardized using (residual) variances.}
}
\value{
  A data.frame containing both unstandardized and standardized model 
  parameters.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
standardizedSolution(fit)
}
