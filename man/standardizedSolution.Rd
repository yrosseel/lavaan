\name{standardizedSolution}
\alias{standardizedSolution}
\alias{standardizedsolution}
\title{Standardized Solution}
\description{
Standardized solution of a latent variable model.}
\usage{
standardizedSolution(object, type = "std.all")
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{If \code{"std.lv"}, the standardized estimates are
on the variances of the (continuous) latent variables only.
If \code{"std.all"}, the standardized estimates are based
on both the variances of both (continuous) observed and latent variables.
If \code{"std.nox"}, the standardized estimates are based
on both the variances of both (continuous) observed and latent variables, 
but not the variances of exogenous covariates.}
}
\value{
  A data.frame containing standardized model parameters.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
standardizedSolution(fit)
}
