\name{standardizedSolution}
\alias{standardizedSolution}
\alias{standardizedsolution}
\title{Standardized Solution}
\description{
Standardized solution of a latent variable model.}
\usage{
standardizedSolution(object, type = "std.all", se = TRUE, 
                     zstat = TRUE, pvalue = TRUE, remove.eq = TRUE, 
                     remove.ineq = TRUE, remove.def = FALSE)
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
\item{se}{Logical. If TRUE, standard errors for the standardized parameters
will be computed, together with a z-statistic and a p-value.}
\item{zstat}{Logical. If \code{TRUE}, an extra column is added containing
the so-called z-statistic, which is simply the value of the estimate divided
by its standard error.}
\item{pvalue}{Logical. If \code{TRUE}, an extra column is added containing
the pvalues corresponding to the z-statistic, evaluated under a standard
normal distribution.}
\item{remove.eq}{Logical. If TRUE, filter the output by removing all rows 
containing equality constraints, if any.}
\item{remove.ineq}{Logical. If TRUE, filter the output by removing all
rows containing inequality constraints, if any.}
\item{remove.def}{Logical. If TRUE, filter the ouitput by removing all
rows containing parameter definitions, if any.}
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
