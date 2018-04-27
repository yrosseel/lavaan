\name{standardizedSolution}
\alias{standardizedSolution}
\alias{standardizedsolution}
\title{Standardized Solution}
\description{
Standardized solution of a latent variable model.}
\usage{
standardizedSolution(object, type = "std.all", se = TRUE, zstat = TRUE, 
                     pvalue = TRUE, ci = TRUE, level = 0.95, cov.std = TRUE, 
                     remove.eq = TRUE, remove.ineq = TRUE, remove.def = FALSE, 
                     partable = NULL, GLIST = NULL, est = NULL)
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
\item{ci}{If \code{TRUE}, simple symmetric confidence intervals are added to 
the output}
\item{level}{The confidence level required.}
\item{cov.std}{Logical. If TRUE, the (residual) observed
covariances are scaled by the square root of the `Theta' diagonal elements, and
the (residual) latent covariances are scaled by the square root of the `Psi'
diagonal elements. If FALSE, the (residual) observed covariances
are scaled by the square root of the diagonal elements of the observed
model-implied covariance matrix (Sigma), and the (residual) latent covariances
are scaled by the square root of diagonal elements of the model-implied
covariance matrix of the latent variables.}
\item{remove.eq}{Logical. If TRUE, filter the output by removing all rows 
containing equality constraints, if any.}
\item{remove.ineq}{Logical. If TRUE, filter the output by removing all
rows containing inequality constraints, if any.}
\item{remove.def}{Logical. If TRUE, filter the ouitput by removing all
rows containing parameter definitions, if any.}
\item{GLIST}{List of model matrices. If provided, they will be used
instead of the GLIST inside the object@Model slot.}
\item{est}{Numeric. Parameter values (as in the `est' column of a
parameter table). If provided, they will be used instead of 
the parameters that can be extract from object.}
\item{partable}{A custom \code{list} or \code{data.frame} in which to store
the standardized parameter values. If provided, it will be used instead of 
the parameter table inside the object@ParTable slot.}
}
\note{
The \code{est}, \code{GLIST}, and \code{partable} arguments are not meant for
everyday users, but for authors of external R packages that depend on
\code{lavaan}. Only to be used with great caution.
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
