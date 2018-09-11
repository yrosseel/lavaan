\name{lavResiduals}
\alias{lavResiduals}
\alias{lavResidual}
\title{Residuals}
\description{
\sQuote{lavResiduals} provides model residuals and standardized residuals from
a fitted lavaan object, as well as various summaries of these residuals.

The \sQuote{residuals()} (and \sQuote{resid()}) methods are just shortcuts to
this function with a limited set of arguments.
}
\usage{
lavResiduals(object, type = "cor.bentler", se = FALSE, zstat = TRUE, 
    summary = TRUE, h1.acov = "unstructured", add.type = TRUE,
    add.labels = TRUE, add.class = TRUE, drop.list.single.group = TRUE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{Character.
If \code{type = "raw"}, this function returns the raw (= unscaled)
difference between the observed and the expected (model-implied) summary
statistics, as well as the standardized version of these residualds.
If \code{type = "cor"}, or \code{type = "cor.bollen"}, the observed and model
implied covariance matrices are first transformed to a correlation matrix
(using \code{cov2cor()}), before the residuals are computed.
If \code{type = "cor.bentler"}, both the observed and model implied covariance
matrices are rescaled by dividing the elements by the square roots of the
corresponding variances of the observed covariance matrix.}
\item{se}{Logical. If \code{TRUE}, show the estimated standard errors
for the residuals.}
\item{zstat}{Logical. If \code{TRUE}, show the standardized residuals, which
are the raw residuals divided by the corresponding (estimated) standard
errors.}
\item{summary}{Logical. If \code{TRUE}, show various summaries of the
(possibly scaled) residuals. When \code{type = "raw"}, we compute the
RMR. When \code{type = "cor.bentler"},
we compute the SRMR. When \code{type = "cor.bollen"}, we compute the CRMR.
An unbiased version of these summaries is also computed, as well as a 
standard error, a z-statistic and a p-value for the test of exact fit
based on these summaries.}
\item{h1.acov}{Character. If \code{"unstructured"}, the observed summary
statistics are used as consistent estimates of the corresponding (unrestricted)
population statistics. If \code{"structured"}, the model-implied summary
statistics are used as consistent estimates of the corresponding (unrestricted)
population statistics. This affects the way the asymptotic variance matrix
of the summary statistics is computed.}
\item{add.type}{Logical. If \code{TRUE}, show the type of residuals
in the output.}
\item{add.labels}{If \code{TRUE}, variable names are added to the vectors
and/or matrices.}
\item{add.class}{If \code{TRUE}, vectors are given the \sQuote{lavaan.vector} 
class; matrices are given the \sQuote{lavaan.matrix} class, and symmetric 
matrices are given the \sQuote{lavaan.matrix.symmetric} class. 
This only affects the way they are printed on the screen.}
\item{drop.list.single.group}{If \code{FALSE}, the results are returned as
a list, where each element corresponds to a group (even if there is only
a single group). If \code{TRUE}, the list will be unlisted if there is
only a single group.}
}
\value{
If \code{drop.list.single.group = TRUE}, a list of (residualized) summary
statistics, including type, standardized residuals, and summaries. If
\code{drop.list.single.group = FALSE}, the list of summary statistics is nested
within a list for each group.
}
\references{
Bentler, P.M. and Dijkstra, T. (1985). Efficient estimation via linearization
in structural models. In Krishnaiah, P.R. (Ed.), Multivariate analysis - VI,
(pp. 9-42). New York, NY: Elsevier.

Ogasawara, H. (2001). Standard errors of fit indices using residuals in
structural equation modeling. Psychometrika, 66(3), 421-436.

Maydeu-Olivares, A. (2017). Assessing the size of model misfit in structural
equation models. Psychometrika, 82(3), 533-558.

Standardized Residuals in Mplus. Document retrieved from URL
http://www.statmodel.com/download/StandardizedResiduals.pdf
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
lavResiduals(fit) 
}
