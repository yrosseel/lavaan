\name{parameterEstimates}
\alias{parameterEstimates}
\alias{parameterestimates}
\title{Parameter Estimates}
\description{
Parameter estimates of a latent variable model.}
\usage{
parameterEstimates(object, se = TRUE, zstat = TRUE, pvalue = TRUE,
                   ci = TRUE, level = 0.95, boot.ci.type = "perc",
                   standardized = FALSE, cov.std = TRUE, fmi = FALSE,
                   remove.system.eq = TRUE, remove.eq = TRUE,
                   remove.ineq = TRUE, remove.def = FALSE,
                   rsquare = FALSE, add.attributes = FALSE, header = TRUE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{se}{Logical. If \code{TRUE}, include column containing the standard
errors. If \code{FALSE}, this implies \code{zstat} and \code{pvalue} and
\code{ci} are also \code{FALSE}.}
\item{zstat}{Logical. If \code{TRUE}, an extra column is added containing
the so-called z-statistic, which is simply the value of the estimate divided
by its standard error.}
\item{pvalue}{Logical. If \code{TRUE}, an extra column is added containing
the pvalues corresponding to the z-statistic, evaluated under a standard
normal distribution.}
\item{ci}{If \code{TRUE}, confidence intervals are added to the output}
\item{level}{The confidence level required.}
\item{boot.ci.type}{If bootstrapping was used, the type of interval required.
  The value should be one of \code{"norm"}, \code{"basic"}, \code{"perc"},
  or \code{"bca.simple"}. For the first three options, see the help page of
  the \code{boot.ci} function in the boot package. The
  \code{"bca.simple"} option produces intervals using the adjusted bootstrap
  percentile (BCa) method, but with no correction for acceleration (only for
  bias).}
\item{standardized}{Logical. If \code{TRUE}, standardized estimates are
added to the output. Note that \emph{SE}s and tests are still based on
unstandardized estimates. Use \code{\link{standardizedSolution}} to obtain
\emph{SE}s and test statistics for standardized estimates.}
\item{cov.std}{Logical. If TRUE, the (residual) observed
covariances are scaled by the square root of the `Theta' diagonal elements, and
the (residual) latent covariances are scaled by the square root of the `Psi'
diagonal elements. If FALSE, the (residual) observed covariances
are scaled by the square root of the diagonal elements of the observed
model-implied covariance matrix (Sigma), and the (residual) latent covariances
are scaled by the square root of diagonal elements of the model-implied
covariance matrix of the latent variables.}
\item{fmi}{Logical. If \code{TRUE}, an extra column is added containing the
fraction of missing information for each estimated parameter. Only
available if
\code{estimator="ML"}, \code{missing="(fi)ml"}, and \code{se="standard"}.
See references for more information.}
\item{remove.eq}{Logical. If \code{TRUE}, filter the output by removing all
rows containing user-specified equality constraints, if any.}
\item{remove.system.eq}{Logical. If \code{TRUE}, filter the output by
removing all rows containing system-generated equality constraints, if any.}
\item{remove.ineq}{Logical. If \code{TRUE}, filter the output by removing all
rows containing inequality constraints, if any.}
\item{remove.def}{Logical. If \code{TRUE}, filter the ouitput by removing all
rows containing parameter definitions, if any.}
\item{rsquare}{Logical. If \code{TRUE}, add additional rows containing
the rsquare values (in the \code{est} column) of all endogenous variables
in the model. Both the \code{lhs} and \code{rhs} column contain the
name of the endogenous variable, while the code{op} column contains \code{r2},
to indicate that the values in the \code{est} column are rsquare values.}
\item{add.attributes}{Logical. If \code{TRUE}, add a class attribute
(class \code{lavaan.parameterEstimates}) and other attributes to be used by
the print function for this class (\code{print.lavaan.parameterEstimates}).
This is used by the \code{summary()} function, to prettify the output.}
\item{header}{Logical. Only used if \code{add.attributes = TRUE}. If
\code{TRUE}, print a header at the top of the parameter list. This header
contains information about the information matrix, if saturated (h1) model
is structured or unstructured, and which type of standard errors are shown
in the output.}
}
\value{
  A data.frame containing the estimated parameters,
  parameters, standard errors, and (by default) z-values , p-values, and
  the lower and upper values of the confidence intervals. If requested,
  extra columns are added with standardized versions of the parameter
  estimates.
}
\references{
Savalei, V. & Rhemtulla, M. (2012). On obtaining estimates of the fraction of
missing information from FIML. Structural Equation Modeling: A
Multidisciplinary Journal, 19(3), 477-494.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
parameterEstimates(fit)
}
