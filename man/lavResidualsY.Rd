\name{lavResidualsY}
\alias{lavResidualsY}
\title{Compute residuals for the y-variables given the values of the x-variables}
\description{
This function can be used to compute the (case-wise) residuals for the
(observed) y-variables in a structural equation model, defined as the
difference between the observed values and the (model-based) predicted
values that are returned by \code{\link{lavPredictY}}.}
\usage{
lavResidualsY(object, newdata = NULL,
              ynames = lav_object_vnames(object, "ov.y"),
              xnames = lav_object_vnames(object, "ov.x"),
              method = "conditional.mean",
              label = TRUE, assemble = TRUE,
              force.zero.mean = FALSE,
              lambda = 0)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{newdata}{An optional data.frame, containing the same variables as
the data.frame that was used when fitting the model in \code{object}.
This data.frame should also include
the y-variables (their observed values are used to compute the residuals).
Note that if no meanstructure was used in the original fit, we will
use the saturated sample means of the original fit as substitutes for the
model-implied means. Alternatively, refit the model using \code{meanstructure = TRUE}.}
\item{ynames}{The names of the observed variables that should be treated as the
y-variables. It is for these variables that the function will compute
the (model-based) residuals for each observation. Can also be a list to allow
for a separate set of variable names per group (or block).}
\item{xnames}{The names of the observed variables that should be treated as the
x-variables. Can also be a list to allow
for a separate set of variable names per group (or block).}
\item{method}{A character string. The only available option for now is
\code{"conditional.mean"}. See Details.}
\item{label}{Logical. If TRUE, the columns of the output are labeled.}
\item{assemble}{Logical. If TRUE,
the residuals for the separate groups in the output are
reassembled into a single data.frame
with a group column, having the same dimensions as the original
(or newdata) dataset.}
\item{force.zero.mean}{Logical. Only relevant if there is no mean structure.
If \code{TRUE}, the (model-implied) mean vector is set to the zero vector.
If \code{FALSE}, the (model-implied) mean vector is set to the (unrestricted)
sample mean vector.}
\item{lambda}{Numeric. A lambda regularization penalty term.}
}
\details{
The residuals are computed as \code{observed - predicted}, where the
predicted values are obtained by \code{\link{lavPredictY}} using the same
arguments. See the help page of \code{\link{lavPredictY}} for more details
about how the predictions are computed.

These residuals can be used, for example, to diagnose the (SEM-based)
out-of-sample prediction of outcome (y) variables, or to check the
distributional assumptions (for example, multivariate normality) for the
y-variables given the x-variables.
}
\seealso{
\code{\link{lavPredictY}} to compute the (model-based) predicted values for
the y-variables given the x-variables.
}
\references{
de Rooij, M., Karch, J.D., Fokkema, M., Bakk, Z., Pratiwi, B.C, and
Kelderman, H. (2022) SEM-Based Out-of-Sample Predictions,
Structural Equation Modeling: A Multidisciplinary Journal.
\doi{10.1080/10705511.2022.2061494}
}
\examples{
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fit <- sem(model, data = PoliticalDemocracy)

lavResidualsY(fit, ynames = c("y5", "y6", "y7", "y8"),
                   xnames = c("x1", "x2", "x3", "y1", "y2", "y3", "y4"))
}
