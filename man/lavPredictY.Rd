\name{lavPredictY}
\alias{lavPredictY}
\title{Predict the values of y-variables given the values of x-variables}
\description{
This function can be used to predict the values of (observed) y-variables
given the values of (observed) x-variables in a structural equation model.}
\usage{
lavPredictY(object, newdata = NULL, 
            ynames = lavNames(object, "ov.y"),
            xnames = lavNames(object, "ov.x"), 
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
the y-variables (although their values will be ignored).
Note that if no meanstructure was used in the original fit, we will
use the saturated sample means of the original fit as substitutes for the
model-implied means. Alternatively, refit the model using \code{meanstructure = TRUE}.}
\item{ynames}{The names of the observed variables that should be treated as the
y-variables. It is for these variables that the function will predict 
the (model-based) values for each observation. Can also be a list to allow 
for a separate set of variable names per group (or block).}
\item{xnames}{The names of the observed variables that should be treated as the
x-variables. Can also be a list to allow 
for a separate set of variable names per group (or block).}
\item{method}{A character string. The only available option for now is
\code{"conditional.mean"}. See Details.}
\item{label}{Logical. If TRUE, the columns of the output are labeled.}
\item{assemble}{Logical. If TRUE, 
the predictions of the separate multiple groups in the output are 
reassembled again to form a single data.frame
with a group column, having the same dimensions as the original 
(or newdata) dataset.}
\item{force.zero.mean}{Logical. Only relevant if there is no mean structure.
If \code{TRUE}, the (model-implied) mean vector is set to the zero vector.
If \code{FALSE}, the (model-implied) mean vector is set to the (unrestricted)
sample mean vector.}
\item{lambda}{Numeric. A lambda regularization penalty term.}
}
\details{
This function can be used for (SEM-based) out-of-sample predictions of
outcome (y) variables, given the values of predictor (x) variables. This
is in contrast to the \code{lavPredict()} function which (historically)
only `predicts' the (factor) scores for latent variables, ignoring the
structural part of the model.

When \code{method = "conditional.mean"}, predictions (for y given x) 
are based on the (joint y and x) model-implied variance-covariance (Sigma) 
matrix and  mean vector (Mu), and the standard expression for the 
conditional mean of a multivariate normal distribution. Note that if the
model is saturated (and hence df = 0), the SEM-based predictions are identical
to ordinary least squares predictions.

Lambda is a regularization penalty term to improve prediction accuracy that can
be determined using the \code{lavPredictY_cv} function.
}
\seealso{
\code{\link{lavPredict}} to compute scores for latent variables.

\code{\link{lavPredictY_cv}} to determine an optimal lambda to increase
prediction accuracy.
}
\references{
de Rooij, M., Karch, J.D., Fokkema, M., Bakk, Z., Pratiwi, B.C, and 
Kelderman, H. (2022) SEM-Based Out-of-Sample Predictions,
Structural Equation Modeling: A Multidisciplinary Journal.
DOI:10.1080/10705511.2022.2061494

Molina, M. D., Molina, L., & Zappaterra, M. W. (2024). Aspects of Higher
Consciousness: A Psychometric Validation and Analysis of a New Model of Mystical
Experience. \doi{https://doi.org/10.31219/osf.io/cgb6e}
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

lavPredictY(fit, ynames = c("y5", "y6", "y7", "y8"),
                 xnames = c("x1", "x2", "x3", "y1", "y2", "y3", "y4"))
}
