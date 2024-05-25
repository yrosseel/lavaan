\name{lavPredictY_cv}
\alias{lavPredictY_cv}
\title{Determine an optimal lambda penalty value through cross-validation}
\description{
This function can be used to determine an optimal lambda value for the
\code{lavPredictY} function. based on cross-validation.}
\usage{
lavPredictY_cv(object, data = NULL,
              xnames = lavNames(object, "ov.x"),
              ynames = lavNames(object, "ov.y"),
              n.folds = 10L,
              lambda.seq = seq(0, 1, 0.1))
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{data}{A data.frame, containing the same variables as the data.frame that
was used when fitting the model in \code{object}.} 
\item{xnames}{The names of the observed variables that should be treated as the
x-variables. Can also be a list to allow for a separate set of variable names
per group (or block).}
\item{ynames}{The names of the observed variables that should be treated as the
y-variables. It is for these variables that the function will predict the
(model-based) values for each observation. Can also be a list to allow for a
separate set of variable names per group (or block).}
\item{n.folds}{Integer. The number of folds to be used during
cross-validation.}
\item{lambda.seq}{An R \code{seq()} containing the range of lambda penalty
values to be tested during cross-validation.}
}
\details{
This function is used to generate an optimal lambda value for
\code{lavPredictY} predictions to improve prediction accuracy.
}
\seealso{
\code{\link{lavPredictY}} to predict the values of (observed) y-variables given
the values of (observed) x-variables in a structural equation model.
}
\references{
de Rooij, M., Karch, J.D., Fokkema, M., Bakk, Z., Pratiwi, B.C, and 
Kelderman, H. (2022) SEM-Based Out-of-Sample Predictions,
Structural Equation Modeling: A Multidisciplinary Journal.
DOI:10.1080/10705511.2022.2061494

Molina, M. D., Molina, L., & Zappaterra, M. W. (2024). Aspects of Higher
Consciousness: A Psychometric Validation and Analysis of a New Model of
Mystical Experience. \doi{https://doi.org/10.31219/osf.io/cgb6e}
}
\examples{
colnames(PoliticalDemocracy) <- c("z1", "z2", "z3", "z4", 
                                  "y1", "y2", "y3", "y4", 
                                  "x1", "x2", "x3")

model <- '
  # latent variable definitions
  ind60 =~ x1 + x2 + x3
  dem60 =~ z1 + z2 + z3 + z4
  dem65 =~ y1 + y2 + y3 + y4
  # regressions
  dem60 ~ ind60
  dem65 ~ ind60 + dem60
  # residual correlations
  z1 ~~ y1
  z2 ~~ z4 + y2
  z3 ~~ y3
  z4 ~~ y4
  y2 ~~ y4
'
fit <- sem(model, data = PoliticalDemocracy, meanstructure = TRUE)

percent <- 0.5
nobs <- lavInspect(fit, "ntotal")
idx <- sort(sample(x = nobs, size = floor(percent*nobs)))

xnames = c("z1", "z2", "z3", "z4", "x1", "x2", "x3")
ynames = c("y1", "y2", "y3", "y4")

reg.results <- lavPredictY_cv(
    fit,
    PoliticalDemocracy[-idx, ],
    xnames = xnames,
    ynames = ynames,
    n.folds = 10L,
    lambda.seq = seq(from = .6, to = 2.5, by = .1)
)
lam <- reg.results$lambda.min

lavPredictY(fit, newdata = PoliticalDemocracy[idx,],
                 ynames  = ynames,
                 xnames  = xnames,
                 lambda  = lam)
}
