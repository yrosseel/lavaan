\name{cvLavPredictY}
\alias{cvLavPredictY}
\title{Determine an optimal lambda penalty value through cross-validation}
\description{
This function can be used to determine an optimal lambda value for the
\code{lavPredictY} function. based on cross-validation.}
\usage{
cvLavPredictY(model, data, 
            ynames = lavNames(object, "ov.y"),
            xnames = lavNames(object, "ov.x"),
            n.folds = 10,
            lambda.seq)
}
\arguments{
\item{model}{A description of the user-specified model using the lavaan model syntax. }
\item{data}{A data.frame, containing the same variables as
the data.frame that was used when fitting the model in \code{object}. 
This data.frame should also include the y-variables.}
\item{ynames}{The names of the observed variables that should be treated as the
y-variables. It is for these variables that the function will predict 
the (model-based) values for each observation. Can also be a list to allow 
for a separate set of variable names per group (or block).}
\item{xnames}{The names of the observed variables that should be treated as the
x-variables. Can also be a list to allow 
for a separate set of variable names per group (or block).}
\item{n.folds}{Numeric. The number of folds to be used during cross-validation.}
\item{lambda.seq}{An R \code{seq()} containing the range of lambda penalty
values to be tested during cross-validation.}
}
\details{
This function is used to generate an optimal lambda value for \code{lavPredictY}
predictions to improve prediction accuracy.
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

Molina, M. D., Molina, L., & Zappaterra, M. W. (2024). Aspects of Higher Consciousness: A Psychometric Validation and Analysis of a New Model of Mystical Experience. DOI:
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
data <- PoliticalDemocracy
colnames(data) = c("z1", "z2", "z3", "z4", "y1", "y2", "y3", "y4", "x1", "x2", "x3")
fit <- sem(model, data = data)

percent <- .9
idx <- sort(sample(x = nrow(data), size = floor(percent*nrow(data)), replace = F))

reg.results <- cvLavPredictY(
    model,
    PoliticalDemocracy[-idx, ],
    xnames,
    ynames,
    n.folds = 10,
    lambda.seq = seq(from = .6, to = 2.5, by = .1)
)
lam <- reg.results$lambda.min

lavPredictY(fit, newdata = PoliticalDemocracy[idx,]
                  ynames = colnames(PoliticalDemocracy)[-c(5,6,7,8)],
                  xnames = colnames(PoliticalDemocracy)[c(5,6,7,8)],
                  lambda = lam)
}
