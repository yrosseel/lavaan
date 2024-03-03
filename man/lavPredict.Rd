\name{lavPredict}
\alias{lavPredict}
\alias{lavpredict}
\title{Predict the values of latent variables (and their indicators).}
\description{
The main purpose of the \code{lavPredict()} function is to compute (or
`predict') estimated values for the latent variables in the model
(`factor scores'). NOTE: the goal of this
function is NOT to predict future values of dependent variables as in the
regression framework! (For models with only continuous observed variables, the function \code{lavPredictY()} supports this.)}
\usage{
lavPredict(object, newdata = NULL, type = "lv", method = "EBM",
           transform = FALSE, se = "none", acov = "none", 
           label = TRUE, fsm = FALSE, mdist = FALSE,
           append.data = FALSE, assemble = FALSE,
           level = 1L, optim.method = "bfgs", ETA = NULL,
           drop.list.single.group = TRUE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{newdata}{An optional data.frame, containing the same variables as
the data.frame used when fitting the model in object.}
\item{type}{A character string. If \code{"lv"}, estimated values for the latent
variables in the model are computed. If \code{"ov"}, model predicted values for
the indicators of the latent variables in the model are computed. If
\code{"yhat"}, the estimated value for the observed indicators, given
user-specified values for the latent variables provided by de \code{ETA}
argument.  If \code{"fy"}, densities (or probabilities) for each observed
indicator, given user-specified values for the latent variables provided by de
\code{ETA} argument.}
\item{method}{A character string. In the linear case (when the indicators are
continuous), the possible options are \code{"regression"} or \code{"Bartlett"}.
In the categorical case, the two options are \code{"EBM"} for
the Empirical Bayes Modal approach, and \code{"ML"} for the maximum
likelihood approach.}
\item{transform}{Logical. If \code{TRUE}, transform the factor scores (per
group) so that their mean and variance-covariance matrix matches the
model-implied mean and variance-covariance matrix. This may be useful if the
individual factor scores will be used in a follow-up (regression) analysis.
Note: the standard errors (if requested) and the factor score matrix (if
requested) are not transformed (yet).}
\item{se}{Character. If \code{"none"}, no standard errors are computed.
If \code{"standard"}, naive standard errors are computed (assuming the
parameters of the measurement model are known). The standard errors are
returned as an attribute. Currently only available for complete continuous
data.}
\item{acov}{Similar to the \code{"se"} argument, but optionally returns the full
sampling covariance matrix of factor scores as an attribute. Currently 
only available for complete continuous data.}
\item{label}{Logical. If TRUE, the columns in the output are labeled.}
\item{fsm}{Logical. If TRUE, return the factor score matrix as an attribute.
           Only for numeric data.}
\item{mdist}{Logical. If TRUE, the (squared)
Mahalanobis distances of the factor scores (if \code{type = "lv"}) or
the casewise residuals (if \code{type = "resid"}) are returned as an 
attribute.}
\item{append.data}{Logical. Only used when \code{type = "lv"}. If TRUE, 
the original data (or the data provided
in the newdata argument) is appended to the factor scores.}
\item{assemble}{Logical. If TRUE, 
the separate multiple groups are reassembled again to form a single data.frame
with a group column, having the same dimensions are the original (or
newdata) dataset.}
\item{level}{Integer. Only used in a multilevel SEM.
If \code{level = 1}, only factor scores for latent variable
defined at the first (within) level are computed; if \code{level = 2},
only factor scores for latent variables defined at the second (between) level
are computed.}
\item{optim.method}{Character string. Only used in the categorical case.
If \code{"nlminb"} (the default in 0.5), the \code{"nlminb()"} function is used
for the optimization. If \code{"bfgs"} or \code{"BFGS"} (the default in 0.6),
the \code{"optim()"} function is used with the BFGS method.}
\item{ETA}{An optional matrix or list, containing latent variable values
  for each observation. Used for computations when \code{type = "ov"}.}
\item{drop.list.single.group}{Logical. If \code{FALSE}, the results are 
    returned as
    a list, where each element corresponds to a group (even if there is only
    a single group). If \code{TRUE}, the list will be unlisted if there is
    only a single group.}
}
\details{
The \code{predict()} function calls the \code{lavPredict()} function
with its default options.

If there are no latent variables in the model, \code{type = "ov"} will
simply return the values of the observed variables. Note that this function
can not be used to `predict' values of dependent variables, given the
values of independent values (in the regression sense). In other words,
the structural component is completely ignored (for now).
}
\seealso{
\code{\link{lavPredictY}} to predict y-variables given x-variables.
}
\examples{
data(HolzingerSwineford1939)

## fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
head(lavPredict(fit))
head(lavPredict(fit, type = "ov"))


## ------------------------------------------
## merge factor scores to original data.frame
## ------------------------------------------

idx <- lavInspect(fit, "case.idx")
fscores <- lavPredict(fit)
## loop over factors
for (fs in colnames(fscores)) {
  HolzingerSwineford1939[idx, fs] <- fscores[ , fs]
}
head(HolzingerSwineford1939)


## multigroup models return a list of factor scores (one per group)
data(HolzingerSwineford1939)
mgfit <- update(fit, group = "school", group.equal = c("loadings","intercepts"))

idx <- lavInspect(mgfit, "case.idx") # list: 1 vector per group
fscores <- lavPredict(mgfit)         # list: 1 matrix per group
## loop over groups and factors
for (g in seq_along(fscores)) {
  for (fs in colnames(fscores[[g]])) {
    HolzingerSwineford1939[ idx[[g]], fs] <- fscores[[g]][ , fs]
  }
}
head(HolzingerSwineford1939)

## -------------------------------------
## Use factor scores in susequent models
## -------------------------------------

## see Examples in semTools package: ?plausibleValues
}
