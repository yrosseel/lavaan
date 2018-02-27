\name{lavPredict}
\alias{lavPredict}
\alias{lavpredict}
\title{Predict the values of latent variables (and their indicators).}
\description{
The \code{lavPredict()} function can be used to compute (or `predict') 
estimated values for latent variables, and given these values, the model-implied
values for the indicators of these latent variables.}
\usage{
lavPredict(object, type = "lv", newdata = NULL, method = "EBM",
           se.fit = FALSE, label = TRUE, fsm = FALSE, 
           optim.method = "bfgs", ETA = NULL)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{A character string. If \code{"lv"}, estimated values for
the latent variables in the model are computed. If \code{"ov"}, model
predicted values for the indicators of the latent variables in the model 
are computed.}
\item{newdata}{An optional data.frame, containing the same variables as
the data.frame used when fitting the model in object.}
\item{method}{A character string. In the linear case (when the indicators are
continuous), the possible options are \code{"regression"} or \code{"Bartlett"}.
In the categorical case, the only option (for now) is \code{"EBM"} for
the Empirical Bayes Modal approach.}
\item{se.fit}{Not used yet.}
\item{label}{Logical. If TRUE, the columns are labeled.}
\item{fsm}{Logical. If TRUE, return the factor score matrix as an attribute. 
           Only for numeric data.}
\item{optim.method}{Character string. Only used in the categorical case. 
If \code{"nlminb"} (the default in 0.5), the \code{"nlminb()"} function is used 
for the optimization. If \code{"bfgs"} or \code{"BFGS"} (the default in 0.6), 
the \code{"optim()"} function is used with the BFGS method.}
\item{ETA}{An optional matrix or list, containing latent variable values
  for each observation. Used for computations when \code{type = "ov"}.}
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
\code{\link{lavaan}}
}
\examples{
# fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
     
fit <- cfa(HS.model, data=HolzingerSwineford1939)
head(lavPredict(fit))
head(lavPredict(fit, type = "ov"))
}
