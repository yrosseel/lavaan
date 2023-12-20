\name{lav_model}
\alias{lav_model_get_parameters}
\alias{lav_model_set_parameters}
\alias{lav_model_implied}
\alias{lav_model_vcov_se}
\title{lavaan model functions}
\description{Utility functions related to internal model representation (lavmodel)}
\usage{
# set/get free parameters
lav_model_set_parameters(lavmodel, x = NULL)
lav_model_get_parameters(lavmodel, GLIST = NULL, type = "free",
                         extra = TRUE)

# compute model-implied statistics
lav_model_implied(lavmodel, GLIST = NULL, delta = TRUE)

# compute standard errors
lav_model_vcov_se(lavmodel, lavpartable, VCOV = NULL, BOOT = NULL)
}
\arguments{
\item{lavmodel}{An internal representation of a lavaan model.}
\item{x}{Numeric. A vector containing the values of all the free model parameters.}
\item{GLIST}{List. A list of model matrices, similar to the output of 
\code{lavInspect(object, "est")}.}
\item{type}{Character string. If \code{"free"}, only return the free model
parameters. If \code{"user"}, return all the parameters (free and fixed) as
they appear in the user-specified parameter table.}
\item{extra}{Logical. If \code{TRUE}, also include values for rows in
the parameter table where the operator is one of \code{":="}, \code{"=="},
\code{"<"} or \code{">"}.}
\item{delta}{Logical. If \code{TRUE}, and a Delta matrix is present in GLIST,
use the (diagonal) values of the Delta matrix to rescale the covariance matrix.
This is usually needed in the categorical setting to convert covariances to
correlations.}
\item{lavpartable}{A parameter table.}
\item{VCOV}{Numeric matrix containing an estimate of the variance covariance matrix of the free model parameters.}
\item{BOOT}{Numeric matrix containing the bootstrap based parameter estimates (in the columns) for each bootstrap sample (in the rows).}
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
lavmodel <- fit@Model

est <- lav_model_get_parameters(lavmodel)
est
}

