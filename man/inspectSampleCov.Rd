\name{lavInspectSampleCov}
\alias{lavInspectSampleCov}
\alias{inspectSampleCov}
\title{Observed Variable Covariance Matrix from a Model and Data}
\description{
The lavaan model syntax describes a latent variable model. Often, the user wants to see the sample covariance matrix of the variables in their model for diagnostic purposes. However, their data may contain many more columns than the variables used in the model.}
\usage{
lavInspectSampleCov(model, data, ...)
inspectSampleCov(model, data, ...)
}
\arguments{
\item{model}{The model that will be fit by lavaan.}
\item{data}{The data frame being used to fit the model.}
\item{...}{Other arguments to \code{\link{sem}} for how to deal with 
  multiple groups, missing values, etc.}
}
\author{Jarrett Byrnes}
\details{
One must supply both a model, coded with proper \code{\link{model.syntax}}, and
a data frame from which a covariance matrix will be calculated. This function
essentially calls \code{\link{sem}} without fitting the model, and then uses
\code{\link{lavInspect}} to obtain the sample covariance matrix and the
meanstructure.
}

\section{See also}{\code{\link{sem}}, \code{\link{lavInspect}} }
