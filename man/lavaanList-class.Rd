\name{lavaanList-class}
\docType{class}
\alias{lavaanList-class}
\alias{coef,lavaanList-method}
\alias{summary,lavaanList-method}

\title{Class For Representing A List of (Fitted) Latent Variable Models}
\description{The \code{lavaanList} class represents a collection of (fitted)
  latent variable models, for a (potentially large) number of datasets. 
  It contains information about the model (which is always the same),
  and for every dataset a set of (user-specified) slots from a regular
  lavaan object.}
\section{Objects from the Class}{
Objects can be created via the
  \code{\link{cfaList}}, \code{\link{semList}}, or
  \code{\link{lavaanList}} functions.
}
\section{Slots}{
  \describe{
    \item{\code{version}:}{The lavaan package version used to create this objects}
    \item{\code{call}:}{The function call as returned by \code{match.call()}.}
    \item{\code{Options}:}{Named list of options that were provided by
      the user, or filled-in automatically.}
    \item{\code{ParTable}:}{Named list describing the model parameters. Can be coerced to a data.frame. In the documentation, this is called the `parameter table'.}
    \item{\code{pta}:}{Named list containing parameter table attributes.}
    \item{\code{Data}:}{Object of internal class \code{"Data"}: information
about the data.}
    \item{\code{Model}:}{Object of internal class \code{"Model"}: the 
      internal (matrix) representation of the model}
    \item{\code{meta}:}{List containing additional flags. For internal use only.}
    \item{\code{timingList}:}{List. Timing slot per dataset.}
    \item{\code{ParTableList}:}{List. ParTable slot per dataset.}
    \item{\code{DataList}:}{List. Data slot per dataset.}
    \item{\code{SampleStatsList}:}{List. SampleStats slot per dataset.}
    \item{\code{CacheList}:}{List. Cache slot per dataset.}
    \item{\code{vcovList}:}{List. vcov slot per dataset.}
    \item{\code{testList}:}{List. test slot per dataset.}
    \item{\code{optimList}:}{List. optim slot per dataset.}
    \item{\code{impliedList}:}{List. implied slot per dataset.}
    \item{\code{h1List}:}{List. h1 slot per dataset.}
    \item{\code{loglikList}:}{List. loglik slot per dataset.}
    \item{\code{baselineList}:}{List. baseline slot per dataset.}
    \item{\code{funList}:}{List. fun slot per dataset.}
    \item{\code{internalList}:}{List. internal slot per dataset.}
    \item{\code{external}:}{List. Empty slot to be used by add-on packages.}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "lavaanList", type = "free")}: Returns 
      the estimates of the parameters in the model as the columns in a matrix;
      each column corresponds to a different dataset.
      If \code{type="free"}, only the free parameters are returned.
      If \code{type="user"}, all parameters listed in the parameter table
      are returned, including constrained and fixed parameters.}
    \item{summary}{\code{signature(object = "lavaanList", header = TRUE, estimates = TRUE, nd = 3L)}: 
      Print a summary of the collection of fitted models.
      If \code{header = TRUE}, the header section is
      printed.  
      If \code{estimates = TRUE}, print the parameter estimates section.
      The argument \code{nd} determines the number of digits after the
      decimal point to be printed (currently only in the parameter estimates
      section.)
      Nothing is returned (use 
      \code{parameterEstimates} or another extractor function 
      to extract information from this object).}
  }
}
\seealso{
\code{\link{cfaList}}, \code{\link{semList}}, \code{\link{lavaanList}}
}
