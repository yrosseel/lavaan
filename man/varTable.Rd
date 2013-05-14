\name{varTable}
\alias{varTable}
\alias{vartable}
\alias{variableTable}
\alias{variabletable}
\title{Variable Table}
\description{
Summary information about the variables included in either a data.frame, 
or a fitted lavaan object.}
\usage{
varTable(object, ov.names = names(object), ov.names.x = NULL, 
    ordered = NULL, factor = NULL, as.data.frame. = TRUE)
}
\arguments{
\item{object}{Either a data.frame, or an object of class 
\code{\linkS4class{lavaan}}.}
\item{ov.names}{Only used if object is a data.frame. A character vector containing the variables that need to be summarized.}
\item{ov.names.x}{Only used if object is a data.frame. A character vector containing additional variables that need to be summarized.}
\item{ordered}{Character vector. Which variables should be treated as ordered
factors}
\item{factor}{Character vector. Which variables should be treated as (unordered) factors?}
\item{as.data.frame.}{If TRUE, return the list as a data.frame.}
}
\value{
  A \code{list} or \code{data.frame} containing summary information about
  variables in a data.frame. If \code{object} is a fitted lavaan object,
  it displays the summary information about the observed variables that are
  included in the model. The summary information includes
  variable type (numeric, ordered, \ldots), the number of non-missing values,
  the mean and variance for numeric variables, the number of levels of 
  ordered variables, and the labels for ordered variables.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
varTable(fit)
}
