\name{varTable}
\alias{varTable}
\alias{vartable}
\alias{variableTable}
\alias{variabletable}
\title{Variable Table}
\description{
Show the variable table of the observed variables used in the model.}
\usage{
varTable(object)
vartable(object)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
}
\value{
  A \code{data.frame} containing information about the (observed) variables
  used in the model. In particular, it provides information about the
  variable type (numeric, ordered, \ldots), the number of non-missing values,
  the mean and variance for numeric variables, the number of levels of 
  ordered variables, and its labels.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
varTable(fit)
}
