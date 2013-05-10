\name{parTable}
\alias{parameterTable}
\alias{parametertable}
\alias{parTable}
\alias{partable}
\title{Parameter Table}
\description{
Show the parameter table of a fitted model.}
\usage{
parameterTable(object)
parTable(object)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
}
\value{
  A \code{data.frame} containing the model parameters. This is
  simply the output of the \code{\link{lavaanify}} function 
  coerced to a \code{data.frame} (with \code{stringsAsFactors = FALSE}).
}
\seealso{\code{\link{lavaanify}}.}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
parTable(fit)
}
