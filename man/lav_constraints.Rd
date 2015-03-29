\name{lav_constraints}
\alias{lav_constraints_parse}
\title{Utility Functions: Constraints}
\description{Utility functions for equality and inequality constraints.}
\usage{
lav_constraints_parse(partable = NULL, constraints = NULL, theta = NULL, 
                     debug = FALSE)
}
\arguments{
\item{partable}{A lavaan parameter table.}
\item{constraints}{A character string containing the constraints.}
\item{theta}{A numeric vector. Optional vector with values for the
model parameters in the parameter table.}
\item{debug}{Logical. If TRUE, show debugging information.}
}
\details{
This is a collection of lower-level constraint related functions
that are used in the lavaan code. They are made public per
request of package developers. Below is a brief description of what
they do:

The \code{lav_constraints_parse} function parses the constraints
specification (provided as a string, see examples), and generates
a list with useful information about the constraints.
}
\examples{
myModel <- 'x1 ~ a*x2 + b*x3 + c*x4'
myParTable <- lavaanify(myModel, as.data.frame. = FALSE)
con <- ' a == 2*b
         b - c == 5 '
conInfo <- lav_constraints_parse(myParTable, constraints = con)
}
