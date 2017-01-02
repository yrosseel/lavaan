\name{lav_constraints}
\alias{lav_constraints_parse}
\alias{lav_partable_constraints_ceq}
\alias{lav_partable_constraints_ciq}
\alias{lav_partable_constraints_def}
\title{Utility Functions: Constraints}
\description{Utility functions for equality and inequality constraints.}
\usage{
lav_constraints_parse(partable = NULL, constraints = NULL, theta = NULL, 
                     debug = FALSE)
lav_partable_constraints_ceq(partable, con = NULL, debug = FALSE, 
                             txtOnly = FALSE)
lav_partable_constraints_ciq(partable, con = NULL, debug = FALSE, 
                             txtOnly = FALSE)
lav_partable_constraints_def(partable, con = NULL, debug = FALSE, 
                             txtOnly = FALSE)
}
\arguments{
\item{partable}{A lavaan parameter table.}
\item{constraints}{A character string containing the constraints.}
\item{theta}{A numeric vector. Optional vector with values for the
model parameters in the parameter table.}
\item{debug}{Logical. If TRUE, show debugging information.}
\item{con}{An optional partable where the operator is one of `==',
`>', `<' or `:='}
\item{txtOnly}{Logical. If TRUE, only the body of the function is returned as a character string. If FALSE, a function is returned.}
}
\details{
This is a collection of lower-level constraint related functions
that are used in the lavaan code. They are made public per
request of package developers. Below is a brief description of what
they do:

The \code{lav_constraints_parse} function parses the constraints
specification (provided as a string, see example), and generates
a list with useful information about the constraints.

The \code{lav_partable_constraints_ceq} function creates a function 
which takes the (unconstrained) parameter vector as input, and
returns the slack values for each equality constraint. If the equality
constraints hold perfectly, this function returns zeroes.

The \code{lav_partable_constraints_ciq} function creates a function
which takes the (unconstrained) parameter vector as input, and
returns the slack values for each inequality constraint. 

The \code{lav_partable_constraints_def} function creates a function
which takes the (unconstrained) parameter vector as input, and
returns the computed values of the defined parameters.
}
\examples{
myModel <- 'x1 ~ a*x2 + b*x3 + c*x4'
myParTable <- lavaanify(myModel, as.data.frame. = FALSE)
con <- ' a == 2*b
         b - c == 5 '
conInfo <- lav_constraints_parse(myParTable, constraints = con)

myModel2 <- 'x1 ~ a*x2 + b*x3 + c*x4
             a == 2*b
             b - c == 5 '
ceq <- lav_partable_constraints_ceq(partable = lavaanify(myModel2))
ceq( c(2,3,4) )
}

