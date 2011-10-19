\name{bootLavaan}
\alias{bootLavaan}
\title{Bootstrap a Lavaan Model}
\description{
Bootstraps a lavaan model using the \code{\link[boot]{boot}} function
in the boot package.}
\usage{
bootLavaan(object, data=NULL, R=1000, ..., verbose=FALSE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{data}{The data frame being used to fit the model.}
\item{R}{The number of bootstrap replicates.}
\item{...}{Other named arguments which are passed to the 
\code{\link[boot]{boot}} function in the boot package.}
\item{verbose}{If \code{TRUE}, information is shown at each bootstrap draw.}
}
\details{
\code{bootLavaan} implements the nonparameteric ordinary bootstrap, by
calling the \code{\link[boot]{boot}} function in the boot package. See
the documentation in the boot package for more information.
}
\value{
\code{bootLavaan} returns an object of class \code{lavaan.boot}, which inherits
    from class \code{boot}. A \code{summary} method will print out a data frame
    containing the original estimates, bootstrapped standard errors, and
    confidence intervals for each free parameter.
}

\seealso{\code{\link[boot]{boot}}}
