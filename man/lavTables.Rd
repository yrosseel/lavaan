\name{lavTables}
\alias{lavTables}
\title{lavaan Tables}
\description{
Pairwise tables for categorical variables.}
\usage{
lavTables(object, categorical = NULL, std.resid = TRUE, 
          min.std.resid = 0.0, average = FALSE, collapse = FALSE) 
}
\arguments{
\item{object}{Either a \code{data.frame}, or an object of class 
\code{\linkS4class{lavaan}}.}
\item{categorical}{Only used if \code{object} is a \code{data.frame}. Specify
variables that need to be treated as categorical.}
\item{std.resid}{If \code{TRUE} and object is a lavaan object, add a column
containing the standardized residuals per cell (see equation 35 
in the reference).}
\item{min.std.resid}{Numeric. If larger than zero, show only those cells
which have a standardized residual larger than \code{min.std.resid}.}
\item{average}{If \code{TRUE}, add columns for the average standardized
residual per table (\code{str.average}), and for the number 
(\code{str.nlarge}) and proportion (\code{str.plarge}) of cells for
which the standardized residual is larger than \code{min.std.resid}.
Is set to \code{TRUE} if \code{collapse} is no longer \code{FALSE}.}
\item{collapse}{Either a logical, or the string \code{"matrix"}. If
\code{TRUE}, only one row per table is shown. If \code{"matrix"}, the
output is a symmetric matrix; the elements of the matrix are the
standardized residuals per pairwise table.}
}
\value{
  If \code{collapse = FALSE}, a \code{data.frame} where each row corresponds
  to a cell of a pairwise table. 
  If \code{collapse = TRUE}, the \code{data.frame} contains only one row
  per table. If \code{collapse = "matrix"}, a symmetric matrix 
  where the elements are the values for the standardized residual for each
  pairwise table.
}
\references{
Joreskog, K.G. & Moustaki, I. (2001). Factor analysis of ordinal variables: A
comparison of three approaches. Multivariate Behavioral Research, 36, 347-387.
}
\seealso{\code{\link{varTable}}.}
\examples{
HS9 <- HolzingerSwineford1939[,c("x1","x2","x3","x4","x5",
                                 "x6","x7","x8","x9")]
HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HSbinary, ordered=names(HSbinary))
lavTables(fit, min.std.resid=1.0)
lavTables(fit, collapse="matrix")
}
