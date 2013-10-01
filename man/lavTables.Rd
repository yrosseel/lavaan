\name{lavTables}
\alias{lavTables1D}
\alias{lavTables}
\title{lavaan Tables}
\description{
Pairwise tables for categorical variables.}
\usage{
lavTables(object, dimension = 2L, type = "cells", categorical = NULL,
    group = NULL, statistic = "default", LR.min = 3, GF.min = 3,
    p.value = FALSE, output = "data.frame", patternAsString = TRUE)
}
\arguments{
\item{object}{Either a \code{data.frame}, or an object of class 
\code{\linkS4class{lavaan}}.}
\item{dimension}{Integer. If 0L, display all response patterns. If 1L, 
display one-dimensional tables; if 2L, display two-dimensional (pairwise) 
tables.}
\item{type}{If \code{"cells"}, display information for each cell in the\
(one-way or two-way) table. If \code{"table"}, display information per
table. If \code{"pattern"}, display response patterns (implying \code{"dimension = 0L"}).}
\item{categorical}{Only used if \code{object} is a \code{data.frame}. Specify
variables that need to be treated as categorical.}
\item{group}{Only used if \code{object} is a \code{data.frame}. Specify
a grouping variable.}
\item{statistic}{Either a character string, or a vector of character
strings requesting various statistics for each cell, pattern or table.}
\item{LR.min}{Numeric. All cells with a LR statistic larger than this number
are considered `large', as reflected in the \code{"LR.plarge"} and
\code{"LR.nlarge"} columns.}
\item{GF.min}{Numeric. All cells with a GF statistic larger than this number
are considered `large', as reflected in the (optional) \code{"GF.plarge"} and
\code{"GF.nlarge"} columns.}
\item{p.value}{Logical. If \code{"TRUE"}, p-values are computed for 
requested statistics (eg LR or GF) if possible.}
\item{output}{If \code{"data.frame"}, the output is presented as a data.frame
where each row is either a cell, a table, or a response pattern, depending
on the \code{"type"} argument. If \code{"table"}, the output is presented
as a table (or matrix) or a list of tables.}
\item{patternAsString}{Logical. Only used for response patterns (dimension = 0L). If \code{"TRUE"}, response patterns are displayed as a compact string. 
If \code{"FALSE"}, as many columns as observed variables are displayed.}
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
lavTables(fit, 1L)
lavTables(fit, 2L, type="cells")
lavTables(fit, 2L, type="table")
lavTables(fit, 2L, type="table", output="table", stat="GF")
}
