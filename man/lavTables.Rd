\name{lavTables}
\alias{lavTables}
\title{lavaan frequency tables}
\description{
Frequency tables for categorical variables and related statistics.}
\usage{
lavTables(object, dimension = 2L, type = "cells", categorical = NULL,
    group = NULL, statistic = "default", G2.min = 3, X2.min = 3,
    p.value = FALSE, output = "data.frame", patternAsString = TRUE)
}
\arguments{
\item{object}{Either a \code{data.frame}, or an object of class 
\code{\linkS4class{lavaan}}.}
\item{dimension}{Integer. If 0L, display all response patterns. If 1L, 
display one-dimensional (one-way) tables; if 2L, display two-dimensional 
(two-way or pairwise) tables. For the latter, we can change the information
per row: if \code{type = "cells"}, each row is a cell in a pairwise table;
if \code{type = "table"}, each row is a table.}
\item{type}{If \code{"cells"}, display information for each cell in the
(one-way or two-way) table. If \code{"table"}, display information per table.
If \code{"pattern"}, display response patterns (implying 
\code{"dimension = 0L"}).}
\item{categorical}{Only used if \code{object} is a \code{data.frame}. Specify
variables that need to be treated as categorical.}
\item{group}{Only used if \code{object} is a \code{data.frame}. Specify
a grouping variable.}
\item{statistic}{Either a character string, or a vector of character strings
requesting one or more statistics for each cell, pattern or table.  Always
available are \code{X2} and \code{G2} for the Pearson and LRT based
goodness-of-fit statistics. A distinction is made between the unrestricted and
restricted model. The statistics based on the former have an extension
\code{*.un}, as in \code{X2.un} and \code{G2.un}.  If object is a
\code{data.frame}, the unrestricted versions of the statistics are the only
ones available.  For one-way tables, additional statistics are the thresholds
(\code{th.un} and \code{th}). For two-way tables and \code{type = "table"}, the
following statistics are available: \code{X2}, \code{G2}, \code{cor}
(polychoric correlation), \code{RMSEA} and the corresponding unrestricted
versions (\code{X2.un} etc). Additional statistics are \code{G2.average},
\code{G2.nlarge} and \code{G2.plarge} statistics based on the cell values 
\code{G2}:
\code{G2.average} is the average of the \code{G2} values in each cell of the
two-way table; \code{G2.nlarge} is the number of cells with a \code{G2} value
larger than \code{G2.min}, and \code{G2.plarge} is the proportion of cells with
a \code{G2} value larger than \code{G2.min}. A similar set of statistics based
on \code{X2} is also available. If \code{"default"}, the selection of
statistics (if any) depends on the \code{dim} and \code{type} arguments, and if
the object is a \code{data.frame} or a fitted lavaan object.}
\item{G2.min}{Numeric. All cells with a G2 statistic larger than this number
are considered `large', as reflected in the (optional) \code{"G2.plarge"} and
\code{"G2.nlarge"} columns.}
\item{X2.min}{Numeric. All cells with a X2 statistic larger than this number
are considered `large', as reflected in the (optional) \code{"X2.plarge"} and
\code{"X2.nlarge"} columns.}
\item{p.value}{Logical. If \code{"TRUE"}, p-values are computed for 
requested statistics (eg G2 or X2) if possible.}
\item{output}{If \code{"data.frame"}, the output is presented as a data.frame
where each row is either a cell, a table, or a response pattern, depending
on the \code{"type"} argument. If \code{"table"}, the output is presented
as a table (or matrix) or a list of tables. Only a single statistic can be
shown in this case, and if the \code{statistic} is empty, the observed 
frequencies are shown.}
\item{patternAsString}{Logical. Only used for response patterns (dimension = 0L). If \code{"TRUE"}, response patterns are displayed as a compact string. 
If \code{"FALSE"}, as many columns as observed variables are displayed.}
}
\value{
  If \code{output = "data.frame"}, the output is presented as a data.frame
  where each row is either a cell, a table, or a response pattern, depending on
  the \code{"type"} argument.
  If \code{output = "table"} (only for two-way tables),
  a list of tables (if \code{type = "cells"}) where each list element 
  corresponds to a pairwise table, or if \code{type = "table"}, a single table
  (per group). In both cases, the table entries are determined by the
  (single) \code{statistic} argument.
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

# using the data only
lavTables(HSbinary, dim = 0L, categorical = names(HSbinary))
lavTables(HSbinary, dim = 1L, categorical = names(HSbinary), stat=c("th.un"))
lavTables(HSbinary, dim = 2L, categorical = names(HSbinary), type = "table")

# fit a model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HSbinary, ordered=names(HSbinary))


lavTables(fit, 1L)
lavTables(fit, 2L, type="cells")
lavTables(fit, 2L, type="table", stat=c("cor.un", "G2", "cor"))
lavTables(fit, 2L, type="table", output="table", stat="X2")
}
