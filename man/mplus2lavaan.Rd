\name{mplus2lavaan}
\alias{mplus2lavaan}
\alias{lavImport}
\title{mplus to lavaan converter}
\description{
Read in an Mplus input file, convert it to lavaan syntax, and fit
the model.}
\usage{
mplus2lavaan(inpfile, run = TRUE)
}
\arguments{
\item{inpfile}{The filename (including a full path) of the Mplus input file.
The data (as referred to in the Mplus input file) should be in the same
directory as the Mplus input file.}
\item{run}{Whether to run the specified Mplus input syntax (\code{TRUE}) or
only to parse and convert the syntax (\code{FALSE}).} 
}
\value{
A \code{lavaan} object with the fitted results of the Mplus model. The parsed
and converted Mplus syntax is preserved in the \code{@external} slot of the \code{lavaan}
object in the \code{$mplus.inp} element. If \code{run} is \code{FALSE}, a \code{list} of converted
syntax is returned.
}
\author{Michael Hallquist}
\seealso{\code{\link{lavExport}}}
\examples{
\dontrun{
out <- mplus2lavaan("ex5.1.inp")
summary(out)
}
}
