\name{mplus2lavaan}
\alias{mplus2lavaan}
\alias{lavImport}
\title{mplus to lavaan converter}
\description{
Read in an Mplus input file, convert it to lavaan syntax, and fit
the model.}
\usage{
mplus2lavaan(inpfile)
}
\arguments{
\item{inpfile}{The filename (including a full path) of the Mplus input file.
The data (as referred to in the Mplus input file) should be in the same
directory as the Mplus input file.}
}
\value{
A list with two elements: \code{mplus.inp} contains the input data, a title,
the variable names, and the converted (lavaan) model syntax; \code{lav.out}
contains the fitted lavaan object.
}
\author{Michael Hallquist}
\seealso{\code{\link{lavExport}}.}
\examples{
\dontrun{
out <- mplus2lavaan("ex5.1.inp")
summary(out$lav.out)
}
}
