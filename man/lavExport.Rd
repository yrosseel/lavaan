\name{lavExport}
\alias{lavExport}
\title{lavaan Export}
\description{
Export a fitted lavaan object to an external program.}
\usage{

lavExport(object, target = "lavaan", prefix = "sem", dir.name = "lavExport", 
          export = TRUE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{target}{The target program. Current options are \code{"lavaan"} and
\code{"Mplus"}.}
\item{prefix}{The prefix used to create the input files; the name of the
input file has the pattern \sQuote{prefix dot target dot in}; the name of the
data file has the pattern \sQuote{prefix dot target dot raw}.}
\item{dir.name}{The directory name (including a full path) where the 
input files will be written.}
\item{export}{If \code{TRUE}, the files are written to the output directory
(\code{dir.name}). If \code{FALSE}, only the syntax is generated as a
character string.}
}
\value{
If \code{export = TRUE}, a directory (called \code{lavExport} by default) will
be created, typically containing a data file, and an input file so that the
same analysis can be run using an external program. If \code{export = FALSE}, a
character string containing the model syntax only for the target program.
}
\details{
This function was mainly created to quickly generate an Mplus syntax file to
compare the results between Mplus and lavaan. The target \code{"lavaan"} can
be useful to create a full model syntax as needed for the \code{lavaan()}
function. More targets (perhaps for \code{LISREL} or \code{EQS}) will be added
in future releases.
}
\seealso{\code{\link{lavaanify}}, \code{\link{mplus2lavaan}}}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
out <- lavExport(fit, target = "Mplus", export=FALSE)
cat(out)
}
