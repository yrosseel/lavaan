\name{mplus2lavaan.modelSyntax}
\alias{mplus2lavaan.modelSyntax}
\title{Convert Mplus model syntax to lavaan}
\description{
Converts Mplus model syntax into lavaan model syntax.}
\usage{
mplus2lavaan.modelSyntax(syntax)
}
\arguments{
\item{syntax}{A character vector containing Mplus model syntax to be
  converted to lavaan model syntax. Note that parsing Mplus syntax often
  requires correct usage of newline characters. If \code{syntax} is a
  vector of multiple strings, these will be joined with newlines prior
  to conversion. Alternatively, \code{\\n} characters can be included
  inline in \code{syntax}.}  
}
\value{
A character string of converted \code{lavaan} model syntax.
}
\author{Michael Hallquist}
\seealso{\code{\link{mplus2lavaan}}}
\examples{
\dontrun{
syntax <- '
    f1 BY x1*1 x2 x3;
    x1 WITH x2;
    x3 (1);
    x2 (1);
'
lavSyntax <- mplus2lavaan.modelSyntax(syntax)
cat(lavSyntax)
}
}
