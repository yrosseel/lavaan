\name{lav_use_lavaanC}
\alias{lav_use_lavaanC}
\alias{uselavaanC}
\title{use C++ code in lavaanC}
\description{
Get or set switch to use the C++ code in lavaanC.}
\usage{
lav_use_lavaanC(x)
}
\arguments{
\item{x}{Optional scalar logical to set the switch ON or OFF.}
}
\value{
If \code{x} is missing the current value of the switch (logical). 
If \code{x} is provided invisible(NULL).
}
\details{
The default value of the switch is ON (TRUE) if the \code{lavaanC} is available and
OFF (FALSE) if it is not. The switch cannot be set to TRUE if \code{lavaanC} is not available
(a warning message is issued).
}
\examples{
lav_use_lavaanC(FALSE)
lav_use_lavaanC()
lav_use_lavaanC(TRUE)
lav_use_lavaanC()
}
