\name{lav_use_lavaanC}
\alias{lav_use_lavaanC}
\alias{uselavaanC}
\title{Use C++ Code From The lavaanC Package In lavaan}
\description{
This function provides a global switch to switch on/off the use
of C++ code (from the lavaanC package) in lavaan.}
\usage{
lav_use_lavaanC(x)
}
\arguments{
\item{x}{Logical. Optional logical to set allow for the use of C++ code within lavaan}
}
\value{
If \code{x} is missing the current value of the switch (logical). 
If \code{x} is provided invisible(NULL).
}
\details{
The default value of the switch is ON (TRUE) if the \code{lavaanC} is 
installed and available and OFF (FALSE) if it is not. The switch cannot 
be set to TRUE if \code{lavaanC} is not available and will result in a
warning message.
}
\examples{
lav_use_lavaanC(FALSE)
lav_use_lavaanC()
lav_use_lavaanC(TRUE)
lav_use_lavaanC()
}
