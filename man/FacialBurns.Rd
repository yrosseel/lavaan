\name{FacialBurns}
\alias{FacialBurns}
\docType{data}
\title{Dataset for illustrating the InformativeTesting function.}

\description{
A dataset from the Dutch burn center (http://www.adbc.nl). 
The data were used to examine psychosocial functioning in patients with 
facial burn wounds. Psychosocial functioning was measured by 
Anxiety and depression symptoms (HADS), and self-esteem 
(Rosenberg's self-esteem scale).} 

\usage{data(FacialBurns)}

\format{
A data frame of 77 observations of 6 variables.
\describe{
    \item{\code{Selfesteem}}{Rosenberg's self-esteem scale}
    \item{\code{HADS}}{Anxiety and depression scale}
    \item{\code{Age}}{Age measured in years, control variable}
    \item{\code{TBSA}}{Total Burned Surface Area}
    \item{\code{RUM}}{Rumination, control variable}
    \item{\code{Sex}}{Gender, grouping variable}
  }

}
\examples{
head(FacialBurns)
}
