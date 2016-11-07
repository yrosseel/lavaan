\name{HolzingerSwineford1939}
\alias{HolzingerSwineford1939}
\docType{data}
\title{
Holzinger and Swineford Dataset (9 Variables)
}
\description{
The classic Holzinger and Swineford (1939) dataset consists of mental
ability test scores of seventh- and eighth-grade children from two
different schools (Pasteur and Grant-White). In the original dataset
(available in the \code{MBESS} package), there are scores for 26 tests.
However, a smaller subset with 9 variables is more widely used in the 
literature (for example in Joreskog's 1969 paper, which also uses the 145 
subjects from the Grant-White school only). 
}
\usage{data(HolzingerSwineford1939)}
\format{
  A data frame with 301 observations of 15 variables.
  \describe{
    \item{\code{id}}{Identifier}
    \item{\code{sex}}{Gender}
    \item{\code{ageyr}}{Age, year part}
    \item{\code{agemo}}{Age, month part}
    \item{\code{school}}{School (Pasteur or Grant-White)}
    \item{\code{grade}}{Grade}
    \item{\code{x1}}{Visual perception}
    \item{\code{x2}}{Cubes}
    \item{\code{x3}}{Lozenges}
    \item{\code{x4}}{Paragraph comprehension}
    \item{\code{x5}}{Sentence completion}
    \item{\code{x6}}{Word meaning}
    \item{\code{x7}}{Speeded addition}
    \item{\code{x8}}{Speeded counting of dots}
    \item{\code{x9}}{Speeded discrimination straight and curved capitals}
  }
}
\source{
This dataset was originally retrieved from \verb{http://web.missouri.edu/~kolenikovs/stata/hs-cfa.dta} (link no longer active) and converted to an R dataset.
}
\references{
Holzinger, K., and Swineford, F. (1939). A study in factor analysis: The 
stability of a bifactor solution. Supplementary Educational Monograph, no. 48. 
Chicago: University of Chicago Press. 

Joreskog, K. G. (1969). A general approach to confirmatory maximum likelihood
factor analysis. \emph{Psychometrika}, 34, 183-202.
}
\seealso{
\code{\link{cfa}}
}
\examples{
head(HolzingerSwineford1939)
}
