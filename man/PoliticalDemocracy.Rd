\name{PoliticalDemocracy}
\alias{PoliticalDemocracy}
\docType{data}
\title{
Industrialization And Political Democracy Dataset
}
\description{
The `famous' Industrialization and Political Democracy dataset. This dataset is
used throughout Bollen's 1989 book (see pages 12, 17, 36 in chapter 2, pages
228 and following in chapter 7, pages 321 and following in chapter 8).
The dataset contains various measures of political democracy and
industrialization in developing countries.
}
\usage{data(PoliticalDemocracy)}
\format{
A data frame of 75 observations of 11 variables.
\describe{
    \item{\code{y1}}{Expert ratings of the freedom of the press in 1960}
    \item{\code{y2}}{The freedom of political opposition in 1960}
    \item{\code{y3}}{The fairness of elections in 1960}
    \item{\code{y4}}{The effectiveness of the elected legislature in 1960}
    \item{\code{y5}}{Expert ratings of the freedom of the press in 1965}
    \item{\code{y6}}{The freedom of political opposition in 1965}
    \item{\code{y7}}{The fairness of elections in 1965}
    \item{\code{y8}}{The effectiveness of the elected legislature in 1965}
    \item{\code{x1}}{The gross national product (GNP) per capita in 1960}
    \item{\code{x2}}{The inanimate energy consumption per capita in 1960}
    \item{\code{x3}}{The percentage of the labor force in industry in 1960}
  }

}
\source{
The dataset was retrieved from \verb{http://web.missouri.edu/~kolenikovs/Stat9370/democindus.txt} (link no longer valid; see discussion on SEMNET 18 Jun 2009)
}
\references{
Bollen, K. A. (1989). \emph{Structural Equations with Latent Variables.} 
Wiley Series in Probability and Mathematical Statistics. New York: Wiley.

Bollen, K. A. (1979). Political democracy and the timing of development.
\emph{American Sociological Review}, 44, 572-587.

Bollen, K. A. (1980). Issues in the comparative measurement of political
democracy. \emph{American Sociological Review}, 45, 370-390.
}
\examples{
head(PoliticalDemocracy)
}
