\name{Demo.twolevel}
\alias{Demo.twolevel}
\docType{data}
\title{
Demo dataset for a illustrating a multilevel CFA.
}
\description{
A toy dataset containing measures on 6 items (y1-y6), 3 within-level 
covariates (x1-x3) and 2 between-level covariates (w1-w2). The data
is clustered (200 clusters of size 5, 10, 15 and 20), and the cluster 
variable is \dQuote{cluster}.
}
\usage{data(Demo.twolevel)}
\format{
A data frame of 2500 observations of 12 variables.
clusters.
\describe{
    \item{\code{y1}}{item 1}
    \item{\code{y2}}{item 2}
    \item{\code{y3}}{item 3}
    \item{\code{y4}}{item 4}
    \item{\code{y5}}{item 5}
    \item{\code{y6}}{item 6}
    \item{\code{x1}}{within-level covariate 1}
    \item{\code{x2}}{within-level covariate 2}
    \item{\code{x3}}{within-level covariate 3}
    \item{\code{w1}}{between-level covariate 1}
    \item{\code{w2}}{between-level covariate 2}
    \item{\code{cluster}}{cluster variable}
  }
}
\examples{
head(Demo.twolevel)

model <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'
fit <- sem(model, data = Demo.twolevel, cluster = "cluster")
summary(fit)
}
