\name{Demo.growth}
\alias{Demo.growth}
\docType{data}
\title{
Demo dataset for a illustrating a linear growth model.
}
\description{
A toy dataset containing measures on 4 time points (t1,t2, t3 and t4),
two predictors (x1 and x2) influencing the random intercept and slope, and 
a time-varying covariate (c1, c2, c3 and c4).
}
\usage{data(Demo.growth)}
\format{
A data frame of 400 observations of 10 variables.
\describe{
    \item{\code{t1}}{Measured value at time point 1}
    \item{\code{t2}}{Measured value at time point 2}
    \item{\code{t3}}{Measured value at time point 3}
    \item{\code{t4}}{Measured value at time point 4}
    \item{\code{x1}}{Predictor 1 influencing intercept and slope}
    \item{\code{x2}}{Predictor 2 influencing intercept and slope}
    \item{\code{c1}}{Time-varying covariate time point 1}
    \item{\code{c2}}{Time-varying covariate time point 2}
    \item{\code{c3}}{Time-varying covariate time point 3}
    \item{\code{c4}}{Time-varying covariate time point 4}
  }

}
\seealso{
\code{\link{growth}}
}
\examples{
head(Demo.growth)
}
