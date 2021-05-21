\name{lav_samplestats}
\alias{lav_samplestats_from_data}
\title{lavaan samplestats functions}
\description{Utility functions related to the sample statistics}
\usage{
# generate samplestats object from full data
lav_samplestats_from_data(lavdata = NULL, lavoptions = NULL, 
                          WLS.V = NULL, NACOV = NULL)
}
\arguments{
\item{lavdata}{A lavdata object.}
\item{lavoptions}{A named list. The Options lsot from a lavaan object.}
\item{WLS.V}{A user provided weight matrix.}
\item{NACOV}{A user provided matrix containing the elements of (N times)
    the asymptotic variance-covariance matrix of the sample statistics.
    For a multiple group analysis, a list with an asymptotic
    variance-covariance matrix for each group.}
}
\examples{
# generate syntax for an independence model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)

# extract data slot and options
lavdata <- fit@Data
lavoptions <- lavInspect(fit, "options")

# generate sample statistics object
sampleStats <- lav_samplestats_from_data(lavdata = lavdata, 
                                         lavoptions = lavoptions)
}

