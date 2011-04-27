\name{fitMeasures}
\alias{fitMeasures}
\alias{fitmeasures}
\alias{fitindices}
\title{Fit Measures for a Latent Variable Model}
\description{
This function computes a variety of fit measures to assess the global
fit of a latent variable model.}
\usage{
fitMeasures(object, fit.measures = "all")
fitmeasures(object, fit.measures = "all")
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{fit.measures}{If \code{"all"}, all fit measures available will be
returned. If only a single or a few fit measures are specified by name,
only those are computed and returned.}
}
\value{
  An named numeric vector of fit measures.
}
%\details{
%The following fit measures are available:
%
%}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
fitMeasures(fit)
fitMeasures(fit, "cfi")
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))
}
