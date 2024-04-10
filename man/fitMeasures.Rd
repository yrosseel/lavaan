\name{fitMeasures}
\alias{fitMeasures}
\alias{fitmeasures}
\alias{fitMeasures,lavaan-method}
\alias{fitmeasures,lavaan-method}
\alias{fitindices}
\title{Fit Measures for a Latent Variable Model}
\description{
This function computes a variety of fit measures to assess the global
fit of a latent variable model.}
\usage{
fitMeasures(object, fit.measures = "all",
            baseline.model = NULL, h1.model = NULL,
            fm.args = list(standard.test     = "default",
                           scaled.test       = "default",
                           rmsea.ci.level    = 0.90,
                           rmsea.close.h0    = 0.05,
                           rmsea.notclose.h0 = 0.08,
                           robust            = TRUE,
                           cat.check.pd      = TRUE),
            output = "vector", ...)
fitmeasures(object, fit.measures = "all",
            baseline.model = NULL, h1.model = NULL,
            fm.args = list(standard.test     = "default",
                           scaled.test       = "default",
                           rmsea.ci.level    = 0.90,
                           rmsea.close.h0    = 0.05,
                           rmsea.notclose.h0 = 0.08,
                           robust            = TRUE,
                           cat.check.pd      = TRUE),
            output = "vector", ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{fit.measures}{If \code{"all"}, all fit measures available will be
  returned. If only a single or a few fit measures are specified by name,
  only those are computed and returned.}
\item{baseline.model}{If not NULL, an object of class
  \code{\linkS4class{lavaan}}, representing a user-specified baseline model.
  If a baseline model is provided, all fit indices relying on a
  baseline model (eg. CFI or TLI) will use the test statistics from
  this user-specified baseline model, instead of the default baseline model.}
\item{h1.model}{If not NULL, an object of class \code{\linkS4class{lavaan}},
  representing a user-specified alternative to the default unrestricted model.
  If \code{h1.model} is provided, all fit indices calculated from chi-squared
  will use the chi-squared \emph{difference} test statistics from 
  \code{\link{lavTestLRT}}, which compare the user-provided \code{h1.model} 
  to \code{object}.}
\item{fm.args}{List. Additional options for certain fit measures. The
  \code{standard.test} element determines the main test statistic (chi-square
  value) that will be used to compute all the fit measures that depend on this
  test statistic. Usually this is \code{"standard"}. The \code{scaled.test}
  element determines which scaling method is to be used for the scaled fit
  measures (in case multiple scaling methods were requested). The
  \code{rmsea.ci.level} element determines the level of the confidence
  interval for the rmsea value. The \code{rmsea.close.h0} element
  is the rmsea value
  that is used under the null hypothesis that \code{rmsea <= rmsea.close.h0}.
  The \code{rmsea.notclose.h0} element is the rmsea value that is used under the
  null hypothesis that \code{rsmsea >= rmsea.notclose.h0}. The \code{robust} 
  element can be set to \code{FALSE} to avoid computing the so-called robust
  rmsea/cfi measures (for example if the computations take too long).
  The \code{cat.check.pd} element is only used when data is categorical. If
  \code{TRUE}, robust values for RMSEA and CFI are only computed if the input
  correlation matrix is positive-definite (for all groups).}
\item{output}{Character. If \code{"vector"} (the default), display the
  output as a named (lavaan-formatted) vector. If \code{"matrix"}, display
  the output as a 1-column matrix. If \code{"text"}, display the
  output using subsections and verbose descriptions. The latter is used
  in the summary output, and does not print the chi-square test by default.
  In addition, \code{fit.measures} should contain the main ingredient
  (for example \code{"rmsea"}) if related fit measures are requested
  (for example \code{"rmsea.ci.lower"}). Otherwise, nothing will be 
  printed in that section.
  See the examples how to add the chi-square test in the text output.}
\item{...}{Further arguments passed to or from other methods. Not currently
  used for \code{lavaan} objects.}
}
\value{
  A named numeric vector of fit measures.
}
%\details{
%The following fit measures are available:
%
%}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
fitMeasures(fit)
fitMeasures(fit, "cfi")
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea"))
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea"), 
            output = "matrix")
fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "rmsea"),
            output = "text")

## fit a more restricted model
fit0 <- cfa(HS.model, data = HolzingerSwineford1939, orthogonal = TRUE)
## Calculate RMSEA_D (Savalei et al., 2023)
## See https://psycnet.apa.org/doi/10.1037/met0000537
fitMeasures(fit0, "rmsea", h1.model = fit)
}
