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
fitMeasures(object, fit_measures = "all",
            baseline_model = NULL, h1_model = NULL,
            fm_args = list(standard.test     = "default",
                           scaled.test       = "default",
                           rmsea.ci.level    = 0.90,
                           rmsea.close.h0    = 0.05,
                           rmsea.notclose.h0 = 0.08,
                           robust            = TRUE,
                           cat.nonpd         = "na"),
            output = "vector", level = NULL, ...)
fitmeasures(object, fit_measures = "all",
            baseline_model = NULL, h1_model = NULL,
            fm_args = list(standard.test     = "default",
                           scaled.test       = "default",
                           rmsea.ci.level    = 0.90,
                           rmsea.close.h0    = 0.05,
                           rmsea.notclose.h0 = 0.08,
                           robust            = TRUE,
                           cat.nonpd         = "na"),
            output = "vector", level = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{fit_measures}{If \code{"all"}, all fit measures available will be
  returned. If only a single or a few fit measures are specified by name,
  only those are computed and returned.
  If additional options for certain fit measures are needed, fit_measures
  should be a list whose \code{fit.measures} element holds the measures
  wanted, with the options supplied as further named elements of the list. The
  \code{standard.test} option determines the main test statistic (chi-square
  value) that will be used to compute all the fit measures that depend on this
  test statistic. Usually this is \code{"standard"}. FMG p-value tests can be
  selected here using names such as \code{"peba4"}, \code{"pols3"},
  \code{"pall"}, or \code{"all"}. The \code{scaled.test}
  option determines which scaling method is to be used for the scaled fit
  measures (in case multiple scaling methods were requested). The
  \code{rmsea.ci.level} option, default 0.9, determines the level of the confidence
  interval for the rmsea value. The \code{rmsea.close.h0} option, default 0.05,
  is the rmsea value that is used under the null hypothesis that \code{rmsea <= rmsea.close.h0}.
  The \code{rmsea.notclose.h0} option, default 0.08, is the rmsea value that is used under the
  null hypothesis that \code{rmsea >= rmsea.notclose.h0}. The
  \code{gfi.ci.level} option, default 0.9, determines the level of the
  confidence interval for the (new) goodness of fit index \code{gfi}. The
  \code{gfi} index (Maydeu-Olivares et al., 2024) is computed from the
  reweighted least squares (RLS) statistic; the version based on the
  (standard) likelihood ratio statistic is also available, and the two can
  be requested by name as \code{gfi_lrt} and \code{gfi_rls} (the latter being
  identical to \code{gfi}), together with their \code{.ci.lower},
  \code{.ci.upper} and \code{.robust} variants. The
  \code{robust}
  option, default \code{TRUE}, can be set to \code{FALSE} to avoid computing the so-called robust
  rmsea/cfi/gfi measures (for example if the computations take too long).
  The \code{cat.nonpd} option, default \code{"na"}, is only used when data is
  categorical, and determines what happens to the robust values of RMSEA, CFI
  and GFI when an input (polychoric/tetrachoric) correlation matrix is not
  positive-definite (for at least one group): if \code{"na"}, the robust
  values are \code{NA}; if \code{"refit"}, the input correlation matrix is
  smoothed to be positive-definite, and the model parameters are re-estimated
  using \code{estimator = "catML"}; if \code{"smooth"}, the input correlation
  matrix is smoothed, but the original (DWLS) parameter estimates are retained
  (this was the behavior in lavaan 0.6-13 and earlier). The older
  \code{cat.check.pd = FALSE} option is a deprecated alias for
  \code{cat.nonpd = "refit"}.}
\item{baseline_model}{If not NULL, an object of class
  \code{\linkS4class{lavaan}}, representing a user-specified baseline model.
  If a baseline model is provided, all fit indices relying on a
  baseline model (e.g., CFI or TLI) will use the test statistics from
  this user-specified baseline model, instead of the default baseline model.}
\item{h1_model}{If not NULL, an object of class \code{\linkS4class{lavaan}},
  representing a user-specified alternative to the default unrestricted model.
  If \code{h1_model} is provided, all fit indices calculated from chi-squared
  will use the chi-squared \emph{difference} test statistics from 
  \code{\link{lavTestLRT}}, which compare the user-provided \code{h1_model} 
  to \code{object}.}
\item{fm_args}{List. Additional options for certain fit measures. DEPRECATED.
  The options should now be specified in the \code{"fit_measures"} argument!}
\item{output}{Character. If \code{"vector"} (the default), display the
  output as a named (lavaan-formatted) vector. If \code{"matrix"}, display
  the output as a 1-column matrix. If \code{"text"}, display the
  output using subsections and verbose descriptions. The latter is used
  in the summary output, and does not print the chi-square test by default.
  In addition, \code{fit_measures} should contain the main ingredient
  (for example \code{"rmsea"}) if related fit measures are requested
  (for example \code{"rmsea.ci.lower"}). Otherwise, nothing will be 
  printed in that section.
  See the examples for how to add the chi-square test in the text output.}
\item{level}{Only used for multilevel models. If not \code{NULL}, either
  a level number (\code{1} or \code{2}) or a level name (\code{"within"} or
  \code{"between"}, or the level label used in the model syntax). The
  returned fit measures then reflect the fit of the model at that level
  only, while the model at the other level is saturated. This uses the
  so-called partially saturated model approach (Ryu & West, 2009): any
  misfit of the partially saturated model can be attributed to the target
  level. The incremental fit indices (e.g., CFI and TLI) are computed with
  respect to a level-specific baseline model: the independence model at
  the target level, again combined with a saturated model at the other
  level. By default, these partially saturated models are fitted (and
  stored) along with the original model; see the \code{fit.by.level}
  option in \code{\link{lavOptions}}. If they were not stored, they are
  re-fitted on the fly.}
\item{...}{Further arguments passed to or from other methods. Not currently
  used for \code{lavaan} objects.}
}
\value{
  A named numeric vector of fit measures.
}
\details{
When a scaled (or robust) test statistic is requested (for example, by using
\code{test = "satorra.bentler"}), the function will also return fit indices
based on the scaled chi-square statistic, rather than the standard version.
These scaled versions of fit measures, such as CFI and RMSEA, are calculated in
the same way as their standard counterparts, with the key difference being that
the scaled chi-square statistic is used in place of the regular one. In the
output of \code{fitMeasures()}, these appear with the \code{.scaled} suffix, 
or in the \code{Scaled} column of the \code{summary()} output.

However, this substitution-based approach---used in SEM software for many
years---has since been shown to be incorrect. Improved versions of robust fit
indices have been proposed, offering better theoretical properties. Although
still under development and not yet implemented for all estimation settings,
these improved robust fit measures are provided when available. They appear
with a \code{.robust} suffix in the output of \code{fitMeasures()}, or in the
\code{Scaled} column of the \code{summary()} output on a row labeled
\code{Robust}. As a general recommendation, these newer robust versions should
be used whenever available, in preference to the older scaled ones. See the
references below for more details.

It is also worth noting that, for models involving ordered categorical data,
robust fit indices are only computed if the underlying matrix of tetrachoric or
polychoric correlations is positive definite. If this condition is not
met---which is not uncommon in small samples---the robust measures are reported 
as \code{NA}.

Finally, in some situations (especially when the data contains missing values),
computing these robust fit indices may be computationally intensive. To avoid
long runtimes, the calculation of robust fit measures can be disabled by
setting the \code{robust} argument to \code{FALSE} in the \code{fit_measures} list.

FMG p-values can be selected through the existing \code{standard.test}
mechanism, for example
\code{fitMeasures(fit, "pvalue", fm_args = list(standard.test = "peba4"))}.
No separate FMG-specific fit-measure names are needed.

For models fitted with \code{\link{sam}} using a local \code{sam.method}
(\code{"local"}, \code{"fsr"} or \code{"cfsr"}), the fit measures are
computed for the step-2 structural model only, conditional on the (fixed)
measurement model of step 1; loglikelihood-based measures (aic, bic, ...)
are not available. If \code{sam.method = "global"}, the fit measures are
based on the joint model; since the parameters were estimated using a
two-step procedure, they should be interpreted with caution.
}
\references{
Ryu, E., & West, S. G. (2009). Level-specific evaluation of model fit in
multilevel structural equation modeling. Structural Equation Modeling, 16(4),
583--601. \doi{10.1080/10705510903203466}

Brosseau-Liard, P. E., Savalei, V., & Li, L. (2012). An investigation of the
sample performance of two nonnormality corrections for RMSEA. Multivariate
behavioral research, 47(6), 904-930.
\doi{10.1080/00273171.2012.715252}

Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit
indices for nonnormality. Multivariate behavioral research, 49(5), 460-470.
\doi{10.1080/00273171.2014.933697}

Savalei, V. (2018). On the computation of the RMSEA and CFI from the
mean-and-variance corrected test statistic with nonnormal data in SEM.
Multivariate behavioral research, 53(3), 419-429.
\doi{10.1080/00273171.2018.1455142}

Savalei, V. (2021). Improving fit indices in structural equation modeling with
categorical data. Multivariate Behavioral Research, 56(3), 390-407.
\doi{10.1080/00273171.2020.1717922}

Savalei, V., Brace, J. C., & Fouladi, R. T. (2023). We need to change how we
compute RMSEA for nested model comparisons in structural equation modeling.
Psychological Methods.
\doi{10.1037/met0000537}

Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI following
FIML and TS estimation with missing data. Psychological Methods, 28(2),
263-283. \doi{10.1037/met0000445}

Foldnes, N., Moss, J., & Gronneberg, S. (2024). Improved goodness of fit
procedures for structural equation models. \emph{Structural Equation Modeling:
A Multidisciplinary Journal}, 1-13. \doi{10.1080/10705511.2024.2372028}

Foldnes, N., Gronneberg, S., & Moss, J. (2026). Penalized eigenvalue block
averaging: Extension to nested model comparison and Monte Carlo evaluations.
\emph{Behavior Research Methods, 58}(4). \doi{10.3758/s13428-026-02968-4}
}

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
fitMeasures(fit, "pvalue", fm_args = list(standard.test = "peba4"))
## specify another threshold for RMSEA confidence interval
fitMeasures(fit, list(
  fit.measures = c("cfi", "rmsea"),
  rmsea.ci.level = 0.95))

## fit a more restricted model
fit0 <- cfa(HS.model, data = HolzingerSwineford1939, orthogonal = TRUE)
## Calculate RMSEA_D (Savalei et al., 2023)
## See https://psycnet.apa.org/doi/10.1037/met0000537
fitMeasures(fit0, "rmsea", h1_model = fit)

## level-specific fit measures for a multilevel model (Ryu & West, 2009)
model <- '
    level: 1
        fw =~ y1 + y2 + y3
        fw ~ x1 + x2 + x3
    level: 2
        fb =~ y1 + y2 + y3
        fb ~ w1 + w2
'
fit2l <- sem(model, data = Demo.twolevel, cluster = "cluster")
fitMeasures(fit2l, c("chisq", "df", "pvalue", "cfi", "rmsea",
                     "srmr_within"), level = "within")
fitMeasures(fit2l, c("chisq", "df", "pvalue", "cfi", "rmsea",
                     "srmr_between"), level = "between")
}
