\name{mimic}
\alias{mimic}
\title{The mimic= option}
\description{
The \code{mimic=} argument of the \code{\link{lavaan}} function (and its
wrappers \code{\link{cfa}}, \code{\link{sem}} and \code{\link{growth}})
controls a collection of primarily technical default settings, with the aim of
making lavaan's output as similar as possible to that of other SEM software.
This applies in particular to a few key results, such as the value of the
Satorra-Bentler scaled test statistic.

The \code{mimic=} option was originally introduced to demonstrate that
differences between lavaan and other SEM programs were often due to technical
implementation details rather than errors in lavaan. Nevertheless, the scope of
the \code{mimic=} option is intentionally limited. It is designed to reproduce
a small number of important results and does not attempt to replicate every
quantity reported by alternative SEM software packages.

This help page documents which options are affected by the
\code{mimic=} argument.
}
\details{
The \code{mimic=} argument does not change the model that is fitted, nor the
estimator that is used. Instead, it overwrites the default values of a number
of (mostly technical) options that are otherwise estimator-dependent. The
user can always override these defaults by setting the corresponding option
explicitly in the \code{\link{lavaan}} call.

The argument is case-insensitive and accepts a number of aliases. After
normalization, only four values remain:
\describe{
  \item{\code{"lavaan"}:}{the native lavaan defaults (also selected by the
    value \code{"default"}); this is the default.}
  \item{\code{"Mplus"}:}{mimic the \code{Mplus} program (also selected by
    the value \code{"mplus"}).}
  \item{\code{"EQS"}:}{mimic the \code{EQS} program (also selected by the
    values \code{"eqs"}, \code{"LISREL"} and \code{"lisrel"}).}
  \item{\code{"lm"}:}{a minimal setting used internally for plain regression
    (also selected by the value \code{"regression"}); not intended for
    general use.}
}

\strong{mimic = "lavaan" (the default)}

This is the reference setting. For maximum likelihood (ML) estimation, the
following defaults are used:
\describe{
  \item{\code{likelihood = "normal"}:}{the sample covariance matrix is
    rescaled by a factor \code{(N-1)/N}, so the loglikelihood and standard
    errors are based on a division by \code{N} (rather than \code{N-1}).}
  \item{\code{conditional.x = FALSE}:}{the exogenous \code{x} covariates are
    regressed out first, and a conditional (on x) model is fitted.}
  \item{\code{fixed.x = TRUE}:}{the means, variances and covariances of the
    exogenous \code{x} covariates are fixed to their sample values (for the
    ML, MML and IV estimator groups, unless \code{start = "simple"}).}
  \item{\code{zero.keep.margins = TRUE}:}{for categorical data, the (one-way)
    margins are preserved when adding a small constant to empty cells.}
}
In addition, the various Mplus-specific technical options
(\code{information.expected.mplus}, \code{gamma.vcov.mplus},
\code{gamma.wls.mplus}, \code{gls.v11.mplus} and
\code{rmsea.scaled.mplus}) are all set to \code{FALSE}.

\strong{mimic = "Mplus"}

Starting from the lavaan defaults, the following options are changed in an
attempt to reproduce the output of \code{Mplus}:
\describe{
  \item{\code{information.expected.mplus = TRUE}:}{use the Mplus variant of
    the expected information matrix. This affects the standard errors (when
    \code{information = "expected"}) and the scaled (Satorra-Bentler and
    Yuan-Bentler) test statistics.}
  \item{\code{gamma.vcov.mplus = TRUE}:}{use the Mplus way of computing the
    \eqn{\Gamma} matrix (the asymptotic covariance matrix of the sample
    statistics) that enters the robust (sandwich) standard errors.}
  \item{\code{gamma.wls.mplus = TRUE}:}{use the Mplus way of computing the
    \eqn{\Gamma} matrix used as the weight matrix in (D)WLS estimation.}
  \item{\code{gls.v11.mplus = TRUE}:}{use the Mplus variant of the weight
    matrix used in GLS estimation (when there is a meanstructure and
    \code{conditional.x = FALSE}).}
  \item{\code{rmsea.scaled.mplus = TRUE}:}{use the Mplus reference
    distribution for the \emph{scaled} RMSEA confidence interval and
    p-values. lavaan refers the \emph{unscaled} test statistic to a
    chi-square distribution with \eqn{tr[U\Gamma]} degrees of freedom, while
    Mplus refers the \emph{scaled} statistic to a chi-square distribution
    with the degrees of freedom of the scaled test. Both match the mean of
    the null distribution, but on a different scale, so they yield different
    intervals and p-values. This only affects the Satorra-Bentler and
    Yuan-Bentler family (\code{MLM}, \code{MLR}, \code{MLMVS},
    \code{WLSM}, ...); for the shifted family (\code{MLMV}, \code{WLSMV},
    \code{ULSMV}) lavaan already uses the Mplus convention.}
  \item{\code{missing = "ml"}:}{for continuous data and the ML or MLR
    estimator, the default missing-data handling becomes full information
    maximum likelihood (FIML), rather than listwise deletion.}
  \item{\code{meanstructure = TRUE}:}{a meanstructure is added by default
    (unless the estimator is PML).}
  \item{\code{group.equal}:}{in a multiple-group analysis, if no equality
    constraints are requested, a default set is imposed: \code{"loadings"}
    and \code{"thresholds"} for categorical data;
    \code{"loadings"} for continuous data without a meanstructure; and
    \code{"loadings"} and \code{"intercepts"} for continuous data with a
    meanstructure.}
  \item{\code{baseline.conditional.x.free.slopes = FALSE}:}{when
    \code{conditional.x = TRUE}, the baseline (independence) model does not
    free the slope structure.}
}
The \code{likelihood}, \code{conditional.x}, \code{fixed.x} and
\code{zero.keep.margins} defaults are the same as for
\code{mimic = "lavaan"}.

\strong{mimic = "EQS"}

Starting from the lavaan defaults, the following options are changed in an
attempt to reproduce the output of \code{EQS} (the value \code{"LISREL"} is
treated as an alias for \code{"EQS"}):
\describe{
  \item{\code{likelihood = "wishart"}:}{for ML estimation, the Wishart
    likelihood is used. The sample covariance matrix is divided by
    \code{N-1} (instead of \code{N}), which affects the chi-square test
    statistic, the loglikelihood and the standard errors.}
  \item{\code{baseline.fixed.x.free.cov = FALSE}:}{when \code{fixed.x = TRUE},
    the baseline (independence) model does not free the (co)variances of the
    exogenous covariates.}
}
In addition, for the MLR estimator, the robust test statistic defaults to
\code{"yuan.bentler"} (instead of the \code{"yuan.bentler.mplus"} variant
that is used for \code{mimic = "lavaan"} and \code{mimic = "Mplus"}).

Note that the \code{mimic=} argument only sets \emph{defaults}: any option
that is set explicitly in the \code{\link{lavaan}} call takes precedence. The
mimicking is also necessarily incomplete; small numerical differences with
the target program may remain.
}
\seealso{
\code{\link{lavOptions}}, \code{\link{lavaan}}
}
\examples{
# default (mimic = "lavaan")
fit1 <- cfa(" visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 ",
            test = "satorra.bentler",
            missing = "listwise",
            data = HolzingerSwineford1939)
fit1

# mimic the Mplus program
fit2 <- cfa(" visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 ",
            test = "satorra.bentler",
            missing = "listwise",
            data = HolzingerSwineford1939, mimic = "Mplus")
fit2

# mimic the EQS program
fit3 <- cfa(" visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 ",
            test = "satorra.bentler",
            missing = "listwise",
            data = HolzingerSwineford1939, mimic = "EQS")
fit3
}
