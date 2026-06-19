\name{estimator_iv}
\alias{estimator_iv}
\alias{iv}
\alias{IV}
\alias{miiv}
\alias{MIIV}
\alias{miivsem}
\alias{MIIVsem}
\alias{instrument}
\alias{instruments}
\alias{sargan}
\alias{Sargan}
\title{Instrumental Variable (MIIV-2SLS) Estimation}
\description{
Estimate a structural equation model using model-implied instrumental
variables and two-stage least squares (MIIV-2SLS), a non-iterative,
equation-by-equation alternative to maximum likelihood (Bollen, 1996).
This page documents how to invoke the estimator, the options that can be
passed via the \code{estimator = list()} argument, how to specify
instruments manually with the \code{|~} operator, and the \code{lavInspect()}
fields that are relevant for IV estimation.
}
\details{
Instead of minimizing a single discrepancy function for the whole model,
MIIV-2SLS estimates each measurement and structural equation separately by
two-stage least squares (2SLS). For each equation, the instruments are other
observed variables in the model that are (under the hypothesized model)
uncorrelated with the equation's composite error term: the \emph{model-implied
instrumental variables} (MIIVs). Because the equations are estimated
separately, a misspecification in one equation does not necessarily propagate
to the others, the method is non-iterative, and no starting values are
required.

The estimator is selected by setting \code{estimator = "IV"} (the alias
\code{"MIIV"} is also accepted; names are case-insensitive), for example
\preformatted{    fit <- sem(model, data = Data, estimator = "IV")}

Estimation proceeds in two stages. In the first stage, the \emph{directed}
effects (factor loadings and regression coefficients) are estimated
equation-by-equation by 2SLS. In the second stage, the \emph{undirected}
effects (residual variances and covariances) are estimated by a weighted
least-squares step applied to the residual moments. Standard errors are
available for both stages (see the options below).

Both continuous and ordered-categorical data are supported (for ordered data,
the polychoric-based PIV estimator of Bollen & Maydeu-Olivares, 2007 is used).
Simple equality constraints (e.g. equal factor loadings) and general linear
equality constraints (e.g. \code{a == 2*b}) are honored. For now, the IV
estimator is restricted to single-group, single-level models.

\subsection{Estimator options}{

Options specific to the IV estimator are passed as named elements of a list
to the \code{estimator} argument, alongside the estimator name:
\preformatted{    fit <- sem(model, data = Data,
               estimator = list(estimator = "IV",
                                iv_varcov_method = "2RLS",
                                iv_sargan_adjust = "BH"))}

The option names use snake_case. For backward compatibility, dot.case names
(e.g. \code{iv.varcov.method}) are also accepted and silently converted.

\describe{
  \item{\code{iv_method}:}{Character. The instrumental-variable method.
    Currently only \code{"2SLS"} (two-stage least squares) is available.
    Default is \code{"2SLS"}.}
  \item{\code{iv_samplestats}:}{Logical. If \code{TRUE} (the default), the
    equations are estimated from the sample moments (covariances and means)
    rather than from the raw data. This is required for categorical data and
    when only sample moments are provided as input.}
  \item{\code{iv_varcov_method}:}{Character. The second-stage method used to
    estimate the residual variances and covariances. One of \code{"ULS"},
    \code{"GLS"}, \code{"2RLS"}, \code{"RLS"} or \code{"NONE"}. Default is
    \code{"RLS"} (reweighted least squares). For categorical data, \code{"ULS"}
    is used. If \code{"NONE"}, the variance/covariance parameters are not
    estimated.}
  \item{\code{iv_sargan}:}{Logical. If \code{TRUE} (the default), a Sargan
    overidentification test is computed for each overidentified equation
    (Sargan, 1958). The results can be retrieved with
    \code{lavInspect(fit, "sargan")}.}
  \item{\code{iv_sargan_adjust}:}{Character. A multiple-comparison adjustment
    applied to the per-equation Sargan p-values, using any method accepted by
    \code{\link[stats]{p.adjust}} (\code{"holm"}, \code{"hochberg"},
    \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"},
    \code{"fdr"}, or \code{"none"}). Default is \code{"none"}. When a method
    other than \code{"none"} is requested, an extra \code{pvalue.adj} column is
    added to the table returned by \code{lavInspect(fit, "sargan")}.}
  \item{\code{iv_weak}:}{Character. How to respond to weak instruments,
    diagnosed per equation by the first-stage F-statistic (Staiger & Stock,
    1997). One of:
    \describe{
      \item{\code{"warn"}:}{(the default) report the affected equations and
        their first-stage F, and suggest supplying stronger instruments
        manually; the estimates are not changed.}
      \item{\code{"prune"}:}{greedily drop the weakest excess instruments from
        weak, overidentified equations (never below just-identified) and report
        what was dropped.}
      \item{\code{"none"}:}{skip the weak-instrument check.}
    }}
  \item{\code{iv_weak_threshold}:}{Numeric. The first-stage F-statistic below
    which instruments are deemed weak. Default is \code{10}.}
  \item{\code{iv_vcov_stage1}:}{Character. The standard-error method for the
    first-stage (directed) coefficients. One of \code{"lm.vcov"},
    \code{"lm.vcov.dfres"}, \code{"gamma"} or \code{"none"}. Default is
    \code{"lm.vcov.dfres"} for continuous data and \code{"gamma"} for
    categorical data. The \code{"gamma"} method requires
    \code{iv_samplestats = TRUE}.}
  \item{\code{iv_vcov_stage2}:}{Character. The standard-error method for the
    second-stage (undirected) parameters. One of \code{"h2"}, \code{"delta"}
    or \code{"none"}. Default is \code{"h2"}. The \code{"h2"} method requires
    \code{iv_samplestats = TRUE}.}
  \item{\code{iv_vcov_gamma_modelbased}:}{Logical. If \code{TRUE} (the
    default), the normal-theory weight matrix (gamma) used for the standard
    errors is based on the model-implied covariance matrix; if \code{FALSE},
    it is based on the unrestricted (H1) sample covariance matrix.}
  \item{\code{iv_vcov_jack_numerical}, \code{iv_vcov_jaca_numerical},
    \code{iv_vcov_jacb_numerical}:}{Logical. If \code{TRUE}, the corresponding
    Jacobians used in the standard-error computation are obtained numerically
    rather than from analytic expressions. The default is \code{FALSE}; these
    are mainly intended for checking the analytic derivatives.}
}
}

\subsection{Specifying instruments manually (the \code{|~} operator)}{

By default lavaan determines the model-implied instruments for each equation
automatically. The \code{|~} operator overrides this choice for a specific
equation: the left-hand side is the dependent variable of the equation (for a
latent variable, its scaling indicator), and the right-hand side lists the
instruments to use.

\preformatted{    model <- '
        ind60 =~ x1 + x2 + x3
        dem60 =~ y1 + y2 + y3 + y4
        dem60 ~ ind60
        # use only x2 and x3 as instruments for the dem60 ~ ind60 equation
        y1 |~ x2 + x3
    '
    fit <- sem(model, data = PoliticalDemocracy, estimator = "IV")}

This is useful when the automatically selected instruments are weak (see
\code{iv_weak} above) or when subject-matter knowledge suggests a different
set of instruments.
}

\subsection{Missing data}{

For continuous data with missing values, the IV estimator defaults to a
two-stage (FIML) approach: the saturated mean vector and covariance matrix are
estimated by the EM algorithm, the model is then estimated from these moments,
and the standard errors are corrected for the additional uncertainty stemming
from the EM estimates (Savalei & Bentler, 2009). This default is used when the
data contain missing values and the user does not set the \code{missing}
argument. Set \code{missing = "robust.two.stage"} for sandwich-corrected
standard errors (Savalei & Falk, 2014), or \code{missing = "listwise"} to
delete incomplete cases instead.
}

\subsection{Inspecting the instruments and the Sargan test}{

Several \code{\link{lavInspect}} fields are specific to IV estimation:
\describe{
  \item{\code{lavInspect(fit, "iv")}:}{(aliases \code{"miiv"},
    \code{"instruments"}) a per-equation table with the dependent variable
    (\code{lhs}), the regressors (\code{rhs}), the observed variables actually
    used (\code{lhs_new}, \code{rhs_new}), the equation type (\code{type}:
    \code{"miiv"}, \code{"ols"} or \code{"user"}), and the instruments
    (\code{iv}) for that equation.}
  \item{\code{lavInspect(fit, "sargan")}:}{a per-equation table with the Sargan
    overidentification test: the statistic (\code{sargan.stat}), its degrees of
    freedom (\code{df}) and p-value (\code{pvalue}), plus an adjusted p-value
    (\code{pvalue.adj}) when \code{iv_sargan_adjust} is set.}
  \item{\code{lavInspect(fit, "eqs")}:}{the underlying per-equation
    information.}
}
}
}
\references{
Bollen, K. A. (1996). An alternative two stage least squares (2SLS) estimator
for latent variable equations. \emph{Psychometrika}, 61(1), 109-121.

Bollen, K. A., & Maydeu-Olivares, A. (2007). A polychoric instrumental variable
(PIV) estimator for structural equation models with categorical variables.
\emph{Psychometrika}, 72(3), 309-326.

Sargan, J. D. (1958). The estimation of economic relationships using
instrumental variables. \emph{Econometrica}, 26(3), 393-415.

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing data:
Theory and application to auxiliary variables.
\emph{Structural Equation Modeling}, 16(3), 477-497.

Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms robust
full information maximum likelihood with incomplete nonnormal data.
\emph{Structural Equation Modeling}, 21(2), 280-302.

Staiger, D., & Stock, J. H. (1997). Instrumental variables regression with weak
instruments. \emph{Econometrica}, 65(3), 557-586.
}
\seealso{
\code{\link{lavaan}}, \code{\link{sem}}, \code{\link{lavOptions}},
\code{\link{lavInspect}}.
}
\examples{
## The classic Bollen (1989) Political Democracy example
model <- '
  # measurement model
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
'

## fit using MIIV-2SLS instead of maximum likelihood
fit <- sem(model, data = PoliticalDemocracy, estimator = "IV")
summary(fit)

## the model-implied instruments used for each equation
lavInspect(fit, "iv")

## the per-equation Sargan overidentification tests
lavInspect(fit, "sargan")

## pass estimator options via estimator = list(...)
fit2 <- sem(model, data = PoliticalDemocracy,
            estimator = list(estimator = "IV",
                             iv_varcov_method = "2RLS",
                             iv_sargan_adjust = "BH",
                             iv_weak = "warn"))
lavInspect(fit2, "sargan")

## specify instruments manually with the |~ operator:
## use only x2, x3 as instruments for the dem60 ~ ind60 equation
model.iv <- '
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem60 ~ ind60
    y1 |~ x2 + x3
'
fit3 <- sem(model.iv, data = PoliticalDemocracy, estimator = "IV")
lavInspect(fit3, "iv")

## equality constraints are honored (here: equal loadings across factors)
model.eq <- '
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + a*y2 + b*y3 + c*y4
    dem65 =~ y5 + a*y6 + b*y7 + c*y8
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
'
fit4 <- sem(model.eq, data = PoliticalDemocracy, estimator = "IV")
coef(fit4)
}
\keyword{models}
\keyword{regression}
