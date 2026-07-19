\name{estimator_js}
\alias{estimator_js}
\alias{js}
\alias{JS}
\alias{jsa}
\alias{JSA}
\alias{james.stein}
\alias{james.stein.aggregated}
\title{James-Stein (Conditional Expectation) Estimation}
\description{
Estimate a structural equation model by replacing the latent variables in
each equation by an estimate of their conditional expectation given
observed proxies, followed by least squares (Burghgraeve, De Neve &
Rosseel, 2021). Two variants are available: the scaling James-Stein
estimator (\code{estimator = "JS"}), which conditions on the
marker/scaling indicators, and the aggregated James-Stein estimator
(\code{estimator = "JSA"}), which conditions on reliability-weighted
aggregates of the indicators. This page documents how to invoke the
estimators, the options that can be passed via the
\code{estimator = list()} argument, and the current restrictions.
}
\details{
Like the instrumental-variable estimator (see
\code{\link{estimator_iv}}), the James-Stein estimators work
equation-by-equation: each measurement and structural equation is
estimated separately, the method is non-iterative, and no starting values
are needed. The two approaches differ in how the latent variables are
replaced by observable quantities. Where MIIV-2SLS projects on
instrumental variables, the James-Stein estimators regress on an estimate
of the conditional expectation of the latent variables given observed
proxies: \code{E(eta | y1) = R y1 + (I - R) E(eta)}, where \code{R} is
the estimated reliability matrix of the proxies (a James-Stein type
shrinkage estimator). For the scaling variant the proxies are the
marker/scaling indicators; for the aggregated variant they are weighted
sums of the indicators, with reliability-maximizing (alphamax; Bentler,
1968) weights. Regressing on these conditional expectations is
algebraically identical to correcting the least-squares normal equations
for the measurement-error variance of the proxies, so the estimates are a
closed-form function of the sample moments (no raw data are needed: the
model can also be fitted from \code{sample.cov} and \code{sample.nobs}).
The measurement-error variances of the indicators are estimated by a
Spearman-type communality estimator, applied per factor. Factors with
three or more indicators use the within-factor tetrads; for 2-indicator
factors the tetrads are formed with external variables (observed model
variables outside the factor), which requires the factor to be
correlated with other variables in the model. For a single-indicator
factor the reliability must come from a residual variance that is fixed
in the model (as the default lavaan rules do: the residual variance of a
single indicator is fixed to zero, so the latent variable coincides with
its indicator), or the residual variances can be supplied directly (see
\code{js_theta} below).

The estimator is selected by setting \code{estimator = "JS"} (scaling) or
\code{estimator = "JSA"} (aggregated), for example
\preformatted{    fit <- sem(model, data = Data, estimator = "JS")}

Estimation proceeds in two stages. In the first stage, the \emph{directed}
effects (factor loadings and regression coefficients) are estimated
equation-by-equation by least squares on the conditional expectations. In
the second stage, the \emph{undirected} effects (residual variances and
covariances) are estimated by a weighted least-squares step applied to the
residual moments, exactly as for the IV estimator (see
\code{js_varcov_method}). When \code{meanstructure = TRUE}, the
marker-intercepts-zero parameterization is used
(\code{marker.int.zero = TRUE}: the latent means are free and equal the
population means of their scaling indicators), and all free mean
parameters are re-estimated by a joint GLS solve (see
\code{js_mean_structure}).

Simple equality constraints (e.g. equal factor loadings) and general
linear equality constraints (e.g. \code{a == 2*b}) among the directed
coefficients are honored by a pooled (system) solve. Multiple groups are
supported, including cross-group equality constraints for measurement
invariance testing (e.g. \code{group.equal = c("loadings",
"intercepts")}). Without cross-group constraints, the per-group estimates
and standard errors are identical to those of separate single-group fits.

\subsection{Estimator options}{

Options specific to the JS/JSA estimators are passed as named elements of
a list to the \code{estimator} argument, alongside the estimator name:
\preformatted{    fit <- sem(model, data = Data,
               estimator = list(estimator = "JS",
                                js_varcov_method = "ULS",
                                js_gamma = "adf"))}

The option names use snake_case. For backward compatibility, dot.case
names (e.g. \code{js.varcov.method}) are also accepted and silently
converted.

\describe{
  \item{\code{js_small_sample}:}{Logical. If \code{TRUE} (the default),
    the Efron-Morris small-sample factor \code{(n-3)/(n-1)} is applied to
    the estimated measurement-error variances in the shrinkage step, as in
    Burghgraeve et al. (2021). If \code{FALSE}, the plain attenuation
    correction is used.}
  \item{\code{js_theta}:}{Character. How the measurement-error (residual)
    variances of the indicators are estimated. \code{"spearman"} (the
    default) uses the Spearman-type communality estimator per factor
    (within-factor tetrads for 3+ indicators; external tetrads for
    2-indicator factors; model-fixed residual variances for factors with
    fewer than 3 indicators); \code{"user"} takes the values from
    \code{js_theta_values}.}
  \item{\code{js_theta_values}:}{Named numeric vector. The residual
    variances of the indicators (on the same scale as the analyzed
    sample moments), used when \code{js_theta = "user"}. This makes it
    possible to plug in externally known reliabilities.}
  \item{\code{js_theta_bounds}:}{Character. Bounds applied to the
    Spearman communality estimates: \code{"wide"} (the default),
    \code{"standard"} or \code{"none"}. In addition, the resulting
    residual variances are always kept within \code{[0, var(y)]}.}
  \item{\code{js_varcov_method}:}{Character. The second-stage method used
    to estimate the residual variances and covariances. One of
    \code{"ULS"}, \code{"GLS"}, \code{"2RLS"}, \code{"RLS"} or
    \code{"NONE"}. Default is \code{"RLS"} (reweighted least squares). If
    \code{"NONE"}, the variance/covariance parameters are not estimated.}
  \item{\code{js_gamma}:}{Character. The asymptotic covariance matrix of
    the sample moments (Gamma) used in the standard errors: \code{"nt"}
    (normal-theory; the default) or \code{"adf"} (asymptotically
    distribution-free; requires complete raw data).}
  \item{\code{js_vcov_gamma_modelbased}:}{Logical. If \code{TRUE} (the
    default), the normal-theory Gamma is evaluated at the model-implied
    covariance matrix; if \code{FALSE} (or when the model-implied matrix
    is not positive definite), the unrestricted (h1) sample covariance
    matrix is used.}
  \item{\code{js_mean_structure}:}{Character. How the free mean
    parameters (observed intercepts and latent means) are estimated when
    \code{meanstructure = TRUE}: \code{"wls"} (the default) re-estimates
    them jointly by GLS, which also provides the means of the exogenous
    latent variables; \code{"moments"} keeps the per-equation intercepts.}
  \item{\code{js_jacobian}:}{Character. How the Jacobian of the
    estimation map (needed for the delta-method standard errors) is
    computed: \code{"analytic"} (the default) or \code{"numeric"}
    (numerical differentiation of the complete map; much slower for
    larger models). A few configurations are not covered by the analytic
    Jacobian (general linear equality constraints among the
    variance/covariance parameters, a non-positive-definite or
    non-converged second-stage weight matrix); these fall back to the
    numerical Jacobian automatically.}
}
}

\subsection{Standard errors and test statistics}{

By default (\code{se = "standard"}), delta-method standard errors are
computed by differentiating the complete estimation map with respect to
the sample moments, and combining the resulting Jacobian with the
asymptotic covariance matrix of the moments (\code{js_gamma}). The
Jacobian is computed analytically in a single chain-rule sweep (see
\code{js_jacobian}), so the standard errors are non-iterative as well.
Because the Jacobian runs through all estimation steps, the additional
variability caused by estimating the reliability matrix (Theorem 2 in
Burghgraeve et al., 2021) -- and, for JSA, the aggregation weights -- is
accounted for automatically. The full covariance matrix of the parameter
estimates is available, so defined parameters (\code{:=}) obtain proper
delta-method standard errors. Alternatively, \code{se = "bootstrap"}
computes bootstrap standard errors, and \code{se = "none"} skips the
standard errors.

There is no discrepancy-function chi-square statistic; the default test
statistic is Browne's (1984) residual-based test
(\code{test = "browne.residual.nt"}; the ADF version can be requested
instead), and the fit measures reported by \code{\link{fitMeasures}()}
are based on it. Likelihood-based quantities (logl, AIC, BIC) are not
available. Modification indices are not available (as for the IV
estimator).
}

\subsection{Missing data}{

By default, cases with missing values are removed listwise. When the data
contain missing values, \code{missing = "two.stage"} is selected
automatically (it can also be requested explicitly, as can
\code{missing = "robust.two.stage"}): the saturated (EM) moments are used
for point estimation, and the standard errors are based on the two-stage
moment covariance (Savalei & Bentler, 2009; Savalei & Falk, 2014). The
test statistic then uses the model-based variant
(\code{"browne.residual.nt.model"}).
}

\subsection{Residual covariances}{

Residual covariances among the indicators (including covariances that
involve a scaling indicator) are supported: each equation conditions on a
\emph{clean} proxy set, i.e. observed variables whose residuals are
uncorrelated with the residual of the equation's dependent variable (and
with each other). When a default proxy (the scaling indicator, or an
aggregate member) is invalidated by a residual covariance, it is replaced
by (an aggregate of) clean pure indicators of the same factor, scaled by
their estimated loadings -- the James-Stein analog of dropping
contaminated instruments in MIIV estimation. Contaminated correlations
are likewise excluded from the Spearman reliability estimates. The
residual covariances themselves are estimated in the second stage. This
requires enough clean indicators: every needed reliability must have at
least one uncontaminated tetrad, and every equation needs at least one
clean proxy per latent regressor; an informative error is raised
otherwise.
}

\subsection{Restrictions}{

The JS/JSA estimators currently require: continuous data (no ordered
categorical variables), single-level models, a marker/scaling indicator
for every latent variable (\code{std.lv = TRUE} is not supported), a
recursive structural part (no feedback loops), and no higher-order
factors. The disturbance of a dependent variable may not be correlated
with (an ancestor of) one of its regressors. Exogenous observed
covariates are supported (they are part of the conditioning set of each
equation; note that \code{fixed.x} is set to \code{FALSE}).
}
}
\references{
Burghgraeve, E., De Neve, J., & Rosseel, Y. (2021). Estimating structural
equation models using James-Stein type shrinkage estimators.
\emph{Psychometrika}, 86(1), 96-130.

Bentler, P. M. (1968). Alpha-maximized factor analysis (alphamax): Its
relation to alpha and canonical factor analysis. \emph{Psychometrika},
33(3), 335-345.

Bollen, K. A. (1996). An alternative two stage least squares (2SLS)
estimator for latent variable equations. \emph{Psychometrika}, 61(1),
109-121.

Browne, M. W. (1984). Asymptotically distribution-free methods for the
analysis of covariance structures. \emph{British Journal of Mathematical
and Statistical Psychology}, 37(1), 62-83.

Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing
data: Theory and application to auxiliary variables.
\emph{Structural Equation Modeling}, 16(3), 477-497.

Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
robust full information maximum likelihood with incomplete nonnormal
data. \emph{Structural Equation Modeling}, 21(2), 280-302.
}
\seealso{
\code{\link{lavaan}}, \code{\link{sem}}, \code{\link{estimator_iv}},
\code{\link{lavOptions}}, \code{\link{lavInspect}}.
}
\examples{
## a CFA example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
fit <- cfa(HS.model, data = HolzingerSwineford1939, estimator = "JS")
summary(fit, fit.measures = TRUE)

## the aggregated variant
fit.a <- cfa(HS.model, data = HolzingerSwineford1939, estimator = "JSA")
coef(fit.a)

## a full structural equation model
model <- '
  # measurement model
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
'
fit2 <- sem(model, data = PoliticalDemocracy, estimator = "JS")
summary(fit2)

## pass estimator options via estimator = list(...)
fit3 <- sem(model, data = PoliticalDemocracy,
            estimator = list(estimator = "JS",
                             js_varcov_method = "ULS",
                             js_gamma = "adf"))

\dontrun{
## multiple groups: metric measurement invariance (equal loadings)
fit.metric <- cfa(HS.model, data = HolzingerSwineford1939,
                  estimator = "JS",
                  group = "school", group.equal = "loadings")
summary(fit.metric)
}
}
\keyword{models}
\keyword{regression}
