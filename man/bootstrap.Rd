\name{lavBootstrap}
\alias{bootstrapLavaan}
\alias{lavBootstrap}
\title{Bootstrapping a Lavaan Model}
\description{Bootstrap the LRT, or any other statistic (or vector
of statistics) you can extract from a fitted lavaan object.}
\usage{
lavBootstrap(object, r = 1000L, type = "ordinary", verbose = FALSE, 
                fun = "coef", keep_idx = FALSE,
                parallel = c("no", "multicore", "snow"),
                ncpus = min(2L, max(1L, parallel::detectCores() - 2L, na.rm = TRUE)),
                cl = NULL, iseed = NULL, h0_rmsea = NULL, method = "refit",
                guard = TRUE, guard_tol = 0.25, ...)
bootstrapLavaan(object, r = 1000L, type = "ordinary", verbose = FALSE, 
                fun = "coef", keep_idx = FALSE,
                parallel = c("no", "multicore", "snow"),
                ncpus = min(2L, max(1L, parallel::detectCores() - 2L, na.rm = TRUE)),
                cl = NULL, iseed = NULL, h0_rmsea = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{r}{Integer. The number of bootstrap draws.}
\item{type}{If \code{"ordinary"} or \code{"nonparametric"}, the usual (naive)
bootstrap method is used. If \code{"bollen.stine"}, the data is first
transformed such that the null hypothesis holds exactly in the resampling
space. If \code{"yuan"}, the data is first transformed by combining data and
theory (model), such that the resampling space is closer to the population
space. Note that both \code{"bollen.stine"} and \code{"yuan"} require the data
to be continuous. They will not work with ordinal data.  If
\code{"parametric"}, the parametric bootstrap approach is used; currently, this
is only valid for continuous data following a multivariate normal distribution.
See references for more details. If the fitted model uses clustered or
multilevel data (that is, a \code{cluster=} argument was used, either for a
two-level model or to obtain cluster-robust standard errors), only the
\code{"ordinary"} (nonparametric) bootstrap is available, and whole clusters
(rather than individual observations) are resampled; see Details.}
\item{fun}{A function which when applied to the \code{\linkS4class{lavaan}} 
object returns a vector containing the statistic(s) of interest.
The default is \code{fun="coef"}, returning the estimated values of the
free parameters in the model.}
\item{\dots}{Other named arguments for \code{fun} which are passed
    unchanged each time it is called.}
\item{verbose}{If \code{TRUE}, show information for each bootstrap draw.}
\item{keep_idx}{If \code{TRUE}, store the indices of each bootstrap run (i.e.,
the observations that were used for this bootstrap run) as an attribute.}
\item{parallel}{The type of parallel operation to be used (if any).  If
missing, the default is \code{"no"}.  }
\item{ncpus}{Integer: number of processes to be used in parallel operation.
By default, this is 2 (following CRAN policy). You may set this to a higher
value (for example, the number of cores as detected by
\code{parallel::detectCores()} minus one) if you wish to use more processes.}
\item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
\code{parallel = "snow"}.  If not supplied, a cluster on the local machine is
created for the duration of the \code{lavBootstrap} call.}
\item{iseed}{An integer to set the seed. Or NULL if no reproducible results are
  needed. This works for both serial (non-parallel) and parallel settings.
  Internally, \code{RNGkind()} is set to \code{"L'Ecuyer-CMRG"} if
  \code{parallel = "multicore"}. If \code{parallel = "snow"} (under windows),
  \code{parallel::clusterSetRNGStream()} is called which
  automatically switches to \code{"L'Ecuyer-CMRG"}. When \code{iseed} is not 
  NULL, \code{.Random.seed} (if it exists) in the global environment is 
  left untouched.}
\item{h0_rmsea}{Only used if \code{type="yuan"}. Allows one to do the Yuan
bootstrap under the hypothesis that the population RMSEA equals a specified
value.}
\item{method}{Character. If \code{"refit"} (the default), the model is refitted
for every bootstrap sample. If \code{"ij"}, the bootstrap draws are
approximated by the (first-order) infinitesimal jackknife, without refitting
the model. If \code{"hoij2"}, the second-order infinitesimal jackknife is
used instead; see Details. The approximations are only available for the
ordinary bootstrap and \code{fun = "coef"}.}
\item{guard}{Logical. Only used if \code{method = "hoij2"}. If \code{TRUE}
(the default), the model is refitted for those bootstrap samples where the
approximation is not accurate enough; see Details.}
\item{guard_tol}{Numeric. Only used if \code{method = "hoij2"} and
\code{guard = TRUE}. The model is refitted for a bootstrap sample if the
estimated distance between the approximate draw and the exact solution
(expressed in standard error units) exceeds this value. The default is
0.25.}
}
\author{Yves Rosseel. Ed Merkle contributed
Yuan's bootstrap. Improvements to Yuan's bootstrap were contributed
by Hao Wu and Chuchu Cheng. The handling of iseed was contributed by
Shu Fai Cheung.}
\value{For \code{lavBootstrap()}, the bootstrap distribution of the value(s)
returned by \code{fun}, when the object can be simplified to a vector.
}
\details{The fun function can return either a scalar or a numeric vector.
This function can be an existing function (for example \code{coef}) or
can be a custom defined function. For example:
\preformatted{
myFUN <- function(x) {
    # require(lavaan)
    modelImpliedCov <- fitted(x)$cov
    vech(modelImpliedCov)
}
}
If \code{parallel="snow"}, it is imperative that the \code{require(lavaan)}
is included in the custom function.

When the fitted model was estimated with clustered or multilevel data (i.e., a
\code{cluster=} argument was used), the bootstrap automatically switches to a
\emph{cluster bootstrap}: instead of resampling individual observations, whole
(level-2) clusters are resampled with replacement, keeping all level-1
observations within a selected cluster together. The number of clusters is held
fixed (equal to the number of clusters in the original data), while the total
number of observations may vary from one bootstrap sample to the next. A cluster
that is drawn more than once is treated as several distinct clusters. This
resampling scheme preserves the within-cluster dependence and is the
appropriate bootstrap both for two-level models and for single-level models with
cluster-robust standard errors. It is currently only available for a single
cluster (level-2) variable, and the \code{"bca"} confidence interval type
(in \code{\link{parameterEstimates}}) is not supported in this case.

When the fitted model contains rotated exploratory factors (EFA or ESEM; that
is, the model syntax uses the \code{efa()} modifier together with a
\code{rotation=} method other than \code{"none"}), bootstrap standard errors and
confidence intervals (requested via \code{se = "bootstrap"} when fitting the
model) are computed for the \emph{rotated} solution. Each bootstrap sample is
refit and rotated, and its factors are then aligned (reordered and reflected) to
the original (full-sample) rotated solution before their variability is used.
This alignment resolves the rotational indeterminacy (the sign and order of the
exploratory factors), which would otherwise inflate the standard errors. Note
that the raw bootstrap coefficients returned by \code{lavBootstrap(object,
fun = "coef")} are \emph{not} aligned to a reference solution in this way.

If \code{method = "ij"}, the model is not refitted. A bootstrap sample is
equivalent to giving each observation a (multinomial) weight, and the
estimates are a smooth function of these weights. The infinitesimal jackknife
(IJ) approximates this function by a first-order Taylor expansion around the
original (unit) weights, using the casewise scores and the observed
information matrix of the original fit. The same bootstrap samples are used
as for \code{method = "refit"} (given the same \code{iseed}, and no
parallelization). The covariance matrix of the IJ draws equals (up to Monte
Carlo error) the Huber-White (sandwich) covariance matrix; the draws are
therefore mostly useful to obtain percentile-type confidence intervals for
(nonlinear) functions of the parameters. Because the approximation is linear
in the weights, it does not capture any skewness that is due to the
estimator itself (for example, for variance parameters in small samples).
The IJ method is currently only available for single-level models that are
estimated by ML using continuous (complete or incomplete) data.

If \code{method = "hoij2"}, a second-order term is added to the expansion
(the second-order infinitesimal jackknife, or HOIJ-2). This term accounts for
the change in the curvature of the objective function, both due to the
reweighting of the observations and due to the change in the parameter
values. As a result, the draws can also reflect the skewness of the
estimator itself. The second-order term needs four evaluations of the
gradient for every bootstrap sample, which is still much cheaper than
refitting the model. Note that HOIJ-2 is a local (quadratic) approximation:
it is usually closer to the refitted estimates than IJ for typical bootstrap
samples, but it can be less accurate for the most extreme bootstrap samples,
and these samples determine the limits of percentile-type confidence
intervals. For this reason, a guard is used by default
(\code{guard = TRUE}). For every bootstrap sample, the gradient of the
objective function is evaluated at the approximate draw, using the sample
statistics of that bootstrap sample. If the draw were exact, this gradient
would be zero. The size of the gradient is converted to an estimate of the
distance between the draw and the exact solution, expressed in standard
error units (the root mean square over the parameters). If this distance
exceeds \code{guard_tol}, the model is refitted for this bootstrap sample,
and the approximate draw is replaced by the refitted one. In large samples
and for well-behaved models, (almost) no bootstrap samples need a refit. In
small samples, or if the model is only weakly identified, many bootstrap
samples may need a refit, and the gain in speed is smaller; a note is given
if this happens for more than half of the bootstrap samples. The distances
and the indices of the refitted draws are available as the attributes
\code{"dist"} and \code{"guard"} of the returned matrix. The guard is not
perfect: the distance is itself an estimate, and a few inaccurate draws may
not be detected. The HOIJ-2 method is currently only
available for single-level models that are estimated by ML using complete
continuous data, possibly with multiple groups and (linear or nonlinear)
equality constraints (but without inequality constraints or bounds).

Both methods can be requested when fitting the model, by using
\code{se = "bootstrap"} together with, for example,
\code{bootstrap = list(R = 1000, method = "ij")} or
\code{bootstrap = list(R = 1000, method = "hoij2", guard_tol = 0.1)}.}
\references{
Bollen, K. and Stine, r. (1992) Bootstrapping Goodness of Fit Measures in
        Structural Equation Models. Sociological Methods and Research, 21,
        205--229.

Yuan, K.-H., Hayashi, K., & Yanagihara, H. (2007). A class of population
covariance matrices in the bootstrap approach to covariance structure analysis.
Multivariate Behavioral Research, 42, 261--281.

Field, C. A., & Welsh, A. H. (2007). Bootstrapping clustered data. Journal of
the Royal Statistical Society: Series B, 69, 369--390.
}

\examples{
# fit the Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939, se="none")

# get the test statistic for the original sample
T.orig <- fitMeasures(fit, "chisq")

# bootstrap to get bootstrap test statistics
# we only generate 10 bootstrap samples in this example; in practice
# you may wish to use a much higher number
T.boot <- lavBootstrap(fit, r=10, type="bollen.stine",
                          fun = fitMeasures, fit_measures = "chisq")

# compute a bootstrap based p-value
pvalue.boot <- length(which(T.boot > T.orig))/length(T.boot)
}
