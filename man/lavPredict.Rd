\name{lavPredict}
\alias{lavPredict}
\alias{lavpredict}
\title{Predict the values of latent variables (and their indicators)}
\description{
The main purpose of the \code{lavPredict()} function is to compute (or
`predict') individual scores for the latent variables in the model
(`factor scores'). NOTE: the goal of this
function is NOT to predict future values of dependent variables as in the
regression framework! (For models with only continuous observed variables, the function \code{lavPredictY()} supports this.)}
\usage{
lavPredict(object, newdata = NULL, type = "lv", method = "EBM",
           transform = FALSE, se = "none", acov = "none",
           label = TRUE, fsm = FALSE, mdist = FALSE, rel = FALSE,
           append_data = FALSE, assemble = FALSE,
           level = 1L, optim_method = "bfgs", eta = NULL,
           parallel = c("auto", "no", "multicore", "snow"),
           ncpus = NULL, cl = NULL,
           drop_list_single_group = TRUE,
           mdist_draws = 2000L,
           ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{newdata}{An optional data.frame, containing the same variables as
the data.frame used when fitting the model in object. For multilevel
models, the cluster variable(s) must be included; the cluster structure
and (for \code{missing = "ml"}) the missing-data patterns are rebuilt
from newdata, so the predictions are conditional on the newdata only
(the cluster labels need not match those of the original data).}
\item{type}{A character string. If \code{"lv"}, estimated values for the latent
variables in the model are computed. If \code{"ov"}, model predicted values for
the indicators of the latent variables in the model are computed. If
\code{"yhat"}, the estimated values for the observed indicators are computed,
given the user-specified values for the latent variables provided by the
\code{eta} argument.  If \code{"fy"}, densities (or probabilities) are computed
for each observed indicator, given the user-specified values for the latent
variables provided by the \code{eta} argument.}
\item{method}{A character string. In the linear case (when the indicators are
continuous), the possible options are \code{"regression"} or \code{"Bartlett"}.
In the categorical case, the two options are \code{"EBM"} for
the Empirical Bayes Modal approach, and \code{"ML"} for the maximum
likelihood approach. For higher-order factors (which have no observed
indicators of their own), the \code{"Bartlett"} scores are computed from the
collapsed measurement model (regressing the observed indicators of the
lower-order factors on the higher-order factors), as in \code{\link{sam}}.}
\item{transform}{Logical. If \code{TRUE}, transform the factor scores (per
group) so that their mean and variance-covariance matrix matches the
model-implied mean and variance-covariance matrix. This may be useful if the
individual factor scores will be used in a follow-up (regression) analysis.
Note: the standard errors (if requested) are not transformed (yet). The resulting
factor scores are often called correlation-preserving factor scores.}
\item{se}{Character. If \code{"none"}, no standard errors are computed.
If \code{"standard"}, standard errors are computed (assuming the
parameters of the measurement model are known). The standard errors are
returned as an attribute. For continuous data, the (naive) standard errors are
the same for every observation. For categorical data, where the factor scores
are obtained by numerical optimization (\code{method = "EBM"} or
\code{method = "ML"}), the standard errors are based on the curvature of the
objective at the optimum (the inverse of its Hessian, a Laplace/observed
-information approximation) and therefore differ from one response pattern to
the next; the result is an (nobs x nfactor) matrix per group.}
\item{acov}{Similar to the \code{"se"} argument, but optionally returns the full
sampling covariance matrix of the factor scores as an attribute. For
continuous data this is a single matrix per group; for categorical data it is
a list with one (nfactor x nfactor) matrix per observation.}
\item{label}{Logical. If TRUE, the columns in the output are labeled.}
\item{fsm}{Logical. If TRUE, return the factor score matrix as an attribute.
           Only for numeric data.}
\item{mdist}{Logical. If TRUE, the (squared)
Mahalanobis distances of the factor scores (if \code{type = "lv"}) or
the casewise residuals (if \code{type = "resid"}) are returned as an
attribute. If (some of) the data is ordered categorical, the distances
are generalized as in Mansolf and Reise (2017): the continuous distance is
replaced by its expected value over the region of the (multivariate
normally distributed) latent response vector that is consistent with the
observed response pattern; this expectation is approximated by Monte Carlo
integration (using \code{mdist_draws} random draws per distinct response
pattern; the result therefore depends on the state of the random number
generator, and \code{set.seed()} can be used for reproducibility).
Observed continuous variables are conditioned on, and missing values are
integrated out. In this categorical setting, \code{type = "resid"} is
allowed, and returns the expected casewise residuals of the latent
responses (a by-product of the same Monte Carlo integration).}
\item{rel}{Logical. Only used if \code{type = "lv"}. If TRUE,
the factor reliabilities are returned as an attribute. (The squared
values are often called the factor determinacies.)}
\item{append_data}{Logical. Only used when \code{type = "lv"}. If TRUE,
the original data (or the data provided
in the newdata argument) is appended to the factor scores.}
\item{assemble}{Logical. If TRUE,
the separate groups are reassembled into a single data.frame
with a group column, having the same dimensions as the original (or
newdata) dataset.}
\item{level}{Integer. Only used in a multilevel SEM.
If \code{level = 1}, only factor scores for latent variable
defined at the first (within) level are computed; if \code{level = 2},
only factor scores for latent variables defined at the second (between) level
are computed. For two-level models with \code{missing = "ml"}, the
posterior means of the cluster-level components and of the missing
values (given all the observed data; computed as in the E-step of the
EM algorithm) are used, so that the factor scores are the exact
posterior means given all the observed data (new in 0.7-1).}
\item{optim_method}{Character string. Only used in the categorical case.
If \code{"nlminb"} (the default in 0.5), the \code{"nlminb()"} function is used
for the optimization. If \code{"bfgs"} or \code{"BFGS"} (the default in 0.6),
the \code{"optim()"} function is used with the BFGS method.}
\item{eta}{An optional matrix or list, containing latent variable values
  for each observation. Used for computations when \code{type = "ov"}.}
\item{parallel}{Character. Only used in the categorical case, where factor
  scores are computed by a per-observation numerical optimization. The options
  are \code{"no"} (serial), \code{"multicore"} (fork-based, not available on
  Windows) and \code{"snow"} (a PSOCK cluster). The default, \code{"auto"},
  uses \code{"multicore"} automatically (on non-Windows platforms) when a
  large number of distinct optimizations is required, and runs serially
  otherwise. Because the computation is deterministic, the results do not
  depend on the value of this argument.}
\item{ncpus}{Integer. The number of processes to use in parallel computation.
  The default is to use all available cores, minus two.}
\item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use when
  \code{parallel = "snow"}. If not supplied, a cluster on the local machine is
  created for the duration of the call.}
\item{drop_list_single_group}{Logical. If \code{FALSE}, the results are
    returned as
    a list, where each element corresponds to a group (even if there is only
    a single group). If \code{TRUE}, the list will be unlisted if there is
    only a single group.}
\item{mdist_draws}{Integer. Only used if \code{mdist = TRUE} and (some of)
    the data is ordered categorical: the number of Monte Carlo draws per
    distinct response pattern used to approximate the expected Mahalanobis
    distances.}
\item{...}{To support old argument names.}
}
\details{
The \code{predict()} function calls the \code{lavPredict()} function
with its default options.

If there are no latent variables in the model, \code{type = "ov"} will
simply return the values of the observed variables. Note that this function
can not be used to `predict' values of dependent variables, given the
values of independent variables (in the regression sense). In other words,
the structural component is completely ignored (for now).
}
\references{
For an overview (and evaluation) of the various factor score methods, see:

Grice, J. W. (2001). Computing and evaluating factor scores.
\emph{Psychological Methods, 6}(4), 430-450. \doi{10.1037/1082-989X.6.4.430}

For the (continuous) regression and Bartlett methods, see:

Bartlett, M. S. (1937). The statistical conception of mental factors.
\emph{British Journal of Psychology, 28}, 97-104.
\doi{10.1111/j.2044-8295.1937.tb00863.x}

Bentler, P. M., & Yuan, K.-H. (1997). Optimal conditionally unbiased
equivariant factor score estimators. In M. Berkane (Ed.), \emph{Latent
variable modeling and applications to causality} (pp. 259-281). New York:
Springer-Verlag. \doi{10.1007/978-1-4612-1842-5_14}

For the (categorical) Empirical Bayes Modal (EBM) and Maximum Likelihood (ML)
methods, see:

Skrondal, A., & Rabe-Hesketh, S. (2004). \emph{Generalized latent variable
modeling: Multilevel, longitudinal, and structural equation models.} Boca
Raton, FL: Chapman & Hall/CRC.

For the Mahalanobis distances with ordered categorical data
(\code{mdist = TRUE}), see:

Mansolf, M., & Reise, S. P. (2017). Case diagnostics for factor analysis of
ordered categorical data with applications to person-fit measurement.
\emph{Structural Equation Modeling: A Multidisciplinary Journal, 25}(1),
86-100. \doi{10.1080/10705511.2017.1367926}

For the correlation-preserving factor scores (\code{transform = TRUE}), see:

ten Berge, J. M. F., Krijnen, W. P., Wansbeek, T., & Shapiro, A. (1999). Some
new results on correlation-preserving factor scores prediction methods.
\emph{Linear Algebra and its Applications, 289}(1-3), 311-318.
\doi{10.1016/S0024-3795(97)10007-6}
}
\seealso{
\code{\link{lavPredictY}} to predict y-variables given x-variables.
}
\examples{
data(HolzingerSwineford1939)

## fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
head(lavPredict(fit))
head(lavPredict(fit, type = "ov"))


## ---------------------------------------------
## standard errors for the factor scores (se =)
## ---------------------------------------------

## the standard errors are returned as the "se" attribute (a list, one
## (nobs x nfactor) matrix per group)

## for continuous indicators, the (naive) standard errors are the same
## for every observation
fscores <- lavPredict(fit, se = "standard")
attr(fscores, "se")[[1]]

## for categorical indicators, the factor scores are obtained by numerical
## optimization, and their standard errors differ from one response pattern
## (observation) to the next

## (a two-factor model is used here, as the numerical optimization and the
##  Monte Carlo integration below both get expensive as the number of
##  ordered indicators grows)
HS.model.ord <- ' visual  =~ x1 + x2 + x3
                  textual =~ x4 + x5 + x6 '
HS.ord <- HolzingerSwineford1939
HS.ord[ , paste0("x", 1:6)] <-
    lapply(HS.ord[ , paste0("x", 1:6)], function(x) ordered(cut(x, 3)))
fit.ord <- cfa(HS.model.ord, data = HS.ord, ordered = paste0("x", 1:6))
fscores <- lavPredict(fit.ord, se = "standard")
head(attr(fscores, "se")[[1]])


## -------------------------------------------------
## casewise Mahalanobis distances (case diagnostics)
## -------------------------------------------------

## continuous data: the squared distances are exact
fscores <- lavPredict(fit, mdist = TRUE)
head(attr(fscores, "mdist")[[1]])

## ordered categorical data: expected (squared) distances of the latent
## responses, obtained by Monte Carlo integration (Mansolf & Reise, 2017);
## set the seed for reproducibility
set.seed(123)
resid.ord <- lavPredict(fit.ord, type = "resid", mdist = TRUE)
## expected casewise residuals of the latent responses:
head(resid.ord)
## expected squared residual-based distances (person fit):
head(attr(resid.ord, "mdist")[[1]])


## ------------------------------------------
## merge factor scores to original data.frame
## ------------------------------------------

idx <- lavInspect(fit, "case.idx")
fscores <- lavPredict(fit)
## loop over factors
for (fs in colnames(fscores)) {
  HolzingerSwineford1939[idx, fs] <- fscores[ , fs]
}
head(HolzingerSwineford1939)


## multigroup models return a list of factor scores (one per group)
data(HolzingerSwineford1939)
mgfit <- update(fit, group = "school", group.equal = c("loadings","intercepts"))

idx <- lavInspect(mgfit, "case.idx") # list: 1 vector per group
fscores <- lavPredict(mgfit)         # list: 1 matrix per group
## loop over groups and factors
for (g in seq_along(fscores)) {
  for (fs in colnames(fscores[[g]])) {
    HolzingerSwineford1939[ idx[[g]], fs] <- fscores[[g]][ , fs]
  }
}
head(HolzingerSwineford1939)

## -------------------------------------
## Use factor scores in subsequent models
## -------------------------------------

## see Examples in semTools package: ?plausibleValues
}
