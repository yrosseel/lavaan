\name{sam}
\alias{sam}
\alias{fsr}
\title{Fit Structural Equation Models using the SAM approach}
\description{
Fit a Structural Equation Model (SEM) using the Structural After Measurement
(SAM) approach.}
\usage{
sam(model = NULL, data = NULL, aux = NULL, cmd = "sem", se = "twostep",
    mm_list = NULL, mm_args = list(bounds = "wide.zerovar"), 
    struc_args = list(estimator = "ML"), 
    sam_method = "local", ..., 
    local_options = list(M.method = "ML", lambda.correction = TRUE, 
                         alpha.correction = 0L, twolevel.method = "h1"), 
    global_options = list(), 
    bootstrap = list(R = 1000L, type = "ordinary", show.progress = FALSE),
    output = "lavaan",
    bootstrap_args = bootstrap
)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See 
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (e.g., the output of the \code{lavParTable()} function) is also
    accepted.}
\item{data}{A data frame containing the observed variables used in
    the model.}
\item{aux}{Character vector. Names of auxiliary observed variables, forwarded
    to the underlying measurement and structural model fits. See the
    \code{\link{lavaan}} function for more details.}
\item{cmd}{Character. Which command is used to run the sem models. The possible
    choices are \code{"sem"}, \code{"cfa"} or \code{"lavaan"}, determining how
    we deal with default options.}
\item{se}{Character. The type of standard errors that are used in the
    final (structural) model. If \code{"twostep"} (the default), the
    standard errors take the estimation uncertainty of the first (measurement)
    stage into account. If \code{"standard"}, this uncertainty is ignored,
    and we treat the measurement information as known. If \code{"none"},
    no standard errors are computed. If the data are clustered (a
    \code{cluster=} argument is provided, without multilevel model syntax),
    \code{se = "twostep"} is automatically replaced by
    \code{se = "twostep.robust"}, and the standard errors (also for
    \code{se = "local"} and \code{se = "bootstrap"}, which resamples whole
    clusters) account for the clustering. For two-level models,
    \code{se = "local"} is available for the single-group continuous
    complete-data setting; it uses the cluster-sandwich covariance of the
    saturated (h1) estimates, and also enables the corrected structural
    test statistic.}
\item{mm_list}{List. Define the measurement blocks. 
     Each element of the list should
     be either a single name of a latent variable, or a vector of latent
     variable names.  If omitted, a separate measurement
     block is used for each latent variable.}
\item{mm_args}{List. Optional arguments for the fitting
    function(s) of the measurement block(s) only. See \code{\link{lavOptions}} 
    for a complete list.}
\item{struc_args}{List. Optional arguments for the fitting function of the
    structural part only. See \code{\link{lavOptions}} for a complete list.}
\item{sam_method}{Character. Can be set to \code{"local"}, \code{"global"}
    or \code{"fsr"}. In the latter case, the results are the same as if
    Bartlett factor scores were used, without any bias correction.}
\item{...}{Many more options can be specified, using 'name = value'.
    See \code{\link{lavOptions}} for a complete list. These options affect
    both the measurement blocks and the structural part.}
\item{local_options}{List. Options specific for local SAM method (these
     options may change over time). If \code{lambda.correction = TRUE},
     we ensure that the variance matrix of the latent variables (VETA)
     is positive definite. The \code{alpha.correction} option must be an
     integer. Acceptable values are in the range 0 to
     N-1. If zero (the default), no small sample correction is performed,
     and the bias-correction is the same as with local SAM. When equal to N-1,
     the bias-correction is eliminated, and the results are the same as
     naive FSR. Typical values are 0, P+1 (where P is the number of
     predictors in the structural model), P+5, and (N-1)/2.}
\item{global_options}{List. Options specific for global SAM method (not
     used for now).}
\item{bootstrap}{List. Only used when \code{se = "bootstrap"}. Typical
     elements of this list are \code{R}: the number of bootstrap samples, and
     \code{type}, which can be set to \code{"ordinary"} (the default) or
     \code{"parametric"}.}
\item{output}{Character. If \code{"lavaan"}, a lavaan object is returned. If
    \code{"list"}, a list is returned with all the ingredients from the 
    different stages.}
\item{bootstrap_args}{List. Deprecated, use \code{bootstrap} instead.}
}
\details{
    The \code{sam} function automates the SAM approach by first
    estimating the measurement part of the model and then the structural
    part of the model. See the reference for more details.

    Note that in the current implementation, all indicators of latent
    variables must be observed. As a result, second-order factor structures
    are not yet supported.

    If \code{sam.method} is \code{"local"}, \code{"fsr"} or \code{"cfsr"},
    the structural part of the model is fitted to the estimated latent
    (co)variances (and means) from step 1. This step-2 structural model is
    stored in the \code{@internal} slot of the returned object, and can be
    retrieved by \code{lavInspect(fit, "sam.struc.fit.object")}. The
    functions \code{\link{fitMeasures}} (and
    \code{summary(fit, fit.measures = TRUE)}), \code{\link{lavResiduals}}
    (and \code{summary(fit, residuals = TRUE)}), \code{\link{modindices}}
    and \code{\link{lavTestLRT}} (when all compared models are local SAM
    models) operate on this structural model: the reported quantities
    describe the fit of the structural part only, conditional on the (fixed)
    measurement model of step 1. When available, the reported (scaled) test
    statistic and the robust fit measures are corrected for the estimation
    of the measurement model in step 1.

    If \code{sam.method} is \code{"global"}, the returned object is an
    ordinary (joint) lavaan object, and \code{fitMeasures} returns the usual
    fit measures of the joint model. Note, however, that the parameters were
    estimated using a two-step procedure; the behavior of fit measures in
    this setting has not been studied, and they should be interpreted with
    caution.
}
\section{Exogenous covariates and conditional.x}{
    Observed exogenous covariates may enter the structural part of the
    model. With continuous data, they are simply treated as (pseudo) latent
    variables (the default \code{conditional.x = FALSE}), but
    \code{conditional.x = TRUE} is supported as well (single-group,
    single-level models): the structural part is then fitted to the latent
    intercepts, the latent-on-covariate slopes, and the residual latent
    (co)variances, all conditional on the covariates.

    With categorical (ordered) data, \code{conditional.x = TRUE} is the
    lavaan default (when exogenous covariates are present), and this is also
    the recommended setting: the covariates are conditioned on (as in a
    probit regression), rather than being absorbed into the joint polychoric
    correlation matrix -- which would (wrongly) treat them as jointly normal,
    a real concern for binary/dummy covariates. In this setting, the
    measurement blocks of step 1 are estimated with
    \code{conditional.x = TRUE} as well (each block latent variable is
    regressed on all exogenous covariates, MIMIC style), so that the
    measurement parameters live on the same conditional latent-response
    scale (\code{Var(y*|x) = 1}, delta parameterization) as the joint model.
    An \emph{ordered endogenous} variable in the structural part is
    supported too (with or without covariates): it enters the structural
    model on its latent-response scale, while its thresholds are fixed at
    their sample-based values.

    Not (yet) available in combination with \code{conditional.x = TRUE}:
    multiple groups (for the local standard errors and the corrected test
    statistic), multilevel models, latent interactions, and the Yuan-Chan
    scaled global test statistic of \code{sam_method = "global"} (the
    standard, unscaled test is reported instead, with a warning).
}
\value{
  If \code{output = "lavaan"}, an object of class
  \code{\linkS4class{lavaan}}, for which several methods
  are available, including a \code{summary} method. If \code{output = "list"},
  a list.
}
\references{
Rosseel and Loh (2021). A structural-after-measurement approach to Structural
Equation Modeling. Psychological Methods. Advance online publication.
https://dx.doi.org/10.1037/met0000503
}
\seealso{
\code{\link{lavaan}}
}
\examples{
## The industrialization and Political Democracy Example 
## Bollen (1989), page 332
model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit.sam <- sam(model, data = PoliticalDemocracy,
               mm.list = list(ind = "ind60", dem = c("dem60", "dem65")))
summary(fit.sam)
}
