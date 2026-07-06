\name{lavSimulateData}
\alias{lavSimulateData}
\alias{lav_data_simulate_old}
\alias{simulateData}
\title{Simulate Data From a Lavaan Model Syntax}
\description{Simulate data starting from a lavaan model syntax.}
\usage{
lavSimulateData(model = NULL, model_type = "sem", meanstructure = FALSE,
    int_ov_free = TRUE, int_lv_free = FALSE,
    marker_int_zero = FALSE, conditional_x = FALSE,
    composites = TRUE, fixed_x = FALSE,
    orthogonal = FALSE, std_lv = TRUE, auto_fix_first = FALSE,
    auto_fix_single = FALSE, auto_var = TRUE, auto_cov_lv_x = TRUE,
    auto_cov_y = TRUE, ..., sample_nobs = 500L, ov_var = NULL,
    group_label = NULL, skewness = NULL,
    kurtosis = NULL, cluster_idx = NULL, seed = NULL, empirical = FALSE,
    mass = FALSE, ordered_center = TRUE, return_type = "data.frame",
    return_fit = FALSE, debug = FALSE, standardized = FALSE)
simulateData(model = NULL, model_type = "sem", meanstructure = FALSE,
    int_ov_free = TRUE, int_lv_free = FALSE,
    marker_int_zero = FALSE, conditional_x = FALSE,
    composites = TRUE, fixed_x = FALSE,
    orthogonal = FALSE, std_lv = TRUE, auto_fix_first = FALSE,
    auto_fix_single = FALSE, auto_var = TRUE, auto_cov_lv_x = TRUE,
    auto_cov_y = TRUE, ..., sample_nobs = 500L, ov_var = NULL,
    group_label = NULL, skewness = NULL,
    kurtosis = NULL, cluster_idx = NULL, seed = NULL, empirical = FALSE,
    mass = FALSE, ordered_center = TRUE, return_type = "data.frame",
    return_fit = FALSE, debug = FALSE, standardized = FALSE)
lav_data_simulate_old(..., ordered_center = FALSE)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (e.g., the output of the \code{lavParTable()} function) is also
    accepted.}
\item{model_type}{Set the model type: possible values
    are \code{"cfa"}, \code{"sem"} or \code{"growth"}. This may affect
    how starting values are computed, and may be used to alter the terminology
    used in the summary output, or the layout of path diagrams that are
    based on a fitted lavaan object.}
\item{meanstructure}{If \code{TRUE}, the means of the observed
    variables enter the model. If \code{"default"}, the value is set based
    on the user-specified model, and/or the values of other arguments.}
\item{int_ov_free}{If \code{FALSE}, the intercepts of the observed variables
    are fixed to zero.}
\item{int_lv_free}{If \code{FALSE}, the intercepts of the latent variables
    are fixed to zero.}
\item{marker_int_zero}{Logical. Only relevant if the metric of each latent
    variable is set by fixing the first factor loading to unity.
    If \code{TRUE}, it implies \code{meanstructure = TRUE} and 
    \code{std_lv = FALSE}, and it fixes the intercepts of the marker
    indicators to zero, while freeing the means/intercepts of the latent
    variables. Only works correctly for single group, single level models.}
\item{conditional_x}{If \code{TRUE}, we set up the model conditional on
    the exogenous `x' covariates; the model-implied sample statistics
    only include the non-x variables. If \code{FALSE}, the exogenous `x'
    variables are modeled jointly with the other variables, and the
    model-implied statistics reflect both sets of variables. If
    \code{"default"}, the value is set depending on the estimator, and
    whether or not the model involves categorical endogenous variables.}
\item{composites}{If \code{TRUE}, use the new (0.6-20) approach
    for handling composites.}
\item{fixed_x}{If \code{TRUE}, the exogenous `x' covariates are considered
    fixed variables and the means, variances and covariances of these variables
    are fixed to their sample values. If \code{FALSE}, they are considered
    random, and the means, variances and covariances are free parameters. If
    \code{"default"}, the value is set depending on the mimic option.}
\item{orthogonal}{If \code{TRUE}, the exogenous latent variables
    are assumed to be uncorrelated.}
\item{std_lv}{If \code{TRUE}, the metric of each latent variable is
    determined by fixing their variances to 1.0. If \code{FALSE}, the metric
    of each latent variable is determined by fixing the factor loading of the
    first indicator to 1.0.}
\item{auto_fix_first}{If \code{TRUE}, the factor loading of the first indicator
    is set to 1.0 for every latent variable.}
\item{auto_fix_single}{If \code{TRUE}, the residual variance (if included)
    of an observed indicator is set to zero if it is the only indicator of a
    latent variable.}
\item{auto_var}{If \code{TRUE}, the (residual) variances of both observed
    and latent variables are set free.}
\item{auto_cov_lv_x}{If \code{TRUE}, the covariances of exogenous latent
    variables are included in the model and set free.}
\item{auto_cov_y}{If \code{TRUE}, the covariances of dependent variables
    (both observed and latent) are included in the model and set free.}
\item{\dots}{additional arguments passed to the \code{\link{lavaan}}
    function.}
\item{sample_nobs}{Number of observations. If a vector, multiple datasets
    are created. If \code{return_type = "matrix"} or 
    \code{return_type = "cov"}, a list of \code{length(sample_nobs)} 
    is returned, with either the data or covariance matrices, each one
    based on the number of observations as specified in \code{sample_nobs}.
    If \code{return_type = "data.frame"}, all datasets are merged and 
    a \code{group} variable is added to mimic a multiple group dataset.}
\item{ov_var}{The user-specified variances of the observed variables.}
\item{group_label}{The group labels that should be used if multiple
    groups are created.}
\item{skewness}{Numeric vector. The skewness values for the observed variables. Defaults to zero.}
\item{kurtosis}{Numeric vector. The kurtosis values for the observed variables. Defaults to zero.}
\item{cluster_idx}{Optional. Only used (and only available via
    \code{lavSimulateData()} and \code{simulateData()}) for multilevel models.
    An integer vector (or a list of integer vectors, one per group) giving the
    cluster membership of each level-1 unit. If supplied, it determines both
    the number of level-1 units and the number of clusters. If not supplied for
    a multilevel model, balanced clusters are created based on \code{sample_nobs}.
    Providing this argument (or a multilevel model syntax with \code{level:}
    blocks) routes the request to the multilevel data-generation engine.}
\item{seed}{Set random seed.}
\item{empirical}{Logical. If \code{TRUE}, the implied moments (Mu and Sigma)
    specify the empirical not population mean and covariance matrix.}
\item{mass}{Logical. If \code{TRUE}, use the mvrnorm() function from the
    MASS package (the default before 0.7-1) instead of the internal
    lav_mvrnorm() function.}
\item{ordered_center}{Logical. Only relevant for categorical data, which is
    generated by the (new) multilevel-aware engine. If \code{TRUE} (the
    default), each ordered variable is centered (its sample mean is subtracted)
    before it is cut at the model-implied thresholds; this also works when the
    model-implied mean is nonzero. If \code{FALSE}, the variable is cut at the
    thresholds directly (the behaviour of the older engine).}
\item{return_type}{If \code{"data.frame"}, a data.frame is returned. If
    \code{"matrix"}, a numeric matrix is returned (without any variable names).
    If \code{"cov"}, a covariance matrix is returned (without any variable 
    names).}
\item{return_fit}{If \code{TRUE}, return the fitted model that has been used
    to generate the data as an attribute (called \code{"fit"}); this 
    may be useful for inspection.}
\item{debug}{If \code{TRUE}, debugging information is displayed.}
\item{standardized}{If \code{TRUE}, the residual variances of the observed
    variables are set in such a way that the model implied variances
    are unity. This allows regression coefficients and factor loadings 
    (involving observed variables) to be specified in a standardized metric.}
}
\details{Model parameters can be specified by fixed values in the lavaan model
syntax. If no fixed values are specified, the value zero will be assumed,
except for factor loadings and variances, which are set to 0.7 and 1.0
respectively. By default, multivariate normal data are generated. However, by
providing skewness and/or kurtosis values, nonnormal multivariate data can be
generated, using the Vale & Maurelli (1983) method.

There is a single data-simulation engine. Multilevel data (model syntax
containing \code{level:} blocks, or when the \code{cluster_idx} argument is
provided) are generated by an internal multilevel worker; all single-level data
(continuous and categorical, including the \code{ov_var}, \code{skewness},
\code{kurtosis}, \code{standardized} and \code{mass} options) are generated by
the historical single-level worker, so that the continuous, single-level output
remains byte-identical to previous versions. For multilevel data, the
\code{ov_var}, \code{skewness}, \code{kurtosis}, \code{standardized} and
\code{mass} arguments are not (yet) supported and are ignored (with a warning).

\code{lav_data_simulate_old()} is a deprecated wrapper kept for backward
compatibility; it forwards to the unified engine (with \code{ordered_center =
FALSE} by default, reproducing the historical threshold cut).}
\value{The generated data. Either as a data.frame 
(if \code{return_type="data.frame"}), 
a numeric matrix (if \code{return_type="matrix"}),
or a covariance matrix (if \code{return_type="cov"}).}
\examples{
# specify population model
population.model <- ' f1 =~ x1 + 0.8*x2 + 1.2*x3
                      f2 =~ x4 + 0.5*x5 + 1.5*x6
                      f3 =~ x7 + 0.1*x8 + 0.9*x9

                      f3 ~ 0.5*f1 + 0.6*f2
                    '

# generate data
set.seed(1234)
myData <- lavSimulateData(population.model, sample_nobs = 100L)

# population moments
fitted(sem(population.model))

# sample moments
round(cov(myData), 3)
round(colMeans(myData), 3)

# fit model
myModel <- ' f1 =~ x1 + x2 + x3
             f2 =~ x4 + x5 + x6
             f3 =~ x7 + x8 + x9
             f3 ~ f1 + f2 '
fit <- sem(myModel, data=myData)
summary(fit)
}
