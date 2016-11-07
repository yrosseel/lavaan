\name{simulateData}
\alias{simulateData}
\title{Simulate Data From a Lavaan Model Syntax}
\description{Simulate data starting from a lavaan model syntax.}
\usage{
simulateData(model = NULL, model.type = "sem", meanstructure = FALSE, 
    int.ov.free = TRUE, int.lv.free = FALSE, conditional.x = FALSE,
    fixed.x = FALSE, 
    orthogonal = FALSE, std.lv = TRUE, auto.fix.first = FALSE, 
    auto.fix.single = FALSE, auto.var = TRUE, auto.cov.lv.x = TRUE, 
    auto.cov.y = TRUE, ..., sample.nobs = 500L, ov.var = NULL, 
    group.label = paste("G", 1:ngroups, sep = ""), skewness = NULL, 
    kurtosis = NULL, seed = NULL, empirical = FALSE, 
    return.type = "data.frame", return.fit = FALSE,
    debug = FALSE, standardized = FALSE)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (eg. the output of the \code{lavaanify()} function) is also
    accepted.}
\item{model.type}{Set the model type: possible values
    are \code{"cfa"}, \code{"sem"} or \code{"growth"}. This may affect
    how starting values are computed, and may be used to alter the terminology
    used in the summary output, or the layout of path diagrams that are
    based on a fitted lavaan object.}
\item{meanstructure}{If \code{TRUE}, the means of the observed
    variables enter the model. If \code{"default"}, the value is set based
    on the user-specified model, and/or the values of other arguments.}
\item{int.ov.free}{If \code{FALSE}, the intercepts of the observed variables
    are fixed to zero.}
\item{int.lv.free}{If \code{FALSE}, the intercepts of the latent variables
    are fixed to zero.}
\item{conditional.x}{If \code{TRUE}, we set up the model conditional on
    the exogenous `x' covariates; the model-implied sample statistics
    only include the non-x variables. If \code{FALSE}, the exogenous `x'
    variables are modeled jointly with the other variables, and the
    model-implied statistics refect both sets of variables. If
    \code{"default"}, the value is set depending on the estimator, and
    whether or not the model involves categorical endogenous variables.}
\item{fixed.x}{If \code{TRUE}, the exogenous `x' covariates are considered
    fixed variables and the means, variances and covariances of these variables
    are fixed to their sample values. If \code{FALSE}, they are considered
    random, and the means, variances and covariances are free parameters. If
    \code{"default"}, the value is set depending on the mimic option.}
\item{orthogonal}{If \code{TRUE}, the exogenous latent variables
    are assumed to be uncorrelated.}
\item{std.lv}{If \code{TRUE}, the metric of each latent variable is
    determined by fixing their variances to 1.0. If \code{FALSE}, the metric
    of each latent variable is determined by fixing the factor loading of the
    first indicator to 1.0.}
\item{auto.fix.first}{If \code{TRUE}, the factor loading of the first indicator
    is set to 1.0 for every latent variable.}
\item{auto.fix.single}{If \code{TRUE}, the residual variance (if included)
    of an observed indicator is set to zero if it is the only indicator of a
    latent variable.}
\item{auto.var}{If \code{TRUE}, the residual variances and the variances
    of exogenous latent variables are included in the model and set free.}
\item{auto.cov.lv.x}{If \code{TRUE}, the covariances of exogenous latent
    variables are included in the model and set free.}
\item{auto.cov.y}{If \code{TRUE}, the covariances of dependent variables
    (both observed and latent) are included in the model and set free.}
\item{\dots}{additional arguments passed to the \code{\link{lavaan}}
    function.}
\item{sample.nobs}{Number of observations. If a vector, multiple datasets
    are created. If \code{return.type = "matrix"} or 
    \code{return.type = "cov"}, a list of \code{length(sample.nobs)} 
    is returned, with either the data or covariance matrices, each one
    based on the number of observations as specified in \code{sample.nobs}.
    If \code{return.type = "data.frame"}, all datasets are merged and 
    a \code{group} variable is added to mimic a multiple group dataset.}
\item{ov.var}{The user-specified variances of the observed variables.}
\item{group.label}{The group labels that should be used if multiple
    groups are created.}
\item{skewness}{Numeric vector. The skewness values for the observed variables. Defaults to zero.}
\item{kurtosis}{Numeric vector. The kurtosis values for the observed variables. Defaults to zero.}
\item{seed}{Set random seed.}
\item{empirical}{Logical. If \code{TRUE}, the implied moments (Mu and Sigma)
    specify the empirical not population mean and covariance matrix.}
\item{return.type}{If \code{"data.frame"}, a data.frame is returned. If
    \code{"matrix"}, a numeric matrix is returned (without any variable names).
    If \code{"cov"}, a covariance matrix is returned (without any variable 
    names).}
\item{return.fit}{If \code{TRUE}, return the fitted model that has been used
    to generate the data as an attribute (called \code{"fit"}); this 
    may be useful for inspection.}
\item{debug}{If \code{TRUE}, debugging information is displayed.}
\item{standardized}{If \code{TRUE}, the residual variances of the observed
    variables are set in such a way such that the model implied variances 
    are unity. This allows regression coefficients and factor loadings 
    (involving observed variables) to be specified in a standardized metric.}
}
\details{Model parameters can be specified by fixed values in the lavaan
model syntax. If no fixed values are specified, the value zero will be 
assumed, except for factor loadings and variances, which are set to unity 
by default. By default, multivariate normal data are generated. However,
by providing skewness and/or kurtosis values, nonnormal multivariate data
can be generated, using the Vale & Maurelli (1983) method.}
\value{The generated data. Either as a data.frame 
(if \code{return.type="data.frame"}), 
a numeric matrix (if \code{return.type="matrix"}),
or a covariance matrix (if \code{return.type="cov"}).}
\examples{
# specify population model
population.model <- ' f1 =~ x1 + 0.8*x2 + 1.2*x3
                      f2 =~ x4 + 0.5*x5 + 1.5*x6
                      f3 =~ x7 + 0.1*x8 + 0.9*x9

                      f3 ~ 0.5*f1 + 0.6*f2
                    '

# generate data
set.seed(1234)
myData <- simulateData(population.model, sample.nobs=100L)

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
