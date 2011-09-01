\name{lavaan}
\alias{lavaan}
\title{Fit a Latent Variable Model}
\description{
Fit a latent variable model.}
\usage{
lavaan(model.syntax = "", model.type = "sem", meanstructure = "default", 
    int.ov.free = FALSE, int.lv.free = FALSE, fixed.x = "default", 
    orthogonal = FALSE, std.lv = FALSE, auto.fix.first = FALSE, 
    auto.fix.single = FALSE, auto.var = FALSE, auto.cov.lv.x = FALSE, 
    auto.cov.y = FALSE, data = NULL, std.ov = FALSE, missing = "default", 
    sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL, 
    group = NULL, group.equal = "", group.partial = "",
    constraints = '', estimator = "default", likelihood = "default", 
    information = "default", se = "default",  test = "default", 
    bootstrap = 1000L, mimic = "default", representation = "default", 
    do.fit = TRUE, control = list(), start = "default",
    verbose = FALSE, warn = TRUE, debug = FALSE)
}
\arguments{
\item{model.syntax}{The model syntax specifying the model. See 
    \code{\link{model.syntax}} for more information. Or a parsed model
    list as outputted by the lavaanify() function.}
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
\item{data}{An optional data frame containing the observed variables used in
    the model.}
\item{std.ov}{If \code{TRUE}, all observed variables are standardized
    before entering the analysis.}
\item{missing}{If \code{"listwise"}, cases with missing values are removed
    listwise from the data frame before analysis. If \code{"direct"} or
    \code{"ml"} or \code{"fiml"} and the estimator is maximum likelihood,
    Full Information Maximum Likelihood (FIML) estimation is used using all
    available data in the data frame. This is only valid if the data are
    missing completely at random (MCAR) or missing at random (MAR). If
    \code{"default"}, the value is set depending on the estimator and the
    mimic option.}
\item{sample.cov}{Numeric matrix. A sample variance-covariance matrix.
    The rownames must contain the observed variable names.
    For a multiple group analysis, a list with a variance-covariance matrix
    for each group.}
\item{sample.mean}{A sample mean vector. For a multiple group analysis,
    a list with a mean vector for each group.}
\item{sample.nobs}{Number of observations if the full data frame is missing
    and only sample moments are given. For a multiple group analysis, a list
    or a vector with the number of observations for each group.}
\item{group}{A variable name in the data frame defining the groups in a
    multiple group analysis.}
\item{group.equal}{A vector of character strings. Only used in
    a multiple group analysis. Can be one or more of the following: 
    \code{"loadings"}, \code{"intercepts"}, \code{"means"}, 
    \code{"regressions"}, \code{"residuals"},
    \code{"residual.covariances"}, \code{"lv.variances"} or
    \code{"lv.covariances"}, specifying the pattern of equality
    constraints across multiple groups.}
\item{group.partial}{A vector of character strings containing the labels
    of the parameters which should be free in all groups (thereby
    overriding the group.equal argument for some specific parameters).}
\item{constraints}{Additional (in)equality constraints not yet included in the
    model syntax. See \code{\link{model.syntax}} for more information.}
\item{estimator}{The estimator to be used. Can be one of the following:
    \code{"ML"} for maximum likelihood, \code{"GLS"} for generalized least 
    squares, \code{"WLS"} for weighted least squares (sometimes called ADF 
    estimation), \code{"MLM"} for maximum likelihood estimation with robust
    standard errors and a Satorra-Bentler scaled test statistic,
    \code{"MLF"} for maximum likelihood estimation with standard errors 
    based on first-order derivatives and a conventional test statistic,
    \code{"MLR"} for maximum likelihood estimation with robust `Huber-White' 
    standard errors and a scaled test statistic which is asymptotically 
    equivalent to the Yuan-Bentler T2-star test statistic. Note that the 
    \code{"MLM"}, \code{"MLF"} and \code{"MLR"} choices only affect the
    standard errors and the test statistic. They also imply 
    \code{mimic="Mplus"}.}
\item{likelihood}{Only relevant for ML estimation. If \code{"wishart"},
    the wishart likelihood approach is used. In this approach, the covariance
    matrix has been divided by N-1, and both standard errors and test 
    statistics are based on N-1. 
    If \code{"normal"}, the normal likelihood approach is used. Here,
    the covariance matrix has been divided by N, and both standard errors 
    and test statistics are based on N. If \code{"default"}, it depends
    on the mimic option: if \code{mimic="Mplus"}, normal likelihood is used;
    otherwise, wishart likelihood is used.}
\item{information}{If \code{"expected"}, the expected information matrix
    is used (to compute the standard errors). If \code{"observed"}, the
    observed information matrix is used. If \code{"default"}, the value is 
    set depending on the estimator and the mimic option.}
\item{se}{If \code{"standard"}, conventional standard errors 
    are computed based on inverting the (expected or observed) information 
    matrix. If \code{"first.order"}, standard errors are computed based on 
    first-order derivatives. If \code{"robust.mlm"}, conventional robust 
    standard errors are computed.
    If \code{"robust.mlr"},
    standard errors are computed based on the `mlr' (aka pseudo ML, 
    Huber-White) approach.
    If \code{"robust"}, either \code{"robust.mlm"} or \code{"robust.mlr"} is
    used depending on the estimator, the mimic option, and whether the data
    are complete or not. 
    If \code{"boot"} or \code{"bootstrap"}, bootstrap standard errors are
    computed using standard bootstrapping.
    If \code{"none"}, no standard errors are computed.}
\item{test}{If \code{"standard"}, a conventional chi-square test is computed.
    If \code{"Satorra-Bentler"}, a Satorra-Bentler scaled test statistic is
    computed. If \code{"Yuan-Bentler"}, a Yuan-Bentler scaled test statistic
    is computed. If \code{"default"}, the value depends on the
    values of other arguments.}
\item{bootstrap}{Number of bootstrap draws, if bootstrapping is used.}
\item{mimic}{If \code{"Mplus"}, an attempt is made to mimic the Mplus 
    program. If \code{"EQS"}, an attempt is made to mimic the EQS program.
    If \code{"default"}, the value is (currently) set to \code{"Mplus"}.}
\item{representation}{If \code{"LISREL"} the classical LISREL matrix 
    representation is used to represent the model (using the all-y variant).}
\item{do.fit}{If \code{FALSE}, the model is not fit, and the current
    starting values of the model parameters are preserved.}
\item{control}{A list containing control parameters passed to the optimizer. 
    By default, lavaan uses \code{"nlminb"}. See the manpage of 
    \code{\link{nlminb}} for an overview of the control parameters. 
    A different optimizer can be chosen by setting the value of 
    \code{optim.method}. For unconstrained optimization (the model syntax 
    does not include any "==", ">" or "<" operators),
    the available options are \code{"nlminb"} (the default), \code{"BFGS"} and 
    \code{"L-BFGS-B"}. See the manpage of the \code{\link{optim}} function for
    the control parameters of the latter two options. For constrained
    optimization, the only available option is \code{"nlminb.constr"}.}
\item{start}{If it is a character string,
    the two options are currently \code{"simple"} and \code{"Mplus"}. 
    In the first
    case, all parameter values are set to zero, except the factor loadings
    (set to one), the variances of latent variables (set to 0.05), and
    the residual variances of observed variables (set to half the observed 
    variance).
    If \code{"Mplus"}, we use a similar scheme, but the factor loadings are
    estimated using the fabin3 estimator (tsls) per factor. 
    If \code{start} is a fitted
    object of class \code{\linkS4class{lavaan}}, the estimated values of
    the corresponding parameters will be extracted. If it is a model list, 
    for example the output of the \code{paramaterEstimates()} function, 
    the values of the \code{est} or \code{start} or \code{ustart} column 
    (whichever is found first) will be extracted.}
\item{verbose}{If \code{TRUE}, the function value is printed out during
    each iteration.}
\item{warn}{If \code{TRUE}, some (possibly harmless) warnings are printed
    out during the iterations.}
\item{debug}{If \code{TRUE}, debugging information is printed out.}
}
\value{
  An object of class \code{\linkS4class{lavaan}}, for which several methods
  are available, including a \code{summary} method.
}
\seealso{
\code{\link{cfa}}, \code{\link{sem}}, \code{\link{growth}}
}
\examples{
# The Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=TRUE,
              auto.cov.lv.x=TRUE)
summary(fit, fit.measures=TRUE)
}
