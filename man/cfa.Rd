\name{cfa}
\alias{cfa}
\title{Fit Confirmatory Factor Analysis Models}
\description{
Fit a Confirmatory Factor Analysis (CFA) model.}
\usage{
cfa(model = NULL, data = NULL,
    meanstructure = "default", 
    conditional.x = "default", fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, 
    parameterization = "default", std.ov = FALSE,
    missing = "default", ordered = NULL, 
    sample.cov = NULL, sample.cov.rescale = "default",
    sample.mean = NULL, sample.nobs = NULL, 
    ridge = 1e-05, group = NULL, 
    group.label = NULL, group.equal = "", group.partial = "", 
    group.w.free = FALSE, cluster = NULL, constraints = '', 
    estimator = "default", likelihood = "default", link = "default",
    information = "default", se = "default", test = "default",
    bootstrap = 1000L, mimic = "default", representation = "default", 
    do.fit = TRUE, control = list(), WLS.V = NULL, NACOV = NULL,
    zero.add = "default", zero.keep.margins = "default",
    zero.cell.warn = TRUE,
    start = "default", verbose = FALSE, warn = TRUE, debug = FALSE)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See                 
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (eg. the output of the \code{lavaanify()} function) is also
    accepted.}
\item{data}{An optional data frame containing the observed variables used in
    the model. If some variables are declared as ordered factors, lavaan will
    treat them as ordinal variables.}
\item{meanstructure}{If \code{TRUE}, the means of the observed
    variables enter the model. If \code{"default"}, the value is set based
    on the user-specified model, and/or the values of other arguments.}
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
    determined by fixing their (residual) variances to 1.0. If \code{FALSE}, 
    the metric
    of each latent variable is determined by fixing the factor loading of the
    first indicator to 1.0.}
\item{parameterization}{Currently only used if data is categorical. If 
    \code{"delta"}, the delta parameterization is used. If \code{"theta"},
    the theta parameterization is used.}
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
\item{ordered}{Character vector. Only used if the data is in a data.frame. 
    Treat these variables as ordered (ordinal) variables, if they are 
    endogenous in the model. Importantly, all other variables will be treated 
    as numeric (unless they are declared as ordered in the original data.frame.)}
\item{sample.cov}{Numeric matrix. A sample variance-covariance matrix.
    The rownames and/or colnames must contain the observed variable names.
    For a multiple group analysis, a list with a variance-covariance matrix
    for each group. Note that if maximum likelihood estimation is used and
    \code{likelihood="normal"}, the user provided covariance matrix is
    internally rescaled by multiplying it with a factor (N-1)/N, to ensure
    that the covariance matrix has been divided by N. This can be turned off
    by setting the \code{sample.cov.rescale} argument to \code{FALSE}.}
\item{sample.cov.rescale}{If \code{TRUE}, the sample covariance matrix provided
    by the user is internally rescaled by multiplying it with a factor (N-1)/N.
    If \code{"default"}, the value is set depending on the estimator and the
    likelihood option: it is set to \code{TRUE} if maximum likelihood 
    estimation is used and \code{likelihood="normal"}, and \code{FALSE} 
    otherwise.}
\item{sample.mean}{A sample mean vector. For a multiple group analysis,
    a list with a mean vector for each group.}
\item{sample.nobs}{Number of observations if the full data frame is missing
    and only sample moments are given. For a multiple group analysis, a list
    or a vector with the number of observations for each group.}
\item{ridge}{Numeric. Small constant used for ridging. Only used if the sample covariance matrix is non positive definite.}
\item{group}{A variable name in the data frame defining the groups in a
    multiple group analysis.}
\item{group.label}{A character vector. The user can specify which group (or
factor) levels need to be selected from the grouping variable, and in which
order. If \code{NULL} (the default), all grouping levels are selected, in the
order as they appear in the data.}
\item{group.equal}{A vector of character strings. Only used in
    a multiple group analysis. Can be one or more of the following:
    \code{"loadings"}, \code{"intercepts"}, \code{"means"}, \code{"thresholds"},
    \code{"regressions"}, \code{"residuals"},
    \code{"residual.covariances"}, \code{"lv.variances"} or
    \code{"lv.covariances"}, specifying the pattern of equality
    constraints across multiple groups.}
\item{group.partial}{A vector of character strings containing the labels
    of the parameters which should be free in all groups (thereby
    overriding the group.equal argument for some specific parameters).}
\item{group.w.free}{Logical. If \code{TRUE}, the group frequencies are
    considered to be free parameters in the model. In this case, a 
    Poisson model is fitted to estimate the group frequencies. If
    \code{FALSE} (the default), the group frequencies are fixed to their
    observed values.}
\item{cluster}{Not used yet.}
\item{constraints}{Additional (in)equality constraints not yet included in the
    model syntax. See \code{\link{model.syntax}} for more information.}
\item{estimator}{The estimator to be used. Can be one of the following:
    \code{"ML"} for maximum likelihood, \code{"GLS"} for generalized least 
    squares, \code{"WLS"} for weighted least squares (sometimes called ADF 
    estimation), \code{"ULS"} for unweighted least squares and \code{"DWLS"} for
    diagonally weighted least squares. These are the main options that affect
    the estimation. For convenience, the \code{"ML"} option can be extended
    as \code{"MLM"}, \code{"MLMV"}, \code{"MLMVS"}, \code{"MLF"}, and 
    \code{"MLR"}. The estimation will still be plain \code{"ML"}, but now
    with robust standard errors and a robust (scaled) test statistic. For
    \code{"MLM"}, \code{"MLMV"}, \code{"MLMVS"}, classic robust standard 
    errors are used (\code{se="robust.sem"}); for \code{"MLF"}, standard
    errors are based on first-order derivatives (\code{se="first.order"});
    for \code{"MLR"}, `Huber-White' robust standard errors are used 
    (\code{se="robust.huber.white"}). In addition, \code{"MLM"} will compute
    a Satorra-Bentler scaled (mean adjusted) test statistic 
    (\code{test="satorra.bentler"}) , \code{"MLMVS"} will compute a 
    mean and variance adjusted test statistic (Satterthwaite style) 
    (\code{test="mean.var.adjusted"}), \code{"MLMV"} will compute a mean 
    and variance adjusted test statistic (scaled and shifted) 
    (\code{test="scaled.shifted"}), and \code{"MLR"} will
    compute a test statistic which is asymptotically 
    equivalent to the Yuan-Bentler T2-star test statistic. Analogously,
    the estimators \code{"WLSM"} and \code{"WLSMV"} imply the \code{"DWLS"}
    estimator (not the \code{"WLS"} estimator) with robust standard errors
    and a mean or mean and variance adjusted test statistic. Estimators
    \code{"ULSM"} and \code{"ULSMV"} imply the \code{"ULS"}
    estimator with robust standard errors
    and a mean or mean and variance adjusted test statistic.}  
\item{likelihood}{Only relevant for ML estimation. If \code{"wishart"},
    the wishart likelihood approach is used. In this approach, the covariance
    matrix has been divided by N-1, and both standard errors and test 
    statistics are based on N-1. 
    If \code{"normal"}, the normal likelihood approach is used. Here,
    the covariance matrix has been divided by N, and both standard errors 
    and test statistics are based on N. If \code{"default"}, it depends
    on the mimic option: if \code{mimic="lavaan"} or \code{mimic="Mplus"}, 
    normal likelihood is used; otherwise, wishart likelihood is used.}
\item{link}{Currently only used if estimator is MML. If \code{"logit"},
    a logit link is used for binary and ordered observed variables.
    If \code{"probit"}, a probit link is used. If \code{"default"}, 
    it is currently set to \code{"probit"} (but this may change).}
\item{information}{If \code{"expected"}, the expected information matrix
    is used (to compute the standard errors). If \code{"observed"}, the
    observed information matrix is used. If \code{"default"}, the value is 
    set depending on the estimator and the mimic option.}
\item{se}{If \code{"standard"}, conventional standard errors 
    are computed based on inverting the (expected or observed) information 
    matrix. If \code{"first.order"}, standard errors are computed based on 
    first-order derivatives. If \code{"robust.sem"}, conventional robust 
    standard errors are computed.  If \code{"robust.huber.white"},
    standard errors are computed based on the `mlr' (aka pseudo ML, 
    Huber-White) approach.
    If \code{"robust"}, either \code{"robust.sem"} or 
    \code{"robust.huber.white"} is used depending on the estimator, 
    the mimic option, and whether the data are complete or not. 
    If \code{"boot"} or \code{"bootstrap"}, bootstrap standard errors are
    computed using standard bootstrapping (unless Bollen-Stine bootstrapping
    is requested for the test statistic; in this case bootstrap standard 
    errors are computed using model-based bootstrapping).
    If \code{"none"}, no standard errors are computed.}
\item{test}{If \code{"standard"}, a conventional chi-square test is computed.
    If \code{"Satorra.Bentler"}, a Satorra-Bentler scaled test statistic is
    computed. If \code{"Yuan.Bentler"}, a Yuan-Bentler scaled test statistic
    is computed. If \code{"mean.var.adjusted"} or \code{"Satterthwaite"}, a 
    mean and variance adjusted test statistic is compute. 
    If \code{"scaled.shifted"}, an alternative mean and variance adjusted test
    statistic is computed (as in Mplus version 6 or higher). 
    If \code{"boot"} or \code{"bootstrap"} or 
    \code{"Bollen.Stine"}, the Bollen-Stine bootstrap is used to compute
    the bootstrap probability value of the test statistic.
    If \code{"default"}, the value depends on the
    values of other arguments.}
\item{bootstrap}{Number of bootstrap draws, if bootstrapping is used.}
\item{mimic}{If \code{"Mplus"}, an attempt is made to mimic the Mplus 
    program. If \code{"EQS"}, an attempt is made to mimic the EQS program.
    If \code{"default"}, the value is (currently) set to \code{"lavaan"},
    which is very close to \code{"Mplus"}.}
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
\item{WLS.V}{A user provided weight matrix to be used by estimator \code{"WLS"};
    if the estimator is \code{"DWLS"}, only the diagonal of this matrix will be 
    used. For a multiple group analysis, a list with a weight matrix
    for each group. The elements of the weight matrix should be in the
    following order (if all data is continuous): first the means (if a 
    meanstructure is involved), then the lower triangular elements of the 
    covariance matrix including the diagonal, ordered column by column. In
    the categorical case: first the thresholds (including the means for 
    continuous variables), then the slopes (if any), the variances of 
    continuous variables (if any), and finally the lower triangular elements 
    of the correlation/covariance matrix excluding the diagonal, ordered 
    column by column.}
\item{NACOV}{A user provided matrix containing the elements of (N times)
    the asymptotic variance-covariance matrix of the sample statistics.
    For a multiple group analysis, a list with an asymptotic 
    variance-covariance matrix for each group. See the \code{WLS.V} 
    argument for information about the order of the elements.}
\item{zero.add}{A numeric vector containing two values. These values affect the
    calculation of polychoric correlations when some frequencies in the
    bivariate table are zero.
    The first value only applies for 2x2 tables. The second value for larger
    tables. This value is added to the zero frequency in the bivariate table.
    If \code{"default"}, the value is set depending on the \code{"mimic"}
    option. By default, lavaan uses \code{zero.add = c(0.5. 0.0)}.}
\item{zero.keep.margins}{Logical. This argument only affects the computation
    of polychoric correlations for 2x2 tables with an empty cell, and where a
    value is added to the empty cell. If \code{TRUE}, the other values of the
    frequency table are adjusted so that all margins are unaffected. If
    \code{"default"}, the value is set depending on the \code{"mimic"}. The
    default is \code{TRUE}.}
\item{zero.cell.warn}{Logical. Only used if some observed endogenous variables
    are categorical. If \code{TRUE}, give a warning if one or more cells
    of a bivariate frequency table are empty.}
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
\details{
    The \code{cfa} function is a wrapper for the more general 
    \code{\link{lavaan}} function, using the following default arguments:
    \code{int.ov.free = TRUE}, \code{int.lv.free = FALSE}, 
    \code{auto.fix.first = TRUE} (unless \code{std.lv = TRUE}), 
    \code{auto.fix.single = TRUE}, \code{auto.var = TRUE}, 
    \code{auto.cov.lv.x = TRUE}, 
    \code{auto.th = TRUE}, \code{auto.delta = TRUE},
    and \code{auto.cov.y = TRUE}.
}
\value{
  An object of class \code{\linkS4class{lavaan}}, for which several methods
  are available, including a \code{summary} method.
}
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. URL
http://www.jstatsoft.org/v48/i02/.}
\seealso{
\code{\link{lavaan}}
}
\examples{
## The famous Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
summary(fit, fit.measures=TRUE)
}
