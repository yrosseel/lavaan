\name{lavInspect}
\alias{lavInspect}
\alias{inspect}
\alias{lavTech}
\title{Inspect or extract information from a fitted lavaan object}
\description{
The \code{lavInspect()} and \code{lavTech()} functions can be used to
inspect/extract information that is stored inside (or can be computed from) a
fitted lavaan object. Note: the (older) \code{inspect()} function is
now simply a shortcut for \code{lavInspect()} with default arguments.
}
\usage{
lavInspect(object, what = "free", add.labels = TRUE, add.class = TRUE,
           list.by.group = TRUE,
           drop.list.single.group = TRUE)

lavTech(object, what = "free", add.labels = FALSE, add.class = FALSE,
        list.by.group = FALSE,
        drop.list.single.group = FALSE)

inspect(object, what = "free", ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{what}{Character. What needs to be inspected/extracted? See Details for a
full list. Note: the \code{what} argument is not case-sensitive
(everything is converted to lower case.)}
\item{add.labels}{If \code{TRUE}, variable names are added to the vectors
and/or matrices.}
\item{add.class}{If \code{TRUE}, vectors are given the `lavaan.vector' class;
matrices are given the `lavaan.matrix' class, and symmetric matrices are
given the `lavaan.matrix.symmetric' class. This only affects the way they
are printed on the screen.}
\item{list.by.group}{Logical. Only used when the output are model matrices.
If \code{TRUE}, the model matrices are nested within groups. If \code{FALSE},
a flattened list is returned containing all model matrices, with repeated
names for multiple groups.}
\item{drop.list.single.group}{If \code{FALSE}, the results are returned as
    a list, where each element corresponds to a group (even if there is only
    a single group). If \code{TRUE}, the list will be unlisted if there is
    only a single group.}
\item{...}{Additional arguments. Not used by lavaan, but by other packages.}
}
\details{
The \code{lavInspect()} and \code{lavTech()} functions only differ in the way
they return the results. The \code{lavInspect()} function will prettify the
output by default, while the \code{lavTech()} will not attempt to prettify the
output by default. The (older) \code{inspect()} function is a simplified
version of \code{lavInspect()} with only the first two arguments.

Below is a list of possible values for the \code{what} argument, organized
in several sections:

Model matrices:

\describe{
    \item{\code{"free"}:}{A list of model matrices. The non-zero integers
        represent the free parameters. The numbers themselves correspond
        to the position of the free parameter in the parameter vector.
        This determines the order of the model parameters in the output
        of for example \code{coef()} and \code{vcov()}.}
    \item{\code{"partable"}:}{A list of model matrices. The non-zero integers
        represent both the fixed parameters (for example, factor loadings
        fixed at 1.0), and the free parameters if we ignore any equality
        constraints. They correspond with all entries (fixed or free)
        in the parameter table. See \code{\link{parTable}}.}
    \item{\code{"se"}:}{A list of model matrices. The non-zero numbers
        represent the standard errors for the free parameters in the model.
        If two parameters are constrained to be equal, they will have the
        same standard error for both parameters.
        Aliases: \code{"std.err"} and \code{"standard.errors"}.}
    \item{\code{"start"}:}{A list of model matrices. The values represent
        the starting values for all model parameters.
        Alias: \code{"starting.values"}.}
    \item{\code{"est"}:}{A list of model matrices. The values represent
        the estimated model parameters. Aliases:
        \code{"estimates"}, and \code{"x"}.}
    \item{\code{"dx.free"}:}{A list of model matrices. The values represent
        the gradient (first derivative) values of the model parameters.
        If two parameters are constrained to be equal, they will have the
        same gradient value.}
    \item{\code{"dx.all"}:}{A list of model matrices. The values represent
        the first derivative with respect to all possible matrix elements.
        Currently, this is only available when the estimator is \code{"ML"}
        or \code{"GLS"}.}
    \item{\code{"std"}:}{A list of model matrices. The values represent
        the (completely) standardized model parameters (the variances of
        both the observed and the latent variables are set to unity).
        Aliases: \code{"std.all"}, \code{"standardized"}.}
    \item{\code{"std.lv"}:}{A list of model matrices. The values represent
        the standardized model parameters (only the variances of the
        latent variables are set to unity.)}
    \item{\code{"std.nox"}:}{A list of model matrices. The values represent
        the (completely) standardized model parameters (the variances of
        both the observed and the latent variables are set to unity;
        however, the variances of any observed exogenous variables are not
        set to unity; hence no-x.)}
}

Information about the data (including missing patterns):

\describe{
    \item{\code{"data"}:}{A matrix containing the observed variables
        that have been used to fit the model. No column/row names are provided.
        Column names correspond to the output of \code{lavNames(object)},
        while the rows correspond to the output of
        \code{lavInspect(object, "case.idx"}.}
    \item{\code{"group"}:}{A character string. The group variable in
        the data.frame (if any).}
    \item{\code{"ngroups"}:}{Integer. The number of groups.}
    \item{\code{"group.label"}:}{A character vector. The group labels.}
    \item{\code{"cluster"}:}{A character vector. The cluster variable(s)
       in the data.frame (if any).}
    \item{\code{"nlevels"}:}{Integer. The number of levels.}
    \item{\code{"ordered"}:}{A character vector. The ordered variables.}
    \item{\code{"nobs"}:}{Integer vector. The number of observations
        in each group that were used in the analysis.}
    \item{\code{"norig"}:}{Integer vector. The original number of observations
        in each group.}
    \item{\code{"ntotal"}:}{Integer. The total number of observations that
        were used in the analysis. If there is just a single group, this
        is the same as the \code{"nobs"} option; if there are multiple groups,
        this is the sum of the \code{"nobs"} numbers for each group.}
    \item{\code{"case.idx"}:}{The case/observation numbers that were used
        in the analysis.
        In the case of multiple groups: a list of numbers.}
    \item{\code{"empty.idx"}:}{The case/observation numbers of those
        cases/observations that contained missing values only
        (at least for the observed variables that were included in the model).
        In the case of multiple groups: a list of numbers.}
    \item{\code{"patterns"}:}{A binary matrix. The rows of the matrix
        are the missing data patterns where 1 and 0 denote non-missing
        and missing values for the corresponding observed variables
        respectively (or
        \code{TRUE} and \code{FALSE} if \code{lavTech()} is used.)
        If the data is complete (no missing values), there will be only
        a single pattern. In the case of multiple groups: a list of
        pattern matrices.}
    \item{\code{"coverage"}:}{A symmetric matrix where each element contains
        the proportion of observed datapoints for the corresponding
        pair of observed variables.
        In the case of multiple groups: a list of coverage matrices.}
}

Observed sample statistics:

\describe{
    \item{\code{"sampstat"}:}{Observed sample statistics. Aliases:
         \code{"obs"}, \code{"observed"}, \code{"samp"}, \code{"sample"},
         \code{"samplestatistics"}. Since
         0.6-3, we always check if an h1 slot is available (the estimates
         for the unrestricted model); if present, we extract the sample
         statistics from this slot. This implies that if variables are
         continuous, and \code{missing = "ml"} (or \code{"fiml"}), we
         return the covariance matrix (and mean vector) as computed by
         the EM algorithm under the unrestricted (h1) model. If the h1 is
         not present (perhaps, because the model was fitted with
         \code{h1 = FALSE}), we return the sample statistics from the
         SampleStats slot. Here, pairwise deletion is used for the elements
         of the covariance matrix (or correlation matrix), and
         listwise deletion for all univariate statistics (means, intercepts
         and thresholds).}
    \item{\code{"sampstat.h1"}:}{Deprecated. Do not use any longer.}
    \item{\code{"wls.obs"}:}{The observed sample statistics (covariance
        elements, intercepts/thresholds, etc.) in a single vector.}
    \item{\code{"wls.v"}:}{The weight vector as used in weighted least
        squares estimation.}
    \item{\code{"gamma"}:}{N times the asymptotic variance matrix of the
        sample statistics. Alias: \code{"sampstat.nacov"}.}
}

Model features:

\describe{
    \item{\code{"meanstructure"}:}{Logical. \code{TRUE} if a meanstructure
        was included in the model.}
    \item{\code{"categorical"}:}{Logical. \code{TRUE} if categorical endogenous
        variables were part of the model.}
    \item{\code{"fixed.x"}:}{Logical. \code{TRUE} if the exogenous x-covariates
        are treated as fixed.}
    \item{\code{"parameterization"}:}{Character. Either \code{"delta"} or
        \code{"theta"}.}
}


Model-implied sample statistics:

\describe{
    \item{\code{"implied"}:}{The model-implied summary statistics.
        Alias: \code{"fitted"}, \code{"expected"}, \code{"exp"}.}
    \item{\code{"resid"}:}{The difference between observed and model-implied 
        summary statistics.
        Alias: \code{"residuals"}, \code{"residual"}, \code{"res"}.}
    \item{\code{"cov.lv"}:}{The model-implied variance-covariance matrix
        of the latent variables. Alias: \code{"veta"} [for V(eta)].}
    \item{\code{"cor.lv"}:}{The model-implied correlation matrix of the
        latent variables.}
    \item{\code{"mean.lv"}:}{The model-implied mean vector of the latent
        variables. Alias: \code{"eeta"} [for E(eta)].}
    \item{\code{"cov.ov"}:}{The model-implied variance-covariance matrix
        of the observed variables.
        Aliases: \code{"sigma"}, \code{"sigma.hat"}.}
    \item{\code{"cor.ov"}:}{The model-implied correlation matrix
        of the observed variables.}
    \item{\code{"mean.ov"}:}{The model-implied mean vector of the observed
        variables. Aliases: \code{"mu"}, \code{"mu.hat"}.}
    \item{\code{"cov.all"}:}{The model-implied variance-covariance matrix
        of both the observed and latent variables.}
    \item{\code{"cor.all"}:}{The model-implied correlation matrix
        of both the observed and latent variables.}
    \item{\code{"th"}:}{The model-implied thresholds.
        Alias: \code{"thresholds"}.}
    \item{\code{"wls.est"}:}{The model-implied sample statistics (covariance
        elements, intercepts/thresholds, etc.) in a single vector.}
    \item{\code{"vy"}:}{The model-implied unconditional variances of the
        observed variables.}
    \item{\code{"rsquare"}:}{The R-square value for all endogenous variables.
        Aliases: \code{"r-square"}, \code{"r2"}.}
}

Optimizer information:

\describe{
    \item{\code{"converged"}:}{Logical. \code{TRUE} if the optimizer has
        converged; \code{FALSE} otherwise.}
    \item{\code{"iteratons"}:}{Integer. The number of iterations used by the
        optimizer.}
    \item{\code{"optim"}:}{List. All available information regarding the
        optimization results.}
    \item{\code{"npar"}:}{Integer. Number of free parameters (ignoring constraints).}
    \item{\code{"coef"}:}{Numeric. The estimated parameter vector.}
}

Gradient, Hessian, observed, expected and first.order information matrices:

\describe{
    \item{\code{"gradient"}:}{Numeric vector containing the first derivatives
        of the discrepancy function with respect to the (free) model
        parameters.}
    \item{\code{"hessian"}:}{Matrix containing the second derivatives
        of the discrepancy function with respect to the (free) model
        parameters.}
    \item{\code{"information"}:}{Matrix containing either the observed or
        the expected information matrix (depending on the information option
        of the fitted model). This is unit-information, not total-information.}
    \item{\code{"information.expected"}:}{Matrix containing the expected
        information matrix for the free model parameters.}
    \item{\code{"information.observed"}:}{Matrix containing the observed
        information matrix for the free model parameters.}
    \item{\code{"information.first.order"}:}{Matrix containing the first.order
        information matrix for the free model parameters. This is the
        outer product of the gradient elements (the first derivative of
        the discrepancy function with respect to the (free) model parameters).
        Alias: \code{"first.order"}.}
    \item{\code{"augmented.information"}:}{Matrix containing either the
        observed or the expected augmented (or bordered) information
        matrix (depending on the information option of the fitted model.
        Only relevant if constraints have been used in the model.}
    \item{\code{"augmented.information.expected"}:}{Matrix containing the
        expected augmented (or bordered) information matrix.
        Only relevant if constraints have been used in the model.}
    \item{\code{"augmented.information.observed"}:}{Matrix containing the
        observed augmented (or bordered) information matrix.
        Only relevant if constraints have been used in the model.}
    \item{\code{"augmented.information.first.order"}:}{Matrix containing
        the first.order augmented (or bordered) information matrix.
        Only relevant if constraints have been used in the model.}
    \item{\code{"inverted.information"}:}{Matrix containing either the
        observed or the expected inverted information matrix
        (depending on the information option of the fitted model.}
    \item{\code{"inverted.information.expected"}:}{Matrix containing the
        inverted expected information matrix for the free model parameters.}
    \item{\code{"inverted.information.observed"}:}{Matrix containing the
        inverted observed information matrix for the free model parameters.}
    \item{\code{"inverted.information.first.order"}:}{Matrix containing the
        inverted first.order information matrix for the free model parameters.}
}

Variance covariance matrix of the model parameters:

\describe{
    \item{\code{"vcov"}:}{Matrix containing the variance covariance matrix
        of the estimated model parameters.}
    \item{\code{"vcov.std.all"}:}{Matrix containing the variance covariance
        matrix of the standardized estimated model parameters. Standardization
        is done with respect to both observed and latent variables.}
    \item{\code{"vcov.std.lv"}:}{Matrix containing the variance covariance
        matrix of the standardized estimated model parameters. Standardization
        is done with respect to the latent variables only.}
    \item{\code{"vcov.std.nox"}:}{Matrix containing the variance covariance
        matrix of the standardized estimated model parameters. Standardization
        is done with respect to both observed and latent variables, but
        ignoring any exogenous observed covariates.}
    \item{\code{"vcov.def"}:}{Matrix containing the variance covariance matrix
        of the user-defined (using the := operator) parameters.}
    \item{\code{"vcov.def.std.all"}:}{Matrix containing the variance covariance
        matrix of the standardized user-defined parameters. Standardization
        is done with respect to both observed and latent variables.}
    \item{\code{"vcov.def.std.lv"}:}{Matrix containing the variance covariance
        matrix of the standardized user-defined parameters. Standardization
        is done with respect to the latent variables only.}
    \item{\code{"vcov.def.std.nox"}:}{Matrix containing the variance covariance
        matrix of the standardized user-defined parameters. Standardization
        is done with respect to both observed and latent variables, but
        ignoring any exogenous observed covariates.}
    \item{\code{"vcov.def.joint"}:}{Matrix containing the joint variance 
        covariance matrix of both the estimated model parameters and
        the defined (using the := operator) parameters.}
    \item{\code{"vcov.def.joint.std.all"}:}{Matrix containing the joint
        variance covariance matrix of both the standardized model parameters
        and the user-defined parameters. Standardization
        is done with respect to both observed and latent variables.}
    \item{\code{"vcov.def.joint.std.lv"}:}{Matrix containing the joint 
        variance covariance matrix of both the standardized model parameters
        and the user-defined parameters. Standardization
        is done with respect to the latent variables only.}   
    \item{\code{"vcov.def.joint.std.nox"}:}{Matrix containing the joint 
        variance covariance matrix of both the standardized model parameters
        and the user-defined parameters. Standardization
        is done with respect to both observed and latent variables, but
        ignoring any exogenous observed covariates.}
}

Miscellaneous:

\describe{
    \item{\code{"UGamma"}:}{Matrix containing the product of 'U' and 'Gamma'
        matrices as used by the Satorra-Bentler correction. The trace of
        this matrix, divided by the degrees of freedom, gives the scaling
        factor.}
    \item{\code{"UfromUGamma"}:}{Matrix containing the 'U' matrix
        as used by the Satorra-Bentler correction. Alias: \code{"U"}.}
    \item{\code{"list"}:}{The parameter table. The same output as given
        by \code{parTable()}.}
    \item{\code{"fit"}:}{The fit measures. Aliases: \code{"fitmeasures"},
        \code{"fit.measures"}, \code{"fit.indices"}. The same output as
        given by \code{fitMeasures()}.}
    \item{\code{"mi"}:}{The modification indices. Alias: \code{"modindices"},
        \code{"modification.indices"}. The same output as given
        by \code{modindices()}.}
    \item{\code{"options"}:}{List. The option list.}
    \item{\code{"call"}:}{List. The call as returned by match.call, coerced to
        a list.}
    \item{\code{"timing"}:}{List. The timing (in milliseconds) of various
        lavaan subprocedures.}
    \item{\code{"test"}:}{List. All available information regarding the
        (goodness-of-fit) test statistic(s).}
    \item{\code{"post.check"}:}{Post-fitting check if the solution is
        admissible. A warning is raised if negative variances are found, or if
        either \code{lavInspect(fit, "cov.lv")} or
        \code{lavInspect(fit, "theta")} return a non-positive definite matrix.}
    \item{\code{"zero.cell.tables"}:}{List. List of bivariate frequency tables
        where at least one cell is empty.}
    \item{\code{"version"}:}{The lavaan version number that was used to 
        construct the fitted lavaan object.}
}

}
\seealso{
\code{\link{lavaan}}
}
\examples{
# fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

     
fit <- cfa(HS.model, data = HolzingerSwineford1939, group = "school")

# extract information
lavInspect(fit, "sampstat")
lavTech(fit, "sampstat")
}
