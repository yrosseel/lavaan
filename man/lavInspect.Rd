\name{lavInspect}
\alias{lavInspect}
\alias{inspect}
\alias{lavTech}
\alias{lavInspect.lavaan}
\alias{inspect.lavaan}
\alias{lavTech.lavaan}
\title{Inspect or extract information from a fitted lavaan object}
\description{
The \code{lavInspect()} and \code{lavTech()} functions can be used to
inspect/extract information that is stored inside (or can be computed from) a
fitted lavaan object. Note: the (older) \code{inspect()} function is
now simply a shortcut for \code{lavInspect()} with default arguments.
}
\usage{
\method{lavInspect}{lavaan}(object, what = "free", add_labels = TRUE, add_class = TRUE,
           list_by_group = TRUE,
           drop_list_single_group = TRUE, ...)

\method{lavTech}{lavaan}(object, what = "free", add_labels = FALSE, add_class = FALSE,
        list_by_group = FALSE,
        drop_list_single_group = FALSE, ...)

\method{inspect}{lavaan}(object, what = "free", ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{what}{Character. What needs to be inspected/extracted? See Details for a
full list. Note: the \code{what} argument is not case-sensitive
(everything is converted to lower case.)}
\item{add_labels}{If \code{TRUE}, variable names are added to the vectors
and/or matrices.}
\item{add_class}{If \code{TRUE}, vectors are given the `lavaan.vector' class;
matrices are given the `lavaan.matrix' class, and symmetric matrices are
given the `lavaan.matrix.symmetric' class. This only affects the way they
are printed on the screen.}
\item{list_by_group}{Logical. Only used when the output are model matrices.
If \code{TRUE}, the model matrices are nested within groups. If \code{FALSE},
a flattened list is returned containing all model matrices, with repeated
names for multiple groups.}
\item{drop_list_single_group}{If \code{FALSE}, the results are returned as
    a list, where each element corresponds to a group (even if there is only
    a single group). If \code{TRUE}, the list will be unlisted if there is
    only a single group.}
\item{...}{For \code{inspect()} additional arguments can be specified, not used by lavaan, but by other packages.
For the functions \code{lavInspect()} and \code{lavTech()} this argument is used to support old names of arguments.}
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
        of, for example, \code{coef()} and \code{vcov()}.}
    \item{\code{"partable"}:}{A list of model matrices. The non-zero integers
        represent both the fixed parameters (for example, factor loadings
        fixed at 1.0) and the free parameters (ignoring any equality
        constraints). They correspond to all entries (fixed or free)
        in the parameter table. See \code{\link{parTable}}.}
    \item{\code{"se"}:}{A list of model matrices. The non-zero numbers
        represent the standard errors for the free parameters in the model.
        If two parameters are constrained to be equal, they will share the
        same standard error.
        Aliases: \code{"std.err"} and \code{"standard.errors"}.}
    \item{\code{"se.std"}:}{A list of model matrices. The non-zero numbers
        represent the standard errors for the (completely) standardized free
        parameters in the model. Alias: \code{"std.se"}.}
    \item{\code{"start"}:}{A list of model matrices. The values represent
        the starting values for all model parameters.
        Alias: \code{"starting.values"}.}
    \item{\code{"est"}:}{A list of model matrices. The values represent
        the estimated model parameters. Aliases:
        \code{"estimates"}, and \code{"x"}.}
    \item{\code{"est.unrotated"}:}{A list of model matrices. The values
        represent the estimated model parameters before rotation was applied.
        Only relevant if rotation was used (for example in EFA); otherwise
        this is identical to \code{"est"}.}
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
    \item{\code{"mm.lambda"}:}{The estimated factor loading matrix (per
        group), as extracted from the \code{@GLIST} slot.
        Alias: \code{"lambda"}.}
    \item{\code{"mm.delta"}:}{The estimated delta (scaling) matrix (per group),
        as extracted from the \code{@GLIST} slot. Only available for
        categorical data or when a correlation structure is used. (Note: this
        is the delta model matrix, not the delta \emph{jacobian} that is
        returned by \code{"delta"}.)}
}

To extract a single model matrix (from the \code{@GLIST} slot), use the
\code{"mm.*"} options (for example \code{"mm.lambda"} for the factor loading
matrix). For backward compatibility, the (older) bare names are still accepted
as aliases.

Information about the data:

\describe{
    \item{\code{"data"}:}{A matrix containing the observed variables
        that were used to fit the model. The matrix has no row or column
        names. The columns correspond to the output of
        \code{lav_object_vnames(object)}, while the rows correspond to the
        output of \code{lavInspect(object, "case.idx")}.}
    \item{\code{"ordered"}:}{A character vector. The ordered variables.}
    \item{\code{"nobs"}:}{Numeric vector. The effective number of observations
        in each group that were used in the analysis. If sampling weights were
        provided, this is the (per group) sum of the weights, consistent with
        the value returned by the \code{nobs()} method. Use \code{"norig"} to
        obtain the number of rows in the data.}
    \item{\code{"norig"}:}{Integer vector. The original number of observations
        (rows in the data) in each group.}
    \item{\code{"ntotal"}:}{Numeric. The total effective number of observations
        that were used in the analysis. If there is just a single group, this
        is the same as the \code{"nobs"} option; if there are multiple groups,
        this is the sum of the \code{"nobs"} numbers for each group. If sampling
        weights were provided, this is the sum of the weights.}
    \item{\code{"case.idx"}:}{Integer vector. The case/observation numbers
        that were used in the analysis.
        In the case of multiple groups: a list of numbers.}
    \item{\code{"empty.idx"}:}{The case/observation numbers of those
        cases/observations that contained missing values only
        (at least for the observed variables that were included in the model).
        In the case of multiple groups: a list of numbers.}
    \item{\code{"th.idx"}:}{Integer vector. For categorical data, the index of
        the observed variable that each threshold belongs to.
        In the case of multiple groups: a list of integer vectors.}
    \item{\code{"patterns"}:}{A binary matrix. The rows of the matrix
        are the missing data patterns where 1 and 0 denote non-missing
        and missing values for the corresponding observed variables
        respectively (or
        \code{TRUE} and \code{FALSE} if \code{lavTech()} is used.)
        If the data is complete (no missing values), there will be only
        a single pattern. In the case of multiple groups: a list of
        pattern matrices.}
    \item{\code{"patterns.freq"}:}{Integer vector. The number of cases
        (rows in the data) per missing data pattern, in the same order as
        the rows of the \code{"patterns"} matrix. If the data is complete,
        a single number (the number of observations). In the case of
        multiple groups: a list of integer vectors.
        Alias: \code{"pattern.freq"}.}
    \item{\code{"patterns.stats"}:}{The sample statistics per missing data
        pattern, as used by the full information maximum likelihood (FIML)
        and two-stage estimation methods. Only available if
        \code{missing = "ml"}, \code{"ml.x"}, \code{"two.stage"} or
        \code{"robust.two.stage"}. A list with one element per missing
        data pattern; each element is a list with the following
        components: \code{cov} (the covariance matrix of the observed
        variables in this pattern, divided by the number of cases),
        \code{mean} (the mean vector of the observed variables in this
        pattern), \code{freq} (the number of cases in this pattern; if
        sampling weights were used, the sum of the weights), and
        \code{var_idx} (a logical vector indicating which observed
        variables are observed in this pattern). In the case of multiple
        groups: a list of such lists. Aliases: \code{"pattern.stats"},
        \code{"missing.stats"}.}
    \item{\code{"response.patterns"}:}{The response patterns of the
        ordered (categorical) variables. A list with the following
        components: \code{pat} (a matrix with one row per observed
        response pattern, sorted from the most to the least frequent
        pattern; the columns correspond to the ordered variables),
        \code{freq} (the number of cases per response pattern),
        \code{npatterns} (the number of observed response patterns),
        \code{total_patterns} (the number of possible response patterns),
        and \code{empty_patterns} (the number of possible response
        patterns that were not observed). Only available if the model
        contains ordered variables and raw data was used. In the case of
        multiple groups: a list of such lists. Alias:
        \code{"resp.patterns"}.}
    \item{\code{"coverage"}:}{A symmetric matrix in which each element
        contains the proportion of jointly observed data points for the
        corresponding pair of observed variables.
        In the case of multiple groups: a list of coverage matrices.}
    \item{\code{"sampling.weights"}:}{Numeric vector. The (normalized)
        sampling weights that were used in the analysis, one per case
        (in the order of \code{"case.idx"}). Depending on the
        \code{sampling_weights_normalization} option, the weights are
        rescaled so that they sum to the number of observations within
        each group (\code{"group"}, the default), to the total number of
        observations (\code{"total"}), or not at all (\code{"none"}).
        Only available if sampling weights were provided.
        In the case of multiple groups: a list of numeric vectors.
        Alias: \code{"weights"}.}
    \item{\code{"group"}:}{A character string. The group variable in
        the data.frame (if any).}
    \item{\code{"ngroups"}:}{Integer. The number of groups.}
    \item{\code{"group.label"}:}{A character vector. The group labels.}
    \item{\code{"level.label"}:}{A character vector. The level labels.}
    \item{\code{"cluster"}:}{A character vector. The cluster variable(s)
        in the data.frame (if any).}
    \item{\code{"nlevels"}:}{Integer. The number of levels.}
    \item{\code{"nclusters"}:}{Integer. The number of clusters
        that were used in the analysis.}
    \item{\code{"ncluster.size"}:}{Integer. The number of different cluster
        sizes.}
    \item{\code{"cluster.size"}:}{Integer vector. The number of observations
        within each cluster. For multigroup multilevel models, a list of
        integer vectors, indicating cluster sizes within each group.}
    \item{\code{"cluster.id"}:}{Integer vector. The cluster IDs identifying
        the clusters. For multigroup multilevel models, a list of
        integer vectors, indicating cluster IDs within each group.}
    \item{\code{"cluster.idx"}:}{Integer vector. The cluster index for
        each observation. The cluster index ranges from 1 to the
        number of clusters. For multigroup multilevel models, a list of
        integer vectors, indicating cluster indices within each group.}
    \item{\code{"cluster.label"}:}{Integer vector. The cluster ID for
        each observation. For multigroup multilevel models, a list of
        integer vectors, indicating the cluster ID for each observation
        within each group.}
    \item{\code{"cluster.sizes"}:}{Integer vector. The different cluster
        sizes that were used in the analysis. For multigroup multilevel
        models, a list of integer vectors, indicating the different cluster
        sizes within each group.}
    \item{\code{"average.cluster.size"}:}{Integer. The average cluster
        size (using the formula
        \code{s = (N^2 - sum(cluster.size^2)) / (N*(nclusters - 1L))}).
        For multigroup multilevel
        models, a list containing the average cluster size per group.}
}

Observed sample statistics:

\describe{
    \item{\code{"sampstat"}:}{Observed sample statistics. Aliases:
         \code{"obs"}, \code{"observed"}, \code{"samp"}, \code{"sample"},
         \code{"samplestatistics"}. Since
         0.6-3, we always check if an h1 slot is available (the estimates
         for the unrestricted model); if present, we extract the sample
         statistics from this slot. This implies that if the variables are
         continuous and \code{missing = "ml"} (or \code{"fiml"}), we
         return the covariance matrix (and mean vector) as computed by
         the EM algorithm under the unrestricted (h1) model. If the h1 slot
         is not present (for example, because the model was fitted with
         \code{h1 = FALSE}), we return the sample statistics from the
         SampleStats slot. In that case, pairwise deletion is used for the
         elements of the covariance matrix (or correlation matrix), and
         listwise deletion for all univariate statistics (means, intercepts,
         and thresholds).}
    \item{\code{"sampstat.std"}:}{Standardized observed sample statistics.
         The covariance matrix is rescaled to a correlation matrix. Aliases:
         \code{"obs.std"}, \code{"observed.std"}, \code{"samp.std"},
         \code{"sample.std"}, \code{"samplestatistics.std"}.}
    \item{\code{"sampstat.h1"}:}{Deprecated. Do not use any longer.}
    \item{\code{"var"}:}{Numeric vector. The observed variances of the
        observed variables: the diagonal of the covariance matrix returned
        by \code{"sampstat"} (that is, the h1-based estimates if the h1
        slot is available). If \code{conditional_x = TRUE}, the variances
        of the joint (y, x) sample covariance matrix are returned; use
        \code{"res.var"} to obtain the residual variances. In the case
        of multiple groups (or levels): a list of vectors.
        Alias: \code{"sampstat.var"}.}
    \item{\code{"res.var"}:}{Numeric vector. The residual variances of the
        (non-exogenous) observed variables, after conditioning on the
        exogenous covariates: the diagonal of the \code{res.cov} matrix
        returned by \code{"sampstat"}. Only available if
        \code{conditional_x = TRUE}.}
    \item{\code{"icov"}:}{The inverse of the observed sample covariance
        matrix (as returned by \code{"sampstat"}). In the case of
        multiple groups: a list of matrices. Not available for multilevel
        models. Alias: \code{"sampstat.icov"}.}
    \item{\code{"res.icov"}:}{The inverse of the residual sample covariance
        matrix (the \code{res.cov} matrix returned by \code{"sampstat"}).
        Only available if \code{conditional_x = TRUE}.}
    \item{\code{"cov.log.det"}:}{Numeric. The logarithm of the determinant
        of the observed sample covariance matrix (as returned by
        \code{"sampstat"}). In the case of multiple groups: a numeric
        vector with one value per group. Not available for multilevel
        models. Alias: \code{"sampstat.cov.log.det"}.}
    \item{\code{"res.cov.log.det"}:}{Numeric. The logarithm of the
        determinant of the residual sample covariance matrix. Only
        available if \code{conditional_x = TRUE}.}
    \item{\code{"res.th.nox"}:}{Numeric vector. For categorical data, the
        sample thresholds of the ordered variables, ignoring the
        exogenous covariates (that is, the thresholds of the marginal
        distribution of each ordered variable). If
        \code{conditional_x = TRUE}, these differ from the (residual)
        thresholds returned by \code{"sampstat"}; if no covariates are
        used, they are identical to the \code{th} element of
        \code{"sampstat"}. In the case of multiple groups: a list of
        vectors. Alias: \code{"th.nox"}.}
    \item{\code{"wls.obs"}:}{The observed sample statistics (covariance
        elements, intercepts/thresholds, etc.) in a single vector.}
    \item{\code{"wls.v"}:}{The weight matrix as used in weighted least
        squares estimation. If the estimator is \code{"DWLS"} or
        \code{"ULS"}, only the diagonal is used, but a full (diagonal)
        matrix is returned.}
    \item{\code{"wls.vd"}:}{Numeric vector. The diagonal of the
        \code{"wls.v"} weight matrix. For the \code{"DWLS"} and
        \code{"ULS"} estimators, this vector contains the (only) weights
        that are used in the estimation. In the case of multiple groups:
        a list of vectors.}
    \item{\code{"gamma"}:}{N times the asymptotic variance matrix of the
        sample statistics. Alias: \code{"sampstat.nacov"}.}
    \item{\code{"group.w"}:}{Numeric vector. The observed group weights:
        the logarithm of the (effective) number of observations in each
        group. This is the same metric as the \code{group.w} element
        of the \code{"implied"} and \code{"sampstat"} output (which is
        only shown if \code{group_w_free = TRUE}). Alias:
        \code{"group.weights"}.}
    \item{\code{"logl"}:}{Numeric vector. The loglikelihood of the fitted
        model, per group (\code{NA} if the estimator does not provide a
        loglikelihood). The sum over the groups is the value returned by
        \code{logLik()}. Alias: \code{"loglik"}.}
    \item{\code{"h1.logl"}:}{Numeric vector. The loglikelihood of the
        unrestricted (h1) model, per group. The sum over the groups is the
        \code{unrestricted.logl} fit measure. Only available if the h1 slot
        is available (\code{h1 = TRUE}, the default). Alias:
        \code{"h1.loglik"}.}
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
    \item{\code{"mm.theta"}:}{The (residual) variance-covariance
        matrix of the observed variables (the theta model matrix).
        Aliases: \code{"theta"}, \code{"theta.cov"}.}
    \item{\code{"theta.cor"}:}{The (residual) variance-covariance
        matrix of the observed variables, expressed in correlation metric.}
    \item{\code{"wls.est"}:}{The model-implied sample statistics (covariance
        elements, intercepts/thresholds, etc.) in a single vector.}
    \item{\code{"vy"}:}{The model-implied unconditional variances of the
        observed variables.}
    \item{\code{"rsquare"}:}{The R-square value for all endogenous variables.
        Aliases: \code{"r-square"}, \code{"r2"}.}
	\item{\code{"fs_determinacy"}:}{The factor determinacies (based on
	   regression factor scores). They represent the (estimated) correlation
       between the factor scores and the latent variable scores.}
	\item{\code{"fs.reliability"}:}{The factor reliabilities (based on
       regression factor scores). They are the square of the factor
       determinacies.}
    \item{\code{"fs_determinacy_Bartlett"}:}{The factor determinacies (based on
       Bartlett factor scores). They represent the (estimated) correlation
       between the factor scores and the latent variable scores.}
    \item{\code{"fs_reliability_bartlett"}:}{The factor reliabilities (based on
       Bartlett factor scores). They are the square of the factor
       determinacies.}
    \item{\code{"icc"}:}{The intraclass correlation coefficients for clustered
       data with multilevel models. Computed from the model-implied (H1)
       within-group and between-group covariance matrices. Only available
       for models with clustered data and multiple levels.}
    \item{\code{"ranef"}:}{The random effects (empirical Bayes predictions of
       the cluster-level latent variables). Only available for clustered data
       (in the long format) with multiple levels.}
    \item{\code{"sampstat.2l"}:}{The intermediate (sufficient) sample
       statistics of the two-level data, as used by the multilevel
       estimation routines. A list with the following components:
       \code{sigma_w} and \code{mu_w} (the pooled within-cluster covariance
       matrix, divided by the number of observations minus the number of
       clusters, and the within-level means), \code{sigma_b} and
       \code{mu_b} (the moment-based between-level covariance matrix and
       means, as used for starting values), \code{s_b} (the between
       covariance matrix of the cluster means), \code{s} (the average
       cluster size, see \code{"average.cluster.size"}), \code{y1y1} (the
       crossproduct matrix of the level-1 data), \code{y2} (the cluster
       means, one row per cluster), \code{mean_d} and \code{cov_d} (the
       mean vectors and covariance matrices of the cluster means, per
       cluster size, in the order of \code{"cluster.sizes"}; the
       between-only variables come first), and \code{loglik_x} (the
       loglikelihood contribution of the fixed exogenous covariates).
       Only available for multilevel models. In the case of multiple
       groups: a list of such lists. Aliases: \code{"sampstat.twolevel"},
       \code{"sampstat.cluster"}.}
}

Diagnostics:

\describe{
    \item{\code{"mdist2.fs"}:}{The squared Mahalanobis distances for the (Bartlett) factor scores.}
    \item{\code{"mdist.fs"}:}{The Mahalanobis distances for the (Bartlett) factor scores.}
    \item{\code{"mdist2.resid"}:}{The squared Mahalanobis distances for the (Bartlett-based) casewise residuals.}
    \item{\code{"mdist.resid"}:}{The Mahalanobis distances for the (Bartlett-based) casewise residuals.}
}

If (some of) the data is ordered categorical, the Mahalanobis distances are
generalized as in Mansolf and Reise (2017): each distance is replaced by its
expected value over the region of the (multivariate normally distributed)
latent response vector that is consistent with the observed response pattern,
approximated by Monte Carlo integration (\code{set.seed()} can be used for
reproducibility). Observed continuous variables are conditioned on, and
missing values are integrated out. Note that in this case
\code{"mdist.fs"}/\code{"mdist.resid"} return the expected distance
\eqn{E[d]} (the person-fit indices \eqn{d_f^*} and \eqn{d_r^*} of Mansolf
and Reise), while \code{"mdist2.fs"}/\code{"mdist2.resid"} return the
expected squared distance \eqn{E[d^2]}, which is not the square of the
former.

Optimizer information:

\describe{
    \item{\code{"converged"}:}{Logical. \code{TRUE} if the optimizer has
        converged; \code{FALSE} otherwise.}
    \item{\code{"iterations"}:}{Integer. The number of iterations used by the
        optimizer.}
    \item{\code{"optim"}:}{List. All available information regarding the
        optimization results.}
    \item{\code{"npar"}:}{Integer. Number of free parameters (ignoring constraints).}
}

Gradient, Hessian, observed, expected and first.order information matrices:

\describe{
    \item{\code{"gradient"}:}{Numeric vector containing the first derivatives
        of the discrepancy function with respect to the (free) model
        parameters.}
    \item{\code{"gradient.logl"}:}{Numeric vector containing the first
        derivatives of the loglikelihood with respect to the (free) model
        parameters.}
    \item{\code{"optim.gradient"}:}{Numeric vector containing the first
        derivatives of the objective function (as seen by the optimizer)
        with respect to the (free) model parameters.}
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
        matrix (depending on the information option of the fitted model).
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
        (depending on the information option of the fitted model).}
    \item{\code{"inverted.information.expected"}:}{Matrix containing the
        inverted expected information matrix for the free model parameters.}
    \item{\code{"inverted.information.observed"}:}{Matrix containing the
        inverted observed information matrix for the free model parameters.}
    \item{\code{"inverted.information.first.order"}:}{Matrix containing the
        inverted first.order information matrix for the free model parameters.}
   \item{\code{"h1.information"}:}{Matrix containing either the observed,
        expected or first.order information matrix (depending on the
        information option of the fitted model) of the unrestricted h1
        model. This is unit-information, not total-information.}
    \item{\code{"h1.information.expected"}:}{Matrix containing the expected
        information matrix for the unrestricted h1 model.}
    \item{\code{"h1.information.observed"}:}{Matrix containing the observed
        information matrix for the unrestricted h1 model.}
    \item{\code{"h1.information.first.order"}:}{Matrix containing the
        first.order information matrix for the unrestricted h1 model.
        Alias: \code{"h1.first.order"}.}
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

Constraint information (equality/inequality constraints, bounds). All
matrices and vectors are reconstructed from the fitted model (the
constraint functions are evaluated at the final parameter estimates),
and are labeled so that rows correspond to individual constraints and
columns to free parameters:

\describe{
    \item{\code{"constraints"}:}{List. All constraint information for the
        fitted model: \code{$ceq.jac} (the Jacobian of the equality
        constraints, one row per constraint, evaluated at the final
        estimates), \code{$ceq.rhs} (right-hand side constants),
        \code{$ceq.resid} (the values of the equality constraint
        functions at the solution; all (near) zero), \code{$ceq.rank}
        (the number of linearly independent equality restrictions; this
        is the number subtracted from the parameter count in the degrees
        of freedom), \code{$ceq.lambda} (Lagrange multipliers, when
        available), \code{$cin.jac}, \code{$cin.rhs}, \code{$cin.slack}
        (the slack of each inequality constraint at the solution;
        (near) zero if the constraint is binding), \code{$cin.active}
        (logical; \code{TRUE} for inequality constraints that were binding
        at the solution), \code{$cin.lambda}, \code{$k} (see
        \code{"constraints.nullspace"} below), \code{$k0} (offset
        vector; only non-\code{NULL} when the model contains only linear
        equality constraints, in which case any admissible parameter
        vector can be written as \code{theta = k0 + K \%*\% theta.reduced}),
        and \code{$ceq.simple} (logical; \code{TRUE} when simple equality
        constraints (\code{a == b}) were absorbed into the parameter
        vector (see the \code{ceq_simple} option), so that the free
        parameter vector already satisfies them and \code{$ceq.jac} is
        empty). The inequality rows appear in this order: finite upper
        bounds, finite lower bounds, explicit \code{"<"}/\code{">"}
        constraints; all are written in slack form (the row value is
        nonnegative at an admissible point). Alias: \code{"con"}.}
    \item{\code{"constraints.jacobian"}:}{Matrix. The Jacobian of all
        constraints: the equality rows stacked above the inequality rows,
        evaluated at the final estimates. The attribute \code{"active"} is
        a logical vector flagging, per row, whether the constraint was
        binding at the solution (equality rows are always active).
        Alias: \code{"con.jac"}.}
    \item{\code{"constraints.nullspace"}:}{Matrix. An orthonormal basis
        \code{K} of the null space of the equality-constraint Jacobian,
        as used by the standard error and (robust) test statistic
        machinery: post-multiplying the Jacobian of the model-implied
        sample statistics (\code{"delta"}) by \code{K} restricts it to
        the directions that satisfy the (linearized) equality
        constraints. For example, the information matrix of the
        constrained model in reduced form is
        \code{t(K) \%*\% I \%*\% K}, and its constrained
        (Moore-Penrose) inverse in the full parameter space is
        \code{K \%*\% solve(t(K) \%*\% I \%*\% K) \%*\% t(K)}. If the model
        contains no equality constraints, the identity matrix is returned.
        The row names indicate the parameter space: they correspond to the
        free parameters, with one row per parameter-table row when simple
        equality constraints were absorbed (\code{ceq_simple = TRUE}).
        Aliases: \code{"constraints.k"}, \code{"eq.constraints.k"},
        \code{"con.k"}.}
}

Miscellaneous:

\describe{
    \item{\code{"coef.boot"}:}{Matrix containing estimated model parameters for each bootstrap sample. Only relevant when bootstrapping was used.}
    \item{\code{"monte.carlo"}:}{Matrix containing the Monte Carlo draws of the
        model parameters, as used for the Monte Carlo confidence intervals of
        defined parameters (Preacher & Selig, 2012). Only relevant when these
        draws were requested. Aliases: \code{"mc"}, \code{"mc.coef"},
        \code{"coef.mc"}.}
    \item{\code{"UGamma"}:}{Matrix containing the product of 'U' and 'Gamma'
        matrices as used by the Satorra-Bentler correction. The trace of
        this matrix, divided by the degrees of freedom, gives the scaling
        factor.}
    \item{\code{"UfromUGamma"}:}{Matrix containing the 'U' matrix
        as used by the Satorra-Bentler correction. Alias: \code{"U"}.}
    \item{\code{"delta"}:}{The delta matrix (per group): the Jacobian of the
        model-implied summary statistics (the rows: means, (co)variances,
        thresholds, \ldots) with respect to the free model parameters (the
        columns). This is the matrix used in the delta method. (Note: this is
        the delta \emph{Jacobian}, not the delta model matrix that is returned
        by \code{"mm.delta"}.) In the case of multiple groups: a list of
        matrices.}
    \item{\code{"delta.rownames"}:}{The row names of the delta matrix (see
        \code{"delta"}): the labels of the model-implied summary statistics.
        In the case of multiple groups: a list of character vectors.}
    \item{\code{"list"}:}{The parameter table. The same output as given
        by \code{parTable()}.}
    \item{\code{"fit"}:}{The fit measures. Aliases: \code{"fitmeasures"},
        \code{"fit.measures"}, \code{"fit.indices"}. The same output as
        given by \code{fitMeasures()}.}
    \item{\code{"mi"}:}{The modification indices. Alias: \code{"modindices"},
        \code{"modification.indices"}. The same output as given
        by \code{modindices()}.}
    \item{\code{"loglik.casewise"}:}{Vector containing the casewise
        loglikelihood contributions. Only available if estimator = \code{"ML"}.}
    \item{\code{"options"}:}{List. The option list.}
    \item{\code{"call"}:}{List. The call as returned by match.call, coerced to
        a list.}
    \item{\code{"timing"}:}{List. The timing (in milliseconds) of various
        lavaan subprocedures.}
    \item{\code{"test"}:}{List. All available information regarding the
        (goodness-of-fit) test statistic(s).}
    \item{\code{"baseline.test"}:}{List. All available information regarding the
        (goodness-of-fit) test statistic(s) of the baseline model.}
    \item{\code{"baseline.partable"}:}{Data.frame. The parameter table of
        the (internal) baseline model.}
    \item{\code{"post.check"}:}{Post-fitting check if the solution is
        admissible. A warning is raised if negative variances are found, or if
        either \code{lavInspect(fit, "cov.lv")} or
        \code{lavInspect(fit, "theta")} return a non-positive definite matrix.}
    \item{\code{"zero.cell.tables"}:}{List. List of bivariate frequency tables
        where at least one cell is empty.}
    \item{\code{"iv"}:}{Data.frame. One row per (transformed) equation, with
        the (model-implied) instrumental variables that were used for each
        equation. Columns are \code{lhs}, \code{rhs}, \code{lhs.new},
        \code{rhs.new}, \code{type}, and \code{instruments}. Only available
        when \code{estimator = "IV"}. Aliases: \code{"ivs"}, \code{"miiv"},
        \code{"miivs"}, \code{"instr"}, \code{"instruments"}. In the case of
        multiple groups: a list of data.frames.}
    \item{\code{"eqs"}:}{List. The (internal) list of equations that was used
        for the (MI)IV-based estimation, including the instruments, the
        coefficients, and the Sargan test per equation. Only available when
        \code{estimator = "IV"}.}
    \item{\code{"sargan"}:}{Data.frame. The per-equation overidentification
        tests, for each overidentified equation. Columns are
        \code{lhs}, \code{rhs}, \code{df}, \code{sargan.stat},
        \code{sargan.pval}, \code{browne.stat}, and \code{browne.pval}. The \code{sargan.*} columns are the classical
        Sargan test (assumes normality; the p-value is \code{NA} for
        categorical data, where it is not valid); the \code{browne.*} columns
        are a robust residual-based test computed with the distribution-free
        (ADF) ACOV for continuous data or the polychoric ACOV for categorical
        data. \code{"hansen"} is an alias. Only available when
        \code{estimator = "IV"}. In the case of multiple groups: a list of
        data.frames.}
    \item{\code{"sam.struc.fit.object"}:}{The (internal) step-2 structural
        model of a model fitted with \code{\link{sam}} using a local
        \code{sam.method} (\code{"local"}, \code{"fsr"} or \code{"cfsr"}):
        a lavaan object, fitted to the estimated latent (co)variances (and
        means) of step 1. \code{fitMeasures()}, \code{lavResiduals()} and
        \code{modindices()} operate on this object when they are called with
        the sam object itself. \code{"sam.struc"} is an alias.}
    \item{\code{"version"}:}{The lavaan version number that was used to
        construct the fitted lavaan object.}
}

}
\seealso{
\code{\link{lavaan}}, \code{\link{estimator_iv}}, \code{\link{estimator_js}}
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
