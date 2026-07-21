\name{lavOptions}
\alias{lavOptions}
\alias{lavoptions}
\title{lavaan Options}
\description{
Show the default options used by the \code{lavaan()} function. The
options can be changed by passing 'name = value' arguments to the
\code{lavaan()} function call, where they are added to the '...'
argument.
}
\usage{
lavOptions(x = NULL, default = NULL, mimic = "lavaan")
}
\arguments{
\item{x}{Character. A character string holding an option name, or a character
string vector holding multiple option names. All option names are converted to
lower case.}
\item{default}{If a single option is specified but not available, this value
is returned.}
\item{mimic}{Not used for now.}
}
\details{
This is the full list of options that are accepted by the \code{lavaan()}
function, organized in several sections:


Model features:

\describe{
    \item{\code{meanstructure}:}{If \code{TRUE}, the means of the observed
      variables enter the model. If \code{"default"}, the value is set based
      on the user-specified model, and/or the values of other arguments.}
    \item{\code{int.ov.free}:}{If \code{FALSE}, the intercepts of the 
      observed variables are fixed to zero.}
    \item{\code{int.lv.free}:}{If \code{FALSE}, the intercepts of the latent 
      variables are fixed to zero.}
    \item{\code{conditional.x}:}{If \code{TRUE}, we set up the model 
      conditional on the exogenous `x' covariates; the model-implied sample 
      statistics only include the non-x variables. If \code{FALSE}, the 
      exogenous `x' variables are modeled jointly with the other variables, and 
      the model-implied statistics reflect both sets of variables. If
      \code{"default"}, the value is set depending on the estimator, and
      whether or not the model involves categorical endogenous variables.}
    \item{\code{fixed.x}:}{If \code{TRUE}, the exogenous `x' covariates are 
      considered fixed variables and the means, variances and covariances of 
      these variables are fixed to their sample values. If \code{FALSE}, they 
      are considered random, and the means, variances and covariances are free 
      parameters. If \code{"default"}, the value is set depending on the mimic 
      option.}
    \item{\code{orthogonal}:}{If \code{TRUE}, all covariances among
      latent variables are set to zero.}
    \item{\code{orthogonal.y}:}{If \code{TRUE}, all covariances among
      endogenous latent variables only are set to zero.}
    \item{\code{orthogonal.x}:}{If \code{TRUE}, all covariances among
      exogenous latent variables only are set to zero.}
    \item{\code{std.lv}:}{If \code{TRUE}, the metric of each latent variable 
      is determined by fixing their (residual) variances to 1.0. If 
      \code{FALSE}, the metric of each latent variable is determined by fixing 
      the factor loading of the first indicator to 1.0. If there are multiple
      groups, \code{std.lv = TRUE} and \code{"loadings"} is included in 
      the \code{group.equal} argument, then only the latent variances
      of the first group will be fixed to 1.0, while the latent
      variances of other groups are set free.}
    \item{\code{effect.coding}:}{Can be logical or character string. If
     logical and \code{TRUE}, this implies 
     \code{effect.coding = c("loadings", "intercepts")}. If logical and 
     \code{FALSE}, it is set equal to the empty string. 
     If \code{"loadings"} is included, equality
     constraints are used so that the average of the factor loadings (per
     latent variable) equals 1. Note that this should not be used 
     together with \code{std.lv = TRUE}. If \code{"intercepts"} is
     included, equality constraints are used so that the sum of the
     intercepts (belonging to the indicators of a single latent variable) 
     equals zero.
     As a result, the latent mean is freely estimated and usually
     equals the average of the means of the involved indicators.}
    \item{\code{ceq.simple}:}{Logical. If \code{TRUE}, and no other 
     general (equality or inequality) constraints are used in the model, 
     simple equality constraints
     are represented in the parameter table as duplicated free parameters
     (instead of extra rows with \code{op = "=="}). The default is
     \code{FALSE}.}
    \item{\code{parameterization}:}{Currently only used if data is 
      categorical. If \code{"delta"}, the delta parameterization is used. 
      If \code{"theta"}, the theta parameterization is used.}
    \item{\code{correlation}:}{Only used for (single-level)
     continuous data. If \code{TRUE}, analyze a correlation matrix (instead
     of a (co)variance matrix). This implies that the residual observed
     variances are no longer free parameters. Instead, they are set to
     values that ensure the model-implied variances are unity. This also
     affects the standard errors. Since 0.7-2, the default estimator is
     ML (before: GLS): the model is fit as the covariance structure
     \eqn{\Sigma = \Delta P(\theta) \Delta} with free scale parameters
     \eqn{\Delta} (the \code{~*~} rows in the output) and unit-diagonal
     \eqn{P(\theta)}, by ordinary maximum likelihood in the full moment
     space. Because this augmented model is exactly scale-equivariant,
     the correlation-metric estimates, standard errors and test
     statistic are identical whether the input is raw data, a
     covariance matrix, or a correlation matrix. The other ML-family
     estimators (\code{"MLM"}, \code{"MLR"}, \ldots) work the same way
     (raw data needed). A meanstructure is supported: the (saturated)
     intercepts live in the original metric and do not affect the
     correlation structure. The ML route requires single-level,
     complete data, \code{conditional.x = FALSE} and no composites; in
     all other settings, a defaulted estimator quietly falls back to
     GLS. Estimator WLS (= ADF weights in the correlation metric)
     covers the non-normal case.
     Both
     \code{fixed.x = FALSE} and \code{fixed.x = TRUE} are supported (the
     latter affects the way the standard errors are computed), and (since
     0.7-2) also \code{conditional.x = TRUE}: the exogenous covariates are
     then standardized along with the outcomes, the input
     \code{sample.cov} (if used) must be the joint (correlation or
     covariance) matrix including the x-variables (it is standardized
     internally), and the unit-variance constraint applies to the
     marginal (not the residual) variances of the outcomes.
     Alternatively, \code{correlation} may be a character vector of
     (observed) variable names (for example
     \code{correlation = c("x1", "x2", "x3")}). In that case only the
     listed variables are standardized to unit variance (their residual
     variances are fixed accordingly), while all remaining observed
     variables keep their original metric (a 'partial' correlation
     structure). This also works with \code{conditional.x = TRUE}: the
     listed exogenous variables are standardized along with the listed
     outcomes, and non-listed outcomes keep a free (residual) variance. If the model contains observed product (interaction)
     terms (for example \code{x1:y1}) and \code{correlation = TRUE},
     the product terms are automatically excluded from the correlation
     set: the product of standardized variables does not have unit
     variance, so their (co)variances remain free parameters. Note that
     product terms are typically highly non-normal, so the normal-theory
     (GLS) standard errors for their (co)variances may be too small;
     consider \code{se = "robust.sem"}.}
    \item{\code{composites.cov}:}{Character string. Only used if the model
     contains composites. Controls whether the composite-indicator
     (co)variances (the elements of the \code{T} matrix) are fixed to their
     sample values or estimated as free parameters. If \code{"fixed"}, the
     \code{T} matrix is fixed to the sample (co)variances. If \code{"free"},
     the \code{T} matrix is estimated as free parameters. If
     \code{"default"}, the value is resolved to \code{"free"} for multilevel
     models (which makes composites at every level identified) and
     \code{"fixed"} otherwise. For a single-level model, freeing \code{T}
     is fit-invariant (the degrees of freedom are unchanged).}
    \item{\code{auto.fix.first}:}{If \code{TRUE}, the factor loading of the
      first indicator is set to 1.0 for every latent variable.}
    \item{\code{bad.marker.crit}:}{Only used if
      \code{auto.fix.first = TRUE}. A single number (default 0.1). If larger
      than zero, lavaan checks -- once the unrestricted (h1) sample statistics
      are available -- whether the first indicator of each latent variable is
      a poor item, in the sense that it correlates weakly with the other
      indicators of the same factor (more precisely: the absolute value of its
      corrected item-total correlation is below \code{bad.marker.crit}). If so
      (and a clearly better indicator is available), a warning is issued and
      the loading of that better indicator is fixed to 1.0 instead, in order
      to set the metric of that latent variable. The main purpose is to avoid
      convergence problems caused by a (very) poor marker item. The check is
      based on the (pooled) unrestricted (h1) covariance matrix, so it behaves
      consistently with or without missing data, categorical data, multiple
      groups, etc. Larger values make the switching more aggressive; a value
      of 0 disables the check and always keeps the first indicator as the
      marker.}
    \item{\code{auto.fix.single}:}{If \code{TRUE}, the residual variance (if
      included) of an observed indicator is set to zero if it is the only 
      indicator of a latent variable.}
    \item{auto.var}{If \code{TRUE}, the (residual) variances of both observed
      and latent variables are set free.}
    \item{\code{auto.cov.lv.x}:}{If \code{TRUE}, the covariances of exogenous
      latent variables are included in the model and set free.}
    \item{\code{auto.cov.x}:}{If \code{TRUE}, the covariances between
      exogenous latent variables and observed exogenous covariates are also
      included in the model and set free (this implies
      \code{auto.cov.lv.x = TRUE}). Ignored if \code{conditional.x = TRUE},
      as these covariances cannot be represented when the model is
      conditioned on the exogenous covariates.}
    \item{\code{auto.cov.y}:}{If \code{TRUE}, the covariances of dependent 
      variables (both observed and latent) are included in the model and set 
      free.}
    \item{\code{auto.th}:}{If \code{TRUE}, thresholds for limited dependent 
      variables are included in the model and set free.}
    \item{\code{auto.delta}:}{If \code{TRUE}, response scaling parameters 
      for limited dependent variables are included in the model and set free.}
    \item{\code{auto.efa}:}{If \code{TRUE}, the necessary constraints are
      imposed to make the (unrotated) exploratory factor analysis blocks
      identifiable: for each block, factor variances are set to 1, factor
      covariances are constrained to be zero, and factor loadings are
      constrained to follow an echelon pattern.}
}

Data options:

\describe{
    \item{\code{std.ov}:}{If \code{TRUE}, observed variables are 
      standardized before entering the analysis. By default, these are
      only the non-exogenous observed variables, unless \code{fixed.x = FALSE}.
      Use this option with caution; it can be used to test whether (for
      example) nonconvergence was due to scaling issues. Note that this is
      still a covariance-based analysis: no constraints are imposed to
      ensure the model-implied (co)variance matrix has unit variances, and
      the standard errors still assume that the input was unstandardized. See
      also the \code{correlation} option.}
    \item{\code{missing}:}{The default setting is \code{"listwise"}: all
         cases with missing values
      are removed listwise from the data before the analysis starts. This is
      only valid if the data are missing completely at random (MCAR).
      Therefore, it may not be the optimal choice, but
      it can be useful for a first run. If the estimator belongs to
      the ML family, another option is \code{"ml"} (alias: \code{"fiml"}
      or \code{"direct"}). This corresponds to the so-called full information
      maximum likelihood approach (fiml), where we compute the likelihood
      case by case, using all available data from that case. Note
      that if the model contains exogenous observed covariates, and
      \code{fixed.x = TRUE} (the default), all cases with any missing values
      on these covariates will be deleted first. The option \code{"ml.x"}
      (alias: \code{"fiml.x"} or \code{"direct.x"}) is similar to \code{"ml"},
      but does not delete any cases with missing values for the exogenous
      covariates, even if \code{fixed.x = TRUE}. (Note: all lavaan versions
      < 0.6 used \code{"ml.x"} instead of \code{"ml"}).
      Since version 0.7-2, \code{"ml"} can also be combined with
      \code{conditional.x = TRUE} (single-level, continuous data): the
      likelihood of the observed non-x variables is then computed case by
      case, conditional on the (complete) exogenous covariates of that
      case. Cases with missing values on the exogenous covariates are
      deleted first, and \code{"ml.x"} is not available in this setting.
      Note that the robust (missing-data corrected) versions of RMSEA/CFI
      are not (yet) available when \code{conditional.x = TRUE}.
      If you wish to use multiple
      imputation, you need to use an external package (e.g., mice) to
      generate imputed datasets, which can then be analyzed using
      the \code{\link{semList}} function. The semTools package contains
      several functions to do this automatically. Another option (with
      continuous data) is to use \code{"two.stage"}
      or \code{"robust.two.stage"}. In this approach, we first estimate
      the sample statistics (mean vector, variance-covariance matrix) using
      an EM algorithm. Then, we use these estimated sample statistics as
      input for a regular analysis (as if the data were complete). The
      standard errors and test statistics
      are adjusted correctly to reflect the two-step procedure. The
      \code{"robust.two.stage"} option produces standard errors and
      a test statistic that are robust against non-normality.
      Both options are available for the ML estimator (since version
      0.7-2 also in combination with \code{conditional.x = TRUE}), and
      also for the
      (continuous-data) least-squares estimators \code{"ULS"}, \code{"GLS"},
      \code{"WLS"} and \code{"DLS"}; for the latter, lavaan automatically
      switches to \code{se = "robust.sem"} and \code{test = "satorra.bentler"},
      using the two-stage variance-covariance matrix of the EM-based sample
      statistics. In addition, for these least-squares estimators, requesting
      \code{missing = "ml"} (or \code{"fiml"}) is interpreted as a generic
      request to handle the missing values, and is silently treated as
      \code{missing = "two.stage"}.
      If (part of) the data is categorical, and the estimator is
      from the (W)LS family, the only option (besides listwise deletion)
      is \code{"pairwise"}. In this three-step approach, missingness is
      only an issue in the first two steps. In the first step, we compute
      thresholds (for categorical variables) and means or intercepts
      (for continuous variables) using univariate information only.
      In this step, we simply ignore
	  the missing values just like in mean(x, na.rm = TRUE). In the second
      step, we compute polychoric/polyserial/pearson correlations using (only)
      two variables at a time. Here we use pairwise deletion: we only keep
      those observations for which both values are observed (not missing),
      and this set may change from pair to pair.
      By default, in the categorical case we use \code{conditional.x = TRUE}. 
      Therefore, any cases
      with missing values on the exogenous covariates will be deleted listwise
      from the data first.
      Finally, if the estimator is \code{"PML"}, the available options are
      \code{"pairwise"}, \code{"available.cases"} and
      \code{"doubly.robust"}. See the PML tutorial on the lavaan website for
      more information about these approaches.}
    \item{\code{sampling.weights.normalization}:}{If \code{"none"}, the
      sampling weights (if provided) will not be transformed. If \code{"total"},
      the sampling weights are normalized by dividing by the total sum of
      the weights, and multiplying again by the total sample size. 
      If \code{"group"}, the sampling weights are normalized per group:
      by dividing by the sum of the weights (in each group), and multiplying
      again by the group size. The default is \code{"group"}. (For a single
      group, \code{"group"} and \code{"total"} are identical; they only differ
      when there are multiple groups, in which case \code{"group"} keeps each
      group's relative contribution proportional to its sample size, matching
      the behavior of Mplus.)}
    \item{\code{sampling.weights.type}:}{Only used when sampling weights are
      provided and the estimator relies on a weight matrix or sandwich variance:
      the least-squares family (GLS/WLS/DWLS/ULS/DLS, both continuous and
      categorical) and PML. Determines how the sampling weights enter the
      asymptotic variance (the Gamma / NACOV matrix, or the PML first-order
      information), and hence the (robust) standard errors and the scaled test
      statistic. For GLS and ULS only the standard errors are affected; for
      WLS/DWLS/DLS the Gamma also feeds the weight matrix, so the point
      estimates shift as well. If
      \code{"design"} (the default), the weights are treated as sampling
      (design) weights and a design-based sandwich is used (the meat is
      weighted by the sum of the squared weights), matching the behavior of
      Mplus. If \code{"frequency"}, the weights are treated as frequency
      (replication) counts (the meat is weighted by the sum of the weights),
      so that the results match a fit on the row-replicated data. Without
      sampling weights the two are identical.}
    \item{\code{samplestats}:}{Logical. If \code{FALSE}, no sample statistics
      will be computed (and no estimation can take place). This can be useful
      when only a dummy lavaan object is requested, without any computations.
      The default is \code{TRUE}.}
}

Data summary options:

\describe{
    \item{\code{sample.cov.rescale}:}{If \code{TRUE}, the sample covariance 
      matrix provided by the user is internally rescaled by multiplying it 
      with a factor (N-1)/N.  If \code{"default"}, the value is set depending 
      on the estimator and the likelihood option: it is set to \code{TRUE} if 
      maximum likelihood estimation is used and \code{likelihood="normal"}, 
      and \code{FALSE} otherwise.}
    \item{\code{ridge}:}{Logical. If \code{TRUE}, a small constant value is
      added to the diagonal elements of the covariance (or correlation)
      matrix before analysis. The value can be set using the
      \code{ridge.constant} option.}
    \item{\code{ridge.constant}:}{Numeric. Small constant used for ridging.
      The default value is 1e-05.}
}

Multiple group options:

\describe{ 
    \item{\code{group.label}:}{A character vector. The user can specify 
      which group (or factor) levels need to be selected from the grouping 
      variable, and in which order. If missing, all grouping levels are 
      selected, in the order as they appear in the data.}
    \item{\code{group.equal}:}{A vector of character strings. Only used in
      a multiple group analysis. Can be one or more of the following:
      \code{"loadings"}, \code{"composite.weights"}, 
      \code{"intercepts"}, \code{"means"},
      \code{"thresholds"}, \code{"regressions"}, \code{"residuals"},
      \code{"residual.covariances"}, \code{"lv.variances"} or
      \code{"lv.covariances"}, specifying the pattern of equality
      constraints across multiple groups. As a shortcut, the single
      value \code{"all"} can be used to constrain all (free) parameters
      to be equal across groups; it is equivalent to listing all of the
      values above. When the model is not the same across groups (for
      example when the \code{group:} modifier is used to specify a
      different model per group), only parameters that are present in two
      or more groups are constrained to be equal; a parameter that occurs
      in a single group only is always left free. The constraints are
      symmetric and do not depend on the first group: a parameter that is
      shared by (say) the second and third group is constrained to be
      equal across those two groups, even if it does not appear in the
      first group.}
    \item{\code{group.partial}:}{A vector of character strings containing 
      the labels of the parameters which should be free in all groups (thereby
      overriding the group.equal argument for some specific parameters).}
    \item{\code{group.w.free}:}{Logical. If \code{TRUE}, the group 
      frequencies are considered to be free parameters in the model. In this 
      case, a Poisson model is fitted to estimate the group frequencies. If
      \code{FALSE} (the default), the group frequencies are fixed to their
      observed values.}
}

Estimation options:

\describe{
    \item{\code{estimator}:}{The estimator to be used. Can be one of the 
      following: \code{"ML"} for maximum likelihood, \code{"GLS"} for 
      (normal theory) generalized least squares, 
      \code{"WLS"} for weighted least squares 
      (sometimes called ADF estimation), \code{"ULS"} for unweighted least 
      squares, \code{"DWLS"} for diagonally weighted least squares,
      and \code{"DLS"} for distributionally-weighted least squares. These
      are the main options that affect the estimation. For convenience, the
      \code{"ML"} option can be extended as \code{"MLM"}, \code{"MLMV"},
      \code{"MLMVS"}, \code{"MLF"}, and \code{"MLR"}.
      The estimation will still be plain \code{"ML"}, but now
      with robust standard errors and a robust (scaled) test statistic. For
      \code{"MLM"}, \code{"MLMV"}, \code{"MLMVS"}, classic robust standard
      errors are used (\code{se="robust.sem"}); for \code{"MLF"}, standard
      errors are based on first-order derivatives 
      (\code{information = "first.order"});
      for \code{"MLR"}, `Huber-White' robust standard errors are used
      (\code{se="robust.huber.white"}). In addition, \code{"MLM"} will compute
      a Satorra-Bentler scaled (mean adjusted) test statistic
      (\code{test="satorra.bentler"}), \code{"MLMVS"} will compute a
      mean and variance adjusted test statistic (Satterthwaite style)
      (\code{test="mean.var.adjusted"}), \code{"MLMV"} will compute a mean
      and variance adjusted test statistic (scaled and shifted)
      (\code{test="scaled.shifted"}), and \code{"MLR"} will
      compute a test statistic which is asymptotically
      equivalent to the Yuan-Bentler T2-star test statistic
      (\code{test="yuan.bentler.mplus"}). Analogously,
      the estimators \code{"WLSM"} and \code{"WLSMV"} imply the \code{"DWLS"}
      estimator (not the \code{"WLS"} estimator) with robust standard errors
      and a mean or mean and variance adjusted test statistic. Estimators
      \code{"ULSM"} and \code{"ULSMV"} imply the \code{"ULS"}
      estimator with robust standard errors
      and a mean or mean and variance adjusted test statistic.
      Finally, \code{"RBM"} requests reduced-bias M-estimation
      (penalized maximum likelihood; continuous-outcome models for now). By
      default sandwich standard errors are used
      (\code{se = "robust.huber.white"}). The type of bias reduction is
      controlled by the \code{rbm.method} element of \code{estimator.args}
      (either \code{"implicit"} (default) or \code{"explicit"}), for example
      \code{estimator = list(estimator = "RBM", rbm.method = "explicit")}.
      The \code{"explicit"} method is much faster for models with many
      parameters (it requires a single ML fit plus a one-step correction),
      and is recommended in that case.}
    \item{\code{likelihood}:}{Only relevant for ML estimation. If 
      \code{"wishart"}, the Wishart likelihood approach is used. In this 
      approach, the covariance matrix has been divided by N-1, and both 
      standard errors and test statistics are based on N-1.
      If \code{"normal"}, the normal likelihood approach is used. Here,
      the covariance matrix has been divided by N, and both standard errors
      and test statistics are based on N. If \code{"default"}, it depends
      on the mimic option: if \code{mimic="lavaan"} or \code{mimic="Mplus"},
      normal likelihood is used; otherwise, Wishart likelihood is used.}
    \item{\code{link}:}{Not used yet. This is just a placeholder until
      the MML estimator is back.}
    \item{\code{information}:}{If \code{"expected"}, the expected 
      information matrix is used (to compute the standard errors). If 
      \code{"observed"}, the observed information matrix is used. 
      If \code{"first.order"}, the information matrix is based on the
      outer product of the casewise scores. See also the options
      \code{"h1.information"} and \code{"observed.information"} for
      further control. If \code{"default"}, the value is set depending 
      on the estimator, the missing argument, and the mimic option. If
      the argument is a vector with two elements, the first element
      is used for the computation of the standard errors, while the
      second element is used for the (robust) test statistic.}
   \item{\code{h1.missing.method}:}{Character. Only used for multilevel
      models with missing data (and \code{missing = "ml"}). If \code{"em"}
      (the default since 0.7-1), the means and (co)variances of the
      unrestricted (h1) model are estimated using the EM algorithm,
      integrating over both the cluster random effects and the missing
      values (this mimics Mplus). If \code{"fiml"}, the unrestricted model
      is estimated using quasi-Newton optimization (nlminb) with a cholesky
      parameterization of the covariance matrices (this was the default
      in versions < 0.7-1). See the \code{em.h1.args} option to control
      the EM algorithm.}
   \item{\code{em.h1.args}:}{List. Options related to the EM algorithm
      that is used to estimate the means and (co)variances of the
      unrestricted (h1) model. This is only relevant for multilevel models,
      and for single-level models with missing data (and
      \code{missing = "ml"} or \code{missing = "two.stage"}). The default
      values are \code{em.h1.args = list(max_iter = "default",
      tol = "default", min_variance = 1e-05, warn = TRUE,
      non_pd_action = "stop", non_pd_tol = 1e-05,
      acceleration = "squarem", fused = TRUE)}. The \code{max_iter}
      element sets the
      maximum number of EM iterations; the \code{tol} element sets the
      tolerance used to determine convergence. If \code{"default"}, they
      are set to 500 and 1e-05 respectively in the single-level setting
      (where the EM algorithm stops when the largest absolute change in
      the parameter values is smaller than \code{tol}), and to 5000 and
      1e-04 respectively in the multilevel setting (where the EM algorithm
      stops when the change in the loglikelihood is smaller than
      \code{tol}). The \code{min_variance} element sets the minimum value
      a variance can take during the EM iterations (multilevel only):
      whenever a within-level variance drops below this value, it is fixed
      to this value (and the corresponding covariances are set to zero).
      If the \code{warn} element is \code{TRUE} (the default),
      a warning is produced if the EM algorithm reaches \code{max_iter}
      iterations without converging. The \code{non_pd_action} element
      controls what happens if an EM estimated variance-covariance matrix
      is (near) singular: if any variance is (near) zero, or if the
      smallest eigenvalue of the matrix is smaller than the
      \code{non_pd_tol} element. This typically happens when a (within or
      between) variance of an observed variable approaches zero, but the
      full matrix may also be (near) singular. If \code{non_pd_action =
      "stop"} (the default), estimation stops with an error message; if
      \code{"warn"}, a warning is produced, and estimation continues; if
      \code{"none"}, estimation continues silently. The
      \code{acceleration} element controls whether the EM iterations are
      accelerated. If \code{"squarem"} (the default), the SQUAREM scheme
      (Varadhan & Roland, 2008) is used: every cycle takes two EM steps,
      extrapolates, and takes one stabilizing EM step from the
      extrapolated point, with a loglikelihood safeguard; this typically
      reduces the number of EM steps considerably (in particular when
      some variances approach zero), and converges to the same solution.
      If \code{"qn"}, a quasi-Newton scheme (Zhou, Alexander & Lange,
      2011; the \code{"qn"} method of the \pkg{turboEM} package) is used
      instead: a low-rank secant approximation of the EM map is used to
      take a quasi-Newton step towards the fixed point, again with a
      loglikelihood safeguard. If \code{"none"}, the plain (unaccelerated)
      EM iterations are used. If the \code{fused} element is \code{TRUE}
      (the default), the
      loglikelihood is computed as a cheap byproduct of the E-step
      quantities (only used in the two-level setting with missing data,
      for the plain EM iterations); if \code{FALSE}, a separate
      (stand-alone) loglikelihood evaluation is used in every iteration;
      both give identical results.}
   \item{\code{h1.information}:}{If \code{"structured"} (the default), the
      unrestricted (h1) information part of the (expected, first.order or 
      observed if h1 is used) information matrix is based on the structured, 
      or model-implied statistics (model-implied covariance matrix, 
      model-implied mean vector, etc.). 
      If \code{"unstructured"}, the unrestricted (h1) information part is
      based on sample-based statistics (observed covariance matrix, observed
      mean vector, etc.) If
      the argument is a vector with two elements, the first element
      is used for the computation of the standard errors, while the
      second element is used for the (robust) test statistic.}
   \item{\code{observed.information}:}{If \code{"hessian"}, the observed
      information matrix is based on the hessian of the objective function.
      If \code{"h1"}, an approximation is used that is based on 
      the observed information matrix of the unrestricted (h1) model. If
      the argument is a vector with two elements, the first element
      is used for the computation of the standard errors, while the
      second element is used for the (robust) test statistic.}
   \item{\code{se}:}{If \code{"standard"}, conventional standard errors
      are computed based on inverting the (expected, observed or first.order) 
      information matrix. If \code{"robust.sem"}, conventional robust
      standard errors are computed.  If \code{"robust.huber.white"},
      standard errors are computed based on the 'mlr' (aka pseudo ML,
      Huber-White) approach.
      If \code{"robust"}, either \code{"robust.sem"} or
      \code{"robust.huber.white"} is used depending on the estimator,
      the mimic option, and whether the data are complete or not.
      For the non-iterative estimators with delta-method standard errors
      over the sample moments (\code{"MGM"}, \code{"JS"}/\code{"JSA"}),
      \code{"robust"} selects the ADF (distribution-free) moment
      covariance, which equals the infinitesimal-jackknife (sandwich)
      covariance of the estimator.
      If \code{"boot"} or \code{"bootstrap"}, bootstrap standard errors are
      computed using standard bootstrapping (unless Bollen-Stine bootstrapping
      is requested for the test statistic; in this case bootstrap standard
      errors are computed using model-based bootstrapping).
      If \code{"none"}, no standard errors are computed.}
   \item{\code{test}:}{Character vector. See the documentation of
      the \code{\link{lavTest}} function for a full list. Multiple names
      of test statistics can be provided. If \code{"default"}, the value 
      depends on the values of other arguments. See also the 
      \code{\link{lavTest}} function to extract (alternative) 
      test statistics from a fitted lavaan object. FMG test names such as
      \code{"peba4"}, \code{"pols3"}, \code{"pall"}, and \code{"all"} can be
      requested here for supported ML-family models. The convenience value
      \code{"fmg"} resolves to \code{"peba4_rls"} for one-model tests.}
   \item{\code{standard.test}:}{Character. Choose the test statistic
    that will be used to compute fit measures (like CFI or RMSEA).
    The default is \code{"standard"}, but it could also be (for example)
    \code{"Browne.residual.nt"}, or an FMG test name such as
    \code{"peba4"} or \code{"fmg"}.}
   \item{\code{scaled.test}:}{Character. Choose the test statistic
    that will be scaled (if a scaled test statistic is requested).
    The default is \code{"standard"}, but it could also be (for example)
    \code{"Browne.residual.nt"}.}
   \item{\code{gamma.n.minus.one}}{Logical. If \code{TRUE}, we divide the
     Gamma matrix by N-1 (instead of the default N).}
   \item{\code{gamma.unbiased}}{Logical. If \code{TRUE}, we compute an
     unbiased version for the Gamma matrix. Only available for single-level
     complete data and when \code{conditional.x = FALSE} and
     \code{fixed.x = FALSE} (for now). Suffixless FMG test names use this
     option to choose between biased and unbiased Gamma.}
   \item{\code{se.delta.second.order}}{Logical. If \code{FALSE} (the
     default), the standard errors of defined parameters (introduced with
     the \code{:=} operator) and of standardized parameters (as returned
     by \code{standardizedSolution()} or \code{lavInspect(., "vcov.std.*")})
     are computed using the first-order delta method. If \code{TRUE}, the
     second-order delta method is used, which adds the term
     \eqn{\frac{1}{2}\textrm{tr}(H V H V)} (where \eqn{H} is the Hessian
     of the transformed parameter and \eqn{V} is the parameter covariance
     matrix) to the first-order variance approximation. This may improve
     accuracy when the transformation is a strongly non-linear function
     of the model parameters. Has no effect when standard errors are
     obtained by bootstrapping or by the Monte Carlo method.}
   \item{\code{se.def}}{Character. Method for computing standard errors
     (and confidence intervals) of defined parameters (introduced with
     the \code{:=} operator). The default is \code{"default"} (or
     equivalently \code{"delta"}), which uses the delta method (respecting
     \code{se.delta.second.order}). Set to \code{"monte.carlo"} to use
     the Monte Carlo method of Preacher and Selig (2012). Under this
     method, \code{R} draws are taken from
     \eqn{\textit{MVN}(\hat{\theta}, \widehat{\textit{Var}}(\hat{\theta}))},
     the defined parameter is evaluated for each draw, and standard
     errors are obtained as the standard deviation of these realizations.
     Confidence intervals (in \code{parameterEstimates()} and
     \code{standardizedSolution()}) are based on percentile bounds of
     the Monte Carlo distribution, and therefore can be asymmetric (just
     like bootstrap intervals). Has no effect when
     \code{se = "bootstrap"}.}
   \item{\code{monte.carlo}}{List of settings for the Monte Carlo method
     (used when \code{se.def = "monte.carlo"}). Currently recognizes:
     \code{R} (integer; number of Monte Carlo draws; default 20000) and
     \code{seed} (integer or \code{NULL}; optional seed for
     reproducibility).}
   \item{\code{bootstrap}:}{Number of bootstrap draws, if bootstrapping is 
     used.}
   \item{\code{do.fit}:}{If \code{FALSE}, the model is not fit, and the 
     current starting values of the model parameters are preserved.}
}

Optimization options:

\describe{
    \item{\code{integration.ngh}:}{Integer. The number of Gauss-Hermite
      quadrature points (per dimension) used for numerical integration.
      The default is 21. Numerical integration is used when
      \code{estimator = "MML"}, and for two-level models with random
      slopes for \emph{latent} covariates (the \code{rv()} modifier;
      only the latent-covariate slopes are integrated numerically, all
      other random effects are handled in closed form).}
    \item{\code{control}:}{A list containing control parameters passed to
      the external optimizer. By default, lavaan uses \code{"nlminb"}. 
      See the manpage of \code{\link{nlminb}} for an overview of the control 
      parameters. If another (external) optimizer is selected, see the
      manpage for that optimizer to see the possible control parameters.}
    \item{\code{optim.method}:}{Character. The optimizer that should be
      used. For unconstrained optimization or models with only linear
      equality constraints (i.e., the model syntax
      does not include any "==", ">" or "<" operators),
      the available options are \code{"nlminb"} (the default), \code{"BFGS"},
      \code{"L-BFGS-B"}. These are all quasi-newton methods. A
      Gauss-Newton (Fisher scoring) method with Levenberg-Marquardt
      damping is also available (\code{optim.method = "GN"}); it can be
      used with the \code{"ML"}, \code{"GLS"}, \code{"WLS"},
      \code{"DWLS"}, \code{"ULS"}, \code{"DLS"}, \code{"NTRLS"} and
      \code{"catML"} estimators, and is the default when
      \code{estimator = "DLS"}. For multilevel models, an EM algorithm
      is also available (\code{optim.method = "em"}); this is the
      default for two-level models with \code{missing = "ml"} (as the
      EM algorithm is more robust in that setting; this also mimics
      Mplus); see the \code{em.args} option to control the EM
      algorithm. For two-level models with random slopes (the
      \code{rv()} modifier), \code{optim.method = "em"} uses a
      dedicated EM algorithm. Note that for any
      \code{optim.method = "em"} fit, variances that are fixed to zero
      are replaced by the (small) value \code{em.args$zerovar_offset}
      (as in Mplus, which uses the same 'minimum variance' device); as
      a result, the reported loglikelihood may differ slightly from an
      \code{optim.method = "nlminb"} fit of the same model (which uses
      the exact zero values).
      For constrained
      optimization, the only available option is \code{"nlminb_constr"},
      which uses an augmented Lagrangian minimization algorithm.}
    \item{\code{optim.force.converged}:}{Logical. If \code{TRUE}, pretend
      the model has converged, no matter what.}
    \item{\code{optim.fix.saturated}:}{Logical. Only used for multilevel
      models. If \code{TRUE} (the default), and one or more blocks (levels)
      of the model are saturated (i.e., all variances/covariances and
      means/intercepts of the observed variables of that block are free
      and unconstrained), the free parameters of these blocks are
      (temporarily) fixed to their h1 (unrestricted model) estimates during
      optimization, as their estimated values are already known. This often
      helps convergence. After the optimization, they are treated as free
      parameters again (e.g., when standard errors are computed).}
    \item{\code{optim.dx.tol}}{Numeric. Tolerance used for checking if
      the elements of the (unscaled) gradient are all zero (in absolute 
      value). The default value is 0.001.}
    \item{\code{gn.args}:}{List. Options related to the Gauss-Newton
      optimizer (\code{optim.method = "GN"}). The default values are
      \code{gn.args = list(max_iter = 200, tol_x = 1e-05,
      tol_g = 1e-06)}. The \code{max_iter} element sets the maximum
      number of Gauss-Newton iterations. Optimization stops when either
      the (scaled) gradient of the discrepancy function is smaller than
      \code{tol_g}, or the root mean square of the difference between
      the old and new parameter values is smaller than \code{tol_x}.}
    \item{\code{bounds}:}{Only used if \code{optim.method = "nlminb"}.
      If logical: \code{FALSE} implies no bounds are imposed on the parameters.
      If \code{TRUE}, this implies \code{bounds = "wide"}. If character,
      possible options are \code{"none"} (the default), \code{"standard"},
      \code{"wide"}, \code{"pos.var"}, \code{"pos.ov.var"}, and 
      \code{"pos.lv.var"}. 
      If \code{bounds = "pos.ov.var"}, the observed variances are forced to be
      nonnegative. If \code{bounds = "pos.lv.var"}, the latent variances are
      forced to be nonnegative. If \code{bounds = "pos.var"}, both observed
      and latent variances are forced to be nonnegative. If 
      \code{bounds = "standard"}, lower and upper bounds are computed for
      observed and latent variances, and factor loadings. If
      \code{bounds = "wide"}, lower and upper bounds are computed for
      observed and latent variances, and factor loadings; but the range of
      the bounds is enlarged (allowing again for slightly negative variances).}
    \item{\code{optim.bounds}:}{List. This can be used instead of the 
      \code{bounds} argument to allow more control. Possible elements of the
      list are \code{lower}, \code{upper}, \code{lower.factor} and 
      \code{upper.factor}. All of these accept a vector. The \code{lower} and
      \code{upper} elements indicate for which type of parameters bounds
      should be computed. Possible choices are \code{"ov.var"}, \code{"lv.var"},
      \code{"loadings"} and \code{"covariances"}. The \code{lower.factor} and 
      \code{upper.factor} elements should have the same length as the 
      \code{lower} and \code{upper} elements respectively. They indicate the
      factor by which the range of the bounds should be enlarged (for 
      example, 1.1 or 1.2; the default is 1.0). Other elements are
      \code{min.reliability.marker} which sets the lower bound for the
      reliability of the marker indicator (if any) of each factor 
      (default is 0.1). Finally, the \code{min.var.lv.endo} element indicates
      the lower bound of the variance of any endogenous latent variable
      (default is 0.0).}
    \item{\code{em.args}:}{List. Options related to the EM algorithm that
      is used to estimate the model parameters when
      \code{optim.method = "em"} (currently only available for multilevel
      models, with complete or missing data). The default values are
      \code{em.args = list(max_iter =
      10000, fx_tol = 1e-08, dx_tol = 1e-04, zerovar_offset = 1e-04,
      acceleration = "squarem", fused = TRUE)}.
      The \code{max_iter} element sets the maximum number of EM
      iterations. The EM algorithm stops when either the change in the
      loglikelihood is smaller than \code{fx_tol}, or the largest absolute
      change in the parameter values is smaller than \code{dx_tol}.
      The \code{zerovar_offset} element is used as the starting value for
      variances that are fixed to zero, as the EM algorithm cannot handle
      exact zero variances. The \code{acceleration} element controls
      whether the EM iterations are accelerated by the SQUAREM scheme,
      and the \code{fused} element controls whether the loglikelihood is
      computed as a byproduct of the E-step (see the \code{em.h1.args}
      option for details on both). If the \code{nlminb_handoff} element
      is \code{TRUE} (the default), and the EM iterations end with a
      non-negligible gradient (larger than \code{dx_tol}; this typically
      happens when the loglikelihood stalls before the gradient has
      converged), the solution is refined by a quasi-Newton (nlminb)
      run, warm-started at the EM values; the refined solution is only
      accepted if nlminb converges to an (at least) equally good
      solution; otherwise, the EM solution is retained. If, in addition,
      the \code{early_handoff} element is \code{TRUE} (the default), the
      EM phase stops early (using a loosened tolerance,
      \code{max(fx_tol, 0.001)}), as it only needs to provide a reliable
      warm start for nlminb; this is usually considerably faster, and
      gives the same solution; if nlminb fails from the early start, a
      safety net resumes the EM iterations with the strict \code{fx_tol}
      (and tries the hand-off once more), reproducing the non-early
      behavior. If \code{verbose = TRUE}, only the EM (or SQUAREM) cycles
      are shown; set the \code{verbose} element to \code{TRUE} to also show
      the optimizer output of each (inner) M-step.}
}

Parallelization options (currently only used for bootstrapping):

\describe{
    \item{parallel}{The type of parallel operation to be used (if any).  If
      missing, the default is \code{"no"}.  }
    \item{ncpus}{Integer: number of processes to be used in parallel operation:
      typically one would set this to the number of available CPUs. By
      default this is 2 (following CRAN policy); you may set this to a higher
      value (for example, the number of cores as detected by
      \code{parallel::detectCores()} minus one) if you wish to use more
      processes.}
    \item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
      \code{parallel = "snow"}. If not supplied, a cluster on the local 
      machine is created for the duration of the \code{lavBootstrap} 
      call.}
    \item{iseed}{An integer to set the seed. Or NULL if no reproducible 
      results are needed. This works for both serial (non-parallel) and 
      parallel settings. Internally, \code{RNGkind()} is set to 
      \code{"L'Ecuyer-CMRG"} if \code{parallel = "multicore"}. If 
      \code{parallel = "snow"} (under windows), 
      \code{parallel::clusterSetRNGStream()} is called which automatically 
      switches to \code{"L'Ecuyer-CMRG"}. When \code{iseed} is not
      NULL, \code{.Random.seed} (if it exists) in the global environment is
      left untouched.}
}

Categorical estimation options:

\describe{
    \item{\code{zero.add}:}{A numeric vector containing two values. These 
      values affect the calculation of polychoric correlations when some 
      frequencies in the bivariate table are zero. The first value only 
      applies for 2x2 tables. The second value for larger tables. This value 
      is added to the zero frequency in the bivariate table. If 
      \code{"default"}, the value is set depending on the \code{"mimic"}
      option. By default, lavaan uses \code{zero.add = c(0.5, 0.0)}.}
    \item{\code{zero.keep.margins}:}{Logical. This argument only affects 
      the computation of polychoric correlations for 2x2 tables with an empty 
      cell, and where a value is added to the empty cell. If \code{TRUE}, the 
      other values of the frequency table are adjusted so that all margins are 
      unaffected. If \code{"default"}, the value is set depending on the
      \code{"mimic"} option. The default is \code{TRUE}.}
    \item{\code{zero.cell.warn}:}{Logical. Only used if some observed 
      endogenous variables are categorical. If \code{TRUE}, give a warning if 
      one or more cells of a bivariate frequency table are empty.}
    \item{\code{allow.empty.cell}:}{Logical. If \code{TRUE}, ignore
      situations where an ordinal variable has fewer categories than
      expected, or where a category is empty in a specific group. This
    argument is currently used by the blavaan package (which
    uses lavaan to set up the model). The argument is not expected to
    salvage a lavaan model that results in errors.}
}

Starting values options:

\describe{
    \item{\code{start}:}{If it is a character string, the two options are 
      currently \code{"simple"} and \code{"default"}. In the first case, all 
      parameter values are set to zero, except the factor loadings and
      (residual) variances, which are set to 0.7 and 1.0 respectively.
      When \code{start} is \code{"default"}, the factor loadings are
      estimated using the fabin3 estimator (tsls) per factor. The
      residual variances of observed variables are set to half the
      observed variance, and all other (residual) variances are set to 0.05.
      The remaining parameters (regression coefficients, covariances) are
      set to zero.
      If \code{start} is a numerical vector, it should contain the starting
      values for all the free parameters (in the same order as the parameter
      table).
      If \code{start} is a fitted object of class \code{\linkS4class{lavaan}},
      the estimated values of the corresponding parameters will be extracted.
      If it is a parameter table, for example the output of the
      \code{parameterEstimates()} function, the values of the \code{est} or
      \code{start} or \code{ustart} column (whichever is found first) will be
      extracted.}
    \item{\code{rstarts}:}{Integer. The number of refits that lavaan should
      try with random starting values. Random starting values are computed
      by drawing random numbers from a uniform distribution. Correlations
      are drawn from the interval [-0.5, +0.5] and then converted to
      covariances. Lower and upper bounds for (residual) variances are
      computed just like the standard bounds in bounded estimation.
      Random starting values are not computed for regression coefficients
      (which are always zero) and factor loadings of higher-order constructs
      (which are always unity). From all the runs that converged, the final
      solution is the one that resulted in the smallest value for the
      discrepancy function.}
}

Check options:

\describe{
    \item{\code{check.start}:}{Logical. If \code{TRUE},
      the starting values are checked for possibly 
      inconsistent values (for example values implying correlations larger 
      than one). If needed, a warning is given.}
    \item{\code{check.sigma.pd}:}{Character. If \code{"chol"}, a Cholesky
      decomposition is used to check if the model implied covariance matrix
      (Sigma) is positive definite. If \code{"eigen"}, an eigen decomposition 
      is used to check if Sigma is positive definite. The latter was the
      default until 0.6-20. From 0.6-21 onwards, \code{"chol"} became the
      default (as it is much faster).}
    \item{\code{check.gradient}:}{Logical. If \code{TRUE}, and the model
      converged, a warning is given if the optimizer reported a (local)
      solution while not all elements of the (unscaled) gradient (as
      seen by the optimizer) are (near) zero, as they should be (the
      tolerance used is 0.001).}
    \item{\code{check.post}:}{Logical. If \code{TRUE}, and the model
      converged, a check is performed after (post) fitting, to verify if 
      the solution is admissible. This implies that all variances are
      non-negative, and all the model-implied covariance matrices are 
      positive (semi-)definite. For the latter test, we tolerate a tiny 
      negative eigenvalue that is smaller than .Machine$double.eps^(3/4), 
      treating it as being zero.}
    \item{\code{check.vcov}:}{Logical. If \code{TRUE}, and the model converged,
      we check if the variance-covariance matrix of the free parameters
      is positive definite. We take into account possible equality and
      active inequality constraints. If needed, a warning is given.}
    \item{\code{check.lv.names}:}{Logical. If \code{TRUE}, and latent variables
      are defined in the model, lavaan will stop with an error message if
      a latent variable name also occurs in the data (implying it is also
      an observed variable).}
}

Verbosity options:

\describe{
    \item{\code{verbose}:}{If \code{TRUE}, show what lavaan is doing. During
      estimation, the function value is printed out 
      during each iteration.}
    \item{\code{warn}:}{If \code{FALSE}, suppress all lavaan-specific 
      warning messages.}
    \item{\code{debug}:}{If \code{TRUE}, debugging information is printed 
      out.}
}

Miscellaneous:

\describe{
    \item{\code{model.type}:}{Set the model type: possible values
      are \code{"cfa"}, \code{"sem"} or \code{"growth"}. This may affect
      how starting values are computed, and may be used to alter the terminology
      used in the summary output, or the layout of path diagrams that are
      based on a fitted lavaan object.}
    \item{\code{mimic}:}{If \code{"Mplus"}, an attempt is made to mimic the
      Mplus program. If \code{"EQS"}, an attempt is made to mimic the EQS
      program. If \code{"default"}, the value is (currently) set to
      \code{"lavaan"}, which is very close to \code{"Mplus"}. See
      \code{\link{mimic}} for the list of options that are affected by each
      setting.}
    \item{\code{representation}:}{If \code{"LISREL"} the classical LISREL 
      matrix representation is used to represent the model (using the all-y 
      variant). No other options are available (for now).}
    \item{\code{implied}:}{Logical. If \code{TRUE}, compute the model-implied
      statistics, and store them in the implied slot.}
    \item{\code{h1}:}{Logical. If \code{TRUE}, compute the unrestricted model
      and store the unrestricted summary statistics (and perhaps a 
      loglikelihood) in the h1 slot.}
    \item{\code{baseline:}}{Logical. If \code{TRUE}, compute a baseline model
      (currently always the independence model, assuming all variables
      are uncorrelated) and store the results in the baseline slot.}
    \item{\code{baseline.conditional.x.free.slopes}:}{Logical. If \code{TRUE},
      and \code{conditional.x = TRUE}, the (default) baseline model will
      allow the slope structure to be unrestricted.}
    \item{\code{fit.by.level}:}{Logical. Only used for multilevel models.
      If \code{TRUE}, fit (and store) the partially saturated models that
      are needed to compute level-specific fit measures: for each level,
      an auxiliary model is fitted where the model at that level is
      retained, while the model at the other level is saturated (in
      addition, level-specific baseline models are fitted as well). See
      the \code{level} argument of \code{\link{fitMeasures}}. Set to
      \code{FALSE} to skip these additional model fits.}
    \item{\code{store.vcov}}{Logical. If \code{TRUE}, and \code{se=} is not
      set to \code{"none"}, store the full variance-covariance matrix of
      the model parameters in the vcov slot of the fitted lavaan object.}
    \item{\code{parser}}{Character. If \code{"new"} (the default), the new
      parser is used to parse the model syntax. If \code{"old"}, the original
      (pre 0.6-18) parser is used.}
}

}

\seealso{
\code{\link{lavaan}}, \code{\link{mimic}}, \code{\link{estimator_iv}}, \code{\link{estimator_js}}
}

\examples{
lavOptions()
lavOptions("std.lv")
lavOptions(c("std.lv", "orthogonal"))
}
