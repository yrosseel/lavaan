\name{lavOptions}
\alias{lavOptions}
\alias{lavoptions}
\title{lavaan Options}
\description{
Show the default options used by the \code{lavaan()} function. The
options can be changed by passing 'name = value' arguments to the 
\code{lavaan()} function call, where they will be added to the '...'
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


Model features (always available):

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
      the model-implied statistics refect both sets of variables. If
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
     As a result, the latent mean will be freely estimated and usually 
     equal the average of the means of the involved indicators.}
    \item{\code{ceq.simple}:}{Logical. If \code{TRUE}, and no other 
     general (equality or inequality) constraints are used in the model, 
     simple equality constraints
     are represented in the parameter table as duplicated free parameters 
     (instead of extra rows with \code{op = "=="}).}
    \item{\code{parameterization}:}{Currently only used if data is 
      categorical. If \code{"delta"}, the delta parameterization is used. 
      If \code{"theta"}, the theta parameterization is used.}
    \item{\code{correlation}:}{(Experimental/incomplete). Only used for 
     continuous data. If \code{TRUE}, analyze a correlation matrix (instead
     of a (co)variance matrix). This implies that the residual observed
     variances are no longer free parameters. Instead, they are set to
     values to ensure the model-implied variances are unity. This also
     affects the standard errors.}
}

Model features (only available for the \code{lavaan()} function):

\describe{
    \item{\code{auto.fix.first}:}{If \code{TRUE}, the factor loading of the 
      first indicator is set to 1.0 for every latent variable.}
    \item{\code{auto.fix.single}:}{If \code{TRUE}, the residual variance (if 
      included) of an observed indicator is set to zero if it is the only 
      indicator of a latent variable.}
    \item{auto.var}{If \code{TRUE}, the (residual) variances of both observed
      and latent variables are set free.}
    \item{\code{auto.cov.lv.x}:}{If \code{TRUE}, the covariances of exogenous 
      latent variables are included in the model and set free.}
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
      Use this option with caution; it can be used to test if (for example)
      nonconvergence was due to scaling issues. But this is still a covariance
      based analysis, in the sense that no constraints are involved (to 
      ensure the model-implied (co)variance matrix has unit variances), and the
      standard errors still assume that the input was unstandardized. See
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
      If you wish to use multiple
      imputation, you need to use an external package (eg. mice) to
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
      those observations for which both values are observed (not-missing).
      And this may change from pair to pair.
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
      again by the group size. The default is \code{"total"}.}
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
    \item{\code{ridge}:}{Logical. If \code{TRUE} a small constant value will
      be added the diagonal elements of the covariance (or correlation) 
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
      \code{"loadings"}, \code{"composite.loadings"}, 
      \code{"intercepts"}, \code{"means"},
      \code{"thresholds"}, \code{"regressions"}, \code{"residuals"},
      \code{"residual.covariances"}, \code{"lv.variances"} or
      \code{"lv.covariances"}, specifying the pattern of equality
      constraints across multiple groups.}
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
      (\code{test="satorra.bentler"}) , \code{"MLMVS"} will compute a
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
      and a mean or mean and variance adjusted test statistic.}
    \item{\code{likelihood}:}{Only relevant for ML estimation. If 
      \code{"wishart"}, the wishart likelihood approach is used. In this 
      approach, the covariance matrix has been divided by N-1, and both 
      standard errors and test statistics are based on N-1.
      If \code{"normal"}, the normal likelihood approach is used. Here,
      the covariance matrix has been divided by N, and both standard errors
      and test statistics are based on N. If \code{"default"}, it depends
      on the mimic option: if \code{mimic="lavaan"} or \code{mimic="Mplus"},
      normal likelihood is used; otherwise, wishart likelihood is used.}
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
      test statistics from a fitted lavaan object.}
   \item{\code{scaled.test}:}{Character. Choose the test statistic
    that will be scaled (if a scaled test statistic is requested).
    The default is \code{"standard"}, but it could also be (for example)
    \code{"Browne.residual.nt"}.}
   \item{\code{gamma.n.minus.one}}{Logical. If \code{TRUE}, we divide the
     Gamma matrix by N-1 (instead of the default N).}
   \item{\code{gamma.unbiased}}{Logical. If \code{TRUE}, we compute an
     unbiased version for the Gamma matrix. Only available for single-level
     complete data and when \code{conditional.x = FALSE} and 
     \code{fixed.x = FALSE} (for now).}
   \item{\code{bootstrap}:}{Number of bootstrap draws, if bootstrapping is 
     used.}
   \item{\code{do.fit}:}{If \code{FALSE}, the model is not fit, and the 
     current starting values of the model parameters are preserved.}
}

Optimization options:

\describe{
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
      \code{"L-BFGS-B"}. These are all quasi-newton methods. A basic 
      implementation of Gauss-Newton is also available
      (\code{optim.method = "GN"}). The latter is the default when 
      \code{estimator = "DLS"}.
      For constrained
      optimization, the only available option is \code{"nlminb.constr"},
      which uses an augmented Lagrangian minimization algorithm.}
    \item{\code{optim.force.converged}:}{Logical. If \code{TRUE}, pretend
      the model has converged, no matter what.}
    \item{\code{optim.dx.tol}}{Numeric. Tolerance used for checking if
      the elements of the (unscaled) gradient are all zero (in absolute 
      value). The default value is 0.001.}
    \item{\code{optim.gn.tol.x}:}{Numeric. Only used when 
      \code{optim.method = "GN"}. Optimization stops when
      the root mean square of the difference between the old and new
      parameter values are smaller than this tolerance value. Default is
      \code{1e-05} for DLS estimation and \code{1e-07} otherwise.}
    \item{\code{optim.gn.iter.max}:}{Integer. Only used when 
      \code{optim.method = "GN"}. The maximum number of GN iterations. 
      The default is 200.}
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
      should be computed. Possible choice are \code{"ov.var"}, \code{"lv.var"},
      \code{"loadings"} and \code{"covariances"}. The \code{lower.factor} and 
      \code{upper.factor} elements should have the same length as the 
      \code{lower} and \code{upper} elements respectively. They indicate the
      factor by which the range of the bounds should be enlarged (for 
      example, 1.1 or 1.2; the default is 1.0). Other elements are
      \code{min.reliability.marker} which sets the lower bound for the
      reliability of the marker indicator (if any) of each factor 
      (default is 0.1). Finally, the \code{min.var.lv.endo} element indicates
      the lower bound of the variance of any endogenous latent variance
      (default is 0.0).}
}

Parallelization options (currently only used for bootstrapping):

\describe{
    \item{parallel}{The type of parallel operation to be used (if any).  If
      missing, the default is \code{"no"}.  }
    \item{ncpus}{Integer: number of processes to be used in parallel operation:
      typically one would chose this to the number of available CPUs. By
      By default this is the number of cores (as detected by 
      \code{parallel::detectCores()}) minus one.}
    \item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
      \code{parallel = "snow"}. If not supplied, a cluster on the local 
      machine is created for the duration of the \code{bootstrapLavaan} 
      or \code{bootstrapLRT} call.}
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
      option. By default, lavaan uses \code{zero.add = c(0.5. 0.0)}.}
    \item{\code{zero.keep.margins}:}{Logical. This argument only affects 
      the computation of polychoric correlations for 2x2 tables with an empty 
      cell, and where a value is added to the empty cell. If \code{TRUE}, the 
      other values of the frequency table are adjusted so that all margins are 
      unaffected. If \code{"default"}, the value is set depending on the 
      \code{"mimic"}. The default is \code{TRUE}.}
    \item{\code{zero.cell.warn}:}{Logical. Only used if some observed 
      endogenous variables are categorical. If \code{TRUE}, give a warning if 
      one or more cells of a bivariate frequency table are empty.}
}

Starting values options:

\describe{
    \item{\code{start}:}{If it is a character string, the two options are 
      currently \code{"simple"} and \code{"Mplus"}. In the first case, all 
      parameter values are set to zero, except the factor loadings and
      (residual) variances, which are set to one.
      When \code{start} is \code{"Mplus"}, the factor loadings are
      estimated using the fabin3 estimator (tsls) per factor. The
      residual variances of observed variables are set tot half the 
      observed variance, and all other (residual) variances are set to 0.05.
      The remaining parameters (regression coefficients, covariances) are
      set to zero.
      If \code{start} is a fitted object of class \code{\linkS4class{lavaan}},
      the estimated values of the corresponding parameters will be extracted.
      If it is a parameter table, for example the output of the
      \code{paramaterEstimates()} function, the values of the \code{est} or
      \code{start} or \code{ustart} column (whichever is found first) will be
      extracted.}
}

Check options:

\describe{
    \item{\code{check.start}:}{Logical. If \code{TRUE},
      the starting values are checked for possibly 
      inconsistent values (for example values implying correlations larger 
      than one). If needed, a warning is given.}
    \item{\code{check.gradient}:}{Logical. If \code{TRUE}, and the model
      converged, a warning
      is given if the optimizer decided that a (local) solution has
      been found, while not all elements of the (unscaled) gradient (as 
      seen by the optimizer) are (near) zero, as
      they should be (the tolerance used is 0.001).}
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
      acitive inequality constraints. If needed, a warning is given.}
    \item{\code{check.lv.names}:}{Logical. If \code{TRUE}, and latent variables
      are defined in the model, lavaan will stop with an error message if
      a latent variable name also occurs in the data (implying it is also
      an observed variable).}
}

Verbosity options:

\describe{
    \item{\code{verbose}:}{If \code{TRUE}, the function value is printed out 
      during each iteration.}
    \item{\code{warn}:}{If \code{TRUE}, some (possibly harmless) warnings 
      are printed out during the iterations.}
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
      program. If \code{"default"}, the value is (currently) set to to 
      \code{"lavaan"}, which is very close to \code{"Mplus"}.}
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
      allow the slopestructure to be unrestricted.}
    \item{\code{store.vcov}}{Logical. If \code{TRUE}, and \code{se=} is not
      set to \code{"none"}, store the full variance-covariance matrix of
      the model parameters in the vcov slot of the fitted lavaan object.}
    \item{\code{parser}}{Character. If \code{"new"} (the default), the new
      parser is used to parse the model syntax. If \code{"old"}, the original
      (pre 0.6-18) parser is used.}
}

}

\seealso{
\code{\link{lavaan}}
}

\examples{
lavOptions()
lavOptions("std.lv")
lavOptions(c("std.lv", "orthogonal"))
}
