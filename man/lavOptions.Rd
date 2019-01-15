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
\item{mimic}{Character. Adjust the option list for this mimic flag.}
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
    \item{\code{std.lv}:}{If \code{TRUE}, the metric of each latent variable 
      is determined by fixing their (residual) variances to 1.0. If 
      \code{FALSE}, the metric of each latent variable is determined by fixing 
      the factor loading of the first indicator to 1.0.}
    \item{\code{parameterization}:}{Currently only used if data is 
      categorical. If \code{"delta"}, the delta parameterization is used. 
      If \code{"theta"}, the theta parameterization is used.}
}

Model features (only available for the \code{lavaan()} function):

\describe{
    \item{\code{auto.fix.first}:}{If \code{TRUE}, the factor loading of the 
      first indicator is set to 1.0 for every latent variable.}
    \item{\code{auto.fix.single}:}{If \code{TRUE}, the residual variance (if 
      included) of an observed indicator is set to zero if it is the only 
      indicator of a latent variable.}
    \item{\code{auto.var}:}{If \code{TRUE}, the residual variances and the 
      variances of exogenous latent variables are included in the model and 
      set free.}
    \item{\code{auto.cov.lv.x}:}{If \code{TRUE}, the covariances of exogenous 
      latent variables are included in the model and set free.}
    \item{\code{auto.cov.y}:}{If \code{TRUE}, the covariances of dependent 
      variables (both observed and latent) are included in the model and set 
      free.}
    \item{\code{auto.th}:}{If \code{TRUE}, thresholds for limited dependent 
      variables are included in the model and set free.}
    \item{\code{auto.delta}:}{If \code{TRUE}, response scaling parameters 
      for limited dependent variables are included in the model and set free.}
}

Data options:

\describe{
    \item{\code{std.ov}:}{If \code{TRUE}, all observed variables are 
      standardized before entering the analysis.}
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
}

Data summary options:

\describe{
    \item{\code{sample.cov.rescale}:}{If \code{TRUE}, the sample covariance 
      matrix provided by the user is internally rescaled by multiplying it 
      with a factor (N-1)/N.  If \code{"default"}, the value is set depending 
      on the estimator and the likelihood option: it is set to \code{TRUE} if 
      maximum likelihood estimation is used and \code{likelihood="normal"}, 
      and \code{FALSE} otherwise.}
    \item{\code{ridge}:}{Numeric. Small constant used for ridging. Only used 
      if the sample covariance matrix is non positive definite.}
}

Multiple group options:

\describe{      
    \item{\code{group.label}:}{A character vector. The user can specify 
      which group (or factor) levels need to be selected from the grouping 
      variable, and in which order. If missing, all grouping levels are 
      selected, in the order as they appear in the data.}
    \item{\code{group.equal}:}{A vector of character strings. Only used in
      a multiple group analysis. Can be one or more of the following:
      \code{"loadings"}, \code{"intercepts"}, \code{"means"},
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
      generalized least squares, \code{"WLS"} for weighted least squares 
      (sometimes called ADF estimation), \code{"ULS"} for unweighted least 
      squares and \code{"DWLS"} for diagonally weighted least squares. These 
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
    \item{\code{link}:}{Currently only used if estimator is MML. 
      If \code{"logit"}, a logit link is used for binary and ordered observed 
      variables. If \code{"probit"}, a probit link is used. If 
      \code{"default"}, it is currently set to \code{"probit"} (but this may 
      change).}
    \item{\code{information}:}{If \code{"expected"}, the expected 
      information matrix is used (to compute the standard errors). If 
      \code{"observed"}, the observed information matrix is used. 
      If \code{"first.order"}, the information matrix is based on the
      outer product of the casewise scores. See also the options
      \code{"h1.information"} and \code{"observed.information"} for
      further control. If \code{"default"}, the value is set depending 
      on the estimator, the missing argument, and the mimic option.}
   \item{\code{h1.information}:}{If \code{"structured"} (the default), the
      unrestricted (h1) information part of the (expected, first.order or 
      observed if h1 is used) information matrix is based on the structured, 
      or model-implied statistics (model-implied covariance matrix, 
      model-implied mean vector, etc.). 
      If \code{"unstructured"}, the unrestricted (h1) information part is
      based on sample-based statistics (observed covariance matrix, observed
      mean vector, etc.)}
   \item{\code{observed.information}:}{If \code{"hessian"}, the observed
      information matrix is based on the hessian of the objective function.
      If \code{"h1"}, an approximation is used that is based on 
      the observed information matrix of the unrestricted (h1) model.}
   \item{\code{se}:}{If \code{"standard"}, conventional standard errors
      are computed based on inverting the (expected, observed or first.order) 
      information matrix. If \code{"robust.sem"}, conventional robust
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
   \item{\code{test}:}{If \code{"standard"}, a conventional chi-square test 
     is computed. If \code{"Satorra.Bentler"}, a Satorra-Bentler scaled test 
     statistic is computed. If \code{"Yuan.Bentler"}, a Yuan-Bentler scaled 
     test statistic is computed. If \code{"Yuan.Bentler.Mplus"}, a 
     test statistic is computed that is asymptotically equal to the
     Yuan-Bentler scaled
     test statistic. If \code{"mean.var.adjusted"} or 
     \code{"Satterthwaite"}, a mean and variance adjusted test statistic is 
     compute. If \code{"scaled.shifted"}, an alternative mean and variance 
     adjusted test statistic is computed (as in Mplus version 6 or higher).
     If \code{"boot"} or \code{"bootstrap"} or \code{"Bollen.Stine"}, the 
     Bollen-Stine bootstrap is used to compute the bootstrap probability value 
     of the test statistic. If \code{"default"}, the value depends on the
     values of other arguments.}
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
      used. For unconstrained optimization (the model syntax
      does not include any "==", ">" or "<" operators),
      the available options are \code{"nlminb"} (the default), \code{"BFGS"} and
      \code{"L-BFGS-B"}. For constrained
      optimization, the only available option is \code{"nlminb.constr"}.}
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
      parameter values are set to zero, except the factor loadings (set to 
      one), the variances of latent variables (set to 0.05), and the residual 
      variances of observed variables (set to half the observed variance).
      If \code{"Mplus"}, we use a similar scheme, but the factor loadings are
      estimated using the fabin3 estimator (tsls) per factor. 
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
      the solution is admissable. This implies that all variances are
      non-negative, and all the model-implied covariance matrices are 
      positive (semi-)definite. For the latter test, we tolerate a tiny 
      negative eigenvalue that is smaller than .Machine$double.eps^(3/4), 
      treating it as being zero.}
    \item{\code{check.vcov}:}{Logical. If \code{TRUE}, and the model converged,
      we check if the variance-covariance matrix of the free parameters
      is positive definite. We take into account possible equality and
      acitive inequality constraints. If needed, a warning is given.}
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
    \item{\code{baseline:}}{Not used yet.}
    \item{\code{baseline.conditional.x.free.slopes}:}{Logical. If \code{TRUE},
      and \code{conditional.x = TRUE}, the (default) baseline model will
      allow the slopestructure to be unrestricted.}
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
