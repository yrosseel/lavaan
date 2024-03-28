\name{lavaan-class}
\docType{class}
\alias{lavaan-class}
\alias{coef,lavaan-method}
\alias{vcov,lavaan-method}
\alias{anova,lavaan-method}
\alias{predict,lavaan-method}
\alias{resid,lavaan-method}
\alias{residuals,lavaan-method}
\alias{fitted,lavaan-method}
\alias{fitted.values,lavaan-method}    %% NOTE: no white space here!
\alias{nobs} %% not imported in < 2.13
\alias{nobs,lavaan-method}
\alias{logLik,lavaan-method}
\alias{update,lavaan-method}
\alias{show,lavaan-method}
\alias{summary,lavaan-method}

\title{Class For Representing A (Fitted) Latent Variable Model}
\description{The \code{lavaan} class represents a (fitted) latent variable
  model. It contains a description of the model as specified by the user,
  a summary of the data, an internal matrix representation, and if the model
  was fitted, the fitting results.}
\section{Objects from the Class}{
Objects can be created via the
  \code{\link{cfa}}, \code{\link{sem}}, \code{\link{growth}} or
  \code{\link{lavaan}} functions.
}
\section{Slots}{
  \describe{
    \item{\code{version}:}{The lavaan package version used to create this objects}
    \item{\code{call}:}{The function call as returned by \code{match.call()}.}
    \item{\code{timing}:}{The elapsed time (user+system) for various parts of
      the program as a list, including the total time.}
    \item{\code{Options}:}{Named list of options that were provided by
      the user, or filled-in automatically.}
    \item{\code{ParTable}:}{Named list describing the model parameters. Can be coerced to a data.frame. In the documentation, this is called the `parameter table'.}
    \item{\code{pta}:}{Named list containing parameter table attributes.}
    \item{\code{Data}:}{Object of internal class \code{"Data"}: information
about the data.}
    \item{\code{SampleStats}:}{Object of internal class \code{"SampleStats"}: sample
      statistics}
    \item{\code{Model}:}{Object of internal class \code{"Model"}: the
      internal (matrix) representation of the model}
    \item{\code{Cache}:}{List using objects that we try to compute only once, and reuse many times.}
    \item{\code{Fit}:}{Object of internal class \code{"Fit"}: the
      results of fitting the model. No longer used.}
    \item{\code{boot}:}{List. Results and information about the bootstrap.}
    \item{\code{optim}:}{List. Information about the optimization.}
    \item{\code{loglik}:}{List. Information about the loglikelihood of the model (if maximum likelihood was used).}
    \item{\code{implied}:}{List. Model implied statistics.}
    \item{\code{vcov}:}{List. Information about the variance matrix (vcov) of the model parameters.}
    \item{\code{test}:}{List. Different test statistics.}
    \item{\code{h1}:}{List. Information about the unrestricted h1 model (if available).}
    \item{\code{baseline}:}{List. Information about a baseline model (often
the independence model) (if available).}
    \item{\code{internal}:}{List. For internal use only.}
    \item{\code{external}:}{List. Empty slot to be used by add-on packages.}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "lavaan", type = "free")}: Returns
      the estimates of the parameters in the model as a named numeric vector.
      If \code{type="free"}, only the free parameters are returned.
      If \code{type="user"}, all parameters listed in the parameter table
      are returned, including constrained and fixed parameters.}
    \item{fitted.values}{\code{signature(object = "lavaan")}: Returns the
      implied moments of the model as a list with two elements (per group):
      \code{cov} for the implied covariance matrix,
      and \code{mean} for the implied mean
      vector. If only the covariance matrix was analyzed, the implied mean
      vector will be zero.}
    \item{fitted}{\code{signature(object = "lavaan")}: an alias for
        \code{fitted.values}.}
    \item{residuals}{\code{signature(object = "lavaan", type="raw")}: 
	  If \code{type = "raw"}, this function returns the raw (= unscaled)
      difference between the observed and the expected (model-implied) summary
      statistics.
	  If \code{type = "cor"}, or \code{type = "cor.bollen"}, the observed and
      model implied covariance matrices are first transformed to a correlation 
      matrix (using \code{cov2cor()}), before the residuals are computed.  
	  If \code{type = "cor.bentler"}, both the observed and model implied
      covariance matrices are rescaled by dividing the elements by the square 
      roots of the corresponding variances of the observed covariance matrix.  
      If \code{type="normalized"}, the residuals are divided by the square 
      root of the asymptotic variance of the corresponding summary statistic 
      (the variance estimate depends on the choice for the \code{se} argument).
      Unfortunately, the corresponding normalized residuals are not entirely
      correct, and this option is only available for historical interest.
      If \code{type="standardized"}, the residuals are divided by the square 
      root of the asymptotic variance of these residuals. The resulting 
      standardized residuals elements can be interpreted as z-scores.  
      If \code{type="standardized.mplus"}, the residuals are divided by the 
      square root of the asymptotic variance of these residuals. However, a 
      simplified formula is used (see the Mplus reference below) which often 
      results in negative estimates for the variances, resulting in many 
      \code{NA} values for the standardized residuals.}
    \item{resid}{\code{signature(object = "lavaan")}: an alias
      for \code{residuals}}
    \item{vcov}{\code{signature(object = "lavaan")}: returns the
      covariance matrix of the estimated parameters.}
    \item{predict}{\code{signature(object = "lavaan")}: compute
      factor scores for all cases that are provided in the data frame. For
      complete data only.}
    \item{anova}{\code{signature(object = "lavaan")}: returns
      model comparison statistics. This method is just a wrapper around
      the function \code{\link{lavTestLRT}}.
      If only a single argument (a fitted
      model) is provided, this model is compared to the unrestricted
      model. If two or more arguments (fitted models) are provided, the models
      are compared in a sequential order. Test statistics are based on the
      likelihood ratio test. For more details and
      further options, see the \code{\link{lavTestLRT}} page.}
    \item{update}{\code{signature(object = "lavaan", model, add, ...,
        evaluate=TRUE)}: update a fitted lavaan object and evaluate it
        (unless \code{evaluate=FALSE}). Note that we use the environment
        that is stored within the lavaan object, which is not necessarily
        the parent frame. The \code{add} argument is analogous to the one
        described in the \code{\link{lavTestScore}} page, and can be used to
        add parameters to the specified model rather than passing an entirely
        new \code{model} argument.}
    \item{nobs}{\code{signature(object = "lavaan")}: returns the effective
      number of observations used when fitting the model. In a multiple group
      analysis, this is the sum of all observations per group.}
    \item{logLik}{\code{signature(object = "lavaan")}:
      returns the log-likelihood of the fitted model, if maximum likelihood estimation
      was used. The \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}}
      methods automatically work via \code{logLik()}.}
    \item{show}{\code{signature(object = "lavaan")}: Print a short summary
      of the model fit}
    \item{summary}{\code{signature(object = "lavaan", header = TRUE,
     fit.measures = FALSE, estimates = TRUE, ci = FALSE, fmi = FALSE,
     standardized = FALSE, std.nox = FALSE,
     remove.step1 = TRUE, remove.unused = TRUE, cov.std = TRUE, rsquare = FALSE,
     modindices = FALSE, ci = FALSE, nd = 3L)}:
      Print a nice summary of the model estimates.
      If \code{header = TRUE}, the header section (including fit measures) is
      printed.
      If \code{fit.measures = TRUE}, additional fit measures are added to the
      header section. The related \code{fm.args} list allows to set options
      related to the fit measures. See \code{\link{fitMeasures}} 
      for more details.
      If \code{estimates = TRUE}, print the parameter estimates section.
      If \code{ci = TRUE}, add confidence intervals to the parameter estimates
      section.
      If \code{fmi = TRUE}, add the fmi (fraction of missing information)
      column, if it is available.
      If \code{standardized=TRUE} or a character vector, the standardized 
      solution is also printed (see \code{\link{parameterEstimates}}).  
      Note that \emph{SE}s and
      tests are still based on unstandardized estimates. Use
      \code{\link{standardizedSolution}} to obtain \emph{SE}s and test
      statistics for standardized estimates.
      The \code{std.nox} argument is deprecated; the \code{standardized}
      argument allows \code{"std.nox"} solution to be specifically requested.
      If \code{remove.step1}, the parameters of the measurement part are not
      shown (only used when using \code{sam()}.)
      If \code{remove.unused}, automatically added parameters that are
      fixed to their default (0 or 1) values are removed.
      If \code{rsquare=TRUE}, the R-Square values for the dependent variables
      in the model are printed.
      If \code{efa = TRUE}, EFA related information is printed. The related
      \code{efa.args} list allows to set options related to the EFA output.
      See \code{\link{summary.efaList}} for more details.
      If \code{modindices=TRUE}, modification indices
      are printed for all fixed parameters.
      The argument \code{nd} determines the number of digits after the
      decimal point to be printed (currently only in the parameter estimates
      section.) Historically, nothing was returned, but since 0.6-12, a
      list is returned of class \code{lavaan.summary} for which is print
      function is available.}
  }
}
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. 
\doi{https://doi.org/10.18637/jss.v048.i02}

Standardized Residuals in Mplus. Document retrieved from URL
https://www.statmodel.com/download/StandardizedResiduals.pdf
}
\seealso{
\code{\link{cfa}}, \code{\link{sem}},
\code{\link{fitMeasures}}, \code{\link{standardizedSolution}},
\code{\link{parameterEstimates}}, \code{\link{lavInspect}},
\code{\link{modindices}}
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)

summary(fit, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
fitted(fit)
coef(fit)
resid(fit, type = "normalized")
}
