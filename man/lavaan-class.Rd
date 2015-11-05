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
\alias{inspect,lavaan-method}
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
      results of fitting the model}
    \item{\code{boot}:}{List. Results and information about the bootstrap.}
    \item{\code{optim}:}{List. Information about the optimization.}
    \item{\code{implied}:}{List. Model implied statistics.}
    \item{\code{vcov}:}{List. Information about the variance matrix (vcov) of the model parameters.}
    \item{\code{test}:}{List. Different test statistics.}
    \item{\code{external}:}{List. Empty slot to be used by add-on packages.}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "lavaan", type = "free")}: Returns 
      the estimates of the parameters in the model as a named numeric vector.
      If \code{type="free"}, only the free parameters are returned. If
      \code{type="unco"}, both free and constrained parameters (simple
      equality constraints only) are returned.
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
      If \code{type="raw"}, this function returns the raw (=unstandardized)
      difference between the implied moments and the observed moments as 
      a list of two elements: \code{cov} for the residual covariance matrix, 
      and \code{mean} for the residual mean vector. 
      If only the covariance matrix was analyzed, the residual mean vector 
      will be zero. 
      If \code{type="cor"}, the observed and model implied covariance matrix 
      is first transformed to a correlation matrix (using \code{cov2cor}), 
      before the residuals are computed.
      If \code{type="normalized"}, the residuals are
      divided by the square root of an asymptotic variance estimate of the
      corresponding sample statistic (the variance estimate depends on the
      choice for the \code{se} argument).
      If \code{type="standardized"}, the residuals are divided by the square
      root of the difference between an asymptotic variance estimate of the
      corresponding sample statistic and an asymptotic variance estimate of
      the corresponding model-implied statistic.
      In the latter case, the residuals have a metric similar
      to z-values. On the other hand, they may often result in \code{NA} values;      for these cases, it may be better to use the normalized residuals. For
      more information about the normalized and standardized residuals, see
      the Mplus reference below.}
    \item{resid}{\code{signature(object = "lavaan")}: an alias
      for \code{residuals}}
    \item{vcov}{\code{signature(object = "lavaan")}: returns the
      covariance matrix of the estimated parameters.}
    \item{predict}{\code{signature(object = "lavaan")}: compute
      factor scores for all cases that are provided in the data frame. For
      complete data only.}
    \item{anova}{\code{signature(object = "lavaan")}: returns 
      model comparison statistics. See \code{\link[stats]{anova}}. At least
      two arguments (fitted models) are required. If the test statistic is
      scaled, an appropriate scaled difference test will be computed.}
    \item{update}{\code{signature(object = "lavaan", model.syntax, ..., 
        evaluate=TRUE)}: update a fitted lavaan object and evaluate it
        (unless \code{evaluate=FALSE}). Note that we use the environment
        that is stored within the lavaan object, which is not necessarily
        the parent frame.}
    \item{nobs}{\code{signature(object = "lavaan")}: returns the effective
      number of observations used when fitting the model. In a multiple group
      analysis, this is the sum of all observations per group.}
    \item{logLik}{\code{signature(object = "lavaan")}:
      returns the log-likelihood of the fitted model, if maximum likelihood estimation
      was used. The \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}}
      methods automatically work via \code{logLik()}.}
    \item{inspect}{\code{signature(object = "lavaan", what = "free")}: This
      method is now a shortcut for the \code{lavInspect()} function. See
      \link{lavInspect} for more details.}
    \item{show}{\code{signature(object = "lavaan")}: Print a short summary 
      of the model fit}
    \item{summary}{\code{signature(object = "lavaan", header = TRUE, fit.measures=FALSE, estimates = TRUE, ci = FALSE, fmi = FALSE, standardized = FALSE, rsquare=FALSE, std.nox = FALSE, modindices=FALSE, ci=FALSE, nd = 3L)}: 
      Print a nice summary of the model estimates. 
      If \code{header = TRUE}, the header section (including fit measures) is
      printed.  
      If \code{fit.measures = TRUE}, additional fit measures are added to the
      header section.
      If \code{estimates = TRUE}, print the parameter estimates section.
      If \code{ci = TRUE}, add confidence intervals to the parameter estimates
      section.
      If \code{fmi = TRUE}, add the fmi (fraction of missing information) 
      column, if it is available.
      If \code{standardized=TRUE},
      the standardized solution is also printed.  
      If \code{rsquare=TRUE}, the R-Square values for the dependent variables
      in the model are printed. 
      If \code{std.nox = TRUE}, the \code{std.all} column contains the 
      the \code{std.nox} column from the parameterEstimates() output.
      If \code{modindices=TRUE}, modification indices
      are printed for all fixed parameters. 
      The argument \code{nd} determines the number of digits after the
      decimal point to be printed (currently only in the parameter estimates
      section.
      Nothing is returned (use 
      \code{inspect} or another extractor function 
      to extract information from a fitted model).}
  }
}
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. URL
http://www.jstatsoft.org/v48/i02/.

Standardized Residuals in Mplus. Document retrieved from URL
http://www.statmodel.com/download/StandardizedResiduals.pdf
}
\seealso{
\code{\link{cfa}}, \code{\link{sem}}, \code{\link{growth}}, 
\code{\link{fitMeasures}}, \code{\link{standardizedSolution}},
\code{\link{parameterEstimates}},
\code{\link{modindices}}
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)

summary(fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)
inspect(fit, "free")
inspect(fit, "start")
inspect(fit, "rsquare")
inspect(fit, "fit")
fitted.values(fit)
coef(fit)
resid(fit, type="normalized")
}
