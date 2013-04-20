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
\alias{inspect}
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
    \item{\code{call}:}{The function call as returned by \code{match.called()}.}
    \item{\code{timing}:}{The elapsed time (user+system) for various parts of 
      the program as a list, including the total time.}
    \item{\code{Options}:}{Named list of options that were provided by
      the user, or filled-in automatically.}
    \item{\code{ParTable}:}{Named list describing the model parameters. Can be coerced to a data.frame. In the documentation, this is called the `parameter table'.}
    \item{\code{Data}:}{Object of internal class \code{"Data"}: information
about the data.}
    \item{\code{SampleStats}:}{Object of internal class \code{"SampleStats"}: sample 
      statistics}
    \item{\code{Model}:}{Object of internal class \code{"Model"}: the 
      internal (matrix) representation of the model}
    \item{\code{Cache}:}{List using objects that we try to compute only once, and reuse many times.}
    \item{\code{Fit}:}{Object of internal class \code{"Fit"}: the 
      results of fitting the model}
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
      If code{type="cor"}, the observed and model implied covariance matrix 
      is first transformed to a correlation matrix (using \code{cov2cor}), 
      before the residuals are computed.
      If \code{type="normalized"}, the residuals are
      normalized. If \code{type="standardized"}, the residuals are
      standardized. In the latter case, the residuals have a metric similar
      to z-values.}
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
      is the main `extractor' function for \code{lavaan} objects.  It
      allows the user to peek into the internal representation of the model. In
      addition, the requested information is returned (typically as a list) and
      can be used for further processing.
      The following values for \code{what} are allowed:
      \describe{
        \item{\code{"free"}:}{A list of model matrices counting the free 
          parameters in the model, typically in the same order as they are 
          specified by the user.}
        \item{\code{"start"}:}{A list of model matrices containing the
          starting values for all parameters.}
        \item{\code{"starting.values"}:}{An alias for \code{"start"}.}
        \item{\code{"sampstat"}:}{The sample statistics used for the analysis.}
        \item{\code{"se"}:}{A list of model matrices containing the 
          estimated standard errors for all free parameters.}
        \item{\code{"std.err"}:}{An alias for \code{"se"}.}
        \item{\code{"standard.errors"}:}{An alias for \code{"se"}.}
        \item{\code{"coef"}:}{A list of model matrices containing the current 
          values of all parameters.}
        \item{\code{"coefficients"}:}{An alias for \code{"coef"}.}
        \item{\code{"parameters"}:}{An alias for \code{"coef"}.}
        \item{\code{"parameter.estimates"}:}{An alias for \code{"coef"}.}
        \item{\code{"parameter.values"}:}{An alias for \code{"coef"}.}
        \item{\code{"estimates"}:}{An alias for \code{"coef"}.}
        \item{\code{"est"}:}{An alias for \code{"coef"}.}
        \item{\code{"x"}:}{An alias for \code{"coef"}.}
        \item{\code{"std.coef"}:}{A data.frame containing both the raw
           and (completely) standardized parameter values.}
        \item{\code{"std"}:}{An alias for \code{"std.coef"}.}
        \item{\code{"standardized"}:}{An alias for \code{"std.coef"}.}
        \item{\code{"standardizedsolution"}:}{An alias for \code{"std.coef"}.}
        \item{\code{"standardized.solution"}:}{An alias for \code{"std.coef"}.}
        \item{\code{"rsquare"}:}{A named vector with the R-Square value
          of the dependent observed and latent variables.}
        \item{\code{"r-square"}:}{An alias for \code{"rsquare"}.}
        \item{\code{"r2"}:}{An alias for \code{"rsquare"}.}
        \item{\code{"dx"}:}{A list of model matrices containing the 
          derivatives of all parameters evaluated at the current (typically 
          minimum) function value.}
        \item{\code{"gradient"}:}{An alias for \code{"dx"}.}
        \item{\code{"derivatives"}:}{An alias for \code{"dx"}.}
        \item{\code{"mi"}:}{A data.frame containing modification indices
          and expected parameter change (EPC) values, both in unstandardized
          and standardized metric.}
        \item{\code{"modindices"}:}{An alias for \code{"mi"}.}
        \item{\code{"modification"}:}{An alias for \code{"mi"}.}
        \item{\code{"modificationindices"}:}{An alias for \code{"mi"}.}
        \item{\code{"modification.indices"}:}{An alias for \code{"mi"}.}
        \item{\code{"converged"}:}{Returns TRUE if the optimization routine
          has converged; FALSE otherwise.}
        \item{\code{"list"}:}{A dataframe showing the internal 
          representation of a lavaan model. 
          Each row corresponds to a model parameter. The 
          columns contain all the information that lavaan stores about 
          these parameters (for example, if it is free of fixed, the
          user-specified starting values, etcetera). This is called 
          the parameter table.}
       }
    }
    \item{show}{\code{signature(object = "lavaan")}: Print a short summary 
      of the model fit}
    \item{summary}{\code{signature(object = "lavaan", standardized=FALSE, 
                                   fit.measures=FALSE, 
                                   rsquare=FALSE, 
                                   modindices=FALSE)}: 
      Print a nice summary of the model estimates. If \code{standardized=TRUE},
      the standardized solution is also printed. If \code{fit.measures=TRUE}, 
      the chi-square statistic is supplemented by several fit measures. 
      If \code{rsquare=TRUE}, the R-Square values for the dependent variables
      in the model are printed. If \code{modindices=TRUE}, modification indices
      are printed for all fixed parameters. Nothing is returned (use 
      \code{inspect} or another extractor function 
      to extract information from a fitted model).}
  }
}
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. URL
http://www.jstatsoft.org/v48/i02/.}
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
