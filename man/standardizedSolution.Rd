\name{standardizedSolution}
\alias{standardizedSolution}
\alias{standardizedsolution}
\title{Standardized Solution}
\description{
Standardized solution of a latent variable model.}
\usage{
standardizedSolution(object, type = "std.all", se = TRUE, zstat = TRUE,
                     pvalue = TRUE, ci = TRUE, level = 0.95,
                     boot_ci_type = "perc", cov_std = TRUE,
                     remove_eq = TRUE, remove_ineq = TRUE, remove_def = FALSE,
                     remove_aux = TRUE,
                     partable = NULL, glist = NULL, est = NULL,
                     output = "data.frame", ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{If \code{"std.lv"}, the standardized estimates are based
on the variances of the (continuous) latent variables only.
If \code{"std.all"}, the standardized estimates are based
on the variances of both (continuous) observed and latent variables.
If \code{"std.nox"}, the standardized estimates are based
on the variances of both (continuous) observed and latent variables,
but not the variances of exogenous covariates.
Alternatively, \code{type} may be a vector of (observed) variable names
(for example \code{type = c("x1", "x2")}); in that case only the parameters
involving these variables are standardized (the other observed variables are
left unstandardized). This is a generalization of \code{"std.nox"}, where the
(observed) exogenous \code{x} variables are the ones left unstandardized.}
\item{se}{Logical. If TRUE, standard errors for the standardized parameters
will be computed, together with a z-statistic and a p-value.}
\item{zstat}{Logical. If \code{TRUE}, an extra column is added containing
the so-called z-statistic, which is simply the value of the estimate divided
by its standard error.}
\item{pvalue}{Logical. If \code{TRUE}, an extra column is added containing
the pvalues corresponding to the z-statistic, evaluated under a standard
normal distribution.}
\item{ci}{If \code{TRUE}, confidence intervals are added to
the output.}
\item{level}{The confidence level required.}
\item{boot_ci_type}{Character. Only used if the model was fitted with
\code{se = "bootstrap"}. In that case, bootstrap standard errors and bootstrap
confidence intervals are computed for the standardized parameters (by
re-standardizing each bootstrap draw), and this argument selects the type of
bootstrap confidence interval, exactly as in
\code{\link{parameterEstimates}}: \code{"norm"}, \code{"basic"},
\code{"perc"} (the default), \code{"bca.simple"} or \code{"bca"}.}
\item{cov_std}{Logical. If TRUE, the (residual) observed
covariances are scaled by the square root of the `Theta' diagonal elements, and
the (residual) latent covariances are scaled by the square root of the `Psi'
diagonal elements. If FALSE, the (residual) observed covariances
are scaled by the square root of the diagonal elements of the observed
model-implied covariance matrix (Sigma), and the (residual) latent covariances
are scaled by the square root of diagonal elements of the model-implied
covariance matrix of the latent variables.}
\item{remove_eq}{Logical. If TRUE, filter the output by removing all rows
containing equality constraints, if any.}
\item{remove_ineq}{Logical. If TRUE, filter the output by removing all
rows containing inequality constraints, if any.}
\item{remove_def}{Logical. If TRUE, filter the output by removing all
rows containing parameter definitions, if any.}
\item{remove_aux}{Logical. If TRUE (the default), filter the output by
removing all rows corresponding to auxiliary (\code{aux=}) variables added to
the model as saturated correlates (only relevant when \code{missing = "ml"}
and the \code{aux=} argument was used).}
\item{glist}{List of model matrices. If provided, they will be used
instead of the GLIST inside the object@Model slot. Only works if the
\code{est} argument is also provided. See Note.}
\item{est}{Numeric. Parameter values (as in the `est' column of a
parameter table). If provided, they will be used instead of
the parameters that can be extracted from object. Only works if the \code{glist}
argument is also provided. See Note.}
\item{partable}{A custom \code{list} or \code{data.frame} in which to store
the standardized parameter values. If provided, it will be used instead of
the parameter table inside the object@ParTable slot.}
\item{output}{Character. If \code{"data.frame"}, the parameter table is
displayed as a standard (albeit lavaan-formatted) data.frame.
If \code{"text"} (or alias \code{"pretty"}), the parameter table is
prettified and displayed with subsections (as used by the summary function).}
\item{...}{To support old argument names.}
}
\details{
The standardized estimates are functions of the (unstandardized) free
parameters and of (a subset of) the model-implied (co)variances used for
scaling. The standard errors reported by \code{standardizedSolution} are
therefore standard errors \emph{for the standardized parameters}, and they
will in general differ from the standard errors of the unstandardized
parameters reported by \code{\link{parameterEstimates}} (or in the
\code{summary()} output). They are also not the same as the standard errors
that would be obtained by simply rescaling the unstandardized standard errors.

How the standard errors are computed depends on how the model was originally
fitted (in particular on the \code{se=} argument of \code{\link{lavaan}},
\code{\link{cfa}}, \code{\link{sem}}, \dots):
\itemize{
  \item By default (\code{se = "standard"}, \code{"robust.sem"},
    \code{"robust.huber.white"}, \dots), the standard errors are obtained with
    the \emph{delta method}: the Jacobian of the standardization function is
    computed (numerically) and combined with the variance-covariance matrix of
    the (unstandardized) parameter estimates. Any robustness present in that
    variance-covariance matrix (for example robust or sandwich-type standard
    errors) is automatically propagated to the standardized solution.
  \item If the model was fitted with \code{se = "bootstrap"} (and the bootstrap
    draws are available), bootstrap standard errors and bootstrap confidence
    intervals are reported instead. These are obtained by re-standardizing each
    bootstrap draw and computing the standard deviation (for the standard
    error) and the requested interval type (see \code{boot_ci_type}) of the
    resulting standardized values. Note that, as in
    \code{\link{parameterEstimates}}, the p-value is still computed by
    referring the z-statistic (standardized estimate divided by its bootstrap
    standard error) to a standard normal distribution.
}

There is no separate argument to choose the \emph{type} of standard error
within \code{standardizedSolution}: the type always follows the \code{se=}
setting that was used when the model was fitted. To obtain, say, robust or
bootstrap standard errors for the standardized solution, refit the model with
the corresponding \code{se=} argument.
}
\note{
The \code{est}, \code{glist}, and \code{partable} arguments are not meant for
everyday users, but for authors of external R packages that depend on
\code{lavaan}. Only to be used with great caution.
}
\value{
  A data.frame containing standardized model parameters.

  If the model was fitted with \code{se = "bootstrap"} (and the bootstrap
  draws are available), the standardized estimates in the bootstrap samples
  are stored as the \code{"boot_std"} attribute of the returned data.frame (a
  matrix with one row per bootstrap draw and one column per row of the
  returned solution). These can be retrieved with
  \code{attr(x, "boot_std")}, analogous to
  \code{lavInspect(object, "coef.boot")} for the unstandardized solution.

  In addition, if \code{boot_ci_type = "bca"}, the empirical-influence design
  matrix used to compute the acceleration constant is stored as the
  \code{"design"} attribute (retrieve with \code{attr(x, "design")}).
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
standardizedSolution(fit)
}
