\name{lavEffects}
\alias{lavEffects}
\alias{laveffects}
\title{Indirect and Total Effects}
\description{
\sQuote{lavEffects} computes various \sQuote{effects} that are functions of
the estimated model parameters of a fitted lavaan object. In this version, the
focus is on the classic (LISREL-style) \emph{total}, \emph{indirect} and
\emph{direct} effects among the variables that appear in the structural part of
the model (that is, the variables that are involved in a regression). Both
observed and latent variables are supported.

This avoids the need to pre-specify all indirect (and total) effects manually
in the model syntax (using the \code{:=} operator), which can be tedious when
there are many of them.
}
\usage{
lavEffects(object, effects = c("total", "indirect"),
    se.def = NULL, level = 0.95, monte.carlo = NULL,
    boot.ci.type = "perc", zstat = TRUE, pvalue = TRUE, ci = TRUE,
    standardized = FALSE, cov.std = TRUE,
    add.class = TRUE, output = "data.frame")
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{effects}{Character vector. One or more of \code{"total"},
    \code{"indirect"} and \code{"direct"}. The default is
    \code{c("total", "indirect")}. The order of the elements does not matter;
    effects are always reported in the order total, indirect, direct.}
\item{se.def}{Character (or \code{NULL}). The method used to compute standard
    errors (and confidence intervals) for the effects. If \code{"monte.carlo"},
    the Monte Carlo method is used: parameters are drawn from a multivariate
    normal distribution with mean the parameter estimates and covariance matrix
    the estimated covariance matrix of the parameters, and the effects are
    recomputed for each draw (percentile confidence intervals). If
    \code{"delta"}, the delta method is used (based on a numerical Jacobian of
    the effects with respect to the free parameters; symmetric normal-theory
    confidence intervals). If \code{"bootstrap"}, the bootstrap draws that are
    stored in the fitted object (when the model was fitted with
    \code{se = "bootstrap"}) are used to compute bootstrap standard errors and
    bootstrap percentile confidence intervals. If \code{"none"}, no standard
    errors are computed.

    If \code{se.def = NULL} (the default), the choice is inherited from the
    fitted object: if the model was fitted with \code{se = "bootstrap"} (and
    bootstrap draws are available), the bootstrap method is used; else if the
    model was fitted with \code{se.def = "monte.carlo"}, the Monte Carlo method
    is used; otherwise the Monte Carlo method is used as the default. An
    explicit value always overrides this inheritance.}
\item{level}{Numeric. The confidence level for the confidence intervals.}
\item{monte.carlo}{List (or \code{NULL}). Settings for the Monte Carlo method
    (only used when \code{se.def = "monte.carlo"}), with elements \code{R} (the
    number of Monte Carlo draws) and \code{seed} (the random seed). If
    \code{NULL} (the default), the \code{monte.carlo} settings that were used
    when the model was fitted are inherited (which themselves default to
    \code{list(R = 20000L, seed = NULL)}). The state of the global random
    number generator is saved and restored, so the seed used here does not
    affect the user's random number stream.}
\item{boot.ci.type}{Character. Only used when \code{se.def = "bootstrap"}. The
    type of bootstrap confidence interval, as in \code{\link{parameterEstimates}}:
    \code{"norm"} (normal-theory using the bootstrap bias and standard error),
    \code{"basic"} (basic bootstrap), \code{"perc"} (the default; percentile
    method), \code{"bca.simple"} (bias-corrected percentile method, without the
    acceleration term) or \code{"bca"} (the adjusted bootstrap percentile
    method, with both bias and acceleration correction).}
\item{zstat}{Logical. If \code{TRUE}, an extra column is added with the
    \eqn{z}-statistic for each effect (the ratio of the estimate to its
    standard error).}
\item{pvalue}{Logical. If \code{TRUE}, an extra column is added with the
    (two-sided) p-value, computed using the standard normal distribution.}
\item{ci}{Logical. If \code{TRUE} (the default), two extra columns
    (\code{ci.lower} and \code{ci.upper}) are added with the lower and upper
    bounds of the confidence interval. If \code{FALSE}, the confidence interval
    bounds are not shown.}
\item{standardized}{Logical, character vector, or vector of (observed) variable
    names, behaving as in \code{\link{parameterEstimates}}. If \code{TRUE},
    standardized effects are added as extra columns (point estimates only):
    \code{std.lv} (only the latent variables are standardized), \code{std.all}
    (both observed and latent variables are standardized) and, if there are
    exogenous observed variables and \code{fixed.x = TRUE}, \code{std.nox}
    (like \code{std.all}, but the exogenous observed variables are not
    standardized). A character vector selects specific types (a subset of
    \code{"std.lv"}, \code{"std.all"}, \code{"std.nox"}); a vector of observed
    variable names standardizes only those observed variables (plus the latent
    variables), and the result is reported in a \code{std.user} column. The
    standardized effect of \eqn{j} on \eqn{i} equals the unstandardized effect
    times the ratio of the model-implied standard deviations
    \eqn{sd_j / sd_i}.}
\item{cov.std}{Logical. See \code{\link{parameterEstimates}}. (Has no influence
    on the standardized regression-type effects reported here.)}
\item{add.class}{Logical. If \code{TRUE}, the result (when
    \code{output = "data.frame"}) is given the class
    \code{c("lavaan.effects", "lavaan.data.frame", "data.frame")}, so that a
    dedicated print method is used.}
\item{output}{Character. If \code{"data.frame"} (the default), the result is a
    (classed) data.frame with one row per effect. If \code{"list"}, the result
    is a list of matrices (one per effect type, and one per block in the
    multilevel/multigroup case), where element \eqn{[i,j]} contains the effect
    of variable \eqn{j} on variable \eqn{i}.}
}
\details{
All effects are derived from the reduced-form matrix
\eqn{(I - B)^{-1}}, where \eqn{B} is the matrix of (direct) regression
coefficients among the structural variables (the \code{beta} matrix in the
LISREL representation that lavaan uses internally). Writing \eqn{B[i,j]} for the
direct effect of variable \eqn{j} on variable \eqn{i}, we have:
\itemize{
  \item the \emph{direct} effect of \eqn{j} on \eqn{i} equals \eqn{B[i,j]};
  \item the \emph{total} effect of \eqn{j} on \eqn{i} equals
        \eqn{\left((I - B)^{-1} - I\right)[i,j]};
  \item the \emph{indirect} effect of \eqn{j} on \eqn{i} equals the total
        effect minus the direct effect.
}
Only effects that are structurally non-zero (that is, for which a directed path
from \eqn{j} to \eqn{i} exists) are reported. For indirect effects, this means
that at least one path of length two or more must exist.

The delta and Monte Carlo methods only require the parameter estimates and the
estimated covariance matrix of the parameters. The delta method produces
(symmetric) normal-theory confidence intervals; the Monte Carlo method produces
percentile-based confidence intervals, which need not be symmetric around the
point estimate. The Monte Carlo method typically behaves better than the delta
method for effects that are products of parameters (such as indirect effects),
in particular in smaller samples. When the model was fitted with
\code{se = "bootstrap"}, the bootstrap draws are reused to obtain bootstrap
standard errors and bootstrap percentile confidence intervals (this requires
no additional model fitting).

Models fitted with \code{conditional.x = TRUE} are also supported: in that
case the structural model is \eqn{eta = B eta + Gamma x + zeta}, where the
exogenous covariates \eqn{x} are stored in the \code{gamma} matrix. The effects
of these covariates on the endogenous variables are then computed as
\eqn{(I - B)^{-1} Gamma} (total), \eqn{Gamma} (direct) and
\eqn{((I - B)^{-1} - I) Gamma} (indirect), and reported alongside the effects
among the \code{beta} variables.
}
\value{
If \code{output = "data.frame"} (the default), a data.frame (of class
\code{lavaan.effects}) with the following columns: \code{effect} (the type of
effect: total, indirect or direct), \code{lhs} (the outcome variable),
\code{op} (always \code{"~"}), \code{rhs} (the predictor variable), and,
depending on the arguments, \code{group}/\code{level}/\code{block},
\code{est}, \code{se}, \code{z}, \code{pvalue}, \code{ci.lower},
\code{ci.upper} and (if \code{standardized} is requested) one or more of
\code{std.lv}, \code{std.all}, \code{std.nox} and \code{std.user}. The effect
in a given row is the effect of \code{rhs} on \code{lhs}.

If \code{output = "list"}, a (possibly nested) list of matrices, where element
\eqn{[i,j]} contains the effect of variable \eqn{j} on variable \eqn{i}.
}
\seealso{
\code{\link{parameterEstimates}}, \code{\link{standardizedSolution}}.
}
\examples{
# a simple mediation model
set.seed(1234)
X <- rnorm(300)
M <- 0.5 * X + rnorm(300)
Y <- 0.4 * M + 0.3 * X + rnorm(300)
Data <- data.frame(X = X, Y = Y, M = M)
model <- ' # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
         '
fit <- sem(model, data = Data)

# total and indirect effects, with Monte Carlo standard errors
lavEffects(fit)

# all effects, using the delta method
lavEffects(fit, effects = c("total", "indirect", "direct"), se.def = "delta")

# add standardized (total and indirect) effects
lavEffects(fit, se.def = "delta", standardized = TRUE)
}
