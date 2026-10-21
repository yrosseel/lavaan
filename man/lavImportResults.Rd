\name{lavImportResults}
\alias{lavImportResults}
\title{Import estimation results from an external backend}
\description{
Build a fitted \code{lavaan} object from a model that lavaan prepared
(with \code{do_fit = FALSE}) and the results of an estimation carried out
elsewhere, for example by an external package that estimates models with
data-defined (\code{":~"}) parameters using a casewise likelihood. The
parameter estimates, standard errors, log-likelihood and convergence
information are placed in the object, so that \code{summary()},
\code{coef()}, \code{vcov()}, \code{parameterEstimates()},
\code{logLik()} and the other post-fit functions can be used.
}
\usage{
lavImportResults(object, est, vcov = NULL, se = NULL, loglik = NULL,
                 converged = TRUE, iterations = NA_integer_,
                 backend = "external", status = NULL, nobs = NULL,
                 loglik_casewise = NULL, scores = NULL,
                 backend_object = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}, typically
    the (unfitted) object returned by \code{lavaan()}, \code{sem()},
    \code{cfa()}, ... with \code{do_fit = FALSE}.}
\item{est}{Numeric vector with the estimates of the free parameters.
    Either a named vector, where the names are the plabels (the
    \code{plabel} column of the parameter table) or the user-specified
    labels of the free parameters, or an unnamed vector of length equal to
    the number of free parameters, in the order of the \code{free} column
    of the parameter table.}
\item{vcov}{Optional variance-covariance matrix of the free parameters,
    one row and column per free parameter. If the matrix has dimnames,
    they are used (as for \code{est}) to put the rows and columns in
    lavaan's order; otherwise, free-parameter order is assumed.}
\item{se}{Optional numeric vector with standard errors (named or in
    free-parameter order, as for \code{est}); only used when
    \code{vcov} is not given. Standard errors of defined (\code{":="})
    parameters are only available when \code{vcov} is given.}
\item{loglik}{Optional single number: the log-likelihood at the
    estimates. Needed for \code{logLik()}, \code{AIC()} and
    \code{BIC()}.}
\item{converged}{Logical. Did the backend converge?}
\item{iterations}{Integer. The number of iterations the backend used (if
    known).}
\item{backend}{Character. The name of the backend; shown in the output.}
\item{status}{Optional: the backend's own status code or message
    (stored verbatim).}
\item{nobs}{Optional: the number of observations the backend used. If
    given, it must match the number of observations lavaan retained (see
    \code{lavInspect(object, "case.idx")}), since the backend must have
    analysed exactly these cases.}
\item{loglik_casewise}{Optional numeric vector with the casewise
    log-likelihood contributions (returned by
    \code{lavInspect(fit, "loglik_casewise")}).}
\item{scores}{Optional matrix with the casewise scores (one row per
    observation, one column per free parameter); stored for future use.}
\item{backend_object}{Optional: the backend's own result object; stored
    in the \code{@external} slot of the returned object.}
\item{\dots}{Not used.}
}
\details{
The object is rebuilt through lavaan's regular pipeline, with the
estimates as (fixed) starting values and without optimization, and the
slots that lavaan can not compute itself are then filled from the
supplied results. Before that, the estimates are checked against the
equality and inequality constraints of the model (a warning is given when
the backend did not honour them), \code{vcov} is checked for symmetry and
positive semi-definiteness, and the standard errors of defined
(\code{":="}) parameters are computed by the delta method.

For models with data-defined (\code{":~"}) parameters, the log-likelihood
is the casewise one supplied by the backend, the moderated (host)
parameters have no scalar estimate (\code{NA} in the parameter table), and
everything derived from a single set of model-implied moments remains
unavailable (see the \sQuote{Data-defined (moderated) parameters} section
in \code{\link{model.syntax}}). For ordinary models, lavaan computes the
model-implied moments, the log-likelihood and the test statistics itself
from the imported estimates; the \code{test} option of \code{object} is
honoured, and since \code{do_fit = FALSE} sets it to \code{"none"} unless
a test was requested explicitly, request the test when preparing the
object (for example \code{test = "standard"}).
}
\value{
An object of class \code{\linkS4class{lavaan}}.
}
\seealso{
\code{\link{lavBackendRegister}} (let lavaan call the backend itself,
through the \code{backend} option), \code{\link{model.syntax}},
\code{\link{lavInspect}}.
}
\examples{
set.seed(1)
n <- 200
z <- rnorm(n); x <- rnorm(n)
y <- 3 + (0.5 + 0.3 * z) * x + rnorm(n)
dat <- data.frame(y = y, x = x, z = z)

model <- ' y ~ b * x
           b :~ b0 + b1 * data(z) '
fit0 <- sem(model, data = dat, do_fit = FALSE)
pt <- parTable(fit0)

# a tiny external 'backend': casewise ML by nlminb()
free_rows <- which(pt$free > 0L)[order(pt$free[pt$free > 0L])]
par_names <- pt$plabel[free_rows] # .p2. (y~~y), .p4. (y~1), .p6. (b0), .p7. (b1)
m2ll <- function(p) {
  mu <- p[2] + (p[3] + p[4] * z) * x
  sum((y - mu)^2) / p[1] + n * log(2 * pi * p[1])
}
opt <- nlminb(c(1, 0, 0, 0), m2ll, lower = c(1e-4, -Inf, -Inf, -Inf))
est <- setNames(opt$par, par_names)
vc <- 2 * solve(optimHess(opt$par, m2ll))

fit <- lavImportResults(fit0, est = est, vcov = vc,
                        loglik = -opt$objective / 2,
                        iterations = opt$iterations, backend = "nlminb")
summary(fit)
logLik(fit)
}
