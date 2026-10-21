\name{lavBackendRegister}
\alias{lavBackendRegister}
\alias{lavBackends}
\title{Register an external estimation backend}
\description{
Register a backend that estimates models lavaan prepares but does not
estimate itself, so that it can be named in the \code{backend} argument
of \code{lavaan()}, \code{sem()}, \code{cfa()}, ... (see
\code{\link{lavOptions}}). \code{lavBackends()} lists the registered
backends.
}
\usage{
lavBackendRegister(name, fun, package = NULL, version = NULL,
                   supports = "all", capabilities = list(),
                   description = "", overwrite = TRUE)
lavBackends()
}
\arguments{
\item{name}{Character. The name of the backend, as used in the
    \code{backend} argument (matched without regard to case). The name
    \code{"lavaan"} is reserved.}
\item{fun}{A function \code{function(object, args)}. It receives the
    (unfitted) lavaan object that lavaan prepared for it, exactly as
    \code{do_fit = FALSE} would return it, and the \code{backend_args}
    list. It must return a list with the results of the estimation: the
    elements are the arguments of \code{\link{lavImportResults}} other
    than \code{object} (\code{est} is required; \code{vcov} or \code{se},
    \code{loglik}, \code{converged}, \code{iterations}, \code{status},
    \code{nobs}, \code{loglik_casewise}, \code{scores} and \code{object}
    (the backend's own result) are optional). lavaan then imports these
    results into the object.}
\item{package}{Character. The package that provides the backend, if
    any. It is loaded (with \code{requireNamespace()}) before the backend
    runs, also on parallel workers. A package registers its backend from
    its \code{.onLoad()} hook; naming a not-yet-loaded package as the
    backend loads it on demand.}
\item{version}{Character. A version string, shown by
    \code{lavBackends()}.}
\item{supports}{Character vector. The model classes the backend handles:
    \code{"all"}, \code{"dv"} (models with data-defined, \code{":~"},
    parameters) and/or \code{"standard"} (ordinary models). lavaan refuses
    a model the backend does not support before anything is estimated.}
\item{capabilities}{List. Free-form capabilities of the backend (for
    example whether it returns a vcov, casewise log-likelihoods, or
    scores); stored with the descriptor for use by lavaan or by the
    backend.}
\item{description}{Character. A one-line description.}
\item{overwrite}{Logical. If \code{FALSE}, registering a name that
    already exists is an error.}
}
\details{
The \code{backend} option is a string resolved through this registry, so
that a call, a saved object, \code{update()} and parallel code all carry
a name rather than an R closure. For interactive use, a bare function
can be passed as \code{backend}; lavaan registers it under an anonymous
name (\code{"function#1"}, ...) for the session.
}
\value{
\code{lavBackendRegister()} invisibly returns the descriptor (a list).
\code{lavBackends()} returns a data frame with one row per registered
backend.
}
\seealso{
\code{\link{lavImportResults}}, \code{\link{lavOptions}},
\code{\link{model.syntax}}.
}
\examples{
# a tiny backend: casewise ML by nlminb(), for the moderated-slope model
# y ~ b*x, b :~ b0 + b1*data(z)
cw_backend <- function(object, args) {
  pt <- parTable(object)
  dat <- lavInspect(object, "data")
  z <- object@Data@aux[[1]][, match("z", object@Data@ov.names.aux[[1]])]
  y <- dat[, "y"]; x <- dat[, "x"]; n <- length(y)
  free_rows <- which(pt$free > 0L)[order(pt$free[pt$free > 0L])]
  par_names <- pt$plabel[free_rows] # y~~y, y~1, b0, b1
  m2ll <- function(p) {
    mu <- p[2] + (p[3] + p[4] * z) * x
    sum((y - mu)^2) / p[1] + n * log(2 * pi * p[1])
  }
  opt <- nlminb(c(1, 0, 0, 0), m2ll, lower = c(1e-4, -Inf, -Inf, -Inf))
  list(est = setNames(opt$par, par_names),
       vcov = 2 * solve(optimHess(opt$par, m2ll)),
       loglik = -opt$objective / 2,
       converged = opt$convergence == 0L,
       iterations = opt$iterations)
}
lavBackendRegister("nlminb_cw", cw_backend, supports = "dv",
                   description = "casewise ML by nlminb (toy)")
lavBackends()

set.seed(1)
n <- 200
z <- rnorm(n); x <- rnorm(n)
y <- 3 + (0.5 + 0.3 * z) * x + rnorm(n)
dat <- data.frame(y = y, x = x, z = z)
fit <- sem(' y ~ b * x
             b :~ b0 + b1 * data(z) ', data = dat, backend = "nlminb_cw")
summary(fit)
}
