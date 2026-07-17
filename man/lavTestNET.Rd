\name{lavTestNET}
\alias{lavTestNET}
\title{Nesting and Equivalence Testing (NET)}
\description{
Checks -- for every pair of fitted models -- whether the more restricted
model is nested within the less restricted model \emph{in the covariance
sense}, using the NET procedure of Bentler & Satorra (2010); for pairs of
models with equal degrees of freedom, it checks whether the models are
equivalent.}
\usage{
lavTestNET(object, ..., crit = 1e-7, npoints = 1L,
           iseed = 12345L, model_names = NULL)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{...}{additional (at least one) objects of class
  \code{\linkS4class{lavaan}}.}
\item{crit}{Numeric. The criterion applied to the minimum of the fit
  function when the less restricted model is fitted to the model-implied
  moments of the more restricted model: values below \code{crit} imply
  nesting (or equivalence, if the degrees of freedom are equal), values
  above \code{10 * crit} imply non-nesting, and values in between are
  inconclusive (following Asparouhov & \enc{Muthén}{Muthen}, 2019).}
\item{npoints}{Integer. The number of parameter points of the restricted
  model where the check is carried out. The first point is always the
  vector of parameter estimates. If \code{npoints > 1L}, additional
  randomly perturbed parameter values are used (respecting any equality
  and inequality constraints). This protects against \sQuote{data
  specific} conclusions: two models may accidentally appear to be nested
  (or equivalent) at the estimated parameter values only. See the
  discussion in Asparouhov & \enc{Muthén}{Muthen} (2019).}
\item{iseed}{Integer. The seed for the random perturbations (only used if
  \code{npoints > 1L}). The state of the random generator of the user is
  left untouched.}
\item{model_names}{Character vector. If provided, use these model names
  in the rows and columns of the result.}
}
\details{
The procedure works as follows: the model-implied (first and) second
order moments of the more restricted model are treated as if they were
sample statistics, and the less restricted model is fitted to these
moments. If (and only if) the minimum of the fit function is zero, the
less restricted model can exactly reproduce these moments, and the more
restricted model is nested within the less restricted model at this
parameter point. If both models have the same degrees of freedom and each
model is nested within the other, the models are equivalent.

Whenever the nesting can be established parametrically (the parameter
table of the restricted model equals the parameter table of the less
restricted model with some free parameters fixed and/or constrained to be
equal), no refitting is needed, and the (numerically exact) verdict is
based on the parameter tables alone.

Note that nesting in the covariance sense does not by itself guarantee
that the usual chi-squared difference test behaves as a chi-square
variate: if the nesting only holds on the boundary of the parameter space
(for example, a variance fixed to zero), the standard regularity
conditions fail.

The procedure is not available for multilevel models, and (currently) not
for categorical models with \code{conditional.x = TRUE}, correlation
structures, or models involving auxiliary variables.
}
\value{
A list of class \code{lavaan.net} with elements \code{nested} (a logical
matrix; entry \eqn{[i, j]} indicates whether the model in row \eqn{i} is
nested within the model in column \eqn{j}; \code{NA} means the pair could
not be checked, or the result was inconclusive), \code{fx} (the
corresponding fit function minima), \code{df} (the degrees of freedom of
the models, in the order used in the matrix), \code{crit} and
\code{npoints}.
}
\references{
Bentler, P. M., & Satorra, A. (2010). Testing model nesting and equivalence.
\emph{Psychological Methods, 15}(2), 111-123. \doi{10.1037/a0019625}

Asparouhov, T., & \enc{Muthén}{Muthen}, B. (2019). Nesting and equivalence
testing for structural equation models. \emph{Structural Equation Modeling:
A Multidisciplinary Journal, 26}(2), 302-309.
\doi{10.1080/10705511.2018.1513795}
}
\seealso{
\code{\link{lavTestLRT}}, which uses this machinery (see its
\code{check_nested} argument) to warn when a chi-squared difference test
is requested for non-nested models.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
# nested: the same model with orthogonal factors
fit0 <- cfa(HS.model, data = HolzingerSwineford1939, orthogonal = TRUE)
# equivalent: the same model with a different scaling of the factors
fit1b <- cfa(HS.model, data = HolzingerSwineford1939, std.lv = TRUE)
# not nested
model2 <- ' f1 =~ x1 + x2 + x3 + x4
            f2 =~ x5 + x6 + x7 + x8 + x9 '
fit2 <- cfa(model2, data = HolzingerSwineford1939)

lavTestNET(fit1, fit0, fit1b, fit2,
           model_names = c("3f", "3f-orth", "3f-std.lv", "2f"))
}
