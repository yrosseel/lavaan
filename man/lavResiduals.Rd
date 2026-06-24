\name{lavResiduals}
\alias{lavResiduals}
\alias{lavResidual}
\title{Residuals}
\description{
\sQuote{lavResiduals} provides model residuals and standardized residuals from
a fitted lavaan object, as well as various summaries of these residuals.

The \sQuote{residuals()} (and \sQuote{resid()}) methods are just shortcuts to
this function with a limited set of arguments.
}
\usage{
lavResiduals(object, type = "cor.bentler", h1 = NULL,
    se = FALSE, zstat = TRUE, summary = TRUE,
    summary.blocks.combined = FALSE, h1.acov = "unstructured",
    add.type = TRUE, add.labels = TRUE, add.class = TRUE,
    drop.list.single.group = TRUE,
    maximum.number = 0L, output = "list")
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{Character.
If \code{type = "raw"}, this function returns the raw (= unscaled)
difference between the observed and the expected (model-implied) summary
statistics, as well as the standardized version of these residuals.
If \code{type = "cor"}, or \code{type = "cor.bollen"}, the observed and model
implied covariance matrices are first transformed to a correlation matrix
(using \code{cov2cor()}), before the residuals are computed.
If \code{type = "cor.bentler"}, both the observed and model implied covariance
matrices are rescaled by dividing the elements by the square roots of the
corresponding variances of the observed covariance matrix.}
\item{h1}{Optional. A user-provided saturated (unrestricted) model supplying
the observed summary statistics against which the model-implied statistics are
compared. It can be a fitted \code{\linkS4class{lavaan}} object (its \code{h1}
slot is used), a \code{lavh1} list (a list with an \code{implied} element), or
an implied-moments list (a list with a \code{cov} or \code{res.cov} element,
per block). It must be compatible with \code{object} (same number of blocks
and same dimensions). If \code{NULL} (the default), the saturated statistics of
\code{object} itself are used.}
\item{se}{Logical. If \code{TRUE}, show the estimated standard errors
for the residuals.}
\item{zstat}{Logical. If \code{TRUE}, show the standardized residuals, which
are the raw residuals divided by the corresponding (estimated) standard
errors.}
\item{summary}{Logical. If \code{TRUE}, show various summaries of the
(possibly scaled) residuals. When \code{type = "raw"}, we compute the
RMR. When \code{type = "cor.bentler"},
we compute the SRMR. When \code{type = "cor.bollen"}, we compute the CRMR.
An unbiased version of these summaries (the URMR, USRMR and UCRMR) is also
computed (Maydeu-Olivares, 2017).

The two versions are reported with different inferential statistics, and this
asymmetry is intentional. The sample RMR/SRMR/CRMR overestimate (are positively
biased for) their population values. For these (biased) statistics, we report a
standard error and a test of \emph{exact} fit (the null hypothesis that the
population value equals zero). This test
relies on the known sampling distribution of the (biased) sample statistic
under exact fit (Maydeu-Olivares, 2017, Eqs. 28--30), where the bias is
explicitly accounted for. For the (approximately) unbiased URMR/USRMR/UCRMR, we
report a standard error, a confidence interval, and a test of \emph{close} fit
(the null hypothesis that the population value equals a small nonzero cutoff,
e.g. 0.05). The confidence interval and the close-fit test require an
(approximately) unbiased and normally distributed estimator, which is precisely
why the unbiased version is used here. Conversely, no test of exact fit is
reported for the unbiased estimator: because it is truncated at zero, its
sampling distribution is degenerate at the exact-fit boundary, so the test of
exact fit is (correctly) only based on the biased statistic.}
\item{summary.blocks.combined}{Logical. Only relevant when multiple blocks are
involved (multiple groups, or multiple levels). If \code{FALSE} (the default),
a separate summary table is shown for each block. If \code{TRUE}, the
block-specific summary tables are replaced by a single overall summary table in
which the residual elements are pooled across all blocks (so, e.g., the SRMR is
computed from the residuals of all groups or all levels at once). The combined
table is returned as the \code{summary} element at the top level of the result.
The pooling correctly accounts for the (block-diagonal) independence of separate
groups and for the (joint) sampling covariance across the levels of a multilevel
model.}
\item{h1.acov}{Character. If \code{"unstructured"}, the observed summary
statistics are used as consistent estimates of the corresponding (unrestricted)
population statistics. If \code{"structured"}, the model-implied summary
statistics are used as consistent estimates of the corresponding (unrestricted)
population statistics. This affects the way the asymptotic variance matrix
of the summary statistics is computed.}
\item{add.type}{Logical. If \code{TRUE}, show the type of residuals
in the output.}
\item{add.labels}{If \code{TRUE}, variable names are added to the vectors
and/or matrices.}
\item{add.class}{If \code{TRUE}, vectors are given the \sQuote{lavaan.vector}
class; matrices are given the \sQuote{lavaan.matrix} class, and symmetric
matrices are given the \sQuote{lavaan.matrix.symmetric} class.
This only affects the way they are printed on the screen.}
\item{drop.list.single.group}{If \code{FALSE}, the results are returned as
a list, where each element corresponds to a group (even if there is only
a single group). If \code{TRUE}, the list will be unlisted if there is
only a single group.}
\item{maximum.number}{Integer. Only used if \code{output = "table"}. If larger
than zero, show only the first maximum.number rows of the data.frame. The
default (\code{0L}) shows all rows.}
\item{output}{Character. By default, \code{output = "list"}, and the output
is a list of elements. If \code{output = "table"}, the residuals are shown in
a single data.frame, one row per residual element, sorted from high (in
absolute value) to low. All residual elements are included: the (co)variances
(\code{lhs~~rhs}), the means/intercepts (\code{lhs~1}), the thresholds
(\code{lhs|t1}) and, when \code{conditional.x = TRUE}, the regression slopes
(\code{lhs~rhs}). Variances (the diagonal) are only shown for the raw-style
metrics (\code{type = "raw"}, \code{"normalized"} or \code{"standardized"}),
not for the correlation metrics where they are structurally zero. Structurally
fixed moments (e.g. saturated means or just-identified thresholds) are omitted.
If \code{se = TRUE}, a \code{se} column with the standard errors is added; if
\code{zstat = TRUE}, a \code{z} column with the z-statistics is added.}
}
\value{
If \code{drop.list.single.group = TRUE}, a list of (residualized) summary
statistics, including type, standardized residuals, and summaries. If
\code{drop.list.single.group = FALSE}, the list of summary statistics is nested
within a list for each group.
}
\references{
Bentler, P.M. and Dijkstra, T. (1985). Efficient estimation via linearization
in structural models. In Krishnaiah, P.R. (Ed.),
\emph{Multivariate analysis - VI}, (pp. 9--42). New York, NY: Elsevier.

Ogasawara, H. (2001). Standard errors of fit indices using residuals in
structural equation modeling. \emph{Psychometrika, 66}(3), 421--436.
doi:10.1007/BF02294443

Maydeu-Olivares, A. (2017). Assessing the size of model misfit in structural
equation models. \emph{Psychometrika, 82}(3), 533--558.
doi:10.1007/s11336-016-9552-7

Standardized Residuals in M\emph{plus}. Document retrieved from URL
http://www.statmodel.com/download/StandardizedResiduals.pdf
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
lavResiduals(fit)
}
