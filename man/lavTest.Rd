\name{lavTest}
\alias{lavTest}
\alias{lavtest}
\title{Test of exact fit}
\description{
Compute a variety of test statistics evaluating the global fit of
the model.}
\usage{
    lavTest(lavobject, test = "standard", scaled_test = "standard",
            output = "list", drop_list_single = TRUE, ...)
}
\arguments{
\item{lavobject}{An object of class \code{\linkS4class{lavaan}}.}
\item{test}{Character vector. Multiple names of test statistics can be provided.
    If \code{"standard"} is included, a conventional chi-square test
    is computed. If \code{"Browne.residual.adf"} is included,
    Browne's residual-based test statistic using ADF theory is computed.
    If \code{"Browne.residual.nt"} is included, Browne's residual-based
    test statistic using normal theory is computed.
    If \code{"Satorra.Bentler"} is included, a Satorra-Bentler scaled test
    statistic is computed. If \code{"Yuan.Bentler"} is included,
    a Yuan-Bentler scaled test statistic is computed.
    If \code{"Yuan.Bentler.Mplus"} is included, a
    test statistic is computed that is asymptotically equivalent to the
    Yuan-Bentler scaled test statistic. If \code{"mean.var.adjusted"} or
    \code{"Satterthwaite"} is included, a mean and variance adjusted test
    statistic is computed. If \code{"scaled.shifted"} is included,
    an alternative mean and variance adjusted test statistic is
    computed (as in Mplus version 6 or higher).
    If \code{"mean.var.adjusted.corrected"} (aliases: \code{"mvc"},
    \code{"hayakawa"}) or \code{"scaled.shifted.corrected"} (alias:
    \code{"ssc"}) is included, the corresponding test statistic is computed
    with the corrected (unbiased) estimator of the trace of the squared
    UGamma matrix (Hayakawa, 2018), which remains accurate when the number
    of observed variables is large relative to the sample size (where the
    naive estimator is severely biased, and the uncorrected adjusted tests
    underreject). These corrected tests are only available for single-group,
    single-level, complete, continuous data and \code{estimator = "ML"}.
    Following Hayakawa (2018), they are best combined with
    \code{scaled_test = "RLS"}.
    If \code{"boot"} or \code{"bootstrap"} or \code{"Bollen.Stine"} is
    included, the
    Bollen-Stine bootstrap is used to compute the bootstrap probability value
    of the (regular) test statistic. FMG p-value tests can be requested using
    names such as \code{"peba4"}, \code{"pols3"}, \code{"pall"}, or
    \code{"all"}. The convenience value \code{"fmg"} resolves to
    \code{"peba4_rls"} for one-model tests. Suffixless FMG names use the
    standard ML chi-square statistic by default; set
    \code{scaled_test = "Browne.residual.nt.model"} or
    \code{scaled_test = "RLS"} to use the Browne residual normal-theory
    model-based statistic. Suffixes such as \code{"_ml"}, \code{"_rls"},
    and \code{"_ug"} can be used for explicit control, for example
    \code{"peba4_ug_rls"} or \code{"pols3_ml"}.}
\item{scaled_test}{Character. Choose the test statistic
    that will be scaled (if a scaled test statistic is requested).
    The default is \code{"standard"}, but it could also be (for example)
    \code{"Browne.residual.nt"}.}
\item{output}{Character. If \code{"list"} (the default), return a list with
    all test statistics. If \code{"text"}, display the output as text with
    verbose descriptions (as in the summary output). If any scaled
    test statistics are included, they are printed first in a two-column
    format, followed by the remaining test statistics in a one-column format.}
\item{drop_list_single}{Logical. Only used when \code{output = "list"}.
    If \code{TRUE} and the list is of length one (i.e. only a single test
    statistic), drop the outer list. If \code{FALSE}, return a nested list
    with one element per test statistic.}
\item{...}{To support old argument names.}
}
\value{
    If \code{output = "list"}: a nested list with test statistics, or if
    only a single test statistic is requested (and
    \code{drop_list_single = TRUE}), a list with details for this test
    statistic. If \code{output = "text"}: the text is printed, and a
    nested list of test statistics (including an info attribute) is
    returned.
}
\references{
Hayakawa, K. (2019). Corrected goodness-of-fit test in covariance
structure analysis. \emph{Psychological Methods, 24}(3), 371--389.
\doi{10.1037/met0000180}

Himeno, T., & Yamada, T. (2014). Estimations for some functions of
covariance matrix in high dimension under non-normality and its
applications. \emph{Journal of Multivariate Analysis, 130}, 27--44.

Srivastava, M. S. (2005). Some tests concerning the covariance matrix in
high dimensional data. \emph{Journal of the Japan Statistical Society,
35}(2), 251--272.
}
\examples{
HS.model <- '
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
'
fit <- cfa(HS.model, data = HolzingerSwineford1939)
lavTest(fit, test = "browne.residual.adf")
lavTest(fit, test = "peba4", output = "text")
lavTest(fit, test = "pols3", scaled_test = "RLS", output = "text")
# corrected adjusted test (Hayakawa, 2018), based on the RLS statistic
lavTest(fit, test = "mean.var.adjusted.corrected", scaled_test = "RLS",
        output = "text")
}
