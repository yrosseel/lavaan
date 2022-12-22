\name{lavTest}
\alias{lavTest}
\alias{lavtest}
\title{Test of exact fit}
\description{
Compute a variety of test statistics evaluating the global fit of
the model.}
\usage{
    lavTest(lavobject, test = "standard", scaled.test = "standard",
            output = "list", drop.list.single = TRUE)
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
    test statistic is computed that is asymptotically equal to the
    Yuan-Bentler scaled test statistic. If \code{"mean.var.adjusted"} or
    \code{"Satterthwaite"} is included, a mean and variance adjusted test 
    statistic is computed. If \code{"scaled.shifted"} is included, 
    an alternative mean and variance adjusted test statistic is 
    computed (as in Mplus version 6 or higher).
    If \code{"boot"} or \code{"bootstrap"} or \code{"Bollen.Stine"} is 
    included, the
    Bollen-Stine bootstrap is used to compute the bootstrap probability value
    of the (regular) test statistic.}
\item{scaled.test}{Character. Choose the test statistic 
    that will be scaled (if a scaled test statistic is requested). 
    The default is \code{"standard"}, but it could also be (for example)
    \code{"Browne.residual.nt"}.}
\item{output}{Character. If \code{"list"} (the default), return a list with
    all test statistics. If \code{"text"}, display the output as text with 
    verbose descriptions (as in the summary output). If any scaled
    test statistics are included, they are printed first in a two-column
    format. Next come the other test statistics in a one-column format.}
\item{drop.list.single}{Logical. Only used when \code{output = "list"}.
    If \code{TRUE} and the list is of length one (i.e. only a single test 
    statistic), drop the outer list. If \code{FALSE}, return a nested list 
    with as many elements as we have test statistics.}
}
\value{
    If \code{output = "list"}: a nested list with test statistics, or if
    only a single test statistic is requested (and 
    \code{drop.list.single = TRUE}), a list with details for this test
    statistic. If \code{output = "text"}: the text is printed, and a
    nested list of test statistics (including an info attribute) is
    returned.
}
\examples{
HS.model <- '
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9
'
fit <- cfa(HS.model, data = HolzingerSwineford1939)
lavTest(fit, test = "browne.residual.adf")
}
