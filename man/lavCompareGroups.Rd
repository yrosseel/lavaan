\name{lavCompareGroups}
\alias{lavCompareGroups}
\title{Compare measurement parameters across groups}
\description{
Side-by-side tables of the measurement parameters of a multiple-group model:
for every factor, its loadings in all groups (indicators in the rows, groups
in the columns), or the intercepts and thresholds of the indicators in all
groups. The range of every parameter over the groups is added, so that
(possibly) non-invariant items stand out. The function is meant for rotated
solutions (ESEM, with or without the agreement rotation, see
\code{\link{esem}}) and aligned models (see \code{\link{alignment}}), but
works for any multiple-group model.}
\usage{
lavCompareGroups(object, what = "loadings", standardized = FALSE,
                 add_range = TRUE, sort = FALSE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}, fitted to
    multiple groups.}
\item{what}{\code{"loadings"} (the default) or \code{"intercepts"}: the
    factor loadings (one table per factor), or the intercepts of the
    indicators together with the thresholds of the ordered indicators
    (one table).}
\item{standardized}{Logical. If \code{FALSE} (the default), the
    estimates in the metric of the model (for a rotated or aligned
    solution: the rotated, aligned loadings). If \code{TRUE}, the completely
    standardized estimates (\code{std.all}).}
\item{add_range}{Logical. If \code{TRUE} (the default), a column
    \code{range} with the difference between the largest and the smallest
    value over the groups is added to every table.}
\item{sort}{Logical. If \code{TRUE}, the rows of every table are sorted by
    decreasing range.}
}
\value{
For \code{what = "loadings"}, a named list with one data frame per factor
(the indicators in the rows, the groups in the columns); for \code{what =
"intercepts"}, a single data frame (the indicators, or the thresholds of
the ordered indicators labelled \code{item|threshold}, in the rows). The
data frames are of class \code{lavaan.data.frame}.
}
\seealso{
\code{\link{lavTestAlignment}} for the invariance tests of an aligned
model; \code{\link{esem}} and \code{\link{alignment}}.
}
\examples{
HS <- HolzingerSwineford1939
model <- ' efa("b")*f1 + efa("b")*f2 + efa("b")*f3 =~
             x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
# the configural (per-group) rotation
fit <- cfa(model, data = HS, group = "school")
lavCompareGroups(fit)
# the agreement rotation: the solutions of the groups are made as similar
# as possible; what remains is (possibly) non-invariance
fit_a <- cfa(model, data = HS, group = "school",
             rotation_args = list(mg_agreement = TRUE))
lavCompareGroups(fit_a, sort = TRUE)
lavCompareGroups(fit_a, standardized = TRUE)
}
