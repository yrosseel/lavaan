\name{modificationIndices}
\alias{modificationIndices}
\alias{modificationindices}
\alias{modindices}
\title{Modification Indices}
\description{
Given a fitted lavaan object, compute the modification indices
(= univariate score tests) for a selected set of fixed-to-zero parameters.
}
\usage{
modificationIndices(object, standardized = TRUE, cov.std = TRUE,
                    information = "expected",
                    power = FALSE, delta = 0.1, alpha = 0.05,
                    high.power = 0.75, sort. = FALSE, minimum.value = 0,
                    maximum.number = nrow(LIST), free.remove = TRUE,
                    na.remove = TRUE, op = NULL)
modindices(object, standardized = TRUE, cov.std = TRUE, information = "expected",
                    power = FALSE, delta = 0.1, alpha = 0.05, high.power = 0.75,
                    sort. = FALSE, minimum.value = 0,
                    maximum.number = nrow(LIST), free.remove = TRUE,
                    na.remove = TRUE, op = NULL)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{standardized}{If \code{TRUE}, two extra columns (sepc.lv and sepc.all)
will contain standardized values for the epc's. In the first column (sepc.lv),
standardizization is based on the variances of the (continuous) latent
variables. In the second column (sepc.all), standardization is based
on both the variances of both (continuous) observed and latent variables.
(Residual) covariances are standardized using (residual) variances.}
\item{cov.std}{Logical. See \code{\link{standardizedSolution}}.}
\item{information}{\code{character} indicating the type of information matrix
to use (check \code{\link{lavInspect}} for available options). \code{"expected"}
information is the default, which provides better control of Type I errors.}
\item{power}{If \code{TRUE}, the (post-hoc) power is computed for each
modification index, using the values of \code{delta} and \code{alpha}.}
\item{delta}{The value of the effect size, as used in the post-hoc power
computation, currently using the unstandardized metric of the epc column.}
\item{alpha}{The significance level used for deciding if the modification
index is statistically significant or not.}
\item{high.power}{If the computed power is higher than this cutoff value,
the power is considered `high'. If not, the power is considered `low'.
This affects the values in the 'decision' column in the output.}
\item{sort.}{Logical. If TRUE, sort the output using the values of
the modification index values. Higher values appear first.}
\item{minimum.value}{Numeric. Filter output and only show rows with a
modification index value equal or higher than this minimum value.}
\item{maximum.number}{Integer. Filter output and only show the first
maximum number rows. Most useful when combined with the \code{sort.} option.}
\item{free.remove}{Logical. If TRUE, filter output by removing all rows
corresponding to free (unconstrained) parameters in the original model.}
\item{na.remove}{Logical. If TRUE, filter output by removing all rows with
NA values for the modification indices.}
\item{op}{Character string. Filter the output by selecting only those rows with
operator \code{op}.}
}
\value{
  A data.frame containing modification indices and EPC's.
}
\details{
  Modification indices are just 1-df (or univariate) score tests. The
  modification index (or score test) for a single parameter reflects
  (approximately) the improvement in model fit (in terms of the chi-square
  test statistic), if we would refit the model but allow this parameter to
  be free.
  This function is a convenience function in the sense that it produces a
  (hopefully sensible) table of currently fixed-to-zero (of fixed to another
  constant) parameters. For each of these parameters, a modification index
  is computed, together with an expected parameter change (epc) value.
  It is important to realize that this function will only consider
  fixed-to-zero parameters. If you have equality constraints in the model,
  and you wish to examine what happens if you release all (or some) of these
  equality constraints, use the \code{\link{lavTestScore}} function.
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
modindices(fit, minimum.value = 10, sort = TRUE)
}
