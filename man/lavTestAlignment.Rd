\name{lavTestAlignment}
\alias{lavTestAlignment}
\title{Invariance analysis of aligned measurement parameters}
\description{
Post-estimation invariance analysis of the measurement parameters of a
multiple-group model with efa blocks that were rotated with the agreement
rotation (\code{rotation_args = list(mg_agreement = TRUE)}), of a
confirmatory model aligned with \code{alignment = "fixed"} or
\code{"free"} (see \code{\link{alignment}}), or of a model with both,
following the
procedure of the alignment method (Asparouhov and Muthen, 2014). For every
aligned loading, and for every aligned intercept (or threshold, for an
ordered indicator) if the intercepts were aligned
(\code{mg_agreement_intercepts}), the groups are divided into an
(approximately) invariant set and the non-invariant groups.}
\usage{
lavTestAlignment(object, p_connect = 0.01, p_set = 0.001, max_iter = 100L)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}: a
    multiple-group model with efa blocks, rotated with
    \code{rotation_args = list(mg_agreement = TRUE)}, or a confirmatory
    model with \code{alignment = "fixed"} or \code{"free"}, with standard
    errors.}
\item{p_connect}{Numeric. Two groups are \sQuote{connected} when the
    p-value of the pairwise test of equality of their parameter exceeds this
    value. The largest connected set of groups is the starting invariant
    set.}
\item{p_set}{Numeric. In the iterations, a group is kept in (or added to)
    the invariant set when the p-value of the test of its parameter against
    the average over the invariant set exceeds this value, and removed
    otherwise.}
\item{max_iter}{Integer. The maximum number of iterations per parameter.}
}
\details{
For every parameter, the procedure is the following. (1) The equality of
the parameter is tested for every pair of groups with a z test, using the
(delta method) covariance matrix of the aligned estimates. Two groups are
connected when the p-value exceeds \code{p_connect}; the starting invariant
set is the largest connected set of groups. (2) The average of the parameter
over the invariant set is computed, and every group is tested against that
average. The groups with a p-value above \code{p_set} form the new invariant
set. This is repeated until the set no longer changes (or a cycle is
detected). (3) The average over the final invariant set, the non-invariant
groups, the contribution of the parameter to the agreement (loss) function
of the rotation, and an R-square are reported. The R-square is the
proportion of the across-group variance of the unaligned parameter (the
loading with a unit-variance factor, or the intercept (threshold) with a zero factor
mean) that is explained by the group differences in factor variances and
means (Asparouhov and Muthen, 2014, eqs. 13 and 14); it is only computed
with three or more groups. A value close to one means that the group
differences in the unaligned parameter are absorbed by the factor variances
and means, a value close to zero that they are not.
}
\value{
A data.frame of class \code{lavaan.data.frame} with one row per parameter:
the parameter (\code{lhs}, \code{op}, \code{rhs}), its \code{average} over
the invariant set, the \code{invariant} set and the \code{noninvariant}
groups (group labels), the R-square (\code{r2}) and the
\code{contribution} of the parameter to the loss function. The attributes
\code{pvalue} (the p-value of every group against the average of the final
invariant set) and \code{invariant} (logical) hold the parameter by group
matrices.
}
\references{
Asparouhov, T., & Muthen, B. (2014). Multiple-group factor analysis
alignment. \emph{Structural Equation Modeling, 21}(4), 495--508.
}
\seealso{
\code{\link{alignment}} for the alignment of a confirmatory model;
\code{\link{esem}} for the rotation options of the agreement rotation,
including the alignment of the loadings and the intercepts of efa blocks.
}
\examples{
\donttest{
HS <- HolzingerSwineford1939
model <- ' efa("b")*f1 + efa("b")*f2 + efa("b")*f3 =~
             x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
# a (partially specified) target towards the confirmatory pattern
target <- matrix(NA, 9, 3)
target[4:9, 1] <- 0; target[c(1:3, 7:9), 2] <- 0; target[1:6, 3] <- 0
# alignment of the loadings and the intercepts
fit <- cfa(model, data = HS, group = "school", meanstructure = TRUE,
           rotation = "target",
           rotation_args = list(target = target, mg_agreement = TRUE,
             mg_agreement_method = "loading.alignment",
             mg_agreement_weight = 1, mg_agreement_scales = "free",
             mg_agreement_rotate = FALSE, mg_agreement_intercepts = "fixed",
             mg_agreement_group_weights = "sqrt.nobs"))
lavTestAlignment(fit)
}
}
