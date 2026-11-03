\name{alignment}
\alias{alignment}
\alias{alignment.args}
\title{Alignment of a multiple-group confirmatory model}
\description{
The \code{alignment} option of \code{\link{lavOptions}} aligns a
multiple-group confirmatory factor model: the factor variances and means of
the configural model (with \code{std.lv = TRUE}: factor variances one and
factor means zero in every group) are chosen to minimize the amount of
measurement non-invariance of the loadings and the intercepts, following
the alignment method of Asparouhov and Muthen (2014). The aligned model has
the same fit as the configural model; it provides comparable group-specific
factor variances and means without assuming measurement invariance, and the
aligned loadings and intercepts can be inspected for (non-)invariance with
\code{\link{lavTestAlignment}}.}
\section{Usage}{
\preformatted{
fit <- cfa(model, data = data, group = "g", std.lv = TRUE,
           meanstructure = TRUE, alignment = "fixed",
           alignment.args = list(method = "loading.alignment",
                                 epsilon = 0.01,
                                 group_weights = "sqrt.nobs",
                                 identification = "reference",
                                 std_ov = FALSE, rstarts = 30L,
                                 reference_group = 1L))
}
}
\section{Options}{
\describe{
  \item{\code{alignment}}{\code{"none"} (the default), \code{"fixed"} or
    \code{"free"}. With \code{"fixed"} the factor means of the first group
    are fixed to zero; with \code{"free"} the means of all groups are
    aligned, which requires more than two groups and some loading
    non-invariance to be identified. The model must be a multiple-group
    configural model: latent variables with indicators, \code{std.lv =
    TRUE}, a mean structure, free loadings and intercepts in every group,
    no equality constraints. The model may also contain efa blocks (see
    \code{\link{esem}}): they are rotated first, with the agreement
    rotation if \code{rotation_args = list(mg_agreement = TRUE)}, and the
    confirmatory factors are then aligned on the rotated model; the two
    transformations act on different factors and compose.}
  \item{\code{alignment.args}}{A list with the settings of the alignment
    loss (snake_case names):
    \code{method} (\code{"loading.alignment"}, the default: the component
    loss \eqn{(d^2 + \epsilon)^{1/4}} summed over all pairs of groups and
    all loadings and intercepts, with \eqn{d} the difference between the
    two parameters; \code{"relative.alignment"}: the same loss on
    relative differences, a loading entering as \eqn{\log(\lambda^2 +
    \epsilon)/2}{log(lambda^2 + epsilon)/2} and an intercept difference divided by the scale of
    the indicator, see \code{mg_agreement_method} in \code{\link{esem}};
    this is the loss that reproduces the alignment of Mplus; or
    \code{"procrustes"}: squared
    differences), \code{epsilon} (0.01), \code{group_weights}
    (\code{"sqrt.nobs"}, the default: a pair of groups is weighted by the
    square root of the product of their sample sizes; \code{"nobs"} or
    \code{"equal"}), \code{identification} (\code{"reference"}, the
    default: the solution of the product identification re-expressed with
    variance one in the first group; \code{"product"}: the product of the
    factor variances is one; \code{"mean"}: their average is one),
    \code{std_ov} (\code{FALSE}; if \code{TRUE}, the loss is evaluated on
    the loadings and intercepts of the indicators standardized over the
    entire population, as the alignment method does), \code{rstarts}
    (30, the number of random starts of the loss minimization in addition
    to the start at equal column sums of squares) and
    \code{reference_group} (1, the group with variance one under the
    \code{"reference"} identification and mean zero under \code{alignment
    = "fixed"}: a group number, a group label, or \code{"auto"}, which
    tries every group as the reference of a fixed alignment and keeps the
    one with the smallest loss, i.e. the most invariant aligned model; with
    \code{alignment = "free"} the loss does not depend on the reference
    group and \code{"auto"} keeps the first group).}
}
}
\section{Details}{
The aligned loadings are the configural loadings divided by the square
roots of the aligned factor variances, the aligned intercepts are the
configural intercepts minus the aligned loadings times the aligned factor
means; all other parameters are transformed so that the model-implied
moments are unchanged (the factor covariances, the regressions among the
factors, the factor intercepts). The aligned factor variances and means
(except those fixed by the identification) are parameters of the aligned
model, with standard errors obtained by the delta method from the
configural model, or by the bootstrap (\code{se = "bootstrap"}: every
bootstrap sample is refit and aligned with the same settings, and the
variability of the aligned parameters is used).
\code{\link{lavInspect}} with \code{"rotation"} returns the alignment
information per group (the aligned variances and means, the loss value,
convergence, the reference group, and the loss value of every random
start in \code{starts}), and \code{"rotation_matrix"} the diagonal
rescaling.

Two diagnostics are reported as warnings. The loss function has several
local optima when the best solution was reached by a single random
start; increase \code{rstarts} and compare the solutions. A free
alignment of the means (\code{alignment = "free"}) may be poorly
identified: the overall level of the means of a factor is identified by
the loss only, and the warning is given when the standard error of the
average of the group means exceeds the largest standard error of a
difference between two group means (which is identified by the data).
This happens with few groups or little non-invariance, as Asparouhov and
Muthen (2014) describe; use the fixed alignment then.
For ordered indicators (\code{ordered}), the thresholds take the place of
the intercepts: every threshold of an indicator is a component of the
alignment loss, together with the loading of the indicator, and the
aligned thresholds are the configural thresholds plus the aligned loadings
times the aligned factor means (the intercepts of ordered indicators stay
fixed to zero); the model-implied thresholds are unchanged. Continuous and
ordered indicators can be mixed. Since the loss is not invariant to the
scale of the indicators, the solution depends on the
\code{parameterization} (\code{"delta"} or \code{"theta"}) of the
configural model.

The model may contain exogenous covariates (\code{f ~ x}), with
\code{conditional.x = TRUE} or \code{FALSE}: the regressions of the
factors on the covariates transform with the factor scales, and the
aligned factor means are the intercepts of the factors given the
covariates (as in Mplus; center the covariates to keep them
interpretable as means). \code{std_ov = TRUE} is not available with
\code{conditional.x = TRUE}. Alignment is not available for multilevel
models, or when an indicator of an aligned factor is itself a predictor
or an outcome in a regression. In a model with efa blocks, the
efa factors are not aligned by this option but by the agreement rotation
(\code{mg_agreement_scales}, \code{mg_agreement_intercepts}); the
rotation information of \code{\link{lavInspect}} then holds both the
efa sets and the \code{alignment} element.
}
\references{
Asparouhov, T., & Muthen, B. (2014). Multiple-group factor analysis
alignment. \emph{Structural Equation Modeling, 21}(4), 495--508.

Asparouhov, T., & Muthen, B. (2023). Multiple group alignment for
exploratory and structural equation models. \emph{Structural Equation
Modeling, 30}(2), 169--191.
}
\seealso{
\code{\link{lavTestAlignment}} for the invariance analysis of the aligned
parameters (efa and confirmatory); \code{\link{esem}} for the alignment
of efa blocks (\code{rotation_args = list(mg_agreement = TRUE, ...)}).
}
\examples{
\donttest{
HS <- HolzingerSwineford1939
model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '
fit <- cfa(model, data = HS, group = "school", std.lv = TRUE,
           meanstructure = TRUE, alignment = "fixed")
lavInspect(fit, "rotation")
lavTestAlignment(fit)
}
}
