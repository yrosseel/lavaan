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
                                 std_ov = FALSE, rstarts = 30L))
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
    no efa blocks, no equality constraints.}
  \item{\code{alignment.args}}{A list with the settings of the alignment
    loss (snake_case names):
    \code{method} (\code{"loading.alignment"}, the default: the component
    loss \eqn{(d^2 + \epsilon)^{1/4}} summed over all pairs of groups and
    all loadings and intercepts; or \code{"procrustes"}: squared
    differences), \code{epsilon} (0.01), \code{group_weights}
    (\code{"sqrt.nobs"}, the default: a pair of groups is weighted by the
    square root of the product of their sample sizes; \code{"nobs"} or
    \code{"equal"}), \code{identification} (\code{"reference"}, the
    default: the solution of the product identification re-expressed with
    variance one in the first group; \code{"product"}: the product of the
    factor variances is one; \code{"mean"}: their average is one),
    \code{std_ov} (\code{FALSE}; if \code{TRUE}, the loss is evaluated on
    the loadings and intercepts of the indicators standardized over the
    entire population, as the alignment method does) and \code{rstarts}
    (30, the number of random starts of the loss minimization in addition
    to the start at equal column sums of squares).}
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
configural model; \code{se = "bootstrap"} is not available.
\code{\link{lavInspect}} with \code{"rotation"} returns the alignment
information per group (the aligned variances and means, the loss value,
convergence), and \code{"rotation_matrix"} the diagonal rescaling.
Alignment is not available for categorical data, \code{conditional.x =
TRUE}, multilevel models, or models in which observed variables appear as
predictors or outcomes in regressions.
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
parameters; \code{\link{esem}} for the alignment of efa blocks
(\code{rotation_args = list(mg_agreement = TRUE, ...)}).
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
