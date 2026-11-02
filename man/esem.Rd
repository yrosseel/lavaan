\name{esem}
\alias{esem}
\alias{ESEM}
\alias{set-ESEM}
\alias{rotation}
\title{Exploratory Structural Equation Modeling}
\description{
Exploratory Structural Equation Modeling (ESEM) combines the flexibility of
exploratory factor analysis (EFA) with the modeling options of structural
equation modeling (SEM). In lavaan, this is done by marking one or more
\emph{sets} of latent variables in the model syntax with the \code{efa()}
modifier. The factors within a set form a block of exploratory factors: all
their factor loadings are (in principle) freely estimated, and the block is
rotated after the model has been fitted. The rest of the model (confirmatory
factors, regressions, covariates, residual covariances, multiple groups,
multiple levels, \ldots) is specified as usual.

This page explains how such sets of factors are defined, what lavaan does to
identify them, gives a brief overview of the available rotation methods and
options, and shows how to use a different target matrix per set in a target
rotation.
}
\details{
\subsection{Defining a set of exploratory factors}{
A set is defined by premultiplying the latent variable names (on the left-hand
side of the \code{=~} operator) by the \code{efa("label")} modifier, where
\code{"label"} is an arbitrary name for the set. All latent variables that
carry the same label belong to the same set:

\preformatted{
  efa("block1")*f1 +
  efa("block1")*f2 +
  efa("block1")*f3 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
}

Each factor in the set should load on the same indicators (here, all of
\code{x1} to \code{x9}). Using different labels defines different sets, for
example one set per measurement occasion or per questionnaire:

\preformatted{
  efa("t1")*f1 + efa("t1")*f2 =~ y1 + y2 + y3 + y4 + y5 + y6
  efa("t2")*f3 + efa("t2")*f4 =~ y7 + y8 + y9 + y10 + y11 + y12
}

Each set is rotated \emph{independently} of the other sets; this is sometimes
called set-ESEM. A set with a single factor is not rotated. Sets can be freely
combined with confirmatory (CFA) factors, with regressions among the factors
(exploratory or not), with observed covariates, and with residual covariances.
For a single set of factors without any further structure, the
\code{\link{efa}} function provides a convenient wrapper.
}

\subsection{Identification of the unrotated solution}{
Within a set, the factors are not identified without further constraints.
By default (\code{auto_efa = TRUE}; see \code{\link{lavOptions}}), lavaan
imposes the constraints that are needed to identify the \emph{unrotated}
solution: the variances of the factors in the set are fixed to 1, their
covariances are fixed to zero, and the factor loadings follow an echelon
pattern (the loading of the second factor on the first indicator is fixed to
zero, the loadings of the third factor on the first two indicators are fixed
to zero, and so on). The unrotated solution can be inspected by setting
\code{rotation = "none"}.

After the model has been fitted, the loadings of each set are rotated, and the
remaining parameters that depend on the factors (factor covariances,
regression coefficients, intercepts) are transformed accordingly, so that the
reported solution is the rotated one. For oblique rotations, the factor
covariances within a set are reported as free parameters. Standard errors for
the rotated solution are obtained by one of the methods described below.

The structural part of the model can be specified in terms of the (rotated)
factors in the usual way, for example
\code{f3 + f4 ~ f1 + f2} or \code{f1 + f2 + f3 ~ age + sex}. Because the
rotation is chosen by the rotation criterion, and not by the structural part
of the model, the rotation does not depend on the regressions.
}

\subsection{Rotation methods}{
The rotation method is selected by the \code{rotation} argument (the default
is \code{"geomin"}). It applies to all sets in the model. By default, the
rotation is oblique; an orthogonal rotation is requested by setting the
\code{orthogonal} element of \code{rotation_args} to \code{TRUE}. The
following criteria are available; most of them are available in both an
orthogonal and an oblique version:
\itemize{
  \item \code{"varimax"}, \code{"quartimax"}, \code{"orthomax"}
    (with \code{orthomax_gamma}), \code{"entropy"}, \code{"mccammon"},
    \code{"infomax"}, \code{"tandem1"}, \code{"tandem2"}: classic
    criteria that are typically used for orthogonal rotation.
  \item \code{"oblimin"} (with \code{oblimin_gamma}), \code{"quartimin"},
    \code{"geomin"} (with \code{geomin_epsilon}), \code{"oblimax"},
    \code{"simplimax"} (making the \eqn{p(m-1)} smallest loadings as
    close to zero as possible, as in GPArotation 2026.8 and later),
    \code{"bentler"}: criteria that are typically used for oblique
    rotation.
  \item \code{"cf"} (Crawford-Ferguson, with \code{cf_gamma}), and its
    special cases \code{"cf-quartimax"}, \code{"cf-varimax"},
    \code{"cf-equamax"}, \code{"cf-parsimax"} and \code{"cf-facparsim"}.
  \item \code{"target.strict"}, \code{"target"} and \code{"pst"}: target
    rotation, where the loadings are rotated towards a user-specified target
    matrix (see below).
  \item \code{"biquartimin"} and \code{"bigeomin"}: bifactor rotation, where
    the first factor of the set is a general factor.
  \item \code{"promax"}: the two-step promax method from the stats package
    (varimax followed by an oblique rotation). No standard errors are
    available for this method.
  \item \code{"none"}: no rotation; the unrotated (echelon) solution is
    reported.
}
The rotation is found by minimizing the criterion with the gradient
projection algorithm (\code{algorithm = "gpa"}, the default) or with a
pairwise algorithm (\code{algorithm = "pairwise"}; one pair of factors at a
time, each plane minimized exactly using the gradient of the criterion, or
by a derivative-free search for simplimax), starting from a number of
random starts (\code{rstarts}, default 30) to avoid local minima. Both
algorithms stop when the projected gradient of the criterion falls below
\code{gpa_tol}. The pairwise algorithm is also used for the multiple-group
agreement rotation (\code{mg_agreement = TRUE}).
}

\subsection{Rotation options}{
Options for the rotation are given as a list in the \code{rotation_args}
argument (they can also be given as additional elements of a list in the
\code{rotation} argument, see \code{\link{efa}}). The most important ones
are:
\describe{
  \item{\code{orthogonal}}{Logical. If \code{TRUE}, the factors within a set
    are kept uncorrelated. Default: \code{FALSE} (oblique rotation).}
  \item{\code{row_weights}}{Standardization of the rows of the loading
    matrix before rotation: \code{"default"}, \code{"kaiser"},
    \code{"cureton.mulaik"} or \code{"none"}. The default is
    \code{"kaiser"} for varimax and promax, and \code{"none"} for all other
    criteria.}
  \item{\code{std_ov}}{Logical. If \code{TRUE} (the default), the rotation
    is carried out on the loadings in the standardized (correlation) metric
    of the indicators; the result is transformed back to the original
    metric.}
  \item{\code{rstarts}}{Integer. The number of random starts for the
    rotation algorithm. Default: \code{30}.}
  \item{\code{algorithm}}{\code{"gpa"} (gradient projection, the default) or
    \code{"pairwise"}.}
  \item{\code{pairwise_accelerate}}{Logical. If \code{TRUE} (the default),
    the pairwise algorithm extrapolates the solution across sweeps (an
    Aitken-type step along the last sweep displacement, accepted only if it
    lowers the criterion), which reduces the number of sweeps on slowly
    converging problems.}
  \item{\code{gpa_algorithm}}{Step size strategy of the gradient projection
    algorithm: \code{"bb"} (Barzilai-Borwein with a non-monotone line
    search, the default), \code{"cayley"} (orthogonal rotation only) or
    \code{"legacy"}. See \code{\link{efa}} for \code{gpa_fwindow}.}
  \item{\code{geomin_epsilon}, \code{orthomax_gamma}, \code{oblimin_gamma},
    \code{cf_gamma}, \code{promax_kappa}, \code{simplimax_k}}{Tuning
    constants of the corresponding criteria. Defaults: \code{0.001},
    \code{1}, \code{0}, \code{0}, \code{4} and \code{0}, where
    \code{simplimax_k = 0} means that the \eqn{p(m-1)} smallest loadings
    are driven to zero.}
  \item{\code{target}, \code{target_mask}}{The target matrix (and optional
    mask) for the target rotation methods; see below.}
  \item{\code{reflect}}{Logical. If \code{TRUE} (the default), the sign of
    each rotated factor is chosen such that the sum of its loadings is
    positive.}
  \item{\code{order_lv_by}}{How the rotated factors within a set are ordered:
    \code{"index"} (the default; the factor whose largest loading is on the
    earliest indicator comes first), \code{"sumofsquares"} or
    \code{"none"}. For target rotation, the order is always \code{"none"}
    (the order of the columns of the target matrix). For bifactor rotation,
    the general factor always stays in the first position, and the ordering
    rule is applied to the specific factors only.}
  \item{\code{mg_agreement}}{Logical. For multiple-group models only. If
    \code{TRUE}, the rotated solutions of the groups are made as similar as
    possible (see also \code{mg_agreement_weight},
    \code{mg_agreement_method} and \code{mg_agreement_algorithm}).}
  \item{\code{mg_agreement_group_weights}}{\code{"equal"} (the default)
    or \code{"nobs"}: whether the groups are weighted by their sample
    sizes in the agreement criterion (the per-group criteria by the group
    weight, and each pair of groups in the agreement term by the product
    of the two weights).}
  \item{\code{gpa_tol}, \code{tol}, \code{max_iter}, \code{warn},
    \code{verbose}}{Convergence tolerances, the maximum number of
    iterations, and reporting options of the rotation algorithm.}
}
A complete list (with defaults) is given in the description of the
\code{rotation} argument of \code{\link{efa}}.
}

\subsection{Target rotation and set-specific targets}{
A target rotation rotates the loadings of a set towards a user-specified
\code{target} matrix. The rows of the target matrix correspond to the
indicators of the set (in the order in which they appear in the model syntax),
and the columns to the factors of the set (in the order in which they appear
in the model syntax). The methods differ in how the elements of the target are
used:
\itemize{
  \item \code{"target.strict"}: all elements of the target matrix are
    used.
  \item \code{"target"} (the default alias for \code{"pst"}, partially
    specified target): only the zero elements of the target matrix are
    used; the non-zero elements are ignored. In other words, the target
    specifies where the (near) zero loadings should be, without specifying
    the size of the other loadings.
  \item \code{"pst"} with an explicit \code{target_mask}: the mask is a
    matrix of the same dimensions as the target, with a \code{1} for the
    elements that should be used, and a \code{0} for the elements that
    should be ignored. If elements of the target are \code{NA}, a mask is
    created automatically (the \code{NA} elements are ignored).
}
The \code{target} (and \code{target_mask}) option can take three forms:
\enumerate{
  \item A single matrix, used for every set and every group. This only works
    if all sets have the same number of indicators and factors.
  \item A \emph{nameless} list of matrices, one per group (the same target
    is then used for all sets in a given group).
  \item A \emph{named} list of matrices, where the names are the set labels
    used in the model syntax, so that each set gets its own target matrix
    (the same targets are then used in all groups). For example, if the
    model contains the sets \code{efa("t1")} and \code{efa("t2")}, use
    \code{rotation_args = list(target = list(t1 = target_t1, t2 =
    target_t2))}. A named list must provide a target for \emph{every} set
    (with two or more factors) in the model; omitting a set results in an
    error.
}
}

\subsection{Standard errors}{
Standard errors for the rotated solution are computed analytically. By
default (\code{rotation_se = "bordered"}), the information matrix is
augmented with the rotation constraints; alternatively,
\code{rotation_se = "delta"} applies the delta method to the (numerical)
Jacobian of the rotation. Both methods give the same results up to numerical
precision. When \code{se = "bootstrap"} is requested, the rotation is applied
to every bootstrap sample, and the rotated factors are aligned (reordered and
reflected) to the original solution; see \code{\link{bootstrapLavaan}}. Robust
standard errors (for example, \code{estimator = "MLR"}) are supported as
well.
}
}
\seealso{
\code{\link{efa}} for exploratory factor analysis (a single set of factors
without a structural part), including the complete list of rotation options;
\code{\link{sem}} and \code{\link{lavaan}} for fitting the model;
\code{\link{lavOptions}} for the \code{auto_efa}, \code{rotation},
\code{rotation_args} and \code{rotation_se} options; \code{\link{bootstrapLavaan}}
for bootstrap standard errors of a rotated solution.
}
\references{
Asparouhov, T., & Muthen, B. (2009). Exploratory structural equation
modeling. \emph{Structural Equation Modeling, 16}(3), 397--438.

Marsh, H. W., Morin, A. J. S., Parker, P. D., & Kaur, G. (2014). Exploratory
structural equation modeling: An integration of the best features of
exploratory and confirmatory factor analysis. \emph{Annual Review of Clinical
Psychology, 10}, 85--110.

Browne, M. W. (2001). An overview of analytic rotation in exploratory factor
analysis. \emph{Multivariate Behavioral Research, 36}(1), 111--150.
}
\examples{
## 1. a single set of three exploratory factors (this is equivalent to
##    efa(data = HolzingerSwineford1939, nfactors = 3, ov_names = ...))
model1 <- '
  efa("block1")*f1 +
  efa("block1")*f2 +
  efa("block1")*f3 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
'
fit1 <- sem(model1, data = HolzingerSwineford1939, rotation = "geomin")
summary(fit1, standardized = TRUE)

# the unrotated (echelon) solution
fit1u <- sem(model1, data = HolzingerSwineford1939, rotation = "none")
lavInspect(fit1u, "est")$lambda

# orthogonal varimax rotation, with a different number of random starts
fit1v <- sem(model1, data = HolzingerSwineford1939, rotation = "varimax",
             rotation_args = list(orthogonal = TRUE, rstarts = 10))
lavInspect(fit1v, "est")$lambda

## 2. ESEM: the exploratory factors are regressed on observed covariates
model2 <- '
  efa("block1")*f1 +
  efa("block1")*f2 +
  efa("block1")*f3 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
  f1 + f2 + f3 ~ ageyr + sex
'
fit2 <- sem(model2, data = HolzingerSwineford1939, rotation = "geomin")
parameterEstimates(fit2, output = "text")

## 3. a set of two exploratory factors combined with a confirmatory factor
model3 <- '
  efa("block1")*f1 +
  efa("block1")*f2 =~ x1 + x2 + x3 + x4 + x5 + x6
  speed =~ x7 + x8 + x9
  speed ~ f1 + f2
'
fit3 <- sem(model3, data = HolzingerSwineford1939, rotation = "oblimin")
summary(fit3)

## 4. two sets of factors (two measurement occasions), with a structural
##    part, and a different target matrix per set
pop_model <- '
  f1 =~ 0.8*y1 + 0.7*y2 + 0.6*y3 + 0.3*y4 + 0.0*y5 + 0.0*y6
  f2 =~ 0.0*y1 + 0.0*y2 + 0.3*y3 + 0.6*y4 + 0.7*y5 + 0.8*y6
  f3 =~ 0.8*y7 + 0.7*y8 + 0.6*y9 + 0.3*y10 + 0.0*y11 + 0.0*y12
  f4 =~ 0.0*y7 + 0.0*y8 + 0.3*y9 + 0.6*y10 + 0.7*y11 + 0.8*y12
  f1 ~~ 0.3*f2
  f3 ~~ 0.3*f4
  f3 ~ 0.6*f1
  f4 ~ 0.6*f2
  y1 ~~ 0.2*y7
  y2 ~~ 0.2*y8
  y3 ~~ 0.2*y9
  y4 ~~ 0.2*y10
  y5 ~~ 0.2*y11
  y6 ~~ 0.2*y12
'
set.seed(1234)
Data <- simulateData(pop_model, sample_nobs = 500)

model4 <- '
  # set "t1": two exploratory factors at time 1
  efa("t1")*f1 +
  efa("t1")*f2 =~ y1 + y2 + y3 + y4 + y5 + y6
  # set "t2": two exploratory factors at time 2
  efa("t2")*f3 +
  efa("t2")*f4 =~ y7 + y8 + y9 + y10 + y11 + y12
  # structural part
  f3 ~ f1 + f2
  f4 ~ f1 + f2
  # residual covariances of the repeated indicators
  y1 ~~ y7
  y2 ~~ y8
  y3 ~~ y9
  y4 ~~ y10
  y5 ~~ y11
  y6 ~~ y12
'
# geomin rotation (each set is rotated independently)
fit4 <- sem(model4, data = Data, rotation = "geomin")
summary(fit4)

# target rotation with a different target per set: the rows of each
# target are the indicators of the set, the columns are its factors
target_t1 <- matrix(0, 6, 2)
target_t1[1:3, 1] <- 1
target_t1[4:6, 2] <- 1
target_t2 <- matrix(0, 6, 2)
target_t2[1:4, 1] <- 1
target_t2[5:6, 2] <- 1
fit4t <- sem(model4, data = Data, rotation = "target",
             rotation_args = list(target = list(t1 = target_t1,
                                                t2 = target_t2)))
lavInspect(fit4t, "est")$lambda

# partially specified target: NA elements are ignored (the cross-loadings
# of y3 and y4 are left unspecified in set "t1")
target_t1_pst <- target_t1
target_t1_pst[3:4, ] <- NA
fit4p <- sem(model4, data = Data, rotation = "pst",
             rotation_args = list(target = list(t1 = target_t1_pst,
                                                t2 = target_t2)))
lavInspect(fit4p, "est")$lambda

## 5. multiple groups, with equal loadings across groups
fit5 <- sem(model1, data = HolzingerSwineford1939, group = "school",
            rotation = "geomin", group_equal = "loadings")
lavInspect(fit5, "est")$Pasteur$lambda
}
