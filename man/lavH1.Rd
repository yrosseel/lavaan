\name{lavH1}
\alias{lavH1}
\title{The unrestricted (H1) model}
\description{
Compute only the summary statistics of the unrestricted (H1) model: the
saturated mean vector and variance-covariance matrix (and thresholds, if
some variables are ordered) that normally end up in the \code{@h1} slot
of a fitted lavaan object. Unlike \code{\link{sem}} or \code{\link{cfa}},
no model is needed. If a model is provided anyway, it is only used to
determine the relevant variables, their role (exogenous, ordered,
per level, \ldots) and the data-related options, so that the result is
exactly the h1 information a fit of that model would use.

If the data are incomplete and \code{missing = "ml"} (fiml), the summary
statistics are the EM estimates of the saturated model; this is a
convenient way to obtain fiml estimates of standard descriptive
statistics. For multilevel data (using the \code{cluster} argument),
the within- and between-level summary statistics (and the intraclass
correlations) are computed.
}
\usage{
lavH1(model = NULL, data = NULL, ordered = NULL, sampling_weights = NULL,
      group = NULL, cluster = NULL, missing = "default",
      ov_names_x = NULL, ov_names_l = list(),
      estimator = "default", se = "default", test = "none", ...,
      output = "list", wls_v = NULL, gamma = NULL, add_extra = TRUE,
      add_labels = TRUE, add_class = TRUE, drop_list_single_group = TRUE)
}
\arguments{
\item{model}{Optional. Either a description of the user-specified model
    (as in \code{\link{sem}}), or a fitted lavaan object. If a model is
    provided, it is only used to determine the variables, their role,
    and the data-related options. For convenience, the first argument
    may also be the data itself (a \code{data.frame}), in which case the
    \code{data} argument should stay empty.}
\item{data}{A \code{data.frame} containing the observed variables. If no
    model is provided, all variables in the data frame are used (except
    the variables given by the \code{group}, \code{cluster} and
    \code{sampling_weights} arguments).}
\item{ordered}{Character vector declaring additional variables as ordered
    (ordinal). Variables that are declared as ordered factors in the data
    frame are always treated as ordinal. If \code{TRUE}, all variables
    are treated as ordinal.}
\item{sampling_weights}{A variable name in the data frame containing
    sampling weight information.}
\item{group}{A variable name in the data frame defining the groups in a
    multiple group analysis.}
\item{cluster}{A variable name in the data frame defining the clusters in
    a two-level dataset.}
\item{missing}{How to handle missing data. If \code{"ml"} (alias:
    \code{"fiml"}), an EM algorithm is used to compute the saturated
    (h1) summary statistics using all available information. If
    \code{"two.stage"} or \code{"robust.two.stage"} (continuous data
    only), the EM saturated summary statistics are returned together with
    the stage-1 \code{NACOV} (the asymptotic covariance matrix of those
    statistics) when \code{output = "lavMoments"}, so that a follow-up
    analysis based on the summary statistics only reproduces the same
    (two-stage) estimates, standard errors and test statistics as when
    the raw data were provided. See \code{\link{lavOptions}} for further
    options.}
\item{ov_names_x}{Character vector. Variables that need to be treated as
    exogenous. Only used if no model is provided.}
\item{ov_names_l}{List. Only used for two-level data when no model is
    provided: the variable names per level. By default, all variables
    are used at both levels.}
\item{estimator}{The estimator (see \code{\link{lavOptions}}); the
    default corresponds to maximum likelihood (\code{"ML"}) for
    continuous data, and to \code{"WLSMV"} when some variables are
    ordered. The estimator determines, among other things, how the
    summary statistics are computed (for example, pairwise maximum
    likelihood for the polychoric correlations when some variables are
    ordered), and whether the h1 loglikelihood is available.}
\item{se}{Only used if \code{output = "lavaan"}: the standard errors for
    the parameters (= the summary statistics) of the unrestricted
    model. By default, standard errors are computed (matching the
    estimator), so that \code{summary()} shows the standard errors,
    Wald (z) statistics and p-values of the summary statistics; use
    \code{se = "none"} to skip them. See \code{\link{lavOptions}} for
    possible options.}
\item{test}{Only used if \code{output = "lavaan"}. See
    \code{\link{lavOptions}} for possible options.}
\item{...}{Optional parameters that are passed to the
    \code{\link{lavaan}} function (for example \code{fixed.x},
    \code{conditional.x}, \code{h1.missing.method}, \ldots).}
\item{output}{Character. If \code{"list"} (the default), a list with the
    h1 summary statistics per group (and per level in the two-level
    setting): the (residual) variance-covariance matrix (\code{$cov}),
    the mean vector (\code{$mean}, if a meanstructure was used), the
    thresholds (\code{$th}, if some variables are ordered), and the
    regression slopes of the exogenous covariates (\code{$slopes}, if
    \code{conditional.x = TRUE}). If \code{"lavaan"}, the unrestricted
    model itself is returned as a fitted lavaan object. If
    \code{"lavMoments"} (alias \code{"moments"}), an object of class
    \code{lavMoments}: a list of summary statistics (including the
    \code{th.idx} attribute for ordered variables, and the residual
    slopes and exogenous-covariate moments when
    \code{conditional.x = TRUE}) structured so that it can be passed
    \emph{directly} as the \code{data=} argument of
    \code{lavaan}/\code{cfa}/\code{sem} to (re)fit a model from summary
    statistics only. For categorical data, the \code{WLS.V} (weight) and
    \code{NACOV} (\code{gamma}) matrices are included by default, so that
    the refit reproduces the estimates, standard errors and test
    statistics exactly.}
\item{wls_v}{Logical. If \code{TRUE}, the \code{WLS.V} weight matrix (the
    h1 expected information; full or diagonal, depending on the
    estimator) is also computed and returned. If \code{NULL} (the
    default), it is set to \code{TRUE} for categorical data when
    \code{output = "lavMoments"}, and to \code{FALSE} otherwise.}
\item{gamma}{Logical. If \code{TRUE}, the \code{gamma} (\code{NACOV})
    matrix (the asymptotic variance matrix of the summary statistics) is
    also computed and returned. If \code{NULL} (the default), it is set
    to \code{TRUE} for categorical data when
    \code{output = "lavMoments"}, and to \code{FALSE} otherwise.}
\item{add_extra}{Logical. If \code{TRUE}, the list output also contains
    the h1 loglikelihood (\code{$logl}, if available), the number of
    observations (\code{$nobs}), the coverage proportions
    (\code{$coverage}, only if the missing data were handled), and the
    intraclass correlations (\code{$icc}, only in the two-level
    setting). Ignored when \code{output = "lavMoments"}.}
\item{add_labels}{Logical. If \code{TRUE}, variable names are attached
    to the returned statistics.}
\item{add_class}{Logical. If \code{TRUE}, lavaan print classes are
    attached to the returned statistics.}
\item{drop_list_single_group}{Logical. If \code{TRUE}, the per-group
    list structure is dropped when there is only a single group.}
}
\value{
If \code{output = "list"} (the default), a list with the h1 summary
statistics (see the \code{output} argument above). If
\code{output = "lavaan"}, a fitted \code{\linkS4class{lavaan}} object for
the unrestricted model; its \code{summary()} shows the standard errors,
Wald (z) statistics and p-values of the summary statistics (for
example, to flag significant correlations when
\code{correlation = TRUE}). The object may also serve as an ingredient
for two-stage estimation procedures. If \code{output = "lavMoments"}, an
object of class \code{lavMoments} that can be passed directly as the
\code{data=} argument of \code{lavaan}/\code{cfa}/\code{sem}.
}
\section{Reproducing a fit from an output = "lavMoments" object}{
A \code{lavMoments} object reproduces the original raw-data fit
\emph{exactly} for the estimators whose estimates, standard errors and
test statistics are functions of the summary statistics only. This
includes maximum likelihood (\code{"ML"}) and generalized least squares
(\code{"GLS"}); the robust variants (\code{"MLM"}, \code{"MLMV"},
\code{"MLMVS"}, \ldots and the Satorra-Bentler / scaled-shifted /
mean-variance-adjusted tests), provided the \code{gamma} (\code{NACOV})
matrix is requested (\code{gamma = TRUE}); the diagonally weighted least
squares estimators for ordered data (\code{"WLSMV"}, \code{"WLSM"},
\code{"DWLS"}, \ldots, which include \code{WLS.V} and \code{NACOV} by
default); and two-stage maximum likelihood for incomplete continuous data
(\code{missing = "two.stage"} or \code{"robust.two.stage"}).

For the estimators whose \emph{weight matrix} is estimator-specific (for
example the identity weight of the \code{"ULS"} family, or the ADF weight
of continuous \code{"WLS"}), pass the SAME \code{estimator} to
\code{lavH1()} as will be used in the refit, so that the matching
\code{WLS.V} is returned.

A raw-data fit cannot be reproduced from summary statistics when the
standard errors or test statistic depend on the individual cases: the
Huber-White ('robust') sandwich (\code{"MLR"}, \code{se =
"robust.huber.white"}), the first-order/outer-product information
(\code{"MLF"}), cluster-robust (\code{cluster =}) standard errors, the
bootstrap (\code{se = "bootstrap"}), pairwise maximum likelihood
(\code{"PML"}), and direct full-information maximum likelihood for
incomplete data (\code{missing = "ml"}; use \code{missing = "two.stage"}
instead).
}
\seealso{
\code{\link{lavCor}} for (polychoric, polyserial and Pearson)
correlations (which uses the same machinery), and
\code{\link{lavInspect}} with \code{what = "h1"} to extract the h1
summary statistics from a fitted lavaan object.
}
\examples{
# saturated summary statistics of the Holzinger & Swineford data
lavH1(HolzingerSwineford1939[, paste0("x", 1:9)])

# fiml (EM) estimates of the descriptive statistics with missing data
HS.missing <- HolzingerSwineford1939
HS.missing$x5[1:10] <- NA
lavH1(HS.missing[, paste0("x", 1:9)], missing = "ml")

# a model is accepted, but only determines the variables and options
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
h1 <- lavH1(HS.model, data = HolzingerSwineford1939, group = "school")

# the unrestricted model as a fitted lavaan object; summary() shows the
# standard errors, z-statistics and p-values of the summary statistics
fit.h1 <- lavH1(HolzingerSwineford1939[, paste0("x", 1:9)],
                output = "lavaan")
summary(fit.h1)

# which correlations are significant?
fit.cor <- lavH1(HolzingerSwineford1939[, paste0("x", 1:9)],
                 output = "lavaan", correlation = TRUE)
summary(fit.cor)

# categorical data: extract all summary statistics (including WLS.V and
# NACOV), then (re)fit a model from the summary statistics only
Data <- HolzingerSwineford1939[, paste0("x", 1:9)]
Data[] <- lapply(Data, cut, breaks = 3, labels = FALSE)
Data[] <- lapply(Data, ordered)
model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '
moments <- lavH1(model, data = Data, output = "lavMoments")
moments
fit <- cfa(model, data = moments)   # identical to cfa(model, data = Data)

# two-stage ML for incomplete (continuous) data: the EM saturated moments
# plus the stage-1 NACOV allow a moments-only refit to reproduce the exact
# two-stage estimates, standard errors and (scaled) test statistics
HS.incomplete <- HolzingerSwineford1939[, paste0("x", 1:9)]
HS.incomplete$x1[1:20] <- HS.incomplete$x5[15:40] <- NA
model2 <- ' visual =~ x1 + x2 + x3
            textual =~ x4 + x5 + x6 '
moments2 <- lavH1(model2, data = HS.incomplete, missing = "two.stage",
                  output = "lavMoments")
fit2 <- cfa(model2, data = moments2)
# same as: cfa(model2, data = HS.incomplete, missing = "two.stage")
}
