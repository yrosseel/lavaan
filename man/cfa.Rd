\name{cfa}
\alias{cfa}
\title{Fit Confirmatory Factor Analysis Models}
\description{
Fit a Confirmatory Factor Analysis (CFA) model.}
\usage{
cfa(model = NULL, data = NULL, ordered = NULL, sampling.weights = NULL,
    sample.cov = NULL, sample.mean = NULL, sample.th = NULL,
    sample.nobs = NULL, group = NULL, cluster = NULL,
    constraints = "", WLS.V = NULL, NACOV = NULL, ov.order = "model",
    ...)

}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See                 
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (eg. the output of the \code{lavaanify()} function) is also
    accepted.}
\item{data}{An optional data frame containing the observed variables used in
    the model. If some variables are declared as ordered factors, lavaan will
    treat them as ordinal variables.}
\item{ordered}{Character vector. Only used if the data is in a data.frame.
    Treat these variables as ordered (ordinal) variables, if they are
    endogenous in the model. Importantly, all other variables will be treated
    as numeric (unless they are declared as ordered in the data.frame.)
    Since 0.6-4, ordered can also be logical. If TRUE, all
    observed endogenous variables are treated as ordered (ordinal).
    If FALSE, all observed endogenous variables are considered to be numeric
    (again, unless they are declared as ordered in the data.frame.)}
\item{sampling.weights}{A variable name in the data frame containing
    sampling weight information. Currently only available for non-clustered 
    data. Depending on the \code{sampling.weights.normalization} option, these
    weights may be rescaled (or not) so that their sum equals the number of 
    observations (total or per group).}
\item{sample.cov}{Numeric matrix. A sample variance-covariance matrix.
    The rownames and/or colnames must contain the observed variable names.
    For a multiple group analysis, a list with a variance-covariance matrix
    for each group.}
\item{sample.mean}{A sample mean vector. For a multiple group analysis,
    a list with a mean vector for each group.}
\item{sample.th}{Vector of sample-based thresholds. For a multiple group 
    analysis, a list with a vector of thresholds for each group.}
\item{sample.nobs}{Number of observations if the full data frame is missing
    and only sample moments are given. For a multiple group analysis, a list
    or a vector with the number of observations for each group.}
\item{group}{Character. A variable name in the data frame defining the groups
in a multiple group analysis.}
\item{cluster}{Character. A (single) variable name in the data frame defining the clusters in a two-level dataset.}
\item{constraints}{Additional (in)equality constraints not yet included in the
    model syntax. See \code{\link{model.syntax}} for more information.}
\item{WLS.V}{A user provided weight matrix to be used by estimator \code{"WLS"};
    if the estimator is \code{"DWLS"}, only the diagonal of this matrix will be 
    used. For a multiple group analysis, a list with a weight matrix
    for each group. The elements of the weight matrix should be in the
    following order (if all data is continuous): first the means (if a 
    meanstructure is involved), then the lower triangular elements of the 
    covariance matrix including the diagonal, ordered column by column. In
    the categorical case: first the thresholds (including the means for 
    continuous variables), then the slopes (if any), the variances of 
    continuous variables (if any), and finally the lower triangular elements 
    of the correlation/covariance matrix excluding the diagonal, ordered 
    column by column.}
\item{NACOV}{A user provided matrix containing the elements of (N times)
    the asymptotic variance-covariance matrix of the sample statistics.
    For a multiple group analysis, a list with an asymptotic 
    variance-covariance matrix for each group. See the \code{WLS.V} 
    argument for information about the order of the elements.}
\item{ov.order}{Character. If \code{"model"} (the default), the order of
    the observed variable names (as reflected for example in the output of
    \code{lavNames()}) is determined by the model syntax. If
    \code{"data"}, the order is determined by the data (either the full
    data.frame or the sample (co)variance matrix). If the \code{WLS.V}
    and/or \code{NACOV} matrices are provided, this argument is currently
    set to \code{"data"}.}
\item{...}{Many more additional options can be defined, using 'name = value'.
    See \code{\link{lavOptions}} for a complete list.}
}
\details{
    The \code{cfa} function is a wrapper for the more general 
    \code{\link{lavaan}} function, using the following default arguments:
    \code{int.ov.free = TRUE}, \code{int.lv.free = FALSE}, 
    \code{auto.fix.first = TRUE} (unless \code{std.lv = TRUE}), 
    \code{auto.fix.single = TRUE}, \code{auto.var = TRUE}, 
    \code{auto.cov.lv.x = TRUE}, \code{auto.efa = TRUE},
    \code{auto.th = TRUE}, \code{auto.delta = TRUE},
    and \code{auto.cov.y = TRUE}.
}
\value{
  An object of class \code{\linkS4class{lavaan}}, for which several methods
  are available, including a \code{summary} method.
}
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36.
\doi{https://doi.org/10.18637/jss.v048.i02}}
\seealso{
\code{\link{lavaan}}
}
\examples{
## The famous Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939)
summary(fit, fit.measures = TRUE)
}
