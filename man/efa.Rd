\name{efa}
\alias{efa}
\alias{rotation}
\title{Exploratory Factor Analysis}
\description{
Fit one or more Exploratory Factor Analysis (EFA) model(s).}
\usage{
efa(data = NULL, nfactors = 1L, sample.cov = NULL, sample.nobs = NULL,
    rotation = "geomin", rotation.args = list(), ov.names = NULL,
    bounds = "pos.var", ..., output = "efa")
}
\arguments{
\item{data}{A data frame containing the observed variables we need for the
    EFA. If only a subset of the observed variables is needed, use the
    \code{ov.names} argument.}
\item{nfactors}{Integer or Integer vector. The desired number of factors
    to extract. Can be a single number, or a vector of numbers
    (e.g., \code{nfactors = 1:4}.), For each different number, a model
    is fitted.}
\item{sample.cov}{Numeric matrix. A sample variance-covariance matrix.
    The rownames and/or colnames must contain the observed variable names.
    Unlike sem and CFA, the matrix may be a correlation matrix.}
\item{sample.nobs}{Number of observations if the full data frame is missing
    and only the sample variance-covariance matrix is given.}
\item{rotation}{Character or a list with first element a character variable.
    The character variable is the rotation method to be used. Possible options
    are varimax, quartimax, orthomax, oblimin, quartimin, geomin, promax,
    entropy, mccammon, infomax, tandem1, tandem2, oblimax, bentler, simplimax, 
    target.strict, target (alias for pst), pst (=partially specified target), 
    cf, crawford-ferguson,
    cf-quartimax,  cf-varimax, cf-equamax,
    cf-parsimax, cf-facparsim, biquartimin, bigeomin. The latter two are
    for bifactor rotation only. The rotation algorithms (except promax
    and target) are similar to those from the GPArotation package, but have 
    been reimplemented for better control. The promax method is taken from the
    stats package. The target.strict method is equal to the target method in
    the GPArotation package. The target method is in fact the pst method where
    all non-zero elements (in the target matrix) are ignored.
    Options for the rotation algorithm can be provided in the list after the method.
    The
   default options (and their alternatives) are \code{orthogonal = FALSE},
   \code{row.weights = "default"} (or \code{"kaiser"},
   \code{"cureton.mulaik"} or \code{"none"}), \code{std.ov = TRUE},
   \code{algorithm = "gpa"} (or \code{"pairwise"}), \code{rstarts = 30},
   \code{gpa.tol = 1e-05}, \code{tol = 1e-08}, \code{max.iter = 10000L},
   \code{warn = FALSE}, \code{verbose = FALSE}, \code{reflect = TRUE},
   \code{order.lv.by = "index"} (or \code{"sumofsquares"} or \code{"none"}).
   Other options are specific for a particular rotation criterion: 
   \code{geomin.epsilon = 0.001}, \code{orthomax.gamma = 1},
   \code{promax.kappa = 4},
   \code{cf.gamma = 0}, and \code{oblimin.gamma = 0}.}
\item{rotation.args}{List. Options related to the rotation algorithm. DEPRECATED.
   The options should now be given in the \code{rotation} argument as a list.}
\item{ov.names}{Character vector. The variables names that are needed for
   the EFA. Should be a subset of the variables names in the data.frame.
   By default (if NULL), all the variables in the data are used.}
\item{bounds}{Per default, \code{bounds = "pos.var"} forces all variances
   of both observed and latent variables to be strictly nonnegative. See
   the entry in \code{\link{lavOptions}} for more options.}
\item{...}{Aditional options to be passed to lavaan, using 'name = value'.
    See \code{\link{lavOptions}} for a complete list.}
\item{output}{Character. If \code{"efa"} (the default), the output mimics
    the typical output of an EFA. If \code{"lavaan"}, a lavaan object returned.
    The latter is only possible if nfactors contains a single (integer) number.}
}
\details{
    The \code{efa} function is essentially a wrapper around the
    \code{lavaan} function. It generates the model syntax (for a given number
    of factors) and then calls \code{lavaan()} treating the factors as 
    a single block that should be rotated. The function only supports
    a single group. Categorical data is handled as usual by first computing
    an appropriate (e.g., tetrachoric or polychoric) correlation matrix,
    which is then used as input for the EFA. 
    There is also (limited) support for
    twolevel data. The same number of factors is then extracted at the
    within and the between level.
    The promax rotation method (taken from the stats package) is only 
    provided for convenience. Because promax is a two-step algorithm (first 
    varimax, then oblique rotation to get simple structure), it does not 
    use the gpa or pairwise rotation algorithms, and as a result, no
    standard errors are provided.
}
\value{
  If \code{output = "lavaan"}, an object of class 
  \code{\linkS4class{lavaan}}. If \code{output = "efa"},
  a list of class \code{efaList} for which a \code{print()}, 
  \code{summary()} and \code{fitMeasures()} method are available. Because
  we added the (standardized) loadings as an extra element, the \code{loadings}
  function (which is not a generic function) from the stats package will
  also work on \code{efaList} objects.
}
\seealso{
    \code{\link{lav_efalist_summary}} for a summary method if the output is
    of class \code{efaList}. 
}
\examples{
## The famous Holzinger and Swineford (1939) example
fit <- efa(data = HolzingerSwineford1939, 
           ov.names = paste("x", 1:9, sep = ""),
           nfactors = 1:3,
           rotation = list("geomin", geomin.epsilon = 0.01, rstarts = 1))
summary(fit, nd = 3L, cutoff = 0.2, dot.cutoff = 0.05)
fitMeasures(fit, fit.measures = "all")

# target rotation
target <- matrix(0, 9, 3)
target[1:3, 1] <- 1
target[4:6, 2] <- 1
target[7:9, 3] <- 1
fit2 <- efa(data = HolzingerSwineford1939,
            ov.names = paste("x", 1:9, sep = ""), 
            nfactors = 3,
            rotation = list("target", target = target))
summary(fit2)
}
