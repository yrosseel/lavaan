\name{efa}
\alias{efa}
\alias{rotation}
\title{Exploratory Factor Analysis}
\description{
Fit one or more Exploratory Factor Analysis (EFA) model(s).}
\usage{
efa(data = NULL, nfactors = 1L, sample.cov = NULL, sample.nobs = NULL,
    rotation = "geomin", rotation.args = list(), ov.names = names(data),
    bounds = "pos.var", baseline = FALSE, ..., output = "efa")
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
\item{rotation}{Character. The rotation method to be used. Possible options
    are varimax, quartimax, orthomax, oblimin, quartimin, geomin, 
    entropy, mccammon, infomax, tandem1, tandem2, oblimax, bentler, simplimax, 
    target, pst (=partially specified target), cf, crawford-ferguson,
    cf-quartimax,  cf-varimax, cf-equamax,
    cf-parsimax, cf-facparsim, biquartimin, bigeomin. The latter two are
    for bifactor rotation only. The rotation algorithms are similar to
    those from the GPArotation package, but have been reimplemented for
    efficiency.}
\item{rotation.args}{List. Options related to the rotation algorithm. The
   default options (and their alternatives) are \code{orthogonal = FALSE},
   \code{row.weights = "default"} (or \code{"kaiser"},
   \code{"cureton.mulaik"} or \code{"none"}), \code{std.ov = TRUE},
   \code{algorithm = "gpa"} (or \code{"pairwise"}), \code{rstarts = 30},
   \code{gpa.tol = 1e-05}, \code{tol = 1e-08}, \code{max.iter = 10000L},
   \code{warn = FALSE}, \code{verbose = FALSE}, \code{reflect = TRUE},
   \code{order.lv.by = "index"} (or \code{"sumofsquares"} or \code{"none"}).
   Other options are specific for a particular rotation criterion: 
   \code{geomin.epsilon = 0.001}, \code{orthomax.gamma = 1},
   \code{cf.gamma = 0}, and \code{oblimin.gamma = 0}.}
\item{ov.names}{Character vector. The variables names that are needed for
   the EFA. Should be a subset of the variables names in the data.frame.
   By default, all the variables in the data are used.}
\item{bounds}{Per default, \code{bounds = "pos.var"} forces all variances
   of both observed and latent variables to be strictly nonnegative. See
   the entry in \code{\link{lavOptions}} for more options.}
\item{baseline}{Logical. The default is \code{FALSE}. This implies that
   fit measures that rely on the baseline model (such as CFI and TLI) cannot
   be computed.}
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
    and use this correlation matrix as input for the EFA. 
    There is also support for
    twolevel data. The same number of factors is then extracted at the
    within and the between level.
}
\value{
  If \code{output = "lavaan"}, an object of class 
  \code{\linkS4class{lavaan}}. If \code{output = "efa"},
  a list of class \code{efaList} for which a \code{print()} and a 
  \code{summary()} method are available.
}
\examples{
## The famous Holzinger and Swineford (1939) example
fit <- efa(data = HolzingerSwineford1939, 
           ov.names = paste("x", 1:9, sep = ""),
           nfactors = 1:3,
           rotation = "geomin",
           rotation.args = list(geomin.epsilon = 0.01, rstarts = 1))
summary(fit, nd = 3L, cutoff = 0.2, dot.cutoff = 0.05)
}
