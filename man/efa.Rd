\name{efa}
\alias{efa}
\alias{rotation}
\title{Exploratory Factor Analysis}
\description{
Fit one or more Exploratory Factor Analysis (EFA) model(s).}
\usage{
efa(data = NULL, nfactors = 1L, sample_cov = NULL, sample_nobs = NULL,
    rotation = "geomin", rotation_args = list(), ov_names = NULL,
    bounds = "pos.var", ..., output = "efa")
}
\arguments{
\item{data}{A data frame containing the observed variables we need for the
    EFA. If only a subset of the observed variables is needed, use the
    \code{ov.names} argument.}
\item{nfactors}{Integer, integer vector, or a list. The desired number of
    factors to extract. Can be a single number, or a vector of numbers
    (e.g., \code{nfactors = 1:4}); for each different number, a model is
    fitted. When the model has multiple blocks (because the data is
    multilevel and/or multigroup), the same number of factors is extracted in
    every block. To extract a different number of factors per block, supply a
    \emph{list}: each element of the list is one fitted model and contains a
    vector with the number of factors per block (or a single number that is
    recycled to all blocks). Blocks are ordered group-major, level-minor, ie
    block \eqn{= (g - 1) \times nlevels + level}. For example, for a twolevel
    model, \code{nfactors = list(c(1, 2))} extracts one factor at the within
    level and two factors at the between level.}
\item{sample_cov}{Numeric matrix. A sample variance-covariance matrix.
    The rownames and/or colnames must contain the observed variable names.
    Unlike sem and CFA, the matrix may be a correlation matrix.}
\item{sample_nobs}{Number of observations if the full data frame is missing
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
    The default options (and their alternatives) are \code{orthogonal = FALSE},
   \code{row_weights = "default"} (or \code{"kaiser"},
   \code{"cureton.mulaik"} or \code{"none"}), \code{std_ov = TRUE},
   \code{algorithm = "gpa"} (or \code{"pairwise"}), \code{rstarts = 30},
   \code{gpa_tol = 1e-05}, \code{tol = 1e-08}, \code{max_iter = 10000L},
   \code{warn = FALSE}, \code{verbose = FALSE}, \code{reflect = TRUE},
   \code{order_lv_by = "index"} (or \code{"sumofsquares"} or \code{"none"}).
   Other options are specific for a particular rotation criterion:
   \code{geomin_epsilon = 0.001}, \code{orthomax_gamma = 1},
   \code{promax_kappa = 4},
   \code{cf_gamma = 0}, and \code{oblimin_gamma = 0}.

   For the target rotation methods (\code{target.strict}, \code{target} and
   \code{pst}), the \code{target} option (and, for \code{pst}, the optional
   \code{target_mask} option) can take three forms: 1) a single matrix,
   which is used for all groups and all efa blocks; 2) a \emph{nameless}
   list of matrices, one per group; or 3) a \emph{named} list of matrices,
   where the names are the efa block labels used in the model syntax (e.g.,
   \code{list(a = target.a, b = target.b)} if the model contains the efa
   blocks \code{efa("a")} and \code{efa("b")}), so that each efa block gets
   its own target matrix (the same targets are then used in all groups).
   The latter is only relevant when the model syntax (of \code{sem()} or
   \code{cfa()}) defines multiple efa blocks (set-ESEM). Note that a named
   list must provide a target matrix for \emph{every} efa block (with two
   or more factors) in the model; omitting a block results in an error, as
   the rotation criterion would be undefined for that block.}
\item{rotation_args}{List. Options related to the rotation algorithm. DEPRECATED.
   The options should now be given in the \code{rotation} argument as a list.}
\item{ov_names}{Character vector. The variables names that are needed for
   the EFA. Should be a subset of the variables names in the data.frame.
   By default (if NULL), all the variables in the data are used.}
\item{bounds}{Per default, \code{bounds = "pos.var"} forces all variances
   of both observed and latent variables to be strictly nonnegative. See
   the entry in \code{\link{lavOptions}} for more options.}
\item{...}{Additional options to be passed to lavaan, using 'name = value'.
    See \code{\link{lavOptions}} for a complete list.}
\item{output}{Character. If \code{"efa"} (the default), the output mimics
    the typical output of an EFA. If \code{"lavaan"}, a lavaan object is returned.
    The latter is only possible if nfactors contains a single (integer) number.}
}
\details{
    The \code{efa} function is essentially a wrapper around the
    \code{lavaan} function. It generates the model syntax (for a given number
    of factors) and then calls \code{lavaan()} treating the factors (in each
    block) as a set that should be rotated. Each block is rotated
    independently. Categorical data is handled as usual by first computing
    an appropriate (e.g., tetrachoric or polychoric) correlation matrix,
    which is then used as input for the EFA.
    Multiple groups (via the \code{group=} argument) and twolevel data
    (via the \code{cluster=} argument) are supported, and these may be
    combined. By default, the same number of factors is extracted in every
    block, but a different number of factors per block can be requested by
    supplying \code{nfactors} as a list (see the \code{nfactors} argument).
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
           rotation = list("geomin", geomin_epsilon = 0.01, rstarts = 1))
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

\dontrun{
# twolevel EFA with a different number of factors per level:
# one factor at the within level, two factors at the between level
fit3 <- efa(data = Demo.twolevel, ov.names = paste0("y", 1:6),
            cluster = "cluster", nfactors = list(c(1, 2)))
summary(fit3)
}
}
