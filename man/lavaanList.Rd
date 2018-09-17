\name{lavaanList}
\alias{lavaanList}
\alias{semList}
\alias{cfaList}
\title{Fit List of Latent Variable Models}
\description{
Fit the same latent variable model, for a (potentially large) number of datasets.}
\usage{
lavaanList(model = NULL, dataList = NULL, dataFunction = NULL, 
    dataFunction.args = list(), ndat = length(dataList), cmd = "lavaan", 
    ..., store.slots = c("partable"), FUN = NULL, show.progress = FALSE,
    store.failed = FALSE, parallel = c("no", "multicore", "snow"), 
    ncpus = max(1L, parallel::detectCores() - 1L), cl = NULL, iseed = NULL)

semList(model = NULL, dataList = NULL, dataFunction = NULL,
    dataFunction.args = list(), ndat = length(dataList),
    ..., store.slots = c("partable"), FUN = NULL, show.progress = FALSE,
    store.failed = FALSE, parallel = c("no", "multicore", "snow"), 
    ncpus = max(1L, parallel::detectCores() - 1L), cl = NULL, iseed = NULL)

cfaList(model = NULL, dataList = NULL, dataFunction = NULL,
    dataFunction.args = list(), ndat = length(dataList),
    ..., store.slots = c("partable"), FUN = NULL, show.progress = FALSE,
    store.failed = FALSE, parallel = c("no", "multicore", "snow"), 
    ncpus = max(1L, parallel::detectCores() - 1L), cl = NULL, iseed = NULL)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See            
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (eg. the output of the \code{lavaanify()} function) is also
    accepted.}
\item{dataList}{List. Each element contains a full data frame containing
      the observed variables used in the model.}
\item{dataFunction}{Function. A function that generated a full data frame
      containing the observed variables used in the model. It can also be
      a matrix, if the columns are named.}
\item{dataFunction.args}{List. Optional list of arguments that are passed
    to the \code{dataFunction} function.}
\item{ndat}{Integer. The number of datasets that should be generated using
    the \code{dataFunction} function.}
\item{cmd}{Character. Which command is used to run the sem models. The possible
    choices are \code{"sem"}, \code{"cfa"} or \code{"lavaan"}, determining how
    we deal with default options.}
\item{\dots}{Other named arguments for \code{lavaan} function.}
\item{store.slots}{Character vector. Which slots (from a lavaan object) 
    should be stored for each dataset? The possible choices are
    \code{"timing"}, \code{"partable"}, \code{"data"}, \code{"samplestats"},
    \code{"vcov"}, \code{"test"}, \code{"optim"}, \code{"h1"}, \code{"loglik"},
    or \code{"implied"}. Finally, \code{"all"} selects all slots.}
\item{FUN}{Function. A function which when applied to the
    \code{\linkS4class{lavaan}} object returns the information of interest.}
\item{store.failed}{Logical. If \code{TRUE}, write (to tempdir()) the dataset
     and (if available) the fitted object when the estimation for a particular
     dataset somehow failed. This will allow posthoc inspection of the problem.}
\item{parallel}{The type of parallel operation to be used (if any). If
missing, the default is \code{"no"}.}
\item{ncpus}{Integer. The number of processes to be used in parallel operation:
typically one would chose this to the number of available CPUs.}
\item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
\code{parallel = "snow"}.  If not supplied, a cluster on the local machine is
created for the duration of the \code{lavaanList} call.}
\item{iseed}{An integer to set the seed. Or NULL if no reproducible seeds are
  needed. To make this work, make sure the first
  RNGkind() element is \code{"L'Ecuyer-CMRG"}. You can check this by
  typing \code{RNGkind()} in the console. You can set it by typing
  \code{RNGkind("L'Ecuyer-CMRG")}, before the lavaanList functions are called.}
\item{show.progress}{If \code{TRUE}, show information for each dataset.}
}
\value{
  An object of class \code{\linkS4class{lavaanList}}, for which several methods
  are available, including a \code{summary} method.
}
\seealso{
class \code{\linkS4class{lavaanList}}
}
\examples{
# The Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

# a data generating function
generateData <- function() simulateData(HS.model, sample.nobs = 100)

set.seed(1234)
fit <- semList(HS.model, dataFunction = generateData, ndat = 5,
               store.slots = "partable")

# show parameter estimates, per dataset
coef(fit)
}
