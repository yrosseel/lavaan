\name{lavaanList}
\alias{lavaanList}
\alias{semList}
\alias{cfaList}
\title{Fit List of Latent Variable Models}
\description{
Fit the same latent variable model to a (potentially large) number of datasets.}
\usage{
lavaanList(model = NULL, data_list = NULL, data_function = NULL, 
    data_function_args = list(), ndat = length(data_list), cmd = "lavaan", 
    ..., store_slots = c("partable"), fun = NULL, show_progress = FALSE,
    store_failed = FALSE, parallel = c("no", "multicore", "snow"), 
    ncpus = min(2L, max(1L, parallel::detectCores() - 1L, na.rm = TRUE)),
    cl = NULL, iseed = NULL)

semList(model = NULL, data_list = NULL, data_function = NULL,
    data_function_args = list(), ndat = length(data_list),
    ..., store_slots = c("partable"), fun = NULL, show_progress = FALSE,
    store_failed = FALSE, parallel = c("no", "multicore", "snow"), 
    ncpus = min(2L, max(1L, parallel::detectCores() - 1L, na.rm = TRUE)),
    cl = NULL, iseed = NULL)

cfaList(model = NULL, data_list = NULL, data_function = NULL,
    data_function_args = list(), ndat = length(data_list),
    ..., store_slots = c("partable"), fun = NULL, show_progress = FALSE,
    store_failed = FALSE, parallel = c("no", "multicore", "snow"), 
    ncpus = min(2L, max(1L, parallel::detectCores() - 1L, na.rm = TRUE)),
    cl = NULL, iseed = NULL)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See 
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (e.g., the output of the \code{lavParTable()} function) is also
    accepted.}
\item{data_list}{List. Each element contains a full data frame containing
      the observed variables used in the model.}
\item{data_function}{Function. A function that generates a full data frame
      containing the observed variables used in the model. It can also be
      a matrix, if the columns are named.}
\item{data_function_args}{List. Optional list of arguments that are passed
    to the \code{data_function} function.}
\item{ndat}{Integer. The number of datasets that should be generated using
    the \code{data_function} function.}
\item{cmd}{Character. The command used to run the SEM models. The possible
    choices are \code{"sem"}, \code{"cfa"} or \code{"lavaan"}, which determine
    how the default options are handled.}
\item{\dots}{Other named arguments for the \code{lavaan} function.}
\item{store_slots}{Character vector. Which slots (from a lavaan object) 
    should be stored for each dataset? The possible choices are
    \code{"timing"}, \code{"partable"}, \code{"data"}, \code{"samplestats"},
    \code{"vcov"}, \code{"test"}, \code{"optim"}, \code{"h1"}, \code{"loglik"},
    or \code{"implied"}. Finally, \code{"all"} selects all slots.}
\item{fun}{Function. A function which when applied to the
    \code{\linkS4class{lavaan}} object returns the information of interest.}
\item{store_failed}{Logical. If \code{TRUE}, write (to tempdir()) the dataset
     and (if available) the fitted object when the estimation for a particular
     dataset failed for some reason. This allows post hoc inspection of the problem.}
\item{parallel}{The type of parallel operation to be used (if any). If
missing, the default is \code{"no"}.}
\item{ncpus}{Integer. The number of processes to be used in parallel operation:
By default, this is 2 (following CRAN policy). You may set this to a higher
value (for example, the number of available CPUs) if you wish to use more
processes.}
\item{cl}{An optional \pkg{parallel} or \pkg{snow} cluster for use if
\code{parallel = "snow"}.  If not supplied, a cluster on the local machine is
created for the duration of the \code{lavaanList} call.}
\item{iseed}{An integer to set the seed. Or NULL if no reproducible seeds are
  needed. To make this work, make sure the first
  RNGkind() element is \code{"L'Ecuyer-CMRG"}. You can check this by
  typing \code{RNGkind()} in the console. You can set it by typing
  \code{RNGkind("L'Ecuyer-CMRG")}, before the lavaanList functions are called.}
\item{show_progress}{If \code{TRUE}, show information for each dataset.}
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
generateData <- function() lavSimulateData(HS.model, sample_nobs = 100)

set.seed(1234)
fit <- semList(HS.model, data_function = generateData, ndat = 5,
               store_slots = "partable")

# show parameter estimates, per dataset
coef(fit)
}
