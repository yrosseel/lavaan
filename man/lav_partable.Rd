\name{lav_partable}
\alias{lav_partable_independence}
\alias{lav_partable_unrestricted}
\alias{lav_partable_df}
\alias{lav_partable_ndat}
\alias{lav_partable_npar}
\alias{lav_partable_labels}
\alias{lav_partable_from_lm}
\alias{lav_partable_complete}
\alias{lav_partable_attributes}
\alias{lav_partable_merge}
\alias{lav_partable_add}
\title{lavaan partable functions}
\description{Utility functions related to the parameter table (partable)}
\usage{
# extract information from a parameter table
lav_partable_df(partable)
lav_partable_ndat(partable)
lav_partable_npar(partable)
lav_partable_attributes(partable, pta = NULL)

# generate parameter labels
lav_partable_labels(partable, blocks = c("group", "level"), 
                    group.equal = "", group.partial = "", type = "user")

# generate parameter table for specific models
lav_partable_independence(lavobject = NULL, lavdata = NULL,
    lavpta = NULL, lavoptions = NULL, lavsamplestats = NULL,
    lavh1 = NULL,
    sample.cov = NULL, sample.mean = NULL, sample.slopes = NULL,
    sample.th = NULL, sample.th.idx = NULL,
    sample.cov.x = NULL, sample.mean.x = NULL)

lav_partable_unrestricted(lavobject = NULL, lavdata = NULL, 
    lavpta = NULL, lavoptions = NULL, lavsamplestats = NULL, 
    lavh1 = NULL,
    sample.cov = NULL, sample.mean = NULL, sample.slopes = NULL, 
    sample.th = NULL, sample.th.idx = NULL,
    sample.cov.x = NULL, sample.mean.x = NULL) 

lav_partable_from_lm(object, est = FALSE, label = FALSE, 
    as.data.frame. = FALSE)

# complete a parameter table only containing a few columns (lhs,op,rhs)
lav_partable_complete(partable = NULL, start = TRUE)

# merge two parameter tables
lav_partable_merge(pt1 = NULL, pt2 = NULL, remove.duplicated = FALSE, 
    fromLast = FALSE, warn = TRUE)

# add a single parameter to an existing parameter table
lav_partable_add(partable = NULL, add = list())
}
\arguments{
\item{partable}{A parameter table. See \code{\link{lavParTable}} 
for more information.}
\item{blocks}{Character vector. Which columns in the parameter table should be
  taken to distinguish between different blocks of parameters (and hence
  be given different labels)? If \code{"blocks"} includes \code{"group"}, a
  suffix \code{".g"} and the group number (or group label) is added for 
  the parameters of  all but the first group. If \code{"blocks"} includes 
  \code{"level"}, a suffix \code{".l"} and the level number is added for 
  the parameters of all but the first level. If \code{"blocks"} includes, 
  say \code{"foo"}, a suffix \code{".foo"} and the corresponding value of
  \code{"foo"} is added to all parameters.}
\item{group.equal}{The same options can be used here as in the fitting functions. Parameters that are constrained to be equal across groups will be given the
same label.}
\item{group.partial}{A vector of character strings containing the labels
    of the parameters which should be free in all groups.}
\item{type}{Character string. Can be either `user' or `free' to select all
entries or only the free parameters from the parameter table respectively.}
\item{lavobject}{An object of class `lavaan'. If this argument is provided,
it should be the only argument. All the values for the other arguments 
are extracted from this object.}
\item{lavdata}{An object of class `lavData'. The Data slot from a lavaan object.}
\item{lavoptions}{A named list. The Options lsot from a lavaan object.}
\item{lavsamplestats}{An object of class `lavSampleStats'. The SampleStats
slot from a lavaan object.}
\item{lavh1}{A named list. The h1 slot from a lavaan object.}
\item{lavpta}{The pta (parameter table attributes) slot from a lavaan object.}
\item{sample.cov}{Optional list of numeric matrices. 
Each list element contains a sample variance-covariance matrix for this group.
If provided, these values will be used as starting values.}
\item{sample.mean}{Optional list of numeric vectors. 
Each list element contains a sample mean vector for this group. 
If provided, these values will be used as starting values.}
\item{sample.slopes}{Optional list of numeric matrices.
Each list element contains the sample slopes for this group (only used
when \code{conditional.x = TRUE}).
If provided, these values will be used as starting values.}
\item{sample.th}{Optional list of numeric vectors. 
Each list element contains a vector of sample thresholds for this group.
If provided (and also sample.th.idx is provided), 
these values will be used as starting values.}
\item{sample.th.idx}{Optional list of integers. Each list contains the
threshold indices for this group.}
\item{sample.cov.x}{Optional list of numeric matrices.  Each list element
contains a sample variance-covariance matrix for the exogenous variables
for this group (only used when \code{conditional.x = TRUE}).  If provided,
these values will be used as starting values.}
\item{sample.mean.x}{Optional list of numeric vectors.
Each list element contains a sample mean vector for the exogenous variables
for this group (only used when \code{conditional.x = TRUE}).
If provided, these values will be used as starting values.}
\item{est}{Logical. If TRUE, include the fitted estimates in the parameter table.}
\item{label}{Logical. If TRUE, include parameter labels in the parameter table.}
\item{as.data.frame.}{Logical. If TRUE, return the parameter table as a data.frame.}
\item{object}{An object of class \code{lm}.}
\item{start}{Logical. If TRUE, include a start column, based on 
the simple method for generating starting values.}
\item{pta}{A list containing parameter attributes.}
\item{pt1}{A parameter table.}
\item{pt2}{A parameter table.}
\item{remove.duplicated}{Logical. If \code{TRUE}, remove duplicated elements
when merging two parameter tables.}
\item{fromLast}{Logical. If \code{TRUE}, duplicated elements are considered
from the bottom of the merged parameter table.}
\item{warn}{Logical. If \code{TRUE}, a warning is produced when duplicated
elements are found, when merging two parameter tables.}
\item{add}{A named list.  A single row of a parameter table as a named list.}
}
\examples{
# generate syntax for an independence model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
lav <- lav_partable_independence(fit)
as.data.frame(lav, stringsAsFactors = FALSE)


# how many free parameters?
lav_partable_npar(lav)

# how many sample statistics?
lav_partable_ndat(lav)
}

