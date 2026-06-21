\name{lav_partable}
\alias{lav_pt_independence}
\alias{lav_pt_unrestricted}
\alias{lav_pt_df}
\alias{lav_pt_ndat}
\alias{lav_pt_npar}
\alias{lav_pt_ngroups}
\alias{lav_pt_labels}
\alias{lav_pt_from_lm}
\alias{lav_pt_complete}
\alias{lav_pt_attributes}
\alias{lav_pt_merge}
\alias{lav_pt_add}
\title{lavaan partable functions}
\description{Utility functions related to the parameter table (partable)}
\usage{
# extract information from a parameter table
lav_pt_df(partable)
lav_pt_ndat(partable)
lav_pt_npar(partable)
lav_pt_ngroups(partable)
lav_pt_attributes(partable, pta = NULL)

# generate parameter labels
lav_pt_labels(partable, blocks = c("group", "level"), 
                    group_equal = "", group_partial = "", type = "user")

# generate parameter table for specific models
lav_pt_independence(lavobject = NULL, lavdata = NULL,
    lavpta = NULL, lavoptions = NULL, lavsamplestats = NULL,
    lavh1 = NULL,
    sample_cov = NULL, sample_mean = NULL, sample_slopes = NULL,
    sample_th = NULL, sample_th_idx = NULL,
    sample_cov_x = NULL, sample_mean_x = NULL)

lav_pt_unrestricted(lavobject = NULL, lavdata = NULL, 
    lavpta = NULL, lavoptions = NULL, lavsamplestats = NULL, 
    lavh1 = NULL,
    sample_cov = NULL, sample_mean = NULL, sample_slopes = NULL, 
    sample_th = NULL, sample_th_idx = NULL,
    sample_cov_x = NULL, sample_mean_x = NULL) 

lav_pt_from_lm(object, est = FALSE, label = FALSE, 
    as_data_frame = FALSE)

# complete a parameter table only containing a few columns (lhs,op,rhs)
lav_pt_complete(partable = NULL, start = TRUE)

# merge two parameter tables
lav_pt_merge(pt1 = NULL, pt2 = NULL, remove_duplicated = FALSE, 
    from_last = FALSE, warn = TRUE)

# add a single parameter to an existing parameter table
lav_pt_add(partable = NULL, add = list())
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
\item{group_equal}{The same options can be used here as in the fitting functions.
Parameters that are constrained to be equal across groups will be given the
same label.}
\item{group_partial}{A vector of character strings containing the labels
    of the parameters which should be free in all groups.}
\item{type}{Character string. Either `user' to select all entries, or
`free' to select only the free parameters from the parameter table.}
\item{lavobject}{An object of class `lavaan'. If this argument is provided,
it should be the only argument. All the values for the other arguments 
are extracted from this object.}
\item{lavdata}{An object of class `lavData'. The Data slot from a lavaan object.}
\item{lavoptions}{A named list. The Options slot from a lavaan object.}
\item{lavsamplestats}{An object of class `lavSampleStats'. The SampleStats
slot from a lavaan object.}
\item{lavh1}{A named list. The h1 slot from a lavaan object.}
\item{lavpta}{The pta (parameter table attributes) slot from a lavaan object.}
\item{sample_cov}{Optional list of numeric matrices. 
Each list element contains a sample variance-covariance matrix for this group.
If provided, these values will be used as starting values.}
\item{sample_mean}{Optional list of numeric vectors. 
Each list element contains a sample mean vector for this group. 
If provided, these values will be used as starting values.}
\item{sample_slopes}{Optional list of numeric matrices.
Each list element contains the sample slopes for this group (only used
when \code{conditional.x = TRUE}).
If provided, these values will be used as starting values.}
\item{sample_th}{Optional list of numeric vectors. 
Each list element contains a vector of sample thresholds for this group.
If provided (and also sample.th.idx is provided), 
these values will be used as starting values.}
\item{sample_th_idx}{Optional list of integers. Each list contains the
threshold indices for this group.}
\item{sample_cov_x}{Optional list of numeric matrices.  Each list element
contains a sample variance-covariance matrix for the exogenous variables
for this group (only used when \code{conditional.x = TRUE}).  If provided,
these values will be used as starting values.}
\item{sample_mean_x}{Optional list of numeric vectors.
Each list element contains a sample mean vector for the exogenous variables
for this group (only used when \code{conditional.x = TRUE}).
If provided, these values will be used as starting values.}
\item{est}{Logical. If TRUE, include the fitted estimates in the parameter table.}
\item{label}{Logical. If TRUE, include parameter labels in the parameter table.}
\item{as_data_frame}{Logical. If TRUE, return the parameter table as a data.frame.}
\item{object}{An object of class \code{lm}.}
\item{start}{Logical. If TRUE, include a start column, based on 
the simple method for generating starting values.}
\item{pta}{A list containing parameter attributes.}
\item{pt1}{A parameter table.}
\item{pt2}{A parameter table.}
\item{remove_duplicated}{Logical. If \code{TRUE}, remove duplicated elements
when merging two parameter tables.}
\item{from_last}{Logical. If \code{TRUE}, duplicated elements are considered
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
lav <- lav_pt_independence(fit)
as.data.frame(lav, stringsAsFactors = FALSE)

# how many free parameters?
lav_pt_npar(lav)

# how many sample statistics?
lav_pt_ndat(lav)
}

