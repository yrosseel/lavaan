\name{lav_partable}
\alias{lav_partable_independence}
\alias{lav_partable_unrestricted}
\alias{lav_partable_df}
\alias{lav_partable_ndat}
\alias{lav_partable_npar}
\alias{lav_partable_labels}
\alias{lav_partable_from_lm}
\title{lavaan partable functions}
\description{Utility functions related to the parameter table (partable)}
\usage{
# extract information from a parameter table
lav_partable_df(partable, group = NULL)
lav_partable_ndat(partable, group = NULL)
lav_partable_npar(partable)

# generate parameter labels
lav_partable_labels(partable, group.equal = "", group.partial = "", 
                    type = "user")

# generate parameter table for specific models
lav_partable_independence(lavobject = NULL, lavdata = NULL, lavoptions = NULL, 
    lavsamplestats = NULL, sample.cov = NULL, sample.mean = NULL, 
    sample.th = NULL, sample.th.idx = NULL)

lav_partable_unrestricted(lavobject = NULL, lavdata = NULL, lavoptions = NULL, 
    lavsamplestats = NULL, sample.cov = NULL, sample.mean = NULL, 
    sample.th = NULL, sample.th.idx = NULL)

lav_partable_from_lm(object, est = FALSE, label = FALSE, as.data.frame. = FALSE)
}
\arguments{
\item{partable}{A parameter table. see \code{\link{lavParTable}} 
for more information.}
\item{group}{Integer. If non-null, only consider this group.}
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
\item{lavoptions}{A names list. The Options lsot from a lavaan object.}
\item{lavsamplestats}{An object of class `lavSampleStats'. The SampleStats
slot from a  lavaan object.}
\item{sample.cov}{Optional list of numeric matrices. 
Each list element contains a sample variance-covariance matrix for this group.
If provided, these values will be used as starting values.}
\item{sample.mean}{Optional list of numeric vectors. 
Each list element contains a sample mean vector for this group. 
If provided, these values will be used as starting values.}
\item{sample.th}{Optional list of numeric vectors. 
Each list element contains a vector of sample thresholds for this group.
If provided (and also sample.th.idx is provided), 
these values will be used as starting values.}
\item{sample.th.idx}{Optional list of integers. Each list contains the
threshold indices for this group.}
\item{est}{Logical. If TRUE, include the fitted estimates in the parameter table.}
\item{label}{Logical. If TRUE, include parameter labels in the parameter table.}
\item{as.data.frame.}{Logical. If TRUE, return the parameter table as a data.frame.}
\item{object}{An object of class \code{lm}.}
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

# how many degrees of freedom?
lav_partable_df(lav)
}

