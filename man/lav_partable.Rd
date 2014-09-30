\name{lav_partable}
\alias{lav_partable_independence}
\alias{lav_partable_unrestricted}
\alias{lav_partable_df}
\alias{lav_partable_ndat}
\alias{lav_partable_npar}
\title{lavaan partable functions}
\description{Utility functions related to the parameter table (partable)}
\usage{
# extract information from a parameter table
lav_partable_df(partable, group = NULL)
lav_partable_ndat(partable, group = NULL)
lav_partable_npar(partable)

# generate parameter table for specific models
lav_partable_independence(lavobject = NULL, 
    ov.names = NULL, ov = NULL, ov.names.x = NULL, 
    sample.cov = NULL, meanstructure = FALSE, sample.mean = NULL, 
    sample.th = NULL, parameterization = "delta", fixed.x = TRUE)

lav_partable_unrestricted(lavobject = NULL,
    ov.names = NULL, ov = NULL, ov.names.x = NULL, 
    sample.cov = NULL, meanstructure = FALSE, sample.mean = NULL, 
    sample.th = NULL, fixed.x = TRUE)
}
\arguments{
\item{partable}{A parameter table. see \code{\link{lavParTable}} 
for more information.}
\item{group}{Integer. If non-null, only consider this group.}
\item{lavobject}{An object of class `lavaan'. If this argument is provided,
it should be the only argument. All the values for the other arguments 
are extracted from this object.}
\item{ov.names}{List of character vectors. Each list element contains the 
names of the observed variables for this group.}
\item{ov}{Variable table. See \code{\link{varTable}} for more information.}
\item{ov.names.x}{List of character vectors. Each list element contains 
the names of exogenous observed variables for this group.}
\item{sample.cov}{Optional list of numeric matrices. 
Each list element contains a sample variance-covariance matrix for this group.
If provided, these values will be used as starting values.}
\item{meanstructure}{Logical. If \code{TRUE}, the means of the observed
    variables enter the model.}
\item{sample.mean}{Optional list of numeric vectors. 
Each list element contains a sample mean vector for this group. 
If provided, these values will be used as starting values.}
\item{sample.th}{Optional list of numeric vectors. 
Each list element contains a vector of sample thresholds for this group.
If provided, these values will be used as starting values.}
\item{parameterization}{Currently only used if data is categorical. If
    \code{"delta"}, the delta parameterization is used. If \code{"theta"},
    the theta parameterization is used.}
\item{fixed.x}{If \code{TRUE}, the exogenous `x' covariates are considered
    fixed variables.}
}
\examples{
# generate independence model for 3 variables with two groups
lav <- lav_partable_independence(ov.names = list( g1 = c("x1","x2","x3"),
                                                  g2 = c("x1","x2","x3") ))

# how many free parameters?
lav_partable_npar(lav)

# how many sample statistics?
lav_partable_ndat(lav)

# how many degrees of freedom?
lav_partable_df(lav)
}

