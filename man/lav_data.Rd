\name{lav_data}
\alias{lav_data_update}
\title{lavaan data functions}
\description{Utility functions related to the Data slot}
\usage{
# update data slot with new data (of the same size)
lav_data_update(lavdata = NULL, newX = NULL, BOOT.idx = NULL, lavoptions = NULL)
}
\arguments{
\item{lavdata}{A lavdata object.}
\item{newX}{A list of (new) data matrices (per group) of the same size. 
They will replace the data stored in the internal dataslot.}
\item{BOOT.idx}{A list of integers. If bootstrapping was used to produce
the data in newX, use these indices to adapt the remaining slots.}
\item{lavoptions}{A named list. The Options lsot from a lavaan object.}
}
\examples{
# generate syntax for an independence model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)

# extract data slot and options
lavdata <- fit@Data
lavoptions <- lavInspect(fit, "options")

# create bootstrap sample
boot.idx <- sample(x = nobs(fit), size = nobs(fit), replace = TRUE)
newX <- list(lavdata@X[[1]][boot.idx,])

# generate update lavdata object
newdata <- lav_data_update(lavdata = lavdata, newX = newX, 
                           lavoptions = lavoptions)
}

