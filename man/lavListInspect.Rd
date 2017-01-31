\name{lavListInspect}
\alias{lavListInspect}
\alias{lavListTech}
\title{Inspect or extract information from a lavaanList object}
\description{
The \code{lavListInspect()} and \code{lavListTech()} functions can be used to
inspect/extract information that is stored inside (or can be computed from) a
lavaanList object.
}
\usage{
lavListInspect(object, what = "free", add.labels = TRUE, 
               add.class = TRUE, list.by.group = TRUE,
               drop.list.single.group = TRUE)

lavListTech(object, what = "free", add.labels = FALSE, 
            add.class = FALSE, list.by.group = FALSE,
            drop.list.single.group = FALSE)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaanList}}.}
\item{what}{Character. What needs to be inspected/extracted? See Details for a
full list. Note: the \code{what} argument is not case-sensitive
(everything is converted to lower case.)}
\item{add.labels}{If \code{TRUE}, variable names are added to the vectors
and/or matrices.}
\item{add.class}{If \code{TRUE}, vectors are given the `lavaan.vector' class;
matrices are given the `lavaan.matrix' class, and symmetric matrices are
given the `lavaan.matrix.symmetric' class. This only affects the way they
are printed on the screen.}
\item{list.by.group}{Logical. Only used when the output are model matrices.
If \code{TRUE}, the model matrices are nested within groups. If \code{FALSE},
a flattened list is returned containing all model matrices, with repeated
names for multiple groups.}
\item{drop.list.single.group}{If \code{FALSE}, the results are returned as
    a list, where each element corresponds to a group (even if there is only
    a single group.) If \code{TRUE}, the list will be unlisted if there is
    only a single group.}
}
\details{
The \code{lavListInspect()} and \code{lavListTech()} functions only differ in
the way they return the results. The \code{lavListInspect()} function will
prettify the output by default, while the \code{lavListTech()} will not attempt
to prettify the output by default.

Below is a list of possible values for the \code{what} argument, organized
in several sections:

Model matrices:

\describe{
    \item{\code{"free"}:}{A list of model matrices. The non-zero integers
        represent the free parameters. The numbers themselves correspond
        to the position of the free parameter in the parameter vector.
        This determines the order of the model parameters in the output
        of for example \code{coef()} and \code{vcov()}.}
    \item{\code{"partable"}:}{A list of model matrices. The non-zero integers
        represent both the fixed parameters (for example, factor loadings
        fixed at 1.0), and the free parameters if we ignore any equality 
        constraints. They correspond with all entries (fixed or free) 
        in the parameter table. See \code{\link{parTable}}.}
    \item{\code{"start"}:}{A list of model matrices. The values represent
        the starting values for all model parameters. 
        Alias: \code{"starting.values"}.}
}

Information about the data (including missing patterns):

\describe{
    \item{\code{"group"}:}{A character string. The group variable in
        the data.frame (if any).}
    \item{\code{"ngroups"}:}{Integer. The number of groups.}
    \item{\code{"group.label"}:}{A character vector. The group labels.}
    \item{\code{"cluster"}:}{A character vector. The cluster variable(s) 
       in the data.frame (if any).}
    \item{\code{"ordered"}:}{A character vector. The ordered variables.}
    \item{\code{"nobs"}:}{Integer vector. The number of observations
        in each group that were used in the analysis (in each dataset).}
    \item{\code{"norig"}:}{Integer vector. The original number of observations
        in each group (in each dataset).}
    \item{\code{"ntotal"}:}{Integer. The total number of observations that
        were used in the analysis. If there is just a single group, this
        is the same as the \code{"nobs"} option; if there are multiple groups,
        this is the sum of the \code{"nobs"} numbers for each group
        (in each dataset).}
}

Model features:

\describe{
    \item{\code{"meanstructure"}:}{Logical. \code{TRUE} if a meanstructure
        was included in the model.}
    \item{\code{"categorical"}:}{Logical. \code{TRUE} if categorical endogenous
        variables were part of the model.}
    \item{\code{"fixed.x"}:}{Logical. \code{TRUE} if the exogenous x-covariates
        are treated as fixed.}
    \item{\code{"parameterization"}:}{Character. Either \code{"delta"} or
        \code{"theta"}.}
}

\describe{
    \item{\code{"list"}:}{The parameter table. The same output as given
        by \code{parTable()}.}
    \item{\code{"options"}:}{List. The option list.}
    \item{\code{"call"}:}{List. The call as returned by match.call, coerced to 
        a list.}
}

}
\seealso{
\code{\link{lavaanList}}
}
\examples{
# fit model
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
     
# a data generating function
generateData <- function() simulateData(HS.model, sample.nobs = 100)
     
set.seed(1234)
fit <- semList(HS.model, dataFunction = generateData, ndat = 5,
               store.slots = "partable")

# extract information
lavListInspect(fit, "free")
lavListTech(fit, "free")
}
