\name{sam}
\alias{sam}
\alias{fsr}
\title{Fit Structural Equation Models using the SAM approach}
\description{
Fit a Structural Equation Model (SEM) using the Structural After Measurement
(SAM) approach.}
\usage{
sam(model = NULL, data = NULL, cmd = "sem", se = "twostep", 
    mm.list = NULL, mm.args = list(bounds = "wide.zerovar"), 
    struc.args = list(estimator = "ML"), 
    sam.method = "local", alpha.correction = 0L, ..., 
    local.options = list(M.method = "ML", veta.force.pd = TRUE, 
    twolevel.method = "h1"), global.options = list(), output = "lavaan")
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See            
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (eg. the output of the \code{lavaanify()} function) is also
    accepted.}
\item{data}{A data frame containing the observed variables used in
    the model.}
\item{cmd}{Character. Which command is used to run the sem models. The possible
    choices are \code{"sem"}, \code{"cfa"} or \code{"lavaan"}, determining how
    we deal with default options.}
\item{se}{Character. The type of standard errors that are used in the
    final (structural) model. If \code{"twostep"} (the default), the
    standard errors take the estimation uncertainty of the first (measurement) 
    stage into account. If \code{"standard"}, this uncertainty is ignored,
    and we treat the measurement information as known. If \code{"none"},
    no standard errors are computed.}
\item{mm.list}{List. Define the measurement blocks. 
     Each element of the list should
     be either a single name of a latent variable, or a vector of latent
     variable names.  If omitted, a separate measurement
     block is used for each latent variable.}
\item{mm.args}{List. Optional arguments for the fitting
    function(s) of the measurement block(s) only. See \code{\link{lavOptions}} 
    for a complete list.}
\item{struc.args}{List. Optional arguments for the fitting function of the
    structural part only. See \code{\link{lavOptions}} for a complete list.}
\item{sam.method}{Character. Can be set to \code{"local"}, \code{"global"}
    or \code{"fsr"}. In the latter case, the results are the same as if
    Bartlett factor scores were used, without any bias correction.}
\item{alpha.correction}{Integer. Acceptable values are in the range 0 till
    N-1. If zero (the default), no small sample correction is performed, 
    and the bias-correction is the same as with local SAM. When equal to N-1, 
    the bias-correction is eliminated, and the results are the same as 
    naive FSR. Typical values are 0, P+1 (where P is the number of
    predictors in the structural model), P+5, and (N-1)/2.}
\item{...}{Many more additional options can be defined, using 'name = value'.
    See \code{\link{lavOptions}} for a complete list. These options affect
    both the measurement blocks and the structural part.}
\item{local.options}{List. Options specific for local SAM method (these
     options may change over time).}
\item{global.options}{List. Options specific for global SAM method (not
     used for now).}
\item{output}{Character. If \code{"lavaan"}, a lavaan object returned. If
    \code{"list"}, a list is returned with all the ingredients from the 
    different stages.}
}
\details{
    The \code{sam} function tries to automate the SAM approach, by first
    estimating the measurement part of the model,
    and then the structural part of the model. See reference for more
    details.
}
\value{
  If \code{output = "lavaan"}, an object of class 
  \code{\linkS4class{lavaan}}, for which several methods
  are available, including a \code{summary} method. If \code{output = "list"},
  a list.
}
\references{
Rosseel and Loh (2021). A structural-after-measurement approach to Structural
Equation Modeling. Psychological Methods. Advance online publication.
https://dx.doi.org/10.1037/met0000503
}
\seealso{
\code{\link{lavaan}}
}
\examples{
## The industrialization and Political Democracy Example 
## Bollen (1989), page 332
model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit.sam <- sam(model, data = PoliticalDemocracy,
               mm.list = list(ind = "ind60", dem = c("dem60", "dem65")))
summary(fit.sam)
}
