\name{lavNames}
\alias{lavNames}
\alias{lavaanNames}
\title{lavaan Names}
\description{
Extract variables names from a fitted lavaan object.}
\usage{
lavNames(object, type = "ov", ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{Character. The type of variables whose names should be extracted. 
See details for a complete list.}
\item{...}{Additional selection variables. For example \code{"group = 2L"}
(in a multiple-group analysis) only considers the variables included
in the model for the second group.}
}
\details{
The order of the variable names, as returned by \code{lavNames} 
determines the order in which the variables are listed in the parameter
table, and therefore also in the summary output.


The following variable types are available:
\itemize{
  \item \code{"ov"}: observed variables

  \item \code{"ov.x"}: (pure) exogenous observed variables (no mediators)

  \item \code{"ov.nox"}: non-exogenous observed variables

  \item \code{"ov.model"}: modelled observed variables (joint vs conditional)

  \item \code{"ov.y"}: (pure) endogenous variables (dependent only) (no mediators)

  \item \code{"ov.num"}: numeric observed variables

  \item \code{"ov.ord"}: ordinal observed variables

  \item \code{"ov.ind"}: observed indicators of latent variables

  \item \code{"ov.orphan"}: lonely observed variables (only intercepts/variancesappear in the model syntax)

  \item \code{"ov.interaction"}: interaction terms (defined by the colon operator)

  \item \code{"th"}: threshold names ordinal variables only

  \item \code{"th.mean"}: threshold names ordinal + numeric variables (if any)

  \item \code{"lv"}: latent variables

  \item \code{"lv.regular"}: latent variables (defined by =~ only)

  \item \code{"lv.formative"}: latent variables (defined by <~ only)

  \item \code{"lv.x"}: (pure) exogenous variables

  \item \code{"lv.y"}: (pure) endogenous variables

  \item \code{"lv.nox"}: non-exogenous latent variables

  \item \code{"lv.nonnormal"}: latent variables with non-normal indicators

  \item \code{"lv.interaction"}:  interaction terms at the latent level

  \item \code{"eqs.y"}: variables that appear as dependent variables in a
                        regression formula (but not indicators of latent
                        variables)
 
  \item \code{"eqs.x"}: variables that appear as independent variables in
                        a regression formula
}
}
\seealso{\code{\link{lavaanify}}, \code{\link{parTable}}}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
lavNames(fit, "ov")
}
