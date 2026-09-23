\name{lavNames}
\alias{lavNames}
\alias{lavaanNames}
\title{lavaan Names}
\description{
Extract variable names from a fitted lavaan object.}
\usage{
lavNames(object, type = "ov", ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{type}{Character. The type of variables whose names should be extracted.
See details for a complete list. Dots and underscores are interchangeable
(\code{"ov_nox"} and \code{"ov.nox"} select the same variables).}
\item{...}{Additional selection criteria. For example, \code{"group = 2L"}
(in a multiple-group analysis) restricts the selection to the variables
included in the model for the second group.}
}
\details{
The order of the variable names, as returned by \code{lav_object_vnames},
determines the order in which the variables are listed in the parameter
table, and therefore also in the summary output.


The following variable types are available:
\itemize{
  \item \code{"ov"}: observed variables

  \item \code{"ov_x"}: (pure) exogenous observed variables (no mediators)

  \item \code{"ov_nox"}: non-exogenous observed variables

  \item \code{"ov_model"}: modeled observed variables (joint vs conditional)

  \item \code{"ov_y"}: (pure) endogenous variables (dependent only) (no mediators)

  \item \code{"ov_num"}: numeric observed variables

  \item \code{"ov_ord"}: ordinal observed variables

  \item \code{"ov_ind"}: observed indicators of latent variables

  \item \code{"ov_orphan"}: isolated observed variables (only their intercepts or variances appear in the model syntax)

  \item \code{"ov_interaction"}: interaction terms (defined by the colon operator)

  \item \code{"th"}: threshold names (ordinal variables only)

  \item \code{"th_mean"}: threshold names (ordinal and numeric variables, if any)

  \item \code{"lv"}: latent variables

  \item \code{"lv_regular"}: measured latent variables (defined by =~ only)

  \item \code{"lv_composite"}: composites (defined by <~ only)

  \item \code{"lv_x"}: (pure) exogenous variables

  \item \code{"lv_y"}: (pure) endogenous variables

  \item \code{"lv_nox"}: non-exogenous latent variables

  \item \code{"lv_nonnormal"}: latent variables with non-normal indicators

  \item \code{"lv_interaction"}:  interaction terms at the latent level

  \item \code{"eqs_y"}: variables that appear as dependent variables in a
                        regression formula (but not indicators of latent
                        variables)
 
  \item \code{"eqs_x"}: variables that appear as independent variables in
                        a regression formula
}
}
\seealso{\code{\link{lavParTable}}, \code{\link{parTable}}}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
lavNames(fit, "ov")
}
