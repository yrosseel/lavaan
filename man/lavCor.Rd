\name{lavCor}
\alias{lavCor}
\title{Polychoric, polyserial and Pearson correlations}
\description{
Fit an unrestricted model to compute polychoric, polyserial and/or Pearson
correlations.}
\usage{
lavCor(object, ordered = NULL, group = NULL, missing = "listwise", 
       ov.names.x = NULL, se = "none", estimator = "two.step", ..., 
       output = "cor")
}
\arguments{
\item{object}{Either a \code{data.frame}, or an object of class
\code{\linkS4class{lavaan}}. If the input is a \code{data.frame}, and
    some variables are declared as ordered factors, lavaan will
    treat them as ordinal variables.}
\item{ordered}{Character vector. Only used if \code{object} is a 
    \code{data.frame}. Treat these variables as ordered (ordinal) variables. 
    Importantly, all other variables will be treated as numeric (unless they 
    are declared as ordered in the original data frame.)}
\item{group}{Only used if \code{object} is a \code{data.frame}. Specify
    a grouping variable.}
\item{missing}{If \code{"listwise"}, cases with missing values are removed
    listwise from the data frame. If \code{"direct"} or
    \code{"ml"} or \code{"fiml"} and the estimator is maximum likelihood,
    an EM algorithm is used to estimate the unrestricted covariance matrix
    (and mean vector). If \code{"pairwise"}, pairwise deletion is used. If
    \code{"default"}, the value is set depending on the estimator and the
    mimic option.}
\item{ov.names.x}{Only used if \code{object} is a \code{data.frame}. Specify
    variables that need to be treated as exogenous. Only used if at least
    one variable is declared as ordered.}
\item{se}{Only used if \code{output} (see below) contains standard errors.
    See the \code{\link{lavaan}} function for possible options.}
\item{estimator}{If \code{"none"} or \code{"two.step"} or \code{"two.stage"},
   only starting values are computed for the correlations (and thresholds),
   without any further estimation. If all variables are continuous, the 
   starting values are the sample covariances (converted to correlations
   if \code{output = "cor"}). If at least one variable is ordered, 
   the thresholds are computed using univariate information only. The
   polychoric and/or polyserial correlations are computed in a second
   stage, keeping the values of the thresholds constant. If an estimator
   (other than \code{"two.step"} or \code{"two.stage"}) is specified
   (for example \code{estimator = "PML"}), these starting values are
   further updated by fitting the unrestricted model using the chosen
   estimator. See the \code{\link{lavaan}} function for alternative 
   estimators.}
\item{...}{Optional parameters that are passed to the \code{\link{lavaan}} 
    function.}
\item{output}{If \code{"cor"}, the function returns the correlation matrix only.   If \code{"cov"}, the function returns the covariance matrix (this only makes
   a difference if at least one variable is numeric). If \code{"th"} or
   \code{"thresholds"}, only the thresholds are returned. If \code{"sampstat"},
   the output equals the result of \code{inspect(fit, "sampstat")} where
   fit is the unrestricted model. If \code{"est"} or \code{"pe"} or 
   \code{"parameterEstimates"}, the output equals the result of
   \code{parameterEstimates(fit)}. Finally, if output is \code{"fit"} or
   \code{"lavaan"}, the function returns an object of class 
   \code{\linkS4class{lavaan}}.}
}
\details{
    This function is a wrapper around the \code{\link{lavaan}} function,
    but where the model is defined as the unrestricted model. The
    following free parameters are included: all covariances/correlations among
    the variables, variances for continuous variables, means for continuous
    variables, thresholds for ordered variables, and if exogenous variables 
    are included (\code{ov.names.x} is not empty) while some variables
    are ordered, also the regression slopes enter the model.
}
\value{
    By default, if \code{output = "cor"} or \code{output = "cov"}, a symmetric
    matrix (of class \code{"lavaan.matrix.symmetric"}, which only affects the
    way the matrix is printed). If \code{output = "th"}, a named vector of
    thresholds. If \code{output = "fit"} or \code{output = "lavaan"}, 
    an object of class  \code{\linkS4class{lavaan}}. 
}
\references{
Olsson, U. (1979). Maximum likelihood estimation of the polychoric correlation 
coefficient. Psychometrika, 44(4), 443-460.

Olsson, U., Drasgow, F., & Dorans, N. J. (1982). The polyserial correlation coefficient. Psychometrika, 47(3), 337-347.
}
\seealso{
\code{\link{lavaan}}
}
\examples{
# Holzinger and Swineford (1939) example
HS9 <- HolzingerSwineford1939[,c("x1","x2","x3","x4","x5",
                                 "x6","x7","x8","x9")]

# Pearson correlations
lavCor(HS9)

# ordinal version, with three categories
HS9ord <- as.data.frame( lapply(HS9, cut, 3, labels=FALSE) )

# polychoric correlations, two-stage estimation
lavCor(HS9ord, ordered=names(HS9ord))

# thresholds only
lavCor(HS9ord, ordered=names(HS9ord), output = "th")

# polychoric correlations, with standard errors
lavCor(HS9ord, ordered=names(HS9ord), se = "standard", output="est")

# polychoric correlations, full output
fit.un <- lavCor(HS9ord, ordered=names(HS9ord), se = "standard", output="fit")
summary(fit.un)
}
