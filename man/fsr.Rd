\name{fsr}
\alias{fsr}
\title{Factor Score Regression}
\description{
Fit a SEM model using factor score regression.}
\usage{
fsr(model = NULL, data = NULL, cmd = "sem", 
    fsr.method = "Croon", fs.method = "Bartlett", 
    fs.scores = FALSE, Gamma.NT = TRUE, lvinfo = FALSE, ...,
    output = "fsr")
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax. See                 
    \code{\link{model.syntax}} for more information. Alternatively, a
    parameter table (eg. the output of the \code{lavaanify()} function) is also
    accepted.}
\item{data}{An optional data frame containing the observed variables used in
    the model. If some variables are declared as ordered factors, lavaan will
    treat them as ordinal variables.}
\item{cmd}{Character string. Which command is used to run the sem models. 
     The possible
    choices are \code{"sem"} or \code{"lavaan"}, determining how
    we deal with default options.}
\item{fsr.method}{Character. Factor score regression method. Possible
options are \code{naive}, \code{Skrondal-Laake}, and \code{Croon}.}
\item{fs.method}{Character. Factor score estimation method. Possible
options are \code{Bartlett} and \code{regression}.}
\item{fs.scores}{Logical. If \code{TRUE}, explicitly compute factor scores; if
\code{FALSE}, only compute the mean vector and variance matrix of the
factor scores.}
\item{Gamma.NT}{Logical. Only needed when \code{se="robust.sem"} and
data is missing; if \code{TRUE}, compute Gamma (N times the variance
matrix of the sample statistics) assuming normality.}
\item{lvinfo}{Logical. If \code{TRUE}, return latent variable information
as an attribute to the output.}
\item{...}{Further arguments that we pass to the \code{"cfa"}, \code{"sem"} or 
           \code{"lavaan"} functions.}
\item{output}{Character. Determines the type of output that is returned. If
\code{"fsr"}, an object of class \code{"fsr"} is returned. If
\code{"lavaan"}, a lavaan object (representing the structural part only) 
is returned.
If \code{"fs.scores"}, only the factor scores are returned. If \code{"fsr.cov"}, the (possibly corrected) covariance matrix of the factor scores is returned. 
If \code{"extra"},  a list with information from the measurement models 
is returned.}

}
\details{
    The \code{fsr} function implements a two-step procedure to estimate
    the parameters of the structural (regression) part of a SEM model.
    In a first step, factor scores are computed for each latent variable.
    In a second step, the latent variables are replaced by the factor
    scores, and a path analysis is used to estimate all remaining 
    model parameters. Special techniques are used in order to ensure
    (approximately) unbiased estimation of point estimates and 
    standard errors.
}
\value{
  An object of class \code{\linkS4class{lavaan}}, for which several methods
  are available, including a \code{summary} method.
}
\references{
Devlieger, I., Mayer, A., & Rosseel, Y. (2015). 
Hypothesis Testing Using Factor
Score Regression: A Comparison of Four Methods. Educational and Psychological
Measurement.
http://epm.sagepub.com/content/early/2015/09/29/0013164415607618.abstract
}
\seealso{
\code{\link{lavaan}}, \code{\link{sem}}, \code{\link{lavPredict}}
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
mm.list <- list(block1 = "ind60",
                block2 = c("dem60", "dem65"))
fit <- fsr(model, data = PoliticalDemocracy, fsr.method = "simple",
           mm.list = mm.list)
summary(fit)
}
