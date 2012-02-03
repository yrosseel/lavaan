\name{InformativeTesting}
\alias{InformativeTesting}
\title{Testing Inequality Constrained Hypotheses in SEM}
\description{Testing inequality constrained Hypotheses in SEM}
\usage{
InformativeTesting(model = NULL, data, constraints = NULL, R = 1000L, 
                   type = "bollen.stine", return.LRT = TRUE, calibrate = FALSE, 
                   calibrate.R = 500L, calibrate.alpha = 0.05, 
                   parallel = c("no", "multicore", "snow"), ncpus = 1L, 
                   cl = NULL, verbose = FALSE, stoptest = NULL, 
                   conclusion = FALSE, ...)
}
\arguments{
\item{model}{Model syntax specifying the model. See \code{\link{model.syntax}} 
    for more information.}
\item{data}{The data frame containing the observed variables being used to 
    fit the model.}
\item{constraints}{The imposed inequality constraints on the model.}
\item{R}{Integer; number of bootstrap draws. Default value is set to 1000.}
\item{type}{If \code{"parametric"}, the parametric bootstrap is used. 
    If \code{"bollen.stine"}, the semi-nonparametric Bollen-Stine boostrap 
    is used. The default value is set to \code{"bollen.stine"}.}
\item{return.LRT}{Logical; if \code{TRUE}, the function returns bootstrapped         
    LRT-values.}
\item{calibrate}{Logical; if \code{TRUE}, a double bootstrap is used to compute 
    an additional set of plug-in p-values for each boostrap sample.}
\item{calibrate.R}{Integer; number of double bootstrap draws. Only used if 
    \code{calibrate=TRUE}. The default value is set to 500.}
\item{calibrate.alpha}{The significance level to compute the adjusted alpha 
    based on the plugin p-values. Only used if \code{calibrate=TRUE}}
\item{parallel}{The type of parallel operation to be used (if any). If missing, 
    the default is "no".}
\item{ncpus}{Integer: number of processes to be used in parallel operation: 
    typically one would chose this to the number of available CPUs.} 
\item{cl}{An optional parallel or snow cluster for use if 
    \code{parallel = "snow"}. If not supplied, a cluster on the local machine 
    is created for the duration of the \code{InformativeTesting} call.}
\item{verbose}{Logical; if \code{TRUE}, information is shown at each bootstrap         
    draw.}
\item{stoptest}{The \code{InformativeTesting()} stops when the plug-in p-value 
    for the null-hypothesis for Type A is larger than the pre-specified value.}
\item{conclusion}{Logical; if \code{TRUE}, a conclusion in words is printed.}
\item{...}{Other named arguments from the lavaan package which are passed to 
    the function. For example \code{"group"} for a multiple group model.}
}
%\details{}
\value{An object of class InformativeTesting for which a summary and a plot 
method is available.}   
%\seealso{}
\examples{
#Multiple regression 

model  <- '
          y1 ~ b1*x1 + b2*x2 + b3*x3 +
          start(1.8)*x1 + start(0.017)*x2 + start(-0.31)*x3
         
          y1~~start(5.7)*y1
          '
      
constraints  <-'
              b2 < b1
              b3 < b1
              '


#We only generate 5 bootstrap samples and 2 double bootstraps; 
#in practice you may wish to use a much higher number.

example <- InformativeTesting(model = model, data = PoliticalDemocracy, 
                              constraints = constraints, R = 5, 
                              type = "bollen.stine",
                              calibrate = TRUE, calibrate.R = 2,                       
                              calibrate.alpha = 0.05)
example
summary(example)
#plot(example)
}

\references{
Van de Schoot, R., Hoijtink, H., & Dekovic, M. (2010). 
Testing inequality constrained hypotheses in SEM models. 
Structural Equation Modeling, 17, 443-463.

Van de Schoot, R., Strohmeier, D. (2011). 
Testing informative hypotheses in SEM increases power: 
An illustration contrasting classical. 
International Journal of Behavioral Development 35(2), 180-190.
} 

\author{
Leonard Vanbrabant
\email{l.g.f.vanbrabant@uu.nl}
}

