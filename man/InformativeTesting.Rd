\name{InformativeTesting}
\alias{InformativeTesting}
\alias{print.InformativeTesting}
\alias{summary.InformativeTesting}
\alias{plot.InformativeTesting}
\title{Testing order Constrained Hypotheses in SEM}
\description{Testing order constrained Hypotheses in SEM}
\usage{
InformativeTesting(model = NULL, data, constraints = NULL, 
                   R = 1000L, type = "bollen.stine",
                   return.LRT = TRUE, 
                   double.bootstrap = "FDB",
                   double.bootstrap.R = 1000L, 
                   double.bootstrap.alpha = 0.05,
                   parallel = c("no", "multicore", "snow"), 
                   ncpus = 1L, cl = NULL, verbose = FALSE, \dots)
}
\arguments{
\item{model}{Model syntax specifying the model. See \code{\link{model.syntax}} 
    for more information.}
\item{data}{The data frame containing the observed variables being used to 
    fit the model.}
\item{constraints}{The imposed inequality constraints on the model.}
\item{R}{Integer; number of bootstrap draws. Default value is set to 1000. 
    If \code{double.bootstrap = "FDB"} the default value is set to 2000.}
\item{type}{If \code{"parametric"}, the parametric bootstrap is used. 
    If \code{"bollen.stine"}, the semi-nonparametric Bollen-Stine boostrap 
    is used. The default is set to \code{"bollen.stine"}.}
\item{return.LRT}{Logical; if \code{TRUE}, the function returns bootstrapped         
    LRT-values.}
\item{double.bootstrap}{If \code{"standard"} the genuine double bootstrap is 
    used to compute an additional set of plug-in p-values for each boostrap       
    sample. If \code{"FDB"}, the fast double bootstrap is used to compute           
    second level LRT-values for each bootstrap sample. If \code{"no"}, no 
    double bootstrap is used. The default is set to \code{"FDB"}.}    
\item{double.bootstrap.R}{Integer; number of double bootstrap draws. Only used         
    if \code{double.bootstrap = "standard"}. The default value is set to 500.}
\item{double.bootstrap.alpha}{The significance level to compute the adjusted 
    alpha based on the plugin p-values. Only used if \code{double.bootstrap = 
    "standard"}. The default value is set to 0.05.}
\item{parallel}{The type of parallel operation to be used (if any). If missing, 
    the default is set "no".}
\item{ncpus}{Integer: number of processes to be used in parallel operation: 
    typically one would chose this to the number of available CPUs.} 
\item{cl}{An optional parallel or snow cluster for use if 
    \code{parallel = "snow"}. If not supplied, a cluster on the local machine 
    is created for the duration of the \code{InformativeTesting} call.}
\item{verbose}{Logical; if \code{TRUE}, information is shown at each bootstrap         
    draw.}
\item{...}{Other named arguments from the lavaan package which are passed to 
    the function. For example \code{"group"} for a multiple group model.}
}
%\details{}
\value{An object of class InformativeTesting for which a summary and a plot 
    method is available.}   
%\seealso{}
\examples{
#Multiple regression model

model  <- '
          y1 ~ b1*x1 + b2*x2 + b3*x3 +
          start(1.8)*x1 + start(0.017)*x2 + start(-0.31)*x3
         
          y1~~start(5.7)*y1
          '
      
constraints  <-'
              b2 < b1
              b3 < b1
              '


#We only generate 5 bootstrap samples; 
#in practice you may wish to use a much higher number.

example <- InformativeTesting(model = model, data = PoliticalDemocracy, 
                              constraints = constraints, R = 5) 
                              
example
#summary(example)

#The original estimates for the unconstrained model.
#summary(example, fit = "hu")

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
\email{l.g.f.vanbrabant@hotmail.com}
}

