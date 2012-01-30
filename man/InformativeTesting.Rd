\name{InformativeTesting}
\alias{InformativeTesting}
\title{Testing Inequality Constrained Hypotheses in SEM}
\description{Testing inequality constrained Hypotheses in SEM}
\usage{
InformativeTesting(model = NULL, data, constraints = NULL, R = 1000L, 
                   type = "bollen.stine", return.LRT = TRUE, calibrate = FALSE, 
                   calibrate.R = 500L, calibrate.alpha = 0.05, 
                   parallel = c("no", "multicore", "snow"), ncpus = 1L, 
                   cl = NULL, verbose = FALSE, stoptest = 1.0, 
                   conclusion = FALSE, ...)
}
\arguments{
\item{model}{Model syntax specifying the model. See \code{\link{model.syntax}} 
    for more information.}
\item{data}{The data frame containing the observed variables being used to 
    fit the model.}
\item{constraints}{The imposed inequality constraints on the model.}
\item{R}{Number of bootstrap draws. Default value is set to 1000.}
\item{type}{If \code{"parametric"}, the parametric bootstrap is used. 
    If \code{"bollen.stine"}, the semi-nonparametric Bollen-Stine boostrap 
    is used. The default is set \code{"bollen.stine"}.}
\item{return.LRT}{If \code{TRUE}, the function returns bootstrapped LRT-values 
    from the observed data.}
\item{calibrate}{If \code{TRUE}, a double bootstrap is performed for adjusting 
    the alpha level.}
\item{calibrate.R}{Number of double bootstrap draws. Only used if 
    \code{calibrate=TRUE}. The default value is set to 500.}
\item{calibrate.alpha}{Sets the x'th quantile of generated plug-in p-value. 
    Default is the 5th quantile (0.05). Only used if \code{calibrate=TRUE}}
\item{parallel}{The type of parallel operation to be used (if any). If missing,
    the default is set \code{"no"}. The value \code{"multicore"} is for 
    Linux/Mac systems and the value \code{"snow"} for Windows systems. 
    If the value is \code{"snow"} the snow package needs to be installed.}
\item{ncpus}{The number of CPU's. Normally one chose this to the number of 
    available CPU's - 1'. The default value is set to 1.} 
\item{cl}{cluster object. For more information see the manual page of the 
    snow package.}
\item{verbose}{If \code{TRUE}, information is shown at each bootstrap draw.}
\item{stoptest}{The procedure stops when the plug-in p-value for Type A is 
    larger than the pre-specified value. The default value is set to 1.}
\item{conclusion}{If \code{TRUE}, a conclusion in words is printed.}
\item{...}{Other named arguments from the lavaan package which are passed to 
    the function. For example \code{group} for a multiple group model.}
}
%\details{}
\value{An object of class InformativeTesting for which a summary method is 
available.}   
%\seealso{}
\examples{
# Multiple regression 

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
                              type = "bollen.stine", return.LRT = FALSE, 
                              calibrate = TRUE, calibrate.R = 2,                       
                              calibrate.alpha = 0.05)
example
summary(example)
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

