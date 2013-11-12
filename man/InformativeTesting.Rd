\name{InformativeTesting}
\alias{InformativeTesting}
\alias{informativetesting}
\alias{print.InformativeTesting}
\title{Testing order Constrained Hypotheses in SEM}
\description{Testing order constrained Hypotheses in SEM}
\usage{
InformativeTesting(model = NULL, data, constraints = NULL, 
                   R = 1000L, type = "bollen.stine",
                   return.LRT = TRUE, 
                   double.bootstrap = "standard",
                   double.bootstrap.R = 500L, 
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
\item{R}{Integer; number of bootstrap draws. The default value is set to 1000.}
\item{type}{If \code{"parametric"}, the parametric bootstrap is used. 
    If \code{"bollen.stine"}, the semi-nonparametric Bollen-Stine bootstrap 
    is used. The default is set to \code{"bollen.stine"}.}
\item{return.LRT}{Logical; if \code{TRUE}, the function returns bootstrapped         
    LRT-values.}
\item{double.bootstrap}{If \code{"standard"} (default) the genuine double bootstrap is 
    used to compute an additional set of plug-in p-values for each bootstrap       
    sample. If \code{"no"}, no double bootstrap is used. If \code{"FDB"}, 
    the fast double bootstrap is used to compute second level LRT-values for 
    each bootstrap sample. Note that the \code{"FDB"} is experimental and should 
    not be used by inexperienced users.}    
\item{double.bootstrap.R}{Integer; number of double bootstrap draws. The default 
    value is set to 500.}
\item{double.bootstrap.alpha}{The significance level to compute the adjusted 
    alpha based on the plugin p-values. Only used if \code{double.bootstrap = "standard"}. 
    The default value is set to 0.05.}
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
    the function. For example \code{"group"} in a multiple group model.}
}
%\details{}
\value{An object of class InformativeTesting for which a plot 
    method is available.}   
%\seealso{}
\examples{
#path model facial burns example

#model syntax with starting values
burns.model <- 'Selfesteem ~ Age + c(m1, f1)*TBSA + HADS +
                           start(-.15, -.25)*TBSA + start(-0.5, -0.5)*HADS 
                HADS ~ Age + c(m2, f2)*TBSA + RUM +
                       start(.15, .25)*TBSA + start(1.0, -0.6)*RUM
                # optional starting values only to get < 5s
                Selfesteem ~ start(37, 35)*1
                HADS       ~ start(-4, 17)*1
                Selfesteem ~~ c(18, 14)*Selfesteem
                HADS       ~~ c(31, 50)*HADS
'



#constraints syntax
burns.constraints <- 'f2 > 0  ; m1 < 0
                      m2 > 0  ; f1 < 0
                      f2 > m2 ; f1 < m1'


# we only generate 1 (double) bootstrap samples in this example; in practice
# you may wish to use a much higher number

example <- InformativeTesting(model = burns.model, data = FacialBurns,
                              R = 1, double.bootstrap.R = 1, 
                              constraints = burns.constraints, group = "Sex")

# dontrun
# plot(example)
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
