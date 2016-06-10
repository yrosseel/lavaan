\name{InformativeTesting}
\alias{InformativeTesting}
\alias{informativetesting}
\title{Testing order/inequality Constrained Hypotheses in SEM}
\description{Testing order/inequality constrained Hypotheses in SEM}
\usage{
InformativeTesting(model = NULL, data, constraints = NULL, 
                   R = 1000L, type = "bollen.stine",
                   return.LRT = TRUE, 
                   double.bootstrap = "standard",
                   double.bootstrap.R = 249L, 
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
    value is set to 249.}
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
\details{The following hypothesis tests are available:
\itemize{
  \item Type A: Test H0: all restriktions with equalities ("=") 
  active against HA: at least one inequality restriktion (">") 
  strictly true.
  \item Type B: Test H0: all restriktions with inequalities (">") 
  (including some equalities ("=")) active against HA: at least 
  one restriktion false (some equality restriktions may be 
  maintained).}
}
\value{An object of class InformativeTesting for which a 
print and a plot method is available.}   
%\seealso{}
\examples{
\dontrun{
#########################
### real data example ###
#########################
# Multiple group path model for facial burns example.

# model syntax with starting values.
  burns.model <- 'Selfesteem ~ Age + c(m1, f1)*TBSA + HADS +
                             start(-.10, -.20)*TBSA  
                 HADS ~ Age + c(m2, f2)*TBSA + RUM +
                        start(.10, .20)*TBSA '
 
 
# constraints syntax
 burns.constraints <- 'f2 > 0  ; m1 < 0
                       m2 > 0  ; f1 < 0
                       f2 > m2 ; f1 < m1'
 
# we only generate 2 bootstrap samples in this example; in practice
# you may wish to use a much higher number. 
# the double bootstrap was switched off; in practice you probably 
# want to set it to "standard".
example1 <- InformativeTesting(model = burns.model, data = FacialBurns,
                               R = 2, constraints = burns.constraints,
                               double.bootstrap = "no", group = "Sex")

example1

##########################
### artificial example ###
##########################
# Simple ANOVA model with 3 groups (N = 20 per group)
set.seed(1234)
Y <- cbind(c(rnorm(20,0,1), rnorm(20,0.5,1), rnorm(20,1,1)))
grp <- c(rep("1", 20), rep("2", 20), rep("3", 20))
Data <- data.frame(Y, grp)

#create model matrix
fit.lm <- lm(Y ~ grp, data = Data)
mfit <- fit.lm$model
mm <- model.matrix(mfit)

Y <- model.response(mfit)
X <- data.frame(mm[,2:3])
names(X) <- c("d1", "d2")
Data.new <- data.frame(Y, X)

# model
model <- 'Y ~ 1 + a1*d1 + a2*d2'

# fit without constraints
fit <- sem(model, data = Data.new)

# constraints syntax: mu1 < mu2 < mu3
constraints <- ' a1 > 0
                 a1 < a2 '

# we only generate 10 bootstrap samples in this example; in practice
# you may wish to use a much higher number, say > 1000. The double 
# bootstrap is not necessary in case of an univariate ANOVA model.
example2 <- InformativeTesting(model = model, data = Data.new, 
                               start = parTable(fit),
                               R = 10L, double.bootstrap = "no",
                               constraints = constraints)
example2
}
}

\references{
Van de Schoot, R., Hoijtink, H., & Dekovic, M. (2010). 
Testing inequality constrained hypotheses in SEM models. 
\emph{Structural Equation Modeling}, \bold{17}, 443-463.

Van de Schoot, R., Strohmeier, D. (2011). 
Testing informative hypotheses in SEM increases power: An 
illustration contrasting classical. \emph{International Journal 
of Behavioral Development}, \bold{35}, 180-190.

Silvapulle, M.J. and Sen, P.K. (2005). \emph{Constrained 
Statistical Inference}. Wiley, New York.

} 

\author{
Leonard Vanbrabant
\email{lgf.vanbrabant@gmail.com}
}
