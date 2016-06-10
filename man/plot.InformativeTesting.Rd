\name{InformativeTesting methods}
\alias{print.InformativeTesting}
\alias{plot.InformativeTesting}
\title{Methods for output InformativeTesting()}
\description{The print function shows the results of hypothesis tests
Type A and Type B. The plot function plots the distributions of 
bootstrapped LRT values and plug-in p-values.}

\usage{
    
\method{print}{InformativeTesting}(x, digits = max(3, getOption("digits") - 3), ...)
  
\method{plot}{InformativeTesting}(x, ..., type = c("lr","ppv"), 
   main = "main", xlab = "xlabel", ylab = "Frequency", freq = TRUE, 
   breaks = 15, cex.main = 1, cex.lab = 1, cex.axis = 1, 
   col = "grey", border = par("fg"), vline = TRUE, 
   vline.col = c("red", "blue"), lty = c(1,2), lwd = 1, 
   legend = TRUE, bty = "o", cex.legend = 1, loc.legend = "topright")
}
 
\arguments{
  \item{x}{object of class "InformativeTesting".}
  \item{digits}{the number of significant digits to use when printing.}
  \item{...}{Currently not used.}
  \item{type}{If \code{"lr"}, a distribution of the first-level 
    bootstrapped LR values is plotted. If \code{"ppv"} a distribution of 
    the bootstrapped plug-in p-values is plotted.}
  \item{main}{The main title(s) for the plot(s).}
  \item{xlab}{A label for the x axis, default depends on input type.}
  \item{ylab}{A label for the y axis.}
  \item{freq}{Logical; if TRUE, the histogram graphic is a representation of 
    frequencies, the counts component of the result; if \code{FALSE}, 
    probability densities, component density, are plotted 
    (so that the histogram has a total area of one). The default is 
    set to \code{TRUE}.}
  \item{breaks}{see \code{\link[graphics]{hist}}} 
  \item{cex.main}{The magnification to be used for main titles relative to the 
    current setting of cex.}
  \item{cex.lab}{The magnification to be used for x and y labels relative to 
    the current setting of cex.}
  \item{cex.axis}{The magnification to be used for axis annotation relative to         
    the current setting of cex.}  
  \item{col}{A colour to be used to fill the bars. The default of NULL 
    yields unfilled bars.}
  \item{border}{Color for rectangle border(s). The default means par("fg").}
  \item{vline}{Logical; if \code{TRUE} a vertical line is drawn 
    at the observed LRT value. If 
    \code{double.bootstrap = "FDB"} a vertical line is drawn at the 1-p* quantile 
    of the second-level LRT values, where p* is the first-level bootstrapped p-value}
  \item{vline.col}{Color(s) for the vline.LRT.}
  \item{lty}{The line type. Line types can either be specified as an integer 
    (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 
     5=longdash, 6=twodash) or as one of the character strings 
    "blank", "solid", "dashed", "dotted", "dotdash", "longdash", or 
    "twodash", where "blank" uses 'invisible lines' 
    (i.e., does not draw them).}
  \item{lwd}{The line width, a positive number, defaulting to 1.}
  \item{legend}{Logical; if \code{TRUE} a legend is added to the plot.}
  \item{bty}{A character string which determined the type of box which is drawn 
    about plots. If bty is one of "o" (the default), "l", "7", "c", "u", or "]" 
    the resulting box resembles the corresponding upper case letter. A value of 
    "n" suppresses the box.}
  \item{cex.legend}{A numerical value giving the amount by which the legend 
    text and symbols should be magnified relative to the default. 
    This starts as 1 when a device is opened, and is reset when the layout 
    is changed.}
  \item{loc.legend}{The location of the legend, specified by a single 
    keyword from the list \code{"bottomright"}, \code{"bottom"}, 
    \code{"bottomleft"}, \code{"left"}, \code{"topleft"}, \code{"top"}, 
    \code{"topright"}, \code{"right"} and \code{"center"}.}
}

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
plot(example1)

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
# plot(example2)
}
}

\author{
Leonard Vanbrabant
\email{lgf.vanbrabant@gmail.com}
}



