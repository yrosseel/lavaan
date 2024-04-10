\name{lavTestLRT}
\alias{lavTestLRT}
\alias{lavtestLRT}
\alias{LRT}
\alias{lavLRTTest}
\alias{lavLRT}
\alias{anova}
\title{LRT test}
\description{
LRT test for comparing (nested) lavaan models.}
\usage{
lavTestLRT(object, ..., method = "default", test = "default",
           A.method = "delta", scaled.shifted = TRUE,
           type = "Chisq", model.names = NULL)
anova(object, ...)
}
\arguments{
\item{object}{An object of class \code{\linkS4class{lavaan}}.}
\item{...}{additional objects of class \code{\linkS4class{lavaan}}.}
\item{method}{Character string. The possible options are 
  \code{"satorra.bentler.2001"}, \code{"satorra.bentler.2010"},
  \code{"satorra.2000"}, and \code{"standard"}. See details.}
\item{test}{Character string specifying which scaled test statistics to use,
  in case multiple scaled \code{test=} options were requested when fitting the
  model(s). See details.}
\item{A.method}{Character string. The possible options are \code{"exact"}
  and \code{"delta"}. This is only used when method = \code{"satorra.2000"}.
  It determines how the Jacobian of the constraint function (the matrix A)
  will be computed. Note that if \code{A.method = "exact"}, the models must
  be nested in the parameter sense, while if \code{A.method = "delta"}, they
  only need to be nested in the covariance matrix sense.}
\item{scaled.shifted}{Logical. Only used when method = \code{"satorra.2000"}.
  If \code{TRUE}, we use a scaled and shifted test statistic; if \code{FALSE},
  we use a mean and variance adjusted (Satterthwaite style) test statistic.}
\item{type}{Character. If \code{"Chisq"}, the test statistic for each
  model is the (scaled or unscaled) model fit test statistic. If \code{"Cf"}, 
  the test statistic for each model is computed by the 
  \code{\link{lavTablesFitCf}} function. If \code{"browne.residual.adf"} (alias
  \code{"browne"}) or \code{"browne.residual.nt"}, the standard chi-squared
  difference is calculated from each model's residual-based statistic.}
\item{model.names}{Character vector. If provided, use these model names in
  the first column of the anova table.}
}
\value{
    An object of class anova. When given a single argument, it simply 
    returns the test statistic of this model. When given a sequence of 
    objects, this function tests the models against one another, after
    reordering the models according to their degrees of freedom.
}
\details{
    The \code{anova} function for lavaan objects simply calls the
    \code{lavTestLRT} function, which has a few additional arguments.
    
    The only \code{test=} options that currently have actual consequences are 
    \code{"satorra.bentler"}, \code{"yuan.bentler"}, or \code{"yuan.bentler.mplus"}
    because \code{"mean.var.adjusted"} and \code{"scaled.shifted"} are
    currently distinguished by the \code{scaled.shifted} argument.
    See \code{\link{lavOptions}} for details about \code{test=} options
    implied by robust \code{estimator=} options). The \code{"default"} is to
    select the first available scaled statistic, if any. To check which test(s)
    were calculated when fitting your model(s), use 
    \code{lavInspect(fit, "options")$test}.

    If \code{type = "Chisq"} and the test statistics are scaled, a
    special scaled difference test statistic is computed. If \code{method} is
    \code{"satorra.bentler.2001"}, a simple approximation is used
    described in Satorra & Bentler (2001). In some settings,
    this can lead to a negative test statistic. To ensure a positive
    test statistic, we can use the method proposed by 
    Satorra & Bentler (2010). Alternatively, when \code{method="satorra.2000"},
    the original formulas of Satorra (2000) are used.  The latter is used for
    model comparison, when \code{...} contain additional (nested) models.
    Even when test statistics are scaled in \code{object} or \code{...},
    users may request the  \code{method="standard"} test statistic,
    without a robust adjustment.
}
\note{
  If there is a \code{\linkS4class{lavaan}} model stored in 
  \code{object@external$h1.model}, it will be added to \code{\dots}
}
\references{
Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample
analysis of moment structures. In Heijmans, R.D.H., Pollock, D.S.G. & Satorra,
A. (eds.), \emph{Innovations in multivariate statistical analysis}:
\emph{A Festschrift for Heinz Neudecker} (pp.233-247). 
London, UK: Kluwer Academic Publishers.     

Satorra, A., & Bentler, P. M. (2001). A scaled difference chi-square test
statistic for moment structure analysis. \emph{Psychometrika, 66}(4), 507-514.
\doi{10.1007/BF02296192}

Satorra, A., & Bentler, P. M. (2010). Ensuring postiveness of the scaled
difference chi-square test statistic. \emph{Psychometrika, 75}(2), 243-248.
\doi{10.1007/s11336-009-9135-y}
}
\examples{
HS.model <- '
    visual  =~ x1 + b1*x2 + x3
    textual =~ x4 + b2*x5 + x6
    speed   =~ x7 + b3*x8 + x9
'
fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
fit0 <- cfa(HS.model, data = HolzingerSwineford1939, 
            orthogonal = TRUE)
lavTestLRT(fit1, fit0)


## When multiple test statistics are selected when the model is fitted,
## use the type= and test= arguments to select a test for comparison.

## refit models, requesting 6 test statistics (in addition to "standard")
t6.1 <- cfa(HS.model, data = HolzingerSwineford1939,
            test = c("browne.residual.adf","scaled.shifted","mean.var.adjusted",
                     "satorra.bentler", "yuan.bentler", "yuan.bentler.mplus"))
t6.0 <- cfa(HS.model, data = HolzingerSwineford1939, orthogonal = TRUE,
            test = c("browne.residual.adf","scaled.shifted","mean.var.adjusted",
                     "satorra.bentler", "yuan.bentler", "yuan.bentler.mplus"))

## By default (test="default", type="Chisq"), the first scaled statistic
## requested will be used. Here, that is "scaled.shifted"
lavTestLRT(t6.1, t6.0)
## But even if "satorra.bentler" were requested first, method="satorra.2000"
## provides the scaled-shifted chi-squared difference test:
lavTestLRT(t6.1, t6.0, method = "satorra.2000")
## == lavTestLRT(update(t6.1, test = "scaled.shifted"), update(t6.0, test = "scaled.shifted"))

## The mean- and variance-adjusted (Satterthwaite) statistic implies
## scaled.shifted = FALSE
lavTestLRT(t6.1, t6.0, method = "satorra.2000", scaled.shifted = FALSE)

## Because "satorra.bentler" is not the first scaled test in the list,
## we MUST request it explicitly:
lavTestLRT(t6.1, t6.0, test = "satorra.bentler") # method="satorra.bentler.2001"
## == lavTestLRT(update(t6.1, test = "satorra.bentler"),
##               update(t6.0, test = "satorra.bentler"))
## The "strictly-positive test" is necessary when the above test is < 0:
lavTestLRT(t6.1, t6.0, test = "satorra.bentler", method = "satorra.bentler.2010")

## Likewise, other scaled statistics can be selected:
lavTestLRT(t6.1, t6.0, test = "yuan.bentler")
## == lavTestLRT(update(t6.1, test = "yuan.bentler"),
##               update(t6.0, test = "yuan.bentler"))
lavTestLRT(t6.1, t6.0, test = "yuan.bentler.mplus")
## == lavTestLRT(update(t6.1, test = "yuan.bentler.mplus"),
##               update(t6.0, test = "yuan.bentler.mplus"))

## To request the difference between Browne's (1984) residual-based statistics,
## rather than statistics based on the fitted model's discrepancy function,
## use the type= argument:
lavTestLRT(t6.1, t6.0, type = "browne.residual.adf")

## Despite requesting multiple robust tests, it is still possible to obtain
## the standard chi-squared difference test (i.e., without a robust correction)
lavTestLRT(t6.1, t6.0, method = "standard")
## == lavTestLRT(update(t6.1, test = "standard"), update(t6.0, test = "standard"))

}
