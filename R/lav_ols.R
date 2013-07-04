# simple wrapper around lm.fit to get scores()
# 
# YR 25 June 2012
#
# NOTES: - X should NOT already contain a column of 1's for the intercept!
#        - weights not used yet

# wrapper function
lavOLS <- function(y, X = NULL, 
                   method = "none", start.values = NULL,
                   control = list(), verbose = FALSE) {

    # initialize
    lavR <- lavRefOLS$new(y = y, X = X)

    # optimize
    if(method != "none") {
        lavR$optimize(method = method, control = control, verbose = verbose,
                      start.values = start.values)
    } else {
        lavR$lik() # initialize
    }

    lavR
}

# lavRefOLS
#
# classic univariate OLS regression
lavRefOLS <- setRefClass("lavOLS",

# inherits
contains = "lavML",

# fields
fields = list(X = "matrix", nexo = "integer", 
              # housekeeping
              int.idx = "integer", slope.idx = "integer", var.idx = "integer",
              # internal
              yhat = "numeric"),

# methods
methods = list(

initialize = function(y, X = NULL, ...) {
    # y
    y <<- y; nobs <<- length(y)

    # X
    if(!is.null(X)) {
        nexo <<- ncol(X)
        X <<- matrix(c(rep.int(1,nobs),X),nobs,nexo+1L)
    } else {
        nexo <<- 0L
        X <<- matrix(1, nobs, 1L)
    }

    # indices of free parameters
      int.idx <<- 1L
    slope.idx <<- seq_len(nexo) + 1L
      var.idx <<- 1L + nexo + 1L
    
    # set up for Optim
    npar  <<- 1L + nexo + 1L # intercept + slopes + var
    start(); theta <<- theta.start
    if(nexo > 0L)
        sl.lab <- paste("beta",seq_len(nexo),sep="")
    else
        sl.lab <- character(0)
    theta.labels <<- c("int", sl.lab, "var.e")
},

start = function() {
    if(nexo > 0L) {
        fit.lm <- lm.fit(y=y, x=X)
        #fit.lm <- lm.wfit(y=y, x=X, w=weights)
        beta.start <- fit.lm$coef
         var.start <- crossprod(fit.lm$residual)/nobs
    } else {
        beta.start <- mean(y, na.rm=TRUE)
         var.start <-  var(y, na.rm=TRUE)*(nobs-1)/nobs # ML
    }
    theta.start <<- c( beta.start, var.start )
},

lik = function(x) {
    if(!missing(x)) theta <<- x
    beta  <- theta[-npar] # not the variance
    e.var <- theta[npar]  # the variance of the error
    if(nexo > 0L) 
        yhat <<- drop(X %*% beta)
    else
        yhat <<- rep(beta[1L], nobs)
    #weights * dnorm(y, mean=yhat, sd=sqrt(e.var))
    dnorm(y, mean=yhat, sd=sqrt(e.var))
},

#gradient = function(x) {
#    if(!missing(x)) logl(x)
#    e.var <- theta[npar]
#
#    # beta
#    if(nexo > 0L) {
#        dx.beta <- 1/e.var * crossprod(X, y - yhat)
#    } else {
#        dx.beta <- 1/e.var * sum(y - yhat)
#    }
#    # var
#    dx.var  <- -nobs/(2*e.var) + 1/(2*e.var^2) * crossprod(y - yhat)
#
#    c(dx.beta, dx.var)
#},

scores = function(x) {
    if(!missing(x)) lik(x)
    e.var <- theta[npar]
    if(length(yhat) == 0L) lik() # not initialized yet

    # beta
    scores.beta <- 1/e.var * X * (y - yhat)
    # var
    scores.var  <- -1/(2*e.var) + 1/(2*e.var^2) * (y - yhat)^2

    cbind(scores.beta, scores.var, deparse.level=0)
},

hessian = function(x) {
    if(!missing(x)) { lik(x); gradient() }
    #cat("hessian num = \n"); print(round(numDeriv:::hessian(func=.self$logl, x=x),3))
    e.var <- theta[npar]

    # beta - beta
    if(nexo > 0L) {
        dx2.beta <- -1/e.var * crossprod(X)
    } else {
        dx2.beta <- -1/e.var * nobs
    }
    # beta - var
    if(nexo > 0L) {
        dx.beta.var <- -1/e.var^2*crossprod(X, y-yhat)
    } else {
        dx.beta.var <- -1/e.var^2* sum(y-yhat)
    }
    # var - var           
    dx2.var <- nobs/(2*e.var^2) - 1/sqrt(e.var)^6*crossprod(y-yhat)

    rbind( cbind(  dx2.beta, dx.beta.var, deparse.level=0),
           cbind(t(dx.beta.var), dx2.var, deparse.level=0), deparse.level=0 )
},

minObjective = function(x) { -logl(x)     },
minGradient  = function(x) { -gradient(x) },
minHessian   = function(x) { -hessian(x)  }

))

