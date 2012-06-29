# simple wrapper around lm.fit to get scores()
# 
# YR 25 June 2012
#
# NOTES: - X should NOT already contain a column of 1's for the intercept!
#        - weights and offset not used yet!!

# wrapper function
lavOLS <- function(y, X = NULL, 
                   weights = rep(1, length(y)),
                   offset  = rep(0, length(y)),
                   method = "none", start = NULL,
                   control = list(), verbose = FALSE) {

    # initialize
    lavR <- lavRefOLS$new(y = y, X = X, weights = weights, offset = offset)

    # optimize
    lavR$optimize(method = method, control = control, verbose = verbose,
                  start = start)

    lavR
}

# lavRefOLS
#
# classic univariate OLS regression
lavRefOLS <- setRefClass("lavOLS",

# inherits
contains = "lavOptim",

# fields
fields = list(y = "numeric", X = "matrix", 
              nobs = "integer", nexo = "integer", 
              weights = "numeric", offset = "numeric",
              beta.idx = "integer", var.idx = "integer",
              # internal
              yhat = "numeric"),

# methods
methods = list(

initialize = function(y, X = NULL,
                      weights = rep(1, length(y)),
                      offset = rep(0, length(y))) {
    # y
    y <<- as.numeric(y); nobs <<- length(y)
    # X
    if(!is.null(X)) {
        nexo <<- ncol(X)
        X <<- cbind(1,X) # add intercept
    } else {
        nexo <<- 0L
    }
    # weights and offset
    weights <<- weights; offset <<- offset

    # indices of free parameters
    beta.idx <<- 1:(nexo + 1L) # note INCLUDES intercept (unlike lavProbit)
     var.idx <<- 1L + nexo + 1L
    
    # set up for Optim
    npar  <<- 1L + nexo + 1L # intercept + beta + var
    theta <<- rep(as.numeric(NA), npar)
},

objective = function(x) {
    if(!missing(x)) theta <<- x
    beta  <- theta[-npar] # not the variance
    e.var <- theta[npar]  # the variance of the error
    if(nexo > 0L) 
        yhat <<- drop(X %*% beta) + offset
    else
        yhat <<- rep(beta[1], nobs)
    log.dy <- dnorm(y, mean=yhat, sd=sqrt(e.var), log=TRUE)
    -sum(weights * log.dy)
},

gradient = function(x) {
    if(!missing(x)) objective(x)
    e.var <- theta[npar]

    # beta
    if(nexo > 0L) {
        dx.beta <- 1/e.var * crossprod(X, y - yhat)
    } else {
        dx.beta <- 1/e.var * sum(y - yhat)
    }
    # var
    dx.var  <- -nobs/(2*e.var) + 1/(2*e.var^2) * crossprod(y - yhat)

    -c(dx.beta, dx.var)
},

scores = function(x) {
    if(!missing(x)) objective(x)
    e.var <- theta[npar]

    # beta
    scores.beta <- 1/e.var * X * (y - yhat)
    # var
    scores.var  <- -1/(2*e.var) + 1/(2*e.var^2) * (y - yhat)^2

    cbind(-scores.beta, -scores.var)
},

hessian = function(x) {
    if(!missing(x)) { objective(x); gradient() }
    #cat("hessian num = \n"); print(round(numDeriv:::hessian(func=.self$objective, x=x),3))
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

    H <- rbind( cbind(dx2.beta, dx.beta.var),
                cbind(t(dx.beta.var), dx2.var) )
    -H
},

start = function() {
    if(nexo > 0L) {
        fit.lm <- lm.fit(y=y, x=X)
        beta.start <- fit.lm$coef
         var.start <- crossprod(fit.lm$residual)/nobs
    } else {
        beta.start <- mean(y, na.rm=TRUE)
         var.start <-  var(y, na.rm=TRUE)*(nobs-1)/nobs # ML
    }
    c( beta.start, var.start )
}

))

