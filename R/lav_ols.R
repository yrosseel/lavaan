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
              missing.values = "logical", missing.idx = "integer",
              # internal
              yhat = "numeric"),

# methods
methods = list(

initialize = function(y, X = NULL, ...) {
    # y
    .self$y <- y; .self$nobs <- length(y)

    # X
    if(!is.null(X)) {
        .self$nexo <- ncol(X)
        .self$X <- matrix(c(rep.int(1,nobs),X),nobs,nexo+1L)
    } else {
        .self$nexo <- 0L
        .self$X <- matrix(1, nobs, 1L)
    }
    if(any(is.na(y)) || (!is.null(X) && any(is.na(X)) )) {
        .self$missing.values <- TRUE
        .self$missing.idx <- which(apply(cbind(y, X), 1, function(x) any(is.na(x))))
    } else {
        .self$missing.values <- FALSE
    }

    # indices of free parameters
      .self$int.idx <- 1L
    .self$slope.idx <- seq_len(nexo) + 1L
      .self$var.idx <- 1L + nexo + 1L

    # set up for Optim
    .self$npar  <- 1L + nexo + 1L # intercept + slopes + var
    start(); .self$theta <- theta.start
    if(nexo > 0L)
        sl.lab <- paste("beta",seq_len(nexo),sep="")
    else
        sl.lab <- character(0)
    .self$theta.labels <- c("int", sl.lab, "var.e")
},

start = function() {
    if(nexo > 0L) {
        if(length(missing.idx) > 0L) {
            fit.lm <- lm.fit(y=y[-missing.idx], x=X[-missing.idx,,drop=FALSE])
        } else {
            fit.lm <- lm.fit(y=y, x=X)
        }
        #fit.lm <- lm.wfit(y=y, x=X, w=weights)
        beta.start <- fit.lm$coef
         var.start <- crossprod(fit.lm$residual)/nobs
    } else {
        beta.start <- mean(y, na.rm=TRUE)
         var.start <-  var(y, na.rm=TRUE)*(nobs-1)/nobs # ML
    }
    .self$theta.start <- c( beta.start, var.start )
},

lik = function(x) {
    if(!missing(x)) .self$theta <- x
    beta  <- theta[-npar] # not the variance
    e.var <- theta[npar]  # the variance of the error
    if(nexo > 0L)
        .self$yhat <- drop(X %*% beta)
    else
        .self$yhat <- rep(beta[1L], nobs)
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
    scores.var  <- -1/(2*e.var) + 1/(2*e.var*e.var) * (y - yhat)*(y - yhat)

    cbind(scores.beta, scores.var, deparse.level=0)
},

hessian = function(x) {
    if(!missing(x)) { lik(x); gradient() }
    #cat("hessian num = \n"); print(round(numDeriv::hessian(func=.self$logl, x=x),3))
    e.var <- theta[npar]

    # beta - beta
    if(nexo > 0L) {
        dx2.beta <- -1/e.var * crossprod(X)
    } else {
        dx2.beta <- -1/e.var * nobs
    }
    # beta - var
    if(nexo > 0L) {
        dx.beta.var <- -1/(e.var*e.var) * crossprod(X, y-yhat)
    } else {
        dx.beta.var <- -1/(e.var*e.var) * sum(y-yhat)
    }
    # var - var
    sq.e.var <- sqrt(e.var)
    sq.e.var6 <- sq.e.var*sq.e.var*sq.e.var*sq.e.var*sq.e.var*sq.e.var
    #dx2.var <- nobs/(2*e.var*e.var) - 1/sq.e.var6 * crossprod(y-yhat)
    dx2.var <- nobs/(2*e.var*e.var) - 1/sq.e.var6 * (e.var * nobs)

    rbind( cbind(  dx2.beta, dx.beta.var, deparse.level=0),
           cbind(t(dx.beta.var), dx2.var, deparse.level=0), deparse.level=0 )
},

minObjective = function(x) { -logl(x)     },
minGradient  = function(x) { -gradient(x) },
minHessian   = function(x) { -hessian(x)  }

))

