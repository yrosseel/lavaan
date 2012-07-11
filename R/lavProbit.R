# ordered probit regression
# 
# YR 21 June 2012
#
# why not using MASS:::polr?
# - it does not deal with binary responses (must use glm.fit instead)
# - we need scores
# - Newton-Raphson is much faster
# - allow for empty X, just to get thresholds (and scores)
# - however, we do NOT force thresholds to be strictly positive!

# NOTE: X should NOT contain a column of 1's for the intercept!!


# wrapper function
lavProbit <- function(y, X=NULL, weights = rep(1, length(y)),
                      offset = rep(0, length(y)), fast=FALSE, 
                      method = "nlminb.hessian", control = list(),
                      verbose = FALSE) {

    # initialize ref class
    lavR <- lavRefProbit$new(y = y, X = X, weights = weights, offset = offset)

    # optimize (only if X)
    if(!is.null(X)) {
        lavR$optimize(method = method, control = control, verbose = verbose)
    } else {
        lavR$theta <- lavR$start()
        lavR$lik() # initialize all elements (z1,z2,Y1,Y2)
    }

    lavR
}

# lavRefProbit
#
# classic probit regression
lavRefProbit <- setRefClass("lavProbit",

# inherits
contains = "lavML",

# fields
fields = list(y = "integer", X = "matrix", 
              nobs = "integer", nexo = "integer", nth = "integer", 
              weights = "numeric", offset = "numeric",
              Y1 = "matrix", Y2 = "matrix",
              th.idx = "integer", slope.idx = "integer",
              # cache:
              # this doesn't result in better speed; but makes the code cleaner
              z1 = "numeric", z2 = "numeric", probits = "numeric",
              p1 = "numeric", p2 = "numeric"),

# methods
methods = list(

initialize = function(y, X=NULL,
                      weights = rep(1, length(y)),
                      offset = rep(0, length(y))) {
    # y
    y <<- as.integer(y); nth <<- length(unique(.self$y)) - 1L
    nobs <<- length(y)
    # X
    if(is.null(X)) { 
        nexo <<- 0L
    } else {
        X <<- unname(X); nexo <<- ncol(X)
    }

    # weights and offset
    weights <<- weights; offset <<- offset
    
    # TH matrices (TRUE/FALSE)
    Y1 <<- matrix(1:nth, nobs, nth, byrow=TRUE) == .self$y
    Y2 <<- matrix(1:nth, nobs, nth, byrow=TRUE) == (.self$y - 1L)

    # indices of free parameters
    th.idx <<- 1:nth
    slope.idx <<- seq_len(nexo) + nth

    # set up for Optim
    npar  <<- nth + nexo
    start(); theta <<- theta.start
    th.lab <- paste("th",  seq_len(nth), sep="")
    sl.lab <- character(0)
    if(nexo > 0L)
        sl.lab <- paste("beta",seq_len(nexo),sep="")
    theta.labels <<- c(th.lab, sl.lab)
},

start = function() {
    th.start <- lavaan:::pc_th(freq=tabulate(y)) # unconditional th's
    beta.start <- rep(0, nexo)
    #Y <- as.numeric(y); range <- 16; Y <- Y*range/nexo
    #fit.ols <- lavOLS(y=Y, X=X, weights=weights, offset=offset)
    #beta.start <- fit.ols$theta[fit.ols$slope.idx[-1L]]
    #print(beta.start)
    theta.start <<- c( th.start, beta.start )
},

# formulas: see Maddala 1983 pages 48 + 49
lik = function(x) {
    if(!missing(x)) theta <<- x
    th <- theta[1:nth]; TH <- c(-Inf, th, +Inf); beta <- theta[-c(1:nth)]
    if(nexo > 0L) 
        eta <- drop(X %*% beta) + offset
    else
        eta <- numeric(nobs)
    z1 <<- pmin( 100, TH[y+1L   ] - eta)
    z2 <<- pmax(-100, TH[y+1L-1L] - eta)
    probits <<- pnorm(z1) - pnorm(z2)
    probits
},

scores = function(x) {
    if(!missing(x)) lik(x)
    if(length(probits) == 0L) lik()
    p1 <<- dnorm(z1); p2 <<- dnorm(z2)

    # th
    scores.th   <- -1 * (Y2*p2 - Y1*p1) * (weights/probits)

    # beta
    scores.beta <- matrix(0, length(p1), 0L)
    if(nexo > 0L)
        scores.beta <- -1 * weights*(p1 - p2)/probits * X

    cbind(scores.th, scores.beta)
},

#gradient = function(x) {
#    if(!missing(x)) objective(x)
#    p1 <<- dnorm(z1); p2 <<- dnorm(z2)
#    
#    # beta
#    dx.beta <- numeric(0L)
#    if(nexo > 0L)
#        dx.beta <- crossprod(X, weights*(p1 - p2)/probits)
#    # th
#    dx.th   <- crossprod(Y2*p2 - Y1*p1, weights/probits)
#
#    c(dx.th, dx.beta)
#},

hessian = function(x) {
    if(!missing(x)) { lik(x); gradient() }
    #cat("hessian num = \n"); print(round(numDeriv:::hessian(func=.self$objective, x=x),3))
    if(length(probits) == 0L) lik(); scores() # not initialized
    gnorm <- function(x) { -x * dnorm(x) }
    wtpr <- weights/probits

    dxa <- Y1*p1 - Y2*p2
    dx2.alpha <- -1 * (crossprod(dxa, (dxa * wtpr / probits)) -
                        ( crossprod(Y1 * gnorm(z1) * wtpr, Y1) -
                          crossprod(Y2 * gnorm(z2) * wtpr, Y2) ) )

    # only for empty X
    if(nexo == 0L) return(dx2.alpha)

    dxb <-  X*p1 - X*p2
    dx2.beta <- -1 * (crossprod(dxb, (dxb * wtpr / probits)) -
                       ( crossprod(X * gnorm(z1) * wtpr, X) - 
                         crossprod(X * gnorm(z2) * wtpr, X) ) )

    dx.ab <- crossprod(dxa, (dxb * wtpr / probits)) -
               ( crossprod(Y1 * gnorm(z1) * wtpr, X) -
                 crossprod(Y2 * gnorm(z2) * wtpr, X) )

    rbind( cbind(dx2.alpha,   dx.ab, deparse.level=0),
           cbind(t(dx.ab), dx2.beta, deparse.level=0)  )
},

minObjective = function(x) { -logl(x)     },
minGradient  = function(x) { -gradient(x) },
minHessian   = function(x) { -hessian(x)  }

))

