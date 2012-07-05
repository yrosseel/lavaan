# ordered probit regression
# 
# YR 21 June 2012
#
# why not using MASS:::polr?
# - it does not deal with binary responses (must use glm.fit instead)
# - we need scores
# - Newton-Raphson is must faster
# - allow for empty X
# - however, we do force thresholds to be strictly positive!

# NOTE: X should NOT contain a column of 1's for the intercept!!


# wrapper function
lavProbit <- function(y, X, weights = rep(1, length(y)),
                      offset = rep(0, length(y)), fast=FALSE, 
                      method = "nlminb.hessian", control = list(),
                      verbose = FALSE) {

    # initialize ref class
    if(fast) {
        lavR <- lavRefProbitB$new(y = y, X = X, weights = weights, 
                                  offset = offset)
    } else {
        lavR <- lavRefProbit$new(y = y, X = X, weights = weights, 
                                 offset = offset)
    }

    # optimize
    lavR$optimize(method = method, control = control, verbose = verbose)

    lavR
}

# lavRefProbit
#
# classic probit regression
lavRefProbit <- setRefClass("lavProbit",

# inherits
contains = "lavOptim",

# fields
fields = list(y = "integer", X = "matrix", 
              nobs = "integer", nexo = "integer", nth = "integer", 
              weights = "numeric", offset = "numeric",
              Y1 = "matrix", Y2 = "matrix",
              th.idx = "integer", beta.idx = "integer",
              # cache:
              # this doesn't result in better speed; but makes the code cleaner
              z1 = "numeric", z2 = "numeric", probits = "numeric",
              p1 = "numeric", p2 = "numeric"),

# methods
methods = list(

initialize = function(y, X,
                      weights = rep(1, length(y)),
                      offset = rep(0, length(y))) {
    # y
    y <<- as.integer(y); nth <<- length(unique(.self$y)) - 1L
    nobs <<- length(y)
    # X
    if(!missing(X)) {
        X <<- X; nexo <<- ncol(X)
    } else {
        nexo <<- 0L
    }
    # weights and offset
    weights <<- weights; offset <<- offset
    
    # TH matrices (TRUE/FALSE)
    Y1 <<- matrix(1:nth, nobs, nth, byrow=TRUE) == .self$y
    Y2 <<- matrix(1:nth, nobs, nth, byrow=TRUE) == (.self$y - 1L)

    # indices of free parameters
    th.idx <<- 1:nth
    beta.idx <<- nth + (1:nexo)

    # set up for Optim
    npar  <<- nth + nexo
    theta <<- rep(as.numeric(NA), npar)
},

# formulas: see Maddala 1983 pages 48 + 49
objective = function(x) {
    if(!missing(x)) theta <<- x
    th <- theta[1:nth]; TH <- c(-Inf, th, +Inf); beta <- theta[-c(1:nth)]
    if(nexo > 0L) 
        eta <- drop(X %*% beta) + offset
    else
        eta <- numeric(nobs)
    z1 <<- pmin( 100, TH[y+1L   ] - eta)
    z2 <<- pmax(-100, TH[y+1L-1L] - eta)
    probits <<- pnorm(z1) - pnorm(z2)
    if(all(probits > 0))
        -sum(weights * log(probits))
    else Inf
},

gradient = function(x) {
    if(!missing(x)) objective(x)
    p1 <<- dnorm(z1); p2 <<- dnorm(z2)
    
    # beta
    dx.beta <- numeric(0L)
    if(nexo > 0L)
        dx.beta <- crossprod(X, weights*(p1 - p2)/probits)
    # th
    dx.th   <- crossprod(Y2*p2 - Y1*p1, weights/probits)

    c(dx.th, dx.beta)
},

scores = function(x) {
    if(!missing(x)) objective(x)
    p1 <<- dnorm(z1); p2 <<- dnorm(z2)

    # beta
    scores.beta <- matrix(0, length(p1), 0L)
    if(nexo > 0L) 
        scores.beta <- weights*(p1 - p2)/probits * X
    # th
    scores.th   <- (Y2*p2 - Y1*p1) * (weights/probits)

    cbind(scores.th, -scores.beta)
},

hessian = function(x) {
    if(!missing(x)) { objective(x); gradient() }
    #cat("hessian num = \n"); print(round(numDeriv:::hessian(func=.self$objective, x=x),3))
    gnorm <- function(x) { -x * dnorm(x) }
    wtpr <- weights/probits

    dxa <- Y1*p1 - Y2*p2
    dx2.alpha <- crossprod(dxa, (dxa * wtpr / probits)) -
                   ( crossprod(Y1 * gnorm(z1) * wtpr, Y1) -
                     crossprod(Y2 * gnorm(z2) * wtpr, Y2) )

    # only for empty X
    if(nexo == 0L) return(dx2.alpha)

    dxb <-  X*p1 - X*p2
    dx2.beta <- crossprod(dxb, (dxb * wtpr / probits)) -
                  ( crossprod(X * gnorm(z1) * wtpr, X) - 
                    crossprod(X * gnorm(z2) * wtpr, X) )

    dx.ab <- crossprod(dxa, (dxb * wtpr / probits)) -
               ( crossprod(Y1 * gnorm(z1) * wtpr, X) -
                 crossprod(Y2 * gnorm(z2) * wtpr, X) )
    dx.ab <- -1*dx.ab 

    H <- rbind( cbind(dx2.alpha,   dx.ab),
                cbind(t(dx.ab), dx2.beta)  )
},

start = function() {
    th.start <- lavaan:::unithord(freq=tabulate(y)) # unconditional th's
    beta.start <- rep(0, nexo)
    #Y <- as.numeric(y); range <- 16; Y <- Y*range/nexo
    #fit.ols <- lavOLS(y=Y, X=X, weights=weights, offset=offset)
    #beta.start <- fit.ols$theta[fit.ols$beta.idx[-1L]]
    print(beta.start)
    c( th.start, beta.start )
}

))



# lavRefProbitB
#
# this version is using the 'B1/B2' trick to save some computations
# this trick is due to Rune Haubo B Christensen <rhbc@imm.dtu.dk>
# see his presentation (slides 35-41) about 'clm' at 
# http://eeecon.uibk.ac.at/psychoco/2012/slides/Christensen.pdf
# see also the 'ordinal' package
#
# in practice, this seems to be about 10-15% faster
#
lavRefProbitB <- setRefClass("lavProbitB",

# inherits
contains = "lavOptim",

# fields
fields = list( nth = "integer", nexo = "integer", freq    = "integer",
    B1 = "matrix", B2 = "matrix", o1 = "numeric", o2 = "numeric",
    probits = "numeric", eta1 = "numeric", eta2 = "numeric",
    weights = "numeric", wtpr = "numeric", dpi.psi = "matrix"),

# methods
methods = list(

initialize = function(y, X, weights = rep(1, length(y)), 
                      offset = rep(0, length(y))) {
    # set up internal fields
    nth <<- nlevels(y) - 1L; freq <<- tabulate(y)
    n <- nrow(X); nexo <<- ncol(X)
    B2.tmp <- 1 * (col(matrix(0, n, nth + 1L)) == c(unclass(y)))
    o1 <<-  c( 1e5 * B2.tmp[, nth + 1L]) - offset
    o2 <<-  c(-1e5 * B2.tmp[,1L]) - offset
    B1 <<-  B2.tmp[, -(nth + 1L), drop = FALSE]
    B2 <<-  B2.tmp[, -1L, drop = FALSE]
    nexo <<- NCOL(X) 
    if(nexo > 0) {
        B1 <<- cbind(B1, -X)
        B2 <<- cbind(B2, -X)
    }
    dimnames(B1) <<- NULL; dimnames(B2) <<- NULL
    probits <<- numeric(length = n)
    weights <<- weights
    npar <<- nth + nexo
    theta <<- rep(as.numeric(NA), npar)
},

objective = function(x) {
    if(!missing(x)) theta <<- x
    eta1 <<- drop(B1 %*% theta) + o1; eta2 <<- drop(B2 %*% theta) + o2
    probits <<- pnorm(eta1) - pnorm(eta2)
    if(all(probits > 0))
          -sum(weights * log(probits))
    else Inf    
},

gradient = function(x) {
    if(!missing(x)) objective(x)
    p1 <- dnorm(eta1); p2 <- dnorm(eta2)
    wtpr <<- weights/probits
    dpi.psi <<- B1 * p1 - B2 * p2
    -crossprod(dpi.psi, wtpr)    
},

hessian = function(x) {
    if(!missing(x)) { objective(x); gradient(x) }
    gnorm <- function(x) { -x * dnorm(x) }
    dg.psi <- crossprod(B1 * gnorm(eta1) * wtpr, B1) -
              crossprod(B2 * gnorm(eta2) * wtpr, B2)
    -dg.psi + crossprod(dpi.psi, (dpi.psi * wtpr / probits))
},

start = function() {
    th.start <- lavaan:::unithord(freq=freq) # unconditional th's
    beta.start <- rep(0, nexo)
    c( th.start, beta.start )
}

))


