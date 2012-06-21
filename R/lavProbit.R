# various flavors of ordered probit regression
# 
# YR 21 June 2012

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
fields = list( nth = "integer", nx = "integer", freq    = "integer",
    y = "integer", Y1 = "matrix", Y2 = "matrix", X = "matrix", 
    probits = "numeric", eta = "numeric",
    wts = "numeric", wtpr = "numeric", offset = "numeric",
    dpi.psi = "matrix"),

# methods
methods = list(

initialize = function(y, X, weights = rep(1, length(y)),
                      offset = rep(0, length(y))) {
    cat("this is lavProbit -- classic!\n")
    # y
    y <<- as.integer(y); nth <<- nlevels(y) - 1L; freq <<- tabulate(y)
    # X
    X <<- X; n <- nrow(X); nx <<- ncol(X); nbeta <- NCOL(X)

    Y <- matrix(0, n, nth)
    Y1 <<- col(Y) == .self$y; Y2 <<- col(Y) == (.self$y - 1L)

    probits <<- numeric(length = n)
    wts <<- weights; offset <<- offset

    # set up for Optim
    npar <<- nth + nbeta
    theta <<- rep(as.numeric(NA), npar)
},

# formulas: see Maddala 1983 pages 48 + 49
objective = function(x) {
    if(!missing(x)) theta <<- x
    th <- theta[1:nth]; TH <- c(-Inf, th, +Inf); beta <- theta[-c(1:nth)]
    eta <<- drop(X %*% beta) + offset
    z1 <- pmin(100,TH[y+1L]  - eta); z2 <- pmax(-100, TH[y+1L-1L] - eta)
    probits <<- pnorm(z1) - pnorm(z2)
    if(all(probits > 0))
        -sum(wts * log(probits))
    else Inf
},

gradient = function(x) {
    if(!missing(x)) objective(x)
    th <- theta[1:nth]; TH <- c(-Inf, th, +Inf); beta <- theta[-c(1:nth)]
    eta <<- drop(X %*% beta) + offset
    z1 <- pmin(100,TH[y+1L]  - eta); z2 <- pmax(-100, TH[y+1L-1L] - eta)
    probits <<- pnorm(z1) - pnorm(z2)
    p1 <- dnorm(z1); p2 <- dnorm(z2)
    
    # beta
    dx.beta <- (wts * (p1 - p2)/probits) %*% X
    # th
    dx.th <- -crossprod(Y1 * p1 - Y2 * p2, wts/probits)

    c(dx.th, dx.beta)
},

hessian = function(x) {
    if(!missing(x)) { objective(x); gradient() }
    #cat("hessian num = \n"); print(round(numDeriv:::hessian(func=.self$objective, x=x),3))
    th <- theta[1:nth]; TH <- c(-Inf, th, +Inf); beta <- theta[-c(1:nth)]
    z1 <- pmin(100,TH[y+1L]  - eta); z2 <- pmax(-100, TH[y+1L-1L] - eta)
    p1 <- dnorm(z1); p2 <- dnorm(z2)
    gnorm <- function(x) { -x * dnorm(x) }

    wtpr <<- wts/probits
    dxb <-  X*p1 - X*p2
    dx2.beta <- -(crossprod(X * gnorm(z1) * wtpr, X) - 
                   crossprod(X * gnorm(z2) * wtpr, X)) +
                 crossprod(dxb, (dxb * wtpr / probits))

    dxa <- Y1*p1 - Y2*p2
    dx2.alpha <- -(crossprod(Y1 * gnorm(z1) * wtpr, Y1) -
                   crossprod(Y2 * gnorm(z2) * wtpr, Y2)) +
                  crossprod(dxa, (dxa * wtpr / probits))

    dx.ab <- -(crossprod(Y1 * gnorm(z1) * wtpr, X) -
               crossprod(Y2 * gnorm(z2) * wtpr, X)) +
              crossprod(dxa, (dxb * wtpr / probits))
    dx.ab <- -1*dx.ab ##??? FIXME?

    H <- rbind( cbind(dx2.alpha, dx.ab),
                cbind(t(dx.ab), dx2.beta) )

    #print(round(H,3))
},

start = function() {
    #th.start <- sort(rnorm(nth))
    th.start <- lavaan:::unithord(freq=freq) # unconditional th's
    beta.start <- rep(0, nx)
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
# in practice, this is only about 15-20% faster
#
# this function will most likely change a lot in the future
# it is one of my first attempts to exploit the setRefClass functionality
# in R
lavRefProbitB <- setRefClass("lavProbitB",

# inherits
contains = "lavOptim",

# fields
fields = list( nth = "integer", nx = "integer", freq    = "integer",
    B1 = "matrix", B2 = "matrix", o1 = "numeric", o2 = "numeric",
    probits = "numeric", eta1 = "numeric", eta2 = "numeric",
    wts = "numeric", wtpr = "numeric", dpi.psi = "matrix"),

# methods
methods = list(

initialize = function(y, X, weights = rep(1, length(y)), 
                      offset = rep(0, length(y))) {
    # set up internal fields
    nth <<- nlevels(y) - 1L; freq <<- tabulate(y)
    n <- nrow(X); nx <<- ncol(X)
    B2.tmp <- 1 * (col(matrix(0, n, nth + 1L)) == c(unclass(y)))
    o1 <<-  c( 1e5 * B2.tmp[, nth + 1L]) - offset
    o2 <<-  c(-1e5 * B2.tmp[,1L]) - offset
    B1 <<-  B2.tmp[, -(nth + 1L), drop = FALSE]
    B2 <<-  B2.tmp[, -1L, drop = FALSE]
    nbeta <- NCOL(X) 
    if(nbeta > 0) {
        B1 <<- cbind(B1, -X)
        B2 <<- cbind(B2, -X)
    }
    dimnames(B1) <<- NULL; dimnames(B2) <<- NULL
    probits <<- numeric(length = n)
    wts <<- weights
    npar <<- nth + nbeta
    theta <<- rep(as.numeric(NA), npar)
},

objective = function(x) {
    if(!missing(x)) theta <<- x
    eta1 <<- drop(B1 %*% theta) + o1; eta2 <<- drop(B2 %*% theta) + o2
    probits <<- pnorm(eta1) - pnorm(eta2)
    if(all(probits > 0))
          -sum(wts * log(probits))
    else Inf    
},

gradient = function(x) {
    if(!missing(x)) objective(x)
    p1 <- dnorm(eta1); p2 <- dnorm(eta2)
    wtpr <<- wts/probits
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
    beta.start <- rep(0, nx)
    c( th.start, beta.start )
}

))


