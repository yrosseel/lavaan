# ordered probit regression
#
# YR 21 June 2012
#
# why not using MASS::polr?
# - it does not deal with binary responses (must use glm.fit instead)
# - we need scores
# - Newton-Raphson is much faster
# - allow for empty X, just to get thresholds (and scores)
# - however, we do NOT force thresholds to be strictly positive!
#   (which is mainly a problem if you allow for 'zero' frequencies, I think)

# NOTE: X should NOT contain a column of 1's for the intercept!!


# wrapper function
lavProbit <- function(y, X=NULL, y.levels=length(tabulate(y)),
                      weights = rep(1, length(y)),
                      offset = rep(0, length(y)), fast=FALSE,
                      method = "nlminb.hessian", control = list(),
                      verbose = FALSE) {

    # sanity check
    y.freq <- tabulate(y, nbins=y.levels)
    if(!missing(y.levels) && y.levels < length(y.freq))
        stop("y.levels smaller than number of categories in y")
    if(y.freq[1L] == 0L)
        warning("first category of y has zero observations")
    if(y.freq[y.levels] == 0L)
        warning("last category of y has zero observations")
    y.middle <- y.freq[-c(1L, y.levels)]
    if(length(y.middle) > 0L && any(y.middle == 0L))
        stop("zero counts in middle categories; please recode")

    # initialize ref class
    lavR <- lavRefProbit$new(y = y, X = X, y.levels=y.levels,
                             weights = weights, offset = offset)

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
              missing.values = "logical", missing.idx = "integer",
              Y1 = "matrix", Y2 = "matrix",
              th.idx = "integer", slope.idx = "integer",
              # cache:
              # this doesn't result in better speed; but makes the code cleaner
              z1 = "numeric", z2 = "numeric", probits = "numeric",
              p1 = "numeric", p2 = "numeric"),

# methods
methods = list(

initialize = function(y, X=NULL, y.levels=length(tabulate(y)),
                      weights = rep(1, length(y)),
                      offset = rep(0, length(y))) {
    # y
    .self$y <- as.integer(y); .self$nth <- as.integer(y.levels - 1L)
    .self$nobs <- length(y)
    # X
    if(is.null(X)) {
        .self$nexo <- 0L
    } else {
        .self$X <- unname(X); .self$nexo <- ncol(X)
    }
    if(any(is.na(y)) || (!is.null(X) && any(is.na(X)) )) {
        .self$missing.values <- TRUE
        .self$missing.idx <- which(apply(cbind(y, X), 1, function(x) any(is.na(x))))
    } else {
        .self$missing.values <- FALSE
    }

    # weights and offset
    .self$weights <- weights; .self$offset <- offset

    # TH matrices (TRUE/FALSE)
    .self$Y1 <- matrix(1:nth, nobs, nth, byrow=TRUE) == .self$y
    .self$Y2 <- matrix(1:nth, nobs, nth, byrow=TRUE) == (.self$y - 1L)

    # indices of free parameters
    .self$th.idx <- 1:nth
    .self$slope.idx <- seq_len(nexo) + nth

    # set up for Optim
    .self$npar  <- nth + nexo
    start(); .self$theta <- theta.start
    th.lab <- paste("th",  seq_len(nth), sep="")
    sl.lab <- character(0)
    if(nexo > 0L)
        sl.lab <- paste("beta",seq_len(nexo),sep="")
    .self$theta.labels <- c(th.lab, sl.lab)
},

start = function() {
    if(nth == 1L && nexo > 0L) {
        th.start <- 0
    } else {
        th.start <- pc_th(freq=tabulate(y, nbins=nth+1L)) # unconditional th's
    }
    beta.start <- rep(0, nexo)
    #Y <- as.numeric(y); range <- 16; Y <- Y*range/nexo
    #fit.ols <- lavOLS(y=Y, X=X, weights=weights, offset=offset)
    #beta.start <- fit.ols$theta[fit.ols$slope.idx[-1L]]
    #print(beta.start)
    .self$theta.start <- c( th.start, beta.start )
},

# formulas: see Maddala 1983 pages 48 + 49
lik = function(x) {
    if(!missing(x)) .self$theta <- x
    th <- theta[1:nth]; TH <- c(-Inf, th, +Inf); beta <- theta[-c(1:nth)]
    if(nexo > 0L) {
        eta <- drop(X %*% beta) + offset
    } else {
        eta <- numeric(nobs)
    }
    .self$z1 <- pmin( 100, TH[y+1L   ] - eta)
    .self$z2 <- pmax(-100, TH[y+1L-1L] - eta)
    .self$probits <- pnorm(z1) - pnorm(z2)
    probits
},

scores = function(x) {
    if(!missing(x)) lik(x)
    if(length(probits) == 0L) lik()
    .self$p1 <- dnorm(z1); .self$p2 <- dnorm(z2)

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
#    .self$p1 <- dnorm(z1); .self$p2 <- dnorm(z2)
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
    #cat("hessian num = \n"); print(round(numDeriv::hessian(func=.self$objective, x=x),3))
    if(length(probits) == 0L) lik(); scores() # not initialized
    gnorm <- function(x) { -x * dnorm(x) }
    wtpr <- weights/probits
    dxa <- Y1*p1 - Y2*p2

    # handle missing values -- FIXME!!! better approach?
    # we could also adapt crossprod, to work pairwise...
    if(missing.values) {
        .probits <- probits[-missing.idx]
        .wtpr <- wtpr[-missing.idx]
        .dxa <- dxa[-missing.idx,,drop=FALSE]
        .Y1 <- Y1[-missing.idx,,drop=FALSE]
        .Y2 <- Y2[-missing.idx,,drop=FALSE]
        .z1 <- z1[-missing.idx]
        .z2 <- z2[-missing.idx]
        .X <- X[-missing.idx,,drop=FALSE]
        .p1 <- p1[-missing.idx]
        .p2 <- p2[-missing.idx]
    } else {
        .probits <- probits
        .wtpr <- wtpr
        .dxa <- dxa
        .Y1 <- Y1
        .Y2 <- Y2
        .z1 <- z1
        .z2 <- z2
        .X <- X
        .p1 <- p1
        .p2 <- p2
    }

    dx2.alpha <- -1 * (crossprod(.dxa, (.dxa * .wtpr / .probits)) -
                        ( crossprod(.Y1 * gnorm(.z1) * .wtpr, .Y1) -
                          crossprod(.Y2 * gnorm(.z2) * .wtpr, .Y2) ) )

    # only for empty X
    if(nexo == 0L) return(dx2.alpha)

    dxb <-  .X*.p1 - .X*.p2
    dx2.beta <- -1 * (crossprod(dxb, (dxb * .wtpr / .probits)) -
                       ( crossprod(.X * gnorm(.z1) * .wtpr, .X) -
                         crossprod(.X * gnorm(.z2) * .wtpr, .X) ) )

    dx.ab <- crossprod(.dxa, (dxb * .wtpr / .probits)) -
               ( crossprod(.Y1 * gnorm(.z1) * .wtpr, .X) -
                 crossprod(.Y2 * gnorm(.z2) * .wtpr, .X) )

    rbind( cbind(dx2.alpha,   dx.ab, deparse.level=0),
           cbind(t(dx.ab), dx2.beta, deparse.level=0)  )
},

minObjective = function(x) { -logl(x)     },
minGradient  = function(x) { -gradient(x) },
minHessian   = function(x) { -hessian(x)  }

))

