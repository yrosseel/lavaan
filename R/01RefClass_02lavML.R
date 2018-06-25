# Maximum Likelihood optimization -- YR 10 july 2012

# super class -- virtual statistical model that needs to be optimized
lavRefML <- setRefClass("lavML",

# inherits
contains = "lavOptim",

# fields
fields = list(
    y = "numeric",        # the (unidimensional) data
    nobs = "integer",     # number of observations
    weights = "numeric"   # weights
),

# methods
methods = list(

logl = function(x) {
    if(!missing(x)) .self$theta <- x
    likelihoods <- lik()
    # FIXME: handle zero/negative/small likelihood values
    sum(log(likelihoods), na.rm=TRUE)
},

lik = function(x) {
    if(!missing(x)) .self$theta <- x
    cat("this is dummy function\n")
    return(rep(as.numeric(NA), nobs))
},

scores = function(x) {
    if(!missing(x)) .self$theta <- x
    cat("this is dummy function\n")
    return(matrix(as.numeric(NA), nobs, npar))
},

gradient = function(x) {
    SCORES <- scores(x)
    apply(SCORES, 2L, base::sum, na.rm=TRUE)
}

))
