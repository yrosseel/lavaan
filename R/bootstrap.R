# main function used by various bootstrap related functions
# this function draws the bootstrap samples, and estimates the
# free parameters for each bootstrap sample
#
# return COEF matrix of size R x npar (R = number of bootstrap samples)
# 
# YR. 9 aug 2011
#
# Notes: - faulty runs are simply ignored (with a warning)
#        - default R=1000


lavaanBootStatistic <- function(data=NULL, boot.idx=NULL, start=NULL,
                                model=NULL, sample=NULL, options=NULL,
                                max.iter=1000L, verbose=FALSE, b=0L) {

    b <- b + 1
    if(verbose) cat("  ... bootstrap draw number: ", b)
    npar <- length(start)

    # 3. construct lavaan Sample (S4) object: description of the data
    bootSample <-
        try( Sample(data          = data[boot.idx,],
                    group         = options$group,
                    sample.cov    = NULL,
                    sample.mean   = NULL,
                    sample.nobs   = NULL,
                    std.ov        = FALSE,

                    ov.names      = sample@ov.names,
                    data.type     = "full",
                    ngroups       = sample@ngroups,
                    group.label   = sample@group.label,
                    estimator     = options$estimator,
                    likelihood    = options$likelihood,
                    mimic         = options$mimic,
                    meanstructure = options$meanstructure,
                    missing       = options$missing,
                    warn          = options$warn,
                    verbose       = options$verbose) )
    if(inherits(bootSample, "try-error")) {
        error.idx <- c(error.idx, b)
        if(verbose) cat("     FAILED: creating sample statistics\n")
        return(rep(NA, npar))
    }

    # 6. estimate free parameters
    lavaanModel <- setModelParameters(model, x = start)
    # switch of verbose in estimateModel
    verbose.old <- options$verbose; options$verbose <- FALSE
    x <- try( estimateModel(lavaanModel,
                            sample  = bootSample,
                            options = options) )
    options$verbose <- verbose.old
    if(inherits(x, "try-error")) {
        error.idx <- c(error.idx, b)
        if(verbose) cat("     FAILED: in estimation\n")
        return(rep(NA, npar))
    } else if(!attr(x, "converged")) {
        error.idx <- c(error.idx, b)
        if(verbose) cat("     FAILED: no convergence\n")
        return(rep(NA, npar))
    } else if(attr(x, "iterations") > max.iter) {
        # FIXME: is this wise? How should we choose max.iter??
        error.idx <- c(error.idx, b)
        if(verbose) cat("     FAILED: too many iterations\n")
        return(rep(NA, npar))
    }

    if(verbose) cat("     ok -- niter = ", attr(x, "iterations"), "\n")

    # strip attributes
    attributes(x) <- NULL

    x
}


# Notes: - faulty runs are REMOVED!
#        - default R=1000
bootstrapParameters.internal <- function(model=NULL, sample=NULL, options=NULL,
                                         data=NULL, R=1000L, verbose=FALSE,
                                         max.iter=1000L) {

    # checks
    stopifnot(!is.null(model), !is.null(sample), !is.null(options), 
              !is.null(data))

    # prepare
    start <- getModelParameters(model, type="free"); npar <- length(start)
    N <- sample@ntotal
    COEF <- matrix(NA, R, npar)

    # run bootstraps
    error.idx <- integer(0)
    for(b in 1:R) {
        # take a bootstrap sample
        boot.idx <- sample(x=N, size=N, replace=TRUE)

        # run bootstrap draw
        x <- lavaanBootStatistic(data=data, boot.idx=boot.idx, start=start,
                                 model=model, sample=sample, options=options,
                                 max.iter=max.iter, verbose=verbose, b=b)

        # catch faulty run
        if(any(is.na(x))) error.idx <- c(error.idx, b)    

        COEF[b,] <- x
    }

    # handle errors
    if(length(error.idx) > 0L) {
        warning("lavaan WARNING: only ", (R-length(error.idx)), " bootstrap draws were successful")
        COEF <- COEF[-error.idx,]
        attr(COEF, "error.idx") <- error.idx
    } else {
        if(verbose) cat("Number of successful bootstrap draws:", 
                        nrow(COEF), "\n")
    }

    COEF
}

# using the 'boot' package
boot.lavaan <- function(object, data=NULL, ...) {

    require("boot", quietly = TRUE)

    boot.out <- boot:::boot(data = data, 
                       statistic = lavaanBootStatistic, 
                       ...,
                       start   = getModelParameters(object@Model, type="free"),
                       model   = object@Model, 
                       sample  = object@Sample, 
                       options = object@Options,
                       # no more than 4 times the number of iterations
                       # than the original run
                       max.iter = max(100, object@Fit@iterations*4), 
                       verbose = FALSE, b = 0L)

    boot.out
}

# bootstrap parameters
bootstrapParameters <- function(object, data = NULL, R = 1000L,
                                verbose = TRUE, summarize = FALSE) {

    COEF <- bootstrapParameters.internal(model   = object@Model, 
                                         sample  = object@Sample, 
                                         options = object@Options, 
                                         data    = data,
                                         R       = R, 
                                         verbose = verbose)

    if(summarize) {
        OUT <- parameterEstimates(object)
        OUT$est.std <- NULL
        OUT$est.std.all <- NULL
     
        boot.mean <- apply(COEF, 2, mean)
        boot.sd   <- apply(COEF, 2, sd) # fix N-1/N?

        GLIST <- x2GLIST(object@Model, x=boot.mean, type="free")
        OUT$est.boot.mean <- 
            getModelParameters(object@Model, GLIST=GLIST, type="user")
        
        GLIST <- x2GLIST(object@Model, x=boot.sd, type="free")
        OUT$est.boot.sd   <- 
            getModelParameters(object@Model, GLIST=GLIST, type="user")
        OUT$est.boot.sd[ which(object@User$free.uncon == 0L) ] <- 0.0
    } else {
        colnames(COEF) <- names( coef(object) )
        OUT <- COEF
    }

    OUT
}

