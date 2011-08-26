# main function used by various bootstrap related functions
# this function draws the bootstrap samples, and estimates the
# free parameters for each bootstrap sample
#
# return BOOT matrix of size R x npar (R = number of bootstrap samples)
# 
# YR. 9 aug 2011
#
# Notes: - faulty runs are simply ignored (with a warning)
#        - default R=1000


lavaanBootStatistic <- function(data=NULL, boot.idx=NULL, start=NULL,
                                model=NULL, sample=NULL, options=NULL,
                                max.iter=1000L, verbose=FALSE, b.iter=0L) {

    # verbose
    if(verbose) {
        if(b.iter == -1L) { 
            cat("  ... bootstrap draw: ")
        } else {
            b.iter <- b.iter + 1
            cat("  ... bootstrap draw number: ", b)
        }
    }

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
        if(verbose) cat("     FAILED: in estimation\n")
        return(rep(NA, npar))
    } else if(!attr(x, "converged")) {
        if(verbose) cat("     FAILED: no convergence\n")
        return(rep(NA, npar))
    } else if(attr(x, "iterations") > max.iter) {
        # FIXME: is this wise? How should we choose max.iter??
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
    BOOT <- matrix(NA, R, npar)

    # run bootstraps
    error.idx <- integer(0)
    for(b in 1:R) {
        # take a bootstrap sample
        boot.idx <- sample(x=N, size=N, replace=TRUE)

        # run bootstrap draw
        x <- lavaanBootStatistic(data=data, boot.idx=boot.idx, start=start,
                                 model=model, sample=sample, options=options,
                                 max.iter=max.iter, verbose=verbose, b.iter=b)

        # catch faulty run
        if(any(is.na(x))) error.idx <- c(error.idx, b)    

        BOOT[b,] <- x
    }

    # handle errors
    if(length(error.idx) > 0L) {
        warning("lavaan WARNING: only ", (R-length(error.idx)), " bootstrap draws were successful")
        BOOT <- BOOT[-error.idx,]
        attr(BOOT, "error.idx") <- error.idx
    } else {
        if(verbose) cat("Number of successful bootstrap draws:", 
                        nrow(BOOT), "\n")
    }

    BOOT
}

# using the 'boot' package
lavaanBoot <- function(object, data=NULL, R=1000, ..., verbose=FALSE) {

    require("boot", quietly = TRUE)

    boot.out <- boot(data = data, 
                     statistic = lavaanBootStatistic, 
                     R = R,
                     ...,
                     start   = getModelParameters(object@Model, type="free"),
                     model   = object@Model, 
                     sample  = object@Sample, 
                     options = object@Options,
                     # no more than 4 times the number of iterations
                     # of the original run
                     max.iter = max(100, object@Fit@iterations*4), 
                     verbose = verbose, b.iter = -1L)

    # add rhs/op/lhs elements
    #free.idx <- which(object@User$free & !duplicated(object@User$free))
    #boot.out$rhs <- object@User$rhs[free.idx]
    #boot.out$op  <-  object@User$op[free.idx]
    #boot.out$lhs <- object@User$lhs[free.idx]
    attr(boot.out$t, "dimnames")[[2]] <-
        getParameterLabels(object@User, type="free")

    # adding defined parameters
    def.idx <- which(object@User$op == ":=")
    if(length(def.idx) > 0L) {
        boot.out$t0 <- c(boot.out$t0, object@Fit@est[def.idx])
        BOOT.def <- apply(boot.out$t, 1, object@Model@def.function)
        if(length(def.idx) == 1L) {
            BOOT.def <- as.matrix(BOOT.def)
        } else {
            BOOT.def <- t(BOOT.def)
        }        
        colnames(BOOT.def) <- object@User$lhs[def.idx]
        boot.out$t <- cbind(boot.out$t, BOOT.def)
        #boot.out$rhs <- c(boot.out$rhs, object@User$rhs[def.idx])
        #boot.out$op  <- c(boot.out$op,   object@User$op[def.idx])
        #boot.out$lhs <- c(boot.out$lhs, object@User$lhs[def.idx])
        
    }

    class(boot.out) <- c("lavaan.boot", "boot")
    boot.out
}


summary.lavaan.boot <- function(object=NULL, ..., type="perc", conf=0.95) {

    # catch 'estimated adjustment 'a' is NA' error
    # if R < nrow(data)
    if(type == "bca" && object$R < nrow(object$data)) {
        stop("lavaan ERROR: number of bootstrap draws must be greater than number of observations if type == \"bca\"")
    }


    # extract ci for one parameter
    ci <- function(index=1L) {
        out <- boot.ci(object, type=type, conf=conf, index=index)
        if(type == "norm") {
            ci <- out[[type, exact=FALSE]][c(2,3)]
        } else {
            ci <- out[[type, exact=FALSE]][c(4,5)]
        }
        ci
    }

    label <- colnames(object$t)
    est <- object$t0
    se <- apply(object$t, 2, sd)
    npar <- length(object$t0)
    CI <- sapply(seq_len(npar), ci)
    lower <- CI[1,]; upper <- CI[2,]

    #LIST <- data.frame(object$lhs, object$op, object$rhs, 
    #                   est, se, lower, upper)
    LIST <- data.frame(label, est, se, lower, upper)
    #names(LIST) <- c("lhs","op","rhs","est", "se", 
    names(LIST) <- c("label", "est", "se",
                     paste(type, ((1 - conf)/2*100), sep=""),
                     paste(type, ((conf+(1-conf)/2)*100), sep=""))
    class(LIST) <- c("lavaan.data.frame", "data.frame")
 
    LIST
}

# bootstrap parameters
bootstrapParameters <- function(object, data = NULL, R = 1000L,
                                verbose = TRUE, summarize = FALSE) {

    BOOT <- bootstrapParameters.internal(model   = object@Model, 
                                         sample  = object@Sample, 
                                         options = object@Options, 
                                         data    = data,
                                         R       = R, 
                                         verbose = verbose)

    # add defined parameters
    def.idx <- which(object@User$op == ":=")
    if(length(def.idx) > 0L) {
        BOOT.def <- apply(BOOT, 1, object@Model@def.function)
        if(length(def.idx) == 1L) {
            BOOT.def <- as.matrix(BOOT.def)
        } else {
            BOOT.def <- t(BOOT.def)
        }
        colnames(BOOT.def) <- object@User$lhs[def.idx]
    }
    

    if(summarize) {
        OUT <- parameterEstimates(object)
        OUT$est.std <- NULL
        OUT$est.std.all <- NULL
     
        boot.mean <- apply(BOOT, 2, mean)
        boot.sd   <- apply(BOOT, 2, sd) # fix N-1/N?

        GLIST <- x2GLIST(object@Model, x=boot.mean, type="free")
        OUT$est.boot.mean <- 
            getModelParameters(object@Model, GLIST=GLIST, type="user")
        
        GLIST <- x2GLIST(object@Model, x=boot.sd, type="free")
        OUT$est.boot.sd   <- 
            getModelParameters(object@Model, GLIST=GLIST, type="user")
        OUT$est.boot.sd[ which(object@User$free.uncon == 0L) ] <- 0.0
    } else {
        colnames(BOOT) <- names( coef(object) )
        OUT <- BOOT
    }

    OUT
}


