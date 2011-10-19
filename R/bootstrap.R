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
                                max.iter=1000L, verbose=FALSE, b.iter=0L,
                                Sigma.hat=NULL, Mu.hat=NULL, test=FALSE) {

    # verbose
    if(verbose) {
        if(b.iter == -1L) { 
            cat("  ... bootstrap draw: ")
        } else {
            cat("  ... bootstrap draw number: ", b.iter)
            b.iter <- b.iter + 1
        }
    }

    npar <- length(start)

    # 3. construct lavaan Sample (S4) object: description of the data
    bootSample <-
        try( Sample(data          = data,
                    group         = options$group,
                    sample.cov    = NULL,
                    sample.mean   = NULL,
                    sample.nobs   = NULL,
                    std.ov        = FALSE,
                    boot.idx      = boot.idx,
                    model.cov     = Sigma.hat,
                    model.mean    = Mu.hat,

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

    if(verbose) cat("     ok -- niter = ", attr(x, "iterations"), 
                    " fx = ", attr(x, "fx"), "\n")

    # strip attributes if not bollen.stine
    if(!test) {
        attributes(x) <- NULL
    } 

    x
}


bootstrapLRT.internal <- function(model1=NULL, model2=NULL, 
                                  sample=NULL, options=NULL,
                                  data=NULL, R=1000L,
                                  type="bollen.stine", 
                                  verbose=FALSE,
                                  max.iter=1000L) {

    # checks
    stopifnot(!is.null(model1), !is.null(model2),
              !is.null(sample), !is.null(options), 
              !is.null(data), 
              type %in% c("bollen.stine", "parametric"))

    # prepare
    LRT <- rep(as.numeric(NA), R)

    # bollen.stine or parametric: we need the Sigma.hat values
    Sigma.hat <- computeSigmaHat(model1)
    Mu.hat <- computeMuHat(model1) 

    # run bootstraps
    error.idx <- integer(0)
    for(b in 1:R) {
        if(type == "bollen.stine") {
            # take a bootstrap sample
            N <- sample@ntotal
            boot.idx <- sample(x=N, size=N, replace=TRUE)
        } else {
            # parametric!
            boot.idx <- NULL
        }

        # verbose
        if(verbose) cat("  ... bootstrap draw number: ", b, "\n")

        # create data
        bootSample <-
            try( Sample(data          = data,
                        group         = options$group,
                        sample.cov    = NULL,
                        sample.mean   = NULL,
                        sample.nobs   = NULL,
                        std.ov        = FALSE,
                        boot.idx      = boot.idx,
                        model.cov     = Sigma.hat,
                        model.mean    = Mu.hat,

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
            error.idx <- c(error.idx, b)
            next
        }

        # estimate model 1
        if(verbose) cat("  ... ... model 1: ")
        start1 <- getModelParameters(model1, type="free")
        lavaanModel1 <- setModelParameters(model1, x = start1)
        # switch of verbose in estimateModel
        verbose.old <- options$verbose; options$verbose <- FALSE
        x1 <- try( estimateModel(lavaanModel1,
                                 sample  = bootSample,
                                 options = options) )
        options$verbose <- verbose.old
        if(inherits(x1, "try-error")) {
            if(verbose) cat("     FAILED: in estimation\n")
            error.idx <- c(error.idx, b)
            next
        } else if(!attr(x1, "converged")) {
            if(verbose) cat("     FAILED: no convergence\n")
            error.idx <- c(error.idx, b)
            next
        } 
        if(verbose) cat("     ok -- niter = ", attr(x1, "iterations"),
                        " fx = ", attr(x1, "fx"), "\n")

        # estimate model 2
        if(verbose) cat("  ... ... model 2: ")
        start2 <- getModelParameters(model2, type="free")
        lavaanModel2 <- setModelParameters(model2, x = start2)
        # switch of verbose in estimateModel
        verbose.old <- options$verbose; options$verbose <- FALSE
        x2 <- try( estimateModel(lavaanModel2,
                                 sample  = bootSample,
                                 options = options) )
        options$verbose <- verbose.old
        if(inherits(x2, "try-error")) {
            if(verbose) cat("     FAILED: in estimation\n")
            error.idx <- c(error.idx, b)
            next
        } else if(!attr(x2, "converged")) {
            if(verbose) cat("     FAILED: no convergence\n")
            error.idx <- c(error.idx, b)
            next
        }
        if(verbose) cat("     ok -- niter = ", attr(x2, "iterations"),
                        " fx = ", attr(x2, "fx"), "\n")


        # store LRT
        NFAC <- 2 * unlist(sample@nobs)
        if(options$estimator == "ML" && options$likelihood == "wishart") {
            # first divide by two
            NFAC <- NFAC / 2
            NFAC <- NFAC - 1
            NFAC <- NFAC * 2
        }
        # model 1
        fx.group1 <- attr(attr(x1, "fx"), "fx.group")
        chisq.group1 <- fx.group1 * NFAC
        chisq.group1[which(chisq.group1 < 0)] <- 0.0
        chisq1 <- sum(chisq.group1)

        # model 2
        fx.group2 <- attr(attr(x2, "fx"), "fx.group")
        chisq.group2 <- fx.group2 * NFAC
        chisq.group2[which(chisq.group2 < 0)] <- 0.0
        chisq2 <- sum(chisq.group2)

        LRT[b] <- abs(chisq1 - chisq2)
    }

    # handle errors
    if(length(error.idx) > 0L) {
        warning("lavaan WARNING: only ", (R-length(error.idx)), " bootstrap draws were successful")
        LRT <- LRT[-error.idx]
        attr(LRT, "error.idx") <- error.idx
    } else {
        if(verbose) cat("Number of successful bootstrap draws:", 
                        (R - length(error.idx)), "\n")
    }

    LRT
}

bootstrap.internal <- function(model=NULL, sample=NULL, options=NULL,
                               data=NULL, R=1000L,
                               type="ordinary",
                               verbose=FALSE,
                               coef=TRUE,
                               test=FALSE,
                               FUN=NULL,
                               ...,
                               max.iter=1000L) {

    # checks
    stopifnot(!is.null(model), !is.null(sample), !is.null(options),
              !is.null(data),
              type %in% c("nonparametric", "ordinary",
                          "bollen.stine", "parametric"))
    if(type == "nonparametric") type <- "ordinary"

    # prepare
    start <- getModelParameters(model, type="free"); npar <- length(start)
    N <- sample@ntotal

    COEF <- TEST <- NULL
    if(coef) {
        COEF <- matrix(NA, R, npar)
    }
    if(test) {
        TEST <- rep(as.numeric(NA), R)
    }

    # bollen.stine or parametric: we need the Sigma.hat values
    if(type == "bollen.stine" || type == "parametric") {
        Sigma.hat <- computeSigmaHat(model)
        Mu.hat <- computeMuHat(model)
    } else {
        Sigma.hat <- NULL
        Mu.hat <- NULL
    }

    # run bootstraps
    error.idx <- integer(0)
    for(b in 1:R) {
        if(type == "bollen.stine" || type == "ordinary") {
            # take a bootstrap sample
            boot.idx <- sample(x=N, size=N, replace=TRUE)
        } else {
            # parametric!
            boot.idx <- NULL
        }

        # run bootstrap draw
        x <- lavaanBootStatistic(data=data, boot.idx=boot.idx, start=start,
                                 model=model, sample=sample, options=options,
                                 max.iter=max.iter, verbose=verbose, b.iter=b,
                                 Sigma.hat=Sigma.hat, Mu.hat=Mu.hat, test=test)

        # catch faulty run
        if(any(is.na(x))) {
            error.idx <- c(error.idx, b)
            next
        }

        if(coef) {
            # store parameters
            COEF[b,] <- x
        }

        # store test statistic
        if(test) {
            NFAC <- 2 * unlist(sample@nobs)
            if(options$estimator == "ML" && options$likelihood == "wishart") {
                # first divide by two
                NFAC <- NFAC / 2
                NFAC <- NFAC - 1
                NFAC <- NFAC * 2
            }
            fx.group <- attr(attr(x, "fx"), "fx.group")
            chisq.group <- fx.group * NFAC
            chisq.group[which(chisq.group < 0)] <- 0.0
            chisq <- sum(chisq.group)
            TEST[b] <- chisq
        }
    }

    # handle errors
    if(length(error.idx) > 0L) {
        warning("lavaan WARNING: only ", (R-length(error.idx)), " bootstrap draws were successful")
        COEF <- COEF[-error.idx,,drop=FALSE]
        #attr(COEF, "error.idx") <- error.idx

        if(test) {
            TEST <- TEST[-error.idx]
        }

    } else {
        if(verbose) cat("Number of successful bootstrap draws:",
                        (R - length(error.idx)), "\n")
    }

    list(coef=COEF, test=TEST)
}

# using the 'internal' function
bootstrapLavaan <- function(object, data=NULL, R=1000, type="ordinary", 
                            verbose=FALSE, coef=TRUE, test=FALSE,
                            FUN=NULL, ...) {

    # three types of information we may want to return:
    # 1) the coefficients (free only)
    # 2) the test statistic (chisq)
    # 3) any other type of information we can extract from a lavaan S4 object

    out <- bootstrap.internal(model=object@Model, 
                              sample=object@Sample,
                              options=object@Options, 
                              data=data, R=R,
                              verbose=verbose,
                              type=type,
                              coef=coef,
                              test=test,
                              FUN=FUN, ...)

    if(coef && !test && is.null(FUN)) {
        # coefficients only
        out <- out$coef
    } else if(!coef && test && is.null(FUN)) {
        out <- out$test
    } else if(coef && test && is.null(FUN)) {
        out <- list(coef=out$coef, test=out$test, FUN=NULL)
    } 

    out
}

bootstrapLRT <- function(h0, h1, data=NULL, R=1000, 
                         type="bollen.stine", verbose=FALSE) {

    out <- bootstrapLRT.internal(model1=h0@Model,
                                 model2=h1@Model,
                                 sample=h0@Sample,
                                 options=h0@Options,
                                 data=data, R=R,
                                 verbose=verbose,
                                 type=type)
    out
}


# using the 'boot' package
bootLavaan <- function(object, data=NULL, R=1000, ..., verbose=FALSE) {

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
        COEF.def <- apply(boot.out$t, 1, object@Model@def.function)
        if(length(def.idx) == 1L) {
            COEF.def <- as.matrix(COEF.def)
        } else {
            COEF.def <- t(COEF.def)
        }        
        colnames(COEF.def) <- object@User$lhs[def.idx]
        boot.out$t <- cbind(boot.out$t, COEF.def)
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
