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
# 
# Updates: - now we have a separate @data slot, we only need to transform once
#            for the bollen.stine bootstrap (13 dec 2011)
#          - bug fix: we need to 'update' the fixed.x variances/covariances
#            for each bootstrap draw!!
#
# Question: if fixed.x=TRUE, should we not keep X fixed, and bootstrap Y
#           only, conditional on X?? How to implement the conditional part?

lavaanBootStatistic <- function(data=NULL, boot.idx=NULL, user=NULL,
                                model=NULL, sample=NULL, options=NULL,
                                verbose=FALSE, b.iter=0L, control=NULL,
                                test=FALSE) {

    # verbose
    if(verbose) {
        if(b.iter == -1L) { 
            cat("  ... bootstrap draw: ")
        } else {
            cat("  ... bootstrap draw number: ", sprintf("%5d", b.iter))
            b.iter <- b.iter + 1
        }
    }

    npar <- model@nx.free

    # 3. construct lavaan Sample (S4) object: description of the data
    # FIXME!!!
    Missing <-  getMissingPatterns(X       = data,
                                   missing = options$missing,
                                   warn    = FALSE,
                                   verbose = FALSE)
    WLS.V <- list()
    if(options$estimator %in% c("GLS", "WLS")) {
        WLS.V <- getWLS.V(X             = data,
                          sample        = NULL,
                          boot.idx      = boot.idx,
                          estimator     = options$estimator,
                          mimic         = options$mimic,
                          meanstructure = options$meanstructure)
    }
    bootSampleStats <- try(getSampleStatsFromData(
                               X           = data,
                               M           = Missing,
                               boot.idx    = boot.idx,
                               rescale     = (options$estimator == "ML" &&
                                              options$likelihood == "normal"),
                               group.label = sample@group.label,
                               WLS.V       = WLS.V)) # fixme!!
    if(inherits(bootSampleStats, "try-error")) {
        if(verbose) cat("     FAILED: creating sample statistics\n")
        return(rep(NA, npar))
    }

    # 5. fill in starting values in Model (free parameters only!)
    #lavaanModel <- setModelParameters(model, x = start)
    lavaanStart <-
        StartingValues(start.method = "default",
                       user         = user,
                       sample       = bootSampleStats,
                       model.type   = options$model.type,
                       mimic        = options$mimic,
                       debug        = options$debug)

    # 5. construct internal model (S4) representation
    lavaanModel <-
        Model(user           = user,
              start          = lavaanStart,
              representation = options$representation,
              debug          = options$debug)

    # 6. estimate free parameters
    # switch of verbose in estimateModel
    #verbose.old <- options$verbose; options$verbose <- FALSE
    x <- try( estimateModel(lavaanModel,
                            sample  = bootSampleStats,
                            control = control,
                            options = options) )
    #options$verbose <- verbose.old
    if(inherits(x, "try-error")) {
        if(verbose) cat("     FAILED: in estimation\n")
        return(rep(NA, npar))
    } else if(!attr(x, "converged")) {
        if(verbose) cat("     FAILED: no convergence\n")

        # DEBUG
        #write.csv(data, file="BDEBUG.TXT", row.names=FALSE)
        #str(options)
        #str(user)
        #str(bootSampleStats)
        #str(lavaanModel)
        #cat("x = \n"); print(x)
        #cat("starting values: ", lavaanStart, "\n")
        #stop("for now")

        return(rep(NA, npar))
    }

    if(verbose) cat("     ok -- niter = ", 
                    sprintf("%3d", attr(x, "iterations")), 
                    " fx = ", sprintf("%13.9f", attr(x, "fx")), "\n")

    # strip attributes if not bollen.stine
    if(!test) {
        attributes(x) <- NULL
    } 

    x
}


bootstrap.internal <- function(model=NULL, sample=NULL, options=NULL,
                               data=NULL, user=NULL, R=1000L,
                               type="ordinary",
                               verbose=FALSE,
                               coef=TRUE,
                               test=FALSE,
                               FUN=NULL,
                               ...,
                               control=list()) {

    # checks
    stopifnot(!is.null(model), !is.null(sample), !is.null(options),
              !is.null(data), !is.null(user),
              type %in% c("nonparametric", "ordinary",
                          "bollen.stine", "parametric"))
    if(type == "nonparametric") type <- "ordinary"

    # prepare
    npar <- model@nx.free
    
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

    # if bollen.stine, transform data here
    if(type == "bollen.stine") {
        for(g in 1:sample@ngroups) {
            sigma.sqrt <- sqrtSymmetricMatrix(  Sigma.hat[[g]])
            S.inv.sqrt <- sqrtSymmetricMatrix(sample@icov[[g]]) 

            # center (needed???)
            X <- scale(data[[g]], center=TRUE, scale=FALSE)
    
            # transform
            X <- X %*% S.inv.sqrt %*% sigma.sqrt

            # add model-based mean
            if(model@meanstructure)
                X <- scale(X, center=(-1*sample@mean[[g]]), scale=FALSE)

            # replace data slot
            data[[g]] <- X
        }
    }

    # run bootstraps
    error.idx <- integer(0)
    for(b in 1:R) {
        if(type == "bollen.stine" || type == "ordinary") {
            # take a bootstrap sample for each group
            boot.idx <- vector("list", length=sample@ngroups)
            for(g in 1:sample@ngroups) 
                boot.idx[[g]] <- sample(x=sample@nobs[[g]], 
                                        size=sample@nobs[[g]], replace=TRUE)
        } else { # parametric!
            boot.idx <- NULL
            for(g in 1:sample@ngroups) {
                    data[[g]] <- MASS.mvrnorm(n     = sample@nobs[[g]],
                                              Sigma = Sigma.hat[[g]],
                                              mu    = Mu.hat[[g]])
            }
        }

        # run bootstrap draw
        x <- lavaanBootStatistic(data=data, boot.idx=boot.idx,
                                 model=model, sample=sample, options=options,
                                 user=user, 
                                 verbose=verbose, b.iter=b, test=test,
                                 control=control)

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
bootstrapLavaan <- function(object, R=1000, type="ordinary", 
                            verbose=FALSE, coef=TRUE, test=FALSE,
                            FUN=NULL, ...) {

    # three types of information we may want to return:
    # 1) the coefficients (free only)
    # 2) the test statistic (chisq)
    # 3) any other type of information we can extract from a lavaan S4 object

    out <- bootstrap.internal(model   = object@Model, 
                              sample  = object@Sample,
                              options = object@Options, 
                              data    = object@Data, 
                              user    = object@User,
                              control = object@Fit@control,
                              R       = R,
                              verbose = verbose,
                              type    = type,
                              coef    = coef,
                              test    = test,
                              FUN     = FUN, 
                              ...)

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
