# bootstrapLRT: bootstrap LR test for nested models

bootstrapLRT <- function(h0, h1, data=NULL, R=1000, 
                         type="bollen.stine", verbose=FALSE) {

    out <- bootstrapLRT.internal(model1  = h0@Model,
                                 model2  = h1@Model,
                                 sample  = h0@Sample,
                                 options = h0@Options,
                                 data    = h0@Data, 
                                 R       = R,
                                 verbose = verbose,
                                 type=type)
    out
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
            if(model1@meanstructure)
                X <- scale(X, center=(-1*sample@mean[[g]]), scale=FALSE)

            # replace data slot
            data[[g]] <- X
        }
    }

    # run bootstraps
    error.idx <- integer(0)
    for(b in 1:R) {
        if(type == "bollen.stine") {
            # take a bootstrap sample for each group
            boot.idx <- vector("list", length=sample@ngroups)
            for(g in 1:sample@ngroups)
                boot.idx[[g]] <- sample(x=sample@nobs[[g]],
                                        size=sample@nobs[[g]], replace=TRUE)
        } else {
            # parametric!
            boot.idx <- NULL
            for(g in 1:sample@ngroups) {
                data[[g]] <- MASS.mvrnorm(n     = sample@nobs[[g]],
                                          mu    = Sigma.hat[[g]],
                                          Sigma = Mu.hat[[g]])
            }
        }

        # verbose
        if(verbose) cat("  ... bootstrap draw number: ", b, "\n")

        # take care of bootstrap sample statistics
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
                                 sample  = bootSampleStats,
                                 # control??
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
                                 sample  = bootSampleStats,
                                 # control??
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

    } # b

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

