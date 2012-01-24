# bootstrapLRT: bootstrap LR test for nested models

# revision YR 23 jan 2012: - use lavaan() for both model fits
#                          - add calibration

bootstrapLRT <- function(h0              = NULL,  # restricted model 
                         h1              = NULL,  # unrestrited model
                         R               = 1000, 
                         type            = "bollen.stine", 
                         verbose         = FALSE,
                         return.LRT      = FALSE,
                         calibrate       = FALSE,
                         calibrate.R     = 1000,
                         calibrate.alpha = 0.05,
                         warn            = 1L,
                        ) {

    # checks
    stopifnot(class(h0) == "lavaan", class(h1) == "lavaan",
              type %in% c("bollen.stine", "parametric"))

    # set warning level
    options(warn = warn)

    # prepare
    LRT <- rep(as.numeric(NA), R)
    LRT.original <- abs(anova(h0, h1)$`Chisq diff`[2L])
    if(calibrate) plugin.pvalues <- numeric(R)
    
    # data
    data <- h0@Data

    # bollen.stine or parametric: we need the Sigma.hat values
    Sigma.hat <- computeSigmaHat(h0@Model)
    Mu.hat    <- computeMuHat(   h0@Model) 

    # if bollen.stine, transform data here
    if(type == "bollen.stine") {
        for(g in 1:h0@Sample@ngroups) {
            sigma.sqrt <- sqrtSymmetricMatrix(     Sigma.hat[[g]])
            S.inv.sqrt <- sqrtSymmetricMatrix(h0@Sample@icov[[g]])

            # center (needed???)
            X <- scale(data[[g]], center=TRUE, scale=FALSE)

            # transform
            X <- X %*% S.inv.sqrt %*% sigma.sqrt

            # add model-based mean
            if(h0@Model@meanstructure)
                X <- scale(X, center=(-1*h0@Sample@mean[[g]]), scale=FALSE)

            # replace data slot
            data[[g]] <- X
        }
    }

    # run bootstraps
    error.idx <- integer(0)
    for(b in 1:R) {
        if(type == "bollen.stine") {
            # take a bootstrap h0@Sample for each group
            boot.idx <- vector("list", length=h0@Sample@ngroups)
            for(g in 1:h0@Sample@ngroups)
                boot.idx[[g]] <- sample(x=h0@Sample@nobs[[g]],
                                        size=h0@Sample@nobs[[g]], replace=TRUE)
        } else {
            # parametric!
            boot.idx <- NULL
            for(g in 1:h0@Sample@ngroups) {
                data[[g]] <- MASS.mvrnorm(n     = h0@Sample@nobs[[g]],
                                          mu    = Mu.hat[[g]],
                                          Sigma = Sigma.hat[[g]])
            }
        }
        colnames(data[[g]]) <- h0@Sample@ov.names[[g]]

        # verbose
        if(verbose) cat("  ... bootstrap draw number: ", b, "\n")

        # take care of bootstrap h0@Sample statistics
        Missing <-  getMissingPatterns(X       = data,
                                       missing = h0@Options$missing,
                                       warn    = FALSE,
                                       verbose = FALSE)
        WLS.V <- list()
        if(h0@Options$estimator %in% c("GLS", "WLS")) {
            WLS.V <- getWLS.V(X             = data,
                              sample        = NULL,
                              boot.idx      = boot.idx,
                              estimator     = h0@Options$estimator,
                              mimic         = h0@Options$mimic,
                              meanstructure = h0@Options$meanstructure)
        }
        bootSampleStats <- try(getSampleStatsFromData(
                               X           = data,
                               M           = Missing,
                               boot.idx    = boot.idx,
                               rescale     = (h0@Options$estimator == "ML" &&
                                              h0@Options$likelihood == "normal"),
                               group.label = h0@Sample@group.label,
                               WLS.V       = WLS.V)) # fixme!!
        if(inherits(bootSampleStats, "try-error")) {
            if(verbose) cat("     FAILED: creating h0@Sample statistics\n")
            error.idx <- c(error.idx, b)
            next
        }


        # estimate model 1
        if(verbose) cat("  ... ... model h0: ")
        h0@Options$verbose <- FALSE
        h0@Options$se <- "none"
        h0@Options$test <- "standard"
        fit.h0 <- lavaan(slotOptions = h0@Options,
                         slotUser    = h0@User,
                        #slotModel   = h0@Model, # only if fixed.x=FALSE???
                         slotSample  = bootSampleStats,
                         slotData    = data)
        if(!fit.h0@Fit@converged) {
            if(verbose) cat("     FAILED: no convergence\n")
            error.idx <- c(error.idx, b)
            next
        } 
        if(verbose) cat("     ok -- niter = ", fit.h0@Fit@iterations,
                        " fx = ", fit.h0@Fit@fx, "\n")

        # estimate model 2
        if(verbose) cat("  ... ... model h1: ")
        h0@Options$verbose <- FALSE
        h0@Options$se <- "none"
        h0@Options$test <- "standard"
        fit.h1 <- lavaan(slotOptions = h0@Options,
                         slotUser    = h1@User,
                        #slotModel   = h0@Model, # only if fixed.x=FALSE???
                         slotSample  = bootSampleStats,
                         slotData    = data)
        if(!fit.h1@Fit@converged) {
            if(verbose) cat("     FAILED: no convergence\n")
            error.idx <- c(error.idx, b)
            next
        } 
        if(verbose) cat("     ok -- niter = ", fit.h1@Fit@iterations,
                        " fx = ", fit.h1@Fit@fx, "\n")


        # store LRT
        LRT[b] <- abs(anova(fit.h1, fit.h0)$`Chisq diff`[2L])
        if(verbose) cat("  ... ... LRT = ", LRT[b],"\n")

        # calibration
        if(calibrate) {
            if(verbose) cat("  ... ... calibrating p.value - ");
            plugin.pvalues[b] <- bootstrapLRT(h0=fit.h0, h1=fit.h1,
                                                 R=calibrate.R,
                                                 type=type, verbose=FALSE,
                                                 calibrate=FALSE,
                                                 return.LRT=FALSE)
            if(verbose) cat(sprintf("%5.3f", plugin.pvalues[b]),"\n")
        }

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

    pvalue <- sum(LRT >= LRT.original)/length(LRT)
    if(return.LRT) { 
        attr(pvalue, "LRT") <- LRT
        attr(pvalue, "LRT.original") <- LRT.original
    }
     
    if(calibrate) {
        adj.alpha <- quantile(plugin.pvalues, alpha)
        attr(pvalue, "plugin.pvalues") <- plugin.pvalues
        attr(pvalue, "adj.alpha") <- adj.alpha
    }

    pvalue
}

