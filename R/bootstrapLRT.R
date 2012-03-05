# bootstrapLRT: bootstrap LR test for nested models

# revision YR 23 jan 2012: - use lavaan() for both model fits
#                          - add calibration

bootstrapLRT <- function(h0              = NULL,  # restricted model 
                         h1              = NULL,  # unrestricted model
                         R               = 1000L, 
                         type            = "bollen.stine", 
                         verbose         = FALSE,
                         return.LRT      = FALSE,
                         calibrate       = FALSE,
                         calibrate.R     = 1000L,
                         calibrate.alpha = 0.05,
                         warn            = -1L,
                         parallel        = c("no", "multicore", "snow"),
                         ncpus           = 1L,
                         cl              = NULL
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

    # parallel stuff (see boot() from the boot package)
    if (missing(parallel)) parallel <- "no"
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!have_mc && !have_snow) ncpus <- 1L
    }

    #if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    #seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
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
            X <- scale(data@X[[g]], center=TRUE, scale=FALSE)

            # transform
            X <- X %*% S.inv.sqrt %*% sigma.sqrt

            # add model-based mean
            if(h0@Model@meanstructure)
                X <- scale(X, center=(-1*h0@Sample@mean[[g]]), 
                              scale=FALSE)

            # replace data slot
            data@X[[g]] <- X
        }
    }

    # run bootstraps
    fn <- function(b) {
        if(type == "bollen.stine") {
            # take a bootstrap h0@Sample for each group
            boot.idx <- vector("list", length=h0@Sample@ngroups)
            for(g in 1:h0@Sample@ngroups) {
                stopifnot(h0@Sample@nobs[[g]] > 1L)
                boot.idx[[g]] <- sample(x=h0@Sample@nobs[[g]],
                                        size=h0@Sample@nobs[[g]], replace=TRUE)
            }
        } else {
            # parametric!
            boot.idx <- NULL
            for(g in 1:h0@Sample@ngroups) {
                data@X[[g]] <- MASS.mvrnorm(n     = h0@Sample@nobs[[g]],
                                            mu    = Mu.hat[[g]],
                                            Sigma = Sigma.hat[[g]])
            }
        }
        # names
        for(g in 1:h0@Sample@ngroups) 
            colnames(data@X[[g]]) <- h0@Data@ov.names[[g]]
 
        # verbose
        if(verbose) cat("  ... bootstrap draw number: ", b, "\n")

        # take care of bootstrap h0@Sample statistics
        Missing <-  getMissingPatterns(Data    = data,
                                       missing = h0@Options$missing,
                                       warn    = FALSE,
                                       verbose = FALSE)
        WLS.V <- list()
        if(h0@Options$estimator %in% c("GLS", "WLS")) {
            WLS.V <- getWLS.V(Data          = data,
                              sample        = NULL,
                              boot.idx      = boot.idx,
                              estimator     = h0@Options$estimator,
                              mimic         = h0@Options$mimic,
                              meanstructure = h0@Options$meanstructure)
        }
        bootSampleStats <- try(getSampleStatsFromData(
                               Data        = data,
                               M           = Missing,
                               boot.idx    = boot.idx,
                               rescale     = (h0@Options$estimator == "ML" &&
                                              h0@Options$likelihood == "normal"),
                               WLS.V       = WLS.V)) # fixme!!
        if(inherits(bootSampleStats, "try-error")) {
            if(verbose) cat("     FAILED: creating h0@Sample statistics\n")
            return(NULL)
        }

        # adjust h0@Model slot if fixed.x variances/covariances
        # have changed:
      ### FIXME #####
        #if(h0@Model@fixed.x && length(vnames(fit0@User, "ov.x")) > 0L) {
        #    for(g in 1:h0@Sample@ngroups) {
        #        
        #    }
        #}

        # just in case we need the dataSlot (lm!)
        if(type == "bollen.stine") {
            data.boot <- data
            for(g in 1:h0@Sample@ngroups) {
                data.boot@X[[g]] <- data@X[[g]][ boot.idx[[g]],,drop=FALSE]
            }
        } else {
            data.boot <- data
        }

        # estimate model 1
        if(verbose) cat("  ... ... model h0: ")
        h0@Options$verbose <- FALSE
        h0@Options$se <- "none"; h0@Options$test <- "standard"
        fit.h0 <- lavaan(slotOptions = h0@Options,
                         slotUser    = h0@User,
                         #slotModel   = h0@Model, # only if fixed.x=FALSE???
                         slotSample  = bootSampleStats,
                         slotData    = data.boot)
        if(!fit.h0@Fit@converged) {
            if(verbose) cat("     FAILED: no convergence\n")
            return(NULL)
        } 
        if(verbose) cat("     ok -- niter = ", fit.h0@Fit@iterations,
                        " fx = ", fit.h0@Fit@fx, "\n")

        # estimate model 2
        if(verbose) cat("  ... ... model h1: ")
        h1@Options$verbose <- FALSE
        h1@Options$se <- "none"; h1@Options$test <- "standard"
        fit.h1 <- lavaan(slotOptions = h1@Options,
                         slotUser    = h1@User,
                         #slotModel   = h1@Model, # only if fixed.x=FALSE???
                         slotSample  = bootSampleStats,
                         slotData    = data.boot)
        if(!fit.h1@Fit@converged) {
            if(verbose) cat("     FAILED: no convergence\n")
            return(NULL)
        } 
        if(verbose) cat("     ok -- niter = ", fit.h1@Fit@iterations,
                        " fx = ", fit.h1@Fit@fx, "\n")


        # store LRT
        lrt.boot <- abs(anova(fit.h1, fit.h0)$`Chisq diff`[2L])
        if(verbose) cat("  ... ... LRT = ", lrt.boot,"\n")

         # calibration
        if(calibrate) {
            if(verbose) cat("  ... ... calibrating p.value - ");
            
            plugin.pvalue <- bootstrapLRT(h0=fit.h0, h1=fit.h1,
                                          R=calibrate.R,
                                          type=type,
                                          verbose=FALSE,
                                          calibrate=FALSE,
                                          return.LRT=FALSE)
            if(verbose) cat(sprintf("%5.3f", plugin.pvalue),"\n")
            attr(lrt.boot, "plugin.pvalue") <- plugin.pvalue
        }

        lrt.boot
    } # b

    # this is from the boot function in package boot
    RR <- sum(R)
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
        if (have_mc) {
            parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
        } else if (have_snow) {
            #list(...) # evaluate any promises
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                if(RNGkind()[1L] == "L'Ecuyer-CMRG")
                    parallel::clusterSetRNGStream(cl)
                res <- parallel::parLapply(cl, seq_len(RR), fn)
                parallel::stopCluster(cl)
                res
            } else parallel::parLapply(cl, seq_len(RR), fn)
        }
    } else lapply(seq_len(RR), fn)

    # handle errors and fill in container
    error.idx <- integer(0)
    for(b in seq_len(RR)) {
        if(!is.null(res[[b]])) {
            LRT[b] <- res[[b]]
            if(calibrate) 
                plugin.pvalues[b] <- attr(res[[b]], "plugin.pvalue")
        } else {
            error.idx <- c(error.idx, b)
        }
    }

    # handle errors
    if(length(error.idx) > 0L) {
        warning("lavaan WARNING: only ", (R-length(error.idx)), " bootstrap draws were successful")
        LRT <- LRT[-error.idx]
        if(calibrate) plugin.pvalues <- plugin.pvalues[-error.idx]
        attr(LRT, "error.idx") <- error.idx
    } else {
        if(verbose) cat("Number of successful bootstrap draws:", 
                        (R - length(error.idx)), "\n")
    }

    pvalue <- sum(LRT > LRT.original)/length(LRT)
    if(return.LRT) { 
        attr(pvalue, "LRT") <- LRT
        attr(pvalue, "LRT.original") <- LRT.original
    }
     
    if(calibrate) {
        adj.alpha <- quantile(plugin.pvalues, calibrate.alpha)
        attr(pvalue, "plugin.pvalues") <- plugin.pvalues
        attr(pvalue, "adj.alpha") <- adj.alpha
    }

    pvalue
}

