bootstrapLRT <- function (h0 = NULL, h1 = NULL, R = NULL, 
                          type = "bollen.stine", verbose = FALSE,
                          return.LRT = TRUE, 
                          double.bootstrap = "FDB",
                          double.bootstrap.R = 500L, 
                          double.bootstrap.alpha = 0.05, 
                          warn = -1L, parallel = c("no", "multicore", "snow"), 
                          ncpus = 1L, cl = NULL) 
{
    # checks
    stopifnot(class(h0) == "lavaan", 
              class(h1) == "lavaan", 
              type %in% c("bollen.stine", "parametric"), 
              double.bootstrap %in% c("no", "FDB", "standard"))
  
    options(warn = warn)

    # prepare
    LRT <- rep(as.numeric(NA), R)
    if((h1@Fit@fx - h0@Fit@fx) > (.Machine$double.eps * 10)) { 
        # restricted fit should not be better!
        cat(" ... h0@Fit@fx = ", h0@Fit@fx, "h1@Fit@fx = ", h1@Fit@fx,
            "h0 should not be better!\n")
        return(NULL)
    }
    LRT.original <- abs(anova(h0, h1)$`Chisq diff`[2L])
    # abs only needed because df may be the same for both models!
  
  
    if(double.bootstrap == "FDB") {
        LRT.2 <- numeric(R)
    } else if(double.bootstrap == "standard") {
        plugin.pvalues <- numeric(R)
    }
  
  
    # prepare for parallel processing
    if(missing(parallel)) parallel <- "no"
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if(parallel != "no" && ncpus > 1L) {
        if(parallel == "multicore") 
            have_mc <- .Platform$OS.type != "windows"
        else if(parallel == "snow") 
            have_snow <- TRUE
        if(!have_mc && !have_snow) 
            ncpus <- 1L
    }

    #data
    data <- h0@Data
  
    #Compute covariance matrix and additional mean vector
    Sigma.hat <- computeSigmaHat(h0@Model)
    Mu.hat    <- computeMuHat(h0@Model)
  
    #Bollen-Stine data transformation
    if (type == "bollen.stine") {
        for (g in 1:h0@Sample@ngroups) {
            sigma.sqrt <- sqrtSymmetricMatrix(     Sigma.hat[[g]])
            S.inv.sqrt <- sqrtSymmetricMatrix(h0@Sample@icov[[g]])

            #center
            X <- scale(data@X[[g]], center = TRUE, scale = FALSE)

            #transform
            X <- X %*% S.inv.sqrt %*% sigma.sqrt

            #add model based mean
            if (h0@Model@meanstructure) 
                X <- scale(X, center = (-1 * h0@Sample@mean[[g]]), scale = FALSE)

            #replace data slot
            data@X[[g]] <- X
        }
    }

    fn <- function(b) {

        #Sampling if Bollen-Stine
        if (type == "bollen.stine") {
            boot.idx <- vector("list", length = h0@Sample@ngroups)
            for (g in 1:h0@Sample@ngroups) {
                stopifnot(h0@Sample@nobs[[g]] > 1L)
                boot.idx[[g]] <- sample(x = h0@Sample@nobs[[g]], 
                                        size = h0@Sample@nobs[[g]], 
                                        replace = TRUE)
            }
        } else {
            #Parametric bootstrap
            boot.idx <- NULL
            for (g in 1:h0@Sample@ngroups) {
                data@X[[g]] <- MASS.mvrnorm(n = h0@Sample@nobs[[g]], 
                                            mu = Mu.hat[[g]], 
                                            Sigma = Sigma.hat[[g]])
            }
        }

        #Add variable names to columns
        for (g in 1:h0@Sample@ngroups) 
            colnames(data@X[[g]]) <- h0@Data@ov.names[[g]]

        if (verbose) 
            cat("  ... bootstrap draw number: ", b, "\n")

        #Missing data handling
        Missing <- getMissingPatterns(Data = data, 
                                      missing = h0@Options$missing, 
                                      warn = FALSE, 
                                      verbose = FALSE)
        WLS.V <- list()

        if (h0@Options$estimator %in% c("GLS", "WLS")) {
            WLS.V <- getWLS.V(Data = data, 
                              sample = NULL, 
                              boot.idx = boot.idx, 
                              estimator = h0@Options$estimator, 
                              mimic = h0@Options$mimic, 
                              meanstructure = h0@Options$meanstructure)
        }

        #Get sample statistics
        bootSampleStats <- 
            try(getSampleStatsFromData(Data = data, 
                                       M = Missing, 
                                       boot.idx = boot.idx, 
                                       rescale = (h0@Options$estimator == "ML" && h0@Options$likelihood =="normal"), 
                                       WLS.V = WLS.V))

        if (inherits(bootSampleStats, "try-error")) {
            if (verbose) 
                cat("     FAILED: creating h0@Sample statistics\n")
            return(NULL)
        }

        #Bollen-Stine sampling
        if (type == "bollen.stine") {
            data.boot <- data
            for (g in 1:h0@Sample@ngroups) {
                data.boot@X[[g]] <- data@X[[g]][boot.idx[[g]], ,drop = FALSE]
            }
        } else {
            data.boot <- data
        }

        if (verbose) cat("  ... ... model h0: ")
        h0@Options$verbose <- FALSE
        h0@Options$se <- "none"
        h0@Options$test <- "standard"

        #Fit h0 model
        fit.h0 <- lavaan(slotOptions = h0@Options,
                         slotUser = h0@User, 
                         slotSample = bootSampleStats, 
                         slotData = data.boot)

        if (!fit.h0@Fit@converged) {
            if (verbose) cat("     FAILED: no convergence\n")
            return(NULL)
        }
        if (verbose) 
            cat("     ok -- niter = ", fit.h0@Fit@iterations, 
                " fx = ", fit.h0@Fit@fx, "\n")
        if (verbose) cat("  ... ... model h1: ")
        h1@Options$verbose <- FALSE
        h1@Options$se <- "none"
        h1@Options$test <- "standard"

        #Fit h1 model
        fit.h1 <- lavaan(slotOptions = h1@Options, 
                         slotUser = h1@User, 
                         slotSample = bootSampleStats, 
                         slotData = data.boot)

        if (!fit.h1@Fit@converged) {
            if (verbose) 
                cat("     FAILED: no convergence -- niter = ", fit.h1@Fit@iterations, 
                    " fx = ", fit.h1@Fit@fx,"\n")
            return(NULL)
        }
        if (verbose) 
            cat("     ok -- niter = ", fit.h1@Fit@iterations, 
                " fx = ", fit.h1@Fit@fx, "\n")

        # store LRT
        if((fit.h1@Fit@fx - fit.h0@Fit@fx) > (.Machine$double.eps * 10)) {
            #if((fit.h1@Fit@fx - fit.h0@Fit@fx) > 0.0) {
            if (verbose)
                cat("  ... ... LRT  = <NA> h0 > h1, delta = ", fit.h1@Fit@fx - fit.h0@Fit@fx, "\n")
            return(NULL)
        } else {
            lrt.boot <- abs(anova(fit.h1, fit.h0)$`Chisq diff`[2L])
            if (verbose)
                cat("  ... ... LRT  = ", lrt.boot, "\n")
        }
     
        #double bootstrap
        if (double.bootstrap == "standard") {
            if (verbose) cat("  ... ... calibrating p.value - ")

            plugin.pvalue <- bootstrapLRT(h0 = fit.h0, h1 = fit.h1, 
                                          R = double.bootstrap.R, 
                                          type = type, 
                                          verbose = FALSE, 
                                          return.LRT = FALSE, #FALSE
                                          warn = warn, 
                                          parallel = parallel, 
                                          ncpus = ncpus, cl = cl,
                                          double.bootstrap = "no")
            if (verbose) cat(sprintf("%5.3f", plugin.pvalue), "\n")
            attr(lrt.boot, "plugin.pvalue") <- plugin.pvalue
        } else if (double.bootstrap == "FDB") {
            #Fast double bootstrap
            plugin.pvalue <- bootstrapLRT(h0 = fit.h0, h1 = fit.h1, 
                                          R = 1, 
                                          type = type, 
                                          verbose = FALSE, 
                                          warn = warn, 
                                          return.LRT = TRUE, #TRUE
                                          parallel = parallel, 
                                          ncpus = ncpus, cl = cl,
                                          double.bootstrap = "no") 

            LRT.2 <- attr(plugin.pvalue, "LRT")

            if (verbose) cat("  ... ... LRT2 = ", LRT.2, "\n")
            attr(lrt.boot, "LRT.2") <- LRT.2
        } 
        lrt.boot
    }
  
    #Parallel processing
    RR <- sum(R)
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
        if (have_mc) {
            parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
        }
        else if (have_snow) {
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus)) #
                if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
                    parallel::clusterSetRNGStream(cl) # 
                res <- parallel::parLapply(cl, seq_len(RR), fn) # 
                parallel::stopCluster(cl) #
                res
            }
            else parallel::parLapply(cl, seq_len(RR), fn)
        }
    } else lapply(seq_len(RR), fn)

    error.idx <- integer(0)
    for (b in seq_len(RR)) {
        if (!is.null(res[[b]])) {
            LRT[b] <- res[[b]]
            if (double.bootstrap == "standard") { 
                plugin.pvalues[b] <- attr(res[[b]], "plugin.pvalue")
            } else if (double.bootstrap == "FDB") {
                LRT.2[b] <- attr(res[[b]], "LRT.2")
            }
        } else {
            error.idx <- c(error.idx, b)
        }
    }


    #Error handling
    if (length(error.idx) > 0L) {
        warning("lavaan WARNING: only ", (R - length(error.idx)), 
                " bootstrap draws were successful")
        LRT <- LRT[-error.idx]
        if(length(LRT) == 0) LRT <- as.numeric(NA)

        if (double.bootstrap == "standard") { 
            plugin.pvalues <- plugin.pvalues[-error.idx]
            attr(LRT, "error.idx") <- error.idx

        }
        if (double.bootstrap == "FDB") { 
            LRT.2 <- LRT.2[-error.idx]
            attr(LRT.2, "error.idx") <- error.idx
        }
    } else {
        if (verbose) 
            cat("Number of successful bootstrap draws:", (R - 
                                                          length(error.idx)), "\n")
    }

    pvalue <- sum(LRT > LRT.original) / length(LRT)

    if (return.LRT) { 
        attr(pvalue, "LRT.original") <- LRT.original
        attr(pvalue, "LRT") <- LRT
    }
    if (double.bootstrap == "FDB") {
        Q <- (1 - pvalue)
        lrt.q <- quantile(LRT.2, Q, na.rm = TRUE)
        adj.pvalue <- sum(LRT > lrt.q) / length(LRT)
        attr(pvalue, "lrt.q") <- lrt.q
        attr(pvalue, "adj.pvalue") <- adj.pvalue

        if (return.LRT) {
            attr(pvalue, "LRT.original") <- LRT.original
            attr(pvalue, "LRT") <- LRT
            attr(pvalue, "LRT2") <- LRT.2
        }  
    } else if (double.bootstrap == "standard") {
        adj.alpha <- quantile(plugin.pvalues, double.bootstrap.alpha,
                              na.rm=TRUE)
        attr(pvalue, "plugin.pvalues") <- plugin.pvalues
        attr(pvalue, "adj.alpha") <- adj.alpha
        # adj.pvalue <- sum(plugin.pvalues < pvalue) / length(R) #adj.pvalue < or <=??
        # attr(pvalue, "adj.pvalue") <- adj.pvalue
    }
    pvalue
}
