bootstrapLRT <- function (h0 = NULL, h1 = NULL, R = 1000L,
                          type = "bollen.stine", verbose = FALSE,
                          return.LRT = FALSE, 
                          double.bootstrap = "no",
                          double.bootstrap.R = 500L, 
                          double.bootstrap.alpha = 0.05, 
                          warn = -1L, parallel = c("no", "multicore", "snow"), 
                          ncpus = 1L, cl = NULL) 
{
    # checks
    type <- tolower(type)
    stopifnot(class(h0) == "lavaan", 
              class(h1) == "lavaan", 
              type %in% c("bollen.stine", "parametric", "yuan"), 
              double.bootstrap %in% c("no", "FDB", "standard"))
  
    old_options <- options(); options(warn = warn)

    # prepare
    LRT <- rep(as.numeric(NA), R)
    if((h1@Fit@fx - h0@Fit@fx) > (.Machine$double.eps * 10)) { 
        # restricted fit should not be better!
        cat(" ... h0@Fit@fx = ", h0@Fit@fx, "h1@Fit@fx = ", h1@Fit@fx,
            "h0 should not be better!\n")
        options(old_options)
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

    # can we use the original data, or do we need to transform it first?
    if(type == "bollen.stine" || type == "yuan") {
        # check if data is complete
        if(h0@Options$missing != "listwise")
            stop("lavaan ERROR: bollen.stine/yuan bootstrap not available for missing data")
        dataX <- vector("list", length=data@ngroups)
    } else {
        dataX <- data@X
    }
  
    #Bollen-Stine data transformation
    if(type == "bollen.stine") {
        for(g in 1:h0@Data@ngroups) {
            sigma.sqrt <- sqrtSymmetricMatrix(     Sigma.hat[[g]])
            S.inv.sqrt <- sqrtSymmetricMatrix(h0@SampleStats@icov[[g]])

            # center
            X <- scale(data@X[[g]], center = TRUE, scale = FALSE)

            # transform
            X <- X %*% S.inv.sqrt %*% sigma.sqrt

            # add model based mean
            if (h0@Model@meanstructure) 
                X <- scale(X, center = (-1 * Mu.hat[[g]]), scale = FALSE)

            # transformed data
            dataX[[g]] <- X
        }
    }

    # Yuan et al data transformation
    if(type == "yuan") {
        # page numbers refer to Yuan et al, 2007      
        # Define a function to find appropriate value of a
        # (p. 272)
        g.a <- function(a, Sigmahat, Sigmahat.inv, S, tau.hat, p){
            S.a <- a*S + (1-a)*Sigmahat
            tmp.term <- S.a %*% Sigmahat.inv
            res <- ((sum(diag(tmp.term)) - log(det(tmp.term)) - p) - tau.hat)^2
            # From p 272
            attr(res, "gradient") <- sum(diag((S - Sigmahat) %*%
                                     (Sigmahat.inv - chol2inv(chol(S.a)))))
            res
        }
      
        # Now use g.a within each group
        for(g in 1:h0@Data@ngroups) {
            S <- h0@SampleStats@cov[[g]]
            # test is in Fit slot
            ghat <- h0@Fit@test[[1]]$stat.group[[g]]
            df <- h0@Fit@test[[1]]$df
            Sigmahat <- Sigma.hat[[g]]
            Sigmahat.inv <- inv.chol(Sigmahat)
            nmv <- nrow(Sigmahat)
            n <- data@nobs[[g]]

            # Calculate tauhat_1, middle p. 267.
            # Yuan et al note that tauhat_1 could be negative;
            # if so, we need to let S.a = Sigmahat. (see middle p 275)
            tau.hat <- (ghat - df)/(n-1)

            if (tau.hat >= 0){
              # Find a to minimize g.a
              a <- optimize(g.a, c(0,1), Sigmahat, Sigmahat.inv,
                            S, tau.hat, nmv)$minimum

              # Calculate S_a (p. 267)
              S.a <- a*S + (1-a)*Sigmahat
            } else {
              S.a <- Sigmahat
            }

            # Transform the data (p. 263)
            S.a.sqrt <- sqrtSymmetricMatrix(S.a)
            S.inv.sqrt <- sqrtSymmetricMatrix(h0@SampleStats@icov[[g]])

            X <- data@X[[g]]
            X <- X %*% S.inv.sqrt %*% S.a.sqrt            

            # transformed data
            dataX[[g]] <- X
        }
    }

    # run bootstraps
    fn <- function(b) {
        if (type == "bollen.stine" || type == "yuan") {
            # take a bootstrap sample for each group
            for(g in 1:h0@Data@ngroups) {
                stopifnot(h0@SampleStats@nobs[[g]] > 1L)
                boot.idx <- sample(x = h0@SampleStats@nobs[[g]], 
                                   size = h0@SampleStats@nobs[[g]], replace = TRUE)
                dataX[[g]] <- dataX[[g]][boot.idx,,drop=FALSE]
            }
        } else { # parametric!
            for(g in 1:h0@Data@ngroups) {
                dataX[[g]] <- MASS::mvrnorm(n     = h0@SampleStats@nobs[[g]], 
                                            mu    = Mu.hat[[g]], 
                                            Sigma = Sigma.hat[[g]])
            }
        }

        # verbose
        if (verbose) cat("  ... bootstrap draw number: ", b, "\n")

        #Get sample statistics
        bootSampleStats <- try(lavSampleStatsFromData(
                               Data     = NULL, 
                               DataX    = dataX,
                               DataOv        = data@ov,
                               DataOvnames   = data@ov.names,
                               missing       = h0@Options$missing,
                               rescale  = (h0@Options$estimator == "ML" && 
                                           h0@Options$likelihood =="normal"), 
                               estimator     = h0@Options$estimator,
                               mimic         = h0@Options$mimic,
                               meanstructure = h0@Options$meanstructure,
                               missing.h1    = TRUE,
                               verbose  = FALSE), silent=TRUE)
        if (inherits(bootSampleStats, "try-error")) {
            if (verbose) cat("     FAILED: creating h0@SampleStats statistics\n")
            options(old_options)
            return(NULL)
        }

        # just in case we need the new X in the data slot (lm!)
        data@X <- dataX

        if (verbose) cat("  ... ... model h0: ")
        h0@Options$verbose <- FALSE
        h0@Options$se <- "none"
        h0@Options$test <- "standard"

        #Fit h0 model
        fit.h0 <- lavaan(slotOptions     = h0@Options,
                         slotParTable    = h0@ParTable, 
                         slotSampleStats = bootSampleStats, 
                         slotData        = data)
        if (!fit.h0@Fit@converged) {
            if (verbose) cat("     FAILED: no convergence\n")
            options(old_options)
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
        fit.h1 <- lavaan(slotOptions     = h1@Options, 
                         slotParTable    = h1@ParTable, 
                         slotSampleStats = bootSampleStats, 
                         slotData        = data)

        if (!fit.h1@Fit@converged) {
            if (verbose) 
                cat("     FAILED: no convergence -- niter = ", fit.h1@Fit@iterations, 
                    " fx = ", fit.h1@Fit@fx,"\n")
            options(old_options)
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
            options(old_options)
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
                                          R = 1L, 
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
        attr(pvalue, "adj.alpha") <- adj.alpha
        adj.pvalue <- sum(plugin.pvalues < pvalue) / length(plugin.pvalues)
        attr(pvalue, "plugin.pvalues") <- plugin.pvalues
        attr(pvalue, "adj.pvalue") <- adj.pvalue
    }

    # restore options
    options(old_options)

    pvalue
}