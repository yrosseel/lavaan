# main function used by various bootstrap related functions
# this function draws the bootstrap samples, and estimates the
# free parameters for each bootstrap sample
#
# return COEF matrix of size R x npar (R = number of bootstrap samples)
# 
# Ed. 9 mar 2012
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



bootstrapLavaan <- function(object, 
                            R           = 1000L, 
                            type        = "ordinary",
                            verbose     = FALSE, 
                            FUN         = "coef",
                            warn        = -1L,
                            return.boot = FALSE,
                            parallel    = c("no", "multicore", "snow"),
                            ncpus       = 1L,
                            cl          = NULL,
                            ...) {

    # checks
    type. <- type # overwritten if nonparametric
    stopifnot(class(object) == "lavaan",
              type. %in% c("nonparametric", "ordinary",
                          "bollen.stine", "parametric", "yuan"))
    if(type. == "nonparametric") type. <- "ordinary"

    # check if options$se is not bootstrap, otherwise, we get an infinite loop
    if(object@Options$se == "bootstrap")
        stop("lavaan ERROR: se == \"bootstrap\"; please refit model with another option for \"se\"")
     
    # check if options$test is not bollen.stine
    if(object@Options$test == "bollen.stine")
        stop("lavaan ERROR: test == \"bollen.stine\"; please refit model with another option for \"test\"")

    bootstrap.internal(object = object,
                       data.       = NULL,
                       model.      = NULL,
                       sample.     = NULL,
                       options.    = NULL,
                       user.       = NULL,
                       R           = R,
                       type        = type.,
                       verbose     = verbose,
                       FUN         = FUN,
                       warn        = warn,
                       return.boot = return.boot,
                       parallel    = parallel,
                       ncpus       = ncpus,
                       cl          = cl,
                       ...)
}

# we need an internal version to be called from VCOV and computeTestStatistic
# when there is no lavaan object yet!
bootstrap.internal <- function(object = NULL,
                               data.       = NULL,
                               model.      = NULL,
                               sample.     = NULL,
                               options.    = NULL,
                               user.       = NULL,
                               R           = 1000L,
                               type        = "ordinary",
                               verbose     = FALSE,
                               FUN         = "coef",
                               warn        = 0L,
                               return.boot = FALSE,
                               parallel    = c("no", "multicore", "snow"),
                               ncpus       = 1L,
                               cl          = NULL,
                               ...) {

    # warning: avoid use of 'options', 'sample' (both are used as functions
    # below...
    # options -> opt
    # sample -> samp


    # object slots
    if(!is.null(object)) {
        data <- object@Data; model <- object@Model; samp <- object@Sample
        opt <- object@Options; user <- object@User
        FUN <- match.fun(FUN)
        t0 <- FUN(object, ...)
        t.star <- matrix(as.numeric(NA), R, length(t0))
        colnames(t.star) <- names(t0)
    } else {
        # internal version!
        data <- data.; model <- model.; samp <- sample.
        opt <- options.; user <- user.
        opt$se <- "none"; opt$test <- "standard"
        opt$verbose <- FALSE
        if(FUN == "coef") {
            t.star <- matrix(as.numeric(NA), R, model@nx.free)
            opt$test <- "none"
        } else if(FUN == "test") {
            t.star <- matrix(as.numeric(NA), R, 1L)
        } else if(FUN == "coeftest") {
            t.star <- matrix(as.numeric(NA), R, model@nx.free + 1L)
        }
    }

    # prepare
    old_options <- options(); options(warn = warn)
    if (missing(parallel)) parallel <- "no"
    # the next 10 lines are borrowed from the boot package
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    ncpus <- ncpus
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!have_mc && !have_snow) ncpus <- 1L
    }

    # only if we return the seed
    #if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    #seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    # bollen.stine, yuan, or parametric: we need the Sigma.hat values
    if(type == "bollen.stine" || type == "parametric" || type == "yuan") {
        Sigma.hat <- computeSigmaHat(model)
        Mu.hat <- computeMuHat(model)
    }

    
    # can we use the original data, or do we need to transform it first?
    if(type == "bollen.stine" || type == "yuan") {
        # check if data is complete
        if(opt$missing != "listwise")
            stop("lavaan ERROR: bollen.stine/yuan bootstrap not available for missing data")
        dataX <- vector("list", length=data@ngroups)
    } else {
        dataX <- data@X
    }

    # if bollen.stine, transform data here
    if(type == "bollen.stine") {
        for(g in 1:samp@ngroups) {
            sigma.sqrt <- sqrtSymmetricMatrix(Sigma.hat[[g]])
            S.inv.sqrt <- sqrtSymmetricMatrix(samp@icov[[g]])

            # center (needed???)
            X <- scale(data@X[[g]], center=TRUE, scale=FALSE)

            # transform
            X <- X %*% S.inv.sqrt %*% sigma.sqrt

            # add model-based mean
            if(model@meanstructure)
                X <- scale(X, center=(-1*Mu.hat[[g]]), scale=FALSE)

            # transformed data
            dataX[[g]] <- X
        }
    }

    # if yuan, transform data here
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
        for(g in 1:samp@ngroups) {
            S <- samp@cov[[g]]
            # test is in Fit slot
            ghat <- object@Fit@test[[1]]$stat.group[[g]]
            df <- object@Fit@test[[1]]$df
            Sigmahat <- Sigma.hat[[g]]
            Sigmahat.inv <- samp@icov[[g]]
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
            S.inv.sqrt <- sqrtSymmetricMatrix(samp@icov[[g]])

            X <- data@X[[g]]
            X <- X %*% S.inv.sqrt %*% S.a.sqrt            

            # transformed data
            dataX[[g]] <- X
        }
    }

    # run bootstraps
    fn <- function(b) {
        if(type == "bollen.stine" || type == "ordinary" || type == "yuan") {
            # take a bootstrap sample for each group
            for(g in 1:samp@ngroups) {
                stopifnot(samp@nobs[[g]] > 1L)
                boot.idx <- sample(x=samp@nobs[[g]],
                                   size=samp@nobs[[g]], replace=TRUE)
                dataX[[g]] <- dataX[[g]][boot.idx,,drop=FALSE]
            }
        } else { # parametric!
            for(g in 1:samp@ngroups) {
                dataX[[g]] <- MASS.mvrnorm(n     = samp@nobs[[g]],
                                           Sigma = Sigma.hat[[g]],
                                           mu    = Mu.hat[[g]])
            }
        }

        # verbose
        if(verbose) cat("  ... bootstrap draw number:", sprintf("%4d", b))

        bootSampleStats <- try(getSampleStatsFromData(
                               Data          = NULL,
                               DataX         = dataX,
                               missing       = opt$missing,
                               rescale       = (opt$estimator == "ML" &&
                                                opt$likelihood == "normal"),
                               WLS.V         = NULL,
                               estimator     = opt$estimator,
                               mimic         = opt$mimic,
                               meanstructure = opt$meanstructure,
                               missing.h1    = (FUN != "coef"),
                               verbose       = FALSE)) 
        if(inherits(bootSampleStats, "try-error")) {
            if(verbose) cat("     FAILED: creating sample statistics\n")
            options(old_options)
            return(NULL)
        }

        # just in case we need the new X in the data slot (lm!)
        data@X <- dataX

        # adjust model slot if fixed.x variances/covariances
        # have changed:
      ### FIXME #####
        #if(model@fixed.x && length(vnames(user, "ov.x")) > 0L) {
        #    for(g in 1:samp@ngroups) {
        #        
        #    }
        #}
        if(model@fixed.x && length(vnames(user, "ov.x")) > 0L) {
            model.boot <- NULL
        } else {
            model.boot <- model
        }

        # fit model on bootstrap sample
        fit.boot <- lavaan(slotOptions = opt,
                           slotUser    = user,
                           slotModel   = model.boot,
                           slotSample  = bootSampleStats,
                           slotData    = data)
        if(!fit.boot@Fit@converged) {
            if(verbose) cat("     FAILED: no convergence\n")
            options(old_options)
            return(NULL)
        } 
        
        # extract information we need
        if(is.null(object)) { # internal use only!
            if(FUN == "coef") {
                out <- fit.boot@Fit@x
            } else if(FUN == "test") { 
                out <- fit.boot@Fit@test[[1L]]$stat
            } else if(FUN == "coeftest") { 
                out <- c(fit.boot@Fit@x, fit.boot@Fit@test[[1L]]$stat)
            } 
        } else { # general use
            out <- try(FUN(fit.boot, ...))
        }
        if(inherits(out, "try-error")) {
            if(verbose) cat("     FAILED: applying FUN to fit.boot\n")
            options(old_options)
            return(NULL)
        } 
        if(verbose) cat("   OK -- niter = ", 
                        sprintf("%3d", fit.boot@Fit@iterations), " fx = ",
                        sprintf("%13.9f", fit.boot@Fit@fx), "\n")
        out
    }

    # this is from the boot function in package boot
    RR <- R
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
        if (have_mc) {
            parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
        } else if (have_snow) {
            list(...) # evaluate any promises
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
        if(!is.null(res[[b]]) && length(res[[b]]) > 0L) {
            t.star[b, ] <- res[[b]]
        } else {
            error.idx <- c(error.idx, b)
        }
    }


    # handle errors
    if(length(error.idx) > 0L) {
        warning("lavaan WARNING: only ", (R-length(error.idx)), " bootstrap draws were successful")
        t.star <- t.star[-error.idx,,drop=FALSE]
        attr(t.star, "error.idx") <- error.idx
    } else {
        if(verbose) cat("Number of successful bootstrap draws:",
                        (R - length(error.idx)), "\n")
    }

    # NOT DONE YET
    if(return.boot) {
        # mimic output boot function
    } 

    # restore options
    options(old_options)

    t.star
}





