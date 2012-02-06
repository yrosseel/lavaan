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
    stopifnot(class(object) == "lavaan",
              type %in% c("nonparametric", "ordinary",
                          "bollen.stine", "parametric"))
    if(type == "nonparametric") type <- "ordinary"

    # check if options$se is not bootstrap, otherwise, we get an infinite loop
    if(object@Options$se == "bootstrap") {
        warning("lavaan WARNING: object@Options$se == \"bootstrap\"; will set this to \"none\"")
        object@Options$se <- "none"
    }
    # check if options$test is not bollen.stine
    if(object@Options$test == "bollen.stine") {
        warning("lavaan WARNING: object@Options$test == \"bollen.stine\"; will set this to \"standard\"")
        object@Options$test <- "standard"
    }

    bootstrap.internal(object = object,
                       data        = NULL,
                       model       = NULL,
                       sample      = NULL,
                       options     = NULL,
                       user        = NULL,
                       R           = R,
                       type        = type,
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
                               data        = NULL,
                               model       = NULL,
                               sample      = NULL,
                               options     = NULL,
                               user        = NULL,
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


    # object slots
    if(!is.null(object)) {
        data <- object@Data; model <- object@Model; sample <- object@Sample
        options <- object@Options; user <- object@User
        FUN <- match.fun(FUN)
        t0 <- FUN(object, ...)
        t.star <- matrix(as.numeric(NA), R, length(t0))
        colnames(t.star) <- names(t0)
    } else {
        # internal version!
        options$se <- "none"; options$test <- "standard"
        options$verbose <- FALSE
        if(FUN == "coef") {
            t.star <- matrix(as.numeric(NA), R, model@nx.free)
            options$test <- "none"
        } else if(FUN == "test") {
            t.star <- matrix(as.numeric(NA), R, 1L)
        } else if(FUN == "coeftest") {
            t.star <- matrix(as.numeric(NA), R, model@nx.free + 1L)
        }
    }

    # prepare
    options(warn = warn)
    if (missing(parallel)) parallel <- "no"
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!have_mc && !have_snow) ncpus <- 1L
    }

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    # bollen.stine or parametric: we need the Sigma.hat values
    if(type == "bollen.stine" || type == "parametric") {
        Sigma.hat <- computeSigmaHat(model)
        Mu.hat <- computeMuHat(model)
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
    fn <- function(b) {
        if(type == "bollen.stine" || type == "ordinary") {
            # take a bootstrap sample for each group
            boot.idx <- vector("list", length=sample@ngroups)
            for(g in 1:sample@ngroups) {
                stopifnot(h0@Sample@nobs[[g]] > 1L)
                boot.idx[[g]] <- sample(x=sample@nobs[[g]],
                                        size=sample@nobs[[g]], replace=TRUE)
            }
        } else { # parametric!
            boot.idx <- NULL
            for(g in 1:sample@ngroups) {
                    data[[g]] <- MASS.mvrnorm(n     = sample@nobs[[g]],
                                              Sigma = Sigma.hat[[g]],
                                              mu    = Mu.hat[[g]])
            }
        }
        # names
        for(g in 1:h0@Sample@ngroups) 
            colnames(data[[g]]) <- sample@ov.names[[g]]

        # verbose
        if(verbose) cat("  ... bootstrap draw number:", sprintf("%4d", b))

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
            return(NULL)
        }

        # just in case we need the dataSlot (lm!)
        if(type == "bollen.stine") {
            data.boot <- data
            for(g in 1:sample@ngroups) {
                data.boot[[g]] <- data[[g]][ boot.idx[[g]],,drop=FALSE]
            }
        } else {
            data.boot <- data
        }

        # adjust model slot if fixed.x variances/covariances
        # have changed:
      ### FIXME #####
        #if(model@fixed.x && length(vnames(user, "ov.x")) > 0L) {
        #    for(g in 1:sample@ngroups) {
        #        
        #    }
        #}
        if(model@fixed.x && length(vnames(user, "ov.x")) > 0L) {
            model.boot <- NULL
        } else {
            model.boot <- model
        }

        # fit model on bootstrap sample
        fit.boot <- lavaan(slotOptions = options,
                           slotUser    = user,
                           slotModel   = model.boot,
                           slotSample  = bootSampleStats,
                           slotData    = data.boot)
        if(!fit.boot@Fit@converged) {
            if(verbose) cat("     FAILED: no convergence\n")
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

    t.star
}





