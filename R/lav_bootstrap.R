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
                            iseed       = NULL,
                            h0.rmsea    = NULL,
                            ...) {

    # checks
    type. <- tolower(type) # overwritten if nonparametric
    stopifnot(inherits(object, "lavaan"),
              type. %in% c("nonparametric", "ordinary",
                          "bollen.stine", "parametric", "yuan"))
    if(type. == "nonparametric") type. <- "ordinary"

    # check if options$se is not bootstrap, otherwise, we get an infinite loop
    if(object@Options$se == "bootstrap")
        stop("lavaan ERROR: se == \"bootstrap\"; please refit model with another option for \"se\"")

    # check if options$test is not bollen.stine
    if(object@Options$test == "bollen.stine")
        stop("lavaan ERROR: test == \"bollen.stine\"; please refit model with another option for \"test\"")

    # check for conditional.x = TRUE
    if(object@Model@conditional.x) {
        stop("lavaan ERROR: this function is not (yet) available if conditional.x = TRUE")
    }

    lavoptions. <- list(parallel = parallel, ncpus = ncpus, cl = cl,
                        iseed = iseed)

    bootstrap.internal(object          = object,
                       lavdata.        = NULL,
                       lavmodel.       = NULL,
                       lavsamplestats. = NULL,
                       lavoptions.     = lavoptions.,
                       lavpartable.    = NULL,
                       R               = R,
                       type            = type.,
                       verbose         = verbose,
                       FUN             = FUN,
                       warn            = warn,
                       return.boot     = return.boot,
                       h0.rmsea        = h0.rmsea,
                       ...)
}

# we need an internal version to be called from VCOV and lav_model_test
# when there is no lavaan object yet!
bootstrap.internal <- function(object          = NULL,
                               lavdata.        = NULL,
                               lavmodel.       = NULL,
                               lavsamplestats. = NULL,
                               lavoptions.     = NULL,
                               lavpartable.    = NULL,
                               R               = 1000L,
                               type            = "ordinary",
                               verbose         = FALSE,
                               FUN             = "coef",
                               warn            = 0L,
                               return.boot     = FALSE,
                               h0.rmsea        = NULL,
                               ...) {

    # warning: avoid use of 'options', 'sample' (both are used as functions
    # below...
    # options -> opt
    # sample -> samp


    # object slots
    if(!is.null(object)) {
        lavdata        <- object@Data
        lavmodel       <- object@Model
        lavsamplestats <- object@SampleStats
        lavoptions     <- object@Options
        lavpartable    <- object@ParTable
        FUN <- match.fun(FUN)
        t0 <- FUN(object, ...)
        t.star <- matrix(as.numeric(NA), R, length(t0))
        colnames(t.star) <- names(t0)
    } else {
        # internal version!
        lavdata        <- lavdata.
        lavmodel       <- lavmodel.
        lavsamplestats <- lavsamplestats.
        lavoptions     <- lavoptions.
        lavpartable    <- lavpartable.
        lavoptions$se <- "none"; lavoptions$test <- "standard"
        lavoptions$verbose <- FALSE
        if(FUN == "coef") {
            t.star <- matrix(as.numeric(NA), R, lavmodel@nx.free)
            lavoptions$test <- "none"
        } else if(FUN == "test") {
            t.star <- matrix(as.numeric(NA), R, 1L)
        } else if(FUN == "coeftest") {
            t.star <- matrix(as.numeric(NA), R, lavmodel@nx.free + 1L)
        }
    }

    parallel <- lavoptions$parallel
    ncpus    <- lavoptions$ncpus
    cl       <- lavoptions$cl
    iseed    <- lavoptions$iseed

    # prepare
    old_options <- options(); options(warn = warn)
    # the next 10 lines are borrowed from the boot package
    have_mc <- have_snow <- FALSE
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
        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
        Mu.hat <- computeMuHat(lavmodel = lavmodel)
    }

    # can we use the original data, or do we need to transform it first?
    if(type == "bollen.stine" || type == "yuan") {
        # check if data is complete
        if(lavoptions$missing != "listwise")
            stop("lavaan ERROR: bollen.stine/yuan bootstrap not available for missing data")
        dataX <- vector("list", length=lavdata@ngroups)
    } else {
        dataX <- lavdata@X
    }
    dataeXo <- lavdata@eXo
    dataWT  <- lavdata@weights
    dataMp  <- lavdata@Mp
    dataRp  <- lavdata@Rp
    dataLp  <- lavdata@Lp

    # if bollen.stine, transform data here
    if(type == "bollen.stine") {
        for(g in 1:lavsamplestats@ngroups) {
            sigma.sqrt <- lav_matrix_symmetric_sqrt(Sigma.hat[[g]])
            S.inv.sqrt <- lav_matrix_symmetric_sqrt(lavsamplestats@icov[[g]])

            # center (needed???)
            X <- scale(lavdata@X[[g]], center=TRUE, scale=FALSE)

            # transform
            X <- X %*% S.inv.sqrt %*% sigma.sqrt

            # add model-based mean
            if(lavmodel@meanstructure)
                X <- scale(X, center=(-1*Mu.hat[[g]]), scale=FALSE)

            # transformed data
            dataX[[g]] <- X
        }
    }

    # if yuan, transform data here
    if(type == "yuan") {
        # page numbers refer to Yuan et al, 2007
        # Define a function to find appropriate value of a
        # (p. 272); code supplied 16 jun 2016 by Cheng & Wu
        search.a <- function(F0, d, p) {
            if (F0 == 0) {
                a0 <- 0
                return(a0)
            }
            max.a <- 1 / (1 - min(d)) - 1e-3
            # starting value; Yuan p. 272
            a0 <- min(sqrt(2 * F0 / sum((d - 1)^2)), max.a)

            # See Yuan p. 280
            for (i in 1:50) {
                dia <- a0 * d + (1 - a0)
                g1 <- -sum(log(dia)) + sum(dia) - p
                dif <- g1 - F0
                if(abs(dif) < 1e-6) return(a0)
                g2 <- a0 * sum((d - 1)^2 / dia)
                a0 <- min(max(a0 - dif/g2, 0), max.a)
            }
            # if search fails to converge in 50 iterations
            warning("lavaan WARNING: yuan bootstrap search for `a` did not converge. h0.rmsea may be too large.")
            a0
        }

        # Now use g.a within each group
        for(g in 1:lavsamplestats@ngroups) {
            S <- lavsamplestats@cov[[g]]
            # test is in Fit slot
            ghat <- object@test[[1]]$stat.group[[g]]
            df <- object@test[[1]]$df
            Sigmahat <- Sigma.hat[[g]]
            nmv <- nrow(Sigmahat)
            n <- lavdata@nobs[[g]]

            # Calculate tauhat_1, middle p. 267.
            # Yuan et al note that tauhat_1 could be negative;
            # if so, we need to let S.a = Sigmahat. (see middle p 275)
            ifelse(length(h0.rmsea)==0,
              tau.hat <- (ghat - df)/(n-1),  # middle p 267
              tau.hat <- df*(h0.rmsea*h0.rmsea))    # middle p 273

            if (tau.hat >= 0) {
                # from Cheng and Wu
                EL <- t(chol(Sigmahat))
                ESE <- forwardsolve(EL, t(forwardsolve(EL, S)))
                d <- eigen(ESE, symmetric = TRUE, only.values = TRUE)$values
                # Find a to minimize g.a
                a <- search.a(tau.hat, d, nmv)
                # Calculate S_a (p. 267)
                S.a <- a*S + (1 - a)*Sigmahat
            } else {
                S.a <- Sigmahat
            }

            # Transform the data (p. 263)
            S.a.sqrt <- lav_matrix_symmetric_sqrt(S.a)
            S.inv.sqrt <- lav_matrix_symmetric_sqrt(lavsamplestats@icov[[g]])

            X <- lavdata@X[[g]]
            X <- X %*% S.inv.sqrt %*% S.a.sqrt

            # transformed data
            dataX[[g]] <- X
        }
    }

    # run bootstraps
    fn <- function(b) {
        if(type == "bollen.stine" || type == "ordinary" || type == "yuan") {
            # take a bootstrap sample for each group
            for(g in 1:lavsamplestats@ngroups) {
                stopifnot(lavsamplestats@nobs[[g]] > 1L)
                boot.idx <- sample(x=lavsamplestats@nobs[[g]],
                                   size=lavsamplestats@nobs[[g]], replace=TRUE)
                dataX[[g]] <- dataX[[g]][boot.idx,,drop=FALSE]
                if(!is.null(dataeXo[[g]])) {
                    dataeXo[[g]] <- dataeXo[[g]][boot.idx,,drop=FALSE]
                }
                if(!is.null(dataWT[[g]])) {
                    dataWT[[g]] <- dataWT[[g]][boot.idx]
                }
                if(lavoptions$missing != "listwise") {
                    dataMp[[g]] <- lav_data_missing_patterns(dataX[[g]],
                                       sort.freq = FALSE, coverage = FALSE)
                }
                if(length(lavdata@ov.names.x[[g]]) == 0L &&
                   all(lavdata@ov.names[[g]] %in%
                       lavdata@ov$name[lavdata@ov$type == "ordered"])) {
                    dataRp[[g]] <- lav_data_resp_patterns(dataX[[g]])
                }
                if(lavdata@nlevels > 1L) {
                    # extract cluster variable(s), for this group
                    clus <- matrix(0, nrow(dataX[[g]]), lavdata@nlevels - 1L)
                    for(l in 2:lavdata@nlevels) {
                        clus[,(l-1L)] <- lavdata@Lp[[g]]$cluster.idx[[l]]
                    }
                    dataLp[[g]] <- lav_data_cluster_patterns(Y = dataX[[g]],
                                        clus = clus,
                                        cluster = lavdata@cluster,
                                        ov.names = lavdata@ov.names[[g]],
                                        ov.names.l = lavdata@ov.names.l[[g]])
                }
            } # g
        } else { # parametric!
            for(g in 1:lavsamplestats@ngroups) {
                dataX[[g]] <- MASS::mvrnorm(n     = lavsamplestats@nobs[[g]],
                                            Sigma = Sigma.hat[[g]],
                                            mu    = Mu.hat[[g]])
            }
        }

        # create new lavdata@ object (new in 0.6-2)
        newData         <- lavdata
        newData@X       <- dataX
        newData@eXo     <- dataeXo
        newData@weights <- dataWT
        newData@Mp      <- dataMp
        newData@Rp      <- dataRp
        newData@Lp      <- dataLp

        # verbose
        if(verbose) cat("  ... bootstrap draw number:", sprintf("%4d", b))
        bootSampleStats <- try(lav_samplestats_from_data(
                               lavdata       = newData,
                               missing       = lavoptions$missing,
                               rescale       = (lavoptions$estimator == "ML" &&
                                                lavoptions$likelihood == "normal"),
                               estimator     = lavoptions$estimator,
                               mimic         = lavoptions$mimic,
                               meanstructure = lavoptions$meanstructure,
                               conditional.x = lavoptions$conditional.x,
                               group.w.free  = lavoptions$group.w.free,
                               #missing.h1    = (FUN != "coef"), # not if fixed.x, otherwise starting values fails!
                               missing.h1    = TRUE,
                               verbose       = FALSE), silent=TRUE)
        if(inherits(bootSampleStats, "try-error")) {
            if(verbose) {
                cat("     FAILED: creating sample statistics\n")
                cat(bootSampleStats[1])
            }
            options(old_options)
            return(NULL)
        }

        if(lavmodel@fixed.x && length(vnames(lavpartable, "ov.x")) > 0L) {
            model.boot <- NULL
        } else {
            model.boot <- lavmodel
        }

        # fit model on bootstrap sample
        fit.boot <- lavaan(slotOptions     = lavoptions,
                           slotParTable    = lavpartable,
                           slotModel       = model.boot,
                           slotSampleStats = bootSampleStats,
                           slotData        = lavdata)
        if(!fit.boot@optim$converged) {
            if(verbose) cat("     FAILED: no convergence\n")
            options(old_options)
            return(NULL)
        }

        # extract information we need
        if(is.null(object)) { # internal use only!
            if(FUN == "coef") {
                out <- fit.boot@optim$x
            } else if(FUN == "test") {
                out <- fit.boot@test[[1L]]$stat
            } else if(FUN == "coeftest") {
                out <- c(fit.boot@optim$x, fit.boot@test[[1L]]$stat)
            }
        } else { # general use
            out <- try(FUN(fit.boot, ...), silent=TRUE)
        }
        if(inherits(out, "try-error")) {
            if(verbose) cat("     FAILED: applying FUN to fit.boot\n")
            options(old_options)
            return(NULL)
        }
        if(verbose) cat("   OK -- niter = ",
                        sprintf("%3d", fit.boot@optim$iterations), " fx = ",
                        sprintf("%13.9f", fit.boot@optim$fx), "\n")
        out
    }

    # this is from the boot function in package boot
    RR <- R
    if(verbose) {
        cat("\n")
    }
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
        if (have_mc) {
            parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
        } else if (have_snow) {
            list(...) # evaluate any promises
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                if(RNGkind()[1L] == "L'Ecuyer-CMRG")
                    parallel::clusterSetRNGStream(cl, iseed = iseed)
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





