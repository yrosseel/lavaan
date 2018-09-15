# compare two nested models, by default using the chi-square
# difference test

# - in 0.5-16, SB.classic = TRUE is the default again (for now)
# - in 0.5-18, SB.classic is replaced by 'method', with the following
#   options:
#     method = "default" (we choose a default method, based on the estimator)
#     method = "Satorra.2000"
#     method = "Satorra.Bentler.2001"
#     method = "Satorra.Bentler.2010"
#     method = "mean.var.adjusted.PLRT"

lavTestLRT <- function(object, ..., method = "default", A.method = "delta",
                       scaled.shifted = TRUE,
                       H1 = TRUE, type = "Chisq", model.names = NULL) {

    if(object@optim$npar > 0L && !object@optim$converged)
        stop("lavaan ERROR: model did not converge")
    type <- tolower(type)
    method <- tolower( gsub("[-_\\.]", "", method ) )

    # NOTE: if we add additional arguments, it is not the same generic
    # anova() function anymore, and match.call will be screwed up

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)

    modp <- if(length(dots))
        sapply(dots, is, "lavaan") else logical(0)

    # some general properties (taken from the first model)
    estimator <- object@Options$estimator
    likelihood <- object@Options$likelihood
    ngroups <- object@Data@ngroups
    nobs <- object@SampleStats@nobs
    ntotal <- object@SampleStats@ntotal

    # shortcut for single argument (just plain LRT)
    if(!any(modp)) {
        if(type == "cf") {
            warning("lavaan WARNING: `type' argument is ignored for a single model")
        }
        aic <- bic <- c(NA, NA)
        if(estimator == "ML") {
            aic <- c(NA, AIC(object))
            bic <- c(NA, BIC(object))
        }

        val <- data.frame(Df = c(0, object@test[[1L]]$df),
                          AIC = aic,
                          BIC = bic,
                          Chisq = c(0, object@test[[1L]]$stat),
                          "Chisq diff" = c(NA, object@test[[1L]]$stat),
                          "Df diff" = c(NA, object@test[[1L]]$df),
                          "Pr(>Chisq)" = c(NA, object@test[[1L]]$pvalue),
                          row.names = c("Saturated", "Model"),
                          check.names = FALSE)
        attr(val, "heading") <- "Chi Square Test Statistic (unscaled)\n"
        class(val) <- c("anova", class(val))
        return(val)
    }

    # list of models
    mods <- c(list(object), dots[modp])
    if(!is.null(model.names)) {
        names(mods) <- model.names
    } else {
        names(mods) <- sapply(as.list(mcall)[which(c(FALSE, TRUE, modp))],
                              deparse)
    }

    ## put them in order (using number of free parameters)
    #nfreepar <- sapply(mods, function(x) x@optim$npar)
    #if(any(duplicated(nfreepar))) { ## FIXME: what to do here?
    #    # what, same number of free parameters?
    #    # maybe, we need to count number of constraints
    #    ncon <- sapply(mods, function(x) { nrow(x@Model@con.jac) })
    #    nfreepar <- nfreepar - ncon
    #}

    # put them in order (using degrees of freedom)
    ndf <- sapply(mods, function(x) x@test[[1]]$df)
    mods <- mods[order(ndf)]

    # here come the checks
    if(TRUE) {
        # 1. same set of observed variables?
        ov.names <- lapply(mods, function(x) { sort(lavNames(x)) })
        OV <- ov.names[[1L]] # the observed variable names of the first model
        if(!all(sapply(ov.names, function(x) identical(x, OV)))) {
            warning("lavaan WARNING: some models are based on a different set of observed variables")
        }
        ## wow FIXME: we may need to reorder the rows/columns first!!
        #COVS <- lapply(mods, function(x) slot(slot(x, "Sample"), "cov")[[1]])
        #if(!all(sapply(COVS, all.equal, COVS[[1]]))) {
        #    stop("lavaan ERROR: models must be fit to the same data")
        #}
        # 2. nested models? *different* npars?

        # TODO!

        # 3. all meanstructure?
        mean.structure <- sapply(mods, inspect, "meanstructure")
        if(sum(mean.structure) > 0L &&
           sum(mean.structure) < length(mean.structure)) {
            warning("lavaan WARNING: not all models have a meanstructure")
        }
    }

    mods.scaled <- unlist( lapply(mods, function(x) {
        any(c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus",
              "mean.var.adjusted", "scaled.shifted") %in%
            unlist(sapply(slot(x, "test"), "[", "test")) ) }))

    if(all(mods.scaled)) {
        scaled <- TRUE
        # which type?
        TEST <- object@test[[2]]$test
    } else if(!any(mods.scaled)) { # thanks to R.M. Bee to fix this
        scaled <- FALSE
        TEST <- "standard"
    } else {
        stop("lavaan ERROR: some models (but not all) have scaled test statistics")
    }

    # which models have used a MEANSTRUCTURE?
    mods.meanstructure <- sapply(mods, function(x) {
                                 unlist(slot(slot(x, "Model"),
                                                     "meanstructure"))})
    if(all(mods.meanstructure)) {
        meanstructure <- "ok"
    } else if(sum(mods.meanstructure) == 0) {
        meanstructure <- "ok"
    } else {
        stop("lavaan ERROR: some models (but not all) have a meanstructure")
    }

    # collect statistics for each model
    if(type == "chisq") {
        Df <- sapply(mods, function(x) slot(x, "test")[[1]]$df)
        STAT <- sapply(mods, function(x) slot(x, "test")[[1]]$stat)
    } else if(type == "cf") {
        tmp <- lapply(mods, lavTablesFitCf)
        STAT <- unlist(tmp)
        Df  <- unlist(lapply(tmp, attr, "DF"))
    } else {
        stop("lavaan ERROR: test type unknown: ", type)
    }


    # difference statistics
    STAT.delta  <- c(NA, diff(STAT))
    Df.delta     <- c(NA, diff(Df))

    # correction for scaled test statistics
    if(type == "chisq" && scaled) {

        # select method
        if(method == "default") {
            if(estimator == "PML") {
                method <- "mean.var.adjusted.PLRT"
            } else if(TEST %in% c("satorra.bentler", "yuan.bentler",
                                  "yuan.bentler.mplus")) {
                method <- "satorra.bentler.2001"
            } else {
                method <- "satorra.2000"
            }
        } else if(method == "meanvaradjustedplrt" ||
                  method == "mean.var.adjusted.PLRT") {
            method <- "mean.var.adjusted.PLRT"
            stopifnot(estimator == "PML")
        } else if(method == "satorra2000") {
            method <- "satorra.2000"
        } else if(method == "satorrabentler2001") {
            method <- "satorra.bentler.2001"
        } else if(method == "satorrabentler2010") {
            method <- "satorra.bentler.2010"
        } else {
            stop("lavaan ERROR: unknown method for scaled difference test: ", method)
        }

        if(method == "satorra.bentler.2001") {
            # use formula from Satorra & Bentler 2001
            for(m in seq_len(length(mods) - 1L)) {
                out <- lav_test_diff_SatorraBentler2001(mods[[m]], mods[[m+1]])
                STAT.delta[m+1] <- out$T.delta
                  Df.delta[m+1] <- out$df.delta
            }
        } else if (method == "mean.var.adjusted.PLRT") {
            for(m in seq_len(length(mods) - 1L)) {
                out <- ctr_pml_plrt_nested(mods[[m]], mods[[m+1]])
                STAT.delta[m+1] <- out$FSMA.PLRT
                  Df.delta[m+1] <- out$adj.df
            }
        } else if(method == "satorra.bentler.2010") {
            for(m in seq_len(length(mods) - 1L)) {
                out <- lav_test_diff_SatorraBentler2010(mods[[m]], mods[[m+1]])
                STAT.delta[m+1] <- out$T.delta
                  Df.delta[m+1] <- out$df.delta
            }
        } else if(method == "satorra.2000") {
            for(m in seq_len(length(mods) - 1L)) {
                if(TEST %in% c("satorra.bentler", "yuan.bentler",
                               "yuan.bentler.mplus")) {
                    Satterthwaite <- FALSE
                } else {
                    Satterthwaite <- TRUE
                }
                out <- lav_test_diff_Satorra2000(mods[[m]], mods[[m+1]],
                                                 H1 = TRUE,
                                                 Satterthwaite = Satterthwaite,
                                                 scaled.shifted = scaled.shifted,
                                                 A.method = A.method)
                STAT.delta[m+1] <- out$T.delta
                  Df.delta[m+1] <- out$df.delta
            }
        }
    }

    # Pvalue
    Pvalue.delta <- pchisq(STAT.delta, Df.delta, lower.tail = FALSE)

    aic <- bic <- rep(NA, length(mods))
    if(estimator == "ML") {
        aic <- sapply(mods, FUN=AIC)
        bic <- sapply(mods, FUN=BIC)
    } else if(estimator == "PML") {
        OUT <- lapply(mods, ctr_pml_aic_bic)
        aic <- sapply(OUT, "[[", "PL_AIC")
        bic <- sapply(OUT, "[[", "PL_BIC")
    }

    if(estimator == "PML") {
        val <- data.frame(Df = Df,
                          PL_AIC = aic,
                          PL_BIC = bic,
                          Chisq = STAT,
                          "Chisq diff" = STAT.delta,
                          "Df diff" = Df.delta,
                          "Pr(>Chisq)" = Pvalue.delta,
                          row.names = names(mods),
                          check.names = FALSE)
    } else {
        val <- data.frame(Df = Df,
                          AIC = aic,
                          BIC = bic,
                          Chisq = STAT,
                          "Chisq diff" = STAT.delta,
                          "Df diff" = Df.delta,
                          "Pr(>Chisq)" = Pvalue.delta,
                          row.names = names(mods),
                          check.names = FALSE)
    }

    # catch Df.delta == 0 cases (reported by Florian Zsok in Zurich)
    # but only if there are no inequality constraints! (0.6-1)
    idx <- which(val[,"Df diff"] == 0)
    if(length(idx) > 0L) {
        # remove models with inequality constraints
        ineq.idx <- which(sapply(lapply(mods, function(x) slot(slot(x, "Model"), "x.cin.idx")), length) > 0L)
        rm.idx <- which(idx %in% ineq.idx)
        if(length(rm.idx) > 0L) {
            idx <- idx[-rm.idx]
        }
    }
    if(length(idx) > 0L) {
        val[idx, "Pr(>Chisq)"] <- as.numeric(NA)
        warning("lavaan WARNING: some models have the same degrees of freedom")
    }

    if(type == "chisq") {

        if(scaled) {
            attr(val, "heading") <-
                paste("Scaled Chi Square Difference Test (method = \"",
                      method, "\")\n", sep="")
        } else {
            attr(val, "heading") <- "Chi Square Difference Test\n"
        }
    } else if(type == "cf") {
        colnames(val)[c(3,4)] <- c("Cf", "Cf diff")
        attr(val, "heading") <- "Cf Difference Test\n"
    }
    class(val) <- c("anova", class(val))

    return(val)

}

