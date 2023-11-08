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
#
# - 0.6-13: RMSEA.D (also known as 'RDR') is added to the table (unless scaled)
# - 0.6-13: fix multiple-group UG^2 bug in Satorra.2000 (reported by
#           Gronneberg, Foldnes and Moss)


lavTestLRT <- function(object, ..., method = "default", A.method = "delta",
                       scaled.shifted = TRUE,
                       H1 = TRUE, type = "Chisq", model.names = NULL) {

    type <- tolower(type)
    method <- tolower( gsub("[-_\\.]", "", method ) )
    if(type %in% c("browne", "browne.residual.adf", "browne.residual.nt")) {
        if(type == "browne") {
            type <- "browne.residual.adf"
        }
        if(!method %in% c("default", "standard")) {
            stop("lavaan ERROR: method cannot be used if type is browne.residual.adf or browne.residual.nt")
        }
        method <- "default"
    }

    # NOTE: if we add additional arguments, it is not the same generic
    # anova() function anymore, and match.call will be screwed up

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    modp <- if(length(dots)) {
        sapply(dots, inherits, "lavaan")
    } else {
         logical(0L)
    }

    # some general properties (taken from the first model)
    estimator <- object@Options$estimator
    likelihood <- object@Options$likelihood
    ngroups <- object@Data@ngroups
    nobs <- object@SampleStats@nobs
    ntotal <- object@SampleStats@ntotal

    # TDJ: check for user-supplied h1 model
    user_h1_exists <- FALSE
    if( !is.null(object@external$h1) ) {
      if(inherits(object@external$h1, "lavaan")) {
        user_h1_exists <- TRUE
      }
    }

    # shortcut for single argument (just plain LRT)
    if(!any(modp) && !user_h1_exists) {
        if(type == "cf") {
            warning("lavaan WARNING: `type' argument is ignored for a single model")
        }
        return(lav_test_lrt_single_model(object))
    }

    # list of models
    mods <- c(list(object), dots[modp])
    if(!is.null(model.names)) {
        names(mods) <- model.names
    } else {
        names(mods) <- sapply(as.list(mcall)[which(c(FALSE, TRUE, modp))],
                              function(x) deparse(x))
    }
    # TDJ: Add user-supplied h1 model, if it exists
    if(user_h1_exists)  mods$h1 <- object@external$h1

    # put them in order (using degrees of freedom)
    ndf <- sapply(mods, function(x) x@test[[1]]$df)
    order.idx <- order(ndf)
    mods <- mods[order.idx]
    ndf <- ndf[order.idx]

    # here come the checks -- eventually, an option may skip this
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

        # 4. all converged?
        if(!all(sapply(mods, lavInspect, "converged"))) {
            warning("lavaan WARNING: not all models converged")
        }
    }

    mods.scaled <- unlist( lapply(mods, function(x) {
        any(c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus",
              "mean.var.adjusted", "scaled.shifted") %in%
            unlist(sapply(slot(x, "test"), "[", "test")) ) }))

    if(all(mods.scaled | ndf == 0) && any(mods.scaled)) {
        # Note: if df=0, test is not really robust, hence the above condition
        scaled <- TRUE
        # which type?
        TEST <- object@test[[2]]$test
    } else if(!any(mods.scaled)) { # thanks to R.M. Bee to fix this
        scaled <- FALSE
        TEST <- "standard"
    } else {
        stop("lavaan ERROR: some models (but not all) have scaled test statistics")
    }
    if(type %in% c("browne.residual.adf", "browne.residual.nt")) {
        scaled <- FALSE
    }

    # select method
    if(method == "default") {
        if(estimator == "PML") {
            method <- "mean.var.adjusted.PLRT"
        } else if(scaled) {
            if(TEST %in% c("satorra.bentler", "yuan.bentler",
                           "yuan.bentler.mplus")) {
                method <- "satorra.bentler.2001"
            } else {
                method <- "satorra.2000"
            }
        } else {
            # nothing to do
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

    # check method if scaled = FALSE
    if(type == "chisq" && !scaled &&
        method %in% c("mean.var.adjusted.PLRT",
                      "satorra.bentler.2001",
                      "satorra.2000",
                      "satorra.bentler.2010")) {

            warning("lavaan WARNING: method = ", dQuote(method),
                    "\n\t but no robust test statistics were used;",
                    "\n\t switching to the standard chi-square difference test")

        method <- "default"
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
    } else if(type == "browne.residual.nt") {
        TESTlist <- lapply(mods,
                           function(x) lavTest(x, test = "browne.residual.nt"))
        Df <- sapply(TESTlist, function(x) x$df)
        STAT <- sapply(TESTlist, function(x) x$stat)
    } else if(type == "browne.residual.adf") {
        TESTlist <- lapply(mods,
                           function(x) lavTest(x, test = "browne.residual.adf"))
        Df <- sapply(TESTlist, function(x) x$df)
        STAT <- sapply(TESTlist, function(x) x$stat)
    } else if(type == "cf") {
        tmp <- lapply(mods, lavTablesFitCf)
        STAT <- unlist(tmp)
        Df  <- unlist(lapply(tmp, attr, "DF"))
    } else {
        stop("lavaan ERROR: test type unknown: ", type)
    }


    # difference statistics
    STAT.delta   <- c(NA, diff(STAT))
    Df.delta     <- c(NA, diff(Df))
    if(method == "satorra.2000" && scaled.shifted) {
        a.delta <- b.delta <- rep(as.numeric(NA), length(STAT))
    }
    # new in 0.6-13
    if(!scaled) {
        RMSEA.delta  <- c(NA, lav_fit_rmsea(X2 = STAT.delta[-1],
                                            df = Df.delta[-1],
                                            N = ntotal, G = ngroups))
    }

    # check for negative values in STAT.delta
    # but with a tolerance (0.6-12)!
    if(any(STAT.delta[-1] < -1*.Machine$double.eps^(1/3))) {
        txt <- c("Some restricted models fit better than less ",
                 "restricted models; either these models are not nested, or ",
                 "the less restricted model failed to reach a global optimum.")
        txt <- c(txt, " Smallest difference = ", min(STAT.delta[-1]))
        warning(lav_txt2message(txt))
    }

    # correction for scaled test statistics
    if(type == "chisq" && scaled) {

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
                out <- lav_test_diff_SatorraBentler2010(mods[[m]], mods[[m+1]],
                                                        H1 = FALSE) # must be F

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
                if(scaled.shifted) {
                    a.delta[m+1] <- out$a
                    b.delta[m+1] <- out$b
                }
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
        if(scaled) {
            val <- data.frame(Df = Df,
                              AIC = aic,
                              BIC = bic,
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
                              "RMSEA" = RMSEA.delta,
                              "Df diff" = Df.delta,
                              "Pr(>Chisq)" = Pvalue.delta,
                              row.names = names(mods),
                              check.names = FALSE)
        }
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
            txt <- paste("The ", dQuote("Chisq"), " column contains standard ",
                      "test statistics, not the robust test that should be ",
                      "reported per model. A robust difference test is a ",
                      "function of two standard (not robust) statistics.",
                      sep = "")
            attr(val, "heading") <-
                paste("\nScaled Chi-Squared Difference Test (method = ",
                      dQuote(method), ")\n\n",
                      lav_txt2message(txt, header = "lavaan NOTE:",
                                           footer = " "), sep = "")
            if(method == "satorra.2000" && scaled.shifted) {
                attr(val, "scale") <- a.delta
                attr(val, "shift") <- b.delta
            }
        } else {
            attr(val, "heading") <- "\nChi-Squared Difference Test\n"
        }
    } else if(type  == "browne.residual.adf") {
        attr(val, "heading") <- "\nChi-Squared Difference Test based on Browne's residual (ADF) Test\n"
    } else if(type == "browne.residual.nt") {
        attr(val, "heading") <- "\nChi-Squared Difference Test based on Browne's residual (NT) Test\n"
    } else if(type == "cf") {
        colnames(val)[c(3,4)] <- c("Cf", "Cf diff")
        attr(val, "heading") <- "\nCf Difference Test\n"
    }
    class(val) <- c("anova", class(val))

    return(val)

}


# anova table for a single model
lav_test_lrt_single_model <- function(object) {

    estimator <- object@Options$estimator

    aic <- bic <- c(NA, NA)
    if(estimator == "ML") {
        aic <- c(NA, AIC(object))
        bic <- c(NA, BIC(object))
    }

    if(length(object@test) > 1L) {
        val <- data.frame(Df = c(0, object@test[[2L]]$df),
                          AIC = aic,
                          BIC = bic,
                          Chisq = c(0, object@test[[2L]]$stat),
                          "Chisq diff" = c(NA, object@test[[2L]]$stat),
                          "Df diff" = c(NA, object@test[[2L]]$df),
                          "Pr(>Chisq)" = c(NA, object@test[[2L]]$pvalue),
                          row.names = c("Saturated", "Model"),
                          check.names = FALSE)
        attr(val, "heading") <- "Chi-Squared Test Statistic (scaled)\n"
    } else {
        val <- data.frame(Df = c(0, object@test[[1L]]$df),
                          AIC = aic,
                          BIC = bic,
                          Chisq = c(0, object@test[[1L]]$stat),
                          "Chisq diff" = c(NA, object@test[[1L]]$stat),
                          "Df diff" = c(NA, object@test[[1L]]$df),
                          "Pr(>Chisq)" = c(NA, object@test[[1L]]$pvalue),
                          row.names = c("Saturated", "Model"),
                          check.names = FALSE)
        attr(val, "heading") <- "Chi-Squared Test Statistic (unscaled)\n"
    }

    class(val) <- c("anova", class(val))

    val
}

