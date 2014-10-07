# compare to nested models, by default using the chi-square
# difference test

# NOTE: in 0.5-16, SB.classic = TRUE is the default again (for now)

lavTestLRT <- function(object, ..., SB.classic = TRUE, SB.H0 = FALSE,
                       type = "Chisq", model.names = NULL) {

    if(object@Fit@npar > 0L && !object@Fit@converged)
        stop("lavaan ERROR: model did not converge")
    type <- tolower(type)

    # NOTE: if we add additional arguments, it is not the same generic
    # anova() function anymore, and match.call will be screwed up

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    #arg.names <- names(dots)
    #arg.idx <- which(nchar(arg.names) > 0L)
    #if(length(arg.idx) > 0L) {
    #    if(!is.null(dots$SB.classic))
    #        SB.classic <- dots$SB.classic
    #    if(!is.null(dots$SB.H0))
    #        SB.H0 <- dots$SB.H0           
    #    dots <- dots[-arg.idx]
    #}
  
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

        val <- data.frame(Df = c(0, object@Fit@test[[1L]]$df),
                          AIC = aic,
                          BIC = bic,
                          Chisq = c(0, object@Fit@test[[1L]]$stat),
                          "Chisq diff" = c(NA, object@Fit@test[[1L]]$stat),
                          "Df diff" = c(NA, object@Fit@test[[1L]]$df),
                          "Pr(>Chisq)" = c(NA, object@Fit@test[[1L]]$pvalue),
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
        names(mods) <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], deparse)
    }

    ## put them in order (using number of free parameters)
    #nfreepar <- sapply(mods, function(x) x@Fit@npar)
    #if(any(duplicated(nfreepar))) { ## FIXME: what to do here?
    #    # what, same number of free parameters?
    #    # maybe, we need to count number of constraints
    #    ncon <- sapply(mods, function(x) { nrow(x@Model@con.jac) })
    #    nfreepar <- nfreepar - ncon
    #}

    # put them in order (using degrees of freedom)
    ndf <- sapply(mods, function(x) x@Fit@test[[1]]$df)    
    mods <- mods[order(ndf)]

    # here come the checks
    if(TRUE) {
        # 1. same set of observed variables?
        ov.names <- lapply(mods, lavNames)
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
        if(!all(mean.structure)) {
            warning("lavaan WARNING: not all models have a meanstructure")
        }
    }

    mods.scaled <- unlist( lapply(mods, function(x) {
        any(c("satorra.bentler", "yuan.bentler", 
              "mean.var.adjusted", "scaled.shifted") %in% 
            unlist(sapply(slot(slot(x, "Fit"), "test"), "[", "test")) ) }))

    if(all(mods.scaled)) {
        scaled <- TRUE
        # which type?
        TEST <- object@Fit@test[[2]]$test
    } else if(!all(mods.scaled)) {
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
        Df <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"), 
                            "test")[[1]]$df))
    } else if(type == "cf") {
        Df <- rep(as.numeric(NA), length(mods))
    } else {
        stop("lavaan ERROR: test type unknown: ", type)
    }
    

    if(type == "chisq") {
        STAT <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"), 
                            "test")[[1]]$stat))
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
        if(SB.classic && TEST %in% c("satorra.bentler", "yuan.bentler")) {
            # use formula from Satorra & Bentler 2001
            scaling.factor <- unlist(lapply(mods, 
                function(x) slot(slot(x, "Fit"), "test")[[2]]$scaling.factor))
 
            # saturated model? (has NA scaling.factor)
            sat.idx <- which(Df == 0)
            if(length(sat.idx) > 0L) {
                scaling.factor[sat.idx] <- 1
            }

            cd1 <- diff(scaling.factor * Df)/diff(Df)


            # check for negative scaling factors
            if(any(cd1 < 0)) {
                warning("lavaan WARNING: some scaling factors are negative: [",
                        paste(round(cd1, 3), collapse=" "),"]; rerun with SB.classic=FALSE")
                cd1[cd1 < 0] <- NA
            }
            cd <- c(NA, cd1)
            STAT.delta <- STAT.delta/cd

            # extract scaled Chisq for each model
            STAT <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"),
                            "test")[[2]]$stat))
        } else {
            # see Mplus Web Note 10 (2006)
            for(m in seq_len(length(mods) - 1L)) {

                if(mods[[m]]@Fit@test[[1]]$df == mods[[m+1]]@Fit@test[[1]]$df) {
                    warnings("lavaan WARNING: some models have the same number of free parameters")
                    next
                }

                if(SB.H0) {
                    # evaluate under H0
                    stop("SB.H0 has not been implemented yet. FIXME!")
                }

                # original M (Satorra)
                Delta1 <- computeDelta(lavmodel = mods[[m]]@Model)
                npar <- ncol(Delta1[[1]])
                WLS.V <- lav_object_inspect_wls_v( mods[[m]] )  ## always H1
                Gamma <- lav_object_inspect_sampstat_nacov( mods[[m]] ) 
                

                # weight WLS.V
                for(g in 1:ngroups) {
                    WLS.V[[g]] <- nobs[[g]]/ntotal * WLS.V[[g]]
                }

                # information matrix
                P1 <- matrix(0, nrow=npar, ncol=npar)
                for(g in 1:ngroups) {
                     P1 <- P1 + (t(Delta1[[g]]) %*% WLS.V[[g]] %*% Delta1[[g]])
                }
                P1.inv <- solve(P1)

                # compute A for these two nested models
                p1 <- mods[[m   ]]@ParTable # partable h1
                p0 <- mods[[m+1L]]@ParTable # partable h0
                af <- lav_partable_constraints_function(p1,p0)
                A <- lav_func_jacobian_complex(func=af, x=mods[[m   ]]@Fit@x)

                trace.UGamma  <- numeric( ngroups )
                trace.UGamma2 <- numeric( ngroups )
                for(g in 1:ngroups) {
                    U <- WLS.V[[g]]  %*% Delta1[[g]] %*%
                         (P1.inv %*% t(A) %*% solve(A %*% P1.inv %*% t(A)) %*% A %*% P1.inv) %*% t(Delta1[[g]]) %*% WLS.V[[g]]
                    UG <- U %*% Gamma[[g]]; tUG <- t(UG)
                    trace.UGamma[g]  <- ntotal/nobs[[g]] * sum( U * Gamma[[g]] )
                    trace.UGamma2[g] <- ntotal/nobs[[g]] * sum( UG * tUG )
                }

                tr.M  <- sum(trace.UGamma)
                tr2.M <- sum(trace.UGamma2)

                # adjust Df.delta?
                if(TEST == "mean.var.adjusted") {
                    # NOT needed for scaled.shifted; see
                    # see 'Simple Second Order Chi-Square Correction' 2010 paper
                    # on www.statmodel.com, section 4
                    # Df.delta[m+1L] <- floor((tr.M1^2 / tr2.M1) + 0.5)
                    Df.delta[m+1L] <- tr.M^2 / tr2.M
                } 

                scaling.factor <- tr.M / Df.delta[m+1L]
                if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
                STAT.delta[m+1L] <- STAT.delta[m+1L]/scaling.factor
            }
        } 
    }

    # Pvalue
    Pvalue.delta <- pchisq(STAT.delta, Df.delta, lower.tail = FALSE)

    aic <- bic <- rep(NA, length(mods))
    if(estimator == "ML") {
        aic <- sapply(mods, FUN=AIC)
        bic <- sapply(mods, FUN=BIC)
    }

    val <- data.frame(Df = Df,
                      AIC = aic,
                      BIC = bic,
                      Chisq = STAT,
                      "Chisq diff" = STAT.delta,
                      "Df diff" = Df.delta,
                      "Pr(>Chisq)" = Pvalue.delta,
                      row.names = names(mods), 
                      check.names = FALSE)

    if(type == "chisq") {
        if(scaled) {
            attr(val, "heading") <- 
                paste("Scaled Chi Square Difference Test (test = ",
                      TEST, ")\n", sep="")
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

