measurementInvariance <- measurementinvariance <- function(..., 
    strict=FALSE, quiet=FALSE) {

    # check for a group.equal argument in ...
    dotdotdot <- list(...)
    if(!is.null(dotdotdot$group.equal))
        stop("lavaan ERROR: group.equal argument should not be used")

    res <- list()
    # base-line model: configural invariance
    res$fit.configural <- cfa(..., group.equal="")

    # fix loadings across groups
    res$fit.loadings <- cfa(..., group.equal=c("loadings"))

    # fix loadings + intercepts across groups
    res$fit.intercepts <- cfa(..., group.equal=c("loadings", 
                                                 "intercepts"))

    if(strict) {
        # fix loadings + intercepts + residuals
        res$fit.residuals <- cfa(..., group.equal=c("loadings", 
                                                    "intercepts", 
                                                    "residuals"))

        # fix loadings + residuals + intercepts + means
        res$fit.means <- cfa(..., group.equal=c("loadings", 
                                                "intercepts", 
                                                "residuals",
                                                "means"))
    } else {
        # fix loadings + intercepts + means
        res$fit.means <- cfa(..., group.equal=c("loadings", 
                                                "intercepts", 
                                                 "means"))
    }

    if(!quiet) {
        cat("\nMeasurement invariance tests:\n")
        cat("\nModel 1: configural invariance:\n")
        printFitLine(res$fit.configural)

        cat("\nModel 2: weak invariance (equal loadings):\n")
        printFitLine(res$fit.loadings)

        cat("\n[Model 1 versus model 2]\n")
        difftest(res$fit.configural, res$fit.loadings)

        cat("\nModel 3: strong invariance (equal loadings + intercepts):\n")
        printFitLine(res$fit.intercepts)
        cat("\n[Model 1 versus model 3]\n")
        difftest(res$fit.configural, res$fit.intercepts)
        cat("\n[Model 2 versus model 3]\n")
        difftest(res$fit.loadings, res$fit.intercepts)


        if(strict) {
            cat("\nModel 4: strict invariance (equal loadings + intercepts + residuals):\n")
            printFitLine(res$fit.residuals)
            cat("\n[Model 1 versus model 4]\n")
            difftest(res$fit.configural, res$fit.residuals)
            cat("\n[Model 3 versus model 4]\n")
            difftest(res$fit.intercepts, res$fit.residuals)
  
            cat("\nModel 5: equal loadings + intercepts + residuals + means:\n")
            printFitLine(res$fit.means,horizontal=TRUE)
            cat("\n[Model 1 versus model 5]\n")
            difftest(res$fit.configural, res$fit.means)
            cat("\n[Model 4 versus model 5]\n")
            difftest(res$fit.residuals, res$fit.means)
        } else {
            cat("\nModel 4: equal loadings + intercepts + means:\n")
            printFitLine(res$fit.means)
            cat("\n[Model 1 versus model 4]\n")
            difftest(res$fit.configural, res$fit.means)
            cat("\n[Model 3 versus model 4]\n")
            difftest(res$fit.intercepts, res$fit.means)
        }
    }
    invisible(res)
}


printFitLine <- function(object, horizontal=TRUE) {

    # which `tests' do we have?
    scaled <- FALSE
    TESTS <- unlist(lapply(object@Fit@test, "[", "test"))
    if(any(c("satorra.bentler", "yuan.bentler") %in% TESTS)) {
        scaled <- TRUE
    }

    if(!scaled) {
        out <- fitMeasures(object, c("chisq", "df", "pvalue",
                                      "cfi", "rmsea", "bic"))
    } else {
        out <- fitMeasures(object, c("chisq.scaled", "df.scaled", 
                                      "pvalue.scaled",
                                      "cfi.scaled", "rmsea.scaled", "bic"))
        names(out) <- c("chisq.scaled", "df", "pvalue",
                        "cfi.scaled", "rmsea.scaled", "bic")
    }
    
    print(out)
}

difftest <- function(model1, model2) {

    # which `tests' do we have for each model?
    model1.scaled <- FALSE
    TESTS <- unlist(lapply(model1@Fit@test, "[", "test"))
    if(any(c("satorra.bentler", "yuan.bentler") %in% TESTS)) {
        model1.scaled <- TRUE
    }
    
    model2.scaled <- FALSE
    TESTS <- unlist(lapply(model2@Fit@test, "[", "test"))
    if(any(c("satorra.bentler", "yuan.bentler") %in% TESTS)) {
        model2.scaled <- TRUE
    }

    if(sum(c(model1.scaled,model2.scaled)) == 2) scaled <- TRUE
    if(sum(c(model1.scaled,model2.scaled)) == 0) scaled <- FALSE
    if(sum(c(model1.scaled,model2.scaled)) == 1) {
        stop("lavaan ERROR: only one of the two models has a scaled test statistic")
    }

    # restricted model is fit0
    if(model1@Fit@test[[1]]$df > model2@Fit@test[[1]]$df) {
        fit0 <- model1
        fit1 <- model2
    } else {
        fit0 <- model2
        fit1 <- model1
    }

    # get fitMeasures
    if(!scaled) {
        fm0 <- fitMeasures(fit0, c("chisq", "df", "cfi"))
        fm1 <- fitMeasures(fit1, c("chisq", "df", "cfi"))
    } else {
        fm0 <- fitMeasures(fit0, c("chisq", "df", "cfi.scaled"))
        fm1 <- fitMeasures(fit1, c("chisq", "df", "cfi.scaled"))
        fm0.scaling <- fit0@Fit@test[[2]]$scaling.factor
        fm1.scaling <- fit1@Fit@test[[2]]$scaling.factor
    }

    result <- numeric(4)
    if(!scaled) {
       # standard chi^2 difference test
       delta.chi     <- (fm0[1] -  fm1[1])
       delta.df      <- (fm0[2] -  fm1[2])
       delta.p.value <- (1 - pchisq(delta.chi, delta.df))
       delta.cfi     <- (fm1[3] -  fm0[3])
       result <- c(delta.chi, delta.df, delta.p.value, delta.cfi)
       names(result) <- c("delta.chisq", "delta.df", "delta.p.value", "delta.cfi")
    } else { # robust!!
        delta.df <- (fm0[2] -  fm1[2])
        naive.chi <- (fm0[1] -  fm1[1])

        # use formula from mplus web note (www.statmodel.com)
        cd <- (fm0[2] * fm0.scaling -
               fm1[2] * fm1.scaling)/delta.df
        delta.chi <- naive.chi/cd

        delta.p.value <- (1 - pchisq(delta.chi, delta.df))
        delta.cfi <- (fm1[3] -  fm0[3])
        result <- c(delta.chi, delta.df, delta.p.value, delta.cfi)
        names(result) <- c("delta.chisq.scaled", "delta.df.scaled", "p.value.scaled", "delta.cfi.scaled")
    }

    print.lavaan.vector(result)

    invisible(result)
}


