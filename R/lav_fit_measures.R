# user-visible function to extract the fit measures
# output can be 1) vector (default), 2) list, 3) matrix, or 4) text
# in the latter case, the result will be of class "lavaan.fitMeasures"
# for which the printing is done by print.lavaan.fitMeasures()

# new in 0.6-13:
# the big families are computed in dedicated functions:
# - lav_fit_rmsea_lavobject
# - lav_fit_cfi_lavobject
# - lav_fit_aic_lavojbect
# - lav_residuals_summary

# Note: fitMeasures/fitmeasures are generic functions; they include a "..."
#       so lavaan.mi can add arguments to pass to lavTestLRT() and
#       lavTestLRT.mi() about how to pool chi-squared.

setMethod("fitMeasures", signature(object = "lavaan"),
function(object, fit.measures = "all", baseline.model = NULL,
         output = "vector", ...) {
    lav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model, output = output)
})

setMethod("fitmeasures", signature(object = "lavaan"),
function(object, fit.measures = "all", baseline.model = NULL,
         output = "vector",  ...) {
    lav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model, output = output)
})

lav_fit_measures <- function(object, fit.measures = "all",
                             baseline.model = NULL, output = "vector") {

    # do we have data? (yep, we had to include this check)
    if(object@Data@data.type == "none") {
        stop("lavaan ERROR: fit measures not available if there is no data.")
    }

    # has the model converged?
    if(object@optim$npar > 0L && !object@optim$converged) {
        stop("lavaan ERROR: fit measures not available if model did not converge")
    }

    # do we have a test statistic?
    TEST <- lavInspect(object, "test")
    if(TEST[[1]]$test == "none") {
        stop("lavaan ERROR: fit measures not available if test = \"none\".")
    }

    # check output argument
    if(output %in% c("vector", "horizontal")) {
        output <- "vector"
    } else if(output %in% c("list")) {
        output <- "list"
    } else if(output %in% c("matrix", "vertical")) {
        output <- "matrix"
    } else if(output %in% c("text", "pretty", "summary")) {
        output <- "text"
    } else {
        stop("lavaan ERROR: output should be ", sQuote("vector"),
             ", ", sQuote("list"),
             ", ", sQuote("matrix"), " or ", sQuote("text"))
    }

    # options
    categorical <- object@Model@categorical
    estimator   <- object@Options$estimator
    test        <- lav_utils_get_test(lavobject = object) # single element!
    scaled      <- lav_utils_get_scaled(lavobject = object) # robust/scaled?

    # basic ingredients
    G <- object@Data@ngroups
    X2 <- TEST[[1]]$stat
    df <- TEST[[1]]$df
    if(scaled) {
        X2.scaled <- TEST[[2]]$stat
        df.scaled <- TEST[[2]]$df
    }
    npar <- lav_utils_get_npar(lavobject = object)
    N    <- lav_utils_get_ntotal(lavobject = object) # N vs N-1


    # define 'sets' of fit measures:
    fit.always <- c("npar")

    # basic chi-square test
    fit.chisq <- c("fmin", "chisq", "df", "pvalue")
    if(scaled) {
        fit.chisq <- c(fit.chisq, "chisq.scaled", "df.scaled", "pvalue.scaled",
                       "chisq.scaling.factor")
    }

    # baseline model
    fit.baseline <- c("baseline.chisq", "baseline.df", "baseline.pvalue")
    if(scaled) {
        fit.baseline <- c(fit.baseline, "baseline.chisq.scaled",
                          "baseline.df.scaled", "baseline.pvalue.scaled",
                          "baseline.chisq.scaling.factor")
    }

    fit.cfi.tli <- c("cfi", "tli")
    if(scaled) {
        fit.cfi.tli <- c(fit.cfi.tli, "cfi.scaled", "tli.scaled",
                                      "cfi.robust", "tli.robust")
    }

    # other incremental fit indices
    fit.cfi.other <- c("nnfi", "rfi", "nfi", "pnfi", "ifi", "rni")
    if(scaled) {
        fit.cfi.other <- c(fit.cfi.other, "nnfi.scaled", "rfi.scaled",
                       "nfi.scaled", "pnfi.scaled", "ifi.scaled", "rni.scaled",
                       "nnfi.robust", "rni.robust")
    }
    fit.cfi <- c(fit.baseline, fit.cfi.tli, fit.cfi.other)

    # likelihood based measures
    if(estimator == "MML") {
        fit.logl <- c("logl", "aic", "bic", "ntotal", "bic2")
    } else {
        fit.logl <- c("logl", "unrestricted.logl", "aic", "bic",
                      "ntotal", "bic2")
    }
    if(scaled && test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
        fit.logl <- c(fit.logl, "scaling.factor.h1", "scaling.factor.h0")
    }

    # rmsea
    fit.rmsea <- c("rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue")
    if(scaled) {
        fit.rmsea <- c(fit.rmsea, "rmsea.scaled", "rmsea.ci.lower.scaled",
                       "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled",
                       "rmsea.robust", "rmsea.ci.lower.robust",
                       "rmsea.ci.upper.robust", "rmsea.pvalue.robust")
    }

    # srmr
    if(categorical) {
        fit.srmr <- c("srmr")
        fit.srmr2 <- c("rmr", "rmr_nomean",
                       "srmr", # per default equal to srmr_bentler_nomean
                       "srmr_bentler", "srmr_bentler_nomean",
                       "crmr", "crmr_nomean",
                       "srmr_mplus", "srmr_mplus_nomean")
    } else {
        if(object@Data@nlevels > 1L) {
            fit.srmr  <- c("srmr","srmr_within", "srmr_between")
            fit.srmr2 <- c("srmr","srmr_within", "srmr_between")
        } else {
            fit.srmr <- c("srmr")
            fit.srmr2 <- c("rmr", "rmr_nomean",
                           "srmr", # the default
                           "srmr_bentler", "srmr_bentler_nomean",
                           "crmr", "crmr_nomean",
                           "srmr_mplus", "srmr_mplus_nomean")
        }
    }

    # various
    if(object@Data@nlevels > 1L) {
        fit.other <- ""
    } else if(estimator == "PML") {
        fit.other <- c("cn_05","cn_01","mfi")
        if(!categorical) { # really needed?
            fit.other <- c(fit.other, "ecvi")
        }
    } else {
        fit.other <- c("cn_05","cn_01","gfi","agfi","pgfi","mfi")
        if(!categorical) { # really needed?
            fit.other <- c(fit.other, "ecvi")
        } else {
            fit.other <- c(fit.other, "wrmr")
        }
    }

    # lower case
    fit.measures <- tolower(fit.measures)

    # select 'default' fit measures
    if(length(fit.measures) == 1L) {
        if(fit.measures == "default") {
            if(estimator == "ML" || estimator == "PML") {
                fit.measures <- c(fit.always, fit.chisq, fit.baseline,
                                  fit.cfi.tli, fit.logl,
                                  fit.rmsea, fit.srmr)
            } else if(estimator == "MML") {
                fit.measures <- c(fit.always, fit.logl)
            } else {
                fit.measures <- c(fit.always,
                                  fit.chisq, fit.baseline, fit.cfi.tli,
                                  fit.rmsea, fit.srmr)
            }
        } else if(fit.measures == "all") {
            if(estimator == "ML") {
                fit.measures <- c(fit.always, fit.chisq,
                                  fit.baseline, fit.cfi.tli, fit.cfi.other,
                                  fit.logl, fit.rmsea, fit.srmr2, fit.other)
            } else {
                fit.measures <- c(fit.always, fit.chisq,
                                  fit.baseline, fit.cfi.tli, fit.cfi.other,
                                  fit.rmsea, fit.srmr2, fit.other)
            }
        }
    }

    # catch empty list
    if(length(fit.measures) == 0L) {
        return(list())
    }

    # main container
    indices <- list()
    indices["npar"]   <- npar
    indices["ntotal"] <- object@SampleStats@ntotal
    indices["fmin"]   <- object@optim$fx # note = 0.5 * fmin if ML

    # CHI-SQUARE TEST
    if(any(fit.chisq %in% fit.measures)) {
	    indices["chisq"] <- X2
        indices["df"] <- df
        indices["pvalue"] <- TEST[[1]]$pvalue
        if(scaled) {
            indices["chisq.scaled"] <- X2.scaled
            indices["df.scaled"] <- df.scaled
            indices["chisq.scaling.factor"] <- TEST[[2]]$scaling.factor
            indices["pvalue.scaled"] <- TEST[[2]]$pvalue
        }
    }

    # BASELINE FAMILY
    if(any(fit.cfi %in% fit.measures)) {
        indices <- c(indices,
                     lav_fit_cfi_lavobject(lavobject    = object,
                                           fit.measures = fit.measures,
                                           scaled       = scaled,
                                           test         = test))
    }

    # INFORMATION CRITERIA
    if(any(fit.logl %in% fit.measures)) {
        indices <- c(indices,
                     lav_fit_aic_lavobject(lavobject    = object,
                                           fit.measures = fit.measures,
                                           scaled       = scaled,
                                           test         = test,
                                           estimator    = estimator))
    }

    # RMSEA and friends
    if(any(fit.rmsea %in% fit.measures)) {
        indices <- c(indices,
                     lav_fit_rmsea_lavobject(lavobject    = object,
                                             fit.measures = fit.measures,
                                             scaled       = scaled,
                                             test         = test))
    }

    # SRMR and friends
    if(any(fit.srmr2 %in% fit.measures)) {
        indices <- c(indices,
                     lav_fit_srmr_lavobject(lavobject    = object,
                                            fit.measures = fit.measures))
    }

    # GFI and friends
    fit.gfi <- c("gfi", "agfi", "pgfi")
    if(any(fit.gfi %in% fit.measures)) {
        indices <- c(indices,
                     lav_fit_gfi_lavobject(lavobject    = object,
                                           fit.measures = fit.measures))
    }

    # various: Hoelter Critical N (CN)
    if(any(c("cn_05", "cn_01") %in% fit.measures)) {
        indices["cn_05"] <- lav_fit_cn(X2 = X2, df = df, N = N, alpha = 0.05)
        indices["cn_01"] <- lav_fit_cn(X2 = X2, df = df, N = N, alpha = 0.01)
    }

    # various: WRMR
    if("wrmr" %in% fit.measures) {
        nel <- length(object@SampleStats@WLS.obs[[1]])
        indices["wrmr"] <- lav_fit_wrmr(X2 = X2, nel = nel)
    }

    # various: MFI
    if("mfi" %in% fit.measures) {
        indices["mfi"] <- lav_fit_mfi(X2 = X2, df = df, N = N)
    }

    # various: ECVI
    if("ecvi" %in% fit.measures) {
        indices["ecvi"] <- lav_fit_ecvi(X2 = X2, npar = npar, N = N)
    }


    # keep only what we need
    out <- indices[fit.measures]

    if(all(is.na(names(out)))) { # perhaps, fit.measures = ""
        # nothing left
        return(numeric(0L))
    }

    # select output type
    if(output == "list") {
        # nothing to do
    } else if(output == "vector") {
        out <- unlist(out)
        class(out) <- c("lavaan.vector", "numeric")
    } else if(output == "matrix") {
        out <- as.matrix(unlist(out))
        colnames(out) <- ""
        class(out) <- c("lavaan.matrix", "matrix")
    } else if(output == "text") {
        out <- unlist(out)
        class(out) <- c("lavaan.fitMeasures", "lavaan.vector", "numeric")
    }

    out
}

# print a nice summary of the fit measures
print.lavaan.fitMeasures <- function(x, ..., nd = 3L, add.h0 = TRUE) {

    names.x <- names(x)

    # scaled?
    scaled <- "chisq.scaled" %in% names.x

    # num format
    num.format  <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

    ## TDJ: optionally add h0 model's fit statistic, for lavaan.mi
    if (add.h0 && "chisq" %in% names.x) {

        cat("\nModel Test User Model:\n\n")

        # container three columns
        c1 <- c2 <- c3 <- character(0L)

        c1 <- c(c1, "Test statistic")
        c2 <- c(c2, sprintf(num.format, x["chisq"]))
        c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["chisq.scaled"]), ""))

        c1 <- c(c1, "Degrees of freedom")
        c2 <- c(c2, x["df"])
        c3 <- c(c3, ifelse(scaled,
                       ifelse(x["df.scaled"]%%1 == 0, x["df.scaled"],
                           sprintf(num.format, x["df.scaled"])), ""))

        c1 <- c(c1, "P-value")
        c2 <- c(c2, sprintf(num.format, x["pvalue"]))
        c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["pvalue.scaled"]), ""))

        if(scaled && "chisq.scaling.factor" %in% names.x) {
            c1 <- c(c1, "Scaling correction factor")
            c2 <- c(c2, "")
            c3 <- c(c3, sprintf(num.format, x["chisq.scaling.factor"]))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 35L)
        c2 <- format(c2, width = 16L + max(0,(nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # independence model
    if("baseline.chisq" %in% names.x) {
        cat("\nModel Test Baseline Model:\n\n")

        c1 <- c2 <- c3 <- character(0L)

        c1 <- c(c1, "Test statistic")
        c2 <- c(c2, sprintf(num.format, x["baseline.chisq"]))
        c3 <- c(c3, ifelse(scaled,
                    sprintf(num.format, x["baseline.chisq.scaled"]), ""))

        c1 <- c(c1, "Degrees of freedom")
        c2 <- c(c2, x["baseline.df"])
        c3 <- c(c3,ifelse(scaled, ifelse(x["baseline.df.scaled"]%%1 == 0,
                                  x["baseline.df.scaled"],
                                  sprintf(num.format, x["baseline.df.scaled"])),
                          ""))

        c1 <- c(c1, "P-value")
        c2 <- c(c2, sprintf(num.format, x["baseline.pvalue"]))
        c3 <- c(c3, ifelse(scaled,
                       sprintf(num.format, x["baseline.pvalue.scaled"]), ""))

        if(scaled && "baseline.chisq.scaling.factor" %in% names.x) {
            c1 <- c(c1, "Scaling correction factor")
            c2 <- c(c2, "")
            c3 <- c(c3, sprintf(num.format, x["baseline.chisq.scaling.factor"]))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 35L)
        c2 <- format(c2, width = 16L + max(0,(nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # cfi/tli
    if(any(c("cfi","tli","nnfi","rfi","nfi","ifi","rni","pnfi") %in% names.x)) {
        cat("\nUser Model versus Baseline Model:\n\n")

        c1 <- c2 <- c3 <- character(0L)

        if("cfi" %in% names.x) {
            c1 <- c(c1, "Comparative Fit Index (CFI)")
            c2 <- c(c2, sprintf(num.format, x["cfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["cfi.scaled"]),
                               ""))
        }

        if("tli" %in% names.x) {
            c1 <- c(c1, "Tucker-Lewis Index (TLI)")
            c2 <- c(c2, sprintf(num.format, x["tli"]))
            c3 <- c(c3, ifelse(scaled,
                               sprintf(num.format, x["tli.scaled"]), ""))
        }

        if("cfi.robust" %in% names.x) {
            c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")
            c1 <- c(c1, "Robust Comparative Fit Index (CFI)")
            c2 <- c(c2, "")
            c3 <- c(c3, sprintf(num.format, x["cfi.robust"]))
        }

        if("tli.robust" %in% names.x) {
            c1 <- c(c1, "Robust Tucker-Lewis Index (TLI)")
            c2 <- c(c2, "")
            c3 <- c(c3, sprintf(num.format, x["tli.robust"]))
        }

        if("nnfi" %in% names.x) {
            c1 <- c(c1, "Bentler-Bonett Non-normed Fit Index (NNFI)")
            c2 <- c(c2, sprintf(num.format, x["nnfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["nnfi.robust"]),
                               ""))
        }

        if("nfi" %in% names.x) {
            c1 <- c(c1, "Bentler-Bonett Normed Fit Index (NFI)")
            c2 <- c(c2, sprintf(num.format, x["nfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["nfi.scaled"]),
                               ""))
        }

        if("pnfi" %in% names.x) {
            c1 <- c(c1, "Parsimony Normed Fit Index (PNFI)")
            c2 <- c(c2, sprintf(num.format, x["pnfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["pnfi.scaled"]),
                               ""))
        }

        if("rfi" %in% names.x) {
            c1 <- c(c1, "Bollen's Relative Fit Index (RFI)")
            c2 <- c(c2, sprintf(num.format, x["rfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["rfi.scaled"]),
                               ""))
        }

        if("ifi" %in% names.x) {
            c1 <- c(c1, "Bollen's Incremental Fit Index (IFI)")
            c2 <- c(c2, sprintf(num.format, x["ifi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["ifi.scaled"]),
                               ""))
        }

        if("rni" %in% names.x) {
            c1 <- c(c1, "Relative Noncentrality Index (RNI)")
            c2 <- c(c2, sprintf(num.format, x["rni"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["rni.robust"]),
                               ""))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 43L)
        c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # likelihood
    if("logl" %in% names.x) {
        cat("\nLoglikelihood and Information Criteria:\n\n")

        c1 <- c2 <- c3 <- character(0L)

        c1 <- c(c1, "Loglikelihood user model (H0)")
        c2 <- c(c2, sprintf(num.format, x["logl"]))
        c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["logl"]), ""))
        if(!is.na(x["scaling.factor.h0"])) {
            c1 <- c(c1, "Scaling correction factor")
            c2 <- c(c2, sprintf("  %10s", ""))
            c3 <- c(c3, sprintf(num.format, x["scaling.factor.h0"]))
            c1 <- c(c1, "    for the MLR correction")
            c2 <- c(c2, ""); c3 <- c(c3, "")
        }

        if("unrestricted.logl" %in% names.x) {
            c1 <- c(c1, "Loglikelihood unrestricted model (H1)")
            c2 <- c(c2, sprintf(num.format, x["unrestricted.logl"]))
            c3 <- c(c3, ifelse(scaled,
                           sprintf(num.format, x["unrestricted.logl"]), ""))
            if(!is.na(x["scaling.factor.h1"])) {
                c1 <- c(c1, "Scaling correction factor")
                c2 <- c(c2, sprintf("  %10s", ""))
                c3 <- c(c3, sprintf(num.format, x["scaling.factor.h1"]))
                c1 <- c(c1, "    for the MLR correction")
                c2 <- c(c2, ""); c3 <- c(c3, "")
            }
        }

        c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")
        #c1 <- c(c1, "Number of free parameters")
        #c2 <- c(c2, sprintf("  %10i", x["npar"]))
        #c3 <- c(c3, ifelse(scaled, sprintf("  %10i", x["npar"]), ""))

        c1 <- c(c1, "Akaike (AIC)")
        c2 <- c(c2, sprintf(num.format, x["aic"]))
        c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["aic"]), ""))
        c1 <- c(c1, "Bayesian (BIC)")
        c2 <- c(c2, sprintf(num.format, x["bic"]))
        c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["bic"]), ""))
        if(!is.na(x["bic2"])) {
            c1 <- c(c1, "Sample-size adjusted Bayesian (SABIC)")
            c2 <- c(c2, sprintf(num.format, x["bic2"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["bic2"]), ""))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 39L)
        c2 <- format(c2, width = 12L + max(0,(nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # RMSEA
    if("rmsea" %in% names.x) {
        cat("\nRoot Mean Square Error of Approximation:\n\n")

        c1 <- c2 <- c3 <- character(0L)

        c1 <- c(c1, "RMSEA")
        c2 <- c(c2, sprintf(num.format, x["rmsea"]))
        c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["rmsea.scaled"]), ""))
        if("rmsea.ci.lower" %in% names.x) {
            c1 <- c(c1, "90 Percent confidence interval - lower")
            c2 <- c(c2, sprintf(num.format, x["rmsea.ci.lower"]))
            c3 <- c(c3, ifelse(scaled,
                        sprintf(num.format, x["rmsea.ci.lower.scaled"]), ""))
            c1 <- c(c1, "90 Percent confidence interval - upper")
            c2 <- c(c2, sprintf(num.format, x["rmsea.ci.upper"]))
            c3 <- c(c3, ifelse(scaled,
                        sprintf(num.format, x["rmsea.ci.upper.scaled"]), ""))
        }
        if("rmsea.pvalue" %in% names.x) {
            c1 <- c(c1, "P-value H_0: RMSEA <= 0.05")
            c2 <- c(c2, sprintf(num.format, x["rmsea.pvalue"]))
            c3 <- c(c3, ifelse(scaled,
                  sprintf(num.format, x["rmsea.pvalue.scaled"]), ""))
        }

        # robust
        if("rmsea.robust" %in% names.x) {
            c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")
            c1 <- c(c1, "Robust RMSEA")
            c2 <- c(c2, "")
            c3 <- c(c3, ifelse(scaled,
                        sprintf(num.format, x["rmsea.robust"]), ""))
        }
        if("rmsea.ci.lower.robust" %in% names.x) {
            c1 <- c(c1, "90 Percent confidence interval - lower")
            c2 <- c(c2, "")
            c3 <- c(c3, sprintf(num.format, x["rmsea.ci.lower.robust"]))
            c1 <- c(c1, "90 Percent confidence interval - upper")
            c2 <- c(c2, "")
            c3 <- c(c3, sprintf(num.format, x["rmsea.ci.upper.robust"]))
        }
        if("rmsea.pvalue.robust" %in% names.x) {
            c1 <- c(c1, "P-value H_0: Robust RMSEA <= 0.05")
            c2 <- c(c2, "")
            c3 <- c(c3, ifelse(scaled,
                            sprintf(num.format, x["rmsea.pvalue.robust"]), ""))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 43L)
        c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # SRMR
    if(any(c("rmr","srmr") %in% names.x) && ! "srmr_within" %in% names.x) {
        cat("\nStandardized Root Mean Square Residual:\n\n")

        c1 <- c2 <- c3 <- character(0L)

        if("rmr" %in% names.x) {
            c1 <- c(c1, "RMR")
            c2 <- c(c2, sprintf(num.format, x["rmr"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["rmr"]), ""))
        }
        if("rmr_nomean" %in% names.x) {
            c1 <- c(c1, "RMR (No Mean)")
            c2 <- c(c2, sprintf(num.format, x["rmr_nomean"]))
            c3 <- c(c3, ifelse(scaled,
                               sprintf(num.format, x["rmr_nomean"]), ""))
        }
        if("srmr" %in% names.x) {
            c1 <- c(c1, "SRMR")
            c2 <- c(c2, sprintf(num.format, x["srmr"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["srmr"]), ""))
        }
        if("srmr_nomean" %in% names.x) {
            c1 <- c(c1, "SRMR (No Mean)")
            c2 <- c(c2, sprintf(num.format, x["srmr_nomean"]))
            c3 <- c(c3, ifelse(scaled,
                               sprintf(num.format, x["srmr_nomean"]), ""))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 43L)
        c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # SRMR -- multilevel
    if(any(c("srmr_within","srmr_between") %in% names.x)) {
        cat("\nStandardized Root Mean Square Residual (corr metric):\n\n")

        c1 <- c2 <- c3 <- character(0L)

        if("srmr_within" %in% names.x) {
            c1 <- c(c1, "SRMR (within covariance matrix)")
            c2 <- c(c2, sprintf(num.format, x["srmr_within"]))
            c3 <- c(c3, ifelse(scaled,
                               sprintf(num.format, x["srmr_within"]), ""))
        }
        if("srmr_between" %in% names.x) {
            c1 <- c(c1, "SRMR (between covariance matrix)")
            c2 <- c(c2, sprintf(num.format, x["srmr_between"]))
            c3 <- c(c3, ifelse(scaled,
                               sprintf(num.format, x["srmr_between"]), ""))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 43L)
        c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # WRMR
    if("wrmr" %in% names.x) {
        cat("\nWeighted Root Mean Square Residual:\n\n")

        c1 <- c2 <- c3 <- character(0L)

        if("wrmr" %in% names.x) {
            c1 <- c(c1, "WRMR")
            c2 <- c(c2, sprintf(num.format, x["wrmr"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["wrmr"]), ""))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 43L)
        c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    # Other
    if(any(c("cn_05","cn_01","gfi","agfi","pgfi","mfi") %in% names.x)) {
        cat("\nOther Fit Indices:\n\n")

        c1 <- c2 <- c3 <- character(0L)

        if("cn_05" %in% names.x) {
            c1 <- c(c1, "Hoelter Critical N (CN) alpha = 0.05")
            c2 <- c(c2, sprintf(num.format, x["cn_05"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["cn_05"]), ""))
        }
        if("cn_01" %in% names.x) {
            c1 <- c(c1, "Hoelter Critical N (CN) alpha = 0.01")
            c2 <- c(c2, sprintf(num.format, x["cn_01"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["cn_01"]), ""))
        }
        if(any(c("cn_05", "cn_01") %in% names.x)) {
            c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")
        }
        if("gfi" %in% names.x) {
            c1 <- c(c1, "Goodness of Fit Index (GFI)")
            c2 <- c(c2, sprintf(num.format, x["gfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["gfi"]), ""))
        }
        if("agfi" %in% names.x) {
            c1 <- c(c1, "Adjusted Goodness of Fit Index (AGFI)")
            c2 <- c(c2, sprintf(num.format, x["agfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["agfi"]), ""))
        }
        if("pgfi" %in% names.x) {
            c1 <- c(c1, "Parsimony Goodness of Fit Index (PGFI)")
            c2 <- c(c2, sprintf(num.format, x["pgfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["pgfi"]), ""))
        }
        if(any(c("gfi","agfi","pgfi") %in% names.x)) {
            c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")
        }
        if("mfi" %in% names.x) {
            c1 <- c(c1, "McDonald Fit Index (MFI)")
            c2 <- c(c2, sprintf(num.format, x["mfi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["mfi"]), ""))
        }
        if("mfi" %in% names.x) {
            c1 <- c(c1, ""); c2 <- c(c2, ""); c3 <- c(c3, "")
        }
        if("ecvi" %in% names.x) {
            c1 <- c(c1, "Expected Cross-Validation Index (ECVI)")
            c2 <- c(c2, sprintf(num.format, x["ecvi"]))
            c3 <- c(c3, ifelse(scaled, sprintf(num.format, x["ecvi"]), ""))
        }

        # format c1/c2/c3
        c1 <- format(c1, width = 43L)
        c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
        c3 <- format(c3, width = 8L + nd, justify = "right")

        # create character matrix
        if(scaled) {
            M <- cbind(c1, c2, c3, deparse.level = 0)
        } else {
            M <- cbind(c1, c2, deparse.level = 0)
        }
        colnames(M) <- rep("",  ncol(M))
        rownames(M) <- rep(" ", nrow(M))

        # print
        write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }

    invisible(x)
}
