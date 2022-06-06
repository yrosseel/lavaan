# TDJ: add "..." to make the method generic, so lavaan.mi can add arguments to
#      pass to lavTestLRT() and lavTestLRT.mi() about how to pool chi-squared.
#      NOT sure this is necessary for the lavaan-method, perhaps only the
#      generic needs "..."?
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

    # has the model converged?
    if(object@optim$npar > 0L && !object@optim$converged) {
        stop("lavaan ERROR: fit measures not available if model did not converge")
    }

    # check output argument
    if(output %in% c("vector", "horizontal")) {
        output <- "vector"
    } else if(output %in% c("matrix", "vertical")) {
        output <- "matrix"
    } else if(output %in% c("text", "pretty")) {
        output <- "text"
    } else {
        stop("lavaan ERROR: output should be ", sQuote("vector"),
             ", ", sQuote("matrix"), " or ", sQuote("text"))
    }

    TEST <- lavInspect(object, "test")

    # do we have a test statistic?
    if(TEST[[1]]$test == "none") {

        # to deal with semTools 0.4-9, we need to check the @Fit@test slot
        #if(object@Fit@test[[1]]$test != "none") {
        #    TEST <- object@Fit@test
        #} else {
            stop("lavaan ERROR: please refit the model with test=\"standard\"")
        #}
    }

    if("all" %in% fit.measures) {
       class.flag <- TRUE
    } else {
       class.flag <- FALSE
    }

    # collect info from the lavaan slots
    GLIST <- object@Model@GLIST

    # N versus N-1
    # this affects BIC, RMSEA, cn_01/05, MFI and ECVI
    # Changed 0.5-15: suggestion by Mark Seeto
    if(object@Options$estimator %in% c("ML","PML","FML") &&
       object@Options$likelihood == "normal") {
        N <- object@SampleStats@ntotal
    } else {
        N <- object@SampleStats@ntotal - object@SampleStats@ngroups
    }

    # Change 0.5-13: take into account explicit equality constraints!!
    # reported by Mark L. Taper (affects AIC and BIC)
    #npar <- object@optim$npar
    npar <- lav_partable_npar(object@ParTable)
    if(nrow(object@Model@con.jac) > 0L) {
        ceq.idx <- attr(object@Model@con.jac, "ceq.idx")
        if(length(ceq.idx) > 0L) {
            neq <- qr(object@Model@con.jac[ceq.idx,,drop=FALSE])$rank
            npar <- npar - neq
        }
    } else if(.hasSlot(object@Model, "ceq.simple.only") &&
              object@Model@ceq.simple.only) {
        npar <- max(object@ParTable$free)
    }

    fx <- object@optim$fx
    fx.group <- object@optim$fx.group
    meanstructure <- object@Model@meanstructure
    categorical   <- object@Model@categorical
    multigroup    <- object@Data@ngroups > 1L
    estimator     <- object@Options$estimator
    test          <- object@Options$test
    # 0.6.5: for now, we make sure that 'test' is a single element
    if(length(test) > 1L) {
        standard.idx <- which(test == "standard")
        if(length(standard.idx) > 0L) {
            test <- test[-standard.idx]
        }
        if(length(test) > 1L) {
            # only retain the first one
            test <- test[1]
        }
    }

    G <- object@Data@ngroups  # number of groups
    X2 <- TEST[[1]]$stat
    df <- TEST[[1]]$df

    # fit stat and df are NA (perhaps test="none"?), try again:
    if(is.na(df)) {
        df <- lav_partable_df(object@ParTable)
        if(nrow(object@Model@con.jac) > 0L) {
            df <- ( df + length(attr(object@Model@con.jac, "ceq.idx")) )
        }
    }
    #if(is.na(X2) && is.finite(fx)) {
    #
    #}

    if(test %in% c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus",
                   "mean.var.adjusted", "scaled.shifted")) {
        scaled <- TRUE
    } else {
        scaled <- FALSE
    }

    # scaled X2
    if(scaled) {
        X2.scaled <- TEST[[2]]$stat
        df.scaled <- TEST[[2]]$df
    }

    # define 'sets' of fit measures:
    fit.always <- c("npar")

    # basic chi-square test
    fit.chisq <- c("fmin", "chisq", "df", "pvalue")
    if(scaled) {
        fit.chisq <- c(fit.chisq, "chisq.scaled", "df.scaled", "pvalue.scaled",
                       "chisq.scaling.factor")
    }

    # basline model
    fit.baseline <- c("baseline.chisq", "baseline.df", "baseline.pvalue")
    if(scaled) {
        fit.baseline <- c(fit.baseline, "baseline.chisq.scaled",
                          "baseline.df.scaled", "baseline.pvalue.scaled",
                          "baseline.chisq.scaling.factor")
    }

    # cfi/tli
    fit.cfi.tli <- c("cfi", "tli")
    if(scaled) {
        fit.cfi.tli <- c(fit.cfi.tli, "cfi.scaled", "tli.scaled",
                                      "cfi.robust", "tli.robust")
    }

    # more incremental fit indices
    fit.incremental <- c("cfi", "tli", "nnfi", "rfi", "nfi", "pnfi",
                         "ifi", "rni")
    if(scaled) {
        fit.incremental <- c(fit.incremental, "cfi.scaled", "tli.scaled",
                             "cfi.robust", "tli.robust",
                             "nnfi.scaled", "nnfi.robust",
                             "rfi.scaled", "nfi.scaled",
                             "ifi.scaled", "rni.scaled", "rni.robust")
    }

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

    # table
    #if(categorical) {
    #    # FIXME: Cp: no exo, all ordinal!
    #    fit.table <- c("c_p", "c_p.df", "c_p.p.value",
    #                   "c_f", "c_f.df", "c_f.p.value",
    #                   "rpat.observed", "rpat.total", "rpat.empty",
    #                   "c_m", "c_m.df", "c_m.p.value")
    #} else {
        fit.table <- character(0L)
    #}

    # various
    if(object@Data@nlevels > 1L) {
        fit.other <- ""
    } else if(estimator == "PML") {
        fit.other <- c("cn_05","cn_01","mfi")
        if(!categorical && G == 1) {
            fit.other <- c(fit.other, "ecvi")
        }
    } else {
        fit.other <- c("cn_05","cn_01","gfi","agfi","pgfi","mfi")
        if(!categorical && G == 1) {
            fit.other <- c(fit.other, "ecvi")
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
                                  fit.rmsea, fit.srmr, fit.table)
                if(object@Options$mimic == "Mplus") {
                    fit.measures <- c(fit.measures, "wrmr")
                }
            }
        } else if(fit.measures == "all") {
            if(estimator == "ML") {
                fit.measures <- c(fit.always,
                                  fit.chisq, fit.baseline, fit.incremental,
                                  fit.logl, fit.rmsea, fit.srmr2, fit.other)
            } else {
                fit.measures <- c(fit.always,
                                  fit.chisq, fit.baseline, fit.incremental,
                                  fit.rmsea, fit.srmr2, fit.other, fit.table)
            }
        }
    }

    # main container
    indices <- list()

    if("npar" %in% fit.measures) {
        indices["npar"] <- npar
    }

    # Chi-square value estimated model (H0)
    if(any(c("fmin", "chisq", "chisq.scaled",
             "chisq.scaling.factor") %in% fit.measures)) {
    indices["fmin"] <- fx
	indices["chisq"] <- X2
        if(scaled) {
            indices["chisq.scaled"] <- X2.scaled
            indices["chisq.scaling.factor"] <- TEST[[2]]$scaling.factor
        }
    }
    if(any(c("df", "df.scaled") %in% fit.measures)) {
        indices["df"] <- df
        if(scaled) {
            indices["df.scaled"] <- df.scaled
        }
    }
    if(any(c("pvalue", "pvalue.scaled") %in% fit.measures)) {
        indices["pvalue"] <- TEST[[1]]$pvalue
        if(scaled) {
            indices["pvalue.scaled"] <- TEST[[2]]$pvalue
        }
    }


    if(any(c("cfi", "cfi.scaled", "cfi.robust",
             "tli", "tli.scaled", "tli.robust",
             "nnfi", "nnfi.scaled", "nnfi.robust",
             "pnfi", "pnfi.scaled",
             "rfi", "rfi.scaled", "nfi", "nfi.scaled",
             "ifi", "ifi.scaled", "rni", "rni.scaled", "rni.robust",
             "baseline.chisq", "baseline.chisq.scaled",
             "baseline.pvalue", "baseline.pvalue.scaled") %in% fit.measures)) {

        baseline.test <- NULL

        # we use the following priority:
        # 1. user-provided baseline model
        # 2. baseline model in @external slot
        # 3. baseline model in @baseline slot
        # 4. nothing -> compute independence model

        # 1. user-provided baseline model
        if( !is.null(baseline.model) ) {
            baseline.test <-
                lav_fit_measures_check_baseline(fit.indep = baseline.model,
                                                object    = object)
        # 2. baseline model in @external slot
        } else if( !is.null(object@external$baseline.model) ) {
            fit.indep <- object@external$baseline.model
            baseline.test <-
                lav_fit_measures_check_baseline(fit.indep = fit.indep,
                                                object    = object)
        # 3. internal @baseline slot
        } else if( .hasSlot(object, "baseline") &&
                   length(object@baseline) > 0L &&
                   !is.null(object@baseline$test) ) {
            baseline.test <- object@baseline$test
        # 4. (re)compute independence model
        } else {
            fit.indep <- try(lav_object_independence(object), silent = TRUE)
            baseline.test <-
                lav_fit_measures_check_baseline(fit.indep = fit.indep,
                                                object    = object)
        }

        if(!is.null(baseline.test)) {
            X2.null <- baseline.test[[1]]$stat
            df.null <- baseline.test[[1]]$df
            if(scaled) {
                X2.null.scaled <- baseline.test[[2]]$stat
                df.null.scaled <- baseline.test[[2]]$df
            }
        } else {
            X2.null <- df.null <- as.numeric(NA)
            X2.null.scaled <- df.null.scaled <- as.numeric(NA)
        }

        # check for NAs
        if(is.na(X2) || is.na(df) || is.na(X2.null) || is.na(df.null)) {
            indices[fit.incremental] <- as.numeric(NA)
        } else {

            if("baseline.chisq" %in% fit.measures) {
                indices["baseline.chisq"] <- X2.null
                if(scaled) {
                    indices["baseline.chisq.scaled"] <- X2.null.scaled
                }
            }
            if("baseline.df" %in% fit.measures) {
                indices["baseline.df"] <- df.null
                if(scaled) {
                    indices["baseline.df.scaled"] <- df.null.scaled
                }
            }
            if("baseline.pvalue" %in% fit.measures) {
                indices["baseline.pvalue"] <- baseline.test[[1]]$pvalue
                if(scaled) {
                    indices["baseline.pvalue.scaled"] <-
                        baseline.test[[2]]$pvalue
                }
            }
            if("baseline.chisq.scaling.factor" %in% fit.measures) {
                indices["baseline.chisq.scaling.factor"] <-
                    baseline.test[[2]]$scaling.factor
            }

            # CFI - comparative fit index (Bentler, 1990)
            if("cfi" %in% fit.measures) {
                t1 <- max( c(X2 - df, 0) )
                t2 <- max( c(X2 - df, X2.null - df.null, 0) )
                if(isTRUE(all.equal(t1,0)) &&
                   isTRUE(all.equal(t2,0))) {
                    indices["cfi"] <- 1
                } else {
                    indices["cfi"] <- 1 - t1/t2
                }
            }
            if("cfi.scaled" %in% fit.measures) {
                t1 <- max( c(X2.scaled - df.scaled, 0) )
                t2 <- max( c(X2.scaled - df.scaled,
                             X2.null.scaled - df.null.scaled, 0) )
                if(is.na(t1) || is.na(t2)) {
                    indices["cfi.scaled"] <- NA
                } else if(isTRUE(all.equal(t1,0)) &&
                          isTRUE(all.equal(t2,0))) {
                    indices["cfi.scaled"] <- 1
                } else {
                    indices["cfi.scaled"] <- 1 - t1/t2
                }
            }
            if("cfi.robust" %in% fit.measures) {

                if(TEST[[2]]$test %in%
                  c("satorra.bentler", "yuan.bentler.mplus", "yuan.bentler")) {

                    # see Brosseau-Liard & Savalei MBR 2014, equation 15

                    # what to do if X2 = 0 and df = 0? in this case,
                    # the scaling factor (ch) will be NA, and we get NA
                    # (instead of 1)
                    if(X2 < .Machine$double.eps && df == 0) {
                        ch <- 0
                    } else {
                        ch <- TEST[[2]]$scaling.factor
                    }
                    cb <- baseline.test[[2]]$scaling.factor

                    t1 <- max( c(X2 - (ch*df), 0) )
                    t2 <- max( c(X2 - (ch*df), X2.null - (cb*df.null), 0) )
                    if(is.na(t1) || is.na(t2)) {
                        indices["cfi.robust"] <- NA
                    } else if(isTRUE(all.equal(t1,0)) &&
                              isTRUE(all.equal(t2,0))) {
                        indices["cfi.robust"] <- 1
                    } else {
                        indices["cfi.robust"] <- 1 - t1/t2
                    }
                } else {
                    indices["cfi.robust"] <- NA
                }
            }

            # RNI - relative noncentrality index (McDonald & Marsh, 1990)
            # same as CFI, but without the max(0,)
            if("rni" %in% fit.measures) {
                t1 <- X2 - df
                t2 <- X2.null - df.null
                if(isTRUE(all.equal(t2,0))) {
                    RNI <- NA
                } else {
                    RNI <- 1 - t1/t2
                }
                indices["rni"] <- RNI
            }
            if("rni.scaled" %in% fit.measures) {
                t1 <- X2.scaled - df.scaled
                t2 <- X2.null.scaled - df.null.scaled
                if(is.na(t1) || is.na(t2)) {
                    RNI <- NA
                } else if(isTRUE(all.equal(t2,0))) {
                    RNI <- NA
                } else {
                    RNI <- 1 - t1/t2
                }
                indices["rni.scaled"] <- RNI
            }
            if("rni.robust" %in% fit.measures) {

                if(TEST[[2]]$test %in%
                   c("satorra.bentler", "yuan.bentler.mplus", "yuan.bentler")) {
                    # see Brosseau-Liard & Savalei MBR 2014, equation 15

                    # what to do if X2 = 0 and df = 0? in this case,
                    # the scaling factor (ch) will be NA, and we get NA
                    # (instead of 1)
                    if(X2 < .Machine$double.eps && df == 0) {
                        ch <- 0
                    } else {
                        ch <- TEST[[2]]$scaling.factor
                    }
                    cb <- baseline.test[[2]]$scaling.factor

                    t1 <- X2 - ch*df
                    t2 <- X2.null - cb*df.null
                    if(is.na(t1) || is.na(t2)) {
                        RNI <- NA
                    } else if(isTRUE(all.equal(t2,0))) {
                        RNI <- NA
                    } else {
                        RNI <- 1 - t1/t2
                    }
                    indices["rni.robust"] <- RNI
                } else {
                    indices["rni.robust"] <- NA
                }
            }

            # TLI - Tucker-Lewis index (Tucker & Lewis, 1973)
            # same as
            # NNFI - nonnormed fit index (NNFI, Bentler & Bonett, 1980)
            if("tli" %in% fit.measures || "nnfi" %in% fit.measures) {
                # note: formula in lavaan <= 0.5-20:
                # t1 <- X2.null/df.null - X2/df
                # t2 <- X2.null/df.null - 1
                # if(t1 < 0 && t2 < 0) {
                #    TLI <- 1
                #} else {
                #    TLI <- t1/t2
                #}
                # note: TLI original formula was in terms of fx/df, not X2/df
                # then, t1 <- fx_0/df.null - fx/df
                #       t2 <- fx_0/df.null - 1/N (or N-1 for wishart)

                # note: in lavaan 0.5-21, we use the alternative formula:
                # TLI <- 1 - ((X2 - df)/(X2.null - df.null) * df.null/df)
                # - this one has the advantage that a 'robust' version
                #   can be derived; this seems non-trivial for the original one
                # - unlike cfi, we do not use 'max(0, )' for t1 and t2
                #   therefore, t1 can go negative, and TLI can be > 1
                t1 <- (X2 - df)*df.null
                t2 <- (X2.null - df.null)*df
                if(df > 0 && abs(t2) > 0) {
                    indices["tli"] <- indices["nnfi"] <- 1 - t1/t2
                } else {
                    indices["tli"] <- indices["nnfi"] <- 1
                }
            }

            if("tli.scaled" %in% fit.measures ||
               "nnfi.scaled" %in% fit.measures) {

                t1 <- (X2.scaled - df.scaled)*df.null.scaled
                t2 <- (X2.null.scaled - df.null.scaled)*df.scaled
                if(is.na(t1) || is.na(t2)) {
                    indices["tli.scaled"] <- indices["nnfi.scaled"] <- NA
                } else if(df > 0 && abs(t2) > 0) {
                    indices["tli.scaled"] <- indices["nnfi.scaled"] <- 1 - t1/t2
                } else {
                    indices["tli.scaled"] <- indices["nnfi.scaled"] <- 1
                }
            }

            if("tli.robust" %in% fit.measures ||
               "nnfi.robust" %in% fit.measures) {

                if(TEST[[2]]$test %in%
                   c("satorra.bentler", "yuan.bentler.mplus", "yuan.bentler")) {
                    #  see Brosseau-Liard & Savalei MBR 2014, equation 16

                    # what to do if X2 = 0 and df = 0? in this case,
                    # the scaling factor (ch) will be NA, and we get NA
                    # (instead of 1)
                    if(X2 < .Machine$double.eps && df == 0) {
                        ch <- 0
                    } else {
                        ch <- TEST[[2]]$scaling.factor
                    }
                    cb <- baseline.test[[2]]$scaling.factor

                    t1 <- (X2 - ch*df)*df.null
                    t2 <- (X2.null - cb*df.null)*df
                    if(is.na(t1) || is.na(t2)) {
                        indices["tli.robust"] <- indices["nnfi.robust"] <- NA
                    } else if(df > 0 && abs(t2) > 0) {
                        indices["tli.robust"] <- indices["nnfi.robust"] <- 1 - t1/t2
                    } else {
                        indices["tli.robust"] <- indices["nnfi.robust"] <- 1
                    }
                } else {
                    indices["tli.robust"] <- indices["nnfi.robust"] <- NA
                }
            }



            # RFI - relative fit index (Bollen, 1986; Joreskog & Sorbom 1993)
            if("rfi" %in% fit.measures) {
                if(df > df.null) {
                    RLI <- as.numeric(NA)
                } else if(df > 0 && df.null > 0) {
                    t1 <- X2.null/df.null - X2/df
                    t2 <- X2.null/df.null
                    if(t1 < 0 || t2 < 0) {
                        RLI <- 1
                    } else {
                        RLI <- t1/t2
                    }
                } else {
                   RLI <- 1
                }
                indices["rfi"] <- RLI
            }
            if("rfi.scaled" %in% fit.measures) {
                if(df > df.null) {
                    RLI <- as.numeric(NA)
                } else if(df > 0) {
                    t1 <- X2.null.scaled/df.null.scaled - X2.scaled/df.scaled
                    t2 <- X2.null.scaled/df.null.scaled
                    if(is.na(t1) || is.na(t2)) {
                        RLI <- NA
                    } else if(t1 < 0 || t2 < 0) {
                        RLI <- 1
                    } else {
                        RLI <- t1/t2
                    }
                } else {
                   RLI <- 1
                }
                indices["rfi.scaled"] <- RLI
            }

            # NFI - normed fit index (Bentler & Bonett, 1980)
            if("nfi" %in% fit.measures) {
                if(df > df.null || isTRUE(all.equal(X2.null,0))) {
                    NFI <- as.numeric(NA)
                } else if(df > 0) {
                    t1 <- X2.null - X2
                    t2 <- X2.null
                    NFI <- t1/t2
                } else {
                    NFI <- 1
                }
                indices["nfi"] <- NFI
            }
            if("nfi.scaled" %in% fit.measures) {
                if(df > df.null || isTRUE(all.equal(X2.null.scaled,0))) {
                    NFI <- as.numeric(NA)
                } else {
                    t1 <- X2.null.scaled - X2.scaled
                    t2 <- X2.null.scaled
                    NFI <- t1/t2
                }
                indices["nfi.scaled"] <- NFI
            }

            # PNFI - Parsimony normed fit index (James, Mulaik & Brett, 1982)
            if("pnfi" %in% fit.measures) {
                if(df.null > 0 && X2.null > 0) {
                    t1 <- X2.null - X2
                    t2 <- X2.null
                    PNFI <- (df/df.null) * t1/t2
                } else {
                    PNFI <- as.numeric(NA)
                }
                indices["pnfi"] <- PNFI
            }
            if("pnfi.scaled" %in% fit.measures) {
                if(df.null > 0 && X2.null.scaled > 0) {
                    t1 <- X2.null.scaled - X2.scaled
                    t2 <- X2.null.scaled
                    PNFI <- (df/df.null) * t1/t2
                } else {
                    PNFI <- as.numeric(NA)
                }
                indices["pnfi.scaled"] <- PNFI
            }

            # IFI - incremental fit index (Bollen, 1989; Joreskog & Sorbom, 1993)
            if("ifi" %in% fit.measures) {
                t1 <- X2.null - X2
                t2 <- X2.null - df
                if(t2 < 0) {
                    IFI <- 1
                } else if(isTRUE(all.equal(t2,0))) {
                    IFI <- as.numeric(NA)
                } else {
                    IFI <- t1/t2
                }
                indices["ifi"] <- IFI
            }
            if("ifi.scaled" %in% fit.measures) {
                t1 <- X2.null.scaled - X2.scaled
                t2 <- X2.null.scaled - df.scaled
                if(is.na(t2)) {
                    IFI <- NA
                } else if(t2 < 0) {
                    IFI <- 1
                } else if(isTRUE(all.equal(t2,0))) {
                    IFI <- as.numeric(NA)
                } else {
                    IFI <- t1/t2
                }
                indices["ifi.scaled"] <- IFI
            }
        }
    }

    if("logl" %in% fit.measures ||
       "unrestricted.logl" %in% fit.measures ||
       "aic" %in% fit.measures ||
       "bic" %in% fit.measures ||
       "bic2" %in% fit.measures) {

        if(estimator == "ML" || estimator == "MML") {

            # do we have a @h1 slot?
            if(.hasSlot(object, "h1") && length(object@h1) > 0L) {
                logl.H1.group <- object@h1$loglik.group
                logl.H1       <- object@h1$loglik
            } else {
                out <- lav_h1_logl(lavdata = object@Data,
                                   lavsamplestats = object@SampleStats,
                                   lavoptions = object@Options)

                logl.H1.group <- out$loglik.group
                logl.H1       <- out$loglik
            }

            if("unrestricted.logl" %in% fit.measures) {
                indices["unrestricted.logl"] <- logl.H1
            }

            # logl H0
            if(.hasSlot(object, "loglik")) {
                logl.H0.group <- object@loglik$loglik.group
                logl.H0       <- object@loglik$loglik
                AIC           <- object@loglik$AIC
                BIC           <- object@loglik$BIC
                BIC2          <- object@loglik$BIC2
            } else {
                out <- lav_model_loglik(lavdata        = object@Data,
                                        lavsamplestats = object@SampleStats,
                                        lavimplied     = object@implied,
                                        lavmodel       = object@Model,
                                        lavoptions     = object@Options)
                logl.H0.group <- out$loglik.group
                logl.H0       <- out$loglik
                AIC           <- out$AIC
                BIC           <- out$BIC
                BIC2          <- out$BIC2
            }

            if("logl" %in% fit.measures) {
                indices["logl"] <- logl.H0
            }

            # AIC
            if("aic" %in% fit.measures) {
                indices["aic"] <- AIC
            }

            # BIC
            if("bic" %in% fit.measures ||
               "bic2" %in% fit.measures) {
                indices["bic"] <- BIC
                indices["bic2"] <- BIC2
            }

            # scaling factor for MLR
            if(test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
                indices["scaling.factor.h1"] <-
                    TEST[[2]]$scaling.factor.h1
                indices["scaling.factor.h0"] <-
                    TEST[[2]]$scaling.factor.h0
            }
        } # ML
        else { # no ML!
            if("logl" %in% fit.measures) {
                indices["logl"] <- as.numeric(NA)
            }
            if("unrestricted.logl" %in% fit.measures) {
                indices["unrestricted.logl"] <- as.numeric(NA)
            }
            if("aic" %in% fit.measures) {
                indices["aic"] <- as.numeric(NA)
            }
            if("bic" %in% fit.measures) {
                indices["bic"] <- as.numeric(NA)
            }
            if("bic2" %in% fit.measures) {
                indices["bic2"] <- as.numeric(NA)
            }
        }
    }

    N.RMSEA <- max(N, X2*4) # FIXME: good strategy??
    if(any(c("rmsea","rmsea.scaled","rmsea.robust") %in% fit.measures)) {
        # RMSEA
        # - RMSEA.scaled replaces X2 by X2.scaled (which is not ok)
        # - RMSEA.robust uses the formula from Broseau-Liard, Savalei & Li
        #   (2012) paper (see eq 8)
        if(is.na(X2) || is.na(df)) {
            RMSEA <- RMSEA.scaled <- RMSEA.robust <- as.numeric(NA)
        } else if(df > 0) {
            if(scaled) {
                d <- sum(TEST[[2]]$trace.UGamma)
                if(is.na(d) || d==0) d <- NA

                # scaling factor
                c.hat <- TEST[[2]]$scaling.factor
            }

            RMSEA <- sqrt( max( c((X2/N)/df - 1/N, 0) ) )
            if(scaled && test != "scaled.shifted") {
                RMSEA.scaled <-
                     sqrt( max( c((X2/N)/d - 1/N, 0) ) )
                RMSEA.robust <-
                     sqrt( max( c((X2/N)/df - c.hat/N, 0) ) )
            } else if(test == "scaled.shifted") {
                RMSEA.scaled <-
                     sqrt( max( c((X2.scaled/N)/df - 1/N, 0)))
                RMSEA.robust <-
                     sqrt( max( c((X2/N)/df - c.hat/N, 0) ) )
            }

            #  multiple group correction
            #  note: recent builds of EQS also use this 'correction'
            #        perhaps we should have an option to obtain the 'old' one
            #if(object@Options$mimic %in% c("Mplus", "lavaan")) {
                RMSEA <- RMSEA * sqrt(G)
                if(scaled) {
                    RMSEA.scaled <- RMSEA.scaled * sqrt(G)
                    RMSEA.robust <- RMSEA.robust * sqrt(G)
                }
            #}

        } else {
            RMSEA <- RMSEA.scaled <- RMSEA.robust <- 0
        }
        indices["rmsea"] <- RMSEA
        if(scaled) {
            indices["rmsea.scaled"] <- RMSEA.scaled
            if(TEST[[2]]$test %in% c("satorra.bentler", "yuan.bentler.mplus",
                                     "yuan.bentler")) {
                indices["rmsea.robust"] <- RMSEA.robust
            } else {
                indices["rmsea.robust"] <- NA
            }
        }
    }

    if("rmsea.ci.lower" %in% fit.measures) {
        lower.lambda <- function(lambda) {
            (pchisq(X2, df=df, ncp=lambda) - 0.95)
        }
        if(is.na(X2) || is.na(df)) {
            indices["rmsea.ci.lower"] <- NA
        } else if(df < 1 || lower.lambda(0) < 0.0) {
            indices["rmsea.ci.lower"] <- 0
        } else {
            lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=X2)$root,
                            silent=TRUE)
            if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
            #if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                GG <- 0
                indices["rmsea.ci.lower"] <-
                    sqrt( lambda.l/((N-GG)*df) ) * sqrt(G)
            #} else {
            #    indices["rmsea.ci.lower"] <- sqrt( lambda.l/(N*df) )
            #}
        }
    }

    if("rmsea.ci.lower.scaled" %in% fit.measures ||
       "rmsea.ci.lower.robust" %in% fit.measures) {
        if(test == "scaled.shifted") {
            XX2 <- X2.scaled
            df2 <- df
        } else {
            XX2 <- X2
            df2 <- sum(TEST[[2]]$trace.UGamma)
        }
        lower.lambda <- function(lambda) {
            (pchisq(XX2, df=df2, ncp=lambda) - 0.95)
        }
        if(is.na(XX2) || is.na(df2)) {
            indices["rmsea.ci.lower.scaled"] <-
            indices["rmsea.ci.lower.robust"] <- NA
        #} else if(df < 1 || df2 < 1 || lower.lambda(0) < 0.0) {
        } else if(df < 1 || lower.lambda(0) < 0.0) { # no longer df2<1 check
            indices["rmsea.ci.lower.scaled"] <-
            indices["rmsea.ci.lower.robust"] <- 0
        } else {
            # 'scaled'
            lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=XX2)$root,
                            silent=TRUE)
            if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
            #if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                indices["rmsea.ci.lower.scaled"] <-
                        sqrt( lambda.l/(N*df2) ) * sqrt(G)
            #} else {
            #    # no multiple group correction
            #    indices["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*df2) )
            #}

            if(TEST[[2]]$test %in% c("satorra.bentler", "yuan.bentler.mplus",
                                     "yuan.bentler")) {
                # robust
                XX2 <- X2.scaled
                df2 <- df
                # scaling factor
                c.hat <- TEST[[2]]$scaling.factor

                lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=XX2)$root,
                                silent=TRUE)
                if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
                #if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                    indices["rmsea.ci.lower.robust"] <-
                            sqrt( (c.hat*lambda.l)/(N*df2) ) * sqrt(G)
                #} else {
                #    # no multiple group correction
                #    indices["rmsea.ci.lower.robust"] <-
                #        sqrt( (c.hat*lambda.l)/(N*df2) )
                #}
            } else {
                indices["rmsea.ci.lower.robust"] <- NA
            }
        }
    }

    if("rmsea.ci.upper" %in% fit.measures) {
        upper.lambda <- function(lambda) {
            (pchisq(X2, df=df, ncp=lambda) - 0.05)
        }
        if(is.na(X2) || is.na(df)) {
            indices["rmsea.ci.upper"] <- NA
        } else if(df < 1 || upper.lambda(N.RMSEA) > 0 || upper.lambda(0) < 0) {
            indices["rmsea.ci.upper"] <- 0
        } else {
            lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=N.RMSEA)$root,
                            silent=TRUE)
            if(inherits(lambda.u, "try-error")) { lambda.u <- NA }
            #if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                GG <- 0
                indices["rmsea.ci.upper"] <-
                    sqrt( lambda.u/((N-GG)*df) ) * sqrt(G)
            #} else {
            #    indices["rmsea.ci.upper"] <- sqrt( lambda.u/(N*df) )
            #}
        }
    }

    if("rmsea.ci.upper.scaled" %in% fit.measures ||
       "rmsea.ci.upper.robust" %in% fit.measures) {
        if(test == "scaled.shifted") {
            XX2 <- X2.scaled
            df2 <- df
        } else {
            XX2 <- X2
            df2 <- sum(TEST[[2]]$trace.UGamma)
        }
        upper.lambda <- function(lambda) {
            (pchisq(XX2, df=df2, ncp=lambda) - 0.05)
        }
        if(is.na(XX2) || is.na(df2)) {
            indices["rmsea.ci.upper.scaled"] <-
            indices["rmsea.ci.upper.robust"] <- NA
        } else if(df < 1 || upper.lambda(N.RMSEA) > 0 ||
                            upper.lambda(0) < 0) {
            indices["rmsea.ci.upper.scaled"] <-
            indices["rmsea.ci.upper.robust"] <- 0
        } else {
            # 'scaled'
            lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=N.RMSEA)$root,
                            silent=TRUE)
            if(inherits(lambda.u, "try-error")) { lambda.u <- NA }
            #if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                indices["rmsea.ci.upper.scaled"] <-
                    sqrt( lambda.u/(N*df2) ) * sqrt(G)
            #} else {
            #    # no multiple group correction
            #    indices["rmsea.ci.upper.scaled"] <-
            #        sqrt( lambda.u/(N*df2) )
            #}

            if(TEST[[2]]$test %in% c("satorra.bentler", "yuan.bentler.mplus",
                                     "yuan.bentler")) {
                # robust
                XX2 <- X2.scaled
                df2 <- df
                # scaling factor
                c.hat <- TEST[[2]]$scaling.factor

                lambda.u <- try(uniroot(f=upper.lambda, lower=0,upper=N.RMSEA)$root,
                                silent=TRUE)
                if(inherits(lambda.u, "try-error")) { lambda.u <- NA }
             #   if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                    indices["rmsea.ci.upper.robust"] <-
                        sqrt( (c.hat*lambda.u)/(N*df2) ) * sqrt(G)
             #   } else {
             #       # no multiple group correction
             #       indices["rmsea.ci.upper.robust"] <-
             #           sqrt( (c.hat*lambda.u)/(N*df2) )
             #   }
            } else {
                indices["rmsea.ci.upper.robust"] <- NA
            }
        }
    }

    if("rmsea.pvalue" %in% fit.measures) {
        if(is.na(X2) || is.na(df)) {
            indices["rmsea.pvalue"] <- as.numeric(NA)
        } else if(df > 0) {
            #if(object@Options$mimic %in% c("lavaan","Mplus")) {
                ncp <- N*df*0.05^2/G
                indices["rmsea.pvalue"] <-
                    1 - pchisq(X2, df=df, ncp=ncp)
            #} else {
            #    indices["rmsea.pvalue"] <-
            #        1 - pchisq(X2, df=df, ncp=(N*df*0.05^2))
            #}
        } else {
            indices["rmsea.pvalue"] <- NA # used to be 1 in < 0.5-21
        }
    }

    if("rmsea.pvalue.scaled" %in% fit.measures ||
       "rmsea.pvalue.robust" %in% fit.measures) {
        if(test == "scaled.shifted") {
            XX2 <- X2.scaled
            df2 <- df
        } else {
            XX2 <- X2
            df2 <- sum(TEST[[2]]$trace.UGamma)
        }
        if(is.na(XX2) || is.na(df2)) {
            indices["rmsea.pvalue.scaled"] <-
            indices["rmsea.pvalue.robust"] <- as.numeric(NA)
        } else if(df > 0) {
            # scaled
            #if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                ncp <- N*df2*0.05^2/G
                indices["rmsea.pvalue.scaled"] <-
                    1 - pchisq(XX2, df=df2, ncp=ncp)
            #} else {
            #    indices["rmsea.pvalue.scaled"] <-
            #        1 - pchisq(XX2, df=df2, ncp=(N*df2*0.05^2))
            #}

            if(TEST[[2]]$test %in% c("satorra.bentler", "yuan.bentler.mplus",
                                     "yuan.bentler")) {
                # robust
                XX2 <- X2.scaled
                df2 <- df
                # scaling factor
                c.hat <- TEST[[2]]$scaling.factor

                #if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                #    ncp <- N*(df2/c.hat)*0.05^2/G
                #    indices["rmsea.pvalue.robust"] <-
                #        1 - pchisq(XX2, df=df2, ncp=ncp)
                #} else {
                #    indices["rmsea.pvalue.robust"] <-
                #        1 - pchisq(XX2, df=df2, ncp=(N*(df2/c.hat)*0.05^2))
                #}
                indices["rmsea.pvalue.robust"] <- NA
            } else {
                indices["rmsea.pvalue.robust"] <- NA
            }
        } else {
                indices["rmsea.pvalue.scaled"] <-
                indices["rmsea.pvalue.robust"] <- NA # used to be 1 in < 0.5-21
        }
    }

    if(any(c("rmr","srmr") %in% fit.measures) &&  object@Data@nlevels == 1L) {
        # RMR and SRMR
        rmr.group <- numeric(G)
        rmr_nomean.group <- numeric(G)
        srmr_bentler.group <- numeric(G)
        srmr_bentler_nomean.group <- numeric(G)
        crmr.group <- numeric(G)
        crmr_nomean.group <- numeric(G)
        srmr_mplus.group <- numeric(G)
        srmr_mplus_nomean.group <- numeric(G)

        for(g in 1:G) {
            # observed
            if(!object@SampleStats@missing.flag) {
                if(object@Model@conditional.x) {
                    S <- object@SampleStats@res.cov[[g]]
                    M <- object@SampleStats@res.int[[g]]
                } else {
                    S <- object@SampleStats@cov[[g]]
                    M <- object@SampleStats@mean[[g]]
                }
            } else {
                # EM estimates
                S <- object@SampleStats@missing.h1[[g]]$sigma
                M <- object@SampleStats@missing.h1[[g]]$mu
            }
            nvar <- ncol(S)

            # estimated
            implied <- object@implied
            lavmodel <- object@Model
            Sigma.hat <- if(lavmodel@conditional.x) implied$res.cov[[g]] else implied$cov[[g]]
            Mu.hat    <- if(lavmodel@conditional.x) implied$res.int[[g]] else implied$mean[[g]]

            # unstandardized residuals
            RR <- (S - Sigma.hat)

            # standardized residual covariance matrix
            # this is the Hu and Bentler definition, not the Bollen one!
            # this one is used by EQS
            # and Mplus, but only if information=expected (god knows why)
            sqrt.d <- 1/sqrt(diag(S))
            D <- diag(sqrt.d, ncol=length(sqrt.d))
            R <- D %*% (S - Sigma.hat) %*% D

            # Bollen approach: simply using cov2cor ('residual correlations')
            S.cor <- cov2cor(S)
            Sigma.cor <- cov2cor(Sigma.hat)
            R.cor <- (S.cor - Sigma.cor)

            if(meanstructure) {
                # standardized residual mean vector
                R.mean <- D %*% (M - Mu.hat) # EQS approach!
                RR.mean <- (M - Mu.hat) # not standardized
                R.cor.mean <- M/sqrt(diag(S)) - Mu.hat/sqrt(diag(Sigma.hat))

                e <- nvar*(nvar+1)/2 + nvar
                srmr_bentler.group[g] <-
                    sqrt( (sum(R[lower.tri(R, diag=TRUE)]^2) +
                           sum(R.mean^2))/ e )
                rmr.group[g] <- sqrt( (sum(RR[lower.tri(RR, diag=TRUE)]^2) +
                                       sum(RR.mean^2))/ e )
                crmr.group[g] <-
                    sqrt( (sum(R.cor[lower.tri(R.cor, diag=TRUE)]^2)  +
                           sum(R.cor.mean^2)) / (e - nvar) )
                # see http://www.statmodel.com/download/SRMR.pdf
                srmr_mplus.group[g] <-
                    sqrt( (sum(R.cor[lower.tri(R.cor, diag=FALSE)]^2)  +
                           sum(R.cor.mean^2) +
                           sum(((diag(S) - diag(Sigma.hat))/diag(S))^2)) / e )

                e <- nvar*(nvar+1)/2
                srmr_bentler_nomean.group[g] <-
                    sqrt(  sum( R[lower.tri( R, diag=TRUE)]^2) / e )
                rmr_nomean.group[g] <-
                    sqrt(  sum(RR[lower.tri(RR, diag=TRUE)]^2) / e )
                if((e - nvar) > 0) {
                    crmr_nomean.group[g] <-
                    sqrt(  sum(R.cor[lower.tri(R.cor, diag=TRUE)]^2) / (e - nvar) )
                } else {
                    crmr_nomean.group[g] <- as.numeric(NA)
                }
                srmr_mplus_nomean.group[g] <-
                    sqrt( (sum(R.cor[lower.tri(R.cor, diag=FALSE)]^2)  +
                           sum(((diag(S) - diag(Sigma.hat))/diag(S))^2)) / e )
            } else {
                e <- nvar*(nvar+1)/2
                srmr_bentler_nomean.group[g] <- srmr_bentler.group[g] <-
                    sqrt( sum(R[lower.tri(R, diag=TRUE)]^2) / e )
                rmr_nomean.group[g] <- rmr.group[g] <-
                    sqrt( sum(RR[lower.tri(RR, diag=TRUE)]^2) / e )
                crmr_nomean.group[g] <- crmr.group[g] <-
                    sqrt(  sum(R.cor[lower.tri(R.cor, diag=TRUE)]^2) / (e - nvar) )
                srmr_mplus_nomean.group[g] <- srmr_mplus.group[g] <-
                    sqrt( (sum(R.cor[lower.tri(R.cor, diag=FALSE)]^2)  +
                           sum(((diag(S) - diag(Sigma.hat))/diag(S))^2)) / e )
            }
        }

        if(G > 1) {
            ## FIXME: get the scaling right
            SRMR_BENTLER <- as.numeric( (unlist(object@SampleStats@nobs) %*% srmr_bentler.group) / object@SampleStats@ntotal )
            SRMR_BENTLER_NOMEAN <- as.numeric( (unlist(object@SampleStats@nobs) %*% srmr_bentler_nomean.group) / object@SampleStats@ntotal )
            CRMR <- as.numeric( (unlist(object@SampleStats@nobs) %*% crmr.group) / object@SampleStats@ntotal )
            CRMR_NOMEAN <- as.numeric( (unlist(object@SampleStats@nobs) %*% crmr_nomean.group) / object@SampleStats@ntotal )
            SRMR_MPLUS <- as.numeric( (unlist(object@SampleStats@nobs) %*% srmr_mplus.group) / object@SampleStats@ntotal )
            SRMR_MPLUS_NOMEAN <- as.numeric( (unlist(object@SampleStats@nobs) %*% srmr_mplus_nomean.group) / object@SampleStats@ntotal )
            RMR <- as.numeric( (unlist(object@SampleStats@nobs) %*% rmr.group) / object@SampleStats@ntotal )
            RMR_NOMEAN <- as.numeric( (unlist(object@SampleStats@nobs) %*% rmr_nomean.group) / object@SampleStats@ntotal )
        } else {
            SRMR_BENTLER <- srmr_bentler.group[1]
            SRMR_BENTLER_NOMEAN <- srmr_bentler_nomean.group[1]
            CRMR <- crmr.group[1]
            CRMR_NOMEAN <- crmr_nomean.group[1]
            SRMR_MPLUS <- srmr_mplus.group[1]
            SRMR_MPLUS_NOMEAN <- srmr_mplus_nomean.group[1]
            RMR <- rmr.group[1]
            RMR_NOMEAN <- rmr_nomean.group[1]
        }

        # the default!
        if(object@Options$mimic %in% c("lavaan", "EQS")) {
            if(categorical) {
                indices["srmr"] <- SRMR_BENTLER_NOMEAN
            } else {
                indices["srmr"] <- SRMR_BENTLER
            }
        } else if(object@Options$mimic == "Mplus") {
            if(object@Options$information[1] == "expected") {
                if(categorical) {
                    indices["srmr"] <- SRMR_BENTLER_NOMEAN
                } else {
                    indices["srmr"] <- SRMR_BENTLER
                }
            } else {
                if(categorical) {
                    indices["srmr"] <- SRMR_MPLUS_NOMEAN
                } else {
                    indices["srmr"] <- SRMR_MPLUS
                }
            }
        }

        # the others
        indices["srmr_bentler"]        <- SRMR_BENTLER
        indices["srmr_bentler_nomean"] <- SRMR_BENTLER_NOMEAN
        indices["crmr"]                <- CRMR
        indices["crmr_nomean"]         <- CRMR_NOMEAN

        # only correct for non-categorical:
        if(object@Model@categorical) {
            # FIXME! Compute Mplus 8.1 way to compute SRMR in the
            #        categorical setting
            #        See 'SRMR in Mplus (2018)' document on Mplus website
            indices["srmr_mplus"]          <- as.numeric(NA)
            indices["srmr_mplus_nomean"]   <- as.numeric(NA)
        } else {
            indices["srmr_mplus"]          <- SRMR_MPLUS
            indices["srmr_mplus_nomean"]   <- SRMR_MPLUS_NOMEAN
        }
        #if(categorical) {
            indices["rmr"]             <- RMR
        #} else {
        #    indices["rmr"]             <- RMR_NOMEAN
        #}
        indices["rmr_nomean"]          <- RMR_NOMEAN
    }

    # multilevel version
    if(any(c("srmr_within", "srmr_between", "srmr") %in% fit.measures) &&
       object@Data@nlevels > 1L) {
        nlevels <-  object@Data@nlevels > 1L
        SRMR.within  <- numeric(G)
        SRMR.between <- numeric(G)
        for(g in 1:G) {

            b.within  <- (g - 1L) * nlevels + 1L
            b.between <- (g - 1L) * nlevels + 2L

            # observed
            S.within  <- object@h1$implied$cov[[  b.within  ]]
            M.within  <- object@h1$implied$mean[[ b.within  ]]
            S.between <- object@h1$implied$cov[[  b.between ]]
            M.between <- object@h1$implied$mean[[ b.between ]]

            # estimated
            implied <- lav_model_implied_cond2uncond(object@implied)
            Sigma.within  <- implied$cov[[  b.within  ]]
            Mu.within     <- implied$mean[[ b.within  ]]
            Sigma.between <- implied$cov[[  b.between ]]
            Mu.between    <- implied$mean[[ b.between ]]

            # force pd for between
            #    S.between <- lav_matrix_symmetric_force_pd(S.between)
            Sigma.between <- lav_matrix_symmetric_force_pd(Sigma.between)

            # Bollen approach: simply using cov2cor ('residual correlations')
            S.within.cor  <- cov2cor(S.within)
            S.between.cor <- cov2cor(S.between)
            Sigma.within.cor <- cov2cor(Sigma.within)
            if(all(diag(Sigma.between) > 0)) {
                Sigma.between.cor <- cov2cor(Sigma.between)
            } else {
                Sigma.between.cor <- matrix(as.numeric(NA),
                                         nrow = nrow(Sigma.between),
                                         ncol = ncol(Sigma.between))
            }
            R.within.cor <- (S.within.cor - Sigma.within.cor)
            R.between.cor <- (S.between.cor - Sigma.between.cor)

            nvar.within <- NCOL(S.within)
            nvar.between <- NCOL(S.between)
            pstar.within <- nvar.within*(nvar.within+1)/2
            pstar.between <- nvar.between*(nvar.between+1)/2

            # SRMR
            SRMR.within[g] <-  sqrt( sum(lav_matrix_vech(R.within.cor)^2) /
                                     pstar.within )
            SRMR.between[g] <- sqrt( sum(lav_matrix_vech(R.between.cor)^2) /
                                     pstar.between )
        }

        if(G > 1) {
            SRMR_WITHIN  <- as.numeric( (unlist(object@SampleStats@nobs) %*%
                              SRMR.within)  / object@SampleStats@ntotal )
            SRMR_BETWEEN <- as.numeric( (unlist(object@SampleStats@nobs) %*%
                              SRMR.between) / object@SampleStats@ntotal )
        } else {
            SRMR_WITHIN  <- SRMR.within[1]
            SRMR_BETWEEN <- SRMR.between[1]
        }

        indices["srmr"] <- SRMR_WITHIN + SRMR_BETWEEN
        indices["srmr_within"]  <- SRMR_WITHIN
        indices["srmr_between"] <- SRMR_BETWEEN
    }


    if(any(c("cn_05", "cn_01") %in% fit.measures)) {
        # catch df=0, X2=0
        if(df == 0 && X2 < .Machine$double.eps) {
            CN_05 <- as.numeric(NA)
            CN_01 <- as.numeric(NA)
        } else {
            CN_05 <- qchisq(p=0.95, df=df)/(X2/N) + 1
            CN_01 <- qchisq(p=0.99, df=df)/(X2/N) + 1
        }
        indices["cn_05"] <- CN_05
        indices["cn_01"] <- CN_01
    }

    if("wrmr" %in% fit.measures) {
        # we use the definition: wrmr = sqrt ( 2*N*F / e )
        e <- length(object@SampleStats@WLS.obs[[1]]) ### only first group???
        WRMR <- sqrt( X2 / e )
        indices["wrmr"] <- WRMR
    }

    if(any(c("gfi","agfi","pgfi") %in% fit.measures)) {
        gfi.group <- numeric(G)
        WLS.obs <- object@SampleStats@WLS.obs
        #WLS.V   <- lav_model_wls_v(lavmodel       = object@Model,
        #                           lavsamplestats = object@SampleStats,
        #                           structured     = TRUE,
        #                           lavdata        = object@Data)
        # WLS.V == h1 expected information
        #WLS.V   <- lav_model_h1_information(lavobject = object)
        WLS.V   <- lav_object_inspect_wls_v(object)
        WLS.est <- lav_object_inspect_wls_est(object)
        for(g in 1:G) {
            wls.obs <- WLS.obs[[g]]
            wls.est <- WLS.est[[g]]
            wls.v   <- WLS.V[[g]]

            if(is.null(wls.v)) {
                gfi.group[g] <- as.numeric(NA)
            } else {
                wls.diff <- wls.obs - wls.est
                if(is.matrix(wls.v)) {
                    # full weight matrix
                    t1 <- crossprod(wls.diff, wls.v) %*% wls.diff
                    t2 <- crossprod(wls.obs, wls.v) %*% wls.obs
                } else {
                    # diagonal weight matrix
                    t1 <- as.numeric(crossprod(wls.diff^2, wls.v))
                    t2 <- as.numeric(crossprod(wls.obs^2, wls.v))
                }
                gfi.group[g] <- 1 - t1/t2
            }
        }
        if(G > 1) {
            ## FIXME: get the scaling right
            GFI <- as.numeric( (unlist(object@SampleStats@nobs) %*% gfi.group) / object@SampleStats@ntotal )
        } else {
            GFI <- gfi.group[1L]
        }
        indices["gfi"] <- GFI
        nel <- length(unlist(WLS.obs)) # total number of modeled sample stats
        if(is.na(df)) {
            indices["agfi"] <- as.numeric(NA)
        } else if(df > 0) {
            indices["agfi"] <- 1 - (nel/df) * (1 - GFI)
        } else {
            indices["agfi"] <- 1
        }
        # LISREL formula (Simplis book 2002, p. 126)
        indices["pgfi"] <- (df/nel)*GFI
    }

    # MFI - McDonald Fit Index (McDonald, 1989)
    if("mfi" %in% fit.measures) {
        #MFI <- exp(-0.5 * (X2 - df)/(N-1)) # Hu & Bentler 1998 Table 1
        MFI <- exp(-0.5 * (X2 - df)/N)
        indices["mfi"] <- MFI
    }

    # ECVI - cross-validation index (Brown & Cudeck, 1989)
    # not defined for multiple groups and/or models with meanstructures
    # TDJ: According to Dudgeon (2004, p. 317), "ECVI requires no adjustment
    #      when a model is fitted simultaneously in multiple samples."
    #      And I think the lack of mean structure in Brown & Cudeck (1989)
    #      was a matter of habitual simplification back then, not necessity.
    if("ecvi" %in% fit.measures) {
        ECVI <- X2/N + (2*npar)/N
        indices["ecvi"] <- ECVI
    }


    # C_p
    if("c_p" %in% fit.measures) {
        out <- lavTablesFitCp(object)
        CpMax <- attr(out, "CpMax")
        indices["c_p"] <- CpMax$LR
        indices["c_p.df"] <- CpMax$df
        indices["c_p.p.value"] <- CpMax$p.value.Bonferroni
    }

    # C_F
    if("c_f" %in% fit.measures) {
        CF <- lavTablesFitCf(object)
        DF <- attr(CF, "DF")
        attributes(CF) <- NULL
        indices["c_f"] <- CF
        indices["c_f.df"] <- DF
        indices["c_f.p.value"] <- pchisq(CF, DF, lower.tail=FALSE)
        indices["rpat.observed"] <- object@Data@Rp[[1L]]$npatterns
        indices["rpat.total"] <- object@Data@Rp[[1L]]$total.patterns
        indices["rpat.empty"] <- object@Data@Rp[[1L]]$empty.patterns
    }


    # C_M
    if("c_m" %in% fit.measures) {
        CM <- lavTablesFitCm(object)
        DF <- attr(CM, "DF")
        attributes(CM) <- NULL
        indices["c_m"] <- CM
        indices["c_m.df"] <- DF
        indices["c_m.p.value"] <- pchisq(CM, DF, lower.tail=FALSE)
    }

    if("ntotal" %in% fit.measures) {
        indices["ntotal"] <- object@SampleStats@ntotal
    }

    # do we have everything that we requested?
    #idx.missing <- which(is.na(match(fit.measures, names(indices))))
    #if(length(idx.missing) > 0L) {
    #    cat("lavaan WARNING: some requested fit measure(s) are not available for this model:\n")
    #    print( fit.measures[ idx.missing ] )
    #    cat("\n")
    #}

    out <- unlist(indices[fit.measures])

    if(length(out) > 0L) {
        if(output == "vector") {
            class(out) <- c("lavaan.vector", "numeric")
        } else if(output == "matrix") {
            out <- as.matrix(out)
            colnames(out) <- ""
            class(out) <- c("lavaan.matrix", "matrix")
        } else if(output == "text") {
            class(out) <- c("lavaan.fitMeasures", "lavaan.vector", "numeric")
        }
    } else {
        return( invisible(numeric(0)) )
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
            c1 <- c(c1, "Sample-size adjusted Bayesian (BIC)")
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
            c1 <- c(c1, "P-value RMSEA <= 0.05")
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

# new in 0.6-5
# internal function to check the (external) baseline model, and
# return baseline 'test' list if everything checks out (and NULL otherwise)
lav_fit_measures_check_baseline <- function(fit.indep = NULL, object = NULL) {

    TEST <- NULL

    # check if everything is in order
    if( inherits(fit.indep, "try-error") ) {
        warning("lavaan WARNING: baseline model estimation failed")
        return(NULL)

    } else if( !inherits(fit.indep, "lavaan") ) {
        warning("lavaan WARNING: (user-provided) baseline model ",
                "is not a fitted lavaan object")
        return(NULL)

    } else if( !fit.indep@optim$converged ) {
        warning("lavaan WARNING: baseline model did not converge")
        return(NULL)

    } else {

        # evaluate if estimator/test matches original object
        # note: we do not need to check for 'se', as it may be 'none'
        sameTest <- all(object@Options$test == fit.indep@Options$test)
        if(!sameTest) {
            warning("lavaan WARNING:\n",
                    "\t Baseline model was using test(s) = ",
                    dQuote(fit.indep@Options$test),
                    "\n\t But original model was using test(s) = ",
                    dQuote(object@Options$test),
                    "\n\t Refitting baseline model!")
        }
        sameEstimator <- ( object@Options$estimator ==
                           fit.indep@Options$estimator )
        if(!sameEstimator) {
            warning("lavaan WARNING:\n",
                    "\t Baseline model was using estimator = ",
                    dQuote(fit.indep@Options$estimator),
                    "\n\t But original model was using estimator = ",
                    dQuote(object@Options$estimator),
                    "\n\t Refitting baseline model!")
        }
        if( !sameTest || !sameEstimator ) {
            lavoptions <- object@Options
            lavoptions$estimator   <- object@Options$estimator
            lavoptions$se          <- "none"
            lavoptions$verbose     <- FALSE
            lavoptions$baseline    <- FALSE
            lavoptions$check.start <- FALSE
            lavoptions$check.post  <- FALSE
            lavoptions$check.vcov  <- FALSE
            lavoptions$test        <- object@Options$test
            fit.indep <- try(lavaan(fit.indep,
                                    slotOptions     = lavoptions,
                                    slotData        = object@Data,
                                    slotSampleStats = object@SampleStats,
                                    sloth1          = object@h1,
                                    slotCache       = object@Cache),
                             silent = TRUE)
            # try again
            TEST <- lav_fit_measures_check_baseline(fit.indep = fit.indep,
                                                    object    = object)

        } else {
            # extract what we need
            TEST <- fit.indep@test
        }

    } # converged lavaan object

    TEST
}


