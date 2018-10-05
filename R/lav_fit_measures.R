# TDJ: add "..." to make the method generic, so lavaan.mi can add arguments to
#      pass to lavTestLRT() and lavTestLRT.mi() about how to pool chi-squared.
#      NOT sure this is necessary for the lavaan-method, perhaps only the
#      generic needs "..."?
setMethod("fitMeasures", signature(object = "lavaan"),
function(object, fit.measures = "all", baseline.model = NULL, ...) {
    lav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model)
})

setMethod("fitmeasures", signature(object = "lavaan"),
function(object, fit.measures = "all", baseline.model = NULL, ...) {
    lav_fit_measures(object = object, fit.measures = fit.measures,
                     baseline.model = baseline.model)
})

lav_fit_measures <- function(object, fit.measures="all",
                             baseline.model = NULL) {

    # has the model converged?
    if(object@optim$npar > 0L && !object@optim$converged) {
        stop("lavaan ERROR: fit measures not available if model did not converge")
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
    npar <- object@optim$npar
    if(nrow(object@Model@con.jac) > 0L) {
        ceq.idx <- attr(object@Model@con.jac, "ceq.idx")
        if(length(ceq.idx) > 0L) {
            neq <- qr(object@Model@con.jac[ceq.idx,,drop=FALSE])$rank
            npar <- npar - neq
        }
    }

    fx <- object@optim$fx
    fx.group <- object@optim$fx.group
    meanstructure <- object@Model@meanstructure
    categorical   <- object@Model@categorical
    multigroup    <- object@Data@ngroups > 1L
    estimator     <- object@Options$estimator
    test          <- object@Options$test
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
    if(scaled && object@Options$test %in%
                 c("yuan.bentler", "yuan.bentler.mplus")) {
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

        # call explicitly independence model
        # this is not strictly needed for ML, but it is for
        # GLS and WLS
        # and MLM and MLR to get the scaling factor(s)!
        if (!is.null(baseline.model) && inherits(baseline.model, "lavaan")) {
            fit.indep <- baseline.model
        } else if (!is.null(object@external$baseline.model) &&
                   inherits(object@external$baseline.model, "lavaan")) {
            fit.indep <- object@external$baseline.model
            ## check baseline converged
            if (!fit.indep@optim$converged) {
                fit.indep <- NULL
            } else {
                ## check test matches and baseline converged
                sameTest <- ( object@Options$test == fit.indep@Options$test )
                sameSE <- ( object@Options$se == fit.indep@Options$se )
                sameEstimator <- ( object@Options$estimator == fit.indep@Options$estimator )
                if (!all(sameTest, sameSE, sameEstimator)) {
                    fit.indep <- try(update(fit.indep,
                                            test = object@Options$test,
                                            se   = object@Options$se,
                                            estimator = object@Options$estimator),
                                     silent = TRUE)
                }
            }
        } else {
            fit.indep <- try(lav_object_independence(object), silent = TRUE)
        }

        if(inherits(fit.indep, "try-error")) {
            X2.null <- df.null <- as.numeric(NA)
            if(scaled) {
                X2.null.scaled <- df.null.scaled <- as.numeric(NA)
            }
        } else {
            if(fit.indep@Data@nlevels > 1L) {
                fit.indep@test[[1]]$stat <- ( -2 * (fit.indep@loglik$loglik -
                                                    object@loglik$loglik) )
                fit.indep@test[[1]]$pvalue <-
                   1 - pchisq(fit.indep@test[[1]]$stat, fit.indep@test[[1]]$df)
            }

            X2.null <- fit.indep@test[[1]]$stat
            df.null <- fit.indep@test[[1]]$df
            if(scaled) {
                X2.null.scaled <- fit.indep@test[[2]]$stat
                df.null.scaled <- fit.indep@test[[2]]$df
            }
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
                indices["baseline.pvalue"] <- fit.indep@test[[1]]$pvalue
                if(scaled) {
                    indices["baseline.pvalue.scaled"] <-
                        fit.indep@test[[2]]$pvalue
                }
            }
            if("baseline.chisq.scaling.factor" %in% fit.measures) {
                indices["baseline.chisq.scaling.factor"] <-
                    fit.indep@test[[2]]$scaling.factor
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
                    cb <- fit.indep@test[[2]]$scaling.factor

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
                    cb <- fit.indep@test[[2]]$scaling.factor

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
                    cb <- fit.indep@test[[2]]$scaling.factor

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
            if(object@Options$test %in%
               c("yuan.bentler", "yuan.bentler.mplus")) {
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
        } else if(df < 1 || df2 < 1 || lower.lambda(0) < 0.0) {
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
        } else if(df < 1 || df2 < 1 || upper.lambda(N.RMSEA) > 0 ||
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
        if(object@Options$mimic == "lavaan") {
            if(categorical) {
                indices["srmr"] <- SRMR_BENTLER_NOMEAN
            } else {
                indices["srmr"] <- SRMR_BENTLER
            }
        } else if(object@Options$mimic == "EQS") {
            indices["srmr"] <- SRMR_BENTLER
        } else if(object@Options$mimic == "Mplus") {
            if(object@Options$information == "expected") {
                indices["srmr"] <- SRMR_BENTLER
            } else {
                indices["srmr"] <- SRMR_MPLUS
            }
        }

        # the others
        indices["srmr_bentler"]        <- SRMR_BENTLER
        indices["srmr_bentler_nomean"] <- SRMR_BENTLER_NOMEAN
        indices["crmr"]                <- CRMR
        indices["crmr_nomean"]         <- CRMR_NOMEAN
        indices["srmr_mplus"]          <- SRMR_MPLUS
        indices["srmr_mplus_nomean"]   <- SRMR_MPLUS_NOMEAN
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
            # observed
            S.within  <- object@h1$implied$cov[[  (g-1)*nlevels + 1 ]]
            M.within  <- object@h1$implied$mean[[ (g-1)*nlevels + 1 ]]
            S.between <- object@h1$implied$cov[[  (g-1)*nlevels + 2 ]]
            M.between <- object@h1$implied$mean[[ (g-1)*nlevels + 2 ]]

            # estimated
            Sigma.within  <- object@implied$cov[[  (g-1)*nlevels + 1 ]]
            Mu.within     <- object@implied$mean[[ (g-1)*nlevels + 1 ]]
            Sigma.between <- object@implied$cov[[  (g-1)*nlevels + 2 ]]
            Mu.between    <- object@implied$mean[[ (g-1)*nlevels + 2 ]]

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
        if(df == 0 && X2 == 0) {
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
        class(out) <- c("lavaan.vector", "numeric")
    } else {
        return( invisible(numeric(0)) )
    }

    out
}

# print a nice summary of the fit measures
print.fit.measures <- function(x) {

   names.x <- names(x)

   # scaled?
   scaled <- "chisq.scaled" %in% names.x

   # table fit measures
   if("C_F" %in% names.x) {
       cat("\nFull response patterns fit statistics:\n\n")

       t0.txt <- sprintf("  %-40s", "Observed response patterns (1st group):")
       t1.txt <- sprintf("  %10i", x["rpat.observed"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Total response patterns (1st group):")
       t1.txt <- sprintf("  %10i", x["rpat.total"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Empty response patterns (1st group):")
       t1.txt <- sprintf("  %10i", x["rpat.empty"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       cat("\n")

       t0.txt <- sprintf("  %-40s", "C_F Test Statistic")
       t1.txt <- sprintf("  %10.3f", x["C_F"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Degrees of freedom")
       t1.txt <- sprintf("  %10i", x["C_F.df"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "P-value")
       t1.txt <- sprintf("  %10.3f", x["C_F.p.value"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       cat("\n")
       t0.txt <- sprintf("  %-40s", "C_M Test Statistic")
       t1.txt <- sprintf("  %10.3f", x["C_M"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Degrees of freedom")
       t1.txt <- sprintf("  %10i", x["C_M.df"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "P-value")
       t1.txt <- sprintf("  %10.3f", x["C_M.p.value"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
   }

   if("C_p" %in% names.x) {
       cat("\nPairwise tables summary statistic:\n\n")
       t0.txt <- sprintf("  %-40s", "C_P Test Statistic")
       t1.txt <- sprintf("  %10.3f", x["C_p"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Degrees of freedom")
       t1.txt <- sprintf("  %10i", x["C_p.df"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Bonferroni corrected P-value")
       t1.txt <- sprintf("  %10.3f", x["C_p.p.value"])
       t2.txt <- ""
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
   }

   # independence model
   if("baseline.chisq" %in% names.x) {
       cat("\nModel test baseline model:\n\n")
       t0.txt <- sprintf("  %-40s", "Minimum Function Test Statistic")
       t1.txt <- sprintf("  %10.3f", x["baseline.chisq"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["baseline.chisq.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

       t0.txt <- sprintf("  %-40s", "Degrees of freedom")
       t1.txt <- sprintf("  %10i", x["baseline.df"])
       t2.txt <- ifelse(scaled,
                        ifelse(round(x["baseline.df.scaled"]) ==
                               x["baseline.df.scaled"],
                        sprintf("  %10i",   x["baseline.df.scaled"]),
                        sprintf("  %10.3f", x["baseline.df.scaled"])),
                        "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

       t0.txt <- sprintf("  %-40s", "P-value")
       t1.txt <- sprintf("  %10.3f", x["baseline.pvalue"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["baseline.pvalue.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    }

    # cfi/tli
    if(any(c("cfi","tli","nnfi","rfi","nfi","ifi","rni","pnfi") %in% names.x)) {
        cat("\nUser model versus baseline model:\n\n")

        if("cfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Comparative Fit Index (CFI)")
            t1.txt <- sprintf("  %10.3f", x["cfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cfi.scaled"]), "")
            #t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cfi.robust"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("tli" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Tucker-Lewis Index (TLI)")
            t1.txt <- sprintf("  %10.3f", x["tli"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["tli.scaled"]), "")
            #t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["tli.robust"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("cfi.robust" %in% names.x) {
            t0.txt <- sprintf("\n  %-40s", "Robust Comparative Fit Index (CFI)")
            t1.txt <- sprintf("  %10s", "")
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cfi.robust"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("tli.robust" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Robust Tucker-Lewis Index (TLI)")
            t1.txt <- sprintf("  %10s", "")
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["tli.robust"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("nnfi" %in% names.x) {
            t0.txt <- sprintf("  %-42s", "Bentler-Bonett Non-normed Fit Index (NNFI)")
            t1.txt <- sprintf("  %8.3f", x["nnfi"])
            #t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["nnfi.scaled"]), "")
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["nnfi.robust"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("nfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Bentler-Bonett Normed Fit Index (NFI)")
            t1.txt <- sprintf("  %10.3f", x["nfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["nfi.scaled"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("nfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Parsimony Normed Fit Index (PNFI)")
            t1.txt <- sprintf("  %10.3f", x["pnfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["pnfi.scaled"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("rfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Bollen's Relative Fit Index (RFI)")
            t1.txt <- sprintf("  %10.3f", x["rfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rfi.scaled"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("ifi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Bollen's Incremental Fit Index (IFI)")
            t1.txt <- sprintf("  %10.3f", x["ifi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["ifi.scaled"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

        if("rni" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Relative Noncentrality Index (RNI)")
            t1.txt <- sprintf("  %10.3f", x["rni"])
            #t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rni.scaled"]), "")
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rni.robust"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
    }

   # likelihood
   if("logl" %in% names.x) {
       cat("\nLoglikelihood and Information Criteria:\n\n")
       t0.txt <- sprintf("  %-40s", "Loglikelihood user model (H0)")
       t1.txt <- sprintf("  %10.3f", x["logl"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["logl"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
       if(!is.na(x["scaling.factor.h0"])) {
           t0.txt <- sprintf("  %-40s", "Scaling correction factor")
           t1.txt <- sprintf("  %10s", "")
           t2.txt <- sprintf("  %10.3f", x["scaling.factor.h0"])
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
           cat("    for the MLR correction\n")
       }

       if("unrestricted.logl" %in% names.x) {
           t0.txt <- sprintf("  %-40s", "Loglikelihood unrestricted model (H1)")
           t1.txt <- sprintf("  %10.3f", x["unrestricted.logl"])
           t2.txt <- ifelse(scaled,
                     sprintf("  %10.3f", x["unrestricted.logl"]), "")
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
           #cat(t0.txt, t1.txt, "\n", sep="")
           if(!is.na(x["scaling.factor.h1"])) {
               t0.txt <- sprintf("  %-40s", "Scaling correction factor")
               t1.txt <- sprintf("  %10s", "")
               t2.txt <- sprintf("  %10.3f", x["scaling.factor.h1"])
               cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
              cat("    for the MLR correction\n")
           }
       }

       cat("\n")
       t0.txt <- sprintf("  %-40s", "Number of free parameters")
       t1.txt <- sprintf("  %10i", x["npar"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10i", x["npar"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

       t0.txt <- sprintf("  %-40s", "Akaike (AIC)")
       t1.txt <- sprintf("  %10.3f", x["aic"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["aic"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Bayesian (BIC)")
       t1.txt <- sprintf("  %10.3f", x["bic"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["bic"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
       if(!is.na(x["bic2"])) {
           t0.txt <- sprintf("  %-40s", "Sample-size adjusted Bayesian (BIC)")
           t1.txt <- sprintf("  %10.3f", x["bic2"])
           t2.txt <- ifelse(scaled,
                     sprintf("  %10.3f", x["bic2"]), "")
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
   }

   # RMSEA
   if("rmsea" %in% names.x) {
       cat("\nRoot Mean Square Error of Approximation:\n\n")
       t0.txt <- sprintf("  %-40s", "RMSEA")
       t1.txt <- sprintf("  %10.3f", x["rmsea"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["rmsea.scaled"]), "")
                 #sprintf("  %10.3f", x["rmsea.robust"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       if("rmsea.ci.lower" %in% names.x) {
           t0.txt <- sprintf("  %-38s", "90 Percent Confidence Interval")
           t1.txt <- sprintf("  %5.3f", x["rmsea.ci.lower"])
           t2.txt <- sprintf("  %5.3f", x["rmsea.ci.upper"])
           t3.txt <- ifelse(scaled,
                     sprintf("       %5.3f  %5.3f", x["rmsea.ci.lower.scaled"],
                                               x["rmsea.ci.upper.scaled"]), "")
                     #sprintf("       %5.3f  %5.3f", x["rmsea.ci.lower.robust"],
                     #                          x["rmsea.ci.upper.robust"]), "")
           cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n", sep="")
       }
       if("rmsea.pvalue" %in% names.x) {
           t0.txt <- sprintf("  %-40s", "P-value RMSEA <= 0.05")
           t1.txt <- sprintf("  %10.3f", x["rmsea.pvalue"])
           t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["rmsea.pvalue.scaled"]), "")
                 #sprintf("  %10.3f", x["rmsea.pvalue.robust"]), "")
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       }

       # robust
       if("rmsea.robust" %in% names.x) {
           t0.txt <- sprintf("\n  %-40s", "Robust RMSEA")
           t1.txt <- sprintf("  %10s", "")
           t2.txt <- ifelse(scaled,
                     #sprintf("  %10.3f", x["rmsea.scaled"]), "")
                     sprintf("  %10.3f", x["rmsea.robust"]), "")
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       }
       if("rmsea.ci.lower.robust" %in% names.x) {
           t0.txt <- sprintf("  %-38s", "90 Percent Confidence Interval")
           t1.txt <- sprintf("  %5.3s", "")
           t2.txt <- sprintf("  %5.3s", "")
           t3.txt <- ifelse(scaled,
                     #sprintf("       %5.3f  %5.3f", x["rmsea.ci.lower.scaled"],
                     #                          x["rmsea.ci.upper.scaled"]), "")
                     sprintf("       %5.3f  %5.3f", x["rmsea.ci.lower.robust"],
                                               x["rmsea.ci.upper.robust"]), "")
           cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n", sep="")
       }
       #if("rmsea.pvalue.robust" %in% names.x) {
       #    t0.txt <- sprintf("  %-40s", "P-value RMSEA <= 0.05")
       #    t1.txt <- sprintf("  %10.3f", x["rmsea.pvalue"])
       #    t2.txt <- ifelse(scaled,
       #          #sprintf("  %10.3f", x["rmsea.pvalue.scaled"]), "")
       #          sprintf("  %10.3f", x["rmsea.pvalue.robust"]), "")
       #    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #}
   }

    # SRMR
    if(any(c("rmr","srmr") %in% names.x) && ! "srmr_within" %in% names.x) {
        cat("\nStandardized Root Mean Square Residual:\n\n")

        if("rmr" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "RMR")
            t1.txt <- sprintf("  %10.3f", x["rmr"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rmr"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("rmr_nomean" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "RMR (No Mean)")
            t1.txt <- sprintf("  %10.3f", x["rmr_nomean"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["rmr_nomean"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("srmr" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "SRMR")
            t1.txt <- sprintf("  %10.3f", x["srmr"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["srmr"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("srmr_nomean" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "SRMR (No Mean)")
            t1.txt <- sprintf("  %10.3f", x["srmr_nomean"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["srmr_nomean"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
    }

    # SRMR -- multilevel
    if(any(c("srmr_within","srmr_between") %in% names.x)) {
        cat("\nStandardized Root Mean Square Residual (corr metric):\n\n")
        if("srmr_within" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "SRMR (within covariance matrix)")
            t1.txt <- sprintf("  %10.3f", x["srmr_within"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["srmr_within"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("srmr_between" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "SRMR (between covariance matrix)")
            t1.txt <- sprintf("  %10.3f", x["srmr_between"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["srmr_between"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
    }

    # WRMR
    if("wrmr" %in% names.x) {
        cat("\nWeighted Root Mean Square Residual:\n\n")

        if("wrmr" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "WRMR")
            t1.txt <- sprintf("  %10.3f", x["wrmr"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["wrmr"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
    }

    # Other
    if(any(c("cn_05","cn_01","gfi","agfi","pgfi","mfi") %in% names.x)) {
        cat("\nOther Fit Indices:\n\n")

        if("cn_05" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Hoelter Critical N (CN) alpha=0.05")
            t1.txt <- sprintf("  %10.3f", x["cn_05"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cn_05"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("cn_01" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Hoelter Critical N (CN) alpha=0.01")
            t1.txt <- sprintf("  %10.3f", x["cn_01"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["cn_01"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if(any(c("cn_05", "cn_01") %in% names.x)) {
            cat("\n")
        }
        if("gfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Goodness of Fit Index (GFI)")
            t1.txt <- sprintf("  %10.3f", x["gfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["gfi"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("agfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Adjusted Goodness of Fit Index (AGFI)")
            t1.txt <- sprintf("  %10.3f", x["agfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["agfi"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("pgfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Parsimony Goodness of Fit Index (PGFI)")
            t1.txt <- sprintf("  %10.3f", x["pgfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["pgfi"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if(any(c("gfi","agfi","pgfi") %in% names.x)) {
            cat("\n")
        }
        if("mfi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "McDonald Fit Index (MFI)")
            t1.txt <- sprintf("  %10.3f", x["mfi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["mfi"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
        if("mfi" %in% names.x) {
            cat("\n")
        }
        if("ecvi" %in% names.x) {
            t0.txt <- sprintf("  %-40s", "Expected Cross-Validation Index (ECVI)")
            t1.txt <- sprintf("  %10.3f", x["ecvi"])
            t2.txt <- ifelse(scaled, sprintf("  %10.3f", x["ecvi"]), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }

    }

    #cat("\n")
}


