# functions related to the RMSEA index of approximate fit

# lower-level functions: no checking of input: just compute the number(s):
# - lav_fit_rmsea
# - lav_fit_rmsea_ci
# - lav_fit_rmsea_closefit
# - lav_fit_rmsea_notclosefit (TODO)

# higher-level functions:
# - lav_fit_rmsea_lavobject

# Y.R. 19 July 2022

# we assume X2 = N * F.val
# lambda = (X2 - df) is the non-centrality parameter

# RMSEA: sqrt( (X2 - df)/(N * df) )
#      = sqrt( lambda/(N * df) )
#      = sqrt( ((N*F.val) - df)/(N * df) )
#      = sqrt( (N.F.val)/(N * df) - df/(N * df) )
#      = sqrt( F.val/df - 1/N )
#      = sqrt( (X2/N)/df - 1/N )


# 'scaled' RMSEA: X2 is replaced by X2-SB (or any other 'scaled' statistic)

# robust RMSEA: sqrt( (X2/N)/df - c.hat/N )
# note:
# - robust RMSEA == scaled RMSEA * sqrt(c.hat)
# - robust RMSEA CI == scaled RMSEA CI * sqrt(c.hat)

# robust RMSEA for MLMV (ie. scaled.shifted):
# - c == a * (df - b) / df
# - robust RMSEA MLM == robust RMSEA MLMV
# - robust RMSEA MLMV == scaled RMSEA MLMV * sqrt(a)
# - robust RMSEA CI MLMV == scaled RMSEA CI MLMV * sqrt(a)

# References:
# Steiger, J. H., & Lind, J. C. (1980, May). Statistically based tests for the
# number of common factors.  Paper presented at the annual meeting of the
# Psychometric Society, Iowa City, IA.

# confidence interval:
# Browne, M. W., & Cudeck, R. (1993). Alternative ways of assessing model fit.
# In K. A. Bollen & J. S. Long (Eds.), Testing structural equation models (pp.
# 136-162). Newbury Park, CA: Sage.

# problems with low df
# Kenny, D. A., Kaniskan, B., & McCoach, D. B. (2015).  The performance of RMSEA
# in models with small degrees of freedom.  Sociological Methods & Research, 44,
# 486-507.

# robust version MLM
# Patricia E. Brosseau-Liard , Victoria Savalei & Libo Li (2012) An
# Investigation of the Sample Performance of Two Nonnormality Corrections for
# RMSEA, Multivariate Behavioral Research, 47:6, 904-930, DOI:
# 10.1080/00273171.2012.715252

# robust version MLMV (scaled.shifted)
# Savalei, V. (2018). On the computation of the RMSEA and CFI from the
# mean-and-variance corrected test statistic with nonnormal data in SEM.
# Multivariate behavioral research, 53(3), 419-429.


# always using N (if a user needs N-1, just replace N by N-1)
# vectorized!
lav_fit_rmsea <- function(X2 = NULL, df = NULL, N = NULL,
                          F.val = NULL, G = 1L, c.hat = 1.0) {

    # did we get a sample size?
    if(missing(N) && !missing(F.val)) {
        # population version
        RMSEA <- sqrt( F.val/df )
    } else {
        nel <- length(X2)
        if(nel == 0) {
            return(as.numeric(NA))
        }
        RMSEA <- ifelse(df > 0,
               # 'standard' way to compute RMSEA
               RMSEA <- sqrt(pmax((X2/N)/df - c.hat/N, rep(0, nel))) * sqrt(G),
               0) # if df == 0
    }

    RMSEA
}

# note: for 'robust' version, X2 should be SB-X2
lav_fit_rmsea_ci <- function(X2 = NULL, df = NULL, N = NULL,
                             G = 1L, c.hat = 1, level = 0.90) {

    if(missing(N) || missing(X2) || missing(df) ||
       !is.finite(X2) || !is.finite(df) || !is.finite(N)) {
        return(list(rmsea.ci.lower = as.numeric(NA),
                    rmsea.ci.upper = as.numeric(NA)))
    }

    if(!is.finite(level) || level < 0 || level > 1.0) {
        warning("lavaan WARNING: invalid level value [",
                        level, "] set to default 0.90.")
        level <- 0.90
    }

    upper.perc <- (1 - (1 - level)/2)
    lower.perc <- (1 - level)/2

    # internal function
    lower.lambda <- function(lambda) {
        (pchisq(X2, df = df, ncp = lambda) - upper.perc)
    }

    upper.lambda <- function(lambda) {
        (pchisq(X2, df = df, ncp = lambda) - lower.perc)
    }


    # lower bound
    if(df < 1 || lower.lambda(0) < 0.0) {
        rmsea.ci.lower <- 0
    } else {
        lambda.l <- try(uniroot(f = lower.lambda, lower = 0, upper = X2)$root,
                        silent = TRUE)
        if(inherits(lambda.l, "try-error")) {
            lambda.l <- as.numeric(NA)
        }
        # lower bound
        rmsea.ci.lower <- sqrt( (c.hat * lambda.l)/(N * df) )

        # multiple groups? -> correction
        if(G > 1L) {
            rmsea.ci.lower <- rmsea.ci.lower * sqrt(G)
        }
    }

    # upper bound
    N.RMSEA <- max(N, X2 * 4)
    if(df < 1 || upper.lambda(N.RMSEA) > 0 || upper.lambda(0) < 0) {
        rmsea.ci.upper <- 0
    } else {
        lambda.u <- try(uniroot(f = upper.lambda, lower = 0,
                                                  upper = N.RMSEA)$root,
                        silent = TRUE)
        if(inherits(lambda.u, "try-error")) {
            lambda.u <- NA
        }
        # upper bound
        rmsea.ci.upper <- sqrt( (c.hat * lambda.u)/(N * df) )

        # multiple groups? -> correction
        if(G > 1L) {
            rmsea.ci.upper <- rmsea.ci.upper * sqrt(G)
        }
    }

    list(rmsea.ci.lower = rmsea.ci.lower,
         rmsea.ci.upper = rmsea.ci.upper)
}

# H_0: RMSEA <= rmsea.h0
lav_fit_rmsea_closefit <- function(X2 = NULL, df = NULL, N = NULL,
                                   G = 1L, c.hat = 1, rmsea.h0 = 0.05) {

    if(missing(N) || missing(X2) || missing(df) ||
       !is.finite(X2) || !is.finite(df) || !is.finite(N)) {
        return(as.numeric(NA))
    }

    rmsea.pvalue <- as.numeric(NA)
    if(df > 0) {
        # see Dudgeon 2004, eq 16 for the 'G' correction
        ncp <- (N * df * 1/c.hat * rmsea.h0^2) / G
        rmsea.pvalue <- 1 - pchisq(X2, df = df, ncp = ncp)
    }

    rmsea.pvalue
}

# MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996).
# H_0: RMSEA >= rmsea.h0
lav_fit_rmsea_notclosefit <- function(X2 = NULL, df = NULL, N = NULL,
                                   G = 1L, c.hat = 1, rmsea.h0 = 0.05) {

    if(missing(N) || missing(X2) || missing(df) ||
       !is.finite(X2) || !is.finite(df) || !is.finite(N)) {
        return(as.numeric(NA))
    }

    rmsea.pvalue <- as.numeric(NA)
    if(df > 0) {
        # see Dudgeon 2004, eq 16 for the 'G' correction
        ncp <- (N * df * 1/c.hat * rmsea.h0^2) / G
        rmsea.pvalue <- pchisq(X2, df = df, ncp = ncp)
    }

    rmsea.pvalue
}


lav_fit_rmsea_lavobject <- function(lavobject = NULL, fit.measures = "rmsea",
                                    standard.test = "standard",
                                    scaled.test   = "none",
                                    ci.level = 0.90,
                                    close.h0 = 0.05, notclose.h0 = 0.08) {

    # check lavobject
    stopifnot(inherits(lavobject, "lavaan"))

    # check for categorical
    categorical <- lavobject@Model@categorical
    if(lavobject@Model@categorical) {
        # 'refit' using estimator = "catML"
        fit.catml <- lav_object_catml(lavobject)
    }

    # tests
    TEST <- lavobject@test
    test.names <- sapply(lavobject@test, "[[", "test")
    if(test.names[1] == "none" || standard.test == "none") {
        return(list())
    }
    test.idx <- which(test.names == standard.test)
    if(length(test.idx) == 0L) {
        return(list())
    }

    scaled.flag <- FALSE
    if(!scaled.test %in% c("none", "standard", "default")) {
        scaled.idx <- which(test.names == scaled.test)
        if(length(scaled.idx) > 0L) {
            scaled.idx <- scaled.idx[1] # only the first one
            scaled.flag <- TRUE
        }
    }

    # supported fit measures in this function
    fit.rmsea <- c("rmsea",
                   "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.ci.level",
                   "rmsea.pvalue", "rmsea.close.h0",
                   "rmsea.notclose.pvalue", "rmsea.notclose.h0")

    if(scaled.flag) {
        fit.rmsea <- c(fit.rmsea, "rmsea.scaled", "rmsea.ci.lower.scaled",
                       "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled",
                       "rmsea.notclose.pvalue.scaled",
                       "rmsea.robust", "rmsea.ci.lower.robust",
                       "rmsea.ci.upper.robust", "rmsea.pvalue.robust",
                       "rmsea.notclose.pvalue.robust")
    }

    # which one do we need?
    if(missing(fit.measures)) {
        # default set
        fit.measures <- fit.rmsea
    } else {
        # remove any not-RMSEA related index from fit.measures
        rm.idx <- which(!fit.measures %in% fit.rmsea)
        if(length(rm.idx) > 0L) {
            fit.measures <- fit.measures[-rm.idx]
        }
        if(length(fit.measures) == 0L) {
            return(list())
        }
    }

    # robust?
    robust.flag <- FALSE
    if(scaled.flag &&
       scaled.test %in% c("satorra.bentler", "yuan.bentler.mplus",
                          "yuan.bentler", "scaled.shifted")) {
        robust.flag <- TRUE
    }

    # basic test statistics
    X2 <- TEST[[test.idx]]$stat
    df <- TEST[[test.idx]]$df
    G <- lavobject@Data@ngroups  # number of groups
    N <- lav_utils_get_ntotal(lavobject = lavobject) # N vs N-1

    # scaled X2/df values
    if(scaled.flag) {
        if(scaled.test == "scaled.shifted") {
            XX2 <- TEST[[scaled.idx]]$stat
            df2 <- df
        } else {
            XX2 <- X2
            df2 <- sum(TEST[[scaled.idx]]$trace.UGamma)
            if(!is.finite(df2) || df2 == 0) {
                df2 <- as.numeric(NA)
            }
        }
        if(robust.flag) {
            if(categorical) {
                XX3 <- fit.catml@test[[1]]$stat
                df3 <- fit.catml@test[[1]]$df
                c.hat3 <- c.hat <- fit.catml@test[[2]]$scaling.factor
                XX3.scaled <- fit.catml@test[[2]]$stat
            } else {
                XX3 <- X2
                df3 <- df
                c.hat <- TEST[[scaled.idx]]$scaling.factor
                if(scaled.test == "scaled.shifted") {
                   # compute c.hat from a and b
                    a <- TEST[[scaled.idx]]$scaling.factor
                    b <- TEST[[scaled.idx]]$shift.parameter
                    c.hat3 <- a * (df - b) / df
                } else {
                    c.hat3 <- c.hat
                }
                XX3.scaled <- TEST[[scaled.idx]]$stat
            }
        }
    }

    # output container
    indices <- list()

    # what do we need?
    rmsea.val.flag <- rmsea.ci.flag <- rmsea.pvalue.flag <- FALSE
    if(any(c("rmsea","rmsea.scaled","rmsea.robust") %in% fit.measures)) {
        rmsea.val.flag <- TRUE
    }
    if(any(c("rmsea.ci.lower", "rmsea.ci.upper", "rmsea.ci.level",
             "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled",
             "rmsea.ci.lower.robust", "rmsea.ci.upper.robust")
       %in% fit.measures)) {
        rmsea.ci.flag <- TRUE
    }
    if(any(c("rmsea.pvalue", "rmsea.pvalue.scaled", "rmsea.pvalue.robust",
             "rmsea.notclose.pvalue", "rmsea.notclose.pvalue.scaled",
             "rmsea.notclose.pvalue.robust",
             "rmsea.close.h0", "rmsea.notclose.h0")
       %in% fit.measures)) {
        rmsea.pvalue.flag <- TRUE
    }


    # 1. RMSEA
    if(rmsea.val.flag) {
        indices["rmsea"] <- lav_fit_rmsea(X2 = X2, df = df, N = N, G = G)
        if(scaled.flag) {
            indices["rmsea.scaled"] <- lav_fit_rmsea(X2 = XX2, df = df2,
                                                     N = N, G = G)
            indices["rmsea.robust"] <- as.numeric(NA)
            if(robust.flag) {
                indices["rmsea.robust"] <-
                    lav_fit_rmsea(X2 = XX3, df = df3, N = N, c.hat = c.hat3,
                                  G = G)
            }
        }
    }

    # 2. RMSEA CI
    if(rmsea.ci.flag) {
        indices["rmsea.ci.level"] <- ci.level
        ci <- lav_fit_rmsea_ci(X2 = X2, df = df, N = N, G = G, level = ci.level)
        indices["rmsea.ci.lower"] <- ci$rmsea.ci.lower
        indices["rmsea.ci.upper"] <- ci$rmsea.ci.upper

        if(scaled.flag) {
            ci.scaled <- lav_fit_rmsea_ci(X2 = XX2, df = df2, N = N, G = G,
                                          level = ci.level)
            indices["rmsea.ci.lower.scaled"] <- ci.scaled$rmsea.ci.lower
            indices["rmsea.ci.upper.scaled"] <- ci.scaled$rmsea.ci.upper
            indices["rmsea.ci.lower.robust"] <- as.numeric(NA)
            indices["rmsea.ci.upper.robust"] <- as.numeric(NA)
            if(robust.flag) {
                # note: input is scaled test statistic!
                ci.robust <- lav_fit_rmsea_ci(X2 = XX3.scaled,
                                 df = df3, N = N, G = G, c.hat = c.hat,
                                 level = ci.level)
                indices["rmsea.ci.lower.robust"] <- ci.robust$rmsea.ci.lower
                indices["rmsea.ci.upper.robust"] <- ci.robust$rmsea.ci.upper
            }
        }
    }

    # 3. RMSEA pvalue
    if(rmsea.pvalue.flag) {
        indices["rmsea.close.h0"]    <- close.h0
        indices["rmsea.notclose.h0"] <- notclose.h0
        indices["rmsea.pvalue"] <-
            lav_fit_rmsea_closefit(X2 = X2, df = df, N = N, G = G,
                                   rmsea.h0 = close.h0)
        indices["rmsea.notclose.pvalue"] <-
            lav_fit_rmsea_notclosefit(X2 = X2, df = df, N = N, G = G,
                                      rmsea.h0 = notclose.h0)
        if(scaled.flag) {
            indices["rmsea.pvalue.scaled"] <-
                lav_fit_rmsea_closefit(X2 = XX2, df = df2, N = N, G = G,
                                       rmsea.h0 = close.h0)
            indices["rmsea.notclose.pvalue.scaled"] <-
                lav_fit_rmsea_notclosefit(X2 = XX2, df = df2, N = N, G = G,
                                          rmsea.h0 = notclose.h0)
            indices["rmsea.pvalue.robust"] <- as.numeric(NA)
            indices["rmsea.notclose.pvalue.robust"] <- as.numeric(NA)
            if(robust.flag) {
                indices["rmsea.pvalue.robust"] <-  # new in 0.6-13
                    lav_fit_rmsea_closefit(X2 = XX3.scaled,
                                           df = df3,
                                           N = N, G = G, c.hat = c.hat,
                                           rmsea.h0 = close.h0)
                indices["rmsea.notclose.pvalue.robust"] <-  # new in 0.6-13
                    lav_fit_rmsea_notclosefit(X2 = XX3.scaled,
                                              df = df3,
                                              N = N, G = G, c.hat = c.hat,
                                              rmsea.h0 = notclose.h0)
            }
        }
    }

    # return only those that were requested
    indices[fit.measures]
}

