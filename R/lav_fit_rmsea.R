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

# categorical data:
# Savalei, V. (2021). Improving fit indices in structural equation modeling with
# categorical data. Multivariate Behavioral Research, 56(3), 390-407. doi:
# 10.1080/00273171.2020.1717922

# missing = "fiml":
# Zhang, X., & Savalei, V. (2022). New computations for RMSEA and CFI following
# FIML and TS estimation with missing data. Psychological Methods.


# always using N (if a user needs N-1, just replace N by N-1)
# vectorized!
lav_fit_rmsea <- function(x2 = NULL, df = NULL, n = NULL,
                          f_val = NULL, g = 1L, c_hat = 1.0) {
  # did we get a sample size?
  if (missing(n) && !missing(f_val)) {
    # population version
    rmsea <- sqrt(f_val / df)
  } else {
    nel <- length(x2)
    if (nel == 0) {
      return(as.numeric(NA))
    }
    rmsea <- ifelse(df > 0,
      # 'standard' way to compute RMSEA
      rmsea <- sqrt(pmax((x2 / n) / df - c_hat / n, rep(0, nel))) * sqrt(g),
      0
    ) # if df == 0
  }

  rmsea
}


# note: for 'robust' version, X2 should be SB-X2
lav_fit_rmsea_ci <- function(x2 = NULL, df = NULL, n = NULL,
                             g = 1L, c_hat = 1, level = 0.90) {
  if (missing(n) || missing(x2) || missing(df) ||
    !is.finite(x2) || !is.finite(df) || !is.finite(n)) {
    return(list(
      rmsea.ci.lower = as.numeric(NA),
      rmsea.ci.upper = as.numeric(NA)
    ))
  }

  if (!is.finite(level) || level < 0 || level > 1.0) {
    lav_msg_warn(gettextf(
      "invalid level value [%s] set to default 0.90.", level
    ))
    level <- 0.90
  }

  upper_perc <- (1 - (1 - level) / 2)
  lower_perc <- (1 - level) / 2

  # internal function
  lower_lambda <- function(lambda) {
    (pchisq(x2, df = df, ncp = lambda) - upper_perc)
  }

  upper_lambda <- function(lambda) {
    (pchisq(x2, df = df, ncp = lambda) - lower_perc)
  }


  # lower bound
  if (df < 1 || lower_lambda(0) < 0.0) {
    rmsea_ci_lower <- 0
  } else {
    lambda_l <- try(uniroot(f = lower_lambda, lower = 0, upper = x2)$root,
      silent = TRUE
    )
    if (inherits(lambda_l, "try-error")) {
      lambda_l <- as.numeric(NA)
    }
    # lower bound
    rmsea_ci_lower <- sqrt((c_hat * lambda_l) / (n * df))

    # multiple groups? -> correction
    if (g > 1L) {
      rmsea_ci_lower <- rmsea_ci_lower * sqrt(g)
    }
  }

  # upper bound
  n_rmsea <- max(n, x2 * 4)
  if (df < 1 || upper_lambda(n_rmsea) > 0 || upper_lambda(0) < 0) {
    rmsea_ci_upper <- 0
  } else {
    lambda_u <- try(
      uniroot(
        f = upper_lambda, lower = 0,
        upper = n_rmsea
      )$root,
      silent = TRUE
    )
    if (inherits(lambda_u, "try-error")) {
      lambda_u <- NA
    }
    # upper bound
    rmsea_ci_upper <- sqrt((c_hat * lambda_u) / (n * df))

    # multiple groups? -> correction
    if (g > 1L) {
      rmsea_ci_upper <- rmsea_ci_upper * sqrt(g)
    }
  }

  list(
    rmsea.ci.lower = rmsea_ci_lower,
    rmsea.ci.upper = rmsea_ci_upper
  )
}

# H_0: RMSEA <= rmsea.h0
lav_fit_rmsea_closefit <- function(x2 = NULL, df = NULL, n = NULL,
                                   g = 1L, c_hat = 1, rmsea_h0 = 0.05) {
  if (missing(n) || missing(x2) || missing(df) ||
    !is.finite(x2) || !is.finite(df) || !is.finite(n)) {
    return(as.numeric(NA))
  }

  rmsea_pvalue <- as.numeric(NA)
  if (df > 0) {
    # see Dudgeon 2004, eq 16 for the 'G' correction
    ncp <- (n * df * 1 / c_hat * rmsea_h0^2) / g
    rmsea_pvalue <- 1 - pchisq(x2, df = df, ncp = ncp)
  }

  rmsea_pvalue
}

# MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996).
# H_0: RMSEA >= rmsea.h0
lav_fit_rmsea_notclosefit <- function(x2 = NULL, df = NULL, n = NULL,
                                      g = 1L, c_hat = 1, rmsea_h0 = 0.05) {
  if (missing(n) || missing(x2) || missing(df) ||
    !is.finite(x2) || !is.finite(df) || !is.finite(n)) {
    return(as.numeric(NA))
  }

  rmsea_pvalue <- as.numeric(NA)
  if (df > 0) {
    # see Dudgeon 2004, eq 16 for the 'G' correction
    ncp <- (n * df * 1 / c_hat * rmsea_h0^2) / g
    rmsea_pvalue <- pchisq(x2, df = df, ncp = ncp)
  }

  rmsea_pvalue
}


lav_fit_rmsea_lavobject <- function(lavobject = NULL, fit_measures = "rmsea",
                                    standard_test = "standard",
                                    scaled_test = "none",
                                    ci_level = 0.90,
                                    close_h0 = 0.05, notclose_h0 = 0.08,
                                    robust = TRUE,
                                    cat_nonpd = "na") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # check for categorical
  categorical_flag <- lavobject@Model@categorical

  # tests
  test <- lavobject@test
  test_names <- sapply(lavobject@test, "[[", "test")
  if (test_names[1] == "none" || standard_test == "none") {
    return(list())
  }
  test_idx <- which(test_names == standard_test)[1]
  if (length(test_idx) == 0L) {
    return(list())
  }

  scaled_flag <- FALSE
  if (!scaled_test %in% c("none", "standard", "default")) {
    scaled_idx <- which(test_names == scaled_test)
    if (length(scaled_idx) > 0L) {
      scaled_idx <- scaled_idx[1] # only the first one
      scaled_flag <- TRUE
    }
  }

  # robust?
  robust_flag <- FALSE
  if (robust && scaled_flag &&
    scaled_test %in% c(
      "satorra.bentler", "yuan.bentler.mplus",
      "yuan.bentler", "yuan.chan", "scaled.shifted",
      "scaled.shifted.corrected"
    )) {
    robust_flag <- TRUE
  }

  # FIML? (the FIML-corrected robust indices are only available for
  # single-level models; for multilevel data, we fall back to the regular
  # scaled/robust computation above, which handles missing data too)
  fiml_flag <- FALSE
  if (robust && lavobject@Options$missing %in% c("ml", "ml.x") &&
    lavobject@Data@nlevels == 1L) {
    fiml_flag <- robust_flag <- TRUE
    # check if we can compute corrected values
    if (scaled_flag) {
      version <- "V3"
    } else {
      version <- "V6"
    }
    fiml <- try(
      lav_fit_fiml_corrected(lavobject,
        baseline_model = NULL,
        version = version
      ),
      silent = TRUE
    )
    if (inherits(fiml, "try-error")) {
      lav_msg_warn(gettext("computation of robust RMSEA failed."))
      fiml <- list(
        XX3 = as.numeric(NA), df3 = as.numeric(NA),
        c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA)
      )
    } else if (anyNA(c(fiml$XX3, fiml$df3, fiml$c.hat3, fiml$XX3.scaled))) {
      lav_msg_warn(gettext(
        "computation of robust RMSEA resulted in NA values."
      ))
    }
  }

  # supported fit measures in this function
  fit_rmsea <- c(
    "rmsea",
    "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.ci.level",
    "rmsea.pvalue", "rmsea.close.h0",
    "rmsea.notclose.pvalue", "rmsea.notclose.h0"
  )
  if (scaled_flag) {
    fit_rmsea <- c(
      fit_rmsea, "rmsea.scaled", "rmsea.ci.lower.scaled",
      "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled",
      "rmsea.notclose.pvalue.scaled"
    )
  }
  if (robust_flag) {
    fit_rmsea <- c(
      fit_rmsea, "rmsea.robust", "rmsea.ci.lower.robust",
      "rmsea.ci.upper.robust", "rmsea.pvalue.robust",
      "rmsea.notclose.pvalue.robust"
    )
  }

  # which one do we need?
  if (missing(fit_measures)) {
    # default set
    fit_measures <- fit_rmsea
  } else {
    # remove any not-RMSEA related index from fit.measures
    rm_idx <- which(!fit_measures %in% fit_rmsea)
    if (length(rm_idx) > 0L) {
      fit_measures <- fit_measures[-rm_idx]
    }
    if (length(fit_measures) == 0L) {
      return(list())
    }
  }


  # basic test statistics
  x2 <- test[[test_idx]]$stat
  df <- test[[test_idx]]$df
  g <- lavobject@Data@ngroups # number of groups
  n <- lav_inspect_ntotal(object = lavobject) # N vs N-1

  # scaled X2/df values
  if (scaled_flag) {
    if (scaled_test %in% c("scaled.shifted", "scaled.shifted.corrected")) {
      xx2 <- test[[scaled_idx]]$stat
      df2 <- df
    } else {
      xx2 <- x2
      df2 <- sum(test[[scaled_idx]]$trace.UGamma)
      if (!is.finite(df2) || df2 == 0) {
        df2 <- as.numeric(NA)
      }
    }
  }

  # robust ingredients
  n3 <- n
  if (robust_flag) {
    if (categorical_flag) {
      out <- try(lav_fit_catml_dwls(lavobject, nonpd = cat_nonpd),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        xx3 <- df3 <- c_hat3 <- xx3_scaled <- as.numeric(NA)
      } else {
        xx3 <- out$XX3
        df3 <- out$df3
        c_hat3 <- c_hat <- out$c.hat3
        xx3_scaled <- out$XX3.scaled
      }
      # new in 0.7-1: the robust quantities (xx3, c.hat3, ...) are all
      # ML-based (chi-square = N * F_ML), so we use N -- not N-1 -- in
      # the denominator; this also matches what a direct
      # estimator = "catML" run would use
      n3 <- lavobject@SampleStats@ntotal
    } else if (fiml_flag) {
      xx3 <- fiml$XX3
      df3 <- fiml$df3
      c_hat3 <- c_hat <- fiml$c.hat3
      xx3_scaled <- fiml$XX3.scaled
    } else {
      xx3 <- x2
      df3 <- df
      c_hat <- test[[scaled_idx]]$scaling.factor
      if (scaled_test %in% c(
        "scaled.shifted", "scaled.shifted.corrected"
      )) {
        # compute c.hat from a and b
        a <- test[[scaled_idx]]$scaling.factor
        b <- test[[scaled_idx]]$shift.parameter
        c_hat3 <- a * (df - b) / df
      } else {
        c_hat3 <- c_hat
      }
      xx3_scaled <- test[[scaled_idx]]$stat
    }
  }

  # output container
  indices <- list()

  # what do we need?
  rmsea_val_flag <- rmsea_ci_flag <- rmsea_pvalue_flag <- FALSE
  if (any(c("rmsea", "rmsea.scaled", "rmsea.robust") %in% fit_measures)) {
    rmsea_val_flag <- TRUE
  }
  if (any(c(
    "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.ci.level",
    "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled",
    "rmsea.ci.lower.robust", "rmsea.ci.upper.robust"
  )
  %in% fit_measures)) {
    rmsea_ci_flag <- TRUE
  }
  if (any(c(
    "rmsea.pvalue", "rmsea.pvalue.scaled", "rmsea.pvalue.robust",
    "rmsea.notclose.pvalue", "rmsea.notclose.pvalue.scaled",
    "rmsea.notclose.pvalue.robust",
    "rmsea.close.h0", "rmsea.notclose.h0"
  )
  %in% fit_measures)) {
    rmsea_pvalue_flag <- TRUE
  }


  # 1. RMSEA
  if (rmsea_val_flag) {
    indices["rmsea"] <- lav_fit_rmsea(x2 = x2, df = df, n = n, g = g)
    if (scaled_flag) {
      indices["rmsea.scaled"] <- lav_fit_rmsea(
        x2 = xx2, df = df2,
        n = n, g = g
      )
    }
    if (robust_flag) {
      indices["rmsea.robust"] <-
        lav_fit_rmsea(x2 = xx3, df = df3, n = n3, c_hat = c_hat3, g = g)
    }
  }

  # 2. RMSEA CI
  if (rmsea_ci_flag) {
    indices["rmsea.ci.level"] <- ci_level
    ci <- lav_fit_rmsea_ci(x2 = x2, df = df, n = n, g = g, level = ci_level)
    indices["rmsea.ci.lower"] <- ci$rmsea.ci.lower
    indices["rmsea.ci.upper"] <- ci$rmsea.ci.upper

    if (scaled_flag) {
      ci_scaled <- lav_fit_rmsea_ci(
        x2 = xx2, df = df2, n = n, g = g,
        level = ci_level
      )
      indices["rmsea.ci.lower.scaled"] <- ci_scaled$rmsea.ci.lower
      indices["rmsea.ci.upper.scaled"] <- ci_scaled$rmsea.ci.upper
    }
    if (robust_flag) {
      # note: input is scaled test statistic!
      ci_robust <- lav_fit_rmsea_ci(
        x2 = xx3_scaled,
        df = df3, n = n3, g = g, c_hat = c_hat,
        level = ci_level
      )
      indices["rmsea.ci.lower.robust"] <- ci_robust$rmsea.ci.lower
      indices["rmsea.ci.upper.robust"] <- ci_robust$rmsea.ci.upper
    }
  }

  # 3. RMSEA pvalue
  if (rmsea_pvalue_flag) {
    indices["rmsea.close.h0"] <- close_h0
    indices["rmsea.notclose.h0"] <- notclose_h0
    indices["rmsea.pvalue"] <-
      lav_fit_rmsea_closefit(
        x2 = x2, df = df, n = n, g = g,
        rmsea_h0 = close_h0
      )
    indices["rmsea.notclose.pvalue"] <-
      lav_fit_rmsea_notclosefit(
        x2 = x2, df = df, n = n, g = g,
        rmsea_h0 = notclose_h0
      )
    if (scaled_flag) {
      indices["rmsea.pvalue.scaled"] <-
        lav_fit_rmsea_closefit(
          x2 = xx2, df = df2, n = n, g = g,
          rmsea_h0 = close_h0
        )
      indices["rmsea.notclose.pvalue.scaled"] <-
        lav_fit_rmsea_notclosefit(
          x2 = xx2, df = df2, n = n, g = g,
          rmsea_h0 = notclose_h0
        )
    }
    if (robust_flag) {
      indices["rmsea.pvalue.robust"] <- # new in 0.6-13
        lav_fit_rmsea_closefit(
          x2 = xx3_scaled,
          df = df3,
          n = n3, g = g, c_hat = c_hat,
          rmsea_h0 = close_h0
        )
      indices["rmsea.notclose.pvalue.robust"] <- # new in 0.6-13
        lav_fit_rmsea_notclosefit(
          x2 = xx3_scaled,
          df = df3,
          n = n3, g = g, c_hat = c_hat,
          rmsea_h0 = notclose_h0
        )
    }
  }

  # return only those that were requested
  indices[fit_measures]
}
