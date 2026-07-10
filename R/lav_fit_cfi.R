# functions related to CFI and other 'incremental' fit indices

# lower-level functions:
# - lav_fit_cfi
# - lav_fit_rni (same as CFI, but without the max(0,))
# - lav_fit_tli/lav_fit_nnfi
# - lav_fit_rfi
# - lav_fit_nfi
# - lav_fit_pnfi
# - lav_fit_ifi

# higher-level functions:
# - lav_fit_cfi_lavobject

# Y.R. 20 July 2022

# CFI - comparative fit index (Bentler, 1990)
# robust version: Brosseau-Liard & Savalei MBR 2014, equation 15
# Brosseau-Liard, P. E., & Savalei, V. (2014). Adjusting incremental fit
# indices for nonnormality. Multivariate behavioral research, 49(5), 460-470.

# robust version MLMV (scaled.shifted)
# Savalei, V. (2018). On the computation of the RMSEA and CFI from the
# mean-and-variance corrected test statistic with nonnormal data in SEM.
# Multivariate behavioral research, 53(3), 419-429. eq 9

# note: robust MLM == robust MLMV

# categorical data:
# Savalei, V. (2021). Improving fit indices in structural equation modeling with
# categorical data. Multivariate Behavioral Research, 56(3), 390-407. doi:
# 10.1080/00273171.2020.1717922

# when missing = "fiml":
# Zhang, X., & Savalei, V. (2023). New computations for RMSEA and CFI following
# FIML and TS estimation with missing data. Psychological Methods, 28(2),
# 263-283. https://doi.org/10.1037/met0000445


lav_fit_cfi <- function(x2 = NULL, df = NULL, x2_null = NULL, df_null = NULL,
                        c_hat = 1, c_hat_null = 1) {
  if (anyNA(c(x2, df, x2_null, df_null, c_hat, c_hat_null))) {
    return(as.numeric(NA))
  }

  # robust?
  if (df > 0 && !missing(c_hat) && !missing(c_hat_null) &&
    c_hat != 1 && c_hat_null != 1) {
    t1 <- max(c(x2 - (c_hat * df), 0))
    t2 <- max(c(x2 - (c_hat * df), x2_null - (c_hat_null * df_null), 0))
  } else {
    t1 <- max(c(x2 - df, 0))
    t2 <- max(c(x2 - df, x2_null - df_null, 0))
  }

  if (isTRUE(all.equal(t1, 0)) && isTRUE(all.equal(t2, 0))) {
    cfi <- 1
  } else {
    cfi <- 1 - t1 / t2
  }

  cfi
}

# RNI - relative noncentrality index (McDonald & Marsh, 1990)
# same as CFI, but without the max(0,)
lav_fit_rni <- function(x2 = NULL, df = NULL, x2_null = NULL, df_null = NULL,
                        c_hat = 1, c_hat_null = 1) {
  if (anyNA(c(x2, df, x2_null, df_null, c_hat, c_hat_null))) {
    return(as.numeric(NA))
  }

  # robust?
  if (df > 0 && !missing(c_hat) && !missing(c_hat_null) &&
    c_hat != 1 && c_hat_null != 1) {
    t1 <- x2 - (c_hat * df)
    t2 <- x2_null - (c_hat_null * df_null)
  } else {
    t1 <- x2 - df
    t2 <- x2_null - df_null
  }

  if (isTRUE(all.equal(t2, 0))) {
    rni <- as.numeric(NA)
  } else if (!is.finite(t1) || !is.finite(t2)) {
    rni <- as.numeric(NA)
  } else {
    rni <- 1 - t1 / t2
  }

  rni
}

# TLI - Tucker-Lewis index (Tucker & Lewis, 1973)
# same as
# NNFI - nonnormed fit index (NNFI, Bentler & Bonett, 1980)
# note: formula in lavaan <= 0.5-20:
# t1 <- X2.null/df.null - X2/df
# t2 <- X2.null/df.null - 1
# if(t1 < 0 && t2 < 0) {
#    TLI <- 1
# } else {
#    TLI <- t1/t2
# }
# note: TLI original formula was in terms of fx/df, not X2/df
# then, t1 <- fx_0/df.null - fx/df
#       t2 <- fx_0/df.null - 1/N (or N-1 for wishart)

# note: in lavaan 0.5-21, we use the alternative formula:
# TLI <- 1 - ((X2 - df)/(X2.null - df.null) * df.null/df)
# - this one has the advantage that a 'robust' version
#   can be derived; this seems non-trivial for the original one
# - unlike cfi, we do not use 'max(0, )' for t1 and t2
#   therefore, t1 can go negative, and TLI can be > 1
lav_fit_tli <- function(x2 = NULL, df = NULL, x2_null = NULL, df_null = NULL,
                        c_hat = 1, c_hat_null = 1) {
  if (anyNA(c(x2, df, x2_null, df_null, c_hat, c_hat_null))) {
    return(as.numeric(NA))
  }

  # robust?
  if (df > 0 && !missing(c_hat) && !missing(c_hat_null) &&
    c_hat != 1 && c_hat_null != 1) {
    t1 <- (x2 - c_hat * df) * df_null
    t2 <- (x2_null - c_hat_null * df_null) * df
  } else {
    t1 <- (x2 - df) * df_null
    t2 <- (x2_null - df_null) * df
  }

  if (df > 0 && abs(t2) > 0) {
    tli <- 1 - t1 / t2
  } else if (!is.finite(t1) || !is.finite(t2)) {
    tli <- as.numeric(NA)
  } else {
    tli <- 1
  }

  tli
}

# alias for nnfi
lav_fit_nnfi <- lav_alias("lav_fit_tli")

# RFI - relative fit index (Bollen, 1986; Joreskog & Sorbom 1993)
lav_fit_rfi <- function(x2 = NULL, df = NULL, x2_null = NULL, df_null = NULL) {
  if (anyNA(c(x2, df, x2_null, df_null))) {
    return(as.numeric(NA))
  }

  if (df > df_null) {
    rli <- as.numeric(NA)
  } else if (df > 0 && df_null > 0) {
    t1 <- x2_null / df_null - x2 / df
    t2 <- x2_null / df_null
    if (!is.finite(t1) || !is.finite(t2)) {
      rli <- as.numeric(NA)
    } else if (t1 < 0 || t2 < 0) {
      rli <- 1
    } else {
      rli <- t1 / t2
    }
  } else {
    rli <- 1
  }

  rli
}

# NFI - normed fit index (Bentler & Bonett, 1980)
lav_fit_nfi <- function(x2 = NULL, df = NULL, x2_null = NULL, df_null = NULL) {
  if (anyNA(c(x2, df, x2_null, df_null))) {
    return(as.numeric(NA))
  }

  if (df > df_null || isTRUE(all.equal(x2_null, 0))) {
    nfi <- as.numeric(NA)
  } else if (df > 0) {
    t1 <- x2_null - x2
    t2 <- x2_null
    nfi <- t1 / t2
  } else {
    nfi <- 1
  }

  nfi
}

# PNFI - Parsimony normed fit index (James, Mulaik & Brett, 1982)
lav_fit_pnfi <- function(x2 = NULL, df = NULL, x2_null = NULL, df_null = NULL) {
  if (anyNA(c(x2, df, x2_null, df_null))) {
    return(as.numeric(NA))
  }

  if (df_null > 0 && x2_null > 0) {
    t1 <- x2_null - x2
    t2 <- x2_null
    pnfi <- (df / df_null) * (t1 / t2)
  } else {
    pnfi <- as.numeric(NA)
  }

  pnfi
}

# IFI - incremental fit index (Bollen, 1989; Joreskog & Sorbom, 1993)
lav_fit_ifi <- function(x2 = NULL, df = NULL, x2_null = NULL, df_null = NULL) {
  if (anyNA(c(x2, df, x2_null, df_null))) {
    return(as.numeric(NA))
  }

  t1 <- x2_null - x2
  t2 <- x2_null - df
  if (!is.finite(t1) || !is.finite(t2)) {
    ifi <- as.numeric(NA)
  } else if (t2 < 0) {
    ifi <- 1
  } else if (isTRUE(all.equal(t2, 0))) {
    ifi <- as.numeric(NA)
  } else {
    ifi <- t1 / t2
  }

  ifi
}

# higher-level function
lav_fit_cfi_lavobject <- function(lavobject = NULL, fit_measures = "cfi",
                                  baseline_model = NULL, h1_model = NULL,
                                  standard_test = "standard",
                                  scaled_test = "none",
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
      "yuan.bentler", "yuan.chan", "scaled.shifted"
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
      lav_fit_fiml_corrected(lavobject, baseline_model,
        version = version
      ),
      silent = TRUE
    )
    if (inherits(fiml, "try-error")) {
      lav_msg_warn(gettext("computation of robust CFI failed."))
      fiml <- list(
        XX3 = as.numeric(NA), df3 = as.numeric(NA),
        c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
        XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
        c.hat3.null = as.numeric(NA)
      )
    } else if (anyNA(c(
      fiml$XX3, fiml$df3, fiml$c.hat3, fiml$XX3.scaled,
      fiml$XX3.null, fiml$df3.null, fiml$c.hat3.null
    ))) {
      lav_msg_warn(gettext("computation of robust CFI resulted in NA values."))
    }
  }

  # supported fit measures in this function
  # baseline model
  fit_baseline <- c("baseline.chisq", "baseline.df", "baseline.pvalue")
  if (scaled_flag) {
    fit_baseline <- c(
      fit_baseline, "baseline.chisq.scaled",
      "baseline.df.scaled", "baseline.pvalue.scaled",
      "baseline.chisq.scaling.factor"
    )
  }

  fit_cfi_tli <- c("cfi", "tli")
  if (scaled_flag) {
    fit_cfi_tli <- c(fit_cfi_tli, "cfi.scaled", "tli.scaled")
  }
  if (robust_flag) {
    fit_cfi_tli <- c(fit_cfi_tli, "cfi.robust", "tli.robust")
  }

  # other incremental fit indices
  fit_cfi_other <- c("nnfi", "rfi", "nfi", "pnfi", "ifi", "rni")
  if (scaled_flag) {
    fit_cfi_other <- c(
      fit_cfi_other, "nnfi.scaled", "rfi.scaled",
      "nfi.scaled", "pnfi.scaled", "ifi.scaled", "rni.scaled"
    )
  }
  if (robust_flag) {
    fit_cfi_other <- c(fit_cfi_other, "nnfi.robust", "rni.robust")
  }

  # which one do we need?
  if (missing(fit_measures)) {
    # default set
    fit_measures <- c(fit_baseline, fit_cfi_tli)
  } else {
    # remove any not-CFI related index from fit.measures
    rm_idx <- which(!fit_measures %in%
      c(fit_baseline, fit_cfi_tli, fit_cfi_other))
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
  # g <- lavobject@Data@ngroups # number of groups
  # n <- lav_inspect_ntotal(object = lavobject) # N vs N-1

  # scaled X2
  if (scaled_flag) {
    x2_scaled <- test[[scaled_idx]]$stat
    df_scaled <- test[[scaled_idx]]$df
  }
  if (robust_flag) {
    xx3 <- x2
    if (categorical_flag) {
      out <- try(lav_fit_catml_dwls(lavobject, nonpd = cat_nonpd),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        xx3 <- df3 <- c_hat <- c_hat3 <- xx3_scaled <- as.numeric(NA)
      } else {
        xx3 <- out$XX3
        df3 <- out$df3
        c_hat3 <- c_hat <- out$c.hat3
        xx3_scaled <- out$XX3.scaled
      }
    } else if (fiml_flag) {
      xx3 <- fiml$XX3
      df3 <- fiml$df3
      c_hat3 <- c_hat <- fiml$c.hat3
      xx3_scaled <- fiml$XX3.scaled
    } else if (scaled_test == "scaled.shifted") {
      # compute c.hat from a and b
      a <- test[[scaled_idx]]$scaling.factor
      b <- test[[scaled_idx]]$shift.parameter
      c_hat <- a * (df - b) / df
    } else {
      c_hat <- test[[scaled_idx]]$scaling.factor
    }
  }

  # output container
  indices <- list()

  # only do what is needed (per groups)
  cfi_baseline_flag <- cfi_tli_flag <- cfi_other_flag <- FALSE
  if (any(fit_baseline %in% fit_measures)) {
    cfi_baseline_flag <- TRUE
  }
  if (any(fit_cfi_tli %in% fit_measures)) {
    cfi_tli_flag <- TRUE
  }
  if (any(fit_cfi_other %in% fit_measures)) {
    cfi_other_flag <- TRUE
  }

  # 1. BASELINE model
  baseline_test <- NULL

  # is the (effective) baseline model the default independence model?
  # (needed below: the robust corrections for categorical data -- Savalei,
  #  2021, eq. 6b -- are only valid for the independence baseline)
  baseline_default_flag <- is.null(baseline_model) &&
    is.null(lavobject@external$baseline.model) &&
    !identical(lavobject@Options$baseline.type, "nested")

  # we use the following priority:
  # 1. user-provided baseline model
  # 2. baseline model in @external slot
  # 3. baseline model in @baseline slot
  # 4. nothing -> compute independence model

  # TDJ: Also check for user-supplied h1.model, using similar priority:
  #        1. user-provided h1 model
  #        2. h1 model in @external slot
  #        3. default h1 model (already in @h1 slot, no update necessary)

  # 1. user-provided h1 model
  if (!is.null(h1_model)) {
    stopifnot(inherits(h1_model, "lavaan"))

    # 2. h1 model in @external slot
  } else if (!is.null(lavobject@external$h1.model)) {
    stopifnot(inherits(lavobject@external$h1.model, "lavaan"))
    h1_model <- lavobject@external$h1.model
  } # else is.null

  # 1. user-provided baseline model
  if (!is.null(baseline_model)) {
    baseline_test <-
      lav_fit_check_baseline(
        fit_indep = baseline_model,
        object = lavobject,
        fit_h1 = h1_model # okay if NULL
      )
    # 2. baseline model in @external slot
  } else if (!is.null(lavobject@external$baseline.model)) {
    fit_indep <- lavobject@external$baseline.model
    baseline_test <-
      lav_fit_check_baseline(
        fit_indep = fit_indep,
        object = lavobject,
        fit_h1 = h1_model # okay if NULL
      )
    # 3. internal @baseline slot
  } else if (length(lavobject@baseline) > 0L &&
    !is.null(lavobject@baseline$test) &&
    ## if there is a custom h1.model, need  _check_baseline() to update @test
    is.null(h1_model)) {
    baseline_test <- lavobject@baseline$test
    # 4. (re)compute independence model
  } else {
    fit_indep <- try(lav_object_independence(lavobject), silent = TRUE)
    baseline_test <-
      lav_fit_check_baseline(
        fit_indep = fit_indep,
        object = lavobject,
        fit_h1 = h1_model # okay if NULL
      )
  }

  # baseline_test_idx
  baseline_test_idx <- which(names(baseline_test) == standard_test)[1]
  baseline_scaled_ok <- FALSE
  if (scaled_flag) {
    baseline_scaled_idx <- which(names(baseline_test) == scaled_test)[1]
    # the baseline model may lack the scaled test entry (eg a user-provided
    # baseline model that was fitted without the scaled test, or a baseline
    # that could not be refitted with the scaled test): degrade to NA for
    # the scaled/robust baseline quantities instead of failing
    baseline_scaled_ok <- !is.na(baseline_scaled_idx)
  }

  if (!is.null(baseline_test)) {
    x2_null <- baseline_test[[baseline_test_idx]]$stat
    df_null <- baseline_test[[baseline_test_idx]]$df
    if (scaled_flag) {
      if (baseline_scaled_ok) {
        x2_null_scaled <- baseline_test[[baseline_scaled_idx]]$stat
        df_null_scaled <- baseline_test[[baseline_scaled_idx]]$df
      } else {
        x2_null_scaled <- df_null_scaled <- as.numeric(NA)
      }
    }
    if (robust_flag) {
      xx3_null <- x2_null
      if (categorical_flag) {
        if (inherits(out, "try-error")) {
          xx3_null <- c_hat_null <- as.numeric(NA)
        } else if (!baseline_default_flag) {
          # the robust corrections for categorical data assume the
          # (default) independence baseline model; for any other baseline
          # model, mixing its test statistic/df with the independence-
          # based correction (tr Gamma) would be wrong -> return NA
          lav_msg_warn(gettext(
            "robust fit indices are currently not available for
             categorical data in combination with a non-default baseline
             model; returning NA."))
          xx3_null <- c_hat_null <- as.numeric(NA)
        } else {
          xx3_null <- out$XX3.null
          c_hat_null <- out$c.hat3.null
        }
      } else if (fiml_flag) {
        xx3_null <- fiml$XX3.null
        c_hat_null <- fiml$c.hat3.null
      } else if (!baseline_scaled_ok) {
        c_hat_null <- as.numeric(NA)
      } else if (scaled_test == "scaled.shifted") {
        # compute c.hat from a and b
        a_null <-
          baseline_test[[baseline_scaled_idx]]$scaling.factor
        b_null <-
          baseline_test[[baseline_scaled_idx]]$shift.parameter
        c_hat_null <- a_null * (df_null - b_null) / df_null
      } else {
        c_hat_null <-
          baseline_test[[baseline_scaled_idx]]$scaling.factor
      }
    }
  } else {
    x2_null <- df_null <- as.numeric(NA)
    x2_null_scaled <- df_null_scaled <- as.numeric(NA)
    c_hat_null <- as.numeric(NA)
  }

  # check for NAs of nonfinite numbers
  if (!is.finite(x2) || !is.finite(df) ||
    !is.finite(x2_null) || !is.finite(df_null)) {
    indices[fit_measures] <- as.numeric(NA)
    return(indices)
  }

  # fill in baseline indices
  if (cfi_baseline_flag) {
    indices["baseline.chisq"] <- x2_null
    indices["baseline.df"] <- df_null
    indices["baseline.pvalue"] <- baseline_test[[baseline_test_idx]]$pvalue
    if (scaled_flag) {
      indices["baseline.chisq.scaled"] <- x2_null_scaled
      indices["baseline.df.scaled"] <- df_null_scaled
      if (baseline_scaled_ok) {
        indices["baseline.pvalue.scaled"] <-
          baseline_test[[baseline_scaled_idx]]$pvalue
        indices["baseline.chisq.scaling.factor"] <-
          baseline_test[[baseline_scaled_idx]]$scaling.factor
      } else {
        indices["baseline.pvalue.scaled"] <- as.numeric(NA)
        indices["baseline.chisq.scaling.factor"] <- as.numeric(NA)
      }
    }
  }

  # 2. CFI and TLI
  if (cfi_tli_flag) {
    indices["cfi"] <- lav_fit_cfi(
      x2 = x2, df = df,
      x2_null = x2_null, df_null = df_null
    )
    indices["tli"] <- lav_fit_tli(
      x2 = x2, df = df,
      x2_null = x2_null, df_null = df_null
    )
    if (scaled_flag) {
      indices["cfi.scaled"] <-
        lav_fit_cfi(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
      indices["tli.scaled"] <-
        lav_fit_tli(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
    }
    if (robust_flag) {
      indices["cfi.robust"] <-
        lav_fit_cfi(
          x2 = xx3, df = df,
          x2_null = xx3_null, df_null = df_null,
          c_hat = c_hat, c_hat_null = c_hat_null
        )
      indices["tli.robust"] <-
        lav_fit_tli(
          x2 = xx3, df = df,
          x2_null = xx3_null, df_null = df_null,
          c_hat = c_hat, c_hat_null = c_hat_null
        )
    }
  }

  # 3. other
  # c("nnfi", "rfi", "nfi", "pnfi", "ifi", "rni")
  if (cfi_other_flag) {
    indices["nnfi"] <-
      lav_fit_nnfi(x2 = x2, df = df, x2_null = x2_null, df_null = df_null)
    indices["rfi"] <-
      lav_fit_rfi(x2 = x2, df = df, x2_null = x2_null, df_null = df_null)
    indices["nfi"] <-
      lav_fit_nfi(x2 = x2, df = df, x2_null = x2_null, df_null = df_null)
    indices["pnfi"] <-
      lav_fit_pnfi(x2 = x2, df = df, x2_null = x2_null, df_null = df_null)
    indices["ifi"] <-
      lav_fit_ifi(x2 = x2, df = df, x2_null = x2_null, df_null = df_null)
    indices["rni"] <-
      lav_fit_rni(x2 = x2, df = df, x2_null = x2_null, df_null = df_null)

    if (scaled_flag) {
      indices["nnfi.scaled"] <-
        lav_fit_nnfi(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
      indices["rfi.scaled"] <-
        lav_fit_rfi(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
      indices["nfi.scaled"] <-
        lav_fit_nfi(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
      indices["pnfi.scaled"] <-
        lav_fit_pnfi(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
      indices["ifi.scaled"] <-
        lav_fit_ifi(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
      indices["rni.scaled"] <-
        lav_fit_rni(
          x2 = x2_scaled, df = df_scaled,
          x2_null = x2_null_scaled, df_null = df_null_scaled
        )
    }
    if (robust_flag) {
      indices["nnfi.robust"] <-
        lav_fit_nnfi(
          x2 = xx3, df = df,
          x2_null = xx3_null, df_null = df_null,
          c_hat = c_hat, c_hat_null = c_hat_null
        )
      indices["rni.robust"] <-
        lav_fit_rni(
          x2 = xx3, df = df,
          x2_null = xx3_null, df_null = df_null,
          c_hat = c_hat, c_hat_null = c_hat_null
        )
    }
  }

  # return only those that were requested
  indices[fit_measures]
}


# new in 0.6-5
# internal function to check the (external) baseline model, and
# return baseline 'test' list if everything checks out (and NULL otherwise)
lav_fit_check_baseline <- function(fit_indep = NULL, object = NULL,
                                            fit_h1 = NULL) {
  test <- NULL

  # check if everything is in order
  if (inherits(fit_indep, "try-error")) {
    lav_msg_warn(gettext("baseline model estimation failed"))
    return(NULL)
  } else if (!inherits(fit_indep, "lavaan")) {
    lav_msg_warn(gettext(
      "(user-provided) baseline model is not a fitted lavaan object"
    ))
    return(NULL)
  } else if (!fit_indep@optim$converged) {
    lav_msg_warn(gettext("baseline model did not converge"))
    return(NULL)
  } else {
    # evaluate if estimator/test matches original object
    # note: we do not need to check for 'se', as it may be 'none'
    same_test <- all(object@Options$test == fit_indep@Options$test)
    if (!same_test) {
      lav_msg_warn(gettextf(
        "Baseline model was using test(s) = %1$s, but original model was using
        test(s) = %2$s. Refitting baseline model!",
        lav_msg_view(fit_indep@Options$test, "none"),
        lav_msg_view(object@Options$test, "none")
      ))
    }
    same_estimator <- (object@Options$estimator ==
      fit_indep@Options$estimator)
    if (!same_estimator) {
      lav_msg_warn(gettextf(
        "Baseline model was using estimator = %1$s, but original model was
        using estimator = %2$s. Refitting baseline model!",
        dQuote(fit_indep@Options$estimator),
        dQuote(object@Options$estimator)
      ))
    }
    if (!same_test || !same_estimator) {
      lavoptions <- object@Options
      lavoptions$estimator <- object@Options$estimator
      lavoptions$se <- "none"
      lavoptions$baseline <- FALSE
      lavoptions$check.start <- FALSE
      lavoptions$check.post <- FALSE
      lavoptions$check.vcov <- FALSE
      lavoptions$test <- object@Options$test
      fit_indep <- try(
        lavaan(fit_indep,
          slot_options = lavoptions,
          slot_data = object@Data,
          slot_sample_stats = object@SampleStats,
          sloth1 = object@h1,
          slot_cache = object@Cache,
          verbose = FALSE
        ),
        silent = TRUE
      )
      # try again
      test <- lav_fit_check_baseline(
        fit_indep = fit_indep,
        object = object
      )
    } else {
      # extract what we need
      test <- fit_indep@test
    }
  } # converged lavaan object



  # TDJ: Check for user-supplied h1.model (here, the fit_h1 == argument)
  #      Similar to BASELINE model, use the following priority:
  #        1. user-provided h1 model
  #        2. h1 model in @external slot
  #        3. default h1 model (already in @h1 slot, no update necessary)
  # FIXME? user-supplied h1 model in object might be in fit.indep, too

  user_h1_exists <- FALSE
  # 1. user-provided h1 model
  if (!is.null(fit_h1)) {
    stopifnot(inherits(fit_h1, "lavaan"))
    user_h1_exists <- TRUE

    # 2. h1 model in @external slot
  } else if (!is.null(object@external$h1.model)) {
    stopifnot(inherits(object@external$h1.model, "lavaan"))
    fit_h1 <- object@external$h1.model
    user_h1_exists <- TRUE
  }

  if (user_h1_exists) {
    ## update @test slot
    test <- lav_update_test_custom_h1(
      lav_obj_h0 = fit_indep,
      lav_obj_h1 = fit_h1
    )@test
  }

  test
}
