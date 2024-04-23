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


lav_fit_cfi <- function(X2 = NULL, df = NULL, X2.null = NULL, df.null = NULL,
                        c.hat = 1, c.hat.null = 1) {
  if (anyNA(c(X2, df, X2.null, df.null, c.hat, c.hat.null))) {
    return(as.numeric(NA))
  }

  # robust?
  if (df > 0 && !missing(c.hat) && !missing(c.hat.null) &&
    c.hat != 1 && c.hat.null != 1) {
    t1 <- max(c(X2 - (c.hat * df), 0))
    t2 <- max(c(X2 - (c.hat * df), X2.null - (c.hat.null * df.null), 0))
  } else {
    t1 <- max(c(X2 - df, 0))
    t2 <- max(c(X2 - df, X2.null - df.null, 0))
  }

  if (isTRUE(all.equal(t1, 0)) && isTRUE(all.equal(t2, 0))) {
    CFI <- 1
  } else {
    CFI <- 1 - t1 / t2
  }

  CFI
}

# RNI - relative noncentrality index (McDonald & Marsh, 1990)
# same as CFI, but without the max(0,)
lav_fit_rni <- function(X2 = NULL, df = NULL, X2.null = NULL, df.null = NULL,
                        c.hat = 1, c.hat.null = 1) {
  if (anyNA(c(X2, df, X2.null, df.null, c.hat, c.hat.null))) {
    return(as.numeric(NA))
  }

  # robust?
  if (df > 0 && !missing(c.hat) && !missing(c.hat.null) &&
    c.hat != 1 && c.hat.null != 1) {
    t1 <- X2 - (c.hat * df)
    t2 <- X2.null - (c.hat.null * df.null)
  } else {
    t1 <- X2 - df
    t2 <- X2.null - df.null
  }

  if (isTRUE(all.equal(t2, 0))) {
    RNI <- as.numeric(NA)
  } else if (!is.finite(t1) || !is.finite(t2)) {
    RNI <- as.numeric(NA)
  } else {
    RNI <- 1 - t1 / t2
  }

  RNI
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
lav_fit_tli <- function(X2 = NULL, df = NULL, X2.null = NULL, df.null = NULL,
                        c.hat = 1, c.hat.null = 1) {
  if (anyNA(c(X2, df, X2.null, df.null, c.hat, c.hat.null))) {
    return(as.numeric(NA))
  }

  # robust?
  if (df > 0 && !missing(c.hat) && !missing(c.hat.null) &&
    c.hat != 1 && c.hat.null != 1) {
    t1 <- (X2 - c.hat * df) * df.null
    t2 <- (X2.null - c.hat.null * df.null) * df
  } else {
    t1 <- (X2 - df) * df.null
    t2 <- (X2.null - df.null) * df
  }

  if (df > 0 && abs(t2) > 0) {
    TLI <- 1 - t1 / t2
  } else if (!is.finite(t1) || !is.finite(t2)) {
    TLI <- as.numeric(NA)
  } else {
    TLI <- 1
  }

  TLI
}

# alias for nnfi
lav_fit_nnfi <- lav_fit_tli

# RFI - relative fit index (Bollen, 1986; Joreskog & Sorbom 1993)
lav_fit_rfi <- function(X2 = NULL, df = NULL, X2.null = NULL, df.null = NULL) {
  if (anyNA(c(X2, df, X2.null, df.null))) {
    return(as.numeric(NA))
  }

  if (df > df.null) {
    RLI <- as.numeric(NA)
  } else if (df > 0 && df.null > 0) {
    t1 <- X2.null / df.null - X2 / df
    t2 <- X2.null / df.null
    if (!is.finite(t1) || !is.finite(t2)) {
      RLI <- as.numeric(NA)
    } else if (t1 < 0 || t2 < 0) {
      RLI <- 1
    } else {
      RLI <- t1 / t2
    }
  } else {
    RLI <- 1
  }

  RLI
}

# NFI - normed fit index (Bentler & Bonett, 1980)
lav_fit_nfi <- function(X2 = NULL, df = NULL, X2.null = NULL, df.null = NULL) {
  if (anyNA(c(X2, df, X2.null, df.null))) {
    return(as.numeric(NA))
  }

  if (df > df.null || isTRUE(all.equal(X2.null, 0))) {
    NFI <- as.numeric(NA)
  } else if (df > 0) {
    t1 <- X2.null - X2
    t2 <- X2.null
    NFI <- t1 / t2
  } else {
    NFI <- 1
  }

  NFI
}

# PNFI - Parsimony normed fit index (James, Mulaik & Brett, 1982)
lav_fit_pnfi <- function(X2 = NULL, df = NULL, X2.null = NULL, df.null = NULL) {
  if (anyNA(c(X2, df, X2.null, df.null))) {
    return(as.numeric(NA))
  }

  if (df.null > 0 && X2.null > 0) {
    t1 <- X2.null - X2
    t2 <- X2.null
    PNFI <- (df / df.null) * (t1 / t2)
  } else {
    PNFI <- as.numeric(NA)
  }

  PNFI
}

# IFI - incremental fit index (Bollen, 1989; Joreskog & Sorbom, 1993)
lav_fit_ifi <- function(X2 = NULL, df = NULL, X2.null = NULL, df.null = NULL) {
  if (anyNA(c(X2, df, X2.null, df.null))) {
    return(as.numeric(NA))
  }

  t1 <- X2.null - X2
  t2 <- X2.null - df
  if (!is.finite(t1) || !is.finite(t2)) {
    IFI <- as.numeric(NA)
  } else if (t2 < 0) {
    IFI <- 1
  } else if (isTRUE(all.equal(t2, 0))) {
    IFI <- as.numeric(NA)
  } else {
    IFI <- t1 / t2
  }

  IFI
}

# higher-level function
lav_fit_cfi_lavobject <- function(lavobject = NULL, fit.measures = "cfi",
                                  baseline.model = NULL, h1.model = NULL,
                                  standard.test = "standard",
                                  scaled.test = "none",
                                  robust = TRUE,
                                  cat.check.pd = TRUE) {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # check for categorical
  categorical.flag <- lavobject@Model@categorical

  # tests
  TEST <- lavobject@test
  test.names <- sapply(lavobject@test, "[[", "test")
  if (test.names[1] == "none" || standard.test == "none") {
    return(list())
  }
  test.idx <- which(test.names == standard.test)[1]
  if (length(test.idx) == 0L) {
    return(list())
  }

  scaled.flag <- FALSE
  if (!scaled.test %in% c("none", "standard", "default")) {
    scaled.idx <- which(test.names == scaled.test)
    if (length(scaled.idx) > 0L) {
      scaled.idx <- scaled.idx[1] # only the first one
      scaled.flag <- TRUE
    }
  }

  # robust?
  robust.flag <- FALSE
  if (robust && scaled.flag &&
    scaled.test %in% c(
      "satorra.bentler", "yuan.bentler.mplus",
      "yuan.bentler", "scaled.shifted"
    )) {
    robust.flag <- TRUE
  }

  # FIML?
  fiml.flag <- FALSE
  if (robust && lavobject@Options$missing %in% c("ml", "ml.x")) {
    fiml.flag <- robust.flag <- TRUE
    # check if we can compute corrected values
    if (scaled.flag) {
      version <- "V3"
    } else {
      version <- "V6"
    }
    fiml <- try(lav_fit_fiml_corrected(lavobject, version = version),
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
  fit.baseline <- c("baseline.chisq", "baseline.df", "baseline.pvalue")
  if (scaled.flag) {
    fit.baseline <- c(
      fit.baseline, "baseline.chisq.scaled",
      "baseline.df.scaled", "baseline.pvalue.scaled",
      "baseline.chisq.scaling.factor"
    )
  }

  fit.cfi.tli <- c("cfi", "tli")
  if (scaled.flag) {
    fit.cfi.tli <- c(fit.cfi.tli, "cfi.scaled", "tli.scaled")
  }
  if (robust.flag) {
    fit.cfi.tli <- c(fit.cfi.tli, "cfi.robust", "tli.robust")
  }

  # other incremental fit indices
  fit.cfi.other <- c("nnfi", "rfi", "nfi", "pnfi", "ifi", "rni")
  if (scaled.flag) {
    fit.cfi.other <- c(
      fit.cfi.other, "nnfi.scaled", "rfi.scaled",
      "nfi.scaled", "pnfi.scaled", "ifi.scaled", "rni.scaled"
    )
  }
  if (robust.flag) {
    fit.cfi.other <- c(fit.cfi.other, "nnfi.robust", "rni.robust")
  }

  # which one do we need?
  if (missing(fit.measures)) {
    # default set
    fit.measures <- c(fit.baseline, fit.cfi.tli)
  } else {
    # remove any not-CFI related index from fit.measures
    rm.idx <- which(!fit.measures %in%
      c(fit.baseline, fit.cfi.tli, fit.cfi.other))
    if (length(rm.idx) > 0L) {
      fit.measures <- fit.measures[-rm.idx]
    }
    if (length(fit.measures) == 0L) {
      return(list())
    }
  }


  # basic test statistics
  X2 <- TEST[[test.idx]]$stat
  df <- TEST[[test.idx]]$df
  G <- lavobject@Data@ngroups # number of groups
  N <- lav_utils_get_ntotal(lavobject = lavobject) # N vs N-1

  # scaled X2
  if (scaled.flag) {
    X2.scaled <- TEST[[scaled.idx]]$stat
    df.scaled <- TEST[[scaled.idx]]$df
  }
  if (robust.flag) {
    XX3 <- X2
    if (categorical.flag) {
      out <- try(lav_fit_catml_dwls(lavobject, check.pd = cat.check.pd),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        XX3 <- df3 <- c.hat <- c.hat3 <- XX3.scaled <- as.numeric(NA)
      } else {
        XX3 <- out$XX3
        df3 <- out$df3
        c.hat3 <- c.hat <- out$c.hat3
        XX3.scaled <- out$XX3.scaled
      }
    } else if (fiml.flag) {
      XX3 <- fiml$XX3
      df3 <- fiml$df3
      c.hat3 <- c.hat <- fiml$c.hat3
      XX3.scaled <- fiml$XX3.scaled
    } else if (scaled.test == "scaled.shifted") {
      # compute c.hat from a and b
      a <- TEST[[scaled.idx]]$scaling.factor
      b <- TEST[[scaled.idx]]$shift.parameter
      c.hat <- a * (df - b) / df
    } else {
      c.hat <- TEST[[scaled.idx]]$scaling.factor
    }
  }

  # output container
  indices <- list()

  # only do what is needed (per groups)
  cfi.baseline.flag <- cfi.tli.flag <- cfi.other.flag <- FALSE
  if (any(fit.baseline %in% fit.measures)) {
    cfi.baseline.flag <- TRUE
  }
  if (any(fit.cfi.tli %in% fit.measures)) {
    cfi.tli.flag <- TRUE
  }
  if (any(fit.cfi.other %in% fit.measures)) {
    cfi.other.flag <- TRUE
  }

  # 1. BASELINE model
  baseline.test <- NULL

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
  if (!is.null(h1.model)) {
    stopifnot(inherits(h1.model, "lavaan"))

    # 2. h1 model in @external slot
  } else if (!is.null(lavobject@external$h1.model)) {
    stopifnot(inherits(lavobject@external$h1.model, "lavaan"))
    h1.model <- lavobject@external$h1.model
  } # else is.null

  # 1. user-provided baseline model
  if (!is.null(baseline.model)) {
    baseline.test <-
      lav_fit_measures_check_baseline(
        fit.indep = baseline.model,
        object = lavobject,
        fit.h1 = h1.model # okay if NULL
      )
    # 2. baseline model in @external slot
  } else if (!is.null(lavobject@external$baseline.model)) {
    fit.indep <- lavobject@external$baseline.model
    baseline.test <-
      lav_fit_measures_check_baseline(
        fit.indep = fit.indep,
        object = lavobject,
        fit.h1 = h1.model # okay if NULL
      )
    # 3. internal @baseline slot
  } else if (.hasSlot(lavobject, "baseline") &&
    length(lavobject@baseline) > 0L &&
    !is.null(lavobject@baseline$test) &&
    ## if there is a custom h1.model, need  _check_baseline() to update @test
    is.null(h1.model)) {
    baseline.test <- lavobject@baseline$test
    # 4. (re)compute independence model
  } else {
    fit.indep <- try(lav_object_independence(lavobject), silent = TRUE)
    baseline.test <-
      lav_fit_measures_check_baseline(
        fit.indep = fit.indep,
        object = lavobject,
        fit.h1 = h1.model # okay if NULL
      )
  }

  # baseline.test.idx
  baseline.test.idx <- which(names(baseline.test) == standard.test)[1]
  if (scaled.flag) {
    baseline.scaled.idx <- which(names(baseline.test) == scaled.test)[1]
  }

  if (!is.null(baseline.test)) {
    X2.null <- baseline.test[[baseline.test.idx]]$stat
    df.null <- baseline.test[[baseline.test.idx]]$df
    if (scaled.flag) {
      X2.null.scaled <- baseline.test[[baseline.scaled.idx]]$stat
      df.null.scaled <- baseline.test[[baseline.scaled.idx]]$df
    }
    if (robust.flag) {
      XX3.null <- X2.null
      if (categorical.flag) {
        if (inherits(out, "try-error")) {
          XX3.null <- c.hat.null <- as.numeric(NA)
        } else {
          XX3.null <- out$XX3.null
          c.hat.null <- out$c.hat3.null
        }
      } else if (fiml.flag) {
        XX3.null <- fiml$XX3.null
        c.hat.null <- fiml$c.hat3.null
      } else if (scaled.test == "scaled.shifted") {
        # compute c.hat from a and b
        a.null <-
          baseline.test[[baseline.scaled.idx]]$scaling.factor
        b.null <-
          baseline.test[[baseline.scaled.idx]]$shift.parameter
        c.hat.null <- a.null * (df.null - b.null) / df.null
      } else {
        c.hat.null <-
          baseline.test[[baseline.scaled.idx]]$scaling.factor
      }
    }
  } else {
    X2.null <- df.null <- as.numeric(NA)
    X2.null.scaled <- df.null.scaled <- as.numeric(NA)
    c.hat.null <- as.numeric(NA)
  }

  # check for NAs of nonfinite numbers
  if (!is.finite(X2) || !is.finite(df) ||
    !is.finite(X2.null) || !is.finite(df.null)) {
    indices[fit.measures] <- as.numeric(NA)
    return(indices)
  }

  # fill in baseline indices
  if (cfi.baseline.flag) {
    indices["baseline.chisq"] <- X2.null
    indices["baseline.df"] <- df.null
    indices["baseline.pvalue"] <- baseline.test[[baseline.test.idx]]$pvalue
    if (scaled.flag) {
      indices["baseline.chisq.scaled"] <- X2.null.scaled
      indices["baseline.df.scaled"] <- df.null.scaled
      indices["baseline.pvalue.scaled"] <-
        baseline.test[[baseline.scaled.idx]]$pvalue
      indices["baseline.chisq.scaling.factor"] <-
        baseline.test[[baseline.scaled.idx]]$scaling.factor
    }
  }

  # 2. CFI and TLI
  if (cfi.tli.flag) {
    indices["cfi"] <- lav_fit_cfi(
      X2 = X2, df = df,
      X2.null = X2.null, df.null = df.null
    )
    indices["tli"] <- lav_fit_tli(
      X2 = X2, df = df,
      X2.null = X2.null, df.null = df.null
    )
    if (scaled.flag) {
      indices["cfi.scaled"] <-
        lav_fit_cfi(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
      indices["tli.scaled"] <-
        lav_fit_tli(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
    }
    if (robust.flag) {
      indices["cfi.robust"] <-
        lav_fit_cfi(
          X2 = XX3, df = df,
          X2.null = XX3.null, df.null = df.null,
          c.hat = c.hat, c.hat.null = c.hat.null
        )
      indices["tli.robust"] <-
        lav_fit_tli(
          X2 = XX3, df = df,
          X2.null = XX3.null, df.null = df.null,
          c.hat = c.hat, c.hat.null = c.hat.null
        )
    }
  }

  # 3. other
  # c("nnfi", "rfi", "nfi", "pnfi", "ifi", "rni")
  if (cfi.other.flag) {
    indices["nnfi"] <-
      lav_fit_nnfi(X2 = X2, df = df, X2.null = X2.null, df.null = df.null)
    indices["rfi"] <-
      lav_fit_rfi(X2 = X2, df = df, X2.null = X2.null, df.null = df.null)
    indices["nfi"] <-
      lav_fit_nfi(X2 = X2, df = df, X2.null = X2.null, df.null = df.null)
    indices["pnfi"] <-
      lav_fit_pnfi(X2 = X2, df = df, X2.null = X2.null, df.null = df.null)
    indices["ifi"] <-
      lav_fit_ifi(X2 = X2, df = df, X2.null = X2.null, df.null = df.null)
    indices["rni"] <-
      lav_fit_rni(X2 = X2, df = df, X2.null = X2.null, df.null = df.null)

    if (scaled.flag) {
      indices["nnfi.scaled"] <-
        lav_fit_nnfi(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
      indices["rfi.scaled"] <-
        lav_fit_rfi(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
      indices["nfi.scaled"] <-
        lav_fit_nfi(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
      indices["pnfi.scaled"] <-
        lav_fit_pnfi(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
      indices["ifi.scaled"] <-
        lav_fit_ifi(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
      indices["rni.scaled"] <-
        lav_fit_rni(
          X2 = X2.scaled, df = df.scaled,
          X2.null = X2.null.scaled, df.null = df.null.scaled
        )
    }
    if (robust.flag) {
      indices["nnfi.robust"] <-
        lav_fit_nnfi(
          X2 = XX3, df = df,
          X2.null = XX3.null, df.null = df.null,
          c.hat = c.hat, c.hat.null = c.hat.null
        )
      indices["rni.robust"] <-
        lav_fit_rni(
          X2 = XX3, df = df,
          X2.null = XX3.null, df.null = df.null,
          c.hat = c.hat, c.hat.null = c.hat.null
        )
    }
  }

  # return only those that were requested
  indices[fit.measures]
}


# new in 0.6-5
# internal function to check the (external) baseline model, and
# return baseline 'test' list if everything checks out (and NULL otherwise)
lav_fit_measures_check_baseline <- function(fit.indep = NULL, object = NULL,
                                            fit.h1 = NULL) {
  TEST <- NULL

  # check if everything is in order
  if (inherits(fit.indep, "try-error")) {
    lav_msg_warn(gettext("baseline model estimation failed"))
    return(NULL)
  } else if (!inherits(fit.indep, "lavaan")) {
    lav_msg_warn(gettext(
      "(user-provided) baseline model is not a fitted lavaan object"))
    return(NULL)
  } else if (!fit.indep@optim$converged) {
    lav_msg_warn(gettext("baseline model did not converge"))
    return(NULL)
  } else {
    # evaluate if estimator/test matches original object
    # note: we do not need to check for 'se', as it may be 'none'
    sameTest <- all(object@Options$test == fit.indep@Options$test)
    if (!sameTest) {
      lav_msg_warn(gettextf(
        "Baseline model was using test(s) = %1$s, but original model was using
        test(s) = %2$s. Refitting baseline model!",
        lav_msg_view(fit.indep@Options$test, "none"),
        lav_msg_view(object@Options$test, "none")))
    }
    sameEstimator <- (object@Options$estimator ==
      fit.indep@Options$estimator)
    if (!sameEstimator) {
      lav_msg_warn(gettextf(
        "Baseline model was using estimator = %1$s, but original model was
        using estimator = %2$s. Refitting baseline model!",
        dQuote(fit.indep@Options$estimator),
        dQuote(object@Options$estimator)))
    }
    if (!sameTest || !sameEstimator) {
      lavoptions <- object@Options
      lavoptions$estimator <- object@Options$estimator
      lavoptions$se <- "none"
      lavoptions$verbose <- FALSE
      lavoptions$baseline <- FALSE
      lavoptions$check.start <- FALSE
      lavoptions$check.post <- FALSE
      lavoptions$check.vcov <- FALSE
      lavoptions$test <- object@Options$test
      fit.indep <- try(
        lavaan(fit.indep,
          slotOptions     = lavoptions,
          slotData        = object@Data,
          slotSampleStats = object@SampleStats,
          sloth1          = object@h1,
          slotCache       = object@Cache
        ),
        silent = TRUE
      )
      # try again
      TEST <- lav_fit_measures_check_baseline(
        fit.indep = fit.indep,
        object = object
      )
    } else {
      # extract what we need
      TEST <- fit.indep@test
    }
  } # converged lavaan object



  # TDJ: Check for user-supplied h1.model (here, the fit.h1= argument)
  #      Similar to BASELINE model, use the following priority:
  #        1. user-provided h1 model
  #        2. h1 model in @external slot
  #        3. default h1 model (already in @h1 slot, no update necessary)
  #FIXME? user-supplied h1 model in object might be in fit.indep, too

  user_h1_exists <- FALSE
  # 1. user-provided h1 model
  if (!is.null(fit.h1)) {
    stopifnot(inherits(fit.h1, "lavaan"))
    user_h1_exists <- TRUE

    # 2. h1 model in @external slot
  } else if (!is.null(object@external$h1.model)) {
    stopifnot(inherits(object@external$h1.model, "lavaan"))
    fit.h1 <- object@external$h1.model
    user_h1_exists <- TRUE
  }

  if (user_h1_exists) {
    ## update @test slot
    TEST <- lav_update_test_custom_h1(lav_obj_h0 = fit.indep,
                                      lav_obj_h1 = fit.h1)@test
  }

  TEST
}
