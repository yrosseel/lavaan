# user-visible function to extract the fit measures
# output can be 1) vector (default), 2) list, 3) matrix, or 4) text
# in the latter case, the result will be of class "lavaan.fitMeasures"
# for which the printing is done by lav_fitmeasures_print()

# new in 0.6-13:
# the big families are computed in dedicated functions:
# - lav_fit_rmsea_lavobject
# - lav_fit_cfi_lavobject
# - lav_fit_aic_lavojbect
# - lav_residuals_summary

# Note: fitMeasures/fitmeasures are generic functions; they include a "..."
#       so lavaan.mi can add arguments to pass to lavTestLRT() and
#       lavTestLRT.mi() about how to pool chi-squared.

setMethod(
  "fitMeasures", signature(object = "lavaan"),
  function(object, fit_measures = "all", baseline_model = NULL, h1_model = NULL,
           fm_args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.nonpd = "na"
           ),
           output = "vector", ...) {
    dotdotdot <- list(...)
    lav_adapt_func(environment(), dotdotdot, NULL)
    if (!is.list(fit_measures))
        fit_measures <- list(fit.measures = fit_measures)
    if (!missing(fm_args)) fit_measures <- c(fit_measures, fm_args)
    lav_fit(
      object = object, fit_measures = fit_measures,
      baseline_model = baseline_model, h1_model = h1_model,
      output = output
    )
  }
)

setMethod(
  "fitmeasures", signature(object = "lavaan"),
  function(object, fit_measures = "all", baseline_model = NULL, h1_model = NULL,
           fm_args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.nonpd = "na"
           ),
           output = "vector", ...) {
    dotdotdot <- list(...)
    lav_adapt_func(environment(), dotdotdot, NULL)
    if (!is.list(fit_measures))
                  fit_measures <- list(fit.measures = fit_measures)
    if (!missing(fm_args)) fit_measures <- c(fit_measures, fm_args)
    lav_fit(
      object = object, fit_measures = fit_measures,
      baseline_model = baseline_model, h1_model = h1_model,
      output = output
    )
  }
)

# S3 method for efaList
lav_efalist_fitmeasures <- function(
    object,
    fit_measures = "all",
    baseline_model = NULL, h1_model = NULL,
    fm_args = list(
      standard.test = "default",
      scaled.test = "default",
      rmsea.ci.level = 0.90,
      rmsea.close.h0 = 0.05,
      rmsea.notclose.h0 = 0.08,
      robust = TRUE,
      cat.nonpd = "na"
    ),
    output = "list", ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, FALSE)

  # kill object$loadings if present
  object[["loadings"]] <- NULL

  # get fit measures for each model
  if (!is.list(fit_measures)) fit_measures <- list(fit.measures = fit_measures)
  if (!missing(fm_args)) fit_measures <- c(fit_measures, fm_args)
  res <- simplify2array(lapply(
    object,
    function(x) {
      lav_fit(
        object = x,
        fit_measures = fit_measures, h1_model = h1_model,
        baseline_model = baseline_model,
        output = "vector"
      )
    }
  )
  )

  # check if res is a matrix
  if (!is.matrix(res)) {
    if (is.numeric(res)) {
      # fit_measures is just 1 element, or only one was correct
      name <- names(res)[1]
      res <- matrix(res, nrow = 1L)
      rownames(res) <- name
    } else { # wrong fit measures?
      # create empty matrix
      res <- matrix(0, nrow = 0L, ncol = length(object))
    }
  }

  # rownames
  nfactors <- sapply(object, function(x) x@pta$nfac[[1]])
  colnames(res) <- paste0("nfactors = ", nfactors)

  # class
  class(res) <- c("lavaan.matrix", "matrix")

  res
}


lav_fit <- function(object, fit_measures = "all",
                             baseline_model = NULL, h1_model = NULL,
                             fm_args = list(
                               standard.test = "default",
                               scaled.test = "default",
                               rmsea.ci.level = 0.90,
                               rmsea.close.h0 = 0.05,
                               rmsea.notclose.h0 = 0.08,
                               robust = TRUE,
                               cat.nonpd = "na"
                             ),
                             output = "vector") {
  # check object
  object <- lav_object_check_version(object)

  # default fm_args
  default_fm_args <- list(
    standard.test = "default",
    scaled.test = "default",
    rmsea.ci.level = 0.90,
    rmsea.close.h0 = 0.05,
    rmsea.notclose.h0 = 0.08,
    robust = TRUE,
    cat.nonpd = "na"
  )
  if (!missing(fm_args)) {
    lav_deprecated_args("fit_measures", "fm_args")
    fm_args <- modifyList(default_fm_args, fm_args)
  } else {
    fm_args <- default_fm_args
  }
  if (is.list(fit_measures)) {
    if (is.null(names(fit_measures)) ||
        is.null(fit_measures$fit.measures)) {
      lav_msg_stop(gettextf(
        "If %s is a list, it must contain a named element %s.",
        "fit_measures", "fit.measures"
      ))
    }
    temp <- fit_measures$fit.measures
    fit_measures$fit.measures <- NULL
    fm_args <- modifyList(default_fm_args, fit_measures)
    fit_measures <- temp
  }

  # resolve the cat.nonpd option (new in 0.7-1): what to do when data are
  # categorical and an input correlation matrix is not positive definite?
  #    "na"     -> robust fit indices are NA (default)
  #    "refit"  -> smooth the matrix, and re-estimate the parameters
  #                using estimator = "catML"
  #    "smooth" -> smooth the matrix, but keep the (DWLS) parameter
  #                estimates (as in lavaan <= 0.6-13)
  # the older cat.check.pd = FALSE is a deprecated alias for "refit"
  if (!is.null(fm_args$cat.check.pd) && !isTRUE(fm_args$cat.check.pd) &&
      identical(fm_args$cat.nonpd, "na")) {
    fm_args$cat.nonpd <- "refit"
  }
  if (!is.character(fm_args$cat.nonpd) ||
      !fm_args$cat.nonpd %in% c("na", "refit", "smooth")) {
    lav_msg_warn(gettextf(
      "invalid cat.nonpd value [%s] set to default \"na\".",
      fm_args$cat.nonpd
    ))
    fm_args$cat.nonpd <- "na"
  }

  # standard test
  if (fm_args$standard.test == "default") {
    fm_args$standard.test <- object@Options$standard.test
    # usually "standard", but could have been changed
    if (is.null(fm_args$standard.test)) { # <older objects
      fm_args$standard.test <- "standard"
    }
  }

  # scaled test
  if (fm_args$scaled.test == "default") {
    fm_args$scaled.test <- object@Options$scaled.test
    # usually "standard", but could have been changed
    if (is.null(fm_args$scaled.test)) { # <older objects
      fm_args$scaled.test <- "standard"
    }
  }

  # do we have data? (yep, we had to include this check)
  if (object@Data@data.type == "none") {
    lav_msg_stop(gettext("fit measures not available if there is no data."))
  }

  # has the model converged?
  if (object@optim$npar > 0L && !object@optim$converged) {
    lav_msg_stop(gettext(
      "fit measures not available if model did not converge"))
  }

  # random slopes (rv() modifier): no test statistics (and no
  # comparable h1 loglikelihood) exist, but the loglikelihood-based
  # measures -- and the H0 scaling correction factor -- are available
  if (length(object@Model@rv.ov) > 0L ||
      length(object@Model@rv.lv) > 0L) {
    return(lav_fit_rv(object, fit_measures = fit_measures,
                      output = output))
  }

  # do we have a test statistic?
  test <- lavInspect(object, "test")
  test_names <- unname(sapply(test, "[[", "test"))
  if (test_names[1] == "none") {
    lav_msg_stop(gettext("fit measures not available if test = \"none\"."))
    #FIXME: allow RMRs, log.likelihoods, info criteria, npar, ntotal
  }

  standard_test <- fm_args$standard.test
  scaled_test <- fm_args$scaled.test

  # check standard.test
  standard_test <- lav_test_rename(standard_test, check = TRUE)[1] # only 1
  fmg_standard_test <- lav_test_fmg_is_fmg(standard_test)

  # check scaled.test
  if (!scaled_test %in% c("none", "default", "standard")) {
    scaled_test <- lav_test_rename(scaled_test, check = TRUE)[1] # only 1
  }

  # which test statistic do we need?
  rerun_lavtest_flag <- FALSE
  if (!standard_test %in% test_names) {
    rerun_lavtest_flag <- TRUE
  }
  if (!scaled_test %in% c("none", "default", "standard") &&
      !scaled_test %in% test_names) {
    rerun_lavtest_flag <- TRUE
  }

  # do we have a scaled test statistic? if so, which one?
  scaled_flag <- FALSE
  if (!fmg_standard_test &&
      scaled_test != "none" &&
      any(test_names %in% c(
        "satorra.bentler",
        "yuan.bentler", "yuan.bentler.mplus", "yuan.chan",
        "mean.var.adjusted", "scaled.shifted"
      ))) {
    scaled_flag <- TRUE
    if (scaled_test %in% c("standard", "default")) {
      tmp_idx <- which(test_names %in% c(
        "satorra.bentler",
        "yuan.bentler", "yuan.bentler.mplus", "yuan.chan",
        "mean.var.adjusted", "scaled.shifted"
      ))
      scaled_test <- test_names[tmp_idx[1]]
    }
  }

  # rerun lavTest?
  if (rerun_lavtest_flag) {
    this_test <- standard_test
    if (scaled_flag) {
      this_test <- unique(this_test, scaled_test)
    }
    lavtest_scaled_test <- standard_test
    if (fmg_standard_test) {
      lavtest_scaled_test <- scaled_test
    }
    test <- lavTest(object,
                    test = this_test, scaled_test = lavtest_scaled_test,
                    drop_list_single = FALSE
    )
    # replace in object, if we pass it to lav_fit_* functions
    object@test <- test
    test_names <- unname(sapply(test, "[[", "test"))
  }


  # TDJ: Check for user-supplied h1 model
  #      Similar to BASELINE model, use the following priority:
  #        1. user-provided h1 model
  #        2. h1 model in @external slot
  #        3. default h1 model (already in @h1 slot, no update necessary)

  user_h1_exists <- FALSE
  # 1. user-provided h1 model
  if (!is.null(h1_model)) {
    stopifnot(inherits(h1_model, "lavaan"))
    user_h1_exists <- TRUE

    # 2. h1 model in @external slot
  } else if (!is.null(object@external$h1.model)) {
    stopifnot(inherits(object@external$h1.model, "lavaan"))
    h1_model <- object@external$h1.model
    user_h1_exists <- TRUE
  }

  ## Update statistics in @test slot?
  if (user_h1_exists) {
    ## update @test slot
    fit <- lav_update_test_custom_h1(lav_obj_h0 = object, lav_obj_h1 = h1_model)

    ## re-assign TEST object that is used below
    object@test <- test <- fit@test
    test_names <- unname(sapply(test, "[[", "test"))
  }

  # get index of standard.test in TEST
  test_idx <- which(test_names == standard_test)[1]

  # get index of scaled test (if any) in TEST
  if (scaled_flag) {
    scaled_idx <- which(test_names == scaled_test)[1]
  }

  # check output argument
  if (output %in% c("vector", "horizontal")) {
    output <- "vector"
  } else if (output %in% c("list")) {
    output <- "list"
  } else if (output %in% c("matrix", "vertical")) {
    output <- "matrix"
  } else if (output %in% c("text", "pretty", "summary")) {
    output <- "text"
  } else {
    lav_msg_stop(gettextf("output should be %s.",
                          lav_msg_view(c("vector", "list", "matrix", "text"),
                          "none", FALSE)
    ))
  }

  # options
  categorical_flag <- object@Model@categorical
  fiml_flag <- (fm_args$robust &&
                  object@Options$missing %in% c("ml", "ml.x"))
  estimator <- object@Options$estimator

  # basic ingredients
  x2 <- test[[test_idx]]$stat
  df <- test[[test_idx]]$df
  if (scaled_flag) {
    x2_scaled <- test[[scaled_idx]]$stat
    df_scaled <- test[[scaled_idx]]$df
  }
  npar <- lav_inspect_npar(object = object, ceq = TRUE)
  n <- lav_inspect_ntotal(object = object) # N vs N-1


  # define 'sets' of fit measures:
  fit_always <- c("npar")

  # basic chi-square test
  fit_chisq <- c("fmin", "chisq", "df", "pvalue")
  if (scaled_flag) {
    fit_chisq <- c(
      fit_chisq, "chisq.scaled", "df.scaled", "pvalue.scaled",
      "chisq.scaling.factor"
    )
  }

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
  if (fm_args$robust && (scaled_flag || categorical_flag || fiml_flag)) {
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
  if (fm_args$robust && (scaled_flag || categorical_flag || fiml_flag)) {
    fit_cfi_other <- c(fit_cfi_other, "nnfi.robust", "rni.robust")
  }
  fit_cfi <- c(fit_baseline, fit_cfi_tli, fit_cfi_other)

  # likelihood based measures
  if (estimator == "MML") {
    fit_logl <- c("logl", "aic", "bic", "ntotal", "bic2")
  } else {
    fit_logl <- c(
      "logl", "unrestricted.logl", "aic", "bic",
      "ntotal", "bic2"
    )
  }
  if (scaled_flag &&
      scaled_test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
    fit_logl <- c(fit_logl, "scaling.factor.h1", "scaling.factor.h0")
  }

  # rmsea
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
  if (fm_args$robust && (scaled_flag || categorical_flag || fiml_flag)) {
    fit_rmsea <- c(
      fit_rmsea, "rmsea.robust", "rmsea.ci.lower.robust",
      "rmsea.ci.upper.robust", "rmsea.pvalue.robust",
      "rmsea.notclose.pvalue.robust"
    )
  }

  # gfi (the 'new' GFI of Maydeu-Olivares et al., 2024, with conf. interval)
  # note: the (printed) 'gfi' is the RLS-based version
  fit_gfi <- c("gfi", "gfi.ci.lower", "gfi.ci.upper", "gfi.ci.level")
  if (fm_args$robust && (scaled_flag || categorical_flag || fiml_flag)) {
    fit_gfi <- c(
      fit_gfi, "gfi.robust",
      "gfi.ci.lower.robust", "gfi.ci.upper.robust"
    )
  }
  # explicit LR-based (gfi_lrt) and RLS-based (gfi_rls) versions; these are not
  # part of the 'default' or 'all' sets, but can be requested by name
  fit_gfi_extra <- c(
    "gfi_lrt", "gfi_lrt.ci.lower", "gfi_lrt.ci.upper",
    "gfi_rls", "gfi_rls.ci.lower", "gfi_rls.ci.upper"
  )
  if (fm_args$robust && (scaled_flag || categorical_flag || fiml_flag)) {
    fit_gfi_extra <- c(
      fit_gfi_extra,
      "gfi_lrt.robust", "gfi_lrt.ci.lower.robust", "gfi_lrt.ci.upper.robust",
      "gfi_rls.robust", "gfi_rls.ci.lower.robust", "gfi_rls.ci.upper.robust"
    )
  }
  # the 'old' (LISREL) GFI and friends
  fit_gfi_lisrel <- c("gfi_lisrel", "agfi_lisrel", "pgfi")

  # srmr
  if (categorical_flag) {
    fit_srmr <- c("srmr")
    fit_srmr2 <- c(
      "rmr", "rmr_nomean",
      "srmr", # per default equal to srmr_bentler_nomean
      "srmr_bentler", "srmr_bentler_nomean",
      "crmr", "crmr_nomean",
      "srmr_mplus", "srmr_mplus_nomean"
    )
  } else {
    if (object@Data@nlevels > 1L) {
      fit_srmr <- c("srmr", "srmr_within", "srmr_between")
      fit_srmr2 <- c("srmr", "srmr_within", "srmr_between")
    } else {
      fit_srmr <- c("srmr")
      fit_srmr2 <- c(
        "rmr", "rmr_nomean",
        "srmr", # the default
        "srmr_bentler", "srmr_bentler_nomean",
        "crmr", "crmr_nomean",
        "srmr_mplus", "srmr_mplus_nomean"
      )
    }
  }

  # various
  if (object@Data@nlevels > 1L) {
    fit_other <- ""
  } else if (estimator == "PML") {
    fit_other <- c("cn_05", "cn_01", "mfi")
    if (!categorical_flag) { # really needed?
      fit_other <- c(fit_other, "ecvi")
    }
  } else {
    fit_other <- c("cn_05", "cn_01", fit_gfi_lisrel, "mfi")
    if (!categorical_flag) { # really needed?
      fit_other <- c(fit_other, "ecvi")
    } else {
      fit_other <- c(fit_other, "wrmr")
    }
  }

  # lower case
  fit_measures <- fit_measures_orig <- tolower(fit_measures)

  # select 'default' fit measures
  if (length(fit_measures) == 1L) {
    if (fit_measures == "default") {
      if (estimator == "ML" || estimator == "PML") {
        fit_measures <- c(
          fit_always, fit_chisq, fit_baseline,
          fit_cfi_tli, fit_logl,
          fit_rmsea, fit_srmr, fit_gfi
        )
      } else if (estimator == "MML") {
        fit_measures <- c(fit_always, fit_logl)
      } else {
        fit_measures <- c(
          fit_always,
          fit_chisq, fit_baseline, fit_cfi_tli,
          fit_rmsea, fit_srmr, fit_gfi
        )
      }
    } else if (fit_measures == "all") {
      if (estimator == "ML") {
        fit_measures <- c(
          fit_always, fit_chisq,
          fit_baseline, fit_cfi_tli, fit_cfi_other,
          fit_logl, fit_rmsea, fit_srmr2, fit_gfi, fit_other
        )
      } else {
        fit_measures <- c(
          fit_always, fit_chisq,
          fit_baseline, fit_cfi_tli, fit_cfi_other,
          fit_rmsea, fit_srmr2, fit_gfi, fit_other
        )
      }
    }
  }

  # catch empty list
  if (length(fit_measures) == 0L) {
    return(list())
  }

  # main container
  indices <- list()
  indices["npar"] <- npar
  indices["ntotal"] <- object@SampleStats@ntotal
  indices["fmin"] <- object@optim$fx # note = 0.5 * fmin if ML

  # CHI-SQUARE TEST
  if (any(fit_chisq %in% fit_measures)) {
    indices["chisq"] <- x2
    indices["df"] <- df
    indices["pvalue"] <- test[[test_idx]]$pvalue
    if (scaled_flag) {
      indices["chisq.scaled"] <- x2_scaled
      indices["df.scaled"] <- df_scaled
      indices["chisq.scaling.factor"] <- test[[scaled_idx]]$scaling.factor
      indices["pvalue.scaled"] <- test[[scaled_idx]]$pvalue
    }
  }

  # BASELINE FAMILY
  if (any(fit_cfi %in% fit_measures)) {

    # rerun baseline?
    if (rerun_lavtest_flag) {
      object@Options$test <- this_test
      fit_indep <- try(lav_object_independence(object), silent = TRUE)
      # override object
      object@baseline$test <- fit_indep@test
    }

    indices <- c(
      indices,
      lav_fit_cfi_lavobject(
        lavobject = object,
        fit_measures = fit_measures,
        baseline_model = baseline_model,
        h1_model = h1_model,
        standard_test = standard_test,
        scaled_test = scaled_test,
        robust = fm_args$robust,
        cat_nonpd = fm_args$cat.nonpd
      )
    )
  }

  # INFORMATION CRITERIA
  if (any(fit_logl %in% fit_measures)) {
    indices <- c(
      indices,
      lav_fit_aic_lavobject(
        lavobject = object,
        fit_measures = fit_measures,
        standard_test = standard_test,
        scaled_test = scaled_test,
        estimator = estimator
      )
    )
  }

  # RMSEA and friends
  if (any(fit_rmsea %in% fit_measures)) {
    # check rmsea options
    rmsea_ci_level <- 0.90
    rmsea_close_h0 <- 0.05
    rmsea_notclose_h0 <- 0.08
    if (!is.null(fm_args$rmsea.ci.level) &&
        is.finite(fm_args$rmsea.ci.level)) {
      rmsea_ci_level <- fm_args$rmsea.ci.level
      if (rmsea_ci_level < 0 || rmsea_ci_level > 1.0) {
        lav_msg_warn(gettextf(
          "invalid rmsea.ci.level value [%s] set to default 0.90.",
          rmsea_ci_level))
        rmsea_ci_level <- 0.90
      }
    }
    if (!is.null(fm_args$rmsea.close.h0) &&
        is.finite(fm_args$rmsea.close.h0)) {
      rmsea_close_h0 <- fm_args$rmsea.close.h0
      if (rmsea_close_h0 < 0) {
        rmsea_close_h0 <- 0
      }
    }
    if (!is.null(fm_args$rmsea.notclose.h0) &&
        is.finite(fm_args$rmsea.notclose.h0)) {
      rmsea_notclose_h0 <- fm_args$rmsea.notclose.h0
      if (rmsea_notclose_h0 < 0) {
        rmsea_notclose_h0 <- 0
      }
    }
    indices <- c(
      indices,
      lav_fit_rmsea_lavobject(
        lavobject = object,
        fit_measures = fit_measures,
        standard_test = standard_test,
        scaled_test = scaled_test,
        ci_level = rmsea_ci_level,
        close_h0 = rmsea_close_h0,
        notclose_h0 = rmsea_notclose_h0,
        robust = fm_args$robust,
        cat_nonpd = fm_args$cat.nonpd
      )
    )
  }

  # SRMR and friends
  if (any(fit_srmr2 %in% fit_measures)) {
    indices <- c(
      indices,
      lav_fit_srmr_lavobject(
        lavobject = object,
        fit_measures = fit_measures
      )
    )
  }

  # GFI and friends
  if (any(c(fit_gfi, fit_gfi_extra, fit_gfi_lisrel) %in% fit_measures)) {
    # check gfi.ci.level option (default 0.90)
    gfi_ci_level <- 0.90
    if (!is.null(fm_args$gfi.ci.level) &&
        is.finite(fm_args$gfi.ci.level)) {
      gfi_ci_level <- fm_args$gfi.ci.level
      if (gfi_ci_level < 0 || gfi_ci_level > 1.0) {
        lav_msg_warn(gettextf(
          "invalid gfi.ci.level value [%s] set to default 0.90.",
          gfi_ci_level))
        gfi_ci_level <- 0.90
      }
    }
    indices <- c(
      indices,
      lav_fit_gfi_lavobject(
        lavobject = object,
        fit_measures = fit_measures,
        standard_test = standard_test,
        scaled_test = scaled_test,
        ci_level = gfi_ci_level,
        robust = fm_args$robust,
        cat_nonpd = fm_args$cat.nonpd
      )
    )
  }

  # various: Hoelter Critical N (CN)
  if (any(c("cn_05", "cn_01") %in% fit_measures)) {
    indices["cn_05"] <- lav_fit_cn(x2 = x2, df = df, n = n, alpha = 0.05)
    indices["cn_01"] <- lav_fit_cn(x2 = x2, df = df, n = n, alpha = 0.01)
  }

  # various: WRMR
  if ("wrmr" %in% fit_measures) {
    nel <- length(object@SampleStats@WLS.obs[[1]])
    indices["wrmr"] <- lav_fit_wrmr(x2 = x2, nel = nel)
  }

  # various: MFI
  if ("mfi" %in% fit_measures) {
    indices["mfi"] <- lav_fit_mfi(x2 = x2, df = df, n = n)
  }

  # various: ECVI
  if ("ecvi" %in% fit_measures) {
    indices["ecvi"] <- lav_fit_ecvi(x2 = x2, npar = npar, n = n)
  }


  # keep only what we need
  out <- indices[fit_measures]

  if (all(is.na(names(out)))) { # perhaps, fit_measures = ""
    # nothing left
    return(numeric(0L))
  }

  # select output type
  if (output == "list") {
    # nothing to do
  } else if (output == "vector") {
    out <- unlist(out)
    class(out) <- c("lavaan.vector", "numeric")
  } else if (output == "matrix") {
    out <- as.matrix(unlist(out))
    colnames(out) <- ""
    class(out) <- c("lavaan.matrix", "matrix")
  } else if (output == "text") {
    out <- unlist(out)
    class(out) <- c("lavaan.fitMeasures", "lavaan.vector", "numeric")
  }

  # attributes?
  # only if fit_measures == "all" or "default"
  if (length(fit_measures_orig) == 1L &&
      fit_measures_orig %in% c("all", "default")) {
    x2_label <- test[[test_idx]]$label # NULL if "standard"
    x2_baseline_label <- object@baseline$test[[test_idx]]$label
    attr(out, "X2.label") <- x2_label
    attr(out, "X2.baseline.label") <- x2_baseline_label
    if (standard_test != "standard") {
      attr(out, "standard.test") <- standard_test
    }
    if (scaled_flag && scaled_test != "standard") {
      attr(out, "scaled.test") <- scaled_test
    }
  }

  out
}

# fit measures for models with random slopes (rv() modifier)
#
# no chi-square test statistic (or fit indices) exist for these
# models (the likelihood is conditional on the covariates carrying
# random slopes: there is no comparable saturated model); we provide
# the loglikelihood-based measures, and the H0 scaling correction
# factor (cfr. the Mplus MLR output) -- the latter is included in
# the default set only when se = "robust.huber.white", but can
# always be requested explicitly
lav_fit_rv <- function(object, fit_measures = "all",
                       output = "vector") {
  loglik <- object@loglik
  robust_flag <- object@Options$se == "robust.huber.white"

  indices <- list()
  indices["npar"] <- loglik$npar
  indices["ntotal"] <- loglik$ntotal
  indices["logl"] <- loglik$loglik
  indices["aic"] <- loglik$AIC
  indices["bic"] <- loglik$BIC
  indices["bic2"] <- loglik$BIC2

  available <- c(names(indices), "scaling.factor.h0")

  if (length(fit_measures) == 1L &&
      fit_measures %in% c("all", "default")) {
    fit_measures <- names(indices)
    if (robust_flag) {
      fit_measures <- c(fit_measures, "scaling.factor.h0")
    }
  } else {
    fit_measures <- tolower(fit_measures)
    bad <- fit_measures[!fit_measures %in% available]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "fit measure(s) not available for models with random slopes:
         %s.", paste(bad, collapse = " ")))
    }
  }

  # H0 scaling correction factor: c = tr(H^{-1} J) / npar
  if ("scaling.factor.h0" %in% fit_measures) {
    rs <- object@Cache[[1]]$rs
    if (is.null(rs)) {
      rs_info <- lav_mvn_cl_rs_info(
        lavmodel = object@Model, lavdata = object@Data
      )
      rs <- list(
        info = rs_info,
        stats = lav_mvn_cl_rs_stats(
          y1 = object@Data@X[[1]], lp = object@Data@Lp[[1]],
          rs_info = rs_info
        )
      )
    }
    indices["scaling.factor.h0"] <-
      lav_mvn_cl_rs_scaling_h0(lavobject = object, rs = rs)
  }

  # keep only what we need (in the requested order)
  out <- indices[fit_measures]

  # select output type
  if (output == "list") {
    # nothing to do
  } else if (output == "vector") {
    out <- unlist(out)
    class(out) <- c("lavaan.vector", "numeric")
  } else if (output == "matrix") {
    out <- as.matrix(unlist(out))
    colnames(out) <- ""
    class(out) <- c("lavaan.matrix", "matrix")
  } else if (output == "text") {
    out <- unlist(out)
    class(out) <- c("lavaan.fitMeasures", "lavaan.vector", "numeric")
  }

  out
}

# print a nice summary of the fit measures
lav_fitmeasures_print <- function(x, ..., nd = 3L, add_h0 = TRUE) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, FALSE)

  names_x <- names(x)

  # scaled?
  scaled_flag <- "chisq.scaled" %in% names_x

  # num format
  num_format <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

  ## TDJ: optionally add h0 model's fit statistic, for lavaan.mi
  if (add_h0 && "chisq" %in% names_x) {
    cat("\nModel Test User Model:\n\n")

    # container three columns
    c1 <- c2 <- c3 <- character(0L)

    # TDJ: Add header used in summary() by lavaan.mi
    if (scaled_flag) {
      c1 <- c("", c1)
      c2 <- c("Standard", c2)
      c3 <- c("Scaled", c3)
    }

    # if test is not standard, add label
    if (!is.null(attr(x, "X2.label"))) {
      c1 <- c(c1, attr(x, "X2.label"))
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }

    c1 <- c(c1, "Test statistic")
    c2 <- c(c2, sprintf(num_format, x["chisq"]))
    c3 <- c(c3, ifelse(scaled_flag,
                       sprintf(num_format, x["chisq.scaled"]), ""
    ))

    c1 <- c(c1, "Degrees of freedom")
    c2 <- c(c2, x["df"])
    c3 <- c(c3, ifelse(scaled_flag,
                       ifelse(x["df.scaled"] %% 1 == 0, x["df.scaled"],
                              sprintf(num_format, x["df.scaled"])
                       ), ""
    ))

    c1 <- c(c1, "P-value")
    c2 <- c(c2, sprintf(num_format, x["pvalue"]))
    c3 <- c(c3, ifelse(scaled_flag,
                       sprintf(num_format, x["pvalue.scaled"]), ""
    ))

    if (scaled_flag && "chisq.scaling.factor" %in% names_x) {
      c1 <- c(c1, "Average scaling correction factor")
      c2 <- c(c2, "")
      c3 <- c(c3, sprintf(num_format, x["chisq.scaling.factor"]))

      ## check for shift parameter
      chisq_shift_parameter <- attr(attr(x, "header"), "shift")
      if (!is.null(chisq_shift_parameter)) {
        c1 <- c(c1, "Average shift parameter",
                "  simple second-order correction")
        c2 <- c(c2, "", "")
        ## This is only provided by the fitMeasures() method for lavaan.mi-class
        c3 <- c(c3, sprintf(num_format, chisq_shift_parameter),
                "")
      }
    }

    # check for lavaan.mi attributes "pool.method" and "pool.robust"
    if (!is.null(attr(x, "pool.method"))) {
      ## extract information from lavaan.mi object about
      ## the method used to pool the test statistic
      pool_method   <- attr(x, "pool.method")
      pool_robust   <- attr(x, "pool.robust")
      standard_test <- attr(attr(x, "header"), "standard.test")
      scaled_test   <- attr(attr(x, "header"),   "scaled.test")

      c1 <- c(c1, "Pooling method")
      c2 <- c(c2, pool_method)
      c3 <- c(c3, "")

      ## (conditionally for D2 method) add other pooling information
      if (scaled_flag) {
        c1 <- c(c1, "  Pooled statistic")
        c2 <- c(c2, ifelse(pool_robust, dQuote(scaled_test),
                                        dQuote(standard_test)))
        c3 <- c(c3, "")

        if (pool_robust && pool_method == "D2") {
          c1 <- c(c1, paste0("  ", dQuote(scaled_test), " correction applied"))
          c2 <- c(c2, "BEFORE")
          c3 <- c(c3, "pooling")
        } else {
          c1 <- c(c1, paste0("  ", dQuote(scaled_test), " correction applied"))
          c2 <- c(c2, "AFTER")
          c3 <- c(c3, "pooling")
        }

      }
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 35L)
    c2 <- format(c2, width = 16L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # print information about standard.test? (new in 0.6-21)
  # (only if "standard.test" is not "standard")
  if (!is.null(attr(x, "standard.test"))) {
    cat("\nNote: fit measures based on the chi-square test statistic\n",
        "     use", attr(x, "X2.label"), "\n")
  }

  # independence model
  if ("baseline.chisq" %in% names_x) {
    cat("\nModel Test Baseline Model:\n\n")

    c1 <- c2 <- c3 <- character(0L)
    # if test is not standard, add label
    if (!is.null(attr(x, "X2.baseline.label"))) {
      c1 <- c(c1, attr(x, "X2.label"))
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }

    c1 <- c(c1, "Test statistic")
    c2 <- c(c2, sprintf(num_format, x["baseline.chisq"]))
    c3 <- c(c3, ifelse(scaled_flag,
                       sprintf(num_format, x["baseline.chisq.scaled"]), ""
    ))

    c1 <- c(c1, "Degrees of freedom")
    c2 <- c(c2, x["baseline.df"])
    c3 <- c(c3, ifelse(scaled_flag,
                  ifelse(x["baseline.df.scaled"] %% 1 == 0,
                         x["baseline.df.scaled"],
                         sprintf(num_format, x["baseline.df.scaled"])
    ),
    ""
    ))

    c1 <- c(c1, "P-value")
    c2 <- c(c2, sprintf(num_format, x["baseline.pvalue"]))
    c3 <- c(c3, ifelse(scaled_flag,
                       sprintf(num_format, x["baseline.pvalue.scaled"]), ""
    ))

    if (scaled_flag && "baseline.chisq.scaling.factor" %in% names_x) {
      c1 <- c(c1, "Scaling correction factor")
      c2 <- c(c2, "")
      c3 <- c(c3, sprintf(num_format, x["baseline.chisq.scaling.factor"]))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 35L)
    c2 <- format(c2, width = 16L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # cfi/tli
  if (any(c("cfi", "tli", "nnfi", "rfi", "nfi", "ifi", "rni", "pnfi") %in%
                                                                 names_x)) {
    cat("\nUser Model versus Baseline Model:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("cfi" %in% names_x) {
      c1 <- c(c1, "Comparative Fit Index (CFI)")
      c2 <- c(c2, sprintf(num_format, x["cfi"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["cfi.scaled"]), ""
      ))
    }

    if ("tli" %in% names_x) {
      c1 <- c(c1, "Tucker-Lewis Index (TLI)")
      c2 <- c(c2, sprintf(num_format, x["tli"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["tli.scaled"]), ""
      ))
    }

    if ("cfi.robust" %in% names_x) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
      c1 <- c(c1, "Robust Comparative Fit Index (CFI)")
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["cfi.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["cfi.robust"]))
        c3 <- c(c3, "")
      }
    }

    if ("tli.robust" %in% names_x) {
      c1 <- c(c1, "Robust Tucker-Lewis Index (TLI)")
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["tli.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["tli.robust"]))
        c3 <- c(c3, "")
      }
    }

    if ("nnfi" %in% names_x) {
      c1 <- c(c1, "Bentler-Bonett Non-normed Fit Index (NNFI)")
      c2 <- c(c2, sprintf(num_format, x["nnfi"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["nnfi.robust"]), ""
      ))
    }

    if ("nfi" %in% names_x) {
      c1 <- c(c1, "Bentler-Bonett Normed Fit Index (NFI)")
      c2 <- c(c2, sprintf(num_format, x["nfi"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["nfi.scaled"]), ""
      ))
    }

    if ("pnfi" %in% names_x) {
      c1 <- c(c1, "Parsimony Normed Fit Index (PNFI)")
      c2 <- c(c2, sprintf(num_format, x["pnfi"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["pnfi.scaled"]), ""
      ))
    }

    if ("rfi" %in% names_x) {
      c1 <- c(c1, "Bollen's Relative Fit Index (RFI)")
      c2 <- c(c2, sprintf(num_format, x["rfi"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["rfi.scaled"]), ""
      ))
    }

    if ("ifi" %in% names_x) {
      c1 <- c(c1, "Bollen's Incremental Fit Index (IFI)")
      c2 <- c(c2, sprintf(num_format, x["ifi"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["ifi.scaled"]), ""
      ))
    }

    if ("rni" %in% names_x) {
      c1 <- c(c1, "Relative Noncentrality Index (RNI)")
      c2 <- c(c2, sprintf(num_format, x["rni"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["rni.robust"]), ""
      ))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # likelihood
  if ("logl" %in% names_x) {
    cat("\nLoglikelihood and Information Criteria:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    c1 <- c(c1, "Loglikelihood user model (H0)")
    c2 <- c(c2, sprintf(num_format, x["logl"]))
    c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["logl"]), ""))
    if (!is.na(x["scaling.factor.h0"])) {
      c1 <- c(c1, "Scaling correction factor")
      if (scaled_flag) {
        c2 <- c(c2, sprintf("  %10s", ""))
        c3 <- c(c3, sprintf(num_format, x["scaling.factor.h0"]))
      } else {
        # no scaled test statistic (eg random slopes): single column
        c2 <- c(c2, sprintf(num_format, x["scaling.factor.h0"]))
        c3 <- c(c3, "")
      }
      c1 <- c(c1, "    for the MLR correction")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }

    if ("unrestricted.logl" %in% names_x) {
      c1 <- c(c1, "Loglikelihood unrestricted model (H1)")
      c2 <- c(c2, sprintf(num_format, x["unrestricted.logl"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["unrestricted.logl"]), ""
      ))
      if (!is.na(x["scaling.factor.h1"])) {
        c1 <- c(c1, "Scaling correction factor")
        c2 <- c(c2, sprintf("  %10s", ""))
        c3 <- c(c3, sprintf(num_format, x["scaling.factor.h1"]))
        c1 <- c(c1, "    for the MLR correction")
        c2 <- c(c2, "")
        c3 <- c(c3, "")
      }
    }

    c1 <- c(c1, "")
    c2 <- c(c2, "")
    c3 <- c(c3, "")
    # c1 <- c(c1, "Number of free parameters")
    # c2 <- c(c2, sprintf("  %10i", x["npar"]))
    # c3 <- c(c3, ifelse(scaled, sprintf("  %10i", x["npar"]), ""))

    c1 <- c(c1, "Akaike (AIC)")
    c2 <- c(c2, sprintf(num_format, x["aic"]))
    c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["aic"]), ""))
    c1 <- c(c1, "Bayesian (BIC)")
    c2 <- c(c2, sprintf(num_format, x["bic"]))
    c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["bic"]), ""))
    if (!is.na(x["bic2"])) {
      c1 <- c(c1, "Sample-size adjusted Bayesian (SABIC)")
      c2 <- c(c2, sprintf(num_format, x["bic2"]))
      c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["bic2"]), ""))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 39L)
    c2 <- format(c2, width = 12L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # RMSEA
  if ("rmsea" %in% names_x) {
    cat("\nRoot Mean Square Error of Approximation:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    c1 <- c(c1, "RMSEA")
    c2 <- c(c2, sprintf(num_format, x["rmsea"]))
    c3 <- c(c3, ifelse(scaled_flag,
                       sprintf(num_format, x["rmsea.scaled"]), ""
    ))

    ci_level <- NULL
    if ("rmsea.ci.level" %in% names_x) {
      ci_level <- x["rmsea.ci.level"]
    }
    if ("rmsea.ci.lower" %in% names_x) {
      if (is.null(ci_level)) {
        c1 <- c(c1, "Confidence interval - lower")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci_level * 100)),
          " Percent confidence interval - lower"
        ))
      }
      c2 <- c(c2, sprintf(num_format, x["rmsea.ci.lower"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["rmsea.ci.lower.scaled"]), ""
      ))
      if (is.null(ci_level)) {
        c1 <- c(c1, "Confidence interval - upper")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci_level * 100)),
          " Percent confidence interval - upper"
        ))
      }
      c2 <- c(c2, sprintf(num_format, x["rmsea.ci.upper"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["rmsea.ci.upper.scaled"]), ""
      ))
    }

    rmsea_close_h0 <- NULL
    if ("rmsea.close.h0" %in% names_x) {
      rmsea_close_h0 <- x["rmsea.close.h0"]
    }
    rmsea_notclose_h0 <- NULL
    if ("rmsea.notclose.h0" %in% names_x) {
      rmsea_notclose_h0 <- x["rmsea.notclose.h0"]
    }
    if ("rmsea.pvalue" %in% names_x) {
      if (is.null(rmsea_close_h0)) {
        c1 <- c(c1, "P-value H_0: RMSEA <= 0.05")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: RMSEA <= ",
          sprintf("%4.3f", rmsea_close_h0)
        ))
      }
      c2 <- c(c2, sprintf(num_format, x["rmsea.pvalue"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["rmsea.pvalue.scaled"]), ""
      ))
    }
    if ("rmsea.notclose.pvalue" %in% names_x) {
      if (is.null(rmsea_notclose_h0)) {
        c1 <- c(c1, "P-value H_0: RMSEA >= 0.080")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: RMSEA >= ",
          sprintf("%4.3f", rmsea_notclose_h0)
        ))
      }
      c2 <- c(c2, sprintf(num_format, x["rmsea.notclose.pvalue"]))
      c3 <- c(c3, ifelse(scaled_flag,
                    sprintf(num_format, x["rmsea.notclose.pvalue.scaled"]), ""
      ))
    }

    # robust
    if ("rmsea.robust" %in% names_x) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
      c1 <- c(c1, "Robust RMSEA")
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["rmsea.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["rmsea.robust"]))
        c3 <- c(c3, "")
      }
    }
    if ("rmsea.ci.lower.robust" %in% names_x) {
      if (is.null(ci_level)) {
        c1 <- c(c1, "Confidence interval - lower")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci_level * 100)),
          " Percent confidence interval - lower"
        ))
      }
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["rmsea.ci.lower.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["rmsea.ci.lower.robust"]))
        c3 <- c(c3, "")
      }

      if (is.null(ci_level)) {
        c1 <- c(c1, "Confidence interval - upper")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci_level * 100)),
          " Percent confidence interval - upper"
        ))
      }
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["rmsea.ci.upper.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["rmsea.ci.upper.robust"]))
        c3 <- c(c3, "")
      }
    }
    if ("rmsea.pvalue.robust" %in% names_x) {
      if (is.null(rmsea_close_h0)) {
        c1 <- c(c1, "P-value H_0: Robust RMSEA <= 0.05")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: Robust RMSEA <= ",
          sprintf("%4.3f", rmsea_close_h0)
        ))
      }
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["rmsea.pvalue.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["rmsea.pvalue.robust"]))
        c3 <- c(c3, "")
      }
    }
    if ("rmsea.notclose.pvalue.robust" %in% names_x) {
      if (is.null(rmsea_notclose_h0)) {
        c1 <- c(c1, "P-value H_0: Robust RMSEA >= 0.080")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: Robust RMSEA >= ",
          sprintf("%4.3f", rmsea_notclose_h0)
        ))
      }
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(
          c3,
          sprintf(num_format, x["rmsea.notclose.pvalue.robust"])
        )
      } else {
        c2 <- c(
          c2,
          sprintf(num_format, x["rmsea.notclose.pvalue.robust"])
        )
        c3 <- c(c3, "")
      }
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # SRMR
  #TODO: add CRMR
  if (any(c("rmr", "srmr") %in% names_x) && !"srmr_within" %in% names_x) {
    cat("\nStandardized Root Mean Square Residual:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("rmr" %in% names_x) {
      c1 <- c(c1, "RMR")
      c2 <- c(c2, sprintf(num_format, x["rmr"]))
      c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["rmr"]), ""))
    }
    if ("rmr_nomean" %in% names_x) {
      c1 <- c(c1, "RMR (No Mean)")
      c2 <- c(c2, sprintf(num_format, x["rmr_nomean"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["rmr_nomean"]), ""
      ))
    }
    if ("srmr" %in% names_x) {
      c1 <- c(c1, "SRMR")
      c2 <- c(c2, sprintf(num_format, x["srmr"]))
      c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["srmr"]), ""))
    }
    if ("srmr_nomean" %in% names_x) {
      c1 <- c(c1, "SRMR (No Mean)")
      c2 <- c(c2, sprintf(num_format, x["srmr_nomean"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["srmr_nomean"]), ""
      ))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # SRMR -- multilevel
  #TODO: add CRMR?
  if (any(c("srmr_within", "srmr_between") %in% names_x)) {
    cat("\nStandardized Root Mean Square Residual (corr metric):\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("srmr_within" %in% names_x) {
      c1 <- c(c1, "SRMR (within covariance matrix)")
      c2 <- c(c2, sprintf(num_format, x["srmr_within"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["srmr_within"]), ""
      ))
    }
    if ("srmr_between" %in% names_x) {
      c1 <- c(c1, "SRMR (between covariance matrix)")
      c2 <- c(c2, sprintf(num_format, x["srmr_between"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["srmr_between"]), ""
      ))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # WRMR
  if ("wrmr" %in% names_x) {
    cat("\nWeighted Root Mean Square Residual:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("wrmr" %in% names_x) {
      c1 <- c(c1, "WRMR")
      c2 <- c(c2, sprintf(num_format, x["wrmr"]))
      c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["wrmr"]), ""))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # GFI (the 'new' GFI of Maydeu-Olivares et al., 2024)
  if ("gfi" %in% names_x) {
    cat("\nGoodness of Fit Index:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    c1 <- c(c1, "Goodness of Fit Index (GFI)")
    c2 <- c(c2, sprintf(num_format, x["gfi"]))
    c3 <- c(c3, "")

    gfi_ci_level <- NULL
    if ("gfi.ci.level" %in% names_x) {
      gfi_ci_level <- x["gfi.ci.level"]
    }
    if ("gfi.ci.lower" %in% names_x) {
      if (is.null(gfi_ci_level)) {
        c1 <- c(c1, "Confidence interval - lower")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(gfi_ci_level * 100)),
          " Percent confidence interval - lower"
        ))
      }
      c2 <- c(c2, sprintf(num_format, x["gfi.ci.lower"]))
      c3 <- c(c3, "")
      if (is.null(gfi_ci_level)) {
        c1 <- c(c1, "Confidence interval - upper")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(gfi_ci_level * 100)),
          " Percent confidence interval - upper"
        ))
      }
      c2 <- c(c2, sprintf(num_format, x["gfi.ci.upper"]))
      c3 <- c(c3, "")
    }

    # robust
    if ("gfi.robust" %in% names_x) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
      c1 <- c(c1, "Robust GFI")
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["gfi.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["gfi.robust"]))
        c3 <- c(c3, "")
      }
    }
    if ("gfi.ci.lower.robust" %in% names_x) {
      if (is.null(gfi_ci_level)) {
        c1 <- c(c1, "Confidence interval - lower")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(gfi_ci_level * 100)),
          " Percent confidence interval - lower"
        ))
      }
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["gfi.ci.lower.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["gfi.ci.lower.robust"]))
        c3 <- c(c3, "")
      }
      if (is.null(gfi_ci_level)) {
        c1 <- c(c1, "Confidence interval - upper")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(gfi_ci_level * 100)),
          " Percent confidence interval - upper"
        ))
      }
      if (scaled_flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num_format, x["gfi.ci.upper.robust"]))
      } else {
        c2 <- c(c2, sprintf(num_format, x["gfi.ci.upper.robust"]))
        c3 <- c(c3, "")
      }
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # Other
  if (any(c("cn_05", "cn_01", "gfi_lisrel", "agfi_lisrel", "pgfi", "mfi")
          %in% names_x)) {
    cat("\nOther Fit Indices:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("cn_05" %in% names_x) {
      c1 <- c(c1, "Hoelter Critical N (CN) alpha = 0.05")
      c2 <- c(c2, sprintf(num_format, x["cn_05"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["cn_05"]), ""
      ))
    }
    if ("cn_01" %in% names_x) {
      c1 <- c(c1, "Hoelter Critical N (CN) alpha = 0.01")
      c2 <- c(c2, sprintf(num_format, x["cn_01"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["cn_01"]), ""
      ))
    }
    if (any(c("cn_05", "cn_01") %in% names_x)) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }
    if ("gfi_lisrel" %in% names_x) {
      c1 <- c(c1, "Goodness of Fit Index LISREL (GFI)")
      c2 <- c(c2, sprintf(num_format, x["gfi_lisrel"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["gfi_lisrel"]), ""
      ))
    }
    if ("agfi_lisrel" %in% names_x) {
      c1 <- c(c1, "Adjusted Goodness of Fit Index LISREL (AGFI)")
      c2 <- c(c2, sprintf(num_format, x["agfi_lisrel"]))
      c3 <- c(c3, ifelse(scaled_flag,
                         sprintf(num_format, x["agfi_lisrel"]), ""
      ))
    }
    if ("pgfi" %in% names_x) {
      c1 <- c(c1, "Parsimony Goodness of Fit Index (PGFI)")
      c2 <- c(c2, sprintf(num_format, x["pgfi"]))
      c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["pgfi"]), ""))
    }
    if (any(c("gfi_lisrel", "agfi_lisrel", "pgfi") %in% names_x)) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }
    if ("mfi" %in% names_x) {
      c1 <- c(c1, "McDonald Fit Index (MFI)")
      c2 <- c(c2, sprintf(num_format, x["mfi"]))
      c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["mfi"]), ""))
    }
    if ("mfi" %in% names_x) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }
    if ("ecvi" %in% names_x) {
      c1 <- c(c1, "Expected Cross-Validation Index (ECVI)")
      c2 <- c(c2, sprintf(num_format, x["ecvi"]))
      c3 <- c(c3, ifelse(scaled_flag, sprintf(num_format, x["ecvi"]), ""))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled_flag) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  invisible(x)
}
