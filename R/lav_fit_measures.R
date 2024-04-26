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

setMethod(
  "fitMeasures", signature(object = "lavaan"),
  function(object, fit.measures = "all", baseline.model = NULL, h1.model = NULL,
           fm.args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.check.pd = TRUE
           ),
           output = "vector", ...) {
    # note: the ... is not used by lavaan
    lav_fit_measures(
      object = object, fit.measures = fit.measures,
      baseline.model = baseline.model, h1.model = h1.model, fm.args = fm.args,
      output = output
    )
  }
)

setMethod(
  "fitmeasures", signature(object = "lavaan"),
  function(object, fit.measures = "all", baseline.model = NULL, h1.model = NULL,
           fm.args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.close.h0 = 0.05,
             rmsea.notclose.h0 = 0.08,
             robust = TRUE,
             cat.check.pd = TRUE
           ),
           output = "vector", ...) {
    # note: the ... is not used by lavaan
    lav_fit_measures(
      object = object, fit.measures = fit.measures,
      baseline.model = baseline.model, h1.model = h1.model, fm.args = fm.args,
      output = output
    )
  }
)

# S3 method for efaList
fitMeasures.efaList <- fitmeasures.efaList <- function(
    object,
    fit.measures = "all",
    baseline.model = NULL, h1.model = NULL,
    fm.args = list(
      standard.test = "default",
      scaled.test = "default",
      rmsea.ci.level = 0.90,
      rmsea.close.h0 = 0.05,
      rmsea.notclose.h0 = 0.08,
      robust = TRUE,
      cat.check.pd = TRUE
    ),
    output = "list", ...) {
  # kill object$loadings if present
  object[["loadings"]] <- NULL

  # get fit measures for each model
  res <- simplify2array(lapply(
    object,
    function(x) {
      lav_fit_measures(
        object = x,
        fit.measures = fit.measures, h1.model = h1.model,
        baseline.model = baseline.model, fm.args = fm.args,
        output = "vector"
      )
    }
  ))

  # check if res is a matrix
  if (!is.matrix(res)) {
    if (is.numeric(res)) {
      # fit.measures is just 1 element, or only one was correct
      NAME <- names(res)[1]
      res <- matrix(res, nrow = 1L)
      rownames(res) <- NAME
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


lav_fit_measures <- function(object, fit.measures = "all",
                             baseline.model = NULL, h1.model = NULL,
                             fm.args = list(
                               standard.test = "default",
                               scaled.test = "default",
                               rmsea.ci.level = 0.90,
                               rmsea.close.h0 = 0.05,
                               rmsea.notclose.h0 = 0.08,
                               robust = TRUE,
                               cat.check.pd = TRUE
                             ),
                             output = "vector") {
  # default fm.args
  default.fm.args <- list(
    standard.test = "default",
    scaled.test = "default",
    rmsea.ci.level = 0.90,
    rmsea.close.h0 = 0.05,
    rmsea.notclose.h0 = 0.08,
    robust = TRUE,
    cat.check.pd = TRUE
  )
  if (!missing(fm.args)) {
    fm.args <- modifyList(default.fm.args, fm.args)
  } else {
    fm.args <- default.fm.args
  }

  # standard test
  if (fm.args$standard.test == "default") {
    fm.args$standard.test <- object@Options$scaled.test
    # usually "standard", but could have been changed
    # the 'scaled' version will be based on the scaled.test!
    if (is.null(fm.args$standard.test)) { # <older objects
      fm.args$standard.test <- "standard"
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

  # do we have a test statistic?
  TEST <- lavInspect(object, "test")
  test.names <- unname(sapply(TEST, "[[", "test"))
  if (test.names[1] == "none") {
    lav_msg_stop(gettext("fit measures not available if test = \"none\"."))
    #FIXME: allow RMRs, log.likelihoods, info criteria, npar, ntotal
  }

  standard.test <- fm.args$standard.test
  scaled.test <- fm.args$scaled.test

  # check standard.test
  standard.test <- lav_test_rename(standard.test, check = TRUE)[1] # only 1

  # check scaled.test
  if (!scaled.test %in% c("none", "default", "standard")) {
    scaled.test <- lav_test_rename(scaled.test, check = TRUE)[1] # only 1
  }

  # which test statistic do we need?
  rerun.lavtest.flag <- FALSE
  if (!standard.test %in% test.names) {
    rerun.lavtest.flag <- TRUE
  }
  if (!scaled.test %in% c("none", "default", "standard") &&
    !scaled.test %in% test.names) {
    rerun.lavtest.flag <- TRUE
  }

  # do we have a scaled test statistic? if so, which one?
  scaled.flag <- FALSE
  if (scaled.test != "none" &&
    any(test.names %in% c(
      "satorra.bentler",
      "yuan.bentler", "yuan.bentler.mplus",
      "mean.var.adjusted", "scaled.shifted"
    ))) {
    scaled.flag <- TRUE
    if (scaled.test %in% c("standard", "default")) {
      tmp.idx <- which(test.names %in% c(
        "satorra.bentler",
        "yuan.bentler", "yuan.bentler.mplus",
        "mean.var.adjusted", "scaled.shifted"
      ))
      scaled.test <- test.names[tmp.idx[1]]
    }
  }

  # rerun lavTest?
  if (rerun.lavtest.flag) {
    this.test <- standard.test
    if (scaled.flag) {
      this.test <- unique(this.test, scaled.test)
    }
    TEST <- lavTest(object,
      test = this.test, scaled.test = standard.test,
      drop.list.single = FALSE
    )
    # replace in object, if we pass it to lav_fit_* functions
    object@test <- TEST
    test.names <- unname(sapply(TEST, "[[", "test"))
  }


  # TDJ: Check for user-supplied h1 model
  #      Similar to BASELINE model, use the following priority:
  #        1. user-provided h1 model
  #        2. h1 model in @external slot
  #        3. default h1 model (already in @h1 slot, no update necessary)

  user_h1_exists <- FALSE
    # 1. user-provided h1 model
  if (!is.null(h1.model)) {
    stopifnot(inherits(h1.model, "lavaan"))
    user_h1_exists <- TRUE

    # 2. h1 model in @external slot
  } else if (!is.null(object@external$h1.model)) {
    stopifnot(inherits(object@external$h1.model, "lavaan"))
    h1.model <- object@external$h1.model
    user_h1_exists <- TRUE
  }

  ## Update statistics in @test slot?
  if (user_h1_exists) {
    ## update @test slot
    FIT <- lav_update_test_custom_h1(lav_obj_h0 = object, lav_obj_h1 = h1.model)

    ## re-assign TEST object that is used below
    object@test <- TEST <- FIT@test
    test.names <- unname(sapply(TEST, "[[", "test"))
  }

  # get index of standard.test in TEST
  test.idx <- which(test.names == standard.test)[1]

  # get index of scaled test (if any) in TEST
  if (scaled.flag) {
    scaled.idx <- which(test.names == scaled.test)[1]
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
      lav_msg_view(c("vector", "list", "matrix", "text"), "none", FALSE)
    ))
  }

  # options
  categorical.flag <- object@Model@categorical
  fiml.flag <- (fm.args$robust &&
    object@Options$missing %in% c("ml", "ml.x"))
  estimator <- object@Options$estimator

  # basic ingredients
  G <- object@Data@ngroups
  X2 <- TEST[[test.idx]]$stat
  df <- TEST[[test.idx]]$df
  if (scaled.flag) {
    X2.scaled <- TEST[[scaled.idx]]$stat
    df.scaled <- TEST[[scaled.idx]]$df
  }
  npar <- lav_utils_get_npar(lavobject = object)
  N <- lav_utils_get_ntotal(lavobject = object) # N vs N-1


  # define 'sets' of fit measures:
  fit.always <- c("npar")

  # basic chi-square test
  fit.chisq <- c("fmin", "chisq", "df", "pvalue")
  if (scaled.flag) {
    fit.chisq <- c(
      fit.chisq, "chisq.scaled", "df.scaled", "pvalue.scaled",
      "chisq.scaling.factor"
    )
  }

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
  if (fm.args$robust && (scaled.flag || categorical.flag || fiml.flag)) {
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
  if (fm.args$robust && (scaled.flag || categorical.flag || fiml.flag)) {
    fit.cfi.other <- c(fit.cfi.other, "nnfi.robust", "rni.robust")
  }
  fit.cfi <- c(fit.baseline, fit.cfi.tli, fit.cfi.other)

  # likelihood based measures
  if (estimator == "MML") {
    fit.logl <- c("logl", "aic", "bic", "ntotal", "bic2")
  } else {
    fit.logl <- c(
      "logl", "unrestricted.logl", "aic", "bic",
      "ntotal", "bic2"
    )
  }
  if (scaled.flag &&
    scaled.test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
    fit.logl <- c(fit.logl, "scaling.factor.h1", "scaling.factor.h0")
  }

  # rmsea
  fit.rmsea <- c(
    "rmsea",
    "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.ci.level",
    "rmsea.pvalue", "rmsea.close.h0",
    "rmsea.notclose.pvalue", "rmsea.notclose.h0"
  )
  if (scaled.flag) {
    fit.rmsea <- c(
      fit.rmsea, "rmsea.scaled", "rmsea.ci.lower.scaled",
      "rmsea.ci.upper.scaled", "rmsea.pvalue.scaled",
      "rmsea.notclose.pvalue.scaled"
    )
  }
  if (fm.args$robust && (scaled.flag || categorical.flag || fiml.flag)) {
    fit.rmsea <- c(
      fit.rmsea, "rmsea.robust", "rmsea.ci.lower.robust",
      "rmsea.ci.upper.robust", "rmsea.pvalue.robust",
      "rmsea.notclose.pvalue.robust"
    )
  }

  # srmr
  if (categorical.flag) {
    fit.srmr <- c("srmr")
    fit.srmr2 <- c(
      "rmr", "rmr_nomean",
      "srmr", # per default equal to srmr_bentler_nomean
      "srmr_bentler", "srmr_bentler_nomean",
      "crmr", "crmr_nomean",
      "srmr_mplus", "srmr_mplus_nomean"
    )
  } else {
    if (object@Data@nlevels > 1L) {
      fit.srmr <- c("srmr", "srmr_within", "srmr_between")
      fit.srmr2 <- c("srmr", "srmr_within", "srmr_between")
    } else {
      fit.srmr <- c("srmr")
      fit.srmr2 <- c(
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
    fit.other <- ""
  } else if (estimator == "PML") {
    fit.other <- c("cn_05", "cn_01", "mfi")
    if (!categorical.flag) { # really needed?
      fit.other <- c(fit.other, "ecvi")
    }
  } else {
    fit.other <- c("cn_05", "cn_01", "gfi", "agfi", "pgfi", "mfi")
    if (!categorical.flag) { # really needed?
      fit.other <- c(fit.other, "ecvi")
    } else {
      fit.other <- c(fit.other, "wrmr")
    }
  }

  # lower case
  fit.measures <- tolower(fit.measures)

  # select 'default' fit measures
  if (length(fit.measures) == 1L) {
    if (fit.measures == "default") {
      if (estimator == "ML" || estimator == "PML") {
        fit.measures <- c(
          fit.always, fit.chisq, fit.baseline,
          fit.cfi.tli, fit.logl,
          fit.rmsea, fit.srmr
        )
      } else if (estimator == "MML") {
        fit.measures <- c(fit.always, fit.logl)
      } else {
        fit.measures <- c(
          fit.always,
          fit.chisq, fit.baseline, fit.cfi.tli,
          fit.rmsea, fit.srmr
        )
      }
    } else if (fit.measures == "all") {
      if (estimator == "ML") {
        fit.measures <- c(
          fit.always, fit.chisq,
          fit.baseline, fit.cfi.tli, fit.cfi.other,
          fit.logl, fit.rmsea, fit.srmr2, fit.other
        )
      } else {
        fit.measures <- c(
          fit.always, fit.chisq,
          fit.baseline, fit.cfi.tli, fit.cfi.other,
          fit.rmsea, fit.srmr2, fit.other
        )
      }
    }
  }

  # catch empty list
  if (length(fit.measures) == 0L) {
    return(list())
  }

  # main container
  indices <- list()
  indices["npar"] <- npar
  indices["ntotal"] <- object@SampleStats@ntotal
  indices["fmin"] <- object@optim$fx # note = 0.5 * fmin if ML

  # CHI-SQUARE TEST
  if (any(fit.chisq %in% fit.measures)) {
    indices["chisq"] <- X2
    indices["df"] <- df
    indices["pvalue"] <- TEST[[test.idx]]$pvalue
    if (scaled.flag) {
      indices["chisq.scaled"] <- X2.scaled
      indices["df.scaled"] <- df.scaled
      indices["chisq.scaling.factor"] <- TEST[[scaled.idx]]$scaling.factor
      indices["pvalue.scaled"] <- TEST[[scaled.idx]]$pvalue
    }
  }

  # BASELINE FAMILY
  if (any(fit.cfi %in% fit.measures)) {
    # rerun baseline?
    if (rerun.lavtest.flag) {
      object@Options$test <- this.test
      fit.indep <- try(lav_object_independence(object), silent = TRUE)
      # override object
      object@baseline$test <- fit.indep@test
    }

    indices <- c(
      indices,
      lav_fit_cfi_lavobject(
        lavobject = object,
        fit.measures = fit.measures,
        baseline.model = baseline.model,
        h1.model = h1.model,
        standard.test = standard.test,
        scaled.test = scaled.test,
        robust = fm.args$robust,
        cat.check.pd = fm.args$cat.check.pd
      )
    )
  }

  # INFORMATION CRITERIA
  if (any(fit.logl %in% fit.measures)) {
    indices <- c(
      indices,
      lav_fit_aic_lavobject(
        lavobject = object,
        fit.measures = fit.measures,
        standard.test = standard.test,
        scaled.test = scaled.test,
        estimator = estimator
      )
    )
  }

  # RMSEA and friends
  if (any(fit.rmsea %in% fit.measures)) {
    # check rmsea options
    rmsea.ci.level <- 0.90
    rmsea.close.h0 <- 0.05
    rmsea.notclose.h0 <- 0.08
    if (!is.null(fm.args$rmsea.ci.level) &&
      is.finite(fm.args$rmsea.ci.level)) {
      rmsea.ci.level <- fm.args$rmsea.ci.level
      if (rmsea.ci.level < 0 || rmsea.ci.level > 1.0) {
        lav_msg_warn(gettextf(
          "invalid rmsea.ci.level value [%s] set to default 0.90.",
          rmsea.ci.level))
        rmsea.ci.level <- 0.90
      }
    }
    if (!is.null(fm.args$rmsea.close.h0) &&
      is.finite(fm.args$rmsea.close.h0)) {
      rmsea.close.h0 <- fm.args$rmsea.close.h0
      if (rmsea.close.h0 < 0) {
        rmsea.close.h0 <- 0
      }
    }
    if (!is.null(fm.args$rmsea.notclose.h0) &&
      is.finite(fm.args$rmsea.notclose.h0)) {
      rmsea.notclose.h0 <- fm.args$rmsea.notclose.h0
      if (rmsea.notclose.h0 < 0) {
        rmsea.notclose.h0 <- 0
      }
    }
    indices <- c(
      indices,
      lav_fit_rmsea_lavobject(
        lavobject = object,
        fit.measures = fit.measures,
        standard.test = standard.test,
        scaled.test = scaled.test,
        ci.level = rmsea.ci.level,
        close.h0 = rmsea.close.h0,
        notclose.h0 = rmsea.notclose.h0,
        robust = fm.args$robust,
        cat.check.pd = fm.args$cat.check.pd
      )
    )
  }

  # SRMR and friends
  if (any(fit.srmr2 %in% fit.measures)) {
    indices <- c(
      indices,
      lav_fit_srmr_lavobject(
        lavobject = object,
        fit.measures = fit.measures
      )
    )
  }

  # GFI and friends
  fit.gfi <- c("gfi", "agfi", "pgfi")
  if (any(fit.gfi %in% fit.measures)) {
    indices <- c(
      indices,
      lav_fit_gfi_lavobject(
        lavobject = object,
        fit.measures = fit.measures
      )
    )
  }

  # various: Hoelter Critical N (CN)
  if (any(c("cn_05", "cn_01") %in% fit.measures)) {
    indices["cn_05"] <- lav_fit_cn(X2 = X2, df = df, N = N, alpha = 0.05)
    indices["cn_01"] <- lav_fit_cn(X2 = X2, df = df, N = N, alpha = 0.01)
  }

  # various: WRMR
  if ("wrmr" %in% fit.measures) {
    nel <- length(object@SampleStats@WLS.obs[[1]])
    indices["wrmr"] <- lav_fit_wrmr(X2 = X2, nel = nel)
  }

  # various: MFI
  if ("mfi" %in% fit.measures) {
    indices["mfi"] <- lav_fit_mfi(X2 = X2, df = df, N = N)
  }

  # various: ECVI
  if ("ecvi" %in% fit.measures) {
    indices["ecvi"] <- lav_fit_ecvi(X2 = X2, npar = npar, N = N)
  }


  # keep only what we need
  out <- indices[fit.measures]

  if (all(is.na(names(out)))) { # perhaps, fit.measures = ""
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

  out
}

# print a nice summary of the fit measures
print.lavaan.fitMeasures <- function(x, ..., nd = 3L, add.h0 = TRUE) {
  names.x <- names(x)

  # scaled?
  scaled.flag <- "chisq.scaled" %in% names.x

  # num format
  num.format <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

  ## TDJ: optionally add h0 model's fit statistic, for lavaan.mi
  if (add.h0 && "chisq" %in% names.x) {
    cat("\nModel Test User Model:\n\n")

    # container three columns
    c1 <- c2 <- c3 <- character(0L)

    # TDJ: Add header used in summary() by lavaan.mi
    if (scaled.flag) {
      c1 <- c("", c1)
      c2 <- c("Standard", c2)
      c3 <- c("Scaled", c3)
    }

    c1 <- c(c1, "Test statistic")
    c2 <- c(c2, sprintf(num.format, x["chisq"]))
    c3 <- c(c3, ifelse(scaled.flag,
      sprintf(num.format, x["chisq.scaled"]), ""
    ))

    c1 <- c(c1, "Degrees of freedom")
    c2 <- c(c2, x["df"])
    c3 <- c(c3, ifelse(scaled.flag,
      ifelse(x["df.scaled"] %% 1 == 0, x["df.scaled"],
        sprintf(num.format, x["df.scaled"])
      ), ""
    ))

    c1 <- c(c1, "P-value")
    c2 <- c(c2, sprintf(num.format, x["pvalue"]))
    c3 <- c(c3, ifelse(scaled.flag,
      sprintf(num.format, x["pvalue.scaled"]), ""
    ))

    if (scaled.flag && "chisq.scaling.factor" %in% names.x) {
      c1 <- c(c1, "Scaling correction factor")
      c2 <- c(c2, "")
      c3 <- c(c3, sprintf(num.format, x["chisq.scaling.factor"]))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 35L)
    c2 <- format(c2, width = 16L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # independence model
  if ("baseline.chisq" %in% names.x) {
    cat("\nModel Test Baseline Model:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    c1 <- c(c1, "Test statistic")
    c2 <- c(c2, sprintf(num.format, x["baseline.chisq"]))
    c3 <- c(c3, ifelse(scaled.flag,
      sprintf(num.format, x["baseline.chisq.scaled"]), ""
    ))

    c1 <- c(c1, "Degrees of freedom")
    c2 <- c(c2, x["baseline.df"])
    c3 <- c(c3, ifelse(scaled.flag, ifelse(x["baseline.df.scaled"] %% 1 == 0,
      x["baseline.df.scaled"],
      sprintf(num.format, x["baseline.df.scaled"])
    ),
    ""
    ))

    c1 <- c(c1, "P-value")
    c2 <- c(c2, sprintf(num.format, x["baseline.pvalue"]))
    c3 <- c(c3, ifelse(scaled.flag,
      sprintf(num.format, x["baseline.pvalue.scaled"]), ""
    ))

    if (scaled.flag && "baseline.chisq.scaling.factor" %in% names.x) {
      c1 <- c(c1, "Scaling correction factor")
      c2 <- c(c2, "")
      c3 <- c(c3, sprintf(num.format, x["baseline.chisq.scaling.factor"]))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 35L)
    c2 <- format(c2, width = 16L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # cfi/tli
  if (any(c("cfi", "tli", "nnfi", "rfi", "nfi", "ifi", "rni", "pnfi") %in% names.x)) {
    cat("\nUser Model versus Baseline Model:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("cfi" %in% names.x) {
      c1 <- c(c1, "Comparative Fit Index (CFI)")
      c2 <- c(c2, sprintf(num.format, x["cfi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["cfi.scaled"]), ""
      ))
    }

    if ("tli" %in% names.x) {
      c1 <- c(c1, "Tucker-Lewis Index (TLI)")
      c2 <- c(c2, sprintf(num.format, x["tli"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["tli.scaled"]), ""
      ))
    }

    if ("cfi.robust" %in% names.x) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
      c1 <- c(c1, "Robust Comparative Fit Index (CFI)")
      if (scaled.flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num.format, x["cfi.robust"]))
      } else {
        c2 <- c(c2, sprintf(num.format, x["cfi.robust"]))
        c3 <- c(c3, "")
      }
    }

    if ("tli.robust" %in% names.x) {
      c1 <- c(c1, "Robust Tucker-Lewis Index (TLI)")
      if (scaled.flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num.format, x["tli.robust"]))
      } else {
        c2 <- c(c2, sprintf(num.format, x["tli.robust"]))
        c3 <- c(c3, "")
      }
    }

    if ("nnfi" %in% names.x) {
      c1 <- c(c1, "Bentler-Bonett Non-normed Fit Index (NNFI)")
      c2 <- c(c2, sprintf(num.format, x["nnfi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["nnfi.robust"]), ""
      ))
    }

    if ("nfi" %in% names.x) {
      c1 <- c(c1, "Bentler-Bonett Normed Fit Index (NFI)")
      c2 <- c(c2, sprintf(num.format, x["nfi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["nfi.scaled"]), ""
      ))
    }

    if ("pnfi" %in% names.x) {
      c1 <- c(c1, "Parsimony Normed Fit Index (PNFI)")
      c2 <- c(c2, sprintf(num.format, x["pnfi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["pnfi.scaled"]), ""
      ))
    }

    if ("rfi" %in% names.x) {
      c1 <- c(c1, "Bollen's Relative Fit Index (RFI)")
      c2 <- c(c2, sprintf(num.format, x["rfi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["rfi.scaled"]), ""
      ))
    }

    if ("ifi" %in% names.x) {
      c1 <- c(c1, "Bollen's Incremental Fit Index (IFI)")
      c2 <- c(c2, sprintf(num.format, x["ifi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["ifi.scaled"]), ""
      ))
    }

    if ("rni" %in% names.x) {
      c1 <- c(c1, "Relative Noncentrality Index (RNI)")
      c2 <- c(c2, sprintf(num.format, x["rni"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["rni.robust"]), ""
      ))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # likelihood
  if ("logl" %in% names.x) {
    cat("\nLoglikelihood and Information Criteria:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    c1 <- c(c1, "Loglikelihood user model (H0)")
    c2 <- c(c2, sprintf(num.format, x["logl"]))
    c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["logl"]), ""))
    if (!is.na(x["scaling.factor.h0"])) {
      c1 <- c(c1, "Scaling correction factor")
      c2 <- c(c2, sprintf("  %10s", ""))
      c3 <- c(c3, sprintf(num.format, x["scaling.factor.h0"]))
      c1 <- c(c1, "    for the MLR correction")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }

    if ("unrestricted.logl" %in% names.x) {
      c1 <- c(c1, "Loglikelihood unrestricted model (H1)")
      c2 <- c(c2, sprintf(num.format, x["unrestricted.logl"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["unrestricted.logl"]), ""
      ))
      if (!is.na(x["scaling.factor.h1"])) {
        c1 <- c(c1, "Scaling correction factor")
        c2 <- c(c2, sprintf("  %10s", ""))
        c3 <- c(c3, sprintf(num.format, x["scaling.factor.h1"]))
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
    c2 <- c(c2, sprintf(num.format, x["aic"]))
    c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["aic"]), ""))
    c1 <- c(c1, "Bayesian (BIC)")
    c2 <- c(c2, sprintf(num.format, x["bic"]))
    c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["bic"]), ""))
    if (!is.na(x["bic2"])) {
      c1 <- c(c1, "Sample-size adjusted Bayesian (SABIC)")
      c2 <- c(c2, sprintf(num.format, x["bic2"]))
      c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["bic2"]), ""))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 39L)
    c2 <- format(c2, width = 12L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # RMSEA
  if ("rmsea" %in% names.x) {
    cat("\nRoot Mean Square Error of Approximation:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    c1 <- c(c1, "RMSEA")
    c2 <- c(c2, sprintf(num.format, x["rmsea"]))
    c3 <- c(c3, ifelse(scaled.flag,
      sprintf(num.format, x["rmsea.scaled"]), ""
    ))

    ci.level <- NULL
    if ("rmsea.ci.level" %in% names.x) {
      ci.level <- x["rmsea.ci.level"]
    }
    if ("rmsea.ci.lower" %in% names.x) {
      if (is.null(ci.level)) {
        c1 <- c(c1, "Confidence interval - lower")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci.level * 100)),
          " Percent confidence interval - lower"
        ))
      }
      c2 <- c(c2, sprintf(num.format, x["rmsea.ci.lower"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["rmsea.ci.lower.scaled"]), ""
      ))
      if (is.null(ci.level)) {
        c1 <- c(c1, "Confidence interval - upper")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci.level * 100)),
          " Percent confidence interval - upper"
        ))
      }
      c2 <- c(c2, sprintf(num.format, x["rmsea.ci.upper"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["rmsea.ci.upper.scaled"]), ""
      ))
    }

    rmsea.close.h0 <- NULL
    if ("rmsea.close.h0" %in% names.x) {
      rmsea.close.h0 <- x["rmsea.close.h0"]
    }
    rmsea.notclose.h0 <- NULL
    if ("rmsea.notclose.h0" %in% names.x) {
      rmsea.notclose.h0 <- x["rmsea.notclose.h0"]
    }
    if ("rmsea.pvalue" %in% names.x) {
      if (is.null(rmsea.close.h0)) {
        c1 <- c(c1, "P-value H_0: RMSEA <= 0.05")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: RMSEA <= ",
          sprintf("%4.3f", rmsea.close.h0)
        ))
      }
      c2 <- c(c2, sprintf(num.format, x["rmsea.pvalue"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["rmsea.pvalue.scaled"]), ""
      ))
    }
    if ("rmsea.notclose.pvalue" %in% names.x) {
      if (is.null(rmsea.notclose.h0)) {
        c1 <- c(c1, "P-value H_0: RMSEA >= 0.080")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: RMSEA >= ",
          sprintf("%4.3f", rmsea.notclose.h0)
        ))
      }
      c2 <- c(c2, sprintf(num.format, x["rmsea.notclose.pvalue"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["rmsea.notclose.pvalue.scaled"]), ""
      ))
    }

    # robust
    if ("rmsea.robust" %in% names.x) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
      c1 <- c(c1, "Robust RMSEA")
      if (scaled.flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num.format, x["rmsea.robust"]))
      } else {
        c2 <- c(c2, sprintf(num.format, x["rmsea.robust"]))
        c3 <- c(c3, "")
      }
    }
    if ("rmsea.ci.lower.robust" %in% names.x) {
      if (is.null(ci.level)) {
        c1 <- c(c1, "Confidence interval - lower")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci.level * 100)),
          " Percent confidence interval - lower"
        ))
      }
      if (scaled.flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num.format, x["rmsea.ci.lower.robust"]))
      } else {
        c2 <- c(c2, sprintf(num.format, x["rmsea.ci.lower.robust"]))
        c3 <- c(c3, "")
      }

      if (is.null(ci.level)) {
        c1 <- c(c1, "Confidence interval - upper")
      } else {
        c1 <- c(c1, paste0(
          sprintf("%2d", round(ci.level * 100)),
          " Percent confidence interval - upper"
        ))
      }
      if (scaled.flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num.format, x["rmsea.ci.upper.robust"]))
      } else {
        c2 <- c(c2, sprintf(num.format, x["rmsea.ci.upper.robust"]))
        c3 <- c(c3, "")
      }
    }
    if ("rmsea.pvalue.robust" %in% names.x) {
      if (is.null(rmsea.close.h0)) {
        c1 <- c(c1, "P-value H_0: Robust RMSEA <= 0.05")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: Robust RMSEA <= ",
          sprintf("%4.3f", rmsea.close.h0)
        ))
      }
      if (scaled.flag) {
        c2 <- c(c2, "")
        c3 <- c(c3, sprintf(num.format, x["rmsea.pvalue.robust"]))
      } else {
        c2 <- c(c2, sprintf(num.format, x["rmsea.pvalue.robust"]))
        c3 <- c(c3, "")
      }
    }
    if ("rmsea.notclose.pvalue.robust" %in% names.x) {
      if (is.null(rmsea.notclose.h0)) {
        c1 <- c(c1, "P-value H_0: Robust RMSEA >= 0.080")
      } else {
        c1 <- c(c1, paste0(
          "P-value H_0: Robust RMSEA >= ",
          sprintf("%4.3f", rmsea.notclose.h0)
        ))
      }
      if (scaled.flag) {
        c2 <- c(c2, "")
        c3 <- c(
          c3,
          sprintf(num.format, x["rmsea.notclose.pvalue.robust"])
        )
      } else {
        c2 <- c(
          c2,
          sprintf(num.format, x["rmsea.notclose.pvalue.robust"])
        )
        c3 <- c(c3, "")
      }
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # SRMR
  #TODO: add CRMR
  if (any(c("rmr", "srmr") %in% names.x) && !"srmr_within" %in% names.x) {
    cat("\nStandardized Root Mean Square Residual:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("rmr" %in% names.x) {
      c1 <- c(c1, "RMR")
      c2 <- c(c2, sprintf(num.format, x["rmr"]))
      c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["rmr"]), ""))
    }
    if ("rmr_nomean" %in% names.x) {
      c1 <- c(c1, "RMR (No Mean)")
      c2 <- c(c2, sprintf(num.format, x["rmr_nomean"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["rmr_nomean"]), ""
      ))
    }
    if ("srmr" %in% names.x) {
      c1 <- c(c1, "SRMR")
      c2 <- c(c2, sprintf(num.format, x["srmr"]))
      c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["srmr"]), ""))
    }
    if ("srmr_nomean" %in% names.x) {
      c1 <- c(c1, "SRMR (No Mean)")
      c2 <- c(c2, sprintf(num.format, x["srmr_nomean"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["srmr_nomean"]), ""
      ))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # SRMR -- multilevel
  #TODO: add CRMR?
  if (any(c("srmr_within", "srmr_between") %in% names.x)) {
    cat("\nStandardized Root Mean Square Residual (corr metric):\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("srmr_within" %in% names.x) {
      c1 <- c(c1, "SRMR (within covariance matrix)")
      c2 <- c(c2, sprintf(num.format, x["srmr_within"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["srmr_within"]), ""
      ))
    }
    if ("srmr_between" %in% names.x) {
      c1 <- c(c1, "SRMR (between covariance matrix)")
      c2 <- c(c2, sprintf(num.format, x["srmr_between"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["srmr_between"]), ""
      ))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # WRMR
  if ("wrmr" %in% names.x) {
    cat("\nWeighted Root Mean Square Residual:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("wrmr" %in% names.x) {
      c1 <- c(c1, "WRMR")
      c2 <- c(c2, sprintf(num.format, x["wrmr"]))
      c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["wrmr"]), ""))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  # Other
  if (any(c("cn_05", "cn_01", "gfi", "agfi", "pgfi", "mfi") %in% names.x)) {
    cat("\nOther Fit Indices:\n\n")

    c1 <- c2 <- c3 <- character(0L)

    if ("cn_05" %in% names.x) {
      c1 <- c(c1, "Hoelter Critical N (CN) alpha = 0.05")
      c2 <- c(c2, sprintf(num.format, x["cn_05"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["cn_05"]), ""
      ))
    }
    if ("cn_01" %in% names.x) {
      c1 <- c(c1, "Hoelter Critical N (CN) alpha = 0.01")
      c2 <- c(c2, sprintf(num.format, x["cn_01"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["cn_01"]), ""
      ))
    }
    if (any(c("cn_05", "cn_01") %in% names.x)) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }
    if ("gfi" %in% names.x) {
      c1 <- c(c1, "Goodness of Fit Index (GFI)")
      c2 <- c(c2, sprintf(num.format, x["gfi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["gfi"]), ""
      ))
    }
    if ("agfi" %in% names.x) {
      c1 <- c(c1, "Adjusted Goodness of Fit Index (AGFI)")
      c2 <- c(c2, sprintf(num.format, x["agfi"]))
      c3 <- c(c3, ifelse(scaled.flag,
        sprintf(num.format, x["agfi"]), ""
      ))
    }
    if ("pgfi" %in% names.x) {
      c1 <- c(c1, "Parsimony Goodness of Fit Index (PGFI)")
      c2 <- c(c2, sprintf(num.format, x["pgfi"]))
      c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["pgfi"]), ""))
    }
    if (any(c("gfi", "agfi", "pgfi") %in% names.x)) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }
    if ("mfi" %in% names.x) {
      c1 <- c(c1, "McDonald Fit Index (MFI)")
      c2 <- c(c2, sprintf(num.format, x["mfi"]))
      c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["mfi"]), ""))
    }
    if ("mfi" %in% names.x) {
      c1 <- c(c1, "")
      c2 <- c(c2, "")
      c3 <- c(c3, "")
    }
    if ("ecvi" %in% names.x) {
      c1 <- c(c1, "Expected Cross-Validation Index (ECVI)")
      c2 <- c(c2, sprintf(num.format, x["ecvi"]))
      c3 <- c(c3, ifelse(scaled.flag, sprintf(num.format, x["ecvi"]), ""))
    }

    # format c1/c2/c3
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (scaled.flag) {
      M <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      M <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(M) <- rep("", ncol(M))
    rownames(M) <- rep(" ", nrow(M))

    # print
    write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)
  }

  invisible(x)
}
