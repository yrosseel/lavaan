# compare two nested models, by default using the chi-square
# difference test

# - in 0.5-16, SB.classic = TRUE is the default again (for now)
# - in 0.5-18, SB.classic is replaced by 'method', with the following
#   options:
#     method = "default" (we choose a default method, based on the estimator)
#     method = "standard" (option to explicitly avoid robust adjustment)
#     method = "Satorra.2000"
#     method = "Satorra.Bentler.2001"
#     method = "Satorra.Bentler.2010"
#     method = "mean.var.adjusted.PLRT"
#
# - 0.6-13: RMSEA.D (also known as 'RDR') is added to the table (also if scaled,
#           since 0.6-20)
# - 0.6-13: fix multiple-group UG^2 bug in Satorra.2000 (reported by
#           Gronneberg, Foldnes and Moss)
#
# - 0.6-18:
#     New option method = "standard" (to explicitly avoid robust adjustment)
#     New test= argument to select scaled stat when
#                           method="satorra.bentler.2001/2010"


lavTestLRT <- function(object, ..., method = "default", test = "default",   # nolint start
                       A.method = "delta", scaled.shifted = TRUE, # only when method="Satorra.2000"
                       type = "Chisq", model.names = NULL) {                # nolint end
  type <- tolower(type[1])
  test <- tolower(test[1])
  method <- tolower(gsub("[-_\\.]", "", method[1]))
  if (type %in% c("browne", "browne.residual.adf", "browne.residual.nt")) {
    if (type == "browne") {
      type <- "browne.residual.adf"
    }
    if (!method %in% c("default", "standard")) {
      lav_msg_stop(gettext(
        "method cannot be used if type is browne.residual.adf or
        browne.residual.nt"))
    }
    method <- "default"
  }

  # NOTE: if we add additional arguments, it is not the same generic
  # anova() function anymore, and match.call will be screwed up

  mcall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  modp <- if (length(dots)) {
    sapply(dots, inherits, "lavaan")
  } else {
    logical(0L)
  }
  # check object
  object <- lav_object_check_version(object)
  # check models in dots
  dots[modp] <- lapply(dots[modp], lav_object_check_version)

  # some general properties (taken from the first model)
  estimator <- object@Options$estimator
  ngroups <- object@Data@ngroups
  ntotal <- object@SampleStats@ntotal

  # TDJ: check for user-supplied h1 model
  user_h1_exists <- FALSE
  if (!is.null(object@external$h1.model)) {
    if (inherits(object@external$h1.model, "lavaan")) {
      user_h1_exists <- TRUE
    }
  }

  # shortcut for single argument (just plain LRT)
  if (!any(modp) && !user_h1_exists) {
    if (type == "cf") {
      lav_msg_warn(gettext("`type' argument is ignored for a single model"))
    }
    return(lav_test_lrt_single_model(object, method = method,
       test = test, type = type))
  }

  # list of models
  mods <- c(list(object), dots[modp])
  if (!is.null(model.names)) {
    names(mods) <- model.names
  } else {
    names(mods) <- sapply(
      as.list(mcall)[which(c(FALSE, TRUE, modp))],
      function(x) deparse(x)
    )
  }
  # TDJ: Add user-supplied h1 model, if it exists
  if (user_h1_exists) mods$user_h1 <- object@external$h1.model

  # put them in order (using degrees of freedom)
  ndf <- sapply(mods, function(x) x@test[[1]]$df)
  order_idx <- order(ndf)
  mods <- mods[order_idx]
  ndf <- ndf[order_idx]

  # here come the checks -- eventually, an option may skip this
  if (TRUE) {
    # 1. same set of observed variables?
    ov_names <- lapply(mods, function(x) {
      sort(lav_object_vnames(x))
    })
    ov <- ov_names[[1L]] # the observed variable names of the first model
    if (!all(sapply(ov_names, function(x) identical(x, ov)))) {
      lav_msg_warn(gettext(
        "some models are based on a different set of observed variables"))
    }
    ## wow FIXME: we may need to reorder the rows/columns first!!
    # COVS <- lapply(mods, function(x) slot(slot(x, "Sample"), "cov")[[1]])
    # if(!all(sapply(COVS, all.equal, COVS[[1]]))) {
    #    stop("lavaan ERROR: models must be fit to the same data")
    # }
    # 2. nested models? *different* npars?

    # TODO!

    # 3. all meanstructure?
    mean_structure <- sapply(mods, lavInspect, "meanstructure")
    if (sum(mean_structure) > 0L &&
      sum(mean_structure) < length(mean_structure)) {
      lav_msg_warn(gettext("not all models have a meanstructure"))
    }

    # 4. all converged?
    if (!all(sapply(mods, lavInspect, "converged"))) {
      lav_msg_warn(gettext("not all models converged"))
    }
  }

  mods_scaled <- unlist(lapply(mods, function(x) {
    any(c(
      "satorra.bentler", "yuan.bentler", "yuan.bentler.mplus",
      "mean.var.adjusted", "scaled.shifted"
    ) %in%
      unlist(sapply(slot(x, "test"), "[[", "test")))
  }))

  if (all(mods_scaled | ndf == 0) && any(mods_scaled)) {
    # Note: if df=0, test is not really robust, hence the above condition
    scaled <- TRUE
    # which test to choose by default?
    # i.e., not determined by method=
    scaled_list <- sapply(mods[[which(ndf > 0)[1]]]@test,
                         # first mod with df>0
                         #FIXME: ? If no mods have df > 0,
                         #         this still yields error
                         function(x) !is.null(x$scaled.test.stat))
    scaled_idx <- which(scaled_list)[[1]]
    default_test <- object@test[[scaled_idx]]$test
    if (test == "default") {
      test_1 <- default_test
    } else if (!test %in% c("satorra.bentler", "yuan.bentler",
     "yuan.bentler.mplus", "mean.var.adjusted", "scaled.shifted")) {
      lav_msg_stop(gettextf(
        "test = %s not found in object. See available tests in
        lavInspect(object, \"options\")$test.", dQuote(test)))
    } else {
      test_1 <- test
    }

    ## is the test available from all models?
    check_scaled <- unlist(lapply(mods, function(x) {
      test_1 %in% unlist(sapply(slot(x, "test"), "[[", "test"))
    }))

    if (any(!check_scaled)) {
      lav_msg_stop(gettextf(
        "test = %1$s not found in model(s): %2$s. Find available tests per model
        using lavInspect(fit, \"options\")$test.", dQuote(test),
        lav_msg_view(names(mods)[which(!check_scaled)], "none")))
    }

  } else if (!any(mods_scaled)) { # thanks to R.M. Bee to fix this
    scaled <- FALSE
    test_1 <- "standard"
    method <- "standard"
  } else {
    lav_msg_stop(gettext(
      "some models (but not all) have scaled test statistics"))
  }
  if (type %in% c("browne.residual.adf", "browne.residual.nt")) {
    scaled <- FALSE
    method <- "standard"
  }
  if (method == "standard") {
    scaled <- FALSE
  }

  # select method
  if (method == "default") {
    if (estimator == "PML") {
      method <- "mean.var.adjusted.PLRT"
    } else if (scaled) {
      if (test_1 %in% c(
        "satorra.bentler", "yuan.bentler",
        "yuan.bentler.mplus"
      )) {
        method <- "satorra.bentler.2001"
      } else {
        method <- "satorra.2000"
      }
    } else {
      # nothing to do
    }
  } else if (method == "meanvaradjustedplrt" ||
    method == "mean.var.adjusted.PLRT") {
    method <- "mean.var.adjusted.PLRT"
    stopifnot(estimator == "PML")
  } else if (method == "satorra2000") {
    method <- "satorra.2000"
  } else if (method == "satorrabentler2001") {
    method <- "satorra.bentler.2001"
  } else if (method == "satorrabentler2010") {
    method <- "satorra.bentler.2010"

    ## only option left:
  } else if (method != "standard") {
    lav_msg_stop(
      gettextf("Unknown method for scaled difference test: %s.", method))
  }

  ## in case users specify method= or test= (but still type="chisq"),
  ## make sure the arguments are consistent for scaled tests
  if (method %in% c("satorra.bentler.2001", "satorra.bentler.2010") && scaled &&
    (!test_1 %in% c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus"))) {
    lav_msg_stop(gettextf(
      "method = %s only available when models are fitted with test =
      \"satorra.bentler\", \"yuan.bentler\", or \"yuan.bentler.mplus\".",
      dQuote(method)))
  } else {
    ## method="satorra.2000" still available when test_1 != scaled.shifted
    ## Or !scaled, so nothing to do.
  }

  # check method if scaled = FALSE
  if (type == "chisq" && !scaled &&
    method %in% c(
      "mean.var.adjusted.PLRT",
      "satorra.bentler.2001",
      "satorra.2000",
      "satorra.bentler.2010"
    )) {
    lav_msg_warn(gettextf(
      "method = %s but no robust test statistics were used; switching to the
      standard chi-squared difference test", dQuote(method)))
    method <- "standard"
  }


  # which models have used a MEANSTRUCTURE?
  mods_meanstructure <- sapply(mods, function(x) {
    unlist(slot(
      slot(x, "Model"),
      "meanstructure"
    ))
  })
  if (!all(mods_meanstructure) && sum(mods_meanstructure) != 0) {
    lav_msg_stop(gettext("some models (but not all) have a meanstructure"))
  }

  # collect statistics for each model
  if (type == "chisq") {
    df_1 <- sapply(mods, function(x) slot(x, "test")[[1]]$df)
    stat_1 <- sapply(mods, function(x) slot(x, "test")[[1]]$stat)
  } else if (type == "browne.residual.nt") {
    testlist <- lapply(
      mods,
      function(x) lavTest(x, test = "browne.residual.nt")
    )
    df_1 <- sapply(testlist, function(x) x$df)
    stat_1 <- sapply(testlist, function(x) x$stat)
  } else if (type == "browne.residual.adf") {
    testlist <- lapply(
      mods,
      function(x) lavTest(x, test = "browne.residual.adf")
    )
    df_1 <- sapply(testlist, function(x) x$df)
    stat_1 <- sapply(testlist, function(x) x$stat)
  } else if (type == "cf") {
    tmp <- lapply(mods, lavTablesFitCf)
    stat_1 <- unlist(tmp)
    df_1 <- unlist(lapply(tmp, attr, "DF"))
  } else {
    lav_msg_stop(gettextf("test type unknown: %s", type))
  }

  # difference statistics
  stat_delta <- stat_delta_orig <- c(NA, diff(stat_1))
  df_delta <- df_delta_orig <- c(NA, diff(df_1))

  # check for negative values in STAT.delta
  # but with a tolerance (0.6-12)!
  if (any(stat_delta[-1] < -1 * .Machine$double.eps^(1 / 3))) {
    lav_msg_warn(gettextf(
      "Some restricted models fit better than less restricted models; either
      these models are not nested, or the less restricted model failed to reach
      a global optimum.Smallest difference = %s.", min(stat_delta[-1])))
  }

  # prepare for scaling versions
  if (method == "satorra.2000" && scaled.shifted) {
    a_delta <- b_delta <- rep(as.numeric(NA), length(stat_1))
    c_delta <- NULL
  } else if (method %in% c("satorra.bentler.2001",
           "satorra.bentler.2010", "satorra.2000")) {
    c_delta <- rep(as.numeric(NA), length(stat_1))
  }

  # correction for scaled test statistics
  if (type == "chisq" && scaled) {
    if (method == "satorra.bentler.2001") {
      # use formula from Satorra & Bentler 2001
      for (m in seq_len(length(mods) - 1L)) {
        out <- lav_test_diff_SatorraBentler2001(mods[[m]], mods[[m + 1]],
                                                # in case not @test[[2]]:
                                                test = test_1)
        stat_delta[m + 1] <- out$T.delta
        df_delta[m + 1] <- out$df.delta
        c_delta[m + 1] <- out$scaling.factor
      }
    } else if (method == "mean.var.adjusted.PLRT") {
      for (m in seq_len(length(mods) - 1L)) {
        out <- lav_pml_test_plrt(mods[[m]], mods[[m + 1]])
        stat_delta[m + 1] <- out$FSMA.PLRT
        df_delta[m + 1] <- out$adj.df
      }
    } else if (method == "satorra.bentler.2010") {
      for (m in seq_len(length(mods) - 1L)) {
        out <- lav_test_diff_SatorraBentler2010(mods[[m]], mods[[m + 1]],
          test = test_1, # in case not @test[[2]]
          H1 = FALSE
        ) # must be F

        stat_delta[m + 1] <- out$T.delta
        df_delta[m + 1] <- out$df.delta
        c_delta[m + 1] <- out$scaling.factor
      }
    } else if (method == "satorra.2000") {
      for (m in seq_len(length(mods) - 1L)) {
        if (test_1 %in% c(
          "satorra.bentler", "yuan.bentler",
          "yuan.bentler.mplus"
        )) {
          satterthwaite <- FALSE
        } else {
          satterthwaite <- TRUE
         }
        out <- lav_test_diff_satorra2000(mods[[m]], mods[[m + 1]],
          h1 = TRUE,
          satterthwaite = satterthwaite,
          scaled_shifted = scaled.shifted,
          a_method = A.method
        )
        stat_delta[m + 1] <- out$T.delta
        df_delta[m + 1] <- out$df.delta
        if (scaled.shifted) {
          a_delta[m + 1] <- out$a
          b_delta[m + 1] <- out$b
        } else {
          c_delta[m + 1] <- out$scaling.factor
        }
      }
    }
  }

  # check if scaled diff failed somehow
  if (scaled &&
     ((method %in% c("satorra.bentler.2001", "satorra.bentler.2010") &&
          is.na(out$scaling.factor)) ||
         (method == "satorra.2000" && scaled.shifted && is.na(out$a)) ||
         (method == "satorra.2000" && !scaled.shifted &&
          is.na(out$scaling.factor)))
     ) {
    scaled <- FALSE
  }


  # unname
  stat_delta <- unname(stat_delta)
  # zap small values (for near-zero values) (anova class does not use rounding)
  stat_delta <- round(stat_delta, 10)
  df_delta <- unname(df_delta)
  stat_delta_orig <- unname(stat_delta_orig)
  df_delta_orig <- unname(df_delta_orig)
  if (scaled && !is.null(c_delta)) {
    c_delta <- unname(c_delta)
  }

  # Pvalue
  pvalue_delta <- pchisq(stat_delta, df_delta, lower.tail = FALSE)

  # new in 0.6-13: RMSEA (RMSEA.D or RDR)
  if (object@Options$missing == "listwise") {
    if (scaled && !is.null(c_delta)) {
      c_hat <- c_delta[-1]
    } else {
      c_hat <- rep(1, length(stat_delta_orig) - 1L)
    }
    rmsea_delta <- c(NA, lav_fit_rmsea(
      X2 = stat_delta_orig[-1],
      df = df_delta_orig[-1],
      N = ntotal,
      G = ngroups,
      c.hat = c_hat
    ))
  }

  # AIC/BIC
  aic <- bic <- rep(NA, length(mods))
  if (estimator == "ML") {
    aic <- sapply(mods, FUN = AIC)
    bic <- sapply(mods, FUN = BIC)
  } else if (estimator == "PML") {
    out_1 <- lapply(mods, lav_pml_object_aic_bic)
    aic <- sapply(out_1, "[[", "PL_AIC")
    bic <- sapply(out_1, "[[", "PL_BIC")
  }

  if (estimator == "PML") {
    val <- data.frame(
      Df = df_1,
      PL_AIC = aic,
      PL_BIC = bic,
      Chisq = stat_1,
      "Chisq diff" = stat_delta,
      "Df diff" = df_delta,
      "Pr(>Chisq)" = pvalue_delta,
      row.names = names(mods),
      check.names = FALSE
    )
  } else if (object@Options$missing == "listwise") {
    val <- data.frame(
      Df = df_1,
      AIC = aic,
      BIC = bic,
      Chisq = stat_1,
      "Chisq diff" = stat_delta,
      "RMSEA" = rmsea_delta,
      "Df diff" = df_delta,
      "Pr(>Chisq)" = pvalue_delta,
      row.names = names(mods),
      check.names = FALSE
    )
  } else {
    val <- data.frame(
      Df = df_1,
      AIC = aic,
      BIC = bic,
      Chisq = stat_1,
      "Chisq diff" = stat_delta,
      #"RMSEA" = RMSEA.delta, # if missing, not yet...
      "Df diff" = df_delta,
      "Pr(>Chisq)" = pvalue_delta,
      row.names = names(mods),
      check.names = FALSE
    )
  }

  # catch Df.delta == 0 cases (reported by Florian Zsok in Zurich)
  # but only if there are no inequality constraints! (0.6-1)
  idx <- which(val[, "Df diff"] == 0)
  if (length(idx) > 0L) {
    # remove models with inequality constraints
    ineq_idx <- which(sapply(lapply(mods, function(x) {
      slot(slot(x, "Model"), "x.cin.idx")
    }), length) > 0L)
    rm_idx <- which(idx %in% ineq_idx)
    if (length(rm_idx) > 0L) {
      idx <- idx[-rm_idx]
    }
  }
  if (length(idx) > 0L) {
    val[idx, "Pr(>Chisq)"] <- as.numeric(NA)
    lav_msg_warn(gettext("some models have the same degrees of freedom"))
  }

  if (type == "chisq") {
    if (scaled) {
      txt <- paste("The ", dQuote("Chisq"), " column contains standard ",
        "test statistics, not the robust test that should be ",
        "reported per model. A robust difference test is a ",
        "function of two standard (not robust) statistics.",
        sep = ""
      )
      attr(val, "heading") <-
        paste("\nScaled Chi-Squared Difference Test (method = ",
          dQuote(method), ")\n\n",
          lav_msg(paste("lavaan NOTE:", txt), showheader = TRUE),
          "\n",
          sep = ""
        )
      if (method == "satorra.2000" && scaled.shifted) {
        attr(val, "scale") <- a_delta
        attr(val, "shift") <- b_delta
      } else if (method %in% c("satorra.bentler.2001",
                "satorra.bentler.2010", "satorra.2000")) {
        attr(val, "scale") <- c_delta
      }
    } else {
      attr(val, "heading") <- "\nChi-Squared Difference Test\n"
    }
  } else if (type == "browne.residual.adf") {
    attr(val, "heading") <-
      "\nChi-Squared Difference Test based on Browne's residual (ADF) Test\n"
  } else if (type == "browne.residual.nt") {
    attr(val, "heading") <-
      "\nChi-Squared Difference Test based on Browne's residual (NT) Test\n"
  } else if (type == "cf") {
    colnames(val)[c(3, 4)] <- c("Cf", "Cf diff")
    attr(val, "heading") <- "\nCf Difference Test\n"
  }
  class(val) <- c("anova", class(val))

  val
}


# anova table for a single model
lav_test_lrt_single_model <- function(object, method = "default",
                                      test = "default", type = "Chisq") {
  estimator <- object@Options$estimator

  aic <- bic <- c(NA, NA)
  if (estimator == "ML") {
    aic <- c(NA, AIC(object))
    bic <- c(NA, BIC(object))
  }

  ## determine which @test element
  tn <- names(object@test)
  if (is.null(tn)) {
    tn <- "standard" # for lavaan <0.6 objects
  }
  if (length(tn) == 1L) {
    test_1 <- 1L # only choice

    ## More than 1.  Cycle through possible user specifications:
  } else if (method[1] == "standard") {
    test_1 <- 1L
  } else if (grepl(pattern = "browne", x = type) && type %in% tn) {
    test_1 <- type
  } else if (test %in% tn) {
    test_1 <- test
  } else {
    ## Nothing explicitly (or validly) requested.
    ## But there is > 1 test, so take the second element (old default)
    test_1 <- 2L
  }

  ## anova table
  val <- data.frame(
    Df = c(0, object@test[[test_1]]$df),
    AIC = aic,
    BIC = bic,
    Chisq = c(0, object@test[[test_1]]$stat),
    "Chisq diff" = c(NA, object@test[[test_1]]$stat),
    "Df diff" = c(NA, object@test[[test_1]]$df),
    "Pr(>Chisq)" = c(NA, object@test[[test_1]]$pvalue),
    row.names = c("Saturated", "Model"),
    check.names = FALSE
  )
  ## scale/shift attributes
  if (!is.null(object@test[[test_1]]$scaling.factor)) {
    attr(val, "scale") <- c(NA, object@test[[test_1]]$scaling.factor)
  }
  if (!is.null(object@test[[test_1]]$shift.parameter)) {
    attr(val, "shift") <- c(NA, object@test[[test_1]]$shift.parameter)
  }

  ## heading
  if (grepl(pattern = "browne", x = test_1)) {
    attr(val, "heading") <- object@test[[test_1]]$label

  } else if (test_1 == 1L) {
    attr(val, "heading") <- "Chi-Squared Test Statistic (unscaled)\n"

  } else {
    label <- object@test[[test_1]]$label
    attr(val, "heading") <- paste0("Chi-Squared Test Statistic (scaled",
                                   ifelse(test_1 == "scaled.shifted",
                                          yes = " and shifted)", no = ")"),
                                   ifelse(is.null(label),
                                          yes = "\n", no = paste("\n ", label)),
                                   "\n")
  }

  class(val) <- c("anova", class(val))

  val
}
