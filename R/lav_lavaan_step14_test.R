lav_lavaan_step14_test <- function(lavoptions = NULL,
                                   lavmodel = NULL,
                                   lavsamplestats = NULL,
                                   lavdata = NULL,
                                   lavpartable = NULL,
                                   lavcache = NULL,
                                   lavimplied = NULL,
                                   lavh1 = NULL,
                                   x = NULL,
                                   VCOV = NULL, # nolint
                                   lavloglik = NULL) {
  # # # # # # # # # # #
  # #  14. lavtest # #
  # # # # # # # # # # #

  # if lavoptions$test != "none" and x converged
  #   compute lavtest via lav_model_test(...)
  # else
  #   lavtest <- list(list(test = "none", stat = NA,
  #                        stat.group = rep(NA, lavdata@ngroups),
  #                        df = NA, refdistr = "unknown", pvalue = NA))
  lavtest <- NULL
  if (!(length(lavoptions$test) == 1L && lavoptions$test == "none") &&
    attr(x, "converged")) {
    if (lavoptions$verbose) {
      cat("computing TEST for test(s) =", lavoptions$test, "...")
    }
    lavtest <- lav_model_test(
      lavmodel = lavmodel,
      lavpartable = lavpartable,
      lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions,
      x = x,
      VCOV = VCOV,
      lavdata = lavdata,
      lavcache = lavcache,
      lavloglik = lavloglik
    )
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  } else {
    lavtest <- list(list(
      test = "none", stat = NA,
      stat.group = rep(NA, lavdata@ngroups), df = NA,
      refdistr = "unknown", pvalue = NA
    ))
  }

  lavtest
}

lav_lavaan_step14_fit <- function(lavpartable = NULL,
                                  lavmodel = NULL,
                                  lavimplied = NULL,
                                  x = NULL,
                                  VCOV = NULL, # nolint
                                  lavtest = NULL) {
  # # # # # # # # # # # #
  # #  14bis. lavfit  # # -> remove if the offending packages are fixed!!
  # # # # # # # # # # # #

  lavfit <- lav_model_fit(
    lavpartable = lavpartable,
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    x = x,
    VCOV = VCOV,
    TEST = lavtest
  )
  # lavfit <- new("Fit")

  lavfit
}
