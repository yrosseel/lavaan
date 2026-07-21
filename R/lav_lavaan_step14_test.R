lav_step14_test <- function(lavoptions = NULL,
                                   lavmodel = NULL,
                                   lavsamplestats = NULL,
                                   lavdata = NULL,
                                   lavpartable = NULL,
                                   lavcache = NULL,
                                   lavimplied = NULL,
                                   lavh1 = NULL,
                                   x = NULL,
                                   vcov = NULL,
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
  # internal rescaling (badly scaled variables): the test statistics
  # were computed by the inner fit in the well-scaled metric; they are
  # invariant under the rescaling, so use them as-is
  if (!is.null(lavoptions$.rescale.test)) {
    return(lavoptions$.rescale.test)
  }

  lavtest <- NULL
  if (!(length(lavoptions$test) == 1L && lavoptions$test == "none") &&
    attr(x, "converged")) {
    if (lav_verbose()) {
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
      vcov_1 = vcov,
      lavdata = lavdata,
      lavcache = lavcache,
      lavloglik = lavloglik
    )
    if (lav_verbose()) {
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

lav_step14_fit <- function(lavpartable = NULL,
                                  lavmodel = NULL,
                                  lavimplied = NULL,
                                  x = NULL,
                                  vcov = NULL,
                                  lavtest = NULL) {
  # # # # # # # # # # # #
  # #  14bis. lavfit  # # -> remove if the offending packages are fixed!!
  # # # # # # # # # # # #

  lavfit <- lav_model_fit(
    lavpartable = lavpartable,
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    x = x,
    vcov_1 = vcov,
    test = lavtest
  )
  # lavfit <- new("Fit")

  lavfit
}
