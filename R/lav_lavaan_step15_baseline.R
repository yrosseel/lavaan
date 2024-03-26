lav_lavaan_step15_baseline <- function(lavoptions = NULL,
                                       lavsamplestats = NULL,
                                       lavdata = NULL,
                                       lavcache = NULL,
                                       lavh1 = NULL,
                                       lavpartable = NULL) {
  # # # # # # # # # # #
  # #  15. baseline # #  (since 0.6-5)
  # # # # # # # # # # #

  # if options$do.fit and options$test not "none" and options$baseline = TRUE
  #   try fit.indep <- lav_object_independence(...)
  #   if not succesfull or not converged
  #     ** warning **
  #     lavbaseline < list()
  #   else
  #     lavbaseline <- list with partable and test of fit.indep
  lavbaseline <- list()
  if (lavoptions$do.fit &&
    !("none" %in% lavoptions$test) &&
    is.logical(lavoptions$baseline) && lavoptions$baseline) {
    if (lavoptions$verbose) {
      cat("lavbaseline ...")
    }
    fit.indep <- try(lav_object_independence(
      object = NULL,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache,
      lavoptions = lavoptions,
      lavpartable = lavpartable,
      lavh1 = lavh1
    ), silent = TRUE)
    if (inherits(fit.indep, "try-error") ||
      !fit.indep@optim$converged) {
      if (lavoptions$warn) {
        lav_msg_warn(gettext("estimation of the baseline model failed."))
      }
      lavbaseline <- list()
      if (lavoptions$verbose) {
        cat(" FAILED.\n")
      }
    } else {
      # store relevant information
      lavbaseline <- list(
        partable = fit.indep@ParTable,
        test = fit.indep@test
      )
      if (lavoptions$verbose) {
        cat(" done.\n")
      }
    }
  }

  lavbaseline
}
