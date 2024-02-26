lav_lavaan_step15_baseline <- function(lavoptions, lavsamplestats, lavdata,
                                     lavcache, lavh1, lavpta) {
  # # # # # # # # # # # 
  # #  15. baseline # #  (since 0.6-5)
  # # # # # # # # # # # 
  lavbaseline <- list()
  if (lavoptions$do.fit &&
      !("none" %in% lavoptions$test) &&
      is.logical(lavoptions$baseline) && lavoptions$baseline) {
    if (lavoptions$verbose) {
      cat("lavbaseline ...")
    }
    fit.indep <- try(lav_object_independence(object = NULL,
                                             lavsamplestats = lavsamplestats,
                                             lavdata        = lavdata,
                                             lavcache       = lavcache,
                                             lavoptions     = lavoptions,
                                             lavpta         = lavpta,
                                             lavh1          = lavh1), silent = TRUE)
    if (inherits(fit.indep, "try-error") ||
        !fit.indep@optim$converged) {
      if (lavoptions$warn) {
        warning("lavaan WARNING: estimation of the baseline model failed.")
      }
      lavbaseline <- list()
      if (lavoptions$verbose) {
        cat(" FAILED.\n")
      }
    } else {
      # store relevant information
      lavbaseline <- list(partable = fit.indep@ParTable,
                          test     = fit.indep@test)
      if (lavoptions$verbose) {
        cat(" done.\n")
      }
    }
  }
  return(lavbaseline)
}