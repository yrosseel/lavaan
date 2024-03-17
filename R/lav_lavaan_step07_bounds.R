lav_lavaan_step07_bounds <- function(lavoptions = NULL,
                                     lavh1 = NULL,
                                     lavdata = NULL,
                                     lavsamplestats = NULL,
                                     lavpartable = NULL) {
  # # # # # # # # # # # # # # #
  # #  7. parameter bounds  # #
  # # # # # # # # # # # # # # #

  # if lavoptions$optim.bounds not NULL and its members lower and upper
  #   have length > 0L
  #   modify lavpartable via lav_partable_add_bounds

  # automatic bounds (new in 0.6-6)
  if (!is.null(lavoptions$optim.bounds) ||
    length(lavoptions$optim.bounds$lower) > 0L ||
    length(lavoptions$optim.bounds$upper) > 0L) {
    if (lavoptions$verbose) {
      cat("lavpartable bounds ...")
    }
    lavpartable <- lav_partable_add_bounds(
      partable = lavpartable,
      lavh1 = lavh1, lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }

  lavpartable
}
