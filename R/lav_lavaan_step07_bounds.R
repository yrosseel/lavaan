lav_lavaan_step07_bounds <- function(lavoptions, lavh1, lavdata, lavsamplestats, lavpartable) {
  
  # # # # # # # # # # # # # # #
  # #  7. parameter bounds # # 
  # # # # # # # # # # # # # # #
  
  # automatic bounds (new in 0.6-6)
  if (!is.null(lavoptions$optim.bounds) ||
      length(lavoptions$optim.bounds$lower) > 0L ||
      length(lavoptions$optim.bounds$upper) > 0L) {
    if (lavoptions$verbose) {
      cat("lavpartable bounds ...")
    }
    lavpartable <- lav_partable_add_bounds(partable = lavpartable,
                                           lavh1 = lavh1, lavdata = lavdata, lavsamplestats = lavsamplestats,
                                           lavoptions = lavoptions)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }
  return(lavpartable)
}