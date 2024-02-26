lav_lavaan_step06_h1 <- function(sloth1, lavoptions, lavsamplestats, lavdata, lavpta) {
  # # # # # # # # # 
  # #  6. lavh1 # # 
  # # # # # # # # # 
  if (!is.null(sloth1)) {
    lavh1 <- sloth1
  } else {
    lavh1 <- list()
    if (is.logical(lavoptions$h1) && lavoptions$h1) {
      if (length(lavsamplestats@ntotal) > 0L) { # lavsamplestats filled in
        if (lavoptions$verbose) {
          cat("lavh1              ... start:\n")
        }
        
        # implied h1 statistics and logl (if available)
        lavh1 <- lav_h1_implied_logl(lavdata        = lavdata,
                                     lavsamplestats = lavsamplestats,
                                     lavpta         = lavpta,
                                     lavoptions     = lavoptions)
        if (lavoptions$debug) {
          print(lavh1)
        }
        if (lavoptions$verbose) {
          cat("lavh1              ... done.\n")
        }
      } else {
        # do nothing for now
      }
    } else {
      if (!is.logical(lavoptions$h1)) {
        stop("lavaan ERROR: argument `h1' must be logical (for now)")
      }
      # TODO: allow h1 to be either a model syntax, a parameter table,
      # or a fitted lavaan object
    }
  }
  return(lavh1)
}