lav_lavaan_step06_h1 <- function(sloth1 = NULL,
                                 lavoptions = NULL,
                                 lavsamplestats = NULL,
                                 lavdata = NULL,
                                 lavpartable = NULL) {
  # # # # # # # # #
  # #  6. lavh1 # #
  # # # # # # # # #

  # if sloth1 not NULL
  #   copy to lavh1
  # else
  #   if lavoptions$h1 TRUE
  #     if length(lavsamplestats$ntotal) > 0
  #       compute lavh1 via lav_h1_implied_logl
  #   else
  #     check lavoptions$h1 is logical, if not *** error ***

  if (!is.null(sloth1)) {
    lavh1 <- sloth1
  } else {
    lavh1 <- list()
    if (is.logical(lavoptions$h1) && lavoptions$h1) {
      if (length(lavsamplestats@ntotal) > 0L ||
	      (!is.null(lavoptions$samplestats) && !lavoptions$samplestats)) {
		  # lavsamplestats filled in
        if (lavoptions$verbose) {
          cat("lavh1              ... start:\n")
        }

        # implied h1 statistics and logl (if available)
        lavh1 <- lav_h1_implied_logl(
          lavdata = lavdata,
          lavsamplestats = lavsamplestats,
          lavpartable = lavpartable,
          lavoptions = lavoptions
        )
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
        lav_msg_stop(gettext("argument `h1' must be logical (for now)"))
      }
      # TODO: allow h1 to be either a model syntax, a parameter table,
      # or a fitted lavaan object
    }
  }

  lavh1
}
