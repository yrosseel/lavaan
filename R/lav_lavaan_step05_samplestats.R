lav_step05_samp <- function(slot_sample_stats = NULL,
                                          lavdata = NULL,
                                          lavoptions = NULL,
                                          wls_v = NULL,
                                          nacov = NULL,
                                          sample_cov = NULL,
                                          sample_mean = NULL,
                                          sample_th = NULL,
                                          sample_nobs = NULL,
                                          ov_names = NULL,
                                          ov_names_x = NULL,
                                          lavpartable = NULL) {
  # # # # # # # # # # # # # #
  # #  5. lavsamplestats  # #
  # # # # # # # # # # # # # #

  # if slot_sample_stats not NULL
  #   copy to lavsamplestats
  # else
  #   if lavdata@data.type == "full"
  #     compute lavsamplestats via lav_samplestats_from_data
  #   else
  #     if lavdata@data.type == "moment"
  #       if lavoptions$meanstructure TRUE but sample_mean is NULL:
  #         ** warning **
  #       compute lavsamplestats via lav_samplestats_from_moments
  #     else
  #       create lavsamplestats object (type lavSampleStats) with data from
  #         lavdata and lavpta

  if (!is.null(slot_sample_stats)) {
    lavsamplestats <- slot_sample_stats
  } else if (lavdata@data.type == "full") {
    if (lav_verbose()) {
      cat("lavsamplestats     ...")
    }
    lavsamplestats <- lav_samp_from_data(
      lavdata       = lavdata,
      lavoptions    = lavoptions,
      wls_v         = wls_v,
      nacov         = nacov
    )
    if (lav_verbose()) {
      cat(" done.\n")
    }
  } else if (lavdata@data.type == "moment") {
    if (lav_verbose()) {
      cat("lavsamplestats ...")
    }
    # check if we have sample_mean and meanstructure = TRUE
    # exception: a FULL correlation structure implies standardized
    # (centered) variables, so zero means are the exact metric convention,
    # not an assumption -- no warning needed (typically triggered by
    # conditional.x = TRUE, which forces meanstructure = TRUE); for a
    # partial structure the warning stays: the non-listed variables keep
    # their original (nonzero) means, which only sample_mean= can supply
    correlation_full <- isTRUE(lavoptions$correlation) &&
      length(lavoptions$.correlation.ov) == 0L
    if (lavoptions$meanstructure && is.null(sample_mean) &&
        !correlation_full) {
      lav_msg_warn(
        gettext("sample_mean= argument is missing, but model contains
                mean/intercept parameters."))
    }
    lavsamplestats <- lav_samp_from_moments(
      sample_cov    = sample_cov,
      sample_mean   = sample_mean,
      sample_th     = sample_th,
      sample_nobs   = sample_nobs,
      ov_names      = ov_names,
      ov_names_x    = ov_names_x,
      wls_v         = wls_v,
      nacov         = nacov,
      lavoptions    = lavoptions
    )
    if (lav_verbose()) {
      cat(" done.\n")
    }
  } else {
    # no data
    lavsamplestats <- new("lavSampleStats",
      ngroups = lavdata@ngroups,
      nobs = as.list(rep(0L, lavdata@ngroups)),
      cov.x = vector("list", length = lavdata@ngroups),
      mean.x = vector("list", length = lavdata@ngroups),
      th.idx = attr(lavpartable, "th.idx"),
      missing.flag = FALSE
    )
  }

  if (lav_debug()) {
    print(str(lavsamplestats))
  }

  lavsamplestats
}
