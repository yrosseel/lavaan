lav_lavaan_step05_samplestats <- function(slotSampleStats = NULL, # nolint
                                          lavdata = NULL,
                                          lavoptions = NULL,
                                          WLS.V = NULL, # nolint
                                          NACOV = NULL, # nolint
                                          sample.cov = NULL,
                                          sample.mean = NULL,
                                          sample.th = NULL,
                                          sample.nobs = NULL,
                                          ov.names = NULL,
                                          ov.names.x = NULL,
                                          lavpartable = NULL) {
  # # # # # # # # # # # # # #
  # #  5. lavsamplestats  # #
  # # # # # # # # # # # # # #

  # if slotSampleStats not NULL
  #   copy to lavsamplestats
  # else
  #   if lavdata@data.type == "full"
  #     compute lavsamplestats via lav_samplestats_from_data
  #   else
  #     if lavdata@data.type == "moment"
  #       if lavoptions$meanstructure TRUE but sample.mean is NULL:
  #         ** warning **
  #       compute lavsamplestats via lav_samplestats_from_moments
  #     else
  #       create lavsamplestats object (type lavSampleStats) with data from
  #         lavdata and lavpta

  if (!is.null(slotSampleStats)) {
    lavsamplestats <- slotSampleStats
  } else if (lavdata@data.type == "full") {
    if (lavoptions$verbose) {
      cat("lavsamplestats     ...")
    }
    lavsamplestats <- lav_samplestats_from_data(
      lavdata       = lavdata,
      lavoptions    = lavoptions,
      WLS.V         = WLS.V,
      NACOV         = NACOV
    )
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  } else if (lavdata@data.type == "moment") {
    if (lavoptions$verbose) {
      cat("lavsamplestats ...")
    }
    # check if we have sample.mean and meanstructure = TRUE
    if (lavoptions$meanstructure && is.null(sample.mean)) {
      lav_msg_warn(
        gettext("sample.mean= argument is missing, but model contains
                mean/intercept parameters."))
    }
    lavsamplestats <- lav_samplestats_from_moments(
      sample.cov    = sample.cov,
      sample.mean   = sample.mean,
      sample.th     = sample.th,
      sample.nobs   = sample.nobs,
      ov.names      = ov.names,
      ov.names.x    = ov.names.x,
      WLS.V         = WLS.V,
      NACOV         = NACOV,
      lavoptions    = lavoptions
    )
    if (lavoptions$verbose) {
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

  if (lavoptions$debug) {
    print(str(lavsamplestats))
  }

  lavsamplestats
}
