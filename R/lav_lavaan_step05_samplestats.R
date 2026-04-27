lav_lavaan_step05_samplestats <- function(slot_sample_stats = NULL, # nolint
                                          lavdata = NULL,
                                          lavoptions = NULL,
                                          wls_v = NULL, # nolint
                                          nacov = NULL, # nolint
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

  if (!is.null(slot_sample_stats)) {
    lavsamplestats <- slot_sample_stats
  } else if (lavdata@data.type == "full") {
    if (lav_verbose()) {
      cat("lavsamplestats     ...")
    }
    lavsamplestats <- lav_samplestats_from_data(
      lavdata       = lavdata,
      lavoptions    = lavoptions,
      WLS.V         = wls_v,
      NACOV         = nacov
    )
    if (lav_verbose()) {
      cat(" done.\n")
    }
  } else if (lavdata@data.type == "moment") {
    if (lav_verbose()) {
      cat("lavsamplestats ...")
    }
    # check if we have sample.mean and meanstructure = TRUE
    if (lavoptions$meanstructure && is.null(sample_mean)) {
      lav_msg_warn(
        gettext("sample.mean= argument is missing, but model contains
                mean/intercept parameters."))
    }
    lavsamplestats <- lav_samplestats_from_moments(
      sample.cov    = sample_cov,
      sample.mean   = sample_mean,
      sample.th     = sample_th,
      sample.nobs   = sample_nobs,
      ov.names      = ov_names,
      ov.names.x    = ov_names_x,
      WLS.V         = wls_v,
      NACOV         = nacov,
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
