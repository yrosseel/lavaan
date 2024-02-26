lav_lavaan_step05_samplestats <- function(slotSampleStats, lavdata, lavoptions, WLS.V, NACOV,
                                          sample.cov, sample.mean, sample.th, sample.nobs,
                                          ov.names, ov.names.x, lavpta) {
  
  # # # # # # # # # # # # # #
  # #  5. lavsamplestats # # 
  # # # # # # # # # # # # # #
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
      NACOV         = NACOV)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  } else if (lavdata@data.type == "moment") {
    if (lavoptions$verbose) {
      cat("lavsamplestats ...")
    }
    # check if we have sample.mean and meanstructure = TRUE
    if (lavoptions$meanstructure && is.null(sample.mean)) {
      txt <- "sample.mean= argument is missing, but model contains mean/intercept parameters."
      warning(lav_txt2message(txt))
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
      lavoptions    = lavoptions)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  } else {
    # no data
    lavsamplestats <- new("lavSampleStats", ngroups = lavdata@ngroups,
                          nobs = as.list(rep(0L, lavdata@ngroups)),
                          cov.x = vector("list", length = lavdata@ngroups),
                          mean.x = vector("list", length = lavdata@ngroups),
                          th.idx = lavpta$th.idx,
                          missing.flag = FALSE)
  }
  if (lavoptions$debug) {
    print(str(lavsamplestats))
  }
  return(lavsamplestats)
}