lav_lavaan_step03_data <- function(slotData = NULL, # nolint
                                   lavoptions = NULL,
                                   ov.names = NULL,
                                   ov.names.y = NULL,
                                   group = NULL,
                                   data = NULL,
                                   cluster = NULL,
                                   ov.names.x = NULL,
                                   ov.names.l = NULL,
                                   ordered = NULL,
                                   sampling.weights = NULL,
                                   sample.cov = NULL,
                                   sample.mean = NULL,
                                   sample.th = NULL,
                                   sample.nobs = NULL,
                                   slotParTable = NULL, # nolint
                                   ngroups = NULL,
                                   dotdotdot = NULL,
                                   flat.model = NULL,
                                   model = NULL,
                                   NACOV = NULL, # nolint
                                   WLS.V = NULL) { # nolint
  # # # # # # # # # # #
  # #  3. lavdata  # #
  # # # # # # # # # # #

  # if slotData not null
  #   copy slotData to lavdata
  # else
  #   create lavdata via function lavData, setting ov.names to ov.names.y
  #     if lavoptions$conditional.x
  # if lavdata$data.type is "none"
  #   set lavoptions$do.fit to FALSE
  #   if flat.model$est not null set lavoptions$start to "est", else set
  #     it to "simple"
  #   set lavoptions$se and lavoptions$test to "none"
  # else
  #   if lavdata$data.type is "moment"
  #     if estimator one of MLM, MLMV, MLR, ULSM, ULSMV, ULSMVS and NACOV
  #       is NULL: *** error ***
  #     if estimator one of WLS, WLSM, WLSMV, WLSMVS, DWLS and WLS.V is
  #       NULL: *** error ***
  #     if lavoptions$se = bootstrap: *** error ***
  # if slotPartable not NULL and model is lavaan-object, check equality
  #   ngroups and lavdata$ngroups
  #                                       --> *** error *** if not

  if (!is.null(slotData)) {
    lavdata <- slotData
  } else {
    if (lavoptions$verbose) {
      cat("lavdata            ...")
    }
    # FIXME: ov.names should always contain both y and x!
    tmp.ov.names <- if (lavoptions$conditional.x) {
      ov.names.y
    } else {
      ov.names
    }
    lavdata <- lavData(
      data = data,
      group = group,
      cluster = cluster,
      ov.names = tmp.ov.names,
      ov.names.x = ov.names.x,
      ov.names.l = ov.names.l,
      ordered = ordered,
      sampling.weights = sampling.weights,
      sample.cov = sample.cov,
      sample.mean = sample.mean,
      sample.th = sample.th,
      sample.nobs = sample.nobs,
      lavoptions = lavoptions
    )

    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }
  # what have we learned from the data?
  if (lavdata@data.type == "none") {
    lavoptions$do.fit <- FALSE
    # check if 'model' was a fitted parameter table
    lavoptions$start <- ifelse(is.null(flat.model$est), "simple", "est")
    lavoptions$se <- "none"
    lavoptions$test <- "none"
  } else if (lavdata@data.type == "moment") {
    # check user-specified options first
    if (!is.null(dotdotdot$estimator)) {
      if (any(dotdotdot$estimator == c(
        "MLM", "MLMV", "MLR", "MLR",
        "ULSM", "ULSMV", "ULSMVS"
      )) &&
        is.null(NACOV)) {
        lav_msg_stop(gettextf(
          "estimator %s requires full data or user-provided NACOV",
          dotdotdot$estimator))
      } else if (any(dotdotdot$estimator == c(
        "WLS", "WLSM", "WLSMV",
        "WLSMVS", "DWLS"
      )) &&
        is.null(WLS.V)) {
        lav_msg_stop(gettextf(
          "estimator %s requires full data or user-provided WLS.V and NACOV",
          dotdotdot$estimator))
      }
    }
    # catch here some options that will not work with moments
    if (lavoptions$se == "bootstrap") {
      lav_msg_stop(gettext("bootstrapping requires full data"))
    }
    # more needed?
  }
  # sanity check
  if (!is.null(slotParTable) || inherits(model, "lavaan")) {
    if (ngroups != lavdata@ngroups) {
      lav_msg_stop(gettext(
        "mismatch between number of groups in data
        and number of groups in model."))
    }
  }
  if (lavoptions$verbose) {
    print(lavdata)
  }
  if (lavoptions$debug) {
    print(str(lavdata))
  }
  # if lavdata@nlevels > 1L, adapt start option (for now)
  # until we figure out how to handle groups+blocks
  # if(lavdata@nlevels > 1L) {
  #   lavoptions$start <- "simple"
  # }

  list(
    lavdata    = lavdata,
    lavoptions = lavoptions
  )
}
