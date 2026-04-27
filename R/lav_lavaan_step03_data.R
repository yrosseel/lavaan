lav_lavaan_step03_data <- function(slot_data = NULL, # nolint
                                   lavoptions = NULL,
                                   ov_names = NULL,
                                   ov_names_y = NULL,
                                   group = NULL,
                                   data = NULL,
                                   cluster = NULL,
                                   ov_names_x = NULL,
                                   ov_names_l = NULL,
                                   ordered = NULL,
                                   sampling_weights = NULL,
                                   sample_cov = NULL,
                                   sample_mean = NULL,
                                   sample_th = NULL,
                                   sample_nobs = NULL,
                                   slot_par_table = NULL, # nolint
                                   ngroups = NULL,
                                   dotdotdot = NULL,
                                   flat_model = NULL,
                                   model = NULL,
                                   nacov = NULL, # nolint
                                   wls_v = NULL) { # nolint
  # # # # # # # # # # #
  # #  3. lavdata  # #
  # # # # # # # # # # #

  # if slotData not null
  #   copy slotData to lavdata
  # else
  #   create lavdata via function lav_lavdata, setting ov.names to ov.names.y
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

  if (!is.null(slot_data)) {
    lavdata <- slot_data
  } else {
    if (lav_verbose()) {
      cat("lavdata            ...")
    }
    # FIXME: ov.names should always contain both y and x!
    tmp_ov_names <- if (lavoptions$conditional.x) {
      ov_names_y
    } else {
      ov_names
    }
    lavdata <- lav_lavdata(
      data = data,
      group = group,
      cluster = cluster,
      ov.names = tmp_ov_names,
      ov.names.x = ov_names_x,
      ov.names.l = ov_names_l,
      ordered = ordered,
      sampling.weights = sampling_weights,
      sample.cov = sample_cov,
      sample.mean = sample_mean,
      sample.th = sample_th,
      sample.nobs = sample_nobs,
      lavoptions = lavoptions
    )

    if (lav_verbose()) {
      cat(" done.\n")
    }
  }
  # what have we learned from the data?
  if (lavdata@data.type == "none") {
    lavoptions$do.fit <- FALSE
    # check if 'model' was a fitted parameter table
    lavoptions$start <- ifelse(is.null(flat_model$est), "simple", "est")
    lavoptions$se <- "none"
    lavoptions$test <- "none"
  } else if (lavdata@data.type == "moment") {
    # check user-specified options first
    if (!is.null(dotdotdot$estimator)) {
      if (is.list(dotdotdot$estimator)) {
        estimator <- toupper(dotdotdot$estimator[[1]])
      } else {
        estimator <- toupper(dotdotdot$estimator)
      }
      if (any(estimator == c(
        "MLM", "MLMV", "MLR", "MLR",
        "ULSM", "ULSMV", "ULSMVS"
      )) &&
        is.null(nacov)) {
        lav_msg_stop(gettextf(
          "estimator %s requires full data or user-provided NACOV", estimator))
      } else if (any(estimator == c(
        "WLS", "WLSM", "WLSMV",
        "WLSMVS", "DWLS"
      )) &&
        is.null(wls_v)) {
        lav_msg_stop(gettextf(
          "estimator %s requires full data or user-provided WLS.V and NACOV",
          estimator))
      }
    }
    # catch here some options that will not work with moments
    if (lavoptions$se == "bootstrap") {
      lav_msg_stop(gettext("bootstrapping requires full data"))
    }
    # more needed?
  }
  # sanity check
  if (!is.null(slot_par_table) || inherits(model, "lavaan")) {
    if (ngroups != lavdata@ngroups) {
      lav_msg_stop(gettext(
        "mismatch between number of groups in data
        and number of groups in model."))
    }
  }
  if (lav_verbose()) {
    print(lavdata)
  }
  if (lav_debug()) {
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
