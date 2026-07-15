lav_step03_data <- function(slot_data = NULL,
                                   lavoptions = NULL,
                                   ov_names = NULL,
                                   ov_names_y = NULL,
                                   group = NULL,
                                   data = NULL,
                                   cluster = NULL,
                                   ov_names_x = NULL,
                                   ov_names_l = NULL,
                                   ov_names_aux = NULL,
                                   ordered = NULL,
                                   sampling_weights = NULL,
                                   sample_cov = NULL,
                                   sample_mean = NULL,
                                   sample_th = NULL,
                                   sample_nobs = NULL,
                                   slot_par_table = NULL,
                                   ngroups = NULL,
                                   dotdotdot = NULL,
                                   flat_model = NULL,
                                   model = NULL,
                                   nacov = NULL,
                                   wls_v = NULL) {
  # # # # # # # # # # #
  # #  3. lavdata  # #
  # # # # # # # # # # #

  # if slot_data not null
  #   copy slot_data to lavdata
  # else
  #   create lavdata via function lav_lavdata, setting ov.names to ov.names.y
  #     if lavoptions$conditional.x
  # if lavdata$data.type is "none"
  #   set lavoptions$do.fit to FALSE
  #   if flat_model$est not null set lavoptions$start to "est", else set
  #     it to "simple"
  #   set lavoptions$se and lavoptions$test to "none"
  # else
  #   if lavdata$data.type is "moment"
  #     if estimator one of MLM, MLMV, MLR, ULSM, ULSMV, ULSMVS and NACOV
  #       is NULL: *** error ***
  #     if estimator one of WLS, WLSM, WLSMV, WLSMVS, DWLS and wls_v is
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
      ov_names = tmp_ov_names,
      ov_names_x = ov_names_x,
      ov_names_l = ov_names_l,
      ov_names_aux = ov_names_aux,
      ordered = ordered,
      sampling_weights = sampling_weights,
      sample_cov = sample_cov,
      sample_mean = sample_mean,
      sample_th = sample_th,
      sample_nobs = sample_nobs,
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
          "estimator %s requires full data or user-provided wls_v and NACOV",
          estimator))
      }
    }
    # catch here some options that will not work with moments
    if (lavoptions$se == "bootstrap") {
      lav_msg_stop(gettext("bootstrapping requires full data"))
    }
    # robust (sandwich) standard errors need the raw-data Gamma/NACOV;
    # robust.sem/robust.sem.nt can be rescued by a user-provided NACOV=,
    # but robust.huber.white/first.order need the casewise scores (raw data)
    if (any(lavoptions$se == c("robust.sem", "robust.sem.nt")) &&
      is.null(nacov)) {
      lav_msg_stop(gettextf(
        "se = %s requires full data or a user-provided NACOV= argument.",
        dQuote(lavoptions$se)))
    } else if (any(lavoptions$se == c("robust.huber.white", "first.order"))) {
      lav_msg_stop(gettextf(
        "se = %s requires full data.", dQuote(lavoptions$se)))
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
