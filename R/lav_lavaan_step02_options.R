lav_lavaan_step02_options <- function(slot_options = NULL, # nolint
                                      slot_data = NULL, # nolint
                                      flat_model = NULL,
                                      ordered = NULL,
                                      ordered_orig = NULL,
                                      sample_cov = NULL,
                                      sample_mean = NULL,
                                      sample_th = NULL,
                                      sample_nobs = NULL,
                                      ov_names_l = NULL,
                                      sampling_weights = NULL,
                                      constraints = NULL,
                                      group = NULL,
                                      ov_names_x = NULL,
                                      ov_names_y = NULL,
                                      dotdotdot = NULL,
                                      cluster = NULL,
                                      data = NULL) {
  # # # # # # # # # # # #
  # #  2. lavoptions  # #
  # # # # # # # # # # # #

  # if slotOptions not NULL
  #   copy to lavoptions and modify categorical/clustered/multilevel
  #     inserting a "." in the first position
  #   if necessary, overwrite with values in dotdotdot and issue a warning
  #   check if all names in dotdotdot are possible options, if not *** error ***
  #   create complete option list (lav_options_default) and substitute values
  #     given in dotdotdot
  #   if data, slotData and sample.cov NULL: opt$bounds = FALSE
  #   if slotData$data.type != "full" or (slotData and data = NULL):
  #     opt$missing = "listwise"
  #   set categorical mode ON if
  #     - an operator "|" (threshold) was used
  #     - data not NULL and one or more elements in ordered parameter
  #     - sample.th provided
  #     - at least one of the non-exogenous observed variables is "ordered"
  #       (ordered factor in R)
  #   if opt$estimator == "catml": set categorical mode OFF
  #       TODO: estimator = "CATML" isn't mentioned in lavOptions / estimator
  #             help text !?
  #   if cluster not NULL, set opt$.clustered TRUE and  *** error *** if
  #     categorical mode is ON
  #   opt$.multilevel = (length(ov.names.l) > 0L &&
  #                      length(ov.names.l[[1]]) > 1L)
  #   if sampling.weights not NULL en categorical mode OFF and opt$estimator
  #     in ("default", "ML", "PML")
  #     set opt$estimator to "MLR"
  #   if constraints present and estimator == "ML", set opt$information to
  #     c("observed", "observed")
  #   if there is an operator "~1" in flat.model and sample.mean not NULL,
  #     set opt$meanstructure TRUE
  #   if there are no exogenous variables but conditional.x explicitly
  #     requested: ** warning **
  #   if there are no exogenous variables set opt$conditional.x FALSE
  #   if there are no exogenous variables and fixed.x not explicitly requested,
  #     set opt$fixed.x to FALSE
  #   if allow.empty.cell and estimator not Bayes, issue a warning

  if (!is.null(slot_options)) {
    lavoptions <- slot_options

    # backwards compatibility
    if (!is.null(lavoptions$categorical)) {
      lavoptions$.categorical <- lavoptions$categorical
      lavoptions$categorical <- NULL
    }
    if (!is.null(lavoptions$clustered)) {
      lavoptions$.clustered <- lavoptions$clustered
      lavoptions$clustered <- NULL
    }
    if (!is.null(lavoptions$multilevel)) {
      lavoptions$.multilevel <- lavoptions$multilevel
      lavoptions$multilevel <- NULL
    }

    # but what if other 'options' are given anyway (eg 'start = ')?
    # give a warning!
    if (length(dotdotdot) > 0L) {
      dot_names <- names(dotdotdot)
      op_idx <- which(dot_names %in% names(slot_options))
      lav_msg_warn(gettext(
        "the following argument(s) override(s) the options in slotOptions:"),
        paste(dot_names[op_idx], collapse = " ")
      )
      lavoptions[dot_names[op_idx]] <- dotdotdot[op_idx]
    }
  } else {
    if (lav_verbose()) {
      cat("lavoptions         ...")
    }

    # load default options
    opt <- lav_options_default()

    # catch unknown options
    ok_names <- names(opt)
    dot_names <- names(dotdotdot)
    wrong_idx <- which(!dot_names %in% ok_names)
    if (length(wrong_idx) > 0L) {
      # stop or warning?? stop for now (there could be more)
      lav_msg_stop(ngettext(length(wrong_idx),
        "unknown argument:", "unknown arguments:"),
        lav_msg_view(dot_names[wrong_idx], "none", FALSE)
      )
    }

    # modifyList
    opt <- modifyList(opt, dotdotdot)

    # extract estimator
    if (is.list(opt$estimator)) {
      estimator <- opt$estimator$estimator
    } else {
      estimator <- opt$estimator
    }

    # no data?
    if (is.null(slot_data) && is.null(data) && is.null(sample_cov)) {
      opt$bounds <- FALSE
    }

    # only sample moments?
    if (!is.null(slot_data) && !slot_data@data.type == "full") {
      opt$missing <- "listwise"
    } else if (is.null(slot_data) && is.null(data)) {
      opt$missing <- "listwise"
    }

    # categorical mode?
    opt$.categorical <- FALSE
    if (any(flat_model$op == "|")) {
      opt$.categorical <- TRUE
    } else if (!is.null(data) && length(ordered) > 0L) {
      opt$.categorical <- TRUE
    } else if (!is.null(sample_th)) {
      opt$.categorical <- TRUE
    } else if (is.data.frame(data)) {
      # first check if we can find ov.names.y in Data
      tmp_ov_names_y <- unique(unlist(ov_names_y))
      # remove possible interaction terms involving an y term
      int_idx <- which(grepl(":", tmp_ov_names_y))
      if (length(int_idx) > 0L) {
        tmp_ov_names_y <- tmp_ov_names_y[-int_idx]
      }
      idx_missing <- which(!(tmp_ov_names_y %in% names(data)))
      if (length(idx_missing)) {
        lav_msg_stop(
          gettext("missing observed variables in dataset:"),
          paste(tmp_ov_names_y[idx_missing], collapse = " ")
        )
      }
      if (any(sapply(data[, tmp_ov_names_y], inherits, "ordered"))) {
        opt$.categorical <- TRUE
      }
    }
    if (tolower(estimator) == "catml") {
      opt$.categorical <- FALSE
    }

    # check for WLSMV estimator with NULL ordered= argument (new in 0.6-22)
    # to avoid unintentional use of WLSMV
    # we do this here, because we need the value of opt$.categorical
    if (!opt$.categorical && is.null(ordered_orig) && !estimator == "default") {
      estimator <- toupper(estimator)
      if (estimator %in% c("WLSM", "WLSMV", "DWLS")) {
        lav_msg_stop(gettextf(
        "Estimator %s is typically used with categorical data only. To use it
         with continuous data, please explicitly set the ordered= argument
         to FALSE.",
        dQuote(lav_options_estimatorgroup(estimator))))
      }
    }

    # clustered?
    if (length(cluster) > 0L) {
      opt$.clustered <- TRUE
      if (opt$.categorical && toupper(estimator) != "PML") {
        lav_msg_stop(gettext("categorical + clustered is not supported yet."))
      }
    } else {
      opt$.clustered <- FALSE
    }

    # multilevel?
    if (length(ov_names_l) > 0L && length(ov_names_l[[1]]) > 1L) {
      opt$.multilevel <- TRUE
    } else {
      opt$.multilevel <- FALSE
    }

    # sampling weights? force MLR
    # HJ 18/10/23: Except for PML
    if (!is.null(sampling_weights) && !opt$.categorical &&
      toupper(estimator) %in% c("DEFAULT", "ML", "PML")) {
      if (opt$se != "none") {
        opt$se <- "robust.huber.white"
      }
      if (opt$se != "none") {
        opt$test <- "yuan.bentler.mplus"
      }
    }

    # constraints
    if (any(nchar(constraints) > 0L) && toupper(estimator) %in% c("ML")) {
      opt$information <- c("observed", "observed")
    }

    # meanstructure
    if (any(flat_model$op == "~1") || !is.null(sample_mean)) {
      opt$meanstructure <- TRUE
    }
    if (!is.null(group) && is.null(dotdotdot$meanstructure)) {
      opt$meanstructure <- TRUE
    }

    # conditional.x
    if ((is.list(ov_names_x) &&
      sum(sapply(ov_names_x, FUN = length)) == 0L) ||
      (is.character(ov_names_x) && length(ov_names_x) == 0L)) {
      # if explicitly set to TRUE, give warning
      if (is.logical(dotdotdot$conditional.x) && dotdotdot$conditional.x) {
        lav_msg_warn(gettext(
          "no exogenous covariates; conditional.x will be set to FALSE"))
      }
      opt$conditional.x <- FALSE
    }

    # fixed.x
    if ((is.list(ov_names_x) &&
      sum(sapply(ov_names_x, FUN = length)) == 0L) ||
      (is.character(ov_names_x) && length(ov_names_x) == 0L)) {
      # if explicitly set to TRUE, give warning
      if (is.logical(dotdotdot$fixed.x) && dotdotdot$fixed.x) {
        # ok, we respect this: keep fixed.x = TRUE
      } else {
        opt$fixed.x <- FALSE
      }
    }

    # allow.empty.cell
    if (opt$allow.empty.cell && opt$do.fit && toupper(estimator) != "BAYES") {
      lav_msg_warn(gettext(
        "allow.empty.cell is not intended to salvage estimation of this model,",
        " see ?lavOptions"))
    }

    # fill in remaining "default" values
    lavoptions <- lav_options_set(opt)

    # store check.sigma.pd in lavaan_cache_env
    assign("opt.check.sigma.pd", opt$check.sigma.pd, lavaan_cache_env)        # nolint

    if (lav_verbose()) {
      cat(" done.\n")
    }
  }

  lavoptions
}
