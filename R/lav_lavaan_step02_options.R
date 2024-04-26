lav_lavaan_step02_options <- function(slotOptions = NULL, # nolint
                                      slotData = NULL, # nolint
                                      flat.model = NULL,
                                      ordered = NULL,
                                      sample.cov = NULL,
                                      sample.mean = NULL,
                                      sample.th = NULL,
                                      sample.nobs = NULL,
                                      ov.names.l = NULL,
                                      sampling.weights = NULL,
                                      constraints = NULL,
                                      group = NULL,
                                      ov.names.x = NULL,
                                      ov.names.y = NULL,
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
  #   if there are no exogene variables but conditional.x explicitly
  #     requested: ** warning **
  #   if there are no exogene variables set opt$conditional.x FALSE
  #   if there are no exogene variables and fixed.x not explicitly requested,
  #     set opt$fixed.x to FALSE

  if (!is.null(slotOptions)) {
    lavoptions <- slotOptions

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
      dot.names <- names(dotdotdot)
      op.idx <- which(dot.names %in% names(slotOptions))
      lav_msg_warn(gettext(
        "the following argument(s) override(s) the options in slotOptions:"),
        paste(dot.names[op.idx], collapse = " ")
      )
      lavoptions[dot.names[op.idx]] <- dotdotdot[op.idx]
    }
  } else {
    if (!is.null(dotdotdot$verbose) && dotdotdot$verbose) {
      cat("lavoptions         ...")
    }

    # load default options
    opt <- lav_options_default()

    # catch unknown options
    ok.names <- names(opt)
    dot.names <- names(dotdotdot)
    wrong.idx <- which(!dot.names %in% ok.names)
    if (length(wrong.idx) > 0L) {
      # stop or warning?? stop for now (there could be more)
      lav_msg_stop(ngettext(length(wrong.idx),
        "unknown argument:", "unknown arguments:"),
        lav_msg_view(dot.names[wrong.idx], "none", FALSE)
      )
    }

    # modifyList
    opt <- modifyList(opt, dotdotdot)

    # no data?
    if (is.null(slotData) && is.null(data) && is.null(sample.cov)) {
      opt$bounds <- FALSE
    }

    # only sample moments?
    if (!is.null(slotData) && !slotData@data.type == "full") {
      opt$missing <- "listwise"
    } else if (is.null(slotData) && is.null(data)) {
      opt$missing <- "listwise"
    }

    # categorical mode?
    opt$.categorical <- FALSE
    if (any(flat.model$op == "|")) {
      opt$.categorical <- TRUE
    } else if (!is.null(data) && length(ordered) > 0L) {
      opt$.categorical <- TRUE
    } else if (!is.null(sample.th)) {
      opt$.categorical <- TRUE
    } else if (is.data.frame(data)) {
      # first check if we can find ov.names.y in Data
      tmp.ov.names.y <- unique(unlist(ov.names.y))
      # remove possible interaction terms involving an y term
      int.idx <- which(grepl(":", tmp.ov.names.y))
      if (length(int.idx) > 0L) {
        tmp.ov.names.y <- tmp.ov.names.y[-int.idx]
      }
      idx.missing <- which(!(tmp.ov.names.y %in% names(data)))
      if (length(idx.missing)) {
        lav_msg_stop(
          gettext("missing observed variables in dataset:"),
          paste(tmp.ov.names.y[idx.missing], collapse = " ")
        )
      }
      if (any(sapply(data[, tmp.ov.names.y], inherits, "ordered"))) {
        opt$.categorical <- TRUE
      }
    }
    if (tolower(opt$estimator) == "catml") {
      opt$.categorical <- FALSE
    }

    # clustered?
    if (length(cluster) > 0L) {
      opt$.clustered <- TRUE
      if (opt$.categorical & opt$estimator != "PML") {
        lav_msg_stop(gettext("categorical + clustered is not supported yet."))
      }
    } else {
      opt$.clustered <- FALSE
    }

    # multilevel?
    if (length(ov.names.l) > 0L && length(ov.names.l[[1]]) > 1L) {
      opt$.multilevel <- TRUE
    } else {
      opt$.multilevel <- FALSE
    }

    # sampling weights? force MLR
    # HJ 18/10/23: Except for PML
    if (!is.null(sampling.weights) && !opt$.categorical &&
      opt$estimator %in% c("default", "ML", "PML")) {
      opt$estimator <- "MLR"
    }

    # constraints
    if (any(nchar(constraints) > 0L) && opt$estimator %in% c("ML")) {
      opt$information <- c("observed", "observed")
    }

    # meanstructure
    if (any(flat.model$op == "~1") || !is.null(sample.mean)) {
      opt$meanstructure <- TRUE
    }
    if (!is.null(group) && is.null(dotdotdot$meanstructure)) {
      opt$meanstructure <- TRUE
    }

    # conditional.x
    if ((is.list(ov.names.x) &&
      sum(sapply(ov.names.x, FUN = length)) == 0L) ||
      (is.character(ov.names.x) && length(ov.names.x) == 0L)) {
      # if explicitly set to TRUE, give warning
      if (is.logical(dotdotdot$conditional.x) && dotdotdot$conditional.x) {
        lav_msg_warn(
          gettext("no exogenous covariates; conditional.x will be set to FALSE"))
      }
      opt$conditional.x <- FALSE
    }

    # fixed.x
    if ((is.list(ov.names.x) &&
      sum(sapply(ov.names.x, FUN = length)) == 0L) ||
      (is.character(ov.names.x) && length(ov.names.x) == 0L)) {
      # if explicitly set to TRUE, give warning
      if (is.logical(dotdotdot$fixed.x) && dotdotdot$fixed.x) {
        # ok, we respect this: keep fixed.x = TRUE
      } else {
        opt$fixed.x <- FALSE
      }
    }

    # fill in remaining "default" values
    lavoptions <- lav_options_set(opt)

    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }

  lavoptions
}
