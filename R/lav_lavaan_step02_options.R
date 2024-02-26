lav_lavaan_step02_options <- function(slotOptions, slotData, flat.model, ordered, sample.cov,
                                      sample.mean, sample.th, sample.nobs,  ov.names.l,
                                      sampling.weights, constraints, group, ov.names.x,
                                      ov.names.y, dotdotdot, cluster, data) {
  # # # # # # # # # # # #
  # #  2. lavoptions # # 
  # # # # # # # # # # # #
  # als slotOptions niet NULL
  #   kopieer naar lavoptions en voorzie eventueel categorical/clustered/multilevel van "." vooraan
  #   overchrijf eventueel met waarden in dotdotdot en geef waarschuwing
  #   kijk of namen in dotdotdot mogelijke opties zijn, zoniet *** error ***
  #   maak complete lijst (lav_options_default) en vervang waarden die meegegeven zijn in dotdotdot
  #   als data, slotData en sample.cov NULL zijn : opt$bounds = FALSE
  #   als slotData$data.type != "full" of slotData en data = NULL: opt$missing = "listwise"
  #   categorical mode als
  #     - er een operator "|" (threshold) gebruikt wordt
  #     - data niet NULL en element(en) in ordered parameter
  #     - sample.th opgegeven
  #     - minstens één van de non-exogenous observed variables is "ordered" (ordered factor in R)
  #   als opt$estimator == "catml" wordt categorical mode off gezet 
  #       TODO: estimator = "CATML" wordt niet vermeld in lavOptions / estimator !?
  #   als cluster opgegeven, zet opt$.clustered TRUE en geef *** error *** als categorical mode aan
  #   opt$.multilevel = (length(ov.names.l) > 0L && length(ov.names.l[[1]]) > 1L)
  #   
  #   
  #     
  #   TODO: hier verder doen
  #   
  #   
  #   
  #   
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
      warning("lavaan WARNING: the following argument(s) override(s) the options in slotOptions:\n\t\t",
              paste(dot.names[op.idx], collapse = " "))
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
      idx <- wrong.idx[1L] # only show first one TODO: why not show all of them?
      # stop or warning?? stop for now (there could be more)
      stop("lavaan ERROR: unknown argument `", dot.names[idx], "'")
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
      OV.names.y <- unique(unlist(ov.names.y))
      # remove possible interaction terms involving an y term
      int.idx <- which(grepl(":", OV.names.y))
      if (length(int.idx) > 0L) {
        OV.names.y <- OV.names.y[-int.idx]
      }
      idx.missing <- which(!(OV.names.y %in% names(data)))
      if (length(idx.missing)) {
        stop("lavaan ERROR: missing observed variables in dataset: ",
             paste(OV.names.y[idx.missing], collapse = " "))
      }
      if (any(sapply(data[, OV.names.y], inherits, "ordered"))) {
        opt$.categorical <- TRUE
      }
    }
    if (tolower(opt$estimator) == "catml") {
      opt$.categorical <- FALSE
    }
    
    # clustered?
    if (length(cluster) > 0L) {
      opt$.clustered <- TRUE
      if (opt$.categorical) {
        stop("lavaan ERROR: categorical + clustered is not supported yet.")
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
        warning("lavaan WARNING: no exogenous covariates; conditional.x will be set to FALSE")
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
  return(lavoptions)
}