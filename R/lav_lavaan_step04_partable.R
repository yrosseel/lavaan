lav_lavaan_step04_partable <- function(slotParTable = NULL, # nolint
                                       model = NULL,
                                       flat.model = NULL,
                                       lavoptions = NULL,
                                       lavdata = NULL,
                                       constraints = NULL) {
  # # # # # # # # # # # #
  # #  4. lavpartable # #
  # # # # # # # # # # # #

  # if slotParTable not null
  #   copy slotParTable to lavpartable
  # else
  #   if model is character or formula
  #     create a temporary variable tmp.data.ov equal to lavdata@ov
  #     if estimator "catML"
  #       set meanstructure to FALSE
  #       set the member type in the temporary variable tmp.data.ov to a
  #         numeric vector with all zeroes
  #     create lavpartable via function lavParTable (=lavaanify)
  #                     using the temporary variable for parameter varTable
  #   else
  #     if model is lavaan object
  #       set lavpartable = parTable(model)
  #     else
  #       if model is a list
  #         set lavpartable to
  #           as.list(lav_partable_complete(as.list(flat.model)))
  #       else
  #         *** error ***
  # if slotParTable is NULL check lavpartable via lav_partable_check
  # if lavoptions$optim.method is "em" and there are variances specified in
  #   partable with free = 0L and
  #    starting value ustart 0, set ustart for these variances to
  #    lavoptions$em.zerovar.offset

  if (!is.null(slotParTable)) {
    lavpartable <- lav_partable_set_cache(slotParTable)
  } else if (is.character(model) ||
    inherits(model, "formula") ||
	# model was already a flat.model
	(is.list(model) && !is.null(model$mod.idx) &&
	 !is.null(attr(model, "modifiers")))) {
    if (lavoptions$verbose) {
      cat("lavpartable        ...")
    }
    # check flat.model before we proceed
    if (lavoptions$debug) {
      print(as.data.frame(flat.model))
    }
    # catch ~~ of fixed.x covariates if fixed.x = TRUE
    # --> done inside lavaanify!

    # if(lavoptions$fixed.x) {
    #    tmp <- lav_partable_vnames(flat.model, type = "ov.x",
    #                               ov.x.fatal = FALSE, warn = TRUE)
    # tmp <- try(vnames(flat.model, type = "ov.x", ov.x.fatal = TRUE),
    #           silent = TRUE)
    # if(inherits(tmp, "try-error")) {
    #    warning("lavaan WARNING: syntax contains parameters involving ",
    #      "exogenous covariates; switching to fixed.x = FALSE")
    #    lavoptions$fixed.x <- FALSE
    # }
    # }
    # if(lavoptions$conditional.x) {
    #    tmp <- vnames(flat.model, type = "ov.x", ov.x.fatal = TRUE)
    # }
    tmp.data.ov <- lavdata@ov
    if (lavoptions$estimator == "catML") {
      lavoptions$meanstructure <- FALSE
      tmp.data.ov$type <- rep("numeric", length(tmp.data.ov$type))
    }
    lavpartable <-
      lavParTable(
        model = flat.model,
        constraints = constraints,
        varTable = tmp.data.ov,
        ngroups = lavdata@ngroups,
        meanstructure = lavoptions$meanstructure,
        int.ov.free = lavoptions$int.ov.free,
        int.lv.free = lavoptions$int.lv.free,
        marker.int.zero = lavoptions$marker.int.zero,
        orthogonal = lavoptions$orthogonal,
        orthogonal.x = lavoptions$orthogonal.x,
        orthogonal.y = lavoptions$orthogonal.y,
        orthogonal.efa = lavoptions$rotation.args$orthogonal,
        conditional.x = lavoptions$conditional.x,
        fixed.x = lavoptions$fixed.x,
        std.lv = lavoptions$std.lv,
        correlation = lavoptions$correlation,
        effect.coding = lavoptions$effect.coding,
        ceq.simple = lavoptions$ceq.simple,
        parameterization = lavoptions$parameterization,
        auto.fix.first = lavoptions$auto.fix.first,
        auto.fix.single = lavoptions$auto.fix.single,
        auto.var = lavoptions$auto.var,
        auto.cov.lv.x = lavoptions$auto.cov.lv.x,
        auto.cov.y = lavoptions$auto.cov.y,
        auto.th = lavoptions$auto.th,
        auto.delta = lavoptions$auto.delta,
        auto.efa = lavoptions$auto.efa,
        group.equal = lavoptions$group.equal,
        group.partial = lavoptions$group.partial,
        group.w.free = lavoptions$group.w.free,
        debug = lavoptions$debug,
        warn = lavoptions$warn,
        as.data.frame. = FALSE
      )
    lavpartable <- lav_partable_set_cache(lavpartable)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  } else if (inherits(model, "lavaan")) {
    lavpartable <- lav_partable_set_cache(as.list(parTable(model)), model@pta)
  } else if (is.list(model)) {
    # we already checked this when creating flat.model
    # but we may need to complete it
    lavpartable <- as.list(flat.model) # in case model is a data.frame
    # complete table
    lavpartable <- as.list(lav_partable_complete(lavpartable))
    lavpartable <- lav_partable_set_cache(lavpartable)
  } else {
    lav_msg_stop(gettextf(
      "model [type = %s] is not of type character or list", class(model)))
    }
  if (lavoptions$debug) {
    print(as.data.frame(lavpartable))
  }

  # at this point, we should check if the partable is complete
  # or not; this is especially relevant if the lavaan() function
  # was used, but the user has forgotten some variances/intercepts...
  if (is.null(slotParTable)) {
    junk <- lav_partable_check(lavpartable,
      categorical = lavoptions$.categorical,
      warn = TRUE
    )
    rm(junk)
  }

  # for EM only (for now), force fixed-to-zero (residual) variances
  # to be slightly larger than zero
  if (lavoptions$optim.method == "em") {
    zero.var.idx <- which(lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs &
      lavpartable$free == 0L &
      lavpartable$ustart == 0)
    if (length(zero.var.idx) > 0L) {
      lavpartable$ustart[zero.var.idx] <- lavoptions$em.zerovar.offset
    }
    lavpartable <- lav_partable_set_cache(lavpartable, NULL, force = TRUE)
  }

  list(
    lavoptions  = lavoptions,
    lavpartable = lavpartable
  )
}
