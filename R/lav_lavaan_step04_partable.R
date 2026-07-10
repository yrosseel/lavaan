lav_step04_pt <- function(slot_par_table = NULL,
                                       model = NULL,
                                       flat_model = NULL,
                                       lavoptions = NULL,
                                       lavdata = NULL,
                                       constraints = NULL,
                                       marker = NULL) {
  # # # # # # # # # # # #
  # #  4. lavpartable # #
  # # # # # # # # # # # #

  # if slot_par_table not null
  #   copy slot_par_table to lavpartable
  # else
  #   if model is character or formula
  #     create a temporary variable tmp.data.ov equal to lavdata@ov
  #     if estimator "catML"
  #       set meanstructure to FALSE
  #       set the member type in the temporary variable tmp.data.ov to a
  #         numeric vector with all zeroes
  #     create lavpartable via function lavParTable (=lav_model_pt)
  #                     using the temporary variable for parameter varTable
  #   else
  #     if model is lavaan object
  #       set lavpartable = parTable(model)
  #     else
  #       if model is a list
  #         set lavpartable to
  #           as.list(lav_pt_complete(as.list(flat_model)))
  #       else
  #         *** error ***
  # if slot_par_table is NULL check lavpartable via lav_partable_check
  # if lavoptions$optim.method is "em" and there are variances specified in
  #   partable with free = 0L and
  #    starting value ustart 0, set ustart for these variances to
  #    lavoptions$em.args$zerovar_offset

  if (!is.null(slot_par_table)) {
    lavpartable <- lav_pt_set_cache(slot_par_table)
  } else if (is.character(model) ||
    inherits(model, "formula") ||
  # model was already a flat_model
  (is.list(model) && !is.null(model$mod.idx) &&
   !is.null(attr(model, "modifiers")))) {
    if (lav_verbose()) {
      cat("lavpartable        ...")
    }
    # check flat_model before we proceed
    if (lav_debug()) {
      print(as.data.frame(flat_model))
    }
    # catch ~~ of fixed.x covariates if fixed.x = TRUE
    # --> done inside lav_model_pt!

    # if(lavoptions$fixed.x) {
    #    tmp <- lav_pt_vnames(flat_model, type = "ov.x",
    #                               ov.x.fatal = FALSE, warn = TRUE)
    # tmp <- try(lav_pt_vnames(flat_model, type = "ov.x",
    #                                         ov.x.fatal = TRUE),
    #           silent = TRUE)
    # if(inherits(tmp, "try-error")) {
    #    warning("lavaan WARNING: syntax contains parameters involving ",
    #      "exogenous covariates; switching to fixed.x = FALSE")
    #    lavoptions$fixed.x <- FALSE
    # }
    # }
    # if(lavoptions$conditional.x) {
    #    tmp <- lav_pt_vnames(flat_model,
    #                  type = "ov.x", ov.x.fatal = TRUE)
    # }
    tmp_data_ov <- lavdata@ov
    if (lavoptions$estimator == "catML") {
      lavoptions$meanstructure <- FALSE
      tmp_data_ov$type <- rep("numeric", length(tmp_data_ov$type))
    }
    lavpartable <-
      lavParTable(
        model = flat_model,
        constraints = constraints,
        var_table = tmp_data_ov,
        ngroups = lavdata@ngroups,
        meanstructure = lavoptions$meanstructure,
        int_ov_free = lavoptions$int.ov.free,
        int_lv_free = lavoptions$int.lv.free,
        marker_int_zero = lavoptions$marker.int.zero,
        orthogonal = lavoptions$orthogonal,
        orthogonal_x = lavoptions$orthogonal.x,
        orthogonal_y = lavoptions$orthogonal.y,
        orthogonal_efa = lavoptions$rotation.args$orthogonal,
        conditional_x = lavoptions$conditional.x,
        fixed_x = lavoptions$fixed.x,
        std_lv = lavoptions$std.lv,
        correlation = if (length(lavoptions$.correlation.ov) > 0L) {
          lavoptions$.correlation.ov
        } else {
          lavoptions$correlation
        },
        composites = lavoptions$composites,
        composites_cov_free = identical(lavoptions$composites.cov, "free"),
        effect_coding = lavoptions$effect.coding,
        ceq_simple = lavoptions$ceq.simple,
        parameterization = lavoptions$parameterization,
        auto_fix_first = lavoptions$auto.fix.first,
        marker = marker,
        auto_fix_single = lavoptions$auto.fix.single,
        auto_var = lavoptions$auto.var,
        auto_cov_lv_x = lavoptions$auto.cov.lv.x,
        auto_cov_x = lavoptions$auto.cov.x,
        auto_cov_y = lavoptions$auto.cov.y,
        auto_th = lavoptions$auto.th,
        auto_delta = lavoptions$auto.delta,
        auto_efa = lavoptions$auto.efa,
        group_equal = lavoptions$group.equal,
        group_partial = lavoptions$group.partial,
        group_w_free = lavoptions$group.w.free,
        as_data_frame = FALSE
      )
    lavpartable <- lav_pt_set_cache(lavpartable)
    if (lav_verbose()) {
      cat(" done.\n")
    }
  } else if (inherits(model, "lavaan")) {
    lavpartable <- lav_pt_set_cache(as.list(parTable(model)), model@pta)
  } else if (is.list(model)) {
    # we already checked this when creating flat_model
    # but we may need to complete it
    lavpartable <- as.list(flat_model) # in case model is a data.frame
    # complete table
    lavpartable <- as.list(lav_pt_complete(lavpartable))
    lavpartable <- lav_pt_set_cache(lavpartable)
  } else {
    lav_msg_stop(gettextf(
      "model [type = %s] is not of type character or list", class(model)))
    }
  if (lav_debug()) {
    print(as.data.frame(lavpartable))
  }

  # at this point, we should check if the partable is complete
  # or not; this is especially relevant if the lavaan() function
  # was used, but the user has forgotten some variances/intercepts...
  if (is.null(slot_par_table)) {
    junk <- lav_pt_check(lavpartable,
      categorical = lavoptions$.categorical
    )
    rm(junk)
  }

  # if optim.method = "em" was set by default (two-level + missing =
  # "ml"; see lav_options_set), quietly fall back to "nlminb" when the
  # EM optimizer does not support the data configuration (single-group
  # two-level only, no conditional.x); this must happen *before* the
  # zero-variance offset tweak below
  if (lavoptions$optim.method == "em" &&
      isTRUE(lavoptions$.optim.em.fallback) &&
      (lavdata@nlevels == 1L ||
       lavdata@ngroups > 1L ||
       lavoptions$conditional.x)) {
    lavoptions$optim.method <- "nlminb"
  }

  # for EM only (for now), force fixed-to-zero (residual) variances
  # to be slightly larger than zero
  if (lavoptions$optim.method == "em") {
    zero_var_idx <- which(lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs &
      lavpartable$free == 0L &
      lavpartable$ustart == 0)
    if (length(zero_var_idx) > 0L) {
      lavpartable$ustart[zero_var_idx] <- lavoptions$em.args$zerovar_offset
    }
    lavpartable <- lav_pt_set_cache(lavpartable, NULL, force = TRUE)
  }

  list(
    lavoptions  = lavoptions,
    lavpartable = lavpartable
  )
}
