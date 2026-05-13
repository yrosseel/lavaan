lav_lavaan_step04_partable_cfa_fast <- function(slot_par_table = NULL, # nolint
                                                flat_model = NULL,
                                                lavoptions = NULL,
                                                lavdata = NULL,
                                                constraints = NULL) {
  if (!is.null(slot_par_table) ||
      !is.null(constraints) && any(nchar(constraints) > 0L) ||
      is.null(flat_model$lhs) ||
      is.null(flat_model$op) ||
      is.null(flat_model$rhs) ||
      is.null(flat_model$mod.idx) ||
      any(flat_model$op != "=~") ||
      any(flat_model$mod.idx != 0L) ||
      any(flat_model$lhs == flat_model$rhs) ||
      any(grepl(":", flat_model$lhs, fixed = TRUE)) ||
      any(grepl(":", flat_model$rhs, fixed = TRUE)) ||
      (!is.null(flat_model$block) && any(flat_model$block != 1L)) ||
      (!is.null(flat_model$level) && any(flat_model$level > 1L)) ||
      (!is.null(flat_model$group) && any(flat_model$group != 1L)) ||
      (!is.null(flat_model$rv) && any(nchar(flat_model$rv) > 0L)) ||
      (!is.null(flat_model$efa) && any(nchar(flat_model$efa) > 0L)) ||
      !isTRUE(lavoptions$std.lv) ||
      isTRUE(lavoptions$meanstructure) ||
      isTRUE(lavoptions$conditional.x) ||
      isTRUE(lavoptions$fixed.x) ||
      isTRUE(lavoptions$correlation) ||
      isTRUE(lavoptions$marker.int.zero) ||
      isTRUE(lavoptions$orthogonal) ||
      isTRUE(lavoptions$orthogonal.x) ||
      isTRUE(lavoptions$orthogonal.y) ||
      isTRUE(lavoptions$rotation.args$orthogonal) ||
      !isTRUE(lavoptions$auto.var) ||
      !isTRUE(lavoptions$auto.cov.lv.x) ||
      !isTRUE(lavoptions$auto.fix.single) ||
      isTRUE(lavoptions$auto.fix.first) ||
      length(lavoptions$group.equal) > 0L ||
      length(lavoptions$group.partial) > 0L ||
      isTRUE(lavoptions$group.w.free) ||
      lavdata@ngroups != 1L ||
      lavdata@nlevels != 1L ||
      length(lavdata@ov.names.x[[1L]]) > 0L ||
      length(lavdata@ordered) > 0L) {
    return(NULL)
  }

  lv_names <- unique(flat_model$lhs)
  ov_names <- unique(flat_model$rhs)
  if (length(lv_names) == 0L ||
      length(ov_names) == 0L ||
      anyNA(lv_names) ||
      anyNA(ov_names) ||
      any(!nzchar(lv_names)) ||
      any(!nzchar(ov_names)) ||
      any(flat_model$rhs %in% lv_names)) {
    return(NULL)
  }

  lv_cov <- if (length(lv_names) > 1L) {
    utils::combn(lv_names, 2L)
  } else {
    matrix(character(0L), nrow = 2L)
  }
  n_load <- length(flat_model$lhs)
  n_ov <- length(ov_names)
  n_lv <- length(lv_names)
  n_lvcov <- ncol(lv_cov)
  n_rows <- n_load + n_ov + n_lv + n_lvcov
  free_load <- seq_len(n_load)
  free_ov <- n_load + seq_len(n_ov)
  free_lvcov <- n_load + n_ov + seq_len(n_lvcov)

  lavpartable <- list(
    id = seq_len(n_rows),
    lhs = c(flat_model$lhs, ov_names, lv_names, lv_cov[1L, ]),
    op = c(
      rep("=~", n_load),
      rep("~~", n_ov + n_lv + n_lvcov)
    ),
    rhs = c(flat_model$rhs, ov_names, lv_names, lv_cov[2L, ]),
    user = c(rep(1L, n_load), rep(0L, n_ov + n_lv + n_lvcov)),
    block = rep(1L, n_rows),
    group = rep(1L, n_rows),
    free = c(free_load, free_ov, rep(0L, n_lv), free_lvcov),
    ustart = c(
      rep(as.numeric(NA), n_load + n_ov),
      rep(1, n_lv),
      rep(as.numeric(NA), n_lvcov)
    ),
    exo = rep(0L, n_rows),
    label = rep("", n_rows),
    plabel = paste0(".p", seq_len(n_rows), ".")
  )

  ov_list <- list(ov_names)
  lv_list <- list(lv_names)
  empty_names <- list(character(0L))
  empty_ov_idx <- list(integer(0L))
  empty_lv_idx <- list(integer(0L))
  lv_marker <- list(stats::setNames(rep("", n_lv), lv_names))
  pta <- list(
    vnames = list(
      ov = ov_list,
      ov.x = empty_names,
      ov.nox = ov_list,
      ov.model = ov_list,
      ov.y = empty_names,
      ov.num = ov_list,
      ov.ord = empty_names,
      ov.ind = ov_list,
      ov.cind = empty_names,
      ov.orphan = empty_names,
      ov.interaction = empty_names,
      ov.efa = empty_names,
      th = empty_names,
      th.mean = ov_list,
      lv = lv_list,
      lv.regular = lv_list,
      lv.formative = empty_names,
      lv.composite = empty_names,
      lv.x = lv_list,
      lv.y = empty_names,
      lv.nox = empty_names,
      lv.nonnormal = empty_names,
      lv.interaction = empty_names,
      lv.efa = empty_names,
      lv.rv = empty_names,
      lv.ind = empty_names,
      lv.ho = empty_names,
      lv.marker = lv_marker,
      eqs.y = empty_names,
      eqs.x = empty_names
    ),
    vidx = list(
      ov = list(seq_len(n_ov)),
      ov.x = empty_ov_idx,
      ov.nox = list(seq_len(n_ov)),
      ov.model = list(seq_len(n_ov)),
      ov.y = empty_ov_idx,
      ov.num = list(seq_len(n_ov)),
      ov.ord = empty_ov_idx,
      ov.ind = list(seq_len(n_ov)),
      ov.cind = empty_ov_idx,
      ov.orphan = empty_ov_idx,
      ov.interaction = empty_ov_idx,
      ov.efa = empty_ov_idx,
      th = empty_ov_idx,
      th.mean = list(seq_len(n_ov)),
      lv = list(seq_len(n_lv)),
      lv.regular = list(seq_len(n_lv)),
      lv.formative = empty_lv_idx,
      lv.composite = empty_lv_idx,
      lv.x = list(seq_len(n_lv)),
      lv.y = empty_lv_idx,
      lv.nox = empty_lv_idx,
      lv.nonnormal = empty_lv_idx,
      lv.interaction = empty_lv_idx,
      lv.efa = empty_lv_idx,
      lv.rv = empty_lv_idx,
      lv.ind = empty_lv_idx,
      lv.ho = empty_lv_idx,
      lv.marker = list(rep(NA_integer_, n_lv)),
      eqs.y = empty_lv_idx,
      eqs.x = empty_lv_idx
    ),
    meanstructure = FALSE,
    nblocks = 1L,
    ngroups = 1L,
    nlevels = 1L,
    nvar = list(n_ov),
    nfac = list(n_lv),
    nfac.nonnormal = list(0L),
    th.idx = list(numeric(n_ov))
  )

  lav_partable_set_cache(lavpartable, pta)
}

lav_lavaan_step04_partable <- function(slot_par_table = NULL, # nolint
                                       model = NULL,
                                       flat_model = NULL,
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
  #     create lavpartable via function lavParTable (=lav_model_partable)
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

  lavpartable <- lav_lavaan_step04_partable_cfa_fast(
    slot_par_table = slot_par_table,
    flat_model = flat_model,
    lavoptions = lavoptions,
    lavdata = lavdata,
    constraints = constraints
  )
  if (!is.null(lavpartable)) {
    return(list(
      lavoptions = lavoptions,
      lavpartable = lavpartable
    ))
  }

  if (!is.null(slot_par_table)) {
    lavpartable <- lav_partable_set_cache(slot_par_table)
  } else if (is.character(model) ||
    inherits(model, "formula") ||
  # model was already a flat.model
  (is.list(model) && !is.null(model$mod.idx) &&
   !is.null(attr(model, "modifiers")))) {
    if (lav_verbose()) {
      cat("lavpartable        ...")
    }
    # check flat.model before we proceed
    if (lav_debug()) {
      print(as.data.frame(flat_model))
    }
    # catch ~~ of fixed.x covariates if fixed.x = TRUE
    # --> done inside lav_model_partable!

    # if(lavoptions$fixed.x) {
    #    tmp <- lav_partable_vnames(flat.model, type = "ov.x",
    #                               ov.x.fatal = FALSE, warn = TRUE)
    # tmp <- try(lav_partable_vnames(flat.model, type = "ov.x",
    #                                         ov.x.fatal = TRUE),
    #           silent = TRUE)
    # if(inherits(tmp, "try-error")) {
    #    warning("lavaan WARNING: syntax contains parameters involving ",
    #      "exogenous covariates; switching to fixed.x = FALSE")
    #    lavoptions$fixed.x <- FALSE
    # }
    # }
    # if(lavoptions$conditional.x) {
    #    tmp <- lav_partable_vnames(flat.model,
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
        varTable = tmp_data_ov,
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
        composites = lavoptions$composites,
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
        as.data.frame. = FALSE
      )
    lavpartable <- lav_partable_set_cache(lavpartable)
    if (lav_verbose()) {
      cat(" done.\n")
    }
  } else if (inherits(model, "lavaan")) {
    lavpartable <- lav_partable_set_cache(as.list(parTable(model)), model@pta)
  } else if (is.list(model)) {
    # we already checked this when creating flat.model
    # but we may need to complete it
    lavpartable <- as.list(flat_model) # in case model is a data.frame
    # complete table
    lavpartable <- as.list(lav_partable_complete(lavpartable))
    lavpartable <- lav_partable_set_cache(lavpartable)
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
    junk <- lav_partable_check(lavpartable,
      categorical = lavoptions$.categorical
    )
    rm(junk)
  }

  # for EM only (for now), force fixed-to-zero (residual) variances
  # to be slightly larger than zero
  if (lavoptions$optim.method == "em") {
    zero_var_idx <- which(lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs &
      lavpartable$free == 0L &
      lavpartable$ustart == 0)
    if (length(zero_var_idx) > 0L) {
      lavpartable$ustart[zero_var_idx] <- lavoptions$em.zerovar.offset
    }
    lavpartable <- lav_partable_set_cache(lavpartable, NULL, force = TRUE)
  }

  list(
    lavoptions  = lavoptions,
    lavpartable = lavpartable
  )
}
