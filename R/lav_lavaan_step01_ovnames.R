lav_lavaan_step01_ovnames_initflat <- function(slotParTable     = NULL, # nolint
                                               model            = NULL,
                                               dotdotdot.parser = "new") {
  # if slotPartable not NULL copy to flat.model
  # else
  #   if model is of type character
  #     parse model to flat.model
  #   else if model is a formula (** warning **)
  #     transform "~ x1 + x2 + x3" to character "f =~ x1 + x2 + x3",
  #         and parse to flat.model
  #     transform "y =~ x1 + x2 + x3" to character and parse to flat.model
  #     something else : *** error ***
  #   else if model is a lavaan object
  #     extract flat.model from @parTable
  #   else if model is a list
  #     if bare minimum present (columns lhs, op, rhs, free)
  #        flat.model = model
  #        replace column block by column group (if present)
  #     else
  #        --> *** error ***
  #   else
  #     --> ***error***
  #

  # 1a. get ov.names and ov.names.x (per group) -- needed for lavData()
  if (!is.null(slotParTable)) {
    flat.model <- slotParTable
  } else if (is.character(model)) {
    if (is.null(dotdotdot.parser)) {
      flat.model <- lavParseModelString(model, parser = "new")
    } else {
      flat.model <- lavParseModelString(model, parser = dotdotdot.parser)
    }
  } else if (inherits(model, "formula")) {
    # two typical cases:
    # 1. regression type formula
    # 2. no quotes, e.g. f =~ x1 + x2 + x3
    #    TODO: this isn't a valid formula !!!
    tmp <- as.character(model)
    if (tmp[1] == "~" && length(tmp) == 2L) {
      # looks like an unquoted single factor model f =~ something
      lav_msg_warn(
        gettext("model seems to be a formula;
                please enclose the model syntax between quotes"))
      # create model and hope for the best
      model.bis <- paste("f =", paste(tmp, collapse = " "), sep = "")
      flat.model <- lavParseModelString(model.bis)
    } else if (tmp[1] == "~" && length(tmp) == 3L) {
      # looks like a (unquoted) regression formula
      lav_msg_warn(
        gettext("model seems to be a formula;
                please enclose the model syntax between quotes"))
      # create model and hope for the best
      model.bis <- paste(tmp[2], tmp[1], tmp[3])
      flat.model <- lavParseModelString(model.bis)
    } else {
      lav_msg_stop(
        gettext("model seems to be a formula;
                please enclose the model syntax between quotes"))
    }
  } else if (inherits(model, "lavaan")) {
    # hm, a lavaan model; let's try to extract the parameter table
    # and see what happens
    flat.model <- parTable(model)
  } else if (is.list(model)) {
    # a list! perhaps a full parameter table, or an initial flat model,
    # or something else...

	# 1. flat.model already (output of lavParseModelString)?
	if (!is.null(model$lhs) && !is.null(model$op) &&
	    !is.null(model$rhs) && !is.null(model$mod.idx) &&
		!is.null(attr(model, "modifiers"))) {
	  flat.model <- model
    }

    # look for the bare minimum columns: lhs - op - rhs
    else if (!is.null(model$lhs) && !is.null(model$op) &&
      !is.null(model$rhs) && !is.null(model$free)) {
      # ok, we have something that looks like a parameter table
      flat.model <- model

      # fix semTools issue here? for auxiliary() which does not use
      # block column yet
      if (!is.null(flat.model$block)) {
        nn <- length(flat.model$lhs)
        if (length(flat.model$block) != nn) {
          flat.model$block <- flat.model$group
        }
        if (any(is.na(flat.model$block))) {
          flat.model$block <- flat.model$group
        }
      } else if (!is.null(flat.model$group)) {
        flat.model$block <- flat.model$group
      }
    } else {
      bare.minimum <- c("lhs", "op", "rhs", "free")
      missing.idx <- is.na(match(bare.minimum, names(model)))
      missing.txt <- paste(bare.minimum[missing.idx], collapse = ", ")
      lav_msg_stop(
        gettextf("model is a list, but not a parameterTable?
                 missing column(s) in parameter table: [%s]",
                  lav_msg_view(bare.minimum[missing.idx], "none")))
    }
  } else {
    lav_msg_stop(gettext("model is NULL or not a valid type for it!"))
  }

  # Ok, we got a flattened model; usually this a flat.model object, but it
  # could also be an already lavaanified parTable, or a bare-minimum list with
  # lhs/op/rhs/free elements
  flat.model
}

lav_lavaan_step01_ovnames_ovorder <- function(flat.model = NULL,       # nolint
                                              ov.order   = "model",
                                              data       = NULL,
                                              sample.cov = NULL,
                                              slotData   = NULL) {     # nolint
  # set ov.order in lowercase, check if it is "data" or "model",
  #  if not *** error ***
  # if ov.order == "data"
  #   try adapt flat.model via lav_partable_ov_from_data
  #   (** warning ** if this fails)

  # new in 0.6-14
  # if ov.order = "data", it would seem we need to intervene here;
  # ldw 1/3/2024:
  # we do this by adding an attribute "ovda" to flat.model and partable
  ov.order <- tolower(ov.order)
  if (ov.order == "data") {
    flat.model.orig <- flat.model
    try(
      flat.model <- lav_partable_ov_from_data(flat.model,
        data = data,
        sample.cov = sample.cov,
        slotData = slotData
      ),
      silent = TRUE
    )
    if (inherits(flat.model, "try-error")) {
      lav_msg_warn(gettext("ov.order = \"data\" setting failed;
                           switching back to ov.order = \"model\""))
      flat.model <- flat.model.orig
    }
  } else if (ov.order != "model") {
    lav_msg_stop(gettext(
      "ov.order= argument should be \"model\" (default) or \"data\""))
  }

  flat.model
}

lav_lavaan_step01_ovnames_group <- function(flat.model = NULL,        # nolint
                                            ngroups    = 1L) {
  # if "group :" appears in flat.model
  #   tmp.group.values: set of names in corresponding right hand sides
  #   copy flat.model without attributes and call lavaanify,
  #     store result in tmp.lav
  #   extract ov.names, ov.names.y, ov.names.x, lv.names from tmp.lav
  #     via lav_partable_vnames
  # else
  #   if flat.model$group not NULL and more then 1 group.value
  #     extract group.values via lav_partable_group_values
  #     extract, for each group.value,
  #       ov.names, ov.names.y, ov.names.x, lv.names from flat.model
  #       via lav_partable_vnames
  #   else
  #     extract ov.names, ov.names.y, ov.names.x, lv.names from flat.model
  #     via lav_partable_vnames
  #
  #     TODO: call lav_partable_vnames only ones and not for each type

  flat.model.2 <- NULL
  tmp.lav <- NULL
  group.values <- NULL
  ov.names  <- character(0L)
  if (any(flat.model$op == ":" & tolower(flat.model$lhs) == "group")) {
    # here, we only need to figure out:
    # - ngroups
    # - ov's per group
    #
    # - FIXME: we need a more efficient way, avoiding
    #          lavaanify/lav_partable_vnames
    #
    group.idx <- which(flat.model$op == ":" &
                       tolower(flat.model$lhs) == "group")
    # replace by 'group' (in case we got 'Group'):
    flat.model$lhs[group.idx] <- "group"
    tmp.group.values <- unique(flat.model$rhs[group.idx])
    tmp.ngroups <- length(tmp.group.values)

    flat.model.2 <- flat.model
    attr(flat.model.2, "modifiers") <- NULL
    attr(flat.model.2, "constraints") <- NULL
    tmp.lav <- lavaanify(flat.model.2, ngroups = tmp.ngroups, warn = FALSE)
    ov.names <- ov.names.y <- ov.names.x <- lv.names <- vector("list",
      length = tmp.ngroups
    )
    attr(tmp.lav, "vnames") <- lav_partable_vnames(tmp.lav, type = "*")
    for (g in seq_len(tmp.ngroups)) {
      ov.names[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
        type = "ov", group = tmp.group.values[g]
      )))
      ov.names.y[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
        type = "ov.nox", group = tmp.group.values[g]
      )))
      ov.names.x[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
        type = "ov.x", group = tmp.group.values[g]
      )))
      lv.names[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
        type = "lv", group = tmp.group.values[g]
      )))
    }
  } else if (!is.null(flat.model$group)) {
    # user-provided full partable with group column!
    attr(flat.model, "vnames") <- lav_partable_vnames(flat.model, type = "*")
    ngroups <- lav_partable_ngroups(flat.model)
    if (ngroups > 1L) {
      group.values <- lav_partable_group_values(flat.model)
      ov.names <- ov.names.y <- ov.names.x <- lv.names <- vector("list",
        length = ngroups
      )
      for (g in seq_len(ngroups)) {
        # collapsed over levels (if any)
        ov.names[[g]] <- unique(unlist(lav_partable_vnames(flat.model,
          type = "ov", group = group.values[g]
        )))
        ov.names.y[[g]] <- unique(unlist(lav_partable_vnames(flat.model,
          type = "ov.nox", group = group.values[g]
        )))
        ov.names.x[[g]] <- unique(unlist(lav_partable_vnames(flat.model,
          type = "ov.x", group = group.values[g]
        )))
        lv.names[[g]] <- unique(unlist(lav_partable_vnames(flat.model,
          type = "lv", group = group.values[g]
        )))
      }
    } else {
      ov.names   <- lav_partable_vnames(flat.model, type = "ov")
      ov.names.y <- lav_partable_vnames(flat.model, type = "ov.nox")
      ov.names.x <- lav_partable_vnames(flat.model, type = "ov.x")
      lv.names   <- lav_partable_vnames(flat.model, type = "lv")
    }
  } else {
    # collapse over levels (if any)
    attr(flat.model, "vnames") <- lav_partable_vnames(flat.model, type = "*")
    ov.names   <- unique(unlist(lav_partable_vnames(flat.model, type = "ov")))
    ov.names.y <- unique(unlist(lav_partable_vnames(flat.model,
      type = "ov.nox")))
    ov.names.x <- unique(unlist(lav_partable_vnames(flat.model, type = "ov.x")))
    lv.names   <- unique(unlist(lav_partable_vnames(flat.model, type = "lv")))
  }

  # sanity check (new in 0.6-8): do we have any ov.names?
  # detect early
  if (length(unlist(ov.names)) == 0L) {
    lav_msg_stop(
      gettext("ov.names is empty: model does not refer to any observed
              variables; check your syntax."))
  }

  list(
    flat.model   = flat.model,
    ov.names     = ov.names,
    ov.names.x   = ov.names.x,
    ov.names.y   = ov.names.y,
    lv.names     = lv.names,
    group.values = group.values,
    ngroups      = ngroups
  )
}

lav_lavaan_step01_ovnames_checklv <- function(                    # nolint
    lv.names    = character(0L),
    data        = NULL,
    sample.cov  = NULL,
    dotdotdot   = NULL,
    slotOptions = NULL) {                                         # nolint
  # latent variables cannot appear in data --> *** error ***
  #   (except when explicitly requested)
  # latent interactions are not supported ---> *** error ***

  # sanity check: ov.names.x should NOT appear in ov.names.y
  # this may happen if 'x' is exogenous in one block, but not in another...
  # endo.idx <- which(ov.names.x %in% ov.names.y)
  # if (length(endo.idx) > 0L) {
  #    # remove from x! (new in 0.6-8)
  #    ov.names.x <- ov.names.x[-endo.idx]
  # }


  # handle for lv.names that are also observed variables (new in 0.6-6)
  lv.lv.names <- unique(unlist(lv.names))
  if (length(lv.lv.names) > 0L) {
    # check for lv.names in data/cov
    if (!is.null(data)) {
      bad.idx <- which(lv.lv.names %in% names(data))
    } else if (!is.null(sample.cov)) {
      bad.idx <- which(lv.lv.names %in% rownames(sample.cov))
    } else {
      bad.idx <- integer(0L)
    }

        # if found, hard stop
    if (length(bad.idx) > 0L) {
      if (!is.null(dotdotdot$check.lv.names) &&
        !dotdotdot$check.lv.names) {
        # ignore it, user switched this check off -- new in 0.6-7
      } else {
        lav_msg_stop(gettext(
          "some latent variable names collide with observed variable names:"),
          paste(lv.lv.names[bad.idx], collapse = " ")
        )
      }

      # rename latent variables (by adding 'lat')
      # flat.model.idx <- which(flat.model$op == "=~" &
      #                  flat.model$lhs %in% lv.names[bad.idx])
      # flat.model$lhs[flat.model.idx] <-
      #   paste(flat.model$lhs[flat.model.idx], "lat", sep = "")

      # add names to ov.names
      # ov.names <- c(ov.names, lv.names[bad.idx])
      # what about ov.names.y and ov.names.x?
    }
  }

  # sanity check: we do not support latent interaction yet (using the :)
  lv.int.idx <- which(grepl(":", lv.lv.names))
  if (length(lv.int.idx) > 0L) {
    if (!is.null(dotdotdot$check.lv.interaction) &&
      !dotdotdot$check.lv.interaction) {
      # ignore, user (or sam) switched this check off - new in 0.6-16
    } else if (!is.null(slotOptions) && !slotOptions$check.lv.interaction) {
      # ignore
    } else {
      lav_msg_stop(gettextf(
        "Interaction terms involving latent variables (%s) are not supported.
        You may consider creating product indicators to define
        the latent interaction term. See the indProd() function
        in the semTools package.", lv.lv.names[lv.int.idx[1]]))
    }
  }

  invisible(NULL)
}

lav_lavaan_step01_ovnames_namesl <- function(data         = NULL,  # nolint
                                             cluster      = NULL,
                                             flat.model   = NULL,
                                             group.values = NULL,
                                             ngroups      = 1L) {
  # if "level :" appears in flat.model
  #   if data not NULL, cluster must not be NULL, if it is: *** error ***
  #   compute tmp.group.values and tmp.level.values from flat.model
  #   there should be at least 2 levels, if not *** error ***
  #   copy flat.model without attributes and lavaanify -> tmp.lav
  #   check at least 2 levels for tmp.lav, if not *** error ***
  #   compute ov.names.l per group and per level (via lav_partable_vnames
  #     on tmp.lav)
  # else
  #   if lav_partable_nlevels(flat.model) > 0
  #   if data not NULL, cluster must not be NULL, if it is: *** error ***
  #     compute ov.names.l per group and per level (via lav_partable_vnames
  #     on flat.model)
  #   else
  #     there are no levels (ov.names.l = list())

  # handle ov.names.l
  if (any(flat.model$op == ":" & tolower(flat.model$lhs) == "level")) {
    # check for cluster argument
    if (!is.null(data) && is.null(cluster)) {
      lav_msg_stop(gettext("cluster argument is missing."))
    }

    # here, we only need to figure out:
    # - nlevels
    # - ov's per level
    # - FIXME: we need a more efficient way, avoiding lavaanify/vnames

    group.idx <- which(flat.model$op == ":" & flat.model$lhs == "group")
    tmp.group.values <- unique(flat.model$rhs[group.idx])
    tmp.ngroups <- max(c(length(tmp.group.values), 1))

    level.idx <- which(flat.model$op == ":" &
                       tolower(flat.model$lhs) == "level")
    # replace by "level" (in case we got 'Level')
    flat.model$lhs[level.idx] <- "level"
    tmp.level.values <- unique(flat.model$rhs[level.idx])
    tmp.nlevels <- length(tmp.level.values)

    # we need at least 2 levels (for now)
    if (tmp.nlevels < 2L) {
      lav_msg_stop(
        gettext("when data is clustered, you must specify a model
                for each level in the model syntax (for now);
                see example(Demo.twolevel)")
      )
    }

    flat.model.2 <- flat.model
    attr(flat.model.2, "modifiers") <- NULL
    attr(flat.model.2, "constraints") <- NULL
    tmp.lav <- lavaanify(flat.model.2, ngroups = tmp.ngroups, warn = FALSE)
    # check for empty levels
    if (max(tmp.lav$level) < 2L) {
      lav_msg_stop(
        gettext("at least one level has no model syntax;
                you must specify a model for each level in the model syntax;
                see example(Demo.twolevel)")
      )
    }
    ov.names.l <- vector("list", length = tmp.ngroups) # per group

    for (g in seq_len(tmp.ngroups)) {
      ov.names.l[[g]] <- vector("list", length = tmp.nlevels)
      for (l in seq_len(tmp.nlevels)) {
        if (tmp.ngroups > 1L) {
          ov.names.l[[g]][[l]] <-
            unique(unlist(lav_partable_vnames(tmp.lav,
              type = "ov",
              group = tmp.group.values[g],
              level = tmp.level.values[l]
            )))
        } else {
          ov.names.l[[g]][[l]] <-
            unique(unlist(lav_partable_vnames(tmp.lav,
              type = "ov",
              level = tmp.level.values[l]
            )))
        }
      } # levels
    } # groups
  } else {
    # perhaps model is already a parameter table
    nlevels <- lav_partable_nlevels(flat.model)
    if (nlevels > 1L) {
      # check for cluster argument (only if we have data)
      if (!is.null(data) && is.null(cluster)) {
        lav_msg_stop(gettext("cluster argument is missing."))
      }

      ngroups <- lav_partable_ngroups(flat.model)
      group.values <- lav_partable_group_values(flat.model)
      ov.names.l <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        # note: lavNames() will return a list if any level:
        ov.names.l[[g]] <- lavNames(flat.model, "ov", group = group.values[g])
      }
    } else {
      # no level: in model syntax
      ov.names.l <- list()
    }
  }

  list(
    flat.model = flat.model,
    ov.names.l = ov.names.l
  )
}

lav_lavaan_step01_ovnames_ordered <- function(ordered    = NULL,  # nolint
                                              flat.model = NULL,
                                              data       = NULL) {
  # interpretation and check ordered parameter, modify if needed

  # sanity check ordered argument (just in case, add lhs variables names)
  if (!is.null(ordered)) { # new in 0.6-4
    if (is.logical(ordered) && ordered) { # ordered = TRUE
      # assume the user means: ordered = names(Data)
      ordered <- lavNames(flat.model, "ov.nox") # new in 0.6-6: changed from ov
    } else if (is.logical(ordered) && !ordered) {
      ordered <- character(0L)
    } else if (!is.character(ordered)) {
      lav_msg_stop(gettext("ordered argument must be a character vector"))
    } else if (length(ordered) == 1L && nchar(ordered) == 0L) {
      ordered <- character(0L)
    } else {
      # check if all names in "ordered" occur in the dataset?
      if (!is.null(data)) {
        if (inherits(data, "data.frame")) {
          data_names <- names(data)
        } else if (inherits(data, "matrix")) {
          data_names <- colnames(data)
        }
        missing.idx <- which(!ordered %in% data_names)
        if (length(missing.idx) > 0L) { # FIXme: warn = FALSE has no eff
          lav_msg_warn(gettextf(
            "ordered variable(s): %s could not be found
            in the data and will be ignored",
             paste(ordered[missing.idx], collapse = " ")))
        }
      }
    }
  }

  # add the variable names that were treated as ordinal
  # in the model syntax
  ordered <- unique(c(ordered, lavNames(flat.model, "ov.ord")))

  ordered
}
