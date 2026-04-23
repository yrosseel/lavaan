# lavaan parameter table
#
# initial version: YR 22/05/2009
#  major revision: YR 02/11/2010: - FLATTEN the model syntax and turn it into a
#                                   data.frame, with a "modifiers" attribute
#                                 - add default elements here
#                                 - check for duplicate elements
#                                 - allow for every possible model...
#                                 - since 0.4-5
#                                 - the end result is a full description of
#                                   a model (but no matrix representation)
# - 14 Jan 2014: merge 02lavaanUser.R with lav_partable.R
#                move syntax-based code to lav_syntax.R
# - 26 April 2016: handle multiple 'blocks' (levels, classes, groups, ...)
# - 24 March 2019: handle efa sets
# - 23 May   2020: support for random slopes

lav_model_partable  <- function(
                                model = NULL,            # nolint start
                                meanstructure = FALSE,
                                int.ov.free = FALSE,
                                int.lv.free = FALSE,
                                marker.int.zero = FALSE,
                                orthogonal = FALSE,
                                orthogonal.y = FALSE,
                                orthogonal.x = FALSE,
                                orthogonal.efa = FALSE,
                                std.lv = FALSE,
                                correlation = FALSE,
                                composites = TRUE,
                                effect.coding = "",
                                conditional.x = FALSE,
                                fixed.x = FALSE,
                                parameterization = "delta",
                                constraints = NULL,
                                ceq.simple = FALSE,
                                auto = FALSE,
                                model.type = "sem",
                                auto.fix.first = FALSE,
                                auto.fix.single = FALSE,
                                auto.var = FALSE,
                                auto.cov.lv.x = FALSE,
                                auto.cov.y = FALSE,
                                auto.th = FALSE,
                                auto.delta = FALSE,
                                auto.efa = FALSE,
                                varTable = NULL,
                                ngroups = 1L,
                                nthresholds = NULL,
                                group.equal = NULL,
                                group.partial = NULL,
                                group.w.free = FALSE,
                                debug = FALSE,
                                warn = TRUE,
                                as.data.frame. = TRUE) {  # nolint end
  if (!missing(debug)) {
    current_debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current_debug), TRUE)
  }
  if (!missing(warn)) {
    current_warn <- lav_warn()
    if (lav_warn(warn))
      on.exit(lav_warn(current_warn), TRUE)
  }

  # copy arguments which are not snake_case and modified within function
  int_ov_free <- int.ov.free
  int_lv_free <- int.lv.free
  effect_coding <- effect.coding
  ceq_simple <- ceq.simple
  auto_fix_first <- auto.fix.first
  auto_fix_single <- auto.fix.single
  auto_var <- auto.var
  auto_cov_lv_x <- auto.cov.lv.x
  auto_cov_y <- auto.cov.y
  auto_th <- auto.th
  auto_delta <- auto.delta
  auto_efa <- auto.efa
  var_table <- varTable

  # check if model is already flat or a full parameter table
  if (is.list(model) && !is.null(model$lhs)) {
    if (is.null(model$mod.idx)) {
      lav_msg_warn(gettext("input already looks like a parameter table"))
      return(lav_partable_set_cache(model))
    } else {
      flat <- model
    }
  } else {
    # parse the model syntax and flatten the user-specified model
    # return a data.frame, where each line is a model element (rhs, op, lhs)
    flat <- lavParseModelString(
      model.syntax = model, debug = FALSE
    )
  }
  # user-specified *modifiers* are returned as an attribute
  tmp_mod <- attr(flat, "modifiers")
  attr(flat, "modifiers") <- NULL
  # user-specified *constraints* are returned as an attribute
  tmp_con <- attr(flat, "constraints")
  attr(flat, "constraints") <- NULL

  # ov_names_data?
  ov_names_data <- attr(flat, "ovda")

  # extra constraints?
  if (!is.null(constraints) && any(nchar(constraints) > 0L)) {
    flat2 <- lavParseModelString(model.syntax = constraints, warn = lav_warn())
    con2 <- attr(flat2, "constraints")
    rm(flat2)
    tmp_con <- c(tmp_con, con2)
  }
  if (length(tmp_con) > 0L) {
    # add 'user' column
    tmp_con <- lapply(tmp_con, function(x) {
      x$user <- 1L
      x
    })
    # any explicit (in)equality constraints? (ignoring := definitions)
    tmp_con_nondef_flag <- (sum(sapply(tmp_con, "[[", "op")
                            %in% c("==", "<", ">")) > 0L)
    # any explicit equality constraints?
    if (tmp_con_nondef_flag) {
      ceq_simple <- FALSE
    }
  }

  if (lav_debug()) {
    cat("[lavaan DEBUG]: flat (flattened user model):\n")
    print(flat)
    cat("[lavaan DEBUG]: tmp_mod (modifiers):\n")
    print(str(tmp_mod))
    cat("[lavaan DEBUG]: tmp_con (constraints):\n")
    print(str(tmp_con))
  }

  # bogus var_table? (if data.type == "none")
  if (!is.null(var_table)) {
    if (!is.list(var_table) || is.null(var_table$name)) {
      lav_msg_stop(gettext(
        "varTable is not a list or does not contain variable names."))
    }
    if (all(var_table$nobs == 0)) {
      var_table <- NULL
    }
  }

  # check for wrongly specified variances/covariances/intercepts
  # of exogenous variables in model syntax (if fixed.x=TRUE)
  if (fixed.x && lav_warn()) { # we ignore the groups here!
    # we only call this function for the warning message
    tmp <- lav_partable_vnames(flat, "ov.x", force.warn = TRUE)
    rm(tmp)
  }

  # check if group.equal is non-empty, but ngroups = 1L
  # fixme: triggers this if mimic="Mplus"!
  # if(ngroups == 1L && length(group.equal) > 0L) {
  #    warning("lavaan WARNING: group.equal= argument",
  #    " has no effect if no groups are specified.")
  # }

  # auto=TRUE?
  if (auto) { # mimic sem/cfa auto behavior
    if (model.type == "sem") {
      int_ov_free <- TRUE
      int_lv_free <- FALSE
      auto_fix_first <- !std.lv
      auto_fix_single <- TRUE
      auto_var <- TRUE
      auto_cov_lv_x <- TRUE
      auto_cov_y <- TRUE
      auto_th <- TRUE
      auto_delta <- TRUE
      auto_efa <- TRUE
    } else if (model.type == "growth") {
      int_ov_free <- FALSE
      int_lv_free <- TRUE
      auto_fix_first <- !std.lv
      auto_fix_single <- TRUE
      auto_var <- TRUE
      auto_cov_lv_x <- TRUE
      auto_cov_y <- TRUE
      auto_th <- TRUE
      auto_delta <- TRUE
      auto_efa <- TRUE
    }
  }

  # check for meanstructure
  if (any(flat$op == "~1")) {
    meanstructure <- TRUE
  }

  # check for block identifiers in the syntax (op = ":")
  n_block_flat <- length(which(flat$op == ":"))
  # this is NOT the number of blocks (eg group 1: level 1: -> 1 block)

  # for each non-empty `block' in n_block_flat, produce a USER
  if (n_block_flat > 0L) {
    # make sure flat is a data.frame
    flat <- as.data.frame(flat, stringsAsFactors = FALSE)

    # what are the block lhs labels?
    blocks_lhs_all <- tolower(flat$lhs[flat$op == ":"])
    tmp_blocks_lhs <- unique(blocks_lhs_all)

    # if we have group and level, check that group comes first!
    if ("group" %in% tmp_blocks_lhs && "level" %in% tmp_blocks_lhs) {
      group_idx <- which(tmp_blocks_lhs == "group")
      level_idx <- which(tmp_blocks_lhs == "level")
      if (group_idx > level_idx) {
        lav_msg_stop(gettext(
          "levels must be nested within groups (not the other way around)."))
      }
    }

    # block op == ":" indices
    block_op_idx <- which(flat$op == ":")

    # check for wrong spelled 'group' lhs
    if (length(grep("group", tmp_blocks_lhs)) > 1L) {
      lav_msg_warn(gettext("ambiguous block identifiers for group:"),
        lav_msg_view(tmp_blocks_lhs[grep("group", tmp_blocks_lhs)], "none"))
    }

    # no empty :rhs fields allowed!
    if (any(nchar(flat$rhs[block_op_idx]) == 0L)) {
      empty_idx <- nchar(flat$rhs[block_op_idx]) == 0L
      txt <- paste(flat$lhs[block_op_idx][empty_idx], ":")
      lav_msg_stop(gettext(
        "syntax contains block identifiers with missing numbers/labels: "), txt)
    }

    # check for ngroups (ngroups is based on the data!)
    if ("group" %in% tmp_blocks_lhs) {
      # how many group blocks?
      group_block_idx <- flat$op == ":" & flat$lhs == "group"
      n_group_flat <- length(unique(flat$rhs[group_block_idx]))

      if (n_group_flat > 0L && n_group_flat != ngroups) {
        lav_msg_stop(gettextf(
          "syntax defines %1$s groups; data (or argument ngroups)
          suggests %2$s groups", n_group_flat, ngroups))
      }
    }

    # figure out how many 'blocks' we have, and store indices/block_labels
    tmp_block_rhs <- rep("0", length(tmp_blocks_lhs))
    block_id <- 0L
    block_info <- vector("list", length = n_block_flat) # too large
    block_op_idx1 <- c(block_op_idx, nrow(flat) + 1L) # add addition row
    for (block_op in seq_len(n_block_flat)) {
      # fill block_rhs value(s)
      block_lhs <- flat$lhs[block_op_idx1[block_op]]
      block_rhs <- flat$rhs[block_op_idx1[block_op]]
      tmp_block_rhs[which(block_lhs == tmp_blocks_lhs)] <- block_rhs

      # another block identifier?
      if (block_op_idx1[block_op + 1L] - block_op_idx1[block_op] == 1L) {
        next
      }

      # we have a 'block'
      block_id <- block_id + 1L

      # select flat rows for this block
      tmp_idx <- seq.int(
        block_op_idx1[block_op] + 1L,
        block_op_idx1[block_op + 1L] - 1L
      )

      # store info in block_info
      block_info[[block_id]] <- list(
        lhs = tmp_blocks_lhs, # always the same
        rhs = tmp_block_rhs, # for this block
        idx = tmp_idx
      )
    }
    block_info <- block_info[seq_len(block_id)]

    # new in 0.6-12
    # check for blocks with the same block_rhs combination
    # (perhaps added later?)
    # - merge the indices
    # - remove the duplicated blocks
    block_labels <- sapply(lapply(block_info, "[[", "rhs"),
      paste,
      collapse = "."
    )
    nblocks <- length(unique(block_labels))
    if (nblocks < length(block_labels)) {
      # it would appear we have duplicated block_labels -> merge
      dup_idx <- which(duplicated(block_labels))
      for (i in seq_along(dup_idx)) {
        this_dup_idx <- dup_idx[i]
        orig_idx <- which(block_labels == block_labels[this_dup_idx])[1]
        block_info[[orig_idx]]$idx <- c(
          block_info[[orig_idx]]$idx,
          block_info[[this_dup_idx]]$idx
        )
      }
      block_info <- block_info[-dup_idx]
    }

    # split the flat data.frame per `block', create tmp_list
    # for each `block', and rbind them together, adding block columns
    for (block in seq_len(nblocks)) {
      tmp_block_rhs <- block_info[[block]]$rhs
      block_lhs <- block_info[[block]]$lhs[length(tmp_blocks_lhs)] # last one
      block_idx <- block_info[[block]]$idx

      flat_block <- flat[block_idx, ]
      # rm 'block' column (if any) in flat_block
      flat_block$block <- NULL

      # new in 0.6-7: check for random slopes, add them here
      if (block_lhs == "level" &&
        block > 1L && # FIXME: multigroup, multilevel
        !is.null(flat$rv) &&
        any(nchar(flat$rv) > 0L)) {
        lv_names_rv <- unique(flat$rv[nchar(flat$rv) > 0L])
        for (i in seq_along(lv_names_rv)) {
          # add phantom latent variable
          tmp <- flat_block[1, ]
          tmp$lhs <- lv_names_rv[i]
          tmp$op <- "=~"
          tmp$rhs <- lv_names_rv[i]
          tmp$mod.idx <- max(flat$mod.idx) + i
          tmp$fixed <- "0"
          tmp$start <- ""
          tmp$lower <- ""
          tmp$upper <- ""
          tmp$label <- ""
          tmp$prior <- ""
          tmp$efa <- ""
          tmp$rv <- lv_names_rv[i]
          flat_block <- rbind(flat_block, tmp, deparse.level = 0L)
          tmp_mod <- c(tmp_mod, list(list(fixed = 0)))
        }
      }

      # new in 0.6-8: if multilevel, use 'global' ov.names.x
      if (fixed.x && block_lhs == "level") {
        tmp_ov_names_x <- lav_partable_vnames(flat, "ov.x") # global
        ov_names_x_block <- lav_partable_vnames(flat_block, "ov.x")
        if (length(ov_names_x_block) > 0L) {
          idx <- which(!ov_names_x_block %in% tmp_ov_names_x)
          if (length(idx) > 0L) {
            # warn!
            lav_msg_warn(gettextf(
              "the variable(s) [%s] are exogenous at one level, but endogenous
              at another level. These variables will be treated as endogenous,
              and their variances/intercepts will be freely estimated.
              To remove this warning, use fixed.x = FALSE.",
              lav_msg_view(ov_names_x_block[idx], "none")))
            ov_names_x_block <- ov_names_x_block[-idx]
          }
        }
      } else {
        ov_names_x_block <- NULL
      }

      # new in 0.6-12: if multilevel and conditional.x, make sure
      # that 'splitted' exogenous covariates become 'y' variables
      if (conditional.x && block_lhs == "level") {
        if (ngroups == 1L) {
          other_block_names <- lav_partable_vnames(flat, "ov",
            block = seq_len(nblocks)[-block]
          )
        } else {
          # TEST ME
          this_group <- ceiling(block / nlevels)
          blocks_within_group <- (this_group - 1L) * nlevels + seq_len(nlevels)
          other_block_names <- lav_partable_vnames(flat, "ov",
            block = blocks_within_group[-block]
          )
        }
        ov_names_x_block <- lav_partable_vnames(flat_block, "ov.x")
        if (length(ov_names_x_block) > 0L) {
          idx <- which(ov_names_x_block %in% other_block_names)
          if (length(idx) > 0L) {
            ov_names_x_block <- ov_names_x_block[-idx]
          }
        }
      } else {
        ov_names_x_block <- NULL
      }

      list_block <- lav_partable_flat(flat_block,
        blocks = tmp_blocks_lhs,
        block.id = block,
        meanstructure = meanstructure,
        int.ov.free = int_ov_free, int.lv.free = int_lv_free,
        orthogonal = orthogonal, orthogonal.y = orthogonal.y,
        orthogonal.x = orthogonal.x, orthogonal.efa = orthogonal.efa,
        std.lv = std.lv, correlation = correlation, composites = composites,
        conditional.x = conditional.x, fixed.x = fixed.x,
        parameterization = parameterization,
        auto.fix.first = auto_fix_first,
        auto.fix.single = auto_fix_single,
        auto.var = auto_var, auto.cov.lv.x = auto_cov_lv_x,
        auto.cov.y = auto_cov_y, auto.th = auto_th,
        auto.delta = auto_delta, auto.efa = auto_efa,
        varTable = var_table, group.equal = NULL,
        group.w.free = group.w.free, ngroups = 1L,
        nthresholds = nthresholds,
        ov.names.x.block = ov_names_x_block
      )
      list_block <- as.data.frame(list_block, stringsAsFactors = FALSE)

      # add block columns with current values in block_rhs
      for (b in seq_along(tmp_blocks_lhs)) {
        block_lhs <- tmp_blocks_lhs[b]
        block_rhs <- tmp_block_rhs[b]
        list_block[block_lhs] <- rep(block_rhs, length(list_block$lhs))
      }

      if (!exists("tmp_list")) {
        tmp_list <- list_block
      } else {
        list_block$id <- list_block$id + max(tmp_list$id)
        tmp_list <- rbind(tmp_list, list_block)
      }
    }
    tmp_list <- as.list(tmp_list)

    # convert block columns to integers if possible
    for (b in seq_along(tmp_blocks_lhs)) {
      block_lhs <- tmp_blocks_lhs[b]
      block_rhs <- tmp_block_rhs[b]
      tmp <- try(scan(
        text = tmp_list[[block_lhs]], what = integer(),
        quiet = TRUE
      ), silent = TRUE)
      if (inherits(tmp, "integer")) {
        tmp_list[[block_lhs]] <- tmp
      }
    }
  } else {
    tmp_list <- lav_partable_flat(flat,
      blocks = "group",
      meanstructure = meanstructure,
      int.ov.free = int_ov_free, int.lv.free = int_lv_free,
      orthogonal = orthogonal, orthogonal.y = orthogonal.y,
      orthogonal.x = orthogonal.x, orthogonal.efa = orthogonal.efa,
      std.lv = std.lv, correlation = correlation, composites = composites,
      conditional.x = conditional.x, fixed.x = fixed.x,
      parameterization = parameterization,
      auto.fix.first = auto_fix_first, auto.fix.single = auto_fix_single,
      auto.var = auto_var, auto.cov.lv.x = auto_cov_lv_x,
      auto.cov.y = auto_cov_y, auto.th = auto_th,
      auto.delta = auto_delta, auto.efa = auto_efa,
      varTable = var_table, group.equal = group.equal,
      group.w.free = group.w.free,
      ngroups = ngroups, nthresholds = nthresholds
    )
  }
  if (lav_debug()) {
    cat("[lavaan DEBUG]: parameter tmp_list without MODIFIERS:\n")
    print(as.data.frame(tmp_list, stringsAsFactors = FALSE))
  }

  # check ordinal variables
  categorical <- FALSE
  ov_ord <- lav_object_vnames(tmp_list, "ov.ord")
  all_ord <- all(lav_object_vnames(tmp_list, "ov") %in% ov_ord)
  if (length(ov_ord) > 0L) {
    categorical <- TRUE
    ord_var_idx <- which(tmp_list$op == "~~" &
                         tmp_list$lhs == tmp_list$rhs &
                         tmp_list$lhs %in% ov_ord &
                         tmp_list$user == 1L)
    if (parameterization == "delta" && length(ord_var_idx) > 0L) {
      lav_msg_warn(gettextf("variances of ordered variables are ignored when
      parameterization = \"delta\"; please remove them from the model syntax
      or use parameterization = \"theta\"; variables involved are: %s",
      paste(tmp_list$lhs[ord_var_idx], collapse = " ")))
      # force them to be nonfree and set ustart to 1 (new in 0.6-20)
      # later, after we have processes the modifiers
    }
  } # ov_ord

  # handle multilevel-specific constraints
  multilevel <- FALSE
  nlevels <- 1L
  if (!is.null(tmp_list$level)) {
    nlevels <- lav_partable_nlevels(tmp_list)
    if (nlevels > 1L) {
      multilevel <- TRUE
    }
  }
  if (multilevel && any(tmp_list$op == "~1")) {
    # fix ov intercepts for all within ov that also appear at level 2
    # FIXME: not tested with > 2 levels
    ov_names <- lav_partable_vnames(tmp_list, "ov") ## all names
    level_values <- lav_partable_level_values(tmp_list)
    other_names <- tmp_list$lhs[tmp_list$op == "~1" &
      tmp_list$level %in% level_values[-1L] &
      tmp_list$lhs %in% ov_names]
    fix_names_idx <- which(tmp_list$op == "~1" &
      tmp_list$level %in% level_values[1L] &
      tmp_list$lhs %in% other_names)
    if (length(fix_names_idx) > 0L) {
      tmp_list$free[fix_names_idx] <- 0L
      tmp_list$ustart[fix_names_idx] <- 0
    }
  }
  if (multilevel && any(tmp_list$op == "|")) {
    # fix ALL thresholds at level 1
    level_values <- lav_partable_level_values(tmp_list)
    th_idx <- which(tmp_list$op == "|" &
      tmp_list$level %in% level_values[1L])
    tmp_list$free[th_idx] <- 0L
    tmp_list$ustart[th_idx] <- 0

    # fix ALL scaling parmaters at higher levels
    scale_idx <- which(tmp_list$op == "~*~" &
      tmp_list$level %in% level_values[-1L])
    tmp_list$free[scale_idx] <- 0L
    tmp_list$ustart[scale_idx] <- 1
  }

  # apply user-specified modifiers
  warn_about_single_label <- FALSE
  if (length(tmp_mod)) {
    for (el in seq_along(tmp_mod)) {
      idx <- which(tmp_list$mod.idx == el) # for each group

      # 0.5-21: check if idx exists
      # perhaps the corresponding element was duplicated, and removed
      if (length(idx) == 0L) {
        next
      }

      tmp_mod_fixed <- tmp_mod[[el]]$fixed
      tmp_mod_start <- tmp_mod[[el]]$start
      tmp_mod_lower <- tmp_mod[[el]]$lower
      tmp_mod_upper <- tmp_mod[[el]]$upper
      tmp_mod_label <- tmp_mod[[el]]$label
      tmp_mod_prior <- tmp_mod[[el]]$prior
      tmp_mod_efa <- tmp_mod[[el]]$efa
      tmp_mod_rv <- tmp_mod[[el]]$rv

      # check for single argument if multiple groups
      if (ngroups > 1L && length(idx) > 1L) {
        # Ok, this is not very consistent:
        # A) here we force same behavior across groups
        if (length(tmp_mod_fixed) == 1L) {
          tmp_mod_fixed <- rep(tmp_mod_fixed, ngroups)
        }
        if (length(tmp_mod_start) == 1L) {
          tmp_mod_start <- rep(tmp_mod_start, ngroups)
        }
        if (length(tmp_mod_lower) == 1L) {
          tmp_mod_lower <- rep(tmp_mod_lower, ngroups)
        }
        if (length(tmp_mod_upper) == 1L) {
          tmp_mod_upper <- rep(tmp_mod_upper, ngroups)
        }
        if (length(tmp_mod_prior) == 1L) {
          tmp_mod_prior <- rep(tmp_mod_prior, ngroups)
        }
        if (length(tmp_mod_efa) == 1L) {
          tmp_mod_efa <- rep(tmp_mod_efa, ngroups)
        }
        if (length(tmp_mod_rv) == 1L) {
          tmp_mod_rv <- rep(tmp_mod_rv, ngroups)
        }

        # new in 0.6-7 (proposal):
        # - always recycle modifiers, including labels
        # - if ngroups > 1 AND group.label= is empty, produce a warning
        #   (as this is a break from < 0.6-6)
        if (length(tmp_mod_label) == 1L) {
          tmp_mod_label <- rep(tmp_mod_label, ngroups)
          if (is.null(group.equal) || length(group.equal) == 0L) {
            warn_about_single_label <- TRUE
          }
        }

        # < 0.6-7 code:
        # B) here we do NOT! otherwise, it would imply an equality
        #                    constraint...
        #    except if group.equal="loadings"!
        # if(length(tmp_mod_label) == 1L) {
        #    if("loadings" %in% group.equal ||
        #       "composite.loadings" %in% group.equal) {
        #        tmp_mod_label <- rep(tmp_mod_label, ngroups)
        #    } else {
        #        tmp_mod_label <- c(tmp_mod_label, rep("", (ngroups-1L)) )
        #    }
        # }
      }

      # check for wrong number of arguments if multiple groups
      nidx <- length(idx)
      if ((!is.null(tmp_mod_fixed) && nidx != length(tmp_mod_fixed)) ||
        (!is.null(tmp_mod_start) && nidx != length(tmp_mod_start)) ||
        (!is.null(tmp_mod_lower) && nidx != length(tmp_mod_lower)) ||
        (!is.null(tmp_mod_upper) && nidx != length(tmp_mod_upper)) ||
        (!is.null(tmp_mod_prior) && nidx != length(tmp_mod_prior)) ||
        (!is.null(tmp_mod_efa) && nidx != length(tmp_mod_efa)) ||
        (!is.null(tmp_mod_rv) && nidx != length(tmp_mod_rv)) ||
        (!is.null(tmp_mod_label) && nidx != length(tmp_mod_label))) {
        el_idx <- which(tmp_list$mod.idx == el)[1L]
        lav_msg_stop(gettextf(
          "wrong number of arguments in modifier (%s) of element",
          lav_msg_view(tmp_mod_label, "none")),
          tmp_list$lhs[el_idx], tmp_list$op[el_idx], tmp_list$rhs[el_idx]
        )
      }

      # apply modifiers
      if (!is.null(tmp_mod_fixed)) {
        # two options: constant or NA
        na_idx <- which(is.na(tmp_mod_fixed))
        not_na_idx <- which(!is.na(tmp_mod_fixed))

        # constant
        tmp_list$ustart[idx][not_na_idx] <- tmp_mod_fixed[not_na_idx]
        tmp_list$free[idx][not_na_idx] <- 0L

        # NA* modifier
        tmp_list$free[idx][na_idx] <- 1L # eg factor loading
        tmp_list$ustart[idx][na_idx] <- as.numeric(NA)
      }
      if (!is.null(tmp_mod_start)) {
        tmp_list$ustart[idx] <- tmp_mod_start
      }
      if (!is.null(tmp_mod_prior)) {
        # do we already have a `prior' column? if not, create one
        if (is.null(tmp_list$prior)) {
          tmp_list$prior <- character(length(tmp_list$lhs))
        }
        tmp_list$prior[idx] <- tmp_mod_prior
      }
      if (!is.null(tmp_mod_efa)) {
        # do we already have a `efa' column? if not, create one
        if (is.null(tmp_list$efa)) {
          tmp_list$efa <- character(length(tmp_list$lhs))
        }
        tmp_list$efa[idx] <- tmp_mod_efa
      }
      if (!is.null(tmp_mod_rv)) {
        # do we already have a `rv' column? if not, create one
        if (is.null(tmp_list$rv)) {
          tmp_list$rv <- character(length(tmp_list$lhs))
        }
        tmp_list$rv[idx] <- tmp_mod_rv

        tmp_list$free[idx] <- 0L
        tmp_list$ustart[idx] <- as.numeric(NA) #
      }
      if (!is.null(tmp_mod_lower)) {
        # do we already have a `lower' column? if not, create one
        if (is.null(tmp_list$lower)) {
          tmp_list$lower <- rep(-Inf, length(tmp_list$lhs))
        }
        tmp_list$lower[idx] <- as.numeric(tmp_mod_lower)
      }
      if (!is.null(tmp_mod_upper)) {
        # do we already have a `upper' column? if not, create one
        if (is.null(tmp_list$upper)) {
          tmp_list$upper <- rep(Inf, length(tmp_list$lhs))
        }
        tmp_list$upper[idx] <- as.numeric(tmp_mod_upper)
      }
      if (!is.null(tmp_mod_label)) {
        tmp_list$label[idx] <- tmp_mod_label
      }
    }
  }
  # remove mod.idx column
  tmp_list$mod.idx <- NULL

  # categorical: check for nonfree variances if parameterization = "delta"
  # we already gave warning; here, we force them to be nonfree
  if (categorical && parameterization == "delta") {
    ord_var_idx <- which(tmp_list$op == "~~" &
                         tmp_list$lhs == tmp_list$rhs &
                         tmp_list$lhs %in% ov_ord &
                         tmp_list$user == 1L)
    if (length(ord_var_idx) > 0L) {
      # force them to be nonfree and set ustart to 1 (new in 0.6-20)
      tmp_list$free[ord_var_idx]   <- rep(0L, length(ord_var_idx))
      tmp_list$ustart[ord_var_idx] <- rep(1,  length(ord_var_idx))
    }
  } # categorical


  # warning about single label in multiple group setting?
  if (warn_about_single_label) {
    lav_msg_warn(gettext(
      "using a single label per parameter in a multiple group setting implies
      imposing equality constraints across all the groups; If this is not
      intended, either remove the label(s), or use a vector of labels (one for
      each group); See the Multiple groups section in the man page of
      model.syntax."
    ))
  }

  # if lower/upper values were added, fix non-free values to ustart values
  # new in 0.6-6
  if (!is.null(tmp_list$lower)) {
    fixed_idx <- which(tmp_list$free == 0L)
    if (length(fixed_idx) > 0L) {
      tmp_list$lower[fixed_idx] <- tmp_list$ustart[fixed_idx]
    }
  }
  if (!is.null(tmp_list$upper)) {
    fixed_idx <- which(tmp_list$free == 0L)
    if (length(fixed_idx) > 0L) {
      tmp_list$upper[fixed_idx] <- tmp_list$ustart[fixed_idx]
    }
  }

  # if rv column is present, add rv_names to ALL rows where they are used
  if (!is.null(tmp_list$rv)) {
    rv_names <- unique(tmp_list$rv[nchar(tmp_list$rv) > 0L])
    for (i in seq_along(rv_names)) {
      lhs_idx <- which(tmp_list$lhs == rv_names[i] &
        tmp_list$op == "=~")
      tmp_list$rv[lhs_idx] <- rv_names[i]
    }
  }

  if (lav_debug()) {
    cat("[lavaan DEBUG]: parameter tmp_list with MODIFIERS:\n")
    print(as.data.frame(tmp_list, stringsAsFactors = FALSE))
  }

  # get 'virtual' parameter labels
  if (n_block_flat > 1L) {
    blocks <- tmp_blocks_lhs
  } else {
    blocks <- "group"
  }
  label <- lav_partable_labels(
    partable = tmp_list,
    blocks = blocks,
    group.equal = group.equal,
    group.partial = group.partial
  )

  if (lav_debug()) {
    cat("[lavaan DEBUG]: parameter tmp_list with LABELS:\n")
    tmp <- tmp_list
    tmp$label <- label
    print(as.data.frame(tmp, stringsAsFactors = FALSE))
  }

  # handle EFA equality constraints
  # YR 14 Jan 2020: 0.6-6 does no longer impose 'explicit' constraints
  #                 if we only need to fix a parameter to 0/1
  # Note: we should also check if they are really needed:
  #       eg., if all the factor-loadings of the 'second' set (time/group)
  #       are constrained to be equal to the factor-loadings of the first
  #       set, no further constraints are needed
  if (auto_efa && !is.null(tmp_list$efa)) {
    tmp_list <- lav_partable_efa_constraints(
      LIST = tmp_list,
      orthogonal.efa = orthogonal.efa,
      group.equal = group.equal
    )
  } # auto_efa

  # handle user-specified equality constraints
  # lavaan 0.6-11:
  #     two settings:
  #       1) simple equality constraints ONLY -> back to basics: only
  #          duplicate 'free' numbers; no longer explicit == rows with plabels
  #       2) mixture of simple and other (explicit) constraints
  #          treat them together as we did in <0.6-11
  tmp_list$plabel <- paste(".p", tmp_list$id, ".", sep = "")
  eq_labels <- unique(label[duplicated(label)])
  eq_id <- integer(length(tmp_list$lhs))
  for (eq_label in eq_labels) {
    tmp_con_idx <- length(tmp_con)
    all_idx <- which(label == eq_label) # all same-label parameters
    ref_idx <- all_idx[1L] # the first one only
    other_idx <- all_idx[-1L] # the others
    eq_id[all_idx] <- ref_idx

    # new in 0.6-6: make sure lower/upper constraints are equal too
    if (!is.null(tmp_list$lower) &&
      length(unique(tmp_list$lower[all_idx])) > 0L) {
      non_inf <- which(is.finite(tmp_list$lower[all_idx]))
      if (length(non_inf) > 0L) {
        smalles_val <- min(tmp_list$lower[all_idx][non_inf])
        tmp_list$lower[all_idx] <- smalles_val
      }
    }
    if (!is.null(tmp_list$upper) &&
      length(unique(tmp_list$upper[all_idx])) > 0L) {
      non_inf <- which(is.finite(tmp_list$upper[all_idx]))
      if (length(non_inf) > 0L) {
        largest_val <- max(tmp_list$upper[all_idx][non_inf])
        tmp_list$upper[all_idx] <- largest_val
      }
    }

    # two possibilities:
    #   1. all_idx contains a fixed parameter: in this case,
    #      we fix them all (hopefully to the same value)
    #   2. all_idx contains only free parameters

    # 1. all_idx contains a fixed parameter
    if (any(tmp_list$free[all_idx] == 0L)) {
      # which one is fixed?
      fixed_all <- all_idx[tmp_list$free[all_idx] == 0L]
      # only pick the first
      fixed_idx <- fixed_all[1]

      # sanity check: are all ustart values equal?
      ustart1 <- tmp_list$ustart[fixed_idx]
      if (all(is.na(tmp_list$ustart[fixed_all]))) {
        # nothing to do; ustart values have not been set yet
      } else if (!all(ustart1 == tmp_list$ustart[fixed_all])) {
        lav_msg_warn(gettext(
          "equality constraints involve fixed parameters with different values;
          only the first one will be used"))
      }

      # make them all fixed
      tmp_list$ustart[all_idx] <- tmp_list$ustart[fixed_idx]
      tmp_list$free[all_idx] <- 0L # not free anymore, since it must
      # be equal to the 'fixed' parameter
      # (Note: Mplus ignores this)
      eq_id[all_idx] <- 0L # remove from eq_id list

      # new in 0.6-8 (for efa + user-specified eq constraints)
      if (any(tmp_list$user[all_idx] %in% c(7L, 77L))) {
        # if involved in an efa block, store in tmp_con anyway
        # we may need it for the rotated solution
        for (o in other_idx) {
          tmp_con_idx <- tmp_con_idx + 1L
          tmp_con[[tmp_con_idx]] <- list(
            op = "==",
            lhs = tmp_list$plabel[ref_idx],
            rhs = tmp_list$plabel[o],
            user = 2L
          )
        }
      }
    } else {
      # 2. all_idx contains only free parameters
      # old system:
      # - add tmp_con entry
      # - in 0.6-11: only if tmp_con is not empty
      if (!ceq_simple) {
        for (o in other_idx) {
          tmp_con_idx <- tmp_con_idx + 1L
          tmp_con[[tmp_con_idx]] <- list(
            op = "==",
            lhs = tmp_list$plabel[ref_idx],
            rhs = tmp_list$plabel[o],
            user = 2L
          )
        }
      } else {
        # new system:
        # - set $free elements to zero, and later to ref id
        tmp_list$free[other_idx] <- 0L # all but the first are non-free
        # but will get a duplicated number
      }

      # just to trick semTools, also add something in the label
      # colum, *if* it is empty
      # update: 0.6-11 we keep this, because it shows the plabels
      #         when eg group.equal = "loadings"
      for (i in all_idx) {
        if (nchar(tmp_list$label[i]) == 0L) {
          tmp_list$label[i] <- tmp_list$plabel[ref_idx]
        }
      }
    } # all free
  } # eq in eq_labels
  if (lav_debug()) {
    print(tmp_con)
  }



  # handle constraints (if any) (NOT per group, but overall - 0.4-11)
  if (length(tmp_con) > 0L) {
    n_con <- length(tmp_con)
    tmp_idx <- length(tmp_list$id) + seq_len(n_con)
    # grow tmp_list with length(tmp_con) extra rows
    tmp_list <- lapply(tmp_list, function(x) {
      if (is.character(x)) {
        c(x, rep("", n_con))
      } else {
        c(x, rep(NA, n_con))
      }
    })

    # fill in some columns
    tmp_list$id[tmp_idx] <- tmp_idx
    tmp_list$lhs[tmp_idx] <- unlist(lapply(tmp_con, "[[", "lhs"))
    tmp_list$op[tmp_idx] <- unlist(lapply(tmp_con, "[[", "op"))
    tmp_list$rhs[tmp_idx] <- unlist(lapply(tmp_con, "[[", "rhs"))
    tmp_list$user[tmp_idx] <- unlist(lapply(tmp_con, "[[", "user"))

    # zero is nicer?
    tmp_list$free[tmp_idx] <- rep(0L, n_con)
    tmp_list$exo[tmp_idx] <- rep(0L, n_con)
    tmp_list$block[tmp_idx] <- rep(0L, n_con)

    if (!is.null(tmp_list$group)) {
      if (is.character(tmp_list$group)) {
        tmp_list$group[tmp_idx] <- rep("", n_con)
      } else {
        tmp_list$group[tmp_idx] <- rep(0L, n_con)
      }
    }
    if (!is.null(tmp_list$level)) {
      if (is.character(tmp_list$level)) {
        tmp_list$level[tmp_idx] <- rep("", n_con)
      } else {
        tmp_list$level[tmp_idx] <- rep(0L, n_con)
      }
    }
    if (!is.null(tmp_list$class)) {
      if (is.character(tmp_list$class)) {
        tmp_list$class[tmp_idx] <- rep("", n_con)
      } else {
        tmp_list$class[tmp_idx] <- rep(0L, n_con)
      }
    }
  }

  # check defined variables (:=)
  def_idx <- which(tmp_list$op == ":=")
  if (length(def_idx) > 0L) {
    # check if the lhs is unique (new in 0.6-20)
    def_lhs <- tmp_list$lhs[def_idx]
    dup_idx <- which(duplicated(def_lhs))
    if (length(dup_idx) > 0L) {
      # warn or stop? warn for now
      lav_msg_warn(gettextf("at least one defined variable (using the :=
                             operator) has been duplicated, and will be
                             overwritten by the last one: %s",
                            paste(def_lhs[dup_idx], collapse = " ")))
    }
    # put lhs of := elements in label column
    tmp_list$label[def_idx] <- def_lhs
  }


  # handle effect_coding related equality constraints
  if (is.logical(effect_coding) && effect_coding) {
    effect_coding <- c("loadings", "intercepts")
  } else if (!is.character(effect_coding)) {
    lav_msg_stop(gettext("effect.coding argument must be a character string"))
  }
  # in ordinal models, do integer coding
  if (all_ord) effect_coding <- c(effect_coding, "thresholds")
  if (any(c("loadings", "intercepts") %in% effect_coding)) {
    tmp <- list()
    # for each block
    nblocks <- lav_partable_nblocks(tmp_list)
    for (b in seq_len(nblocks)) {

      # which group?
      this_group <- floor(b / nlevels + 0.5)

      # lv's for this block/set
      lv_names <- unique(tmp_list$lhs[tmp_list$op == "=~" &
        tmp_list$block == b])

      if (length(lv_names) == 0L) {
        next
      }

      for (lv in lv_names) {
        # ind_names
        ind_names <- tmp_list$rhs[tmp_list$op == "=~" &
          tmp_list$block == b &
          tmp_list$lhs == lv]

        if ("loadings" %in% effect_coding &&
            (!"loadings" %in% group.equal || this_group == 1L)) {
          # factor loadings indicators of this lv
          loadings_idx <- which(tmp_list$op == "=~" &
            tmp_list$block == b &
            tmp_list$rhs %in% ind_names &
            tmp_list$lhs == lv)

          # all free?
          if (length(loadings_idx) > 0L &&
            all(tmp_list$free[loadings_idx] > 0L)) {
            # add eq constraint
            plabel <- tmp_list$plabel[loadings_idx]

            # Note:  we write them as
            # .p1. == 3 - .p2. - .p3.
            # instead of
            # 3 ==  .p1.+.p2.+.p3.
            # as this makes it easier to translate things to
            # JAGS/stan

            tmp_lhs <- plabel[1]
            if (length(loadings_idx) > 1L) {
              tmp_rhs <- paste(length(loadings_idx), "-",
                paste(plabel[-1], collapse = "-"),
                sep = ""
              )
            } else {
              tmp_rhs <- length(loadings_idx)
            }

            tmp$lhs <- c(tmp$lhs, tmp_lhs)
            tmp$op <- c(tmp$op, "==")
            tmp$rhs <- c(tmp$rhs, tmp_rhs)
            tmp$block <- c(tmp$block, 0L)
            tmp$user <- c(tmp$user, 2L)
            tmp$ustart <- c(tmp$ustart, as.numeric(NA))
          }
        } # loadings

        if ("thresholds" %in% effect_coding &&
            (!"thresholds" %in% group.equal || this_group == 1L)) {
          thresholds_idx <- which(tmp_list$op == "|" &
            tmp_list$block == b &
            tmp_list$lhs %in% ind_names)

          # all free?
          if (length(thresholds_idx) > 0L &&
            all(tmp_list$free[thresholds_idx] > 0L)) {
            nlevs <- table(tmp_list$lhs[thresholds_idx]) + 1
            kmax <- max(nlevs)
            for (ind in ind_names) {
              # 1) fix bottom and top per ov
              trows <- which(tmp_list$lhs[thresholds_idx] == ind)
              thisidx <- thresholds_idx[trows[1]]
              tmp_list$free[thisidx] <- 0L
              tmp_list$ustart[thisidx] <- round(.5 + kmax / nlevs[ind], 3)
              tmp_list$user[thisidx] <- 2L

              if (nlevs[ind] > 2L) {
                thisidx <- thresholds_idx[trows[length(trows)]]
                tmp_list$free[thisidx] <- 0L
                tmp_list$ustart[thisidx] <- round(.5 + kmax * (nlevs[ind] - 1) /
                    nlevs[ind], 3)
                tmp_list$user[thisidx] <- 2L
              }

              # 2) free intercepts of variables with > 2 levels,
              #                     only if automatically fixed
              intercept_idx <- which(tmp_list$op == "~1" &
                tmp_list$block == b &
                tmp_list$lhs == ind)
              if (nlevs[ind] > 2L && length(intercept_idx) > 0L &&
                  tmp_list$user[intercept_idx] == 0L) {
                tmp_list$free[intercept_idx] <- 1L
              }

              # 3) free variances/deltas of variables with > 2 levels,
              #                           only if automatically fixed
              varop <- ifelse(parameterization == "delta", "~*~", "~~")
              var_idx <- which(tmp_list$op == varop &
                            tmp_list$block == b &
                            tmp_list$lhs == ind &
                            tmp_list$lhs == tmp_list$rhs)
              if (length(var_idx) > 0L &&
                  tmp_list$user[var_idx] == 0L) {
                tmp_list$free[var_idx] <- 1L
              }
            }
            # 4) fix latent variance if all variables have 2 levels
            if (all(nlevs == 2L)) {
              lv_var_idx <- which(tmp_list$op == "~~" &
                              tmp_list$block == b &
                              tmp_list$lhs == lv &
                              tmp_list$lhs == tmp_list$rhs)
              if (length(lv_var_idx) > 0L && tmp_list$user[lv_var_idx] == 0L) {
                tmp_list$free[lv_var_idx] <- 0L
                tmp_list$ustart[lv_var_idx] <- 1L
              }
            }
          }
        } # thresholds

        if ("intercepts" %in% effect_coding &&
            (!"intercepts" %in% group.equal || this_group == 1L)) {
          # intercepts for indicators of this lv
          intercepts_idx <- which(tmp_list$op == "~1" &
            tmp_list$block == b &
            tmp_list$lhs %in% ind_names)

          # all free?
          if (length(intercepts_idx) > 0L &&
            (all(tmp_list$free[intercepts_idx] > 0L) ||
              "thresholds" %in% effect_coding)) {
            # 1) add eq constraint
            plabel <- tmp_list$plabel[intercepts_idx]

            tmp_lhs <- plabel[1]
            if (length(intercepts_idx) > 1L) {
              tmp_rhs <- paste("0-",
                paste(plabel[-1], collapse = "-"),
                sep = ""
              )
            } else {
              tmp_rhs <- 0L
            }

            tmp$lhs <- c(tmp$lhs, tmp_lhs)
            tmp$op <- c(tmp$op, "==")
            tmp$rhs <- c(tmp$rhs, tmp_rhs)
            tmp$block <- c(tmp$block, 0L)
            tmp$user <- c(tmp$user, 2L)
            tmp$ustart <- c(tmp$ustart, as.numeric(NA))

            # 2) release latent mean
            lv_int_idx <- which(tmp_list$op == "~1" &
              tmp_list$block == b &
              tmp_list$lhs == lv)
            # free only if automatically added
            if (length(lv_int_idx) > 0L &&
              tmp_list$user[lv_int_idx] == 0L) {
              tmp_list$free[lv_int_idx] <- 1L
            }
          }
        } # intercepts
      } # lv
    } # blocks

    tmp_list <- lav_partable_merge(tmp_list, tmp)
  }

  # marker.int.zero
  if (meanstructure && marker.int.zero) {
    # for each block
    nblocks <- lav_partable_nblocks(tmp_list)
    for (b in seq_len(nblocks)) {
      # lv's for this block/set
      lv_names <- lav_partable_vnames(tmp_list,
        type = "lv.regular",
        block = b
      )
      lv_marker <- lav_partable_vnames(tmp_list,
        type = "lv.marker",
        block = b
      )
      ov_num <- lav_partable_vnames(tmp_list,
        type = "ov.num",
        block = b
      )
      # new in 0.6-22: only consider ov_num markers
      lv_marker <- lv_marker[lv_marker %in% ov_num]

      if (length(lv_marker) == 0L) {
        next
      }

      # fix marker intercepts to zero
      marker_idx <- which(tmp_list$op == "~1" &
        tmp_list$lhs %in% lv_marker & tmp_list$block == b &
        tmp_list$user == 0L)
      tmp_list$free[marker_idx] <- 0L
      tmp_list$ustart[marker_idx] <- 0

      # free latent means
      lv_idx <- which(tmp_list$op == "~1" &
        tmp_list$lhs %in% lv_names & tmp_list$block == b &
        tmp_list$user == 0L)
      tmp_list$free[lv_idx] <- 1L
      tmp_list$ustart[lv_idx] <- as.numeric(NA)
    } # block
  }


  # mg.lv.variances
  if (ngroups > 1L && "mg.lv.variances" %in% effect_coding) {
    tmp <- list()

    # do not include 'EFA' lv's
    if (!is.null(tmp_list$efa)) {
      lv_names <- unique(tmp_list$lhs[tmp_list$op == "=~" &
        !nchar(tmp_list$efa) > 0L])
    } else {
      lv_names <- unique(tmp_list$lhs[tmp_list$op == "=~"])
    }
    group_values <- lav_partable_group_values(tmp_list)

    for (lv in lv_names) {
      # factor variances
      lv_var_idx <- which(tmp_list$op == "~~" &
        tmp_list$lhs == lv &
        tmp_list$rhs == tmp_list$lhs &
        tmp_list$lhs == lv)

      # all free (but the first?)
      if (length(lv_var_idx) > 0L &&
        all(tmp_list$free[lv_var_idx][-1] > 0L)) {
        # 1) add eq constraint
        plabel <- tmp_list$plabel[lv_var_idx]

        tmp_lhs <- plabel[1]
        if (length(lv_var_idx) > 1L) {
          tmp_rhs <- paste(length(lv_var_idx), "-",
            paste(plabel[-1], collapse = "-"),
            sep = ""
          )
        } else {
          tmp_rhs <- length(lv_var_idx)
        }

        tmp$lhs <- c(tmp$lhs, tmp_lhs)
        tmp$op <- c(tmp$op, "==")
        tmp$rhs <- c(tmp$rhs, tmp_rhs)
        tmp$block <- c(tmp$block, 0L)
        tmp$user <- c(tmp$user, 2L)
        tmp$ustart <- c(tmp$ustart, as.numeric(NA))

        # 2) free lv variances first group
        lv_var_g1_idx <- which(tmp_list$op == "~~" &
          tmp_list$group == group_values[1] &
          tmp_list$lhs == lv &
          tmp_list$rhs == tmp_list$lhs &
          tmp_list$lhs == lv)
        # free only if automatically added
        if (length(lv_var_g1_idx) > 0L &&
          tmp_list$user[lv_var_g1_idx] == 0L) {
          tmp_list$free[lv_var_g1_idx] <- 1L
        }
      }
    } # lv

    tmp_list <- lav_partable_merge(tmp_list, tmp)
  }

  # mg.lv.efa.variances
  if (ngroups > 1L && "mg.lv.efa.variances" %in% effect_coding) {
    tmp <- list()

    # only 'EFA' lv's
    if (!is.null(tmp_list$efa)) {
      lv_names <- unique(tmp_list$lhs[tmp_list$op == "=~" &
        nchar(tmp_list$efa) > 0L])
    } else {
      lv_names <- character(0L)
    }
    group_values <- lav_partable_group_values(tmp_list)

    for (lv in lv_names) {
      # factor variances
      lv_var_idx <- which(tmp_list$op == "~~" &
        tmp_list$lhs == lv &
        tmp_list$rhs == tmp_list$lhs &
        tmp_list$lhs == lv)

      # all free (but the first?)
      if (length(lv_var_idx) > 0L &&
        all(tmp_list$free[lv_var_idx][-1] > 0L)) {
        # 1) add eq constraint
        plabel <- tmp_list$plabel[lv_var_idx]

        tmp_lhs <- plabel[1]
        if (length(lv_var_idx) > 1L) {
          tmp_rhs <- paste(length(lv_var_idx), "-",
            paste(plabel[-1], collapse = "-"),
            sep = ""
          )
        } else {
          tmp_rhs <- length(lv_var_idx)
        }

        tmp$lhs <- c(tmp$lhs, tmp_lhs)
        tmp$op <- c(tmp$op, "==")
        tmp$rhs <- c(tmp$rhs, tmp_rhs)
        tmp$block <- c(tmp$block, 0L)
        tmp$user <- c(tmp$user, 2L)
        tmp$ustart <- c(tmp$ustart, as.numeric(NA))

        # 2) free lv variances first group
        lv_var_g1_idx <- which(tmp_list$op == "~~" &
          tmp_list$group == group_values[1] &
          tmp_list$lhs == lv &
          tmp_list$rhs == tmp_list$lhs &
          tmp_list$lhs == lv)
        # free only if automatically added
        if (length(lv_var_g1_idx) > 0L &&
          tmp_list$user[lv_var_g1_idx] == 0L) {
          tmp_list$free[lv_var_g1_idx] <- 1L
        }
      }
    } # lv

    tmp_list <- lav_partable_merge(tmp_list, tmp)
  }


  # count free parameters
  idx_free <- which(tmp_list$free > 0L)
  tmp_list$free[idx_free] <- seq_along(idx_free)

  # new in 0.6-11: add free counter to this element (as in < 0.5-18)
  # unless we have other constraints
  if (ceq_simple) {
    idx_equal <- which(eq_id > 0)
    tmp_list$free[idx_equal] <- tmp_list$free[eq_id[idx_equal]]
  }

  # new in 0.6-14: add 'da' entries to reflect data-based order of ov's
  # now via attribute "ovda"
  attr(tmp_list, "ovda") <- ov_names_data

  # backwards compatibility...
  if (!is.null(tmp_list$unco)) {
    tmp_list$unco[idx_free] <- seq_along(sum(tmp_list$free > 0L))
  }

  if (lav_debug()) {
    cat("[lavaan DEBUG] lavParTable\n")
    print(as.data.frame(tmp_list))
  }

  # data.frame?
  if (as.data.frame.) {
    tmp_list <- as.data.frame(tmp_list, stringsAsFactors = FALSE)
    attr(tmp_list, "ovda") <- ov_names_data
  } else {
    tmp_list <- lav_partable_set_cache(tmp_list) # add cached "pta" data
  }

  tmp_list
}
lavParTable <- lav_model_partable   # synonym #nolint
lavaanify <- lav_model_partable     # synonym
