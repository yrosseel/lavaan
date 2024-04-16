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

lavaanify <- lavParTable <- function( # nolint
                                     model = NULL,
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
                                     varTable = NULL, # nolint
                                     ngroups = 1L,
                                     nthresholds = NULL,
                                     group.equal = NULL,
                                     group.partial = NULL,
                                     group.w.free = FALSE,
                                     debug = FALSE,
                                     warn = TRUE,
                                     as.data.frame. = TRUE) { # nolint
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
      model.syntax = model, warn = warn,
      debug = FALSE
    )
  }
  # user-specified *modifiers* are returned as an attribute
  tmp.mod <- attr(flat, "modifiers")
  attr(flat, "modifiers") <- NULL
  # user-specified *constraints* are returned as an attribute
  tmp.con <- attr(flat, "constraints")
  attr(flat, "constraints") <- NULL

  # ov.names.data?
  ov.names.data <- attr(flat, "ovda")

  # extra constraints?
  if (!is.null(constraints) && any(nchar(constraints) > 0L)) {
    flat2 <- lavParseModelString(model.syntax = constraints, warn = warn)
    con2 <- attr(flat2, "constraints")
    rm(flat2)
    tmp.con <- c(tmp.con, con2)
  }
  if (length(tmp.con) > 0L) {
    # add 'user' column
    tmp.con <- lapply(tmp.con, function(x) {
      x$user <- 1L
      x
    })
    # any explicit (in)equality constraints? (ignoring := definitions)
    tmp.con.nondef.flag <- (sum(sapply(tmp.con, "[[", "op")
    %in% c("==", "<", ">")) > 0L)
    # any explicit equality constraints?
    tmp.con.eq.flag <- (sum(sapply(tmp.con, "[[", "op") == "==") > 0L)
    if (tmp.con.nondef.flag) {
      ceq.simple <- FALSE
    }
  }

  if (debug) {
    cat("[lavaan DEBUG]: flat (flattened user model):\n")
    print(flat)
    cat("[lavaan DEBUG]: tmp.mod (modifiers):\n")
    print(str(tmp.mod))
    cat("[lavaan DEBUG]: tmp.con (constraints):\n")
    print(str(tmp.con))
  }

  # bogus varTable? (if data.type == "none")
  if (!is.null(varTable)) {
    if (!is.list(varTable) || is.null(varTable$name)) {
      lav_msg_stop(gettext(
        "varTable is not a list or does not contain variable names."))
    }
    if (all(varTable$nobs == 0)) {
      varTable <- NULL # nolint
    }
  }

  # check for wrongly specified variances/covariances/intercepts
  # of exogenous variables in model syntax (if fixed.x=TRUE)
  if (fixed.x && warn) { # we ignore the groups here!
    # we only call this function for the warning message
    tmp <- lav_partable_vnames(flat, "ov.x", warn = TRUE)
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
      int.ov.free <- TRUE
      int.lv.free <- FALSE
      auto.fix.first <- !std.lv
      auto.fix.single <- TRUE
      auto.var <- TRUE
      auto.cov.lv.x <- TRUE
      auto.cov.y <- TRUE
      auto.th <- TRUE
      auto.delta <- TRUE
      auto.efa <- TRUE
    } else if (model.type == "growth") {
      model.type <- "growth"
      int.ov.free <- FALSE
      int.lv.free <- TRUE
      auto.fix.first <- !std.lv
      auto.fix.single <- TRUE
      auto.var <- TRUE
      auto.cov.lv.x <- TRUE
      auto.cov.y <- TRUE
      auto.th <- TRUE
      auto.delta <- TRUE
      auto.efa <- TRUE
    }
  }

  # check for meanstructure
  if (any(flat$op == "~1")) {
    meanstructure <- TRUE
  }

  # check for block identifiers in the syntax (op = ":")
  n.block.flat <- length(which(flat$op == ":"))
  # this is NOT the number of blocks (eg group 1: level 1: -> 1 block)

  # for each non-empty `block' in n.block.flat, produce a USER
  if (n.block.flat > 0L) {
    # make sure flat is a data.frame
    flat <- as.data.frame(flat, stringsAsFactors = FALSE)

    # what are the block lhs labels?
    blocks.lhs.all <- tolower(flat$lhs[flat$op == ":"])
    tmp.block.lhs <- unique(blocks.lhs.all)

    # if we have group and level, check that group comes first!
    if ("group" %in% tmp.block.lhs && "level" %in% tmp.block.lhs) {
      group.idx <- which(tmp.block.lhs == "group")
      level.idx <- which(tmp.block.lhs == "level")
      if (group.idx > level.idx) {
        lav_msg_stop(gettext(
          "levels must be nested within groups (not the other way around)."))
      }
    }

    # block op == ":" indices
    block.op.idx <- which(flat$op == ":")

    # check for wrong spelled 'group' lhs
    if (length(grep("group", tmp.block.lhs)) > 1L) {
      lav_msg_warn(gettext("ambiguous block identifiers for group:"),
        lav_msg_view(tmp.block.lhs[grep("group", tmp.block.lhs)], "none"))
    }

    # no empty :rhs fields allowed!
    if (any(nchar(flat$rhs[block.op.idx]) == 0L)) {
      empty.idx <- nchar(flat$rhs[block.op.idx]) == 0L
      txt <- paste(flat$lhs[block.op.idx][empty.idx], ":")
      lav_msg_stop(gettext(
        "syntax contains block identifiers with missing numbers/labels: "), txt)
    }

    # check for ngroups (ngroups is based on the data!)
    if ("group" %in% tmp.block.lhs) {
      # how many group blocks?
      group.block.idx <- flat$op == ":" & flat$lhs == "group"
      n.group.flat <- length(unique(flat$rhs[group.block.idx]))

      if (n.group.flat > 0L && n.group.flat != ngroups) {
        lav_msg_stop(gettextf(
          "syntax defines %1$s groups; data (or argument ngroups)
          suggests %2$s groups", n.group.flat, ngroups))
      }
    }

    # figure out how many 'blocks' we have, and store indices/block.labels
    tmp.block.rhs <- rep("0", length(tmp.block.lhs))
    block.id <- 0L
    block.info <- vector("list", length = n.block.flat) # too large
    block.op.idx1 <- c(block.op.idx, nrow(flat) + 1L) # add addition row
    for (block.op in seq_len(n.block.flat)) {
      # fill block.rhs value(s)
      block.lhs <- flat$lhs[block.op.idx1[block.op]]
      block.rhs <- flat$rhs[block.op.idx1[block.op]]
      tmp.block.rhs[which(block.lhs == tmp.block.lhs)] <- block.rhs

      # another block identifier?
      if (block.op.idx1[block.op + 1L] - block.op.idx1[block.op] == 1L) {
        next
      }

      # we have a 'block'
      block.id <- block.id + 1L

      # select flat rows for this block
      tmp.idx <- seq.int(
        block.op.idx1[block.op] + 1L,
        block.op.idx1[block.op + 1L] - 1L
      )

      # store info in block.info
      block.info[[block.id]] <- list(
        lhs = tmp.block.lhs, # always the same
        rhs = tmp.block.rhs, # for this block
        idx = tmp.idx
      )
    }
    block.info <- block.info[seq_len(block.id)]

    # new in 0.6-12
    # check for blocks with the same block.rhs combination
    # (perhaps added later?)
    # - merge the indices
    # - remove the duplicated blocks
    block.labels <- sapply(lapply(block.info, "[[", "rhs"),
      paste,
      collapse = "."
    )
    nblocks <- length(unique(block.labels))
    if (nblocks < length(block.labels)) {
      # it would appear we have duplicated block.labels -> merge
      dup.idx <- which(duplicated(block.labels))
      for (i in seq_along(dup.idx)) {
        this.dup.idx <- dup.idx[i]
        orig.idx <- which(block.labels == block.labels[this.dup.idx])[1]
        block.info[[orig.idx]]$idx <- c(
          block.info[[orig.idx]]$idx,
          block.info[[this.dup.idx]]$idx
        )
      }
      block.info <- block.info[-dup.idx]
    }

    # split the flat data.frame per `block', create tmp.list
    # for each `block', and rbind them together, adding block columns
    for (block in seq_len(nblocks)) {
      tmp.block.rhs <- block.info[[block]]$rhs
      block.lhs <- block.info[[block]]$lhs[length(tmp.block.lhs)] # last one
      block.idx <- block.info[[block]]$idx

      flat.block <- flat[block.idx, ]
      # rm 'block' column (if any) in flat.block
      flat.block$block <- NULL

      # new in 0.6-7: check for random slopes, add them here
      if (block.lhs == "level" &&
        block > 1L && # FIXME: multigroup, multilevel
        !is.null(flat$rv) &&
        any(nchar(flat$rv) > 0L)) {
        lv.names.rv <- unique(flat$rv[nchar(flat$rv) > 0L])
        for (i in seq_along(lv.names.rv)) {
          # add phantom latent variable
          tmp <- flat.block[1, ]
          tmp$lhs <- lv.names.rv[i]
          tmp$op <- "=~"
          tmp$rhs <- lv.names.rv[i]
          tmp$mod.idx <- max(flat$mod.idx) + i
          tmp$fixed <- "0"
          tmp$start <- ""
          tmp$lower <- ""
          tmp$upper <- ""
          tmp$label <- ""
          tmp$prior <- ""
          tmp$efa <- ""
          tmp$rv <- lv.names.rv[i]
          flat.block <- rbind(flat.block, tmp, deparse.level = 0L)
          tmp.mod <- c(tmp.mod, list(list(fixed = 0)))
        }
      }

      # new in 0.6-8: if multilevel, use 'global' ov.names.x
      if (fixed.x && block.lhs == "level") {
        tmp.ov.names.x <- lav_partable_vnames(flat, "ov.x") # global
        ov.names.x.block <- lav_partable_vnames(flat.block, "ov.x")
        if (length(ov.names.x.block) > 0L) {
          idx <- which(!ov.names.x.block %in% tmp.ov.names.x)
          if (length(idx) > 0L) {
            # warn!
            lav_msg_warn(gettextf(
              "the variable(s) [%s] are exogenous at one level, but endogenous
              at another level. These variables will be treated as endogenous,
              and their variances/intercepts will be freely estimated.
              To remove this warning, use fixed.x = FALSE.",
              lav_msg_view(ov.names.x.block[idx], "none")))
            ov.names.x.block <- ov.names.x.block[-idx]
          }
        }
      } else {
        ov.names.x.block <- NULL
      }

      # new in 0.6-12: if multilevel and conditional.x, make sure
      # that 'splitted' exogenous covariates become 'y' variables
      if (conditional.x && block.lhs == "level") {
        if (ngroups == 1L) {
          other.block.names <- lav_partable_vnames(flat, "ov",
            block = seq_len(nblocks)[-block]
          )
        } else {
          # TEST ME
          this.group <- ceiling(block / nlevels)
          blocks.within.group <- (this.group - 1L) * nlevels + seq_len(nlevels)
          other.block.names <- lav_partable_vnames(flat, "ov",
            block = blocks.within.group[-block]
          )
        }
        ov.names.x.block <- lav_partable_vnames(flat.block, "ov.x")
        if (length(ov.names.x.block) > 0L) {
          idx <- which(ov.names.x.block %in% other.block.names)
          if (length(idx) > 0L) {
            ov.names.x.block <- ov.names.x.block[-idx]
          }
        }
      } else {
        ov.names.x.block <- NULL
      }

      list.block <- lav_partable_flat(flat.block,
        blocks = tmp.block.lhs,
        block.id = block,
        meanstructure = meanstructure,
        int.ov.free = int.ov.free, int.lv.free = int.lv.free,
        orthogonal = orthogonal, orthogonal.y = orthogonal.y,
        orthogonal.x = orthogonal.x, orthogonal.efa = orthogonal.efa,
        std.lv = std.lv, correlation = correlation,
        conditional.x = conditional.x, fixed.x = fixed.x,
        parameterization = parameterization,
        auto.fix.first = auto.fix.first,
        auto.fix.single = auto.fix.single,
        auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
        auto.cov.y = auto.cov.y, auto.th = auto.th,
        auto.delta = auto.delta, auto.efa = auto.efa,
        varTable = varTable, group.equal = NULL,
        group.w.free = group.w.free, ngroups = 1L,
        nthresholds = nthresholds,
        ov.names.x.block = ov.names.x.block
      )
      list.block <- as.data.frame(list.block, stringsAsFactors = FALSE)

      # add block columns with current values in block.rhs
      for (b in seq_len(length(tmp.block.lhs))) {
        block.lhs <- tmp.block.lhs[b]
        block.rhs <- tmp.block.rhs[b]
        list.block[block.lhs] <- rep(block.rhs, length(list.block$lhs))
      }

      if (!exists("tmp.list")) {
        tmp.list <- list.block
      } else {
        list.block$id <- list.block$id + max(tmp.list$id)
        tmp.list <- rbind(tmp.list, list.block)
      }
    }
    tmp.list <- as.list(tmp.list)

    # convert block columns to integers if possible
    for (b in seq_len(length(tmp.block.lhs))) {
      block.lhs <- tmp.block.lhs[b]
      block.rhs <- tmp.block.rhs[b]
      tmp <- try(scan(
        text = tmp.list[[block.lhs]], what = integer(),
        quiet = TRUE
      ), silent = TRUE)
      if (inherits(tmp, "integer")) {
        tmp.list[[block.lhs]] <- tmp
      }
    }
  } else {
    tmp.list <- lav_partable_flat(flat,
      blocks = "group",
      meanstructure = meanstructure,
      int.ov.free = int.ov.free, int.lv.free = int.lv.free,
      orthogonal = orthogonal, orthogonal.y = orthogonal.y,
      orthogonal.x = orthogonal.x, orthogonal.efa = orthogonal.efa,
      std.lv = std.lv, correlation = correlation,
      conditional.x = conditional.x, fixed.x = fixed.x,
      parameterization = parameterization,
      auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
      auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
      auto.cov.y = auto.cov.y, auto.th = auto.th,
      auto.delta = auto.delta, auto.efa = auto.efa,
      varTable = varTable, group.equal = group.equal,
      group.w.free = group.w.free,
      ngroups = ngroups, nthresholds = nthresholds
    )
  }
  if (debug) {
    cat("[lavaan DEBUG]: parameter tmp.list without MODIFIERS:\n")
    print(as.data.frame(tmp.list, stringsAsFactors = FALSE))
  }

  # handle multilevel-specific constraints
  multilevel <- FALSE
  if (!is.null(tmp.list$level)) {
    nlevels <- lav_partable_nlevels(tmp.list)
    if (nlevels > 1L) {
      multilevel <- TRUE
    }
  }
  if (multilevel && any(tmp.list$op == "~1")) {
    # fix ov intercepts for all within ov that also appear at level 2
    # FIXME: not tested with > 2 levels
    ov.names <- lav_partable_vnames(tmp.list, "ov") ## all names
    level.values <- lav_partable_level_values(tmp.list)
    other.names <- tmp.list$lhs[tmp.list$op == "~1" &
      tmp.list$level %in% level.values[-1L] &
      tmp.list$lhs %in% ov.names]
    fix.names.idx <- which(tmp.list$op == "~1" &
      tmp.list$level %in% level.values[1L] &
      tmp.list$lhs %in% other.names)
    if (length(fix.names.idx) > 0L) {
      tmp.list$free[fix.names.idx] <- 0L
      tmp.list$ustart[fix.names.idx] <- 0
    }
  }
  if (multilevel && any(tmp.list$op == "|")) {
    # fix ALL thresholds at level 1
    level.values <- lav_partable_level_values(tmp.list)
    th.idx <- which(tmp.list$op == "|" &
      tmp.list$level %in% level.values[1L])
    tmp.list$free[th.idx] <- 0L
    tmp.list$ustart[th.idx] <- 0

    # fix ALL scaling parmaters at higher levels
    scale.idx <- which(tmp.list$op == "~*~" &
      tmp.list$level %in% level.values[-1L])
    tmp.list$free[scale.idx] <- 0L
    tmp.list$ustart[scale.idx] <- 1
  }

  # apply user-specified modifiers
  warn.about.single.label <- FALSE
  if (length(tmp.mod)) {
    for (el in seq_along(tmp.mod)) {
      idx <- which(tmp.list$mod.idx == el) # for each group

      # 0.5-21: check if idx exists
      # perhaps the corresponding element was duplicated, and removed
      if (length(idx) == 0L) {
        next
      }

      tmp.mod.fixed <- tmp.mod[[el]]$fixed
      tmp.mod.start <- tmp.mod[[el]]$start
      tmp.mod.lower <- tmp.mod[[el]]$lower
      tmp.mod.upper <- tmp.mod[[el]]$upper
      tmp.mod.label <- tmp.mod[[el]]$label
      tmp.mod.prior <- tmp.mod[[el]]$prior
      tmp.mod.efa <- tmp.mod[[el]]$efa
      tmp.mod.rv <- tmp.mod[[el]]$rv

      # check for single argument if multiple groups
      if (ngroups > 1L && length(idx) > 1L) {
        # Ok, this is not very consistent:
        # A) here we force same behavior across groups
        if (length(tmp.mod.fixed) == 1L) {
          tmp.mod.fixed <- rep(tmp.mod.fixed, ngroups)
        }
        if (length(tmp.mod.start) == 1L) {
          tmp.mod.start <- rep(tmp.mod.start, ngroups)
        }
        if (length(tmp.mod.lower) == 1L) {
          tmp.mod.lower <- rep(tmp.mod.lower, ngroups)
        }
        if (length(tmp.mod.upper) == 1L) {
          tmp.mod.upper <- rep(tmp.mod.upper, ngroups)
        }
        if (length(tmp.mod.prior) == 1L) {
          tmp.mod.prior <- rep(tmp.mod.prior, ngroups)
        }
        if (length(tmp.mod.efa) == 1L) {
          tmp.mod.efa <- rep(tmp.mod.efa, ngroups)
        }
        if (length(tmp.mod.rv) == 1L) {
          tmp.mod.rv <- rep(tmp.mod.rv, ngroups)
        }

        # new in 0.6-7 (proposal):
        # - always recycle modifiers, including labels
        # - if ngroups > 1 AND group.label= is empty, produce a warning
        #   (as this is a break from < 0.6-6)
        if (length(tmp.mod.label) == 1L) {
          tmp.mod.label <- rep(tmp.mod.label, ngroups)
          if (is.null(group.equal) || length(group.equal) == 0L) {
            warn.about.single.label <- TRUE
          }
        }

        # < 0.6-7 code:
        # B) here we do NOT! otherwise, it would imply an equality
        #                    constraint...
        #    except if group.equal="loadings"!
        # if(length(tmp.mod.label) == 1L) {
        #    if("loadings" %in% group.equal ||
        #       "composite.loadings" %in% group.equal) {
        #        tmp.mod.label <- rep(tmp.mod.label, ngroups)
        #    } else {
        #        tmp.mod.label <- c(tmp.mod.label, rep("", (ngroups-1L)) )
        #    }
        # }
      }

      # check for wrong number of arguments if multiple groups
      nidx <- length(idx)
      if ((!is.null(tmp.mod.fixed) && nidx != length(tmp.mod.fixed)) ||
        (!is.null(tmp.mod.start) && nidx != length(tmp.mod.start)) ||
        (!is.null(tmp.mod.lower) && nidx != length(tmp.mod.lower)) ||
        (!is.null(tmp.mod.upper) && nidx != length(tmp.mod.upper)) ||
        (!is.null(tmp.mod.prior) && nidx != length(tmp.mod.prior)) ||
        (!is.null(tmp.mod.efa) && nidx != length(tmp.mod.efa)) ||
        (!is.null(tmp.mod.rv) && nidx != length(tmp.mod.rv)) ||
        (!is.null(tmp.mod.label) && nidx != length(tmp.mod.label))) {
        el.idx <- which(tmp.list$mod.idx == el)[1L]
        lav_msg_stop(gettextf(
          "wrong number of arguments in modifier (%s) of element",
          lav_msg_view(tmp.mod.label, "none")),
          tmp.list$lhs[el.idx], tmp.list$op[el.idx], tmp.list$rhs[el.idx]
        )
      }

      # apply modifiers
      if (!is.null(tmp.mod.fixed)) {
        # two options: constant or NA
        na.idx <- which(is.na(tmp.mod.fixed))
        not.na.idx <- which(!is.na(tmp.mod.fixed))

        # constant
        tmp.list$ustart[idx][not.na.idx] <- tmp.mod.fixed[not.na.idx]
        tmp.list$free[idx][not.na.idx] <- 0L

        # NA* modifier
        tmp.list$free[idx][na.idx] <- 1L # eg factor loading
        tmp.list$ustart[idx][na.idx] <- as.numeric(NA)
      }
      if (!is.null(tmp.mod.start)) {
        tmp.list$ustart[idx] <- tmp.mod.start
      }
      if (!is.null(tmp.mod.prior)) {
        # do we already have a `prior' column? if not, create one
        if (is.null(tmp.list$prior)) {
          tmp.list$prior <- character(length(tmp.list$lhs))
        }
        tmp.list$prior[idx] <- tmp.mod.prior
      }
      if (!is.null(tmp.mod.efa)) {
        # do we already have a `efa' column? if not, create one
        if (is.null(tmp.list$efa)) {
          tmp.list$efa <- character(length(tmp.list$lhs))
        }
        tmp.list$efa[idx] <- tmp.mod.efa
      }
      if (!is.null(tmp.mod.rv)) {
        # do we already have a `rv' column? if not, create one
        if (is.null(tmp.list$rv)) {
          tmp.list$rv <- character(length(tmp.list$lhs))
        }
        tmp.list$rv[idx] <- tmp.mod.rv

        tmp.list$free[idx] <- 0L
        tmp.list$ustart[idx] <- as.numeric(NA) #
      }
      if (!is.null(tmp.mod.lower)) {
        # do we already have a `lower' column? if not, create one
        if (is.null(tmp.list$lower)) {
          tmp.list$lower <- rep(-Inf, length(tmp.list$lhs))
        }
        tmp.list$lower[idx] <- as.numeric(tmp.mod.lower)
      }
      if (!is.null(tmp.mod.upper)) {
        # do we already have a `upper' column? if not, create one
        if (is.null(tmp.list$upper)) {
          tmp.list$upper <- rep(Inf, length(tmp.list$lhs))
        }
        tmp.list$upper[idx] <- as.numeric(tmp.mod.upper)
      }
      if (!is.null(tmp.mod.label)) {
        tmp.list$label[idx] <- tmp.mod.label
      }
    }
  }
  # remove mod.idx column
  tmp.list$mod.idx <- NULL

  # warning about single label in multiple group setting?
  if (warn.about.single.label) {
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
  if (!is.null(tmp.list$lower)) {
    fixed.idx <- which(tmp.list$free == 0L)
    if (length(fixed.idx) > 0L) {
      tmp.list$lower[fixed.idx] <- tmp.list$ustart[fixed.idx]
    }
  }
  if (!is.null(tmp.list$upper)) {
    fixed.idx <- which(tmp.list$free == 0L)
    if (length(fixed.idx) > 0L) {
      tmp.list$upper[fixed.idx] <- tmp.list$ustart[fixed.idx]
    }
  }

  # if rv column is present, add rv.names to ALL rows where they are used
  if (!is.null(tmp.list$rv)) {
    rv.names <- unique(tmp.list$rv[nchar(tmp.list$rv) > 0L])
    for (i in seq_len(length(rv.names))) {
      lhs.idx <- which(tmp.list$lhs == rv.names[i] &
        tmp.list$op == "=~")
      tmp.list$rv[lhs.idx] <- rv.names[i]
    }
  }

  if (debug) {
    cat("[lavaan DEBUG]: parameter tmp.list with MODIFIERS:\n")
    print(as.data.frame(tmp.list, stringsAsFactors = FALSE))
  }

  # get 'virtual' parameter labels
  if (n.block.flat > 1L) {
    blocks <- tmp.block.lhs
  } else {
    blocks <- "group"
  }
  label <- lav_partable_labels(
    partable = tmp.list,
    blocks = blocks,
    group.equal = group.equal,
    group.partial = group.partial
  )

  if (debug) {
    cat("[lavaan DEBUG]: parameter tmp.list with LABELS:\n")
    tmp <- tmp.list
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
  if (auto.efa && !is.null(tmp.list$efa)) {
    tmp.list <- lav_partable_efa_constraints(
      LIST = tmp.list,
      orthogonal.efa = orthogonal.efa,
      group.equal = group.equal
    )
  } # auto.efa

  # handle user-specified equality constraints
  # lavaan 0.6-11:
  #     two settings:
  #       1) simple equality constraints ONLY -> back to basics: only
  #          duplicate 'free' numbers; no longer explicit == rows with plabels
  #       2) mixture of simple and other (explicit) constraints
  #          treat them together as we did in <0.6-11
  tmp.list$plabel <- paste(".p", tmp.list$id, ".", sep = "")
  eq.labels <- unique(label[duplicated(label)])
  eq.id <- integer(length(tmp.list$lhs))
  for (eq.label in eq.labels) {
    tmp.con.idx <- length(tmp.con)
    all.idx <- which(label == eq.label) # all same-label parameters
    ref.idx <- all.idx[1L] # the first one only
    other.idx <- all.idx[-1L] # the others
    eq.id[all.idx] <- ref.idx

    # new in 0.6-6: make sure lower/upper constraints are equal too
    if (!is.null(tmp.list$lower) &&
      length(unique(tmp.list$lower[all.idx])) > 0L) {
      non.inf <- which(is.finite(tmp.list$lower[all.idx]))
      if (length(non.inf) > 0L) {
        smallest.val <- min(tmp.list$lower[all.idx][non.inf])
        tmp.list$lower[all.idx] <- smallest.val
      }
    }
    if (!is.null(tmp.list$upper) &&
      length(unique(tmp.list$upper[all.idx])) > 0L) {
      non.inf <- which(is.finite(tmp.list$upper[all.idx]))
      if (length(non.inf) > 0L) {
        largest.val <- max(tmp.list$upper[all.idx][non.inf])
        tmp.list$upper[all.idx] <- largest.val
      }
    }

    # two possibilities:
    #   1. all.idx contains a fixed parameter: in this case,
    #      we fix them all (hopefully to the same value)
    #   2. all.idx contains only free parameters

    # 1. all.idx contains a fixed parameter
    if (any(tmp.list$free[all.idx] == 0L)) {
      # which one is fixed?
      fixed.all <- all.idx[tmp.list$free[all.idx] == 0L]
      # only pick the first
      fixed.idx <- fixed.all[1]

      # sanity check: are all ustart values equal?
      ustart1 <- tmp.list$ustart[fixed.idx]
      if (!all(ustart1 == tmp.list$ustart[fixed.all])) {
        lav_msg_warn(gettext(
          "equality constraints involve fixed parameters with different values;
          only the first one will be used"))
      }

      # make them all fixed
      tmp.list$ustart[all.idx] <- tmp.list$ustart[fixed.idx]
      tmp.list$free[all.idx] <- 0L # not free anymore, since it must
      # be equal to the 'fixed' parameter
      # (Note: Mplus ignores this)
      eq.id[all.idx] <- 0L # remove from eq.id list

      # new in 0.6-8 (for efa + user-specified eq constraints)
      if (any(tmp.list$user[all.idx] %in% c(7L, 77L))) {
        # if involved in an efa block, store in tmp.con anyway
        # we may need it for the rotated solution
        for (o in other.idx) {
          tmp.con.idx <- tmp.con.idx + 1L
          tmp.con[[tmp.con.idx]] <- list(
            op = "==",
            lhs = tmp.list$plabel[ref.idx],
            rhs = tmp.list$plabel[o],
            user = 2L
          )
        }
      }
    } else {
      # 2. all.idx contains only free parameters
      # old system:
      # - add tmp.con entry
      # - in 0.6-11: only if tmp.con is not empty
      if (!ceq.simple) {
        for (o in other.idx) {
          tmp.con.idx <- tmp.con.idx + 1L
          tmp.con[[tmp.con.idx]] <- list(
            op = "==",
            lhs = tmp.list$plabel[ref.idx],
            rhs = tmp.list$plabel[o],
            user = 2L
          )
        }
      } else {
        # new system:
        # - set $free elements to zero, and later to ref id
        tmp.list$free[other.idx] <- 0L # all but the first are non-free
        # but will get a duplicated number
      }

      # just to trick semTools, also add something in the label
      # colum, *if* it is empty
      # update: 0.6-11 we keep this, because it shows the plabels
      #         when eg group.equal = "loadings"
      for (i in all.idx) {
        if (nchar(tmp.list$label[i]) == 0L) {
          tmp.list$label[i] <- tmp.list$plabel[ref.idx]
        }
      }
    } # all free
  } # eq in eq.labels
  if (debug) {
    print(tmp.con)
  }



  # handle constraints (if any) (NOT per group, but overall - 0.4-11)
  if (length(tmp.con) > 0L) {
    n.con <- length(tmp.con)
    tmp.idx <- length(tmp.list$id) + seq_len(n.con)
    # grow tmp.list with length(tmp.con) extra rows
    tmp.list <- lapply(tmp.list, function(x) {
      if (is.character(x)) {
        c(x, rep("", n.con))
      } else {
        c(x, rep(NA, n.con))
      }
    })

    # fill in some columns
    tmp.list$id[tmp.idx] <- tmp.idx
    tmp.list$lhs[tmp.idx] <- unlist(lapply(tmp.con, "[[", "lhs"))
    tmp.list$op[tmp.idx] <- unlist(lapply(tmp.con, "[[", "op"))
    tmp.list$rhs[tmp.idx] <- unlist(lapply(tmp.con, "[[", "rhs"))
    tmp.list$user[tmp.idx] <- unlist(lapply(tmp.con, "[[", "user"))

    # zero is nicer?
    tmp.list$free[tmp.idx] <- rep(0L, n.con)
    tmp.list$exo[tmp.idx] <- rep(0L, n.con)
    tmp.list$block[tmp.idx] <- rep(0L, n.con)

    if (!is.null(tmp.list$group)) {
      if (is.character(tmp.list$group)) {
        tmp.list$group[tmp.idx] <- rep("", n.con)
      } else {
        tmp.list$group[tmp.idx] <- rep(0L, n.con)
      }
    }
    if (!is.null(tmp.list$level)) {
      if (is.character(tmp.list$level)) {
        tmp.list$level[tmp.idx] <- rep("", n.con)
      } else {
        tmp.list$level[tmp.idx] <- rep(0L, n.con)
      }
    }
    if (!is.null(tmp.list$class)) {
      if (is.character(tmp.list$class)) {
        tmp.list$class[tmp.idx] <- rep("", n.con)
      } else {
        tmp.list$class[tmp.idx] <- rep(0L, n.con)
      }
    }
  }

  # put lhs of := elements in label column
  def.idx <- which(tmp.list$op == ":=")
  tmp.list$label[def.idx] <- tmp.list$lhs[def.idx]


  # handle effect.coding related equality constraints
  if (is.logical(effect.coding) && effect.coding) {
    effect.coding <- c("loadings", "intercepts")
  } else if (!is.character(effect.coding)) {
    lav_msg_stop(gettext("effect.coding argument must be a character string"))
  }
  if (any(c("loadings", "intercepts") %in% effect.coding)) {
    tmp <- list()
    # for each block
    nblocks <- lav_partable_nblocks(tmp.list)
    for (b in seq_len(nblocks)) {
      # lv's for this block/set
      lv.names <- unique(tmp.list$lhs[tmp.list$op == "=~" &
        tmp.list$block == b])

      if (length(lv.names) == 0L) {
        next
      }

      int.plabel <- character(0L)
      for (lv in lv.names) {
        # ind.names
        ind.names <- tmp.list$rhs[tmp.list$op == "=~" &
          tmp.list$block == b &
          tmp.list$lhs == lv]

        if ("loadings" %in% effect.coding) {
          # factor loadings indicators of this lv
          loadings.idx <- which(tmp.list$op == "=~" &
            tmp.list$block == b &
            tmp.list$rhs %in% ind.names &
            tmp.list$lhs == lv)

          # all free?
          if (length(loadings.idx) > 0L &&
            all(tmp.list$free[loadings.idx] > 0L)) {
            # add eq constraint
            plabel <- tmp.list$plabel[loadings.idx]

            # Note:  we write them as
            # .p1. == 3 - .p2. - .p3.
            # instead of
            # 3 ==  .p1.+.p2.+.p3.
            # as this makes it easier to translate things to
            # JAGS/stan

            tmp.lhs <- plabel[1]
            if (length(loadings.idx) > 1L) {
              tmp.rhs <- paste(length(loadings.idx), "-",
                paste(plabel[-1], collapse = "-"),
                sep = ""
              )
            } else {
              tmp.rhs <- length(loadings.idx)
            }

            tmp$lhs <- c(tmp$lhs, tmp.lhs)
            tmp$op <- c(tmp$op, "==")
            tmp$rhs <- c(tmp$rhs, tmp.rhs)
            tmp$block <- c(tmp$block, 0L)
            tmp$user <- c(tmp$user, 2L)
            tmp$ustart <- c(tmp$ustart, as.numeric(NA))
          }
        } # loadings

        if ("intercepts" %in% effect.coding) {
          # intercepts for indicators of this lv
          intercepts.idx <- which(tmp.list$op == "~1" &
            tmp.list$block == b &
            tmp.list$lhs %in% ind.names)

          # all free?
          if (length(intercepts.idx) > 0L &&
            all(tmp.list$free[intercepts.idx] > 0L)) {
            # 1) add eq constraint
            plabel <- tmp.list$plabel[intercepts.idx]

            tmp.lhs <- plabel[1]
            if (length(intercepts.idx) > 1L) {
              tmp.rhs <- paste("0-",
                paste(plabel[-1], collapse = "-"),
                sep = ""
              )
            } else {
              tmp.rhs <- 0L
            }

            tmp$lhs <- c(tmp$lhs, tmp.lhs)
            tmp$op <- c(tmp$op, "==")
            tmp$rhs <- c(tmp$rhs, tmp.rhs)
            tmp$block <- c(tmp$block, 0L)
            tmp$user <- c(tmp$user, 2L)
            tmp$ustart <- c(tmp$ustart, as.numeric(NA))

            # 2) release latent mean
            lv.int.idx <- which(tmp.list$op == "~1" &
              tmp.list$block == b &
              tmp.list$lhs == lv)
            # free only if automatically added
            if (length(lv.int.idx) > 0L &&
              tmp.list$user[lv.int.idx] == 0L) {
              tmp.list$free[lv.int.idx] <- 1L
            }
          }
        } # intercepts
      } # lv
    } # blocks

    tmp.list <- lav_partable_merge(tmp.list, tmp)
  }

  # marker.int.zero
  if (meanstructure && marker.int.zero) {
    # for each block
    nblocks <- lav_partable_nblocks(tmp.list)
    for (b in seq_len(nblocks)) {
      # lv's for this block/set
      lv.names <- lav_partable_vnames(tmp.list,
        type = "lv.regular",
        block = b
      )
      lv.marker <- lav_partable_vnames(tmp.list,
        type = "lv.regular",
        block = b
      )

      if (length(lv.names) == 0L) {
        next
      }

      # markers for this block
      lv.marker <- lav_partable_vnames(tmp.list,
        type = "lv.marker",
        block = b
      )

      # fix marker intercepts to zero
      marker.idx <- which(tmp.list$op == "~1" &
        tmp.list$lhs %in% lv.marker & tmp.list$block == b &
        tmp.list$user == 0L)
      tmp.list$free[marker.idx] <- 0L
      tmp.list$ustart[marker.idx] <- 0

      # free latent means
      lv.idx <- which(tmp.list$op == "~1" &
        tmp.list$lhs %in% lv.names & tmp.list$block == b &
        tmp.list$user == 0L)
      tmp.list$free[lv.idx] <- 1L
      tmp.list$ustart[lv.idx] <- as.numeric(NA)
    } # block
  }


  # mg.lv.variances
  if (ngroups > 1L && "mg.lv.variances" %in% effect.coding) {
    tmp <- list()

    # do not include 'EFA' lv's
    if (!is.null(tmp.list$efa)) {
      lv.names <- unique(tmp.list$lhs[tmp.list$op == "=~" &
        !nchar(tmp.list$efa) > 0L])
    } else {
      lv.names <- unique(tmp.list$lhs[tmp.list$op == "=~"])
    }
    group.values <- lav_partable_group_values(tmp.list)

    for (lv in lv.names) {
      # factor variances
      lv.var.idx <- which(tmp.list$op == "~~" &
        tmp.list$lhs == lv &
        tmp.list$rhs == tmp.list$lhs &
        tmp.list$lhs == lv)

      # all free (but the first?)
      if (length(lv.var.idx) > 0L &&
        all(tmp.list$free[lv.var.idx][-1] > 0L)) {
        # 1) add eq constraint
        plabel <- tmp.list$plabel[lv.var.idx]

        tmp.lhs <- plabel[1]
        if (length(lv.var.idx) > 1L) {
          tmp.rhs <- paste(length(lv.var.idx), "-",
            paste(plabel[-1], collapse = "-"),
            sep = ""
          )
        } else {
          tmp.rhs <- length(lv.var.idx)
        }

        tmp$lhs <- c(tmp$lhs, tmp.lhs)
        tmp$op <- c(tmp$op, "==")
        tmp$rhs <- c(tmp$rhs, tmp.rhs)
        tmp$block <- c(tmp$block, 0L)
        tmp$user <- c(tmp$user, 2L)
        tmp$ustart <- c(tmp$ustart, as.numeric(NA))

        # 2) free lv variances first group
        lv.var.g1.idx <- which(tmp.list$op == "~~" &
          tmp.list$group == group.values[1] &
          tmp.list$lhs == lv &
          tmp.list$rhs == tmp.list$lhs &
          tmp.list$lhs == lv)
        # free only if automatically added
        if (length(lv.var.g1.idx) > 0L &&
          tmp.list$user[lv.var.g1.idx] == 0L) {
          tmp.list$free[lv.var.g1.idx] <- 1L
        }
      }
    } # lv

    tmp.list <- lav_partable_merge(tmp.list, tmp)
  }

  # mg.lv.efa.variances
  if (ngroups > 1L && "mg.lv.efa.variances" %in% effect.coding) {
    tmp <- list()

    # only 'EFA' lv's
    if (!is.null(tmp.list$efa)) {
      lv.names <- unique(tmp.list$lhs[tmp.list$op == "=~" &
        nchar(tmp.list$efa) > 0L])
    } else {
      lv.names <- character(0L)
    }
    group.values <- lav_partable_group_values(tmp.list)

    for (lv in lv.names) {
      # factor variances
      lv.var.idx <- which(tmp.list$op == "~~" &
        tmp.list$lhs == lv &
        tmp.list$rhs == tmp.list$lhs &
        tmp.list$lhs == lv)

      # all free (but the first?)
      if (length(lv.var.idx) > 0L &&
        all(tmp.list$free[lv.var.idx][-1] > 0L)) {
        # 1) add eq constraint
        plabel <- tmp.list$plabel[lv.var.idx]

        tmp.lhs <- plabel[1]
        if (length(lv.var.idx) > 1L) {
          tmp.rhs <- paste(length(lv.var.idx), "-",
            paste(plabel[-1], collapse = "-"),
            sep = ""
          )
        } else {
          tmp.rhs <- length(lv.var.idx)
        }

        tmp$lhs <- c(tmp$lhs, tmp.lhs)
        tmp$op <- c(tmp$op, "==")
        tmp$rhs <- c(tmp$rhs, tmp.rhs)
        tmp$block <- c(tmp$block, 0L)
        tmp$user <- c(tmp$user, 2L)
        tmp$ustart <- c(tmp$ustart, as.numeric(NA))

        # 2) free lv variances first group
        lv.var.g1.idx <- which(tmp.list$op == "~~" &
          tmp.list$group == group.values[1] &
          tmp.list$lhs == lv &
          tmp.list$rhs == tmp.list$lhs &
          tmp.list$lhs == lv)
        # free only if automatically added
        if (length(lv.var.g1.idx) > 0L &&
          tmp.list$user[lv.var.g1.idx] == 0L) {
          tmp.list$free[lv.var.g1.idx] <- 1L
        }
      }
    } # lv

    tmp.list <- lav_partable_merge(tmp.list, tmp)
  }


  # count free parameters
  idx.free <- which(tmp.list$free > 0L)
  tmp.list$free[idx.free] <- seq_along(idx.free)

  # new in 0.6-11: add free counter to this element (as in < 0.5-18)
  # unless we have other constraints
  if (ceq.simple) {
    idx.equal <- which(eq.id > 0)
    tmp.list$free[idx.equal] <- tmp.list$free[eq.id[idx.equal]]
  }

  # new in 0.6-14: add 'da' entries to reflect data-based order of ov's
  # now via attribute "ovda"
  attr(tmp.list, "ovda") <- ov.names.data

  # backwards compatibility...
  if (!is.null(tmp.list$unco)) {
    tmp.list$unco[idx.free] <- seq_along(sum(tmp.list$free > 0L))
  }

  if (debug) {
    cat("[lavaan DEBUG] lavParTable\n")
    print(as.data.frame(tmp.list))
  }

  # data.frame?
  if (as.data.frame.) {
    tmp.list <- as.data.frame(tmp.list, stringsAsFactors = FALSE)
    attr(tmp.list, "ovda") <- ov.names.data
  } else {
    tmp.list <- lav_partable_set_cache(tmp.list) # add cached "pta" data
  }

  tmp.list
}
