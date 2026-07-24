# support for data-defined ("moderated") parameters using the ":~" operator
#
# a ":~" definition redefines an EXISTING (labeled) free parameter as a
# per-observation function of (new) component parameters and observed
# ("definition") variables, e.g.
#
#     f1 =~ 1*y1 + L2*y2
#     L2 :~ L2_0 + L2_1 * data(z)
#
# rules:
#  - the left-hand side must be the label of an existing free parameter
#    (the "host" parameter); the host row(s) become fixed (free = 0), and
#    are flagged in a (new) $dv column
#  - every definition variable must be wrapped in data(); this makes the
#    classification of the right-hand side tokens context-free:
#      * data(z)                -> definition variable (a column of the data)
#      * an existing label      -> reference to that (free) parameter
#      * another ":~" lhs       -> reference to that definition
#      * a known function       -> function call
#      * anything else          -> a NEW (free) component parameter
#  - component parameters are added to the parameter table as free scalar
#    rows with op "dp" (definition parameter)
#
# note: lavaan itself cannot (yet) estimate these models: they require a
# casewise likelihood. the parameter table + the dv.function in the @Model
# slot provide the complete model specification for external backends.
#
# YR/Claude 24 July 2026: initial version

# collect the variable names inside data(...) calls in expression 'x'
lav_dv_expr_datavars <- function(x) {
  if (is.call(x)) {
    if (identical(x[[1L]], as.name("data"))) {
      if (length(x) != 2L || !is.name(x[[2L]])) {
        lav_msg_stop(gettext(
          "data() in a ':~' expression must contain a single variable name."))
      }
      return(as.character(x[[2L]]))
    }
    if (length(x) < 2L) {
      return(character(0L))
    }
    return(unique(unlist(lapply(as.list(x)[-1L], lav_dv_expr_datavars))))
  }
  character(0L)
}

# replace every data(z) call in expression 'x' by the (mangled) symbol
# `.data.z.' (used when building the body of dv.function; the mangling
# avoids name clashes with parameter labels)
lav_dv_expr_data2sym <- function(x) {
  if (is.call(x)) {
    if (identical(x[[1L]], as.name("data"))) {
      return(as.name(paste0(".data.", as.character(x[[2L]]), ".")))
    }
    if (length(x) >= 2L) {
      for (j in seq.int(2L, length(x))) {
        x[[j]] <- lav_dv_expr_data2sym(x[[j]])
      }
    }
  }
  x
}

# replace every data(z) call in expression 'x' by the (numeric) value
# in values[["z"]] (used to compute representative starting values)
lav_dv_expr_sub_data <- function(x, values) {
  if (is.call(x)) {
    if (identical(x[[1L]], as.name("data"))) {
      return(values[[as.character(x[[2L]])]])
    }
    if (length(x) >= 2L) {
      for (j in seq.int(2L, length(x))) {
        x[[j]] <- lav_dv_expr_sub_data(x[[j]], values)
      }
    }
  }
  x
}

# replace every data(...) call in expression 'x' by a sentinel symbol,
# so that all.vars()/all.names() no longer see the definition variables
lav_dv_expr_strip_data <- function(x) {
  if (is.call(x)) {
    if (identical(x[[1L]], as.name("data"))) {
      return(as.name(".dv.data."))
    }
    if (length(x) >= 2L) {
      for (j in seq.int(2L, length(x))) {
        x[[j]] <- lav_dv_expr_strip_data(x[[j]])
      }
    }
  }
  x
}

# process the ":~" rows of a (still list-based) parameter table:
#  - assign labels to the ":~" rows
#  - classify the right-hand side tokens of each definition
#  - flag the host rows (free = 0, $dv column)
#  - append free "dp" rows for the (new) component parameters
# the free counter is (re)assigned later in lav_model_pt()
lav_pt_dv <- function(partable, var_table = NULL, data_names = NULL) {
  dv_idx <- which(partable$op == ":~")
  if (length(dv_idx) == 0L) {
    return(partable)
  }

  con_ops <- c("==", "<", ">", ":=", ":~")
  dv_lhs <- partable$lhs[dv_idx]

  # each parameter can have only one definition
  dup_idx <- which(duplicated(dv_lhs))
  if (length(dup_idx) > 0L) {
    lav_msg_stop(gettextf(
      "parameter(s) defined (using the :~ operator) more than once: %s.",
      lav_msg_view(unique(dv_lhs[dup_idx]), "none")))
  }

  # put lhs of ":~" rows in the label column (as we do for ":=")
  partable$label[dv_idx] <- dv_lhs

  # model rows and their labels/variable names
  model_row_idx <- which(!partable$op %in% con_ops)
  model_labels <- unique(partable$label[model_row_idx][
    nchar(partable$label[model_row_idx]) > 0L])
  model_vnames <- unique(c(partable$lhs[model_row_idx],
                           partable$rhs[model_row_idx]))
  model_vnames <- model_vnames[nchar(model_vnames) > 0L]
  def_lhs <- unique(partable$lhs[partable$op == ":="])

  # process each definition
  dv_data_vars <- vector("list", length(dv_idx))
  dv_refs <- vector("list", length(dv_idx)) # references to other ":~" lhs
  component_names <- character(0L)
  for (i in seq_along(dv_idx)) {
    this_lhs <- dv_lhs[i]
    this_rhs <- partable$rhs[dv_idx[i]]
    expr <- lav_pt_con_parse(this_rhs, gettext("parameter definition(s)"))
    if (length(expr) != 1L) {
      lav_msg_stop(gettextf(
        "right-hand side of `%s' contains more than one expression.",
        paste(this_lhs, ":~", this_rhs)))
    }
    expr <- expr[[1L]]

    # definition variables: data(...)
    data_vars <- lav_dv_expr_datavars(expr)
    if (length(data_vars) == 0L) {
      lav_msg_stop(gettextf(
        "the ':~' definition `%1$s' does not contain any data() reference;
        use the ':=' operator (or an equality constraint) for parameters
        that are a function of other parameters only.",
        paste(this_lhs, ":~", this_rhs)))
    }
    dv_data_vars[[i]] <- data_vars

    # remaining tokens (data() parts stripped)
    stripped <- lav_dv_expr_strip_data(expr)
    tokens <- setdiff(all.vars(stripped), ".dv.data.")

    # referencing a ':=' defined parameter is not allowed (for now)
    bad <- tokens[tokens %in% def_lhs]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "the ':~' definition `%1$s' refers to parameter(s) defined by the
        ':=' operator: %2$s; this is not supported.",
        paste(this_lhs, ":~", this_rhs), lav_msg_view(bad, "none")))
    }

    # references to other ':~' definitions (evaluated as their
    # per-observation values)
    dv_refs[[i]] <- tokens[tokens %in% dv_lhs]

    # references to existing (free) parameter labels
    ref_labels <- tokens[tokens %in% model_labels & !tokens %in% dv_lhs]
    for (lab in ref_labels) {
      lab_free <- partable$free[model_row_idx][
        partable$label[model_row_idx] == lab]
      if (all(lab_free == 0L)) {
        lav_msg_stop(gettextf(
          "the ':~' definition `%1$s' refers to the FIXED parameter `%2$s';
          only free parameters can be referenced.",
          paste(this_lhs, ":~", this_rhs), lab))
      }
    }

    # everything else is a NEW component parameter
    new_tokens <- setdiff(tokens, c(dv_lhs, ref_labels))

    # safety: a new component parameter should not have the same name as
    # an observed/latent variable in the model
    bad <- new_tokens[new_tokens %in% model_vnames]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "token(s) in the ':~' definition `%1$s' match the name of a
        model variable: %2$s; if you intended a definition variable, use
        data(); otherwise, please choose a different parameter label.",
        paste(this_lhs, ":~", this_rhs), lav_msg_view(bad, "none")))
    }

    # safety: warn if a new component parameter matches a column of the
    # data (most likely a forgotten data() wrapper)
    known_data_names <- unique(c(var_table$name, data_names))
    if (length(known_data_names) > 0L) {
      bad <- new_tokens[new_tokens %in% known_data_names]
      if (length(bad) > 0L) {
        lav_msg_warn(gettextf(
          "token(s) in the ':~' definition `%1$s' match a variable name
          in the data: %2$s; they are treated as NEW parameters -- if you
          intended the data column, use data().",
          paste(this_lhs, ":~", this_rhs), lav_msg_view(bad, "none")))
      }
    }

    component_names <- c(component_names,
                         new_tokens[!new_tokens %in% component_names])

    # check for unknown functions
    bad_fun <- lav_pt_con_undefined_fun(deparse(stripped, width.cutoff = 500L),
                                        environment())
    bad_fun <- setdiff(bad_fun, "data")
    if (length(bad_fun) > 0L) {
      lav_msg_stop(gettextf(
        "could not find function(s) used in the ':~' definition `%1$s': %2$s",
        paste(this_lhs, ":~", this_rhs), lav_msg_view(bad_fun, "none")))
    }
  }

  # check for cycles among cross-referenced definitions: order the
  # definitions topologically; in a cycle-free system, each definition
  # (in sorted order) only refers to definitions that came earlier
  # (note: on a cycle, lav_graph_order_adj_mat() silently returns the
  # original order, so we validate the sorted order explicitly)
  adj_mat <- matrix(0L, nrow = length(dv_idx), ncol = length(dv_idx))
  for (i in seq_along(dv_idx)) {
    adj_mat[match(dv_refs[[i]], dv_lhs), i] <- 1L
  }
  order_idx <- lav_graph_order_adj_mat(adj_mat, warn = FALSE)
  available <- character(0L)
  for (i in order_idx) {
    bad <- dv_refs[[i]][!dv_refs[[i]] %in% available]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "cyclic cross-reference(s) among ':~' definitions: `%1$s'
        refers to %2$s.",
        paste(dv_lhs[i], ":~", partable$rhs[dv_idx[i]]),
        lav_msg_view(bad, "none")))
    }
    available <- c(available, dv_lhs[i])
  }

  # the host rows: all model rows whose label equals a ":~" lhs
  for (i in seq_along(dv_idx)) {
    host_idx <- model_row_idx[partable$label[model_row_idx] == dv_lhs[i]]
    if (length(host_idx) == 0L) {
      lav_msg_stop(gettextf(
        "the left-hand side of the ':~' definition `%1$s' is not a label
        of an existing model parameter.",
        paste(dv_lhs[i], ":~", partable$rhs[dv_idx[i]])))
    }
    if (any(partable$free[host_idx] == 0L)) {
      lav_msg_stop(gettextf(
        "the ':~' operator cannot be used to define the FIXED
        parameter `%s'.", dv_lhs[i]))
    }
    # flag + fix the host rows
    if (is.null(partable$dv)) {
      partable$dv <- character(length(partable$lhs))
    }
    partable$dv[host_idx] <- dv_lhs[i]
    partable$free[host_idx] <- 0L
  }

  # equality constraints between host rows (eg auto-generated plabel
  # constraints when the same label is used in several groups): the
  # constraint is vacuous when both sides are hosts of the SAME
  # definition (each cell holds the same per-observation value), so we
  # drop it; an equality involving a host on only one side is an error
  if (!is.null(partable$plabel)) {
    dv_of_plabel <- partable$dv
    names(dv_of_plabel) <- partable$plabel
    eq_idx <- which(partable$op == "==")
    drop_idx <- integer(0L)
    for (i in eq_idx) {
      lhs_dv <- dv_of_plabel[partable$lhs[i]]
      rhs_dv <- dv_of_plabel[partable$rhs[i]]
      lhs_dv <- if (is.na(lhs_dv)) "" else lhs_dv
      rhs_dv <- if (is.na(rhs_dv)) "" else rhs_dv
      if (nchar(lhs_dv) > 0L && lhs_dv == rhs_dv) {
        drop_idx <- c(drop_idx, i)
      } else if (nchar(lhs_dv) > 0L || nchar(rhs_dv) > 0L) {
        lav_msg_stop(gettextf(
          "equality constraint `%1$s == %2$s' involves the data-defined
          parameter `%3$s'; a moderated parameter can not be constrained
          to be equal to an ordinary parameter.",
          partable$lhs[i], partable$rhs[i],
          if (nchar(lhs_dv) > 0L) lhs_dv else rhs_dv))
      }
    }
    if (length(drop_idx) > 0L) {
      partable <- lapply(partable, function(x) x[-drop_idx])
      # keep the id/position invariant intact
      partable$id <- seq_along(partable$id)
      dv_idx <- which(partable$op == ":~")
    }
  }

  # references from other constraints (==, <, >, :=) to a moderated
  # parameter are not allowed: no single scalar value exists
  other_con_idx <- which(partable$op %in% c("==", "<", ">", ":="))
  for (i in other_con_idx) {
    txt <- c(if (partable$op[i] == ":=") character(0L) else partable$lhs[i],
             partable$rhs[i])
    txt <- txt[nchar(txt) > 0L]
    con_vars <- unique(unlist(lapply(txt, function(s) {
      all.vars(lav_pt_con_parse(s, gettext("constraint(s)")))
    })))
    bad <- con_vars[con_vars %in% dv_lhs]
    if (length(bad) > 0L) {
      lav_msg_stop(gettextf(
        "the constraint/definition `%1$s %2$s %3$s' refers to data-defined
        parameter(s): %4$s; this is not supported.",
        partable$lhs[i], partable$op[i], partable$rhs[i],
        lav_msg_view(bad, "none")))
    }
  }

  # append the component parameters as free scalar rows (op = "dp")
  n_dp <- length(component_names)
  if (n_dp > 0L) {
    dp_idx <- length(partable$id) + seq_len(n_dp)
    partable <- lapply(partable, function(x) {
      if (is.character(x)) {
        c(x, rep("", n_dp))
      } else {
        c(x, rep(NA, n_dp))
      }
    })
    partable$id[dp_idx] <- dp_idx
    partable$lhs[dp_idx] <- component_names
    partable$op[dp_idx] <- rep("dp", n_dp)
    partable$rhs[dp_idx] <- rep("", n_dp)
    partable$user[dp_idx] <- rep(1L, n_dp)
    partable$label[dp_idx] <- component_names
    partable$free[dp_idx] <- rep(1L, n_dp) # renumbered later
    partable$exo[dp_idx] <- rep(0L, n_dp)
    partable$block[dp_idx] <- rep(0L, n_dp)
    for (bl in c("group", "level", "class")) {
      if (!is.null(partable[[bl]])) {
        if (is.character(partable[[bl]])) {
          partable[[bl]][dp_idx] <- rep("", n_dp)
        } else {
          partable[[bl]][dp_idx] <- rep(0L, n_dp)
        }
      }
    }
  }

  partable
}

# starting values for data-defined (":~") parameters:
#  - component ("dp") parameters: ustart if given, 0 otherwise
#  - host rows: a REPRESENTATIVE value = the definition evaluated at the
#    component starting values and the sample means of the definition
#    variables; this keeps the model matrices numerically complete, and
#    provides external backends with a reasonable warm start
# (called after lav_start; needs lavdata for the definition variable means)
lav_start_dv <- function(lavpartable, lavdata = NULL) {
  dv_idx <- which(lavpartable$op == ":~")
  if (length(dv_idx) == 0L) {
    return(lavpartable)
  }
  dv_lhs <- lavpartable$lhs[dv_idx]

  # 1. component parameters
  dp_idx <- which(lavpartable$op == "dp")
  if (length(dp_idx) > 0L) {
    user_idx <- !is.na(lavpartable$ustart[dp_idx])
    lavpartable$start[dp_idx[!user_idx]] <- 0
    lavpartable$start[dp_idx[user_idx]] <-
      lavpartable$ustart[dp_idx[user_idx]]
  }

  # 2. sample means of the definition variables
  data_vars <- unique(unlist(lapply(lavpartable$rhs[dv_idx], function(s) {
    lav_dv_expr_datavars(
      lav_pt_con_parse(s, gettext("parameter definition(s)"))[[1L]])
  })))
  data_means <- vector("list", length(data_vars))
  names(data_means) <- data_vars
  for (v in data_vars) {
    values <- NULL
    if (!is.null(lavdata) && lavdata@data.type == "full") {
      values <- unlist(lapply(seq_len(lavdata@ngroups), function(g) {
        aux_col <- match(v, lavdata@ov.names.aux[[g]])
        if (!is.na(aux_col)) {
          return(lavdata@aux[[g]][, aux_col])
        }
        ov_col <- match(v, lavdata@ov.names[[g]])
        if (!is.na(ov_col)) {
          return(lavdata@X[[g]][, ov_col])
        }
        NULL
      }))
    }
    if (is.null(values)) {
      # no data available (eg standalone call): evaluate at 0
      data_means[[v]] <- 0
    } else {
      data_means[[v]] <- mean(values, na.rm = TRUE)
    }
  }

  # 3. evaluate the definitions (in dependency order) in an environment
  #    holding the starting values of all referenced (labeled) parameters;
  #    for each definition, we first try to improve the (0) start of ONE
  #    'location' component by numerically solving
  #        definition(component, means) == lav_start value of the host
  #    (eg L2_0 becomes the unmoderated loading value, v0 becomes
  #     log(variance start), ...); if this fails, the components stay
  #    at their defaults and the representative value is used as-is
  env <- new.env(parent = parent.frame())
  free_idx <- which(lavpartable$free > 0L & nchar(lavpartable$label) > 0L)
  for (i in free_idx) {
    assign(lavpartable$label[i], lavpartable$start[i], envir = env)
  }
  dp_labels <- lavpartable$lhs[dp_idx]
  # count in how many definitions each component appears (only components
  # unique to a definition are candidates for the location-solve)
  dp_use_count <- integer(length(dp_labels))
  all_refs <- vector("list", length(dv_idx))
  for (i in seq_along(dv_idx)) {
    expr <- lav_pt_con_parse(lavpartable$rhs[dv_idx[i]],
                             gettext("parameter definition(s)"))[[1L]]
    all_refs[[i]] <- setdiff(all.vars(lav_dv_expr_strip_data(expr)),
                             ".dv.data.")
    dp_use_count <- dp_use_count + (dp_labels %in% all_refs[[i]])
  }
  # topological order: definitions referencing other definitions come later
  todo <- seq_along(dv_idx)
  while (length(todo) > 0L) {
    progress <- FALSE
    for (i in todo) {
      refs <- all_refs[[i]][all_refs[[i]] %in% dv_lhs]
      if (any(!refs %in% ls(env))) {
        next # references a definition that has not been evaluated yet
      }
      expr <- lav_pt_con_parse(lavpartable$rhs[dv_idx[i]],
                               gettext("parameter definition(s)"))[[1L]]
      expr <- lav_dv_expr_sub_data(expr, data_means)
      host_idx <- which(!lavpartable$op %in%
                          c("==", "<", ">", ":=", ":~", "dp") &
                        lavpartable$label == dv_lhs[i])
      # try to seed a location component to match the lav_start value
      # of the host (candidate: first component unique to this definition)
      target <- lavpartable$start[host_idx[1L]]
      cand <- all_refs[[i]][all_refs[[i]] %in%
                              dp_labels[dp_use_count == 1L]]
      if (length(cand) > 0L && is.finite(target)) {
        c1 <- cand[1L]
        c1_old <- get(c1, envir = env)
        gfun <- function(t) {
          assign(c1, t, envir = env)
          eval(expr, envir = env) - target
        }
        root <- tryCatch(
          stats::uniroot(gfun, interval = c(-30, 30))$root,
          error = function(e) NULL)
        if (!is.null(root) && is.finite(root)) {
          assign(c1, root, envir = env)
          lavpartable$start[dp_idx[match(c1, dp_labels)]] <- root
        } else {
          assign(c1, c1_old, envir = env)
        }
      }
      value <- tryCatch(eval(expr, envir = env), error = function(e) NA_real_)
      if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
        value <- 0
      }
      assign(dv_lhs[i], value, envir = env)
      # fill in the host rows (and the ":~" row itself, for display)
      lavpartable$start[host_idx] <- value
      lavpartable$start[dv_idx[i]] <- value
      todo <- todo[todo != i]
      progress <- TRUE
    }
    if (!progress) {
      break # cannot happen if lav_pt_dv() validated the system
    }
  }

  lavpartable
}

# build the dv.function from the parameter table (cfr. lav_pt_con_def):
#
#     dv.function(.x., .DATA.) -> n x K matrix
#
# where .x. is the free-parameter vector, .DATA. a data.frame (or named
# list) holding (at least) the definition variable columns, K the number
# of ":~" definitions (columns in parameter-table order, named by the
# defined parameter), and n the number of rows in .DATA.; row i holds the
# per-observation values of the data-defined parameters for observation i
lav_pt_con_dv <- function(partable) {
  # empty function
  dv_function <- function() NULL

  dv_idx <- which(partable$op == ":~")
  if (length(dv_idx) == 0L) {
    return(dv_function)
  }
  dv_lhs <- partable$lhs[dv_idx]

  # parse all expressions; collect (per definition) the data variables,
  # the referenced labels, and the cross-references
  exprs <- vector("list", length(dv_idx))
  data_vars <- vector("list", length(dv_idx))
  ref_vars <- vector("list", length(dv_idx))
  for (i in seq_along(dv_idx)) {
    exprs[[i]] <- lav_pt_con_parse(partable$rhs[dv_idx[i]],
                    gettext("parameter definition(s)"))[[1L]]
    data_vars[[i]] <- lav_dv_expr_datavars(exprs[[i]])
    ref_vars[[i]] <- setdiff(all.vars(lav_dv_expr_strip_data(exprs[[i]])),
                             ".dv.data.")
  }

  # labels (components + plain references, NOT other ":~" lhs) -> .x. idx
  def_labels <- setdiff(unique(unlist(ref_vars)), dv_lhs)
  def_x_idx <- partable$free[match(def_labels, partable$label)]
  if (any(is.na(def_x_idx))) {
    lav_msg_stop(gettext(
      "unknown label(s) in ':~' definition(s):"),
      lav_msg_view(def_labels[which(is.na(def_x_idx))], "none"))
  }
  if (any(def_x_idx == 0)) {
    lav_msg_stop(gettext(
      "non-free parameter(s) in ':~' definition(s):"),
      lav_msg_view(def_labels[which(def_x_idx == 0)], "none"))
  }

  # topological order of the definitions
  adj_mat <- matrix(0L, nrow = length(dv_idx), ncol = length(dv_idx))
  for (i in seq_along(dv_idx)) {
    adj_mat[match(intersect(ref_vars[[i]], dv_lhs), dv_lhs), i] <- 1L
  }
  def_order <- lav_graph_order_adj_mat(adj_mat, warn = FALSE)

  # create function
  formals(dv_function) <- alist(.x. = , .DATA. = , ... = )
  body_txt <- paste("{\n# data-defined (:~) parameters\n\n")

  # 1. bind the referenced labels
  if (length(def_labels) > 0L) {
    body_txt <- paste(body_txt, "# parameter labels\n",
      paste(def_labels, " <- .x.[", def_x_idx, "]", collapse = "\n",
            sep = ""),
      "\n", sep = "")
  }

  # 2. bind the definition variables (mangled names, to avoid clashes)
  all_data_vars <- unique(unlist(data_vars))
  body_txt <- paste(body_txt, "\n# definition variables\n",
    paste(".data.", all_data_vars, ". <- .DATA.[[\"", all_data_vars,
          "\"]]", collapse = "\n", sep = ""),
    "\n", sep = "")

  # 3. write the definitions (in dependency order)
  body_txt <- paste(body_txt, "\n# definitions\n", sep = "")
  for (i in def_order) {
    def_txt <- paste(deparse(lav_dv_expr_data2sym(exprs[[i]]),
                             width.cutoff = 500L), collapse = " ")
    body_txt <- paste(body_txt, dv_lhs[i], " <- ", def_txt, "\n", sep = "")
  }

  # 4. return an n x K matrix (columns in parameter-table order)
  body_txt <- paste(body_txt, "\nout <- cbind(",
    paste(dv_lhs, collapse = ", "), ")\n", sep = "")
  body_txt <- paste(body_txt, "colnames(out) <- c(\"",
    paste(dv_lhs, collapse = "\", \""), "\")\n", sep = "")
  body_txt <- paste(body_txt, "return(out)\n}\n", sep = "")

  body(dv_function) <- parse(file = "", text = body_txt)

  if (lav_debug()) {
    cat("dv.function = \n")
    print(dv_function)
    cat("\n")
  }

  dv_function
}
