# build def function from partable
lav_partable_constraints_def <- function(partable, con = NULL, debug = FALSE, # nolint start
                                         txtOnly = FALSE, warn = TRUE) {      # nolint end
  if (!missing(debug)) {
    current_debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current_debug), TRUE)
  }
  # empty function
  def_function <- function() NULL

  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get := definitions
  def_idx <- which(partable$op == ":=")

  # catch empty def
  if (length(def_idx) == 0L) {
    if (txtOnly) {
      return(character(0L))
    } else {
      return(def_function)
    }
  }

  # sort order of def.idx by dependencies
  deps <- lapply(partable$rhs[def_idx], FUN = function(x) {
                                            all.vars(parse(text = x)) })
  lab_unsorted <- partable$lhs[def_idx]
  adj_mat <- matrix(0L, nrow = length(def_idx), ncol = length(def_idx))
  for (i in seq_along(lab_unsorted)) {
    adj_mat[lab_unsorted %in% deps[[i]], i] <- 1L
  }
  def_idx <- def_idx[lav_graph_order_adj_mat(adj_mat, warn = warn)]

  # create function
  formals(def_function) <- alist(.x. = , ... = )
  if (txtOnly) {
    body_txt <- ""
  } else {
    body_txt <- paste("{\n# parameter definitions\n\n")
  }

  lhs_names <- partable$lhs[def_idx]
  def_labels <- all.vars(parse(file = "", text = partable$rhs[def_idx]))
  # remove the ones in lhs.names
  idx <- which(def_labels %in% lhs_names)
  if (length(idx) > 0L) def_labels <- def_labels[-idx]

  # get corresponding 'x' indices
  def_x_idx <- partable$free[match(def_labels, partable$label)]
  if (any(is.na(def_x_idx))) {
    lav_msg_stop(gettext(
      "unknown label(s) in variable definition(s):"),
      lav_msg_view(def_labels[which(is.na(def_x_idx))], "none")
    )
  }
  if (any(def_x_idx == 0)) {
    lav_msg_stop(gettext(
      "non-free parameter(s) in variable definition(s):"),
      lav_msg_view(def_labels[which(def_x_idx == 0)], "none")
    )
  }
  def_x_lab <- paste(".x.[", def_x_idx, "]", sep = "")
  # put both the labels the function BODY
  if (length(def_x_idx) > 0L) {
    body_txt <- paste(body_txt, "# parameter labels\n",
      paste(def_labels, " <- ", def_x_lab, collapse = "\n"),
      "\n",
      sep = ""
    )
  }

  # write the definitions literally
  body_txt <- paste(body_txt, "\n# parameter definitions\n", sep = "")
  for (i in seq_along(def_idx)) {
    body_txt <- paste(body_txt,
      lhs_names[i], " <- ", partable$rhs[def_idx[i]], "\n",
      sep = ""
    )
  }

  if (txtOnly) {
    return(body_txt)
  }

  # put the results in 'out'
  body_txt <- paste(body_txt, "\nout <- ",
    paste("c(", paste(lhs_names, collapse = ","), ")\n", sep = ""),
    sep = ""
  )
  # what to do with NA values? -> return +Inf???
  body_txt <- paste(body_txt, "out[is.na(out)] <- Inf\n", sep = "")
  body_txt <- paste(body_txt, "names(out) <- ",
    paste("c(\"", paste(lhs_names, collapse = "\",\""), "\")\n", sep = ""),
    sep = ""
  )
  body_txt <- paste(body_txt, "return(out)\n}\n", sep = "")

  body(def_function) <- parse(file = "", text = body_txt)
  if (lav_debug()) {
    cat("def.function = \n")
    print(def_function)
    cat("\n")
  }

  def_function
}

# build ceq function from partable
#     non-trivial equality constraints (linear or nonlinear)
#     convert to 'ceq(x)' function where 'x' is the (free) parameter vector
#     and ceq(x) returns the evaluated equality constraints
#
#     eg. if b1 + b2 == 2 (and b1 correspond to, say,  x[10] and x[17])
#         ceq <- function(x) {
#             out <- rep(NA, 1)
#             b1 = x[10]; b2 = x[17]
#             out[1] <- b1 + b2 - 2
#         }
lav_partable_constraints_ceq <- function(partable, con = NULL, debug = FALSE, # nolint start
                                         txtOnly = FALSE) {                   # nolint end
  if (!missing(debug)) {
    current_debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current_debug), TRUE)
  }
  # empty function
  ceq_function <- function() NULL

  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get equality constraints
  eq_idx <- which(partable$op == "==")

  # catch empty ceq
  if (length(eq_idx) == 0L) {
    if (txtOnly) {
      return(character(0L))
    } else {
      return(ceq_function)
    }
  }

  # create function
  formals(ceq_function) <- alist(.x. = , ... = )
  if (txtOnly) {
    body_txt <- ""
  } else {
    body_txt <- paste("{\nout <- rep(NA, ", length(eq_idx), ")\n", sep = "")
  }

  # first come the variable definitions
  def_txt <- lav_partable_constraints_def(partable, txtOnly = TRUE,
                                          warn = FALSE)
  def_idx <- which(partable$op == ":=")
  body_txt <- paste(body_txt, def_txt, "\n", sep = "")


  # extract labels
  lhs_labels <- all.vars(parse(file = "", text = partable$lhs[eq_idx]))
  rhs_labels <- all.vars(parse(file = "", text = partable$rhs[eq_idx]))
  eq_labels <- unique(c(lhs_labels, rhs_labels))
  # remove def.names from eq.labels
  if (length(def_idx) > 0L) {
    def_names <- as.character(partable$lhs[def_idx])
    d_idx <- which(eq_labels %in% def_names)
    if (length(d_idx) > 0) eq_labels <- eq_labels[-d_idx]
  }
  eq_x_idx <- rep(as.integer(NA), length(eq_labels))
  # get user-labels ids
  ulab_idx <- which(eq_labels %in% partable$label)
  if (length(ulab_idx) > 0L) {
    eq_x_idx[ulab_idx] <- partable$free[match(
      eq_labels[ulab_idx],
      partable$label
    )]
  }
  # get plabels ids
  plab_idx <- which(eq_labels %in% partable$plabel)
  if (length(plab_idx) > 0L) {
    eq_x_idx[plab_idx] <- partable$free[match(
      eq_labels[plab_idx],
      partable$plabel
    )]
  }

  # check if we have found the label
  if (any(is.na(eq_x_idx))) {
    lav_msg_stop(gettext("unknown label(s) in equality constraint(s):"),
      lav_msg_view(eq_labels[which(is.na(eq_x_idx))], "none")
    )
  }
  # check if they are all 'free'
  if (any(eq_x_idx == 0)) {
    fixed_eq_idx <- which(eq_x_idx == 0)
    # FIXME: what should we do here? we used to stop with an error
    # from 0.5.18, we give a warning, and replace the non-free label
    # with its fixed value in ustart
    # warning(
    #    "lavaan WARNING: non-free parameter(s) in equality constraint(s): ",
    #    paste(eq.labels[fixed.eq.idx], collapse=" "))

    fixed_lab_lhs <- eq_labels[fixed_eq_idx]
    fixed_lab_rhs <- numeric(length(fixed_lab_lhs))

    for (i in seq_along(fixed_lab_lhs)) {
      # first try label
      idx <- match(fixed_lab_lhs[i], partable$label)
      # then try plabel
      if (is.na(idx)) {
        idx <- match(fixed_lab_lhs[i], partable$plabel)
      }
      if (is.na(idx)) {
        # hm, not found? fill in zero, or NA?
      } else {
        fixed_lab_rhs[i] <- partable$ustart[idx]
      }
    }

    body_txt <- paste(body_txt, "# non-free parameter labels\n",
      paste(fixed_lab_lhs, "<-", fixed_lab_rhs, collapse = "\n"),
      "\n",
      sep = ""
    )

    eq_x_idx <- eq_x_idx[-fixed_eq_idx]
    eq_labels <- eq_labels[-fixed_eq_idx]
  }

  # put the labels the function BODY
  eq_x_lab <- paste(".x.[", eq_x_idx, "]", sep = "")
  if (length(eq_x_idx) > 0L) {
    body_txt <- paste(body_txt, "# parameter labels\n",
      paste(eq_labels, "<-", eq_x_lab, collapse = "\n"),
      "\n",
      sep = ""
    )
  }

  # write the equality constraints literally
  body_txt <- paste(body_txt, "\n# equality constraints\n", sep = "")
  for (i in seq_along(eq_idx)) {
    lhs <- partable$lhs[eq_idx[i]]
    rhs <- partable$rhs[eq_idx[i]]
    if (rhs == "0") {
      eq_string <- lhs
    } else {
      eq_string <- paste(lhs, " - (", rhs, ")", sep = "")
    }
    body_txt <- paste(body_txt, "out[", i, "] <- ", eq_string, "\n", sep = "")
  }

  if (txtOnly) {
    return(body_txt)
  }

  # put the results in 'out'
  # BODY.txt <- paste(BODY.txt, "\nout <- ",
  #    paste("c(", paste(lhs.names, collapse=","),")\n", sep=""), sep="")

  # what to do with NA values? -> return +Inf???
  body_txt <- paste(body_txt, "\n", "out[is.na(out)] <- Inf\n", sep = "")
  body_txt <- paste(body_txt, "return(out)\n}\n", sep = "")
  body(ceq_function) <- parse(file = "", text = body_txt)
  if (lav_debug()) {
    cat("ceq.function = \n")
    print(ceq_function)
    cat("\n")
  }

  ceq_function
}


# build ciq function from partable
#     non-trivial inequality constraints (linear or nonlinear)
#     convert to 'cin(x)' function where 'x' is the (free) parameter vector
#     and cin(x) returns the evaluated inequality constraints
#
#     eg. if b1 + b2 > 2 (and b1 correspond to, say,  x[10] and x[17])
#         cin <- function(x) {
#             out <- rep(NA, 1)
#             b1 = x[10]; b2 = x[17]
#             out[1] <- b1 + b2 - 2
#         }
#
#     new in 0.6-19: we also add the lower/upper bounds
#
# NOTE: very similar, but not identitical to ceq, because we need to take
#       care of the difference between '<' and '>'
lav_partable_constraints_ciq <- function(partable, con = NULL, debug = FALSE, # nolint start
                                         txtOnly = FALSE) {                   # nolint end
  if (!missing(debug)) {
    current_debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current_debug), TRUE)
  }
  # empty function
  cin_function <- function() NULL

  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get explicit inequality constraints
  ineq_idx <- which(partable$op == ">" | partable$op == "<")

  # get lower/upper bounds
  upper_idx <- integer(0L)
  lower_idx <- integer(0L)
  if (!is.null(partable$upper)) {
    upper_idx <- which(partable$free > 0L  & is.finite(partable$upper))
  }
  if (!is.null(partable$lower)) {
    lower_idx <- which(partable$free > 0L & is.finite(partable$lower))
  }

  # add them to ineq.idx
  ineq_only_idx <- ineq_idx
  ineq_idx <- c(upper_idx, lower_idx, ineq_idx)

  # catch empty ciq
  if (length(ineq_idx) == 0L) {
    if (txtOnly) {
      return(character(0L))
    } else {
      return(cin_function)
    }
  }

  # create function
  formals(cin_function) <- alist(.x. = , ... = )
  if (txtOnly) {
    body_txt <- ""
  } else {
    body_txt <- paste("{\nout <- rep(NA, ", length(ineq_idx), ")\n", sep = "")
  }

  # first come the variable definitions
  def_txt <- lav_partable_constraints_def(partable, txtOnly = TRUE,
                                          warn = FALSE)
  def_idx <- which(partable$op == ":=")
  body_txt <- paste(body_txt, def_txt, "\n", sep = "")

  # extract labels
  lhs_labels <- all.vars(parse(file = "", text = partable$lhs[ineq_only_idx]))
  rhs_labels <- all.vars(parse(file = "", text = partable$rhs[ineq_only_idx]))
  ineq_labels <- unique(c(lhs_labels, rhs_labels))
  # remove def.names from ineq.labels
  if (length(def_idx) > 0L) {
    def_names <- as.character(partable$lhs[def_idx])
    d_idx <- which(ineq_labels %in% def_names)
    if (length(d_idx) > 0) ineq_labels <- ineq_labels[-d_idx]
  }
  ineq_x_idx <- rep(as.integer(NA), length(ineq_labels))
  # get user-labels ids
  ulab_idx <- which(ineq_labels %in% partable$label)
  if (length(ulab_idx) > 0L) {
    ineq_x_idx[ulab_idx] <- partable$free[match(
      ineq_labels[ulab_idx],
      partable$label
    )]
  }
  # get plabels ids
  plab_idx <- which(ineq_labels %in% partable$plabel)
  if (length(plab_idx) > 0L) {
    ineq_x_idx[plab_idx] <- partable$free[match(
      ineq_labels[plab_idx],
      partable$plabel
    )]
  }

  # check if we have found the label
  if (length(ineq_x_idx) > 0L && any(is.na(ineq_x_idx))) {
    lav_msg_stop(gettext("unknown label(s) in inequality constraint(s):"),
      lav_msg_view(ineq_labels[which(is.na(ineq_x_idx))], "none")
    )
  }

  # check if they are all 'free'
  if (length(ineq_x_idx) > 0L && any(ineq_x_idx == 0)) {
    fixed_ineq_idx <- which(ineq_x_idx == 0)
    # FIXME: what should we do here? we used to stop with an error
    # from 0.5.18, we give a warning, and replace the non-free label
    # with its fixed value in ustart
    lav_msg_warn(gettext("non-free parameter(s) in inequality constraint(s):"),
      lav_msg_view(ineq_labels[fixed_ineq_idx], "none")
    )

    fixed_lab_lhs <- ineq_labels[fixed_ineq_idx]
    fixed_lab_rhs <- partable$ustart[match(fixed_lab_lhs, partable$label)]
    body_txt <- paste(body_txt, "# non-free parameter labels\n",
      paste(fixed_lab_lhs, "<-", fixed_lab_rhs, collapse = "\n"),
      "\n",
      sep = ""
    )

    ineq_x_idx <- ineq_x_idx[-fixed_ineq_idx]
    ineq_labels <- ineq_labels[-fixed_ineq_idx]
  }

  # put the labels the function BODY
  if (length(ineq_x_idx) > 0L) {
    ineq_x_lab <- paste(".x.[", ineq_x_idx, "]", sep = "")
    if (length(ineq_x_idx) > 0L) {
      body_txt <- paste(body_txt, "# parameter labels\n",
        paste(ineq_labels, "<-", ineq_x_lab, collapse = "\n"),
        "\n",
        sep = ""
      )
    }
  }

  # write the constraints literally
  body_txt <- paste(body_txt, "\n# inequality constraints\n", sep = "")
  free <- partable$free
  free[free > 0] <- seq_along(free[free > 0])
  for (i in seq_along(ineq_idx)) {
    lhs <- partable$lhs[ineq_idx[i]]
    op <- partable$op[ineq_idx[i]]
    rhs <- partable$rhs[ineq_idx[i]]

    if (ineq_idx[i] %in% ineq_only_idx) {
      # EXPLICIT inequality constraints
      # note,this is different from ==, because we have < AND >
      if (rhs == "0" && op == ">") {
        ineq_string <- lhs
      } else if (rhs == "0" && op == "<") {
        ineq_string <- paste(rhs, " - (", lhs, ")", sep = "")
      } else if (rhs != "0" && op == ">") {
        ineq_string <- paste(lhs, " - (", rhs, ")", sep = "")
      } else if (rhs != "0" && op == "<") {
        ineq_string <- paste(rhs, " - (", lhs, ")", sep = "")
      }
    } else if (ineq_idx[i] %in% upper_idx) {
      # simple upper bound
      val <- partable$upper[ineq_idx[i]]
      xlab <- paste(".x.[", free[ineq_idx[i]], "]", sep = "")
      ineq_string <- paste(val, " - (", xlab, ")", sep = "")
    } else if (ineq_idx[i] %in% lower_idx) {
      # simple lower bound
      val <- partable$lower[ineq_idx[i]]
      xlab <- paste(".x.[", free[ineq_idx[i]], "]", sep = "")
      ineq_string <- paste(xlab, " - (", val, ")", sep = "")
    }

    body_txt <- paste(body_txt, "out[", i, "] <- ", ineq_string, "\n", sep = "")
  }

  if (txtOnly) {
    return(body_txt)
  }

  # put the results in 'out'
  body_txt <- paste(body_txt, "\n", "out[is.na(out)] <- Inf\n", sep = "")
  if (length(upper_idx) > 0L || length(lower_idx) > 0L) {
    bound_idx <- which(ineq_idx %in% c(upper_idx, lower_idx))
    bound_txt <- paste("c(", paste(bound_idx, collapse = ", "), ")", sep = "")
    body_txt <- paste(body_txt, "attr(out, \"bound.idx\") <- ", bound_txt,
                      "\n", sep = "")
  }
  body_txt <- paste(body_txt, "return(out)\n}\n", sep = "")
  body(cin_function) <- parse(file = "", text = body_txt)
  if (lav_debug()) {
    cat("cin.function = \n")
    print(cin_function)
    cat("\n")
  }

  cin_function
}

# return a named vector of the 'free' indices, for the labels that
# are used in a constrained (or optionally a definition)
# (always 0 for definitions)
lav_partable_constraints_label_id <- function(partable, con = NULL, # nolint
                                              def = TRUE) {
  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get constraints
  if (def) {
    con_idx <- which(partable$op %in% c("==", "<", ">", ":="))
  } else {
    con_idx <- which(partable$op %in% c("==", "<", ">"))
  }

  # catch empty con
  if (length(con_idx) == 0L) {
    return(integer(0L))
  }

  def_idx <- which(partable$op == ":=")

  # extract labels
  lhs_labels <- all.vars(parse(file = "", text = partable$lhs[con_idx]))
  rhs_labels <- all.vars(parse(file = "", text = partable$rhs[con_idx]))
  con_labels <- unique(c(lhs_labels, rhs_labels))

  # remove def.names from con.labels (unless def = TRUE)
  if (!def && length(def_idx) > 0L) {
    def_names <- as.character(partable$lhs[def_idx])
    d_idx <- which(con_labels %in% def_names)
    if (length(d_idx) > 0) {
      con_labels <- con_labels[-d_idx]
    }
  }
  con_x_idx <- rep(as.integer(NA), length(con_labels))

  # get user-labels ids
  ulab_idx <- which(con_labels %in% partable$label)
  if (length(ulab_idx) > 0L) {
    con_x_idx[ulab_idx] <- partable$free[match(
      con_labels[ulab_idx],
      partable$label
    )]
  }
  # get plabels ids
  plab_idx <- which(con_labels %in% partable$plabel)
  if (length(plab_idx) > 0L) {
    con_x_idx[plab_idx] <- partable$free[match(
      con_labels[plab_idx],
      partable$plabel
    )]
  }

  # check if we have found the label
  if (any(is.na(con_x_idx))) {
    lav_msg_warn(gettext("unknown label(s) in equality constraint(s):"),
      lav_msg_view(con_labels[which(is.na(con_x_idx))], "none")
    )
  }

  # return named integer vector
  names(con_x_idx) <- con_labels

  con_x_idx
}
