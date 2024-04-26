# build def function from partable
lav_partable_constraints_def <- function(partable, con = NULL, debug = FALSE,
                                         txtOnly = FALSE) {
  # empty function
  def.function <- function() NULL

  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get := definitions
  def.idx <- which(partable$op == ":=")

  # catch empty def
  if (length(def.idx) == 0L) {
    if (txtOnly) {
      return(character(0L))
    } else {
      return(def.function)
    }
  }

  # create function
  formals(def.function) <- alist(.x. = , ... = )
  if (txtOnly) {
    BODY.txt <- ""
  } else {
    BODY.txt <- paste("{\n# parameter definitions\n\n")
  }

  lhs.names <- partable$lhs[def.idx]
  def.labels <- all.vars(parse(file = "", text = partable$rhs[def.idx]))
  # remove the ones in lhs.names
  idx <- which(def.labels %in% lhs.names)
  if (length(idx) > 0L) def.labels <- def.labels[-idx]

  # get corresponding 'x' indices
  def.x.idx <- partable$free[match(def.labels, partable$label)]
  if (any(is.na(def.x.idx))) {
    lav_msg_stop(gettext(
      "unknown label(s) in variable definition(s):"),
      lav_msg_view(def.labels[which(is.na(def.x.idx))], "none")
    )
  }
  if (any(def.x.idx == 0)) {
    lav_msg_stop(gettext(
      "non-free parameter(s) in variable definition(s):"),
      lav_msg_view(def.labels[which(def.x.idx == 0)], "none")
    )
  }
  def.x.lab <- paste(".x.[", def.x.idx, "]", sep = "")
  # put both the labels the function BODY
  if (length(def.x.idx) > 0L) {
    BODY.txt <- paste(BODY.txt, "# parameter labels\n",
      paste(def.labels, " <- ", def.x.lab, collapse = "\n"),
      "\n",
      sep = ""
    )
  }

  # write the definitions literally
  BODY.txt <- paste(BODY.txt, "\n# parameter definitions\n", sep = "")
  for (i in 1:length(def.idx)) {
    BODY.txt <- paste(BODY.txt,
      lhs.names[i], " <- ", partable$rhs[def.idx[i]], "\n",
      sep = ""
    )
  }

  if (txtOnly) {
    return(BODY.txt)
  }

  # put the results in 'out'
  BODY.txt <- paste(BODY.txt, "\nout <- ",
    paste("c(", paste(lhs.names, collapse = ","), ")\n", sep = ""),
    sep = ""
  )
  # what to do with NA values? -> return +Inf???
  BODY.txt <- paste(BODY.txt, "out[is.na(out)] <- Inf\n", sep = "")
  BODY.txt <- paste(BODY.txt, "names(out) <- ",
    paste("c(\"", paste(lhs.names, collapse = "\",\""), "\")\n", sep = ""),
    sep = ""
  )
  BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep = "")

  body(def.function) <- parse(file = "", text = BODY.txt)
  if (debug) {
    cat("def.function = \n")
    print(def.function)
    cat("\n")
  }

  def.function
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
lav_partable_constraints_ceq <- function(partable, con = NULL, debug = FALSE,
                                         txtOnly = FALSE) {
  # empty function
  ceq.function <- function() NULL

  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get equality constraints
  eq.idx <- which(partable$op == "==")

  # catch empty ceq
  if (length(eq.idx) == 0L) {
    if (txtOnly) {
      return(character(0L))
    } else {
      return(ceq.function)
    }
  }

  # create function
  formals(ceq.function) <- alist(.x. = , ... = )
  if (txtOnly) {
    BODY.txt <- ""
  } else {
    BODY.txt <- paste("{\nout <- rep(NA, ", length(eq.idx), ")\n", sep = "")
  }

  # first come the variable definitions
  DEF.txt <- lav_partable_constraints_def(partable, txtOnly = TRUE)
  def.idx <- which(partable$op == ":=")
  BODY.txt <- paste(BODY.txt, DEF.txt, "\n", sep = "")


  # extract labels
  lhs.labels <- all.vars(parse(file = "", text = partable$lhs[eq.idx]))
  rhs.labels <- all.vars(parse(file = "", text = partable$rhs[eq.idx]))
  eq.labels <- unique(c(lhs.labels, rhs.labels))
  # remove def.names from eq.labels
  if (length(def.idx) > 0L) {
    def.names <- as.character(partable$lhs[def.idx])
    d.idx <- which(eq.labels %in% def.names)
    if (length(d.idx) > 0) eq.labels <- eq.labels[-d.idx]
  }
  eq.x.idx <- rep(as.integer(NA), length(eq.labels))
  # get user-labels ids
  ulab.idx <- which(eq.labels %in% partable$label)
  if (length(ulab.idx) > 0L) {
    eq.x.idx[ulab.idx] <- partable$free[match(
      eq.labels[ulab.idx],
      partable$label
    )]
  }
  # get plabels ids
  plab.idx <- which(eq.labels %in% partable$plabel)
  if (length(plab.idx) > 0L) {
    eq.x.idx[plab.idx] <- partable$free[match(
      eq.labels[plab.idx],
      partable$plabel
    )]
  }

  # check if we have found the label
  if (any(is.na(eq.x.idx))) {
    lav_msg_stop(gettext("unknown label(s) in equality constraint(s):"),
      lav_msg_view(eq.labels[which(is.na(eq.x.idx))], "none")
    )
  }
  # check if they are all 'free'
  if (any(eq.x.idx == 0)) {
    fixed.eq.idx <- which(eq.x.idx == 0)
    # FIXME: what should we do here? we used to stop with an error
    # from 0.5.18, we give a warning, and replace the non-free label
    # with its fixed value in ustart
    # warning("lavaan WARNING: non-free parameter(s) in equality constraint(s): ",
    #    paste(eq.labels[fixed.eq.idx], collapse=" "))

    fixed.lab.lhs <- eq.labels[fixed.eq.idx]
    fixed.lab.rhs <- numeric(length(fixed.lab.lhs))

    for (i in 1:length(fixed.lab.lhs)) {
      # first try label
      idx <- match(fixed.lab.lhs[i], partable$label)
      # then try plabel
      if (is.na(idx)) {
        idx <- match(fixed.lab.lhs[i], partable$plabel)
      }
      if (is.na(idx)) {
        # hm, not found? fill in zero, or NA?
      } else {
        fixed.lab.rhs[i] <- partable$ustart[idx]
      }
    }

    BODY.txt <- paste(BODY.txt, "# non-free parameter labels\n",
      paste(fixed.lab.lhs, "<-", fixed.lab.rhs, collapse = "\n"),
      "\n",
      sep = ""
    )

    eq.x.idx <- eq.x.idx[-fixed.eq.idx]
    eq.labels <- eq.labels[-fixed.eq.idx]
  }

  # put the labels the function BODY
  eq.x.lab <- paste(".x.[", eq.x.idx, "]", sep = "")
  if (length(eq.x.idx) > 0L) {
    BODY.txt <- paste(BODY.txt, "# parameter labels\n",
      paste(eq.labels, "<-", eq.x.lab, collapse = "\n"),
      "\n",
      sep = ""
    )
  }

  # write the equality constraints literally
  BODY.txt <- paste(BODY.txt, "\n# equality constraints\n", sep = "")
  for (i in 1:length(eq.idx)) {
    lhs <- partable$lhs[eq.idx[i]]
    rhs <- partable$rhs[eq.idx[i]]
    if (rhs == "0") {
      eq.string <- lhs
    } else {
      eq.string <- paste(lhs, " - (", rhs, ")", sep = "")
    }
    BODY.txt <- paste(BODY.txt, "out[", i, "] <- ", eq.string, "\n", sep = "")
  }

  if (txtOnly) {
    return(BODY.txt)
  }

  # put the results in 'out'
  # BODY.txt <- paste(BODY.txt, "\nout <- ",
  #    paste("c(", paste(lhs.names, collapse=","),")\n", sep=""), sep="")

  # what to do with NA values? -> return +Inf???
  BODY.txt <- paste(BODY.txt, "\n", "out[is.na(out)] <- Inf\n", sep = "")
  BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep = "")
  body(ceq.function) <- parse(file = "", text = BODY.txt)
  if (debug) {
    cat("ceq.function = \n")
    print(ceq.function)
    cat("\n")
  }

  ceq.function
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
# NOTE: very similar, but not identitical to ceq, because we need to take
#       care of the difference between '<' and '>'
lav_partable_constraints_ciq <- function(partable, con = NULL, debug = FALSE,
                                         txtOnly = FALSE) {
  # empty function
  cin.function <- function() NULL

  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get inequality constraints
  ineq.idx <- which(partable$op == ">" | partable$op == "<")

  # catch empty ciq
  if (length(ineq.idx) == 0L) {
    if (txtOnly) {
      return(character(0L))
    } else {
      return(cin.function)
    }
  }

  # create function
  formals(cin.function) <- alist(.x. = , ... = )
  if (txtOnly) {
    BODY.txt <- ""
  } else {
    BODY.txt <- paste("{\nout <- rep(NA, ", length(ineq.idx), ")\n", sep = "")
  }

  # first come the variable definitions
  DEF.txt <- lav_partable_constraints_def(partable, txtOnly = TRUE)
  def.idx <- which(partable$op == ":=")
  BODY.txt <- paste(BODY.txt, DEF.txt, "\n", sep = "")

  # extract labels
  lhs.labels <- all.vars(parse(file = "", text = partable$lhs[ineq.idx]))
  rhs.labels <- all.vars(parse(file = "", text = partable$rhs[ineq.idx]))
  ineq.labels <- unique(c(lhs.labels, rhs.labels))
  # remove def.names from ineq.labels
  if (length(def.idx) > 0L) {
    def.names <- as.character(partable$lhs[def.idx])
    d.idx <- which(ineq.labels %in% def.names)
    if (length(d.idx) > 0) ineq.labels <- ineq.labels[-d.idx]
  }
  ineq.x.idx <- rep(as.integer(NA), length(ineq.labels))
  # get user-labels ids
  ulab.idx <- which(ineq.labels %in% partable$label)
  if (length(ulab.idx) > 0L) {
    ineq.x.idx[ulab.idx] <- partable$free[match(
      ineq.labels[ulab.idx],
      partable$label
    )]
  }
  # get plabels ids
  plab.idx <- which(ineq.labels %in% partable$plabel)
  if (length(plab.idx) > 0L) {
    ineq.x.idx[plab.idx] <- partable$free[match(
      ineq.labels[plab.idx],
      partable$plabel
    )]
  }

  # check if we have found the label
  if (any(is.na(ineq.x.idx))) {
    lav_msg_stop(gettext("unknown label(s) in inequality constraint(s):"),
      lav_msg_view(ineq.labels[which(is.na(ineq.x.idx))], "none")
    )
  }
  # check if they are all 'free'
  if (any(ineq.x.idx == 0)) {
    fixed.ineq.idx <- which(ineq.x.idx == 0)
    # FIXME: what should we do here? we used to stop with an error
    # from 0.5.18, we give a warning, and replace the non-free label
    # with its fixed value in ustart
    lav_msg_warn(gettext("non-free parameter(s) in inequality constraint(s):"),
      lav_msg_view(ineq.labels[fixed.ineq.idx],"none")
    )

    fixed.lab.lhs <- ineq.labels[fixed.ineq.idx]
    fixed.lab.rhs <- partable$ustart[match(fixed.lab.lhs, partable$label)]
    BODY.txt <- paste(BODY.txt, "# non-free parameter labels\n",
      paste(fixed.lab.lhs, "<-", fixed.lab.rhs, collapse = "\n"),
      "\n",
      sep = ""
    )

    ineq.x.idx <- ineq.x.idx[-fixed.ineq.idx]
    ineq.labels <- ineq.labels[-fixed.ineq.idx]
  }

  # put the labels the function BODY
  ineq.x.lab <- paste(".x.[", ineq.x.idx, "]", sep = "")
  if (length(ineq.x.idx) > 0L) {
    BODY.txt <- paste(BODY.txt, "# parameter labels\n",
      paste(ineq.labels, "<-", ineq.x.lab, collapse = "\n"),
      "\n",
      sep = ""
    )
  }

  # write the constraints literally
  BODY.txt <- paste(BODY.txt, "\n# inequality constraints\n", sep = "")
  for (i in 1:length(ineq.idx)) {
    lhs <- partable$lhs[ineq.idx[i]]
    op <- partable$op[ineq.idx[i]]
    rhs <- partable$rhs[ineq.idx[i]]

    # note,this is different from ==, because we have < AND >
    if (rhs == "0" && op == ">") {
      ineq.string <- lhs
    } else if (rhs == "0" && op == "<") {
      ineq.string <- paste(rhs, " - (", lhs, ")", sep = "")
    } else if (rhs != "0" && op == ">") {
      ineq.string <- paste(lhs, " - (", rhs, ")", sep = "")
    } else if (rhs != "0" && op == "<") {
      ineq.string <- paste(rhs, " - (", lhs, ")", sep = "")
    }

    BODY.txt <- paste(BODY.txt, "out[", i, "] <- ", ineq.string, "\n", sep = "")
  }

  if (txtOnly) {
    return(BODY.txt)
  }

  # put the results in 'out'
  # BODY.txt <- paste(BODY.txt, "\nout <- ",
  #    paste("c(", paste(lhs.names, collapse=","),")\n", sep=""), sep="")

  # what to do with NA values? -> return +Inf???
  BODY.txt <- paste(BODY.txt, "\n", "out[is.na(out)] <- Inf\n", sep = "")
  BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep = "")
  body(cin.function) <- parse(file = "", text = BODY.txt)
  if (debug) {
    cat("cin.function = \n")
    print(cin.function)
    cat("\n")
  }

  cin.function
}

# return a named vector of the 'free' indices, for the labels that
# are used in a constrained (or optionally a definition)
# (always 0 for definitions)
lav_partable_constraints_label_id <- function(partable, con = NULL,
                                              def = TRUE,
                                              warn = TRUE) {
  # if 'con', merge partable + con
  if (!is.null(con)) {
    partable$lhs <- c(partable$lhs, con$lhs)
    partable$op <- c(partable$op, con$op)
    partable$rhs <- c(partable$rhs, con$rhs)
  }

  # get constraints
  if (def) {
    con.idx <- which(partable$op %in% c("==", "<", ">", ":="))
  } else {
    con.idx <- which(partable$op %in% c("==", "<", ">"))
  }

  # catch empty con
  if (length(con.idx) == 0L) {
    return(integer(0L))
  }

  def.idx <- which(partable$op == ":=")

  # extract labels
  lhs.labels <- all.vars(parse(file = "", text = partable$lhs[con.idx]))
  rhs.labels <- all.vars(parse(file = "", text = partable$rhs[con.idx]))
  con.labels <- unique(c(lhs.labels, rhs.labels))

  # remove def.names from con.labels (unless def = TRUE)
  if (!def && length(def.idx) > 0L) {
    def.names <- as.character(partable$lhs[def.idx])
    d.idx <- which(con.labels %in% def.names)
    if (length(d.idx) > 0) {
      con.labels <- con.labels[-d.idx]
    }
  }
  con.x.idx <- rep(as.integer(NA), length(con.labels))

  # get user-labels ids
  ulab.idx <- which(con.labels %in% partable$label)
  if (length(ulab.idx) > 0L) {
    con.x.idx[ulab.idx] <- partable$free[match(
      con.labels[ulab.idx],
      partable$label
    )]
  }
  # get plabels ids
  plab.idx <- which(con.labels %in% partable$plabel)
  if (length(plab.idx) > 0L) {
    con.x.idx[plab.idx] <- partable$free[match(
      con.labels[plab.idx],
      partable$plabel
    )]
  }

  # check if we have found the label
  if (any(is.na(con.x.idx)) && warn) {
    lav_msg_warn(gettext("unknown label(s) in equality constraint(s):"),
      lav_msg_view(con.labels[which(is.na(con.x.idx))], "none")
    )
  }

  # return named integer vector
  names(con.x.idx) <- con.labels

  con.x.idx
}
