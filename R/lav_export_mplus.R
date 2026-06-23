# export a fitted lavaan object to Mplus input (+ a data file Mplus can read)
#
# The generated input, when run in Mplus, reproduces the lavaan model:
#  - a single data file is written, with extra columns for the grouping
#    variable (GROUPING) and/or the cluster variable (CLUSTER) when needed
#  - equality constraints are imposed through shared (clean) parameter labels;
#    other constraints / defined parameters go to a MODEL CONSTRAINT block
#  - multigroup models use group-specific MODEL blocks; multilevel models use
#    %WITHIN% / %BETWEEN% sections (nested within group blocks if both apply)

# main entry point: returns the complete Mplus input file as a string
# (data_file gives the name(s) of the data file(s) referenced in the DATA block)
lav_export_mplus <- function(object, data_file = "data.raw") {
  if (inherits(object, "lavaan")) {
    pt      <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    ngroups <- object@Data@ngroups
    nlevels <- object@Data@nlevels
    opt     <- object@Options
  } else {
    pt      <- as.data.frame(lav_export_check(object), stringsAsFactors = FALSE)
    ngroups <- if (is.null(pt$group)) 1L else lav_pt_ngroups(pt)
    nlevels <- if (is.null(pt$level)) 1L else max(pt$level)
    opt     <- list(estimator = "ML", se = "standard", test = "standard",
                    information = "expected", meanstructure = FALSE)
  }
  nn <- length(pt$lhs)
  if (is.null(pt$free))   pt$free   <- rep(0L, nn)
  if (is.null(pt$ustart)) pt$ustart <- if (!is.null(pt$est)) pt$est else rep(NA, nn)
  if (is.null(pt$label))  pt$label  <- rep("", nn)
  if (is.null(pt$plabel)) pt$plabel <- rep("", nn)
  if (is.null(pt$exo))    pt$exo    <- rep(0L, nn)
  if (is.null(pt$block))  pt$block  <- rep(1L, nn)
  if (is.null(pt$group))  pt$group  <- ifelse(pt$block == 0L, 0L, 1L)
  if (is.null(pt$level))  pt$level  <- ifelse(pt$block == 0L, 0L, 1L)
  pt$label[is.na(pt$label)] <- ""

  # NOTE: variable names are sanitized (dots -> underscores) only when emitting
  # model statements (in lav_export_mplus_block), *after* label/constraint
  # resolution, so that plabel references keep matching.

  is.con <- pt$op %in% c("==", "<", ">", ":=")
  display.label <- lav_export_lavaan_labels(pt, is.con)
  mlabel <- vapply(display.label, lav_export_mplus_label, character(1L),
                   USE.NAMES = FALSE)

  ov.ord <- lav_export_mplus_var(lav_pt_vnames(pt, type = "ov.ord"))

  # ---- MODEL block(s) ----
  body <- character(0L)
  for (g in seq_len(ngroups)) {
    header <- if (g == 1L) "MODEL:" else paste0("MODEL g", g, ":")
    body <- c(body, header)
    if (nlevels > 1L) {
      for (l in seq_len(nlevels)) {
        sec <- if (l == 1L) "%WITHIN%" else "%BETWEEN%"
        body <- c(body, paste0("  ", sec))
        rows <- which(!is.con & pt$group == g & pt$level == l)
        # in Mplus, means/intercepts exist only on the between level
        body <- c(body, lav_export_mplus_block(pt, rows, mlabel, ov.ord,
                                               indent = "    ",
                                               drop.means = (l == 1L)))
      }
    } else {
      rows <- which(!is.con & pt$group == g)
      body <- c(body, lav_export_mplus_block(pt, rows, mlabel, ov.ord,
                                            indent = "  "))
    }
  }

  # ---- MODEL CONSTRAINT block ----
  con.df <- lav_export_constraints_df(pt, is.con, display.label)
  con.block <- ""
  if (nrow(con.df)) {
    lines <- character(0L)
    new.par <- con.df$lhs[con.df$op == ":="]
    if (length(new.par)) {
      lines <- c(lines, paste0("  NEW(", paste(new.par, collapse = " "), ");"))
    }
    for (i in seq_len(nrow(con.df))) {
      op <- con.df$op[i]
      mop <- switch(op, "==" = "=", ":=" = "=", op)
      lhs <- lav_export_mplus_expr(con.df$lhs[i])
      rhs <- lav_export_mplus_expr(con.df$rhs[i])
      lines <- c(lines, paste0("  ", lhs, " ", mop, " ", rhs, ";"))
    }
    con.block <- paste0("MODEL CONSTRAINT:\n", paste(lines, collapse = "\n"),
                        "\n")
  }

  # ---- header (TITLE / DATA / VARIABLE / ANALYSIS) ----
  head <- lav_export_mplus_header(object, pt, data_file = data_file,
                                  ngroups = ngroups, nlevels = nlevels,
                                  ov.ord = ov.ord, opt = opt)

  out <- paste0(head,
                paste(body, collapse = "\n"), "\n",
                con.block,
                "OUTPUT:\n  sampstat standardized tech1;\n")
  class(out) <- c("lavaan.character", "character")
  out
}

# emit Mplus MODEL statements for a set of parameter rows
lav_export_mplus_block <- function(pt, rows, mlabel, ov.ord, indent = "  ",
                                   drop.means = FALSE) {
  if (length(rows) == 0L) {
    return(character(0L))
  }
  out <- character(0L)
  for (i in rows) {
    op <- pt$op[i]
    lhs <- lav_export_mplus_var(pt$lhs[i])
    rhs <- lav_export_mplus_var(pt$rhs[i])

    # skip parameters that Mplus handles implicitly
    if (op %in% c("~~", "~1") && pt$exo[i] == 1L) next            # fixed.x exo
    if (op == "~~" && lhs == rhs && lhs %in% ov.ord) next         # ordinal var
    if (op == "~1" && lhs %in% ov.ord) next                       # ordinal int
    if (op == "~1" && drop.means) next             # within level: no means

    lab <- if (nzchar(mlabel[i])) paste0(" (", mlabel[i], ")") else ""
    mod <- lav_export_mplus_mod(pt$free[i], pt$ustart[i])

    stmt <- switch(op,
      "=~" = paste0(lhs, " BY ", rhs, mod, lab, ";"),
      "<~" = paste0(lhs, " BY ", rhs, mod, lab, ";"),
      "~"  = paste0(lhs, " ON ", rhs, mod, lab, ";"),
      "~~" = if (lhs == rhs) {
               paste0(lhs, mod, lab, ";")
             } else {
               paste0(lhs, " WITH ", rhs, mod, lab, ";")
             },
      "~1" = paste0("[", lhs, mod, "]", lab, ";"),
      "|"  = paste0("[", lhs, "$", sub("t", "", rhs), mod, "]", lab, ";"),
      "~*~" = paste0("{", lhs, mod, lab, "};"),
      NULL
    )
    if (!is.null(stmt)) out <- c(out, paste0(indent, stmt))
  }
  out
}

# Mplus value modifier: @value for fixed, * for free
lav_export_mplus_mod <- function(free, ustart) {
  if (free == 0L) {
    if (is.na(ustart)) {
      return("")
    }
    return(paste0("@", lav_export_num(ustart)))
  }
  "*"
}

# sanitize a variable name for Mplus (no dots)
lav_export_mplus_var <- function(x) {
  gsub("\\.", "_", x)
}

# sanitize a parameter label for Mplus (alnum/underscore, must start w/ letter)
lav_export_mplus_label <- function(x) {
  if (!nzchar(x)) {
    return("")
  }
  x <- gsub("[^A-Za-z0-9_]", "", x)
  if (!grepl("^[A-Za-z]", x)) x <- paste0("p", x)
  x
}

# rewrite a constraint expression: sanitize labels embedded in it and convert
# the few operators that differ between lavaan and Mplus (eg ^ -> **)
lav_export_mplus_expr <- function(expr) {
  if (is.na(expr) || !nzchar(expr)) {
    return(expr)
  }
  expr <- gsub("\\^", "**", expr)
  # drop any leftover dots inside plabel-style tokens
  expr <- gsub("\\.", "", expr)
  expr
}

# build the TITLE / DATA / VARIABLE / ANALYSIS header
lav_export_mplus_header <- function(object, pt, data_file, ngroups, nlevels,
                                    ov.ord, opt) {
  estimator <- if (inherits(object, "lavaan")) {
    lav_mplus_estimator(object)
  } else {
    "ML"
  }
  meanstructure <- isTRUE(opt$meanstructure)
  information <- opt$information
  if (is.null(information)) information <- "expected"

  ov.names <- lav_export_mplus_var(lav_pt_vnames(pt, type = "ov"))
  weight.name <- if (inherits(object, "lavaan")) {
    object@Data@sampling.weights
  } else {
    character(0L)
  }
  weight.name <- lav_export_mplus_var(weight.name)

  # all columns present in the data file (must match the written data)
  data.cols <- ov.names
  if (length(weight.name)) data.cols <- c(data.cols, weight.name)
  if (ngroups > 1L) data.cols <- c(data.cols, "GRP")
  if (nlevels > 1L) data.cols <- c(data.cols, "CLUS")

  data.type <- if (inherits(object, "lavaan")) {
    object@Data@data.type
  } else {
    "full"
  }
  listwise <- isTRUE(opt$missing == "listwise")

  # TITLE
  c_title <- "TITLE:\n  this Mplus input file is autogenerated by lavExport()\n"

  # DATA
  c_data <- "DATA:\n"
  if (ngroups > 1L && length(data_file) > 1L) {
    for (g in seq_len(ngroups)) {
      c_data <- paste0(c_data, "  file (g", g, ") is ", data_file[g], ";\n")
    }
  } else {
    c_data <- paste0(c_data, "  file is ", data_file[1L], ";\n")
  }
  if (data.type == "full") {
    c_data <- paste0(c_data, "  type is individual;\n")
    if (listwise) c_data <- paste0(c_data, "  listwise = on;\n")
  } else if (data.type == "moment") {
    c_data <- paste0(c_data, "  type is fullcov;\n")
    c_data <- paste0(c_data, "  nobservations are ",
                     object@Data@nobs[[1L]], ";\n")
  }

  # VARIABLE
  c_var <- paste0("VARIABLE:\n  names are\n",
                  lav_export_mplus_wrap(data.cols), ";\n")
  c_var <- paste0(c_var, "  usevariables are\n",
                  lav_export_mplus_wrap(ov.names), ";\n")
  if (data.type == "full") {
    c_var <- paste0(c_var, "  missing are all (-999999);\n")
  }
  if (length(ov.ord)) {
    c_var <- paste0(c_var, "  categorical are\n",
                    lav_export_mplus_wrap(ov.ord), ";\n")
  }
  if (ngroups > 1L && length(data_file) <= 1L) {
    glab <- paste0(seq_len(ngroups), "=g", seq_len(ngroups))
    c_var <- paste0(c_var, "  grouping is GRP (", paste(glab, collapse = " "),
                    ");\n")
  }
  if (nlevels > 1L) {
    c_var <- paste0(c_var, "  cluster is CLUS;\n")
    # declare variables that are modeled only at the within or between level
    if (inherits(object, "lavaan")) {
      ovl <- object@Data@ov.names.l[[1L]]
      if (length(ovl) >= 2L) {
        within.only  <- lav_export_mplus_var(setdiff(ovl[[1L]], ovl[[2L]]))
        between.only <- lav_export_mplus_var(setdiff(ovl[[2L]], ovl[[1L]]))
        if (length(within.only)) {
          c_var <- paste0(c_var, "  within = ",
                          paste(within.only, collapse = " "), ";\n")
        }
        if (length(between.only)) {
          c_var <- paste0(c_var, "  between = ",
                          paste(between.only, collapse = " "), ";\n")
        }
      }
    }
  }
  if (length(weight.name)) {
    c_var <- paste0(c_var, "  weight is ", weight.name, ";\n")
  }

  # ANALYSIS
  atype <- if (nlevels > 1L) "twolevel" else "general"
  c_ana <- paste0("ANALYSIS:\n  type = ", atype, ";\n")
  c_ana <- paste0(c_ana, "  estimator = ", toupper(estimator), ";\n")
  if (toupper(estimator) %in% c("ML", "MLR", "MLF") && nlevels == 1L) {
    c_ana <- paste0(c_ana, "  information = ", information[1L], ";\n")
  }
  if (!meanstructure && nlevels == 1L) {
    c_ana <- paste0(c_ana, "  model = nomeanstructure;\n")
  }

  paste0(c_title, c_data, c_var, c_ana)
}

# wrap a vector of names over multiple indented lines (6 per line)
lav_export_mplus_wrap <- function(x, per = 6L, indent = "    ") {
  if (length(x) == 0L) {
    return("")
  }
  grp <- (seq_along(x) - 1L) %/% per
  lines <- tapply(x, grp, function(z) paste0(indent, paste(z, collapse = " ")))
  paste(lines, collapse = "\n")
}


# helper functions ----------------------------------------------------------

lav_mplus_estimator <- function(object) {
  estimator <- object@Options$estimator
  if (estimator == "DWLS") {
    estimator <- "WLS"
  }

  # only 1 argument for 'test' is allowed: drop the (default) "standard" test
  # and keep the scaled/robust one if present
  test <- object@Options$test
  test <- test[test != "standard"]
  if (length(test) == 0L) {
    test <- "standard"
  } else if (length(test) > 1L) {
    test <- test[1L]
  }
  object@Options$test <- test

  if (estimator == "ML") {
    if (object@Options$test %in% c("yuan.bentler", "yuan.bentler.mplus")) {
      estimator <- "MLR"
    } else if (object@Options$test == "satorra.bentler") {
      estimator <- "MLM"
    } else if (object@Options$test == "scaled.shifted") {
      estimator <- "MLMV"
    } else if (object@Options$se == "first.order") {
      estimator <- "MLF"
    }
  } else if (estimator %in% c("ULS", "WLS")) {
    if (object@Options$test == "satorra.bentler") {
      estimator <- paste(estimator, "M", sep = "")
    } else if (object@Options$test == "scaled.shifted") {
      estimator <- paste(estimator, "MV", sep = "")
    }
  } else if (estimator == "MML") {
    estimator <- "ML"
  }

  estimator
}
