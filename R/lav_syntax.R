# parse lavaan syntax
# YR 14 Jan 2014: move to lav_syntax.R
# YR 17 Oct 2023: add ldw parser
# YR 23 Oct 2024: switch to "c.r"
# LDW 9 Jan 2026: remove c.r option

lavParseModelString <- function(model.syntax = "", as.data.frame. = FALSE,   # nolint
                                parser = "new", warn = TRUE, debug = FALSE) {
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
  parser <- tolower(parser)
  if (!parser %in% c("old", "new")) {
    lav_msg_stop(gettext("parser= argument should
     be \"old\" or \"new\""))
  }

  if (parser == "old") {
    # original/classic parser
    out <- lav_parse_model_string_orig(
      model_syntax = model.syntax,
      as_data_frame = as.data.frame.
    )
  } else {
    # new parser
    out <- lav_parse_model_string(
      model_syntax = model.syntax,
      as_data_frame = as.data.frame.
      )
  }
  out
}

# the 'original' parser (up to 0.6-17)
lav_parse_model_string_orig <- function(model_syntax = "",
                                        as_data_frame = FALSE) {
  # check for empty syntax
  if (length(model_syntax) == 0) {
    lav_msg_stop(gettextf("lavaan ERROR: empty model syntax"))
  }

  # remove comments prior to split:
  # match from comment character to newline, but don't eliminate newline
  model_syntax <- gsub("[#!].*(?=\n)", "", model_syntax, perl = TRUE)

  # replace semicolons with newlines prior to split
  model_syntax <- gsub(";", "\n", model_syntax, fixed = TRUE)

  # remove all whitespace prior to split
  model_syntax <- gsub("[ \t]+", "", model_syntax, perl = TRUE)
  # remove any occurrence of >= 2 consecutive newlines to eliminate
  # blank statements; this retains a blank newline at the beginning,
  # if such exists, but parser will not choke because of start.idx
  model_syntax <- gsub("\n{2,}", "\n", model_syntax, perl = TRUE)

  # replace 'strange' tildes (in some locales) (new in 0.6-6)
  model_syntax <- gsub(pattern = "\u02dc", replacement = "~", model_syntax)

  # break up in lines
  model <- unlist(strsplit(model_syntax, "\n"))

  # check for multi-line formulas: they contain no operator symbol
  # but before we do that, we remove all strings between double quotes
  # to avoid confusion with for example equal("f1=~x1") statements
  # model.simple <- gsub("\\(.*\\)\\*", "MODIFIER*", model)
  model_simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)

  # start.idx <- grep("[~=<>:|%]", model.simple)
  operators <- c(
    "=~", "<~", "~*~", "~~", "~", "==", "<", ">", ":=",
    ":", "\\|", "%"
  )
  lhs_modifiers <- c("efa")
  operators_extra <- c(operators, lhs_modifiers)
  start_idx <- grep(paste(operators_extra, collapse = "|"), model_simple)

  # check for empty start.idx: no operator found (new in 0.6-1)
  if (length(start_idx) == 0L) {
    lav_msg_stop(gettext("lavaan ERROR: model does not contain lavaan
                          syntax (no operators found)"))
  }

  # check for lonely lhs modifiers (only efa() for now):
  # if found, remove the following start.idx
  efa_idx <- grep("efa\\(", model_simple)
  op_idx <- grep(paste(operators, collapse = "|"), model_simple)
  both_idx <- which(efa_idx %in% op_idx)
  if (length(both_idx) > 0L) {
    efa_idx <- efa_idx[-which(efa_idx %in% op_idx)]
  }
  if (length(efa_idx) > 0L) {
    start_idx <- start_idx[-(match(efa_idx, start_idx) + 1L)]
  }

  # check for non-empty string, without an operator in the first lines
  # (new in 0.6-1)
  if (start_idx[1] > 1L) {
    # two possibilities:
    # - we have an empty line (ok)
    # - the element contains no operator (warn!)
    for (el in 1:(start_idx[1] - 1L)) {
      # not empty?
      if (nchar(model_simple[el]) > 0L) {
        lav_msg_warn(gettextf("lavaan WARNING: no operator found in this
                               syntax line: ", model_simple[el], "\n",
                              "     This syntax line will be ignored!"))
      }
    }
  }

  end_idx <- c(start_idx[-1] - 1, length(model))
  model_orig <- model
  model <- character(length(start_idx))
  for (i in seq_along(start_idx)) {
    model[i] <- paste(model_orig[start_idx[i]:end_idx[i]], collapse = "")
  }

  # ok, in all remaining lines, we should have an operator outside the ""
  model_simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)
  idx_wrong <- which(!grepl(
    paste(operators, collapse = "|"),
    model_simple
  ))
  # idx.wrong <- which(!grepl("[~=<>:|%]", model.simple))
  if (length(idx_wrong) > 0) {
    cat("lavaan: missing operator in formula(s):\n")
    print(model[idx_wrong])
    lav_msg_stop(gettext("lavaan ERROR: syntax error in lavaan model syntax"))
  }

  # but perhaps we have a '+' as the first character?
  idx_wrong <- which(grepl("^\\+", model))
  if (length(idx_wrong) > 0) {
    cat("lavaan: some formula(s) start with a plus (+) sign:\n")
    print(model[idx_wrong])
    lav_msg_stop(gettext("lavaan ERROR: syntax error in lavaan model syntax"))
  }


  # main operation: flatten formulas into single bivariate pieces
  # with a left-hand-side (lhs), an operator (eg "=~"), and a
  # right-hand-side (rhs)
  # both lhs and rhs can have a modifier
  flat_lhs <- character(0)
  flat_op <- character(0)
  flat_rhs <- character(0)
  flat_rhs_mod_idx <- integer(0)
  flat_block <- integer(0) # keep track of groups using ":" operator

  flat_fixed <- character(0) # only for display purposes!
  flat_start <- character(0) # only for display purposes!
  flat_lower <- character(0) # only for display purposes!
  flat_upper <- character(0) # only for display purposes!
  flat_label <- character(0) # only for display purposes!
  flat_prior <- character(0)
  flat_efa <- character(0)
  flat_rv <- character(0)
  flat_idx <- 0L
  mod_idx_1 <- 0L
  con_idx <- 0L
  mod_1 <- vector("list", length = 0L)
  con <- vector("list", length = 0L)
  block <- 1L
  block_op <- FALSE
  for (i in seq_along(model)) {
    x <- model[i]
    if (lav_debug()) {
      cat("formula to parse:\n")
      print(x)
      cat("\n")
    }

    # 1. which operator is used?
    line_simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", x)
    # "=~" operator?
    if (grepl("=~", line_simple, fixed = TRUE)) {
      op <- "=~"
      # "<~" operator?
    } else if (grepl("<~", line_simple, fixed = TRUE)) {
      op <- "<~"
    } else if (grepl("~*~", line_simple, fixed = TRUE)) {
      op <- "~*~"
      # "~~" operator?
    } else if (grepl("~~", line_simple, fixed = TRUE)) {
      op <- "~~"
      # "~" operator?
    } else if (grepl("~", line_simple, fixed = TRUE)) {
      op <- "~"
      # "==" operator?
    } else if (grepl("==", line_simple, fixed = TRUE)) {
      op <- "=="
      # "<" operator?
    } else if (grepl("<", line_simple, fixed = TRUE)) {
      op <- "<"
      # ">" operator?
    } else if (grepl(">", line_simple, fixed = TRUE)) {
      op <- ">"
      # ":=" operator?
    } else if (grepl(":=", line_simple, fixed = TRUE)) {
      op <- ":="
      # ":" operator?
    } else if (grepl(":", line_simple, fixed = TRUE)) {
      op <- ":"
      # "|" operator?
    } else if (grepl("|", line_simple, fixed = TRUE)) {
      op <- "|"
      # "%" operator?
    } else if (grepl("%", line_simple, fixed = TRUE)) {
      op <- "%"
    } else {
      lav_msg_stop(gettext("unknown operator in ", model[i]))
    }

    # 2. split by operator (only the *first* occurence!)
    # check first if equal/label modifier has been used on the LEFT!
    if (substr(x, 1, 6) == "label(") {
      lav_msg_stop(gettext("label modifier can not be used on the",
            "left-hand side of the operator"))
    }
    if (op == "|") {
      op_idx <- regexpr("\\|", x)
    } else if (op == "~*~") {
      op_idx <- regexpr("~\\*~", x)
    } else {
      op_idx <- regexpr(op, x)
    }
    lhs <- substr(x, 1L, op_idx - 1L)

    # right-hand side string
    rhs <- substr(x, op_idx + attr(op_idx, "match.length"), nchar(x))

    # check if first character of rhs is '+'; if so, remove silently
    # (for those who copied multiline R input from a website/pdf)
    if (substr(rhs, 1, 1) == "+") {
      rhs <- substr(rhs, 2, nchar(rhs))
    }

    # 2b. if operator is "==" or "<" or ">" or ":=", put it in CON
    if (op == "==" || op == "<" || op == ">" || op == ":=") {
      # remove quotes, if any
      lhs <- gsub("\\\"", "", lhs)
      rhs <- gsub("\\\"", "", rhs)
      con_idx <- con_idx + 1L
      con[[con_idx]] <- list(op = op, lhs = lhs, rhs = rhs, user = 1L)
      next
    }

    # 2c if operator is ":", put it in BLOCK
    if (op == ":") {
      # check if rhs is empty (new in 0.6-4)
      if (nchar(rhs) == 0L) {
        lav_msg_stop(gettextf(
          "syntax contains block identifier %s with missing number/label. The
          correct syntax is: \"LHS: RHS\", where LHS is a block identifier (eg
          group or level), and RHS is the group/level/block number or label.",
          dQuote(lhs))
        )
      }

      # check lhs (new in 0.6-4) - note: class is for nlsem
      lhs_orig <- lhs
      lhs <- tolower(lhs)
      if (!lhs %in% c("group", "level", "block", "class")) {
        lav_msg_stop(gettextf(
          "unknown block identifier: %s. Block identifier should be
           group, level or block.", dQuote(lhs_orig))
        )
      }

      flat_idx <- flat_idx + 1L
      flat_lhs[flat_idx] <- lhs
      flat_op[flat_idx] <- op
      flat_rhs[flat_idx] <- rhs
      flat_fixed[flat_idx] <- ""
      flat_start[flat_idx] <- ""
      flat_lower[flat_idx] <- ""
      flat_upper[flat_idx] <- ""
      flat_label[flat_idx] <- ""
      flat_prior[flat_idx] <- ""
      flat_efa[flat_idx] <- ""
      flat_rv[flat_idx] <- ""
      flat_rhs_mod_idx[flat_idx] <- 0L
      if (block_op) {
        block <- block + 1L
      }
      flat_block[flat_idx] <- block
      block_op <- TRUE
      next
    }

    # 3. parse left hand

    # new in 0.6-3
    # first check if all lhs names are valid (in R); see ?make.names
    # and ?reserved
    # for example, 'NA' is a reserved keyword, and should not be used
    # this usually only happens for latent variable names
    #
    # check should not come earlier, as we do not need it for :,==,<,>,:=
    lhs_1 <- strsplit(lhs, split = "+", fixed = TRUE)[[1]]
    # remove modifiers
    lhs_1 <- gsub("^\\S*\\*", "", lhs_1)
    if (!all(make.names(lhs_1) == lhs_1)) {
      lav_msg_stop(gettextf(
        "lavaan ERROR: left hand side (lhs) of this formula: %1$s %2$s %3$s
        contains either a reserved word (in R) or an illegal character: %4$s.
        See ?reserved for a list of reserved words in R. Please use a variable
        name that is not a reserved word in R and use only characters, digits,
        or the dot symbol.",
        lhs, op, rhs, dQuote(lhs_1[!make.names(lhs_1) == lhs_1])
      ))
    }

    lhs_formula <- as.formula(paste("~", lhs))
    lhs_out <- lav_syntax_parse_rhs(rhs = lhs_formula[[2L]], op = op)
    lhs_names <- names(lhs_out)


    # new in 0.6-4
    # handle LHS modifiers (if any)
    # if(sum(sapply(lhs.out, length)) > 0L) {
    # warning("lavaan WARNING: left-hand side of formula ",
    #         "below contains modifier:\n", x,"\n")
    # }

    # 4. lav_syntax_parse_rhs (as rhs of a single-sided formula)

    # new 0.5-12: before we do this, replace '0.2?' by 'start(0.2)*'
    # requested by the simsem folks
    rhs <- gsub("\\(?([-]?[0-9]*\\.?[0-9]*)\\)?\\?", "start(\\1)\\*", rhs)




    # new in 0.6-6, check for rhs NAMES that are reserved names
    # like in foo =~ in + out
    rhs_1 <- strsplit(rhs, split = "+", fixed = TRUE)[[1]]
    rhs_names <- gsub("^\\S*\\*", "", rhs_1)
    bad <- c("if", "else", "repeat", "while", "function", "for", "in")
    if (any(rhs_names %in% c(bad, "NA"))) { # "NA" added in 0.6-8
      lav_msg_stop(gettextf(
        "Right hand side (rhs) of this formula:\n    ",
        lhs, " ", op, " ", rhs,
        "\n    contains either a reserved word (in R) ",
        "or an illegal character: ",
        dQuote(rhs_names[which(rhs_names %in% bad)[1]]),
        "\n    See ?reserved for a list of reserved words in R",
        "\n    Please use a variable name that is not a reserved word in R",
        "\n    and use only characters, digits, or the dot symbol."
      ))
    }
    # new in 0.6-6, check for rhs LABELS that are reserved names
    # like in foo =~ in*bar
    rhs_1 <- strsplit(rhs, split = "+", fixed = TRUE)[[1]]
    rhs_labels <- gsub("\\*\\S*$", "", rhs_1)
    if (any(rhs_labels %in% bad)) {
      lav_msg_stop(gettextf(
        "Right hand side (rhs) of this formula:\n    ",
        lhs, " ", op, " ", rhs,
        "\n    contains either a reserved word (in R) ",
        "or an illegal character: ",
        dQuote(rhs_names[which(rhs_labels %in% bad)[1]]),
        "\n    See ?reserved for a list of reserved words in R",
        "\n    Please use a variable name that is not a reserved word in R",
        "\n    and use only characters, digits, or the dot symbol."
      ))
    }
    # new in 0.6-12: check for three-way interaction terms (which we do
    # NOT support)
    if (any(grepl(":", rhs_names))) {
      ncolon <- sapply(gregexpr(":", rhs_names), length)
      if (any(ncolon > 1L)) {
        idx <- which(ncolon > 1L)
        lav_msg_stop(gettext(
        "Three-way or higher-order interaction terms (using
multiple colons) are not supported in the lavaan syntax; please manually
construct the product terms yourself in the data.frame, give them an
appropriate name, and then you can use these interaction variables as any
other (observed) variable in the model syntax. Problematic term is: "),
          rhs_names[idx[1]]
        )
      }
    }


    rhs_formula <- as.formula(paste("~", rhs))
    out <- lav_syntax_parse_rhs(rhs = rhs_formula[[2L]], op = op)

    if (lav_debug()) print(out)

    # for each lhs element
    for (l in seq_along(lhs_names)) {
      # for each rhs element
      for (j in seq_along(out)) {
        # catch intercepts
        if (names(out)[j] == "intercept") {
          if (op == "~") {
            rhs_name <- ""
          } else {
            # either number (1), or reserved name?
            lav_msg_stop(gettextf("Right-hand side of
             formula contains an invalid variable name:\n    %s", x))
          }
        } else if (names(out)[j] == "..zero.." && op == "~") {
          rhs_name <- ""
        } else if (names(out)[j] == "..constant.." && op == "~") {
          rhs_name <- ""
        } else {
          rhs_name <- names(out)[j]
        }

        # move this 'check' to post-parse
        # if(op == "|") {
        #    th.name <- paste("t", j, sep="")
        #    if(names(out)[j] != th.name) {
        #        stop("lavaan ERROR: threshold ", j, " of variable ",
        #             sQuote(lhs.names[1]), " should be named ",
        #             sQuote(th.name), "; found ",
        #             sQuote(names(out)[j]), "\n")
        #    }
        # }

        # catch lhs = rhs and op = "=~"
        if (op == "=~" && lhs_names[l] == names(out)[j]) {
          lav_msg_stop(gettextf(
            "Latent variable `%s` can not be measured by itself",
              lhs_names[l]))
        }

        # check if we not already have this combination (in this group)
        # 1. asymmetric (=~, ~, ~1)
        if (op != "~~") {
          idx <- which(flat_lhs == lhs_names[l] &
            flat_op == op &
            flat_block == block &
            flat_rhs == rhs_name)
          if (length(idx) > 0L) {
            lav_msg_stop(gettextf(
              "Duplicate model element in: %s.", model[i]))
          }
        } else {
          # 2. symmetric (~~)
          idx <- which(flat_lhs == rhs_name &
            flat_op == "~~" &
            flat_block == block &
            flat_rhs == lhs_names[l])
          if (length(idx) > 0L) {
            lav_msg_stop(gettextf(
              "Duplicate model element in: %s.", model[i]))
          }
        }

        # check if we have a self-loop (y ~ y)
        if (op %in% c("~", "<~") && rhs_name == lhs_names[l]) {
          # stop("lavaan ERROR: lhs and rhs are the same in: ",
          #     model[i])
          # this breaks pompom package, example uSEM
          lav_msg_warn(gettextf(
            "Lhs and rhs are the same in: %s.",
            model[i]
          ))
        }


        flat_idx <- flat_idx + 1L
        flat_lhs[flat_idx] <- lhs_names[l]
        flat_op[flat_idx] <- op
        flat_rhs[flat_idx] <- rhs_name
        flat_block[flat_idx] <- block
        flat_fixed[flat_idx] <- ""
        flat_start[flat_idx] <- ""
        flat_label[flat_idx] <- ""
        flat_lower[flat_idx] <- ""
        flat_upper[flat_idx] <- ""
        flat_prior[flat_idx] <- ""
        flat_efa[flat_idx] <- ""
        flat_rv[flat_idx] <- ""

        mod <- list()
        rhs_mod <- 0L
        if (length(lhs_out[[l]]$efa) > 0L) {
          mod$efa <- lhs_out[[l]]$efa
          flat_efa[flat_idx] <- paste(mod$efa, collapse = ";")
          rhs_mod <- 1L # despite being a LHS modifier
        }
        if (length(out[[j]]$fixed) > 0L) {
          mod$fixed <- out[[j]]$fixed
          flat_fixed[flat_idx] <- paste(mod$fixed, collapse = ";")
          rhs_mod <- 1L
        }
        if (length(out[[j]]$start) > 0L) {
          mod$start <- out[[j]]$start
          flat_start[flat_idx] <- paste(mod$start, collapse = ";")
          rhs_mod <- 1L
        }
        if (length(out[[j]]$lower) > 0L) {
          mod$lower <- out[[j]]$lower
          flat_lower[flat_idx] <- paste(mod$lower, collapse = ";")
          rhs_mod <- 1L
        }
        if (length(out[[j]]$upper) > 0L) {
          mod$upper <- out[[j]]$upper
          flat_upper[flat_idx] <- paste(mod$upper, collapse = ";")
          rhs_mod <- 1L
        }
        if (length(out[[j]]$label) > 0L) {
          mod$label <- out[[j]]$label
          flat_label[flat_idx] <- paste(mod$label, collapse = ";")
          rhs_mod <- 1L
        }
        if (length(out[[j]]$rv) > 0L) {
          mod$rv <- out[[j]]$rv
          flat_rv[flat_idx] <- paste(mod$rv, collapse = ";")
          rhs_mod <- 1L
        }
        if (length(out[[j]]$prior) > 0L) {
          mod$prior <- out[[j]]$prior
          flat_prior[flat_idx] <- paste(mod$prior, collapse = ";")
          rhs_mod <- 1L
        }
        # if(op == "~1" && rhs == "0") {
        #    mod$fixed <- 0
        #    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
        #    rhs.mod <- 1L
        # }
        if (op == "=~" && rhs == "0") {
          mod$fixed <- 0
          flat_rhs[flat_idx] <- flat_lhs[flat_idx]
          flat_fixed[flat_idx] <- paste(mod$fixed, collapse = ";")
          rhs_mod <- 1L
        }

        flat_rhs_mod_idx[flat_idx] <- rhs_mod

        if (rhs_mod > 0L) {
          mod_idx_1 <- mod_idx_1 + 1L
          mod_1[[mod_idx_1]] <- mod
        }
      } # rhs elements
    } # lhs elements
  } # model elements

  # enumerate modifier indices
  mod_idx <- which(flat_rhs_mod_idx > 0L)
  flat_rhs_mod_idx[mod_idx] <- seq_along(mod_idx)

  flat <- list(
    lhs = flat_lhs, op = flat_op, rhs = flat_rhs,
    mod.idx = flat_rhs_mod_idx, block = flat_block,
    fixed = flat_fixed, start = flat_start,
    lower = flat_lower, upper = flat_upper,
    label = flat_label, prior = flat_prior,
    efa = flat_efa, rv = flat_rv
  )

  # change op for intercepts (for convenience only)
  int_idx <- which(flat$op == "~" & flat$rhs == "")
  if (length(int_idx) > 0L) {
    flat$op[int_idx] <- "~1"
  }

  # new in 0.6, reorder covariances here!
  flat <- lav_partable_covariance_reorder(flat)

  if (as_data_frame) {
    flat <- as.data.frame(flat, stringsAsFactors = FALSE)
  }

  # new in 0.6-4: check for 'group' within 'level'
  if (any(flat$op == ":")) {
    op_idx <- which(flat$op == ":")
    if (length(op_idx) < 2L) {
      # only 1 block identifier? this is weird -> give warning
      lav_msg_warn(gettextf(
        "Syntax contains only a single block identifier: %s",
         flat$lhs[op_idx]))
    } else {
      first_block <- flat$lhs[op_idx[1L]]
      second_block <- flat$lhs[op_idx[2L]]
      if (first_block == "level" &&
        second_block == "group") {
        lav_msg_stop(gettext(
          "Groups can not be nested within levels."))
      }
    }
  }

  attr(flat, "modifiers") <- mod_1
  attr(flat, "constraints") <- con

  flat
}

lav_syntax_parse_rhs <- function(rhs, op = "") {
  # new version YR 15 dec 2011!
  # - no 'equal' field anymore (only labels!)
  # - every modifier is evaluated
  # - unquoted labels are allowed (eg. x1 + x2 + c(v1,v2,v3)*x3)

  # fill in rhs list
  out <- list()
  repeat {
    if (length(rhs) == 1L) { # last one and only a single element
      out <- c(vector("list", 1L), out)
      name <- all.vars(rhs)
      if (length(name) > 0L) {
        names(out)[1L] <- name
      } else { # intercept or zero?
        if (as.character(rhs) == "1") {
          names(out)[1L] <- "intercept"
        } else if (as.character(rhs) == "0") {
          names(out)[1L] <- "..zero.."
          out[[1L]]$fixed <- 0
        } else {
          names(out)[1L] <- "..constant.."
          out[[1L]]$fixed <- 0
        }
      }
      break
    } else if (rhs[[1L]] == "*") { # last one, but with modifier
      out <- c(vector("list", 1L), out)
      name <- all.vars(rhs[[3L]])

      if (length(name) > 0L) { # not an intercept
        # catch interaction term
        rhs3_names <- all.names(rhs[[3L]])
        if (rhs3_names[1L] == ":") {
          if (length(name) == 1) {
            name <- paste(name[1L], ":", name[1L], sep = "")
          } else {
            name <- paste(name[1L], ":", name[2L], sep = "")
          }
        }
        names(out)[1L] <- name
      } else { # intercept
        names(out)[1L] <- "intercept"
      }
      i_var <- all.vars(rhs[[2L]], unique = FALSE)
      if (length(i_var) > 0L) {
        # modifier are unquoted labels
        out[[1L]]$label <- i_var
      } else {
        # modifer is something else
        out[[1L]] <- lav_syntax_get_modifier(rhs[[2L]])
      }
      break
    } else if (rhs[[1L]] == ":") { # last one, but interaction term
      out <- c(vector("list", 1L), out)
      name <- all.vars(rhs)
      if (length(name) == 1) {
        name <- paste(name[1L], ":", name[1L], sep = "")
      } else {
        name <- paste(name[1L], ":", name[2L], sep = "")
      }
      names(out)[1L] <- name
      break
    } else if (rhs[[1L]] == "+") { # not last one!

      # three possibilities:
      # 1. length(rhs[[3]] == 3), and rhs[[3L]][[1]] == "*" -> modifier
      # 2. length(rhs[[3]] == 3), and rhs[[3L]][[1]] == ":" -> interaction
      # 3. length(rhs[[3]] == 1) -> single element

      out <- c(vector("list", 1L), out)

      # modifier or not?
      if (length(rhs[[3L]]) == 3L && rhs[[3L]][[1]] == "*") {
        # modifier!!
        name <- all.vars(rhs[[3L]][[3]])

        if (length(name) > 0L) { # not an intercept
          # catch interaction term
          rhs3_names <- all.names(rhs[[3L]][[3]])
          if (rhs3_names[1L] == ":") {
            if (length(name) == 1) {
              name <- paste(name[1L], ":", name[1L], sep = "")
            } else {
              name <- paste(name[1L], ":", name[2L], sep = "")
            }
          }
          names(out)[1L] <- name
        } else { # intercept
          names(out)[1L] <- "intercept"
        }
        i_var <- all.vars(rhs[[3]][[2L]], unique = FALSE)
        if (length(i_var) > 0L) {
          # modifier are unquoted labels
          out[[1L]]$label <- i_var
        } else {
          # modifer is something else
          out[[1L]] <- lav_syntax_get_modifier(rhs[[3]][[2L]])
        }

        # interaction term?
      } else if (length(rhs[[3L]]) == 3L && rhs[[3L]][[1]] == ":") {
        # interaction term, without modifier
        name <- all.vars(rhs[[3L]])
        if (length(name) == 1) {
          name <- paste(name[1L], ":", name[1L], sep = "")
        } else {
          name <- paste(name[1L], ":", name[2L], sep = "")
        }
        names(out)[1L] <- name
      } else { # no modifier!!
        name <- all.vars(rhs[[3]])
        if (length(name) > 0L) {
          names(out)[1L] <- name
        } else { # intercept or zero?
          if (as.character(rhs[[3]]) == "1") {
            names(out)[1L] <- "intercept"
          } else if (as.character(rhs[[3]]) == "0") {
            names(out)[1L] <- "..zero.."
            out[[1L]]$fixed <- 0
          } else {
            names(out)[1L] <- "..constant.."
            out[[1L]]$fixed <- 0
          }
        }
      }


      # next element
      rhs <- rhs[[2L]]
    } else {
      lav_msg_stop(gettextf(
        "lavaan ERROR: I'm confused parsing this line: ", rhs, "\n"))
    }
  }

  # if multiple elements, check for duplicated elements and merge if found
  if (length(out) > 1L) {
    rhs_names <- names(out)
    while (!is.na(idx <- which(duplicated(rhs_names))[1L])) {
      dup_name <- rhs_names[idx]
      orig_idx <- match(dup_name, rhs_names)
      merged <- c(out[[orig_idx]], out[[idx]])
      if (!is.null(merged)) { # be careful, NULL will delete element
        out[[orig_idx]] <- merged
      }
      out <- out[-idx]
      rhs_names <- names(out)
    }
  }

  # if thresholds, check order and reorder if necessary
  # if(op == "|") {
  #    t.names <- names(out)
  #    idx <- match(sort(t.names), t.names)
  #    out <- out[idx]
  # }

  out
}


lav_syntax_get_modifier <- function(mod) {
  if (length(mod) == 1L) {
    # three possibilites: 1) numeric, 2) NA, or 3) quoted character
    if (is.numeric(mod)) {
      return(list(fixed = mod))
    }
    if (is.na(mod)) {
      return(list(fixed = as.numeric(NA)))
    }
    if (is.character(mod)) {
      list(label = mod)
    }
  } else if (mod[[1L]] == "start") {
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    list(start = cof)
  } else if (mod[[1L]] == "lower") {
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    list(lower = cof)
  } else if (mod[[1L]] == "upper") {
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    list(upper = cof)
  } else if (mod[[1L]] == "equal") {
    label <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    list(label = label)
  } else if (mod[[1L]] == "label") {
    label <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    label[is.na(label)] <- "" # catch 'NA' elements in a label
    list(label = label)
  } else if (mod[[1L]] == "rv") {
    rv <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    if (anyNA(rv)) {
      lav_msg_stop(gettextf("lavaan ERROR: some rv() labels are NA"))
    }
    list(rv = rv)
  } else if (mod[[1L]] == "prior") {
    prior <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    list(prior = prior)
  } else if (mod[[1L]] == "efa") {
    efa <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    list(efa = efa)
  } else if (mod[[1L]] == "c") {
    # vector: we allow numeric and character only!
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    if (all(is.na(cof))) {
      list(fixed = rep(as.numeric(NA), length(cof)))
    } else if (is.numeric(cof)) {
      list(fixed = cof)
    } else if (is.character(cof)) {
      cof[is.na(cof)] <- "" # catch 'NA' elements in a label
      list(label = cof)
    } else {
      lav_msg_stop(gettextf("lavaan ERROR: can not parse modifier:", mod, "\n"))
    }
  } else {
    # unknown expression
    # as a final attempt, we will evaluate it and coerce it
    # to either a numeric or character (vector)
    cof <- try(eval(mod, envir = NULL, enclos = NULL), silent = TRUE)
    if (inherits(cof, "try-error")) {
      lav_msg_stop(gettextf(
        "lavaan ERROR: evaluating modifier failed: ",
        paste(as.character(mod)[[1]], "()*", sep = ""), "\n"
      ))
    } else if (is.numeric(cof)) {
      list(fixed = cof)
    } else if (is.character(cof)) {
      list(label = cof)
    } else {
      lav_msg_stop(gettextf(
        "lavaan ERROR: can not parse modifier: ",
        paste(as.character(mod)[[1]], "()*", sep = ""), "\n"
      ))
    }
  }
}
