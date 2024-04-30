# parse lavaan syntax
# YR 14 Jan 2014: move to lav_syntax.R
# YR 17 Oct 2023: add ldw parser

lavParseModelString <- function(model.syntax = "", as.data.frame. = FALSE,
                                parser = "new",
                                warn = TRUE, debug = FALSE) {
  parser <- tolower(parser)
  if (!parser %in% c("old", "new")) {
    lav_msg_stop(gettext("parser= argument should be \"old\" or \"new\""))
  }

  if (parser == "old") {
    # original/classic parser
    out <- lav_parse_model_string_orig(
      model.syntax = model.syntax,
      as.data.frame. = as.data.frame.,
      warn = warn,
      debug = debug
    )
  } else {
    # new parser
    out <- ldw_parse_model_string(
      model.syntax = model.syntax,
      as.data.frame. = as.data.frame.,
      warn = warn,
      debug = debug
    )
  }

  out
}

# the 'original' parser (up to 0.6-17)
lav_parse_model_string_orig <- function(model.syntax = "",
                                        as.data.frame. = FALSE,
                                        warn = TRUE, debug = FALSE) {
  # check for empty syntax
  if (length(model.syntax) == 0) {
    stop("lavaan ERROR: empty model syntax")
  }

  # remove comments prior to split:
  # match from comment character to newline, but don't eliminate newline
  model.syntax <- gsub("[#!].*(?=\n)", "", model.syntax, perl = TRUE)

  # replace semicolons with newlines prior to split
  model.syntax <- gsub(";", "\n", model.syntax, fixed = TRUE)

  # remove all whitespace prior to split
  model.syntax <- gsub("[ \t]+", "", model.syntax, perl = TRUE)
  # remove any occurrence of >= 2 consecutive newlines to eliminate
  # blank statements; this retains a blank newline at the beginning,
  # if such exists, but parser will not choke because of start.idx
  model.syntax <- gsub("\n{2,}", "\n", model.syntax, perl = TRUE)

  # replace 'strange' tildes (in some locales) (new in 0.6-6)
  model.syntax <- gsub(pattern = "\u02dc", replacement = "~", model.syntax)

  # break up in lines
  model <- unlist(strsplit(model.syntax, "\n"))

  # check for multi-line formulas: they contain no operator symbol
  # but before we do that, we remove all strings between double quotes
  # to avoid confusion with for example equal("f1=~x1") statements
  # model.simple <- gsub("\\(.*\\)\\*", "MODIFIER*", model)
  model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)

  # start.idx <- grep("[~=<>:|%]", model.simple)
  operators <- c(
    "=~", "<~", "~*~", "~~", "~", "==", "<", ">", ":=",
    ":", "\\|", "%"
  )
  lhs.modifiers <- c("efa")
  operators.extra <- c(operators, lhs.modifiers)
  start.idx <- grep(paste(operators.extra, collapse = "|"), model.simple)

  # check for empty start.idx: no operator found (new in 0.6-1)
  if (length(start.idx) == 0L) {
    stop("lavaan ERROR: model does not contain lavaan syntax (no operators found)")
  }

  # check for lonely lhs modifiers (only efa() for now):
  # if found, remove the following start.idx
  efa.idx <- grep("efa\\(", model.simple)
  op.idx <- grep(paste(operators, collapse = "|"), model.simple)
  both.idx <- which(efa.idx %in% op.idx)
  if (length(both.idx) > 0L) {
    efa.idx <- efa.idx[-which(efa.idx %in% op.idx)]
  }
  if (length(efa.idx) > 0L) {
    start.idx <- start.idx[-(match(efa.idx, start.idx) + 1L)]
  }

  # check for non-empty string, without an operator in the first lines
  # (new in 0.6-1)
  if (start.idx[1] > 1L) {
    # two possibilities:
    # - we have an empty line (ok)
    # - the element contains no operator (warn!)
    for (el in 1:(start.idx[1] - 1L)) {
      # not empty?
      if (nchar(model.simple[el]) > 0L) {
        warning("lavaan WARNING: no operator found in this syntax line: ", model.simple[el], "\n", "                  This syntax line will be ignored!")
      }
    }
  }

  end.idx <- c(start.idx[-1] - 1, length(model))
  model.orig <- model
  model <- character(length(start.idx))
  for (i in 1:length(start.idx)) {
    model[i] <- paste(model.orig[start.idx[i]:end.idx[i]], collapse = "")
  }

  # ok, in all remaining lines, we should have an operator outside the ""
  model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)
  idx.wrong <- which(!grepl(
    paste(operators, collapse = "|"),
    model.simple
  ))
  # idx.wrong <- which(!grepl("[~=<>:|%]", model.simple))
  if (length(idx.wrong) > 0) {
    cat("lavaan: missing operator in formula(s):\n")
    print(model[idx.wrong])
    stop("lavaan ERROR: syntax error in lavaan model syntax")
  }

  # but perhaps we have a '+' as the first character?
  idx.wrong <- which(grepl("^\\+", model))
  if (length(idx.wrong) > 0) {
    cat("lavaan: some formula(s) start with a plus (+) sign:\n")
    print(model[idx.wrong])
    stop("lavaan ERROR: syntax error in lavaan model syntax")
  }


  # main operation: flatten formulas into single bivariate pieces
  # with a left-hand-side (lhs), an operator (eg "=~"), and a
  # right-hand-side (rhs)
  # both lhs and rhs can have a modifier
  FLAT.lhs <- character(0)
  FLAT.op <- character(0)
  FLAT.rhs <- character(0)
  FLAT.rhs.mod.idx <- integer(0)
  FLAT.block <- integer(0) # keep track of groups using ":" operator

  FLAT.fixed <- character(0) # only for display purposes!
  FLAT.start <- character(0) # only for display purposes!
  FLAT.lower <- character(0) # only for display purposes!
  FLAT.upper <- character(0) # only for display purposes!
  FLAT.label <- character(0) # only for display purposes!
  FLAT.prior <- character(0)
  FLAT.efa <- character(0)
  FLAT.rv <- character(0)
  FLAT.idx <- 0L
  MOD.idx <- 0L
  CON.idx <- 0L
  MOD <- vector("list", length = 0L)
  CON <- vector("list", length = 0L)
  BLOCK <- 1L
  BLOCK_OP <- FALSE
  for (i in 1:length(model)) {
    x <- model[i]
    if (debug) {
      cat("formula to parse:\n")
      print(x)
      cat("\n")
    }

    # 1. which operator is used?
    line.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", x)
    # "=~" operator?
    if (grepl("=~", line.simple, fixed = TRUE)) {
      op <- "=~"
      # "<~" operator?
    } else if (grepl("<~", line.simple, fixed = TRUE)) {
      op <- "<~"
    } else if (grepl("~*~", line.simple, fixed = TRUE)) {
      op <- "~*~"
      # "~~" operator?
    } else if (grepl("~~", line.simple, fixed = TRUE)) {
      op <- "~~"
      # "~" operator?
    } else if (grepl("~", line.simple, fixed = TRUE)) {
      op <- "~"
      # "==" operator?
    } else if (grepl("==", line.simple, fixed = TRUE)) {
      op <- "=="
      # "<" operator?
    } else if (grepl("<", line.simple, fixed = TRUE)) {
      op <- "<"
      # ">" operator?
    } else if (grepl(">", line.simple, fixed = TRUE)) {
      op <- ">"
      # ":=" operator?
    } else if (grepl(":=", line.simple, fixed = TRUE)) {
      op <- ":="
      # ":" operator?
    } else if (grepl(":", line.simple, fixed = TRUE)) {
      op <- ":"
      # "|" operator?
    } else if (grepl("|", line.simple, fixed = TRUE)) {
      op <- "|"
      # "%" operator?
    } else if (grepl("%", line.simple, fixed = TRUE)) {
      op <- "%"
    } else {
      stop("unknown operator in ", model[i])
    }

    # 2. split by operator (only the *first* occurence!)
    # check first if equal/label modifier has been used on the LEFT!
    if (substr(x, 1, 6) == "label(") {
      stop("label modifier can not be used on the left-hand side of the operator")
    }
    if (op == "|") {
      op.idx <- regexpr("\\|", x)
    } else if (op == "~*~") {
      op.idx <- regexpr("~\\*~", x)
    } else {
      op.idx <- regexpr(op, x)
    }
    lhs <- substr(x, 1L, op.idx - 1L)

    # right-hand side string
    rhs <- substr(x, op.idx + attr(op.idx, "match.length"), nchar(x))

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
      CON.idx <- CON.idx + 1L
      CON[[CON.idx]] <- list(op = op, lhs = lhs, rhs = rhs, user = 1L)
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
      lhs.orig <- lhs
      lhs <- tolower(lhs)
      if (!lhs %in% c("group", "level", "block", "class")) {
        lav_msg_stop(gettextf(
          "unknown block identifier: %s. Block identifier should be
           group, level or block.", dQuote(lhs.orig))
        )
      }

      FLAT.idx <- FLAT.idx + 1L
      FLAT.lhs[FLAT.idx] <- lhs
      FLAT.op[FLAT.idx] <- op
      FLAT.rhs[FLAT.idx] <- rhs
      FLAT.fixed[FLAT.idx] <- ""
      FLAT.start[FLAT.idx] <- ""
      FLAT.lower[FLAT.idx] <- ""
      FLAT.upper[FLAT.idx] <- ""
      FLAT.label[FLAT.idx] <- ""
      FLAT.prior[FLAT.idx] <- ""
      FLAT.efa[FLAT.idx] <- ""
      FLAT.rv[FLAT.idx] <- ""
      FLAT.rhs.mod.idx[FLAT.idx] <- 0L
      if (BLOCK_OP) {
        BLOCK <- BLOCK + 1L
      }
      FLAT.block[FLAT.idx] <- BLOCK
      BLOCK_OP <- TRUE
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
    LHS <- strsplit(lhs, split = "+", fixed = TRUE)[[1]]
    # remove modifiers
    LHS <- gsub("^\\S*\\*", "", LHS)
    if (!all(make.names(LHS) == LHS)) {
      lav_msg_stop(gettextf(
        "lavaan ERROR: left hand side (lhs) of this formula: %1$s %2$s %3$s
        contains either a reserved word (in R) or an illegal character: %4$s.
        See ?reserved for a list of reserved words in R. Please use a variable
        name that is not a reserved word in R and use only characters, digits,
        or the dot symbol.",
        lhs, op, rhs, dQuote(LHS[!make.names(LHS) == LHS])
      ))
    }

    lhs.formula <- as.formula(paste("~", lhs))
    lhs.out <- lav_syntax_parse_rhs(rhs = lhs.formula[[2L]], op = op)
    lhs.names <- names(lhs.out)


    # new in 0.6-4
    # handle LHS modifiers (if any)
    # if(sum(sapply(lhs.out, length)) > 0L) {
    # warning("lavaan WARNING: left-hand side of formula below contains modifier:\n", x,"\n")
    # }

    # 4. lav_syntax_parse_rhs (as rhs of a single-sided formula)

    # new 0.5-12: before we do this, replace '0.2?' by 'start(0.2)*'
    # requested by the simsem folks
    rhs <- gsub("\\(?([-]?[0-9]*\\.?[0-9]*)\\)?\\?", "start(\\1)\\*", rhs)




    # new in 0.6-6, check for rhs NAMES that are reserved names
    # like in foo =~ in + out
    RHS <- strsplit(rhs, split = "+", fixed = TRUE)[[1]]
    RHS.names <- gsub("^\\S*\\*", "", RHS)
    BAD <- c("if", "else", "repeat", "while", "function", "for", "in")
    if (any(RHS.names %in% c(BAD, "NA"))) { # "NA" added in 0.6-8
      stop(
        "lavaan ERROR: right hand side (rhs) of this formula:\n    ",
        lhs, " ", op, " ", rhs,
        "\n    contains either a reserved word (in R) or an illegal character: ",
        dQuote(RHS.names[which(RHS.names %in% BAD)[1]]),
        "\n    See ?reserved for a list of reserved words in R",
        "\n    Please use a variable name that is not a reserved word in R",
        "\n    and use only characters, digits, or the dot symbol."
      )
    }
    # new in 0.6-6, check for rhs LABELS that are reserved names
    # like in foo =~ in*bar
    RHS <- strsplit(rhs, split = "+", fixed = TRUE)[[1]]
    RHS.labels <- gsub("\\*\\S*$", "", RHS)
    if (any(RHS.labels %in% BAD)) {
      stop(
        "lavaan ERROR: right hand side (rhs) of this formula:\n    ",
        lhs, " ", op, " ", rhs,
        "\n    contains either a reserved word (in R) or an illegal character: ",
        dQuote(RHS.names[which(RHS.labels %in% BAD)[1]]),
        "\n    See ?reserved for a list of reserved words in R",
        "\n    Please use a variable name that is not a reserved word in R",
        "\n    and use only characters, digits, or the dot symbol."
      )
    }
    # new in 0.6-12: check for three-way interaction terms (which we do
    # NOT support)
    if (any(grepl(":", RHS.names))) {
      ncolon <- sapply(gregexpr(":", RHS.names), length)
      if (any(ncolon > 1L)) {
        idx <- which(ncolon > 1L)
        lav_msg_stop(gettext(
        "Three-way or higher-order interaction terms (using
multiple colons) are not supported in the lavaan syntax; please manually
construct the product terms yourself in the data.frame, give them an
appropriate name, and then you can use these interaction variables as any
other (observed) variable in the model syntax. Problematic term is: "),
          RHS.names[idx[1]]
        )
      }
    }


    rhs.formula <- as.formula(paste("~", rhs))
    out <- lav_syntax_parse_rhs(rhs = rhs.formula[[2L]], op = op)

    if (debug) print(out)

    # for each lhs element
    for (l in 1:length(lhs.names)) {
      # for each rhs element
      for (j in 1:length(out)) {
        # catch intercepts
        if (names(out)[j] == "intercept") {
          if (op == "~") {
            rhs.name <- ""
          } else {
            # either number (1), or reserved name?
            stop("lavaan ERROR: right-hand side of formula contains an invalid variable name:\n    ", x)
          }
        } else if (names(out)[j] == "..zero.." && op == "~") {
          rhs.name <- ""
        } else if (names(out)[j] == "..constant.." && op == "~") {
          rhs.name <- ""
        } else {
          rhs.name <- names(out)[j]
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
        if (op == "=~" && lhs.names[l] == names(out)[j]) {
          stop("lavaan ERROR: latent variable `", lhs.names[l], "' can not be measured by itself")
        }

        # check if we not already have this combination (in this group)
        # 1. asymmetric (=~, ~, ~1)
        if (op != "~~") {
          idx <- which(FLAT.lhs == lhs.names[l] &
            FLAT.op == op &
            FLAT.block == BLOCK &
            FLAT.rhs == rhs.name)
          if (length(idx) > 0L) {
            stop("lavaan ERROR: duplicate model element in: ", model[i])
          }
        } else {
          # 2. symmetric (~~)
          idx <- which(FLAT.lhs == rhs.name &
            FLAT.op == "~~" &
            FLAT.block == BLOCK &
            FLAT.rhs == lhs.names[l])
          if (length(idx) > 0L) {
            stop("lavaan ERROR: duplicate model element in: ", model[i])
          }
        }

        # check if we have a self-loop (y ~ y)
        if (op %in% c("~", "<~") && rhs.name == lhs.names[l]) {
          # stop("lavaan ERROR: lhs and rhs are the same in: ",
          #     model[i])
          # this breaks pompom package, example uSEM
          warning(
            "lavaan WARNING: lhs and rhs are the same in: ",
            model[i]
          )
        }


        FLAT.idx <- FLAT.idx + 1L
        FLAT.lhs[FLAT.idx] <- lhs.names[l]
        FLAT.op[FLAT.idx] <- op
        FLAT.rhs[FLAT.idx] <- rhs.name
        FLAT.block[FLAT.idx] <- BLOCK
        FLAT.fixed[FLAT.idx] <- ""
        FLAT.start[FLAT.idx] <- ""
        FLAT.label[FLAT.idx] <- ""
        FLAT.lower[FLAT.idx] <- ""
        FLAT.upper[FLAT.idx] <- ""
        FLAT.prior[FLAT.idx] <- ""
        FLAT.efa[FLAT.idx] <- ""
        FLAT.rv[FLAT.idx] <- ""

        mod <- list()
        rhs.mod <- 0L
        if (length(lhs.out[[l]]$efa) > 0L) {
          mod$efa <- lhs.out[[l]]$efa
          FLAT.efa[FLAT.idx] <- paste(mod$efa, collapse = ";")
          rhs.mod <- 1L # despite being a LHS modifier
        }
        if (length(out[[j]]$fixed) > 0L) {
          mod$fixed <- out[[j]]$fixed
          FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$start) > 0L) {
          mod$start <- out[[j]]$start
          FLAT.start[FLAT.idx] <- paste(mod$start, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$lower) > 0L) {
          mod$lower <- out[[j]]$lower
          FLAT.lower[FLAT.idx] <- paste(mod$lower, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$upper) > 0L) {
          mod$upper <- out[[j]]$upper
          FLAT.upper[FLAT.idx] <- paste(mod$upper, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$label) > 0L) {
          mod$label <- out[[j]]$label
          FLAT.label[FLAT.idx] <- paste(mod$label, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$rv) > 0L) {
          mod$rv <- out[[j]]$rv
          FLAT.rv[FLAT.idx] <- paste(mod$rv, collapse = ";")
          rhs.mod <- 1L
        }
        if (length(out[[j]]$prior) > 0L) {
          mod$prior <- out[[j]]$prior
          FLAT.prior[FLAT.idx] <- paste(mod$prior, collapse = ";")
          rhs.mod <- 1L
        }
        # if(op == "~1" && rhs == "0") {
        #    mod$fixed <- 0
        #    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
        #    rhs.mod <- 1L
        # }
        if (op == "=~" && rhs == "0") {
          mod$fixed <- 0
          FLAT.rhs[FLAT.idx] <- FLAT.lhs[FLAT.idx]
          FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse = ";")
          rhs.mod <- 1L
        }

        FLAT.rhs.mod.idx[FLAT.idx] <- rhs.mod

        if (rhs.mod > 0L) {
          MOD.idx <- MOD.idx + 1L
          MOD[[MOD.idx]] <- mod
        }
      } # rhs elements
    } # lhs elements
  } # model elements

  # enumerate modifier indices
  mod.idx <- which(FLAT.rhs.mod.idx > 0L)
  FLAT.rhs.mod.idx[mod.idx] <- 1:length(mod.idx)

  FLAT <- list(
    lhs = FLAT.lhs, op = FLAT.op, rhs = FLAT.rhs,
    mod.idx = FLAT.rhs.mod.idx, block = FLAT.block,
    fixed = FLAT.fixed, start = FLAT.start,
    lower = FLAT.lower, upper = FLAT.upper,
    label = FLAT.label, prior = FLAT.prior,
    efa = FLAT.efa, rv = FLAT.rv
  )

  # change op for intercepts (for convenience only)
  int.idx <- which(FLAT$op == "~" & FLAT$rhs == "")
  if (length(int.idx) > 0L) {
    FLAT$op[int.idx] <- "~1"
  }

  # new in 0.6, reorder covariances here!
  FLAT <- lav_partable_covariance_reorder(FLAT)

  if (as.data.frame.) {
    FLAT <- as.data.frame(FLAT, stringsAsFactors = FALSE)
  }

  # new in 0.6-4: check for 'group' within 'level'
  if (any(FLAT$op == ":")) {
    op.idx <- which(FLAT$op == ":")
    if (length(op.idx) < 2L) {
      # only 1 block identifier? this is weird -> give warning
      warning("lavaan WARNING: syntax contains only a single block identifier: ", FLAT$lhs[op.idx])
    } else {
      first.block <- FLAT$lhs[op.idx[1L]]
      second.block <- FLAT$lhs[op.idx[2L]]
      if (first.block == "level" &&
        second.block == "group") {
        stop("lavaan ERROR: groups can not be nested within levels")
      }
    }
  }

  attr(FLAT, "modifiers") <- MOD
  attr(FLAT, "constraints") <- CON

  FLAT
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
      NAME <- all.vars(rhs)
      if (length(NAME) > 0L) {
        names(out)[1L] <- NAME
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
      NAME <- all.vars(rhs[[3L]])

      if (length(NAME) > 0L) { # not an intercept
        # catch interaction term
        rhs3.names <- all.names(rhs[[3L]])
        if (rhs3.names[1L] == ":") {
          if (length(NAME) == 1) {
            NAME <- paste(NAME[1L], ":", NAME[1L], sep = "")
          } else {
            NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
          }
        }
        names(out)[1L] <- NAME
      } else { # intercept
        names(out)[1L] <- "intercept"
      }
      i.var <- all.vars(rhs[[2L]], unique = FALSE)
      if (length(i.var) > 0L) {
        # modifier are unquoted labels
        out[[1L]]$label <- i.var
      } else {
        # modifer is something else
        out[[1L]] <- lav_syntax_get_modifier(rhs[[2L]])
      }
      break
    } else if (rhs[[1L]] == ":") { # last one, but interaction term
      out <- c(vector("list", 1L), out)
      NAME <- all.vars(rhs)
      if (length(NAME) == 1) {
        NAME <- paste(NAME[1L], ":", NAME[1L], sep = "")
      } else {
        NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
      }
      names(out)[1L] <- NAME
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
        NAME <- all.vars(rhs[[3L]][[3]])

        if (length(NAME) > 0L) { # not an intercept
          # catch interaction term
          rhs3.names <- all.names(rhs[[3L]][[3]])
          if (rhs3.names[1L] == ":") {
            if (length(NAME) == 1) {
              NAME <- paste(NAME[1L], ":", NAME[1L], sep = "")
            } else {
              NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
            }
          }
          names(out)[1L] <- NAME
        } else { # intercept
          names(out)[1L] <- "intercept"
        }
        i.var <- all.vars(rhs[[3]][[2L]], unique = FALSE)
        if (length(i.var) > 0L) {
          # modifier are unquoted labels
          out[[1L]]$label <- i.var
        } else {
          # modifer is something else
          out[[1L]] <- lav_syntax_get_modifier(rhs[[3]][[2L]])
        }

        # interaction term?
      } else if (length(rhs[[3L]]) == 3L && rhs[[3L]][[1]] == ":") {
        # interaction term, without modifier
        NAME <- all.vars(rhs[[3L]])
        if (length(NAME) == 1) {
          NAME <- paste(NAME[1L], ":", NAME[1L], sep = "")
        } else {
          NAME <- paste(NAME[1L], ":", NAME[2L], sep = "")
        }
        names(out)[1L] <- NAME
      } else { # no modifier!!
        NAME <- all.vars(rhs[[3]])
        if (length(NAME) > 0L) {
          names(out)[1L] <- NAME
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
      stop("lavaan ERROR: I'm confused parsing this line: ", rhs, "\n")
    }
  }

  # if multiple elements, check for duplicated elements and merge if found
  if (length(out) > 1L) {
    rhs.names <- names(out)
    while (!is.na(idx <- which(duplicated(rhs.names))[1L])) {
      dup.name <- rhs.names[idx]
      orig.idx <- match(dup.name, rhs.names)
      merged <- c(out[[orig.idx]], out[[idx]])
      if (!is.null(merged)) { # be careful, NULL will delete element
        out[[orig.idx]] <- merged
      }
      out <- out[-idx]
      rhs.names <- names(out)
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
      return(list(label = mod))
    }
  } else if (mod[[1L]] == "start") {
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    return(list(start = cof))
  } else if (mod[[1L]] == "lower") {
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    return(list(lower = cof))
  } else if (mod[[1L]] == "upper") {
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    return(list(upper = cof))
  } else if (mod[[1L]] == "equal") {
    label <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    return(list(label = label))
  } else if (mod[[1L]] == "label") {
    label <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    label[is.na(label)] <- "" # catch 'NA' elements in a label
    return(list(label = label))
  } else if (mod[[1L]] == "rv") {
    rv <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    if (anyNA(rv)) {
      stop("lavaan ERROR: some rv() labels are NA")
    }
    return(list(rv = rv))
  } else if (mod[[1L]] == "prior") {
    prior <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    return(list(prior = prior))
  } else if (mod[[1L]] == "efa") {
    efa <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    return(list(efa = efa))
  } else if (mod[[1L]] == "c") {
    # vector: we allow numeric and character only!
    cof <- unlist(lapply(as.list(mod)[-1],
      eval,
      envir = NULL, enclos = NULL
    ))
    if (all(is.na(cof))) {
      return(list(fixed = rep(as.numeric(NA), length(cof))))
    } else if (is.numeric(cof)) {
      return(list(fixed = cof))
    } else if (is.character(cof)) {
      cof[is.na(cof)] <- "" # catch 'NA' elements in a label
      return(list(label = cof))
    } else {
      stop("lavaan ERROR: can not parse modifier:", mod, "\n")
    }
  } else {
    # unknown expression
    # as a final attempt, we will evaluate it and coerce it
    # to either a numeric or character (vector)
    cof <- try(eval(mod, envir = NULL, enclos = NULL), silent = TRUE)
    if (inherits(cof, "try-error")) {
      stop(
        "lavaan ERROR: evaluating modifier failed: ",
        paste(as.character(mod)[[1]], "()*", sep = ""), "\n"
      )
    } else if (is.numeric(cof)) {
      return(list(fixed = cof))
    } else if (is.character(cof)) {
      return(list(label = cof))
    } else {
      stop(
        "lavaan ERROR: can not parse modifier: ",
        paste(as.character(mod)[[1]], "()*", sep = ""), "\n"
      )
    }
  }
}
