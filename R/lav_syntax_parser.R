# New version of parser, written by Luc De Wilde in september/october 2023
# Adapt to style guide by LDW on April 4, 2026

# ----------------------- lav_create_enum ------------------------------------ #
# function to create an Enumerable like structure in R
#  usage example  mycolors <- lav_create_enum(c("black", "white",
#                                     "orange", "green", "red", "blue"))
#                 xyz <- mycolors$red
# values are default 1L, ..., number of names, but can be user specified
# ---------------------------------------------------------------------------- #
lav_create_enum <- function(names, values = seq_along(names)) {
  stopifnot(identical(unique(names), names), is.character(names))
  stopifnot(length(names) == length(values))
  res <- as.list(setNames(values, names))
  res$enum_names <- names
  res$enum_values <- values
  res$enum_size <- length(values)
  res <- as.environment(res)
  lockEnvironment(res, bindings = TRUE)
  res
}

# ------------------------ lav_parse_sublist --------------------------------- #
# function to create a list with only some indexes for all members
# ---------------------------------------------------------------------------- #
lav_parse_sublist <- function(inlist, indexes) {
  for (j in seq_along(inlist)) {
    inlist[[j]] <- inlist[[j]][indexes]
  }
  inlist
}

# ------------------------ lav_parse_eval_r_expression ------------------- #
# help function to evaluate the value of an r expression formed by the elements
# with index 'from' to 'to' of a formula 'formul1'
# returns "_error_" if evaluation failed
# used only in lav_parse_modifier
# ---------------------------------------------------------------------------- #
lav_parse_eval_r_expression <- function(formul1, from, to, types) {
  strings <- vapply(seq.int(from, to), function(x) {
    if (formul1$elem_type[x] == types$stringliteral) {
      paste0('"', formul1$elem_text[x], '"')
    } else {
      formul1$elem_text[x]
    }
  }, "")
  txt <- paste(strings, collapse = "")
  result <- try(eval(parse(text = txt),
    envir = NULL,
    enclos = baseenv()
  ), silent = TRUE)
  if (inherits(result, "try-error")) {
    return("_error_")
  }
  result
}
# ------------------------ lav_parse_txtloc  --------------------------------- #
# function which translates a position in the model source string to a
# user friendly locator (=[1L]) and the line with position (=[2L])
# ---------------------------------------------------------------------------- #
lav_parse_txtloc <- function(modelsrc, position) {
  txt <- c("", "")
  if (nchar(modelsrc) >= position && position > 0) {
    newlines <- gregexpr("\n", paste0(modelsrc, "\n"), fixed = TRUE)[[1]]
    lijn <- which(newlines >= position)[1]
    if (lijn == 1L) {
      pos <- position
      lijnchar <- substr(modelsrc, 1L, newlines[1])
    } else {
      pos <- position - newlines[lijn - 1L]
      lijnchar <- substr(modelsrc, newlines[lijn - 1L] + 1L, newlines[lijn])
    }
    if (nchar(lijnchar) == 1L) {
      lijnchar <- ""
    } else {
      lijnchar <- substr(lijnchar, 1L, nchar(lijnchar) - 1)
    }
    # adapt line number when first line blank :
    if (grepl("^[ \t]*\n", modelsrc)) lijn <- lijn - 1L
    txt <- c(
      gettextf(" at line %1$s, pos %2$s", lijn, pos),
      paste(lijnchar, "\n", strrep(" ", pos - 1L), "^\n", sep = "")
    )
  }
  txt
}

# ------------------------ lav_parse_text_tokens ----------------------------- #
# function to split the model source in tokens.
# Returns a list with tokens with their attributes
#   elem_pos  : position in source
#   elem_type : type of token (cf. definition of types
#               in lav_parse_model_string)
#   elem_text : the text of the token
#   elem_formule.number : sequence number of the 'logical'
#                         formula where the token occurs
# the function returns the stored tokens in a list
# ---------------------------------------------------------------------------- #
lav_parse_text_tokens <- function(modelsrc, types) {
  nmax <- nchar(modelsrc)
  elem_pos <- vector("integer", nmax)
  elem_type <- elem_pos
  elem_text <- vector("character", nmax)
  elem_i <- 1L
  modelsrcw <- paste0(modelsrc, "\n") # working model, must end
  # with a newline for tests via regexpr
  stringliterals <- gregexpr("\"[^\"]*?[\"\n]", modelsrcw)[[1L]]
  if (stringliterals[1L] > -1L) {
    stringliteral_lengths <- attr(stringliterals, "match.length")
    for (i in seq_along(stringliterals)) {
      pfpos <- stringliterals[i]
      pflen <- stringliteral_lengths[i]
      substr(modelsrcw, pfpos + 1L, pfpos + pflen - 2L) <-
        strrep(" ", pflen - 2L)
      elem_pos[elem_i] <- pfpos
      elem_text[elem_i] <- substr(modelsrc, pfpos + 1L, pfpos + pflen - 2L)
      elem_type[elem_i] <- types$stringliteral
      elem_i <- elem_i + 1L
    }
  }
  comments <- gregexpr("[#!].*?\n", modelsrcw)[[1L]]
  if (comments[1] > -1L) {
    comment_lengths <- attr(comments, "match.length")
    for (i in seq_along(comments)) {
      substr(modelsrcw, comments[i], comments[i] + comment_lengths[i] - 1L) <-
        strrep(" ", comment_lengths[i] - 1L)
      # check for stringliterals in comment
      str.in.comment <- (elem_pos > comments[i] &
        elem_pos < comments[i] + comment_lengths[i])
      if (any(str.in.comment)) {
        elem_type[str.in.comment] <- 0
      }
    }
  }
  modelsrcw <- gsub("\t", " ", modelsrcw)
  newlines <- gregexpr("[;\n]", modelsrcw)[[1L]]
  if (newlines[1L] > -1L) {
    for (i in seq_along(newlines)) {
      pfpos <- newlines[i]
      substr(modelsrcw, pfpos, pfpos) <- "\n"
      elem_pos[elem_i] <- pfpos
      elem_text[elem_i] <- "\n"
      elem_type[elem_i] <- types$newline
      elem_i <- elem_i + 1L
    }
  }
  # --------------------- handling spaces in operators ----------------------- #
  if (grepl("= +~", modelsrcw)) {
    waar <- regexpr("= +~", modelsrcw)[1]
    modelsrcw <- gsub("=( +)~", "=~\\1", modelsrcw)
    tl <- lav_parse_txtloc(modelsrc, waar)
    lav_msg_warn(gettext("splitting of '=~' deprecated"),
      tl[1L],
      footer = tl[2L]
    )
  }
  if (grepl("[^=~]~ +~", modelsrcw)) {
    waar <- regexpr("[^=~]~ +~", modelsrcw)[1]
    modelsrcw <- gsub("([^=~])~( +)~", "\\1~~\\2", modelsrcw)
    tl <- lav_parse_txtloc(modelsrc, waar + 1L)
    lav_msg_warn(gettext("splitting of '~~' deprecated"),
      tl[1L],
      footer = tl[2L]
    )
  }
  # -------------------------------------------------------------------------- #
  lavops <- gregexpr("=~|<~|~\\*~|~~|~|\\|~|==|<|>|:=|:|\\||%", modelsrcw)[[1]]
  if (lavops[1L] > -1L) {
    lavop_lengths <- attr(lavops, "match.length")
    for (i in seq_along(lavops)) {
      pfpos <- lavops[i]
      pflen <- lavop_lengths[i]
      elem_pos[elem_i] <- pfpos
      elem_text[elem_i] <- substr(modelsrcw, pfpos, pfpos + pflen - 1L)
      elem_type[elem_i] <- types$lavaanoperator
      substr(modelsrcw, pfpos, pfpos + pflen - 1L) <- strrep(" ", pflen)
      elem_i <- elem_i + 1L
    }
  }
  symbols <- gregexpr("[,()/*?^']", modelsrcw)[[1L]] # f1=~x2 + 0.5 ? x3
  symbols1 <- gregexpr("[-+][^.0-9]", modelsrcw)[[1L]] # f1=~x2+x3
  symbols2 <- gregexpr("[._0-9a-df-zA-DF-Z)] *[-+][.0-9]", modelsrcw)[[1L]]
  # f1=~x2+2*x3, len-2 !
  symbols3 <- gregexpr("[^.0-9][eE] *[-+][.0-9]", modelsrcw)[[1L]]
  # f1=~xe+2*x3, len-2 !
  if (symbols1[1L] > -1L) {
    if (symbols[1L] == -1L) {
      symbols <- symbols1
    } else {
      symbols <- c(symbols, symbols1)
    }
  }
  if (symbols2[1L] > -1L) {
    symbols2_lengths <- attr(symbols2, "match.length")
    symbols2 <- symbols2 + symbols2_lengths - 2L
    if (symbols[1L] == -1L) {
      symbols <- symbols2
    } else {
      symbols <- c(symbols, symbols2)
    }
  }
  if (symbols3[1L] > -1L) {
    symbols3_lengths <- attr(symbols3, "match.length")
    symbols3 <- symbols3 + symbols3_lengths - 2L
    if (symbols[1L] == -1L) {
      symbols <- symbols3
    } else {
      symbols <- c(symbols, symbols3)
    }
  }
  if (symbols[1L] > -1L) {
    for (i in seq_along(symbols)) {
      pfpos <- symbols[i]
      substr(modelsrcw, pfpos, pfpos) <- " "
      elem_pos[elem_i] <- pfpos
      elem_text[elem_i] <- substr(modelsrc, pfpos, pfpos)
      elem_type[elem_i] <- types$symbol
      elem_i <- elem_i + 1L
    }
  }

  numliterals <- gregexpr(
    "([ \n][-+][.0-9]|[ \n]\\.[0-9]|[ \n][0-9])[-+\\.0-9eE]*",
    paste0(" ", modelsrcw)
  )[[1]]
  if (numliterals[1L] > -1L) {
    numliteral_lengths <- attr(numliterals, "match.length") - 1L
    for (i in seq_along(numliterals)) {
      pfpos <- numliterals[i]
      pflen <- numliteral_lengths[i]
      substr(modelsrcw, pfpos, pfpos + pflen - 1L) <- strrep(" ", pflen)
      elem_pos[elem_i] <- pfpos
      elem_text[elem_i] <- substr(modelsrc, pfpos, pfpos + pflen - 1L)
      elem_type[elem_i] <- types$numliteral
      elem_i <- elem_i + 1L
    }
  }
  identifiers <- gregexpr("[ \n][_.[:alpha:]][_.[:alnum:]]*",
                 paste0(" ", modelsrcw)
  )[[1]]
  identifier_lengths <- attr(identifiers, "match.length") - 1L
  for (i in seq_along(identifiers)) {
    pfpos <- identifiers[i]
    pflen <- identifier_lengths[i]
    substr(modelsrcw, pfpos, pfpos + pflen - 1L) <- strrep(" ", pflen)
    elem_pos[elem_i] <- pfpos
    elem_text[elem_i] <- substr(modelsrc, pfpos, pfpos + pflen - 1L)
    elem_type[elem_i] <- types$identifier
    elem_i <- elem_i + 1L
  }
  # check for uninterpreted chars
  wrong <- regexpr("[^\"\n ]", modelsrcw)
  if (wrong != -1L) {
    tl <- lav_parse_txtloc(modelsrc, wrong)
    lav_msg_stop(gettext("unexpected character"),
      tl[1L],
      footer = tl[2L]
    )
  }
  # remove unused elements from vectors
  elements <- which(elem_type > 0L)
  elem_pos <- elem_pos[elements]
  elem_type <- elem_type[elements]
  elem_text <- elem_text[elements]
  # order tokens
  token_order <- order(elem_pos)
  elem_pos <- elem_pos[token_order]
  elem_type <- elem_type[token_order]
  elem_text <- elem_text[token_order]

  # concatenate identifiers with only spaces in between - LDW 22/4/2024
  elem_i <- length(elem_pos)
  concatenated <- FALSE
  while (elem_i > 1L) {
    if (any(elem_type[elem_i] == c(types$identifier, types$numliteral)) &&
      elem_type[elem_i - 1L] == types$identifier) {
      spaces_between <- elem_pos[elem_i] - elem_pos[elem_i - 1L] -
        length(elem_text[elem_i - 1L])
      elem_text[elem_i - 1L] <- paste0(
        elem_text[elem_i - 1L],
        strrep(" ", spaces_between),
        elem_text[elem_i]
      )
      elem_type[elem_i] <- 0L
      concatenated <- TRUE
    }
    elem_i <- elem_i - 1L
  }
  if (concatenated) { # remove items with type 0
    elements <- which(elem_type > 0L)
    elem_pos <- elem_pos[elements]
    elem_type <- elem_type[elements]
    elem_text <- elem_text[elements]
  }

  # to set formula number
  elem_formula_number <- rep(0L, length(elem_type))
  frm_number <- 1L
  frm_hasefa <- FALSE
  frm_lastplus <- FALSE
  frm_incremented <- FALSE
  for (i in seq_along(elem_type)) {
    elem_formula_number[i] <- frm_number
    if (elem_type[i] == types$identifier && elem_text[i] == "efa") {
      frm_hasefa <- TRUE
    }
    if (any(elem_text[i] ==
      c("+", "*", "=~", "-", "<~", "~*~", "~~", "~", "|~", "|", "%"))) {
      if (frm_incremented) {
        frm_number <- frm_number - 1L
        elem_formula_number[i] <- frm_number
        frm_incremented <- FALSE
      }
      frm_lastplus <- TRUE
    } else {
      if (any(elem_type[i] == c(
        types$stringliteral, types$identifier, types$numliteral,
        types$stringliteral, types$symbol
      ))) {
        frm_lastplus <- FALSE
      }
      if (i > 1 && elem_type[i] != types$newline &&
        elem_type[i - 1L] == types$lavaanoperator) {
        frm_hasefa <- FALSE
      }
    }
    if (elem_type[i] == types$newline) {
      if (i > 1 && elem_type[i - 1L] != types$newline) { # ignore multiple nl's
        if (!frm_hasefa && !frm_lastplus) {
          frm_number <- frm_number + 1L
          frm_incremented <- TRUE
        } else {
          frm_hasefa <- FALSE
        }
      }
    } else {
      frm_incremented <- FALSE
    }
  }
  list(
    elem_pos = elem_pos, elem_type = elem_type,
    elem_text = elem_text, elem_formula_number = elem_formula_number
  )
}

# ------------------------ lav_parse_tokens_formulas ------------------------- #
# function to group the modellist tokens in 'mono' formulas.
# mono means that the terms (for formulas other then blocks and constraints)
#   are split in seperate formula's, e.g.
#   a1 + a2 =~ b1 + b2  becomes
#     /  a1 =~ b1
#     |  a1 =~ b2
#     |  a2 =~ b1
#     \  a2 =~ b2
# newlines are removed
# the function returns a list of formulas
# ---------------------------------------------------------------------------- #
lav_parse_tokens_formulas <- function(modellist, modelsrc, types) {
  real_operators <- c("=~", "<~", "~*~", "~~", "~", "|~", "|", "%")
  welke <- modellist$elem_type != types$newline
  formula_numbers <- unique(modellist$elem_formula_number[welke])
  formulas <- lapply(formula_numbers, function(s) {
    welkenu <- modellist$elem_formula_number == s & welke
    list(
      elem_pos = modellist$elem_pos[welkenu],
      elem_type = modellist$elem_type[welkenu],
      elem_text = modellist$elem_text[welkenu]
    )
  })
  maxnum <- length(formula_numbers) + sum(modellist$elem_text == "+")
  outval <- vector(mode = "list", length = maxnum)
  realnum <- 0L
  for (i in seq_along(formulas)) {
    formul1 <- formulas[[i]]
    opi <- which(formul1$elem_type == types$lavaanoperator)
    nelem <- length(formul1$elem_type)
    if (length(opi) == 0L) {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[1])
      lav_msg_stop(gettext("formula without valid operator"),
        tl[1L],
        footer = tl[2L]
      )
    }
    if (length(opi) > 1L) opi <- opi[1] # only first operator taken
    if (any(formul1$elem_text[opi] == real_operators) &&
      sum(formul1$elem_text == "+") > 0) {
      # check + symbols outside parentheses in left and right hand side
      lhplusjes <- integer(0)
      openparentheses <- 0L
      for (jj in seq.int(1L, opi - 1L)) {
        if (formul1$elem_text[jj] == "(") {
          openparentheses <- openparentheses + 1L
          next
        }
        if (formul1$elem_text[jj] == ")") {
          openparentheses <- openparentheses - 1L
          next
        }
        if (formul1$elem_text[jj] == "+" && openparentheses == 0L) {
          lhplusjes <- c(lhplusjes, jj)
        }
      }
      lhplusjes <- c(lhplusjes, opi)
      plusjes <- integer(0)
      openparentheses <- 0L
      for (jj in seq.int(opi + 1L, nelem)) {
        if (formul1$elem_text[jj] == "(") {
          openparentheses <- openparentheses + 1L
          next
        }
        if (formul1$elem_text[jj] == ")") {
          openparentheses <- openparentheses - 1L
          next
        }
        if (formul1$elem_text[jj] == "+" && openparentheses == 0L) {
          plusjes <- c(plusjes, jj)
        }
      }
      plusjes <- c(plusjes, nelem + 1)
      # splitting lhs and rhs on '+' signs
      for (j in seq_along(lhplusjes)) {
        j0 <- 1L
        if (j > 1L) j0 <- lhplusjes[j - 1L] + 1L
        j1 <- lhplusjes[j] - 1L
        if (j1 < j0) next # skip empty parts
        for (k in seq_along(plusjes)) {
          k0 <- opi + 1L
          k1 <- plusjes[k] - 1L
          if (k > 1L) k0 <- plusjes[k - 1L] + 1L
          if (k1 < k0) next # skip empty parts
          welke <- c(seq.int(j0, j1), opi, seq.int(k0, k1))
          realnum <- realnum + 1L
          outval[[realnum]] <- lav_parse_sublist(formul1, welke)
        }
      }
    } else {
      realnum <- realnum + 1L
      outval[[realnum]] <- formul1
    }
  }
  outval[seq_len(realnum)]
}
# ------------------------ lav_parse_check_name ------------------------ #
# checks if an element of the elem_text member in a list is a valid r-name
# ---------------------------------------------------------------------------- #
lav_parse_check_name <- function(formul1, ind, modelsrc) {
  # allow spaces, LDW 22/4/2024
  testitem <- gsub(" ", "_", formul1$elem_text[ind], fixed = TRUE)

  if (make.names(testitem) != testitem) {
    tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[ind])
    lav_msg_stop(
      gettext("identifier is either a reserved word (in R) or
              contains an illegal character"),
      tl[1L],
      footer = tl[2L]
    )
  }
  invisible(NULL)
}
# ------------------------ lav_parse_modifier ---------------------------- #
# The function takes a list with tokens belonging to a single 'mono' lavaan
# formula as input. The other arguments are:
#       lhs : check for lhs or rhs modifier
#       opi : index of the lavaan operator in the list-items
#  modelsrc : the model source string (for error messages and warnings)
#     types : the types of tokens
#       rme : index of last element of modifier in formula (*)
#   rmeprev : index of first element of modifier in formula - 1L (*)
# The function return the modifier detected as element of a list
#  with name the modifier type (efa, fixed, start, label, lower, upper, prior or
#   rv) and value an array of values (length > 1 if vector via c(...)) for the
#   modifier value.
# (*) if rme > remprev the rhs is limited to the elements with index
#     rmeprev+1:rme, this is to support multiple modifiers for the same element.
# An error message is produced when no modifier can be determined.
# ---------------------------------------------------------------------------- #
lav_parse_modifier <- function(formul1, lhs, opi, modelsrc, types,
                                   rme = 0L, rmeprev = 0L, modenv) {
  if (rme > rmeprev) {
    welke <- c(seq.int(1L, opi), seq.int(rmeprev + 1L, rme),
               length(formul1$elem_type))
    formul1 <- lav_parse_sublist(formul1, welke)
  }
  nelem <- length(formul1$elem_type)
  # remove unnecessary parentheses (one element between parentheses, previous
  # no identifier)
  check_more <- TRUE
  while (check_more && nelem > 4L) {
    check_more <- FALSE
    for (par.i in seq.int(3L, nelem - 1L)) {
      if (formul1$elem_text[par.i - 1L] == "(" &&
        formul1$elem_text[par.i + 1L] == ")" &&
        formul1$elem_type[par.i - 2L] != types$identifier) {
        formul1$elem_type[par.i - 1L] <- 0L
        formul1$elem_type[par.i + 1L] <- 0L
        check_more <- TRUE
      }
    }
    if (check_more) {
      formul1 <- lav_parse_sublist(formul1, which(formul1$elem_type > 0))
      nelem <- length(formul1$elem_type)
    }
  }
  if (lhs) {
    # modifier on left hand side
    # only 1 possibility : efa ( expression-resulting-in-char ) *
    #                                        identifier operator ... (rhs) ...
    if (formul1$elem_text[1L] == "efa" &&
      formul1$elem_text[2L] == "(" &&
      formul1$elem_text[opi - 3L] == ")" &&
      formul1$elem_text[opi - 2L] == "*") {
      temp <- lav_parse_eval_r_expression(formul1, 3L, opi - 4L, types)
      if (is.character(temp) && temp[1] != "_error_") {
        return(list(efa = temp))
      }
    }
    tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[1L])
    lav_msg_stop(gettext("invalid left hand side modifier"),
      tl[1L],
      footer = tl[2L]
    )
  } else {
    getmodifier <- function(s) {
      if (s %in% c("c", "t")) {
        tmp <- s
        attr(tmp, "tiepe") <- "label"
        return(tmp)
      }
      v <- all.vars(str2expression(s))
      for (v1 in v) {
        assign(v1, v1, modenv)
      }
      tmp <- eval(str2expression(s), modenv)
      rm(list = v, pos = modenv)
      if (is.null(attr(tmp, "tiepe"))) {
        if (is.numeric(tmp) || is.logical(tmp)) {
          attr(tmp, "tiepe") <- "fixed"
        } else {
          attr(tmp, "tiepe") <- "label"
        }
      }
      tmp
    }
    strings <- vapply(seq.int(opi + 1L, nelem - 2L), function(x) {
      if (formul1$elem_type[x] == types$stringliteral) {
        paste0('"', formul1$elem_text[x], '"')
      } else {
        formul1$elem_text[x]
      }
    }, "")
    txt <- paste(strings, collapse = "")
    if (formul1$elem_text[nelem - 1L] == "?") {
      formul1$elem_text[nelem - 1L] <- "*"
      txt <- paste0("start(", txt, ")")
    }
    if (formul1$elem_text[nelem - 1L] != "*") {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[nelem - 1L])
      lav_msg_stop(gettext("invalid modifier symbol (should be '*' or '?')"),
        tl[1L],
        footer = tl[2L]
      )
    }
    modifier <- try(getmodifier(txt), silent = TRUE)
    if (inherits(modifier, "try-error")) {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[opi + 1L])
      lav_msg_stop(gettext("invalid modifier specification"),
        tl[1L],
        footer = tl[2L]
      )
    }
    if (attr(modifier, "tiepe") == "label") {
      if (!is.character(modifier)) {
        tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[opi + 1L])
        lav_msg_stop(gettext("invalid label modifier (should be character)"),
          tl[1L],
          footer = tl[2L]
        )
      }
    } else {
      if (!is.numeric(modifier) && !all(is.na(modifier))) {
        tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[opi + 1L])
        lav_msg_stop(gettext("invalid numeric modifier"),
          tl[1L],
          footer = tl[2L]
        )
      }
    }
    modifierlist <- list(as.vector(modifier))
    names(modifierlist) <- attr(modifier, "tiepe")
    modifierlist
  }
}

# -------------------- main parsing function --------------------------------- #
lav_parse_model_string <- function(model_syntax = "", as_data_frame = FALSE) {
  stopifnot(length(model_syntax) > 0L)
  # replace 'strange' tildes (in some locales) (new in 0.6-6)
  modelsrc <- gsub(
    pattern = "\u02dc",
    replacement = "~",
    paste(unlist(model_syntax), "", collapse = "\n")
  )
  hashstring <- paste0("mdl_", lav_char2hash(paste0(modelsrc, as_data_frame)))
  if (exists(hashstring, envir = lavaan_cache_env)) {
    return(get(hashstring, envir = lavaan_cache_env))
  }
  modenv <- new.env()
assign("label", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "label"
  x
}, modenv)
assign("start", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "start"
  x
}, modenv)
assign("fixed", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "fixed"
  x
}, modenv)
assign("upper", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "upper"
  x
}, modenv)
assign("lower", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "lower"
  x
}, modenv)
assign("rv", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "rv"
  x
}, modenv)
assign("prior", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "prior"
  x
}, modenv)
assign("equal", function(...) {
  x <- do.call("c", args = list(...))
  attr(x, "tiepe") <- "label"
  x
}, modenv)
  types <- lav_create_enum(c(
    "identifier", "numliteral", "stringliteral",
    "symbol", "lavaanoperator", "newline"
  ))
  modellist <- lav_parse_text_tokens(modelsrc, types)
  if (lav_debug()) {
    print(data.frame(
      pos = modellist$elem_pos,
      type = types$enum_names[modellist$elem_type],
      text = modellist$elem_text,
      formula = modellist$elem_formula_number
    ))
  }
  formulalist <- lav_parse_tokens_formulas(modellist, modelsrc, types)
  #---- analyse syntax formulas and put in flat_-----
  max.mono.formulas <- length(formulalist)
  flat_lhs <- character(max.mono.formulas)
  flat_op <- character(max.mono.formulas)
  flat_rhs <- character(max.mono.formulas)
  flat_rhs_mod_idx <- integer(max.mono.formulas)
  flat_block <- integer(max.mono.formulas) # keep track of groups using ":" opr
  flat_fixed <- character(max.mono.formulas) # only for display purposes!
  flat_start <- character(max.mono.formulas) # only for display purposes!
  flat_lower <- character(max.mono.formulas) # only for display purposes!
  flat_upper <- character(max.mono.formulas) # only for display purposes!
  flat_label <- character(max.mono.formulas) # only for display purposes!
  flat_prior <- character(max.mono.formulas)
  flat_efa <- character(max.mono.formulas)
  flat_rv <- character(max.mono.formulas)
  flat_idx <- 0L
  mod_idx <- 0L
  constraints <- list()
  mod <- list()
  block <- 1L
  block_op <- FALSE
  if (lav_debug()) {
    cat("formula to analyse:\n")
  }
  #  operators <- c("=~", "<~", "~*~", "~~", "~", "|~", "==", "<", ">", ":=",
  #                 ":", "|", "%")
  constraint_operators <- c("==", "<", ">", ":=")
  for (s in seq_along(formulalist)) {
    formul1 <- formulalist[[s]]
    if (lav_debug()) {
      cat(vapply(seq_along(formul1$elem_type), function(j) {
        if (formul1$elem_type[j] == types$stringliteral) {
          return(dQuote(formul1$elem_text[j], FALSE))
        }
        formul1$elem_text[j]
      }, ""), "\n")
    }
    nelem <- length(formul1$elem_type)
    # where is the operator
    opi <- which(formul1$elem_type == types$lavaanoperator)
    if (length(opi) > 1L) { # if more then 1 operator skip operators ':'
      opii <- 1L
      while (formul1$elem_text[opi[opii]] == ":" && opii < length(opi)) {
        opii <- opii + 1L
      }
      opi <- opi[opii]
    }
    op <- formul1$elem_text[opi]
    if (any(op == constraint_operators)) { # ----- constraints -------
      lhs <- paste(formul1$elem_text[seq.int(1L, opi - 1L)], collapse = "")
      rhs <- paste(formul1$elem_text[seq.int(opi + 1L, nelem)], collapse = "")
      constraints <- c(
        constraints,
        list(list(
          op = op,
          lhs = lhs,
          rhs = rhs,
          user = 1L
        ))
      )
      next
    }
    if (op == ":") { # ------------------------- block start ----------------- #
      if (opi == 1L) {
        tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[1])
        lav_msg_stop(
          gettext(
          "Missing block identifier. The correct syntax is: \"LHS: RHS\",
          where LHS is a block identifier (eg group or level), and RHS is
          the group/level/block number or label."
          ), tl[1L], footer = tl[2L]
        )
      }
      if (opi > 2L || all(tolower(formul1$elem_text[1]) !=
        c("group", "level", "block", "class"))) {
        tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[1])
        lav_msg_stop(
          gettext(
          "Invalid block identifier. The correct syntax is: \"LHS: RHS\",
          where LHS is a block identifier (eg group or level), and RHS is
          the group/level/block number or label."
          ),
          tl[1L], footer = tl[2L]
        )
      }
      if (nelem != 3 || all(formul1$elem_type[3] !=
        c(types$stringliteral, types$identifier, types$numliteral))) {
        tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[1])
        lav_msg_stop(
          gettext("syntax contains block identifier \"group\" with missing or
                invalid number/label.The correct syntax is: \"LHS: RHS\", where
                LHS is a block identifier (eg group or level), and RHS is the
                group/level/block number or label."),
          tl[1L],
          footer = tl[2L]
        )
      }
      flat_idx <- flat_idx + 1L
      flat_lhs[flat_idx] <- formul1$elem_text[1]
      flat_op[flat_idx] <- op
      flat_rhs[flat_idx] <- formul1$elem_text[3]
      flat_rhs_mod_idx[flat_idx] <- 0L
      if (block_op) {
        block <- block + 1L
      } else {
        if (flat_idx != 1) {
          tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[1])
          lav_msg_warn(
            gettext("First block defined after other formula's"),
            tl[1L],
            footer = tl[2L]
          )
        }
      }
      flat_block[flat_idx] <- block
      block_op <- TRUE
      next
    }
    # ------------------ relational operators -------------------------------- #
    # warn if some identifiers contain spaces
    contsp <- which(formul1$elem_type == types$identifier &
      grepl(" ", formul1$elem_text, fixed = TRUE))
    if (length(contsp) > 0L) {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[contsp[1L]])
      lav_msg_warn(
        gettextf(
          "having identifiers with spaces ('%s') is deprecated",
          formul1$elem_text[contsp[1]]
        ),
        tl[1L],
        footer = tl[2L]
      )
    }
    # checks for valid names in lhs and rhs
    lav_parse_check_name(formul1, opi - 1L, modelsrc) # valid name lhs
    for (j in seq.int(opi + 1L, nelem)) { # valid names rhs
      if (formul1$elem_type[j] == types$identifier &&
        formul1$elem_text[j] != "NA") {
        lav_parse_check_name(formul1, j, modelsrc)
      }
    }
    if (formul1$elem_type[nelem] != types$identifier &&
        (formul1$elem_type[nelem] != types$numliteral ||
         all(op != c("~", "|~", "=~")))) {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[nelem])
      lav_msg_stop(
        gettext("Last element of rhs part expected to be an identifier or,
                for operator ~, |~ or =~, a numeric literal!"),
        tl[1L],
        footer = tl[2L]
      )
    }
    # intercept fixed on 0
    # replace 'lhs ~ 0' => 'lhs ~ 0 * 1' - intercept fixed on zero
    if (formul1$elem_text[nelem] == "0" && op == "~" && opi == nelem - 1L) {
      formul1$elem_type <- c(formul1$elem_type, types$symbol, types$numliteral)
      formul1$elem_text <- c(formul1$elem_text, "*", "1")
      formul1$elem_pos <- c(formul1$elem_pos, rep(formul1$elem_pos[nelem], 2))
      nelem <- length(formul1$elem_type)
    }
    # phantom latent variable
    # replace 'lhs =~ 0' => 'lhs =~ fixed(0)*lhs', 0 can be other numliteral
    #          also, lhs is last element before '=~'
    if (formul1$elem_type[nelem] == types$numliteral && op == "=~") {
      formul1$elem_type <- c(
        formul1$elem_type[seq.int(1L, nelem - 1L)], types$identifier,
        types$symbol, types$numliteral, types$symbol, types$symbol,
        types$identifier
      )
      formul1$elem_text <- c(
        formul1$elem_text[seq.int(1L, nelem - 1L)], "fixed", "(",
        formul1$elem_text[nelem], ")", "*", formul1$elem_text[opi - 1L]
      )
      formul1$elem_pos <- c(
        formul1$elem_pos[seq.int(1L, nelem - 1L)],
        rep(formul1$elem_pos[nelem], 6)
      )
      nelem <- length(formul1$elem_type)
    }
    # handling interaction variable types
    colons <- which(formul1$elem_text[seq.int(1L, nelem - 1L)] == ":" &
      formul1$elem_type[seq.int(2L, nelem)] == types$identifier)
    # check at most 1 colon
    if (length(colons) > 2L ||
        (length(colons) == 2L && (colons[1L] > opi || colons[2L] < opi))) {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[colons[2]])
      lav_msg_stop(
        gettext(
        "Three-way or higher-order interaction terms (using multiple
         colons) are not supported in the lavaan syntax; please manually
         construct the product terms yourself in the data.frame, give
         them an appropriate name, and then you can use these interaction
         variables as any other (observed) variable in the model syntax."
        ), tl[1L], footer = tl[2L]
      )
    }
    if (length(colons) > 0L) {
      # collapse items around colon "a" ":" "b" => "a:b"
      formul1$elem_text[colons[1L] - 1L] <-
        paste(formul1$elem_text[seq.int(colons[1L] - 1L, colons[1L] + 1L)],
          collapse = ""
        )
      formul1 <- lav_parse_sublist(formul1,
          setdiff(seq.int(1L, nelem), seq.int(colons[1L], colons[1L] + 1L)))
      nelem <- length(formul1$elem_type)
      if (colons[1L] < opi) {
        opi <- opi - 2L # is in LHS
        if (length(colons) == 2L) colons[2L] <- colons[2L] - 2L
      }
    }
    if (length(colons) == 2L) {
      # collapse items around colon "a" ":" "b" => "a:b"
      formul1$elem_text[colons[2L] - 1L] <-
        paste(formul1$elem_text[seq.int(colons[2L] - 1L, colons[2L] + 1L)],
              collapse = ""
        )
      formul1 <- lav_parse_sublist(formul1,
          setdiff(seq.int(1L, nelem), seq.int(colons[2L], colons[2L] + 1L)))
      nelem <- length(formul1$elem_type)
    }
    # modifiers
    rhsmodelems <- which(seq_along(formul1$elem_type) > opi &
      formul1$elem_type == types$symbol &
      (formul1$elem_text == "*" | formul1$elem_text == "?"))
    for (j in seq_along(rhsmodelems)) {
      if (sum(formul1$elem_text[seq.int(opi, rhsmodelems[j])] == "(") !=
          sum(formul1$elem_text[seq.int(opi, rhsmodelems[j])] == ")"))
        rhsmodelems[j] <- 0L
    }
    rhsmodelems <- rhsmodelems[rhsmodelems != 0L]
    if (length(rhsmodelems) == 0L) rhsmodelems <- opi
    lhs <- formul1$elem_text[opi - 1L]
    rhs <- formul1$elem_text[nelem]
    for (rmei in seq_along(rhsmodelems)) {
      rme <- rhsmodelems[rmei]
      rmeprev <- if (rmei == 1L) opi else rhsmodelems[rmei - 1L]
      already <- which(flat_lhs == lhs & flat_op == op & flat_block == block &
        (flat_rhs == rhs | (flat_rhs == "" & (op == "~" | op == "|~") &
          formul1$elem_type[nelem] == types$numliteral)))
      if (length(already) == 1L) {
        idx <- already
      } else {
        flat_idx <- flat_idx + 1L
        idx <- flat_idx
        flat_lhs[idx] <- lhs
        flat_op[idx] <- op
        flat_rhs[idx] <- rhs
        flat_block[idx] <- block
        if (formul1$elem_type[nelem] == types$numliteral) {
          if (op == "~" || op == "|~") flat_rhs[idx] <- ""
        }
      }
      lhsmod <- list()
      if (opi > 2 && rmei == 1L) {
        lhsmod <- lav_parse_modifier(
          formul1,
          TRUE, opi, modelsrc, types, modenv = modenv
        )
      }
      rhsmod <- list()
      if (nelem - opi > 1) {
        rhsmod <- lav_parse_modifier(
          formul1,
          FALSE, opi, modelsrc, types, rme, rmeprev, modenv = modenv
        )
      }
      flat_fixed[idx] <- if (is.null(rhsmod$fixed)) {
        flat_fixed[idx]
      } else {
        paste(rhsmod$fixed, collapse = ";")
      }
      flat_start[idx] <- if (is.null(rhsmod$start)) {
        flat_start[idx]
      } else {
        paste(rhsmod$start, collapse = ";")
      }
      flat_label[idx] <- if (is.null(rhsmod$label)) {
        flat_label[idx]
      } else {
        paste(rhsmod$label, collapse = ";")
      }
      flat_lower[idx] <- if (is.null(rhsmod$lower)) {
        flat_lower[idx]
      } else {
        paste(rhsmod$lower, collapse = ";")
      }
      flat_upper[idx] <- if (is.null(rhsmod$upper)) {
        flat_upper[idx]
      } else {
        paste(rhsmod$upper, collapse = ";")
      }
      flat_prior[idx] <- if (is.null(rhsmod$prior)) {
        flat_prior[idx]
      } else {
        paste(rhsmod$prior, collapse = ";")
      }
      flat_efa[idx] <- if (is.null(lhsmod$efa)) {
        flat_efa[idx]
      } else {
        paste(lhsmod$efa, collapse = ";")
      }
      flat_rv[idx] <- if (is.null(rhsmod$rv)) {
        flat_rv[idx]
      } else {
        paste(rhsmod$rv, collapse = ";")
      }
      modnu <- c(lhsmod, rhsmod)
      if (length(modnu) > 0L) { # there is a modifier here
        if (length(already) == 0) { # unknown element
          mod_idx <- mod_idx + 1L
          cur_mod_idx <- mod_idx
          mod[[cur_mod_idx]] <- modnu
          flat_rhs_mod_idx[idx] <- cur_mod_idx
        } else { # known element
          if (flat_rhs_mod_idx[idx] == 0) { # not yet modifier
            mod_idx <- mod_idx + 1L
            cur_mod_idx <- mod_idx
            mod[[cur_mod_idx]] <- modnu
            flat_rhs_mod_idx[idx] <- cur_mod_idx
          } else { # use existing modifier index
            cur_mod_idx <- flat_rhs_mod_idx[idx]
            overwrite <- names(modnu)[names(modnu) %in%
              names(mod[[cur_mod_idx]])]
            if (length(overwrite) > 0) {
              tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[rmeprev + 1L])
              lav_msg_warn(
                gettextf(
                  "modifier %s specified multiple times, overwritten",
                  overwrite[1L]
                ), tl[1L],
                footer = tl[2L]
              )
            }
            mod[[cur_mod_idx]] <- modifyList(mod[[cur_mod_idx]], modnu)
          }
        }
      }
    }
    # check for variable regressed on itself
    if (formul1$elem_text[opi] == "~" &&
        formul1$elem_text[opi - 1L] == formul1$elem_text[nelem]) {
      if (!grepl("^0\\.?0*$", flat_fixed[idx])) {
        tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[opi])
        lav_msg_stop(
          gettext("a variable cannot be regressed on itself"),
          tl[1L],
          footer = tl[2L]
        )
      }
    }
  }
  # create flat (omit items without operator)
  filled_ones <- which(flat_op != "")
  flat <- list(
    lhs = flat_lhs[filled_ones],
    op = flat_op[filled_ones],
    rhs = flat_rhs[filled_ones],
    mod.idx = flat_rhs_mod_idx[filled_ones],
    block = flat_block[filled_ones],
    fixed = flat_fixed[filled_ones],
    start = flat_start[filled_ones],
    lower = flat_lower[filled_ones],
    upper = flat_upper[filled_ones],
    label = flat_label[filled_ones],
    prior = flat_prior[filled_ones],
    efa = flat_efa[filled_ones],
    rv = flat_rv[filled_ones]
  )
  # change op for intercepts (for convenience only)
  int_idx <- which(flat_op == "~" & flat_rhs == "")
  if (length(int_idx) > 0L) {
    flat$op[int_idx] <- "~1"
  }
  # change op for ininstruments (for convenience only)
  int_idx <- which(flat_op == "|~" & flat_rhs == "")
  if (length(int_idx) > 0L) {
    flat$op[int_idx] <- "|~1"
  }
  # if there are constraints that are simple lower or upper limits, put
  # them in these members, add a modifier and remove constraint
  aantal <- length(constraints)
  if (aantal > 0) {
    for (j in aantal:1) {
      if (any(flat$label == constraints[[j]]$lhs) &&
          any(constraints[[j]]$op == c("<", ">"))) {
          rhslang <- str2lang(constraints[[j]]$rhs)
          numbound <- NA_real_
          if (mode(rhslang) == "numeric") {
            numbound <- as.numeric(constraints[[j]]$rhs)
          } else {
            if (mode(rhslang) == "call") {
              if (is.numeric(tryCatch(eval(rhslang),
                                      error = function(e) "error"))) {
                numbound <- eval(rhslang)
              }
            }
          }
          if (!is.na(numbound)) {
            nrs <- which(flat$label == constraints[[j]]$lhs)
            for (nr in nrs) {
              nrm <- length(mod) + 1L
              if (flat$mod.idx[nr] > 0L) {
                nrm <- flat$mod.idx[nr]
              } else {
                flat$mod.idx[nr] <- nrm
                mod <- c(mod, list(label = constraints[[j]]$lhs))
              }
              if (constraints[[j]]$op == "<") {
                flat$upper[nr] <- as.character(numbound)
                mod[[nrm]]$upper <- numbound
              } else {
                flat$lower[nr] <- as.character(numbound)
                mod[[nrm]]$lower <- numbound
              }
            }
            constraints <- constraints[-j]
          }
        }
    }
  }
  # new in 0.6, reorder covariances here!
  flat <- lav_partable_covariance_reorder(flat)
  if (as_data_frame) {
    flat <- as.data.frame(flat, stringsAsFactors = FALSE)
  }
  # new in 0.6-4: check for 'group' within 'level'
  if (any(flat_op == ":")) {
    op_idx <- which(flat_op == ":")
    if (length(op_idx) < 2L) {
      # only 1 block identifier? this is weird -> give warning
      lav_msg_warn(gettext("syntax contains only a single block identifier!"))
    } else {
      first_block <- flat_lhs[op_idx[1L]]
      second_block <- flat_lhs[op_idx[2L]]
      if (first_block == "level" && second_block == "group") {
        lav_msg_stop(gettext("groups can not be nested within levels!"))
      }
    }
  }
  # create output
  attr(flat, "modifiers") <- mod
  attr(flat, "constraints") <- constraints
  assign(hashstring, flat, envir = lavaan_cache_env)
  flat
}
