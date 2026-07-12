# ---------------------lav_parse_options ------------------------------------- #
# function to configurate settings which conduct the parser function           #
# can be used to add operators, modifiers and grouping classes                 #
# the function returns the status of the parser config                         #
# ---------------------------------------------------------------------------- #
lav_parse_options <- function(
    pkgname   = NULL, # must be given if any of the following arguments not NULL
    operators = NULL, # char vector
    modifiers = NULL, # data.frame(mod = ., side = ., expect = .)
    groupings  = NULL # char vector
) {
  stopifnot(is.null(pkgname) || (is.character(pkgname) && pkgname != "lavaan"))
  stopifnot(is.null(operators) ||
    (is.character(operators) && !is.null(pkgname)))
  stopifnot(is.null(modifiers) ||
    (is.data.frame(modifiers) && !is.null(pkgname)))
  stopifnot(is.null(groupings) ||
    (is.character(groupings) && !is.null(pkgname)))
  config <- get0("mdl_config", lavaan_cache_env)
  if (is.null(config)) {
    ops = data.frame(
      op = c("=~", "<~", "~*~", "~~", "~", "|~", "==",
      "<", ">", ":=", ":", "|", "%", ":~"),
      pkg = "lavaan"
    )
    mods = data.frame(
      mod = c("efa", "fixed", "label", "start",
      "lower", "upper", "rv", "equal"),
      side = c("l", "r", "r", "r",
        "r", "r", "r", "r"),
      expect = c("char", "num", "char", "num",
        "num", "num", "char", "char"),
      pkg = "lavaan"
    )
    groups = data.frame(
      group = c("group", "level", "class", "block"),
      pkg = "lavaan"
    )
    config <- list(operators = ops,
      modifiers = mods, groupings = groups)
    assign("mdl_config", config, envir = lavaan_cache_env)
  }
  if (is.null(pkgname)) return(config)
  if (!exists("ops", environment())) {
    ops    <- config$operators
    mods   <- config$modifiers
    groups <- config$groupings
  }
  # remove existing values
  ops <- ops[ops$pkg != pkgname, ]
  mods <- mods[mods$pkg != pkgname, ]
  groups <- groups[groups$pkg != pkgname, ]
  # add supplied values
  if (!is.null(operators)) {
    if (any(!grepl("^%.+%$", operators))) {
      lav_msg_stop(gettext("Operators must have the format '%.+%'."))
    }
    if (any(operators %in% ops$op)) {
      lav_msg_stop(gettext("Some operators already defined."))
    }
    ops <- rbind(ops,
    data.frame(
      op = operators,
      pkg = pkgname
    ))
  }
  if (!is.null(modifiers)) {
    if (is.null(modifiers$mod)) {
      lav_msg_stop(gettext("Modifiers must have a column 'mod'."))
    }
    if (any(!grepl("^[a-z_]+$", modifiers$mod))) {
      lav_msg_stop(gettext("Modifiers$mod must have the format '[a-z_]+'."))
    }
    if (is.null(modifiers$side)) {
      lav_msg_stop(gettext("Modifiers must have a column 'side'."))
    }
    if (any(!grepl("^l$|^r$", modifiers$side))) {
      lav_msg_stop(gettext("Modifiers$side must contain 'l' or 'r'."))
    }
    if (is.null(modifiers$expect)) {
      lav_msg_stop(gettext("Modifiers must have a column 'expect'."))
    }
    if (any(!grepl("^char$|^num$|^raw$", modifiers$expect))) {
      lav_msg_stop(gettext(
        "Modifiers$expect must contain 'char', 'num' or 'raw'."))
    }
    forbidden <- c("lhs", "op", "rhs", "block")
    if (any(modifiers %in% forbidden)) {
      lav_msg_stop(gettextf("Modifiers cannot be %s.",
                    lav_msg_view(forbidden, "or")))
    }
    if (any(modifiers %in% mods$mod)) {
      lav_msg_stop(gettext("Some modifiers already defined."))
    }
    mods <- rbind(mods,
    data.frame(
      mod = modifiers$mod,
      side = modifiers$side,
      expect = modifiers$expect,
      pkg = pkgname
    ))
  }
  if (!is.null(groupings)) {
    if (any(!grepl("^[a-z_]+$", groupings))) {
      lav_msg_stop(gettext("groupings must have the format '[a-z_]+'."))
    }
    if (any(groupings %in% groups$group)) {
      lav_msg_stop(gettext("Some groupings already defined."))
    }
    groups <- rbind(groups,
    data.frame(
      group = groupings,
      pkg = pkgname
    ))
  }
  config <- list(operators = ops,
    modifiers = mods, groupings = groups)
  assign("mdl_config", config, envir = lavaan_cache_env)
  config
}

# ------------------------ lav_parse_tokens_open ----------------------------- #
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
lav_parse_tokens_open <- function(modelsrc, types) {
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
      str_in_comment <- (elem_pos > comments[i] &
                           elem_pos < comments[i] + comment_lengths[i])
      if (any(str_in_comment)) {
        elem_type[str_in_comment] <- 0
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
  config <- lav_parse_options()
  for (opr in config$operators$op[nchar(config$operators$op) == 2L]) {
    oprre <- paste0("([^=~%])[", substr(opr, 1L, 1L), "]( +)[",
                                            substr(opr, 2L, 2L), "]([^=~%])")
    oprsub <- paste0("\\1", opr, "\\2\\3")
    if (grepl(oprre, modelsrcw)) {
      waar <- regexpr(oprre, modelsrcw)[1] + 1L
      modelsrcw <- gsub(oprre, oprsub, modelsrcw)
      tl <- lav_parse_txtloc(modelsrc, waar)
      lav_msg_warn(gettextf("splitting of operator '%s' deprecated", opr),
                  tl[1L],
                  footer = tl[2L]
      )
    }
  }
  # -------------------------------------------------------------------------- #
  lavopsre <- paste(
    gsub("\t", "\\\\", gsub("([.|()[{^$*+?])", "\t\\1",
             config$operators$op)), collapse = "|")
  lavops <- gregexpr(lavopsre, modelsrcw)[[1]]
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
  if (identifiers[1L] > -1L) {
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

  list_before_numbering <-  list(
    elem_pos = elem_pos, elem_type = elem_type,
    elem_text = elem_text
  )

  lav_parse_formula_numbers(list_before_numbering, types)
}


# ------------------------ lav_parse_formula_numbers ------------------------- #
# help function to add formula numbers to tokens
# ---------------------------------------------------------------------------- #
lav_parse_formula_numbers <- function(list_before_numbering, types) {
  elem_pos <- list_before_numbering$elem_pos
  elem_type <- list_before_numbering$elem_type
  elem_text <- list_before_numbering$elem_text
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
      c("+", "*", "=~", "-", "<~", "~*~", "~~", "~", "|~", "|", "%", ":~"))) {
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

# ------------------------ lav_parse_formulas_open ------------------------- #
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
lav_parse_formulas_open <- function(modellist, modelsrc, types) {
  config <- lav_parse_options()
  group_operator <- ":"
  define_operator <- ":="
  conditional_operators <- c("==", "<", ">")
  non_relational <- config$operators$op %in%
                      c(group_operator, define_operator, conditional_operators)
  relational_operators <- config$operators$op[!non_relational]
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
    if (any(formul1$elem_text[opi] == relational_operators) &&
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

# ------------------------ lav_parse_modifier_open --------------------------- #
# The function takes a list with tokens belonging to a single 'mono' lavaan
# formula as input. The other arguments are:
#       lhs : check for lhs or rhs modifier
#       opi : index of the lavaan operator in the list-items
#  modelsrc : the model source string (for error messages and warnings)
#     types : the types of tokens
#       rme : index of last element of modifier in formula (*)
#   rmeprev : index of first element of modifier in formula - 1L (*)
# The function returns the modifier detected as element of a list
#  with name the modifier type and value an array of values
#   (length > 1 if vector via c(...)) for the modifier value.
# The lhs or rhs is limited to the elements with index
#   rmeprev+1:rme, this is to support multiple modifiers for the same element.
# An error message is produced when no modifier can be determined.
# ---------------------------------------------------------------------------- #
lav_parse_modifier_open <- function(formul, opi, modelsrc, types,
                                    rme, rmeprev, modenv, config) {
  if (rme < opi) { # lhs modifier
    lhs <- TRUE
    welke <- c(seq.int(rmeprev + 1L, rme),
               seq.int(opi - 1L, length(formul$elem_type)))
  } else {
    lhs <- FALSE
    welke <- c(seq.int(1L, opi), seq.int(rmeprev + 1L, rme),
               length(formul$elem_type))
  }
  formul1 <- lav_parse_sublist(formul, welke)
  nelem <- length(formul1$elem_type)
  opi <- which(formul1$elem_type == types$lavaanoperator)[1L]
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
  getmodifier <- function(s) {
    if (s %in% c("c", "t")) {
      tmp <- s
      attr(tmp, "tiepe") <- "label"
      return(tmp)
    }
    v <- all.vars(str2expression(s))
    for (v1 in v) {
      if (v1 != "pi") assign(v1, v1, modenv)
    }
    tmp <- eval(str2expression(s), modenv)
    suppressWarnings(rm(list = v, pos = modenv))
    if (is.null(attr(tmp, "tiepe"))) {
      if (is.numeric(tmp) || is.logical(tmp)) {
        attr(tmp, "tiepe") <- "fixed"
      } else {
        attr(tmp, "tiepe") <- "label"
      }
    }
    tmp
  }
  strings <- vapply(seq.int(
    if (lhs) 1L else opi + 1L,
    if (lhs) opi - 3L else nelem - 2L
    ), function(x) {
    if (formul1$elem_type[x] == types$stringliteral) {
      paste0('"', formul1$elem_text[x], '"')
    } else {
      formul1$elem_text[x]
    }
  }, "")
  txt <- paste(strings, collapse = "")
  if (!lhs && formul1$elem_text[nelem - 1L] == "?") {
    formul1$elem_text[nelem - 1L] <- "*"
    txt <- paste0("start(", txt, ")")
  }
   modifier <- tryCatch(getmodifier(txt), error = function(e) "__error__")
  if (!is.na(modifier[1L]) && modifier[1L] == "__error__") {
    tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[opi + 1L])
    lav_msg_stop(gettext("invalid modifier specification"),
                  tl[1L],
                  footer = tl[2L]
    )
  }
  mdtype <- attr(modifier, "tiepe")
  mdind <- which(config$modifiers$mod == mdtype)[1L]
  if (config$modifiers$expect[mdind] == "char") {
    if (!is.character(modifier)) {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[opi + 1L])
      lav_msg_stop(gettextf("invalid %s modifier (should be character)",
        mdtype), tl[1L], footer = tl[2L])
    }
  } else if (config$modifiers$expect[mdind] == "num") {
    if (!is.numeric(modifier) && !all(is.na(modifier))) {
      tl <- lav_parse_txtloc(modelsrc, formul1$elem_pos[opi + 1L])
      lav_msg_stop(gettextf("invalid %s modifier (should be numeric)",
        mdtype), tl[1L], footer = tl[2L])
    }
  }
  modifierlist <- list(as.vector(modifier))
  names(modifierlist) <- attr(modifier, "tiepe")
  modifierlist
}

# -------------------- create environment for evaluating modifiers ----------- #
lav_parse_modenv <- function(config, side) {
  modenv <- new.env()
  for (j in which(config$modifiers$side == side)) {
    mod1 <- config$modifiers$mod[j]
    src <- paste0("assign('", mod1, "', function(...) {\n",
    if (config$modifiers$expect[j] == "raw") {
      "x <- sub('^.*?\\\\(', '', sub('\\\\)$', '', deparse(sys.call())))\n"
    } else {
      "x <- do.call('c', args = list(...))\n"
    },
    "attr(x, 'tiepe') <- '", mod1, "'\n",
    "x\n}, modenv)")
    source(exprs = str2expression(src), local = TRUE)
  }
  # equal is synonym for label
  if (side == "r") assign("equal", modenv$label, modenv)
  modenv
}

# -------------------- main parsing function --------------------------------- #
lav_parse_model_string_open <- function(model_syntax = "",
                                        as_data_frame = FALSE) {
  # replace 'strange' tildes (in some locales) (new in 0.6-6)
  modelsrc <- gsub(
    pattern = "\u02dc",
    replacement = "~",
    paste(unlist(model_syntax), "", collapse = "\n")
  )
  # check for empty, whitespace-only or comment-only model syntax
  if (length(model_syntax) == 0L ||
      !grepl("[^[:space:];]", gsub("[#!][^\n]*", "", modelsrc))) {
    lav_msg_stop(gettext("Model syntax is empty."))
  }
  hashstring <- paste0("mdl_o_", lav_char2hash(paste0(modelsrc, as_data_frame)))
  if (exists(hashstring, envir = lavaan_cache_env)) {
    return(get(hashstring, envir = lavaan_cache_env))
  }
  config <- lav_parse_options()
  modenvl <- lav_parse_modenv(config, "l")
  modenvr <- lav_parse_modenv(config, "r")
  types <- lav_create_enum(c(
    "identifier", "numliteral", "stringliteral",
    "symbol", "lavaanoperator", "newline"
  ))
  modellist <- lav_parse_tokens_open(modelsrc, types)
  if (lav_debug()) {
    print(data.frame(
      pos = modellist$elem_pos,
      type = types$enum_names[modellist$elem_type],
      text = modellist$elem_text,
      formula = modellist$elem_formula_number
    ))
  }
  formulalist <- lav_parse_formulas_open(modellist, modelsrc, types)
  #---- analyse syntax formulas and put in flat_-----
  max_mono_formulas <- length(formulalist)
  flat <- list()
  flat$lhs <- character(max_mono_formulas)
  flat$op <- character(max_mono_formulas)
  flat$rhs <- character(max_mono_formulas)
  flat$mod_idx <- integer(max_mono_formulas)
  flat$block <- integer(max_mono_formulas) # keep track of groups using ":" opr
  for (mod1 in config$modifiers$mod) {
    flat[[mod1]] <- character(max_mono_formulas) # only for display purposes!
  }
  tmplist <- list(
    flat = flat,
    flat_idx = 0L,
    mod_idx = 0L,
    constraints = list(),
    modlist = list(),
    block = 1L,
    block_op = FALSE
  )
  if (lav_debug()) {
    cat("formula to analyse:\n")
  }
  for (s in seq_along(formulalist)) {
    formule <- formulalist[[s]]
    tmplist <- lav_parse_handle_formule(formule, tmplist, types, modelsrc,
                                        modenvl, modenvr, config)
  }
  lav_parse_final_operations(tmplist$flat, tmplist$modlist, config,
    tmplist$constraints, as_data_frame, hashstring)
}

lav_parse_handle_formule <- function(formule, tmplist, types, modelsrc,
                                     modenvl, modenvr, config) {
  flat <- tmplist$flat
  flat_idx <- tmplist$flat_idx
  mod_idx <- tmplist$mod_idx
  constraints <- tmplist$constraints
  modlist <- tmplist$modlist
  block <- tmplist$block
  block_op <- tmplist$block_op
  constraint_operators <- c("==", "<", ">", ":=")
  if (lav_debug()) {
    cat(vapply(seq_along(formule$elem_type), function(j) {
      if (formule$elem_type[j] == types$stringliteral) {
        return(dQuote(formule$elem_text[j], FALSE))
      }
      formule$elem_text[j]
    }, ""), "\n")
  }
  nelem <- length(formule$elem_type)
  # where is the operator
  opi <- which(formule$elem_type == types$lavaanoperator)
  if (length(opi) > 1L) { # if more then 1 operator skip operators ':'
    opii <- 1L
    while (formule$elem_text[opi[opii]] == ":" && opii < length(opi)) {
      opii <- opii + 1L
    }
    opi <- opi[opii]
  }
  op <- formule$elem_text[opi]
  if (any(op == constraint_operators)) {         # ----- constraints ------- #
    # a constraint/definition needs a non-empty left- and right-hand side
    if (opi <= 1L || opi >= nelem) {
      tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[opi])
      lav_msg_stop(
        gettextf("Missing or empty %1$s of operator '%2$s'.",
                 if (opi <= 1L) gettext("left-hand side")
                 else gettext("right-hand side"),
                 op),
        tl[1L], footer = tl[2L]
      )
    }
    lhs <- paste(formule$elem_text[seq.int(1L, opi - 1L)], collapse = "")
    rhs <- paste(formule$elem_text[seq.int(opi + 1L, nelem)], collapse = "")
    constraints <- c(
      constraints,
      list(list(
        op = op,
        lhs = lhs,
        rhs = rhs,
        user = 1L
      ))
    )
    tmplist$constraints <- constraints
    return(tmplist)
  }
  if (op == ":") {           # --------------- block start ----------------- #
    correct_syntax <- gettextf("The correct syntax is: \"LHS : RHS\", where
    LHS is  block identifier (%s) and RHS is the block number or label.",
                             lav_msg_view(config$groupings$group, "or"))
    if (opi == 1L) {
      tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[1])
      lav_msg_stop(
        gettextf(
          "Missing block identifier. %s", correct_syntax
        ), tl[1L], footer = tl[2L]
      )
    }
    if (opi > 2L || all(tolower(formule$elem_text[1]) !=
                        config$groupings$group)) {
      tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[1])
      lav_msg_stop(
        gettextf(
          "Invalid block identifier. %s", correct_syntax
        ),
        tl[1L], footer = tl[2L]
      )
    }
    if (nelem != 3 || all(formule$elem_type[3] !=
            c(types$stringliteral, types$identifier, types$numliteral))) {
      tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[1])
      lav_msg_stop(
        gettextf("Syntax contains block identifier with missing or
              invalid number/label. %s", correct_syntax),
        tl[1L],
        footer = tl[2L]
      )
    }
    flat_idx <- flat_idx + 1L
    flat$lhs[flat_idx] <- formule$elem_text[1]
    flat$op[flat_idx] <- op
    flat$rhs[flat_idx] <- formule$elem_text[3]
    flat$mod_idx[flat_idx] <- 0L
    if (block_op) {
      block <- block + 1L
    } else {
      if (flat_idx != 1) {
        tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[1])
        lav_msg_warn(
          gettext("First block defined after other formula's"),
          tl[1L],
          footer = tl[2L]
        )
      }
    }
    flat$block[flat_idx] <- block
    tmplist$flat <- flat
    tmplist$flat_idx <- flat_idx
    tmplist$block <- block
    tmplist$block_op <- TRUE
    return(tmplist)
  }
  # ------------------ relational operators -------------------------------- #
  lav_parse_check_relational(formule, modelsrc, types, opi, nelem, op)
  tmpenv <- list2env(list(
    formule = formule,
    opi = opi,
    nelem = nelem
  ))
  lav_parse_update_relational(tmpenv, modelsrc, types, op)
  formule <- get("formule", tmpenv)
  opi <- get("opi", tmpenv)
  nelem <- get("nelem", tmpenv)
  rm(tmpenv)
  lhs <- formule$elem_text[opi - 1L]
  rhs <- formule$elem_text[nelem]
  already <- which(flat$lhs == lhs & flat$op == op & flat$block == block &
            (flat$rhs == rhs | (flat$rhs == "" &
              (op == "~" | op == "|~" | op == ":~") &
              formule$elem_type[nelem] == types$numliteral)))
  if (length(already) == 1L) {
    idx <- already
  } else {
    flat_idx <- flat_idx + 1L
    idx <- flat_idx
    flat$lhs[idx] <- lhs
    flat$op[idx] <- op
    flat$rhs[idx] <- rhs
    flat$block[idx] <- block
    if (formule$elem_type[nelem] == types$numliteral) {
      if (op == "~" || op == "|~"| op == ":~") flat$rhs[idx] <- ""
    }
  }
  # modifiers always come before an asterix or a question mark
  modelems <- which(formule$elem_type == types$symbol &
              (formule$elem_text == "*" | formule$elem_text == "?"))
  for (j in seq_along(modelems)) {
    if (sum(formule$elem_text[seq.int(opi, modelems[j])] == "(") !=
        sum(formule$elem_text[seq.int(opi, modelems[j])] == ")"))
      modelems[j] <- 0L
  }
  modelems <- modelems[modelems != 0L]
  if (length(modelems) != 0L) {
    modleftfirst <- TRUE
    modrightfirst <- TRUE
    for (rmei in seq_along(modelems)) {
      rme <- modelems[rmei]
      if (rme < opi) {
        rmeprev <- if (modleftfirst) 0L else modelems[rmei - 1L]
        modleftfirst <- FALSE
      } else {
        rmeprev <- if (modrightfirst) opi else modelems[rmei - 1L]
        modrightfirst <- FALSE
      }
      modnu <- lav_parse_modifier_open(
        formule,
        opi, modelsrc, types, rme, rmeprev,
        modenv = if (rme < opi) modenvl else modenvr,
        config
      )
      for (j in seq_along(config$modifiers$mod)) {
        mod1 <- config$modifiers$mod[j]
        flat[[mod1]][idx] <-
          if (is.null(modnu[[mod1]])) {
            flat[[mod1]][idx]
          } else {
            paste(modnu[[mod1]], collapse = ";")
          }
      }
      if (flat$mod_idx[idx] == 0) { # not yet modifier
        mod_idx <- mod_idx + 1L
        cur_mod_idx <- mod_idx
        modlist[[cur_mod_idx]] <- modnu
        flat$mod_idx[idx] <- cur_mod_idx
      } else { # use existing modifier index
        cur_mod_idx <- flat$mod_idx[idx]
        overwrite <- names(modnu)[names(modnu) %in%
                                    names(modlist[[cur_mod_idx]])]
        if (length(overwrite) > 0) {
          tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[rmeprev + 1L])
          lav_msg_warn(
            gettextf(
              "modifier %s specified multiple times, overwritten",
              overwrite[1L]
            ), tl[1L],
            footer = tl[2L]
          )
        }
        modlist[[cur_mod_idx]] <- modifyList(modlist[[cur_mod_idx]], modnu)
      }
    }
  }
  # check for variable regressed on itself
  if (formule$elem_text[opi] %in% c("~", ":~") &&
      formule$elem_text[opi - 1L] == formule$elem_text[nelem]) {
    if (!grepl("^0\\.?0*$", flat$fixed[idx])) {
      tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[opi])
      lav_msg_stop(
        gettext("a variable cannot be regressed on itself"),
        tl[1L],
        footer = tl[2L]
      )
    }
  }
  list(
    flat = flat,
    flat_idx = flat_idx,
    mod_idx = mod_idx,
    constraints = constraints,
    modlist = modlist,
    block = block,
    block_op = block_op
  )
}

lav_parse_check_relational <- function(formule, modelsrc, types, opi,
                                       nelem, op) {
  # warn if some identifiers contain spaces
  contsp <- which(formule$elem_type == types$identifier &
                    grepl(" ", formule$elem_text, fixed = TRUE))
  if (length(contsp) > 0L) {
    tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[contsp[1L]])
    lav_msg_warn(
      gettextf(
        "having identifiers with spaces ('%s') is deprecated",
        formule$elem_text[contsp[1]]
      ),
      tl[1L],
      footer = tl[2L]
    )
  }
  # checks for valid names in lhs and rhs
  lav_parse_check_name(formule, opi - 1L, modelsrc) # valid name lhs
  for (j in seq.int(opi + 1L, nelem)) { # valid names rhs
    if (formule$elem_type[j] == types$identifier &&
        formule$elem_text[j] != "NA") {
      lav_parse_check_name(formule, j, modelsrc)
    }
  }
  if (formule$elem_type[nelem] != types$identifier &&
      (formule$elem_type[nelem] != types$numliteral ||
        all(op != c("~", "|~", "=~", ":~")))) {
    tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[nelem])
    lav_msg_stop(
      gettext("Last element of rhs part expected to be an identifier or,
              for operator ~, |~, =~ or :~, a numeric literal!"),
      tl[1L],
      footer = tl[2L]
    )
  }
   if (nelem < opi + 1L &&
     !(formule$elem_text[nelem - 1L] %in% c("*", "?"))) {
    tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[nelem - 1L])
    lav_msg_stop(gettext("invalid modifier symbol (should be '*' or '?')"),
                  tl[1L],
                  footer = tl[2L]
    )
  }
  invisible(NULL)
}

lav_parse_update_relational <- function(tmpenv, modelsrc, types, op) {
  formule <- tmpenv$formule
  opi <- tmpenv$opi
  nelem <- tmpenv$nelem
  # intercept fixed on 0
  # replace 'lhs ~ 0' => 'lhs ~ 0 * 1' - intercept fixed on zero
  if (formule$elem_text[nelem] == "0" && op == "~" && opi == nelem - 1L) {
    formule$elem_type <- c(formule$elem_type, types$symbol, types$numliteral)
    formule$elem_text <- c(formule$elem_text, "*", "1")
    formule$elem_pos <- c(formule$elem_pos, rep(formule$elem_pos[nelem], 2))
    nelem <- length(formule$elem_type)
  }
  # phantom latent variable
  # replace 'lhs =~ 0' => 'lhs =~ fixed(0)*lhs', 0 can be other numliteral
  #          also, lhs is last element before '=~'
  if (formule$elem_type[nelem] == types$numliteral && op == "=~") {
    formule$elem_type <- c(
      formule$elem_type[seq.int(1L, nelem - 1L)], types$identifier,
      types$symbol, types$numliteral, types$symbol, types$symbol,
      types$identifier
    )
    formule$elem_text <- c(
      formule$elem_text[seq.int(1L, nelem - 1L)], "fixed", "(",
      formule$elem_text[nelem], ")", "*", formule$elem_text[opi - 1L]
    )
    formule$elem_pos <- c(
      formule$elem_pos[seq.int(1L, nelem - 1L)],
      rep(formule$elem_pos[nelem], 6)
    )
    nelem <- length(formule$elem_type)
  }
  # handling interaction variable types
  colons <- which(formule$elem_text[seq.int(1L, nelem - 1L)] == ":" &
                    formule$elem_type[seq.int(2L, nelem)] == types$identifier)
  # check at most 1 colon
  if (length(colons) > 2L ||
      (length(colons) == 2L && (colons[1L] > opi || colons[2L] < opi))) {
    tl <- lav_parse_txtloc(modelsrc, formule$elem_pos[colons[2]])
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
    formule$elem_text[colons[1L] - 1L] <-
      paste(formule$elem_text[seq.int(colons[1L] - 1L, colons[1L] + 1L)],
            collapse = ""
      )
    formule <- lav_parse_sublist(formule,
          setdiff(seq.int(1L, nelem), seq.int(colons[1L], colons[1L] + 1L)))
    nelem <- length(formule$elem_type)
    if (colons[1L] < opi) {
      opi <- opi - 2L # is in LHS
      if (length(colons) == 2L) colons[2L] <- colons[2L] - 2L
    }
  }
  if (length(colons) == 2L) {
    # collapse items around colon "a" ":" "b" => "a:b"
    formule$elem_text[colons[2L] - 1L] <-
      paste(formule$elem_text[seq.int(colons[2L] - 1L, colons[2L] + 1L)],
            collapse = ""
      )
    formule <- lav_parse_sublist(formule,
          setdiff(seq.int(1L, nelem), seq.int(colons[2L], colons[2L] + 1L)))
    nelem <- length(formule$elem_type)
  }
  assign("formule", formule, tmpenv)
  assign("opi", opi, tmpenv)
  assign("nelem", nelem, tmpenv)
  invisible(NULL)
}

lav_parse_final_operations <- function(flat, modlist, config,
                        constraints, as_data_frame, hashstring) {
  # update flat (omit items without operator)
  filled_ones <- which(flat$op != "")
  flat$lhs <- flat$lhs[filled_ones]
  flat$op <- flat$op[filled_ones]
  flat$rhs <- flat$rhs[filled_ones]
  flat$mod.idx <- flat$mod_idx[filled_ones]
  flat$mod_idx <- NULL
  flat$block <- flat$block[filled_ones]

  # change op for intercepts (for convenience only)
  int_idx <- which(flat$op == "~" & flat$rhs == "")
  if (length(int_idx) > 0L) {
    flat$op[int_idx] <- "~1"
  }
  # change op for instruments (for convenience only)
  int_idx <- which(flat$op == "|~" & flat$rhs == "")
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
            nrm <- length(modlist) + 1L
            if (flat$mod.idx[nr] > 0L) {
              nrm <- flat$mod.idx[nr]
            } else {
              flat$mod.idx[nr] <- nrm
              modlist <- c(modlist, list(label = constraints[[j]]$lhs))
            }
            if (constraints[[j]]$op == "<") {
              flat$upper[nr] <- as.character(numbound)
              modlist[[nrm]]$upper <- numbound
            } else {
              flat$lower[nr] <- as.character(numbound)
              modlist[[nrm]]$lower <- numbound
            }
          }
          constraints <- constraints[-j]
        }
      }
    }
  }
  # new in 0.6, reorder covariances here!
  flat <- lav_pt_covariance_reorder(flat)
  # remove empty modifier columns
  for (mod1 in config$modifiers$mod) {
    flat[[mod1]] <- flat[[mod1]][filled_ones]
    if (all(flat[[mod1]] == "")) flat[[mod1]] <- NULL
  }
  # transform to data.frame is asked
  if (as_data_frame) {
    flat <- as.data.frame(flat, stringsAsFactors = FALSE)
  }
  # new in 0.6-4: check for 'group' within 'level'
  if (any(flat$op == ":")) {
    op_idx <- which(flat$op == ":")
    if (length(op_idx) < 2L) {
      # only 1 block identifier? this is weird -> give warning
      lav_msg_warn(gettext("syntax contains only a single block identifier!"))
    } else {
      first_block <- flat$lhs[op_idx[1L]]
      second_block <- flat$lhs[op_idx[2L]]
      if (first_block == "level" && second_block == "group") {
        lav_msg_stop(gettext("groups can not be nested within levels!"))
      }
    }
  }
  # create output
  attr(flat, "modifiers") <- modlist
  attr(flat, "constraints") <- constraints
  assign(hashstring, flat, envir = lavaan_cache_env)
  flat
}
