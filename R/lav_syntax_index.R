# ------------------------ indexed names ------------------------------------- #
# YR 17 Oct 2026: indexed names for models of repeated measures (wide data)
#
#   wx[t in 2:4] ~ a*wx[t-1] + b*wy[t-1]
#
# A name followed by an index in square brackets either BINDS an index to a
# range of integers (x[t in 1:4]) or REFERS to an index bound elsewhere in the
# same statement (x[t], x[t-1], x[t+2]). A statement with bindings stands for
# one copy of itself per value (per combination of values) of its indices;
# the name is written with the value after it (x[t-1] -> x2 for t = 3).
#
# This is syntactic sugar: the copies are ordinary statements, and the
# expansion happens on the tokens, before the '++' expansion. Every generated
# name token remembers its index and offset (elem_index, elem_offset), so that
# the parser can store, for each element of the parameter table, the lag
# between the left- and the right-hand side (the 'lag' column).
#
# Rules (see ?model.syntax, section 'Indexed names'):
# - a range lo:hi is ascending; a bound is an integer, an index bound earlier
#   in the statement, or (index + k) / (index - k) in parentheses
# - two integer bounds with lo > hi: error; a dependent range with lo > hi is
#   empty (that copy disappears)
# - a leading zero in an integer bound sets the width (y[t in 03:10] -> y03)
# - a bound index used as a number (in a modifier or an expression) stands
#   for its value: s =~ (t-1)*y[t in 1:5]; a bare index as the whole modifier
#   (t*x[t]) is an error (it reads as a label)
# - no space between a name and its index; no text after the index
# ---------------------------------------------------------------------------- #

lav_expand_index <- function(modellist, modelsrc, types) {
  is_open <- modellist$elem_type == types$symbol &
    modellist$elem_text %in% c("[", "]")
  if (!any(is_open)) {
    return(modellist)
  }
  ntok <- length(modellist$elem_type)
  modellist$elem_index <- character(ntok)
  modellist$elem_offset <- integer(ntok)

  out <- list(
    elem_pos = integer(0L), elem_type = integer(0L),
    elem_text = character(0L), elem_formula_number = integer(0L),
    elem_index = character(0L), elem_offset = integer(0L)
  )
  add_tokens <- function(out, tok, frm) {
    out$elem_pos <- c(out$elem_pos, tok$elem_pos)
    out$elem_type <- c(out$elem_type, tok$elem_type)
    out$elem_text <- c(out$elem_text, tok$elem_text)
    out$elem_formula_number <- c(out$elem_formula_number,
                                 rep(frm, length(tok$elem_type)))
    out$elem_index <- c(out$elem_index, tok$elem_index)
    out$elem_offset <- c(out$elem_offset, tok$elem_offset)
    out
  }
  newline_token <- function(pos) {
    list(elem_pos = pos, elem_type = types$newline, elem_text = "\n",
         elem_index = "", elem_offset = 0L)
  }

  frm_numbers <- unique(modellist$elem_formula_number)
  new_frm <- 0L
  for (frm in frm_numbers) {
    idx <- which(modellist$elem_formula_number == frm)
    stmt <- lav_parse_sublist(modellist, idx)
    has_bracket <- any(stmt$elem_type == types$symbol &
                         stmt$elem_text %in% c("[", "]"))
    if (!has_bracket) {
      new_frm <- new_frm + 1L
      out <- add_tokens(out, stmt, new_frm)
      next
    }
    # the statement without its line breaks
    keep <- stmt$elem_type != types$newline
    trailing_pos <- stmt$elem_pos[length(stmt$elem_pos)]
    stmt <- lav_parse_sublist(stmt, which(keep))
    copies <- lav_expand_index_statement(stmt, modelsrc, types)
    if (length(copies) == 0L) {
      tl <- lav_parse_txtloc(modelsrc, stmt$elem_pos[1L])
      lav_msg_warn(
        gettext("this statement with indexed names generates no elements
                (all index ranges are empty)"),
        tl[1L], footer = tl[2L]
      )
    }
    for (cp in copies) {
      new_frm <- new_frm + 1L
      out <- add_tokens(out, cp, new_frm)
      out <- add_tokens(out, newline_token(trailing_pos), new_frm)
    }
  }
  out
}

# expand one statement (a token list without line breaks) into a list of
# token lists, one per combination of index values
lav_expand_index_statement <- function(stmt, modelsrc, types) {
  ntok <- length(stmt$elem_type)
  txt <- stmt$elem_text
  tp <- stmt$elem_type
  pos <- stmt$elem_pos
  is_sym <- function(i, s) tp[i] == types$symbol && txt[i] == s
  where <- function(i) lav_parse_txtloc(modelsrc, pos[i])
  is_int <- function(s) grepl("^[0-9]+$", s)

  # 1. locate the indexed names: name '[' ... ']'
  groups <- list()
  in_bracket <- logical(ntok)
  i <- 1L
  while (i <= ntok) {
    if (is_sym(i, "]")) {
      tl <- where(i)
      lav_msg_stop(gettext("unexpected ']' (an index starts with '[')"),
                   tl[1L], footer = tl[2L])
    }
    if (!is_sym(i, "[")) {
      i <- i + 1L
      next
    }
    if (i == 1L || tp[i - 1L] != types$identifier) {
      tl <- where(i)
      lav_msg_stop(gettext("only a name can have an index"),
                   tl[1L], footer = tl[2L])
    }
    if (pos[i - 1L] + nchar(txt[i - 1L]) != pos[i]) {
      tl <- where(i)
      lav_msg_stop(gettext("there can be no space between a name and its
                         index"),
                   tl[1L], footer = tl[2L])
    }
    j <- i + 1L
    while (j <= ntok && !is_sym(j, "]")) {
      if (tp[j] == types$lavaanoperator && txt[j] != ":") {
        tl <- where(i)
        lav_msg_stop(gettext("index without closing ']'"),
                     tl[1L], footer = tl[2L])
      }
      if (is_sym(j, "[")) {
        tl <- where(j)
        lav_msg_stop(gettext("an index can not contain another index"),
                     tl[1L], footer = tl[2L])
      }
      j <- j + 1L
    }
    if (j > ntok) {
      tl <- where(i)
      lav_msg_stop(gettext("index without closing ']'"),
                   tl[1L], footer = tl[2L])
    }
    if (j == i + 1L) {
      tl <- where(i)
      lav_msg_stop(gettext("empty index"),
                   tl[1L], footer = tl[2L])
    }
    if (j < ntok && tp[j + 1L] %in% c(types$identifier, types$numliteral) &&
        pos[j + 1L] == pos[j] + 1L) {
      tl <- where(j + 1L)
      lav_msg_stop(gettext("text after an index is not supported; the
                              index must come at the end of the name"),
                   tl[1L], footer = tl[2L])
    }
    groups[[length(groups) + 1L]] <- lav_expand_index_parse(
      i - 1L, i, j, stmt, types, where, is_int
    )
    in_bracket[seq.int(i, j)] <- TRUE
    i <- j + 1L
  }

  # 2. the bindings, in order of appearance
  bind_flag <- vapply(groups, function(g) g$binding, logical(1L))
  bindings <- groups[bind_flag]
  if (length(bindings) == 0L) {
    tl <- where(groups[[1L]]$name_i)
    lav_msg_stop(gettextf("the index '%s' is not bound in this statement (bind it
                     once, as in %s[%s in 1:4])", groups[[1L]]$index,
                     txt[groups[[1L]]$name_i], groups[[1L]]$index),
                 tl[1L], footer = tl[2L])
  }
  bnames <- vapply(bindings, function(g) g$index, character(1L))
  if (anyDuplicated(bnames)) {
    dup <- bindings[[anyDuplicated(bnames)]]
    tl <- where(dup$name_i)
    lav_msg_stop(gettextf("the index '%s' is bound twice in this
                                 statement", dup$index),
                 tl[1L], footer = tl[2L])
  }
  for (b in seq_along(bindings)) {
    for (bnd in list(bindings[[b]]$lo, bindings[[b]]$hi)) {
      if (!is.null(bnd$ref) && !bnd$ref %in% bnames[seq_len(b - 1L)]) {
        tl <- where(bindings[[b]]$name_i)
        lav_msg_stop(gettextf("the bound '%s' of a range must be an index bound
                         earlier in this statement", bnd$ref),
                     tl[1L], footer = tl[2L])
      }
    }
    if (is.null(bindings[[b]]$lo$ref) && is.null(bindings[[b]]$hi$ref) &&
        bindings[[b]]$lo$value > bindings[[b]]$hi$value) {
      tl <- where(bindings[[b]]$name_i)
      lav_msg_stop(gettext("an index range must be ascending (lo:hi with
                      lo <= hi)"),
                   tl[1L], footer = tl[2L])
    }
  }
  for (g in groups[!bind_flag]) {
    if (!g$index %in% bnames) {
      tl <- where(g$name_i)
      lav_msg_stop(gettextf("the index '%s' is not bound in this
                                 statement", g$index),
                   tl[1L], footer = tl[2L])
    }
  }
  width <- vapply(bindings, function(g) g$width, integer(1L))
  names(width) <- bnames

  # 3. what kind of statement? (the first operator outside the brackets,
  #    other than ':' (block or interaction))
  ops <- which(tp == types$lavaanoperator & !in_bracket & txt != ":")
  op <- if (length(ops) > 0L) txt[ops[1L]] else ""
  expression_stmt <- op %in% c(":=", ":~", "==", "<", ">")

  # 4. bare indices (an index used as a number)
  bare <- which(tp == types$identifier & !in_bracket & txt %in% bnames)
  name_tokens <- vapply(groups, function(g) g$name_i, integer(1L))
  bare <- setdiff(bare, name_tokens)
  bare <- bare[!(bare < ntok & vapply(bare + 1L, function(k)
    k <= ntok && is_sym(k, "("), logical(1L)))] # a function name
  if (!expression_stmt) {
    term_end <- c("+", "++")
    for (k in bare) {
      nxt_is_end <- k == ntok ||
        (tp[k + 1L] == types$symbol && txt[k + 1L] %in% term_end) ||
        tp[k + 1L] == types$lavaanoperator
      if (nxt_is_end) {
        tl <- where(k)
        lav_msg_stop(gettextf("the index '%s' can not be used as a variable;
                            a name with an index is written as x[%s]",
                            txt[k], txt[k]),
                     tl[1L], footer = tl[2L])
      }
      prv_is_start <- k == 1L ||
        (tp[k - 1L] == types$symbol &&
           txt[k - 1L] %in% c("+", "++", "*", "?")) ||
        tp[k - 1L] == types$lavaanoperator
      nxt_is_mod <- tp[k + 1L] == types$symbol && txt[k + 1L] %in% c("*", "?")
      if (prv_is_start && nxt_is_mod) {
        tl <- where(k)
        lav_msg_stop(gettextf("'%1$s' is an index here, and can not be a label;
                            use fixed(%1$s) for its value as a fixed value",
                            txt[k]),
                     tl[1L], footer = tl[2L])
      }
    }
  }

  # 5. the copies
  copies <- list()
  lav_expand_index_emit <- function(env) {
    tok <- list(elem_pos = integer(0L), elem_type = integer(0L),
                elem_text = character(0L), elem_index = character(0L),
                elem_offset = integer(0L))
    put <- function(tok, p, ty, te, ix = "", of = 0L) {
      tok$elem_pos <- c(tok$elem_pos, p)
      tok$elem_type <- c(tok$elem_type, ty)
      tok$elem_text <- c(tok$elem_text, te)
      tok$elem_index <- c(tok$elem_index, ix)
      tok$elem_offset <- c(tok$elem_offset, of)
      tok
    }
    gi <- 1L
    k <- 1L
    while (k <= ntok) {
      if (gi <= length(groups) && k == groups[[gi]]$name_i) {
        g <- groups[[gi]]
        v <- env[[g$index]] + g$offset
        if (v < 0L) {
          tl <- where(k)
          lav_msg_stop(gettextf("the index value %1$s for '%2$s' is negative",
                              v, txt[k]),
                       tl[1L], footer = tl[2L])
        }
        vtxt <- if (width[[g$index]] > 0L) {
          formatC(v, width = width[[g$index]], flag = "0")
        } else {
          as.character(v)
        }
        tok <- put(tok, pos[k], types$identifier, paste0(txt[k], vtxt),
                   g$index, g$offset)
        k <- g$close_i + 1L
        gi <- gi + 1L
        next
      }
      if (k %in% bare) {
        tok <- put(tok, pos[k], types$numliteral, as.character(env[[txt[k]]]))
      } else {
        tok <- put(tok, pos[k], tp[k], txt[k])
      }
      k <- k + 1L
    }
    tok
  }
  bound_value <- function(bnd, env) {
    if (is.null(bnd$ref)) bnd$value else env[[bnd$ref]] + bnd$value
  }
  recurse <- function(b, env) {
    if (b > length(bindings)) {
      copies[[length(copies) + 1L]] <<- lav_expand_index_emit(env)
      return(invisible(NULL))
    }
    lo <- bound_value(bindings[[b]]$lo, env)
    hi <- bound_value(bindings[[b]]$hi, env)
    if (lo > hi) {
      return(invisible(NULL)) # an empty dependent range
    }
    for (v in seq.int(lo, hi)) {
      env[[bnames[b]]] <- v
      recurse(b + 1L, env)
    }
    invisible(NULL)
  }
  recurse(1L, list())
  copies
}

# parse the tokens of one index: name '[' ... ']'
# returns list(name_i, close_i, index, binding, offset, lo, hi, width)
lav_expand_index_parse <- function(name_i, open_i, close_i, stmt, types,
                                   where, is_int) {
  txt <- stmt$elem_text
  tp <- stmt$elem_type
  inner <- seq.int(open_i + 1L, close_i - 1L)
  n <- length(inner)
  t1 <- inner[1L]
  if (tp[t1] != types$identifier) {
    tl <- where(t1)
    lav_msg_stop(gettext("an index is a name (x[t]), possibly with an
                        integer offset (x[t-1]), or a binding (x[t in 1:4])"),
                 tl[1L], footer = tl[2L])
  }
  res <- list(name_i = name_i, close_i = close_i, index = txt[t1],
              binding = FALSE, offset = 0L, lo = NULL, hi = NULL,
              width = 0L)
  if (n == 1L) {
    return(res)
  }
  t2 <- inner[2L]
  # reference with an offset: t + k, t - k (the sign may be glued to k)
  if (tp[t2] == types$symbol && txt[t2] %in% c("+", "-")) {
    if (n != 3L || tp[inner[3L]] != types$numliteral ||
        !is_int(txt[inner[3L]])) {
      tl <- where(t2)
      lav_msg_stop(gettext("the offset of an index must be an integer
                          (x[t-1], x[t+2])"),
                   tl[1L], footer = tl[2L])
    }
    k <- as.integer(txt[inner[3L]])
    res$offset <- if (txt[t2] == "-") -k else k
    return(res)
  }
  if (tp[t2] == types$numliteral && grepl("^[-+][0-9]+$", txt[t2]) &&
      n == 2L) {
    res$offset <- as.integer(txt[t2])
    return(res)
  }
  # binding: t in lo : hi
  if (tp[t2] != types$identifier || txt[t2] != "in") {
    tl <- where(t2)
    lav_msg_stop(gettext("an index is a name (x[t]), possibly with an
                        integer offset (x[t-1]), or a binding (x[t in 1:4])"),
                 tl[1L], footer = tl[2L])
  }
  colon <- inner[tp[inner] == types$lavaanoperator & txt[inner] == ":"]
  if (length(colon) != 1L) {
    tl <- where(t2)
    lav_msg_stop(gettext("the range of an index is written lo:hi
                        (x[t in 1:4])"),
                 tl[1L], footer = tl[2L])
  }
  lav_expand_index_bound <- function(from, to) {
    if (from > to) {
      tl <- where(colon)
      lav_msg_stop(gettext("the range of an index is written lo:hi
                             (x[t in 1:4])"),
                   tl[1L], footer = tl[2L])
    }
    m <- to - from + 1L
    if (m == 1L && tp[from] == types$numliteral && is_int(txt[from])) {
      w <- if (nchar(txt[from]) > 1L && substr(txt[from], 1L, 1L) == "0") {
        nchar(txt[from])
      } else {
        0L
      }
      return(list(ref = NULL, value = as.integer(txt[from]), width = w))
    }
    if (m == 1L && tp[from] == types$identifier) {
      return(list(ref = txt[from], value = 0L, width = 0L))
    }
    # ( s + k ) or ( s - k )
    if (m == 5L && tp[from] == types$symbol && txt[from] == "(" &&
        tp[from + 1L] == types$identifier &&
        tp[from + 2L] == types$symbol && txt[from + 2L] %in% c("+", "-") &&
        tp[from + 3L] == types$numliteral && is_int(txt[from + 3L]) &&
        tp[from + 4L] == types$symbol && txt[from + 4L] == ")") {
      k <- as.integer(txt[from + 3L])
      return(list(ref = txt[from + 1L],
                  value = if (txt[from + 2L] == "-") -k else k, width = 0L))
    }
    tl <- where(from)
    lav_msg_stop(gettext("a bound of an index range is an integer, an index
                          bound earlier in the statement, or such an index
                          plus or minus an integer, in parentheses:
                          x[t in (s+1):4]"),
                 tl[1L], footer = tl[2L])
  }
  res$binding <- TRUE
  res$lo <- lav_expand_index_bound(inner[3L], colon - 1L)
  res$hi <- lav_expand_index_bound(colon + 1L, close_i - 1L)
  res$width <- max(res$lo$width, res$hi$width)
  res
}
