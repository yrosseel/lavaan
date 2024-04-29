# Displays a message (... concatenated with spaces in between) with header
# 'lavaan(function):' and formatted to have a maximum line length of 'txt.width'
# while all but the first line start with 'indent' spaces. The message is shown
# via R function 'message()'.
lav_msg_note <- function(..., showheader = FALSE) {
  wat <- unlist(list(...), use.names = FALSE)
  if (!showheader) wat <- c("lavaan NOTE: ___", wat)
  message(lav_msg(wat, showheader = showheader), domain = NA)
}

# Displays a message with header and formatted as
# above via R function 'warning()'.
lav_msg_warn <- function(...) {
  wat <- unlist(list(...), use.names = FALSE)
  warning(lav_msg(wat), call. = FALSE, domain = NA)
}

# Displays a message with header and formatted as
# above via R function 'stop()'.
lav_msg_stop <- function(...) {
  wat <- unlist(list(...), use.names = FALSE)
  stop(lav_msg(wat), call. = FALSE, domain = NA)
}

# Displays a message with header and formatted as
# above via R function 'stop()', where the message is prepended with "FIXME:",
# to indicate an internal error, e.g. an error condition which was supposed
# to be handled in the calling functions. Such error message do not have to
# be created by [n]gettext[f] because they don't have to be translated!!!
lav_msg_fixme <- function(...) {
  wat <- c("FIXME: ", unlist(list(...), use.names = FALSE))
  stop(lav_msg(wat), call. = FALSE, domain = NA)
}

# subroutine for above functions
lav_msg <- function(wat, txt.width = getOption("width", 80L),
                    indent = 4L, showheader = TRUE) {
  if (showheader) {
    ignore.in.stack <- c(
      "^eval$", "^try", "^doTryCatch", "^lav_msg", "^stop$", "^warning$",
      "^which$", "^unique$", "^as\\.", "^unlist$", "^message$",
      "^source$", "^withVisible$", "^tryCatch.W.E$", "^withCallingHandlers$",
      "^do.call$", "^paste"
    )
    sc <- sys.calls()
    sc.i <- length(sc)
    sc.naam <- ""
    while (sc.i > 0L) {
      x <- tryCatch(
        as.character(sc[[sc.i]][[1L]]),
        error = function(e) {"unknown"}
      )
      if (length(x) == 3L) {
        # needed if a function specified in namespace, e.g.
        # as.character(str2lang("lavaan::sem(m, d)")[[1L]])
        x <- x[[3L]]
      }
      skip <- FALSE
      for (re in ignore.in.stack) {
        if (grepl(re, x)) {
          skip <- TRUE
          break
        }
      }
      if (!skip) {
        sc.naam <- x
        break
      }
      sc.i <- sc.i - 1L
    }
    if (sc.naam == "") {
      header <- "lavaan:"
    } else {
      header <- paste0("lavaan->", sc.naam, "():")
    }
  } else {
    header <- ""
  }
  # make sure we only have a single string
  txt <- paste(wat, collapse = " ")
  # split the txt in little chunks
  chunks <- strsplit(paste(header, txt), "\\s+", fixed = FALSE)[[1]]

  # chunk size (number of characters)
  chunk.size <- nchar(chunks)

  # remove empty chunk in position 1 (if txt starts with whitespace)
  if (chunk.size[1L] == 0L) {
    chunks <- chunks[-1L]
    chunk.size <- chunk.size[-1]
  }

  nstart <- 1L
  nstop <- 1L
  corr.line1 <- 7L # first line possibly contains "error: "
  while (nstart <= length(chunks)) {
    while (nstop < length(chunks) &&
           sum(chunk.size[seq.int(nstart, nstop + 1L)]) + corr.line1 +
           nstop - nstart + indent < txt.width && chunks[nstop + 1L] != "___") {
      nstop <- nstop + 1
    }
    corr.line1 <- 0L
    if (nstop < length(chunks) && chunks[nstop + 1L] == "___") {
      # forced line break
      chunks[nstop + 1L] <-  ""
      nstop <- nstop + 1L
    }
    if (nstop < length(chunks)) {
      chunks[nstop + 1L] <- paste0(
        "\n", strrep(" ", indent),
        chunks[nstop + 1L]
      )
    }
    nstart <- nstop + 1L
    nstop <- nstart
  }
  paste(chunks, collapse = " ")
}

# Transforms a value to a character representation for use in messages
# if logsep = "array" (default), letters[1:3] is transformed to ("a", "b", "c")
# if logsep = "none", c("x", "y", "z") is transformed to "x", "y", "z"
# if logsep = "and", 1:3 is transformed to 1, 2 and 3
# if logsep = "or", c("a", "b", "c") is transformed to "a", "b" or "c"
# The type of quote can be modified via parameter qd (default = TRUE).
# If the object has names, the names will be prepended with a colon before the
# value, e.g. c(x = 2.3, y = 4) --> (x : 2.3, y : 4).
lav_msg_view <- function(x,
                         log.sep = c("array", "none", "and", "or"),
                         qd = TRUE) {
  if (missing(x)) {
    return("NULL")
  }
  log.sep <- match.arg(log.sep)
  xn <- names(x)
  if (is.list(x)) {
    xx <- sapply(x, lav_msg_view)
  } else {
    if (is.character(x)) {
      if (qd) {
        xx <- dQuote(x, q = FALSE)
      } else {
        xx <- sQuote(x, q = FALSE)
      }
    } else {
      xx <- as.character(x)
    }
    xx[is.na(x)] <- "NA"
  }
  if (!is.null(xn)) xx <- paste(xn, ":", xx)
  if (length(xx) == 1) {
    rv <- xx
  } else {
    if (log.sep == "array") rv <- paste0("(", paste(xx, collapse = ", "), ")")
    if (log.sep == "none") rv <- paste(xx, collapse = ", ")
    if (log.sep == "and") rv <- paste(paste(xx[-length(xx)], collapse = ", "), gettext("and"), xx[length(xx)])
    if (log.sep == "or") rv <- paste(paste(xx[-length(xx)], collapse = ", "), gettext("or"), xx[length(xx)])
  }
  rv
}
#  ---------------  examples of use ----------------------
# # warning if argument x is missing
#   lav_msg_warn(gettextf(
#     "argument %1$s is missing, using %2$s.",
#     x, lav_msg_view(usedvalue)
#   ))
#
# # warning if length of an argument x is greater then 1 and cannot be
#   lav_msg_warn(gettextf("%1$s argument should be a single character string.
#   Only the first one (%2$s) is used.", xname, x[[1]]))
#
# # error if argument is unknown (show value)
#   lav_msg_stop(gettextf(
#     "%1$s argument unknown: %2$s",
#     xname, lav_msg_view(xvalue)
#   ))
#
# # error if argument isn't one of the allowed values, show values allowed
#   if (length(allowed) == 2L) {
#     lav_msg_stop(gettextf(
#       "%1$s argument must be either %2$s",
#       x, lav_msg_view(allowed, "or")
#     ))
#   } else {
#     lav_msg_stop(gettextf(
#       "%1$s argument must be one of %2$s",
#       x, lav_msg_view(allowed, "or")
#     ))
#   }
#
# # error if argument isn't one of the allowed values (show invalid ones)
#   lav_msg_stop(sprintf(
#     ngettext(
#       length(invalids),
#       "invalid value in %1$s argument: %2$s.",
#       "invalid values in %1$s argument: %2$s."
#     ),
#     x, lav_msg_view(invalids, log.sep = "none")
#   ))
