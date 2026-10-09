lav_expand_strings <- function(str1, str2) {
  expanded_strings <- character(0)
  if (str1 == str2) return(expanded_strings)
  getdigit <- function(i) {if (i < 48L || i > 57L) -1 else i - 48L}
  byte1 <- as.integer(charToRaw(str1))
  byte2 <- as.integer(charToRaw(str2))
  aantal <- if (length(byte1) < length(byte2)) length(byte1) else length(byte2)
  numeric_flags <- logical(aantal)
  leading_zeroes <- integer(aantal)
  start_values <- integer(aantal)
  stop_values <- integer(aantal)
  numeric_difference <- 0L
  nb_of_items <- 1L
  j1 <- 1L
  j2 <- 1L
  while (j1 <= length(byte1) && j2 <= length(byte2)) {
    if (byte1[j1] == byte2[j2]) {
      start_values[nb_of_items] <- stop_values[nb_of_items] <- byte1[j1]
      nb_of_items <- nb_of_items + 1L
      j1 <- j1 + 1L
      j2 <- j2 + 1L
      next
    }
    if (getdigit(byte1[j1]) >= 0 && getdigit(byte2[j2]) >= 0) {
      ledz <- (getdigit(byte1[j1]) == 0 || getdigit(byte2[j2]) == 0)
      inc1 <- 1L
      sval <- getdigit(byte1[j1])
      while ((j1 + inc1) <= length(byte1) &&
              getdigit(byte1[j1 + inc1]) >= 0) {
          sval <- 10L * sval + getdigit(byte1[j1 + inc1])
          inc1 <- inc1 + 1L
      }
      start_values[nb_of_items] <- sval
      j1 <- j1 + inc1
      inc2 <- 1L
      sval <- getdigit(byte2[j2])
      while ((j2 + inc2) <= length(byte2) &&
              getdigit(byte2[j2 + inc2]) >= 0) {
          sval <- 10L * sval + getdigit(byte2[j2 + inc2])
          inc2 <- inc2 + 1L
      }
      stop_values[nb_of_items] <- sval
      j2 <- j2 + inc2
      if (ledz)
        leading_zeroes[nb_of_items] <- if (inc1 < inc2) inc2 else inc1
      numeric_flags[nb_of_items] <- TRUE
      if (numeric_difference == 0) {
          numeric_difference <- stop_values[nb_of_items] -
                                start_values[nb_of_items]
      } else if (numeric_difference !=
                (stop_values[nb_of_items] - start_values[nb_of_items])) {
          attr(expanded_strings, "error") <- 1L # distances not conform for expanding
          return(expanded_strings)
      }
      nb_of_items <- nb_of_items + 1L
    } else {
      start_values[nb_of_items] <- byte1[j1]
      stop_values[nb_of_items] <- byte2[j2]
      if (numeric_difference == 0) {
        numeric_difference <- stop_values[nb_of_items] -
                              start_values[nb_of_items]
      } else if (numeric_difference !=
            (stop_values[nb_of_items] - start_values[nb_of_items])) {
        attr(expanded_strings, "error") <- 1L # distances not conform for expanding
        return(expanded_strings)
      }
      nb_of_items <- nb_of_items + 1L
      j1 <- j1 + 1L
      j2 <- j2 + 1L
    }
  }
  if (j1 <= length(byte1) || j2 <= length(byte2)) {
    attr(expanded_strings, "error") <- 2L # string lengths not conform for expanding
    return(expanded_strings)
  }
  nb_of_items <- nb_of_items - 1L
  teken <- sign(numeric_difference)
  if (numeric_difference < -1L || numeric_difference > 1L) {
    for (i in seq(teken, numeric_difference - teken, teken)) {
      new_string_text <- ""
      for (k in seq(1, nb_of_items, 1)) {
          new_value <- start_values[k]
          if (start_values[k] != stop_values[k]) {
              new_value <- new_value + i
          }
          if (numeric_flags[k]) {
              formaat <- "%d"
              if (leading_zeroes[k]) formaat <-
                                  sprintf("%%0%dd", leading_zeroes[k])
              new_value_str <- sprintf(formaat, new_value)
              new_string_text <- paste0(new_string_text, new_value_str)
          } else {
              new_string_text <- paste0(new_string_text,
                                            rawToChar(as.raw(new_value)))
          }
      }
      expanded_strings <- c(expanded_strings, new_string_text)
    }
  }
  expanded_strings
  }

lav_expand_tokens <- function(modellist,
  modelsrc,
  types,
  ppi, # index of '++'
  starti, # index of first item of starting tokens
  stopi # index of last item of ending tokens
  ) {
  expanded_tokens <- list()
  if (ppi - starti != stopi - ppi) {
    tl <- lav_parse_txtloc(modelsrc, modellist$elem_pos[ppi])
    lav_msg_stop(gettext("tokens before and after ++ not conform"),
                   tl[1L],
                   footer = tl[2L]
    )
  }
  aantal = 0L
  for (i in seq(starti, ppi - 1L, 1L)) {
    i2 = ppi + i - starti + 1L
    if (modellist$elem_text[i] != modellist$elem_text[i2]) {
        sub_expanded_strings = lav_expand_strings(modellist$elem_text[i], modellist$elem_text[i2])
        if (!is.null(attr(sub_expanded_strings, "error", TRUE))) {
          the_error <- attr(sub_expanded_strings, "error", TRUE)
          tl <- lav_parse_txtloc(modelsrc, modellist$elem_pos[ppi])
          if (the_error == 1L) {
            lav_msg_stop(gettext("substring distances not conform for expanding"),
                    tl[1L],
                    footer = tl[2L])
          } else {
            lav_msg_stop(gettext("string lengths not conform for expanding"),
                    tl[1L],
                    footer = tl[2L])
          }
        }
        if (aantal == 0) {
            aantal = length(sub_expanded_strings)
        } else if (length(sub_expanded_strings) != aantal) {
          tl <- lav_parse_txtloc(modelsrc, modellist$elem_pos[ppi])
          lav_msg_stop(gettext("substring distances not conform for expanding"),
                   tl[1L],
                   footer = tl[2L]
          )
        }
    }
  }
  expnum <- 1L
  expanded_tokens$elem_text[expnum] <- "+"
  expanded_tokens$elem_type[expnum] <- types$symbol
  expanded_tokens$elem_pos[expnum] <- modellist$elem_pos[ppi]
  expanded_tokens$elem_formula_number[expnum] <- modellist$elem_formula_number[ppi]
  expnum <- expnum + 1L
  if (aantal > 0L) {
    for (j in seq(1L, aantal, 1L)) {
      for (i in seq(starti, ppi - 1L, 1L)) {
          i2 = ppi + i - starti + 1
          if (modellist$elem_text[i] != modellist$elem_text[i2]) {
              sub_expanded_strings = lav_expand_strings(modellist$elem_text[i], modellist$elem_text[i2])
              expanded_tokens$elem_text[expnum] <- sub_expanded_strings[j]
          } else {
              expanded_tokens$elem_text[expnum] <- modellist$elem_text[i]
          }
          expanded_tokens$elem_type[expnum] <- modellist$elem_type[i]
          expanded_tokens$elem_pos[expnum] <- modellist$elem_pos[i]
          expanded_tokens$elem_formula_number[expnum] <- modellist$elem_formula_number[ppi]
          expnum <- expnum + 1L
      }
      expanded_tokens$elem_text[expnum] <- "+"
      expanded_tokens$elem_type[expnum] <- types$symbol
      expanded_tokens$elem_pos[expnum] <- modellist$elem_pos[ppi]
      expanded_tokens$elem_formula_number[expnum] <- modellist$elem_formula_number[ppi]
      expnum <- expnum + 1L
    }
  }
  return(expanded_tokens)
  }

lav_expand_plusplus <- function(modellist, modelsrc, types) {
    ppi = -1L
    for (i in seq_along(modellist$elem_type)) {
        if (modellist$elem_text[i] == "++") {
            ppi = i
            break
        }
    }
    if (ppi == -1) return(modellist)  # no '++' to handle
    newvec <- modellist
    while (ppi > -1) {
      # find starting item
      starti = ppi - 1L
      balance = 0L
      while (starti > 0L && ((
          newvec$elem_text[starti] != "+" &&
          newvec$elem_type[starti] != types$newline &&
          newvec$elem_type[starti] != types$lavaanoperator) ||
          balance != 0L)
        ) {
          if (newvec$elem_text[starti] == "(") balance <- balance + 1L
          if (newvec$elem_text[starti] == ")") balance <- balance - 1L
          starti <- starti - 1L
        }
        starti <- starti + 1L
        # find ending item
        balance = 0L
        stopi = ppi + 1L
        while (stopi <= length(newvec$elem_type) && ((
          newvec$elem_text[stopi] != "+" &&
          newvec$elem_text[stopi] != "++" &&
          newvec$elem_type[stopi] != types$newline &&
          newvec$elem_type[stopi] != types$lavaanoperator) ||
          balance != 0L)
      ) {
          if (newvec$elem_text[stopi] == "(") balance <- balance + 1L
          if (newvec$elem_text[stopi] == ")") balance <- balance - 1L
          stopi <- stopi + 1L
      }
      stopi <- stopi - 1L
      # get added tokens
      toinsert = lav_expand_tokens(newvec, modelsrc, types, ppi, starti, stopi)
      rv <- list()
      rvi <- 1L
      for (i in seq(1L, ppi - 1L, 1L)) {
        rv$elem_text[rvi] <- newvec$elem_text[i]
        rv$elem_type[rvi] <- newvec$elem_type[i]
        rv$elem_pos[rvi] <- newvec$elem_pos[i]
        rv$elem_formula_number[rvi] <- newvec$elem_formula_number[i]
        rvi <- rvi + 1L
      }
      for (i in seq_along(toinsert$elem_type)) {
        rv$elem_text[rvi] <- toinsert$elem_text[i]
        rv$elem_type[rvi] <- toinsert$elem_type[i]
        rv$elem_pos[rvi] <- toinsert$elem_pos[i]
        rv$elem_formula_number[rvi] <- toinsert$elem_formula_number[i]
        rvi <- rvi + 1L
      }
      for (i in seq(ppi + 1L, length(newvec$elem_typ), 1L)) {
        rv$elem_text[rvi] <- newvec$elem_text[i]
        rv$elem_type[rvi] <- newvec$elem_type[i]
        rv$elem_pos[rvi] <- newvec$elem_pos[i]
        rv$elem_formula_number[rvi] <- newvec$elem_formula_number[i]
        rvi <- rvi + 1L
      }
      # check for more '++'s
      ppi = -1
      for (i in seq_along(rv$elem_type)) {
        if (rv$elem_text[i] == "++") {
            ppi = i
            break
        }
      }
      newvec <- rv
    }
    rv
}
