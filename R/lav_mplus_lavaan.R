# this code is written by Michael Hallquist
# First draft of parser to convert Mplus model syntax to lavaan model syntax

# idea: build parTable and run model from mplus syntax
# then perhaps write export function: parTable2Mplus
# and/or parTable2lavaan

lav_mplus_trim <- function(string) {
  string_trim <- sapply(string, function(x) {
    x <- sub("^\\s*", "", x, perl = TRUE)
    x <- sub("\\s*$", "", x, perl = TRUE)
    x
  }, USE.NAMES = FALSE)
  string_trim
}

# small utility function to join strings in a regexp loop
lav_mplus_join_regex <- function(cmd, arg_expand, matches,
                  iterator, match_length = "match.length") {
  if (iterator == 1 && matches[iterator] > 1) {
    pre <- substr(cmd, 1, matches[iterator] - 1)
  } else {
    pre <- ""
  }

  # if this is not the final match, then get sub-string between the end of
  #                               this match and the beginning of the next
  # otherwise, match to the end of the command
  post_end <- ifelse(iterator < length(matches),
               matches[iterator + 1] - 1, nchar(cmd))
  post <- substr(cmd, matches[iterator] + attr(matches,
               match_length)[iterator], post_end)

  cmd_expand <- paste(pre, arg_expand, post, sep = "")
  cmd_expand
}

# expand Mplus hyphen syntax (will also expand constraints with hyphens)
lav_mplus_expand_cmd <- function(cmd, alpha_start = TRUE) {
  # use negative lookahead and negative lookbehind to eliminate possibility of
  #  hyphen being used as a negative starting value (e.g., x*-1)
  # also avoid match of anything that includes a decimal point, such as a
  #  floating-point starting value -10.5*x1

  # if alphaStart==TRUE, then require that the matches before and after hyphens
  #  begin with alpha character
  # this is used for variable names, whereas the more generic expansion works
  #  for numeric constraints and such

  # need to do a better job of this so that u1-u20* is supported...
  #         I don't think the regexp below is general enough

  # if (alphaStart) {
  #  hyphens <- gregexpr("[_A-Za-z]+\\w*\\s*-\\s*[_A-Za-z]+\\w*",
  #                      cmd, perl=TRUE)[[1]]
  # } else {
  #  hyphens <- gregexpr(
  # "(?!<(\\*|\\.))\\w+(?!(\\*|\\.))\\s*-\\s*(?!<(\\*|\\.))\\w+(?!(\\*|\\.))",
  #  cmd, perl=TRUE)[[1]]
  # }

  # hyphens <- gregexpr(
  # "(?!<(\\*|\\.))\\w+(?!(\\*|\\.))\\s*-\\s*(?!<(\\*|\\.))\\w+(?!(\\*|\\.))",
  #  cmd, perl=TRUE)[[1]]

  # support trailing @XXX. Still still fail on Trait1-Trait3*XXX
  hyphens <- gregexpr(
    "(?!<(\\*|\\.))\\w+(?!(\\*|\\.))\\s*-\\s*(?!<(\\*|\\.))\\w+(?!(\\*|\\.))(@[\\d\\.\\-]+)?", # nolint
     cmd, perl = TRUE)[[1]]

  # Promising, but this is still failing in the case of x3*1 -4.25*x4
  # On either side of a hyphen, require alpha character followed
  # by alphanumeric
  # This enforces that neither side of the hyphen can be a number
  # Alternatively, match digits on either side alone
  # hyphens <- gregexpr(
  # "([A-z]+\\w*\\s*-\\s*[A-z]+\\w*(@[\\d\\.-]+)?|\\d+\\s*-\\s*\\d+)",
  #  cmd, perl=TRUE)[[1]]

  if (hyphens[1L] > 0) {
    cmd_expand <- c()
    ep <- 1

    for (v in seq_along(hyphens)) {
      # match one keyword before and after hyphen
      argsplit <- strsplit(substr(cmd, hyphens[v], hyphens[v] +
        attr(hyphens, "match.length")[v] - 1), "\\s*-\\s*", perl = TRUE)[[1]]

      v_pre <- argsplit[1]
      v_post <- argsplit[2]

      v_post_suffix <-                   # will be empty string if not present
        sub("^([^@]+)(@[\\d\\-.]+)?$", "\\2", v_post, perl = TRUE)
      v_post <- sub("@[\\d\\-.]+$", "", v_post, perl = TRUE) # trim @ suffix

      # If v_pre and v_post contain leading alpha characters,
      #   verify that these prefixes match.
      # Otherwise, there is nothing to expand, as in the case of
      #   MODEL CONSTRAINT: e1e2=e1-e2_n.
      v_pre_alpha <- sub("\\d+$", "", v_pre, perl = TRUE)
      v_post_alpha <- sub("\\d+$", "", v_post, perl = TRUE)

      # only enforce prefix match if we have leading alpha characters
      #  (i.e., not simple numeric 1 - 3 syntax)
      if (length(v_pre_alpha) > 0L && length(v_post_alpha) > 0L) {
        # if alpha prefixes do match, assume that the hyphen is not for
        #  expansion (e.g., in subtraction case)
        if (v_pre_alpha != v_post_alpha) {
          return(cmd)
        }
      }

      # the basic positive lookbehind blows up with pure numeric
      #   constraints (1 - 3) because no alpha char precedes digit
      # can use an non-capturing alternation grouping to allow for digits only
      #  or the final digits after alphas (as in v_post.num)
      v_pre_num <- as.integer(sub(
        "\\w*(?<=[A-Za-z_])(\\d+)$", "\\1", v_pre, perl = TRUE))
         # use positive lookbehind to avoid greedy \w+
         # match -- capture all digits

      v_post_match <- regexpr("^(?:\\w*(?<=[A-Za-z_])(\\d+)|(\\d+))$",
                      v_post, perl = TRUE)
      stopifnot(v_post_match[1L] > 0)

      # match mat be under capture[1] or capture[2]
      # because of alternation above
      which_capture <- which(attr(v_post_match, "capture.start") > 0)

      v_post_num <- as.integer(substr(v_post, attr(v_post_match,
         "capture.start")[which_capture],
         attr(v_post_match, "capture.start")[which_capture] +
           attr(v_post_match, "capture.length")[which_capture] - 1))
      v_post_prefix <- substr(v_post, 1, attr(v_post_match,
         "capture.start")[which_capture] - 1)
         # just trusting that pre and post match

      if (is.na(v_pre_num) || is.na(v_post_num)) lav_msg_stop(
        gettext("Cannot expand variables:"), v_pre, ", ", v_post)
      v_expand <- paste(v_post_prefix, v_pre_num:v_post_num, v_post_suffix,
         sep = "", collapse = " ")

      # for first hyphen, there may be non-hyphenated syntax
      #    preceding the initial match
      cmd_expand[ep] <- lav_mplus_join_regex(cmd, v_expand, hyphens, v)

      # This won't really work because the cmd.expand element may contain
      # other variables that are at the beginning or end, prior to hyphen stuff
      # This is superseded by logic above where @ is included in hyphen match,
      # then trapped as suffix I don't think it will work yet for this Mplus
      # syntax: y1-y10*5 -- the 5 wouldn't propagate
      # handle the case of @ fixed values or * starting values used in a list
      # example: Trait1-Trait3@1
      ## if (grepl("@|\\*", cmd.expand[ep], perl=TRUE)) {
      ##   exp_split <- strsplit(cmd.expand[ep], "\\s+", perl=TRUE)[[1]]
      ##   suffixes <- sub("^([^@\\*]+)([@*][\\d\\.-]+)?$", "\\2",
      #                exp_split, perl=TRUE)
      ##   variables <- sub("^([^@\\*]+)([@*][\\d\\.-]+)?$", "\\1",
      #                exp_split, perl=TRUE)
      ##   suffixes <- suffixes[suffixes != ""]
      ##   if (length(unique(suffixes)) > 1L) {
      ##     browser()

      ##     #stop("Don't know how to interpret syntax: ", cmd)
      ##   } else {
      ##     variables <- paste0(variables, suffixes[1])
      ##     cmd.expand[ep] <- paste(variables, collapse=" ")
      ##   }
      ## }

      ep <- ep + 1
    }
    paste(cmd_expand, collapse = "")
  } else {
    cmd # no hyphens to expand
  }
}


# handle starting values and fixed parameters on rhs
lav_mplus_cmd_fix_start <- function(cmd) {
  cmd_parse <- c()
  ep <- 1L

  # support ESEM-like syntax: F BY a1* a2*
  # The easy path: putting in 1s before we proceed on parsing
  # Mar2023 bugfix: support parenthesis after * in case a parameter
  #                      constraint comes next
  cmd <- gsub("([A-z]+\\w*)\\s*\\*(?=\\s+\\(?[A-z]+|\\s*$)", "\\1*1",
             cmd, perl = TRUE)

  if ((fixed_starts <- gregexpr("[\\w\\.\\-$]+\\s*([@*])\\s*[\\w\\.\\-]+",
                                cmd, perl = TRUE)[[1]])[1L] > 0) {
       # shouldn't it be \\*, not * ?! Come back to this.
    for (f in seq_along(fixed_starts)) {
      # capture above obtains the fixed/start character (@ or *),
      # whereas match obtains the full regex match
      opchar <- substr(cmd, attr(fixed_starts, "capture.start")[f],
      attr(fixed_starts, "capture.start")[f] + attr(fixed_starts,
         "capture.length")[f] - 1)

      # match arguments around asterisk/at symbol
      argsplit <- strsplit(substr(cmd, fixed_starts[f], fixed_starts[f] +
        attr(fixed_starts, "match.length")[f] - 1), paste0("\\s*",
        ifelse(opchar == "*", "\\*", opchar), "\\s*"), perl = TRUE)[[1]]
      v_pre <- argsplit[1]
      v_post <- argsplit[2]

      if (suppressWarnings(is.na(as.numeric(v_pre)))) {
        # fixed.starts value post-multiplier
        var <- v_pre
        val <- v_post
      } else if (suppressWarnings(is.na(as.numeric(v_post)))) {
        # starting value pre-multiplier
        var <- v_post
        val <- v_pre
      } else {
        lav_msg_stop(
          gettext("Cannot parse Mplus fixed/starts values specification:"),
          v_pre, v_post)
      }

      if (opchar == "@") {
        cmd_parse[ep] <- lav_mplus_join_regex(cmd,
          paste0(val, "*", var, sep = ""), fixed_starts, f)
        ep <- ep + 1L
      } else {
        cmd_parse[ep] <- lav_mplus_join_regex(cmd,
          paste0("start(", val, ")*", var, sep = ""), fixed_starts, f)
        ep <- ep + 1L
      }
    }
    paste(cmd_parse, collapse = "")
  } else {
    cmd
  }
}

lav_mplus_cmd_constraints <- function(cmd) {
  # Allow cmd to have newlines embedded. In this case, split on newlines, and
  #  loop over and parse each chunk
  # Dump leading and trailing newlines, which contain no information about
  #  constraints, but may add dummy elements to vector after strsplit
  # Maybe return LHS and RHS parsed command where constraints only appear on
  #  the RHS, whereas the LHS contains only parameters.
  # Example: LHS is v1 v2 v3 and RHS is con1*v1 con2*v2 con3*v3

  cmd_split <- strsplit(cmd, "\n")[[1]]

  # drop empty lines (especially leading newline)
  cmd_split <- if (length(empty_pos <- which(cmd_split == "")) > 0L) {
    cmd_split[-1 * empty_pos]
  } else {
    cmd_split
  }

  # Create a version of the command with no modifiers (constraints,
  #   starting values, etc.) specifications.
  # This is useful for syntax that uses the params on the LHS and with a
  #   modified RHS. Example: v1 ~~ conB*v1
  cmd_nomodifiers <- paste0(gsub("(start\\([^\\)]+\\)\\*|[\\d\\-\\.]+\\*)",
   "", cmd_split, perl = TRUE), collapse = " ") # peel off premultiplication
  cmd_nomodifiers <- gsub("\\([^\\)]+\\)", "", cmd_nomodifiers, perl = TRUE)

  cmd_tojoin <- c()
  # will store all chunks divided by newlines, which will be joined at the end.

  # iterate over each newline segment
  for (n in seq_along(cmd_split)) {
    # in principle, now that we respect newlines, parens should only be of
    # length 1, since Mplus syntax dictates newlines for each use of
    # parentheses for constraints
    if ((parens <- gregexpr("(?<!start)\\(([^\\)]+)\\)", cmd_split[n],
     perl = TRUE)[[1L]])[1L] > 0) { # match parentheses, but not start()
      # the syntax chunk after all parentheses have been matched
      cmd_expand <- c()

      for (p in seq_along(parens)) {
        # string within the constraint parentheses
        constraints <- substr(cmd_split[n], attr(parens, "capture.start")[p],
         attr(parens,
          "capture.start")[p] + attr(parens, "capture.length")[p] - 1)

        # Divide constraints on spaces to determine number of constraints to
        # parse. Use lav_mplus_trim to avoid problem of user including
        # leading/trailing spaces within parentheses.
        con_split <- strsplit(lav_mplus_trim(constraints), "\\s+",
         perl = TRUE)[[1]]

        # if Mplus uses a purely numeric constraint, then add ".con" prefix
        #   to be consistent with R naming.
        con_split <- sapply(con_split, function(x) {
          if (!suppressWarnings(is.na(as.numeric(x)))) {
            make.names(paste0(".con", x))
          } else {
            x
          }
        })

        # determine the parameters that precede the parentheses (either first
        #  character for p == 1 or character after preceding parentheses)
        prestr_start <- ifelse(p > 1, attr(parens, "capture.start")[p - 1] +
          attr(parens, "capture.length")[p - 1] + 1, 1)

        # obtain the parameters that precede the parentheses, divide into
        #  arguments on spaces
        # use lav_mplus_trim here because first char after prestrStart for
        #  p > 1 will probably be a space
        precmd_split <- strsplit(lav_mplus_trim(substr(cmd_split[n],
           prestr_start, parens[p] - 1)), "\\s+", perl = TRUE)[[1]]

        # peel off any potential LHS arguments, such as F1 BY
        precmd_lhsop <- which(tolower(precmd_split) %in% c("by", "with", "on"))
        if (any(precmd_lhsop)) {
          lhsop <- paste0(precmd_split[1:precmd_lhsop[1L]], " ",
                          collapse = " ")
          # join lhs and op as a single string, add trailing space so that
          #  paste with expanded RHS is right.
          rhs <- precmd_split[(precmd_lhsop + 1):length(precmd_split)]
        } else {
          lhsop <- ""
          rhs <- precmd_split
        }

        if (length(con_split) > 1L) {
          # several constraints listed within parentheses.
          #  Example: F1 BY X1 X2 X3 X4 (C2 C3 C4)
          # thus, backwards match the constraints to parameters

          # restrict parameters to backwards match to be of the same length
          #  as number of constraints
          rhs_backmatch <- rhs[(length(rhs) -
                           length(con_split) + 1):length(rhs)]

          rhs_expand <- c()

          # check that no mean or scale markers are part of the rhs param
          #  to expand
          if ((pre_mark_match <- regexpr("^\\s*[\\[\\{]", rhs_backmatch[1L],
                  perl = TRUE))[1L] > 0) {
            pre_mark <- substr(rhs_backmatch[1L], pre_mark_match[1L],
              pre_mark_match[1L] + attr(pre_mark_match,
                "match.length")[1L] - 1)
            rhs_backmatch[1L] <- substr(rhs_backmatch[1L],
              pre_mark_match[1L] + attr(pre_mark_match,
                 "match.length")[1L], nchar(rhs_backmatch[1L]))
          } else {
            pre_mark <- ""
          }

          if ((post_mark_match <- regexpr("[\\]\\}]\\s*$",
             rhs_backmatch[length(rhs_backmatch)], perl = TRUE))[1L] > 0) {
            post_mark <- substr(rhs_backmatch[length(rhs_backmatch)],
             post_mark_match[1L], nchar(rhs_backmatch[length(rhs_backmatch)]))
            rhs_backmatch[length(rhs_backmatch)] <-
              substr(rhs_backmatch[length(rhs_backmatch)], 1,
                       post_mark_match[1L] - 1)
          } else {
            post_mark <- ""
          }


          # pre-multiply each parameter with each corresponding constraint
          for (i in seq_along(rhs_backmatch)) {
            rhs_expand[i] <- paste0(con_split[i], "*", rhs_backmatch[i])
          }

          # join rhs as string and add back in mean/scale operator, if present
          rhs_expand <- paste0(pre_mark, paste(rhs_expand,
             collapse = " "), post_mark)

          # if there were params that preceded the backwards match,
          #  then add these back to the syntax
          # append this syntax to the parsed command, cmd.expand
          if (length(rhs) - length(con_split) > 0L) {
            cmd_expand <- c(cmd_expand, paste(lhsop, paste(rhs[1:(length(rhs) -
                            length(con_split))], collapse = " "), rhs_expand))
          } else {
            cmd_expand <- c(cmd_expand, paste0(lhsop, rhs_expand))
          }
        } else {
          # should be able to reduce redundancy with above

          # all parameters on the right hand side are to be equated
          # thus, pre-multiply each parameter by the constraint

          # check that no mean or scale markers are part of the
          #   rhs param to expand
          # DUPE CODE FROM ABOVE. Make Function?!
          if ((pre_mark_match <-
            regexpr("^\\s*[\\[\\{]", rhs[1L], perl = TRUE))[1L] > 0) {
            pre_mark <- substr(rhs[1L], pre_mark_match[1L], pre_mark_match[1L] +
              attr(pre_mark_match, "match.length")[1L] - 1)
            rhs[1L] <- substr(rhs[1L], pre_mark_match[1L] + attr(pre_mark_match,
               "match.length")[1L], nchar(rhs[1L]))
          } else {
            pre_mark <- ""
          }

          if ((post_mark_match <-
             regexpr("[\\]\\}]\\s*$", rhs[length(rhs)], perl = TRUE))[1L] > 0) {
            post_mark <- substr(rhs[length(rhs)], post_mark_match[1L],
             nchar(rhs[length(rhs)]))
            rhs[length(rhs)] <- substr(rhs[length(rhs)], 1,
             post_mark_match[1L] - 1)
          } else {
            post_mark <- ""
          }


          rhs_expand <- c()
          for (i in seq_along(rhs)) {
            rhs_expand[i] <- paste0(con_split[1L], "*", rhs[i])
          }

          # join rhs as string
          rhs_expand <- paste0(pre_mark,
            paste(rhs_expand, collapse = " "), post_mark)

          cmd_expand <- c(cmd_expand, paste0(lhsop, rhs_expand))
        }
      }

      cmd_tojoin[n] <- paste(cmd_expand, collapse = " ")
    } else {
      cmd_tojoin[n] <- cmd_split[n]
    } # no parens
  }

  # eliminate newlines in this function so that they
  # don't mess up \\s+ splits downstream
  to_return <- paste(cmd_tojoin, collapse = " ")
  attr(to_return, "noModifiers") <- cmd_nomodifiers

  to_return
}

lav_mplus_cmd_growth <- function(cmd) {
  # can assume that any spaces between tscore and variable were stripped by
  #  lav_mplus_cmd_fix_start

  # verify that this is not a random slope
  if (any(tolower(strsplit(cmd, "\\s+", perl = TRUE)[[1]]) %in%
                   c("on", "at"))) {
    lav_msg_stop(gettext(
      "lavaan does not support random slopes or individually varying
      growth model time scores"))
  }

  cmd_split <- strsplit(cmd, "\\s*\\|\\s*", perl = TRUE)[[1]]
  if (!length(cmd_split) == 2) {
    lav_msg_stop(gettext("Unknown growth syntax:"), cmd)
  }

  lhs <- cmd_split[1]
  lhs_split <- strsplit(lhs, "\\s+", perl = TRUE)[[1]]

  rhs <- cmd_split[2]
  rhs_split <- strsplit(rhs, "(\\*|\\s+)", perl = TRUE)[[1]]

  if (length(rhs_split) %% 2 != 0) {
    lav_msg_stop(gettext(
      "Number of variables and number of tscores does not match:"), rhs)
  }
  tscores <- as.numeric(rhs_split[seq_along(rhs_split) %% 2 != 0])
             # pre-multipliers

  vars <- rhs_split[seq_along(rhs_split) %% 2 == 0]

  cmd_expand <- c()

  for (p in 0:(length(lhs_split) - 1)) {
    if (p == 0) {
      # intercept
      cmd_expand <- c(cmd_expand, paste(lhs_split[(p + 1)], "=~",
       paste("1*", vars, sep = "", collapse = " + ")))
    } else {
      cmd_expand <- c(cmd_expand, paste(lhs_split[(p + 1)], "=~",
       paste(tscores^p, "*", vars, sep = "", collapse = " + ")))
    }
  }

  cmd_expand
}

# function to wrap long lines at a certain width,
#  splitting on + symbols to be consistent with R syntax
lav_mplus_cmd_wrap <- function(cmd, width = 90, exdent = 5) {
  result <- lapply(cmd, function(line) {
    if (nchar(line) > width) {
      split <- c()
      spos <- 1L

      plus_match <- gregexpr("+", line, fixed = TRUE)[[1]]
#      mpos <- 1L

      if (plus_match[1L] > 0L) {
        # split after plus symbol
        chars_remain <- nchar(line)
        while (chars_remain > 0L) {
#          to_process <- substr(line, nchar(line) - chars_remain + 1,
#           nchar(line))
          offset <- nchar(line) - chars_remain + 1

          if (nchar(remainder <- substr(line, offset, nchar(line))) <=
                              (width - exdent)) {
            # remainder of line fits within width -- no need
            #            to continue wrapping
            split[spos] <- remainder
            chars_remain <- 0
          } else {
            wrap_at <- which(plus_match < (width + offset - exdent))
            wrap_at <- wrap_at[length(wrap_at)] # at the final +

            split[spos] <- substr(line, offset, plus_match[wrap_at])
            chars_remain <- chars_remain - nchar(split[spos])
            spos <- spos + 1
          }
        }

        # remove leading and trailing chars
        split <- lav_mplus_trim(split)

        # handle exdent
        split <- sapply(seq_along(split), function(x) {
          if (x > 1) {
            paste0(paste(rep(" ", exdent), collapse = ""), split[x])
          } else {
            split[x]
          }
        })

        split
      } else {
        strwrap(line, width = width, exdent = exdent)
           # convention strwrap when no + present
      }
    } else {
      line
    }
  })

  # bind together multi-line expansions into single vector
  unname(do.call(c, result))
}

lav_mplus_syntax_constraints <- function(syntax) {
  # should probably pass in model syntax along with some tracking
  #  of which parameter labels are defined.

  # convert MODEL CONSTRAINT section to lavaan model syntax
  syntax <- paste(lapply(lav_mplus_trim(strsplit(syntax, "\n")), function(x) {
    if (length(x) == 0L && is.character(x)) "" else x
  }), collapse = "\n")

  # replace ! with # for comment lines.
  #  Also strip newline and replace with semicolon
  syntax <- gsub("(\\s*)!(.+)\n", "\\1#\\2;", syntax, perl = TRUE)

  # split into vector of strings
  # need to peel off leading or trailing newlines -- leads to parsing
  #  confusion downstream otherwise
  syntax_split <- gsub("(^\n|\n$)", "",
                unlist(strsplit(syntax, ";")), perl = TRUE)

  constraint_out <- c()

  # TODO: Handle PLOT and LOOP syntax for model constraints.
  # TODO: Handle DO loop convention

  # first parse new parameters defined in MODEL CONSTRAINT into a vector
  new_parameters <- c() # parameters that are defined in constraint section
  if (length(new_con_lines <- grep("^\\s*NEW\\s*\\([^\\)]+\\)",
        syntax_split, perl = TRUE, ignore.case = TRUE)) > 0L) {
    for (cmd in syntax_split[new_con_lines]) {
      # process new constraint definition
      new_con <- regexpr("^\\s*NEW\\s*\\(([^\\)]+)\\)", cmd,
        perl = TRUE, ignore.case = TRUE)
      if (new_con[1L] == -1)
        lav_msg_stop(gettext("Unable to parse names of new contraints"))
      new_con <- substr(cmd, attr(new_con, "capture.start"),
       attr(new_con, "capture.start") + attr(new_con, "capture.length") - 1L)
      new_con <- lav_mplus_expand_cmd(new_con) # allow for hyphen expansion
      new_parameters <- c(new_parameters, strsplit(lav_mplus_trim(new_con),
          "\\s+", perl = TRUE)[[1L]])
    }

    syntax_split <- syntax_split[-1L * new_con_lines] # drop out these lines
    parameters_undefined <- new_parameters
    # to be used below to handle ambiguity of equation versus definition
  }

  for (cmd in syntax_split) {
    if (grepl("^\\s*#", cmd, perl = TRUE)) { # comment line
      constraint_out <- c(constraint_out, gsub("\n", "", cmd, fixed = TRUE))
      # drop any newlines
    } else if (grepl("^\\s+$", cmd, perl = TRUE)) {
      # do nothing, just a space line
    } else {
      # constraint proper
      cmd <- gsub("**", "^", cmd, fixed = TRUE) # handle exponent

      # lower case the math operations supported by Mplus to be
      #        consistent with R
      # match all math operators, then lower case each and rejoin string
      maths <- gregexpr(
        "(SQRT|LOG|LOG10|EXP|ABS|SIN|COS|TAN|ASIN|ACOS|ATAN)\\s*\\(",
         cmd, perl = TRUE)[[1L]]
      if (maths[1L] > 0) {
        maths_replace <- c()
        ep <- 1

        for (i in seq_along(maths)) {
          operator <- tolower(substr(cmd, attr(maths, "capture.start")[i],
           attr(maths, "capture.start")[i] +
             attr(maths, "capture.length")[i] - 1))
          maths_replace[ep] <- lav_mplus_join_regex(cmd, operator, maths, i,
             match_length = "capture.length")
          # only match operator, not opening (
          ep <- ep + 1
        }
        cmd <- paste(maths_replace, collapse = "")
      }

      # equating some lhs and rhs: could reflect definition of new parameter
      if ((equals <- regexpr("=", cmd, fixed = TRUE))[1L] > 0) {
        lhs <- lav_mplus_trim(substr(cmd, 1, equals - 1))
        rhs <- lav_mplus_trim(substr(cmd, equals +
          attr(equals, "match.length"), nchar(cmd)))

        # possibility of lhs or rhs containing the single variable to be equated
        if (regexpr("\\s+", lhs, perl = TRUE)[1L] > 0L) {
          def <- rhs
          body <- lhs
        } else if (regexpr("\\s+", rhs, perl = TRUE)[1L] > 0L) {
          def <- lhs
          body <- rhs
        } else {
          # warning("Can't figure out which side of constraint
          #         defines a parameter")
          # this would occur for simple rel5 = rel2 sort of syntax
          def <- lhs
          body <- rhs
        }

        # must decide whether this is a new parameter (:=)
        #      or equation of exising labels (==)
        # alternatively, could be zero, as in  0 = x + y
        # this is tricky, because mplus doesn't differentiate
        #  definition from equation
        # consequently, could confuse the issue as in ex5.20
        # NEW(rel2 rel5 stan3 stan6);
        # rel2 = lam2**2*vf1/(lam2**2*vf1 + ve2);
        # rel5 = lam5**2*vf2/(lam5**2*vf2 + ve5);
        # rel5 = rel2;

        # for now, only define a new constraint if it's not already defined
        # otherwise equate
        if (def %in% new_parameters && def %in% parameters_undefined) {
          constraint_out <- c(constraint_out, paste(def, ":=", body))
          parameters_undefined <-
            parameters_undefined[!parameters_undefined == def]
        } else {
          constraint_out <- c(constraint_out, paste(def, "==", body))
        }
      } else {
        # inequality constraints -- paste as is
        constraint_out <- c(constraint_out, cmd)
      }
    }
  }

  wrap <- paste(lav_mplus_cmd_wrap(constraint_out, width = 90, exdent = 5),
     collapse = "\n")
  wrap
}

lav_mplus_syntax_model <- function(syntax) {
  if (is.character(syntax)) {
    if (length(syntax) > 1L) {
      syntax <- paste(syntax, collapse = "\n")
    } # concatenate into a long string separated by newlines
  } else {
    lav_msg_stop(gettext(
      "lav_mplus_syntax_model accepts a single character string or
      character vector containing all model syntax"))
  }

  # because this is now exposed as a function in the package, handle the case of
  #   the user passing in full .inp file as text
  # we should only be interested in the MODEL and MODEL CONSTRAINT sections
  by_line <- strsplit(syntax, "\r?\n", perl = TRUE)[[1]]
  input_headers <- grep(
    "^\\s*(title:|data.*:|variable:|define:|analysis:|model.*:|output:|savedata:|plot:|montecarlo:)", # nolint
     by_line, ignore.case = TRUE, perl = TRUE)
  con_syntax <- c()
  if (length(input_headers) > 0L) {
    # warning("lav_mplus_syntax_model is intended to accept only the model
    #  section, not an entire .inp file. For the .inp file case,
    #  use lav_mplus_lavaan")
    parsed_syntax <- lav_mplus_text_sections(by_line, "local")

    # handle model constraint
    if ("model.constraint" %in% names(parsed_syntax)) {
      con_syntax <- strsplit(lav_mplus_syntax_constraints(
                  parsed_syntax$model.constraint), "\n")[[1]]
    }

    # just keep model syntax before continuing
    syntax <- parsed_syntax$model
  }

  # initial strip of leading/trailing whitespace, which can interfere with
  #   splitting on spaces
  # strsplit generates character(0) for empty strings, which causes problems in
  #   paste because paste actually includes it as a literal
  # example: paste(list(character(0), "asdf", character(0)), collapse=" ")
  # thus, use lapply to convert these to empty strings first
  syntax <- paste(lapply(lav_mplus_trim(strsplit(syntax, "\n")), function(x) {
    if (length(x) == 0L && is.character(x)) "" else x
  }), collapse = "\n")

  # replace ! with # for comment lines. Also strip newline and replace
  #  with semicolon
  syntax <- gsub("(\\s*)!(.+)\n*", "\\1#\\2;", syntax, perl = TRUE)

  # new direction: retain newlines in parsed syntax until after
  #  constraints have been parsed

  # delete newlines
  # syntax <- gsub("\n", "", syntax, fixed=TRUE)

  # replace semicolons with newlines prior to split (divide into commands)
  # syntax <- gsub(";", "\n", syntax, fixed=TRUE)

  # split into vector of strings
  # syntax.split <- unlist( strsplit(syntax, "\n") )
  syntax_split <- lav_mplus_trim(unlist(strsplit(syntax, ";")))

  # format of parTable to mimic.
  # 'data.frame':	34 obs. of  12 variables:
  #  $ id    : int  1 2 3 4 5 6 7 8 9 10 ...
  #  $ lhs   : chr  "ind60" "ind60" "ind60" "dem60" ...
  #  $ op    : chr  "=~" "=~" "=~" "=~" ...
  #  $ rhs   : chr  "x1" "x2" "x3" "y1" ...
  #  $ user  : int  1 1 1 1 1 1 1 1 1 1 ...
  #  $ group : int  1 1 1 1 1 1 1 1 1 1 ...
  #  $ free  : int  0 1 2 0 3 4 5 0 6 7 ...
  #  $ ustart: num  1 NA NA 1 NA NA NA 1 NA NA ...
  #  $ exo   : int  0 0 0 0 0 0 0 0 0 0 ...
  #  $ label : chr  "" "" "" "" ...
  #  $ eq.id : int  0 0 0 0 0 0 0 0 0 0 ...
  #  $ unco  : int  0 1 2 0 3 4 5 0 6 7 ...

  # vector of lavaan syntax
  lavaan_out <- c()

  for (cmd in syntax_split) {
    if (grepl("^\\s*#", cmd, perl = TRUE)) { # comment line
      lavaan_out <- c(lavaan_out, gsub("\n", "", cmd, fixed = TRUE))
       # drop any newlines (otherwise done by lav_mplus_cmd_constraints)
    } else if (grepl("^\\s*$", cmd, perl = TRUE)) {
      # do nothing, just a space or blank line
    } else {
      # hyphen expansion
      cmd <- lav_mplus_expand_cmd(cmd)

      # parse fixed parameters and starting values
      cmd <- lav_mplus_cmd_fix_start(cmd)

      # parse any constraints here (avoid weird logic below)
      cmd <- lav_mplus_cmd_constraints(cmd)

      if ((op <- regexpr("\\s+(by|on|with|pwith)\\s+", cmd, ignore.case = TRUE,
       perl = TRUE))[1L] > 0) { # regressions, factors, covariances

        lhs <- substr(cmd, 1, op - 1)
        # using op takes match.start which will omit spaces before operator
        rhs <- substr(cmd, op + attr(op, "match.length"), nchar(cmd))
        operator <- tolower(substr(cmd, attr(op, "capture.start"),
           attr(op, "capture.start") + attr(op, "capture.length") - 1))

        if (operator == "by") {
          lav_operator <- "=~"
        } else if (operator == "with" || operator == "pwith") {
          lav_operator <- "~~"
        } else if (operator == "on") {
          lav_operator <- "~"
        }

        # handle parameter combinations
        lhs_split <- strsplit(lhs, "\\s+")[[1]] # lav_mplus_trim(

        # handle pwith syntax
        if (operator == "pwith") {
          # TODO: Figure out if pwith can be paired with constraints?

          rhs_split <- strsplit(rhs, "\\s+")[[1]] # lav_mplus_trim(
          if (length(lhs_split) != length(rhs_split)) {
            browser()
            lav_msg_stop(gettext(
              "PWITH command does not have the same number of arguments on
              the left and right sides."))
          }

          cmd <- sapply(seq_along(lhs_split), function(i) {
            paste(lhs_split[i], lav_operator, rhs_split[i])
          })
        } else {
          # insert plus signs on the rhs as long as it isn't preceded or
          #    followed by a plus already
          rhs <- gsub("(?<!\\+)\\s+(?!\\+)", " + ", rhs, perl = TRUE)

          if (length(lhs_split) > 1L) {
            # expand using possible combinations
            cmd <- sapply(lhs_split, function(larg) {
              pair <- paste(larg, lav_operator, rhs)
              pair
            })
          } else {
            cmd <- paste(lhs, lav_operator, rhs)
          }
        }
          } else if ((means_scales <-
        regexpr("^\\s*([\\[\\{])([^\\]\\}]+)[\\]\\}]\\s*$", cmd,
         ignore.case = TRUE, perl = TRUE))[1L] > 0) {
        # intercepts/means or scales
        # first capture is the operator: [ or {
        operator <- substr(cmd, attr(means_scales, "capture.start")[1L],
        attr(means_scales, "capture.start")[1L] + attr(means_scales,
           "capture.length")[1L] - 1)

        params <- substr(cmd, attr(means_scales, "capture.start")[2L],
         attr(means_scales, "capture.start")[2L] +
           attr(means_scales, "capture.length")[2L] - 1)

        # obtain parameters with no modifiers specified for LHS
        params_no_modifiers <- sub(
          "^\\s*[\\[\\{]([^\\]\\}]+)[\\]\\}]\\s*$", "\\1",
           attr(cmd, "noModifiers"), perl = TRUE)

        means_scales_split <- strsplit(params, "\\s+")[[1]] # lav_mplus_trim(
        means_scales_no_modifiers_split <-                                 # nolint
          strsplit(params_no_modifiers, "\\s+")[[1]] # lav_mplus_trim(

        if (operator == "[") {
          # Tricky syntax shift (and corresponding kludge). For means, need to
          #    put constraint on RHS as pre-multiplier of 1 (e.g., x1 ~ 5*1).
          # But lav_mplus_cmd_constraints returns constraints
          #     multiplied by parameters
          cmd <- sapply(means_scales_split, function(v) {
            # shift pre-multiplier
            if ((premult <- regexpr("([^\\*]+\\*[^\\*]+)\\*([^\\*]+)",
             v, perl = TRUE))[1L] > 0) { # double modifier: label and constraint
              modifier <- substr(v, attr(premult, "capture.start")[1L],
               attr(premult, "capture.start")[1L] +
                 attr(premult, "capture.length")[1L] - 1)
              param_name <- substr(v, attr(premult, "capture.start")[2L],
               attr(premult, "capture.start")[2L] +
                 attr(premult, "capture.length")[2L] - 1)
              paste0(param_name, " ~ ", modifier, "*1")
            } else if ((premult <-
              regexpr("([^\\*]+)\\*([^\\*]+)", v, perl = TRUE))[1L] > 0) {
              modifier <- substr(v, attr(premult, "capture.start")[1L],
                attr(premult, "capture.start")[1L] +
                  attr(premult, "capture.length")[1L] - 1)
              param_name <- substr(v, attr(premult, "capture.start")[2L],
                attr(premult, "capture.start")[2L] +
                   attr(premult, "capture.length")[2L] - 1)
              paste0(param_name, " ~ ", modifier, "*1")
            } else {
              paste(v, "~ 1")
            }
          })
        } else if (operator == "{") {
          # only include constraints on RHS
          cmd <- sapply(seq_along(means_scales_split),
                 function(v) {
                   paste(means_scales_no_modifiers_split[v],
                   "~*~", means_scales_split[v])
                 })
        } else {
          lav_msg_stop(gettext("What's the operator?!"))
        }
                 } else if (grepl("|", cmd, fixed = TRUE)) {
        # expand growth modeling language
        cmd <- lav_mplus_cmd_growth(cmd)
      } else { # no operator, no means, must be variance.
        # cat("assuming vars: ", cmd, "\n")

        vars_lhs <- strsplit(attr(cmd, "noModifiers"), "\\s+")[[1]]
        # lav_mplus_trim(
        vars_rhs <- strsplit(cmd, "\\s+")[[1]] # lav_mplus_trim(

        cmd <- sapply(seq_along(vars_lhs),
             function(v) paste(vars_lhs[v], "~~", vars_rhs[v]))
      }

      # handle threshold substitution: $ -> |
      cmd <- gsub("$", "|", cmd, fixed = TRUE)

      # if we have both starting/fixed values and constraints, these must be
      #    handled by separate commands.
      # starting and fixed values are already handled in the pipeline by this
      #    point, so should be evident in the command
      # bfi BY lab1*start(1)*bfi_1 ==> bfi BY lab1*bfi_1 + start(1)*bfi_1
      double_asterisks <- grepl(
        "\\s*[\\w\\(\\)\\.]+\\*[\\w\\(\\)\\.]+\\*[\\w\\(\\)\\.]+",
         cmd, perl = TRUE)

      if (isTRUE(double_asterisks[1])) {
        ss <- strsplit(cmd, "*", fixed = TRUE)[[1]]
        if (length(ss) != 3) {
          lav_msg_warn(gettext("problem interpreting double asterisk syntax:"),
                       cmd) # sanity check on my logic
        } else {
          cmd <- paste0(ss[1], "*", ss[3], " + ", ss[2], "*", ss[3])
        }
      }

      lavaan_out <- c(lavaan_out, cmd)
    }
  }

  # new threshold syntax shifts things to the form:
  # VAR | t1 + t2 + t3 (left to write ordering)
  # Parameter labels, fixed values, and starting values are tacked on in
  # the usual way, like
  # VAR | 5*t1 + start(1.5)*t2 + par_label*t3 (left to write ordering)

  thresh_lines <- grep("^\\s*[A-z]+\\w*\\|\\d+", lavaan_out, perl = TRUE)
  if (length(thresh_lines) > 0L) {
    thresh_vars <- unname(sub("^\\s*([A-z]+\\w*).*", "\\1",
     lavaan_out[thresh_lines], perl = TRUE))
    thresh_split <- split(thresh_lines, thresh_vars)
    drop_elements <- c()
    for (i in seq_along(thresh_split)) {
      this_set <- lavaan_out[thresh_split[[i]]]
      tnum <- as.integer(sub("^\\s*[A-z]+\\w*\\|(\\d+)\\s*.*", "\\1", this_set))
      this_set <- this_set[order(tnum)]
      # ensure that threshold numbering matches ascending order
      this_set <- sub("[^~]+\\s*~\\s*", "",
      this_set, perl = TRUE) # drop variable and ~

      # convert to new t1, t2 syntax by combining modifiers
      #    with threshold numbers
      this_set <- sapply(seq_along(this_set), function(j) {
        # gsub("[^~]+\\s*~\\s*([\\w\\.\\-]+\\*)*1",
        #   paste0("\\1t", j), this_set[j], perl=TRUE)
        gsub("([\\w\\.\\-]+\\*)*1", paste0("\\1t", j), this_set[j], perl = TRUE)
      })

      new_str <- paste(names(thresh_split)[i], "|",
       paste(this_set, collapse = " + "))
      # replace in model string on the first line having relevant syntax
      lavaan_out[thresh_split[[i]][1]] <- new_str
      drop_elements <- c(drop_elements, thresh_split[[i]][-1])
    }
    lavaan_out <- lavaan_out[-drop_elements]
  }


  # tack on constraint syntax, if included
  lavaan_out <- c(lavaan_out, con_syntax)

  # for now, include a final lav_mplus_trim call since some arguments have
  #   leading/trailing space stripped.
  wrap <- paste(lav_mplus_cmd_wrap(lavaan_out, width = 90, exdent = 5),
   collapse = "\n") # lav_mplus_trim(
  wrap
}

lav_mplus_lavaan <- function(inpfile, run = TRUE) {
  stopifnot(length(inpfile) == 1L)
  stopifnot(grepl("\\.inp$", inpfile, ignore.case = TRUE))
  if (!file.exists(inpfile)) {
    lav_msg_stop(gettext("Could not find file:"), inpfile)
  }

  # for future consideration. For now, require a .inp file
  #  if (length(inpfile) == 1L && grepl("\\.inp$", inpfile)) {
  #    if (!file.exists(inpfile)) { stop("Could not find file: ", inpfile) }
  #    inpfile.text <- scan(inpfile, what="character", sep="\n",
  #         strip.white=FALSE, blank.lines.skip=FALSE, quiet=TRUE)
  #  } else {
  #    #assume that inpfile itself is syntax (e.g., in a character vector)
  #    inpfile.text <- inpfile
  #  }

  inpfile_text <- scan(inpfile, what = "character", sep = "\n",
           strip.white = FALSE, blank.lines.skip = FALSE, quiet = TRUE)
  sections <- lav_mplus_text_sections(inpfile_text, inpfile)

  mplus_inp <- list()

  mplus_inp$title <- lav_mplus_trim(paste(sections$title, collapse = " "))
  mplus_inp$data <- lav_mplus_text_fields(sections$data, required = "file")
  mplus_inp$variable <- lav_mplus_text_fields(sections$variable,
       required = "names")
  mplus_inp$analysis <- lav_mplus_text_fields(sections$analysis)

  meanstructure <- "default" # lavaan default
  if (!is.null(mplus_inp$analysis$model)) {
    if (tolower(mplus_inp$analysis$model) == "nomeanstructure") {
      meanstructure <- FALSE
    } # explicitly disable mean structure
  }

  information <- "default" # lavaan default
  if (!is.null(mplus_inp$analysis$information)) {
    information <- tolower(mplus_inp$analysis$information)
  }

  estimator <- "default"
  if (!is.null(est <- mplus_inp$analysis$estimator)) {
    # no memory of what this is up to....
    if (toupper(est) == "MUML") lav_msg_warn(gettext(
      "Mplus does not support MUML estimator. Using default instead."))
    estimator <- est

    # march 2013: handle case where categorical data are specified,
    #    but ML-based estimator requested.
    # use WLSMV instead
    if (!is.null(mplus_inp$variable$categorical) &&
        toupper(substr(mplus_inp$analysis$estimator, 1, 2)) == "ML") {
      lav_msg_warn(gettext(
        "Lavaan does not yet support ML-based estimation for categorical data.
        Reverting to WLSMV"))
      estimator <- "WLSMV"
    }
  }

  # expand hyphens in variable names and split into vector
  #   that will be the names for read.table
  mplus_inp$variable$names <-
    strsplit(lav_mplus_expand_cmd(mplus_inp$variable$names), "\\s+",
             perl = TRUE)[[1]]

  # expand hyphens in categorical declaration
  if (!is.null(mplus_inp$variable$categorical))
    mplus_inp$variable$categorical <- strsplit(
      lav_mplus_expand_cmd(mplus_inp$variable$categorical), "\\s+",
      perl = TRUE)[[1]]

  # convert mplus syntax to lavaan syntax
  mplus_inp$model <- lav_mplus_syntax_model(sections$model)

  # handle model constraint
  if ("model.constraint" %in% names(sections)) {
    mplus_inp$model.constraint <-
      lav_mplus_syntax_constraints(sections$model.constraint)
    mplus_inp$model <-
      paste(mplus_inp$model, mplus_inp$model.constraint, sep = "\n")
  }

  # read mplus data (and handle missing spec)
  mplus_inp$data <- lav_mplus_path_data(mplus_inp, inpfile)

  # handle bootstrapping specification
  se <- "default"
  bootstrap <- 1000L
  test <- "default"
  if (!is.null(mplus_inp$analysis$bootstrap)) {
    boot_type <- "standard"
    # check whether standard versus residual bootstrap is specified
    if ((boot_match <- regexpr("\\((\\w+)\\)", mplus_inp$analysis$bootstrap,
     perl = TRUE)) > 0L) {
      boot_type <- tolower(substr(mplus_inp$analysis$bootstrap,
        attr(boot_match, "capture.start"), attr(boot_match, "capture.start") +
          attr(boot_match, "capture.length") - 1L))
    }

    if (boot_type == "residual") test <- "Bollen.Stine"

    se <- "bootstrap"

    if ((nboot_match <- regexpr("^\\s*(\\d+)",
    mplus_inp$analysis$bootstrap, perl = TRUE)) > 0L) {
      bootstrap <- as.numeric(substr(mplus_inp$analysis$bootstrap,
        attr(nboot_match, "capture.start"), attr(nboot_match, "capture.start") +
          attr(nboot_match, "capture.length") - 1L))
    }
  }

  if (run) {
    fit <- sem(mplus_inp$model, data = mplus_inp$data, meanstructure =
      meanstructure, mimic = "Mplus", estimator = estimator, test = test,
      se = se, bootstrap = bootstrap, information = information)
    fit@external <- list(mplus.inp = mplus_inp)
  } else {
    fit <- mplus_inp # just return the syntax outside of a lavaan object
  }

  fit
}


lav_mplus_text_fields <- function(section_text, required) {
  if (is.null(section_text)) {
    return(NULL)
  }

  # The parser breaks down when there is a line with a trailing comment because
  #    then splitting on semicolon will combine it with the following line
  # Thus, trim off trailing comments before initial split
  section_text <- gsub("\\s*!.*$", "", section_text, perl = TRUE)
  section_split <- strsplit(paste(section_text, collapse = " "), ";",
       fixed = TRUE)[[1]] # split on semicolons
  section_divide <- list()

  for (cmd in section_split) {
    if (grepl("^\\s*!.*", cmd, perl = TRUE)) next # skip comment lines
    if (grepl("^\\s+$", cmd, perl = TRUE)) next # skip blank lines

    # mplus is apparently tolerant of specifications that don't include IS/ARE/=
    # example: usevariables x1-x10;
    # thus, split on spaces and assume that first element is lhs,
    #    drop second element if IS/ARE/=, and assume remainder is rhs

    # but if user uses equals sign, then spaces will not always
    #           be present (e.g., usevariables=x1-x10)
    if ((leading_equals <- regexpr(
      "^\\s*[A-Za-z]+[A-Za-z_-]*\\s*(=)", cmd[1L], perl = TRUE))[1L] > 0) {
      cmd_name <- lav_mplus_trim(substr(cmd[1L], 1,
         attr(leading_equals, "capture.start") - 1))
      cmd_args <- lav_mplus_trim(substr(cmd[1L],
        attr(leading_equals, "capture.start") + 1, nchar(cmd[1L])))
    } else {
      cmd_spacesplit <-
        strsplit(lav_mplus_trim(cmd[1L]), "\\s+", perl = TRUE)[[1L]]

      if (length(cmd_spacesplit) < 2L) {
        # for future: make room for this function to prase things
        #     like just TECH13 (no rhs)
      } else {
        cmd_name <- lav_mplus_trim(cmd_spacesplit[1L])
        if (length(cmd_spacesplit) > 2L && tolower(cmd_spacesplit[2L]) %in%
                    c("is", "are")) {
          cmd_args <- paste(cmd_spacesplit[3L:length(cmd_spacesplit)],
               collapse = " ") # remainder, removing is/are
        } else {
          cmd_args <- paste(cmd_spacesplit[2L:length(cmd_spacesplit)],
               collapse = " ") # is/are not used, so just join rhs
        }
      }
    }

    section_divide[[make.names(tolower(cmd_name))]] <- cmd_args
  }

  if (!missing(required)) {
    stopifnot(all(required %in% names(section_divide)))
  }
  section_divide
}

# helper function
lav_mplus_path_splitted <- function(abspath) {
  # function to split path into path and filename
  # code adapted from R.utils filePath command
  if (!is.character(abspath)) lav_msg_stop(gettext(
    "Path not a character string"))
  if (nchar(abspath) < 1 || is.na(abspath)) lav_msg_stop(gettext(
    "Path is missing or of zero length"))

  components <- strsplit(abspath, split = "[\\/]")[[1]]
  lcom <- length(components)

  stopifnot(lcom > 0)

  # the file is the last element in the list.
  # In the case of length == 1, this will extract the only element.
  rel_filename <- components[lcom]
  absolute <- FALSE

  if (lcom == 1) {
    dirpart <- NA_character_
  } else if (lcom > 1) {
    # drop the file from the list (the last element)
    components <- components[-lcom]
    dirpart <- do.call("file.path", as.list(components))

    # if path begins with C:, /, //, or \\, then treat as absolute
    if (grepl("^([A-Z]{1}:|/|//|\\\\)+.*$", dirpart, perl = TRUE))
        absolute <- TRUE
  }

  list(directory = dirpart, filename = rel_filename, absolute = absolute)
}

lav_mplus_path_data <- function(mplus_inp, inpfile) {
  # handle issue of lav_mplus_lavaan being called with an absolute path,
  #   whereas mplus has only a local data file
  inpfile_split <- lav_mplus_path_splitted(inpfile)
  datfile_split <- lav_mplus_path_splitted(mplus_inp$data$file)

  # if inp file target directory is non-empty, but mplus data is without
  # directory, then append inp file directory to mplus data. This ensures that
  # R need not be in the working directory to read the dat file. But if mplus
  # data has an absolute directory, don't append

  # if mplus data directory is present and absolute, or if no directory
  #    in input file, just use filename as is
  if (!is.na(datfile_split$directory) && datfile_split$absolute) {
    dat_file <- mplus_inp$data$file
    # else: just use mplus data filename if it has absolute path
  } else if (is.na(inpfile_split$directory)) {
    dat_file <- mplus_inp$data$file
  } else {
    # just use mplus data filename if inp file is missing path (working dir)
    dat_file <- file.path(inpfile_split$directory, mplus_inp$data$file)
  } # dat file path is relative or absent, and inp file directory is present

  if (!file.exists(dat_file)) {
    lav_msg_warn(gettext("Cannot find data file:"), dat_file)
    return(NULL)
  }

  # handle missing is/are:
  miss_list <- NULL
  if (!is.null(miss_spec <- mplus_inp$variable$missing)) {
    expand_miss_vec <- function(miss_str) {
      # sub-function to obtain a vector of all missing values
      #   within a set of parentheses
      miss_split <- strsplit(miss_str, "\\s+")[[1L]]
      miss_vals <- c()
      for (f in miss_split) {
        if ((hyphen_pos <- regexpr("\\d+(-)\\d+", f, perl = TRUE))[1L] > -1L) {
          # expand hyphen
          pre_hyphen <- substr(f, 1, attr(hyphen_pos, "capture.start") - 1)
          post_hyphen <- substr(f,
              attr(hyphen_pos, "capture.start") + 1, nchar(f))
          miss_vals <- c(miss_vals, as.character(seq(pre_hyphen, post_hyphen)))
        } else {
          # append to vector
          miss_vals <- c(miss_vals, f)
        }
      }
      as.numeric(miss_vals)
    }

    if (miss_spec == "." || miss_spec == "*") { # case 1: MISSING ARE|=|IS .;
      na_strings <- miss_spec
    } else if ((all_match <- regexpr("\\s*ALL\\s*\\(([^\\)]+)\\)", miss_spec,
     perl = TRUE))[1L] > -1L) { # case 2: use of ALL with parens
      miss_str <- lav_mplus_trim(substr(miss_spec,
        attr(all_match, "capture.start"),
        attr(all_match, "capture.start") +
          attr(all_match, "capture.length") - 1L))
      na_strings <- expand_miss_vec(miss_str)
    } else { # case 3: specific missing values per variable
      # process each element
      miss_blocks <- gregexpr("(?:(\\w+)\\s+\\(([^\\)]+)\\))+", miss_spec,
         perl = TRUE)[[1]]
      miss_list <- list()

      if (miss_blocks[1L] > -1L) {
        for (i in seq_along(miss_blocks)) {
          vname <- substr(miss_spec, attr(miss_blocks,
            "capture.start")[i, 1L], attr(miss_blocks, "capture.start")[i, 1L] +
               attr(miss_blocks, "capture.length")[i, 1L] - 1L)
          vmiss <- substr(miss_spec, attr(miss_blocks, "capture.start")[i, 2L],
           attr(miss_blocks, "capture.start")[i, 2L] + attr(miss_blocks,
              "capture.length")[i, 2L] - 1L)

          vname_hyphen <- regexpr("(\\w+)-(\\w+)", vname, perl = TRUE)[1L]
          if (vname_hyphen > -1L) {
            # lookup against variable names
            vstart <- which(mplus_inp$variable$names == substr(vname,
               attr(vname_hyphen, "capture.start")[1L], attr(vname_hyphen,
                                   "capture.start")[1L]
                 + attr(vname_hyphen, "capture.length")[1L] - 1L))
            vend <- which(mplus_inp$variable$names ==
              substr(vname, attr(vname_hyphen, "capture.start")[2L],
            attr(vname_hyphen, "capture.start")[2L] +
              attr(vname_hyphen, "capture.length")[2L] - 1L))
            if (length(vstart) == 0L || length(vend) == 0L) {
              lav_msg_stop(gettext("Unable to lookup missing variable list: "),
                           vname)
            }
            # I suppose start or finish could be mixed up
            if (vstart > vend) {
              vstart_orig <- vstart
              vstart <- vend
              vend <- vstart_orig
            }
            vname <- mplus_inp$variable$names[vstart:vend]
          }

          miss_vals <- expand_miss_vec(vmiss)

          for (j in seq_along(vname)) {
            miss_list[[vname[j]]] <- miss_vals
          }
        }
      } else {
        lav_msg_stop(gettext("I don't understand this missing specification:"),
                     miss_spec)
      }
    }
  } else {
    na_strings <- "NA"
  }

  if (!is.null(miss_list)) {
    dat <- read.table(dat_file, header = FALSE,
      col.names = mplus_inp$variable$names, colClasses = "numeric")
    # loop over variables in missList and set missing values to NA
    dat[, names(miss_list)] <- lapply(names(miss_list), function(vmiss) {
      dat[which(dat[, vmiss] %in% miss_list[[vmiss]]), vmiss] <- NA
      dat[, vmiss]
    })

    names(dat) <- mplus_inp$variable$names # loses these from the lapply
  } else {
    dat <- read.table(dat_file, header = FALSE,
      col.names = mplus_inp$variable$names, na.strings = na_strings,
      colClasses = "numeric")
  }


  # TODO: support covariance/mean+cov inputs

  # store categorical variables as ordered factors
  if (!is.null(mplus_inp$variable$categorical)) {
    dat[, c(mplus_inp$variable$categorical)] <-
      lapply(dat[, c(mplus_inp$variable$categorical), drop = FALSE], ordered)
  }

  dat
}


lav_mplus_text_sections <- function(inpfile_text, filename) {
  input_headers <- grep(
    "^\\s*(title:|data.*:|variable:|define:|analysis:|model.*:|output:|savedata:|plot:|montecarlo:)", # nolint
     inpfile_text, ignore.case = TRUE, perl = TRUE)

  stopifnot(length(input_headers) > 0L)

  mplus_sections <- list()

  for (h in seq_along(input_headers)) {
    section_end <- ifelse(h < length(input_headers), input_headers[h + 1] - 1,
        length(inpfile_text))
    section <- inpfile_text[input_headers[h]:section_end]
    section_name <- lav_mplus_trim(sub("^([^:]+):.*$", "\\1", section[1L],
     perl = TRUE)) # obtain text before the colon

    # dump section name from input syntax
    section[1L] <- sub("^[^:]+:(.*)$", "\\1", section[1L], perl = TRUE)

    mplus_sections[[make.names(tolower(section_name))]] <- section
  }

  mplus_sections
}
