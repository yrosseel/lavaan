# this code is written by Michael Hallquist
# First draft of parser to convert Mplus model syntax to lavaan model syntax

# idea: build parTable and run model from mplus syntax
# then perhaps write export function: parTable2Mplus
# and/or parTable2lavaan

trimSpace <- function(string) {
  stringTrim <- sapply(string, function(x) {
    x <- sub("^\\s*", "", x, perl = TRUE)
    x <- sub("\\s*$", "", x, perl = TRUE)
    return(x)
  }, USE.NAMES = FALSE)
  return(stringTrim)
}

# small utility function to join strings in a regexp loop
joinRegexExpand <- function(cmd, argExpand, matches, iterator, matchLength = "match.length") {
  if (iterator == 1 && matches[iterator] > 1) {
    pre <- substr(cmd, 1, matches[iterator] - 1)
  } else {
    pre <- ""
  }

  # if this is not the final match, then get sub-string between the end of this match and the beginning of the next
  # otherwise, match to the end of the command
  post.end <- ifelse(iterator < length(matches), matches[iterator + 1] - 1, nchar(cmd))
  post <- substr(cmd, matches[iterator] + attr(matches, matchLength)[iterator], post.end)

  cmd.expand <- paste(pre, argExpand, post, sep = "")
  return(cmd.expand)
}

# expand Mplus hyphen syntax (will also expand constraints with hyphens)
expandCmd <- function(cmd, alphaStart = TRUE) {
  # use negative lookahead and negative lookbehind to eliminate possibility of hyphen being used as a negative starting value (e.g., x*-1)
  # also avoid match of anything that includes a decimal point, such as a floating-point starting value -10.5*x1

  # if alphaStart==TRUE, then require that the matches before and after hyphens begin with alpha character
  # this is used for variable names, whereas the more generic expansion works for numeric constraints and such

  # need to do a better job of this so that u1-u20* is supported... I don't think the regexp below is general enough

  # if (alphaStart) {
  #  hyphens <- gregexpr("[_A-Za-z]+\\w*\\s*-\\s*[_A-Za-z]+\\w*", cmd, perl=TRUE)[[1]]
  # } else {
  #  hyphens <- gregexpr("(?!<(\\*|\\.))\\w+(?!(\\*|\\.))\\s*-\\s*(?!<(\\*|\\.))\\w+(?!(\\*|\\.))", cmd, perl=TRUE)[[1]]
  # }

  # hyphens <- gregexpr("(?!<(\\*|\\.))\\w+(?!(\\*|\\.))\\s*-\\s*(?!<(\\*|\\.))\\w+(?!(\\*|\\.))", cmd, perl=TRUE)[[1]]

  # support trailing @XXX. Still still fail on Trait1-Trait3*XXX
  hyphens <- gregexpr("(?!<(\\*|\\.))\\w+(?!(\\*|\\.))\\s*-\\s*(?!<(\\*|\\.))\\w+(?!(\\*|\\.))(@[\\d\\.\\-]+)?", cmd, perl = TRUE)[[1]]

  # Promising, but this is still failing in the case of x3*1 -4.25*x4
  # On either side of a hyphen, require alpha character followed by alphanumeric
  # This enforces that neither side of the hyphen can be a number
  # Alternatively, match digits on either side alone
  # hyphens <- gregexpr("([A-z]+\\w*\\s*-\\s*[A-z]+\\w*(@[\\d\\.-]+)?|\\d+\\s*-\\s*\\d+)", cmd, perl=TRUE)[[1]]

  if (hyphens[1L] > 0) {
    cmd.expand <- c()
    ep <- 1

    for (v in 1:length(hyphens)) {
      # match one keyword before and after hyphen
      argsplit <- strsplit(substr(cmd, hyphens[v], hyphens[v] + attr(hyphens, "match.length")[v] - 1), "\\s*-\\s*", perl = TRUE)[[1]]

      v_pre <- argsplit[1]
      v_post <- argsplit[2]

      v_post.suffix <- sub("^([^@]+)(@[\\d\\-.]+)?$", "\\2", v_post, perl = TRUE) # will be empty string if not present
      v_post <- sub("@[\\d\\-.]+$", "", v_post, perl = TRUE) # trim @ suffix

      # If v_pre and v_post contain leading alpha characters, verify that these prefixes match.
      # Otherwise, there is nothing to expand, as in the case of MODEL CONSTRAINT: e1e2=e1-e2_n.
      v_pre.alpha <- sub("\\d+$", "", v_pre, perl = TRUE)
      v_post.alpha <- sub("\\d+$", "", v_post, perl = TRUE)

      # only enforce prefix match if we have leading alpha characters (i.e., not simple numeric 1 - 3 syntax)
      if (length(v_pre.alpha) > 0L && length(v_post.alpha) > 0L) {
        # if alpha prefixes do match, assume that the hyphen is not for expansion (e.g., in subtraction case)
        if (v_pre.alpha != v_post.alpha) {
          return(cmd)
        }
      }

      # the basic positive lookbehind blows up with pure numeric constraints (1 - 3) because no alpha char precedes digit
      # can use an non-capturing alternation grouping to allow for digits only or the final digits after alphas (as in v_post.num)
      v_pre.num <- as.integer(sub("\\w*(?<=[A-Za-z_])(\\d+)$", "\\1", v_pre, perl = TRUE)) # use positive lookbehind to avoid greedy \w+ match -- capture all digits

      v_post.match <- regexpr("^(?:\\w*(?<=[A-Za-z_])(\\d+)|(\\d+))$", v_post, perl = TRUE)
      stopifnot(v_post.match[1L] > 0)

      # match mat be under capture[1] or capture[2] because of alternation above
      whichCapture <- which(attr(v_post.match, "capture.start") > 0)

      v_post.num <- as.integer(substr(v_post, attr(v_post.match, "capture.start")[whichCapture], attr(v_post.match, "capture.start")[whichCapture] + attr(v_post.match, "capture.length")[whichCapture] - 1))
      v_post.prefix <- substr(v_post, 1, attr(v_post.match, "capture.start")[whichCapture] - 1) # just trusting that pre and post match

      if (is.na(v_pre.num) || is.na(v_post.num)) lav_msg_stop(
        gettext("Cannot expand variables:"), v_pre, ", ", v_post)
      v_expand <- paste(v_post.prefix, v_pre.num:v_post.num, v_post.suffix, sep = "", collapse = " ")

      # for first hyphen, there may be non-hyphenated syntax preceding the initial match
      cmd.expand[ep] <- joinRegexExpand(cmd, v_expand, hyphens, v)

      # This won't really work because the cmd.expand element may contain other variables
      # that are at the beginning or end, prior to hyphen stuff
      # This is superseded by logic above where @ is included in hyphen match, then trapped as suffix
      # I don't think it will work yet for this Mplus syntax: y1-y10*5 -- the 5 wouldn't propagate
      # handle the case of @ fixed values or * starting values used in a list
      # example: Trait1-Trait3@1
      ## if (grepl("@|\\*", cmd.expand[ep], perl=TRUE)) {
      ##   exp_split <- strsplit(cmd.expand[ep], "\\s+", perl=TRUE)[[1]]
      ##   suffixes <- sub("^([^@\\*]+)([@*][\\d\\.-]+)?$", "\\2", exp_split, perl=TRUE)
      ##   variables <- sub("^([^@\\*]+)([@*][\\d\\.-]+)?$", "\\1", exp_split, perl=TRUE)
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
    return(paste(cmd.expand, collapse = ""))
  } else {
    return(cmd) # no hyphens to expand
  }
}


# handle starting values and fixed parameters on rhs
parseFixStart <- function(cmd) {
  cmd.parse <- c()
  ep <- 1L

  # support ESEM-like syntax: F BY a1* a2*
  # The easy path: putting in 1s before we proceed on parsing
  # Mar2023 bugfix: support parenthesis after * in case a parameter constraint comes next
  cmd <- gsub("([A-z]+\\w*)\\s*\\*(?=\\s+\\(?[A-z]+|\\s*$)", "\\1*1", cmd, perl = TRUE)

  if ((fixed.starts <- gregexpr("[\\w\\.\\-$]+\\s*([@*])\\s*[\\w\\.\\-]+", cmd, perl = TRUE)[[1]])[1L] > 0) { # shouldn't it be \\*, not * ?! Come back to this.
    for (f in 1:length(fixed.starts)) {
      # capture above obtains the fixed/start character (@ or *), whereas match obtains the full regex match
      opchar <- substr(cmd, attr(fixed.starts, "capture.start")[f], attr(fixed.starts, "capture.start")[f] + attr(fixed.starts, "capture.length")[f] - 1)

      # match arguments around asterisk/at symbol
      argsplit <- strsplit(substr(cmd, fixed.starts[f], fixed.starts[f] + attr(fixed.starts, "match.length")[f] - 1), paste0("\\s*", ifelse(opchar == "*", "\\*", opchar), "\\s*"), perl = TRUE)[[1]]
      v_pre <- argsplit[1]
      v_post <- argsplit[2]

      if (suppressWarnings(is.na(as.numeric(v_pre)))) { # fixed.starts value post-multiplier
        var <- v_pre
        val <- v_post
      } else if (suppressWarnings(is.na(as.numeric(v_post)))) { # starting value pre-multiplier
        var <- v_post
        val <- v_pre
      } else {
        lav_msg_stop(
          gettext("Cannot parse Mplus fixed/starts values specification:"),
          v_pre, v_post)
      }

      if (opchar == "@") {
        cmd.parse[ep] <- joinRegexExpand(cmd, paste0(val, "*", var, sep = ""), fixed.starts, f)
        ep <- ep + 1L
      } else {
        cmd.parse[ep] <- joinRegexExpand(cmd, paste0("start(", val, ")*", var, sep = ""), fixed.starts, f)
        ep <- ep + 1L
      }
    }
    return(paste(cmd.parse, collapse = ""))
  } else {
    return(cmd)
  }
}

parseConstraints <- function(cmd) {
  # Allow cmd to have newlines embedded. In this case, split on newlines, and loop over and parse each chunk
  # Dump leading and trailing newlines, which contain no information about constraints, but may add dummy elements to vector after strsplit
  # Maybe return LHS and RHS parsed command where constraints only appear on the RHS, whereas the LHS contains only parameters.
  # Example: LHS is v1 v2 v3 and RHS is con1*v1 con2*v2 con3*v3

  cmd.split <- strsplit(cmd, "\n")[[1]]

  # drop empty lines (especially leading newline)
  cmd.split <- if (length(emptyPos <- which(cmd.split == "")) > 0L) {
    cmd.split[-1 * emptyPos]
  } else {
    cmd.split
  }

  # Create a version of the command with no modifiers (constraints, starting values, etc.) specifications.
  # This is useful for syntax that uses the params on the LHS and with a modified RHS. Example: v1 ~~ conB*v1
  cmd.nomodifiers <- paste0(gsub("(start\\([^\\)]+\\)\\*|[\\d\\-\\.]+\\*)", "", cmd.split, perl = TRUE), collapse = " ") # peel off premultiplication
  cmd.nomodifiers <- gsub("\\([^\\)]+\\)", "", cmd.nomodifiers, perl = TRUE)

  cmd.tojoin <- c() # will store all chunks divided by newlines, which will be joined at the end.

  # iterate over each newline segment
  for (n in 1:length(cmd.split)) {
    # in principle, now that we respect newlines, parens should only be of length 1, since Mplus syntax dictates newlines for each use of parentheses for constraints
    if ((parens <- gregexpr("(?<!start)\\(([^\\)]+)\\)", cmd.split[n], perl = TRUE)[[1L]])[1L] > 0) { # match parentheses, but not start()
      # the syntax chunk after all parentheses have been matched
      cmd.expand <- c()

      for (p in 1:length(parens)) {
        # string within the constraint parentheses
        constraints <- substr(cmd.split[n], attr(parens, "capture.start")[p], attr(parens, "capture.start")[p] + attr(parens, "capture.length")[p] - 1)

        # Divide constraints on spaces to determine number of constraints to parse. Use trimSpace to avoid problem of user including leading/trailing spaces within parentheses.
        con.split <- strsplit(trimSpace(constraints), "\\s+", perl = TRUE)[[1]]

        # if Mplus uses a purely numeric constraint, then add ".con" prefix to be consistent with R naming.
        con.split <- sapply(con.split, function(x) {
          if (!suppressWarnings(is.na(as.numeric(x)))) {
            make.names(paste0(".con", x))
          } else {
            x
          }
        })

        # determine the parameters that precede the parentheses (either first character for p == 1 or character after preceding parentheses)
        prestrStart <- ifelse(p > 1, attr(parens, "capture.start")[p - 1] + attr(parens, "capture.length")[p - 1] + 1, 1)

        # obtain the parameters that precede the parentheses, divide into arguments on spaces
        # use trimSpace here because first char after prestrStart for p > 1 will probably be a space
        precmd.split <- strsplit(trimSpace(substr(cmd.split[n], prestrStart, parens[p] - 1)), "\\s+", perl = TRUE)[[1]]

        # peel off any potential LHS arguments, such as F1 BY
        precmdLHSOp <- which(tolower(precmd.split) %in% c("by", "with", "on"))
        if (any(precmdLHSOp)) {
          lhsop <- paste0(precmd.split[1:precmdLHSOp[1L]], " ", collapse = " ") # join lhs and op as a single string, add trailing space so that paste with expanded RHS is right.
          rhs <- precmd.split[(precmdLHSOp + 1):length(precmd.split)]
        } else {
          lhsop <- ""
          rhs <- precmd.split
        }

        if (length(con.split) > 1L) {
          # several constraints listed within parentheses. Example: F1 BY X1 X2 X3 X4 (C2 C3 C4)
          # thus, backwards match the constraints to parameters

          # restrict parameters to backwards match to be of the same length as number of constraints
          rhs.backmatch <- rhs[(length(rhs) - length(con.split) + 1):length(rhs)]

          rhs.expand <- c()

          # check that no mean or scale markers are part of the rhs param to expand
          if ((preMark.match <- regexpr("^\\s*[\\[\\{]", rhs.backmatch[1L], perl = TRUE))[1L] > 0) {
            preMark <- substr(rhs.backmatch[1L], preMark.match[1L], preMark.match[1L] + attr(preMark.match, "match.length")[1L] - 1)
            rhs.backmatch[1L] <- substr(rhs.backmatch[1L], preMark.match[1L] + attr(preMark.match, "match.length")[1L], nchar(rhs.backmatch[1L]))
          } else {
            preMark <- ""
          }

          if ((postMark.match <- regexpr("[\\]\\}]\\s*$", rhs.backmatch[length(rhs.backmatch)], perl = TRUE))[1L] > 0) {
            postMark <- substr(rhs.backmatch[length(rhs.backmatch)], postMark.match[1L], nchar(rhs.backmatch[length(rhs.backmatch)]))
            rhs.backmatch[length(rhs.backmatch)] <- substr(rhs.backmatch[length(rhs.backmatch)], 1, postMark.match[1L] - 1)
          } else {
            postMark <- ""
          }


          # pre-multiply each parameter with each corresponding constraint
          for (i in 1:length(rhs.backmatch)) {
            rhs.expand[i] <- paste0(con.split[i], "*", rhs.backmatch[i])
          }

          # join rhs as string and add back in mean/scale operator, if present
          rhs.expand <- paste0(preMark, paste(rhs.expand, collapse = " "), postMark)

          # if there were params that preceded the backwards match, then add these back to the syntax
          # append this syntax to the parsed command, cmd.expand
          if (length(rhs) - length(con.split) > 0L) {
            cmd.expand <- c(cmd.expand, paste(lhsop, paste(rhs[1:(length(rhs) - length(con.split))], collapse = " "), rhs.expand))
          } else {
            cmd.expand <- c(cmd.expand, paste0(lhsop, rhs.expand))
          }
        } else {
          # should be able to reduce redundancy with above

          # all parameters on the right hand side are to be equated
          # thus, pre-multiply each parameter by the constraint

          # check that no mean or scale markers are part of the rhs param to expand
          # DUPE CODE FROM ABOVE. Make Function?!
          if ((preMark.match <- regexpr("^\\s*[\\[\\{]", rhs[1L], perl = TRUE))[1L] > 0) {
            preMark <- substr(rhs[1L], preMark.match[1L], preMark.match[1L] + attr(preMark.match, "match.length")[1L] - 1)
            rhs[1L] <- substr(rhs[1L], preMark.match[1L] + attr(preMark.match, "match.length")[1L], nchar(rhs[1L]))
          } else {
            preMark <- ""
          }

          if ((postMark.match <- regexpr("[\\]\\}]\\s*$", rhs[length(rhs)], perl = TRUE))[1L] > 0) {
            postMark <- substr(rhs[length(rhs)], postMark.match[1L], nchar(rhs[length(rhs)]))
            rhs[length(rhs)] <- substr(rhs[length(rhs)], 1, postMark.match[1L] - 1)
          } else {
            postMark <- ""
          }


          rhs.expand <- c()
          for (i in 1:length(rhs)) {
            rhs.expand[i] <- paste0(con.split[1L], "*", rhs[i])
          }

          # join rhs as string
          rhs.expand <- paste0(preMark, paste(rhs.expand, collapse = " "), postMark)

          cmd.expand <- c(cmd.expand, paste0(lhsop, rhs.expand))
        }
      }

      cmd.tojoin[n] <- paste(cmd.expand, collapse = " ")
    } else {
      cmd.tojoin[n] <- cmd.split[n]
    } # no parens
  }

  # eliminate newlines in this function so that they don't mess up \\s+ splits downstream
  toReturn <- paste(cmd.tojoin, collapse = " ")
  attr(toReturn, "noModifiers") <- cmd.nomodifiers

  return(toReturn)
}

expandGrowthCmd <- function(cmd) {
  # can assume that any spaces between tscore and variable were stripped by parseFixStart

  # verify that this is not a random slope
  if (any(tolower(strsplit(cmd, "\\s+", perl = TRUE)[[1]]) %in% c("on", "at"))) {
    lav_msg_stop(gettext(
      "lavaan does not support random slopes or individually varying
      growth model time scores"))
  }

  cmd.split <- strsplit(cmd, "\\s*\\|\\s*", perl = TRUE)[[1]]
  if (!length(cmd.split) == 2) {
    lav_msg_stop(gettext("Unknown growth syntax:"), cmd)
  }

  lhs <- cmd.split[1]
  lhs.split <- strsplit(lhs, "\\s+", perl = TRUE)[[1]]

  rhs <- cmd.split[2]
  rhs.split <- strsplit(rhs, "(\\*|\\s+)", perl = TRUE)[[1]]

  if (length(rhs.split) %% 2 != 0) {
    lav_msg_stop(gettext(
      "Number of variables and number of tscores does not match:"), rhs)
  }
  tscores <- as.numeric(rhs.split[1:length(rhs.split) %% 2 != 0]) # pre-multipliers

  vars <- rhs.split[1:length(rhs.split) %% 2 == 0]

  cmd.expand <- c()

  for (p in 0:(length(lhs.split) - 1)) {
    if (p == 0) {
      # intercept
      cmd.expand <- c(cmd.expand, paste(lhs.split[(p + 1)], "=~", paste("1*", vars, sep = "", collapse = " + ")))
    } else {
      cmd.expand <- c(cmd.expand, paste(lhs.split[(p + 1)], "=~", paste(tscores^p, "*", vars, sep = "", collapse = " + ")))
    }
  }

  return(cmd.expand)
}

# function to wrap long lines at a certain width, splitting on + symbols to be consistent with R syntax
wrapAfterPlus <- function(cmd, width = 90, exdent = 5) {
  result <- lapply(cmd, function(line) {
    if (nchar(line) > width) {
      split <- c()
      spos <- 1L

      plusMatch <- gregexpr("+", line, fixed = TRUE)[[1]]
      mpos <- 1L

      if (plusMatch[1L] > 0L) {
        # split after plus symbol
        charsRemain <- nchar(line)
        while (charsRemain > 0L) {
          toProcess <- substr(line, nchar(line) - charsRemain + 1, nchar(line))
          offset <- nchar(line) - charsRemain + 1

          if (nchar(remainder <- substr(line, offset, nchar(line))) <= (width - exdent)) {
            # remainder of line fits within width -- no need to continue wrapping
            split[spos] <- remainder
            charsRemain <- 0
          } else {
            wrapAt <- which(plusMatch < (width + offset - exdent))
            wrapAt <- wrapAt[length(wrapAt)] # at the final +

            split[spos] <- substr(line, offset, plusMatch[wrapAt])
            charsRemain <- charsRemain - nchar(split[spos])
            spos <- spos + 1
          }
        }

        # remove leading and trailing chars
        split <- trimSpace(split)

        # handle exdent
        split <- sapply(1:length(split), function(x) {
          if (x > 1) {
            paste0(paste(rep(" ", exdent), collapse = ""), split[x])
          } else {
            split[x]
          }
        })

        return(split)
      } else {
        return(strwrap(line, width = width, exdent = exdent)) # convention strwrap when no + present
      }
    } else {
      return(line)
    }
  })

  # bind together multi-line expansions into single vector
  return(unname(do.call(c, result)))
}

mplus2lavaan.constraintSyntax <- function(syntax) {
  # should probably pass in model syntax along with some tracking of which parameter labels are defined.

  # convert MODEL CONSTRAINT section to lavaan model syntax
  syntax <- paste(lapply(trimSpace(strsplit(syntax, "\n")), function(x) {
    if (length(x) == 0L && is.character(x)) "" else x
  }), collapse = "\n")

  # replace ! with # for comment lines. Also strip newline and replace with semicolon
  syntax <- gsub("(\\s*)!(.+)\n", "\\1#\\2;", syntax, perl = TRUE)

  # split into vector of strings
  # need to peel off leading or trailing newlines -- leads to parsing confusion downstream otherwise
  syntax.split <- gsub("(^\n|\n$)", "", unlist(strsplit(syntax, ";")), perl = TRUE)

  constraint.out <- c()

  # TODO: Handle PLOT and LOOP syntax for model constraints.
  # TODO: Handle DO loop convention

  # first parse new parameters defined in MODEL CONSTRAINT into a vector
  new.parameters <- c() # parameters that are defined in constraint section
  if (length(new.con.lines <- grep("^\\s*NEW\\s*\\([^\\)]+\\)", syntax.split, perl = TRUE, ignore.case = TRUE)) > 0L) {
    for (cmd in syntax.split[new.con.lines]) {
      # process new constraint definition
      new.con <- regexpr("^\\s*NEW\\s*\\(([^\\)]+)\\)", cmd, perl = TRUE, ignore.case = TRUE)
      if (new.con[1L] == -1)
        lav_msg_stop(gettext("Unable to parse names of new contraints"))
      new.con <- substr(cmd, attr(new.con, "capture.start"), attr(new.con, "capture.start") + attr(new.con, "capture.length") - 1L)
      new.con <- expandCmd(new.con) # allow for hyphen expansion
      new.parameters <- c(new.parameters, strsplit(trimSpace(new.con), "\\s+", perl = TRUE)[[1L]])
    }

    syntax.split <- syntax.split[-1L * new.con.lines] # drop out these lines
    parameters.undefined <- new.parameters # to be used below to handle ambiguity of equation versus definition
  }

  for (cmd in syntax.split) {
    if (grepl("^\\s*#", cmd, perl = TRUE)) { # comment line
      constraint.out <- c(constraint.out, gsub("\n", "", cmd, fixed = TRUE)) # drop any newlines
    } else if (grepl("^\\s+$", cmd, perl = TRUE)) {
      # do nothing, just a space line
    } else {
      # constraint proper
      cmd <- gsub("**", "^", cmd, fixed = TRUE) # handle exponent

      # lower case the math operations supported by Mplus to be consistent with R
      # match all math operators, then lower case each and rejoin string
      maths <- gregexpr("(SQRT|LOG|LOG10|EXP|ABS|SIN|COS|TAN|ASIN|ACOS|ATAN)\\s*\\(", cmd, perl = TRUE)[[1L]]
      if (maths[1L] > 0) {
        maths.replace <- c()
        ep <- 1

        for (i in 1:length(maths)) {
          operator <- tolower(substr(cmd, attr(maths, "capture.start")[i], attr(maths, "capture.start")[i] + attr(maths, "capture.length")[i] - 1))
          maths.replace[ep] <- joinRegexExpand(cmd, operator, maths, i, matchLength = "capture.length") # only match operator, not opening (
          ep <- ep + 1
        }
        cmd <- paste(maths.replace, collapse = "")
      }

      # equating some lhs and rhs: could reflect definition of new parameter
      if ((equals <- regexpr("=", cmd, fixed = TRUE))[1L] > 0) {
        lhs <- trimSpace(substr(cmd, 1, equals - 1))
        rhs <- trimSpace(substr(cmd, equals + attr(equals, "match.length"), nchar(cmd)))

        # possibility of lhs or rhs containing the single variable to be equated
        if (regexpr("\\s+", lhs, perl = TRUE)[1L] > 0L) {
          def <- rhs
          body <- lhs
        } else if (regexpr("\\s+", rhs, perl = TRUE)[1L] > 0L) {
          def <- lhs
          body <- rhs
        } else {
          # warning("Can't figure out which side of constraint defines a parameter")
          # this would occur for simple rel5 = rel2 sort of syntax
          def <- lhs
          body <- rhs
        }

        # must decide whether this is a new parameter (:=) or equation of exising labels (==)
        # alternatively, could be zero, as in  0 = x + y
        # this is tricky, because mplus doesn't differentiate definition from equation
        # consequently, could confuse the issue as in ex5.20
        # NEW(rel2 rel5 stan3 stan6);
        # rel2 = lam2**2*vf1/(lam2**2*vf1 + ve2);
        # rel5 = lam5**2*vf2/(lam5**2*vf2 + ve5);
        # rel5 = rel2;

        # for now, only define a new constraint if it's not already defined
        # otherwise equate
        if (def %in% new.parameters && def %in% parameters.undefined) {
          constraint.out <- c(constraint.out, paste(def, ":=", body))
          parameters.undefined <- parameters.undefined[!parameters.undefined == def]
        } else {
          constraint.out <- c(constraint.out, paste(def, "==", body))
        }
      } else {
        # inequality constraints -- paste as is
        constraint.out <- c(constraint.out, cmd)
      }
    }
  }

  wrap <- paste(wrapAfterPlus(constraint.out, width = 90, exdent = 5), collapse = "\n")
  return(wrap)
}

mplus2lavaan.modelSyntax <- function(syntax) {
  if (is.character(syntax)) {
    if (length(syntax) > 1L) {
      syntax <- paste(syntax, collapse = "\n")
    } # concatenate into a long string separated by newlines
  } else {
    lav_msg_stop(gettext(
      "mplus2lavaan.modelSyntax accepts a single character string or
      character vector containing all model syntax"))
  }

  # because this is now exposed as a function in the package, handle the case of the user passing in full .inp file as text
  # we should only be interested in the MODEL and MODEL CONSTRAINT sections
  by_line <- strsplit(syntax, "\r?\n", perl = TRUE)[[1]]
  inputHeaders <- grep("^\\s*(title:|data.*:|variable:|define:|analysis:|model.*:|output:|savedata:|plot:|montecarlo:)", by_line, ignore.case = TRUE, perl = TRUE)
  con_syntax <- c()
  if (length(inputHeaders) > 0L) {
    # warning("mplus2lavaan.modelSyntax is intended to accept only the model section, not an entire .inp file. For the .inp file case, use mplus2lavaan")
    parsed_syntax <- divideInputIntoSections(by_line, "local")

    # handle model constraint
    if ("model.constraint" %in% names(parsed_syntax)) {
      con_syntax <- strsplit(mplus2lavaan.constraintSyntax(parsed_syntax$model.constraint), "\n")[[1]]
    }

    # just keep model syntax before continuing
    syntax <- parsed_syntax$model
  }

  # initial strip of leading/trailing whitespace, which can interfere with splitting on spaces
  # strsplit generates character(0) for empty strings, which causes problems in paste because paste actually includes it as a literal
  # example: paste(list(character(0), "asdf", character(0)), collapse=" ")
  # thus, use lapply to convert these to empty strings first
  syntax <- paste(lapply(trimSpace(strsplit(syntax, "\n")), function(x) {
    if (length(x) == 0L && is.character(x)) "" else x
  }), collapse = "\n")

  # replace ! with # for comment lines. Also strip newline and replace with semicolon
  syntax <- gsub("(\\s*)!(.+)\n*", "\\1#\\2;", syntax, perl = TRUE)

  # new direction: retain newlines in parsed syntax until after constraints have been parsed

  # delete newlines
  # syntax <- gsub("\n", "", syntax, fixed=TRUE)

  # replace semicolons with newlines prior to split (divide into commands)
  # syntax <- gsub(";", "\n", syntax, fixed=TRUE)

  # split into vector of strings
  # syntax.split <- unlist( strsplit(syntax, "\n") )
  syntax.split <- trimSpace(unlist(strsplit(syntax, ";")))

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
  lavaan.out <- c()

  for (cmd in syntax.split) {
    if (grepl("^\\s*#", cmd, perl = TRUE)) { # comment line
      lavaan.out <- c(lavaan.out, gsub("\n", "", cmd, fixed = TRUE)) # drop any newlines (otherwise done by parseConstraints)
    } else if (grepl("^\\s*$", cmd, perl = TRUE)) {
      # do nothing, just a space or blank line
    } else {
      # hyphen expansion
      cmd <- expandCmd(cmd)

      # parse fixed parameters and starting values
      cmd <- parseFixStart(cmd)

      # parse any constraints here (avoid weird logic below)
      cmd <- parseConstraints(cmd)

      if ((op <- regexpr("\\s+(by|on|with|pwith)\\s+", cmd, ignore.case = TRUE, perl = TRUE))[1L] > 0) { # regressions, factors, covariances

        lhs <- substr(cmd, 1, op - 1) # using op takes match.start which will omit spaces before operator
        rhs <- substr(cmd, op + attr(op, "match.length"), nchar(cmd))
        operator <- tolower(substr(cmd, attr(op, "capture.start"), attr(op, "capture.start") + attr(op, "capture.length") - 1))

        if (operator == "by") {
          lav.operator <- "=~"
        } else if (operator == "with" || operator == "pwith") {
          lav.operator <- "~~"
        } else if (operator == "on") {
          lav.operator <- "~"
        }

        # handle parameter combinations
        lhs.split <- strsplit(lhs, "\\s+")[[1]] # trimSpace(

        # handle pwith syntax
        if (operator == "pwith") {
          # TODO: Figure out if pwith can be paired with constraints?

          rhs.split <- strsplit(rhs, "\\s+")[[1]] # trimSpace(
          if (length(lhs.split) != length(rhs.split)) {
            browser()
            lav_msg_stop(gettext(
              "PWITH command does not have the same number of arguments on
              the left and right sides."))
          }

          cmd <- sapply(1:length(lhs.split), function(i) paste(lhs.split[i], lav.operator, rhs.split[i]))
        } else {
          # insert plus signs on the rhs as long as it isn't preceded or followed by a plus already
          rhs <- gsub("(?<!\\+)\\s+(?!\\+)", " + ", rhs, perl = TRUE)

          if (length(lhs.split) > 1L) {
            # expand using possible combinations
            cmd <- sapply(lhs.split, function(larg) {
              pair <- paste(larg, lav.operator, rhs)
              return(pair)
            })
          } else {
            cmd <- paste(lhs, lav.operator, rhs)
          }
        }
      } else if ((means.scales <- regexpr("^\\s*([\\[\\{])([^\\]\\}]+)[\\]\\}]\\s*$", cmd, ignore.case = TRUE, perl = TRUE))[1L] > 0) { # intercepts/means or scales
        # first capture is the operator: [ or {
        operator <- substr(cmd, attr(means.scales, "capture.start")[1L], attr(means.scales, "capture.start")[1L] + attr(means.scales, "capture.length")[1L] - 1)

        params <- substr(cmd, attr(means.scales, "capture.start")[2L], attr(means.scales, "capture.start")[2L] + attr(means.scales, "capture.length")[2L] - 1)

        # obtain parameters with no modifiers specified for LHS
        params.noModifiers <- sub("^\\s*[\\[\\{]([^\\]\\}]+)[\\]\\}]\\s*$", "\\1", attr(cmd, "noModifiers"), perl = TRUE)

        means.scales.split <- strsplit(params, "\\s+")[[1]] # trimSpace(
        means.scales.noModifiers.split <- strsplit(params.noModifiers, "\\s+")[[1]] # trimSpace(

        if (operator == "[") {
          # Tricky syntax shift (and corresponding kludge). For means, need to put constraint on RHS as pre-multiplier of 1 (e.g., x1 ~ 5*1).
          # But parseConstraints returns constraints multiplied by parameters
          cmd <- sapply(means.scales.split, function(v) {
            # shift pre-multiplier
            if ((premult <- regexpr("([^\\*]+\\*[^\\*]+)\\*([^\\*]+)", v, perl = TRUE))[1L] > 0) { # double modifier: label and constraint
              modifier <- substr(v, attr(premult, "capture.start")[1L], attr(premult, "capture.start")[1L] + attr(premult, "capture.length")[1L] - 1)
              paramName <- substr(v, attr(premult, "capture.start")[2L], attr(premult, "capture.start")[2L] + attr(premult, "capture.length")[2L] - 1)
              paste0(paramName, " ~ ", modifier, "*1")
            } else if ((premult <- regexpr("([^\\*]+)\\*([^\\*]+)", v, perl = TRUE))[1L] > 0) {
              modifier <- substr(v, attr(premult, "capture.start")[1L], attr(premult, "capture.start")[1L] + attr(premult, "capture.length")[1L] - 1)
              paramName <- substr(v, attr(premult, "capture.start")[2L], attr(premult, "capture.start")[2L] + attr(premult, "capture.length")[2L] - 1)
              paste0(paramName, " ~ ", modifier, "*1")
            } else {
              paste(v, "~ 1")
            }
          })
        } else if (operator == "{") {
          # only include constraints on RHS
          cmd <- sapply(1:length(means.scales.split), function(v) paste(means.scales.noModifiers.split[v], "~*~", means.scales.split[v]))
        } else {
          lav_msg_stop(gettext("What's the operator?!"))
        }
      } else if (grepl("|", cmd, fixed = TRUE)) {
        # expand growth modeling language
        cmd <- expandGrowthCmd(cmd)
      } else { # no operator, no means, must be variance.
        # cat("assuming vars: ", cmd, "\n")

        vars.lhs <- strsplit(attr(cmd, "noModifiers"), "\\s+")[[1]] # trimSpace(
        vars.rhs <- strsplit(cmd, "\\s+")[[1]] # trimSpace(

        cmd <- sapply(1:length(vars.lhs), function(v) paste(vars.lhs[v], "~~", vars.rhs[v]))
      }

      # handle threshold substitution: $ -> |
      cmd <- gsub("$", "|", cmd, fixed = TRUE)

      # if we have both starting/fixed values and constraints, these must be handled by separate commands.
      # starting and fixed values are already handled in the pipeline by this point, so should be evident in the command
      # bfi BY lab1*start(1)*bfi_1 ==> bfi BY lab1*bfi_1 + start(1)*bfi_1
      double_asterisks <- grepl("\\s*[\\w\\(\\)\\.]+\\*[\\w\\(\\)\\.]+\\*[\\w\\(\\)\\.]+", cmd, perl = TRUE)

      if (isTRUE(double_asterisks[1])) {
        ss <- strsplit(cmd, "*", fixed = TRUE)[[1]]
        if (length(ss) != 3) {
          lav_msg_warn(gettext("problem interpreting double asterisk syntax:"),
                       cmd) # sanity check on my logic
        } else {
          cmd <- paste0(ss[1], "*", ss[3], " + ", ss[2], "*", ss[3])
        }
      }

      lavaan.out <- c(lavaan.out, cmd)
    }
  }

  # new threshold syntax shifts things to the form:
  # VAR | t1 + t2 + t3 (left to write ordering)
  # Parameter labels, fixed values, and starting values are tacked on in the usual way, like
  # VAR | 5*t1 + start(1.5)*t2 + par_label*t3 (left to write ordering)

  thresh_lines <- grep("^\\s*[A-z]+\\w*\\|\\d+", lavaan.out, perl = TRUE)
  if (length(thresh_lines) > 0L) {
    thresh_vars <- unname(sub("^\\s*([A-z]+\\w*).*", "\\1", lavaan.out[thresh_lines], perl = TRUE))
    thresh_split <- split(thresh_lines, thresh_vars)
    drop_elements <- c()
    for (i in seq_along(thresh_split)) {
      this_set <- lavaan.out[thresh_split[[i]]]
      tnum <- as.integer(sub("^\\s*[A-z]+\\w*\\|(\\d+)\\s*.*", "\\1", this_set))
      this_set <- this_set[order(tnum)] # ensure that threshold numbering matches ascending order
      this_set <- sub("[^~]+\\s*~\\s*", "", this_set, perl = T) # drop variable and ~

      # convert to new t1, t2 syntax by combining modifiers with threshold numbers
      this_set <- sapply(seq_along(this_set), function(j) {
        # gsub("[^~]+\\s*~\\s*([\\w\\.\\-]+\\*)*1", paste0("\\1t", j), this_set[j], perl=TRUE)
        gsub("([\\w\\.\\-]+\\*)*1", paste0("\\1t", j), this_set[j], perl = TRUE)
      })

      new_str <- paste(names(thresh_split)[i], "|", paste(this_set, collapse = " + "))
      # replace in model string on the first line having relevant syntax
      lavaan.out[thresh_split[[i]][1]] <- new_str
      drop_elements <- c(drop_elements, thresh_split[[i]][-1])
    }
    lavaan.out <- lavaan.out[-drop_elements]
  }


  # tack on constraint syntax, if included
  lavaan.out <- c(lavaan.out, con_syntax)

  # for now, include a final trimSpace call since some arguments have leading/trailing space stripped.
  wrap <- paste(wrapAfterPlus(lavaan.out, width = 90, exdent = 5), collapse = "\n") # trimSpace(
  return(wrap)
}

mplus2lavaan <- function(inpfile, run = TRUE) {
  stopifnot(length(inpfile) == 1L)
  stopifnot(grepl("\\.inp$", inpfile, ignore.case = TRUE))
  if (!file.exists(inpfile)) {
    lav_msg_stop(gettext("Could not find file:"), inpfile)
  }

  # for future consideration. For now, require a .inp file
  #  if (length(inpfile) == 1L && grepl("\\.inp$", inpfile)) {
  #    if (!file.exists(inpfile)) { stop("Could not find file: ", inpfile) }
  #    inpfile.text <- scan(inpfile, what="character", sep="\n", strip.white=FALSE, blank.lines.skip=FALSE, quiet=TRUE)
  #  } else {
  #    #assume that inpfile itself is syntax (e.g., in a character vector)
  #    inpfile.text <- inpfile
  #  }

  inpfile.text <- scan(inpfile, what = "character", sep = "\n", strip.white = FALSE, blank.lines.skip = FALSE, quiet = TRUE)
  sections <- divideInputIntoSections(inpfile.text, inpfile)

  mplus.inp <- list()

  mplus.inp$title <- trimSpace(paste(sections$title, collapse = " "))
  mplus.inp$data <- divideIntoFields(sections$data, required = "file")
  mplus.inp$variable <- divideIntoFields(sections$variable, required = "names")
  mplus.inp$analysis <- divideIntoFields(sections$analysis)

  meanstructure <- "default" # lavaan default
  if (!is.null(mplus.inp$analysis$model)) {
    if (tolower(mplus.inp$analysis$model) == "nomeanstructure") {
      meanstructure <- FALSE
    } # explicitly disable mean structure
  }

  information <- "default" # lavaan default
  if (!is.null(mplus.inp$analysis$information)) {
    information <- tolower(mplus.inp$analysis$information)
  }

  estimator <- "default"
  if (!is.null(est <- mplus.inp$analysis$estimator)) {
    # no memory of what this is up to....
    if (toupper(est) == "MUML") lav_msg_warn(gettext(
      "Mplus does not support MUML estimator. Using default instead."))
    estimator <- est

    # march 2013: handle case where categorical data are specified, but ML-based estimator requested.
    # use WLSMV instead
    if (!is.null(mplus.inp$variable$categorical) && toupper(substr(mplus.inp$analysis$estimator, 1, 2)) == "ML") {
      lav_msg_warn(gettext(
        "Lavaan does not yet support ML-based estimation for categorical data.
        Reverting to WLSMV"))
      estimator <- "WLSMV"
    }
  }

  # expand hyphens in variable names and split into vector that will be the names for read.table
  mplus.inp$variable$names <- strsplit(expandCmd(mplus.inp$variable$names), "\\s+", perl = TRUE)[[1]]

  # expand hyphens in categorical declaration
  if (!is.null(mplus.inp$variable$categorical)) mplus.inp$variable$categorical <- strsplit(expandCmd(mplus.inp$variable$categorical), "\\s+", perl = TRUE)[[1]]

  # convert mplus syntax to lavaan syntax
  mplus.inp$model <- mplus2lavaan.modelSyntax(sections$model)

  # handle model constraint
  if ("model.constraint" %in% names(sections)) {
    mplus.inp$model.constraint <- mplus2lavaan.constraintSyntax(sections$model.constraint)
    mplus.inp$model <- paste(mplus.inp$model, mplus.inp$model.constraint, sep = "\n")
  }

  # read mplus data (and handle missing spec)
  mplus.inp$data <- readMplusInputData(mplus.inp, inpfile)

  # handle bootstrapping specification
  se <- "default"
  bootstrap <- 1000L
  test <- "default"
  if (!is.null(mplus.inp$analysis$bootstrap)) {
    boot.type <- "standard"
    # check whether standard versus residual bootstrap is specified
    if ((boot.match <- regexpr("\\((\\w+)\\)", mplus.inp$analysis$bootstrap, perl = TRUE)) > 0L) {
      boot.type <- tolower(substr(mplus.inp$analysis$bootstrap, attr(boot.match, "capture.start"), attr(boot.match, "capture.start") + attr(boot.match, "capture.length") - 1L))
    }

    if (boot.type == "residual") test <- "Bollen.Stine"

    se <- "bootstrap"

    if ((nboot.match <- regexpr("^\\s*(\\d+)", mplus.inp$analysis$bootstrap, perl = TRUE)) > 0L) {
      bootstrap <- as.numeric(substr(mplus.inp$analysis$bootstrap, attr(nboot.match, "capture.start"), attr(nboot.match, "capture.start") + attr(nboot.match, "capture.length") - 1L))
    }
  }

  if (run) {
    fit <- sem(mplus.inp$model, data = mplus.inp$data, meanstructure = meanstructure, mimic = "Mplus", estimator = estimator, test = test, se = se, bootstrap = bootstrap, information = information)
    fit@external <- list(mplus.inp = mplus.inp)
  } else {
    fit <- mplus.inp # just return the syntax outside of a lavaan object
  }

  return(fit)
}


divideIntoFields <- function(section.text, required) {
  if (is.null(section.text)) {
    return(NULL)
  }

  # The parser breaks down when there is a line with a trailing comment because then splitting on semicolon will combine it with the following line
  # Thus, trim off trailing comments before initial split
  section.text <- gsub("\\s*!.*$", "", section.text, perl = TRUE)
  section.split <- strsplit(paste(section.text, collapse = " "), ";", fixed = TRUE)[[1]] # split on semicolons
  section.divide <- list()

  for (cmd in section.split) {
    if (grepl("^\\s*!.*", cmd, perl = TRUE)) next # skip comment lines
    if (grepl("^\\s+$", cmd, perl = TRUE)) next # skip blank lines

    # mplus is apparently tolerant of specifications that don't include IS/ARE/=
    # example: usevariables x1-x10;
    # thus, split on spaces and assume that first element is lhs, drop second element if IS/ARE/=, and assume remainder is rhs

    # but if user uses equals sign, then spaces will not always be present (e.g., usevariables=x1-x10)
    if ((leadingEquals <- regexpr("^\\s*[A-Za-z]+[A-Za-z_-]*\\s*(=)", cmd[1L], perl = TRUE))[1L] > 0) {
      cmdName <- trimSpace(substr(cmd[1L], 1, attr(leadingEquals, "capture.start") - 1))
      cmdArgs <- trimSpace(substr(cmd[1L], attr(leadingEquals, "capture.start") + 1, nchar(cmd[1L])))
    } else {
      cmd.spacesplit <- strsplit(trimSpace(cmd[1L]), "\\s+", perl = TRUE)[[1L]]

      if (length(cmd.spacesplit) < 2L) {
        # for future: make room for this function to prase things like just TECH13 (no rhs)
      } else {
        cmdName <- trimSpace(cmd.spacesplit[1L])
        if (length(cmd.spacesplit) > 2L && tolower(cmd.spacesplit[2L]) %in% c("is", "are")) {
          cmdArgs <- paste(cmd.spacesplit[3L:length(cmd.spacesplit)], collapse = " ") # remainder, removing is/are
        } else {
          cmdArgs <- paste(cmd.spacesplit[2L:length(cmd.spacesplit)], collapse = " ") # is/are not used, so just join rhs
        }
      }
    }

    section.divide[[make.names(tolower(cmdName))]] <- cmdArgs
  }

  if (!missing(required)) {
    stopifnot(all(required %in% names(section.divide)))
  }
  return(section.divide)
}

# helper function
splitFilePath <- function(abspath) {
  # function to split path into path and filename
  # code adapted from R.utils filePath command
  if (!is.character(abspath)) lav_msg_stop(gettext(
    "Path not a character string"))
  if (nchar(abspath) < 1 || is.na(abspath)) lav_msg_stop(gettext(
    "Path is missing or of zero length"))

  components <- strsplit(abspath, split = "[\\/]")[[1]]
  lcom <- length(components)

  stopifnot(lcom > 0)

  # the file is the last element in the list. In the case of length == 1, this will extract the only element.
  relFilename <- components[lcom]
  absolute <- FALSE

  if (lcom == 1) {
    dirpart <- NA_character_
  } else if (lcom > 1) {
    # drop the file from the list (the last element)
    components <- components[-lcom]
    dirpart <- do.call("file.path", as.list(components))

    # if path begins with C:, /, //, or \\, then treat as absolute
    if (grepl("^([A-Z]{1}:|/|//|\\\\)+.*$", dirpart, perl = TRUE)) absolute <- TRUE
  }

  return(list(directory = dirpart, filename = relFilename, absolute = absolute))
}

readMplusInputData <- function(mplus.inp, inpfile) {
  # handle issue of mplus2lavaan being called with an absolute path, whereas mplus has only a local data file
  inpfile.split <- splitFilePath(inpfile)
  datfile.split <- splitFilePath(mplus.inp$data$file)

  # if inp file target directory is non-empty, but mplus data is without directory, then append
  # inp file directory to mplus data. This ensures that R need not be in the working directory
  # to read the dat file. But if mplus data has an absolute directory, don't append

  # if mplus data directory is present and absolute, or if no directory in input file, just use filename as is
  if (!is.na(datfile.split$directory) && datfile.split$absolute) {
    datFile <- mplus.inp$data$file
  } # just use mplus data filename if it has absolute path
  else if (is.na(inpfile.split$directory)) {
    datFile <- mplus.inp$data$file
  } # just use mplus data filename if inp file is missing path (working dir)
  else {
    datFile <- file.path(inpfile.split$directory, mplus.inp$data$file)
  } # dat file path is relative or absent, and inp file directory is present

  if (!file.exists(datFile)) {
    lav_msg_warn(gettext("Cannot find data file:"), datFile)
    return(NULL)
  }

  # handle missing is/are:
  missList <- NULL
  if (!is.null(missSpec <- mplus.inp$variable$missing)) {
    expandMissVec <- function(missStr) {
      # sub-function to obtain a vector of all missing values within a set of parentheses
      missSplit <- strsplit(missStr, "\\s+")[[1L]]
      missVals <- c()
      for (f in missSplit) {
        if ((hyphenPos <- regexpr("\\d+(-)\\d+", f, perl = TRUE))[1L] > -1L) {
          # expand hyphen
          preHyphen <- substr(f, 1, attr(hyphenPos, "capture.start") - 1)
          postHyphen <- substr(f, attr(hyphenPos, "capture.start") + 1, nchar(f))
          missVals <- c(missVals, as.character(seq(preHyphen, postHyphen)))
        } else {
          # append to vector
          missVals <- c(missVals, f)
        }
      }
      return(as.numeric(missVals))
    }

    if (missSpec == "." || missSpec == "*") { # case 1: MISSING ARE|=|IS .;
      na.strings <- missSpec
    } else if ((allMatch <- regexpr("\\s*ALL\\s*\\(([^\\)]+)\\)", missSpec, perl = TRUE))[1L] > -1L) { # case 2: use of ALL with parens
      missStr <- trimSpace(substr(missSpec, attr(allMatch, "capture.start"), attr(allMatch, "capture.start") + attr(allMatch, "capture.length") - 1L))
      na.strings <- expandMissVec(missStr)
    } else { # case 3: specific missing values per variable
      # process each element
      missBlocks <- gregexpr("(?:(\\w+)\\s+\\(([^\\)]+)\\))+", missSpec, perl = TRUE)[[1]]
      missList <- list()

      if (missBlocks[1L] > -1L) {
        for (i in 1:length(missBlocks)) {
          vname <- substr(missSpec, attr(missBlocks, "capture.start")[i, 1L], attr(missBlocks, "capture.start")[i, 1L] + attr(missBlocks, "capture.length")[i, 1L] - 1L)
          vmiss <- substr(missSpec, attr(missBlocks, "capture.start")[i, 2L], attr(missBlocks, "capture.start")[i, 2L] + attr(missBlocks, "capture.length")[i, 2L] - 1L)

          vnameHyphen <- regexpr("(\\w+)-(\\w+)", vname, perl = TRUE)[1L]
          if (vnameHyphen > -1L) {
            # lookup against variable names
            vstart <- which(mplus.inp$variable$names == substr(vname, attr(vnameHyphen, "capture.start")[1L], attr(vnameHyphen, "capture.start")[1L] + attr(vnameHyphen, "capture.length")[1L] - 1L))
            vend <- which(mplus.inp$variable$names == substr(vname, attr(vnameHyphen, "capture.start")[2L], attr(vnameHyphen, "capture.start")[2L] + attr(vnameHyphen, "capture.length")[2L] - 1L))
            if (length(vstart) == 0L || length(vend) == 0L) {
              lav_msg_stop(gettext("Unable to lookup missing variable list: "),
                           vname)
            }
            # I suppose start or finish could be mixed up
            if (vstart > vend) {
              vstart.orig <- vstart
              vstart <- vend
              vend <- vstart.orig
            }
            vname <- mplus.inp$variable$names[vstart:vend]
          }

          missVals <- expandMissVec(vmiss)

          for (j in 1:length(vname)) {
            missList[[vname[j]]] <- missVals
          }
        }
      } else {
        lav_msg_stop(gettext("I don't understand this missing specification:"),
                     missSpec)
      }
    }
  } else {
    na.strings <- "NA"
  }

  if (!is.null(missList)) {
    dat <- read.table(datFile, header = FALSE, col.names = mplus.inp$variable$names, colClasses = "numeric")
    # loop over variables in missList and set missing values to NA
    dat[, names(missList)] <- lapply(names(missList), function(vmiss) {
      dat[which(dat[, vmiss] %in% missList[[vmiss]]), vmiss] <- NA
      return(dat[, vmiss])
    })

    names(dat) <- mplus.inp$variable$names # loses these from the lapply
  } else {
    dat <- read.table(datFile, header = FALSE, col.names = mplus.inp$variable$names, na.strings = na.strings, colClasses = "numeric")
  }


  # TODO: support covariance/mean+cov inputs

  # store categorical variables as ordered factors
  if (!is.null(mplus.inp$variable$categorical)) {
    dat[, c(mplus.inp$variable$categorical)] <- lapply(dat[, c(mplus.inp$variable$categorical), drop = FALSE], ordered)
  }

  return(dat)
}


divideInputIntoSections <- function(inpfile.text, filename) {
  inputHeaders <- grep("^\\s*(title:|data.*:|variable:|define:|analysis:|model.*:|output:|savedata:|plot:|montecarlo:)", inpfile.text, ignore.case = TRUE, perl = TRUE)

  stopifnot(length(inputHeaders) > 0L)

  mplus.sections <- list()

  for (h in 1:length(inputHeaders)) {
    sectionEnd <- ifelse(h < length(inputHeaders), inputHeaders[h + 1] - 1, length(inpfile.text))
    section <- inpfile.text[inputHeaders[h]:sectionEnd]
    sectionName <- trimSpace(sub("^([^:]+):.*$", "\\1", section[1L], perl = TRUE)) # obtain text before the colon

    # dump section name from input syntax
    section[1L] <- sub("^[^:]+:(.*)$", "\\1", section[1L], perl = TRUE)

    mplus.sections[[make.names(tolower(sectionName))]] <- section
  }

  return(mplus.sections)
}
