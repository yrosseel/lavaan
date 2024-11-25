# -------------------- main parsing function in C or R ----------------------- #
ldw_parse_model_string_cr <- function(model.syntax = "",
                                      as.data.frame. = FALSE) {
  stopifnot(length(model.syntax) > 0L)
  # replace 'strange' tildes (in some locales) (new in 0.6-6)
  modelsrc <- gsub(
    pattern = "\u02dc",
    replacement = "~",
    paste(unlist(model.syntax), "", collapse = "\n")
  )

  if (requireNamespace("lavaanC", quietly = TRUE)) {
    flat <- lavaanC::lav_parse_model_string_c(modelsrc)
  } else {
    flat <- lav_parse_model_string_r(modelsrc)
  }

  if (!is.list(flat)) {
    if (!is.integer(flat) || length(flat) != 2L) {
      lav_msg_stop(gettext("internal error in lav_parse_model_string_c!"))
    } else if (flat[1] == 1L) { # SPE_MALLOC
      lav_msg_stop(gettextf("Error allocating memory on line %d of C code.",
                            flat[2L]))
    } else if (flat[1] == 2L) { # SPE_PROGERROR
      lav_msg_stop(gettextf("Internal program error on line %d of C code.",
                            flat[2L]))
    } else {
      tl <- ldw_txtloc(modelsrc, flat[2L] + 1L)
      if (flat[1L] == 21L) { # SPE_ILLNUMLIT
        lav_msg_stop(gettext("Illegal numeric literal"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 22L) { # SPE_EMPTYMODEL
        lav_msg_stop(gettext("The model specification seems empty!"))
      } else if (flat[1L] == 23L) { # SPE_FORMUL1
        lav_msg_stop(gettext("Formula with only 1 token!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 24L) { # SPE_ILLCHAR
        lav_msg_stop(gettext("Illegal character in model source!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 31L) { # SPE_NOOPERATOR
        lav_msg_stop(gettext("formula without valid operator"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 32L) { # SPE_PARENTHESES
        lav_msg_stop(gettext("Formula with non-matching parentheses!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 33L) { # SPE_3WAYINTERACTION
        lav_msg_stop(gettext(
        "Three-way or higher-order interaction terms (using multiple
         colons) are not supported in the lavaan syntax; please manually
         construct the product terms yourself in the data.frame, give
         them an appropriate name, and then you can use these interaction
         variables as any other (observed) variable in the model syntax."),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 34L) { # SPE_AUTOREGRESS
        lav_msg_stop(gettext("a variable cannot be regressed on itself"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 41L) { # SPE_INVALIDNAME
        lav_msg_stop(gettext("identifier is either a reserved word (in R) or
              contains an illegal character"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 42L) { # SPE_INVALIDLHS
        lav_msg_stop(gettext("invalid left hand side modifier"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 43L) { # SPE_INVALIDVECTOR
        lav_msg_stop(gettext("invalid vector specification"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 44L) { # SPE_MODNOLITORID
        lav_msg_stop(gettext(
          "invalid value (should be numeric, identifier or string)"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 45L) { # SPE_MODNONUM
        lav_msg_stop(gettext("invalid value (should be numeric)"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 46L) { # SPE_MODNOSTR
        lav_msg_stop(gettext("invalid value (should be string)"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 47L) { # SPE_INVALIDBLOCK
        lav_msg_stop(gettext(
        "Invalid block identifier. The correct syntax is: \"LHS: RHS\",
          where LHS is a block identifier (eg group or level), and RHS is
          the group/level/block number or label."),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 48L) { # SPE_INVALIDLAST
        lav_msg_stop(gettext(
          "Last element of rhs part expected to be an identifier or,
                for operator ~ or =~, a numeric literal!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 49L) { # SPE_INVALIDMODSMB
        lav_msg_stop(gettext(
          "Invalid modifier symbol (should be '*' or '?')!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 50L) { # SPE_INVALIDEXPR
        lav_msg_stop(gettext("Invalid expression!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 51L) { # SPE_INVALIDEXPRTYP
        lav_msg_stop(gettext("Incorrrect type for expression!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (flat[1L] == 52L) { # SPE_LVLGRP
        lav_msg_stop(gettext("groups can not be nested within levels!"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else {
        lav_msg_stop(gettext(
          "internal error (code unknown) in lav_parse_model_string_cr!"))
      }
    }
    return(NULL)
  }
  constraints <- attr(flat, "constraints", TRUE)
  mod <- attr(flat, "modifiers", TRUE)
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

  if (!is.null(attr(flat, "warns"))) {
    warns <- attr(flat, "warns")
    attr(flat, "warns") <- NULL
    for (w in warns) {
      tl <- ldw_txtloc(modelsrc, w[2L] + 1L)
      if (w[1L] == 101L) { # SPW_OPERATORBLANKS
        lav_msg_warn(gettext("Blancs in lavaan operators are deprecated."),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (w[1L] == 102L) { # SPW_IDENTIFIERBLANKS
        lav_msg_warn(gettext("Blancs in identifiers are deprecated."),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (w[1L] == 103L) { # SPW_FIRSTBLK
        lav_msg_warn(gettext("First block defined after other formula's"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (w[1L] == 104L) { # SPW_MODMULTIPLE
        lav_msg_warn(gettext("modifier specified multiple times, overwritten"),
                     tl[1L],
                     footer = tl[2L]
        )
      } else if (w[1L] == 105L) { # SPW_MODMULTIPLE
        lav_msg_warn(gettext("syntax contains only a single block identifier!"),
                     tl[1L],
                     footer = tl[2L]
        )
      }
    }
  }

  if (as.data.frame.) {
    flat <- as.data.frame(flat, stringsAsFactors = FALSE)
  }
  # create output
  attr(flat, "modifiers") <- mod
  attr(flat, "constraints") <- constraints
  flat
}
