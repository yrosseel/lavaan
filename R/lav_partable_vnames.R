# lav_partable_names
#
# YR. 29 june 2013
#  - as separate file; used to be in utils-user.R
#  - lav_partable_names (aka 'vnames') allows multiple options in 'type'
#    returning them all as a list (or just a vector if only 1 type is needed)

# public version
lavNames <- function(object, type = "ov", ...) { # nolint
  if (inherits(object, "lavaan") || inherits(object, "lavaanList")) {
    partable <- object@ParTable
  } else if (inherits(object, "list") ||
    inherits(object, "data.frame")) {
    partable <- object
  } else if (inherits(object, "character")) {
    # just a model string?
    partable <- lavParseModelString(object)
  }
  lav_partable_vnames(partable, type = type, ...)
}

# alias for backwards compatibility
lavaanNames <- lavNames # nolint

# return variable names in a partable
# - the 'type' argument determines the status of the variable: observed,
#   latent, endo/exo/...; default = "ov", but most used is type = "all"
#   LDW 30/1/24: there is no default and most used is a single type: "all" is
#                rarely used in other code
# - the 'group' argument either selects a single group (if group is an integer)
#   or returns a list per group
# - the 'level' argument either selects a single level (if level is an integer)
#   or returns a list per level
#   LDW 30/1/24: 'level' argument not explicitly tested !?
# - the 'block' argument either selects a single block (if block is an integer)
#   or returns a list per block
#   LDW 30/1/24: 'block' argument not explicitly tested !?
#   LDW 29/2/24: ov.order = "data" via attribute "ovda"
lav_partable_vnames <- function(partable, type = NULL, ..., # nolint
                                warn = FALSE, ov.x.fatal = FALSE) {
  # This function derives the names of some types of variable (as specified
  # in type) from a 'partable'. The 'warn' parameter needs no explanation.
  # The ov.x.fatal parameter implies, when set to TRUE, that the function
  # issues a 'stop' when there are exogenious variables present in variance/
  # covariance or intercept formulas.
  # The call of this function can also contain extra parameters (...) which
  # have to be the name(s) of blockvariable(s) to be used to select names.
  # If more than 1 blockvariable given, all must be satisfied to select!
  # The 'partable' must be a list with minimum members lhs, op, rhs.
  # Other members of 'partable' used (if present): block, ustart, free, exo,
  #                                                user, efa, rv, 'blockname'
  # If the 'partable' contains an attribute vnames and the type is not "*",
  # the 'type' elements of this attribute are used.

  # ----- lav_partable_vnames ---- common ----------------------------------
  # sanity check
  stopifnot(is.list(partable), !missing(type))
  # check for empty table
  if (length(partable$lhs) == 0) {
    return(character(0L))
  }
  # dotdotdot
  dotdotdot <- list(...)
  type.list <- c(
    "ov", # observed variables (ov)
    "ov.x", # (pure) exogenous observed variables
    "ov.nox", # non-exogenous observed variables
    "ov.model", # modeled observed variables (joint vs cond)
    "ov.y", # (pure) endogenous variables (dependent only)
    "ov.num", # numeric observed variables
    "ov.ord", # ordinal observed variables
    "ov.ind", # observed indicators of latent variables
    "ov.orphan", # lonely observed intercepts/variances
    "ov.interaction", # interaction terms (with colon)
    "ov.efa", # indicators involved in efa

    "th", # thresholds ordinal only
    "th.mean", # thresholds ordinal + numeric variables

    "lv", # latent variables
    "lv.regular", # latent variables (defined by =~ only)
    "lv.formative", # latent variables (defined by <~ only)
    "lv.x", # (pure) exogenous variables
    "lv.y", # (pure) endogenous variables
    "lv.nox", # non-exogenous latent variables
    "lv.nonnormal", # latent variables with non-normal indicators
    "lv.interaction", # interaction terms
    "lv.efa", # latent variables involved in efa
    "lv.rv", # random slopes, random variables
    "lv.ind", # latent indicators (higher-order cfa)
    "lv.marker", # marker indicator per lv

    "eqs.y", # y's in regression
    "eqs.x" # x's in regression
  )
  if (type[1L] != "all" && type[1L] != "*" && !all(type %in% type.list)) {
    wrongtypes <- type[!(type %in% type.list)]
    lav_msg_stop(sprintf(
      ngettext(length(wrongtypes),
      "type = %s is not a valid option",
      "type = %s are not valid options"),
      lav_msg_view(wrongtypes, "none", FALSE)))
  }
  return.value <- NULL
  if (type[1L] != "*" && !is.null(attr(partable, "vnames"))) {
    # ----- lav_partable_vnames ---- cached data --------------------------
    # uncomment/comment following line to enable/disable trace
    # ldw_trace(paste("cached:", paste(type, collapse = ",")))

    if (type[1L] == "all") {
      return.value <- attr(partable, "vnames")
    } else {
      return.value <- attr(partable, "vnames")[type]
    }
  }
  # ----- lav_partable_vnames ---- common ----------------------------------
  if (type[1L] == "all" || type[1L] == "*") {
    type <- type.list
  }
  # ALWAYS need `block' column -- create one if missing
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$lhs))
  }
  # per default, use full partable
  block.select <- lav_partable_block_values(partable)
  # check for ... selection argument(s)
  ndotdotdot <- length(dotdotdot)
  if (ndotdotdot > 0L) {
    dot.names <- names(dotdotdot)
    row.select <- rep(TRUE, length(partable$lhs))
    for (dot in seq_len(ndotdotdot)) {
      # selection variable?
      block.var <- dot.names[dot]
      block.val <- dotdotdot[[block.var]]
      # do we have this 'block.var' in partable?
      if (is.null(partable[[block.var]])) {
        # for historical reasons, treat "group = 1" special
        if (block.var == "group" && block.val == 1L) {
          partable$group <- rep(1L, length(partable$lhs))
          # remove block == 0
          idx <- which(partable$block == 0L)
          if (length(idx) > 0L) {
            partable$group[idx] <- 0L
          }
          row.select <- (row.select &
            partable[[block.var]] %in% block.val)
        } else {
          lav_msg_stop(gettextf(
            "selection variable '%s' not found in the parameter table.",
            block.var))
        }
      } else {
        if (!all(block.val %in% partable[[block.var]])) {
          lav_msg_stop(gettextf(
            "%1$s column does not contain value `%2$s'", block.var, block.val))
        }
        row.select <- (row.select &
          !partable$op %in% c("==", "<", ">", ":=") &
          partable[[block.var]] %in% block.val)
      }
    } # dot
    block.select <- unique(partable$block[row.select])
    if (length(block.select) == 0L) {
      lav_msg_warn(gettext("no blocks selected."))
    }
  }
  if (is.null(return.value)) {
    # ----- lav_partable_vnames ---- no cache ------------------------

    # uncomment/comment following line to enable/disable trace
    # ldw_trace(paste("computed:", paste(type, collapse = ",")))

    # random slope names, if any (new in 0.6-7)
    if (!is.null(partable$rv) && any(nchar(partable$rv) > 0L)) {
      rv.names <- unique(partable$rv[nchar(partable$rv) > 0L])
    } else {
      rv.names <- character(0L)
    }

    # output: list per block
    return.value <- lapply(type, function(x) {
      vector("list", length = length(block.select))
    })
    names(return.value) <- type

    for (b in block.select) {
      # indices for this block
      block.ind <- partable$block == b

      # always compute lv.names
      lv.names <- unique(partable$lhs[block.ind &
        (partable$op == "=~" | partable$op == "<~")])
      # including random slope names
      lv.names2 <- unique(c(lv.names, rv.names))

      # determine lv interactions
      int.names <- unique(partable$rhs[block.ind &
        grepl(":", partable$rhs, fixed = TRUE)])
      n.int <- length(int.names)
      if (n.int > 0L) {
        ok.idx <- logical(n.int)
        for (iv in seq_len(n.int)) {
          tmp.names <- strsplit(int.names[iv], ":", fixed = TRUE)[[1L]]
          # three scenario's:
          # - both variables are latent (ok)
          # - both variables are observed (ignore)
          # - only one latent (warn??) -> upgrade observed to latent
          # thus if at least one is in lv.names, we treat it as a
          # latent interaction
          if (any(tmp.names %in% lv.names)) {
            ok.idx[iv] <- TRUE
          }
        }
        lv.interaction <- int.names[ok.idx]
        lv.names <- c(lv.names, lv.interaction)
        lv.names2 <- c(lv.names2, lv.interaction)
      } else {
        lv.interaction <- character(0L)
      }
      if (length(type) == 1L) {
        # ----- lav_partable_vnames ---- no cache ----- 1 type -----------
        # store lv
        if ("lv" == type) {
          # check if FLAT for random slopes
          # if (!is.null(partable$rv) && any(nchar(partable$rv) > 0L) &&
          #    !is.null(partable$block)) {
          #    return.value$lv[[b]] <- lv.names2
          # } else {
          # here, they will be 'defined' at level 2 as regular =~ lvs
          return.value$lv[[b]] <- lv.names
          # }
          next
        }

        # regular latent variables ONLY (ie defined by =~ only)
        if ("lv.regular" == type) {
          out <- unique(partable$lhs[block.ind &
            partable$op == "=~" &
            !partable$lhs %in% rv.names])
          return.value$lv.regular[[b]] <- out
          next
        }

        # interaction terms involving latent variables (only)
        if ("lv.interaction" == type) {
          return.value$lv.interaction[[b]] <- lv.interaction
          next
        }

        # formative latent variables ONLY (ie defined by <~ only)
        if ("lv.formative" == type) {
          out <- unique(partable$lhs[block.ind &
            partable$op == "<~"])
          return.value$lv.formative[[b]] <- out
          next
        }

        # lv's involved in efa
        if (any(type == c("lv.efa", "ov.efa"))) {
          if (is.null(partable$efa)) {
            out <- character(0L)
          } else {
            set.names <- lav_partable_efa_values(partable)
            out <- unique(partable$lhs[partable$op == "=~" &
              block.ind &
              partable$efa %in% set.names])
          }
          if (type == "ov.efa") ov_efa <- out
          if (type == "lv.efa") {
            return.value$lv.efa[[b]] <- out
            next
          }
        }

        # lv's that are random slopes
        if ("lv.rv" == type) {
          if (is.null(partable$rv)) {
            out <- character(0L)
          } else {
            out <- unique(partable$lhs[partable$op == "=~" &
              block.ind &
              partable$lhs %in% rv.names])
          }
          return.value$lv.rv[[b]] <- out
          next
        }

        # lv's that are indicators of a higher-order factor
        if ("lv.ind" == type) {
          out <- unique(partable$rhs[block.ind &
            partable$op == "=~" &
            partable$rhs %in% lv.names])
          return.value$lv.ind[[b]] <- out
          next
        }

        # eqs.y
        if (!(any(type == c("lv", "lv.regular")))) {
          eqs.y <- unique(partable$lhs[block.ind &
            partable$op == "~"])
        }

        # store eqs.y
        if ("eqs.y" == type) {
          return.value$eqs.y[[b]] <- eqs.y
          next
        }

        # eqs.x
        if (!(any(type == c("lv", "lv.regular", "lv.x")))) {
          eqs.x <- unique(partable$rhs[block.ind &
            (partable$op == "~" |
              partable$op == "<~")])
        }

        # store eqs.x
        if ("eqs.x" == type) {
          return.value$eqs.x[[b]] <- eqs.x
          next
        }

        # v.ind -- indicators of latent variables
        if (!(any(type == c("lv", "lv.regular")))) {
          v.ind <- unique(partable$rhs[block.ind &
            partable$op == "=~"])
        }

        # ov.*
        if (!(any(type == c("lv", "lv.regular", "lv.x", "lv.y")))) {
          # 1. indicators, which are not latent variables themselves
          ov.ind <- v.ind[!v.ind %in% lv.names2]
          # 2. dependent ov's
          ov.y <- eqs.y[!eqs.y %in% c(lv.names2, ov.ind)]
          # 3. independent ov's
          if (lav_partable_nlevels(partable) > 1L && b > 1L) {
            # NEW in 0.6-8: if an 'x' was an 'y' in a previous level,
            #               treat it as 'y'
            tmp.eqs.y <- unique(partable$lhs[partable$op == "~"]) # all blocks
            ov.x <- eqs.x[!eqs.x %in% c(lv.names2, ov.ind, tmp.eqs.y)]
          } else {
            ov.x <- eqs.x[!eqs.x %in% c(lv.names2, ov.ind, ov.y)]
          }
          # new in 0.6-12: if we have interaction terms in ov.x, check
          # if some terms are in eqs.y; if so, remove the interaction term
          # from ov.x
          int.idx <- which(grepl(":", ov.x, fixed = TRUE))
          bad.idx <- integer(0L)
          for (iv in int.idx) {
            tmp.names <- strsplit(ov.x[iv], ":", fixed = TRUE)[[1L]]
            if (any(tmp.names %in% eqs.y)) {
              bad.idx <- c(bad.idx, iv)
            }
          }
          if (length(bad.idx) > 0L) {
            ov.y <- unique(c(ov.y, ov.x[bad.idx]))
            # it may be removed later, but needed to construct ov.names
            ov.x <- ov.x[-bad.idx]
          }
        }

        # observed variables
        # easy approach would be: everything that is not in lv.names,
        # but the main purpose here is to 'order' the observed variables
        # according to 'type' (indicators, ov.y, ov.x, orphans)
        if (!(any(type == c("lv", "lv.regular", "lv.x", "lv.y")))) {
          # 4. orphaned covariances
          ov.cov <- c(
            partable$lhs[block.ind &
              partable$op == "~~" &
              !partable$lhs %in% lv.names2],
            partable$rhs[block.ind &
              partable$op == "~~" &
              !partable$rhs %in% lv.names2]
          )
          # 5. orphaned intercepts/thresholds
          ov.int <- partable$lhs[block.ind &
            (partable$op == "~1" |
              partable$op == "|") &
            !partable$lhs %in% lv.names2]

          ov.tmp <- c(ov.ind, ov.y, ov.x)
          ov.extra <- unique(c(ov.cov, ov.int)) # must be in this order!
          # so that
          # lav_partable_independence
          # retains the same order
          ov.names <- c(ov.tmp, ov.extra[!ov.extra %in% ov.tmp])
        }

        # store ov?
        if ("ov" == type) {
          return.value$ov[[b]] <- ov.names
          next
        }

        if ("ov.ind" == type) {
          return.value$ov.ind[[b]] <- ov.ind
          next
        }

        if ("ov.interaction" == type) {
          ov.int.names <- ov.names[grepl(":", ov.names, fixed = TRUE)]
          n.int <- length(ov.int.names)
          if (n.int > 0L) {
            ov.names.noint <- ov.names[!ov.names %in% ov.int.names]

            ok.idx <- logical(n.int)
            for (iv in seq_len(n.int)) {
              tmp.names <- strsplit(ov.int.names[iv], ":", fixed = TRUE)[[1L]]

              # two scenario's:
              # - both variables are in ov.names.noint (ok)
              # - at least one variables is NOT in ov.names.noint (ignore)
              if (all(tmp.names %in% ov.names.noint)) {
                ok.idx[iv] <- TRUE
              }
            }
            ov.interaction <- ov.int.names[ok.idx]
          } else {
            ov.interaction <- character(0L)
          }

          return.value$ov.interaction[[b]] <- ov.interaction
          next
        }

        if ("ov.efa" == type) {
          ov.efa <- partable$rhs[partable$op == "=~" &
            block.ind &
            partable$rhs %in% ov.ind &
            partable$lhs %in% ov_efa]
          return.value$ov.efa[[b]] <- unique(ov.efa)
          next
        }


        # exogenous `x' covariates
        if (any(type == c(
          "ov.x", "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          # correction: is any of these ov.names.x mentioned as a variance,
          #             covariance, or intercept?
          # this should trigger a warning in lavaanify()
          if (is.null(partable$user)) { # FLAT!
            partable$user <- rep(1L, length(partable$lhs))
          }
          vars <- c(
            partable$lhs[block.ind &
              partable$op == "~1" &
              partable$user == 1],
            partable$lhs[block.ind &
              partable$op == "~~" &
              partable$user == 1],
            partable$rhs[block.ind &
              partable$op == "~~" &
              partable$user == 1]
          )
          idx.no.x <- which(ov.x %in% vars)
          if (length(idx.no.x)) {
            if (ov.x.fatal) {
              lav_msg_stop(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s];  Please remove them
                and try again.", lav_msg_view(ov.x[idx.no.x], "none")))
            }
            if (warn) {
              lav_msg_warn(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s]; these variables will
                now be treated as random introducing additional free parameters.
                If you wish to treat those variables as fixed, remove these
                formulas from the model syntax. Otherwise, consider adding the
                fixed.x = FALSE option.", lav_msg_view(ov.x[idx.no.x], "none")))
            }
            ov.x <- ov.x[-idx.no.x]
          }
          ov.tmp.x <- ov.x

          # extra
          if (!is.null(partable$exo)) {
            ov.cov <- c(
              partable$lhs[block.ind &
                partable$op == "~~" &
                partable$exo == 1L],
              partable$rhs[block.ind &
                partable$op == "~~" &
                partable$exo == 1L]
            )
            ov.int <- partable$lhs[block.ind &
              partable$op == "~1" &
              partable$exo == 1L]
            ov.extra <- unique(c(ov.cov, ov.int))
            ov.tmp.x <- c(ov.tmp.x, ov.extra[!ov.extra %in% ov.tmp.x])
          }

          ov.names.x <- ov.tmp.x
        }

        # store ov.x?
        if ("ov.x" == type) {
          return.value$ov.x[[b]] <- ov.names.x
          next
        }

        # story ov.orphan?
        if ("ov.orphan" == type) {
          return.value$ov.orphan[[b]] <- ov.extra
          next
        }

        # ov's withouth ov.x
        if (any(type == c(
          "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          ov.names.nox <- ov.names[!ov.names %in% ov.names.x]
        }

        # store ov.nox
        if ("ov.nox" == type) {
          return.value$ov.nox[[b]] <- ov.names.nox
          next
        }

        # store ov.model
        if ("ov.model" == type) {
          # if no conditional.x, this is just ov
          # else, this is ov.nox
          if (any(block.ind & partable$op == "~" &
            partable$exo == 1L)) {
            return.value$ov.model[[b]] <- ov.names.nox
          } else {
            return.value$ov.model[[b]] <- ov.names
          }
          next
        }

        # ov's strictly ordered
        if (any(type == c(
          "ov.ord", "th", "th.mean",
          "ov.num", "lv.nonnormal"
        ))) {
          tmp <- unique(partable$lhs[block.ind &
            partable$op == "|"])
          ord.names <- ov.names[ov.names %in% tmp]
        }

        if ("ov.ord" == type) {
          return.value$ov.ord[[b]] <- ord.names
          next
        }

        # ov's strictly numeric
        if (any(type == c("ov.num", "lv.nonnormal"))) {
          ov.num <- ov.names[!ov.names %in% ord.names]
        }

        if ("ov.num" == type) {
          return.value$ov.num[[b]] <- ov.num
          next
        }

        # nonnormal lv's
        if ("lv.nonnormal" == type) {
          # regular lv's
          lv.reg <- unique(partable$lhs[block.ind &
            partable$op == "=~"])
          if (length(lv.reg) > 0L) {
            out <- unlist(lapply(lv.reg, function(x) {
              # get indicators for this lv
              tmp.ind <- unique(partable$rhs[block.ind &
                partable$op == "=~" &
                partable$lhs == x])
              if (!all(tmp.ind %in% ov.num)) {
                return(x)
              } else {
                return(character(0))
              }
            }), use.names = FALSE)
            return.value$lv.nonnormal[[b]] <- out
          } else {
            return.value$lv.nonnormal[[b]] <- character(0)
          }
          next
        }

        if (any(c("th", "th.mean") == type)) {
          tmp.th.lhs <- partable$lhs[block.ind &
            partable$op == "|"]
          tmp.th.rhs <- partable$rhs[block.ind &
            partable$op == "|"]
        }

        # threshold
        if ("th" == type) {
          if (length(ord.names) > 0L) {
            # return in the right order (following ord.names!)
            out <- unlist(lapply(ord.names, function(x) {
              idx <- which(x == tmp.th.lhs)
              tmp.th <- unique(paste(tmp.th.lhs[idx], "|",
                tmp.th.rhs[idx],
                sep = ""
              ))
              # make sure the th's are in increasing order
              # sort(tmp.th)
              # NO!, don't do that; t10 will be before t2
              # fixed in 0.6-1 (bug report from Myrsini)
              # in 0.6-12, we do this anyway like this:

              # get var name
              tmp.th1 <- sapply(
                strsplit(tmp.th, split = "\\|t"),
                "[[", 1
              )
              # get number, and sort
              tmp.th2 <- as.character(sort(as.integer(sapply(
                strsplit(tmp.th, split = "\\|t"), "[[", 2
              ))))
              # paste back togehter in the right order
              paste(tmp.th1, tmp.th2, sep = "|t")
            }), use.names = FALSE)
          } else {
            out <- character(0L)
          }
          return.value$th[[b]] <- out
          next
        }

        # thresholds and mean/intercepts of numeric variables
        if ("th.mean" == type) {
          # if fixed.x -> use ov.names.nox
          # else -> use ov.names
          if (is.null(partable$exo) || all(partable$exo == 0L)) {
            tmp.ov.names <- ov.names
          } else {
            tmp.ov.names <- ov.names.nox
          }
          if (length(tmp.ov.names) > 0L) {
            # return in the right order (following ov.names.nox!)
            out <- unlist(lapply(tmp.ov.names, function(x) {
              if (x %in% ord.names) {
                idx <- which(x == tmp.th.lhs)
                tmp.th <- unique(
                  paste(tmp.th.lhs[idx], "|",
                    tmp.th.rhs[idx],
                    sep = ""
                  ),
                  use.names = FALSE
                )
                # make sure the th's are in increasing order
                # get var name
                tmp.th1 <- sapply(
                  strsplit(tmp.th, split = "\\|t"),
                  "[[", 1
                )
                # get number, and sort
                tmp.th2 <- as.character(sort(as.integer(sapply(
                  strsplit(tmp.th, split = "\\|t"), "[[", 2
                ))))
                # paste back togehter in the right order
                paste(tmp.th1, tmp.th2, sep = "|t")
              } else {
                x
              }
            }))
          } else {
            out <- character(0L)
          }
          return.value$th.mean[[b]] <- out
          next
        }

        # exogenous lv's
        if (any(c("lv.x", "lv.nox") == type)) {
          tmp <- lv.names[!lv.names %in% c(v.ind, eqs.y)]
          lv.names.x <- lv.names[lv.names %in% tmp]
        }

        if ("lv.x" == type) {
          return.value$lv.x[[b]] <- lv.names.x
          next
        }

        # dependent ov (but not also indicator or x)
        if ("ov.y" == type) {
          tmp <- eqs.y[!eqs.y %in% c(v.ind, eqs.x, lv.names)]
          return.value$ov.y[[b]] <- ov.names[ov.names %in% tmp]
          next
        }

        # dependent lv (but not also indicator or x)
        if ("lv.y" == type) {
          tmp <- eqs.y[!eqs.y %in% c(v.ind, eqs.x) &
            eqs.y %in% lv.names]
          return.value$lv.y[[b]] <- lv.names[lv.names %in% tmp]
          next
        }

        # non-exogenous latent variables
        if ("lv.nox" == type) {
          return.value$lv.nox[[b]] <- lv.names[!lv.names %in% lv.names.x]
          next
        }

        # marker indicator (if any) for each lv
        if ("lv.marker" == type) {
          # default: "" per lv
          out <- character(length(lv.names))
          names(out) <- lv.names
          for (l in seq_len(length(lv.names))) {
            this.lv.name <- lv.names[l]
            # try to see if we can find a 'marker' indicator for this factor
            marker.idx <- which(block.ind &
              partable$op == "=~" &
              partable$lhs == this.lv.name &
              partable$rhs %in% v.ind &
              partable$ustart == 1L &
              partable$free == 0L)
            if (length(marker.idx) == 1L) { # unique only!!
              # check if 'other' loadings are fixed to zero
              other.idx <- which(block.ind &
                partable$op == "=~" &
                partable$lhs != this.lv.name &
                partable$rhs == partable$rhs[marker.idx] &
                partable$free == 0L)
              if (length(other.idx) == 0L) {
                out[l] <- partable$rhs[marker.idx]
              } else if (all(partable$ustart[other.idx] == 0)) {
                out[l] <- partable$rhs[marker.idx]
              }
            }
          }
          return.value$lv.marker[[b]] <- out
        }
      } else {
        # ----- lav_partable_vnames ---- no cache ----- more than 1 type ------
        # store lv
        if (any("lv" == type)) {
          # check if FLAT for random slopes
          # if (!is.null(partable$rv) && any(nchar(partable$rv) > 0L) &&
          #    !is.null(partable$block)) {
          #    return.value$lv[[b]] <- lv.names2
          # } else {
          # here, they will be 'defined' at level 2 as regular =~ lvs
          return.value$lv[[b]] <- lv.names
          # }
        }

        # regular latent variables ONLY (ie defined by =~ only)
        if (any("lv.regular" == type)) {
          out <- unique(partable$lhs[block.ind &
            partable$op == "=~" &
            !partable$lhs %in% rv.names])
          return.value$lv.regular[[b]] <- out
        }


        # interaction terms involving latent variables (only)
        if (any("lv.interaction" == type)) {
          return.value$lv.interaction[[b]] <- lv.interaction
        }

        # formative latent variables ONLY (ie defined by <~ only)
        if (any("lv.formative" == type)) {
          out <- unique(partable$lhs[block.ind &
            partable$op == "<~"])
          return.value$lv.formative[[b]] <- out
        }

        # lv's involved in efa
        if (any(type %in% c("lv.efa", "ov.efa"))) {
          if (is.null(partable$efa)) {
            out <- character(0L)
          } else {
            set.names <- lav_partable_efa_values(partable)
            out <- unique(partable$lhs[partable$op == "=~" &
              block.ind &
              partable$efa %in% set.names])
          }
          if (any(type == "ov.efa")) ov_efa <- out
          if (any(type == "lv.efa")) return.value$lv.efa[[b]] <- out
        }

        # lv's that are random slopes
        if (any("lv.rv" == type)) {
          if (is.null(partable$rv)) {
            out <- character(0L)
          } else {
            out <- unique(partable$lhs[partable$op == "=~" &
              block.ind &
              partable$lhs %in% rv.names])
          }
          return.value$lv.rv[[b]] <- out
        }

        # lv's that are indicators of a higher-order factor
        if (any("lv.ind" == type)) {
          out <- unique(partable$rhs[block.ind &
            partable$op == "=~" &
            partable$rhs %in% lv.names])
          return.value$lv.ind[[b]] <- out
        }

        # eqs.y
        eqs.y <- unique(partable$lhs[block.ind &
          partable$op == "~"])

        # store eqs.y
        if (any("eqs.y" == type)) {
          return.value$eqs.y[[b]] <- eqs.y
        }

        # eqs.x
        eqs.x <- unique(partable$rhs[block.ind &
          (partable$op == "~" |
            partable$op == "<~")])

        # store eqs.x
        if (any("eqs.x" == type)) {
          return.value$eqs.x[[b]] <- eqs.x
        }

        # v.ind -- indicators of latent variables
        v.ind <- unique(partable$rhs[block.ind &
          partable$op == "=~"])

        # ov.*
        # 1. indicators, which are not latent variables themselves
        ov.ind <- v.ind[!v.ind %in% lv.names2]
        # 2. dependent ov's
        ov.y <- eqs.y[!eqs.y %in% c(lv.names2, ov.ind)]
        # 3. independent ov's
        if (lav_partable_nlevels(partable) > 1L && b > 1L) {
          # NEW in 0.6-8: if an 'x' was an 'y' in a previous level,
          #               treat it as 'y'
          tmp.eqs.y <- unique(partable$lhs[partable$op == "~"]) # all blocks
          ov.x <- eqs.x[!eqs.x %in% c(lv.names2, ov.ind, tmp.eqs.y)]
        } else {
          ov.x <- eqs.x[!eqs.x %in% c(lv.names2, ov.ind, ov.y)]
        }
        # new in 0.6-12: if we have interaction terms in ov.x, check
        # if some terms are in eqs.y; if so, remove the interaction term
        # from ov.x
        int.idx <- which(grepl(":", ov.x, fixed = TRUE))
        bad.idx <- integer(0L)
        for (iv in int.idx) {
          tmp.names <- strsplit(ov.x[iv], ":", fixed = TRUE)[[1L]]
          if (any(tmp.names %in% eqs.y)) {
            bad.idx <- c(bad.idx, iv)
          }
        }
        if (length(bad.idx) > 0L) {
          ov.y <- unique(c(ov.y, ov.x[bad.idx]))
          # it may be removed later, but needed to construct ov.names
          ov.x <- ov.x[-bad.idx]
        }

        # observed variables
        # easy approach would be: everything that is not in lv.names,
        # but the main purpose here is to 'order' the observed variables
        # according to 'type' (indicators, ov.y, ov.x, orphans)

        # 4. orphaned covariances
        ov.cov <- c(
          partable$lhs[block.ind &
            partable$op == "~~" &
            !partable$lhs %in% lv.names2],
          partable$rhs[block.ind &
            partable$op == "~~" &
            !partable$rhs %in% lv.names2]
        )
        # 5. orphaned intercepts/thresholds
        ov.int <- partable$lhs[block.ind &
          (partable$op == "~1" |
            partable$op == "|") &
          !partable$lhs %in% lv.names2]

        ov.tmp <- c(ov.ind, ov.y, ov.x)
        ov.extra <- unique(c(ov.cov, ov.int)) # must be in this order!
        # so that
        # lav_partable_independence
        # retains the same order
        ov.names <- c(ov.tmp, ov.extra[!ov.extra %in% ov.tmp])

        # store ov?
        if (any("ov" == type)) {
          return.value$ov[[b]] <- ov.names
        }

        if (any("ov.ind" == type)) {
          return.value$ov.ind[[b]] <- ov.ind
        }

        if (any("ov.interaction" == type)) {
          ov.int.names <- ov.names[grepl(":", ov.names, fixed = TRUE)]
          n.int <- length(ov.int.names)
          if (n.int > 0L) {
            ov.names.noint <- ov.names[!ov.names %in% ov.int.names]

            ok.idx <- logical(n.int)
            for (iv in seq_len(n.int)) {
              tmp.names <- strsplit(ov.int.names[iv], ":", fixed = TRUE)[[1L]]

              # two scenario's:
              # - both variables are in ov.names.noint (ok)
              # - at least one variables is NOT in ov.names.noint (ignore)
              if (all(tmp.names %in% ov.names.noint)) {
                ok.idx[iv] <- TRUE
              }
            }
            ov.interaction <- ov.int.names[ok.idx]
          } else {
            ov.interaction <- character(0L)
          }

          return.value$ov.interaction[[b]] <- ov.interaction
        }

        if (any("ov.efa" == type)) {
          ov.efa <- partable$rhs[partable$op == "=~" &
            block.ind &
            partable$rhs %in% ov.ind &
            partable$lhs %in% ov_efa]
          return.value$ov.efa[[b]] <- unique(ov.efa)
        }


        # exogenous `x' covariates
        if (any(type %in% c(
          "ov.x", "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          # correction: is any of these ov.names.x mentioned as a variance,
          #             covariance, or intercept?
          # this should trigger a warning in lavaanify()
          if (is.null(partable$user)) { # FLAT!
            partable$user <- rep(1L, length(partable$lhs))
          }
          vars <- c(
            partable$lhs[block.ind &
              partable$op == "~1" &
              partable$user == 1],
            partable$lhs[block.ind &
              partable$op == "~~" &
              partable$user == 1],
            partable$rhs[block.ind &
              partable$op == "~~" &
              partable$user == 1]
          )
          idx.no.x <- which(ov.x %in% vars)
          if (length(idx.no.x)) {
            if (ov.x.fatal) {
              lav_msg_stop(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s];  Please remove them
                and try again.", lav_msg_view(ov.x[idx.no.x], "none")))
            }
            if (warn) {
              lav_msg_warn(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s]; these variables will
                now be treated as random introducing additional free parameters.
                If you wish to treat those variables as fixed, remove these
                formulas from the model syntax. Otherwise, consider adding the
                fixed.x = FALSE option.", lav_msg_view(ov.x[idx.no.x], "none")))
            }
            ov.x <- ov.x[-idx.no.x]
          }
          ov.tmp.x <- ov.x

          # extra
          if (!is.null(partable$exo)) {
            ov.cov <- c(
              partable$lhs[block.ind &
                partable$op == "~~" &
                partable$exo == 1L],
              partable$rhs[block.ind &
                partable$op == "~~" &
                partable$exo == 1L]
            )
            ov.int <- partable$lhs[block.ind &
              partable$op == "~1" &
              partable$exo == 1L]
            ov.extra <- unique(c(ov.cov, ov.int))
            ov.tmp.x <- c(ov.tmp.x, ov.extra[!ov.extra %in% ov.tmp.x])
          }

          ov.names.x <- ov.tmp.x
        }

        # store ov.x?
        if (any("ov.x" == type)) {
          return.value$ov.x[[b]] <- ov.names.x
        }

        # story ov.orphan?
        if (any("ov.orphan" == type)) {
          return.value$ov.orphan[[b]] <- ov.extra
        }

        # ov's withouth ov.x
        if (any(type %in% c(
          "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          ov.names.nox <- ov.names[!ov.names %in% ov.names.x]
        }

        # store ov.nox
        if (any("ov.nox" == type)) {
          return.value$ov.nox[[b]] <- ov.names.nox
        }

        # store ov.model
        if (any("ov.model" == type)) {
          # if no conditional.x, this is just ov
          # else, this is ov.nox
          if (any(block.ind & partable$op == "~" &
            partable$exo == 1L)) {
            return.value$ov.model[[b]] <- ov.names.nox
          } else {
            return.value$ov.model[[b]] <- ov.names
          }
        }

        # ov's strictly ordered
        if (any(type %in% c(
          "ov.ord", "th", "th.mean",
          "ov.num", "lv.nonnormal"
        ))) {
          tmp <- unique(partable$lhs[block.ind &
            partable$op == "|"])
          ord.names <- ov.names[ov.names %in% tmp]
        }

        if (any("ov.ord" == type)) {
          return.value$ov.ord[[b]] <- ord.names
        }

        # ov's strictly numeric
        if (any(type %in% c("ov.num", "lv.nonnormal"))) {
          ov.num <- ov.names[!ov.names %in% ord.names]
        }

        if (any("ov.num" == type)) {
          return.value$ov.num[[b]] <- ov.num
        }

        # nonnormal lv's
        if (any("lv.nonnormal" == type)) {
          # regular lv's
          lv.reg <- unique(partable$lhs[block.ind &
            partable$op == "=~"])
          if (length(lv.reg) > 0L) {
            out <- unlist(lapply(lv.reg, function(x) {
              # get indicators for this lv
              tmp.ind <- unique(partable$rhs[block.ind &
                partable$op == "=~" &
                partable$lhs == x])
              if (!all(tmp.ind %in% ov.num)) {
                return(x)
              } else {
                return(character(0))
              }
            }), use.names = FALSE)
            return.value$lv.nonnormal[[b]] <- out
          } else {
            return.value$lv.nonnormal[[b]] <- character(0)
          }
        }

        if (any(c("th", "th.mean") %in% type)) {
          tmp.th.lhs <- partable$lhs[block.ind &
            partable$op == "|"]
          tmp.th.rhs <- partable$rhs[block.ind &
            partable$op == "|"]
        }

        # threshold
        if (any("th" == type)) {
          if (length(ord.names) > 0L) {
            # return in the right order (following ord.names!)
            out <- unlist(lapply(ord.names, function(x) {
              idx <- which(x == tmp.th.lhs)
              tmp.th <- unique(paste(tmp.th.lhs[idx], "|",
                tmp.th.rhs[idx],
                sep = ""
              ))
              # make sure the th's are in increasing order
              # sort(tmp.th)
              # NO!, don't do that; t10 will be before t2
              # fixed in 0.6-1 (bug report from Myrsini)
              # in 0.6-12, we do this anyway like this:

              # get var name
              tmp.th1 <- sapply(
                strsplit(tmp.th, split = "\\|t"),
                "[[", 1
              )
              # get number, and sort
              tmp.th2 <- as.character(sort(as.integer(sapply(
                strsplit(tmp.th, split = "\\|t"), "[[", 2
              ))))
              # paste back togehter in the right order
              paste(tmp.th1, tmp.th2, sep = "|t")
            }), use.names = FALSE)
          } else {
            out <- character(0L)
          }
          return.value$th[[b]] <- out
        }

        # thresholds and mean/intercepts of numeric variables
        if (any("th.mean" == type)) {
          # if fixed.x -> use ov.names.nox
          # else -> use ov.names
          if (is.null(partable$exo) || all(partable$exo == 0L)) {
            tmp.ov.names <- ov.names
          } else {
            tmp.ov.names <- ov.names.nox
          }
          if (length(tmp.ov.names) > 0L) {
            # return in the right order (following ov.names.nox!)
            out <- unlist(lapply(tmp.ov.names, function(x) {
              if (x %in% ord.names) {
                idx <- which(x == tmp.th.lhs)
                tmp.th <- unique(paste(tmp.th.lhs[idx], "|",
                  tmp.th.rhs[idx],
                  sep = ""
                ))
                # make sure the th's are in increasing order
                # get var name
                tmp.th1 <- sapply(
                  strsplit(tmp.th, split = "\\|t"),
                  "[[", 1
                )
                # get number, and sort
                tmp.th2 <- as.character(sort(as.integer(sapply(
                  strsplit(tmp.th, split = "\\|t"), "[[", 2
                ))))
                # paste back togehter in the right order
                paste(tmp.th1, tmp.th2, sep = "|t")
              } else {
                x
              }
            }), use.names = FALSE)
          } else {
            out <- character(0L)
          }
          return.value$th.mean[[b]] <- out
        }

        # exogenous lv's
        if (any(c("lv.x", "lv.nox") %in% type)) {
          tmp <- lv.names[!lv.names %in% c(v.ind, eqs.y)]
          lv.names.x <- lv.names[lv.names %in% tmp]
        }

        if (any("lv.x" == type)) {
          return.value$lv.x[[b]] <- lv.names.x
        }

        # dependent ov (but not also indicator or x)
        if (any("ov.y" == type)) {
          tmp <- eqs.y[!eqs.y %in% c(v.ind, eqs.x, lv.names)]
          return.value$ov.y[[b]] <- ov.names[ov.names %in% tmp]
        }

        # dependent lv (but not also indicator or x)
        if (any("lv.y" == type)) {
          tmp <- eqs.y[!eqs.y %in% c(v.ind, eqs.x) &
            eqs.y %in% lv.names]
          return.value$lv.y[[b]] <- lv.names[lv.names %in% tmp]
        }

        # non-exogenous latent variables
        if (any("lv.nox" == type)) {
          return.value$lv.nox[[b]] <- lv.names[!lv.names %in% lv.names.x]
        }

        # marker indicator (if any) for each lv
        if (any("lv.marker" == type)) {
          # default: "" per lv
          out <- character(length(lv.names))
          names(out) <- lv.names
          for (l in seq_len(length(lv.names))) {
            this.lv.name <- lv.names[l]
            # try to see if we can find a 'marker' indicator for this factor
            marker.idx <- which(block.ind &
              partable$op == "=~" &
              partable$lhs == this.lv.name &
              partable$rhs %in% v.ind &
              partable$ustart == 1L &
              partable$free == 0L)
            if (length(marker.idx) == 1L) { # unique only!!
              # check if 'other' loadings are fixed to zero
              other.idx <- which(block.ind &
                partable$op == "=~" &
                partable$lhs != this.lv.name &
                partable$rhs == partable$rhs[marker.idx] &
                partable$free == 0L)
              if (length(other.idx) == 0L) {
                out[l] <- partable$rhs[marker.idx]
              } else if (all(partable$ustart[other.idx] == 0)) {
                out[l] <- partable$rhs[marker.idx]
              }
            }
          }
          return.value$lv.marker[[b]] <- out
        }
      }
    } # b

    # ----- lav_partable_vnames ---- no cache ------------------------
    # new in 0.6-14: if 'da' operator, change order! (for ov.order = "data")
    # now via attribute "ovda"
    ov.names.data <- attr(partable, "ovda")
    if (!is.null(ov.names.data)) {
      return.value <- lapply(return.value, function(x) {
        for (b in seq_len(length(x))) {
          m <- match(x[[b]], ov.names.data)
          target.idx <- which(!is.na(m))
          if (length(target.idx) > 1L) {
            x[[b]][target.idx] <- x[[b]][target.idx][order(m[target.idx])]
          }
        }
        x
      })
    }
  }
  # ----- lav_partable_vnames ---- common ------ 1 type --------
  # to mimic old behaviour, if length(type) == 1L
  if (length(type) == 1L) {
    return.value <- return.value[[type]]
    # to mimic old behaviour, if specific block is requested
    if (ndotdotdot == 0L) {
      if (type == "lv.marker") {
        return.value <- unlist(return.value)
        # no unique(), as unique() drops attributes, and reduces
        # c("", "", "") to a single ""
        # (but, say for 2 groups, you get 2 copies)
        # as this is only for 'display', we leave it like that
      } else {
        return.value <- unique(unlist(return.value))
      }
    } else if (length(block.select) == 1L) {
      return.value <- return.value[[block.select]]
    } else {
      return.value <- return.value[block.select]
    }
  }
  # ----- lav_partable_vnames ---- common ------------------------
  return.value
}

# alias for backward compatibility
vnames <- lav_partable_vnames
