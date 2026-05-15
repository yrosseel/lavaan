# lav_partable_names
#
# YR. 29 june 2013
#  - as separate file; used to be in utils-user.R
#  - lav_partable_names (aka 'vnames') allows multiple options in 'type'
#    returning them all as a list (or just a vector if only 1 type is needed)

# public version
lav_object_vnames <- function(object, type = "ov", ...) { # nolint
  if (inherits(object, "lavaan") || inherits(object, "lavaanList")) {
    # check object
    object <- lav_object_check_version(object)
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
lavNames <- lav_object_vnames    # synonym #nolint

lav_partable_vnames_cached_block_select <- function(partable, dotdotdot) {
  dot_names <- names(dotdotdot)
  if (is.null(dot_names) ||
      any(!nzchar(dot_names)) ||
      !all(dot_names %in% c("block", "group", "level"))) {
    return(NULL)
  }

  if (is.null(partable$block)) {
    partable_block <- rep(1L, length(partable$lhs))
  } else {
    partable_block <- partable$block
  }
  valid_row <- partable_block > 0L &
    !partable$op %in% c("==", "<", ">", ":=")
  block_select <- unique(na.omit(partable_block[valid_row]))
  if (length(block_select) == 0L) {
    lav_msg_warn(gettext("no blocks selected."))
    return(block_select)
  }

  block_row <- match(block_select, partable_block)
  for (dot in seq_along(dotdotdot)) {
    block_var <- dot_names[dot]
    block_val <- dotdotdot[[dot]]

    if (block_var == "block") {
      if (!all(block_val %in% partable_block)) {
        lav_msg_stop(gettextf(
          "%1$s column does not contain value `%2$s'", block_var, block_val))
      }
      block_select <- block_select[block_select %in% block_val]
      block_row <- match(block_select, partable_block)
      next
    }

    block_var_values <- partable[[block_var]]
    if (is.null(block_var_values) || length(block_var_values) == 0L) {
      if (block_var == "group" &&
          length(block_val) == 1L &&
          !is.na(block_val) &&
          block_val == 1L) {
        next
      }
      return(NULL)
    }
    if (!all(block_val %in% block_var_values)) {
      lav_msg_stop(gettextf(
        "%1$s column does not contain value `%2$s'", block_var, block_val))
    }

    block_select <- block_select[block_var_values[block_row] %in% block_val]
    block_row <- match(block_select, partable_block)
  }

  if (length(block_select) == 0L) {
    lav_msg_warn(gettext("no blocks selected."))
  }
  block_select
}

lav_partable_vnames_cached_block_subset <- function(cached_values,
                                                    block_select) {
  if (length(block_select) == 0L) {
    return(vector("list", length = 0L))
  }

  out <- vector("list", length = max(block_select))
  out[block_select] <- cached_values[block_select]
  out
}

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
# - YR 11/02/25: add type = "lv.ho" for higher-order latent variables
lav_partable_vnames <- function(partable, type = NULL, ..., # nolint
                                force_warn = FALSE, ov_x_fatal = FALSE) {
  # This function derives the names of some types of variable (as specified
  # in type) from a 'partable'. The 'warn' parameter needs no explanation.
  # The ov.x.fatal parameter implies, when set to TRUE, that the function
  # issues a 'stop' when there are exogenous variables present in variance/
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
  # this is a special function where the default is to suppress warnings,
  # overwritten if parameter force.warn TRUE (used in lav_partable function)
  current_warn <- lav_warn()
  if (force_warn) {
    if (lav_warn(TRUE)) on.exit(lav_warn(current_warn))
  } else {
    if (lav_warn(FALSE)) on.exit(lav_warn(current_warn))
  }
  # check for empty table
  if (length(partable$lhs) == 0) {
    return(character(0L))
  }
  # dotdotdot
  dotdotdot <- list(...)
  ndotdotdot <- length(dotdotdot)
  type_list <- c(
    "ov", # observed variables (ov)
    "ov.x", # (pure) exogenous observed variables
    "ov.nox", # non-exogenous observed variables
    "ov.model", # modeled observed variables (joint vs cond)
    "ov.y", # (pure) endogenous variables (dependent only)
    "ov.num", # numeric observed variables
    "ov.ord", # ordinal observed variables
    "ov.ind",  # observed indicators of latent variables
    "ov.cind", # observed indicators of composites (new in 0.6-20)
    "ov.orphan", # lonely observed intercepts/variances
    "ov.interaction", # interaction terms (with colon)
    "ov.efa", # indicators involved in efa

    "th", # thresholds ordinal only
    "th.mean", # thresholds ordinal + numeric variables

    "lv", # latent variables
    "lv.regular", # latent variables (defined by =~ only)
    "lv.formative", # latent variables (defined by <~ only) (old style)
    "lv.composite", # latent variables (defined by <~ only) (new style)
    "lv.x", # (pure) exogenous variables
    "lv.y", # (pure) endogenous variables
    "lv.nox", # non-exogenous latent variables
    "lv.nonnormal", # latent variables with non-normal indicators
    "lv.interaction", # interaction terms
    "lv.efa", # latent variables involved in efa
    "lv.rv", # random slopes, random variables
    "lv.ind", # latent indicators (higher-order cfa)
    "lv.ho",  # higher-order latent variables
    "lv.marker", # marker indicator per lv

    "eqs.y", # y's in regression
    "eqs.x" # x's in regression
  )
  if (type[1L] != "all" && type[1L] != "*" && !all(type %in% type_list)) {
    wrongtypes <- type[!(type %in% type_list)]
    lav_msg_stop(sprintf(
      ngettext(length(wrongtypes),
      "type = %s is not a valid option",
      "type = %s are not valid options"),
      lav_msg_view(wrongtypes, "none", FALSE)))
  }
  return_value <- NULL
  if (type[1L] != "*" && !is.null(attr(partable, "vnames"))) {
    # ----- lav_partable_vnames ---- cached data --------------------------
    # uncomment/comment following line to enable/disable trace
    # lav_trace(paste("cached:", paste(type, collapse = ",")))

    if (type[1L] == "all") {
      return_value <- attr(partable, "vnames")
    } else {
      return_value <- attr(partable, "vnames")[type]
    }
    if (ndotdotdot == 0L) {
      if (type[1L] == "all") {
        return(return_value)
      } else if (length(type) == 1L) {
        if (type == "lv.marker") {
          return(unlist(return_value[[type]]))
        } else {
          return(unique(unlist(return_value[[type]])))
        }
      } else {
        return(return_value)
      }
    } else {
      block_select <- lav_partable_vnames_cached_block_select(
        partable,
        dotdotdot
      )
      if (!is.null(block_select)) {
        if (type[1L] == "all") {
          return(lapply(return_value,
            lav_partable_vnames_cached_block_subset,
            block_select = block_select
          ))
        } else if (length(type) == 1L) {
          return_value <- return_value[[type]]
          if (length(block_select) == 1L) {
            return(return_value[[block_select]])
          } else {
            return(lav_partable_vnames_cached_block_subset(
              return_value,
              block_select
            ))
          }
        } else {
          return(lapply(return_value,
            lav_partable_vnames_cached_block_subset,
            block_select = block_select
          ))
        }
      }
    }
  }
  # ----- lav_partable_vnames ---- common ----------------------------------
  if (type[1L] == "all" || type[1L] == "*") {
    type <- type_list
  }
  # ALWAYS need `block' column -- create one if missing
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$lhs))
  }
  # per default, use full partable
  block_select <- lav_partable_block_values(partable)
  # check for ... selection argument(s)
  if (ndotdotdot > 0L) {
    dot_names <- names(dotdotdot)
    row_select <- rep(TRUE, length(partable$lhs))
    for (dot in seq_len(ndotdotdot)) {
      # selection variable?
      block_var <- dot_names[dot]
      block_val <- dotdotdot[[block_var]]
      # do we have this 'block.var' in partable?
      if (is.null(partable[[block_var]])) {
        # for historical reasons, treat "group = 1" special
        if (block_var == "group" && block_val == 1L) {
          partable$group <- rep(1L, length(partable$lhs))
          # remove block == 0
          idx <- which(partable$block == 0L)
          if (length(idx) > 0L) {
            partable$group[idx] <- 0L
          }
          row_select <- (row_select &
            partable[[block_var]] %in% block_val)
        } else {
          lav_msg_stop(gettextf(
            "selection variable '%s' not found in the parameter table.",
            block_var))
        }
      } else {
        if (!all(block_val %in% partable[[block_var]])) {
          lav_msg_stop(gettextf(
            "%1$s column does not contain value `%2$s'", block_var, block_val))
        }
        row_select <- (row_select &
          !partable$op %in% c("==", "<", ">", ":=") &
          partable[[block_var]] %in% block_val)
      }
    } # dot
    block_select <- unique(partable$block[row_select])
    if (length(block_select) == 0L) {
      lav_msg_warn(gettext("no blocks selected."))
    }
  }
  if (is.null(return_value)) {
    # ----- lav_partable_vnames ---- no cache ------------------------

    # uncomment/comment following line to enable/disable trace
    # lav_trace(paste("computed:", paste(type, collapse = ",")))

    # random slope names, if any (new in 0.6-7)
    if (!is.null(partable$rv) && any(nchar(partable$rv) > 0L)) {
      rv_names <- unique(partable$rv[nchar(partable$rv) > 0L])
    } else {
      rv_names <- character(0L)
    }

    # output: list per block
    return_value <- lapply(type, function(x) {
      vector("list", length = length(block_select))
    })
    names(return_value) <- type

    for (b in block_select) {
      # indices for this block
      block_ind <- partable$block == b

      # always compute lv.names
      lv_names <- unique(partable$lhs[block_ind &
        (partable$op == "=~" | partable$op == "<~")])
      # including random slope names
      lv_names2 <- unique(c(lv_names, rv_names))

      # determine lv interactions
      int_names <- unique(partable$rhs[block_ind &
        grepl(":", partable$rhs, fixed = TRUE)])
      n_int <- length(int_names)
      if (n_int > 0L) {
        ok_idx <- logical(n_int)
        for (iv in seq_len(n_int)) {
          tmp_names <- strsplit(int_names[iv], ":", fixed = TRUE)[[1L]]
          # three scenario's:
          # - both variables are latent (ok)
          # - both variables are observed (ignore)
          # - only one latent (warn??) -> upgrade observed to latent
          # thus if at least one is in lv.names, we treat it as a
          # latent interaction
          if (any(tmp_names %in% lv_names)) {
            ok_idx[iv] <- TRUE
          }
        }
        lv_interaction <- int_names[ok_idx]
        lv_names <- c(lv_names, lv_interaction)
        lv_names2 <- c(lv_names2, lv_interaction)
      } else {
        lv_interaction <- character(0L)
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
          return_value$lv[[b]] <- lv_names
          # }
          next
        }

        # regular latent variables ONLY (ie defined by =~ only)
        if ("lv.regular" == type) {
          out <- unique(partable$lhs[block_ind &
            partable$op == "=~" &
            partable$lhs != partable$rhs & # no phantom
            !partable$lhs %in% rv_names])
          return_value$lv.regular[[b]] <- out
          next
        }

        # interaction terms involving latent variables (only)
        if ("lv.interaction" == type) {
          return_value$lv.interaction[[b]] <- lv_interaction
          next
        }

        # formative latent variables: phantom lv + zero residual
        if ("lv.formative" == type) {
          out <- unique(partable$lhs[block_ind &
            partable$lhs == partable$rhs &
            partable$op == "=~"])
          rm_idx <- integer(0L)
          for (i in seq_along(out)) {
            var_idx <- which(block_ind &
                             partable$op == "~~" &
                             partable$lhs == out[i] &
                             partable$lhs == partable$rhs)
            if (length(var_idx) == 0L) {
              next
            }
            if (!is.null(partable$mod.idx) && !is.null(partable$fixed)) {
              if (partable$fixed[var_idx] != "0") {
                rm_idx <- c(rm_idx, i)
              }
            } else if (!is.null(partable$free) && !is.null(partable$ustart)) {
              if (partable$free[var_idx] > 0L ||
                  partable$ustart[var_idx] != 0) {
                rm_idx <- c(rm_idx, i)
              }
            }
          }
          if (length(rm_idx) > 0L) {
            out <- out[-rm_idx]
          }
          return_value$lv.formative[[b]] <- out
          next
        }

        # composites defined by "<~"
        if ("lv.composite" == type) {
          out <- unique(partable$lhs[block_ind &
            partable$op == "<~"])
          return_value$lv.composite[[b]] <- out
          next
        }

        # lv's involved in efa
        if (any(type == c("lv.efa", "ov.efa"))) {
          if (is.null(partable$efa)) {
            out <- character(0L)
          } else {
            set_names <- lav_partable_efa_values(partable)
            out <- unique(partable$lhs[partable$op == "=~" &
              block_ind &
              partable$efa %in% set_names])
          }
          if (type == "ov.efa") ov_efa <- out
          if (type == "lv.efa") {
            return_value$lv.efa[[b]] <- out
            next
          }
        }

        # lv's that are random slopes
        if ("lv.rv" == type) {
          if (is.null(partable$rv)) {
            out <- character(0L)
          } else {
            out <- unique(partable$lhs[partable$op == "=~" &
              block_ind &
              partable$lhs %in% rv_names])
          }
          return_value$lv.rv[[b]] <- out
          next
        }

        # lv's that are indicators of a higher-order factor
        if ("lv.ind" == type) {
          out <- unique(partable$rhs[block_ind &
            partable$op == "=~" &
            partable$rhs %in% lv_names &
            partable$lhs != partable$rhs]) # no phantom lv's
          return_value$lv.ind[[b]] <- out
          next
        }

        # higher-order latent variables
        if ("lv.ho" == type) {
          out_ind <- unique(partable$rhs[block_ind &
            partable$op == "=~" &
            partable$rhs %in% lv_names &
            partable$lhs != partable$rhs]) # no phantom lv's
          out <- unique(partable$lhs[block_ind &
            partable$op == "=~" &
            partable$lhs %in% lv_names &
            partable$rhs %in% out_ind])
          return_value$lv.ho[[b]] <- out
          next
        }

        # eqs.y
        if (!(any(type == c("lv", "lv.regular")))) {
          eqs_y <- unique(partable$lhs[block_ind &
            partable$op == "~"])
        }

        # store eqs.y
        if ("eqs.y" == type) {
          return_value$eqs.y[[b]] <- eqs_y
          next
        }

        # eqs.x
        if (!(any(type == c("lv", "lv.regular", "lv.x")))) {
          eqs_x <- unique(partable$rhs[block_ind &
            (partable$op == "~" |
              partable$op == "<~")])
        }

        # store eqs.x
        if ("eqs.x" == type) {
          return_value$eqs.x[[b]] <- eqs_x
          next
        }

        # v.ind -- indicators of latent variables
        if (!(any(type == c("lv", "lv.regular")))) {
          v_ind <- unique(partable$rhs[block_ind &
            partable$op == "=~"])
        }

        # v.cind -- indicators of composites
        if (!(any(type == c("lv", "lv.regular")))) {
          v_cind <- unique(partable$rhs[block_ind &
            partable$op == "<~"])
        }

        # ov.*
        if (!(any(type == c("lv", "lv.regular", "lv.x", "lv.y")))) {
          # 1. indicators, which are not latent variables themselves
          ov_ind <- v_ind[!v_ind %in% lv_names2]
          # 1b. indicator of composites
          ov_cind <- v_cind[!v_cind %in% lv_names2]
          # 2. dependent ov's
          ov_y <- eqs_y[!eqs_y %in% c(lv_names2, ov_ind, ov_cind)]
          # 3. independent ov's
          if (lav_partable_nlevels(partable) > 1L && b > 1L) {
            # NEW in 0.6-8: if an 'x' was an 'y' in a previous level,
            #               treat it as 'y'
            tmp_eqs_y <- unique(partable$lhs[partable$op == "~"]) # all blocks
            ov_x <- eqs_x[!eqs_x %in% c(lv_names2, ov_ind, ov_cind, tmp_eqs_y)]
          } else {
            ov_x <- eqs_x[!eqs_x %in% c(lv_names2, ov_ind, ov_cind, ov_y)]
          }
          # new in 0.6-12: if we have interaction terms in ov.x, check
          # if some terms are in eqs.y; if so, remove the interaction term
          # from ov.x
          int_idx <- which(grepl(":", ov_x, fixed = TRUE))
          bad_idx <- integer(0L)
          for (iv in int_idx) {
            tmp_names <- strsplit(ov_x[iv], ":", fixed = TRUE)[[1L]]
            if (any(tmp_names %in% eqs_y)) {
              bad_idx <- c(bad_idx, iv)
            }
          }
          if (length(bad_idx) > 0L) {
            ov_y <- unique(c(ov_y, ov_x[bad_idx]))
            # it may be removed later, but needed to construct ov.names
            ov_x <- ov_x[-bad_idx]
          }
        }

        # observed variables
        # easy approach would be: everything that is not in lv.names,
        # but the main purpose here is to 'order' the observed variables
        # according to 'type' (indicators, ov.y, ov.x, orphans)
        if (!(any(type == c("lv", "lv.regular", "lv.x", "lv.y")))) {
          # 4. orphaned covariances
          ov_cov <- c(
            partable$lhs[block_ind &
              partable$op == "~~" &
              !partable$lhs %in% lv_names2],
            partable$rhs[block_ind &
              partable$op == "~~" &
              !partable$rhs %in% lv_names2]
          )
          # 5. orphaned intercepts/thresholds
          ov_int <- partable$lhs[block_ind &
            (partable$op == "~1" |
              partable$op == "|") &
            !partable$lhs %in% lv_names2]

          ov_tmp <- c(ov_ind, ov_cind, ov_y, ov_x)
          ov_extra <- unique(c(ov_cov, ov_int)) # must be in this order!
          # so that
          # lav_partable_independence
          # retains the same order
          ov_names <- c(ov_tmp, ov_extra[!ov_extra %in% ov_tmp])
        }

        # store ov?
        if ("ov" == type) {
          return_value$ov[[b]] <- ov_names
          next
        }

        if ("ov.ind" == type) {
          return_value$ov.ind[[b]] <- ov_ind
          next
        }

        if ("ov.cind" == type) {
          return_value$ov.cind[[b]] <- ov_cind
          next
        }

        if ("ov.interaction" == type) {
          ov_int_names <- ov_names[grepl(":", ov_names, fixed = TRUE)]
          n_int <- length(ov_int_names)
          if (n_int > 0L) {
            ov_names_noint <- ov_names[!ov_names %in% ov_int_names]

            ok_idx <- logical(n_int)
            for (iv in seq_len(n_int)) {
              tmp_names <- strsplit(ov_int_names[iv], ":", fixed = TRUE)[[1L]]

              # two scenario's:
              # - both variables are in ov.names.noint (ok)
              # - at least one variable is NOT in ov.names.noint (ignore)
              if (all(tmp_names %in% ov_names_noint)) {
                ok_idx[iv] <- TRUE
              }
            }
            ov_interaction <- ov_int_names[ok_idx]
          } else {
            ov_interaction <- character(0L)
          }

          return_value$ov.interaction[[b]] <- ov_interaction
          next
        }

        if ("ov.efa" == type) {
          ov_efa_1 <- partable$rhs[partable$op == "=~" &
            block_ind &
            partable$rhs %in% ov_ind &
            partable$lhs %in% ov_efa]
          return_value$ov.efa[[b]] <- unique(ov_efa_1)
          next
        }


        # exogenous `x' covariates
        if (any(type == c(
          "ov.x", "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          # correction: is any of these ov.names.x mentioned as a variance,
          #             covariance, or intercept?
          # this should trigger a warning in lav_model_partable()
          if (is.null(partable$user)) { # FLAT!
            partable$user <- rep(1L, length(partable$lhs))
          }
          vars <- c(
            partable$lhs[block_ind &
              partable$op == "~1" &
              partable$user == 1],
            partable$lhs[block_ind &
              partable$op == "~~" &
              partable$user == 1],
            partable$rhs[block_ind &
              partable$op == "~~" &
              partable$user == 1]
          )
          idx_no_x <- which(ov_x %in% vars)
          if (length(idx_no_x)) {
            if (ov_x_fatal) {
              lav_msg_stop(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s];  Please remove them
                and try again.", lav_msg_view(ov_x[idx_no_x], "none")))
            }
            lav_msg_warn(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s]; these variables will
                now be treated as random introducing additional free parameters.
                If you wish to treat those variables as fixed, remove these
                formulas from the model syntax. Otherwise, consider adding the
                fixed.x = FALSE option.", lav_msg_view(ov_x[idx_no_x], "none")))
            ov_x <- ov_x[-idx_no_x]
          }
          ov_tmp_x <- ov_x

          # extra
          if (!is.null(partable$exo)) {
            ov_cov <- c(
              partable$lhs[block_ind &
                partable$op == "~~" &
                partable$exo == 1L],
              partable$rhs[block_ind &
                partable$op == "~~" &
                partable$exo == 1L]
            )
            ov_int <- partable$lhs[block_ind &
              partable$op == "~1" &
              partable$exo == 1L]
            ov_extra <- unique(c(ov_cov, ov_int))
            ov_tmp_x <- c(ov_tmp_x, ov_extra[!ov_extra %in% ov_tmp_x])
          }

          ov_names_x <- ov_tmp_x
        }

        # store ov.x?
        if ("ov.x" == type) {
          return_value$ov.x[[b]] <- ov_names_x
          next
        }

        # story ov.orphan?
        if ("ov.orphan" == type) {
          return_value$ov.orphan[[b]] <- ov_extra
          next
        }

        # ov's without ov.x
        if (any(type == c(
          "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          ov_names_nox <- ov_names[!ov_names %in% ov_names_x]
        }

        # store ov.nox
        if ("ov.nox" == type) {
          return_value$ov.nox[[b]] <- ov_names_nox
          next
        }

        # store ov.model
        if ("ov.model" == type) {
          # if no conditional.x, this is just ov
          # else, this is ov.nox
          if (any(block_ind & partable$op == "~" &
            partable$exo == 1L)) {
            return_value$ov.model[[b]] <- ov_names_nox
          } else {
            return_value$ov.model[[b]] <- ov_names
          }
          next
        }

        # ov's strictly ordered
        if (any(type == c(
          "ov.ord", "th", "th.mean",
          "ov.num", "lv.nonnormal"
        ))) {
          tmp <- unique(partable$lhs[block_ind &
            partable$op == "|"])
          ord_names <- ov_names[ov_names %in% tmp]
        }

        if ("ov.ord" == type) {
          return_value$ov.ord[[b]] <- ord_names
          next
        }

        # ov's strictly numeric
        if (any(type == c("ov.num", "lv.nonnormal"))) {
          ov_num <- ov_names[!ov_names %in% ord_names]
        }

        if ("ov.num" == type) {
          return_value$ov.num[[b]] <- ov_num
          next
        }

        # nonnormal lv's
        if ("lv.nonnormal" == type) {
          # regular lv's
          lv_reg <- unique(partable$lhs[block_ind &
                           partable$op == "=~" &
                           partable$lhs != partable$rhs])
          if (length(lv_reg) > 0L) {
            out <- unlist(lapply(lv_reg, function(x) {
              # get indicators for this lv
              tmp_ind <- unique(partable$rhs[block_ind &
                partable$op == "=~" &
                partable$lhs == x])
              if (!all(tmp_ind %in% ov_num)) {
                x
              } else {
                character(0)
              }
            }), use.names = FALSE)
            return_value$lv.nonnormal[[b]] <- out
          } else {
            return_value$lv.nonnormal[[b]] <- character(0)
          }
          next
        }

        if (any(c("th", "th.mean") == type)) {
          tmp_th_lhs <- partable$lhs[block_ind &
            partable$op == "|"]
          tmp_th_rhs <- partable$rhs[block_ind &
            partable$op == "|"]
        }

        # threshold
        if ("th" == type) {
          if (length(ord_names) > 0L) {
            # return in the right order (following ord.names!)
            out <- unlist(lapply(ord_names, function(x) {
              idx <- which(x == tmp_th_lhs)
              tmp_th <- unique(paste(tmp_th_lhs[idx], "|",
                tmp_th_rhs[idx],
                sep = ""
              ))
              # make sure the th's are in increasing order
              # sort(tmp.th)
              # NO!, don't do that; t10 will be before t2
              # fixed in 0.6-1 (bug report from Myrsini)
              # in 0.6-12, we do this anyway like this:

              # get var name
              tmp_th1 <- sapply(
                strsplit(tmp_th, split = "\\|t"),
                "[[", 1
              )
              # get number, and sort
              tmp_th2 <- as.character(sort(as.integer(sapply(
                strsplit(tmp_th, split = "\\|t"), "[[", 2
              ))))
              # paste back together in the right order
              paste(tmp_th1, tmp_th2, sep = "|t")
            }), use.names = FALSE)
          } else {
            out <- character(0L)
          }
          return_value$th[[b]] <- out
          next
        }

        # thresholds and mean/intercepts of numeric variables
        if ("th.mean" == type) {
          # if fixed.x -> use ov.names.nox
          # else -> use ov.names
          if (is.null(partable$exo) || all(partable$exo == 0L)) {
            tmp_ov_names <- ov_names
          } else {
            tmp_ov_names <- ov_names_nox
          }
          if (length(tmp_ov_names) > 0L) {
            # return in the right order (following ov.names.nox!)
            out <- unlist(lapply(tmp_ov_names, function(x) {
              if (x %in% ord_names) {
                idx <- which(x == tmp_th_lhs)
                tmp_th <- unique(
                  paste(tmp_th_lhs[idx], "|",
                    tmp_th_rhs[idx],
                    sep = ""
                  ),
                  use.names = FALSE
                )
                # make sure the th's are in increasing order
                # get var name
                tmp_th1 <- sapply(
                  strsplit(tmp_th, split = "\\|t"),
                  "[[", 1
                )
                # get number, and sort
                tmp_th2 <- as.character(sort(as.integer(sapply(
                  strsplit(tmp_th, split = "\\|t"), "[[", 2
                ))))
                # paste back together in the right order
                paste(tmp_th1, tmp_th2, sep = "|t")
              } else {
                x
              }
            }))
          } else {
            out <- character(0L)
          }
          return_value$th.mean[[b]] <- out
          next
        }

        # exogenous lv's
        if (any(c("lv.x", "lv.nox") == type)) {
          tmp <- lv_names[!lv_names %in% c(v_ind, eqs_y)]
          lv_names_x <- lv_names[lv_names %in% tmp]
        }

        if ("lv.x" == type) {
          return_value$lv.x[[b]] <- lv_names_x
          next
        }

        # dependent ov (but not also indicator or x)
        if ("ov.y" == type) {
          tmp <- eqs_y[!eqs_y %in% c(v_ind, eqs_x, lv_names)]
          return_value$ov.y[[b]] <- ov_names[ov_names %in% tmp]
          next
        }

        # dependent lv (but not also indicator or x)
        if ("lv.y" == type) {
          tmp <- eqs_y[!eqs_y %in% c(v_ind, eqs_x) &
            eqs_y %in% lv_names]
          return_value$lv.y[[b]] <- lv_names[lv_names %in% tmp]
          next
        }

        # non-exogenous latent variables
        if ("lv.nox" == type) {
          return_value$lv.nox[[b]] <- lv_names[!lv_names %in% lv_names_x]
          next
        }

        # marker indicator (if any) for each lv
        if ("lv.marker" == type) {
          # default: "" per lv
          out <- character(length(lv_names))
          names(out) <- lv_names
          for (l in seq_along(lv_names)) {
            this_lv_name <- lv_names[l]
            # try to see if we can find a 'marker' indicator for this factor
            # here defined as: has fixed-to-one factor loading
            marker_idx <- which(block_ind &
              partable$op == "=~" &
              partable$lhs == this_lv_name &
              partable$rhs %in% v_ind &
              partable$ustart == 1L &
              partable$free == 0L)
            if (length(marker_idx) > 0L) {
              # we may have multiple potential markers
              potential_markers <- partable$rhs[marker_idx]
              valid_markers <- character(0L)
              for (m in seq_along(potential_markers)) {
                this_marker_idx <- marker_idx[m]
                # check if 'other' loadings are fixed to zero
                other_idx <- which(block_ind &
                  partable$op == "=~" &
                  partable$lhs != this_lv_name &
                  partable$rhs == partable$rhs[this_marker_idx] &
                  partable$free == 0L)
                if (length(other_idx) == 0L) {
                  # simple structure, or one factor
                  valid_markers <- c(valid_markers,
                                     partable$rhs[this_marker_idx])
                } else if (all(partable$ustart[other_idx] == 0)) {
                  valid_markers <- c(valid_markers,
                                     partable$rhs[this_marker_idx])
                }
              }
              if (length(valid_markers) > 0L) {
                out[l] <- valid_markers[1L] # pick the first one
              }
            }
          } # l
          return_value$lv.marker[[b]] <- out
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
          return_value$lv[[b]] <- lv_names
          # }
        }

        # regular latent variables ONLY (ie defined by =~ only)
        if (any("lv.regular" == type)) {
          out <- unique(partable$lhs[block_ind &
            partable$op == "=~" &
            partable$lhs != partable$rhs & # no phantom
            !partable$lhs %in% rv_names])
          return_value$lv.regular[[b]] <- out
        }


        # interaction terms involving latent variables (only)
        if (any("lv.interaction" == type)) {
          return_value$lv.interaction[[b]] <- lv_interaction
        }

        # formative latent variables: phantom lv + zero residual
        if (any("lv.formative" == type)) {
          out <- unique(partable$lhs[block_ind &
            partable$lhs == partable$rhs &
            partable$op == "=~"])
          rm_idx <- integer(0L)
          for (i in seq_along(out)) {
            var_idx <- which(block_ind &
                             partable$op == "~~" &
                             partable$lhs == out[i] &
                             partable$lhs == partable$rhs)
            if (length(var_idx) == 0L) {
              next
            }
            if (!is.null(partable$mod.idx) && !is.null(partable$fixed)) {
              if (partable$fixed[var_idx] != "0") {
                rm_idx <- c(rm_idx, i)
              }
            } else if (!is.null(partable$free) && !is.null(partable$ustart)) {
              if (partable$free[var_idx] > 0L ||
                  partable$ustart[var_idx] != 0) {
                rm_idx <- c(rm_idx, i)
              }
            }
          }
          if (length(rm_idx) > 0L) {
            out <- out[-rm_idx]
          }
          return_value$lv.formative[[b]] <- out
        }

        # composite latent variables (ie defined by <~)
        if (any("lv.composite" == type)) {
          out <- unique(partable$lhs[block_ind &
            partable$op == "<~"])
          return_value$lv.composite[[b]] <- out
        }

        # lv's involved in efa
        if (any(type %in% c("lv.efa", "ov.efa"))) {
          if (is.null(partable$efa)) {
            out <- character(0L)
          } else {
            set_names <- lav_partable_efa_values(partable)
            out <- unique(partable$lhs[partable$op == "=~" &
              block_ind &
              partable$efa %in% set_names])
          }
          if (any(type == "ov.efa")) ov_efa <- out
          if (any(type == "lv.efa")) return_value$lv.efa[[b]] <- out
        }

        # lv's that are random slopes
        if (any("lv.rv" == type)) {
          if (is.null(partable$rv)) {
            out <- character(0L)
          } else {
            out <- unique(partable$lhs[partable$op == "=~" &
              block_ind &
              partable$lhs %in% rv_names])
          }
          return_value$lv.rv[[b]] <- out
        }

        # lv's that are indicators of a higher-order factor
        if (any("lv.ind" == type)) {
          out <- unique(partable$rhs[block_ind &
            partable$op == "=~" &
            partable$rhs %in% lv_names &
            partable$lhs != partable$rhs]) # no phantom lv's
          return_value$lv.ind[[b]] <- out
        }

        # higher-order latent variables
        if (any("lv.ho" == type)) {
          out_ind <- unique(partable$rhs[block_ind &
            partable$op == "=~" &
            partable$rhs %in% lv_names &
            partable$lhs != partable$rhs]) # no phantom lv's
          out <- unique(partable$lhs[block_ind &
            partable$op == "=~" &
            partable$lhs %in% lv_names &
            partable$rhs %in% out_ind])
          return_value$lv.ho[[b]] <- out
        }

        # eqs.y
        eqs_y <- unique(partable$lhs[block_ind &
          partable$op == "~"])

        # store eqs.y
        if (any("eqs.y" == type)) {
          return_value$eqs.y[[b]] <- eqs_y
        }

        # eqs.x
        eqs_x <- unique(partable$rhs[block_ind &
          (partable$op == "~" |
            partable$op == "<~")])

        # store eqs.x
        if (any("eqs.x" == type)) {
          return_value$eqs.x[[b]] <- eqs_x
        }

        # v.ind -- indicators of latent variables
        v_ind <- unique(partable$rhs[block_ind &
          partable$op == "=~"])

        # v.cind -- indicators of composites
        v_cind <- unique(partable$rhs[block_ind &
          partable$op == "<~"])


        # ov.*
        # 1. indicators, which are not latent variables themselves
        ov_ind <- v_ind[!v_ind %in% lv_names2]
        # 1a. indicators of composites
        ov_cind <- v_cind[!v_cind %in% lv_names2]
        # 2. dependent ov's
        ov_y <- eqs_y[!eqs_y %in% c(lv_names2, ov_ind, ov_cind)]
        # 3. independent ov's
        if (lav_partable_nlevels(partable) > 1L && b > 1L) {
          # NEW in 0.6-8: if an 'x' was an 'y' in a previous level,
          #               treat it as 'y'
          tmp_eqs_y <- unique(partable$lhs[partable$op == "~"]) # all blocks
          ov_x <- eqs_x[!eqs_x %in% c(lv_names2, ov_ind, ov_cind, tmp_eqs_y)]
        } else {
          ov_x <- eqs_x[!eqs_x %in% c(lv_names2, ov_ind, ov_cind, ov_y)]
        }
        # new in 0.6-12: if we have interaction terms in ov.x, check
        # if some terms are in eqs.y; if so, remove the interaction term
        # from ov.x
        int_idx <- which(grepl(":", ov_x, fixed = TRUE))
        bad_idx <- integer(0L)
        for (iv in int_idx) {
          tmp_names <- strsplit(ov_x[iv], ":", fixed = TRUE)[[1L]]
          if (any(tmp_names %in% eqs_y)) {
            bad_idx <- c(bad_idx, iv)
          }
        }
        if (length(bad_idx) > 0L) {
          ov_y <- unique(c(ov_y, ov_x[bad_idx]))
          # it may be removed later, but needed to construct ov.names
          ov_x <- ov_x[-bad_idx]
        }

        # observed variables
        # easy approach would be: everything that is not in lv.names,
        # but the main purpose here is to 'order' the observed variables
        # according to 'type' (indicators, ov.y, ov.x, orphans)

        # 4. orphaned covariances
        ov_cov <- c(
          partable$lhs[block_ind &
            partable$op == "~~" &
            !partable$lhs %in% lv_names2],
          partable$rhs[block_ind &
            partable$op == "~~" &
            !partable$rhs %in% lv_names2]
        )
        # 5. orphaned intercepts/thresholds
        ov_int <- partable$lhs[block_ind &
          (partable$op == "~1" |
            partable$op == "|") &
          !partable$lhs %in% lv_names2]

        ov_tmp <- c(ov_ind, ov_cind, ov_y, ov_x)
        ov_extra <- unique(c(ov_cov, ov_int)) # must be in this order!
        # so that
        # lav_partable_independence
        # retains the same order
        ov_names <- c(ov_tmp, ov_extra[!ov_extra %in% ov_tmp])

        # store ov?
        if (any("ov" == type)) {
          return_value$ov[[b]] <- ov_names
        }

        if (any("ov.ind" == type)) {
          return_value$ov.ind[[b]] <- ov_ind
        }

        if (any("ov.cind" == type)) {
          return_value$ov.cind[[b]] <- ov_cind
        }

        if (any("ov.interaction" == type)) {
          ov_int_names <- ov_names[grepl(":", ov_names, fixed = TRUE)]
          n_int <- length(ov_int_names)
          if (n_int > 0L) {
            ov_names_noint <- ov_names[!ov_names %in% ov_int_names]

            ok_idx <- logical(n_int)
            for (iv in seq_len(n_int)) {
              tmp_names <- strsplit(ov_int_names[iv], ":", fixed = TRUE)[[1L]]

              # two scenario's:
              # - both variables are in ov.names.noint (ok)
              # - at least one variable is NOT in ov.names.noint (ignore)
              if (all(tmp_names %in% ov_names_noint)) {
                ok_idx[iv] <- TRUE
              }
            }
            ov_interaction <- ov_int_names[ok_idx]
          } else {
            ov_interaction <- character(0L)
          }

          return_value$ov.interaction[[b]] <- ov_interaction
        }

        if (any("ov.efa" == type)) {
          ov_efa_1 <- partable$rhs[partable$op == "=~" &
            block_ind &
            partable$rhs %in% ov_ind &
            partable$lhs %in% ov_efa]
          return_value$ov.efa[[b]] <- unique(ov_efa_1)
        }


        # exogenous `x' covariates
        if (any(type %in% c(
          "ov.x", "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          # correction: is any of these ov.names.x mentioned as a variance,
          #             covariance, or intercept?
          # this should trigger a warning in lav_model_partable()
          if (is.null(partable$user)) { # FLAT!
            partable$user <- rep(1L, length(partable$lhs))
          }
          vars <- c(
            partable$lhs[block_ind &
              partable$op == "~1" &
              partable$user == 1],
            partable$lhs[block_ind &
              partable$op == "~~" &
              partable$user == 1],
            partable$rhs[block_ind &
              partable$op == "~~" &
              partable$user == 1]
          )
          idx_no_x <- which(ov_x %in% vars)
          if (length(idx_no_x)) {
            if (ov_x_fatal) {
              lav_msg_stop(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s];  Please remove them
                and try again.", lav_msg_view(ov_x[idx_no_x], "none")))
            }
              lav_msg_warn(gettextf(
                "model syntax contains variance/covariance/intercept formulas
                involving (an) exogenous variable(s): [%s]; these variables will
                now be treated as random introducing additional free parameters.
                If you wish to treat those variables as fixed, remove these
                formulas from the model syntax. Otherwise, consider adding the
                fixed.x = FALSE option.", lav_msg_view(ov_x[idx_no_x], "none")))
            ov_x <- ov_x[-idx_no_x]
          }
          ov_tmp_x <- ov_x

          # extra
          if (!is.null(partable$exo)) {
            ov_cov <- c(
              partable$lhs[block_ind &
                partable$op == "~~" &
                partable$exo == 1L],
              partable$rhs[block_ind &
                partable$op == "~~" &
                partable$exo == 1L]
            )
            ov_int <- partable$lhs[block_ind &
              partable$op == "~1" &
              partable$exo == 1L]
            ov_extra <- unique(c(ov_cov, ov_int))
            ov_tmp_x <- c(ov_tmp_x, ov_extra[!ov_extra %in% ov_tmp_x])
          }

          ov_names_x <- ov_tmp_x
        }

        # store ov.x?
        if (any("ov.x" == type)) {
          return_value$ov.x[[b]] <- ov_names_x
        }

        # story ov.orphan?
        if (any("ov.orphan" == type)) {
          return_value$ov.orphan[[b]] <- ov_extra
        }

        # ov's without ov.x
        if (any(type %in% c(
          "ov.nox", "ov.model",
          "th.mean", "lv.nonnormal"
        ))) {
          ov_names_nox <- ov_names[!ov_names %in% ov_names_x]
        }

        # store ov.nox
        if (any("ov.nox" == type)) {
          return_value$ov.nox[[b]] <- ov_names_nox
        }

        # store ov.model
        if (any("ov.model" == type)) {
          # if no conditional.x, this is just ov
          # else, this is ov.nox
          if (any(block_ind & partable$op == "~" &
            partable$exo == 1L)) {
            return_value$ov.model[[b]] <- ov_names_nox
          } else {
            return_value$ov.model[[b]] <- ov_names
          }
        }

        # ov's strictly ordered
        if (any(type %in% c(
          "ov.ord", "th", "th.mean",
          "ov.num", "lv.nonnormal"
        ))) {
          tmp <- unique(partable$lhs[block_ind &
            partable$op == "|"])
          ord_names <- ov_names[ov_names %in% tmp]
        }

        if (any("ov.ord" == type)) {
          return_value$ov.ord[[b]] <- ord_names
        }

        # ov's strictly numeric
        if (any(type %in% c("ov.num", "lv.nonnormal"))) {
          ov_num <- ov_names[!ov_names %in% ord_names]
        }

        if (any("ov.num" == type)) {
          return_value$ov.num[[b]] <- ov_num
        }

        # nonnormal lv's
        if (any("lv.nonnormal" == type)) {
          # regular lv's
          lv_reg <- unique(partable$lhs[block_ind &
                           partable$op == "=~" &
                           partable$lhs != partable$rhs])
          if (length(lv_reg) > 0L) {
            out <- unlist(lapply(lv_reg, function(x) {
              # get indicators for this lv
              tmp_ind <- unique(partable$rhs[block_ind &
                partable$op == "=~" &
                partable$lhs == x])
              if (!all(tmp_ind %in% ov_num)) {
                x
              } else {
                character(0)
              }
            }), use.names = FALSE)
            return_value$lv.nonnormal[[b]] <- out
          } else {
            return_value$lv.nonnormal[[b]] <- character(0)
          }
        }

        if (any(c("th", "th.mean") %in% type)) {
          tmp_th_lhs <- partable$lhs[block_ind &
            partable$op == "|"]
          tmp_th_rhs <- partable$rhs[block_ind &
            partable$op == "|"]
        }

        # threshold
        if (any("th" == type)) {
          if (length(ord_names) > 0L) {
            # return in the right order (following ord.names!)
            out <- unlist(lapply(ord_names, function(x) {
              idx <- which(x == tmp_th_lhs)
              tmp_th <- unique(paste(tmp_th_lhs[idx], "|",
                tmp_th_rhs[idx],
                sep = ""
              ))
              # make sure the th's are in increasing order
              # sort(tmp.th)
              # NO!, don't do that; t10 will be before t2
              # fixed in 0.6-1 (bug report from Myrsini)
              # in 0.6-12, we do this anyway like this:

              # get var name
              tmp_th1 <- sapply(
                strsplit(tmp_th, split = "\\|t"),
                "[[", 1
              )
              # get number, and sort
              tmp_th2 <- as.character(sort(as.integer(sapply(
                strsplit(tmp_th, split = "\\|t"), "[[", 2
              ))))
              # paste back together in the right order
              paste(tmp_th1, tmp_th2, sep = "|t")
            }), use.names = FALSE)
          } else {
            out <- character(0L)
          }
          return_value$th[[b]] <- out
        }

        # thresholds and mean/intercepts of numeric variables
        if (any("th.mean" == type)) {
          # if fixed.x -> use ov.names.nox
          # else -> use ov.names
          if (is.null(partable$exo) || all(partable$exo == 0L)) {
            tmp_ov_names <- ov_names
          } else {
            tmp_ov_names <- ov_names_nox
          }
          if (length(tmp_ov_names) > 0L) {
            # return in the right order (following ov.names.nox!)
            out <- unlist(lapply(tmp_ov_names, function(x) {
              if (x %in% ord_names) {
                idx <- which(x == tmp_th_lhs)
                tmp_th <- unique(paste(tmp_th_lhs[idx], "|",
                  tmp_th_rhs[idx],
                  sep = ""
                ))
                # make sure the th's are in increasing order
                # get var name
                tmp_th1 <- sapply(
                  strsplit(tmp_th, split = "\\|t"),
                  "[[", 1
                )
                # get number, and sort
                tmp_th2 <- as.character(sort(as.integer(sapply(
                  strsplit(tmp_th, split = "\\|t"), "[[", 2
                ))))
                # paste back together in the right order
                paste(tmp_th1, tmp_th2, sep = "|t")
              } else {
                x
              }
            }), use.names = FALSE)
          } else {
            out <- character(0L)
          }
          return_value$th.mean[[b]] <- out
        }

        # exogenous lv's
        if (any(c("lv.x", "lv.nox") %in% type)) {
          tmp <- lv_names[!lv_names %in% c(v_ind, eqs_y)]
          lv_names_x <- lv_names[lv_names %in% tmp]
        }

        if (any("lv.x" == type)) {
          return_value$lv.x[[b]] <- lv_names_x
        }

        # dependent ov (but not also indicator or x)
        if (any("ov.y" == type)) {
          tmp <- eqs_y[!eqs_y %in% c(v_ind, eqs_x, lv_names)]
          return_value$ov.y[[b]] <- ov_names[ov_names %in% tmp]
        }

        # dependent lv (but not also indicator or x)
        if (any("lv.y" == type)) {
          tmp <- eqs_y[!eqs_y %in% c(v_ind, eqs_x) &
            eqs_y %in% lv_names]
          return_value$lv.y[[b]] <- lv_names[lv_names %in% tmp]
        }

        # non-exogenous latent variables
        if (any("lv.nox" == type)) {
          return_value$lv.nox[[b]] <- lv_names[!lv_names %in% lv_names_x]
        }

        # marker indicator (if any) for each lv
        if (any("lv.marker" == type)) {
          # default: "" per lv
          out <- character(length(lv_names))
          names(out) <- lv_names
          for (l in seq_along(lv_names)) {
            this_lv_name <- lv_names[l]
            # try to see if we can find a 'marker' indicator for this factor
            # here defined as: has a fixed-to-one factor loading
            marker_idx <- which(block_ind &
              partable$op == "=~" &
              partable$lhs == this_lv_name &
              partable$rhs %in% v_ind &
              partable$ustart == 1L &
              partable$free == 0L)
            if (length(marker_idx) > 0L) {
              # we may have multiple potential markers
              potential_markers <- partable$rhs[marker_idx]
              valid_markers <- character(0L)
              for (m in seq_along(potential_markers)) {
                this_marker_idx <- marker_idx[m]
                # check if 'other' loadings are fixed to zero
                other_idx <- which(block_ind &
                  partable$op == "=~" &
                  partable$lhs != this_lv_name &
                  partable$rhs == partable$rhs[this_marker_idx] &
                  partable$free == 0L)
                if (length(other_idx) == 0L) {
                  # simple structure, or one factor
                  valid_markers <- c(valid_markers,
                                     partable$rhs[this_marker_idx])
                } else if (all(partable$ustart[other_idx] == 0)) {
                  valid_markers <- c(valid_markers,
                                     partable$rhs[this_marker_idx])
                }
              }
              if (length(valid_markers) > 0L) {
                out[l] <- valid_markers[1L] # pick the first one
              }
            }
          } # l
          return_value$lv.marker[[b]] <- out
        }
      }
    } # b

    # ----- lav_partable_vnames ---- no cache ------------------------
    # new in 0.6-14: if 'da' operator, change order! (for ov.order = "data")
    # now via attribute "ovda"
    ov_names_data <- attr(partable, "ovda")
    if (!is.null(ov_names_data)) {
      return_value <- lapply(return_value, function(x) {
        for (b in seq_along(x)) {
          m <- match(x[[b]], ov_names_data)
          target_idx <- which(!is.na(m))
          if (length(target_idx) > 1L) {
            x[[b]][target_idx] <- x[[b]][target_idx][order(m[target_idx])]
          }
        }
        x
      })
    }
  }
  # ----- lav_partable_vnames ---- common ------ 1 type --------
  # to mimic old behaviour, if length(type) == 1L
  if (length(type) == 1L) {
    return_value <- return_value[[type]]
    # to mimic old behaviour, if specific block is requested
    if (ndotdotdot == 0L) {
      if (type == "lv.marker") {
        return_value <- unlist(return_value)
        # no unique(), as unique() drops attributes, and reduces
        # c("", "", "") to a single ""
        # (but, say for 2 groups, you get 2 copies)
        # as this is only for 'display', we leave it like that
      } else {
        return_value <- unique(unlist(return_value))
      }
    } else if (length(block_select) == 1L) {
      return_value <- return_value[[block_select]]
    } else {
      return_value <- return_value[block_select]
    }
  }
  # ----- lav_partable_vnames ---- common ------------------------
  return_value
}
