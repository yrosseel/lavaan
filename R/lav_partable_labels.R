# generate labels for each parameter
lav_partable_labels <- function(partable,                                # nolint start
                                blocks = c("group", "level"),
                                group.equal = "", group.partial = "",
                                type = "user") {                         # nolint end
  # catch empty partable
  if (length(partable$lhs) == 0L) {
    return(character(0L))
  }

  # default labels
  label <- paste(partable$lhs, partable$op, partable$rhs, sep = "")

  # handle multiple groups
  if ("group" %in% blocks) {
    if (is.character(partable$group)) {
      group_label <- unique(partable$group)
      group_label <- group_label[nchar(group_label) > 0L]
      ngroups <- length(group_label)
    } else {
      ngroups <- lav_partable_ngroups(partable)
      group_label <- 1:ngroups
    }
    if (ngroups > 1L) {
      for (g in 2:ngroups) {
        label[partable$group == group_label[g]] <-
          paste(label[partable$group == group_label[g]],
            ".g", g,
            sep = ""
          )
      }
    }
  } else {
    ngroups <- 1L
  }

  # cat("DEBUG: label start:\n"); print(label); cat("\n")
  # cat("group.equal = ", group.equal, "\n")
  # cat("group.partial = ", group.partial, "\n")

  # use group.equal so that equal sets of parameters get the same label
  if (ngroups > 1L && length(group.equal) > 0L) {
    if ("intercepts" %in% group.equal ||
      "residuals" %in% group.equal ||
      "residual.covariances" %in% group.equal) {
      ov_names_nox <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        ov_names_nox[[g]] <- unique(unlist(
          lav_partable_vnames(partable, "ov.nox", group = g)))
      }
    }
    if ("thresholds" %in% group.equal) {
      ov_names_ord <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        ov_names_ord[[g]] <- unique(unlist(
          lav_partable_vnames(partable, "ov.ord", group = g)))
      }
    }
    if ("means" %in% group.equal ||
      "lv.variances" %in% group.equal ||
      "lv.covariances" %in% group.equal) {
      lv_names <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        lv_names[[g]] <- unique(unlist(
          lav_partable_vnames(partable, "lv", group = g)))
      }
    }

    # g1.flag: TRUE if included, FALSE if not
    g1_flag <- logical(length(partable$lhs))

    # LOADINGS
    if ("loadings" %in% group.equal) {
      g1_flag[partable$op == "=~" & partable$group == 1L] <- TRUE
    }
    # COMPOSITE LOADINGS (new in 0.6-4)
    if ("composite.loadings" %in% group.equal) {
      # new setting (0.6-20): <~
      if (any(partable$op == "<~" & partable$group == 1L)) {
        lav_msg_warn(gettext("composite.loadings are in fact composite weights;
                              better use composite.weights"))
        g1_flag[partable$op == "<~" & partable$group == 1L] <- TRUE
      } else {
        # old school: composites are phantom constructs with zero residual...
        lv_f_names <- unique(unlist(
          lav_partable_vnames(partable, "lv.formative")))
        g1_flag[partable$op == "~" &
                partable$lhs %in% lv_f_names &
                partable$group == 1L] <- TRUE
      }
    }
    # COMPOSITE WEIGHTS (new in 0.6-20) # same as 'loadings'...
    if ("composite.weights" %in% group.equal) {
      g1_flag[partable$op == "<~" & partable$group == 1L] <- TRUE
    }
    # INTERCEPTS (OV)
    if ("intercepts" %in% group.equal) {
      g1_flag[partable$op == "~1" & partable$group == 1L &
        partable$lhs %in% ov_names_nox[[1L]]] <- TRUE
    }
    # THRESHOLDS (OV-ORD)
    if ("thresholds" %in% group.equal) {
      g1_flag[partable$op == "|" & partable$group == 1L &
        partable$lhs %in% ov_names_ord[[1L]]] <- TRUE
    }
    # MEANS (LV)
    if ("means" %in% group.equal) {
      g1_flag[partable$op == "~1" & partable$group == 1L &
        partable$lhs %in% lv_names[[1L]]] <- TRUE
    }
    # REGRESSIONS
    if ("regressions" %in% group.equal) {
      g1_flag[partable$op == "~" & partable$group == 1L] <- TRUE
    }
    # RESIDUAL variances (FIXME: OV ONLY!)
    if ("residuals" %in% group.equal) {
      g1_flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% ov_names_nox[[1L]] &
        partable$lhs == partable$rhs] <- TRUE
    }
    # RESIDUAL covariances (FIXME: OV ONLY!)
    if ("residual.covariances" %in% group.equal) {
      g1_flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% ov_names_nox[[1L]] &
        partable$lhs != partable$rhs] <- TRUE
    }
    # LV VARIANCES
    if ("lv.variances" %in% group.equal) {
      g1_flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% lv_names[[1L]] &
        partable$lhs == partable$rhs] <- TRUE
    }
    # LV COVARIANCES
    if ("lv.covariances" %in% group.equal) {
      g1_flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% lv_names[[1L]] &
        partable$lhs != partable$rhs] <- TRUE
    }

    # if group.partial, set corresponding flag to FALSE
    if (length(group.partial) > 0L) {
      g1_flag[label %in% group.partial &
        partable$group == 1L] <- FALSE
    }

    # for each (constrained) parameter in 'group 1', find a similar one
    # in the other groups (we assume here that the models need
    # NOT be the same across groups!
    g1_idx <- which(g1_flag)
    for (i in seq_along(g1_idx)) {
      ref_idx <- g1_idx[i]
      idx <- which(partable$lhs == partable$lhs[ref_idx] &
        partable$op == partable$op[ref_idx] &
        partable$rhs == partable$rhs[ref_idx] &
        partable$group > 1L)
      label[idx] <- label[ref_idx]
    }
  }

  # cat("DEBUG: g1.idx = ", g1.idx, "\n")
  # cat("DEBUG: label after group.equal:\n"); print(label); cat("\n")

  # handle other block identifier (not 'group')
  for (block in blocks) {
    if (block == "group") {
      next
    } else if (block == "level" && !is.null(partable[[block]])) {
      # all but first level
      lev_vals <- lav_partable_level_values(partable)
      idx <- which(partable[[block]] != lev_vals[1])
      label[idx] <- paste(label[idx], ".", "l",
        partable[[block]][idx],
        sep = ""
      )
    } else if (!is.null(partable[[block]])) {
      label <- paste(label, ".", block, partable[[block]], sep = "")
    }
  }

  # user-specified labels -- override everything!!
  user_idx <- which(nchar(partable$label) > 0L)
  label[user_idx] <- partable$label[user_idx]

  # cat("DEBUG: user.idx = ", user.idx, "\n")
  # cat("DEBUG: label after user.idx:\n"); print(label); cat("\n")

  # which labels do we need?
  if (type == "user") {
    idx <- seq_along(label)
  } else if (type == "free") {
    # idx <- which(partable$free > 0L & !duplicated(partable$free))
    idx <- which(partable$free > 0L)
    # } else if(type == "unco") {
    #    idx <- which(partable$unco > 0L & !duplicated(partable$unco))
  } else {
    lav_msg_stop(gettext("argument `type' must be one of free or user"))
  }

  label[idx]
}
