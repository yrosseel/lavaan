# generate labels for each parameter
lav_pt_labels <- function(partable,
                                blocks = c("group", "level"),
                                group_equal = "", group_partial = "",
                                type = "user") {
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
      ngroups <- lav_pt_ngroups(partable)
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
  # cat("group_equal = ", group_equal, "\n")
  # cat("group_partial = ", group_partial, "\n")

  # use group_equal so that equal sets of parameters get the same label
  if (ngroups > 1L && length(group_equal) > 0L) {
    if ("intercepts" %in% group_equal ||
      "residuals" %in% group_equal ||
      "residual.covariances" %in% group_equal) {
      ov_names_nox <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        ov_names_nox[[g]] <- unique(unlist(
          lav_pt_vnames(partable, "ov.nox", group = g)))
      }
    }
    if ("thresholds" %in% group_equal) {
      ov_names_ord <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        ov_names_ord[[g]] <- unique(unlist(
          lav_pt_vnames(partable, "ov.ord", group = g)))
      }
    }
    if ("means" %in% group_equal ||
      "lv.variances" %in% group_equal ||
      "lv.covariances" %in% group_equal) {
      lv_names <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        lv_names[[g]] <- unique(unlist(
          lav_pt_vnames(partable, "lv", group = g)))
      }
    }

    # eq_flag: TRUE if the parameter is a candidate for a (cross-group)
    # equality constraint. Note that the models need NOT be the same across
    # groups (eg when the 'group:' syntax is used): a parameter qualifies in
    # a given group if it belongs to the requested set *in that group*.
    eq_flag <- logical(length(partable$lhs))

    # helper: TRUE for a row if its 'lhs' is in 'names_list[[group]]' of the
    # group that the row belongs to
    in_group_names <- function(names_list) {
      out <- logical(length(partable$lhs))
      for (g in seq_len(ngroups)) {
        out[partable$group == g &
          partable$lhs %in% names_list[[g]]] <- TRUE
      }
      out
    }

    # LOADINGS
    if ("loadings" %in% group_equal) {
      eq_flag[partable$op == "=~"] <- TRUE
    }
    # COMPOSITE LOADINGS (new in 0.6-4)
    if ("composite.loadings" %in% group_equal) {
      # new setting (0.6-20): <~
      if (any(partable$op == "<~")) {
        lav_msg_warn(gettext("composite.loadings are in fact composite weights;
                              better use composite.weights"))
        eq_flag[partable$op == "<~"] <- TRUE
      } else {
        # old school: composites are phantom constructs with zero residual...
        lv_f_names <- unique(unlist(
          lav_pt_vnames(partable, "lv.formative")))
        eq_flag[partable$op == "~" &
                partable$lhs %in% lv_f_names] <- TRUE
      }
    }
    # COMPOSITE WEIGHTS (new in 0.6-20) # same as 'loadings'...
    if ("composite.weights" %in% group_equal) {
      eq_flag[partable$op == "<~"] <- TRUE
    }
    # INTERCEPTS (OV)
    if ("intercepts" %in% group_equal) {
      eq_flag[partable$op == "~1" & in_group_names(ov_names_nox)] <- TRUE
    }
    # THRESHOLDS (OV-ORD)
    if ("thresholds" %in% group_equal) {
      eq_flag[partable$op == "|" & in_group_names(ov_names_ord)] <- TRUE
    }
    # MEANS (LV)
    if ("means" %in% group_equal) {
      eq_flag[partable$op == "~1" & in_group_names(lv_names)] <- TRUE
    }
    # REGRESSIONS
    if ("regressions" %in% group_equal) {
      eq_flag[partable$op == "~"] <- TRUE
    }
    # RESIDUAL variances (FIXME: OV ONLY!)
    if ("residuals" %in% group_equal) {
      eq_flag[partable$op == "~~" & in_group_names(ov_names_nox) &
        partable$lhs == partable$rhs] <- TRUE
    }
    # RESIDUAL covariances (FIXME: OV ONLY!)
    if ("residual.covariances" %in% group_equal) {
      eq_flag[partable$op == "~~" & in_group_names(ov_names_nox) &
        partable$lhs != partable$rhs] <- TRUE
    }
    # LV VARIANCES
    if ("lv.variances" %in% group_equal) {
      eq_flag[partable$op == "~~" & in_group_names(lv_names) &
        partable$lhs == partable$rhs] <- TRUE
    }
    # LV COVARIANCES
    if ("lv.covariances" %in% group_equal) {
      eq_flag[partable$op == "~~" & in_group_names(lv_names) &
        partable$lhs != partable$rhs] <- TRUE
    }

    # if group_partial, do NOT constrain these parameters (in any group)
    if (length(group_partial) > 0L) {
      base_label <- paste(partable$lhs, partable$op, partable$rhs, sep = "")
      eq_flag[base_label %in% group_partial] <- FALSE
    }

    # Give all candidate parameters that share the same lhs/op/rhs the same
    # label, so they become equality-constrained -- across *every* group that
    # contains them, independently of group 1. Parameters that appear in a
    # single group only keep their unique label (and stay free). This way the
    # constraints are symmetric: a parameter present in groups 2 and 3 (but
    # not group 1) is constrained just as one present in groups 1 and 2.
    eq_idx <- which(eq_flag)
    if (length(eq_idx) > 0L) {
      eq_key <- paste(partable$lhs, partable$op, partable$rhs,
        sep = "")[eq_idx]
      for (key in unique(eq_key)) {
        same_idx <- eq_idx[eq_key == key]
        if (length(same_idx) > 1L) {
          label[same_idx] <- label[same_idx[1L]]
        }
      }
    }
  }

  # cat("DEBUG: g1.idx = ", g1.idx, "\n")
  # cat("DEBUG: label after group_equal:\n"); print(label); cat("\n")

  # handle other block identifier (not 'group')
  for (block in blocks) {
    if (block == "group") {
      next
    } else if (block == "level" && !is.null(partable[[block]])) {
      # all but first level
      lev_vals <- lav_pt_level_values(partable)
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
