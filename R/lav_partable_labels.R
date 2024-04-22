# generate labels for each parameter
lav_partable_labels <- function(partable,
                                blocks = c("group", "level"),
                                group.equal = "", group.partial = "",
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
      group.label <- unique(partable$group)
      group.label <- group.label[nchar(group.label) > 0L]
      ngroups <- length(group.label)
    } else {
      ngroups <- lav_partable_ngroups(partable)
      group.label <- 1:ngroups
    }
    if (ngroups > 1L) {
      for (g in 2:ngroups) {
        label[partable$group == group.label[g]] <-
          paste(label[partable$group == group.label[g]],
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
      ov.names.nox <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        ov.names.nox[[g]] <- unique(unlist(lav_partable_vnames(partable, "ov.nox", group = g)))
      }
    }
    if ("thresholds" %in% group.equal) {
      ov.names.ord <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        ov.names.ord[[g]] <- unique(unlist(lav_partable_vnames(partable, "ov.ord", group = g)))
      }
    }
    if ("means" %in% group.equal ||
      "lv.variances" %in% group.equal ||
      "lv.covariances" %in% group.equal) {
      lv.names <- vector("list", length = ngroups)
      for (g in 1:ngroups) {
        lv.names[[g]] <- unique(unlist(lav_partable_vnames(partable, "lv", group = g)))
      }
    }

    # g1.flag: TRUE if included, FALSE if not
    g1.flag <- logical(length(partable$lhs))

    # LOADINGS
    if ("loadings" %in% group.equal) {
      g1.flag[partable$op == "=~" & partable$group == 1L] <- TRUE
    }
    # COMPOSITE LOADINGS (new in 0.6-4)
    if ("composite.loadings" %in% group.equal) {
      g1.flag[partable$op == "<~" & partable$group == 1L] <- TRUE
    }
    # INTERCEPTS (OV)
    if ("intercepts" %in% group.equal) {
      g1.flag[partable$op == "~1" & partable$group == 1L &
        partable$lhs %in% ov.names.nox[[1L]]] <- TRUE
    }
    # THRESHOLDS (OV-ORD)
    if ("thresholds" %in% group.equal) {
      g1.flag[partable$op == "|" & partable$group == 1L &
        partable$lhs %in% ov.names.ord[[1L]]] <- TRUE
    }
    # MEANS (LV)
    if ("means" %in% group.equal) {
      g1.flag[partable$op == "~1" & partable$group == 1L &
        partable$lhs %in% lv.names[[1L]]] <- TRUE
    }
    # REGRESSIONS
    if ("regressions" %in% group.equal) {
      g1.flag[partable$op == "~" & partable$group == 1L] <- TRUE
    }
    # RESIDUAL variances (FIXME: OV ONLY!)
    if ("residuals" %in% group.equal) {
      g1.flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% ov.names.nox[[1L]] &
        partable$lhs == partable$rhs] <- TRUE
    }
    # RESIDUAL covariances (FIXME: OV ONLY!)
    if ("residual.covariances" %in% group.equal) {
      g1.flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% ov.names.nox[[1L]] &
        partable$lhs != partable$rhs] <- TRUE
    }
    # LV VARIANCES
    if ("lv.variances" %in% group.equal) {
      g1.flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% lv.names[[1L]] &
        partable$lhs == partable$rhs] <- TRUE
    }
    # LV COVARIANCES
    if ("lv.covariances" %in% group.equal) {
      g1.flag[partable$op == "~~" & partable$group == 1L &
        partable$lhs %in% lv.names[[1L]] &
        partable$lhs != partable$rhs] <- TRUE
    }

    # if group.partial, set corresponding flag to FALSE
    if (length(group.partial) > 0L) {
      g1.flag[label %in% group.partial &
        partable$group == 1L] <- FALSE
    }

    # for each (constrained) parameter in 'group 1', find a similar one
    # in the other groups (we assume here that the models need
    # NOT be the same across groups!
    g1.idx <- which(g1.flag)
    for (i in 1:length(g1.idx)) {
      ref.idx <- g1.idx[i]
      idx <- which(partable$lhs == partable$lhs[ref.idx] &
        partable$op == partable$op[ref.idx] &
        partable$rhs == partable$rhs[ref.idx] &
        partable$group > 1L)
      label[idx] <- label[ref.idx]
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
  user.idx <- which(nchar(partable$label) > 0L)
  label[user.idx] <- partable$label[user.idx]

  # cat("DEBUG: user.idx = ", user.idx, "\n")
  # cat("DEBUG: label after user.idx:\n"); print(label); cat("\n")

  # which labels do we need?
  if (type == "user") {
    idx <- 1:length(label)
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
