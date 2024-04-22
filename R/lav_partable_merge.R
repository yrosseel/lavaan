# merge two parameter tables
# - but allow different number of columns
lav_partable_merge <- function(pt1 = NULL, pt2 = NULL,
                               remove.duplicated = FALSE,
                               fromLast = FALSE,
                               warn = TRUE) {
  # check for empty pt2
  if (is.null(pt2) || length(pt2) == 0L) {
    return(pt1)
  }

  pt1 <- as.data.frame(pt1, stringsAsFactors = FALSE)
  pt2 <- as.data.frame(pt2, stringsAsFactors = FALSE)

  # check minimum requirements: lhs, op, rhs
  stopifnot(
    !is.null(pt1$lhs), !is.null(pt1$op), !is.null(pt1$rhs),
    !is.null(pt2$lhs), !is.null(pt2$op), !is.null(pt2$rhs)
  )

  # both should have block (or not)
  if (is.null(pt1$block) && is.null(pt2$block)) {
    pt1$block <- rep(1L, length(pt1$lhs))
    pt2$block <- rep(1L, length(pt2$lhs))
    TMP <- rbind(
      pt1[, c("lhs", "op", "rhs", "block")],
      pt2[, c("lhs", "op", "rhs", "block")]
    )
  } else {
    if (is.null(pt1$block) && !is.null(pt2$block)) {
      pt1$block <- rep(1L, length(pt1$lhs))
    } else if (is.null(pt2$block) && !is.null(pt1$block)) {
      pt2$block <- rep(1L, length(pt2$lhs))
    }
    TMP <- rbind(
      pt1[, c("lhs", "op", "rhs", "block")],
      pt2[, c("lhs", "op", "rhs", "block")]
    )
  }

  # if missing columns, provide default values of the right type
  # (numeric/integer/character)

  # group
  if (is.null(pt1$group) && !is.null(pt2$group)) {
    pt1$group <- rep(1L, length(pt1$lhs))
  } else if (is.null(pt2$group) && !is.null(pt1$group)) {
    pt2$group <- rep(1L, length(pt2$lhs))
  }

  # level
  if (is.null(pt1$level) && !is.null(pt2$level)) {
    pt1$level <- rep(1L, length(pt1$lhs))
  } else if (is.null(pt2$level) && !is.null(pt1$level)) {
    pt2$level <- rep(1L, length(pt2$lhs))
  }

  # user
  if (is.null(pt1$user) && !is.null(pt2$user)) {
    pt1$user <- rep(0L, length(pt1$lhs))
  } else if (is.null(pt2$user) && !is.null(pt1$user)) {
    pt2$user <- rep(0L, length(pt2$lhs))
  }

  # free
  if (is.null(pt1$free) && !is.null(pt2$free)) {
    pt1$free <- rep(0L, length(pt1$lhs))
  } else if (is.null(pt2$free) && !is.null(pt1$free)) {
    pt2$free <- rep(0L, length(pt2$lhs))
  }

  # ustart -- set to zero!!
  if (is.null(pt1$ustart) && !is.null(pt2$ustart)) {
    pt1$ustart <- rep(0, length(pt1$lhs))
  } else if (is.null(pt2$ustart) && !is.null(pt1$ustart)) {
    pt2$ustart <- rep(0, length(pt2$lhs))
  }

  # exo
  if (is.null(pt1$exo) && !is.null(pt2$exo)) {
    pt1$exo <- rep(0L, length(pt1$lhs))
  } else if (is.null(pt2$exo) && !is.null(pt1$exo)) {
    pt2$exo <- rep(0L, length(pt2$lhs))
  }

  # label
  if (is.null(pt1$label) && !is.null(pt2$label)) {
    pt1$label <- rep("", length(pt1$lhs))
  } else if (is.null(pt2$label) && !is.null(pt1$label)) {
    pt2$label <- rep("", length(pt2$lhs))
  }

  # plabel
  if (is.null(pt1$plabel) && !is.null(pt2$plabel)) {
    pt1$plabel <- rep("", length(pt1$lhs))
  } else if (is.null(pt2$plabel) && !is.null(pt1$plabel)) {
    pt2$plabel <- rep("", length(pt2$lhs))
  }

  # efa
  if (is.null(pt1$efa) && !is.null(pt2$efa)) {
    pt1$efa <- rep("", length(pt1$lhs))
  } else if (is.null(pt2$efa) && !is.null(pt1$efa)) {
    pt2$efa <- rep("", length(pt2$lhs))
  }

  # start
  if (is.null(pt1$start) && !is.null(pt2$start)) {
    pt1$start <- rep(as.numeric(NA), length(pt1$lhs))
  } else if (is.null(pt2$start) && !is.null(pt1$start)) {
    pt2$start <- rep(as.numeric(NA), length(pt2$lhs))
  }

  # est
  if (is.null(pt1$est) && !is.null(pt2$est)) {
    pt1$est <- rep(0, length(pt1$lhs))
  } else if (is.null(pt2$est) && !is.null(pt1$est)) {
    pt2$est <- rep(0, length(pt2$lhs))
  }


  # check for duplicated elements
  if (remove.duplicated) {
    # if fromLast = TRUE, idx is in pt1
    # if fromLast = FALSE, idx is in pt2
    idx <- which(duplicated(TMP, fromLast = fromLast))

    if (length(idx)) {
      if (warn) {
        lav_msg_warn(
          gettext("duplicated parameters are ignored:"),
          paste(apply(TMP[idx, c("lhs", "op", "rhs")], 1,
            paste,
            collapse = " "
          ), collapse = "\n")
        )
      }
      if (fromLast) {
        pt1 <- pt1[-idx, ]
      } else {
        idx <- idx - nrow(pt1)
        pt2 <- pt2[-idx, ]
      }
    }
  } else if (!is.null(pt1$start) && !is.null(pt2$start)) {
    # copy start values from pt1 to pt2
    for (i in 1:length(pt1$lhs)) {
      idx <- which(pt2$lhs == pt1$lhs[i] &
        pt2$op == pt1$op[i] &
        pt2$rhs == pt1$rhs[i] &
        pt2$block == pt1$block[i])

      pt2$start[idx] <- pt1$start[i]
    }
  }

  # nicely merge, using 'id' column (if it comes first)
  if (is.null(pt1$id) && !is.null(pt2$id)) {
    nid <- max(pt2$id)
    pt1$id <- (nid + 1L):(nid + nrow(pt1))
  } else if (is.null(pt2$id) && !is.null(pt1$id)) {
    nid <- max(pt1$id)
    pt2$id <- (nid + 1L):(nid + nrow(pt2))
  }

  NEW <- base::merge(pt1, pt2, all = TRUE, sort = FALSE)

  # make sure group/block/level are zero (or "") if
  # op %in% c("==", "<", ">", ":=")
  op.idx <- which(NEW$op %in% c("==", "<", ">", ":="))
  if (length(op.idx) > 0L) {
    if (!is.null(NEW$block)) {
      # ALWAYS integer
      NEW$block[op.idx] <- 0L
    }
    if (!is.null(NEW$group)) {
      if (is.character(NEW$level)) {
        NEW$group[op.idx] <- ""
      } else {
        NEW$group[op.idx] <- 0L
      }
    }
    if (!is.null(NEW$level)) {
      if (is.character(NEW$level)) {
        NEW$level[op.idx] <- ""
      } else {
        NEW$level[op.idx] <- 0L
      }
    }
  }

  NEW
}
