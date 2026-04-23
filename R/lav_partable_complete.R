# handle bare-minimum partables
# add some additional columns
lav_partable_complete <- function(partable = NULL, start = TRUE) { # nolint
  # check if we hava a data.frame
  # if so, check for columns that are 'factor' and convert them to 'character'
  ovda <- attr(partable, "ovda")
  if (is.data.frame(partable)) {
    fac_idx <- sapply(partable, is.factor)
    partable[fac_idx] <- lapply(partable[fac_idx], as.character)
  }

  # check if we have lhs, op, rhs
  stopifnot(!is.null(partable$lhs),
    !is.null(partable$op),
    !is.null(partable$rhs))

  # number of elements
  tmp_n <- length(partable$lhs)
  if (!is.data.frame(partable)) {
    # check for equal column length
    nel <- sapply(partable, length)
    short_idx <- which(nel < tmp_n)
    long_idx <- which(nel > tmp_n)
    if (length(long_idx) > 0L) {
      lav_msg_warn(gettext("partable columns have unequal length"))
    }
    if (length(short_idx) > 0L) {
      # try to extend them in a 'natural' way
      for (i in short_idx) {
        too_short <- tmp_n - nel[i]
        if (is.integer(partable[[i]])) {
          partable[[i]] <- c(partable[[i]],
            integer(too_short))
        } else if (is.numeric(partable[[i]])) {
          partable[[i]] <- c(partable[[i]],
            numeric(too_short))
        } else {
          partable[[i]] <- c(partable[[i]],
            character(too_short))
        }
      }
    }
  }

  # create new id column
  # if(is.null(partable$id)) {
  partable$id <- seq_len(tmp_n)
  # }

  # add user column
  if (is.null(partable$user)) {
    partable$user <- rep(1L, tmp_n)
  } else {
    partable$user <- as.integer(partable$user)
  }

  # add block column
  if (is.null(partable$block)) {
    partable$block <- rep(1L, tmp_n)
  } else {
    partable$block <- as.integer(partable$block)
  }

  # add group column
  if (is.null(partable$group)) {
    partable$group <- rep(1L, tmp_n)
  } else {
    # partable$group <- as.integer(partable$group) # maybe labels?
  }

  # add free column
  if (is.null(partable$free)) {
    partable$free <- seq_len(tmp_n)
    # 0.6-11: check for simple equality constraints
    #         note: this is perhaps only a subset (eg SAM!) of a larger
    #         table, and we have to renumber the 'free' column
  } else if (is.integer(partable$free) &&
    any(partable$free > 0L)   &&
    !any(partable$op == "==") &&
    !is.null(partable$label)  &&
    !is.null(partable$plabel) &&
    any(duplicated(partable$free[partable$free > 0L]))) {
    dup_idx <- which(partable$free > 0L & duplicated(partable$free))
    all_idx <- which(partable$free %in% unique(partable$free[dup_idx]))
    eq_labels <- unique(partable$free[all_idx])
    eq_id <- integer(length(partable$lhs))
    eq_id[all_idx] <- partable$free[all_idx]
    partable$free[dup_idx] <- 0L
    idx_free <- which(partable$free > 0L)
    partable$free <- rep(0L, tmp_n)
    partable$free[idx_free] <- seq_along(idx_free)
    for (eq_label in eq_labels) {
      all_idx <- which(eq_id == eq_label)
      ref_idx <- all_idx[1L]
      other_idx <- all_idx[-1L]
      partable$free[other_idx] <- partable$free[ref_idx]
    }
  } else {
    # treat non-zero as 'free'
    free_idx <- which(as.logical(partable$free))
    partable$free <- rep(0L, tmp_n)
    if (length(free_idx) > 0L) {
      partable$free[free_idx] <- seq_along(free_idx)
    }
  }

  # add ustart column
  if (is.null(partable$ustart)) {
    # do we have something else? start? est?
    if (!is.null(partable$start)) {
      partable$ustart <- as.numeric(partable$start)
    } else if (!is.null(partable$est)) {
      partable$ustart <- as.numeric(partable$est)
    } else {
      partable$ustart <- rep(as.numeric(NA), tmp_n)
      non_free <- which(!partable$free)
      if (length(non_free)) {
        partable$ustart[non_free] <- 0
      }
    }
  } else {
    partable$ustart <- as.numeric(partable$ustart)
  }

  # add exo column
  if (is.null(partable$exo)) {
    partable$exo <- rep(0, tmp_n)
  } else {
    partable$exo <- as.integer(partable$exo)
  }

  # add label column
  if (is.null(partable$label)) {
    partable$label <- rep("", tmp_n)
  } else {
    partable$label <- as.character(partable$label)
  }

  # order them nicely: id lhs op rhs group
  idx <- match(c("id", "lhs", "op", "rhs", "user", "block", "group",
    "free", "ustart", "exo", "label"),
  names(partable))
  tmp <- partable[idx]
  partable <- c(tmp, partable[-idx])

  # add start column
  if (start) {
    if (is.null(partable$start)) {
      partable$start <- lav_start(start.method = "simple",
        lavpartable = partable)
    }
  }
  attr(partable, "ovda") <- ovda
  attr(partable, "vnames") <- lav_partable_vnames(partable, "*")
  partable
}
