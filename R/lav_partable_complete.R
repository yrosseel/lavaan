# handle bare-minimum partables
# add some additional columns
lav_partable_complete <- function(partable = NULL, start = TRUE) { # nolint
  # check if we hava a data.frame
  # if so, check for columns that are 'factor' and convert them to 'character'
  ovda <- attr(partable, "ovda")
  if (is.data.frame(partable)) {
    fac.idx <- sapply(partable, is.factor)
    partable[fac.idx] <- lapply(partable[fac.idx], as.character)
  }

  # check if we have lhs, op, rhs
  stopifnot(!is.null(partable$lhs),
    !is.null(partable$op),
    !is.null(partable$rhs))

  # number of elements
  tmp.n <- length(partable$lhs)
  if (!is.data.frame(partable)) {
    # check for equal column length
    nel <- sapply(partable, length)
    short.idx <- which(nel < tmp.n)
    long.idx <- which(nel > tmp.n)
    if (length(long.idx) > 0L) {
      lav_msg_warn(gettext("partable columns have unequal length"))
    }
    if (length(short.idx) > 0L) {
      # try to extend them in a 'natural' way
      for (i in short.idx) {
        too.short <- tmp.n - nel[i]
        if (is.integer(partable[[i]])) {
          partable[[i]] <- c(partable[[i]],
            integer(too.short))
        } else if (is.numeric(partable[[i]])) {
          partable[[i]] <- c(partable[[i]],
            numeric(too.short))
        } else {
          partable[[i]] <- c(partable[[i]],
            character(too.short))
        }
      }
    }
  }

  # create new id column
  # if(is.null(partable$id)) {
  partable$id <- seq_len(tmp.n)
  # }

  # add user column
  if (is.null(partable$user)) {
    partable$user <- rep(1L, tmp.n)
  } else {
    partable$user <- as.integer(partable$user)
  }

  # add block column
  if (is.null(partable$block)) {
    partable$block <- rep(1L, tmp.n)
  } else {
    partable$block <- as.integer(partable$block)
  }

  # add group column
  if (is.null(partable$group)) {
    partable$group <- rep(1L, tmp.n)
  } else {
    # partable$group <- as.integer(partable$group) # maybe labels?
  }

  # add free column
  if (is.null(partable$free)) {
    partable$free <- seq_len(tmp.n)
    # 0.6-11: check for simple equality constraints
    #         note: this is perhaps only a subset (eg SAM!) of a larger
    #         table, and we have to renumber the 'free' column
  } else if (is.integer(partable$free) &&
    any(partable$free > 0L)   &&
    !any(partable$op == "==") &&
    !is.null(partable$label)  &&
    !is.null(partable$plabel) &&
    any(duplicated(partable$free[partable$free > 0L]))) {
    dup.idx <- which(partable$free > 0L & duplicated(partable$free))
    all.idx <- which(partable$free %in% unique(partable$free[dup.idx]))
    eq.labels <- unique(partable$free[all.idx])
    eq.id <- integer(length(partable$lhs))
    eq.id[all.idx] <- partable$free[all.idx]
    partable$free[dup.idx] <- 0L
    idx.free <- which(partable$free > 0L)
    partable$free <- rep(0L, tmp.n)
    partable$free[idx.free] <- seq_along(idx.free)
    for (eq.label in eq.labels) {
      all.idx <- which(eq.id == eq.label)
      ref.idx <- all.idx[1L]
      other.idx <- all.idx[-1L]
      partable$free[other.idx] <- partable$free[ref.idx]
    }
  } else {
    # treat non-zero as 'free'
    free.idx <- which(as.logical(partable$free))
    partable$free <- rep(0L, tmp.n)
    if (length(free.idx) > 0L) {
      partable$free[free.idx] <- seq_len(length(free.idx))
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
      partable$ustart <- rep(as.numeric(NA), tmp.n)
      non.free <- which(!partable$free)
      if (length(non.free)) {
        partable$ustart[non.free] <- 0
      }
    }
  } else {
    partable$ustart <- as.numeric(partable$ustart)
  }

  # add exo column
  if (is.null(partable$exo)) {
    partable$exo <- rep(0, tmp.n)
  } else {
    partable$exo <- as.integer(partable$exo)
  }

  # add label column
  if (is.null(partable$label)) {
    partable$label <- rep("", tmp.n)
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
