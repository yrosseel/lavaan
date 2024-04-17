# what are the block values (not necessarily integers)
lav_partable_block_values <- function(partable) {
  if (is.null(partable$block)) {
    block.values <- 1L
  } else {
    # always integers
    tmp <- partable$block[partable$block > 0L & # non-zero only
      !partable$op %in% c("==", "<", ">", ":=")]
    block.values <- unique(na.omit(tmp)) # could be, eg, '2' only
  }

  block.values
}

# guess number of blocks from a partable
lav_partable_nblocks <- function(partable) {
  length(lav_partable_block_values(partable))
}

# what are the group values (not necessarily integers)
lav_partable_group_values <- function(partable) {
  # FLAT?
  if (any(partable$op == ":")) {
    colon.idx <- which(partable$op == ":" &
      tolower(partable$lhs) == "group")
    if (length(colon.idx) > 0L) {
      group.values <- unique(partable$rhs[colon.idx])
    }
    # regular partable
  } else if (is.null(partable$group)) {
    group.values <- 1L
  } else if (is.numeric(partable$group)) {
    tmp <- partable$group[partable$group > 0L &
      !partable$op %in% c("==", "<", ">", ":=")]
    group.values <- unique(na.omit(tmp))
  } else { # character
    tmp <- partable$group[nchar(partable$group) > 0L &
      !partable$op %in% c("==", "<", ">", ":=")]
    group.values <- unique(na.omit(tmp))
  }

  group.values
}

# guess number of groups from a partable
lav_partable_ngroups <- function(partable) {
  length(lav_partable_group_values(partable))
}

# what are the level values (not necessarily integers)
lav_partable_level_values <- function(partable) {
  # FLAT?
  if (any(partable$op == ":")) {
    colon.idx <- which(partable$op == ":" &
      tolower(partable$lhs) == "level")
    level.values <- integer(0L)
    if (length(colon.idx) > 0L) {
      level.values <- unique(partable$rhs[colon.idx])
    }
    # regular partable
  } else if (is.null(partable$level)) {
    level.values <- 1L
  } else if (is.numeric(partable$level)) {
    tmp <- partable$level[partable$level > 0L &
      !partable$op %in% c("==", "<", ">", ":=")]
    level.values <- unique(na.omit(tmp))
  } else { # character
    tmp <- partable$level[nchar(partable$level) > 0L &
      !partable$op %in% c("==", "<", ">", ":=")]
    level.values <- unique(na.omit(tmp))
  }

  level.values
}

# guess number of levels from a partable
lav_partable_nlevels <- function(partable) {
  length(lav_partable_level_values(partable))
}

# efa sets values
lav_partable_efa_values <- function(partable) {
  if (is.null(partable$efa)) {
    efa.values <- character(0L)
  } else { # should be character
    tmp.efa <- as.character(partable$efa)
    tmp <- tmp.efa[nchar(tmp.efa) > 0L &
      !partable$op %in% c("==", "<", ">", ":=")]
    efa.values <- unique(na.omit(tmp))
  }

  efa.values
}

# number of efa sets from a partable
lav_partable_nefa <- function(partable) {
  length(lav_partable_efa_values(partable))
}




# number of sample statistics per block
lav_partable_ndat <- function(partable) {
  # global
  meanstructure <- any(partable$op == "~1")
  fixed.x <- any(partable$exo > 0L & partable$free == 0L)
  conditional.x <- any(partable$exo > 0L & partable$op == "~")
  categorical <- any(partable$op == "|")
  correlation <- any(partable$op == "~*~")
  if (categorical) {
    meanstructure <- TRUE
  }

  # blocks
  nblocks <- lav_partable_nblocks(partable)
  nlevels <- lav_partable_nlevels(partable)
  ndat <- integer(nblocks)

  for (b in seq_len(nblocks)) {
    # how many observed variables in this block?
    if (conditional.x) {
      ov.names <- lav_partable_vnames(partable, "ov.nox", block = b)
    } else {
      ov.names <- lav_partable_vnames(partable, "ov", block = b)
    }
    nvar <- length(ov.names)

    # pstar
    pstar <- nvar * (nvar + 1) / 2
    if (meanstructure) {
      pstar <- pstar + nvar

      # no meanstructure if within level, except ov.x which is not
      # decomposed
      if (nlevels > 1L && (b %% nlevels) == 1L) {
        # all zero
        pstar <- pstar - nvar

        # except within-only 'y'
        ov.names.y <- lav_partable_vnames(partable, "ov.nox", block = b)
        ov.names.y2 <- unlist(lav_partable_vnames(partable, "ov",
          block = seq_len(nblocks)[-b]
        ))
        ov.names.y <- ov.names.y[!ov.names.y %in% ov.names.y2]
        if (length(ov.names.y) > 0L) {
          pstar <- pstar + length(ov.names.y)
        }

        # except within-only 'x' (unless fixed.x)
        ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
        ov.names.x2 <- unlist(lav_partable_vnames(partable, "ov",
          block = seq_len(nblocks)[-b]
        ))
        ov.names.x <- ov.names.x[!ov.names.x %in% ov.names.x2]
        if (!fixed.x && length(ov.names.x) > 0L) {
          pstar <- pstar + length(ov.names.x)
        }
      }
    }

    ndat[b] <- pstar

    # correction for fixed.x?
    if (!conditional.x && fixed.x) {
      ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
      nvar.x <- length(ov.names.x)
      pstar.x <- nvar.x * (nvar.x + 1) / 2
      if (meanstructure) {
        if (nlevels > 1L && (b %% nlevels) == 1L) {
          # do nothing, they are already removed
        } else {
          pstar.x <- pstar.x + nvar.x
        }
      }
      ndat[b] <- ndat[b] - pstar.x
    }

    # correction for ordinal data?
    if (categorical) {
      ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
      nexo <- length(ov.names.x)
      ov.ord <- lav_partable_vnames(partable, "ov.ord", block = b)
      nvar.ord <- length(ov.ord)
      th <- lav_partable_vnames(partable, "th", block = b)
      nth <- length(th)
      # no variances
      ndat[b] <- ndat[b] - nvar.ord
      # no means
      ndat[b] <- ndat[b] - nvar.ord
      # but additional thresholds
      ndat[b] <- ndat[b] + nth
      # add slopes
      if (conditional.x) {
        ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
        nexo <- length(ov.names.x)
        ndat[b] <- ndat[b] + (nvar * nexo)
      }
    }

    # correction for correlation not categorical
    if (correlation && !categorical) {
      ndat[b] <- ndat[b] - nvar
    }

    # correction for conditional.x not categorical
    if (conditional.x && !categorical) {
      ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
      nexo <- length(ov.names.x)
      # add slopes
      ndat[b] <- ndat[b] + (nvar * nexo)
    }

    # correction for group proportions?
    group.idx <- which(partable$lhs == "group" &
      partable$op == "%" &
      partable$block == b)
    if (length(group.idx) > 0L) {
      # ndat <- ndat + (length(group.idx) - 1L) # G - 1 (sum to one)
      ndat[b] <- ndat[b] + 1L # poisson: each cell a parameter
    }
  } # blocks

  # sum over all blocks
  sum(ndat)
}

# total number of free parameters (ignoring equality constraints)
lav_partable_npar <- function(partable) {
  # we only assume non-zero values
  npar <- length(which(partable$free > 0L))
  npar
}

# global degrees of freedom: ndat - npar
# ignoring constraints! (not very useful)
#
# we need to find the rank of con.jac to find the exact amount
# of non-redundant equality constraints (this is done in lav_test.R)
lav_partable_df <- function(partable) {
  npar <- lav_partable_npar(partable)
  ndat <- lav_partable_ndat(partable)

  # degrees of freedom
  df <- ndat - npar

  as.integer(df)
}

# check order of covariances: we only fill the upper.tri
# therefore, we 'switch' lhs & rhs if they appear in the wrong order
lav_partable_covariance_reorder <- function(partable, # nolint
                                            ov.names = NULL,
                                            lv.names = NULL) {
  # shortcut
  cov.idx <- which(partable$op == "~~" & partable$lhs != partable$rhs)
  if (length(cov.idx) == 0L) {
    # nothing to do
    return(partable)
  }

  # get names
  if (is.null(ov.names)) {
    ov.names <- lav_partable_vnames(partable, "ov")
  } else {
    ov.names <- unlist(ov.names)
  }
  if (is.null(lv.names)) {
    lv.names <- lav_partable_vnames(partable, "lv")
    # add random slopes (if any)
    if (!is.null(partable$rv) && any(nchar(partable$rv) > 0L)) {
      rv.names <- unique(partable$rv[nchar(partable$rv) > 0L])
      lv.names <- c(lv.names, rv.names)
    }
  } else {
    lv.names <- unlist(lv.names)
  }
  lv.ov.names <- c(lv.names, ov.names)

  # identify wrong ordering
  lhs.idx <- match(partable$lhs[cov.idx], lv.ov.names)
  rhs.idx <- match(partable$rhs[cov.idx], lv.ov.names)
  swap.idx <- cov.idx[lhs.idx > rhs.idx]

  if (length(swap.idx) == 0L) {
    # nothing to do
    return(partable)
  }

  # swap!
  tmp <- partable$lhs[swap.idx]
  partable$lhs[swap.idx] <- partable$rhs[swap.idx]
  partable$rhs[swap.idx] <- tmp

  partable
}

# add a single parameter to an existing parameter table
lav_partable_add <- function(partable = NULL, add = list()) {
  # treat partable as list, not as a data.frame
  partable <- as.list(partable)

  # number of elements
  nel <- length(partable$lhs)

  # add copy of last row
  for (c in seq_len(length(partable))) {
    if (is.integer(partable[[c]][[1]])) {
      if (partable[[c]][nel] == 0L) {
        partable[[c]][nel + 1] <- 0L
      } else if (partable[[c]][nel] == 1L) {
        partable[[c]][nel + 1] <- 1L
      } else {
        partable[[c]][nel + 1] <- partable[[c]][nel] + 1L
      }
    } else if (is.character(partable[[c]][[1]])) {
      partable[[c]][nel + 1] <- ""
    } else if (is.numeric(partable[[c]][[1]])) {
      partable[[c]][nel + 1] <- 0
    } else {
      partable[[c]][nel + 1] <- partable[[c]][nel]
    }

    # replace
    if (names(partable)[c] %in% names(add)) {
      partable[[c]][nel + 1] <- add[[names(partable)[c]]]
    }
  }

  partable
}


# look for p2-row-idx of p1 elements
# p1 is usually a subset of p2
# return NA if not found
lav_partable_map_id_p1_in_p2 <- function(p1, p2, stopifnotfound = TRUE,
                                         exclude.nonpar = TRUE) {
  # check if we have a 'block' column (in both p1 and p2)
  if (is.null(p1$block)) {
    if (is.null(p1$group)) {
      p1$block <- rep.int(1L, length(p1$lhs))
    } else {
      p1$block <- p1$group
    }
  }
  if (is.null(p2$block)) {
    if (is.null(p2$group)) {
      p2$block <- rep.int(1L, length(p2$lhs))
    } else {
      p2$block <- p2$group
    }
  }

  # ALL rows from p1, or only 'parameters'?
  if (exclude.nonpar) {
    # get all parameters that have a '.p*' plabel
    # (they exclude "==", "<", ">", ":=")
    p1.idx <- which(grepl("\\.p", p1$plabel))
  } else {
    # all of it
    # note: block should be '0' in both p1 and p2
    p1.idx <- seq_len(length(p1$lhs))
  }
  np1 <- length(p1.idx)

  # return p2.id
  p2.id <- integer(np1)

  # check every parameter in p1
  for (i in seq_len(np1)) {
    # identify parameter in p1
    lhs <- p1$lhs[i]
    op <- p1$op[i]
    rhs <- p1$rhs[i]
    block <- p1$block[i]

    # search for corresponding parameter in p2
    p2.idx <- which(p2$lhs == lhs & p2$op == op & p2$rhs == rhs &
      p2$block == block)

    # found?
    if (length(p2.idx) == 0L) {
      if (stopifnotfound) {
        lav_msg_stop(gettext("parameter in p1 not found in p2:"),
          paste(lhs, op, rhs, "(block = ", block, ")", sep = " ")
        )
      } else {
        p2.id[i] <- as.integer(NA)
      }
    } else {
      p2.id[i] <- p2.idx
    }
  }

  p2.id
}
lav_partable_da2ovda <- function(partable) {
  # convert handling of ov.order = "data" with "da-operator elements"
  #         to "ovda attribute"
  if (any(partable$op == "da")) {
    da.idx <- which(partable$op == "da")
    ov.names.data <- partable$lhs[da.idx]
    temp <- lapply(partable, function(x) x[-da.idx])
    # names(temp) <- names(partable)
    attr(temp, "ovda") <- ov.names.data
    return(temp)
  }
  return(partable)
}
