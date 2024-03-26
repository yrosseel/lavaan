# create `full' parameter table, containing (almost) all parameters
# that could be free
#
# main motivation: univariate scores tests (modification indices)
#
lav_partable_full <- function(partable = NULL,
                              strict.exo = FALSE,
                              free = FALSE, start = FALSE) {
  # check minimum requirements: lhs, op, rhs
  stopifnot(
    !is.null(partable$lhs),
    !is.null(partable$op),
    !is.null(partable$rhs)
  )

  # lavpta?
  lavpta <- lav_partable_attributes(partable)

  # meanstructure
  if (!is.null(lavpta$meanstructure)) {
    meanstructure <- lavpta$meanstructure
  } else {
    # old object
    meanstructure <- any(partable$op == "~1")
  }

  # number of blocks
  nblocks <- lavpta$nblocks
  ngroups <- lavpta$ngroups
  nlevels <- lavpta$nlevels

  lhs <- rhs <- op <- character(0L)
  block <- group <- level <- integer(0L)

  # new in 0.6-3:
  GROUP.values <- lav_partable_group_values(partable)
  LEVEL.values <- lav_partable_level_values(partable)
  if (is.character(GROUP.values[1])) {
    group <- character(0L)
  }
  if (is.character(LEVEL.values[1L])) {
    level <- character(0L)
  }

  # block number
  b <- 0L
  for (g in 1:ngroups) {
    for (l in 1:nlevels) {
      # block
      b <- b + 1L

      ov.names <- lavpta$vnames$ov[[b]]
      ov.names.nox <- lavpta$vnames$ov.nox[[b]]
      ov.names.x <- lavpta$vnames$ov.x[[b]]
      ov.names.ind <- lavpta$vnames$ov.ind[[b]]
      ov.names.ord <- lavpta$vnames$ov.ord[[b]]

      lv.names <- lavpta$vnames$lv[[b]]

      # eqs.y, eqs.x
      eqs.names <- unique(c(
        lavpta$vnames$eqs.y[[b]],
        lavpta$vnames$eqs.x[[b]]
      ))
      if (length(eqs.names) > 0L) {
        eqs.y <- eqs.names
        if (strict.exo) {
          x.idx <- which(eqs.names %in% ov.names.x)
          if (length(x.idx) > 0L) {
            eqs.y <- eqs.names[-x.idx]
          }
        }
        eqs.x <- eqs.names
      } else {
        eqs.y <- character(0L)
        eqs.x <- character(0L)
      }

      # 1 "=~"
      l.lhs <- rep(lv.names, each = length(ov.names.nox))
      l.rhs <- rep(ov.names.nox, times = length(lv.names))

      # remove factor ~ eqs.y combinations, if any
      # because they also appear as a regression
      bad.idx <- which(l.lhs %in% lv.names &
        l.rhs %in% eqs.y)
      if (length(bad.idx) > 0L) {
        l.lhs <- l.lhs[-bad.idx]
        l.rhs <- l.rhs[-bad.idx]
      }
      l.op <- rep("=~", length(l.lhs))

      # 2a. "~~" ov ## FIXME: ov.names.nox or ov.names??
      # if(strict.exo) {
      OV <- ov.names.nox
      # } else {
      #    OV <- ov.names
      # }
      nx <- length(OV)
      idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
      ov.lhs <- rep(OV, each = nx)[idx] # fill upper.tri
      ov.rhs <- rep(OV, times = nx)[idx]
      ov.op <- rep("~~", length(ov.lhs))

      # remove dummy indicators that correlate with 'proper'
      # indicators; new in 0.6-14; fixed in 0.6-16
      ov.other <- ov.names[!ov.names %in% c(
        ov.names.ind, ov.names.x,
        eqs.x, eqs.y
      )]
      if (length(ov.other) > 0L) {
        bad.idx <- which((ov.lhs %in% ov.names &
          ov.rhs %in% ov.other) |
          (ov.lhs %in% ov.other &
            ov.rhs %in% ov.names))
        if (length(bad.idx) > 0L) {
          ov.lhs <- ov.lhs[-bad.idx]
          ov.rhs <- ov.rhs[-bad.idx]
          ov.op <- ov.op[-bad.idx]
        }
      }

      # exo ~~
      if (!strict.exo && length(ov.names.x) > 0L) {
        OV <- ov.names.x
        nx <- length(OV)
        idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
        more.lhs <- rep(OV, each = nx)[idx] # fill upper.tri
        more.rhs <- rep(OV, times = nx)[idx]
        ov.lhs <- c(ov.lhs, more.lhs)
        ov.rhs <- c(ov.rhs, more.rhs)
        ov.op <- c(ov.op, rep("~~", length(more.lhs)))
      }

      # 2b. "~~" lv
      nx <- length(lv.names)
      idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
      lv.lhs <- rep(lv.names, each = nx)[idx] # fill upper.tri
      lv.rhs <- rep(lv.names, times = nx)[idx]
      lv.op <- rep("~~", length(lv.lhs))

      # 3 regressions?
      r.lhs <- r.rhs <- r.op <- character(0)
      if (length(eqs.names) > 0L) {
        r.lhs <- rep(eqs.y, each = length(eqs.x))
        r.rhs <- rep(eqs.x, times = length(eqs.y))

        # remove self-arrows
        idx <- which(r.lhs == r.rhs)
        if (length(idx) > 0L) {
          r.lhs <- r.lhs[-idx]
          r.rhs <- r.rhs[-idx]
        }

        # remove indicator ~ factor if they exist
        bad.idx <- which(r.lhs %in% ov.names.ind &
          r.rhs %in% lv.names)
        if (length(bad.idx) > 0L) {
          r.lhs <- r.lhs[-bad.idx]
          r.rhs <- r.rhs[-bad.idx]
        }

        r.op <- rep("~", length(r.rhs))
      }

      # 4. intercepts
      int.lhs <- int.rhs <- int.op <- character(0)
      if (meanstructure) {
        if (strict.exo) {
          int.lhs <- c(ov.names.nox, lv.names)
        } else {
          int.lhs <- c(ov.names, lv.names)
        }
        int.rhs <- rep("", length(int.lhs))
        int.op <- rep("~1", length(int.lhs))
      }

      # 5. thresholds
      th.lhs <- th.rhs <- th.op <- character(0)
      if (length(ov.names.ord) > 0L) {
        th.names <- lavpta$vnames$th[[b]]
        tmp <- strsplit(th.names, "\\|")
        th.lhs <- sapply(tmp, function(x) x[1])
        th.rhs <- sapply(tmp, function(x) x[2])
        th.op <- rep("|", length(th.lhs))
      }

      # 6. scaling parameters
      delta.lhs <- delta.rhs <- delta.op <- character(0)
      if (ngroups > 1L && length(ov.names.ord) > 0L) {
        delta.lhs <- ov.names.ord
        delta.rhs <- ov.names.ord
        delta.op <- rep("~*~", length(delta.lhs))
      }

      # combine
      this.lhs <- c(
        l.lhs, ov.lhs, lv.lhs, r.lhs, int.lhs, th.lhs,
        delta.lhs
      )
      this.rhs <- c(
        l.rhs, ov.rhs, lv.rhs, r.rhs, int.rhs, th.rhs,
        delta.rhs
      )
      this.op <- c(
        l.op, ov.op, lv.op, r.op, int.op, th.op,
        delta.op
      )
      n.el <- length(this.lhs)

      lhs <- c(lhs, this.lhs)
      rhs <- c(rhs, this.rhs)
      op <- c(op, this.op)
      block <- c(block, rep(b, n.el))
      group <- c(group, rep(GROUP.values[g], n.el))
      level <- c(level, rep(LEVEL.values[l], n.el))
    } # level
  } # group

  LIST <- data.frame(
    lhs = lhs, op = op, rhs = rhs, block = block,
    group = group, level = level, stringsAsFactors = FALSE
  )

  if (free) {
    LIST$free <- rep(0L, nrow(LIST))
  }

  if (start) {
    LIST$start <- rep(0, nrow(LIST))
  }

  LIST
}
