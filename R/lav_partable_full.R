# create `full' parameter table, containing (almost) all parameters
# that could be free
#
# main motivation: univariate scores tests (modification indices)
#
lav_partable_full <- function(partable = NULL,
                              strict_exo = FALSE,
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
  # nblocks <- lavpta$nblocks
  ngroups <- lavpta$ngroups
  nlevels <- lavpta$nlevels

  lhs <- rhs <- op <- character(0L)
  block <- group <- level <- integer(0L)

  # new in 0.6-3:
  group_values <- lav_partable_group_values(partable)
  level_values <- lav_partable_level_values(partable)
  if (is.character(group_values[1])) {
    group <- character(0L)
  }
  if (is.character(level_values[1L])) {
    level <- character(0L)
  }

  # block number
  b <- 0L
  for (g in 1:ngroups) {
    for (l in 1:nlevels) {
      # block
      b <- b + 1L

      ov_names <- lavpta$vnames$ov[[b]]
      ov_names_nox <- lavpta$vnames$ov.nox[[b]]
      ov_names_x <- lavpta$vnames$ov.x[[b]]
      ov_names_ind <- lavpta$vnames$ov.ind[[b]]
      ov_names_ord <- lavpta$vnames$ov.ord[[b]]

      lv_names <- lavpta$vnames$lv[[b]]

      # eqs.y, eqs.x
      eqs_names <- unique(c(
        lavpta$vnames$eqs.y[[b]],
        lavpta$vnames$eqs.x[[b]]
      ))
      if (length(eqs_names) > 0L) {
        eqs_y <- eqs_names
        if (strict_exo) {
          x_idx <- which(eqs_names %in% ov_names_x)
          if (length(x_idx) > 0L) {
            eqs_y <- eqs_names[-x_idx]
          }
        }
        eqs_x <- eqs_names
      } else {
        eqs_y <- character(0L)
        eqs_x <- character(0L)
      }

      # 1 "=~"
      l_lhs <- rep(lv_names, each = length(ov_names_nox))
      l_rhs <- rep(ov_names_nox, times = length(lv_names))

      # remove factor ~ eqs.y combinations, if any
      # because they also appear as a regression
      bad_idx <- which(l_lhs %in% lv_names &
        l_rhs %in% eqs_y)
      if (length(bad_idx) > 0L) {
        l_lhs <- l_lhs[-bad_idx]
        l_rhs <- l_rhs[-bad_idx]
      }
      l_op <- rep("=~", length(l_lhs))

      # 2a. "~~" ov ## FIXME: ov.names.nox or ov.names??
      # if(strict.exo) {
      ov <- ov_names_nox
      # } else {
      #    OV <- ov.names
      # }
      nx <- length(ov)
      idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
      ov_lhs <- rep(ov, each = nx)[idx] # fill upper.tri
      ov_rhs <- rep(ov, times = nx)[idx]
      ov_op <- rep("~~", length(ov_lhs))

      # remove dummy indicators that correlate with 'proper'
      # indicators; new in 0.6-14; fixed in 0.6-16
      ov_other <- ov_names[!ov_names %in% c(
        ov_names_ind, ov_names_x,
        eqs_x, eqs_y
      )]
      if (length(ov_other) > 0L) {
        bad_idx <- which((ov_lhs %in% ov_names &
          ov_rhs %in% ov_other) |
          (ov_lhs %in% ov_other &
            ov_rhs %in% ov_names))
        if (length(bad_idx) > 0L) {
          ov_lhs <- ov_lhs[-bad_idx]
          ov_rhs <- ov_rhs[-bad_idx]
          ov_op <- ov_op[-bad_idx]
        }
      }

      # exo ~~
      if (!strict_exo && length(ov_names_x) > 0L) {
        ov <- ov_names_x
        nx <- length(ov)
        idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
        more_lhs <- rep(ov, each = nx)[idx] # fill upper.tri
        more_rhs <- rep(ov, times = nx)[idx]
        ov_lhs <- c(ov_lhs, more_lhs)
        ov_rhs <- c(ov_rhs, more_rhs)
        ov_op <- c(ov_op, rep("~~", length(more_lhs)))
      }

      # 2b. "~~" lv
      nx <- length(lv_names)
      idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
      lv_lhs <- rep(lv_names, each = nx)[idx] # fill upper.tri
      lv_rhs <- rep(lv_names, times = nx)[idx]
      lv_op <- rep("~~", length(lv_lhs))

      # 3 regressions?
      r_lhs <- r_rhs <- r_op <- character(0)
      if (length(eqs_names) > 0L) {
        r_lhs <- rep(eqs_y, each = length(eqs_x))
        r_rhs <- rep(eqs_x, times = length(eqs_y))

        # remove self-arrows
        idx <- which(r_lhs == r_rhs)
        if (length(idx) > 0L) {
          r_lhs <- r_lhs[-idx]
          r_rhs <- r_rhs[-idx]
        }

        # remove indicator ~ factor if they exist
        bad_idx <- which(r_lhs %in% ov_names_ind &
          r_rhs %in% lv_names)
        if (length(bad_idx) > 0L) {
          r_lhs <- r_lhs[-bad_idx]
          r_rhs <- r_rhs[-bad_idx]
        }

        r_op <- rep("~", length(r_rhs))
      }

      # 4. intercepts
      int_lhs <- int_rhs <- int_op <- character(0)
      if (meanstructure) {
        if (strict_exo) {
          int_lhs <- c(ov_names_nox, lv_names)
        } else {
          int_lhs <- c(ov_names, lv_names)
        }
        int_rhs <- rep("", length(int_lhs))
        int_op <- rep("~1", length(int_lhs))
      }

      # 5. thresholds
      th_lhs <- th_rhs <- th_op <- character(0)
      if (length(ov_names_ord) > 0L) {
        th_names <- lavpta$vnames$th[[b]]
        tmp <- strsplit(th_names, "\\|")
        th_lhs <- sapply(tmp, function(x) x[1])
        th_rhs <- sapply(tmp, function(x) x[2])
        th_op <- rep("|", length(th_lhs))
      }

      # 6. scaling parameters
      delta_lhs <- delta_rhs <- delta_op <- character(0)
      if (ngroups > 1L && length(ov_names_ord) > 0L) {
        delta_lhs <- ov_names_ord
        delta_rhs <- ov_names_ord
        delta_op <- rep("~*~", length(delta_lhs))
      }

      # combine
      this_lhs <- c(
        l_lhs, ov_lhs, lv_lhs, r_lhs, int_lhs, th_lhs,
        delta_lhs
      )
      this_rhs <- c(
        l_rhs, ov_rhs, lv_rhs, r_rhs, int_rhs, th_rhs,
        delta_rhs
      )
      this_op <- c(
        l_op, ov_op, lv_op, r_op, int_op, th_op,
        delta_op
      )
      n_el <- length(this_lhs)

      lhs <- c(lhs, this_lhs)
      rhs <- c(rhs, this_rhs)
      op <- c(op, this_op)
      block <- c(block, rep(b, n_el))
      group <- c(group, rep(group_values[g], n_el))
      level <- c(level, rep(level_values[l], n_el))
    } # level
  } # group

  list_1 <- data.frame(
    lhs = lhs, op = op, rhs = rhs, block = block,
    group = group, level = level, stringsAsFactors = FALSE
  )

  if (free) {
    list_1$free <- rep(0L, nrow(list_1))
  }

  if (start) {
    list_1$start <- rep(0, nrow(list_1))
  }

  list_1
}
