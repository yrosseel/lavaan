
# handle EFA equality constraints
# YR 14 Jan 2020: 0.6-6 does no longer impose 'explicit' constraints
#                 if we only need to fix a parameter to 0/1
# Note: we should also check if they are really needed:
#       eg., if all the factor-loadings of the 'second' set (time/group)
#       are constrained to be equal to the factor-loadings of the first
#       set, no further constraints are needed
lav_partable_efa_constraints <- function(LIST = NULL,
                                         orthogonal.efa = FALSE,
                                         group.equal = character(0L)) {
  # for each set, for each block
  nblocks <- lav_partable_nblocks(LIST)
  set.names <- lav_partable_efa_values(LIST)
  nsets <- length(set.names)

  for (b in seq_len(nblocks)) {
    for (s in seq_len(nsets)) {
      # lv's for this block/set
      lv.nam.efa <- unique(LIST$lhs[LIST$op == "=~" &
        LIST$block == b &
        LIST$efa == set.names[s]])
      if (length(lv.nam.efa) == 1L) {
        # nothing to do (warn?)
        next
      }

      # equality constraints on ALL factor loadings in this set?
      # two scenario's:
      # 1. eq constraints within the same block, perhaps time1/time2/
      # 2. eq constraints across groups (group.equal = "loadings")
      # --> no constraints are needed

      # store labels (if any)
      fix.to.zero <- TRUE

      # 1. within block/group
      if (s == 1L) {
        set.idx <- which(LIST$op == "=~" &
          LIST$block == b &
          LIST$lhs %in% lv.nam.efa)
        LABEL.set1 <- LIST$label[set.idx]
      } else {
        # collect labels for this set, if any
        set.idx <- which(LIST$op == "=~" &
          LIST$block == b &
          LIST$lhs %in% lv.nam.efa)

        # user-provided labels (if any)
        this.label.set <- LIST$label[set.idx]

        # same as in reference set?
        if (all(nchar(this.label.set) > 0L) &&
          all(this.label.set %in% LABEL.set1)) {
          fix.to.zero <- FALSE
        }
      }

      # 2. across groups
      if (b == 1L) {
        set.idx <- which(LIST$op == "=~" &
          LIST$block == b &
          LIST$lhs %in% lv.nam.efa)
        LABEL.group1 <- LIST$label[set.idx]
      } else {
        if ("loadings" %in% group.equal) {
          fix.to.zero <- FALSE
        } else {
          # collect labels for this set, if any
          set.idx <- which(LIST$op == "=~" &
            LIST$block == b &
            LIST$lhs %in% lv.nam.efa)

          # user-provided labels (if any)
          this.label.set <- LIST$label[set.idx]

          # same as in reference set?
          if (all(nchar(this.label.set) > 0L) &&
            all(this.label.set %in% LABEL.group1)) {
            fix.to.zero <- FALSE
          }
        }
      }

      # 1. echelon pattern
      nfac <- length(lv.nam.efa)
      for (f in seq_len(nfac)) {
        if (f == 1L) {
          next
        }
        nzero <- (f - 1L)
        ind.idx <- which(LIST$op == "=~" &
          LIST$block == b &
          LIST$lhs %in% lv.nam.efa[f])
        if (length(ind.idx) < nzero) {
          lav_msg_stop(gettextf(
            "efa factor %s has not enough indicators for echelon pattern",
            lv.nam.efa[f]))
        }

        # fix to zero
        if (fix.to.zero) {
          LIST$free[ind.idx[seq_len(nzero)]] <- 0L
          LIST$ustart[ind.idx[seq_len(nzero)]] <- 0
          LIST$user[ind.idx[seq_len(nzero)]] <- 7L
        } else {
          LIST$user[ind.idx[seq_len(nzero)]] <- 77L
        }
      }

      # 2. covariances constrained to zero (only if oblique rotation)
      if (!orthogonal.efa) {
        # skip if user == 1 (user-override!)
        cov.idx <- which(LIST$op == "~~" &
          LIST$block == b &
          LIST$user == 0L &
          LIST$lhs %in% lv.nam.efa &
          LIST$rhs %in% lv.nam.efa &
          LIST$lhs != LIST$rhs)

        # fix to zero
        if (fix.to.zero) {
          LIST$free[cov.idx] <- 0L
          LIST$ustart[cov.idx] <- 0
          LIST$user[cov.idx] <- 7L
        } else {
          LIST$user[cov.idx] <- 77L
        }
      }
    } # sets
  } # blocks

  LIST
}
