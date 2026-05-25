
# handle EFA equality constraints
# YR 14 Jan 2020: 0.6-6 does no longer impose 'explicit' constraints
#                 if we only need to fix a parameter to 0/1
# Note: we should also check if they are really needed:
#       eg., if all the factor-loadings of the 'second' set (time/group)
#       are constrained to be equal to the factor-loadings of the first
#       set, no further constraints are needed
lav_pt_efa_con <- function(list_1 = NULL,
                                         orthogonal_efa = FALSE,
                                         group_equal = character(0L)) {
  # for each set, for each block
  nblocks <- lav_pt_nblocks(list_1)
  set_names <- lav_pt_efa_values(list_1)
  nsets <- length(set_names)

  for (b in seq_len(nblocks)) {
    for (s in seq_len(nsets)) {
      # lv's for this block/set
      lv_nam_efa <- unique(list_1$lhs[list_1$op == "=~" &
        list_1$block == b &
        list_1$efa == set_names[s]])
      if (length(lv_nam_efa) == 1L) {
        # nothing to do (warn?)
        next
      }

      # equality constraints on ALL factor loadings in this set?
      # two scenario's:
      # 1. eq constraints within the same block, perhaps time1/time2/
      # 2. eq constraints across groups (group.equal = "loadings")
      # --> no constraints are needed

      # store labels (if any)
      fix_to_zero <- TRUE

      # 1. within block/group
      if (s == 1L) {
        set_idx <- which(list_1$op == "=~" &
          list_1$block == b &
          list_1$lhs %in% lv_nam_efa)
        label_set1 <- list_1$label[set_idx]
      } else {
        # collect labels for this set, if any
        set_idx <- which(list_1$op == "=~" &
          list_1$block == b &
          list_1$lhs %in% lv_nam_efa)

        # user-provided labels (if any)
        this_label_set <- list_1$label[set_idx]

        # same as in reference set?
        if (all(nchar(this_label_set) > 0L) &&
          all(this_label_set %in% label_set1)) {
          fix_to_zero <- FALSE
        }
      }

      # 2. across groups
      if (b == 1L) {
        set_idx <- which(list_1$op == "=~" &
          list_1$block == b &
          list_1$lhs %in% lv_nam_efa)
        label_group1 <- list_1$label[set_idx]
      } else {
        if ("loadings" %in% group_equal) {
          fix_to_zero <- FALSE
        } else {
          # collect labels for this set, if any
          set_idx <- which(list_1$op == "=~" &
            list_1$block == b &
            list_1$lhs %in% lv_nam_efa)

          # user-provided labels (if any)
          this_label_set <- list_1$label[set_idx]

          # same as in reference set?
          if (all(nchar(this_label_set) > 0L) &&
            all(this_label_set %in% label_group1)) {
            fix_to_zero <- FALSE
          }
        }
      }

      # 1. echelon pattern
      nfac <- length(lv_nam_efa)
      for (f in seq_len(nfac)) {
        if (f == 1L) {
          next
        }
        nzero <- (f - 1L)
        ind_idx <- which(list_1$op == "=~" &
          list_1$block == b &
          list_1$lhs %in% lv_nam_efa[f])
        if (length(ind_idx) < nzero) {
          lav_msg_stop(gettextf(
            "efa factor %s has not enough indicators for echelon pattern",
            lv_nam_efa[f]))
        }

        # fix to zero
        if (fix_to_zero) {
          list_1$free[ind_idx[seq_len(nzero)]] <- 0L
          list_1$ustart[ind_idx[seq_len(nzero)]] <- 0
          list_1$user[ind_idx[seq_len(nzero)]] <- 7L
        } else {
          list_1$user[ind_idx[seq_len(nzero)]] <- 77L
        }
      }

      # 2. covariances constrained to zero (only if oblique rotation)
      if (!orthogonal_efa) {
        # skip if user == 1 (user-override!)
        cov_idx <- which(list_1$op == "~~" &
          list_1$block == b &
          list_1$user == 0L &
          list_1$lhs %in% lv_nam_efa &
          list_1$rhs %in% lv_nam_efa &
          list_1$lhs != list_1$rhs)

        # fix to zero
        if (fix_to_zero) {
          list_1$free[cov_idx] <- 0L
          list_1$ustart[cov_idx] <- 0
          list_1$user[cov_idx] <- 7L
        } else {
          list_1$user[cov_idx] <- 77L
        }
      }
    } # sets
  } # blocks

  list_1
}
