
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

  # group.equal = "loadings": the loadings of the efa factors are shared
  # across groups, and a single rotation (based on the first group) is used
  # for all the groups. The model is only invariant to this rotation if the
  # factor covariance matrix of the other groups is either completely free
  # (the default), or completely equal to the one of the first group (both
  # "lv.variances" and "lv.covariances" in group.equal)
  group_values <- unique(list_1$group[list_1$group > 0L])
  mg_loadings_flag <- (length(group_values) > 1L &&
    "loadings" %in% group_equal)
  free_lv_var_flag <- !("lv.variances" %in% group_equal)
  free_lv_cov_flag <- !("lv.covariances" %in% group_equal)
  if (mg_loadings_flag && free_lv_var_flag != free_lv_cov_flag &&
    nsets > 0L) {
    lav_msg_warn(gettext(
      "group_equal contains \"loadings\" and only one of \"lv.variances\" and
      \"lv.covariances\": the model for the efa factors is not rotation
      invariant, and the results will depend on the (arbitrary) unrotated
      solution. Use either both or none of them."))
  }

  for (b in seq_len(nblocks)) {
    # is this block part of the first group?
    first_group_flag <- TRUE
    if (length(group_values) > 1L) {
      first_group_flag <-
        (list_1$group[list_1$block == b][1L] == group_values[1L])
    }
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

      # 3. group.equal = "loadings": free the factor variances and
      #    covariances in all but the first group (also if the rotation is
      #    orthogonal: only the factors of the first group are orthogonal),
      #    unless they are constrained to be equal across groups
      #    (skip if user == 1: user-override!)
      if (mg_loadings_flag && !first_group_flag) {
        if (free_lv_var_flag) {
          var_idx <- which(list_1$op == "~~" &
            list_1$block == b &
            list_1$user == 0L &
            list_1$lhs %in% lv_nam_efa &
            list_1$lhs == list_1$rhs)
          list_1$free[var_idx] <- 1L
          list_1$ustart[var_idx] <- as.numeric(NA)
        }
        if (free_lv_cov_flag) {
          cov_idx <- which(list_1$op == "~~" &
            list_1$block == b &
            list_1$user %in% c(0L, 77L) &
            list_1$lhs %in% lv_nam_efa &
            list_1$rhs %in% lv_nam_efa &
            list_1$lhs != list_1$rhs)
          list_1$free[cov_idx] <- 1L
          list_1$ustart[cov_idx] <- as.numeric(NA)
          list_1$user[cov_idx] <- 77L
        }
      }
    } # sets
  } # blocks

  list_1
}
