# YR 11 feb 2017: initial version

# given a parameter table (PT), extract a part of the model:
# eg.:
# - only the measurement model (with saturated latent variables)
# - only the structural part
# - a single measurement block
# ...

# YR 25 June 2021: - add.exo.cov = TRUE for structural model
#                  - fixed.x = FALSE/TRUE -> exo flags


# FIXME:
# - but fixed-to-zero covariances may not be present in PT...
# - if indicators are regressed on exogenous covariates, should we
#   add them here? (no for now, unless add.ind.predictors = TRUE)
# - new in 0.6-20: - check for 2nd, 3rd order lv.names...
#                  - allow for conditional.x (global SAM)
lav_partable_subset_measurement_model <- function(pt_1 = NULL,       # nolint
                                                  lv_names = NULL,
                                                  add_lv_cov = TRUE,
                                                  add_ind_predictors = FALSE,
                                                  add_idx = FALSE,
                                                  idx_only = FALSE) {
  # PT
  pt_1 <- as.data.frame(pt_1, stringsAsFactors = FALSE)

  # check if we have free.unrotated column (rotation!)
  if (!is.null(pt_1$free.unrotated)) {
    pt_1$free.orig <- pt_1$free
    pt_1$free <- pt_1$free.unrotated
    pt_1$est.unrotated <- NULL
    pt_1$free.unrotated <- NULL
    pt_1$est.std <- NULL
  }

  # lavpta
  lavpta <- lav_partable_attributes(pt_1)

  # nblocks
  nblocks <- lavpta$nblocks
  block_values <- lav_partable_block_values(pt_1)

  # lv.names: list with element per block
  if (is.null(lv_names)) {
    lv_names <- lavpta$vnames$lv.regular
  } else if (!is.list(lv_names)) {
    lv_names <- rep(list(lv_names), nblocks)
  }

  # if lv.names contains a higher-order latent variable,
  # add its (latent) indicators
  for (g in 1:nblocks) {
    if (length(lavpta$vnames$lv.ind[[g]]) == 0L) {
      next
    }
    # ALL lv names
    lv_names_1 <- lavpta$vnames$lv[[g]]

    # lv.names in this block
    this_lv <- lv_names[[g]]

    all_ind_are_observed_flag <- FALSE
    new_lv <- character(0L)
    while (!all_ind_are_observed_flag) {
      # check indicators of all this.lv
      rhs <- pt_1$rhs[which(pt_1$op == "=~" & pt_1$block == g &
                          pt_1$lhs %in% this_lv)]
      lv_idx <- which(rhs %in% lv_names_1)
      if (length(lv_idx) > 0L) {
        # add these 'new' ones to new.lv
        new_lv <- c(new_lv, rhs[lv_idx])
        this_lv <- rhs[lv_idx]
      } else {
        all_ind_are_observed_flag <- TRUE
      }
    }

    # update lv.names for this block
    lv_names[[g]] <- unique(c(lv_names[[g]], new_lv))
  }

  # keep rows idx
  keep_idx <- integer(0L)

  # remove not-needed measurement models
  for (g in 1:nblocks) {
    # indicators for latent variables we keep
    ind_idx <- which(pt_1$op == "=~" &
      pt_1$lhs %in% lv_names[[g]] &
      pt_1$block == block_values[g])
    ind <- pt_1$rhs[ind_idx]
    # ind_plabel <- pt_1$plabel[ind_idx]

    # keep =~
    keep_idx <- c(keep_idx, ind_idx)

    # new in 0.6-17: indicators regressed on predictors
    if (add_ind_predictors) {
      pred_idx <- which(pt_1$op == "~" &
        pt_1$lhs %in% ind &
        pt_1$block == block_values[g])
      extra <- unique(pt_1$rhs[pred_idx])
      keep_idx <- c(keep_idx, pred_idx)
      # add them to IND, so we include their variances/intercepts
      ind <- c(ind, extra)
    }

    # keep ~~
    ov_var_idx <- which(pt_1$op == "~~" &
      pt_1$lhs %in% ind &
      pt_1$rhs %in% ind &
      pt_1$block == block_values[g])
    keep_idx <- c(keep_idx, ov_var_idx)

    lv_var_idx <- which(pt_1$op == "~~" &
      pt_1$lhs %in% lv_names[[g]] &
      pt_1$rhs %in% lv_names[[g]] &
      pt_1$block == block_values[g])
    keep_idx <- c(keep_idx, lv_var_idx)

    # intercepts indicators
    ov_int_idx <- which(pt_1$op == "~1" &
      pt_1$lhs %in% ind &
      pt_1$block == block_values[g])
    keep_idx <- c(keep_idx, ov_int_idx)

    # intercepts latent variables
    lv_int_idx <- which(pt_1$op == "~1" &
      pt_1$lhs %in% lv_names[[g]] &
      pt_1$block == block_values[g])
    keep_idx <- c(keep_idx, lv_int_idx)

    # thresholds
    th_idx <- which(pt_1$op == "|" &
      pt_1$lhs %in% ind &
      pt_1$block == block_values[g])
    keep_idx <- c(keep_idx, th_idx)

    # scaling factors
    sc_idx <- which(pt_1$op == "~*~" &
      pt_1$lhs %in% ind &
      pt_1$block == block_values[g])
    keep_idx <- c(keep_idx, sc_idx)

    # defined/constraints
    if (any(pt_1$op %in% c("==", "<", ">", ":="))) {
      # get the 'id' numbers and the labels involved in def/constraints
      pt2 <- pt_1
      pt2$free <- pt_1$id # use 'id' numbers instead of 'free' indices
      id <- lav_partable_constraints_label_id(pt2, def = TRUE)
      label <- names(id)

      # what are the row indices that we currently keep?
      free_id <- pt_1$id[keep_idx]
    }

    # defined parameters
    def_idx <- which(pt_1$op == ":=")
    if (length(def_idx) > 0L) {
      def_keep <- logical(length(def_idx))
      for (def in seq_along(def_idx)) {
        # rhs
        rhs_labels <- all.vars(as.formula(paste(
          "~",
          pt_1[def_idx[def], "rhs"]
        )))
        if (length(rhs_labels) > 0L) {
          # par id
          rhs_freeid <- id[match(rhs_labels, label)]

          # keep?
          if (all(rhs_freeid %in% free_id)) {
            def_keep[def] <- TRUE
          }
        } else { # only constants?
          def_keep[def] <- TRUE
        }
      }
      keep_idx <- c(keep_idx, def_idx[def_keep])
      # add 'id' numbers of := definitions that we keep
      free_id <- c(free_id, pt_1$id[def_idx[def_keep]])
    }

    # (in)equality constraints
    con_idx <- which(pt_1$op %in% c("==", "<", ">"))
    if (length(con_idx) > 0L) {
      con_keep <- logical(length(con_idx))
      for (con in seq_along(con_idx)) {
        lhs_keep <- FALSE
        rhs_keep <- FALSE

        # lhs
        lhs_labels <- all.vars(as.formula(paste(
          "~",
          pt_1[con_idx[con], "lhs"]
        )))
        if (length(lhs_labels) > 0L) {
          # par id
          lhs_freeid <- id[match(lhs_labels, label)]

          # keep?
          if (all(lhs_freeid %in% free_id)) {
            lhs_keep <- TRUE
          }
        } else {
          lhs_keep <- TRUE
        }

        # rhs
        rhs_labels <- all.vars(as.formula(paste(
          "~",
          pt_1[con_idx[con], "rhs"]
        )))
        if (length(rhs_labels) > 0L) {
          # par id
          rhs_freeid <- id[match(rhs_labels, label)]

          # keep?
          if (all(rhs_freeid %in% free_id)) {
            rhs_keep <- TRUE
          }
        } else {
          rhs_keep <- TRUE
        }

        if (lhs_keep && rhs_keep) {
          con_keep[con] <- TRUE
        }
      }

      keep_idx <- c(keep_idx, con_idx[con_keep])
    } # con
  } # block

  # remove any duplicated (only with higher-order factors)
  keep_idx <- keep_idx[!duplicated(keep_idx)]

  if (idx_only) {
    return(keep_idx)
  }

  pt_1 <- pt_1[keep_idx, , drop = FALSE]

  # check if we have enough indicators?
  # TODO

  # add covariances among latent variables?
  if (add_lv_cov) {
    pt_1 <- lav_partable_add_lv_cov(
      pt_1 = pt_1,
      lv_names = lv_names
    )
  }

  # clean up
  pt_1 <- lav_partable_complete(pt_1)

  if (add_idx) {
    attr(pt_1, "idx") <- keep_idx
  }

  pt_1
}

# NOTE: only within same level
lav_partable_add_lv_cov <- function(pt_1, lv_names = NULL) {
  # PT
  if (!is.data.frame(pt_1)) {
    pt_1 <- as.data.frame(pt_1, stringsAsFactors = FALSE)
  }

  # lavpta
  lavpta <- lav_partable_attributes(pt_1)

  # nblocks
  nblocks <- lavpta$nblocks
  block_values <- lav_partable_block_values(pt_1)

  # lv.names: list with element per block
  if (is.null(lv_names)) {
    lv_names <- lavpta$vnames$lv.regular
  } else if (!is.list(lv_names)) {
    lv_names <- rep(list(lv_names), nblocks)
  }

  # check for higher-order models (new in 0.6-20)
  for (b in seq_len(nblocks)) {
    if (length(lavpta$vnames$lv.ind[[b]]) > 0L) {
      ind_idx <- which(lv_names[[b]] %in% lavpta$vnames$lv.ind[[b]])
      if (length(ind_idx) > 0L) {
        lv_names[[b]] <- lv_names[[b]][-ind_idx]
      }
    }
  }

  # remove lv.names if not present at same level/block
  if (nblocks > 1L) {
    for (b in seq_len(nblocks)) {
      rm_idx <- which(!lv_names[[b]] %in% lavpta$vnames$lv.regular[[b]])
      if (length(rm_idx) > 0L) {
        lv_names[[b]] <- lv_names[[b]][-rm_idx]
      }
    } # b
  }

  # add covariances among latent variables
  for (b in seq_len(nblocks)) {
    if (length(lv_names[[b]]) > 1L) {
      tmp <- utils::combn(lv_names[[b]], 2L)
      for (i in seq_len(ncol(tmp))) {
        # already present?
        cov1_idx <- which(pt_1$op == "~~" &
          pt_1$block == block_values[b] &
          pt_1$lhs == tmp[1, i] & pt_1$rhs == tmp[2, i])
        cov2_idx <- which(pt_1$op == "~~" &
          pt_1$block == block_values[b] &
          pt_1$lhs == tmp[2, i] & pt_1$rhs == tmp[1, i])

        # if not, add
        if (length(c(cov1_idx, cov2_idx)) == 0L) {
          add <- list(
            lhs = tmp[1, i],
            op = "~~",
            rhs = tmp[2, i],
            user = 3L,
            free = max(pt_1$free) + 1L,
            block = b
          )
          # add group column
          if (!is.null(pt_1$group)) {
            add$group <- unique(pt_1$block[pt_1$block == b])
          }
          # add level column
          if (!is.null(pt_1$level)) {
            add$level <- unique(pt_1$level[pt_1$block == b])
          }
          # add lower column
          if (!is.null(pt_1$lower)) {
            add$lower <- as.numeric(-Inf)
          }
          # add upper column
          if (!is.null(pt_1$upper)) {
            add$upper <- as.numeric(+Inf)
          }
          pt_1 <- lav_partable_add(pt_1, add = add)
        }
      }
    } # lv.names
  } # blocks

  pt_1
}

# this function takes a 'full' SEM (measurement models + structural part)
# and returns only the structural part
#
# - what to do if we have no regressions among the latent variables?
#   -> we return all covariances among the latent variables
#
lav_partable_subset_structural_model <- function(pt_1 = NULL,          # nolint
                                                 add_idx = FALSE,
                                                 idx_only = FALSE,
                                                 add_exo_cov = FALSE,
                                                 fixed_x = FALSE,
                                                 conditional_x = FALSE,
                                                 free_fixed_var = FALSE,
                                                 meanstructure = FALSE) {
  # PT
  pt_1 <- as.data.frame(pt_1, stringsAsFactors = FALSE)

  # remove any EFA related information -- new in 0.6-18
  if (!is.null(pt_1$efa)) {
    pt_1$efa <- NULL
    pt_1$est.unrotated <- NULL
    seven_idx <- which(pt_1$user == 7L & pt_1$op == "~~")
    if (length(seven_idx) > 0L) {
      pt_1$user[seven_idx] <- 0L
      pt_1$free[seven_idx] <- 1L
      pt_1$ustart[seven_idx] <- as.numeric(NA)
      pt_1$est[seven_idx] <- pt_1$est.std[seven_idx]
    }
    pt_1$est.std <- NULL
  }

  # lavpta
  lavpta <- lav_partable_attributes(pt_1)

  # nblocks
  nblocks <- lavpta$nblocks
  block_values <- lav_partable_block_values(pt_1)

  # eqs.names
  # eqs_x_names <- lavpta$vnames$eqs.x
  # eqs_y_names <- lavpta$vnames$eqs.y
  lv_names <- lavpta$vnames$lv.regular

  # keep rows idx
  keep_idx <- integer(0L)

  # remove not-needed measurement models
  for (g in 1:nblocks) {

    # eqs.names
    eqs_names <- unique(c(
      lavpta$vnames$eqs.x[[g]],
      lavpta$vnames$eqs.y[[g]]
    ))
    if (length(eqs_names) == 0L) { # no structural model
      eqs_names <- lv_names[[g]]
    }
    # all.names <- unique(c(
    #   eqs.names,
    #   lavpta$vnames$lv.regular[[g]]
    # ))

    # regressions
    reg_idx <- which(pt_1$op == "~" & pt_1$block == block_values[g] &
      pt_1$lhs %in% eqs_names &
      pt_1$rhs %in% eqs_names)

    # the variances
    var_idx <- which(pt_1$op == "~~" & pt_1$block == block_values[g] &
      pt_1$lhs %in% eqs_names &
      pt_1$rhs %in% eqs_names &
      pt_1$lhs == pt_1$rhs)

    # optionally covariances (exo!)
    cov_idx <- which(pt_1$op == "~~" & pt_1$block == block_values[g] &
      pt_1$lhs %in% eqs_names &
      pt_1$rhs %in% eqs_names &
      pt_1$lhs != pt_1$rhs)

    # means/intercepts
    int_idx <- which(pt_1$op == "~1" & pt_1$block == block_values[g] &
      pt_1$lhs %in% eqs_names)

    keep_idx <- c(
      keep_idx, reg_idx, var_idx, cov_idx, int_idx
    )

    # defined/constraints
    if (any(pt_1$op %in% c("==", "<", ">", ":="))) {
      # get the 'id' numbers and the labels involved in def/constraints
      pt2 <- pt_1
      pt2$free <- pt_1$id # use 'id' numbers instead of 'free' indices
      id <- lav_partable_constraints_label_id(pt2, def = TRUE)
      label <- names(id)

      # what are the row indices that we currently keep?
      free_id <- pt_1$id[keep_idx]
    }

    # defined parameters
    def_idx <- which(pt_1$op == ":=")
    if (length(def_idx) > 0L) {
      def_keep <- logical(length(def_idx))
      for (def in seq_along(def_idx)) {
        # rhs
        rhs_labels <- all.vars(as.formula(paste(
          "~",
          pt_1[def_idx[def], "rhs"]
        )))
        if (length(rhs_labels) > 0L) {
          # par id
          rhs_freeid <- id[match(rhs_labels, label)]

          # keep?
          if (all(rhs_freeid %in% free_id)) {
            def_keep[def] <- TRUE
          }
        } else { # only constants?
          def_keep[def] <- TRUE
        }
      }
      keep_idx <- c(keep_idx, def_idx[def_keep])
      # add 'id' numbers of := definitions that we keep
      free_id <- c(free_id, pt_1$id[def_idx[def_keep]])
    }

    # (in)equality constraints
    con_idx <- which(pt_1$op %in% c("==", "<", ">"))
    if (length(con_idx) > 0L) {
      con_keep <- logical(length(con_idx))
      for (con in seq_along(con_idx)) {
        lhs_keep <- FALSE
        rhs_keep <- FALSE

        # lhs
        lhs_labels <- all.vars(as.formula(paste(
          "~",
          pt_1[con_idx[con], "lhs"]
        )))
        if (length(lhs_labels) > 0L) {
          # par id
          lhs_freeid <- id[match(lhs_labels, label)]

          # keep?
          if (all(lhs_freeid %in% free_id)) {
            lhs_keep <- TRUE
          }
        } else {
          lhs_keep <- TRUE
        }

        # rhs
        rhs_labels <- all.vars(as.formula(paste(
          "~",
          pt_1[con_idx[con], "rhs"]
        )))
        if (length(rhs_labels) > 0L) {
          # par id
          rhs_freeid <- id[match(rhs_labels, label)]

          # keep?
          if (all(rhs_freeid %in% free_id)) {
            rhs_keep <- TRUE
          }
        } else {
          rhs_keep <- TRUE
        }

        if (lhs_keep && rhs_keep) {
          con_keep[con] <- TRUE
        }
      }

      keep_idx <- c(keep_idx, con_idx[con_keep])
    } # con
  } # block

  if (idx_only) {
    return(keep_idx)
  }

  pt_1 <- pt_1[keep_idx, , drop = FALSE]

  # add any missing covariances among exogenous variables
  if (add_exo_cov) {
    if (conditional_x) {
      pt_1 <- lav_partable_add_exo_cov(pt_1, ov_names_x = lavpta$vnames$ov.x)
    } else {
      pt_1 <- lav_partable_add_exo_cov(pt_1)
    }
  }

  # if meanstructure, 'free' user=0 intercepts
  if (meanstructure) {
    int_idx <- which(pt_1$op == "~1" & pt_1$user == 0L & pt_1$free == 0L)
    if (length(int_idx) > 0L) {
      pt_1$free[int_idx] <- max(pt_1$free) + seq_along(int_idx)
      pt_1$ustart[int_idx] <- as.numeric(NA)
      pt_1$user[int_idx] <- 3L
    }
  }

  # if fixed.x = FALSE, remove all remaining (free) exo=1 elements
  if (!fixed_x) {
    exo_idx <- which(pt_1$exo != 0L)
    if (length(exo_idx) > 0L) {
      pt_1$exo[exo_idx] <- 0L
      pt_1$user[exo_idx] <- 3L
      pt_1$free[exo_idx] <- max(pt_1$free) + seq_along(exo_idx)
    }

  } else if (conditional_x) {
    # don't change exo column!
  } else {
    # first, wipe out exo column
    exo_idx <- which(pt_1$exo != 0L)
    if (length(exo_idx) > 0L) {
      pt_1$exo[exo_idx] <- 0L
    }

    # keep ov.x as in global model!
    for (g in 1:nblocks) {
      # ov.names.x <- lav_partable_vnames(PT,
      #   type = "ov.x",
      #   block = block.values[g]
      # )
      ov_names_x <- lavpta$vnames$ov.x[[g]]
      if (length(ov_names_x) == 0L) {
        next
      }

      # 1. variances/covariances
      exo_var_idx <- which(
        pt_1$op == "~~" &
          pt_1$block == block_values[g] &
          pt_1$rhs %in% ov_names_x &
          pt_1$lhs %in% ov_names_x &
          pt_1$user %in% c(0L, 3L)
      )
      if (length(exo_var_idx) > 0L) {
        pt_1$ustart[exo_var_idx] <- as.numeric(NA) # to be overridden
        pt_1$free[exo_var_idx] <- 0L
        pt_1$exo[exo_var_idx] <- 1L
        pt_1$user[exo_var_idx] <- 3L
      }

      # 2. intercepts
      exo_int_idx <- which(
        pt_1$op == "~1" &
          pt_1$block == block_values[g] &
          pt_1$lhs %in% ov_names_x &
          pt_1$user == 0L
      )
      if (length(exo_int_idx) > 0L) {
        pt_1$ustart[exo_int_idx] <- as.numeric(NA) # to be overridden
        pt_1$free[exo_int_idx] <- 0L
        pt_1$exo[exo_int_idx] <- 1L
        pt_1$user[exo_int_idx] <- 3L
      }
    } # blocks
  } # fixed.x

  # if conditional.x, check if we have 'additional' ov.x variables
  # that were not ov.x in the global model
  if (conditional_x) {
    for (b in 1:nblocks) {
      global_ov_x <- lavpta$vnames$ov.x[[b]]
      local_ov_x <- lav_partable_vnames(pt_1, type = "ov.x",
        block = block_values[b]
      )
      extra_idx <- which(!local_ov_x %in% global_ov_x)
      if (length(extra_idx) > 0L) {
        extra_ov_names <- local_ov_x[extra_idx]
        for (i in seq_along(extra_idx)) {
          add <- list(
            lhs = extra_ov_names[i],
            op = "~",
            rhs = global_ov_x[1],
            user = 3L,
            free = 0L,
            block = b,
            ustart = 0,
            exo = 1L
          )
          # add group column
          if (!is.null(pt_1$group)) {
            add$group <- unique(pt_1$block[pt_1$block == b])
          }
          # add level column
          if (!is.null(pt_1$level)) {
            add$level <- unique(pt_1$level[pt_1$block == b])
          }
          # add lower column
          if (!is.null(pt_1$lower)) {
            add$lower <- as.numeric(-Inf)
          }
          # add upper column
          if (!is.null(pt_1$upper)) {
            add$upper <- as.numeric(+Inf)
          }
          pt_1 <- lav_partable_add(pt_1, add = add)
        } # i
      } # extra.idx
    } # b
  } # conditional.x

  # if free.fixed.var, free up all 'fixed (to unity)' variances
  if (free_fixed_var) {
    fixed_var_idx <- which(pt_1$op == "~~" & pt_1$lhs == pt_1$rhs &
                       pt_1$free == 0 & pt_1$user == 0L & pt_1$ustart == 1)
    if (length(fixed_var_idx) > 0L) {
      pt_1$free[fixed_var_idx] <- max(pt_1$free) + seq_along(fixed_var_idx)
      pt_1$ustart[fixed_var_idx] <- as.numeric(NA)
    }
  }

  # clean up
  pt_1 <- lav_partable_complete(pt_1)

  if (add_idx) {
    attr(pt_1, "idx") <- keep_idx
  }

  pt_1
}

# NOTE: only within same level
lav_partable_add_exo_cov <- function(pt_1, ov_names_x = NULL) {
  # PT
  pt_1 <- as.data.frame(pt_1, stringsAsFactors = FALSE)

  # lavpta
  lavpta <- lav_partable_attributes(pt_1)

  # nblocks
  nblocks <- lavpta$nblocks
  block_values <- lav_partable_block_values(pt_1)


  # ov.names.x: list with element per block
  if (is.null(ov_names_x)) {
    ov_names_x <- lavpta$vnames$ov.x
  } else if (!is.list(ov_names_x)) {
    ov_names_x <- rep(list(ov_names_x), nblocks)
  }

  # remove ov.names.x if not present at same level/block
  if (nblocks > 1L) {
    for (b in seq_len(nblocks)) {
      rm_idx <- which(!ov_names_x[[b]] %in% lavpta$vnames$ov.x[[b]])
      if (length(rm_idx) > 0L) {
        ov_names_x[[b]] <- ov_names_x[[b]][-rm_idx]
      }
    } # b
  }

  # add covariances among latent variables
  for (b in seq_len(nblocks)) {
    if (length(ov_names_x[[b]]) > 1L) {
      tmp <- utils::combn(ov_names_x[[b]], 2L)
      for (i in seq_len(ncol(tmp))) {
        # already present?
        cov1_idx <- which(pt_1$op == "~~" &
          pt_1$block == block_values[b] &
          pt_1$lhs == tmp[1, i] & pt_1$rhs == tmp[2, i])
        cov2_idx <- which(pt_1$op == "~~" &
          pt_1$block == block_values[b] &
          pt_1$lhs == tmp[2, i] & pt_1$rhs == tmp[1, i])

        # if not, add
        if (length(c(cov1_idx, cov2_idx)) == 0L) {
          add <- list(
            lhs = tmp[1, i],
            op = "~~",
            rhs = tmp[2, i],
            user = 3L,
            free = max(pt_1$free) + 1L,
            block = b,
            ustart = as.numeric(NA)
          )
          # add group column
          if (!is.null(pt_1$group)) {
            add$group <- unique(pt_1$block[pt_1$block == b])
          }
          # add level column
          if (!is.null(pt_1$level)) {
            add$level <- unique(pt_1$level[pt_1$block == b])
          }
          # add lower column
          if (!is.null(pt_1$lower)) {
            add$lower <- as.numeric(-Inf)
          }
          # add upper column
          if (!is.null(pt_1$upper)) {
            add$upper <- as.numeric(+Inf)
          }
          pt_1 <- lav_partable_add(pt_1, add = add)
        }
      }
    } # ov.names.x
  } # blocks

  pt_1
}
