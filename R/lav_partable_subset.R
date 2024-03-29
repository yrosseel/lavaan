# YR 11 feb 2017: initial version

# given a parameter table (PT), extract a part of the model:
# eg.:
# - only the measurement model (with saturated latent variables)
# - only the stuctural part
# - a single measurement block
# ...

# YR 25 June 2021: - add.exo.cov = TRUE for structural model
#                  - fixed.x = FALSE/TRUE -> exo flags


# FIXME:
# - but fixed-to-zero covariances may not be present in PT...
# - if indicators are regressed on exogenous covariates, should we
#   add them here? (no for now, unless add.ind.predictors = TRUE)
lav_partable_subset_measurement_model <- function(PT = NULL,
                                                  lv.names = NULL,
                                                  add.lv.cov = TRUE,
                                                  add.ind.predictors = FALSE,
                                                  add.idx = FALSE,
                                                  idx.only = FALSE) {
  # PT
  PT <- as.data.frame(PT, stringsAsFactors = FALSE)

  # lavpta
  lavpta <- lav_partable_attributes(PT)

  # nblocks
  nblocks <- lavpta$nblocks
  block.values <- lav_partable_block_values(PT)

  # lv.names: list with element per block
  if (is.null(lv.names)) {
    lv.names <- lavpta$vnames$lv.regular
  } else if (!is.list(lv.names)) {
    lv.names <- rep(list(lv.names), nblocks)
  }

  # keep rows idx
  keep.idx <- integer(0L)

  # remove not-needed measurement models
  for (g in 1:nblocks) {
    # indicators for latent variables we keep
    IND.idx <- which(PT$op == "=~" &
      PT$lhs %in% lv.names[[g]] &
      PT$block == block.values[g])
    IND <- PT$rhs[IND.idx]
    IND.plabel <- PT$plabel[IND.idx]

    # keep =~
    keep.idx <- c(keep.idx, IND.idx)

    # new in 0.6-17: indicators regressed on predictors
    if (add.ind.predictors) {
      PRED.idx <- which(PT$op == "~" &
        PT$lhs %in% IND &
        PT$block == block.values[g])
      EXTRA <- unique(PT$rhs[PRED.idx])
      keep.idx <- c(keep.idx, PRED.idx)
      # add them to IND, so we include their variances/intercepts
      IND <- c(IND, EXTRA)
    }


    # keep ~~
    OV.VAR.idx <- which(PT$op == "~~" &
      PT$lhs %in% IND &
      PT$rhs %in% IND &
      PT$block == block.values[g])
    keep.idx <- c(keep.idx, OV.VAR.idx)

    LV.VAR.idx <- which(PT$op == "~~" &
      PT$lhs %in% lv.names[[g]] &
      PT$rhs %in% lv.names[[g]] &
      PT$block == block.values[g])
    keep.idx <- c(keep.idx, LV.VAR.idx)

    # intercepts indicators
    OV.INT.idx <- which(PT$op == "~1" &
      PT$lhs %in% IND &
      PT$block == block.values[g])
    keep.idx <- c(keep.idx, OV.INT.idx)

    # intercepts latent variables
    LV.INT.idx <- which(PT$op == "~1" &
      PT$lhs %in% lv.names[[g]] &
      PT$block == block.values[g])
    keep.idx <- c(keep.idx, LV.INT.idx)

    # thresholds
    TH.idx <- which(PT$op == "|" &
      PT$lhs %in% IND &
      PT$block == block.values[g])
    keep.idx <- c(keep.idx, TH.idx)

    # scaling factors
    SC.idx <- which(PT$op == "~*~" &
      PT$lhs %in% IND &
      PT$block == block.values[g])
    keep.idx <- c(keep.idx, SC.idx)

    # defined/constraints
    if (any(PT$op %in% c("==", "<", ">", ":="))) {
      # get the 'id' numbers and the labels involved in def/constraints
      PT2 <- PT
      PT2$free <- PT$id # us 'id' numbers instead of 'free' indices
      ID <- lav_partable_constraints_label_id(PT2, def = TRUE)
      LABEL <- names(ID)

      # what are the row indices that we currently keep?
      FREE.id <- PT$id[keep.idx]
    }

    # defined parameters
    def.idx <- which(PT$op == ":=")
    if (length(def.idx) > 0L) {
      def.keep <- logical(length(def.idx))
      for (def in seq_len(length(def.idx))) {
        # rhs
        RHS.labels <- all.vars(as.formula(paste(
          "~",
          PT[def.idx[def], "rhs"]
        )))
        if (length(RHS.labels) > 0L) {
          # par id
          RHS.freeid <- ID[match(RHS.labels, LABEL)]

          # keep?
          if (all(RHS.freeid %in% FREE.id)) {
            def.keep[def] <- TRUE
          }
        } else { # only constants?
          def.keep[def] <- TRUE
        }
      }
      keep.idx <- c(keep.idx, def.idx[def.keep])
      # add 'id' numbers of := definitions that we keep
      FREE.id <- c(FREE.id, PT$id[def.idx[def.keep]])
    }

    # (in)equality constraints
    con.idx <- which(PT$op %in% c("==", "<", ">"))
    if (length(con.idx) > 0L) {
      con.keep <- logical(length(con.idx))
      for (con in seq_len(length(con.idx))) {
        lhs.keep <- FALSE
        rhs.keep <- FALSE

        # lhs
        LHS.labels <- all.vars(as.formula(paste(
          "~",
          PT[con.idx[con], "lhs"]
        )))
        if (length(LHS.labels) > 0L) {
          # par id
          LHS.freeid <- ID[match(LHS.labels, LABEL)]

          # keep?
          if (all(LHS.freeid %in% FREE.id)) {
            lhs.keep <- TRUE
          }
        } else {
          lhs.keep <- TRUE
        }

        # rhs
        RHS.labels <- all.vars(as.formula(paste(
          "~",
          PT[con.idx[con], "rhs"]
        )))
        if (length(RHS.labels) > 0L) {
          # par id
          RHS.freeid <- ID[match(RHS.labels, LABEL)]

          # keep?
          if (all(RHS.freeid %in% FREE.id)) {
            rhs.keep <- TRUE
          }
        } else {
          rhs.keep <- TRUE
        }

        if (lhs.keep && rhs.keep) {
          con.keep[con] <- TRUE
        }
      }

      keep.idx <- c(keep.idx, con.idx[con.keep])
    } # con
  } # block

  if (idx.only) {
    return(keep.idx)
  }

  PT <- PT[keep.idx, , drop = FALSE]

  # check if we have enough indicators?
  # TODO

  # add covariances among latent variables?
  if (add.lv.cov) {
    PT <- lav_partable_add_lv_cov(
      PT = PT,
      lv.names = lv.names
    )
  }

  # clean up
  PT <- lav_partable_complete(PT)

  if (add.idx) {
    attr(PT, "idx") <- keep.idx
  }

  PT
}

# NOTE: only within same level
lav_partable_add_lv_cov <- function(PT, lv.names = NULL) {
  # PT
  PT <- as.data.frame(PT, stringsAsFactors = FALSE)

  # lavpta
  lavpta <- lav_partable_attributes(PT)

  # nblocks
  nblocks <- lavpta$nblocks
  block.values <- lav_partable_block_values(PT)


  # lv.names: list with element per block
  if (is.null(lv.names)) {
    lv.names <- lavpta$vnames$lv.regular
  } else if (!is.list(lv.names)) {
    lv.names <- rep(list(lv.names), nblocks)
  }

  # remove lv.names if not present at same level/block
  if (nblocks > 1L) {
    for (b in seq_len(nblocks)) {
      rm.idx <- which(!lv.names[[b]] %in% lavpta$vnames$lv.regular[[b]])
      if (length(rm.idx) > 0L) {
        lv.names[[b]] <- lv.names[[b]][-rm.idx]
      }
    } # b
  }

  # add covariances among latent variables
  for (b in seq_len(nblocks)) {
    if (length(lv.names[[b]]) > 1L) {
      tmp <- utils::combn(lv.names[[b]], 2L)
      for (i in seq_len(ncol(tmp))) {
        # already present?
        cov1.idx <- which(PT$op == "~~" &
          PT$block == block.values[b] &
          PT$lhs == tmp[1, i] & PT$rhs == tmp[2, i])
        cov2.idx <- which(PT$op == "~~" &
          PT$block == block.values[b] &
          PT$lhs == tmp[2, i] & PT$rhs == tmp[1, i])

        # if not, add
        if (length(c(cov1.idx, cov2.idx)) == 0L) {
          ADD <- list(
            lhs = tmp[1, i],
            op = "~~",
            rhs = tmp[2, i],
            user = 3L,
            free = max(PT$free) + 1L,
            block = b
          )
          # add group column
          if (!is.null(PT$group)) {
            ADD$group <- unique(PT$block[PT$block == b])
          }
          # add level column
          if (!is.null(PT$level)) {
            ADD$level <- unique(PT$level[PT$block == b])
          }
          # add lower column
          if (!is.null(PT$lower)) {
            ADD$lower <- as.numeric(-Inf)
          }
          # add upper column
          if (!is.null(PT$upper)) {
            ADD$upper <- as.numeric(+Inf)
          }
          PT <- lav_partable_add(PT, add = ADD)
        }
      }
    } # lv.names
  } # blocks

  PT
}


# this function takes a 'full' SEM (measurement models + structural part)
# and returns only the structural part
#
# - what to do if we have no regressions among the latent variables?
#   we return all covariances among the latent variables
#
# - also, we should check if we have any 'higher' order factors
#
lav_partable_subset_structural_model <- function(PT = NULL,
                                                 add.idx = FALSE,
                                                 idx.only = FALSE,
                                                 add.exo.cov = FALSE,
                                                 fixed.x = FALSE,
                                                 meanstructure = FALSE) {
  # PT
  PT <- as.data.frame(PT, stringsAsFactors = FALSE)

  # remove any EFA related information -- new in 0.6-18
  if (!is.null(PT$efa)) {
    PT$efa <- NULL
    PT$est.unrotated <- NULL
    seven.idx <- which(PT$user == 7L & PT$op == "~~")
    if (length(seven.idx) > 0L) {
      PT$user[seven.idx] <- 0L
      PT$free[seven.idx] <- 1L
      PT$ustart[seven.idx] <- as.numeric(NA)
      PT$est[seven.idx] <- PT$est.std[seven.idx]
    }
    PT$est.std <- NULL
  }

  # lavpta
  lavpta <- lav_partable_attributes(PT)

  # nblocks
  nblocks <- lavpta$nblocks
  block.values <- lav_partable_block_values(PT)

  # eqs.names
  eqs.x.names <- lavpta$vnames$eqs.x
  eqs.y.names <- lavpta$vnames$eqs.y
  lv.names <- lavpta$vnames$lv.regular

  # keep rows idx
  keep.idx <- integer(0L)

  # remove not-needed measurement models
  for (g in 1:nblocks) {
    # higher-order factor loadings
    fac.idx <- which(PT$op == "=~" & PT$block == block.values[g] &
      PT$lhs %in% lavpta$vnames$lv.regular[[g]] &
      PT$rhs %in% lavpta$vnames$lv.regular[[g]])

    # eqs.names
    eqs.names <- unique(c(
      lavpta$vnames$eqs.x[[g]],
      lavpta$vnames$eqs.y[[g]]
    ))
    all.names <- unique(c(
      eqs.names,
      lavpta$vnames$lv.regular[[g]]
    ))

    # regressions
    reg.idx <- which(PT$op == "~" & PT$block == block.values[g] &
      PT$lhs %in% eqs.names &
      PT$rhs %in% eqs.names)

    # the variances
    var.idx <- which(PT$op == "~~" & PT$block == block.values[g] &
      PT$lhs %in% all.names &
      PT$rhs %in% all.names &
      PT$lhs == PT$rhs)

    # optionally covariances (exo!)
    cov.idx <- which(PT$op == "~~" & PT$block == block.values[g] &
      PT$lhs %in% all.names &
      PT$rhs %in% all.names &
      PT$lhs != PT$rhs)

    # means/intercepts
    int.idx <- which(PT$op == "~1" & PT$block == block.values[g] &
      PT$lhs %in% all.names)

    keep.idx <- c(
      keep.idx, reg.idx, var.idx, cov.idx, int.idx,
      fac.idx
    )

    # defined/constraints
    if (any(PT$op %in% c("==", "<", ">", ":="))) {
      # get the 'id' numbers and the labels involved in def/constraints
      PT2 <- PT
      PT2$free <- PT$id # us 'id' numbers instead of 'free' indices
      ID <- lav_partable_constraints_label_id(PT2, def = TRUE)
      LABEL <- names(ID)

      # what are the row indices that we currently keep?
      FREE.id <- PT$id[keep.idx]
    }

    # defined parameters
    def.idx <- which(PT$op == ":=")
    if (length(def.idx) > 0L) {
      def.keep <- logical(length(def.idx))
      for (def in seq_len(length(def.idx))) {
        # rhs
        RHS.labels <- all.vars(as.formula(paste(
          "~",
          PT[def.idx[def], "rhs"]
        )))
        if (length(RHS.labels) > 0L) {
          # par id
          RHS.freeid <- ID[match(RHS.labels, LABEL)]

          # keep?
          if (all(RHS.freeid %in% FREE.id)) {
            def.keep[def] <- TRUE
          }
        } else { # only constants?
          def.keep[def] <- TRUE
        }
      }
      keep.idx <- c(keep.idx, def.idx[def.keep])
      # add 'id' numbers of := definitions that we keep
      FREE.id <- c(FREE.id, PT$id[def.idx[def.keep]])
    }

    # (in)equality constraints
    con.idx <- which(PT$op %in% c("==", "<", ">"))
    if (length(con.idx) > 0L) {
      con.keep <- logical(length(con.idx))
      for (con in seq_len(length(con.idx))) {
        lhs.keep <- FALSE
        rhs.keep <- FALSE

        # lhs
        LHS.labels <- all.vars(as.formula(paste(
          "~",
          PT[con.idx[con], "lhs"]
        )))
        if (length(LHS.labels) > 0L) {
          # par id
          LHS.freeid <- ID[match(LHS.labels, LABEL)]

          # keep?
          if (all(LHS.freeid %in% FREE.id)) {
            lhs.keep <- TRUE
          }
        } else {
          lhs.keep <- TRUE
        }

        # rhs
        RHS.labels <- all.vars(as.formula(paste(
          "~",
          PT[con.idx[con], "rhs"]
        )))
        if (length(RHS.labels) > 0L) {
          # par id
          RHS.freeid <- ID[match(RHS.labels, LABEL)]

          # keep?
          if (all(RHS.freeid %in% FREE.id)) {
            rhs.keep <- TRUE
          }
        } else {
          rhs.keep <- TRUE
        }

        if (lhs.keep && rhs.keep) {
          con.keep[con] <- TRUE
        }
      }

      keep.idx <- c(keep.idx, con.idx[con.keep])
    } # con
  } # block

  if (idx.only) {
    return(keep.idx)
  }

  PT <- PT[keep.idx, , drop = FALSE]

  # add any missing covariances among exogenous variables
  if (add.exo.cov) {
    PT <- lav_partable_add_exo_cov(PT)
  }

  # if meanstructure, 'free' user=0 intercepts
  if (meanstructure) {
    int.idx <- which(PT$op == "~1" & PT$user == 0L & PT$free == 0L)
    if (length(int.idx) > 0L) {
      PT$free[int.idx] <- max(PT$free) + seq_len(length(int.idx))
      PT$ustart[int.idx] <- as.numeric(NA)
      PT$user[int.idx] <- 3L
    }
  }

  # if fixed.x = FALSE, remove all remaining (free) exo=1 elements
  if (!fixed.x) {
    exo.idx <- which(PT$exo != 0L)
    if (length(exo.idx) > 0L) {
      PT$exo[exo.idx] <- 0L
      PT$user[exo.idx] <- 3L
      PT$free[exo.idx] <- max(PT$free) + seq_len(length(exo.idx))
    }

    # if fixed.x = TRUE, check/set all exo elements
  } else {
    # redefine ov.x for the structural part only; set exo flag
    for (g in 1:nblocks) {
      ov.names.x <- lav_partable_vnames(PT,
        type = "ov.x",
        block = block.values[g]
      )
      if (length(ov.names.x) == 0L) {
        next
      }

      # 1. variances/covariances
      exo.var.idx <- which(
        PT$op == "~~" &
          PT$block == block.values[g] &
          PT$rhs %in% ov.names.x &
          PT$lhs %in% ov.names.x &
          PT$user %in% c(0L, 3L)
      )
      if (length(exo.var.idx) > 0L) {
        PT$ustart[exo.var.idx] <- as.numeric(NA) # to be overriden
        PT$free[exo.var.idx] <- 0L
        PT$exo[exo.var.idx] <- 1L
        PT$user[exo.var.idx] <- 3L
      }

      # 2. intercepts
      exo.int.idx <- which(
        PT$op == "~1" &
          PT$block == block.values[g] &
          PT$lhs %in% ov.names.x &
          PT$user == 0L
      )
      if (length(exo.int.idx) > 0L) {
        PT$ustart[exo.int.idx] <- as.numeric(NA) # to be overriden
        PT$free[exo.int.idx] <- 0L
        PT$exo[exo.int.idx] <- 1L
        PT$user[exo.var.idx] <- 3L
      }
    } # blocks
  } # fixed.x

  # clean up
  PT <- lav_partable_complete(PT)

  if (add.idx) {
    attr(PT, "idx") <- keep.idx
  }

  PT
}

# NOTE: only within same level
lav_partable_add_exo_cov <- function(PT, ov.names.x = NULL) {
  # PT
  PT <- as.data.frame(PT, stringsAsFactors = FALSE)

  # lavpta
  lavpta <- lav_partable_attributes(PT)

  # nblocks
  nblocks <- lavpta$nblocks
  block.values <- lav_partable_block_values(PT)


  # ov.names.x: list with element per block
  if (is.null(ov.names.x)) {
    ov.names.x <- lavpta$vnames$ov.x
  } else if (!is.list(ov.names.x)) {
    ov.names.x <- rep(list(ov.names.x), nblocks)
  }

  # remove ov.names.x if not present at same level/block
  if (nblocks > 1L) {
    for (b in seq_len(nblocks)) {
      rm.idx <- which(!ov.names.x[[b]] %in% lavpta$vnames$ov.x[[b]])
      if (length(rm.idx) > 0L) {
        ov.names.x[[b]] <- ov.names.x[[b]][-rm.idx]
      }
    } # b
  }

  # add covariances among latent variables
  for (b in seq_len(nblocks)) {
    if (length(ov.names.x[[b]]) > 1L) {
      tmp <- utils::combn(ov.names.x[[b]], 2L)
      for (i in seq_len(ncol(tmp))) {
        # already present?
        cov1.idx <- which(PT$op == "~~" &
          PT$block == block.values[b] &
          PT$lhs == tmp[1, i] & PT$rhs == tmp[2, i])
        cov2.idx <- which(PT$op == "~~" &
          PT$block == block.values[b] &
          PT$lhs == tmp[2, i] & PT$rhs == tmp[1, i])

        # if not, add
        if (length(c(cov1.idx, cov2.idx)) == 0L) {
          ADD <- list(
            lhs = tmp[1, i],
            op = "~~",
            rhs = tmp[2, i],
            user = 3L,
            free = max(PT$free) + 1L,
            block = b,
            ustart = as.numeric(NA)
          )
          # add group column
          if (!is.null(PT$group)) {
            ADD$group <- unique(PT$block[PT$block == b])
          }
          # add level column
          if (!is.null(PT$level)) {
            ADD$level <- unique(PT$level[PT$block == b])
          }
          # add lower column
          if (!is.null(PT$lower)) {
            ADD$lower <- as.numeric(-Inf)
          }
          # add upper column
          if (!is.null(PT$upper)) {
            ADD$upper <- as.numeric(+Inf)
          }
          PT <- lav_partable_add(PT, add = ADD)
        }
      }
    } # ov.names.x
  } # blocks

  PT
}
