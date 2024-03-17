# update lavdata object
# - new dataset (lav_data_update)
# - only subset of data (lav_data_update_subset) (for sam())

# YR - 18 Jan 2021 (so we don't need to export lav_data_*_patterns functions)
#    - 28 May 2023 lav_data_update_subset()

# update lavdata object with new dataset
# - assuming everything else stays the same
# - optionally, also provide boot.idx (per group) to adapt internal slots
lav_data_update <- function(lavdata = NULL, newX = NULL, BOOT.idx = NULL,
                            lavoptions = NULL) {
  stopifnot(length(newX) == lavdata@ngroups)
  stopifnot(!is.null(lavoptions))
  newdata <- lavdata

  # replace data 'X' slot for each group
  for (g in 1:lavdata@ngroups) {
    # replace raw data
    newdata@X[[g]] <- newX[[g]]

    # Mp + nobs
    if (lavoptions$missing != "listwise") {
      newdata@Mp[[g]] <- lav_data_missing_patterns(newX[[g]],
        sort.freq = FALSE, coverage = FALSE
      )
      newdata@nobs[[g]] <-
        (nrow(newdata@X[[g]]) - length(newdata@Mp[[g]]$empty.idx))
    }

    # Rp
    if (length(lavdata@ov.names.x[[g]]) == 0L &&
      all(lavdata@ov.names[[g]] %in%
        lavdata@ov$name[lavdata@ov$type == "ordered"])) {
      newdata@Rp[[g]] <- lav_data_resp_patterns(newX[[g]])
    }

    # Lp
    if (lavdata@nlevels > 1L) {
      # CHECKME!
      # extract cluster variable(s), for this group
      clus <- matrix(0, nrow(newX[[g]]), lavdata@nlevels - 1L)
      for (l in 2:lavdata@nlevels) {
        clus[, (l - 1L)] <- lavdata@Lp[[g]]$cluster.idx[[l]]
      }
      newdata@Lp[[g]] <- lav_data_cluster_patterns(
        Y = newX[[g]],
        clus = clus,
        cluster = lavdata@cluster,
        ov.names = lavdata@ov.names[[g]],
        ov.names.l = lavdata@ov.names.l[[g]]
      )
    }
  }

  # if boot.idx if provided, also adapt eXo and WT
  if (!is.null(BOOT.idx)) {
    boot.idx <- BOOT.idx[[g]]

    # eXo
    if (!is.null(lavdata@eXo[[g]])) {
      newdata@eXo[[g]] <- lavdata@eXo[[g]][boot.idx, , drop = FALSE]
    }

    # sampling weights
    if (!is.null(lavdata@weights[[g]])) {
      newdata@weights[[g]] <- lavdata@weights[[g]][boot.idx]
    }
  } # g

  # return update data object
  newdata
}

# update lavdata, keeping only a subset of the observed variables
# (assuming everything else stays the same)
lav_data_update_subset <- function(lavdata = NULL, ov.names = NULL) {
  stopifnot(length(ov.names) == length(lavdata@ov.names))
  newdata <- lavdata

  # replace ov.names
  newdata@ov.names <- ov.names

  # ordered?
  if (length(lavdata@ordered) > 0L) {
    newdata@ordered <- lavdata@ordered[lavdata@ordered %in% ov.names]
  }

  # replace/update slots for each group
  for (g in 1:lavdata@ngroups) {
    # sanity check:
    if (all(lavdata@ov.names[[g]] %in% ov.names[[g]])) {
      # nothing to do
      next
    }

    # replace ov.names.x
    if (length(lavdata@ov.names.x[[g]]) > 0L) {
      newdata@ov.names.x[[g]] <- lavdata@ov.names.x[[g]][lavdata@ov.names.x[[g]] %in% ov.names[[g]]]
    }

    # replace ov.names.l
    if (newdata@nlevels > 1L) {
      for (l in 1:newdata@nlevels) {
        newdata@ov.names.l[[g]][[l]] <- lavdata@ov.names.l[[g]][[l]][lavdata@ov.names.l[[g]][[l]] %in% ov.names[[g]]]
      }
    }

    # ov table
    keep.idx <- which(lavdata@ov$name %in% unlist(ov.names))
    newdata@ov <- lapply(lavdata@ov, "[", keep.idx)

    # replace raw data
    newdata@X[[g]] <- lavdata@X[[g]][, lavdata@ov.names[[g]] %in% ov.names[[g]], drop = FALSE]

    # eXo
    if (length(newdata@ov.names.x[[g]]) == 0L) {
      newdata@eXo[g] <- list(NULL)
    } else {
      newdata@eXo[[g]] <- lavdata@eXo[[g]][, lavdata@ov.names.x[[g]] %in% ov.names[[g]], drop = FALSE]
    }

    # Mp + nobs
    if (lavdata@missing != "listwise") {
      newdata@Mp[[g]] <- lav_data_missing_patterns(newdata@X[[g]],
        sort.freq = FALSE, coverage = FALSE
      )
      newdata@nobs[[g]] <-
        (nrow(newdata@X[[g]]) - length(newdata@Mp[[g]]$empty.idx))
    }

    # Rp
    if (length(newdata@ordered) == 0L) {
      # nothing to do
    } else if (length(newdata@ov.names.x[[g]]) == 0L &&
      all(newdata@ov.names[[g]] %in%
        newdata@ov$name[newdata@ov$type == "ordered"])) {
      newdata@Rp[[g]] <- lav_data_resp_patterns(newdata@X[[g]])
    }

    # Lp
    if (length(newdata@cluster) > 0L) {
      # extract cluster variable(s), for this group
      clus <- matrix(0, nrow(newdata@X[[g]]), lavdata@nlevels - 1L)
      for (l in 2:lavdata@nlevels) {
        clus[, (l - 1L)] <- lavdata@Lp[[g]]$cluster.idx[[l]]
      }
      if (newdata@nlevels > 1L) {
        multilevel <- TRUE
      } else {
        multilevel <- FALSE
      }
      OV.NAMES <- unique(c(ov.names[[g]], newdata@ov.names.x[[g]]))
      newdata@Lp[[g]] <- lav_data_cluster_patterns(
        Y = newdata@X[[g]],
        clus = clus,
        cluster = newdata@cluster,
        multilevel = multilevel,
        ov.names = OV.NAMES,
        ov.names.x = newdata@ov.names.x[[g]],
        ov.names.l = newdata@ov.names.l[[g]]
      )
    }
  } # g

  # return update data object
  newdata
}
