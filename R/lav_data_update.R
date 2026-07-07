# update lavdata object
# - new dataset (lav_data_update)
# - only subset of data (lav_data_update_subset) (for sam())

# YR - 18 Jan 2021 (so we don't need to export lav_data_*_patterns functions)
#    - 28 May 2023 lav_data_update_subset()

# update lavdata object with new dataset
# - assuming everything else stays the same
# - optionally, also provide boot_idx (per group) to adapt internal slots
#   (eXo, sampling weights); the row order of new_x must match boot_idx
# - optionally, also provide boot_clus (per group): a matrix of (relabeled)
#   cluster ids for the new data. This is used by the cluster bootstrap so
#   that a cluster that is resampled more than once is treated as several
#   distinct clusters. If NULL, the original cluster ids are reused.
lav_data_update <- function(lavdata = NULL,
                            new_x = NULL,
                            boot_idx = NULL,
                            boot_clus = NULL,
                            lavoptions = NULL,
                            ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  stopifnot(length(new_x) == lavdata@ngroups)
  stopifnot(!is.null(lavoptions))
  newdata <- lavdata

  # replace data 'X' slot for each group
  for (g in 1:lavdata@ngroups) {
    # replace raw data
    newdata@X[[g]] <- new_x[[g]]

    # if boot_idx is provided, also adapt eXo and sampling weights, using
    # the same (per-group) row ordering as new_x
    if (!is.null(boot_idx)) {
      boot_idx_g <- boot_idx[[g]]

      # eXo
      if (!is.null(lavdata@eXo[[g]])) {
        newdata@eXo[[g]] <- lavdata@eXo[[g]][boot_idx_g, , drop = FALSE]
      }

      # sampling weights
      if (!is.null(lavdata@weights[[g]])) {
        newdata@weights[[g]] <- lavdata@weights[[g]][boot_idx_g]
      }
    }

    # Mp + nobs
    if (lavoptions$missing != "listwise") {
      newdata@Mp[[g]] <- lav_data_mi_patterns(new_x[[g]],
        sort_freq = FALSE, coverage = TRUE
      )
      newdata@nobs[[g]] <-
        (nrow(newdata@X[[g]]) - length(newdata@Mp[[g]]$empty.idx))
    } else {
      # the number of rows may have changed (e.g., cluster bootstrap, where
      # the total N varies from sample to sample); keep nobs in sync
      newdata@nobs[[g]] <- nrow(newdata@X[[g]])
    }

    # Rp
    if (length(lavdata@ov.names.x[[g]]) == 0L &&
      all(lavdata@ov.names[[g]] %in%
        lavdata@ov$name[lavdata@ov$type == "ordered"])) {
      newdata@Rp[[g]] <- lav_data_resp_patterns(new_x[[g]])
    }

    # Lp (cluster structure): rebuild whenever there is a cluster variable.
    # This covers both two-level models (nlevels > 1) and single-level models
    # with cluster-robust standard errors (cluster = , nlevels == 1).
    if (length(lavdata@cluster) > 0L) {
      multilevel <- lavdata@nlevels > 1L
      n_clus_vars <- length(lavdata@cluster)
      # cluster id(s) for the new data
      if (!is.null(boot_clus)) {
        # cluster bootstrap: use the (relabeled) cluster ids, so that a
        # cluster that was resampled more than once is treated as several
        # distinct clusters
        clus <- boot_clus[[g]]
      } else {
        # reuse the original cluster ids (cluster.idx is stored from level 2)
        clus <- matrix(0L, nrow(new_x[[g]]), n_clus_vars)
        for (cc in seq_len(n_clus_vars)) {
          clus[, cc] <- lavdata@Lp[[g]]$cluster.idx[[cc + 1L]]
        }
      }
      # ALWAYS add ov.names.x at the end (as in lav_data_full())
      ov_names_2 <- unique(c(lavdata@ov.names[[g]], lavdata@ov.names.x[[g]]))
      newdata@Lp[[g]] <- lav_data_cl_patterns(
        y = new_x[[g]],
        clus = clus,
        cluster = lavdata@cluster,
        multilevel = multilevel,
        ov_names = ov_names_2,
        ov_names_x = lavdata@ov.names.x[[g]],
        ov_names_l = lavdata@ov.names.l[[g]]
      )
    }
  } # g

  # return updated data object
  newdata
}

# update lavdata, keeping only a subset of the observed variables
# (assuming everything else stays the same)
lav_data_update_subset <- function(lavdata = NULL, ov_names = NULL) {
  stopifnot(length(ov_names) == length(lavdata@ov.names))
  newdata <- lavdata

  # replace ov.names
  newdata@ov.names <- ov_names

  # ordered?
  if (length(lavdata@ordered) > 0L) {
    newdata@ordered <- lavdata@ordered[lavdata@ordered %in% ov_names]
  }

  # replace/update slots for each group
  for (g in 1:lavdata@ngroups) {
    # sanity check:
    if (all(lavdata@ov.names[[g]] %in% ov_names[[g]])) {
      # nothing to do
      next
    }

    # replace ov.names.x
    if (length(lavdata@ov.names.x[[g]]) > 0L) {
      newdata@ov.names.x[[g]] <-
        lavdata@ov.names.x[[g]][lavdata@ov.names.x[[g]] %in%
          ov_names[[g]]]
    }

    # replace ov.names.l
    if (newdata@nlevels > 1L) {
      for (l in 1:newdata@nlevels) {
        newdata@ov.names.l[[g]][[l]] <-
          lavdata@ov.names.l[[g]][[l]][lavdata@ov.names.l[[g]][[l]] %in%
            ov_names[[g]]]
      }
    }

    # ov table
    keep_idx <- which(lavdata@ov$name %in% unlist(ov_names))
    newdata@ov <- lapply(lavdata@ov, "[", keep_idx)

    # replace raw data
    newdata@X[[g]] <- lavdata@X[[g]][, lavdata@ov.names[[g]] %in%
      ov_names[[g]], drop = FALSE]

    # eXo
    if (length(newdata@ov.names.x[[g]]) == 0L) {
      newdata@eXo[g] <- list(NULL)
    } else {
      newdata@eXo[[g]] <- lavdata@eXo[[g]][, lavdata@ov.names.x[[g]] %in%
        ov_names[[g]], drop = FALSE]
    }

    # Mp + nobs
    if (lavdata@missing != "listwise") {
      newdata@Mp[[g]] <- lav_data_mi_patterns(newdata@X[[g]],
        sort_freq = FALSE, coverage = TRUE
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

    # Lp: rebuild whenever there is a cluster variable; this covers both
    # two-level models (nlevels > 1) and single-level models with
    # cluster-robust standard errors (cluster = , nlevels == 1)
    if (length(newdata@cluster) > 0L) {
      # extract cluster variable(s), for this group
      # (cluster.idx is stored from level 2 onwards)
      n_clus_vars <- length(newdata@cluster)
      clus <- matrix(0L, nrow(newdata@X[[g]]), n_clus_vars)
      for (cc in seq_len(n_clus_vars)) {
        clus[, cc] <- lavdata@Lp[[g]]$cluster.idx[[cc + 1L]]
      }
      multilevel <- newdata@nlevels > 1L
      ov_names_1 <- unique(c(ov_names[[g]], newdata@ov.names.x[[g]]))
      newdata@Lp[[g]] <- lav_data_cl_patterns(
        y = newdata@X[[g]],
        clus = clus,
        cluster = newdata@cluster,
        multilevel = multilevel,
        ov_names = ov_names_1,
        ov_names_x = newdata@ov.names.x[[g]],
        ov_names_l = newdata@ov.names.l[[g]]
      )
    }
  } # g

  # return update data object
  newdata
}
