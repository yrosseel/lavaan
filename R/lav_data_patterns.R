# get missing patterns
lav_data_mi_patterns <- function(y, sort_freq = FALSE, coverage = FALSE,
                                      lp = NULL) {
  # handle two-level data
  if (!is.null(lp)) {
    y_orig <- y
    z <- NULL
    if (length(lp$between.idx[[2]]) > 0L) {
      y <- y[, -lp$between.idx[[2]], drop = FALSE]
      z_idx <- which(!duplicated(lp$cluster.idx[[2]]))
      z <- y_orig[z_idx, lp$between.idx[[2]], drop = FALSE]
    }
  }

  # construct TRUE/FALSE matrix: TRUE if value is observed
  obs <- !is.na(y)

  # empty cases
  empty_idx <- which(rowSums(obs) == 0L)

  # pattern of observed values per observation
  case_id <- apply(1L * obs, 1L, paste, collapse = "")

  # remove empty patterns
  if (length(empty_idx)) {
    case_id_nonempty <- case_id[-empty_idx]
  } else {
    case_id_nonempty <- case_id
  }

  # sort non-empty patterns (from high occurrence to low occurrence)
  if (sort_freq) {
    table_1 <- sort(table(case_id_nonempty), decreasing = TRUE)
  } else {
    table_1 <- table(case_id_nonempty)
  }

  # unique pattern ids
  pat_id <- names(table_1)

  # number of patterns
  pat_npatterns <- length(pat_id)

  # case idx per pattern
  pat_case_idx <- lapply(
    seq_len(pat_npatterns),
    function(p) which(case_id == pat_id[p])
  )

  # unique pattern frequencies
  pat_freq <- as.integer(table_1)

  # first occurrence of each pattern
  pat_first <- match(pat_id, case_id)

  # TRUE/FALSE for each pattern
  pat_obs <- obs[pat_first, , drop = FALSE] # observed per pattern

  mp <- list(
    npatterns = pat_npatterns, id = pat_id, freq = pat_freq,
    case.idx = pat_case_idx, pat = pat_obs, empty.idx = empty_idx,
    nel = sum(obs)
  )

  if (coverage) {
    # FIXME: if we have empty cases, include them in N?
    # no for now
    mp$coverage <- crossprod(obs) / sum(pat_freq)
    # Mp$coverage <- crossprod(OBS) / NROW(Y)
  }

  # additional info if we have two-level data
  if (!is.null(lp)) {
    mp$j.idx <- lapply(
      seq_len(pat_npatterns),
      function(p) lp$cluster.idx[[2]][mp$case.idx[[p]]]
    )
    mp$j1.idx <- lapply(
      seq_len(pat_npatterns),
      function(p) sort(unique.default(mp$j.idx[[p]])) # new in 0.6-19: sort!!
                                                    # because table/tabulate
                            # also sorts...
    )
    mp$j.freq <- lapply(
      seq_len(pat_npatterns),
      function(p) as.integer(unname(table(sort(mp$j.idx[[p]])))) # sort needed?
    )

    # between-level patterns
    if (!is.null(z)) {
      mp$Zp <- lav_data_mi_patterns(z,
        sort_freq = FALSE,
        coverage = FALSE, lp = NULL
      )
    }
  }

  mp
}

# get response patterns (ignore empty cases!)
lav_data_resp_patterns <- function(y) {
  # construct TRUE/FALSE matrix: TRUE if value is observed
  obs <- !is.na(y)

  # empty cases
  empty_idx <- which(rowSums(obs) == 0L)

  # remove empty cases
  if (length(empty_idx) > 0L) {
    y <- y[-empty_idx, , drop = FALSE]
  }

  # ntotal <- nrow(y)
  # nvar <- ncol(y)

  # identify, label and sort response patterns
  id <- apply(y, MARGIN = 1, paste, collapse = "")

  # sort patterns (from high occurrence to low occurrence)
  table_1 <- sort(table(id), decreasing = TRUE)
  order <- names(table_1)
  npatterns <- length(table_1)
  pat <- y[match(order, id), , drop = FALSE]
  row.names(pat) <- as.character(table_1)

  # handle NA?
  y[is.na(y)] <- -9
  total_patterns <- prod(apply(y, 2, function(x) length(unique(x))))
  empty_patterns <- total_patterns - npatterns
  # return a list
  # out <- list(nobs=ntotal, nvar=nvar,
  #            id=id, npatterns=npatterns,
  #            order=order, pat=pat)

  # only return pat
  out <- list(
    npatterns = npatterns, pat = pat, total.patterns = total_patterns,
    empty.patterns = empty_patterns
  )

  out
}

# get cluster information
# - cluster can be a vector!
# - clus can contain multiple columns!
lav_data_cl_patterns <- function(y = NULL,
                                      clus = NULL, # the cluster ids
                                      cluster = NULL, # the cluster 'names'
                                      multilevel = FALSE,
                                      ov_names = NULL,
                                      ov_names_x = NULL,
                                      ov_names_l = NULL) {
  # how many levels?
  nlevels <- length(cluster) + 1L

  # did we get any data (or is this just for lav_data_simulate_old)
  if (!is.null(y) && !is.null(clus)) {
    have_data <- TRUE
  } else {
    have_data <- FALSE
  }

  # check clus
  if (have_data) {
    stopifnot(ncol(clus) == (nlevels - 1L), nrow(y) == nrow(clus))
  }

  cluster_size    <- vector("list", length = nlevels)
  cluster_id      <- vector("list", length = nlevels)
  cluster_idx     <- vector("list", length = nlevels)
  nclusters       <- vector("list", length = nlevels)
  cluster_sizes   <- vector("list", length = nlevels)
  ncluster_sizes  <- vector("list", length = nlevels)
  cluster_size_ns <- vector("list", length = nlevels)

  ov_idx          <- vector("list", length = nlevels)
  ov_x_idx        <- vector("list", length = nlevels)
  ov_y_idx        <- vector("list", length = nlevels)
  both_idx        <- vector("list", length = nlevels)
  within_idx      <- vector("list", length = nlevels)
  within_x_idx    <- vector("list", length = nlevels)
  within_y_idx    <- vector("list", length = nlevels)
  between_idx     <- vector("list", length = nlevels)
  between_x_idx   <- vector("list", length = nlevels)
  between_y_idx   <- vector("list", length = nlevels)
  both_names      <- vector("list", length = nlevels)
  within_names    <- vector("list", length = nlevels)
  within_x_names  <- vector("list", length = nlevels)
  within_y_names  <- vector("list", length = nlevels)
  between_names   <- vector("list", length = nlevels)
  between_x_names <- vector("list", length = nlevels)
  between_y_names <- vector("list", length = nlevels)

  # level-1 is special
  if (have_data) {
    nclusters[[1]] <- NROW(y)
  }

  # higher levels:
  for (l in 2:nlevels) {
    if (have_data) {
      clus_1 <- clus[, (l - 1L)]

    # cluster.id: original cluster identifier
      cluster_id[[l]] <- unique(clus_1)

    # cluster.idx: internal cluster idx (always an integer from 1 to J
    # for every observation;
    # the first cluster id we observe is '1', the second '2', etc...)
      cluster_idx[[l]] <- match(clus_1, cluster_id[[l]])

      cluster_size[[l]] <- tabulate(cluster_idx[[l]])
      nclusters[[l]] <- length(cluster_size[[l]])
      # check if we have more observations than clusters
      if (nclusters[[1]] == nclusters[[l]]) {
        lav_msg_stop(gettext("every cluster contains only one observation."))
      }
      mean_cluster_size <- mean(cluster_size[[l]])
      if (mean_cluster_size < 1.5) {
        lav_msg_warn(gettextf(
          "mean cluster size is %s. This means that many clusters only
          contain a single observation.", mean_cluster_size))
        }
      cluster_sizes[[l]] <- unique(cluster_size[[l]])
      ncluster_sizes[[l]] <- length(cluster_sizes[[l]])
      cluster_size_ns[[l]] <- as.integer(table(factor(cluster_size[[l]],
        levels = as.character(cluster_sizes[[l]])
      )))
    } else {
      cluster_id[[l]] <- integer(0L)
      cluster_idx[[l]] <- integer(0L)
      cluster_size[[l]] <- integer(0L)
      nclusters[[l]] <- integer(0L)
      cluster_sizes[[l]] <- integer(0L)
      ncluster_sizes[[l]] <- integer(0L)
      cluster_size_ns[[l]] <- integer(0L)
    }
  }

  # for all levels:
  if (multilevel) {
    for (l in 1:nlevels) {
      # index of ov.names for this level
      ov_idx[[l]] <- match(ov_names_l[[l]], ov_names)

      # new in 0.6-12: always preserve the order of ov.idx[[l]]
      idx <- which(ov_names %in% ov_names_l[[1]] &
        ov_names %in% ov_names_l[[2]])
      both_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]

      idx <- which(ov_names %in% ov_names_l[[1]] &
        !ov_names %in% ov_names_l[[2]])
      within_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]
      # backwards compatibility: also store in within.idx[[2]]
      if (l == 2) {
        within_idx[[l]] <- within_idx[[1]]
      }

      idx <- which(!ov_names %in% ov_names_l[[1]] &
        ov_names %in% ov_names_l[[2]])
      between_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]

      # names
      # both.names[[l]]     <- ov.names[ ov.names %in% ov.names.l[[1]] &
      #                                 ov.names %in% ov.names.l[[2]] ]
      # within.names[[l]]   <- ov.names[ ov.names %in% ov.names.l[[1]] &
      #                                !ov.names %in% ov.names.l[[2]] ]
      # between.names[[l]]  <- ov.names[!ov.names %in% ov.names.l[[1]] &
      #                                 ov.names %in% ov.names.l[[2]] ]
      both_names[[l]] <- ov_names[both_idx[[l]]]
      within_names[[l]] <- ov_names[within_idx[[l]]]
      between_names[[l]] <- ov_names[between_idx[[l]]]
    }
  }

  # fixed.x wrt variable index
  if (multilevel && length(ov_names_x) > 0L) {
    for (l in 1:nlevels) {
      # some ov.names.x could be 'split', and end up in both.names
      # they should NOT be part ov.x.idx (as they become latent variables)
      idx <- which(ov_names %in% ov_names_x &
        ov_names %in% ov_names_l[[l]] &
        !ov_names %in% unlist(both_names))
      ov_x_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]

      # not any longer, we split them, but still treat them as 'fixed'
      # ov.x.idx[[l]] <- which( ov.names %in% ov.names.x &
      #                        ov.names %in% ov.names.l[[l]] )

      # if some ov.names.x have been 'split', and end up in both.names,
      # they should become part of ov.y.idx (despite being exogenous)
      # as they are now latent variables
      idx <- which(ov_names %in% ov_names_l[[l]] &
        !ov_names %in% ov_names_x[!ov_names_x %in% unlist(both_names)])
      ov_y_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]

      # not any longer, ov.x stays ov.x (even if we split)
      # ov.y.idx[[l]] <- which( ov.names %in% ov.names.l[[l]] &
      #                       !ov.names %in% ov.names.x )


      # if(l == 1L) {
      #    next
      # }
      # below, we only fill in the [[2]] element (and higher)

      idx <- which(ov_names %in% ov_names_l[[1]] &
        !ov_names %in% ov_names_l[[2]] &
        ov_names %in% ov_names_x)
      within_x_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]
      # backwards compatibility: also store in within.x.idx[[2]]
      if (l == 2) {
        within_x_idx[[l]] <- within_x_idx[[1]]
      }

      idx <- which(ov_names %in% ov_names_l[[1]] &
        !ov_names %in% ov_names_l[[2]] &
        !ov_names %in% ov_names_x)
      within_y_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]
      # backwards compatibility: also store in within.y.idx[[2]]
      if (l == 2) {
        within_y_idx[[l]] <- within_y_idx[[1]]
      }

      idx <- which(!ov_names %in% ov_names_l[[1]] &
        ov_names %in% ov_names_l[[2]] &
        ov_names %in% ov_names_x)
      between_x_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]

      idx <- which(!ov_names %in% ov_names_l[[1]] &
        ov_names %in% ov_names_l[[2]] &
        !ov_names %in% ov_names_x)
      between_y_idx[[l]] <- ov_idx[[l]][ov_idx[[l]] %in% idx]

      # within.x.names[[l]]   <- ov.names[ ov.names %in% ov.names.l[[1]] &
      #                                   ov.names %in% ov.names.x &
      #                                  !ov.names %in% ov.names.l[[2]] ]
      # within.y.names[[l]]   <- ov.names[ ov.names %in% ov.names.l[[1]] &
      #                                  !ov.names %in% ov.names.x &
      #                                  !ov.names %in% ov.names.l[[2]] ]
      # between.x.names[[l]]  <- ov.names[!ov.names %in% ov.names.l[[1]] &
      #                                   ov.names %in% ov.names.x &
      #                                   ov.names %in% ov.names.l[[2]] ]
      # between.y.names[[l]]  <- ov.names[!ov.names %in% ov.names.l[[1]] &
      #                                  !ov.names %in% ov.names.x &
      #                                   ov.names %in% ov.names.l[[2]] ]
      within_x_names[[l]] <- ov_names[within_x_idx[[l]]]
      within_y_names[[l]] <- ov_names[within_y_idx[[l]]]
      between_x_names[[l]] <- ov_names[between_x_idx[[l]]]
      between_y_names[[l]] <- ov_names[between_y_idx[[l]]]
    }
  } else {
    ov_y_idx <- ov_idx
  }

  out <- list(
    ov.names = ov_names, ov.names.x = ov_names_x, # for this group
    cluster = cluster, # clus = clus,
    # per level
    nclusters = nclusters,
    cluster.size = cluster_size, cluster.id = cluster_id,
    cluster.idx = cluster_idx, cluster.sizes = cluster_sizes,
    ncluster.sizes = ncluster_sizes,
    cluster.size.ns = cluster_size_ns,
    ov.idx = ov_idx, ov.x.idx = ov_x_idx, ov.y.idx = ov_y_idx,
    both.idx = both_idx, within.idx = within_idx,
    within.x.idx = within_x_idx, within.y.idx = within_y_idx,
    between.idx = between_idx,
    between.x.idx = between_x_idx, between.y.idx = between_y_idx,
    both.names = both_names, within.names = within_names,
    within.x.names = within_x_names,
    within.y.names = within_y_names,
    between.names = between_names,
    between.x.names = between_x_names,
    between.y.names = between_y_names
  )

  out
}
