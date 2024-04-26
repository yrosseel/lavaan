# get missing patterns
lav_data_missing_patterns <- function(Y, sort.freq = FALSE, coverage = FALSE,
                                      Lp = NULL) {
  # handle two-level data
  if (!is.null(Lp)) {
    Y.orig <- Y
    Z <- NULL
    if (length(Lp$between.idx[[2]]) > 0L) {
      Y <- Y[, -Lp$between.idx[[2]], drop = FALSE]
      z.idx <- which(!duplicated(Lp$cluster.idx[[2]]))
      Z <- Y.orig[z.idx, Lp$between.idx[[2]], drop = FALSE]
    }
  }

  # construct TRUE/FALSE matrix: TRUE if value is observed
  OBS <- !is.na(Y)

  # empty cases
  empty.idx <- which(rowSums(OBS) == 0L)

  # pattern of observed values per observation
  case.id <- apply(1L * OBS, 1L, paste, collapse = "")

  # remove empty patterns
  if (length(empty.idx)) {
    case.id.nonempty <- case.id[-empty.idx]
  } else {
    case.id.nonempty <- case.id
  }

  # sort non-empty patterns (from high occurence to low occurence)
  if (sort.freq) {
    TABLE <- sort(table(case.id.nonempty), decreasing = TRUE)
  } else {
    TABLE <- table(case.id.nonempty)
  }

  # unique pattern ids
  pat.id <- names(TABLE)

  # number of patterns
  pat.npatterns <- length(pat.id)

  # case idx per pattern
  pat.case.idx <- lapply(
    seq_len(pat.npatterns),
    function(p) which(case.id == pat.id[p])
  )

  # unique pattern frequencies
  pat.freq <- as.integer(TABLE)

  # first occurrence of each pattern
  pat.first <- match(pat.id, case.id)

  # TRUE/FALSE for each pattern
  pat.obs <- OBS[pat.first, , drop = FALSE] # observed per pattern

  Mp <- list(
    npatterns = pat.npatterns, id = pat.id, freq = pat.freq,
    case.idx = pat.case.idx, pat = pat.obs, empty.idx = empty.idx,
    nel = sum(OBS)
  )

  if (coverage) {
    # FIXME: if we have empty cases, include them in N?
    # no for now
    Mp$coverage <- crossprod(OBS) / sum(pat.freq)
    # Mp$coverage <- crossprod(OBS) / NROW(Y)
  }

  # additional info in we have two-level data
  if (!is.null(Lp)) {
    Mp$j.idx <- lapply(
      seq_len(pat.npatterns),
      function(p) Lp$cluster.idx[[2]][Mp$case.idx[[p]]]
    )
    Mp$j1.idx <- lapply(
      seq_len(pat.npatterns),
      function(p) unique.default(Mp$j.idx[[p]])
    )
    Mp$j.freq <- lapply(
      seq_len(pat.npatterns),
      function(p) as.integer(unname(table(Mp$j.idx[[p]])))
    )

    # between-level patterns
    if (!is.null(Z)) {
      Mp$Zp <- lav_data_missing_patterns(Z,
        sort.freq = FALSE,
        coverage = FALSE, Lp = NULL
      )
    }
  }

  Mp
}

# get response patterns (ignore empty cases!)
lav_data_resp_patterns <- function(Y) {
  # construct TRUE/FALSE matrix: TRUE if value is observed
  OBS <- !is.na(Y)

  # empty cases
  empty.idx <- which(rowSums(OBS) == 0L)

  # removeYempty cases
  if (length(empty.idx) > 0L) {
    Y <- Y[-empty.idx, , drop = FALSE]
  }

  ntotal <- nrow(Y)
  nvar <- ncol(Y)

  # identify, label and sort response patterns
  id <- apply(Y, MARGIN = 1, paste, collapse = "")

  # sort patterns (from high occurence to low occurence)
  TABLE <- sort(table(id), decreasing = TRUE)
  order <- names(TABLE)
  npatterns <- length(TABLE)
  pat <- Y[match(order, id), , drop = FALSE]
  row.names(pat) <- as.character(TABLE)

  # handle NA?
  Y[is.na(Y)] <- -9
  total.patterns <- prod(apply(Y, 2, function(x) length(unique(x))))
  empty.patterns <- total.patterns - npatterns
  # return a list
  # out <- list(nobs=ntotal, nvar=nvar,
  #            id=id, npatterns=npatterns,
  #            order=order, pat=pat)

  # only return pat
  out <- list(
    npatterns = npatterns, pat = pat, total.patterns = total.patterns,
    empty.patterns = empty.patterns
  )

  out
}

# get cluster information
# - cluster can be a vector!
# - clus can contain multiple columns!
lav_data_cluster_patterns <- function(Y = NULL,
                                      clus = NULL, # the cluster ids
                                      cluster = NULL, # the cluster 'names'
                                      multilevel = FALSE,
                                      ov.names = NULL,
                                      ov.names.x = NULL,
                                      ov.names.l = NULL) {
  # how many levels?
  nlevels <- length(cluster) + 1L

  # did we get any data (or is this just for simulateData)
  if (!is.null(Y) && !is.null(clus)) {
    haveData <- TRUE
  } else {
    haveData <- FALSE
  }

  # check clus
  if (haveData) {
    stopifnot(ncol(clus) == (nlevels - 1L), nrow(Y) == nrow(clus))
  }

  cluster.size <- vector("list", length = nlevels)
  cluster.id <- vector("list", length = nlevels)
  cluster.idx <- vector("list", length = nlevels)
  nclusters <- vector("list", length = nlevels)
  cluster.sizes <- vector("list", length = nlevels)
  ncluster.sizes <- vector("list", length = nlevels)
  cluster.size.ns <- vector("list", length = nlevels)
  ov.idx <- vector("list", length = nlevels)
  ov.x.idx <- vector("list", length = nlevels)
  ov.y.idx <- vector("list", length = nlevels)
  both.idx <- vector("list", length = nlevels)
  within.idx <- vector("list", length = nlevels)
  within.x.idx <- vector("list", length = nlevels)
  within.y.idx <- vector("list", length = nlevels)
  between.idx <- vector("list", length = nlevels)
  between.x.idx <- vector("list", length = nlevels)
  between.y.idx <- vector("list", length = nlevels)
  both.names <- vector("list", length = nlevels)
  within.names <- vector("list", length = nlevels)
  within.x.names <- vector("list", length = nlevels)
  within.y.names <- vector("list", length = nlevels)
  between.names <- vector("list", length = nlevels)
  between.x.names <- vector("list", length = nlevels)
  between.y.names <- vector("list", length = nlevels)

  # level-1 is special
  if (haveData) {
    nclusters[[1]] <- NROW(Y)
  }

  # higher levels:
  for (l in 2:nlevels) {
    if (haveData) {
      CLUS <- clus[, (l - 1L)]
      cluster.id[[l]] <- unique(CLUS)
      cluster.idx[[l]] <- match(CLUS, cluster.id[[l]])
      cluster.size[[l]] <- tabulate(cluster.idx[[l]])
      nclusters[[l]] <- length(cluster.size[[l]])
      # check if we have more observations than clusters
      if (nclusters[[1]] == nclusters[[l]]) {
        lav_msg_stop(gettext("every cluster contains only one observation."))
      }
      mean.cluster.size <- mean(cluster.size[[l]])
      if (mean.cluster.size < 1.5) {
        lav_msg_warn(gettextf(
          "mean cluster size is %s. This means that many clusters only
          contain a single observation.", mean.cluster.size))
        }
      cluster.sizes[[l]] <- unique(cluster.size[[l]])
      ncluster.sizes[[l]] <- length(cluster.sizes[[l]])
      cluster.size.ns[[l]] <- as.integer(table(factor(cluster.size[[l]],
        levels = as.character(cluster.sizes[[l]])
      )))
    } else {
      cluster.id[[l]] <- integer(0L)
      cluster.idx[[l]] <- integer(0L)
      cluster.size[[l]] <- integer(0L)
      nclusters[[l]] <- integer(0L)
      cluster.sizes[[l]] <- integer(0L)
      ncluster.sizes[[l]] <- integer(0L)
      cluster.size.ns[[l]] <- integer(0L)
    }
  }

  # for all levels:
  if (multilevel) {
    for (l in 1:nlevels) {
      # index of ov.names for this level
      ov.idx[[l]] <- match(ov.names.l[[l]], ov.names)

      # new in 0.6-12: always preserve the order of ov.idx[[l]]
      idx <- which(ov.names %in% ov.names.l[[1]] &
        ov.names %in% ov.names.l[[2]])
      both.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]

      idx <- which(ov.names %in% ov.names.l[[1]] &
        !ov.names %in% ov.names.l[[2]])
      within.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]
      # backwards compatibility: also store in within.idx[[2]]
      if (l == 2) {
        within.idx[[l]] <- within.idx[[1]]
      }

      idx <- which(!ov.names %in% ov.names.l[[1]] &
        ov.names %in% ov.names.l[[2]])
      between.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]

      # names
      # both.names[[l]]     <- ov.names[ ov.names %in% ov.names.l[[1]] &
      #                                 ov.names %in% ov.names.l[[2]] ]
      # within.names[[l]]   <- ov.names[ ov.names %in% ov.names.l[[1]] &
      #                                !ov.names %in% ov.names.l[[2]] ]
      # between.names[[l]]  <- ov.names[!ov.names %in% ov.names.l[[1]] &
      #                                 ov.names %in% ov.names.l[[2]] ]
      both.names[[l]] <- ov.names[both.idx[[l]]]
      within.names[[l]] <- ov.names[within.idx[[l]]]
      between.names[[l]] <- ov.names[between.idx[[l]]]
    }
  }

  # fixed.x wrt variable index
  if (multilevel && length(ov.names.x) > 0L) {
    for (l in 1:nlevels) {
      # some ov.names.x could be 'splitted', and end up in both.names
      # they should NOT be part ov.x.idx (as they become latent variables)
      idx <- which(ov.names %in% ov.names.x &
        ov.names %in% ov.names.l[[l]] &
        !ov.names %in% unlist(both.names))
      ov.x.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]

      # not any longer, we split them, but still treat them as 'fixed'
      # ov.x.idx[[l]] <- which( ov.names %in% ov.names.x &
      #                        ov.names %in% ov.names.l[[l]] )

      # if some ov.names.x have been 'splitted', and end up in both.names,
      # they should become part of ov.y.idx (despite being exogenous)
      # as they are now latent variables
      idx <- which(ov.names %in% ov.names.l[[l]] &
        !ov.names %in% ov.names.x[!ov.names.x %in% unlist(both.names)])
      ov.y.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]

      # not any longer, ov.x stays ov.x (even if we split)
      # ov.y.idx[[l]] <- which( ov.names %in% ov.names.l[[l]] &
      #                       !ov.names %in% ov.names.x )


      # if(l == 1L) {
      #    next
      # }
      # below, we only fill in the [[2]] element (and higher)

      idx <- which(ov.names %in% ov.names.l[[1]] &
        !ov.names %in% ov.names.l[[2]] &
        ov.names %in% ov.names.x)
      within.x.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]
      # backwards compatibility: also store in within.x.idx[[2]]
      if (l == 2) {
        within.x.idx[[l]] <- within.x.idx[[1]]
      }

      idx <- which(ov.names %in% ov.names.l[[1]] &
        !ov.names %in% ov.names.l[[2]] &
        !ov.names %in% ov.names.x)
      within.y.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]
      # backwards compatibility: also store in within.y.idx[[2]]
      if (l == 2) {
        within.y.idx[[l]] <- within.y.idx[[1]]
      }

      idx <- which(!ov.names %in% ov.names.l[[1]] &
        ov.names %in% ov.names.l[[2]] &
        ov.names %in% ov.names.x)
      between.x.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]

      idx <- which(!ov.names %in% ov.names.l[[1]] &
        ov.names %in% ov.names.l[[2]] &
        !ov.names %in% ov.names.x)
      between.y.idx[[l]] <- ov.idx[[l]][ov.idx[[l]] %in% idx]

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
      within.x.names[[l]] <- ov.names[within.x.idx[[l]]]
      within.y.names[[l]] <- ov.names[within.y.idx[[l]]]
      between.x.names[[l]] <- ov.names[between.x.idx[[l]]]
      between.y.names[[l]] <- ov.names[between.y.idx[[l]]]
    }
  } else {
    ov.y.idx <- ov.idx
  }

  out <- list(
    ov.names = ov.names, ov.names.x = ov.names.x, # for this group
    cluster = cluster, # clus = clus,
    # per level
    nclusters = nclusters,
    cluster.size = cluster.size, cluster.id = cluster.id,
    cluster.idx = cluster.idx, cluster.sizes = cluster.sizes,
    ncluster.sizes = ncluster.sizes,
    cluster.size.ns = cluster.size.ns,
    ov.idx = ov.idx, ov.x.idx = ov.x.idx, ov.y.idx = ov.y.idx,
    both.idx = both.idx, within.idx = within.idx,
    within.x.idx = within.x.idx, within.y.idx = within.y.idx,
    between.idx = between.idx,
    between.x.idx = between.x.idx, between.y.idx = between.y.idx,
    both.names = both.names, within.names = within.names,
    within.x.names = within.x.names,
    within.y.names = within.y.names,
    between.names = between.names,
    between.x.names = between.x.names,
    between.y.names = between.y.names
  )

  out
}
