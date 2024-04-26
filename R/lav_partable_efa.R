# generate a parameter table for an EFA model
#
# YR 20 Sept 2022: initial verion
lav_partable_generate_efa <- function(ov.names = NULL,
                                      nfactors = 1L,
                                      meanstructure = FALSE,
                                      varTable = NULL) {
  # currently, we support only a single block (but we plan for more)
  nblocks <- 1L
  ov.names <- list(ov.names)

  # currently, we only support continuous data (but ordered is planned)
  stopifnot(is.null(ordered))

  lhs <- rhs <- op <- character(0)
  block <- free <- integer(0)
  ustart <- numeric(0)

  # create factor names
  lv.names <- paste("f", 1:nfactors, sep = "")

  # block number
  for (b in seq_len(nblocks)) {
    # ov.names for this block
    OV.NAMES <- ov.names[[b]]
    nvar <- length(OV.NAMES)
    nel <- nvar * nfactors

    # get 'ordered' variables from varTable
    categorical <- FALSE
    if (!is.null(varTable)) {
      ov.names.ord <-
        as.character(varTable$name[varTable$type == "ordered"])
      # remove those that do appear in the model syntax
      idx <- which(!ov.names.ord %in% OV.NAMES)
      if (length(idx) > 0L) {
        ov.names.ord <- ov.names.ord[-idx]
      }
      if (length(ov.names.ord) > 0L) {
        ov.names.ord <- OV.NAMES[OV.NAMES %in% ov.names.ord]
        categorical <- TRUE
      }
    }

    # a) factor loadings
    lhs <- c(lhs, rep(lv.names, each = nvar))
    op <- c(op, rep("=~", nel))
    rhs <- c(rhs, rep(OV.NAMES, times = nfactors))
    block <- c(block, rep(b, nel))
    # group <- c(group,  rep(1L,   nel))
    # level <- c(level,  rep(1L,   nel))
    free <- c(free, rep(1L, nel)) # for now
    ustart <- c(ustart, rep(as.numeric(NA), nel))

    # b) ov variances
    lhs <- c(lhs, OV.NAMES)
    op <- c(op, rep("~~", nvar))
    rhs <- c(rhs, OV.NAMES)
    block <- c(block, rep(b, nvar))
    # group <- c(group, rep(1L,   nvar))
    # level <- c(level, rep(1L,   nvar))
    free <- c(free, rep(1L, nvar))
    ustart <- c(ustart, rep(as.numeric(NA), nvar))

    # c) lv variances
    lhs <- c(lhs, lv.names)
    op <- c(op, rep("~~", nfactors))
    rhs <- c(rhs, lv.names)
    block <- c(block, rep(b, nfactors))
    # group <- c(group,  rep(1L,   nfactors))
    # level <- c(level,  rep(1L,   nfactors))
    free <- c(free, rep(0L, nfactors)) # standardized!
    ustart <- c(ustart, rep(1, nfactors))

    # d) lv covariances
    pstar <- nfactors * (nfactors - 1) / 2
    if (pstar > 0L) { # only if more than 1 variable
      tmp <- utils::combn(lv.names, 2)
      lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
      op <- c(op, rep("~~", pstar))
      rhs <- c(rhs, tmp[2, ])
      block <- c(block, rep(b, pstar))
      # group  <- c(group,  rep(g,    pstar))
      # level  <- c(level,  rep(l,    pstar))
      free <- c(free, rep(1L, pstar)) # to be changed...
      ustart <- c(ustart, rep(as.numeric(NA), pstar))
    }

    if (meanstructure) {
      # e) ov means/intercepts
      lhs <- c(lhs, OV.NAMES)
      op <- c(op, rep("~1", nvar))
      rhs <- c(rhs, rep("", nvar))
      block <- c(block, rep(b, nvar))
      # group <- c(group,  rep(1L,   nvar))
      # level <- c(level,  rep(1L,   nvar))
      free <- c(free, rep(1L, nvar))
      ustart <- c(ustart, rep(as.numeric(NA), nvar))

      # f) lv means/intercepts
      lhs <- c(lhs, lv.names)
      op <- c(op, rep("~1", nfactors))
      rhs <- c(rhs, rep("", nfactors))
      block <- c(block, rep(b, nfactors))
      # group <- c(group,  rep(1L,   nfactors))
      # level <- c(level,  rep(1L,   nfactors))
      free <- c(free, rep(0L, nfactors))
      ustart <- c(ustart, rep(0, nfactors))
    } # meanstructure
  } # blocks

  # create LIST
  LIST <- list(
    id = 1:length(lhs),
    lhs = lhs,
    op = op,
    rhs = rhs,
    user = rep(0L, length(lhs)), # all system-generated
    block = block,
    group = rep(1L, length(lhs)),
    # level       = level,
    free = free,
    ustart = ustart,
    exo = rep(0L, length(lhs)),
    label = rep("", length(lhs)),
    efa = rep("", length(lhs))
  )

  # add 'efa' column with a single block string (i.e., "efa")
  LIST$efa[LIST$op == "=~"] <- "efa"

  # take care of EFA constraints
  LIST <- lav_partable_efa_constraints(LIST)

  # free counter
  idx.free <- which(LIST$free > 0)
  LIST$free[idx.free] <- seq_along(idx.free)

  # needed?
  LIST <- lav_partable_complete(LIST)

  LIST
}

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
