
# what are the block values (not necessarly 1..nb)
lav_partable_block_values <- function(partable) {

    if(is.null(partable$block)) {
        block.values <- 1L
    } else {
        # always integers
        tmp <- partable$block[ partable$block > 0L ] # non-zero only
        block.values <- unique(na.omit(tmp)) # could be, eg, '2' only
    }

    block.values
}

# guess number of blocks from a partable
lav_partable_nblocks <- function(partable) {
    length( lav_partable_block_values(partable) )
}

# what are the group values (not necessarily integers)
lav_partable_group_values <- function(partable) {

    if(is.null(partable$group)) {
        group.values <- 1L
    } else if(is.numeric(partable$group)) {
        tmp <- partable$group[ partable$group > 0L ]
        group.values <- unique(na.omit(tmp))
    } else { # character
        tmp <- partable$group[nchar(partable$group) > 0L]
        group.values <- unique(na.omit(tmp))
    }

    group.values
}

# guess number of groups from a partable
lav_partable_ngroups <- function(partable) {
    length( lav_partable_group_values(partable) )
}

# what are the level values (not necessarily integers)
lav_partable_level_values <- function(partable) {

    if(is.null(partable$level)) {
        level.values <- 1L
    } else if(is.numeric(partable$level)) {
        tmp <- partable$level[ partable$level > 0L ]
        level.values <- unique(na.omit(tmp))
    } else { # character
        tmp <- partable$level[nchar(partable$level) > 0L]
        level.values <- unique(na.omit(tmp))
    }

    level.values
}

# guess number of levels from a partable
lav_partable_nlevels <- function(partable) {
    length( lav_partable_level_values(partable) )
}






# number of sample statistics per block
lav_partable_ndat <- function(partable) {

    # global
    meanstructure <- any(partable$op == "~1")
    fixed.x <- any(partable$exo > 0L & partable$free == 0L)
    conditional.x <- any(partable$exo > 0L & partable$op == "~")
    categorical <- any(partable$op == "|")
    if(categorical) { 
        meanstructure <- TRUE
    }

    # blocks
    nblocks <- lav_partable_nblocks(partable)
    nlevels <- lav_partable_nlevels(partable)
    ndat <- integer(nblocks)

    for(b in seq_len(nblocks)) {

        # how many observed variables in this block?
        if(conditional.x) {
            ov.names <- lav_partable_vnames(partable, "ov.nox", block = b)
        } else {
            ov.names <- lav_partable_vnames(partable, "ov", block = b)
        }
        nvar <- length(ov.names)

        # pstar
        pstar <- nvar*(nvar+1)/2
        if(meanstructure) {
            pstar <- pstar + nvar

            # no meanstructure if within level, except ov.x which is not
            # decomposed
            if(nlevels > 1L && b == 1L) {
                pstar <- pstar - nvar

                ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
                ov.names.x2 <- unlist(lav_partable_vnames(partable, "ov.x", 
                                      block = seq_len(nblocks)[-b]))
                ov.names.x <- ov.names.x[ !ov.names.x %in% ov.names.x2 ]
                if(!fixed.x && length(ov.names.x) > 0L) {
                    pstar <- pstar + length(ov.names.x)
                }
            }
        }

        ndat[b] <- pstar

        # correction for fixed.x?
        if(!conditional.x && fixed.x) {
            ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
            nvar.x <- length(ov.names.x)
            pstar.x <- nvar.x * (nvar.x + 1) / 2
            if(meanstructure) {
                if(nlevels > 1L && b == 1L) {
                    # do nothing, they are already removed
                } else {
                    pstar.x <- pstar.x + nvar.x
                }
            }
            ndat[b] <- ndat[b] - pstar.x
        }

        # correction for ordinal data?
        if(categorical) {
            ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
            nexo     <- length(ov.names.x)
            ov.ord   <- lav_partable_vnames(partable, "ov.ord", block = b)
            nvar.ord <- length(ov.ord)
            th       <- lav_partable_vnames(partable, "th", block = b)
            nth      <- length(th)
            # no variances
            ndat[b] <- ndat[b] - nvar.ord
            # no means
            ndat[b] <- ndat[b] - nvar.ord
            # but additional thresholds
            ndat[b] <- ndat[b] + nth
            # add slopes
            ndat[b] <- ndat[b] + (nvar * nexo)
        }

        # correction for conditional.x not categorical
        if(conditional.x && !categorical) {
            ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
            nexo <- length(ov.names.x)
            # add slopes
            ndat[b] <- ndat[b] + (nvar * nexo)
        }

        # correction for group proportions?
        group.idx <- which(partable$lhs == "group" &
                           partable$op == "%" &
                           partable$block == b)
        if(length(group.idx) > 0L) {
            # ndat <- ndat + (length(group.idx) - 1L) # G - 1 (sum to one)
            ndat[b] <- ndat[b] + 1L # poisson: each cell a parameter
        }
    } # blocks

    # sum over all blocks
    sum(ndat)
}

# total number of free parameters
lav_partable_npar <- function(partable) {

    # we only assume non-zero values
    npar <- length( which(partable$free > 0L) )
    npar
}

# global degrees of freedom: ndat - npar
lav_partable_df <- function(partable) {

    npar <- lav_partable_npar(partable)
    ndat <- lav_partable_ndat(partable)

    # degrees of freedom
    df <- ndat - npar

    as.integer(df)
}

# check order of covariances: we only fill the upper.tri
# therefore, we 'switch' lhs & rhs if they appear in the wrong order
lav_partable_covariance_reorder <- function(partable, 
                                            ov.names = NULL,
                                            lv.names = NULL) {

    # shortcut
    cov.idx <- which(partable$op == "~~" & partable$lhs != partable$rhs)
    if(length(cov.idx) == 0L) {
        # nothing to do
        return(partable)
    }

    # get names
    if(is.null(ov.names)) {
        ov.names <- lav_partable_vnames(partable, "ov")
    } else {
        ov.names <- unlist(ov.names)
    }
    if(is.null(lv.names)) {
        lv.names <- lav_partable_vnames(partable, "lv")
    } else {
        lv.names <- unlist(lv.names)
    }
    lv.ov.names <- c(lv.names, ov.names)

    # identify wrong ordering
    lhs.idx <- match(partable$lhs[ cov.idx ], lv.ov.names)
    rhs.idx <- match(partable$rhs[ cov.idx ], lv.ov.names)
    swap.idx <- which(lhs.idx > rhs.idx)

    if(length(swap.idx) == 0L) {
        # nothing to do
        return(partable)
    }

    # swap!
    tmp <- partable$lhs[ swap.idx ]
    partable$lhs[ swap.idx ] <- partable$rhs[ swap.idx ]
    partable$rhs[ swap.idx ] <- tmp
  
    partable
}


