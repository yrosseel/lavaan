# guess number of blocks from a partable
lav_partable_nblocks <- function(partable) {

    if(is.null(partable$block)) {
        nblocks <- 1L
    } else {
        # always integers
        tmp <- partable$block[ partable$block > 0L ] # non-zero only
        nblocks <- length(unique(na.omit(tmp)))      # all combinations
    }

    nblocks
}

# guess number of groups from a partable
lav_partable_ngroups <- function(partable) {

    if(is.null(partable$group)) {
        ngroups <- 1L
    } else if(is.numeric(partable$group)) {
        tmp <- partable$group[ partable$group > 0L ]
        ngroups <- length(unique(na.omit(tmp)))
    } else { # character
        tmp <- partable$group[nchar(partable$group) > 0L]
        ngroups <- length(unique(na.omit(tmp)))
    }

    ngroups
}

# guess number of levels from a partable
lav_partable_nlevels <- function(partable) {

    if(is.null(partable$level)) {
        nlevels. <- 1L
    } else if(is.numeric(partable$level)) {
        tmp <- partable$level[ partable$level > 0L ]
        nlevels. <- length(unique(na.omit(tmp)))
    } else { # character
        tmp <- partable$level[nchar(partable$level) > 0L]
        nlevels. <- length(unique(na.omit(tmp)))
    }

    nlevels.
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
        }

        ndat[b] <- pstar

        # correction for fixed.x?
        if(!conditional.x && fixed.x) {
            ov.names.x <- lav_partable_vnames(partable, "ov.x", block = b)
            nvar.x <- length(ov.names.x)
            pstar.x <- nvar.x * (nvar.x + 1) / 2
            if(meanstructure) {
                pstar.x <- pstar.x + nvar.x
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

