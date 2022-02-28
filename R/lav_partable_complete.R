# handle bare-minimum partables
# add some additional columns
lav_partable_complete <- function(partable = NULL, start = TRUE) {

    # check if we hava a data.frame
    # if so, check for columns that are 'factor' and convert them to 'character'
    if(is.data.frame(partable)) {
        fac.idx <- sapply(partable, is.factor)
        partable[fac.idx] <- lapply(partable[fac.idx], as.character)
    }

    # check if we have lhs, op, rhs
    stopifnot(!is.null(partable$lhs),
              !is.null(partable$op),
              !is.null(partable$rhs))

    # number of elements
    N <- length(partable$lhs)
    if(!is.data.frame(partable)) {
         # check for equal column length
         nel <- sapply(partable, length)
         short.idx <- which(nel < N)
         long.idx <- which(nel > N)
         if(length(long.idx) > 0L) {
             warning("lavaan WARNING: partable columns have unequal length")
         }
         if(length(short.idx) > 0L) {
             # try to extend them in a 'natural' way
             for(i in short.idx) {
                 too.short <- N - nel[i]
                 if(is.integer(partable[[i]])) {
                     partable[[i]] <- c(partable[[i]],
                                           integer( too.short ))
                 } else if(is.numeric(partable[[i]])) {
                     partable[[i]] <- c(partable[[i]],
                                           numeric( too.short ))
                 } else {
                     partable[[i]] <- c(partable[[i]],
                                           character( too.short ))
                 }
             }
         }
    }

    # create new id column
    #if(is.null(partable$id)) {
        partable$id <- seq_len(N)
    #}

    # add block column
    if(is.null(partable$block)) {
        partable$block <- rep(1L, N)
    } else {
        partable$block <- as.integer(partable$block)
    }

    # add user column
    if(is.null(partable$user)) {
        partable$user <- rep(1L, N)
    } else {
         partable$user <- as.integer( partable$user )
    }

    # add free column
    if(is.null(partable$free)) {
        partable$free <- seq_len(N)
    # 0.6-11: check for simple equality constraints
    #         note: this is perhaps only a subset (eg SAM!) of a larger
    #         table, and we have to renumber the 'free' column
    } else if( is.integer(partable$free) &&
               any(partable$free > 0L)   &&
               !any(partable$op == "==") &&
               !is.null(partable$label)  &&
               !is.null(partable$plabel) &&
               any(duplicated(partable$free[partable$free > 0L])) ) {
        dup.idx <- which(partable$free > 0L & duplicated(partable$free))
        all.idx <- which(partable$free %in% unique(partable$free[dup.idx]))
        eq.LABELS <- unique(partable$free[all.idx])
        eq.id <- integer(length(partable$lhs))
        eq.id[all.idx] <- partable$free[all.idx]
        partable$free[dup.idx] <- 0L
        idx.free <- which(partable$free > 0L)
        partable$free <- rep(0L, N)
        partable$free[idx.free] <- seq_along(idx.free)
        for(eq.label in eq.LABELS) {
            all.idx <- which(eq.id == eq.label)
            ref.idx <- all.idx[1L]
            other.idx <- all.idx[-1L]
            partable$free[other.idx] <- partable$free[ref.idx]
        }
    } else {
        # treat non-zero as 'free'
        free.idx <- which(as.logical(partable$free))
        partable$free <- rep(0L, N)
        if(length(free.idx) > 0L) {
            partable$free[free.idx] <- seq_len(length(free.idx))
        }
    }

    # add ustart column
    if(is.null(partable$ustart)) {
        # do we have something else? start? est?
        if(!is.null(partable$start)) {
            partable$ustart <- as.numeric(partable$start)
        } else if(!is.null(partable$est)) {
            partable$ustart <- as.numeric(partable$est)
        } else {
            partable$ustart <- rep(as.numeric(NA), N)
            non.free <- which(!partable$free)
            if(length(non.free)) {
                partable$ustart[non.free] <- 0
            }
        }
    } else {
        partable$ustart <- as.numeric(partable$ustart)
    }

    # add exo column
    if(is.null(partable$exo)) {
        partable$exo <- rep(0, N)
    } else {
        partable$exo <- as.integer( partable$exo )
    }

    # add label column
    if(is.null(partable$label)) {
        partable$label <- rep("", N)
    } else {
        partable$label <- as.character( partable$label )
    }

    # order them nicely: id lhs op rhs group
    idx <- match(c("id", "lhs","op","rhs", "block","user",
                   "free","ustart","exo","label"),
                 names(partable))
    tmp <- partable[idx]
    partable <- c(tmp, partable[-idx])

    # add start column
    if(start) {
        if(is.null(partable$start)) {
            partable$start <- lav_start(start.method = "simple",
                                        lavpartable = partable)
        }
    }

    partable
}
