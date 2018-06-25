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

    # add eq.id column
    #if(is.null(partable$eq.id)) {
    #    partable$eq.id <- rep(0, N)
    #}

    # add unco column
    #if(is.null(partable$unco)) {
    #    partable$unco <- partable$free
    #}

    # order them nicely: id lhs op rhs group
    idx <- match(c("id", "lhs","op","rhs", "block","user",
                   "free","ustart","exo","label"),
                 names(partable))

    # order them nicely: id lhs op rhs group
    #idx <- match(c("id", "lhs","op","rhs", "group","user",
    #               "free","ustart","exo","label","eq.id","unco"),
    #             names(partable))
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
