# handle bare-minimum partables
# add some additional columns
lav_partable_complete <- function(lavpartable = NULL, start = TRUE) {

    # check if we hava a data.frame
    # if so, check for columns that are 'factor' and convert them to 'character'
    if(is.data.frame(lavpartable)) {
        fac.idx <- sapply(lavpartable, is.factor)
        lavpartable[fac.idx] <- lapply(lavpartable[fac.idx], as.character)
    }

    # check if we have lhs, op, rhs
    stopifnot(!is.null(lavpartable$lhs),
              !is.null(lavpartable$op),
              !is.null(lavpartable$rhs))
        
    # number of elements
    N <- length(lavpartable$lhs)
    if(!is.data.frame(lavpartable)) {
         # check for equal column length
         nel <- sapply(lavpartable, length)
         short.idx <- which(nel < N)
         long.idx <- which(nel > N)
         if(length(long.idx) > 0L) {
             warning("lavaan WARNING: partable columns have unequal length")
         }
         if(length(short.idx) > 0L) {
             # try to extend them in a 'natural' way
             for(i in short.idx) {
                 too.short <- N - nel[i]
                 if(is.integer(lavpartable[[i]])) {
                     lavpartable[[i]] <- c(lavpartable[[i]], 
                                           integer( too.short ))
                 } else if(is.numeric(lavpartable[[i]])) {
                     lavpartable[[i]] <- c(lavpartable[[i]], 
                                           numeric( too.short ))
                 } else {
                     lavpartable[[i]] <- c(lavpartable[[i]],
                                           character( too.short ))
                 }
             }
         }
    }

    # add id column
    if(is.null(lavpartable$id)) {
        lavpartable$id <- seq_len(N)
    }

    # add group column
    if(is.null(lavpartable$group)) {
        lavpartable$group <- rep(1L, N)
    }

    # add user column
    if(is.null(lavpartable$user)) {
        lavpartable$user <- rep(1L, N)
    }

    # add free column
    if(is.null(lavpartable$free)) {
        lavpartable$free <- seq_len(N)
    } else {
        # treat non-zero as 'free'
        free.idx <- which(as.logical(lavpartable$free))
        lavpartable$free <- rep(0L, N)
        if(length(free.idx) > 0L) {
            lavpartable$free[free.idx] <- seq_len(length(free.idx))
        }
    }

    # add ustart column
    if(is.null(lavpartable$ustart)) {
        # do we have something else? start? est?
        if(!is.null(lavpartable$start)) {
            lavpartable$ustart <- lavpartable$start
        } else if(!is.null(lavpartable$est)) {
            lavpartable$ustart <- lavpartable$est
        } else {
            lavpartable$ustart <- rep(as.numeric(NA), N)
            non.free <- which(!lavpartable$free)
            if(length(non.free)) {
                lavpartable$ustart[non.free] <- 0
            }
        }
    }

    # add exo column
    if(is.null(lavpartable$exo)) {
        lavpartable$exo <- rep(0, N)
    }

    # add label column
    if(is.null(lavpartable$label)) {
        lavpartable$label <- rep("", N)
    }

    # add eq.id column 
    #if(is.null(lavpartable$eq.id)) {
    #    lavpartable$eq.id <- rep(0, N)
    #}

    # add unco column
    #if(is.null(lavpartable$unco)) {
    #    lavpartable$unco <- lavpartable$free
    #}

    # order them nicely: id lhs op rhs group
    idx <- match(c("id", "lhs","op","rhs", "group","user",
                   "free","ustart","exo","label"),
                 names(lavpartable))

    # order them nicely: id lhs op rhs group
    #idx <- match(c("id", "lhs","op","rhs", "group","user",
    #               "free","ustart","exo","label","eq.id","unco"), 
    #             names(lavpartable))
    tmp <- lavpartable[idx]
    lavpartable <- c(tmp, lavpartable[-idx])

    # add start column
    if(start) {
        if(is.null(lavpartable$start)) {
            lavpartable$start <- lav_start(start.method = "simple",
                                           lavpartable = lavpartable)
        }
    }

    lavpartable
}
