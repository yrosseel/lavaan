
# handle ov.order = "data" by add dummy rhs/op/lhs entries to trick lavNames()
lav_partable_ov_from_data <- function(FLAT = NULL, data = NULL,
                                      sample.cov = NULL, slotData = NULL) {

    # store original FLAT
    FLAT.orig <- FLAT
    ATTR <- attributes(FLAT)

    # current model-based ov.names
    ov.names <- lav_partable_vnames(FLAT, type = "ov")

    # get data-based ov.names
    DATA.NAMES <- NULL
    if(!is.null(data)) {
        DATA.NAMES <- names(data)
    } else if(!is.null(sample.cov)) {
        # multiple group/blocks?
        if(is.list(sample.cov)) {
            DATA.NAMES <- unique(unlist(lapply(sample.cov, colnames)))
            if(is.null(DATA.NAMES)) {
                # try again with rows
                DATA.NAMES <- unique(unlist(lapply(sample.cov, rownames)))
            }
        } else {
            DATA.NAMES <- colnames(sample.cov)
            if(is.null(DATA.NAMES)) {
                # try again with rows
                DATA.NAMES <- rownames(sample.cov)
            }
        }
    } else if(!is.null(slotData)) {
        DATA.NAMES <- unique(unlist(slotData@ov.names))
    }

    if(is.null(DATA.NAMES) || length(DATA.NAMES) == 0L) {
        stop("lavaan ERROR: could not find variable names in data/sample.cov")
    }

    # extract needed ov.names in the same order as the data
    ov.names.data <- DATA.NAMES[ DATA.NAMES %in% ov.names ]

    # check if we have all of them
    if(length(ov.names.data) != length(ov.names)) {
        idx.missing <- which(!(ov.names %in% ov.names.data))
        stop("lavaan ERROR: some (observed) variables specified in the model are not found in the datat: ", paste(ov.names[idx.missing], collapse=" "))
    }

    # check if the order is the same
    if(identical(ov.names, ov.names.data)) {
        # nothing to do!
        return(FLAT.orig)
    }

    # ok, do we have a regular FLAT object?
    #if(!is.null(FLAT$mod.idx) && !is.null(FLAT$fixed)) {
    #    attr(FLAT, "ov.names.data") <- ov.names.data
    #    return(FLAT)
    #}
    # if FLAT is full/partial partable, append "rhs da lhs" entries

    # nvar
    nvar <- length(ov.names.data)

    # add all ov.names.data to lhs/op/rhs
    FLAT <- as.list(FLAT)
    FLAT$lhs <- c(FLAT$lhs, ov.names.data)
    FLAT$op  <- c(FLAT$op,  rep("da", nvar))
    FLAT$rhs <- c(FLAT$rhs, ov.names.data)

    # enlarge all other list elements
    n.old <- length(FLAT.orig$lhs)
    n.new <- n.old + nvar
    FLAT <- lapply(FLAT, function(x) {
                       if(length(x) != n.new) {
                           if(inherits(x, "character")) {
                               x <- c(x, rep("", nvar))
                           } else if(inherits(x, "integer")) {
                               x <- c(x, rep(0L, nvar))
                           } else if(inherits(x, "numeric")) {
                               x <- c(x, rep(0,  nvar))
                           } else {
                               stop("lavaan ERROR: unknown class [",
                                    class(x), "] in FLAT object")
                           }
                       }
                       x
                   })


    # add attributes
    attributes(FLAT) <- ATTR

    # return
    FLAT
}
