
# handle ov.order = "data" by adding attribute "ovda" to FLAT
lav_partable_ov_from_data <- function(FLAT = NULL, data = NULL,
                                      sample.cov = NULL, slotData = NULL) {
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
        stop("lavaan ERROR: some (observed) variables specified in the model ", 
             "are not found in the data: ", 
             paste(ov.names[idx.missing], collapse=" "))
    }

    # check if the order is the same
    if (!identical(ov.names, ov.names.data)) {
        attr(FLAT, "ovda") <- ov.names.data
    }
    return(FLAT)
}
