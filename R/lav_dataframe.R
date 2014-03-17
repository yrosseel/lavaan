# data.frame utilities, or how to avoid copying (parts) of the data
# Y.R. 11 April 2013

# this is to replace sapply(frame, function(x) class(x)[1])
# try (in R3.0.0):
# N <- 100000
# frame <- data.frame(a=factor(sample(1:5, size=N, replace=TRUE)), 
#                     b=factor(sample(1:5, size=N, replace=TRUE)),
#                     c=rnorm(N))
# system.time(replicate(1000, sapply(frame, function(x) class(x)[1])))
# #   user  system elapsed 
# #  1.223   0.000   1.222
# system.time(replicate(1000, lav_dataframe_check_vartype(frame)))
# #   user  system elapsed 
# #  0.093   0.000   0.093
lav_dataframe_check_vartype <- function(frame = NULL, ov.names = character(0)) {
    if(missing(ov.names)) {
        var.idx <- seq_len(ncol(frame))
    } else {
        var.idx <- match(unlist(ov.names, use.names = FALSE), names(frame))
    }
    nvar <- length(var.idx)
    out <- character(nvar)
    for(i in seq_len(nvar)) {
        out[i] <- class(frame[[var.idx[i]]])[1L]
        # watch out for matrix type with 1 column
        if(out[i] == "matrix" && ncol(frame[[var.idx[i]]]) == 1L) {
            out[i] <- "numeric"
        }
    }
    out
}

# check if any of the variables in frame are ordered or not
lav_dataframe_check_ordered <- function(frame = NULL, ov.names = character(0)) {
    if(missing(ov.names)) {
        var.idx <- seq_len(ncol(frame))
    } else {
        var.idx <- match(unlist(ov.names, use.names = FALSE), names(frame))
    }
    nvar <- length(var.idx)
    for(i in seq_len(nvar)) {
        if(class(frame[[var.idx[i]]])[1L] == "ordered")
            return(TRUE)
    }
    FALSE
}

# construct vartable, but allow 'ordered/factor' argument to intervene
# we do NOT change the data.frame
lav_dataframe_vartable <- function(frame = NULL, ov.names = NULL, 
                                   ov.names.x = NULL,
                                   ordered = NULL,
                                   factor = NULL,
                                   as.data.frame. = FALSE) {

    if(missing(ov.names)) {
        var.names <- names(frame)
    } else {
        ov.names <- unlist(ov.names, use.names=FALSE)
        ov.names.x <- unlist(ov.names.x, use.names=FALSE)
        var.names <- unique(c(ov.names, ov.names.x))
    }
    nvar <- length(var.names)
    var.idx <- match(var.names, names(frame))


    nobs <- integer(nvar)
    type <- character(nvar)
    user <- integer(nvar)
    exo  <- ifelse(var.names %in% ov.names.x, 1L, 0L)
    mean <- numeric(nvar); var <- numeric(nvar)
    nlev <- integer(nvar); lnam <- character(nvar)
    for(i in seq_len(nvar)) {
        x <- frame[[var.idx[i]]]
        
        type.x <- class(x)[1L]

        # correct for matrix with 1 column
        if(type.x == "matrix" && ncol(x) == 1L) {
            type.x <- "numeric"
        }

        # correct for integers
        if(type.x == "integer") {
            type.x <- "numeric"
        }


        # handle ordered/factor
        if(!is.null(ordered) && var.names[i] %in% ordered) {
            type.x <- "ordered"
            lev <- sort(unique(x)) # we assume integers!
            nlev[i] <- length(lev)
            lnam[i] <- paste(lev, collapse="|")   
            user[i] <- 1L
        } else if(!is.null(factor) && var.names[i] %in% factor) {
            type.x <- "factor"
            lev <- sort(unique(x)) # we assume integers!
            nlev[i] <- length(lev)
            lnam[i] <- paste(lev, collapse="|")
            user[i] <- 1L
        } else {
            nlev[i] <- nlevels(x)
            lnam[i] <- paste(levels(x), collapse="|")
        }

        type[i] <- type.x
        nobs[i] <- sum(!is.na(x))
        mean[i] <- ifelse(type.x == "numeric", mean(x, na.rm=TRUE), 
                          as.numeric(NA))
        var[i]  <- ifelse(type.x == "numeric", var(x, na.rm=TRUE), 
                          as.numeric(NA))
    }

    VAR <- list(name=var.names, idx=var.idx, nobs=nobs, type=type, exo=exo,
                user=user, mean=mean, var=var, nlev=nlev, lnam=lnam)

    if(as.data.frame.) {
        VAR <- as.data.frame(VAR, stringsAsFactors=FALSE,
                             row.names=1:length(VAR$name))
        class(VAR) <- c("lavaan.data.frame", "data.frame")
    }

    VAR
}

