# methods
setMethod("show", "lavaanList",
function(object) {
    # show only basic information
    lav_lavaanList_short_summary(object)

})

lav_lavaanList_short_summary <- function(object) {
    cat(sprintf("lavaanList (%s) -- based on %d datasets (%d converged)\n",
                 packageDescription("lavaan", fields="Version"),
                 object@meta$ndat,
                 sum(object@meta$ok)))
}

setMethod("summary", "lavaanList",
function(object, header       = TRUE,
                 nd = 3L) {
    lav_lavaanList_summary(object, header = header, nd = nd)
})

lav_lavaanList_summary <- function(object,
                                   header    = TRUE,
                                   estimates = TRUE,
                                   nd = 3L) {

    if(header) {
        # show only basic information
        lav_lavaanList_short_summary(object)
    }

    if(estimates && "partable" %in% object@meta$store.slots) {
        PE <- parameterEstimates(object, zstat = FALSE, pvalue = FALSE,
                                 ci = FALSE, standardized = FALSE,
                                 add.attributes = TRUE)
        # 'default' is average, so better change the name
        PE$est.av <- PE$est
        PE$est <- NULL
        print(PE, nd = nd)
    } else {
        cat("available slots (per dataset) are:\n")
        print(object@meta$store.slots)
    }
}

setMethod("coef", "lavaanList",
function(object, type = "free", labels = TRUE) {
    lav_lavaanList_coef(object = object, type = type, labels = labels)
})

lav_lavaanList_coef <- function(object, type = "free", labels = TRUE) {
    if("partable" %in% object@meta$store.slots) {
        COF <- sapply(object@ParTableList, "[[", "est")
    } else {
        stop("lavaan ERROR: no ParTable slot stored in lavaanList object")
    }

    if(type == "user" || type == "all") {
        type <- "user"
        idx <- 1:length( object@ParTable$lhs )
    } else if(type == "free") {
        idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
    } else {
        stop("lavaan ERROR: argument `type' must be one of free or user")
    }

    COF <- COF[idx, , drop = FALSE]

    if(labels) {
        rownames(COF) <- lav_partable_labels(object@ParTable, type = type)
    }

    COF
}

