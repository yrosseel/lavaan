# for blavaan
# TDJ: add "..." to make the generic actually generic, for lavaan.mi objects
setGeneric("fitMeasures",
    function(object, fit.measures = "all", baseline.model = NULL, ...)
    standardGeneric("fitMeasures"))
setGeneric("fitmeasures",
    function(object, fit.measures = "all", baseline.model = NULL, ...)
    standardGeneric("fitmeasures"))


# S3 generics
inspect <- function(object, what = "free", ...) {
    UseMethod("inspect", object)
}

lavInspect <- function(object, what                   = "free",
                               add.labels             = TRUE,
                               add.class              = TRUE,
                               list.by.group          = TRUE,
                               drop.list.single.group = TRUE) {
    UseMethod("lavInspect", object)
}

lavTech <- function(object, what                   = "free",
                            add.labels             = FALSE,
                            add.class              = FALSE,
                            list.by.group          = FALSE,
                            drop.list.single.group = FALSE) {
    UseMethod("lavTech", object)
}
