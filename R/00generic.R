# for simsem: uses 'inspect' in exportMethods in NAMESPACE
setGeneric("inspect", function(object, ...) standardGeneric("inspect"))

# for blavaan
setGeneric("fitMeasures", 
    function(object, fit.measures = "all", baseline.model = NULL) 
    standardGeneric("fitMeasures"))
setGeneric("fitmeasures",
    function(object, fit.measures = "all", baseline.model = NULL) 
    standardGeneric("fitmeasures"))


# S3 generics
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
