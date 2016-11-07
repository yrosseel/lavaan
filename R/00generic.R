# for simsem: uses 'inspect' in exportMethods in NAMESPACE
setGeneric("inspect", function(object, ...) standardGeneric("inspect"))

# for blavaan
setGeneric("fitMeasures", 
    function(object, fit.measures = "all", baseline.model = NULL) 
    standardGeneric("fitMeasures"))
setGeneric("fitmeasures",
    function(object, fit.measures = "all", baseline.model = NULL) 
    standardGeneric("fitmeasures"))

