# generic methods
#
# initial version: YR 25/03/2009

setGeneric("computeMuHat", function(object, ...) standardGeneric("computeMuHat"))

setGeneric("computeSigmaHat", function(object, ...) standardGeneric("computeSigmaHat"))

setGeneric("computeObjective", function(object, ...) standardGeneric("computeObjective"))

setGeneric("computeGradient", function(object, ...) standardGeneric("computeGradient"))

#setGeneric("computeExpectedInformation", function(object, ...) standardGeneric("computeExpectedInformation"))

#setGeneric("computeObservedInformation", function(object, ...) standardGeneric("computeObservedInformation"))

#setGeneric("vcov.internal", function(object, ...) standardGeneric("vcov.internal"))

#setGeneric("computeStdErrors", function(object, ...) standardGeneric("computeStdErrors"))

setGeneric("estimateModel", function(object, ...) standardGeneric("estimateModel"))

#setGeneric("estimateVCOV", function(object, ...) standardGeneric("estimateVCOV"))

setGeneric("inspect", function(object, ...) standardGeneric("inspect"))
