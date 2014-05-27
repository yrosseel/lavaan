# backwards compatibility (for semTools!) - YR 14 Jan 2014 (dev 0.5-16)
#
# 1) computeVY has samplestats argument
# 2) computeObjective has been renamed, see below
# 3) computeTestStatistic has been renamed, see below

computeObjective <- function(object, GLIST=NULL,
                             samplestats=NULL, X = NULL,
                             cache=NULL,
                             estimator="ML", verbose=FALSE, forcePD=TRUE,
                             debug=FALSE) {
    lav_model_objective(lavmodel = object, GLIST = GLIST,
                        lavsamplestats = samplestats,
                        lavdata = NULL, # X??? not always needed
                        lavcache = cache,
                        estimator = estimator, verbose=verbose, forcePD=forcePD,
                        debug=debug)
}


computeTestStatistic <- function(object, partable=NULL, samplestats=NULL,
                                 options=NULL, x=NULL, VCOV=NULL, cache=NULL,
                                 data=NULL, control=list()) {

    lav_model_test(lavmodel = object, lavpartable = partable, 
                   lavsamplestats = samplestats, lavoptions = options,
                   x = x, VCOV = VCOV, lavcache = cache, 
                   lavdata = data, control = control)
}

Fit <- function(partable=NULL, start, model, x=NULL, VCOV=NULL, TEST=NULL) {

    lav_model_fit(lavpartable = partable, start = start, lavmodel = model,
                  x = x, VCOV = VCOV, TEST = TEST)
}

# for simsem: uses 'inspect' in exportMethods in NAMESPACE
# setGeneric("inspect", function(object, ...) standardGeneric("inspect"))
