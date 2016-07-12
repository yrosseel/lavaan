# lavSimulate: fit the *same* model, on simulated datasets
# YR - 4 July 2016

lavSimulate <- function(pop.model     = NULL,             # population model
                        model         = NULL,             # user model
                        dataFunction  = simulateData,
                        dataFunction.args = list(model = pop.model,
                                                 sample.nobs = 1000L),   
                        ndat          = 1000L,
                        cmd           = "sem",
                        cmd.pop       = "sem",
                        ...,
                        store.slots   = c("partable"),
                        FUN           = NULL,
                        show.progress = FALSE,
                        parallel      = c("no", "multicore", "snow"),
                        ncpus         = max(1L, parallel::detectCores() - 1L),
                        cl            = NULL) {

    # 'fit' population model, to get 'true' parameters
    fit.orig <- do.call(cmd.pop, args = list(model = pop.model))

    # create (empty) 'model' object
    if(is.null(model)) {
        model <- fit.orig
    } else {
        # this is more tricky!! we must 'embed' the model in the original
        # model, otherwise, the order of the parameters will be different!!

        # TODO
        stop("model argument not supported yet")
    }

    # dotdotdot
    dotdotdot <- list()

    # per default, use 'true' values as starting values
    if(is.null(dotdotdot$start)) {
        dotdotdot$start = fit.orig
    }

    # generate simulations
    fit <- do.call("lavaanList", args = c(list(model = model, 
                   dataFunction = dataFunction,
                   dataFunction.args = dataFunction.args,
                   ndat = ndat, cmd = cmd,
                   store.slots = store.slots, FUN = FUN,
                   show.progress = show.progress,
                   parallel = parallel, ncpus = ncpus, cl = cl), dotdotdot))

    # store 'true' parameters in meta$est.true
    fit@meta$lavSimulate <- TRUE
    fit@meta$est.true <- fit.orig@ParTable$est

    fit
}
