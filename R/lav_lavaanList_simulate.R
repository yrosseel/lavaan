# lavSimulate: fit the *same* model, on simulated datasets
# YR - 4 July 2016

lavSimulate <- function(pop.model     = NULL,             # population model
                        model         = NULL,             # user model
                        dataFunction  = simulateData,
                        dataFunction.args = list(model = pop.model,
                                                 sample.nobs = 1000L),   
                        ndat          = 1000L,
                        cmd           = "sem",
                        ...,
                        store.slots   = c("partable"),
                        FUN           = NULL,
                        show.progress = FALSE,
                        parallel      = c("no", "multicore", "snow"),
                        ncpus         = max(1L, parallel::detectCores() - 1L),
                        cl            = NULL) {

    if(is.null(model)) {
        model <- pop.model
    }

    # 'fit' population model, to get 'true' parameters
    fit.orig <- do.call(cmd, args = list(model = pop.model))

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
