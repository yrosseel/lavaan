# lavSimulate: fit the *same* model, on simulated datasets
# YR - 4 July 2016: initial version
# YR - 15 Oct 2024: add iseed (similar to bootstrapLavaan)
# YR - 26 Oct 2024: rm pop.model, add est.true argument

lavSimulate <- function(model = NULL, # user model
                        dataFunction = simulateData,
                        dataFunction.args = list(),
                        est.true = NULL,
                        ndat = 1000L,
                        cmd = "sem",
                        ...,
                        store.slots = c("partable"),
                        FUN = NULL,
                        show.progress = FALSE,
                        store.failed = FALSE,
                        parallel = c("no", "multicore", "snow"),
                        ncpus = max(1L, parallel::detectCores() - 1L),
                        cl = NULL,
                        iseed = NULL) {
  # dotdotdot
  dotdotdot <- list(...)

  # dotdotdot for fit.pop
  dotdotdot.pop <- dotdotdot
  dotdotdot.pop$verbose <- FALSE
  dotdotdot.pop$debug <- FALSE
  dotdotdot.pop$data <- NULL
  dotdotdot.pop$sample.cov <- NULL

  # 'fit' model without data, check 'true' parameters
  fit.pop <- do.call(cmd,
    args = c(list(model = model), dotdotdot.pop)
  )

  # check est.true= argument
  stopifnot(!missing(est.true), is.numeric(est.true),
            length(fit.pop@ParTable$lhs) == length(est.true))
  fit.pop@ParTable$est <- est.true
  fit.pop@ParTable$start <- est.true

  # per default, use 'true' values as starting values
  if (is.null(dotdotdot$start)) {
    dotdotdot$start <- fit.pop
  }

  # no warnings during/after the simulations
  # add 'warn = FALSE' to args

  # generate simulations
  fit <- do.call("lavaanList", args = c(list(
    model = model,
    dataFunction = dataFunction,
    dataFunction.args = dataFunction.args,
    ndat = ndat, cmd = cmd,
    store.slots = store.slots, FUN = FUN,
    show.progress = show.progress,
    store.failed = store.failed,
    parallel = parallel, ncpus = ncpus,
    cl = cl, iseed = iseed
  ), dotdotdot))

  # flag this is a simulation
  fit@meta$lavSimulate <- TRUE

  # store 'true' parameters in meta$est.true
  fit@meta$est.true <- est.true

  fit
}
