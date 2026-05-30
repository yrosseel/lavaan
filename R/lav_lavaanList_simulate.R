# lavSimulate: fit the *same* model, on simulated datasets
# YR - 4 July 2016: initial version
# YR - 15 Oct 2024: add iseed (similar to lavBootstrap)
# YR - 26 Oct 2024: rm pop.model, add est.true argument

lavSimulate <- function(model = NULL, # user model
                        data_function = lav_data_simulate_old,
                        data_function_args = list(),
                        est_true = NULL,
                        ndat = 1000L,
                        cmd = "sem",
                        ...,
                        store_slots = c("partable"),
                        fun = NULL,
                        show_progress = FALSE,
                        store_failed = FALSE,
                        parallel = c("no", "multicore", "snow"),
                        ncpus = max(1L, parallel::detectCores() - 1L),
                        cl = NULL,
                        iseed = NULL) {
  # dotdotdot
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot)

  # dotdotdot for fit.pop
  dotdotdot_pop <- dotdotdot
  dotdotdot_pop$verbose <- FALSE
  dotdotdot_pop$debug <- FALSE
  dotdotdot_pop$data <- NULL
  dotdotdot_pop$sample.cov <- NULL

  # 'fit' model without data, check 'true' parameters
  fit_pop <- do.call(cmd,
    args = c(list(model = model), dotdotdot_pop)
  )

  # check est.true= argument
  stopifnot(!missing(est_true), is.numeric(est_true),
            length(fit_pop@ParTable$lhs) == length(est_true))
  fit_pop@ParTable$est <- est_true
  fit_pop@ParTable$start <- est_true

  # per default, use 'true' values as starting values
  if (is.null(dotdotdot$start)) {
    dotdotdot$start <- fit_pop
  }

  # no warnings during/after the simulations
  # add 'warn = FALSE' to args

  # generate simulations
  fit <- do.call("lavaanList", args = c(list(
    model = model,
    dataFunction = data_function,
    dataFunction.args = data_function_args,
    ndat = ndat, cmd = cmd,
    store.slots = store_slots, FUN = fun,
    show.progress = show_progress,
    store.failed = store_failed,
    parallel = parallel, ncpus = ncpus,
    cl = cl, iseed = iseed
  ), dotdotdot))

  # flag this is a simulation
  fit@meta$lavSimulate <- TRUE

  # store 'true' parameters in meta$est.true
  fit@meta$est.true <- est_true

  fit
}
