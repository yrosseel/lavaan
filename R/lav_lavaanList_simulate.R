# lavSimulate: fit the *same* model, on simulated datasets
# YR - 4 July 2016

lavSimulate <- function(pop.model = NULL, # population model
                        model = NULL, # user model
                        dataFunction = simulateData,
                        dataFunction.args = list(
                          model = pop.model,
                          sample.nobs = 1000L
                        ),
                        ndat = 1000L,
                        cmd = "sem",
                        cmd.pop = "sem",
                        ...,
                        store.slots = c("partable"),
                        FUN = NULL,
                        show.progress = FALSE,
                        store.failed = FALSE,
                        parallel = c("no", "multicore", "snow"),
                        ncpus = max(1L, parallel::detectCores() - 1L),
                        cl = NULL) {
  # dotdotdot
  dotdotdot <- list(...)

  # dotdotdot for fit.pop
  dotdotdot.pop <- dotdotdot
  dotdotdot.pop$verbose <- FALSE
  dotdotdot.pop$debug <- FALSE
  dotdotdot.pop$data <- NULL
  dotdotdot.pop$sample.cov <- NULL

  # 'fit' population model without data, to get 'true' parameters
  fit.pop <- do.call(cmd.pop,
    args = c(list(model = pop.model), dotdotdot.pop)
  )

  # check model object
  if (is.null(model)) {
    model <- fit.pop@ParTable
  }

  # per default, use 'true' values as starting values
  if (is.null(dotdotdot$start)) {
    dotdotdot$start <- fit.pop
  }

  # no warnings during/after the simulations
  # dotdotdot$warn <- FALSE

  # generate simulations
  fit <- do.call("lavaanList", args = c(list(
    model = model,
    dataFunction = dataFunction,
    dataFunction.args = dataFunction.args,
    ndat = ndat, cmd = cmd,
    store.slots = store.slots, FUN = FUN,
    show.progress = show.progress,
    store.failed = store.failed,
    parallel = parallel, ncpus = ncpus, cl = cl
  ), dotdotdot))

  # flag this is a simulation
  fit@meta$lavSimulate <- TRUE

  # NOTE!!!
  # if the model != pop.model, we may need to 'reorder' the
  # 'true' parameters, so they correspond to the 'model' parameters
  p2.id <- lav_partable_map_id_p1_in_p2(
    p1 = fit@ParTable,
    p2 = fit.pop@ParTable,
    stopifnotfound = FALSE
  )
  est1 <- fit@ParTable$est
  na.idx <- which(is.na(p2.id))
  if (length(na.idx) > 0L) {
    lav_msg_warn(gettext(
      "some estimated parameters were not mentioned in the population model;
      partable user model idx = ", lav_msg_view(na.idx, "none")))

    # replace NA by '1' (override later!)
    p2.id[na.idx] <- 1L
  }
  est.pop <- fit.pop@ParTable$est[p2.id]

  # by default, the 'unknown' population values are set to 0.0
  if (length(na.idx) > 0L) {
    est.pop[na.idx] <- 0
  }

  # store 'true' parameters in meta$est.true
  fit@meta$est.true <- est.pop

  fit
}
