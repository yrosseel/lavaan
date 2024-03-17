# lavMultipleImputation: fit the *same* model, on a set of imputed datasets
# YR - 11 July 2016

lavMultipleImputation <-
  function(model = NULL,
           dataList = NULL,
           ndat = length(dataList),
           cmd = "sem",
           ...,
           store.slots = c("partable"),
           FUN = NULL,
           show.progress = FALSE,
           parallel = c("no", "multicore", "snow"),
           ncpus = max(1L, parallel::detectCores() - 1L),
           cl = NULL) {
    # dotdotdot
    dotdotdot <- list()

    # fit multiple times
    fit <- do.call("lavaanList", args = c(list(
      model = model,
      dataList = dataList, ndat = ndat, cmd = cmd,
      store.slots = store.slots, FUN = FUN,
      show.progress = show.progress,
      parallel = parallel, ncpus = ncpus, cl = cl
    ), dotdotdot))

    # flag multiple imputation
    fit@meta$lavMultipleImputation <- TRUE

    fit
  }
