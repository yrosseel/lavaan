lav_lavaan_step17_lavaan <- function(lavmc = NULL,
                                     timing = NULL,
                                     lavoptions = NULL,
                                     lavpartable = NULL,
                                     lavdata = NULL,
                                     lavsamplestats = NULL,
                                     lavmodel = NULL,
                                     lavcache = NULL,
                                     lavfit = NULL,
                                     lavboot = NULL,
                                     lavoptim = NULL,
                                     lavimplied = NULL,
                                     lavloglik = NULL,
                                     lavvcov = NULL,
                                     lavtest = NULL,
                                     lavh1 = NULL,
                                     lavbaseline = NULL,
                                     start.time0 = NULL) {
  # # # # # # # # # #
  # #  17. lavaan # #
  # # # # # # # # # #

  # stop timer
  # create lavaan object
  # if lavmodel@nefa > 0
  #   compute standardizedSolution and store in ParTable slot in lavaan object
  # if post-checking demanded and converged, execute
  #   lavInspect(lavaan, "post.check")
  #
  timing$total <- (proc.time()[3] - start.time0)
  timing$start.time <- NULL
  lavpta <- lav_partable_attributes(lavpartable)
  lavpartable <- lav_partable_remove_cache(lavpartable)
  lavaan <- new("lavaan", # type_of_slot - where created or modified ?
    # ------------   ------------------------- -
    version = packageDescription("lavaan", fields = "Version"),
    call = lavmc, # match.call - ldw_adapt_match_call
    timing = timing, # list - ldw_add_timing
    Options = lavoptions, # list - options (2) / data (3) / partable (4)
    ParTable = lavpartable,
    # list - partable/bounds/start/model/estoptim/vcovboot/rotation
    pta = lavpta, # list - lav_partable_attributes
    Data = lavdata, # S4 class - data (3)
    SampleStats = lavsamplestats, # S4 class - samplestats (5)
    Model = lavmodel, # S4 class - model (9) / estoptim (11) / vcovboot (13)
    Cache = lavcache, # list - cache (10)
    Fit = lavfit, # S4 class - lav_model_fit (14bis)
    boot = lavboot, # list - vcovboot (13)
    optim = lavoptim, # list - estoptim (11)
    implied = lavimplied, # list - lav_model_implied (12)
    loglik = lavloglik, # list - lav_model_loglik (12)
    vcov = lavvcov, # list - vcovboot (13)
    test = lavtest, # list - test (14)
    h1 = lavh1, # list - h1 (6)
    baseline = lavbaseline, # list - baseline (15)
    internal = list(), # empty list
    external = list() # empty list
  )

  # if model.type = "efa", add standardized solution to partable
  if ((.hasSlot(lavmodel, "nefa")) && (lavmodel@nefa > 0L)) {
    if (lavoptions$verbose) {
      cat("computing standardized solution ... ")
    }
    std <- standardizedSolution(lavaan,
      remove.eq = FALSE,
      remove.ineq = FALSE, remove.def = FALSE
    )
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
    lavaan@ParTable$est.std <- std$est.std
    lavaan@ParTable$se.std <- std$se
  }

  # post-fitting check of parameters
  if (!is.null(lavoptions$check.post) && lavoptions$check.post &&
    lavTech(lavaan, "converged")) {
    if (lavoptions$verbose) {
      cat("post check  ...")
    }
    lavInspect(lavaan, "post.check")
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }

  lavaan
}
