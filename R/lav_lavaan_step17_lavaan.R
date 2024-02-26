lav_lavaan_step17_lavaan <- function(mc, timing, lavoptions, lavpartable,
                                   lavpta, lavdata, lavsamplestats, lavmodel, lavcache, lavfit, lavboot,
                                   lavoptim, lavimplied, lavloglik, lavvcov, lavtest, lavh1, lavbaseline,
                                   start.time0) {
  # # # # # # # # # # 
  # #  17. lavaan # # 
  # # # # # # # # # # 
  
  # stop timer
  timing$total <- (proc.time()[3] - start.time0)
  timing$start.time <- NULL
  
  lavaan <- new("lavaan",                           # type_of_slot - where created or modified ?
                                                    # ------------   ------------------------- -
                version      = as.character(packageVersion("lavaan")),
                call         = mc,                  # match.call - ldw_adapt_match_call
                timing       = timing,              # list - ldw_add_timing
                Options      = lavoptions,          # list - options (2) / data (3) / partable (4)
                ParTable     = lavpartable,         # list - partable/bounds/start/model/estoptim/vcovboot/rotation
                pta          = lavpta,              # list - lav_partable_attributes
                Data         = lavdata,             # S4 class - data (3)
                SampleStats  = lavsamplestats,      # S4 class - samplestats (5)
                Model        = lavmodel,            # S4 class - model (9) / estoptim (11) / vcovboot (13)
                Cache        = lavcache,            # list - cache (10)
                Fit          = lavfit,              # S4 class - lav_model_fit (14bis)
                boot         = lavboot,             # list - vcovboot (13)
                optim        = lavoptim,            # list - estoptim (11)
                implied      = lavimplied,          # list - lav_model_implied (12)
                loglik       = lavloglik,           # list - lav_model_loglik (12)
                vcov         = lavvcov,             # list - vcovboot (13)
                test         = lavtest,             # list - test (14)
                h1           = lavh1,               # list - h1 (6)
                baseline     = lavbaseline,         # list - baseline (15)
                internal     = list(),              # empty list
                external     = list()               # empty list
  )
  
  # if model.type = "efa", add standardized solution to partable
  if ((.hasSlot(lavmodel, "nefa")) && (lavmodel@nefa > 0L)) {
    if (lavoptions$verbose) {
      cat("computing standardized solution ... ")
    }
    STD <- standardizedSolution(lavaan, remove.eq = FALSE,
                                remove.ineq = FALSE, remove.def = FALSE)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
    lavaan@ParTable$est.std <- STD$est.std
    lavaan@ParTable$se.std  <- STD$se
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
  return(lavaan)
}
