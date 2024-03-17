lav_lavaan_step12_implied <- function(lavoptions = NULL,
                                      lavmodel = NULL) {
  # # # # # # # # # # # #
  # #  12. lavimplied # #
  # # # # # # # # # # # #

  # if lavoptions$implied compute lavimplied via lav_model_implied
  lavimplied <- list()
  if (lavoptions$implied) {
    if (lavoptions$verbose) {
      cat("lavimplied  ...")
    }
    lavimplied <- lav_model_implied(lavmodel)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }

  lavimplied
}

lav_lavaan_step12_loglik <- function(lavoptions = NULL,
                                     lavdata = NULL,
                                     lavsamplestats = NULL,
                                     lavimplied = NULL,
                                     lavmodel = NULL) {
  # # # # # # # # # # # #
  # #  12. lavloglik  # #
  # # # # # # # # # # # #

  # if lavoptions$loglik compute lavloglik via lav_model_loglik
  lavloglik <- list()
  if (lavoptions$loglik) {
    if (lavoptions$verbose) {
      cat("lavloglik   ...")
    }
    lavloglik <- lav_model_loglik(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      lavmodel = lavmodel,
      lavoptions = lavoptions
    )
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }

  lavloglik
}
