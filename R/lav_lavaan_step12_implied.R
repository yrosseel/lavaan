lav_lavaan_step12_implied <- function(lavoptions = NULL,
                                      lavmodel = NULL) {
  # # # # # # # # # # # #
  # #  12. lavimplied # #
  # # # # # # # # # # # #

  # if lavoptions$implied compute lavimplied via lav_model_implied
  lavimplied <- list()
  if (lavoptions$implied) {
    if (lav_verbose()) {
      cat("lavimplied  ...")
    }
    lavimplied <- lav_model_implied(lavmodel)
    if (lav_verbose()) {
      cat(" done.\n")
    }
  }

  lavimplied
}

lav_lavaan_step12_loglik <- function(lavoptions = NULL,
                                     lavdata = NULL,
                                     lavsamplestats = NULL,
                                     lavh1 = NULL,
                                     lavimplied = NULL,
                                     lavmodel = NULL) {
  # # # # # # # # # # # #
  # #  12. lavloglik  # #
  # # # # # # # # # # # #

  # only when missing = "ml" and zero coverage
  if (lavoptions$missing == "ml" && lavoptions$model.type == "unrestricted" &&
      length(lavh1) == 0L) {
    lavh1 <- list(implied = lavimplied,
                  logl = list())
  }

  # if lavoptions$loglik compute lavloglik via lav_model_loglik
  lavloglik <- list()
  if (lavoptions$loglik) {
    if (lav_verbose()) {
      cat("lavloglik   ...")
    }
    lavloglik <- lav_model_loglik(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavh1 = lavh1,
      lavimplied = lavimplied,
      lavmodel = lavmodel,
      lavoptions = lavoptions
    )
    if (lav_verbose()) {
      cat(" done.\n")
    }
  }

  lavloglik
}
