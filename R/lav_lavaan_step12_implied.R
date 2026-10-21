lav_step12_implied <- function(lavoptions = NULL,
                                      lavmodel = NULL) {
  # # # # # # # # # # # #
  # #  12. lavimplied # #
  # # # # # # # # # # # #

  # if lavoptions$implied compute lavimplied via lav_model_implied
  lavimplied <- list()
  # data-defined (":~") parameters: the model-implied moments vary across
  # observations; evaluating them at the representative values of the
  # moderated cells would give misleading numbers, so the slot stays empty
  if (lavoptions$implied && !lav_model_dv_flag(lavmodel)) {
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

lav_step12_loglik <- function(lavoptions = NULL,
                                     lavdata = NULL,
                                     lavsamplestats = NULL,
                                     lavh1 = NULL,
                                     lavimplied = NULL,
                                     lavmodel = NULL,
                                     lavcache = NULL) {
  # # # # # # # # # # # #
  # #  12. lavloglik  # #
  # # # # # # # # # # # #

  # only when missing = "ml" (or "ml.x") and zero coverage
  if (lavoptions$missing %in% c("ml", "ml.x") &&
      lavoptions$model.type == "unrestricted" &&
      length(lavh1) == 0L) {
    lavh1 <- list(implied = lavimplied,
                  logl = list())
  }

  # if lavoptions$loglik compute lavloglik via lav_model_loglik
  lavloglik <- list()
  # (data-defined parameters: no summary-statistics loglikelihood; the
  #  slot stays empty until an external backend supplies the casewise one)
  if (lavoptions$loglik && !lav_model_dv_flag(lavmodel)) {
    if (lav_verbose()) {
      cat("lavloglik   ...")
    }
    lavloglik <- lav_model_loglik(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavh1 = lavh1,
      lavimplied = lavimplied,
      lavmodel = lavmodel,
      lavoptions = lavoptions,
      lavcache = lavcache
    )
    if (lav_verbose()) {
      cat(" done.\n")
    }
  }

  lavloglik
}
