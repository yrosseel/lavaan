lav_lavaan_step12_implied <- function(lavoptions, lavmodel) {
  # # # # # # # # # # # #  
  # #  12. lavimplied # # 
  # # # # # # # # # # # #  
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
  return(lavimplied)
}
lav_lavaan_step12_loglik <- function(lavoptions, lavdata, lavsamplestats, lavimplied, lavmodel) {
  # # # # # # # # # # # # 
  # #  12. lavloglik  # # 
  # # # # # # # # # # # # 
  lavloglik <- list()
  if (lavoptions$loglik) {
    if (lavoptions$verbose) {
      cat("lavloglik   ...")
    }
    lavloglik <- lav_model_loglik(lavdata        = lavdata,
                                  lavsamplestats = lavsamplestats,
                                  lavimplied     = lavimplied,
                                  lavmodel       = lavmodel,
                                  lavoptions     = lavoptions)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }
  return(lavloglik)
}
