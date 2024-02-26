lav_lavaan_step13_vcov_boot <- function(lavoptions, lavmodel, lavsamplestats, lavdata, lavpartable,
                                     lavcache, lavimplied, lavh1, x) {
  # # # # # # # # # # # # # # # #
  # #  13. lavvcov + lavboot # # 
  # # # # # # # # # # # # # # # #
  VCOV <- NULL
  if (lavoptions$se != "none" && lavoptions$se != "external" &&
      lavoptions$se != "twostep" &&
      #(.hasSlot(lavmodel, "nefa") &&
      #    (lavmodel@nefa == 0L ||
      #       (lavmodel@nefa > 0L && lavoptions$rotation == "none") ||
      #       (lavmodel@nefa > 0L && lavoptions$rotation.se == "delta")
      #)
      #) &&
      lavmodel@nx.free > 0L && (attr(x, "converged") ||
                                lavoptions$optim.method == "none")) {
    
    if (lavoptions$verbose) {
      cat("computing VCOV for      se =", lavoptions$se, "...")
    }
    VCOV <- lav_model_vcov(lavmodel        = lavmodel,
                           lavsamplestats  = lavsamplestats,
                           lavoptions      = lavoptions,
                           lavdata         = lavdata,
                           lavpartable     = lavpartable,
                           lavcache        = lavcache,
                           lavimplied      = lavimplied,
                           lavh1           = lavh1)
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
    
  } # VCOV
  
  # extract bootstrap results (if any)
  if (!is.null(attr(VCOV, "BOOT.COEF"))) {
    lavboot <- list()
    lavboot$coef <- attr(VCOV, "BOOT.COEF")
  } else {
    lavboot <- list()
  }
  
  # store VCOV in vcov
  # strip all attributes but 'dim'
  tmp.attr <- attributes(VCOV)
  VCOV1 <- VCOV
  attributes(VCOV1) <- tmp.attr["dim"]
  # store vcov? new in 0.6-6
  if (!is.null(lavoptions$store.vcov) && !is.null(VCOV1)) {
    if (is.logical(lavoptions$store.vcov) && !lavoptions$store.vcov) {
      VCOV1 <- NULL
    }
    if (is.character(lavoptions$store.vcov) &&
        lavoptions$rotation == "none" &&
        lavoptions$store.vcov == "default" &&
        ncol(VCOV1) > 200L) {
      VCOV1 <- NULL
    }
  }
  lavvcov <- list(se = lavoptions$se, information = lavoptions$information,
                  vcov = VCOV1)
  
  # store se in partable
  if (lavoptions$se == "external") {
    if (is.null(lavpartable$se)) {
      lavpartable$se <- lav_model_vcov_se(lavmodel = lavmodel,
                                          lavpartable = lavpartable,
                                          VCOV = NULL, BOOT = NULL)
      warning("lavaan WARNING: se = \"external\" but parameter table does not contain a `se' column")
    }
  } else if (lavoptions$se %in% c("none", "twostep")) {
    # do nothing
  } else {
    lavpartable$se <- lav_model_vcov_se(lavmodel = lavmodel,
                                        lavpartable = lavpartable,
                                        VCOV = VCOV,
                                        BOOT = lavboot$coef)
  }
  return(list(
    lavpartable = lavpartable,
    lavvcov = lavvcov,
    VCOV = VCOV,
    lavmodel = lavmodel,
    lavboot = lavboot
  ))
}
