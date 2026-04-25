lav_lavaan_step13_vcov_boot <- function(lavoptions = NULL,
                                        lavmodel = NULL,
                                        lavsamplestats = NULL,
                                        lavdata = NULL,
                                        lavpartable = NULL,
                                        lavcache = NULL,
                                        lavimplied = NULL,
                                        lavh1 = NULL,
                                        x = NULL) {
  # # # # # # # # # # # # # # # #
  # #  13. lavvcov + lavboot # #
  # # # # # # # # # # # # # # # #

  # set VCOV to NULL
  # if lavoptions$se not "none", "external", "twostep" and
  #   lavmodel@nx.free > 0L and x converged or optim.method == "none"
  #   compute VCOV via lav_model_vcov
  # if attribute BOOT.COEFF of VCOV not NULL, store it in lavboot$coef
  # lavvcov <- list(se = lavoptions$se, information = lavoptions$information,
  #   vcov = VCOV1)
  #   where VCOV1 = VCOV without attributes (except dim) or ...
  #                 NULL if lavoptions$store.vcov FALSE or
  #                 store.vcov=="default" and rotation="none"
  # if lavoptions$se == "external"
  #   if lavpartable$se NULL
  #     lavpartable$se <- lav_model_vcov_se(..., VCOV=NULL, BOOT=NULL)
  #       + ** warning **
  # if lavpartable not "external" or "none" or "twostep"
  #     lavpartable$se <- lav_model_vcov_se(...)

  vcov_1 <- NULL # nolint
  if (lavoptions$se != "none" && lavoptions$se != "external" &&
    lavoptions$se != "twostep" &&
    (lavmodel@nefa == 0L ||
     (lavmodel@nefa > 0L && lavoptions$rotation == "none") ||
     (lavmodel@nefa > 0L && lavoptions$rotation.se == "delta")
    ) &&
    lavmodel@nx.free > 0L && (attr(x, "converged") ||
    lavoptions$optim.method == "none")) {
    if (lav_verbose()) {
      cat("computing VCOV for      se =", lavoptions$se, "...")
    }
    # special case: estimator = "IV"
    if (lavoptions$estimator %in% "IV" && !is.null(attr(x, "eqs"))) {
      vcov_1 <- lav_sem_miiv_vcov(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavoptions = lavoptions,
        lavpartable = lavpartable,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        eqs = attr(x, "eqs")
      )
    } else {
      # everything else:
      vcov_1 <- lav_model_vcov( # nolint
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavoptions = lavoptions,
        lavdata = lavdata,
        lavpartable = lavpartable,
        lavcache = lavcache,
        lavimplied = lavimplied,
        lavh1 = lavh1
      )
    }
    if (lav_verbose()) {
      cat(" done.\n")
    }
  } # VCOV

  # extract bootstrap results (if any)
  if (!is.null(attr(vcov_1, "BOOT.COEF"))) {
    lavboot <- list()
    lavboot$coef <- attr(vcov_1, "BOOT.COEF")
  } else {
    lavboot <- list()
  }

  # store VCOV in vcov
  # strip all attributes but 'dim'
  tmp_attr <- attributes(vcov_1)
  vcov1 <- vcov_1 # nolint
  attributes(vcov1) <- tmp_attr["dim"] # nolint
  # store vcov? new in 0.6-6
  if (!is.null(lavoptions$store.vcov) && !is.null(vcov1)) {
    if (is.logical(lavoptions$store.vcov) && !lavoptions$store.vcov) {
      vcov1 <- NULL # nolint
    }
    if (is.character(lavoptions$store.vcov) &&
      lavoptions$rotation == "none" &&
      lavoptions$store.vcov == "default" &&
      ncol(vcov1) > 200L) {
      vcov1 <- NULL # nolint
    }
  }
  lavvcov <- list(
    se = lavoptions$se,
    information = lavoptions$information[1],
    vcov = vcov1
  )

  # store se in partable
  if (lavoptions$se == "external") {
    if (is.null(lavpartable$se)) {
      lavpartable$se <- lav_model_vcov_se(
        lavmodel = lavmodel,
        lavpartable = lavpartable,
        VCOV = NULL, BOOT = NULL
      )
      lav_msg_warn(gettext(
        "se = \"external\" but parameter table does not contain a `se' column"))
    }
  } else if (lavoptions$se %in% c("none", "twostep")) {
    # do nothing
  } else {
    lavpartable$se <- lav_model_vcov_se(
      lavmodel = lavmodel,
      lavpartable = lavpartable,
      VCOV = vcov_1,
      BOOT = lavboot$coef
    )
  }

  list(
    lavpartable = lavpartable,
    lavvcov     = lavvcov,
    VCOV        = vcov_1,
    lavmodel    = lavmodel,
    lavboot     = lavboot
  )
}
