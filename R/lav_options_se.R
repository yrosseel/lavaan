# in separate file LDW 06/04/2024

# rename names of se values, and check for invalid values
lav_options_check_se <- function(opt = NULL) {
  # unlike test=, se= should be a single character string
  if (length(opt$se) > 1L) opt$se <- opt$se[1]

  # backwards compatibility (0.4 -> 0.5)
  if (opt$se == "first.order") {
    opt$se <- "standard"
    opt$information[1] <- "first.order"
    if (length(opt$information) > 1L &&
        opt$information[2] == "default") {
      opt$information[2] <- "first.order"
    }
  } else if (opt$se == "observed") {
    opt$se <- "standard"
    opt$information[1] <- "observed"
    if (length(opt$information) > 1L &&
        opt$information[2] == "default") {
      opt$information[2] <- "observed"
    }
  } else if (opt$se == "expected") {
    opt$se <- "standard"
    opt$information[1] <- "expected"
    if (length(opt$information) > 1L &&
        opt$information[2] == "default") {
      opt$information[2] <- "expected"
    }
  }

  # handle generic 'robust' (except clustered/multilvel)
  # else if(opt$se == "robust" && !opt$.clustered && !opt$.multilevel) {
  #    if(opt$missing %in% c("ml", "ml.x")) {
  #        opt$se <- "robust.huber.white"
  #    } else if(opt$missing == "two.stage") {
  #        opt$se <- "two.stage"
  #    } else if(opt$missing == "robust.two.stage") {
  #        opt$se <- "robust.two.stage"
  #    } else {
  #        # depends on estimator!
  #        opt$se <- "robust.sem"
  #    }
  # }

  # GLS, NTRLS, FML, UMN
  ok.flag <- TRUE
  if (any(opt$estimator == c("gls", "ntrls", "fml"))) {
    ok.flag <- any(opt$se == c(
      "default", "none", "standard",
      "bootstrap", "external"
    ))
  }

  # WLS, DLS, DWLS, WLSM, WLSMV, WLSMVS, ULS, ULSM, ULSMV, ULSMVS
  else if (any(opt$estimator == c(
    "wls", "dls",
    "dwls", "wlsm", "wlsmv", "wlsmvs",
    "uls", "ulsm", "ulsmv", "ulsmvs"
  ))) {
    ok.flag <- any(opt$se == c(
      "default", "none", "standard",
      "bootstrap", "external",
      "robust", "robust.sem"
    ))
  }

  # PML
  else if (opt$estimator == "pml") {
    ok.flag <- any(opt$se == c(
      "default", "none", "standard",
      "bootstrap", "external",
      "robust.huber.white"
    ))
  }

  # FABIN, GUTTMAN1952, BENTLER1982, ...
  else if (any(opt$estimator == c(
    "fabin2", "fabin3", "mgm"))) {
    ok.flag <- any(opt$se == c("default", "none", "bootstrap", "external"))
  }

  # OTHERS
  else if (any(opt$estimator == c("fml", "mml", "reml"))) {
    ok.flag <- any(opt$se == c("default", "none", "standard", "external"))
  }

  if (!ok.flag) {
    lav_msg_stop(gettextf(
      "invalid value (%1$s) in se= argument for estimator %2$s.",
      opt$se, toupper(opt$estimator)))
  }

  opt
}
