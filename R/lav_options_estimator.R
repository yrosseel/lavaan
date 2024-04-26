# lav_options specific per estimator in separate functions LDW 06/04/2024

lav_options_est_ml <- function(opt) {
  # ML and friends: MLF, MLM, MLMV, MLMVS, MLR                     ####
  # se
  if (opt$se == "bootstrap" &&
      opt$estimator %in% c("mlf", "mlm", "mlmv", "mlmvs", "mlr")) {
    lav_msg_stop(gettext("use ML estimator for bootstrap"))
  } else if (opt$se == "default") {
    if (opt$estimator %in% c("ml", "mlf")) {
      opt$se <- "standard"
    } else if (opt$estimator %in% c("mlm", "mlmv", "mlmvs")) {
      opt$se <- "robust.sem"
    } else if (opt$estimator == "mlr") {
      opt$se <- "robust.huber.white"
    }
  } else if (opt$se == "robust") {
    if (opt$missing %in% c("ml", "ml.x")) {
      opt$se <- "robust.huber.white"
    } else if (opt$missing == "two.stage") { # needed?
      opt$se <- "two.stage"
    } else if (opt$missing == "robust.two.stage") { # needed?
      opt$se <- "robust.two.stage"
    } else {
      opt$se <- "robust.sem"
    }
  }
  # information
  if (opt$estimator == "mlf") {
    if (opt$information[1] == "default") {
      opt$information[1] <- "first.order"
    }
    if (length(opt$information) > 1L &&
        opt$information[2] == "default") {
      opt$information[2] <- "first.order"
    }
  }
  # test
  if (!opt$test[1] == "none") {
    if (opt$estimator %in% c("ml", "mlf")) {
      if (opt$test[1] == "default") {
        opt$test <- "standard"
      } # else {
      #    opt$test <- union("standard", opt$test)
      # }
    } else if (opt$estimator == "mlm") {
      if (opt$test[1] == "default") {
        opt$test <- "satorra.bentler"
      } else {
        opt$test <- union("satorra.bentler", opt$test)
      }
    } else if (opt$estimator == "mlmv") {
      if (opt$test[1] == "default") {
        opt$test <- "scaled.shifted"
      } else {
        opt$test <- union("scaled.shifted", opt$test)
      }
    } else if (opt$estimator == "mlmvs") {
      if (opt$test[1] == "default") {
        opt$test <- "mean.var.adjusted"
      } else {
        opt$test <- union("mean.var.adjusted", opt$test)
      }
    }
  }
  opt
}

lav_options_est_gls <- function(opt) {
  # GLS                                                            ####
  # FIXME: catch categorical, clustered, ...
  # se
  if (opt$se == "default") {
    opt$se <- "standard"
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "standard"
  }
  bad.idx <- which(!opt$test %in% c(
    "standard", "none",
    "browne.residual.nt", # == standard
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  if (length(bad.idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is GLS: %s.",
      paste(opt$test[bad.idx], collapse = " ")))
  }
  # missing
  opt$missing <- "listwise" # also pairwise?
  opt
}

lav_options_est_ntrls <- function(opt) {
  # NTRLS (experimental)                                           ####
  # optim.gradient
  opt$optim.gradien <- "numerical"
  # se
  if (opt$se == "default") {
    opt$se <- "standard"
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "standard"
  }
  bad.idx <- which(!opt$test %in% c(
    "standard", "none",
    "browne.residual.nt",
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  if (length(bad.idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is NTRLS: %s.",
      paste(opt$test[bad.idx], collapse = " ")))
  }
  # missing
  opt$missing <- "listwise"
  opt
}

lav_options_est_catml <- function(opt) {
  # catML (experimental)                                           ####
  # optim.gradient
  opt$optim.gradient <- "numerical" # for now
  # force correlation = TRUE, and categorical = FALSE
  opt$correlation <- TRUE
  opt$.categorical <- FALSE # we 'pretend' to have continuous data!
  # se
  if (opt$se == "default") {
    opt$se <- "robust.sem" # for now
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "satorra.bentler"
  }
  # missing
  if (opt$missing %in% c("listwise", "pairwise")) {
    # nothing to do
  } else if (opt$missing == "default") {
    opt$missing <- "listwise"
  } else {
    lav_msg_stop(gettext(
      "missing argument should be listwise or pairwise if estimator is catML"))
  }
  opt
}

lav_options_est_wls <- function(opt) {
  # WLS                                                            ####
  # se
  if (opt$se == "default") {
    opt$se <- "standard"
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "standard"
  }
  bad.idx <- which(!opt$test %in% c(
    "standard", "none",
    "browne.residual.nt",
    "browne.residual.nt.model",
    "browne.residual.adf", # == standard
    "browne.residual.adf.model"
  ))
  if (length(bad.idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is WLS: %s.",
      paste(opt$test[bad.idx], collapse = " ")))
  }
  # missing
  # opt$missing <- "listwise" (could be pairwise)
  opt
}

lav_options_est_dls <- function(opt) {
  # DLS                                                            ####
  # se
  if (opt$se == "default") {
    opt$se <- "robust.sem"
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "satorra.bentler"
  }
  bad.idx <- which(!opt$test %in% c(
    "standard", "none",
    "satorra.bentler",
    "browne.residual.nt", # == standard
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  if (length(bad.idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is DLS: %s.",
      paste(opt$test[bad.idx], collapse = " ")))
  }
  # always include "satorra.bentler"
  if (opt$test[1] %in% c(
    "browne.residual.nt", "browne.residual.adf",
    "browne.residual.nt.model",
    "browne.residual.adf.model"
  )) {
    opt$test <- union("satorra.bentler", opt$test)
  }
  # missing
  opt$missing <- "listwise"
  # estimator.args
  if (is.null(opt$estimator.args)) {
    opt$estimator.args <- list(
      dls.a = 1.0, dls.GammaNT = "model",
      dls.FtimesNmin1 = FALSE
    )
  } else {
    if (is.null(opt$estimator.args$dls.a)) {
      opt$estimator.args$dls.a <- 1.0
    } else {
      stopifnot(is.numeric(opt$estimator.args$dls.a))
      if (opt$estimator.args$dls.a < 0.0 ||
          opt$estimator.args$dls.a > 1.0) {
        lav_msg_stop(gettext(
          "dls.a value in estimator.args must be between 0 and 1."))
      }
    }
    if (is.null(opt$estimator.args$dls.GammaNT)) {
      opt$estimator.args$dls.GammaNT <- "model"
    } else {
      stopifnot(is.character(opt$estimator.args$dls.GammaNT))
      opt$estimator.args$dls.GammaNT <-
        tolower(opt$estimator.args$dls.GammaNT)
      if (!opt$estimator.args$dls.GammaNT %in% c("sample", "model")) {
        lav_msg_stop(gettextf(
          "dls.GammaNT value in estimator.args must be either %s.",
          lav_msg_view(c("sample", "model"), log.sep = "or")))
      }
    }
    if (is.null(opt$estimator.args$dls.FtimesNminus1)) {
      opt$estimator.args$dls.FtimesNminus1 <- FALSE
    } else {
      stopifnot(is.logical(opt$estimator.args$dls.FtimesNminus1))
    }
  }
  if (opt$estimator.args$dls.GammaNT == "sample") {
    if (opt$optim.method %in% c("nlminb", "gn")) {
      # nothing to do
    } else if (opt$optim.method == "default") {
      opt$optim.method <- "gn"
    } else {
      lav_msg_stop(gettext(
        "optim.method must be either nlminb or gn if estimator is DLS."))
    }
  } else {
    if (opt$optim.method %in% c("gn")) {
      # nothing to do
    } else if (opt$optim.method == "default") {
      opt$optim.method <- "gn"
    } else if (opt$optim.method == "nlminb") {
      opt$optim.gradient <- "numerical"
    } else {
      lav_msg_stop(gettext(
        "optim.method must be either nlminb or gn if estimator is DLS."))
    }
  }
  opt
}

lav_options_est_dwls <- function(opt) {
  # DWLS, WLSM, WLSMV, WLSMVS                                      ####
  # new in 0.6-17: if !categorical, give a warning
  if (!opt$.categorical) {
    lav_msg_warn(gettextf(
      "estimator %s is not recommended for continuous data.
      Did you forget to set the ordered= argument?",
      dQuote(lav_options_estimatorgroup(opt$estimator))))
  }
  # se
  if (opt$se == "bootstrap" &&
      opt$estimator %in% c("wlsm", "wlsmv", "wlsmvs")) {
    lav_msg_stop(gettext("use (D)WLS estimator for bootstrap"))
  } else if (opt$se == "default") {
    if (opt$estimator == "dwls") {
      opt$se <- "standard"
    } else {
      opt$se <- "robust.sem"
    }
  } else if (opt$se == "robust") {
    opt$se <- "robust.sem"
  }
  # test
  if (!opt$test[1] == "none") {
    if (opt$estimator == "dwls") {
      if (opt$test[1] == "default") {
        opt$test <- "standard"
      } # else {
      #  opt$test <- union("standard", opt$test)
      # }
    } else if (opt$estimator == "wlsm") {
      if (opt$test[1] == "default") {
        opt$test <- "satorra.bentler"
      } else {
        opt$test <- union("satorra.bentler", opt$test)
      }
    } else if (opt$estimator == "wlsmv") {
      if (opt$test[1] == "default") {
        opt$test <- "scaled.shifted"
      } else {
        opt$test <- union("scaled.shifted", opt$test)
      }
    } else if (opt$estimator == "wlsmvs") {
      if (opt$test[1] == "default") {
        opt$test <- "mean.var.adjusted"
      } else {
        opt$test <- union("mean.var.adjusted", opt$test)
      }
    }
  }
  opt
}

lav_options_est_uls <- function(opt) {
  # ULS, ULSM, ULSMV, ULSMVS                                       ####
  # se
  if (opt$se == "bootstrap" &&
      opt$estimator %in% c("ulsm", "ulsmv", "ulsmvs")) {
    lav_msg_stop(gettext("use ULS estimator for bootstrap"))
  } else if (opt$se == "default") {
    if (opt$estimator == "uls") {
      opt$se <- "standard"
    } else {
      opt$se <- "robust.sem"
    }
  } else if (opt$se == "robust") {
    opt$se <- "robust.sem"
  }
  # test
  if (!opt$test[1] == "none") {
    if (opt$estimator == "uls") {
      if (opt$test[1] == "default") {
        opt$test <- "standard"
      } # else {
      #  opt$test <- union("standard", opt$test)
      # }
    } else if (opt$estimator == "ulsm") {
      if (opt$test[1] == "default") {
        opt$test <- "satorra.bentler"
      } else {
        opt$test <- union("satorra.bentler", opt$test)
      }
    } else if (opt$estimator == "ulsmv") {
      if (opt$test[1] == "default") {
        opt$test <- "scaled.shifted"
      } else {
        opt$test <- union("scaled.shifted", opt$test)
      }
    } else if (opt$estimator == "ulsmvs") {
      if (opt$test[1] == "default") {
        opt$test <- "mean.var.adjusted"
      } else {
        opt$test <- union("mean.var.adjusted", opt$test)
      }
    }
  }
  opt
}

lav_options_est_pml <- function(opt) {
  # PML                                                            ####
  # se
  if (opt$se == "default") {
    opt$se <- "robust.huber.white"
  }

  # information
  opt$information[1] <- "observed"
  if (length(opt$information) > 1L &&
      opt$information[2] == "default") {
    opt$information[2] <- "observed"
  }
  if (length(opt$observed.information) > 1L &&
      opt$observed.information[2] == "default") {
    opt$observed.information[2] <- "hessian"
  }

  # test
  if (length(opt$test) > 1L) {
    lav_msg_stop(gettext(
      "only a single test statistic is allow when estimator is PML."))
  }
  if (!opt$test[1] == "none") {
    opt$test <- "mean.var.adjusted"
  }
  opt
}

lav_options_est_fml <- function(opt) {
  # FML - UMN                                                      ####
  # optim.gradient
  opt$optim.gradient <- "numerical"
  # se
  if (opt$se == "default") {
    opt$se <- "standard"
  }
  # information
  opt$information[1] <- "observed"
  if (length(opt$information) > 1L &&
      opt$information[2] == "default") {
    opt$information[2] <- "observed"
  }
  # test
  if (!opt$test[1] == "none") {
    opt$test <- "standard"
  }
  opt
}

lav_options_est_reml <- function(opt) {
  # REML                                                           ####
  # se
  if (opt$se == "default") {
    opt$se <- "standard"
  }
  # information
  opt$information[1] <- "observed"
  if (length(opt$information) > 1L &&
      opt$information[2] == "default") {
    opt$information[2] <- "observed"
  }
  # test
  if (!opt$test[1] == "none") {
    opt$test <- "standard"
  }
  # missing
  opt$missing <- "listwise"
  opt
}

lav_options_est_mml <- function(opt) {
  # MML                                                            ####
  # se
  if (opt$se == "default") {
    opt$se <- "standard"
  }
  # information
  opt$information[1] <- "observed"
  opt$meanstructure <- TRUE
  if (length(opt$information) > 1L &&
      opt$information[2] == "default") {
    opt$information[2] <- "observed"
  }
  # test
  opt$test <- "none"
  # link
  if (opt$link == "default") {
    # opt$link <- "logit"
    opt$link <- "probit"
  } else if (opt$link %in% c("logit", "probit")) {
    # nothing to do
  } else {
    lav_msg_stop(gettext("link must be `logit' or `probit'"))
  }
  # parameterization
  if (opt$parameterization == "default") {
    opt$parameterization <- "mml"
  } else {
    lav_msg_stop(gettext(
      "parameterization argument is ignored if estimator = MML"))
  }
  opt
}

lav_options_est_fabin <- function(opt) {
  # FABIN, MULTIPLE-GROUP-METHOD (MGM), BENTLER1982, ...           ####
  # experimental, for cfa or sam only
  # se
  if (opt$se == "default") {
    opt$se <- "none"
  }
  # bounds
  if (!is.null(opt$bounds) && opt$bounds == "default" &&
      length(opt$optim.bounds) == 0L) {
    opt$bounds <- "standard"
  }
  # test
  if (opt$test == "default") {
    opt$test <- "none" # for now
  }
  # missing
  opt$missing <- "listwise" # for now (until we have two-stage working)
  # options for fabin
  if (lav_options_estimatorgroup(opt$estimator) %in% c("FABIN2", "FABIN3")) {
    if (is.null(opt$estimator.args)) {
      opt$estimator.args <- list(thetapsi.method = "GLS")
    } else {
      if (is.null(opt$estimator.args$thetapsi.method)) {
        opt$estimator.args$thetapsi.method <- "GLS"
      } else {
        opt$estimator.args$thetapsi.method <-
          toupper(opt$estimator.args$thetapsi.method)
        if (opt$estimator.args$thetapsi.method %in% c(
          "ULS",
          "GLS", "WLS", "ULS.ML", "GLS.ML", "WLS.ML"
        )) {
          if (opt$estimator.args$thetapsi.method == "WLS") {
            opt$estimator.args$thetapsi.method <- "GLS"
          }
          if (opt$estimator.args$thetapsi.method == "WLS.ML") {
            opt$estimator.args$thetapsi.method <- "GLS.ML"
          }
        } else {
          lav_msg_stop(gettextf(
            "unknown value for estimator.args$thetapsi.method option: %s.",
            opt$estimator.args$thetapsi.method))
        }
      }
    }
  }
  # options for Bentler
  if (lav_options_estimatorgroup(opt$estimator) == "BENTLER1982") {
    if (is.null(opt$estimator.args)) {
      opt$estimator.args <- list(GLS = FALSE, quadprog = FALSE)
    } else {
      if (is.null(opt$estimator.args$GLS)) {
        opt$estimator.args$GLS <- FALSE
      }
      if (is.null(opt$estimator.args$quadprog)) {
        opt$estimator.args$quadprog <- FALSE
      }
    }
  }
  # options for guttman1952 multiple group method
  if (lav_options_estimatorgroup(opt$estimator) == "MGM") {
    if (is.null(opt$estimator.args)) {
      opt$estimator.args <- list(
        psi.mapping = FALSE,
        quadprog = FALSE
      )
    } else {
      if (is.null(opt$estimator.args$psi.mapping)) {
        opt$estimator.args$psi.mapping <- FALSE
      }
      if (is.null(opt$estimator.args$quadprog)) {
        opt$estimator.args$quadprog <- FALSE
      }
    }
  }
  # brute-force override
  opt$optim.method <- "noniter"
  opt$start <- "simple"
  opt
}

lav_options_est_miiv <- function(opt) {
  # MIIV-2SLS and friends                                          ####
  # se
  if (opt$se == "default") {
    opt$se <- "none" # for now
  }
  # bounds
  if (!is.null(opt$bounds) && opt$bounds == "default" &&
      length(opt$optim.bounds) == 0L) {
    opt$bounds <- "standard"
  }
  # test
  if (opt$test == "default") {
    opt$test <- "none" # for now
  }
  # missing
  opt$missing <- "listwise" # for now
  # estimator options
  if (is.null(opt$estimator.args)) {
    opt$estimator.args <- list(method = "2SLS")
  } else {
    if (is.null(opt$estimator.args$method)) {
      opt$estimator.args$method <- "2SLS"
    }
  }
  # brute-force override
  opt$optim.method <- "noniter"
  opt$start <- "simple"
  opt
}

lav_options_est_none <- function(opt) {
  # NONE                                                           ####
  # se
  if (opt$se == "default") {
    opt$se <- "none"
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "none"
  }
  opt
}