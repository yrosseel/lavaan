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
  # FIXME: catch clustered, ...
  if (opt$.categorical) {
    lav_msg_stop(gettext(
      "ordered categorical data is not supported when estimator is GLS."))
  }
  # two-stage missing data: estimation uses the (EM) saturated moments, while
  # the SEs use the robust.sem sandwich (with the two-stage NACOV) and the
  # scaled (satorra.bentler) test; se/test/missing are set in lav_options()
  two_stage <- any(opt$missing == c("two.stage", "robust.two.stage"))
  if (two_stage) {
    return(opt)
  }
  # se
  if (opt$se == "default") {
    # GLS uses a normal-theory weight matrix, so the 'standard' se/test are
    # only valid for normal, i.i.d. data. With sampling weights the moments
    # follow the sampling design, so default to the (ADF) sandwich se and the
    # ADF-based test instead (cf. ULS).
    if (isTRUE(opt$.sampling.weights)) {
      opt$se <- "robust.sem"
    } else {
      opt$se <- "standard"
    }
  } else if (opt$se == "robust") {
    opt$se <- "robust.sem"
  }
  # test
  if (opt$test[1] == "default") {
    if (isTRUE(opt$.sampling.weights)) {
      opt$test <- "browne.residual.adf"
    } else {
      opt$test <- "standard"
    }
  }
  bad_idx <- which(!opt$test %in% c(
    "standard", "none",
    "browne.residual.nt", # == standard
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  if (length(bad_idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is GLS: %s.",
      paste(opt$test[bad_idx], collapse = " ")))
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
  bad_idx <- which(!opt$test %in% c(
    "standard", "none",
    "browne.residual.nt",
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  if (length(bad_idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is NTRLS: %s.",
      paste(opt$test[bad_idx], collapse = " ")))
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
  # two-stage missing data: see lav_options_est_gls()
  if (any(opt$missing == c("two.stage", "robust.two.stage"))) {
    return(opt)
  }
  # se
  if (opt$se == "default") {
    opt$se <- "standard"
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "standard"
  }
  bad_idx <- which(!opt$test %in% c(
    "standard", "none",
    "browne.residual.nt",
    "browne.residual.nt.model",
    "browne.residual.adf", # == standard
    "browne.residual.adf.model"
  ))
  if (length(bad_idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is WLS: %s.",
      paste(opt$test[bad_idx], collapse = " ")))
  }
  # missing
  # opt$missing <- "listwise" (could be pairwise)
  opt
}

lav_options_est_dls <- function(opt) {
  # DLS                                                            ####
  # no categorical
  if (opt$.categorical) {
    lav_msg_stop(gettext(
      "ordered categorical data is not supported when estimator is DLS."))
  }
  # se
  if (opt$se == "default") {
    opt$se <- "robust.sem"
  }
  # test
  if (opt$test[1] == "default") {
    opt$test <- "satorra.bentler"
  }
  bad_idx <- which(!opt$test %in% c(
    "standard", "none",
    "satorra.bentler",
    "browne.residual.nt", # == standard
    "browne.residual.nt.model",
    "browne.residual.adf",
    "browne.residual.adf.model"
  ))
  if (length(bad_idx) > 0L) {
    lav_msg_stop(gettextf(
      "invalid value(s) in test= argument when estimator is DLS: %s.",
      paste(opt$test[bad_idx], collapse = " ")))
  }
  # always include "satorra.bentler"
  if (opt$test[1] %in% c(
    "browne.residual.nt", "browne.residual.adf",
    "browne.residual.nt.model",
    "browne.residual.adf.model"
  )) {
    opt$test <- union("satorra.bentler", opt$test)
  }
  # missing (two-stage missing data keeps the EM moments; see
  # lav_options_est_gls())
  if (!any(opt$missing == c("two.stage", "robust.two.stage"))) {
    opt$missing <- "listwise"
  }
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
          lav_msg_view(c("sample", "model"), log_sep = "or")))
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
  # changed in 0.7-1: catch this earlier
  #if (!opt$.categorical) {
  #  lav_msg_warn(gettextf(
  #    "estimator %s is not recommended for continuous data.
  #    Did you forget to set the ordered= argument?",
  #    dQuote(lav_options_estimatorgroup(opt$estimator))))
  #}
  # se
  if (opt$se == "bootstrap" &&
      opt$estimator %in% c("wlsm", "wlsmv", "wlsmvs")) {
    lav_msg_stop(gettext("use (D)WLS estimator for bootstrap"))
  } else if (opt$se == "default") {
    if (opt$estimator == "dwls" && !opt$.categorical) {
      # opt$se <- "standard"
      # with continuous data, DWLS follows the same defaults as ULS: with
      # sampling weights, use the ADF-based se (and, below, the ADF test)
      if (isTRUE(opt$.sampling.weights)) {
        opt$se <- "robust.sem"
      } else {
        opt$se <- "robust.sem.nt" # new in 0.6-21
      }
    } else {
      opt$se <- "robust.sem"
    }
  } else if (opt$se == "robust") {
    opt$se <- "robust.sem"
  }
  # test
  if (!opt$test[1] == "none") {
    if (opt$estimator == "dwls") {
      if (opt$test[1] == "default" && !opt$.categorical) {
        if (opt$se == "robust.sem") { # user-specified?
          opt$test <- "browne.residual.adf" # new in 0.7-1
          opt$standard.test <- "browne.residual.adf"
          opt$scaled.test <- "browne.residual.adf"
        } else {
          opt$test <- "browne.residual.nt" # new in 0.6-21
          opt$standard.test <- "browne.residual.nt"
          opt$scaled.test <- "browne.residual.nt"
        }
      } else if (opt$test[1] == "default" && opt$.categorical) {
        opt$test <- "standard" # bad choice!
      } else {
        opt$test <- union("standard", opt$test)
      }
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
  # two-stage missing data: see lav_options_est_gls()
  if (any(opt$missing == c("two.stage", "robust.two.stage"))) {
    return(opt)
  }
  # se
  if (opt$se == "bootstrap" &&
      opt$estimator %in% c("ulsm", "ulsmv", "ulsmvs")) {
    lav_msg_stop(gettext("use ULS estimator for bootstrap"))
  } else if (opt$se == "default") {
    if (opt$estimator == "uls" && !opt$.categorical) {
      #opt$se <- "standard"
      # with sampling weights, the (ADF) full Gamma is computed, so we can
      # default to the ADF-based robust SE (and, below, the ADF-based test)
      # instead of the normal-theory versions
      if (isTRUE(opt$.sampling.weights)) {
        opt$se <- "robust.sem"
      } else {
        opt$se <- "robust.sem.nt" # new in 0.6-21
      }
    } else {
      opt$se <- "robust.sem"
    }
  } else if (opt$se == "robust") {
    opt$se <- "robust.sem"
  }
  # test
  if (!opt$test[1] == "none") {
    if (opt$estimator == "uls") {
      if (opt$test[1] == "default" && !opt$.categorical) {
        if (opt$se == "robust.sem") { # user-specified?
          opt$test <- "browne.residual.adf" # new in 0.7-1
          opt$standard.test <- "browne.residual.adf"
          opt$scaled.test <- "browne.residual.adf"
        } else {
          opt$test <- "browne.residual.nt" # new in 0.6-21
          opt$standard.test <- "browne.residual.nt"
          opt$scaled.test <- "browne.residual.nt"
        }
      } else if (opt$test[1] == "default" && opt$.categorical) {
        opt$test <- "standard"
      } else {
        opt$test <- union("standard", opt$test)
      }
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
      "only a single test statistic is allowed when estimator is PML."))
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
        zero.after.efa = TRUE,
        psi.mapping = FALSE,
        quadprog = FALSE
      )
    } else {
      if (is.null(opt$estimator.args$zero.after.efa)) {
        opt$estimator.args$zero.after.efa <- TRUE
      }
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

lav_options_est_iv <- function(opt) {
  # (MI)IV-2SLS and friends                                          ####

  # the IV estimator.args use snake_case names (iv_samplestats, iv_vcov_stage1,
  # ...); normalize any user-supplied dot.case / CamelCase names so both styles
  # are accepted during the transition to snake_case
  if (is.list(opt$estimator.args) && length(opt$estimator.args) > 0L) {
    opt$estimator.args <- lav_snake_case(opt$estimator.args)
  }

  # brute-force override
  opt$optim.method <- "noniter"
  opt$marker.int.zero <- TRUE

  # treat simple equality constraints (eg equal loadings) as shared
  # free parameters; this lets the equation-by-equation estimator honor
  # them via a pooled (system) solve -- see lav_sem_miiv_pool_directed()
  opt$ceq.simple <- TRUE

  # missing data
  # - two.stage / robust.two.stage: use the (saturated) EM moments for point
  #   estimation, with standard errors corrected for the EM uncertainty (see
  #   lav_sem_miiv_vcov()); requires the sample-statistics path
  # - anything else falls back to listwise deletion
  two_stage <- any(opt$missing == c("two.stage", "robust.two.stage"))
  if (two_stage && isTRUE(opt$.categorical)) {
    lav_msg_warn(gettext(
      "[IV] two-stage missing-data handling is not available for categorical
       data; using listwise deletion instead."))
    two_stage <- FALSE
  }
  if (!two_stage) {
    opt$missing <- "listwise" # for now
  }

  # se: the general missing = two.stage handling (in lav_options()) sets se to
  # (robust.)two.stage, but the IV estimator uses its own SE machinery and
  # reads opt$missing to pick the two-stage moment covariance; reset to
  # "standard" here
  if (any(opt$se == c("default", "two.stage", "robust.two.stage"))) {
    opt$se <- "standard" # for now
  }
  # bounds
  if (!is.null(opt$bounds) && opt$bounds == "default" &&
      length(opt$optim.bounds) == 0L) {
    opt$bounds <- "standard"
  }
  # test (the general two.stage handling forces satorra.bentler; reset to the
  # IV default Browne residual test). With missing data the sample-based
  # version is not available, so use the model-based variant.
  if ((length(opt$test) == 1L && opt$test == "default") || two_stage) {
    if (opt$.categorical) {
      opt$test <- if (two_stage) {
        "browne.residual.adf.model"
      } else {
        "browne.residual.adf" # always sample-based
      }
    } else {
      opt$test <- if (two_stage) {
        "browne.residual.nt.model"
      } else {
        "browne.residual.nt" # sample-based (especially for baseline)
      }                       # model-based Sigma is here diagonal!
    }
  }
  opt$standard.test <- opt$test[1]

  # fixed.x
  opt$fixed.x <- FALSE # for now

  # sample.icov not needed
  opt$sample.icov <- FALSE

  # no loglik
  opt$loglik <- FALSE

  # no rescale of vcov
  opt$sample.cov.rescale <- FALSE

  # estimator options
  if (is.null(opt$estimator.args)) {
    # create default list
    opt$estimator.args <- list(iv_method = "2SLS",
                               iv_samplestats = TRUE,
                               iv_varcov_method = "RLS",
                               iv_sargan = TRUE,
                               iv_sargan_adjust = "none",
                               iv_weak = "warn",
                               iv_weak_threshold = 10,
                               iv_vcov_stage1 = "lm.vcov.dfres",
                               iv_vcov_stage2 = "h2",
                               iv_vcov_gamma_modelbased = TRUE,
                               iv_vcov_jack_numerical = FALSE,
                               iv_vcov_jaca_numerical = FALSE,
                               iv_vcov_jacb_numerical = FALSE)
  } else {
    if (is.null(opt$estimator.args$iv_method)) {
      opt$estimator.args$iv_method <- "2SLS"
    } else if (!opt$estimator.args$iv_method %in% "2SLS") {
      lav_msg_stop(gettext("iv_method should be 2SLS (for now)."))
    }
    if (is.null(opt$estimator.args$iv_samplestats)) {
      opt$estimator.args$iv_samplestats <- TRUE
    }
    # [[ ]] exact match: $ would partial-match iv_sargan to iv_sargan_adjust
    if (is.null(opt$estimator.args[["iv_sargan"]])) {
      opt$estimator.args$iv_sargan <- TRUE
    }
    # multiple-comparison adjustment for the per-equation Sargan p-values
    if (is.null(opt$estimator.args[["iv_sargan_adjust"]])) {
      opt$estimator.args$iv_sargan_adjust <- "none"
    } else if (!opt$estimator.args[["iv_sargan_adjust"]] %in%
               stats::p.adjust.methods) {
      lav_msg_stop(gettextf(
        "iv_sargan_adjust should be one of %s.",
        lav_msg_view(stats::p.adjust.methods, log_sep = "or")))
    }
    # NOTE: use [[ ]] (exact match); $ would partial-match iv_weak to
    # iv_weak_threshold
    if (is.null(opt$estimator.args[["iv_weak"]])) {
      opt$estimator.args$iv_weak <- "warn"
    } else if (!tolower(opt$estimator.args[["iv_weak"]]) %in%
               c("warn", "prune", "none")) {
      lav_msg_stop(gettext("iv_weak should be warn, prune, or none."))
    }
    if (is.null(opt$estimator.args[["iv_weak_threshold"]])) {
      opt$estimator.args$iv_weak_threshold <- 10
    } else if (!is.numeric(opt$estimator.args[["iv_weak_threshold"]]) ||
               opt$estimator.args[["iv_weak_threshold"]] < 0) {
      lav_msg_stop(gettext(
        "iv_weak_threshold should be a non-negative number."))
    }
    if (is.null(opt$estimator.args$iv_vcov_stage1)) {
      if (opt$.categorical) {
        opt$estimator.args$iv_vcov_stage1 <- "gamma"
      } else {
        opt$estimator.args$iv_vcov_stage1 <- "lm.vcov.dfres"
      }
      if (tolower(opt$se[1]) == "none") {
        opt$estimator.args$iv_vcov_stage1 <- "none"
      }
    } else if (!tolower(opt$estimator.args$iv_vcov_stage1) %in%
               c("lm.vcov", "lm.vcov.dfres", "gamma", "none")) {
      lav_msg_stop(gettext("iv_vcov_stage1 should be lm.vcov, lm.vcov.dfres,
                            gamma or none."))
    } else if (tolower(opt$estimator.args$iv_vcov_stage1) == "gamma" &&
               !opt$estimator.args$iv_samplestats) {
      lav_msg_stop(gettext("iv_vcov_stage1 cannot be gamma if
                            iv_samplestats is FALSE"))
    }
    if (is.null(opt$estimator.args$iv_vcov_stage2)) {
      if (tolower(opt$se[1]) == "none") {
        opt$estimator.args$iv_vcov_stage2 <- "none"
      } else {
        opt$estimator.args$iv_vcov_stage2 <- "h2"
      }
    } else if (!tolower(opt$estimator.args$iv_vcov_stage2) %in%
               c("h2", "delta", "none")) {
      lav_msg_stop(gettext("iv_vcov_stage2 should be h2, delta, or none."))
    } else if (tolower(opt$estimator.args$iv_vcov_stage2) == "h2" &&
               !opt$estimator.args$iv_samplestats) {
      lav_msg_stop(gettext("iv_vcov_stage2 should be delta (or none) if
                   iv_samplestats = FALSE."))
    }
    if (opt$.categorical) {
      opt$estimator.args$iv_samplestats <- TRUE
      opt$estimator.args$iv_varcov_method = "ULS"
    }
    if (is.null(opt$estimator.args$iv_varcov_method)) {
      opt$estimator.args$iv_varcov_method <- "RLS"
    } else if (!toupper(opt$estimator.args$iv_varcov_method) %in%
               c("ULS", "GLS", "2RLS", "RLS", "NONE")) {
      lav_msg_stop(gettext("iv_varcov_method should be ULS, GLS, 2RLS, RLS
                            or NONE."))
    }
    if (is.null(opt$estimator.args$iv_vcov_gamma_modelbased)) {
      opt$estimator.args$iv_vcov_gamma_modelbased <- TRUE
    }
    if (is.null(opt$estimator.args$iv_mean_structure)) {
      opt$estimator.args$iv_mean_structure <- "wls"
    } else if (!tolower(opt$estimator.args$iv_mean_structure) %in%
               c("moments", "wls")) {
      lav_msg_stop(gettext("iv_mean_structure should be moments or wls."))
    }
    if (is.null(opt$estimator.args$iv_vcov_jack_numerical)) {
      opt$estimator.args$iv_vcov_jack_numerical <- FALSE
    }
    if (is.null(opt$estimator.args$iv_vcov_jaca_numerical)) {
      opt$estimator.args$iv_vcov_jaca_numerical <- FALSE
    }
    if (is.null(opt$estimator.args$iv_vcov_jacb_numerical)) {
      opt$estimator.args$iv_vcov_jacb_numerical <- FALSE
    }
  }

  # two-stage missing data needs the (EM) sample moments, and the standard
  # errors must run through the gamma path so the moment covariance can be
  # replaced by the two-stage one (see lav_sem_miiv_vcov())
  if (two_stage) {
    opt$estimator.args$iv_samplestats <- TRUE
    if (tolower(opt$se) != "none" &&
        tolower(opt$estimator.args$iv_vcov_stage1) %in%
          c("lm.vcov", "lm.vcov.dfres")) {
      opt$estimator.args$iv_vcov_stage1 <- "gamma"
    }
  }

  # override test if iv_varcov_method is "none"
  if (tolower(opt$estimator.args$iv_varcov_method) == "none") {
    opt$test <- "none"
    opt$standard.test <- "none"
  }

  opt
}

# normalize the rbm.method estimator.arg (short forms allowed)
# ("none" is an undocumented value that falls back to plain ML; kept for
#  debugging / parity checks)
lav_options_est_rbm_method <- function(x) {
  if (length(x) == 0L) {
    return("implicit")
  }
  x <- tolower(as.character(x)[1L])
  if (x %in% c("none", "ml", "false")) {
    "none"
  } else if (startsWith(x, "e")) {
    "explicit" # explicit, erbm, e
  } else if (startsWith(x, "i")) {
    "implicit" # implicit, irbm, i
  } else {
    lav_msg_stop(gettextf(
      "estimator.args$rbm.method must be one of %s.",
      lav_msg_view(c("implicit", "explicit"), log_sep = "or")))
  }
}

lav_options_est_rbm <- function(opt) {
  # RBM reduced-bias M-estimation (continuous SEM, for now)  # FIXME
  #
  # penalized-ML point estimation; the discrepancy, implied moments, loglik,
  # test and information matrices are all the ML ones, so we keep the estimator
  # as "ml" internally. The dedicated fit function is triggered in step 11 by
  # the presence of estimator.args$rbm.method (the optimizer itself is nlminb).

  # estimator.args
  if (is.null(opt$estimator.args)) {
    opt$estimator.args <- list(rbm.method = "implicit")
  } else {
    opt$estimator.args$rbm.method <-
      lav_options_est_rbm_method(opt$estimator.args$rbm.method)
  }

  # information: the penalty/bias need the observed information
  opt$information[1] <- "observed"
  if (length(opt$information) > 1L &&
      opt$information[2] == "default") {
    opt$information[2] <- "observed"
  }

  # se: sandwich by default
  if (opt$se == "default") {
    opt$se <- "robust.huber.white"
  } else if (opt$se == "robust") {
    opt$se <- "robust.huber.white"
  }

  # test
  if (opt$test[1] == "default") {
    opt$test <- "standard"
  }

  # missing
  if (opt$missing == "default") {
    opt$missing <- "listwise" # for now
  }

  # the actual optimizer is nlminb (the rbm fit is triggered in step 11 by
  # estimator.args$rbm.method)
  opt$optim.method <- "nlminb"

  # internally, treat as ML for all downstream machinery
  opt$estimator <- "ml"

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
