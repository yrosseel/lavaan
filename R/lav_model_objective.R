# model objective

lav_model_objective <- function(lavmodel = NULL,
                                GLIST = NULL,
                                lavsamplestats = NULL,
                                lavdata = NULL,
                                lavcache = NULL,
                                verbose = FALSE,
                                debug = FALSE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  # shortcut for data.type == "none" or estimator == "none"
  if (lavmodel@estimator == "none" || length(lavsamplestats@cov) == 0L) {
    fx <- as.numeric(NA)
    attr(fx, "fx.group") <- rep(as.numeric(NA), lavsamplestats@ngroups)
    return(fx)
  }

  meanstructure <- lavmodel@meanstructure
  estimator <- lavmodel@estimator
  categorical <- lavmodel@categorical
  if (.hasSlot(lavmodel, "correlation")) {
    correlation <- lavmodel@correlation
  } else {
    correlation <- FALSE
  }
  group.w.free <- lavmodel@group.w.free
  fixed.x <- lavmodel@fixed.x
  conditional.x <- lavmodel@conditional.x
  num.idx <- lavmodel@num.idx
  th.idx <- lavmodel@th.idx
  if (.hasSlot(lavmodel, "estimator.args")) {
    estimator.args <- lavmodel@estimator.args
  } else {
    estimator.args <- list()
  }


  # do we need WLS.est?
  if (estimator %in% c("ULS", "WLS", "DWLS", "NTRLS", "DLS")) {
    lavimplied <- lav_model_implied(lavmodel, GLIST = GLIST)
    # check for COV with negative diagonal elements?
    for (g in 1:lavsamplestats@ngroups) {
      COV <- if (lavmodel@conditional.x) {
        lavimplied$res.cov[[g]]
      } else {
        lavimplied$cov[[g]]
      }
      dCOV <- diag(COV)
      if (anyNA(COV) || any(dCOV < 0)) {
        # return NA
        fx <- as.numeric(NA)
        attr(fx, "fx.group") <- rep(as.numeric(NA), lavsamplestats@ngroups)
        return(fx)
      }
    }
    WLS.est <- lav_model_wls_est(
      lavmodel = lavmodel, GLIST = GLIST,
      lavimplied = lavimplied
    ) # ,
    # cov.x = lavsamplestats@cov.x)
    if (estimator == "NTRLS") {
      Sigma.hat <- computeSigmaHat(
        lavmodel = lavmodel, GLIST = GLIST,
        extra = TRUE
      )
      Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
    }
    if (estimator == "DLS" && estimator.args$dls.GammaNT == "model") {
      Sigma.hat <- computeSigmaHat(
        lavmodel = lavmodel, GLIST = GLIST,
        extra = FALSE
      )
      Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
    }
    if (debug) print(WLS.est)
  } else if (estimator %in% c("ML", "GLS", "PML", "FML", "REML", "catML") &&
    lavdata@nlevels == 1L) {
    # compute moments for all groups
    # if(conditional.x) {
    #    Sigma.hat <- computeSigmaHatJoint(lavmodel = lavmodel,
    #                     GLIST = GLIST, lavsamplestats = lavsamplestats,
    #                     extra = (estimator %in% c("ML", "REML","NTRLS")))
    # } else {
    Sigma.hat <- computeSigmaHat(
      lavmodel = lavmodel, GLIST = GLIST,
      extra = (estimator %in% c(
        "ML", "REML",
        "NTRLS", "catML"
      ))
    )
    # }

    if (estimator == "REML") {
      LAMBDA <- computeLAMBDA(lavmodel = lavmodel, GLIST = GLIST)
    }

    # ridge?
    if (lavsamplestats@ridge > 0.0) {
      for (g in 1:lavsamplestats@ngroups) {
        diag(Sigma.hat[[g]]) <- diag(Sigma.hat[[g]]) +
          lavsamplestats@ridge
      }
    }
    if (debug) print(Sigma.hat)

    if (meanstructure) {
      # if(conditional.x) {
      #    Mu.hat <- computeMuHatJoint(lavmodel = lavmodel, GLIST = GLIST,
      #                           lavsamplestats = lavsamplestats)
      # } else {
      Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
      # }
    }

    if (categorical) {
      TH <- computeTH(lavmodel = lavmodel, GLIST = GLIST)
    }

    if (conditional.x) {
      PI <- computePI(lavmodel = lavmodel, GLIST = GLIST)
    }

    if (group.w.free) {
      GW <- computeGW(lavmodel = lavmodel, GLIST = GLIST)
    }
  } else if (estimator == "MML") {
    TH <- computeTH(lavmodel = lavmodel, GLIST = GLIST)
    THETA <- computeTHETA(lavmodel = lavmodel, GLIST = GLIST)
    GW <- computeGW(lavmodel = lavmodel, GLIST = GLIST)
  }

  fx <- 0.0
  fx.group <- numeric(lavsamplestats@ngroups)
  logl.group <- rep(as.numeric(NA), lavsamplestats@ngroups)

  for (g in 1:lavsamplestats@ngroups) {
    # incomplete data and fiml?
    if (lavsamplestats@missing.flag && estimator != "Bayes") {
      if (estimator == "ML" && lavdata@nlevels == 1L) {
        # FIML
        if (!attr(Sigma.hat[[g]], "po")) {
          return(Inf)
        }
        group.fx <- estimator.FIML(
          Sigma.hat = Sigma.hat[[g]],
          Mu.hat = Mu.hat[[g]],
          Yp = lavsamplestats@missing[[g]],
          h1 = lavsamplestats@missing.h1[[g]]$h1, N = lavsamplestats@nobs[[g]]
        )
      } else if (estimator == "ML" && lavdata@nlevels > 1L) {
        # FIML twolevel
        group.fx <- estimator.2L(
          lavmodel = lavmodel,
          GLIST = GLIST,
          Y1 = lavdata@X[[g]],
          Lp = lavdata@Lp[[g]],
          Mp = lavdata@Mp[[g]],
          lavsamplestats = lavsamplestats,
          group = g
        )
      } else {
        lav_msg_stop(gettextf(
          "this estimator: `%s' can not be used with incomplete data and
          the missing=\"ml\" option", estimator))
      }
    } else if (estimator == "ML" || estimator == "Bayes" ||
      estimator == "catML") {
      # complete data
      # ML and friends
      if (lavdata@nlevels > 1L) {
        if (estimator %in% c("catML", "Bayes")) {
          lav_msg_stop(gettext("multilevel data not supported for estimator"),
                       estimator)
          }
        group.fx <- estimator.2L(
          lavmodel = lavmodel,
          GLIST = GLIST,
          Lp = lavdata@Lp[[g]],
          Mp = NULL, # complete data
          lavsamplestats = lavsamplestats,
          group = g
        )
      } else if (conditional.x) {
        group.fx <- estimator.ML_res(
          Sigma.hat        = Sigma.hat[[g]],
          Mu.hat           = Mu.hat[[g]],
          PI               = PI[[g]],
          res.cov          = lavsamplestats@res.cov[[g]],
          res.int          = lavsamplestats@res.int[[g]],
          res.slopes       = lavsamplestats@res.slopes[[g]],
          res.cov.log.det  = lavsamplestats@res.cov.log.det[[g]],
          cov.x            = lavsamplestats@cov.x[[g]],
          mean.x           = lavsamplestats@mean.x[[g]]
        )
      } else {
        group.fx <- estimator.ML(
          Sigma.hat        = Sigma.hat[[g]],
          Mu.hat           = Mu.hat[[g]],
          data.cov         = lavsamplestats@cov[[g]],
          data.mean        = lavsamplestats@mean[[g]],
          data.cov.log.det = lavsamplestats@cov.log.det[[g]],
          meanstructure    = meanstructure
        )
      }


      ### GLS #### (0.6-10: nog using WLS function any longer)
    } else if (estimator == "GLS") {
      group.fx <- estimator.GLS(
        Sigma.hat        = Sigma.hat[[g]],
        Mu.hat           = Mu.hat[[g]],
        data.cov         = lavsamplestats@cov[[g]],
        data.cov.inv     = lavsamplestats@icov[[g]],
        data.mean        = lavsamplestats@mean[[g]],
        meanstructure    = meanstructure
      )
    } else if (estimator == "WLS" ||
      estimator == "DLS" ||
      estimator == "NTRLS") {
      # full weight matrix
      if (estimator == "WLS") {
        WLS.V <- lavsamplestats@WLS.V[[g]]
      } else if (estimator == "DLS") {
        if (estimator.args$dls.GammaNT == "sample") {
          WLS.V <- lavsamplestats@WLS.V[[g]]
        } else {
          dls.a <- estimator.args$dls.a
          GammaNT <- lav_samplestats_Gamma_NT(
            COV            = Sigma.hat[[g]],
            MEAN           = Mu.hat[[g]],
            rescale        = FALSE,
            x.idx          = lavsamplestats@x.idx[[g]],
            fixed.x        = lavmodel@fixed.x,
            conditional.x  = lavmodel@conditional.x,
            meanstructure  = lavmodel@meanstructure,
            slopestructure = lavmodel@conditional.x
          )
          W.DLS <- (1 - dls.a) * lavsamplestats@NACOV[[g]] + dls.a * GammaNT
          WLS.V <- lav_matrix_symmetric_inverse(W.DLS)
        }
      } else if (estimator == "NTRLS") {
        # WLS.V <- lav_samplestats_Gamma_inverse_NT(
        #             ICOV = attr(Sigma.hat[[g]],"inv")[,,drop=FALSE],
        #             COV            = Sigma.hat[[g]][,,drop=FALSE],
        #             MEAN           = Mu.hat[[g]],
        #             x.idx          = c(10000,10001), ### FIXME!!!!
        #             fixed.x        = fixed.x,
        #             conditional.x  = conditional.x,
        #             meanstructure  = meanstructure,
        #             slopestructure = conditional.x)
        WLS.V <- lav_mvnorm_information_expected(
          Sigma = Sigma.hat[[g]],
          x.idx = lavsamplestats@x.idx[[g]],
          meanstructure = lavmodel@meanstructure
        )
        # DEBUG!!!!
        # WLS.V <- 2*WLS.V
      }

      group.fx <- estimator.WLS(
        WLS.est = WLS.est[[g]],
        WLS.obs = lavsamplestats@WLS.obs[[g]],
        WLS.V = WLS.V
      )
      attr(group.fx, "WLS.est") <- WLS.est[[g]]
    } else if (estimator == "DWLS" || estimator == "ULS") {
      # diagonal weight matrix
      group.fx <- estimator.DWLS(
        WLS.est = WLS.est[[g]],
        WLS.obs = lavsamplestats@WLS.obs[[g]],
        WLS.VD = lavsamplestats@WLS.VD[[g]]
      )
      attr(group.fx, "WLS.est") <- WLS.est[[g]]
    } else if (estimator == "PML") {
      # Pairwise maximum likelihood
      if (lavdata@nlevels > 1L) {
        # group.fx <- estimator.PML.2L(lavmodel       = lavmodel,
        #                             GLIST          = GLIST,
        #                             Lp             = lavdata@Lp[[g]],
        #                             lavsamplestats = lavsamplestats,
        #                             group          = g)
        group.fx <- 0 # for now
        attr(group.fx, "logl") <- 0
      } else if (conditional.x) {
        group.fx <- estimator.PML(
          Sigma.hat = Sigma.hat[[g]],
          Mu.hat = Mu.hat[[g]],
          TH = TH[[g]],
          PI = PI[[g]],
          th.idx = th.idx[[g]],
          num.idx = num.idx[[g]],
          X = lavdata@X[[g]],
          eXo = lavdata@eXo[[g]],
          wt = lavdata@weights[[g]],
          lavcache = lavcache[[g]],
          missing = lavdata@missing
        )
      } else {
        group.fx <- estimator.PML(
          Sigma.hat = Sigma.hat[[g]],
          Mu.hat = Mu.hat[[g]],
          TH = TH[[g]],
          PI = NULL,
          th.idx = th.idx[[g]],
          num.idx = num.idx[[g]],
          X = lavdata@X[[g]],
          eXo = NULL,
          wt = lavdata@weights[[g]],
          lavcache = lavcache[[g]],
          missing = lavdata@missing
        )
      }
      logl.group[g] <- attr(group.fx, "logl")
    } else if (estimator == "FML") {
      # Full maximum likelihood (underlying multivariate normal)
      group.fx <- estimator.FML(
        Sigma.hat = Sigma.hat[[g]],
        TH = TH[[g]],
        th.idx = th.idx[[g]],
        num.idx = num.idx[[g]],
        X = lavdata@X[[g]],
        lavcache = lavcache[[g]]
      )
    } else if (estimator == "MML") {
      # marginal maximum likelihood
      group.fx <- estimator.MML(
        lavmodel = lavmodel,
        GLIST = GLIST,
        THETA = THETA[[g]],
        TH = TH[[g]],
        group = g,
        lavdata = lavdata,
        sample.mean = lavsamplestats@mean[[g]],
        sample.mean.x = lavsamplestats@mean.x[[g]],
        lavcache = lavcache
      )
    } else if (estimator == "REML") {
      # restricted/residual maximum likelihood
      group.fx <- estimator.REML(
        Sigma.hat = Sigma.hat[[g]],
        Mu.hat = Mu.hat[[g]],
        data.cov = lavsamplestats@cov[[g]],
        data.mean = lavsamplestats@mean[[g]],
        data.cov.log.det = lavsamplestats@cov.log.det[[g]],
        meanstructure = meanstructure,
        group = g,
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata
      )
    } else {
      lav_msg_stop(gettext("unsupported estimator:"), estimator)
    }

    if (estimator %in% c("ML", "REML", "NTRLS", "catML")) {
      if (lavdata@nlevels == 1L) {
        group.fx <- 0.5 * group.fx ## FIXME
      }
    } else if (estimator == "PML" || estimator == "FML" ||
      estimator == "MML") {
      # do nothing
    } else if (estimator == "DLS") {
      if (estimator.args$dls.FtimesNminus1) {
        group.fx <- 0.5 * (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@nobs[[g]] * group.fx
      } else {
        group.fx <- 0.5 * group.fx
      }
    } else {
      group.fx <- 0.5 * (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@nobs[[g]] * group.fx
    }

    fx.group[g] <- group.fx
  } # g

  if (lavsamplestats@ngroups > 1) {
    ## FIXME: if group.w.free, should we use group.w or nobs???
    ##  - if we use estimated group.w, gradient changes!!!!
    ##  - but, if group models are misspecified, the group weights
    ##    will be affected too... which is unwanted (I think)
    # if(group.w.free) {
    # nobs <- unlist(GW) * lavsamplestats@ntotal
    #   nobs <- exp(unlist(GW))
    # } else {
    if (estimator == "PML") {
      # no weighting needed! (since N_g is part of the logl per group)
      fx <- sum(fx.group)
    } else {
      nobs <- unlist(lavsamplestats@nobs)
      # }
      fx <- weighted.mean(fx.group, w = nobs)
    }
  } else { # single group
    fx <- fx.group[1]
  }

  # penalty for group.w + ML
  if (group.w.free && estimator %in% c(
    "ML", "MML", "FML", "PML",
    "REML", "catML"
  )) {
    # obs.prop <- unlist(lavsamplestats@group.w)
    # est.prop <- unlist(GW)
    # if(estimator %in% c("WLS", "GLS", ...) {
    #    # X2 style discrepancy measures (aka GLS/WLS!!)
    #    fx.w <- sum ( (obs.prop-est.prop)^2/est.prop )
    # } else {
    #    # G2 style discrepancy measures (aka ML)
    #    # deriv is here -2 * (obs.prop - est.prop)
    # fx.w <- sum(obs.prop * log(obs.prop/est.prop) )
    # }

    # poisson kernel
    obs.freq <- unlist(lavsamplestats@group.w) * lavsamplestats@ntotal
    est.freq <- exp(unlist(GW))
    fx.w <- -1 * sum(obs.freq * log(est.freq) - est.freq)
    # divide by N (to be consistent with the rest of lavaan)
    fx.w <- fx.w / lavsamplestats@ntotal

    fx.sat <- sum(obs.freq * log(obs.freq) - obs.freq)
    fx.sat <- fx.sat / lavsamplestats@ntotal

    # saturated - poisson
    # fx.w <- sum(obs.freq * log(obs.freq/est.freq))
    # does not work without constraints? --> need lagrange multiplier

    fx <- fx + (fx.w + fx.sat)
  }

  fx.value <- as.numeric(fx)

  attr(fx, "fx.group") <- fx.group
  if (estimator == "PML") {
    attr(fx, "logl.group") <- logl.group
    attr(fx, "fx.pml") <- fx.value
  }

  fx
}
