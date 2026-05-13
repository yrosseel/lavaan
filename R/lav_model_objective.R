# model objective

lav_model_objective <- function(lavmodel = NULL,
                                glist = NULL,
                                lavsamplestats = NULL,
                                lavdata = NULL,
                                lavcache = NULL) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  # shortcut for data.type == "none" or estimator == "none"
  if (lavmodel@estimator == "none" || length(lavsamplestats@cov) == 0L) {
    fx <- as.numeric(NA)
    attr(fx, "fx.group") <- rep(as.numeric(NA), lavsamplestats@ngroups)
    return(fx)
  }

  meanstructure <- lavmodel@meanstructure
  estimator <- lavmodel@estimator
  categorical <- lavmodel@categorical
  correlation <- lavmodel@correlation
  group_w_free <- lavmodel@group.w.free
  # fixed.x <- lavmodel@fixed.x
  conditional_x <- lavmodel@conditional.x
  num_idx <- lavmodel@num.idx
  th_idx <- lavmodel@th.idx
  estimator_args <- lavmodel@estimator.args

  # do we need WLS.est?
  if (estimator %in% c("ULS", "WLS", "DWLS", "NTRLS", "DLS")) {
    lavimplied <- lav_model_implied(lavmodel, GLIST = glist)
    # check for COV with negative diagonal elements?
    for (g in 1:lavsamplestats@ngroups) {
      cov_1 <- if (lavmodel@conditional.x) {
        lavimplied$res.cov[[g]]
      } else {
        lavimplied$cov[[g]]
      }
      d_cov <- diag(cov_1)
      if (anyNA(cov_1) || any(d_cov < 0)) {
        # return NA
        fx <- as.numeric(NA)
        attr(fx, "fx.group") <- rep(as.numeric(NA), lavsamplestats@ngroups)
        return(fx)
      }
    }
    wls_est <- lav_model_wls_est(
      lavmodel = lavmodel, glist = glist,
      lavimplied = lavimplied
    ) # ,
    # cov.x = lavsamplestats@cov.x)
    if (estimator == "NTRLS") {
      implied_fast <- lav_model_implied_fast(
        lavmodel = lavmodel, glist = glist,
        need_sigma = TRUE,
        need_mu = TRUE,
        extra = TRUE
      )
      sigma_hat <- implied_fast$sigma
      mu_hat <- implied_fast$mu
    }
    if (estimator == "DLS" && estimator_args$dls.GammaNT == "model") {
      implied_fast <- lav_model_implied_fast(
        lavmodel = lavmodel, glist = glist,
        need_sigma = TRUE,
        need_mu = TRUE,
        extra = FALSE
      )
      sigma_hat <- implied_fast$sigma
      mu_hat <- implied_fast$mu
    }
    if (lav_debug()) print(wls_est)
  } else if (estimator %in% c("ML", "GLS", "PML", "FML", "REML", "catML") &&
    lavdata@nlevels == 1L) {
    # compute moments for all groups
    # if(conditional.x) {
    #    Sigma.hat <- lav_model_cond2joint_sigma(lavmodel = lavmodel,
    #                     GLIST = GLIST, lavsamplestats = lavsamplestats,
    #                     extra = (estimator %in% c("ML", "REML","NTRLS")))
    # } else {
    implied_fast <- lav_model_implied_fast(
      lavmodel = lavmodel, glist = glist,
      need_sigma = TRUE,
      need_mu = meanstructure,
      need_th = categorical,
      need_pi = conditional_x,
      extra = (estimator %in% c(
        "ML", "REML",
        "NTRLS", "catML"
      ))
    )
    sigma_hat <- implied_fast$sigma
    # }

    # if (estimator == "REML") {
    #   LAMBDA <- lav_model_lambda(lavmodel = lavmodel, GLIST = GLIST)
    # }

    # ridge?
    if (lavsamplestats@ridge > 0.0) {
      for (g in 1:lavsamplestats@ngroups) {
        diag(sigma_hat[[g]]) <- diag(sigma_hat[[g]]) +
          lavsamplestats@ridge
      }
    }
    if (lav_debug()) print(sigma_hat)

    if (meanstructure) {
      # if(conditional.x) {
      #    Mu.hat <- lav_model_cond2joint_mu(lavmodel = lavmodel, GLIST = GLIST,
      #                           lavsamplestats = lavsamplestats)
      # } else {
      mu_hat <- implied_fast$mu
      # }
    }

    if (categorical) {
      th <- implied_fast$th
    }

    if (conditional_x) {
      pi0 <- implied_fast$pi
    }

    if (group_w_free) {
      gw <- lav_model_gw(lavmodel = lavmodel, glist = glist)
    }
  } else if (estimator == "MML") {
    th <- lav_model_th(lavmodel = lavmodel, glist = glist)
    mm_theta <- lav_model_theta(lavmodel = lavmodel, glist = glist)
    gw <- lav_model_gw(lavmodel = lavmodel, glist = glist)
  }

  fx <- 0.0
  fx_group <- numeric(lavsamplestats@ngroups)
  logl_group <- rep(as.numeric(NA), lavsamplestats@ngroups)

  for (g in 1:lavsamplestats@ngroups) {
    # incomplete data and fiml?
    if (lavsamplestats@missing.flag && estimator != "Bayes") {
      if (estimator == "ML" && lavdata@nlevels == 1L) {
        # FIML
        if (!attr(sigma_hat[[g]], "po")) {
          fx <- as.numeric(Inf)
          attr(fx, "fx.group") <- rep(as.numeric(Inf), lavsamplestats@ngroups)
          return(fx)
        }
        # check if h1 is defined (eg zero coverage)
        if (is.null(lavsamplestats@missing.h1[[g]]$h1)) {
          #this.h1 <- lav_mvnorm_missing_loglik_samplestats(
          #  Yp = lavsamplestats@missing[[g]],
          #  Mu = Mu.hat[[g]], Sigma = Sigma.hat[[g]],
          #  log2pi = FALSE, minus.two = TRUE) / lavsamplestats@nobs[[g]]
          #this.h1 <- this.h1 * 0.9999999 # avoid perfect fit
          this_h1 <- NULL #for now
        } else {
          this_h1 <- lavsamplestats@missing.h1[[g]]$h1
        }
        group_fx <- lav_model_objective_fiml(
          sigma_hat = sigma_hat[[g]],
          mu_hat = mu_hat[[g]],
          yp = lavsamplestats@missing[[g]],
          h1 = this_h1, n = lavsamplestats@nobs[[g]]
        )
      } else if (estimator == "ML" && lavdata@nlevels > 1L) {
        # FIML twolevel
        group_fx <- lav_model_objective_2l(
          lavmodel = lavmodel,
          glist = glist,
          y1 = lavdata@X[[g]],
          lp = lavdata@Lp[[g]],
          mp = lavdata@Mp[[g]],
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
        group_fx <- lav_model_objective_2l(
          lavmodel = lavmodel,
          glist = glist,
          lp = lavdata@Lp[[g]],
          mp = NULL, # complete data
          lavsamplestats = lavsamplestats,
          group = g
        )
      } else if (conditional_x) {
        group_fx <- lav_model_objective_ml_res(
          sigma_hat        = sigma_hat[[g]],
          mu_hat           = mu_hat[[g]],
          pi0              = pi0[[g]],
          res_cov          = lavsamplestats@res.cov[[g]],
          res_int          = lavsamplestats@res.int[[g]],
          res_slopes       = lavsamplestats@res.slopes[[g]],
          res_cov_log_det  = lavsamplestats@res.cov.log.det[[g]],
          cov_x            = lavsamplestats@cov.x[[g]],
          mean_x           = lavsamplestats@mean.x[[g]]
        )
      } else {
        group_fx <- lav_model_objective_ml(
          sigma_hat        = sigma_hat[[g]],
          mu_hat           = mu_hat[[g]],
          data_cov         = lavsamplestats@cov[[g]],
          data_mean        = lavsamplestats@mean[[g]],
          data_cov_log_det = lavsamplestats@cov.log.det[[g]],
          meanstructure    = meanstructure
        )
      }


      ### GLS #### (0.6-10: not using WLS function any longer)
    } else if (estimator == "GLS") {
      group_fx <- lav_model_objective_gls(
        sigma_hat        = sigma_hat[[g]],
        mu_hat           = mu_hat[[g]],
        data_cov         = lavsamplestats@cov[[g]],
        data_cov_inv     = lavsamplestats@icov[[g]],
        data_mean        = lavsamplestats@mean[[g]],
        meanstructure    = meanstructure,
    correlation      = correlation
      )
    } else if (estimator == "WLS" ||
      estimator == "DLS" ||
      estimator == "NTRLS") {
      # full weight matrix
      if (estimator == "WLS") {
        wls_v <- lavsamplestats@WLS.V[[g]]
      } else if (estimator == "DLS") {
        if (estimator_args$dls.GammaNT == "sample") {
          wls_v <- lavsamplestats@WLS.V[[g]]
        } else {
          dls_a <- estimator_args$dls.a
          gamma_nt <- lav_samplestats_gamma_nt(
            m_cov          = sigma_hat[[g]],
            m_mean         = mu_hat[[g]],
            x_idx          = lavsamplestats@x.idx[[g]],
            fixed_x        = lavmodel@fixed.x,
            conditional_x  = lavmodel@conditional.x,
            meanstructure  = lavmodel@meanstructure,
            slopestructure = lavmodel@conditional.x
          )
          w_dls <- (1 - dls_a) * lavsamplestats@NACOV[[g]] + dls_a * gamma_nt
          wls_v <- lav_matrix_symmetric_inverse(w_dls)
        }
      } else if (estimator == "NTRLS") {
        # WLS.V <- lav_samplestats_gamma_inverse_nt(
        #             m_icov = attr(Sigma.hat[[g]],"inv")[,,drop=FALSE],
        #             m_cov            = Sigma.hat[[g]][,,drop=FALSE],
        #             m_mean           = Mu.hat[[g]],
        #             x_idx          = c(10000,10001), ### FIXME!!!!
        #             fixed_x        = fixed.x,
        #             conditional_x  = conditional.x,
        #             meanstructure  = meanstructure,
        #             slopestructure = conditional.x)
        wls_v <- lav_mvnorm_information_expected(
          sigma_1 = sigma_hat[[g]],
          x_idx = lavsamplestats@x.idx[[g]],
          meanstructure = lavmodel@meanstructure
        )
        # DEBUG!!!!
        # WLS.V <- 2*WLS.V
      }

      group_fx <- lav_model_objective_wls(
        wls_est = wls_est[[g]],
        wls_obs = lavsamplestats@WLS.obs[[g]],
        wls_v = wls_v
      )
      attr(group_fx, "WLS.est") <- wls_est[[g]]
    } else if (estimator == "DWLS" || estimator == "ULS") {
      # diagonal weight matrix
      group_fx <- lav_model_objective_dwls(
        wls_est = wls_est[[g]],
        wls_obs = lavsamplestats@WLS.obs[[g]],
        wls_vd = lavsamplestats@WLS.VD[[g]]
      )
      attr(group_fx, "WLS.est") <- wls_est[[g]]
    } else if (estimator == "PML") {
      # Pairwise maximum likelihood
      if (lavdata@nlevels > 1L) {
        # group.fx <- lav_model_objective_pml.2L(lavmodel       = lavmodel,
        #                             GLIST          = GLIST,
        #                             Lp             = lavdata@Lp[[g]],
        #                             lavsamplestats = lavsamplestats,
        #                             group          = g)
        group_fx <- 0 # for now
        attr(group_fx, "logl") <- 0
      } else if (conditional_x) {
        group_fx <- lav_model_objective_pml(
          sigma_hat = sigma_hat[[g]],
          mu_hat = mu_hat[[g]],
          th = th[[g]],
          pi0 = pi0[[g]],
          th_idx = th_idx[[g]],
          num_idx = num_idx[[g]],
          x = lavdata@X[[g]],
          exo = lavdata@eXo[[g]],
          wt = lavdata@weights[[g]],
          lavcache = lavcache[[g]],
          missing = lavdata@missing
        )
      } else {
        group_fx <- lav_model_objective_pml(
          sigma_hat = sigma_hat[[g]],
          mu_hat = mu_hat[[g]],
          th = th[[g]],
          pi0 = NULL,
          th_idx = th_idx[[g]],
          num_idx = num_idx[[g]],
          x = lavdata@X[[g]],
          exo = NULL,
          wt = lavdata@weights[[g]],
          lavcache = lavcache[[g]],
          missing = lavdata@missing
        )
      }
      logl_group[g] <- attr(group_fx, "logl")
    } else if (estimator == "FML") {
      # Full maximum likelihood (underlying multivariate normal)
      group_fx <- lav_model_objective_fml(
        sigma_hat = sigma_hat[[g]],
        th = th[[g]],
        th_idx = th_idx[[g]],
        num_idx = num_idx[[g]],
        x = lavdata@X[[g]],
        lavcache = lavcache[[g]]
      )
    } else if (estimator == "MML") {
      # marginal maximum likelihood
      group_fx <- lav_model_objective_mml(
        lavmodel = lavmodel,
        glist = glist,
        mm_theta = mm_theta[[g]],
        th = th[[g]],
        group = g,
        lavdata = lavdata,
        sample_mean = lavsamplestats@mean[[g]],
        sample_mean_x = lavsamplestats@mean.x[[g]],
        lavcache = lavcache
      )
    } else if (estimator == "REML") {
      # restricted/residual maximum likelihood
      group_fx <- lav_model_objective_reml(
        sigma_hat = sigma_hat[[g]],
        mu_hat = mu_hat[[g]],
        data_cov = lavsamplestats@cov[[g]],
        data_mean = lavsamplestats@mean[[g]],
        data_cov_log_det = lavsamplestats@cov.log.det[[g]],
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
        group_fx <- 0.5 * group_fx ## FIXME
      }
    } else if (estimator == "PML" || estimator == "FML" ||
      estimator == "MML") {
      # do nothing
    } else if (estimator == "DLS") {
      if (estimator_args$dls.FtimesNminus1) {
        group_fx <- 0.5 * (lavsamplestats@nobs[[g]] - 1) /
            lavsamplestats@nobs[[g]] * group_fx
      } else {
        group_fx <- 0.5 * group_fx
      }
    } else {
      group_fx <- 0.5 * (lavsamplestats@nobs[[g]] - 1) /
          lavsamplestats@nobs[[g]] * group_fx
    }

    fx_group[g] <- group_fx
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
      fx <- sum(fx_group)
  } else if (lavdata@nlevels > 1L) {
    # no weighting needed! (implicit in obj, which is based on loglik)
    fx <- sum(fx_group)
    } else {
      nobs <- unlist(lavsamplestats@nobs)
      # }
      fx <- weighted.mean(fx_group, w = nobs)
    }
  } else { # single group
    fx <- fx_group[1]
  }

  # penalty for group.w + ML
  if (group_w_free && estimator %in% c(
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
    obs_freq <- unlist(lavsamplestats@group.w) * lavsamplestats@ntotal
    est_freq <- exp(unlist(gw))
    fx_w <- -1 * sum(obs_freq * log(est_freq) - est_freq)
    # divide by N (to be consistent with the rest of lavaan)
    fx_w <- fx_w / lavsamplestats@ntotal

    fx_sat <- sum(ifelse(obs_freq > 0, obs_freq * log(obs_freq), 0) - obs_freq)
    fx_sat <- fx_sat / lavsamplestats@ntotal

    # saturated - poisson
    # fx.w <- sum(obs.freq * log(obs.freq/est.freq))
    # does not work without constraints? --> need lagrange multiplier

    fx <- fx + (fx_w + fx_sat)
  }

  fx_value <- as.numeric(fx)

  attr(fx, "fx.group") <- fx_group
  if (estimator == "PML") {
    attr(fx, "logl.group") <- logl_group
    attr(fx, "fx.pml") <- fx_value
  }

  fx
}
