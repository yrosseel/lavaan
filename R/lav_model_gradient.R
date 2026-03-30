# model gradient

lav_model_gradient <- function(lavmodel = NULL,
                               GLIST = NULL,
                               lavsamplestats = NULL,
                               lavdata = NULL,
                               lavcache = NULL,
                               type = "free",
                               group.weight = TRUE,
                               Delta = NULL,
                               m.el.idx = NULL,
                               x.el.idx = NULL,
                               ceq.simple = FALSE) {
  nmat <- lavmodel@nmat
  estimator <- lavmodel@estimator
  representation <- lavmodel@representation
  meanstructure <- lavmodel@meanstructure
  categorical <- lavmodel@categorical
  group.w.free <- lavmodel@group.w.free
  fixed.x <- lavmodel@fixed.x
  conditional.x <- lavmodel@conditional.x
  num.idx <- lavmodel@num.idx
  th.idx <- lavmodel@th.idx
  nx.free <- lavmodel@nx.free
  estimator.args <- lavmodel@estimator.args

  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  if (estimator == "REML") lav_msg_warn(gettext(
    "analytical gradient not implement; use numerical approximation"))

  # group.weight
  # FIXME --> block.weight
  if (group.weight) {
    if (estimator %in% c("ML", "PML", "FML", "MML", "REML", "NTRLS", "catML")) {
      group.w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)
    } else if (estimator == "DLS") {
      if (estimator.args$dls.FtimesNminus1) {
        group.w <- ((unlist(lavsamplestats@nobs) - 1) / lavsamplestats@ntotal)
      } else {
        group.w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)
      }
    } else {
      # FIXME: double check!
      group.w <- ((unlist(lavsamplestats@nobs) - 1) / lavsamplestats@ntotal)
    }
  } else {
    group.w <- rep(1.0, lavmodel@nblocks)
  }

  # do we need WLS.est?
  if (estimator %in% c("WLS", "DWLS", "ULS", "GLS", "NTRLS", "DLS")) {
    # always compute WLS.est
    WLS.est <- lav_model_wls_est(lavmodel = lavmodel, GLIST = GLIST) # ,
    # cov.x = lavsamplestats@cov.x)
  }

  if (estimator %in% c("ML", "PML", "FML", "REML", "NTRLS", "catML")) {
    # compute moments for all groups
    # if(conditional.x) {
    #    Sigma.hat <- lav_model_cond2joint_sigma(lavmodel = lavmodel,
    #                     GLIST = GLIST,
    #                     extra = (estimator %in% c("ML", "REML","NTRLS")))
    # } else {
    Sigma.hat <- lav_model_sigma(
      lavmodel = lavmodel, GLIST = GLIST,
      extra = (estimator %in% c(
        "ML", "REML",
        "NTRLS", "catML"
      ))
    )
    # }

    if (meanstructure) {
      # if(conditional.x) {
      #    Mu.hat <- lav_model_mu(lavmodel = lavmodel, GLIST = GLIST)
      # } else {
      Mu.hat <- lav_model_mu(lavmodel = lavmodel, GLIST = GLIST)
      # }
    }

    if (categorical) {
      TH <- lav_model_th(lavmodel = lavmodel, GLIST = GLIST)
    }

    if (conditional.x) {
      PI <- lav_model_pi(lavmodel = lavmodel, GLIST = GLIST)
    } else if (estimator == "PML") {
      PI <- vector("list", length = lavmodel@nblocks)
    }

    if (group.w.free) {
      GW <- lav_model_gw(lavmodel = lavmodel, GLIST = GLIST)
    }
  } else if (estimator == "DLS" && estimator.args$dls.GammaNT == "model") {
    Sigma.hat <- lav_model_sigma(
      lavmodel = lavmodel, GLIST = GLIST,
      extra = FALSE
    )
    Mu.hat <- lav_model_mu(lavmodel = lavmodel, GLIST = GLIST)
  } else if (estimator == "MML") {
    TH <- lav_model_th(lavmodel = lavmodel, GLIST = GLIST)
    THETA <- lav_model_theta(lavmodel = lavmodel, GLIST = GLIST)
    GW <- lav_model_gw(lavmodel = lavmodel, GLIST = GLIST)
  }

  # four approaches (FIXME!!!! merge this!)
  # - ML approach: using Omega (and Omega.mu)
  #    Omega = 'POST' = Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)
  #   (still 2x faster than Delta method)
  # - WLS/DWLS/GLS: using Delta + WLS.V; support for fixed.x, conditional.x
  # - (ML)/NTRLS: using Delta, no support for fixed.x, conditional.x
  # - PML/FML/MML: custom

  # composites?
  composites.flag <- lavmodel@composites

  # 1. ML approach
  if ((estimator == "ML" || estimator == "REML" || estimator == "catML") &&
    lavdata@nlevels == 1L && !composites.flag &&
    !lavmodel@conditional.x) {
	correlation <- lavmodel@correlation
    if (meanstructure) {
      Omega <- lav_model_omega(
        Sigma.hat = Sigma.hat, Mu.hat = Mu.hat,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = TRUE,
        conditional.x = conditional.x,
		correlation = correlation
      )
      Omega.mu <- attr(Omega, "mu")
    } else {
      Omega <- lav_model_omega(
        Sigma.hat = Sigma.hat, Mu.hat = NULL,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = FALSE,
        conditional.x = conditional.x,
		correlation = correlation
      )
      Omega.mu <- vector("list", length = lavmodel@nblocks)
    }

    # compute DX (for all elements in every model matrix)
    DX <- vector("list", length = length(GLIST))

    for (g in 1:lavmodel@nblocks) {
      # which mm belong to group g?
      mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
      mm.names <- names(GLIST[mm.in.group])

      if (representation == "LISREL") {
        DX.group <- lav_lisrel_df_dmlist(
          GLIST[mm.in.group],
          Omega[[g]],
          Omega.mu[[g]]
        )

        # FIXME!!!
        # add empty gamma
        if (lavmodel@conditional.x) {
          DX.group$gamma <- lavmodel@GLIST$gamma
        }

        # only save what we need
        DX[mm.in.group] <- DX.group[mm.names]
      } else if (representation == "RAM") {
        DX.group <- lav_ram_df(
          GLIST[mm.in.group],
          Omega[[g]],
          Omega.mu[[g]]
        )
        # only save what we need
        DX[mm.in.group] <- DX.group[mm.names]
      } else {
        lav_msg_stop(gettext(
          "only LISREL and RAM representation has been implemented for now"))
      }

      # weight by group
      if (lavmodel@nblocks > 1L) {
        for (mm in mm.in.group) {
          DX[[mm]] <- group.w[g] * DX[[mm]]
        }
      }
    }

    # extract free parameters

    if (type == "free") {
      if (lavmodel@ceq.simple.only) { # new in 0.6-11
        dx <- numeric(lavmodel@nx.unco)
        for (g in 1:lavmodel@nblocks) {
          mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
          for (mm in mm.in.group) {
            m.free.idx <- lavmodel@m.free.idx[[mm]]
            x.unco.idx <- lavmodel@x.unco.idx[[mm]]
            dx[x.unco.idx] <- DX[[mm]][m.free.idx]
          }
        }
        if (ceq.simple) {
          dx <- drop(crossprod(lavmodel@ceq.simple.K, dx))
        }
      } else {
        dx <- numeric(nx.free)
        for (g in 1:lavmodel@nblocks) {
          mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
          for (mm in mm.in.group) {
            m.free.idx <- lavmodel@m.free.idx[[mm]]
            x.free.idx <- lavmodel@x.free.idx[[mm]]
            dx[x.free.idx] <- DX[[mm]][m.free.idx]
          }
        }
      }
    } else {
      dx <- DX
      # handle equality constraints
      ### FIXME!!!! TODO!!!!
    }
  } else # ML

  # 2. using Delta - *LS family
  if (estimator %in% c("WLS", "DWLS", "ULS", "GLS", "NTGLS", "DLS")) {
    if (type != "free") {
      if (is.null(Delta)) {
        lav_msg_fixme("Delta should be given if type != free")
      }
      # stop("FIXME: WLS gradient with type != free needs fixing!")
    } else {
      Delta <- lav_model_delta(
        lavmodel = lavmodel, GLIST. = GLIST,
        ceq.simple = ceq.simple
      )
    }

    for (g in 1:lavmodel@nblocks) {
      # diff <- as.matrix(lavsamplestats@WLS.obs[[g]]  - WLS.est[[g]])
      # group.dx <- -1 * ( t(Delta[[g]]) %*% lavsamplestats@WLS.V[[g]] %*% diff)
      # 0.5-17: use crossprod twice; treat DWLS/ULS special
      if (estimator == "WLS" ||
        estimator == "GLS" ||
        estimator == "DLS" ||
        estimator == "NTRLS") {
        # full weight matrix
        diff <- lavsamplestats@WLS.obs[[g]] - WLS.est[[g]]

        # full weight matrix
        if (estimator == "GLS" || estimator == "WLS") {
          WLS.V <- lavsamplestats@WLS.V[[g]]
          group.dx <- -1 * crossprod(
            Delta[[g]],
            crossprod(WLS.V, diff)
          )
        } else if (estimator == "DLS") {
          if (estimator.args$dls.GammaNT == "sample") {
            WLS.V <- lavsamplestats@WLS.V[[g]] # for now
          } else {
            dls.a <- estimator.args$dls.a
            GammaNT <- lav_samplestats_gamma_nt(
              m_cov          = Sigma.hat[[g]],
              m_mean         = Mu.hat[[g]],
              rescale        = FALSE,
              x_idx          = lavsamplestats@x.idx[[g]],
              fixed_x        = lavmodel@fixed.x,
              conditional_x  = lavmodel@conditional.x,
              meanstructure  = lavmodel@meanstructure,
              slopestructure = lavmodel@conditional.x
            )
            W.DLS <- (1 - dls.a) * lavsamplestats@NACOV[[g]] + dls.a * GammaNT
            WLS.V <- lav_matrix_symmetric_inverse(W.DLS)
          }
          group.dx <- -1 * crossprod(
            Delta[[g]],
            crossprod(WLS.V, diff)
          )
        } else if (estimator == "NTRLS") {
          stopifnot(!conditional.x)
          # WLS.V <- lav_samplestats_gamma_inverse_nt(
          #         m_icov = attr(Sigma.hat[[g]],"inv")[,,drop=FALSE],
          #         m_cov            = Sigma.hat[[g]][,,drop=FALSE],
          #         m_mean           = Mu.hat[[g]],
          #         x_idx          = lavsamplestats@x.idx[[g]],
          #         fixed_x        = fixed.x,
          #         conditional_x  = conditional.x,
          #         meanstructure  = meanstructure,
          #         slopestructure = conditional.x)

          S <- lavsamplestats@cov[[g]]
          Sigma <- Sigma.hat[[g]]
          Sigma.inv <- attr(Sigma, "inv")
          nvar <- NROW(Sigma)

          if (meanstructure) {
            MEAN <- lavsamplestats@mean[[g]]
            Mu <- Mu.hat[[g]]
            POST.Sigma <- lav_matrix_duplication_pre(
              matrix(
                (Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)) %*%
                  (diag(nvar) + (S - Sigma) %*% Sigma.inv) +
                  (Sigma.inv %*% tcrossprod(MEAN - Mu) %*% Sigma.inv),
                ncol = 1
              )
            )
            POST.Mu <- as.numeric(2 * Sigma.inv %*% (MEAN - Mu))
            POST <- c(POST.Mu, POST.Sigma)
          } else {
            POST <- lav_matrix_duplication_pre(
              matrix((Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)) %*%
                (diag(nvar) + (S - Sigma) %*% Sigma.inv), ncol = 1)
            )
          }

          group.dx <- as.numeric(-1 * crossprod(Delta[[g]], POST))
        }
      } else if (estimator == "DWLS" || estimator == "ULS") {
        # diagonal weight matrix
        diff <- lavsamplestats@WLS.obs[[g]] - WLS.est[[g]]
        group.dx <- -1 * crossprod(
          Delta[[g]],
          lavsamplestats@WLS.VD[[g]] * diff
        )
      }

      group.dx <- group.w[g] * group.dx
      if (g == 1) {
        dx <- group.dx
      } else {
        dx <- dx + group.dx
      }
    } # g

    if (type == "free") {
      # nothing to do
    } else {
      # make a GLIST
      dx <- lav_model_x2glist(
        lavmodel = lavmodel, x = dx,
        type = "custom", setDelta = FALSE,
        m.el.idx = m.el.idx,
        x.el.idx = x.el.idx
      )
    }
  } # WLS

  # ML + conditional.x
  else if (estimator %in% c("ML", "catML") && lavmodel@conditional.x &&
    lavdata@nlevels == 1L) {
    if (type != "free") {
      if (is.null(Delta)) {
        lav_msg_fixme("Delta should be given if type != free")
      }
      # stop("FIXME: WLS gradient with type != free needs fixing!")
    } else {
      Delta <- lav_model_delta(
        lavmodel = lavmodel, GLIST. = GLIST,
        ceq.simple = ceq.simple
      )
    }

    for (g in 1:lavmodel@nblocks) {
      # augmented mean.x + cov.x matrix
      mean.x <- lavsamplestats@mean.x[[g]]
      cov.x <- lavsamplestats@cov.x[[g]]
      C3 <- rbind(
        c(1, mean.x),
        cbind(mean.x, cov.x + tcrossprod(mean.x))
      )

      Sigma <- Sigma.hat[[g]]
      Mu.g <- Mu.hat[[g]]
      PI.g <- PI[[g]]
      Sigma.inv <- attr(Sigma, "inv")
      nvar <- NROW(Sigma)
      S <- lavsamplestats@res.cov[[g]]

      # beta
      OBS <- t(cbind(
        lavsamplestats@res.int[[g]],
        lavsamplestats@res.slopes[[g]]
      ))
      EST <- t(cbind(Mu.g, PI.g))
      # obs.beta <- c(lavsamplestats@res.int[[g]],
      #              lav_matrix_vec(lavsamplestats@res.slopes[[g]]))
      # est.beta <- c(Mu.g,  lav_matrix_vec(PI.g))
      # beta.COV <- C3 %x% Sigma.inv

      # a <- t(obs.beta - est.beta)
      # b <- as.matrix(obs.beta - est.beta)
      # K <- lav_matrix_commutation(m = nvar, n = nvar)
      # AB <- (K %x% diag(NROW(C3)*NROW(C3))) %*%
      #          (diag(nvar) %x% lav_matrix_vec(C3) %x% diag(nvar))
      # K <- lav_matrix_commutation(m = nvar, n = NROW(C3))
      # AB <- ( diag(NROW(C3)) %x% K %x% diag(nvar) ) %*%
      #        (lav_matrix_vec(C3) %x% diag( nvar * nvar) )

      # POST.beta <- 2 *  beta.COV %*% (obs.beta - est.beta)
      d.BETA <- C3 %*% (OBS - EST) %*% Sigma.inv
      # NOTE: the vecr here, unlike lav_mvreg_dlogl_beta
      #       this is because DELTA has used vec(t(BETA)),
      #       instead of vec(BETA)
      # POST.beta <- 2 * lav_matrix_vecr(d.BETA)
      # NOT any longer, since 0.6-1!!!
      POST.beta <- 2 * lav_matrix_vec(d.BETA)

      # POST.sigma1 <- lav_matrix_duplication_pre(
      #        (Sigma.inv %x% Sigma.inv)  %*% t(AB)  %*% (t(a) %x% b) )

      # Sigma
      # POST.sigma2 <- lav_matrix_duplication_pre(
      #                 matrix( lav_matrix_vec(
      #          Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)), ncol = 1L))
      W.tilde <- S + t(OBS - EST) %*% C3 %*% (OBS - EST)
      d.SIGMA <- (Sigma.inv - Sigma.inv %*% W.tilde %*% Sigma.inv)
      d.vechSigma <- as.numeric(lav_matrix_duplication_pre(
        as.matrix(lav_matrix_vec(d.SIGMA))
      ))
      POST.sigma <- -1 * d.vechSigma

      # POST <- c(POST.beta, POST.sigma1 + POST.sigma2)
      POST <- c(POST.beta, POST.sigma)

      group.dx <- as.numeric(-1 * crossprod(Delta[[g]], POST))

      # because we still use obj/2, we need to divide by 2!
      group.dx <- group.dx / 2 # fixed in 0.6-1

      group.dx <- group.w[g] * group.dx
      if (g == 1) {
        dx <- group.dx
      } else {
        dx <- dx + group.dx
      }
    } # g

    if (type == "free") {
      # nothing to do
    } else {
      # make a GLIST
      dx <- lav_model_x2glist(
        lavmodel = lavmodel, x = dx,
        type = "custom", setDelta = FALSE,
        m.el.idx = m.el.idx,
        x.el.idx = x.el.idx
      )
    }
  } # ML + conditional.x

  else if (estimator == "ML" && lavdata@nlevels > 1L) {
    if (type != "free") {
      lav_msg_fixme("type != free in lav_model_gradient for
                    estimator ML for nlevels > 1")
    } else {
      Delta <- lav_model_delta(
        lavmodel = lavmodel, GLIST. = GLIST,
        ceq.simple = ceq.simple
      )
    }

    # for each upper-level group....
    for (g in 1:lavmodel@ngroups) {
      if (!lavsamplestats@missing.flag) { # complete data
        if (lavmodel@conditional.x) {
          DX <- lav_mvreg_cluster_dlogl_2l_samplestats(
            YLp = lavsamplestats@YLp[[g]],
            Lp = lavdata@Lp[[g]],
            Res.Sigma.W = Sigma.hat[[(g - 1) * 2 + 1]],
            Res.Int.W = Mu.hat[[(g - 1) * 2 + 1]],
            Res.Pi.W = PI[[(g - 1) * 2 + 1]],
            Res.Sigma.B = Sigma.hat[[(g - 1) * 2 + 2]],
            Res.Int.B = Mu.hat[[(g - 1) * 2 + 2]],
            Res.Pi.B = PI[[(g - 1) * 2 + 2]],
            Sinv.method = "eigen"
          )
        } else {
          DX <- lav_mvnorm_cluster_dlogl_2l_samplestats(
            YLp = lavsamplestats@YLp[[g]],
            Lp = lavdata@Lp[[g]],
            Mu.W = Mu.hat[[(g - 1) * 2 + 1]],
            Sigma.W = Sigma.hat[[(g - 1) * 2 + 1]],
            Mu.B = Mu.hat[[(g - 1) * 2 + 2]],
            Sigma.B = Sigma.hat[[(g - 1) * 2 + 2]],
            Sinv.method = "eigen"
          )
        }
      } else {
        # missing data
        if (lavmodel@conditional.x) {
          lav_msg_stop(gettext("gradient for twolevel + conditional.x + fiml
                  is not ready; use optim.gradient = \"numerical\""))
        } else {
          DX <- lav_mvnorm_cluster_missing_dlogl_2l_samplestats(
            Y1 = lavdata@X[[g]],
            Y2 = lavsamplestats@YLp[[g]][[2]]$Y2,
            Lp = lavdata@Lp[[g]],
            Mp = lavdata@Mp[[g]],
            Mu.W = Mu.hat[[(g - 1) * 2 + 1]],
            Sigma.W = Sigma.hat[[(g - 1) * 2 + 1]],
            Mu.B = Mu.hat[[(g - 1) * 2 + 2]],
            Sigma.B = Sigma.hat[[(g - 1) * 2 + 2]],
            Sinv.method = "eigen"
          )
        }
      }

      group.dx <- as.numeric(DX %*% Delta[[g]])

      # group weights (if any)
      group.dx <- group.w[g] * group.dx
      if (g == 1) {
        dx <- group.dx
      } else {
        dx <- dx + group.dx
      }
    } # g

    # divide by 2 * N
    dx <- dx / (2 * lavsamplestats@ntotal)

    # cat("dx1 (numerical) = \n"); print( zapsmall(dx1) )
    # cat("dx  (analytic)  = \n"); print( zapsmall(dx ) )
  } # ML + two-level

  else if (estimator == "PML" || estimator == "FML" ||
    estimator == "MML") {
    if (type != "free") {
      lav_msg_fixme("type != free in lav_model_gradient for estimator PML")
    } else {
      Delta <- lav_model_delta(
        lavmodel = lavmodel, GLIST. = GLIST,
        ceq.simple = ceq.simple
      )
    }

    for (g in 1:lavmodel@nblocks) {
      # print(GLIST)
      # print(lav_model_get_parameters(lavmodel = lavmodel, GLIST = GLIST))
      # print(Sigma.hat[[g]])
      # print(TH[[g]])
      # cat("*****\n")

      # compute partial derivative of logLik with respect to
      # thresholds/means, slopes, variances, correlations
      if (estimator == "PML") {
        if (lavdata@nlevels > 1L) {
          lav_msg_stop(gettext(
            "PL gradient + multilevel not implemented;
            try optim.gradient = \"numerical\""))
        } else if (conditional.x) {
          d1 <- lav_pml_dploglik_dimplied(
            Sigma.hat = Sigma.hat[[g]],
            Mu.hat = Mu.hat[[g]],
            TH = TH[[g]],
            th.idx = th.idx[[g]],
            num.idx = num.idx[[g]],
            X = lavdata@X[[g]],
            lavcache = lavcache[[g]],
            eXo = lavdata@eXo[[g]],
            wt = lavdata@weights[[g]],
            PI = PI[[g]],
            missing = lavdata@missing
          )
        } else {
          d1 <- lav_pml_dploglik_dimplied(
            Sigma.hat = Sigma.hat[[g]],
            Mu.hat = Mu.hat[[g]],
            TH = TH[[g]],
            th.idx = th.idx[[g]],
            num.idx = num.idx[[g]],
            X = lavdata@X[[g]],
            lavcache = lavcache[[g]],
            eXo = NULL,
            wt = lavdata@weights[[g]],
            PI = NULL,
            missing = lavdata@missing
          )
        } # not conditional.x

        # chain rule (fmin)
        group.dx <-
          as.numeric(t(d1) %*% Delta[[g]])
      } # PML

      else if (estimator == "FML") {
        d1 <- lav_pml_fml_dploglik_dimplied(
          Sigma.hat = Sigma.hat[[g]],
          TH = TH[[g]],
          th.idx = th.idx[[g]],
          num.idx = num.idx[[g]],
          X = lavdata@X[[g]],
          lavcache = lavcache[[g]]
        )

        # chain rule (fmin)
        group.dx <-
          as.numeric(t(d1) %*% Delta[[g]]) / lavsamplestats@nobs[[g]]
      } else if (estimator == "MML") {
        group.dx <-
          lav_model_gradient_mml(
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
      }

      # group weights (if any)
      group.dx <- group.w[g] * group.dx
      if (g == 1) {
        dx <- group.dx
      } else {
        dx <- dx + group.dx
      }
    } # g
  } else {
    lav_msg_stop(gettext(
      "no analytical gradient available for estimator"), estimator)
  }


  # group.w.free for ML
  if (lavmodel@group.w.free &&
    estimator %in% c("ML", "MML", "FML", "PML", "REML", "catML")) {
    # est.prop <- unlist( lav_model_gw(lavmodel = lavmodel, GLIST = GLIST) )
    # obs.prop <- unlist(lavsamplestats@group.w)
    # FIXME: G2 based -- ML and friends only!!
    # dx.GW <- - (obs.prop - est.prop)

    # poisson version
    est.freq <- exp(unlist(lav_model_gw(lavmodel = lavmodel, GLIST = GLIST)))
    obs.freq <- unlist(lavsamplestats@group.w) * lavsamplestats@ntotal
    dx.GW <- -(obs.freq - est.freq)
    # divide by N (to be consistent with the rest of lavaan)
    dx.GW <- dx.GW / lavsamplestats@ntotal

    # remove last element (fixed LAST group to zero)
    # dx.GW <- dx.GW[-length(dx.GW)]

    # fill in in dx
    gw.mat.idx <- which(names(lavmodel@GLIST) == "gw")
    gw.x.idx <- unlist(lavmodel@x.free.idx[gw.mat.idx])
    dx[gw.x.idx] <- dx.GW
  }

  # dx is 1xnpar matrix of LIST (type != "free")
  if (is.matrix(dx)) {
    dx <- as.numeric(dx)
  }

  dx
}

# for testing purposes only
lav_model_delta_numerical <- function(lavmodel = NULL, GLIST = NULL, g = 1L) {

  # state or final?
  if(is.null(GLIST)) GLIST <- lavmodel@GLIST

  compute.moments <- function(x) {
    GLIST <- lav_model_x2glist(lavmodel = lavmodel, x = x, type = "free")
    Sigma.hat <- lav_model_sigma(lavmodel = lavmodel, GLIST = GLIST)
    S.vec <- lav_matrix_vech(Sigma.hat[[g]])
    if(lavmodel@meanstructure) {
      Mu.hat <- lav_model_mu(lavmodel = lavmodel, GLIST=GLIST)
      out <- c(Mu.hat[[g]], S.vec)
    } else {
      out <- S.vec
    }
    out
  }

  x <- lav_model_get_parameters(lavmodel = lavmodel,
                                GLIST = GLIST, type = "free")
  Delta <- lav_func_jacobian_complex(func = compute.moments, x = x)

  Delta
}


### FIXME: should we here also:
###        - weight for groups? (no, for now)
###        - handle equality constraints? (yes, for now)
lav_model_delta <- function(lavmodel = NULL, GLIST. = NULL,
                         m.el.idx. = NULL, x.el.idx. = NULL,
                         ceq.simple = FALSE,
                         force.conditional.x.false = FALSE) {

  # temporary fix
  if (lavmodel@composites) {
    Delta <- vector("list", length = lavmodel@nblocks)
    for (b in seq_len(lavmodel@nblocks)) {
      Delta[[b]] <- lav_model_delta_numerical(lavmodel = lavmodel, g = b)  
    }  
    return(Delta)
  }

  representation <- lavmodel@representation
  categorical <- lavmodel@categorical
  correlation <- lavmodel@correlation
  conditional.x <- lavmodel@conditional.x
  group.w.free <- lavmodel@group.w.free
  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  nvar <- lavmodel@nvar
  num.idx <- lavmodel@num.idx
  th.idx <- lavmodel@th.idx
  nexo <- lavmodel@nexo
  parameterization <- lavmodel@parameterization

  # number of thresholds per group (if any)
  nth <- sapply(th.idx, function(x) sum(x > 0L))

  # state or final?
  if (is.null(GLIST.)) {
    GLIST <- lavmodel@GLIST
  } else {
    GLIST <- GLIST.
  }

  # type = "free" or something else?
  type <- "nonfree"
  m.el.idx <- m.el.idx.
  x.el.idx <- x.el.idx.
  if (is.null(m.el.idx) && is.null(x.el.idx)) {
    type <- "free"
  }

  # number of rows in DELTA.group
  pstar <- integer(nblocks)
  for (g in 1:nblocks) {
    pstar[g] <- as.integer(nvar[g] * (nvar[g] + 1) / 2)
    if (lavmodel@meanstructure) {
      pstar[g] <- nvar[g] + pstar[g] # first the means, then sigma
    }
    if (categorical) {
      pstar[g] <- pstar[g] - nvar[g] # remove variances
      pstar[g] <- pstar[g] - nvar[g] # remove means

      pstar[g] <- pstar[g] + nth[g] # add thresholds
      pstar[g] <- pstar[g] + length(num.idx[[g]]) # add num means
      pstar[g] <- pstar[g] + length(num.idx[[g]]) # add num vars
    } else if (correlation) {
      pstar[g] <- pstar[g] - nvar[g] # remove variances
    }
    if (conditional.x && nexo[g] > 0L) {
      pstar[g] <- pstar[g] + (nvar[g] * nexo[g]) # add slopes
    }
    if (group.w.free) {
      pstar[g] <- pstar[g] + 1L # add group weight
    }
  }


  # number of columns in DELTA + m.el.idx/x.el.idx
  if (type == "free") {
    if (lavmodel@ceq.simple.only) {
      NCOL <- lavmodel@nx.unco
    } else {
      NCOL <- lavmodel@nx.free
    }
    m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
    for (mm in 1:length(GLIST)) {
      m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
      if (lavmodel@ceq.simple.only) {
        x.el.idx[[mm]] <- lavmodel@x.unco.idx[[mm]]
      } else {
        x.el.idx[[mm]] <- lavmodel@x.free.idx[[mm]]
      }
      # handle symmetric matrices
      if (lavmodel@isSymmetric[mm]) {
        # since we use 'x.free.idx', only symmetric elements
        # are duplicated (not the equal ones, only in x.free.free)
        dix <- duplicated(x.el.idx[[mm]])
        if (any(dix)) {
          m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
          x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
        }
      }
    }
  } else {
    ## FIXME: this does *not* take into account symmetric
    ##        matrices; hence NCOL will be too large, and empty
    ##        columns will be added
    ##        this is ugly, but it doesn't hurt
    ## alternative could be:
    ## NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
    # NCOL <- sum(unlist(lapply(m.el.idx, length)))
    NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
    # sanity check
    # nx <- sum(unlist(lapply(x.el.idx, length)))
    # stopifnot(NCOL == nx)
  }


  # compute Delta
  Delta <- vector("list", length = nblocks)
  for (g in 1:nblocks) {
    Delta.group <- matrix(0, nrow = pstar[g], ncol = NCOL)

    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]

    # label rows of Delta.group --- FIXME!!!
    # if(categorical) {
    #    # 1. th (means interleaved?)
    #    # 2. pi
    #    # 3. var num + cor
    # } else {
    #    if(meanstructure) {
    #    }
    # }
    # if(group.w.free) {
    # }

    # if theta, do some preparation
    if (representation == "LISREL" && parameterization == "theta") {
      sigma.hat <- lav_lisrel_sigma(
        MLIST = GLIST[mm.in.group],
        delta = FALSE
      )
      dsigma <- diag(sigma.hat)
      # dcor/dcov for sigma
      R <- lav_deriv_cov2cor(sigma.hat, num.idx = lavmodel@num.idx[[g]])
      theta.var.idx <- lav_matrix_diagh_idx(nvar[g])
    }

    for (mm in mm.in.group) {
      mname <- names(lavmodel@GLIST)[mm]

      # skip empty ones
      if (!length(m.el.idx[[mm]])) next

      # get Delta columns for this model matrix
      if (representation == "LISREL") {
        # Sigma
        DELTA <- dxSigma <-
          lav_lisrel_dsigma_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]],
            delta = parameterization == "delta"
          )
        if (categorical && parameterization == "theta") {
          DELTA <- R %*% DELTA
        }

        if (categorical) {
          # reorder: first variances (of numeric), then covariances
          cov.idx <- lav_matrix_vech_idx(nvar[g])
          covd.idx <- lav_matrix_vech_idx(nvar[g], diagonal = FALSE)

          var.idx <- which(is.na(match(
            cov.idx,
            covd.idx
          )))[num.idx[[g]]]
          cor.idx <- match(covd.idx, cov.idx)

          DELTA <- rbind(
            DELTA[var.idx, , drop = FALSE],
            DELTA[cor.idx, , drop = FALSE]
          )
        }

        # correlation structure?
        if (!categorical && correlation) {
          rm.idx <- lav_matrix_diagh_idx(nvar[g])
          DELTA <- DELTA[-rm.idx, , drop = FALSE]
        }

        if (!categorical) {
          if (conditional.x) {
            # means/intercepts
            DELTA.mu <- lav_lisrel_dmu_dx(
              MLIST = GLIST[mm.in.group],
              m = mname,
              idx = m.el.idx[[mm]]
            )

            # slopes
            if (lavmodel@nexo[g] > 0L) {
              DELTA.pi <- lav_lisrel_dpi_dx(
                 MLIST = GLIST[mm.in.group],
                m = mname,
                idx = m.el.idx[[mm]]
              )

              if (lavmodel@multilevel) {
                DELTA <- rbind(DELTA.mu, DELTA.pi, DELTA)
              } else {
                # ATTENTION: we need to change the order here
                # lav_mvreg_scores_* uses 'Beta' where the
                # the intercepts are just the first row
                # using the col-major approach, we need to
                # interweave the intercepts with the slopes!

                nEls <- NROW(DELTA.mu) + NROW(DELTA.pi)
                # = (nexo + 1 int) * nvar

                # intercepts on top
                tmp <- rbind(DELTA.mu, DELTA.pi)
                # change row index
                row.idx <- lav_matrix_vec(matrix(seq.int(nEls),
                  nrow = lavmodel@nexo[g] + 1L,
                  ncol = lavmodel@nvar[g], byrow = TRUE
                ))
                DELTA.beta <- tmp[row.idx, , drop = FALSE]
                DELTA <- rbind(DELTA.beta, DELTA)
              }
            } else {
              DELTA <- rbind(DELTA.mu, DELTA)
            }
          } else if (!conditional.x && lavmodel@meanstructure) {
            DELTA.mu <- lav_lisrel_dmu_dx(
               MLIST = GLIST[mm.in.group],
              m = mname,
              idx = m.el.idx[[mm]]
            )
            DELTA <- rbind(DELTA.mu, DELTA)
          }
        } else if (categorical) {
          DELTA.th <- lav_lisrel_dth_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]],
            th.idx = th.idx[[g]],
            delta = TRUE
          )
          if (parameterization == "theta") {
            # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
            dDelta.dx <-
              (dxSigma[theta.var.idx, , drop = FALSE] *
                -0.5 / (dsigma * sqrt(dsigma)))
            dth.dDelta <-
              lav_lisrel_dth_dx(
                MLIST = GLIST[mm.in.group],
                m = "delta",
                idx = 1:nvar[g],
                th.idx = th.idx[[g]]
              )
            # add dth.dDelta %*% dDelta.dx
            no.num.idx <- which(th.idx[[g]] > 0)
            DELTA.th[no.num.idx, ] <-
              DELTA.th[no.num.idx, , drop = FALSE] +
              (dth.dDelta %*% dDelta.dx)[no.num.idx, , drop = FALSE]
          }
          if (conditional.x && lavmodel@nexo[g] > 0L) {
            DELTA.pi <-
              lav_lisrel_dpi_dx(
                MLIST = GLIST[mm.in.group],
                m = mname,
                idx = m.el.idx[[mm]]
              )
            if (parameterization == "theta") {
              dpi.dDelta <-
                lav_lisrel_dpi_dx(
                  MLIST = GLIST[mm.in.group],
                  m = "delta",
                  idx = 1:nvar[g]
                )
              # add dpi.dDelta %*% dDelta.dx
              no.num.idx <-
                which(!seq.int(1L, nvar[g]) %in% num.idx[[g]])
              no.num.idx <- rep(seq.int(0, nexo[g] - 1) * nvar[g],
                each = length(no.num.idx)
              ) + no.num.idx
              DELTA.pi[no.num.idx, ] <-
                DELTA.pi[no.num.idx, , drop = FALSE] +
                (dpi.dDelta %*% dDelta.dx)[no.num.idx, , drop = FALSE]
            }
            DELTA <- rbind(DELTA.th, DELTA.pi, DELTA)
          } else {
            DELTA <- rbind(DELTA.th, DELTA)
          }
        }
        if (group.w.free) {
          DELTA.gw <- lav_lisrel_dgw_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
          DELTA <- rbind(DELTA.gw, DELTA)
        }
      } else if (representation == "RAM") {
        DELTA <- dxSigma <-
          lav_ram_dsigma(
            m = mname,
            idx = m.el.idx[[mm]],
            MLIST = GLIST[mm.in.group]
          )
        if (lavmodel@meanstructure) {
          DELTA.mu <- lav_ram_dmu(
            m = mname,
            idx = m.el.idx[[mm]],
            MLIST = GLIST[mm.in.group]
          )
          DELTA <- rbind(DELTA.mu, DELTA)
        }
      } else {
        lav_msg_stop(gettextf("representation %s not implemented yet",
                              representation))
      }
      Delta.group[, x.el.idx[[mm]]] <- DELTA
    } # mm

    # if type == "free" take care of equality constraints
    if (type == "free" && ceq.simple && lavmodel@ceq.simple.only) {
      Delta.group <- Delta.group %*% lavmodel@ceq.simple.K
    }

    Delta[[g]] <- Delta.group
  } # g

  # if multilevel, rbind levels within group
  if (lavmodel@multilevel) {
    DELTA <- vector("list", length = lavmodel@ngroups)
    for (g in 1:lavmodel@ngroups) {
      DELTA[[g]] <- rbind(
        Delta[[(g - 1) * 2 + 1]],
        Delta[[(g - 1) * 2 + 2]]
      )
    }
    Delta <- DELTA
  }

  Delta
}

lav_model_ddelta_dx <- function(lavmodel = NULL, GLIST = NULL, target = "lambda",
                           ceq.simple = FALSE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  representation <- lavmodel@representation
  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  th.idx <- lavmodel@th.idx

  # number of columns in DELTA + m.el.idx/x.el.idx
  type <- "free"
  # if(type == "free") {
  if (lavmodel@ceq.simple.only) {
    NCOL <- lavmodel@nx.unco
  } else {
    NCOL <- lavmodel@nx.free
  }
  m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
  for (mm in 1:length(GLIST)) {
    m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
    if (lavmodel@ceq.simple.only) {
      x.el.idx[[mm]] <- lavmodel@x.unco.idx[[mm]]
    } else {
      x.el.idx[[mm]] <- lavmodel@x.free.idx[[mm]]
    }
    # handle symmetric matrices
    if (lavmodel@isSymmetric[mm]) {
      # since we use 'x.free.idx', only symmetric elements
      # are duplicated (not the equal ones, only in x.free.free)
      dix <- duplicated(x.el.idx[[mm]])
      if (any(dix)) {
        m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
        x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
      }
    }
  }
  # } else {
  #    NCOL <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
  # }

  # compute Delta per group
  Delta <- vector("list", length = nblocks)
  for (g in 1:nblocks) {
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    Delta.group <- NULL
    for (mm in mm.in.group) {
      mname <- names(lavmodel@GLIST)[mm]

      # skip empty ones
      if (!length(m.el.idx[[mm]])) next

      # get Delta columns for this model matrix
      if (representation == "LISREL") {
        if (target == "lambda") {
         DELTA <- lav_lisrel_dlambda_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "th") {
          DELTA <- lav_lisrel_dth_dx(
             MLIST = GLIST[mm.in.group], m = mname, th.idx = th.idx[[g]],
            idx = m.el.idx[[mm]],
            delta = TRUE
          )
        } else if (target == "mu") {
          DELTA <- lav_lisrel_dmu_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "nu") {
          DELTA <- lav_lisrel_dnu_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "tau") {
          DELTA <- lav_lisrel_dtau_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "theta") {
          DELTA <- lav_lisrel_dtheta_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "gamma") {
          DELTA <- lav_lisrel_dgamma_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "beta") {
          DELTA <- lav_lisrel_dbeta_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "alpha") {
          DELTA <- lav_lisrel_dalpha_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "psi") {
          DELTA <- lav_lisrel_dpsi_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]]
          )
        } else if (target == "sigma") {
          DELTA <- lav_lisrel_dsigma_dx(
            MLIST = GLIST[mm.in.group],
            m = mname,
            idx = m.el.idx[[mm]],
            delta = TRUE
          )
        } else {
          lav_msg_stop(gettextf("target %s not implemented yet", target))
        }

        # initialize?
        if (is.null(Delta.group)) {
          Delta.group <- matrix(0, nrow = nrow(DELTA), ncol = NCOL)
        }
        Delta.group[, x.el.idx[[mm]]] <- DELTA
      }
    } # mm

    if (type == "free" && ceq.simple && lavmodel@ceq.simple.only) {
      Delta.group <- Delta.group %*% lavmodel@ceq.simple.K
    }

    Delta[[g]] <- Delta.group
  } # g

  Delta
}

# single block, continuous data/LISREL only
# new approach (0.6-22):
# - only the free elements per model matrix
# - in one go
lav_model_delta_lisrel <- function(lavmodel, block = 1L) {

  stopifnot(lavmodel@representation == "LISREL",
            !lavmodel@conditional.x,
            lavmodel@parameterization == "delta")

  # model matrices for this block
  mm.in.group <-
    seq_len(lavmodel@nmat[block]) + cumsum(c(0, lavmodel@nmat))[block]
  MLIST <- lavmodel@GLIST[mm.in.group]

  delta.flag <- beta.flag <- meanstructure <- categorical <- FALSE
  if (!is.null(MLIST$delta) && any(MLIST$delta[,1] != 1)) {
    delta.flag <- TRUE
  }
  if (!is.null(MLIST$beta)) {
    beta.flag <- TRUE
  }
  if (lavmodel@meanstructure) {
    meanstructure <- lavmodel@meanstructure
  }
  if (lavmodel@categorical) {
    categorical <- lavmodel@categorical
  }

  # model matrices in this block
  mnames <- names(MLIST)

  mm.lambda.idx <- which(mnames == "lambda")
  x.lambda.idx <- lavmodel@x.free.idx[[mm.lambda.idx]]
  m.lambda.idx <- lavmodel@m.free.idx[[mm.lambda.idx]]
  n_lam <- length(m.lambda.idx)

  mm.psi.idx <- which(mnames == "psi")
  tmp <- lavmodel@x.free.idx[[mm.psi.idx]]
  x.psi.idx <- x.psi.idx <- tmp[!duplicated(tmp)]
  m.psi.idx <- lavmodel@m.free.idx[[mm.psi.idx]][!duplicated(tmp)]
  n_psi <- length(m.psi.idx)


  mm.theta.idx <- which(mnames == "theta")
  tmp <- lavmodel@x.free.idx[[mm.theta.idx]]
  x.theta.idx <- tmp[!duplicated(tmp)]
  m.theta.idx <- lavmodel@m.free.idx[[mm.theta.idx]][!duplicated(tmp)]
  n_the <- length(m.theta.idx)

  n_bet <- 0L
  x.beta.idx <- integer(0L)
  if (beta.flag) {
    mm.beta.idx <- which(mnames == "beta")
    x.beta.idx <- lavmodel@x.free.idx[[mm.beta.idx]]
    m.beta.idx <- lavmodel@m.free.idx[[mm.beta.idx]]
    n_bet <- length(m.beta.idx)
  }

  n_del <- 0L
  x.delta.idx <- integer(0L)
  if (delta.flag) {
    mm.delta.idx <- which(mnames == "delta")
    x.delta.idx <- lavmodel@x.free.idx[[mm.delta.idx]]
    m.delta.idx <- lavmodel@m.free.idx[[mm.delta.idx]]
    n_del <- length(m.delta.idx)
  }

  # in case !meanstructure or !categorical
  n_nu  <- 0L;  m.nu.idx    <- integer(0L);  x.nu.idx    <- integer(0L)
  n_alp <- 0L;  m.alpha.idx <- integer(0L);  x.alpha.idx <- integer(0L)

  if (meanstructure) {
    mm.nu.idx <- which(mnames == "nu")
    x.nu.idx <- lavmodel@x.free.idx[[mm.nu.idx]]
    m.nu.idx <- lavmodel@m.free.idx[[mm.nu.idx]]
    n_nu <- length(m.nu.idx)

    mm.alpha.idx <- which(mnames == "alpha")
    x.alpha.idx <- lavmodel@x.free.idx[[mm.alpha.idx]]
    m.alpha.idx <- lavmodel@m.free.idx[[mm.alpha.idx]]
    n_alp <- length(m.alpha.idx)
  }

  if (categorical) {
    mm.th.idx <- which(mnames == "tau")
    x.th.idx <- lavmodel@x.free.idx[[mm.th.idx]]
    m.th.idx <- lavmodel@m.free.idx[[mm.th.idx]]
    n_th <- length(m.th.idx)
  }

  nvar <- nrow(MLIST$lambda)
  nfac <- ncol(MLIST$lambda)
  pstar <- nvar * (nvar + 1L) / 2L

  # precompute
  if (beta.flag) {
    A <- lav_lisrel_ibinv(MLIST)
    M <- MLIST$lambda %*% A
    G <- A %*% MLIST$psi %*% t(A)
    LG <- MLIST$lambda %*% G
  } else {
    M <- MLIST$lambda
    G <- MLIST$psi
    LG <- MLIST$lambda %*% MLIST$psi
  }

  if (delta.flag) {
    delta  <- as.vector(MLIST$delta)
  }

  # vech structure for Sigma
  r_s <- lav_matrix_vech_row_idx(nvar)
  c_s <- lav_matrix_vech_col_idx(nvar)
  sigma_lut <- lav_matrix_vech_reverse(seq_len(pstar))

  if (delta.flag) {
    # scaling weights: delta[r] * delta[s] for each vech position
    delta_weight <- delta[r_s] * delta[c_s]
  }

  # helper: vec index -> (row, col)
  vec2rc <- function(idx, nr) {
    cbind(row = (idx - 1L) %% nr + 1L,
          col = (idx - 1L) %/% nr + 1L)
  }

  # prepare output matrix for sigma
  n_free <- n_lam + n_bet + n_psi + n_the + n_del
  jac_sigma <- matrix(0, pstar, n_free)

  col <- 1L
  # Lambda [i,j]: dS0[r,s] = LG[s,j]*I(r==i) + LG[r,j]*I(s==i)
  if (n_lam > 0L) {
    rc  <- vec2rc(m.lambda.idx, nvar);  i_v <- rc[, 1L];  j_v <- rc[, 2L]

    T1 <- LG[c_s, j_v, drop = FALSE] * outer(r_s, i_v, `==`)
    T2 <- LG[r_s, j_v, drop = FALSE] * outer(c_s, i_v, `==`)

    DX <- T1 + T2
    if (delta.flag) {
      DX <- DX * delta_weight
    }
    jac_sigma[, col:(col + n_lam - 1L)] <- DX
    col <- col + n_lam
  }

  # Beta [i,j]: dS0[r,s] = M[r,i]*LG[s,j] + LG[r,j]*M[s,i]
  if (n_bet > 0L) {
    rc  <- vec2rc(m.beta.idx, nfac);  i_v <- rc[, 1L];  j_v <- rc[, 2L]

    T1 <- M[r_s,  i_v, drop = FALSE] * LG[c_s, j_v, drop = FALSE]
    T2 <- LG[r_s, j_v, drop = FALSE] * M[c_s,  i_v, drop = FALSE]

    DX <- T1 + T2
    if (delta.flag) {
      DX <- DX * delta_weight
    }
    jac_sigma[, col:(col + n_bet - 1L)] <- DX
    col <- col + n_bet
  }

  # Psi [k,l] symmetric: dS0[r,s] = M[r,k]*M[s,l] + M[r,l]*M[s,k]
  # diagonal (k==l): halve
  if (n_psi > 0L) {
    rc  <- vec2rc(m.psi.idx, nfac)
    k_v <- pmax(rc[, 1L], rc[, 2L])
    l_v <- pmin(rc[, 1L], rc[, 2L])

    T1 <- M[r_s, k_v, drop = FALSE] * M[c_s, l_v, drop = FALSE]
    T2 <- M[r_s, l_v, drop = FALSE] * M[c_s, k_v, drop = FALSE]
    DX <- T1 + T2

    diag_mask <- (k_v == l_v)
    if (any(diag_mask)) {
      DX[, diag_mask] <- DX[, diag_mask] * 0.5
    }

    if (delta.flag) {
      DX <- DX * delta_weight
    }
    jac_sigma[, col:(col + n_psi - 1L)] <- DX
    col <- col + n_psi
  }

  # Theta [k,l] symmetric: dS = Delta %*% dTheta %*% Delta
  #  -> unit vector at vech pos (k,l)
  if (n_the > 0L) {
    rc  <- vec2rc(m.theta.idx, nvar)
    k_v <- pmax(rc[, 1L], rc[, 2L])
    l_v <- pmin(rc[, 1L], rc[, 2L])

    Jt <- matrix(0, pstar, n_the)
    vech_pos <- sigma_lut[cbind(k_v, l_v)]
    if (delta.flag) {
      Jt[cbind(vech_pos, seq_len(n_the))] <- delta_weight[vech_pos]
    } else {
      Jt[cbind(vech_pos, seq_len(n_the))] <- 1
    }

    jac_sigma[, col:(col + n_the - 1L)] <- Jt
    col <- col + n_the
  }

  # Delta[k] (diagonal only):
  # dS[r,s] = I(r==k)*Sigma0[k,s]*delta[s] + delta[r]*Sigma0[r,k]*I(s==k)
  if (n_del > 0L) {
    Sigma0 <- M %*% MLIST$psi %*% t(M) + MLIST$theta
    k_v <- m.delta.idx
    T1 <- Sigma0[c_s, k_v, drop = FALSE] * delta[c_s] * outer(r_s, k_v, `==`)
    T2 <- delta[r_s] * Sigma0[r_s, k_v, drop = FALSE] * outer(c_s, k_v, `==`)
    jac_sigma[, col:(col + n_del - 1L)] <- T1 + T2
  }

  # categorical: reorder
  if (categorical) {
    # reorder: first variances (of numeric), then covariances
    cov.idx <- lav_matrix_vech_idx(lavmodel@nvar[block])
    covd.idx <- lav_matrix_vech_idx(lavmodel@nvar[block], diagonal = FALSE)

    var.idx <- which(is.na(match(
      cov.idx,
      covd.idx
    )))[lavmodel@num.idx[[block]]]
    cor.idx <- match(covd.idx, cov.idx)

    jac_sigma <- rbind(
      jac_sigma[var.idx, , drop = FALSE],
      jac_sigma[cor.idx, , drop = FALSE]
    )
  }

  # correlation structure
  if (!categorical && lavmodel@correlation) {
    rm.idx <- lav_matrix_diagh_idx(lavmodel@nvar[block])
    jac_sigma <- jac_sigma[-rm.idx, , drop = FALSE]
  }

  jac_th   <- NULL
  nth_full <- 0L
  if (!categorical) {
    if (meanstructure) {
      # precompute: IB.inv * alpha
      if (beta.flag) {
        a <- drop(A %*% MLIST$alpha)
      } else {
        a <- MLIST$alpha
      }

      n_free <- n_nu + n_lam + n_bet + n_alp
      jac_mean <- matrix(0, nvar, n_free)

      col <- 1L
      # nu[i]:  dmu = e_i
      if (n_nu > 0L) {
        Jn <- matrix(0, nvar, n_nu)
        Jn[cbind(m.nu.idx, seq_len(n_nu))] <- 1
        jac_mean[, col:(col + n_nu - 1L)] <- Jn
        col <- col + n_nu
      }

      # Lambda[i,j]:  dmu[r] = a[j] * I(r == i)
      if (n_lam > 0L) {
        rc  <- vec2rc(m.lambda.idx, nvar); i_v <- rc[, 1L];  j_v <- rc[, 2L]

        Jl <- matrix(0, nvar, n_lam)
        Jl[cbind(i_v, seq_len(n_lam))] <- a[j_v]
        jac_mean[, col:(col + n_lam - 1L)] <- Jl
        col <- col + n_lam
      }

      # Beta[i,j]:  dmu[r] = M[r,i] * a[j]
      if (n_bet > 0L) {
        rc  <- vec2rc(m.beta.idx, nfac); i_v <- rc[, 1L];  j_v <- rc[, 2L]

        jac_mean[, col:(col + n_bet - 1L)] <-
          M[, i_v, drop = FALSE] * rep(a[j_v], each = nvar)
        col <- col + n_bet
      }

      # alpha[k]:  dmu = M[,k]
      if (n_alp > 0L) {
        jac_mean[, col:(col + n_alp - 1L)] <- M[, m.alpha.idx, drop = FALSE]
      }
    } # meanstructure

  # categorical
  } else if (categorical) {

    th_idx_g <- lavmodel@th.idx[[block]]
    nth_full  <- length(th_idx_g)

    # v_slot[t] = observed-variable index for TH element t
    nlev_v <- tabulate(th_idx_g, nbins = nvar)
    nlev_v[nlev_v == 0L] <- 1L
    v_slot <- rep(seq_len(nvar), times = nlev_v)

    # per-slot delta scale
    if (delta.flag) {
      delta_star <- delta[v_slot]
    } else {
      delta_star <- rep(1, nth_full)
    }

    # ord_slots[k] = position in TH vector that corresponds to tau row k
    ord_slots <- which(th_idx_g > 0L)

    # a_vec = IB.inv * alpha  (for lambda/beta derivatives)
    if (!is.null(MLIST$alpha)) {
      a_vec <- if (beta.flag) {
        as.vector(A %*% MLIST$alpha)
      } else {
        as.vector(MLIST$alpha)
      }
    } else {
      a_vec <- rep(0, nfac)
    }

    # nu_vec (intercepts; zero for purely ordinal models)
    nu_vec <- if (!is.null(MLIST$nu)) as.vector(MLIST$nu) else rep(0, nvar)

    # pi0 = nu + Lambda * IB.inv * alpha = nu + M * alpha
    if (!is.null(MLIST$alpha)) {
      pi0 <- nu_vec + as.vector(M %*% as.vector(MLIST$alpha))
    } else {
      pi0 <- nu_vec
    }

    # tau_full: tau values at ordinal TH slots, 0 at numeric TH slots
    tau_full <- numeric(nth_full)
    if (!is.null(MLIST$tau) && n_th > 0L) {
      tau_full[ord_slots] <- as.vector(MLIST$tau)
    }

    # unscaled TH: TH0[t] = tau_full[t] - pi0[v(t)]
    TH0 <- tau_full - pi0[v_slot]

    # build jac_th  (nth_full rows, n_free_th columns)
    # column order: tau | delta | nu | lambda | beta | alpha
    n_free_th <- n_th + n_del + n_nu + n_lam + n_bet + n_alp
    jac_th <- matrix(0, nth_full, n_free_th)
    col_th <- 1L

    # tau[k]: delta_star[t] * I(t == ord_slots[k])
    if (n_th > 0L) {
      t_slots <- ord_slots[m.th.idx]
      Jt <- matrix(0, nth_full, n_th)
      Jt[cbind(t_slots, seq_len(n_th))] <- delta_star[t_slots]
      jac_th[, col_th:(col_th + n_th - 1L)] <- Jt
      col_th <- col_th + n_th
    }

    # delta[k]: I(v_slot[t]==k) * TH0[t]
    if (n_del > 0L) {
      jac_th[, col_th:(col_th + n_del - 1L)] <-
        outer(v_slot, m.delta.idx, `==`) * TH0
      col_th <- col_th + n_del
    }

    # nu[i]: -delta_star[t] * I(v_slot[t]==i)
    if (n_nu > 0L) {
      jac_th[, col_th:(col_th + n_nu - 1L)] <-
        -outer(v_slot, m.nu.idx, `==`) * delta_star
      col_th <- col_th + n_nu
    }

    # lambda[i,j]: -delta_star[t]*I(v_slot[t]==i)*a_vec[j]
    if (n_lam > 0L) {
      rc  <- vec2rc(m.lambda.idx, nvar); i_v <- rc[, 1L]; j_v <- rc[, 2L]
      jac_th[, col_th:(col_th + n_lam - 1L)] <-
        -(outer(v_slot, i_v, `==`) * delta_star) *
         matrix(a_vec[j_v], nth_full, n_lam, byrow = TRUE)
      col_th <- col_th + n_lam
    }

    # beta[i,j]: -delta_star[t]*M[v_slot[t],i]*a_vec[j]
    if (n_bet > 0L) {
      rc  <- vec2rc(m.beta.idx, nfac); i_v <- rc[, 1L]; j_v <- rc[, 2L]
      jac_th[, col_th:(col_th + n_bet - 1L)] <-
        -delta_star * M[v_slot, i_v, drop = FALSE] *
         matrix(a_vec[j_v], nth_full, n_bet, byrow = TRUE)
      col_th <- col_th + n_bet
    }

    # alpha[k]: -delta_star[t] * M[v_slot[t], k]
    if (n_alp > 0L) {
     jac_th[, col_th:(col_th + n_alp - 1L)] <-
       -delta_star * M[v_slot, m.alpha.idx, drop = FALSE]
    }
  }

  # sigma
  out <- matrix(0, nrow = nrow(jac_sigma), ncol = lavmodel@nx.free)
  el.idx <- c(x.lambda.idx, x.beta.idx, x.psi.idx, x.theta.idx, x.delta.idx)
  out[, el.idx] <- jac_sigma

  # meanstructure
  if (!categorical && meanstructure) {
    # right order
    el.idx <- c(x.nu.idx, x.lambda.idx, x.beta.idx, x.alpha.idx)
    outm <- matrix(0, nrow = nrow(jac_mean), ncol = lavmodel@nx.free)
    outm[, el.idx] <- jac_mean
    out <- rbind(outm, out)
  } else if (categorical && !is.null(jac_th)) {
    el.idx_th <- c(
      if (n_th  > 0L) x.th.idx     else integer(0L),
      if (n_del > 0L) x.delta.idx  else integer(0L),
      if (n_nu  > 0L) x.nu.idx     else integer(0L),
      if (n_lam > 0L) x.lambda.idx else integer(0L),
      if (n_bet > 0L) x.beta.idx   else integer(0L),
      if (n_alp > 0L) x.alpha.idx  else integer(0L)
    )
    out_th <- matrix(0, nth_full, lavmodel@nx.free)
    out_th[, el.idx_th] <- jac_th
    out <- rbind(out_th, out)
  }

  out
}


lav_model_omega <- function(Sigma.hat = NULL, Mu.hat = NULL,
                         lavsamplestats = NULL, estimator = "ML",
                         meanstructure = FALSE, conditional.x = FALSE,
						 correlation = FALSE) {
  # nblocks
  nblocks <- length(Sigma.hat)

  Omega <- vector("list", length = nblocks)
  Omega.mu <- vector("list", length = nblocks)

  for (g in 1:nblocks) {
    # ML
    if (estimator %in% c("ML", "REML", "catML")) {
      if (attr(Sigma.hat[[g]], "po") == FALSE) {
        # FIXME: WHAT IS THE BEST THING TO DO HERE??
        # CURRENTLY: stop
        lav_msg_warn(gettext(
          "lav_model_gradient: Sigma.hat is not positive definite\n"))
        Sigma.hat.inv <- MASS::ginv(Sigma.hat[[g]])
      } else {
        Sigma.hat.inv <- attr(Sigma.hat[[g]], "inv")
      }

      if (!lavsamplestats@missing.flag) { # complete data
        if (meanstructure) {
          if (conditional.x) {
            diff <- lavsamplestats@res.int[[g]] - Mu.hat[[g]]
            W.tilde <- lavsamplestats@res.cov[[g]] + tcrossprod(diff)
          } else {
            diff <- lavsamplestats@mean[[g]] - Mu.hat[[g]]
            W.tilde <- lavsamplestats@cov[[g]] + tcrossprod(diff)
          }
          # Browne 1995 eq 4.55
          Omega.mu[[g]] <- t(t(diff) %*% Sigma.hat.inv)
          Omega[[g]] <-
            (Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*%
              Sigma.hat.inv)
        } else {
          if (conditional.x) {
            W.tilde <- lavsamplestats@res.cov[[g]]
          } else {
            W.tilde <- lavsamplestats@cov[[g]]
          }
          Omega[[g]] <-
            (Sigma.hat.inv %*% (W.tilde - Sigma.hat[[g]]) %*%
              Sigma.hat.inv)
        }
      } else { # missing data
        M <- lavsamplestats@missing[[g]]

        nvar <- ncol(lavsamplestats@cov[[g]])
        OMEGA <- matrix(0, nvar, nvar)
        OMEGA.MU <- matrix(0, nvar, 1)

        for (p in 1:length(M)) {
          SX <- M[[p]][["SY"]]
          MX <- M[[p]][["MY"]]
          nobs <- M[[p]][["freq"]]
          var.idx <- M[[p]][["var.idx"]]

          Sigma.inv <- try(chol2inv(chol(Sigma.hat[[g]][var.idx, var.idx])),
                          silent = TRUE)
          if (inherits(Sigma.inv, "try-error")) {
            Sigma.inv <- MASS::ginv(Sigma.hat[[g]][var.idx, var.idx])
          }
          Mu <- Mu.hat[[g]][var.idx]
          W.tilde <- SX + tcrossprod(MX - Mu)

          OMEGA.MU[var.idx, 1] <-
            (OMEGA.MU[var.idx, 1] + nobs / lavsamplestats@ntotal *
              t(t(MX - Mu) %*% Sigma.inv))

          OMEGA[var.idx, var.idx] <-
            (OMEGA[var.idx, var.idx] + nobs / lavsamplestats@ntotal *
              (Sigma.inv %*%
                (W.tilde - Sigma.hat[[g]][var.idx, var.idx]) %*%
                Sigma.inv))
        }
        Omega.mu[[g]] <- OMEGA.MU
        Omega[[g]] <- OMEGA
      } # missing

      # GLS
    } else if (estimator == "GLS") {
      W.inv <- lavsamplestats@icov[[g]]
      W <- lavsamplestats@cov[[g]]
      Omega[[g]] <- (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@nobs[[g]] *
        (W.inv %*% (W - Sigma.hat[[g]]) %*% W.inv)
      if (meanstructure) {
        diff <- as.matrix(lavsamplestats@mean[[g]] - Mu.hat[[g]])
        Omega.mu[[g]] <- t(t(diff) %*% W.inv)
      }
    }

    # new in 0.6-18
	if(correlation) {
	    diag(Omega[[g]]) <- 0
	}
  } # g

  if (meanstructure) attr(Omega, "mu") <- Omega.mu

  Omega
}

lav_model_gradient_dd <- function(lavmodel, g_list = NULL, group = 1L) {
  if (is.null(g_list)) g_list <- lavmodel@GLIST

  #### FIX th + mu!!!!!
  delta_lambda <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "lambda")[[group]]
  delta_tau <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "tau")[[group]]
  delta_nu <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "nu")[[group]]
  delta_theta <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "theta")[[group]]
  delta_beta <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "beta")[[group]]
  delta_psi <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "psi")[[group]]
  delta_alpha <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "alpha")[[group]]
  delta_gamma <-
    lav_model_ddelta_dx(lavmodel, GLIST = g_list, target = "gamma")[[group]]

  ov_y_dummy_ov_idx <- lavmodel@ov.y.dummy.ov.idx[[group]]
  ov_x_dummy_ov_idx <- lavmodel@ov.x.dummy.ov.idx[[group]]
  ov_y_dummy_lv_idx <- lavmodel@ov.y.dummy.lv.idx[[group]]
  ov_x_dummy_lv_idx <- lavmodel@ov.x.dummy.lv.idx[[group]]
  ov_dummy_idx <- c(ov_y_dummy_ov_idx, ov_x_dummy_ov_idx)
  lv_dummy_idx <- c(ov_y_dummy_lv_idx, ov_x_dummy_lv_idx)
  num_idx <- lavmodel@num.idx[[group]]

  # fix Delta's...
  mm_in_group <- 1:lavmodel@nmat[group] + cumsum(c(0, lavmodel@nmat))[group]
  m_list <- g_list[mm_in_group]

  dd <- list()
  nvar <- lavmodel@nvar
  nfac <- ncol(m_list$lambda) - length(lv_dummy_idx)

  # dd$theta
  theta_idx <- lav_matrix_diagh_idx(nvar)
  dd$theta <- delta_theta[theta_idx, , drop = FALSE]
  if (length(ov_dummy_idx) > 0L) {
    psi_idx <- lav_matrix_diagh_idx(ncol(m_list$psi))[lv_dummy_idx]
    dd$theta[ov_dummy_idx, ] <- delta_psi[psi_idx, , drop = FALSE]
  }
  # num only? FIXME or just all of them?
  dd$theta <- dd$theta[num_idx, , drop = FALSE]

  # dd$nu
  dd$nu <- delta_nu
  if (length(ov_dummy_idx) > 0L) {
    dd$nu[ov_dummy_idx, ] <- delta_alpha[lv_dummy_idx, ]
  }
  dd$nu <- dd$nu[num_idx, , drop = FALSE] # needed?

  # dd$lambda
  nr <- nvar
  nc <- nfac
  lambda_idx <- nr * ((1:nc) - 1L) + rep(1:nvar, each = nc)
  dd$lambda <- delta_lambda[lambda_idx, , drop = FALSE]
  if (length(ov_dummy_idx) > 0L) {
    nr <- nrow(m_list$beta)
    nc <- nfac # only the first 1:nfac columns
    # beta_idx <- rep(nr*((1:nc) - 1L), each=length(lv_dummy_idx)) +
    #                                 rep(lv_dummy_idx, times=nc) ## FIXME
    beta_idx <- rep(nr * ((1:nc) - 1L),
           times = length(lv_dummy_idx)) + rep(lv_dummy_idx, each = nc)

    # l_idx <- inr*((1:nc) - 1L) + rep(ov_dummy_idx, each=nc) ## FIXME
    # l_idx <- rep(nr*((1:nc) - 1L), each=length(ov_dummy_idx)) +
    #                                rep(ov_dummy_idx, times=nc)
    l_idx <- rep(nr * ((1:nc) - 1L),
                times = length(ov_dummy_idx)) + rep(ov_dummy_idx, each = nc)
    dd$lambda[match(l_idx, lambda_idx), ] <-
                delta_beta[beta_idx, , drop = FALSE]
  }

  # dd$KAPPA
  dd$kappa <- delta_gamma
  if (length(ov_dummy_idx) > 0L) {
    nr <- nrow(m_list$gamma)
    nc <- ncol(m_list$gamma)
    kappa_idx <- nr * ((1:nc) - 1L) + rep(lv_dummy_idx, each = nc)
    dd$kappa <- dd$kappa[kappa_idx, , drop = FALSE]
  }

  # dd$GAMMA
  if (!is.null(m_list$gamma)) {
    nr <- nrow(m_list$gamma)
    nc <- ncol(m_list$gamma)
    lv_idx <- 1:nfac
    # MUST BE ROWWISE!
    gamma_idx <- rep(nr * ((1:nc) - 1L),
                 times = length(lv_idx)) + rep(lv_idx, each = nc)
    dd$gamma <- delta_gamma[gamma_idx, , drop = FALSE]
  }

  # dd$BETA
  if (!is.null(m_list$beta)) {
    nr <- nc <- nrow(m_list$beta)
    lv_idx <- 1:nfac
    # MUST BE ROWWISE!
    beta_idx <- rep(nr * ((1:nfac) - 1L),
                times = nfac) + rep(lv_idx, each = nfac)
    dd$beta <- delta_beta[beta_idx, , drop = FALSE]
  }

  ## dd$psi
  dd$psi <- delta_psi
  if (length(lv_dummy_idx) > 0L) {
    nr <- nc <- nrow(m_list$psi)
    lv_idx <- 1:nfac
    # MUST BE ROWWISE!
    psi_idx <- rep(nr * ((1:nfac) - 1L),
               times = nfac) + rep(lv_idx, each = nfac)

    dd$psi <- dd$psi[psi_idx, , drop = FALSE]
  }

  ## dd$tau
  if (!is.null(m_list$tau)) {
    dd$tau <- delta_tau
  }

  dd
}
