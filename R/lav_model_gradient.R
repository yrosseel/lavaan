# model gradient

lav_model_gradient <- function(lavmodel = NULL,
                               GLIST = NULL,
                               lavsamplestats = NULL,
                               lavdata = NULL,
                               lavcache = NULL,
                               type = "free",
                               verbose = FALSE,
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
  if (.hasSlot(lavmodel, "estimator.args")) {
    estimator.args <- lavmodel@estimator.args
  } else {
    estimator.args <- list()
  }

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
    #    Sigma.hat <- computeSigmaHatJoint(lavmodel = lavmodel,
    #                     GLIST = GLIST,
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

    if (meanstructure) {
      # if(conditional.x) {
      #    Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
      # } else {
      Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
      # }
    }

    if (categorical) {
      TH <- computeTH(lavmodel = lavmodel, GLIST = GLIST)
    }

    if (conditional.x) {
      PI <- computePI(lavmodel = lavmodel, GLIST = GLIST)
    } else if (estimator == "PML") {
      PI <- vector("list", length = lavmodel@nblocks)
    }

    if (group.w.free) {
      GW <- computeGW(lavmodel = lavmodel, GLIST = GLIST)
    }
  } else if (estimator == "DLS" && estimator.args$dls.GammaNT == "model") {
    Sigma.hat <- computeSigmaHat(
      lavmodel = lavmodel, GLIST = GLIST,
      extra = FALSE
    )
    Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
  } else if (estimator == "MML") {
    TH <- computeTH(lavmodel = lavmodel, GLIST = GLIST)
    THETA <- computeTHETA(lavmodel = lavmodel, GLIST = GLIST)
    GW <- computeGW(lavmodel = lavmodel, GLIST = GLIST)
  }

  # four approaches (FIXME!!!! merge this!)
  # - ML approach: using Omega (and Omega.mu)
  #    Omega = 'POST' = Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)
  #   (still 2x faster than Delta method)
  # - WLS/DWLS/GLS: using Delta + WLS.V; support for fixed.x, conditional.x
  # - (ML)/NTRLS: using Delta, no support for fixed.x, conditional.x
  # - PML/FML/MML: custom

  # 1. ML approach
  if ((estimator == "ML" || estimator == "REML" || estimator == "catML") &&
    lavdata@nlevels == 1L &&
    !lavmodel@conditional.x) {
    if (meanstructure) {
      Omega <- computeOmega(
        Sigma.hat = Sigma.hat, Mu.hat = Mu.hat,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = TRUE,
        conditional.x = conditional.x
      )
      Omega.mu <- attr(Omega, "mu")
    } else {
      Omega <- computeOmega(
        Sigma.hat = Sigma.hat, Mu.hat = NULL,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = FALSE,
        conditional.x = conditional.x
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
        DX.group <- derivative.F.LISREL(
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
      if (.hasSlot(lavmodel, "ceq.simple.only") &&
        lavmodel@ceq.simple.only) { # new in 0.6-11
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
      Delta <- computeDelta(
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
          group.dx <- -1 * crossprod(
            Delta[[g]],
            crossprod(WLS.V, diff)
          )
        } else if (estimator == "NTRLS") {
          stopifnot(!conditional.x)
          # WLS.V <- lav_samplestats_Gamma_inverse_NT(
          #         ICOV = attr(Sigma.hat[[g]],"inv")[,,drop=FALSE],
          #         COV            = Sigma.hat[[g]][,,drop=FALSE],
          #         MEAN           = Mu.hat[[g]],
          #         x.idx          = lavsamplestats@x.idx[[g]],
          #         fixed.x        = fixed.x,
          #         conditional.x  = conditional.x,
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
      dx <- lav_model_x2GLIST(
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
      Delta <- computeDelta(
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
      dx <- lav_model_x2GLIST(
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
      Delta <- computeDelta(
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
      Delta <- computeDelta(
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
          d1 <- pml_deriv1(
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
          d1 <- pml_deriv1(
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
        d1 <- fml_deriv1(
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
    # est.prop <- unlist( computeGW(lavmodel = lavmodel, GLIST = GLIST) )
    # obs.prop <- unlist(lavsamplestats@group.w)
    # FIXME: G2 based -- ML and friends only!!
    # dx.GW <- - (obs.prop - est.prop)

    # poisson version
    est.freq <- exp(unlist(computeGW(lavmodel = lavmodel, GLIST = GLIST)))
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
# computeDeltaNumerical <- function(lavmodel = NULL, GLIST = NULL, g = 1L) {
#
#     # state or final?
#    if(is.null(GLIST)) GLIST <- lavmodel@GLIST
#
#    compute.moments <- function(x) {
#        GLIST <- lav_model_x2GLIST(lavmodel = lavmodel, x=x, type="free")
#        Sigma.hat <- computeSigmaHat(lavmodel = lavmodel, GLIST = GLIST)
#         S.vec <- lav_matrix_vech(Sigma.hat[[g]])
#         if(lavmodel@meanstructure) {
#             Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST=GLIST)
#             out <- c(Mu.hat[[g]], S.vec)
#         } else {
#             out <- S.vec
#         }
#         out
#     }
#
#     x <- lav_model_get_parameters(lavmodel = lavmodel, GLIST=GLIST, type="free")
#     Delta <- lav_func_jacobian_complex(func=compute.moments, x = x)
#
#     Delta
# }


### FIXME: should we here also:
###        - weight for groups? (no, for now)
###        - handle equality constraints? (yes, for now)
computeDelta <- function(lavmodel = NULL, GLIST. = NULL,
                         m.el.idx. = NULL, x.el.idx. = NULL,
                         ceq.simple = FALSE,
                         force.conditional.x.false = FALSE) {
  representation <- lavmodel@representation
  categorical <- lavmodel@categorical
  if (.hasSlot(lavmodel, "correlation")) {
    correlation <- lavmodel@correlation
  } else {
    correlation <- FALSE
  }
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
    if (.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      NCOL <- lavmodel@nx.unco
    } else {
      NCOL <- lavmodel@nx.free
    }
    m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
    for (mm in 1:length(GLIST)) {
      m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
      if (.hasSlot(lavmodel, "ceq.simple.only") &&
        lavmodel@ceq.simple.only) {
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
      sigma.hat <- computeSigmaHat.LISREL(
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
          derivative.sigma.LISREL(
            m = mname,
            idx = m.el.idx[[mm]],
            MLIST = GLIST[mm.in.group],
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
            DELTA.mu <- derivative.mu.LISREL(
              m = mname,
              idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
            )

            # slopes
            if (lavmodel@nexo[g] > 0L) {
              DELTA.pi <- derivative.pi.LISREL(
                m = mname,
                idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
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
            DELTA.mu <- derivative.mu.LISREL(
              m = mname,
              idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
            )
            DELTA <- rbind(DELTA.mu, DELTA)
          }
        } else if (categorical) {
          DELTA.th <- derivative.th.LISREL(
            m = mname,
            idx = m.el.idx[[mm]],
            th.idx = th.idx[[g]],
            MLIST = GLIST[mm.in.group],
            delta = TRUE
          )
          if (parameterization == "theta") {
            # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
            dDelta.dx <-
              (dxSigma[theta.var.idx, , drop = FALSE] *
                -0.5 / (dsigma * sqrt(dsigma)))
            dth.dDelta <-
              derivative.th.LISREL(
                m = "delta",
                idx = 1:nvar[g],
                MLIST = GLIST[mm.in.group],
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
              derivative.pi.LISREL(
                m = mname,
                idx = m.el.idx[[mm]],
                MLIST = GLIST[mm.in.group]
              )
            if (parameterization == "theta") {
              dpi.dDelta <-
                derivative.pi.LISREL(
                  m = "delta",
                  idx = 1:nvar[g],
                  MLIST = GLIST[mm.in.group]
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
          DELTA.gw <- derivative.gw.LISREL(
            m = mname,
            idx = m.el.idx[[mm]],
            MLIST = GLIST[mm.in.group]
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
    if (type == "free" && ceq.simple &&
      .hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      Delta.group <- Delta.group %*% lavmodel@ceq.simple.K
    }

    Delta[[g]] <- Delta.group
  } # g

  # if multilevel, rbind levels within group
  if (.hasSlot(lavmodel, "multilevel") && lavmodel@multilevel) {
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

computeDeltaDx <- function(lavmodel = NULL, GLIST = NULL, target = "lambda",
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
  if (.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
    NCOL <- lavmodel@nx.unco
  } else {
    NCOL <- lavmodel@nx.free
  }
  m.el.idx <- x.el.idx <- vector("list", length = length(GLIST))
  for (mm in 1:length(GLIST)) {
    m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
    if (.hasSlot(lavmodel, "ceq.simple.only") &&
      lavmodel@ceq.simple.only) {
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
          DELTA <- derivative.lambda.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "th") {
          DELTA <- derivative.th.LISREL(
            m = mname, th.idx = th.idx[[g]],
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group],
            delta = TRUE
          )
        } else if (target == "mu") {
          DELTA <- derivative.mu.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "nu") {
          DELTA <- derivative.nu.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "tau") {
          DELTA <- derivative.tau.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "theta") {
          DELTA <- derivative.theta.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "gamma") {
          DELTA <- derivative.gamma.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "beta") {
          DELTA <- derivative.beta.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "alpha") {
          DELTA <- derivative.alpha.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "psi") {
          DELTA <- derivative.psi.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group]
          )
        } else if (target == "sigma") {
          DELTA <- derivative.sigma.LISREL(
            m = mname,
            idx = m.el.idx[[mm]], MLIST = GLIST[mm.in.group],
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

    if (type == "free" && ceq.simple &&
      .hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      Delta.group <- Delta.group %*% lavmodel@ceq.simple.K
    }

    Delta[[g]] <- Delta.group
  } # g

  Delta
}

computeOmega <- function(Sigma.hat = NULL, Mu.hat = NULL,
                         lavsamplestats = NULL, estimator = "ML",
                         meanstructure = FALSE, conditional.x = FALSE) {
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

          Sigma.inv <- inv.chol(Sigma.hat[[g]][var.idx, var.idx],
            logdet = FALSE
          )
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
  } # g

  if (meanstructure) attr(Omega, "mu") <- Omega.mu

  Omega
}

lav_model_gradient_DD <- function(lavmodel, GLIST = NULL, group = 1L) {
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  #### FIX th + mu!!!!!
  Delta.lambda <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "lambda")[[group]]
  Delta.tau <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "tau")[[group]]
  Delta.nu <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "nu")[[group]]
  Delta.theta <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "theta")[[group]]
  Delta.beta <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "beta")[[group]]
  Delta.psi <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "psi")[[group]]
  Delta.alpha <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "alpha")[[group]]
  Delta.gamma <- computeDeltaDx(lavmodel, GLIST = GLIST, target = "gamma")[[group]]

  ov.y.dummy.ov.idx <- lavmodel@ov.y.dummy.ov.idx[[group]]
  ov.x.dummy.ov.idx <- lavmodel@ov.x.dummy.ov.idx[[group]]
  ov.y.dummy.lv.idx <- lavmodel@ov.y.dummy.lv.idx[[group]]
  ov.x.dummy.lv.idx <- lavmodel@ov.x.dummy.lv.idx[[group]]
  ov.dummy.idx <- c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)
  th.idx <- lavmodel@th.idx[[group]]
  num.idx <- lavmodel@num.idx[[group]]
  ord.idx <- unique(th.idx[th.idx > 0L])

  # fix Delta's...
  mm.in.group <- 1:lavmodel@nmat[group] + cumsum(c(0, lavmodel@nmat))[group]
  MLIST <- GLIST[mm.in.group]

  DD <- list()
  nvar <- lavmodel@nvar
  nfac <- ncol(MLIST$lambda) - length(lv.dummy.idx)

  # DD$theta
  theta.idx <- lav_matrix_diagh_idx(nvar)
  DD$theta <- Delta.theta[theta.idx, , drop = FALSE]
  if (length(ov.dummy.idx) > 0L) {
    psi.idx <- lav_matrix_diagh_idx(ncol(MLIST$psi))[lv.dummy.idx]
    DD$theta[ov.dummy.idx, ] <- Delta.psi[psi.idx, , drop = FALSE]
  }
  # num only? FIXME or just all of them?
  DD$theta <- DD$theta[num.idx, , drop = FALSE]

  # DD$nu
  DD$nu <- Delta.nu
  if (length(ov.dummy.idx) > 0L) {
    DD$nu[ov.dummy.idx, ] <- Delta.alpha[lv.dummy.idx, ]
  }
  DD$nu <- DD$nu[num.idx, , drop = FALSE] # needed?

  # DD$lambda
  nr <- nvar
  nc <- nfac
  lambda.idx <- nr * ((1:nc) - 1L) + rep(1:nvar, each = nc)
  DD$lambda <- Delta.lambda[lambda.idx, , drop = FALSE]
  if (length(ov.dummy.idx) > 0L) {
    nr <- nrow(MLIST$beta)
    nc <- nfac # only the first 1:nfac columns
    # beta.idx <- rep(nr*((1:nc) - 1L), each=length(lv.dummy.idx)) + rep(lv.dummy.idx, times=nc) ## FIXME
    beta.idx <- rep(nr * ((1:nc) - 1L), times = length(lv.dummy.idx)) + rep(lv.dummy.idx, each = nc)

    # l.idx <- inr*((1:nc) - 1L) + rep(ov.dummy.idx, each=nc) ## FIXME
    # l.idx <- rep(nr*((1:nc) - 1L), each=length(ov.dummy.idx)) + rep(ov.dummy.idx, times=nc)
    l.idx <- rep(nr * ((1:nc) - 1L), times = length(ov.dummy.idx)) + rep(ov.dummy.idx, each = nc)
    DD$lambda[match(l.idx, lambda.idx), ] <- Delta.beta[beta.idx, , drop = FALSE]
  }

  # DD$KAPPA
  DD$kappa <- Delta.gamma
  if (length(ov.dummy.idx) > 0L) {
    nr <- nrow(MLIST$gamma)
    nc <- ncol(MLIST$gamma)
    kappa.idx <- nr * ((1:nc) - 1L) + rep(lv.dummy.idx, each = nc)
    DD$kappa <- DD$kappa[kappa.idx, , drop = FALSE]
  }

  # DD$GAMMA
  if (!is.null(MLIST$gamma)) {
    nr <- nrow(MLIST$gamma)
    nc <- ncol(MLIST$gamma)
    lv.idx <- 1:nfac
    # MUST BE ROWWISE!
    gamma.idx <- rep(nr * ((1:nc) - 1L), times = length(lv.idx)) + rep(lv.idx, each = nc)
    DD$gamma <- Delta.gamma[gamma.idx, , drop = FALSE]
  }

  # DD$BETA
  if (!is.null(MLIST$beta)) {
    nr <- nc <- nrow(MLIST$beta)
    lv.idx <- 1:nfac
    # MUST BE ROWWISE!
    beta.idx <- rep(nr * ((1:nfac) - 1L), times = nfac) + rep(lv.idx, each = nfac)
    DD$beta <- Delta.beta[beta.idx, , drop = FALSE]
  }

  ## DD$psi
  DD$psi <- Delta.psi
  if (length(lv.dummy.idx) > 0L) {
    nr <- nc <- nrow(MLIST$psi)
    lv.idx <- 1:nfac
    # MUST BE ROWWISE!
    psi.idx <- rep(nr * ((1:nfac) - 1L), times = nfac) + rep(lv.idx, each = nfac)

    DD$psi <- DD$psi[psi.idx, , drop = FALSE]
  }

  ## DD$tau
  if (!is.null(MLIST$tau)) {
    DD$tau <- Delta.tau
  }

  DD
}
