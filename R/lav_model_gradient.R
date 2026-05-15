# model gradient

lav_model_gradient <- function(lavmodel = NULL,
                               glist = NULL,
                               lavsamplestats = NULL,
                               lavdata = NULL,
                               lavcache = NULL,
                               type = "free",
                               group_weight = TRUE,
                               delta = NULL,
                               m_el_idx = NULL,
                               x_el_idx = NULL,
                               ceq_simple = FALSE,
                               implied = NULL) {
  nmat <- lavmodel@nmat
  mm_idx <- lav_model_get_mm_idx(lavmodel)
  estimator <- lavmodel@estimator
  representation <- lavmodel@representation
  meanstructure <- lavmodel@meanstructure
  categorical <- lavmodel@categorical
  group_w_free <- lavmodel@group.w.free
  # fixed_x <- lavmodel@fixed.x
  conditional_x <- lavmodel@conditional.x
  num_idx <- lavmodel@num.idx
  th_idx <- lavmodel@th.idx
  nx_free <- lavmodel@nx.free
  estimator_args <- lavmodel@estimator.args

  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  if (estimator == "REML") lav_msg_warn(gettext(
    "analytical gradient not implemented; use numerical approximation"))

  # group.weight
  # FIXME --> block.weight
  if (group_weight) {
    if (estimator %in% c("ML", "PML", "FML", "MML", "REML", "NTRLS", "catML")) {
      group_w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)
    } else if (estimator == "DLS") {
      if (estimator_args$dls.FtimesNminus1) {
        group_w <- ((unlist(lavsamplestats@nobs) - 1) / lavsamplestats@ntotal)
      } else {
        group_w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)
      }
    } else {
      # FIXME: double check!
      group_w <- ((unlist(lavsamplestats@nobs) - 1) / lavsamplestats@ntotal)
    }
  } else {
    group_w <- rep(1.0, lavmodel@nblocks)
  }

  # do we need WLS.est?
  if (estimator %in% c("WLS", "DWLS", "ULS", "GLS", "NTRLS", "DLS")) {
    # always compute WLS.est
    wls_est <- lav_model_wls_est(lavmodel = lavmodel, glist = glist) # ,
    # cov.x = lavsamplestats@cov.x)
  }

  if (estimator %in% c("ML", "PML", "FML", "REML", "NTRLS", "catML")) {
    # compute moments for all groups
    # if(conditional.x) {
    #    Sigma.hat <- lav_model_cond2joint_sigma(lavmodel = lavmodel,
    #                     GLIST = GLIST,
    #                     extra = (estimator %in% c("ML", "REML","NTRLS")))
    # } else {
    implied_fast <- lav_model_implied_fast_state(
      lavmodel = lavmodel, glist = glist,
      implied = implied,
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

    if (meanstructure) {
      # if(conditional.x) {
      #    Mu.hat <- lav_model_mu(lavmodel = lavmodel, GLIST = GLIST)
      # } else {
      mu_hat <- implied_fast$mu
      # }
    }

    if (categorical) {
      th <- implied_fast$th
    }

    if (conditional_x) {
      pi0 <- implied_fast$pi
    } else if (estimator == "PML") {
      pi0 <- vector("list", length = lavmodel@nblocks)
    }

    # if (group_w_free) {
    #   gw <- lav_model_gw(lavmodel = lavmodel, glist = glist)
    # }
  } else if (estimator == "DLS" && estimator_args$dls.GammaNT == "model") {
    implied_fast <- lav_model_implied_fast_state(
      lavmodel = lavmodel, glist = glist,
      implied = implied,
      need_sigma = TRUE,
      need_mu = TRUE,
      extra = FALSE
    )
    sigma_hat <- implied_fast$sigma
    mu_hat <- implied_fast$mu
  } else if (estimator == "MML") {
    th <- lav_model_th(lavmodel = lavmodel, glist = glist)
    mm_theta <- lav_model_theta(lavmodel = lavmodel, glist = glist)
    # gw <- lav_model_gw(lavmodel = lavmodel, glist = glist)
  }

  # four approaches (FIXME!!!! merge this!)
  # - ML approach: using Omega (and Omega.mu)
  #    Omega = 'POST' = Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)
  #   (still 2x faster than Delta method)
  # - WLS/DWLS/GLS: using Delta + WLS.V; support for fixed.x, conditional.x
  # - (ML)/NTRLS: using Delta, no support for fixed.x, conditional.x
  # - PML/FML/MML: custom

  # composites?
  composites_flag <- lavmodel@composites

  # 1. ML approach
  if ((estimator == "ML" || estimator == "REML" || estimator == "catML") &&
    lavdata@nlevels == 1L && !composites_flag &&
    !lavmodel@conditional.x) {
  correlation <- lavmodel@correlation
    if (meanstructure) {
      omega <- lav_model_omega(
        sigma_hat = sigma_hat, mu_hat = mu_hat,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = TRUE,
        conditional_x = conditional_x,
    correlation = correlation
      )
      omega_mu <- attr(omega, "mu")
    } else {
      omega <- lav_model_omega(
        sigma_hat = sigma_hat, mu_hat = NULL,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = FALSE,
        conditional_x = conditional_x,
    correlation = correlation
      )
      omega_mu <- vector("list", length = lavmodel@nblocks)
    }

    # compute DX (for all elements in every model matrix)
    dx_1 <- vector("list", length = length(glist))

    for (g in 1:lavmodel@nblocks) {
      # which mm belong to group g?
      mm_in_group <- mm_idx[[g]]
      mm_names <- names(glist[mm_in_group])

      if (representation == "LISREL") {
        dx_group <- lav_lisrel_df_dmlist(
          glist[mm_in_group],
          omega[[g]],
          omega_mu[[g]]
        )

        # FIXME!!!
        # add empty gamma
        if (lavmodel@conditional.x) {
          dx_group$gamma <- lavmodel@GLIST$gamma
        }

        # only save what we need
        dx_1[mm_in_group] <- dx_group[mm_names]
      } else if (representation == "RAM") {
        dx_group <- lav_ram_df(
          glist[mm_in_group],
          omega[[g]],
          omega_mu[[g]]
        )
        # only save what we need
        dx_1[mm_in_group] <- dx_group[mm_names]
      } else {
        lav_msg_stop(gettext(
          "only LISREL and RAM representation has been implemented for now"))
      }

      # weight by group
      if (lavmodel@nblocks > 1L) {
        for (mm in mm_in_group) {
          dx_1[[mm]] <- group_w[g] * dx_1[[mm]]
        }
      }
    }

    # extract free parameters

    if (type == "free") {
      if (lavmodel@ceq.simple.only) { # new in 0.6-11
        dx <- numeric(lavmodel@nx.unco)
        for (g in 1:lavmodel@nblocks) {
          mm_in_group <- mm_idx[[g]]
          for (mm in mm_in_group) {
            m_free_idx <- lavmodel@m.free.idx[[mm]]
            x_unco_idx <- lavmodel@x.unco.idx[[mm]]
            dx[x_unco_idx] <- dx_1[[mm]][m_free_idx]
          }
        }
        if (ceq_simple) {
          dx <- drop(crossprod(lavmodel@ceq.simple.K, dx))
        }
      } else {
        dx <- numeric(nx_free)
        for (g in 1:lavmodel@nblocks) {
          mm_in_group <- mm_idx[[g]]
          for (mm in mm_in_group) {
            m_free_idx <- lavmodel@m.free.idx[[mm]]
            x_free_idx <- lavmodel@x.free.idx[[mm]]
            dx[x_free_idx] <- dx_1[[mm]][m_free_idx]
          }
        }
      }
    } else {
      dx <- dx_1
      # handle equality constraints
      ### FIXME!!!! TODO!!!!
    }
  } else # ML

  # 1a. ML approach for composites: use Delta directly. The existing
  # per-matrix lav_lisrel_df_dmlist() path does not handle composites
  # (it ignores wmat / WTW.inv). We reuse the standard ML "POST" trick:
  #   dF/dphi = crossprod(Delta, POST)  where
  #     POST = c(-2 * Omega.mu, -D' vec(Omega))
  # This works for any model where lav_model_delta() can produce an
  # analytical Jacobian -- including composites since 0.6-23.
  if ((estimator == "ML" || estimator == "REML" || estimator == "catML") &&
      lavdata@nlevels == 1L && composites_flag &&
      !lavmodel@conditional.x) {
    correlation <- lavmodel@correlation
    if (meanstructure) {
      omega <- lav_model_omega(
        sigma_hat = sigma_hat, mu_hat = mu_hat,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = TRUE,
        conditional_x = conditional_x,
        correlation = correlation
      )
      omega_mu <- attr(omega, "mu")
    } else {
      omega <- lav_model_omega(
        sigma_hat = sigma_hat, mu_hat = NULL,
        lavsamplestats = lavsamplestats,
        estimator = estimator,
        meanstructure = FALSE,
        conditional_x = conditional_x,
        correlation = correlation
      )
      omega_mu <- vector("list", length = lavmodel@nblocks)
    }

    delta <- lav_model_delta(lavmodel = lavmodel, glist = glist)

    if (lavmodel@ceq.simple.only) {
      dx <- numeric(lavmodel@nx.unco)
    } else {
      dx <- numeric(nx_free)
    }

    # lavaan's ML objective is 0.5 * F_ML (lav_model_objective.R: group.fx
     # <- 0.5 * group.fx for ML/REML/NTRLS/catML), so dF/dphi has a 0.5
     # factor relative to the textbook ML loss.
    for (g in 1:lavmodel@nblocks) {
      post_sigma_1 <- -0.5 * lav_matrix_duplication_pre(
        matrix(omega[[g]], ncol = 1L))
      if (meanstructure) {
        post_mu <- -1 * as.numeric(omega_mu[[g]])
        post <- c(post_mu, post_sigma_1)
      } else {
        post <- post_sigma_1
      }
      # Delta rows: [gw (if group.w.free) | mu (if meanstructure) | vech(Sigma)]
      # group weight is overwritten below via the explicit Poisson term, so
      # pad POST with a leading zero when needed (gradient there is set later)
      if (lavmodel@group.w.free) {
        post <- c(0, post)
      }
      group_dx <- as.numeric(crossprod(delta[[g]], post))
      dx <- dx + group_w[g] * group_dx
    }

    if (lavmodel@ceq.simple.only && ceq_simple) {
      dx <- drop(crossprod(lavmodel@ceq.simple.K, dx))
    }
  } else # ML + composites

  # 2. using Delta - *LS family
  if (estimator %in% c("WLS", "DWLS", "ULS", "GLS", "NTGLS", "DLS")) {
    if (type != "free") {
      if (is.null(delta)) {
        lav_msg_fixme("Delta should be given if type != free")
      }
      # stop("FIXME: WLS gradient with type != free needs fixing!")
    } else {
      delta <- lav_model_delta(lavmodel = lavmodel, glist = glist)
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
        diff <- lavsamplestats@WLS.obs[[g]] - wls_est[[g]]

        # full weight matrix
        if (estimator == "GLS" || estimator == "WLS") {
          wls_v <- lavsamplestats@WLS.V[[g]]
          group_dx <- -1 * crossprod(
            delta[[g]],
            crossprod(wls_v, diff)
          )
        } else if (estimator == "DLS") {
          if (estimator_args$dls.GammaNT == "sample") {
            wls_v <- lavsamplestats@WLS.V[[g]] # for now
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
          group_dx <- -1 * crossprod(
            delta[[g]],
            crossprod(wls_v, diff)
          )
        } else if (estimator == "NTRLS") {
          stopifnot(!conditional_x)
          # WLS.V <- lav_samplestats_gamma_inverse_nt(
          #         m_icov = attr(Sigma.hat[[g]],"inv")[,,drop=FALSE],
          #         m_cov            = Sigma.hat[[g]][,,drop=FALSE],
          #         m_mean           = Mu.hat[[g]],
          #         x_idx          = lavsamplestats@x.idx[[g]],
          #         fixed_x        = fixed.x,
          #         conditional_x  = conditional.x,
          #         meanstructure  = meanstructure,
          #         slopestructure = conditional.x)

          s <- lavsamplestats@cov[[g]]
          sigma_1 <- sigma_hat[[g]]
          sigma_inv <- attr(sigma_1, "inv")
          nvar <- NROW(sigma_1)

          if (meanstructure) {
            mean_1 <- lavsamplestats@mean[[g]]
            mu <- mu_hat[[g]]
            post_sigma_1 <- lav_matrix_duplication_pre(
              matrix(
                (sigma_inv %*% (s - sigma_1) %*% t(sigma_inv)) %*%
                  (diag(nvar) + (s - sigma_1) %*% sigma_inv) +
                  (sigma_inv %*% tcrossprod(mean_1 - mu) %*% sigma_inv),
                ncol = 1
              )
            )
            post_mu <- as.numeric(2 * sigma_inv %*% (mean_1 - mu))
            post <- c(post_mu, post_sigma_1)
          } else {
            post <- lav_matrix_duplication_pre(
              matrix((sigma_inv %*% (s - sigma_1) %*% t(sigma_inv)) %*%
                (diag(nvar) + (s - sigma_1) %*% sigma_inv), ncol = 1)
            )
          }

          group_dx <- as.numeric(-1 * crossprod(delta[[g]], post))
        }
      } else if (estimator == "DWLS" || estimator == "ULS") {
        # diagonal weight matrix
        diff <- lavsamplestats@WLS.obs[[g]] - wls_est[[g]]
        group_dx <- -1 * crossprod(
          delta[[g]],
          lavsamplestats@WLS.VD[[g]] * diff
        )
      }

      group_dx <- group_w[g] * group_dx
      if (g == 1) {
        dx <- group_dx
      } else {
        dx <- dx + group_dx
      }
    } # g

    if (type == "free") {
      # nothing to do
    } else {
      # make a GLIST
      dx <- lav_model_x2glist(
        lavmodel = lavmodel, x = dx,
        type = "custom", set_delta = FALSE,
        m_el_idx = m_el_idx,
        x_el_idx = x_el_idx
      )
    }

  # WLS
  # ML + conditional.x
  } else if (estimator %in% c("ML", "catML") && lavmodel@conditional.x &&
    lavdata@nlevels == 1L) {
    if (type != "free") {
      if (is.null(delta)) {
        lav_msg_fixme("Delta should be given if type != free")
      }
      # stop("FIXME: WLS gradient with type != free needs fixing!")
    } else {
      delta <- lav_model_delta(lavmodel = lavmodel, glist = glist)
    }

    for (g in 1:lavmodel@nblocks) {
      # augmented mean.x + cov.x matrix
      mean_x <- lavsamplestats@mean.x[[g]]
      cov_x <- lavsamplestats@cov.x[[g]]
      c3 <- rbind(
        c(1, mean_x),
        cbind(mean_x, cov_x + tcrossprod(mean_x))
      )

      sigma_1 <- sigma_hat[[g]]
      mu_g <- mu_hat[[g]]
      pi_g <- pi0[[g]]
      sigma_inv <- attr(sigma_1, "inv")
      nvar <- NROW(sigma_1)
      s <- lavsamplestats@res.cov[[g]]

      # beta
      obs <- t(cbind(
        lavsamplestats@res.int[[g]],
        lavsamplestats@res.slopes[[g]]
      ))
      est <- t(cbind(mu_g, pi_g))
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
      d_beta <- c3 %*% (obs - est) %*% sigma_inv
      # NOTE: the vecr here, unlike lav_mvreg_dlogl_beta
      #       this is because DELTA has used vec(t(BETA)),
      #       instead of vec(BETA)
      # POST.beta <- 2 * lav_matrix_vecr(d.BETA)
      # NOT any longer, since 0.6-1!!!
      post_beta <- 2 * lav_matrix_vec(d_beta)

      # POST.sigma1 <- lav_matrix_duplication_pre(
      #        (Sigma.inv %x% Sigma.inv)  %*% t(AB)  %*% (t(a) %x% b) )

      # Sigma
      # POST.sigma2 <- lav_matrix_duplication_pre(
      #                 matrix( lav_matrix_vec(
      #          Sigma.inv %*% (S - Sigma) %*% t(Sigma.inv)), ncol = 1L))
      w_tilde <- s + t(obs - est) %*% c3 %*% (obs - est)
      d_sigma <- (sigma_inv - sigma_inv %*% w_tilde %*% sigma_inv)
      d_vech_sigma <- as.numeric(lav_matrix_duplication_pre(
        as.matrix(lav_matrix_vec(d_sigma))
      ))
      post_sigma <- -1 * d_vech_sigma

      # POST <- c(POST.beta, POST.sigma1 + POST.sigma2)
      post <- c(post_beta, post_sigma)

      group_dx <- as.numeric(-1 * crossprod(delta[[g]], post))

      # because we still use obj/2, we need to divide by 2!
      group_dx <- group_dx / 2 # fixed in 0.6-1

      group_dx <- group_w[g] * group_dx
      if (g == 1) {
        dx <- group_dx
      } else {
        dx <- dx + group_dx
      }
    } # g

    if (type == "free") {
      # nothing to do
    } else {
      # make a GLIST
      dx <- lav_model_x2glist(
        lavmodel = lavmodel, x = dx,
        type = "custom", set_delta = FALSE,
        m_el_idx = m_el_idx,
        x_el_idx = x_el_idx
      )
    }

  # ML + conditional.x
  } else if (estimator == "ML" && lavdata@nlevels > 1L) {
    if (type != "free") {
      lav_msg_fixme("type != free in lav_model_gradient for
                    estimator ML for nlevels > 1")
    } else {
      delta <- lav_model_delta(lavmodel = lavmodel, glist = glist)
    }

    # for each upper-level group....
    for (g in 1:lavmodel@ngroups) {
      if (!lavsamplestats@missing.flag) { # complete data
        if (lavmodel@conditional.x) {
          dx_1 <- lav_mvreg_cluster_dlogl_2l_samplestats(
            ylp = lavsamplestats@YLp[[g]],
            lp = lavdata@Lp[[g]],
            res_sigma_w = sigma_hat[[(g - 1) * 2 + 1]],
            res_int_w = mu_hat[[(g - 1) * 2 + 1]],
            res_pi_w = pi0[[(g - 1) * 2 + 1]],
            res_sigma_b = sigma_hat[[(g - 1) * 2 + 2]],
            res_int_b = mu_hat[[(g - 1) * 2 + 2]],
            res_pi_b = pi0[[(g - 1) * 2 + 2]],
            sinv_method = "eigen"
          )
        } else {
          dx_1 <- lav_mvnorm_cluster_dlogl_2l_samplestats(
            ylp = lavsamplestats@YLp[[g]],
            lp = lavdata@Lp[[g]],
            mu_w = mu_hat[[(g - 1) * 2 + 1]],
            sigma_w = sigma_hat[[(g - 1) * 2 + 1]],
            mu_b = mu_hat[[(g - 1) * 2 + 2]],
            sigma_b = sigma_hat[[(g - 1) * 2 + 2]],
            sinv_method = "eigen"
          )
        }
      } else {
        # missing data
        if (lavmodel@conditional.x) {
          lav_msg_stop(gettext("gradient for twolevel + conditional.x + fiml
                  is not ready; use optim.gradient = \"numerical\""))
        } else {
          dx_1 <- lav_mvnorm_cluster_missing_dlogl_2l_samplestats(
            y1 = lavdata@X[[g]],
            y2 = lavsamplestats@YLp[[g]][[2]]$Y2,
            lp = lavdata@Lp[[g]],
            mp = lavdata@Mp[[g]],
            mu_w = mu_hat[[(g - 1) * 2 + 1]],
            sigma_w = sigma_hat[[(g - 1) * 2 + 1]],
            mu_b = mu_hat[[(g - 1) * 2 + 2]],
            sigma_b = sigma_hat[[(g - 1) * 2 + 2]],
            sinv_method = "eigen"
          )
        }
      }

      group_dx <- as.numeric(dx_1 %*% delta[[g]])

      # group weights (if any)
      group_dx <- group_w[g] * group_dx
      if (g == 1) {
        dx <- group_dx
      } else {
        dx <- dx + group_dx
      }
    } # g

    # divide by 2 * N
    dx <- dx / (2 * lavsamplestats@ntotal)

  # cat("dx1 (numerical) = \n"); print( zapsmall(dx1) )
  # cat("dx  (analytic)  = \n"); print( zapsmall(dx ) )
  # ML + two-level
  } else if (estimator == "PML" || estimator == "FML" ||
    estimator == "MML") {
    if (type != "free") {
      lav_msg_fixme("type != free in lav_model_gradient for estimator PML")
    } else {
      delta <- lav_model_delta(lavmodel = lavmodel, glist = glist)
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
        } else if (conditional_x) {
          d1 <- lav_pml_dploglik_dimplied(
            sigma_hat = sigma_hat[[g]],
            mu_hat = mu_hat[[g]],
            th = th[[g]],
            th_idx = th_idx[[g]],
            num_idx = num_idx[[g]],
            x = lavdata@X[[g]],
            lavcache = lavcache[[g]],
            exo = lavdata@eXo[[g]],
            wt = lavdata@weights[[g]],
            pi0 = pi0[[g]],
            missing = lavdata@missing
          )
        } else {
          d1 <- lav_pml_dploglik_dimplied(
            sigma_hat = sigma_hat[[g]],
            mu_hat = mu_hat[[g]],
            th = th[[g]],
            th_idx = th_idx[[g]],
            num_idx = num_idx[[g]],
            x = lavdata@X[[g]],
            lavcache = lavcache[[g]],
            exo = NULL,
            wt = lavdata@weights[[g]],
            pi0 = NULL,
            missing = lavdata@missing
          )
        } # not conditional.x

        # chain rule (fmin)
        group_dx <-
          as.numeric(t(d1) %*% delta[[g]])
    # PML
    } else if (estimator == "FML") {
        d1 <- lav_pml_fml_dploglik_dimplied(
          sigma_hat = sigma_hat[[g]],
          th = th[[g]],
          th_idx = th_idx[[g]],
          num_idx = num_idx[[g]],
          x = lavdata@X[[g]],
          lavcache = lavcache[[g]]
        )

        # chain rule (fmin)
        group_dx <-
          as.numeric(t(d1) %*% delta[[g]]) / lavsamplestats@nobs[[g]]
      } else if (estimator == "MML") {
        group_dx <-
          lav_model_gradient_mml(
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
      }

      # group weights (if any)
      group_dx <- group_w[g] * group_dx
      if (g == 1) {
        dx <- group_dx
      } else {
        dx <- dx + group_dx
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
    est_freq <- exp(unlist(lav_model_gw(lavmodel = lavmodel, glist = glist)))
    obs_freq <- unlist(lavsamplestats@group.w) * lavsamplestats@ntotal
    dx_gw <- -(obs_freq - est_freq)
    # divide by N (to be consistent with the rest of lavaan)
    dx_gw <- dx_gw / lavsamplestats@ntotal

    # remove last element (fixed LAST group to zero)
    # dx.GW <- dx.GW[-length(dx.GW)]

    # fill in in dx
    gw_mat_idx <- which(names(lavmodel@GLIST) == "gw")
    gw_x_idx <- unlist(lavmodel@x.free.idx[gw_mat_idx])
    dx[gw_x_idx] <- dx_gw
  }

  # dx is 1xnpar matrix of LIST (type != "free")
  if (is.matrix(dx)) {
    dx <- as.numeric(dx)
  }

  dx
}

# for testing purposes only
lav_model_delta_numerical <- function(lavmodel = NULL, glist = NULL, g = 1L) {

  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  compute_moments <- function(x) {
    glist <- lav_model_x2glist(lavmodel = lavmodel, x = x, type = "free")
    sigma_hat <- lav_model_sigma(lavmodel = lavmodel, glist = glist)
    s_vec <- lav_matrix_vech(sigma_hat[[g]])
    if (lavmodel@meanstructure) {
      mu_hat <- lav_model_mu(lavmodel = lavmodel, glist = glist)
      out <- c(mu_hat[[g]], s_vec)
    } else {
      out <- s_vec
    }
    out
  }

  x <- lav_model_get_parameters(lavmodel = lavmodel,
                                GLIST = glist, type = "free")
  delta <- lav_func_jacobian_complex(func = compute_moments, x = x)

  delta
}


lav_model_ddelta_dx <- function(lavmodel = NULL, glist = NULL,
   target = "lambda", ceq_simple = FALSE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  representation <- lavmodel@representation
  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  th_idx <- lavmodel@th.idx

  # number of columns in DELTA + m.el.idx/x.el.idx
  type <- "free"
  # if(type == "free") {
  if (lavmodel@ceq.simple.only) {
    n_col <- lavmodel@nx.unco
  } else {
    n_col <- lavmodel@nx.free
  }
  m_el_idx <- x_el_idx <- vector("list", length = length(glist))
  mm_idx <- lav_model_get_mm_idx(lavmodel)
  for (mm in seq_along(glist)) {
    m_el_idx[[mm]] <- lavmodel@m.free.idx[[mm]]
    if (lavmodel@ceq.simple.only) {
      x_el_idx[[mm]] <- lavmodel@x.unco.idx[[mm]]
    } else {
      x_el_idx[[mm]] <- lavmodel@x.free.idx[[mm]]
    }
    # handle symmetric matrices
    if (lavmodel@isSymmetric[mm]) {
      # since we use 'x.free.idx', only symmetric elements
      # are duplicated (not the equal ones, only in x.free.free)
      dix <- duplicated(x_el_idx[[mm]])
      if (any(dix)) {
        m_el_idx[[mm]] <- m_el_idx[[mm]][!dix]
        x_el_idx[[mm]] <- x_el_idx[[mm]][!dix]
      }
    }
  }
  # } else {
  #    n_col <- sum(unlist(lapply(x.el.idx, function(x) length(unique(x)))))
  # }

  # compute Delta per group
  delta <- vector("list", length = nblocks)
  for (g in 1:nblocks) {
    mm_in_group <- mm_idx[[g]]
    delta_group <- NULL
    for (mm in mm_in_group) {
      mname <- names(lavmodel@GLIST)[mm]

      # skip empty ones
      if (!length(m_el_idx[[mm]])) next

      # get Delta columns for this model matrix
      if (representation == "LISREL") {
        if (target == "lambda") {
         mm_delta <- lav_lisrel_dlambda_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "th") {
          mm_delta <- lav_lisrel_dth_dx(
             mlist = glist[mm_in_group], m = mname, th_idx = th_idx[[g]],
            idx = m_el_idx[[mm]],
            delta = TRUE
          )
        } else if (target == "mu") {
          mm_delta <- lav_lisrel_dmu_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "nu") {
          mm_delta <- lav_lisrel_dnu_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "tau") {
          mm_delta <- lav_lisrel_dtau_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "theta") {
          mm_delta <- lav_lisrel_dtheta_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "gamma") {
          mm_delta <- lav_lisrel_dgamma_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "beta") {
          mm_delta <- lav_lisrel_dbeta_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "alpha") {
          mm_delta <- lav_lisrel_dalpha_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "psi") {
          mm_delta <- lav_lisrel_dpsi_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]]
          )
        } else if (target == "sigma") {
          mm_delta <- lav_lisrel_dsigma_dx(
            mlist = glist[mm_in_group],
            m = mname,
            idx = m_el_idx[[mm]],
            delta = TRUE
          )
        } else {
          lav_msg_stop(gettextf("target %s not implemented yet", target))
        }

        # initialize?
        if (is.null(delta_group)) {
          delta_group <- matrix(0, nrow = nrow(mm_delta), ncol = n_col)
        }
        delta_group[, x_el_idx[[mm]]] <- mm_delta
      }
    } # mm

    if (type == "free" && ceq_simple && lavmodel@ceq.simple.only) {
      delta_group <- delta_group %*% lavmodel@ceq.simple.K
    }

    delta[[g]] <- delta_group
  } # g

  delta
}


# Compute the Delta matrix: Jacobian of the model-implied moments wrt the
# free parameters, per block.
#
# Supported:
#   - LISREL representation, full surface: continuous + categorical +
#     correlation, conditional.x on/off, parameterization "delta" or
#     "theta", composites, group.w.free. Both the standard "free" indexing
#     and caller-supplied m.el.idx./x.el.idx. (e.g. lav_pml_test_plrt(),
#     modindices()) are handled.
#   - RAM representation, classical mean- and covariance structures only.
#     RAM with categorical / correlation / conditional.x / composites is
#     not supported and produces an error.
#
# This is the housekeeping wrapper: it loops over blocks, extracts the
# block-specific MLIST and free-element indices, and dispatches the actual
# computation to the representation-specific worker.
lav_model_delta <- function(lavmodel = NULL, glist = NULL,
                            m_el_idx = NULL, x_el_idx = NULL,
                            ceq_simple = FALSE) {

  representation <- lavmodel@representation
  if (representation == "RAM") {
    if (lavmodel@conditional.x || lavmodel@categorical ||
        lavmodel@correlation   || lavmodel@composites) {
      lav_msg_stop(gettext(
        "lav_model_delta(): RAM representation does not support
         conditional.x, categorical, correlation, or composites."))
    }
  } else if (representation != "LISREL") {
    lav_msg_stop(gettextf(
      "lav_model_delta(): representation %s not supported.",
      representation))
  }

  # state or final?
  if (is.null(glist)) {
    glist_1 <- lavmodel@GLIST
  } else {
    glist_1 <- glist
  }

  # which (model element, x column) lists do we use?
  #  - 'free' (default): use lavmodel@m.free.idx + x.free.idx (or x.unco.idx
  #     when ceq.simple.only); n_col = nx.free (or nx.unco). The worker writes
  #     Delta with nx.unco columns when ceq.simple.only; the caller maps to
  #     nx.free via ceq.simple.K.
  #  - 'custom': caller passes m.el.idx./x.el.idx. directly (e.g. modindices
  #     after lav_object_extended()); n_col = max column index.
  use_custom <- !is.null(m_el_idx) || !is.null(x_el_idx)
  if (use_custom) {
    m_idx <- m_el_idx
    x_idx <- x_el_idx
    nx    <- max(unlist(x_idx), 0L)
  } else if (lavmodel@ceq.simple.only) {
    m_idx <- lavmodel@m.free.idx
    x_idx <- lavmodel@x.unco.idx
    nx    <- lavmodel@nx.unco
  } else {
    m_idx <- lavmodel@m.free.idx
    x_idx <- lavmodel@x.free.idx
    nx    <- lavmodel@nx.free
  }

  delta <- vector("list", length = lavmodel@nblocks)
  mm_idx <- lav_model_get_mm_idx(lavmodel)
  for (b in seq_len(lavmodel@nblocks)) {
    mm_in_group <- mm_idx[[b]]

    if (representation == "LISREL") {
      delta[[b]] <- lav_lisrel_dimplied_dx(
        mlist            = glist_1[mm_in_group],
        m_free_idx       = m_idx[mm_in_group],
        x_free_idx       = x_idx[mm_in_group],
        nx_free          = nx,
        meanstructure    = lavmodel@meanstructure,
        categorical      = lavmodel@categorical,
        correlation      = lavmodel@correlation,
        conditional_x    = lavmodel@conditional.x,
        num_idx          = lavmodel@num.idx[[b]],
        th_idx           = lavmodel@th.idx[[b]],
        group_w_free     = lavmodel@group.w.free,
        parameterization = lavmodel@parameterization
      )
    } else { # RAM
      delta[[b]] <- lav_ram_dimplied_dx(
        mlist         = glist_1[mm_in_group],
        m_free_idx    = m_idx[mm_in_group],
        x_free_idx    = x_idx[mm_in_group],
        nx_free       = nx,
        meanstructure = lavmodel@meanstructure,
        group_w_free  = lavmodel@group.w.free
      )
    }
  }

  # if multilevel, rbind levels within group
  if (lavmodel@multilevel) {
    delta_tmp <- vector("list", length = lavmodel@ngroups)
    for (g in 1:lavmodel@ngroups) {
      delta_tmp[[g]] <- rbind(
        delta[[(g - 1) * 2 + 1]],
        delta[[(g - 1) * 2 + 2]]
      )
    }
    delta <- delta_tmp
  }

  delta
}



lav_model_omega <- function(sigma_hat = NULL, mu_hat = NULL,
                         lavsamplestats = NULL, estimator = "ML",
                         meanstructure = FALSE, conditional_x = FALSE,
             correlation = FALSE) {
  # nblocks
  nblocks <- length(sigma_hat)

  omega <- vector("list", length = nblocks)
  omega_mu <- vector("list", length = nblocks)

  for (g in 1:nblocks) {
    # ML
    if (estimator %in% c("ML", "REML", "catML")) {
      if (attr(sigma_hat[[g]], "po") == FALSE) {
        # FIXME: WHAT IS THE BEST THING TO DO HERE??
        # CURRENTLY: stop
        lav_msg_warn(gettext(
          "lav_model_gradient: Sigma.hat is not positive definite\n"))
        sigma_hat_inv <- MASS::ginv(sigma_hat[[g]])
      } else {
        sigma_hat_inv <- attr(sigma_hat[[g]], "inv")
      }

      if (!lavsamplestats@missing.flag) { # complete data
        if (meanstructure) {
          if (conditional_x) {
            diff <- lavsamplestats@res.int[[g]] - mu_hat[[g]]
            w_tilde <- lavsamplestats@res.cov[[g]] + tcrossprod(diff)
          } else {
            diff <- lavsamplestats@mean[[g]] - mu_hat[[g]]
            w_tilde <- lavsamplestats@cov[[g]] + tcrossprod(diff)
          }
          # Browne 1995 eq 4.55
          omega_mu[[g]] <- t(t(diff) %*% sigma_hat_inv)
          omega[[g]] <-
            (sigma_hat_inv %*% (w_tilde - sigma_hat[[g]]) %*%
              sigma_hat_inv)
        } else {
          if (conditional_x) {
            w_tilde <- lavsamplestats@res.cov[[g]]
          } else {
            w_tilde <- lavsamplestats@cov[[g]]
          }
          omega[[g]] <-
            (sigma_hat_inv %*% (w_tilde - sigma_hat[[g]]) %*%
              sigma_hat_inv)
        }
      } else { # missing data
        m <- lavsamplestats@missing[[g]]

        nvar <- ncol(lavsamplestats@cov[[g]])
        omega_1 <- matrix(0, nvar, nvar)
        omega_mu_1 <- matrix(0, nvar, 1)

        for (p in seq_along(m)) {
          sx <- m[[p]][["SY"]]
          mx <- m[[p]][["MY"]]
          nobs <- m[[p]][["freq"]]
          var_idx <- m[[p]][["var.idx"]]

          sigma_inv <- try(chol2inv(chol(sigma_hat[[g]][var_idx, var_idx])),
                          silent = TRUE)
          if (inherits(sigma_inv, "try-error")) {
            sigma_inv <- MASS::ginv(sigma_hat[[g]][var_idx, var_idx])
          }
          mu <- mu_hat[[g]][var_idx]
          w_tilde <- sx + tcrossprod(mx - mu)

          omega_mu_1[var_idx, 1] <-
            (omega_mu_1[var_idx, 1] + nobs / lavsamplestats@ntotal *
              t(t(mx - mu) %*% sigma_inv))

          omega_1[var_idx, var_idx] <-
            (omega_1[var_idx, var_idx] + nobs / lavsamplestats@ntotal *
              (sigma_inv %*%
                (w_tilde - sigma_hat[[g]][var_idx, var_idx]) %*%
                sigma_inv))
        }
        omega_mu[[g]] <- omega_mu_1
        omega[[g]] <- omega_1
      } # missing

      # GLS
    } else if (estimator == "GLS") {
      w_inv <- lavsamplestats@icov[[g]]
      m_w <- lavsamplestats@cov[[g]]
      omega[[g]] <- (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@nobs[[g]] *
        (w_inv %*% (m_w - sigma_hat[[g]]) %*% w_inv)
      if (meanstructure) {
        diff <- as.matrix(lavsamplestats@mean[[g]] - mu_hat[[g]])
        omega_mu[[g]] <- t(t(diff) %*% w_inv)
      }
    }

    # new in 0.6-18
  if (correlation) {
      diag(omega[[g]]) <- 0
  }
  } # g

  if (meanstructure) attr(omega, "mu") <- omega_mu

  omega
}

lav_model_gradient_dd <- function(lavmodel, g_list = NULL, group = 1L) {
  if (is.null(g_list)) g_list <- lavmodel@GLIST

  #### FIX th + mu!!!!!
  delta_lambda <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "lambda")[[group]]
  delta_tau <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "tau")[[group]]
  delta_nu <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "nu")[[group]]
  delta_theta <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "theta")[[group]]
  delta_beta <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "beta")[[group]]
  delta_psi <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "psi")[[group]]
  delta_alpha <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "alpha")[[group]]
  delta_gamma <-
    lav_model_ddelta_dx(lavmodel, glist = g_list, target = "gamma")[[group]]

  ov_y_dummy_ov_idx <- lavmodel@ov.y.dummy.ov.idx[[group]]
  ov_x_dummy_ov_idx <- lavmodel@ov.x.dummy.ov.idx[[group]]
  ov_y_dummy_lv_idx <- lavmodel@ov.y.dummy.lv.idx[[group]]
  ov_x_dummy_lv_idx <- lavmodel@ov.x.dummy.lv.idx[[group]]
  ov_dummy_idx <- c(ov_y_dummy_ov_idx, ov_x_dummy_ov_idx)
  lv_dummy_idx <- c(ov_y_dummy_lv_idx, ov_x_dummy_lv_idx)
  num_idx <- lavmodel@num.idx[[group]]

  # fix Delta's...
  mm_idx <- lav_model_get_mm_idx(lavmodel)
  mm_in_group <- mm_idx[[group]]
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
