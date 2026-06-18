lav_model_grad_mml <- function(lavmodel = NULL,
                                   mm_theta = NULL,
                                   th = NULL,
                                   glist = NULL,
                                   group = 1L,
                                   lavdata = NULL,
                                   sample_mean = NULL,
                                   sample_mean_x = NULL,
                                   lavcache = NULL) {
  if (lavmodel@link == "logit") {
    lav_msg_stop(gettext("logit link not implemented yet; use probit"))
  }

  # shortcut
  ov_y_dummy_ov_idx <- lavmodel@ov.y.dummy.ov.idx[[group]]
  ov_x_dummy_ov_idx <- lavmodel@ov.x.dummy.ov.idx[[group]]
  ov_y_dummy_lv_idx <- lavmodel@ov.y.dummy.lv.idx[[group]]
  ov_x_dummy_lv_idx <- lavmodel@ov.x.dummy.lv.idx[[group]]
  # ov_dummy_idx <- c(ov_y_dummy_ov_idx, ov_x_dummy_ov_idx)
  lv_dummy_idx <- c(ov_y_dummy_lv_idx, ov_x_dummy_lv_idx)
  th_idx <- lavmodel@th.idx[[group]]
  num_idx <- lavmodel@num.idx[[group]]
  ord_idx <- unique(th_idx[th_idx > 0L])


  # data for this group
  x_1 <- lavdata@X[[group]]
  nobs <- nrow(x_1)
  nvar <- ncol(x_1)
  exo <- lavdata@eXo[[group]]

  # MLIST (for veta and yhat)
  mm_in_group <- 1:lavmodel@nmat[group] + cumsum(c(0, lavmodel@nmat))[group]
  mlist <- glist[mm_in_group]

  # quadrature points
  gh <- lavcache[[group]]$GH
  n_gh <- nrow(gh$x)
  nfac <- ncol(gh$x)

  # compute VETAx (latent lv only)
  # VETAx <- lav_lisrel_vetax(MLIST = MLIST, lv.dummy.idx = lv.dummy.idx)
  vetax <- lav_lisrel_vetax(mlist = mlist)
  # check for negative values?
  if (any(diag(vetax) < 0)) {
    lav_msg_warn(gettext("--- VETAx contains negative values"))
    print(vetax)
    return(0)
  }

  # cholesky?
  # if(is.null(lavmodel@control$cholesky)) {
  cholesky <- TRUE
  # } else {
  #    CHOLESKY <- as.logical(lavmodel@control$cholesky)
  # if(nfac > 1L && !CHOLESKY) {
  #    warning("lavaan WARNING: CHOLESKY is OFF but nfac > 1L")
  # }
  # }

  if (!cholesky) {
    # we should still 'scale' the factors, if std.lv=FALSE
    eta_sd <- sqrt(diag(vetax))
  } else {
    # cholesky takes care of scaling
    eta_sd <- rep(1, nfac)
    tchol_veta <- try(chol(vetax), silent = TRUE)
    if (inherits(tchol_veta, "try-error")) {
      lav_msg_warn(gettext("--- VETAx not positive definite"))
      print(vetax)
      return(0)
    }
    if (!is.null(mlist$alpha) || !is.null(mlist$gamma)) {
      eetax <- lav_lisrel_eetax(
        mlist = mlist, exo = exo, n = nobs,
        sample_mean = sample_mean,
        ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
        ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
        ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
        ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
      )
      # if(length(lv.dummy.idx) > 0L) {
      #    EETAx <- EETAx[,-lv.dummy.idx,drop=FALSE]
      # }
    }
  }

  # prepare common stuff
  # fix Lambda?
  mm_lambda <- lav_lisrel_lambda(
    mlist = mlist,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  # fix ALPHA
  mm_alpha <- mlist$alpha
  if (is.null(mm_alpha)) {
    mm_alpha <- numeric(nfac)
  } else if (length(lv_dummy_idx)) {
    mm_alpha <- mm_alpha[-lv_dummy_idx, , drop = FALSE]
  }

  # Beta?
  mm_beta <- mlist$beta
  if (is.null(mm_beta)) {
    lambda__ib_inv <- mm_lambda
  } else {
    tmp <- -mm_beta
    nr <- nrow(mm_beta)
    i <- seq_len(nr)
    tmp[cbind(i, i)] <- 1
    ib_inv <- solve(tmp)
    lambda__ib_inv <- mlist$lambda %*% ib_inv ## no need to FIX???
    if (length(lv_dummy_idx) > 0L) {
      lambda__ib_inv <- lambda__ib_inv[, -lv_dummy_idx, drop = FALSE]
    }

    # fix BETA
    if (length(lv_dummy_idx)) {
      mm_beta <- mlist$beta[-lv_dummy_idx, -lv_dummy_idx, drop = FALSE]
    }
    tmp <- -mm_beta
    nr <- nrow(mm_beta)
    i <- seq_len(nr)
    tmp[cbind(i, i)] <- 1
    ib_inv <- solve(tmp)
  }

  # fix GAMMA
  mm_gamma <- mlist$gamma
  if (is.null(mm_gamma)) {
    alpha_gamma_e_xo <- matrix(as.numeric(mm_alpha), nobs, nfac, byrow = TRUE)
  } else if (length(lv_dummy_idx)) {
    mm_gamma <- mm_gamma[-lv_dummy_idx, , drop = FALSE]
    alpha_gamma_e_xo <- sweep(exo %*% t(mm_gamma),
      MARGIN = 2, STATS = as.numeric(mm_alpha), FUN = "+"
    )
  }

  # Delta
  ## DD <- lavcache[[group]]$DD
  dd <- lav_model_grad_dd(lavmodel, g_list = glist, group = group)

  ## FIXME!!! do this analytically...
  x <- lav_model_get_parameters(lavmodel = lavmodel, glist = mlist)
  d_vetadx <- function(x, lavmodel = lavmodel, g = 1L) {
    glist <- lav_model_x2glist(lavmodel, x = x, type = "free")
    vetax <- lav_model_vetax(lavmodel, glist = glist)[[g]]
    if (cholesky) {
      s <- chol(vetax) ### FIXME or t(chol())????
    } else {
      s <- diag(sqrt(diag(vetax)))
    }
    s
  }
  delta_s <- lav_func_jacobian_simple(func = d_vetadx,
               x = x, lavmodel = lavmodel, g = group)
  dd$S <- delta_s

  # compute dL/dx for each node
  # dLdx <-  matrix(0, nGH, lavmodel@nx.free)
  d_fyp <- matrix(0, nobs, lavmodel@nx.free)
  sum_log_fy <- matrix(0, nrow = n_gh, ncol = nobs)
  for (q in 1:n_gh) {
    # contribution to dFYp for this q
    d_fyp_q <- matrix(0, nobs, lavmodel@nx.free)

    # current value(s) for ETA
    eta <- ksi <- gh$x[q, , drop = FALSE]

    # rescale/unwhiten
    if (cholesky) {
      eta <- eta %*% tchol_veta
    } else {
      # no unit scale? (un-standardize)
      eta <- sweep(eta, MARGIN = 2, STATS = eta_sd, FUN = "*")
    }

    # eta_i = alpha + BETA eta_i + GAMMA eta_i + error
    #
    # - direct effect of BETA is already in VETAx, and hence tchol.VETA
    # - need to add alpha, and GAMMA eta_i
    if (!is.null(mlist$alpha) || !is.null(mlist$gamma)) {
      eta <- sweep(eetax, MARGIN = 2, STATS = eta, FUN = "+")
    }

    # again, compute yhat for this node (eta)
    if (lavmodel@conditional.x) {
      yhat <- lav_lisrel_eyetax(
        mlist = mlist, exo = exo,
        eta = eta, sample_mean = sample_mean,
        ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
        ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
        ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
        ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
      )
    } else {
      yhat <- lav_lisrel_eyetax3(
        mlist = mlist,
        eta = eta, sample_mean = sample_mean,
        mean_x = sample_mean_x,
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[group]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[group]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[group]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[group]]
      )
    }

    # compute fy.var, for this node (eta): P(Y_i =  y_i | eta_i, x_i)
    log_fy_var <- lav_predict_fy_internal(
      x = x_1, yhat = yhat,
      th = th, mm_theta = mm_theta,
      num_idx = num_idx, th_idx = th_idx,
      link = lavmodel@link, log_1 = TRUE
    )

    # if log, fy is just the sum of log.fy.var
    log_fy <- apply(log_fy_var, 1L, sum)

    # store log likelihoods for this node
    sum_log_fy[q, ] <- log_fy

    # FY
    fy <- exp(log_fy_var) ### FIXME log/exp/log/...
    lik_eta <- apply(fy, 1, prod)
    # fyp <- LIK.eta * GH$w[q]

    ######### dFY_p ###########################################
    # note, dFYp is actually 1/FY[,p] * dFYp

    pre <- matrix(0, nobs, nvar)
    if (length(num_idx) > 0L) {
      tmp <- x_1[, num_idx, drop = FALSE] - yhat[, num_idx, drop = FALSE]
      theta_var <- diag(mm_theta)[num_idx]
      pre[, num_idx] <- sweep(tmp, MARGIN = 2, STATS = 1 / theta_var, FUN = "*")
    }

    if (length(ord_idx) > 0L) {
      for (p in ord_idx) {
        # just in case we need theta[v,v] after all...
        sd_v_inv <- 1 / sqrt(mm_theta[p, p])

        # lav_probit
        y <- x_1[, p]
        th_y <- th[th_idx == p]
        th_y_1 <- c(-Inf, th_y, Inf)
        ncat <- length(th_y) + 1L
        nth <- ncat - 1L
        y1 <- matrix(1:nth, nobs, nth, byrow = TRUE) == y
        y2 <- matrix(1:nth, nobs, nth, byrow = TRUE) == (y - 1L)
        z1 <- pmin(100, th_y_1[y + 1L] - yhat[, p])
        z2 <- pmax(-100, th_y_1[y + 1L - 1L] - yhat[, p])
        p1 <- dnorm(z1)
        p2 <- dnorm(z2)
        # probits = p1 - p2

        pre[, p] <- -1 * (p1 - p2) * sd_v_inv * (1 / fy[, p])

        # [nobx * n.th]
        # dth <- -1 * (Y2*p2 - Y1*p1) * sd.v.inv
        dth <- -1 * (y2 * p2 - y1 * p1) * sd_v_inv * (1 / fy[, p])
        d_fyp_q <- d_fyp_q +
          (dth %*% dd$tau[which(th_idx == p), , drop = FALSE])
      }
    }

    if (length(num_idx) > 0L) {
      # THETA (num only)
      dsigma2 <- sweep(0.5 * pre[, num_idx, drop = FALSE] *
        pre[, num_idx, drop = FALSE],
        MARGIN = 2,
        STATS = 1 / (2 * theta_var), FUN = "-"
      )
      d_fyp_q <- d_fyp_q + (dsigma2 %*% dd$theta)

      # NU (num only)
      dnu <- pre[, num_idx, drop = FALSE]
      d_fyp_q <- d_fyp_q + (dnu %*% dd$nu)
    }

    # LAMBDA
    if (nrow(eta) == 1L) {
      dlambda <- pre %*% eta
      ### FIXME!!!!!
    } else {
      dlambda <- matrix(apply(pre, 2, function(x) x * eta), nobs, )
      # dlambda <- sweep(PRE, MARGIN=1, STATS=eta, FUN="*")
    }
    d_fyp_q <- d_fyp_q + (dlambda %*% dd$lambda)

    # PSI
    # if(nrow(ksi) == 1L) {
    dpsi <- pre %*% kronecker(mm_lambda[, , drop = FALSE], ksi)
    # } else {
    #    dpsi <- PRE * kronecker(LAMBDA[,,drop=FALSE], ksi)
    # }
    d_fyp_q <- d_fyp_q + (dpsi %*% dd$S)

    # KAPPA
    if (length(ov_y_dummy_ov_idx) > 0L) {
      dkappa <- matrix(apply(
        pre[, ov_y_dummy_ov_idx, drop = FALSE], 2,
        function(x) x * exo
      ), nobs, )
      d_fyp_q <- d_fyp_q + (dkappa %*% dd$kappa)
    }

    # GAMMA
    if (!is.null(exo)) {
      dgamma <- matrix(apply(
        pre %*% lambda__ib_inv, 2,
        function(x) x * exo
      ), nobs, )
      d_fyp_q <- d_fyp_q + (dgamma %*% dd$gamma)
    }

    # BETA
    if (!is.null(mm_beta)) {
      # tmp <- kronecker(LAMBDA, ALPHA.GAMMA.eXo) %*%
      #         t( kronecker(t(IB.inv), IB.inv) )
      # dbeta <- apply(matrix(as.numeric(PRE) * tmp, nobs, ), 1, sum)
      dbeta <- matrix(apply(
        pre %*% lambda__ib_inv, 2,
        function(x) x * alpha_gamma_e_xo
      ), nobs, )
      d_fyp_q <- d_fyp_q + (dbeta %*% dd$beta)
    }

    d_fyp <- d_fyp + ((lik_eta * gh$w[q]) * d_fyp_q)
  }

  lik <- as.numeric(t(gh$w) %*% exp(sum_log_fy))
  # avoid underflow
  idx <- which(lik < exp(-600))
  if (length(idx) > 0L) {
    lik[idx] <- exp(-600)
  }

  d_fyp <- 1 / lik * d_fyp

  dx <- apply(d_fyp, 2, sum)

  # integration
  # dx <- apply(as.numeric(GH$w) * dLdx, 2, sum)

  # minimize
  dx <- -1 * dx

  dx
}
