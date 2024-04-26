lav_model_gradient_mml <- function(lavmodel = NULL,
                                   THETA = NULL,
                                   TH = NULL,
                                   GLIST = NULL,
                                   group = 1L,
                                   lavdata = NULL,
                                   sample.mean = NULL,
                                   sample.mean.x = NULL,
                                   lavcache = NULL) {
  if (lavmodel@link == "logit") {
    lav_msg_stop(gettext("logit link not implemented yet; use probit"))
  }

  # shortcut
  ov.y.dummy.ov.idx <- lavmodel@ov.y.dummy.ov.idx[[group]]
  ov.x.dummy.ov.idx <- lavmodel@ov.x.dummy.ov.idx[[group]]
  ov.y.dummy.lv.idx <- lavmodel@ov.y.dummy.lv.idx[[group]]
  ov.x.dummy.lv.idx <- lavmodel@ov.x.dummy.lv.idx[[group]]
  ov.dummy.idx <- c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)
  th.idx <- lavmodel@th.idx[[group]]
  num.idx <- lavmodel@num.idx[[group]]
  ord.idx <- unique(th.idx[th.idx > 0L])


  # data for this group
  X <- lavdata@X[[group]]
  nobs <- nrow(X)
  nvar <- ncol(X)
  eXo <- lavdata@eXo[[group]]

  # MLIST (for veta and yhat)
  mm.in.group <- 1:lavmodel@nmat[group] + cumsum(c(0, lavmodel@nmat))[group]
  MLIST <- GLIST[mm.in.group]

  # quadrature points
  GH <- lavcache[[group]]$GH
  nGH <- nrow(GH$x)
  nfac <- ncol(GH$x)

  # compute VETAx (latent lv only)
  # VETAx <- computeVETAx.LISREL(MLIST = MLIST, lv.dummy.idx = lv.dummy.idx)
  VETAx <- computeVETAx.LISREL(MLIST = MLIST)
  # check for negative values?
  if (any(diag(VETAx) < 0)) {
    lav_msg_warn(gettext("--- VETAx contains negative values"))
    print(VETAx)
    return(0)
  }

  # cholesky?
  # if(is.null(lavmodel@control$cholesky)) {
  CHOLESKY <- TRUE
  # } else {
  #    CHOLESKY <- as.logical(lavmodel@control$cholesky)
  # if(nfac > 1L && !CHOLESKY) {
  #    warning("lavaan WARNING: CHOLESKY is OFF but nfac > 1L")
  # }
  # }

  if (!CHOLESKY) {
    # we should still 'scale' the factors, if std.lv=FALSE
    ETA.sd <- sqrt(diag(VETAx))
  } else {
    # cholesky takes care of scaling
    ETA.sd <- rep(1, nfac)
    tchol.VETA <- try(chol(VETAx), silent = TRUE)
    if (inherits(tchol.VETA, "try-error")) {
      lav_msg_warn(gettext("--- VETAx not positive definite"))
      print(VETAx)
      return(0)
    }
    if (!is.null(MLIST$alpha) || !is.null(MLIST$gamma)) {
      EETAx <- computeEETAx.LISREL(
        MLIST = MLIST, eXo = eXo, N = nobs,
        sample.mean = sample.mean,
        ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
        ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
        ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
        ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
      )
      # if(length(lv.dummy.idx) > 0L) {
      #    EETAx <- EETAx[,-lv.dummy.idx,drop=FALSE]
      # }
    }
  }

  # prepare common stuff
  # fix Lambda?
  LAMBDA <- computeLAMBDA.LISREL(
    MLIST = MLIST,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # fix ALPHA
  ALPHA <- MLIST$alpha
  if (is.null(ALPHA)) {
    ALPHA <- numeric(nfac)
  } else if (length(lv.dummy.idx)) {
    ALPHA <- ALPHA[-lv.dummy.idx, , drop = FALSE]
  }

  # Beta?
  BETA <- MLIST$beta
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    tmp <- -BETA
    nr <- nrow(BETA)
    i <- seq_len(nr)
    tmp[cbind(i, i)] <- 1
    IB.inv <- solve(tmp)
    LAMBDA..IB.inv <- MLIST$lambda %*% IB.inv ## no need to FIX???
    if (length(lv.dummy.idx) > 0L) {
      LAMBDA..IB.inv <- LAMBDA..IB.inv[, -lv.dummy.idx, drop = FALSE]
    }

    # fix BETA
    if (length(lv.dummy.idx)) {
      BETA <- MLIST$beta[-lv.dummy.idx, -lv.dummy.idx, drop = FALSE]
    }
    tmp <- -BETA
    nr <- nrow(BETA)
    i <- seq_len(nr)
    tmp[cbind(i, i)] <- 1
    IB.inv <- solve(tmp)
  }

  # fix GAMMA
  GAMMA <- MLIST$gamma
  if (is.null(GAMMA)) {
    ALPHA.GAMMA.eXo <- matrix(as.numeric(ALPHA), nobs, nfac, byrow = TRUE)
  } else if (length(lv.dummy.idx)) {
    GAMMA <- GAMMA[-lv.dummy.idx, , drop = FALSE]
    ALPHA.GAMMA.eXo <- sweep(eXo %*% t(GAMMA),
      MARGIN = 2, STATS = as.numeric(ALPHA), FUN = "+"
    )
  }

  # Delta
  ## DD <- lavcache[[group]]$DD
  DD <- lav_model_gradient_DD(lavmodel, GLIST = GLIST, group = group)

  ## FIXME!!! do this analytically...
  x <- lav_model_get_parameters(lavmodel = lavmodel, GLIST = MLIST)
  dVetadx <- function(x, lavmodel = lavmodel, g = 1L) {
    GLIST <- lav_model_x2GLIST(lavmodel, x = x, type = "free")
    VETAx <- computeVETAx(lavmodel, GLIST = GLIST)[[g]]
    if (CHOLESKY) {
      S <- chol(VETAx) ### FIXME or t(chol())????
    } else {
      S <- diag(sqrt(diag(VETAx)))
    }
    S
  }
  Delta.S <- lav_func_jacobian_simple(func = dVetadx, x = x, lavmodel = lavmodel, g = group)
  DD$S <- Delta.S

  # compute dL/dx for each node
  # dLdx <-  matrix(0, nGH, lavmodel@nx.free)
  dFYp <- matrix(0, nobs, lavmodel@nx.free)
  SUM.LOG.FY <- matrix(0, nrow = nGH, ncol = nobs)
  for (q in 1:nGH) {
    # contribution to dFYp for this q
    dFYp.q <- matrix(0, nobs, lavmodel@nx.free)

    # current value(s) for ETA
    eta <- ksi <- GH$x[q, , drop = FALSE]

    # rescale/unwhiten
    if (CHOLESKY) {
      eta <- eta %*% tchol.VETA
    } else {
      # no unit scale? (un-standardize)
      eta <- sweep(eta, MARGIN = 2, STATS = ETA.sd, FUN = "*")
    }

    # eta_i = alpha + BETA eta_i + GAMMA eta_i + error
    #
    # - direct effect of BETA is already in VETAx, and hence tchol.VETA
    # - need to add alpha, and GAMMA eta_i
    if (!is.null(MLIST$alpha) || !is.null(MLIST$gamma)) {
      eta <- sweep(EETAx, MARGIN = 2, STATS = eta, FUN = "+")
    }

    # again, compute yhat for this node (eta)
    if (lavmodel@conditional.x) {
      yhat <- computeEYetax.LISREL(
        MLIST = MLIST, eXo = eXo,
        ETA = eta, sample.mean = sample.mean,
        ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
        ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
        ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
        ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
      )
    } else {
      yhat <- computeEYetax3.LISREL(
        MLIST = MLIST,
        ETA = eta, sample.mean = sample.mean,
        mean.x = sample.mean.x,
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]]
      )
    }

    # compute fy.var, for this node (eta): P(Y_i =  y_i | eta_i, x_i)
    log.fy.var <- lav_predict_fy_internal(
      X = X, yhat = yhat,
      TH = TH, THETA = THETA,
      num.idx = num.idx, th.idx = th.idx,
      link = lavmodel@link, log. = TRUE
    )

    # if log, fy is just the sum of log.fy.var
    log.fy <- apply(log.fy.var, 1L, sum)

    # store log likelihoods for this node
    SUM.LOG.FY[q, ] <- log.fy

    # FY
    FY <- exp(log.fy.var) ### FIXME log/exp/log/...
    LIK.eta <- apply(FY, 1, prod)
    # fyp <- LIK.eta * GH$w[q]

    ######### dFY_p ###########################################
    # note, dFYp is actually 1/FY[,p] * dFYp

    PRE <- matrix(0, nobs, nvar)
    if (length(num.idx) > 0L) {
      tmp <- X[, num.idx, drop = FALSE] - yhat[, num.idx, drop = FALSE]
      theta.var <- diag(THETA)[num.idx]
      PRE[, num.idx] <- sweep(tmp, MARGIN = 2, STATS = 1 / theta.var, FUN = "*")
    }

    if (length(ord.idx) > 0L) {
      for (p in ord.idx) {
        # just in case we need theta[v,v] after all...
        sd.v.inv <- 1 / sqrt(THETA[p, p])

        # lav_probit
        y <- X[, p]
        th.y <- TH[th.idx == p]
        TH.Y <- c(-Inf, th.y, Inf)
        ncat <- length(th.y) + 1L
        nth <- ncat - 1L
        Y1 <- matrix(1:nth, nobs, nth, byrow = TRUE) == y
        Y2 <- matrix(1:nth, nobs, nth, byrow = TRUE) == (y - 1L)
        z1 <- pmin(100, TH.Y[y + 1L] - yhat[, p])
        z2 <- pmax(-100, TH.Y[y + 1L - 1L] - yhat[, p])
        p1 <- dnorm(z1)
        p2 <- dnorm(z2)
        # probits = p1 - p2

        PRE[, p] <- -1 * (p1 - p2) * sd.v.inv * (1 / FY[, p])

        # [nobx * n.th]
        # dth <- -1 * (Y2*p2 - Y1*p1) * sd.v.inv
        dth <- -1 * (Y2 * p2 - Y1 * p1) * sd.v.inv * (1 / FY[, p])
        dFYp.q <- dFYp.q +
          (dth %*% DD$tau[which(th.idx == p), , drop = FALSE])
      }
    }

    if (length(num.idx) > 0L) {
      # THETA (num only)
      dsigma2 <- sweep(0.5 * PRE[, num.idx] * PRE[, num.idx],
        MARGIN = 2,
        STATS = 1 / (2 * theta.var), FUN = "-"
      )
      dFYp.q <- dFYp.q + (dsigma2 %*% DD$theta)

      # NU (num only)
      dnu <- PRE[, num.idx]
      dFYp.q <- dFYp.q + (dnu %*% DD$nu)
    }

    # LAMBDA
    if (nrow(eta) == 1L) {
      dlambda <- PRE %*% eta
      ### FIXME!!!!!
    } else {
      dlambda <- matrix(apply(PRE, 2, function(x) x * eta), nobs, )
      # dlambda <- sweep(PRE, MARGIN=1, STATS=eta, FUN="*")
    }
    dFYp.q <- dFYp.q + (dlambda %*% DD$lambda)

    # PSI
    # if(nrow(ksi) == 1L) {
    dpsi <- PRE %*% kronecker(LAMBDA[, , drop = FALSE], ksi)
    # } else {
    #    dpsi <- PRE * kronecker(LAMBDA[,,drop=FALSE], ksi)
    # }
    dFYp.q <- dFYp.q + (dpsi %*% DD$S)

    # KAPPA
    if (length(ov.y.dummy.ov.idx) > 0L) {
      dkappa <- matrix(apply(
        PRE[, ov.y.dummy.ov.idx, drop = FALSE], 2,
        function(x) x * eXo
      ), nobs, )
      dFYp.q <- dFYp.q + (dkappa %*% DD$kappa)
    }

    # GAMMA
    if (!is.null(eXo)) {
      dgamma <- matrix(apply(
        PRE %*% LAMBDA..IB.inv, 2,
        function(x) x * eXo
      ), nobs, )
      dFYp.q <- dFYp.q + (dgamma %*% DD$gamma)
    }

    # BETA
    if (!is.null(BETA)) {
      # tmp <- kronecker(LAMBDA, ALPHA.GAMMA.eXo) %*%
      #         t( kronecker(t(IB.inv), IB.inv) )
      # dbeta <- apply(matrix(as.numeric(PRE) * tmp, nobs, ), 1, sum)
      dbeta <- matrix(apply(
        PRE %*% LAMBDA..IB.inv, 2,
        function(x) x * ALPHA.GAMMA.eXo
      ), nobs, )
      dFYp.q <- dFYp.q + (dbeta %*% DD$beta)
    }

    dFYp <- dFYp + ((LIK.eta * GH$w[q]) * dFYp.q)
  }

  lik <- as.numeric(t(GH$w) %*% exp(SUM.LOG.FY))
  # avoid underflow
  idx <- which(lik < exp(-600))
  if (length(idx) > 0L) {
    lik[idx] <- exp(-600)
  }

  dFYp <- 1 / lik * dFYp

  dx <- apply(dFYp, 2, sum)

  # integration
  # dx <- apply(as.numeric(GH$w) * dLdx, 2, sum)

  # minimize
  dx <- -1 * dx

  dx
}
