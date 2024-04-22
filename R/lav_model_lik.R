# casewise likelihoods

# closed-form marginal likelihood
# - classic SEM models, continous observed variables only
lav_model_lik_ml <- function(lavmodel = NULL,
                             GLIST = NULL,
                             lavdata = NULL,
                             lavsamplestats = NULL) {


}


# marginal ML
lav_model_lik_mml <- function(lavmodel = NULL,
                              THETA = NULL,
                              TH = NULL,
                              GLIST = NULL,
                              group = 1L,
                              lavdata = NULL,
                              sample.mean = NULL,
                              sample.mean.x = NULL,
                              lavcache = NULL) {
  conditional.x <- lavmodel@conditional.x

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
  lv.dummy.idx <- c(
    lavmodel@ov.y.dummy.lv.idx[[group]],
    lavmodel@ov.x.dummy.lv.idx[[group]]
  )
  VETAx <- computeVETAx.LISREL(
    MLIST = MLIST,
    lv.dummy.idx = lv.dummy.idx
  )
  # VETAx <- computeVETAx.LISREL(MLIST = MLIST)
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
    tchol.VETA <- try(chol(VETAx), silent = TRUE)
    if (inherits(tchol.VETA, "try-error")) {
      lav_msg_warn(gettext("--- VETAx not positive definite"))
      print(VETAx)
      return(0)
    }
    if (!is.null(MLIST$alpha) || !is.null(MLIST$gamma)) {
      if (conditional.x) {
        EETAx <- computeEETAx.LISREL(
          MLIST = MLIST, eXo = eXo, N = nobs,
          sample.mean = sample.mean,
          ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
          ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
          ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
          ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]]
        )
      } else {
        EETA <- computeEETA.LISREL(
          MLIST = MLIST,
          mean.x = sample.mean.x,
          sample.mean = sample.mean,
          ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
          ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
          ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
          ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]]
        )
      }
      # if(length(lv.dummy.idx) > 0L) {
      #    EETAx <- EETAx[,-lv.dummy.idx,drop=FALSE]
      # }
    }
  }

  # compute (log)lik for each node, for each observation
  SUM.LOG.FY <- matrix(0, nrow = nGH, ncol = nobs)
  for (q in 1:nGH) {
    # current value(s) for ETA
    # eta <- matrix(0, nrow = 1, ncol = ncol(MLIST$lambda))

    # non-dummy elements -> quadrature points
    # eta[1L, -lv.dummy.idx] <- GH$x[q,,drop=FALSE]
    XQ <- GH$x[q, , drop = FALSE]

    # rescale/unwhiten
    if (CHOLESKY) {
      # un-orthogonalize
      XQ <- XQ %*% tchol.VETA
    } else {
      # no unit scale? (un-standardize)
      XQ <- sweep(XQ, MARGIN = 2, STATS = ETA.sd, FUN = "*")
    }

    eta <- matrix(0, nrow = 1, ncol = ncol(MLIST$lambda))
    if (length(lv.dummy.idx) > 0L) {
      eta[, -lv.dummy.idx] <- XQ
    } else {
      eta <- XQ
    }

    # eta_i = alpha + BETA eta_i + GAMMA eta_i + error
    #
    # - direct effect of BETA is already in VETAx, and hence tchol.VETA
    # - need to add alpha, and GAMMA eta_i
    if (!is.null(MLIST$alpha) || !is.null(MLIST$gamma)) {
      if (conditional.x) {
        eta <- sweep(EETAx, MARGIN = 2, STATS = eta, FUN = "+")
      } else {
        eta <- eta + EETA
      }
    }

    # compute yhat for this node (eta)
    if (lavmodel@conditional.x) {
      yhat <- computeEYetax.LISREL(
        MLIST = MLIST, eXo = eXo,
        ETA = eta, sample.mean = sample.mean,
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[group]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[group]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[group]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[group]]
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
      num.idx = lavmodel@num.idx[[group]],
      th.idx = lavmodel@th.idx[[group]],
      link = lavmodel@link, log. = TRUE
    )

    # if log, fy is just the sum of log.fy.var
    log.fy <- apply(log.fy.var, 1L, sum)

    # store log likelihoods for this node
    SUM.LOG.FY[q, ] <- log.fy
  }

  # integration
  lik <- as.numeric(t(GH$w) %*% exp(SUM.LOG.FY))

  # avoid underflow
  idx <- which(lik < exp(-600))
  if (length(idx) > 0L) {
    lik[idx] <- exp(-600)
  }

  lik
}
