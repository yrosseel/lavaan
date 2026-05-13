

# marginal ML
lav_model_lik_mml <- function(lavmodel = NULL,
                              mm_theta = NULL,
                              th = NULL,
                              glist = NULL,
                              group = 1L,
                              lavdata = NULL,
                              sample_mean = NULL,
                              sample_mean_x = NULL,
                              lavcache = NULL) {
  conditional_x <- lavmodel@conditional.x

  # data for this group
  x <- lavdata@X[[group]]
  nobs <- nrow(x)
  # nvar <- ncol(x)
  exo <- lavdata@eXo[[group]]

  # MLIST (for veta and yhat)
  mm_in_group <- 1:lavmodel@nmat[group] + cumsum(c(0, lavmodel@nmat))[group]
  mlist <- glist[mm_in_group]

  # quadrature points
  gh <- lavcache[[group]]$GH
  n_gh <- nrow(gh$x)
  # nfac <- ncol(gh$x)

  # compute VETAx (latent lv only)
  lv_dummy_idx <- c(
    lavmodel@ov.y.dummy.lv.idx[[group]],
    lavmodel@ov.x.dummy.lv.idx[[group]]
  )
  vetax <- lav_lisrel_vetax(
    mlist = mlist,
    lv_dummy_idx = lv_dummy_idx
  )
  # VETAx <- lav_lisrel_vetax(MLIST = MLIST)
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
    tchol_veta <- try(chol(vetax), silent = TRUE)
    if (inherits(tchol_veta, "try-error")) {
      lav_msg_warn(gettext("--- VETAx not positive definite"))
      print(vetax)
      return(0)
    }
    if (!is.null(mlist$alpha) || !is.null(mlist$gamma)) {
      if (conditional_x) {
        eetax <- lav_lisrel_eetax(
          mlist = mlist, exo = exo, n = nobs,
          sample_mean = sample_mean,
          ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[group]],
          ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[group]],
          ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[group]],
          ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[group]]
        )
      } else {
        eeta <- lav_lisrel_eeta(
          mlist = mlist,
          mean_x = sample_mean_x,
          sample_mean = sample_mean,
          ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[group]],
          ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[group]],
          ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[group]],
          ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[group]]
        )
      }
      # if(length(lv.dummy.idx) > 0L) {
      #    EETAx <- EETAx[,-lv.dummy.idx,drop=FALSE]
      # }
    }
  }

  # compute (log)lik for each node, for each observation
  sum_log_fy <- matrix(0, nrow = n_gh, ncol = nobs)
  for (q in 1:n_gh) {
    # current value(s) for ETA
    # eta <- matrix(0, nrow = 1, ncol = ncol(MLIST$lambda))

    # non-dummy elements -> quadrature points
    # eta[1L, -lv.dummy.idx] <- GH$x[q,,drop=FALSE]
    xq <- gh$x[q, , drop = FALSE]

    # rescale/unwhiten
    if (cholesky) {
      # un-orthogonalize
      xq <- xq %*% tchol_veta
    } else {
      # no unit scale? (un-standardize)
      xq <- sweep(xq, MARGIN = 2, STATS = eta_sd, FUN = "*")
    }

    eta <- matrix(0, nrow = 1, ncol = ncol(mlist$lambda))
    if (length(lv_dummy_idx) > 0L) {
      eta[, -lv_dummy_idx] <- xq
    } else {
      eta <- xq
    }

    # eta_i = alpha + BETA eta_i + GAMMA eta_i + error
    #
    # - direct effect of BETA is already in VETAx, and hence tchol.VETA
    # - need to add alpha, and GAMMA eta_i
    if (!is.null(mlist$alpha) || !is.null(mlist$gamma)) {
      if (conditional_x) {
        eta <- sweep(eetax, MARGIN = 2, STATS = eta, FUN = "+")
      } else {
        eta <- eta + eeta
      }
    }

    # compute yhat for this node (eta)
    if (lavmodel@conditional.x) {
      yhat <- lav_lisrel_eyetax(
        mlist = mlist, exo = exo,
        eta = eta, sample_mean = sample_mean,
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[group]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[group]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[group]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[group]]
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
      X = x, yhat = yhat,
      TH = th, THETA = mm_theta,
      num.idx = lavmodel@num.idx[[group]],
      th.idx = lavmodel@th.idx[[group]],
      link = lavmodel@link, log. = TRUE
    )

    # if log, fy is just the sum of log.fy.var
    log_fy <- apply(log_fy_var, 1L, sum)

    # store log likelihoods for this node
    sum_log_fy[q, ] <- log_fy
  }

  # integration
  lik <- as.numeric(t(gh$w) %*% exp(sum_log_fy))

  # avoid underflow
  idx <- which(lik < exp(-600))
  if (length(idx) > 0L) {
    lik[idx] <- exp(-600)
  }

  lik
}
