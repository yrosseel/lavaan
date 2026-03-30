# YR 21 March 2015
# new approach to compute 'Gamma': the asymptotic variance matrix of
#                                  sqrt{N} times the
#                                  observed sample statistics (means + varcov)
#
# Gamma = N x ACOV[ ybar, vech(S) ]
#       = NACOV[ ybar, vech(S) ]
#
# - one single function for mean + cov
# - handle 'fixed.x' exogenous covariates


# - YR  3 Dec 2015: allow for conditional.x = TRUE
# - YR 22 Jan 2023: add model.based= argument (if object is lavaan object)
# - YR 30 May 2024: add lav_samplestats_cor_gamma(_NT)

# generic public function (not exported yet)
# input for lavGamma can be lavobject, lavdata, data.frame, or matrix
lavGamma <- function(object, group = NULL, missing = "listwise",
                     ov.names.x = NULL, fixed.x = FALSE, conditional.x = FALSE,
                     meanstructure = FALSE, slopestructure = FALSE,
                     gamma.n.minus.one = FALSE, gamma.unbiased = FALSE,
                     ADF = TRUE, model.based = FALSE, NT.rescale = FALSE,
                     Mplus.WLS = FALSE, add.labels = FALSE) {
  # check object

  # 1. object is lavaan object
  if (inherits(object, "lavaan")) {
    object <- lav_object_check_version(object)
    lav_object_gamma(
      lavobject = object, ADF = ADF,
      model.based = model.based, Mplus.WLS = Mplus.WLS
    )
  } else if (inherits(object, "lavData")) {
    lavdata <- object
    model.based <- FALSE
  } else if (inherits(object, "data.frame") ||
    inherits(object, "matrix")) {
    model.based <- FALSE
    NAMES <- names(object)
    if (!is.null(NAMES) && !is.null(group)) {
      NAMES <- NAMES[-match(group, NAMES)]
    }
    lavdata <- lav_lavdata(
      data = object, group = group,
      ov.names = NAMES, ordered = NULL,
      ov.names.x = ov.names.x,
      lavoptions = list(
        warn = FALSE,
        missing = missing
      )
    )
  } else {
    lav_msg_stop(
      gettextf("lavGamma can not handle objects of class %s",
                lav_msg_view(class(object)))
    )
  }

  # extract data
  Y <- lavdata@X
  if (conditional.x) {
    eXo <- lavdata@eXo
    for (g in seq_len(lavdata@ngroups)) {
      Y[[g]] <- cbind(Y[[g]], eXo[[g]])
    }
  }

  # x.idx
  x.idx <- lapply(
    seq_len(lavdata@ngroups),
    function(g) {
      match(
        lavdata@ov.names.x[[g]],
        lavdata@ov.names[[g]]
      )
    }
  )

  OUT <- lapply(seq_len(lavdata@ngroups), function(g) {
    if (length(lavdata@cluster) > 0L) {
      cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
    } else {
      cluster.idx <- NULL
    }
    if (ADF) {
      out <- lav_samplestats_gamma(
        m_y = Y[[g]],
        m_mu = NULL,
        m_sigma = NULL,
        x_idx = x.idx[[g]],
        cluster_idx = cluster.idx,
        fixed_x = fixed.x,
        conditional_x = conditional.x,
        meanstructure = meanstructure,
        slopestructure = conditional.x,
        gamma_n_minus_one = gamma.n.minus.one,
        unbiased = gamma.unbiased,
        mplus_wls = Mplus.WLS
      )
    } else {
      out <- lav_samplestats_gamma_nt(
        m_y = Y[[g]],
        wt = NULL, # for now
        rescale = NT.rescale,
        x_idx = x.idx[[g]],
        fixed_x = fixed.x,
        conditional_x = conditional.x,
        meanstructure = meanstructure,
        slopestructure = conditional.x
      )
    }
    out
  })

  # todo: labels

  OUT
}



# for internal use -- lavobject or internal slots
lav_object_gamma <- function(lavobject = NULL,
                             # or individual slots
                             lavdata = NULL,
                             lavoptions = NULL,
                             lavsamplestats = NULL,
                             lavh1 = NULL,
                             lavimplied = NULL,
                             # other options
                             ADF = TRUE, model.based = FALSE,
                             Mplus.WLS = FALSE) {
  # extract slots
  if (!is.null(lavobject)) {
    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
    lavsamplestats <- lavobject@SampleStats
    lavh1 <- lavobject@h1
    lavimplied <- lavobject@implied
  }

  missing <- lavoptions$missing
  if (!missing %in% c("listwise", "pairwise")) {
    model.based <- TRUE
  }
  fixed.x <- lavoptions$fixed.x
  conditional.x <- lavoptions$conditional.x
  meanstructure <- lavoptions$meanstructure
  gamma.n.minus.one <- lavoptions$gamma.n.minus.one
  gamma.unbiased <- lavoptions$gamma.unbiased

  if (ADF && model.based && conditional.x) {
    lav_msg_stop(gettext(
      "ADF + model.based + conditional.x is not supported yet."))
  }

  # output container
  OUT <- vector("list", length = lavdata@ngroups)

  # compute Gamma matrix for each group
  for (g in seq_len(lavdata@ngroups)) {
    x.idx <- lavsamplestats@x.idx[[g]]
    COV <- MEAN <- NULL
    if (!ADF || model.based) {
      implied <- lavh1$implied # saturated/unstructured
      if (model.based) {
        implied <- lavimplied # model-based/structured
      }
      if (conditional.x) {
        # convert to joint COV/MEAN
        res.S <- implied$res.cov[[g]]
        res.slopes <- implied$res.slopes[[g]]
        res.int <- implied$res.int[[g]]
        S.xx <- implied$cov.x[[g]]
        M.x <- implied$mean.x[[g]]

        S.yy <- res.S + res.slopes %*% S.xx %*% t(res.slopes)
        S.yx <- res.slopes %*% S.xx
        S.xy <- S.xx %*% t(res.slopes)
        M.y <- res.int + res.slopes %*% M.x

        COV <- rbind(cbind(S.yy, S.yx), cbind(S.xy, S.xx))
        MEAN <- c(M.y, M.x)
      } else {
        # not conditional.x
        COV <- implied$cov[[g]]
        MEAN <- implied$mean[[g]]
      }
    } # COV/MEAN

    if (ADF) {
      if (conditional.x) {
        Y <- cbind(lavdata@X[[g]], lavdata@eXo[[g]])
      } else {
        Y <- lavdata@X[[g]]
      }
      if (length(lavdata@cluster) > 0L) {
        cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster.idx <- NULL
      }
      OUT[[g]] <- lav_samplestats_gamma(
        m_y = Y,
        m_mu = MEAN,
        m_sigma = COV,
        x_idx = x.idx,
        cluster_idx = cluster.idx,
        fixed_x = fixed.x,
        conditional_x = conditional.x,
        meanstructure = meanstructure,
        slopestructure = conditional.x,
        gamma_n_minus_one = gamma.n.minus.one,
        unbiased = gamma.unbiased,
        mplus_wls = Mplus.WLS
      )
    } else {
      OUT[[g]] <- lav_samplestats_gamma_nt(
        m_cov = COV, # joint!
        m_mean = MEAN, # joint!
        x_idx = x.idx,
        fixed_x = fixed.x,
        conditional_x = conditional.x,
        meanstructure = meanstructure,
        slopestructure = conditional.x
      )
    }

    # group.w.free
    if (lavoptions$group.w.free) {
      # checkme!
      OUT[[g]] <- lav_matrix_bdiag(matrix(1, 1, 1), OUT[[g]])
    }

  } # g

  OUT
}





# NOTE:
#  - three types:
#       1) plain         (conditional.x = FALSE, fixed.x = FALSE)
#       2) fixed.x       (conditional.x = FALSE, fixed.x = TRUE)
#       3) conditional.x (conditional.x = TRUE)
#  - if conditional.x = TRUE, we ignore fixed.x (can be TRUE or FALSE)

# NORMAL-THEORY
lav_samplestats_gamma_nt <- function(m_y = NULL, # should include
                                     # eXo if
                                     # conditional.x=TRUE
                                     wt = NULL,
                                     m_cov = NULL, # joint!
                                     m_mean = NULL, # joint!
                                     rescale = FALSE,
                                     x_idx = integer(0L),
                                     fixed_x = FALSE,
                                     conditional_x = FALSE,
                                     meanstructure = FALSE,
                                     slopestructure = FALSE) {
  # check arguments
  if (length(x_idx) == 0L) {
    conditional_x <- FALSE
    fixed_x <- FALSE
  }

  # compute m_cov from m_y
  if (is.null(m_cov)) {
    stopifnot(!is.null(m_y))

    # coerce to matrix
    m_y <- unname(as.matrix(m_y))
    n <- nrow(m_y)
    if (is.null(wt)) {
      m_cov <- cov(m_y)
      if (rescale) {
        m_cov <- m_cov * (n - 1) / n # (normal) ML version
      }
    } else {
      out <- stats::cov.wt(m_y, wt = wt, method = "ML")
      m_cov <- out$cov
    }
  } else {
    if (!missing(rescale)) {
      lav_msg_warn(gettext("rescale= argument has no effect if m_cov is given"))
    }
    if (!missing(wt)) {
      lav_msg_warn(gettext("wt= argument has no effect if m_cov is given"))
    }
  }

  # if needed, compute m_mean from m_y
  if (conditional_x && length(x_idx) > 0L && is.null(m_mean) &&
    (meanstructure || slopestructure)) {
    stopifnot(!is.null(m_y))
    if (is.null(wt)) {
      m_mean <- colMeans(m_y, na.rm = TRUE)
    } else {
      m_mean <- out$center
    }
  }

  # rename
  m_s <- m_cov
  m_m <- m_mean

  # unconditional
  if (!conditional_x) {
    # unconditional - stochastic x
    if (!fixed_x) {
      # if (lav_use_lavaanC()) {
      #   m_gamma <- lavaanC::m_kronecker_dup_ginv_pre_post(m_s,
      #                                              multiplicator = 2.0)
      # } else {
        m_gamma <- 2 * lav_matrix_duplication_ginv_pre_post(m_s %x% m_s)
      # }
      if (meanstructure) {
        m_gamma <- lav_matrix_bdiag(m_s, m_gamma)
      }

      # unconditional - fixed x
    } else {
      # handle fixed_x = TRUE

      # cov(m_y|X) = m_a - m_b m_c^{-1} m_b'
      # where m_a = cov(m_y), m_b = cov(m_y,X), m_c = cov(X)
      m_a <- m_s[-x_idx, -x_idx, drop = FALSE]
      m_b <- m_s[-x_idx, x_idx, drop = FALSE]
      m_c <- m_s[x_idx, x_idx, drop = FALSE]
      ybarx <- m_a - m_b %*% solve(m_c, t(m_b))

      # reinsert ybarx in m_y+X (residual) covariance matrix
      ybarx_aug <- matrix(0, nrow = NROW(m_s), ncol = NCOL(m_s))
      ybarx_aug[-x_idx, -x_idx] <- ybarx

      # take difference
      m_r <- m_s - ybarx_aug

      # if (lav_use_lavaanC()) {
      #   gamma_s <- lavaanC::m_kronecker_dup_ginv_pre_post(m_s,
      #                                               multiplicator = 2.0)
      #   gamma_r <- lavaanC::m_kronecker_dup_ginv_pre_post(m_r,
      #                                               multiplicator = 2.0)
      # } else {
        gamma_s <- 2 * lav_matrix_duplication_ginv_pre_post(m_s %x% m_s)
        gamma_r <- 2 * lav_matrix_duplication_ginv_pre_post(m_r %x% m_r)
      # }
      m_gamma <- gamma_s - gamma_r

      if (meanstructure) {
        m_gamma <- lav_matrix_bdiag(ybarx_aug, m_gamma)
      }
    }
  } else {
    # conditional_x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    # regress m_y on X, and compute covariance of residuals 'm_r'
    m_a <- m_s[-x_idx, -x_idx, drop = FALSE]
    m_b <- m_s[-x_idx, x_idx, drop = FALSE]
    m_c <- m_s[x_idx, x_idx, drop = FALSE]
    cov_ybarx <- m_a - m_b %*% solve(m_c) %*% t(m_b)
    # if (lav_use_lavaanC()) {
    #   m_gamma <- lavaanC::m_kronecker_dup_ginv_pre_post(cov_ybarx,
    #                                                  multiplicator = 2.0)
    # } else {
      m_gamma <- 2 *
               lav_matrix_duplication_ginv_pre_post(cov_ybarx %x% cov_ybarx)
    # }

    if (meanstructure || slopestructure) {
      m_mx <- m_m[x_idx]
      c3 <- rbind(
        c(1, m_mx),
        cbind(m_mx, m_c + tcrossprod(m_mx))
      )
      # B3 <- cbind(m_my, m_b + tcrossprod(m_my,m_mx))
    }

    if (meanstructure) {
      if (slopestructure) {
        # a11 <- solve(c3) %x% cov_ybarx
        a11 <- cov_ybarx %x% solve(c3)
      } else {
        # a11 <- solve(c3)[1, 1, drop=FALSE] %x% cov_ybarx
        a11 <- cov_ybarx %x% solve(c3)[1, 1, drop = FALSE]
      }
    } else {
      if (slopestructure) {
        # a11 <- solve(c3)[-1, -1, drop=FALSE] %x% cov_ybarx
        a11 <- cov_ybarx %x% solve(c3)[-1, -1, drop = FALSE]
      } else {
        a11 <- matrix(0, 0, 0)
      }
    }

    m_gamma <- lav_matrix_bdiag(a11, m_gamma)
  }

  m_gamma
}

# NOTE:
#  - three types:
#       1) plain         (conditional_x = FALSE, fixed_x = FALSE)
#       2) fixed_x       (conditional_x = FALSE, fixed_x = TRUE)
#       3) conditional_x (conditional_x = TRUE)
#  - if conditional_x = TRUE, we ignore fixed_x (can be TRUE or FALSE)
#
#  - new in 0.6-1: if Mu/Sigma is provided, compute 'model-based' Gamma
#                  (only if conditional_x = FALSE, for now)
#  - new in 0.6-2: if cluster.idx is not NULL, correct for clustering
#  - new in 0.6-13: add unbiased = TRUE (for the 'plain' setting only)

# ADF THEORY
lav_samplestats_gamma <- function(m_y, # Y+X if cond!
                                  m_mu = NULL,
                                  m_sigma = NULL,
                                  x_idx = integer(0L),
                                  cluster_idx = NULL,
                                  fixed_x = FALSE,
                                  conditional_x = FALSE,
                                  meanstructure = FALSE,
                                  slopestructure = FALSE,
                                  gamma_n_minus_one = FALSE,
                                  unbiased = FALSE,
                                  mplus_wls = FALSE) {
  # coerce to matrix
  m_y <- unname(as.matrix(m_y))
  n <- nrow(m_y)
  p <- ncol(m_y)

  # unbiased?
  if (unbiased) {
    if (conditional_x || fixed_x || !is.null(m_sigma) ||
      !is.null(cluster_idx)) {
      lav_msg_stop(
        gettext("unbiased Gamma only available for the simple
        (not conditional_x or fixed_x or model-based or clustered) setting."))
    } else {
    # data really should be complete
      m_cov <- stats::cov(m_y, use = "pairwise.complete.obs")
      m_cov <- m_cov * (n - 1) / n
      cov_vech <- lav_matrix_vech(m_cov)
    }
  }

  # model-based?
  if (!is.null(m_sigma)) {
    stopifnot(!conditional_x)
    model_based <- TRUE
    if (meanstructure) {
      stopifnot(!is.null(m_mu))
      sigma <- c(as.numeric(m_mu), lav_matrix_vech(m_sigma))
    } else {
      m_mu <- colMeans(m_y, na.rm = TRUE) # for centering!
      sigma <- lav_matrix_vech(m_sigma)
    }
  } else {
    model_based <- FALSE
  }

  # denominator
  if (gamma_n_minus_one) {
    n <- n - 1
  }

  # check arguments
  if (length(x_idx) == 0L) {
    conditional_x <- FALSE
    fixed_x <- FALSE
  }
  if (mplus_wls) {
    stopifnot(!conditional_x, !fixed_x)
  }

  if (!conditional_x && !fixed_x) {
    # center, so we can use crossprod instead of cov
    if (model_based) {
      m_yc <- t(t(m_y) - as.numeric(m_mu))
    } else {
      m_yc <- t(t(m_y) - colMeans(m_y, na.rm = TRUE))
    }

    # create z where the rows_i contain the following elements:
    #  - Y_i (if meanstructure is TRUE)
    #  - vech(Yc_i' %*% Yc_i) where Yc_i are the residuals
    idx1 <- lav_matrix_vech_col_idx(p)
    idx2 <- lav_matrix_vech_row_idx(p)
    if (meanstructure) {
      z <- cbind(m_y, m_yc[, idx1, drop = FALSE] *
        m_yc[, idx2, drop = FALSE])
    } else {
      z <- (m_yc[, idx1, drop = FALSE] *
        m_yc[, idx2, drop = FALSE])
    }

    if (model_based) {
      if (meanstructure) {
        stopifnot(!is.null(m_mu))
        sigma <- c(as.numeric(m_mu), lav_matrix_vech(m_sigma))
      } else {
        sigma <- lav_matrix_vech(m_sigma)
      }
      zc <- t(t(z) - sigma)
    } else {
      zc <- t(t(z) - colMeans(z, na.rm = TRUE))
    }

    # clustered?
    if (length(cluster_idx) > 0L) {
      zc <- rowsum(zc, cluster_idx)
    }

    if (anyNA(zc)) {
      m_gamma <- lav_matrix_crossprod(zc) / n
    } else {
      m_gamma <- base::crossprod(zc) / n
    }
  } else if (!conditional_x && fixed_x) {
    if (model_based) {
      m_yc <- t(t(m_y) - as.numeric(m_mu))
      y_bar <- colMeans(m_y, na.rm = TRUE)
      res_cov <- (m_sigma[-x_idx, -x_idx, drop = FALSE] -
        m_sigma[-x_idx, x_idx, drop = FALSE] %*%
        solve(m_sigma[x_idx, x_idx, drop = FALSE]) %*%
        m_sigma[x_idx, -x_idx, drop = FALSE])
      res_slopes <- (solve(m_sigma[x_idx, x_idx, drop = FALSE]) %*%
        m_sigma[x_idx, -x_idx, drop = FALSE])
      res_int <- (y_bar[-x_idx] -
        as.numeric(colMeans(m_y[, x_idx, drop = FALSE],
          na.rm = TRUE
        ) %*% res_slopes))
      x_bar <- y_bar[x_idx]
      yhat__bar <- as.numeric(res_int + as.numeric(x_bar) %*% res_slopes)
      yhat_bar <- numeric(p)
      yhat_bar[-x_idx] <- yhat__bar
      yhat_bar[x_idx] <- x_bar
      yhat_cov <- m_sigma
      yhat_cov[-x_idx, -x_idx] <- m_sigma[-x_idx, -x_idx] - res_cov


      yhat <- cbind(1, m_y[, x_idx]) %*% rbind(res_int, res_slopes)
      y_hat <- m_y
      y_hat[, -x_idx] <- yhat
      # y_hat <- cbind(yhat, m_y[,x_idx])
      y_hatc <- t(t(y_hat) - yhat_bar)
      idx1 <- lav_matrix_vech_col_idx(p)
      idx2 <- lav_matrix_vech_row_idx(p)
      if (meanstructure) {
        z <- (cbind(m_y, m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE]) -
          cbind(y_hat, y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE]))
        sigma1 <- c(m_mu, lav_matrix_vech(m_sigma))
        sigma2 <- c(yhat_bar, lav_matrix_vech(yhat_cov))
      } else {
        z <- (m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE] -
          y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE])
        sigma1 <- lav_matrix_vech(m_sigma)
        sigma2 <- lav_matrix_vech(yhat_cov)
      }
      zc <- t(t(z) - (sigma1 - sigma2))
    } else {
      m_qr <- qr(cbind(1, m_y[, x_idx, drop = FALSE]))
      yhat <- qr.fitted(m_qr, m_y[, -x_idx, drop = FALSE])
      # y_hat <- cbind(yhat, m_y[,x_idx])
      y_hat <- m_y
      y_hat[, -x_idx] <- yhat

      m_yc <- t(t(m_y) - colMeans(m_y, na.rm = TRUE))
      y_hatc <- t(t(y_hat) - colMeans(y_hat, na.rm = TRUE))
      idx1 <- lav_matrix_vech_col_idx(p)
      idx2 <- lav_matrix_vech_row_idx(p)
      if (meanstructure) {
        z <- (cbind(m_y, m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE]) -
          cbind(y_hat, y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE]))
      } else {
        z <- (m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE] -
          y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE])
      }
      zc <- t(t(z) - colMeans(z, na.rm = TRUE))
    }

    # clustered?
    if (length(cluster_idx) > 0L) {
      zc <- rowsum(zc, cluster_idx)
    }

    if (anyNA(zc)) {
      m_gamma <- lav_matrix_crossprod(zc) / n
    } else {
      m_gamma <- base::crossprod(zc) / n
    }
  } else {
    # conditional_x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    # regress m_y on m_x, and compute residuals
    m_x <- cbind(1, m_y[, x_idx, drop = FALSE])
    m_qr <- qr(m_x)
    m_res <- qr.resid(m_qr, m_y[, -x_idx, drop = FALSE])
    p <- ncol(m_res)

    idx1 <- lav_matrix_vech_col_idx(p)
    idx2 <- lav_matrix_vech_row_idx(p)

    if (meanstructure || slopestructure) {
      xtx_inv <- unname(solve(crossprod(m_x)))
      xi <- (m_x %*% xtx_inv) * n ## FIXME, shorter way?
      nc_x <- NCOL(m_x)
      nc_y <- NCOL(m_res)
    }

    if (meanstructure) {
      if (slopestructure) {
        # xi_idx <- rep(seq_len(nc_x), each  = nc_y)
        # res_idx <- rep(seq_len(nc_y), times = nc_x)
        xi_idx <- rep(seq_len(nc_x), times = nc_y)
        res_idx <- rep(seq_len(nc_y), each = nc_x)
        z <- cbind(
          xi[, xi_idx, drop = FALSE] *
            m_res[, res_idx, drop = FALSE],
          m_res[, idx1, drop = FALSE] *
            m_res[, idx2, drop = FALSE]
        )
      } else {
        xi_idx <- rep(1L, each = nc_y)
        z <- cbind(
          xi[, xi_idx, drop = FALSE] *
            m_res,
          m_res[, idx1, drop = FALSE] *
            m_res[, idx2, drop = FALSE]
        )
      }
    } else {
      if (slopestructure) {
        # xi_idx <- rep(seq_len(nc_x), each  = nc_y)
        # xi_idx <- xi_idx[ -seq_len(nc_y) ]
        xi_idx <- rep(seq(2, nc_x), times = nc_y)
        # res_idx <- rep(seq_len(nc_y), times = (nc_x - 1L))
        res_idx <- rep(seq_len(nc_y), each = (nc_x - 1L))
        z <- cbind(
          xi[, xi_idx, drop = FALSE] *
            m_res[, res_idx, drop = FALSE],
          m_res[, idx1, drop = FALSE] *
            m_res[, idx2, drop = FALSE]
        )
      } else {
        z <- m_res[, idx1, drop = FALSE] * m_res[, idx2, drop = FALSE]
      }
    }

    if (model_based) {
      zc <- t(t(z) - sigma)
    } else {
      zc <- t(t(z) - colMeans(z, na.rm = TRUE))
    }

    # clustered?
    if (length(cluster_idx) > 0L) {
      zc <- rowsum(zc, cluster_idx)
    }

    if (anyNA(zc)) {
      m_gamma <- lav_matrix_crossprod(zc) / n
    } else {
      m_gamma <- base::crossprod(zc) / n
    }
  }


  # only to mimic Mplus when estimator = "WLS"
  if (mplus_wls && !fixed_x && !conditional_x) {
    # adjust G_22 (the varcov part)
    m_s <- cov(m_y, use = "pairwise")
    w <- lav_matrix_vech(m_s)
    w_biased <- (n - 1) / n * w
    diff <- outer(w, w) - outer(w_biased, w_biased)
    if (meanstructure) {
      m_gamma[-seq_len(p), -seq_len(p)] <-
        m_gamma[-seq_len(p), -seq_len(p), drop = FALSE] - diff
    } else {
      m_gamma <- m_gamma - diff
    }

    if (meanstructure) {
      # adjust G_12/G_21 (third-order)
      # strange rescaling?
      n1 <- (n - 1) / n
      m_gamma[seq_len(p), -seq_len(p)] <- m_gamma[seq_len(p), -seq_len(p)] * n1
      m_gamma[-seq_len(p), seq_len(p)] <- m_gamma[-seq_len(p), seq_len(p)] * n1
    }
  }

  # clustered?
  if (length(cluster_idx) > 0L) {
    nc <- nrow(zc)
    m_gamma <- m_gamma * nc / (nc - 1)
  }

  # unbiased?
  if (unbiased) {
    # normal-theory m_gamma (cov only)
    # if (lav_use_lavaanC()) {
    #   gammant_cov <- lavaanC::m_kronecker_dup_ginv_pre_post(m_cov,
    #                                               multiplicator = 2.0)
    # } else {
      gammant_cov <- 2 * lav_matrix_duplication_ginv_pre_post(m_cov %x% m_cov)
    # }

    if (meanstructure) {
      gamma_cov <- m_gamma[-(1:p), -(1:p), drop = FALSE]
      gamma_mean_cov <- m_gamma[1:p, -(1:p), drop = FALSE]
    } else {
      gamma_cov <- m_gamma
    }

    # Browne's unbiased DF estimator (m_cov part)
    gamma_u <- (n * (n - 1) / (n - 2) / (n - 3) * gamma_cov -
      n / (n - 2) / (n - 3) * (gammant_cov -
        2 / (n - 1) * tcrossprod(cov_vech)))
    if (meanstructure) {
      m_gamma <- lav_matrix_bdiag(m_cov, gamma_u)

      # 3-rd order:
      m_gamma[1:p, (p + 1):ncol(m_gamma)] <- gamma_mean_cov * n / (n - 2)
      m_gamma[(p + 1):ncol(m_gamma), 1:p] <- t(gamma_mean_cov * n / (n - 2))
    } else {
      m_gamma <- gamma_u
    }
  } # unbiased

  m_gamma
}

# ADF Gamma for correlations
#
# 30 May 2024: basic version: fixed_x=FALSE, conditional_x=FALSE, ...
lav_samplestats_cor_gamma <- function(m_y, meanstructure = FALSE) {

  # coerce to matrix
  m_y <- unname(as.matrix(m_y))
  n <- nrow(m_y)
  p <- ncol(m_y)

  # compute m_s and m_r
  m_s <- cov(m_y) * (n - 1) / n
  m_r <- cov2cor(m_s)

  # create z-scores
  s_sd <- sqrt(diag(m_s))
  yz <- t((t(m_y) - colMeans(m_y)) / s_sd)

  # create squared z-scores
  yz2 <- yz * yz

  # find indices so we avoid 1) double subscripts (diagonal!), and
  #                          2) duplicated subscripts (symmetric!)
  idx1 <- lav_matrix_vech_col_idx(p, diagonal = FALSE)
  idx2 <- lav_matrix_vech_row_idx(p, diagonal = FALSE)

  zr1 <- (yz[, idx1, drop = FALSE] * yz[, idx2, drop = FALSE])
  zr2 <- (yz2[, idx1, drop = FALSE] + yz2[, idx2, drop = FALSE])
  zr2 <- t(t(zr2) * lav_matrix_vech(m_r, diagonal = FALSE))
  zrr <- zr1 - 0.5 * zr2
  if (meanstructure) {
      zrr <- cbind(yz, zrr)
  }
  m_gamma <- crossprod(zrr) / n

  m_gamma
}

# normal theory version
# 30 May 2024: basic version: fixed_x=FALSE, conditional_x=FALSE, ...
lav_samplestats_cor_gamma_nt <- function(m_y = NULL,
                                         wt = NULL,
                                         m_cov = NULL, # joint!
                                         m_mean = NULL, # joint!
                                         rescale = FALSE,
                                         x_idx = integer(0L),
                                         fixed_x = FALSE,
                                         conditional_x = FALSE,
                                         meanstructure = FALSE,
                                         slopestructure = FALSE) {
  # check arguments
  if (length(x_idx) == 0L) {
    conditional_x <- FALSE
    fixed_x <- FALSE
  } else {
    lav_msg_stop(gettext("x_idx not supported (yet) for correlations; use
                          fixed_x = FALSE (for now)"))
  }
  if (conditional_x) {
    lav_msg_stop(gettext("conditional_x = TRUE not supported (yet) for
                          correlations"))
  }

  # compute m_cov from m_y
  if (is.null(m_cov)) {
    stopifnot(!is.null(m_y))

    # coerce to matrix
    m_y <- unname(as.matrix(m_y))
    n <- nrow(m_y)
    if (is.null(wt)) {
      m_cov <- cov(m_y)
      if (rescale) {
        m_cov <- m_cov * (n - 1) / n # (normal) ML version
      }
    } else {
      out <- stats::cov.wt(m_y, wt = wt, method = "ML")
      m_cov <- out$cov
    }
  } else {
    if (!missing(rescale)) {
      lav_msg_warn(gettext("rescale= argument has no effect if m_cov is given"))
    }
    if (!missing(wt)) {
      lav_msg_warn(gettext("wt= argument has no effect if m_cov is given"))
    }
  }

  # if needed, compute m_mean from m_y
  if (conditional_x && length(x_idx) > 0L && is.null(m_mean) &&
    (meanstructure || slopestructure)) {
    stopifnot(!is.null(m_y))
    if (is.null(wt)) {
      m_mean <- colMeans(m_y, na.rm = TRUE)
    } else {
      m_mean <- out$center
    }
  }

  # rename
  m_s <- m_cov
  m_r <- cov2cor(m_s)
  s_p <- nrow(m_s)

  # unconditional
  if (!conditional_x) {
    # unconditional - stochastic x
    if (!fixed_x) {
      m_ip <- diag(s_p) %x% m_r
      m_rr <- m_r %x% m_r
      gamma_z_nt <- m_rr + lav_matrix_commutation_pre(m_rr)
      tmp <- (m_ip + lav_matrix_commutation_pre(m_ip)) / 2
      zero_idx <- seq_len(s_p * s_p)[-lav_matrix_diag_idx(s_p)]
      tmp[, zero_idx] <- 0
      m_a <- -tmp
      diag(m_a) <- 1 - diag(tmp)
      gamma_nt_big <- m_a %*% gamma_z_nt %*% t(m_a)
      r_idx <- lav_matrix_vech_idx(s_p, diagonal = FALSE)
      v_gamma <- gamma_nt_big[r_idx, r_idx]

      if (meanstructure) {
        v_gamma <- lav_matrix_bdiag(m_r, v_gamma)
      }

      # unconditional - fixed x
    } else {
      # TODO
    }
  } else {
    # conditional_x
    # TODO
  }

  v_gamma
}
