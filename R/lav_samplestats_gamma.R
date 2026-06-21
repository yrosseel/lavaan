


# for internal use -- lavobject or internal slots
lav_object_gamma <- function(lavobject = NULL,
                             # or individual slots
                             lavdata = NULL,
                             lavoptions = NULL,
                             lavsamplestats = NULL,
                             lavh1 = NULL,
                             lavimplied = NULL,
                             # other options
                             adf = TRUE, model_based = FALSE,
                             mplus_wls = FALSE) {
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
    model_based <- TRUE
  }
  fixed_x <- lavoptions$fixed.x
  conditional_x <- lavoptions$conditional.x
  meanstructure <- lavoptions$meanstructure
  gamma_n_minus_one <- lavoptions$gamma.n.minus.one
  gamma_unbiased <- lavoptions$gamma.unbiased

  if (adf && model_based && conditional_x) {
    lav_msg_stop(gettext(
      "ADF + model.based + conditional.x is not supported yet."))
  }

  # output container
  out <- vector("list", length = lavdata@ngroups)

  # compute Gamma matrix for each group
  for (g in seq_len(lavdata@ngroups)) {
    x_idx <- lavsamplestats@x.idx[[g]]
    cov_1 <- mean_1 <- NULL
    if (!adf || model_based) {
      implied <- lavh1$implied # saturated/unstructured
      if (model_based) {
        implied <- lavimplied # model-based/structured
      }
      if (conditional_x) {
        # convert to joint COV/MEAN
        res_s <- implied$res.cov[[g]]
        res_slopes <- implied$res.slopes[[g]]
        res_int <- implied$res.int[[g]]
        s_xx <- implied$cov.x[[g]]
        m_x <- implied$mean.x[[g]]

        s_yy <- res_s + res_slopes %*% s_xx %*% t(res_slopes)
        s_yx <- res_slopes %*% s_xx
        s_xy <- s_xx %*% t(res_slopes)
        m_y <- res_int + res_slopes %*% m_x

        cov_1 <- rbind(cbind(s_yy, s_yx), cbind(s_xy, s_xx))
        mean_1 <- c(m_y, m_x)
      } else {
        # not conditional.x
        cov_1 <- implied$cov[[g]]
        mean_1 <- implied$mean[[g]]
      }
    } # COV/MEAN

    if (adf) {
      if (conditional_x) {
        y <- cbind(lavdata@X[[g]], lavdata@eXo[[g]])
      } else {
        y <- lavdata@X[[g]]
      }
      if (length(lavdata@cluster) > 0L) {
        cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster_idx <- NULL
      }
      out[[g]] <- lav_samp_gamma(
        m_y = y,
        m_mu = mean_1,
        m_sigma = cov_1,
        x_idx = x_idx,
        cluster_idx = cluster_idx,
        fixed_x = fixed_x,
        conditional_x = conditional_x,
        meanstructure = meanstructure,
        slopestructure = conditional_x,
        gamma_n_minus_one = gamma_n_minus_one,
        unbiased = gamma_unbiased,
        mplus_wls = mplus_wls
      )
    } else {
      out[[g]] <- lav_samp_gamma_nt(
        m_cov = cov_1, # joint!
        m_mean = mean_1, # joint!
        x_idx = x_idx,
        fixed_x = fixed_x,
        conditional_x = conditional_x,
        meanstructure = meanstructure,
        slopestructure = conditional_x
      )
    }

    # group.w.free
    if (lavoptions$group.w.free) {
      # checkme!
      out[[g]] <- lav_mat_bdiag(matrix(1, 1, 1), out[[g]])
    }

  } # g

  out
}





# NOTE:
#  - three types:
#       1) plain         (conditional.x = FALSE, fixed.x = FALSE)
#       2) fixed.x       (conditional.x = FALSE, fixed.x = TRUE)
#       3) conditional.x (conditional.x = TRUE)
#  - if conditional.x = TRUE, we ignore fixed.x (can be TRUE or FALSE)

# NORMAL-THEORY
lav_samp_gamma_nt <- function(m_y = NULL, # should include
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
        m_gamma <- 2 * lav_mat_dup_ginv_pre_post(m_s %x% m_s)
      # }
      if (meanstructure) {
        m_gamma <- lav_mat_bdiag(m_s, m_gamma)
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
        gamma_s <- 2 * lav_mat_dup_ginv_pre_post(m_s %x% m_s)
        gamma_r <- 2 * lav_mat_dup_ginv_pre_post(m_r %x% m_r)
      # }
      m_gamma <- gamma_s - gamma_r

      if (meanstructure) {
        m_gamma <- lav_mat_bdiag(ybarx_aug, m_gamma)
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
               lav_mat_dup_ginv_pre_post(cov_ybarx %x% cov_ybarx)
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

    m_gamma <- lav_mat_bdiag(a11, m_gamma)
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
lav_samp_gamma <- function(m_y, # Y+X if cond!
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
                                  mplus_wls = FALSE,
                                  wt = NULL,
                                  sampling_weights_type = "design") {
  # coerce to matrix
  m_y <- unname(as.matrix(m_y))
  n <- nrow(m_y)
  p <- ncol(m_y)

  # sampling weights?
  # Gamma is the asymptotic (co)variance of the sample moments, which for
  # the saturated (h1) model equals the sandwich I^{-1} J I^{-1}, where I is
  # the expected and J the first-order h1 information matrix. We reuse the
  # (weight-aware) expected information I, and construct the meat J according
  # to how the sampling weights are interpreted (sampling.weights.type):
  #  - "design"    : J = crossprod(wt * sc) / sum(wt), i.e. sum(wt^2)-weighted,
  #                  the design-based sandwich; the parameter vcov is then
  #                  invariant to the overall scale of the weights, matching
  #                  the continuous ML robust (huber.white) SEs.
  #  - "frequency" : J = crossprod(sqrt(wt) * sc) / sum(wt), i.e. sum(wt)-
  #                  weighted; treats the weights as frequencies, so the result
  #                  equals the Gamma from physically replicating each case wt
  #                  times. (NB: this also feeds the WLS/DWLS/DLS weight matrix,
  #                  so it shifts those point estimates relative to "design".)
  # Under normality J = I and this reduces to the normal-theory Gamma.
  if (!is.null(wt)) {
    if (conditional_x) {
      lav_msg_stop(gettext(
        "Gamma with sampling weights is not available (yet) if
         conditional.x = TRUE."))
    }
    if (length(cluster_idx) > 0L) {
      lav_msg_stop(gettext(
        "Gamma with sampling weights is not available (yet) for
         clustered data."))
    }
    out <- stats::cov.wt(m_y, wt = wt, method = "ML")
    i_mat <- lav_mvn_h1_info_expected(
      y = m_y, wt = wt, sample_cov = out$cov, x_idx = x_idx,
      meanstructure = meanstructure
    )
    # unit (unweighted) casewise scores of the saturated model, evaluated
    # at the weighted sample moments
    sc <- lav_mvn_sc_mu_sigma(
      y = m_y, mu = out$center, sigma_1 = out$cov, x_idx = x_idx
    )
    if (!meanstructure) {
      sc <- sc[, -seq_len(p), drop = FALSE]
    }
    if (identical(sampling_weights_type, "frequency")) {
      j_mat <- crossprod(sqrt(wt) * sc) / sum(wt)
    } else {
      j_mat <- crossprod(wt * sc) / sum(wt)
    }
    # fixed.x zeroes the x-blocks of I, so use a pseudo-inverse there
    if (length(x_idx) > 0L) {
      i_inv <- MASS::ginv(i_mat)
    } else {
      i_inv <- solve(i_mat)
    }
    return(i_inv %*% j_mat %*% i_inv)
  }

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
      cov_vech <- lav_mat_vech(m_cov)
    }
  }

  # model-based?
  if (!is.null(m_sigma)) {
    stopifnot(!conditional_x)
    model_based <- TRUE
    if (meanstructure) {
      stopifnot(!is.null(m_mu))
      sigma <- c(as.numeric(m_mu), lav_mat_vech(m_sigma))
    } else {
      m_mu <- colMeans(m_y, na.rm = TRUE) # for centering!
      sigma <- lav_mat_vech(m_sigma)
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
    idx1 <- lav_mat_vech_col_idx(p)
    idx2 <- lav_mat_vech_row_idx(p)
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
        sigma <- c(as.numeric(m_mu), lav_mat_vech(m_sigma))
      } else {
        sigma <- lav_mat_vech(m_sigma)
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
      m_gamma <- lav_mat_crossprod(zc) / n
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
      idx1 <- lav_mat_vech_col_idx(p)
      idx2 <- lav_mat_vech_row_idx(p)
      if (meanstructure) {
        z <- (cbind(m_y, m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE]) -
          cbind(y_hat, y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE]))
        sigma1 <- c(m_mu, lav_mat_vech(m_sigma))
        sigma2 <- c(yhat_bar, lav_mat_vech(yhat_cov))
      } else {
        z <- (m_yc[, idx1, drop = FALSE] *
          m_yc[, idx2, drop = FALSE] -
          y_hatc[, idx1, drop = FALSE] *
            y_hatc[, idx2, drop = FALSE])
        sigma1 <- lav_mat_vech(m_sigma)
        sigma2 <- lav_mat_vech(yhat_cov)
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
      idx1 <- lav_mat_vech_col_idx(p)
      idx2 <- lav_mat_vech_row_idx(p)
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
      m_gamma <- lav_mat_crossprod(zc) / n
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

    idx1 <- lav_mat_vech_col_idx(p)
    idx2 <- lav_mat_vech_row_idx(p)

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
      m_gamma <- lav_mat_crossprod(zc) / n
    } else {
      m_gamma <- base::crossprod(zc) / n
    }
  }


  # only to mimic Mplus when estimator = "WLS"
  if (mplus_wls && !fixed_x && !conditional_x) {
    # adjust G_22 (the varcov part)
    m_s <- cov(m_y, use = "pairwise")
    w <- lav_mat_vech(m_s)
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
      gammant_cov <- 2 * lav_mat_dup_ginv_pre_post(m_cov %x% m_cov)
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
      m_gamma <- lav_mat_bdiag(m_cov, gamma_u)

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
# ADF (sandwich) Gamma for a *partial* correlation structure: only the
# variables in 'cor_idx' are standardized to unit variance, the others keep
# their original metric. The partial moment vector is
#   [ var(num_idx), off-diagonal (co)variances ]   (and the means, if any)
# where num_idx = complement of cor_idx. Each partial moment is a smooth
# function of the raw moments [mean, vech(S)], so the ADF Gamma follows from
# the delta method: Gamma_partial = J %*% Gamma_full %*% t(J), with J the
# Jacobian of the partial moments wrt the raw moments. (For cor_idx = all
# variables this reduces to lav_samp_cor_gamma().)
lav_samp_partial_cor_gamma <- function(m_y, cor_idx = NULL,
                                       meanstructure = FALSE) {
  m_y <- unname(as.matrix(m_y))
  n <- nrow(m_y)
  p <- ncol(m_y)
  m_s <- cov(m_y) * (n - 1) / n

  is_cor <- logical(p)
  is_cor[cor_idx] <- TRUE
  num_idx <- which(!is_cor)
  # scaling: sd for correlation variables, 1 for the others
  d <- rep(1.0, p)
  d[is_cor] <- sqrt(diag(m_s)[is_cor])

  # raw ADF Gamma over [mean(p)?, vech(S)]
  gamma_full <- lav_samp_gamma(m_y, meanstructure = meanstructure)

  # vech bookkeeping (full, with diagonal)
  vr <- lav_mat_vech_row_idx(p)
  vc <- lav_mat_vech_col_idx(p)
  pstar <- length(vr)
  diag_pos <- which(vr == vc) # variances, in variable order
  off_pos <- which(vr != vc) # covariances
  var_pos <- diag_pos[num_idx] # variances we keep (non-correlation vars)

  # Jacobian (covariance part): rows = [var(num_idx), offdiag], cols = vech(S)
  nkeep <- length(var_pos) + length(off_pos)
  j_cov <- matrix(0, nkeep, pstar)
  if (length(var_pos) > 0L) {
    j_cov[cbind(seq_along(var_pos), var_pos)] <- 1
  }
  off0 <- length(var_pos)
  for (i in seq_along(off_pos)) {
    pos <- off_pos[i]
    k <- vr[pos]
    l <- vc[pos]
    m_kl <- m_s[k, l] / (d[k] * d[l])
    row <- off0 + i
    j_cov[row, pos] <- 1 / (d[k] * d[l])
    if (is_cor[k]) {
      j_cov[row, diag_pos[k]] <- j_cov[row, diag_pos[k]] -
        0.5 * m_kl / m_s[k, k]
    }
    if (is_cor[l]) {
      j_cov[row, diag_pos[l]] <- j_cov[row, diag_pos[l]] -
        0.5 * m_kl / m_s[l, l]
    }
  }

  if (meanstructure) {
    # partial moments: [mean(p), var(num_idx), offdiag]; the mean of a
    # standardized variable is scaled by 1/sd (first order), others by 1
    j_mat <- matrix(0, p + nkeep, p + pstar)
    j_mat[cbind(seq_len(p), seq_len(p))] <- 1 / d
    j_mat[(p + 1L):(p + nkeep), (p + 1L):(p + pstar)] <- j_cov
  } else {
    j_mat <- j_cov
  }

  j_mat %*% gamma_full %*% t(j_mat)
}

lav_samp_cor_gamma <- function(m_y, meanstructure = FALSE) {

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
  idx1 <- lav_mat_vech_col_idx(p, diagonal = FALSE)
  idx2 <- lav_mat_vech_row_idx(p, diagonal = FALSE)

  zr1 <- (yz[, idx1, drop = FALSE] * yz[, idx2, drop = FALSE])
  zr2 <- (yz2[, idx1, drop = FALSE] + yz2[, idx2, drop = FALSE])
  zr2 <- t(t(zr2) * lav_mat_vech(m_r, diagonal = FALSE))
  zrr <- zr1 - 0.5 * zr2
  if (meanstructure) {
      zrr <- cbind(yz, zrr)
  }
  m_gamma <- crossprod(zrr) / n

  m_gamma
}

# normal theory version
# 30 May 2024: basic version: fixed_x=FALSE, conditional_x=FALSE, ...
lav_samp_cor_gamma_nt <- function(m_y = NULL,
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
      gamma_z_nt <- m_rr + lav_mat_com_pre(m_rr)
      tmp <- (m_ip + lav_mat_com_pre(m_ip)) / 2
      zero_idx <- seq_len(s_p * s_p)[-lav_mat_diag_idx(s_p)]
      tmp[, zero_idx] <- 0
      m_a <- -tmp
      diag(m_a) <- 1 - diag(tmp)
      gamma_nt_big <- m_a %*% gamma_z_nt %*% t(m_a)
      r_idx <- lav_mat_vech_idx(s_p, diagonal = FALSE)
      v_gamma <- gamma_nt_big[r_idx, r_idx, drop = FALSE]

      if (meanstructure) {
        v_gamma <- lav_mat_bdiag(m_r, v_gamma)
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
