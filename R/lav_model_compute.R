lav_model_group_mm_indices <- function(nmat) {
  nmat_cumsum <- cumsum(c(0L, nmat))
  lapply(seq_along(nmat), function(g) {
    seq_len(nmat[g]) + nmat_cumsum[g]
  })
}

lav_model_get_mm_idx <- function(lavmodel) {
  if (.hasSlot(lavmodel, "mm.idx") && length(lavmodel@mm.idx) > 0L) {
    return(lavmodel@mm.idx)
  }

  lav_model_group_mm_indices(lavmodel@nmat)
}

lav_model_sigma_extra <- function(sigma_hat_g = NULL,
                                   check_sigma_pd = "chol") {
  # check if matrix is positive definite
  if (check_sigma_pd == "chol") {
    # fast path: try Cholesky directly
    c_s <- tryCatch(chol(sigma_hat_g), error = function(e) NULL)
    is_pd <- !is.null(c_s)
  } else {
    # slow path: eigenvalue-based PD check
    ev <- eigen(sigma_hat_g, symmetric = TRUE, only.values = TRUE)$values
    is_pd <- !(any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0)
  }

  if (!is_pd) {
    sigma_hat_inv <- MASS::ginv(sigma_hat_g)
    sigma_hat_log_det <- log(.Machine$double.eps)
    attr(sigma_hat_g, "po") <- FALSE
    attr(sigma_hat_g, "inv") <- sigma_hat_inv
    attr(sigma_hat_g, "log.det") <- sigma_hat_log_det
  } else {
    if (check_sigma_pd != "chol") {
      c_s <- chol(sigma_hat_g)
    }
    sigma_hat_inv <- chol2inv(c_s)
    d <- diag(c_s)
    sigma_hat_log_det <- 2 * sum(log(d))
    attr(sigma_hat_g, "po") <- TRUE
    attr(sigma_hat_g, "inv") <- sigma_hat_inv
    attr(sigma_hat_g, "log.det") <- sigma_hat_log_det
  }

  sigma_hat_g
}

lav_model_implied_fast <- function(lavmodel = NULL, glist = NULL,
                                   need_sigma = FALSE,
                                   need_mu = FALSE,
                                   need_th = FALSE,
                                   need_pi = FALSE,
                                   extra = FALSE,
                                   delta = TRUE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  # check.sigma.pd -- new in 0.6-21
  check_sigma_pd <- get0("opt_check_sigma_pd", lavaan_cache_env,
                         ifnotfound = "chol")

  nblocks <- lavmodel@nblocks
  representation <- lavmodel@representation
  meanstructure <- lavmodel@meanstructure
  categorical <- lavmodel@categorical
  conditional_x <- lavmodel@conditional.x
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  out <- list()
  if (need_sigma) out$sigma <- vector("list", length = nblocks)
  if (need_mu) out$mu <- vector("list", length = nblocks)
  if (need_th) out$th <- vector("list", length = nblocks)
  if (need_pi) out$pi <- vector("list", length = nblocks)

  for (g in seq_len(nblocks)) {
    mlist <- glist[mm_idx[[g]]]

    block_need_mu <- need_mu && meanstructure
    block_need_th <- need_th && categorical && length(lavmodel@th.idx[[g]]) > 0L
    block_need_pi <- need_pi && conditional_x

    if (representation == "LISREL") {
      block <- lav_lisrel_implied_fast(
        mlist = mlist,
        th_idx = lavmodel@th.idx[[g]],
        need_sigma = need_sigma,
        need_mu = block_need_mu,
        need_th = block_need_th,
        need_pi = block_need_pi,
        delta = delta
      )
    } else if (representation == "RAM") {
      block <- lav_ram_implied_fast(
        mlist = mlist,
        need_sigma = need_sigma,
        need_mu = block_need_mu,
        delta = delta
      )
      if (block_need_th || block_need_pi) {
        lav_msg_stop(gettext(
          "only LISREL representation supports thresholds and slopes"))
      }
    } else {
      lav_msg_stop(gettext(
        "only LISREL and RAM representation has been implemented for now"))
    }

    if (need_sigma) {
      out$sigma[[g]] <- block$sigma
      if (lav_debug()) print(out$sigma[[g]])
      if (extra) {
        out$sigma[[g]] <- lav_model_sigma_extra(
          sigma_hat_g = out$sigma[[g]],
          check_sigma_pd = check_sigma_pd
        )
      }
    }

    if (need_mu) {
      if (meanstructure) {
        out$mu[[g]] <- block$mu
      } else {
        out$mu[[g]] <- numeric(lavmodel@nvar[g])
      }

      # new in 0.6-20: if a variable is ordinal, set its mean to zero
      # (even if NU is not all zero, as in a multiple group analysis with
      # group.equal = "thresholds")
      if (categorical) {
        ord_idx <- unique(lavmodel@th.idx[[g]][lavmodel@th.idx[[g]] > 0L])
        out$mu[[g]][ord_idx] <- 0
      }
    }

    if (need_th) {
      if (length(lavmodel@th.idx[[g]]) == 0L) {
        out$th[[g]] <- numeric(0L)
      } else {
        out$th[[g]] <- block$th
      }
    }

    if (need_pi) {
      if (conditional_x) {
        out$pi[[g]] <- block$pi
      } else {
        out$pi[[g]] <- numeric(lavmodel@nvar[g])
      }
    }
  } # nblocks

  out
}

lav_model_implied_fast_state <- function(lavmodel = NULL, glist = NULL,
                                         implied = NULL,
                                         need_sigma = FALSE,
                                         need_mu = FALSE,
                                         need_th = FALSE,
                                         need_pi = FALSE,
                                         extra = FALSE,
                                         delta = TRUE) {
  spec <- list(
    need_sigma = need_sigma,
    need_mu = need_mu,
    need_th = need_th,
    need_pi = need_pi,
    extra = extra,
    delta = delta
  )

  if (is.environment(implied)) {
    if (!is.null(implied$implied) &&
        identical(implied$implied_spec, spec)) {
      return(implied$implied)
    }

    implied$implied <- lav_model_implied_fast(
      lavmodel = lavmodel,
      glist = glist,
      need_sigma = need_sigma,
      need_mu = need_mu,
      need_th = need_th,
      need_pi = need_pi,
      extra = extra,
      delta = delta
    )
    implied$implied_spec <- spec
    return(implied$implied)
  }

  lav_model_implied_fast(
    lavmodel = lavmodel,
    glist = glist,
    need_sigma = need_sigma,
    need_mu = need_mu,
    need_th = need_th,
    need_pi = need_pi,
    extra = extra,
    delta = delta
  )
}

lav_model_sigma <- function(lavmodel = NULL, glist = NULL, extra = FALSE,
                            delta = TRUE) {
  lav_model_implied_fast(
    lavmodel = lavmodel,
    glist = glist,
    need_sigma = TRUE,
    extra = extra,
    delta = delta
  )$sigma
}

## only if conditional.x = TRUE
## compute the (larger) unconditional 'joint' covariance matrix (y,x)
##
## Sigma (Joint ) = [ (S11, S12),
##                    (S21, S22) ] where
##     S11 = Sigma.res + PI %*% cov.x %*% t(PI)
##     S12 = PI %*% cov.x
##     S21 = cov.x %*% t(PI)
##     S22 = cov.x
lav_model_cond2joint_sigma <- function(lavmodel = NULL, glist = NULL,
                                       extra = FALSE, delta = TRUE) {
  stopifnot(lavmodel@conditional.x)

  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  # check.sigma.pd -- new in 0.6-21
  check_sigma_pd <- get0("opt_check_sigma_pd", lavaan_cache_env,
                         ifnotfound = "chol")

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  sigma_hat <- vector("list", length = nblocks)

  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      res_sigma <- lav_lisrel_sigma(mlist = mlist, delta = delta)
      # res.int <- lav_lisrel_mu(MLIST = MLIST)
      res_slopes <- lav_lisrel_pi(mlist = mlist)
      s_xx <- mlist$cov.x

      s_yy <- res_sigma + res_slopes %*% s_xx %*% t(res_slopes)
      s_yx <- res_slopes %*% s_xx
      s_xy <- s_xx %*% t(res_slopes)

      sigma_hat[[g]] <- rbind(cbind(s_yy, s_yx), cbind(s_xy, s_xx))
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }
    if (lav_debug()) print(sigma_hat[[g]])

    if (extra) {
      sigma_hat[[g]] <- lav_model_sigma_extra(
        sigma_hat_g = sigma_hat[[g]],
        check_sigma_pd = check_sigma_pd
      )
    }
  } # nblocks

  sigma_hat
}

lav_model_mu <- function(lavmodel = NULL, glist = NULL) {
  lav_model_implied_fast(
    lavmodel = lavmodel,
    glist = glist,
    need_mu = TRUE
  )$mu
}

## only if conditional.x = TRUE
## compute the (larger) unconditional 'joint' mean vector (y,x)
##
## Mu (Joint ) = [ Mu.y, Mu.x ] where
##     Mu.y = res.int + PI %*% M.x
##     Mu.x = M.x
lav_model_cond2joint_mu <- function(lavmodel = NULL, glist = NULL) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  representation <- lavmodel@representation
  meanstructure <- lavmodel@meanstructure
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  mu_hat <- vector("list", length = nblocks)

  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]

    if (!meanstructure) {
      mu_hat[[g]] <- numeric(lavmodel@nvar[g])
    } else if (representation == "LISREL") {
      mlist <- glist[mm_in_group]
      res_int <- lav_lisrel_mu(mlist = mlist)
      res_slopes <- lav_lisrel_pi(mlist = mlist)
      m_x <- mlist$mean.x

      mu_y <- res_int + res_slopes %*% m_x
      mu_x <- m_x
      mu_hat[[g]] <- c(mu_y, mu_x)
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }
  } # nblocks

  mu_hat
}

# TH.star = DELTA.star * (th.star - pi0.star)
# see Muthen 1984 eq 11
lav_model_th <- function(lavmodel = NULL, glist = NULL, delta = TRUE) {
  lav_model_implied_fast(
    lavmodel = lavmodel,
    glist = glist,
    need_th = TRUE,
    delta = delta
  )$th
}

# PI = slope structure
# see Muthen 1984 eq 12
lav_model_pi <- function(lavmodel = NULL, glist = NULL, delta = TRUE) {
  lav_model_implied_fast(
    lavmodel = lavmodel,
    glist = glist,
    need_pi = TRUE,
    delta = delta
  )$pi
}


# GW = group weight
lav_model_gw <- function(lavmodel = NULL, glist = NULL) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  group_w_free <- lavmodel@group.w.free
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  gw <- vector("list", length = nblocks)

  # compute GW for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (!group_w_free) {
      gw_g <- 0.0 # FIXME
    } else if (representation == "LISREL") {
      gw_g <- as.numeric(mlist$gw[1, 1])
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    gw[[g]] <- gw_g
  }

  # transform to proportions
  # gw <- unlist(GW)
  # gw <- exp(gw) / sum(exp(gw))
  # for(g in 1:nblocks) {
  #    GW[[g]] <- gw[g]
  # }

  gw
}

# *unconditional* variance/covariance matrix of Y
#  - same as Sigma.hat if all Y are continuous)
#  - if also Gamma, cov.x is used (only if categorical)
lav_model_vy <- function(lavmodel = NULL, glist = NULL, diagonal_only = FALSE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  vy <- vector("list", length = nblocks)

  # compute TH for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      vy_g <- lav_lisrel_vy(mlist = mlist)
    } else if (representation == "RAM") {
      # does not work for categorical setting yet
      stopifnot(!lavmodel@categorical)
      # does not work if conditional.x = TRUE
      stopifnot(!lavmodel@conditional.x)
      vy_g <- lav_ram_sigmahat(mlist = mlist)
    } else {
      lav_msg_stop(gettext(
        "only RAM and LISREL representation has been implemented for now"))
    }

    if (diagonal_only) {
      vy[[g]] <- diag(vy_g)
    } else {
      vy[[g]] <- vy_g
    }
  }

  vy
}

# V(ETA): latent variances variances/covariances
lav_model_veta <- function(lavmodel = NULL, glist = NULL,
                           remove_dummy_lv = FALSE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  veta <- vector("list", length = nblocks)

  # compute VETA for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      veta_g <- lav_lisrel_veta(mlist = mlist)

      if (remove_dummy_lv) {
        # remove all dummy latent variables
        lv_idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (!is.null(lv_idx)) {
          veta_g <- veta_g[-lv_idx, -lv_idx, drop = FALSE]
        }
      }
    } else if (representation == "RAM") {
      veta_g <- lav_ram_veta(mlist = mlist)
    } else {
      lav_msg_stop(gettext(
        "only LISREL and RAM representation has been implemented for now"))
    }

    veta[[g]] <- veta_g
  }

  veta
}

# V(ETA|x_i): latent variances variances/covariances, conditional on x_
# - this is always (I-B)^-1 PSI (I-B)^-T, after REMOVING lv dummies
lav_model_vetax <- function(lavmodel = NULL, glist = NULL) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  eta <- vector("list", length = nblocks)

  # compute ETA for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      lv_idx <- c(
        lavmodel@ov.y.dummy.lv.idx[[g]],
        lavmodel@ov.x.dummy.lv.idx[[g]]
      )
      eta_g <- lav_lisrel_vetax(
        mlist = mlist,
        lv_dummy_idx = lv_idx
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    eta[[g]] <- eta_g
  }

  eta
}

# COV: observed+latent variances variances/covariances
lav_model_cov_both <- function(lavmodel = NULL, glist = NULL,
                               remove_dummy_lv = FALSE, delta = TRUE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  cov_1 <- vector("list", length = nblocks)

  # compute COV for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      cov_g <- lav_lisrel_cov_both(mlist = mlist, delta = delta)

      if (remove_dummy_lv) {
        # remove all dummy latent variables
        lv_idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (!is.null(lv_idx)) {
          # offset for ov
          lambda_names <-
            lavmodel@dimNames[[which(names(glist) == "lambda")[g]]][[1L]]
          lv_idx <- lv_idx + length(lambda_names)
          cov_g <- cov_g[-lv_idx, -lv_idx, drop = FALSE]
        }
      }
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    cov_1[[g]] <- cov_g
  }

  cov_1
}


# E(ETA): expectation (means) of latent variables (return vector)
lav_model_eeta <- function(lavmodel = NULL, glist = NULL, lavsamplestats = NULL,
                           remove_dummy_lv = FALSE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  eeta <- vector("list", length = nblocks)

  # compute E(ETA) for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      eeta_g <- lav_lisrel_eeta(mlist,
        mean_x = lavsamplestats@mean.x[[g]],
        sample_mean = lavsamplestats@mean[[g]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]]
      )
      if (remove_dummy_lv) {
        # remove dummy
        lv_dummy_idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (length(lv_dummy_idx) > 0L) {
          eeta_g <- eeta_g[-lv_dummy_idx]
        }
      }
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    eeta[[g]] <- eeta_g
  }

  eeta
}

# E(ETA|x_i): conditional expectation (means) of latent variables
# for a given value of x_i (instead of E(x_i))
lav_model_eetax <- function(lavmodel = NULL, glist = NULL,
                            lavsamplestats = NULL, exo = NULL,
                            nobs = NULL, remove_dummy_lv = FALSE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  eetax <- vector("list", length = nblocks)

  # compute E(ETA) for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    exo_1 <- exo[[g]]
    if (is.null(exo_1)) {
      # create empty matrix
      exo_1 <- matrix(0, nobs[[g]], 0L)
    }

    if (representation == "LISREL") {
      eetax_g <- lav_lisrel_eetax(mlist,
        exo = exo_1, n = nobs[[g]],
        sample_mean = lavsamplestats@mean[[g]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]]
      )

      if (remove_dummy_lv) {
        # remove dummy
        lv_dummy_idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (length(lv_dummy_idx) > 0L) {
          eetax_g <- eetax_g[, -lv_dummy_idx, drop = FALSE]
        }
      }
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    eetax[[g]] <- eetax_g
  }

  eetax
}

# return 'regular' LAMBDA
lav_model_lambda <- function(lavmodel = NULL, glist = NULL,
                             handle_dummy_lv = TRUE,
                             remove_dummy_lv = FALSE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  mm_lambda <- vector("list", length = nblocks)

  # compute LAMBDA for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      if (handle_dummy_lv) {
        ov_y_dummy_ov_idx <- lavmodel@ov.y.dummy.ov.idx[[g]]
        ov_x_dummy_ov_idx <- lavmodel@ov.x.dummy.ov.idx[[g]]
        ov_y_dummy_lv_idx <- lavmodel@ov.y.dummy.lv.idx[[g]]
        ov_x_dummy_lv_idx <- lavmodel@ov.x.dummy.lv.idx[[g]]
      } else {
        ov_y_dummy_ov_idx <- NULL
        ov_x_dummy_ov_idx <- NULL
        ov_y_dummy_lv_idx <- NULL
        ov_x_dummy_lv_idx <- NULL
      }
      lambda_g <- lav_lisrel_lambda(
        mlist = mlist,
        ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
        ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
        ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
        ov_x_dummy_lv_idx = ov_x_dummy_lv_idx,
        remove_dummy_lv = remove_dummy_lv
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    mm_lambda[[g]] <- lambda_g
  }

  mm_lambda
}

# THETA: observed (residual) variances
lav_model_theta <- function(lavmodel = NULL, glist = NULL, fix = TRUE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  mm_theta <- vector("list", length = nblocks)

  # compute THETA for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      if (fix) {
        theta_g <- lav_lisrel_theta(
          mlist = mlist,
          ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
          ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
          ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
          ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]]
        )
      } else {
        theta_g <- lav_lisrel_theta(mlist = mlist)
      }
    } else if (representation == "RAM") {
      ov_idx <- as.integer(mlist$ov.idx[1, ])
      theta_g <- mlist$S[ov_idx, ov_idx, drop = FALSE]
    } else {
      lav_msg_stop(gettext(
        "only LISREL and RAM representation has been implemented for now"))
    }

    mm_theta[[g]] <- theta_g
  }

  mm_theta
}

# E(Y): expectation (mean) of observed variables
# returns vector 1 x nvar
lav_model_ey <- function(lavmodel = NULL, glist = NULL, lavsamplestats = NULL,
                         delta = TRUE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  ey <- vector("list", length = nblocks)

  # compute E(Y) for each group
  for (g in seq_len(nblocks)) {
    # which mm belong to group g?
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (representation == "LISREL") {
      ey_g <- lav_lisrel_ey(
        mlist = mlist,
        mean_x = lavsamplestats@mean.x[[g]],
        sample_mean = lavsamplestats@mean[[g]],
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        delta = delta
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    ey[[g]] <- ey_g
  }

  ey
}


# E(Y | ETA, x_i): conditional expectation (means) of observed variables
# for a given value of x_i AND eta_i
lav_model_yhat <- function(lavmodel = NULL, glist = NULL, lavsamplestats = NULL,
                           exo = NULL, nobs = NULL, eta = NULL,
                           duplicate = FALSE, delta = TRUE) {
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST

  # ngroups, not nblocks!
  ngroups <- lavsamplestats@ngroups
  mm_idx <- lav_model_get_mm_idx(lavmodel)

  # return a list
  yhat <- vector("list", length = ngroups)

  # compute YHAT for each group
  for (g in seq_len(ngroups)) {
    # which mm belong to group g?
    # FIXME: what if more than g blocks???
    mm_in_group <- mm_idx[[g]]
    mlist <- glist[mm_in_group]

    if (is.null(exo[[g]]) && duplicate) {
      nobs_1 <- nobs[[g]]
    } else {
      nobs_1 <- 1L
    }

    if (lavmodel@representation == "LISREL") {
      if (lavmodel@conditional.x) {
        yhat[[g]] <- lav_lisrel_eyetax(
          mlist = mlist,
          exo = exo[[g]], eta = eta[[g]], n = nobs_1,
          sample_mean = lavsamplestats@mean[[g]],
          ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
          ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
          ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
          ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
          delta = delta
        )
      } else {
        # unconditional case
        yhat[[g]] <- lav_lisrel_eyetax3(
          mlist = mlist,
          eta = eta[[g]],
          sample_mean = lavsamplestats@mean[[g]],
          mean_x = lavsamplestats@mean.x[[g]],
          ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
          ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
          ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
          ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
          delta = delta
        )
        # impute back ov.y values that are NOT indicators
      }
    } else {
      lav_msg_stop(gettextf("representation %s not supported yet.",
                            lavmodel@representation))
    }
  }

  yhat
}
