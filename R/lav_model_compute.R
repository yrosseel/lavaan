lav_model_sigma <- function(lavmodel = NULL, GLIST = NULL, extra = FALSE,
                            delta = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  # check.sigma.pd -- new in 0.6-21
  check.sigma.pd <- get0("opt.check.sigma.pd", lavaan_cache_env,
                         ifnotfound = "chol")

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  representation <- lavmodel@representation

  # return a list
  Sigma.hat <- vector("list", length = nblocks)

  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      Sigma.hat[[g]] <- lav_lisrel_sigma(
        MLIST = MLIST,
        delta = delta
      )
    } else if (representation == "RAM") {
      Sigma.hat[[g]] <- lav_ram_sigmahat(MLIST = MLIST, delta = delta)
    } else {
      lav_msg_stop(gettext(
        "only LISREL and RAM representation has been implemented for now"))
    }
    if (lav_debug()) print(Sigma.hat[[g]])

    if (extra) {
      # check if matrix is positive definite
      if (check.sigma.pd == "chol") {
        # fast path: try Cholesky directly
        cS <- tryCatch(chol(Sigma.hat[[g]]), error = function(e) NULL)
        is_pd <- !is.null(cS)
      } else {
        # slow path: eigenvalue-based PD check
        ev <- eigen(Sigma.hat[[g]], symmetric = TRUE, only.values = TRUE)$values
        is_pd <- !(any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0)
      }
      if (!is_pd) {
        Sigma.hat.inv <- MASS::ginv(Sigma.hat[[g]])
        Sigma.hat.log.det <- log(.Machine$double.eps)
        attr(Sigma.hat[[g]], "po") <- FALSE
        attr(Sigma.hat[[g]], "inv") <- Sigma.hat.inv
        attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
      } else {
        if (check.sigma.pd != "chol") {
          cS <- chol(Sigma.hat[[g]])
        }
        Sigma.hat.inv <- chol2inv(cS)
        d <- diag(cS)
        Sigma.hat.log.det <- 2 * sum(log(d))
        attr(Sigma.hat[[g]], "po") <- TRUE
        attr(Sigma.hat[[g]], "inv") <- Sigma.hat.inv
        attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
      }
    }
  } # nblocks
  Sigma.hat
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
lav_model_cond2joint_sigma <- function(lavmodel = NULL, GLIST = NULL,
                                       extra = FALSE, delta = TRUE) {
  stopifnot(lavmodel@conditional.x)

  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  # check.sigma.pd -- new in 0.6-21
  check.sigma.pd <- get0("opt.check.sigma.pd", lavaan_cache_env,
                         ifnotfound = "chol")

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  representation <- lavmodel@representation

  # return a list
  Sigma.hat <- vector("list", length = nblocks)

  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      res.Sigma <- lav_lisrel_sigma(MLIST = MLIST, delta = delta)
      # res.int <- lav_lisrel_mu(MLIST = MLIST)
      res.slopes <- lav_lisrel_pi(MLIST = MLIST)
      S.xx <- MLIST$cov.x

      S.yy <- res.Sigma + res.slopes %*% S.xx %*% t(res.slopes)
      S.yx <- res.slopes %*% S.xx
      S.xy <- S.xx %*% t(res.slopes)

      Sigma.hat[[g]] <- rbind(cbind(S.yy, S.yx), cbind(S.xy, S.xx))
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }
    if (lav_debug()) print(Sigma.hat[[g]])

    if (extra) {
      # check if matrix is positive definite
      if (check.sigma.pd == "chol") {
        # fast path: try Cholesky directly
        cS <- tryCatch(chol(Sigma.hat[[g]]), error = function(e) NULL)
        is_pd <- !is.null(cS)
      } else {
        # slow path: eigenvalue-based PD check
        ev <- eigen(Sigma.hat[[g]], symmetric = TRUE, only.values = TRUE)$values
        is_pd <- !(any(ev < sqrt(.Machine$double.eps)) || sum(ev) == 0)
      }
      if (!is_pd) {
        Sigma.hat.inv <- MASS::ginv(Sigma.hat[[g]])
        Sigma.hat.log.det <- log(.Machine$double.eps)
        attr(Sigma.hat[[g]], "po") <- FALSE
        attr(Sigma.hat[[g]], "inv") <- Sigma.hat.inv
        attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
      } else {
        if (check.sigma.pd != "chol") {
          cS <- chol(Sigma.hat[[g]])
        }
        Sigma.hat.inv <- chol2inv(cS)
        d <- diag(cS)
        Sigma.hat.log.det <- 2 * sum(log(d))
        attr(Sigma.hat[[g]], "po") <- TRUE
        attr(Sigma.hat[[g]], "inv") <- Sigma.hat.inv
        attr(Sigma.hat[[g]], "log.det") <- Sigma.hat.log.det
      }
    }
  } # nblocks

  Sigma.hat
}

lav_model_mu <- function(lavmodel = NULL, GLIST = NULL) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  representation <- lavmodel@representation
  meanstructure <- lavmodel@meanstructure

  # return a list
  Mu.hat <- vector("list", length = nblocks)

  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (!meanstructure) {
      Mu.hat[[g]] <- numeric(lavmodel@nvar[g])
    } else if (representation == "LISREL") {
      Mu.hat[[g]] <- lav_lisrel_mu(MLIST = MLIST)
    } else if (representation == "RAM") {
      Mu.hat[[g]] <- lav_ram_muhat(MLIST = MLIST)
    } else {
      lav_msg_stop(gettext(
        "only RAM and LISREL representation has been implemented for now"))
    }

    # new in 0.6-20: if a variable is ordinal, set its mean to zero
    # (even if NU is not all zero, as in a multiple group analysis with
    # group.equal = "thresholds")
    #
    # the logic is: Mu.hat is about 'y', not 'y-star'
    # the non-free intercepts (in TAU) are used when computing the
    # model-implied thresholds, but the do not say anything about the
    # 'observed' mean of 'y'
    if (lavmodel@categorical) {
      ord.idx <- unique(lavmodel@th.idx[[g]][lavmodel@th.idx[[g]] > 0L])
      Mu.hat[[g]][ord.idx] <- 0
    }
  } # nblocks

  Mu.hat
}

## only if conditional.x = TRUE
## compute the (larger) unconditional 'joint' mean vector (y,x)
##
## Mu (Joint ) = [ Mu.y, Mu.x ] where
##     Mu.y = res.int + PI %*% M.x
##     Mu.x = M.x
lav_model_cond2joint_mu <- function(lavmodel = NULL, GLIST = NULL) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  # check.sigma.pd -- new in 0.6-21
  check.sigma.pd <- get0("opt.check.sigma.pd", lavaan_cache_env,
                         ifnotfound = "chol")

  nmat <- lavmodel@nmat
  nblocks <- lavmodel@nblocks
  representation <- lavmodel@representation
  meanstructure <- lavmodel@meanstructure

  # return a list
  Mu.hat <- vector("list", length = nblocks)

  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]

    if (!meanstructure) {
      Mu.hat[[g]] <- numeric(lavmodel@nvar[g])
    } else if (representation == "LISREL") {
      MLIST <- GLIST[mm.in.group]
      res.int <- lav_lisrel_mu(MLIST = MLIST)
      res.slopes <- lav_lisrel_pi(MLIST = MLIST)
      M.x <- MLIST$mean.x

      Mu.y <- res.int + res.slopes %*% M.x
      Mu.x <- M.x
      Mu.hat[[g]] <- c(Mu.y, Mu.x)
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }
  } # nblocks

  Mu.hat
}

# TH.star = DELTA.star * (th.star - pi0.star)
# see Muthen 1984 eq 11
lav_model_th <- function(lavmodel = NULL, GLIST = NULL, delta = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  th.idx <- lavmodel@th.idx

  # return a list
  TH <- vector("list", length = nblocks)

  # compute TH for each group
  for (g in 1:nblocks) {
    if (length(th.idx[[g]]) == 0) {
      TH[[g]] <- numeric(0L)
      next
    }

    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]

    if (representation == "LISREL") {
      TH[[g]] <- lav_lisrel_th(
        MLIST = GLIST[mm.in.group],
        th.idx = th.idx[[g]], delta = delta
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }
  }

  TH
}

# PI = slope structure
# see Muthen 1984 eq 12
lav_model_pi <- function(lavmodel = NULL, GLIST = NULL, delta = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  conditional.x <- lavmodel@conditional.x

  # return a list
  PI <- vector("list", length = nblocks)

  # compute TH for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (!conditional.x) {
      PI.g <- numeric(lavmodel@nvar[g])
    } else if (representation == "LISREL") {
      PI.g <- lav_lisrel_pi(MLIST = MLIST, delta = delta)
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    PI[[g]] <- PI.g
  }

  PI
}


# GW = group weight
lav_model_gw <- function(lavmodel = NULL, GLIST = NULL) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation
  group.w.free <- lavmodel@group.w.free

  # return a list
  GW <- vector("list", length = nblocks)

  # compute GW for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (!group.w.free) {
      GW.g <- 0.0 # FIXME
    } else if (representation == "LISREL") {
      GW.g <- as.numeric(MLIST$gw[1, 1])
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    GW[[g]] <- GW.g
  }

  # transform to proportions
  # gw <- unlist(GW)
  # gw <- exp(gw) / sum(exp(gw))
  # for(g in 1:nblocks) {
  #    GW[[g]] <- gw[g]
  # }

  GW
}

# *unconditional* variance/covariance matrix of Y
#  - same as Sigma.hat if all Y are continuous)
#  - if also Gamma, cov.x is used (only if categorical)
lav_model_vy <- function(lavmodel = NULL, GLIST = NULL, diagonal.only = FALSE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  VY <- vector("list", length = nblocks)

  # compute TH for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      VY.g <- lav_lisrel_vy(MLIST = MLIST)
    } else if (representation == "RAM") {
      # does not work for categorical setting yet
      stopifnot(!lavmodel@categorical)
      # does not work if conditional.x = TRUE
      stopifnot(!lavmodel@conditional.x)
      VY.g <- lav_ram_sigmahat(MLIST = MLIST)
    } else {
      lav_msg_stop(gettext(
        "only RAM and LISREL representation has been implemented for now"))
    }

    if (diagonal.only) {
      VY[[g]] <- diag(VY.g)
    } else {
      VY[[g]] <- VY.g
    }
  }

  VY
}

# V(ETA): latent variances variances/covariances
lav_model_veta <- function(lavmodel = NULL, GLIST = NULL,
                           remove.dummy.lv = FALSE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  VETA <- vector("list", length = nblocks)

  # compute VETA for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      VETA.g <- lav_lisrel_veta(MLIST = MLIST)

      if (remove.dummy.lv) {
        # remove all dummy latent variables
        lv.idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (!is.null(lv.idx)) {
          VETA.g <- VETA.g[-lv.idx, -lv.idx, drop = FALSE]
        }
      }
    } else if (representation == "RAM") {
      VETA.g <- lav_ram_veta(MLIST = MLIST)
    } else {
      lav_msg_stop(gettext(
        "only LISREL and RAM representation has been implemented for now"))
    }

    VETA[[g]] <- VETA.g
  }

  VETA
}

# V(ETA|x_i): latent variances variances/covariances, conditional on x_
# - this is always (I-B)^-1 PSI (I-B)^-T, after REMOVING lv dummies
lav_model_vetax <- function(lavmodel = NULL, GLIST = NULL) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  ETA <- vector("list", length = nblocks)

  # compute ETA for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      lv.idx <- c(
        lavmodel@ov.y.dummy.lv.idx[[g]],
        lavmodel@ov.x.dummy.lv.idx[[g]]
      )
      ETA.g <- lav_lisrel_vetax(
        MLIST = MLIST,
        lv.dummy.idx = lv.idx
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    ETA[[g]] <- ETA.g
  }

  ETA
}

# COV: observed+latent variances variances/covariances
lav_model_cov_both <- function(lavmodel = NULL, GLIST = NULL,
                               remove.dummy.lv = FALSE, delta = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  COV <- vector("list", length = nblocks)

  # compute COV for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      COV.g <- lav_lisrel_cov_both(MLIST = MLIST, delta = delta)

      if (remove.dummy.lv) {
        # remove all dummy latent variables
        lv.idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (!is.null(lv.idx)) {
          # offset for ov
          lambda.names <-
            lavmodel@dimNames[[which(names(GLIST) == "lambda")[g]]][[1L]]
          lv.idx <- lv.idx + length(lambda.names)
          COV.g <- COV.g[-lv.idx, -lv.idx, drop = FALSE]
        }
      }
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    COV[[g]] <- COV.g
  }

  COV
}


# E(ETA): expectation (means) of latent variables (return vector)
lav_model_eeta <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL,
                           remove.dummy.lv = FALSE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  EETA <- vector("list", length = nblocks)

  # compute E(ETA) for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      EETA.g <- lav_lisrel_eeta(MLIST,
        mean.x = lavsamplestats@mean.x[[g]],
        sample.mean = lavsamplestats@mean[[g]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]]
      )
      if (remove.dummy.lv) {
        # remove dummy
        lv.dummy.idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (length(lv.dummy.idx) > 0L) {
          EETA.g <- EETA.g[-lv.dummy.idx]
        }
      }
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    EETA[[g]] <- EETA.g
  }

  EETA
}

# E(ETA|x_i): conditional expectation (means) of latent variables
# for a given value of x_i (instead of E(x_i))
lav_model_eetax <- function(lavmodel = NULL, GLIST = NULL,
                            lavsamplestats = NULL, eXo = NULL,
                            nobs = NULL, remove.dummy.lv = FALSE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  EETAx <- vector("list", length = nblocks)

  # compute E(ETA) for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    EXO <- eXo[[g]]
    if (is.null(EXO)) {
      # create empty matrix
      EXO <- matrix(0, nobs[[g]], 0L)
    }

    if (representation == "LISREL") {
      EETAx.g <- lav_lisrel_eetax(MLIST,
        eXo = EXO, N = nobs[[g]],
        sample.mean = lavsamplestats@mean[[g]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]]
      )

      if (remove.dummy.lv) {
        # remove dummy
        lv.dummy.idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[g]],
          lavmodel@ov.x.dummy.lv.idx[[g]]
        )
        if (length(lv.dummy.idx) > 0L) {
          EETAx.g <- EETAx.g[, -lv.dummy.idx, drop = FALSE]
        }
      }
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    EETAx[[g]] <- EETAx.g
  }

  EETAx
}

# return 'regular' LAMBDA
lav_model_lambda <- function(lavmodel = NULL, GLIST = NULL,
                             handle.dummy.lv = TRUE,
                             remove.dummy.lv = FALSE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  LAMBDA <- vector("list", length = nblocks)

  # compute LAMBDA for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      if (handle.dummy.lv) {
        ov.y.dummy.ov.idx <- lavmodel@ov.y.dummy.ov.idx[[g]]
        ov.x.dummy.ov.idx <- lavmodel@ov.x.dummy.ov.idx[[g]]
        ov.y.dummy.lv.idx <- lavmodel@ov.y.dummy.lv.idx[[g]]
        ov.x.dummy.lv.idx <- lavmodel@ov.x.dummy.lv.idx[[g]]
      } else {
        ov.y.dummy.ov.idx <- NULL
        ov.x.dummy.ov.idx <- NULL
        ov.y.dummy.lv.idx <- NULL
        ov.x.dummy.lv.idx <- NULL
      }
      LAMBDA.g <- lav_lisrel_lambda(
        MLIST = MLIST,
        ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
        ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
        ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
        ov.x.dummy.lv.idx = ov.x.dummy.lv.idx,
        remove.dummy.lv = remove.dummy.lv
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    LAMBDA[[g]] <- LAMBDA.g
  }

  LAMBDA
}

# THETA: observed (residual) variances
lav_model_theta <- function(lavmodel = NULL, GLIST = NULL, fix = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  THETA <- vector("list", length = nblocks)

  # compute THETA for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      if (fix) {
        THETA.g <- lav_lisrel_theta(
          MLIST = MLIST,
          ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
          ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
          ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
          ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]]
        )
      } else {
        THETA.g <- lav_lisrel_theta(MLIST = MLIST)
      }
    } else if (representation == "RAM") {
      ov.idx <- as.integer(MLIST$ov.idx[1, ])
      THETA.g <- MLIST$S[ov.idx, ov.idx, drop = FALSE]
    } else {
      lav_msg_stop(gettext(
        "only LISREL and RAM representation has been implemented for now"))
    }

    THETA[[g]] <- THETA.g
  }

  THETA
}

# NU: observed intercepts
lav_model_nu <- function(lavmodel = NULL, GLIST = NULL,
                         lavsamplestats = NULL) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  NU <- vector("list", length = nblocks)

  # compute NU for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      NU.g <- lav_lisrel_nu(
        MLIST = MLIST,
        sample.mean = lavsamplestats@mean[[g]],
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]]
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    NU[[g]] <- as.matrix(NU.g)
  }

  NU
}

# E(Y): expectation (mean) of observed variables
# returns vector 1 x nvar
lav_model_ey <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL,
                         delta = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  EY <- vector("list", length = nblocks)

  # compute E(Y) for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      EY.g <- lav_lisrel_ey(
        MLIST = MLIST,
        mean.x = lavsamplestats@mean.x[[g]],
        sample.mean = lavsamplestats@mean[[g]],
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        delta = delta
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    EY[[g]] <- EY.g
  }

  EY
}

# E(Y|x_i): conditional expectation (mean) of observed variables
# returns matrix N x nvar
lav_model_eyx <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL,
                          eXo = NULL, delta = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  nblocks <- lavmodel@nblocks
  nmat <- lavmodel@nmat
  representation <- lavmodel@representation

  # return a list
  EYx <- vector("list", length = nblocks)

  # compute E(Y) for each group
  for (g in 1:nblocks) {
    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (representation == "LISREL") {
      EYx.g <- lav_lisrel_eyx(
        MLIST = MLIST,
		eXo = eXo[[g]],
        sample.mean = lavsamplestats@mean[[g]],
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        delta = delta
      )
    } else {
      lav_msg_stop(gettext(
        "only representation LISREL has been implemented for now"))
    }

    EYx[[g]] <- EYx.g
  }

  EYx
}


# E(Y | ETA, x_i): conditional expectation (means) of observed variables
# for a given value of x_i AND eta_i
lav_model_yhat <- function(lavmodel = NULL, GLIST = NULL, lavsamplestats = NULL,
                           eXo = NULL, nobs = NULL, ETA = NULL,
                           duplicate = FALSE, delta = TRUE) {
  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  # ngroups, not nblocks!
  ngroups <- lavsamplestats@ngroups

  # return a list
  YHAT <- vector("list", length = ngroups)

  # compute YHAT for each group
  for (g in seq_len(ngroups)) {
    # which mm belong to group g?
    # FIXME: what if more than g blocks???
    mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0L, lavmodel@nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (is.null(eXo[[g]]) && duplicate) {
      Nobs <- nobs[[g]]
    } else {
      Nobs <- 1L
    }

    if (lavmodel@representation == "LISREL") {
      if (lavmodel@conditional.x) {
        YHAT[[g]] <- lav_lisrel_eyetax(
          MLIST = MLIST,
          eXo = eXo[[g]], ETA = ETA[[g]], N = Nobs,
          sample.mean = lavsamplestats@mean[[g]],
          ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
          ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
          ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
          ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
          delta = delta
        )
      } else {
        # unconditional case
        YHAT[[g]] <- lav_lisrel_eyetax3(
          MLIST = MLIST,
          ETA = ETA[[g]],
          sample.mean = lavsamplestats@mean[[g]],
          mean.x = lavsamplestats@mean.x[[g]],
          ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
          ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
          ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
          ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
          delta = delta
        )
        # impute back ov.y values that are NOT indicators
      }
    } else {
      lav_msg_stop(gettextf("representation %s not supported yet.",
                            lavmodel@representation))
    }
  }

  YHAT
}
