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
    lavdata <- lavData(
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
      out <- lav_samplestats_Gamma(
        Y = Y[[g]],
        Mu = NULL,
        Sigma = NULL,
        x.idx = x.idx[[g]],
        cluster.idx = cluster.idx,
        fixed.x = fixed.x,
        conditional.x = conditional.x,
        meanstructure = meanstructure,
        slopestructure = conditional.x,
        gamma.n.minus.one = gamma.n.minus.one,
        unbiased = gamma.unbiased,
        Mplus.WLS = Mplus.WLS
      )
    } else {
      out <- lav_samplestats_Gamma_NT(
        Y = Y[[g]],
        wt = NULL, # for now
        rescale = NT.rescale,
        x.idx = x.idx[[g]],
        fixed.x = fixed.x,
        conditional.x = conditional.x,
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
    if (.hasSlot(lavobject, "h1")) {
      lavh1 <- lavobject@h1
    } else {
      # only for <0.6
      out <- lav_h1_implied_logl(
        lavdata = lavdata,
        lavsamplestats = lavsamplestats,
        lavpartable = lavobject@ParTable,
        lavoptions = lavoptions
      )
      h1.implied <- out$implied
      h1.loglik <- out$logl$loglik
      h1.loglik.group <- out$logl$loglik.group
      lavh1 <- list(
        implied = h1.implied,
        loglik = h1.loglik,
        loglik.group = h1.loglik.group
      )
      lavoptions$gamma.n.minus.one <- FALSE
      lavoptions$gamma.unbiased <- FALSE
    }
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
      OUT[[g]] <- lav_samplestats_Gamma(
        Y = Y,
        Mu = MEAN,
        Sigma = COV,
        x.idx = x.idx,
        cluster.idx = cluster.idx,
        fixed.x = fixed.x,
        conditional.x = conditional.x,
        meanstructure = meanstructure,
        slopestructure = conditional.x,
        gamma.n.minus.one = gamma.n.minus.one,
        unbiased = gamma.unbiased,
        Mplus.WLS = Mplus.WLS
      )
    } else {
      OUT[[g]] <- lav_samplestats_Gamma_NT(
        COV = COV, # joint!
        MEAN = MEAN, # joint!
        x.idx = x.idx,
        fixed.x = fixed.x,
        conditional.x = conditional.x,
        meanstructure = meanstructure,
        slopestructure = conditional.x
      )
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
lav_samplestats_Gamma_NT <- function(Y = NULL, # should include
                                     # eXo if
                                     # conditional.x=TRUE
                                     wt = NULL,
                                     COV = NULL, # joint!
                                     MEAN = NULL, # joint!
                                     rescale = FALSE,
                                     x.idx = integer(0L),
                                     fixed.x = FALSE,
                                     conditional.x = FALSE,
                                     meanstructure = FALSE,
                                     slopestructure = FALSE) {
  # check arguments
  if (length(x.idx) == 0L) {
    conditional.x <- FALSE
    fixed.x <- FALSE
  }

  # compute COV from Y
  if (is.null(COV)) {
    stopifnot(!is.null(Y))

    # coerce to matrix
    Y <- unname(as.matrix(Y))
    N <- nrow(Y)
    if (is.null(wt)) {
      COV <- cov(Y)
      if (rescale) {
        COV <- COV * (N - 1) / N # (normal) ML version
      }
    } else {
      out <- stats::cov.wt(Y, wt = wt, method = "ML")
      COV <- out$cov
    }
  } else {
    if (!missing(rescale)) {
      lav_msg_warn(gettext("rescale= argument has no effect if COV is given"))
    }
    if (!missing(wt)) {
      lav_msg_warn(gettext("wt= argument has no effect if COV is given"))
    }
  }

  # if needed, compute MEAN from Y
  if (conditional.x && length(x.idx) > 0L && is.null(MEAN) &&
    (meanstructure || slopestructure)) {
    stopifnot(!is.null(Y))
    if (is.null(wt)) {
      MEAN <- colMeans(Y, na.rm = TRUE)
    } else {
      MEAN <- out$center
    }
  }

  # rename
  S <- COV
  M <- MEAN

  # unconditional
  if (!conditional.x) {
    # unconditional - stochastic x
    if (!fixed.x) {
      Gamma <- 2 * lav_matrix_duplication_ginv_pre_post(S %x% S)
      if (meanstructure) {
        Gamma <- lav_matrix_bdiag(S, Gamma)
      }

      # unconditional - fixed x
    } else {
      # handle fixed.x = TRUE

      # cov(Y|X) = A - B C^{-1} B'
      # where A = cov(Y), B = cov(Y,X), C = cov(X)
      A <- S[-x.idx, -x.idx, drop = FALSE]
      B <- S[-x.idx, x.idx, drop = FALSE]
      C <- S[x.idx, x.idx, drop = FALSE]
      YbarX <- A - B %*% solve(C, t(B))

      # reinsert YbarX in Y+X (residual) covariance matrix
      YbarX.aug <- matrix(0, nrow = NROW(S), ncol = NCOL(S))
      YbarX.aug[-x.idx, -x.idx] <- YbarX

      # take difference
      R <- S - YbarX.aug

      Gamma.S <- 2 * lav_matrix_duplication_ginv_pre_post(S %x% S)
      Gamma.R <- 2 * lav_matrix_duplication_ginv_pre_post(R %x% R)
      Gamma <- Gamma.S - Gamma.R

      if (meanstructure) {
        Gamma <- lav_matrix_bdiag(YbarX.aug, Gamma)
      }
    }
  } else {
    # conditional.x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    # regress Y on X, and compute covariance of residuals 'R'
    A <- S[-x.idx, -x.idx, drop = FALSE]
    B <- S[-x.idx, x.idx, drop = FALSE]
    C <- S[x.idx, x.idx, drop = FALSE]
    Cov.YbarX <- A - B %*% solve(C) %*% t(B)
    Gamma <- 2 * lav_matrix_duplication_ginv_pre_post(Cov.YbarX %x% Cov.YbarX)

    if (meanstructure || slopestructure) {
      MY <- M[-x.idx]
      MX <- M[x.idx]
      C3 <- rbind(
        c(1, MX),
        cbind(MX, C + tcrossprod(MX))
      )
      # B3 <- cbind(MY, B + tcrossprod(MY,MX))
    }

    if (meanstructure) {
      if (slopestructure) {
        # A11 <- solve(C3) %x% Cov.YbarX
        A11 <- Cov.YbarX %x% solve(C3)
      } else {
        # A11 <- solve(C3)[1, 1, drop=FALSE] %x% Cov.YbarX
        A11 <- Cov.YbarX %x% solve(C3)[1, 1, drop = FALSE]
      }
    } else {
      if (slopestructure) {
        # A11 <- solve(C3)[-1, -1, drop=FALSE] %x% Cov.YbarX
        A11 <- Cov.YbarX %x% solve(C3)[-1, -1, drop = FALSE]
      } else {
        A11 <- matrix(0, 0, 0)
      }
    }

    Gamma <- lav_matrix_bdiag(A11, Gamma)
  }

  Gamma
}

# NOTE:
#  - three types:
#       1) plain         (conditional.x = FALSE, fixed.x = FALSE)
#       2) fixed.x       (conditional.x = FALSE, fixed.x = TRUE)
#       3) conditional.x (conditional.x = TRUE)
#  - if conditional.x = TRUE, we ignore fixed.x (can be TRUE or FALSE)
#
#  - new in 0.6-1: if Mu/Sigma is provided, compute 'model-based' Gamma
#                  (only if conditional.x = FALSE, for now)
#  - new in 0.6-2: if cluster.idx is not NULL, correct for clustering
#  - new in 0.6-13: add unbiased = TRUE (for the 'plain' setting only)

# ADF THEORY
lav_samplestats_Gamma <- function(Y, # Y+X if cond!
                                  Mu = NULL,
                                  Sigma = NULL,
                                  x.idx = integer(0L),
                                  cluster.idx = NULL,
                                  fixed.x = FALSE,
                                  conditional.x = FALSE,
                                  meanstructure = FALSE,
                                  slopestructure = FALSE,
                                  gamma.n.minus.one = FALSE,
                                  unbiased = FALSE,
                                  Mplus.WLS = FALSE) {
  # coerce to matrix
  Y <- unname(as.matrix(Y))
  N <- nrow(Y)
  p <- ncol(Y)

  # unbiased?
  if (unbiased) {
    if (conditional.x || fixed.x || !is.null(Sigma) ||
      !is.null(cluster.idx)) {
      lav_msg_stop(
        gettext("unbiased Gamma only available for the simple
        (not conditional.x or fixed.x or model-based or clustered) setting."))
    } else {
      COV <- COV.unbiased <- cov(Y)
      COV <- COV * (N - 1) / N
      cov.vech <- lav_matrix_vech(COV)
    }
  }

  # model-based?
  if (!is.null(Sigma)) {
    stopifnot(!conditional.x)
    model.based <- TRUE
    if (meanstructure) {
      stopifnot(!is.null(Mu))
      sigma <- c(as.numeric(Mu), lav_matrix_vech(Sigma))
    } else {
      Mu <- colMeans(Y, na.rm = TRUE) # for centering!
      sigma <- lav_matrix_vech(Sigma)
    }
  } else {
    model.based <- FALSE
  }

  # denominator
  if (gamma.n.minus.one) {
    N <- N - 1
  }

  # check arguments
  if (length(x.idx) == 0L) {
    conditional.x <- FALSE
    fixed.x <- FALSE
  }
  if (Mplus.WLS) {
    stopifnot(!conditional.x, !fixed.x)
  }

  if (!conditional.x && !fixed.x) {
    # center, so we can use crossprod instead of cov
    if (model.based) {
      Yc <- t(t(Y) - as.numeric(Mu))
    } else {
      Yc <- t(t(Y) - colMeans(Y, na.rm = TRUE))
    }

    # create Z where the rows_i contain the following elements:
    #  - Y_i (if meanstructure is TRUE)
    #  - vech(Yc_i' %*% Yc_i) where Yc_i are the residuals
    idx1 <- lav_matrix_vech_col_idx(p)
    idx2 <- lav_matrix_vech_row_idx(p)
    if (meanstructure) {
      Z <- cbind(Y, Yc[, idx1, drop = FALSE] *
        Yc[, idx2, drop = FALSE])
    } else {
      Z <- (Yc[, idx1, drop = FALSE] *
        Yc[, idx2, drop = FALSE])
    }

    if (model.based) {
      if (meanstructure) {
        stopifnot(!is.null(Mu))
        sigma <- c(as.numeric(Mu), lav_matrix_vech(Sigma))
      } else {
        sigma <- lav_matrix_vech(Sigma)
      }
      Zc <- t(t(Z) - sigma)
    } else {
      Zc <- t(t(Z) - colMeans(Z, na.rm = TRUE))
    }

    # clustered?
    if (length(cluster.idx) > 0L) {
      Zc <- rowsum(Zc, cluster.idx)
    }

    if (anyNA(Zc)) {
      Gamma <- lav_matrix_crossprod(Zc) / N
    } else {
      Gamma <- base::crossprod(Zc) / N
    }
  } else if (!conditional.x && fixed.x) {
    if (model.based) {
      Yc <- t(t(Y) - as.numeric(Mu))
      Y.bar <- colMeans(Y, na.rm = TRUE)
      res.cov <- (Sigma[-x.idx, -x.idx, drop = FALSE] -
        Sigma[-x.idx, x.idx, drop = FALSE] %*%
        solve(Sigma[x.idx, x.idx, drop = FALSE]) %*%
        Sigma[x.idx, -x.idx, drop = FALSE])
      res.slopes <- (solve(Sigma[x.idx, x.idx, drop = FALSE]) %*%
        Sigma[x.idx, -x.idx, drop = FALSE])
      res.int <- (Y.bar[-x.idx] -
        as.numeric(colMeans(Y[, x.idx, drop = FALSE],
          na.rm = TRUE
        ) %*% res.slopes))
      x.bar <- Y.bar[x.idx]
      yhat.bar <- as.numeric(res.int + as.numeric(x.bar) %*% res.slopes)
      YHAT.bar <- numeric(p)
      YHAT.bar[-x.idx] <- yhat.bar
      YHAT.bar[x.idx] <- x.bar
      YHAT.cov <- Sigma
      YHAT.cov[-x.idx, -x.idx] <- Sigma[-x.idx, -x.idx] - res.cov


      yhat <- cbind(1, Y[, x.idx]) %*% rbind(res.int, res.slopes)
      YHAT <- Y
      YHAT[, -x.idx] <- yhat
      # YHAT <- cbind(yhat, Y[,x.idx])
      YHATc <- t(t(YHAT) - YHAT.bar)
      idx1 <- lav_matrix_vech_col_idx(p)
      idx2 <- lav_matrix_vech_row_idx(p)
      if (meanstructure) {
        Z <- (cbind(Y, Yc[, idx1, drop = FALSE] *
          Yc[, idx2, drop = FALSE]) -
          cbind(YHAT, YHATc[, idx1, drop = FALSE] *
            YHATc[, idx2, drop = FALSE]))
        sigma1 <- c(Mu, lav_matrix_vech(Sigma))
        sigma2 <- c(YHAT.bar, lav_matrix_vech(YHAT.cov))
      } else {
        Z <- (Yc[, idx1, drop = FALSE] *
          Yc[, idx2, drop = FALSE] -
          YHATc[, idx1, drop = FALSE] *
            YHATc[, idx2, drop = FALSE])
        sigma1 <- lav_matrix_vech(Sigma)
        sigma2 <- lav_matrix_vech(YHAT.cov)
      }
      Zc <- t(t(Z) - (sigma1 - sigma2))
    } else {
      QR <- qr(cbind(1, Y[, x.idx, drop = FALSE]))
      yhat <- qr.fitted(QR, Y[, -x.idx, drop = FALSE])
      # YHAT <- cbind(yhat, Y[,x.idx])
      YHAT <- Y
      YHAT[, -x.idx] <- yhat

      Yc <- t(t(Y) - colMeans(Y, na.rm = TRUE))
      YHATc <- t(t(YHAT) - colMeans(YHAT, na.rm = TRUE))
      idx1 <- lav_matrix_vech_col_idx(p)
      idx2 <- lav_matrix_vech_row_idx(p)
      if (meanstructure) {
        Z <- (cbind(Y, Yc[, idx1, drop = FALSE] *
          Yc[, idx2, drop = FALSE]) -
          cbind(YHAT, YHATc[, idx1, drop = FALSE] *
            YHATc[, idx2, drop = FALSE]))
      } else {
        Z <- (Yc[, idx1, drop = FALSE] *
          Yc[, idx2, drop = FALSE] -
          YHATc[, idx1, drop = FALSE] *
            YHATc[, idx2, drop = FALSE])
      }
      Zc <- t(t(Z) - colMeans(Z, na.rm = TRUE))
    }

    # clustered?
    if (length(cluster.idx) > 0L) {
      Zc <- rowsum(Zc, cluster.idx)
    }

    if (anyNA(Zc)) {
      Gamma <- lav_matrix_crossprod(Zc) / N
    } else {
      Gamma <- base::crossprod(Zc) / N
    }
  } else {
    # conditional.x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    # regress Y on X, and compute residuals
    X <- cbind(1, Y[, x.idx, drop = FALSE])
    QR <- qr(X)
    RES <- qr.resid(QR, Y[, -x.idx, drop = FALSE])
    p <- ncol(RES)

    idx1 <- lav_matrix_vech_col_idx(p)
    idx2 <- lav_matrix_vech_row_idx(p)

    if (meanstructure || slopestructure) {
      XtX.inv <- unname(solve(crossprod(X)))
      Xi <- (X %*% XtX.inv) * N ## FIXME, shorter way?
      ncX <- NCOL(X)
      ncY <- NCOL(RES)
    }

    if (meanstructure) {
      if (slopestructure) {
        # Xi.idx <- rep(seq_len(ncX), each  = ncY)
        # Res.idx <- rep(seq_len(ncY), times = ncX)
        Xi.idx <- rep(seq_len(ncX), times = ncY)
        Res.idx <- rep(seq_len(ncY), each = ncX)
        Z <- cbind(
          Xi[, Xi.idx, drop = FALSE] *
            RES[, Res.idx, drop = FALSE],
          RES[, idx1, drop = FALSE] *
            RES[, idx2, drop = FALSE]
        )
      } else {
        Xi.idx <- rep(1L, each = ncY)
        Z <- cbind(
          Xi[, Xi.idx, drop = FALSE] *
            RES,
          RES[, idx1, drop = FALSE] *
            RES[, idx2, drop = FALSE]
        )
      }
    } else {
      if (slopestructure) {
        # Xi.idx <- rep(seq_len(ncX), each  = ncY)
        # Xi.idx <- Xi.idx[ -seq_len(ncY) ]
        Xi.idx <- rep(seq(2, ncX), times = ncY)
        # Res.idx <- rep(seq_len(ncY), times = (ncX - 1L))
        Res.idx <- rep(seq_len(ncY), each = (ncX - 1L))
        Z <- cbind(
          Xi[, Xi.idx, drop = FALSE] *
            RES[, Res.idx, drop = FALSE],
          RES[, idx1, drop = FALSE] *
            RES[, idx2, drop = FALSE]
        )
      } else {
        Z <- RES[, idx1, drop = FALSE] * RES[, idx2, drop = FALSE]
      }
    }

    if (model.based) {
      Zc <- t(t(Z) - sigma)
    } else {
      Zc <- t(t(Z) - colMeans(Z, na.rm = TRUE))
    }

    # clustered?
    if (length(cluster.idx) > 0L) {
      Zc <- rowsum(Zc, cluster.idx)
    }

    if (anyNA(Zc)) {
      Gamma <- lav_matrix_crossprod(Zc) / N
    } else {
      Gamma <- base::crossprod(Zc) / N
    }
  }


  # only to mimic Mplus when estimator = "WLS"
  if (Mplus.WLS && !fixed.x && !conditional.x) {
    # adjust G_22 (the varcov part)
    S <- cov(Y, use = "pairwise")
    w <- lav_matrix_vech(S)
    w.biased <- (N - 1) / N * w
    diff <- outer(w, w) - outer(w.biased, w.biased)
    if (meanstructure) {
      Gamma[-seq_len(p), -seq_len(p)] <-
        Gamma[-seq_len(p), -seq_len(p), drop = FALSE] - diff
    } else {
      Gamma <- Gamma - diff
    }

    if (meanstructure) {
      # adjust G_12/G_21 (third-order)
      # strange rescaling?
      N1 <- (N - 1) / N
      Gamma[seq_len(p), -seq_len(p)] <- Gamma[seq_len(p), -seq_len(p)] * N1
      Gamma[-seq_len(p), seq_len(p)] <- Gamma[-seq_len(p), seq_len(p)] * N1
    }
  }

  # clustered?
  if (length(cluster.idx) > 0L) {
    nC <- nrow(Zc)
    Gamma <- Gamma * nC / (nC - 1)
  }

  # unbiased?
  if (unbiased) {
    # normal-theory Gamma (cov only)
    GammaNT.cov <- 2 * lav_matrix_duplication_ginv_pre_post(COV %x% COV)

    if (meanstructure) {
      Gamma.cov <- Gamma[-(1:p), -(1:p), drop = FALSE]
      Gamma.mean.cov <- Gamma[1:p, -(1:p), drop = FALSE]
    } else {
      Gamma.cov <- Gamma
    }

    # Browne's unbiased DF estimator (COV part)
    Gamma.u <- (N * (N - 1) / (N - 2) / (N - 3) * Gamma.cov -
      N / (N - 2) / (N - 3) * (GammaNT.cov -
        2 / (N - 1) * tcrossprod(cov.vech)))
    if (meanstructure) {
      Gamma <- lav_matrix_bdiag(COV, Gamma.u)

      # 3-rd order:
      Gamma[1:p, (p + 1):ncol(Gamma)] <- Gamma.mean.cov * N / (N - 2)
      Gamma[(p + 1):ncol(Gamma), 1:p] <- t(Gamma.mean.cov * N / (N - 2))
    } else {
      Gamma <- Gamma.u
    }
  } # unbiased

  Gamma
}
