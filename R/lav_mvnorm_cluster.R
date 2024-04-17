# loglikelihood clustered/twolevel data

# YR: first version around Feb 2017


# take model-implied mean+variance matrices, and reorder/augment them
# to facilitate computing of (log)likelihood in the two-level case

# when conditional.x = FALSE:
# - sigma.w and sigma.b: same dimensions, level-1 variables only
# - sigma.zz: level-2 variables only
# - sigma.yz: cov(level-1, level-2)
# - mu.y: level-1 variables only (mu.w + mu.b)
# - mu.w: y within  part
# - mu.b: y between part
# - mu.z: level-2 variables only
lav_mvnorm_cluster_implied22l <- function(Lp = NULL,
                                          implied = NULL,
                                          Mu.W = NULL,
                                          Mu.B = NULL,
                                          Sigma.W = NULL,
                                          Sigma.B = NULL) {
  if (!is.null(implied)) {
    # FIXME: only for single-group analysis!
    Sigma.W <- implied$cov[[1]]
    Mu.W <- implied$mean[[1]]

    Sigma.B <- implied$cov[[2]]
    Mu.B <- implied$mean[[2]]
  }

  # within/between.idx
  between.idx <- Lp$between.idx[[2]]
  within.idx <- Lp$within.idx[[2]]
  both.idx <- Lp$both.idx[[2]]

  # ov.idx per level
  ov.idx <- Lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))

  # Sigma.W.tilde
  Sigma.W.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.W.tilde[ov.idx[[1]], ov.idx[[1]]] <- Sigma.W

  # Sigma.B.tilde
  Sigma.B.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.B.tilde[ov.idx[[2]], ov.idx[[2]]] <- Sigma.B

  # Mu.W.tilde
  Mu.W.tilde <- numeric(p.tilde)
  Mu.W.tilde[ov.idx[[1]]] <- Mu.W

  # Mu.B.tilde
  Mu.B.tilde <- numeric(p.tilde)
  Mu.B.tilde[ov.idx[[2]]] <- Mu.B


  # add Mu.W[within.idx] to Mu.B
  Mu.WB.tilde <- numeric(p.tilde)
  Mu.WB.tilde[within.idx] <- Mu.W.tilde[within.idx]
  Mu.WB.tilde[both.idx] <- (Mu.B.tilde[both.idx] +
    Mu.W.tilde[both.idx])

  # map to matrices needed for loglik
  if (length(within.idx) > 0L) {
    Mu.B.tilde[within.idx] <- 0
  }
  if (length(between.idx) > 0L) {
    mu.z <- Mu.B.tilde[between.idx]
    mu.y <- Mu.WB.tilde[-between.idx]
    mu.w <- Mu.W.tilde[-between.idx]
    mu.b <- Mu.B.tilde[-between.idx]
    sigma.zz <- Sigma.B.tilde[between.idx, between.idx, drop = FALSE]
    sigma.yz <- Sigma.B.tilde[-between.idx, between.idx, drop = FALSE]
    sigma.b <- Sigma.B.tilde[-between.idx, -between.idx, drop = FALSE]
    sigma.w <- Sigma.W.tilde[-between.idx, -between.idx, drop = FALSE]
  } else {
    mu.z <- numeric(0L)
    mu.y <- Mu.WB.tilde
    mu.w <- Mu.W.tilde
    mu.b <- Mu.B.tilde
    sigma.zz <- matrix(0, 0L, 0L)
    sigma.yz <- matrix(0, nrow(Sigma.B.tilde), 0L)
    sigma.b <- Sigma.B.tilde
    sigma.w <- Sigma.W.tilde
  }

  list(
    sigma.w = sigma.w, sigma.b = sigma.b, sigma.zz = sigma.zz,
    sigma.yz = sigma.yz, mu.z = mu.z, mu.y = mu.y, mu.w = mu.w,
    mu.b = mu.b
  )
}

lav_mvnorm_cluster_2l2implied <- function(Lp,
                                          sigma.w = NULL,
                                          sigma.b = NULL,
                                          sigma.zz = NULL,
                                          sigma.yz = NULL,
                                          mu.z = NULL,
                                          mu.y = NULL,
                                          mu.w = NULL,
                                          mu.b = NULL) {
  # between.idx
  between.idx <- Lp$between.idx[[2]]
  within.idx <- Lp$within.idx[[2]]

  # ov.idx per level
  ov.idx <- Lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))

  # if we have mu.y, convert to mu.w and mu.b
  if (!is.null(mu.y)) {
    mu.b <- mu.y
    mu.w.tilde <- numeric(p.tilde)
    mu.w.tilde[ov.idx[[1]]] <- mu.y

    # NO NEED TO SET THIS TO ZERO!
    # otherwise, we get non-symmetric Hessian!! 0.6-5

    # if(length(within.idx) > 0L) {
    #    mu.w.tilde[  -within.idx ] <- 0
    # } else {
    #    mu.w.tilde[] <- 0
    # }
    mu.w <- mu.w.tilde[ov.idx[[1]]]
  }

  Mu.W.tilde <- numeric(p.tilde)
  ########## DEBUG ##############
  # if(length(within.idx) > 0) {
  Mu.W.tilde[ov.idx[[1]]] <- mu.w
  # }
  ###############################
  Mu.W <- Mu.W.tilde[ov.idx[[1]]]

  # Mu.B
  Mu.B.tilde <- numeric(p.tilde)
  Mu.B.tilde[ov.idx[[1]]] <- mu.b
  Mu.B.tilde[between.idx] <- mu.z
  if (length(within.idx) > 0) {
    Mu.B.tilde[within.idx] <- 0
  }
  Mu.B <- Mu.B.tilde[ov.idx[[2]]]

  # Sigma.W
  Sigma.W <- sigma.w

  # Sigma.B
  Sigma.B.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.B.tilde[ov.idx[[1]], ov.idx[[1]]] <- sigma.b
  Sigma.B.tilde[ov.idx[[1]], between.idx] <- sigma.yz
  Sigma.B.tilde[between.idx, ov.idx[[1]]] <- t(sigma.yz)
  Sigma.B.tilde[between.idx, between.idx] <- sigma.zz
  Sigma.B <- Sigma.B.tilde[ov.idx[[2]], ov.idx[[2]], drop = FALSE]

  list(Mu.W = Mu.W, Mu.B = Mu.B, Sigma.W = Sigma.W, Sigma.B = Sigma.B)
}


# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
# (not yet reordered)
lav_mvnorm_cluster_loglik_samplestats_2l <- function(YLp = NULL,
                                                     Lp = NULL,
                                                     Mu.W = NULL,
                                                     Sigma.W = NULL,
                                                     Mu.B = NULL,
                                                     Sigma.B = NULL,
                                                     Sinv.method = "eigen",
                                                     log2pi = FALSE,
                                                     minus.two = TRUE) {
  # map implied to 2l matrices
  out <- lav_mvnorm_cluster_implied22l(
    Lp = Lp, Mu.W = Mu.W, Mu.B = Mu.B,
    Sigma.W = Sigma.W, Sigma.B = Sigma.B
  )
  mu.y <- out$mu.y
  mu.z <- out$mu.z
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz

  # Lp
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  between.idx <- Lp$between.idx[[2]]
  cluster.sizes <- Lp$cluster.sizes[[2]]
  ncluster.sizes <- Lp$ncluster.sizes[[2]]
  cluster.size.ns <- Lp$cluster.size.ns[[2]]

  # Y1 samplestats
  if (length(between.idx) > 0L) {
    S.PW <- YLp[[2]]$Sigma.W[-between.idx, -between.idx, drop = FALSE]
  } else {
    S.PW <- YLp[[2]]$Sigma.W
  }

  # Y2 samplestats
  cov.d <- YLp[[2]]$cov.d
  mean.d <- YLp[[2]]$mean.d

  # common parts:
  sigma.w.inv <- lav_matrix_symmetric_inverse(
    S = sigma.w,
    logdet = TRUE, Sinv.method = Sinv.method
  )
  sigma.w.logdet <- attr(sigma.w.inv, "logdet")
  attr(sigma.w.inv, "logdet") <- NULL

  if (length(between.idx) > 0L) {
    sigma.zz.inv <- lav_matrix_symmetric_inverse(
      S = sigma.zz,
      logdet = TRUE, Sinv.method = Sinv.method
    )
    sigma.zz.logdet <- attr(sigma.zz.inv, "logdet")
    attr(sigma.zz.inv, "logdet") <- NULL
    sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
    sigma.zi.zy <- t(sigma.yz.zi)
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy
  } else {
    sigma.zz.logdet <- 0
    sigma.b.z <- sigma.b
  }

  # min 2* logliklihood
  L <- numeric(ncluster.sizes) # logdet
  B <- numeric(ncluster.sizes) # between qf
  for (clz in seq_len(ncluster.sizes)) {
    # cluster size
    nj <- cluster.sizes[clz]

    # data between
    Y2Yc <- (cov.d[[clz]] + tcrossprod(mean.d[[clz]] - c(mu.z, mu.y)))

    # FIXME: avoid reorder/b.idx, so we can use between.idx
    if (length(between.idx) > 0L) {
      b.idx <- seq_len(length(Lp$between.idx[[2]]))
      Y2Yc.zz <- Y2Yc[b.idx, b.idx, drop = FALSE]
      Y2Yc.yz <- Y2Yc[-b.idx, b.idx, drop = FALSE]
      Y2Yc.yy <- Y2Yc[-b.idx, -b.idx, drop = FALSE]
    } else {
      Y2Yc.yy <- Y2Yc
    }

    # construct sigma.j
    sigma.j <- (nj * sigma.b.z) + sigma.w
    sigma.j.inv <- lav_matrix_symmetric_inverse(
      S = sigma.j,
      logdet = TRUE, Sinv.method = Sinv.method
    )
    sigma.j.logdet <- attr(sigma.j.inv, "logdet")
    attr(sigma.j.inv, "logdet") <- NULL

    # check: what if sigma.j is non-pd? should not happen
    if (is.na(sigma.j.logdet)) {
      # stop, and return NA right away
      # return(as.numeric(NA))
      # FORCE?
      # sigma.j <- lav_matrix_symmetric_force_pd(sigma.j)
      # sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j,
      #           logdet = TRUE, Sinv.method = Sinv.method)
      # sigma.j.logdet <- attr(sigma.j.inv, "logdet")
      # attr(sigma.j.inv, "logdet") <- NULL
    }

    # logdet -- between only
    L[clz] <- (sigma.zz.logdet + sigma.j.logdet)

    if (length(between.idx) > 0L) {
      # part 1 -- zz
      sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
      Vinv.11 <- sigma.zz.inv + nj * (sigma.zi.zy %*% sigma.ji.yz.zi)
      q.zz <- sum(Vinv.11 * Y2Yc.zz)

      # part 2 -- yz
      q.yz <- -nj * sum(sigma.ji.yz.zi * Y2Yc.yz)
    } else {
      q.zz <- q.yz <- 0
    }

    # part 5 -- yyc
    q.yyc <- -nj * sum(sigma.j.inv * Y2Yc.yy)

    # qf -- between only
    B[clz] <- q.zz + 2 * q.yz - q.yyc
  }
  # q.yya + q.yyb
  # the reason why we multiply the trace by 'N - nclusters' is
  # S.PW has been divided by 'N - nclusters'
  q.W <- sum(cluster.size - 1) * sum(sigma.w.inv * S.PW)
  # logdet within part
  L.W <- sum(cluster.size - 1) * sigma.w.logdet

  # -2*times logl (without the constant)
  loglik <- sum(L * cluster.size.ns) + sum(B * cluster.size.ns) + q.W + L.W

  # functions below compute -2 * logl
  if (!minus.two) {
    loglik <- loglik / (-2)
  }

  # constant
  # Note: total 'N' = (nobs * #within vars) + (nclusters * #between vars)
  if (log2pi) {
    LOG.2PI <- log(2 * pi)
    nWithin <- length(c(Lp$both.idx[[2]], Lp$within.idx[[2]]))
    nBetween <- length(Lp$between.idx[[2]])
    P <- Lp$nclusters[[1]] * nWithin + Lp$nclusters[[2]] * nBetween
    constant <- -(P * LOG.2PI) / 2
    loglik <- loglik + constant
  }

  # loglik.x (only if loglik is requested)
  if (length(unlist(Lp$ov.x.idx)) > 0L && log2pi && !minus.two) {
    loglik <- loglik - YLp[[2]]$loglik.x
  }

  loglik
}


# first derivative -2*logl wrt Mu.W, Mu.B, Sigma.W, Sigma.B
lav_mvnorm_cluster_dlogl_2l_samplestats <- function(YLp = NULL,
                                                    Lp = NULL,
                                                    Mu.W = NULL,
                                                    Sigma.W = NULL,
                                                    Mu.B = NULL,
                                                    Sigma.B = NULL,
                                                    return.list = FALSE,
                                                    Sinv.method = "eigen") {
  # map implied to 2l matrices
  out <- lav_mvnorm_cluster_implied22l(
    Lp = Lp, Mu.W = Mu.W, Mu.B = Mu.B,
    Sigma.W = Sigma.W, Sigma.B = Sigma.B
  )
  mu.y <- out$mu.y
  mu.z <- out$mu.z
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz

  # Lp
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.sizes <- Lp$cluster.sizes[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]
  between.idx <- Lp$between.idx[[2]]
  ncluster.sizes <- Lp$ncluster.sizes[[2]]
  cluster.size.ns <- Lp$cluster.size.ns[[2]]

  # Y1
  if (length(between.idx) > 0L) {
    S.PW <- YLp[[2]]$Sigma.W[-between.idx, -between.idx, drop = FALSE]
  } else {
    S.PW <- YLp[[2]]$Sigma.W
  }

  # Y2
  cov.d <- YLp[[2]]$cov.d
  mean.d <- YLp[[2]]$mean.d

  # common parts:
  sigma.w.inv <- lav_matrix_symmetric_inverse(
    S = sigma.w,
    logdet = FALSE, Sinv.method = Sinv.method
  )

  # both level-1 and level-2
  G.muy <- matrix(0, ncluster.sizes, length(mu.y))
  G.Sigma.w <- matrix(0, ncluster.sizes, length(lav_matrix_vech(sigma.w)))
  G.Sigma.b <- matrix(0, ncluster.sizes, length(lav_matrix_vech(sigma.b)))

  if (length(between.idx) > 0L) {
    G.muz <- matrix(0, ncluster.sizes, length(mu.z))
    G.Sigma.zz <- matrix(
      0, ncluster.sizes,
      length(lav_matrix_vech(sigma.zz))
    )
    G.Sigma.yz <- matrix(0, ncluster.sizes, length(lav_matrix_vec(sigma.yz)))

    sigma.zz.inv <- lav_matrix_symmetric_inverse(
      S = sigma.zz,
      logdet = FALSE, Sinv.method = Sinv.method
    )
    sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
    sigma.zi.zy <- t(sigma.yz.zi)
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy


    for (clz in seq_len(ncluster.sizes)) {
      # cluster size
      nj <- cluster.sizes[clz]

      # level-2 vectors
      b.idx <- seq_len(length(Lp$between.idx[[2]]))
      zyc <- mean.d[[clz]] - c(mu.z, mu.y)
      yc <- zyc[-b.idx]
      zc <- zyc[b.idx]

      # level-2 crossproducts
      Y2Yc <- (cov.d[[clz]] + tcrossprod(mean.d[[clz]] - c(mu.z, mu.y)))
      b.idx <- seq_len(length(Lp$between.idx[[2]]))
      Y2Yc.zz <- Y2Yc[b.idx, b.idx, drop = FALSE]
      Y2Yc.yz <- Y2Yc[-b.idx, b.idx, drop = FALSE]
      Y2Yc.yy <- Y2Yc[-b.idx, -b.idx, drop = FALSE]

      # construct sigma.j
      sigma.j <- (nj * sigma.b.z) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(
        S = sigma.j,
        logdet = FALSE, Sinv.method = Sinv.method
      )
      sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
      sigma.zi.zy.ji <- t(sigma.ji.yz.zi)

      # common parts
      jYZj <- nj * (sigma.j.inv %*%
        (sigma.yz.zi %*% Y2Yc.zz %*% t(sigma.yz.zi)
          - Y2Yc.yz %*% t(sigma.yz.zi)
          - t(Y2Yc.yz %*% t(sigma.yz.zi)) + Y2Yc.yy)
        %*% sigma.j.inv)

      Z1 <- Y2Yc.zz %*% t(sigma.ji.yz.zi) %*% sigma.yz
      YZ1 <- t(Y2Yc.yz) %*% sigma.j.inv %*% sigma.yz


      # Mu.Z
      G.muz[clz, ] <- -2 * as.numeric(
        (sigma.zz.inv + nj * (sigma.zi.zy.ji %*% sigma.yz.zi)) %*% zc
          - nj * sigma.zi.zy.ji %*% yc
      )

      # MU.Y
      G.muy[clz, ] <- 2 * nj * as.numeric(zc %*% sigma.zi.zy.ji -
        yc %*% sigma.j.inv)

      # SIGMA.W (between part)
      g.sigma.w <- sigma.j.inv - jYZj
      tmp <- g.sigma.w * 2
      diag(tmp) <- diag(g.sigma.w)
      G.Sigma.w[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g.sigma.b <- nj * (sigma.j.inv - jYZj)
      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.Sigma.b[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZZ
      g.sigma.zz <- (sigma.zz.inv + nj * sigma.zz.inv %*% (
        t(sigma.yz) %*% (sigma.j.inv - jYZj) %*% sigma.yz
          - (1 / nj * Y2Yc.zz + t(Z1) + Z1 - t(YZ1) - YZ1)) %*%
        sigma.zz.inv)

      tmp <- g.sigma.zz * 2
      diag(tmp) <- diag(g.sigma.zz)
      G.Sigma.zz[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZY
      g.sigma.yz <- 2 * nj * (
        (sigma.j.inv %*%
          (sigma.yz.zi %*% Y2Yc.zz - sigma.yz - Y2Yc.yz)
          + jYZj %*% sigma.yz) %*% sigma.zz.inv)

      G.Sigma.yz[clz, ] <- lav_matrix_vec(g.sigma.yz)
    }

    # level-1
    d.mu.y <- colSums(G.muy * cluster.size.ns)
    d.sigma.w1 <- lav_matrix_vech_reverse(colSums(G.Sigma.w *
      cluster.size.ns))
    d.sigma.b <- lav_matrix_vech_reverse(colSums(G.Sigma.b *
      cluster.size.ns))

    # level-2
    d.mu.z <- colSums(G.muz * cluster.size.ns)
    d.sigma.zz <- lav_matrix_vech_reverse(colSums(G.Sigma.zz *
      cluster.size.ns))
    d.sigma.yz <- matrix(
      colSums(G.Sigma.yz * cluster.size.ns),
      nrow(sigma.yz), ncol(sigma.yz)
    )
  } # between.idx

  else { # no level-2 variables

    for (clz in seq_len(ncluster.sizes)) {
      # cluster size
      nj <- cluster.sizes[clz]

      # level-2 vectors
      yc <- mean.d[[clz]] - mu.y

      # level-2 crossproducts
      Y2Yc.yy <- (cov.d[[clz]] + tcrossprod(mean.d[[clz]] - mu.y))

      # construct sigma.j
      sigma.j <- (nj * sigma.b) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(
        S = sigma.j,
        logdet = FALSE, Sinv.method = Sinv.method
      )
      # common part
      jYYj <- nj * sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv

      # MU.Y
      G.muy[clz, ] <- -2 * nj * as.numeric(yc %*% sigma.j.inv)

      # SIGMA.W (between part)
      g.sigma.w <- sigma.j.inv - jYYj
      tmp <- g.sigma.w * 2
      diag(tmp) <- diag(g.sigma.w)
      G.Sigma.w[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g.sigma.b <- nj * (sigma.j.inv - jYYj)
      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.Sigma.b[clz, ] <- lav_matrix_vech(tmp)
    }

    # level-1
    d.mu.y <- colSums(G.muy * cluster.size.ns)
    d.sigma.w1 <- lav_matrix_vech_reverse(colSums(G.Sigma.w *
      cluster.size.ns))
    d.sigma.b <- lav_matrix_vech_reverse(colSums(G.Sigma.b *
      cluster.size.ns))
    # level-2
    d.mu.z <- numeric(0L)
    d.sigma.zz <- matrix(0, 0L, 0L)
    d.sigma.yz <- matrix(0, 0L, 0L)
  }

  # Sigma.W (bis)
  d.sigma.w2 <- (Lp$nclusters[[1]] - nclusters) * (sigma.w.inv
  - sigma.w.inv %*% S.PW %*% sigma.w.inv)
  tmp <- d.sigma.w2 * 2
  diag(tmp) <- diag(d.sigma.w2)
  d.sigma.w2 <- tmp

  d.sigma.w <- d.sigma.w1 + d.sigma.w2

  # rearrange
  dout <- lav_mvnorm_cluster_2l2implied(
    Lp = Lp,
    sigma.w = d.sigma.w, sigma.b = d.sigma.b,
    sigma.yz = d.sigma.yz, sigma.zz = d.sigma.zz,
    mu.y = d.mu.y, mu.z = d.mu.z
  )

  if (return.list) {
    out <- dout
  } else {
    out <- c(
      dout$Mu.W, lav_matrix_vech(dout$Sigma.W),
      dout$Mu.B, lav_matrix_vech(dout$Sigma.B)
    )
  }

  out
}

# cluster-wise scores -2*logl wrt Mu.W, Mu.B, Sigma.W, Sigma.B
lav_mvnorm_cluster_scores_2l <- function(Y1 = NULL,
                                         YLp = NULL,
                                         Lp = NULL,
                                         Mu.W = NULL,
                                         Sigma.W = NULL,
                                         Mu.B = NULL,
                                         Sigma.B = NULL,
                                         Sinv.method = "eigen") {
  # map implied to 2l matrices
  out <- lav_mvnorm_cluster_implied22l(
    Lp = Lp, Mu.W = Mu.W, Mu.B = Mu.B,
    Sigma.W = Sigma.W, Sigma.B = Sigma.B
  )
  mu.y <- out$mu.y
  mu.z <- out$mu.z
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz

  # Lp
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]
  between.idx <- Lp$between.idx[[2]]

  # Y1
  if (length(between.idx) > 0L) {
    Y1w <- Y1[, -Lp$between.idx[[2]], drop = FALSE]
  } else {
    Y1w <- Y1
  }
  Y1w.cm <- t(t(Y1w) - mu.y)

  # Y2
  Y2 <- YLp[[2]]$Y2
  # NOTE: ORDER mu.b must match Y2
  mu.b <- numeric(ncol(Y2))
  if (length(between.idx) > 0L) {
    mu.b[-Lp$between.idx[[2]]] <- mu.y
    mu.b[Lp$between.idx[[2]]] <- mu.z
  } else {
    mu.b <- mu.y
  }
  Y2.cm <- t(t(Y2) - mu.b)

  # common parts:
  sigma.w.inv <- lav_matrix_symmetric_inverse(
    S = sigma.w,
    logdet = FALSE, Sinv.method = Sinv.method
  )

  # both level-1 and level-2
  G.muy <- matrix(0, nclusters, length(mu.y))
  G.Sigma.w <- matrix(0, nclusters, length(lav_matrix_vech(sigma.w)))
  G.Sigma.b <- matrix(0, nclusters, length(lav_matrix_vech(sigma.b)))
  G.muz <- matrix(0, nclusters, length(mu.z))
  G.Sigma.zz <- matrix(0, nclusters, length(lav_matrix_vech(sigma.zz)))
  G.Sigma.yz <- matrix(0, nclusters, length(lav_matrix_vec(sigma.yz)))

  if (length(between.idx) > 0L) {
    sigma.zz.inv <- lav_matrix_symmetric_inverse(
      S = sigma.zz,
      logdet = FALSE, Sinv.method = Sinv.method
    )
    sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
    sigma.zi.zy <- t(sigma.yz.zi)
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy


    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster.size[cl]

      # data within for the cluster (centered by mu.y)
      Y1m <- Y1w.cm[cluster.idx == cl, , drop = FALSE]
      yc <- Y2.cm[cl, -Lp$between.idx[[2]]]
      zc <- Y2.cm[cl, Lp$between.idx[[2]]]

      # data between
      Y2Yc <- tcrossprod(Y2.cm[cl, ])
      Y2Yc.zz <- Y2Yc[Lp$between.idx[[2]],
        Lp$between.idx[[2]],
        drop = FALSE
      ]
      Y2Yc.yz <- Y2Yc[-Lp$between.idx[[2]],
        Lp$between.idx[[2]],
        drop = FALSE
      ]
      Y2Yc.yy <- Y2Yc[-Lp$between.idx[[2]],
        -Lp$between.idx[[2]],
        drop = FALSE
      ]

      # construct sigma.j
      sigma.j <- (nj * sigma.b.z) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(
        S = sigma.j,
        logdet = FALSE, Sinv.method = Sinv.method
      )
      sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
      sigma.zi.zy.ji <- t(sigma.ji.yz.zi)

      # common parts
      jYZj <- nj * (sigma.j.inv %*%
        (sigma.yz.zi %*% Y2Yc.zz %*% t(sigma.yz.zi)
          - Y2Yc.yz %*% t(sigma.yz.zi)
          - t(Y2Yc.yz %*% t(sigma.yz.zi)) + Y2Yc.yy)
        %*% sigma.j.inv)

      Z1 <- Y2Yc.zz %*% t(sigma.ji.yz.zi) %*% sigma.yz
      YZ1 <- t(Y2Yc.yz) %*% sigma.j.inv %*% sigma.yz


      # Mu.Z
      G.muz[cl, ] <- -2 * as.numeric(
        (sigma.zz.inv + nj * (sigma.zi.zy.ji %*% sigma.yz.zi)) %*% zc
          - nj * sigma.zi.zy.ji %*% yc
      )

      # MU.Y
      G.muy[cl, ] <- 2 * nj * as.numeric(zc %*% sigma.zi.zy.ji -
        yc %*% sigma.j.inv)

      # SIGMA.W
      g.sigma.w <- ((nj - 1) * sigma.w.inv
        - sigma.w.inv %*% (crossprod(Y1m) - nj * Y2Yc.yy) %*% sigma.w.inv
        + sigma.j.inv - jYZj)

      tmp <- g.sigma.w * 2
      diag(tmp) <- diag(g.sigma.w)
      G.Sigma.w[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g.sigma.b <- nj * (sigma.j.inv - jYZj)

      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.Sigma.b[cl, ] <- lav_matrix_vech(tmp)


      # SIGMA.ZZ
      g.sigma.zz <- (sigma.zz.inv + nj * sigma.zz.inv %*% (
        t(sigma.yz) %*% (sigma.j.inv - jYZj) %*% sigma.yz
          - (1 / nj * Y2Yc.zz + t(Z1) + Z1 - t(YZ1) - YZ1)) %*%
        sigma.zz.inv)

      tmp <- g.sigma.zz * 2
      diag(tmp) <- diag(g.sigma.zz)
      G.Sigma.zz[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZY
      g.sigma.yz <- 2 * nj * (
        (sigma.j.inv %*%
          (sigma.yz.zi %*% Y2Yc.zz - sigma.yz - Y2Yc.yz)
          + jYZj %*% sigma.yz) %*% sigma.zz.inv)

      G.Sigma.yz[cl, ] <- lav_matrix_vec(g.sigma.yz)
    }
  } # between.idx

  else { # no level-2 variables

    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster.size[cl]

      # data within for the cluster (centered by mu.y)
      Y1m <- Y1w.cm[cluster.idx == cl, , drop = FALSE]
      yc <- Y2.cm[cl, ]

      # data between
      Y2Yc.yy <- tcrossprod(Y2.cm[cl, ])

      # construct sigma.j
      sigma.j <- (nj * sigma.b) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(
        S = sigma.j,
        logdet = FALSE, Sinv.method = Sinv.method
      )
      # common part
      jYYj <- nj * sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv

      # MU.Y
      G.muy[cl, ] <- -2 * nj * as.numeric(yc %*% sigma.j.inv)

      # SIGMA.W
      g.sigma.w <- ((nj - 1) * sigma.w.inv
        - sigma.w.inv %*% (crossprod(Y1m) - nj * Y2Yc.yy) %*% sigma.w.inv
        + sigma.j.inv - jYYj)
      tmp <- g.sigma.w * 2
      diag(tmp) <- diag(g.sigma.w)
      G.Sigma.w[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g.sigma.b <- nj * (sigma.j.inv - jYYj)
      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.Sigma.b[cl, ] <- lav_matrix_vech(tmp)
    }
  }

  # rearrange columns to Mu.W, Mu.B, Sigma.W, Sigma.B
  ov.idx <- Lp$ov.idx
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))

  # Mu.W (for within-only)
  Mu.W.tilde <- matrix(0, nclusters, p.tilde)
  Mu.W.tilde[, ov.idx[[1]]] <- G.muy
  Mu.W.tilde[, Lp$both.idx[[2]]] <- 0 # ZERO!!!
  Mu.W <- Mu.W.tilde[, ov.idx[[1]], drop = FALSE]

  # Mu.B
  Mu.B.tilde <- matrix(0, nclusters, p.tilde)
  Mu.B.tilde[, ov.idx[[1]]] <- G.muy
  if (length(between.idx) > 0L) {
    Mu.B.tilde[, between.idx] <- G.muz
  }
  Mu.B <- Mu.B.tilde[, ov.idx[[2]], drop = FALSE]

  # Sigma.W
  Sigma.W <- G.Sigma.w

  # Sigma.B
  if (length(between.idx) > 0L) {
    p.tilde.star <- p.tilde * (p.tilde + 1) / 2
    B.tilde <- lav_matrix_vech_reverse(seq_len(p.tilde.star))

    Sigma.B.tilde <- matrix(0, nclusters, p.tilde.star)

    col.idx <- lav_matrix_vech(B.tilde[ov.idx[[1]], ov.idx[[1]],
      drop = FALSE
    ])
    Sigma.B.tilde[, col.idx] <- G.Sigma.b

    col.idx <- lav_matrix_vec(B.tilde[ov.idx[[1]], between.idx,
      drop = FALSE
    ])
    Sigma.B.tilde[, col.idx] <- G.Sigma.yz

    col.idx <- lav_matrix_vech(B.tilde[between.idx, between.idx,
      drop = FALSE
    ])
    Sigma.B.tilde[, col.idx] <- G.Sigma.zz

    col.idx <- lav_matrix_vech(B.tilde[ov.idx[[2]], ov.idx[[2]],
      drop = FALSE
    ])
    Sigma.B <- Sigma.B.tilde[, col.idx, drop = FALSE]
  } else {
    p.tilde.star <- p.tilde * (p.tilde + 1) / 2
    B.tilde <- lav_matrix_vech_reverse(seq_len(p.tilde.star))

    Sigma.B.tilde <- matrix(0, nclusters, p.tilde.star)

    col.idx <- lav_matrix_vech(B.tilde[ov.idx[[1]], ov.idx[[1]],
      drop = FALSE
    ])
    Sigma.B.tilde[, col.idx] <- G.Sigma.b

    col.idx <- lav_matrix_vech(B.tilde[ov.idx[[2]], ov.idx[[2]],
      drop = FALSE
    ])
    Sigma.B <- Sigma.B.tilde[, col.idx, drop = FALSE]
    # Sigma.B <- G.Sigma.b
  }

  SCORES <- cbind(Mu.W, Sigma.W, Mu.B, Sigma.B)

  SCORES
}


# first-order information: outer crossprod of scores per cluster
lav_mvnorm_cluster_information_firstorder <- function(Y1 = NULL,
                                                      YLp = NULL,
                                                      Lp = NULL,
                                                      Mu.W = NULL,
                                                      Sigma.W = NULL,
                                                      Mu.B = NULL,
                                                      Sigma.B = NULL,
                                                      x.idx = NULL,
                                                      divide.by.two = FALSE,
                                                      Sinv.method = "eigen") {
  N <- NROW(Y1)

  SCORES <- lav_mvnorm_cluster_scores_2l(
    Y1 = Y1,
    YLp = YLp,
    Lp = Lp,
    Mu.W = Mu.W,
    Sigma.W = Sigma.W,
    Mu.B = Mu.B,
    Sigma.B = Sigma.B,
    Sinv.method = Sinv.method
  )

  # divide by 2 (if we want scores wrt objective function)
  if (divide.by.two) {
    SCORES <- SCORES / 2
  }

  # unit information
  information <- crossprod(SCORES) / Lp$nclusters[[2]]

  # if x.idx, set rows/cols to zero
  if (length(x.idx) > 0L) {
    nw <- length(as.vector(Mu.W))
    nw.star <- nw * (nw + 1) / 2
    nb <- length(as.vector(Mu.B))
    ov.idx <- Lp$ov.idx

    x.idx.w <- which(ov.idx[[1]] %in% x.idx)
    if (length(x.idx.w) > 0L) {
      xw.idx <- c(
        x.idx.w,
        nw + lav_matrix_vech_which_idx(n = nw, idx = x.idx.w)
      )
    } else {
      xw.idx <- integer(0L)
    }
    x.idx.b <- which(ov.idx[[2]] %in% x.idx)
    if (length(x.idx.b) > 0L) {
      xb.idx <- c(
        x.idx.b,
        nb + lav_matrix_vech_which_idx(n = nb, idx = x.idx.b)
      )
    } else {
      xb.idx <- integer(0L)
    }

    all.idx <- c(xw.idx, nw + nw.star + xb.idx)

    information[all.idx, ] <- 0
    information[, all.idx] <- 0
  }

  information
}

# expected information 'h1' model
# order: mu.w within, vech(sigma.w) within, mu.b between, vech(sigma.b) between
# mu.w rows/cols that are splitted within/between are forced to zero
lav_mvnorm_cluster_information_expected <- function(Lp = NULL,
                                                    Mu.W = NULL,
                                                    Sigma.W = NULL,
                                                    Mu.B = NULL,
                                                    Sigma.B = NULL,
                                                    x.idx = integer(0L),
                                                    Sinv.method = "eigen") {
  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    Lp = Lp,
    Mu.W = Mu.W, Mu.B = Mu.B,
    Sigma.W = Sigma.W, Sigma.B = Sigma.B
  )
  mu.y <- out$mu.y
  mu.z <- out$mu.z
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz

  # create Delta.W.tilde, Delta.B.tilde
  ov.idx <- Lp$ov.idx
  nw <- length(ov.idx[[1]])
  nb <- length(ov.idx[[2]])
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))
  p.tilde.star <- p.tilde * (p.tilde + 1) / 2
  npar <- p.tilde + p.tilde.star
  B.tilde <- lav_matrix_vech_reverse(seq_len(p.tilde.star))
  w.idx <- lav_matrix_vech(B.tilde[ov.idx[[1]], ov.idx[[1]], drop = FALSE])
  b.idx <- lav_matrix_vech(B.tilde[ov.idx[[2]], ov.idx[[2]], drop = FALSE])

  Delta.W.tilde <- matrix(0, npar, npar)
  Delta.B.tilde <- matrix(0, npar, npar)
  Delta.W.tilde[
    c(ov.idx[[1]], w.idx + p.tilde),
    c(ov.idx[[1]], w.idx + p.tilde)
  ] <- diag(nw + nw * (nw + 1) / 2)
  Delta.B.tilde[
    c(ov.idx[[2]], b.idx + p.tilde),
    c(ov.idx[[2]], b.idx + p.tilde)
  ] <- diag(nb + nb * (nb + 1) / 2)
  Delta.W.tilde <- cbind(Delta.W.tilde, matrix(0, npar, npar))
  Delta.B.tilde <- cbind(matrix(0, npar, npar), Delta.B.tilde)

  nobs <- Lp$nclusters[[1]]
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.sizes <- Lp$cluster.sizes[[2]]
  ncluster.sizes <- Lp$ncluster.sizes[[2]]
  n.s <- Lp$cluster.size.ns[[2]]
  between.idx <- Lp$between.idx[[2]]

  information.j <- matrix(0, npar * 2, npar * 2)
  for (clz in seq_len(ncluster.sizes)) {
    # cluster size
    nj <- cluster.sizes[clz]

    # Delta.j -- changes per cluster(size)
    # this is why we can not write info = t(delta) info.sat delta
    Delta.j <- Delta.B.tilde + 1 / nj * Delta.W.tilde

    # compute Sigma.j
    sigma.j <- sigma.w + nj * sigma.b
    if (length(between.idx) > 0L) {
      omega.j <- matrix(0, p.tilde, p.tilde)
      omega.j[-between.idx, -between.idx] <- 1 / nj * sigma.j
      omega.j[-between.idx, between.idx] <- sigma.yz
      omega.j[between.idx, -between.idx] <- t(sigma.yz)
      omega.j[between.idx, between.idx] <- sigma.zz
      # omega.j <- rbind( cbind(sigma.zz, t(sigma.yz)),
      #                  cbind(sigma.yz, 1/nj * sigma.j) )
    } else {
      omega.j <- 1 / nj * sigma.j
    }
    omega.j.inv <- solve(omega.j)

    I11.j <- omega.j.inv
    I22.j <- 0.5 * lav_matrix_duplication_pre_post(omega.j.inv %x% omega.j.inv)
    I.j <- lav_matrix_bdiag(I11.j, I22.j)
    info.j <- t(Delta.j) %*% I.j %*% Delta.j

    information.j <- information.j + n.s[clz] * info.j
  }

  Sigma.W.inv <- lav_matrix_symmetric_inverse(
    S = Sigma.W, logdet = FALSE,
    Sinv.method = Sinv.method
  )
  # create Sigma.W.inv.tilde
  Sigma.W.inv.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.W.inv.tilde[ov.idx[[1]], ov.idx[[1]]] <- Sigma.W.inv

  I11.w <- Sigma.W.inv.tilde
  I22.w <- 0.5 * lav_matrix_duplication_pre_post(Sigma.W.inv.tilde %x% Sigma.W.inv.tilde)
  I.w <- lav_matrix_bdiag(I11.w, I22.w)
  information.w <- (nobs - nclusters) *
    (t(Delta.W.tilde) %*% I.w %*% Delta.W.tilde)

  # unit information
  information.tilde <- 1 / Lp$nclusters[[2]] * (information.w + information.j)

  # force zero for means both.idx in within part
  information.tilde[Lp$both.idx[[2]], ] <- 0
  information.tilde[, Lp$both.idx[[2]]] <- 0

  # if x.idx, set rows/cols to zero
  if (length(x.idx) > 0L) {
    xw.idx <- c(
      x.idx,
      p.tilde + lav_matrix_vech_which_idx(n = p.tilde, idx = x.idx)
    )
    xb.idx <- npar + xw.idx
    all.idx <- c(xw.idx, xb.idx)
    information.tilde[all.idx, ] <- 0
    information.tilde[, all.idx] <- 0
  }

  # remove redundant rows/cols
  ok.idx <- c(
    ov.idx[[1]],
    w.idx + p.tilde,
    npar + ov.idx[[2]],
    npar + b.idx + p.tilde
  )

  information <- information.tilde[ok.idx, ok.idx]

  information
}


# expected information -- delta
# for non-saturated models only
lav_mvnorm_cluster_information_expected_delta <- function(Lp = NULL,
                                                          Delta = NULL,
                                                          Mu.W = NULL,
                                                          Sigma.W = NULL,
                                                          Mu.B = NULL,
                                                          Sigma.B = NULL,
                                                          Sinv.method = "eigen") {
  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    Lp = Lp,
    Mu.W = Mu.W, Mu.B = Mu.B,
    Sigma.W = Sigma.W, Sigma.B = Sigma.B
  )
  mu.y <- out$mu.y
  mu.z <- out$mu.z
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz

  # Delta -- this group
  npar <- NCOL(Delta)

  # create Delta.W.tilde, Delta.B.tilde
  ov.idx <- Lp$ov.idx
  nw <- length(ov.idx[[1]])
  nw.star <- nw * (nw + 1) / 2
  nb <- length(ov.idx[[2]])

  Delta.W <- Delta[1:(nw + nw.star), , drop = FALSE]
  Delta.B <- Delta[-(1:(nw + nw.star)), , drop = FALSE]

  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))
  p.tilde.star <- p.tilde * (p.tilde + 1) / 2
  Delta.W.tilde.Mu <- matrix(0, p.tilde, npar)
  Delta.W.tilde.Sigma <- matrix(0, p.tilde.star, npar)
  Delta.B.tilde.Mu <- matrix(0, p.tilde, npar)
  Delta.B.tilde.Sigma <- matrix(0, p.tilde.star, npar)

  Delta.W.tilde.Mu[ov.idx[[1]], ] <- Delta.W[1:nw, ]
  Delta.B.tilde.Mu[ov.idx[[2]], ] <- Delta.B[1:nb, ]

  # correct Delta to reflect Mu.W[ both.idx ] is added to Mu.B[ both.idx ]
  # changed in 0.6-5
  Delta.B.tilde.Mu[Lp$both.idx[[2]], ] <-
    (Delta.B.tilde.Mu[Lp$both.idx[[2]], ] +
      Delta.W.tilde.Mu[Lp$both.idx[[2]], ])
  Delta.W.tilde.Mu[Lp$both.idx[[2]], ] <- 0


  B.tilde <- lav_matrix_vech_reverse(seq_len(p.tilde.star))
  w.idx <- lav_matrix_vech(B.tilde[ov.idx[[1]], ov.idx[[1]], drop = FALSE])
  b.idx <- lav_matrix_vech(B.tilde[ov.idx[[2]], ov.idx[[2]], drop = FALSE])
  Delta.W.tilde.Sigma[w.idx, ] <- Delta.W[-(1:nw), ]
  Delta.B.tilde.Sigma[b.idx, ] <- Delta.B[-(1:nb), ]

  Delta.W.tilde <- rbind(Delta.W.tilde.Mu, Delta.W.tilde.Sigma)
  Delta.B.tilde <- rbind(Delta.B.tilde.Mu, Delta.B.tilde.Sigma)

  nobs <- Lp$nclusters[[1]]
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.sizes <- Lp$cluster.sizes[[2]]
  ncluster.sizes <- Lp$ncluster.sizes[[2]]
  n.s <- Lp$cluster.size.ns[[2]]
  between.idx <- Lp$between.idx[[2]]

  information.j <- matrix(0, npar, npar)
  for (clz in seq_len(ncluster.sizes)) {
    # cluster size
    nj <- cluster.sizes[clz]

    # Delta.j -- changes per cluster(size)
    # this is why we can not write info = t(delta) info.sat delta
    Delta.j <- Delta.B.tilde + 1 / nj * Delta.W.tilde

    # compute Sigma.j
    sigma.j <- sigma.w + nj * sigma.b
    if (length(between.idx) > 0L) {
      omega.j <- matrix(0, p.tilde, p.tilde)
      omega.j[-between.idx, -between.idx] <- 1 / nj * sigma.j
      omega.j[-between.idx, between.idx] <- sigma.yz
      omega.j[between.idx, -between.idx] <- t(sigma.yz)
      omega.j[between.idx, between.idx] <- sigma.zz
      # omega.j <- rbind( cbind(sigma.zz, t(sigma.yz)),
      #                  cbind(sigma.yz, 1/nj * sigma.j) )
    } else {
      omega.j <- 1 / nj * sigma.j
    }
    omega.j.inv <- solve(omega.j)

    I11.j <- omega.j.inv
    I22.j <- 0.5 * lav_matrix_duplication_pre_post(omega.j.inv %x% omega.j.inv)
    I.j <- lav_matrix_bdiag(I11.j, I22.j)
    info.j <- t(Delta.j) %*% I.j %*% Delta.j

    information.j <- information.j + n.s[clz] * info.j
  }


  Sigma.W.inv <- lav_matrix_symmetric_inverse(
    S = sigma.w, logdet = FALSE,
    Sinv.method = Sinv.method
  )
  I11.w <- Sigma.W.inv
  I22.w <- 0.5 * lav_matrix_duplication_pre_post(Sigma.W.inv %x% Sigma.W.inv)
  I.w <- lav_matrix_bdiag(I11.w, I22.w)

  # force zero for means both.idx in within part
  # changed in 0.6-5
  I.w[Lp$both.idx[[2]], ] <- 0
  I.w[, Lp$both.idx[[2]]] <- 0

  information.w <- (nobs - nclusters) * (t(Delta.W) %*% I.w %*% Delta.W)

  # unit information
  information <- 1 / Lp$nclusters[[2]] * (information.w + information.j)


  information
}


# observed information
# order: mu.w within, vech(sigma.w) within, mu.b between, vech(sigma.b) between
# mu.w rows/cols that are splitted within/between are forced to zero
#
# numerical approximation (for now)
lav_mvnorm_cluster_information_observed <- function(Lp = NULL,
                                                    YLp = NULL,
                                                    Mu.W = NULL,
                                                    Sigma.W = NULL,
                                                    Mu.B = NULL,
                                                    Sigma.B = NULL,
                                                    x.idx = integer(0L),
                                                    Sinv.method = "eigen") {
  nobs <- Lp$nclusters[[1]]

  nw <- length(as.vector(Mu.W))
  nw.star <- nw * (nw + 1) / 2
  nb <- length(as.vector(Mu.B))
  nb.star <- nb * (nb + 1) / 2

  ov.idx <- Lp$ov.idx
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))

  # Mu.W (for within-only)
  Mu.W.tilde <- numeric(p.tilde)
  Mu.W.tilde[ov.idx[[1]]] <- Mu.W

  # local function -- gradient
  GRAD <- function(x) {
    # Mu.W (for within-only)
    Mu.W.tilde2 <- numeric(p.tilde)
    Mu.W.tilde2[ov.idx[[1]]] <- x[1:nw]
    Mu.W.tilde2[Lp$both.idx[[2]]] <- Mu.W.tilde[Lp$both.idx[[2]]]
    Mu.W2 <- Mu.W.tilde2[ov.idx[[1]]]

    Sigma.W2 <- lav_matrix_vech_reverse(x[nw + 1:nw.star])
    Mu.B2 <- x[nw + nw.star + 1:nb]
    Sigma.B2 <- lav_matrix_vech_reverse(x[nw + nw.star + nb + 1:nb.star])

    dx <- lav_mvnorm_cluster_dlogl_2l_samplestats(
      YLp = YLp,
      Lp = Lp, Mu.W = Mu.W2, Sigma.W = Sigma.W2,
      Mu.B = Mu.B2, Sigma.B = Sigma.B2,
      return.list = FALSE,
      Sinv.method = Sinv.method
    )

    # dx is for -2*logl
    -1 / 2 * dx
  }

  # start.x
  start.x <- c(
    as.vector(Mu.W), lav_matrix_vech(Sigma.W),
    as.vector(Mu.B), lav_matrix_vech(Sigma.B)
  )

  # total information
  information <- -1 * numDeriv::jacobian(func = GRAD, x = start.x)

  # unit information
  information <- information / Lp$nclusters[[2]]

  # if x.idx, set rows/cols to zero
  if (length(x.idx) > 0L) {
    x.idx.w <- which(ov.idx[[1]] %in% x.idx)
    if (length(x.idx.w) > 0L) {
      xw.idx <- c(
        x.idx.w,
        nw + lav_matrix_vech_which_idx(n = nw, idx = x.idx.w)
      )
    } else {
      xw.idx <- integer(0L)
    }
    x.idx.b <- which(ov.idx[[2]] %in% x.idx)
    if (length(x.idx.b) > 0L) {
      xb.idx <- c(
        x.idx.b,
        nb + lav_matrix_vech_which_idx(n = nb, idx = x.idx.b)
      )
    } else {
      xb.idx <- integer(0L)
    }

    all.idx <- c(xw.idx, nw + nw.star + xb.idx)

    information[all.idx, ] <- 0
    information[, all.idx] <- 0
  }

  information
}

# estimate ML estimates of Mu.W, Mu.B, Sigma.W, Sigma.B
# using the EM algorithm
#
# per cluster-SIZE
#
lav_mvnorm_cluster_em_sat <- function(YLp = NULL,
                                      Lp = NULL,
                                      verbose = TRUE,
                                      tol = 1e-04,
                                      max.iter = 5000,
                                      min.variance = 1e-05) {
  # lavdata
  between.idx <- Lp$between.idx[[2]]
  within.idx <- Lp$within.idx[[2]]
  Y2 <- YLp[[2]]$Y2

  # starting values for Sigma
  ov.idx <- Lp$ov.idx
  # COVT <- lavsamplestats@cov[[1]]
  # Sigma.W <- diag( diag(COVT)[ov.idx[[1]]] )
  # Sigma.B <- diag( diag(COVT)[ov.idx[[2]]] )
  Sigma.W <- diag(length(ov.idx[[1]]))
  Sigma.B <- diag(length(ov.idx[[2]]))
  Mu.W <- numeric(length(ov.idx[[1]]))
  Mu.B <- numeric(length(ov.idx[[2]]))
  # Mu.W.tilde <- YLp[[2]]$Mu.W
  # Mu.B.tilde <- YLp[[2]]$Mu.B
  # if(length(between.idx) > 0) {
  #    Mu.W <- Mu.W.tilde[-between.idx]
  # } else {
  #    Mu.W <- Mu.W.tilde
  # }
  # if(length(within.idx) > 0) {
  #    Mu.B <- Mu.B.tilde[-within.idx]
  # } else {
  #    Mu.B <- Mu.B.tilde
  # }

  # report initial fx
  fx <- lav_mvnorm_cluster_loglik_samplestats_2l(
    YLp = YLp, Lp = Lp,
    Mu.W = Mu.W, Sigma.W = Sigma.W,
    Mu.B = Mu.B, Sigma.B = Sigma.B,
    Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE
  )

  # if verbose, report
  if (verbose) {
    cat(
      "EM iter:", sprintf("%3d", 0),
      " fx =", sprintf("%17.10f", fx),
      "\n"
    )
  }

  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    Lp = Lp,
    Mu.W = Mu.W, Sigma.W = Sigma.W, Mu.B = Mu.B, Sigma.B = Sigma.B
  )
  mu.y <- out$mu.y
  mu.z <- out$mu.z
  mu.w <- out$mu.w
  mu.b <- out$mu.b
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz

  # mu.z and sigma.zz can be computed beforehand
  if (length(between.idx) > 0L) {
    Z <- Y2[, between.idx, drop = FALSE]
    mu.z <- colMeans(Z, na.rm = TRUE)
    sigma.zz <- cov(Z, use = "pairwise.complete.obs") * (Lp$nclusters[[2]] - 1L) / Lp$nclusters[[2]]
    # sigma.zz <- 1/Lp$nclusters[[2]] * crossprod(Z) - tcrossprod(mu.z)
    # Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop=FALSE]
  }

  # EM iterations
  fx.old <- fx
  for (i in 1:max.iter) {
    # E-step
    estep <- lav_mvnorm_cluster_em_estepb( # Y1 = Y1,
      YLp = YLp,
      Lp = Lp,
      sigma.w = sigma.w,
      sigma.b = sigma.b,
      mu.w = mu.w,
      mu.b = mu.b,
      sigma.yz = sigma.yz,
      sigma.zz = sigma.zz,
      mu.z = mu.z
    )

    # mstep
    sigma.w <- estep$sigma.w
    sigma.b <- estep$sigma.b
    sigma.yz <- estep$sigma.yz
    mu.w <- estep$mu.w
    mu.b <- estep$mu.b

    implied2 <- lav_mvnorm_cluster_2l2implied(
      Lp = Lp,
      sigma.w = estep$sigma.w, sigma.b = estep$sigma.b,
      sigma.zz = sigma.zz, sigma.yz = estep$sigma.yz,
      mu.z = mu.z,
      mu.y = NULL, mu.w = estep$mu.w, mu.b = estep$mu.b
    )

    # check for (near-zero) variances at the within level, and set
    # them to min.variance
    Sigma.W <- implied2$Sigma.W
    zero.var <- which(diag(Sigma.W) < min.variance)
    if (length(zero.var) > 0L) {
      Sigma.W[, zero.var] <- sigma.w[, zero.var] <- 0
      Sigma.W[zero.var, ] <- sigma.w[zero.var, ] <- 0
      diag(Sigma.W)[zero.var] <- diag(sigma.w)[zero.var] <- min.variance
    }

    fx <- lav_mvnorm_cluster_loglik_samplestats_2l(
      YLp = YLp,
      Lp = Lp, Mu.W = implied2$Mu.W, Sigma.W = Sigma.W,
      Mu.B = implied2$Mu.B, Sigma.B = implied2$Sigma.B,
      Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE
    )

    # fx.delta
    fx.delta <- fx - fx.old

    # what if fx.delta is negative?
    if (fx.delta < 0) {
      lav_msg_warn(gettext(
        "logl decreased during EM steps of the saturated (H1) model"))
    }

    if (verbose) {
      cat(
        "EM iter:", sprintf("%3d", i),
        " fx =", sprintf("%17.10f", fx),
        " fx.delta =", sprintf("%9.8f", fx.delta),
        "\n"
      )
    }

    # convergence check
    if (fx.delta < tol) {
      break
    } else {
      fx.old <- fx
    }
  } # EM iterations

  list(
    Sigma.W = implied2$Sigma.W, Sigma.B = implied2$Sigma.B,
    Mu.W = implied2$Mu.W, Mu.B = implied2$Mu.B, logl = fx
  )
}


# based on lav_mvnorm_cluster_em_estep
lav_mvnorm_cluster_em_h0 <- function(lavsamplestats = NULL,
                                     lavdata = NULL,
                                     lavimplied = NULL,
                                     lavpartable = NULL,
                                     lavmodel = NULL,
                                     lavoptions = NULL,
                                     verbose = FALSE,
                                     verbose.x = FALSE,
                                     fx.tol = 1e-08,
                                     dx.tol = 1e-05,
                                     max.iter = 5000,
                                     mstep.iter.max = 10000L,
                                     mstep.rel.tol = 1e-10) {
  # single group only for now
  stopifnot(lavdata@ngroups == 1L)

  # lavdata
  Lp <- lavdata@Lp[[1]] # first group only (for now)
  ov.names.l <- lavdata@ov.names.l[[1]] # first group only (for now)
  Y1 <- lavdata@X[[1]] # first group only
  YLp <- lavsamplestats@YLp[[1]] # first group only

  between.idx <- Lp$between.idx[[2]]
  Y2 <- YLp[[2]]$Y2

  # initial values
  x.current <- lav_model_get_parameters(lavmodel)

  # implied
  if (is.null(lavimplied)) {
    lavimplied <- lav_model_implied(lavmodel)
  }

  # TODO: what if current 'starting' parameters imply a non-pd sigma.b?

  # report initial fx
  fx <- lav_mvnorm_cluster_loglik_samplestats_2l(
    YLp = YLp, Lp = Lp,
    Mu.W = lavimplied$mean[[1]], Sigma.W = lavimplied$cov[[1]],
    Mu.B = lavimplied$mean[[2]], Sigma.B = lavimplied$cov[[2]],
    Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE
  )

  # if verbose, report
  if (verbose) {
    cat(
      "EM iter:", sprintf("%3d", 0),
      " fx =", sprintf("%17.10f", fx),
      "\n"
    )
  }

  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    Lp = Lp,
    Mu.W = lavimplied$mean[[1]], Sigma.W = lavimplied$cov[[1]],
    Mu.B = lavimplied$mean[[2]], Sigma.B = lavimplied$cov[[2]]
  )
  mu.y <- out$mu.y
  mu.z <- out$mu.z
  mu.w <- out$mu.w
  mu.b <- out$mu.b
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz

  # mu.z and sigma.zz can be computed beforehand
  if (length(between.idx) > 0L) {
    Z <- Y2[, between.idx, drop = FALSE]
    mu.z <- colMeans(Y2)[between.idx]
    sigma.zz <- cov(Z) * (Lp$nclusters[[2]] - 1L) / Lp$nclusters[[2]]
    # sigma.zz <- 1/Lp$nclusters[[2]] * crossprod(Z) - tcrossprod(mu.z)
    # Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop=FALSE]
  }

  # EM iterations
  fx.old <- fx
  fx2.old <- 0
  REL <- numeric(max.iter)
  for (i in 1:max.iter) {
    # E-step
    estep <- lav_mvnorm_cluster_em_estepb(
      YLp = YLp,
      Lp = Lp,
      sigma.w = sigma.w,
      sigma.b = sigma.b,
      mu.w = mu.w,
      mu.b = mu.b,
      sigma.yz = sigma.yz,
      sigma.zz = sigma.zz,
      mu.z = mu.z
    )

    # back to model-implied dimensions
    implied <- lav_mvnorm_cluster_2l2implied(
      Lp = Lp,
      sigma.w = estep$sigma.w, sigma.b = estep$sigma.b,
      sigma.zz = sigma.zz, sigma.yz = estep$sigma.yz,
      mu.z = mu.z,
      mu.y = NULL, mu.w = estep$mu.w, mu.b = estep$mu.b
    )
    rownames(implied$Sigma.W) <- ov.names.l[[1]]
    rownames(implied$Sigma.B) <- ov.names.l[[2]]

    # M-step

    # fit two-group model
    local.partable <- lavpartable
    # if a group column exists, delete it (it will be overriden anyway)
    local.partable$group <- NULL
    level.idx <- which(names(local.partable) == "level")
    names(local.partable)[level.idx] <- "group"
    local.partable$est <- NULL
    local.partable$se <- NULL

    # give current values as starting values
    free.idx <- which(lavpartable$free > 0L)
    local.partable$ustart[free.idx] <- x.current

    local.fit <- lavaan(local.partable,
      sample.cov = list(
        within = implied$Sigma.W,
        between = implied$Sigma.B
      ),
      sample.mean = list(
        within = implied$Mu.W,
        between = implied$Mu.B
      ),
      sample.nobs = Lp$nclusters,
      sample.cov.rescale = FALSE,
      control = list(
        iter.max = mstep.iter.max,
        rel.tol = mstep.rel.tol
      ),
      fixed.x = any(lavpartable$exo == 1L),
      estimator = "ML",
      warn = FALSE, # no warnings
      check.start = FALSE,
      check.post = FALSE,
      check.gradient = FALSE,
      check.vcov = FALSE,
      baseline = FALSE,
      h1 = FALSE,
      se = "none",
      test = "none"
    )

    # end of M-step

    implied2 <- local.fit@implied
    fx <- lav_mvnorm_cluster_loglik_samplestats_2l(
      YLp = YLp,
      Lp = Lp, Mu.W = implied2$mean[[1]], Sigma.W = implied2$cov[[1]],
      Mu.B = implied2$mean[[2]], Sigma.B = implied2$cov[[2]],
      Sinv.method = "eigen", log2pi = TRUE, minus.two = FALSE
    )

    # fx.delta
    fx.delta <- fx - fx.old

    # derivatives
    lavmodel <- lav_model_set_parameters(lavmodel, x = local.fit@optim$x)
    dx <- lav_model_gradient(lavmodel,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats
    )
    max.dx <- max(abs(dx))

    if (verbose) {
      cat(
        "EM iter:", sprintf("%3d", i),
        " fx =", sprintf("%17.10f", fx),
        " fx.delta =", sprintf("%9.8f", fx.delta),
        " mstep.iter =", sprintf(
          "%3d",
          lavInspect(local.fit, "iterations")
        ),
        " max.dx = ", sprintf("%9.8f", max.dx),
        "\n"
      )
    }

    # stopping rule check
    if (fx.delta < fx.tol) {
      if (verbose) {
        cat("EM stopping rule reached: fx.delta < ", fx.tol, "\n")
      }
      break
    } else {
      fx.old <- fx
      x.current <- local.fit@optim$x
      if (verbose.x) {
        print(round(x.current, 3))
      }
    }

    # second stopping rule check -- derivatives
    if (max.dx < dx.tol) {
      if (verbose) {
        cat("EM stopping rule reached: max.dx < ", dx.tol, "\n")
      }
      break
    }

    # translate to internal matrices
    out <- lav_mvnorm_cluster_implied22l(
      Lp = Lp,
      Mu.W = implied2$mean[[1]], Sigma.W = implied2$cov[[1]],
      Mu.B = implied2$mean[[2]], Sigma.B = implied2$cov[[2]]
    )
    mu.y <- out$mu.y
    mu.z <- out$mu.z
    mu.w <- out$mu.w
    mu.b <- out$mu.b
    sigma.w <- out$sigma.w
    sigma.b <- out$sigma.b
    sigma.zz <- out$sigma.zz
    sigma.yz <- out$sigma.yz
  } # EM iterations

  x <- local.fit@optim$x

  # add attributes
  if (i < max.iter) {
    attr(x, "converged") <- TRUE
    attr(x, "warn.txt") <- ""
  } else {
    attr(x, "converged") <- FALSE
    attr(x, "warn.txt") <- paste("maxmimum number of iterations (",
      max.iter, ") ",
      "was reached without convergence.\n",
      sep = ""
    )
  }
  attr(x, "iterations") <- i
  attr(x, "control") <- list(
    em.iter.max = max.iter,
    em.fx.tol = fx.tol,
    em.dx.tol = dx.tol
  )
  attr(fx, "fx.group") <- fx # single group for now
  attr(x, "fx") <- fx

  x
}

# get the random effects (here: expected values for cluster means)
# and optionally a standard error
lav_mvnorm_cluster_em_estep_ranef <- function(YLp = NULL,
                                              Lp = NULL,
                                              sigma.w = NULL,
                                              sigma.b = NULL,
                                              sigma.yz = NULL,
                                              sigma.zz = NULL,
                                              mu.z = NULL,
                                              mu.w = NULL,
                                              mu.b = NULL,
                                              se = FALSE) {
  # sample stats
  nobs <- Lp$nclusters[[1]]
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  between.idx <- Lp$between.idx[[2]]

  Y2 <- YLp[[2]]$Y2

  nvar.y <- ncol(sigma.w)
  nvar.z <- ncol(sigma.zz)

  MB.j <- matrix(0, nrow = nclusters, ncol = nvar.y)
  SE.j <- matrix(0, nrow = nclusters, ncol = nvar.y)

  mu.y <- mu.w + mu.b

  if (length(between.idx) > 0L) {
    sigma.1 <- cbind(sigma.yz, sigma.b)
    mu <- c(mu.z, mu.y)
  } else {
    sigma.1 <- sigma.b
    mu <- mu.y
  }

  # E-step
  for (cl in seq_len(nclusters)) {
    nj <- cluster.size[cl]

    # data
    if (length(between.idx) > 0L) {
      # z comes first!
      b.j <- c(
        Y2[cl, between.idx],
        Y2[cl, -between.idx]
      )
      ybar.j <- Y2[cl, -between.idx]
    } else {
      ybar.j <- b.j <- Y2[cl, ]
    }

    sigma.j <- sigma.w + nj * sigma.b
    if (length(between.idx) > 0L) {
      omega.j <- rbind(
        cbind(sigma.zz, t(sigma.yz)),
        cbind(sigma.yz, 1 / nj * sigma.j)
      )
    } else {
      omega.j <- 1 / nj * sigma.j
    }
    omega.j.inv <- solve(omega.j)

    # E(v|y)
    Ev <- as.numeric(mu.b + (sigma.1 %*% omega.j.inv %*% (b.j - mu)))
    MB.j[cl, ] <- Ev

    if (se) {
      # Cov(v|y)
      Covv <- sigma.b - (sigma.1 %*% omega.j.inv %*% t(sigma.1))

      # force symmetry
      Covv <- (Covv + t(Covv)) / 2

      Covv.diag <- diag(Covv)
      nonzero.idx <- which(Covv.diag > 0)

      SE.j[cl, ] <- numeric(length(Covv.diag))
      SE.j[cl, nonzero.idx] <- sqrt(Covv.diag[nonzero.idx])
    }
  }

  if (se) {
    attr(MB.j, "se") <- SE.j
  }

  MB.j
}

# per cluster
lav_mvnorm_cluster_em_estep <- function( # Y1           = NULL,
                                        YLp = NULL,
                                        Lp = NULL,
                                        sigma.w = NULL,
                                        sigma.b = NULL,
                                        sigma.yz = NULL,
                                        sigma.zz = NULL,
                                        mu.z = NULL,
                                        mu.w = NULL,
                                        mu.b = NULL) {
  # sample stats
  nobs <- Lp$nclusters[[1]]
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]
  within.idx <- Lp$within.idx[[2]]
  between.idx <- Lp$between.idx[[2]]
  both.idx <- Lp$both.idx[[2]]

  Y2 <- YLp[[2]]$Y2
  Y1Y1 <- YLp[[2]]$Y1Y1

  nvar.y <- ncol(sigma.w)
  nvar.z <- ncol(sigma.zz)

  CW2.j <- matrix(0, nrow = nvar.y, ncol = nvar.y)
  CB.j <- matrix(0, nrow = nvar.y, ncol = nvar.y)
  MW.j <- matrix(0, nrow = nclusters, ncol = nvar.y)
  MB.j <- matrix(0, nrow = nclusters, ncol = nvar.y)
  ZY.j <- matrix(0, nrow = nvar.z, ncol = nvar.y)

  mu.y <- mu.w + mu.b

  if (length(between.idx) > 0L) {
    sigma.1 <- cbind(sigma.yz, sigma.b)
    mu <- c(mu.z, mu.y)
    Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop = FALSE]
  } else {
    sigma.1 <- sigma.b
    mu <- mu.y
  }

  # E-step
  for (cl in seq_len(nclusters)) {
    nj <- cluster.size[cl]

    # data
    if (length(between.idx) > 0L) {
      # z comes first!
      b.j <- c(
        Y2[cl, between.idx],
        Y2[cl, -between.idx]
      )
      ybar.j <- Y2[cl, -between.idx]
    } else {
      ybar.j <- b.j <- Y2[cl, ]
    }

    sigma.j <- sigma.w + nj * sigma.b
    if (length(between.idx) > 0L) {
      omega.j <- rbind(
        cbind(sigma.zz, t(sigma.yz)),
        cbind(sigma.yz, 1 / nj * sigma.j)
      )
    } else {
      omega.j <- 1 / nj * sigma.j
    }
    omega.j.inv <- solve(omega.j)

    # E(v|y)
    Ev <- as.numeric(mu.b + (sigma.1 %*% omega.j.inv %*% (b.j - mu)))

    # Cov(v|y)
    Covv <- sigma.b - (sigma.1 %*% omega.j.inv %*% t(sigma.1))

    # force symmetry
    Covv <- (Covv + t(Covv)) / 2

    # E(vv|y) = Cov(v|y) + E(v|y)E(v|y)^T
    Evv <- Covv + tcrossprod(Ev)

    # store for this cluster
    MW.j[cl, ] <- ybar.j - Ev
    MB.j[cl, ] <- Ev
    CW2.j <- CW2.j + nj * (Evv - tcrossprod(ybar.j, Ev)
      - tcrossprod(Ev, ybar.j))
    CB.j <- CB.j + Evv

    # between only
    if (length(between.idx) > 0L) {
      ZY.j <- ZY.j + tcrossprod(Y2[cl, between.idx], Ev)
    }
  }

  M.w <- 1 / nobs * colSums(MW.j * cluster.size)
  M.b <- 1 / nclusters * colSums(MB.j)
  C.b <- 1 / nclusters * CB.j
  C.w <- 1 / nobs * (Y1Y1 + CW2.j)
  # end of E-step

  # make symmetric (not needed here?)
  # C.b <- (C.b + t(C.b))/2
  # C.w <- (C.w + t(C.w))/2

  # between only
  if (length(between.idx) > 0L) {
    A <- 1 / nclusters * ZY.j - tcrossprod(mu.z, M.b)
  }

  sigma.w <- C.w - tcrossprod(M.w)
  sigma.b <- C.b - tcrossprod(M.b)
  mu.w <- M.w
  mu.b <- M.b

  if (length(between.idx) > 0L) {
    sigma.yz <- t(A)
  }

  list(
    sigma.w = sigma.w, sigma.b = sigma.b, mu.w = mu.w, mu.b = mu.b,
    sigma.yz = sigma.yz, sigma.zz = sigma.zz, mu.z = mu.z
  )
}

# per cluster SIZE
lav_mvnorm_cluster_em_estepb <- function( # Y1           = NULL, # not used!
                                         YLp = NULL,
                                         Lp = NULL,
                                         sigma.w = NULL,
                                         sigma.b = NULL,
                                         sigma.yz = NULL,
                                         sigma.zz = NULL,
                                         mu.z = NULL,
                                         mu.w = NULL,
                                         mu.b = NULL) {
  # sample stats
  nobs <- Lp$nclusters[[1]]
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]
  between.idx <- Lp$between.idx[[2]]
  cluster.sizes <- Lp$cluster.sizes[[2]]
  ncluster.sizes <- Lp$ncluster.sizes[[2]]
  n.s <- Lp$cluster.size.ns[[2]]

  Y2 <- YLp[[2]]$Y2
  Y1Y1 <- YLp[[2]]$Y1Y1

  nvar.y <- ncol(sigma.w)
  nvar.z <- ncol(sigma.zz)

  mu.y <- mu.w + mu.b

  if (length(between.idx) > 0L) {
    sigma.1 <- cbind(sigma.yz, sigma.b)
    mu <- c(mu.z, mu.y)
    Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop = FALSE]
  } else {
    sigma.1 <- sigma.b
    mu <- mu.y
  }

  # per cluster SIZE
  CW2.s <- matrix(0, nrow = nvar.y, ncol = nvar.y)
  CB.s <- matrix(0, nrow = nvar.y, ncol = nvar.y)
  MW.s <- matrix(0, nrow = ncluster.sizes, ncol = nvar.y)
  MB.s <- matrix(0, nrow = ncluster.sizes, ncol = nvar.y)
  ZY.s <- matrix(0, nvar.z, nvar.y)

  # E-step
  for (clz in seq_len(ncluster.sizes)) {
    # cluster size
    nj <- cluster.sizes[clz]

    # data
    if (length(between.idx) > 0L) {
      # z comes first!
      b.j <- cbind(
        Y2[cluster.size == nj, between.idx, drop = FALSE],
        Y2[cluster.size == nj, -between.idx, drop = FALSE]
      )
      ybar.j <- Y2[cluster.size == nj, -between.idx, drop = FALSE]
    } else {
      ybar.j <- b.j <- Y2[cluster.size == nj, , drop = FALSE]
    }

    sigma.j <- sigma.w + nj * sigma.b
    if (length(between.idx) > 0L) {
      omega.j <- rbind(
        cbind(sigma.zz, t(sigma.yz)),
        cbind(sigma.yz, 1 / nj * sigma.j)
      )
    } else {
      omega.j <- 1 / nj * sigma.j
    }
    omega.j.inv <- solve(omega.j)
    sigma.1.j.inv <- sigma.1 %*% omega.j.inv

    # E(v|y)
    b.jc <- t(t(b.j) - mu)
    tmp <- b.jc %*% t(sigma.1.j.inv)
    Ev <- t(t(tmp) + mu.b)

    # Cov(v|y)
    Covv <- n.s[clz] * (sigma.b - (sigma.1.j.inv %*% t(sigma.1)))

    # force symmetry
    Covv <- (Covv + t(Covv)) / 2

    # E(vv|y) = Cov(v|y) + E(v|y)E(v|y)^T
    Evv <- Covv + crossprod(Ev)

    # store for this cluster SIZE
    MW.s[clz, ] <- nj * colSums(ybar.j - Ev)
    MB.s[clz, ] <- colSums(Ev)
    CW2.s <- CW2.s + nj * (Evv - crossprod(ybar.j, Ev)
      - crossprod(Ev, ybar.j))
    CB.s <- CB.s + Evv

    # between only
    if (length(between.idx) > 0L) {
      ZY.s <- ZY.s + crossprod(Y2[cluster.size == nj, between.idx,
        drop = FALSE
      ], Ev)
    }
  } # cluster-sizes

  M.ws <- 1 / nobs * colSums(MW.s)
  M.bs <- 1 / nclusters * colSums(MB.s)
  C.bs <- 1 / nclusters * CB.s
  C.ws <- 1 / nobs * (Y1Y1 + CW2.s)

  # between only
  if (length(between.idx) > 0L) {
    As <- 1 / nclusters * ZY.s - tcrossprod(mu.z, M.bs)
  }

  sigma.w <- C.ws - tcrossprod(M.ws)
  sigma.b <- C.bs - tcrossprod(M.bs)
  mu.w <- M.ws
  mu.b <- M.bs
  if (length(between.idx) > 0L) {
    sigma.yz <- t(As)
  }

  list(
    sigma.w = sigma.w, sigma.b = sigma.b, mu.w = mu.w, mu.b = mu.b,
    sigma.yz = sigma.yz, sigma.zz = sigma.zz, mu.z = mu.z
  )
}
