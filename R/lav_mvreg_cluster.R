# loglikelihood clustered/twolevel data -- conditional.x = TRUE

# YR: first version around Sept 2021

# take model-implied mean+variance matrices, and reorder/augment them
# to facilitate computing of (log)likelihood in the two-level case

# when conditional.x = TRUE:
# - sigma.w and sigma.b: same dimensions, level-1 'Y' variables only
# - sigma.zz: level-2 variables only
# - sigma.yz: cov(level-1, level-2)
# - beta.w: beta y within part
# - beta.b: beta y between part
# - beta.z: beta z (between-only)
lav_mvreg_cluster_implied22l <- function(Lp = NULL,
                                         implied = NULL,
                                         Res.Int.W = NULL,
                                         Res.Int.B = NULL,
                                         Res.Pi.W = NULL,
                                         Res.Pi.B = NULL,
                                         Res.Sigma.W = NULL,
                                         Res.Sigma.B = NULL) {
  if (!is.null(implied)) {
    # FIXME: only for single-group analysis!
    Res.Sigma.W <- implied$res.cov[[1]]
    Res.Int.W <- implied$res.int[[1]]
    Res.Pi.W <- implied$res.slopes[[1]]

    Res.Sigma.B <- implied$res.cov[[2]]
    Res.Int.B <- implied$res.int[[2]]
    Res.Pi.B <- implied$res.slopes[[2]]
  }

  # within/between idx
  within.x.idx <- Lp$within.x.idx[[1]]
  between.y.idx <- Lp$between.y.idx[[2]]
  between.x.idx <- Lp$between.x.idx[[2]]

  # ov.idx per level
  ov.idx <- Lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))

  # only 'y'
  ov.y.idx <- Lp$ov.y.idx

  # two levels only (for now)
  ov.y.idx1 <- ov.y.idx[[1]]
  ov.y.idx2 <- ov.y.idx[[2]]

  # Sigma.W.tilde
  Sigma.W.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.W.tilde[ov.y.idx1, ov.y.idx1] <- Res.Sigma.W

  # INT.W.tilde
  INT.W.tilde <- matrix(0, p.tilde, 1L)
  INT.W.tilde[ov.y.idx1, 1L] <- Res.Int.W

  # PI.W.tilde
  PI.W.tilde <- matrix(0, p.tilde, ncol(Res.Pi.W))
  PI.W.tilde[ov.y.idx1, ] <- Res.Pi.W

  BETA.W.tilde <- rbind(t(INT.W.tilde), t(PI.W.tilde))



  # Sigma.B.tilde
  Sigma.B.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.B.tilde[ov.y.idx2, ov.y.idx2] <- Res.Sigma.B

  # INT.B.tilde
  INT.B.tilde <- matrix(0, p.tilde, 1L)
  INT.B.tilde[ov.y.idx2, 1L] <- Res.Int.B

  # PI.B.tilde
  PI.B.tilde <- matrix(0, p.tilde, ncol(Res.Pi.B))
  PI.B.tilde[ov.y.idx2, ] <- Res.Pi.B

  BETA.B.tilde <- rbind(t(INT.B.tilde), t(PI.B.tilde))

  if (length(between.y.idx) > 0L) {
    rm.idx <- c(within.x.idx, between.x.idx, between.y.idx) # between AND x
    beta.z <- BETA.B.tilde[, between.y.idx, drop = FALSE]
    beta.b <- BETA.B.tilde[, -rm.idx, drop = FALSE]
    beta.w <- BETA.W.tilde[, -rm.idx, drop = FALSE]
    sigma.zz <- Sigma.B.tilde[between.y.idx, between.y.idx, drop = FALSE]
    sigma.yz <- Sigma.B.tilde[-rm.idx, between.y.idx, drop = FALSE]
    sigma.b <- Sigma.B.tilde[-rm.idx, -rm.idx, drop = FALSE]
    sigma.w <- Sigma.W.tilde[-rm.idx, -rm.idx, drop = FALSE]
  } else {
    rm.idx <- c(within.x.idx, between.x.idx) # all 'x'
    beta.z <- matrix(0, 0L, 0L)
    sigma.zz <- matrix(0, 0L, 0L)
    beta.b <- BETA.B.tilde[, -rm.idx, drop = FALSE]
    beta.w <- BETA.W.tilde[, -rm.idx, drop = FALSE]
    sigma.b <- Sigma.B.tilde[-rm.idx, -rm.idx, drop = FALSE]
    sigma.w <- Sigma.W.tilde[-rm.idx, -rm.idx, drop = FALSE]
    sigma.yz <- matrix(0, nrow(sigma.w), 0L)
  }


  # beta.wb # FIXme: not correct if some 'x' are splitted (overlap)
  # but because we ALWAYS treat splitted-x as 'y', this is not a problem
  beta.wb <- rbind(beta.w, beta.b[-1, , drop = FALSE])
  beta.wb[1, ] <- beta.wb[1, , drop = FALSE] + beta.b[1, , drop = FALSE]

  list(
    sigma.w = sigma.w, sigma.b = sigma.b, sigma.zz = sigma.zz,
    sigma.yz = sigma.yz, beta.w = beta.w, beta.b = beta.b, beta.z = beta.z,
    beta.wb = beta.wb
  )
}


# recreate implied matrices from 2L matrices
lav_mvreg_cluster_2l2implied <- function(Lp,
                                         sigma.w = NULL,
                                         sigma.b = NULL,
                                         sigma.zz = NULL,
                                         sigma.yz = NULL,
                                         beta.w = NULL,
                                         beta.b = NULL,
                                         beta.z = NULL) {
  # within/between idx
  within.x.idx <- Lp$within.x.idx[[1]]
  between.y.idx <- Lp$between.y.idx[[2]]
  between.x.idx <- Lp$between.x.idx[[2]]

  # ov.idx per level
  ov.idx <- Lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))

  # only 'y'
  ov.y.idx <- Lp$ov.y.idx

  # two levels only (for now)
  ov.y.idx1 <- ov.y.idx[[1]]
  ov.y.idx2 <- ov.y.idx[[2]]

  # Sigma.W.tilde
  Sigma.W.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.W.tilde[ov.y.idx1, ov.y.idx1] <- sigma.w

  # INT.W.tilde
  INT.W.tilde <- matrix(0, p.tilde, 1L)
  INT.W.tilde[ov.y.idx1, 1L] <- beta.w[1L, ]

  # PI.W.tilde
  PI.W.tilde <- matrix(0, p.tilde, nrow(beta.w) - 1L)
  PI.W.tilde[ov.y.idx1, ] <- t(beta.w[-1L, ])

  # Sigma.B.tilde
  Sigma.B.tilde <- matrix(0, p.tilde, p.tilde)
  Sigma.B.tilde[ov.y.idx1, ov.y.idx1] <- sigma.b

  # INT.B.tilde
  INT.B.tilde <- matrix(0, p.tilde, 1L)
  INT.B.tilde[ov.y.idx1, 1L] <- beta.b[1L, ]

  # PI.B.tilde
  PI.B.tilde <- matrix(0, p.tilde, nrow(beta.b) - 1L)
  PI.B.tilde[ov.y.idx1, ] <- t(beta.b[-1L, ])

  if (length(between.y.idx) > 0L) {
    INT.B.tilde[between.y.idx, 1L] <- beta.z[1L, ]
    PI.B.tilde[between.y.idx, ] <- t(beta.z[-1L, ])
    Sigma.B.tilde[between.y.idx, between.y.idx] <- sigma.zz
    Sigma.B.tilde[ov.y.idx1, between.y.idx] <- sigma.yz
    Sigma.B.tilde[between.y.idx, ov.y.idx1] <- t(sigma.yz)
  }

  Res.Sigma.W <- Sigma.W.tilde[ov.y.idx1, ov.y.idx1, drop = FALSE]
  Res.Int.W <- INT.W.tilde[ov.y.idx1, , drop = FALSE]
  Res.Pi.W <- PI.W.tilde[ov.y.idx1, , drop = FALSE]

  Res.Sigma.B <- Sigma.B.tilde[ov.y.idx2, ov.y.idx2, drop = FALSE]
  Res.Int.B <- INT.B.tilde[ov.y.idx2, , drop = FALSE]
  Res.Pi.B <- PI.B.tilde[ov.y.idx2, , drop = FALSE]

  implied <- list(
    res.cov = list(Res.Sigma.W, Res.Sigma.B),
    res.int = list(Res.Int.W, Res.Int.B),
    res.slopes = list(Res.Pi.W, Res.Pi.B)
  )

  # Note: cov.x and mean.x must be added by the caller
  implied
}

lav_mvreg_cluster_loglik_samplestats_2l <- function(YLp = NULL,
                                                    Lp = NULL,
                                                    Res.Sigma.W = NULL,
                                                    Res.Int.W = NULL,
                                                    Res.Pi.W = NULL,
                                                    Res.Sigma.B = NULL,
                                                    Res.Int.B = NULL,
                                                    Res.Pi.B = NULL,
                                                    out = NULL, # 2l
                                                    Sinv.method = "eigen",
                                                    log2pi = FALSE,
                                                    minus.two = TRUE) {
  # map implied to 2l matrices
  if (is.null(out)) {
    out <- lav_mvreg_cluster_implied22l(
      Lp = Lp, implied = NULL,
      Res.Sigma.W = Res.Sigma.W,
      Res.Int.W = Res.Int.W, Res.Pi.W = Res.Pi.W,
      Res.Sigma.B = Res.Sigma.B,
      Res.Int.B = Res.Int.B, Res.Pi.B = Res.Pi.B
    )
  }
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz
  beta.w <- out$beta.w
  beta.b <- out$beta.b
  beta.z <- out$beta.z
  beta.wb <- out$beta.wb

  # check for beta.wb
  if (is.null(out$beta.wb)) {
    beta.wb <- rbind(beta.w, beta.b[-1, , drop = FALSE])
    beta.wb[1, ] <- beta.wb[1, , drop = FALSE] + beta.b[1, , drop = FALSE]
  }

  # log 2*pi
  LOG.2PI <- log(2 * pi)

  # Lp
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.sizes <- Lp$cluster.sizes[[2]]
  ncluster.sizes <- Lp$ncluster.sizes[[2]]
  n.s <- Lp$cluster.size.ns[[2]]

  # dependent 'y' level-2 ('Z') only variables?
  between.y.idx <- Lp$between.y.idx[[2]]

  # extract (the many) sample statistics from YLp
  sample.wb <- YLp[[2]]$sample.wb
  sample.YYres.wb1 <- YLp[[2]]$sample.YYres.wb1
  sample.XX.wb1 <- YLp[[2]]$sample.XX.wb1
  sample.wb2 <- YLp[[2]]$sample.wb2
  sample.YYres.wb2 <- YLp[[2]]$sample.YYres.wb2
  sample.YresX.wb2 <- YLp[[2]]$sample.YresX.wb2
  sample.XX.wb2 <- YLp[[2]]$sample.XX.wb2
  sample.clz.Y2.res <- YLp[[2]]$sample.clz.Y2.res
  sample.clz.Y2.XX <- YLp[[2]]$sample.clz.Y2.XX
  sample.clz.Y2.B <- YLp[[2]]$sample.clz.Y2.B
  if (length(between.y.idx) > 0L) {
    sample.clz.ZZ.res <- YLp[[2]]$sample.clz.ZZ.res
    sample.clz.ZZ.XX <- YLp[[2]]$sample.clz.ZZ.XX
    sample.clz.ZZ.B <- YLp[[2]]$sample.clz.ZZ.B
    sample.clz.YZ.res <- YLp[[2]]$sample.clz.YZ.res
    sample.clz.YZ.XX <- YLp[[2]]$sample.clz.YZ.XX
    sample.clz.YresXZ <- YLp[[2]]$sample.clz.YresXZ # zero?
    sample.clz.XWZres <- YLp[[2]]$sample.clz.XWZres
  }

  # reconstruct S.PW
  wb1.diff <- sample.wb - beta.wb
  Y1Y1.wb.res <- (sample.YYres.wb1 +
    t(wb1.diff) %*% sample.XX.wb1 %*% (wb1.diff))

  # this one is weighted -- not the same as crossprod(Y2w.res)
  wb2.diff <- sample.wb2 - beta.wb
  Y2Y2w.res <- (sample.YYres.wb2 +
    sample.YresX.wb2 %*% (wb2.diff) +
    t(wb2.diff) %*% t(sample.YresX.wb2) +
    t(wb2.diff) %*% sample.XX.wb2 %*% (wb2.diff))
  S.PW <- (Y1Y1.wb.res - Y2Y2w.res) / sum(cluster.size - 1)

  # common parts:
  sigma.w.inv <- lav_matrix_symmetric_inverse(
    S = sigma.w,
    logdet = TRUE
  )
  sigma.w.logdet <- attr(sigma.w.inv, "logdet")
  if (length(between.y.idx) > 0L) {
    sigma.zz.inv <- lav_matrix_symmetric_inverse(
      S = sigma.zz,
      logdet = TRUE
    )
    sigma.zz.logdet <- attr(sigma.zz.inv, "logdet")
    sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
    sigma.zi.zy <- t(sigma.yz.zi)
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy
  } else {
    sigma.b.z <- sigma.b
  }

  # min 2* logliklihood
  DIST <- numeric(ncluster.sizes)
  LOGDET <- numeric(ncluster.sizes)
  CONST <- numeric(ncluster.sizes)
  for (clz in seq_len(ncluster.sizes)) {
    # cluster size
    nj <- cluster.sizes[clz]

    # data between
    nj.idx <- which(cluster.size == nj)
    y2.diff <- sample.clz.Y2.B[[clz]] - beta.wb
    Y2Yc.yy <- (sample.clz.Y2.res[[clz]] +
      t(y2.diff) %*% sample.clz.Y2.XX[[clz]] %*% (y2.diff))
    if (length(between.y.idx) > 0L) {
      zz.diff <- sample.clz.ZZ.B[[clz]] - beta.z
      Y2Yc.zz <- (sample.clz.ZZ.res[[clz]] +
        t(zz.diff) %*% sample.clz.ZZ.XX[[clz]] %*% (zz.diff))
      Y2Yc.yz <- (sample.clz.YZ.res[[clz]] +
        sample.clz.YresXZ[[clz]] %*% zz.diff + # zero?
        t(y2.diff) %*% sample.clz.XWZres[[clz]] +
        t(y2.diff) %*% sample.clz.YZ.XX[[clz]] %*% zz.diff)
    }

    # construct sigma.j
    sigma.j <- (nj * sigma.b.z) + sigma.w
    sigma.j.inv <- lav_matrix_symmetric_inverse(
      S = sigma.j,
      logdet = TRUE
    )
    sigma.j.logdet <- attr(sigma.j.inv, "logdet")

    if (length(between.y.idx) > 0L) {
      sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi

      # part 1 -- zz
      Vinv.11 <- sigma.zz.inv + nj * (sigma.zi.zy %*% sigma.ji.yz.zi)
      q.zz <- sum(Vinv.11 * Y2Yc.zz)

      # part 2 -- yz
      q.yz <- -nj * sum(sigma.ji.yz.zi * Y2Yc.yz)
    } else {
      q.zz <- q.yz <- sigma.zz.logdet <- 0
    }

    # part 5 -- yyc
    q.yyc <- -nj * sum(sigma.j.inv * Y2Yc.yy)

    if (log2pi) {
      P <- nj * nrow(sigma.w) + nrow(sigma.zz)
      CONST[clz] <- P * LOG.2PI
    }
    LOGDET[clz] <- sigma.zz.logdet + sigma.j.logdet
    DIST[clz] <- q.zz + 2 * q.yz - q.yyc
  }
  # q.yya + q.yyb
  q.W <- sum(cluster.size - 1) * sum(sigma.w.inv * S.PW)
  # logdet within part
  L.W <- sum(cluster.size - 1) * sigma.w.logdet

  # -2*times logl (without the constant) (for optimization)
  loglik <- sum(LOGDET * n.s) + sum(DIST) + q.W + L.W

  if (log2pi) {
    loglik <- loglik + sum(CONST * n.s)
  }

  # functions below compute -2 * logl
  if (!minus.two) {
    loglik <- loglik / (-2)
  }

  loglik
}

# first derivative -2*logl wrt Beta.W, Beta.B, Sigma.W, Sigma.B
lav_mvreg_cluster_dlogl_2l_samplestats <- function(YLp = NULL,
                                                   Lp = NULL,
                                                   Res.Sigma.W = NULL,
                                                   Res.Int.W = NULL,
                                                   Res.Pi.W = NULL,
                                                   Res.Sigma.B = NULL,
                                                   Res.Int.B = NULL,
                                                   Res.Pi.B = NULL,
                                                   out = NULL, # 2l
                                                   return.list = FALSE,
                                                   Sinv.method = "eigen") {
  # map implied to 2l matrices
  if (is.null(out)) {
    out <- lav_mvreg_cluster_implied22l(
      Lp = Lp, implied = NULL,
      Res.Sigma.W = Res.Sigma.W,
      Res.Int.W = Res.Int.W, Res.Pi.W = Res.Pi.W,
      Res.Sigma.B = Res.Sigma.B,
      Res.Int.B = Res.Int.B, Res.Pi.B = Res.Pi.B
    )
  }
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz
  beta.w <- out$beta.w
  beta.b <- out$beta.b
  beta.z <- out$beta.z
  beta.wb <- out$beta.wb

  # check for beta.wb
  if (is.null(out$beta.wb)) {
    beta.wb <- rbind(beta.w, beta.b[-1, , drop = FALSE])
    beta.wb[1, ] <- beta.wb[1, , drop = FALSE] + beta.b[1, , drop = FALSE]
  }

  # Lp
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.sizes <- Lp$cluster.sizes[[2]]
  ncluster.sizes <- Lp$ncluster.sizes[[2]]
  n.s <- Lp$cluster.size.ns[[2]]

  within.x.idx <- Lp$within.x.idx[[1]]
  between.y.idx <- Lp$between.y.idx[[2]]

  w1.idx <- seq_len(length(within.x.idx) + 1L)
  b1.idx <- c(1L, seq_len(nrow(beta.wb))[-w1.idx])

  # extract (the many) sample statistics from YLp
  sample.wb <- YLp[[2]]$sample.wb
  sample.YYres.wb1 <- YLp[[2]]$sample.YYres.wb1
  sample.XX.wb1 <- YLp[[2]]$sample.XX.wb1
  sample.wb2 <- YLp[[2]]$sample.wb2
  sample.YYres.wb2 <- YLp[[2]]$sample.YYres.wb2
  sample.YresX.wb2 <- YLp[[2]]$sample.YresX.wb2
  sample.XX.wb2 <- YLp[[2]]$sample.XX.wb2
  sample.clz.Y2.res <- YLp[[2]]$sample.clz.Y2.res
  sample.clz.Y2.XX <- YLp[[2]]$sample.clz.Y2.XX
  sample.clz.Y2.B <- YLp[[2]]$sample.clz.Y2.B
  if (length(between.y.idx) > 0L) {
    sample.clz.ZZ.res <- YLp[[2]]$sample.clz.ZZ.res
    sample.clz.ZZ.XX <- YLp[[2]]$sample.clz.ZZ.XX
    sample.clz.ZZ.B <- YLp[[2]]$sample.clz.ZZ.B
    sample.clz.YZ.res <- YLp[[2]]$sample.clz.YZ.res
    sample.clz.YZ.XX <- YLp[[2]]$sample.clz.YZ.XX
    sample.clz.YresXZ <- YLp[[2]]$sample.clz.YresXZ # zero?
    sample.clz.XWZres <- YLp[[2]]$sample.clz.XWZres
  }

  # reconstruct S.PW
  wb1.diff <- sample.wb - beta.wb
  Y1Y1.wb.res <- (sample.YYres.wb1 +
    t(wb1.diff) %*% sample.XX.wb1 %*% (wb1.diff))

  # this one is weighted -- not the same as crossprod(Y2w.res)
  wb2.diff <- sample.wb2 - beta.wb
  Y2Y2w.res <- (sample.YYres.wb2 +
    sample.YresX.wb2 %*% (wb2.diff) +
    t(wb2.diff) %*% t(sample.YresX.wb2) +
    t(wb2.diff) %*% sample.XX.wb2 %*% (wb2.diff))
  S.PW <- (Y1Y1.wb.res - Y2Y2w.res) / sum(cluster.size - 1)

  # common parts:
  sigma.w.inv <- lav_matrix_symmetric_inverse(S = sigma.w)

  G.beta.w <- matrix(0, ncluster.sizes, length(beta.w))
  G.beta.b <- matrix(0, ncluster.sizes, length(beta.b))
  G.beta.wb <- matrix(0, ncluster.sizes, length(beta.wb))
  G.sigma.w1 <- matrix(0, ncluster.sizes, length(lav_matrix_vech(sigma.w)))
  G.sigma.b <- matrix(0, ncluster.sizes, length(lav_matrix_vech(sigma.b)))

  if (length(between.y.idx) > 0L) {
    G.beta.z <- matrix(0, ncluster.sizes, length(beta.z))
    G.sigma.zz <- matrix(0, ncluster.sizes, length(lav_matrix_vech(sigma.zz)))
    G.sigma.yz <- matrix(0, ncluster.sizes, length(sigma.yz))

    sigma.zz.inv <- lav_matrix_symmetric_inverse(S = sigma.zz)
    sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
    sigma.zi.zy <- t(sigma.yz.zi)
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy

    for (clz in seq_len(ncluster.sizes)) {
      # cluster size
      nj <- cluster.sizes[clz]

      y2.diff <- sample.clz.Y2.B[[clz]] - beta.wb
      XX.y2.diff <- sample.clz.Y2.XX[[clz]] %*% y2.diff
      Y2Yc.yy <- sample.clz.Y2.res[[clz]] + crossprod(y2.diff, XX.y2.diff)

      zz.diff <- sample.clz.ZZ.B[[clz]] - beta.z
      Y2Yc.zz <- (sample.clz.ZZ.res[[clz]] +
        t(zz.diff) %*% sample.clz.ZZ.XX[[clz]] %*% (zz.diff))
      Y2Yc.yz <- (sample.clz.YZ.res[[clz]] +
        sample.clz.YresXZ[[clz]] %*% zz.diff + # zero?
        t(y2.diff) %*% sample.clz.XWZres[[clz]] +
        t(y2.diff) %*% sample.clz.YZ.XX[[clz]] %*% zz.diff)

      # construct sigma.j
      sigma.j <- (nj * sigma.b.z) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j)
      sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
      sigma.zi.zy.ji <- t(sigma.ji.yz.zi)
      sigma.ji.yz <- sigma.j.inv %*% sigma.yz
      ns.sigma.j.inv <- n.s[clz] * sigma.j.inv
      ns.sigma.zz.inv <- n.s[clz] * sigma.zz.inv
      ns.sigma.yz <- n.s[clz] * sigma.yz
      ns.sigma.ji.yz.zi <- n.s[clz] * sigma.ji.yz.zi

      # common parts
      ZZ.zi.yz.ji <- Y2Yc.zz %*% sigma.zi.zy.ji
      ji.YZ.zi <- sigma.j.inv %*% Y2Yc.yz %*% sigma.zz.inv

      jYZj.yy <- sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv
      jYZj.yz <- tcrossprod(ji.YZ.zi, sigma.ji.yz)
      jYZj.zz <- sigma.ji.yz.zi %*% ZZ.zi.yz.ji

      jYZj <- nj * (jYZj.yy + jYZj.zz - jYZj.yz - t(jYZj.yz))

      # SIGMA.W (between part)
      g.sigma.w1 <- ns.sigma.j.inv - jYZj
      tmp <- g.sigma.w1 * 2
      diag(tmp) <- diag(g.sigma.w1)
      G.sigma.w1[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g.sigma.b <- nj * g.sigma.w1
      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.sigma.b[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZZ
      YZ1 <- ZZ.zi.yz.ji %*% sigma.yz
      YZ2 <- crossprod(Y2Yc.yz, sigma.ji.yz)
      tmp <- (t(sigma.yz) %*% g.sigma.w1 %*% sigma.yz
        - 1 / nj * Y2Yc.zz - t(YZ1) - YZ1 + t(YZ2) + YZ2)
      g.sigma.zz <- (ns.sigma.zz.inv +
        nj * sigma.zz.inv %*% tmp %*% sigma.zz.inv)
      tmp <- g.sigma.zz * 2
      diag(tmp) <- diag(g.sigma.zz)
      G.sigma.zz[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZY
      tmp1 <- crossprod(ZZ.zi.yz.ji, sigma.zz.inv)
      tmp2 <- ns.sigma.ji.yz.zi
      tmp3 <- ji.YZ.zi
      tmp4 <- jYZj %*% sigma.yz.zi
      g.sigma.yz <- 2 * nj * (tmp1 - tmp2 - tmp3 + tmp4)
      G.sigma.yz[clz, ] <- lav_matrix_vec(g.sigma.yz)

      # BETA.Z
      A <- (sigma.zz.inv + nj * (sigma.zi.zy.ji %*% sigma.yz.zi)) # symm!
      B <- nj * (sigma.zi.zy.ji)
      tmp.z <- (sample.clz.ZZ.XX[[clz]] %*% zz.diff %*% A -
        (t(sample.clz.YresXZ[[clz]]) +
          t(sample.clz.YZ.XX[[clz]]) %*% y2.diff) %*% t(B))
      G.beta.z[clz, ] <- as.vector(-2 * tmp.z)

      # BETA.W (between part only) + BETA.B
      tmp <- (sample.clz.XWZres[[clz]] +
        sample.clz.YZ.XX[[clz]] %*% zz.diff)
      out.b <- tmp %*% sigma.zi.zy.ji - XX.y2.diff %*% sigma.j.inv
      out.w <- out.b + XX.y2.diff %*% sigma.w.inv
      tmp.b <- out.b[b1.idx, , drop = FALSE]
      tmp.w <- out.w[w1.idx, , drop = FALSE]
      G.beta.b[clz, ] <- as.vector(2 * nj * tmp.b)
      G.beta.w[clz, ] <- as.vector(2 * nj * tmp.w)
    } # clz

    d.beta.w1 <- matrix(colSums(G.beta.w), nrow(beta.w), ncol(beta.w))
    d.beta.b <- matrix(colSums(G.beta.b), nrow(beta.b), ncol(beta.b))
    d.sigma.w1 <- lav_matrix_vech_reverse(colSums(G.sigma.w1))
    d.sigma.b <- lav_matrix_vech_reverse(colSums(G.sigma.b))

    # z
    d.beta.z <- matrix(colSums(G.beta.z), nrow(beta.z), ncol(beta.z))
    d.sigma.zz <- lav_matrix_vech_reverse(colSums(G.sigma.zz))
    d.sigma.yz <- matrix(colSums(G.sigma.yz), nrow(sigma.yz), ncol(sigma.yz))
  } # between.y.idx

  else { # no beween.y.idx

    for (clz in seq_len(ncluster.sizes)) {
      # cluster size
      nj <- cluster.sizes[clz]

      y2.diff <- sample.clz.Y2.B[[clz]] - beta.wb
      XX.y2.diff <- sample.clz.Y2.XX[[clz]] %*% y2.diff
      Y2Yc.yy <- sample.clz.Y2.res[[clz]] + crossprod(y2.diff, XX.y2.diff)

      # construct sigma.j
      sigma.j <- (nj * sigma.b) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j)

      # common part
      jYYj <- nj * sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv

      # SIGMA.W (between part)
      g.sigma.w1 <- (n.s[clz] * sigma.j.inv) - jYYj
      tmp <- g.sigma.w1 * 2
      diag(tmp) <- diag(g.sigma.w1)
      G.sigma.w1[clz, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g.sigma.b <- nj * g.sigma.w1
      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.sigma.b[clz, ] <- lav_matrix_vech(tmp)

      # BETA.W (between part only) + BETA.B
      out.b <- -1 * XX.y2.diff %*% sigma.j.inv
      out.w <- out.b + XX.y2.diff %*% sigma.w.inv
      tmp.b <- out.b[b1.idx, , drop = FALSE]
      tmp.w <- out.w[w1.idx, , drop = FALSE]
      G.beta.b[clz, ] <- as.vector(2 * nj * tmp.b)
      G.beta.w[clz, ] <- as.vector(2 * nj * tmp.w)
    } # cl

    d.beta.w1 <- matrix(colSums(G.beta.w), nrow(beta.w), ncol(beta.w))
    d.beta.b <- matrix(colSums(G.beta.b), nrow(beta.b), ncol(beta.b))
    d.sigma.w1 <- lav_matrix_vech_reverse(colSums(G.sigma.w1))
    d.sigma.b <- lav_matrix_vech_reverse(colSums(G.sigma.b))

    # z
    d.beta.z <- matrix(0, 0L, 0L)
    d.sigma.zz <- matrix(0, 0L, 0L)
    d.sigma.yz <- matrix(0, 0L, 0L)
  } # no-between-y

  # Sigma.W (bis)
  d.sigma.w2 <- sum(cluster.size - 1) * (sigma.w.inv
  - sigma.w.inv %*% S.PW %*% sigma.w.inv)
  tmp <- d.sigma.w2 * 2
  diag(tmp) <- diag(d.sigma.w2)
  d.sigma.w2 <- tmp

  d.sigma.w <- d.sigma.w1 + d.sigma.w2

  # beta.w (bis)
  d.beta.w2 <- -2 * (sample.XX.wb1 %*% (sample.wb - beta.wb))[w1.idx, , drop = FALSE] %*% sigma.w.inv

  d.beta.w <- d.beta.w1 + d.beta.w2

  # rearrange
  dimplied <- lav_mvreg_cluster_2l2implied(Lp,
    sigma.w = d.sigma.w, sigma.b = d.sigma.b,
    sigma.zz = d.sigma.zz, sigma.yz = d.sigma.yz,
    beta.w = d.beta.w, beta.b = d.beta.b, beta.z = d.beta.z
  )

  if (return.list) {
    return(dimplied)
  }

  # as a single vector
  out <- c(
    drop(dimplied$res.int[[1]]),
    lav_matrix_vec(dimplied$res.slopes[[1]]),
    lav_matrix_vech(dimplied$res.cov[[1]]),
    drop(dimplied$res.int[[2]]),
    lav_matrix_vec(dimplied$res.slopes[[2]]),
    lav_matrix_vech(dimplied$res.cov[[2]])
  )
  out
}

# cluster-wise scores -2*logl wrt Beta.W, Beta.B, Sigma.W, Sigma.B
lav_mvreg_cluster_scores_2l <- function(Y1 = NULL,
                                        YLp = NULL,
                                        Lp = NULL,
                                        Res.Sigma.W = NULL,
                                        Res.Int.W = NULL,
                                        Res.Pi.W = NULL,
                                        Res.Sigma.B = NULL,
                                        Res.Int.B = NULL,
                                        Res.Pi.B = NULL,
                                        out = NULL, # 2l
                                        Sinv.method = "eigen") {
  # map implied to 2l matrices
  if (is.null(out)) {
    out <- lav_mvreg_cluster_implied22l(
      Lp = Lp, implied = NULL,
      Res.Sigma.W = Res.Sigma.W,
      Res.Int.W = Res.Int.W, Res.Pi.W = Res.Pi.W,
      Res.Sigma.B = Res.Sigma.B,
      Res.Int.B = Res.Int.B, Res.Pi.B = Res.Pi.B
    )
  }
  sigma.w <- out$sigma.w
  sigma.b <- out$sigma.b
  sigma.zz <- out$sigma.zz
  sigma.yz <- out$sigma.yz
  beta.w <- out$beta.w
  beta.b <- out$beta.b
  beta.z <- out$beta.z
  beta.wb <- out$beta.wb

  # check for beta.wb
  if (is.null(out$beta.wb)) {
    beta.wb <- rbind(beta.w, beta.b[-1, , drop = FALSE])
    beta.wb[1, ] <- beta.wb[1, , drop = FALSE] + beta.b[1, , drop = FALSE]
  }

  # Lp
  nclusters <- Lp$nclusters[[2]]
  cluster.size <- Lp$cluster.size[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]

  within.x.idx <- Lp$within.x.idx[[1]]
  between.idx <- Lp$between.idx[[2]]
  between.y.idx <- Lp$between.y.idx[[2]]
  between.x.idx <- Lp$between.x.idx[[2]]

  y1.idx <- Lp$ov.y.idx[[1]]
  x1.idx <- c(within.x.idx, between.x.idx) # in that order

  # residuals for 'Y'
  Y1.wb <- Y1[, y1.idx, drop = FALSE]

  if (length(x1.idx) > 0L) {
    EXO.wb <- cbind(1, Y1[, x1.idx, drop = FALSE])
    Y1.wb.hat <- EXO.wb %*% beta.wb
    Y1.wb.res <- Y1.wb - Y1.wb.hat
  } else {
    Y1.wb.res <- Y1.wb
  }

  # residuals 'Y' (level 2)
  Y2 <- YLp[[2]]$Y2
  if (length(x1.idx) > 0L) {
    EXO.wb2 <- cbind(1, Y2[, x1.idx, drop = FALSE])
    Y2w.res <- Y2[, y1.idx, drop = FALSE] - EXO.wb2 %*% beta.wb
  } else {
    EXO.wb2 <- matrix(1, nrow(Y2), 1L)
    Y2w.res <- Y2[, y1.idx, drop = FALSE]
  }

  # residual 'Z' (level 2)
  if (length(between.y.idx) > 0L) {
    if (length(between.x.idx) > 0L) {
      EXO.z <- cbind(1, Y2[, between.x.idx, drop = FALSE])
      Y2.z <- Y2[, between.y.idx, drop = FALSE]
      Y2z.res <- Y2.z - EXO.z %*% beta.z
      # sample.z
      # XX.z <- crossprod(EXO.z)
      # sample.z <- try(solve(XX.z, crossprod(EXO.z, Y2.z)))
      # if(inherits(sample.z, "try-error")) {
      #    sample.z <- MASS::ginv(XX.z) %*% crossprod(EXO.z, Y2.z)
      # }

      # sample.wb2
      # sample.wb2 <- YLp[[2]]$sample.wb2
    } else {
      Y2z.res <- Y2[, between.y.idx, drop = FALSE]
    }
  }

  # common parts:
  sigma.w.inv <- lav_matrix_symmetric_inverse(S = sigma.w)

  G.beta.w1 <- matrix(0, nclusters, length(beta.w))
  G.beta.b <- matrix(0, nclusters, length(beta.b))
  G.beta.wb <- matrix(0, nclusters, length(beta.wb))
  G.sigma.w1 <- matrix(0, nclusters, length(lav_matrix_vech(sigma.w)))
  G.sigma.b <- matrix(0, nclusters, length(lav_matrix_vech(sigma.b)))


  if (length(between.y.idx) > 0L) {
    G.beta.z <- matrix(0, nclusters, length(beta.z))
    G.sigma.zz <- matrix(0, nclusters, length(lav_matrix_vech(sigma.zz)))
    G.sigma.yz <- matrix(0, nclusters, length(sigma.yz))

    sigma.zz.inv <- lav_matrix_symmetric_inverse(S = sigma.zz)
    sigma.yz.zi <- sigma.yz %*% sigma.zz.inv
    sigma.zi.zy <- t(sigma.yz.zi)
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy

    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster.size[cl]

      # data within for the cluster (centered)
      Y1m <- Y1.wb.res[cluster.idx == cl, , drop = FALSE]
      yc <- Y2w.res[cl, ]

      # data between
      zc <- Y2z.res[cl, ]
      Y2Yc.yy <- tcrossprod(Y2w.res[cl, ])
      Y2Yc.zz <- tcrossprod(Y2z.res[cl, ])
      Y2Yc.yz <- tcrossprod(Y2w.res[cl, ], Y2z.res[cl, ])

      # construct sigma.j
      sigma.j <- (nj * sigma.b.z) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j)
      sigma.ji.yz.zi <- sigma.j.inv %*% sigma.yz.zi
      sigma.zi.zy.ji <- t(sigma.ji.yz.zi)
      sigma.ji.yz <- sigma.j.inv %*% sigma.yz

      # common parts
      ZZ.zi.yz.ji <- Y2Yc.zz %*% sigma.zi.zy.ji
      ji.YZ.zi <- sigma.j.inv %*% Y2Yc.yz %*% sigma.zz.inv

      jYZj.yy <- sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv
      jYZj.yz <- tcrossprod(ji.YZ.zi, sigma.ji.yz)
      jYZj.zz <- sigma.ji.yz.zi %*% ZZ.zi.yz.ji

      jYZj <- nj * (jYZj.yy + jYZj.zz - jYZj.yz - t(jYZj.yz))

      # SIGMA.W (between part)
      g.sigma.w1 <- sigma.j.inv - jYZj
      tmp <- g.sigma.w1 * 2
      diag(tmp) <- diag(g.sigma.w1)
      G.sigma.w1[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.W (within part)
      # g.sigma.w2 <- ( (nj-1) * sigma.w.inv
      #    - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv )
      # tmp <- g.sigma.w2*2; diag(tmp) <- diag(g.sigma.w2)
      # G.sigma.w2[cl,] <- lav_matrix_vech(tmp)
      # G.sigma.w[cl,] <- G.sigma.w1[cl,] + G.sigma.w2[cl,]

      # SIGMA.B
      g.sigma.b <- nj * g.sigma.w1
      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.sigma.b[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZZ
      YZ1 <- ZZ.zi.yz.ji %*% sigma.yz
      YZ2 <- crossprod(Y2Yc.yz, sigma.ji.yz)
      tmp <- (t(sigma.yz) %*% g.sigma.w1 %*% sigma.yz
        - (1 / nj * Y2Yc.zz + t(YZ1) + YZ1 - t(YZ2) - YZ2))
      g.sigma.zz <- (sigma.zz.inv +
        nj * sigma.zz.inv %*% tmp %*% sigma.zz.inv)
      tmp <- g.sigma.zz * 2
      diag(tmp) <- diag(g.sigma.zz)
      G.sigma.zz[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZY
      # g.sigma.yz <- 2 * nj * (
      #              (sigma.j.inv %*%
      #                  (sigma.yz.zi %*% Y2Yc.zz - sigma.yz - Y2Yc.yz)
      #                   + jYZj %*% sigma.yz) %*% sigma.zz.inv )
      tmp1 <- crossprod(ZZ.zi.yz.ji, sigma.zz.inv)
      tmp2 <- sigma.ji.yz.zi
      tmp3 <- ji.YZ.zi
      tmp4 <- jYZj %*% sigma.yz.zi
      g.sigma.yz <- 2 * nj * (tmp1 - tmp2 - tmp3 + tmp4)
      G.sigma.yz[cl, ] <- lav_matrix_vec(g.sigma.yz)

      # BETA.Z
      # here, we avoid the (sample.z - beta.z) approach
      exo.z <- cbind(1, Y2[cl, between.x.idx, drop = FALSE])
      tmp1 <- (sigma.zz.inv + nj * (sigma.zi.zy.ji %*% sigma.yz.zi)) %*% zc
      tmp2 <- nj * (sigma.zi.zy.ji) %*% yc
      tmp.z <- crossprod(exo.z, drop(tmp1 - tmp2))
      G.beta.z[cl, ] <- as.vector(-2 * tmp.z)

      # BETA.W
      #            exo.w <- cbind(1,
      #                           Y1[cluster.idx == cl, within.x.idx, drop = FALSE])
      #            G.beta.w[cl,] <- as.vector( 2 * t(exo.w) %*% (
      #                                 matrix(1, nj, 1) %x% (zc %*% sigma.zi.zy.ji -
      #                                                       yc %*% sigma.j.inv +
      #                                                       yc %*% sigma.w.inv) -
      #                                 Y1m %*% sigma.w.inv) )

      # BETA.W (between part only)
      exo2.w <- cbind(1, Y2[cl, within.x.idx, drop = FALSE])
      tmp2 <- (zc %*% sigma.zi.zy.ji -
        yc %*% sigma.j.inv +
        yc %*% sigma.w.inv)
      G.beta.w1[cl, ] <- as.vector(2 * nj * crossprod(exo2.w, tmp2))

      # BETA.W (within part only)
      # exo.w <- cbind(1,
      #               Y1[cluster.idx == cl, within.x.idx, drop = FALSE])
      # tmp1 <- - Y1m %*% sigma.w.inv
      # G.beta.ww <- as.vector( 2 * crossprod(exo.w, tmp1) )
      # G.beta.w[cl,] <- G.beta.w1 + G.beta.ww
      # G.beta.w2[cl,] <- G.beta.ww

      # BETA.B
      exo2.b <- cbind(1, Y2[cl, between.x.idx, drop = FALSE])
      tmp <- (zc %*% sigma.zi.zy.ji - yc %*% sigma.j.inv)
      G.beta.b[cl, ] <- as.vector(2 * nj * crossprod(exo2.b, tmp))
    } # cl
  } # between.y.idx

  else { # no beween.y.idx

    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster.size[cl]

      # data within for the cluster (centered)
      Y1m <- Y1.wb.res[cluster.idx == cl, , drop = FALSE]
      yc <- Y2w.res[cl, ]

      # data between
      Y2Yc.yy <- tcrossprod(Y2w.res[cl, ])

      # construct sigma.j
      sigma.j <- (nj * sigma.b) + sigma.w
      sigma.j.inv <- lav_matrix_symmetric_inverse(S = sigma.j)

      # common part
      jYYj <- nj * sigma.j.inv %*% Y2Yc.yy %*% sigma.j.inv

      # SIGMA.W
      # g.sigma.w <- ( (nj-1) * sigma.w.inv
      #    - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv
      #    + sigma.j.inv - jYYj )
      # tmp <- g.sigma.w*2; diag(tmp) <- diag(g.sigma.w)
      # G.sigma.w[cl,] <- lav_matrix_vech(tmp)

      # SIGMA.W (between part)
      g.sigma.w1 <- sigma.j.inv - jYYj
      tmp <- g.sigma.w1 * 2
      diag(tmp) <- diag(g.sigma.w1)
      G.sigma.w1[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g.sigma.b <- nj * (sigma.j.inv - jYYj)
      tmp <- g.sigma.b * 2
      diag(tmp) <- diag(g.sigma.b)
      G.sigma.b[cl, ] <- lav_matrix_vech(tmp)

      # BETA.W (between part only)
      exo2.w <- cbind(1, Y2[cl, within.x.idx, drop = FALSE])
      tmp2 <- (-yc %*% sigma.j.inv + yc %*% sigma.w.inv)
      G.beta.w1[cl, ] <- as.vector(2 * nj * crossprod(exo2.w, tmp2))

      # BETA.B
      exo2.b <- cbind(1, Y2[cl, between.x.idx, drop = FALSE])
      tmp <- -yc %*% sigma.j.inv
      G.beta.b[cl, ] <- as.vector(2 * nj * crossprod(exo2.b, tmp))
    } # cl
  } # no-between-y

  # beta.w (bis)

  #    d.beta.w2 <- -2 * t(EXO.wb[,1:(length(within.x.idx) + 1L), drop = FALSE]) %*% Y1.wb.res %*% sigma.w.inv

  Y1.wb.res.i <- Y1.wb.res %*% sigma.w.inv
  w1.idx <- seq_len(length(within.x.idx) + 1L)
  a1.idx <- rep(w1.idx, times = ncol(Y1.wb.res.i))
  b1.idx <- rep(seq_len(ncol(Y1.wb.res.i)), each = length(w1.idx))
  TMP <- EXO.wb[, a1.idx, drop = FALSE] * Y1.wb.res.i[, b1.idx, drop = FALSE]
  G.beta.w2 <- -2 * rowsum.default(TMP, cluster.idx,
    reorder = FALSE,
    na.rm = TRUE
  )
  G.beta.w <- G.beta.w1 + G.beta.w2

  # Sigma.W (bis)
  # d.sigma.w2 <- sum(cluster.size - 1) * ( sigma.w.inv
  #               - sigma.w.inv %*% S.PW %*% sigma.w.inv )
  # tmp <- d.sigma.w2*2; diag(tmp) <- diag(d.sigma.w2)
  # d.sigma.w2 <- tmp

  # g.sigma.w2 <- ( (nj-1) * sigma.w.inv
  #        - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv )

  Y1a.res <- Y1.wb.res - Y2w.res[cluster.idx, , drop = FALSE]
  Y1a.res.i <- Y1a.res %*% sigma.w.inv
  idx1 <- lav_matrix_vech_col_idx(nrow(sigma.w))
  idx2 <- lav_matrix_vech_row_idx(nrow(sigma.w))
  SW2 <- matrix(lav_matrix_vech(sigma.w.inv),
    nrow = nclusters,
    length(lav_matrix_vech(sigma.w.inv)), byrow = TRUE
  )
  SW2 <- SW2 * (cluster.size - 1)
  TMP <- Y1a.res.i[, idx1, drop = FALSE] * Y1a.res.i[, idx2, drop = FALSE]
  TMP2 <- rowsum.default(TMP, cluster.idx, reorder = FALSE, na.rm = TRUE)
  G.sigma.w2 <- 2 * (SW2 - TMP2)
  diagh.idx <- lav_matrix_diagh_idx(nrow(sigma.w))
  G.sigma.w2[, diagh.idx] <- G.sigma.w2[, diagh.idx, drop = FALSE] / 2
  G.sigma.w <- G.sigma.w1 + G.sigma.w2



  # rearrange columns to Res.Int.W, Res.Pi.W, Res.Sigma.W,
  #                      Res.Int.B, Res.Pi.B, Res.Sigma.B

  # ov.idx per level
  ov.idx <- Lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p.tilde <- length(unique(c(ov.idx[[1]], ov.idx[[2]])))
  p.tilde.star <- p.tilde * (p.tilde + 1) / 2
  B.tilde <- lav_matrix_vech_reverse(seq_len(p.tilde.star))

  # only 'y'
  ov.y.idx <- Lp$ov.y.idx

  # two levels only (for now)
  ov.y.idx1 <- ov.y.idx[[1]]
  ov.y.idx2 <- ov.y.idx[[2]]

  # WITHIN (is easy)
  BETA.W.idx <- matrix(seq_len(length(beta.w)), nrow(beta.w), ncol(beta.w))
  BETA.B.idx <- matrix(seq_len(length(beta.b)), nrow(beta.b), ncol(beta.b))

  Res.Int.W <- G.beta.w[, BETA.W.idx[1L, ], drop = FALSE]
  Res.Pi.W <- G.beta.w[, lav_matrix_vecr(BETA.W.idx[-1L, ]), drop = FALSE]
  Res.Sigma.W <- G.sigma.w

  # Sigma.B
  Sigma.B.tilde <- matrix(0, nclusters, p.tilde.star)
  col.idx <- lav_matrix_vech(B.tilde[ov.y.idx1, ov.y.idx1, drop = FALSE])
  Sigma.B.tilde[, col.idx] <- G.sigma.b

  # Int.B
  BETA.B.tilde <- matrix(seq_len(nrow(beta.b) * p.tilde), nrow(beta.b), p.tilde)
  Int.B <- matrix(0, nclusters, p.tilde)
  Int.B[, ov.y.idx1] <- G.beta.b[, BETA.B.idx[1L, ]]

  # Pi.B
  Pi.B <- matrix(0, nclusters, p.tilde * (nrow(beta.b) - 1L))
  col.idx <- lav_matrix_vecr(BETA.B.tilde[-1L, ov.y.idx1, drop = FALSE])
  Pi.B[, col.idx] <- G.beta.b[, lav_matrix_vecr(BETA.B.idx[-1L, ]), drop = FALSE]

  if (length(between.y.idx) > 0L) {
    # Sigma.B: add yz/zz parts
    col.idx <- lav_matrix_vec(B.tilde[ov.y.idx1, between.y.idx, drop = FALSE])
    Sigma.B.tilde[, col.idx] <- G.sigma.yz
    col.idx <- lav_matrix_vech(B.tilde[between.y.idx, between.y.idx,
      drop = FALSE
    ])
    Sigma.B.tilde[, col.idx] <- G.sigma.zz

    # Int.B: add z-part
    BETA.Z.idx <- matrix(seq_len(length(beta.z)), nrow(beta.z), ncol(beta.z))
    Int.B[, between.y.idx] <- G.beta.z[, BETA.Z.idx[1L, ], drop = FALSE]

    # Pi.B: add beta.z
    col.idx <- lav_matrix_vecr(BETA.B.tilde[-1L, between.y.idx, drop = FALSE])
    Pi.B[, col.idx] <-
      G.beta.z[, lav_matrix_vecr(BETA.Z.idx[-1L, ]), drop = FALSE]
  }

  # only extract ov.y.idx2 for BETWEEN
  col.idx <- lav_matrix_vech(B.tilde[ov.y.idx2, ov.y.idx2, drop = FALSE])
  Res.Sigma.B <- Sigma.B.tilde[, col.idx, drop = FALSE]

  Res.Int.B <- Int.B[, ov.y.idx2, drop = FALSE]

  col.idx <- lav_matrix_vecr(BETA.B.tilde[-1, ov.y.idx2])
  Res.Pi.B <- Pi.B[, col.idx, drop = FALSE]

  SCORES <- cbind(
    Res.Int.W, Res.Pi.W, Res.Sigma.W,
    Res.Int.B, Res.Pi.B, Res.Sigma.B
  )

  SCORES
}

# first-order information: outer crossprod of scores per cluster
lav_mvreg_cluster_information_firstorder <- function(Y1 = NULL,
                                                     YLp = NULL,
                                                     Lp = NULL,
                                                     Res.Sigma.W = NULL,
                                                     Res.Int.W = NULL,
                                                     Res.Pi.W = NULL,
                                                     Res.Sigma.B = NULL,
                                                     Res.Int.B = NULL,
                                                     Res.Pi.B = NULL,
                                                     divide.by.two = FALSE,
                                                     Sinv.method = "eigen") {
  N <- NROW(Y1)

  SCORES <- lav_mvreg_cluster_scores_2l(
    Y1 = Y1,
    YLp = YLp,
    Lp = Lp,
    Res.Sigma.W = Res.Sigma.W,
    Res.Int.W = Res.Int.W,
    Res.Pi.W = Res.Pi.W,
    Res.Sigma.B = Res.Sigma.B,
    Res.Int.B = Res.Int.B,
    Res.Pi.B = Res.Pi.B,
    Sinv.method = Sinv.method
  )

  # divide by 2 (if we want scores wrt objective function)
  if (divide.by.two) {
    SCORES <- SCORES / 2
  }

  # unit information
  information <- crossprod(SCORES) / Lp$nclusters[[2]]

  information
}
