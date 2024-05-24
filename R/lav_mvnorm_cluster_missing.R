# loglikelihood clustered/twolevel data in the presence of missing data

# YR:
# - objective function: first version around March 2021 (see Psych paper)
# - analytic gradient: first version around May 2021


# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
lav_mvnorm_cluster_missing_loglik_samplestats_2l <- function(Y1 = NULL,
                                                             Y2 = NULL,
                                                             Lp = NULL,
                                                             Mp = NULL,
                                                             Mu.W = NULL,
                                                             Sigma.W = NULL,
                                                             Mu.B = NULL,
                                                             Sigma.B = NULL,
                                                             Sinv.method = "eigen",
                                                             log2pi = FALSE,
                                                             loglik.x = 0,
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
  between.idx <- Lp$between.idx[[2]]
  both.idx <- Lp$both.idx[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]

  # sanity checks
  if (any(diag(sigma.w) < 0) || any(diag(sigma.b) < 0)) {
    return(+Inf)
  }

  # check is both.idx part of sigma.b is 'too' negative; if so, return +Inf
  ev <- eigen(sigma.b[both.idx, both.idx, drop = FALSE],
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (any(ev < -0.05)) {
    return(+Inf)
  }

  # cat("sigma.w = \n"); print(sigma.w)
  # cat("sigma.b = \n"); print(sigma.b)
  # cat("mu.y = \n"); print(mu.y)

  # global
  sigma.w.inv <- solve.default(sigma.w)
  sigma.w.logdet <- log(det(sigma.w))
  sigma.b <- sigma.b[both.idx, both.idx] # only both part

  # y
  ny <- ncol(sigma.w)
  if (length(between.idx) > 0L) {
    Y1w <- Y1[, -between.idx, drop = FALSE]
  } else {
    Y1w <- Y1
  }
  Y1w.c <- t(t(Y1w) - mu.y)
  PIJ <- matrix(0, nrow(Y1w.c), ny)

  # z
  nz <- length(between.idx)
  if (nz > 0L) {
    # check is sigma.zz is PD; if not, return +Inf
    ev <- eigen(sigma.zz, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < sqrt(.Machine$double.eps))) {
      return(+Inf)
    }

    Z <- Y2[, between.idx, drop = FALSE]
    Z.c <- t(t(Z) - mu.z)

    sigma.yz <- sigma.yz[both.idx, , drop = FALSE] # only both part
    sigma.zy <- t(sigma.yz)
    sigma.zz.inv <- solve.default(sigma.zz)
    sigma.zz.logdet <- log(det(sigma.zz))
    sigma.zi.zy <- sigma.zz.inv %*% sigma.zy
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy
    GZ <- Z.c %*% sigma.zz.inv # for complete cases only
  }

  # containters per cluster
  q.yy.b <- q.zy <- q.zz.b <- numeric(nclusters)
  IBZA.j.logdet <- numeric(nclusters)
  ALIST <- rep(list(matrix(
    0, length(both.idx),
    length(both.idx)
  )), nclusters)

  # Z per missing pattern
  if (nz > 0L) {
    Zp <- Mp$Zp
    ZPAT2J <- integer(nclusters) # which sigma.b.z per cluster
    SIGMA.B.Z <- vector("list", length = Zp$npatterns + 1L)

    sigma.j.zz.logdet <- q.zz.a <- 0
    for (p in seq_len(Zp$npatterns)) {
      freq <- Zp$freq[p]
      z.na.idx <- which(!Zp$pat[p, ])
      j.idx <- Zp$case.idx[[p]] # cluster indices with this pattern
      ZPAT2J[j.idx] <- p

      if (length(z.na.idx) > 0L) {
        zp <- sigma.zz[-z.na.idx, -z.na.idx, drop = FALSE]
        zp.inv <- lav_matrix_symmetric_inverse_update(
          S.inv = sigma.zz.inv, rm.idx = z.na.idx,
          logdet = TRUE, S.logdet = sigma.zz.logdet
        )
        zp.logdet <- attr(zp.inv, "logdet")
        sigma.j.zz.logdet <- sigma.j.zz.logdet + (zp.logdet * freq)

        GZ[j.idx, -z.na.idx] <- Z.c[j.idx, -z.na.idx] %*% zp.inv

        yziy <- (sigma.yz[, -z.na.idx, drop = FALSE] %*% zp.inv %*%
          sigma.zy[-z.na.idx, , drop = FALSE])
        SIGMA.B.Z[[p]] <- (sigma.b - yziy)
      } else {
        # complete case
        sigma.j.zz.logdet <-
          sigma.j.zz.logdet + (sigma.zz.logdet * freq)
        SIGMA.B.Z[[p]] <- sigma.b.z
      }
    } # p

    # add empty patterns (if any)
    if (length(Zp$empty.idx) > 0L) {
      ZPAT2J[Zp$empty.idx] <- p + 1L
      SIGMA.B.Z[[p + 1L]] <- sigma.b
    }

    q.zz.a <- sum(GZ * Z.c, na.rm = TRUE)
    GZ0 <- GZ
    GZ0[is.na(GZ0)] <- 0
    GJ <- GZ0 %*% sigma.zy # only both part
  }

  # Y per missing pattern
  W.logdet <- 0
  MPi <- integer(nrow(Y1))
  for (p in seq_len(Mp$npatterns)) {
    freq <- Mp$freq[p]
    na.idx <- which(!Mp$pat[p, ])
    j.idx <- Mp$j.idx[[p]]
    j1.idx <- Mp$j1.idx[[p]]
    TAB <- integer(nclusters)
    TAB[j1.idx] <- Mp$j.freq[[p]]

    # compute sigma.w.inv for this pattern
    if (length(na.idx) > 0L) {
      MPi[Mp$case.idx[[p]]] <- p
      wp <- sigma.w[-na.idx, -na.idx, drop = FALSE]
      wp.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = sigma.w.inv, rm.idx = na.idx,
        logdet = TRUE, S.logdet = sigma.w.logdet
      )
      wp.logdet <- attr(wp.inv, "logdet")
      W.logdet <- W.logdet + (wp.logdet * freq)

      PIJ[Mp$case.idx[[p]], -na.idx] <-
        Y1w.c[Mp$case.idx[[p]], -na.idx] %*% wp.inv

      A.j <- matrix(0, ny, ny)
      A.j[-na.idx, -na.idx] <- wp.inv
      for (j in j1.idx) {
        ALIST[[j]] <- ALIST[[j]] + (A.j[both.idx, both.idx] * TAB[j])
      }
      # WIP[[p]][-na.idx, -na.idx] <- wp.inv
    } else {
      # complete case
      W.logdet <- W.logdet + (sigma.w.logdet * freq)
      PIJ[Mp$case.idx[[p]], ] <-
        Y1w.c[Mp$case.idx[[p]], ] %*% sigma.w.inv
      for (j in j1.idx) {
        ALIST[[j]] <-
          ALIST[[j]] + (sigma.w.inv[both.idx, both.idx] * TAB[j])
      }
    }
  } # p
  q.yy.a <- sum(PIJ * Y1w.c, na.rm = TRUE)
  PJ <- rowsum.default(PIJ[, both.idx], cluster.idx,
    reorder = FALSE,
    na.rm = TRUE
  ) # only both part is needed

  # per cluster
  both.diag.idx <- lav_matrix_diag_idx(length(both.idx))
  for (j in seq_len(nclusters)) {
    # we only need the 'both.idx' part of A.j, sigma.b.z, p.j, g.j ,...
    A.j <- ALIST[[j]]
    p.j <- PJ[j, ]
    if (nz > 0L) {
      sigma.b.z <- SIGMA.B.Z[[ZPAT2J[j]]]
    } else {
      sigma.b.z <- sigma.b
    }
    IBZA.j <- sigma.b.z %*% A.j
    IBZA.j[both.diag.idx] <- IBZA.j[both.diag.idx] + 1
    # logdet IBZA.j
    tmp <- determinant.matrix(IBZA.j, logarithm = TRUE)
    IBZA.j.logdet[j] <- tmp$modulus * tmp$sign
    # IBZA.j.inv.BZ.p
    IBZA.j.inv.BZ.p <- solve.default(IBZA.j, drop(sigma.b.z %*% p.j))
    q.yy.b[j] <- sum(p.j * IBZA.j.inv.BZ.p)

    if (nz > 0L) {
      g.j <- GJ[j, ]
      IBZA.j.inv.g <- solve.default(IBZA.j, g.j)
      A.IBZA.j.inv.g <- A.j %*% IBZA.j.inv.g

      q.zz.b[j] <- sum(g.j * A.IBZA.j.inv.g)
      q.zy[j] <- -sum(p.j * IBZA.j.inv.g)
    }
  }


  if (nz > 0L) {
    P <- Mp$nel + Zp$nel
    DIST <- (q.yy.a - sum(q.yy.b)) + 2 * sum(q.zy) + (q.zz.a + sum(q.zz.b))
    LOGDET <- W.logdet + sum(IBZA.j.logdet) + sigma.j.zz.logdet
  } else {
    P <- Mp$nel
    DIST <- (q.yy.a - sum(q.yy.b))
    LOGDET <- W.logdet + sum(IBZA.j.logdet)
  }

  # loglik?
  if (log2pi && !minus.two) {
    LOG.2PI <- log(2 * pi)
    loglik <- -(P * LOG.2PI + LOGDET + DIST) / 2
  } else {
    loglik <- DIST + LOGDET
  }

  # loglik.x (only if loglik is requested)
  if (length(unlist(Lp$ov.x.idx)) > 0L && log2pi && !minus.two) {
    loglik <- loglik - loglik.x
  }

  loglik
}

# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
lav_mvnorm_cluster_missing_dlogl_2l_samplestats <- function(
    Y1 = NULL,
    Y2 = NULL,
    Lp = NULL,
    Mp = NULL,
    Mu.W = NULL,
    Sigma.W = NULL,
    Mu.B = NULL,
    Sigma.B = NULL,
    Sinv.method = "eigen",
    return.list = FALSE) {
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

  # containers for dx
  dx.mu.y <- numeric(length(mu.y))
  dx.mu.z <- numeric(length(mu.z))
  dx.sigma.zz <- matrix(0, nrow(sigma.zz), ncol(sigma.zz))
  dx.sigma.yz <- matrix(0, nrow(sigma.yz), ncol(sigma.yz))
  dx.sigma.b <- matrix(0, nrow(sigma.b), ncol(sigma.b))
  dx.sigma.w <- matrix(0, nrow(sigma.w), ncol(sigma.w))

  # Lp
  nclusters <- Lp$nclusters[[2]]
  between.idx <- Lp$between.idx[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]
  both.idx <- Lp$both.idx[[2]]

  # sigma.w
  sigma.w.inv <- solve.default(sigma.w)
  sigma.b <- sigma.b[both.idx, both.idx] # only both part

  # y
  ny <- ncol(sigma.w)
  if (length(between.idx) > 0L) {
    Y1w <- Y1[, -between.idx, drop = FALSE]
  } else {
    Y1w <- Y1
  }
  Y1w.c <- t(t(Y1w) - mu.y)
  PIJ <- matrix(0, nrow(Y1w.c), ny)

  # z
  nz <- length(between.idx)
  if (nz > 0L) {
    Z <- Y2[, between.idx, drop = FALSE]
    Z.c <- t(t(Z) - mu.z)
    sigma.yz <- sigma.yz[both.idx, , drop = FALSE] # only both part
    sigma.zy <- t(sigma.yz)
    sigma.zz.inv <- solve.default(sigma.zz)
    sigma.zz.logdet <- log(det(sigma.zz))
    sigma.zi.zy <- sigma.zz.inv %*% sigma.zy
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy
    GZ <- Z.c %*% sigma.zz.inv # for complete cases only
  }

  # containters per cluster
  # ALIST <- rep(list(matrix(0, length(both.idx),
  #                             length(both.idx))), nclusters)
  ALIST <- rep(list(matrix(0, ny, ny)), nclusters)

  # Z per missing pattern
  if (nz > 0L) {
    Zp <- Mp$Zp
    ZPAT2J <- integer(nclusters) # which pattern per cluster
    SIGMA.B.Z <- vector("list", length = Zp$npatterns + 1L) # +1 for empty
    ZIZY <- rep(list(matrix(
      0, nrow(sigma.zy),
      ncol(sigma.zy)
    )), Zp$npatterns + 1L)
    ZIP <- rep(list(matrix(
      0, nrow(sigma.zz),
      ncol(sigma.zz)
    )), Zp$npatterns + 1L)
    for (p in seq_len(Zp$npatterns)) {
      freq <- Zp$freq[p]
      z.na.idx <- which(!Zp$pat[p, ])
      j.idx <- Zp$case.idx[[p]] # cluster indices with this pattern
      ZPAT2J[j.idx] <- p

      if (length(z.na.idx) > 0L) {
        zp <- sigma.zz[-z.na.idx, -z.na.idx, drop = FALSE]
        zp.inv <- lav_matrix_symmetric_inverse_update(
          S.inv = sigma.zz.inv, rm.idx = z.na.idx,
          logdet = FALSE
        )
        ZIP[[p]][-z.na.idx, -z.na.idx] <- zp.inv
        GZ[j.idx, -z.na.idx] <- Z.c[j.idx, -z.na.idx] %*% zp.inv
        Z.G.ZY <- zp.inv %*% sigma.zy[-z.na.idx, , drop = FALSE]
        ZIZY[[p]][-z.na.idx, ] <-
          zp.inv %*% sigma.zy[-z.na.idx, , drop = FALSE]
        yziy <- sigma.yz[, -z.na.idx, drop = FALSE] %*% Z.G.ZY
        SIGMA.B.Z[[p]] <- (sigma.b - yziy)
      } else {
        # complete case
        ZIZY[[p]] <- sigma.zi.zy
        ZIP[[p]] <- sigma.zz.inv
        SIGMA.B.Z[[p]] <- sigma.b.z
      }
    } # p

    # add empty patterns (if any)
    if (length(Zp$empty.idx) > 0L) {
      ZPAT2J[Zp$empty.idx] <- p + 1L
      SIGMA.B.Z[[p + 1L]] <- sigma.b
    }

    GZ[is.na(GZ)] <- 0
    GJ <- GZ %*% sigma.zy
  }

  # Y per missing pattern
  WIP <- rep(list(matrix(0, ny, ny)), Mp$npatterns)
  MPi <- integer(nrow(Y1))
  for (p in seq_len(Mp$npatterns)) {
    freq <- Mp$freq[p]
    na.idx <- which(!Mp$pat[p, ])
    j.idx <- Mp$j.idx[[p]]
    j1.idx <- Mp$j1.idx[[p]]
    TAB <- integer(nclusters)
    TAB[j1.idx] <- Mp$j.freq[[p]]

    if (length(na.idx) > 0L) {
      MPi[Mp$case.idx[[p]]] <- p
      wp.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = sigma.w.inv, rm.idx = na.idx,
        logdet = FALSE
      )
      WIP[[p]][-na.idx, -na.idx] <- wp.inv
      PIJ[Mp$case.idx[[p]], -na.idx] <-
        Y1w.c[Mp$case.idx[[p]], -na.idx] %*% wp.inv

      for (j in j1.idx) {
        ALIST[[j]] <-
          # ALIST[[j]] + (WIP[[p]][both.idx, both.idx] * TAB[j])
          ALIST[[j]] + (WIP[[p]] * TAB[j])
      }
    } else {
      # complete case
      PIJ[Mp$case.idx[[p]], ] <-
        Y1w.c[Mp$case.idx[[p]], ] %*% sigma.w.inv
      WIP[[p]] <- sigma.w.inv
      for (j in j1.idx) {
        ALIST[[j]] <-
          # ALIST[[j]] + (sigma.w.inv[both.idx, both.idx] * TAB[j])
          ALIST[[j]] + (sigma.w.inv * TAB[j])
      }
    }
  } # p
  PJ <- rowsum.default(PIJ[, , drop = FALSE],
    cluster.idx,
    reorder = FALSE, na.rm = TRUE
  )

  # per cluster
  both.diag.idx <- lav_matrix_diag_idx(length(both.idx))
  for (j in seq_len(nclusters)) {
    A.j.full <- ALIST[[j]]
    A.j <- A.j.full[both.idx, both.idx, drop = FALSE]
    p.j <- as.matrix(PJ[j, ])
    pb.j <- as.matrix(PJ[j, both.idx]) # only both.idx part
    if (nz > 0L) {
      sigma.b.z <- SIGMA.B.Z[[ZPAT2J[j]]]
    } else {
      sigma.b.z <- sigma.b
    }

    IBZA.j <- sigma.b.z %*% A.j
    IBZA.j[both.diag.idx] <- IBZA.j[both.diag.idx] + 1

    IBZA.j.inv.BZ <- solve.default(IBZA.j, sigma.b.z)
    IBZA.j.inv.BZ.p <- IBZA.j.inv.BZ %*% pb.j
    A.IBZA.j.inv.BZ <- A.j %*% IBZA.j.inv.BZ
    A.IBZA.j.inv.BZ.p <- A.IBZA.j.inv.BZ %*% pb.j

    IBZA.j.inv <- solve.default(IBZA.j)
    A.IBZA.j.inv <- A.j %*% IBZA.j.inv
    p.IBZA.j.inv <- t(crossprod(pb.j, IBZA.j.inv))

    # only if we have between-only variables
    if (nz > 0L) {
      g.j <- as.matrix(GJ[j, ])
      zij <- as.matrix(GZ[j, ])
      zizy <- ZIZY[[ZPAT2J[j]]]
      zip <- ZIP[[ZPAT2J[j]]]

      IBZA.j.inv.zizy <- solve.default(IBZA.j, t(zizy))
      IBZA.j.inv.g <- IBZA.j.inv %*% g.j
      IBZA.j.inv.p <- IBZA.j.inv %*% pb.j
      A.IBZA.j.inv.g <- A.j %*% IBZA.j.inv.g
      A.IBZA.j.inv.zizy <- A.j %*% IBZA.j.inv.zizy
      zizy.A.IBZA.j.inv.g <- zizy %*% A.IBZA.j.inv.g
      p.IBZA.j.inv.zizy <- crossprod(pb.j, IBZA.j.inv.zizy)
      ggbzpp <- 2 * A.IBZA.j.inv.g + A.IBZA.j.inv.BZ.p - pb.j
      ZIJzizyp <- (2 * zij - zizy %*% pb.j)

      ###########
      # dx.mu.z #
      ###########
      tmp <- 2 * (t(p.IBZA.j.inv.zizy) - zij - zizy.A.IBZA.j.inv.g)
      dx.mu.z <- dx.mu.z + drop(tmp)

      ###############
      # dx.sigma.zz #
      ###############
      tmp1 <- (zip + zizy %*% A.IBZA.j.inv.zizy # logdet
        - tcrossprod(zij) # ZA
        - tcrossprod(zizy.A.IBZA.j.inv.g)) # ZB-1

      d <- (t((2 * zizy.A.IBZA.j.inv.g + zizy %*% A.IBZA.j.inv.BZ.p)
      %*% p.IBZA.j.inv.zizy)
      + ZIJzizyp %*% p.IBZA.j.inv.zizy
        - 2 * tcrossprod(zizy.A.IBZA.j.inv.g, zij))
      tmp2 <- (d + t(d)) / 2
      tmp <- tmp1 + tmp2
      # symmetry correction
      ZZ <- 2 * tmp
      diag(ZZ) <- diag(tmp)
      dx.sigma.zz <- dx.sigma.zz + ZZ

      ###############
      # dx.sigma.yz #
      ###############
      t0 <- -2 * A.IBZA.j.inv.zizy

      t1 <- (-2 * tcrossprod(p.IBZA.j.inv, g.j)
        - 1 * tcrossprod(p.IBZA.j.inv, sigma.b.z %*% pb.j)
        + 2 * tcrossprod(A.IBZA.j.inv.g, g.j)) %*% A.IBZA.j.inv.zizy
      t2 <- -ggbzpp %*% p.IBZA.j.inv.zizy
      t3 <- -tcrossprod(p.IBZA.j.inv, ZIJzizyp)
      t4 <- 2 * tcrossprod(A.IBZA.j.inv.g, zij)
      tmp <- t0 + t1 + t2 + t3 + t4
      dx.sigma.yz[both.idx, ] <- dx.sigma.yz[both.idx, , drop = FALSE] + tmp

      ##############
      # dx.sigma.b #
      ##############
      c <- tcrossprod(ggbzpp, p.IBZA.j.inv)
      tmp <- t(A.IBZA.j.inv) - tcrossprod(A.IBZA.j.inv.g) + (c + t(c)) / 2
      # symmetry correction
      ZZ <- 2 * tmp
      diag(ZZ) <- diag(tmp)
      dx.sigma.b[both.idx, both.idx] <-
        dx.sigma.b[both.idx, both.idx, drop = FALSE] + ZZ

      # for dx.sigma.w
      PART1.b <- -1 * (IBZA.j.inv.g %*%
        (2 * t(IBZA.j.inv.BZ.p) + t(g.j) -
          t(g.j) %*% A.IBZA.j.inv.BZ)
        + IBZA.j.inv.BZ + tcrossprod(IBZA.j.inv.BZ.p))
      PART2.b <- 2 * (IBZA.j.inv.g + IBZA.j.inv.BZ.p) # vector
    } else {
      ##############
      # dx.sigma.b #
      ##############
      bzpp <- A.IBZA.j.inv.BZ.p - pb.j
      c <- tcrossprod(bzpp, p.IBZA.j.inv)
      tmp <- t(A.IBZA.j.inv) + (c + t(c)) / 2
      # symmetry correction
      ZZ <- 2 * tmp
      diag(ZZ) <- diag(tmp)
      dx.sigma.b[both.idx, both.idx] <-
        dx.sigma.b[both.idx, both.idx, drop = FALSE] + ZZ

      PART1.b <- -1 * (IBZA.j.inv.BZ + tcrossprod(IBZA.j.inv.BZ.p))
      PART2.b <- 2 * IBZA.j.inv.BZ.p # vector
    }

    ##############
    # dx.sigma.w #
    ##############

    PART1 <- matrix(0, ny, ny)
    PART1[both.idx, both.idx] <- PART1.b

    PART2 <- matrix(0, ny, 1L)
    PART2[both.idx, 1L] <- PART2.b

    ij.index <- which(cluster.idx == j)
    pij <- PIJ[ij.index, , drop = FALSE]

    which.compl <- which(MPi[ij.index] == 0L)
    which.incompl <- which(MPi[ij.index] != 0L)

    AP2 <- rep(list(sigma.w.inv %*% PART2), length(ij.index))
    AP1A.a <- AP1A.b <- matrix(0, ny, ny)
    # A.j.full <- matrix(0, ny, ny)
    if (length(which.compl) > 0L) {
      tmp <- (sigma.w.inv %*% PART1 %*% sigma.w.inv)
      AP1A.a <- tmp * length(which.compl)
      # A.j.full <- A.j.full + sigma.w.inv * length(which.compl)
    }
    if (length(which.incompl) > 0L) {
      p.idx <- MPi[ij.index][which.incompl]
      tmp <- lapply(WIP[p.idx], function(x) {
        x %*% PART1 %*% x
      })
      AP1A.b <- Reduce("+", tmp)
      AP2[which.incompl] <-
        lapply(WIP[p.idx], function(x) {
          x %*% PART2
        })
      # A.j.full <- A.j.full + Reduce("+", WIP[ p.idx ])
    }
    t1 <- AP1A.a + AP1A.b
    t2 <- (do.call("cbind", AP2) - t(pij)) %*% pij

    AA.wj <- t1 + t2

    tmp <- A.j.full + (AA.wj + t(AA.wj)) / 2
    # symmetry correction
    ZZ <- 2 * tmp
    diag(ZZ) <- diag(tmp)
    dx.sigma.w <- dx.sigma.w + ZZ

    ###########
    # dx.mu.y #
    ###########
    tmp <- numeric(ny)
    if (nz > 0L) {
      tmp[both.idx] <- IBZA.j.inv.g + IBZA.j.inv.BZ.p
    } else {
      tmp[both.idx] <- IBZA.j.inv.BZ.p
    }
    gbzpp <- A.j.full %*% tmp - p.j
    dx.mu.y <- dx.mu.y + drop(2 * gbzpp)
  } # j

  # rearrange
  dout <- lav_mvnorm_cluster_2l2implied(
    Lp = Lp,
    sigma.w = dx.sigma.w, sigma.b = dx.sigma.b,
    sigma.yz = dx.sigma.yz, sigma.zz = dx.sigma.zz,
    mu.y = dx.mu.y, mu.z = dx.mu.z
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
lav_mvnorm_cluster_missing_scores_2l <- function(
    Y1 = NULL,
    Y2 = NULL,
    Lp = NULL,
    Mp = NULL,
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
  between.idx <- Lp$between.idx[[2]]
  cluster.idx <- Lp$cluster.idx[[2]]
  both.idx <- Lp$both.idx[[2]]

  # sigma.w
  sigma.w.inv <- solve.default(sigma.w)
  sigma.b <- sigma.b[both.idx, both.idx] # only both part

  # y
  ny <- ncol(sigma.w)
  if (length(between.idx) > 0L) {
    Y1w <- Y1[, -between.idx, drop = FALSE]
  } else {
    Y1w <- Y1
  }
  Y1w.c <- t(t(Y1w) - mu.y)
  PIJ <- matrix(0, nrow(Y1w.c), ny)

  # z
  nz <- length(between.idx)
  if (nz > 0L) {
    Z <- Y2[, between.idx, drop = FALSE]
    Z.c <- t(t(Z) - mu.z)
    sigma.yz <- sigma.yz[both.idx, , drop = FALSE] # only both part
    sigma.zy <- t(sigma.yz)
    sigma.zz.inv <- solve.default(sigma.zz)
    sigma.zz.logdet <- log(det(sigma.zz))
    sigma.zi.zy <- sigma.zz.inv %*% sigma.zy
    sigma.b.z <- sigma.b - sigma.yz %*% sigma.zi.zy
    GZ <- Z.c %*% sigma.zz.inv # for complete cases only
  }

  # containters per cluster
  # ALIST <- rep(list(matrix(0, length(both.idx),
  #                             length(both.idx))), nclusters)
  ALIST <- rep(list(matrix(0, ny, ny)), nclusters)

  # both level-1 and level-2
  G.muy <- matrix(0, nclusters, length(mu.y))
  G.Sigma.w <- matrix(0, nclusters, length(lav_matrix_vech(sigma.w)))
  G.Sigma.b <- matrix(0, nclusters, length(lav_matrix_vech(out$sigma.b)))
  G.muz <- matrix(0, nclusters, length(mu.z))
  G.Sigma.zz <- matrix(0, nclusters, length(lav_matrix_vech(sigma.zz)))
  G.Sigma.yz <- matrix(0, nclusters, length(lav_matrix_vec(out$sigma.yz)))

  # Z per missing pattern
  if (nz > 0L) {
    Zp <- Mp$Zp
    ZPAT2J <- integer(nclusters) # which pattern per cluster
    SIGMA.B.Z <- vector("list", length = Zp$npatterns + 1L) # +1 for empty
    ZIZY <- rep(list(matrix(
      0, nrow(sigma.zy),
      ncol(sigma.zy)
    )), Zp$npatterns + 1L)
    ZIP <- rep(list(matrix(
      0, nrow(sigma.zz),
      ncol(sigma.zz)
    )), Zp$npatterns + 1L)
    for (p in seq_len(Zp$npatterns)) {
      freq <- Zp$freq[p]
      z.na.idx <- which(!Zp$pat[p, ])
      j.idx <- Zp$case.idx[[p]] # cluster indices with this pattern
      ZPAT2J[j.idx] <- p

      if (length(z.na.idx) > 0L) {
        zp <- sigma.zz[-z.na.idx, -z.na.idx, drop = FALSE]
        zp.inv <- lav_matrix_symmetric_inverse_update(
          S.inv = sigma.zz.inv, rm.idx = z.na.idx,
          logdet = FALSE
        )
        ZIP[[p]][-z.na.idx, -z.na.idx] <- zp.inv
        GZ[j.idx, -z.na.idx] <- Z.c[j.idx, -z.na.idx] %*% zp.inv
        Z.G.ZY <- zp.inv %*% sigma.zy[-z.na.idx, , drop = FALSE]
        ZIZY[[p]][-z.na.idx, ] <-
          zp.inv %*% sigma.zy[-z.na.idx, , drop = FALSE]
        yziy <- sigma.yz[, -z.na.idx, drop = FALSE] %*% Z.G.ZY
        SIGMA.B.Z[[p]] <- (sigma.b - yziy)
      } else {
        # complete case
        ZIZY[[p]] <- sigma.zi.zy
        ZIP[[p]] <- sigma.zz.inv
        SIGMA.B.Z[[p]] <- sigma.b.z
      }
    } # p

    # add empty patterns (if any)
    if (length(Zp$empty.idx) > 0L) {
      ZPAT2J[Zp$empty.idx] <- p + 1L
      SIGMA.B.Z[[p + 1L]] <- sigma.b
    }

    GZ[is.na(GZ)] <- 0
    GJ <- GZ %*% sigma.zy
  }

  # Y per missing pattern
  WIP <- rep(list(matrix(0, ny, ny)), Mp$npatterns)
  MPi <- integer(nrow(Y1))
  for (p in seq_len(Mp$npatterns)) {
    freq <- Mp$freq[p]
    na.idx <- which(!Mp$pat[p, ])
    j.idx <- Mp$j.idx[[p]]
    j1.idx <- Mp$j1.idx[[p]]
    TAB <- integer(nclusters)
    TAB[j1.idx] <- Mp$j.freq[[p]]

    if (length(na.idx) > 0L) {
      MPi[Mp$case.idx[[p]]] <- p
      wp.inv <- lav_matrix_symmetric_inverse_update(
        S.inv = sigma.w.inv, rm.idx = na.idx,
        logdet = FALSE
      )
      WIP[[p]][-na.idx, -na.idx] <- wp.inv
      PIJ[Mp$case.idx[[p]], -na.idx] <-
        Y1w.c[Mp$case.idx[[p]], -na.idx] %*% wp.inv

      for (j in j1.idx) {
        ALIST[[j]] <-
          # ALIST[[j]] + (WIP[[p]][both.idx, both.idx] * TAB[j])
          ALIST[[j]] + (WIP[[p]] * TAB[j])
      }
    } else {
      # complete case
      PIJ[Mp$case.idx[[p]], ] <-
        Y1w.c[Mp$case.idx[[p]], ] %*% sigma.w.inv
      WIP[[p]] <- sigma.w.inv
      for (j in j1.idx) {
        ALIST[[j]] <-
          # ALIST[[j]] + (sigma.w.inv[both.idx, both.idx] * TAB[j])
          ALIST[[j]] + (sigma.w.inv * TAB[j])
      }
    }
  } # p

  PJ <- rowsum.default(PIJ[, , drop = FALSE],
    cluster.idx,
    reorder = FALSE, na.rm = TRUE
  )

  # per cluster
  both.diag.idx <- lav_matrix_diag_idx(length(both.idx))
  for (j in seq_len(nclusters)) {
    A.j.full <- ALIST[[j]]
    A.j <- A.j.full[both.idx, both.idx, drop = FALSE]
    p.j <- as.matrix(PJ[j, ])
    pb.j <- as.matrix(PJ[j, both.idx]) # only both.idx part
    if (nz > 0L) {
      sigma.b.z <- SIGMA.B.Z[[ZPAT2J[j]]]
    } else {
      sigma.b.z <- sigma.b
    }

    IBZA.j <- sigma.b.z %*% A.j
    IBZA.j[both.diag.idx] <- IBZA.j[both.diag.idx] + 1

    IBZA.j.inv.BZ <- solve.default(IBZA.j, sigma.b.z)
    IBZA.j.inv.BZ.p <- IBZA.j.inv.BZ %*% pb.j
    A.IBZA.j.inv.BZ <- A.j %*% IBZA.j.inv.BZ
    A.IBZA.j.inv.BZ.p <- A.IBZA.j.inv.BZ %*% pb.j

    IBZA.j.inv <- solve.default(IBZA.j)
    A.IBZA.j.inv <- A.j %*% IBZA.j.inv
    p.IBZA.j.inv <- t(crossprod(pb.j, IBZA.j.inv))

    # only if we have between-only variables
    if (nz > 0L) {
      g.j <- as.matrix(GJ[j, ])
      zij <- as.matrix(GZ[j, ])
      zizy <- ZIZY[[ZPAT2J[j]]]
      zip <- ZIP[[ZPAT2J[j]]]

      IBZA.j.inv.zizy <- solve.default(IBZA.j, t(zizy))
      IBZA.j.inv.g <- IBZA.j.inv %*% g.j
      IBZA.j.inv.p <- IBZA.j.inv %*% pb.j
      A.IBZA.j.inv.g <- A.j %*% IBZA.j.inv.g
      A.IBZA.j.inv.zizy <- A.j %*% IBZA.j.inv.zizy
      zizy.A.IBZA.j.inv.g <- zizy %*% A.IBZA.j.inv.g
      p.IBZA.j.inv.zizy <- crossprod(pb.j, IBZA.j.inv.zizy)
      ggbzpp <- 2 * A.IBZA.j.inv.g + A.IBZA.j.inv.BZ.p - pb.j
      ZIJzizyp <- (2 * zij - zizy %*% pb.j)

      ###########
      # dx.mu.z #
      ###########
      tmp <- 2 * (t(p.IBZA.j.inv.zizy) - zij - zizy.A.IBZA.j.inv.g)
      G.muz[j, ] <- drop(tmp)

      ###############
      # dx.sigma.zz #
      ###############
      tmp1 <- (zip + zizy %*% A.IBZA.j.inv.zizy # logdet
        - tcrossprod(zij) # ZA
        - tcrossprod(zizy.A.IBZA.j.inv.g)) # ZB-1

      d <- (t((2 * zizy.A.IBZA.j.inv.g + zizy %*% A.IBZA.j.inv.BZ.p)
      %*% p.IBZA.j.inv.zizy)
      + ZIJzizyp %*% p.IBZA.j.inv.zizy
        - 2 * tcrossprod(zizy.A.IBZA.j.inv.g, zij))
      tmp2 <- (d + t(d)) / 2
      tmp <- tmp1 + tmp2
      # symmetry correction
      ZZ <- 2 * tmp
      diag(ZZ) <- diag(tmp)
      G.Sigma.zz[j, ] <- lav_matrix_vech(ZZ)

      ###############
      # dx.sigma.yz #
      ###############
      t0 <- -2 * A.IBZA.j.inv.zizy

      t1 <- (-2 * tcrossprod(p.IBZA.j.inv, g.j)
        - 1 * tcrossprod(p.IBZA.j.inv, sigma.b.z %*% pb.j)
        + 2 * tcrossprod(A.IBZA.j.inv.g, g.j)) %*% A.IBZA.j.inv.zizy
      t2 <- -ggbzpp %*% p.IBZA.j.inv.zizy
      t3 <- -tcrossprod(p.IBZA.j.inv, ZIJzizyp)
      t4 <- 2 * tcrossprod(A.IBZA.j.inv.g, zij)
      tmp <- t0 + t1 + t2 + t3 + t4
      tmp2 <- matrix(0, nrow(out$sigma.yz), ncol(out$sigma.yz))
      tmp2[both.idx, ] <- tmp
      G.Sigma.yz[j, ] <- lav_matrix_vec(tmp2)

      ##############
      # dx.sigma.b #
      ##############
      c <- tcrossprod(ggbzpp, p.IBZA.j.inv)
      tmp <- t(A.IBZA.j.inv) - tcrossprod(A.IBZA.j.inv.g) + (c + t(c)) / 2
      # symmetry correction
      ZZ <- 2 * tmp
      diag(ZZ) <- diag(tmp)
      ZZ2 <- matrix(0, nrow(out$sigma.b), ncol(out$sigma.b))
      ZZ2[both.idx, both.idx] <- ZZ
      G.Sigma.b[j, ] <- lav_matrix_vech(ZZ2)
      # dx.sigma.b[both.idx, both.idx] <-
      #  dx.sigma.b[both.idx, both.idx, drop = FALSE] + ZZ

      # for dx.sigma.w
      PART1.b <- -1 * (IBZA.j.inv.g %*%
        (2 * t(IBZA.j.inv.BZ.p) + t(g.j) -
          t(g.j) %*% A.IBZA.j.inv.BZ)
        + IBZA.j.inv.BZ + tcrossprod(IBZA.j.inv.BZ.p))
      PART2.b <- 2 * (IBZA.j.inv.g + IBZA.j.inv.BZ.p) # vector
    } else {
      ##############
      # dx.sigma.b #
      ##############
      bzpp <- A.IBZA.j.inv.BZ.p - pb.j
      c <- tcrossprod(bzpp, p.IBZA.j.inv)
      tmp <- t(A.IBZA.j.inv) + (c + t(c)) / 2
      # symmetry correction
      ZZ <- 2 * tmp
      diag(ZZ) <- diag(tmp)
      ZZ2 <- matrix(0, nrow(out$sigma.b), ncol(out$sigma.b))
      ZZ2[both.idx, both.idx] <- ZZ
      G.Sigma.b[j, ] <- lav_matrix_vech(ZZ2)
      # dx.sigma.b[both.idx, both.idx] <-
      #  dx.sigma.b[both.idx, both.idx, drop = FALSE] + ZZ

      PART1.b <- -1 * (IBZA.j.inv.BZ + tcrossprod(IBZA.j.inv.BZ.p))
      PART2.b <- 2 * IBZA.j.inv.BZ.p # vector
    }

    ##############
    # dx.sigma.w #
    ##############

    PART1 <- matrix(0, ny, ny)
    PART1[both.idx, both.idx] <- PART1.b

    PART2 <- matrix(0, ny, 1L)
    PART2[both.idx, 1L] <- PART2.b

    ij.index <- which(cluster.idx == j)
    pij <- PIJ[ij.index, , drop = FALSE]

    which.compl <- which(MPi[ij.index] == 0L)
    which.incompl <- which(MPi[ij.index] != 0L)

    AP2 <- rep(list(sigma.w.inv %*% PART2), length(ij.index))
    AP1A.a <- AP1A.b <- matrix(0, ny, ny)
    # A.j.full <- matrix(0, ny, ny)
    if (length(which.compl) > 0L) {
      tmp <- (sigma.w.inv %*% PART1 %*% sigma.w.inv)
      AP1A.a <- tmp * length(which.compl)
      # A.j.full <- A.j.full + sigma.w.inv * length(which.compl)
    }
    if (length(which.incompl) > 0L) {
      p.idx <- MPi[ij.index][which.incompl]
      tmp <- lapply(WIP[p.idx], function(x) {
        x %*% PART1 %*% x
      })
      AP1A.b <- Reduce("+", tmp)
      AP2[which.incompl] <-
        lapply(WIP[p.idx], function(x) {
          x %*% PART2
        })
      # A.j.full <- A.j.full + Reduce("+", WIP[ p.idx ])
    }
    t1 <- AP1A.a + AP1A.b
    t2 <- (do.call("cbind", AP2) - t(pij)) %*% pij

    AA.wj <- t1 + t2

    tmp <- A.j.full + (AA.wj + t(AA.wj)) / 2
    # symmetry correction
    ZZ <- 2 * tmp
    diag(ZZ) <- diag(tmp)
    G.Sigma.w[j, ] <- lav_matrix_vech(ZZ)
    # dx.sigma.w <- dx.sigma.w + ZZ

    ###########
    # dx.mu.y #
    ###########
    tmp <- numeric(ny)
    if (nz > 0L) {
      tmp[both.idx] <- IBZA.j.inv.g + IBZA.j.inv.BZ.p
    } else {
      tmp[both.idx] <- IBZA.j.inv.BZ.p
    }
    gbzpp <- A.j.full %*% tmp - p.j
    # dx.mu.y <- dx.mu.y + drop(2 * gbzpp)
    G.muy[j, ] <- drop(2 * gbzpp)
  } # j

  # browser()

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
lav_mvnorm_cluster_missing_information_firstorder <- function(
    Y1 = NULL,
    Y2 = NULL,
    Lp = NULL,
    Mp = NULL,
    Mu.W = NULL,
    Sigma.W = NULL,
    Mu.B = NULL,
    Sigma.B = NULL,
    x.idx = NULL,
    divide.by.two = FALSE,
    Sinv.method = "eigen") {
  N <- NROW(Y1)

  SCORES <- lav_mvnorm_cluster_missing_scores_2l(
    Y1 = Y1,
    Y2 = Y2,
    Lp = Lp,
    Mp = Mp,
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

# observed information
# order: mu.w within, vech(sigma.w) within, mu.b between, vech(sigma.b) between
# mu.w rows/cols that are splitted within/between are forced to zero
#
# numerical approximation (for now)
lav_mvnorm_cluster_missing_information_observed <- function(
    Y1 = NULL,
    Y2 = NULL,
    Lp = NULL,
    Mp = NULL,
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

    dx <- lav_mvnorm_cluster_missing_dlogl_2l_samplestats(
      Y1 = Y1, Y2 = Y2, Lp = Lp, Mp = Mp,
      Mu.W = Mu.W2, Sigma.W = Sigma.W2,
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
