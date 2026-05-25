# loglikelihood clustered/twolevel data in the presence of missing data

# YR:
# - objective function: first version around March 2021 (see Psych paper)
# - analytic gradient: first version around May 2021


# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
lav_mvn_cl_mi_loglik_samp_2l <- function(y1 = NULL,
                                                         y2 = NULL,
                                                         lp = NULL,
                                                         mp = NULL,
                                                         mu_w = NULL,
                                                         sigma_w = NULL,
                                                         mu_b = NULL,
                                                         sigma_b = NULL,
                                                         sinv_method = "eigen",
                                                         log2pi = FALSE,
                                                         loglik_x = 0,
                                                         minus_two = TRUE) {
  # map implied to 2l matrices
  out <- lav_mvn_cl_implied22l(
    lp = lp, mu_w = mu_w, mu_b = mu_b,
    sigma_w = sigma_w, sigma_b = sigma_b
  )
  mu_y <- out$mu.y
  mu_z <- out$mu.z
  sigma_w_1 <- out$sigma.w
  sigma_b_1 <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz

  # Lp
  nclusters <- lp$nclusters[[2]]
  between_idx <- lp$between.idx[[2]]
  both_idx <- lp$both.idx[[2]]
  cluster_idx <- lp$cluster.idx[[2]]

  # sanity checks
  if (any(diag(sigma_w_1) < 0) || any(diag(sigma_b_1) < 0)) {
    return(+Inf)
  }

  # check is both.idx part of sigma.b is 'too' negative; if so, return +Inf
  ev <- eigen(sigma_b_1[both_idx, both_idx, drop = FALSE],
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
  sigma_w_inv <- solve.default(sigma_w_1)
  sigma_w_logdet <- log(det(sigma_w_1))
  sigma_b_1 <- sigma_b_1[both_idx, both_idx] # only both part

  # y
  ny <- ncol(sigma_w_1)
  if (length(between_idx) > 0L) {
    y1w <- y1[, -between_idx, drop = FALSE]
  } else {
    y1w <- y1
  }
  y1w_c <- t(t(y1w) - mu_y)
  pij <- matrix(0, nrow(y1w_c), ny)

  # z
  nz <- length(between_idx)
  if (nz > 0L) {
    # check is sigma_zz is PD; if not, return +Inf
    ev <- eigen(sigma_zz, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < sqrt(.Machine$double.eps))) {
      return(+Inf)
    }

    z <- y2[, between_idx, drop = FALSE]
    z_c <- t(t(z) - mu_z)

    sigma_yz <- sigma_yz[both_idx, , drop = FALSE] # only both part
    sigma_zy <- t(sigma_yz)
    sigma_zz_inv <- solve.default(sigma_zz)
    sigma_zz_logdet <- log(det(sigma_zz))
    sigma_zi_zy <- sigma_zz_inv %*% sigma_zy
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy
    gz <- z_c %*% sigma_zz_inv # for complete cases only
  }

  # containters per cluster
  q_yy_b <- q_zy <- q_zz_b <- numeric(nclusters)
  ibza_j_logdet <- numeric(nclusters)
  alist_1 <- rep(list(matrix(
    0, length(both_idx),
    length(both_idx)
  )), nclusters)

  # Z per missing pattern
  if (nz > 0L) {
    zp_1 <- mp$Zp
    zpat2j <- integer(nclusters) # which sigma.b.z per cluster
    sigma_b_z_1 <- vector("list", length = zp_1$npatterns + 1L)

    sigma_j_zz_logdet <- q_zz_a <- 0
    for (p in seq_len(zp_1$npatterns)) {
      freq <- zp_1$freq[p]
      z_na_idx <- which(!zp_1$pat[p, ])
      j_idx <- zp_1$case.idx[[p]] # cluster indices with this pattern
      zpat2j[j_idx] <- p

      if (length(z_na_idx) > 0L) {
        # zp <- sigma_zz[-z_na_idx, -z_na_idx, drop = FALSE]
        zp_inv <- lav_mat_sym_inverse_update(
          s_inv = sigma_zz_inv, rm_idx = z_na_idx,
          logdet = TRUE, s_logdet = sigma_zz_logdet
        )
        zp_logdet <- attr(zp_inv, "logdet")
        sigma_j_zz_logdet <- sigma_j_zz_logdet + (zp_logdet * freq)

        gz[j_idx, -z_na_idx] <- z_c[j_idx, -z_na_idx] %*% zp_inv

        yziy <- (sigma_yz[, -z_na_idx, drop = FALSE] %*% zp_inv %*%
          sigma_zy[-z_na_idx, , drop = FALSE])
        sigma_b_z_1[[p]] <- (sigma_b_1 - yziy)
      } else {
        # complete case
        sigma_j_zz_logdet <-
          sigma_j_zz_logdet + (sigma_zz_logdet * freq)
        sigma_b_z_1[[p]] <- sigma_b_z
      }
    } # p

    # add empty patterns (if any)
    if (length(zp_1$empty.idx) > 0L) {
      zpat2j[zp_1$empty.idx] <- p + 1L
      sigma_b_z_1[[p + 1L]] <- sigma_b_1
    }

    q_zz_a <- sum(gz * z_c, na.rm = TRUE)
    gz0 <- gz
    gz0[is.na(gz0)] <- 0
    gj <- gz0 %*% sigma_zy # only both part
  }

  # Y per missing pattern
  w_logdet <- 0
  #MPi <- integer(nrow(Y1))
  for (p in seq_len(mp$npatterns)) {
    freq <- mp$freq[p]
    na_idx <- which(!mp$pat[p, ])
    j_idx <- mp$j.idx[[p]]
    j1_idx <- mp$j1.idx[[p]]
    tab <- integer(nclusters)
    tab[j1_idx] <- mp$j.freq[[p]]

    # compute sigma.w.inv for this pattern
    if (length(na_idx) > 0L) {
      #MPi[Mp$case.idx[[p]]] <- p
      # wp <- sigma_w_1[-na_idx, -na_idx, drop = FALSE]
      wp_inv <- lav_mat_sym_inverse_update(
        s_inv = sigma_w_inv, rm_idx = na_idx,
        logdet = TRUE, s_logdet = sigma_w_logdet
      )
      wp_logdet <- attr(wp_inv, "logdet")
      w_logdet <- w_logdet + (wp_logdet * freq)

      pij[mp$case.idx[[p]], -na_idx] <-
        y1w_c[mp$case.idx[[p]], -na_idx] %*% wp_inv

      a_j <- matrix(0, ny, ny)
      a_j[-na_idx, -na_idx] <- wp_inv
      for (j in j1_idx) {
        alist_1[[j]] <- alist_1[[j]] + (a_j[both_idx, both_idx] * tab[j])
      }
      # WIP[[p]][-na.idx, -na.idx] <- wp.inv
    } else {
      # complete case
      w_logdet <- w_logdet + (sigma_w_logdet * freq)
      pij[mp$case.idx[[p]], ] <-
        y1w_c[mp$case.idx[[p]], ] %*% sigma_w_inv
      for (j in j1_idx) {
        alist_1[[j]] <-
          alist_1[[j]] + (sigma_w_inv[both_idx, both_idx] * tab[j])
      }
    }
  } # p
  q_yy_a <- sum(pij * y1w_c, na.rm = TRUE)
  pj <- rowsum.default(pij[, both_idx], cluster_idx,
    reorder = FALSE,
    na.rm = TRUE
  ) # only both part is needed

  # per cluster
  both_diag_idx <- lav_mat_diag_idx(length(both_idx))
  for (j in seq_len(nclusters)) {
    # we only need the 'both.idx' part of A.j, sigma.b.z, p.j, g.j ,...
    a_j <- alist_1[[j]]
    p_j <- pj[j, ]
    if (nz > 0L) {
      sigma_b_z <- sigma_b_z_1[[zpat2j[j]]]
    } else {
      sigma_b_z <- sigma_b_1
    }
    ibza_j <- sigma_b_z %*% a_j
    ibza_j[both_diag_idx] <- ibza_j[both_diag_idx] + 1
    # logdet IBZA.j
    tmp <- determinant.matrix(ibza_j, logarithm = TRUE)
    ibza_j_logdet[j] <- tmp$modulus * tmp$sign
    # IBZA.j.inv.BZ.p
    ibza_j_inv_bz_p <- solve.default(ibza_j, drop(sigma_b_z %*% p_j))
    q_yy_b[j] <- sum(p_j * ibza_j_inv_bz_p)

    if (nz > 0L) {
      g_j <- gj[j, ]
      ibza_j_inv_g <- solve.default(ibza_j, g_j)
      a_ibza_j_inv_g <- a_j %*% ibza_j_inv_g

      q_zz_b[j] <- sum(g_j * a_ibza_j_inv_g)
      q_zy[j] <- -sum(p_j * ibza_j_inv_g)
    }
  }


  if (nz > 0L) {
    p_1 <- mp$nel + zp_1$nel
    dist_1 <- (q_yy_a - sum(q_yy_b)) + 2 * sum(q_zy) + (q_zz_a + sum(q_zz_b))
    logdet <- w_logdet + sum(ibza_j_logdet) + sigma_j_zz_logdet
  } else {
    p_1 <- mp$nel
    dist_1 <- (q_yy_a - sum(q_yy_b))
    logdet <- w_logdet + sum(ibza_j_logdet)
  }

  # loglik?
  if (log2pi && !minus_two) {
    log_2pi <- log(2 * pi)
    loglik <- -(p_1 * log_2pi + logdet + dist_1) / 2
  } else {
    loglik <- dist_1 + logdet
  }

  # loglik.x (only if loglik is requested)
  if (length(unlist(lp$ov.x.idx)) > 0L && log2pi && !minus_two) {
    loglik <- loglik - loglik_x
  }

  loglik
}

# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
lav_mvn_cl_mi_dlogl_2l_samp <- function(
    y1 = NULL,
    y2 = NULL,
    lp = NULL,
    mp = NULL,
    mu_w = NULL,
    sigma_w = NULL,
    mu_b = NULL,
    sigma_b = NULL,
    sinv_method = "eigen",
    return_list = FALSE) {
  # map implied to 2l matrices
  out <- lav_mvn_cl_implied22l(
    lp = lp, mu_w = mu_w, mu_b = mu_b,
    sigma_w = sigma_w, sigma_b = sigma_b
  )
  mu_y <- out$mu.y
  mu_z <- out$mu.z
  sigma_w_1 <- out$sigma.w
  sigma_b_1 <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz

  # containers for dx
  dx_mu_y <- numeric(length(mu_y))
  dx_mu_z <- numeric(length(mu_z))
  dx_sigma_zz <- matrix(0, nrow(sigma_zz), ncol(sigma_zz))
  dx_sigma_yz <- matrix(0, nrow(sigma_yz), ncol(sigma_yz))
  dx_sigma_b <- matrix(0, nrow(sigma_b_1), ncol(sigma_b_1))
  dx_sigma_w <- matrix(0, nrow(sigma_w_1), ncol(sigma_w_1))

  # Lp
  nclusters <- lp$nclusters[[2]]
  between_idx <- lp$between.idx[[2]]
  cluster_idx <- lp$cluster.idx[[2]]
  both_idx <- lp$both.idx[[2]]

  # sigma.w
  sigma_w_inv <- solve.default(sigma_w_1)
  sigma_b_1 <- sigma_b_1[both_idx, both_idx] # only both part

  # y
  ny <- ncol(sigma_w_1)
  if (length(between_idx) > 0L) {
    y1w <- y1[, -between_idx, drop = FALSE]
  } else {
    y1w <- y1
  }
  y1w_c <- t(t(y1w) - mu_y)
  pij_1 <- matrix(0, nrow(y1w_c), ny)

  # z
  nz <- length(between_idx)
  if (nz > 0L) {
    z <- y2[, between_idx, drop = FALSE]
    z_c <- t(t(z) - mu_z)
    sigma_yz <- sigma_yz[both_idx, , drop = FALSE] # only both part
    sigma_zy <- t(sigma_yz)
    sigma_zz_inv <- solve.default(sigma_zz)
    # sigma_zz_logdet <- log(det(sigma_zz))
    sigma_zi_zy <- sigma_zz_inv %*% sigma_zy
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy
    gz <- z_c %*% sigma_zz_inv # for complete cases only
  }

  # containters per cluster
  # ALIST <- rep(list(matrix(0, length(both.idx),
  #                             length(both.idx))), nclusters)
  alist_1 <- rep(list(matrix(0, ny, ny)), nclusters)

  # Z per missing pattern
  if (nz > 0L) {
    zp_1 <- mp$Zp
    zpat2j <- integer(nclusters) # which pattern per cluster
    sigma_b_z_1 <- vector("list", length = zp_1$npatterns + 1L) # +1 for empty
    zizy_1 <- rep(list(matrix(
      0, nrow(sigma_zy),
      ncol(sigma_zy)
    )), zp_1$npatterns + 1L)
    zip_1 <- rep(list(matrix(
      0, nrow(sigma_zz),
      ncol(sigma_zz)
    )), zp_1$npatterns + 1L)
    for (p in seq_len(zp_1$npatterns)) {
      # freq <- zp_1$freq[p]
      z_na_idx <- which(!zp_1$pat[p, ])
      j_idx <- zp_1$case.idx[[p]] # cluster indices with this pattern
      zpat2j[j_idx] <- p

      if (length(z_na_idx) > 0L) {
        # zp <- sigma_zz[-z_na_idx, -z_na_idx, drop = FALSE]
        zp_inv <- lav_mat_sym_inverse_update(
          s_inv = sigma_zz_inv, rm_idx = z_na_idx,
          logdet = FALSE
        )
        zip_1[[p]][-z_na_idx, -z_na_idx] <- zp_inv
        gz[j_idx, -z_na_idx] <- z_c[j_idx, -z_na_idx] %*% zp_inv
        z_g_zy <- zp_inv %*% sigma_zy[-z_na_idx, , drop = FALSE]
        zizy_1[[p]][-z_na_idx, ] <-
          zp_inv %*% sigma_zy[-z_na_idx, , drop = FALSE]
        yziy <- sigma_yz[, -z_na_idx, drop = FALSE] %*% z_g_zy
        sigma_b_z_1[[p]] <- (sigma_b_1 - yziy)
      } else {
        # complete case
        zizy_1[[p]] <- sigma_zi_zy
        zip_1[[p]] <- sigma_zz_inv
        sigma_b_z_1[[p]] <- sigma_b_z
      }
    } # p

    # add empty patterns (if any)
    if (length(zp_1$empty.idx) > 0L) {
      zpat2j[zp_1$empty.idx] <- p + 1L
      sigma_b_z_1[[p + 1L]] <- sigma_b_1
    }

    gz[is.na(gz)] <- 0
    gj <- gz %*% sigma_zy
  }

  # Y per missing pattern
  wip <- rep(list(matrix(0, ny, ny)), mp$npatterns)
  mpi <- integer(nrow(y1))
  for (p in seq_len(mp$npatterns)) {
    # freq <- mp$freq[p]
    na_idx <- which(!mp$pat[p, ])
    j_idx <- mp$j.idx[[p]]
    j1_idx <- mp$j1.idx[[p]]
    tab <- integer(nclusters)
    tab[j1_idx] <- mp$j.freq[[p]]

    if (length(na_idx) > 0L) {
      mpi[mp$case.idx[[p]]] <- p
      wp_inv <- lav_mat_sym_inverse_update(
        s_inv = sigma_w_inv, rm_idx = na_idx,
        logdet = FALSE
      )
      wip[[p]][-na_idx, -na_idx] <- wp_inv
      pij_1[mp$case.idx[[p]], -na_idx] <-
        y1w_c[mp$case.idx[[p]], -na_idx] %*% wp_inv

      for (j in j1_idx) {
        alist_1[[j]] <-
          # ALIST[[j]] + (WIP[[p]][both.idx, both.idx] * TAB[j])
          alist_1[[j]] + (wip[[p]] * tab[j])
      }
    } else {
      # complete case
      pij_1[mp$case.idx[[p]], ] <-
        y1w_c[mp$case.idx[[p]], ] %*% sigma_w_inv
      wip[[p]] <- sigma_w_inv
      for (j in j1_idx) {
        alist_1[[j]] <-
          # ALIST[[j]] + (sigma.w.inv[both.idx, both.idx] * TAB[j])
          alist_1[[j]] + (sigma_w_inv * tab[j])
      }
    }
  } # p
  pj <- rowsum.default(pij_1[, , drop = FALSE],
    cluster_idx,
    reorder = FALSE, na.rm = TRUE
  )

  # per cluster
  both_diag_idx <- lav_mat_diag_idx(length(both_idx))
  for (j in seq_len(nclusters)) {
    a_j_full <- alist_1[[j]]
    a_j <- a_j_full[both_idx, both_idx, drop = FALSE]
    p_j <- as.matrix(pj[j, ])
    pb_j <- as.matrix(pj[j, both_idx]) # only both.idx part
    if (nz > 0L) {
      sigma_b_z <- sigma_b_z_1[[zpat2j[j]]]
    } else {
      sigma_b_z <- sigma_b_1
    }

    ibza_j <- sigma_b_z %*% a_j
    ibza_j[both_diag_idx] <- ibza_j[both_diag_idx] + 1

    ibza_j_inv_bz <- solve.default(ibza_j, sigma_b_z)
    ibza_j_inv_bz_p <- ibza_j_inv_bz %*% pb_j
    a_ibza_j_inv_bz <- a_j %*% ibza_j_inv_bz
    a_ibza_j_inv_bz_p <- a_ibza_j_inv_bz %*% pb_j

    ibza_j_inv <- solve.default(ibza_j)
    a_ibza_j_inv <- a_j %*% ibza_j_inv
    p_ibza_j_inv <- t(crossprod(pb_j, ibza_j_inv))

    # only if we have between-only variables
    if (nz > 0L) {
      g_j <- as.matrix(gj[j, ])
      zij <- as.matrix(gz[j, ])
      zizy <- zizy_1[[zpat2j[j]]]
      zip <- zip_1[[zpat2j[j]]]

      ibza_j_inv_zizy <- solve.default(ibza_j, t(zizy))
      ibza_j_inv_g <- ibza_j_inv %*% g_j
      # ibza_j_inv_p <- ibza_j_inv %*% pb_j
      a_ibza_j_inv_g <- a_j %*% ibza_j_inv_g
      a_ibza_j_inv_zizy <- a_j %*% ibza_j_inv_zizy
      zizy_a_ibza_j_inv_g <- zizy %*% a_ibza_j_inv_g
      p_ibza_j_inv_zizy <- crossprod(pb_j, ibza_j_inv_zizy)
      ggbzpp <- 2 * a_ibza_j_inv_g + a_ibza_j_inv_bz_p - pb_j
      zijzizyp <- (2 * zij - zizy %*% pb_j)

      ###########
      # dx.mu.z #
      ###########
      tmp <- 2 * (t(p_ibza_j_inv_zizy) - zij - zizy_a_ibza_j_inv_g)
      dx_mu_z <- dx_mu_z + drop(tmp)

      ###############
      # dx.sigma.zz #
      ###############
      tmp1 <- (zip + zizy %*% a_ibza_j_inv_zizy # logdet
        - tcrossprod(zij) # ZA
        - tcrossprod(zizy_a_ibza_j_inv_g)) # ZB-1

      d <- (t((2 * zizy_a_ibza_j_inv_g + zizy %*% a_ibza_j_inv_bz_p)
      %*% p_ibza_j_inv_zizy)
      + zijzizyp %*% p_ibza_j_inv_zizy
        - 2 * tcrossprod(zizy_a_ibza_j_inv_g, zij))
      tmp2 <- (d + t(d)) / 2
      tmp <- tmp1 + tmp2
      # symmetry correction
      zz <- 2 * tmp
      diag(zz) <- diag(tmp)
      dx_sigma_zz <- dx_sigma_zz + zz

      ###############
      # dx.sigma.yz #
      ###############
      t0 <- -2 * a_ibza_j_inv_zizy

      t1 <- (-2 * tcrossprod(p_ibza_j_inv, g_j)
        - 1 * tcrossprod(p_ibza_j_inv, sigma_b_z %*% pb_j)
        + 2 * tcrossprod(a_ibza_j_inv_g, g_j)) %*% a_ibza_j_inv_zizy
      t2 <- -ggbzpp %*% p_ibza_j_inv_zizy
      t3 <- -tcrossprod(p_ibza_j_inv, zijzizyp)
      t4 <- 2 * tcrossprod(a_ibza_j_inv_g, zij)
      tmp <- t0 + t1 + t2 + t3 + t4
      dx_sigma_yz[both_idx, ] <- dx_sigma_yz[both_idx, , drop = FALSE] + tmp

      ##############
      # dx.sigma.b #
      ##############
      c <- tcrossprod(ggbzpp, p_ibza_j_inv)
      tmp <- t(a_ibza_j_inv) - tcrossprod(a_ibza_j_inv_g) + (c + t(c)) / 2
      # symmetry correction
      zz <- 2 * tmp
      diag(zz) <- diag(tmp)
      dx_sigma_b[both_idx, both_idx] <-
        dx_sigma_b[both_idx, both_idx, drop = FALSE] + zz

      # for dx.sigma.w
      part1_b <- -1 * (ibza_j_inv_g %*%
        (2 * t(ibza_j_inv_bz_p) + t(g_j) -
          t(g_j) %*% a_ibza_j_inv_bz)
        + ibza_j_inv_bz + tcrossprod(ibza_j_inv_bz_p))
      part2_b <- 2 * (ibza_j_inv_g + ibza_j_inv_bz_p) # vector
    } else {
      ##############
      # dx.sigma.b #
      ##############
      bzpp <- a_ibza_j_inv_bz_p - pb_j
      c <- tcrossprod(bzpp, p_ibza_j_inv)
      tmp <- t(a_ibza_j_inv) + (c + t(c)) / 2
      # symmetry correction
      zz <- 2 * tmp
      diag(zz) <- diag(tmp)
      dx_sigma_b[both_idx, both_idx] <-
        dx_sigma_b[both_idx, both_idx, drop = FALSE] + zz

      part1_b <- -1 * (ibza_j_inv_bz + tcrossprod(ibza_j_inv_bz_p))
      part2_b <- 2 * ibza_j_inv_bz_p # vector
    }

    ##############
    # dx.sigma.w #
    ##############

    part1 <- matrix(0, ny, ny)
    part1[both_idx, both_idx] <- part1_b

    part2 <- matrix(0, ny, 1L)
    part2[both_idx, 1L] <- part2_b

    ij_index <- which(cluster_idx == j)
    pij <- pij_1[ij_index, , drop = FALSE]

    which_compl <- which(mpi[ij_index] == 0L)
    which_incompl <- which(mpi[ij_index] != 0L)

    ap2 <- rep(list(sigma_w_inv %*% part2), length(ij_index))
    ap1a_a <- ap1a_b <- matrix(0, ny, ny)
    # A.j.full <- matrix(0, ny, ny)
    if (length(which_compl) > 0L) {
      tmp <- (sigma_w_inv %*% part1 %*% sigma_w_inv)
      ap1a_a <- tmp * length(which_compl)
      # A.j.full <- A.j.full + sigma.w.inv * length(which.compl)
    }
    if (length(which_incompl) > 0L) {
      p_idx <- mpi[ij_index][which_incompl]
      tmp <- lapply(wip[p_idx], function(x) {
        x %*% part1 %*% x
      })
      ap1a_b <- Reduce("+", tmp)
      ap2[which_incompl] <-
        lapply(wip[p_idx], function(x) {
          x %*% part2
        })
      # A.j.full <- A.j.full + Reduce("+", WIP[ p.idx ])
    }
    t1 <- ap1a_a + ap1a_b
    t2 <- (do.call("cbind", ap2) - t(pij)) %*% pij

    aa_wj <- t1 + t2

    tmp <- a_j_full + (aa_wj + t(aa_wj)) / 2
    # symmetry correction
    zz <- 2 * tmp
    diag(zz) <- diag(tmp)
    dx_sigma_w <- dx_sigma_w + zz

    ###########
    # dx.mu.y #
    ###########
    tmp <- numeric(ny)
    if (nz > 0L) {
      tmp[both_idx] <- ibza_j_inv_g + ibza_j_inv_bz_p
    } else {
      tmp[both_idx] <- ibza_j_inv_bz_p
    }
    gbzpp <- a_j_full %*% tmp - p_j
    dx_mu_y <- dx_mu_y + drop(2 * gbzpp)
  } # j

  # rearrange
  dout <- lav_mvn_cl_2l2implied(
    lp = lp,
    sigma_w = dx_sigma_w, sigma_b = dx_sigma_b,
    sigma_yz = dx_sigma_yz, sigma_zz = dx_sigma_zz,
    mu_y = dx_mu_y, mu_z = dx_mu_z
  )

  if (return_list) {
    out <- dout
  } else {
    out <- c(
      dout$Mu.W, lav_mat_vech(dout$Sigma.W),
      dout$Mu.B, lav_mat_vech(dout$Sigma.B)
    )
  }

  out
}

# cluster-wise scores -2*logl wrt Mu.W, Mu.B, Sigma.W, Sigma.B
lav_mvn_cl_mi_sc_2l <- function(
    y1 = NULL,
    y2 = NULL,
    lp = NULL,
    mp = NULL,
    mu_w = NULL,
    sigma_w = NULL,
    mu_b = NULL,
    sigma_b = NULL,
    sinv_method = "eigen") {
  # map implied to 2l matrices
  out <- lav_mvn_cl_implied22l(
    lp = lp, mu_w = mu_w, mu_b = mu_b,
    sigma_w = sigma_w, sigma_b = sigma_b
  )
  mu_y <- out$mu.y
  mu_z <- out$mu.z
  sigma_w_1 <- out$sigma.w
  sigma_b_1 <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz

  # Lp
  nclusters <- lp$nclusters[[2]]
  between_idx <- lp$between.idx[[2]]
  cluster_idx <- lp$cluster.idx[[2]]
  both_idx <- lp$both.idx[[2]]

  # sigma.w
  sigma_w_inv <- solve.default(sigma_w_1)
  sigma_b_1 <- sigma_b_1[both_idx, both_idx] # only both part

  # y
  ny <- ncol(sigma_w_1)
  if (length(between_idx) > 0L) {
    y1w <- y1[, -between_idx, drop = FALSE]
  } else {
    y1w <- y1
  }
  y1w_c <- t(t(y1w) - mu_y)
  pij_1 <- matrix(0, nrow(y1w_c), ny)

  # z
  nz <- length(between_idx)
  if (nz > 0L) {
    z <- y2[, between_idx, drop = FALSE]
    z_c <- t(t(z) - mu_z)
    sigma_yz <- sigma_yz[both_idx, , drop = FALSE] # only both part
    sigma_zy <- t(sigma_yz)
    sigma_zz_inv <- solve.default(sigma_zz)
    # sigma_zz_logdet <- log(det(sigma_zz))
    sigma_zi_zy <- sigma_zz_inv %*% sigma_zy
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy
    gz <- z_c %*% sigma_zz_inv # for complete cases only
  }

  # containters per cluster
  # ALIST <- rep(list(matrix(0, length(both.idx),
  #                             length(both.idx))), nclusters)
  alist_1 <- rep(list(matrix(0, ny, ny)), nclusters)

  # both level-1 and level-2
  g_muy <- matrix(0, nclusters, length(mu_y))
  g_sigma_w <- matrix(0, nclusters, length(lav_mat_vech(sigma_w_1)))
  g_sigma_b <- matrix(0, nclusters, length(lav_mat_vech(out$sigma.b)))
  g_muz <- matrix(0, nclusters, length(mu_z))
  g_sigma_zz <- matrix(0, nclusters, length(lav_mat_vech(sigma_zz)))
  g_sigma_yz <- matrix(0, nclusters, length(lav_mat_vec(out$sigma.yz)))

  # Z per missing pattern
  if (nz > 0L) {
    zp_1 <- mp$Zp
    zpat2j <- integer(nclusters) # which pattern per cluster
    sigma_b_z_1 <- vector("list", length = zp_1$npatterns + 1L) # +1 for empty
    zizy_1 <- rep(list(matrix(
      0, nrow(sigma_zy),
      ncol(sigma_zy)
    )), zp_1$npatterns + 1L)
    zip_1 <- rep(list(matrix(
      0, nrow(sigma_zz),
      ncol(sigma_zz)
    )), zp_1$npatterns + 1L)
    for (p in seq_len(zp_1$npatterns)) {
      # freq <- zp_1$freq[p]
      z_na_idx <- which(!zp_1$pat[p, ])
      j_idx <- zp_1$case.idx[[p]] # cluster indices with this pattern
      zpat2j[j_idx] <- p

      if (length(z_na_idx) > 0L) {
        # zp <- sigma_zz[-z_na_idx, -z_na_idx, drop = FALSE]
        zp_inv <- lav_mat_sym_inverse_update(
          s_inv = sigma_zz_inv, rm_idx = z_na_idx,
          logdet = FALSE
        )
        zip_1[[p]][-z_na_idx, -z_na_idx] <- zp_inv
        gz[j_idx, -z_na_idx] <- z_c[j_idx, -z_na_idx] %*% zp_inv
        z_g_zy <- zp_inv %*% sigma_zy[-z_na_idx, , drop = FALSE]
        zizy_1[[p]][-z_na_idx, ] <-
          zp_inv %*% sigma_zy[-z_na_idx, , drop = FALSE]
        yziy <- sigma_yz[, -z_na_idx, drop = FALSE] %*% z_g_zy
        sigma_b_z_1[[p]] <- (sigma_b_1 - yziy)
      } else {
        # complete case
        zizy_1[[p]] <- sigma_zi_zy
        zip_1[[p]] <- sigma_zz_inv
        sigma_b_z_1[[p]] <- sigma_b_z
      }
    } # p

    # add empty patterns (if any)
    if (length(zp_1$empty.idx) > 0L) {
      zpat2j[zp_1$empty.idx] <- p + 1L
      sigma_b_z_1[[p + 1L]] <- sigma_b_1
    }

    gz[is.na(gz)] <- 0
    gj <- gz %*% sigma_zy
  }

  # Y per missing pattern
  wip <- rep(list(matrix(0, ny, ny)), mp$npatterns)
  mpi <- integer(nrow(y1))
  for (p in seq_len(mp$npatterns)) {
    # freq <- mp$freq[p]
    na_idx <- which(!mp$pat[p, ])
    j_idx <- mp$j.idx[[p]]
    j1_idx <- mp$j1.idx[[p]]
    tab <- integer(nclusters)
    tab[j1_idx] <- mp$j.freq[[p]]

    if (length(na_idx) > 0L) {
      mpi[mp$case.idx[[p]]] <- p
      wp_inv <- lav_mat_sym_inverse_update(
        s_inv = sigma_w_inv, rm_idx = na_idx,
        logdet = FALSE
      )
      wip[[p]][-na_idx, -na_idx] <- wp_inv
      pij_1[mp$case.idx[[p]], -na_idx] <-
        y1w_c[mp$case.idx[[p]], -na_idx] %*% wp_inv

      for (j in j1_idx) {
        alist_1[[j]] <-
          # ALIST[[j]] + (WIP[[p]][both.idx, both.idx] * TAB[j])
          alist_1[[j]] + (wip[[p]] * tab[j])
      }
    } else {
      # complete case
      pij_1[mp$case.idx[[p]], ] <-
        y1w_c[mp$case.idx[[p]], ] %*% sigma_w_inv
      wip[[p]] <- sigma_w_inv
      for (j in j1_idx) {
        alist_1[[j]] <-
          # ALIST[[j]] + (sigma.w.inv[both.idx, both.idx] * TAB[j])
          alist_1[[j]] + (sigma_w_inv * tab[j])
      }
    }
  } # p

  pj <- rowsum.default(pij_1[, , drop = FALSE],
    cluster_idx,
    reorder = FALSE, na.rm = TRUE
  )

  # per cluster
  both_diag_idx <- lav_mat_diag_idx(length(both_idx))
  for (j in seq_len(nclusters)) {
    a_j_full <- alist_1[[j]]
    a_j <- a_j_full[both_idx, both_idx, drop = FALSE]
    p_j <- as.matrix(pj[j, ])
    pb_j <- as.matrix(pj[j, both_idx]) # only both.idx part
    if (nz > 0L) {
      sigma_b_z <- sigma_b_z_1[[zpat2j[j]]]
    } else {
      sigma_b_z <- sigma_b_1
    }

    ibza_j <- sigma_b_z %*% a_j
    ibza_j[both_diag_idx] <- ibza_j[both_diag_idx] + 1

    ibza_j_inv_bz <- solve.default(ibza_j, sigma_b_z)
    ibza_j_inv_bz_p <- ibza_j_inv_bz %*% pb_j
    a_ibza_j_inv_bz <- a_j %*% ibza_j_inv_bz
    a_ibza_j_inv_bz_p <- a_ibza_j_inv_bz %*% pb_j

    ibza_j_inv <- solve.default(ibza_j)
    a_ibza_j_inv <- a_j %*% ibza_j_inv
    p_ibza_j_inv <- t(crossprod(pb_j, ibza_j_inv))

    # only if we have between-only variables
    if (nz > 0L) {
      g_j <- as.matrix(gj[j, ])
      zij <- as.matrix(gz[j, ])
      zizy <- zizy_1[[zpat2j[j]]]
      zip <- zip_1[[zpat2j[j]]]

      ibza_j_inv_zizy <- solve.default(ibza_j, t(zizy))
      ibza_j_inv_g <- ibza_j_inv %*% g_j
      # ibza_j_inv_p <- ibza_j_inv %*% pb_j
      a_ibza_j_inv_g <- a_j %*% ibza_j_inv_g
      a_ibza_j_inv_zizy <- a_j %*% ibza_j_inv_zizy
      zizy_a_ibza_j_inv_g <- zizy %*% a_ibza_j_inv_g
      p_ibza_j_inv_zizy <- crossprod(pb_j, ibza_j_inv_zizy)
      ggbzpp <- 2 * a_ibza_j_inv_g + a_ibza_j_inv_bz_p - pb_j
      zijzizyp <- (2 * zij - zizy %*% pb_j)

      ###########
      # dx.mu.z #
      ###########
      tmp <- 2 * (t(p_ibza_j_inv_zizy) - zij - zizy_a_ibza_j_inv_g)
      g_muz[j, ] <- drop(tmp)

      ###############
      # dx.sigma.zz #
      ###############
      tmp1 <- (zip + zizy %*% a_ibza_j_inv_zizy # logdet
        - tcrossprod(zij) # ZA
        - tcrossprod(zizy_a_ibza_j_inv_g)) # ZB-1

      d <- (t((2 * zizy_a_ibza_j_inv_g + zizy %*% a_ibza_j_inv_bz_p)
      %*% p_ibza_j_inv_zizy)
      + zijzizyp %*% p_ibza_j_inv_zizy
        - 2 * tcrossprod(zizy_a_ibza_j_inv_g, zij))
      tmp2 <- (d + t(d)) / 2
      tmp <- tmp1 + tmp2
      # symmetry correction
      zz <- 2 * tmp
      diag(zz) <- diag(tmp)
      g_sigma_zz[j, ] <- lav_mat_vech(zz)

      ###############
      # dx.sigma.yz #
      ###############
      t0 <- -2 * a_ibza_j_inv_zizy

      t1 <- (-2 * tcrossprod(p_ibza_j_inv, g_j)
        - 1 * tcrossprod(p_ibza_j_inv, sigma_b_z %*% pb_j)
        + 2 * tcrossprod(a_ibza_j_inv_g, g_j)) %*% a_ibza_j_inv_zizy
      t2 <- -ggbzpp %*% p_ibza_j_inv_zizy
      t3 <- -tcrossprod(p_ibza_j_inv, zijzizyp)
      t4 <- 2 * tcrossprod(a_ibza_j_inv_g, zij)
      tmp <- t0 + t1 + t2 + t3 + t4
      tmp2 <- matrix(0, nrow(out$sigma.yz), ncol(out$sigma.yz))
      tmp2[both_idx, ] <- tmp
      g_sigma_yz[j, ] <- lav_mat_vec(tmp2)

      ##############
      # dx.sigma.b #
      ##############
      c <- tcrossprod(ggbzpp, p_ibza_j_inv)
      tmp <- t(a_ibza_j_inv) - tcrossprod(a_ibza_j_inv_g) + (c + t(c)) / 2
      # symmetry correction
      zz <- 2 * tmp
      diag(zz) <- diag(tmp)
      zz2 <- matrix(0, nrow(out$sigma.b), ncol(out$sigma.b))
      zz2[both_idx, both_idx] <- zz
      g_sigma_b[j, ] <- lav_mat_vech(zz2)
      # dx.sigma.b[both.idx, both.idx] <-
      #  dx.sigma.b[both.idx, both.idx, drop = FALSE] + ZZ

      # for dx.sigma.w
      part1_b <- -1 * (ibza_j_inv_g %*%
        (2 * t(ibza_j_inv_bz_p) + t(g_j) -
          t(g_j) %*% a_ibza_j_inv_bz)
        + ibza_j_inv_bz + tcrossprod(ibza_j_inv_bz_p))
      part2_b <- 2 * (ibza_j_inv_g + ibza_j_inv_bz_p) # vector
    } else {
      ##############
      # dx.sigma.b #
      ##############
      bzpp <- a_ibza_j_inv_bz_p - pb_j
      c <- tcrossprod(bzpp, p_ibza_j_inv)
      tmp <- t(a_ibza_j_inv) + (c + t(c)) / 2
      # symmetry correction
      zz <- 2 * tmp
      diag(zz) <- diag(tmp)
      zz2 <- matrix(0, nrow(out$sigma.b), ncol(out$sigma.b))
      zz2[both_idx, both_idx] <- zz
      g_sigma_b[j, ] <- lav_mat_vech(zz2)
      # dx.sigma.b[both.idx, both.idx] <-
      #  dx.sigma.b[both.idx, both.idx, drop = FALSE] + ZZ

      part1_b <- -1 * (ibza_j_inv_bz + tcrossprod(ibza_j_inv_bz_p))
      part2_b <- 2 * ibza_j_inv_bz_p # vector
    }

    ##############
    # dx.sigma.w #
    ##############

    part1 <- matrix(0, ny, ny)
    part1[both_idx, both_idx] <- part1_b

    part2 <- matrix(0, ny, 1L)
    part2[both_idx, 1L] <- part2_b

    ij_index <- which(cluster_idx == j)
    pij <- pij_1[ij_index, , drop = FALSE]

    which_compl <- which(mpi[ij_index] == 0L)
    which_incompl <- which(mpi[ij_index] != 0L)

    ap2 <- rep(list(sigma_w_inv %*% part2), length(ij_index))
    ap1a_a <- ap1a_b <- matrix(0, ny, ny)
    # A.j.full <- matrix(0, ny, ny)
    if (length(which_compl) > 0L) {
      tmp <- (sigma_w_inv %*% part1 %*% sigma_w_inv)
      ap1a_a <- tmp * length(which_compl)
      # A.j.full <- A.j.full + sigma.w.inv * length(which.compl)
    }
    if (length(which_incompl) > 0L) {
      p_idx <- mpi[ij_index][which_incompl]
      tmp <- lapply(wip[p_idx], function(x) {
        x %*% part1 %*% x
      })
      ap1a_b <- Reduce("+", tmp)
      ap2[which_incompl] <-
        lapply(wip[p_idx], function(x) {
          x %*% part2
        })
      # A.j.full <- A.j.full + Reduce("+", WIP[ p.idx ])
    }
    t1 <- ap1a_a + ap1a_b
    t2 <- (do.call("cbind", ap2) - t(pij)) %*% pij

    aa_wj <- t1 + t2

    tmp <- a_j_full + (aa_wj + t(aa_wj)) / 2
    # symmetry correction
    zz <- 2 * tmp
    diag(zz) <- diag(tmp)
    g_sigma_w[j, ] <- lav_mat_vech(zz)
    # dx.sigma.w <- dx.sigma.w + ZZ

    ###########
    # dx.mu.y #
    ###########
    tmp <- numeric(ny)
    if (nz > 0L) {
      tmp[both_idx] <- ibza_j_inv_g + ibza_j_inv_bz_p
    } else {
      tmp[both_idx] <- ibza_j_inv_bz_p
    }
    gbzpp <- a_j_full %*% tmp - p_j
    # dx.mu.y <- dx.mu.y + drop(2 * gbzpp)
    g_muy[j, ] <- drop(2 * gbzpp)
  } # j

  # browser()

  # rearrange columns to Mu.W, Mu.B, Sigma.W, Sigma.B
  ov_idx <- lp$ov.idx
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))

  # Mu.W (for within-only)
  mu_w_tilde <- matrix(0, nclusters, p_tilde)
  mu_w_tilde[, ov_idx[[1]]] <- g_muy
  mu_w_tilde[, lp$both.idx[[2]]] <- 0 # ZERO!!!
  mu_w <- mu_w_tilde[, ov_idx[[1]], drop = FALSE]

  # Mu.B
  mu_b_tilde <- matrix(0, nclusters, p_tilde)
  mu_b_tilde[, ov_idx[[1]]] <- g_muy
  if (length(between_idx) > 0L) {
    mu_b_tilde[, between_idx] <- g_muz
  }
  mu_b <- mu_b_tilde[, ov_idx[[2]], drop = FALSE]

  # Sigma.W
  sigma_w <- g_sigma_w

  # Sigma.B
  if (length(between_idx) > 0L) {
    p_tilde_star <- p_tilde * (p_tilde + 1) / 2
    b_tilde <- lav_mat_vech_rev(seq_len(p_tilde_star))

    sigma_b_tilde <- matrix(0, nclusters, p_tilde_star)

    col_idx <- lav_mat_vech(b_tilde[ov_idx[[1]], ov_idx[[1]],
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_b

    col_idx <- lav_mat_vec(b_tilde[ov_idx[[1]], between_idx,
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_yz

    col_idx <- lav_mat_vech(b_tilde[between_idx, between_idx,
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_zz

    col_idx <- lav_mat_vech(b_tilde[ov_idx[[2]], ov_idx[[2]],
      drop = FALSE
    ])
    sigma_b <- sigma_b_tilde[, col_idx, drop = FALSE]
  } else {
    p_tilde_star <- p_tilde * (p_tilde + 1) / 2
    b_tilde <- lav_mat_vech_rev(seq_len(p_tilde_star))

    sigma_b_tilde <- matrix(0, nclusters, p_tilde_star)

    col_idx <- lav_mat_vech(b_tilde[ov_idx[[1]], ov_idx[[1]],
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_b

    col_idx <- lav_mat_vech(b_tilde[ov_idx[[2]], ov_idx[[2]],
      drop = FALSE
    ])
    sigma_b <- sigma_b_tilde[, col_idx, drop = FALSE]
    # Sigma.B <- G.Sigma.b
  }

  scores <- cbind(mu_w, sigma_w, mu_b, sigma_b)

  scores
}

# first-order information: outer crossprod of scores per cluster
lav_mvn_cl_mi_info_firstorder <- function(
    y1 = NULL,
    y2 = NULL,
    lp = NULL,
    mp = NULL,
    mu_w = NULL,
    sigma_w = NULL,
    mu_b = NULL,
    sigma_b = NULL,
    x_idx = NULL,
    divide_by_two = FALSE,
    sinv_method = "eigen") {

  scores <- lav_mvn_cl_mi_sc_2l(
    y1 = y1,
    y2 = y2,
    lp = lp,
    mp = mp,
    mu_w = mu_w,
    sigma_w = sigma_w,
    mu_b = mu_b,
    sigma_b = sigma_b,
    sinv_method = sinv_method
  )

  # divide by 2 (if we want scores wrt objective function)
  if (divide_by_two) {
    scores <- scores / 2
  }

  # unit information
  information <- crossprod(scores) / lp$nclusters[[2]]


  # if x.idx, set rows/cols to zero
  if (length(x_idx) > 0L) {
    nw <- length(as.vector(mu_w))
    nw_star <- nw * (nw + 1) / 2
    nb <- length(as.vector(mu_b))
    ov_idx <- lp$ov.idx

    x_idx_w <- which(ov_idx[[1]] %in% x_idx)
    if (length(x_idx_w) > 0L) {
      xw_idx <- c(
        x_idx_w,
        nw + lav_mat_vech_which_idx(n = nw, idx = x_idx_w)
      )
    } else {
      xw_idx <- integer(0L)
    }
    x_idx_b <- which(ov_idx[[2]] %in% x_idx)
    if (length(x_idx_b) > 0L) {
      xb_idx <- c(
        x_idx_b,
        nb + lav_mat_vech_which_idx(n = nb, idx = x_idx_b)
      )
    } else {
      xb_idx <- integer(0L)
    }

    all_idx <- c(xw_idx, nw + nw_star + xb_idx)

    information[all_idx, ] <- 0
    information[, all_idx] <- 0
  }

  information
}

# observed information
# order: mu.w within, vech(sigma.w) within, mu.b between, vech(sigma.b) between
# mu.w rows/cols that are splitted within/between are forced to zero
#
# numerical approximation (for now)
lav_mvn_cl_mi_info_observed <- function(
    y1 = NULL,
    y2 = NULL,
    lp = NULL,
    mp = NULL,
    ylp = NULL,
    mu_w = NULL,
    sigma_w = NULL,
    mu_b = NULL,
    sigma_b = NULL,
    x_idx = integer(0L),
    sinv_method = "eigen") {

  nw <- length(as.vector(mu_w))
  nw_star <- nw * (nw + 1) / 2
  nb <- length(as.vector(mu_b))
  nb_star <- nb * (nb + 1) / 2

  ov_idx <- lp$ov.idx
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))

  # Mu.W (for within-only)
  mu_w_tilde <- numeric(p_tilde)
  mu_w_tilde[ov_idx[[1]]] <- mu_w

  # local function -- gradient
  grad <- function(x) {
    # Mu.W (for within-only)
    mu_w_tilde2 <- numeric(p_tilde)
    mu_w_tilde2[ov_idx[[1]]] <- x[1:nw]
    mu_w_tilde2[lp$both.idx[[2]]] <- mu_w_tilde[lp$both.idx[[2]]]
    mu_w2 <- mu_w_tilde2[ov_idx[[1]]]

    sigma_w2 <- lav_mat_vech_rev(x[nw + 1:nw_star])
    mu_b2 <- x[nw + nw_star + 1:nb]
    sigma_b2 <- lav_mat_vech_rev(x[nw + nw_star + nb + 1:nb_star])

    dx <- lav_mvn_cl_mi_dlogl_2l_samp(
      y1 = y1, y2 = y2, lp = lp, mp = mp,
      mu_w = mu_w2, sigma_w = sigma_w2,
      mu_b = mu_b2, sigma_b = sigma_b2,
      return_list = FALSE,
      sinv_method = sinv_method
    )

    # dx is for -2*logl
    -1 / 2 * dx
  }

  # start.x
  start_x <- c(
    as.vector(mu_w), lav_mat_vech(sigma_w),
    as.vector(mu_b), lav_mat_vech(sigma_b)
  )

  # total information
  information <- -1 * numDeriv::jacobian(func = grad, x = start_x)

  # unit information
  information <- information / lp$nclusters[[2]]


  # if x.idx, set rows/cols to zero
  if (length(x_idx) > 0L) {
    x_idx_w <- which(ov_idx[[1]] %in% x_idx)
    if (length(x_idx_w) > 0L) {
      xw_idx <- c(
        x_idx_w,
        nw + lav_mat_vech_which_idx(n = nw, idx = x_idx_w)
      )
    } else {
      xw_idx <- integer(0L)
    }
    x_idx_b <- which(ov_idx[[2]] %in% x_idx)
    if (length(x_idx_b) > 0L) {
      xb_idx <- c(
        x_idx_b,
        nb + lav_mat_vech_which_idx(n = nb, idx = x_idx_b)
      )
    } else {
      xb_idx <- integer(0L)
    }

    all_idx <- c(xw_idx, nw + nw_star + xb_idx)

    information[all_idx, ] <- 0
    information[, all_idx] <- 0
  }

  information
}
