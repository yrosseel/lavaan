# loglikelihood clustered/twolevel data

# YR: first version around Feb 2017
# YR 28 Oct 2024: [EM steps:] if fx.delta is NA, check if we have (near)-perfect
#                 correlations, and provide an informative warning


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
lav_mvnorm_cluster_implied22l <- function(lp = NULL,
                                          implied = NULL,
                                          mu_w = NULL,
                                          mu_b = NULL,
                                          sigma_w = NULL,
                                          sigma_b = NULL) {
  if (!is.null(implied)) {
    # FIXME: only for single-group analysis!
    sigma_w <- implied$cov[[1]]
    mu_w <- implied$mean[[1]]

    sigma_b <- implied$cov[[2]]
    mu_b <- implied$mean[[2]]
  }

  # within/between.idx
  between_idx <- lp$between.idx[[2]]
  within_idx <- lp$within.idx[[2]]
  both_idx <- lp$both.idx[[2]]

  # ov.idx per level
  ov_idx <- lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))

  # Sigma.W.tilde
  sigma_w_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_w_tilde[ov_idx[[1]], ov_idx[[1]]] <- sigma_w

  # Sigma.B.tilde
  sigma_b_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_b_tilde[ov_idx[[2]], ov_idx[[2]]] <- sigma_b

  # Mu.W.tilde
  mu_w_tilde <- numeric(p_tilde)
  mu_w_tilde[ov_idx[[1]]] <- mu_w

  # Mu.B.tilde
  mu_b_tilde <- numeric(p_tilde)
  mu_b_tilde[ov_idx[[2]]] <- mu_b

  # add Mu.W[within.idx] to Mu.B
  mu_wb_tilde <- numeric(p_tilde)
  mu_wb_tilde[within_idx] <- mu_w_tilde[within_idx]
  mu_wb_tilde[both_idx] <- (mu_b_tilde[both_idx] +
    mu_w_tilde[both_idx])

  # set Mu.W[both.idx] to zero (after we added to WB)
  mu_w_tilde[both_idx] <- 0
  # get Mu.B[both.idx[ from WB
  mu_b_tilde[both_idx] <- mu_wb_tilde[both_idx]

  # map to matrices needed for loglik
  if (length(within_idx) > 0L) {
    mu_b_tilde[within_idx] <- 0
  }
  if (length(between_idx) > 0L) {
    mu_z <- mu_b_tilde[between_idx]
    mu_y <- mu_wb_tilde[-between_idx]
    mu_w_1 <- mu_w_tilde[-between_idx]
    mu_b_1 <- mu_b_tilde[-between_idx]
    sigma_zz <- sigma_b_tilde[between_idx, between_idx, drop = FALSE]
    sigma_yz <- sigma_b_tilde[-between_idx, between_idx, drop = FALSE]
    sigma_b_1 <- sigma_b_tilde[-between_idx, -between_idx, drop = FALSE]
    sigma_w_1 <- sigma_w_tilde[-between_idx, -between_idx, drop = FALSE]
  } else {
    mu_z <- numeric(0L)
    mu_y <- mu_wb_tilde
    mu_w_1 <- mu_w_tilde
    mu_b_1 <- mu_b_tilde
    sigma_zz <- matrix(0, 0L, 0L)
    sigma_yz <- matrix(0, nrow(sigma_b_tilde), 0L)
    sigma_b_1 <- sigma_b_tilde
    sigma_w_1 <- sigma_w_tilde
  }

  list(
    sigma.w = sigma_w_1, sigma.b = sigma_b_1, sigma.zz = sigma_zz,
    sigma.yz = sigma_yz, mu.z = mu_z, mu.y = mu_y, mu.w = mu_w_1,
    mu.b = mu_b_1
  )
}

lav_mvnorm_cluster_2l2implied <- function(lp,
                                          sigma_w = NULL,
                                          sigma_b = NULL,
                                          sigma_zz = NULL,
                                          sigma_yz = NULL,
                                          mu_z = NULL,
                                          mu_y = NULL,
                                          mu_w = NULL,
                                          mu_b = NULL) {
  # between.idx
  between_idx <- lp$between.idx[[2]]
  within_idx <- lp$within.idx[[2]]
  # both_idx <- lp$both.idx[[2]]

  # ov.idx per level
  ov_idx <- lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))

  # if we have mu.y, convert to mu.w and mu.b
  if (!is.null(mu_y)) {
    mu_b <- mu_y
    mu_w_tilde <- numeric(p_tilde)
    mu_w_tilde[ov_idx[[1]]] <- mu_y

    # NO NEED TO SET THIS TO ZERO!
    # otherwise, we get non-symmetric Hessian!! 0.6-5

    # if(length(within.idx) > 0L) {
    #    mu.w.tilde[  -within.idx ] <- 0
    # } else {
    #    mu.w.tilde[] <- 0
    # }
    mu_w <- mu_w_tilde[ov_idx[[1]]]
  }

  # new in 0.6-18: ensure mu.w[both.idx] is zero?
  # NO: we get Hessian is not fully symmetric again!!
  # only do this at the very end (post-estimation)

  # Mu.W.tilde <- numeric(p.tilde)
  # Mu.B.tilde <- numeric(p.tilde)
  # Mu.W.tilde[ov.idx[[1]]] <- mu.w
  # Mu.B.tilde[ov.idx[[2]]] <- mu.b
  # Mu.B.tilde[between.idx] <- mu.z
  # if (length(within.idx) > 0) {
  #   Mu.B.tilde[within.idx] <- 0
  # }
  # Mu.B.tilde[both.idx] <- Mu.W.tilde[both.idx] + Mu.B.tilde[both.idx]
  # Mu.W.tilde[both.idx] <- 0
  # Mu.W <- Mu.W.tilde[ov.idx[[1]]]
  # Mu.B <- Mu.B.tilde[ov.idx[[2]]]

  mu_w_tilde_1 <- numeric(p_tilde)
  ###### DEBUG ##############
  #if (length(within.idx) > 0) {
      mu_w_tilde_1[ov_idx[[1]]] <- mu_w
  #}
  ###########################
  mu_w_1 <- mu_w_tilde_1[ov_idx[[1]]]

  # Mu.B
  mu_b_tilde <- numeric(p_tilde)
  mu_b_tilde[ov_idx[[1]]] <- mu_b
  mu_b_tilde[between_idx] <- mu_z
  if (length(within_idx) > 0) {
      mu_b_tilde[within_idx] <- 0
  }
  mu_b_1 <- mu_b_tilde[ov_idx[[2]]]


  # Sigma.W
  sigma_w_1 <- sigma_w

  # Sigma.B
  sigma_b_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_b_tilde[ov_idx[[1]], ov_idx[[1]]] <- sigma_b
  sigma_b_tilde[ov_idx[[1]], between_idx] <- sigma_yz
  sigma_b_tilde[between_idx, ov_idx[[1]]] <- t(sigma_yz)
  sigma_b_tilde[between_idx, between_idx] <- sigma_zz
  sigma_b_1 <- sigma_b_tilde[ov_idx[[2]], ov_idx[[2]], drop = FALSE]

  list(Mu.W = mu_w_1, Mu.B = mu_b_1, Sigma.W = sigma_w_1, Sigma.B = sigma_b_1)
}


# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics
# (not yet reordered)
lav_mvnorm_cluster_loglik_samplestats_2l <- function(ylp = NULL, # nolint
                                                     lp = NULL,
                                                     mu_w = NULL,
                                                     sigma_w = NULL,
                                                     mu_b = NULL,
                                                     sigma_b = NULL,
                                                     sinv_method = "eigen",
                                                     log2pi = FALSE,
                                                     minus_two = TRUE) {
  # map implied to 2l matrices
  out <- lav_mvnorm_cluster_implied22l(
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
  # nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]
  between_idx <- lp$between.idx[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  cluster_size_ns <- lp$cluster.size.ns[[2]]

  # Y1 samplestats
  if (length(between_idx) > 0L) {
    s_pw <- ylp[[2]]$Sigma.W[-between_idx, -between_idx, drop = FALSE]
  } else {
    s_pw <- ylp[[2]]$Sigma.W
  }

  # Y2 samplestats
  cov_d <- ylp[[2]]$cov.d
  mean_d <- ylp[[2]]$mean.d

  # common parts:
  sigma_w_inv <- lav_matrix_symmetric_inverse(
    s = sigma_w_1,
    logdet = TRUE, sinv_method = sinv_method
  )
  sigma_w_logdet <- attr(sigma_w_inv, "logdet")
  attr(sigma_w_inv, "logdet") <- NULL

  if (length(between_idx) > 0L) {
    sigma_zz_inv <- lav_matrix_symmetric_inverse(
      s = sigma_zz,
      logdet = TRUE, sinv_method = sinv_method
    )
    sigma_zz_logdet <- attr(sigma_zz_inv, "logdet")
    attr(sigma_zz_inv, "logdet") <- NULL
    sigma_yz_zi <- sigma_yz %*% sigma_zz_inv
    sigma_zi_zy <- t(sigma_yz_zi)
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy
  } else {
    sigma_zz_logdet <- 0
    sigma_b_z <- sigma_b_1
  }

  # min 2* logliklihood
  l <- numeric(ncluster_sizes) # logdet
  m_b <- numeric(ncluster_sizes) # between qf
  for (clz in seq_len(ncluster_sizes)) {
    # cluster size
    nj <- cluster_sizes[clz]

    # data between
    y2yc <- (cov_d[[clz]] + tcrossprod(mean_d[[clz]] - c(mu_z, mu_y)))

    # FIXME: avoid reorder/b.idx, so we can use between.idx
    if (length(between_idx) > 0L) {
      b_idx <- seq_along(lp$between.idx[[2]])
      y2yc_zz <- y2yc[b_idx, b_idx, drop = FALSE]
      y2yc_yz <- y2yc[-b_idx, b_idx, drop = FALSE]
      y2yc_yy <- y2yc[-b_idx, -b_idx, drop = FALSE]
    } else {
      y2yc_yy <- y2yc
    }

    # construct sigma.j
    sigma_j <- (nj * sigma_b_z) + sigma_w_1
    sigma_j_inv <- lav_matrix_symmetric_inverse(
      s = sigma_j,
      logdet = TRUE, sinv_method = sinv_method
    )
    sigma_j_logdet <- attr(sigma_j_inv, "logdet")
    attr(sigma_j_inv, "logdet") <- NULL

    # check: what if sigma.j is non-pd? should not happen
    if (is.na(sigma_j_logdet)) {
      # stop, and return NA right away
      # return(as.numeric(NA))
      # FORCE?
      # sigma.j <- lav_matrix_symmetric_force_pd(sigma.j)
      # sigma.j.inv <- lav_matrix_symmetric_inverse(s = sigma.j,
      #           logdet = TRUE, sinv_method = Sinv.method)
      # sigma.j.logdet <- attr(sigma.j.inv, "logdet")
      # attr(sigma.j.inv, "logdet") <- NULL
    }

    # logdet -- between only
    l[clz] <- (sigma_zz_logdet + sigma_j_logdet)

    if (length(between_idx) > 0L) {
      # part 1 -- zz
      sigma_ji_yz_zi <- sigma_j_inv %*% sigma_yz_zi
      vinv_11 <- sigma_zz_inv + nj * (sigma_zi_zy %*% sigma_ji_yz_zi)
      q_zz <- sum(vinv_11 * y2yc_zz)

      # part 2 -- yz
      q_yz <- -nj * sum(sigma_ji_yz_zi * y2yc_yz)
    } else {
      q_zz <- q_yz <- 0
    }

    # part 5 -- yyc
    q_yyc <- -nj * sum(sigma_j_inv * y2yc_yy)

    # qf -- between only
    m_b[clz] <- q_zz + 2 * q_yz - q_yyc
  }
  # q.yya + q.yyb
  # the reason why we multiply the trace by 'N - nclusters' is
  # S.PW has been divided by 'N - nclusters'
  q_w <- sum(cluster_size - 1) * sum(sigma_w_inv * s_pw)
  # logdet within part
  l_w <- sum(cluster_size - 1) * sigma_w_logdet

  # -2*times logl (without the constant)
  loglik <- sum(l * cluster_size_ns) + sum(m_b * cluster_size_ns) + q_w + l_w

  # functions below compute -2 * logl
  if (!minus_two) {
    loglik <- loglik / (-2)
  }

  # constant
  # Note: total 'N' = (nobs * #within vars) + (nclusters * #between vars)
  if (log2pi) {
    log_2pi <- log(2 * pi)
    n_within <- length(c(lp$both.idx[[2]], lp$within.idx[[2]]))
    n_between <- length(lp$between.idx[[2]])
    p <- lp$nclusters[[1]] * n_within + lp$nclusters[[2]] * n_between
    constant <- -(p * log_2pi) / 2
    loglik <- loglik + constant
  }

  # loglik.x (only if loglik is requested)
  if (length(unlist(lp$ov.x.idx)) > 0L && log2pi && !minus_two) {
    loglik <- loglik - ylp[[2]]$loglik.x
  }

  loglik
}


# first derivative -2*logl wrt Mu.W, Mu.B, Sigma.W, Sigma.B
lav_mvnorm_cluster_dlogl_2l_samplestats <- function(ylp = NULL, # nolint
                                                    lp = NULL,
                                                    mu_w = NULL,
                                                    sigma_w = NULL,
                                                    mu_b = NULL,
                                                    sigma_b = NULL,
                                                    return_list = FALSE,
                                                    sinv_method = "eigen") {
  # map implied to 2l matrices
  out <- lav_mvnorm_cluster_implied22l(
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
  # cluster_size <- lp$cluster.size[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  # cluster_idx <- lp$cluster.idx[[2]]
  between_idx <- lp$between.idx[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  cluster_size_ns <- lp$cluster.size.ns[[2]]

  # Y1
  if (length(between_idx) > 0L) {
    s_pw <- ylp[[2]]$Sigma.W[-between_idx, -between_idx, drop = FALSE]
  } else {
    s_pw <- ylp[[2]]$Sigma.W
  }

  # Y2
  cov_d <- ylp[[2]]$cov.d
  mean_d <- ylp[[2]]$mean.d

  # common parts:
  sigma_w_inv <- lav_matrix_symmetric_inverse(
    s = sigma_w_1,
    logdet = FALSE, sinv_method = sinv_method
  )

  # weighted level-1 and level-2 derivative accumulators
  d_mu_y <- numeric(length(mu_y))
  d_sigma_w1_vech <- numeric(length(lav_matrix_vech(sigma_w_1)))
  d_sigma_b_vech <- numeric(length(lav_matrix_vech(sigma_b_1)))

  if (length(between_idx) > 0L) {
    d_mu_z <- numeric(length(mu_z))
    d_sigma_zz_vech <- numeric(length(lav_matrix_vech(sigma_zz)))
    d_sigma_yz_vec <- numeric(length(lav_matrix_vec(sigma_yz)))

    sigma_zz_inv <- lav_matrix_symmetric_inverse(
      s = sigma_zz,
      logdet = FALSE, sinv_method = sinv_method
    )
    sigma_yz_zi <- sigma_yz %*% sigma_zz_inv
    sigma_zi_zy <- t(sigma_yz_zi)
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy

    b_idx <- seq_along(lp$between.idx[[2]])

    for (clz in seq_len(ncluster_sizes)) {
      # cluster size
      nj <- cluster_sizes[clz]
      cluster_weight <- cluster_size_ns[clz]

      # level-2 vectors
      zyc <- mean_d[[clz]] - c(mu_z, mu_y)
      yc <- zyc[-b_idx]
      zc <- zyc[b_idx]

      # level-2 crossproducts
      y2yc <- (cov_d[[clz]] + tcrossprod(mean_d[[clz]] - c(mu_z, mu_y)))
      y2yc_zz <- y2yc[b_idx, b_idx, drop = FALSE]
      y2yc_yz <- y2yc[-b_idx, b_idx, drop = FALSE]
      y2yc_yy <- y2yc[-b_idx, -b_idx, drop = FALSE]

      # construct sigma.j
      sigma_j <- (nj * sigma_b_z) + sigma_w_1
      sigma_j_inv <- lav_matrix_symmetric_inverse(
        s = sigma_j,
        logdet = FALSE, sinv_method = sinv_method
      )
      sigma_ji_yz_zi <- sigma_j_inv %*% sigma_yz_zi
      sigma_zi_zy_ji <- t(sigma_ji_yz_zi)

      # common parts
      j_yzj <- nj * (sigma_j_inv %*%
        (sigma_yz_zi %*% y2yc_zz %*% t(sigma_yz_zi)
          - y2yc_yz %*% t(sigma_yz_zi)
          - t(y2yc_yz %*% t(sigma_yz_zi)) + y2yc_yy)
        %*% sigma_j_inv)

      z1 <- y2yc_zz %*% t(sigma_ji_yz_zi) %*% sigma_yz
      yz1 <- t(y2yc_yz) %*% sigma_j_inv %*% sigma_yz


      # Mu.Z
      d_mu_z <- d_mu_z + cluster_weight * (-2 * as.numeric(
        (sigma_zz_inv + nj * (sigma_zi_zy_ji %*% sigma_yz_zi)) %*% zc
          - nj * sigma_zi_zy_ji %*% yc
      ))

      # MU.Y
      d_mu_y <- d_mu_y + cluster_weight * (
        2 * nj * as.numeric(zc %*% sigma_zi_zy_ji - yc %*% sigma_j_inv)
      )

      # SIGMA.W (between part)
      g_sigma_w <- sigma_j_inv - j_yzj
      tmp <- g_sigma_w * 2
      diag(tmp) <- diag(g_sigma_w)
      d_sigma_w1_vech <- d_sigma_w1_vech +
        cluster_weight * lav_matrix_vech(tmp)

      # SIGMA.B
      g_sigma_b <- nj * (sigma_j_inv - j_yzj)
      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      d_sigma_b_vech <- d_sigma_b_vech +
        cluster_weight * lav_matrix_vech(tmp)

      # SIGMA.ZZ
      g_sigma_zz <- (sigma_zz_inv + nj * sigma_zz_inv %*% (
        t(sigma_yz) %*% (sigma_j_inv - j_yzj) %*% sigma_yz
          - (1 / nj * y2yc_zz + t(z1) + z1 - t(yz1) - yz1)) %*%
        sigma_zz_inv)

      tmp <- g_sigma_zz * 2
      diag(tmp) <- diag(g_sigma_zz)
      d_sigma_zz_vech <- d_sigma_zz_vech +
        cluster_weight * lav_matrix_vech(tmp)

      # SIGMA.ZY
      g_sigma_yz <- 2 * nj * (
        (sigma_j_inv %*%
          (sigma_yz_zi %*% y2yc_zz - sigma_yz - y2yc_yz)
          + j_yzj %*% sigma_yz) %*% sigma_zz_inv)

      d_sigma_yz_vec <- d_sigma_yz_vec +
        cluster_weight * lav_matrix_vec(g_sigma_yz)
    }

    # level-1
    d_sigma_w1 <- lav_matrix_vech_reverse(d_sigma_w1_vech)
    d_sigma_b <- lav_matrix_vech_reverse(d_sigma_b_vech)

    # level-2
    d_sigma_zz <- lav_matrix_vech_reverse(d_sigma_zz_vech)
    d_sigma_yz <- matrix(
      d_sigma_yz_vec,
      nrow(sigma_yz), ncol(sigma_yz)
    )
  # between.idx
  } else { # no level-2 variables

    for (clz in seq_len(ncluster_sizes)) {
      # cluster size
      nj <- cluster_sizes[clz]
      cluster_weight <- cluster_size_ns[clz]

      # level-2 vectors
      yc <- mean_d[[clz]] - mu_y

      # level-2 crossproducts
      y2yc_yy <- (cov_d[[clz]] + tcrossprod(mean_d[[clz]] - mu_y))

      # construct sigma.j
      sigma_j <- (nj * sigma_b_1) + sigma_w_1
      sigma_j_inv <- lav_matrix_symmetric_inverse(
        s = sigma_j,
        logdet = FALSE, sinv_method = sinv_method
      )
      # common part
      j_yyj <- nj * sigma_j_inv %*% y2yc_yy %*% sigma_j_inv

      # MU.Y
      d_mu_y <- d_mu_y +
        cluster_weight * (-2 * nj * as.numeric(yc %*% sigma_j_inv))

      # SIGMA.W (between part)
      g_sigma_w <- sigma_j_inv - j_yyj
      tmp <- g_sigma_w * 2
      diag(tmp) <- diag(g_sigma_w)
      d_sigma_w1_vech <- d_sigma_w1_vech +
        cluster_weight * lav_matrix_vech(tmp)

      # SIGMA.B
      g_sigma_b <- nj * (sigma_j_inv - j_yyj)
      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      d_sigma_b_vech <- d_sigma_b_vech +
        cluster_weight * lav_matrix_vech(tmp)
    }

    # level-1
    d_sigma_w1 <- lav_matrix_vech_reverse(d_sigma_w1_vech)
    d_sigma_b <- lav_matrix_vech_reverse(d_sigma_b_vech)
    # level-2
    d_mu_z <- numeric(0L)
    d_sigma_zz <- matrix(0, 0L, 0L)
    d_sigma_yz <- matrix(0, 0L, 0L)
  }

  # Sigma.W (bis)
  d_sigma_w2 <- (lp$nclusters[[1]] - nclusters) * (sigma_w_inv
  - sigma_w_inv %*% s_pw %*% sigma_w_inv)
  tmp <- d_sigma_w2 * 2
  diag(tmp) <- diag(d_sigma_w2)
  d_sigma_w2 <- tmp

  d_sigma_w <- d_sigma_w1 + d_sigma_w2

  # rearrange
  dout <- lav_mvnorm_cluster_2l2implied(
    lp = lp,
    sigma_w = d_sigma_w, sigma_b = d_sigma_b,
    sigma_yz = d_sigma_yz, sigma_zz = d_sigma_zz,
    mu_y = d_mu_y, mu_z = d_mu_z
  )

  if (return_list) {
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
lav_mvnorm_cluster_scores_2l <- function(y1 = NULL,
                                         ylp = NULL,
                                         lp = NULL,
                                         mu_w = NULL,
                                         sigma_w = NULL,
                                         mu_b = NULL,
                                         sigma_b = NULL,
                                         sinv_method = "eigen") {
  # map implied to 2l matrices
  out <- lav_mvnorm_cluster_implied22l(
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
  cluster_size <- lp$cluster.size[[2]]
  cluster_idx <- lp$cluster.idx[[2]]
  between_idx <- lp$between.idx[[2]]

  # Y1
  if (length(between_idx) > 0L) {
    y1w <- y1[, -lp$between.idx[[2]], drop = FALSE]
  } else {
    y1w <- y1
  }
  y1w_cm <- t(t(y1w) - mu_y)

  # Y2
  y2 <- ylp[[2]]$Y2
  # NOTE: ORDER mu.b must match Y2
  mu_b_1 <- numeric(ncol(y2))
  if (length(between_idx) > 0L) {
    mu_b_1[-lp$between.idx[[2]]] <- mu_y
    mu_b_1[lp$between.idx[[2]]] <- mu_z
  } else {
    mu_b_1 <- mu_y
  }
  y2_cm <- t(t(y2) - mu_b_1)

  # common parts:
  sigma_w_inv <- lav_matrix_symmetric_inverse(
    s = sigma_w_1,
    logdet = FALSE, sinv_method = sinv_method
  )

  # both level-1 and level-2
  g_muy <- matrix(0, nclusters, length(mu_y))
  g_sigma_w_1 <- matrix(0, nclusters, length(lav_matrix_vech(sigma_w_1)))
  g_sigma_b_1 <- matrix(0, nclusters, length(lav_matrix_vech(sigma_b_1)))
  g_muz <- matrix(0, nclusters, length(mu_z))
  g_sigma_zz_1 <- matrix(0, nclusters, length(lav_matrix_vech(sigma_zz)))
  g_sigma_yz_1 <- matrix(0, nclusters, length(lav_matrix_vec(sigma_yz)))

  if (length(between_idx) > 0L) {
    sigma_zz_inv <- lav_matrix_symmetric_inverse(
      s = sigma_zz,
      logdet = FALSE, sinv_method = sinv_method
    )
    sigma_yz_zi <- sigma_yz %*% sigma_zz_inv
    sigma_zi_zy <- t(sigma_yz_zi)
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy


    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster_size[cl]

      # data within for the cluster (centered by mu.y)
      y1m <- y1w_cm[cluster_idx == cl, , drop = FALSE]
      yc <- y2_cm[cl, -lp$between.idx[[2]]]
      zc <- y2_cm[cl, lp$between.idx[[2]]]

      # data between
      y2yc <- tcrossprod(y2_cm[cl, ])
      y2yc_zz <- y2yc[lp$between.idx[[2]],
        lp$between.idx[[2]],
        drop = FALSE
      ]
      y2yc_yz <- y2yc[-lp$between.idx[[2]],
        lp$between.idx[[2]],
        drop = FALSE
      ]
      y2yc_yy <- y2yc[-lp$between.idx[[2]],
        -lp$between.idx[[2]],
        drop = FALSE
      ]

      # construct sigma.j
      sigma_j <- (nj * sigma_b_z) + sigma_w_1
      sigma_j_inv <- lav_matrix_symmetric_inverse(
        s = sigma_j,
        logdet = FALSE, sinv_method = sinv_method
      )
      sigma_ji_yz_zi <- sigma_j_inv %*% sigma_yz_zi
      sigma_zi_zy_ji <- t(sigma_ji_yz_zi)

      # common parts
      j_yzj <- nj * (sigma_j_inv %*%
        (sigma_yz_zi %*% y2yc_zz %*% t(sigma_yz_zi)
          - y2yc_yz %*% t(sigma_yz_zi)
          - t(y2yc_yz %*% t(sigma_yz_zi)) + y2yc_yy)
        %*% sigma_j_inv)

      z1 <- y2yc_zz %*% t(sigma_ji_yz_zi) %*% sigma_yz
      yz1 <- t(y2yc_yz) %*% sigma_j_inv %*% sigma_yz


      # Mu.Z
      g_muz[cl, ] <- -2 * as.numeric(
        (sigma_zz_inv + nj * (sigma_zi_zy_ji %*% sigma_yz_zi)) %*% zc
          - nj * sigma_zi_zy_ji %*% yc
      )

      # MU.Y
      g_muy[cl, ] <- 2 * nj * as.numeric(zc %*% sigma_zi_zy_ji -
        yc %*% sigma_j_inv)

      # SIGMA.W
      g_sigma_w <- ((nj - 1) * sigma_w_inv
        - sigma_w_inv %*% (crossprod(y1m) - nj * y2yc_yy) %*% sigma_w_inv
        + sigma_j_inv - j_yzj)

      tmp <- g_sigma_w * 2
      diag(tmp) <- diag(g_sigma_w)
      g_sigma_w_1[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g_sigma_b <- nj * (sigma_j_inv - j_yzj)

      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      g_sigma_b_1[cl, ] <- lav_matrix_vech(tmp)


      # SIGMA.ZZ
      g_sigma_zz <- (sigma_zz_inv + nj * sigma_zz_inv %*% (
        t(sigma_yz) %*% (sigma_j_inv - j_yzj) %*% sigma_yz
          - (1 / nj * y2yc_zz + t(z1) + z1 - t(yz1) - yz1)) %*%
        sigma_zz_inv)

      tmp <- g_sigma_zz * 2
      diag(tmp) <- diag(g_sigma_zz)
      g_sigma_zz_1[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.ZY
      g_sigma_yz <- 2 * nj * (
        (sigma_j_inv %*%
          (sigma_yz_zi %*% y2yc_zz - sigma_yz - y2yc_yz)
          + j_yzj %*% sigma_yz) %*% sigma_zz_inv)

      g_sigma_yz_1[cl, ] <- lav_matrix_vec(g_sigma_yz)
    }
  # between.idx
  } else { # no level-2 variables

    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster_size[cl]

      # data within for the cluster (centered by mu.y)
      y1m <- y1w_cm[cluster_idx == cl, , drop = FALSE]
      yc <- y2_cm[cl, ]

      # data between
      y2yc_yy <- tcrossprod(y2_cm[cl, ])

      # construct sigma.j
      sigma_j <- (nj * sigma_b_1) + sigma_w_1
      sigma_j_inv <- lav_matrix_symmetric_inverse(
        s = sigma_j,
        logdet = FALSE, sinv_method = sinv_method
      )
      # common part
      j_yyj <- nj * sigma_j_inv %*% y2yc_yy %*% sigma_j_inv

      # MU.Y
      g_muy[cl, ] <- -2 * nj * as.numeric(yc %*% sigma_j_inv)

      # SIGMA.W
      g_sigma_w <- ((nj - 1) * sigma_w_inv
        - sigma_w_inv %*% (crossprod(y1m) - nj * y2yc_yy) %*% sigma_w_inv
        + sigma_j_inv - j_yyj)
      tmp <- g_sigma_w * 2
      diag(tmp) <- diag(g_sigma_w)
      g_sigma_w_1[cl, ] <- lav_matrix_vech(tmp)

      # SIGMA.B
      g_sigma_b <- nj * (sigma_j_inv - j_yyj)
      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      g_sigma_b_1[cl, ] <- lav_matrix_vech(tmp)
    }
  }

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
  sigma_w <- g_sigma_w_1

  # Sigma.B
  if (length(between_idx) > 0L) {
    p_tilde_star <- p_tilde * (p_tilde + 1) / 2
    b_tilde <- lav_matrix_vech_reverse(seq_len(p_tilde_star))

    sigma_b_tilde <- matrix(0, nclusters, p_tilde_star)

    col_idx <- lav_matrix_vech(b_tilde[ov_idx[[1]], ov_idx[[1]],
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_b_1

    col_idx <- lav_matrix_vec(b_tilde[ov_idx[[1]], between_idx,
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_yz_1

    col_idx <- lav_matrix_vech(b_tilde[between_idx, between_idx,
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_zz_1

    col_idx <- lav_matrix_vech(b_tilde[ov_idx[[2]], ov_idx[[2]],
      drop = FALSE
    ])
    sigma_b <- sigma_b_tilde[, col_idx, drop = FALSE]
  } else {
    p_tilde_star <- p_tilde * (p_tilde + 1) / 2
    b_tilde <- lav_matrix_vech_reverse(seq_len(p_tilde_star))

    sigma_b_tilde <- matrix(0, nclusters, p_tilde_star)

    col_idx <- lav_matrix_vech(b_tilde[ov_idx[[1]], ov_idx[[1]],
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_b_1

    col_idx <- lav_matrix_vech(b_tilde[ov_idx[[2]], ov_idx[[2]],
      drop = FALSE
    ])
    sigma_b <- sigma_b_tilde[, col_idx, drop = FALSE]
    # Sigma.B <- G.Sigma.b
  }

  scores <- cbind(mu_w, sigma_w, mu_b, sigma_b)

  scores
}


# first-order information: outer crossprod of scores per cluster
lav_mvnorm_cluster_information_firstorder <- function(y1 = NULL,  # nolint
                                                      ylp = NULL,
                                                      lp = NULL,
                                                      mu_w = NULL,
                                                      sigma_w = NULL,
                                                      mu_b = NULL,
                                                      sigma_b = NULL,
                                                      x_idx = NULL,
                                                      divide_by_two = FALSE,
                                                      sinv_method = "eigen") {
  # n <- NROW(y1)

  scores <- lav_mvnorm_cluster_scores_2l(
    y1 = y1,
    ylp = ylp,
    lp = lp,
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
        nw + lav_matrix_vech_which_idx(n = nw, idx = x_idx_w)
      )
    } else {
      xw_idx <- integer(0L)
    }
    x_idx_b <- which(ov_idx[[2]] %in% x_idx)
    if (length(x_idx_b) > 0L) {
      xb_idx <- c(
        x_idx_b,
        nb + lav_matrix_vech_which_idx(n = nb, idx = x_idx_b)
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

# expected information 'h1' model
# order: mu.w within, vech(sigma.w) within, mu.b between, vech(sigma.b) between
# mu.w rows/cols that are splitted within/between are forced to zero
lav_mvnorm_cluster_information_expected <- function(lp = NULL,  # nolint
                                                    mu_w = NULL,
                                                    sigma_w = NULL,
                                                    mu_b = NULL,
                                                    sigma_b = NULL,
                                                    x_idx = integer(0L),
                                                    sinv_method = "eigen") {
  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    lp = lp,
    mu_w = mu_w, mu_b = mu_b,
    sigma_w = sigma_w, sigma_b = sigma_b
  )
  # mu_y <- out$mu.y
  # mu_z <- out$mu.z
  sigma_w_1 <- out$sigma.w
  sigma_b_1 <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz

  # create Delta.W.tilde, Delta.B.tilde
  ov_idx <- lp$ov.idx
  nw <- length(ov_idx[[1]])
  nb <- length(ov_idx[[2]])
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))
  p_tilde_star <- p_tilde * (p_tilde + 1) / 2
  npar <- p_tilde + p_tilde_star
  b_tilde <- lav_matrix_vech_reverse(seq_len(p_tilde_star))
  w_idx <- lav_matrix_vech(b_tilde[ov_idx[[1]], ov_idx[[1]], drop = FALSE])
  b_idx <- lav_matrix_vech(b_tilde[ov_idx[[2]], ov_idx[[2]], drop = FALSE])

  delta_w_tilde <- matrix(0, npar, npar)
  delta_b_tilde <- matrix(0, npar, npar)
  delta_w_tilde[
    c(ov_idx[[1]], w_idx + p_tilde),
    c(ov_idx[[1]], w_idx + p_tilde)
  ] <- diag(nw + nw * (nw + 1) / 2)
  delta_b_tilde[
    c(ov_idx[[2]], b_idx + p_tilde),
    c(ov_idx[[2]], b_idx + p_tilde)
  ] <- diag(nb + nb * (nb + 1) / 2)
  delta_w_tilde <- cbind(delta_w_tilde, matrix(0, npar, npar))
  delta_b_tilde <- cbind(matrix(0, npar, npar), delta_b_tilde)

  nobs <- lp$nclusters[[1]]
  nclusters <- lp$nclusters[[2]]
  # cluster_size <- lp$cluster.size[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  n_s <- lp$cluster.size.ns[[2]]
  between_idx <- lp$between.idx[[2]]

  information_j <- matrix(0, npar * 2, npar * 2)
  for (clz in seq_len(ncluster_sizes)) {
    # cluster size
    nj <- cluster_sizes[clz]

    # Delta.j -- changes per cluster(size)
    # this is why we can not write info = t(delta) info.sat delta
    delta_j <- delta_b_tilde + 1 / nj * delta_w_tilde

    # compute Sigma.j
    sigma_j <- sigma_w_1 + nj * sigma_b_1
    if (length(between_idx) > 0L) {
      omega_j <- matrix(0, p_tilde, p_tilde)
      omega_j[-between_idx, -between_idx] <- 1 / nj * sigma_j
      omega_j[-between_idx, between_idx] <- sigma_yz
      omega_j[between_idx, -between_idx] <- t(sigma_yz)
      omega_j[between_idx, between_idx] <- sigma_zz
      # omega.j <- rbind( cbind(sigma.zz, t(sigma.yz)),
      #                  cbind(sigma.yz, 1/nj * sigma.j) )
    } else {
      omega_j <- 1 / nj * sigma_j
    }
    omega_j_inv <- solve(omega_j)

    i11_j <- omega_j_inv
    # if (lav_use_lavaanC()) {
    #   I22.j <-
    #     lavaanC::m_kronecker_dup_pre_post(omega.j.inv, multiplicator = 0.5)
    # } else {
      i22_j <- 0.5 *
        lav_matrix_duplication_pre_post(omega_j_inv %x% omega_j_inv)
    # }
    i_j <- lav_matrix_bdiag(i11_j, i22_j)
    info_j <- t(delta_j) %*% i_j %*% delta_j

    information_j <- information_j + n_s[clz] * info_j
  }

  sigma_w_inv <- lav_matrix_symmetric_inverse(
    s = sigma_w, logdet = FALSE,
    sinv_method = sinv_method
  )
  # create Sigma.W.inv.tilde
  sigma_w_inv_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_w_inv_tilde[ov_idx[[1]], ov_idx[[1]]] <- sigma_w_inv

  i11_w <- sigma_w_inv_tilde
  # if (lav_use_lavaanC()) {
  #   I22.W <- lavaanC::m_kronecker_dup_pre_post(Sigma.W.inv.tilde,
  #                                                    multiplicator = 0.5)
  # } else {
    i22_w <- 0.5 *
      lav_matrix_duplication_pre_post(sigma_w_inv_tilde %x% sigma_w_inv_tilde)
  # }
  i_w <- lav_matrix_bdiag(i11_w, i22_w)
  information_w <- (nobs - nclusters) *
    (t(delta_w_tilde) %*% i_w %*% delta_w_tilde)

  # unit information
  information_tilde <- 1 / lp$nclusters[[2]] * (information_w + information_j)

  # force zero for means both.idx in within part
  information_tilde[lp$both.idx[[2]], ] <- 0
  information_tilde[, lp$both.idx[[2]]] <- 0

  # if x.idx, set rows/cols to zero
  if (length(x_idx) > 0L) {
    xw_idx <- c(
      x_idx,
      p_tilde + lav_matrix_vech_which_idx(n = p_tilde, idx = x_idx)
    )
    xb_idx <- npar + xw_idx
    all_idx <- c(xw_idx, xb_idx)
    information_tilde[all_idx, ] <- 0
    information_tilde[, all_idx] <- 0
  }

  # remove redundant rows/cols
  ok_idx <- c(
    ov_idx[[1]],
    w_idx + p_tilde,
    npar + ov_idx[[2]],
    npar + b_idx + p_tilde
  )

  information <- information_tilde[ok_idx, ok_idx]

  information
}


# expected information -- delta
# for non-saturated models only
lav_mvnorm_cluster_information_expected_delta <- function(lp = NULL, # nolint
                                                      delta = NULL,
                                                      mu_w = NULL,
                                                      sigma_w = NULL,
                                                      mu_b = NULL,
                                                      sigma_b = NULL,
                                                      sinv_method = "eigen") {
  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    lp = lp,
    mu_w = mu_w, mu_b = mu_b,
    sigma_w = sigma_w, sigma_b = sigma_b
  )
  # mu_y <- out$mu.y
  # mu_z <- out$mu.z
  sigma_w_1 <- out$sigma.w
  sigma_b_1 <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz

  # Delta -- this group
  npar <- NCOL(delta)

  # create Delta.W.tilde, Delta.B.tilde
  ov_idx <- lp$ov.idx
  nw <- length(ov_idx[[1]])
  nw_star <- nw * (nw + 1) / 2
  nb <- length(ov_idx[[2]])

  delta_w <- delta[1:(nw + nw_star), , drop = FALSE]
  delta_b <- delta[-(1:(nw + nw_star)), , drop = FALSE]

  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))
  p_tilde_star <- p_tilde * (p_tilde + 1) / 2
  delta_w_tilde_mu <- matrix(0, p_tilde, npar)
  delta_w_tilde_sigma <- matrix(0, p_tilde_star, npar)
  delta_b_tilde_mu <- matrix(0, p_tilde, npar)
  delta_b_tilde_sigma <- matrix(0, p_tilde_star, npar)

  delta_w_tilde_mu[ov_idx[[1]], ] <- delta_w[1:nw, ]
  delta_b_tilde_mu[ov_idx[[2]], ] <- delta_b[1:nb, ]

  # correct Delta to reflect Mu.W[ both.idx ] is added to Mu.B[ both.idx ]
  # changed in 0.6-5
  delta_b_tilde_mu[lp$both.idx[[2]], ] <-
    (delta_b_tilde_mu[lp$both.idx[[2]], ] +
      delta_w_tilde_mu[lp$both.idx[[2]], ])
  delta_w_tilde_mu[lp$both.idx[[2]], ] <- 0


  b_tilde <- lav_matrix_vech_reverse(seq_len(p_tilde_star))
  w_idx <- lav_matrix_vech(b_tilde[ov_idx[[1]], ov_idx[[1]], drop = FALSE])
  b_idx <- lav_matrix_vech(b_tilde[ov_idx[[2]], ov_idx[[2]], drop = FALSE])
  delta_w_tilde_sigma[w_idx, ] <- delta_w[-(1:nw), ]
  delta_b_tilde_sigma[b_idx, ] <- delta_b[-(1:nb), ]

  delta_w_tilde <- rbind(delta_w_tilde_mu, delta_w_tilde_sigma)
  delta_b_tilde <- rbind(delta_b_tilde_mu, delta_b_tilde_sigma)

  nobs <- lp$nclusters[[1]]
  nclusters <- lp$nclusters[[2]]
  # cluster_size <- lp$cluster.size[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  n_s <- lp$cluster.size.ns[[2]]
  between_idx <- lp$between.idx[[2]]

  information_j <- matrix(0, npar, npar)
  for (clz in seq_len(ncluster_sizes)) {
    # cluster size
    nj <- cluster_sizes[clz]

    # Delta.j -- changes per cluster(size)
    # this is why we can not write info = t(delta) info.sat delta
    delta_j <- delta_b_tilde + 1 / nj * delta_w_tilde

    # compute Sigma.j
    sigma_j <- sigma_w_1 + nj * sigma_b_1
    if (length(between_idx) > 0L) {
      omega_j <- matrix(0, p_tilde, p_tilde)
      omega_j[-between_idx, -between_idx] <- 1 / nj * sigma_j
      omega_j[-between_idx, between_idx] <- sigma_yz
      omega_j[between_idx, -between_idx] <- t(sigma_yz)
      omega_j[between_idx, between_idx] <- sigma_zz
      # omega.j <- rbind( cbind(sigma.zz, t(sigma.yz)),
      #                  cbind(sigma.yz, 1/nj * sigma.j) )
    } else {
      omega_j <- 1 / nj * sigma_j
    }
    omega_j_inv <- solve(omega_j)

    i11_j <- omega_j_inv
    # if (lav_use_lavaanC()) {
    #   I22.j <-
    #      lavaanC::m_kronecker_dup_pre_post(omega.j.inv, multiplicator = 0.5)
    # } else {
      i22_j <- 0.5 *
             lav_matrix_duplication_pre_post(omega_j_inv %x% omega_j_inv)
    # }
    i_j <- lav_matrix_bdiag(i11_j, i22_j)
    info_j <- t(delta_j) %*% i_j %*% delta_j

    information_j <- information_j + n_s[clz] * info_j
  }


  sigma_w_inv <- lav_matrix_symmetric_inverse(
    s = sigma_w_1, logdet = FALSE,
    sinv_method = sinv_method
  )
  i11_w <- sigma_w_inv
  # if (lav_use_lavaanC()) {
  #   I22.w <- lavaanC::m_kronecker_dup_pre_post(Sigma.W.inv,
  #                                                       multiplicator = 0.5)
  # } else {
    i22_w <- 0.5 * lav_matrix_duplication_pre_post(sigma_w_inv %x% sigma_w_inv)
  # }
  i_w <- lav_matrix_bdiag(i11_w, i22_w)

  # force zero for means both.idx in within part
  # changed in 0.6-5
  i_w[lp$both.idx[[2]], ] <- 0
  i_w[, lp$both.idx[[2]]] <- 0

  information_w <- (nobs - nclusters) * (t(delta_w) %*% i_w %*% delta_w)

  # unit information
  information <- 1 / lp$nclusters[[2]] * (information_w + information_j)


  information
}


# observed information
# order: mu.w within, vech(sigma.w) within, mu.b between, vech(sigma.b) between
# mu.w rows/cols that are splitted within/between are forced to zero
#
# numerical approximation (for now)
lav_mvnorm_cluster_information_observed <- function(lp = NULL,     # nolint
                                                    ylp = NULL,
                                                    mu_w = NULL,
                                                    sigma_w = NULL,
                                                    mu_b = NULL,
                                                    sigma_b = NULL,
                                                    x_idx = integer(0L),
                                                    sinv_method = "eigen") {
  # nobs <- lp$nclusters[[1]]

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

    sigma_w2 <- lav_matrix_vech_reverse(x[nw + 1:nw_star])
    mu_b2 <- x[nw + nw_star + 1:nb]
    sigma_b2 <- lav_matrix_vech_reverse(x[nw + nw_star + nb + 1:nb_star])

    dx <- lav_mvnorm_cluster_dlogl_2l_samplestats(
      ylp = ylp,
      lp = lp, mu_w = mu_w2, sigma_w = sigma_w2,
      mu_b = mu_b2, sigma_b = sigma_b2,
      return_list = FALSE,
      sinv_method = sinv_method
    )

    # dx is for -2*logl
    -1 / 2 * dx
  }

  # start.x
  start_x <- c(
    as.vector(mu_w), lav_matrix_vech(sigma_w),
    as.vector(mu_b), lav_matrix_vech(sigma_b)
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
        nw + lav_matrix_vech_which_idx(n = nw, idx = x_idx_w)
      )
    } else {
      xw_idx <- integer(0L)
    }
    x_idx_b <- which(ov_idx[[2]] %in% x_idx)
    if (length(x_idx_b) > 0L) {
      xb_idx <- c(
        x_idx_b,
        nb + lav_matrix_vech_which_idx(n = nb, idx = x_idx_b)
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

# estimate ML estimates of Mu.W, Mu.B, Sigma.W, Sigma.B
# using the EM algorithm
#
# per cluster-SIZE
#
lav_mvnorm_cluster_em_sat <- function(ylp = NULL,
                                      lp = NULL,
                                      tol = 1e-04,
                                      max_iter = 5000,
                                      min_variance = 1e-05) {
  # lavdata
  between_idx <- lp$between.idx[[2]]
  # within_idx <- lp$within.idx[[2]]
  y2 <- ylp[[2]]$Y2

  # starting values for Sigma
  ov_idx <- lp$ov.idx
  # COVT <- lavsamplestats@cov[[1]]
  # Sigma.W <- diag( diag(COVT)[ov.idx[[1]]] )
  # Sigma.B <- diag( diag(COVT)[ov.idx[[2]]] )
  sigma_w_1 <- diag(length(ov_idx[[1]]))
  sigma_b_1 <- diag(length(ov_idx[[2]]))
  mu_w_1 <- numeric(length(ov_idx[[1]]))
  mu_b_1 <- numeric(length(ov_idx[[2]]))
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
    ylp = ylp, lp = lp,
    mu_w = mu_w_1, sigma_w = sigma_w_1,
    mu_b = mu_b_1, sigma_b = sigma_b_1,
    sinv_method = "eigen", log2pi = TRUE, minus_two = FALSE
  )

  # if verbose, report
  if (lav_verbose()) {
    cat(
      "EM iter:", sprintf("%3d", 0),
      " fx =", sprintf("%17.10f", fx),
      "\n"
    )
  }

  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    lp = lp,
    mu_w = mu_w_1, sigma_w = sigma_w_1, mu_b = mu_b_1, sigma_b = sigma_b_1
  )
  # mu_y <- out$mu.y
  mu_z <- out$mu.z
  mu_w <- out$mu.w
  mu_b <- out$mu.b
  sigma_w <- out$sigma.w
  sigma_b <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz

  # mu.z and sigma.zz can be computed beforehand
  if (length(between_idx) > 0L) {
    z <- y2[, between_idx, drop = FALSE]
    mu_z <- colMeans(z, na.rm = TRUE)
    sigma_zz <- cov(z, use = "pairwise.complete.obs") *
                              (lp$nclusters[[2]] - 1L) / lp$nclusters[[2]]
    # sigma.zz <- 1/Lp$nclusters[[2]] * crossprod(Z) - tcrossprod(mu.z)
    # Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop=FALSE]
  }

  # EM iterations
  fx_old <- fx
  for (i in 1:max_iter) {
    # E-step
    estep <- lav_mvnorm_cluster_em_estepb( # Y1 = Y1,
      ylp = ylp,
      lp = lp,
      sigma_w = sigma_w,
      sigma_b = sigma_b,
      mu_w = mu_w,
      mu_b = mu_b,
      sigma_yz = sigma_yz,
      sigma_zz = sigma_zz,
      mu_z = mu_z
    )

    # mstep
    sigma_w <- estep$sigma.w
    sigma_b <- estep$sigma.b
    sigma_yz <- estep$sigma.yz
    mu_w <- estep$mu.w
    mu_b <- estep$mu.b

    implied2 <- lav_mvnorm_cluster_2l2implied(
      lp = lp,
      sigma_w = estep$sigma.w, sigma_b = estep$sigma.b,
      sigma_zz = sigma_zz, sigma_yz = estep$sigma.yz,
      mu_z = mu_z,
      mu_y = NULL, mu_w = estep$mu.w, mu_b = estep$mu.b
    )

    # check for (near-zero) variances at the within level, and set
    # them to min.variance
    sigma_w_1 <- implied2$Sigma.W
    zero_var <- which(diag(sigma_w_1) < min_variance)
    if (length(zero_var) > 0L) {
      sigma_w_1[, zero_var] <- sigma_w[, zero_var] <- 0
      sigma_w_1[zero_var, ] <- sigma_w[zero_var, ] <- 0
      diag(sigma_w_1)[zero_var] <- diag(sigma_w)[zero_var] <- min_variance
    }

    fx <- lav_mvnorm_cluster_loglik_samplestats_2l(
      ylp = ylp,
      lp = lp, mu_w = implied2$Mu.W, sigma_w = sigma_w_1,
      mu_b = implied2$Mu.B, sigma_b = implied2$Sigma.B,
      sinv_method = "eigen", log2pi = TRUE, minus_two = FALSE
    )

    # fx.delta
    fx_delta <- fx - fx_old

    # check if fx.delta is finite
    if (!is.finite(fx_delta)) {
      # not good ... something is very wrong; perhaps near-singular
      # matrices?
      cat("\n")
      cat("FATAL problem: dumping Sigma.W and Sigma.B matrices:\n\n")
      cat("Sigma.W:\n")
      print(sigma_w_1)
      cat("\n")
      cat("Sigma.B:\n")
      print(implied2$Sigma.B)
      cat("\n")
      lav_msg_stop(gettext(
        "EM steps of the saturated (H1) model failed; some matrices may
         be singular; please check your data for (near-)perfect correlations."))
    }

    # what if fx.delta is negative?
    if (fx_delta < 0) {
      lav_msg_warn(gettext(
        "logl decreased during EM steps of the saturated (H1) model"))
    }

    if (lav_verbose()) {
      cat(
        "EM iter:", sprintf("%3d", i),
        " fx =", sprintf("%17.10f", fx),
        " fx.delta =", sprintf("%9.8f", fx_delta),
        "\n"
      )
    }

    # convergence check
    if (fx_delta < tol) {
      break
    } else {
      fx_old <- fx
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
                                     verbose_x = FALSE,
                                     fx_tol = 1e-08,
                                     dx_tol = 1e-05,
                                     max_iter = 5000,
                                     mstep_iter_max = 10000L,
                                     mstep_rel_tol = 1e-10) {
  # single group only for now
  stopifnot(lavdata@ngroups == 1L)

  # lavdata
  lp <- lavdata@Lp[[1]] # first group only (for now)
  ov_names_l <- lavdata@ov.names.l[[1]] # first group only (for now)
  # y1 <- lavdata@X[[1]] # first group only
  ylp <- lavsamplestats@YLp[[1]] # first group only

  between_idx <- lp$between.idx[[2]]
  y2 <- ylp[[2]]$Y2

  # initial values
  x_current <- lav_model_get_parameters(lavmodel)

  # implied
  if (is.null(lavimplied)) {
    lavimplied <- lav_model_implied(lavmodel)
  }

  # TODO: what if current 'starting' parameters imply a non-pd sigma.b?

  # report initial fx
  fx <- lav_mvnorm_cluster_loglik_samplestats_2l(
    ylp = ylp, lp = lp,
    mu_w = lavimplied$mean[[1]], sigma_w = lavimplied$cov[[1]],
    mu_b = lavimplied$mean[[2]], sigma_b = lavimplied$cov[[2]],
    sinv_method = "eigen", log2pi = TRUE, minus_two = FALSE
  )

  # if verbose, report
  if (lav_verbose()) {
    cat(
      "EM iter:", sprintf("%3d", 0),
      " fx =", sprintf("%17.10f", fx),
      "\n"
    )
  }

  # translate to internal matrices
  out <- lav_mvnorm_cluster_implied22l(
    lp = lp,
    mu_w = lavimplied$mean[[1]], sigma_w = lavimplied$cov[[1]],
    mu_b = lavimplied$mean[[2]], sigma_b = lavimplied$cov[[2]]
  )
  # mu_y <- out$mu.y
  mu_z <- out$mu.z
  mu_w <- out$mu.w
  mu_b <- out$mu.b
  sigma_w <- out$sigma.w
  sigma_b <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz

  # mu.z and sigma.zz can be computed beforehand
  if (length(between_idx) > 0L) {
    z <- y2[, between_idx, drop = FALSE]
    mu_z <- colMeans(y2)[between_idx]
    sigma_zz <- cov(z) * (lp$nclusters[[2]] - 1L) / lp$nclusters[[2]]
    # sigma.zz <- 1/Lp$nclusters[[2]] * crossprod(Z) - tcrossprod(mu.z)
    # Y1Y1 <- Y1Y1[-between.idx, -between.idx, drop=FALSE]
  }

  # EM iterations
  fx_old <- fx
  # fx2_old <- 0
  # rel <- numeric(max_iter)
  for (i in 1:max_iter) {
    # E-step
    estep <- lav_mvnorm_cluster_em_estepb(
      ylp = ylp,
      lp = lp,
      sigma_w = sigma_w,
      sigma_b = sigma_b,
      mu_w = mu_w,
      mu_b = mu_b,
      sigma_yz = sigma_yz,
      sigma_zz = sigma_zz,
      mu_z = mu_z
    )

    # back to model-implied dimensions
    implied <- lav_mvnorm_cluster_2l2implied(
      lp = lp,
      sigma_w = estep$sigma.w, sigma_b = estep$sigma.b,
      sigma_zz = sigma_zz, sigma_yz = estep$sigma.yz,
      mu_z = mu_z,
      mu_y = NULL, mu_w = estep$mu.w, mu_b = estep$mu.b
    )
    rownames(implied$Sigma.W) <- ov_names_l[[1]]
    rownames(implied$Sigma.B) <- ov_names_l[[2]]

    # M-step

    # fit two-group model
    local_partable <- lavpartable
    # if a group column exists, delete it (it will be overriden anyway)
    local_partable$group <- NULL
    level_idx <- which(names(local_partable) == "level")
    names(local_partable)[level_idx] <- "group"
    local_partable$est <- NULL
    local_partable$se <- NULL

    # give current values as starting values
    free_idx <- which(lavpartable$free > 0L)
    local_partable$ustart[free_idx] <- x_current

    local_fit <- lavaan(local_partable,
      sample.cov = list(
        within = implied$Sigma.W,
        between = implied$Sigma.B
      ),
      sample.mean = list(
        within = implied$Mu.W,
        between = implied$Mu.B
      ),
      sample.nobs = lp$nclusters,
      sample.cov.rescale = FALSE,
      control = list(
        iter.max = mstep_iter_max,
        rel.tol = mstep_rel_tol
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

    implied2 <- local_fit@implied
    fx <- lav_mvnorm_cluster_loglik_samplestats_2l(
      ylp = ylp,
      lp = lp, mu_w = implied2$mean[[1]], sigma_w = implied2$cov[[1]],
      mu_b = implied2$mean[[2]], sigma_b = implied2$cov[[2]],
      sinv_method = "eigen", log2pi = TRUE, minus_two = FALSE
    )

    # fx.delta
    fx_delta <- fx - fx_old

    # derivatives
    lavmodel <- lav_model_set_parameters(lavmodel, x = local_fit@optim$x)
    dx <- lav_model_gradient(lavmodel,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats
    )
    max_dx <- max(abs(dx))

    if (lav_verbose()) {
      cat(
        "EM iter:", sprintf("%3d", i),
        " fx =", sprintf("%17.10f", fx),
        " fx.delta =", sprintf("%9.8f", fx_delta),
        " mstep.iter =", sprintf(
          "%3d",
          lavInspect(local_fit, "iterations")
        ),
        " max.dx = ", sprintf("%9.8f", max_dx),
        "\n"
      )
    }

    # stopping rule check
    if (fx_delta < fx_tol) {
      if (lav_verbose()) {
        cat("EM stopping rule reached: fx.delta < ", fx_tol, "\n")
      }
      break
    } else {
      fx_old <- fx
      x_current <- local_fit@optim$x
      if (verbose_x) {
        print(round(x_current, 3))
      }
    }

    # second stopping rule check -- derivatives
    if (max_dx < dx_tol) {
      if (lav_verbose()) {
        cat("EM stopping rule reached: max.dx < ", dx_tol, "\n")
      }
      break
    }

    # translate to internal matrices
    out <- lav_mvnorm_cluster_implied22l(
      lp = lp,
      mu_w = implied2$mean[[1]], sigma_w = implied2$cov[[1]],
      mu_b = implied2$mean[[2]], sigma_b = implied2$cov[[2]]
    )
    # mu_y <- out$mu.y
    mu_z <- out$mu.z
    mu_w <- out$mu.w
    mu_b <- out$mu.b
    sigma_w <- out$sigma.w
    sigma_b <- out$sigma.b
    sigma_zz <- out$sigma.zz
    sigma_yz <- out$sigma.yz
  } # EM iterations

  x <- local_fit@optim$x

  # add attributes
  if (i < max_iter) {
    attr(x, "converged") <- TRUE
    attr(x, "warn.txt") <- ""
  } else {
    attr(x, "converged") <- FALSE
    attr(x, "warn.txt") <- paste("maxmimum number of iterations (",
      max_iter, ") ",
      "was reached without convergence.\n",
      sep = ""
    )
  }
  attr(x, "iterations") <- i
  attr(x, "control") <- list(
    em.iter.max = max_iter,
    em.fx.tol = fx_tol,
    em.dx.tol = dx_tol
  )
  attr(fx, "fx.group") <- fx # single group for now
  attr(x, "fx") <- fx

  x
}

# get the random effects (here: expected values for cluster means)
# and optionally a standard error
lav_mvnorm_cluster_em_estep_ranef <- function(ylp = NULL, # nolint
                                              lp = NULL,
                                              sigma_w = NULL,
                                              sigma_b = NULL,
                                              sigma_yz = NULL,
                                              sigma_zz = NULL,
                                              mu_z = NULL,
                                              mu_w = NULL,
                                              mu_b = NULL,
                                              se = FALSE) {
  # sample stats
  # nobs <- lp$nclusters[[1]]
  nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]
  between_idx <- lp$between.idx[[2]]

  y2 <- ylp[[2]]$Y2

  nvar_y <- ncol(sigma_w)
  # nvar_z <- ncol(sigma_zz)

  mb_j <- matrix(0, nrow = nclusters, ncol = nvar_y)
  se_j <- matrix(0, nrow = nclusters, ncol = nvar_y)

  mu_y <- mu_w + mu_b

  if (length(between_idx) > 0L) {
    sigma_1 <- cbind(sigma_yz, sigma_b)
    mu <- c(mu_z, mu_y)
  } else {
    sigma_1 <- sigma_b
    mu <- mu_y
  }

  # E-step
  for (cl in seq_len(nclusters)) {
    nj <- cluster_size[cl]

    # data
    if (length(between_idx) > 0L) {
      # z comes first!
      b_j <- c(
        y2[cl, between_idx],
        y2[cl, -between_idx]
      )
      # ybar_j <- y2[cl, -between_idx]
    } else {
      # ybar_j <- b_j <- y2[cl, ]
    }

    sigma_j <- sigma_w + nj * sigma_b
    if (length(between_idx) > 0L) {
      omega_j <- rbind(
        cbind(sigma_zz, t(sigma_yz)),
        cbind(sigma_yz, 1 / nj * sigma_j)
      )
    } else {
      omega_j <- 1 / nj * sigma_j
    }
    omega_j_inv <- solve(omega_j)

    # E(v|y)
    ev <- as.numeric(mu_b + (sigma_1 %*% omega_j_inv %*% (b_j - mu)))
    mb_j[cl, ] <- ev

    if (se) {
      # Cov(v|y)
      covv <- sigma_b - (sigma_1 %*% omega_j_inv %*% t(sigma_1))

      # force symmetry
      covv <- (covv + t(covv)) / 2

      covv_diag <- diag(covv)
      nonzero_idx <- which(covv_diag > 0)

      se_j[cl, ] <- numeric(length(covv_diag))
      se_j[cl, nonzero_idx] <- sqrt(covv_diag[nonzero_idx])
    }
  }

  if (se) {
    attr(mb_j, "se") <- se_j
  }

  mb_j
}

# per cluster
lav_mvnorm_cluster_em_estep <- function( # Y1           = NULL,
                                        ylp = NULL,
                                        lp = NULL,
                                        sigma_w = NULL,
                                        sigma_b = NULL,
                                        sigma_yz = NULL,
                                        sigma_zz = NULL,
                                        mu_z = NULL,
                                        mu_w = NULL,
                                        mu_b = NULL) {
  # sample stats
  nobs <- lp$nclusters[[1]]
  nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]
  # cluster_idx <- lp$cluster.idx[[2]]
  # within_idx <- lp$within.idx[[2]]
  between_idx <- lp$between.idx[[2]]
  # both_idx <- lp$both.idx[[2]]

  y2 <- ylp[[2]]$Y2
  y1y1 <- ylp[[2]]$Y1Y1

  nvar_y <- ncol(sigma_w)
  nvar_z <- ncol(sigma_zz)

  cw2_j <- matrix(0, nrow = nvar_y, ncol = nvar_y)
  cb_j <- matrix(0, nrow = nvar_y, ncol = nvar_y)
  mw_j <- matrix(0, nrow = nclusters, ncol = nvar_y)
  mb_j <- matrix(0, nrow = nclusters, ncol = nvar_y)
  zy_j <- matrix(0, nrow = nvar_z, ncol = nvar_y)

  mu_y <- mu_w + mu_b

  if (length(between_idx) > 0L) {
    sigma_1 <- cbind(sigma_yz, sigma_b)
    mu <- c(mu_z, mu_y)
    y1y1 <- y1y1[-between_idx, -between_idx, drop = FALSE]
  } else {
    sigma_1 <- sigma_b
    mu <- mu_y
  }

  # E-step
  for (cl in seq_len(nclusters)) {
    nj <- cluster_size[cl]

    # data
    if (length(between_idx) > 0L) {
      # z comes first!
      b_j <- c(
        y2[cl, between_idx],
        y2[cl, -between_idx]
      )
      ybar_j <- y2[cl, -between_idx]
    } else {
      ybar_j <- b_j <- y2[cl, ]
    }

    sigma_j <- sigma_w + nj * sigma_b
    if (length(between_idx) > 0L) {
      omega_j <- rbind(
        cbind(sigma_zz, t(sigma_yz)),
        cbind(sigma_yz, 1 / nj * sigma_j)
      )
    } else {
      omega_j <- 1 / nj * sigma_j
    }
    omega_j_inv <- solve(omega_j)

    # E(v|y)
    ev <- as.numeric(mu_b + (sigma_1 %*% omega_j_inv %*% (b_j - mu)))

    # Cov(v|y)
    covv <- sigma_b - (sigma_1 %*% omega_j_inv %*% t(sigma_1))

    # force symmetry
    covv <- (covv + t(covv)) / 2

    # E(vv|y) = Cov(v|y) + E(v|y)E(v|y)^T
    evv <- covv + tcrossprod(ev)

    # store for this cluster
    mw_j[cl, ] <- ybar_j - ev
    mb_j[cl, ] <- ev
    cw2_j <- cw2_j + nj * (evv - tcrossprod(ybar_j, ev)
      - tcrossprod(ev, ybar_j))
    cb_j <- cb_j + evv

    # between only
    if (length(between_idx) > 0L) {
      zy_j <- zy_j + tcrossprod(y2[cl, between_idx], ev)
    }
  }

  m_w <- 1 / nobs * colSums(mw_j * cluster_size)
  m_b <- 1 / nclusters * colSums(mb_j)
  c_b <- 1 / nclusters * cb_j
  c_w <- 1 / nobs * (y1y1 + cw2_j)
  # end of E-step

  # make symmetric (not needed here?)
  # C.b <- (C.b + t(C.b))/2
  # C.w <- (C.w + t(C.w))/2

  # between only
  if (length(between_idx) > 0L) {
    a <- 1 / nclusters * zy_j - tcrossprod(mu_z, m_b)
  }

  sigma_w <- c_w - tcrossprod(m_w)
  sigma_b <- c_b - tcrossprod(m_b)
  mu_w <- m_w
  mu_b <- m_b

  if (length(between_idx) > 0L) {
    sigma_yz <- t(a)
  }

  list(
    sigma.w = sigma_w, sigma.b = sigma_b, mu.w = mu_w, mu.b = mu_b,
    sigma.yz = sigma_yz, sigma.zz = sigma_zz, mu.z = mu_z
  )
}

# per cluster SIZE
lav_mvnorm_cluster_em_estepb <- function( # Y1           = NULL, # not used!
                                         ylp = NULL,
                                         lp = NULL,
                                         sigma_w = NULL,
                                         sigma_b = NULL,
                                         sigma_yz = NULL,
                                         sigma_zz = NULL,
                                         mu_z = NULL,
                                         mu_w = NULL,
                                         mu_b = NULL) {
  # sample stats
  nobs <- lp$nclusters[[1]]
  nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]
  # cluster_idx <- lp$cluster.idx[[2]]
  between_idx <- lp$between.idx[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  n_s <- lp$cluster.size.ns[[2]]

  y2 <- ylp[[2]]$Y2
  y1y1 <- ylp[[2]]$Y1Y1

  nvar_y <- ncol(sigma_w)
  nvar_z <- ncol(sigma_zz)

  mu_y <- mu_w + mu_b

  if (length(between_idx) > 0L) {
    sigma_1 <- cbind(sigma_yz, sigma_b)
    mu <- c(mu_z, mu_y)
    y1y1 <- y1y1[-between_idx, -between_idx, drop = FALSE]
  } else {
    sigma_1 <- sigma_b
    mu <- mu_y
  }

  # per cluster SIZE
  cw2_s <- matrix(0, nrow = nvar_y, ncol = nvar_y)
  cb_s <- matrix(0, nrow = nvar_y, ncol = nvar_y)
  mw_s <- matrix(0, nrow = ncluster_sizes, ncol = nvar_y)
  mb_s <- matrix(0, nrow = ncluster_sizes, ncol = nvar_y)
  zy_s <- matrix(0, nvar_z, nvar_y)

  # E-step
  for (clz in seq_len(ncluster_sizes)) {
    # cluster size
    nj <- cluster_sizes[clz]

    # data
    if (length(between_idx) > 0L) {
      # z comes first!
      b_j <- cbind(
        y2[cluster_size == nj, between_idx, drop = FALSE],
        y2[cluster_size == nj, -between_idx, drop = FALSE]
      )
      ybar_j <- y2[cluster_size == nj, -between_idx, drop = FALSE]
    } else {
      ybar_j <- b_j <- y2[cluster_size == nj, , drop = FALSE]
    }

    sigma_j <- sigma_w + nj * sigma_b
    if (length(between_idx) > 0L) {
      omega_j <- rbind(
        cbind(sigma_zz, t(sigma_yz)),
        cbind(sigma_yz, 1 / nj * sigma_j)
      )
    } else {
      omega_j <- 1 / nj * sigma_j
    }
    omega_j_inv <- solve(omega_j)
    sigma_1_j_inv <- sigma_1 %*% omega_j_inv

    # E(v|y)
    b_jc <- t(t(b_j) - mu)
    tmp <- b_jc %*% t(sigma_1_j_inv)
    ev <- t(t(tmp) + mu_b)

    # Cov(v|y)
    covv <- n_s[clz] * (sigma_b - (sigma_1_j_inv %*% t(sigma_1)))

    # force symmetry
    covv <- (covv + t(covv)) / 2

    # E(vv|y) = Cov(v|y) + E(v|y)E(v|y)^T
    evv <- covv + crossprod(ev)

    # store for this cluster SIZE
    mw_s[clz, ] <- nj * colSums(ybar_j - ev)
    mb_s[clz, ] <- colSums(ev)
    cw2_s <- cw2_s + nj * (evv - crossprod(ybar_j, ev)
      - crossprod(ev, ybar_j))
    cb_s <- cb_s + evv

    # between only
    if (length(between_idx) > 0L) {
      zy_s <- zy_s + crossprod(y2[cluster_size == nj, between_idx,
        drop = FALSE
      ], ev)
    }
  } # cluster-sizes

  m_ws <- 1 / nobs * colSums(mw_s)
  m_bs <- 1 / nclusters * colSums(mb_s)
  c_bs <- 1 / nclusters * cb_s
  c_ws <- 1 / nobs * (y1y1 + cw2_s)

  # between only
  if (length(between_idx) > 0L) {
    as_1 <- 1 / nclusters * zy_s - tcrossprod(mu_z, m_bs)
  }

  sigma_w <- c_ws - tcrossprod(m_ws)
  sigma_b <- c_bs - tcrossprod(m_bs)
  mu_w <- m_ws
  mu_b <- m_bs
  if (length(between_idx) > 0L) {
    sigma_yz <- t(as_1)
  }

  list(
    sigma.w = sigma_w, sigma.b = sigma_b, mu.w = mu_w, mu.b = mu_b,
    sigma.yz = sigma_yz, sigma.zz = sigma_zz, mu.z = mu_z
  )
}
