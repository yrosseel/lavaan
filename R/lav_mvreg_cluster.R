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
lav_mvreg_cl_implied22l <- function(lp = NULL,
                                         implied = NULL,
                                         res_int_w = NULL,
                                         res_int_b = NULL,
                                         res_pi_w = NULL,
                                         res_pi_b = NULL,
                                         res_sigma_w = NULL,
                                         res_sigma_b = NULL) {
  if (!is.null(implied)) {
    # FIXME: only for single-group analysis!
    res_sigma_w <- implied$res.cov[[1]]
    res_int_w <- implied$res.int[[1]]
    res_pi_w <- implied$res.slopes[[1]]

    res_sigma_b <- implied$res.cov[[2]]
    res_int_b <- implied$res.int[[2]]
    res_pi_b <- implied$res.slopes[[2]]
  }

  # within/between idx
  within_x_idx <- lp$within.x.idx[[1]]
  between_y_idx <- lp$between.y.idx[[2]]
  between_x_idx <- lp$between.x.idx[[2]]

  # ov.idx per level
  ov_idx <- lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))

  # only 'y'
  ov_y_idx <- lp$ov.y.idx

  # two levels only (for now)
  ov_y_idx1 <- ov_y_idx[[1]]
  ov_y_idx2 <- ov_y_idx[[2]]

  # Sigma.W.tilde
  sigma_w_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_w_tilde[ov_y_idx1, ov_y_idx1] <- res_sigma_w

  # INT.W.tilde
  int_w_tilde <- matrix(0, p_tilde, 1L)
  int_w_tilde[ov_y_idx1, 1L] <- res_int_w

  # PI.W.tilde
  pi_w_tilde <- matrix(0, p_tilde, ncol(res_pi_w))
  pi_w_tilde[ov_y_idx1, ] <- res_pi_w

  beta_w_tilde <- rbind(t(int_w_tilde), t(pi_w_tilde))



  # Sigma.B.tilde
  sigma_b_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_b_tilde[ov_y_idx2, ov_y_idx2] <- res_sigma_b

  # INT.B.tilde
  int_b_tilde <- matrix(0, p_tilde, 1L)
  int_b_tilde[ov_y_idx2, 1L] <- res_int_b

  # PI.B.tilde
  pi_b_tilde <- matrix(0, p_tilde, ncol(res_pi_b))
  pi_b_tilde[ov_y_idx2, ] <- res_pi_b

  beta_b_tilde <- rbind(t(int_b_tilde), t(pi_b_tilde))

  if (length(between_y_idx) > 0L) {
    rm_idx <- c(within_x_idx, between_x_idx, between_y_idx) # between AND x
    beta_z <- beta_b_tilde[, between_y_idx, drop = FALSE]
    beta_b <- beta_b_tilde[, -rm_idx, drop = FALSE]
    beta_w <- beta_w_tilde[, -rm_idx, drop = FALSE]
    sigma_zz <- sigma_b_tilde[between_y_idx, between_y_idx, drop = FALSE]
    sigma_yz <- sigma_b_tilde[-rm_idx, between_y_idx, drop = FALSE]
    sigma_b <- sigma_b_tilde[-rm_idx, -rm_idx, drop = FALSE]
    sigma_w <- sigma_w_tilde[-rm_idx, -rm_idx, drop = FALSE]
  } else {
    rm_idx <- c(within_x_idx, between_x_idx) # all 'x'
    beta_z <- matrix(0, 0L, 0L)
    sigma_zz <- matrix(0, 0L, 0L)
    beta_b <- beta_b_tilde[, -rm_idx, drop = FALSE]
    beta_w <- beta_w_tilde[, -rm_idx, drop = FALSE]
    sigma_b <- sigma_b_tilde[-rm_idx, -rm_idx, drop = FALSE]
    sigma_w <- sigma_w_tilde[-rm_idx, -rm_idx, drop = FALSE]
    sigma_yz <- matrix(0, nrow(sigma_w), 0L)
  }


  # beta.wb # FIXME: not correct if some 'x' are split (overlap)
  # but because we ALWAYS treat split-x as 'y', this is not a problem
  beta_wb <- rbind(beta_w, beta_b[-1, , drop = FALSE])
  beta_wb[1, ] <- beta_wb[1, , drop = FALSE] + beta_b[1, , drop = FALSE]

  list(
    sigma.w = sigma_w, sigma.b = sigma_b, sigma.zz = sigma_zz,
    sigma.yz = sigma_yz, beta.w = beta_w, beta.b = beta_b, beta.z = beta_z,
    beta.wb = beta_wb
  )
}


# recreate implied matrices from 2L matrices
lav_mvreg_cl_2l2implied <- function(lp,
                                         sigma_w = NULL,
                                         sigma_b = NULL,
                                         sigma_zz = NULL,
                                         sigma_yz = NULL,
                                         beta_w = NULL,
                                         beta_b = NULL,
                                         beta_z = NULL) {
  # within/between idx
  # within_x_idx <- lp$within.x.idx[[1]]
  between_y_idx <- lp$between.y.idx[[2]]
  # between_x_idx <- lp$between.x.idx[[2]]

  # ov.idx per level
  ov_idx <- lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))

  # only 'y'
  ov_y_idx <- lp$ov.y.idx

  # two levels only (for now)
  ov_y_idx1 <- ov_y_idx[[1]]
  ov_y_idx2 <- ov_y_idx[[2]]

  # Sigma.W.tilde
  sigma_w_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_w_tilde[ov_y_idx1, ov_y_idx1] <- sigma_w

  # INT.W.tilde
  int_w_tilde <- matrix(0, p_tilde, 1L)
  int_w_tilde[ov_y_idx1, 1L] <- beta_w[1L, ]

  # PI.W.tilde
  pi_w_tilde <- matrix(0, p_tilde, nrow(beta_w) - 1L)
  pi_w_tilde[ov_y_idx1, ] <- t(beta_w[-1L, ])

  # Sigma.B.tilde
  sigma_b_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_b_tilde[ov_y_idx1, ov_y_idx1] <- sigma_b

  # INT.B.tilde
  int_b_tilde <- matrix(0, p_tilde, 1L)
  int_b_tilde[ov_y_idx1, 1L] <- beta_b[1L, ]

  # PI.B.tilde
  pi_b_tilde <- matrix(0, p_tilde, nrow(beta_b) - 1L)
  pi_b_tilde[ov_y_idx1, ] <- t(beta_b[-1L, ])

  if (length(between_y_idx) > 0L) {
    int_b_tilde[between_y_idx, 1L] <- beta_z[1L, ]
    pi_b_tilde[between_y_idx, ] <- t(beta_z[-1L, ])
    sigma_b_tilde[between_y_idx, between_y_idx] <- sigma_zz
    sigma_b_tilde[ov_y_idx1, between_y_idx] <- sigma_yz
    sigma_b_tilde[between_y_idx, ov_y_idx1] <- t(sigma_yz)
  }

  res_sigma_w <- sigma_w_tilde[ov_y_idx1, ov_y_idx1, drop = FALSE]
  res_int_w <- int_w_tilde[ov_y_idx1, , drop = FALSE]
  res_pi_w <- pi_w_tilde[ov_y_idx1, , drop = FALSE]

  res_sigma_b <- sigma_b_tilde[ov_y_idx2, ov_y_idx2, drop = FALSE]
  res_int_b <- int_b_tilde[ov_y_idx2, , drop = FALSE]
  res_pi_b <- pi_b_tilde[ov_y_idx2, , drop = FALSE]

  implied <- list(
    res.cov = list(res_sigma_w, res_sigma_b),
    res.int = list(res_int_w, res_int_b),
    res.slopes = list(res_pi_w, res_pi_b)
  )

  # Note: cov.x and mean.x must be added by the caller
  implied
}

lav_mvreg_cl_loglik_samp_2l <- function(ylp = NULL,
                                                    lp = NULL,
                                                    res_sigma_w = NULL,
                                                    res_int_w = NULL,
                                                    res_pi_w = NULL,
                                                    res_sigma_b = NULL,
                                                    res_int_b = NULL,
                                                    res_pi_b = NULL,
                                                    out = NULL, # 2l
                                                    sinv_method = "eigen",
                                                    log2pi = FALSE,
                                                    minus_two = TRUE) {
  # map implied to 2l matrices
  if (is.null(out)) {
    out <- lav_mvreg_cl_implied22l(
      lp = lp, implied = NULL,
      res_sigma_w = res_sigma_w,
      res_int_w = res_int_w, res_pi_w = res_pi_w,
      res_sigma_b = res_sigma_b,
      res_int_b = res_int_b, res_pi_b = res_pi_b
    )
  }
  sigma_w <- out$sigma.w
  sigma_b <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz
  beta_w <- out$beta.w
  beta_b <- out$beta.b
  beta_z <- out$beta.z
  beta_wb <- out$beta.wb

  # check for beta.wb
  if (is.null(out$beta.wb)) {
    beta_wb <- rbind(beta_w, beta_b[-1, , drop = FALSE])
    beta_wb[1, ] <- beta_wb[1, , drop = FALSE] + beta_b[1, , drop = FALSE]
  }

  # log 2*pi
  log_2pi <- log(2 * pi)

  # Lp
  # nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  n_s <- lp$cluster.size.ns[[2]]

  # dependent 'y' level-2 ('Z') only variables?
  between_y_idx <- lp$between.y.idx[[2]]

  # extract (the many) sample statistics from YLp
  sample_wb <- ylp[[2]]$sample.wb
  sample_yyres_wb1 <- ylp[[2]]$sample.YYres.wb1
  sample_xx_wb1 <- ylp[[2]]$sample.XX.wb1
  sample_wb2 <- ylp[[2]]$sample.wb2
  sample_yyres_wb2 <- ylp[[2]]$sample.YYres.wb2
  sample_yres_x_wb2 <- ylp[[2]]$sample.YresX.wb2
  sample_xx_wb2 <- ylp[[2]]$sample.XX.wb2
  sample_clz_y2_res <- ylp[[2]]$sample.clz.Y2.res
  sample_clz_y2_xx <- ylp[[2]]$sample.clz.Y2.XX
  sample_clz_y2_b <- ylp[[2]]$sample.clz.Y2.B
  if (length(between_y_idx) > 0L) {
    sample_clz_zz_res <- ylp[[2]]$sample.clz.ZZ.res
    sample_clz_zz_xx <- ylp[[2]]$sample.clz.ZZ.XX
    sample_clz_zz_b <- ylp[[2]]$sample.clz.ZZ.B
    sample_clz_yz_res <- ylp[[2]]$sample.clz.YZ.res
    sample_clz_yz_xx <- ylp[[2]]$sample.clz.YZ.XX
    sample_clz_yres_xz <- ylp[[2]]$sample.clz.YresXZ # zero?
    sample_clz_xwzres <- ylp[[2]]$sample.clz.XWZres
  }

  # reconstruct S.PW
  wb1_diff <- sample_wb - beta_wb
  y1y1_wb_res <- (sample_yyres_wb1 +
    t(wb1_diff) %*% sample_xx_wb1 %*% (wb1_diff))

  # this one is weighted -- not the same as crossprod(Y2w.res)
  wb2_diff <- sample_wb2 - beta_wb
  y2y2w_res <- (sample_yyres_wb2 +
    sample_yres_x_wb2 %*% (wb2_diff) +
    t(wb2_diff) %*% t(sample_yres_x_wb2) +
    t(wb2_diff) %*% sample_xx_wb2 %*% (wb2_diff))
  s_pw <- (y1y1_wb_res - y2y2w_res) / sum(cluster_size - 1)

  # common parts:
  sigma_w_inv <- lav_mat_sym_inverse(
    s = sigma_w,
    logdet = TRUE
  )
  sigma_w_logdet <- attr(sigma_w_inv, "logdet")
  if (length(between_y_idx) > 0L) {
    sigma_zz_inv <- lav_mat_sym_inverse(
      s = sigma_zz,
      logdet = TRUE
    )
    sigma_zz_logdet <- attr(sigma_zz_inv, "logdet")
    sigma_yz_zi <- sigma_yz %*% sigma_zz_inv
    sigma_zi_zy <- t(sigma_yz_zi)
    sigma_b_z <- sigma_b - sigma_yz %*% sigma_zi_zy
  } else {
    sigma_b_z <- sigma_b
  }

  # min 2* logliklihood
  dist_1 <- numeric(ncluster_sizes)
  logdet <- numeric(ncluster_sizes)
  const <- numeric(ncluster_sizes)
  for (clz in seq_len(ncluster_sizes)) {
    # cluster size
    nj <- cluster_sizes[clz]

    # data between
    # nj_idx <- which(cluster_size == nj)
    y2_diff <- sample_clz_y2_b[[clz]] - beta_wb
    y2yc_yy <- (sample_clz_y2_res[[clz]] +
      t(y2_diff) %*% sample_clz_y2_xx[[clz]] %*% (y2_diff))
    if (length(between_y_idx) > 0L) {
      zz_diff <- sample_clz_zz_b[[clz]] - beta_z
      y2yc_zz <- (sample_clz_zz_res[[clz]] +
        t(zz_diff) %*% sample_clz_zz_xx[[clz]] %*% (zz_diff))
      y2yc_yz <- (sample_clz_yz_res[[clz]] +
        sample_clz_yres_xz[[clz]] %*% zz_diff + # zero?
        t(y2_diff) %*% sample_clz_xwzres[[clz]] +
        t(y2_diff) %*% sample_clz_yz_xx[[clz]] %*% zz_diff)
    }

    # construct sigma.j
    sigma_j <- (nj * sigma_b_z) + sigma_w
    sigma_j_inv <- lav_mat_sym_inverse(
      s = sigma_j,
      logdet = TRUE
    )
    sigma_j_logdet <- attr(sigma_j_inv, "logdet")

    if (length(between_y_idx) > 0L) {
      sigma_ji_yz_zi <- sigma_j_inv %*% sigma_yz_zi

      # part 1 -- zz
      vinv_11 <- sigma_zz_inv + nj * (sigma_zi_zy %*% sigma_ji_yz_zi)
      q_zz <- sum(vinv_11 * y2yc_zz)

      # part 2 -- yz
      q_yz <- -nj * sum(sigma_ji_yz_zi * y2yc_yz)
    } else {
      q_zz <- q_yz <- sigma_zz_logdet <- 0
    }

    # part 5 -- yyc
    q_yyc <- -nj * sum(sigma_j_inv * y2yc_yy)

    if (log2pi) {
      p <- nj * nrow(sigma_w) + nrow(sigma_zz)
      const[clz] <- p * log_2pi
    }
    logdet[clz] <- sigma_zz_logdet + sigma_j_logdet
    dist_1[clz] <- q_zz + 2 * q_yz - q_yyc
  }
  # q.yya + q.yyb
  q_w <- sum(cluster_size - 1) * sum(sigma_w_inv * s_pw)
  # logdet within part
  l_w <- sum(cluster_size - 1) * sigma_w_logdet

  # -2*times logl (without the constant) (for optimization)
  loglik <- sum(logdet * n_s) + sum(dist_1) + q_w + l_w

  if (log2pi) {
    loglik <- loglik + sum(const * n_s)
  }

  # functions below compute -2 * logl
  if (!minus_two) {
    loglik <- loglik / (-2)
  }

  loglik
}

# first derivative -2*logl wrt Beta.W, Beta.B, Sigma.W, Sigma.B
lav_mvreg_cl_dlogl_2l_samp <- function(ylp = NULL,
                                                   lp = NULL,
                                                   res_sigma_w = NULL,
                                                   res_int_w = NULL,
                                                   res_pi_w = NULL,
                                                   res_sigma_b = NULL,
                                                   res_int_b = NULL,
                                                   res_pi_b = NULL,
                                                   out = NULL, # 2l
                                                   return_list = FALSE,
                                                   sinv_method = "eigen") {
  # map implied to 2l matrices
  if (is.null(out)) {
    out <- lav_mvreg_cl_implied22l(
      lp = lp, implied = NULL,
      res_sigma_w = res_sigma_w,
      res_int_w = res_int_w, res_pi_w = res_pi_w,
      res_sigma_b = res_sigma_b,
      res_int_b = res_int_b, res_pi_b = res_pi_b
    )
  }
  sigma_w <- out$sigma.w
  sigma_b <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz
  beta_w <- out$beta.w
  beta_b <- out$beta.b
  beta_z <- out$beta.z
  beta_wb <- out$beta.wb

  # check for beta.wb
  if (is.null(out$beta.wb)) {
    beta_wb <- rbind(beta_w, beta_b[-1, , drop = FALSE])
    beta_wb[1, ] <- beta_wb[1, , drop = FALSE] + beta_b[1, , drop = FALSE]
  }

  # Lp
  # nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  n_s <- lp$cluster.size.ns[[2]]

  within_x_idx <- lp$within.x.idx[[1]]
  between_y_idx <- lp$between.y.idx[[2]]

  w1_idx <- seq_len(length(within_x_idx) + 1L)
  b1_idx <- c(1L, seq_len(nrow(beta_wb))[-w1_idx])

  # extract (the many) sample statistics from YLp
  sample_wb <- ylp[[2]]$sample.wb
  sample_yyres_wb1 <- ylp[[2]]$sample.YYres.wb1
  sample_xx_wb1 <- ylp[[2]]$sample.XX.wb1
  sample_wb2 <- ylp[[2]]$sample.wb2
  sample_yyres_wb2 <- ylp[[2]]$sample.YYres.wb2
  sample_yres_x_wb2 <- ylp[[2]]$sample.YresX.wb2
  sample_xx_wb2 <- ylp[[2]]$sample.XX.wb2
  sample_clz_y2_res <- ylp[[2]]$sample.clz.Y2.res
  sample_clz_y2_xx <- ylp[[2]]$sample.clz.Y2.XX
  sample_clz_y2_b <- ylp[[2]]$sample.clz.Y2.B
  if (length(between_y_idx) > 0L) {
    sample_clz_zz_res <- ylp[[2]]$sample.clz.ZZ.res
    sample_clz_zz_xx <- ylp[[2]]$sample.clz.ZZ.XX
    sample_clz_zz_b <- ylp[[2]]$sample.clz.ZZ.B
    sample_clz_yz_res <- ylp[[2]]$sample.clz.YZ.res
    sample_clz_yz_xx <- ylp[[2]]$sample.clz.YZ.XX
    sample_clz_yres_xz <- ylp[[2]]$sample.clz.YresXZ # zero?
    sample_clz_xwzres <- ylp[[2]]$sample.clz.XWZres
  }

  # reconstruct S.PW
  wb1_diff <- sample_wb - beta_wb
  y1y1_wb_res <- (sample_yyres_wb1 +
    t(wb1_diff) %*% sample_xx_wb1 %*% (wb1_diff))

  # this one is weighted -- not the same as crossprod(Y2w.res)
  wb2_diff <- sample_wb2 - beta_wb
  y2y2w_res <- (sample_yyres_wb2 +
    sample_yres_x_wb2 %*% (wb2_diff) +
    t(wb2_diff) %*% t(sample_yres_x_wb2) +
    t(wb2_diff) %*% sample_xx_wb2 %*% (wb2_diff))
  s_pw <- (y1y1_wb_res - y2y2w_res) / sum(cluster_size - 1)

  # common parts:
  sigma_w_inv <- lav_mat_sym_inverse(s = sigma_w)

  g_beta_w <- matrix(0, ncluster_sizes, length(beta_w))
  g_beta_b <- matrix(0, ncluster_sizes, length(beta_b))
  # g_beta_wb <- matrix(0, ncluster_sizes, length(beta_wb))
  g_sigma_w1_1 <- matrix(0, ncluster_sizes, length(lav_mat_vech(sigma_w)))
  g_sigma_b_1 <- matrix(0, ncluster_sizes, length(lav_mat_vech(sigma_b)))

  if (length(between_y_idx) > 0L) {
    g_beta_z <- matrix(0, ncluster_sizes, length(beta_z))
    g_sigma_zz_1 <- matrix(0, ncluster_sizes, length(lav_mat_vech(sigma_zz)))
    g_sigma_yz_1 <- matrix(0, ncluster_sizes, length(sigma_yz))

    sigma_zz_inv <- lav_mat_sym_inverse(s = sigma_zz)
    sigma_yz_zi <- sigma_yz %*% sigma_zz_inv
    sigma_zi_zy <- t(sigma_yz_zi)
    sigma_b_z <- sigma_b - sigma_yz %*% sigma_zi_zy

    for (clz in seq_len(ncluster_sizes)) {
      # cluster size
      nj <- cluster_sizes[clz]

      y2_diff <- sample_clz_y2_b[[clz]] - beta_wb
      xx_y2_diff <- sample_clz_y2_xx[[clz]] %*% y2_diff
      y2yc_yy <- sample_clz_y2_res[[clz]] + crossprod(y2_diff, xx_y2_diff)

      zz_diff <- sample_clz_zz_b[[clz]] - beta_z
      y2yc_zz <- (sample_clz_zz_res[[clz]] +
        t(zz_diff) %*% sample_clz_zz_xx[[clz]] %*% (zz_diff))
      y2yc_yz <- (sample_clz_yz_res[[clz]] +
        sample_clz_yres_xz[[clz]] %*% zz_diff + # zero?
        t(y2_diff) %*% sample_clz_xwzres[[clz]] +
        t(y2_diff) %*% sample_clz_yz_xx[[clz]] %*% zz_diff)

      # construct sigma.j
      sigma_j <- (nj * sigma_b_z) + sigma_w
      sigma_j_inv <- lav_mat_sym_inverse(s = sigma_j)
      sigma_ji_yz_zi <- sigma_j_inv %*% sigma_yz_zi
      sigma_zi_zy_ji <- t(sigma_ji_yz_zi)
      sigma_ji_yz <- sigma_j_inv %*% sigma_yz
      ns_sigma_j_inv <- n_s[clz] * sigma_j_inv
      ns_sigma_zz_inv <- n_s[clz] * sigma_zz_inv
      # ns_sigma_yz <- n_s[clz] * sigma_yz
      ns_sigma_ji_yz_zi <- n_s[clz] * sigma_ji_yz_zi

      # common parts
      zz_zi_yz_ji <- y2yc_zz %*% sigma_zi_zy_ji
      ji_yz_zi <- sigma_j_inv %*% y2yc_yz %*% sigma_zz_inv

      j_yzj_yy <- sigma_j_inv %*% y2yc_yy %*% sigma_j_inv
      j_yzj_yz <- tcrossprod(ji_yz_zi, sigma_ji_yz)
      j_yzj_zz <- sigma_ji_yz_zi %*% zz_zi_yz_ji

      j_yzj <- nj * (j_yzj_yy + j_yzj_zz - j_yzj_yz - t(j_yzj_yz))

      # SIGMA.W (between part)
      g_sigma_w1 <- ns_sigma_j_inv - j_yzj
      tmp <- g_sigma_w1 * 2
      diag(tmp) <- diag(g_sigma_w1)
      g_sigma_w1_1[clz, ] <- lav_mat_vech(tmp)

      # SIGMA.B
      g_sigma_b <- nj * g_sigma_w1
      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      g_sigma_b_1[clz, ] <- lav_mat_vech(tmp)

      # SIGMA.ZZ
      yz1 <- zz_zi_yz_ji %*% sigma_yz
      yz2 <- crossprod(y2yc_yz, sigma_ji_yz)
      tmp <- (t(sigma_yz) %*% g_sigma_w1 %*% sigma_yz
        - 1 / nj * y2yc_zz - t(yz1) - yz1 + t(yz2) + yz2)
      g_sigma_zz <- (ns_sigma_zz_inv +
        nj * sigma_zz_inv %*% tmp %*% sigma_zz_inv)
      tmp <- g_sigma_zz * 2
      diag(tmp) <- diag(g_sigma_zz)
      g_sigma_zz_1[clz, ] <- lav_mat_vech(tmp)

      # SIGMA.ZY
      tmp1 <- crossprod(zz_zi_yz_ji, sigma_zz_inv)
      tmp2 <- ns_sigma_ji_yz_zi
      tmp3 <- ji_yz_zi
      tmp4 <- j_yzj %*% sigma_yz_zi
      g_sigma_yz <- 2 * nj * (tmp1 - tmp2 - tmp3 + tmp4)
      g_sigma_yz_1[clz, ] <- lav_mat_vec(g_sigma_yz)

      # BETA.Z
      a <- (sigma_zz_inv + nj * (sigma_zi_zy_ji %*% sigma_yz_zi)) # symm!
      m_b <- nj * (sigma_zi_zy_ji)
      tmp_z <- (sample_clz_zz_xx[[clz]] %*% zz_diff %*% a -
        (t(sample_clz_yres_xz[[clz]]) +
          t(sample_clz_yz_xx[[clz]]) %*% y2_diff) %*% t(m_b))
      g_beta_z[clz, ] <- as.vector(-2 * tmp_z)

      # BETA.W (between part only) + BETA.B
      tmp <- (sample_clz_xwzres[[clz]] +
        sample_clz_yz_xx[[clz]] %*% zz_diff)
      out_b <- tmp %*% sigma_zi_zy_ji - xx_y2_diff %*% sigma_j_inv
      out_w <- out_b + xx_y2_diff %*% sigma_w_inv
      tmp_b <- out_b[b1_idx, , drop = FALSE]
      tmp_w <- out_w[w1_idx, , drop = FALSE]
      g_beta_b[clz, ] <- as.vector(2 * nj * tmp_b)
      g_beta_w[clz, ] <- as.vector(2 * nj * tmp_w)
    } # clz

    d_beta_w1 <- matrix(colSums(g_beta_w), nrow(beta_w), ncol(beta_w))
    d_beta_b <- matrix(colSums(g_beta_b), nrow(beta_b), ncol(beta_b))
    d_sigma_w1 <- lav_mat_vech_rev(colSums(g_sigma_w1_1))
    d_sigma_b <- lav_mat_vech_rev(colSums(g_sigma_b_1))

    # z
    d_beta_z <- matrix(colSums(g_beta_z), nrow(beta_z), ncol(beta_z))
    d_sigma_zz <- lav_mat_vech_rev(colSums(g_sigma_zz_1))
    d_sigma_yz <- matrix(colSums(g_sigma_yz_1), nrow(sigma_yz), ncol(sigma_yz))
  # between.y.idx

  } else { # no between.y.idx

    for (clz in seq_len(ncluster_sizes)) {
      # cluster size
      nj <- cluster_sizes[clz]

      y2_diff <- sample_clz_y2_b[[clz]] - beta_wb
      xx_y2_diff <- sample_clz_y2_xx[[clz]] %*% y2_diff
      y2yc_yy <- sample_clz_y2_res[[clz]] + crossprod(y2_diff, xx_y2_diff)

      # construct sigma.j
      sigma_j <- (nj * sigma_b) + sigma_w
      sigma_j_inv <- lav_mat_sym_inverse(s = sigma_j)

      # common part
      j_yyj <- nj * sigma_j_inv %*% y2yc_yy %*% sigma_j_inv

      # SIGMA.W (between part)
      g_sigma_w1 <- (n_s[clz] * sigma_j_inv) - j_yyj
      tmp <- g_sigma_w1 * 2
      diag(tmp) <- diag(g_sigma_w1)
      g_sigma_w1_1[clz, ] <- lav_mat_vech(tmp)

      # SIGMA.B
      g_sigma_b <- nj * g_sigma_w1
      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      g_sigma_b_1[clz, ] <- lav_mat_vech(tmp)

      # BETA.W (between part only) + BETA.B
      out_b <- -1 * xx_y2_diff %*% sigma_j_inv
      out_w <- out_b + xx_y2_diff %*% sigma_w_inv
      tmp_b <- out_b[b1_idx, , drop = FALSE]
      tmp_w <- out_w[w1_idx, , drop = FALSE]
      g_beta_b[clz, ] <- as.vector(2 * nj * tmp_b)
      g_beta_w[clz, ] <- as.vector(2 * nj * tmp_w)
    } # cl

    d_beta_w1 <- matrix(colSums(g_beta_w), nrow(beta_w), ncol(beta_w))
    d_beta_b <- matrix(colSums(g_beta_b), nrow(beta_b), ncol(beta_b))
    d_sigma_w1 <- lav_mat_vech_rev(colSums(g_sigma_w1_1))
    d_sigma_b <- lav_mat_vech_rev(colSums(g_sigma_b_1))

    # z
    d_beta_z <- matrix(0, 0L, 0L)
    d_sigma_zz <- matrix(0, 0L, 0L)
    d_sigma_yz <- matrix(0, 0L, 0L)
  } # no-between-y

  # Sigma.W (bis)
  d_sigma_w2 <- sum(cluster_size - 1) * (sigma_w_inv
  - sigma_w_inv %*% s_pw %*% sigma_w_inv)
  tmp <- d_sigma_w2 * 2
  diag(tmp) <- diag(d_sigma_w2)
  d_sigma_w2 <- tmp

  d_sigma_w <- d_sigma_w1 + d_sigma_w2

  # beta.w (bis)
  d_beta_w2 <- -2 * (sample_xx_wb1 %*%
      (sample_wb - beta_wb))[w1_idx, , drop = FALSE] %*% sigma_w_inv

  d_beta_w <- d_beta_w1 + d_beta_w2

  # rearrange
  dimplied <- lav_mvreg_cl_2l2implied(lp,
    sigma_w = d_sigma_w, sigma_b = d_sigma_b,
    sigma_zz = d_sigma_zz, sigma_yz = d_sigma_yz,
    beta_w = d_beta_w, beta_b = d_beta_b, beta_z = d_beta_z
  )

  if (return_list) {
    return(dimplied)
  }

  # as a single vector
  out <- c(
    drop(dimplied$res.int[[1]]),
    lav_mat_vec(dimplied$res.slopes[[1]]),
    lav_mat_vech(dimplied$res.cov[[1]]),
    drop(dimplied$res.int[[2]]),
    lav_mat_vec(dimplied$res.slopes[[2]]),
    lav_mat_vech(dimplied$res.cov[[2]])
  )
  out
}

# cluster-wise scores -2*logl wrt Beta.W, Beta.B, Sigma.W, Sigma.B
lav_mvreg_cl_sc_2l <- function(y1 = NULL,
                                        ylp = NULL,
                                        lp = NULL,
                                        res_sigma_w = NULL,
                                        res_int_w = NULL,
                                        res_pi_w = NULL,
                                        res_sigma_b = NULL,
                                        res_int_b = NULL,
                                        res_pi_b = NULL,
                                        out = NULL, # 2l
                                        sinv_method = "eigen") {
  # map implied to 2l matrices
  if (is.null(out)) {
    out <- lav_mvreg_cl_implied22l(
      lp = lp, implied = NULL,
      res_sigma_w = res_sigma_w,
      res_int_w = res_int_w, res_pi_w = res_pi_w,
      res_sigma_b = res_sigma_b,
      res_int_b = res_int_b, res_pi_b = res_pi_b
    )
  }
  sigma_w <- out$sigma.w
  sigma_b <- out$sigma.b
  sigma_zz <- out$sigma.zz
  sigma_yz <- out$sigma.yz
  beta_w <- out$beta.w
  beta_b <- out$beta.b
  beta_z <- out$beta.z
  beta_wb <- out$beta.wb

  # check for beta.wb
  if (is.null(out$beta.wb)) {
    beta_wb <- rbind(beta_w, beta_b[-1, , drop = FALSE])
    beta_wb[1, ] <- beta_wb[1, , drop = FALSE] + beta_b[1, , drop = FALSE]
  }

  # Lp
  nclusters <- lp$nclusters[[2]]
  cluster_size <- lp$cluster.size[[2]]
  cluster_idx <- lp$cluster.idx[[2]]

  within_x_idx <- lp$within.x.idx[[1]]
  # between_idx <- lp$between.idx[[2]]
  between_y_idx <- lp$between.y.idx[[2]]
  between_x_idx <- lp$between.x.idx[[2]]

  y1_idx <- lp$ov.y.idx[[1]]
  x1_idx <- c(within_x_idx, between_x_idx) # in that order

  # residuals for 'Y'
  y1_wb <- y1[, y1_idx, drop = FALSE]

  if (length(x1_idx) > 0L) {
    exo_wb <- cbind(1, y1[, x1_idx, drop = FALSE])
    y1_wb_hat <- exo_wb %*% beta_wb
    y1_wb_res <- y1_wb - y1_wb_hat
  } else {
    y1_wb_res <- y1_wb
  }

  # residuals 'Y' (level 2)
  y2 <- ylp[[2]]$Y2
  if (length(x1_idx) > 0L) {
    exo_wb2 <- cbind(1, y2[, x1_idx, drop = FALSE])
    y2w_res <- y2[, y1_idx, drop = FALSE] - exo_wb2 %*% beta_wb
  } else {
    exo_wb2 <- matrix(1, nrow(y2), 1L)
    y2w_res <- y2[, y1_idx, drop = FALSE]
  }

  # residual 'Z' (level 2)
  if (length(between_y_idx) > 0L) {
    if (length(between_x_idx) > 0L) {
      exo_z_1 <- cbind(1, y2[, between_x_idx, drop = FALSE])
      y2_z <- y2[, between_y_idx, drop = FALSE]
      y2z_res <- y2_z - exo_z_1 %*% beta_z
      # sample.z
      # XX.z <- crossprod(EXO.z)
      # sample.z <- try(solve(XX.z, crossprod(EXO.z, Y2.z)))
      # if(inherits(sample.z, "try-error")) {
      #    sample.z <- MASS::ginv(XX.z) %*% crossprod(EXO.z, Y2.z)
      # }

      # sample.wb2
      # sample.wb2 <- YLp[[2]]$sample.wb2
    } else {
      y2z_res <- y2[, between_y_idx, drop = FALSE]
    }
  }

  # common parts:
  sigma_w_inv <- lav_mat_sym_inverse(s = sigma_w)

  g_beta_w1 <- matrix(0, nclusters, length(beta_w))
  g_beta_b <- matrix(0, nclusters, length(beta_b))
  # g_beta_wb <- matrix(0, nclusters, length(beta_wb))
  g_sigma_w1_1 <- matrix(0, nclusters, length(lav_mat_vech(sigma_w)))
  g_sigma_b_1 <- matrix(0, nclusters, length(lav_mat_vech(sigma_b)))


  if (length(between_y_idx) > 0L) {
    g_beta_z <- matrix(0, nclusters, length(beta_z))
    g_sigma_zz_1 <- matrix(0, nclusters, length(lav_mat_vech(sigma_zz)))
    g_sigma_yz_1 <- matrix(0, nclusters, length(sigma_yz))

    sigma_zz_inv <- lav_mat_sym_inverse(s = sigma_zz)
    sigma_yz_zi <- sigma_yz %*% sigma_zz_inv
    sigma_zi_zy <- t(sigma_yz_zi)
    sigma_b_z <- sigma_b - sigma_yz %*% sigma_zi_zy

    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster_size[cl]

      # data within for the cluster (centered)
      # y1m <- y1_wb_res[cluster_idx == cl, , drop = FALSE]
      yc <- y2w_res[cl, ]

      # data between
      zc <- y2z_res[cl, ]
      y2yc_yy <- tcrossprod(y2w_res[cl, ])
      y2yc_zz <- tcrossprod(y2z_res[cl, ])
      y2yc_yz <- tcrossprod(y2w_res[cl, ], y2z_res[cl, ])

      # construct sigma.j
      sigma_j <- (nj * sigma_b_z) + sigma_w
      sigma_j_inv <- lav_mat_sym_inverse(s = sigma_j)
      sigma_ji_yz_zi <- sigma_j_inv %*% sigma_yz_zi
      sigma_zi_zy_ji <- t(sigma_ji_yz_zi)
      sigma_ji_yz <- sigma_j_inv %*% sigma_yz

      # common parts
      zz_zi_yz_ji <- y2yc_zz %*% sigma_zi_zy_ji
      ji_yz_zi <- sigma_j_inv %*% y2yc_yz %*% sigma_zz_inv

      j_yzj_yy <- sigma_j_inv %*% y2yc_yy %*% sigma_j_inv
      j_yzj_yz <- tcrossprod(ji_yz_zi, sigma_ji_yz)
      j_yzj_zz <- sigma_ji_yz_zi %*% zz_zi_yz_ji

      j_yzj <- nj * (j_yzj_yy + j_yzj_zz - j_yzj_yz - t(j_yzj_yz))

      # SIGMA.W (between part)
      g_sigma_w1 <- sigma_j_inv - j_yzj
      tmp <- g_sigma_w1 * 2
      diag(tmp) <- diag(g_sigma_w1)
      g_sigma_w1_1[cl, ] <- lav_mat_vech(tmp)

      # SIGMA.W (within part)
      # g.sigma.w2 <- ( (nj-1) * sigma.w.inv
      #    - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv )
      # tmp <- g.sigma.w2*2; diag(tmp) <- diag(g.sigma.w2)
      # G.sigma.w2[cl,] <- lav_mat_vech(tmp)
      # G.sigma.w[cl,] <- G.sigma.w1[cl,] + G.sigma.w2[cl,]

      # SIGMA.B
      g_sigma_b <- nj * g_sigma_w1
      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      g_sigma_b_1[cl, ] <- lav_mat_vech(tmp)

      # SIGMA.ZZ
      yz1 <- zz_zi_yz_ji %*% sigma_yz
      yz2 <- crossprod(y2yc_yz, sigma_ji_yz)
      tmp <- (t(sigma_yz) %*% g_sigma_w1 %*% sigma_yz
        - (1 / nj * y2yc_zz + t(yz1) + yz1 - t(yz2) - yz2))
      g_sigma_zz <- (sigma_zz_inv +
        nj * sigma_zz_inv %*% tmp %*% sigma_zz_inv)
      tmp <- g_sigma_zz * 2
      diag(tmp) <- diag(g_sigma_zz)
      g_sigma_zz_1[cl, ] <- lav_mat_vech(tmp)

      # SIGMA.ZY
      # g.sigma.yz <- 2 * nj * (
      #              (sigma.j.inv %*%
      #                  (sigma.yz.zi %*% Y2Yc.zz - sigma.yz - Y2Yc.yz)
      #                   + jYZj %*% sigma.yz) %*% sigma.zz.inv )
      tmp1 <- crossprod(zz_zi_yz_ji, sigma_zz_inv)
      tmp2 <- sigma_ji_yz_zi
      tmp3 <- ji_yz_zi
      tmp4 <- j_yzj %*% sigma_yz_zi
      g_sigma_yz <- 2 * nj * (tmp1 - tmp2 - tmp3 + tmp4)
      g_sigma_yz_1[cl, ] <- lav_mat_vec(g_sigma_yz)

      # BETA.Z
      # here, we avoid the (sample.z - beta.z) approach
      exo_z <- cbind(1, y2[cl, between_x_idx, drop = FALSE])
      tmp1 <- (sigma_zz_inv + nj * (sigma_zi_zy_ji %*% sigma_yz_zi)) %*% zc
      tmp2 <- nj * (sigma_zi_zy_ji) %*% yc
      tmp_z <- crossprod(exo_z, drop(tmp1 - tmp2))
      g_beta_z[cl, ] <- as.vector(-2 * tmp_z)

      # BETA.W
      #       exo.w <- cbind(1,
      #                Y1[cluster.idx == cl, within.x.idx, drop = FALSE])
      #       G.beta.w[cl,] <- as.vector( 2 * t(exo.w) %*% (
      #                        matrix(1, nj, 1) %x% (zc %*% sigma.zi.zy.ji -
      #                                              yc %*% sigma.j.inv +
      #                                              yc %*% sigma.w.inv) -
      #                            Y1m %*% sigma.w.inv) )

      # BETA.W (between part only)
      exo2_w <- cbind(1, y2[cl, within_x_idx, drop = FALSE])
      tmp2 <- (zc %*% sigma_zi_zy_ji -
        yc %*% sigma_j_inv +
        yc %*% sigma_w_inv)
      g_beta_w1[cl, ] <- as.vector(2 * nj * crossprod(exo2_w, tmp2))

      # BETA.W (within part only)
      # exo.w <- cbind(1,
      #               Y1[cluster.idx == cl, within.x.idx, drop = FALSE])
      # tmp1 <- - Y1m %*% sigma.w.inv
      # G.beta.ww <- as.vector( 2 * crossprod(exo.w, tmp1) )
      # G.beta.w[cl,] <- G.beta.w1 + G.beta.ww
      # G.beta.w2[cl,] <- G.beta.ww

      # BETA.B
      exo2_b <- cbind(1, y2[cl, between_x_idx, drop = FALSE])
      tmp <- (zc %*% sigma_zi_zy_ji - yc %*% sigma_j_inv)
      g_beta_b[cl, ] <- as.vector(2 * nj * crossprod(exo2_b, tmp))
    } # cl
  # between.y.idx

  } else { # no between.y.idx

    for (cl in seq_len(nclusters)) {
      # cluster size
      nj <- cluster_size[cl]

      # data within for the cluster (centered)
      # y1m <- y1_wb_res[cluster_idx == cl, , drop = FALSE]
      yc <- y2w_res[cl, ]

      # data between
      y2yc_yy <- tcrossprod(y2w_res[cl, ])

      # construct sigma.j
      sigma_j <- (nj * sigma_b) + sigma_w
      sigma_j_inv <- lav_mat_sym_inverse(s = sigma_j)

      # common part
      j_yyj <- nj * sigma_j_inv %*% y2yc_yy %*% sigma_j_inv

      # SIGMA.W
      # g.sigma.w <- ( (nj-1) * sigma.w.inv
      #    - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv
      #    + sigma.j.inv - jYYj )
      # tmp <- g.sigma.w*2; diag(tmp) <- diag(g.sigma.w)
      # G.sigma.w[cl,] <- lav_mat_vech(tmp)

      # SIGMA.W (between part)
      g_sigma_w1 <- sigma_j_inv - j_yyj
      tmp <- g_sigma_w1 * 2
      diag(tmp) <- diag(g_sigma_w1)
      g_sigma_w1_1[cl, ] <- lav_mat_vech(tmp)

      # SIGMA.B
      g_sigma_b <- nj * (sigma_j_inv - j_yyj)
      tmp <- g_sigma_b * 2
      diag(tmp) <- diag(g_sigma_b)
      g_sigma_b_1[cl, ] <- lav_mat_vech(tmp)

      # BETA.W (between part only)
      exo2_w <- cbind(1, y2[cl, within_x_idx, drop = FALSE])
      tmp2 <- (-yc %*% sigma_j_inv + yc %*% sigma_w_inv)
      g_beta_w1[cl, ] <- as.vector(2 * nj * crossprod(exo2_w, tmp2))

      # BETA.B
      exo2_b <- cbind(1, y2[cl, between_x_idx, drop = FALSE])
      tmp <- -yc %*% sigma_j_inv
      g_beta_b[cl, ] <- as.vector(2 * nj * crossprod(exo2_b, tmp))
    } # cl
  } # no-between-y

  # beta.w (bis)

  #    d.beta.w2 <- -2 * t(EXO.wb[,1:(length(within.x.idx) + 1L), drop = FALSE])
  #                                                %*% Y1.wb.res %*% sigma.w.inv

  y1_wb_res_i <- y1_wb_res %*% sigma_w_inv
  w1_idx <- seq_len(length(within_x_idx) + 1L)
  a1_idx <- rep(w1_idx, times = ncol(y1_wb_res_i))
  b1_idx <- rep(seq_len(ncol(y1_wb_res_i)), each = length(w1_idx))
  tmp_1 <- exo_wb[, a1_idx, drop = FALSE] * y1_wb_res_i[, b1_idx, drop = FALSE]
  g_beta_w2 <- -2 * rowsum.default(tmp_1, cluster_idx,
    reorder = FALSE,
    na.rm = TRUE
  )
  g_beta_w <- g_beta_w1 + g_beta_w2

  # Sigma.W (bis)
  # d.sigma.w2 <- sum(cluster.size - 1) * ( sigma.w.inv
  #               - sigma.w.inv %*% S.PW %*% sigma.w.inv )
  # tmp <- d.sigma.w2*2; diag(tmp) <- diag(d.sigma.w2)
  # d.sigma.w2 <- tmp

  # g.sigma.w2 <- ( (nj-1) * sigma.w.inv
  #        - sigma.w.inv %*% (crossprod(Y1m) - nj*Y2Yc.yy) %*% sigma.w.inv )

  y1a_res <- y1_wb_res - y2w_res[cluster_idx, , drop = FALSE]
  y1a_res_i <- y1a_res %*% sigma_w_inv
  idx1 <- lav_mat_vech_col_idx(nrow(sigma_w))
  idx2 <- lav_mat_vech_row_idx(nrow(sigma_w))
  sw2 <- matrix(lav_mat_vech(sigma_w_inv),
    nrow = nclusters,
    length(lav_mat_vech(sigma_w_inv)), byrow = TRUE
  )
  sw2 <- sw2 * (cluster_size - 1)
  tmp_1 <- y1a_res_i[, idx1, drop = FALSE] * y1a_res_i[, idx2, drop = FALSE]
  tmp2_1 <- rowsum.default(tmp_1, cluster_idx, reorder = FALSE, na.rm = TRUE)
  g_sigma_w2 <- 2 * (sw2 - tmp2_1)
  diagh_idx <- lav_mat_diagh_idx(nrow(sigma_w))
  g_sigma_w2[, diagh_idx] <- g_sigma_w2[, diagh_idx, drop = FALSE] / 2
  g_sigma_w <- g_sigma_w1_1 + g_sigma_w2



  # rearrange columns to Res.Int.W, Res.Pi.W, Res.Sigma.W,
  #                      Res.Int.B, Res.Pi.B, Res.Sigma.B

  # ov.idx per level
  ov_idx <- lp$ov.idx

  # 'tilde' matrices: ALL variables within and between
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))
  p_tilde_star <- p_tilde * (p_tilde + 1) / 2
  b_tilde <- lav_mat_vech_rev(seq_len(p_tilde_star))

  # only 'y'
  ov_y_idx <- lp$ov.y.idx

  # two levels only (for now)
  ov_y_idx1 <- ov_y_idx[[1]]
  ov_y_idx2 <- ov_y_idx[[2]]

  # WITHIN (is easy)
  beta_w_idx <- matrix(seq_along(beta_w), nrow(beta_w), ncol(beta_w))
  beta_b_idx <- matrix(seq_along(beta_b), nrow(beta_b), ncol(beta_b))

  res_int_w <- g_beta_w[, beta_w_idx[1L, ], drop = FALSE]
  res_pi_w <- g_beta_w[, lav_mat_vecr(beta_w_idx[-1L, ]), drop = FALSE]
  res_sigma_w <- g_sigma_w

  # Sigma.B
  sigma_b_tilde <- matrix(0, nclusters, p_tilde_star)
  col_idx <- lav_mat_vech(b_tilde[ov_y_idx1, ov_y_idx1, drop = FALSE])
  sigma_b_tilde[, col_idx] <- g_sigma_b_1

  # Int.B
  beta_b_tilde <- matrix(seq_len(nrow(beta_b) * p_tilde), nrow(beta_b), p_tilde)
  int_b <- matrix(0, nclusters, p_tilde)
  int_b[, ov_y_idx1] <- g_beta_b[, beta_b_idx[1L, ]]

  # Pi.B
  pi_b <- matrix(0, nclusters, p_tilde * (nrow(beta_b) - 1L))
  col_idx <- lav_mat_vecr(beta_b_tilde[-1L, ov_y_idx1, drop = FALSE])
  pi_b[, col_idx] <- g_beta_b[, lav_mat_vecr(beta_b_idx[-1L, ]),
                                                                drop = FALSE]

  if (length(between_y_idx) > 0L) {
    # Sigma.B: add yz/zz parts
    col_idx <- lav_mat_vec(b_tilde[ov_y_idx1, between_y_idx, drop = FALSE])
    sigma_b_tilde[, col_idx] <- g_sigma_yz_1
    col_idx <- lav_mat_vech(b_tilde[between_y_idx, between_y_idx,
      drop = FALSE
    ])
    sigma_b_tilde[, col_idx] <- g_sigma_zz_1

    # Int.B: add z-part
    beta_z_idx <- matrix(seq_along(beta_z), nrow(beta_z), ncol(beta_z))
    int_b[, between_y_idx] <- g_beta_z[, beta_z_idx[1L, ], drop = FALSE]

    # Pi.B: add beta.z
    col_idx <- lav_mat_vecr(beta_b_tilde[-1L, between_y_idx, drop = FALSE])
    pi_b[, col_idx] <-
      g_beta_z[, lav_mat_vecr(beta_z_idx[-1L, ]), drop = FALSE]
  }

  # only extract ov.y.idx2 for BETWEEN
  col_idx <- lav_mat_vech(b_tilde[ov_y_idx2, ov_y_idx2, drop = FALSE])
  res_sigma_b <- sigma_b_tilde[, col_idx, drop = FALSE]

  res_int_b <- int_b[, ov_y_idx2, drop = FALSE]

  col_idx <- lav_mat_vecr(beta_b_tilde[-1, ov_y_idx2])
  res_pi_b <- pi_b[, col_idx, drop = FALSE]

  scores <- cbind(
    res_int_w, res_pi_w, res_sigma_w,
    res_int_b, res_pi_b, res_sigma_b
  )

  scores
}

# first-order information: outer crossprod of scores per cluster
lav_mvreg_cl_info_firstorder <- function(y1 = NULL,
                                                     ylp = NULL,
                                                     lp = NULL,
                                                     res_sigma_w = NULL,
                                                     res_int_w = NULL,
                                                     res_pi_w = NULL,
                                                     res_sigma_b = NULL,
                                                     res_int_b = NULL,
                                                     res_pi_b = NULL,
                                                     divide_by_two = FALSE,
                                                     sinv_method = "eigen") {
  # n <- NROW(y1)

  scores <- lav_mvreg_cl_sc_2l(
    y1 = y1,
    ylp = ylp,
    lp = lp,
    res_sigma_w = res_sigma_w,
    res_int_w = res_int_w,
    res_pi_w = res_pi_w,
    res_sigma_b = res_sigma_b,
    res_int_b = res_int_b,
    res_pi_b = res_pi_b,
    sinv_method = sinv_method
  )

  # divide by 2 (if we want scores wrt objective function)
  if (divide_by_two) {
    scores <- scores / 2
  }

  # unit information
  information <- crossprod(scores) / lp$nclusters[[2]]

  information
}
