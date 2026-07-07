# shared computational engines for the two-level (clustered) mvnorm/mvreg
# machinery
#
# the gradient of -2*logl and the cluster-wise scores share (almost) all
# of their algebra; the engines below compute both, controlled by a
# score_mode= switch:
#   - score_mode = FALSE: one 'unit' per cluster SIZE, aggregate (dlogl)
#   - score_mode = TRUE : one 'unit' per cluster, cluster-wise rows (scores)
#
# YR/CC 03 July 2026: first version (refactoring; no change in behavior)


# 0) shared scatter/gather core of the implied <-> 2l conversions
#    (see the representation invariants at the top of lav_mvnorm_cluster.R)
#
#    the sigma blocks of lav_mvn_cl_implied22l/lav_mvreg_cl_implied22l and
#    their inverses only differ in which index sets partition the tilde
#    universe:
#      mvnorm: l1_idx/l2_idx = ov.idx[[1]]/[[2]], z_idx = between.idx,
#              rm_idx = between.idx
#      mvreg : l1_idx/l2_idx = ov.y.idx[[1]]/[[2]], z_idx = between.y.idx,
#              rm_idx = within.x.idx + between.x.idx (+ between.y.idx)

# implied -> 2l: scatter the per-level covariance matrices into the tilde
# universe and gather the y/z blocks
lav_cl_sigma_22l <- function(sigma_w = NULL,
                             sigma_b = NULL,
                             p_tilde = NULL,
                             l1_idx = NULL,
                             l2_idx = NULL,
                             z_idx = NULL,
                             rm_idx = NULL) {
  # scatter into the tilde universe
  sigma_w_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_w_tilde[l1_idx, l1_idx] <- sigma_w
  sigma_b_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_b_tilde[l2_idx, l2_idx] <- sigma_b

  # gather the y/z blocks
  y_idx <- seq_len(p_tilde)
  if (length(rm_idx) > 0L) {
    y_idx <- y_idx[-rm_idx]
  }
  if (length(z_idx) > 0L) {
    sigma_zz <- sigma_b_tilde[z_idx, z_idx, drop = FALSE]
    sigma_yz <- sigma_b_tilde[y_idx, z_idx, drop = FALSE]
  } else {
    sigma_zz <- matrix(0, 0L, 0L)
    sigma_yz <- matrix(0, length(y_idx), 0L)
  }

  list(
    sigma.w = sigma_w_tilde[y_idx, y_idx, drop = FALSE],
    sigma.b = sigma_b_tilde[y_idx, y_idx, drop = FALSE],
    sigma.zz = sigma_zz, sigma.yz = sigma_yz
  )
}

# 2l -> implied: rebuild the level-2 covariance matrix from its y/z blocks
# (the within matrix needs no reassembly: it is the y block itself)
lav_cl_sigma_2l2 <- function(sigma_b = NULL,
                             sigma_zz = NULL,
                             sigma_yz = NULL,
                             p_tilde = NULL,
                             y_idx = NULL,
                             z_idx = NULL,
                             l2_idx = NULL) {
  sigma_b_tilde <- matrix(0, p_tilde, p_tilde)
  sigma_b_tilde[y_idx, y_idx] <- sigma_b
  if (length(z_idx) > 0L) {
    sigma_b_tilde[y_idx, z_idx] <- sigma_yz
    sigma_b_tilde[z_idx, y_idx] <- t(sigma_yz)
    sigma_b_tilde[z_idx, z_idx] <- sigma_zz
  }
  sigma_b_tilde[l2_idx, l2_idx, drop = FALSE]
}

# vech/vec positions (in the tilde vech space) of the per-level and y/z
# blocks; shared by the scores/information rearrangement code below
lav_cl_vech_blocks <- function(p_tilde = NULL,
                               l1_idx = NULL,
                               l2_idx = NULL,
                               z_idx = integer(0L)) {
  b_tilde <- lav_mat_vech_rev(seq_len(p_tilde * (p_tilde + 1) / 2))
  out <- list(
    l1 = lav_mat_vech(b_tilde[l1_idx, l1_idx, drop = FALSE]),
    l2 = lav_mat_vech(b_tilde[l2_idx, l2_idx, drop = FALSE])
  )
  if (length(z_idx) > 0L) {
    out$yz <- lav_mat_vec(b_tilde[l1_idx, z_idx, drop = FALSE])
    out$zz <- lav_mat_vech(b_tilde[z_idx, z_idx, drop = FALSE])
  }
  out
}


# 1) per-cluster-size cache of Sigma.j = nj*Sigma.b.z + Sigma.w quantities
#
#    Sigma.j only depends on the cluster SIZE nj, so the scores engine
#    (looping over clusters) can reuse one inverse per unique size
lav_mvn_cl_sigma_j_cache <- function(cluster_sizes = NULL,
                                     sigma_b_z = NULL,
                                     sigma_w_1 = NULL,
                                     sigma_yz_zi = NULL,
                                     sinv_method = "eigen") {
  lapply(cluster_sizes, function(nj) {
    sigma_j <- (nj * sigma_b_z) + sigma_w_1
    sigma_j_inv <- lav_mat_sym_inverse(
      s = sigma_j,
      logdet = FALSE, sinv_method = sinv_method
    )
    if (!is.null(sigma_yz_zi)) {
      sigma_ji_yz_zi <- sigma_j_inv %*% sigma_yz_zi
      list(
        sigma_j_inv = sigma_j_inv,
        sigma_ji_yz_zi = sigma_ji_yz_zi,
        sigma_zi_zy_ji = t(sigma_ji_yz_zi)
      )
    } else {
      list(sigma_j_inv = sigma_j_inv)
    }
  })
}

# 2) rearrange cluster-wise score blocks (g_muy, g_muz, g_sigma_*) to the
#    Mu.W, vech(Sigma.W), Mu.B, vech(Sigma.B) column layout
#    (shared by the complete-data and missing-data scores)
lav_mvn_cl_scores_2implied <- function(lp = NULL,
                                       g_muy = NULL,
                                       g_muz = NULL,
                                       g_sigma_w = NULL,
                                       g_sigma_b = NULL,
                                       g_sigma_yz = NULL,
                                       g_sigma_zz = NULL) {
  nclusters <- NROW(g_muy)
  between_idx <- lp$between.idx[[2]]
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
  p_tilde_star <- p_tilde * (p_tilde + 1) / 2
  vidx <- lav_cl_vech_blocks(
    p_tilde = p_tilde, l1_idx = ov_idx[[1]], l2_idx = ov_idx[[2]],
    z_idx = between_idx
  )

  sigma_b_tilde <- matrix(0, nclusters, p_tilde_star)
  sigma_b_tilde[, vidx$l1] <- g_sigma_b
  if (length(between_idx) > 0L) {
    sigma_b_tilde[, vidx$yz] <- g_sigma_yz
    sigma_b_tilde[, vidx$zz] <- g_sigma_zz
  }
  sigma_b <- sigma_b_tilde[, vidx$l2, drop = FALSE]

  cbind(mu_w, sigma_w, mu_b, sigma_b)
}

# 3) zero out the rows/cols of a (mu.w, vech(sigma.w), mu.b, vech(sigma.b))
#    information matrix that correspond to fixed exogenous variables
lav_mvn_cl_zero_x_idx <- function(information,
                                  lp = NULL,
                                  mu_w = NULL,
                                  mu_b = NULL,
                                  x_idx = NULL) {
  if (length(x_idx) == 0L) {
    return(information)
  }

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

  information
}

# 4) the per-cluster-size Omega.j loop of the expected information:
#    sum over cluster sizes of n_s * t(Delta.j) I(Omega.j) Delta.j
#    (shared by lav_mvn_cl_info_expected and _info_expected_delta)
lav_mvn_cl_info_omega_j_loop <- function(lp = NULL,
                                         sigma_w_1 = NULL,
                                         sigma_b_1 = NULL,
                                         sigma_zz = NULL,
                                         sigma_yz = NULL,
                                         delta_w_tilde = NULL,
                                         delta_b_tilde = NULL) {
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  n_s <- lp$cluster.size.ns[[2]]
  between_idx <- lp$between.idx[[2]]
  ov_idx <- lp$ov.idx
  p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))
  npar <- NCOL(delta_w_tilde)

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
    } else {
      omega_j <- 1 / nj * sigma_j
    }
    omega_j_inv <- solve(omega_j)

    i11_j <- omega_j_inv
    i22_j <- lav_mvn_kron_dup_half(omega_j_inv)
    i_j <- lav_mat_bdiag(i11_j, i22_j)
    info_j <- t(delta_j) %*% i_j %*% delta_j

    information_j <- information_j + n_s[clz] * info_j
  }

  information_j
}

# 5) numerical observed information: jacobian of an analytic gradient
#    function (shared by the complete-data and missing-data versions)
#
#    dlogl_fn is called with (lp, mu_w, sigma_w, mu_b, sigma_b, sinv_method,
#    return_list = FALSE) plus the extra arguments in dlogl_args, and must
#    return the gradient of -2*logl
lav_mvn_cl_info_obs_engine <- function(lp = NULL,
                                       mu_w = NULL,
                                       sigma_w = NULL,
                                       mu_b = NULL,
                                       sigma_b = NULL,
                                       x_idx = integer(0L),
                                       sinv_method = "eigen",
                                       dlogl_fn = NULL,
                                       dlogl_args = list()) {
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

    dx <- do.call(dlogl_fn, c(
      list(
        lp = lp, mu_w = mu_w2, sigma_w = sigma_w2,
        mu_b = mu_b2, sigma_b = sigma_b2,
        return_list = FALSE, sinv_method = sinv_method
      ),
      dlogl_args
    ))

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
  lav_mvn_cl_zero_x_idx(information,
    lp = lp, mu_w = mu_w, mu_b = mu_b,
    x_idx = x_idx
  )
}


# 6) gradient/scores engine, complete data
#
#    score_mode = FALSE: gradient of -2*logl wrt Mu.W, vech(Sigma.W),
#                        Mu.B, vech(Sigma.B) (loop over cluster sizes)
#    score_mode = TRUE : cluster-wise scores, one row per cluster
lav_mvn_cl_grad_engine <- function(ylp = NULL,
                                   y1 = NULL,
                                   lp = NULL,
                                   mu_w = NULL,
                                   sigma_w = NULL,
                                   mu_b = NULL,
                                   sigma_b = NULL,
                                   score_mode = FALSE,
                                   return_list = FALSE,
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
  cluster_size <- lp$cluster.size[[2]]
  cluster_idx <- lp$cluster.idx[[2]]
  between_idx <- lp$between.idx[[2]]
  cluster_sizes <- lp$cluster.sizes[[2]]
  ncluster_sizes <- lp$ncluster.sizes[[2]]
  cluster_size_ns <- lp$cluster.size.ns[[2]]

  # Y1 (pooled-within sample statistics)
  if (length(between_idx) > 0L) {
    s_pw <- ylp[[2]]$Sigma.W[-between_idx, -between_idx, drop = FALSE]
  } else {
    s_pw <- ylp[[2]]$Sigma.W
  }

  # units: cluster sizes (gradient) or clusters (scores)
  if (score_mode) {
    nunits <- nclusters
    unit_nj <- cluster_size
    unit_clz <- match(cluster_size, cluster_sizes)

    # Y1 within part (centered)
    if (length(between_idx) > 0L) {
      y1w <- y1[, -between_idx, drop = FALSE]
    } else {
      y1w <- y1
    }
    y1w_cm <- t(t(y1w) - mu_y)

    # Y2
    y2 <- ylp[[2]]$Y2
    # NOTE: ORDER mu.b must match Y2
    mu_b_2 <- numeric(ncol(y2))
    if (length(between_idx) > 0L) {
      mu_b_2[-between_idx] <- mu_y
      mu_b_2[between_idx] <- mu_z
    } else {
      mu_b_2 <- mu_y
    }
    y2_cm <- t(t(y2) - mu_b_2)
  } else {
    nunits <- ncluster_sizes
    unit_nj <- cluster_sizes
    unit_clz <- seq_len(ncluster_sizes)

    # Y2 samplestats (aggregated per cluster size)
    cov_d <- ylp[[2]]$cov.d
    mean_d <- ylp[[2]]$mean.d
  }

  # common parts:
  sigma_w_inv <- lav_mat_sym_inverse(
    s = sigma_w_1,
    logdet = FALSE, sinv_method = sinv_method
  )

  # unit-wise gradient blocks
  g_muy <- matrix(0, nunits, length(mu_y))
  g_sigma_w_1 <- matrix(0, nunits, length(lav_mat_vech(sigma_w_1)))
  g_sigma_b_1 <- matrix(0, nunits, length(lav_mat_vech(sigma_b_1)))
  g_muz <- matrix(0, nunits, length(mu_z))
  g_sigma_zz_1 <- matrix(0, nunits, length(lav_mat_vech(sigma_zz)))
  g_sigma_yz_1 <- matrix(0, nunits, length(lav_mat_vec(sigma_yz)))

  if (length(between_idx) > 0L) {
    sigma_zz_inv <- lav_mat_sym_inverse(
      s = sigma_zz,
      logdet = FALSE, sinv_method = sinv_method
    )
    sigma_yz_zi <- sigma_yz %*% sigma_zz_inv
    sigma_zi_zy <- t(sigma_yz_zi)
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy

    # Sigma.j quantities per unique cluster size
    cache <- lav_mvn_cl_sigma_j_cache(
      cluster_sizes = cluster_sizes,
      sigma_b_z = sigma_b_z, sigma_w_1 = sigma_w_1,
      sigma_yz_zi = sigma_yz_zi, sinv_method = sinv_method
    )

    b_idx <- seq_along(between_idx) # z-first layout (gradient mode)
    for (u in seq_len(nunits)) {
      # cluster size
      nj <- unit_nj[u]

      # Sigma.j quantities (cached per cluster size)
      cj <- cache[[unit_clz[u]]]
      sigma_j_inv <- cj$sigma_j_inv
      sigma_ji_yz_zi <- cj$sigma_ji_yz_zi
      sigma_zi_zy_ji <- cj$sigma_zi_zy_ji

      # level-2 vectors + crossproducts
      if (score_mode) {
        yc <- y2_cm[u, -between_idx]
        zc <- y2_cm[u, between_idx]

        y2yc <- tcrossprod(y2_cm[u, ])
        y2yc_zz <- y2yc[between_idx, between_idx, drop = FALSE]
        y2yc_yz <- y2yc[-between_idx, between_idx, drop = FALSE]
        y2yc_yy <- y2yc[-between_idx, -between_idx, drop = FALSE]
      } else {
        zyc <- mean_d[[u]] - c(mu_z, mu_y)
        yc <- zyc[-b_idx]
        zc <- zyc[b_idx]

        y2yc <- (cov_d[[u]] + tcrossprod(mean_d[[u]] - c(mu_z, mu_y)))
        y2yc_zz <- y2yc[b_idx, b_idx, drop = FALSE]
        y2yc_yz <- y2yc[-b_idx, b_idx, drop = FALSE]
        y2yc_yy <- y2yc[-b_idx, -b_idx, drop = FALSE]
      }

      # common parts
      j_yzj <- nj * (sigma_j_inv %*%
        (sigma_yz_zi %*% y2yc_zz %*% t(sigma_yz_zi)
          - y2yc_yz %*% t(sigma_yz_zi)
          - t(y2yc_yz %*% t(sigma_yz_zi)) + y2yc_yy)
        %*% sigma_j_inv)

      z1 <- y2yc_zz %*% t(sigma_ji_yz_zi) %*% sigma_yz
      yz1 <- t(y2yc_yz) %*% sigma_j_inv %*% sigma_yz


      # Mu.Z
      g_muz[u, ] <- -2 * as.numeric(
        (sigma_zz_inv + nj * (sigma_zi_zy_ji %*% sigma_yz_zi)) %*% zc
          - nj * sigma_zi_zy_ji %*% yc
      )

      # MU.Y
      g_muy[u, ] <- 2 * nj * as.numeric(zc %*% sigma_zi_zy_ji -
        yc %*% sigma_j_inv)

      # SIGMA.W
      if (score_mode) {
        # data within for the cluster (centered by mu.y)
        y1m <- y1w_cm[cluster_idx == u, , drop = FALSE]
        g_sigma_w <- ((nj - 1) * sigma_w_inv
          - sigma_w_inv %*% (crossprod(y1m) - nj * y2yc_yy) %*% sigma_w_inv
          + sigma_j_inv - j_yzj)
      } else {
        # between part only; the pooled-within part is added below
        g_sigma_w <- sigma_j_inv - j_yzj
      }
      g_sigma_w_1[u, ] <- lav_mat_vech_dd(g_sigma_w)

      # SIGMA.B
      g_sigma_b <- nj * (sigma_j_inv - j_yzj)
      g_sigma_b_1[u, ] <- lav_mat_vech_dd(g_sigma_b)

      # SIGMA.ZZ
      g_sigma_zz <- (sigma_zz_inv + nj * sigma_zz_inv %*% (
        t(sigma_yz) %*% (sigma_j_inv - j_yzj) %*% sigma_yz
          - (1 / nj * y2yc_zz + t(z1) + z1 - t(yz1) - yz1)) %*%
        sigma_zz_inv)
      g_sigma_zz_1[u, ] <- lav_mat_vech_dd(g_sigma_zz)

      # SIGMA.ZY
      g_sigma_yz <- 2 * nj * (
        (sigma_j_inv %*%
          (sigma_yz_zi %*% y2yc_zz - sigma_yz - y2yc_yz)
          + j_yzj %*% sigma_yz) %*% sigma_zz_inv)

      g_sigma_yz_1[u, ] <- lav_mat_vec(g_sigma_yz)
    }
  # between.idx
  } else { # no level-2 variables

    # Sigma.j quantities per unique cluster size
    cache <- lav_mvn_cl_sigma_j_cache(
      cluster_sizes = cluster_sizes,
      sigma_b_z = sigma_b_1, sigma_w_1 = sigma_w_1,
      sinv_method = sinv_method
    )

    for (u in seq_len(nunits)) {
      # cluster size
      nj <- unit_nj[u]

      # Sigma.j (cached per cluster size)
      sigma_j_inv <- cache[[unit_clz[u]]]$sigma_j_inv

      # level-2 vectors + crossproducts
      if (score_mode) {
        yc <- y2_cm[u, ]
        y2yc_yy <- tcrossprod(y2_cm[u, ])
      } else {
        yc <- mean_d[[u]] - mu_y
        y2yc_yy <- (cov_d[[u]] + tcrossprod(mean_d[[u]] - mu_y))
      }

      # common part
      j_yyj <- nj * sigma_j_inv %*% y2yc_yy %*% sigma_j_inv

      # MU.Y
      g_muy[u, ] <- -2 * nj * as.numeric(yc %*% sigma_j_inv)

      # SIGMA.W
      if (score_mode) {
        # data within for the cluster (centered by mu.y)
        y1m <- y1w_cm[cluster_idx == u, , drop = FALSE]
        g_sigma_w <- ((nj - 1) * sigma_w_inv
          - sigma_w_inv %*% (crossprod(y1m) - nj * y2yc_yy) %*% sigma_w_inv
          + sigma_j_inv - j_yyj)
      } else {
        # between part only; the pooled-within part is added below
        g_sigma_w <- sigma_j_inv - j_yyj
      }
      g_sigma_w_1[u, ] <- lav_mat_vech_dd(g_sigma_w)

      # SIGMA.B
      g_sigma_b <- nj * (sigma_j_inv - j_yyj)
      g_sigma_b_1[u, ] <- lav_mat_vech_dd(g_sigma_b)
    }
  }

  # scores: rearrange columns to Mu.W, Mu.B, Sigma.W, Sigma.B
  if (score_mode) {
    return(lav_mvn_cl_scores_2implied(
      lp = lp,
      g_muy = g_muy, g_muz = g_muz,
      g_sigma_w = g_sigma_w_1, g_sigma_b = g_sigma_b_1,
      g_sigma_yz = g_sigma_yz_1, g_sigma_zz = g_sigma_zz_1
    ))
  }

  # gradient: weighted sum over the cluster-size groups

  # level-1
  d_mu_y <- colSums(g_muy * cluster_size_ns)
  d_sigma_w1 <- lav_mat_vech_rev(colSums(g_sigma_w_1 *
    cluster_size_ns))
  d_sigma_b <- lav_mat_vech_rev(colSums(g_sigma_b_1 *
    cluster_size_ns))

  # level-2
  if (length(between_idx) > 0L) {
    d_mu_z <- colSums(g_muz * cluster_size_ns)
    d_sigma_zz <- lav_mat_vech_rev(colSums(g_sigma_zz_1 *
      cluster_size_ns))
    d_sigma_yz <- matrix(
      colSums(g_sigma_yz_1 * cluster_size_ns),
      nrow(sigma_yz), ncol(sigma_yz)
    )
  } else {
    d_mu_z <- numeric(0L)
    d_sigma_zz <- matrix(0, 0L, 0L)
    d_sigma_yz <- matrix(0, 0L, 0L)
  }

  # Sigma.W (bis): pooled-within part
  d_sigma_w2 <- (lp$nclusters[[1]] - nclusters) * (sigma_w_inv
  - sigma_w_inv %*% s_pw %*% sigma_w_inv)
  tmp <- d_sigma_w2 * 2
  diag(tmp) <- diag(d_sigma_w2)
  d_sigma_w2 <- tmp

  d_sigma_w <- d_sigma_w1 + d_sigma_w2

  # rearrange
  dout <- lav_mvn_cl_2l2implied(
    lp = lp,
    sigma_w = d_sigma_w, sigma_b = d_sigma_b,
    sigma_yz = d_sigma_yz, sigma_zz = d_sigma_zz,
    mu_y = d_mu_y, mu_z = d_mu_z
  )

  if (return_list) {
    dout
  } else {
    c(
      dout$Mu.W, lav_mat_vech(dout$Sigma.W),
      dout$Mu.B, lav_mat_vech(dout$Sigma.B)
    )
  }
}


# 7) gradient/scores engine, missing data (FIML)
#
#    score_mode = FALSE: gradient of -2*logl wrt Mu.W, vech(Sigma.W),
#                        Mu.B, vech(Sigma.B)
#    score_mode = TRUE : cluster-wise scores, one row per cluster
#
#    (the gradient equals the column sums of the cluster-wise scores;
#     both used to be separate, almost identical, functions)
lav_mvn_cl_mi_grad_engine <- function(y1 = NULL,
                                      y2 = NULL,
                                      lp = NULL,
                                      mp = NULL,
                                      mu_w = NULL,
                                      sigma_w = NULL,
                                      mu_b = NULL,
                                      sigma_b = NULL,
                                      score_mode = FALSE,
                                      return_list = FALSE,
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
  sigma_b_1 <- sigma_b_1[both_idx, both_idx, drop = FALSE] # only both part

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
    sigma_zi_zy <- sigma_zz_inv %*% sigma_zy
    sigma_b_z <- sigma_b_1 - sigma_yz %*% sigma_zi_zy
    gz <- z_c %*% sigma_zz_inv # for complete cases only
  }

  # containters per cluster
  alist_1 <- rep(list(matrix(0, ny, ny)), nclusters)

  # cluster-wise gradient blocks
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
      z_na_idx <- which(!zp_1$pat[p, ])
      j_idx <- zp_1$case.idx[[p]] # cluster indices with this pattern
      zpat2j[j_idx] <- p

      if (length(z_na_idx) > 0L) {
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
          alist_1[[j]] + (wip[[p]] * tab[j])
      }
    } else {
      # complete case
      pij_1[mp$case.idx[[p]], ] <-
        y1w_c[mp$case.idx[[p]], ] %*% sigma_w_inv
      wip[[p]] <- sigma_w_inv
      for (j in j1_idx) {
        alist_1[[j]] <-
          alist_1[[j]] + (sigma_w_inv * tab[j])
      }
    }
  } # p

  # fully-missing-within cases (empty.idx) belong to a cluster but
  # contribute nothing to the loglik; flag them (-1L) so they are not
  # counted as complete cases when computing dx.sigma.w below
  mpi[mp$empty.idx] <- -1L

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

    # one inverse per cluster (used to be 3 solve.default() calls)
    ibza_j_inv <- solve.default(ibza_j)

    ibza_j_inv_bz <- ibza_j_inv %*% sigma_b_z
    ibza_j_inv_bz_p <- ibza_j_inv_bz %*% pb_j
    a_ibza_j_inv_bz <- a_j %*% ibza_j_inv_bz
    a_ibza_j_inv_bz_p <- a_ibza_j_inv_bz %*% pb_j

    a_ibza_j_inv <- a_j %*% ibza_j_inv
    p_ibza_j_inv <- t(crossprod(pb_j, ibza_j_inv))

    # only if we have between-only variables
    if (nz > 0L) {
      g_j <- as.matrix(gj[j, ])
      zij <- as.matrix(gz[j, ])
      zizy <- zizy_1[[zpat2j[j]]]
      zip <- zip_1[[zpat2j[j]]]

      ibza_j_inv_zizy <- ibza_j_inv %*% t(zizy)
      ibza_j_inv_g <- ibza_j_inv %*% g_j
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
      g_sigma_zz[j, ] <- lav_mat_vech_dd(tmp)

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
    which_incompl <- which(mpi[ij_index] > 0L)

    ap2 <- rep(list(sigma_w_inv %*% part2), length(ij_index))
    ap1a_a <- ap1a_b <- matrix(0, ny, ny)
    if (length(which_compl) > 0L) {
      tmp <- (sigma_w_inv %*% part1 %*% sigma_w_inv)
      ap1a_a <- tmp * length(which_compl)
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
    }
    t1 <- ap1a_a + ap1a_b
    t2 <- (do.call("cbind", ap2) - t(pij)) %*% pij

    aa_wj <- t1 + t2

    tmp <- a_j_full + (aa_wj + t(aa_wj)) / 2
    # symmetry correction
    g_sigma_w[j, ] <- lav_mat_vech_dd(tmp)

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
    g_muy[j, ] <- drop(2 * gbzpp)
  } # j

  # scores: rearrange columns to Mu.W, Mu.B, Sigma.W, Sigma.B
  if (score_mode) {
    return(lav_mvn_cl_scores_2implied(
      lp = lp,
      g_muy = g_muy, g_muz = g_muz,
      g_sigma_w = g_sigma_w, g_sigma_b = g_sigma_b,
      g_sigma_yz = g_sigma_yz, g_sigma_zz = g_sigma_zz
    ))
  }

  # gradient: sum over the clusters
  dx_mu_y <- colSums(g_muy)
  dx_sigma_w <- lav_mat_vech_rev(colSums(g_sigma_w))
  dx_sigma_b <- lav_mat_vech_rev(colSums(g_sigma_b))
  if (nz > 0L) {
    dx_mu_z <- colSums(g_muz)
    dx_sigma_zz <- lav_mat_vech_rev(colSums(g_sigma_zz))
    dx_sigma_yz <- matrix(
      colSums(g_sigma_yz),
      nrow(out$sigma.yz), ncol(out$sigma.yz)
    )
  } else {
    dx_mu_z <- numeric(0L)
    dx_sigma_zz <- matrix(0, 0L, 0L)
    dx_sigma_yz <- matrix(0, nrow(out$sigma.yz), 0L)
  }

  # rearrange
  dout <- lav_mvn_cl_2l2implied(
    lp = lp,
    sigma_w = dx_sigma_w, sigma_b = dx_sigma_b,
    sigma_yz = dx_sigma_yz, sigma_zz = dx_sigma_zz,
    mu_y = dx_mu_y, mu_z = dx_mu_z
  )

  if (return_list) {
    dout
  } else {
    c(
      dout$Mu.W, lav_mat_vech(dout$Sigma.W),
      dout$Mu.B, lav_mat_vech(dout$Sigma.B)
    )
  }
}
