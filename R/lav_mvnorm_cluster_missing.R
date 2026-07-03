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
  sigma_b_1 <- sigma_b_1[both_idx, both_idx, drop = FALSE] # only both part

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
  lav_mvn_cl_mi_grad_engine(
    y1 = y1, y2 = y2, lp = lp, mp = mp,
    mu_w = mu_w, sigma_w = sigma_w,
    mu_b = mu_b, sigma_b = sigma_b,
    score_mode = FALSE, return_list = return_list,
    sinv_method = sinv_method
  )
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
  lav_mvn_cl_mi_grad_engine(
    y1 = y1, y2 = y2, lp = lp, mp = mp,
    mu_w = mu_w, sigma_w = sigma_w,
    mu_b = mu_b, sigma_b = sigma_b,
    score_mode = TRUE,
    sinv_method = sinv_method
  )
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
  information <- lav_mvn_cl_zero_x_idx(information,
    lp = lp, mu_w = mu_w, mu_b = mu_b, x_idx = x_idx
  )

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
  lav_mvn_cl_info_obs_engine(
    lp = lp,
    mu_w = mu_w, sigma_w = sigma_w,
    mu_b = mu_b, sigma_b = sigma_b,
    x_idx = x_idx, sinv_method = sinv_method,
    dlogl_fn = lav_mvn_cl_mi_dlogl_2l_samp,
    dlogl_args = list(y1 = y1, y2 = y2, mp = mp)
  )
}
