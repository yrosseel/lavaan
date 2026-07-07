# loglikelihood clustered/twolevel data in the presence of missing data

# YR:
# - objective function: first version around March 2021 (see Psych paper)
# - analytic gradient: first version around May 2021

# NOTE: see the top of lav_mvnorm_cluster.R for the 'implied' vs '2l'
# representations and the invariants of the conversion functions


# Mu.W, Mu.B, Sigma.W, Sigma.B are the model-implied statistics; a
# pre-converted 2l list (as returned by lav_mvn_cl_implied22l) may be
# passed as out= instead, to avoid re-converting
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
                                         minus_two = TRUE,
                                         out = NULL) {
  # map implied to 2l matrices (unless the caller already did)
  if (is.null(out)) {
    out <- lav_mvn_cl_implied22l(
      lp = lp, mu_w = mu_w, mu_b = mu_b,
      sigma_w = sigma_w, sigma_b = sigma_b
    )
  }
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
  for (p in seq_len(mp$npatterns)) {
    freq <- mp$freq[p]
    na_idx <- which(!mp$pat[p, ])
    j_idx <- mp$j.idx[[p]]
    j1_idx <- mp$j1.idx[[p]]
    tab <- integer(nclusters)
    tab[j1_idx] <- mp$j.freq[[p]]

    # compute sigma.w.inv for this pattern
    if (length(na_idx) > 0L) {
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

# estimate ML estimates of Mu.W, Mu.B, Sigma.W, Sigma.B in the presence
# of missing data, using the EM algorithm (new in 0.7-1)
#
# this extends lav_mvn_cl_em_sat() to the missing data setting: in the
# E-step, we integrate over both the (latent) cluster random effects and
# the missing values; this mimics the approach used by Mplus (H1
# estimation with TYPE = TWOLEVEL and missing data)
#
# notation (same '2l' parameterization as the loglik function above):
# - y_ij = mu.y + u_j + e_ij, where mu.y = mu.w + mu.b,
#   u_j ~ N(0, sigma.b) (only the 'both' part has nonzero variance), and
#   e_ij ~ N(0, sigma.w)
# - z_j are the between-only variables: z_j = mu.z + (z-part), with
#   Var(z) = sigma.zz and Cov(u, z) = sigma.yz (only the 'both' rows)
#
# E-step, per cluster j (writing beta = u_j['both'] and centering all
# observed data):
# - the observed units contribute the precision A.j = sum_i T_i' W_i^{-1} T_i
#   and the information vector p.j = sum_i T_i' W_i^{-1} (y_i.obs - mu.y),
#   where W_i is the observed block of sigma.w for the pattern of unit i
#   (these are the same quantities used in the loglik function above)
# - the prior of w = (beta, z.mis) given z.obs is N(m0, C0) (standard
#   conditioning of the joint normal); the posterior is then
#       m1 = (I + C0 %*% Atilde)^{-1} (m0 + C0 %*% ptilde)
#       C1 = (I + C0 %*% Atilde)^{-1} C0
#   with Atilde = blkdiag(A.j, 0) and ptilde = c(p.j, 0); this form remains
#   valid even if C0 is singular
# - the missing y-values of unit i, given (beta, y_i.obs), follow the usual
#   within-pattern regression: mean B_p (y.obs - mu.y[obs] - beta[obs]),
#   covariance C_p
#
# the M-step is the standard saturated update based on the expected
# complete-data sufficient statistics
#
# the E-step machinery is also used by the h0 EM optimizer
# (lav_mvn_cl_em_h0, when missing = "ml"); it is therefore exposed as an
# 'engine': a list of closures over the data and the missing patterns:
# - theta0:  packed (available-case) starting values
# - pack():  pack the internal '2l' matrices into a parameter vector
# - step():  one EM step (E-step + saturated M-step update)
# - logl():  observed-data loglikelihood at theta
# - implied(): the implied-space matrices (Sigma.W, Sigma.B, Mu.W, Mu.B)
lav_mvn_cl_mi_em_engine <- function(y1 = NULL,
                                    y2 = NULL,
                                    lp = NULL,
                                    mp = NULL,
                                    loglik_x = 0,
                                    min_variance = 1e-05) {
  if (is.null(loglik_x)) {
    loglik_x <- 0
  }

  # dimensions and indices
  nobs <- lp$nclusters[[1]]
  nclusters <- lp$nclusters[[2]]
  cluster_idx <- lp$cluster.idx[[2]]
  between_idx <- lp$between.idx[[2]]
  both_idx <- lp$both.idx[[2]] # between-only variables always come last,
                               # so these indices are also valid in y-space
  nz <- length(between_idx)
  nb <- length(both_idx)
  zp <- mp$Zp

  # y data (level-1 variables only) and z data (one row per cluster)
  if (nz > 0L) {
    y1w <- y1[, -between_idx, drop = FALSE]
    z <- y2[, between_idx, drop = FALSE]
  } else {
    y1w <- y1
    z <- NULL
  }
  ny <- ncol(y1w)

  # starting values (available-case)
  mu_y0 <- colMeans(y1w, na.rm = TRUE)
  mu_y0[!is.finite(mu_y0)] <- 0
  var_y0 <- apply(y1w, 2L, var, na.rm = TRUE)
  var_y0[!is.finite(var_y0) | var_y0 < min_variance] <- 1
  mu_w <- mu_y0
  mu_w[both_idx] <- 0
  mu_b <- numeric(ny)
  mu_b[both_idx] <- mu_y0[both_idx]
  sigma_w <- diag(var_y0, nrow = ny)
  sigma_w[both_idx, both_idx] <- diag(var_y0[both_idx] / 2,
                                      nrow = nb)
  sigma_b <- matrix(0, ny, ny)
  sigma_b[both_idx, both_idx] <- diag(var_y0[both_idx] / 2, nrow = nb)
  if (nz > 0L) {
    mu_z <- colMeans(z, na.rm = TRUE)
    mu_z[!is.finite(mu_z)] <- 0
    var_z0 <- apply(z, 2L, var, na.rm = TRUE)
    var_z0[!is.finite(var_z0) | var_z0 < min_variance] <- 1
    sigma_zz <- diag(var_z0, nrow = nz)
    sigma_yz <- matrix(0, ny, nz)
  } else {
    mu_z <- numeric(0L)
    sigma_zz <- matrix(0, 0L, 0L)
    sigma_yz <- matrix(0, ny, 0L)
  }

  # z-values per cluster: use the between-level missing patterns (Zp)
  # to decide which entries are observed
  if (nz > 0L) {
    zpat2j <- integer(nclusters)
    for (q in seq_len(zp$npatterns)) {
      zpat2j[zp$case.idx[[q]]] <- q
    }
    if (length(zp$empty.idx) > 0L) {
      zpat2j[zp$empty.idx] <- zp$npatterns + 1L
    }
  }

  # constants for the fused loglikelihood (a byproduct of the E-step;
  # see em_step below)
  p_1 <- mp$nel
  if (nz > 0L) {
    p_1 <- p_1 + zp$nel
  }
  has_x <- length(unlist(lp$ov.x.idx)) > 0L
  log_2pi <- log(2 * pi)

  # pack/unpack the internal parameter vector (structural zeros stay zero
  # under linear extrapolation, so we simply pack the full matrices)
  em_pack <- function(mu_w, mu_b, sigma_w, sigma_b, mu_z, sigma_zz,
                      sigma_yz) {
    c(mu_w, mu_b, lav_mat_vech(sigma_w), lav_mat_vech(sigma_b), mu_z,
      lav_mat_vech(sigma_zz), as.vector(sigma_yz))
  }
  em_unpack <- function(theta) {
    i <- 0L
    mu_w <- theta[i + seq_len(ny)]; i <- i + ny
    mu_b <- theta[i + seq_len(ny)]; i <- i + ny
    sigma_w <- lav_mat_vech_rev(theta[i + seq_len(ny * (ny + 1) / 2)])
    i <- i + ny * (ny + 1) / 2
    sigma_b <- lav_mat_vech_rev(theta[i + seq_len(ny * (ny + 1) / 2)])
    i <- i + ny * (ny + 1) / 2
    mu_z <- theta[i + seq_len(nz)]; i <- i + nz
    sigma_zz <- if (nz > 0L) {
      lav_mat_vech_rev(theta[i + seq_len(nz * (nz + 1) / 2)])
    } else {
      matrix(0, 0L, 0L)
    }
    i <- i + nz * (nz + 1) / 2
    sigma_yz <- matrix(theta[i + seq_len(ny * nz)], ny, nz)
    list(mu_w = mu_w, mu_b = mu_b, sigma_w = sigma_w, sigma_b = sigma_b,
         mu_z = mu_z, sigma_zz = sigma_zz, sigma_yz = sigma_yz)
  }

  # one EM step: theta -> theta'; if logl = TRUE, the observed-data
  # loglikelihood at the *input* theta -- a cheap byproduct of the E-step
  # quantities -- is attached to the result as the "logl" attribute (it
  # matches the value of lav_mvn_cl_mi_loglik_samp_2l() exactly)
  em_step <- function(theta, logl = FALSE) {
    th <- em_unpack(theta)
    mu_w <- th$mu_w
    mu_b <- th$mu_b
    sigma_w <- th$sigma_w
    sigma_b <- th$sigma_b
    mu_z <- th$mu_z
    sigma_zz <- th$sigma_zz
    sigma_yz <- th$sigma_yz

    # current 'both' blocks
    sb <- sigma_b[both_idx, both_idx, drop = FALSE]
    if (nz > 0L) {
      syz <- sigma_yz[both_idx, , drop = FALSE]
    }
    mu_y <- mu_w + mu_b
    y1w_c <- t(t(y1w) - mu_y)

    sigma_w_inv <- try(solve.default(sigma_w), silent = TRUE)
    if (inherits(sigma_w_inv, "try-error")) {
      lav_msg_stop(gettext(
        "EM steps of the saturated (H1) model failed; some matrices may
         be singular; please check your data for (near-)perfect
         correlations."))
    }
    sigma_w_logdet <- 0
    if (logl) {
      sigma_w_logdet <- c(determinant.matrix(sigma_w,
                                             logarithm = TRUE)$modulus)
      # accumulators for the fused loglikelihood
      w_logdet <- zz_logdet <- ibza_logdet <- 0
      q_yy_a <- q_zz_a <- q_corr <- 0
    }

    # 1. per-pattern quantities + A.j and p.j per cluster
    winv_p <- breg_p <- ccond_p <- vector("list", mp$npatterns)
    alist_1 <- rep(list(matrix(0, nb, nb)), nclusters)
    pij <- matrix(0, nrow(y1w), ny)
    for (p in seq_len(mp$npatterns)) {
      o_idx <- which(mp$pat[p, ])
      na_idx <- which(!mp$pat[p, ])
      case_idx <- mp$case.idx[[p]]

      if (length(na_idx) > 0L) {
        wp_inv <- lav_mat_sym_inverse_update(
          s_inv = sigma_w_inv, rm_idx = na_idx,
          logdet = logl, s_logdet = sigma_w_logdet
        )
        if (logl) {
          w_logdet <- w_logdet + (attr(wp_inv, "logdet") * mp$freq[p])
        }
        breg_p[[p]] <- sigma_w[na_idx, o_idx, drop = FALSE] %*% wp_inv
        ccond_p[[p]] <- sigma_w[na_idx, na_idx, drop = FALSE] -
          breg_p[[p]] %*% sigma_w[o_idx, na_idx, drop = FALSE]
      } else {
        wp_inv <- sigma_w_inv
        if (logl) {
          w_logdet <- w_logdet + (sigma_w_logdet * mp$freq[p])
        }
      }
      winv_p[[p]] <- wp_inv

      pij[case_idx, o_idx] <-
        y1w_c[case_idx, o_idx, drop = FALSE] %*% wp_inv

      # contribution to A.j ('both' block only)
      a_full <- matrix(0, ny, ny)
      a_full[o_idx, o_idx] <- wp_inv
      a_bb <- a_full[both_idx, both_idx, drop = FALSE]
      j1_idx <- mp$j1.idx[[p]]
      j_freq <- mp$j.freq[[p]]
      for (k in seq_along(j1_idx)) {
        alist_1[[j1_idx[k]]] <- alist_1[[j1_idx[k]]] + (a_bb * j_freq[k])
      }
    } # patterns
    pj <- rowsum.default(pij[, both_idx, drop = FALSE], cluster_idx,
      reorder = FALSE, na.rm = TRUE
    )
    if (logl) {
      q_yy_a <- sum(pij * y1w_c, na.rm = TRUE)
    }

    # 2. prior of (beta, z.mis) given z.obs, per z-pattern
    if (nz > 0L) {
      zc <- t(t(z) - mu_z)
      nqz <- zp$npatterns + (length(zp$empty.idx) > 0L)
      c0_q <- m0coef_q <- zmis_q <- vector("list", nqz)
      for (q in seq_len(nqz)) {
        if (q > zp$npatterns) { # empty pattern: no z observed
          zo_idx <- integer(0L)
          zm_idx <- seq_len(nz)
        } else {
          zo_idx <- which(zp$pat[q, ])
          zm_idx <- which(!zp$pat[q, ])
        }
        zmis_q[[q]] <- zm_idx
        k_q <- rbind(
          cbind(sb, syz[, zm_idx, drop = FALSE]),
          cbind(
            t(syz[, zm_idx, drop = FALSE]),
            sigma_zz[zm_idx, zm_idx, drop = FALSE]
          )
        )
        if (length(zo_idx) > 0L) {
          l_q <- rbind(
            syz[, zo_idx, drop = FALSE],
            sigma_zz[zm_idx, zo_idx, drop = FALSE]
          )
          szz_o <- sigma_zz[zo_idx, zo_idx, drop = FALSE]
          soo_inv <- solve.default(szz_o)
          m0coef_q[[q]] <- l_q %*% soo_inv
          c0_q[[q]] <- k_q - m0coef_q[[q]] %*% t(l_q)

          if (logl) {
            # fused loglikelihood: marginal of the observed z-values
            zz_logdet <- zz_logdet +
              (c(determinant.matrix(szz_o, logarithm = TRUE)$modulus) *
               zp$freq[q])
            zc_o <- zc[zp$case.idx[[q]], zo_idx, drop = FALSE]
            q_zz_a <- q_zz_a + sum((zc_o %*% soo_inv) * zc_o)
          }
        } else {
          m0coef_q[[q]] <- matrix(0, nb + length(zm_idx), 0L)
          c0_q[[q]] <- k_q
        }
      }
    }

    # 3. per-cluster posterior of (beta, z.mis); accumulate the
    #    between-level (and z) expected sufficient statistics
    eb <- matrix(0, nclusters, nb) # E(beta | data)
    vb_list <- vector("list", nclusters) # Cov(beta | data)
    tb1 <- numeric(nb)
    tb2 <- matrix(0, nb, nb)
    if (nz > 0L) {
      tz1 <- numeric(nz)
      tz2 <- matrix(0, nz, nz)
      tzv <- matrix(0, nz, nb)
    }
    for (j in seq_len(nclusters)) {
      a_j <- alist_1[[j]]
      p_j <- pj[j, ]
      if (nz > 0L) {
        q <- zpat2j[j]
        c0 <- c0_q[[q]]
        zm_idx <- zmis_q[[q]]
        nzm <- length(zm_idx)
        if (q > zp$npatterns) {
          m0 <- numeric(nb + nzm)
        } else {
          zo_idx <- which(zp$pat[q, ])
          m0 <- drop(m0coef_q[[q]] %*% zc[j, zo_idx])
        }
      } else {
        c0 <- sb
        nzm <- 0L
        m0 <- numeric(nb)
      }

      # posterior of w = (beta, z.mis)
      nw <- nb + nzm
      m_j <- c0[, seq_len(nb), drop = FALSE] %*% a_j # C0 %*% Atilde
      m_j <- cbind(m_j, matrix(0, nw, nzm))
      m_j[lav_mat_diag_idx(nw)] <- m_j[lav_mat_diag_idx(nw)] + 1
      rhs <- cbind(m0 + drop(c0[, seq_len(nb), drop = FALSE] %*% p_j), c0)
      sol <- try(solve.default(m_j, rhs), silent = TRUE)
      if (inherits(sol, "try-error")) {
        lav_msg_stop(gettext(
          "EM steps of the saturated (H1) model failed; some matrices may
           be singular; please check your data for (near-)perfect
           correlations."))
      }
      m1 <- sol[, 1L]
      c1 <- sol[, -1L, drop = FALSE]
      c1 <- (c1 + t(c1)) / 2

      ebeta <- m1[seq_len(nb)]
      vb <- c1[seq_len(nb), seq_len(nb), drop = FALSE]
      eb[j, ] <- ebeta
      vb_list[[j]] <- vb

      if (logl) {
        # fused loglikelihood: m_j is block lower-triangular, so
        # det(m_j) = det(I + C0[beta, beta] %*% A.j); the second term is
        # the correction from integrating out w = (beta, z.mis)
        ldm <- determinant.matrix(m_j, logarithm = TRUE)
        ibza_logdet <- ibza_logdet + (c(ldm$modulus) * ldm$sign)
        m0_b <- m0[seq_len(nb)]
        q_corr <- q_corr + sum(p_j * (m0_b + ebeta)) -
          sum(m0_b * (a_j %*% ebeta))
      }

      # E(u_j | data), uncentered ('both' part)
      ev_b <- mu_b[both_idx] + ebeta
      tb1 <- tb1 + ev_b
      tb2 <- tb2 + vb + tcrossprod(ev_b)

      if (nz > 0L) {
        ez <- z[j, ]
        if (nzm > 0L) {
          ez[zm_idx] <- mu_z[zm_idx] + m1[nb + seq_len(nzm)]
        }
        tz1 <- tz1 + ez
        ezz <- tcrossprod(ez)
        ezv <- tcrossprod(ez, ev_b)
        if (nzm > 0L) {
          ezz[zm_idx, zm_idx] <- ezz[zm_idx, zm_idx, drop = FALSE] +
            c1[nb + seq_len(nzm), nb + seq_len(nzm), drop = FALSE]
          ezv[zm_idx, ] <- ezv[zm_idx, , drop = FALSE] +
            c1[nb + seq_len(nzm), seq_len(nb), drop = FALSE]
        }
        tz2 <- tz2 + ezz
        tzv <- tzv + ezv
      }
    } # clusters

    # 4. within-level expected sufficient statistics
    #    d_i = y_i - u_j; E(d_i) and E(d_i d_i') per pattern
    tw1 <- numeric(ny)
    tw2 <- matrix(0, ny, ny)
    eb_full <- matrix(0, nclusters, ny)
    eb_full[, both_idx] <- eb
    for (p in seq_len(mp$npatterns)) {
      o_idx <- which(mp$pat[p, ])
      na_idx <- which(!mp$pat[p, ])
      case_idx <- mp$case.idx[[p]]
      j_idx <- mp$j.idx[[p]]

      # E(d_i): observed part = y - mu.b - E(beta); missing part via B_p
      f_o <- y1w_c[case_idx, o_idx, drop = FALSE] -
        eb_full[j_idx, o_idx, drop = FALSE]
      e_full <- matrix(0, length(case_idx), ny)
      e_full[, o_idx] <- t(t(f_o) + mu_w[o_idx])
      if (length(na_idx) > 0L) {
        e_full[, na_idx] <- t(t(f_o %*% t(breg_p[[p]])) + mu_w[na_idx])
      }
      tw1 <- tw1 + colSums(e_full)
      tw2 <- tw2 + crossprod(e_full)

      # covariance contributions, per cluster in this pattern
      ob <- which(o_idx %in% both_idx) # both-var positions within o_idx
      b_pos <- match(o_idx[ob], both_idx)
      j1_idx <- mp$j1.idx[[p]]
      j_freq <- mp$j.freq[[p]]
      b_p <- breg_p[[p]]
      for (k in seq_along(j1_idx)) {
        v_oo <- matrix(0, length(o_idx), length(o_idx))
        if (length(ob) > 0L) {
          v_oo[ob, ob] <- vb_list[[j1_idx[k]]][b_pos, b_pos, drop = FALSE]
        }
        m_cov <- matrix(0, ny, ny)
        m_cov[o_idx, o_idx] <- v_oo
        if (length(na_idx) > 0L) {
          bv <- b_p %*% v_oo
          m_cov[na_idx, o_idx] <- bv
          m_cov[o_idx, na_idx] <- t(bv)
          m_cov[na_idx, na_idx] <- tcrossprod(bv, b_p) + ccond_p[[p]]
        }
        tw2 <- tw2 + (j_freq[k] * m_cov)
      }
    } # patterns

    # empty units (all level-1 variables missing): E(d) = mu.w,
    # Cov(d) = sigma.w -- no data contribution
    n_empty <- length(mp$empty.idx)
    if (n_empty > 0L) {
      tw1 <- tw1 + (n_empty * mu_w)
      tw2 <- tw2 + n_empty * (sigma_w + tcrossprod(mu_w))
    }

    # 5. M-step (saturated updates)
    mu_w <- tw1 / nobs
    sigma_w <- tw2 / nobs - tcrossprod(mu_w)
    sigma_w <- (sigma_w + t(sigma_w)) / 2
    m_b <- tb1 / nclusters
    mu_b <- numeric(ny)
    mu_b[both_idx] <- m_b
    sigma_b <- matrix(0, ny, ny)
    sigma_b[both_idx, both_idx] <- tb2 / nclusters - tcrossprod(m_b)
    if (nz > 0L) {
      mu_z <- tz1 / nclusters
      sigma_zz <- tz2 / nclusters - tcrossprod(mu_z)
      sigma_zz <- (sigma_zz + t(sigma_zz)) / 2
      sigma_yz <- matrix(0, ny, nz)
      sigma_yz[both_idx, ] <- t(tzv / nclusters - tcrossprod(mu_z, m_b))
    }

    # check for (near-zero) variances at the within level, and set
    # them to min.variance
    zero_var <- which(diag(sigma_w) < min_variance)
    if (length(zero_var) > 0L) {
      sigma_w[, zero_var] <- 0
      sigma_w[zero_var, ] <- 0
      diag(sigma_w)[zero_var] <- min_variance
    }

    theta_new <- em_pack(mu_w, mu_b, sigma_w, sigma_b, mu_z, sigma_zz,
                         sigma_yz)
    if (logl) {
      # assemble the loglikelihood at the input theta (note: all the
      # accumulators were computed *before* the M-step overwrote the
      # parameter matrices)
      fx_in <- -(p_1 * log_2pi + (w_logdet + ibza_logdet + zz_logdet) +
                 (q_yy_a + q_zz_a - q_corr)) / 2
      if (has_x) {
        fx_in <- fx_in - loglik_x
      }
      attr(theta_new, "logl") <- fx_in
    }
    theta_new
  } # em_step

  # loglikelihood at theta -- the state is already in 2l form, so we can
  # skip the 2l -> implied -> 2l round trip (mu.y = mu.w + mu.b, with
  # mu.b[within-only] structurally zero; see the invariants at the top of
  # lav_mvnorm_cluster.R)
  em_logl <- function(theta) {
    th <- em_unpack(theta)
    out <- list(
      sigma.w = th$sigma_w, sigma.b = th$sigma_b,
      sigma.zz = th$sigma_zz, sigma.yz = th$sigma_yz,
      mu.z = th$mu_z, mu.y = th$mu_w + th$mu_b,
      mu.w = th$mu_w, mu.b = th$mu_b
    )
    lav_mvn_cl_mi_loglik_samp_2l(
      y1 = y1, y2 = y2, lp = lp, mp = mp,
      loglik_x = loglik_x, log2pi = TRUE, minus_two = FALSE,
      out = out
    )
  }

  # implied-space matrices at theta
  em_implied <- function(theta) {
    th <- em_unpack(theta)
    lav_mvn_cl_2l2implied(
      lp = lp, sigma_w = th$sigma_w, sigma_b = th$sigma_b,
      sigma_zz = th$sigma_zz, sigma_yz = th$sigma_yz, mu_z = th$mu_z,
      mu_y = NULL, mu_w = th$mu_w, mu_b = th$mu_b
    )
  }

  list(
    theta0 = em_pack(mu_w, mu_b, sigma_w, sigma_b, mu_z, sigma_zz,
                     sigma_yz),
    pack = em_pack,
    step = em_step,
    step_logl = function(theta) em_step(theta, logl = TRUE),
    logl = em_logl,
    implied = em_implied
  )
}

lav_mvn_cl_mi_em_sat <- function(y1 = NULL,
                                 y2 = NULL,
                                 lp = NULL,
                                 mp = NULL,
                                 loglik_x = 0,
                                 tol = 1e-04, # = em.h1.args$tol
                                 max_iter = 5000L, # = em.h1.args$max_iter
                                 min_variance = 1e-05,
                                 acceleration = "none",
                                 fused = TRUE) { # = em.h1.args$fused
  if (is.null(fused)) {
    fused <- TRUE
  }
  engine <- lav_mvn_cl_mi_em_engine(
    y1 = y1, y2 = y2, lp = lp, mp = mp,
    loglik_x = loglik_x, min_variance = min_variance
  )

  theta <- engine$theta0

  # EM iterations
  if (acceleration %in% c("squarem", "qn")) {
    # report initial fx
    fx <- engine$logl(theta)
    if (lav_verbose()) {
      cat(
        "EM iter:", sprintf("%3d", 0),
        " fx =", sprintf("%17.10f", fx),
        "\n"
      )
    }
    accel_fn <- if (acceleration == "qn") lav_em_qn else lav_em_squarem
    out_em <- accel_fn(
      theta = theta, step_fn = engine$step, logl_fn = engine$logl,
      fx0 = fx, tol = tol, max_iter = max_iter
    )
    theta <- out_em$theta
    fx <- out_em$fx
    converged <- out_em$converged
  } else {
    # plain EM iterations; if fused = TRUE, the loglikelihood at theta is
    # obtained as a cheap byproduct of engine$step_logl(theta) (returned
    # as the "logl" attribute), so no separate loglikelihood evaluations
    # are needed; we merely look one EM step ahead, and at the end
    # (theta, fx) still form an exact pair; if fused = FALSE, the
    # stand-alone loglikelihood function is used instead (same iterates,
    # same stopping rule, same result)
    if (fused) {
      theta_next <- engine$step_logl(theta)
      fx <- attr(theta_next, "logl") # logl at theta0
    } else {
      fx <- engine$logl(theta)
    }
    if (lav_verbose()) {
      cat(
        "EM iter:", sprintf("%3d", 0),
        " fx =", sprintf("%17.10f", fx),
        "\n"
      )
    }
    fx_old <- fx
    converged <- FALSE
    for (i in seq_len(max_iter)) {
      if (fused) {
        theta <- theta_next
        theta_next <- engine$step_logl(theta)
        fx <- attr(theta_next, "logl") # logl at the current theta
      } else {
        theta <- engine$step(theta)
        fx <- engine$logl(theta)
      }

      # fx.delta
      fx_delta <- fx - fx_old

      # check if fx.delta is finite
      if (!is.finite(fx_delta)) {
        implied2 <- engine$implied(theta)
        cat("\n")
        cat("FATAL problem: dumping Sigma.W and Sigma.B matrices:\n\n")
        cat("Sigma.W:\n")
        print(implied2$Sigma.W)
        cat("\n")
        cat("Sigma.B:\n")
        print(implied2$Sigma.B)
        cat("\n")
        lav_msg_stop(gettext(
          "EM steps of the saturated (H1) model failed; some matrices may
           be singular; please check your data for (near-)perfect
           correlations."))
      }

      # what if fx.delta is negative?
      if (fx_delta < -1e-06) {
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
      if (abs(fx_delta) < tol) {
        converged <- TRUE
        break
      } else {
        fx_old <- fx
      }
    }
  } # EM iterations

  # warning?
  if (!converged) {
    lav_msg_warn(gettext(
      "Maximum number of iterations reached when computing the sample
       moments of the saturated (H1) model using EM; increase the max_iter
       element of the em.h1.args= argument to increase the number of
       iterations"))
  }

  # final implied matrices
  implied2 <- engine$implied(theta)

  list(
    Sigma.W = implied2$Sigma.W, Sigma.B = implied2$Sigma.B,
    Mu.W = implied2$Mu.W, Mu.B = implied2$Mu.B, logl = fx
  )
}

# the posterior distribution of the cluster-level unknowns, given the
# observed data and the model-implied statistics (Mu.W, Sigma.W, Mu.B,
# Sigma.B): per cluster j, w_j = (beta_j, ztilde_j[mis]) has posterior
# N(m1_j, C1_j) (see the E-step of lav_mvn_cl_mi_em_engine); also
# returns the per-pattern sweep regressions (breg, ccond) needed to
# impute the missing level-1 values, and various intermediate
# quantities (used by lav_mvn_cl_mi_scov_louis and
# lav_mvn_cl_mi_estep_ranef)
lav_mvn_cl_mi_posterior <- function(y1 = NULL,
                                    y2 = NULL,
                                    lp = NULL,
                                    mp = NULL,
                                    mu_w = NULL, # implied Mu.W
                                    sigma_w = NULL, # implied Sigma.W
                                    mu_b = NULL, # implied Mu.B
                                    sigma_b = NULL, # implied Sigma.B
                                    out = NULL) { # pre-converted 2l list
  # map implied to 2l matrices (unless the caller already did)
  if (is.null(out)) {
    out <- lav_mvn_cl_implied22l(
      lp = lp, mu_w = mu_w, mu_b = mu_b,
      sigma_w = sigma_w, sigma_b = sigma_b
    )
  }
  mu_y <- out$mu.y
  mu_z <- out$mu.z
  sw <- out$sigma.w
  sb <- out$sigma.b[lp$both.idx[[2]], lp$both.idx[[2]], drop = FALSE]
  sigma_zz <- out$sigma.zz
  syz <- out$sigma.yz[lp$both.idx[[2]], , drop = FALSE]

  # dimensions and indices
  nclusters <- lp$nclusters[[2]]
  cluster_idx <- lp$cluster.idx[[2]]
  between_idx <- lp$between.idx[[2]]
  both_idx <- lp$both.idx[[2]]
  nz <- length(between_idx)
  nb <- length(both_idx)
  zp <- mp$Zp
  if (nz > 0L) {
    y1w <- y1[, -between_idx, drop = FALSE]
    z <- y2[, between_idx, drop = FALSE]
    zc <- t(t(z) - mu_z)
  } else {
    y1w <- y1
    z <- NULL
    zc <- NULL
  }
  ny <- ncol(y1w)
  y1w_c <- t(t(y1w) - mu_y)

  # per-pattern quantities + A.j and p.j per cluster
  sigma_w_inv <- solve.default(sw)
  winv_p <- breg_p <- ccond_p <- vector("list", mp$npatterns)
  alist_1 <- rep(list(matrix(0, nb, nb)), nclusters)
  pij <- matrix(0, nrow(y1w), ny)
  for (p in seq_len(mp$npatterns)) {
    o_idx <- which(mp$pat[p, ])
    na_idx <- which(!mp$pat[p, ])
    case_idx <- mp$case.idx[[p]]
    if (length(na_idx) > 0L) {
      wp_inv <- lav_mat_sym_inverse_update(
        s_inv = sigma_w_inv, rm_idx = na_idx
      )
      breg_p[[p]] <- sw[na_idx, o_idx, drop = FALSE] %*% wp_inv
      ccond_p[[p]] <- sw[na_idx, na_idx, drop = FALSE] -
        breg_p[[p]] %*% sw[o_idx, na_idx, drop = FALSE]
    } else {
      wp_inv <- sigma_w_inv
    }
    winv_p[[p]] <- wp_inv
    pij[case_idx, o_idx] <- y1w_c[case_idx, o_idx, drop = FALSE] %*% wp_inv
    a_full <- matrix(0, ny, ny)
    a_full[o_idx, o_idx] <- wp_inv
    a_bb <- a_full[both_idx, both_idx, drop = FALSE]
    j1_idx <- mp$j1.idx[[p]]
    j_freq <- mp$j.freq[[p]]
    for (k in seq_along(j1_idx)) {
      alist_1[[j1_idx[k]]] <- alist_1[[j1_idx[k]]] + (a_bb * j_freq[k])
    }
  }
  pj <- rowsum.default(pij[, both_idx, drop = FALSE], cluster_idx,
    reorder = FALSE, na.rm = TRUE
  )

  # z-patterns: prior of (beta, z.mis) given z.obs
  zmis_q <- zpat2j <- NULL
  if (nz > 0L) {
    zpat2j <- integer(nclusters)
    for (q in seq_len(zp$npatterns)) {
      zpat2j[zp$case.idx[[q]]] <- q
    }
    if (length(zp$empty.idx) > 0L) {
      zpat2j[zp$empty.idx] <- zp$npatterns + 1L
    }
    nqz <- zp$npatterns + (length(zp$empty.idx) > 0L)
    c0_q <- m0coef_q <- zmis_q <- vector("list", nqz)
    for (q in seq_len(nqz)) {
      if (q > zp$npatterns) {
        zo_idx <- integer(0L)
        zm_idx <- seq_len(nz)
      } else {
        zo_idx <- which(zp$pat[q, ])
        zm_idx <- which(!zp$pat[q, ])
      }
      zmis_q[[q]] <- zm_idx
      k_q <- rbind(
        cbind(sb, syz[, zm_idx, drop = FALSE]),
        cbind(
          t(syz[, zm_idx, drop = FALSE]),
          sigma_zz[zm_idx, zm_idx, drop = FALSE]
        )
      )
      if (length(zo_idx) > 0L) {
        l_q <- rbind(
          syz[, zo_idx, drop = FALSE],
          sigma_zz[zm_idx, zo_idx, drop = FALSE]
        )
        soo_inv <- solve.default(sigma_zz[zo_idx, zo_idx, drop = FALSE])
        m0coef_q[[q]] <- l_q %*% soo_inv
        c0_q[[q]] <- k_q - m0coef_q[[q]] %*% t(l_q)
      } else {
        m0coef_q[[q]] <- matrix(0, nb + length(zm_idx), 0L)
        c0_q[[q]] <- k_q
      }
    }
  }

  # per-cluster posterior of w = (beta, z.mis): N(m1, C1)
  m1_list <- c1_list <- vector("list", nclusters)
  for (j in seq_len(nclusters)) {
    a_j <- alist_1[[j]]
    p_j <- pj[j, ]
    if (nz > 0L) {
      q <- zpat2j[j]
      c0 <- c0_q[[q]]
      zm_idx <- zmis_q[[q]]
      nzm <- length(zm_idx)
      if (q > zp$npatterns) {
        m0 <- numeric(nb + nzm)
      } else {
        zo_idx <- which(zp$pat[q, ])
        m0 <- drop(m0coef_q[[q]] %*% zc[j, zo_idx])
      }
    } else {
      c0 <- sb
      nzm <- 0L
      m0 <- numeric(nb)
    }
    nw <- nb + nzm
    m_j <- c0[, seq_len(nb), drop = FALSE] %*% a_j
    m_j <- cbind(m_j, matrix(0, nw, nzm))
    m_j[lav_mat_diag_idx(nw)] <- m_j[lav_mat_diag_idx(nw)] + 1
    rhs <- cbind(m0 + drop(c0[, seq_len(nb), drop = FALSE] %*% p_j), c0)
    sol <- solve.default(m_j, rhs)
    m1_list[[j]] <- sol[, 1L]
    c1 <- sol[, -1L, drop = FALSE]
    c1_list[[j]] <- (c1 + t(c1)) / 2
  }

  list(
    m1 = m1_list, c1 = c1_list,
    breg = breg_p, ccond = ccond_p, winv = winv_p,
    alist = alist_1, pj = pj,
    y1w = y1w, y1w_c = y1w_c, z = z, zc = zc,
    mu_y = mu_y, mu_z = mu_z,
    sb = sb, syz = syz, sigma_zz = sigma_zz, sigma_w = sw,
    mu_b_2l = out$mu.b, mu_w_2l = out$mu.w,
    both_idx = both_idx, between_idx = between_idx,
    zmis_q = zmis_q, zpat2j = zpat2j
  )
}

# the E-step 'random effects': the posterior means of the cluster-level
# components nu_j = (mu_b + beta_j), given ALL the observed data; this
# is the missing-data counterpart of lav_mvn_cl_em_estep_ranef()
#
# returns a nclusters x ny matrix (ny = number of level-1 variables) in
# the same convention as the complete-data version: for both-level
# variables, the (uncentered, total-mean scale) posterior mean of nu_j;
# zero for within-only variables
#
# if se = TRUE, the posterior standard deviations are returned in the
# "se" attribute
#
# if impute = TRUE, two extra attributes are returned:
# - "y1w.imputed": the level-1 data (level-1 variables only), with the
#   missing values replaced by their posterior means
#   E[y.mis | y.obs, E(beta)] (exact, as the conditional mean is linear
#   in beta)
# - "z.imputed": the cluster-level data (between-only variables, one
#   row per cluster), with the missing values replaced by their
#   posterior means
lav_mvn_cl_mi_estep_ranef <- function(y1 = NULL,
                                      y2 = NULL,
                                      lp = NULL,
                                      mp = NULL,
                                      mu_w = NULL, # implied Mu.W
                                      sigma_w = NULL, # implied Sigma.W
                                      mu_b = NULL, # implied Mu.B
                                      sigma_b = NULL, # implied Sigma.B
                                      se = FALSE,
                                      impute = FALSE) {
  post <- lav_mvn_cl_mi_posterior(
    y1 = y1, y2 = y2, lp = lp, mp = mp,
    mu_w = mu_w, sigma_w = sigma_w, mu_b = mu_b, sigma_b = sigma_b
  )
  nclusters <- lp$nclusters[[2]]
  both_idx <- post$both_idx
  nb <- length(both_idx)
  nz <- length(post$between_idx)
  ny <- ncol(post$y1w)

  # E(nu_j | all observed data), in the y1w space (cf. the complete-data
  # version: mu_b_2l[both] is the total mean of the both-level variables)
  eb <- matrix(0, nclusters, nb)
  for (j in seq_len(nclusters)) {
    eb[j, ] <- post$m1[[j]][seq_len(nb)]
  }
  mb_j <- matrix(0, nrow = nclusters, ncol = ny)
  mb_j[, both_idx] <- rep(post$mu_b_2l[both_idx], each = nclusters) + eb

  # posterior standard deviations
  if (se) {
    se_j <- matrix(0, nrow = nclusters, ncol = ny)
    for (j in seq_len(nclusters)) {
      c1_diag <- diag(post$c1[[j]])[seq_len(nb)]
      c1_diag[c1_diag < 0] <- 0
      se_j[j, both_idx] <- sqrt(c1_diag)
    }
    attr(mb_j, "se") <- se_j
  }

  # impute the missing values by their posterior means
  if (impute) {
    # level-1 variables: E[y.mis | y.obs, beta] is linear in beta, so
    # plugging in E[beta] gives the exact posterior mean
    y1w_i <- post$y1w
    eb_full <- matrix(0, nclusters, ny)
    eb_full[, both_idx] <- eb
    for (p in seq_len(mp$npatterns)) {
      na_idx <- which(!mp$pat[p, ])
      if (length(na_idx) == 0L) {
        next
      }
      o_idx <- which(mp$pat[p, ])
      case_idx <- mp$case.idx[[p]]
      j_idx <- mp$j.idx[[p]]
      f_o <- post$y1w_c[case_idx, o_idx, drop = FALSE] -
        eb_full[j_idx, o_idx, drop = FALSE]
      y1w_i[case_idx, na_idx] <-
        t(t(f_o %*% t(post$breg[[p]])) + post$mu_y[na_idx]) +
        eb_full[j_idx, na_idx, drop = FALSE]
    }
    # empty units (all level-1 variables missing)
    if (length(mp$empty.idx) > 0L) {
      cl_empty <- lp$cluster.idx[[2]][mp$empty.idx]
      y1w_i[mp$empty.idx, ] <-
        t(t(eb_full[cl_empty, , drop = FALSE]) + post$mu_y)
    }
    attr(mb_j, "y1w.imputed") <- y1w_i

    # between-only variables (one row per cluster)
    if (nz > 0L) {
      z_i <- post$z
      for (j in seq_len(nclusters)) {
        zm_idx <- post$zmis_q[[post$zpat2j[j]]]
        if (length(zm_idx) > 0L) {
          z_i[j, zm_idx] <- post$mu_z[zm_idx] +
            post$m1[[j]][nb + seq_along(zm_idx)]
        }
      }
      attr(mb_j, "z.imputed") <- z_i
    }
  }

  mb_j
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
