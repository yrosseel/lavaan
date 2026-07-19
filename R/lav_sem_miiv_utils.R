# varcov solution is:
# theta_2 = solve(t(Delta2) %*% W2 %*% Delta2) %*% t(Delta2) %*% W2 %*% s_vech

# we need to compute the Jacobian of theta_2 wrt the elements of s_vech
#
# four options:
# - ULS
# - GLS
# - 2RLS
# - RLS

# based on lav_utils_wls_linearization(), but only for ULS and H
lav_sem_miiv_utils_jacb_uls <- function(sample_cov = NULL, delta2 = NULL) {
  nvar <- nrow(sample_cov)
  pstar <- nvar * (nvar + 1L) / 2L

  jac_cov <- delta2
  # nc <- ncol(jac_cov)
  meanstructure <- FALSE
  if (nrow(delta2) > pstar) {
    meanstructure <- TRUE
  }

  # split delta2/svec
  if (meanstructure) {
    mean_idx <- seq_len(nvar)
    jac_mean <- delta2[mean_idx, , drop = FALSE]
    jac_cov <- delta2[-mean_idx, , drop = FALSE]
  }

  # vech
  w <- rep(1.0, pstar)
  w[lav_mat_diagh_idx(nvar)] <- 0.5

  q_cov <- w * jac_cov
  a <- crossprod(jac_cov, q_cov)
  if (meanstructure) {
    a <- a + crossprod(jac_mean)
  }

  # tQ only needed if H is requested
  t_q <- t(q_cov)
  if (meanstructure) t_q <- cbind(t(jac_mean), t_q)

  r <- chol(a)
  ainv  <- chol2inv(r)
  h <- ainv %*% t_q

  h
}

# this is the GLS version where W2 = t(D) %*% (S.inv %x% S.inv) %*% D
lav_sem_miiv_utils_jacb_gls <- function(sample_cov = NULL, delta2 = NULL,
                                        theta2 = NULL) {
  nvar  <- nrow(sample_cov)
  pstar <- nvar * (nvar + 1L) / 2L
  s_vech <- lav_mat_vech(sample_cov)

  # vech weights
  w <- rep(1.0, pstar)
  w[lav_mat_diagh_idx(nvar)] <- 0.5

  # only covariance block
  meanstructure_flag <- FALSE
  if (nrow(delta2) > pstar) {
    meanstructure_flag <- TRUE
    delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
  }

  nc <- ncol(delta2)
  c_s <- chol(sample_cov)
  si <- chol2inv(c_s)

  # q_cov[,j] = w * vech(Si %*% Vj %*% Si),  Vj = vech_inv(delta2[,j])
  # equivalent to (1/2) W2 %*% delta2
  q_cov <- matrix(0.0, pstar, nc)
  for (j in seq_len(nc)) {
    v_j <- lav_mat_vech_rev(delta2[, j])
    si_vj <- si %*% v_j
    q_cov[, j] <- w * lav_mat_vech(si_vj %*% si)
  }

  # M_half = (1/2) M = t(delta2) %*% q_cov
  m_half <- crossprod(delta2, q_cov)
  b_half <- drop(crossprod(q_cov, s_vech))

  # theta2: GLS estimate
  # (scaling cancels: M_half^{-1} b_half = M^{-1} t(delta2) W2 s_vech)
  r <- chol(m_half)
  if (is.null(theta2)) {
    theta2 <- backsolve(r, forwardsolve(t(r), b_half))
  }

  # error term and GLS correction q_mat = Si %*% e_mat %*% Si
  e_mat <- lav_mat_vech_rev(s_vech - delta2 %*% theta2)
  q_mat <- si %*% e_mat %*% si

  # correction_cov[,j] = w * vech(Tj + t(Tj)),  Tj = Si %*% Vj %*% q_mat
  # equivalent to (1/2) correction %*% delta2
  correction_cov <- matrix(0.0, pstar, nc)
  for (j in seq_len(nc)) {
    v_j              <- lav_mat_vech_rev(delta2[, j])
    t_j              <- si %*% v_j %*% q_mat
    correction_cov[, j] <- w * lav_mat_vech(t_j + t(t_j))
  }

  # out = M^{-1} %*% t(delta2) %*% (W2 - correction)
  #     = M_half^{-1} %*% t(q_cov - correction_cov)
  out <- backsolve(r, forwardsolve(t(r), t(q_cov - correction_cov)))

  # meanstructure: prepend zero block for mean parameters
  if (meanstructure_flag) {
    out <- cbind(matrix(0.0, nrow = nrow(out), ncol = nvar), out)
  }

  out
}


# 2RLS solution, where W2 = t(D) %*% (Sigma.inv %x% Sigma.inv) %*% D
# and Sigma is based on the ULS estimate of theta2
lav_sem_miiv_utils_jacb_2rls <- function(sample_cov = NULL, delta2 = NULL,
                                         theta2 = NULL) {

  nvar  <- nrow(sample_cov)
  pstar <- nvar * (nvar + 1L) / 2L
  s_vech <- lav_mat_vech(sample_cov)

  # vech weights: 1 for off-diagonal positions, 0.5 for diagonal positions
  w <- rep(1.0, pstar)
  w[lav_mat_diagh_idx(nvar)] <- 0.5

  # only covariance block
  meanstructure_flag <- FALSE
  if (nrow(delta2) > pstar) {
    meanstructure_flag <- TRUE
    delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
  }
  nc <- ncol(delta2)

  # step 1: ULS estimate of theta_uls
  q_uls <- w * delta2
  m_uls <- crossprod(delta2, q_uls)
  r_uls <- chol(m_uls)
  theta_uls <-
    backsolve(r_uls, forwardsolve(t(r_uls), drop(crossprod(q_uls, s_vech))))

  # sigma to be used for GLS estimate
  sigma <- lav_mat_vech_rev(delta2 %*% theta_uls)
  c_si <- chol(sigma)
  sigma_inv <- chol2inv(c_si)

  # step 2: GLS estimate of theta using sigma_inv
  # q_cov[,j] = w * vech(Si Vj Si)  =  (1/2) w2 delta2[,j]
  q_cov <- matrix(0.0, pstar, nc)
  for (j in seq_len(nc)) {
    v_j <- lav_mat_vech_rev(delta2[, j])
    q_cov[, j] <- w * lav_mat_vech(sigma_inv %*% v_j %*% sigma_inv)
  }
  m_sigma_half <- crossprod(delta2, q_cov)
  r_sigma <- chol(m_sigma_half)
  if (is.null(theta2)) {
    theta2 <- backsolve(r_sigma,
      forwardsolve(t(r_sigma), drop(crossprod(q_cov, s_vech))))
  }

  # step 3: correction matrix (as column operations, no c_mat formed)
  # c_cov[,j] = w * vech(Tj + t(Tj)),  Tj = Si Vj q_mat
  #           = (1/2) c_mat delta2[,j]
  e_mat <- lav_mat_vech_rev(s_vech - delta2 %*% theta2)
  q_mat <- sigma_inv %*% e_mat %*% sigma_inv

  c_cov <- matrix(0.0, pstar, nc)
  for (j in seq_len(nc)) {
    v_j <- lav_mat_vech_rev(delta2[, j])
    t_j <- sigma_inv %*% v_j %*% q_mat
    c_cov[, j] <- w * lav_mat_vech(t_j + t(t_j))
  }

  # step 4: Jacobian
  #   out = M_sigma_half^{-1} (t(q_cov) - t(c_cov) delta2 m_uls^{-1} t(q_uls))
  # build rhs = t(q_cov) - t(c_cov) delta2 m_uls^{-1} t(q_uls)  (nc x pstar)
  # Step A: CtD = t(c_cov) delta2  (nc x nc)
  ct_d  <- crossprod(c_cov, delta2)
  mi_ct_d <- backsolve(r_uls, forwardsolve(t(r_uls), t(ct_d)))
  # Step B: rhs = t(q_cov) - MiCtD^T t(q_uls)  (nc x pstar)
  rhs <- t(q_cov) - t(mi_ct_d) %*% t(q_uls)
  # Step C: out = M_sigma_half^{-1} rhs
  out <- backsolve(r_sigma, forwardsolve(t(r_sigma), rhs))

  # meanstructure: prepend zero block for mean parameters
  if (meanstructure_flag) {
    out <- cbind(matrix(0.0, nrow = nrow(out), ncol = nvar), out)
  }

  out
}

# RLS solution
lav_sem_miiv_utils_jacb_rls <- function(sample_cov = NULL, delta2 = NULL,
                                        theta2 = NULL) {
  nvar  <- nrow(sample_cov)
  pstar <- nvar * (nvar + 1L) / 2L
  s_vech <- lav_mat_vech(sample_cov)

  # vech weights: 0.5 for diagonal positions, 1 for off-diagonal positions
  w <- rep(1.0, pstar)
  w[lav_mat_diagh_idx(nvar)] <- 0.5

  # only covariance block
  meanstructure_flag <- FALSE
  if (nrow(delta2) > pstar) {
    meanstructure_flag <- TRUE
    delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
  }

  nc <- ncol(delta2)

  # Pre-compute vech_reverse of each column of delta2 once; used in the
  # final block (and in the loop when theta2 is not supplied).
  v_list <- vector("list", nc)
  for (j in seq_len(nc)) v_list[[j]] <- lav_mat_vech_rev(delta2[, j])

  q_cov <- matrix(0.0, pstar, nc)

  # using final theta2 estimate
  new_sigma <- lav_mat_vech_rev(delta2 %*% theta2)
  c_si       <- chol(new_sigma)
  sigma_inv <- chol2inv(c_si)
  for (j in seq_len(nc)) {
    q_cov[, j] <- w * lav_mat_vech(sigma_inv %*% v_list[[j]] %*% sigma_inv)
  }
  m_mat <- crossprod(delta2, q_cov)
  r_mat <- chol(m_mat)

  # correction columns: c_cov[,j] = w * vech(Tj + t(Tj)),  Tj = Si Vj q_mat
  # where q_mat = Si e_mat Si.
  # JW_theta2 = -c_cov  (dW_e2_from_v(delta2[,j]) = -c_cov[,j], derived above)
  e_mat <- lav_mat_vech_rev(s_vech - delta2 %*% theta2)
  q_mat <- sigma_inv %*% e_mat %*% sigma_inv
  c_cov <- matrix(0.0, pstar, nc)
  for (j in seq_len(nc)) {
    t_j <- sigma_inv %*% v_list[[j]] %*% q_mat
    c_cov[, j] <- w * lav_mat_vech(t_j + t(t_j))
  }

  # Jacobian:
  #   A           = m_mat_inv t(delta2) JW_theta2
  #               = -solve(m_mat, crossprod(delta2, c_cov))
  #   IminusA_inv = solve(I - A)
  #   out         = IminusA_inv %*% solve(m_mat, t(q_cov))
  dt_c <- crossprod(delta2, c_cov)
  a <- -backsolve(r_mat, forwardsolve(t(r_mat), dt_c))
  iminus_a_inv <- solve(diag(nc) - a)
  minv_t_q <- backsolve(r_mat, forwardsolve(t(r_mat), t(q_cov)))
  out     <- iminus_a_inv %*% minv_t_q

  # meanstructure: prepend zero block for mean parameters
  if (meanstructure_flag) {
    out <- cbind(matrix(0.0, nrow = nrow(out), ncol = nvar), out)
  }

  out
}


# we also need to compute the Jacobian of theta_2 wrt the elements theta_1
#
# four options:
# - ULS
# - GLS
# - 2RLS
# - RLS

# ULS and GLS: only delta2 depends on theta1
lav_sem_miiv_utils_jaca_uls_gls <- function(lavmodel = NULL,  # nolint
                                            lavpartable = NULL,
                                            lavh1 = NULL,
                                            free_directed_idx = integer(0L),
                                            free_undirected_idx = integer(0L),
                                            iv_varcov_method = "ULS") {
  nblocks <- lavmodel@nblocks

  # delta across all blocks
  delta_block <- lav_model_delta(lavmodel = lavmodel)

  # augment lavpartable to include model matrices and row/col indices
  mm_info <- lav_lisrel(lavpartable, target = NULL, extra = FALSE)

  jac_a <- matrix(0, length(free_undirected_idx), length(free_directed_idx))
  for (b in seq_len(nblocks)) {
    # s_vech
    sample_cov <- lavh1$implied$cov[[b]]
    nvar <- nrow(sample_cov)
    s_vech <- lav_mat_vech(sample_cov)

    # delta2
    delta2 <- delta_block[[b]][, free_undirected_idx, drop = FALSE]
    if (lavmodel@meanstructure) {
      delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
    }

    # w2
    if (iv_varcov_method == "GLS") {
      s_inv <- solve(sample_cov)
    } else {
      s_inv <- diag(1, nrow = nrow(sample_cov))
    }
    w2_22 <- 0.5 * lav_mat_dup_pre_post(s_inv %x% s_inv)
    w2 <- w2_22

    # MLIST
    mm_in_group <- seq_len(lavmodel@nmat[b]) + cumsum(c(0, lavmodel@nmat))[b]
    mlist <- lavmodel@GLIST[mm_in_group]
    if (!is.null(mlist$beta)) {
      mlist$IB.inv <- lav_lisrel_ibinv(mlist)
    }

    m_mat <- t(delta2) %*% w2 %*% delta2
    theta2 <- solve(m_mat, t(delta2) %*% w2 %*% s_vech)
    m_inv <- solve(m_mat)
    w2e <- w2 %*% (s_vech - delta2 %*% theta2)
    delta2tw2 <- t(delta2) %*% w2

    # container for Jacobian for this block
    for (k in seq_along(free_directed_idx)) {
      d_idx <- match(free_directed_idx[k], lavpartable$free)[1]
      # correct block?
      if (lavpartable$block[d_idx] != b) {
        next
      }
      d_mat <- mm_info$mat[d_idx]
      d_row <- mm_info$row[d_idx]
      d_col <- mm_info$col[d_idx]
      d_deltak <- matrix(0, nrow(delta2), ncol(delta2))
      # we only need to fill the columns of dDeltak that correspond to the
      # psi elements in free.undirected.idx
      for (i in seq_along(free_undirected_idx)) {
        un_idx <- match(free_undirected_idx[i], lavpartable$free)[1]
        un_mat <- mm_info$mat[un_idx]
        if (un_mat == "theta") {
          next
        }
        un_row <- mm_info$row[un_idx]
        un_col <- mm_info$col[un_idx]

        if (d_mat == "lambda") {
          tmp <- lav_lisrel_d2sigma_lambda_psi(
            mlist = mlist,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        } else if (d_mat == "beta") {
          tmp <- lav_lisrel_d2sigma_beta_psi(
            mlist = mlist,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        }
        d_deltak[, i] <- lav_mat_vech(tmp)
      }
      jac_a[, k] <-
        m_inv %*% (t(d_deltak) %*% w2e - delta2tw2 %*% d_deltak %*% theta2)
    }
  }
  jac_a
}

# 2RLS
# - first step: only Delta_2 depends on theta1
# - second step: W_2 also depends on theta1!
lav_sem_miiv_utils_jaca_2rls <- function(lavmodel = NULL,
                                         lavpartable = NULL,
                                         lavh1 = NULL,
                                         free_directed_idx = integer(0L),
                                         free_undirected_idx = integer(0L)) {
  nblocks <- lavmodel@nblocks

  # delta across all blocks
  delta_block <- lav_model_delta(lavmodel = lavmodel)

  # augment lavpartable to include model matrices and row/col indices
  mm_info <- lav_lisrel(lavpartable, target = NULL, extra = FALSE)

  # first get jac_a for ULS (step 1)
  jac_uls <- lav_sem_miiv_utils_jaca_uls_gls(
    lavmodel = lavmodel,
    lavpartable = lavpartable, lavh1 = lavh1,
    free_directed_idx = free_directed_idx,
    free_undirected_idx = free_undirected_idx, iv_varcov_method = "ULS"
  )

  jac_a <- matrix(0, length(free_undirected_idx), length(free_directed_idx))
  for (b in seq_len(nblocks)) {
    # s_vech
    sample_cov <- lavh1$implied$cov[[b]]
    nvar <- nrow(sample_cov)
    s_vech <- lav_mat_vech(sample_cov)

    # delta2
    delta2 <- delta_block[[b]][, free_undirected_idx, drop = FALSE]
    if (lavmodel@meanstructure) {
      delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
    }
    # w2
    s_inv <- diag(1, nrow = nrow(sample_cov))
    w2_uls <- 0.5 * lav_mat_dup_pre_post(s_inv %x% s_inv)

    # MLIST
    mm_in_group <- seq_len(lavmodel@nmat[b]) + cumsum(c(0, lavmodel@nmat))[b]
    mlist <- lavmodel@GLIST[mm_in_group]
    if (!is.null(mlist$beta)) {
      mlist$IB.inv <- lav_lisrel_ibinv(mlist)
    }

    # step 1
    m_uls <- t(delta2) %*% w2_uls %*% delta2
    m_uls_inv <- solve(m_uls)
    theta_uls <- drop(m_uls_inv %*% t(delta2) %*% w2_uls %*% s_vech)

    # step 2
    new_sigma <- lav_mat_vech_rev(delta2 %*% theta_uls)
    sigma_inv <- solve(new_sigma)
    w2 <- 0.5 * lav_mat_dup_pre_post(sigma_inv %x% sigma_inv)
    m_mat <- t(delta2) %*% w2 %*% delta2
    m_mat_inv <- solve(m_mat)
    theta2 <- drop(m_mat_inv %*% t(delta2) %*% w2 %*% s_vech)
    e2 <- s_vech - delta2 %*% theta2
    e_mat <- lav_mat_vech_rev(e2)

    # pre-compute
    minv_delta2t <- m_mat_inv %*% t(delta2)
    w2_e2 <- w2 %*% e2
    delta2t_w2 <- t(delta2) %*% w2
    sigma_inv_e <- sigma_inv %*% e_mat

    for (k in seq_along(free_directed_idx)) {
      d_idx <- match(free_directed_idx[k], lavpartable$free)[1]
      # correct block?
      if (lavpartable$block[d_idx] != b) {
        next
      }
      d_mat <- mm_info$mat[d_idx]
      d_row <- mm_info$row[d_idx]
      d_col <- mm_info$col[d_idx]
      d_deltak <- matrix(0, nrow(delta2), ncol(delta2))
      # we only need to fill the columns of dDeltak that correspond to the
      # psi elements in free.undirected.idx
      for (i in seq_along(free_undirected_idx)) {
        un_idx <- match(free_undirected_idx[i], lavpartable$free)[1]
        un_mat <- mm_info$mat[un_idx]
        if (un_mat == "theta") {
          next
        }
        un_row <- mm_info$row[un_idx]
        un_col <- mm_info$col[un_idx]

        if (d_mat == "lambda") {
          tmp <- lav_lisrel_d2sigma_lambda_psi(
            mlist = mlist,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        } else if (d_mat == "beta") {
          tmp <- lav_lisrel_d2sigma_beta_psi(
            mlist = mlist,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        }
        d_deltak[, i] <- lav_mat_vech(tmp)
      }

      # d vech(Sigma)/d theta1[k]: product rule on unvech(Delta2 theta_uls)
      v_k <- d_deltak %*% theta_uls + delta2 %*% jac_uls[, k]
      d_sigma_k <- lav_mat_vech_rev(v_k)

      # d_sigma_inv = -sigma.inv dSigma_k sigma.inv
      d_sigma_inv <- -sigma_inv %*% d_sigma_k %*% sigma_inv

      # the tricky one:
      # (dw2_k) e2 = (1/2)
      # D^T vec(sigma_inv e_mat dsigma_inv_k + dsigma_invk e_mat sigma_inv)
      sym_part <- sigma_inv_e %*% d_sigma_inv + d_sigma_inv %*% t(sigma_inv_e)
      dw2_e2 <-
        0.5 * drop(lav_mat_dup_pre(as.matrix(as.vector(sym_part))))

      # term 1: direct delta2 variation in final WLS
      term1 <-
        m_mat_inv %*% (t(d_deltak) %*% w2_e2 -
                                           delta2t_w2 %*% d_deltak %*% theta2)

      # term 2: W2 variation via Sigma = unvech(delta2 theta_uls)
      term2 <- minv_delta2t %*% dw2_e2

      jac_a[, k] <- term1 + term2
    }
  }
  jac_a
}

# RLS
# - first step: only Delta_2 depends on theta1
# - second step: iteratively updating Sigma -> W_2 also depends on theta1!
lav_sem_miiv_utils_jaca_rls <- function(lavmodel = NULL,
                                        lavpartable = NULL,
                                        lavh1 = NULL,
                                        free_directed_idx = integer(0L),
                                        free_undirected_idx = integer(0L)) {
  nblocks <- lavmodel@nblocks

  # delta across all blocks
  delta_block <- lav_model_delta(lavmodel = lavmodel)

  # augment lavpartable to include model matrices and row/col indices
  mm_info <- lav_lisrel(lavpartable, target = NULL, extra = FALSE)

  jac_a <- matrix(0, length(free_undirected_idx), length(free_directed_idx))
  for (b in seq_len(nblocks)) {
    # s_vech
    sample_cov <- lavh1$implied$cov[[b]]
    nvar <- nrow(sample_cov)
    s_vech <- lav_mat_vech(sample_cov)

    # delta2
    delta2 <- delta_block[[b]][, free_undirected_idx, drop = FALSE]
    if (lavmodel@meanstructure) {
      delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
    }

    # w2
    s_inv <- diag(1, nrow = nrow(sample_cov))
    w2_uls <- 0.5 * lav_mat_dup_pre_post(s_inv %x% s_inv)

    # MLIST
    mm_in_group <- seq_len(lavmodel@nmat[b]) + cumsum(c(0, lavmodel@nmat))[b]
    mlist <- lavmodel@GLIST[mm_in_group]
    if (!is.null(mlist$beta)) {
      mlist$IB.inv <- lav_lisrel_ibinv(mlist)
    }

    # step 1
    m_uls <- t(delta2) %*% w2_uls %*% delta2
    m_uls_inv <- solve(m_uls)
    theta_uls <- drop(m_uls_inv %*% t(delta2) %*% w2_uls %*% s_vech)

    # step 2
    theta2 <- theta_uls
    for (i in seq_len(200)) {
      old_x <- theta2
      new_sigma <- lav_mat_vech_rev(delta2 %*% theta2)
      sigma_inv <- solve(new_sigma)
      w2 <- 0.5 * lav_mat_dup_pre_post(sigma_inv %x% sigma_inv)
      m_mat <- t(delta2) %*% w2 %*% delta2
      theta2 <- drop(solve(m_mat, t(delta2) %*% w2 %*% s_vech))
      if (sum((old_x - theta2)^2) < 1e-12 * (1 + sum(theta2^2))) break
    }

    # final quantities
    new_sigma <- lav_mat_vech_rev(delta2 %*% theta2)
    sigma_inv <- solve(new_sigma)
    w2 <- 0.5 * lav_mat_dup_pre_post(sigma_inv %x% sigma_inv)
    m_mat <- t(delta2) %*% w2 %*% delta2
    m_mat_inv <- solve(m_mat)
    e2 <- s_vech - delta2 %*% theta2
    e_mat <- lav_mat_vech_rev(e2)
    sigma_inv_e <- sigma_inv %*% e_mat

    minv_delta2t <- m_mat_inv %*% t(delta2)
    w2_e2 <- w2 %*% e2
    delta2t_w2 <- t(delta2) %*% w2

    # helper: compute (dW)e2 given d vech(Sigma) as a vector v
    d_w_e2_from_v <- function(v) {
      d_sigma <- lav_mat_vech_rev(v)
      d_sigma_inv <- -sigma_inv %*% d_sigma %*% sigma_inv
      sym_part <- sigma_inv_e %*% d_sigma_inv + d_sigma_inv %*% t(sigma_inv_e)
      dw2_e2 <-
        0.5 * drop(lav_mat_dup_pre(as.matrix(as.vector(sym_part))))
      dw2_e2
    }

    # A = M2^{-1} delta2^T J_W^(theta2): q2 x q2
    # column j of J_W^(theta2): dSigma = unvech(delta2[:,j])
    jw_theta2 <- matrix(0.0, nrow(delta2), ncol(delta2))
    for (j in seq_len(ncol(delta2))) {
      jw_theta2[, j] <- d_w_e2_from_v(delta2[, j])
    }
    a <- minv_delta2t %*% jw_theta2
    iminus_a_inv <- solve(diag(ncol(delta2)) - a)

    for (k in seq_along(free_directed_idx)) {
      d_idx <- match(free_directed_idx[k], lavpartable$free)[1]
      # correct block?
      if (lavpartable$block[d_idx] != b) {
        next
      }
      d_mat <- mm_info$mat[d_idx]
      d_row <- mm_info$row[d_idx]
      d_col <- mm_info$col[d_idx]
      d_deltak <- matrix(0, nrow(delta2), ncol(delta2))
      # we only need to fill the columns of dDeltak that correspond to the
      # psi elements in free.undirected.idx
      for (i in seq_along(free_undirected_idx)) {
        un_idx <- match(free_undirected_idx[i], lavpartable$free)[1]
        un_mat <- mm_info$mat[un_idx]
        if (un_mat == "theta") {
          next
        }
        un_row <- mm_info$row[un_idx]
        un_col <- mm_info$col[un_idx]

        if (d_mat == "lambda") {
          tmp <- lav_lisrel_d2sigma_lambda_psi(
            mlist = mlist,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        } else if (d_mat == "beta") {
          tmp <- lav_lisrel_d2sigma_beta_psi(
            mlist = mlist,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        }
        d_deltak[, i] <- lav_mat_vech(tmp)
      }

      # term 1: direct delta2 variation in final WLS
      term1 <-
        m_mat_inv %*% (t(d_deltak) %*% w2_e2 -
                                            delta2t_w2 %*% d_deltak %*% theta2)

      # term 2: W variation from dSigma^(Delta2) = unvech(dDelta2_k %*% theta2)
      # (the dSigma^(theta2) part is handled implicitly by (I-A)^{-1})
      v_k <- as.vector(d_deltak %*% theta2)
      term2 <- minv_delta2t %*% d_w_e2_from_v(v_k)

      jac_a[, k] <- iminus_a_inv %*% (term1 + term2)
    }
  }
  jac_a
}


# jac_k: d theta1 / d svec
# see Fisher & Bollen (2020) section 2.4

# d(bvec - amat %*% beta) / d vech(S) for one equation of the moments
# engine, at an ARBITRARY slope vector beta (the expression is linear in
# beta). For an OLS equation (no instruments) bvec - amat beta =
# s_xy - s_xx beta; for a 2SLS equation it is S_xz S_zz^{-1} (s_zy -
# s_zx beta). Rows = ALL slopes of the equation (free and fixed);
# columns = vech(S) of the block, diagonal included.
#
# Evaluated at the equation's own unconstrained solution this is
# amat %*% k_mat (the per-equation moment Jacobian of Fisher & Bollen,
# 2020, section 2.4); evaluated at the POOLED solution it provides the
# right-hand side of the pooled (equality-constrained) Jacobian chain
# (see lav_sem_miiv_jack_analytic).
lav_sem_miiv_eq_dresid_dvech <- function(sample_cov = NULL,
                                         x_idx = NULL, y_idx = NULL,
                                         i_idx = integer(0L),
                                         beta = NULL) {
  nvar <- nrow(sample_cov)
  pstar <- nvar * (nvar + 1L) / 2L
  nx <- length(x_idx)
  if (nx == 0L) {
    return(matrix(0, nrow = 0L, ncol = pstar))
  }
  iv_flag <- length(i_idx) > 0L

  # Ex[k, j] = 1 iff x_idx[k] == j
  ex <- matrix(0, nrow = nx, ncol = nvar)
  ex[cbind(seq_len(nx), x_idx)] <- 1.0
  a_vec <- lav_mat_vech_row_idx(nvar)
  b_vec <- lav_mat_vech_col_idx(nvar)
  offdiag <- (a_vec > b_vec)

  if (!iv_flag) {
    # d(s_xy - s_xx beta) / d vech(S)
    g <- -as.vector(crossprod(ex, beta))
    g[y_idx] <- g[y_idx] + 1.0
    m <- sweep(ex[, b_vec, drop = FALSE], 2L, g[a_vec], `*`) +
      sweep(ex[, a_vec, drop = FALSE], 2L, g[b_vec] * offdiag, `*`)
  } else {
    # d(S_xz S_zz^{-1} (s_zy - s_zx beta)) / d vech(S)
    nz <- length(i_idx)
    s_zx <- sample_cov[i_idx, x_idx, drop = FALSE]
    s_zz <- sample_cov[i_idx, i_idx, drop = FALSE]
    s_zy <- sample_cov[i_idx, y_idx, drop = FALSE]
    m_w <- lav_mat_sym_solve_spd(s_zz, s_zx)
    r_z <- drop(s_zy) - s_zx %*% beta
    u <- drop(lav_mat_sym_solve_spd(s_zz, r_z))
    ez <- matrix(0, nrow = nz, ncol = nvar)
    ez[cbind(seq_len(nz), i_idx)] <- 1.0
    ew <- crossprod(m_w, ez)
    cu <- as.vector(crossprod(ez, u))
    cb_x <- as.vector(crossprod(ex, beta))
    h <- -cu
    h[y_idx] <- h[y_idx] + 1.0
    m <-
      sweep(ex[, b_vec, drop = FALSE], 2L, cu[a_vec], `*`) +
      sweep(ew[, a_vec, drop = FALSE], 2L, h[b_vec], `*`) -
      sweep(ew[, b_vec, drop = FALSE], 2L, cb_x[a_vec], `*`) +
      sweep(ex[, a_vec, drop = FALSE], 2L, cu[b_vec] * offdiag, `*`) +
      sweep(ew[, b_vec, drop = FALSE], 2L, h[a_vec] * offdiag, `*`) -
      sweep(ew[, a_vec, drop = FALSE], 2L, cb_x[b_vec] * offdiag, `*`)
  }
  m
}


# per equation, one block only!
lav_sem_miiv_utils_jack_eqs <- function(eqs = NULL, # one block only
                                        block = 1L,
                                        lavmodel = NULL,
                                        lavpartable = NULL,
                                        free_directed_idx = integer(0L)) {
  # continuous data only (for now)
  nvar <- lavmodel@nvar[block]
  if (lavmodel@categorical) {
    pstar <- nvar * (nvar - 1L) / 2
  } else {
    pstar <- nvar * (nvar + 1L) / 2
  }
  ntheta1 <- length(free_directed_idx)

  b <- block

  # K matrix
  if (lavmodel@categorical) {
    nc <- ncol(eqs[[1]][[1]]$k_mat)
    k_mat <- matrix(0.0, nrow = ntheta1, ncol = nc)
  } else if (lavmodel@meanstructure) {
    k_mat <- matrix(0.0, nrow = ntheta1, ncol = (nvar + pstar))
  } else {
    k_mat <- matrix(0.0, nrow = ntheta1, ncol = pstar)
  }

  # collect k_mat matrices from each equation
  for (j in seq_along(eqs[[b]])) {
    # this equation
    eq <- eqs[[b]][[j]]
    tmp <- lavpartable$free[eq$pt] # free index
    tmp <- tmp[!tmp == 0]
    free_idx <- match(tmp, free_directed_idx)
    nx <- nrow(eq$k_mat)
    if (length(free_idx) > 0L && nx > 0L) {
      if (all(free_idx > 0L)) {
        k_mat[free_idx, ] <- eq$k_mat
      } else {
        # remove non-free elements
        zero_idx <- which(free_idx == 0L)
        free_idx <- free_idx[-zero_idx]
        if (length(free_idx) > 0) {
          k_mat[free_idx, ] <- eq$k_mat[-zero_idx, , drop = FALSE]
        }
      }
    }
    if (lavmodel@meanstructure) {
      tmp <- lavpartable$free[eq$ptint]
      tmp <- tmp[!tmp == 0]
      free_int_idx <- match(tmp, free_directed_idx)
      if (length(free_int_idx) > 0L && free_int_idx > 0L) {
        k_mat[free_int_idx, ] <- eq$k_mat_int
      }
    }
  } # eq

  k_mat
}


# analytic directed Jacobian d theta1 / d svec over ALL blocks (the
# stacked moment vector), including the pooled solve for equality
# constraints among the directed coefficients.
#
# - without pooling (no shared slope columns, no shared intercepts, no
#   general linear constraints): block-diagonal assembly of the
#   per-equation Jacobians (eq$k_mat / eq$k_mat_int, Fisher & Bollen
#   2020) -- this also covers categorical data
# - with pooling (continuous moments engine only): the pooled solve is
#   theta = D^{-1} sum_e C_e' b_e with D = sum_e C_e' A_e C_e, so
#     d theta = D^{-1} sum_e C_e' d(b_e - A_e C_e theta)/d svec
#   where the per-equation term is the Fisher-Bollen M matrix evaluated
#   at the POOLED coefficients (lav_sem_miiv_eq_dresid_dvech is linear
#   in beta); general linear constraints add the usual KKT bordering,
#   and a label repeated within one equation collapses C'AC-style
#   (matching lav_sem_miiv_pool_directed). Pooled intercepts are the
#   nobs-weighted averages of ybar - xbar' beta(pooled), and their rows
#   chain through d theta.
#
# Returns ntheta1 x sum(md), or NULL when the configuration is not
# covered (pooling + categorical data); external instruments and the
# raw-data engine (no eq$k_mat) must be excluded by the caller.
lav_sem_miiv_jack_analytic <- function(eqs = NULL, lavmodel = NULL,
                                       lavpartable = NULL, lavdata = NULL,
                                       lavh1 = NULL,
                                       free_directed_idx = integer(0L),
                                       con = NULL, x = NULL,
                                       md = NULL, moff = NULL) {
  nblocks <- lavmodel@nblocks
  ntheta1 <- length(free_directed_idx)
  ntot <- sum(md)
  meanstructure <- lavmodel@meanstructure

  # pooling active? (mirror lav_sem_miiv_apply_directed_pool)
  all_gcol <- integer(0L)
  all_int <- integer(0L)
  for (b in seq_len(nblocks)) {
    for (eq in eqs[[b]]) {
      if (is.null(eq$slope_block)) {
        next
      }
      all_gcol <- c(all_gcol, eq$slope_block$gcol)
      if (!is.null(eq$ptint) && length(eq$ptint) == 1L) {
        fint <- lavpartable$free[eq$ptint]
        if (fint > 0L) {
          all_int <- c(all_int, fint)
        }
      }
    }
  }
  pool_active <- anyDuplicated(all_gcol) > 0L || !is.null(con) ||
    anyDuplicated(all_int) > 0L

  if (!pool_active) {
    # block-diagonal assembly of the per-equation Jacobians
    jac <- matrix(0, nrow = ntheta1, ncol = ntot)
    for (b in seq_len(nblocks)) {
      jac[, moff[b] + seq_len(md[b])] <- lav_sem_miiv_utils_jack_eqs(
        eqs = eqs, block = b, lavmodel = lavmodel,
        lavpartable = lavpartable, free_directed_idx = free_directed_idx
      )
    }
    return(jac)
  }

  # pooled chain: continuous moments engine only
  if (lavmodel@categorical) {
    return(NULL)
  }

  # 1. the pooled system D (same accumulation as
  #    lav_sem_miiv_pool_directed, including the within-equation
  #    duplicate-label collapse)
  free_slope_idx <- sort(unique(all_gcol))
  npool <- length(free_slope_idx)
  dmat <- matrix(0, nrow = npool, ncol = npool)
  for (b in seq_len(nblocks)) {
    for (eq in eqs[[b]]) {
      sb <- eq$slope_block
      if (is.null(sb)) {
        next
      }
      p <- match(sb$gcol, free_slope_idx)
      if (anyDuplicated(p) > 0L) {
        a_c <- rowsum(sb$amat, group = p)
        a_c <- t(rowsum(t(a_c), group = p))
        pu <- sort(unique(p))
        dmat[pu, pu] <- dmat[pu, pu] + a_c
      } else {
        dmat[p, p] <- dmat[p, p] + sb$amat
      }
    }
  }

  # 2. right-hand side: sum over equations of C' M(beta_pooled)
  rhs <- matrix(0, nrow = npool, ncol = ntot)
  beta_full_list <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    ov_names_b <- lavdata@ov.names[[b]]
    cov_b <- lavh1$implied$cov[[b]]
    nvar_b <- lavmodel@nvar[b]
    pstar_b <- nvar_b * (nvar_b + 1L) / 2L
    mean_off <- if (meanstructure) nvar_b else 0L
    cov_cols <- moff[b] + mean_off + seq_len(pstar_b)
    beta_full_list[[b]] <- vector("list", length(eqs[[b]]))
    for (j in seq_along(eqs[[b]])) {
      eq <- eqs[[b]][[j]]
      sb <- eq$slope_block
      if (is.null(sb)) {
        next
      }
      # instruments? (mirror the estimation loop)
      iv_flag <- TRUE
      if (!is.null(eq$iv_type)) {
        iv_flag <- eq$iv_type != "ols"
      } else if (identical(eq$rhs_new, eq$miiv)) {
        iv_flag <- FALSE
      }
      y_idx <- match(eq$lhs_new, ov_names_b)
      i_idx <- if (iv_flag) match(eq$iv, ov_names_b) else integer(0L)
      if (anyNA(y_idx) || anyNA(i_idx)) {
        # instruments outside the model moments (external instruments):
        # not covered by the pooled chain
        return(NULL)
      }
      # full slope vector at the POOLED estimates (fixed values kept)
      eq_free_idx <- lavpartable$free[eq$pt]
      beta_full <- lavpartable$ustart[eq$pt]
      beta_full[is.na(beta_full)] <- 0
      fp <- which(eq_free_idx > 0L)
      beta_full[fp] <- x[eq_free_idx[fp]]
      beta_full_list[[b]][[j]] <- beta_full
      m_star <- lav_sem_miiv_eq_dresid_dvech(
        sample_cov = cov_b, x_idx = sb$x_idx, y_idx = y_idx,
        i_idx = i_idx, beta = beta_full
      )
      m_free <- m_star[fp, , drop = FALSE]
      pw <- match(sb$gcol, free_slope_idx)
      mm <- rowsum(m_free, group = pw)
      pu <- sort(unique(pw))
      rhs[pu, cov_cols] <- rhs[pu, cov_cols] + mm
    }
  }

  # 3. d theta (KKT bordering for general linear constraints)
  if (!is.null(con) && nrow(con$jac) > 0L) {
    ncon <- nrow(con$jac)
    kkt <- rbind(
      cbind(dmat, t(con$jac)),
      cbind(con$jac, matrix(0, ncon, ncon))
    )
    dtheta <- solve(kkt, rbind(rhs, matrix(0, ncon, ntot)))
    dtheta <- dtheta[seq_len(npool), , drop = FALSE]
  } else {
    dtheta <- solve(dmat, rhs)
  }

  # 4. assemble the slope rows
  jac <- matrix(0, nrow = ntheta1, ncol = ntot)
  ridx <- match(free_slope_idx, free_directed_idx)
  ok <- which(!is.na(ridx))
  jac[ridx[ok], ] <- dtheta[ok, , drop = FALSE]

  # 5. intercept rows (meanstructure): per-equation rows in estimation
  #    order (last write wins, as in the estimation loop); equations with
  #    a slope block are then overwritten by the pooled nobs-weighted
  #    average (as in lav_sem_miiv_apply_directed_pool)
  if (meanstructure) {
    int_acc <- list() # fint -> list(num = row-sum, den = weight-sum)
    for (b in seq_len(nblocks)) {
      ov_names_b <- lavdata@ov.names[[b]]
      nvar_b <- lavmodel@nvar[b]
      mean_cols <- moff[b] + seq_len(nvar_b)
      for (j in seq_along(eqs[[b]])) {
        eq <- eqs[[b]][[j]]
        if (is.null(eq$ptint) || length(eq$ptint) != 1L) {
          next
        }
        fint <- lavpartable$free[eq$ptint]
        if (fint <= 0L) {
          next
        }
        rrow <- match(fint, free_directed_idx)
        if (is.na(rrow)) {
          next
        }
        y_idx <- match(eq$lhs_new, ov_names_b)
        sb <- eq$slope_block
        if (identical(eq$rhs_new, "1")) {
          # intercept-only equation: beta0 = ybar(y)
          row <- numeric(ntot)
          row[mean_cols[y_idx]] <- 1
          jac[rrow, ] <- row
        } else if (is.null(sb)) {
          # slopes all fixed: beta0 = ybar - xbar' b_fix
          x_idx <- match(eq$rhs_new, ov_names_b)
          b_fix <- lavpartable$ustart[eq$pt]
          b_fix[is.na(b_fix)] <- 0
          row <- numeric(ntot)
          row[mean_cols[y_idx]] <- 1
          row[mean_cols[x_idx]] <- row[mean_cols[x_idx]] - b_fix
          jac[rrow, ] <- row
        } else {
          # beta0 = ybar - xbar' beta(pooled):
          #   d = e_y - sum_l beta_l e_{x_l} - xbar' d beta
          x_idx <- sb$x_idx
          beta_full <- beta_full_list[[b]][[j]]
          eq_free_idx <- lavpartable$free[eq$pt]
          row <- numeric(ntot)
          row[mean_cols[y_idx]] <- 1
          row[mean_cols[x_idx]] <- row[mean_cols[x_idx]] - beta_full
          xbar <- sb$x_bar
          fp <- which(eq_free_idx > 0L)
          if (length(fp) > 0L) {
            dbeta <- dtheta[match(sb$gcol, free_slope_idx), ,
                            drop = FALSE]
            row <- row - drop(crossprod(dbeta, xbar[fp]))
          }
          w_e <- if (!is.null(eq$nobs)) eq$nobs else 1
          key <- as.character(fint)
          if (is.null(int_acc[[key]])) {
            int_acc[[key]] <- list(num = w_e * row, den = w_e)
          } else {
            int_acc[[key]]$num <- int_acc[[key]]$num + w_e * row
            int_acc[[key]]$den <- int_acc[[key]]$den + w_e
          }
        }
      }
    }
    for (key in names(int_acc)) {
      rrow <- match(as.integer(key), free_directed_idx)
      jac[rrow, ] <- int_acc[[key]]$num / int_acc[[key]]$den
    }
  }

  jac
}


# per-equation overidentification test based on Browne's residual-based
# statistic, using the asymptotic covariance (ACOV/Gamma) of the sample
# (co)variances or correlations. This replaces the standard Sargan test when
# the latter is not valid, in particular when the moments are polychoric
# correlations (categorical data), but also for non-normal continuous data.
#
# For an equation  y = beta' x + e  with instruments z (the moment conditions
# being cov(z, e) = 0), the overidentifying restrictions are
#   g = s_zy - S_zx %*% beta = 0   (a q-vector, q = #instruments)
# The (reduced) residual-based statistic is
#   T = N * g' [W - W D (D'W D)^{-1} D'W] g,   W = Gamma_g^{-1},  D = S_zx
# where Gamma_g is the ACOV of g, obtained from the ACOV of the relevant
# sample moments via the delta method (g is linear in those moments with
# d g_k / d s_{z_k y} = 1 and d g_k / d s_{z_k x_j} = -beta_j). This is
# algebraically equivalent to Browne's residual test applied to the single
# equation viewed as a covariance/correlation structure that is saturated
# except for the q - p overidentifying restrictions. Under the classical
# homoskedastic-normal weight Gamma_g = resvar * S_zz, T reduces exactly to
# the standard Sargan statistic (and D'W g = 0, so the projection drops out).
#
# Arguments:
#   y_idx, x_idx, i_idx : column indices (into 'sample_cov') of the outcome,
#                         the regressors x, and the instruments z
#   beta                : the (estimated) directed slopes (length = #x)
#   sample_cov          : the (saturated) covariance/correlation matrix over
#                         all variables in this block
#   nacov               : the ACOV of the sample statistics (Gamma); for
#                         categorical data this is the polychoric ACOV over
#                         [thresholds, correlations]
#   nth                 : number of leading entries in 'nacov' that precede the
#                         (co)variance/correlation block (e.g. thresholds);
#                         0 for the continuous case
#   moment_idx_mat      : symmetric integer matrix mapping a variable pair
#                         (a, b) to its position within the moment block of
#                         'nacov' (vech order, off-diagonal for correlations)
#   nobs                : number of observations for this block
lav_sem_miiv_browne_test <- function(y_idx = integer(0L), x_idx = integer(0L),
                                     i_idx = integer(0L), beta = numeric(0L),
                                     sample_cov = NULL, nacov = NULL,
                                     nth = 0L, moment_idx_mat = NULL,
                                     nobs = 0L) {
  q <- length(i_idx)
  p <- length(x_idx)
  df <- q - p
  out <- c(stat = as.numeric(NA), df = df, pvalue = as.numeric(NA))
  if (df <= 0L) {
    return(out)
  }

  s_zy <- sample_cov[i_idx, y_idx, drop = FALSE]
  s_zx <- sample_cov[i_idx, x_idx, drop = FALSE]
  g <- drop(s_zy - s_zx %*% beta)

  # indices into 'nacov' of the moments that enter g: the q instrument-outcome
  # moments, followed by the q * p instrument-regressor moments
  p_y <- nth + moment_idx_mat[cbind(i_idx, rep.int(y_idx, q))]
  p_x <- nth + as.vector(moment_idx_mat[cbind(rep(i_idx, times = p),
                                              rep(x_idx, each = q))])
  p_all <- c(p_y, p_x)
  gam_p <- nacov[p_all, p_all, drop = FALSE]

  # delta-method Jacobian of g w.r.t. those moments
  jmat <- matrix(0, q, length(p_all))
  jmat[seq_len(q), seq_len(q)] <- diag(q)
  for (j in seq_len(p)) {
    cols <- q + (j - 1L) * q + seq_len(q)
    jmat[cbind(seq_len(q), cols)] <- -beta[j]
  }
  gamma_g <- jmat %*% gam_p %*% t(jmat)

  # T = N g' [W - W D (D'W D)^{-1} D'W] g
  w <- lav_mat_sym_inverse(gamma_g)
  wg <- drop(w %*% g)
  q1 <- sum(g * wg)
  q2 <- 0
  if (p > 0L) {
    dtwg <- crossprod(s_zx, wg)
    dtwd <- crossprod(s_zx, w %*% s_zx)
    q2 <- drop(crossprod(dtwg, lav_mat_sym_solve_spd(dtwd, dtwg)))
  }
  stat <- as.numeric(nobs * (q1 - q2))

  out["stat"] <- stat
  out["pvalue"] <- pchisq(stat, df, lower.tail = FALSE)
  out
}
