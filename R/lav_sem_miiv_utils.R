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
