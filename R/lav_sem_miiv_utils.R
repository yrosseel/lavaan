# varcov solution is:
# theta_2 = solve(t(Delta2) %*% W2 %*% Delta2) %*% t(Delta2) %*% W2 %*% s_vech

# we need to compute the Jacobian of theta_2 wrt the elements of s_vech
#
# four options:
# - ULS
# - GLS
# - 2RLS
# - RLS

# this is the GLS version where W2 = t(D) %*% (S.inv %x% S.inv) %*% D
lav_sem_miiv_utils_jacb_gls <- function(sample_cov = NULL, delta2 = NULL) {
  nvar <- nrow(sample_cov)
  s_vech <- lav_matrix_vech(sample_cov)
  s_inv <- solve(sample_cov)

  # only covariance block
  meanstructure.flag <- FALSE
  if (nrow(delta2) > length(s_vech)) {
    meanstructure.flag <- TRUE
    delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
  }
  w2 <- lav_matrix_duplication_pre_post(s_inv %x% s_inv)

  # shortcut
  M <- t(delta2) %*% w2 %*% delta2

  # result
  theta2 <- solve(M, t(delta2) %*% w2 %*% s_vech)

  # error term
  e_mat <- lav_matrix_vech_reverse(s_vech - delta2 %*% theta2)

  # correction term
  q_mat <- s_inv %*% e_mat %*% s_inv
  correction <-
    lav_matrix_duplication_pre_post(q_mat %x% s_inv + s_inv %x% q_mat)

  out <- solve(M, t(delta2)) %*% (w2 - correction)

  # meanstructure?
  if (meanstructure.flag) {
    out <- cbind(matrix(0, nrow = nrow(out), ncol = nvar), out)
  }

  out
}

# 2RLS solution, where W2 = t(D) %*% (Sigma.inv %x% Sigma.inv) %*% D
# and Sigma is based on the ULS estimate of theta2
lav_sem_miiv_utils_jacb_2rls <- function(sample_cov = NULL, delta2 = NULL) {
  nvar <- nrow(sample_cov)
  s_vech <- lav_matrix_vech(sample_cov)

  # only covariance block
  meanstructure.flag <- FALSE
  if (nrow(delta2) > length(s_vech)) {
    meanstructure.flag <- TRUE
    delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
  }
  ip <- diag(nvar)
  w_uls <- 0.5 * lav_matrix_duplication_pre_post(ip %x% ip)

  # Step 1: initial theta2, Sigma
  m_uls <- t(delta2) %*% w_uls %*% delta2
  theta2 <- solve(m_uls, t(delta2) %*% w_uls %*% s_vech)
  sigma <- lav_matrix_vech_reverse(delta2 %*% theta2)
  sigma_inv <- solve(sigma)

  # Step 2: W2, M, theta final
  w2 <- lav_matrix_duplication_pre_post(sigma_inv %x% sigma_inv)
  m_sigma <- t(delta2) %*% w2 %*% delta2
  theta <- solve(m_sigma, t(delta2) %*% w2 %*% s_vech)

  # Step 3: correction matrix c_mat
  e_mat <- lav_matrix_vech_reverse(s_vech - delta2 %*% theta)
  q_mat <- sigma_inv %*% e_mat %*% sigma_inv
  c_mat <-
    lav_matrix_duplication_pre_post(q_mat %x% sigma_inv + sigma_inv %x% q_mat)

  # Step 4: Jacobian
  J_theta2 <- solve(m_uls, t(delta2) %*% w_uls)
  out <- solve(m_sigma, t(delta2)) %*% (w2 - c_mat %*% delta2 %*% J_theta2)

  # meanstructure?
  if (meanstructure.flag) {
    out <- cbind(matrix(0, nrow = nrow(out), ncol = nvar), out)
  }

  out
}

# RLS solution
lav_sem_miiv_utils_jacb_rls <- function(sample_cov = NULL, delta2 = NULL) {
  nvar <- nrow(sample_cov)
  s_vech <- lav_matrix_vech(sample_cov)

  # only covariance block
  meanstructure.flag <- FALSE
  if (nrow(delta2) > length(s_vech)) {
    meanstructure.flag <- TRUE
    delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
  }
  ip <- diag(nvar)
  w_uls <- 0.5 * lav_matrix_duplication_pre_post(ip %x% ip)

  # step 1
  m_uls <- t(delta2) %*% w_uls %*% delta2
  theta2_uls <- solve(m_uls, t(delta2) %*% w_uls %*% s_vech)

  # step 2 (iterative)
  theta2 <- theta2_uls
  for (i in seq_len(200L)) {
    old_x <- theta2
    new_sigma <- lav_matrix_vech_reverse(delta2 %*% theta2)
    sigma_inv <- solve(new_sigma)
    w2 <- 0.5 * lav_matrix_duplication_pre_post(sigma_inv %x% sigma_inv)
    m_mat <- t(delta2) %*% w2 %*% delta2
    theta2 <- drop(solve(m_mat, t(delta2) %*% w2 %*% s_vech))
    if (sum((old_x - theta2)^2) < 1e-12 * (1 + sum(theta2^2))) break
  }

  # final quantities
  new_sigma <- lav_matrix_vech_reverse(delta2 %*% theta2)
  sigma_inv <- solve(new_sigma)
  w2 <- 0.5 * lav_matrix_duplication_pre_post(sigma_inv %x% sigma_inv)
  m_mat <- t(delta2) %*% w2 %*% delta2
  m_mat_inv <- solve(m_mat)
  e2 <- s_vech - delta2 %*% theta2
  e_mat <- lav_matrix_vech_reverse(e2)
  sigma_inv_e <- sigma_inv %*% e_mat
  minv_delta2t <- m_mat_inv %*% t(delta2)

  # helper: compute (dW)e2 given d vech(Sigma) as a vector v
  dW_e2_from_v <- function(v) {
    dSigma <- lav_matrix_vech_reverse(v)
    d_sigma_inv <- -sigma_inv %*% dSigma %*% sigma_inv
    sym_part <- sigma_inv_e %*% d_sigma_inv + d_sigma_inv %*% t(sigma_inv_e)
    dw2_e2 <-
      0.5 * drop(lav_matrix_duplication_pre(as.matrix(as.vector(sym_part))))
    dw2_e2
  }

  # A = M2^{-1} delta2^T J_W^(theta2): q2 x q2
  # column j of J_W^(theta2): dSigma = unvech(delta2[:,j])
  JW_theta2 <- matrix(0.0, nrow(delta2), ncol(delta2))
  for (j in seq_len(ncol(delta2))) {
    JW_theta2[, j] <- dW_e2_from_v(delta2[, j])
  }
  A <- minv_delta2t %*% JW_theta2
  IminusA_inv <- solve(diag(ncol(delta2)) - A)

  # Jacobian: (I - A)^{-1} M2^{-1} Delta2^T W
  # (for fixed W case this would just be M2^{-1} Delta2^T W = M2inv_D2t %*% W2)
  out <- IminusA_inv %*% minv_delta2t %*% w2

  # meanstructure?
  if (meanstructure.flag) {
    out <- cbind(matrix(0, nrow = nrow(out), ncol = nvar), out)
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
lav_sem_miiv_utils_jaca_uls_gls <- function(lavmodel = NULL,
                                            lavpartable = NULL,
                                            lavh1 = NULL,
                                            free.directed.idx = integer(0L),
                                            free.undirected.idx = integer(0L),
                                            iv.varcov.method = "ULS") {
  nblocks <- lavmodel@nblocks

  # delta across all blocks
  delta_block <- lav_model_delta(
    lavmodel = lavmodel,
    ceq.simple = lavmodel@ceq.simple.only
  )

  # augment lavpartable to include model matrices and row/col indices
  mm_info <- lav_lisrel(lavpartable, target = NULL, extra = FALSE)

  jac_a <- matrix(0, length(free.undirected.idx), length(free.directed.idx))
  for (b in seq_len(nblocks)) {
    # s_vech
    sample_cov <- lavh1$implied$cov[[b]]
    nvar <- nrow(sample_cov)
    s_vech <- lav_matrix_vech(sample_cov)

    # delta2
    delta2 <- delta_block[[b]][, free.undirected.idx, drop = FALSE]
    if (lavmodel@meanstructure) {
      delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
    }

    # w2
    if (iv.varcov.method == "GLS") {
      s_inv <- solve(sample_cov)
    } else {
      s_inv <- diag(1, nrow = nrow(sample_cov))
    }
    w2_22 <- 0.5 * lav_matrix_duplication_pre_post(s_inv %x% s_inv)
    w2 <- w2_22

    # MLIST
    mm.in.group <- seq_len(lavmodel@nmat[b]) + cumsum(c(0, lavmodel@nmat))[b]
    MLIST <- lavmodel@GLIST[mm.in.group]
    if (!is.null(MLIST$beta)) {
      MLIST$IB.inv <- lav_lisrel_ibinv(MLIST)
    }

    m_mat <- t(delta2) %*% w2 %*% delta2
    theta2 <- solve(m_mat, t(delta2) %*% w2 %*% s_vech)
    m_inv <- solve(m_mat)
    w2e <- w2 %*% (s_vech - delta2 %*% theta2)
    delta2tw2 <- t(delta2) %*% w2

    # container for Jacobian for this block
    for (k in seq_along(free.directed.idx)) {
      d_idx <- match(free.directed.idx[k], lavpartable$free)[1]
      # correct block?
      if (lavpartable$block[d_idx] != b) {
        next
      }
      d_mat <- mm_info$mat[d_idx]
      d_row <- mm_info$row[d_idx]
      d_col <- mm_info$col[d_idx]
      dDeltak <- matrix(0, nrow(delta2), ncol(delta2))
      # we only need to fill the columns of dDeltak that correspond to the
      # psi elements in free.undirected.idx
      for (i in seq_along(free.undirected.idx)) {
        un_idx <- match(free.undirected.idx[i], lavpartable$free)[1]
        un_mat <- mm_info$mat[un_idx]
        if (un_mat == "theta") {
          next
        }
        un_row <- mm_info$row[un_idx]
        un_col <- mm_info$col[un_idx]

        if (d_mat == "lambda") {
          tmp <- lav_lisrel_d2sigma_lambda_psi(
            MLIST = MLIST,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        } else if (d_mat == "beta") {
          tmp <- lav_lisrel_d2sigma_beta_psi(
            MLIST = MLIST,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        }
        dDeltak[, i] <- lav_matrix_vech(tmp)
      }
      jac_a[, k] <-
        m_inv %*% (t(dDeltak) %*% w2e - delta2tw2 %*% dDeltak %*% theta2)
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
                                         free.directed.idx = integer(0L),
                                         free.undirected.idx = integer(0L)) {
  nblocks <- lavmodel@nblocks

  # delta across all blocks
  delta_block <- lav_model_delta(
    lavmodel = lavmodel,
    ceq.simple = lavmodel@ceq.simple.only
  )

  # augment lavpartable to include model matrices and row/col indices
  mm_info <- lav_lisrel(lavpartable, target = NULL, extra = FALSE)

  # first get jac_a for ULS (step 1)
  jac_uls <- lav_sem_miiv_utils_jaca_uls_gls(
    lavmodel = lavmodel,
    lavpartable = lavpartable, lavh1 = lavh1,
    free.directed.idx = free.directed.idx,
    free.undirected.idx = free.undirected.idx, iv.varcov.method = "ULS"
  )

  jac_a <- matrix(0, length(free.undirected.idx), length(free.directed.idx))
  for (b in seq_len(nblocks)) {
    # s_vech
    sample_cov <- lavh1$implied$cov[[b]]
    nvar <- nrow(sample_cov)
    s_vech <- lav_matrix_vech(sample_cov)

    # delta2
    delta2 <- delta_block[[b]][, free.undirected.idx, drop = FALSE]
    if (lavmodel@meanstructure) {
      delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
    }
    # w2
    s_inv <- diag(1, nrow = nrow(sample_cov))
    w2_uls <- 0.5 * lav_matrix_duplication_pre_post(s_inv %x% s_inv)

    # MLIST
    mm.in.group <- seq_len(lavmodel@nmat[b]) + cumsum(c(0, lavmodel@nmat))[b]
    MLIST <- lavmodel@GLIST[mm.in.group]
    if (!is.null(MLIST$beta)) {
      MLIST$IB.inv <- lav_lisrel_ibinv(MLIST)
    }

    # step 1
    m_uls <- t(delta2) %*% w2_uls %*% delta2
    m_uls_inv <- solve(m_uls)
    theta_uls <- drop(m_uls_inv %*% t(delta2) %*% w2_uls %*% s_vech)

    # step 2
    new_sigma <- lav_matrix_vech_reverse(delta2 %*% theta_uls)
    sigma_inv <- solve(new_sigma)
    w2 <- 0.5 * lav_matrix_duplication_pre_post(sigma_inv %x% sigma_inv)
    m_mat <- t(delta2) %*% w2 %*% delta2
    m_mat_inv <- solve(m_mat)
    theta2 <- drop(m_mat_inv %*% t(delta2) %*% w2 %*% s_vech)
    e2 <- s_vech - delta2 %*% theta2
    e_mat <- lav_matrix_vech_reverse(e2)

    # pre-compute
    minv_delta2t <- m_mat_inv %*% t(delta2)
    w2_e2 <- w2 %*% e2
    delta2t_w2 <- t(delta2) %*% w2
    sigma_inv_e <- sigma_inv %*% e_mat

    for (k in seq_along(free.directed.idx)) {
      d_idx <- match(free.directed.idx[k], lavpartable$free)[1]
      # correct block?
      if (lavpartable$block[d_idx] != b) {
        next
      }
      d_mat <- mm_info$mat[d_idx]
      d_row <- mm_info$row[d_idx]
      d_col <- mm_info$col[d_idx]
      dDeltak <- matrix(0, nrow(delta2), ncol(delta2))
      # we only need to fill the columns of dDeltak that correspond to the
      # psi elements in free.undirected.idx
      for (i in seq_along(free.undirected.idx)) {
        un_idx <- match(free.undirected.idx[i], lavpartable$free)[1]
        un_mat <- mm_info$mat[un_idx]
        if (un_mat == "theta") {
          next
        }
        un_row <- mm_info$row[un_idx]
        un_col <- mm_info$col[un_idx]

        if (d_mat == "lambda") {
          tmp <- lav_lisrel_d2sigma_lambda_psi(
            MLIST = MLIST,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        } else if (d_mat == "beta") {
          tmp <- lav_lisrel_d2sigma_beta_psi(
            MLIST = MLIST,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        }
        dDeltak[, i] <- lav_matrix_vech(tmp)
      }

      # d vech(Sigma)/d theta1[k]: product rule on unvech(Delta2 theta_uls)
      v_k <- dDeltak %*% theta_uls + delta2 %*% jac_uls[, k]
      dSigma_k <- lav_matrix_vech_reverse(v_k)

      # d_sigma_inv = -sigma.inv dSigma_k sigma.inv
      d_sigma_inv <- -sigma_inv %*% dSigma_k %*% sigma_inv

      # the tricky one:
      # (dw2_k) e2 = (1/2)
      # D^T vec(sigma_inv e_mat dsigma_inv_k + dsigma_invk e_mat sigma_inv)
      sym_part <- sigma_inv_e %*% d_sigma_inv + d_sigma_inv %*% t(sigma_inv_e)
      dw2_e2 <-
        0.5 * drop(lav_matrix_duplication_pre(as.matrix(as.vector(sym_part))))

      # term 1: direct delta2 variation in final WLS
      term1 <-
        m_mat_inv %*% (t(dDeltak) %*% w2_e2 - delta2t_w2 %*% dDeltak %*% theta2)

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
                                        free.directed.idx = integer(0L),
                                        free.undirected.idx = integer(0L)) {
  nblocks <- lavmodel@nblocks

  # delta across all blocks
  delta_block <- lav_model_delta(
    lavmodel = lavmodel,
    ceq.simple = lavmodel@ceq.simple.only
  )

  # augment lavpartable to include model matrices and row/col indices
  mm_info <- lav_lisrel(lavpartable, target = NULL, extra = FALSE)

  jac_a <- matrix(0, length(free.undirected.idx), length(free.directed.idx))
  for (b in seq_len(nblocks)) {
    # s_vech
    sample_cov <- lavh1$implied$cov[[b]]
    nvar <- nrow(sample_cov)
    s_vech <- lav_matrix_vech(sample_cov)

    # delta2
    delta2 <- delta_block[[b]][, free.undirected.idx, drop = FALSE]
    if (lavmodel@meanstructure) {
      delta2 <- delta2[-seq_len(nvar), , drop = FALSE]
    }

    # w2
    s_inv <- diag(1, nrow = nrow(sample_cov))
    w2_uls <- 0.5 * lav_matrix_duplication_pre_post(s_inv %x% s_inv)

    # MLIST
    mm.in.group <- seq_len(lavmodel@nmat[b]) + cumsum(c(0, lavmodel@nmat))[b]
    MLIST <- lavmodel@GLIST[mm.in.group]
    if (!is.null(MLIST$beta)) {
      MLIST$IB.inv <- lav_lisrel_ibinv(MLIST)
    }

    # step 1
    m_uls <- t(delta2) %*% w2_uls %*% delta2
    m_uls_inv <- solve(m_uls)
    theta_uls <- drop(m_uls_inv %*% t(delta2) %*% w2_uls %*% s_vech)

    # step 2
    theta2 <- theta_uls
    for (i in seq_len(200)) {
      old_x <- theta2
      new_sigma <- lav_matrix_vech_reverse(delta2 %*% theta2)
      sigma_inv <- solve(new_sigma)
      w2 <- 0.5 * lav_matrix_duplication_pre_post(sigma_inv %x% sigma_inv)
      m_mat <- t(delta2) %*% w2 %*% delta2
      theta2 <- drop(solve(m_mat, t(delta2) %*% w2 %*% s_vech))
      if (sum((old_x - theta2)^2) < 1e-12 * (1 + sum(theta2^2))) break
    }

    # final quantities
    new_sigma <- lav_matrix_vech_reverse(delta2 %*% theta2)
    sigma_inv <- solve(new_sigma)
    w2 <- 0.5 * lav_matrix_duplication_pre_post(sigma_inv %x% sigma_inv)
    m_mat <- t(delta2) %*% w2 %*% delta2
    m_mat_inv <- solve(m_mat)
    e2 <- s_vech - delta2 %*% theta2
    e_mat <- lav_matrix_vech_reverse(e2)
    sigma_inv_e <- sigma_inv %*% e_mat

    minv_delta2t <- m_mat_inv %*% t(delta2)
    w2_e2 <- w2 %*% e2
    delta2t_w2 <- t(delta2) %*% w2

    # helper: compute (dW)e2 given d vech(Sigma) as a vector v
    dW_e2_from_v <- function(v) {
      dSigma <- lav_matrix_vech_reverse(v)
      d_sigma_inv <- -sigma_inv %*% dSigma %*% sigma_inv
      sym_part <- sigma_inv_e %*% d_sigma_inv + d_sigma_inv %*% t(sigma_inv_e)
      dw2_e2 <-
        0.5 * drop(lav_matrix_duplication_pre(as.matrix(as.vector(sym_part))))
      dw2_e2
    }

    # A = M2^{-1} delta2^T J_W^(theta2): q2 x q2
    # column j of J_W^(theta2): dSigma = unvech(delta2[:,j])
    JW_theta2 <- matrix(0.0, nrow(delta2), ncol(delta2))
    for (j in seq_len(ncol(delta2))) {
      JW_theta2[, j] <- dW_e2_from_v(delta2[, j])
    }
    A <- minv_delta2t %*% JW_theta2
    IminusA_inv <- solve(diag(ncol(delta2)) - A)

    for (k in seq_along(free.directed.idx)) {
      d_idx <- match(free.directed.idx[k], lavpartable$free)[1]
      # correct block?
      if (lavpartable$block[d_idx] != b) {
        next
      }
      d_mat <- mm_info$mat[d_idx]
      d_row <- mm_info$row[d_idx]
      d_col <- mm_info$col[d_idx]
      dDeltak <- matrix(0, nrow(delta2), ncol(delta2))
      # we only need to fill the columns of dDeltak that correspond to the
      # psi elements in free.undirected.idx
      for (i in seq_along(free.undirected.idx)) {
        un_idx <- match(free.undirected.idx[i], lavpartable$free)[1]
        un_mat <- mm_info$mat[un_idx]
        if (un_mat == "theta") {
          next
        }
        un_row <- mm_info$row[un_idx]
        un_col <- mm_info$col[un_idx]

        if (d_mat == "lambda") {
          tmp <- lav_lisrel_d2sigma_lambda_psi(
            MLIST = MLIST,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        } else if (d_mat == "beta") {
          tmp <- lav_lisrel_d2sigma_beta_psi(
            MLIST = MLIST,
            i = d_row, j = d_col, k = un_row, l = un_col
          )
        }
        dDeltak[, i] <- lav_matrix_vech(tmp)
      }

      # term 1: direct delta2 variation in final WLS
      term1 <-
        m_mat_inv %*% (t(dDeltak) %*% w2_e2 - delta2t_w2 %*% dDeltak %*% theta2)

      # term 2: W variation from dSigma^(Delta2) = unvech(dDelta2_k %*% theta2)
      # (the dSigma^(theta2) part is handled implicitly by (I-A)^{-1})
      v_k <- as.vector(dDeltak %*% theta2)
      term2 <- minv_delta2t %*% dW_e2_from_v(v_k)

      jac_a[, k] <- IminusA_inv %*% (term1 + term2)
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
                                        free.directed.idx = integer(0L)) {
  # continuous data only (for now)
  nvar <- lavmodel@nvar[block]
  pstar <- nvar * (nvar + 1L) / 2
  ntheta1 <- length(free.directed.idx)

  b <- block

  # K matrix
  if (lavmodel@meanstructure) {
    K_mat <- matrix(0.0, nrow = ntheta1, ncol = (nvar + pstar))
  } else {
    K_mat <- matrix(0.0, nrow = ntheta1, ncol = pstar)
  }

  # collect k_mat matrices from each equation
  for (j in seq_along(eqs[[b]])) {
    # this equation
    eq <- eqs[[b]][[j]]
    free.idx <- match(lavpartable$free[eq$pt], free.directed.idx)
    nx <- nrow(eq$k_mat)
    if (nx > 0L) {
      if (all(free.idx > 0L)) {
        K_mat[free.idx,] <- eq$k_mat
      } else {
        # remove non-free elements
        zero.idx <- which(free.idx == 0L)
        free.idx <- free.idx[-zero.idx]
        if (length(free.idx) > 0) {
          K_mat[free.idx,] <- eq$k_mat[-zero.idx, , drop = FALSE]
        }
      }
    }
    if (lavmodel@meanstructure) {
      free.int.idx <- match(lavpartable$free[eq$ptint], free.directed.idx)
      if (free.int.idx > 0L) {
        K_mat[free.int.idx,] <- eq$k_mat_int
      }
    }
  } # eq

  K_mat
}
