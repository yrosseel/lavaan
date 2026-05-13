# RAM representation
#
# initial version: YR 2021-10-04

lav_ram <- function(partable = NULL,
                    target = NULL,
                    extra = FALSE,
                    remove_nonexisting = TRUE) {
  # prepare target list
  if (is.null(target)) target <- partable

  stopifnot(!is.null(target$block))

  # not for categorical data (yet)
  if (any(partable$op == "|")) {
    lav_msg_stop(gettext("RAM representation is not (yet) supported for
                         categorical endogenous variables."))
  }

  # not for conditional.x = TRUE yet
  conditional_x <- any(partable$exo > 0L & partable$op == "~")
  if (conditional_x) {
    lav_msg_stop(gettext("RAM representation is not (yet) supported
                         if conditional.x = TRUE"))
  }

  # prepare output
  n <- length(target$lhs)
  tmp_mat <- character(n)
  tmp_row <- integer(n)
  tmp_col <- integer(n)

  # global settings
  meanstructure <- any(partable$op == "~1")
  # categorical <- any(partable$op == "|")
  # group_w_free <- any(partable$lhs == "group" & partable$op == "%")

  # number of blocks
  nblocks <- lav_partable_nblocks(partable)

  # always return ov.idx
  ov_idx <- vector("list", nblocks)
  ov_dummy_names_nox <- vector("list", nblocks)
  ov_dummy_names_x <- vector("list", nblocks)
  if (extra) {
    rep_mm_names <- vector("list", nblocks)
    rep_mm_number <- vector("list", nblocks)
    rep_mm_rows <- vector("list", nblocks)
    rep_mm_cols <- vector("list", nblocks)
    rep_mm_dim_names <- vector("list", nblocks)
    rep_mm_symmetric <- vector("list", nblocks)
  }

  for (g in 1:nblocks) {
    # info from user model per block
    ov_names <- lav_partable_vnames(partable, "ov", block = g)
    nvar <- length(ov_names)
    ov_idx[[g]] <- seq_len(nvar)
    ov_dummy_names_nox[[g]] <- character(0)
    ov_dummy_names_x[[g]] <- character(0)

    lv_names <- lav_partable_vnames(partable, "lv", block = g)
    both_names <- c(ov_names, lv_names)
    nboth <- length(both_names)

    # 1. "=~" indicators
    idx <- which(target$block == g & target$op == "=~")
    tmp_mat[idx] <- "A"
    tmp_row[idx] <- match(target$rhs[idx], both_names)
    tmp_col[idx] <- match(target$lhs[idx], both_names)

    # 2. "~" regressions
    idx <- which(target$block == g & (target$op == "~" |
      target$op == "<~"))
    tmp_mat[idx] <- "A"
    tmp_row[idx] <- match(target$lhs[idx], both_names)
    tmp_col[idx] <- match(target$rhs[idx], both_names)

    # 3. "~~" variances/covariances
    idx <- which(target$block == g & target$op == "~~")
    tmp_mat[idx] <- "S"
    tmp_row[idx] <- match(target$lhs[idx], both_names)
    tmp_col[idx] <- match(target$rhs[idx], both_names)

    # catch lower-elements in theta/psi
    idx_lower <- which(tmp_mat == "S" & tmp_row > tmp_col)
    if (length(idx_lower) > 0L) {
      tmp <- tmp_row[idx_lower]
      tmp_row[idx_lower] <- tmp_col[idx_lower]
      tmp_col[idx_lower] <- tmp
    }

    # 4. "~1" means/intercepts
    idx <- which(target$block == g & target$op == "~1")
    tmp_mat[idx] <- "m"
    tmp_row[idx] <- match(target$lhs[idx], both_names)
    tmp_col[idx] <- 1L

    # 5. "|" th
    # not used yet

    # 6. "~*~" scales
    # not used yet

    # 7. group weights
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "%")
    tmp_mat[idx] <- "gw"
    tmp_row[idx] <- 1L
    tmp_col[idx] <- 1L

    # 8. instruments
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "|~")
    tmp_mat[idx] <- "miiv"
    tmp_row[idx] <- 0L
    tmp_col[idx] <- 0L

    if (extra) {
      # mRows
      mm_rows <- list(
        ov.idx = 1L,
        A = nboth,
        S = nboth,
        m = nboth,
        gw = 1L
      )

      # mCols
      mm_cols <- list(
        ov.idx = nvar,
        A = nboth,
        S = nboth,
        m = 1L,
        gw = 1L
      )

      # dimNames for LISREL model matrices
      mm_dim_names <- list(
        ov.idx = list("ov.idx", ov_names),
        A = list(both_names, both_names),
        S = list(both_names, both_names),
        m = list(both_names, "intercept"),
        gw = list("group", "weight")
      )
      # isSymmetric
      mm_symmetric <- list(
        ov.idx = FALSE,
        A = FALSE,
        S = TRUE,
        m = FALSE,
        gw = FALSE
      )

      # which mm's do we need? (always include ov.idx, A and S)
      idx_1 <- which(target$block == g)
      mm_names <- c("ov.idx", "A", "S")
      if (meanstructure) {
        mm_names <- c(mm_names, "m")
      }
      if ("gw" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "gw")
      }

      rep_mm_names[[g]] <- mm_names
      rep_mm_number[[g]] <- length(mm_names)
      rep_mm_rows[[g]] <- unlist(mm_rows[mm_names])
      rep_mm_cols[[g]] <- unlist(mm_cols[mm_names])
      rep_mm_dim_names[[g]] <- mm_dim_names[mm_names]
      rep_mm_symmetric[[g]] <- unlist(mm_symmetric[mm_names])
    } # extra
  } # nblocks

  rep_1 <- list(
    mat = tmp_mat,
    row = tmp_row,
    col = tmp_col
  )

  # always return ov.idx attribute
  attr(rep_1, "ov.idx") <- ov_idx
  attr(rep_1, "ov.dummy.names.nox") <- ov_dummy_names_nox
  attr(rep_1, "ov.dummy.names.x") <- ov_dummy_names_x

  if (extra) {
    attr(rep_1, "mmNames") <- rep_mm_names
    attr(rep_1, "mmNumber") <- rep_mm_number
    attr(rep_1, "mmRows") <- rep_mm_rows
    attr(rep_1, "mmCols") <- rep_mm_cols
    attr(rep_1, "mmDimNames") <- rep_mm_dim_names
    attr(rep_1, "mmSymmetric") <- rep_mm_symmetric
  }

  rep_1
}

# the model-implied variance/covariance matrix of the observed variables
lav_ram_sigmahat <- function(mlist = NULL, delta = NULL) {
  ov_idx <- as.integer(mlist$ov.idx[1, ])
  a <- mlist$A
  s <- mlist$S

  # get (I-A)^{-1}
  ia_inv <- lav_matrix_inverse_iminus(a)

  # compute Sigma for all ov and lv
  vyeta <- tcrossprod(ia_inv %*% s, ia_inv)

  # select only observed part
  vy <- vyeta[ov_idx, ov_idx, drop = FALSE]

  # if delta, scale
  if (!is.null(mlist$delta) && delta) {
    nvar <- ncol(vy)
    mm_delta <- diag(mlist$delta[, 1L], nrow = nvar, ncol = nvar)
    vy <- mm_delta %*% vy %*% mm_delta
  }

  vy
}

# VETA: the variance/covariance matrix of the latent variables only
lav_ram_veta <- function(mlist = NULL) {
  ov_idx <- as.integer(mlist$ov.idx[1, ])
  a <- mlist$A
  s <- mlist$S

  # get (I-A)^{-1}
  ia_inv <- lav_matrix_inverse_iminus(a)

  # compute Sigma for all ov and lv
  vyeta <- tcrossprod(ia_inv %*% s, ia_inv)

  # select only latent part
  veta <- vyeta[-ov_idx, -ov_idx, drop = FALSE]

  veta
}

# MuHat: the model-implied means/intercepts
lav_ram_muhat <- function(mlist = NULL) {
  ov_idx <- as.integer(mlist$ov.idx[1, ])
  a <- mlist$A
  m <- mlist$m

  # shortcut
  if (is.null(m)) {
    return(matrix(0, nrow = length(ov_idx), 1L))
  }

  # get (I-A)^{-1}
  ia_inv <- lav_matrix_inverse_iminus(a)

  # all means/intercepts
  eyeta <- ia_inv %*% m

  # select observed only
  muhat <- eyeta[ov_idx, , drop = FALSE]

  muhat
}

lav_ram_implied_fast <- function(mlist = NULL,
                                 need_sigma = FALSE,
                                 need_mu = FALSE,
                                 delta = TRUE) {
  ov_idx <- as.integer(mlist$ov.idx[1, ])
  a <- mlist$A
  s <- mlist$S

  out <- list()

  if (need_sigma || need_mu) {
    ia_inv <- lav_matrix_inverse_iminus(a)
  }

  if (need_sigma) {
    # compute Sigma for all ov and lv
    vyeta <- tcrossprod(ia_inv %*% s, ia_inv)

    # select only observed part
    vy <- vyeta[ov_idx, ov_idx, drop = FALSE]

    # if delta, scale
    if (!is.null(mlist$delta) && delta) {
      nvar <- ncol(vy)
      mm_delta <- diag(mlist$delta[, 1L], nrow = nvar, ncol = nvar)
      vy <- mm_delta %*% vy %*% mm_delta
    }

    out$sigma <- vy
  }

  if (need_mu) {
    m <- mlist$m

    # shortcut
    if (is.null(m)) {
      out$mu <- matrix(0, nrow = length(ov_idx), 1L)
    } else {
      # all means/intercepts
      eyeta <- ia_inv %*% m

      # select observed only
      out$mu <- eyeta[ov_idx, , drop = FALSE]
    }
  }

  out
}

# derivative of 'Sigma' wrt the (freel) elements in A and/or S
lav_ram_dsigma <- function(m = "A",
                           idx = seq_along(mlist[[m]]),
                           mlist = NULL,
                           vech = TRUE) {
  ov_idx <- as.integer(mlist$ov.idx[1, ])
  a <- mlist$A
  s <- mlist$S

  nvar <- length(ov_idx)
  nboth <- nrow(a)

  # shortcut for ov.idx, m, ...
  if (!m %in% c("A", "S")) {
    pstar <- nvar * (nvar + 1) / 2
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # get (I-A)^{-1}
  ia_inv <- lav_matrix_inverse_iminus(a)

  if (m == "A") {
    l1 <- (ia_inv %*% s %*% t(ia_inv))[ov_idx, , drop = FALSE]
    kol_idx <- matrix(1:(nboth * nboth), nboth, nboth, byrow = TRUE)[idx]
    dx <- (l1 %x% ia_inv[ov_idx, , drop = FALSE])[, idx, drop = FALSE] +
      (ia_inv[ov_idx, , drop = FALSE] %x% l1)[, kol_idx, drop = FALSE]
    # this is not really needed (because we select idx=m.el.idx)
    # but just in case we need all elements of beta...
    dx[, which(idx %in% lav_matrix_diag_idx(nboth))] <- 0.0
  } else if (m == "S") {
    dx <- (ia_inv[ov_idx, , drop = FALSE] %x% ia_inv[ov_idx, , drop = FALSE])
    # symmetry correction, but keeping all duplicated elements
    # since we depend on idx=m.el.idx
    lower_idx <- lav_matrix_vech_idx(nboth, diagonal = FALSE)
    upper_idx <- lav_matrix_vechru_idx(nboth, diagonal = FALSE)
    offdiag_sum <- dx[, lower_idx] + dx[, upper_idx]
    dx[, c(lower_idx, upper_idx)] <- cbind(offdiag_sum, offdiag_sum)
    dx <- dx[, idx, drop = FALSE]
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  # vech?
  if (vech) {
    v_idx <- lav_matrix_vech_idx(nvar)
    dx <- dx[v_idx, , drop = FALSE]
  }

  dx
}

# derivative of 'Mu' wrt the (free) elements in A and/or m
lav_ram_dmu <- function(m = "A",
                        idx = seq_along(mlist[[m]]),
                        mlist = NULL,
                        vech = TRUE) {
  ov_idx <- as.integer(mlist$ov.idx[1, ])
  a <- mlist$A
  # s <- mlist$S

  nvar <- length(ov_idx)
  # nboth <- nrow(a)

  # shortcut for ov.idx, m, ...
  if (!m %in% c("A", "m")) {
    return(matrix(0.0, nrow = nvar, ncol = length(idx)))
  }

  # get (I-A)^{-1}
  ia_inv <- lav_matrix_inverse_iminus(a)

  if (m == "A") {
    dx <- (t(ia_inv %*% mlist$m) %x% ia_inv)[ov_idx, idx, drop = FALSE]
  } else if (m == "m") {
    dx <- ia_inv[ov_idx, idx, drop = FALSE]
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  dx
}

# derivative of ML/GLS objective function F wrt the free parameters
lav_ram_df <- function(mlist = NULL, omega = NULL, omega_mu = NULL) {
  ov_idx <- as.integer(mlist$ov.idx[1, ])
  a <- mlist$A
  s <- mlist$S

  # nvar <- length(ov_idx)
  # nboth <- nrow(a)

  # get (I-A)^{-1}
  ia_inv <- lav_matrix_inverse_iminus(a)

  # meanstructure?
  meanstructure <- FALSE
  if (!is.null(omega_mu)) meanstructure <- TRUE

  # pre-compute
  t_ia_inv <- t(ia_inv)
  omega__ia_inv__s__t_ia_inv <- (omega %*% ia_inv[ov_idx, , drop = FALSE] %*%
                                                             s %*% t(ia_inv))

  # 1. A
  if (meanstructure) {
    a_deriv <-
      -1.0 * ((t(ia_inv)[, ov_idx, drop = FALSE] %*%
                                    (omega_mu %*% t(mlist$m)) %*% t(ia_inv)) +
        (t_ia_inv[, ov_idx, drop = FALSE] %*% omega__ia_inv__s__t_ia_inv))
  } else {
    a_deriv <- -1.0 * (t_ia_inv[, ov_idx, drop = FALSE] %*%
                                                    omega__ia_inv__s__t_ia_inv)
  }

  # 2. S
  s_deriv <- -1.0 * (t_ia_inv[, ov_idx, drop = FALSE] %*% omega %*%
                                                   ia_inv[ov_idx, , drop = FALSE])
  diag(s_deriv) <- 0.5 * diag(s_deriv)

  if (meanstructure) {
    m_deriv <- -1.0 * t(t(omega_mu) %*% ia_inv[ov_idx, , drop = FALSE])
  } else {
    m_deriv <- NULL
  }

  list(
    A = a_deriv,
    S = s_deriv,
    m = m_deriv
  )
}


# Jacobian of the model-implied moments wrt the free parameters
# (single block, RAM representation, classical continuous case only:
# no categorical, no correlation structure, no conditional.x).
#
# Inputs:
#   MLIST         : named list of model matrices for the block
#                   (ov.idx, A, S, [m])
#   m.free.idx    : list of length(MLIST); m.free.idx[[mm]] holds the
#                   matrix-element indices of the free elements of MLIST[[mm]]
#   x.free.idx    : list of length(MLIST); x.free.idx[[mm]] holds the
#                   parameter-vector positions of those free elements
#   nx.free       : total number of free parameters (column dim of out)
#   meanstructure : logical
#   group.w.free  : logical; if TRUE, prepend a row for the free group weight
#
# Output: a matrix with nx.free columns; rows are
# [group weight (if group.w.free) | Mu (if meanstructure)] then [vech(Sigma)].
lav_ram_dimplied_dx <- function(mlist         = NULL,
                                m_free_idx    = NULL,
                                x_free_idx    = NULL,
                                nx_free       = NULL,
                                meanstructure = FALSE,
                                group_w_free  = FALSE) {

  ov_idx <- as.integer(mlist$ov.idx[1L, ])
  nvar   <- length(ov_idx)
  nboth  <- nrow(mlist$A)
  pstar  <- nvar * (nvar + 1L) / 2L

  mnames <- names(mlist)

  # A (asymmetric)
  mm_a_idx <- which(mnames == "A")
  x_a_idx  <- x_free_idx[[mm_a_idx]]
  m_a_idx  <- m_free_idx[[mm_a_idx]]
  n_a      <- length(m_a_idx)

  # S (symmetric)
  mm_s_idx <- which(mnames == "S")
  tmp      <- x_free_idx[[mm_s_idx]]
  x_s_idx  <- tmp[!duplicated(tmp)]
  m_s_idx  <- m_free_idx[[mm_s_idx]][!duplicated(tmp)]
  n_s      <- length(m_s_idx)

  # m (only if meanstructure)
  n_m     <- 0L
  x_m_idx <- integer(0L)
  m_m_idx <- integer(0L)
  if (meanstructure) {
    mm_m_idx <- which(mnames == "m")
    x_m_idx  <- x_free_idx[[mm_m_idx]]
    m_m_idx  <- m_free_idx[[mm_m_idx]]
    n_m      <- length(m_m_idx)
  }

  # precompute
  ia_inv <- lav_matrix_inverse_iminus(mlist$A)            # nboth x nboth
  fobs   <- ia_inv[ov_idx, , drop = FALSE]                # nvar  x nboth
  m_w      <- ia_inv %*% mlist$S %*% t(fobs)                # nboth x nvar
  wt     <- t(m_w)                                          # nvar  x nboth
  if (meanstructure) {
    a <- as.vector(ia_inv %*% mlist$m)                    # nboth
  }

  # vech structure for Sigma
  r_s <- lav_matrix_vech_row_idx(nvar)
  c_s <- lav_matrix_vech_col_idx(nvar)

  # helper: vec index -> (row, col)
  vec2rc <- function(idx, nr) {
    cbind(row = (idx - 1L) %% nr + 1L,
          col = (idx - 1L) %/% nr + 1L)
  }

  # ---- jac_sigma ----
  n_free    <- n_a + n_s
  jac_sigma <- matrix(0, pstar, n_free)
  col       <- 1L

  # A[k,l]: dSigma[r,s] = Fobs[r,k]*W[l,s] + Fobs[s,k]*W[l,r]
  if (n_a > 0L) {
    rc  <- vec2rc(m_a_idx, nboth)
    k_v <- rc[, 1L]
    l_v <- rc[, 2L]

    t1 <- fobs[r_s, k_v, drop = FALSE] * wt[c_s, l_v, drop = FALSE]
    t2 <- fobs[c_s, k_v, drop = FALSE] * wt[r_s, l_v, drop = FALSE]

    jac_sigma[, col:(col + n_a - 1L)] <- t1 + t2
    col <- col + n_a
  }

  # S[k,l] symmetric: dSigma[r,s] = Fobs[r,k]*Fobs[s,l] + Fobs[r,l]*Fobs[s,k]
  # diagonal (k==l): halve (formula double-counts the single parameter)
  if (n_s > 0L) {
    rc  <- vec2rc(m_s_idx, nboth)
    k_v <- pmax(rc[, 1L], rc[, 2L])
    l_v <- pmin(rc[, 1L], rc[, 2L])

    t1 <- fobs[r_s, k_v, drop = FALSE] * fobs[c_s, l_v, drop = FALSE]
    t2 <- fobs[r_s, l_v, drop = FALSE] * fobs[c_s, k_v, drop = FALSE]
    dx <- t1 + t2

    diag_mask <- (k_v == l_v)
    if (any(diag_mask)) {
      dx[, diag_mask] <- dx[, diag_mask] * 0.5
    }

    jac_sigma[, col:(col + n_s - 1L)] <- dx
  }

  # ---- jac_mean (if meanstructure) ----
  if (meanstructure) {
    n_free_mu <- n_a + n_m
    jac_mean  <- matrix(0, nvar, n_free_mu)
    col       <- 1L

    # A[k,l]:  dMu[r] = Fobs[r,k] * a[l]
    if (n_a > 0L) {
      rc  <- vec2rc(m_a_idx, nboth)
      k_v <- rc[, 1L]
      l_v <- rc[, 2L]
      jac_mean[, col:(col + n_a - 1L)] <-
        fobs[, k_v, drop = FALSE] * rep(a[l_v], each = nvar)
      col <- col + n_a
    }

    # m[k]:  dMu[r] = Fobs[r, k]
    if (n_m > 0L) {
      jac_mean[, col:(col + n_m - 1L)] <- fobs[, m_m_idx, drop = FALSE]
    }
  }

  # ---- assemble output ----
  out    <- matrix(0, nrow = pstar, ncol = nx_free)
  el_idx <- c(x_a_idx, x_s_idx)
  out[, el_idx] <- jac_sigma

  if (meanstructure) {
    el_idx_mu      <- c(x_a_idx, x_m_idx)
    outm           <- matrix(0, nrow = nvar, ncol = nx_free)
    outm[, el_idx_mu] <- jac_mean
    out <- rbind(outm, out)
  }

  # group weight: prepend a row with 1.0 at the gw column
  if (group_w_free) {
    mm_gw_idx <- which(mnames == "gw")
    if (length(mm_gw_idx) > 0L) {
      x_gw_idx <- x_free_idx[[mm_gw_idx]]
      out_gw <- matrix(0, nrow = 1L, ncol = nx_free)
      out_gw[1L, x_gw_idx] <- 1.0
      out <- rbind(out_gw, out)
    }
  }

  out
}
