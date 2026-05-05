# RAM representation
#
# initial version: YR 2021-10-04

lav_ram <- function(partable = NULL,
                    target = NULL,
                    extra = FALSE,
                    remove.nonexisting = TRUE) {
  # prepare target list
  if (is.null(target)) target <- partable

  stopifnot(!is.null(target$block))

  # not for categorical data (yet)
  if (any(partable$op == "|")) {
    lav_msg_stop(gettext("RAM representation is not (yet) supported for
                         categorical endogenous variables."))
  }

  # not for conditional.x = TRUE yet
  conditional.x <- any(partable$exo > 0L & partable$op == "~")
  if (conditional.x) {
    lav_msg_stop(gettext("RAM representation is not (yet) supported
                         if conditional.x = TRUE"))
  }

  # prepare output
  N <- length(target$lhs)
  tmp.mat <- character(N)
  tmp.row <- integer(N)
  tmp.col <- integer(N)

  # global settings
  meanstructure <- any(partable$op == "~1")
  categorical <- any(partable$op == "|")
  group.w.free <- any(partable$lhs == "group" & partable$op == "%")

  # number of blocks
  nblocks <- lav_partable_nblocks(partable)

  # always return ov.idx
  ov.idx <- vector("list", nblocks)
  ov.dummy.names.nox <- vector("list", nblocks)
  ov.dummy.names.x <- vector("list", nblocks)
  if (extra) {
    REP.mmNames <- vector("list", nblocks)
    REP.mmNumber <- vector("list", nblocks)
    REP.mmRows <- vector("list", nblocks)
    REP.mmCols <- vector("list", nblocks)
    REP.mmDimNames <- vector("list", nblocks)
    REP.mmSymmetric <- vector("list", nblocks)
  }

  for (g in 1:nblocks) {
    # info from user model per block
    ov.names <- lav_partable_vnames(partable, "ov", block = g)
    nvar <- length(ov.names)
    ov.idx[[g]] <- seq_len(nvar)
    ov.dummy.names.nox[[g]] <- character(0)
    ov.dummy.names.x[[g]] <- character(0)

    lv.names <- lav_partable_vnames(partable, "lv", block = g)
    both.names <- c(ov.names, lv.names)
    nboth <- length(both.names)

    # 1. "=~" indicators
    idx <- which(target$block == g & target$op == "=~")
    tmp.mat[idx] <- "A"
    tmp.row[idx] <- match(target$rhs[idx], both.names)
    tmp.col[idx] <- match(target$lhs[idx], both.names)

    # 2. "~" regressions
    idx <- which(target$block == g & (target$op == "~" |
      target$op == "<~"))
    tmp.mat[idx] <- "A"
    tmp.row[idx] <- match(target$lhs[idx], both.names)
    tmp.col[idx] <- match(target$rhs[idx], both.names)

    # 3. "~~" variances/covariances
    idx <- which(target$block == g & target$op == "~~")
    tmp.mat[idx] <- "S"
    tmp.row[idx] <- match(target$lhs[idx], both.names)
    tmp.col[idx] <- match(target$rhs[idx], both.names)

    # catch lower-elements in theta/psi
    idx.lower <- which(tmp.mat == "S" & tmp.row > tmp.col)
    if (length(idx.lower) > 0L) {
      tmp <- tmp.row[idx.lower]
      tmp.row[idx.lower] <- tmp.col[idx.lower]
      tmp.col[idx.lower] <- tmp
    }

    # 4. "~1" means/intercepts
    idx <- which(target$block == g & target$op == "~1")
    tmp.mat[idx] <- "m"
    tmp.row[idx] <- match(target$lhs[idx], both.names)
    tmp.col[idx] <- 1L

    # 5. "|" th
    # not used yet

    # 6. "~*~" scales
    # not used yet

    # 7. group weights
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "%")
    tmp.mat[idx] <- "gw"
    tmp.row[idx] <- 1L
    tmp.col[idx] <- 1L

    # 8. instruments
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "|~")
    tmp.mat[idx] <- "miiv"
    tmp.row[idx] <- 0L
    tmp.col[idx] <- 0L

    if (extra) {
      # mRows
      mmRows <- list(
        ov.idx = 1L,
        A = nboth,
        S = nboth,
        m = nboth,
        gw = 1L
      )

      # mCols
      mmCols <- list(
        ov.idx = nvar,
        A = nboth,
        S = nboth,
        m = 1L,
        gw = 1L
      )

      # dimNames for LISREL model matrices
      mmDimNames <- list(
        ov.idx = list("ov.idx", ov.names),
        A = list(both.names, both.names),
        S = list(both.names, both.names),
        m = list(both.names, "intercept"),
        gw = list("group", "weight")
      )
      # isSymmetric
      mmSymmetric <- list(
        ov.idx = FALSE,
        A = FALSE,
        S = TRUE,
        m = FALSE,
        gw = FALSE
      )

      # which mm's do we need? (always include ov.idx, A and S)
      IDX <- which(target$block == g)
      mmNames <- c("ov.idx", "A", "S")
      if (meanstructure) {
        mmNames <- c(mmNames, "m")
      }
      if ("gw" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "gw")
      }

      REP.mmNames[[g]] <- mmNames
      REP.mmNumber[[g]] <- length(mmNames)
      REP.mmRows[[g]] <- unlist(mmRows[mmNames])
      REP.mmCols[[g]] <- unlist(mmCols[mmNames])
      REP.mmDimNames[[g]] <- mmDimNames[mmNames]
      REP.mmSymmetric[[g]] <- unlist(mmSymmetric[mmNames])
    } # extra
  } # nblocks

  REP <- list(
    mat = tmp.mat,
    row = tmp.row,
    col = tmp.col
  )

  # always return ov.idx attribute
  attr(REP, "ov.idx") <- ov.idx
  attr(REP, "ov.dummy.names.nox") <- ov.dummy.names.nox
  attr(REP, "ov.dummy.names.x") <- ov.dummy.names.x

  if (extra) {
    attr(REP, "mmNames") <- REP.mmNames
    attr(REP, "mmNumber") <- REP.mmNumber
    attr(REP, "mmRows") <- REP.mmRows
    attr(REP, "mmCols") <- REP.mmCols
    attr(REP, "mmDimNames") <- REP.mmDimNames
    attr(REP, "mmSymmetric") <- REP.mmSymmetric
  }

  REP
}

# the model-implied variance/covariance matrix of the observed variables
lav_ram_sigmahat <- function(MLIST = NULL, delta = NULL) {
  ov.idx <- as.integer(MLIST$ov.idx[1, ])
  A <- MLIST$A
  S <- MLIST$S

  # get (I-A)^{-1}
  IA.inv <- lav_matrix_inverse_iminus(A)

  # compute Sigma for all ov and lv
  VYeta <- tcrossprod(IA.inv %*% S, IA.inv)

  # select only observed part
  VY <- VYeta[ov.idx, ov.idx, drop = FALSE]

  # if delta, scale
  if (!is.null(MLIST$delta) && delta) {
    nvar <- ncol(VY)
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VY <- DELTA %*% VY %*% DELTA
  }

  VY
}

# VETA: the variance/covariance matrix of the latent variables only
lav_ram_veta <- function(MLIST = NULL) {
  ov.idx <- as.integer(MLIST$ov.idx[1, ])
  A <- MLIST$A
  S <- MLIST$S

  # get (I-A)^{-1}
  IA.inv <- lav_matrix_inverse_iminus(A)

  # compute Sigma for all ov and lv
  VYeta <- tcrossprod(IA.inv %*% S, IA.inv)

  # select only latent part
  VETA <- VYeta[-ov.idx, -ov.idx, drop = FALSE]

  VETA
}

# MuHat: the model-implied means/intercepts
lav_ram_muhat <- function(MLIST = NULL) {
  ov.idx <- as.integer(MLIST$ov.idx[1, ])
  A <- MLIST$A
  m <- MLIST$m

  # shortcut
  if (is.null(m)) {
    return(matrix(0, nrow = length(ov.idx), 1L))
  }

  # get (I-A)^{-1}
  IA.inv <- lav_matrix_inverse_iminus(A)

  # all means/intercepts
  EYeta <- IA.inv %*% m

  # select observed only
  muhat <- EYeta[ov.idx, , drop = FALSE]

  muhat
}

# derivative of 'Sigma' wrt the (freel) elements in A and/or S
lav_ram_dsigma <- function(m = "A",
                           idx = seq_len(length(MLIST[[m]])),
                           MLIST = NULL,
                           vech = TRUE) {
  ov.idx <- as.integer(MLIST$ov.idx[1, ])
  A <- MLIST$A
  S <- MLIST$S

  nvar <- length(ov.idx)
  nboth <- nrow(A)

  # shortcut for ov.idx, m, ...
  if (!m %in% c("A", "S")) {
    pstar <- nvar * (nvar + 1) / 2
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # get (I-A)^{-1}
  IA.inv <- lav_matrix_inverse_iminus(A)

  if (m == "A") {
    L1 <- (IA.inv %*% S %*% t(IA.inv))[ov.idx, , drop = FALSE]
    KOL.idx <- matrix(1:(nboth * nboth), nboth, nboth, byrow = TRUE)[idx]
    DX <- (L1 %x% IA.inv[ov.idx, , drop = FALSE])[, idx, drop = FALSE] +
      (IA.inv[ov.idx, , drop = FALSE] %x% L1)[, KOL.idx, drop = FALSE]
    # this is not really needed (because we select idx=m.el.idx)
    # but just in case we need all elements of beta...
    DX[, which(idx %in% lav_matrix_diag_idx(nboth))] <- 0.0
  } else if (m == "S") {
    DX <- (IA.inv[ov.idx, , drop = FALSE] %x% IA.inv[ov.idx, , drop = FALSE])
    # symmetry correction, but keeping all duplicated elements
    # since we depend on idx=m.el.idx
    lower.idx <- lav_matrix_vech_idx(nboth, diagonal = FALSE)
    upper.idx <- lav_matrix_vechru_idx(nboth, diagonal = FALSE)
    offdiagSum <- DX[, lower.idx] + DX[, upper.idx]
    DX[, c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
    DX <- DX[, idx, drop = FALSE]
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  # vech?
  if (vech) {
    v.idx <- lav_matrix_vech_idx(nvar)
    DX <- DX[v.idx, , drop = FALSE]
  }

  DX
}

# derivative of 'Mu' wrt the (free) elements in A and/or m
lav_ram_dmu <- function(m = "A",
                        idx = seq_len(length(MLIST[[m]])),
                        MLIST = NULL,
                        vech = TRUE) {
  ov.idx <- as.integer(MLIST$ov.idx[1, ])
  A <- MLIST$A
  S <- MLIST$S

  nvar <- length(ov.idx)
  nboth <- nrow(A)

  # shortcut for ov.idx, m, ...
  if (!m %in% c("A", "m")) {
    return(matrix(0.0, nrow = nvar, ncol = length(idx)))
  }

  # get (I-A)^{-1}
  IA.inv <- lav_matrix_inverse_iminus(A)

  if (m == "A") {
    DX <- (t(IA.inv %*% MLIST$m) %x% IA.inv)[ov.idx, idx, drop = FALSE]
  } else if (m == "m") {
    DX <- IA.inv[ov.idx, idx, drop = FALSE]
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  DX
}

# derivative of ML/GLS objective function F wrt the free parameters
lav_ram_df <- function(MLIST = NULL, Omega = NULL, Omega.mu = NULL) {
  ov.idx <- as.integer(MLIST$ov.idx[1, ])
  A <- MLIST$A
  S <- MLIST$S

  nvar <- length(ov.idx)
  nboth <- nrow(A)

  # get (I-A)^{-1}
  IA.inv <- lav_matrix_inverse_iminus(A)

  # meanstructure?
  meanstructure <- FALSE
  if (!is.null(Omega.mu)) meanstructure <- TRUE

  # pre-compute
  tIA.inv <- t(IA.inv)
  Omega..IA.inv..S..tIA.inv <- (Omega %*% IA.inv[ov.idx, , drop = FALSE] %*% S %*% t(IA.inv))

  # 1. A
  if (meanstructure) {
    A.deriv <-
      -1.0 * ((t(IA.inv)[, ov.idx, drop = FALSE] %*% (Omega.mu %*% t(MLIST$m)) %*% t(IA.inv)) +
        (tIA.inv[, ov.idx, drop = FALSE] %*% Omega..IA.inv..S..tIA.inv))
  } else {
    A.deriv <- -1.0 * (tIA.inv[, ov.idx, drop = FALSE] %*% Omega..IA.inv..S..tIA.inv)
  }

  # 2. S
  S.deriv <- -1.0 * (tIA.inv[, ov.idx, drop = FALSE] %*% Omega %*% IA.inv[ov.idx, , drop = FALSE])
  diag(S.deriv) <- 0.5 * diag(S.deriv)

  if (meanstructure) {
    m.deriv <- -1.0 * t(t(Omega.mu) %*% IA.inv[ov.idx, , drop = FALSE])
  } else {
    m.deriv <- NULL
  }

  list(
    A = A.deriv,
    S = S.deriv,
    m = m.deriv
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
lav_ram_dimplied_dx <- function(MLIST         = NULL,
                                m.free.idx    = NULL,
                                x.free.idx    = NULL,
                                nx.free       = NULL,
                                meanstructure = FALSE,
                                group.w.free  = FALSE) {

  ov.idx <- as.integer(MLIST$ov.idx[1L, ])
  nvar   <- length(ov.idx)
  nboth  <- nrow(MLIST$A)
  pstar  <- nvar * (nvar + 1L) / 2L

  mnames <- names(MLIST)

  # A (asymmetric)
  mm.A.idx <- which(mnames == "A")
  x.A.idx  <- x.free.idx[[mm.A.idx]]
  m.A.idx  <- m.free.idx[[mm.A.idx]]
  n_A      <- length(m.A.idx)

  # S (symmetric)
  mm.S.idx <- which(mnames == "S")
  tmp      <- x.free.idx[[mm.S.idx]]
  x.S.idx  <- tmp[!duplicated(tmp)]
  m.S.idx  <- m.free.idx[[mm.S.idx]][!duplicated(tmp)]
  n_S      <- length(m.S.idx)

  # m (only if meanstructure)
  n_m     <- 0L
  x.m.idx <- integer(0L)
  m.m.idx <- integer(0L)
  if (meanstructure) {
    mm.m.idx <- which(mnames == "m")
    x.m.idx  <- x.free.idx[[mm.m.idx]]
    m.m.idx  <- m.free.idx[[mm.m.idx]]
    n_m      <- length(m.m.idx)
  }

  # precompute
  IA.inv <- lav_matrix_inverse_iminus(MLIST$A)            # nboth x nboth
  Fobs   <- IA.inv[ov.idx, , drop = FALSE]                # nvar  x nboth
  W      <- IA.inv %*% MLIST$S %*% t(Fobs)                # nboth x nvar
  WT     <- t(W)                                          # nvar  x nboth
  if (meanstructure) {
    a <- as.vector(IA.inv %*% MLIST$m)                    # nboth
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
  n_free    <- n_A + n_S
  jac_sigma <- matrix(0, pstar, n_free)
  col       <- 1L

  # A[k,l]: dSigma[r,s] = Fobs[r,k]*W[l,s] + Fobs[s,k]*W[l,r]
  if (n_A > 0L) {
    rc  <- vec2rc(m.A.idx, nboth);  k_v <- rc[, 1L];  l_v <- rc[, 2L]

    T1 <- Fobs[r_s, k_v, drop = FALSE] * WT[c_s, l_v, drop = FALSE]
    T2 <- Fobs[c_s, k_v, drop = FALSE] * WT[r_s, l_v, drop = FALSE]

    jac_sigma[, col:(col + n_A - 1L)] <- T1 + T2
    col <- col + n_A
  }

  # S[k,l] symmetric: dSigma[r,s] = Fobs[r,k]*Fobs[s,l] + Fobs[r,l]*Fobs[s,k]
  # diagonal (k==l): halve (formula double-counts the single parameter)
  if (n_S > 0L) {
    rc  <- vec2rc(m.S.idx, nboth)
    k_v <- pmax(rc[, 1L], rc[, 2L])
    l_v <- pmin(rc[, 1L], rc[, 2L])

    T1 <- Fobs[r_s, k_v, drop = FALSE] * Fobs[c_s, l_v, drop = FALSE]
    T2 <- Fobs[r_s, l_v, drop = FALSE] * Fobs[c_s, k_v, drop = FALSE]
    DX <- T1 + T2

    diag_mask <- (k_v == l_v)
    if (any(diag_mask)) {
      DX[, diag_mask] <- DX[, diag_mask] * 0.5
    }

    jac_sigma[, col:(col + n_S - 1L)] <- DX
  }

  # ---- jac_mean (if meanstructure) ----
  if (meanstructure) {
    n_free_mu <- n_A + n_m
    jac_mean  <- matrix(0, nvar, n_free_mu)
    col       <- 1L

    # A[k,l]:  dMu[r] = Fobs[r,k] * a[l]
    if (n_A > 0L) {
      rc  <- vec2rc(m.A.idx, nboth);  k_v <- rc[, 1L];  l_v <- rc[, 2L]
      jac_mean[, col:(col + n_A - 1L)] <-
        Fobs[, k_v, drop = FALSE] * rep(a[l_v], each = nvar)
      col <- col + n_A
    }

    # m[k]:  dMu[r] = Fobs[r, k]
    if (n_m > 0L) {
      jac_mean[, col:(col + n_m - 1L)] <- Fobs[, m.m.idx, drop = FALSE]
    }
  }

  # ---- assemble output ----
  out    <- matrix(0, nrow = pstar, ncol = nx.free)
  el.idx <- c(x.A.idx, x.S.idx)
  out[, el.idx] <- jac_sigma

  if (meanstructure) {
    el.idx_mu      <- c(x.A.idx, x.m.idx)
    outm           <- matrix(0, nrow = nvar, ncol = nx.free)
    outm[, el.idx_mu] <- jac_mean
    out <- rbind(outm, out)
  }

  # group weight: prepend a row with 1.0 at the gw column
  if (group.w.free) {
    mm.gw.idx <- which(mnames == "gw")
    if (length(mm.gw.idx) > 0L) {
      x.gw.idx <- x.free.idx[[mm.gw.idx]]
      out_gw <- matrix(0, nrow = 1L, ncol = nx.free)
      out_gw[1L, x.gw.idx] <- 1.0
      out <- rbind(out_gw, out)
    }
  }

  out
}
