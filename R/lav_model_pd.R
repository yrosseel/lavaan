# pd_lv = TRUE: keep the variance-covariance matrix of the latent variables
# positive definite during optimization
#
# YR 24 Nov 2026: first version
#
# The (residual) variances and covariances of the latent variables live in
# the PSI matrix. Because VETA = (I-B)^{-1} PSI (I-B)^{-T}, the model-implied
# variance-covariance matrix of the latent variables is positive definite
# if (and only if) PSI is positive definite.
#
# To keep (a block of) PSI positive definite, the optimizer does not see the
# free elements of PSI, but the elements of its Cholesky factor L instead
# (PSI = L L'), with the diagonal elements of L on a log scale:
#
#   - free variance psi[i,i]            <->   log(L[i,i])
#   - free covariance psi[i,j] (i > j)  <->   L[i,j]
#   - psi[i,j] fixed to zero (i > j)    ->    L[i,j] is NOT a parameter: it
#       is computed (row by row) from the elements of L that precede it:
#           L[i,j] = -( sum_{k<j} L[i,k] * L[j,k] ) / L[j,j]
#       so that psi[i,j] is exactly zero
#
# This is a smooth one-to-one mapping between the free elements of the PSI
# block and the set of all positive definite matrices with the same pattern
# of fixed zeros (whatever the ordering of the variables). The number of
# parameters does not change, and the transformation is only used during
# optimization (see lav_model_est()): everything else (the objective
# function, the gradient, the information matrix, ...) keeps using the
# original parameterization.
#
# PSI is split in 'blocks': the connected components of the pattern of its
# (free or fixed) nonzero off-diagonal elements. Only the blocks that
# contain at least one (regular) latent variable and at least one free
# parameter are reparameterized; these blocks must be 'unconstrained':
#   - all variances must be free
#   - the covariances are either free, or fixed to zero
#   - the free elements are not involved in (in)equality constraints
#   - the free elements have no (user-specified) bounds
# All other blocks (eg variances fixed to zero, the fixed.x block, blocks
# with observed variables only) are left alone.
#
# Note: the automatic bounds (bounds = "standard", "pos.var", ...) can not
# be imposed on the elements of the Cholesky factor: they are ignored for
# the parameters that are reparameterized.

# setup: find the PSI blocks that must be reparameterized, and check if we
# can handle them (if not, an error is raised)
#
# returns NULL if there is nothing to do, or a list with elements:
#   - blocks: a list with one element per PSI block
#   - x_idx:  the positions (in the vector of free parameters) of all the
#             parameters that are reparameterized
lav_model_pd_setup <- function(lavmodel = NULL, lavpartable = NULL,
                               lavoptions = NULL) {
  if (lavmodel@representation != "LISREL") {
    lav_msg_stop(gettext(
      "pd_lv = TRUE is only available for the LISREL representation."))
  }

  # how many times is each free parameter used (over all model matrices)?
  # (more than once: simple equality constraints)
  glist <- lavmodel@GLIST
  x_all <- integer(0L)
  for (mm in seq_along(glist)) {
    if (length(lavmodel@m.free.idx[[mm]]) == 0L) {
      next
    }
    if (lavmodel@isSymmetric[mm]) {
      id <- matrix(0L, nrow(glist[[mm]]), ncol(glist[[mm]]))
      id[lavmodel@m.free.idx[[mm]]] <- lavmodel@x.free.idx[[mm]]
      id <- id[lower.tri(id, diag = TRUE)]
      x_all <- c(x_all, id[id > 0L])
    } else {
      x_all <- c(x_all, lavmodel@x.free.idx[[mm]])
    }
  }
  x_dup <- unique(x_all[duplicated(x_all)])

  # free parameters that are involved in explicit (in)equality constraints
  x_con <- integer(0L)
  con_idx <- which(lavpartable$op %in% c("==", "<", ">"))
  if (length(con_idx) > 0L) {
    con_txt <- c(lavpartable$lhs[con_idx], lavpartable$rhs[con_idx])
    con_ref <- unique(unlist(lapply(con_txt, function(txt) {
      out <- try(all.vars(parse(text = txt)), silent = TRUE)
      if (inherits(out, "try-error")) txt else out
    })))
    lab_flag <- (nzchar(lavpartable$label) &
                 lavpartable$label %in% con_ref)
    if (!is.null(lavpartable$plabel)) {
      lab_flag <- lab_flag | (nzchar(lavpartable$plabel) &
                              lavpartable$plabel %in% con_ref)
    }
    x_con <- unique(lavpartable$free[lab_flag & lavpartable$free > 0L])
  }

  # free parameters with a user-specified bound (eg "a > 0.1"); a
  # non-positive lower bound for a variance is harmless; if automatic
  # bounds are used (the bounds= option), we can not tell them apart
  x_bnd <- integer(0L)
  if (is.null(lavoptions$bounds) || lavoptions$bounds == "none") {
    npt <- length(lavpartable$lhs)
    lower <- lavpartable$lower
    if (is.null(lower)) {
      lower <- rep(-Inf, npt)
    }
    upper <- lavpartable$upper
    if (is.null(upper)) {
      upper <- rep(+Inf, npt)
    }
    var_flag <- lavpartable$op == "~~" & lavpartable$lhs == lavpartable$rhs
    bnd_flag <- is.finite(upper) | (is.finite(lower) & !(var_flag & lower <= 0))
    x_bnd <- unique(lavpartable$free[bnd_flag & lavpartable$free > 0L])
  }

  blocks <- list()
  psi_idx <- which(names(glist) == "psi")
  for (b in seq_along(psi_idx)) {
    mm <- psi_idx[b]
    psi <- glist[[mm]]
    n <- nrow(psi)
    if (n == 0L || length(lavmodel@m.free.idx[[mm]]) == 0L) {
      next
    }
    psi_names <- lavmodel@dimNames[[mm]][[1L]]
    lv_names <- lav_pt_vnames(lavpartable, type = "lv.regular", block = b)

    # free parameter id's (0 = fixed)
    id <- matrix(0L, n, n)
    id[lavmodel@m.free.idx[[mm]]] <- lavmodel@x.free.idx[[mm]]

    # connected components of the nonzero (free or fixed) pattern
    adj <- (id > 0L) | (psi != 0)
    adj <- adj | t(adj)
    diag(adj) <- TRUE
    comp <- lav_model_pd_components(adj)

    for (cc in unique(comp)) {
      v <- which(comp == cc)
      id_c <- id[v, v, drop = FALSE]
      # nothing free, or no latent variables: leave this block alone
      if (!any(id_c > 0L) || !any(psi_names[v] %in% lv_names)) {
        next
      }
      psi_c <- psi[v, v, drop = FALSE]
      v_names <- psi_names[v]
      where_txt <- ""
      if (length(psi_idx) > 1L) {
        where_txt <- gettextf(" (block %s)", b)
      }

      # 1. all variances must be free
      bad <- which(diag(id_c) == 0L)
      if (length(bad) > 0L) {
        lav_msg_stop(gettextf(
          "pd_lv = TRUE is (currently) only available if the variances of
          the (correlated) latent variables are free parameters; the
          following variances are fixed%1$s: %2$s. (Perhaps std_lv = TRUE
          was used?)", where_txt,
          lav_msg_view(v_names[bad], "none", FALSE)))
      }

      # 2. covariances: free, or fixed to zero
      low <- lower.tri(id_c)
      bad <- which(low & id_c == 0L & psi_c != 0, arr.ind = TRUE)
      if (nrow(bad) > 0L) {
        lav_msg_stop(gettextf(
          "pd_lv = TRUE is (currently) only available if the covariances
          among the latent variables are either free, or fixed to zero;
          the following covariances are fixed to a nonzero value%1$s: %2$s.",
          where_txt, lav_msg_view(paste0(v_names[bad[, 2L]], " ~~ ",
                                         v_names[bad[, 1L]]),
                                  "none", FALSE)))
      }

      # 3. no (in)equality constraints
      low_d <- lower.tri(id_c, diag = TRUE)
      bad <- which(low_d & (id_c %in% x_dup | id_c %in% x_con),
                   arr.ind = TRUE)
      if (nrow(bad) > 0L) {
        lav_msg_stop(gettextf(
          "pd_lv = TRUE is (currently) not available if the variances or
          covariances of the latent variables are involved in equality or
          inequality constraints; the following parameters are
          constrained%1$s: %2$s.",
          where_txt, lav_msg_view(paste0(v_names[bad[, 2L]], " ~~ ",
                                         v_names[bad[, 1L]]),
                                  "none", FALSE)))
      }

      # 4. no (user-specified) bounds
      bad <- which(low_d & id_c %in% x_bnd, arr.ind = TRUE)
      if (nrow(bad) > 0L) {
        lav_msg_stop(gettextf(
          "pd_lv = TRUE is (currently) not available if the variances or
          covariances of the latent variables have (user-specified) lower
          or upper bounds; the following parameters have bounds%1$s: %2$s.",
          where_txt, lav_msg_view(paste0(v_names[bad[, 2L]], " ~~ ",
                                         v_names[bad[, 1L]]),
                                  "none", FALSE)))
      }

      # lower-triangular elements: free or (fixed to) zero
      free_pos <- which(low & id_c > 0L, arr.ind = TRUE)
      zero_cols <- lapply(seq_along(v), function(i) {
        which(id_c[i, seq_len(i - 1L)] == 0L)
      })
      blocks[[length(blocks) + 1L]] <- list(
        n = length(v),
        block = b,
        names = v_names,
        diag_x = diag(id_c),
        low_pos = free_pos,
        low_x = id_c[free_pos],
        zero_cols = zero_cols,
        zero_flag = any(lengths(zero_cols) > 0L)
      )
    } # components
  } # psi matrices

  if (length(blocks) == 0L) {
    return(NULL)
  }

  list(
    blocks = blocks,
    x_idx = unlist(lapply(blocks, function(blk) c(blk$diag_x, blk$low_x)))
  )
}

# connected components of a (symmetric, logical) adjacency matrix
# returns an integer vector with the component number of each row
lav_model_pd_components <- function(adj) {
  n <- nrow(adj)
  comp <- integer(n)
  cc <- 0L
  for (i in seq_len(n)) {
    if (comp[i] > 0L) {
      next
    }
    cc <- cc + 1L
    todo <- i
    while (length(todo) > 0L) {
      comp[todo] <- cc
      nb <- which(colSums(adj[todo, , drop = FALSE]) > 0L)
      todo <- nb[comp[nb] == 0L]
    }
  }
  comp
}

# the Cholesky factor L of a single block, given the transformed
# parameters (x_pd: the full parameter vector in the 'Cholesky' metric)
lav_model_pd_block_l <- function(x_pd, blk) {
  n <- blk$n
  l_mat <- matrix(0, n, n)
  l_mat[blk$low_pos] <- x_pd[blk$low_x]
  diag(l_mat) <- exp(x_pd[blk$diag_x])

  # elements that are implied by the fixed zeros in PSI
  if (blk$zero_flag) {
    for (i in seq_len(n)[-1L]) {
      for (j in blk$zero_cols[[i]]) {
        k <- seq_len(j - 1L)
        l_mat[i, j] <- -sum(l_mat[i, k] * l_mat[j, k]) / l_mat[j, j]
      }
    }
  }

  l_mat
}

# 'Cholesky' metric -> original (model) metric
lav_model_pd_to_model <- function(x_pd, pd) {
  x <- x_pd
  for (blk in pd$blocks) {
    l_mat <- lav_model_pd_block_l(x_pd, blk)
    psi <- tcrossprod(l_mat)
    x[blk$diag_x] <- diag(psi)
    x[blk$low_x] <- psi[blk$low_pos]
  }
  x
}

# the PSI block, given the parameters in the original (model) metric
lav_model_pd_block_psi <- function(x, blk) {
  psi <- matrix(0, blk$n, blk$n)
  psi[blk$low_pos] <- x[blk$low_x]
  psi <- psi + t(psi)
  diag(psi) <- x[blk$diag_x]
  psi
}

# original (model) metric -> 'Cholesky' metric
#
# this is only used for the starting values; if a PSI block is not
# (comfortably) positive definite, the starting values of the covariances
# are shrunken towards zero until it is
lav_model_pd_to_cholesky <- function(x, pd) {
  x_pd <- x
  for (blk in pd$blocks) {
    psi <- lav_model_pd_block_psi(x, blk)
    psi_var <- diag(psi)
    bad <- which(!is.finite(psi_var) | psi_var <= 0)
    if (length(bad) > 0L) {
      psi_var[bad] <- 0.05
    }
    psi[!is.finite(psi)] <- 0
    diag(psi) <- psi_var

    l_mat <- NULL
    for (attempt in 1:20) {
      l_try <- try(t(chol(psi)), silent = TRUE)
      # 'comfortably': conditional sd > 5% of the marginal sd
      if (!inherits(l_try, "try-error") &&
          all(diag(l_try) > 0.05 * sqrt(psi_var))) {
        l_mat <- l_try
        break
      }
      psi <- psi * 0.5
      diag(psi) <- psi_var
    }
    if (is.null(l_mat)) {
      l_mat <- diag(sqrt(psi_var), nrow = blk$n)
    }
    x_pd[blk$diag_x] <- log(diag(l_mat))
    x_pd[blk$low_x] <- l_mat[blk$low_pos]
  }
  x_pd
}

# chain rule: the gradient with respect to the free parameters in the
# original (model) metric (dx) -> the gradient in the 'Cholesky' metric
lav_model_pd_gradient <- function(dx, x_pd, pd) {
  dx_pd <- dx
  for (blk in pd$blocks) {
    n <- blk$n
    l_mat <- lav_model_pd_block_l(x_pd, blk)

    # dF/dPSI, with psi[i,j] and psi[j,i] as two distinct elements
    g_mat <- matrix(0, n, n)
    g_mat[blk$low_pos] <- dx[blk$low_x] / 2
    g_mat <- g_mat + t(g_mat)
    diag(g_mat) <- dx[blk$diag_x]

    # dF/dL, all (lower-triangular) elements of L being independent
    l_bar <- 2 * g_mat %*% l_mat

    # reverse sweep: the elements of L implied by the fixed zeros in PSI
    # are functions of the elements that precede them
    if (blk$zero_flag) {
      for (i in rev(seq_len(n)[-1L])) {
        for (j in rev(blk$zero_cols[[i]])) {
          a <- l_bar[i, j] / l_mat[j, j]
          k <- seq_len(j - 1L)
          l_bar[i, k] <- l_bar[i, k] - a * l_mat[j, k]
          l_bar[j, k] <- l_bar[j, k] - a * l_mat[i, k]
          l_bar[j, j] <- l_bar[j, j] - a * l_mat[i, j]
        }
      }
    }

    dx_pd[blk$low_x] <- l_bar[blk$low_pos]
    dx_pd[blk$diag_x] <- diag(l_bar) * diag(l_mat) # log scale
  }
  dx_pd
}

# parameter scaling (optim.parscale): the scaling factors for the
# parameters in the 'Cholesky' metric, given the scaling factors of the
# original parameters (variance: s_i^2, covariance: s_i * s_j)
lav_model_pd_parscale <- function(parscale, pd) {
  for (blk in pd$blocks) {
    s2 <- parscale[blk$diag_x]
    if (any(!is.finite(s2)) || any(s2 <= 0)) {
      s2 <- rep(1.0, blk$n)
    }
    # log(L[i,i]): the scaling is a shift, not a factor
    parscale[blk$diag_x] <- 1.0
    # L[i,j]: scaled like the sd of variable i
    parscale[blk$low_x] <- sqrt(s2)[blk$low_pos[, 1L]]
  }
  parscale
}

# is the solution (x: original metric) on the boundary of the positive
# definite region? The optimizer can only approach this boundary, so we
# use a tolerance: a PSI block is flagged if a variable is (almost)
# perfectly predicted by the other variables of the block (tol_cor), or if
# its variance has (almost) vanished (tol_var) relative to var_ref: a
# reference variance for each block (group/level) of the model (typically
# the largest model-implied variance of the observed variables)
#
# returns the names of the variables in the offending blocks
lav_model_pd_boundary <- function(x, pd, var_ref = NULL,
                                  tol_cor = 1e-06, tol_var = 1e-06) {
  out <- character(0L)
  for (blk in pd$blocks) {
    psi <- lav_model_pd_block_psi(x, blk)
    psi_var <- diag(psi)
    flag <- FALSE
    if (any(!is.finite(psi)) || any(psi_var <= 0)) {
      flag <- TRUE
    } else if (blk$n > 1L) {
      ev <- eigen(cov2cor(psi), symmetric = TRUE, only.values = TRUE)$values
      flag <- min(ev) < tol_cor
    }
    if (!flag && length(var_ref) >= blk$block &&
        isTRUE(var_ref[blk$block] > 0)) {
      flag <- any(psi_var < tol_var * var_ref[blk$block])
    }
    if (flag) {
      out <- c(out, blk$names)
    }
  }
  unique(out)
}
