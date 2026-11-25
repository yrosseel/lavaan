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
#   - variance psi[i,i] fixed to a (positive) value v (eg std_lv = TRUE):
#       row i of L must have length sqrt(v), and log(L[i,i]) is NOT a
#       parameter; instead, we use an unnormalized row w with w[i,i] = 1,
#           free covariance psi[i,j] (i > j)  <->   asinh(w[i,j])
#       (the elements implied by fixed zeros are computed as above), and
#       row i of L is sqrt(v) * w / ||w||, so that psi[i,i] is exactly v;
#       the boundary (a singular PSI) corresponds to w[i,j] -> +/- infinity,
#       and the asinh() scale makes sure that the optimizer can get close
#       (just like the log scale for the free variances)
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
# contain at least one latent variable (a regular latent variable, or a
# random slope) and at least one free parameter are reparameterized; these
# blocks must be 'unconstrained':
#   - the variances are either free, or fixed to a positive value
#   - the covariances are either free, or fixed to zero
#   - the free elements are not involved in (in)equality constraints
#   - the free elements have no (user-specified) bounds
# All other blocks (eg variances fixed to zero, the fixed.x block, blocks
# with observed variables only) are left alone.
#
# EM (two-level models): the M-step is itself a lavaan() fit of a
# two-group model to the expected complete-data moments; with pd_lv = TRUE,
# this M-step is also done with pd_lv = TRUE (a constrained M-step), so
# that all EM iterates are positive definite. In the M-step model of the
# random-slope EM, the random slopes are observed variables of the second
# group: their names are passed on via 'pd_lv_extra_names' in
# lavaan_cache_env (see lav_model_pd_extra_names()), and their variances
# and covariances may then also be found in the THETA matrix.
#
# pd_lv = "auto": the same, but a block that can not be handled (see the
# list above) is silently left alone instead of raising an error; the
# names of the variables in these blocks are returned in 'skipped'. The
# same holds for an optimizer that can not be used (gn, noniter).
#
# Note: the automatic bounds (bounds = "standard", "pos.var", ...) can not
# be imposed on the elements of the Cholesky factor: they are ignored for
# the parameters that are reparameterized.

# setup: find the PSI blocks that must be reparameterized, and check if we
# can handle them (if not, an error is raised)
#
# returns NULL if there is nothing to do, or a list with elements:
#   - blocks:  a list with one element per PSI block
#   - x_idx:   the positions (in the vector of free parameters) of all the
#              parameters that are reparameterized
#   - skipped: (pd_lv = "auto" only) the names of the variables in the
#              blocks that could not be handled
lav_model_pd_setup <- function(lavmodel = NULL, lavpartable = NULL,
                               lavoptions = NULL) {
  # pd_lv = "auto": skip what we can not handle
  auto <- identical(lavoptions$pd.lv, "auto")
  skipped <- character(0L)

  if (lavmodel@representation != "LISREL") {
    if (auto) {
      return(NULL)
    }
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

  # extra variables that must be treated as latent variables (the random
  # slopes in the M-step model of the random-slope EM)
  extra_names <- lav_model_pd_extra_names()

  blocks <- list()
  psi_idx <- which(names(glist) == "psi")
  # THETA is only needed for the extra variables
  theta_idx <- integer(0L)
  if (length(extra_names) > 0L) {
    theta_idx <- which(names(glist) == "theta")
  }
  mat_idx <- c(psi_idx, theta_idx)
  mat_block <- c(seq_along(psi_idx), seq_along(theta_idx))
  mat_psi <- c(rep(TRUE, length(psi_idx)), rep(FALSE, length(theta_idx)))
  for (m in seq_along(mat_idx)) {
    mm <- mat_idx[m]
    b <- mat_block[m]
    psi <- glist[[mm]]
    n <- nrow(psi)
    if (n == 0L || length(lavmodel@m.free.idx[[mm]]) == 0L) {
      next
    }
    psi_names <- lavmodel@dimNames[[mm]][[1L]]
    lv_names <- extra_names
    if (mat_psi[m]) {
      lv_names <- c(
        lv_names,
        lav_pt_vnames(lavpartable, type = "lv.regular", block = b),
        lav_pt_vnames(lavpartable, type = "lv.rv", block = b)
      )
    }

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

      # 1. variances: free, or fixed to a positive value
      var_free <- diag(id_c) > 0L
      bad <- which(!var_free & !(diag(psi_c) > 0))
      if (length(bad) > 0L) {
        if (auto) {
          skipped <- c(skipped, v_names)
          next
        }
        lav_msg_stop(gettextf(
          "pd_lv = TRUE is not available if the variance of a latent
          variable that has free covariances is fixed to zero (or to a
          negative value); the following variances are not positive%1$s:
          %2$s.", where_txt,
          lav_msg_view(v_names[bad], "none", FALSE)))
      }

      # 2. covariances: free, or fixed to zero
      low <- lower.tri(id_c)
      bad <- which(low & id_c == 0L & psi_c != 0, arr.ind = TRUE)
      if (nrow(bad) > 0L) {
        if (auto) {
          skipped <- c(skipped, v_names)
          next
        }
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
        if (auto) {
          skipped <- c(skipped, v_names)
          next
        }
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
        if (auto) {
          skipped <- c(skipped, v_names)
          next
        }
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
        var_free = var_free,
        var_fixed = ifelse(var_free, NA_real_, diag(psi_c)),
        diag_x = diag(id_c), # 0 = fixed variance
        low_row = unname(free_pos[, 1L]),
        low_sinh = !var_free[free_pos[, 1L]], # row with a fixed variance
        low_pos = free_pos,
        low_x = id_c[free_pos],
        zero_cols = zero_cols,
        zero_flag = any(lengths(zero_cols) > 0L)
      )
    } # components
  } # psi matrices

  if (auto && length(skipped) > 0L && lav_verbose()) {
    cat("pd_lv = \"auto\": skipped the block(s) with:",
        paste(unique(skipped), collapse = " "), "\n")
  }

  if (length(blocks) == 0L) {
    return(NULL)
  }

  list(
    skipped = unique(skipped),
    blocks = blocks,
    x_idx = unlist(lapply(blocks, function(blk) {
      c(blk$diag_x[blk$var_free], blk$low_x)
    }))
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
#
# returns a list with L, and the 'unnormalized' rows W: for a variable with
# a free variance, the rows of W and L coincide; for a variable with a
# fixed variance v, row i of L is sqrt(v) * W[i, ] / ||W[i, ]||
lav_model_pd_block_l <- function(x_pd, blk) {
  n <- blk$n
  w_mat <- matrix(0, n, n)
  w_low <- x_pd[blk$low_x]
  w_low[blk$low_sinh] <- sinh(w_low[blk$low_sinh])
  w_mat[blk$low_pos] <- w_low
  w_diag <- rep(1.0, n)
  w_diag[blk$var_free] <- exp(x_pd[blk$diag_x[blk$var_free]])
  diag(w_mat) <- w_diag

  # nothing special: L = W
  if (!blk$zero_flag && all(blk$var_free)) {
    return(list(l = w_mat, w = w_mat, w_norm = rep(1.0, n)))
  }

  l_mat <- w_mat
  w_norm <- rep(1.0, n)
  for (i in seq_len(n)) {
    # elements that are implied by the fixed zeros in PSI
    for (j in blk$zero_cols[[i]]) {
      k <- seq_len(j - 1L)
      w_mat[i, j] <- -sum(w_mat[i, k] * l_mat[j, k]) / l_mat[j, j]
    }
    # fixed variance: rescale
    if (blk$var_free[i]) {
      l_mat[i, ] <- w_mat[i, ]
    } else {
      w_norm[i] <- sqrt(sum(w_mat[i, ]^2))
      l_mat[i, ] <- sqrt(blk$var_fixed[i]) * w_mat[i, ] / w_norm[i]
    }
  }

  list(l = l_mat, w = w_mat, w_norm = w_norm)
}

# 'Cholesky' metric -> original (model) metric
lav_model_pd_to_model <- function(x_pd, pd) {
  x <- x_pd
  for (blk in pd$blocks) {
    l_mat <- lav_model_pd_block_l(x_pd, blk)$l
    psi <- tcrossprod(l_mat)
    x[blk$diag_x[blk$var_free]] <- diag(psi)[blk$var_free]
    x[blk$low_x] <- psi[blk$low_pos]
  }
  x
}

# the PSI block, given the parameters in the original (model) metric
lav_model_pd_block_psi <- function(x, blk) {
  psi <- matrix(0, blk$n, blk$n)
  psi[blk$low_pos] <- x[blk$low_x]
  psi <- psi + t(psi)
  psi_var <- blk$var_fixed
  psi_var[blk$var_free] <- x[blk$diag_x[blk$var_free]]
  diag(psi) <- psi_var
  psi
}

# original (model) metric -> 'Cholesky' metric
#
# this is used for the starting values (and for the EM accelerators); if
# a PSI block is not positive definite, the covariances are shrunken
# towards zero (as little as possible) until it is
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

    # shrink the covariances (gently at first) until the block is
    # positive definite
    psi_off <- psi
    diag(psi_off) <- 0
    l_mat <- NULL
    for (shrink in c(1, 1 - 1e-08, 1 - 1e-06, 1 - 1e-04, 0.99, 0.5^(1:20))) {
      psi_try <- psi_off * shrink
      diag(psi_try) <- psi_var
      l_try <- try(t(chol(psi_try)), silent = TRUE)
      if (!inherits(l_try, "try-error") && all(is.finite(l_try)) &&
          all(diag(l_try) > 0)) {
        l_mat <- l_try
        break
      }
    }
    if (is.null(l_mat)) {
      l_mat <- diag(sqrt(psi_var), nrow = blk$n)
    }
    # fixed variances: unnormalized rows, with a unit diagonal
    w_mat <- l_mat
    fixed <- which(!blk$var_free)
    if (length(fixed) > 0L) {
      w_mat[fixed, ] <- l_mat[fixed, , drop = FALSE] / diag(l_mat)[fixed]
    }
    x_pd[blk$diag_x[blk$var_free]] <- log(diag(l_mat)[blk$var_free])
    w_low <- w_mat[blk$low_pos]
    w_low[blk$low_sinh] <- asinh(w_low[blk$low_sinh])
    x_pd[blk$low_x] <- w_low
  }
  x_pd
}

# chain rule: the gradient with respect to the free parameters in the
# original (model) metric (dx) -> the gradient in the 'Cholesky' metric
lav_model_pd_gradient <- function(dx, x_pd, pd) {
  dx_pd <- dx
  for (blk in pd$blocks) {
    n <- blk$n
    tmp <- lav_model_pd_block_l(x_pd, blk)
    l_mat <- tmp$l
    w_mat <- tmp$w

    # dF/dPSI, with psi[i,j] and psi[j,i] as two distinct elements
    # (fixed variances: zero)
    g_mat <- matrix(0, n, n)
    g_mat[blk$low_pos] <- dx[blk$low_x] / 2
    g_mat <- g_mat + t(g_mat)
    g_diag <- numeric(n)
    g_diag[blk$var_free] <- dx[blk$diag_x[blk$var_free]]
    diag(g_mat) <- g_diag

    # dF/dL, all (lower-triangular) elements of L being independent
    l_bar <- 2 * g_mat %*% l_mat
    l_bar[upper.tri(l_bar)] <- 0

    # reverse sweep (last row first): when we arrive at row i, l_bar[i, ]
    # is complete; go back from L[i, ] to W[i, ] (fixed variances), and
    # then to the elements of W[i, ] that are parameters (the elements
    # implied by the fixed zeros in PSI are functions of the elements that
    # precede them, and of the previous rows of L)
    w_bar <- l_bar
    if (blk$zero_flag || !all(blk$var_free)) {
      for (i in rev(seq_len(n))) {
        wb <- l_bar[i, ]
        if (!blk$var_free[i]) {
          w_i <- w_mat[i, ]
          w_n <- tmp$w_norm[i]
          wb <- (sqrt(blk$var_fixed[i]) / w_n) *
            (wb - sum(wb * w_i) * w_i / (w_n * w_n))
        }
        for (j in rev(blk$zero_cols[[i]])) {
          a <- wb[j] / l_mat[j, j]
          k <- seq_len(j - 1L)
          wb[k] <- wb[k] - a * l_mat[j, k]
          l_bar[j, k] <- l_bar[j, k] - a * w_mat[i, k]
          l_bar[j, j] <- l_bar[j, j] - a * w_mat[i, j]
        }
        w_bar[i, ] <- wb
      }
    }

    w_low <- w_bar[blk$low_pos]
    # fixed variances: asinh scale
    w_low[blk$low_sinh] <- w_low[blk$low_sinh] *
      cosh(x_pd[blk$low_x[blk$low_sinh]])
    dx_pd[blk$low_x] <- w_low
    # free variances: log scale
    dx_pd[blk$diag_x[blk$var_free]] <-
      (diag(w_bar) * diag(w_mat))[blk$var_free]
  }
  dx_pd
}

# parameter scaling (optim.parscale): the scaling factors for the
# parameters in the 'Cholesky' metric, given the scaling factors of the
# original parameters (variance: s_i^2, covariance: s_i * s_j)
lav_model_pd_parscale <- function(parscale, pd) {
  for (blk in pd$blocks) {
    s2 <- rep(1.0, blk$n)
    s2[blk$var_free] <- parscale[blk$diag_x[blk$var_free]]
    if (any(!is.finite(s2)) || any(s2 <= 0)) {
      s2 <- rep(1.0, blk$n)
    }
    # log(L[i,i]): the scaling is a shift, not a factor
    parscale[blk$diag_x[blk$var_free]] <- 1.0
    # L[i,j]: scaled like the sd of variable i; fixed variance: W[i,j] is
    # a ratio, on an asinh scale (no scaling)
    parscale[blk$low_x] <- sqrt(s2)[blk$low_row]
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
# The tolerances (1e-04) are deliberately not too strict: when the
# likelihood is flat (eg a random slope without variance, few clusters),
# the optimizer stops at a relative variance of 1e-05 or so; if such a
# solution is not flagged, the information matrix is evaluated at a point
# that is not a stationary point in the original metric, and the standard
# errors may not exist (a Monte Carlo study with random-slope models:
# NaN standard errors in up to 2.4% of the replications with 1e-06, none
# of them flagged).
#
# For a flagged block, we also return the jacobian of the constraints that
# are active at the solution. If the columns of U span the null space of
# the (singular) PSI block, the matrices in the neighbourhood of PSI that
# have the same rank are (to first order) PSI + D with U' D U = 0: this is
# the tangent space of the boundary. Each pair of null vectors (u_a, u_b)
# therefore gives one linear constraint on the free elements of PSI:
#     u_a' D u_b = sum_i u_ai u_bi d_ii
#                + sum_{i>j} (u_ai u_bj + u_aj u_bi) d_ij = 0
# (a single null vector if the rank deficiency is one; for a variance that
# has vanished, the null vector is a unit vector, and the constraint simply
# fixes that variance; for a correlation of one between two standardized
# variables, the constraint fixes their covariance). These rows are added
# to the active (inequality) constraints of the model, so that the
# standard errors are computed on this tangent space, exactly as for a
# parameter that sits on a (simple) bound.
#
# Second-order term: the boundary is curved (unless a variance vanishes),
# and the hessian of the objective function alone does not describe the
# curvature of the objective function along the boundary. If LAMBDA is the
# (positive semi-definite) matrix of Lagrange multipliers, dF/dPSI = LAMBDA
# = U M U' at the solution, the quadratic form that matters for a tangent
# direction D is
#     (1/2) D' H D + tr(LAMBDA D PSI^+ D)
# (PSI^+ is the Moore-Penrose inverse of PSI; this is the 'sigma term' of
# the second-order conditions of semidefinite programming). The second
# term is returned as a matrix ('hess', to be added to the hessian of the
# objective function when the observed information is computed); it needs
# the gradient of the objective function with respect to the free
# parameters (dx), from which M is obtained by least squares. Without
# this term, the observed information restricted to the tangent space can
# be indefinite (no standard errors).
#
# returns a list with elements:
#   - names: the names of the variables in the offending blocks
#   - jac:   the jacobian of the active constraints (one row per
#            constraint, one column per free parameter)
#   - hess:  the second-order term (NULL if dx is not given)
lav_model_pd_boundary <- function(x, pd, var_ref = NULL, dx = NULL,
                                  tol_cor = 1e-04, tol_var = 1e-04) {
  out_names <- character(0L)
  out_jac <- matrix(0, 0L, length(x))
  out_hess <- NULL
  if (!is.null(dx) && length(dx) == length(x) && all(is.finite(dx))) {
    out_hess <- matrix(0, length(x), length(x))
  }

  for (blk in pd$blocks) {
    n <- blk$n
    psi <- lav_model_pd_block_psi(x, blk)
    psi_var <- diag(psi)
    if (any(!is.finite(psi))) {
      out_names <- c(out_names, blk$names)
      next
    }

    # null vectors (columns of u_mat)
    u_mat <- matrix(0, n, 0L)

    # 1. variances that have vanished
    zero_idx <- which(psi_var <= 0)
    if (length(var_ref) >= blk$block && isTRUE(var_ref[blk$block] > 0)) {
      zero_idx <- which(psi_var < tol_var * var_ref[blk$block])
    }
    if (length(zero_idx) > 0L) {
      u_mat <- diag(n)[, zero_idx, drop = FALSE]
    }

    # 2. the other variables: (almost) singular correlation matrix
    rest_idx <- setdiff(seq_len(n), zero_idx)
    if (length(rest_idx) > 1L) {
      sd_rest <- sqrt(psi_var[rest_idx])
      ev <- eigen(cov2cor(psi[rest_idx, rest_idx, drop = FALSE]),
                  symmetric = TRUE)
      null_idx <- which(ev$values < tol_cor)
      if (length(null_idx) > 0L) {
        u_rest <- ev$vectors[, null_idx, drop = FALSE] / sd_rest
        u_rest <- t(t(u_rest) / sqrt(colSums(u_rest * u_rest)))
        u_new <- matrix(0, n, length(null_idx))
        u_new[rest_idx, ] <- u_rest
        u_mat <- cbind(u_mat, u_new)
      }
    }

    if (ncol(u_mat) == 0L) {
      next
    }
    out_names <- c(out_names, blk$names)

    # one constraint per pair of null vectors
    vf <- which(blk$var_free)
    li <- blk$low_pos[, 1L]
    lj <- blk$low_pos[, 2L]
    blk_x <- c(blk$diag_x[vf], blk$low_x)
    pair_rows <- matrix(0, 0L, length(blk_x))
    pair_a <- pair_b <- integer(0L)
    for (a in seq_len(ncol(u_mat))) {
      for (b in seq_len(a)) {
        row <- numeric(length(x))
        row[blk$diag_x[vf]] <- u_mat[vf, a] * u_mat[vf, b]
        row[blk$low_x] <- u_mat[li, a] * u_mat[lj, b] +
          u_mat[lj, a] * u_mat[li, b]
        if (sqrt(sum(row * row)) > 1e-10) {
          out_jac <- rbind(out_jac, row, deparse.level = 0L)
          pair_rows <- rbind(pair_rows, row[blk_x], deparse.level = 0L)
          pair_a <- c(pair_a, a)
          pair_b <- c(pair_b, b)
        }
      }
    }

    # second-order term
    if (!is.null(out_hess) && nrow(pair_rows) > 0L) {
      # multipliers: dx = sum_c coef_c * row_c (least squares)
      coef_c <- try(drop(MASS::ginv(t(pair_rows)) %*% dx[blk_x]),
                    silent = TRUE)
      if (!inherits(coef_c, "try-error") && all(is.finite(coef_c))) {
        k <- ncol(u_mat)
        m_mat <- matrix(0, k, k)
        for (cc in seq_along(coef_c)) {
          a <- pair_a[cc]
          b <- pair_b[cc]
          if (a == b) {
            m_mat[a, a] <- coef_c[cc]
          } else {
            m_mat[a, b] <- m_mat[b, a] <- coef_c[cc] / 2
          }
        }
        # LAMBDA must be positive semi-definite
        m_ev <- eigen(m_mat, symmetric = TRUE)
        m_val <- pmax(m_ev$values, 0)
        lambda <- u_mat %*% (m_ev$vectors %*% (m_val * t(m_ev$vectors))) %*%
          t(u_mat)
        # PSI^+, on the orthogonal complement of the null space
        u_orth <- qr.Q(qr(u_mat))
        proj <- diag(n) - tcrossprod(u_orth)
        p_ev <- eigen(proj %*% psi %*% proj, symmetric = TRUE)
        keep <- which(p_ev$values > 1e-08 * max(p_ev$values))
        if (any(m_val > 0) && length(keep) > 0L) {
          psi_plus <- p_ev$vectors[, keep, drop = FALSE] %*%
            (t(p_ev$vectors[, keep, drop = FALSE]) / p_ev$values[keep])
          # basis matrices E_p (dPSI / dx_p) of the free elements
          npb <- length(blk_x)
          le <- vector("list", npb) # LAMBDA %*% E_p
          pe <- vector("list", npb) # PSI^+ %*% E_p
          for (p in seq_len(npb)) {
            e_mat <- matrix(0, n, n)
            if (p <= length(vf)) {
              e_mat[vf[p], vf[p]] <- 1
            } else {
              q <- p - length(vf)
              e_mat[li[q], lj[q]] <- e_mat[lj[q], li[q]] <- 1
            }
            le[[p]] <- lambda %*% e_mat
            pe[[p]] <- psi_plus %*% e_mat
          }
          h_blk <- matrix(0, npb, npb)
          for (p in seq_len(npb)) {
            for (q in seq_len(npb)) {
              # tr(LAMBDA E_p PSI^+ E_q)
              h_blk[p, q] <- sum(le[[p]] * t(pe[[q]]))
            }
          }
          out_hess[blk_x, blk_x] <- out_hess[blk_x, blk_x] +
            h_blk + t(h_blk)
        }
      }
    }
  }

  # remove redundant rows (if any)
  if (nrow(out_jac) > 1L) {
    qr_jac <- qr(t(out_jac))
    out_jac <- out_jac[qr_jac$pivot[seq_len(qr_jac$rank)], , drop = FALSE]
  }

  list(names = unique(out_names), jac = out_jac, hess = out_hess)
}

# extra variables that must be treated as latent variables when the PSI
# (and THETA) blocks are selected; only used by the M-step of the
# random-slope EM, where the random slopes are observed variables: the
# names are stored in lavaan_cache_env during the M-step fit
#
# if 'value' is given, the names are set (NULL: removed), and the previous
# value is returned
lav_model_pd_extra_names <- function(value) {
  old <- get0("pd_lv_extra_names", lavaan_cache_env,
              ifnotfound = character(0L))
  if (!missing(value)) {
    if (length(value) == 0L) {
      if (exists("pd_lv_extra_names", envir = lavaan_cache_env,
                 inherits = FALSE)) {
        rm(list = "pd_lv_extra_names", envir = lavaan_cache_env)
      }
    } else {
      assign("pd_lv_extra_names", value, envir = lavaan_cache_env)
    }
    return(invisible(old))
  }
  old
}

# are all the PSI blocks (x: original metric) positive semi-definite?
# (up to a tolerance, relative to the size of the variances)
lav_model_pd_admissible <- function(x, pd, tol = 1e-08) {
  for (blk in pd$blocks) {
    psi <- lav_model_pd_block_psi(x, blk)
    if (any(!is.finite(psi))) {
      return(FALSE)
    }
    ev <- eigen(psi, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -tol * max(abs(diag(psi)), .Machine$double.eps)) {
      return(FALSE)
    }
  }
  TRUE
}

# EM acceleration (squarem, qn): the accelerators extrapolate in the
# parameter space, and an extrapolated point may well be outside the
# admissible region; with pd_lv = TRUE, the extrapolation is done in the
# 'Cholesky' metric instead, where every point is admissible
#
# returns a list with the wrapped functions, and the two conversions
lav_model_pd_em_wrap <- function(pd, step_fn, logl_fn, conv_fn) {
  to_model <- function(theta) lav_model_pd_to_model(theta, pd = pd)
  to_chol <- function(x) {
    x_attr <- attributes(x)
    out <- lav_model_pd_to_cholesky(as.numeric(x), pd = pd)
    # keep the "logl" attribute (a byproduct of the E-step)
    if (!is.null(x_attr$logl)) {
      attr(out, "logl") <- x_attr$logl
    }
    out
  }
  list(
    to_model = to_model,
    to_chol = to_chol,
    step_fn = function(theta, ...) to_chol(step_fn(to_model(theta), ...)),
    logl_fn = function(theta) logl_fn(to_model(theta)),
    conv_fn = function(theta_old, theta, fx_old, fx) {
      conv_fn(to_model(theta_old), to_model(theta), fx_old, fx)
    }
  )
}

# the M-step parameter table: remove the bounds of the parameters that are
# reparameterized (the automatic bounds of the original model would be
# mistaken for user-specified bounds in the M-step fit)
lav_model_pd_mstep_partable <- function(partable, pd) {
  idx <- which(partable$free %in% pd$x_idx)
  if (length(idx) > 0L) {
    if (!is.null(partable$lower)) {
      partable$lower[idx] <- -Inf
    }
    if (!is.null(partable$upper)) {
      partable$upper[idx] <- +Inf
    }
  }
  partable
}

# is the pd_lv option switched on (TRUE, or "auto")?
lav_model_pd_on <- function(lavoptions = NULL) {
  isTRUE(lavoptions$pd.lv) || identical(lavoptions$pd.lv, "auto")
}
