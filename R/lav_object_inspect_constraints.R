# lavInspect() workers for constraint information:
#   "constraints"            (alias "con")      -- annotated list
#   "constraints.jacobian"   (alias "con.jac")  -- stacked constraint Jacobian
#   "constraints.nullspace"  (aliases "constraints.k", "eq.constraints.k",
#                             "con.k")          -- null-space basis 'K'
#
# Everything is RECONSTRUCTED from the fitted model (constraint functions
# evaluated at the final estimates, lav_con_eq_basis() for 'K'), never
# copied from the internal packing slots: which slot holds the constraint
# information depends on the constraint regime (eq.constraints /
# ceq.simple.only / general), and some slots (eg ceq.JAC) are evaluated at
# starting values.
#
# Two parameter spaces show up in the output. Normally they coincide, but
# when simple equality constraints are absorbed into the parameter vector
# (ceq.simple = TRUE), the 'free' vector that the constraint machinery
# sees is the reduced (packed) one, while bound constraints operate on the
# expanded vector with one entry per free parameter-table row. The column
# names of the returned matrices always identify the space.

# gather all constraint information; returns a plain (unlabeled) list
lav_inspect_con_info <- function(object) {
  lavmodel <- object@Model
  partable <- object@ParTable

  # all free rows; when simple equality constraints are absorbed
  # (ceq.simple = TRUE), several rows share the same 'free' index
  free_row_idx <- which(partable$free > 0L)
  free_idx <- partable$free[free_row_idx]
  npar <- lavmodel@nx.free
  absorbed <- any(duplicated(free_idx))

  # parameter labels/values in the expanded (one per free row) and the
  # reduced (one per free index) space; identical unless absorbed
  labels_unco <- lav_pt_labels(partable, type = "free")
  labels_user <- lav_pt_labels(partable, type = "user")
  est_user <- lav_inspect_est(object)
  x_unco <- est_user[free_row_idx]
  first_idx <- match(seq_len(npar), free_idx)
  labels_free <- labels_unco[first_idx]
  x_free <- x_unco[first_idx]

  # ---------------- equality constraints ----------------
  # the ceq function returns lhs - rhs per declared "==" row, and indexes
  # the reduced parameter vector
  ceq_fun <- lavmodel@ceq.function
  if (!is.null(body(ceq_fun))) {
    ceq_resid <- as.numeric(ceq_fun(x_free))
    ceq_jac <- try(lav_func_jacobian_complex(func = ceq_fun, x = x_free),
      silent = TRUE)
    if (inherits(ceq_jac, "try-error")) { # eg. pnorm()
      ceq_jac <- lav_func_jacobian_simple(func = ceq_fun, x = x_free)
    }
    ceq_rhs <- try(-1 * as.numeric(ceq_fun(numeric(npar))), silent = TRUE)
    if (inherits(ceq_rhs, "try-error")) {
      ceq_rhs <- rep(as.numeric(NA), length(ceq_resid))
    }
  } else {
    ceq_jac <- matrix(0, nrow = 0L, ncol = npar)
    ceq_resid <- numeric(0L)
    ceq_rhs <- numeric(0L)
  }
  n_ceq_fun <- nrow(ceq_jac)

  # row labels: the declared "==" rows, in parameter-table order
  eq_pt_idx <- which(partable$op == "==")
  if (length(eq_pt_idx) == 0L && n_ceq_fun == 0L) {
    ceq_names <- character(0L)
  } else if (length(eq_pt_idx) == n_ceq_fun) {
    ceq_names <- paste(
      partable$lhs[eq_pt_idx], "==", partable$rhs[eq_pt_idx])
  } else {
    # unused (all-zero) constraint rows may have been pruned internally
    ceq_names <- paste0("eq", seq_len(n_ceq_fun))
  }

  # EFA models: the rotation (criterion stationarity) constraints are not
  # part of ceq.function, but they do restrict the parameters and enter
  # the df/vcov bookkeeping; append them, clearly flagged
  efa_jac <- lavmodel@ceq.efa.JAC
  if (nrow(efa_jac) > 0L && ncol(efa_jac) == ncol(ceq_jac)) {
    ceq_jac <- rbind(ceq_jac, efa_jac)
    ceq_names <- c(ceq_names, paste0("efa", seq_len(nrow(efa_jac))))
    ceq_resid <- c(ceq_resid, rep(as.numeric(NA), nrow(efa_jac)))
    ceq_rhs <- c(ceq_rhs, rep(as.numeric(NA), nrow(efa_jac)))
  }

  # number of (linearly independent) equality restrictions
  if (nrow(ceq_jac) > 0L) {
    ceq_rank <- qr(ceq_jac)$rank
  } else {
    ceq_rank <- 0L
  }

  # ---------------- inequality constraints ----------------
  # the cin function returns the slack (>= 0 at a feasible point) for the
  # finite upper bounds, then the finite lower bounds, then the explicit
  # "<"/">" rows; the bound rows index the expanded parameter vector
  cin_fun <- lavmodel@cin.function
  if (absorbed) {
    x_cin <- x_unco
    cin_par_labels <- labels_unco
  } else {
    x_cin <- x_free
    cin_par_labels <- labels_free
  }
  if (!is.null(body(cin_fun))) {
    cin_slack <- as.numeric(cin_fun(x_cin))
    cin_jac <- try(lav_func_jacobian_complex(func = cin_fun, x = x_cin),
      silent = TRUE)
    if (inherits(cin_jac, "try-error")) { # eg. pnorm()
      cin_jac <- lav_func_jacobian_simple(func = cin_fun, x = x_cin)
    }
    cin_rhs <- try(-1 * as.numeric(cin_fun(numeric(length(x_cin)))),
      silent = TRUE)
    if (inherits(cin_rhs, "try-error")) {
      cin_rhs <- rep(as.numeric(NA), length(cin_slack))
    }
  } else {
    cin_jac <- matrix(0, nrow = 0L, ncol = length(x_cin))
    cin_slack <- numeric(0L)
    cin_rhs <- numeric(0L)
  }
  n_cin <- nrow(cin_jac)

  # row labels, in the same order as the cin function output
  upper_pt_idx <- lower_pt_idx <- integer(0L)
  if (!is.null(partable$upper)) {
    upper_pt_idx <- which(partable$free > 0L & is.finite(partable$upper))
  }
  if (!is.null(partable$lower)) {
    lower_pt_idx <- which(partable$free > 0L & is.finite(partable$lower))
  }
  ineq_pt_idx <- which(partable$op %in% c("<", ">"))
  if (length(upper_pt_idx) + length(lower_pt_idx) +
      length(ineq_pt_idx) == n_cin) {
    # note: guard each part against zero-length indices (paste() would
    # recycle them up to length one)
    cin_names <- character(0L)
    if (length(upper_pt_idx) > 0L) {
      cin_names <- c(cin_names, paste(labels_user[upper_pt_idx], "<=",
        partable$upper[upper_pt_idx]))
    }
    if (length(lower_pt_idx) > 0L) {
      cin_names <- c(cin_names, paste(labels_user[lower_pt_idx], ">=",
        partable$lower[lower_pt_idx]))
    }
    if (length(ineq_pt_idx) > 0L) {
      cin_names <- c(cin_names, paste(partable$lhs[ineq_pt_idx],
        partable$op[ineq_pt_idx], partable$rhs[ineq_pt_idx]))
    }
  } else {
    cin_names <- paste0("ineq", seq_len(n_cin))
  }

  # ---------------- activity + Lagrange multipliers ----------------
  # the optimizer records (in con.jac/con.lambda) which inequality rows
  # were binding at the solution; when that record is available and its
  # layout matches the declared rows, use it -- otherwise fall back to a
  # conservative slack-only classification
  cin_active <- cin_slack < 1e-05
  cin_lambda <- rep(as.numeric(NA), n_cin)
  ceq_lambda <- rep(as.numeric(NA), n_ceq_fun)
  con_jac_slot <- lavmodel@con.jac
  if (nrow(con_jac_slot) > 0L) {
    ceq_idx_attr <- attr(con_jac_slot, "ceq.idx")
    cin_idx_attr <- attr(con_jac_slot, "cin.idx")
    inactive_attr <- attr(con_jac_slot, "inactive.idx")
    lambda <- lavmodel@con.lambda
    lambda_ok <- length(lambda) == nrow(con_jac_slot)
    if (length(cin_idx_attr) == n_cin && n_cin > 0L) {
      active <- rep(TRUE, n_cin)
      if (length(inactive_attr) > 0L) {
        active[match(inactive_attr, cin_idx_attr)] <- FALSE
      }
      cin_active <- active
      if (lambda_ok) {
        cin_lambda <- lambda[cin_idx_attr]
      }
    }
    if (length(ceq_idx_attr) == n_ceq_fun && n_ceq_fun > 0L && lambda_ok) {
      ceq_lambda <- lambda[ceq_idx_attr]
    }
  }
  if (nrow(ceq_jac) > n_ceq_fun) {
    # appended EFA rows
    ceq_lambda <- c(ceq_lambda,
      rep(as.numeric(NA), nrow(ceq_jac) - n_ceq_fun))
  }

  # ---------------- null-space basis 'K' ----------------
  # the basis actually used by the test/vcov machinery, valid in every
  # constraint regime; post-multiplying the (expanded) Delta matrix by K
  # restricts it to the directions that satisfy the equality constraints
  k <- lav_con_eq_basis(lavmodel)
  if (is.null(k)) {
    # no equality constraints: the constrained tangent space is the full
    # parameter space
    k <- diag(npar)
  }
  if (nrow(k) == length(labels_unco)) {
    k_par_labels <- labels_unco
  } else if (nrow(k) == npar) {
    k_par_labels <- labels_free
  } else {
    k_par_labels <- NULL
  }

  # offset vector: only available in the linear-equality-only regime,
  # where theta = k0 + K %*% theta.reduced
  if (lavmodel@eq.constraints) {
    k0 <- lavmodel@eq.constraints.k0
  } else {
    k0 <- NULL
  }

  list(
    ceq.jac = ceq_jac, ceq.names = ceq_names, ceq.rhs = ceq_rhs,
    ceq.resid = ceq_resid, ceq.rank = ceq_rank, ceq.lambda = ceq_lambda,
    cin.jac = cin_jac, cin.names = cin_names, cin.rhs = cin_rhs,
    cin.slack = cin_slack, cin.active = cin_active,
    cin.lambda = cin_lambda,
    k = k, k.par.labels = k_par_labels, k0 = k0,
    ceq.simple = absorbed,
    labels.free = labels_free, cin.par.labels = cin_par_labels
  )
}

# what = "constraints" (alias "con")
lav_inspect_con <- function(object, add_labels = FALSE,
                                   add_class = FALSE) {
  info <- lav_inspect_con_info(object)

  ceq_jac <- info$ceq.jac
  cin_jac <- info$cin.jac
  ceq_rhs <- info$ceq.rhs
  ceq_resid <- info$ceq.resid
  ceq_lambda <- info$ceq.lambda
  cin_rhs <- info$cin.rhs
  cin_slack <- info$cin.slack
  cin_active <- info$cin.active
  cin_lambda <- info$cin.lambda
  k <- info$k
  k0 <- info$k0

  if (add_labels) {
    dimnames(ceq_jac) <- list(info$ceq.names, info$labels.free)
    dimnames(cin_jac) <- list(info$cin.names, info$cin.par.labels)
    names(ceq_rhs) <- names(ceq_resid) <- info$ceq.names
    names(ceq_lambda) <- info$ceq.names
    names(cin_rhs) <- names(cin_slack) <- info$cin.names
    names(cin_active) <- names(cin_lambda) <- info$cin.names
    rownames(k) <- info$k.par.labels
    if (!is.null(k0)) {
      names(k0) <- info$labels.free
    }
  }

  if (add_class) {
    class(ceq_jac) <- c("lavaan.matrix", "matrix")
    class(cin_jac) <- c("lavaan.matrix", "matrix")
    class(k) <- c("lavaan.matrix", "matrix")
    class(ceq_rhs) <- class(ceq_resid) <- c("lavaan.vector", "numeric")
    class(ceq_lambda) <- c("lavaan.vector", "numeric")
    class(cin_rhs) <- class(cin_slack) <- c("lavaan.vector", "numeric")
    class(cin_lambda) <- c("lavaan.vector", "numeric")
    if (!is.null(k0)) {
      class(k0) <- c("lavaan.vector", "numeric")
    }
  }

  list(
    ceq.jac = ceq_jac,
    ceq.rhs = ceq_rhs,
    ceq.resid = ceq_resid,
    ceq.rank = info$ceq.rank,
    ceq.lambda = ceq_lambda,
    cin.jac = cin_jac,
    cin.rhs = cin_rhs,
    cin.slack = cin_slack,
    cin.active = cin_active,
    cin.lambda = cin_lambda,
    k = k,
    k0 = k0,
    ceq.simple = info$ceq.simple
  )
}

# what = "constraints.jacobian" (alias "con.jac"): the equality rows
# stacked above the inequality rows, evaluated at the final estimates,
# with an "active" attribute (equality rows are always active)
lav_inspect_con_jac <- function(object, add_labels = FALSE,
                                       add_class = FALSE) {
  info <- lav_inspect_con_info(object)

  n_ceq <- nrow(info$ceq.jac)
  n_cin <- nrow(info$cin.jac)
  if (n_cin == 0L) {
    out <- info$ceq.jac
    par_labels <- info$labels.free
  } else if (n_ceq == 0L) {
    out <- info$cin.jac
    par_labels <- info$cin.par.labels
  } else if (ncol(info$ceq.jac) == ncol(info$cin.jac)) {
    out <- rbind(info$ceq.jac, info$cin.jac)
    par_labels <- info$labels.free
  } else {
    # different parameter spaces (simple equalities absorbed + bounds in
    # the expanded space); return the equality part only
    lav_msg_warn(gettext(
      "equality and inequality constraints live in different parameter spaces for this model; returning the equality constraints only (use what = \"constraints\" to see both)."))
    out <- info$ceq.jac
    par_labels <- info$labels.free
    n_cin <- 0L
  }

  active <- c(rep(TRUE, n_ceq), if (n_cin > 0L) info$cin.active)
  ceq_idx <- seq_len(n_ceq)
  cin_idx <- n_ceq + seq_len(n_cin)

  if (add_labels) {
    row_names <- c(if (n_ceq > 0L) info$ceq.names,
      if (n_cin > 0L) info$cin.names)
    dimnames(out) <- list(row_names, par_labels)
  }
  if (add_class) {
    class(out) <- c("lavaan.matrix", "matrix")
  }
  attr(out, "active") <- active
  attr(out, "ceq.idx") <- ceq_idx
  attr(out, "cin.idx") <- cin_idx

  out
}

# what = "constraints.nullspace" (aliases "constraints.k",
# "eq.constraints.k", "con.k"): orthonormal basis of the null space of the
# equality-constraint Jacobian, as used by the test statistic and vcov
# machinery; the identity matrix if the model has no equality constraints
lav_inspect_con_nullspace <- function(object, add_labels = FALSE,
                                             add_class = FALSE) {
  info <- lav_inspect_con_info(object)

  k <- info$k
  if (add_labels) {
    rownames(k) <- info$k.par.labels
  }
  if (add_class) {
    class(k) <- c("lavaan.matrix", "matrix")
  }

  k
}
