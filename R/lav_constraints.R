# Evaluate 'expr' using a local, fixed RNG state, and restore the global
# '.Random.seed' afterwards. This keeps RNG-using helpers (eg the jacobian
# linearity checks below) from disturbing the global stream -- which would,
# for example, change the seed picked by a subsequent se = "bootstrap" run.
lav_with_local_seed <- function(expr, seed = 1234L) {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv,
                   inherits = FALSE), add = TRUE)
  } else {
    # no '.Random.seed' existed: remove the one set.seed() creates below
    on.exit(rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
            add = TRUE)
  }
  set.seed(seed)
  force(expr)
}

lav_con_parse <- function(partable = NULL, constraints = NULL,
                                  theta = NULL,
                                  debug = FALSE) {
  if (!missing(debug)) {
    current_debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current_debug), TRUE)
  }
  # just in case we do not have a $free column in partable
  if (is.null(partable$free)) {
    partable$free <- seq_along(partable$lhs)
  }

  # from the partable: free parameters
  if (!is.null(theta)) {
    # nothing to do
  } else if (!is.null(partable$est)) {
    theta <- partable$est[partable$free > 0L]
  } else if (!is.null(partable$start)) {
    theta <- partable$start[partable$free > 0L]
  } else {
    theta <- rep(0, length(partable$lhs))
  }

  # note: avoid zero values in theta (so jac is not all null)
  # (eg, regression coefficients are often zero in partable$start)
  # we perturb with small random numbers, but isolate the RNG state
  # (see lav_with_local_seed() above)
  zero_theta_idx <- which(theta == 0)
  if (length(zero_theta_idx) > 0L) {
    theta[zero_theta_idx] <- lav_with_local_seed(
      rnorm(length(zero_theta_idx), 0, 0.1))
  }

  # number of free (but possibly constrained) parameters
  npar <- length(theta)

  # parse the constraints
  if (is.null(constraints)) {
    list_1 <- NULL
  } else if (!is.character(constraints)) {
    lav_msg_stop(gettext("constraints should be a string"))
  } else {
    flat <- lavParseModelString(constraints)
    con <- attr(flat, "constraints")
    list_1 <- list()
    if (length(con) > 0L) {
      lhs <- unlist(lapply(con, "[[", "lhs"))
      op <- unlist(lapply(con, "[[", "op"))
      rhs <- unlist(lapply(con, "[[", "rhs"))
      list_1$lhs <- c(list_1$lhs, lhs)
      list_1$op <- c(list_1$op, op)
      list_1$rhs <- c(list_1$rhs, rhs)
    } else {
      lav_msg_stop(gettext("no constraints found in constraints argument"))
    }
  }

  # simple equality constraints?
  ceq_simple <- FALSE
  if (!is.null(partable$unco)) {
    ceq_simple <- TRUE
  }

  # simple inequality constraints?
  cin_simple <- TRUE
  if (any(partable$op %in% c(">", "<"))) {
    cin_simple <- FALSE
  }
  if (!is.null(list_1) && any(list_1$op %in% c(">", "<"))) {
    cin_simple <- FALSE
  }
  if (is.null(partable$upper) && is.null(partable$lower)) {
    cin_simple <- FALSE
  } else {
    # only check free parameters!
    free_idx <- which(partable$free > 0L)
    if (all(partable$lower[free_idx] == -Inf) &&
        all(partable$upper[free_idx] == +Inf)) {
      cin_simple <- FALSE
    }
  }

  # variable definitions
  def_function <- lav_pt_con_def(partable,
    con = list_1,
    debug = debug
  )

  # construct ceq/ciq functions
  ceq_function <- lav_pt_con_ceq(partable,
    con = list_1,
    debug = debug
  )
  # linear or nonlinear?
  ceq_linear_idx <- lav_con_linear_idx(
    func = ceq_function,
    npar = npar
  )
  ceq_nonlinear_idx <- lav_con_nonlinear_idx(
    func = ceq_function,
    npar = npar
  )

  # inequalities
  cin_function <- lav_pt_con_ciq(partable,
    con = list_1,
    debug = debug
  )

  # linear or nonlinear?
  cin_linear_idx <- lav_con_linear_idx(
    func = cin_function,
    npar = npar
  )
  cin_nonlinear_idx <- lav_con_nonlinear_idx(
    func = cin_function,
    npar = npar
  )

  # Jacobians
  if (!is.null(body(ceq_function))) {
    ceq_jac <- try(lav_func_jacobian_complex(
      func = ceq_function,
      x = theta
    ), silent = TRUE)
    if (inherits(ceq_jac, "try-error")) { # eg. pnorm()
      ceq_jac <- lav_func_jacobian_simple(func = ceq_function, x = theta)
    }

    # constants
    # do we have a non-zero 'rhs' elements? FIXME!!! is this reliable??
    ceq_rhs <- -1 * ceq_function(numeric(npar))

    # evaluate constraints
    ceq_theta <- ceq_function(theta)
  } else {
    ceq_jac <- matrix(0, nrow = 0L, ncol = npar)
    ceq_rhs <- numeric(0L)
    ceq_theta <- numeric(0L)
  }

  if (!is.null(body(cin_function))) {
    cin_jac <- try(lav_func_jacobian_complex(
      func = cin_function,
      x = theta
    ), silent = TRUE)
    if (inherits(cin_jac, "try-error")) { # eg. pnorm()
      cin_jac <- lav_func_jacobian_simple(func = cin_function, x = theta)
    }

    # constants
    # do we have a non-zero 'rhs' elements? FIXME!!! is this reliable??
    cin_rhs <- -1 * cin_function(numeric(npar))

    # evaluate constraints
    cin_theta <- cin_function(theta)
  } else {
    cin_jac <- matrix(0, nrow = 0L, ncol = npar)
    cin_rhs <- numeric(0L)
    cin_theta <- numeric(0L)
  }

  # check for empty/unused constraints (new in 0.7-1)
  if (nrow(ceq_jac) > 0L) {
    if (all(ceq_jac == 0)) {
      ceq_jac <- matrix(0, nrow = 0L, ncol = npar)
      ceq_rhs <- numeric(0L)
      ceq_theta <- numeric(0L)
      ceq_linear_idx <- integer(0L)
      ceq_nonlinear_idx <- integer(0L)
    } else {
      zero_idx <- which(apply(ceq_jac, 1, function(x) all(x == 0)))
      if (length(zero_idx) > 0L) {
        ceq_jac <- ceq_jac[-zero_idx, , drop = FALSE]
        ceq_rhs <- ceq_rhs[-zero_idx]
        ceq_theta <- ceq_theta[-zero_idx]
        # hm, how to handle these? indices no longer match rows of ceq.JAC!
        ceq_linear_idx <- ceq_linear_idx[!ceq_linear_idx %in% zero_idx]
        ceq_nonlinear_idx <- ceq_nonlinear_idx[!ceq_nonlinear_idx %in% zero_idx]
      }
    }
  }
  if (nrow(cin_jac) > 0L) {
    if (all(cin_jac == 0)) {
      cin_jac <- matrix(0, nrow = 0L, ncol = npar)
      cin_rhs <- numeric(0L)
      cin_theta <- numeric(0L)
      cin_linear_idx <- integer(0L)
      cin_nonlinear_idx <- integer(0L)
    } else {
      zero_idx <- which(apply(cin_jac, 1, function(x) all(x == 0)))
      if (length(zero_idx) > 0L) {
        cin_jac <- cin_jac[-zero_idx, , drop = FALSE]
        cin_rhs <- cin_rhs[-zero_idx]
        cin_theta <- cin_theta[-zero_idx]
        # hm, how to handle these? indices no longer match rows of cin.JAC!
        cin_linear_idx <- cin_linear_idx[!cin_linear_idx %in% zero_idx]
        cin_nonlinear_idx <- cin_nonlinear_idx[!cin_nonlinear_idx %in% zero_idx]
      }
    }
  }

  # shortcut flags
  ceq_linear_flag <- length(ceq_linear_idx) > 0L
  ceq_nonlinear_flag <- length(ceq_nonlinear_idx) > 0L
  ceq_flag <- ceq_linear_flag || ceq_nonlinear_flag

  cin_linear_flag <- length(cin_linear_idx) > 0L
  cin_nonlinear_flag <- length(cin_nonlinear_idx) > 0L
  cin_flag <- cin_linear_flag || cin_nonlinear_flag
  if (cin_simple) {
    cin_flag <- FALSE
  }

  # ceq_only_flag <- ceq_flag && !cin_flag
  cin_only_flag <- cin_flag && !ceq_flag

  ceq_linear_only_flag <- (ceq_linear_flag &&
    !ceq_nonlinear_flag &&
    !cin_flag && !cin_simple)

  ceq_simple_only <- ceq_simple && !ceq_flag && !cin_flag
  cin_simple_only <- cin_simple && !ceq_linear_flag

  # additional info if ceq.linear.flag
  if (ceq_linear_only_flag) {
    ## NEW: 18 nov 2014: handle general *linear* constraints
    ##
    ## see Nocedal & Wright (2006) 15.3
    ## - from x to x.red:
    ##       x.red <- MASS::ginv(Q2) %*% (x - Q1 %*% solve(t(R)) %*% b)
    ##   or
    ##       x.red <- as.numeric((x - b %*% qr.coef(QR,diag(npar))) %*% Q2)
    ##
    ## - from x.red to x
    ##       x <- as.numeric(Q1 %*% solve(t(R)) %*% b + Q2 %*% x.red)
    ##   or
    ##       x <- as.numeric(b %*% qr.coef(QR, diag(npar))) +
    ##                       as.numeric(Q2 %*% x.red)
    ##
    ## we write eq.constraints.K = Q2
    ##          eq.constraints.k0 = b %*% qr.coef(QR, diag(npar)))

    # compute range+null space of the jacobian (JAC) of the constraint
    # matrix
    # JAC <- lav_func_jacobian_complex(func = ceq.function,
    #           x = lavpartable$start[lavpartable$free > 0L]
    qr_1 <- qr(t(ceq_jac))
    ran_k <- qr_1$rank
    q_1 <- qr.Q(qr_1, complete = TRUE)
    # Q1 <- Q[,1:ranK, drop = FALSE]         # range space
    # Q2 <- Q[,-seq_len(ranK), drop = FALSE] # null space
    # R <- qr.R(QR)
    ceq_jac_null <- q_1[, -seq_len(ran_k), drop = FALSE]

    if (all(ceq_rhs == 0)) {
      ceq_rhs_null <- numeric(npar)
    } else {
      tmp <- qr.coef(qr_1, diag(npar))
      na_idx <- which(is.na(rowSums(tmp))) # catch NAs
      if (length(na_idx) > 0L) {
        tmp[na_idx, ] <- 0
      }
      ceq_rhs_null <- as.numeric(ceq_rhs %*% tmp)
    }
  } else {
    ceq_jac_null <- matrix(0, 0L, 0L)
    ceq_rhs_null <- numeric(0L)
  }

  # if simple equalities only, create 'K' matrix
  ceq_simple_k <- matrix(0, 0, 0)
  if (ceq_simple_only) {
    n_unco <- max(partable$unco)
    n_free <- max(partable$free)
    ceq_simple_k <- matrix(0, nrow = n_unco, ncol = n_free)
    #####
    #####     FIXME !
    #####
    idx_free <- partable$free[partable$free > 0]
    for (k in 1:n_unco) {
      c <- idx_free[k]
      ceq_simple_k[k, c] <- 1
    }
  }

  # dummy jacobian 'function'
  ceq_jacobian <- function() NULL
  cin_jacobian <- function() NULL


  out <- list(
    def.function = def_function,
    ceq.function = ceq_function,
    ceq.JAC = ceq_jac,
    ceq.jacobian = ceq_jacobian,
    ceq.rhs = ceq_rhs,
    ceq.theta = ceq_theta,
    ceq.linear.idx = ceq_linear_idx,
    ceq.nonlinear.idx = ceq_nonlinear_idx,
    ceq.linear.flag = ceq_linear_flag,
    ceq.nonlinear.flag = ceq_nonlinear_flag,
    ceq.flag = ceq_flag,
    ceq.linear.only.flag = ceq_linear_only_flag,
    ceq.JAC.NULL = ceq_jac_null,
    ceq.rhs.NULL = ceq_rhs_null,
    ceq.simple.only = ceq_simple_only,
    ceq.simple.K = ceq_simple_k,
    cin.function = cin_function,
    cin.JAC = cin_jac,
    cin.jacobian = cin_jacobian,
    cin.rhs = cin_rhs,
    cin.theta = cin_theta,
    cin.linear.idx = cin_linear_idx,
    cin.nonlinear.idx = cin_nonlinear_idx,
    cin.linear.flag = cin_linear_flag,
    cin.nonlinear.flag = cin_nonlinear_flag,
    cin.flag = cin_flag,
    cin.only.flag = cin_only_flag,
    cin.simple.only = cin_simple_only
  )

  out
}

lav_con_linear_idx <- function(func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) {
    return(integer(0L))
  }

  # compare the jacobian at two (distinct) evaluation points; isolate the
  # RNG state so we do not disturb the global stream (see lav_with_local_seed)
  a0min_a1 <- lav_with_local_seed({
    a0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar)) # point 1
    a1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar)) # point 2
    a0 - a1
  })
  linear <- apply(a0min_a1, 1, function(x) all(x == 0))
  which(linear)
}

lav_con_nonlinear_idx <- function(func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) {
    return(integer(0L))
  }

  # compare the jacobian at two (distinct) evaluation points; isolate the
  # RNG state so we do not disturb the global stream (see lav_with_local_seed)
  a0min_a1 <- lav_with_local_seed({
    a0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar)) # point 1
    a1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar)) # point 2
    a0 - a1
  })
  linear <- apply(a0min_a1, 1, function(x) all(x == 0))
  which(!linear)
}

# check if the equality constraints are 'simple' (a == b)
lav_con_check_simple <- function(lavmodel = NULL) {
  ones <- (lavmodel@ceq.JAC == 1 | lavmodel@ceq.JAC == -1)
  simple <- all(lavmodel@ceq.rhs == 0) &&
    all(apply(lavmodel@ceq.JAC != 0, 1, sum) == 2) &&
    all(apply(ones, 1, sum) == 2) &&
    length(lavmodel@ceq.nonlinear.idx) == 0

  # TRUE or FALSE
  simple
}

# orthonormal basis 'K' of the equality-constrained (tangent) space, in the
# column space of lav_model_delta(): post-multiplying Delta by K restricts it
# to the directions that satisfy the equality constraints.
#
# This plays the role of eq.constraints.K, but is also available when
# equality constraints coexist with inequality constraints or with parameter
# bounds: in that case BOTH lavmodel@eq.constraints and
# lavmodel@ceq.simple.only are FALSE (they are *packing* flags, meaning
# "only equality constraints"), eq.constraints.K is not computed, and code
# that gates the K-reduction on those flags would silently drop the equality
# constraints (wrong statistics with unchanged df).
#
# Returns NULL when the model has no equality constraints. For nonlinear
# equality constraints the basis is the null space of the Jacobian
# (a local, asymptotically correct linearization).
lav_con_eq_basis <- function(lavmodel = NULL) {
  if (lavmodel@eq.constraints) {
    # only linear equality constraints: precomputed orthonormal null space
    out <- lavmodel@eq.constraints.K
  } else if (lavmodel@ceq.simple.only) {
    # ceq.simple.K is a 0/1 duplication matrix (t(K) %*% K != I);
    # orthonormalize so that K / t(K) round-trips are consistent
    out <- qr.Q(qr(lavmodel@ceq.simple.K))
  } else if (nrow(lavmodel@ceq.JAC) > 0L) {
    # equality constraints together with inequality constraints/bounds,
    # or nonlinear equality constraints
    jac <- lavmodel@ceq.JAC
    # ceq.JAC is evaluated at the STARTING values; when nonlinear equality
    # constraints are present, the linearization must be taken at the
    # solution, so use the optimizer's final Jacobian (con.jac) instead
    # (for all-linear constraints both are identical, and ceq.JAC is kept
    # as it is computed with higher numerical accuracy)
    if (length(lavmodel@ceq.nonlinear.idx) > 0L &&
        nrow(lavmodel@con.jac) > 0L) {
      ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
      if (length(ceq_idx) > 0L) {
        jac <- lavmodel@con.jac[ceq_idx, , drop = FALSE]
      }
    }
    out <- lav_mat_ortho_complement(t(jac))
  } else {
    out <- NULL
  }
  out
}

# classify the inequality constraints at a converged solution: which rows
# are INACTIVE (not binding)? A constraint is only binding when it is
# satisfied with (near) equality AND its Lagrange multiplier exerts actual
# force on the objective gradient (strict complementarity). Both tests are
# needed: an estimate that merely lands within the slack of a bound, with a
# vanishing multiplier, is an interior solution -- treating the constraint
# as a binding equality would (wrongly) collapse its standard error to
# zero, with a discontinuous jump to the free standard error just across
# the slack threshold.
#
# The multiplier is compared on the scale of the force it exerts on the
# gradient, |lambda| * ||jac row||, which is invariant to the scaling of
# the constraint function (for c * g(x) >= 0, lambda scales with 1/c and
# the jacobian row with c). At an interior optimum the (post-hoc)
# multiplier is gradient-level noise, many orders of magnitude below the
# force of even a weakly binding constraint (which is proportional to how
# deep the unconstrained optimum violates the bound).
#
# arguments:
# - con0: the constraint values c(ceq, cin) at the solution
# - lambda: the Lagrange multipliers, aligned with con0
# - jac: the constraint Jacobian, rows aligned with con0
# - cin_flag: logical, TRUE for the inequality rows
lav_con_cin_inactive_idx <- function(con0 = NULL, lambda = NULL,
                                     jac = NULL, cin_flag = NULL,
                                     slack = 1e-05) {
  # measured post-hoc multiplier noise at interior optima is <= ~1e-7
  # (gradient-level), while the force of even a weakly binding constraint
  # is proportional to how deep the unconstrained optimum violates the
  # bound (~1e-4 for a violation of 1e-4); 1e-6 separates the two with a
  # comfortable margin on both sides
  force_tol <- 1e-06
  # clearly interior: strictly satisfied beyond the slack
  interior <- con0 > slack
  # no force: multiplier (in scale-free force units) vanishing; this also
  # covers a (wrong-sided) negative multiplier
  row_norm <- sqrt(rowSums(jac * jac))
  no_force <- (lambda * row_norm) < force_tol
  which(cin_flag & (interior | no_force))
}

lav_con_r2k <- function(lavmodel = NULL) {
  # constraint matrix
  m_r <- NULL
  if (!is.null(lavmodel)) {
    m_r <- lavmodel@ceq.JAC
  }
  stopifnot(!is.null(m_r))

  npar_full <- NCOL(m_r)
  # npar_red <-
  npar_full - NROW(m_r)

  m_k <- diag(npar_full)
  for (i in seq_len(NROW(m_r))) {
    idx1 <- which(m_r[i, ] == 1)
    idx2 <- which(m_r[i, ] == -1)
    m_k[idx2, idx1] <- 1
  }

  # remove redundant columns
  neg_idx <- which(colSums(m_r) < 0)
  m_k <- m_k[, -neg_idx]

  m_k
}

lav_con_lambda_pre <- function(lavobject = NULL, method = "Don") {
  # compute factor 'pre' so that pre %*% g = lambda
  method <- tolower(method)

  # note: drop = FALSE, so a single-row constraint matrix stays a matrix
  r <- lavobject@Model@con.jac[, , drop = FALSE]
  if (is.null(r) || length(r) == 0L) {
    return(numeric(0L))
  }

  info <- lavTech(lavobject, "information.first.order")
  npar <- nrow(info)

  # Don 1985
  if (method == "don") {
    r_plus <- MASS::ginv(r)

    # construct augmented matrix
    z <- rbind(
      cbind(info, t(r)),
      cbind(r, matrix(0, nrow = nrow(r), ncol = nrow(r)))
    )
    z_plus <- MASS::ginv(z)
    p_star <- z_plus[1:npar, 1:npar]
    pre <- t(r_plus) %*% (diag(npar) - info %*% p_star)

    # Bentler EQS manual
  } else if (method == "bentler") {
    info_inv <- solve(info)
    pre <- solve(r %*% info_inv %*% t(r)) %*% r %*% info_inv
  }

  pre
}
