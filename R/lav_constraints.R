lav_constraints_parse <- function(partable = NULL, constraints = NULL,
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

  no_explicit_constraints <- is.null(list_1) &&
    !any(partable$op %in% c(":=", "==", "<", ">")) &&
    !ceq_simple
  if (no_explicit_constraints) {
    def_function <- function() NULL
    ceq_function <- function() NULL
    cin_function <- function() NULL
    ceq_jacobian <- function() NULL
    cin_jacobian <- function() NULL
    ceq_jac <- matrix(0, nrow = 0L, ncol = npar)
    ceq_rhs <- numeric(0L)
    ceq_theta <- numeric(0L)
    ceq_linear_idx <- integer(0L)
    ceq_nonlinear_idx <- integer(0L)
    ceq_jac_null <- matrix(0, 0L, 0L)
    ceq_rhs_null <- numeric(0L)
    ceq_simple_k <- matrix(0, 0, 0)

    upper_idx <- if (is.null(partable$upper)) {
      integer(0L)
    } else {
      which(partable$free > 0L & is.finite(partable$upper))
    }
    lower_idx <- if (is.null(partable$lower)) {
      integer(0L)
    } else {
      which(partable$free > 0L & is.finite(partable$lower))
    }
    bound_idx <- c(upper_idx, lower_idx)
    if (length(bound_idx) > 0L) {
      free <- partable$free
      free[free > 0L] <- seq_along(free[free > 0L])
      cin_jac <- matrix(0, nrow = length(bound_idx), ncol = npar)
      upper_rows <- seq_along(upper_idx)
      lower_rows <- length(upper_idx) + seq_along(lower_idx)
      if (length(upper_idx) > 0L) {
        cin_jac[cbind(upper_rows, free[upper_idx])] <- -1
      }
      if (length(lower_idx) > 0L) {
        cin_jac[cbind(lower_rows, free[lower_idx])] <- 1
      }
      cin_rhs <- c(-partable$upper[upper_idx], partable$lower[lower_idx])
      cin_theta <- c(
        partable$upper[upper_idx] - theta[free[upper_idx]],
        theta[free[lower_idx]] - partable$lower[lower_idx]
      )
      attr(cin_rhs, "bound.idx") <- seq_along(bound_idx)
      attr(cin_theta, "bound.idx") <- seq_along(bound_idx)
      cin_linear_idx <- seq_along(bound_idx)
      cin_function <- local({
        upper_free <- free[upper_idx]
        lower_free <- free[lower_idx]
        upper <- partable$upper[upper_idx]
        lower <- partable$lower[lower_idx]
        function(.x., ...) {
          out <- c(upper - .x.[upper_free], .x.[lower_free] - lower)
          attr(out, "bound.idx") <- seq_along(out)
          out
        }
      })
    } else {
      cin_jac <- matrix(0, nrow = 0L, ncol = npar)
      cin_rhs <- numeric(0L)
      cin_theta <- numeric(0L)
      cin_linear_idx <- integer(0L)
    }

    return(list(
      def.function = def_function,
      ceq.function = ceq_function,
      ceq.JAC = ceq_jac,
      ceq.jacobian = ceq_jacobian,
      ceq.rhs = ceq_rhs,
      ceq.theta = ceq_theta,
      ceq.linear.idx = ceq_linear_idx,
      ceq.nonlinear.idx = ceq_nonlinear_idx,
      ceq.linear.flag = FALSE,
      ceq.nonlinear.flag = FALSE,
      ceq.flag = FALSE,
      ceq.linear.only.flag = FALSE,
      ceq.JAC.NULL = ceq_jac_null,
      ceq.rhs.NULL = ceq_rhs_null,
      ceq.simple.only = FALSE,
      ceq.simple.K = ceq_simple_k,
      cin.function = cin_function,
      cin.JAC = cin_jac,
      cin.jacobian = cin_jacobian,
      cin.rhs = cin_rhs,
      cin.theta = cin_theta,
      cin.linear.idx = cin_linear_idx,
      cin.nonlinear.idx = integer(0L),
      cin.linear.flag = length(cin_linear_idx) > 0L,
      cin.nonlinear.flag = FALSE,
      cin.flag = FALSE,
      cin.only.flag = FALSE,
      cin.simple.only = length(cin_linear_idx) > 0L
    ))
  }

  # variable definitions
  def_function <- lav_partable_constraints_def(partable,
    con = list_1,
    debug = debug
  )

  # construct ceq/ciq functions
  ceq_function <- lav_partable_constraints_ceq(partable,
    con = list_1,
    debug = debug
  )
  # linear or nonlinear?
  ceq_linear_idx <- lav_constraints_linear_idx(
    func = ceq_function,
    npar = npar
  )
  ceq_nonlinear_idx <- lav_constraints_nonlinear_idx(
    func = ceq_function,
    npar = npar
  )

  # inequalities
  cin_function <- lav_partable_constraints_ciq(partable,
    con = list_1,
    debug = debug
  )

  # linear or nonlinear?
  cin_linear_idx <- lav_constraints_linear_idx(
    func = cin_function,
    npar = npar
  )
  cin_nonlinear_idx <- lav_constraints_nonlinear_idx(
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

  # check for empty/unused constraints (new in 0.6-22)
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

lav_constraints_linear_idx <- function(func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) {
    return(integer(0L))
  }

  # seed 1: rnorm
  a0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  # seed 2: rnorm
  a1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  a0min_a1 <- a0 - a1
  linear <- apply(a0min_a1, 1, function(x) all(x == 0))
  which(linear)
}

lav_constraints_nonlinear_idx <- function(func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) {
    return(integer(0L))
  }

  # seed 1: rnorm
  a0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  # seed 2: rnorm
  a1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  a0min_a1 <- a0 - a1
  linear <- apply(a0min_a1, 1, function(x) all(x == 0))
  which(!linear)
}

# check if the equality constraints are 'simple' (a == b)
lav_constraints_check_simple <- function(lavmodel = NULL) {
  ones <- (lavmodel@ceq.JAC == 1 | lavmodel@ceq.JAC == -1)
  simple <- all(lavmodel@ceq.rhs == 0) &&
    all(apply(lavmodel@ceq.JAC != 0, 1, sum) == 2) &&
    all(apply(ones, 1, sum) == 2) &&
    length(lavmodel@ceq.nonlinear.idx) == 0

  # TRUE or FALSE
  simple
}

lav_constraints_r2k <- function(lavmodel = NULL) {
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

lav_constraints_lambda_pre <- function(lavobject = NULL, method = "Don") {
  # compute factor 'pre' so that pre %*% g = lambda
  method <- tolower(method)

  r <- lavobject@Model@con.jac[, ]
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
