lav_constraints_parse <- function(partable = NULL, constraints = NULL,
                                  theta = NULL,
                                  debug = FALSE) {
  # just in case we do not have a $free column in partable
  if (is.null(partable$free)) {
    partable$free <- seq_len(length(partable$lhs))
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

  # number of free (but possibliy constrained) parameters
  npar <- length(theta)

  # parse the constraints
  if (is.null(constraints)) {
    LIST <- NULL
  } else if (!is.character(constraints)) {
    lav_msg_stop(gettext("constraints should be a string"))
  } else {
    FLAT <- lavParseModelString(constraints)
    CON <- attr(FLAT, "constraints")
    LIST <- list()
    if (length(CON) > 0L) {
      lhs <- unlist(lapply(CON, "[[", "lhs"))
      op <- unlist(lapply(CON, "[[", "op"))
      rhs <- unlist(lapply(CON, "[[", "rhs"))
      LIST$lhs <- c(LIST$lhs, lhs)
      LIST$op <- c(LIST$op, op)
      LIST$rhs <- c(LIST$rhs, rhs)
    } else {
      lav_msg_stop(gettext("no constraints found in constraints argument"))
    }
  }

  # simple equality constraints?
  ceq.simple <- FALSE
  if (!is.null(partable$unco)) {
    ceq.simple <- TRUE
  }

  # variable definitions
  def.function <- lav_partable_constraints_def(partable,
    con = LIST,
    debug = debug
  )

  # construct ceq/ciq functions
  ceq.function <- lav_partable_constraints_ceq(partable,
    con = LIST,
    debug = debug
  )
  # linear or nonlinear?
  ceq.linear.idx <- lav_constraints_linear_idx(
    func = ceq.function,
    npar = npar
  )
  ceq.nonlinear.idx <- lav_constraints_nonlinear_idx(
    func = ceq.function,
    npar = npar
  )

  # inequalities
  cin.function <- lav_partable_constraints_ciq(partable,
    con = LIST,
    debug = debug
  )

  # linear or nonlinear?
  cin.linear.idx <- lav_constraints_linear_idx(
    func = cin.function,
    npar = npar
  )
  cin.nonlinear.idx <- lav_constraints_nonlinear_idx(
    func = cin.function,
    npar = npar
  )

  # Jacobians
  if (!is.null(body(ceq.function))) {
    ceq.JAC <- try(lav_func_jacobian_complex(
      func = ceq.function,
      x = theta
    ), silent = TRUE)
    if (inherits(ceq.JAC, "try-error")) { # eg. pnorm()
      ceq.JAC <- lav_func_jacobian_simple(func = ceq.function, x = theta)
    }

    # constants
    # do we have a non-zero 'rhs' elements? FIXME!!! is this reliable??
    ceq.rhs <- -1 * ceq.function(numeric(npar))

    # evaluate constraints
    ceq.theta <- ceq.function(theta)
  } else {
    ceq.JAC <- matrix(0, nrow = 0L, ncol = npar)
    ceq.rhs <- numeric(0L)
    ceq.theta <- numeric(0L)
  }

  if (!is.null(body(cin.function))) {
    cin.JAC <- try(lav_func_jacobian_complex(
      func = cin.function,
      x = theta
    ), silent = TRUE)
    if (inherits(cin.JAC, "try-error")) { # eg. pnorm()
      cin.JAC <- lav_func_jacobian_simple(func = cin.function, x = theta)
    }

    # constants
    # do we have a non-zero 'rhs' elements? FIXME!!! is this reliable??
    cin.rhs <- -1 * cin.function(numeric(npar))

    # evaluate constraints
    cin.theta <- cin.function(theta)
  } else {
    cin.JAC <- matrix(0, nrow = 0L, ncol = npar)
    cin.rhs <- numeric(0L)
    cin.theta <- numeric(0L)
  }

  # shortcut flags
  ceq.linear.flag <- length(ceq.linear.idx) > 0L
  ceq.nonlinear.flag <- length(ceq.nonlinear.idx) > 0L
  ceq.flag <- ceq.linear.flag || ceq.nonlinear.flag

  cin.linear.flag <- length(cin.linear.idx) > 0L
  cin.nonlinear.flag <- length(cin.nonlinear.idx) > 0L
  cin.flag <- cin.linear.flag || cin.nonlinear.flag

  ceq.only.flag <- ceq.flag && !cin.flag
  cin.only.flag <- cin.flag && !ceq.flag

  ceq.linear.only.flag <- (ceq.linear.flag &&
    !ceq.nonlinear.flag &&
    !cin.flag)

  ceq.simple.only <- ceq.simple && !ceq.flag && !cin.flag

  # additional info if ceq.linear.flag
  if (ceq.linear.flag) {
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

    # compute range+null space of the jacobion (JAC) of the constraint
    # matrix
    # JAC <- lav_func_jacobian_complex(func = ceq.function,
    #           x = lavpartable$start[lavpartable$free > 0L]
    QR <- qr(t(ceq.JAC))
    ranK <- QR$rank
    Q <- qr.Q(QR, complete = TRUE)
    # Q1 <- Q[,1:ranK, drop = FALSE]         # range space
    # Q2 <- Q[,-seq_len(ranK), drop = FALSE] # null space
    # R <- qr.R(QR)
    ceq.JAC.NULL <- Q[, -seq_len(ranK), drop = FALSE]

    if (all(ceq.rhs == 0)) {
      ceq.rhs.NULL <- numeric(npar)
    } else {
      tmp <- qr.coef(QR, diag(npar))
      NA.idx <- which(is.na(rowSums(tmp))) # catch NAs
      if (length(NA.idx) > 0L) {
        tmp[NA.idx, ] <- 0
      }
      ceq.rhs.NULL <- as.numeric(ceq.rhs %*% tmp)
    }
  } else {
    ceq.JAC.NULL <- matrix(0, 0L, 0L)
    ceq.rhs.NULL <- numeric(0L)
  }

  # if simple equalities only, create 'K' matrix
  ceq.simple.K <- matrix(0, 0, 0)
  if (ceq.simple.only) {
    n.unco <- max(partable$unco)
    n.free <- max(partable$free)
    ceq.simple.K <- matrix(0, nrow = n.unco, ncol = n.free)
    #####
    #####     FIXME !
    #####
    idx.free <- partable$free[partable$free > 0]
    for (k in 1:n.unco) {
      c <- idx.free[k]
      ceq.simple.K[k, c] <- 1
    }
  }

  # dummy jacobian 'function'
  ceq.jacobian <- function() NULL
  cin.jacobian <- function() NULL


  OUT <- list(
    def.function = def.function,
    ceq.function = ceq.function,
    ceq.JAC = ceq.JAC,
    ceq.jacobian = ceq.jacobian,
    ceq.rhs = ceq.rhs,
    ceq.theta = ceq.theta,
    ceq.linear.idx = ceq.linear.idx,
    ceq.nonlinear.idx = ceq.nonlinear.idx,
    ceq.linear.flag = ceq.linear.flag,
    ceq.nonlinear.flag = ceq.nonlinear.flag,
    ceq.flag = ceq.flag,
    ceq.linear.only.flag = ceq.linear.only.flag,
    ceq.JAC.NULL = ceq.JAC.NULL,
    ceq.rhs.NULL = ceq.rhs.NULL,
    ceq.simple.only = ceq.simple.only,
    ceq.simple.K = ceq.simple.K,
    cin.function = cin.function,
    cin.JAC = cin.JAC,
    cin.jacobian = cin.jacobian,
    cin.rhs = cin.rhs,
    cin.theta = cin.theta,
    cin.linear.idx = cin.linear.idx,
    cin.nonlinear.idx = cin.nonlinear.idx,
    cin.linear.flag = cin.linear.flag,
    cin.nonlinear.flag = cin.nonlinear.flag,
    cin.flag = cin.flag,
    cin.only.flag = cin.only.flag
  )

  OUT
}

lav_constraints_linear_idx <- function(func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) {
    return(integer(0L))
  }

  # seed 1: rnorm
  A0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  # seed 2: rnorm
  A1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  A0minA1 <- A0 - A1
  linear <- apply(A0minA1, 1, function(x) all(x == 0))
  which(linear)
}

lav_constraints_nonlinear_idx <- function(func = NULL, npar = NULL) {
  if (is.null(func) || is.null(body(func))) {
    return(integer(0L))
  }

  # seed 1: rnorm
  A0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  # seed 2: rnorm
  A1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

  A0minA1 <- A0 - A1
  linear <- apply(A0minA1, 1, function(x) all(x == 0))
  which(!linear)
}


# FIXME: is there a more elegant/robust way to do this??
lav_constraints_check_linear <- function(model) {
  # seed 1: rnorm
  A.ceq <- A.cin <- matrix(0, model@nx.free, 0)
  if (!is.null(body(model@ceq.function))) {
    A.ceq <- t(lav_func_jacobian_complex(func = model@ceq.function, x = rnorm(model@nx.free)))
  }
  if (!is.null(body(model@cin.function))) {
    A.cin <- t(lav_func_jacobian_complex(func = model@cin.function, x = rnorm(model@nx.free)))
  }
  A0 <- cbind(A.ceq, A.cin)

  # seed 2: rnorm
  A.ceq <- A.cin <- matrix(0, model@nx.free, 0)
  if (!is.null(body(model@ceq.function))) {
    A.ceq <- t(lav_func_jacobian_complex(func = model@ceq.function, x = rnorm(model@nx.free)))
  }
  if (!is.null(body(model@cin.function))) {
    A.cin <- t(lav_func_jacobian_complex(func = model@cin.function, x = rnorm(model@nx.free)))
  }
  A1 <- cbind(A.ceq, A.cin)

  A0minA1 <- all.equal(A0, A1)
  if (is.logical(A0minA1) && A0minA1 == TRUE) {
    return(TRUE)
  } else {
    return(FALSE)
  }
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

lav_constraints_R2K <- function(lavmodel = NULL, R = NULL) {
  # constraint matrix
  if (!is.null(lavmodel)) {
    R <- lavmodel@ceq.JAC
  }
  stopifnot(!is.null(R))

  npar.full <- NCOL(R)
  npar.red <- npar.full - NROW(R)

  K <- diag(npar.full)
  for (i in 1:NROW(R)) {
    idx1 <- which(R[i, ] == 1)
    idx2 <- which(R[i, ] == -1)
    K[idx2, idx1] <- 1
  }

  # remove redundant columns
  neg.idx <- which(colSums(R) < 0)
  K <- K[, -neg.idx]

  K
}

lav_constraints_lambda_pre <- function(lavobject = NULL, method = "Don") {
  # compute factor 'pre' so that pre %*% g = lambda
  method <- tolower(method)

  R <- lavobject@Model@con.jac[, ]
  if (is.null(R) || length(R) == 0L) {
    return(numeric(0L))
  }

  INFO <- lavTech(lavobject, "information.first.order")
  npar <- nrow(INFO)

  # Don 1985
  if (method == "don") {
    R.plus <- MASS::ginv(R)

    # construct augmented matrix
    Z <- rbind(
      cbind(INFO, t(R)),
      cbind(R, matrix(0, nrow = nrow(R), ncol = nrow(R)))
    )
    Z.plus <- MASS::ginv(Z)
    P.star <- Z.plus[1:npar, 1:npar]
    PRE <- t(R.plus) %*% (diag(npar) - INFO %*% P.star)

    # Bentler EQS manual
  } else if (method == "bentler") {
    INFO.inv <- solve(INFO)
    PRE <- solve(R %*% INFO.inv %*% t(R)) %*% R %*% INFO.inv
  }

  PRE
}
