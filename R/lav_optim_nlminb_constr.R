# constrained optimization
# - references: * Nocedal & Wright (2006) Chapter 17
#               * Optimization with constraints by Madsen, Nielsen & Tingleff
#               * original papers: Powell, 1969 and Rockafeller, 1974
# - using 'nlminb' for the unconstrained subproblem
# - convergence scheme is based on the auglag function in the alabama package
nlminb_constr <- function(start, objective, gradient = NULL, hessian = NULL,
                          ..., scale = 1, control = list(),
                          lower = -Inf, upper = Inf,
                          ceq = NULL, ceq_jac = NULL,
                          cin = NULL, cin_jac = NULL,
                          control_outer = list()) {
  # we need a gradient
  stopifnot(!is.null(gradient))

  # if no 'ceq' or 'cin' function, we create a dummy one
  if (is.null(ceq)) {
    ceq <- function(x, ...) {
      numeric(0)
    }
  }
  if (is.null(cin)) {
    cin <- function(x, ...) {
      numeric(0)
    }
  }

  # if no user-supplied jacobian functions, create them
  if (is.null(ceq_jac)) {
    if (is.null(ceq)) {
      ceq_jac <- function(x, ...) {
        matrix(0, nrow = 0L, ncol = length(x))
      }
    } else {
      ceq_jac <- function(x, ...) {
        numDeriv::jacobian(func = ceq, x = x, ...)
      }
    }
  }
  if (is.null(cin_jac)) {
    if (is.null(cin)) {
      cin_jac <- function(x, ...) {
        matrix(0, nrow = 0L, ncol = length(x))
      }
    } else {
      cin_jac <- function(x, ...) {
        numDeriv::jacobian(func = cin, x = x, ...)
      }
    }
  }

  # how many ceq and cin constraints?
  nceq <- length(ceq(start))
  ncin <- length(cin(start))
  ncon <- nceq + ncin
  ceq_idx <- cin_idx <- integer(0)
  if (nceq > 0L) ceq_idx <- 1:nceq
  if (ncin > 0L) cin_idx <- nceq + 1:ncin
  cin_flag <- rep(FALSE, length(ncon))
  if (ncin > 0L) cin_flag[cin_idx] <- TRUE

  # control outer default values
  control_outer_default <- list(
    mu0 = 100,
    lambda0 = 10,
    tol = 1e-06, # changed this in 0.4-12
    itmax = 100L,
    verbose = FALSE
  )
  control_outer <- modifyList(control_outer_default, control_outer)


  # construct augmented lagrangian function
  auglag <- function(x, ...) {
    # apply constraints
    ceq0 <- ceq(x, ...)
    cin0 <- cin(x, ...)
    con0 <- c(ceq0, cin0)
    # 'release' inactive constraints
    if (ncin > 0L) {
      slack <- lambda / mu
      inactive_idx <- which(cin_flag & con0 > slack)
      con0[inactive_idx] <- slack[inactive_idx]
    }
    objective(x, ...) - sum(lambda * con0) + (mu / 2) * sum(con0 * con0)
  }

  fgrad <- function(x, ...) {
    # apply constraints
    ceq0 <- ceq(x, ...)
    cin0 <- cin(x, ...)
    con0 <- c(ceq0, cin0)
    # jacobian
    jac <- rbind(ceq_jac(x, ...), cin_jac(x, ...))
    lambda_jac <- lambda * jac

    # handle inactive constraints
    if (ncin > 0L) {
      slack <- lambda / mu
      inactive_idx <- which(cin_flag & con0 > slack)
      if (length(inactive_idx) > 0L) {
        jac <- jac[-inactive_idx, , drop = FALSE]
        lambda_jac <- lambda_jac[-inactive_idx, , drop = FALSE]
        con0 <- con0[-inactive_idx]
      }
    }

    if (nrow(jac) > 0L) {
      (gradient(x, ...) - colSums(lambda_jac) +
        mu * as.numeric(t(jac) %*% con0))
    } else {
      gradient(x, ...)
    }
  }


  # initialization
  ceq0 <- ceq(start, ...)
  cin0 <- cin(start, ...)
  con0 <- c(ceq0, cin0)
  lambda <- rep(control_outer$lambda0, length(con0))
  mu <- control_outer$mu0
  inactive_idx <- integer(0)
  if (ncin > 0L) {
    slack <- lambda / mu
    inactive_idx <- which(cin_flag & con0 > slack)
    con0[inactive_idx] <- slack[inactive_idx]
  }
  m_k <- max(abs(con0))
  if (control_outer$verbose) {
    cat("init cin0 values: ", cin0, "\n")
    cat("init ceq0 values: ", ceq0, "\n")
    cat("init slack values: ", lambda / mu, "\n")
    cat("init inactive idx: ", inactive_idx, "\n")
    cat("init con0 values: ", con0, "\n")
    cat("K = max con0: ", m_k, "\n")
  }

  r <- objective(start, ...) # r := obj
  feval <- 0L
  geval <- 0L
  niter <- 0L
  ilack <- 0L
  kprev <- m_k
  mu0 <- control_outer$mu0 / kprev
  if (is.infinite(mu0)) mu0 <- 1.0
  mu <- mu0

  m_k <- Inf
  x_par <- start
  for (i in 1:control_outer$itmax) {
    x_old <- x_par
    r_old <- r
    ############################################################
    if (control_outer$verbose) {
      cat("\nStarting inner optimization [", i, "]:\n")
      cat("lambda: ", lambda, "\n")
      cat("mu: ", mu, "\n")
    }
    optim_out <- nlminb(
      start = x_par, objective = auglag,
      gradient = fgrad, control = control,
      lower = lower, upper = upper,
      scale = scale, ...
    )
    ############################################################
    x_par <- optim_out$par
    r <- optim_out$objective
    feval <- feval + optim_out$evaluations[1]
    geval <- geval + optim_out$evaluations[2]
    niter <- niter + optim_out$iterations

    # check constraints
    ceq0 <- ceq(x_par, ...)
    cin0 <- cin(x_par, ...)
    con0 <- c(ceq0, cin0)
    if (ncin > 0L) {
      slack <- lambda / mu
      inactive_idx <- which(cin_flag & con0 > slack)
      con0[inactive_idx] <- slack[inactive_idx]
    }
    m_k <- max(abs(con0))
    if (control_outer$verbose) {
      cat("cin0 values: ", cin0, "\n")
      cat("ceq0 values: ", ceq0, "\n")
      cat("active threshold: ", lambda / mu, "\n")
      cat("inactive idx: ", inactive_idx, "\n")
      cat("con0 values: ", con0, "\n")
      cat("K = max con0: ", m_k, " Kprev = ", kprev, "\n")
    }

    # update K or mu (see Powell, 1969)
    if (m_k <= kprev / 4) {
      lambda <- lambda - (mu * con0)
      kprev <- m_k
    } else {
      mu <- 10 * mu
    }

    # check convergence
    pconv <- max(abs(x_par - x_old))
    if (pconv < control_outer$tol) {
      ilack <- ilack + 1L
    } else {
      ilack <- 0L
    }

    if ((is.finite(r) && is.finite(r_old) &&
      abs(r - r_old) < control_outer$tol && m_k < control_outer$tol) ||
      ilack >= 3) {
      break
    }
  }

  # output
  a <- list()

  if (i == control_outer$itmax) {
    a$convergence <- 10L
    a$message <- "nlminb_constr ran out of iterations and did not converge"
  } else if (m_k > control_outer$tol) {
    a$convergence <- 11L
    a$message <- "Convergence due to lack of progress in parameter updates"
  } else {
    a$convergence <- 0L
    a$message <- "converged"
  }
  a$par <- optim_out$par
  a$outer.iterations <- i
  a$lambda <- lambda
  a$mu <- mu
  # a$value <- objective(a$start, ...)
  # a$cin <- cin(a$start, ...)
  # a$ceq <- ceq(a$start, ...)
  a$evaluations <- c(feval, geval)
  a$iterations <- niter
  # a$kkt1 <- max(abs(a$fgrad)) <= 0.01 * (1 + abs(a$value))
  # a$kkt2 <- any(eigen(a$hessian)$value * control.optim$objectivescale> 0)

  # jacobian of ceq and 'active' cin
  ceq0 <- ceq(a$par, ...)
  cin0 <- cin(a$par, ...)
  con0 <- c(ceq0, cin0)
  jac <- rbind(ceq_jac(a$par, ...), cin_jac(a$par, ...))
  inactive_idx <- integer(0L)
  cin_idx <- which(cin_flag)
  # ceq.idx <- which(!cin.flag)
  if (ncin > 0L) {
    # FIXME: slack value not too strict??
    slack <- 1e-05
    # cat("DEBUG:\n"); print(con0)
    inactive_idx <- which(cin_flag & abs(con0) > slack)
    # if(length(inactive.idx) > 0L) {
    #    JAC        <-        JAC[-inactive.idx,,drop=FALSE]
    # }
  }
  attr(jac, "inactive.idx") <- inactive_idx
  attr(jac, "cin.idx") <- cin_idx
  attr(jac, "ceq.idx") <- ceq_idx
  a$con.jac <- jac

  a
}
