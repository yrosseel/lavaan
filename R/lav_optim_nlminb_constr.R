# constrained optimization
# - references: * Nocedal & Wright (2006) Chapter 17
#               * Optimization with constraints by Madsen, Nielsen & Tingleff
#               * original papers: Powell, 1969 and Rockafeller, 1974
# - using 'nlminb' for the unconstrained subproblem
# - convergence scheme is based on the auglag function in the alabama package
nlminb.constr <- function(start, objective, gradient = NULL, hessian = NULL,
                          ..., scale = 1, control = list(),
                          lower = -Inf, upper = Inf,
                          ceq = NULL, ceq.jac = NULL,
                          cin = NULL, cin.jac = NULL,
                          control.outer = list()) {
  # we need a gradient
  stopifnot(!is.null(gradient))

  # if no 'ceq' or 'cin' function, we create a dummy one
  if (is.null(ceq)) {
    ceq <- function(x, ...) {
      return(numeric(0))
    }
  }
  if (is.null(cin)) {
    cin <- function(x, ...) {
      return(numeric(0))
    }
  }

  # if no user-supplied jacobian functions, create them
  if (is.null(ceq.jac)) {
    if (is.null(ceq)) {
      ceq.jac <- function(x, ...) {
        matrix(0, nrow = 0L, ncol = length(x))
      }
    } else {
      ceq.jac <- function(x, ...) {
        numDeriv::jacobian(func = ceq, x = x, ...)
      }
    }
  }
  if (is.null(cin.jac)) {
    if (is.null(cin)) {
      cin.jac <- function(x, ...) {
        matrix(0, nrow = 0L, ncol = length(x))
      }
    } else {
      cin.jac <- function(x, ...) {
        numDeriv::jacobian(func = cin, x = x, ...)
      }
    }
  }

  # how many ceq and cin constraints?
  nceq <- length(ceq(start))
  ncin <- length(cin(start))
  ncon <- nceq + ncin
  ceq.idx <- cin.idx <- integer(0)
  if (nceq > 0L) ceq.idx <- 1:nceq
  if (ncin > 0L) cin.idx <- nceq + 1:ncin
  cin.flag <- rep(FALSE, length(ncon))
  if (ncin > 0L) cin.flag[cin.idx] <- TRUE

  # control outer default values
  control.outer.default <- list(
    mu0 = 100,
    lambda0 = 10,
    tol = 1e-06, # changed this in 0.4-12
    itmax = 100L,
    verbose = FALSE
  )
  control.outer <- modifyList(control.outer.default, control.outer)


  # construct augmented lagrangian function
  auglag <- function(x, ...) {
    # apply constraints
    ceq0 <- ceq(x, ...)
    cin0 <- cin(x, ...)
    con0 <- c(ceq0, cin0)
    # 'release' inactive constraints
    if (ncin > 0L) {
      slack <- lambda / mu
      inactive.idx <- which(cin.flag & con0 > slack)
      con0[inactive.idx] <- slack[inactive.idx]
    }
    objective(x, ...) - sum(lambda * con0) + (mu / 2) * sum(con0 * con0)
  }

  fgrad <- function(x, ...) {
    # apply constraints
    ceq0 <- ceq(x, ...)
    cin0 <- cin(x, ...)
    con0 <- c(ceq0, cin0)
    # jacobian
    JAC <- rbind(ceq.jac(x, ...), cin.jac(x, ...))
    lambda.JAC <- lambda * JAC

    # handle inactive constraints
    if (ncin > 0L) {
      slack <- lambda / mu
      inactive.idx <- which(cin.flag & con0 > slack)
      if (length(inactive.idx) > 0L) {
        JAC <- JAC[-inactive.idx, , drop = FALSE]
        lambda.JAC <- lambda.JAC[-inactive.idx, , drop = FALSE]
        con0 <- con0[-inactive.idx]
      }
    }

    if (nrow(JAC) > 0L) {
      (gradient(x, ...) - colSums(lambda.JAC) +
        mu * as.numeric(t(JAC) %*% con0))
    } else {
      gradient(x, ...)
    }
  }


  # initialization
  ceq0 <- ceq(start, ...)
  cin0 <- cin(start, ...)
  con0 <- c(ceq0, cin0)
  lambda <- rep(control.outer$lambda0, length(con0))
  mu <- control.outer$mu0
  inactive.idx <- integer(0)
  if (ncin > 0L) {
    slack <- lambda / mu
    inactive.idx <- which(cin.flag & con0 > slack)
    con0[inactive.idx] <- slack[inactive.idx]
  }
  K <- max(abs(con0))
  if (control.outer$verbose) {
    cat("init cin0 values: ", cin0, "\n")
    cat("init ceq0 values: ", ceq0, "\n")
    cat("init slack values: ", lambda / mu, "\n")
    cat("init inactive idx: ", inactive.idx, "\n")
    cat("init con0 values: ", con0, "\n")
    cat("K = max con0: ", K, "\n")
  }

  r <- obj <- objective(start, ...)
  feval <- 0L
  geval <- 0L
  niter <- 0L
  ilack <- 0L
  Kprev <- K
  mu0 <- control.outer$mu0 / Kprev
  if (is.infinite(mu0)) mu0 <- 1.0
  mu <- mu0

  K <- Inf
  x.par <- start
  for (i in 1:control.outer$itmax) {
    x.old <- x.par
    r.old <- r
    ############################################################
    if (control.outer$verbose) {
      cat("\nStarting inner optimization [", i, "]:\n")
      cat("lambda: ", lambda, "\n")
      cat("mu: ", mu, "\n")
    }
    optim.out <- nlminb(
      start = x.par, objective = auglag,
      gradient = fgrad, control = control,
      lower = lower, upper = upper,
      scale = scale, ...
    )
    ############################################################
    x.par <- optim.out$par
    r <- optim.out$objective
    feval <- feval + optim.out$evaluations[1]
    geval <- geval + optim.out$evaluations[2]
    niter <- niter + optim.out$iterations

    # check constraints
    ceq0 <- ceq(x.par, ...)
    cin0 <- cin(x.par, ...)
    con0 <- c(ceq0, cin0)
    if (ncin > 0L) {
      slack <- lambda / mu
      inactive.idx <- which(cin.flag & con0 > slack)
      con0[inactive.idx] <- slack[inactive.idx]
    }
    K <- max(abs(con0))
    if (control.outer$verbose) {
      cat("cin0 values: ", cin0, "\n")
      cat("ceq0 values: ", ceq0, "\n")
      cat("active threshold: ", lambda / mu, "\n")
      cat("inactive idx: ", inactive.idx, "\n")
      cat("con0 values: ", con0, "\n")
      cat("K = max con0: ", K, " Kprev = ", Kprev, "\n")
    }

    # update K or mu (see Powell, 1969)
    if (K <= Kprev / 4) {
      lambda <- lambda - (mu * con0)
      Kprev <- K
    } else {
      mu <- 10 * mu
    }

    # check convergence
    pconv <- max(abs(x.par - x.old))
    if (pconv < control.outer$tol) {
      ilack <- ilack + 1L
    } else {
      ilack <- 0L
    }

    if ((is.finite(r) && is.finite(r.old) &&
      abs(r - r.old) < control.outer$tol && K < control.outer$tol) |
      ilack >= 3) {
      break
    }
  }

  # output
  a <- list()

  if (i == control.outer$itmax) {
    a$convergence <- 10L
    a$message <- "nlminb.constr ran out of iterations and did not converge"
  } else if (K > control.outer$tol) {
    a$convergence <- 11L
    a$message <- "Convergence due to lack of progress in parameter updates"
  } else {
    a$convergence <- 0L
    a$message <- "converged"
  }
  a$par <- optim.out$par
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
  JAC <- rbind(ceq.jac(a$par, ...), cin.jac(a$par, ...))
  inactive.idx <- integer(0L)
  cin.idx <- which(cin.flag)
  # ceq.idx <- which(!cin.flag)
  if (ncin > 0L) {
    # FIXME: slack value not too strict??
    slack <- 1e-05
    # cat("DEBUG:\n"); print(con0)
    inactive.idx <- which(cin.flag & con0 > slack)
    # if(length(inactive.idx) > 0L) {
    #    JAC        <-        JAC[-inactive.idx,,drop=FALSE]
    # }
  }
  attr(JAC, "inactive.idx") <- inactive.idx
  attr(JAC, "cin.idx") <- cin.idx
  attr(JAC, "ceq.idx") <- ceq.idx
  a$con.jac <- JAC

  a
}
