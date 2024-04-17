# Gauss-Newton style optimization
#
# Initial version needed for DLS - model based
# YR - 19 Jan 2021
#
# TODo:
# - what to do if the function value goes up?
# - handle general (nonlinear) equality constraints
# - handle general (nonlinear) inequality constraints
# - better approach for simple bounds
# ...

# YR - 04 Nov 2023: add huber = TRUE option to get 'outlier-robust' estimates
#      (see Yuan and Zhong 2008, where they call this IRLS_r)


# objective function, plus 'extra' information
# needed for a Gauss Newton step
lav_objective_GN <- function(x, lavsamplestats = NULL, lavmodel = NULL,
                             lavoptions = NULL, lavdata = NULL,
                             extra = FALSE, lambda = NULL) {
  # evaluate objective function
  lavmodel <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
  obj <- lav_model_objective(
    lavmodel = lavmodel, lavdata = lavdata,
    lavsamplestats = lavsamplestats
  )
  attributes(obj) <- NULL

  # monitoring obj only
  if (!extra) {
    # handle linear equality constraints
    if (lavmodel@eq.constraints) {
      hx <- lavmodel@ceq.function(x)
      obj <- obj + max(abs(lambda)) * sum(abs(hx))
    }
    return(list(obj = obj, U.invQ = NULL, lambda = lambda))
  }

  # model implied statistics
  lavimplied <- lav_model_implied(lavmodel = lavmodel)
  wls.est <- lav_model_wls_est(lavmodel = lavmodel, lavimplied = lavimplied)

  # observed statistics
  wls.obs <- lavsamplestats@WLS.obs

  # always use expected information
  A1 <- lav_model_h1_information_expected(
    lavobject = NULL,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavimplied = lavimplied,
    lavh1 = NULL,
    lavcache = NULL
  )
  # Delta
  Delta <- computeDelta(lavmodel = lavmodel)

  # first group
  g <- 1L
  if (lavmodel@estimator == "DWLS") {
    PRE.g <- t(Delta[[g]] * A1[[g]])
  } else {
    PRE.g <- t(Delta[[g]]) %*% A1[[g]]
  }
  Q.g <- PRE.g %*% (wls.obs[[g]] - wls.est[[g]])
  U.g <- PRE.g %*% Delta[[g]]

  # additional groups (if any)
  if (lavsamplestats@ngroups > 1L) {
    fg <- lavsamplestats@nobs[[1]] / lavsamplestats@ntotal
    Q <- fg * Q.g
    U <- fg * U.g

    for (g in 2:lavsamplestats@ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      if (lavmodel@estimator == "DWLS") {
        PRE.g <- t(Delta[[g]] * A1[[g]])
      } else {
        PRE.g <- t(Delta[[g]]) %*% A1[[g]]
      }
      Q.g <- PRE.g %*% (wls.obs[[g]] - wls.est[[g]])
      U.g <- PRE.g %*% Delta[[g]]

      Q <- Q + fg * Q.g
      U <- U + fg * U.g
    }
  } else {
    Q <- Q.g
    U <- U.g
  }

  # handle equality constraints
  # this can be made more efficient; see Jamshidian & Bentler 1993
  # where instead of inverting a p+r matrix, they use a p-r matrix
  # (if the eq constraints are linear)
  if (lavmodel@eq.constraints) {
    hx <- lavmodel@ceq.function(x)
    npar <- nrow(U)
    H <- lavmodel@con.jac
    U <- U + crossprod(H)
    U <- rbind(
      cbind(U, t(H)),
      cbind(H, matrix(0, nrow(H), nrow(H)))
    )
    Q <- rbind(Q, matrix(-hx, nrow(H), 1))
  }

  # compute step
  # note, we could use U + k*I for a given scalar 'k' (Levenberg, 1944)
  #                 or U + k*(diag(U) (Marquardt, 1963)
  U.invQ <- drop(solve(U, Q))
  if (lavmodel@eq.constraints) {
    # merit function
    lambda <- U.invQ[-seq_len(npar)]
    obj <- obj + max(abs(lambda)) * sum(abs(hx))
  } else {
    lambda <- NULL
  }

  list(obj = obj, U.invQ = U.invQ, lambda = lambda)
}

lav_optim_gn <- function(lavmodel = NULL, lavsamplestats = NULL,
                         lavpartable = NULL,
                         lavdata = NULL, lavoptions = NULL) {
  # no support (yet) for nonlinear constraints
  nonlinear.idx <- c(
    lavmodel@ceq.nonlinear.idx,
    lavmodel@cin.nonlinear.idx
  )
  if (length(nonlinear.idx) > 0L) {
    lav_msg_stop(gettext(
      "nonlinear constraints not supported (yet) with optim.method = \"GN\"."))
  }

  # no support (yet) for inequality constraints
  if (!is.null(body(lavmodel@cin.function))) {
    lav_msg_stop(gettext(
      "inequality constraints not supported (yet) with optim.method = \"GN\"."))
  }

  # extract current set of free parameters
  x <- lav_model_get_parameters(lavmodel)
  npar <- length(x)

  # extract bounds (if any)
  lb <- ub <- NULL
  if (!is.null(lavpartable) && !is.null(lavpartable$lower)) {
    lb <- lavpartable$lower[lavpartable$free > 0]
    stopifnot(length(x) == length(lb))
    lb.idx <- which(x < lb)
    if (length(lb.idx) > 0L) {
      x[lb.idx] <- lb[lb.idx]
    }
  }
  if (!is.null(lavpartable) && !is.null(lavpartable$upper)) {
    ub <- lavpartable$upper[lavpartable$free > 0]
    stopifnot(length(x) == length(ub))
    ub.idx <- which(x > ub)
    if (length(ub.idx) > 0L) {
      x[ub.idx] <- ub[ub.idx]
    }
  }

  # options
  verbose <- lavoptions$verbose
  iter.max <- lavoptions$optim.gn.iter.max
  tol.x <- lavoptions$optim.gn.tol.x
  stephalf.max <- as.integer(lavoptions$optim.gn.stephalf.max)
  if (stephalf.max < 0L) {
    stephalf.max <- 0L
  }

  # initialize
  iter <- 0
  alpha <- 1.0
  old.x <- x

  # start Gauss-Newton steps
  for (iter in seq_len(iter.max)) {
    old.out <- lav_objective_GN(
      x = old.x, lavsamplestats = lavsamplestats,
      lavoptions = lavoptions, lavdata = lavdata,
      lavmodel = lavmodel, extra = TRUE
    )
    old.obj <- old.out$obj
    U.invQ <- old.out$U.invQ

    # only the first time
    if (verbose && iter == 1L) {
      cat("iteration = ", sprintf("%2d", iter - 1L),
        ": objective = ", sprintf("%11.9f", old.obj), "\n",
        sep = ""
      )
    }

    # update
    alpha <- 1.0
    step <- U.invQ[seq_len(npar)]
    # TODO: if step-halving fails, we could also
    # allow the steps to be negative
    for (h in 1:max(1L, stephalf.max)) {
      new.x <- old.x + (alpha * step)

      # apply simple bounds (if any)
      if (!is.null(lb)) {
        lb.idx <- which(new.x < lb)
        if (length(lb.idx) > 0L) {
          new.x[lb.idx] <- lb[lb.idx]
        }
      }
      if (!is.null(ub)) {
        ub.idx <- which(new.x > ub)
        if (length(ub.idx) > 0L) {
          new.x[ub.idx] <- ub[ub.idx]
        }
      }

      new.obj <- lav_objective_GN(
        x = new.x,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata, lavoptions = lavoptions,
        lavmodel = lavmodel, extra = FALSE,
        lambda = old.out$lambda
      )$obj

      if (is.finite(new.obj) && new.obj < old.obj) {
        break
      } else if (stephalf.max == 0L) { # no step-halving!
        break
      } else {
        # step-halving
        alpha <- alpha / 2.0
        # if(verbose) {
        #    cat(" -- step halving -- : alpha = ", alpha, "\n")
        # }
      }
    }

    # TODO - if this fails, we need to recover somehow
    # negative steps:
    if (stephalf.max != 0L && h == stephalf.max) {
      if (verbose) {
        cat(" -- step halving failed; function value may increase.\n")
      }
      # forcing step with alpha = 1
      new.x <- old.x + (1 * step)
    }

    rms.x <- sqrt(mean((old.x - new.x) * (old.x - new.x)))

    # verbose?
    if (verbose) {
      cat("iteration = ", sprintf("%2d", iter),
        ": objective = ", sprintf("%11.9f", new.obj),
        " alpha = ", sprintf("%6.5f", alpha),
        " rms.x = ", sprintf("%9.9f", rms.x), "\n",
        sep = ""
      )
      # print(new.x)
    }

    # check for convergence
    if (rms.x < tol.x) {
      old.x <- new.x
      old.obj <- new.obj
      if (verbose) {
        cat("Gauss-Newton algorithm converged: rms.x = ",
          sprintf("%12.12f", rms.x), " < ",
          sprintf("%12.12f", tol.x), "\n",
          sep = ""
        )
      }
      break
    } else {
      old.x <- new.x
      old.obj <- new.obj
    }
  } # iter

  x <- new.x

  # one last evaluation, to get fx.group attribute
  lavmodel <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
  fx <- lav_model_objective(
    lavmodel = lavmodel, lavdata = lavdata,
    lavsamplestats = lavsamplestats
  )

  # add attributes
  if (iter < iter.max) {
    attr(x, "converged") <- TRUE
    attr(x, "warn.txt") <- ""
  } else {
    attr(x, "converged") <- FALSE
    attr(x, "warn.txt") <- paste("maxmimum number of iterations (",
      iter.max, ") ",
      "was reached without convergence.\n",
      sep = ""
    )
  }
  attr(x, "iterations") <- iter
  attr(x, "control") <- list(
    iter.max = iter.max,
    tol.x = tol.x
  )
  attr(x, "fx") <- fx

  x
}
