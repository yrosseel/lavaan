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
lav_objective_gn <- function(x, lavsamplestats = NULL, lavmodel = NULL,
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
  wls_est <- lav_model_wls_est(lavmodel = lavmodel, lavimplied = lavimplied)

  # observed statistics
  wls_obs <- lavsamplestats@WLS.obs

  # always use expected information
  a1 <- lav_model_h1_information_expected(
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
  delta <- lav_model_delta(lavmodel = lavmodel)

  # first group
  g <- 1L
  if (lavmodel@estimator == "DWLS") {
    pre_g <- t(delta[[g]] * a1[[g]])
  } else {
    pre_g <- t(delta[[g]]) %*% a1[[g]]
  }
  q_g <- pre_g %*% (wls_obs[[g]] - wls_est[[g]])
  u_g <- pre_g %*% delta[[g]]

  # additional groups (if any)
  if (lavsamplestats@ngroups > 1L) {
    fg <- lavsamplestats@nobs[[1]] / lavsamplestats@ntotal
    q_1 <- fg * q_g
    u <- fg * u_g

    for (g in 2:lavsamplestats@ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      if (lavmodel@estimator == "DWLS") {
        pre_g <- t(delta[[g]] * a1[[g]])
      } else {
        pre_g <- t(delta[[g]]) %*% a1[[g]]
      }
      q_g <- pre_g %*% (wls_obs[[g]] - wls_est[[g]])
      u_g <- pre_g %*% delta[[g]]

      q_1 <- q_1 + fg * q_g
      u <- u + fg * u_g
    }
  } else {
    q_1 <- q_g
    u <- u_g
  }

  # handle equality constraints
  # this can be made more efficient; see Jamshidian & Bentler 1993
  # where instead of inverting a p+r matrix, they use a p-r matrix
  # (if the eq constraints are linear)
  if (lavmodel@eq.constraints) {
    hx <- lavmodel@ceq.function(x)
    npar <- nrow(u)
    h <- lavmodel@con.jac
    u <- u + crossprod(h)
    u <- rbind(
      cbind(u, t(h)),
      cbind(h, matrix(0, nrow(h), nrow(h)))
    )
    q_1 <- rbind(q_1, matrix(-hx, nrow(h), 1))
  }

  # compute step
  # note, we could use U + k*I for a given scalar 'k' (Levenberg, 1944)
  #                 or U + k*(diag(U) (Marquardt, 1963)
  u_inv_q <- drop(solve(u, q_1))
  if (lavmodel@eq.constraints) {
    # merit function
    lambda <- u_inv_q[-seq_len(npar)]
    obj <- obj + max(abs(lambda)) * sum(abs(hx))
  } else {
    lambda <- NULL
  }

  list(obj = obj, U.invQ = u_inv_q, lambda = lambda)
}

lav_optim_gn <- function(lavmodel = NULL, lavsamplestats = NULL,
                         lavpartable = NULL,
                         lavdata = NULL, lavoptions = NULL) {
  # no support (yet) for nonlinear constraints
  nonlinear_idx <- c(
    lavmodel@ceq.nonlinear.idx,
    lavmodel@cin.nonlinear.idx
  )
  if (length(nonlinear_idx) > 0L) {
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
    lb_idx <- which(x < lb)
    if (length(lb_idx) > 0L) {
      x[lb_idx] <- lb[lb_idx]
    }
  }
  if (!is.null(lavpartable) && !is.null(lavpartable$upper)) {
    ub <- lavpartable$upper[lavpartable$free > 0]
    stopifnot(length(x) == length(ub))
    ub_idx <- which(x > ub)
    if (length(ub_idx) > 0L) {
      x[ub_idx] <- ub[ub_idx]
    }
  }

  # options
  iter_max <- lavoptions$optim.gn.iter.max
  tol_x <- lavoptions$optim.gn.tol.x
  stephalf_max <- as.integer(lavoptions$optim.gn.stephalf.max)
  if (stephalf_max < 0L) {
    stephalf_max <- 0L
  }

  # initialize
  iter <- 0
  alpha <- 1.0
  old_x <- x

  # start Gauss-Newton steps
  for (iter in seq_len(iter_max)) {
    old_out <- lav_objective_gn(
      x = old_x, lavsamplestats = lavsamplestats,
      lavoptions = lavoptions, lavdata = lavdata,
      lavmodel = lavmodel, extra = TRUE
    )
    old_obj <- old_out$obj
    u_inv_q <- old_out$U.invQ

    # only the first time
    if (lav_verbose() && iter == 1L) {
      cat("iteration = ", sprintf("%2d", iter - 1L),
        ": objective = ", sprintf("%11.9f", old_obj), "\n",
        sep = ""
      )
    }

    # update
    alpha <- 1.0
    step <- u_inv_q[seq_len(npar)]
    # TODO: if step-halving fails, we could also
    # allow the steps to be negative
    for (h in 1:max(1L, stephalf_max)) {
      new_x <- old_x + (alpha * step)

      # apply simple bounds (if any)
      if (!is.null(lb)) {
        lb_idx <- which(new_x < lb)
        if (length(lb_idx) > 0L) {
          new_x[lb_idx] <- lb[lb_idx]
        }
      }
      if (!is.null(ub)) {
        ub_idx <- which(new_x > ub)
        if (length(ub_idx) > 0L) {
          new_x[ub_idx] <- ub[ub_idx]
        }
      }

      new_obj <- lav_objective_gn(
        x = new_x,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata, lavoptions = lavoptions,
        lavmodel = lavmodel, extra = FALSE,
        lambda = old_out$lambda
      )$obj

      if (is.finite(new_obj) && new_obj < old_obj) {
        break
      } else if (stephalf_max == 0L) { # no step-halving!
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
    if (stephalf_max != 0L && h == stephalf_max) {
      if (lav_verbose()) {
        cat(" -- step halving failed; function value may increase.\n")
      }
      # forcing step with alpha = 1
      new_x <- old_x + (1 * step)
    }

    rms_x <- sqrt(mean((old_x - new_x) * (old_x - new_x)))

    # verbose?
    if (lav_verbose()) {
      cat("iteration = ", sprintf("%2d", iter),
        ": objective = ", sprintf("%11.9f", new_obj),
        " alpha = ", sprintf("%6.5f", alpha),
        " rms.x = ", sprintf("%9.9f", rms_x), "\n",
        sep = ""
      )
      # print(new.x)
    }

    # check for convergence
    if (rms_x < tol_x) {
      old_x <- new_x
      old_obj <- new_obj
      if (lav_verbose()) {
        cat("Gauss-Newton algorithm converged: rms.x = ",
          sprintf("%12.12f", rms_x), " < ",
          sprintf("%12.12f", tol_x), "\n",
          sep = ""
        )
      }
      break
    } else {
      old_x <- new_x
      old_obj <- new_obj
    }
  } # iter

  x <- new_x

  # one last evaluation, to get fx.group attribute
  lavmodel <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
  fx <- lav_model_objective(
    lavmodel = lavmodel, lavdata = lavdata,
    lavsamplestats = lavsamplestats
  )

  # add attributes
  if (iter < iter_max) {
    attr(x, "converged") <- TRUE
    attr(x, "warn.txt") <- ""
  } else {
    attr(x, "converged") <- FALSE
    attr(x, "warn.txt") <- paste("maximum number of iterations (",
      iter_max, ") ",
      "was reached without convergence.\n",
      sep = ""
    )
  }
  attr(x, "iterations") <- iter
  attr(x, "control") <- list(
    iter.max = iter_max,
    tol.x = tol_x
  )
  attr(x, "fx") <- fx

  x
}
