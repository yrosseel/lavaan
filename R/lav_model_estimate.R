# overview NLMINB (default) versus CONSTR (=constrained optimization)

#                         | cin.simple | nonlinear | no cin
# ----------------------------------------------------------
# eq.constraints (linear) | CONSTR     | CONSTR    | NLMINB
# ceq.nonlinear           | CONSTR     | CONSTR    | CONSTR
# ceq.simple              | NLMINB     | CONSTR    | NLMINB
# no ceq                  | NLMINB     | CONSTR    | NLMINB



# model estimation
lav_model_est <- function(lavmodel = NULL,
                               lavpartable = NULL, # for parscale = "stand"
                               lavh1 = NULL, # for multilevel + parsc
                               lavsamplestats = NULL,
                               lavdata = NULL,
                               lavoptions = NULL,
                               lavcache = list(),
                               start = "model",
                               do_fit = TRUE) {
  lavpartable <- lav_pt_set_cache(lavpartable)
  estimator <- lavoptions$estimator
  verbose <- lav_verbose()
  debug <- lav_debug()
  ngroups <- lavsamplestats@ngroups

  if (lavsamplestats@missing.flag || estimator == "PML" ||
      lavdata@nlevels > 1L) {
    group_weight <- FALSE
  } else {
    group_weight <- TRUE
  }

  # backwards compatibility < 0.6-11
  if (is.null(lavoptions$optim.partrace)) {
    lavoptions$optim.partrace <- FALSE
  }

  if (lavoptions$optim.partrace) {
    # fx + parameter values
    penv <- new.env()
    penv$PARTRACE <- matrix(NA, nrow = 0, ncol = lavmodel@nx.free + 1L)
  }

  # starting values (ignoring equality constraints)
  x_unpack <- lav_model_get_parameters(lavmodel)

  # override? use simple instead? (new in 0.6-7)
  if (start == "simple") {
    start_1 <- numeric(length(lavpartable$lhs))
    # set loadings to 0.7
    loadings_idx <- which(lavpartable$free > 0L &
      lavpartable$op == "=~")
    if (length(loadings_idx) > 0L) {
      start_1[loadings_idx] <- 0.7
    }
    # set (only) variances to 1
    var_idx <- which(lavpartable$free > 0L &
      lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs)
    if (length(var_idx) > 0L) {
      start_1[var_idx] <- 1
    }

    if (lavmodel@ceq.simple.only) {
      x_unpack <- start_1[lavpartable$free > 0L &
        !duplicated(lavpartable$free)]
    } else {
      x_unpack <- start_1[lavpartable$free > 0L]
    }

    # override? use random starting values instead? (new in 0.6-18)
  } else if (start == "random") {
    start_1 <- lav_pt_random(
      lavpartable = lavpartable,
      # needed if we still need to compute bounds:
      lavh1 = lavh1,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )

    if (lavmodel@ceq.simple.only) {
      x_unpack <- start_1[lavpartable$free > 0L &
        !duplicated(lavpartable$free)]
    } else {
      x_unpack <- start_1[lavpartable$free > 0L]
    }
  }

  # new in 0.7-1: multilevel models where a whole block (level) is
  # saturated: the h1 estimates for this block are already available;
  # temporarily fix these parameters at their h1 values, so they do not
  # enter the optimization again as free parameters (this often helps
  # convergence); after optimization, they are free parameters again
  h1_sat_idx <- integer(0L)
  h1_sat_ok <- lavoptions$optim.fix.saturated
  if (is.null(h1_sat_ok)) { # backwards compatibility
    h1_sat_ok <- TRUE
  }
  if (h1_sat_ok && lavdata@nlevels > 1L && estimator == "ML" &&
      # only if the optimizer honours box constraints (NLMINB family)
      length(lavmodel@ceq.nonlinear.idx) == 0L &&
      (lavmodel@cin.simple.only ||
       (length(lavmodel@cin.linear.idx) == 0L &&
        length(lavmodel@cin.nonlinear.idx) == 0L)) &&
      (is.null(lavoptions$optim.method) ||
       toupper(lavoptions$optim.method) %in%
         c("NLMINB", "NLMINB0", "NLMINB1", "L.BFGS.B"))) {
    tmp <- lav_model_est_h1_saturated(
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavdata = lavdata, lavh1 = lavh1
    )
    if (length(tmp$x.idx) > 0L) {
      h1_sat_idx <- tmp$x.idx
      x_unpack[h1_sat_idx] <- tmp$value
    }
  }

  # 1. parameter scaling (new in 0.6-2, rewritten in 0.7-2)
  #
  # optim.parscale = "standardized": the optimizer works with the scaled
  # (packed) parameters u = p * scale, where p is the packed parameter
  # vector and 'scale' holds the factors that transform each free
  # parameter to its standardized-solution metric (as if the data were
  # standardized); see lav_model_est_parscale(). The reparameterization
  # is applied AFTER the equality-constraint packing (below), and undone
  # BEFORE the objective, the gradient and the constraint functions are
  # evaluated, so the constraints are always imposed on the original
  # parameters. (Up to 0.7-1, the free parameters were scaled BEFORE
  # packing, which distorted all but a == b equality constraints.)

  # for < 0.6 compatibility
  if (is.null(lavoptions$optim.parscale)) {
    lavoptions$optim.parscale <- "none"
  }

  # scaling factors in the metric of the free parameters
  parscale <- rep(1.0, length(x_unpack))
  if (lavoptions$optim.parscale != "none") {
    parscale <- lav_model_est_parscale(
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavh1 = lavh1, lavoptions = lavoptions
    )
    if (length(parscale) != length(x_unpack)) {
      parscale <- rep(1.0, length(x_unpack))
    }

    # repair degenerate variance starts: some starting values are
    # raw-metric constants (e.g., the 0.05 default for latent variances);
    # when the marker has a large variance, such a start is (almost) zero
    # in the standardized metric, which puts the start on the boundary of
    # the positive-definite region and can derail the optimizer; give
    # those the standardized-world default start (0.05) instead
    if (lavmodel@ceq.simple.only) {
      keep <- lavpartable$free > 0L & !duplicated(lavpartable$free)
    } else {
      keep <- lavpartable$free > 0L
    }
    is_var <- (lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs)[keep]
    small_idx <- which(is_var & abs(x_unpack * parscale) < 0.01)
    small_idx <- setdiff(small_idx, h1_sat_idx)
    if (length(small_idx) > 0L) {
      x_unpack[small_idx] <- 0.05 / parscale[small_idx]
    }
  }
  if (debug) {
    cat("parscale = ", parscale, "\n")
  }

  # 2. pack (apply equality constraints)
  if (lavmodel@eq.constraints && ncol(lavmodel@eq.constraints.K) > 0L) {
    x_pack <- as.numeric((x_unpack - lavmodel@eq.constraints.k0) %*%
      lavmodel@eq.constraints.K)
  } else {
    x_pack <- x_unpack
  }

  # 3. scale (in the packed metric)
  scale <- lav_model_est_parscale_pack(
    parscale = parscale, lavmodel = lavmodel
  )
  if (length(scale) != length(x_pack)) {
    scale <- rep(1.0, length(x_pack))
  }
  scaling <- any(scale != 1.0)

  # final starting values for optimizer
  start_x <- x_pack * scale
  if (debug) {
    cat("start.x = ", start_x, "\n")
  }

  # user-specified bounds? (new in 0.6-2)
  if (is.null(lavpartable$lower)) {
    lower <- -Inf
  } else {
    if (lavmodel@ceq.simple.only) {
      free_idx <- which(lavpartable$free > 0L &
        !duplicated(lavpartable$free))
      lower <- lavpartable$lower[free_idx]
    } else if (lavmodel@eq.constraints) {
      # bounds have no effect any longer....
      # 0.6-19 -> we switch to constrained estimation
      #lav_msg_warn(gettext(
      #  "bounds have no effect in the presence of linear
      #          equality constraints"))
      lower <- -Inf
    } else {
      lower <- lavpartable$lower[lavpartable$free > 0L]
    }
  }
  if (is.null(lavpartable$upper)) {
    upper <- +Inf
  } else {
    if (lavmodel@ceq.simple.only) {
      free_idx <- which(lavpartable$free > 0L &
        !duplicated(lavpartable$free))
      upper <- lavpartable$upper[free_idx]
    } else if (lavmodel@eq.constraints) {
      # bounds have no effect any longer....
      if (is.null(lavpartable$lower)) {
        # bounds have no effect any longer....
        # 0.6-19 -> we switch to constrained estimation
        #lav_msg_warn(gettext(
        # "bounds have no effect in the presence of linear
        #          equality constraints"))
      }
      upper <- +Inf
    } else {
      upper <- lavpartable$upper[lavpartable$free > 0L]
    }
  }

  # the optimizer iterates in the u = p * scale metric, so the box
  # constraints must be transformed as well (scale > 0, so the bounds
  # keep their orientation)
  if (scaling) {
    if (length(lower) %in% c(1L, length(scale))) {
      lower <- lower * scale
    }
    if (length(upper) %in% c(1L, length(scale))) {
      upper <- upper * scale
    }
  }

  # check for inconsistent lower/upper bounds
  # this may happen if we have equality constraints; qr() may switch
  # the sign...
  bad_idx <- which(lower > upper)
  if (length(bad_idx) > 0L) {
    # switch
    # tmp <- lower[bad.idx]
    # lower[bad.idx] <- upper[bad.idx]
    # upper[bad.idx] <- tmp
    lower[bad_idx] <- -Inf
    upper[bad_idx] <- +Inf
  }

  # new in 0.7-1: fix the parameters of saturated blocks at their h1
  # values during optimization, by forcing lower == upper == start
  if (length(h1_sat_idx) > 0L) {
    lower <- rep(lower, length.out = length(start_x))
    upper <- rep(upper, length.out = length(start_x))
    lower[h1_sat_idx] <- upper[h1_sat_idx] <- start_x[h1_sat_idx]
  }

  # function to be minimized
  objective_function <- function(x, verbose = FALSE, inf_to_max = FALSE,
                                 debug = FALSE) {
    # 3. standard deviations to variances
    # WARNING: x is still packed here!
    # if(lavoptions$optim.var.transform == "sqrt" &&
    #   length(lavmodel@x.free.var.idx) > 0L) {
    #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
    #    x.var <- x[lavmodel@x.free.var.idx]
    #    x.var.sign <- sign(x.var)
    #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
    # }

    # 3. unscale
    x <- x / scale

    # 2. unpack
    if (lavmodel@eq.constraints) {
      x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
        lavmodel@eq.constraints.k0
    }

    # update GLIST (change `state') and make a COPY!
    glist <- lav_model_x2glist(lavmodel, x = x)

    fx <- lav_model_objective(
      lavmodel = lavmodel,
      glist = glist,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache
    )

    # only for PML: divide by N (to speed up convergence)
    if (estimator == "PML") {
      fx <- fx / lavsamplestats@ntotal
    }



    if (debug || verbose) {
      cat("  objective function  = ",
        sprintf("%18.16f", fx), "\n",
        sep = ""
      )
    }
    if (debug) {
      # cat("Current unconstrained parameter values =\n")
      # tmp.x <- lav_model_get_parameters(lavmodel, glist=GLIST, type="unco")
      # print(tmp.x); cat("\n")
      cat("Current free parameter values =\n")
      print(x)
      cat("\n")
    }

    if (lavoptions$optim.partrace) {
      penv$PARTRACE <- rbind(penv$PARTRACE, c(fx, x))
    }

    # for L-BFGS-B
    # if(infToMax && is.infinite(fx)) fx <- 1e20
    if (!is.finite(fx)) {
      fx_group <- attr(fx, "fx.group")
      fx <- 1e20
      attr(fx, "fx.group") <- fx_group # only for lav_model_fit()
    }

    fx
  }

  gradient_function <- function(x, verbose = FALSE, inf_to_max = FALSE,
                                debug = FALSE) {
    # transform variances back
    # if(lavoptions$optim.var.transform == "sqrt" &&
    #   length(lavmodel@x.free.var.idx) > 0L) {
    #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
    #    x.var <- x[lavmodel@x.free.var.idx]
    #    x.var.sign <- sign(x.var)
    #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
    # }

    # 3. unscale
    x <- x / scale

    # 2. unpack
    if (lavmodel@eq.constraints) {
      x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
        lavmodel@eq.constraints.k0
    }

    # update GLIST (change `state') and make a COPY!
    glist <- lav_model_x2glist(lavmodel, x = x)

    dx <- lav_model_grad(
      lavmodel = lavmodel,
      glist = glist,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache,
      type = "free",
      group_weight = group_weight, ### check me!!
      ceq_simple = lavmodel@ceq.simple.only
    )

    if (debug) {
      cat("Gradient function (analytical) =\n")
      print(dx)
      cat("\n")
    }

    # 2. pack
    if (lavmodel@eq.constraints) {
      dx <- as.numeric(dx %*% lavmodel@eq.constraints.K)
    }

    # 3. scale (note: divide, not multiply!)
    dx <- dx / scale

    # 3. transform variances back
    # if(lavoptions$optim.var.transform == "sqrt" &&
    #   length(lavmodel@x.free.var.idx) > 0L) {
    #    x.var <- x[lavmodel@x.free.var.idx] # here in 'var' metric
    #    x.var.sign <- sign(x.var)
    #    x.var <- abs(x.var)
    #    x.sd <- sqrt(x.var)
    #    dx[lavmodel@x.free.var.idx] <-
    #        ( 2 * x.var.sign * dx[lavmodel@x.free.var.idx] * x.sd )
    # }

    # only for PML: divide by N (to speed up convergence)
    if (estimator == "PML") {
      dx <- dx / lavsamplestats@ntotal
    }

    if (debug) {
      cat("Gradient function (analytical, after eq.constraints.K) =\n")
      print(dx)
      cat("\n")
    }

    dx
  }

  gradient_function_numerical <- function(x, verbose = FALSE, debug = FALSE) {
    # NOTE: no need to 'transform' anything here (var/eq)
    # this is done anyway in objective_function

    # numerical approximation using the Richardson method
    npar <- length(x)
    h <- 10e-6
    dx <- numeric(npar)

    ## FIXME: call lav_model_objective directly!!
    for (i in 1:npar) {
      x_left <- x_left2 <- x_right <- x_right2 <- x
      x_left[i] <- x[i] - h
      x_left2[i] <- x[i] - 2 * h
      x_right[i] <- x[i] + h
      x_right2[i] <- x[i] + 2 * h
      fx_left <- objective_function(x_left, verbose = FALSE, debug = FALSE)
      fx_left2 <- objective_function(x_left2, verbose = FALSE, debug = FALSE)
      fx_right <- objective_function(x_right, verbose = FALSE, debug = FALSE)
      fx_right2 <- objective_function(x_right2, verbose = FALSE, debug = FALSE)
      dx[i] <- (fx_left2 - 8 * fx_left + 8 * fx_right - fx_right2) / (12 * h)
    }

    # dx <- lavGradientC(func=objective_function, x=x)
    # does not work if pnorm is involved... (eg PML)

    if (debug) {
      cat("Gradient function (numerical) =\n")
      print(dx)
      cat("\n")
    }

    dx
  }

  gradient_function_numerical_complex <- function(x, verbose = FALSE, debug = FALSE) { # nolint
    dx <- Re(lav_func_grad_complex(
      func = objective_function, x = x,
      h = sqrt(.Machine$double.eps)
    ))
    # does not work if pnorm is involved... (eg PML)

    if (debug) {
      cat("Gradient function (numerical complex) =\n")
      print(dx)
      cat("\n")
    }

    dx
  }


  # check if the initial values produce a positive definite Sigma
  # to begin with -- but only for estimator="ML"
  if (estimator %in% c("ML", "FML", "MML")) {
    sigma_hat <- lav_model_sigma(lavmodel, extra = TRUE)
    for (g in 1:ngroups) {
      if (!attr(sigma_hat[[g]], "po")) {
        group_txt <-
          if (ngroups > 1) gettextf(" in group %s.", g) else "."
        if (debug) {
          print(sigma_hat[[g]][, ])
        }
        lav_msg_warn(gettext(
          "initial model-implied matrix (Sigma) is not positive definite;
          check your model and/or starting parameters"), group_txt)
        x <- start_x
        fx <- as.numeric(NA)
        attr(fx, "fx.group") <- rep(as.numeric(NA), ngroups)
        attr(x, "converged") <- FALSE
        attr(x, "iterations") <- 0L
        attr(x, "control") <- lavoptions$control
        attr(x, "fx") <- fx
        return(x)
      }
    }
  }


  # nlminb's own scale heuristic (unchanged since 0.5): parameters that
  # start above 1 (in absolute value) are scaled by 1/|start|
  scale_1 <- rep(1.0, length(start_x))
  idx <- which(abs(start_x) > 1.0)
  if (length(idx) > 0L) {
    scale_1[idx] <- abs(1.0 / start_x[idx])
  }
  if (debug) {
    cat("SCALE = ", scale_1, "\n")
  }


  # first try: check if starting values return a finite value
  fx <- objective_function(start_x, verbose = verbose, debug = debug)
  if (!is.finite(fx)) {
    # emergency change of start.x
    if (length(h1_sat_idx) > 0L) {
      # do not touch the (fixed) h1 values of saturated blocks
      start_x[-h1_sat_idx] <- start_x[-h1_sat_idx] / 10
    } else {
      start_x <- start_x / 10
    }
  }



  # first some nelder mead steps? (default = FALSE)
  init_nelder_mead <- lavoptions$optim.init_nelder_mead

  # gradient: analytic, numerical or NULL?
  if (is.character(lavoptions$optim.gradient)) {
    if (lavoptions$optim.gradient %in% c("analytic", "analytical")) {
      gradient <- gradient_function
    } else if (lavoptions$optim.gradient %in% c("numerical", "numeric")) {
      gradient <- gradient_function_numerical
    } else if (lavoptions$optim.gradient %in% c("numeric.complex", "complex")) {
      gradient <- gradient_function_numerical_complex
    } else if (lavoptions$optim.gradient %in% c("NULL", "null")) {
      gradient <- NULL
    } else {
      lav_msg_warn(gettext("gradient should be analytic, numerical or NULL"))
    }
  } else if (is.logical(lavoptions$optim.gradient)) {
    if (lavoptions$optim.gradient) {
      gradient <- gradient_function
    } else {
      gradient <- NULL
    }
  } else if (is.null(lavoptions$optim.gradient)) {
    gradient <- gradient_function
  }


  # default optimizer
  if (length(lavmodel@ceq.nonlinear.idx) == 0L &&
      (lavmodel@cin.simple.only || (length(lavmodel@cin.linear.idx)    == 0L &&
                                    length(lavmodel@cin.nonlinear.idx) == 0L))
     ) {
    if (is.null(lavoptions$optim.method)) {
      optimizer <- "NLMINB"
      # OPTIMIZER <- "BFGS"  # slightly slower, no bounds; better scaling!
      # OPTIMIZER <- "L-BFGS-B"  # trouble with Inf values for fx!
    } else {
      optimizer <- toupper(lavoptions$optim.method)
      stopifnot(optimizer %in% c(
        "NLMINB0", "NLMINB1", "NLMINB2",
        "NLMINB", "BFGS", "L.BFGS.B", "NONE"
      ))
      if (optimizer == "NLMINB1") {
        optimizer <- "NLMINB"
      }
    }
  } else {
    if (is.null(lavoptions$optim.method)) {
      optimizer <- "NLMINB.CONSTR"
    } else {
      optimizer <- toupper(lavoptions$optim.method)
      stopifnot(optimizer %in% c("NLMINB.CONSTR", "NLMINB", "NONE"))
    }
    if (optimizer == "NLMINB") {
      optimizer <- "NLMINB.CONSTR"
    }
  }

  if (init_nelder_mead) {
    if (verbose) cat("  initial Nelder-Mead step:\n")
    # trace <- 0L
    # if (verbose) trace <- 1L
    optim_out <- optim(
      par = start_x,
      fn = objective_function,
      method = "Nelder-Mead",
      # control=list(maxit=10L,
      #             parscale=SCALE,
      #             trace=trace),
      hessian = FALSE,
      verbose = verbose, debug = debug
    )
    cat("\n")
    start_x <- optim_out$par
  }



  if (optimizer == "NLMINB0") {
    if (verbose)
      cat("  quasi-Newton steps using NLMINB0 (no analytic gradient):\n")
    # if(debug) control$trace <- 1L;
    control_nlminb <- list(
      eval.max = 20000L,
      iter.max = 10000L,
      trace = 0L,
      # abs.tol=1e-20, ### important!! fx never negative
      abs.tol = (.Machine$double.eps * 10),
      rel.tol = 1e-10,
      # step.min=2.2e-14, # in =< 0.5-12
      step.min = 1.0, # 1.0 in < 0.5-21
      step.max = 1.0,
      x.tol = 1.5e-8,
      xf.tol = 2.2e-14
    )
    control_nlminb <- modifyList(control_nlminb, lavoptions$control)
    control <- control_nlminb[c(
      "eval.max", "iter.max", "trace",
      "step.min", "step.max",
      "abs.tol", "rel.tol", "x.tol", "xf.tol"
    )]
    # cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
    optim_out <- nlminb(
      start = start_x,
      objective = objective_function,
      gradient = NULL,
      lower = lower,
      upper = upper,
      control = control,
      scale = scale_1,
      verbose = verbose, debug = debug
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  nlminb message says: ", optim_out$message, "\n")
      cat("  number of iterations: ", optim_out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$evaluations, "\n"
      )
    }

    # try again
    if (optim_out$convergence != 0L) {
      optim_out <- nlminb(
        start = start_x,
        objective = objective_function,
        gradient = NULL,
        lower = lower,
        upper = upper,
        control = control,
        scale = scale_1,
        verbose = verbose, debug = debug
      )
    }

    iterations <- optim_out$iterations
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "NLMINB") {
    if (verbose) cat("  quasi-Newton steps using NLMINB:\n")
    # if(debug) control$trace <- 1L;
    control_nlminb <- list(
      eval.max = 20000L,
      iter.max = 10000L,
      trace = 0L,
      # abs.tol=1e-20, ### important!! fx never negative
      abs.tol = (.Machine$double.eps * 10),
      rel.tol = 1e-10,
      # step.min=2.2e-14, # in =< 0.5-12
      step.min = 1.0, # 1.0 in < 0.5-21
      step.max = 1.0,
      x.tol = 1.5e-8,
      xf.tol = 2.2e-14
    )
    control_nlminb <- modifyList(control_nlminb, lavoptions$control)
    control <- control_nlminb[c(
      "eval.max", "iter.max", "trace",
      "step.min", "step.max",
      "abs.tol", "rel.tol", "x.tol", "xf.tol"
    )]
    # cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
    optim_out <- nlminb(
      start = start_x,
      objective = objective_function,
      gradient = gradient,
      lower = lower,
      upper = upper,
      control = control,
      scale = scale_1,
      verbose = verbose, debug = debug
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  nlminb message says: ", optim_out$message, "\n")
      cat("  number of iterations: ", optim_out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$evaluations, "\n"
      )
    }

    iterations <- optim_out$iterations
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "BFGS") {
    # warning: Bollen example with estimator=GLS does NOT converge!
    # (but WLS works!)
    # - BB.ML works too

    # note: optim()'s parscale is the typical SIZE of a parameter
    # (optimization is performed on par/parscale), the inverse of
    # nlminb's scale (up to 0.7-1, scale_1 itself was passed here)
    control_bfgs <- list(
      trace = 0L, fnscale = 1,
      parscale = 1 / scale_1,
      ndeps = 1e-3,
      maxit = 10000,
      abstol = 1e-20,
      reltol = 1e-10,
      REPORT = 1L
    )
    control_bfgs <- modifyList(control_bfgs, lavoptions$control)
    control <- control_bfgs[c(
      "trace", "fnscale", "parscale", "ndeps",
      "maxit", "abstol", "reltol", "REPORT"
    )]
    # trace <- 0L; if(verbose) trace <- 1L
    optim_out <- optim(
      par = start_x,
      fn = objective_function,
      gr = gradient,
      method = "BFGS",
      control = control,
      hessian = FALSE,
      verbose = verbose, debug = debug
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  optim BFGS message says: ", optim_out$message, "\n")
      # cat("number of iterations: ", optim.out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$counts, "\n"
      )
    }

    # iterations <- optim.out$iterations
    iterations <- optim_out$counts[1]
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "L.BFGS.B") {
    # warning, does not cope with Inf values!!

    control_lbfgsb <- list(
      trace = 0L, fnscale = 1,
      parscale = 1 / scale_1, # see BFGS above
      ndeps = 1e-3,
      maxit = 10000,
      REPORT = 1L,
      lmm = 5L,
      factr = 1e7,
      pgtol = 0
    )
    control_lbfgsb <- modifyList(control_lbfgsb, lavoptions$control)
    control <- control_lbfgsb[c(
      "trace", "fnscale", "parscale",
      "ndeps", "maxit", "REPORT", "lmm",
      "factr", "pgtol"
    )]
    optim_out <- optim(
      par = start_x,
      fn = objective_function,
      gr = gradient,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = control,
      hessian = FALSE,
      verbose = verbose, debug = debug,
      inf_to_max = TRUE
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  optim L-BFGS-B message says: ", optim_out$message, "\n")
      # cat("number of iterations: ", optim.out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$counts, "\n"
      )
    }

    # iterations <- optim.out$iterations
    iterations <- optim_out$counts[1]
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "NLMINB.CONSTR") {
    ocontrol <- list(verbose = verbose)
    if (!is.null(lavoptions$control$control.outer)) {
      ocontrol <- c(lavoptions$control$control.outer, verbose = verbose)
    }
    control_nlminb <- list(
      eval.max = 20000L,
      iter.max = 10000L,
      trace = 0L,
      # abs.tol=1e-20,
      abs.tol = (.Machine$double.eps * 10),
      rel.tol = 1e-9, # 1e-10 seems 'too strict'
      step.min = 1.0, # 1.0 in < 0.5-21
      step.max = 1.0,
      x.tol = 1.5e-8,
      xf.tol = 2.2e-14
    )
    control_nlminb <- modifyList(control_nlminb, lavoptions$control)
    control <- control_nlminb[c(
      "eval.max", "iter.max", "trace",
      "abs.tol", "rel.tol"
    )]
    cin <- cin_jac <- ceq <- ceq_jac <- NULL
    if (!is.null(body(lavmodel@cin.function))) cin <- lavmodel@cin.function
    if (!is.null(body(lavmodel@cin.jacobian))) cin_jac <- lavmodel@cin.jacobian
    if (!is.null(body(lavmodel@ceq.function))) ceq <- lavmodel@ceq.function
    if (!is.null(body(lavmodel@ceq.jacobian))) ceq_jac <- lavmodel@ceq.jacobian
    # parameter scaling: the constraint functions (and their jacobians)
    # are defined for the original parameters (there is no
    # equality-constraint packing in this branch), while the optimizer
    # works with u = x * scale
    if (scaling) {
      scale_con <- function(fun) {
        if (is.null(fun)) {
          return(NULL)
        }
        function(x, ...) fun(x / scale, ...)
      }
      scale_con_jac <- function(fun) {
        if (is.null(fun)) {
          return(NULL)
        }
        function(x, ...) {
          jac <- fun(x / scale, ...)
          if (!is.matrix(jac)) {
            jac <- matrix(jac, ncol = length(scale))
          }
          sweep(jac, 2L, scale, "/")
        }
      }
      cin <- scale_con(cin)
      ceq <- scale_con(ceq)
      cin_jac <- scale_con_jac(cin_jac)
      ceq_jac <- scale_con_jac(ceq_jac)
    }
    trace <- FALSE
    if (verbose) trace <- TRUE
    optim_out <- nlminb_constr(
      start = start_x,
      objective = objective_function,
      gradient = gradient,
      control = control,
      scale = scale_1,
      verbose = verbose, debug = debug,
      lower = lower,
      upper = upper,
      cin = cin, cin_jac = cin_jac,
      ceq = ceq, ceq_jac = ceq_jac,
      control_outer = ocontrol
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  nlminb_constr message says: ", optim_out$message, "\n")
      cat("  number of outer iterations: ", optim_out$outer.iterations, "\n")
      cat("  number of inner iterations: ", optim_out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$evaluations, "\n"
      )
    }

    iterations <- optim_out$iterations
    x <- optim_out$par
    if (optim_out$convergence == 0) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
    # the jacobian of the constraints back in the original metric
    if (scaling && is.matrix(optim_out$con.jac)) {
      jac_attr <- attributes(optim_out$con.jac)
      optim_out$con.jac <- sweep(optim_out$con.jac, 2L, scale, "*")
      attributes(optim_out$con.jac) <- jac_attr
    }
  } else if (optimizer == "NONE") {
    x <- start_x
    iterations <- 0L
    converged <- TRUE
    control <- list()

    # if inequality constraints, add con.jac/lambda
    # needed for df!
    if (length(lavmodel@ceq.nonlinear.idx) == 0L &&
        (lavmodel@cin.simple.only ||
         (length(lavmodel@cin.linear.idx) == 0L &&
          length(lavmodel@cin.nonlinear.idx) == 0L))) {
      optim_out <- list()
    } else {
      # if inequality constraints, add con.jac/lambda
      # needed for df!

      optim_out <- list()
      x_con <- start_x / scale
      if (is.null(body(lavmodel@ceq.function))) {
        ceq <- function(x, ...) {
          numeric(0)
        }
      } else {
        ceq <- lavmodel@ceq.function
      }
      if (is.null(body(lavmodel@cin.function))) {
        cin <- function(x, ...) {
          numeric(0)
        }
      } else {
        cin <- lavmodel@cin.function
      }
      ceq0 <- ceq(x_con)
      cin0 <- cin(x_con)
      con0 <- c(ceq0, cin0)
      jac <- rbind(
        numDeriv::jacobian(ceq, x = x_con),
        numDeriv::jacobian(cin, x = x_con)
      )
      nceq <- length(ceq(x_con))
      ncin <- length(cin(x_con))
      ncon <- nceq + ncin
      ceq_idx <- cin_idx <- integer(0)
      if (nceq > 0L) ceq_idx <- 1:nceq
      if (ncin > 0L) cin_idx <- nceq + 1:ncin
      cin_flag <- rep(FALSE, ncon)
      if (ncin > 0L) cin_flag[cin_idx] <- TRUE

      inactive_idx <- integer(0L)
      cin_idx <- which(cin_flag)
      if (ncin > 0L) {
        # optimizer == "NONE": no multiplier information is available for
        # the externally supplied solution, so a constraint that touches
        # the boundary is conservatively treated as binding (this only
        # feeds the df bookkeeping); no strict-complementarity test here
        slack <- 1e-05
        inactive_idx <- which(cin_flag & con0 > slack)
      }
      attr(jac, "inactive.idx") <- inactive_idx
      attr(jac, "cin.idx") <- cin_idx
      attr(jac, "ceq.idx") <- ceq_idx

      optim_out$con.jac <- jac
      optim_out$lambda <- rep(0, ncon)
    }
  }

  # new in 0.6-19
  # if NLMINB() + cin.simple.only, add con.jac and lambda to optim.out
  if (optimizer %in% c("NLMINB", "NLMINB0", "L.BFGS.B") &&
      (lavmodel@cin.simple.only || lavmodel@ceq.simple.only)) {

    if (lavmodel@cin.simple.only && nrow(lavmodel@cin.JAC) > 0L) {
      # JAC
      cin_jac_1 <- lavmodel@cin.JAC

      # lambda (post-hoc): the bound rows involve one parameter each, so
      # the rows of cin.JAC are mutually orthogonal and the least-squares
      # multiplier reduces to the plain inner product with the gradient.
      # NOTE: for ceq.simple models both the parameter vector and the
      # gradient live in the reduced (packed) space, while cin.function
      # and cin.JAC operate on the full (unco) space -- unpack BOTH
      # (evaluating cin.function on the packed x returned garbage
      # constraint values, silently corrupting the old classification)
      # back to the packed metric (the optimizer works with u = p * scale)
      dx <- gradient(x) * scale
      x_p <- x / scale
      if (lavmodel@ceq.simple.only) {
        unpack_idx <- lavpartable$free[lavpartable$free > 0]
        x_unpack <- x_p[unpack_idx]
        dx_unpack <- dx[unpack_idx]
      } else if (lavmodel@eq.constraints) {
        # unreachable via the standard pipeline: eq.constraints is a
        # packing flag that is only TRUE when equality constraints are the
        # ONLY constraints, which contradicts cin.simple.only; kept as a
        # safety net
        x_unpack <- as.numeric(lavmodel@eq.constraints.K %*% x_p) +
          lavmodel@eq.constraints.k0
        dx_unpack <- as.numeric(lavmodel@eq.constraints.K %*% dx)
      } else {
        x_unpack <- x_p
        dx_unpack <- dx
      }
      con0 <- lavmodel@cin.function(x_unpack)
      cin_lambda <- drop(cin_jac_1 %*% dx_unpack)

      # strict complementarity: a bound is only ACTIVE (binding) when the
      # solution sits at the bound AND the multiplier exerts force on the
      # gradient (see lav_con_cin_inactive_idx()); an interior solution
      # that merely grazes a bound keeps its ordinary standard error
      inactive_idx <- lav_con_cin_inactive_idx(
        con0 = con0, lambda = cin_lambda, jac = cin_jac_1,
        cin_flag = rep(TRUE, nrow(cin_jac_1))
      )
      cin_lambda[inactive_idx] <- 0

      # remove all inactive rows
      #if (length(inactive.idx) > 0L) {
      #  cin.JAC <- cin.JAC[-inactive.idx, , drop = FALSE]
      #  cin.lambda <- cin.lambda[-inactive.idx]
      #  inactive.idx <- integer(0L)
      #}
    } else {
      npar <- length(lavpartable$free[lavpartable$free > 0])
      cin_jac_1 <- matrix(0, nrow = 0L, ncol = npar)
      inactive_idx <- integer(0L)
      cin_lambda <- numeric(0L)
    }

    if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L) {
      ceq_jac_1 <- t(lav_mat_ortho_complement(lavmodel@ceq.simple.K))
      ceq_lambda <- numeric(nrow(ceq_jac_1))
    } else {
      npar <- length(lavpartable$free[lavpartable$free > 0])
      ceq_jac_1 <- matrix(0, nrow = 0L, ncol = npar)
      ceq_lambda <- numeric(0L)
    }

    # combine
    jac <- rbind(cin_jac_1, ceq_jac_1)
    attr(jac, "inactive.idx") <- inactive_idx
    attr(jac, "cin.idx") <- seq_len(nrow(cin_jac_1))
    attr(jac, "ceq.idx") <- nrow(cin_jac_1) + seq_len(nrow(ceq_jac_1))
    lambda <- c(cin_lambda, ceq_lambda)

    optim_out$con.jac <- jac
    optim_out$lambda <- lambda
  }


  fx <- objective_function(x) # to get "fx.group" attribute

  # check convergence
  warn_txt <- ""
  if (converged) {
    # check.gradient
    if (!is.null(gradient) &&
      optimizer %in% c("NLMINB", "BFGS", "L.BFGS.B")) {
      # gradient in the z (optimizer) metric: when parameter scaling is
      # active, this is the right metric for an absolute tolerance --
      # the z parameters live on a standardized-world scale, while the
      # RAW gradient components of (say) a variance of order 1e-6 are
      # amplified by the inverse scale and could never pass a fixed
      # cutoff, even at a machine-precision optimum
      dx <- gradient(x)

      if (converged && lavoptions$check.gradient &&
        any(abs(dx) > lavoptions$optim.dx.tol)) {
        # ok, identify the non-zero elements
        non_zero <- which(abs(dx) > lavoptions$optim.dx.tol)

        # which ones are 'boundary' points, defined by lower/upper?
        bound_idx <- integer(0L)
        if (!is.null(lavpartable$lower)) {
          bound_idx <- c(bound_idx, which(lower == x))
        }
        if (!is.null(lavpartable$upper)) {
          bound_idx <- c(bound_idx, which(upper == x))
        }
        # parameters fixed at their h1 values (saturated blocks) were
        # not part of the search
        if (length(h1_sat_idx) > 0L) {
          bound_idx <- c(bound_idx, h1_sat_idx)
        }
        if (length(bound_idx) > 0L) {
          non_zero <- non_zero[-which(non_zero %in% bound_idx)]
        }

        # this has many implications ... so should be careful to
        # avoid false alarm
        if (length(non_zero) > 0L) {
          converged <- FALSE
          warn_txt <- paste("the optimizer claimed the model converged,\n",
            "       but not all elements of the gradient are (near) zero;\n",
            "       the optimizer may not have found a local solution;\n",
            "       use check_gradient = FALSE to skip this check.",
            sep = ""
          )
        }
      }
    } else {
      dx <- numeric(0L)
    }
  } else {
    dx <- numeric(0L)
    warn_txt <- "the optimizer warns that a solution has NOT been found!"
  }

  # transform back
  # 3.
  # if(lavoptions$optim.var.transform == "sqrt" &&
  #       length(lavmodel@x.free.var.idx) > 0L) {
  #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
  #    x.var <- x[lavmodel@x.free.var.idx]
  #    x.var.sign <- sign(x.var)
  #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
  # }

  # 3. unscale
  x <- x / scale

  # 2. unpack
  if (lavmodel@eq.constraints) {
    x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
      lavmodel@eq.constraints.k0
  }

  # runaway solution? (residual variance far more negative than the
  # observed variance: a drift towards an infimum at infinity, not a
  # stationary point); the caller may decide to try again
  runaway <- NULL
  if (converged && lavoptions$estimator != "PML") {
    runaway <- lav_model_est_runaway(
      x = x, lavpartable = lavpartable, lavsamplestats = lavsamplestats,
      lavh1 = lavh1, lavdata = lavdata, lavoptions = lavoptions
    )
    if (!is.null(runaway)) {
      warn_txt <- paste0(
        gettext("the optimizer claimed the model converged,\n"),
        gettext("       but the solution seems to have run away:\n"),
        paste0("       ", gettextf(
          "the estimated residual variance of %1$s is %2$s, while its observed variance is only %3$s",
          runaway$name, formatC(runaway$est, digits = 4, format = "g"),
          formatC(runaway$obs, digits = 4, format = "g")
        ), collapse = ";\n"),
        gettext(";\n       the objective function may not have a minimum (an over-factored or otherwise unidentified model?); consider adding bounds (e.g., bounds = \"pos.var\").")
      )
    }
  }

  attr(x, "converged") <- converged
  attr(x, "runaway") <- !is.null(runaway)
  attr(x, "start") <- start_x
  attr(x, "warn.txt") <- warn_txt
  attr(x, "iterations") <- iterations
  attr(x, "control") <- control
  attr(x, "fx") <- fx
  attr(x, "dx") <- dx
  attr(x, "parscale") <- parscale
  attr(x, "parscale_packed") <- scale
  if (!is.null(optim_out$con.jac)) attr(x, "con.jac") <- optim_out$con.jac
  if (!is.null(optim_out$lambda)) attr(x, "con.lambda") <- optim_out$lambda
  if (lavoptions$optim.partrace) {
    attr(x, "partrace") <- penv$PARTRACE
  }

  x
}

# backwards compatibility
# estimateModel <- lav_model_estimate

# new in 0.7-1: identify free parameters that belong to a 'saturated'
# block of a multilevel model
#
# a block is considered to be saturated if it only contains (free,
# unconstrained) variances/covariances -- and means/intercepts -- of the
# observed variables of that block; in that case, the 'estimated' values of
# these parameters are already available in lavh1$implied, and they do not
# need to enter the optimization again as free parameters: they are
# (temporarily) fixed at their h1 values during optimization only; once
# model estimation is over, they are treated as free parameters again
# (eg to compute standard errors)
#
# returns a list with two elements:
#   - x.idx: positions of these parameters in the (packed) parameter vector
#   - value: the corresponding h1 estimates
lav_model_est_h1_saturated <- function(lavmodel = NULL,
                                       lavpartable = NULL,
                                       lavdata = NULL,
                                       lavh1 = NULL) {
  empty <- list(x.idx = integer(0L), value = numeric(0L))

  # multilevel only (for now), and we need the h1 estimates
  if (lavdata@nlevels == 1L || length(lavh1) == 0L ||
      is.null(lavh1$implied$cov)) {
    return(empty)
  }

  # general linear equality constraints: the optimizer operates in a
  # reduced (packed) parameter space; we cannot fix individual parameters
  if (lavmodel@eq.constraints) {
    return(empty)
  }

  nlevels <- lavdata@nlevels
  ngroups <- lavdata@ngroups

  # position of each free parameter in the (packed) parameter vector
  free_id <- lavpartable$free[lavpartable$free > 0L]
  uid <- unique(free_id) # ceq.simple: duplicated ids are packed out
  dup_id <- unique(free_id[duplicated(free_id)])

  # labels that show up in (in)equality constraints
  con_idx <- which(lavpartable$op %in% c("==", "<", ">"))
  con_ref <- unique(c(lavpartable$lhs[con_idx], lavpartable$rhs[con_idx]))

  x_idx <- integer(0L)
  value <- numeric(0L)

  for (g in seq_len(ngroups)) {
    for (l in seq_len(nlevels)) {
      b <- (g - 1L) * nlevels + l
      ovn <- lavdata@ov.names.l[[g]][[l]]
      p <- length(ovn)
      cov_h1 <- lavh1$implied$cov[[b]]
      mean_h1 <- lavh1$implied$mean[[b]]
      if (p == 0L || is.null(cov_h1) || nrow(cov_h1) != p ||
          length(mean_h1) != p) {
        next
      }

      row_idx <- which(lavpartable$block == b)
      if (length(row_idx) == 0L) {
        next
      }
      op <- lavpartable$op[row_idx]
      lhs <- lavpartable$lhs[row_idx]
      rhs <- lavpartable$rhs[row_idx]
      free <- lavpartable$free[row_idx]

      # only variances/covariances and intercepts of observed variables
      # (no latent variables, no regressions)
      if (any(!op %in% c("~~", "~1")) ||
          any(!lhs %in% ovn) ||
          any(op == "~~" & !rhs %in% ovn)) {
        next
      }

      # all p*(p+1)/2 variances/covariances must be present and free
      cov_flag <- op == "~~"
      if (any(free[cov_flag] == 0L)) {
        next
      }
      i1 <- match(lhs[cov_flag], ovn)
      i2 <- match(rhs[cov_flag], ovn)
      pair_id <- pmin(i1, i2) + p * pmax(i1, i2)
      if (length(unique(pair_id)) != p * (p + 1L) / 2L) {
        next
      }

      # intercepts/means: all present; either free, or fixed at their h1
      # value (eg the zero within-level means of shared variables)
      int_flag <- op == "~1"
      if (lavmodel@meanstructure) {
        if (sum(int_flag) != p || any(duplicated(lhs[int_flag]))) {
          next
        }
        int_fixed <- which(int_flag & free == 0L)
        if (length(int_fixed) > 0L) {
          fixed_val <- lavpartable$start[row_idx[int_fixed]]
          h1_val <- mean_h1[match(lhs[int_fixed], ovn)]
          if (any(abs(fixed_val - h1_val) > 1e-10)) {
            next
          }
        }
      } else if (any(int_flag)) {
        next
      }

      # none of the free parameters may be involved in equality or
      # inequality constraints
      block_id <- free[free > 0L]
      if (any(block_id %in% dup_id)) {
        next
      }
      lab <- c(lavpartable$label[row_idx], lavpartable$plabel[row_idx])
      lab <- lab[nzchar(lab)]
      if (length(con_ref) > 0L && any(lab %in% con_ref)) {
        next
      }

      # this block is saturated: collect positions + h1 values
      b_free <- which(free > 0L)
      b_value <- numeric(length(b_free))
      for (i in seq_along(b_free)) {
        r <- b_free[i]
        if (op[r] == "~~") {
          b_value[i] <- cov_h1[match(lhs[r], ovn), match(rhs[r], ovn)]
        } else {
          b_value[i] <- mean_h1[match(lhs[r], ovn)]
        }
      }
      x_idx <- c(x_idx, match(free[b_free], uid))
      value <- c(value, b_value)
    } # levels
  } # groups

  list(x.idx = x_idx, value = value)
}


# detect a 'runaway' solution
#
# the optimizer may report convergence at a point that is not a stationary
# point at all, but lies on a path towards an infimum at infinity; the
# classic case is an over-factored (efa) block where a factor collapses
# onto a single indicator: its loading grows without bound, the residual
# variance of that indicator goes to minus infinity (so that the implied
# variance stays put), and the objective creeps towards a limit it never
# reaches; because the drift is so slow, nlminb eventually satisfies its
# relative tolerance, and at a rather arbitrary point
#
# we flag a solution as 'runaway' when a free residual variance of a
# (continuous) observed variable is more negative than the observed
# variance of that variable itself (ratio < -tol); a proper local minimum
# with a mild Heywood case (a slightly negative residual variance) is not
# affected
#
# returns NULL if nothing is wrong, otherwise a data.frame with the
# offending variables, their estimated residual variances, and their
# observed variances
lav_model_est_runaway <- function(x = NULL, lavpartable = NULL,
                                  lavsamplestats = NULL, lavh1 = NULL,
                                  lavdata = NULL, lavoptions = NULL,
                                  tol = 1) {
  # only single-level (multilevel: TODO)
  if (is.null(lavdata) || lavdata@nlevels > 1L || is.null(lavsamplestats) ||
      is.null(lavpartable$free)) {
    return(NULL)
  }

  if (is.null(lavpartable$group)) {
    lavpartable$group <- rep(1L, length(lavpartable$lhs))
  }
  group_values <- lav_pt_group_values(lavpartable)
  ngroups <- length(group_values)
  if (ngroups != lavsamplestats@ngroups) {
    return(NULL)
  }

  # continuous observed variables only
  ov_cont <- lavdata@ov$name[lavdata@ov$type == "numeric"]

  out <- NULL
  for (g in seq_len(ngroups)) {
    ov_names <- lavdata@ov.names[[g]]

    # observed variances for this group (same ordering as ov.names)
    ov_var <- NULL
    if (lavsamplestats@missing.flag) {
      if (!is.null(lavh1$implied$cov[[g]])) {
        ov_var <- diag(lavh1$implied$cov[[g]])
      } else if (!is.null(lavsamplestats@missing.h1[[g]]$sigma)) {
        ov_var <- diag(lavsamplestats@missing.h1[[g]]$sigma)
      }
    } else if (isTRUE(lavoptions$conditional.x)) {
      ov_var <- diag(lavsamplestats@res.cov[[g]])
    } else {
      ov_var <- diag(lavsamplestats@cov[[g]])
    }
    if (is.null(ov_var) || length(ov_var) != length(ov_names)) {
      next
    }

    # free residual variances of continuous observed variables
    row_idx <- which(lavpartable$op == "~~" &
                     lavpartable$lhs == lavpartable$rhs &
                     lavpartable$free > 0L &
                     lavpartable$group == group_values[g] &
                     lavpartable$lhs %in% ov_names &
                     lavpartable$lhs %in% ov_cont)
    if (length(row_idx) == 0L) {
      next
    }
    est <- x[lavpartable$free[row_idx]]
    obs <- ov_var[match(lavpartable$lhs[row_idx], ov_names)]
    bad <- which(is.finite(est) & is.finite(obs) & obs > 0 & est < -tol * obs)
    if (length(bad) > 0L) {
      out <- rbind(out, data.frame(
        name = lavpartable$lhs[row_idx[bad]],
        group = rep(g, length(bad)),
        est = est[bad],
        obs = obs[bad],
        stringsAsFactors = FALSE
      ))
    }
  }

  out
}
