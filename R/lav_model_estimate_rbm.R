# Reduced-bias M-estimation (RBM) for continuous-outcome SEM.
#
# RBM is a penalized maximum likelihood method (Firth/Kosmidis-style). The
# discrepancy function, implied moments, loglik, test statistic and information
# matrices are all the ML ones, so the estimator is kept as "ML" internally
# (see lav_options_est_rbm) and only the point estimate differs:
#
#   penalty  P(theta) = -1/2 * tr( J(theta)^-1 E(theta) )
#     J = observed information, E = first-order (outer-product) information.
#
#   implicit (iRBM): maximize  loglik(theta) + P(theta)
#   explicit (eRBM): maximize  loglik(theta), then one-step bias-correct
#                    theta <- theta - bias(theta), bias = -J^-1 dP/dtheta
#
# Standard errors use the sandwich estimator (se = "robust.huber.white"),
# computed by lav_model_vcov() in the usual way.
#
# Reference: 
# - Jamil, H., Rosseel, Y., Kemp, O., & Kosmidis, I. (2026). Bias-Reduced 
#   Estimation of Structural Equation Models. Structural Equation Modeling: A 
#   Multidisciplinary Journal, 33(3), 376–392. 
#   https://doi.org/10.1080/10705511.2025.2610462
# - ported from the brlavaan package (github.com/haziqj/brlavaan)

# check whether a matrix is non-positive-definite or has NAs
lav_model_rbm_check_mat <- function(mat) {
  if (anyNA(mat)) {
    return(TRUE)
  }
  eig <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  any(eig < -1e-06 * eig[1])
}

# total (N-scaled) information matrix in the (free) parameter space
# kind: "observed" (J) or "first.order" (E)
lav_model_rbm_information <- function(x, lavmodel, lavsamplestats, lavdata,
                                      lavoptions, kind = "observed") {
  lavmodel <- lav_model_set_parameters(lavmodel, x)
  if (kind == "observed") {
    info <- lav_model_info_observed(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats,
      lavdata = lavdata, lavoptions = lavoptions
    )
  } else {
    info <- lav_model_info_firstorder(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats,
      lavdata = lavdata, lavoptions = lavoptions
    )
  }
  # the info functions return unit (per-observation) information
  info <- lavsamplestats@ntotal * info

  # reduce to the optimizer space for simple equality constraints
  if (lavmodel@ceq.simple.only) {
    k <- lavmodel@ceq.simple.K
    info <- t(k) %*% info %*% k
  }
  info
}

# RBM penalty P(theta) = -1/2 * tr( J^-1 E ), evaluated on the loglik scale
lav_model_rbm_penalty <- function(x, lavmodel, lavsamplestats, lavdata,
                                  lavoptions, fallback = -1e8) {
  e <- lav_model_rbm_information(x, lavmodel, lavsamplestats, lavdata,
    lavoptions,
    kind = "first.order"
  )
  j <- lav_model_rbm_information(x, lavmodel, lavsamplestats, lavdata,
    lavoptions,
    kind = "observed"
  )
  if (lav_model_rbm_check_mat(j)) {
    return(fallback)
  }
  jinv <- try(solve(j), silent = TRUE)
  if (inherits(jinv, "try-error")) {
    fallback
  } else {
    # sum(jinv * e) == tr(jinv %*% e) for symmetric jinv
    -sum(jinv * e) / 2
  }
}

# one-step bias for the explicit method: bias = -J^-1 dP/dtheta
lav_model_rbm_bias <- function(x, lavmodel, lavsamplestats, lavdata,
                               lavoptions) {
  j <- lav_model_rbm_information(x, lavmodel, lavsamplestats, lavdata,
    lavoptions,
    kind = "observed"
  )
  if (lav_model_rbm_check_mat(j)) {
    return(rep(as.numeric(NA), length(x)))
  }
  jinv <- try(solve(j), silent = TRUE)
  if (inherits(jinv, "try-error")) {
    return(rep(as.numeric(NA), length(x)))
  }
  a <- numDeriv::grad(
    func = lav_model_rbm_penalty, x = x,
    lavmodel = lavmodel, lavsamplestats = lavsamplestats,
    lavdata = lavdata, lavoptions = lavoptions
  )
  -drop(jinv %*% a)
}

# RBM estimation; returns the free parameter vector 'x' with the same
# attributes as lav_model_est() so that lav_step11_estoptim() can use it.
lav_model_est_rbm <- function(lavmodel = NULL,
                              lavpartable = NULL,
                              lavsamplestats = NULL,
                              lavdata = NULL,
                              lavoptions = NULL,
                              lavcache = list()) {
  rbm_method <- lavoptions$estimator.args$rbm.method
  if (is.null(rbm_method)) {
    rbm_method <- "implicit"
  }

  # general (non-simple) equality constraints are not supported yet
  if (lavmodel@eq.constraints) {
    lav_msg_stop(gettext(
      "estimator \"rbm\" does not support general equality constraints yet;
       use simple (label-based) equality constraints instead."))
  }

  ntotal <- lavsamplestats@ntotal

  # ----- explicit / none: reuse the standard ML optimizer -----
  if (rbm_method %in% c("explicit", "none")) {
    # lav_model_est() dispatches on optim.method; restore the ML optimizer
    lavoptions_ml <- lavoptions
    lavoptions_ml$optim.method <- "nlminb"
    x_ml <- lav_model_est(
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavoptions = lavoptions_ml, lavcache = lavcache
    )
    if (rbm_method == "none") {
      return(x_ml)
    }

    # explicit: one-step bias correction
    b <- lav_model_rbm_bias(
      as.numeric(x_ml), lavmodel, lavsamplestats, lavdata, lavoptions
    )
    x <- as.numeric(x_ml) - b

    lavmodel_x <- lav_model_set_parameters(lavmodel, x)
    fx <- lav_model_objective(
      lavmodel = lavmodel_x, lavsamplestats = lavsamplestats,
      lavdata = lavdata, lavcache = lavcache
    )

    attr(x, "converged") <- attr(x_ml, "converged")
    attr(x, "iterations") <- attr(x_ml, "iterations")
    attr(x, "control") <- attr(x_ml, "control")
    attr(x, "warn.txt") <- attr(x_ml, "warn.txt")
    attr(x, "fx") <- fx
    attr(x, "dx") <- numeric(0L)
    attr(x, "parscale") <- attr(x_ml, "parscale")
    return(x)
  }

  # ----- implicit: penalized-ML optimization -----

  # starting values (in the optimizer/free space)
  start_x <- lav_model_get_parameters(lavmodel)

  # bounds (mirror lav_model_est; eq.constraints excluded above)
  if (is.null(lavpartable$lower)) {
    lower <- -Inf
  } else if (lavmodel@ceq.simple.only) {
    free_idx <- which(lavpartable$free > 0L & !duplicated(lavpartable$free))
    lower <- lavpartable$lower[free_idx]
  } else {
    lower <- lavpartable$lower[lavpartable$free > 0L]
  }
  if (is.null(lavpartable$upper)) {
    upper <- +Inf
  } else if (lavmodel@ceq.simple.only) {
    free_idx <- which(lavpartable$free > 0L & !duplicated(lavpartable$free))
    upper <- lavpartable$upper[free_idx]
  } else {
    upper <- lavpartable$upper[lavpartable$free > 0L]
  }
  bad_idx <- which(lower > upper)
  if (length(bad_idx) > 0L) {
    lower[bad_idx] <- -Inf
    upper[bad_idx] <- +Inf
  }

  # theta-independent part of the Gaussian loglik (loglik = const - N * fx_ML).
  # Optimizing the penalized negative loglik on its natural scale (rather than
  # the small N * fx_ML - P) is what lets nlminb reach 'relative convergence'
  # (matching brlavaan); the constant does not affect the estimate.
  log2pi <- log(2 * pi)
  loglik_const <- 0
  for (g in seq_len(lavsamplestats@ngroups)) {
    ng <- lavsamplestats@nobs[[g]]
    pg <- NROW(lavsamplestats@cov[[g]])
    loglik_const <- loglik_const -
      0.5 * ng * (pg * log2pi + lavsamplestats@cov.log.det[[g]] + pg)
  }
  if (!is.finite(loglik_const)) {
    loglik_const <- 0
  }

  # penalized objective: minimize  -(loglik + P) = N * fx_ML - P - const
  objective_function <- function(x) {
    lavmodel_x <- lav_model_set_parameters(lavmodel, x)
    fx <- lav_model_objective(
      lavmodel = lavmodel_x, lavsamplestats = lavsamplestats,
      lavdata = lavdata, lavcache = lavcache
    )
    if (!is.finite(fx)) {
      return(1e20)
    }
    pen <- lav_model_rbm_penalty(
      x, lavmodel, lavsamplestats, lavdata, lavoptions
    )
    val <- ntotal * as.numeric(fx) - pen - loglik_const
    if (!is.finite(val)) 1e20 else val
  }

  # nlminb control (mirror lav_model_est)
  control_nlminb <- list(
    eval.max = 20000L, iter.max = 10000L, trace = 0L,
    abs.tol = (.Machine$double.eps * 10), rel.tol = 1e-10,
    step.min = 1.0, step.max = 1.0, x.tol = 1.5e-8, xf.tol = 2.2e-14
  )
  if (!is.null(lavoptions$control)) {
    control_nlminb <- modifyList(control_nlminb, lavoptions$control)
  }
  control <- control_nlminb[c(
    "eval.max", "iter.max", "trace", "step.min", "step.max",
    "abs.tol", "rel.tol", "x.tol", "xf.tol"
  )]

  optim_out <- nlminb(
    start = start_x, objective = objective_function, gradient = NULL,
    lower = lower, upper = upper, control = control
  )

  x <- optim_out$par

  # ML objective (with fx.group attribute) at the solution, for steps 12-14
  lavmodel_x <- lav_model_set_parameters(lavmodel, x)
  fx <- lav_model_objective(
    lavmodel = lavmodel_x, lavsamplestats = lavsamplestats,
    lavdata = lavdata, lavcache = lavcache
  )

  converged <- optim_out$convergence == 0L
  attr(x, "converged") <- converged
  attr(x, "iterations") <- optim_out$iterations
  attr(x, "control") <- control
  attr(x, "warn.txt") <- if (converged) {
    ""
  } else {
    "the optimizer warns that a solution has NOT been found!"
  }
  attr(x, "fx") <- fx
  attr(x, "dx") <- numeric(0L)
  attr(x, "parscale") <- rep(1.0, length(x))
  x
}
