# Gauss-Newton optimization: Fisher scoring with Levenberg-Marquardt
# globalization
#
# The search direction solves the damped normal equations
#
#     (H + mu * D) d = -g
#
# where g is the (scoring) gradient of the discrepancy function,
# H = t(Delta) %*% A1 %*% Delta is the expected (Gauss-Newton)
# information, and D = diag(H) is the Marquardt scaling. The damping
# parameter mu is adapted from the gain ratio (actual versus predicted
# reduction): mu -> 0 recovers the undamped Gauss-Newton step (fast
# local convergence), while a large mu turns the step into a short
# (scaled) steepest-descent step (globally safe). A step is only taken
# when it decreases the objective (or, for estimators whose scoring
# gradient is not the exact gradient of the monitored objective, when it
# decreases the scoring residual; see below).
#
# constraints:
# - linear equality constraints are eliminated exactly by optimizing in
#   the packed space x = K %*% z + k0 (the same route as nlminb); as in
#   lav_model_est, bounds are ignored in this case
# - simple equality constraints (ceq.simple) are handled by reducing the
#   gradient/information via ceq.simple.K
# - simple bounds (the lower/upper columns of the parameter table) are
#   honored with projected steps plus an active set: parameters pinned at
#   a bound with an outward-pointing gradient are frozen, the damped
#   system is solved for the remaining parameters, and the trial point is
#   projected onto the box
#
# estimator = "DLS" with dls.GammaNT = "model": the weight matrix is a
# function of the model-implied moments, and the (IRLS) scoring gradient
# omits the dV/dtheta terms of the monitored objective. The estimator is
# defined by the scoring fixed point g(theta) = 0, so near convergence
# the objective value may (slightly) increase along the scoring
# direction; such steps are accepted whenever they shrink the scoring
# residual ||g||.

# YR - 19 Jan 2021: initial version, needed for DLS
# YR - 10 Jul 2026: rewritten: Levenberg-Marquardt damping, packed
#      equality constraints, active-set bounds, Delta/A1 computed once
#      per iteration (shared by gradient and information), streamed
#      ML/GLS kernels via lav_model_grad/lav_model_info_expected,
#      gradient-based convergence test

# evaluate the objective at x; when objective_only = FALSE, also return
# the scoring gradient g and the expected (Gauss-Newton) information H,
# both in the free-parameter space (nx.free)
lav_optim_gn_eval <- function(x, lavmodel = NULL, lavsamplestats = NULL,
                              lavdata = NULL, lavoptions = NULL,
                              objective_only = TRUE,
                              fast_ls = FALSE, group_w = NULL,
                              group_weight = TRUE) {
  lavmodel2 <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
  fx <- lav_model_objective(
    lavmodel = lavmodel2, lavdata = lavdata,
    lavsamplestats = lavsamplestats
  )
  fx <- as.numeric(fx)[1]

  if (objective_only) {
    return(list(fx = fx, g = NULL, H = NULL))
  }

  if (fast_ls) {
    # single-level moment-based estimators (WLS/DWLS/ULS/DLS): compute
    # Delta and the weight matrix A1 once, and share them between the
    # gradient and the information; for DLS with dls.GammaNT = "model"
    # this builds the (expensive) theta-dependent weight matrix a single
    # time per iteration
    lavimplied <- lav_model_implied(lavmodel = lavmodel2)
    delta <- lav_model_delta(lavmodel = lavmodel2)
    a1 <- lav_model_h1_info_expected(
      lavobject = NULL,
      lavmodel = lavmodel2,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavimplied = lavimplied,
      lavh1 = NULL,
      lavcache = NULL
    )
    wls_est <- lav_model_wls_est(
      lavmodel = lavmodel2,
      lavimplied = lavimplied
    )
    diagonal_v <- lavmodel@estimator %in% c("DWLS", "ULS")

    g <- H <- NULL
    for (gr in seq_len(lavsamplestats@ngroups)) {
      diff_g <- lavsamplestats@WLS.obs[[gr]] - wls_est[[gr]]
      if (diagonal_v) {
        # t(Delta) %*% diag(a1)
        pre_g <- t(a1[[gr]] * delta[[gr]])
      } else {
        pre_g <- crossprod(delta[[gr]], a1[[gr]])
      }
      g_gr <- group_w[gr] * drop(pre_g %*% diff_g)
      h_gr <- group_w[gr] * (pre_g %*% delta[[gr]])
      if (gr == 1L) {
        g <- -1 * g_gr
        H <- h_gr
      } else {
        g <- g - g_gr
        H <- H + h_gr
      }
    }

    # simple equality constraints: reduce from nx.unco to nx.free
    if (lavmodel@ceq.simple.only &&
      length(g) == lavmodel@nx.unco &&
      lavmodel@nx.unco > lavmodel@nx.free) {
      k_s <- lavmodel@ceq.simple.K
      g <- drop(crossprod(k_s, g))
      H <- crossprod(k_s, H %*% k_s)
    }
  } else {
    # generic route: the same analytic gradient as the nlminb path, and
    # the expected information; both use the streamed ML/GLS kernels
    # when available (the pstar x pstar weight matrix is never formed)
    g <- lav_model_grad(
      lavmodel = lavmodel2,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      type = "free",
      group_weight = group_weight,
      ceq_simple = TRUE
    )
    g <- as.numeric(g)
    H <- lav_model_info_expected(
      lavmodel = lavmodel2,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    # simple equality constraints: reduce from nx.unco to nx.free
    # (the gradient is already reduced via ceq_simple = TRUE)
    if (lavmodel@ceq.simple.only &&
      NCOL(H) == lavmodel@nx.unco &&
      lavmodel@nx.unco > lavmodel@nx.free) {
      k_s <- lavmodel@ceq.simple.K
      H <- crossprod(k_s, H %*% k_s)
    }
  }

  # enforce symmetry (guards against round-off in the group sums)
  H <- (H + t(H)) / 2

  list(fx = fx, g = g, H = H)
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

  # no support (yet) for general inequality constraints. Simple ones are
  # fine: when cin.simple.only is TRUE, every inequality constraint (and
  # every bounds= bound, which is encoded as a cin.function too) has been
  # lowered into the lower/upper columns of the parameter table, which the
  # projected Gauss-Newton iterations below honor natively.
  if (!is.null(body(lavmodel@cin.function)) && !lavmodel@cin.simple.only) {
    lav_msg_stop(gettext(
      "inequality constraints not supported (yet) with optim.method = \"GN\"."))
  }

  # only for estimators with an analytic (scoring) gradient and an
  # expected information in moment space
  if (!lavmodel@estimator %in% c(
    "ML", "GLS", "WLS", "DWLS", "ULS", "DLS", "NTRLS", "catML"
  )) {
    lav_msg_stop(gettextf(
      "optim.method = \"GN\" is not available for estimator %s.",
      dQuote(lavmodel@estimator)))
  }

  # options (gn.args sublist)
  max_iter <- as.integer(lavoptions$gn.args$max_iter)
  tol_x <- lavoptions$gn.args$tol_x
  tol_g <- lavoptions$gn.args$tol_g

  # current set of free parameters
  x <- lav_model_get_parameters(lavmodel)
  npar <- length(x)

  # nothing to do?
  if (npar == 0L) {
    fx <- lav_model_objective(
      lavmodel = lavmodel, lavdata = lavdata,
      lavsamplestats = lavsamplestats
    )
    attr(x, "converged") <- TRUE
    attr(x, "warn.txt") <- ""
    attr(x, "iterations") <- 0L
    attr(x, "control") <- list(iter.max = max_iter, tol.x = tol_x,
                               tol.g = tol_g)
    attr(x, "fx") <- fx
    return(x)
  }

  # linear (non-simple) equality constraints: optimize in the packed
  # space x = K %*% z + k0 (the constraints then hold exactly along the
  # whole path); this mirrors the nlminb route in lav_model_est
  eq_flag <- (lavmodel@eq.constraints &&
    ncol(lavmodel@eq.constraints.K) > 0L)
  if (eq_flag) {
    eq_k <- lavmodel@eq.constraints.K
    eq_k0 <- lavmodel@eq.constraints.k0
  }
  unpack <- function(z) {
    if (eq_flag) {
      as.numeric(eq_k %*% z) + eq_k0
    } else {
      z
    }
  }

  # simple bounds; as in lav_model_est, bounds are dropped in the
  # presence of (non-simple) linear equality constraints (the packed
  # coordinates are linear combinations of the original parameters)
  lb <- rep(-Inf, npar)
  ub <- rep(+Inf, npar)
  if (!eq_flag && !is.null(lavpartable) &&
    (!is.null(lavpartable$lower) || !is.null(lavpartable$upper))) {
    if (lavmodel@ceq.simple.only) {
      free_idx <- which(lavpartable$free > 0L &
        !duplicated(lavpartable$free))
    } else {
      free_idx <- which(lavpartable$free > 0L)
    }
    if (!is.null(lavpartable$lower)) {
      lb <- lavpartable$lower[free_idx]
    }
    if (!is.null(lavpartable$upper)) {
      ub <- lavpartable$upper[free_idx]
    }
    stopifnot(length(lb) == npar, length(ub) == npar)
    # inconsistent bounds may occur with equality constraints (qr() may
    # switch the sign); disable them
    bad_idx <- which(lb > ub)
    if (length(bad_idx) > 0L) {
      lb[bad_idx] <- -Inf
      ub[bad_idx] <- +Inf
    }
    # make the starting values feasible
    x <- pmin(pmax(x, lb), ub)
  }
  bounds_flag <- any(is.finite(lb)) || any(is.finite(ub))

  # group weights for the fast *LS route, matching the analytic gradient
  # (lav_model_grad) so that g is the exact gradient of the objective
  ngroups <- lavsamplestats@ngroups
  group_w <- numeric(ngroups)
  for (g in seq_len(ngroups)) {
    if (lavmodel@estimator == "DLS" &&
      !isTRUE(lavmodel@estimator.args$dls.FtimesNminus1)) {
      group_w[g] <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    } else {
      group_w[g] <- (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@ntotal
    }
  }

  # single-level moment-based estimators: share Delta/A1 between the
  # gradient and the information (see lav_optim_gn_eval)
  fast_ls <- (lavdata@nlevels == 1L &&
    lavmodel@estimator %in% c("WLS", "DWLS", "ULS", "DLS") &&
    !lavmodel@group.w.free &&
    !lavsamplestats@missing.flag)

  # generic route: same group-weight switch as lav_model_est
  group_weight <- !(lavsamplestats@missing.flag || lavdata@nlevels > 1L)

  # is the scoring gradient the exact gradient of the objective? not for
  # DLS with a model-based weight matrix (see the top of this file)
  scoring_flag <- (lavmodel@estimator == "DLS" &&
    lavmodel@estimator.args$dls.GammaNT == "model")

  eval_point <- function(x_1, objective_only = TRUE) {
    lav_optim_gn_eval(
      x = x_1, lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavoptions = lavoptions, objective_only = objective_only,
      fast_ls = fast_ls, group_w = group_w, group_weight = group_weight
    )
  }

  # move the gradient/information to the packed space (if needed)
  pack_point <- function(fit) {
    if (eq_flag && !is.null(fit$g)) {
      fit$g <- as.numeric(fit$g %*% eq_k)
      fit$H <- crossprod(eq_k, fit$H %*% eq_k)
    }
    fit
  }

  # working coordinates
  if (eq_flag) {
    z <- as.numeric((x - eq_k0) %*% eq_k)
    lb_z <- rep(-Inf, length(z))
    ub_z <- rep(+Inf, length(z))
  } else {
    z <- x
    lb_z <- lb
    ub_z <- ub
  }
  nz <- length(z)

  # initial evaluation
  cur <- pack_point(eval_point(unpack(z), objective_only = FALSE))
  warn_txt <- ""
  converged <- FALSE
  iter <- 0L

  if (!is.finite(cur$fx) || anyNA(cur$g) || anyNA(cur$H)) {
    x <- unpack(z)
    fx <- lav_model_objective(
      lavmodel = lav_model_set_parameters(lavmodel = lavmodel, x = x),
      lavdata = lavdata, lavsamplestats = lavsamplestats
    )
    attr(x, "converged") <- FALSE
    attr(x, "warn.txt") <- gettext(
      "objective or gradient could not be evaluated at the starting values.")
    attr(x, "iterations") <- 0L
    attr(x, "control") <- list(iter.max = max_iter, tol.x = tol_x,
                               tol.g = tol_g)
    attr(x, "fx") <- fx
    return(x)
  }

  if (lav_verbose()) {
    cat("iter = ", sprintf("%3d", 0L),
      ": obj = ", sprintf("%13.10f", cur$fx), "\n",
      sep = ""
    )
  }

  # scaled norm of the projected gradient (Dennis & Schnabel style
  # relative gradient); the projection removes the components that point
  # outward at an active bound
  relgrad <- function(z_1, fit) {
    pg <- fit$g
    if (bounds_flag) {
      pg[z_1 <= lb_z & pg > 0] <- 0
      pg[z_1 >= ub_z & pg < 0] <- 0
    }
    max(abs(pg) * pmax(abs(z_1), 1)) / max(abs(fit$fx), 1)
  }

  # Levenberg-Marquardt damping state
  mu <- 1e-04
  nu <- 2
  mu_max <- 1e+12
  reject_max <- 25L
  bound_eps <- 1e-12

  while (iter < max_iter) {
    # gradient-based convergence test
    if (relgrad(z, cur) < tol_g) {
      converged <- TRUE
      if (lav_verbose()) {
        cat("Gauss-Newton: relative gradient is below tolerance (",
          sprintf("%g", tol_g), ")\n",
          sep = ""
        )
      }
      break
    }

    iter <- iter + 1L

    # active set: parameters at a bound with an outward-pointing gradient
    # stay put this iteration
    if (bounds_flag) {
      active <- (z <= lb_z + bound_eps & cur$g > 0) |
        (z >= ub_z - bound_eps & cur$g < 0)
    } else {
      active <- logical(nz)
    }
    free_set <- which(!active)
    if (length(free_set) == 0L) {
      # all parameters pinned at a bound; nothing can improve
      converged <- TRUE
      break
    }

    h_ff <- cur$H[free_set, free_set, drop = FALSE]
    g_f <- cur$g[free_set]
    # Marquardt scaling, with a floor for (near-)zero curvature
    d_diag <- diag(h_ff)
    d_floor <- max(abs(d_diag), 1) * 1e-08
    d_diag <- pmax(abs(d_diag), d_floor)

    # inner loop: increase the damping until an acceptable step is found
    accepted <- FALSE
    new_full <- NULL
    n_reject <- 0L
    rho <- as.numeric(NA)
    while (!accepted && n_reject <= reject_max && mu <= mu_max) {
      m_damp <- h_ff
      diag(m_damp) <- diag(m_damp) + mu * d_diag
      r_chol <- try(chol(m_damp), silent = TRUE)
      if (inherits(r_chol, "try-error")) {
        # not positive definite at this damping level: increase and retry
        mu <- mu * nu
        nu <- 2 * nu
        n_reject <- n_reject + 1L
        next
      }
      d_step <- numeric(nz)
      d_step[free_set] <- -1 * backsolve(
        r_chol,
        backsolve(r_chol, g_f, transpose = TRUE)
      )

      # trial point (projected onto the box)
      z_new <- pmin(pmax(z + d_step, lb_z), ub_z)
      d_eff <- z_new - z
      # predicted reduction of the quadratic model along the actual
      # (projected) step
      pred <- -1 * (sum(cur$g * d_eff) +
        0.5 * sum(d_eff * drop(cur$H %*% d_eff)))

      trial <- eval_point(unpack(z_new), objective_only = TRUE)
      ared <- cur$fx - trial$fx

      if (is.finite(trial$fx) && ared > 0) {
        # objective decreased: accept, and relax the damping according
        # to the gain ratio (Nielsen)
        accepted <- TRUE
        if (is.finite(pred) && pred > 0) {
          rho <- ared / pred
        } else {
          rho <- 1
        }
        mu <- max(mu * max(1 / 3, 1 - (2 * rho - 1)^3), 1e-12)
        nu <- 2
      } else if (scoring_flag && is.finite(trial$fx)) {
        # scoring fallback: the (IRLS) fixed point is not the exact
        # minimizer of the monitored objective; accept the step if it
        # shrinks the scoring residual
        new_full <- pack_point(eval_point(unpack(z_new),
          objective_only = FALSE
        ))
        if (!anyNA(new_full$g) && !anyNA(new_full$H) &&
          sqrt(sum(new_full$g^2)) <
            (1 - 1e-04) * sqrt(sum(cur$g^2))) {
          accepted <- TRUE
          new_full$fx <- trial$fx
          rho <- as.numeric(NA)
          nu <- 2
        } else {
          new_full <- NULL
        }
      }

      if (!accepted) {
        mu <- mu * nu
        nu <- 2 * nu
        n_reject <- n_reject + 1L
      }
    } # inner

    if (!accepted) {
      # no acceptable step, even with (very) heavy damping: we are either
      # at a (numerically) stationary point, or stuck
      if (relgrad(z, cur) < sqrt(tol_g)) {
        converged <- TRUE
      } else {
        warn_txt <- gettext(
          "the optimizer could no longer improve the objective
          (Gauss-Newton step rejected at maximum damping).")
      }
      break
    }

    # accepted step: update the state
    x_old <- unpack(z)
    z <- z_new
    x_new <- unpack(z)
    rms_x <- sqrt(mean((x_new - x_old)^2))

    if (is.null(new_full)) {
      new_full <- pack_point(eval_point(x_new, objective_only = FALSE))
      new_full$fx <- trial$fx
    }
    cur <- new_full

    if (lav_verbose()) {
      # (kept within 80 columns)
      cat("iter = ", sprintf("%3d", iter),
        ": obj = ", sprintf("%13.10f", cur$fx),
        " mu = ", sprintf("%7.1e", mu),
        if (is.finite(rho)) {
          paste0(" rho = ", sprintf("%6.3f", rho))
        } else {
          " (scoring)"
        },
        " rms.x = ", sprintf("%12.10f", rms_x), "\n",
        sep = ""
      )
    }

    if (anyNA(cur$g) || anyNA(cur$H)) {
      warn_txt <- gettext(
        "the gradient could not be evaluated at the updated parameter
        values.")
      break
    }

    # step-based convergence test (backward compatible with < 0.7)
    if (rms_x < tol_x) {
      converged <- TRUE
      if (lav_verbose()) {
        cat("Gauss-Newton algorithm converged: rms.x = ",
          sprintf("%12.12f", rms_x), " < ",
          sprintf("%12.12f", tol_x), "\n",
          sep = ""
        )
      }
      break
    }
  } # outer

  x <- unpack(z)

  # one last evaluation, to get the fx.group attribute
  lavmodel <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
  fx <- lav_model_objective(
    lavmodel = lavmodel, lavdata = lavdata,
    lavsamplestats = lavsamplestats
  )

  # add attributes
  if (!converged && iter >= max_iter && !nzchar(warn_txt)) {
    warn_txt <- paste("maximum number of iterations (",
      max_iter, ") ",
      "was reached without convergence.\n",
      sep = ""
    )
  }
  attr(x, "converged") <- converged
  attr(x, "warn.txt") <- warn_txt
  attr(x, "iterations") <- iter
  attr(x, "control") <- list(
    iter.max = max_iter,
    tol.x = tol_x,
    tol.g = tol_g
  )
  attr(x, "fx") <- fx

  x
}
