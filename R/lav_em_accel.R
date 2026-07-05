# acceleration of EM-type fixed-point iterations
#
# YR + Claude, July 2026
#
# SQUAREM-1 (Varadhan & Roland, 2008, Scand. J. Statist. 35, 335-353):
# each cycle takes two EM steps, extrapolates along the 'squared' step,
# and then takes one stabilizing EM step from the extrapolated point;
# a loglikelihood safeguard falls back to the plain (double) EM step
# whenever the extrapolation fails or decreases the loglikelihood too much

# theta     = packed parameter vector
# step_fn   = the EM map: theta -> theta' (one E-step + M-step); the map
#             must return a valid parameter vector for any theta it accepts
# logl_fn   = loglikelihood at theta (to be maximized); used to safeguard
#             the extrapolation and to monitor convergence
# fx0       = loglikelihood at the initial theta (recomputed if NULL)
# conv      = convergence criterion (per cycle): "logl" = |change in logl|
#             < tol (as in the plain multilevel EM iterations); "param" =
#             largest absolute change in the parameter values < tol (as in
#             the plain single-level EM iterations)
# conv_fn   = optional user-supplied convergence check, overriding conv;
#             a function(theta_old, theta, fx_old, fx) returning TRUE
#             (converged) or FALSE (continue); used by the h0 EM optimizer
#             to combine the logl criterion with a gradient check
# tol       = convergence tolerance
# max_iter  = maximum number of EM map evaluations (not cycles), so that
#             the em.h1.args$max_iter option keeps the same meaning as for
#             the non-accelerated EM iterations
# step_max0 = initial upper bound for the extrapolation step length
# mstep_factor = growth/shrink factor for the step length bound
# logl_dec  = maximum allowed *decrease* in the loglikelihood for an
#             extrapolated point to be accepted (mild non-monotonicity
#             speeds up convergence; see Varadhan & Roland); if even the
#             plain (fallback) EM step decreases the loglikelihood by more
#             than logl_dec, the iterations stop at the previous point
#             (returned with stalled = TRUE)
lav_em_squarem <- function(theta = NULL,
                           step_fn = NULL,
                           logl_fn = NULL,
                           fx0 = NULL,
                           conv = "logl",
                           conv_fn = NULL,
                           tol = 1e-04,
                           max_iter = 5000L,
                           step_max0 = 1,
                           mstep_factor = 4,
                           logl_dec = 1) {
  if (is.null(fx0)) {
    fx0 <- logl_fn(theta)
  }
  fx <- fx_old <- fx0
  fpeval <- 0L
  cycle <- 0L
  converged <- FALSE
  stalled <- FALSE
  step_max <- step_max0

  while (fpeval < max_iter) {
    cycle <- cycle + 1L

    # two plain EM steps
    p0 <- theta
    p1 <- step_fn(p0)
    fpeval <- fpeval + 1L
    q1 <- p1 - p0
    p2 <- step_fn(p1)
    fpeval <- fpeval + 1L
    q2 <- p2 - p1
    v <- q2 - q1
    sv2 <- sum(v * v)

    # extrapolate: theta' = p0 + 2*alpha*q1 + alpha^2*(q2 - q1)
    # with steplength alpha = sqrt(|q1|^2 / |q2 - q1|^2) ("SqS3"),
    # bounded by step_max; note that alpha = 1 gives theta' = p2, the
    # plain (double) EM step
    alpha <- 1
    if (sv2 > .Machine$double.eps) {
      alpha <- max(1, min(step_max, sqrt(sum(q1 * q1) / sv2)))
    }
    accepted <- FALSE
    if (alpha > 1.01) {
      p_new <- p0 + 2 * alpha * q1 + alpha * alpha * v
      # stabilize: one extra EM step from the extrapolated point; this
      # also repairs (mildly) infeasible extrapolations, as the output
      # of an M-step is always a valid parameter vector
      p_try <- try(step_fn(p_new), silent = TRUE)
      fpeval <- fpeval + 1L
      if (!inherits(p_try, "try-error") &&
          all(is.finite(p_try))) {
        fx_try <- try(logl_fn(p_try), silent = TRUE)
        if (!inherits(fx_try, "try-error") && is.finite(fx_try) &&
            fx_try >= fx_old - logl_dec) {
          theta <- p_try
          fx <- fx_try
          accepted <- TRUE
        }
      }
      if (!accepted) {
        # rejected: shrink the steplength bound (if we hit it)
        if (alpha >= step_max) {
          step_max <- max(step_max0, step_max / mstep_factor)
        }
        alpha <- 1
      }
    }
    if (!accepted) {
      # plain (double) EM step
      theta <- p2
      fx <- logl_fn(theta)
      if (!is.finite(fx)) {
        lav_msg_stop(gettext(
          "The EM iterations failed; some matrices may be singular;
           please check your data for (near-)perfect correlations."))
      }
      if (fx < fx_old - logl_dec) {
        # even the plain (double) EM step decreased the loglikelihood
        # substantially: the map is not behaving like a monotone EM map
        # at this point (this can happen when the M-step is only solved
        # approximately, as in the h0 EM optimizer); revert and stop
        theta <- p0
        fx <- fx_old
        stalled <- TRUE
        break
      }
      accepted <- TRUE
    }
    # if we hit the bound, allow longer steps in the next cycle
    if (accepted && alpha >= step_max) {
      step_max <- mstep_factor * step_max
    }

    fx_delta <- fx - fx_old
    if (lav_verbose()) {
      cat(
        "EM cycle (squarem):", sprintf("%3d", cycle),
        " fx =", sprintf("%17.10f", fx),
        " fx.delta =", sprintf("%9.8f", fx_delta),
        " alpha =", sprintf("%6.2f", alpha),
        "\n"
      )
    }

    # convergence check (per cycle)
    if (!is.null(conv_fn)) {
      converged <- isTRUE(conv_fn(p0, theta, fx_old, fx))
    } else if (conv == "param") {
      converged <- max(abs(theta - p0)) < tol
    } else {
      converged <- abs(fx_delta) < tol
    }
    if (converged) {
      break
    }
    fx_old <- fx
  }

  list(theta = theta, fx = fx, converged = converged, stalled = stalled,
       fpeval = fpeval, cycles = cycle)
}
