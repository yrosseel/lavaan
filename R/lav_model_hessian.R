# numeric approximation of the Hessian
# using an analytic gradient
#
# method = "richardson" (the default): 4-point Richardson
#   extrapolation (4 x npar gradient evaluations), with a fixed
#   (absolute) step size h
# method = "central": plain central differences (2 x npar gradient
#   evaluations), with a *relative* step size h * max(1, |x_j|);
#   cheaper, and the appropriate choice when the gradient itself
#   contains a numerical layer (eg random-slope models)
#
# `NA' semantics: if the gradient cannot be computed (or is not
# finite) at one of the perturbed parameter values, an all-NA matrix
# is returned (callers can check with anyNA()), rather than a hard
# stop or a silently wrong Hessian
lav_model_hessian <- function(lavmodel = NULL,
                              lavsamplestats = NULL,
                              lavdata = NULL,
                              lavoptions = NULL,
                              lavcache = NULL,
                              group_weight = TRUE,
                              ceq_simple = FALSE,
                              h = 1e-06,
                              method = "richardson") {
  # estimator <- lavmodel@estimator

  # catch numerical gradient
  if (lavoptions$optim.gradient == "numerical") {
    obj_f <- function(x) {
      lavmodel2 <- lav_model_set_parameters(lavmodel, x = x)
      lav_model_objective(
        lavmodel = lavmodel2,
        lavsamplestats = lavsamplestats, lavdata = lavdata,
        lavcache = lavcache
      )[1]
    }
    x <- lav_model_get_parameters(lavmodel = lavmodel)
    hessian <- try(numDeriv::hessian(func = obj_f, x = x),
                   silent = TRUE)
    if (inherits(hessian, "try-error")) {
      hessian <- matrix(as.numeric(NA), length(x), length(x))
    }
    return(hessian)
  }

  # two-level + missing data: use Louis's (1982) method (one analytic
  # pass) instead of numerically differentiating the gradient (which
  # needs 4 x npar full-data gradient passes); new in 0.7-1
  if (!is.null(lavdata) && lavdata@nlevels > 1L &&
      lavdata@ngroups == 1L &&
      lavsamplestats@missing.flag &&
      lavmodel@estimator == "ML" &&
      !lavmodel@conditional.x &&
      length(lavmodel@rv.ov) == 0L && length(lavmodel@rv.lv) == 0L &&
      is.null(lavdata@weights[[1]])) {
    hessian <- try(
      lav_mvn_cl_mi_h_louis(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        ceq_simple = ceq_simple
      ),
      silent = TRUE
    )
    if (!inherits(hessian, "try-error") && all(is.finite(hessian))) {
      return(hessian)
    }
    # if the Louis computation fails, fall through to the numerical
    # approximation below
  }

  # computing the Richardson extrapolation (or central differences)
  if (!ceq_simple && lavmodel@ceq.simple.only) {
    npar <- lavmodel@nx.unco
    type_glist <- "unco"
  } else {
    npar <- lavmodel@nx.free
    type_glist <- "free"
  }
  hessian <- matrix(0, npar, npar)
  x <- lav_model_get_parameters(lavmodel = lavmodel)
  if (!ceq_simple && lavmodel@ceq.simple.only) {
    # unpack
    x <- drop(x %*% t(lavmodel@ceq.simple.K))
  }

  # the gradient at a perturbed x; NULL if it cannot be computed
  grad_at <- function(x_p) {
    out <- try(
      lav_model_grad(
        lavmodel = lavmodel,
        glist = lav_model_x2glist(
          lavmodel = lavmodel, type = type_glist, x = x_p
        ),
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        type = "free",
        group_weight = group_weight,
        ceq_simple = ceq_simple
      ),
      silent = TRUE
    )
    if (inherits(out, "try-error") || anyNA(out) ||
        any(!is.finite(out))) {
      return(NULL)
    }
    out
  }
  hessian_na <- matrix(as.numeric(NA), npar, npar)

  for (j in seq_len(npar)) {
    if (method == "central") {
      # central differences, relative step size
      h_j <- h * max(1, abs(x[j]))
      x_left <- x_right <- x
      x_left[j] <- x[j] - h_j
      x_right[j] <- x[j] + h_j

      g_left <- grad_at(x_left)
      g_right <- grad_at(x_right)
      if (is.null(g_left) || is.null(g_right)) {
        return(hessian_na)
      }

      hessian[, j] <- (g_right - g_left) / (2 * h_j)
    } else {
      # Richardson extrapolation
      # FIXME: the number below should vary as a function of 'x[j]'
      h_j <- h
      x_left <- x_left2 <- x_right <- x_right2 <- x
      x_left[j] <- x[j] - h_j
      x_left2[j] <- x[j] - 2 * h_j
      x_right[j] <- x[j] + h_j
      x_right2[j] <- x[j] + 2 * h_j

      g_left <- grad_at(x_left)
      g_left2 <- grad_at(x_left2)
      g_right <- grad_at(x_right)
      g_right2 <- grad_at(x_right2)
      if (is.null(g_left) || is.null(g_left2) ||
          is.null(g_right) || is.null(g_right2)) {
        return(hessian_na)
      }

      hessian[, j] <- (g_left2 - 8 * g_left + 8 * g_right - g_right2) /
        (12 * h_j)
    }
  }

  # check if Hessian is (almost) symmetric, as it should be
  max_diff <- max(abs(hessian - t(hessian)))
  if (max_diff > 1e-05 * max(diag(hessian))) {
    # hm, Hessian is not symmetric -> WARNING!
    lav_msg_warn(gettextf(
    "Hessian is not fully symmetric. Max diff = %1$s (Max diag Hessian = %2$s)",
    max_diff, max(diag(hessian))))
    # FIXME: use numDeriv::hessian instead?
  }
  hessian <- (hessian + t(hessian)) / 2.0

  hessian
}
