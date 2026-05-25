# numeric approximation of the Hessian
# using an analytic gradient
lav_model_hessian <- function(lavmodel = NULL,
                              lavsamplestats = NULL,
                              lavdata = NULL,
                              lavoptions = NULL,
                              lavcache = NULL,
                              group_weight = TRUE,
                              ceq_simple = FALSE,
                              h = 1e-06) {
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
    hessian <- numDeriv::hessian(func = obj_f, x = x)
    return(hessian)
  }

  # computing the Richardson extrapolation
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
  for (j in seq_len(npar)) {
    # FIXME: the number below should vary as a function of 'x[j]'
    h_j <- h
    x_left <- x_left2 <- x_right <- x_right2 <- x
    x_left[j] <- x[j] - h_j
    x_left2[j] <- x[j] - 2 * h_j
    x_right[j] <- x[j] + h_j
    x_right2[j] <- x[j] + 2 * h_j

    g_left <-
      lav_model_grad(
        lavmodel = lavmodel,
        glist = lav_model_x2glist(
          lavmodel =
            lavmodel, type = type_glist,
          x_left
        ),
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        type = "free",
        group_weight = group_weight,
        ceq_simple = ceq_simple
      )
    g_left2 <-
      lav_model_grad(
        lavmodel = lavmodel,
        glist = lav_model_x2glist(
          lavmodel =
            lavmodel, type = type_glist,
          x_left2
        ),
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        type = "free",
        group_weight = group_weight,
        ceq_simple = ceq_simple
      )

    g_right <-
      lav_model_grad(
        lavmodel = lavmodel,
        glist = lav_model_x2glist(
          lavmodel =
            lavmodel, type = type_glist,
          x_right
        ),
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        type = "free",
        group_weight = group_weight,
        ceq_simple = ceq_simple
      )

    g_right2 <-
      lav_model_grad(
        lavmodel = lavmodel,
        glist = lav_model_x2glist(
          lavmodel =
            lavmodel, type = type_glist,
          x_right2
        ),
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        type = "free",
        group_weight = group_weight,
        ceq_simple = ceq_simple
      )

    hessian[, j] <- (g_left2 - 8 * g_left + 8 * g_right - g_right2) / (12 * h_j)
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
