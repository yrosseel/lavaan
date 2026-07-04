lav_step11_estoptim <- function(lavdata = NULL,
                                       lavmodel = NULL,
                                       lavcache = NULL,
                                       lavsamplestats = NULL,
                                       lavh1 = NULL,
                                       lavoptions = NULL,
                                       lavpartable = NULL) {
  # # # # # # # # # # # # # #
  # #  11. est + lavoptim # #
  # # # # # # # # # # # # # #

  # if lavoptions$do.fit and lavoptions$estimator not "none" and
  #   lavmodel$nx.free > 0
  #   select case lavoptions$optim.method
  #   case "noniter"
  #     try x <- lav_optim_noniter(...)
  #   case "em"
  #     if nlevels < 2L *** error ***
  #     try x <- lav_mvn_cl_em_h0(...)
  #   case "gn"
  #     try x <- lav_optim_gn(...)
  #   case else
  #     set 1 in lavoptions$optim.attempts if it wasn't specified
  #     try x <- lav_model_est(...)
  #     if not successful and optim.attempts > 1L
  #       try x <- lav_optim_estimate(...) with
  #         options$optim.parscale = "standardized"
  #       if not successful and optim.attempts > 2L
  #         try x <- lav_optim_estimate(...) with start = "simple"
  #         if not successful and optim.attempts > 3L
  #           try x <- lav_optim_estimate(...) with
  #             options$optim.parscale = "standardized" and start = "simple"
  #   end select
  #   if x not successfully computed
  #     ** warning **
  #     set starting values and appropriate attributes in x
  #   in case of non-linear constraints: store final con.jac and
  #     con.lambda in lavmodel
  #   store parameters in lavmodel
  #   store parameters in partable$est
  # else
  #   initialize x and attributes (iterations, converged, warn.txt,
  #     control, dx) of x
  #   try fx <- lav_model_objective
  #   if not successful
  #     fx = NA_real_
  #     attribute fx.group of fx = NA_real_
  #   store fx in attribute "fx" of x
  #   set lavpartable$est to starting values
  # if lavoptions$optim.force.converged set attribute converged of x to TRUE
  # store optimization info in lavoptim

  x <- NULL
  if (lavoptions$do.fit && lavoptions$estimator != "none" &&
    lavmodel@nx.free > 0L) {
    if (lav_verbose()) {
      cat("lavoptim           ... start:\n")
    }

    # the EM optimizer only supports two-level models with a single
    # group, an ML-family estimator, no conditional.x, and no
    # (nonlinear or inequality) constraints; if optim.method = "em" was
    # merely the default (missing = "ml"; see lav_options_set), quietly
    # fall back to "nlminb"; if the user asked for it explicitly, warn
    if (lavoptions$optim.method == "em" &&
        (lavdata@nlevels == 1L ||
         lavdata@ngroups > 1L ||
         lavmodel@conditional.x ||
         lavmodel@estimator != "ML" ||
         length(lavmodel@ceq.nonlinear.idx) > 0L ||
         length(lavmodel@cin.linear.idx) > 0L ||
         length(lavmodel@cin.nonlinear.idx) > 0L)) {
      if (!isTRUE(lavoptions$.optim.em.fallback)) {
        lav_msg_warn(gettext(
          "optim.method = \"em\" is not available for this model;
           using optim.method = \"nlminb\" instead."))
      }
      lavoptions$optim.method <- "nlminb"
    }

    # non-iterative methods (fabin, miiv, ...)
    if (lavoptions$optim.method == "noniter") {
      x <- try(
        lav_optim_noniter(
          lavmodel = lavmodel,
          lavsamplestats = lavsamplestats,
          lavh1 = lavh1,
          lavdata = lavdata,
          lavpartable = lavpartable,
          lavoptions = lavoptions
        ),
        silent = FALSE
      )
      # EM for multilevel models
    } else if (lavoptions$optim.method == "em") {
      # multilevel only for now
      stopifnot(lavdata@nlevels > 1L)
      x <- try(
        lav_mvn_cl_em_h0(
          lavsamplestats = lavsamplestats,
          lavdata = lavdata,
          lavimplied = NULL,
          lavpartable = lavpartable,
          lavmodel = lavmodel,
          lavoptions = lavoptions,
          fx_tol = lavoptions$em.args$fx_tol,
          dx_tol = lavoptions$em.args$dx_tol,
          max_iter = lavoptions$em.args$iter_max,
          acceleration = lavoptions$em.args$acceleration,
          fused = lavoptions$em.args$fused
        ),
        silent = TRUE
      )

      # EM-stall -> nlminb hand-off (new in 0.7-2): the EM iterations
      # typically stop because the loglikelihood stalls (fx_tol) before
      # the gradient becomes small (dx_tol); in that case, refine the
      # solution with nlminb, warm-started at the EM values; if nlminb
      # gets lost in parameter space (non-finite values, a worse
      # solution) or does not converge, we keep the EM solution
      if (!inherits(x, "try-error") &&
          !isFALSE(lavoptions$em.args$nlminb_handoff)) {
        em_dx_tol <- lavoptions$em.args$dx_tol
        if (!is.numeric(em_dx_tol)) {
          em_dx_tol <- 1e-04
        }
        lavmodel_em <- lav_model_set_parameters(lavmodel,
          x = as.numeric(x))
        dx_em <- try(lav_model_grad(
          lavmodel = lavmodel_em,
          lavsamplestats = lavsamplestats,
          lavdata = lavdata
        ), silent = TRUE)
        if (!inherits(dx_em, "try-error") && all(is.finite(dx_em)) &&
            max(abs(dx_em)) >= em_dx_tol) {
          if (lav_verbose()) {
            cat("EM stalled (max|grad| = ", max(abs(dx_em)),
                "); handing off to nlminb\n")
          }
          # objective at the EM solution (for the accept/reject check)
          fx_em <- try(lav_model_objective(
            lavmodel = lavmodel_em,
            lavsamplestats = lavsamplestats,
            lavdata = lavdata,
            lavcache = lavcache
          ), silent = TRUE)
          lavoptions2 <- lavoptions
          lavoptions2$optim.method <- "nlminb"
          # do NOT pin saturated blocks at their h1 values: the EM
          # values are (at least as) good, and pinning would overwrite
          # the warm start
          lavoptions2$optim.fix.saturated <- FALSE
          x2 <- try(
            lav_model_est(
              lavmodel = lavmodel_em, # warm start at the EM solution
              lavpartable = lavpartable,
              lavh1 = lavh1,
              lavsamplestats = lavsamplestats,
              lavdata = lavdata,
              lavoptions = lavoptions2,
              lavcache = lavcache
            ),
            silent = TRUE
          )
          accept <- (!inherits(x2, "try-error") &&
            isTRUE(attr(x2, "converged")) &&
            all(is.finite(x2)) &&
            is.finite(attr(x2, "fx")[1]) &&
            (inherits(fx_em, "try-error") || !is.finite(fx_em[1]) ||
             attr(x2, "fx")[1] <= fx_em[1] + 1e-10))
          if (accept) {
            x <- x2
          }
          if (lav_verbose()) {
            if (accept) {
              cat("nlminb hand-off accepted (fx = ",
                  attr(x2, "fx")[1], ")\n")
            } else {
              cat("nlminb hand-off rejected; keeping the EM solution\n")
            }
          }
        }
      }

      # Gauss-Newton
    } else if (lavoptions$optim.method == "gn") {
      # only tested for DLS (for now)
      x <- try(
        lav_optim_gn(
          lavmodel = lavmodel,
          lavsamplestats = lavsamplestats,
          lavdata = lavdata,
          lavpartable = lavpartable,
          lavoptions = lavoptions
        ),
        silent = TRUE
      )

      # reduced-bias M-estimation (RBM); optim.method is nlminb, so this is
      # keyed on the estimator.args marker instead
    } else if (!is.null(lavoptions$estimator.args$rbm.method)) {
      x <- try(
        lav_model_est_rbm(
          lavmodel = lavmodel,
          lavpartable = lavpartable,
          lavsamplestats = lavsamplestats,
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavcache = lavcache
        ),
        silent = TRUE
      )

      # Quasi-Newton
    } else {
      # for backwards compatibility (<0.6)
      if (is.null(lavoptions$optim.attempts)) {
        lavoptions$optim.attempts <- 1L
      }

      # try 1
      if (lav_verbose()) {
        cat("attempt 1 -- default options\n")
      }
      x <- try(
        lav_model_est(
          lavmodel = lavmodel,
          lavpartable = lavpartable,
          lavh1 = lavh1,
          lavsamplestats = lavsamplestats,
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavcache = lavcache
        ),
        silent = TRUE
      )

      # try 2: optim.parscale = "standardize" (new in 0.6-7)
      if (lavoptions$optim.attempts > 1L &&
        lavoptions$rstarts == 0L &&
        (inherits(x, "try-error") || !attr(x, "converged"))) {
        lavoptions2 <- lavoptions
        lavoptions2$optim.parscale <- "standardized"
        if (lav_verbose()) {
          str(x)
          cat("attempt 2 -- optim.parscale = \"standardized\"\n")
        }
        x <- try(
          lav_model_est(
            lavmodel = lavmodel,
            lavpartable = lavpartable,
            lavh1 = lavh1,
            lavsamplestats = lavsamplestats,
            lavdata = lavdata,
            lavoptions = lavoptions2,
            lavcache = lavcache
          ),
          silent = TRUE
        )
      }

      # try 3: start = "simple"
      if (lavoptions$optim.attempts > 2L &&
        lavoptions$rstarts == 0L &&
        (inherits(x, "try-error") || !attr(x, "converged"))) {
        if (lav_verbose()) {
          str(x)
          cat("attempt 3 -- start = \"simple\"\n")
        }
        x <- try(
          lav_model_est(
            lavmodel = lavmodel,
            lavpartable = lavpartable,
            lavh1 = lavh1,
            lavsamplestats = lavsamplestats,
            lavdata = lavdata,
            lavoptions = lavoptions,
            start = "simple",
            lavcache = lavcache
          ),
          silent = TRUE
        )
      }

      # try 4: start = "simple" + optim.parscale = "standardize"
      if (lavoptions$optim.attempts > 3L &&
        lavoptions$rstarts == 0L &&
        (inherits(x, "try-error") || !attr(x, "converged"))) {
        lavoptions2 <- lavoptions
        lavoptions2$optim.parscale <- "standardized"
        if (lav_verbose()) {
          str(x)
          cat(
            "attempt 4 -- optim.parscale = \"standardized\" + ",
            "start = \"simple\"\n"
          )
        }
        x <- try(
          lav_model_est(
            lavmodel = lavmodel,
            lavpartable = lavpartable,
            lavh1 = lavh1,
            lavsamplestats = lavsamplestats,
            lavdata = lavdata,
            lavoptions = lavoptions2,
            start = "simple",
            lavcache = lavcache
          ),
          silent = TRUE
        )
      }


      # random starts? -- new in 0.6-18
      # run this even if we already have a converged solution
      # perhaps we find a better solution?
      if (lavoptions$rstarts > 0L) {
        x_rstarts <- vector("list", length = lavoptions$rstarts)
        if (lav_verbose()) {
          str(x)
          cat("trying again with random starts (", lavoptions$rstarts,
            " in total):\n",
            sep = ""
          )
        }
        for (i in seq_len(lavoptions$rstarts)) {
          if (lav_verbose()) {
            cat("-- random start run: ", i, "\n")
          }
          x_rstarts[[i]] <-
            try(
              lav_model_est(
                lavmodel = lavmodel,
                lavpartable = lavpartable,
                lavh1 = lavh1,
                lavsamplestats = lavsamplestats,
                lavdata = lavdata,
                lavoptions = lavoptions,
                start = "random",
                lavcache = lavcache
              ),
              silent = TRUE
            )
        }

        # pick best solution (if any)
        x_converged <- vector("list", length = 0L)
        fx_rstarts <- numeric(0L)
        ok_flag <- sapply(x_rstarts, function(x) {
            if (inherits(x, "try-error")) {
              FALSE
            } else {
              attr(x, "converged")
            }
          })
        if (sum(ok_flag) > 0L) {
          x_converged <- x_rstarts[ok_flag]
        }
        if (length(x_converged) > 0L) {
          fx_rstarts <- sapply(x_converged, "attr", "fx")
          x_best <- x_converged[[which.min(fx_rstarts)]]
          fx_best <- attr(x_best, "fx")[1]

          # if we did not find a converged solution, use x.best
          if (inherits(x, "try-error") || !attr(x, "converged")) {
            x <- x_best


            # if we already had a converged solution, only replace
            # if fx.best is better than attr(x, "fx")[1]
          } else {
            if (fx_best < attr(x, "fx")[1]) {
              x <- x_best
            }
          }
        }

        attr(x, "x.rstarts") <- x_rstarts
      } # random starts
    }

    # optimization failed with error
    if (inherits(x, "try-error")) {
      warn_txt <- gettext("Model estimation FAILED! Returning starting values.")
      x <- lav_model_get_parameters(
        lavmodel = lavmodel,
        type = "free"
      ) # starting values
      attr(x, "iterations") <- 0L
      attr(x, "converged") <- FALSE
      attr(x, "warn.txt") <- warn_txt
      attr(x, "control") <- lavoptions$control
      attr(x, "dx") <- numeric(0L)
      fx <- as.numeric(NA)
      attr(fx, "fx.group") <- as.numeric(NA)
      attr(x, "fx") <- fx
    }

    # if a warning was produced, say it here
    warn_txt <- attr(x, "warn.txt")
    if (!is.null(warn_txt) && nchar(warn_txt) > 0L) {
      lav_msg_warn(gettext(warn_txt))
    }

    # in case of non-linear constraints: store final con.jac and con.lambda
    # in lavmodel
    if (!is.null(attr(x, "con.jac"))) {
      lavmodel@con.jac <- attr(x, "con.jac")
    }
    if (!is.null(attr(x, "con.lambda"))) {
      lavmodel@con.lambda <- attr(x, "con.lambda")
    }

    # store parameters in lavmodel
    lavmodel <- lav_model_set_parameters(lavmodel, x = as.numeric(x))

    # store parameters in @ParTable$est
    lavpartable$est <- lav_model_get_parameters(
      lavmodel = lavmodel,
      type = "user", extra = TRUE
    )

    if (lav_verbose()) {
      cat("lavoptim    ... done.\n")
    }
  } else {
    x <- numeric(0L)
    attr(x, "iterations") <- 0L
    attr(x, "converged") <- FALSE
    attr(x, "warn.txt") <- ""
    attr(x, "control") <- lavoptions$control
    attr(x, "dx") <- numeric(0L)
    fx <- try(lav_model_objective(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavcache = lavcache
    ), silent = TRUE)
    if (!inherits(fx, "try-error")) {
      attr(x, "fx") <- fx
    } else {
      fx <- as.numeric(NA)
      attr(fx, "fx.group") <- as.numeric(NA)
      attr(x, "fx") <- fx
    }

    lavpartable$est <- lavpartable$start
  }

  # should we fake/force convergence? (eg. to enforce the
  # computation of a test statistic)
  if (lavoptions$optim.force.converged) {
    attr(x, "converged") <- TRUE
  }

  # store optimization info in lavoptim
  lavoptim <- list()
  x2 <- x
  attributes(x2) <- NULL
  lavoptim$x <- x2
  lavoptim$dx <- attr(x, "dx")
  lavoptim$npar <- length(x)
  lavoptim$iterations <- attr(x, "iterations")
  lavoptim$converged <- attr(x, "converged")
  lavoptim$warn.txt <- attr(x, "warn.txt")
  lavoptim$parscale <- attr(x, "parscale")
  lavoptim$partrace <- attr(x, "partrace")
  fx_copy <- fx <- attr(x, "fx")
  attributes(fx) <- NULL
  lavoptim$fx <- fx
  lavoptim$fx.group <- attr(fx_copy, "fx.group")
  if (!is.null(attr(fx_copy, "logl.group"))) {
    lavoptim$logl.group <- attr(fx_copy, "logl.group")
    lavoptim$logl <- sum(lavoptim$logl.group)
  } else {
    lavoptim$logl.group <- as.numeric(NA)
    lavoptim$logl <- as.numeric(NA)
  }
  lavoptim$control <- attr(x, "control")
  if (!is.null(attr(x, "x.rstarts"))) {
    lavoptim$x.rstarts <- attr(x, "x.rstarts")
  }

  lavpartable <- lav_pt_set_cache(lavpartable, force = TRUE)
  list(
    lavoptim = lavoptim, lavmodel = lavmodel, lavpartable = lavpartable,
    x = x
  )
}
