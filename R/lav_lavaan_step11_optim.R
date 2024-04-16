lav_lavaan_step11_estoptim <- function(lavdata = NULL, # nolint
                                       lavmodel = NULL,
                                       lavcache = NULL,
                                       lavsamplestats = NULL,
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
  #     try x <- lav_mvnorm_cluster_em_h0(...)
  #   case "gn"
  #     try x <- lav_optim_gn(...)
  #   case else
  #     set 1 in lavoptions$optim.attempts is it wasn't specified
  #     try x <- lav_model_estimate(...)
  #     if not successfull and optim.attempts > 1L
  #       try x <- lav_optim_estimate(...) with
  #         options$optim.parscale = "standardized"
  #       if not successfull and optim.attempts > 2L
  #         try x <- lav_optim_estimate(...) with start = "simple"
  #         if not successfull and optim.attempts > 3L
  #           try x <- lav_optim_estimate(...) with
  #             options$optim.parscale = "standardized" and start = "simple"
  #   end select
  #   if x not succesfully computed
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
  #   if not successfull
  #     fx = NA_real_
  #     attribute fx.group of fx = NA_real_
  #   store fx in attribute "fx" of x
  #   set lavpartable$est to starting values
  # if lavoptions$optim.force.converged set attribute converged of x to TRUE
  # store optimization info in lavoptim

  x <- NULL
  if (lavoptions$do.fit && lavoptions$estimator != "none" &&
    lavmodel@nx.free > 0L) {
    if (lavoptions$verbose) {
      cat("lavoptim           ... start:\n")
    }

    # non-iterative methods (fabin, ...)
    if (lavoptions$optim.method == "noniter") {
      x <- try(
        lav_optim_noniter(
          lavmodel = lavmodel,
          lavsamplestats = lavsamplestats,
          lavdata = lavdata,
          lavpartable = lavpartable,
          lavoptions = lavoptions
        ),
        silent = TRUE
      )
      # EM for multilevel models
    } else if (lavoptions$optim.method == "em") {
      # multilevel only for now
      stopifnot(lavdata@nlevels > 1L)
      x <- try(
        lav_mvnorm_cluster_em_h0(
          lavsamplestats = lavsamplestats,
          lavdata = lavdata,
          lavimplied = NULL,
          lavpartable = lavpartable,
          lavmodel = lavmodel,
          lavoptions = lavoptions,
          verbose = lavoptions$verbose,
          fx.tol = lavoptions$em.fx.tol,
          dx.tol = lavoptions$em.dx.tol,
          max.iter = lavoptions$em.iter.max
        ),
        silent = TRUE
      )
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

      # Quasi-Newton
    } else {
      # for backwards compatibility (<0.6)
      if (is.null(lavoptions$optim.attempts)) {
        lavoptions$optim.attempts <- 1L
      }

      # try 1
      if (lavoptions$verbose) {
        cat("attempt 1 -- default options\n")
      }
      x <- try(
        lav_model_estimate(
          lavmodel = lavmodel,
          lavpartable = lavpartable,
          lavsamplestats = lavsamplestats,
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavcache = lavcache
        ),
        silent = TRUE
      )
      # store first attempt
      # x.first <- x

      # try 2: optim.parscale = "standardize" (new in 0.6-7)
      if (lavoptions$optim.attempts > 1L &&
        lavoptions$rstarts == 0L &&
        (inherits(x, "try-error") || !attr(x, "converged"))) {
        lavoptions2 <- lavoptions
        lavoptions2$optim.parscale <- "standardized"
        if (lavoptions$verbose) {
          cat("attempt 2 -- optim.parscale = \"standardized\"\n")
        }
        x <- try(
          lav_model_estimate(
            lavmodel = lavmodel,
            lavpartable = lavpartable,
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
        if (lavoptions$verbose) {
          cat("attempt 3 -- start = \"simple\"\n")
        }
        x <- try(
          lav_model_estimate(
            lavmodel = lavmodel,
            lavpartable = lavpartable,
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
        if (lavoptions$verbose) {
          cat(
            "attempt 4 -- optim.parscale = \"standardized\" + ",
            "start = \"simple\"\n"
          )
        }
        x <- try(
          lav_model_estimate(
            lavmodel = lavmodel,
            lavpartable = lavpartable,
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
        x.rstarts <- vector("list", length = lavoptions$rstarts)
        if (lavoptions$verbose) {
          cat("trying again with random starts (", lavoptions$rstarts,
            " in total):\n",
            sep = ""
          )
        }
        for (i in seq_len(lavoptions$rstarts)) {
          if (lavoptions$verbose) {
            cat("-- random start run: ", i, "\n")
          }
          x.rstarts[[i]] <-
            try(
              lav_model_estimate(
                lavmodel = lavmodel,
                lavpartable = lavpartable,
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
        x.noerror <- x.rstarts[!sapply(
          x.rstarts,
          inherits, "try-error"
        )]
        x.converged <- vector("list", length = 0L)
        fx.rstarts <- numeric(0L)
        if (length(x.noerror) > 0L) {
          x.converged <-
            x.noerror[sapply(x.rstarts, "attr", "converged")]
        }
        if (length(x.converged) > 0L) {
          fx.rstarts <- sapply(x.converged, "attr", "fx")
          x.best <- x.converged[[which.min(fx.rstarts)]]
          fx.best <- attr(x.best, "fx")[1]

          # if we did not find a converged solution, use x.best
          if (inherits(x, "try-error") || !attr(x, "converged")) {
            x <- x.best


            # if we already had a converged solution, only replace
            # if fx.best is better than attr(x, "fx")[1]
          } else {
            if (fx.best < attr(x, "fx")[1]) {
              x <- x.best
            }
          }
        }

        attr(x, "x.rstarts") <- x.rstarts
      } # random starts
    }

    # optimization failed with error
    if (inherits(x, "try-error")) {
      warn.txt <- gettext("Model estimation FAILED! Returning starting values.")
      x <- lav_model_get_parameters(
        lavmodel = lavmodel,
        type = "free"
      ) # starting values
      attr(x, "iterations") <- 0L
      attr(x, "converged") <- FALSE
      attr(x, "warn.txt") <- warn.txt
      attr(x, "control") <- lavoptions$control
      attr(x, "dx") <- numeric(0L)
      fx <- as.numeric(NA)
      attr(fx, "fx.group") <- as.numeric(NA)
      attr(x, "fx") <- fx
    }

    # if a warning was produced, say it here
    warn.txt <- attr(x, "warn.txt")
    if (lavoptions$warn && nchar(warn.txt) > 0L) {
      lav_msg_warn(
        gettext("Model estimation FAILED! Returning starting values."))
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

    if (lavoptions$verbose) {
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
  fx.copy <- fx <- attr(x, "fx")
  attributes(fx) <- NULL
  lavoptim$fx <- fx
  lavoptim$fx.group <- attr(fx.copy, "fx.group")
  if (!is.null(attr(fx.copy, "logl.group"))) {
    lavoptim$logl.group <- attr(fx.copy, "logl.group")
    lavoptim$logl <- sum(lavoptim$logl.group)
  } else {
    lavoptim$logl.group <- as.numeric(NA)
    lavoptim$logl <- as.numeric(NA)
  }
  lavoptim$control <- attr(x, "control")
  if (!is.null(attr(x, "x.rstarts"))) {
    lavoptim$x.rstarts <- attr(x, "x.rstarts")
  }

  lavpartable <- lav_partable_set_cache(lavpartable, force = TRUE)
  list(
    lavoptim = lavoptim, lavmodel = lavmodel, lavpartable = lavpartable,
    x = x
  )
}
