# lavSimulate: fit the *same* model, on simulated datasets
# YR - 4 July 2016: initial version
# YR - 15 Oct 2024: add iseed (similar to lavBootstrap)
# YR - 26 Oct 2024: rm pop.model, add est.true argument

lavSimulate <- function(model = NULL, # user model
                        data_function = lav_data_simulate_old,
                        data_function_args = list(),
                        est_true = NULL,
                        ndat = 1000L,
                        cmd = "sem",
                        ...,
                        store_slots = c("partable"),
                        fun = NULL,
                        show_progress = FALSE,
                        store_failed = FALSE,
                        parallel = c("no", "multicore", "snow"),
                        ncpus = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
                        cl = NULL,
                        iseed = NULL) {
  # dotdotdot
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, FALSE)

  # dotdotdot for fit.pop
  dotdotdot_pop <- dotdotdot
  dotdotdot_pop$verbose <- FALSE
  dotdotdot_pop$debug <- FALSE
  dotdotdot_pop$data <- NULL
  dotdotdot_pop$sample_cov <- NULL

  # 'fit' model without data, check 'true' parameters
  fit_pop <- do.call(cmd,
    args = c(list(model = model), dotdotdot_pop)
  )

  # check est.true= argument
  stopifnot(!missing(est_true), is.numeric(est_true),
            length(fit_pop@ParTable$lhs) == length(est_true))
  fit_pop@ParTable$est <- est_true
  fit_pop@ParTable$start <- est_true

  # per default, use 'true' values as starting values
  if (is.null(dotdotdot$start)) {
    dotdotdot$start <- fit_pop
  }

  # honor fixed.x = TRUE: when the exogenous covariates are treated as fixed,
  # the analytic standard errors are *conditional* on these covariates (and
  # their sample (co)variances are fixed to the observed values). For the
  # simulation to match, the exogenous covariates must therefore be identical
  # in every replication. We hold the exogenous data fixed and regenerate the
  # endogenous variables conditionally: each replication is generated with the
  # user's data_function, but the endogenous part is re-expressed at the fixed
  # exogenous design, using the population regression of the (observed)
  # endogenous on the (observed) exogenous variables and the residuals of the
  # freshly generated dataset. This keeps the residual distribution intact
  # (e.g. non-normal) while fixing the exogenous (co)variances.
  fixed_x_flag <- lavInspect(fit_pop, "options")$fixed.x
  ov_x <- lavNames(fit_pop, "ov.x")
  if (isTRUE(fixed_x_flag) && length(ov_x) > 0L &&
      lavInspect(fit_pop, "ngroups") == 1L && is.function(data_function)) {
    ov_all <- lavNames(fit_pop)
    ov_y <- ov_all[!ov_all %in% ov_x]
    if (length(ov_y) > 0L) {
      # reproducible population/design generation
      if (!is.null(iseed)) {
        set.seed(iseed)
      }
      orig_data_function <- data_function
      orig_args <- data_function_args

      # estimate the (raw) population moments by stacking generated datasets
      d_one <- as.data.frame(do.call(orig_data_function, orig_args))
      if (!all(c(ov_x, ov_y) %in% names(d_one))) {
        lav_msg_stop(gettext(
          "fixed.x = TRUE: the data_function output must contain all observed
           variables of the model (as named columns)."))
      }
      n_one <- nrow(d_one)
      nrep_pop <- max(1L, min(200L, ceiling(1e5 / n_one)))
      if (nrep_pop > 1L) {
        d_more <- lapply(seq_len(nrep_pop - 1L), function(b) {
          as.data.frame(do.call(orig_data_function, orig_args))
        })
        d_big <- do.call("rbind", c(list(d_one), d_more))
      } else {
        d_big <- d_one
      }
      mu_pop <- colMeans(d_big[, ov_all, drop = FALSE])
      s_pop <- stats::cov(d_big[, ov_all, drop = FALSE])
      xi <- match(ov_x, ov_all)
      yi <- match(ov_y, ov_all)
      beta_yx <- s_pop[yi, xi, drop = FALSE] %*%
        solve(s_pop[xi, xi, drop = FALSE])
      mu_x <- mu_pop[xi]
      mu_y <- mu_pop[yi]

      # the fixed exogenous design (identical in every replication). Recolor
      # it so that its (co)variances match the population exogenous moments
      # exactly; otherwise this single realization would differ from the
      # population by O(1/sqrt(n)), and the conditional 'true' values would
      # not match est.true (which is defined at the population level).
      x_raw <- as.matrix(d_one[, ov_x, drop = FALSE])
      sx_pop <- s_pop[xi, xi, drop = FALSE]
      w_x <- solve(chol(stats::cov(x_raw))) %*% chol(sx_pop)
      x_fixed <- sweep(sweep(x_raw, 2L, colMeans(x_raw)) %*% w_x,
                       2L, mu_x, "+")
      colnames(x_fixed) <- ov_x

      data_function <- function(...) {
        d <- as.data.frame(do.call(orig_data_function, orig_args))
        if (nrow(d) != nrow(x_fixed)) {
          lav_msg_stop(gettext(
            "fixed.x = TRUE: the data_function must return a constant number
             of observations across replications."))
        }
        x_r <- as.matrix(d[, ov_x, drop = FALSE])
        y_r <- as.matrix(d[, ov_y, drop = FALSE])
        # residuals of the freshly generated endogenous data (under the
        # population regression)
        e_r <- sweep(y_r, 2L, mu_y) - sweep(x_r, 2L, mu_x) %*% t(beta_yx)
        # endogenous, conditional on the *fixed* exogenous design
        y_new <- sweep(sweep(x_fixed, 2L, mu_x) %*% t(beta_yx),
                       2L, mu_y, "+") + e_r
        d[, ov_y] <- y_new
        d[, ov_x] <- x_fixed
        d
      }
      data_function_args <- list() # original args are captured in the closure
    }
  }

  # no warnings during/after the simulations
  # add 'warn = FALSE' to args

  # generate simulations
  fit <- do.call("lavaanList", args = c(list(
    model = model,
    dataFunction = data_function,
    dataFunction.args = data_function_args,
    ndat = ndat, cmd = cmd,
    store.slots = store_slots, FUN = fun,
    show.progress = show_progress,
    store.failed = store_failed,
    parallel = parallel, ncpus = ncpus,
    cl = cl, iseed = iseed
  ), dotdotdot))

  # flag this is a simulation
  fit@meta$lavSimulate <- TRUE

  # store 'true' parameters in meta$est.true
  fit@meta$est.true <- est_true

  fit
}
