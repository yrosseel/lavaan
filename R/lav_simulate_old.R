# simulate data starting from a user-specified model
#
# initial version: YR 24 jan 2011
# revision for 0.4-11: YR 21 okt 2011
#
#
#
lav_data_simulate_old <- function( # user-specified model    # nolint start
                         model = NULL,
                         model.type = "sem",
                         # model modifiers
                         meanstructure = FALSE,
                         int.ov.free = TRUE,
                         int.lv.free = FALSE,
                         marker.int.zero = FALSE,
                         conditional.x = FALSE,
                         composites = TRUE,
                         fixed.x = FALSE,
                         orthogonal = FALSE,
                         std.lv = TRUE,
                         auto.fix.first = FALSE,
                         auto.fix.single = FALSE,
                         auto.var = TRUE,
                         auto.cov.lv.x = TRUE,
                         auto.cov.y = TRUE,
                         ...,
                         # data properties
                         sample.nobs = 500L,
                         ov.var = NULL,
                         group.label = paste("G", 1:ngroups, sep = ""),
                         skewness = NULL,
                         kurtosis = NULL,
                         # control
                         seed = NULL,
                         empirical = FALSE,
                         mass = FALSE,
                         return.type = "data.frame",
                         return.fit = FALSE,
                         debug = FALSE,
                         standardized = FALSE) {           # nolint end
  if (!missing(debug)) {
    current_debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current_debug), TRUE)
  }
  if (!is.null(seed)) set.seed(seed)
  # if(!exists(".Random.seed", envir = .GlobalEnv))
  #    runif(1)               # initialize the RNG if necessary
  # RNGstate <- .Random.seed

  # lav_model_pt
  if (is.list(model)) {
    # two possibilities: either model is already lavaanified
    # or it is something else...
    if (!is.null(model$lhs) && !is.null(model$op) &&
      !is.null(model$rhs) && !is.null(model$free)) {
      lav <- model

      # until 0.6-5, we only used the 'ustart' column
      # but what if 'lav' is a fitted lavaan object -> use 'est'
      if (!is.null(lav$est)) {
        lav$ustart <- lav$est
        lav$se <- NULL
        lav$est <- NULL
        lav$start <- NULL
      }
    } else if (is.character(model[[1]])) {
      lav_msg_stop(gettext("model is a list, but not a parameterTable?"))
    }
  } else {
    lav <- lav_model_pt(
      model = model,
      meanstructure = meanstructure,
      int_ov_free = int.ov.free,
      int_lv_free = int.lv.free,
      marker_int_zero = marker.int.zero,
      composites = composites,
      conditional_x = conditional.x,
      fixed_x = fixed.x,
      orthogonal = orthogonal,
      std_lv = std.lv,
      auto_fix_first = auto.fix.first,
      auto_fix_single = auto.fix.single,
      auto_var = auto.var,
      auto_cov_lv_x = auto.cov.lv.x,
      auto_cov_y = auto.cov.y,
      ngroups = length(sample.nobs)
    )
  }

  group_values <- lav_pt_group_values(lav)
  if (lav_debug()) {
    cat("initial lav\n")
    print(as.data.frame(lav))
  }

  # fill in any remaining NA values (needed for unstandardize)
  # 1 for variances and (unstandardized) factor loadings, 0 otherwise
  idx <- which(lav$op == "=~" & is.na(lav$ustart))
  if (length(idx) > 0L) {
    if (standardized) {
      lav$ustart[idx] <- 0.7
    } else {
      lav$ustart[idx] <- 1.0
    }
  }

  idx <- which(lav$op == "~~" & is.na(lav$ustart) & lav$lhs == lav$rhs)
  if (length(idx) > 0L) lav$ustart[idx] <- 1.0

  idx <- which(lav$op == "~" & is.na(lav$ustart))
  if (length(idx) > 0L) {
    lav_msg_warn(gettext(
      "some regression coefficients are unspecified and will be set to zero"))
  }

  idx <- which(is.na(lav$ustart))
  if (length(idx) > 0L) lav$ustart[idx] <- 0.0

  if (lav_debug()) {
    cat("lav + default values\n")
    print(as.data.frame(lav))
  }

  # set residual variances to enforce a standardized solution
  # but only if no *residual* variances have been specified in the syntax

  if (standardized) {
    # check if factor loadings are smaller than 1.0
    lambda_idx <- which(lav$op == "=~")
    if (any(lav$ustart[lambda_idx] >= 1.0)) {
      lav_msg_warn(gettext("standardized=TRUE but factor loadings are >= 1.0"))
    }

    # check if regression coefficients are smaller than 1.0
    reg_idx <- which(lav$op == "~")
    if (any(lav$ustart[reg_idx] >= 1.0)) {
      lav_msg_warn(gettext(
        "standardized=TRUE but regression coefficients are >= 1.0"))
    }

    # for ordered observed variables, we will get '0.0', but that is ok
    # so there is no need to make a distinction between numeric/ordered
    # here??
    ngroups <- lav_pt_ngroups(lav)
    ov_names <- lav_pt_vnames(lav, "ov")
    ov_nox <- lav_pt_vnames(lav, "ov.nox")
    # lv_names <- lav_pt_vnames(lav, "lv")
    # lv_y <- lav_pt_vnames(lav, "lv.y")
    lv_nox <- lav_pt_vnames(lav, "lv.nox")
    ov_var_idx <- which(lav$op == "~~" & lav$lhs %in% ov_nox &
      lav$rhs == lav$lhs)
    lv_var_idx <- which(lav$op == "~~" & lav$lhs %in% lv_nox &
      lav$rhs == lav$lhs)
    if (any(lav$user[c(ov_var_idx, lv_var_idx)] > 0L)) {
      lav_msg_warn(gettext(
        "if residual variances are specified, please use standardized=FALSE"))
    }

    # new in 0.6-20: - use lav_lisrel_residual_variances
    #                - use lav_lisrel_comp_set_intresvar
    dotdotdot <- list(...)
    dotdotdot$sample.nobs <- sample.nobs
    dotdotdot$fixed.x <- FALSE # for now
    dotdotdot$representation <- "LISREL"
    dotdotdot$composites <- composites
    dotdotdot$correlation <- TRUE # this is the trick
    tmp_fit <- do.call("lavaan", args = c(list(model = lav), dotdotdot))
    # set/get parameters to invoke lav_lisrel_residual_variances
    tmp_lav <- tmp_fit@ParTable
    tmp_x <- lav_model_get_parameters(tmp_fit@Model)
    tmp_model <- lav_model_set_parameters(tmp_fit@Model, x = tmp_x)
    tmp_lav$ustart <- lav_model_get_parameters(tmp_model, type = "user")

    # copy residual values to lav (without assuming parameter tables look the
    # same)
    res_idx <- c(ov_var_idx, lv_var_idx)
    for (i in seq_along(res_idx)) {
      # lookup this parameter in tmp.lav
      idx_in_lav <- res_idx[i]
      this_lhs <- lav$lhs[idx_in_lav]
      idx_in_tmp <- which(tmp_lav$op == "~~" & tmp_lav$lhs == this_lhs &
                          tmp_lav$rhs == tmp_lav$lhs)
      if (length(idx_in_tmp) == 0L) {
        # hm, not found? Give a warning?
      } else {
        vals <- tmp_lav$ustart[idx_in_tmp]
        # check if we have unfortunate values
        bad_idx <- which(!is.finite(vals) | vals < 0)
        vals[bad_idx] <- 1.0 # not pretty, but safe
        lav$ustart[idx_in_lav] <- vals
      }
    }

    # this is what we did <0.6-20
    # fit <- lavaan(model = lav, sample.nobs = sample.nobs, ...)
    # Sigma.hat <- lav_model_sigma(lavmodel = fit@Model)
    # ETA <- lav_model_veta(lavmodel = fit@Model)

    # if (lav_debug()) {
    #   cat("Sigma.hat:\n")
    #   print(Sigma.hat)
    #   cat("Eta:\n")
    #   print(ETA)
    # }

    # # stage 1: standardize LV
    # if (length(lv.nox) > 0L) {
    #   for (g in 1:ngroups) {
    #     var.group <- which(lav$op == "~~" & lav$lhs %in% lv.nox &
    #       lav$rhs == lav$lhs &
    #       lav$group == group.values[g])
    #     eta.idx <- match(lv.nox, lv.names)
    #     lav$ustart[var.group] <- 1 - diag(ETA[[g]])[eta.idx]
    #   }
    # }
    # # refit
    # fit <- lavaan(model = lav, sample.nobs = sample.nobs, ...)
    # Sigma.hat <- lav_model_sigma(lavmodel = fit@Model)

    # if (lav_debug()) {
    #   cat("after stage 1:\n")
    #   cat("Sigma.hat:\n")
    #   print(Sigma.hat)
    # }

    # # stage 2: standardize OV
    # for (g in 1:ngroups) {
    #   var.group <- which(lav$op == "~~" & lav$lhs %in% ov.nox &
    #     lav$rhs == lav$lhs &
    #     lav$group == group.values[g])
    #   ov.idx <- match(ov.nox, ov.names)
    #   lav$ustart[var.group] <- 1 - diag(Sigma.hat[[g]])[ov.idx]
    # }

    # if (lav_debug()) {
    #   cat("after standardisation lav\n")
    #   print(as.data.frame(lav))
    # }
  }


  # unstandardize
  if (!is.null(ov.var)) {
    # FIXME: if ov.var is named, check the order of the elements

    # 1. unstandardize observed variables
    lav$ustart <- lav_unstandardize_ov(partable = lav, ov_var = ov.var)

    # 2. unstandardized latent variables

    if (lav_debug()) {
      cat("after unstandardisation lav\n")
      print(as.data.frame(lav))
    }
  }

  # fit the model without data
  fit <- lavaan(model = lav, sample.nobs = sample.nobs, ...)

  # the model-implied moments for the population
  sigma_hat <- lav_model_sigma(lavmodel = fit@Model)
  mu_hat <- lav_model_mu(lavmodel = fit@Model)
  if (fit@Model@categorical) {
    th <- lav_model_th(lavmodel = fit@Model)
  }

  if (lav_debug()) {
    cat("\nModel-implied moments (before Vale-Maurelli):\n")
    print(sigma_hat)
    print(mu_hat)
    if (exists("TH")) print(th)
  }

  # ngroups
  ngroups <- length(sample.nobs)

  # prepare
  x <- vector("list", length = ngroups)
  # out <- vector("list", length = ngroups)

  for (g in 1:ngroups) {
    cov_1 <- sigma_hat[[g]]

    # if empirical = TRUE, rescale by N/(N-1), so that estimator=ML
    # returns exact results
    if (empirical) {
      cov_1 <- cov_1 * sample.nobs[g] / (sample.nobs[g] - 1)
    }

    # Using sign-invariant method for cross-machine reproducibility
    if (is.null(skewness) && is.null(kurtosis)) {
      if (mass) {
        x[[g]] <- MASS::mvrnorm(
          n = sample.nobs[g],
          mu = mu_hat[[g]],
          Sigma = cov_1,
          empirical = empirical
        )
      } else {
        x[[g]] <- lav_mvrnorm(
          n = sample.nobs[g],
          mu = mu_hat[[g]],
          sigma_1 = cov_1,
          empirical = empirical
        )
      }
    } else {
      # first generate Z
      z <- lav_data_valemaurelli1983(
        n = sample.nobs[g],
        cor_1 = cov2cor(cov_1),
        skewness = skewness, # FIXME: per group?
        kurtosis = kurtosis, mass = mass
      )
      # rescale
      # Note: 'scale()' will first center, and then scale
      # but we need to first scale, and then center...
      # this was reported by Jordan Brace (9 may 2014)
      # X[[g]] <- scale(Z, center = -Mu.hat[[g]],
      #                   scale  = 1/sqrt(diag(COV)))

      # first, we scale
      tmp <- scale(z,
        center = FALSE,
        scale = 1 / sqrt(diag(cov_1))
      )[, , drop = FALSE]

      # then, we center
      x[[g]] <- sweep(tmp, MARGIN = 2, STATS = mu_hat[[g]], FUN = "+")
    }

    # any categorical variables?
    ov_ord <- lav_pt_vnames(lav, type = "ov.ord", group = group_values[g])
    if (length(ov_ord) > 0L) {
      ov_names <- lav_pt_vnames(lav, type = "ov", group = group_values[g])
      # use thresholds to cut
      for (o in ov_ord) {
        o_idx <- which(o == ov_names)
        th_idx <- which(lav$op == "|" & lav$lhs == o &
          lav$group == group_values[g])
        th_val <- c(-Inf, sort(lav$ustart[th_idx]), +Inf)
        x[[g]][, o_idx] <- as.integer(cut(x[[g]][, o_idx], th_val))
      }
    }

    if (return.type == "data.frame") x[[g]] <- as.data.frame(x[[g]])
  }

  if (return.type == "matrix") {
    if (ngroups == 1L) {
      x[[1L]]
    } else {
      x
    }
  } else if (return.type == "data.frame") {
    data_1 <- x[[1L]]

    # if multiple groups, add group column
    if (ngroups > 1L) {
      for (g in 2:ngroups) {
        data_1 <- rbind(data_1, x[[g]])
      }
      data_1$group <- rep(1:ngroups, times = sample.nobs)
    }
    var_names <- lav_pt_vnames(fit@ParTable, type = "ov", group = 1L)
    if (ngroups > 1L) var_names <- c(var_names, "group")
    names(data_1) <- var_names
    if (return.fit) {
      attr(data_1, "fit") <- fit
    }
    data_1
  } else if (return.type == "cov") {
    if (ngroups == 1L) {
      cov(x[[1L]])
    } else {
      cov_list <- lapply(x, cov)
      cov_list
    }
  }
}
lavSimulateData <- lav_data_simulate_old  # synonym #nolint



lav_data_valemaurelli1983 <- function(n = 100L, cor_1, skewness, kurtosis,
                                      mass = FALSE) {
  fleishman1978_abcd <- function(skewness, kurtosis) {
    system_function <- function(x, skewness, kurtosis) {
      b <- x[1L]
      c_1 <- x[2L]
      d <- x[3L]
      eq1 <- b * b + 6 * b * d + 2 * c_1 * c_1 + 15 * d * d - 1
      eq2 <- 2 * c_1 * (b * b + 24 * b * d + 105 * d * d + 2) - skewness
      eq3 <- 24 * (b * d + c_1 * c_1 * (1 + b * b + 28 * b * d) +
        d * d * (12 + 48 * b * d + 141 * c_1 * c_1 + 225 * d * d)) - kurtosis
      eq <- c(eq1, eq2, eq3)
      sum(eq * eq) ## SS
    }

    out <- nlminb(
      start = c(1, 0, 0), objective = system_function,
      scale = 10,
      control = list(trace = 0),
      skewness = skewness, kurtosis = kurtosis
    )
    if (out$convergence != 0 || out$objective > 1e-5) {
      lav_msg_warn(gettext("lav_data_valemaurelli1983 method did not converge,
                   or it did not find the roots"))
    }
    b <- out$par[1L]
    c_1 <- out$par[2L]
    d <- out$par[3L]
    a <- -c_1
    c(a, b, c_1, d)
  }

  get_icov <- function(b1, c1, d1, b2, c2, d2, r) {
    objective_function <- function(x, b1, c1, d1, b2, c2, d2, r) {
      rho <- x[1L]
      eq <- rho * (b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2) +
        rho * rho * (2 * c1 * c2) + rho * rho * rho * (6 * d1 * d2) - r
      eq * eq
    }

    # gradientFunction <- function(x, bcd1, bcd2, R) {
    #
    # }

    out <- nlminb(
      start = r, objective = objective_function,
      scale = 10, control = list(trace = 0),
      b1 = b1, c1 = c1, d1 = d1, b2 = b2, c2 = c2, d2 = d2, r = r
    )
    if (out$convergence != 0 || out$objective > 1e-5)
      lav_msg_warn(gettext("no convergence"))
    rho <- out$par[1L]
    rho
  }

  # number of variables
  nvar <- ncol(cor_1)
  # check skewness
  if (is.null(skewness)) {
    sk <- rep(0, nvar)
  } else if (length(skewness) == nvar) {
    sk <- skewness
  } else if (length(skewness) == 1L) {
    sk <- rep(skewness, nvar)
  } else {
    lav_msg_stop(gettext("skewness has wrong length"))
  }

  if (is.null(kurtosis)) {
    ku <- rep(0, nvar)
  } else if (length(kurtosis) == nvar) {
    ku <- kurtosis
  } else if (length(kurtosis) == 1L) {
    ku <- rep(kurtosis, nvar)
  } else {
    lav_msg_stop(gettext("kurtosis has wrong length"))
  }

  # create Fleishman table
  ftable_1 <- matrix(0, nvar, 4L)
  for (i in 1:nvar) {
    ftable_1[i, ] <- fleishman1978_abcd(skewness = sk[i], kurtosis = ku[i])
  }

  # compute intermediate correlations between all pairs
  icor <- diag(nvar)
  for (j in 1:(nvar - 1L)) {
    for (i in (j + 1):nvar) {
      if (cor_1[i, j] == 0) next
      icor[i, j] <- icor[j, i] <-
        get_icov(ftable_1[i, 2], ftable_1[i, 3], ftable_1[i, 4],
          ftable_1[j, 2], ftable_1[j, 3], ftable_1[j, 4],
          r = cor_1[i, j]
        )
    }
  }

  if (lav_debug()) {
    cat("\nOriginal correlations (for Vale-Maurelli):\n")
    print(cor_1)
    cat("\nIntermediate correlations (for Vale-Maurelli):\n")
    print(icor)
    cat("\nEigen values ICOR:\n")
    print(eigen(icor)$values)
  }

  # generate Z (using sign-invariant method for cross-machine reproducibility)
  if (mass) {
    x <- z <- MASS::mvrnorm(n = n, mu = rep(0, nvar), Sigma = icor)
  } else {
    x <- z <- lav_mvrnorm(n = n, mu = rep(0, nvar), sigma_1 = icor)
  }

  # transform Z using Fleishman constants
  for (i in 1:nvar) {
    x[, i] <- ftable_1[i, 1L] + ftable_1[i, 2L] * z[, i] +
      ftable_1[i, 3L] * z[, i] * z[, i] +
      ftable_1[i, 4L] * z[, i] * z[, i] * z[, i]
  }

  x
}
