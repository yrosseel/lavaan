# simulate data starting from a user-specified model
#
# initial version: YR 24 jan 2011
# revision for 0.4-11: YR 21 okt 2011
#
#
#
lav_data_simulate_old <- function( # user-specified model
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
                         return.type = "data.frame",
                         return.fit = FALSE,
                         debug = FALSE,
                         standardized = FALSE) {
  if (!missing(debug)) {
    current.debug <- lav_debug()
    if (lav_debug(debug))
      on.exit(lav_debug(current.debug), TRUE)
  }
  if (!is.null(seed)) set.seed(seed)
  # if(!exists(".Random.seed", envir = .GlobalEnv))
  #    runif(1)               # initialize the RNG if necessary
  # RNGstate <- .Random.seed

  # lav_model_partable
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
    lav <- lav_model_partable(
      model = model,
      meanstructure = meanstructure,
      int.ov.free = int.ov.free,
      int.lv.free = int.lv.free,
      marker.int.zero = marker.int.zero,
      composites = composites,
      conditional.x = conditional.x,
      fixed.x = fixed.x,
      orthogonal = orthogonal,
      std.lv = std.lv,
      auto.fix.first = auto.fix.first,
      auto.fix.single = auto.fix.single,
      auto.var = auto.var,
      auto.cov.lv.x = auto.cov.lv.x,
      auto.cov.y = auto.cov.y,
      ngroups = length(sample.nobs)
    )
  }

  group.values <- lav_partable_group_values(lav)
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
    lambda.idx <- which(lav$op == "=~")
    if (any(lav$ustart[lambda.idx] >= 1.0)) {
      lav_msg_warn(gettext("standardized=TRUE but factor loadings are >= 1.0"))
    }

    # check if regression coefficients are smaller than 1.0
    reg.idx <- which(lav$op == "~")
    if (any(lav$ustart[reg.idx] >= 1.0)) {
      lav_msg_warn(gettext(
        "standardized=TRUE but regression coefficients are >= 1.0"))
    }

    # for ordered observed variables, we will get '0.0', but that is ok
    # so there is no need to make a distinction between numeric/ordered
    # here??
    ngroups <- lav_partable_ngroups(lav)
    ov.names <- lav_partable_vnames(lav, "ov")
    ov.nox <- lav_partable_vnames(lav, "ov.nox")
    lv.names <- lav_partable_vnames(lav, "lv")
    lv.y <- lav_partable_vnames(lav, "lv.y")
    lv.nox <- lav_partable_vnames(lav, "lv.nox")
    ov.var.idx <- which(lav$op == "~~" & lav$lhs %in% ov.nox &
      lav$rhs == lav$lhs)
    lv.var.idx <- which(lav$op == "~~" & lav$lhs %in% lv.nox &
      lav$rhs == lav$lhs)
    if (any(lav$user[c(ov.var.idx, lv.var.idx)] > 0L)) {
      lav_msg_warn(gettext(
        "if residual variances are specified, please use standardized=FALSE"))
    }

    # new in 0.6-20: - use lav_lisrel_residual_variances
    #                - use lav_lisrel_composites_variances
    dotdotdot <- list(...)
    dotdotdot$sample.nobs <- sample.nobs
    dotdotdot$fixed.x <- FALSE # for now
    dotdotdot$representation <- "LISREL"
    dotdotdot$composites <- composites
    dotdotdot$correlation <- TRUE # this is the trick
    tmp.fit <- do.call("lavaan", args = c(list(model = lav), dotdotdot))
    # set/get parameters to invoke lav_lisrel_residual_variances
    tmp.lav <- tmp.fit@ParTable
    tmp.x <- lav_model_get_parameters(tmp.fit@Model)
    tmp.model <- lav_model_set_parameters(tmp.fit@Model, x = tmp.x)
    tmp.lav$ustart <- lav_model_get_parameters(tmp.model, type = "user")

    # copy residual values to lav (without assuming parameter tables look the
    # same)
    res.idx <- c(ov.var.idx, lv.var.idx)
    for (i in seq_len(length(res.idx))) {
      # lookup this parameter in tmp.lav
      idx.in.lav <- res.idx[i]
      this.lhs <- lav$lhs[idx.in.lav]
      idx.in.tmp <- which(tmp.lav$op == "~~" & tmp.lav$lhs == this.lhs &
                          tmp.lav$rhs == tmp.lav$lhs)
      if (length(idx.in.tmp) == 0L) {
        # hm, not found? Give a warning?
      } else {
        vals <- tmp.lav$ustart[idx.in.tmp]
        # check if we have unfortunate values
        bad.idx <- which(!is.finite(vals) | vals < 0)
        vals[bad.idx] <- 1.0 # not pretty, but safe
        lav$ustart[idx.in.lav] <- vals
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
    lav$ustart <- lav_unstandardize_ov(partable = lav, ov.var = ov.var)

    # 2. unstandardized latent variables

    if (lav_debug()) {
      cat("after unstandardisation lav\n")
      print(as.data.frame(lav))
    }
  }

  # fit the model without data
  fit <- lavaan(model = lav, sample.nobs = sample.nobs, ...)

  # the model-implied moments for the population
  Sigma.hat <- lav_model_sigma(lavmodel = fit@Model)
  Mu.hat <- lav_model_mu(lavmodel = fit@Model)
  if (fit@Model@categorical) {
    TH <- lav_model_th(lavmodel = fit@Model)
  }

  if (lav_debug()) {
    cat("\nModel-implied moments (before Vale-Maurelli):\n")
    print(Sigma.hat)
    print(Mu.hat)
    if (exists("TH")) print(TH)
  }

  # ngroups
  ngroups <- length(sample.nobs)

  # prepare
  X <- vector("list", length = ngroups)
  out <- vector("list", length = ngroups)

  for (g in 1:ngroups) {
    COV <- Sigma.hat[[g]]

    # if empirical = TRUE, rescale by N/(N-1), so that estimator=ML
    # returns exact results
    if (empirical) {
      COV <- COV * sample.nobs[g] / (sample.nobs[g] - 1)
    }

    # FIXME: change to rmvnorm once we include the library?
    if (is.null(skewness) && is.null(kurtosis)) {
      X[[g]] <- MASS::mvrnorm(
        n = sample.nobs[g],
        mu = Mu.hat[[g]],
        Sigma = COV,
        empirical = empirical
      )
    } else {
      # first generate Z
      Z <- lav_data_valemaurelli1983(
        n = sample.nobs[g],
        COR = cov2cor(COV),
        skewness = skewness, # FIXME: per group?
        kurtosis = kurtosis
      )
      # rescale
      # Note: 'scale()' will first center, and then scale
      # but we need to first scale, and then center...
      # this was reported by Jordan Brace (9 may 2014)
      # X[[g]] <- scale(Z, center = -Mu.hat[[g]],
      #                   scale  = 1/sqrt(diag(COV)))

      # first, we scale
      TMP <- scale(Z,
        center = FALSE,
        scale = 1 / sqrt(diag(COV))
      )[, , drop = FALSE]

      # then, we center
      X[[g]] <- sweep(TMP, MARGIN = 2, STATS = Mu.hat[[g]], FUN = "+")
    }

    # any categorical variables?
    ov.ord <- lav_partable_vnames(lav, type = "ov.ord", group = group.values[g])
    if (length(ov.ord) > 0L) {
      ov.names <- lav_partable_vnames(lav, type = "ov", group = group.values[g])
      # use thresholds to cut
      for (o in ov.ord) {
        o.idx <- which(o == ov.names)
        th.idx <- which(lav$op == "|" & lav$lhs == o &
          lav$group == group.values[g])
        th.val <- c(-Inf, sort(lav$ustart[th.idx]), +Inf)
        X[[g]][, o.idx] <- as.integer(cut(X[[g]][, o.idx], th.val))
      }
    }

    if (return.type == "data.frame") X[[g]] <- as.data.frame(X[[g]])
  }

  if (return.type == "matrix") {
    if (ngroups == 1L) {
      return(X[[1L]])
    } else {
      return(X)
    }
  } else if (return.type == "data.frame") {
    Data <- X[[1L]]

    # if multiple groups, add group column
    if (ngroups > 1L) {
      for (g in 2:ngroups) {
        Data <- rbind(Data, X[[g]])
      }
      Data$group <- rep(1:ngroups, times = sample.nobs)
    }
    var.names <- lav_partable_vnames(fit@ParTable, type = "ov", group = 1L)
    if (ngroups > 1L) var.names <- c(var.names, "group")
    names(Data) <- var.names
    if (return.fit) {
      attr(Data, "fit") <- fit
    }
    return(Data)
  } else if (return.type == "cov") {
    if (ngroups == 1L) {
      return(cov(X[[1L]]))
    } else {
      cov.list <- lapply(X, cov)
      return(cov.list)
    }
  }
}

lav_skewness <- function(x., N1 = TRUE) {
  x <- x.
  x <- x[!is.na(x)]
  N <- length(x)
  mean.x <- mean(x)
  xc <- x - mean.x
  var.x <- var(x)
  if (!N1) var.x <- var.x * (N - 1) / N
  sd.x <- sqrt(var.x)
  sk <- sum(xc * xc * xc) / (sd.x * sd.x * sd.x)
  skewness <- N * sk / ((N - 1) * (N - 2))
  skewness
}

lav_kurtosis <- function(x., N1 = TRUE) {
  x <- x.
  x <- x[!is.na(x)]
  N <- length(x)
  mean.x <- mean(x)
  xc <- x - mean.x
  var.x <- var(x)
  if (!N1) var.x <- var.x * (N - 1) / N
  k <- sum(xc * xc * xc * xc) / (var.x * var.x)
  kurtosis <- N * (N + 1) * k / ((N - 1) * (N - 2) * (N - 3)) - 3 * (N - 1) * (N - 1) / ((N - 2) * (N - 3))
  kurtosis
}



lav_data_valemaurelli1983 <- function(n = 100L, COR, skewness, kurtosis) {
  fleishman1978_abcd <- function(skewness, kurtosis) {
    system.function <- function(x, skewness, kurtosis) {
      b. <- x[1L]
      c. <- x[2L]
      d. <- x[3L]
      eq1 <- b. * b. + 6 * b. * d. + 2 * c. * c. + 15 * d. * d. - 1
      eq2 <- 2 * c. * (b. * b. + 24 * b. * d. + 105 * d. * d. + 2) - skewness
      eq3 <- 24 * (b. * d. + c. * c. * (1 + b. * b. + 28 * b. * d.) +
        d. * d. * (12 + 48 * b. * d. + 141 * c. * c. + 225 * d. * d.)) - kurtosis
      eq <- c(eq1, eq2, eq3)
      sum(eq * eq) ## SS
    }

    out <- nlminb(
      start = c(1, 0, 0), objective = system.function,
      scale = 10,
      control = list(trace = 0),
      skewness = skewness, kurtosis = kurtosis
    )
    if (out$convergence != 0 || out$objective > 1e-5) {
      lav_msg_warn(gettext("lav_data_valemaurelli1983 method did not convergence,
                   or it did not find the roots"))
    }
    b. <- out$par[1L]
    c. <- out$par[2L]
    d. <- out$par[3L]
    a. <- -c.
    c(a., b., c., d.)
  }

  getICOV <- function(b1, c1, d1, b2, c2, d2, R) {
    objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {
      rho <- x[1L]
      eq <- rho * (b1 * b2 + 3 * b1 * d2 + 3 * d1 * b2 + 9 * d1 * d2) +
        rho * rho * (2 * c1 * c2) + rho * rho * rho * (6 * d1 * d2) - R
      eq * eq
    }

    # gradientFunction <- function(x, bcd1, bcd2, R) {
    #
    # }

    out <- nlminb(
      start = R, objective = objectiveFunction,
      scale = 10, control = list(trace = 0),
      b1 = b1, c1 = c1, d1 = d1, b2 = b2, c2 = c2, d2 = d2, R = R
    )
    if (out$convergence != 0 || out$objective > 1e-5)
      lav_msg_warn(gettext("no convergence"))
    rho <- out$par[1L]
    rho
  }

  # number of variables
  nvar <- ncol(COR)
  # check skewness
  if (is.null(skewness)) {
    SK <- rep(0, nvar)
  } else if (length(skewness) == nvar) {
    SK <- skewness
  } else if (length(skewness) == 1L) {
    SK <- rep(skewness, nvar)
  } else {
    lav_msg_stop(gettext("skewness has wrong length"))
  }

  if (is.null(kurtosis)) {
    KU <- rep(0, nvar)
  } else if (length(kurtosis) == nvar) {
    KU <- kurtosis
  } else if (length(kurtosis) == 1L) {
    KU <- rep(kurtosis, nvar)
  } else {
    lav_msg_stop(gettext("kurtosis has wrong length"))
  }

  # create Fleishman table
  FTable <- matrix(0, nvar, 4L)
  for (i in 1:nvar) {
    FTable[i, ] <- fleishman1978_abcd(skewness = SK[i], kurtosis = KU[i])
  }

  # compute intermediate correlations between all pairs
  ICOR <- diag(nvar)
  for (j in 1:(nvar - 1L)) {
    for (i in (j + 1):nvar) {
      if (COR[i, j] == 0) next
      ICOR[i, j] <- ICOR[j, i] <-
        getICOV(FTable[i, 2], FTable[i, 3], FTable[i, 4],
          FTable[j, 2], FTable[j, 3], FTable[j, 4],
          R = COR[i, j]
        )
    }
  }

  if (lav_debug()) {
    cat("\nOriginal correlations (for Vale-Maurelli):\n")
    print(COR)
    cat("\nIntermediate correlations (for Vale-Maurelli):\n")
    print(ICOR)
    cat("\nEigen values ICOR:\n")
    print(eigen(ICOR)$values)
  }

  # generate Z ## FIXME: replace by rmvnorm once we use that package
  X <- Z <- MASS::mvrnorm(n = n, mu = rep(0, nvar), Sigma = ICOR)

  # transform Z using Fleishman constants
  for (i in 1:nvar) {
    X[, i] <- FTable[i, 1L] + FTable[i, 2L] * Z[, i] + FTable[i, 3L] * Z[, i] * Z[, i] +
      FTable[i, 4L] * Z[, i] * Z[, i] * Z[, i]
  }

  X
}
