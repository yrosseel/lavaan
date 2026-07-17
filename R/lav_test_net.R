# NET: nesting and equivalence testing for SEM models
# see Bentler & Satorra (2010, Psychological Methods), and
# Asparouhov & Muthen (2019, SEM journal, 26:2, 302-309)
#
# core idea: model '0' (the more restricted model) is nested within model '1'
# (the less restricted model) if -- and only if -- model 1 can perfectly
# reproduce ANY set of (first and) second order moments that model 0 can
# produce. To check this, we compute the model-implied moments of model 0
# (by default at the parameter estimates; optionally also at randomly
# perturbed parameter values), feed these moments to model 1 as if they were
# sample statistics, and inspect the minimum of the (ML or ULS) fit function:
# if F0 is (numerically) zero, model 0 is nested within model 1; if in
# addition both models have the same degrees of freedom, they are equivalent.
#
# following Asparouhov & Muthen (2019), we use the fit function value F0
# itself (not a chi-square statistic) as the criterion: F0 < crit (1e-7 by
# default) implies nesting, F0 > 10*crit implies non-nesting, and values in
# between are inconclusive.
#
# note: nesting alone does not guarantee that the chi-squared difference
# test behaves as a chi-square variate: if the nesting holds only on the
# boundary of the parameter space (eg a variance fixed to zero), the usual
# regularity conditions fail (Self & Liang, 1987).
#
# YR 17 July 2026: first version

# exported user function: pairwise nesting/equivalence matrix, in the
# spirit of semTools::net()
lavTestNET <- function(object, ..., crit = 1e-7, npoints = 1L,          # nolint
                       iseed = 12345L, model_names = NULL) {
  mcall <- match.call(expand.dots = TRUE)
  dotdotdot <- list(...)
  modp <- if (length(dotdotdot)) {
    sapply(dotdotdot, inherits, "lavaan")
  } else {
    logical(0L)
  }
  if (!inherits(object, "lavaan")) {
    lav_msg_stop(gettext("object is not a lavaan object."))
  }
  if (sum(modp) == 0L) {
    lav_msg_stop(gettext("please provide at least two lavaan objects."))
  }

  # check (and if needed upgrade) all objects
  object <- lav_object_check_version(object)
  dotdotdot[modp] <- lapply(dotdotdot[modp], lav_object_check_version)

  # list of models + names
  mods <- c(list(object), dotdotdot[modp])
  if (!is.null(model_names)) {
    names(mods) <- model_names
  } else {
    names(mods) <- sapply(
      as.list(mcall)[which(c(FALSE, TRUE, modp))],
      function(x) deparse(x)
    )
  }
  names(mods) <- make.unique(names(mods))

  # order by degrees of freedom (least restricted first)
  ndf <- sapply(mods, function(x) x@test[[1]]$df)
  order_idx <- order(ndf)
  mods <- mods[order_idx]
  ndf <- ndf[order_idx]
  nmods <- length(mods)

  # entry [i, j] (i > j): is model i (more df, more restricted) nested
  # within model j (fewer df, less restricted)? for pairs with equal df,
  # both directions are tested: nested both ways means equivalent
  nested <- matrix(as.logical(NA), nmods, nmods)
  fx_mat <- matrix(as.numeric(NA), nmods, nmods)
  diag(nested) <- TRUE
  diag(fx_mat) <- 0
  for (i in seq_len(nmods)[-1L]) {
    for (j in seq_len(i - 1L)) {
      out <- try(
        lav_test_net_pair(
          m1 = mods[[j]], m0 = mods[[i]],
          crit = crit, npoints = npoints, iseed = iseed
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) {
        next
      }
      nested[i, j] <- out$nested
      fx_mat[i, j] <- out$fx
      if (ndf[i] == ndf[j]) {
        # same df: also test the reverse direction (equivalence)
        out2 <- try(
          lav_test_net_pair(
            m1 = mods[[i]], m0 = mods[[j]],
            crit = crit, npoints = npoints, iseed = iseed
          ),
          silent = TRUE
        )
        if (!inherits(out2, "try-error")) {
          nested[j, i] <- out2$nested
          fx_mat[j, i] <- out2$fx
        }
      }
    }
  }

  dimnames(nested) <- dimnames(fx_mat) <-
    list(paste0(names(mods), " (df = ", ndf, ")"), names(mods))

  out <- list(
    nested = nested, fx = fx_mat, df = ndf,
    crit = crit, npoints = npoints
  )
  class(out) <- c("lavaan.net", "list")
  out
}

# print method for class "lavaan.net"
lav_net_print <- function(x, ...) {
  cat("\nNesting and Equivalence Testing (NET)\n\n")
  nested <- x$nested
  nmods <- nrow(nested)
  m_char <- matrix("", nmods, nmods, dimnames = dimnames(nested))
  for (i in seq_len(nmods)) {
    for (j in seq_len(nmods)) {
      if (i == j) {
        m_char[i, j] <- "-"
      } else if (is.na(nested[i, j])) {
        m_char[i, j] <- if (i > j || x$df[i] == x$df[j]) "?" else ""
      } else if (nested[i, j]) {
        m_char[i, j] <- if (isTRUE(nested[j, i])) "equivalent" else "nested"
      } else {
        m_char[i, j] <- "no"
      }
    }
  }
  print(m_char, quote = FALSE, right = TRUE)
  cat("\nEntry [i, j]: is the model in row i nested within the model in",
      "column j?\n")
  cat("(criterion: fit function value <", format(x$crit),
      "computed at", x$npoints, "parameter point(s))\n")
  invisible(x)
}

# workhorse: is model 'm0' (more restricted) nested within model 'm1'
# (less restricted)?
#
# returns a list:
#   nested   -- TRUE / FALSE / NA (NA = could not check, or inconclusive)
#   fx       -- worst (largest) fit function value over all points checked
#   method   -- "partable" (proven parametrically) or "net" (moment refits)
#   reason   -- why the check was skipped (if nested is NA)
lav_test_net_pair <- function(m1 = NULL, m0 = NULL, crit = 1e-7,
                              npoints = 1L, iseed = 12345L,
                              fast = TRUE) {
  out <- list(
    nested = as.logical(NA), fx = as.numeric(NA),
    method = "net", reason = ""
  )

  # can we check this pair at all?
  reason <- lav_test_net_pair_check(m1 = m1, m0 = m0)
  if (nchar(reason) > 0L) {
    out$reason <- reason
    return(out)
  }

  # fast path: if the parameter table of m0 equals the parameter table of
  # m1 with some free parameters fixed and/or constrained to be equal, the
  # nesting holds algebraically, and no refitting is needed
  if (fast && lav_test_net_pt_nested(m0 = m0, m1 = m1)) {
    out$nested <- TRUE
    out$fx <- 0
    out$method <- "partable"
    return(out)
  }

  # perturbation needs random numbers: do not disturb the user's stream
  if (npoints > 1L) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv),
              add = TRUE)
    } else {
      on.exit(suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)),
              add = TRUE)
    }
    set.seed(iseed)
  }

  # the parameter point(s) of m0 where the check is performed
  x_hat <- lav_model_get_parameters(m0@Model)
  nonlin_flag <- length(m0@Model@ceq.nonlinear.idx) > 0L

  fx_points <- numeric(0L)
  conv_points <- logical(0L)
  for (k in seq_len(npoints)) {
    if (k == 1L) {
      x_k <- x_hat
    } else if (nonlin_flag) {
      # we cannot (easily) generate perturbed parameter values that still
      # satisfy nonlinear equality constraints; check the estimates only
      break
    } else {
      x_k <- lav_test_net_perturb(lavmodel = m0@Model, x_hat = x_hat)
      if (is.null(x_k)) {
        next
      }
    }

    # model-implied moments of m0 at this parameter point
    implied <- lav_test_net_implied(lavmodel = m0@Model, x = x_k)
    if (is.null(implied)) {
      next
    }

    # refit the (partable of the) less restricted model m1 to these moments
    refit <- lav_test_net_refit(m1 = m1, m0 = m0, implied = implied)
    if (is.null(refit)) {
      next
    }
    if (refit$fx >= crit && !refit$converged) {
      # a large minimum from a non-converged refit is no evidence of
      # non-nesting; try again with default starting values
      refit2 <- lav_test_net_refit(
        m1 = m1, m0 = m0, implied = implied,
        default_start = TRUE
      )
      if (!is.null(refit2) && refit2$fx < refit$fx) {
        refit <- refit2
      }
    }
    fx_points <- c(fx_points, refit$fx)
    conv_points <- c(conv_points, refit$converged)
  }

  if (length(fx_points) == 0L) {
    out$reason <- gettext("all moment refits failed")
    return(out)
  }
  out$fx <- max(fx_points)

  # verdict, following Asparouhov & Muthen (2019): F0 below crit at every
  # point -> nested; F0 above 10*crit at some (converged) point -> not
  # nested; anything else -> inconclusive
  if (all(fx_points < crit)) {
    out$nested <- TRUE
  } else if (any(fx_points > 10 * crit & conv_points)) {
    out$nested <- FALSE
  } # else: NA (inconclusive)

  out
}

# guards: return "" if the pair can be checked, or a reason why not
lav_test_net_pair_check <- function(m1 = NULL, m0 = NULL) {
  for (m in list(m1, m0)) {
    if (m@Data@nlevels > 1L) {
      return(gettext("multilevel models are not supported"))
    }
    if (m@Model@group.w.free) {
      return(gettext("group.w.free models are not supported"))
    }
    if (.hasSlot(m@Model, "correlation") && m@Model@correlation) {
      return(gettext("correlation structures are not supported"))
    }
    if (m@Model@estimator %in% c("MML", "FML")) {
      return(gettext("marginal ML estimation is not supported"))
    }
    if (m@Data@data.type == "none") {
      return(gettext("no data (or sample statistics)"))
    }
    if (length(unlist(m@Data@ov.names.aux)) > 0L) {
      return(gettext("models with auxiliary variables are not supported"))
    }
  }
  if (m1@Data@ngroups != m0@Data@ngroups) {
    return(gettext("unequal number of groups"))
  }
  if (!identical(m1@Data@group.label, m0@Data@group.label)) {
    return(gettext("the groups (or their order) differ"))
  }
  ov1 <- lapply(m1@Data@ov.names, sort)
  ov0 <- lapply(m0@Data@ov.names, sort)
  if (!identical(ov1, ov0)) {
    return(gettext("different sets of observed variables"))
  }
  if (m1@Model@meanstructure != m0@Model@meanstructure) {
    return(gettext("only one of the models has a meanstructure"))
  }
  if (m1@Model@categorical != m0@Model@categorical) {
    return(gettext("only one of the models is categorical"))
  }
  if (m0@Model@categorical) {
    if (m0@Model@conditional.x || m1@Model@conditional.x) {
      return(gettext(
        "categorical models with conditional.x are not supported"))
    }
    if (!identical(m1@Data@ov.names, m0@Data@ov.names)) {
      return(gettext(
        "categorical models with a different variable order are not
        supported"))
    }
  }
  ""
}

# model-implied (unconditional) moments of 'lavmodel' with free parameters
# set to 'x'; returns NULL if a model-implied covariance matrix is not
# positive definite
lav_test_net_implied <- function(lavmodel = NULL, x = NULL) {
  lavmodel <- lav_model_set_parameters(lavmodel, x = x)
  implied <- lav_model_implied(lavmodel)
  if (lavmodel@conditional.x) {
    implied <- lav_model_implied_cond2uncond(implied)
  }
  # all covariance matrices must be positive definite
  for (g in seq_along(implied$cov)) {
    ev <- eigen(implied$cov[[g]], symmetric = TRUE,
                only.values = TRUE)$values
    if (ev[length(ev)] < max(abs(ev)) * .Machine$double.eps^0.75) {
      return(NULL)
    }
  }
  implied
}

# a randomly perturbed parameter vector for the restricted model,
# respecting linear equality constraints, inequality constraints, and
# positive definiteness of the implied covariance matrices
lav_test_net_perturb <- function(lavmodel = NULL, x_hat = NULL,
                                 magnitude = 0.1, max_tries = 8L) {
  lineq_flag <- lavmodel@eq.constraints &&
    ncol(lavmodel@eq.constraints.K) > 0L
  cin_flag <- !is.null(body(lavmodel@cin.function)) &&
    length(lavmodel@x.cin.idx) > 0L

  for (i in seq_len(max_tries)) {
    if (lineq_flag) {
      # perturb in the reduced (unconstrained) coordinates
      k_mat <- lavmodel@eq.constraints.K
      k0 <- lavmodel@eq.constraints.k0
      u <- qr.solve(k_mat, x_hat - k0)
      u <- u + runif(length(u), -magnitude, magnitude) * pmax(1, abs(u))
      x_new <- as.numeric(k_mat %*% u) + k0
    } else {
      x_new <- x_hat +
        runif(length(x_hat), -magnitude, magnitude) * pmax(1, abs(x_hat))
    }

    # respect inequality constraints
    if (cin_flag) {
      con <- try(lavmodel@cin.function(x_new), silent = TRUE)
      if (inherits(con, "try-error") || any(con < 0)) {
        magnitude <- magnitude / 2
        next
      }
    }

    # the implied covariance matrices must be positive definite
    implied <- try(lav_test_net_implied(lavmodel = lavmodel, x = x_new),
                   silent = TRUE)
    if (inherits(implied, "try-error") || is.null(implied)) {
      magnitude <- magnitude / 2
      next
    }
    return(x_new)
  }
  NULL
}

# fit the parameter table of m1 to the implied moments of m0; returns
# list(fx = minimum fit function value, converged = flag), or NULL if the
# refit failed
lav_test_net_refit <- function(m1 = NULL, m0 = NULL, implied = NULL,
                               default_start = FALSE) {
  pt <- parTable(m1)
  if (default_start) {
    pt$ustart[pt$free > 0L] <- as.numeric(NA)
  } else {
    # the estimated values of m1 are (usually) excellent starting values
    pt$ustart[pt$free > 0L] <- pt$est[pt$free > 0L]
  }
  pt$est <- pt$se <- pt$start <- NULL

  ngroups <- m0@Data@ngroups
  ov_names <- m0@Data@ov.names
  if (m0@Model@conditional.x) {
    # lav_model_implied_cond2uncond() returns the joint moments in
    # [y, x] order
    for (g in seq_len(ngroups)) {
      ov_x <- m0@Data@ov.names.x[[g]]
      ov_names[[g]] <- c(setdiff(ov_names[[g]], ov_x), ov_x)
    }
  }
  categorical <- m0@Model@categorical
  meanstructure <- m1@Model@meanstructure

  # sample statistics (= the model-implied moments of m0)
  cov_list <- implied$cov
  for (g in seq_len(ngroups)) {
    dimnames(cov_list[[g]]) <- list(ov_names[[g]], ov_names[[g]])
  }
  if (length(m0@Data@group.label) > 0L) {
    names(cov_list) <- m0@Data@group.label
  }
  mean_list <- NULL
  if (meanstructure) {
    mean_list <- implied$mean
    for (g in seq_len(ngroups)) {
      mean_list[[g]] <- as.numeric(mean_list[[g]])
      names(mean_list[[g]]) <- ov_names[[g]]
    }
  }
  th_list <- NULL
  if (categorical) {
    th_list <- lapply(implied$th, as.numeric)
    attr(th_list, "th.idx") <- m0@SampleStats@th.idx
  }
  nobs <- unlist(m0@Data@nobs)

  # if m1 was fitted with conditional.x = TRUE, lavaan cannot derive the
  # conditional sample statistics from joint moments: decompose the joint
  # moments here, with respect to the exogenous set of m1, and pass them
  # via the res.slopes/cov.x/mean.x attributes of sample.cov
  if (m1@Model@conditional.x) {
    res_slopes <- cov_x <- mean_x <- vector("list", ngroups)
    for (g in seq_len(ngroups)) {
      ov_all <- ov_names[[g]]
      ov_x <- m1@Data@ov.names.x[[g]]
      y_idx <- which(!ov_all %in% ov_x)
      x_idx <- match(ov_x, ov_all)
      s_full <- cov_list[[g]]
      s_xx <- s_full[x_idx, x_idx, drop = FALSE]
      s_yx <- s_full[y_idx, x_idx, drop = FALSE]
      beta <- s_yx %*% solve(s_xx)
      cov_list[[g]] <- s_full[y_idx, y_idx, drop = FALSE] -
        beta %*% t(s_yx)
      dimnames(cov_list[[g]]) <- list(ov_all[y_idx], ov_all[y_idx])
      dimnames(beta) <- list(ov_all[y_idx], ov_x)
      res_slopes[[g]] <- beta
      cov_x[[g]] <- s_xx
      dimnames(cov_x[[g]]) <- list(ov_x, ov_x)
      mu <- if (is.null(mean_list)) {
        numeric(length(ov_all))
      } else {
        mean_list[[g]]
      }
      mean_x[[g]] <- mu[x_idx]
      names(mean_x[[g]]) <- ov_x
      if (!is.null(mean_list)) {
        mean_list[[g]] <- as.numeric(mu[y_idx] - beta %*% mu[x_idx])
        names(mean_list[[g]]) <- ov_all[y_idx]
      }
    }
    attr(cov_list, "res.slopes") <- res_slopes
    attr(cov_list, "cov.x") <- cov_x
    attr(cov_list, "mean.x") <- mean_x
  }

  # single group: plain matrix/vector input
  if (ngroups == 1L) {
    extra <- attributes(cov_list)[c("res.slopes", "cov.x", "mean.x")]
    cov_list <- cov_list[[1]]
    for (a in names(extra)) {
      if (!is.null(extra[[a]])) {
        attr(cov_list, a) <- extra[[a]][[1]]
      }
    }
    if (!is.null(mean_list)) {
      mean_list <- mean_list[[1]]
    }
    if (!is.null(th_list)) {
      th_list <- th_list[[1]]
      attr(th_list, "th.idx") <- m0@SampleStats@th.idx[[1]]
    }
  }

  # the estimator is irrelevant for the nesting question (the minimum is
  # zero if and only if the moments are reproduced exactly); we use ML for
  # continuous data (matching Asparouhov & Muthen, 2019) and ULS whenever
  # ML is not available
  estimator <- if (categorical) "ULS" else "ML"

  args <- list(
    model = pt,
    sample_cov = cov_list,
    sample_nobs = nobs,
    "sample.cov.rescale" = FALSE,
    estimator = estimator,
    likelihood = "normal",
    meanstructure = meanstructure,
    "fixed.x" = m1@Options$fixed.x,
    # the moments we pass are always the *joint* moments; if m1 was fitted
    # with conditional.x = TRUE, we keep that setting (lavaan derives the
    # conditional sample statistics from the joint moments itself)
    "conditional.x" = m1@Model@conditional.x,
    missing = "listwise",
    se = "none",
    test = "none",
    baseline = FALSE,
    h1 = FALSE,
    bounds = "none"
  )
  if (!is.null(mean_list)) {
    args$sample_mean <- mean_list
  }
  if (categorical) {
    args$sample_th <- th_list
    args$parameterization <- m1@Options$parameterization
  }
  if (default_start) {
    args$start <- "simple"
  }

  fit_new <- try(suppressWarnings(do.call(lavaan, args)), silent = TRUE)
  if (inherits(fit_new, "try-error")) {
    return(NULL)
  }
  fx <- as.numeric(fit_new@optim$fx)
  if (!is.finite(fx)) {
    return(NULL)
  }
  list(fx = fx, converged = isTRUE(fit_new@optim$converged))
}

# fast path: try to PROVE that m0 is parametrically nested within m1, by
# comparing the parameter tables: m0 must equal m1 with some free
# parameters fixed and/or constrained to be equal. Parametric nesting
# implies covariance nesting, so a TRUE verdict is definitive; FALSE only
# means that the (more expensive) moment refits are needed.
lav_test_net_pt_nested <- function(m0 = NULL, m1 = NULL) {
  p0 <- m0@ParTable
  p1 <- m1@ParTable

  # bail out on anything beyond plain (in)equality constraints via labels
  # or plabels
  if (any(p1$op %in% c("<", ">"))) {
    return(FALSE) # inequality constraints restrict the space of m1
  }
  ok_ops <- c("=~", "~", "~~", "~1", "|", "~*~", "==", ":=", "<", ">")
  if (!all(p0$op %in% ok_ops) || !all(p1$op %in% ok_ops)) {
    return(FALSE)
  }
  # bounds also restrict the parameter space of m1
  if (!is.null(p1$lower) && any(is.finite(p1$lower[p1$free > 0L]))) {
    return(FALSE)
  }
  if (!is.null(p1$upper) && any(is.finite(p1$upper[p1$free > 0L]))) {
    return(FALSE)
  }

  # canonical description of one parameter table:
  #   key    -> (free, value, class id)
  # returns NULL if the table is too complex for this fast path
  describe <- function(p) {
    par_idx <- which(p$op %in% c("=~", "~", "~~", "~1", "|", "~*~"))
    lhs <- p$lhs[par_idx]
    rhs <- p$rhs[par_idx]
    op <- p$op[par_idx]
    # canonicalize symmetric (co)variances
    swap_idx <- which(op == "~~" & lhs > rhs)
    if (length(swap_idx) > 0L) {
      tmp <- lhs[swap_idx]
      lhs[swap_idx] <- rhs[swap_idx]
      rhs[swap_idx] <- tmp
    }
    key <- paste0(lhs, op, rhs, ".b", p$block[par_idx])
    if (any(duplicated(key))) {
      return(NULL)
    }
    free <- p$free[par_idx] > 0L
    value <- ifelse(free, as.numeric(NA), p$est[par_idx])
    # a parameter that is absent from a table is implicitly fixed to zero,
    # but only for these operators (a missing threshold or scaling factor
    # is never zero)
    zero_ok <- op %in% c("=~", "~", "~~", "~1")

    # equality classes: shared labels, and simple "==" rows in which both
    # sides are (p)labels of parameters
    class_id <- rep(as.integer(NA), length(key))
    all_labels <- character(length(key))
    if (!is.null(p$label)) {
      all_labels <- p$label[par_idx]
    }
    plabels <- character(length(key))
    if (!is.null(p$plabel)) {
      plabels <- p$plabel[par_idx]
    }
    # start: rows sharing a non-empty label form one class
    class_list <- list()
    lab_idx <- which(nchar(all_labels) > 0L)
    for (lab in unique(all_labels[lab_idx])) {
      class_list <- c(class_list, list(which(all_labels == lab)))
    }
    # add "==" rows
    eq_idx <- which(p$op == "==")
    for (i in eq_idx) {
      lhs_match <- which(plabels == p$lhs[i] | all_labels == p$lhs[i])
      rhs_match <- which(plabels == p$rhs[i] | all_labels == p$rhs[i])
      if (length(lhs_match) == 0L || length(rhs_match) == 0L) {
        # a constraint involving := definitions or expressions
        return(NULL)
      }
      class_list <- c(class_list, list(c(lhs_match, rhs_match)))
    }
    # merge overlapping classes
    if (length(class_list) > 0L) {
      changed <- TRUE
      while (changed) {
        changed <- FALSE
        i <- 1L
        while (i < length(class_list)) {
          j <- i + 1L
          while (j <= length(class_list)) {
            if (length(intersect(class_list[[i]], class_list[[j]])) > 0L) {
              class_list[[i]] <- union(class_list[[i]], class_list[[j]])
              class_list[[j]] <- NULL
              changed <- TRUE
            } else {
              j <- j + 1L
            }
          }
          i <- i + 1L
        }
      }
      for (ci in seq_along(class_list)) {
        class_id[class_list[[ci]]] <- ci
      }
    }
    list(key = key, free = free, value = value, class_id = class_id,
         zero_ok = zero_ok)
  }

  d0 <- describe(p0)
  d1 <- describe(p1)
  if (is.null(d0) || is.null(d1)) {
    return(FALSE)
  }

  # every parameter of m1 must be present in m0 in a (weakly) more
  # restricted form
  tol <- 1e-10
  idx0 <- match(d1$key, d0$key)
  for (i in seq_along(d1$key)) {
    j <- idx0[i]
    if (is.na(j)) {
      # absent in m0: only allowed if implicitly fixed to zero there, and
      # m1 must then allow the zero value (free, or fixed to zero)
      if (!d1$zero_ok[i]) {
        return(FALSE)
      }
      if (!d1$free[i] && abs(d1$value[i]) > tol) {
        return(FALSE)
      }
      next
    }
    if (!d1$free[i]) {
      # fixed in m1: must be fixed to the same value in m0
      if (d0$free[j]) {
        return(FALSE)
      }
      if (is.na(d1$value[i]) || is.na(d0$value[j]) ||
          abs(d1$value[i] - d0$value[j]) > tol) {
        return(FALSE)
      }
    }
    # free in m1: m0 may do anything with this parameter
  }
  # parameters of m0 that are absent from m1 must be fixed to zero in m0
  # (a parameter absent from BOTH tables is fixed to zero in both)
  extra0 <- which(!d0$key %in% d1$key)
  for (j in extra0) {
    if (!d0$zero_ok[j] || d0$free[j] || is.na(d0$value[j]) ||
        abs(d0$value[j]) > tol) {
      return(FALSE)
    }
  }

  # every equality class of m1 must be preserved in m0: its members must
  # all be in a single class in m0, or all fixed to the same value
  cls1 <- unique(d1$class_id[!is.na(d1$class_id)])
  for (ci in cls1) {
    members <- which(d1$class_id == ci)
    j_members <- idx0[members]
    if (any(is.na(j_members))) {
      return(FALSE)
    }
    if (all(!d0$free[j_members])) {
      vals <- d0$value[j_members]
      if (any(is.na(vals)) || max(vals) - min(vals) > tol) {
        return(FALSE)
      }
    } else if (all(d0$free[j_members])) {
      cid <- unique(d0$class_id[j_members])
      if (length(cid) != 1L || is.na(cid[1])) {
        return(FALSE)
      }
    } else {
      return(FALSE)
    }
  }

  TRUE
}
