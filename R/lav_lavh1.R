# lavH1: compute only the unrestricted ('H1') model: the saturated summary
# statistics -- and, when available, the corresponding loglikelihood --
# that normally end up in the @h1 slot of a fitted lavaan object.
#
# Unlike sem/cfa/lavaan, no model is needed; but if one is provided
# (model syntax, parameter table, or a fitted lavaan object), it is used
# to determine the relevant variables, their role (exogenous, ordered,
# per-level, ...) and the data-related options -- so the result is
# exactly the h1 information a fit of that model would produce.
#
# The workhorse lav_h1_lavobject() is shared with lavCor().
#
# YR 19 July 2026: first released version, based on a June 2024 prototype
#                  and suggestions by Terrence Jorgensen

lavH1 <- function(model = NULL, # nolint
                  data = NULL,
                  ordered = NULL,
                  sampling_weights = NULL,
                  group = NULL,
                  cluster = NULL,
                  missing = "default",
                  ov_names_x = NULL,
                  ov_names_l = list(),
                  estimator = "default",
                  se = "default",
                  test = "none",
                  ...,
                  output = "list",
                  wls_v = NULL,
                  gamma = NULL,
                  add_extra = TRUE,
                  add_labels = TRUE,
                  add_class = TRUE,
                  drop_list_single_group = TRUE) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, TRUE)

  # output
  output <- tolower(output)
  if (output == "h1") {
    output <- "list"
  }
  if (output == "moments") {
    output <- "lavmoments"
  }
  if (!output %in% c("list", "lavaan", "lavmoments")) {
    lav_msg_stop(gettextf(
      "%1$s argument unknown: %2$s",
      "output", lav_msg_view(output)
    ))
  }

  # if the WLS.V and/or gamma (NACOV) matrices are requested, the h1
  # object needs an estimator that defines them: bump a "none" estimator
  # to "default" (which resolves to DWLS/WLSMV for categorical data, and
  # to ML for continuous data)
  need_moment_stats <- (output == "lavmoments") ||
    isTRUE(wls_v) || isTRUE(gamma)
  if (need_moment_stats && tolower(estimator) == "none") {
    estimator <- "default"
  }

  # these options are set by lavH1() itself
  fixed_names <- c("h1", "do.fit", "baseline")
  drop_idx <- which(names(dotdotdot) %in% fixed_names)
  if (length(drop_idx) > 0L) {
    lav_msg_warn(gettextf(
      "the following argument(s) are set by lavH1() and will be ignored: %s",
      lav_msg_view(names(dotdotdot)[drop_idx], "none")
    ))
    dotdotdot <- dotdotdot[-drop_idx]
  }

  # standard errors/test statistics only make sense for the (fitted)
  # unrestricted model itself; for output = "lavaan", the default se
  # ("default") resolves to standard errors for the summary statistics,
  # so that summary() shows their se/z/p-values
  if (output != "lavaan") {
    se <- "none"
    test <- "none"
  }

  # convenience: the first argument may be the data
  if (inherits(model, "data.frame") || inherits(model, "matrix") ||
      inherits(model, "lavData")) {
    if (!is.null(data)) {
      lav_msg_stop(gettext(
        "both the model= and the data= argument seem to contain data."))
    }
    data <- model
    model <- NULL
  }

  # case 1: model is a fitted lavaan object: reuse its (data) slots
  if (inherits(model, "lavaan")) {
    lavobject <- lav_object_check_version(model)
    if (output == "lavaan") {
      lavobject <- lav_h1_lavobject(
        object = lavobject@Data, missing = "default",
        estimator = estimator, se = se, test = test,
        dotdotdot = dotdotdot
      )
    }
    return(lav_h1_output(lavobject,
      output = output, wls_v = wls_v, gamma = gamma,
      add_extra = add_extra,
      add_labels = add_labels, add_class = add_class,
      drop_list_single_group = drop_list_single_group
    ))
  }

  # case 2: model syntax or parameter table: run the standard workflow
  # (without fitting the model) to process the data and compute h1
  if (!is.null(model)) {
    fit0 <- do.call(lavaan, c(list(
      model = model, data = data, ordered = ordered,
      sampling_weights = sampling_weights,
      group = group, cluster = cluster,
      missing = missing, estimator = estimator,
      model.type = "sem", cmd = "sem",
      se = "none", test = "none", baseline = FALSE,
      h1 = TRUE, do.fit = FALSE
    ), dotdotdot))
    if (output == "lavaan") {
      fit0 <- lav_h1_lavobject(
        object = fit0@Data, missing = "default",
        estimator = estimator, se = se, test = test,
        dotdotdot = dotdotdot
      )
    }
    return(lav_h1_output(fit0,
      output = output, wls_v = wls_v, gamma = gamma,
      add_extra = add_extra,
      add_labels = add_labels, add_class = add_class,
      drop_list_single_group = drop_list_single_group
    ))
  }

  # case 3: no model: use all variables in the data
  if (is.null(data)) {
    lav_msg_stop(gettext("both the model= and data= arguments are empty."))
  }
  lavobject <- lav_h1_lavobject(
    object = data,
    ordered = ordered, group = group, cluster = cluster,
    missing = missing, ov_names_x = ov_names_x, ov_names_l = ov_names_l,
    sampling_weights = sampling_weights,
    estimator = estimator, se = se, test = test, baseline = FALSE,
    do_fit = (output == "lavaan"),
    dotdotdot = dotdotdot
  )

  lav_h1_output(lavobject,
    output = output, wls_v = wls_v, gamma = gamma,
    add_extra = add_extra,
    add_labels = add_labels, add_class = add_class,
    drop_list_single_group = drop_list_single_group
  )
}

# workhorse: construct -- and unless do_fit = FALSE, fit -- the
# unrestricted (H1) model as a lavaan object, starting from a raw
# data.frame/matrix or a lavData object; shared by lavH1() and lavCor()
lav_h1_lavobject <- function(object,
                             ordered = NULL,
                             group = NULL,
                             cluster = NULL,
                             missing = "default",
                             ov_names_x = NULL,
                             ov_names_l = list(),
                             sampling_weights = NULL,
                             estimator = "none",
                             se = "none",
                             test = "none",
                             baseline = FALSE,
                             do_fit = TRUE,
                             dotdotdot = list()) {
  # a user-provided do.fit in the dots overrides the do_fit argument
  if (!is.null(dotdotdot$do.fit)) {
    do_fit <- dotdotdot$do.fit
    dotdotdot$do.fit <- NULL
  }

  # extract sampling.weights.normalization from dots (for lav_lavdata() call)
  dots <- dotdotdot
  sampling_weights_normalization <- "group"
  if (!is.null(dots$sampling.weights.normalization)) {
    sampling_weights_normalization <- dots$sampling.weights.normalization
  }

  # check object class
  if (inherits(object, "lavData")) {
    lav_data <- object
    if (lav_data@data.type != "full") {
      lav_msg_stop(gettext(
        "the unrestricted model can only be constructed from full data (not from summary statistics)."))
    }
    # adopt the missing method stored in the lavData object (unless the
    # user asked for a specific one)
    if (missing == "default") {
      missing <- lav_data@missing
    }
  } else if (inherits(object, "data.frame") ||
    inherits(object, "matrix")) {
    object <- as.data.frame(object)
    names_1 <- names(object)
    if (!is.null(group)) {
      names_1 <- names_1[-match(group, names_1)]
    }
    if (!is.null(cluster)) {
      names_1 <- names_1[-match(cluster, names_1)]
    }
    if (!is.null(sampling_weights)) {
      names_1 <- names_1[-match(sampling_weights, names_1)]
    }
    if (is.logical(ordered)) {
      ordered_flag <- ordered
      if (ordered_flag) {
        ordered <- names_1
        if (length(ov_names_x) > 0L) {
          ordered <- ordered[-which(ordered %in% ov_names_x)]
        }
      } else {
        ordered <- character(0L)
      }
    } else if (is.null(ordered)) {
      ordered <- character(0L)
    } else if (!is.character(ordered)) {
      lav_msg_stop(gettext("ordered argument must be a character vector"))
    } else if (length(ordered) == 1L && nchar(ordered) == 0L) {
      ordered <- character(0L)
    } else {
      # check if all names in "ordered" occur in the dataset?
      missing_idx <- which(!ordered %in% names_1)
      if (length(missing_idx) > 0L) {
        lav_msg_warn(gettextf(
          "ordered variable(s): %s could not be found
          in the data and will be ignored",
          lav_msg_view(ordered[missing_idx])))
      }
    }

    # multilevel: which variables per level? (default: all variables at
    # both levels, unless the user provided ov.names.l)
    if (!is.null(cluster)) {
      ngroups_1 <- 1L
      if (!is.null(group)) {
        ngroups_1 <- length(unique(object[[group]]))
      }
      if (is.null(ov_names_l) || length(ov_names_l) == 0L) {
        ov_names_l <- rep(list(rep(list(names_1), 2L)), ngroups_1)
      } else if (!is.list(ov_names_l)) {
        # a single vector of names: use it at both levels
        ov_names_l <- rep(list(rep(list(ov_names_l), 2L)), ngroups_1)
      } else if (length(ov_names_l) == ngroups_1 &&
                 is.list(ov_names_l[[ngroups_1]])) {
        # already a list (per group) of lists (per level)
      } else {
        # a list (per level): recycle over the groups
        ov_names_l <- rep(list(ov_names_l), ngroups_1)
      }
    }

    lav_data <- lav_lavdata(
      data = object, group = group, cluster = cluster,
      ov_names = names_1, ordered = ordered,
      sampling_weights = sampling_weights,
      ov_names_x = ov_names_x,
      ov_names_l = ov_names_l,
      lavoptions = list(
        missing = missing,
        sampling.weights.normalization = sampling_weights_normalization
      )
    )
  } else {
    lav_msg_stop(
      gettext("cannot handle objects of class"),
      paste(class(object), collapse = " ")
    )
  }

  # set default estimator if se != "none"
  categorical <- any(lav_data@ov$type == "ordered")
  if (tolower(se) != "none" && tolower(estimator) == "none") {
    if (categorical) {
      estimator <- "WLSMV"
    } else {
      estimator <- "ML"
    }
  }

  # extract more partable options from dots
  meanstructure <- FALSE
  fixed_x <- FALSE
  mimic <- "lavaan"
  conditional_x <- FALSE
  correlation <- FALSE
  if (!is.null(dots$meanstructure)) {
    meanstructure <- dots$meanstructure
  }
  if (!is.null(dots$correlation)) {
    # correlation structure: the unrestricted partable must fix the
    # variances (and add the ~*~ scaling rows); without this, the
    # information matrix is singular and all standard errors are NA
    correlation <- dots$correlation
  }
  if (lav_data@ngroups > 1L || lav_data@nlevels > 1L || categorical ||
      tolower(missing) %in% c("ml", "fiml", "direct")) {
    meanstructure <- TRUE
  }
  if (!is.null(dots$fixed.x)) {
    fixed_x <- dots$fixed.x
  }
  if (!is.null(dots$mimic)) {
    mimic <- dots$mimic
  }
  if (!is.null(dots$conditional.x)) {
    conditional_x <- dots$conditional.x
  }

  # generate partable for unrestricted model
  pt_un <- lav_pt_unrestricted(
    lavobject = NULL,
    lavdata = lav_data,
    lavoptions = list(
      meanstructure = meanstructure,
      fixed.x = fixed_x,
      conditional.x = conditional_x,
      correlation = correlation,
      group.w.free = FALSE,
      missing = missing,
      estimator = estimator,
      mimic = mimic
    ),
    sample_cov = NULL,
    sample_mean = NULL,
    sample_th = NULL
  )

  do.call(lavaan, c(list(
    slot_par_table = pt_un, slot_data = lav_data,
    model.type = "unrestricted",
    missing = missing,
    baseline = baseline, h1 = TRUE, # must be TRUE!
    do.fit = do_fit,
    se = se, test = test, estimator = estimator
  ), dotdotdot))
}

# extract the h1 information from a (fitted or unfitted) lavaan object
lav_h1_output <- function(lavobject,
                          output = "list",
                          wls_v = NULL,
                          gamma = NULL,
                          add_extra = TRUE,
                          add_labels = TRUE,
                          add_class = TRUE,
                          drop_list_single_group = TRUE) {
  if (output == "lavaan") {
    return(lavobject)
  }

  # resolve the wls_v/gamma defaults: for output = "lavMoments", both are
  # needed to refit categorical data from summary statistics only, so
  # default them to TRUE when the data are categorical (and FALSE
  # otherwise, since a continuous ML refit needs neither); for the plain
  # list output, they are opt-in (default FALSE)
  categorical <- lavobject@Model@categorical
  if (is.null(wls_v)) {
    wls_v <- (output == "lavmoments") && categorical
  }
  if (is.null(gamma)) {
    gamma <- (output == "lavmoments") && categorical
  }

  # output = "lavMoments": a lavMoments-class list that can be passed
  # directly as the data= argument of lavaan()/cfa()/sem()
  if (output == "lavmoments") {
    return(lav_h1_moments(lavobject,
      wls_v = wls_v, gamma = gamma,
      add_labels = add_labels, add_class = add_class
    ))
  }

  # make sure the @h1 slot is filled in (the object may have been fitted
  # with h1 = FALSE); the computation only needs the data/samplestats
  if (length(lavobject@h1) == 0L) {
    lavobject@h1 <- lav_h1_implied_logl(
      lavdata = lavobject@Data,
      lavsamplestats = lavobject@SampleStats,
      lavpartable = lavobject@ParTable,
      lavoptions = lavobject@Options
    )
  }

  # the h1 summary statistics (cov, mean, th, ... per group/block)
  out <- lav_inspect_sampstat(lavobject,
    h1 = TRUE,
    add_labels = add_labels, add_class = add_class,
    drop_list_single_group = drop_list_single_group
  )

  if (!add_extra && !wls_v && !gamma) {
    return(out)
  }

  if (add_extra) {
    # loglikelihood of the unrestricted (h1) model (if available)
    logl <- lavobject@h1$logl
    if (!is.null(logl$loglik) && is.finite(logl$loglik)) {
      out$logl <- logl$loglik
      if (lavobject@Data@ngroups > 1L) {
        logl_group <- logl$loglik.group
        if (add_labels) {
          names(logl_group) <- lavobject@Data@group.label
        }
        if (add_class) {
          class(logl_group) <- c("lavaan.vector", "numeric")
        }
        out$logl.group <- logl_group
      }
    }

    # number of observations
    nobs <- unlist(lavobject@SampleStats@nobs)
    if (lavobject@Data@ngroups > 1L && add_labels) {
      names(nobs) <- lavobject@Data@group.label
    }
    out$nobs <- nobs

    # missing data: coverage proportions
    if (!all(sapply(lavobject@Data@Mp, is.null))) {
      out$coverage <- lav_inspect_mi_coverage(lavobject,
        add_labels = add_labels, add_class = add_class,
        drop_list_single_group = drop_list_single_group
      )
    }

    # multilevel: intraclass correlations
    if (lavobject@Data@nlevels > 1L) {
      out$icc <- lav_inspect_icc(lavobject,
        add_labels = add_labels, add_class = add_class,
        drop_list_single_group = drop_list_single_group
      )
    }
  }

  # optionally: the WLS.V (weight) and gamma (NACOV) matrices, so that the
  # model can be refitted from summary statistics only (see also
  # output = "lavMoments")
  if (wls_v) {
    out$WLS.V <- lav_inspect_wls_v(lavobject,
      add_labels = add_labels, add_class = add_class,
      drop_list_single_group = drop_list_single_group
    )
  }
  if (gamma) {
    out$gamma <- lav_inspect_sampstat_gamma(lavobject,
      add_labels = add_labels, add_class = add_class,
      drop_list_single_group = drop_list_single_group
    )
  }

  out
}

# build a 'lavMoments' object: a list of summary statistics -- structured
# exactly as the sample.cov/sample.mean/sample.th/sample.nobs (and,
# optionally, WLS.V/NACOV) arguments of lavaan() -- that can be passed
# directly as the data= argument to refit a model from summary statistics
# only. Handles multiple groups and conditional.x (where the residual
# moments and the exogenous-covariate statistics travel as attributes of
# sample.cov, following the lavaan() convention).
lav_h1_moments <- function(lavobject,
                           wls_v = TRUE,
                           gamma = TRUE,
                           add_labels = TRUE,
                           add_class = TRUE) {
  lavmodel <- lavobject@Model
  ngroups <- lavobject@Data@ngroups
  conditional_x <- lavmodel@conditional.x
  categorical <- lavmodel@categorical
  meanstructure <- lavmodel@meanstructure

  if (lavobject@Data@nlevels > 1L) {
    lav_msg_stop(gettext(
      "output = \"lavMoments\" is not available for multilevel data."))
  }

  # per-group summary statistics (always keep the per-group list here)
  stat <- lav_inspect_sampstat(lavobject,
    h1 = TRUE, add_labels = add_labels, add_class = add_class,
    drop_list_single_group = FALSE
  )
  th_idx <- NULL
  if (categorical) {
    th_idx <- lav_inspect_th_idx(lavobject,
      add_labels = add_labels, drop_list_single_group = FALSE
    )
  }

  # assemble the sample.cov / sample.mean / sample.th pieces per group.
  # note: a mean/th component may be NULL (no meanstructure / continuous);
  # `list[[g]] <- NULL` would *delete* the slot, so build with lapply and
  # keep the whole component NULL when it is empty for every group
  keys <- if (conditional_x) {
    list(cov = "res.cov", mean = "res.int", th = "res.th")
  } else {
    list(cov = "cov", mean = "mean", th = "th")
  }
  pick <- function(field) {
    got <- lapply(stat, function(sg) sg[[field]])
    if (all(vapply(got, is.null, logical(1)))) NULL else got
  }
  cov_list <- pick(keys$cov)
  mean_list <- pick(keys$mean)
  th_list <- pick(keys$th)
  if (conditional_x) {
    slopes_list <- pick("res.slopes")
    covx_list <- pick("cov.x")
    meanx_list <- pick("mean.x")
  }

  # single group: unwrap the per-group lists; multiple groups: keep them.
  # For conditional.x, the residual slopes and the exogenous-covariate
  # moments travel as attributes of sample.cov (the lavaan() convention);
  # note that for a single group they attach to the covariance *matrix*,
  # but for several groups to the *outer* list, as lists-per-group
  if (ngroups == 1L) {
    sample_cov <- cov_list[[1]]
    sample_mean <- if (!is.null(mean_list)) mean_list[[1]] else NULL
    sample_th <- if (!is.null(th_list)) th_list[[1]] else NULL
    if (conditional_x) {
      attr(sample_cov, "res.slopes") <- slopes_list[[1]]
      attr(sample_cov, "cov.x") <- covx_list[[1]]
      attr(sample_cov, "mean.x") <- meanx_list[[1]]
    }
    if (categorical) {
      attr(sample_th, "th.idx") <- th_idx[[1]]
    }
  } else {
    sample_cov <- cov_list
    sample_mean <- mean_list # NULL if no meanstructure
    sample_th <- th_list # NULL if not categorical
    if (conditional_x) {
      attr(sample_cov, "res.slopes") <- slopes_list
      attr(sample_cov, "cov.x") <- covx_list
      attr(sample_cov, "mean.x") <- meanx_list
    }
    if (categorical) {
      # the th.idx travels as an attribute of the *outer* list (lavaan
      # convention, see the sample.th argument of lavaan())
      attr(sample_th, "th.idx") <- th_idx
    }
  }

  # number of observations (scalar for one group, vector for several)
  sample_nobs <- unlist(lavobject@SampleStats@nobs)
  if (ngroups > 1L && add_labels) {
    names(sample_nobs) <- lavobject@Data@group.label
  }

  # name the per-group pieces (cosmetic, but keeps the group identity)
  if (ngroups > 1L && length(lavobject@Data@group.label) == ngroups) {
    g_lab <- lavobject@Data@group.label
    names(sample_cov) <- g_lab
    if (!is.null(sample_mean)) names(sample_mean) <- g_lab
    if (!is.null(sample_th)) names(sample_th) <- g_lab
  }

  # the options that the refit must reproduce (only those that change the
  # interpretation of the summary statistics; step00 copies these into
  # the call unless the user overrides them). The covariance matrix is the
  # (divided-by-N) maximum-likelihood version that lavaan uses internally,
  # so the refit must NOT rescale it by (N-1)/N
  lav_options <- list(
    conditional.x = conditional_x,
    sample.cov.rescale = FALSE
  )
  if (conditional_x) {
    lav_options$fixed.x <- lavmodel@fixed.x
  }
  if (isTRUE(lavobject@Options$correlation)) {
    lav_options$correlation <- TRUE
  }
  if (ngroups > 1L && length(lavobject@Data@group.label) == ngroups) {
    lav_options$group.label <- lavobject@Data@group.label
  }

  out <- list(sample.cov = sample_cov, sample.nobs = sample_nobs)
  if (!is.null(sample_mean)) {
    out$sample.mean <- sample_mean
  }
  if (!is.null(sample_th)) {
    out$sample.th <- sample_th
  }
  if (wls_v) {
    out$WLS.V <- lav_inspect_wls_v(lavobject,
      add_labels = add_labels, add_class = add_class,
      drop_list_single_group = TRUE
    )
  }
  if (gamma) {
    out$NACOV <- lav_inspect_sampstat_gamma(lavobject,
      add_labels = add_labels, add_class = add_class,
      drop_list_single_group = TRUE
    )
  }
  out$lavOptions <- lav_options

  class(out) <- c("lavMoments", "list")
  out
}

# compact print method for a 'lavMoments' object: describe what it holds
# (the full matrices are usually large) and how to use it
lav_moments_print <- function(x, ...) {
  # fall back to a plain list print if this is not the structure we expect
  if (!("sample.cov" %in% names(x))) {
    y <- x
    class(y) <- "list"
    print(y, ...)
    return(invisible(x))
  }

  cov1 <- if (is.list(x$sample.cov)) x$sample.cov[[1]] else x$sample.cov
  ngroups <- if (is.list(x$sample.cov)) length(x$sample.cov) else 1L
  nvar <- NROW(cov1)
  th_idx <- attr(x$sample.th, "th.idx")
  ncat <- 0L
  if (!is.null(th_idx)) {
    ti <- if (is.list(th_idx)) th_idx[[1]] else th_idx
    ncat <- length(unique(ti[ti > 0L])) # distinct ordered variables
  }
  components <- intersect(
    c("sample.cov", "sample.mean", "sample.th", "WLS.V", "NACOV"),
    names(x)
  )

  cat("lavaan summary statistics (class 'lavMoments')\n\n")
  cat(sprintf("  %-32s %s\n", "Number of groups", ngroups))
  cat(sprintf("  %-32s %s\n", "Number of observations",
    paste(x$sample.nobs, collapse = " ")))
  cat(sprintf("  %-32s %d\n", "Number of variables", nvar))
  if (ncat > 0L) {
    cat(sprintf("  %-32s %s\n", "Ordered (categorical) variables",
      paste0("yes (", ncat, ")")))
  }
  if (isTRUE(x$lavOptions$conditional.x)) {
    cat(sprintf("  %-32s %s\n", "Conditional.x", "yes"))
  }
  cat(sprintf("  %-32s %s\n", "Components", paste(components, collapse = " ")))
  cat("\n")
  cat("Pass this object as the data= argument of lavaan()/cfa()/sem()\n")
  cat("to (re)fit a model from these summary statistics.\n")

  invisible(x)
}
