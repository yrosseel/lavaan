# This file will (eventually) contain functions that can be used to
# 'predict' the values of outcome variables (y), given the values of
# input variables (x).

# first version YR 2 Nov 2022

# method = "conditional.mean" is based on the following article:

# Mark de Rooij, Julian D. Karch, Marjolein Fokkema, Zsuzsa Bakk, Bunga Citra
# Pratiwi & Henk Kelderman (2022) SEM-Based Out-of-Sample Predictions,
# StructuralEquation Modeling: A Multidisciplinary Journal
# DOI:10.1080/10705511.2022.2061494

# YR 31 Jan 2023: we always 'force' meanstructure = TRUE (for now)

# internal function: prepare the data and the y/x indices that are needed
# by both lavPredictY() and lavResidualsY()
lav_predict_y_prepare <- function(object,
                                  newdata = NULL,
                                  ynames = lav_object_vnames(object, "ov.y"),
                                  xnames = lav_object_vnames(object, "ov.x")) {

  lavmodel <- object@Model
  lavdata <- object@Data
  lavimplied <- object@implied

  # check meanstructure
  if (!lavmodel@meanstructure) {
    lavimplied$mean <- lapply(object@SampleStats@mean, as.matrix)
  }

  new_data <- NULL
  newdata_bygroup <- FALSE
  group_idx <- seq_len(lavdata@ngroups)

  # need full data set
  if (is.null(newdata)) {
    # use internal copy:
    if (lavdata@data.type != "full") {
      lav_msg_stop(gettext(
        "sample statistics were used for fitting and newdata is empty"
      ))
    } else if (is.null(lavdata@X[[1]])) {
      lav_msg_stop(gettext("no local copy of data; FIXME!"))
    } else {
      data_obs <- lavdata@X
      ov_names <- lavdata@ov.names
    }
    # eXo <- lavdata@eXo
  } else if (is.list(newdata) && !is.data.frame(newdata)) {
    # newdata is provided per group, as a list with one element per group.
    # Elements may be NULL to skip a group. This is useful when the groups
    # have a different set of variables, so that ynames= and xnames= can be
    # specified per group (see GitHub issue #369). Each element is processed
    # as a single-group dataset, using the ov.names of the corresponding
    # group in the fitted model.
    newdata_bygroup <- TRUE
    if (length(newdata) != lavdata@ngroups) {
      lav_msg_stop(gettextf(
        "newdata is a list of length %1$d, but the model has %2$d group(s);
        please provide one data.frame per group (use NULL to skip a group).",
        length(newdata), lavdata@ngroups
      ))
    }
    ov <- lavdata@ov
    ordered_names <- ov$name[ov$type == "ordered"]
    data_obs <- ov_names <- vector("list", lavdata@ngroups)
    for (g in seq_len(lavdata@ngroups)) {
      if (is.null(newdata[[g]])) {
        next # skip this group; data_obs[[g]] stays NULL
      }
      new_data_g <- lav_lavdata(
        data = newdata[[g]],
        group = NULL, # treat as a single group
        ov_names = lavdata@ov.names[[g]],
        ov_names_x = lavdata@ov.names.x[[g]],
        ordered = ordered_names,
        lavoptions = list(
          std.ov = lavdata@std.ov,
          group.label = character(0L),
          missing = "ml.x", # always!
          warn = TRUE
        ),
        allow_single_case = TRUE
      )
      data_obs[[g]] <- new_data_g@X[[1L]]
      ov_names[[g]] <- new_data_g@ov.names[[1L]]
    }
    group_idx <- which(!vapply(data_obs, is.null, logical(1)))
    if (length(group_idx) == 0L) {
      lav_msg_stop(gettext(
        "newdata is a list, but all its elements are NULL; please provide the
        data for at least one group."
      ))
    }
  } else {
    # newdata is given (a single data.frame containing all groups)!

    # create lavData object
    ov <- lavdata@ov
    new_data <- lav_lavdata(
      data = newdata,
      group = lavdata@group,
      ov_names = lavdata@ov.names,
      ov_names_x = lavdata@ov.names.x,
      ordered = ov$name[ov$type == "ordered"],
      lavoptions = list(
        std.ov = lavdata@std.ov,
        group.label = lavdata@group.label,
        missing = "ml.x", # always!
        warn = TRUE
      ),
      allow_single_case = TRUE
    )
    # if ordered, check if number of levels is still the same (new in 0.6-7)
    if (lavmodel@categorical) {
      orig_ordered_idx <- which(lavdata@ov$type == "ordered")
      orig_ordered_lev <- lavdata@ov$nlev[orig_ordered_idx]
      match_new_idx <- match(
        lavdata@ov$name[orig_ordered_idx],
        new_data@ov$name
      )
      new_ordered_lev <- new_data@ov$nlev[match_new_idx]
      if (any(orig_ordered_lev - new_ordered_lev != 0)) {
        lav_msg_stop(gettext(
          "mismatch number of categories for some ordered variables in
          newdata compared to original data."
        ))
      }
    }

    data_obs <- new_data@X
    # eXo <- newData@eXo
    ov_names <- new_data@ov.names
  } # newdata

  # check ynames
  if (length(ynames) == 0L) {
    lav_msg_stop(gettext(
      "please specify the y-variables in the ynames= argument"
    ))
  }

  if (anyDuplicated((ynames))) {
    lav_msg_stop(gettext(
      "ynames contains duplicate variable names"
    ))
  }

  if (!is.list(ynames)) {
    ynames <- rep(list(ynames), lavdata@ngroups)
  }

  # check xnames
  if (length(xnames) == 0L) {
    lav_msg_stop(gettext(
      "please specify the x-variables in the xnames= argument"
    ))
  }

  if (anyDuplicated((xnames))) {
    lav_msg_stop(gettext(
      "ynames contains duplicate variable names"
    ))
  }

  if (!is.list(xnames)) {
    xnames <- rep(list(xnames), lavdata@ngroups)
  }

  # create y.idx and x.idx (only for the groups we will predict for)
  y_idx <- x_idx <- vector("list", lavdata@ngroups)
  for (g in group_idx) {
    # ynames in ov.names for this group?
    missing_idx <- which(!ynames[[g]] %in% ov_names[[g]])
    if (length(missing_idx) > 0L) {
      lav_msg_stop(
        gettext(
          "some variable names in ynames do not appear in the dataset:"
        ),
        lav_msg_view(ynames[[g]][missing_idx], "none")
      )
    } else {
      y_idx[[g]] <- match(ynames[[g]], ov_names[[g]])
    }

    # xnames in ov.names for this group?
    missing_idx <- which(!xnames[[g]] %in% ov_names[[g]])
    if (length(missing_idx) > 0L) {
      lav_msg_stop(
        gettext(
          "some variable names in xnames do not appear in the dataset:"
        ),
        lav_msg_view(xnames[[g]][missing_idx], "none")
      )
    } else {
      x_idx[[g]] <- match(xnames[[g]], ov_names[[g]])
    }
  }

  list(
    lavmodel = lavmodel, lavdata = lavdata, lavimplied = lavimplied,
    new_data = new_data, data_obs = data_obs, ov_names = ov_names,
    ynames = ynames, xnames = xnames, y_idx = y_idx, x_idx = x_idx,
    newdata_bygroup = newdata_bygroup, group_idx = group_idx
  )
}


# internal function: label + (optionally) assemble the per-group output
# matrices (predictions or residuals) into the user-facing result. This is
# shared by lavPredictY() and lavResidualsY().
lav_predict_y_assemble <- function(out, ynames, lavdata,
                                   newdata = NULL, new_data = NULL,
                                   label = TRUE, assemble = TRUE) {
  # label?
  if (label) {
    # column names
    for (g in seq_len(lavdata@ngroups)) {
      colnames(out[[g]]) <- ynames[[g]]
    }

    # group.labels
    if (lavdata@ngroups > 1L) {
      names(out) <- lavdata@group.label
    }
  }

  # lavaan.matrix
  out <- lapply(out, "class<-", c("lavaan.matrix", "matrix"))

  if (lavdata@ngroups == 1L) {
    res <- out[[1L]]
  } else {
    res <- out
  }

  # assemble multiple groups into a single data.frame?
  if (lavdata@ngroups > 1L && assemble) {
    if (!is.null(newdata)) {
      lavdata <- new_data
    }
    data_1 <- matrix(as.numeric(NA),
      nrow = sum(unlist(lavdata@norig)),
      ncol = ncol(out[[1L]])
    ) # assume == per g
    colnames(data_1) <- colnames(out[[1L]])
    for (g in seq_len(lavdata@ngroups)) {
      data_1[lavdata@case.idx[[g]], ] <- out[[g]]
    }
    data_1 <- as.data.frame(data_1, stringsAsFactors = FALSE)

    if (!is.null(newdata)) {
      data_1[, lavdata@group] <- newdata[, lavdata@group]
    } else {
      # add group
      data_1[, lavdata@group] <- rep(as.character(NA), nrow(data_1))
      if (lavdata@missing == "listwise") {
        # we will lose the group label of omitted variables!
        data_1[unlist(lavdata@case.idx), lavdata@group] <-
          rep(lavdata@group.label, unlist(lavdata@nobs))
      } else {
        data_1[unlist(lavdata@case.idx), lavdata@group] <-
          rep(lavdata@group.label, unlist(lavdata@norig))
      }
    }

    res <- data_1
  }

  res
}


# internal function: label the per-group output matrices when newdata was
# provided as a per-group list (see lav_predict_y_prepare()). Only the
# groups in group_idx are returned; there is no single-data.frame layout to
# reassemble into, so a (labeled) list is returned -- or a single matrix if
# only one group was requested.
lav_predict_y_assemble_bygroup <- function(out, ynames, lavdata, group_idx,
                                           label = TRUE) {
  active <- out[group_idx]
  if (label) {
    for (i in seq_along(group_idx)) {
      colnames(active[[i]]) <- ynames[[group_idx[i]]]
    }
    if (lavdata@ngroups > 1L) {
      names(active) <- lavdata@group.label[group_idx]
    }
  }
  active <- lapply(active, "class<-", c("lavaan.matrix", "matrix"))
  if (length(active) == 1L) {
    active[[1L]]
  } else {
    active
  }
}


# main function
lavPredictY <- function(object,                                    # nolint
                        newdata = NULL,
                        ynames = lav_object_vnames(object, "ov.y"),
                        xnames = lav_object_vnames(object, "ov.x"),
                        method = "conditional.mean",
                        label = TRUE,
                        assemble = TRUE,
                        force_zero_mean = FALSE,
                        lambda = 0,
                      ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  stopifnot(inherits(object, "lavaan"))
  # check object
  object <- lav_object_check_version(object)

  # prepare data + y/x indices
  prep <- lav_predict_y_prepare(
    object = object, newdata = newdata,
    ynames = ynames, xnames = xnames
  )

  # prediction method
  method <- tolower(method)
  if (method == "conditional.mean") {
    out <- lav_predict_y_conditional_mean(
      lavobject = NULL,
      lavmodel = prep$lavmodel, lavdata = prep$lavdata,
      lavimplied = prep$lavimplied,
      data_obs = prep$data_obs, y_idx = prep$y_idx, x_idx = prep$x_idx,
      force_zero_mean = force_zero_mean,
      lambda = lambda
    )
  } else {
    lav_msg_stop(gettext("method must be \"conditional.mean\" (for now)."))
  }

  # label + assemble
  if (prep$newdata_bygroup) {
    lav_predict_y_assemble_bygroup(
      out = out, ynames = prep$ynames, lavdata = prep$lavdata,
      group_idx = prep$group_idx, label = label
    )
  } else {
    lav_predict_y_assemble(
      out = out, ynames = prep$ynames, lavdata = prep$lavdata,
      newdata = newdata, new_data = prep$new_data,
      label = label, assemble = assemble
    )
  }
}


# compute residuals (observed - predicted) for the y-variables, given the
# x-variables. See lavPredictY() for the prediction part. First version:
# YR 2 July 2026 (feature request in GitHub issue #269).
lavResidualsY <- function(object,                                  # nolint
                          newdata = NULL,
                          ynames = lav_object_vnames(object, "ov.y"),
                          xnames = lav_object_vnames(object, "ov.x"),
                          method = "conditional.mean",
                          label = TRUE,
                          assemble = TRUE,
                          force_zero_mean = FALSE,
                          lambda = 0,
                        ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  stopifnot(inherits(object, "lavaan"))
  # check object
  object <- lav_object_check_version(object)

  # prepare data + y/x indices (identical to lavPredictY)
  prep <- lav_predict_y_prepare(
    object = object, newdata = newdata,
    ynames = ynames, xnames = xnames
  )

  # prediction method
  method <- tolower(method)
  if (method == "conditional.mean") {
    ypred <- lav_predict_y_conditional_mean(
      lavobject = NULL,
      lavmodel = prep$lavmodel, lavdata = prep$lavdata,
      lavimplied = prep$lavimplied,
      data_obs = prep$data_obs, y_idx = prep$y_idx, x_idx = prep$x_idx,
      force_zero_mean = force_zero_mean,
      lambda = lambda
    )
  } else {
    lav_msg_stop(gettext("method must be \"conditional.mean\" (for now)."))
  }

  # residuals = observed - predicted (only for the groups we predicted for)
  out <- vector("list", length = prep$lavdata@ngroups)
  for (g in prep$group_idx) {
    yobs_g <- prep$data_obs[[g]][, prep$y_idx[[g]], drop = FALSE]
    out[[g]] <- yobs_g - ypred[[g]]
  }

  # label + assemble
  if (prep$newdata_bygroup) {
    lav_predict_y_assemble_bygroup(
      out = out, ynames = prep$ynames, lavdata = prep$lavdata,
      group_idx = prep$group_idx, label = label
    )
  } else {
    lav_predict_y_assemble(
      out = out, ynames = prep$ynames, lavdata = prep$lavdata,
      newdata = newdata, new_data = prep$new_data,
      label = label, assemble = assemble
    )
  }
}


# method = "conditional.mean"
lav_predict_y_conditional_mean <- function(
    lavobject = NULL, # for convenience
    # object ingredients
    lavmodel = NULL,
    lavdata = NULL,
    lavimplied = NULL,
    # new data
    data_obs = NULL,
    # y and x
    y_idx = NULL,
    x_idx = NULL,
    # options
    force_zero_mean = FALSE,
    lambda = lambda,
    level = 1L) { # not used for now

  # full object?
  if (inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    # lavsamplestats <- lavobject@SampleStats
    lavimplied <- lavobject@implied
  } else {
    stopifnot(
      !is.null(lavmodel), !is.null(lavdata),
      # !is.null(lavsamplestats),
      !is.null(lavimplied)
    )
  }

  # data.obs?
  if (is.null(data_obs)) {
    data_obs <- lavdata@X
  }

  # checks
  if (lavmodel@categorical) {
    lav_msg_stop(gettext("no support for categorical data (yet)."))
  }
  if (lavdata@nlevels > 1L) {
    lav_msg_stop(gettext("no support for multilevel data (yet)."))
  }

  # conditional.x?
  if (lavmodel@conditional.x) {
    sigma_hat <- lav_model_cond2joint_sigma(lavmodel)
    if (lavmodel@meanstructure) {
      mu_hat <- lav_model_cond2joint_mu(lavmodel)
    }
  } else {
    sigma_hat <- lavimplied$cov
    mu_hat <- lavimplied$mean
  }

  # output container
  ypred <- vector("list", length = lavdata@ngroups)

  # run over all groups
  for (g in 1:lavdata@ngroups) {
    # multiple levels?
    if (lavdata@nlevels > 1L) {
      # TODO!
      lav_msg_stop(gettext("no support for multilevel data (yet)!"))
    } else {
      data_obs_g <- data_obs[[g]]

      # no (new)data for this group? (e.g., per-group newdata list with a
      # NULL entry); leave ypred[[g]] as NULL and move on
      if (is.null(data_obs_g)) {
        next
      }

      # model-implied variance-covariance matrix for this group
      cov_g <- sigma_hat[[g]]

      # model-implied mean vector for this group
      if (force_zero_mean) {
        mean_g <- rep(0, ncol(data_obs_g))
      } else {
        mean_g <- as.numeric(mu_hat[[g]])
      }

      # indices (in ov.names)
      y_idx_g <- y_idx[[g]]
      x_idx_g <- x_idx[[g]]

      # partition y/x
      sxx <- cov_g[x_idx_g, x_idx_g, drop = FALSE]
      sxy <- cov_g[x_idx_g, y_idx_g, drop = FALSE]

      # x-data only
      xtest <- data_obs_g[, x_idx_g, drop = FALSE]

      # mx/my
      mx <- mean_g[x_idx_g]
      my <- mean_g[y_idx_g]

      # center using mx
      xtest <- t(t(xtest) - mx)

      # Apply regularization
      sxx <- sxx + lambda * diag(nrow(sxx))

      # prediction rule
      tmp <- xtest %*% solve(sxx, sxy)
      ypred[[g]] <- t(t(tmp) + my)
    } # single level
  } # g

  ypred
}


# Takes a sequence of lambdas and performs k-fold cross-validation to determine
# the best lambda
lavPredictY_cv <- function(                    # nolint
    object,
    data = NULL,
    xnames = lav_object_vnames(object, "ov.x"),
    ynames = lav_object_vnames(object, "ov.y"),
    n_folds = 10L,
    lambda_seq = seq(0, 1, 0.1),
  ...) {

  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  # object should be (or inherit from) a lavaan object
  stopifnot(inherits(object, "lavaan"))
  # check object
  object <- lav_object_check_version(object)

  # results container
  results <- matrix(as.numeric(NA),
    nrow = length(lambda_seq) * n_folds,
    ncol = 2L
  )
  colnames(results) <- c("mse", "lambda")

  # shuffle folds
  folds <- sample(rep(1:n_folds, length.out = nrow(data)))

  # extract Y-data
  y <- as.matrix(data[, ynames, drop = FALSE])

  j <- 0L
  for (i in 1:n_folds) {
    indis <- which(folds == i)
    fold_fit <- try(update(object,
      data = data[-indis, , drop = FALSE],
      warn = FALSE
    ), silent = TRUE)
    if (inherits(fold_fit, "try-error")) {
      lav_msg_warn(gettextf("failed fit in fold %s", i))
      next
    }
    for (l in lambda_seq) {
      j <- j + 1L
      yhat <- lavPredictY(
        fold_fit,
        newdata = data[indis, , drop = FALSE],
        xnames = xnames,
        ynames = ynames,
        lambda = l
      )
      y_error <- y[indis, , drop = FALSE] - yhat
      mse <- mean(y_error * y_error)
      results[j, ] <- c(mse, l)
    }
  }

  # Group by lambda and determine average MSE per group
  avg <- aggregate(results[, "mse"],
    by = list(results[, "lambda"]),
    FUN = mean, na.rm = TRUE
  )
  avg <- avg[order(avg[, 2]), ]
  names(avg) <- c("lambda", "mse")
  lambda_min <- avg[1L, "lambda"]

  list(results = avg, lambda.min = lambda_min)
}
