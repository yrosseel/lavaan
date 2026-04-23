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

# main function
lavPredictY <- function(object,                                    # nolint start
                        newdata = NULL,
                        ynames = lav_object_vnames(object, "ov.y"),
                        xnames = lav_object_vnames(object, "ov.x"),
                        method = "conditional.mean",
                        label = TRUE,
                        assemble = TRUE,
                        force.zero.mean = FALSE,
                        lambda = 0) {                              # nolint end
  stopifnot(inherits(object, "lavaan"))
  # check object
  object <- lav_object_check_version(object)

  lavmodel <- object@Model
  lavdata <- object@Data
  lavimplied <- object@implied

  # check meanstructure
  if (!lavmodel@meanstructure) {
    lavimplied$mean <- lapply(object@SampleStats@mean, as.matrix)
  }

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
  } else {
    # newdata is given!

    # create lavData object
    ov <- lavdata@ov
    new_data <- lav_lavdata(
      data = newdata,
      group = lavdata@group,
      ov.names = lavdata@ov.names,
      ov.names.x = lavdata@ov.names.x,
      ordered = ov$name[ov$type == "ordered"],
      lavoptions = list(
        std.ov = lavdata@std.ov,
        group.label = lavdata@group.label,
        missing = "ml.x", # always!
        warn = TRUE
      ),
      allow.single.case = TRUE
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

  # create y.idx and x.idx
  y_idx <- x_idx <- vector("list", lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
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

  # prediction method
  method <- tolower(method)
  if (method == "conditional.mean") {
    out <- lav_predict_y_conditional_mean(
      lavobject = NULL,
      lavmodel = lavmodel, lavdata = lavdata,
      lavimplied = lavimplied,
      data_obs = data_obs, y_idx = y_idx, x_idx = x_idx,
      force_zero_mean = force.zero.mean,
      lambda = lambda
    )
  } else {
    lav_msg_stop(gettext("method must be \"conditional.mean\" (for now)."))
  }

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
        # we will loose the group label of omitted variables!
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
lavPredictY_cv <- function(                    # nolint start
    object,
    data = NULL,
    xnames = lav_object_vnames(object, "ov.x"),
    ynames = lav_object_vnames(object, "ov.y"),
    n.folds = 10L,
    lambda.seq = seq(0, 1, 0.1)) {             # nolint end

  # object should be (or inherit from) a lavaan object
  stopifnot(inherits(object, "lavaan"))
  # check object
  object <- lav_object_check_version(object)

  # results container
  results <- matrix(as.numeric(NA),
    nrow = length(lambda.seq) * n.folds,
    ncol = 2L
  )
  colnames(results) <- c("mse", "lambda")

  # shuffle folds
  folds <- sample(rep(1:n.folds, length.out = nrow(data)))

  # extract Y-data
  y <- as.matrix(data[, ynames, drop = FALSE])

  j <- 0L
  for (i in 1:n.folds) {
    indis <- which(folds == i)
    fold_fit <- try(update(object,
      data = data[-indis, , drop = FALSE],
      warn = FALSE
    ), silent = TRUE)
    if (inherits(fold_fit, "try-error")) {
      lav_msg_warn(gettextf("failed fit in fold %s", i))
      next
    }
    for (l in lambda.seq) {
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
