# EFA: exploratory factor analysis
#
# EFA is implemented as a special version of ESEM
# - it is therefore a wrapper around the lavaan() function to simplify
#   the input
# - a lavaan model is generated with a single 'block' that can be rotated
# - the 'default' output produces output that is more in line with traditional
#   EFA software (in R) like factanal() and fa() from the psych package

# YR 20 Sept 2022 - first version

efa <- function(data = NULL,
                nfactors = 1L,
                sample_cov = NULL,
                sample_nobs = NULL,
                rotation = "geomin",
                rotation_args = list(),
                ov_names = NULL,
                bounds = "pos.var",
                ...,
                output = "efa") {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, TRUE)
  # rotation.args deprecation handling
  if (!missing(rotation_args)) {
    lav_deprecated_args("rotation", "rotation_args")
  }

  if (is.list(rotation)) {
    rotation_args <- modifyList(list(), rotation)
    rotation_args <- lav_snake_case(rotation_args)
    rotation <- rotation[[1L]]
  }

  # twolevel?
  twolevel_flag <- !is.null(dotdotdot$cluster)

  # multiple groups?
  group_flag <- !is.null(dotdotdot$group)

  # sampling weights?
  sampling_weights_flag <- !is.null(dotdotdot$sampling.weights)

  # if data= argument is used, convert to data.frame (eg matrix, tibble, ...)
  if (!is.null(data) && !inherits(data, "lavMoments")) {
    data <- as.data.frame(data)
  }

  # handle ov_names
  if (!is.null(data) && inherits(data, "lavMoments")) {
    if ("sample.cov" %in% names(data)) {
      ov_names <- rownames(data$sample.cov)
      if (is.null(ov_names)) {
        ov_names <- colnames(data$sample.cov)
      }
    } else {
      lav_msg_stop(gettext(
        "When data= is of class lavMoments, it must contain sample.cov"))
    }

  } else if (!is.null(data) && inherits(data, "data.frame")) {
    if (length(ov_names) > 0L) {
      names_1 <- ov_names
      if (twolevel_flag) {
        names_1 <- c(names_1, dotdotdot$cluster)
      }
      if (group_flag) {
        names_1 <- c(names_1, dotdotdot$group)
      }
      if (sampling_weights_flag) {
        names_1 <- c(names_1, dotdotdot$sampling.weights)
      }
      data <- data[, names_1, drop = FALSE]
    } else {
      ov_names <- names(data)
      if (twolevel_flag) {
        if (!all(dotdotdot$cluster %in% ov_names)) {
          lav_msg_stop(gettextf(
            "cluster variable(s) %s not found in data.",
            lav_msg_view(dotdotdot$cluster[!dotdotdot$cluster %in% ov_names])))
        }
        ov_names <- setdiff(ov_names, dotdotdot$cluster)
      }
      if (group_flag) {
        if (!all(dotdotdot$group %in% ov_names)) {
          lav_msg_stop(gettextf(
            "group variable(s) %s not found in data.",
            lav_msg_view(dotdotdot$group[!dotdotdot$group %in% ov_names])))
        }
        ov_names <- setdiff(ov_names, dotdotdot$group)
      }
      if (sampling_weights_flag) {
        if (!all(dotdotdot$sampling.weights %in% ov_names)) {
          lav_msg_stop(gettextf(
            "sampling.weights variable(s) %s not found in data.",
            lav_msg_view(dotdotdot$sampling.weights[
              !dotdotdot$sampling.weights %in% ov_names])))
        }
        ov_names <- setdiff(ov_names, dotdotdot$sampling.weights)
      }
    }
  } else if (!is.null(sample_cov)) {
    ov_names <- rownames(sample_cov)
    if (is.null(ov_names)) {
      ov_names <- colnames(sample_cov)
    }
  }
  # ov_names?
  if (length(ov_names) == 0L) {
    lav_msg_stop(gettext(
      "could not extract variable names from data or sample_cov"))
  }

  # determine the block structure (blocks are ordered group-major, level-minor,
  # ie block = (g - 1L) * nlevels + level)
  nlevels <- if (twolevel_flag) 2L else 1L
  if (group_flag && !is.null(data) && inherits(data, "data.frame")) {
    ngroups <- length(unique(data[[dotdotdot$group]]))
  } else if (!is.null(sample_cov) && is.list(sample_cov) &&
    !inherits(sample_cov, "lavMoments")) {
    ngroups <- length(sample_cov)
  } else {
    ngroups <- 1L
  }
  nblocks <- ngroups * nlevels

  # normalize 'nfactors' into a list of per-block (integer) vectors, one element
  # per candidate model; a scalar or plain vector keeps its original meaning,
  # ie 'candidate models, the same factor count in every block'
  if (is.list(nfactors)) {
    nfactors_list <- lapply(nfactors, function(x) {
      x <- as.integer(x)
      if (length(x) == 1L) {
        x <- rep.int(x, nblocks)
      }
      if (length(x) != nblocks) {
        lav_msg_stop(gettextf("each element of nfactors must have length 1 or
                              the number of blocks (= %1$s); found length %2$s.",
          nblocks, length(x)))
      }
      x
    })
  } else {
    nfactors_list <- lapply(as.integer(nfactors), function(k) {
      rep.int(k, nblocks)
    })
  }
  # labels for the resulting fits ('nf2', or 'nf1_2' if it varies per block)
  nfactors_labels <- vapply(nfactors_list, function(x) {
    if (length(unique(x)) == 1L) {
      paste0("nf", x[1L])
    } else {
      paste0("nf", paste(x, collapse = "_"))
    }
  }, character(1L))

  # check nfactors
  all_nfactors <- unlist(nfactors_list)
  if (any(all_nfactors < 1L)) {
    lav_msg_stop(gettext("nfactors must be greater than zero."))
  } else {
    # check for maximum number of factors
    # Fixme: can we do this more efficiently? also holds for categorical?
    nvar <- length(ov_names)
    p_star <- nvar * (nvar + 1) / 2
    nfac_max <- 0L
    for (nfac in seq_len(nvar)) {
      # compute number of free parameters
      npar <- nfac * nvar + nfac * (nfac + 1L) / 2 + nvar - nfac^2
      if (npar > p_star) {
        nfac_max <- nfac - 1L
        break
      }
    }
    if (any(all_nfactors > nfac_max)) {
      lav_msg_stop(gettextf("when nvar = %1$s the maximum number of factors
                            is %2$s", nvar, nfac_max))
    }
  }

  # output
  output <- tolower(output)
  if (!output %in% c("lavaan", "efa")) {
    lav_msg_stop(gettext("output= must be either \"lavaan\" or \"efa\""))
  }
  if (output == "lavaan" && length(nfactors_list) > 1L) {
    lav_msg_stop(gettext("when output = \"lavaan\", nfactors must be a
                         single (integer) number."))
  }

  # fit models
  nfits <- length(nfactors_list)
  out <- vector("list", length = nfits)
  for (f in seq_len(nfits)) {
    # generate model syntax
    model_syntax <- lav_syntax_efa(
      ov_names = ov_names,
      nfactors = nfactors_list[[f]],
      nlevels = nlevels,
      ngroups = ngroups
    )
    # call lavaan (using sem())
    fit <- do.call("sem",
      args = c(
        list(
          model = model_syntax,
          data = data,
          sample_cov = sample_cov,
          sample_nobs = sample_nobs,
          rotation = rotation,
          rotation.args = rotation_args,
          bounds = bounds,
          cmd = "efa"
        ),
        dotdotdot
      )
    )

    if (output == "efa") {
      fit@Options$model.type <- "efa"
    }

    out[[f]] <- fit
  }

  # class
  if (nfits == 1L && output == "lavaan") {
    out <- out[[1]]
  } else {
    names(out) <- nfactors_labels
    # add loadings element to the end of the list
    # so we an use the non-generic but useful loadings() function
    # from the stats package
    out$loadings <- lav_efa_get_loadings(out)
    class(out) <- c("efaList", "list")
  }

  out
}
