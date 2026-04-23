# EFA: exploratory factor analysis
#
# EFA is implemented as a special version of ESEM
# - it is therefore a wrapper around the lavaan() function to simplify
#   the input
# - a lavaan model is generated with a single 'block' that can be rotated
# - the 'default' output produces output that is more in line with traditional
#   EFA software (in R) like factanal() and fa() from the psych package

# YR 20 Sept 2022 - first version

efa <- function(data = NULL,              # nolint start
                nfactors = 1L,
                sample.cov = NULL,
                sample.nobs = NULL,
                rotation = "geomin",
                rotation.args = list(),
                ov.names = NULL,
                bounds = "pos.var",
                ...,
                output = "efa") {        # nolint end
  # rotation.args deprecation handling
  if (!missing(rotation.args)) {
    lav_deprecated_args("rotation", "rotation.args")
  }
  rotation_args <- rotation.args
  ov_names <- ov.names

  if (is.list(rotation)) {
    rotation_args <- modifyList(list(), rotation)
    rotation <- rotation[[1L]]
  }

  # handle dotdotdot
  dotdotdot <- list(...)

  # twolevel?
  twolevel_flag <- !is.null(dotdotdot$cluster)

  # sampling weights?
  sampling_weights_flag <- !is.null(dotdotdot$sampling.weights)

  # check for unallowed arguments
  if (!is.null(dotdotdot$group)) {
    lav_msg_stop(gettext("efa has no support for multiple groups (for now);
                          consider using the cfa() function in combination
                          with the efa() modifier."))
  }
  #if (!is.null(dotdotdot$sampling.weights)) {
  #  lav_msg_stop(gettext("efa has no support for sampling weights (for now);
  #                        consider using the cfa() function in combination
  #                        with the efa() modifier."))
  #}

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
      if (sampling_weights_flag) {
        names_1 <- c(names_1, dotdotdot$sampling.weights)
      }
      data <- data[, names_1, drop = FALSE]
    } else {
      ov_names <- names(data)
      if (twolevel_flag) {
        ov_names <- ov_names[-which(ov_names == dotdotdot$cluster)]
      }
      if (sampling_weights_flag) {
        ov_names <- ov_names[-which(ov_names == dotdotdot$sampling.weights)]
      }
    }
  } else if (!is.null(sample.cov)) {
    ov_names <- rownames(sample.cov)
    if (is.null(ov_names)) {
      ov_names <- colnames(sample.cov)
    }
  }
  # ov_names?
  if (length(ov_names) == 0L) {
    lav_msg_stop(gettext(
      "could not extract variable names from data or sample.cov"))
  }

  # check nfactors
  if (any(nfactors < 1L)) {
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
    if (any(nfactors > nfac_max)) {
      lav_msg_stop(gettextf("when nvar = %1$s the maximum number of factors
                            is %2$s", nvar, nfac_max))
    }
  }

  # output
  output <- tolower(output)
  if (!output %in% c("lavaan", "efa")) {
    lav_msg_stop(gettext("output= must be either \"lavaan\" or \"efa\""))
  }
  if (output == "lavaan" && length(nfactors) > 1L) {
    lav_msg_stop(gettext("when output = \"lavaan\", nfactors must be a
                         single (integer) number."))
  }

  # fit models
  nfits <- length(nfactors)
  out <- vector("list", length = nfits)
  for (f in seq_len(nfits)) {
    # generate model syntax
    model_syntax <- lav_syntax_efa(
      ov.names = ov_names,
      nfactors = nfactors[f],
      twolevel = twolevel_flag
    )
    # call lavaan (using sem())
    fit <- do.call("sem",
      args = c(
        list(
          model = model_syntax,
          data = data,
          sample.cov = sample.cov,
          sample.nobs = sample.nobs,
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
    names(out) <- paste0("nf", nfactors)
    # add loadings element to the end of the list
    # so we an use the non-generic but useful loadings() function
    # from the stats package
    out$loadings <- lav_efa_get_loadings(out)
    class(out) <- c("efaList", "list")
  }

  out
}
