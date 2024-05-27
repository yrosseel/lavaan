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
                sample.cov = NULL,
                sample.nobs = NULL,
                rotation = "geomin",
                rotation.args = list(),
                ov.names = names(data),
                bounds = "pos.var",
                ...,
                output = "efa") {
  # handle dotdotdot
  dotdotdot <- list(...)

  # twolevel?
  twolevel.flag <- !is.null(dotdotdot$cluster)

  # check for unallowed arguments
  if (!is.null(dotdotdot$group)) {
    lav_msg_stop(gettext("efa has no support for multiple groups (for now)"))
  }

  # handle ov.names
  if (!is.null(data) && inherits(data, "lavMoments")) {
    if ("sample.cov" %in% names(data)) {
      ov.names <- rownames(data$sample.cov)
      if (is.null(ov.names)) {
        ov.names <- colnames(data$sample.cov)
      }
    } else {
      lav_msg_stop(gettext(
        "When data= is of class lavMoments, it must contain sample.cov"))
    }

  } else if (!is.null(data) && inherits(data, "data.frame")) {
    if (length(ov.names) > 0L) {
      if (twolevel.flag) {
        data <- data[, c(ov.names, dotdotdot$cluster)]
      } else {
        data <- data[, ov.names, drop = FALSE]
      }
    } else {
      ov.names <- names(data)
    }
  } else if (!is.null(sample.cov)) {
    ov.names <- rownames(sample.cov)
    if (is.null(ov.names)) {
      ov.names <- colnames(sample.cov)
    }
  }
  # ov.names?
  if (length(ov.names) == 0L) {
    lav_msg_stop(gettext(
      "could not extract variable names from data or sample.cov"))
  }

  # check nfactors
  if (any(nfactors < 1L)) {
    lav_msg_stop(gettext("nfactors must be greater than zero."))
  } else {
    # check for maximum number of factors
    # Fixme: can we do this more efficiently? also holds for categorical?
    nvar <- length(ov.names)
    p.star <- nvar * (nvar + 1) / 2
    nfac.max <- 0L
    for (nfac in seq_len(nvar)) {
      # compute number of free parameters
      npar <- nfac * nvar + nfac * (nfac + 1L) / 2 + nvar - nfac^2
      if (npar > p.star) {
        nfac.max <- nfac - 1L
        break
      }
    }
    if (any(nfactors > nfac.max)) {
      lav_msg_stop(gettextf("when nvar = %1$s the maximum number of factors
                            is %2$s", nvar, nfac.max))
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
    model.syntax <- lav_syntax_efa(
      ov.names = ov.names,
      nfactors = nfactors[f],
      twolevel = twolevel.flag
    )
    # call lavaan (using sem())
    FIT <- do.call("sem",
      args = c(
        list(
          model = model.syntax,
          data = data,
          sample.cov = sample.cov,
          sample.nobs = sample.nobs,
          rotation = rotation,
          rotation.args = rotation.args,
          bounds = bounds
        ),
        dotdotdot
      )
    )

    if (output == "efa") {
      FIT@Options$model.type <- "efa"
    }

    out[[f]] <- FIT
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
