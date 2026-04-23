# user-visible routine to
# compute polychoric/polyserial/... correlations
#
# YR 17 Sept 2013
#
# - YR 26 Nov 2013: big change - make it a wrapper around lavaan()
#                   estimator = "none" means two.step (starting values)

lav_object_cor <- function(object,                 # nolint start
                   # lav.data options
                   ordered = NULL,
                   group = NULL,
                   missing = "listwise",
                   ov.names.x = NULL,
                   sampling.weights = NULL,
                   # lavaan options
                   se = "none",
                   test = "none",
                   estimator = "two.step",
                   baseline = FALSE,
                   # other options (for lavaan)
                   ...,
                   cor.smooth = FALSE,
                   cor.smooth.tol = 1e-04, # was 1e-06 in <0.6-14
                   output = "cor") {             # nolint end
  # shortcut if object = lavaan object
  if (inherits(object, "lavaan")) {
    # check object
    object <- lav_object_check_version(object)
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          "function lav_object_cor for lavaan-object")
        )
      }
    }
    out <- lav_object_cor_output(object, output = output)
    return(out)
  }

  # check estimator
  estimator <- tolower(estimator)
  if (estimator %in% c("two.step", "two.stage")) {
    estimator <- "none"
  }

  # se?
  se <- tolower(se)
  output <- tolower(output)
  if (se != "none") {
    if (output %in% c("cor", "cov", "sampstat", "th", "thresholds")) {
      lav_msg_warn(gettext("argument `se' is ignored since standard errors
                           are not needed for the requested `output'"))
      se <- "none"
    }
  }

  # extract sampling.weights.normalization from dots (for lav_lavdata() call)
  dots <- list(...)
  sampling_weights_normalization <- "total"
  if (!is.null(dots$sampling.weights.normalization)) {
    sampling_weights_normalization <- dots$sampling.weights.normalization
  }


  # check object class
  if (inherits(object, "lavData")) {
    lav_data <- object
  } else if (inherits(object, "data.frame") ||
    inherits(object, "matrix")) {
    object <- as.data.frame(object)
    names_1 <- names(object)
    if (!is.null(group)) {
      names_1 <- names_1[-match(group, names_1)]
    }
    if (!is.null(sampling.weights)) {
      names_1 <- names_1[-match(sampling.weights, names_1)]
    }
    if (is.logical(ordered)) {
      ordered_flag <- ordered
      if (ordered_flag) {
        ordered <- names_1
        if (length(ov.names.x) > 0L) {
          ordered <- ordered[-which(ordered %in% ov.names.x)]
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
    lav_data <- lav_lavdata(
      data = object, group = group,
      ov.names = names_1, ordered = ordered,
      sampling.weights = sampling.weights,
      ov.names.x = ov.names.x,
      lavoptions = list(
        missing = missing,
        sampling.weights.normalization = sampling_weights_normalization
      )
    )
  } else {
    lav_msg_stop(gettext("lav_object_cor can not handle objects of class"),
      paste(class(object), collapse = " ")
    )
  }

  # set default estimator if se != "none"
  categorical <- any(lav_data@ov$type == "ordered")
  if (se != "none" && estimator == "none") {
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
  if (!is.null(dots$meanstructure)) {
    meanstructure <- dots$meanstructure
  }
  if (lav_data@ngroups > 1L || categorical ||
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

  # override, only for backwards compatibility (eg moments() in JWileymisc)
  # if(missing %in% c("ml", "fiml")) {
  #    meanstructure = TRUE
  # }

  # generate partable for unrestricted model
  pt_un <-
    lav_partable_unrestricted(
      lavobject = NULL,
      lavdata = lav_data,
      lavoptions = list(
        meanstructure = meanstructure,
        fixed.x = fixed_x,
        conditional.x = conditional_x,
        # sampling.weights.normalization = sampling.weights.normalization,
        group.w.free = FALSE,
        missing = missing,
        estimator = estimator,
        mimic = mimic
      ),
      sample.cov = NULL,
      sample.mean = NULL,
      sample.th = NULL
    )


  fit <- lavaan(
    slotParTable = pt_un, slotData = lav_data,
    model.type = "unrestricted",
    missing = missing,
    baseline = baseline, h1 = TRUE, # must be TRUE!
    se = se, test = test, estimator = estimator, ...
  )

  out <- lav_object_cor_output(fit, output = output)

  # smooth correlation matrix? (only if output = "cor")
  if (output == "cor" && cor.smooth) {
    tmp_attr <- attributes(out)
    out <- cov2cor(lav_matrix_symmetric_force_pd(out, tol = cor.smooth.tol))
    # we lost most of the attributes
    attributes(out) <- tmp_attr
  }

  out
}
lavCor <- lav_object_cor  # synonym #nolint

lav_object_cor_output <- function(object, output = "cor") {
  # check output
  if (output %in% c("cor", "cov")) {
    out <- lavInspect(object, "sampstat")
    if (object@Data@ngroups == 1L) {
      if (object@Model@conditional.x) {
        out <- out$res.cov
      } else {
        out <- out$cov
      }
      if (output == "cor") {
        out <- cov2cor(out)
      }
    } else {
      if (object@Model@conditional.x) {
        out <- lapply(out, "[[", "res.cov")
      } else {
        out <- lapply(out, "[[", "cov")
      }
      if (output == "cor") {
        out <- lapply(out, cov2cor)
      }
    }
  } else if (output %in% c("th", "thresholds")) {
    out <- lavInspect(object, "sampstat")
    if (object@Data@ngroups == 1L) {
      if (object@Model@conditional.x) {
        out <- out$res.th
      } else {
        out <- out$th
      }
    } else {
      if (object@Model@conditional.x) {
        out <- lapply(out, "[[", "res.th")
      } else {
        out <- lapply(out, "[[", "th")
      }
    }
  } else if (output %in% c("sampstat")) {
    out <- lavInspect(object, "sampstat")
  } else if (output %in% c(
    "parameterEstimates", "pe",
    "parameterestimates", "est"
  )) {
    out <- standardizedSolution(object)
  } else {
    out <- object
  }

  out
}
