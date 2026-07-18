# user-visible routine to
# compute polychoric/polyserial/... correlations
#
# YR 17 Sept 2013
#
# - YR 26 Nov 2013: big change - make it a wrapper around lavaan()
#                   estimator = "none" means two.step (starting values)

lav_object_cor <- function(object,
                   # lav.data options
                   ordered = NULL,
                   group = NULL,
                   missing = "listwise",
                   ov_names_x = NULL,
                   sampling_weights = NULL,
                   # lavaan options
                   se = "none",
                   test = "none",
                   estimator = "two.step",
                   baseline = FALSE,
                   # other options (for lavaan)
                   ...,
                   cor_smooth = FALSE,
                   cor_smooth_tol = 1e-04, # was 1e-06 in <0.6-14
                   output = "cor") {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, TRUE)
  # shortcut if object = lavaan object
  if (inherits(object, "lavaan")) {
    # check object
    object <- lav_object_check_version(object)
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

  # construct (and possibly fit) the unrestricted model; the workhorse
  # lav_h1_lavobject() is shared with lavH1()
  fit <- lav_h1_lavobject(
    object = object,
    ordered = ordered, group = group,
    missing = missing,
    ov_names_x = ov_names_x,
    sampling_weights = sampling_weights,
    estimator = estimator, se = se, test = test,
    baseline = baseline,
    dotdotdot = dotdotdot
  )

  out <- lav_object_cor_output(fit, output = output)

  # smooth correlation matrix? (only if output = "cor")
  if (output == "cor" && cor_smooth) {
    tmp_attr <- attributes(out)
    out <- cov2cor(lav_mat_sym_force_pd(out, tol = cor_smooth_tol))
    # we lost most of the attributes
    attributes(out) <- tmp_attr
  }

  out
}
lavCor <- lav_alias("lav_object_cor") # synonym #nolint

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
