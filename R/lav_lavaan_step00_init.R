lav_lavaan_step00_parameters <- function(matchcall = NULL,
                                         syscall = NULL,
                                         dotdotdot = NULL) {
  # 1. to resolve a problem where parameter 'cl' is matched to 'cluster'
  #    and shouldn't
  # 2. to apply defaut options for cfa/sem/growth functions
  # 3. if dotdotdot$control present, copy to dotdotdot$... for
  #   optim.method,
  #   optim.force.converged,
  #   gradient -> optim.gradient !!! overwritten by dotdotdot$gradient
  #     if present
  #   init_nelder_mead -> optim.init_nelder_mead

  mc <- matchcall
  sc <- syscall
  ddd <- dotdotdot

  # catch partial matching of 'cl' (expanded to cluster)
  if (!is.null(sc[["cl"]]) &&
    is.null(sc[["cluster"]]) &&
    !is.null(mc[["cluster"]])) {
    mc[["cl"]] <- mc[["cluster"]]
    mc[["cluster"]] <- NULL
    ddd$cl <- sc[["cl"]]
  }
  if (!is.null(mc$cluster)) mc$cluster <- eval(mc$cluster, parent.frame(2))

  # default options
  if (any(ddd$model.type == c("sem", "cfa", "growth"))) {
    # default options for sem/cfa or growth
    defaults <- list(
      int.ov.free     = ddd$model.type != "growth",
      int.lv.free     = ddd$model.type == "growth",
      auto.fix.first  = TRUE, # (re)set in lav_options_set
      auto.fix.single = TRUE,
      auto.var        = TRUE,
      auto.cov.lv.x   = TRUE,
      auto.cov.y      = TRUE,
      auto.th         = TRUE,
      auto.delta      = TRUE,
      auto.efa        = TRUE
    )
    for (dflt.i in seq_along(defaults)) {
      argname <- names(defaults)[dflt.i]
      if (is.null(mc[[argname]])) {
        mc[[argname]] <- defaults[[dflt.i]]
        ddd[[argname]] <- defaults[[dflt.i]]
      }
    }
  }

  # backwards compatibility, control= argument (<0.5-23)
  if (!is.null(ddd$control)) {
    # optim.method
    if (!is.null(ddd$control$optim.method)) {
      ddd$optim.method <- ddd$control$optim.method
    }
    # cor.optim.method
    if (!is.null(ddd$control$cor.optim.method)) {
      # ignore it silently
    }
    # control$optim.force.converged
    if (!is.null(ddd$control$optim.force.converged)) {
      ddd$optim.force.converged <- ddd$control$optim.force.converged
    }
    # gradient
    if (!is.null(ddd$control$gradient)) {
      ddd$optim.gradient <- ddd$control$gradient
    }
    if (!is.null(ddd$gradient)) {
      ddd$optim.gradient <- ddd$gradient
    }
    # init_nelder_mead
    if (!is.null(ddd$control$init_nelder_mead)) {
      ddd$optim.init_nelder_mead <- ddd$control$init_nelder_mead
    }
  }

  list(mc = mc, dotdotdot = ddd)
}

lav_lavaan_step00_checkdata <- function(data = NULL,
                                        dotdotdot = NULL,
                                        sample.cov = NULL,
                                        sample.nobs = NULL,
                                        sample.mean = NULL,
                                        sample.th = NULL,
                                        NACOV = NULL, # nolint
                                        WLS.V = NULL, # nolint
                                        ov.order = NULL) {
  # if data not NULL:
  #   if it is an 'enriched' data.frame (e.g. a tibble), simplify to an
  #   ordinary data.frame
  #   if class is 'lavMoments':
  #         check if it contains sample.cov and sample.nobs (***error*** if not)
  #         if sample.mean, sample.th, NACOV and/or WLS.V present,
  #                   copy to corresponding arguments of the function
  #         if lavOptions present in data, copy those that are not provided
  #         in the function call to dotdotdot$...
  #             set data to NULL
  #   if it is a function --> ***error***
  #   TODO: other tests are present in lavData, should we copy them here ???
  # if NACOV or WLS.V not NULL, set ov.order to "data"

  if (!is.null(data)) {
    if (inherits(data, "data.frame")) {
      # just in case it is not a traditional data.frame
      data <- as.data.frame(data)
    } else if (inherits(data, "lavMoments")) {
      # This object must contain summary statistics
      # e.g., created by lavaan.mi::poolSat

      # set required-data arguments
      if ("sample.cov" %in% names(data)) {
        sample.cov <- data$sample.cov
      } else {
        lav_msg_stop(gettext(
          "When data= is of class lavMoments, it must contain sample.cov"))
      }

      if ("sample.nobs" %in% names(data)) {
        sample.nobs <- data$sample.nobs
      } else {
        lav_msg_stop(gettext(
          "When data= is of class lavMoments, it must contain sample.nobs"))
      }

      # check for optional-data arguments
      if ("sample.mean" %in% names(data)) sample.mean <- data$sample.mean
      if ("sample.th" %in% names(data)) sample.th <- data$sample.th
      if ("NACOV" %in% names(data)) NACOV <- data$NACOV # nolint
      if ("WLS.V" %in% names(data)) WLS.V <- data$WLS.V # nolint

      # set other args not included in dotdotdot
      if (length(data$lavOptions)) {
        newdots <- setdiff(names(data$lavOptions), names(dotdotdot))
        if (length(newdots)) {
          for (dd in newdots) dotdotdot[[dd]] <- data$lavOptions[[dd]]
        }
      }

      # FIXME: Should WLS.V be an I(dentity) matrix when ULS is requested?
      #       Unused for point estimates, but still used to scale/shift test
      # if (!is.null(dotdotdot$estimator)) {
      #   if (grepl(pattern = "ULS", x = toupper(dotdotdot$estimator[1L])) &&
      #       !is.null(WLS.V)) {
      #     #  set to diagonal
      #     if (is.list(WLS.V)) {
      #       WLS.V <- lapply(WLS.V, function(w) {diag(w) <- 1 ; return(w) })
      #     } else diag(WLS.V) <- 1
      #   }
      # }

      # get rid of data= argument
      data <- NULL
    }

    if (is.function(data)) {
      lav_msg_stop(gettext("data is a function; it should be a data.frame"))
    }
  }
  # new in 0.6-14: if NACOV and/or WLS.V are provided, we force
  # ov.order="data" for now
  # until we have reliable code to re-arrange/select col/rows for
  # of NACOV/WLS.V based on the model-based ov.names
  if (!is.null(NACOV) || !is.null(WLS.V)) {
    ov.order <- "data"
  }

  list(
    data = data, dotdotdot = dotdotdot, sample.cov = sample.cov,
    sample.nobs = sample.nobs, sample.mean = sample.mean,
    sample.th = sample.th, NACOV = NACOV, WLS.V = WLS.V, ov.order = ov.order
  )
}
