# only compute 'H1' (=unrestricted/saturated) sample statistics
# (which usually go into the lavh1 slot) (and perhaps some extras, like icc)
#
# YR 14 June 2024: initial version

lavH1 <- function(object, # nolint
                  # lavdata options
                  ordered = NULL,
                  sampling.weights = NULL,
                  group = NULL,
                  cluster = NULL,
                  ov.names.x = character(0L),
                  ov.names.l = list(),
                  # other options (for lavaan)
                  ...,
                  output = "h1",
                  add.extra = TRUE,
                  # lavInspect() options
                  add.labels = TRUE,
                  add.class = TRUE,
                  drop.list.single.group = TRUE) {

  # shortcut if object = lavaan object
  if (inherits(object, "lavaan")) {
    lavh1 <- lavInspect(object, "h1")
    # if multilevel: add icc
    if (add.extra && object@Data@nlevels > 1L) {
      lavh1$icc <- lavInspect(object, "icc",
        add.labels = add.labels,
        add.class = add.class, drop.list.single.group = drop.list.single.group
      )
    }
    return(lavh1)
  }

  # output
  output <- tolower(output)

  # dotdotdot
  dotdotdot <- list(...)

  # verbose?
  if (!is.null(dotdotdot$verbose)) {
    if (dotdotdot$verbose) {
      lav_verbose(TRUE)
    }
    dotdotdot$verbose <- NULL
  }

  # check object class
  if (inherits(object, "lavData")) {
    lavdata <- object
  } else if (inherits(object, "data.frame") ||
    inherits(object, "matrix")) {
    object <- as.data.frame(object)
  } else {
    lav_msg_stop(
      gettext("lavH1 can not handle objects of class"),
      paste(class(object), collapse = " ")
    )
  }

  # prepare to create lavdata object
  if (inherits(object, "lavData")) {
    ov.names <- unique(unlist(lavdata@ov.names))
    ov.names.x <- unique(unlist(lavdata@ov.names.x))
    ordered <- lavdata@ordered
    ov.names.l <- lavdata@ov.names.l
    group <- lavdata@group
    cluster <- lavdata@cluster
    dotdotdot@missing <- lavdata@missing

    # just a regular data.frame
  } else {
    # ov.names
    ov.names <- names(object)
    if (!is.null(group)) {
      ov.names <- ov.names[-match(group, ov.names)]
    }
    if (!is.null(cluster)) {
      ov.names <- ov.names[-match(cluster, ov.names)]
    }
    if (!is.null(sampling.weights)) {
      ov.names <- ov.names[-match(sampling.weights, ov.names)]
    }
    if (!is.null(dotdotdot$conditional.x) && dotdotdot$conditional.x) {
      if (length(ov.names.x) == 0L) {
        lav_msg_stop(
          gettext("conditional.x = TRUE, but ov.names.x is empty")
        )
      }
      ov.names.x <- unique(unlist(ov.names.x))
      ov.names <- ov.names[-match(ov.names.x, ov.names)]
    }

    # ordered
    if (is.logical(ordered) && ordered) {
      ordered <- ov.names
    } else if (is.null(ordered) ||
      (length(ordered) == 1L && nchar(ordered) == 0L)) {
      # check data
      ordered <- ov.names[sapply(object, is.ordered)]
    } else if (!is.character(ordered)) {
      lav_msg_stop(gettext("ordered argument must be a character vector"))
    }

    if (length(ordered) > 0L) {
      # check if all names in "ordered" occur in the dataset?
      missing.idx <- which(!ordered %in% ov.names)
      if (length(missing.idx) > 0L) {
        lav_msg_warn(gettextf("ordered variable(s): %s could not be found
          in the data and will be ignored", lav_msg_view(ordered[missing.idx])))
      }
      # if ov.names.x, remove them from ordered (sanity check)
      if (length(ov.names.x) > 0L) {
        rm.idx <- which(ordered %in% ov.names.x)
        if (length(rm.idx) > 0L) {
          ordered <- ordered[-rm.idx]
        }
      }
    }

    # group/cluster
    if (!is.null(cluster)) {
      ngroups <- 1L
      if (!is.null(group)) {
        ngroups <- length(unique(object[, group]))
      }
      if (is.null(ov.names.l) || length(ov.names.l) == 0L) {
        ov.names.l <- vector("list", ngroups)
        # only two levels for now
        for (g in seq_len(ngroups)) {
          ov.names.l[[g]] <- vector("list", 2L)
          ov.names.l[[g]][[1]] <- ov.names
          ov.names.l[[g]][[2]] <- ov.names
        }
      } else {
        # check ov.names.l
        if (!is.list(ov.names.l)) { # just a single vector of names?
          ov.names.l <- rep(list(ov.names.l), 2L) # twolevels
        } else if (length(ov.names.l) == ngroups &&
          is.list(ov.names.l[[ngroups]])) {
        } else {
          ov.names.l <- rep(list(ov.names.l), ngroups)
        }
        tmp <- unique(unlist(ov.names.l))
        if (!all(tmp %in% ov.names)) {
          lav_msg_stop(gettextf(
            "some variable names in ov.names.l
            are missing: %1$s",
            lav_msg_view(tmp[!tmp %in% ov.names])
          ))
        }
      }
    }
  } # object = data

  # lavoptions
  lavoptions <- lav_lavaan_step02_options(
    slotOptions = NULL, # nolint
    slotData = NULL, # nolint
    flat.model = NULL,
    ordered = ordered,
    sample.cov = NULL,
    sample.mean = NULL,
    sample.th = NULL,
    sample.nobs = NULL,
    ov.names.l = ov.names.l,
    sampling.weights = sampling.weights,
    constraints = NULL,
    group = group,
    ov.names.x = ov.names.x,
    ov.names.y = ov.names[!ov.names %in% ov.names.x],
    dotdotdot = dotdotdot,
    cluster = cluster,
    data = "not.null"
  )

  # annoying: ov.names should not contain ov.names.x if categorical
  if (length(ordered) > 0L && lavoptions$conditional.x) {
    rm.idx <- which(ov.names %in% ov.names.x)
    if (length(rm.idx) > 0L) {
      ov.names <- ov.names[-rm.idx]
    }
  }

  # create lavdata (if needed)
  if (!inherits(object, "lavdata")) {
    lavdata <- lavData(
      data = object, group = group, cluster = cluster,
      ov.names = ov.names, ordered = ordered,
      sampling.weights = sampling.weights,
      ov.names.x = ov.names.x, ov.names.l = ov.names.l,
      lavoptions = lavoptions
    )
  }

  # create SampleStats (needed for multilevel - YLp slot)
  lavsamplestats <- lav_samplestats_from_data(
    lavdata = lavdata,
    lavoptions = lavoptions
  )

  # generate partable for unrestricted model
  lavpartable.un <- lav_partable_unrestricted(
    lavobject = NULL,
    lavdata = lavdata,
    lavsamplestats = lavsamplestats,
    lavh1 = NULL,
    lavoptions = lavoptions
  )

  # override some options
  lavoptions$h1 <- TRUE
  lavoptions$model.type <- "unrestricted"
  lavoptions$start <- "simple" # we have a full ustart column in lavpartable.un
  if (output != "lavaan") {
    lavoptions$implied <- FALSE
    lavoptions$loglik <- FALSE
    lavoptions$postcheck <- FALSE
    lavoptions$bounds <- "none"
    lavoptions$optim.method <- "none"
    lavoptions$estimator <- "none"
    lavoptions$se <- "none"
    lavoptions$test <- "none"
    lavoptions$baseline <- FALSE
  }

  # create lavaan object (usually without fitting)
  lavobject <- lavaan(
    slotParTable = lavpartable.un, slotData = lavdata,
    slotSampleStats = lavsamplestats, slotOptions = lavoptions
  )

  output <- tolower(output)
  if (output == "lavaan") {
    return(lavobject)
  } else {
    lavh1 <- lav_object_inspect_sampstat(
      object = lavobject, h1 = TRUE,
      add.labels = add.labels, add.class = add.class,
      drop.list.single.group = drop.list.single.group
    )

    # add additional info?

    # if multilevel: add icc
    if (add.extra && lavobject@Data@nlevels > 1L) {
      lavh1$icc <- lavInspect(lavobject, "icc",
        add.labels = add.labels,
        add.class = add.class, drop.list.single.group = drop.list.single.group
      )
    }
  }

  lavh1
}
