#
# the lavData class describes how the data looks like
#  - do we have a full data frame, or only sample statistics?
#    (TODO: allow for patterns + freq, if data is categorical)
#  - variable type ("numeric", "ordered", ...)
#  - how many groups, how many observations, ...
#  - what about missing patterns?
#
# initial version: YR 14 April 2012

# YR 23 Feb 2017: blocks/levels/groups, but everything is group-based!

# FIXME: if nlevels > 1L, and ngroups > 1L, we should check that
# group is at the upper-level

# YR 08 May  2019: sampling weights normalization -> different options

# extract the data we need for this particular model
lavData <- function(data = NULL, # data.frame
                    group = NULL, # multiple groups?
                    cluster = NULL, # clusters?
                    ov.names = NULL, # variables in model
                    ov.names.x = character(0), # exo variables
                    ov.names.l = list(), # names per level
                    ordered = NULL, # ordered variables
                    sampling.weights = NULL, # sampling weights
                    sample.cov = NULL, # sample covariance(s)
                    sample.mean = NULL, # sample mean vector(s)
                    sample.th = NULL, # sample thresholds
                    sample.nobs = NULL, # sample nobs
                    lavoptions = lavOptions(), # lavoptions
                    allow.single.case = FALSE # for newdata in predict
) {
  # get info from lavoptions

  # group.labels
  group.label <- lavoptions$group.label
  if (is.null(group.label)) {
    group.label <- character(0L)
  }

  # level.labels
  level.label <- lavoptions$level.label
  if (is.null(level.label)) {
    level.label <- character(0L)
  }

  # block.labels
  block.label <- character(0L)
  if (length(group.label) > 0L && length(level.label) == 0L) {
    block.label <- group.label
  } else if (length(level.label) > 0L && length(group.label) == 0L) {
    block.label <- level.label
  } else if (length(group.label) > 0L &&
    length(level.label) > 0L) {
    block.label <- paste(rep(group.label, each = length(level.label)),
      rep(level.label, times = length(group.label)),
      sep = "."
    )
  }

  # std.ov?
  std.ov <- lavoptions$std.ov
  if (is.null(std.ov)) {
    std.ov <- FALSE
  }

  # missing?
  missing <- lavoptions$missing
  if (is.null(missing) || missing == "default") {
    missing <- "listwise"
  }

  # warn?
  warn <- lavoptions$warn
  if (is.null(warn)) {
    warn <- TRUE
  }
  if (allow.single.case) { # eg, in lavPredict
    warn <- FALSE
  }

  # four scenarios:
  #    0) data is already a lavData object: do nothing
  #    1) data is full data.frame (or a matrix)
  #    2) data are sample statistics only
  #    3) no data at all

  # 1) full data
  if (!is.null(data)) {
    # catch lavaan/lavData objects
    if (inherits(data, "lavData")) {
      return(data)
    } else if (inherits(data, "lavaan")) {
      return(data@Data)
    }

    # catch matrix
    if (!is.data.frame(data)) {
      # is it a matrix?
      if (is.matrix(data)) {
        if (nrow(data) == ncol(data)) {
          # perhaps it is a covariance matrix?
          if (isSymmetric(data)) {
            lav_msg_warn(
              gettext("data argument looks like a covariance matrix;
              please use the sample.cov argument instead"))
          }
        }
        # or perhaps it is a data matrix?
        ### FIXME, we should avoid as.data.frame() and handle
        ### data matrices directly
        data <- as.data.frame(data, stringsAsFactors = FALSE)
      } else {
        lav_msg_stop(gettextf(
          "data= argument is not a data.frame, but of class ",
          sQuote(class(data))))
      }
    }

    # no ov.names?
    if (is.null(ov.names)) {
      ov.names <- names(data)
      # remove group variable, if provided
      if (length(group) > 0L) {
        group.idx <- which(ov.names == group)
        ov.names <- ov.names[-group.idx]
      }
      # remove cluster variable, if provided
      if (length(cluster) > 0L) {
        cluster.idx <- which(ov.names == cluster)
        ov.names <- ov.names[-cluster.idx]
      }
    }

    lavData <- lav_data_full(
      data = data,
      group = group,
      cluster = cluster,
      group.label = group.label,
      level.label = level.label,
      block.label = block.label,
      ov.names = ov.names,
      ordered = ordered,
      sampling.weights = sampling.weights,
      sampling.weights.normalization =
        lavoptions$sampling.weights.normalization,
      ov.names.x = ov.names.x,
      ov.names.l = ov.names.l,
      std.ov = std.ov,
      missing = missing,
      warn = warn,
      allow.single.case = allow.single.case
    )
    sample.cov <- NULL # not needed, but just in case
  }


  # 2) sample moments
  if (is.null(data) && !is.null(sample.cov)) {
    # for now: no levels!!
    nlevels <- 1L

    # we also need the number of observations (per group)
    if (is.null(sample.nobs)) {
      lav_msg_stop(gettext("please specify number of observations"))
    }

    # if a 'group' argument was provided, keep it -- new in 0.6-4
    if (is.null(group)) {
      group <- character(0L)
    } else if (is.character(group)) {
      # nothing to do, just store it
    } else {
      lav_msg_stop(gettext("group argument should be a string"))
    }

    # list?
    if (is.list(sample.cov)) {
      # multiple groups, multiple cov matrices
      if (!is.null(sample.mean)) {
        stopifnot(length(sample.mean) == length(sample.cov))
      }
      if (!is.null(sample.th)) {
        stopifnot(length(sample.th) == length(sample.cov))
      }
      # multiple groups, multiple cov matrices
      ngroups <- length(sample.cov)
      LABEL <- names(sample.cov)
      if (is.null(group.label) || length(group.label) == 0L) {
        if (is.null(LABEL)) {
          group.label <- paste("Group ", 1:ngroups, sep = "")
        } else {
          group.label <- LABEL
        }
      } else {
        if (is.null(LABEL)) {
          stopifnot(length(group.label) == ngroups)
        } else {
          # FIXME!!!!
          # check if they match
        }
      }
    } else {
      ngroups <- 1L
      group.label <- character(0)
      if (!is.matrix(sample.cov)) {
        lav_msg_stop(gettext(
          "sample.cov must be a matrix or a list of matrices"))
      }
      sample.cov <- list(sample.cov)
    }

    # get ov.names
    if (is.null(ov.names)) {
      ov.names <- lapply(sample.cov, row.names)
    } else if (!is.list(ov.names)) {
      # duplicate ov.names for each group
      tmp <- ov.names
      ov.names <- vector("list", length = ngroups)
      ov.names[1:ngroups] <- list(tmp)
    } else {
      if (length(ov.names) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names assumes %1$s groups; data contains %2$s groups",
          length(ov.names), ngroups))
      }
      # nothing to do
    }

    # handle ov.names.x
    if (!is.list(ov.names.x)) {
      tmp <- ov.names.x
      ov.names.x <- vector("list", length = ngroups)
      ov.names.x[1:ngroups] <- list(tmp)
    } else {
      if (length(ov.names.x) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names.x assumes %1$s groups; data contains %2$s groups",
          length(ov.names.x), ngroups))
      }
    }

    ov <- list()
    ov$name <- unique(unlist(c(ov.names, ov.names.x)))
    nvar <- length(ov$name)
    ov$idx <- rep(NA, nvar)
    ov$nobs <- rep(sum(unlist(sample.nobs)), nvar)
    ov$type <- rep("numeric", nvar)
    ov$nlev <- rep(0, nvar)
    # check for categorical
    if (!is.null(sample.th)) {
      th.idx <- attr(sample.th, "th.idx")
      if (is.list(th.idx)) {
        th.idx <- th.idx[[1]] ## FIRST group only (assuming same ths!)
      }
      if (any(th.idx > 0)) {
        TAB <- table(th.idx[th.idx > 0])
        ord.idx <- as.numeric(names(TAB))
        nlev <- as.integer(unname(TAB) + 1)
        ov$type[ord.idx] <- "ordered"
        ov$nlev[ord.idx] <- nlev
      }
    }

    # if std.ov = TRUE, give a warning (suggested by Peter Westfall)
    if (std.ov && warn) {
      lav_msg_warn(gettext(
        "std.ov argument is ignored if only sample statistics are provided."))
    }

    # check variances (new in 0.6-7)
    for (g in seq_len(ngroups)) {
      VAR <- diag(sample.cov[[g]])
      # 1. finite?
      if (!all(is.finite(VAR))) {
        lav_msg_stop(gettext(
        "at least one variance in the sample covariance matrix is not finite."))
        }
      # 2. near zero (or negative)?
      if (any(VAR < .Machine$double.eps)) {
        lav_msg_stop(
          gettext("at least one variance in the sample covariance matrix is
          (near) zero or negative."))
      }
      # 3. very large?
      max.var <- max(VAR)
      if (max.var > 1000000) {
        lav_msg_warn(
          gettext("some observed variances in the sample covariance matrix
          are larger than 1000000."))
      }
    }

    # block.labels
    block.label <- character(0L)
    if (length(group.label) > 0L && length(level.label) == 0L) {
      block.label <- group.label
    } else if (length(level.label) > 0L && length(group.label) == 0L) {
      block.label <- level.label
    } else if (length(group.label) > 0L &&
      length(level.label) > 0L) {
      block.label <- paste(rep(group.label, each = length(level.label)),
        rep(level.label, times = length(group.label)),
        sep = "."
      )
    }

    # construct lavData object
    lavData <- new("lavData",
      data.type = "moment",
      ngroups = ngroups,
      group = group,
      nlevels = 1L, # for now
      cluster = character(0L),
      group.label = group.label,
      level.label = character(0L),
      block.label = block.label,
      nobs = as.list(sample.nobs),
      norig = as.list(sample.nobs),
      ov.names = ov.names,
      ov.names.x = ov.names.x,
      ov.names.l = ov.names.l,
      ordered = as.character(ordered),
      weights = vector("list", length = ngroups),
      sampling.weights = character(0L),
      ov = ov,
      std.ov = FALSE,
      missing = "listwise",
      case.idx = vector("list", length = ngroups),
      Mp = vector("list", length = ngroups),
      Rp = vector("list", length = ngroups),
      Lp = vector("list", length = ngroups),
      X = vector("list", length = ngroups),
      eXo = vector("list", length = ngroups)
    )
  }

  # 3) data.type = "none":  both data and sample.cov are NULL
  if (is.null(data) && is.null(sample.cov)) {
    # clustered/multilevel? --> ov.names.l should be filled in
    if (length(ov.names.l) > 0L) {
      nlevels <- length(ov.names.l[[1]]) # we assume the same number
      # of levels in each group!

      # do we have a cluster argument? if not, create one
      if (is.null(cluster)) {
        if (nlevels == 2L) {
          cluster <- "cluster"
        } else {
          cluster <- paste0("cluster", seq_len(nlevels - 1L))
        }
      }

      # default level.labels
      if (length(level.label) == 0L) {
        level.label <- c("within", cluster)
      } else {
        # check if length(level.label) = 1 + length(cluster)
        if (length(level.label) != length(cluster) + 1L) {
          lav_msg_stop(gettext("length(level.label) != length(cluster) + 1L"))
        }
        # nothing to do
      }
    } else {
      nlevels <- 1L
      cluster <- character(0L)
      level.label <- character(0L)
    }

    # ngroups: ov.names (when group: is used), or sample.nobs
    if (is.null(ov.names)) {
      lav_msg_warn(gettext("ov.names is NULL"))
      ov.names <- character(0L)
      if (is.null(sample.nobs)) {
        ngroups <- 1L
        sample.nobs <- rep(list(0L), ngroups)
      } else {
        sample.nobs <- as.list(sample.nobs)
        ngroups <- length(sample.nobs)
      }
    } else if (!is.list(ov.names)) {
      if (is.null(sample.nobs)) {
        ngroups <- 1L
        sample.nobs <- rep(list(0L), ngroups)
      } else {
        sample.nobs <- as.list(sample.nobs)
        ngroups <- length(sample.nobs)
      }
      ov.names <- rep(list(ov.names), ngroups)
    } else if (is.list(ov.names)) {
      ngroups <- length(ov.names)
      if (is.null(sample.nobs)) {
        sample.nobs <- rep(list(0L), ngroups)
      } else {
        sample.nobs <- as.list(sample.nobs)
        if (length(sample.nobs) != ngroups) {
          lav_msg_stop(gettextf(
            "length(sample.nobs) = %1$s but syntax implies ngroups = %2$s",
            length(sample.nobs), ngroups))
        }
      }
    }


    # group.label
    if (ngroups > 1L) {
      if (is.null(group)) {
        group <- "group"
      }
      group.label <- paste("Group", 1:ngroups, sep = "")
    } else {
      group <- character(0L)
      group.label <- character(0L)
    }

    # handle ov.names.x
    if (!is.list(ov.names.x)) {
      ov.names.x <- rep(list(ov.names.x), ngroups)
    }

    ov <- list()
    ov$name <- unique(unlist(c(ov.names, ov.names.x)))
    nvar <- length(ov$name)
    ov$idx <- rep(NA, nvar)
    ov$nobs <- rep(0L, nvar)
    ov$type <- rep("numeric", nvar)
    ov$nlev <- rep(0L, nvar)

    # collect information per upper-level group
    Lp <- vector("list", length = ngroups)
    for (g in 1:ngroups) {
      if (nlevels > 1L) {
        # ALWAYS add ov.names.x at the end, even if conditional.x
        OV.NAMES <- unique(c(ov.names[[g]], ov.names.x[[g]]))
        Lp[[g]] <- lav_data_cluster_patterns(
          Y = NULL, clus = NULL,
          cluster = cluster,
          multilevel = TRUE,
          ov.names = OV.NAMES,
          ov.names.x = ov.names.x[[g]],
          ov.names.l = ov.names.l[[g]]
        )
      }
    } # g

    # block.labels
    block.label <- character(0L)
    if (length(group.label) > 0L && length(level.label) == 0L) {
      block.label <- group.label
    } else if (length(level.label) > 0L && length(group.label) == 0L) {
      block.label <- level.label
    } else if (length(group.label) > 0L &&
      length(level.label) > 0L) {
      block.label <- paste(rep(group.label, each = length(level.label)),
        rep(level.label, times = length(group.label)),
        sep = "."
      )
    }

    # construct lavData object
    lavData <- new("lavData",
      data.type = "none",
      ngroups = ngroups,
      group = group,
      nlevels = nlevels,
      cluster = cluster,
      group.label = group.label,
      level.label = level.label,
      block.label = block.label,
      nobs = sample.nobs,
      norig = sample.nobs,
      ov.names = ov.names,
      ov.names.x = ov.names.x,
      ov.names.l = ov.names.l,
      ordered = as.character(ordered),
      weights = vector("list", length = ngroups),
      sampling.weights = character(0L),
      ov = ov,
      missing = "listwise",
      case.idx = vector("list", length = ngroups),
      Mp = vector("list", length = ngroups),
      Rp = vector("list", length = ngroups),
      Lp = Lp,
      X = vector("list", length = ngroups),
      eXo = vector("list", length = ngroups)
    )
  }

  lavData
}


# handle full data
lav_data_full <- function(data = NULL, # data.frame
                          group = NULL, # multiple groups?
                          cluster = NULL, # clustered?
                          group.label = NULL, # custom group labels?
                          level.label = NULL,
                          block.label = NULL,
                          ov.names = NULL, # variables needed
                          # in model
                          ordered = NULL, # ordered variables
                          sampling.weights = NULL, # sampling weights
                          sampling.weights.normalization = "none",
                          ov.names.x = character(0L), # exo variables
                          ov.names.l = list(), # var per level
                          std.ov = FALSE, # standardize ov's?
                          missing = "listwise", # remove missings?
                          warn = TRUE, # produce warnings?
                          allow.single.case = FALSE # allow single case?
) {
  # number of groups and group labels
  if (!is.null(group) && length(group) > 0L) {
    if (!(group %in% names(data))) {
      lav_msg_stop(gettextf(
        "grouping variable %1$s not found; variable names
        found in data frame are: %2$s",
        sQuote(group),paste(names(data), collapse = " ")))
    }
    # note: by default, we use the order as in the data;
    # not as in levels(data[,group])
    if (length(group.label) == 0L) {
      group.label <- unique(as.character(data[[group]]))
      if (warn && any(is.na(group.label))) {
        lav_msg_warn(gettextf("group variable %s contains missing values",
                              sQuote(group)))
      }
      group.label <- group.label[!is.na(group.label)]
    } else {
      group.label <- unique(as.character(group.label))
      # check if user-provided group labels exist
      LABEL <- unique(as.character(data[[group]]))
      idx <- match(group.label, LABEL)
      if (warn && any(is.na(idx))) {
        lav_msg_warn(gettextf(
          "some group.labels do not appear in the grouping variable: %s"),
          lav_msg_view(group.label[which(is.na(idx))], log.sep = "none")
        )
      }
      group.label <- group.label[!is.na(idx)]
      # any groups left?
      if (length(group.label) == 0L) {
        lav_msg_stop(gettext(
          "no group levels left; check the group.label argument"))
      }
    }
    ngroups <- length(group.label)
  } else {
    if (warn && length(group.label) > 0L) {
      lav_msg_warn(gettext(
       "`group.label' argument will be ignored if `group' argument is missing"))
    }
    ngroups <- 1L
    group.label <- character(0L)
    group <- character(0L)
  }

  # sampling weights
  if (!is.null(sampling.weights)) {
    if (is.character(sampling.weights)) {
      if (!(sampling.weights %in% names(data))) {
        lav_msg_stop(
          gettextf("sampling weights variable %1$s not found;
                   variable names found in data frame are: %2$s",
          sQuote(sampling.weights)), paste(names(data), collapse = " "))
      }
      # check for missing values in sampling weight variable
      if (any(is.na(data[[sampling.weights]]))) {
        lav_msg_stop(
          gettextf("sampling.weights variable %s contains missing values",
          sQuote(sampling.weights)))
      }
    } else {
      lav_msg_stop(gettext(
        "sampling weights argument should be a variable name in the data.frame"
        ))
    }
  }

  # clustered?
  if (!is.null(cluster) && length(cluster) > 0L) {
    # cluster variable in data?
    if (!all(cluster %in% names(data))) {
      # which one did we not find?
      not.ok <- which(!cluster %in% names(data))
      lav_msg_stop(gettextf(
        "cluster variable(s) %1$s not found;
        variable names found in data frame are: %2$s",
        sQuote(cluster[not.ok]), paste(names(data), collapse = " ")))
    }

    # check for missing values in cluster variable(s)
    for (cl in 1:length(cluster)) {
      if (warn && anyNA(data[[cluster[cl]]])) {
        lav_msg_warn(gettextf("cluster variable %s contains missing values",
          sQuote(cluster[cl])))
      }
    }

    # multilevel?
    if (length(ov.names.l) > 0L) {
      # default level.labels
      if (length(level.label) == 0L) {
        level.label <- c("within", cluster)
      } else {
        # check if length(level.label) = 1 + length(cluster)
        if (length(level.label) != length(cluster) + 1L) {
          lav_msg_stop(gettext("length(level.label) != length(cluster) + 1L"))
        }
        # nothing to do
      }
      nlevels <- length(level.label)
    } else {
      # just clustered data, but no random effects
      nlevels <- 1L
      level.label <- character(0L)
    }
  } else {
    if (warn && length(level.label) > 0L) {
      lav_msg_warn(gettext(
       "`level.label' argument will be ignored if `cluster' argument is missing"
      ))
    }
    nlevels <- 1L
    level.label <- character(0L)
    cluster <- character(0L)
  }

  # check ov.names vs ngroups
  if (ngroups > 1L) {
    if (is.list(ov.names)) {
      if (length(ov.names) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names assumes %1$s groups; data contains %2$s groups",
          length(ov.names), ngroups))
      }
    } else {
      tmp <- ov.names
      ov.names <- vector("list", length = ngroups)
      ov.names[1:ngroups] <- list(tmp)
    }
    if (is.list(ov.names.x)) {
      if (length(ov.names.x) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names.x assumes %1$s groups; data contains %2$s groups",
          length(ov.names.x), ngroups))
      }
    } else {
      tmp <- ov.names.x
      ov.names.x <- vector("list", length = ngroups)
      ov.names.x[1:ngroups] <- list(tmp)
    }
  } else {
    if (is.list(ov.names)) {
      if (length(ov.names) > 1L) {
        lav_msg_stop(gettext(
          "model syntax defines multiple groups; data suggests a single group"))
      }
    } else {
      ov.names <- list(ov.names)
    }
    if (is.list(ov.names.x)) {
      if (length(ov.names.x) > 1L) {
        lav_msg_stop(gettext(
          "model syntax defines multiple groups; data suggests a single group"))
      }
    } else {
      ov.names.x <- list(ov.names.x)
    }
  }

  # check if all ov.names can be found in the data.frame
  for (g in 1:ngroups) {
    # does the data contain all the observed variables
    # needed in the user-specified model for this group
    ov.all <- unique(c(ov.names[[g]], ov.names.x[[g]])) # no overlap if categ

    # handle interactions
    ov.int.names <- ov.all[grepl(":", ov.all)]
    n.int <- length(ov.int.names)
    if (n.int > 0L) {
      ov.names.noint <- ov.all[!ov.all %in% ov.int.names]
      for (iv in seq_len(n.int)) {
        NAMES <- strsplit(ov.int.names[iv], ":", fixed = TRUE)[[1L]]
        if (all(NAMES %in% ov.names.noint)) {
          # add this interaction term to the data.frame, unless
          # it already exists
          if (is.null(data[[ov.int.names[iv]]])) {
            data[[ov.int.names[iv]]] <-
              data[[NAMES[1L]]] * data[[NAMES[2L]]]
          }
        }
      }
    }

    # check for missing observed variables
    idx.missing <- which(!(ov.all %in% names(data)))

    if (length(idx.missing)) {
      lav_msg_stop(
        gettextf("some (observed) variables specified in the model are
        not found in the dataset: %s",
        paste(ov.all[idx.missing], collapse = " ")))
    }
  }


  # here, we know for sure all ov.names exist in the data.frame
  # create varTable
  # FIXME: should we add the 'group'/'cluster' variable (no for now)
  ov <- lav_dataframe_vartable(
    frame = data, ov.names = ov.names,
    ov.names.x = ov.names.x, ordered = ordered,
    as.data.frame. = FALSE
  )

  # do some checking
  # check for unordered factors (but only if nlev > 2)
  if ("factor" %in% ov$type) {
    f.names <- ov$name[ov$type == "factor" & ov$nlev > 2L]
    f.names.all <- ov$name[ov$type == "factor"]
    OV.names <- unlist(ov.names)
    OV.names.x <- unlist(ov.names.x)
    OV.names.nox <- OV.names[!OV.names %in% OV.names.x]
    if (any(f.names %in% OV.names.x)) {
      lav_msg_stop(
        gettext("unordered factor(s) with more than 2 levels detected
            as exogenous covariate(s): "),
        paste(f.names, collapse = " "))
    } else if (any(f.names.all %in% OV.names.nox)) {
      lav_msg_stop(
        gettext("unordered factor(s) detected; make them numeric or ordered:"),
        paste(f.names.all, collapse = " "))
    }
  }
  # check for ordered exogenous variables
  if ("ordered" %in% ov$type[ov$name %in% unlist(ov.names.x)]) {
    f.names <- ov$name[ov$type == "ordered" &
      ov$name %in% unlist(ov.names.x)]
    if (warn && any(f.names %in% unlist(ov.names.x))) {
      lav_msg_warn(gettextf(
        "exogenous variable(s) declared as ordered in data: %s",
        lav_msg_view(f.names, log.sep = "none")))
    }
  }
  # check for ordered endogenous variables with more than 12 levels
  if ("ordered" %in% ov$type[!ov$name %in% unlist(ov.names.x)]) {
    f.names <- ov$name[ov$type == "ordered" &
      !ov$name %in% unlist(ov.names.x) &
      ov$nlev > 12L]
    if (warn && length(f.names) > 0L) {
      lav_msg_warn(gettextf(
        "some ordered categorical variable(s) have more than 12 levels: %s",
        lav_msg_view(f.names, log.sep = "none")))
    }
  }
  # check for zero-cases
  idx <- which(ov$nobs == 0L | ov$var == 0)
  if (!allow.single.case && length(idx) > 0L) {
    OV <- as.data.frame(ov)
    rn <- rownames(OV)
    rn[idx] <- paste(rn[idx], "***", sep = "")
    rownames(OV) <- rn
    print(OV)
    lav_msg_stop(gettext(
      "some variables have no values (only missings) or no variance"))
  }
  # check for single cases (no variance!)
  idx <- which(ov$nobs == 1L | (ov$type == "numeric" & !is.finite(ov$var)))
  if (!allow.single.case && length(idx) > 0L) {
    OV <- as.data.frame(ov)
    rn <- rownames(OV)
    rn[idx] <- paste(rn[idx], "***", sep = "")
    rownames(OV) <- rn
    print(OV)
    lav_msg_stop(gettext(
      "some variables have only 1 observation or no finite variance"))
  }
  # check for ordered variables with only 1 level
  idx <- which(ov$type == "ordered" & ov$nlev == 1L)
  if (!allow.single.case && length(idx) > 0L) {
    OV <- as.data.frame(ov)
    rn <- rownames(OV)
    rn[idx] <- paste(rn[idx], "***", sep = "")
    rownames(OV) <- rn
    print(OV)
    lav_msg_stop(gettext("ordered variable(s) has/have only 1 level"))
  }
  # check for mix small/large variances (NOT including exo variables)
  if (!std.ov && !allow.single.case && warn && any(ov$type == "numeric")) {
    num.idx <- which(ov$type == "numeric" & ov$exo == 0L)
    if (length(num.idx) > 0L) {
      min.var <- min(ov$var[num.idx])
      max.var <- max(ov$var[num.idx])
      rel.var <- max.var / min.var
      if (warn && rel.var > 1000) {
        lav_msg_warn(
          gettext("some observed variances are (at least) a factor 1000 times
           larger than others; use varTable(fit) to investigate"))
      }
    }
  }
  # check for really large variances (perhaps -999999 for missing?)
  if (!std.ov && warn && any(ov$type == "numeric")) {
    num.idx <- which(ov$type == "numeric" & ov$exo == 0L)
    if (length(num.idx) > 0L) {
      max.var <- max(ov$var[num.idx])
      if (warn && max.var > 1000000) {
        lav_msg_warn(
          gettext("some observed variances are larger than 1000000
          use varTable(fit) to investigate"))
      }
    }
  }
  # check for all-exogenous variables (eg in f <~ x1 + x2 + x3)
  if (warn && all(ov$exo == 1L)) {
    lav_msg_warn(gettext(
      "all observed variables are exogenous; model may not be identified"))
  }

  # prepare empty lists

  # group-based
  case.idx <- vector("list", length = ngroups)
  Mp <- vector("list", length = ngroups)
  Rp <- vector("list", length = ngroups)
  norig <- vector("list", length = ngroups)
  nobs <- vector("list", length = ngroups)
  X <- vector("list", length = ngroups)
  eXo <- vector("list", length = ngroups)
  Lp <- vector("list", length = ngroups)
  weights <- vector("list", length = ngroups)

  # collect information per upper-level group
  for (g in 1:ngroups) {
    # extract variables in correct order
    if (nlevels > 1L) {
      # keep 'joint' (Y,X) matrix in @X if multilevel (or always?)
      # yes for multilevel (for now); no for clustered only
      OV.NAMES <- unique(c(ov.names[[g]], ov.names.x[[g]]))
      ov.idx <- ov$idx[match(OV.NAMES, ov$name)]
    } else {
      ov.idx <- ov$idx[match(ov.names[[g]], ov$name)]
    }
    exo.idx <- ov$idx[match(ov.names.x[[g]], ov$name)]
    all.idx <- unique(c(ov.idx, exo.idx))

    # extract cases per group
    if (ngroups > 1L || length(group.label) > 0L) {
      if (missing == "listwise") {
        case.idx[[g]] <- which(data[[group]] == group.label[g] &
          complete.cases(data[all.idx]))
        nobs[[g]] <- length(case.idx[[g]])
        norig[[g]] <- length(which(data[[group]] == group.label[g]))
        # } else if(missing == "pairwise" && length(exo.idx) > 0L) {
        #    case.idx[[g]] <- which(data[[group]] == group.label[g] &
        #                           complete.cases(data[exo.idx]))
        #    nobs[[g]] <- length(case.idx[[g]])
        #    norig[[g]] <- length(which(data[[group]] == group.label[g]))
      } else if (length(exo.idx) > 0L && missing != "ml.x") {
        case.idx[[g]] <- which(data[[group]] == group.label[g] &
          complete.cases(data[exo.idx]))
        nobs[[g]] <- length(case.idx[[g]])
        norig[[g]] <- length(which(data[[group]] == group.label[g]))
        if (warn && (nobs[[g]] < norig[[g]])) {
          lav_msg_warn(gettextf(
              "%1$s cases were deleted in group %2$s  due to missing values
              in  exogenous variable(s), while fixed.x = TRUE.",
              (norig[[g]] - nobs[[g]]), group.label[g]))
        }
      } else {
        case.idx[[g]] <- which(data[[group]] == group.label[g])
        nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
      }
    } else {
      if (missing == "listwise") {
        case.idx[[g]] <- which(complete.cases(data[all.idx]))
        nobs[[g]] <- length(case.idx[[g]])
        norig[[g]] <- nrow(data)
        # } else if(missing == "pairwise" && length(exo.idx) > 0L) {
        #    case.idx[[g]] <- which(complete.cases(data[exo.idx]))
        #    nobs[[g]] <- length(case.idx[[g]])
        #    norig[[g]] <- nrow(data)
      } else if (length(exo.idx) > 0L && missing != "ml.x") {
        case.idx[[g]] <- which(complete.cases(data[exo.idx]))
        nobs[[g]] <- length(case.idx[[g]])
        norig[[g]] <- nrow(data)
        if (warn && (nobs[[g]] < norig[[g]])) {
          lav_msg_warn(
            gettextf("%s cases were deleted due to missing values in
                     exogenous variable(s), while fixed.x = TRUE.",
                     (norig[[g]] - nobs[[g]])))
        }
      } else {
        case.idx[[g]] <- 1:nrow(data)
        nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
      }
    }

    # extract data
    X[[g]] <- data.matrix(data[case.idx[[g]], ov.idx, drop = FALSE])
    dimnames(X[[g]]) <- NULL ### copy?

    # sampling weights (but no normalization yet)
    if (!is.null(sampling.weights)) {
      WT <- data[[sampling.weights]][case.idx[[g]]]
      if (any(WT < 0)) {
        lav_msg_stop(gettext("some sampling weights are negative"))
      }

      # check for missing values in sampling weight variable
      if (any(is.na(WT))) {
        lav_msg_stop(gettextf(
          "sampling.weights variable %s contains missing values",
          sQuote(sampling.weights)))
      }

      weights[[g]] <- WT
    }

    # construct integers for user-declared 'ordered' factors
    # FIXME: is this really (always) needed???
    #  (but still better than doing lapply(data[,idx], ordered) which
    #   generated even more copies)
    user.ordered.names <- ov$name[ov$type == "ordered" & ov$user == 1L]
    user.ordered.idx <- which(ov.names[[g]] %in% user.ordered.names)
    if (length(user.ordered.idx) > 0L) {
      for (i in user.ordered.idx) {
        X[[g]][, i][is.na(X[[g]][, i])] <- NA # change NaN to NA
        X[[g]][, i] <- as.numeric(as.factor(X[[g]][, i]))
        # possible alternative to the previous two lines:
        # X[[g]][,i] <- as.numeric(factor(X[[g]][,i], exclude = c(NA, NaN)))
      }
    }

    ## FIXME:
    ## - why also in X? (for samplestats, for now)
    if (length(exo.idx) > 0L) {
      eXo[[g]] <- data.matrix(data[case.idx[[g]], exo.idx, drop = FALSE])
      dimnames(eXo[[g]]) <- NULL
    } else {
      eXo[g] <- list(NULL)
    }

    # standardize observed variables? numeric only!
    if (std.ov) {
      num.idx <- which(ov$name %in% ov.names[[g]] &
        ov$type == "numeric" & ov$exo == 0L)
      if (length(num.idx) > 0L) {
        X[[g]][, num.idx] <-
          scale(X[[g]][, num.idx, drop = FALSE])[, , drop = FALSE]
        # three copies are made!!!!!
      }
      if (length(exo.idx) > 0L) {
        eXo[[g]] <- scale(eXo[[g]])[, , drop = FALSE]
      }
    }

    # response patterns (ordered variables only)
    ord.idx <- which(ov.names[[g]] %in% ov$name[ov$type == "ordered"])
    if (length(ord.idx) > 0L) {
      Rp[[g]] <- lav_data_resp_patterns(X[[g]][, ord.idx, drop = FALSE])
    }

    # warn if we have a small number of observations (but NO error!)
    if (!allow.single.case && warn &&
      nobs[[g]] < (nvar <- length(ov.idx))) {
      txt <- ""
      if (ngroups > 1L) txt <- gettextf("in group %s", g)
      lav_msg_warn(
        gettextf("small number of observations (nobs < nvar) %1$s:
                 nobs = %2$s  nvar = %3$s", txt, nobs[[g]], nvar))
    }
    # check variances per group (if we have multiple groups)
    # to catch zero-variance variables within a group (new in 0.6-8)
    if (ngroups > 1L) {
      # X
      group.var <- apply(X[[g]], 2, var, na.rm = TRUE)
      zero.var <- which(group.var < .Machine$double.eps)
      if (length(zero.var) == 0L) {
        # all is good
      } else {
        # some zero variances!
        gtxt <- if (ngroups > 1L) {
          gettextf("in group %s", g)
        } else {
          ""
        }
       lav_msg_stop(
          gettext("some variables have no variance"), gtxt,
          ":", paste(ov.names[[g]][zero.var], collapse = " "))
      }

      # eXo (if conditional.x = TRUE)...
      if (length(exo.idx) > 0L) {
        group.var <- apply(eXo[[g]], 2, var, na.rm = TRUE)
        zero.var <- which(group.var < .Machine$double.eps)
        if (length(zero.var) == 0L) {
          # all is good
        } else {
          # some zero variances!
          gtxt <- if (ngroups > 1L) {
            gettextf("in group %s", g)
          } else {
            ""
          }
          lav_msg_stop(
            gettext("some exogenous variables have no variance"), gtxt,
            ":", paste(ov.names.x[[g]][zero.var], collapse = " ")
          )
        }
      }
    }

    # cluster information
    if (length(cluster) > 0L) {
      # extract cluster variable(s), for this group
      clus <- data.matrix(data[case.idx[[g]], cluster])
      if (nlevels > 1L) {
        multilevel <- TRUE
      } else {
        multilevel <- FALSE
      }
      # ALWAYS add ov.names.x at the end, even if conditional.x (0.6-7)
      OV.NAMES <- unique(c(ov.names[[g]], ov.names.x[[g]]))
      Lp[[g]] <- lav_data_cluster_patterns(
        Y = X[[g]], clus = clus,
        cluster = cluster,
        multilevel = multilevel,
        ov.names = OV.NAMES,
        ov.names.x = ov.names.x[[g]],
        ov.names.l = ov.names.l[[g]]
      )

      # new in 0.6-4
      # check for 'level-1' variables with zero within variance
      l1.idx <- c(
        Lp[[g]]$within.idx[[2]], # within only
        Lp[[g]]$both.idx[[2]]
      )
      for (v in l1.idx) {
        within.var <- tapply(X[[g]][, v], Lp[[g]]$cluster.idx[[2]],
          FUN = var, na.rm = TRUE
        )
        # ignore singletons
        singleton.idx <- which(Lp[[g]]$cluster.size[[2]] == 1L)
        if (length(singleton.idx) > 0L) {
          within.var[singleton.idx] <- 10 # non-zero variance
        }
        zero.var <- which(within.var < .Machine$double.eps)
        if (length(zero.var) == 0L) {
          # all is good
        } else if (length(zero.var) == length(within.var)) {
          # all zero! possibly a between-level variable
          gtxt <- if (ngroups > 1L) {
            gettextf("in group %s", g)
          } else {
            ""
          }
          lav_msg_warn(
            gettextf("Level-1 variable %1$s has no variance at the within
                     level %2$s. The variable appears to be a between-level
                     variable. Please remove this variable from the level 1
                     section in the model syntax.",
                     dQuote(ov.names[[g]][v])), gtxt)
        } else {
          # some zero variances!
          gtxt <- if (ngroups > 1L) {
            gettextf("in group %s", g)
          } else {
            ""
          }
          lav_msg_warn(gettextf(
          "Level-1 variable %1$s has no variance within some clusters %2$s.
          The cluster ids with zero within variance are: %3$s.",
          dQuote(ov.names[[g]][v])), gtxt,
          lav_msg_view(Lp[[g]]$cluster.id[[2]][zero.var], "none"))
        }
      }

      # new in 0.6-4
      # check for 'level-2' only variables with non-zero within variance
      l2.idx <- Lp[[g]]$between.idx[[2]] # between only
      error.flag <- FALSE
      for (v in l2.idx) {
        within.var <- tapply(X[[g]][, v], Lp[[g]]$cluster.idx[[2]],
          FUN = var, na.rm = TRUE
        )
        non.zero.var <- which(unname(within.var) > .Machine$double.eps)
        if (length(non.zero.var) == 0L) {
          # all is good
        } else if (length(non.zero.var) == 1L) {
          # just one
          gtxt <- if (ngroups > 1L) {
            gettextf("in group %s.", g)
          } else {
            "."
          }
          lav_msg_warn(gettextf(
            "Level-2 variable %1$ss has non-zero variance at the within
            level %2$s in one cluster with id: %3$ss. Please double-check
            if this is a between only variable.",
            dQuote(ov.names[[g]][v])), gtxt,
            Lp[[g]]$cluster.id[[2]][non.zero.var])
        } else {
          error.flag <- TRUE
          # several
          gtxt <- if (ngroups > 1L) {
            gettextf("in group %s", g)
          } else {
            ""
          }
          lav_msg_warn(gettextf(
            "Level-2 variable %1$s has non-zero variance at the within level
            %2$s. The cluster ids with non-zero within variance are: %3$s",
            dQuote(ov.names[[g]][v])), gtxt,
            lav_msg_view(Lp[[g]]$cluster.id[[2]][non.zero.var], "none"))
        }
      }
      if (error.flag) {
        lav_msg_stop(
          gettext("Some between-level (only) variables have
          non-zero variance at the within-level.
          Please double-check your data.")
        )
      }
    } # clustered data

    # missing data
    if (missing != "listwise") {
      if (length(cluster) > 0L) {
        # get missing patterns
        Mp[[g]] <- lav_data_missing_patterns(X[[g]],
          sort.freq = TRUE, coverage = TRUE,
          Lp = Lp[[g]]
        )
      } else {
        # get missing patterns
        Mp[[g]] <- lav_data_missing_patterns(X[[g]],
          sort.freq = TRUE, coverage = TRUE,
          Lp = NULL
        )
      }

      # checking!
      if (length(Mp[[g]]$empty.idx) > 0L) {
        # new in 0.6-4: return 'original' index in full data.frame
        empty.case.idx <- case.idx[[g]][Mp[[g]]$empty.idx]
        if (warn) {
          lav_msg_warn(gettextf(
            "some cases are empty and will be ignored: %s.",
            paste(empty.case.idx, collapse = " ")))
        }
      }
      if (warn && any(Mp[[g]]$coverage == 0)) {
        lav_msg_warn(gettext(
          "due to missing values, some pairwise combinations have 0% coverage;
          use lavInspect(fit, \"coverage\") to investigate."))
      } else if (warn && any(Mp[[g]]$coverage < 0.1)) {
        lav_msg_warn(gettext(
          "due to missing values, some pairwise combinations have less than
          10% coverage; use lavInspect(fit, \"coverage\") to investigate."))
      }
      # in case we had observations with only missings
      nobs[[g]] <- NROW(X[[g]]) - length(Mp[[g]]$empty.idx)
    } # missing
  } # groups, at first level

  # sampling weigths, again
  if (is.null(sampling.weights)) {
    sampling.weights <- character(0L)
  } else {
    # check if we need normalization
    if (sampling.weights.normalization == "none") {
      # nothing to do
    } else if (sampling.weights.normalization == "total") {
      sum.weights <- sum(unlist(weights))
      ntotal <- sum(unlist(nobs))
      for (g in 1:ngroups) {
        WT <- weights[[g]]
        WT2 <- WT / sum.weights * ntotal
        weights[[g]] <- WT2
      }
    } else if (sampling.weights.normalization == "group") {
      for (g in 1:ngroups) {
        WT <- weights[[g]]
        WT2 <- WT / sum(WT) * nobs[[g]]
        weights[[g]] <- WT2
      }
    } else {
      lav_msg_stop(gettext(
        "sampling.weights.normalization should be total, group or none."))
    }
  }

  # block.labels
  block.label <- character(0L)
  if (length(group.label) > 0L && length(level.label) == 0L) {
    block.label <- group.label
  } else if (length(level.label) > 0L && length(group.label) == 0L) {
    block.label <- level.label
  } else if (length(group.label) > 0L &&
    length(level.label) > 0L) {
    block.label <- paste(rep(group.label, each = length(level.label)),
      rep(level.label, times = length(group.label)),
      sep = "."
    )
  }


  lavData <- new("lavData",
    data.type = "full",
    ngroups = ngroups,
    group = group,
    nlevels = nlevels,
    cluster = cluster,
    group.label = group.label,
    level.label = level.label,
    block.label = block.label,
    std.ov = std.ov,
    nobs = nobs,
    norig = norig,
    ov.names = ov.names,
    ov.names.x = ov.names.x,
    ov.names.l = ov.names.l,
    # ov.types        = ov.types,
    # ov.idx          = ov.idx,
    ordered = as.character(ordered),
    weights = weights,
    sampling.weights = sampling.weights,
    ov = ov,
    case.idx = case.idx,
    missing = missing,
    X = X,
    eXo = eXo,
    Mp = Mp,
    Rp = Rp,
    Lp = Lp
  )
  lavData
}
