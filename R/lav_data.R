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
lav_lavdata <- function(data = NULL, # data.frame
                    group = NULL, # multiple groups?
                    cluster = NULL, # clusters?
                    ov_names = NULL, # variables in model
                    ov_names_x = character(0), # exo variables
                    ov_names_l = list(), # names per level
                    ov_names_aux = character(0), # auxiliary variables
                    ordered = NULL, # ordered variables
                    sampling_weights = NULL, # sampling weights
                    sample_cov = NULL, # sample covariance(s)
                    sample_mean = NULL, # sample mean vector(s)
                    sample_th = NULL, # sample thresholds
                    sample_nobs = NULL, # sample nobs
                    lavoptions = lavOptions(), # lavoptions
                    allow_single_case = FALSE # for newdata in predict
) {
  # get info from lavoptions

  # group.labels
  group_label <- lavoptions$group.label
  if (is.null(group_label)) {
    group_label <- character(0L)
  }

  # level.labels
  level_label <- lavoptions$level.label
  if (is.null(level_label)) {
    level_label <- character(0L)
  }

  # allow empty categories of ordinal variable
  allow_empty_cell <- lavoptions$allow.empty.cell

  # block.labels
  block_label <- character(0L)
  if (length(group_label) > 0L && length(level_label) == 0L) {
    block_label <- group_label
  } else if (length(level_label) > 0L && length(group_label) == 0L) {
    block_label <- level_label
  } else if (length(group_label) > 0L &&
    length(level_label) > 0L) {
    block_label <- paste(rep(group_label, each = length(level_label)),
      rep(level_label, times = length(group_label)),
      sep = "."
    )
  }

  # std.ov?
  std_ov <- lavoptions$std.ov
  if (is.null(std_ov)) {
    std_ov <- FALSE
  }
  # partial correlation structure: only these (observed) variables are
  # standardized to unit variance (empty -> all, the usual std.ov behavior)
  correlation_ov <- lavoptions$.correlation.ov
  if (is.null(correlation_ov)) {
    correlation_ov <- character(0L)
  }
  if (length(correlation_ov) > 0L && !is.null(ov_names)) {
    all_ov <- unique(unlist(ov_names))
    bad <- correlation_ov[!(correlation_ov %in% all_ov)]
    if (length(bad) > 0L) {
      lav_msg_warn(gettextf(
        "correlation= argument contains unknown observed variable(s): %s",
        paste(bad, collapse = ", ")))
    }
  }

  # missing? (lav_object_cor() does not parse options before calling lavdata...)
  missing <- tolower(lavoptions$missing)
  if (is.null(missing) || missing == "default") {
    missing <- "listwise"
  } else if (missing %in% c("ml", "fiml", "direct")) {
    missing <- "ml"
  } else if (missing %in% c("ml.x", "fiml.x", "direct.x")) {
    missing <- "ml.x"
  }

  # warn?
  if (allow_single_case) { # eg, in lavPredict
    current_warn <- lav_warn()
    if (lav_warn(FALSE))
        on.exit(lav_warn(current_warn), TRUE)
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
              please use the sample_cov= argument instead"))
          }
        }
        # or perhaps it is a data matrix?
        ### FIXME, we should avoid as.data.frame() and handle
        ### data matrices directly
        data <- as.data.frame(data, stringsAsFactors = FALSE)
      } else {
        lav_msg_stop(gettextf(
          "data= argument is not a data.frame, but of class %s",
          sQuote(class(data))))
      }
    }

    # no ov.names?
    if (is.null(ov_names)) {
      ov_names <- names(data)
      # remove group variable, if provided
      if (length(group) > 0L) {
        group_idx <- which(ov_names == group)
        ov_names <- ov_names[-group_idx]
      }
      # remove cluster variable, if provided
      if (length(cluster) > 0L) {
        cluster_idx <- which(ov_names == cluster)
        ov_names <- ov_names[-cluster_idx]
      }
    }

    lav_data <- lav_data_full(
      data = data,
      group = group,
      cluster = cluster,
      group_label = group_label,
      level_label = level_label,
      block_label = block_label,
      ov_names = ov_names,
      ordered = ordered,
      sampling_weights = sampling_weights,
      sampling_weights_normalization =
        lavoptions$sampling.weights.normalization,
      ov_names_x = ov_names_x,
      ov_names_l = ov_names_l,
      ov_names_aux = ov_names_aux,
      std_ov = std_ov,
      correlation_ov = correlation_ov,
      missing = missing,
      allow_single_case = allow_single_case,
      allow_empty_cell = allow_empty_cell
    )
    sample_cov <- NULL # not needed, but just in case
  }


  # 2) sample moments
  if (is.null(data) && !is.null(sample_cov)) {
    # for now: no levels!!
    nlevels <- 1L

    # we also need the number of observations (per group)
    if (is.null(sample_nobs)) {
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
    if (is.list(sample_cov)) {
      # multiple groups, multiple cov matrices
      if (!is.null(sample_mean)) {
        stopifnot(length(sample_mean) == length(sample_cov))
      }
      if (!is.null(sample_th)) {
        stopifnot(length(sample_th) == length(sample_cov))
      }
      # multiple groups, multiple cov matrices
      ngroups <- length(sample_cov)
      label <- names(sample_cov)
      if (is.null(group_label) || length(group_label) == 0L) {
        if (is.null(label)) {
          group_label <- paste("Group ", 1:ngroups, sep = "")
        } else {
          group_label <- label
        }
      } else {
        if (is.null(label)) {
          stopifnot(length(group_label) == ngroups)
        } else {
          # FIXME!!!!
          # check if they match
        }
      }
    } else {
      ngroups <- 1L
      group_label <- character(0)
      if (!is.matrix(sample_cov)) {
        lav_msg_stop(gettext(
          "sample_cov must be a matrix or a list of matrices"))
      }
      sample_cov <- list(sample_cov)
    }

    # get ov.names
    if (is.null(ov_names)) {
      ov_names <- lapply(sample_cov, row.names)
    } else if (!is.list(ov_names)) {
      # duplicate ov.names for each group
      tmp <- ov_names
      ov_names <- vector("list", length = ngroups)
      ov_names[1:ngroups] <- list(tmp)
    } else {
      if (length(ov_names) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names assumes %1$s groups; data contains %2$s groups",
          length(ov_names), ngroups))
      }
      # nothing to do
    }

    # handle ov.names.x
    if (!is.list(ov_names_x)) {
      tmp <- ov_names_x
      ov_names_x <- vector("list", length = ngroups)
      ov_names_x[1:ngroups] <- list(tmp)
    } else {
      if (length(ov_names_x) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names.x assumes %1$s groups; data contains %2$s groups",
          length(ov_names_x), ngroups))
      }
    }

    ov <- list()
    ov$name <- unique(unlist(c(ov_names, ov_names_x)))
    nvar <- length(ov$name)
    ov$idx <- rep(NA, nvar)
    ov$nobs <- rep(sum(unlist(sample_nobs)), nvar)
    ov$type <- rep("numeric", nvar)
    ov$nlev <- rep(0, nvar)
    # check for categorical
    if (!is.null(sample_th)) {
      th_idx <- attr(sample_th, "th.idx")
      if (is.list(th_idx)) {
        th_idx <- th_idx[[1]] ## FIRST group only (assuming same ths!)
      }
      if (any(th_idx > 0)) {
        tab <- table(th_idx[th_idx > 0])
        ord_idx <- as.numeric(names(tab))
        nlev <- as.integer(unname(tab) + 1)
        ov$type[ord_idx] <- "ordered"
        ov$nlev[ord_idx] <- nlev
      }
    }

    # if std.ov = TRUE, give a warning (suggested by Peter Westfall)
    if (std_ov && !lavoptions$correlation) {
      lav_msg_warn(gettext(
        "std.ov argument is ignored if only sample statistics are provided."))
    }

    # check variances (new in 0.6-7)
    if (!allow_single_case) {
      for (g in seq_len(ngroups)) {
        var_1 <- diag(sample_cov[[g]])
        # 1. finite?
        if (!all(is.finite(var_1))) {
          lav_msg_stop(gettext(
          "at least one variance in the sample
          covariance matrix is not finite."))
          }
        # 2. near zero (or negative)?
        if (any(var_1 < .Machine$double.eps)) {
          lav_msg_stop(
            gettext("at least one variance in the sample covariance matrix is
            (near) zero or negative."))
        }
        # 3. very large?
        max_var <- max(var_1)
        if (max_var > 1000000) {
          lav_msg_warn(
            gettext("some observed variances in the sample covariance matrix
            are larger than 1000000."))
        }
      }
    }

    # block.labels
    block_label <- character(0L)
    if (length(group_label) > 0L && length(level_label) == 0L) {
      block_label <- group_label
    } else if (length(level_label) > 0L && length(group_label) == 0L) {
      block_label <- level_label
    } else if (length(group_label) > 0L &&
      length(level_label) > 0L) {
      block_label <- paste(rep(group_label, each = length(level_label)),
        rep(level_label, times = length(group_label)),
        sep = "."
      )
    }

    # construct lavData object
    lav_data <- new("lavData",
      data.type = "moment",
      ngroups = ngroups,
      group = group,
      nlevels = 1L, # for now
      cluster = character(0L),
      group.label = group_label,
      level.label = character(0L),
      block.label = block_label,
      nobs = as.list(sample_nobs),
      norig = as.list(sample_nobs),
      ov.names = ov_names,
      ov.names.x = ov_names_x,
      ov.names.l = ov_names_l,
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

  # 3) data.type = "none":  both data and sample_cov are NULL
  if (is.null(data) && is.null(sample_cov)) {
    # clustered/multilevel? --> ov.names.l should be filled in
    if (length(ov_names_l) > 0L) {
      nlevels <- length(ov_names_l[[1]]) # we assume the same number
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
      if (length(level_label) == 0L) {
        level_label <- c("within", cluster)
      } else {
        # check if length(level.label) = 1 + length(cluster)
        if (length(level_label) != length(cluster) + 1L) {
          lav_msg_stop(gettext("length(level.label) != length(cluster) + 1L"))
        }
        # nothing to do
      }
    } else {
      nlevels <- 1L
      cluster <- character(0L)
      level_label <- character(0L)
    }

    # ngroups: ov.names (when group: is used), or sample_nobs
    if (is.null(ov_names)) {
      lav_msg_warn(gettext("ov.names is NULL"))
      ov_names <- character(0L)
      if (is.null(sample_nobs)) {
        ngroups <- 1L
        sample_nobs <- rep(list(0L), ngroups)
      } else {
        sample_nobs <- as.list(sample_nobs)
        ngroups <- length(sample_nobs)
      }
    } else if (!is.list(ov_names)) {
      if (is.null(sample_nobs)) {
        ngroups <- 1L
        sample_nobs <- rep(list(0L), ngroups)
      } else {
        sample_nobs <- as.list(sample_nobs)
        ngroups <- length(sample_nobs)
      }
      ov_names <- rep(list(ov_names), ngroups)
    } else if (is.list(ov_names)) {
      ngroups <- length(ov_names)
      if (is.null(sample_nobs)) {
        sample_nobs <- rep(list(0L), ngroups)
      } else {
        sample_nobs <- as.list(sample_nobs)
        if (length(sample_nobs) != ngroups) {
          lav_msg_stop(gettextf(
            "length(sample_nobs) = %1$s but syntax implies ngroups = %2$s",
            length(sample_nobs), ngroups))
        }
      }
    }


    # group.label
    if (ngroups > 1L) {
      if (is.null(group)) {
        group <- "group"
      }
      group_label <- paste("Group", 1:ngroups, sep = "")
    } else {
      group <- character(0L)
      group_label <- character(0L)
    }

    # handle ov.names.x
    if (!is.list(ov_names_x)) {
      ov_names_x <- rep(list(ov_names_x), ngroups)
    }

    # handle ov.names.l: we need a model for each group
    if (length(ov_names_l) == 1L && ngroups > 1L) {
      lav_msg_stop(gettextf(
        "this analysis assumes %s groups, but the multilevel model syntax is
        specified for a single group only; combining multiple groups with
        multiple levels requires explicit group: blocks in the model
        syntax, one for each group (for now).", ngroups))
    }

    ov <- list()
    ov$name <- unique(unlist(c(ov_names, ov_names_x)))
    nvar <- length(ov$name)
    ov$idx <- rep(NA, nvar)
    ov$nobs <- rep(0L, nvar)
    ov$type <- rep("numeric", nvar)
    ov$nlev <- rep(0L, nvar)

    # collect information per upper-level group
    lp <- vector("list", length = ngroups)
    for (g in 1:ngroups) {
      if (nlevels > 1L) {
        # ALWAYS add ov.names.x at the end, even if conditional.x
        ov_names_1 <- unique(c(ov_names[[g]], ov_names_x[[g]]))
        lp[[g]] <- lav_data_cl_patterns(
          y = NULL, clus = NULL,
          cluster = cluster,
          multilevel = TRUE,
          ov_names = ov_names_1,
          ov_names_x = ov_names_x[[g]],
          ov_names_l = ov_names_l[[g]]
        )
      }
    } # g

    # block.labels
    block_label <- character(0L)
    if (length(group_label) > 0L && length(level_label) == 0L) {
      block_label <- group_label
    } else if (length(level_label) > 0L && length(group_label) == 0L) {
      block_label <- level_label
    } else if (length(group_label) > 0L &&
      length(level_label) > 0L) {
      block_label <- paste(rep(group_label, each = length(level_label)),
        rep(level_label, times = length(group_label)),
        sep = "."
      )
    }

    # construct lavData object
    lav_data <- new("lavData",
      data.type = "none",
      ngroups = ngroups,
      group = group,
      nlevels = nlevels,
      cluster = cluster,
      group.label = group_label,
      level.label = level_label,
      block.label = block_label,
      nobs = sample_nobs,
      norig = sample_nobs,
      ov.names = ov_names,
      ov.names.x = ov_names_x,
      ov.names.l = ov_names_l,
      ordered = as.character(ordered),
      weights = vector("list", length = ngroups),
      sampling.weights = character(0L),
      ov = ov,
      missing = "listwise",
      case.idx = vector("list", length = ngroups),
      Mp = vector("list", length = ngroups),
      Rp = vector("list", length = ngroups),
      Lp = lp,
      X = vector("list", length = ngroups),
      eXo = vector("list", length = ngroups)
    )
  }

  lav_data
}


# validate the aux= argument and return the names of the auxiliary variables
# that can actually be used in the current configuration
#
# auxiliary variables are observed variables that are NOT part of the model;
# they are never part of the model-implied summary statistics. For now, they
# are only used to improve the EM-based h1 (saturated) moments when
# missing = "ml" (continuous data). In any other configuration, they are
# ignored (with a note).
#
# binary or ordered/categorical auxiliary variables are not supported and are
# removed (with a warning).
lav_data_aux_check <- function(aux        = NULL,
                               data       = NULL,
                               ov_names   = NULL,
                               ordered    = NULL,
                               lavoptions = NULL) {
  # nothing to do?
  if (is.null(aux)) {
    return(character(0L))
  }
  if (!is.character(aux)) {
    lav_msg_stop(gettext("aux= argument must be a character vector."))
  }
  aux <- unique(aux[nchar(aux) > 0L])
  if (length(aux) == 0L) {
    return(character(0L))
  }

  # we can only use auxiliary variables if we have raw data, the model is
  # continuous, and missing is one of the EM-based options
  missing <- tolower(lavoptions$missing)
  aux_missing <- c("ml", "ml.x", "two.stage", "robust.two.stage")
  if (is.null(data) || !(is.data.frame(data) || is.matrix(data)) ||
      isTRUE(lavoptions$categorical) ||
      !missing %in% aux_missing) {
    lav_msg_note(gettext(
      "auxiliary (aux) variables are currently only used (to improve the
       EM-based h1/saturated summary statistics) when missing = \"ml\",
       \"two.stage\" or \"robust.two.stage\" with continuous data;
       they will be ignored here."))
    return(character(0L))
  }

  # drop aux names that are not usable (not in data, model variables, or
  # binary/ordered-categorical); shared with the FIML saturated-correlates path
  aux <- lav_aux_clean_names(
    aux = aux, data = data, ov_names = ov_names, ordered = ordered
  )
  if (length(aux) == 0L) {
    return(character(0L))
  }

  # note: for the two-stage (and robust two-stage) standard errors, the
  # auxiliary variables are accounted for via the augmented stage-1 ACOV
  # (Savalei & Bentler, 2009; see lav_model_nvcov_two_stage()), except when
  # fixed-x covariates are present, in which case the standard errors fall
  # back to the model-only ACOV

  aux
}


# handle full data
lav_data_full <- function(data = NULL, # data.frame
                          group = NULL, # multiple groups?
                          cluster = NULL, # clustered?
                          group_label = NULL, # custom group labels?
                          level_label = NULL,
                          block_label = NULL,
                          ov_names = NULL, # variables needed
                          # in model
                          ordered = NULL, # ordered variables
                          sampling_weights = NULL, # sampling weights
                          sampling_weights_normalization = "none",
                          ov_names_x = character(0L), # exo variables
                          ov_names_l = list(), # var per level
                          ov_names_aux = character(0L), # auxiliary variables
                          std_ov = FALSE, # standardize ov's?
                          correlation_ov = character(0L), # std only these
                          missing = "listwise", # remove missings?
                          allow_single_case = FALSE, # allow single case?
                          allow_empty_cell = FALSE
) {
  # number of groups and group labels
  if (!is.null(group) && length(group) > 0L) {
    if (!(group %in% names(data))) {
      lav_msg_stop(gettextf(
        "grouping variable %1$s not found; variable names
        found in data frame are: %2$s",
        sQuote(group), paste(names(data), collapse = " ")))
    }
    # note: by default, we use the order as in the data;
    # not as in levels(data[,group])
    if (length(group_label) == 0L) {
      group_label <- unique(as.character(data[[group]]))
      if (any(is.na(group_label))) {
        lav_msg_warn(gettextf("group variable %s contains missing values",
                              sQuote(group)))
      }
      group_label <- group_label[!is.na(group_label)]
    } else {
      group_label <- unique(as.character(group_label))
      # check if user-provided group labels exist
      label <- unique(as.character(data[[group]]))
      idx <- match(group_label, label)
      if (any(is.na(idx))) {
        lav_msg_warn(gettextf(
          "some group.labels do not appear in the grouping variable: %s",
          lav_msg_view(group_label[which(is.na(idx))], log_sep = "none"))
        )
      }
      group_label <- group_label[!is.na(idx)]
      # any groups left?
      if (length(group_label) == 0L) {
        lav_msg_stop(gettext(
          "no group levels left; check the group.label argument"))
      }
    }
    ngroups <- length(group_label)
  } else {
    if (length(group_label) > 0L) {
      lav_msg_warn(gettext(
       "`group.label' argument will be ignored if `group' argument is missing"))
    }
    ngroups <- 1L
    group_label <- character(0L)
    group <- character(0L)
  }

  # ensure allow.empty.cell is logical
  if (is.null(allow_empty_cell)) allow_empty_cell <- FALSE

  # sampling weights
  if (!is.null(sampling_weights)) {
    if (is.character(sampling_weights)) {
      if (!(sampling_weights %in% names(data))) {
        lav_msg_stop(
          gettextf("sampling weights variable %1$s not found;
                   variable names found in data frame are: %2$s",
          sQuote(sampling_weights), paste(names(data), collapse = " ")))
      }
      # check for missing values in sampling weight variable
      if (any(is.na(data[[sampling_weights]]))) {
        lav_msg_stop(
          gettextf("sampling.weights variable %s contains missing values",
          sQuote(sampling_weights)))
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
      not_ok <- which(!cluster %in% names(data))
      lav_msg_stop(gettextf(
        "cluster variable(s) %1$s not found;
        variable names found in data frame are: %2$s",
        sQuote(cluster[not_ok]), paste(names(data), collapse = " ")))
    }

    # check for missing values in cluster variable(s)
    for (cl in seq_along(cluster)) {
      if (anyNA(data[[cluster[cl]]])) {
        lav_msg_warn(gettextf("cluster variable %s contains missing values",
          sQuote(cluster[cl])))
      }
    }

    # multilevel?
    if (length(ov_names_l) > 0L) {
      # default level.labels
      if (length(level_label) == 0L) {
        level_label <- c("within", cluster)
      } else {
        # check if length(level.label) = 1 + length(cluster)
        if (length(level_label) != length(cluster) + 1L) {
          lav_msg_stop(gettext("length(level.label) != length(cluster) + 1L"))
        }
        # nothing to do
      }
      nlevels <- length(level_label)
    } else {
      # just clustered data, but no random effects
      nlevels <- 1L
      level_label <- character(0L)
    }
  } else {
    if (length(level_label) > 0L) {
      lav_msg_warn(gettext(
       "`level.label' argument will be ignored if `cluster' argument is missing"
      ))
    }
    nlevels <- 1L
    level_label <- character(0L)
    cluster <- character(0L)
  }

  # check ov.names vs ngroups
  if (ngroups > 1L) {
    if (is.list(ov_names)) {
      if (length(ov_names) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names assumes %1$s groups; data contains %2$s groups",
          length(ov_names), ngroups))
      }
    } else {
      tmp <- ov_names
      ov_names <- vector("list", length = ngroups)
      ov_names[1:ngroups] <- list(tmp)
    }
    if (is.list(ov_names_x)) {
      if (length(ov_names_x) != ngroups) {
        lav_msg_stop(gettextf(
          "ov.names.x assumes %1$s groups; data contains %2$s groups",
          length(ov_names_x), ngroups))
      }
    } else {
      tmp <- ov_names_x
      ov_names_x <- vector("list", length = ngroups)
      ov_names_x[1:ngroups] <- list(tmp)
    }
    # multilevel: we need a model (with level: blocks) for each group
    if (length(ov_names_l) > 0L && length(ov_names_l) != ngroups) {
      if (length(ov_names_l) == 1L) {
        lav_msg_stop(gettextf(
          "the data contains %s groups, but the multilevel model syntax is
          specified for a single group only; combining multiple groups with
          multiple levels requires explicit group: blocks in the model
          syntax, one for each group (for now).", ngroups))
      } else {
        lav_msg_stop(gettextf(
          "ov.names.l assumes %1$s groups; data contains %2$s groups",
          length(ov_names_l), ngroups))
      }
    }
  } else {
    if (is.list(ov_names)) {
      if (length(ov_names) > 1L) {
        lav_msg_stop(gettext(
          "model syntax defines multiple groups; data suggests a single group"))
      }
    } else {
      ov_names <- list(ov_names)
    }
    if (is.list(ov_names_x)) {
      if (length(ov_names_x) > 1L) {
        lav_msg_stop(gettext(
          "model syntax defines multiple groups; data suggests a single group"))
      }
    } else {
      ov_names_x <- list(ov_names_x)
    }
  }

  # check if all ov.names can be found in the data.frame
  for (g in 1:ngroups) {
    # does the data contain all the observed variables
    # needed in the user-specified model for this group
    ov_all <- unique(c(ov_names[[g]], ov_names_x[[g]])) # no overlap if categ

    # handle interactions
    ov_int_names <- ov_all[grepl(":", ov_all)]
    n_int <- length(ov_int_names)
    if (n_int > 0L) {
      ov_names_noint <- ov_all[!ov_all %in% ov_int_names]
      for (iv in seq_len(n_int)) {
        names_1 <- strsplit(ov_int_names[iv], ":", fixed = TRUE)[[1L]]
        if (all(names_1 %in% ov_names_noint)) {
          # add this interaction term to the data.frame, unless
          # it already exists
          if (is.null(data[[ov_int_names[iv]]])) {
            data[[ov_int_names[iv]]] <-
              data[[names_1[1L]]] * data[[names_1[2L]]]
          }
        }
      }
    }

    # check for missing observed variables
    idx_missing <- which(!(ov_all %in% names(data)))

    if (length(idx_missing)) {
      lav_msg_stop(
        gettextf("some (observed) variables specified in the model are
        not found in the dataset: %s",
        paste(ov_all[idx_missing], collapse = " ")))
    }
  }


  # here, we know for sure all ov.names exist in the data.frame
  # create varTable
  # FIXME: should we add the 'group'/'cluster' variable (no for now)
  ov <- lav_dataframe_vartable(
    frame = data, ov_names = ov_names,
    ov_names_x = ov_names_x, ordered = ordered,
    as_data_frame = FALSE, allow_empty_cell = allow_empty_cell
  )

  # variable table for the auxiliary/instrument variables (their type and
  # number of levels are needed, eg by the categorical IV estimator, to compute
  # the augmented polychoric sample statistics)
  ov_aux <- list()
  if (length(ov_names_aux) > 0L) {
    ov_aux <- lav_dataframe_vartable(
      frame = data, ov_names = ov_names_aux,
      ordered = ordered, as_data_frame = FALSE,
      allow_empty_cell = allow_empty_cell
    )
  }

  # do some checking
  # check for unordered factors (but only if nlev > 2)
  if ("factor" %in% ov$type) {
    f_names <- ov$name[ov$type == "factor" & ov$nlev > 2L]
    f_names_all <- ov$name[ov$type == "factor"]
    ov_names_1 <- unlist(ov_names)
    ov_names_x_1 <- unlist(ov_names_x)
    ov_names_nox <- ov_names_1[!ov_names_1 %in% ov_names_x_1]
    if (any(f_names %in% ov_names_x_1)) {
      lav_msg_stop(
        gettext("unordered factor(s) with more than 2 levels detected
            as exogenous covariate(s): "),
        paste(f_names, collapse = " "))
    } else if (any(f_names_all %in% ov_names_nox)) {
      lav_msg_stop(
        gettext("unordered factor(s) detected; make them numeric or ordered:"),
        paste(f_names_all, collapse = " "))
    }
  }
  # check for ordered exogenous variables
  if ("ordered" %in% ov$type[ov$name %in% unlist(ov_names_x)]) {
    f_names <- ov$name[ov$type == "ordered" &
      ov$name %in% unlist(ov_names_x)]
    if (any(f_names %in% unlist(ov_names_x))) {
      lav_msg_warn(gettextf(
        "exogenous variable(s) declared as ordered in data: %s",
        lav_msg_view(f_names, log_sep = "none")))
    }
  }
  # check for ordered endogenous variables with more than 12 levels
  if ("ordered" %in% ov$type[!ov$name %in% unlist(ov_names_x)]) {
    f_names <- ov$name[ov$type == "ordered" &
      !ov$name %in% unlist(ov_names_x) &
      ov$nlev > 12L]
    if (length(f_names) > 0L) {
      lav_msg_warn(gettextf(
        "some ordered categorical variable(s) have more than 12 levels: %s",
        lav_msg_view(f_names, log_sep = "none")))
    }
  }
  # check for zero-cases
  idx <- which(ov$nobs == 0L | ov$var == 0)
  if (!allow_single_case && length(idx) > 0L) {
    ov_1 <- as.data.frame(ov)
    rn <- rownames(ov_1)
    rn[idx] <- paste(rn[idx], "***", sep = "")
    rownames(ov_1) <- rn
    print(ov_1)
    lav_msg_stop(gettext(
      "some variables have no values (only missings) or no variance"))
  }
  # check for single cases (no variance!)
  idx <- which(ov$nobs == 1L | (ov$type == "numeric" & !is.finite(ov$var)))
  if (!allow_single_case && length(idx) > 0L) {
    ov_1 <- as.data.frame(ov)
    rn <- rownames(ov_1)
    rn[idx] <- paste(rn[idx], "***", sep = "")
    rownames(ov_1) <- rn
    print(ov_1)
    lav_msg_stop(gettext(
      "some variables have only 1 observation or no finite variance"))
  }
  # check for ordered variables with only 1 level
  idx <- which(ov$type == "ordered" & ov$nlev == 1L)
  if (!allow_single_case && length(idx) > 0L) {
    ov_1 <- as.data.frame(ov)
    rn <- rownames(ov_1)
    rn[idx] <- paste(rn[idx], "***", sep = "")
    rownames(ov_1) <- rn
    print(ov_1)
    lav_msg_stop(gettext("ordered variable(s) has/have only 1 level"))
  }
  # check for mix small/large variances (NOT including exo variables)
  if (!std_ov && !allow_single_case && any(ov$type == "numeric")) {
    num_idx <- which(ov$type == "numeric" & ov$exo == 0L)
    if (length(num_idx) > 0L) {
      min_var <- min(ov$var[num_idx])
      max_var <- max(ov$var[num_idx])
      rel_var <- max_var / min_var
      if (rel_var > 1000) {
        lav_msg_warn(
          gettext("some observed variances are (at least) a factor 1000 times
           larger than others; use varTable(fit) to investigate"))
      }
    }
  }
  # check for really large variances (perhaps -999999 for missing?)
  if (!allow_single_case && !std_ov && any(ov$type == "numeric")) {
    num_idx <- which(ov$type == "numeric" & ov$exo == 0L)
    if (length(num_idx) > 0L) {
      max_var <- max(ov$var[num_idx])
      if (max_var > 1000000) {
        lav_msg_warn(
          gettext("some observed variances are larger than 1000000
          use varTable(fit) to investigate"))
      }
    }
  }
  # check for all-exogenous variables (eg in f <~ x1 + x2 + x3)
  if (all(ov$exo == 1L)) {
    lav_msg_warn(gettext(
      "all observed variables are exogenous; model may not be identified"))
  }
  # check for perfect correlations (NOT including exo variables)
  if (!allow_single_case && any(ov$type == "numeric")) {
    num_idx <- which(ov$type == "numeric" & ov$exo == 0L)
    cor_1 <- try(cor(data[, ov$idx[num_idx]], use = "pairwise.complete.obs"),
               silent = TRUE)
    # replace any NAs by 0 (as we only wish to detect perfect correlations)
    cor_1[is.na(cor_1)] <- 0
    if (!inherits(cor_1, "try-error") &&
        any(lav_mat_vech(cor_1, diagonal = FALSE) == 1)) {
      cor_1[upper.tri(cor_1, diag = TRUE)] <- 0
      idx <- which(cor_1 == 1)
      this_names <- ov$name[num_idx]
      bad_names <- this_names[sort(unique(c(row(cor_1)[idx],
                              col(cor_1)[idx]))) ]
      # should we make this a hard stop? things will most likely fail later...
      lav_msg_warn(gettextf(
      "some observed variables are perfectly correlated;
       please check your data; variables involved are: %s",
       paste(bad_names, collapse = " ")))
    }
  }

  # prepare empty lists

  # group-based
  case_idx <- vector("list", length = ngroups)
  mp <- vector("list", length = ngroups)
  rp <- vector("list", length = ngroups)
  norig <- vector("list", length = ngroups)
  nobs <- vector("list", length = ngroups)
  x <- vector("list", length = ngroups)
  exo <- vector("list", length = ngroups)
  aux <- vector("list", length = ngroups)
  lp <- vector("list", length = ngroups)
  weights <- vector("list", length = ngroups)

  # collect information per upper-level group
  # datam <- data.matrix(data) # YR: not yet, this breaks the stuart package!
  for (g in 1:ngroups) {
    # extract variables in correct order
    if (nlevels > 1L) {
      # keep 'joint' (Y,X) matrix in @X if multilevel (or always?)
      # yes for multilevel (for now); no for clustered only
      ov_names_2 <- unique(c(ov_names[[g]], ov_names_x[[g]]))
      ov_idx <- ov$idx[match(ov_names_2, ov$name)]
    } else {
      ov_idx <- ov$idx[match(ov_names[[g]], ov$name)]
    }
    exo_idx <- ov$idx[match(ov_names_x[[g]], ov$name)]
    all_idx <- unique(c(ov_idx, exo_idx))

    # extract cases per group
    if (ngroups > 1L || length(group_label) > 0L) {
      if (missing == "listwise") {
        case_idx[[g]] <- which(data[[group]] == group_label[g] &
          complete.cases(data[all_idx]))
        nobs[[g]] <- length(case_idx[[g]])
        norig[[g]] <- length(which(data[[group]] == group_label[g]))
        # check for empty data
        if (nobs[[g]] == 0L) {
          lav_msg_stop(gettextf(
            "all observations were deleted due to missing
             data after listwise deletion in group [%s];
             check your data or consider a different option
             for the missing= argument.", group_label[g]))
        }
        # } else if(missing == "pairwise" && length(exo.idx) > 0L) {
        #    case.idx[[g]] <- which(data[[group]] == group.label[g] &
        #                           complete.cases(data[exo.idx]))
        #    nobs[[g]] <- length(case.idx[[g]])
        #    norig[[g]] <- length(which(data[[group]] == group.label[g]))
      } else if (length(exo_idx) > 0L && missing != "ml.x") {
        case_idx[[g]] <- which(data[[group]] == group_label[g] &
          complete.cases(data[exo_idx]))
        nobs[[g]] <- length(case_idx[[g]])
        norig[[g]] <- length(which(data[[group]] == group_label[g]))
        if ((nobs[[g]] < norig[[g]])) {
          lav_msg_warn(gettextf(
              "%1$s cases were deleted in group %2$s  due to missing values
              in  exogenous variable(s), while fixed.x = TRUE.",
              (norig[[g]] - nobs[[g]]), group_label[g]))
        }
      } else {
        case_idx[[g]] <- which(data[[group]] == group_label[g])
        nobs[[g]] <- norig[[g]] <- length(case_idx[[g]])
      }
    } else {
      if (missing == "listwise") {
        case_idx[[g]] <- which(complete.cases(data[all_idx]))
        nobs[[g]] <- length(case_idx[[g]])
        norig[[g]] <- nrow(data)
        if (nobs[[g]] == 0L) {
          lav_msg_stop(gettext(
            "all observations were deleted due to missing
             data after listwise deletion;
             check your data or consider a different option
             for the missing= argument."))
        }
        # } else if(missing == "pairwise" && length(exo.idx) > 0L) {
        #    case.idx[[g]] <- which(complete.cases(data[exo.idx]))
        #    nobs[[g]] <- length(case.idx[[g]])
        #    norig[[g]] <- nrow(data)
      } else if (length(exo_idx) > 0L && missing != "ml.x") {
        case_idx[[g]] <- which(complete.cases(data[exo_idx]))
        nobs[[g]] <- length(case_idx[[g]])
        norig[[g]] <- nrow(data)
        if ((nobs[[g]] < norig[[g]])) {
          lav_msg_warn(
            gettextf("%s cases were deleted due to missing values in
                     exogenous variable(s), while fixed.x = TRUE.",
                     (norig[[g]] - nobs[[g]])))
        }
      } else {
        case_idx[[g]] <- seq_len(nrow(data))
        nobs[[g]] <- norig[[g]] <- length(case_idx[[g]])
      }
    }

    # extract data
    #X[[g]] <- datam[case.idx[[g]], ov.idx, drop = FALSE]
    x[[g]] <- data.matrix(data[case_idx[[g]], ov_idx, drop = FALSE])
    dimnames(x[[g]]) <- NULL ### copy?

    # sampling weights (but no normalization yet)
    if (!is.null(sampling_weights)) {
      wt <- data[[sampling_weights]][case_idx[[g]]]
      if (any(wt < 0)) {
        lav_msg_stop(gettext("some sampling weights are negative"))
      }

      # check for missing values in sampling weight variable
      if (any(is.na(wt))) {
        lav_msg_stop(gettextf(
          "sampling.weights variable %s contains missing values",
          sQuote(sampling_weights)))
      }

      weights[[g]] <- wt
    }

    # construct integers for user-declared 'ordered' factors
    # FIXME: is this really (always) needed???
    #  (but still better than doing lapply(data[,idx], ordered) which
    #   generated even more copies)
    user_ordered_names <- ov$name[ov$type == "ordered" & ov$user == 1L]
    user_ordered_idx <- which(ov_names[[g]] %in% user_ordered_names)
    if (length(user_ordered_idx) > 0L) {
      for (i in user_ordered_idx) {
        x[[g]][, i][is.na(x[[g]][, i])] <- NA # change NaN to NA
        if (!allow_empty_cell) x[[g]][, i] <- as.numeric(as.factor(x[[g]][, i]))
        # possible alternative to the previous two lines:
        # X[[g]][,i] <- as.numeric(factor(X[[g]][,i], exclude = c(NA, NaN)))
      }
    }

    ## FIXME:
    ## - why also in X? (for samplestats, for now)
    if (length(exo_idx) > 0L) {
      exo[[g]] <- data.matrix(data[case_idx[[g]], exo_idx, drop = FALSE])
      dimnames(exo[[g]]) <- NULL
    } else {
      exo[g] <- list(NULL)
    }

    # auxiliary variables (not part of the model); kept here only so that
    # downstream code (eg the EM-based h1 moments) can use them; they never
    # enter the model-implied summary statistics. Rows are aligned with
    # case.idx[[g]], hence with x[[g]].
    if (length(ov_names_aux) > 0L) {
      aux[[g]] <- data.matrix(data[case_idx[[g]], ov_names_aux, drop = FALSE])
      dimnames(aux[[g]]) <- NULL
    } else {
      aux[g] <- list(NULL)
    }

    # standardize observed variables? numeric only!
    if (std_ov) {
      num_idx <- which(ov$name %in% ov_names[[g]] &
        ov$type == "numeric" & ov$exo == 0L)
      # partial correlation structure: only standardize the requested
      # variables, leave the others in their original metric
      if (length(correlation_ov) > 0L && length(num_idx) > 0L) {
        num_idx <- num_idx[ov$name[num_idx] %in% correlation_ov]
      }
      if (length(num_idx) > 0L) {
        x[[g]][, num_idx] <-
          scale(x[[g]][, num_idx, drop = FALSE])[, , drop = FALSE]
        # three copies are made!!!!!
      }
      if (length(exo_idx) > 0L) {
        exo[[g]] <- scale(exo[[g]])[, , drop = FALSE]
      }
    }

    # response patterns (ordered variables only)
    ord_idx <- which(ov_names[[g]] %in% ov$name[ov$type == "ordered"])
    if (length(ord_idx) > 0L) {
      rp[[g]] <- lav_data_resp_patterns(x[[g]][, ord_idx, drop = FALSE])
    }

    # warn if we have a small number of observations (but NO error!)
    if (!allow_single_case &&
      nobs[[g]] < (nvar <- length(ov_idx))) {
      txt <- ""
      if (ngroups > 1L) txt <- gettextf("in group %s", g)
      lav_msg_warn(
        gettextf("small number of observations (nobs < nvar) %1$s:
                 nobs = %2$s  nvar = %3$s", txt, nobs[[g]], nvar))
    }
    # check variances per group (if we have multiple groups)
    # to catch zero-variance variables within a group (new in 0.6-8)
    if (ngroups > 1L && !allow_empty_cell) {
      # X
      group_var <- apply(x[[g]], 2, var, na.rm = TRUE)
      zero_var <- which(group_var < .Machine$double.eps)
      if (length(zero_var) == 0L) {
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
          ":", paste(ov_names[[g]][zero_var], collapse = " "))
      }

      # eXo (if conditional.x = TRUE)...
      if (length(exo_idx) > 0L) {
        group_var <- apply(exo[[g]], 2, var, na.rm = TRUE)
        zero_var <- which(group_var < .Machine$double.eps)
        if (length(zero_var) == 0L) {
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
            ":", paste(ov_names_x[[g]][zero_var], collapse = " ")
          )
        }
      }
    }

    # cluster information
    if (length(cluster) > 0L) {
      # extract cluster variable(s), for this group
      clus <- data.matrix(data[case_idx[[g]], cluster])
      if (nlevels > 1L) {
        multilevel <- TRUE
      } else {
        multilevel <- FALSE
      }
      # ALWAYS add ov.names.x at the end, even if conditional.x (0.6-7)
      ov_names_2 <- unique(c(ov_names[[g]], ov_names_x[[g]]))
      lp[[g]] <- lav_data_cl_patterns(
        y = x[[g]], clus = clus,
        cluster = cluster,
        multilevel = multilevel,
        ov_names = ov_names_2,
        ov_names_x = ov_names_x[[g]],
        ov_names_l = ov_names_l[[g]]
      )

      # new in 0.6-4
      # check for 'level-1' variables with zero within variance
      l1_idx <- c(
        lp[[g]]$within.idx[[2]], # within only
        lp[[g]]$both.idx[[2]]
      )
      for (v in l1_idx) {
        within_var <- tapply(x[[g]][, v], lp[[g]]$cluster.idx[[2]],
          FUN = var, na.rm = TRUE
        )
        # ignore singletons
        singleton_idx <- which(lp[[g]]$cluster.size[[2]] == 1L)
        if (length(singleton_idx) > 0L) {
          within_var[singleton_idx] <- 10 # non-zero variance
        }
        zero_var <- which(within_var < .Machine$double.eps)
        if (length(zero_var) == 0L) {
          # all is good
        } else if (length(zero_var) == length(within_var)) {
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
                     dQuote(ov_names[[g]][v]), gtxt))
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
          dQuote(ov_names[[g]][v]), gtxt,
          lav_msg_view(lp[[g]]$cluster.id[[2]][zero_var], "none")))
        }
      }

      # new in 0.6-4
      # check for 'level-2' only variables with non-zero within variance
      l2_idx <- lp[[g]]$between.idx[[2]] # between only
      error_flag <- FALSE
      for (v in l2_idx) {
        within_var <- tapply(x[[g]][, v], lp[[g]]$cluster.idx[[2]],
          FUN = var, na.rm = TRUE
        )
        non_zero_var <- which(unname(within_var) > .Machine$double.eps)
        if (length(non_zero_var) == 0L) {
          # all is good
        } else if (length(non_zero_var) == 1L) {
          # just one
          gtxt <- if (ngroups > 1L) {
            gettextf("in group %s.", g)
          } else {
            "."
          }
          lav_msg_warn(gettextf(
            "Level-2 variable %1$s has non-zero variance at the within
            level %2$s in one cluster with id: %3$s. Please double-check
            if this is a between only variable.",
            dQuote(ov_names[[g]][v]), gtxt,
            lp[[g]]$cluster.id[[2]][non_zero_var]))
        } else {
          error_flag <- TRUE
          # several
          gtxt <- if (ngroups > 1L) {
            gettextf("in group %s", g)
          } else {
            ""
          }
          lav_msg_warn(gettextf(
            "Level-2 variable %1$s has non-zero variance at the within level
            %2$s. The cluster ids with non-zero within variance are: %3$s",
            dQuote(ov_names[[g]][v]), gtxt,
            lav_msg_view(lp[[g]]$cluster.id[[2]][non_zero_var], "none")))
        }
      }
      if (error_flag) {
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
        mp[[g]] <- lav_data_mi_patterns(x[[g]],
          sort_freq = TRUE, coverage = TRUE,
          lp = lp[[g]]
        )
      } else {
        # get missing patterns
        mp[[g]] <- lav_data_mi_patterns(x[[g]],
          sort_freq = TRUE, coverage = TRUE,
          lp = NULL
        )
      }

      # checking!
      if (length(mp[[g]]$empty.idx) > 0L) {
        # new in 0.6-4: return 'original' index in full data.frame
        empty_case_idx <- case_idx[[g]][mp[[g]]$empty.idx]
        lav_msg_warn(gettextf(
          "some cases are empty and will be ignored: %s.",
          paste(empty_case_idx, collapse = " ")))
      }
    if (any(mp[[g]]$coverage < 0.1)) {
      coverage_vech <- lav_mat_vech(mp[[g]]$coverage, diagonal = FALSE)
      small_idx <- which(coverage_vech < 0.1)
      if (all(coverage_vech[small_idx] == 0)) {
      # 0.6-18: no warning --> this could be due to missing by design
          # 0.6-20: give warning anyway (as EM is ignoring this)
          lav_msg_warn(gettext(
            "due to missing values, some pairwise combinations have zero
             coverage; the corresponding covariances are not identified;
             use lavInspect(fit, \"coverage\") to investigate."))
    } else {
          lav_msg_warn(gettext(
            "due to missing values, some pairwise combinations have less than
            10% coverage; use lavInspect(fit, \"coverage\") to investigate."))
      }
      }
      # in case we had observations with only missings
      nobs[[g]] <- NROW(x[[g]]) - length(mp[[g]]$empty.idx)
    } # missing
  } # groups, at first level

  # sampling weights, again
  if (is.null(sampling_weights)) {
    sampling_weights <- character(0L)
  } else {
    # check if we need normalization
    if (sampling_weights_normalization == "none") {
      # nothing to do
    } else if (sampling_weights_normalization == "total") {
      sum_weights <- sum(unlist(weights))
      ntotal <- sum(unlist(nobs))
      for (g in 1:ngroups) {
        wt <- weights[[g]]
        wt2 <- wt / sum_weights * ntotal
        weights[[g]] <- wt2
      }
    } else if (sampling_weights_normalization == "group") {
      for (g in 1:ngroups) {
        wt <- weights[[g]]
        wt2 <- wt / sum(wt) * nobs[[g]]
        weights[[g]] <- wt2
      }
    } else {
      lav_msg_stop(gettext(
        "sampling.weights.normalization should be total, group or none."))
    }
  }

  # block.labels
  block_label <- character(0L)
  if (length(group_label) > 0L && length(level_label) == 0L) {
    block_label <- group_label
  } else if (length(level_label) > 0L && length(group_label) == 0L) {
    block_label <- level_label
  } else if (length(group_label) > 0L &&
    length(level_label) > 0L) {
    block_label <- paste(rep(group_label, each = length(level_label)),
      rep(level_label, times = length(group_label)),
      sep = "."
    )
  }


  lav_data <- new("lavData",
    data.type = "full",
    ngroups = ngroups,
    group = group,
    nlevels = nlevels,
    cluster = cluster,
    group.label = group_label,
    level.label = level_label,
    block.label = block_label,
    std.ov = std_ov,
    nobs = nobs,
    norig = norig,
    ov.names = ov_names,
    ov.names.x = ov_names_x,
    ov.names.l = ov_names_l,
    ov.names.aux = if (length(ov_names_aux) > 0L) {
      rep(list(as.character(ov_names_aux)), ngroups)
    } else {
      vector("list", length = ngroups)
    },
    # ov.types        = ov.types,
    # ov.idx          = ov.idx,
    ordered = as.character(ordered),
    weights = weights,
    sampling.weights = sampling_weights,
    ov = ov,
    ov.aux = ov_aux,
    case.idx = case_idx,
    missing = missing,
    X = x,
    eXo = exo,
    aux = aux,
    Mp = mp,
    Rp = rp,
    Lp = lp
  )
  lav_data
}
