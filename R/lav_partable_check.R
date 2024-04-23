# check if the partable is complete/consistent
# we may have added intercepts/variances (user = 0), fixed to zero
lav_partable_check <- function(partable, categorical = FALSE, warn = TRUE) {
  check <- TRUE

  # check for empy table - or should we WARN?
  if (length(partable$lhs) == 0) {
    return(check)
  }

  # get observed/latent variables
  ov.names <- vnames(partable, "ov.nox") # no need to specify exo??
  lv.names <- vnames(partable, "lv")
  all.names <- c(ov.names, lv.names)
  ov.names.ord <- vnames(partable, "ov.ord")

  nlevels <- lav_partable_nlevels(partable)

  # if categorical, we should have some ov.names.ord
  if (categorical && length(ov.names.ord) == 0L) {
    check <- FALSE
    if (warn) {
      lav_msg_warn(gettext("parameter table does not contain thresholds"))
    }
  }

  # we should have a (residual) variance for *each* ov/lv
  # note: if lavaanify() has been used, this is always TRUE
  var.idx <- which(partable$op == "~~" &
    partable$lhs == partable$rhs)
  missing.idx <- which(is.na(match(all.names, partable$lhs[var.idx])))
  if (length(missing.idx) > 0L) {
    check <- FALSE
    if (warn) {
      lav_msg_warn(gettextf(
        "parameter table does not contain (residual) variances for
        one or more variables: %s",
        lav_msg_view(all.names[missing.idx])))
    }
  }

  # meanstructure?
  meanstructure <- any(partable$op == "~1")

  # if meanstructure, check for missing intercepts
  # note if lavaanify() has been used, this is always TRUE
  if (meanstructure) {
    # we should have a intercept for *each* ov/lv
    int.idx <- which(partable$op == "~1")
    missing.idx <- which(is.na(match(all.names, partable$lhs[int.idx])))
    if (length(missing.idx) > 0L) {
      check <- FALSE
      if (warn) {
        lav_msg_warn(gettextf(
          "parameter table does not contain intercepts
          for one or more variables: %s",
          lav_msg_view(all.names[missing.idx])))
      }
    }
  }

  # ok, now the 'real' checks

  # do we have added (residual) variances (user = 0) that are fixed to zero?
  # this is not necessarily problematic!
  # eg. in latent change score models
  # therefore, we do NOT give a warning

  # var.fixed <- which(partable$op == "~~" &
  #                   partable$lhs == partable$rhs &
  #                   partable$user == 0 &
  #                   partable$free == 0)
  # if(length(var.fixed) > 0L) {
  #    check <- FALSE
  #     if(warn) {
  #        warning("lavaan WARNING: missing (residual) variances are set to",
  #        " zero: [", paste(partable$lhs[var.fixed],  collapse = " "), "]")
  #    }
  # }

  # do we have added intercepts (user = 0) that are fixed to zero?
  # this is not necessarily problematic; perhaps only for
  # exogenous variables?
  ov.ind <- unique(partable$rhs[partable$op == "=~"])
  lv.names <- unique(partable$lhs[partable$op == "=~"])
  int.fixed <- which(partable$op == "~1" &
    partable$user == 0L &
    partable$free == 0L &
    partable$ustart == 0L &
    # ignore block/group 1 -- typically within level exo
    !(partable$block %% nlevels == 1L) &
    # do not include factors
    !partable$lhs %in% lv.names &
    # do not include ordered variables
    !partable$lhs %in% ov.names.ord &
    # do not include indicators
    !partable$lhs %in% ov.ind)

  if (length(int.fixed) > 0L) {
    check <- FALSE
    if (warn) {
      lav_msg_warn(gettext("automatically added intercepts are set to zero:"),
                   lav_msg_view(partable$lhs[int.fixed]))
    }
  }

  # return check code
  check
}
