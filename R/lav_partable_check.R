# check if the partable is complete/consistent
# for example, we may have added intercepts/variances (user = 0), fixed to zero
lav_pt_check <- function(partable, categorical = FALSE) {
  check <- TRUE

  # check for empty table - or should we WARN?
  if (length(partable$lhs) == 0) {
    return(check)
  }

  # get observed/latent variables
  ov_names <- lav_pt_vnames(partable, "ov.nox") # no need to specify exo??
  lv_names <- lav_pt_vnames(partable, "lv")
  lv_names_c <- lav_pt_vnames(partable, "lv.composite")
  lv_names_noc <- lv_names[!lv_names %in% lv_names_c]
  all_names <- c(ov_names, lv_names_noc)
  ov_names_ord <- lav_pt_vnames(partable, "ov.ord")

  nlevels <- lav_pt_nlevels(partable)

  # composites: within a single block, an indicator may not serve a composite
  # (<~) and a common factor (=~) at the same time (each indicator relates to
  # one construct only; see Schamberger et al.). Catch this here with a clear
  # message, otherwise it crashes later (non-positive-definite Sigma / NA in
  # loglik). NOTE: the check is per-block, because in a multilevel (or
  # multigroup) model the same observed variable may legitimately be a composite
  # indicator in one block and a common-factor indicator in another.
  if (length(lv_names_c) > 0L) {
    blk <- if (!is.null(partable$block)) partable$block else
      rep(1L, length(partable$lhs))
    for (b in unique(blk)) {
      cind <- unique(partable$rhs[partable$op == "<~" & blk == b])
      find <- unique(partable$rhs[partable$op == "=~" & blk == b])
      both <- cind[cind %in% find]
      if (length(both) > 0L) {
        lav_msg_stop(gettextf(
          "indicator(s) used by both a composite (<~) and a common factor (=~) %s: %s",
          "in the same block", lav_msg_view(both)))
      }
    }
  }

  # if categorical, we should have some ov.names.ord
  if (categorical && length(ov_names_ord) == 0L) {
    check <- FALSE
    lav_msg_warn(gettext("parameter table does not contain thresholds"))
  }

  # we should have a (residual) variance for *each* ov/lv
  # note: if lav_model_pt() has been used, this is always TRUE
  var_idx <- which(partable$op == "~~" &
    partable$lhs == partable$rhs & !partable$lhs %in% lv_names_c)
  missing_idx <- which(is.na(match(all_names, partable$lhs[var_idx])))
  if (length(missing_idx) > 0L) {
    check <- FALSE
    lav_msg_warn(gettextf(
      "parameter table does not contain (residual) variances for
      one or more variables: %s",
      lav_msg_view(all_names[missing_idx])))
  }

  # meanstructure?
  meanstructure <- any(partable$op == "~1")

  # if meanstructure, check for missing intercepts
  # note if lav_model_pt() has been used, this is always TRUE
  if (meanstructure) {
    # we should have an intercept for *each* ov/lv
    int_idx <- which(partable$op == "~1")
    missing_idx <- which(is.na(match(all_names, partable$lhs[int_idx])))
    if (length(missing_idx) > 0L) {
      check <- FALSE
      lav_msg_warn(gettextf(
        "parameter table does not contain intercepts
        for one or more variables: %s",
        lav_msg_view(all_names[missing_idx])))
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
  ov_ind <- unique(partable$rhs[partable$op == "=~"])
  lv_names <- unique(partable$lhs[partable$op == "=~"])
  int_fixed <- which(partable$op == "~1" &
    partable$user == 0L &
    partable$free == 0L &
    partable$ustart == 0L &
    # ignore block/group 1 -- typically within level exo
    !(partable$block %% nlevels == 1L) &
    # do not include factors
    !partable$lhs %in% lv_names &
    # do not include ordered variables
    !partable$lhs %in% ov_names_ord &
    # do not include indicators
    !partable$lhs %in% ov_ind &
    # do not include lv composites
    !partable$lhs %in% lv_names_c)

  if (length(int_fixed) > 0L) {
    check <- FALSE
    lav_msg_warn(gettext("automatically added intercepts are set to zero:"),
                 lav_msg_view(partable$lhs[int_fixed]))
  }

  # return check code
  check
}
