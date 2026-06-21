lav_pt_flat <- function(flat = NULL,
                              blocks = "group",
                              block_id = NULL,
                              meanstructure = FALSE,
                              int_ov_free = FALSE,
                              int_lv_free = FALSE,
                              orthogonal = FALSE,
                              orthogonal_y = FALSE,
                              orthogonal_x = FALSE,
                              orthogonal_efa = FALSE,
                              std_lv = FALSE,
                              correlation = FALSE,
                              composites = TRUE,
                              conditional_x = FALSE,
                              fixed_x = TRUE,
                              parameterization = "delta",
                              auto_fix_first = FALSE,
                              marker = NULL,
                              auto_fix_single = FALSE,
                              auto_var = FALSE,
                              auto_cov_lv_x = FALSE,
                              auto_cov_y = FALSE,
                              auto_th = FALSE,
                              auto_delta = FALSE,
                              auto_efa = FALSE,
                              var_table = NULL,
                              group_equal = NULL,
                              group_w_free = FALSE,
                              ngroups = 1L,
                              nthresholds = NULL,
                              ov_names_x_block = NULL) {
  categorical <- FALSE

  ### tmp.default elements: parameters that are typically not specified by
  ###                   users, but should typically be considered,
  ###                   either free or fixed

  # extract `names' of various types of variables:
  lv_names <- lav_pt_vnames(flat, type = "lv") # latent variables
  # lv.names.r   <- lav_pt_vnames(FLAT, type="lv.regular")
  # regular latent variables
  if (composites) {
    lv_names_f <- character(0L)
    lv_names_c <- lav_pt_vnames(flat, type = "lv.composite")
    ov_ind_c <- lav_pt_vnames(flat, type = "ov.cind")
    lv_names_noc <- lv_names[!lv_names %in% lv_names_c]
  } else {
    lv_names_c <- character(0L)
    ov_ind_c <- character(0L)
    lv_names_f <- lav_pt_vnames(flat, type = "lv.formative")
    lv_names_noc <- lv_names
  }

  # formative latent variables
  ov_names <- lav_pt_vnames(flat, type = "ov")
  # observed variables
  ov_names_x <- lav_pt_vnames(flat, type = "ov.x")
  # exogenous x covariates
  lv_names_int <- lav_pt_vnames(flat, type = "lv.interaction")
  # lv interactions

  if (is.null(ov_names_x_block)) {
    ov_names_x_block <- ov_names_x
  }
  ov_names_nox <- lav_pt_vnames(flat, type = "ov.nox")
  lv_names_x <- lav_pt_vnames(flat, type = "lv.x") # exogenous lv
  ov_names_y <- lav_pt_vnames(flat, type = "ov.y") # dependent ov
  lv_names_y <- lav_pt_vnames(flat, type = "lv.y") # dependent lv
  lv_names_efa <- lav_pt_vnames(flat, type = "lv.efa")
  # lvov.names.y <- c(ov.names.y, lv.names.y)
  lvov_names_y <- c(lv_names_y, ov_names_y)

  # get 'ordered' variables, either from FLAT or var_table
  ov_names_ord1 <- lav_pt_vnames(flat, type = "ov.ord")
  # check if we have "|" for exogenous variables
  if (length(ov_names_ord1) > 0L) {
    idx <- which(ov_names_ord1 %in% ov_names_x)
    if (length(idx) > 0L) {
      lav_msg_warn(gettext("thresholds are defined for exogenous variables:"),
                   lav_msg_view(ov_names_ord1[idx], "none"))
    }
  }

  # check data
  if (!is.null(var_table)) {
    ov_names_ord2 <-
      as.character(var_table$name[var_table$type == "ordered"])
    # remove fixed.x variables
    idx <- which(ov_names_ord2 %in% ov_names_x)
    if (length(idx) > 0L) {
      ov_names_ord2 <- ov_names_ord2[-idx]
    }

    # remove those that do appear in the model syntax
    idx <- which(!ov_names_ord2 %in% ov_names)
    if (length(idx) > 0L) {
      ov_names_ord2 <- ov_names_ord2[-idx]
    }
  } else {
    ov_names_ord2 <- character(0L)
  }

  # check nthresholds, if it is a named vector
  ov_names_ord3 <- character(0L)
  if (!is.null(nthresholds)) {
    if (!is.null(var_table)) {
      lav_msg_stop(gettext(
        "the var_table and nthresholds arguments should not be used together."))
    }
    if (!is.numeric(nthresholds)) {
      lav_msg_stop(gettext("nthresholds should be a named vector of integers."))
    }
    nth_names <- names(nthresholds)
    if (!is.null(nth_names)) {
      ov_names_ord3 <- nth_names
    } else {
      # if nthresholds is just a number, all is good; otherwise it
      # should be a names vector
      if (length(nthresholds) > 1L) {
        lav_msg_warn(gettext("nthresholds must be a named vector of integers."))
      }
      # just a single number -> assume ALL y variables are ordered
      ov_names_ord3 <- ov_names_nox
    }
  }

  # final ov.names.ord
  tmp <- unique(c(ov_names_ord1, ov_names_ord2, ov_names_ord3))
  ov_names_ord <- ov_names[ov_names %in% tmp]

  # if we have the "|" in the model syntax, check the number of thresholds
  # if(!is.null(var_table) && length(ov.names.ord1) > 0L) {
  #    for(o in ov.names.ord1) {
  #        nth <- var_table$nlev[ var_table$name == o ] - 1L
  #        nth.in.partable <- sum(FLAT$op == "|" & FLAT$lhs == o)
  #        if(nth != nth.in.partable) {
  #            stop("lavaan ERROR: expected ", max(0,nth),
  #                 " threshold(s) for variable ",
  #                 sQuote(o), "; syntax contains ", nth.in.partable, "\n")
  #        }
  #    }
  # }

  if (length(ov_names_ord) > 0L) {
    categorical <- TRUE
  }

  # do we have any EFA lv's? they need special treatment if auto.efa = TRUE
  if (!is.null(flat$efa) && auto_efa) {
    lv_names_efa <- unique(flat$lhs[flat$op == "=~" &
      nchar(flat$efa) > 0L])
    # remove them from lv.names.x
    # if(length(lv.names.x) > 0L) {
    #    both.idx <- which(lv.names.x %in% lv.names.efa)
    #    if(length(both.idx) > 0L) {
    #        lv.names.x <- lv.names.x[ -both.idx ]
    #    }
    # }

    # remove them from lvov.names.y
    # if(length(lvov.names.y) > 0L) {
    #    both.idx <- which(lvov.names.y %in% lv.names.efa)
    #    if(length(both.idx) > 0L) {
    #        lvov.names.y <- lvov.names.y[ -both.idx ]
    #    }
    # }
  } else {
    lv_names_efa <- character(0)
  }

  lhs <- rhs <- character(0)

  # 1. THRESHOLDS (based on var_table)
  #    NOTE: - new in 0.5-18: ALWAYS include threshold parameters in partable,
  #            but only free them if auto.th = TRUE
  #          - [only ov.names.ord2, because ov.names.ord1 are already
  #            in tmp.user and we only need to add 'default' parameters here]
  #            (not any longer: we create them for ALL ordered var (0.6-12)
  nth <- 0L
  # if(auto.th && length(ov.names.ord2) > 0L) {
  # if(length(ov.names.ord2) > 0L) {
  if (length(ov_names_ord) > 0L) {
    # for(o in ov.names.ord2) {
    for (o in ov_names_ord) {
      if (!is.null(var_table)) {
        nth <- var_table$nlev[var_table$name == o] - 1L
      } else if (!is.null(nthresholds)) {
        if (length(nthresholds) == 1L && is.null(nth_names)) {
          nth <- nthresholds
        } else {
          # we can assume nthresholds is a named vector
          nth <- unname(nthresholds[o])
          if (is.na(nth)) {
            lav_msg_stop(gettextf("ordered variable %s not found in the
                                  named vector nthresholds.", o))
          }
        }
      }
      if (nth < 1L) next
      lhs <- c(lhs, rep(o, nth))
      rhs <- c(rhs, paste("t", seq_len(nth), sep = ""))
    }
    nth <- length(lhs)
  }

  # 2. default (residual) variances and covariances

  # a) (residual) VARIANCES (all ov's except exo, and all lv's)
  # NOTE: change since 0.5-17: we ALWAYS include the vars in the
  #       parameter table; but only if auto.var = TRUE, we set them free
  # if(auto.var) {
  ov_var <- ov_names_nox
  # auto-remove ordinal variables
  # idx <- match(ov.names.ord, ov.var)
  # if(length(idx)) ov.var <- ov.var[-idx]
  lhs <- c(lhs, ov_var, lv_names)
  rhs <- c(rhs, ov_var, lv_names)
  # }

  # b) `independent` latent variable COVARIANCES (lv.names.x)
  if (auto_cov_lv_x && length(lv_names_x) > 1L) {
    tmp <- utils::combn(lv_names_x, 2)
    lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
    rhs <- c(rhs, tmp[2, ])
  }

  # c) `dependent` latent variables COVARIANCES (lv.y.idx + ov.y.lv.idx)
  if (auto_cov_y && length(lvov_names_y) > 1L) {
    tmp <- utils::combn(lvov_names_y, 2L)
    lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
    rhs <- c(rhs, tmp[2, ])
  }

  # d) exogenous x covariates: VARIANCES + COVARIANCES
  if ((nx <- length(ov_names_x)) > 0L) {
    if (conditional_x) {
      # new in 0.6-12: we make a distinction between ov.names.x and
      # ov.names.x.block: we treat them 'separately' (with no covariances
      # among them)
      # but we add 'regressions' instead (see below)
      ov_names_x1 <- ov_names_x[!ov_names_x %in% ov_names_x_block]
      ov_names_x2 <- ov_names_x_block
      nx1 <- length(ov_names_x1) # split x
      nx2 <- length(ov_names_x2) # regular  x
      if (nx1 > 0L) {
        idx <- lower.tri(matrix(0, nx1, nx1), diag = TRUE)
        lhs <- c(lhs, rep(ov_names_x1, each = nx1)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov_names_x1, times = nx1)[idx])
      }
      if (nx2 > 0L) {
        idx <- lower.tri(matrix(0, nx2, nx2), diag = TRUE)
        lhs <- c(lhs, rep(ov_names_x2, each = nx2)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov_names_x2, times = nx2)[idx])
      }
    } else {
      idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
      lhs <- c(lhs, rep(ov_names_x, each = nx)[idx]) # fill upper.tri
      rhs <- c(rhs, rep(ov_names_x, times = nx)[idx])
    }
  }

  # e) indicators of composites: COVARIANCES
  #    but only within/intra blocks
  if ((ncx <- length(ov_ind_c)) > 0L) {
    # create W1
    w1 <- matrix(0, length(ov_ind_c), length(lv_names_c))
    c_idx <- which(flat$op == "<~")
    w1[cbind(match(flat$rhs[c_idx], ov_ind_c),
             match(flat$lhs[c_idx], lv_names_c))] <- 1
    w1w1 <- tcrossprod(w1)
    w1w1[upper.tri(w1w1, diag = TRUE)] <- 0 # keep lower.tri only
    if (ncx > 1L) {
      lhs <- c(lhs, ov_ind_c[col(w1w1)[as.logical(w1w1)]])
      rhs <- c(rhs, ov_ind_c[row(w1w1)[as.logical(w1w1)]])
    }
  }

  # f) efa latent variables COVARIANCES; only needed for 'mediators'
  #    (not in lv.names.x, not in lv.names.y) -- added in 0.6-18
  if (auto_efa && length(lv_names_efa) > 1L) {
    efa_values <- lav_pt_efa_values(flat)
    for (set in efa_values) {
      # correlated factors within each set
      this_set_lv <- unique(flat$lhs[flat$op == "=~" &
        !flat$lhs %in% lv_names_x &
        !flat$lhs %in% lv_names_y &
        flat$efa == set])
      if (length(this_set_lv) > 0L) {
        tmp <- utils::combn(this_set_lv, 2)
        lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
        rhs <- c(rhs, tmp[2, ])
      }
    }
  }

  # create 'op' (thresholds come first, then variances)
  op <- rep("~~", length(lhs))
  op[seq_len(nth)] <- "|"

  # LATENT RESPONSE SCALES (DELTA)
  #    NOTE: - new in 0.5-19: ALWAYS include scaling parameters in partable,
  #            but only free them if auto.delta = TRUE (and parameterization
  #            is "delta"
  # if(auto.delta && auto.th && length(ov.names.ord) > 0L &&
  #   # length(lv.names) > 0L &&
  #   (ngroups > 1L || any(FLAT$op == "~*~") || parameterization == "theta")) {
  if (length(ov_names_ord) > 0L) {
    lhs <- c(lhs, ov_names_ord)
    rhs <- c(rhs, ov_names_ord)
    op <- c(op, rep("~*~", length(ov_names_ord)))
  }

  # same for correlation structures, but now for ALL variables
  if (!categorical && correlation) {
    lhs <- c(lhs, ov_names)
    rhs <- c(rhs, ov_names)
    op <- c(op, rep("~*~", length(ov_names)))
  }

  # 3. INTERCEPTS
  if (meanstructure) {
    # if(conditional.x) {
    #    ov.int <- ov.names.nox
    # } else {
    ov_int <- ov_names
    # }
    # auto-remove ordinal variables
    # idx <- which(ov.int %in% ov.names.ord)
    # if(length(idx)) ov.int <- ov.int[-idx]

    int_lhs <- c(ov_int, lv_names)
    lhs <- c(lhs, int_lhs)
    rhs <- c(rhs, rep("", length(int_lhs)))
    op <- c(op, rep("~1", length(int_lhs)))
  }

  # 4. REGRESSIONS
  if (conditional_x) {
    # new in 0.6-12: we make a distinction between ov.names.x and
    # ov.names.x.block: we treat them 'separately' (with no covariances
    # among them)
    # but we add 'regressions' instead!
    ov_names_x1 <- ov_names_x[!ov_names_x %in% ov_names_x_block]
    ov_names_x2 <- ov_names_x_block
    nx1 <- length(ov_names_x1) # split x
    nx2 <- length(ov_names_x2) # regular  x
    if (nx1 > 0L && nx2 > 0L) {
      # add regressions for split-x ~ regular-x
      lhs <- c(lhs, rep(ov_names_x1, times = nx2))
      op <- c(op, rep("~", nx2 * nx1))
      rhs <- c(rhs, rep(ov_names_x2, each = nx1))
    }
  }

  # free group weights
  if (group_w_free) {
    lhs <- c(lhs, "group")
    rhs <- c(rhs, "w")
    op <- c(op, "%")
  }

  tmp_default <- data.frame(
    lhs = lhs, op = op, rhs = rhs,
    mod.idx = rep(0L, length(lhs)),
    stringsAsFactors = FALSE
  )


  # 4. USER: user-specified elements
  lhs <- flat$lhs
  op <- flat$op
  rhs <- flat$rhs
  mod_idx <- flat$mod.idx

  lv_names <- lav_pt_vnames(flat, type = "lv") # latent variables
  ov_names <- lav_pt_vnames(flat, type = "ov") # observed variables
  tmp_user <- data.frame(
    lhs = lhs, op = op, rhs = rhs, mod.idx = mod_idx,
    stringsAsFactors = FALSE
  )

  # check for duplicated elements in tmp.user
  tmp_tmp <- tmp_user[, 1:3]
  idx <- which(duplicated(tmp_tmp))
  if (length(idx) > 0L) {
    txt <- sapply(seq_along(idx), function(i) {
      paste(
        "    ", tmp_tmp[idx[i], "lhs"],
        tmp_tmp[idx[i], "op"],
        tmp_tmp[idx[i], "rhs"]
      )
    })
    lav_msg_warn(gettext(
      "duplicated elements in model syntax have been ignored:"),
      lav_msg_view(txt, "none"))
    tmp_user <- tmp_user[-idx, ]
  }

  # check for duplicated elements in tmp.default
  # - FIXME: can we not avoid this somehow??
  # - for example, if the user model includes 'x1 ~~ x1'
  #   or 'x1 ~ 1'
  # - remove them from tmp.default
  tmp_tmp <- rbind(tmp_default[, 1:3], tmp_user[, 1:3])
  idx <- which(duplicated(tmp_tmp, fromLast = TRUE))
  # idx should be in tmp.default
  if (length(idx)) {
    for (i in idx) {
      flat_idx <- which(tmp_user$lhs == tmp_default$lhs[i] &
        tmp_user$op == tmp_default$op[i] &
        tmp_user$rhs == tmp_default$rhs[i])
      if (length(flat_idx) != 1L) {
        cat("[lavaan DEBUG] idx in tmp.tmp: i = ", i, "\n")
        print(tmp_tmp[i, ])
        cat("[lavaan DEBUG] idx in tmp.default: i = ", i, "\n")
        print(tmp_default[i, ])
        cat("[lavaan DEBUG] flat.idx:")
        print(flat_idx)
      }
    }
    tmp_default <- tmp_default[-idx, ]
  }

  # now that we have removed all duplicated elements, we can construct
  # the tmp.list for a single group/block
  lhs <- c(tmp_user$lhs, tmp_default$lhs)
  op <- c(tmp_user$op, tmp_default$op)
  rhs <- c(tmp_user$rhs, tmp_default$rhs)
  user <- c(
    rep(1L, length(tmp_user$lhs)),
    rep(0L, length(tmp_default$lhs))
  )
  mod_idx <- c(tmp_user$mod.idx, tmp_default$mod.idx)

  # by default: everything is free!
  free <- rep(1L, length(lhs))
  ustart <- rep(as.numeric(NA), length(lhs))
  # label   <- paste(lhs, op, rhs, sep="")
  label <- rep(character(1), length(lhs))
  exo <- rep(0L, length(lhs))

  # 0a. if auto.th = FALSE, set fix the thresholds
  if (!auto_th) {
    th_idx <- which(op == "|" & user == 0L)
    free[th_idx] <- 0L
  }

  # 0b. if auto.var = FALSE, set the unspecified variances to zero
  if (!auto_var) {
    var_idx <- which(op == "~~" &
      lhs == rhs &
      !lhs %in% ov_ind_c &
      user == 0L)
    ustart[var_idx] <- 0.0
    free[var_idx] <- 0L
  } else if (length(lv_names_f) > 0L) {
    # 'formative' (residual) variances are set to zero by default
    var_idx <- which(op == "~~" &
      lhs == rhs &
      lhs %in% lv_names_f &
      user == 0L)
    ustart[var_idx] <- 0.0
    free[var_idx] <- 0L
  }

  # 0c. for the ~~ for composite indicators: currently ALWAYS fixed
  #     todo: create an option to free them anyway
  if (length(ov_ind_c) > 0) {
    var_idx <- which(op == "~~" & lhs %in% ov_ind_c)
    ustart[var_idx] <- as.numeric(NA)
    free[var_idx] <- 0L
  }

  # 0d. variances for composites: ALWAYS fixed (should be set later
  #     by lav_lisrel_comp_set_intresvar
  if (length(lv_names_c) > 0) {
    var_idx <- which(op == "~~" & lhs %in% lv_names_c & lhs == rhs)
    ustart[var_idx] <- as.numeric(NA)
    free[var_idx] <- 0L
  }

  # 0e. instruments: ALWAYS nonfree
  iv_idx <- which(op == "|~")
  ustart[iv_idx] <- 0.0
  free[iv_idx] <- 0L

  # 1. fix metric of regular latent variables
  if (std_lv) {
    # fix metric by fixing the variance of the latent variable
    lv_var_idx <- which(op == "~~" &
      lhs %in% lv_names & lhs == rhs)
    ustart[lv_var_idx] <- 1.0
    free[lv_var_idx] <- 0L
  }
  if (auto_efa && length(lv_names_efa) > 0L) {
    # fix lv variances of efa blocks to unity
    lv_var_idx <- which(op == "~~" &
      lhs %in% lv_names_efa & lhs == rhs)
    ustart[lv_var_idx] <- 1.0
    free[lv_var_idx] <- 0L
  }
  if (auto_fix_first) {
    # fix metric by fixing the loading of the first indicator
    # (but not for efa factors)
    #
    # if 'marker' is provided (a named vector lv -> indicator), fix the
    # loading of that indicator instead of the first one; this is used by
    # the bad.marker.crit mechanism to switch to another marker if the first
    # indicator turns out to be a poor item (see lav_pt_marker_adapt())
    mm_idx <- which(op == "=~" & !(lhs %in% lv_names_efa))
    if (is.null(marker)) {
      first_idx <- mm_idx[which(!duplicated(lhs[mm_idx]))]
      ustart[first_idx] <- 1.0
      free[first_idx] <- 0L
    } else {
      for (lv in unique(lhs[mm_idx])) {
        lv_rows <- mm_idx[lhs[mm_idx] == lv]
        marker_row <- lv_rows[1L] # default: first indicator
        if (lv %in% names(marker) && !is.na(marker[[lv]])) {
          tmp_row <- lv_rows[rhs[lv_rows] == marker[[lv]]]
          if (length(tmp_row) > 0L) {
            marker_row <- tmp_row[1L]
          }
        }
        ustart[marker_row] <- 1.0
        free[marker_row] <- 0L
      }
    }
    if (composites && length(lv_names_c) > 0L) {
      mm_idx <- which(op == "<~")
      first_idx <- mm_idx[which(!duplicated(lhs[mm_idx]))]
      ustart[first_idx] <- 1.0
      free[first_idx] <- 0L
    }
  }

  # 2. fix residual variance of single indicators to zero
  if (auto_var && auto_fix_single) {
    mm_idx <- which(op == "=~")
    tmp_t <- table(lhs[mm_idx])
    if (any(tmp_t == 1L)) {
      # ok, we have a LV with only a single indicator
      lv_names_single <- names(tmp_t)[tmp_t == 1L]
      # get corresponding indicator if unique
      lhs_mm <- lhs[mm_idx]
      rhs_mm <- rhs[mm_idx]
      single_ind <- rhs_mm[which(lhs_mm %in% lv_names_single &
        lhs_mm != rhs_mm & # exclude phantom
        !(duplicated(rhs_mm) |
          duplicated(rhs_mm, fromLast = TRUE)))]
      # is the indicator unique?
      if (length(single_ind) > 0L) {
        var_idx <- which(op == "~~" & lhs %in% single_ind &
          rhs %in% single_ind &
          lhs == rhs &
          user == 0L)
        ustart[var_idx] <- 0.0
        free[var_idx] <- 0L
      }
    }
  }

  # 3. orthogonal = TRUE?
  if (orthogonal) {
    lv_cov_idx <- which(op == "~~" &
      lhs %in% lv_names &
      rhs %in% lv_names &
      lhs != rhs &
      user == 0L)
    ustart[lv_cov_idx] <- 0.0
    free[lv_cov_idx] <- 0L
  }
  # 3b. orthogonal.y = TRUE?
  if (orthogonal_y) {
    lv_cov_idx <- which(op == "~~" &
      lhs %in% lv_names_y &
      rhs %in% lv_names_y &
      lhs != rhs &
      user == 0L)
    ustart[lv_cov_idx] <- 0.0
    free[lv_cov_idx] <- 0L
  }
  # 3c. orthogonal.x = TRUE?
  if (orthogonal_x) {
    lv_cov_idx <- which(op == "~~" &
      lhs %in% lv_names_x &
      rhs %in% lv_names_x &
      lhs != rhs &
      user == 0L)
    ustart[lv_cov_idx] <- 0.0
    free[lv_cov_idx] <- 0L
  }
  # 3d. orthogonal.efa = TRUE?
  if (orthogonal_efa) {
    lv_cov_idx <- which(op == "~~" &
      lhs %in% lv_names_efa &
      rhs %in% lv_names_efa &
      lhs != rhs &
      user == 0L)
    ustart[lv_cov_idx] <- 0.0
    free[lv_cov_idx] <- 0L
  }

  # 4. intercepts
  if (meanstructure) {
    if (categorical) {
      # zero intercepts/means ordinal variables
      ov_int_idx <- which(op == "~1" &
        lhs %in% ov_names_ord &
        user == 0L)
      ustart[ov_int_idx] <- 0.0
      free[ov_int_idx] <- 0L
    }
    if (int_ov_free == FALSE) {
      # zero intercepts/means observed variables
      ov_int_idx <- which(op == "~1" &
        lhs %in% ov_names &
        user == 0L)
      ustart[ov_int_idx] <- 0.0
      free[ov_int_idx] <- 0L
    }
    if (int_lv_free == FALSE) {
      # zero intercepts/means latent variables
      lv_int_idx <- which(op == "~1" &
        lhs %in% lv_names &
        user == 0L)
      ustart[lv_int_idx] <- 0.0
      free[lv_int_idx] <- 0L
    }
    # 4b. fixed effect (only if we have random slopes)
    if (!is.null(flat$rv) && any(nchar(flat$rv) > 0L)) {
      lv_names_rv <- lav_pt_vnames(flat, "lv.rv")
      lv_rv_idx <- which(op == "~1" &
        lhs %in% lv_names_rv &
        user == 0L)
      ustart[lv_rv_idx] <- as.numeric(NA)
      free[lv_rv_idx] <- 1L
    }
    if (length(lv_names_int) > 0L) {
      lv_int_idx <- which(op == "~1" &
        lhs %in% lv_names_int &
        user == 0L)
      ustart[lv_int_idx] <- as.numeric(NA)
      free[lv_int_idx] <- 1L
    }
    # composites: always non-free, but with ustart = NA; value should be
    # filled in later as a function of the other parameters
    if (length(lv_names_c) > 0L) {
      c_int_idx <- which(op == "~1" & lhs %in% lv_names_c & user == 0L)
      ustart[c_int_idx] <- as.numeric(NA)
      free[c_int_idx] <- 0L
    }
  }

  # 4b. fixed effect (only if we have random slopes)
  # if(!is.null(FLAT$rv)) {
  #    }

  # 5. handle exogenous `x' covariates
  # usually, ov.names.x.block == ov.names.x
  # except if multilevel, where 'split' ov.x are treated as endogenous

  # 5a conditional.x = FALSE
  if (!conditional_x && fixed_x && length(ov_names_x_block) > 0) {
    # 1. variances/covariances
    exo_var_idx <- which(op == "~~" &
      rhs %in% ov_names_x_block &
      lhs %in% ov_names_x_block &
      user == 0L)
    ustart[exo_var_idx] <- as.numeric(NA) # should be overridden later!
    free[exo_var_idx] <- 0L
    exo[exo_var_idx] <- 1L

    # 2. intercepts
    exo_int_idx <- which(op == "~1" &
      lhs %in% ov_names_x_block &
      user == 0L)
    ustart[exo_int_idx] <- as.numeric(NA) # should be overridden later!
    free[exo_int_idx] <- 0L
    exo[exo_int_idx] <- 1L
  }

  # 5a-bis. conditional.x = TRUE
  if (conditional_x && length(ov_names_x) > 0L) {
    # 1. variances/covariances
    exo_var_idx <- which(op == "~~" &
      rhs %in% ov_names_x &
      lhs %in% ov_names_x &
      user == 0L)
    if (fixed_x) {
      ustart[exo_var_idx] <- as.numeric(NA) # should be overridden later!
      free[exo_var_idx] <- 0L
    }
    exo[exo_var_idx] <- 1L

    # 2. intercepts
    exo_int_idx <- which(op == "~1" &
      lhs %in% ov_names_x &
      user == 0L)
    if (fixed_x) {
      ustart[exo_int_idx] <- as.numeric(NA) # should be overridden later!
      free[exo_int_idx] <- 0L
    }
    exo[exo_int_idx] <- 1L

    # 3. regressions ov + lv
    exo_reg_idx <- which(op %in% c("~", "<~") &
      lhs %in% c(lv_names, ov_names_nox) &
      rhs %in% ov_names_x)
    exo[exo_reg_idx] <- 1L

    # 3b regression split.x ~ regular.x
    exo_reg2_idx <- which(op %in% c("~", "<~") &
      lhs %in% ov_names_x &
      rhs %in% ov_names_x)
    if (fixed_x) {
      ustart[exo_reg2_idx] <- as.numeric(NA) # should be overridden later!
      free[exo_reg2_idx] <- 0L
    }
    exo[exo_reg2_idx] <- 1L
  }

  # 5b. residual variances of ordinal variables?
  if (length(ov_names_ord) > 0L) {
    ord_idx <- which(lhs %in% ov_names_ord &
      op == "~~" &
      user == 0L & ## New in 0.6-1
      lhs == rhs)
    ustart[ord_idx] <- 1L ## FIXME!! or 0?? (0 breaks ex3.12)
    free[ord_idx] <- 0L
  }

  # 5c latent response scales of ordinal variables?
  #    by default, all fixed to 1.0
  if (length(ov_names_ord) > 0L) {
    delta_idx <- which(op == "~*~" &
      user == 0L) ## New in 0.6-1
    ustart[delta_idx] <- 1.0
    free[delta_idx] <- 0L
  }

  # correlation structure (new in 0.6-13)
  if (correlation) {
    var_idx <- which(lhs %in% ov_names &
      op == "~~" &
      user == 0L &
      lhs == rhs)
    ustart[var_idx] <- 1L
    free[var_idx] <- 0L

    delta_idx <- which(op == "~*~" &
      user == 0L)
    ustart[delta_idx] <- 1.0
    free[delta_idx] <- 0L
  }

  # group proportions (group 1L)
  if (group_w_free) {
    group_idx <- which(lhs == "group" & op == "%")
    # if(ngroups > 1L) {
    free[group_idx] <- 1L
    ustart[group_idx] <- as.numeric(NA)
    # } else {
    #      free[ group.idx ] <- 0L
    #    ustart[ group.idx ] <- 0.0 # last group
    # }
  }

  # 6. multiple groups?
  group <- rep(1L, length(lhs))
  if (ngroups > 1) {

    # only if "loadings" in group.equal and !std.lv:
    # construct temporary tmp.list to obtain lv.marker
    if (!std_lv && "loadings" %in% group_equal) {
      tmp_list <- list(
        id          = seq_along(lhs),
        lhs         = lhs,
        op          = op,
        rhs         = rhs,
        free        = free,
        ustart      = ustart,
        block       = rep(1, length(rhs)))
      lv_marker <- lav_pt_vnames(tmp_list, "lv.marker")
    }


    group <- rep(1:ngroups, each = length(lhs))
    user <- rep(user, times = ngroups)
    lhs <- rep(lhs, times = ngroups)
    op <- rep(op, times = ngroups)
    rhs <- rep(rhs, times = ngroups)
    free <- rep(free, times = ngroups)
    ustart <- rep(ustart, times = ngroups)
    mod_idx <- rep(mod_idx, times = ngroups)
    label <- rep(label, times = ngroups)
    exo <- rep(exo, times = ngroups)

    # specific changes per group
    for (g in 2:ngroups) {
      # free/fix intercepts latent variables
      if (meanstructure) {
        int_idx <- which(op == "~1" &
          lhs %in% lv_names_noc &
          user == 0L &
          group == g)
        if (int_lv_free == FALSE && g > 1 &&
          ("intercepts" %in% group_equal) &&
          !("means" %in% group_equal)) {
          free[int_idx] <- 1L
          ustart[int_idx] <- as.numeric(NA)
        }
      }

      # free intercept indicators if equal thresholds (new in 0.6-20)
      if (meanstructure && length(ov_names_ord) > 0L) {
        ord_idx <- which(op == "~1" &
          lhs %in% ov_names_ord &
          user == 0L &
          group == g)
        if (int_lv_free == FALSE && g > 1 &&
           "thresholds" %in% group_equal) {
          free[ord_idx] <- 1L
          ustart[ord_idx] <- as.numeric(NA)
        }
      }

      # latent variances if std.lv = TRUE (new in 0.6-4)
      if (std_lv && "loadings" %in% group_equal &&
        !"lv.variances" %in% group_equal) {
        lv_var_idx <- which(op == "~~" &
          lhs %in% lv_names &
          !lhs %in% lv_names_efa &
          lhs == rhs &
          user == 0L &
          group == g)
        if (length(lv_var_idx) > 0L) {
          free[lv_var_idx] <- 1L
          ustart[lv_var_idx] <- as.numeric(NA)
        }
      # marker indicator if std.lv = FALSE (new in 0.6-20)
      } else if (!std_lv && "loadings" %in% group_equal) {
        marker_idx <- which(op == "=~" &
          rhs %in% lv_marker &
          free == 0L &
          ustart == 1L &
          group == g)
        if (length(marker_idx) > 0L) {
          free[marker_idx] <- 1L
          ustart[marker_idx] <- as.numeric(NA)
        }
      }

      # latent variances if efa = TRUE (new in 0.6-5)
      if (length(lv_names_efa) > 0L &&
          auto_efa && "loadings" %in% group_equal &&
        !"lv.variances" %in% group_equal) {
        lv_var_idx <- which(op == "~~" &
          lhs %in% lv_names_efa &
          lhs == rhs &
          user == 0L &
          group == g)
        if (length(lv_var_idx) > 0L) {
          free[lv_var_idx] <- 1L
          ustart[lv_var_idx] <- as.numeric(NA)
        }
      }

      # latent response scaling -- categorical only
      # - if thresholds are equal -> free scalings/residual variances
      # - but not for binary indicators!
      if (length(ov_names_ord) > 0L) {
        nth <- sapply(ov_names_ord,
          function(x) sum(lhs == x & op == "|" & group == 1L))
        ov_names_ord_notbinary <- ov_names_ord[nth > 1L]
        if (auto_delta && parameterization == "delta") {
          if (any(op == "~*~" & group == g) &&
            ("thresholds" %in% group_equal)) {
            delta_idx <- which(op == "~*~" & group == g &
              lhs %in% ov_names_ord_notbinary)
            free[delta_idx] <- 1L
            ustart[delta_idx] <- as.numeric(NA)
          }
        } else if (parameterization == "theta") {
          if (any(op == "~*~" & group == g) &&
            ("thresholds" %in% group_equal)) {
            var_ord_idx <- which(op == "~~" & group == g &
              lhs %in% ov_names_ord_notbinary & lhs == rhs)
            free[var_ord_idx] <- 1L
            ustart[var_ord_idx] <- as.numeric(NA)
          }
        }
      }

      # group proportions
      if (group_w_free) {
        group_idx <- which(lhs == "group" & op == "%" & group == g)
        # if(g == ngroups) {
        #      free[ group.idx ] <- 0L
        #    ustart[ group.idx ] <- 0.0 # last group
        # } else {
        free[group_idx] <- 1L
        ustart[group_idx] <- as.numeric(NA)
        # }
      }
    } # g
  } # ngroups

  # construct tmp.list
  tmp_list <- list(
    id          = seq_along(lhs),
    lhs         = lhs,
    op          = op,
    rhs         = rhs,
    user        = user
  )

  # add block column (before group/level columns)
  if (!is.null(block_id)) {
    # only one block
    tmp_list$block <- rep(block_id, length(lhs))
  } else {
    # block is a combination of at least group, level, ...
    # for now, only group
    tmp_list$block <- group
  }

  # block columns (typically only group)
  for (block in blocks) {
    if (block == "group") {
      tmp_list[[block]] <- group
    } else {
      tmp_list[[block]] <- rep(0L, length(lhs))
    }
  }

  # other columns
  tmp_list2 <- list(
    mod.idx     = mod_idx,
    free        = free,
    ustart      = ustart,
    exo         = exo,
    label       = label
  )

  tmp_list <- c(tmp_list, tmp_list2)
}
