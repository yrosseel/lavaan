lav_partable_flat <- function(FLAT = NULL, # nolint
                              blocks = "group",
                              block.id = NULL,
                              meanstructure = FALSE,
                              int.ov.free = FALSE,
                              int.lv.free = FALSE,
                              orthogonal = FALSE,
                              orthogonal.y = FALSE,
                              orthogonal.x = FALSE,
                              orthogonal.efa = FALSE,
                              std.lv = FALSE,
                              correlation = FALSE,
                              conditional.x = FALSE,
                              fixed.x = TRUE,
                              parameterization = "delta",
                              auto.fix.first = FALSE,
                              auto.fix.single = FALSE,
                              auto.var = FALSE,
                              auto.cov.lv.x = FALSE,
                              auto.cov.y = FALSE,
                              auto.th = FALSE,
                              auto.delta = FALSE,
                              auto.efa = FALSE,
                              varTable = NULL, # nolint
                              group.equal = NULL,
                              group.w.free = FALSE,
                              ngroups = 1L,
                              nthresholds = NULL,
                              ov.names.x.block = NULL) {
  categorical <- FALSE

  ### tmp.default elements: parameters that are typically not specified by
  ###                   users, but should typically be considered,
  ###                   either free or fixed

  # extract `names' of various types of variables:
  lv.names <- lav_partable_vnames(FLAT, type = "lv") # latent variables
  # lv.names.r   <- lav_partable_vnames(FLAT, type="lv.regular")
  # regular latent variables
  lv.names.f <- lav_partable_vnames(FLAT, type = "lv.formative")
  # formative latent variables
  ov.names <- lav_partable_vnames(FLAT, type = "ov")
  # observed variables
  ov.names.x <- lav_partable_vnames(FLAT, type = "ov.x")
  # exogenous x covariates
  lv.names.int <- lav_partable_vnames(FLAT, type = "lv.interaction")
  # lv interactions

  if (is.null(ov.names.x.block)) {
    ov.names.x.block <- ov.names.x
  }
  ov.names.nox <- lav_partable_vnames(FLAT, type = "ov.nox")
  lv.names.x <- lav_partable_vnames(FLAT, type = "lv.x") # exogenous lv
  ov.names.y <- lav_partable_vnames(FLAT, type = "ov.y") # dependent ov
  lv.names.y <- lav_partable_vnames(FLAT, type = "lv.y") # dependent lv
  lv.names.efa <- lav_partable_vnames(FLAT, type = "lv.efa")
  # lvov.names.y <- c(ov.names.y, lv.names.y)
  lvov.names.y <- c(lv.names.y, ov.names.y)

  # get 'ordered' variables, either from FLAT or varTable
  ov.names.ord1 <- lav_partable_vnames(FLAT, type = "ov.ord")
  # check if we have "|" for exogenous variables
  if (length(ov.names.ord1) > 0L) {
    idx <- which(ov.names.ord1 %in% ov.names.x)
    if (length(idx) > 0L) {
      lav_msg_warn(gettext("thresholds are defined for exogenous variables:"),
                   lav_msg_view(ov.names.ord1[idx], "none"))
    }
  }

  # check data
  if (!is.null(varTable)) {
    ov.names.ord2 <-
      as.character(varTable$name[varTable$type == "ordered"])
    # remove fixed.x variables
    idx <- which(ov.names.ord2 %in% ov.names.x)
    if (length(idx) > 0L) {
      ov.names.ord2 <- ov.names.ord2[-idx]
    }

    # remove those that do appear in the model syntax
    idx <- which(!ov.names.ord2 %in% ov.names)
    if (length(idx) > 0L) {
      ov.names.ord2 <- ov.names.ord2[-idx]
    }
  } else {
    ov.names.ord2 <- character(0L)
  }

  # check nthresholds, if it is a named vector
  ov.names.ord3 <- character(0L)
  if (!is.null(nthresholds)) {
    if (!is.null(varTable)) {
      lav_msg_stop(gettext(
        "the varTable and nthresholds arguments should not be used together."))
    }
    if (!is.numeric(nthresholds)) {
      lav_msg_stop(gettext("nthresholds should be a named vector of integers."))
    }
    nth.names <- names(nthresholds)
    if (!is.null(nth.names)) {
      ov.names.ord3 <- nth.names
    } else {
      # if nthresholds is just a number, all is good; otherwise it
      # should be a names vector
      if (length(nthresholds) > 1L) {
        lav_msg_warn(gettext("nthresholds must be a named vector of integers."))
      }
      # just a single number -> assume ALL y variables are ordered
      ov.names.ord3 <- ov.names.nox
    }
  }

  # final ov.names.ord
  tmp <- unique(c(ov.names.ord1, ov.names.ord2, ov.names.ord3))
  ov.names.ord <- ov.names[ov.names %in% tmp]

  # if we have the "|" in the model syntax, check the number of thresholds
  # if(!is.null(varTable) && length(ov.names.ord1) > 0L) {
  #    for(o in ov.names.ord1) {
  #        nth <- varTable$nlev[ varTable$name == o ] - 1L
  #        nth.in.partable <- sum(FLAT$op == "|" & FLAT$lhs == o)
  #        if(nth != nth.in.partable) {
  #            stop("lavaan ERROR: expected ", max(0,nth),
  #                 " threshold(s) for variable ",
  #                 sQuote(o), "; syntax contains ", nth.in.partable, "\n")
  #        }
  #    }
  # }

  if (length(ov.names.ord) > 0L) {
    categorical <- TRUE
  }

  # std.lv = TRUE, group.equal includes "loadings"
  # if(ngroups > 1L && std.lv && "loadings" %in% group.equal) {
  # suggested by Michael Hallquist
  # in 0.6.3, we gave a warning,
  # warning("lavaan WARNING: std.lv = TRUE forces all variances to be unity",
  # " in all groups, despite group.equal = \"loadings\"")
  # in >0.6.4, we free the lv variances in all but the first group,
  # }


  # do we have any EFA lv's? they need special treatment if auto.efa = TRUE
  if (!is.null(FLAT$efa) && auto.efa) {
    lv.names.efa <- unique(FLAT$lhs[FLAT$op == "=~" &
      nchar(FLAT$efa) > 0L])
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
    lv.names.efa <- character(0)
  }

  lhs <- rhs <- character(0)

  # 1. THRESHOLDS (based on varTable)
  #    NOTE: - new in 0.5-18: ALWAYS include threshold parameters in partable,
  #            but only free them if auto.th = TRUE
  #          - [only ov.names.ord2, because ov.names.ord1 are already
  #            in tmp.user and we only need to add 'default' parameters here]
  #            (not any longer: we create them for ALL ordered var (0.6-12)
  nth <- 0L
  # if(auto.th && length(ov.names.ord2) > 0L) {
  # if(length(ov.names.ord2) > 0L) {
  if (length(ov.names.ord) > 0L) {
    # for(o in ov.names.ord2) {
    for (o in ov.names.ord) {
      if (!is.null(varTable)) {
        nth <- varTable$nlev[varTable$name == o] - 1L
      } else if (!is.null(nthresholds)) {
        if (length(nthresholds) == 1L && is.null(nth.names)) {
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
  ov.var <- ov.names.nox
  # auto-remove ordinal variables
  # idx <- match(ov.names.ord, ov.var)
  # if(length(idx)) ov.var <- ov.var[-idx]
  lhs <- c(lhs, ov.var, lv.names)
  rhs <- c(rhs, ov.var, lv.names)
  # }

  # b) `independent` latent variable COVARIANCES (lv.names.x)
  if (auto.cov.lv.x && length(lv.names.x) > 1L) {
    tmp <- utils::combn(lv.names.x, 2)
    lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
    rhs <- c(rhs, tmp[2, ])
  }

  # c) `dependent` latent variables COVARIANCES (lv.y.idx + ov.y.lv.idx)
  if (auto.cov.y && length(lvov.names.y) > 1L) {
    tmp <- utils::combn(lvov.names.y, 2L)
    lhs <- c(lhs, tmp[1, ]) # to fill upper.tri
    rhs <- c(rhs, tmp[2, ])
  }

  # d) exogenous x covariates: VARIANCES + COVARIANCES
  if ((nx <- length(ov.names.x)) > 0L) {
    if (conditional.x) {
      # new in 0.6-12: we make a distinction between ov.names.x and
      # ov.names.x.block: we treat them 'separately' (with no covariances
      # among them)
      # but we add 'regressions' instead (see below)
      ov.names.x1 <- ov.names.x[!ov.names.x %in% ov.names.x.block]
      ov.names.x2 <- ov.names.x.block
      nx1 <- length(ov.names.x1) # splitted x
      nx2 <- length(ov.names.x2) # regular  x
      if (nx1 > 0L) {
        idx <- lower.tri(matrix(0, nx1, nx1), diag = TRUE)
        lhs <- c(lhs, rep(ov.names.x1, each = nx1)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov.names.x1, times = nx1)[idx])
      }
      if (nx2 > 0L) {
        idx <- lower.tri(matrix(0, nx2, nx2), diag = TRUE)
        lhs <- c(lhs, rep(ov.names.x2, each = nx2)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov.names.x2, times = nx2)[idx])
      }
    } else {
      idx <- lower.tri(matrix(0, nx, nx), diag = TRUE)
      lhs <- c(lhs, rep(ov.names.x, each = nx)[idx]) # fill upper.tri
      rhs <- c(rhs, rep(ov.names.x, times = nx)[idx])
    }
  }

  # e) efa latent variables COVARIANCES; only needed for 'mediators'
  #    (not in lv.names.x, not in lv.names.y) -- added in 0.6-18
  if (auto.efa && length(lv.names.efa) > 1L) {
    efa.values <- lav_partable_efa_values(FLAT)
    for (set in efa.values) {
      # correlated factors within each set
      this.set.lv <- unique(FLAT$lhs[FLAT$op == "=~" &
        !FLAT$lhs %in% lv.names.x &
        !FLAT$lhs %in% lv.names.y &
        FLAT$efa == set])
      if (length(this.set.lv) > 0L) {
        tmp <- utils::combn(this.set.lv, 2)
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
  if (length(ov.names.ord) > 0L) {
    lhs <- c(lhs, ov.names.ord)
    rhs <- c(rhs, ov.names.ord)
    op <- c(op, rep("~*~", length(ov.names.ord)))
  }

  # same for correlation structures, but now for ALL variables
  if (!categorical && correlation) {
    lhs <- c(lhs, ov.names)
    rhs <- c(rhs, ov.names)
    op <- c(op, rep("~*~", length(ov.names)))
  }

  # 3. INTERCEPTS
  if (meanstructure) {
    # if(conditional.x) {
    #    ov.int <- ov.names.nox
    # } else {
    ov.int <- ov.names
    # }
    # auto-remove ordinal variables
    # idx <- which(ov.int %in% ov.names.ord)
    # if(length(idx)) ov.int <- ov.int[-idx]

    int.lhs <- c(ov.int, lv.names)
    lhs <- c(lhs, int.lhs)
    rhs <- c(rhs, rep("", length(int.lhs)))
    op <- c(op, rep("~1", length(int.lhs)))
  }

  # 4. REGRESSIONS
  if (conditional.x) {
    # new in 0.6-12: we make a distinction between ov.names.x and
    # ov.names.x.block: we treat them 'separately' (with no covariances
    # among them)
    # but we add 'regressions' instead!
    ov.names.x1 <- ov.names.x[!ov.names.x %in% ov.names.x.block]
    ov.names.x2 <- ov.names.x.block
    nx1 <- length(ov.names.x1) # splitted x
    nx2 <- length(ov.names.x2) # regular  x
    if (nx1 > 0L && nx2 > 0L) {
      # add regressions for splitted-x ~ regular-x
      lhs <- c(lhs, rep(ov.names.x1, times = nx2))
      op <- c(op, rep("~", nx2 * nx1))
      rhs <- c(rhs, rep(ov.names.x2, each = nx1))
    }
  }

  # free group weights
  if (group.w.free) {
    lhs <- c(lhs, "group")
    rhs <- c(rhs, "w")
    op <- c(op, "%")
  }

  tmp.default <- data.frame(
    lhs = lhs, op = op, rhs = rhs,
    mod.idx = rep(0L, length(lhs)),
    stringsAsFactors = FALSE
  )


  # 4. USER: user-specified elements
  lhs <- FLAT$lhs
  op <- FLAT$op
  rhs <- FLAT$rhs
  mod.idx <- FLAT$mod.idx

  lv.names <- lav_partable_vnames(FLAT, type = "lv") # latent variables
  ov.names <- lav_partable_vnames(FLAT, type = "ov") # observed variables
  tmp.user <- data.frame(
    lhs = lhs, op = op, rhs = rhs, mod.idx = mod.idx,
    stringsAsFactors = FALSE
  )

  # check for duplicated elements in tmp.user
  tmp.tmp <- tmp.user[, 1:3]
  idx <- which(duplicated(tmp.tmp))
  if (length(idx) > 0L) {
    txt <- sapply(seq_along(idx), function(i) {
      paste(
        "    ", tmp.tmp[idx[i], "lhs"],
        tmp.tmp[idx[i], "op"],
        tmp.tmp[idx[i], "rhs"]
      )
    })
    lav_msg_warn(gettext(
      "duplicated elements in model syntax have been ignored:"),
      lav_msg_view(txt, "none"))
    tmp.user <- tmp.user[-idx, ]
  }

  # check for duplicated elements in tmp.default
  # - FIXME: can we not avoid this somehow??
  # - for example, if the user model includes 'x1 ~~ x1'
  #   or 'x1 ~ 1'
  # - remove them from tmp.default
  tmp.tmp <- rbind(tmp.default[, 1:3], tmp.user[, 1:3])
  idx <- which(duplicated(tmp.tmp, fromLast = TRUE))
  # idx should be in tmp.default
  if (length(idx)) {
    for (i in idx) {
      flat.idx <- which(tmp.user$lhs == tmp.default$lhs[i] &
        tmp.user$op == tmp.default$op[i] &
        tmp.user$rhs == tmp.default$rhs[i])
      if (length(flat.idx) != 1L) {
        cat("[lavaan DEBUG] idx in tmp.tmp: i = ", i, "\n")
        print(tmp.tmp[i, ])
        cat("[lavaan DEBUG] idx in tmp.default: i = ", i, "\n")
        print(tmp.default[i, ])
        cat("[lavaan DEBUG] flat.idx:")
        print(flat.idx)
      }
    }
    tmp.default <- tmp.default[-idx, ]
  }

  # now that we have removed all duplicated elements, we can construct
  # the tmp.list for a single group/block
  lhs <- c(tmp.user$lhs, tmp.default$lhs)
  op <- c(tmp.user$op, tmp.default$op)
  rhs <- c(tmp.user$rhs, tmp.default$rhs)
  user <- c(
    rep(1L, length(tmp.user$lhs)),
    rep(0L, length(tmp.default$lhs))
  )
  mod.idx <- c(tmp.user$mod.idx, tmp.default$mod.idx)
  free <- rep(1L, length(lhs))
  ustart <- rep(as.numeric(NA), length(lhs))
  # label   <- paste(lhs, op, rhs, sep="")
  label <- rep(character(1), length(lhs))
  exo <- rep(0L, length(lhs))

  # 0a. if auto.th = FALSE, set fix the thresholds
  if (!auto.th) {
    th.idx <- which(op == "|" & user == 0L)
    free[th.idx] <- 0L
  }

  # 0b. if auto.var = FALSE, set the unspecified variances to zero
  if (!auto.var) {
    var.idx <- which(op == "~~" &
      lhs == rhs &
      user == 0L)
    ustart[var.idx] <- 0.0
    free[var.idx] <- 0L
  } else if (length(lv.names.f) > 0L) {
    # 'formative' (residual) variances are set to zero by default
    var.idx <- which(op == "~~" &
      lhs == rhs &
      lhs %in% lv.names.f &
      user == 0L)
    ustart[var.idx] <- 0.0
    free[var.idx] <- 0L
  }


  # 1. fix metric of regular latent variables
  if (std.lv) {
    # fix metric by fixing the variance of the latent variable
    lv.var.idx <- which(op == "~~" &
      lhs %in% lv.names & lhs == rhs)
    ustart[lv.var.idx] <- 1.0
    free[lv.var.idx] <- 0L
  }
  if (auto.efa && length(lv.names.efa) > 0L) {
    # fix lv variances of efa blocks to unity
    lv.var.idx <- which(op == "~~" &
      lhs %in% lv.names.efa & lhs == rhs)
    ustart[lv.var.idx] <- 1.0
    free[lv.var.idx] <- 0L
  }
  if (auto.fix.first) {
    # fix metric by fixing the loading of the first indicator
    # (but not for efa factors)
    mm.idx <- which(op == "=~" & !(lhs %in% lv.names.efa))
    first.idx <- mm.idx[which(!duplicated(lhs[mm.idx]))]
    ustart[first.idx] <- 1.0
    free[first.idx] <- 0L
  }

  # 2. fix residual variance of single indicators to zero
  if (auto.var && auto.fix.single) {
    mm.idx <- which(op == "=~")
    tmp.t <- table(lhs[mm.idx])
    if (any(tmp.t == 1L)) {
      # ok, we have a LV with only a single indicator
      lv.names.single <- names(tmp.t)[tmp.t == 1L]
      # get corresponding indicator if unique
      lhs.mm <- lhs[mm.idx]
      rhs.mm <- rhs[mm.idx]
      single.ind <- rhs.mm[which(lhs.mm %in% lv.names.single &
        lhs.mm != rhs.mm & # exclude phantom
        !(duplicated(rhs.mm) |
          duplicated(rhs.mm, fromLast = TRUE)))]
      # is the indicator unique?
      if (length(single.ind) > 0L) {
        var.idx <- which(op == "~~" & lhs %in% single.ind &
          rhs %in% single.ind &
          lhs == rhs &
          user == 0L)
        ustart[var.idx] <- 0.0
        free[var.idx] <- 0L
      }
    }
  }

  # 3. orthogonal = TRUE?
  if (orthogonal) {
    lv.cov.idx <- which(op == "~~" &
      lhs %in% lv.names &
      rhs %in% lv.names &
      lhs != rhs &
      user == 0L)
    ustart[lv.cov.idx] <- 0.0
    free[lv.cov.idx] <- 0L
  }
  # 3b. orthogonal.y = TRUE?
  if (orthogonal.y) {
    lv.cov.idx <- which(op == "~~" &
      lhs %in% lv.names.y &
      rhs %in% lv.names.y &
      lhs != rhs &
      user == 0L)
    ustart[lv.cov.idx] <- 0.0
    free[lv.cov.idx] <- 0L
  }
  # 3c. orthogonal.x = TRUE?
  if (orthogonal.x) {
    lv.cov.idx <- which(op == "~~" &
      lhs %in% lv.names.x &
      rhs %in% lv.names.x &
      lhs != rhs &
      user == 0L)
    ustart[lv.cov.idx] <- 0.0
    free[lv.cov.idx] <- 0L
  }
  # 3d. orthogonal.efa = TRUE?
  if (orthogonal.efa) {
    lv.cov.idx <- which(op == "~~" &
      lhs %in% lv.names.efa &
      rhs %in% lv.names.efa &
      lhs != rhs &
      user == 0L)
    ustart[lv.cov.idx] <- 0.0
    free[lv.cov.idx] <- 0L
  }

  # 4. intercepts
  if (meanstructure) {
    if (categorical) {
      # zero intercepts/means ordinal variables
      ov.int.idx <- which(op == "~1" &
        lhs %in% ov.names.ord &
        user == 0L)
      ustart[ov.int.idx] <- 0.0
      free[ov.int.idx] <- 0L
    }
    if (int.ov.free == FALSE) {
      # zero intercepts/means observed variables
      ov.int.idx <- which(op == "~1" &
        lhs %in% ov.names &
        user == 0L)
      ustart[ov.int.idx] <- 0.0
      free[ov.int.idx] <- 0L
    }
    if (int.lv.free == FALSE) {
      # zero intercepts/means latent variables
      lv.int.idx <- which(op == "~1" &
        lhs %in% lv.names &
        user == 0L)
      ustart[lv.int.idx] <- 0.0
      free[lv.int.idx] <- 0L
    }
    # 4b. fixed effect (only if we have random slopes)
    if (!is.null(FLAT$rv) && any(nchar(FLAT$rv) > 0L)) {
      lv.names.rv <- lav_partable_vnames(FLAT, "lv.rv")
      lv.rv.idx <- which(op == "~1" &
        lhs %in% lv.names.rv &
        user == 0L)
      ustart[lv.rv.idx] <- as.numeric(NA)
      free[lv.rv.idx] <- 1L
    }
    if (length(lv.names.int) > 0L) {
      lv.int.idx <- which(op == "~1" &
        lhs %in% lv.names.int &
        user == 0L)
      ustart[lv.int.idx] <- as.numeric(NA)
      free[lv.int.idx] <- 1L
    }
  }

  # 4b. fixed effect (only if we have random slopes)
  # if(!is.null(FLAT$rv)) {
  #    }

  # 5. handle exogenous `x' covariates
  # usually, ov.names.x.block == ov.names.x
  # except if multilevel, where 'splitted' ov.x are treated as endogenous

  # 5a conditional.x = FALSE
  if (!conditional.x && fixed.x && length(ov.names.x.block) > 0) {
    # 1. variances/covariances
    exo.var.idx <- which(op == "~~" &
      rhs %in% ov.names.x.block &
      lhs %in% ov.names.x.block &
      user == 0L)
    ustart[exo.var.idx] <- as.numeric(NA) # should be overriden later!
    free[exo.var.idx] <- 0L
    exo[exo.var.idx] <- 1L

    # 2. intercepts
    exo.int.idx <- which(op == "~1" &
      lhs %in% ov.names.x.block &
      user == 0L)
    ustart[exo.int.idx] <- as.numeric(NA) # should be overriden later!
    free[exo.int.idx] <- 0L
    exo[exo.int.idx] <- 1L
  }

  # 5a-bis. conditional.x = TRUE
  if (conditional.x && length(ov.names.x) > 0L) {
    # 1. variances/covariances
    exo.var.idx <- which(op == "~~" &
      rhs %in% ov.names.x &
      lhs %in% ov.names.x &
      user == 0L)
    if (fixed.x) {
      ustart[exo.var.idx] <- as.numeric(NA) # should be overriden later!
      free[exo.var.idx] <- 0L
    }
    exo[exo.var.idx] <- 1L

    # 2. intercepts
    exo.int.idx <- which(op == "~1" &
      lhs %in% ov.names.x &
      user == 0L)
    if (fixed.x) {
      ustart[exo.int.idx] <- as.numeric(NA) # should be overriden later!
      free[exo.int.idx] <- 0L
    }
    exo[exo.int.idx] <- 1L

    # 3. regressions ov + lv
    exo.reg.idx <- which(op %in% c("~", "<~") &
      lhs %in% c(lv.names, ov.names.nox) &
      rhs %in% ov.names.x)
    exo[exo.reg.idx] <- 1L

    # 3b regression splitted.x ~ regular.x
    exo.reg2.idx <- which(op %in% c("~", "<~") &
      lhs %in% ov.names.x &
      rhs %in% ov.names.x)
    if (fixed.x) {
      ustart[exo.reg2.idx] <- as.numeric(NA) # should be overriden later!
      free[exo.reg2.idx] <- 0L
    }
    exo[exo.reg2.idx] <- 1L
  }

  # 5b. residual variances of ordinal variables?
  if (length(ov.names.ord) > 0L) {
    ord.idx <- which(lhs %in% ov.names.ord &
      op == "~~" &
      user == 0L & ## New in 0.6-1
      lhs == rhs)
    ustart[ord.idx] <- 1L ## FIXME!! or 0?? (0 breaks ex3.12)
    free[ord.idx] <- 0L
  }

  # 5c latent response scales of ordinal variables?
  #    by default, all fixed to 1.0
  if (length(ov.names.ord) > 0L) {
    delta.idx <- which(op == "~*~" &
      user == 0L) ## New in 0.6-1
    ustart[delta.idx] <- 1.0
    free[delta.idx] <- 0L
  }

  # correlation structure (new in 0.6-13)
  if (correlation) {
    var.idx <- which(lhs %in% ov.names &
      op == "~~" &
      user == 0L &
      lhs == rhs)
    ustart[var.idx] <- 1L
    free[var.idx] <- 0L

    delta.idx <- which(op == "~*~" &
      user == 0L)
    ustart[delta.idx] <- 1.0
    free[delta.idx] <- 0L
  }

  # group proportions (group 1L)
  if (group.w.free) {
    group.idx <- which(lhs == "group" & op == "%")
    # if(ngroups > 1L) {
    free[group.idx] <- 1L
    ustart[group.idx] <- as.numeric(NA)
    # } else {
    #      free[ group.idx ] <- 0L
    #    ustart[ group.idx ] <- 0.0 # last group
    # }
  }

  # 6. multiple groups?
  group <- rep(1L, length(lhs))
  if (ngroups > 1) {
    group <- rep(1:ngroups, each = length(lhs))
    user <- rep(user, times = ngroups)
    lhs <- rep(lhs, times = ngroups)
    op <- rep(op, times = ngroups)
    rhs <- rep(rhs, times = ngroups)
    free <- rep(free, times = ngroups)
    ustart <- rep(ustart, times = ngroups)
    mod.idx <- rep(mod.idx, times = ngroups)
    label <- rep(label, times = ngroups)
    exo <- rep(exo, times = ngroups)

    # specific changes per group
    for (g in 2:ngroups) {
      # label
      # label[group == g] <- paste(label[group == 1], ".g", g, sep="")

      # free/fix intercepts
      if (meanstructure) {
        int.idx <- which(op == "~1" &
          lhs %in% lv.names &
          user == 0L &
          group == g)
        if (int.lv.free == FALSE && g > 1 &&
          ("intercepts" %in% group.equal ||
            "thresholds" %in% group.equal) &&
          !("means" %in% group.equal)) {
          free[int.idx] <- 1L
          ustart[int.idx] <- as.numeric(NA)
        }
      }

      # latent variances if std.lv = TRUE (new in 0.6-4)
      if (std.lv && "loadings" %in% group.equal &&
        !"lv.variances" %in% group.equal) {
        lv.var.idx <- which(op == "~~" &
          lhs %in% lv.names &
          !lhs %in% lv.names.efa &
          lhs == rhs &
          user == 0L &
          group == g)
        if (length(lv.var.idx) > 0L) {
          free[lv.var.idx] <- 1L
          ustart[lv.var.idx] <- as.numeric(NA)
        }
      }

      # latent variances if efa = TRUE (new in 0.6-5)
      if (auto.efa && "loadings" %in% group.equal &&
        !"lv.variances" %in% group.equal) {
        lv.var.idx <- which(op == "~~" &
          lhs %in% lv.names.efa &
          lhs == rhs &
          user == 0L &
          group == g)
        if (length(lv.var.idx) > 0L) {
          free[lv.var.idx] <- 1L
          ustart[lv.var.idx] <- as.numeric(NA)
        }
      }

      # latent response scaling
      if (auto.delta && parameterization == "delta") {
        if (any(op == "~*~" & group == g) &&
          ("thresholds" %in% group.equal)) {
          delta.idx <- which(op == "~*~" & group == g)
          free[delta.idx] <- 1L
          ustart[delta.idx] <- as.numeric(NA)
        }
      } else if (parameterization == "theta") {
        if (any(op == "~*~" & group == g) &&
          ("thresholds" %in% group.equal)) {
          var.ord.idx <- which(op == "~~" & group == g &
            lhs %in% ov.names.ord & lhs == rhs)
          free[var.ord.idx] <- 1L
          ustart[var.ord.idx] <- as.numeric(NA)
        }
      }

      # group proportions
      if (group.w.free) {
        group.idx <- which(lhs == "group" & op == "%" & group == g)
        # if(g == ngroups) {
        #      free[ group.idx ] <- 0L
        #    ustart[ group.idx ] <- 0.0 # last group
        # } else {
        free[group.idx] <- 1L
        ustart[group.idx] <- as.numeric(NA)
        # }
      }
    } # g
  } # ngroups

  # construct tmp.list
  tmp.list <- list(
    id          = seq_along(lhs),
    lhs         = lhs,
    op          = op,
    rhs         = rhs,
    user        = user
  )

  # add block column (before group/level columns)
  if (!is.null(block.id)) {
    # only one block
    tmp.list$block <- rep(block.id, length(lhs))
  } else {
    # block is a combination of at least group, level, ...
    # for now, only group
    tmp.list$block <- group
  }

  # block columns (typically only group)
  for (block in blocks) {
    if (block == "group") {
      tmp.list[[block]] <- group
    } else {
      tmp.list[[block]] <- rep(0L, length(lhs))
    }
  }

  # other columns
  tmp.list2 <- list(
    mod.idx     = mod.idx,
    free        = free,
    ustart      = ustart,
    exo         = exo,
    label       = label
  )

  tmp.list <- c(tmp.list, tmp.list2)
}
