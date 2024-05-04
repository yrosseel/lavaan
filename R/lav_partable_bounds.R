# add parameter bounds to the parameter table
# lavoptions$optim.bounds
lav_partable_add_bounds <- function(partable = NULL,
                                    lavh1 = NULL,
                                    lavdata = NULL,
                                    lavsamplestats = NULL,
                                    lavoptions = NULL) {
  # no support (yet) for multilevel
  if (lav_partable_nlevels(partable) > 1L) {
    return(partable)
  }

  # check optim.bounds
  if (is.null(lavoptions$optim.bounds)) {
    # <0.6-6 version
    return(partable)
  } else if (!is.null(lavoptions$samplestats) && !lavoptions$samplestats) {
    # no sample statistics
    return(partable)
  } else {
    if (!is.null(lavoptions$bounds) && lavoptions$bounds == "none") {
      # no bounds needed
      return(partable)
    }

    # no support from effect.coding (for now)
    if (!is.null(lavoptions$effect.coding) &&
      nchar(lavoptions$effect.coding[1L]) > 0L) {
      lav_msg_warn(gettext(
        "automatic bounds not available (yet) if effect.coding is used"
        ))
      return(partable)
    }

    optim.bounds <- lavoptions$optim.bounds

    # check the elements
    if (is.null(optim.bounds$lower)) {
      optim.bounds$lower <- character(0L)
    } else {
      optim.bounds$lower <- as.character(optim.bounds$lower)
    }
    if (is.null(optim.bounds$upper)) {
      optim.bounds$upper <- character(0L)
    } else {
      optim.bounds$upper <- as.character(optim.bounds$upper)
    }

    if (is.null(optim.bounds$min.reliability.marker)) {
      optim.bounds$min.reliability.marker <- 0.0
    } else {
      if (optim.bounds$min.reliability.marker < 0 ||
        optim.bounds$min.reliability.marker > 1.0) {
        lav_msg_stop(gettextf(
          "optim.bounds$min.reliability.marker is out of range: %s",
          optim.bounds$min.reliability.marker
        ))
      }
    }

    if (is.null(optim.bounds$min.var.ov)) {
      optim.bounds$min.var.ov <- -Inf
    }

    if (is.null(optim.bounds$min.var.lv.exo)) {
      optim.bounds$min.var.lv.exo <- 0.0
    }

    if (is.null(optim.bounds$min.var.lv.endo)) {
      optim.bounds$min.var.lv.endo <- 0.0
    }

    if (is.null(optim.bounds$max.r2.lv.endo)) {
      optim.bounds$max.r2.lv.endo <- 1.0
    }

    if (is.null(optim.bounds$lower.factor)) {
      optim.bounds$lower.factor <- rep(1.0, length(optim.bounds$lower))
    } else {
      if (length(optim.bounds$lower.factor) == 1L &&
        is.numeric(optim.bounds$lower.factor)) {
        optim.bounds$lower.factor <- rep(
          optim.bounds$lower.factor,
          length(optim.bounds$lower)
        )
      } else if (length(optim.bounds$lower.factor) !=
        length(optim.bounds$lower)) {
        lav_msg_stop(
          gettext("length(optim.bounds$lower.factor) is not equal to
                  length(optim.bounds$lower)")
        )
      }
    }
    lower.factor <- optim.bounds$lower.factor

    if (is.null(optim.bounds$upper.factor)) {
      optim.bounds$upper.factor <- rep(1.0, length(optim.bounds$upper))
    } else {
      if (length(optim.bounds$upper.factor) == 1L &&
        is.numeric(optim.bounds$upper.factor)) {
        optim.bounds$upper.factor <- rep(
          optim.bounds$upper.factor,
          length(optim.bounds$upper)
        )
      } else if (length(optim.bounds$upper.factor) !=
        length(optim.bounds$upper)) {
        lav_msg_stop(
          gettext("length(optim.bounds$lower.factor) is not equal to
                  length(optim.bounds$upper)")
        )
      }
    }
    upper.factor <- optim.bounds$upper.factor
  }

  # new in 0.6-17: check if we have theta parameterization
  theta.parameterization.flag <- FALSE
  if (any(partable$op == "~*~") && lavoptions$parameterization == "theta") {
    # some fixed-to-1 theta elements?
    ov.scaled <- partable$lhs[partable$op == "~*~"]
    ov.var.idx <- which(partable$op == "~~" &
      partable$lhs %in% ov.scaled &
      partable$free == 0L &
      partable$ustart == 1)
    if (length(ov.var.idx) > 0L) {
      theta.parameterization.flag <- TRUE
      theta.parameterization.names <- partable$lhs[ov.var.idx]
    }
  }

  # shortcut
  REL <- optim.bounds$min.reliability.marker

  # nothing to do
  if (length(optim.bounds$lower) == 0L &&
    length(optim.bounds$upper) == 0L) {
    return(partable)
  } else {
    # we compute ALL bounds, then we select what we need
    # (otherwise, we can not use the 'factor')

    if (!is.null(partable$lower)) {
      lower.user <- partable$lower
    } else {
      partable$lower <- lower.user <- rep(-Inf, length(partable$lhs))
    }
    if (!is.null(partable$upper)) {
      upper.user <- partable$upper
    } else {
      partable$upper <- upper.user <- rep(+Inf, length(partable$lhs))
    }

    # the 'automatic' bounds
    lower.auto <- rep(-Inf, length(partable$lhs))
    upper.auto <- rep(+Inf, length(partable$lhs))
  }

  lavpta <- lav_partable_attributes(partable)

  # check blocks
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$lhs))
  }
  block.values <- lav_partable_block_values(partable)

  # check groups
  if (is.null(partable$group)) {
    partable$group <- rep(1L, length(partable$lhs))
  }
  group.values <- lav_partable_group_values(partable)
  ngroups <- length(group.values)

  # compute bounds per group ### TODO: add levels/classes/...
  b <- 0L
  for (g in seq_len(ngroups)) {
    # next block
    b <- b + 1L

    # for this block
    ov.names <- lavpta$vnames$ov[[b]]
    lv.names <- lavpta$vnames$lv[[b]]
    lv.names.x <- lavpta$vnames$lv.x[[b]]
    if (length(lv.names.x) > 0L) {
      lv.names.endo <- lv.names[!lv.names %in% lv.names.x]
    } else {
      lv.names.endo <- lv.names
    }
    lv.marker <- lavpta$vnames$lv.marker[[b]]

    # OV.VAR for this group
    if (lavsamplestats@missing.flag && lavdata@nlevels == 1L) {
      OV.VAR <- diag(lavsamplestats@missing.h1[[g]]$sigma)
    } else {
      if (lavoptions$conditional.x) {
        OV.VAR <- diag(lavsamplestats@res.cov[[g]])
      } else {
        OV.VAR <- diag(lavsamplestats@cov[[g]])
      }
    }

    # new in 0.6-17: increase observed variances for 'scaled' parameters
    # if theta parameterization
    if (theta.parameterization.flag) {
      sc.idx <- match(theta.parameterization.names, ov.names)
      OV.VAR[sc.idx] <- OV.VAR[sc.idx] / REL
    }


    # we 'process' the parameters per 'type', so we can choose
    # to apply (or not) upper/lower bounds for each type separately

    ################################
    ## 1. (residual) ov variances ##
    ################################
    par.idx <- which(partable$group == group.values[g] &
      partable$op == "~~" &
      partable$lhs %in% ov.names &
      partable$lhs == partable$rhs)

    if (length(par.idx) > 0L) {
      # lower == 0
      lower.auto[par.idx] <- 0

      # upper == var(ov)
      var.idx <- match(partable$lhs[par.idx], ov.names)
      upper.auto[par.idx] <- OV.VAR[var.idx]

      # if reliability > 0, adapt marker indicators only
      if (REL > 0) {
        marker.idx <- which(partable$group == group.values[g] &
          partable$op == "~~" &
          partable$lhs %in% lv.marker &
          partable$lhs == partable$rhs)
        marker.var.idx <- match(partable$lhs[marker.idx], ov.names)

        # upper = (1-REL)*OVAR
        upper.auto[marker.idx] <- (1 - REL) * OV.VAR[marker.var.idx]
      }

      # range
      bound.range <- upper.auto[par.idx] - pmax(lower.auto[par.idx], 0)

      # enlarge lower?
      if ("ov.var" %in% optim.bounds$lower) {
        factor <- lower.factor[which(optim.bounds$lower == "ov.var")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          diff <- abs(new.range - bound.range)
          lower.auto[par.idx] <- lower.auto[par.idx] - diff
        }
      }

      # enlarge upper?
      if ("ov.var" %in% optim.bounds$upper) {
        factor <- upper.factor[which(optim.bounds$upper == "ov.var")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          diff <- abs(new.range - bound.range)
          upper.auto[par.idx] <- upper.auto[par.idx] + diff
        }
      }

      # min.var.ov?
      min.idx <- which(lower.auto[par.idx] < optim.bounds$min.var.ov)
      if (length(min.idx) > 0L) {
        lower.auto[par.idx[min.idx]] <- optim.bounds$min.var.ov
      }

      # requested?
      if ("ov.var" %in% optim.bounds$lower) {
        partable$lower[par.idx] <- lower.auto[par.idx]
      }
      if ("ov.var" %in% optim.bounds$upper) {
        partable$upper[par.idx] <- upper.auto[par.idx]
      }
    } # (res) ov variances

    ################################
    ## 2. (residual) lv variances ##
    ################################

    # first collect lower/upper bounds for TOTAL variances in lv.names
    LV.VAR.LB <- numeric(length(lv.names))
    LV.VAR.UB <- numeric(length(lv.names))

    if (lavoptions$std.lv) {
      LV.VAR.LB <- rep(1.0, length(lv.names))
      LV.VAR.UB <- rep(1.0, length(lv.names))
    } else {
      for (i in seq_len(length(lv.names))) {
        this.lv.name <- lv.names[i]
        this.lv.marker <- lv.marker[i]

        if (nchar(this.lv.marker) > 0L && this.lv.marker %in% ov.names) {
          marker.var <- OV.VAR[match(this.lv.marker, ov.names)]
          LOWER <- marker.var - (1 - REL) * marker.var
          LV.VAR.LB[i] <- max(LOWER, optim.bounds$min.var.lv.exo)
          # LV.VAR.UB[i] <- marker.var - REL*marker.var
          LV.VAR.UB[i] <- marker.var

          # new in 0.6-17
          if (theta.parameterization.flag) {
            LV.VAR.LB[i] <- REL
          }
        } else {
          LV.VAR.LB[i] <- optim.bounds$min.var.lv.exo
          LV.VAR.UB[i] <- max(OV.VAR)
        }
      }
    }

    # use these bounds for the free parameters
    par.idx <- which(partable$group == group.values[g] &
      partable$op == "~~" &
      partable$lhs %in% lv.names &
      partable$lhs == partable$rhs)

    if (length(par.idx) > 0L) {
      # adjust for endogenenous lv
      LV.VAR.LB2 <- LV.VAR.LB
      endo.idx <- which(lv.names %in% lv.names.endo)
      if (length(endo.idx) > 0L) {
        LV.VAR.LB2[endo.idx] <- optim.bounds$min.var.lv.endo
        if (optim.bounds$max.r2.lv.endo != 1) {
          LV.VAR.LB2[endo.idx] <- (1 - optim.bounds$max.r2.lv.endo) * LV.VAR.UB[endo.idx]
        }
      }
      exo.idx <- which(!lv.names %in% lv.names.endo)
      if (length(exo.idx) > 0L && optim.bounds$min.var.lv.exo != 0) {
        LV.VAR.LB2[exo.idx] <- optim.bounds$min.var.lv.exo
      }

      lower.auto[par.idx] <- LV.VAR.LB2[match(
        partable$lhs[par.idx],
        lv.names
      )]
      upper.auto[par.idx] <- LV.VAR.UB[match(
        partable$lhs[par.idx],
        lv.names
      )]

      # range
      bound.range <- upper.auto[par.idx] - pmax(lower.auto[par.idx], 0)

      # enlarge lower?
      if ("lv.var" %in% optim.bounds$lower) {
        factor <- lower.factor[which(optim.bounds$lower == "lv.var")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          diff <- abs(new.range - bound.range)
          lower.auto[par.idx] <- lower.auto[par.idx] - diff
        }
      }

      # enlarge upper?
      if ("lv.var" %in% optim.bounds$upper) {
        factor <- upper.factor[which(optim.bounds$upper == "lv.var")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          diff <- abs(new.range - bound.range)
          upper.auto[par.idx] <- upper.auto[par.idx] + diff
        }
      }

      # requested?
      if ("lv.var" %in% optim.bounds$lower) {
        partable$lower[par.idx] <- lower.auto[par.idx]
      }
      if ("lv.var" %in% optim.bounds$upper) {
        partable$upper[par.idx] <- upper.auto[par.idx]
      }
    } # lv variances


    #############################################
    ## 3. factor loadings (ov indicators only) ##
    #############################################

    # lambda_p^(u) = sqrt( upper(res.var.indicators_p) /
    #                      lower(var.factor) )

    ov.ind.names <- lavpta$vnames$ov.ind[[b]]
    par.idx <- which(partable$group == group.values[g] &
      partable$op == "=~" &
      partable$lhs %in% lv.names &
      partable$rhs %in% ov.ind.names)

    if (length(par.idx) > 0L) {
      # if negative LV variances are allowed (due to factor > 1)
      # make them equal to zero
      LV.VAR.LB[LV.VAR.LB < 0] <- 0.0

      var.all <- OV.VAR[match(partable$rhs[par.idx], ov.names)]
      tmp <- LV.VAR.LB[match(partable$lhs[par.idx], lv.names)]
      tmp[is.na(tmp)] <- 0 # just in case...
      lower.auto[par.idx] <- -1 * sqrt(var.all / tmp) # -Inf if tmp==0
      upper.auto[par.idx] <- +1 * sqrt(var.all / tmp) # +Inf if tmp==0

      # if std.lv = TRUE, force 'first' loading to be positive?
      # if(lavoptions$std.lv) {
      #    # get index 'first' indicators
      #    first.idx <- which(!duplicated(partable$lhs[par.idx]))
      #    lower.auto[par.idx][first.idx] <- 0
      # }

      # range
      bound.range <- upper.auto[par.idx] - lower.auto[par.idx]

      # enlarge lower?
      if ("loadings" %in% optim.bounds$lower) {
        factor <- lower.factor[which(optim.bounds$lower == "loadings")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          ok.idx <- is.finite(new.range)
          if (length(ok.idx) > 0L) {
            diff <- abs(new.range[ok.idx] - bound.range[ok.idx])
            lower.auto[par.idx][ok.idx] <-
              lower.auto[par.idx][ok.idx] - diff
          }
        }
      }

      # enlarge upper?
      if ("loadings" %in% optim.bounds$upper) {
        factor <- upper.factor[which(optim.bounds$upper == "loadings")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          ok.idx <- is.finite(new.range)
          if (length(ok.idx) > 0L) {
            diff <- abs(new.range[ok.idx] - bound.range[ok.idx])
            upper.auto[par.idx][ok.idx] <-
              upper.auto[par.idx][ok.idx] + diff
          }
        }
      }



      # requested?
      if ("loadings" %in% optim.bounds$lower) {
        partable$lower[par.idx] <- lower.auto[par.idx]
      }
      if ("loadings" %in% optim.bounds$upper) {
        partable$upper[par.idx] <- upper.auto[par.idx]
      }
    } # lambda


    ####################
    ## 4. covariances ##
    ####################

    # | sqrt(var(x)) sqrt(var(y)) | <= cov(x,y)

    par.idx <- which(partable$group == group.values[g] &
      partable$op == "~~" &
      partable$lhs != partable$rhs)

    if (length(par.idx) > 0L) {
      for (i in seq_len(length(par.idx))) {
        # this lhs/rhs
        this.lhs <- partable$lhs[par.idx[i]]
        this.rhs <- partable$rhs[par.idx[i]]

        # 2 possibilities:
        # - variances are free parameters
        # - variances are fixed (eg std.lv = TRUE)

        # var idx
        lhs.var.idx <- which(partable$group == group.values[g] &
          partable$op == "~~" &
          partable$lhs == this.lhs &
          partable$lhs == partable$rhs)
        rhs.var.idx <- which(partable$group == group.values[g] &
          partable$op == "~~" &
          partable$lhs == this.rhs &
          partable$lhs == partable$rhs)
        # upper bounds
        lhs.upper <- upper.auto[lhs.var.idx]
        rhs.upper <- upper.auto[rhs.var.idx]

        # compute upper bounds for this cov (assuming >0 vars)
        if (is.finite(lhs.upper) && is.finite(rhs.upper)) {
          upper.cov <- sqrt(lhs.upper) * sqrt(rhs.upper)
          upper.auto[par.idx[i]] <- +1 * upper.cov
          lower.auto[par.idx[i]] <- -1 * upper.cov
        }
      }

      # range
      bound.range <- upper.auto[par.idx] - lower.auto[par.idx]

      # enlarge lower?
      if ("covariances" %in% optim.bounds$lower) {
        factor <-
          lower.factor[which(optim.bounds$lower == "covariances")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          ok.idx <- is.finite(new.range)
          if (length(ok.idx) > 0L) {
            diff <- new.range[ok.idx] - bound.range[ok.idx]
            lower.auto[par.idx][ok.idx] <-
              lower.auto[par.idx][ok.idx] - diff
          }
        }
      }

      # enlarge upper?
      if ("covariances" %in% optim.bounds$upper) {
        factor <-
          upper.factor[which(optim.bounds$upper == "covariances")]
        if (is.finite(factor) && factor != 1.0) {
          new.range <- bound.range * factor
          ok.idx <- is.finite(new.range)
          if (length(ok.idx) > 0L) {
            diff <- new.range[ok.idx] - bound.range[ok.idx]
            upper.auto[par.idx][ok.idx] <-
              upper.auto[par.idx][ok.idx] + diff
          }
        }
      }

      # requested?
      if ("covariances" %in% optim.bounds$lower) {
        partable$lower[par.idx] <- lower.auto[par.idx]
      }
      if ("covariances" %in% optim.bounds$upper) {
        partable$upper[par.idx] <- upper.auto[par.idx]
      }
    } # covariances
  } # g

  # overwrite with lower.user (except -Inf)
  not.inf.idx <- which(lower.user > -Inf)
  if (length(not.inf.idx) > 0L) {
    partable$lower[not.inf.idx] <- lower.user[not.inf.idx]
  }

  # overwrite with upper.user (except +Inf)
  not.inf.idx <- which(upper.user < +Inf)
  if (length(not.inf.idx) > 0L) {
    partable$upper[not.inf.idx] <- upper.user[not.inf.idx]
  }

  # non-free
  non.free.idx <- which(partable$free == 0L)
  if (length(non.free.idx) > 0L && !is.null(partable$ustart)) {
    partable$lower[non.free.idx] <- partable$ustart[non.free.idx]
    partable$upper[non.free.idx] <- partable$ustart[non.free.idx]
  }

  partable
}
