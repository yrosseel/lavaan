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

    optim_bounds <- lavoptions$optim.bounds

    # check the elements
    if (is.null(optim_bounds$lower)) {
      optim_bounds$lower <- character(0L)
    } else {
      optim_bounds$lower <- as.character(optim_bounds$lower)
    }
    if (is.null(optim_bounds$upper)) {
      optim_bounds$upper <- character(0L)
    } else {
      optim_bounds$upper <- as.character(optim_bounds$upper)
    }

    if (is.null(optim_bounds$min.reliability.marker)) {
      optim_bounds$min.reliability.marker <- 0.0
    } else {
      if (optim_bounds$min.reliability.marker < 0 ||
        optim_bounds$min.reliability.marker > 1.0) {
        lav_msg_stop(gettextf(
          "optim.bounds$min.reliability.marker is out of range: %s",
          optim_bounds$min.reliability.marker
        ))
      }
    }

    if (is.null(optim_bounds$min.var.ov)) {
      optim_bounds$min.var.ov <- -Inf
    }

    if (is.null(optim_bounds$min.var.lv.exo)) {
      optim_bounds$min.var.lv.exo <- 0.0
    }

    if (is.null(optim_bounds$min.var.lv.endo)) {
      optim_bounds$min.var.lv.endo <- 0.0
    }

    if (is.null(optim_bounds$max.r2.lv.endo)) {
      optim_bounds$max.r2.lv.endo <- 1.0
    }

    if (is.null(optim_bounds$lower.factor)) {
      optim_bounds$lower.factor <- rep(1.0, length(optim_bounds$lower))
    } else {
      if (length(optim_bounds$lower.factor) == 1L &&
        is.numeric(optim_bounds$lower.factor)) {
        optim_bounds$lower.factor <- rep(
          optim_bounds$lower.factor,
          length(optim_bounds$lower)
        )
      } else if (length(optim_bounds$lower.factor) !=
        length(optim_bounds$lower)) {
        lav_msg_stop(
          gettext("length(optim.bounds$lower.factor) is not equal to
                  length(optim.bounds$lower)")
        )
      }
    }
    lower_factor <- optim_bounds$lower.factor

    if (is.null(optim_bounds$upper.factor)) {
      optim_bounds$upper.factor <- rep(1.0, length(optim_bounds$upper))
    } else {
      if (length(optim_bounds$upper.factor) == 1L &&
        is.numeric(optim_bounds$upper.factor)) {
        optim_bounds$upper.factor <- rep(
          optim_bounds$upper.factor,
          length(optim_bounds$upper)
        )
      } else if (length(optim_bounds$upper.factor) !=
        length(optim_bounds$upper)) {
        lav_msg_stop(
          gettext("length(optim.bounds$lower.factor) is not equal to
                  length(optim.bounds$upper)")
        )
      }
    }
    upper_factor <- optim_bounds$upper.factor
  }

  # new in 0.6-17: check if we have theta parameterization
  theta_parameterization_flag <- FALSE
  if (any(partable$op == "~*~") && lavoptions$parameterization == "theta") {
    # some fixed-to-1 theta elements?
    ov_scaled <- partable$lhs[partable$op == "~*~"]
    ov_var_idx <- which(partable$op == "~~" &
      partable$lhs %in% ov_scaled &
      partable$free == 0L &
      partable$ustart == 1)
    if (length(ov_var_idx) > 0L) {
      theta_parameterization_flag <- TRUE
      theta_parameterization_names <- partable$lhs[ov_var_idx]
    }
  }

  # shortcut
  rel <- optim_bounds$min.reliability.marker

  # nothing to do
  if (length(optim_bounds$lower) == 0L &&
    length(optim_bounds$upper) == 0L) {
    return(partable)
  } else {
    # we compute ALL bounds, then we select what we need
    # (otherwise, we can not use the 'factor')

    if (!is.null(partable$lower)) {
      lower_user <- partable$lower
    } else {
      partable$lower <- lower_user <- rep(-Inf, length(partable$lhs))
    }
    if (!is.null(partable$upper)) {
      upper_user <- partable$upper
    } else {
      partable$upper <- upper_user <- rep(+Inf, length(partable$lhs))
    }

    # the 'automatic' bounds
    lower_auto <- rep(-Inf, length(partable$lhs))
    upper_auto <- rep(+Inf, length(partable$lhs))
  }

  lavpta <- lav_partable_attributes(partable)

  # check blocks
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$lhs))
  }
  # block_values <- lav_partable_block_values(partable)

  # check groups
  if (is.null(partable$group)) {
    partable$group <- rep(1L, length(partable$lhs))
  }
  group_values <- lav_partable_group_values(partable)
  ngroups <- length(group_values)

  # compute bounds per group ### TODO: add levels/classes/...
  b <- 0L
  for (g in seq_len(ngroups)) {
    # next block
    b <- b + 1L

    # for this block
    ov_names <- lavpta$vnames$ov[[b]]
    lv_names <- lavpta$vnames$lv[[b]]
    lv_names_x <- lavpta$vnames$lv.x[[b]]
    if (length(lv_names_x) > 0L) {
      lv_names_endo <- lv_names[!lv_names %in% lv_names_x]
    } else {
      lv_names_endo <- lv_names
    }
    lv_marker <- lavpta$vnames$lv.marker[[b]]

    # OV.VAR for this group
    if (lavsamplestats@missing.flag && lavdata@nlevels == 1L) {
      if (!is.null(lavh1$implied$cov[[g]])) {
        ov_var <- diag(lavh1$implied$cov[[g]])
      } else {
        ov_var <- diag(lavsamplestats@missing.h1[[g]]$sigma)
      }
    } else {
      if (lavoptions$conditional.x) {
        ov_var <- diag(lavsamplestats@res.cov[[g]])
      } else {
        ov_var <- diag(lavsamplestats@cov[[g]])
      }
    }

    # new in 0.6-17: increase observed variances for 'scaled' parameters
    # if theta parameterization
    if (theta_parameterization_flag) {
      sc_idx <- match(theta_parameterization_names, ov_names)
      ov_var[sc_idx] <- ov_var[sc_idx] / rel
    }


    # we 'process' the parameters per 'type', so we can choose
    # to apply (or not) upper/lower bounds for each type separately

    ################################
    ## 1. (residual) ov variances ##
    ################################
    par_idx <- which(partable$group == group_values[g] &
      partable$op == "~~" &
      partable$lhs %in% ov_names &
      partable$lhs == partable$rhs)

    if (length(par_idx) > 0L) {
      # lower == 0
      lower_auto[par_idx] <- 0

      # upper == var(ov)
      var_idx <- match(partable$lhs[par_idx], ov_names)
      upper_auto[par_idx] <- ov_var[var_idx]

      # if reliability > 0, adapt marker indicators only
      if (rel > 0) {
        marker_idx <- which(partable$group == group_values[g] &
          partable$op == "~~" &
          partable$lhs %in% lv_marker &
          partable$lhs == partable$rhs)
        marker_var_idx <- match(partable$lhs[marker_idx], ov_names)

        # upper = (1-REL)*OVAR
        upper_auto[marker_idx] <- (1 - rel) * ov_var[marker_var_idx]
      }

      # range
      bound_range <- upper_auto[par_idx] - pmax(lower_auto[par_idx], 0)

      # enlarge lower?
      if ("ov.var" %in% optim_bounds$lower) {
        factor <- lower_factor[which(optim_bounds$lower == "ov.var")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          diff <- abs(new_range - bound_range)
          lower_auto[par_idx] <- lower_auto[par_idx] - diff
        }
      }

      # enlarge upper?
      if ("ov.var" %in% optim_bounds$upper) {
        factor <- upper_factor[which(optim_bounds$upper == "ov.var")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          diff <- abs(new_range - bound_range)
          upper_auto[par_idx] <- upper_auto[par_idx] + diff
        } else if (is.finite(factor) && factor == 1.0) {
          # new in 0.6-20
          # enlarge anyway, but only with 0.5%
          # this is in particular useful for exogenous variances, otherwise,
          # they will always end up on the boundary
          new_range <- bound_range * 1.005
          diff <- abs(new_range - bound_range)
          upper_auto[par_idx] <- upper_auto[par_idx] + diff
        }
      }

      # min.var.ov?
      min_idx <- which(lower_auto[par_idx] < optim_bounds$min.var.ov)
      if (length(min_idx) > 0L) {
        lower_auto[par_idx[min_idx]] <- optim_bounds$min.var.ov
      }

      # requested?
      if ("ov.var" %in% optim_bounds$lower) {
        partable$lower[par_idx] <- lower_auto[par_idx]
      }
      if ("ov.var" %in% optim_bounds$upper) {
        partable$upper[par_idx] <- upper_auto[par_idx]
      }
    } # (res) ov variances

    ################################
    ## 2. (residual) lv variances ##
    ################################

    # first collect lower/upper bounds for TOTAL variances in lv.names
    lv_var_lb <- numeric(length(lv_names))
    lv_var_ub <- numeric(length(lv_names))

    if (lavoptions$std.lv) {
      lv_var_lb <- rep(1.0, length(lv_names))
      lv_var_ub <- rep(1.0, length(lv_names))
    } else {
      for (i in seq_along(lv_names)) {
        # this_lv_name <- lv_names[i]
        this_lv_marker <- lv_marker[i]

        if (nchar(this_lv_marker) > 0L && this_lv_marker %in% ov_names) {
          marker_var <- ov_var[match(this_lv_marker, ov_names)]
          lower <- marker_var - (1 - rel) * marker_var
          lv_var_lb[i] <- max(lower, optim_bounds$min.var.lv.exo)
          # LV.VAR.UB[i] <- marker.var - REL*marker.var
          lv_var_ub[i] <- marker_var

          # new in 0.6-17
          if (theta_parameterization_flag) {
            lv_var_lb[i] <- rel
          }
        } else {
          lv_var_lb[i] <- optim_bounds$min.var.lv.exo
          lv_var_ub[i] <- max(ov_var)
        }
      }
    }

    # use these bounds for the free parameters
    par_idx <- which(partable$group == group_values[g] &
      partable$op == "~~" &
      partable$lhs %in% lv_names &
      partable$lhs == partable$rhs)

    if (length(par_idx) > 0L) {
      # adjust for endogenous lv
      lv_var_lb2 <- lv_var_lb
      endo_idx <- which(lv_names %in% lv_names_endo)
      if (length(endo_idx) > 0L) {
        lv_var_lb2[endo_idx] <- optim_bounds$min.var.lv.endo
        if (optim_bounds$max.r2.lv.endo != 1) {
          lv_var_lb2[endo_idx] <-
            (1 - optim_bounds$max.r2.lv.endo) * lv_var_ub[endo_idx]
        }
      }
      exo_idx <- which(!lv_names %in% lv_names_endo)
      if (length(exo_idx) > 0L && optim_bounds$min.var.lv.exo != 0) {
        lv_var_lb2[exo_idx] <- optim_bounds$min.var.lv.exo
      }

      lower_auto[par_idx] <- lv_var_lb2[match(
        partable$lhs[par_idx],
        lv_names
      )]
      upper_auto[par_idx] <- lv_var_ub[match(
        partable$lhs[par_idx],
        lv_names
      )]

      # range
      bound_range <- upper_auto[par_idx] - pmax(lower_auto[par_idx], 0)

      # enlarge lower?
      if ("lv.var" %in% optim_bounds$lower) {
        factor <- lower_factor[which(optim_bounds$lower == "lv.var")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          diff <- abs(new_range - bound_range)
          lower_auto[par_idx] <- lower_auto[par_idx] - diff
        }
      }

      # enlarge upper?
      if ("lv.var" %in% optim_bounds$upper) {
        factor <- upper_factor[which(optim_bounds$upper == "lv.var")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          diff <- abs(new_range - bound_range)
          upper_auto[par_idx] <- upper_auto[par_idx] + diff
        }
      }

      # requested?
      if ("lv.var" %in% optim_bounds$lower) {
        partable$lower[par_idx] <- lower_auto[par_idx]
      }
      if ("lv.var" %in% optim_bounds$upper) {
        partable$upper[par_idx] <- upper_auto[par_idx]
      }
    } # lv variances


    #############################################
    ## 3. factor loadings (ov indicators only) ##
    #############################################

    # lambda_p^(u) = sqrt( upper(res.var.indicators_p) /
    #                      lower(var.factor) )

    ov_ind_names <- lavpta$vnames$ov.ind[[b]]
    par_idx <- which(partable$group == group_values[g] &
      partable$op == "=~" &
      partable$lhs %in% lv_names &
      partable$rhs %in% ov_ind_names)

    if (length(par_idx) > 0L) {
      # if negative LV variances are allowed (due to factor > 1)
      # make them equal to zero
      lv_var_lb[lv_var_lb < 0] <- 0.0

      var_all <- ov_var[match(partable$rhs[par_idx], ov_names)]
      tmp <- lv_var_lb[match(partable$lhs[par_idx], lv_names)]
      tmp[is.na(tmp)] <- 0 # just in case...
      lower_auto[par_idx] <- -1 * sqrt(var_all / tmp) # -Inf if tmp==0
      upper_auto[par_idx] <- +1 * sqrt(var_all / tmp) # +Inf if tmp==0

      # if std.lv = TRUE, force 'first' loading to be positive?
      # if(lavoptions$std.lv) {
      #    # get index 'first' indicators
      #    first.idx <- which(!duplicated(partable$lhs[par.idx]))
      #    lower.auto[par.idx][first.idx] <- 0
      # }

      # range
      bound_range <- upper_auto[par_idx] - lower_auto[par_idx]

      # enlarge lower?
      if ("loadings" %in% optim_bounds$lower) {
        factor <- lower_factor[which(optim_bounds$lower == "loadings")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          ok_idx <- is.finite(new_range)
          if (length(ok_idx) > 0L) {
            diff <- abs(new_range[ok_idx] - bound_range[ok_idx])
            lower_auto[par_idx][ok_idx] <-
              lower_auto[par_idx][ok_idx] - diff
          }
        }
      }

      # enlarge upper?
      if ("loadings" %in% optim_bounds$upper) {
        factor <- upper_factor[which(optim_bounds$upper == "loadings")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          ok_idx <- is.finite(new_range)
          if (length(ok_idx) > 0L) {
            diff <- abs(new_range[ok_idx] - bound_range[ok_idx])
            upper_auto[par_idx][ok_idx] <-
              upper_auto[par_idx][ok_idx] + diff
          }
        }
      }



      # requested?
      if ("loadings" %in% optim_bounds$lower) {
        partable$lower[par_idx] <- lower_auto[par_idx]
      }
      if ("loadings" %in% optim_bounds$upper) {
        partable$upper[par_idx] <- upper_auto[par_idx]
      }
    } # lambda


    ####################
    ## 4. covariances ##
    ####################

    # | sqrt(var(x)) sqrt(var(y)) | <= cov(x,y)

    par_idx <- which(partable$group == group_values[g] &
      partable$op == "~~" &
      partable$lhs != partable$rhs)

    if (length(par_idx) > 0L) {
      for (i in seq_along(par_idx)) {
        # this lhs/rhs
        this_lhs <- partable$lhs[par_idx[i]]
        this_rhs <- partable$rhs[par_idx[i]]

        # 2 possibilities:
        # - variances are free parameters
        # - variances are fixed (eg std.lv = TRUE)

        # var idx
        lhs_var_idx <- which(partable$group == group_values[g] &
          partable$op == "~~" &
          partable$lhs == this_lhs &
          partable$lhs == partable$rhs)
        rhs_var_idx <- which(partable$group == group_values[g] &
          partable$op == "~~" &
          partable$lhs == this_rhs &
          partable$lhs == partable$rhs)
        # upper bounds
        lhs_upper <- upper_auto[lhs_var_idx]
        rhs_upper <- upper_auto[rhs_var_idx]

        # compute upper bounds for this cov (assuming >0 vars)
        if (is.finite(lhs_upper) && is.finite(rhs_upper)) {
          upper_cov <- sqrt(lhs_upper) * sqrt(rhs_upper)
          upper_auto[par_idx[i]] <- +1 * upper_cov
          lower_auto[par_idx[i]] <- -1 * upper_cov
        }
      }

      # range
      bound_range <- upper_auto[par_idx] - lower_auto[par_idx]

      # enlarge lower?
      if ("covariances" %in% optim_bounds$lower) {
        factor <-
          lower_factor[which(optim_bounds$lower == "covariances")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          ok_idx <- is.finite(new_range)
          if (length(ok_idx) > 0L) {
            diff <- new_range[ok_idx] - bound_range[ok_idx]
            lower_auto[par_idx][ok_idx] <-
              lower_auto[par_idx][ok_idx] - diff
          }
        }
      }

      # enlarge upper?
      if ("covariances" %in% optim_bounds$upper) {
        factor <-
          upper_factor[which(optim_bounds$upper == "covariances")]
        if (is.finite(factor) && factor != 1.0) {
          new_range <- bound_range * factor
          ok_idx <- is.finite(new_range)
          if (length(ok_idx) > 0L) {
            diff <- new_range[ok_idx] - bound_range[ok_idx]
            upper_auto[par_idx][ok_idx] <-
              upper_auto[par_idx][ok_idx] + diff
          }
        }
      }

      # requested?
      if ("covariances" %in% optim_bounds$lower) {
        partable$lower[par_idx] <- lower_auto[par_idx]
      }
      if ("covariances" %in% optim_bounds$upper) {
        partable$upper[par_idx] <- upper_auto[par_idx]
      }
    } # covariances
  } # g

  # overwrite with lower.user (except -Inf)
  not_inf_idx <- which(lower_user > -Inf)
  if (length(not_inf_idx) > 0L) {
    partable$lower[not_inf_idx] <- lower_user[not_inf_idx]
  }

  # overwrite with upper.user (except +Inf)
  not_inf_idx <- which(upper_user < +Inf)
  if (length(not_inf_idx) > 0L) {
    partable$upper[not_inf_idx] <- upper_user[not_inf_idx]
  }

  # non-free
  non_free_idx <- which(partable$free == 0L)
  if (length(non_free_idx) > 0L && !is.null(partable$ustart)) {
    partable$lower[non_free_idx] <- partable$ustart[non_free_idx]
    partable$upper[non_free_idx] <- partable$ustart[non_free_idx]
  }

  partable
}
