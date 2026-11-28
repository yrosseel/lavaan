# add parameter bounds to the parameter table
# lavoptions$optim.bounds
lav_pt_add_bounds <- function(partable = NULL,
                                    lavh1 = NULL,
                                    lavdata = NULL,
                                    lavsamplestats = NULL,
                                    lavoptions = NULL) {
  # multilevel?
  nlevels <- lav_pt_nlevels(partable)

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

    # multilevel: the bounds are computed per level, using the variances
    # of the unrestricted (h1) model for that level; this is only possible
    # for random-intercept models
    if (nlevels > 1L &&
      (length(lavoptions$optim.bounds$lower) > 0L ||
        length(lavoptions$optim.bounds$upper) > 0L)) {
      if (!is.null(partable$rv) && any(nchar(partable$rv) > 0L)) {
        lav_msg_warn(gettext(
          "automatic bounds not available (yet) for multilevel models with
           random slopes"
        ))
        return(partable)
      }
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

    if (is.null(optim_bounds$min_reliability_marker)) {
      optim_bounds$min_reliability_marker <- 0.0
    } else {
      if (optim_bounds$min_reliability_marker < 0 ||
        optim_bounds$min_reliability_marker > 1.0) {
        lav_msg_stop(gettextf(
          "optim_bounds$min_reliability_marker is out of range: %s",
          optim_bounds$min_reliability_marker
        ))
      }
    }

    if (is.null(optim_bounds$min_var_ov)) {
      optim_bounds$min_var_ov <- -Inf
    }

    if (is.null(optim_bounds$min_var_lv_exo)) {
      optim_bounds$min_var_lv_exo <- 0.0
    }

    if (is.null(optim_bounds$min_var_lv_endo)) {
      optim_bounds$min_var_lv_endo <- 0.0
    }

    if (is.null(optim_bounds$max_r2_lv_endo)) {
      optim_bounds$max_r2_lv_endo <- 1.0
    }

    # multilevel only: the factor by which the between-level variances (of
    # the h1 model) are multiplied before they are used to compute upper
    # bounds; unlike the (total) variances in a single-level model, the
    # model-implied between-level variances can be (much) larger than the
    # between-level variances of the h1 model, in particular if the
    # intraclass correlations are small and/or the number of clusters is
    # small
    if (is.null(optim_bounds$upper_factor_between)) {
      optim_bounds$upper_factor_between <- 1.5
    } else if (!is.numeric(optim_bounds$upper_factor_between) ||
      length(optim_bounds$upper_factor_between) != 1L ||
      !is.finite(optim_bounds$upper_factor_between) ||
      optim_bounds$upper_factor_between < 1.0) {
      lav_msg_stop(gettext(
        "optim_bounds$upper_factor_between should be a single number that
         is equal to or larger than 1."))
    }

    if (is.null(optim_bounds$lower_factor)) {
      optim_bounds$lower_factor <- rep(1.0, length(optim_bounds$lower))
    } else {
      if (length(optim_bounds$lower_factor) == 1L &&
        is.numeric(optim_bounds$lower_factor)) {
        optim_bounds$lower_factor <- rep(
          optim_bounds$lower_factor,
          length(optim_bounds$lower)
        )
      } else if (length(optim_bounds$lower_factor) !=
        length(optim_bounds$lower)) {
        lav_msg_stop(
          gettext("length(optim_bounds$lower_factor) is not equal to
                  length(optim_bounds$lower)")
        )
      }
    }
    lower_factor <- optim_bounds$lower_factor

    if (is.null(optim_bounds$upper_factor)) {
      optim_bounds$upper_factor <- rep(1.0, length(optim_bounds$upper))
    } else {
      if (length(optim_bounds$upper_factor) == 1L &&
        is.numeric(optim_bounds$upper_factor)) {
        optim_bounds$upper_factor <- rep(
          optim_bounds$upper_factor,
          length(optim_bounds$upper)
        )
      } else if (length(optim_bounds$upper_factor) !=
        length(optim_bounds$upper)) {
        lav_msg_stop(
          gettext("length(optim_bounds$lower_factor) is not equal to
                  length(optim_bounds$upper)")
        )
      }
    }
    upper_factor <- optim_bounds$upper_factor
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
  rel <- optim_bounds$min_reliability_marker

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

  lavpta <- lav_pt_attributes(partable)

  # check blocks
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$lhs))
  }

  # check groups
  if (is.null(partable$group)) {
    partable$group <- rep(1L, length(partable$lhs))
  }
  group_values <- lav_pt_group_values(partable)
  ngroups <- length(group_values)

  # compute bounds per block (group, or level within group)
  for (b in seq_len(ngroups * nlevels)) {
    # the group of this block
    g <- (b - 1L) %/% nlevels + 1L

    # the parameters of this block
    if (nlevels > 1L) {
      block_flag <- partable$block == b
    } else {
      block_flag <- partable$group == group_values[g]
    }

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
    if (nlevels > 1L) {
      # multilevel: the (within or between) variances of the h1 model
      ov_var <- lav_pt_bounds_ov_var_level(
        b = b, g = g, level = (b - 1L) %% nlevels + 1L,
        ov_names = ov_names, lavh1 = lavh1, lavdata = lavdata,
        lavsamplestats = lavsamplestats
      )
      if (is.null(ov_var)) {
        # no information about the variances per level
        partable$lower <- lower_user
        partable$upper <- upper_user
        return(partable)
      }
      # the h1 estimate of a between variance can be (close to) zero;
      # bounds that are based on such a variance would (almost) fix the
      # parameters of this variable
      ov_var[!is.finite(ov_var) | ov_var < 0] <- 0
    } else if (lavsamplestats@missing.flag && lavdata@nlevels == 1L) {
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

    # the variances that are used for the upper bounds: the same, unless
    # this is a between-level block (see upper_factor_between above)
    ov_var_ub <- ov_var
    if (nlevels > 1L && (b - 1L) %% nlevels + 1L > 1L) {
      ov_var_ub <- ov_var * optim_bounds$upper_factor_between
    }


    # we 'process' the parameters per 'type', so we can choose
    # to apply (or not) upper/lower bounds for each type separately

    ################################
    ## 1. (residual) ov variances ##
    ################################
    par_idx <- which(block_flag &
      partable$op == "~~" &
      partable$lhs %in% ov_names &
      partable$lhs == partable$rhs)

    if (length(par_idx) > 0L) {
      # lower == 0
      lower_auto[par_idx] <- 0

      # upper == var(ov)
      var_idx <- match(partable$lhs[par_idx], ov_names)
      upper_auto[par_idx] <- ov_var_ub[var_idx]

      # if reliability > 0, adapt marker indicators only
      if (rel > 0) {
        marker_idx <- which(block_flag &
          partable$op == "~~" &
          partable$lhs %in% lv_marker &
          partable$lhs == partable$rhs)
        marker_var_idx <- match(partable$lhs[marker_idx], ov_names)

        # upper = (1-REL)*OVAR
        upper_auto[marker_idx] <- (1 - rel) * ov_var_ub[marker_var_idx]
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

      # min_var_ov?
      min_idx <- which(lower_auto[par_idx] < optim_bounds$min_var_ov)
      if (length(min_idx) > 0L) {
        lower_auto[par_idx[min_idx]] <- optim_bounds$min_var_ov
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
          lv_var_lb[i] <- max(lower, optim_bounds$min_var_lv_exo)
          # LV.VAR.UB[i] <- marker.var - REL*marker.var
          lv_var_ub[i] <- ov_var_ub[match(this_lv_marker, ov_names)]

          # new in 0.6-17
          if (theta_parameterization_flag) {
            lv_var_lb[i] <- rel
          }
        } else {
          lv_var_lb[i] <- optim_bounds$min_var_lv_exo
          lv_var_ub[i] <- max(ov_var_ub)
        }
      }
    }

    # use these bounds for the free parameters
    par_idx <- which(block_flag &
      partable$op == "~~" &
      partable$lhs %in% lv_names &
      partable$lhs == partable$rhs)

    if (length(par_idx) > 0L) {
      # adjust for endogenous lv
      lv_var_lb2 <- lv_var_lb
      endo_idx <- which(lv_names %in% lv_names_endo)
      if (length(endo_idx) > 0L) {
        lv_var_lb2[endo_idx] <- optim_bounds$min_var_lv_endo
        if (optim_bounds$max_r2_lv_endo != 1) {
          lv_var_lb2[endo_idx] <-
            (1 - optim_bounds$max_r2_lv_endo) * lv_var_ub[endo_idx]
        }
      }
      exo_idx <- which(!lv_names %in% lv_names_endo)
      if (length(exo_idx) > 0L && optim_bounds$min_var_lv_exo != 0) {
        lv_var_lb2[exo_idx] <- optim_bounds$min_var_lv_exo
      }

      lower_auto[par_idx] <- lv_var_lb2[match(
        partable$lhs[par_idx],
        lv_names
      )]
      upper_auto[par_idx] <- lv_var_ub[match(
        partable$lhs[par_idx],
        lv_names
      )]

      # new in 0.6-20: if std.lv = TRUE, the latent variances are normally
      # fixed to 1.0; but in some settings (eg multiple groups +
      # group.equal = "loadings") the latent variances in groups 2+ are
      # free parameters. For those free parameters, the [1, 1] bounds are
      # wrong (they would fix the variance to 1.0). Use data-driven bounds
      # instead (std.lv has no marker indicator, so no marker-based bounds).
      if (lavoptions$std.lv) {
        std_free_idx <- which(partable$free[par_idx] > 0L)
        if (length(std_free_idx) > 0L) {
          fidx <- par_idx[std_free_idx]
          lower_auto[fidx] <- optim_bounds$min_var_lv_exo
          upper_auto[fidx] <- max(ov_var_ub)
          # endogenous lv: respect min_var_lv_endo / max_r2_lv_endo
          endo_f <- which(partable$lhs[fidx] %in% lv_names_endo)
          if (length(endo_f) > 0L) {
            lower_auto[fidx[endo_f]] <- optim_bounds$min_var_lv_endo
            if (optim_bounds$max_r2_lv_endo != 1) {
              lower_auto[fidx[endo_f]] <-
                (1 - optim_bounds$max_r2_lv_endo) * upper_auto[fidx[endo_f]]
            }
          }
        }
      }

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
    par_idx <- which(block_flag &
      partable$op == "=~" &
      partable$lhs %in% lv_names &
      partable$rhs %in% ov_ind_names)

    if (length(par_idx) > 0L) {
      # if negative LV variances are allowed (due to factor > 1)
      # make them equal to zero
      lv_var_lb[lv_var_lb < 0] <- 0.0

      var_all <- ov_var_ub[match(partable$rhs[par_idx], ov_names)]
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

    par_idx <- which(block_flag &
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
        lhs_var_idx <- which(block_flag &
          partable$op == "~~" &
          partable$lhs == this_lhs &
          partable$lhs == partable$rhs)
        rhs_var_idx <- which(block_flag &
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

# multilevel: the variances of the observed variables of block 'b' (level
# 'level' in group 'g'), to be used for the automatic bounds
#
# - if available, we use the unrestricted (h1) model
# - if not (eg h1 = FALSE), we use the (moment-based) estimates of the
#   within and between covariance matrices that are also used to compute
#   the starting values
# - if both are missing, NULL is returned
lav_pt_bounds_ov_var_level <- function(b = 1L, g = 1L, level = 1L,
                                       ov_names = character(0L),
                                       lavh1 = NULL,
                                       lavdata = NULL,
                                       lavsamplestats = NULL) {
  nvar <- length(ov_names)

  # 1. h1
  cov_b <- NULL
  if (length(lavh1$implied$cov) >= b) {
    cov_b <- lavh1$implied$cov[[b]]
  }
  if (!is.null(cov_b) && nrow(cov_b) == nvar) {
    return(diag(cov_b))
  }

  # 2. moment-based
  ylp <- NULL
  if (!is.null(lavsamplestats) && length(lavsamplestats@YLp) >= g &&
    length(lavsamplestats@YLp[[g]]) >= 2L) {
    ylp <- lavsamplestats@YLp[[g]][[2L]]
  }
  if (level > 2L || is.null(ylp) || is.null(ylp$Sigma.W) ||
    is.null(ylp$Sigma.B)) {
    return(NULL)
  }
  ov_idx <- match(ov_names, lavdata@ov.names[[g]])
  if (anyNA(ov_idx)) {
    return(NULL)
  }
  if (level == 1L) {
    out <- diag(ylp$Sigma.W)[ov_idx]
  } else {
    out <- diag(ylp$Sigma.B)[ov_idx]
  }

  out
}
