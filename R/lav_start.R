# lav_start.R: provide starting values for model parameters
#
# YR 30/11/2010: initial version
# YR 08/06/2011: add fabin3 start values for factor loadings
# YR 14 Jan 2014: moved to lav_start.R

# fill in the 'ustart' column in a User data.frame with reasonable
# starting values, using the sample data

lav_start <- function(start_method = "default",
                      lavpartable = NULL,
                      lavsamplestats = NULL,
                      lavh1 = NULL, # fixme: only use lavh1?
                      model_type = "sem",
                      reflect = FALSE, # rotation only
                      samplestats_flag = TRUE,
                      order_lv_by = "none" # rotation only
                      ) {
  # check arguments
  stopifnot(is.list(lavpartable))

  # categorical?
  categorical <- any(lavpartable$op == "|")

  # correlation structure?
  correlation <- any(lavpartable$op == "~*~")

  # composites?
  composites <- any(lavpartable$op == "<~")

  # conditional.x?
  conditional_x <- any(lavpartable$exo == 1L &
    lavpartable$op %in% c("~", "<~"))
  # ord.names <- unique(lavpartable$lhs[ lavpartable$op == "|" ])

  # nlevels?
  nlevels <- lav_pt_nlevels(lavpartable)

  # reflect/order.lv.by
  if (is.null(reflect)) {
    reflect <- FALSE
  }
  if (is.null(order_lv_by)) {
    order_lv_by <- "index"
  }

  # check start.method
  #if (mimic == "lavaan") {
  start_initial <- "lavaan"
  #} else if (mimic == "Mplus") {
  #  start.initial <- "mplus"
  #} else {
  #  # FIXME: use LISREL/EQS/AMOS/.... schemes
  #  start.initial <- "lavaan"
  #}

  # start.method
  start_user <- NULL
  if (is.character(start_method)) {
    start_method_lc <- tolower(start_method)
  if (start_method_lc != "simple" && !samplestats_flag) {
      start_method_lc <- start_method <- "simple"
    }
    if (start_method_lc == "default") {
      # nothing to do
    } else if (start_method == "simple") {
      start <- numeric(length(lavpartable$ustart))
      # if(categorical || correlation) {
      start[which(lavpartable$op == "=~")] <- 0.7
      start[which(lavpartable$op == "<~")] <- 1
      # } else {
      #    start[ which(lavpartable$op == "=~") ] <- 1.0
      # }
      start[which(lavpartable$op == "~*~")] <- 1.0
      ov_names_ord <- lav_pt_vnames(lavpartable, "ov.ord")
      var_idx <- which(lavpartable$op == "~~" &
        lavpartable$lhs == lavpartable$rhs &
        !(lavpartable$lhs %in% ov_names_ord))
      start[var_idx] <- 1.0
      user_idx <- which(!is.na(lavpartable$ustart))
      start[user_idx] <- lavpartable$ustart[user_idx]
      return(start) # assuming fixed.x = FALSE!
    } else if (start_method == "est") {
      return(lavpartable$est)
    } else if (start_method_lc %in% c("simple", "lavaan")) {
      start_initial <- start_method_lc
    } else {
     lav_msg_stop(gettext("unknown value for start argument"))
    }
  } else if (is.list(start_method)) {
    start_user <- start_method
  } else if (is.numeric(start_method)) {
    nx_free <- sum(lavpartable$free > 0L)
    if (length(start_method) != nx_free) {
     lav_msg_stop(gettextf(
       "start argument contains %1$s elements; but parameter table
       expects %2$s free parameters.", length(start_method), nx_free))
    }
    lavpartable$ustart[lavpartable$free > 0L] <- start_method
  } else if (inherits(start_method, "lavaan")) {
    start_user <- parTable(start_method)
  }

  # check model list elements, if provided
  if (!is.null(start_user)) {
    if (is.null(start_user$lhs) ||
      is.null(start_user$op) ||
      is.null(start_user$rhs)) {
     lav_msg_stop(gettext(
       "problem with start argument: model list does not contain
       all elements: lhs/op/rhs"))
    }
    if (!is.null(start_user$est)) {
      # excellent, we got an est column; nothing to do
    } else if (!is.null(start_user$start)) {
      # no est column, but we use the start column
      start_user$est <- start_user$start
    } else if (!is.null(start_user$ustart)) {
      # not ideal, but better than nothing
      start_user$est <- start_user$ustart
    } else {
     lav_msg_stop(gettext(
       "problem with start argument: could not find est/start column
       in model list"))
    }
  }


  # global settings
  # 0. everything is zero
  start <- numeric(length(lavpartable$ustart))

  # 1. =~ factor loadings:
  if (categorical || correlation) {
    # if std.lv=TRUE, 0.8 is too large
    start[which(lavpartable$op == "=~")] <- 0.7
  } else {
    start[which(lavpartable$op == "=~")] <- 1.0
  }

  # 2. (residual) lv variances for latent variables
  lv_names <- lav_pt_vnames(lavpartable, "lv") # all groups
  lv_var_idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lv_names &
    lavpartable$lhs == lavpartable$rhs)
  start[lv_var_idx] <- 0.05
  # start[lv.var.idx] <- 0.5 # new in 0.6-2? (for optim.parscale = "stand")

  # 3. latent response scales (if any)
  delta_idx <- which(lavpartable$op == "~*~")
  start[delta_idx] <- 1.0


  # group-specific settings
  ngroups <- lav_pt_ngroups(lavpartable)

  # for now, if no group column, add one (again), until we rewrite
  # this function to handle block/group hybrid settings
  if (is.null(lavpartable$group) && ngroups == 1L) {
    lavpartable$group <- rep(1L, length(lavpartable$lhs))
    lavpartable$group[lavpartable$block == 0L] <- 0L
  }

  # group values
  group_values <- lav_pt_group_values(lavpartable)

  for (g in 1:ngroups) {

    # info from user model for this group
    if (conditional_x) {
      ov_names <- lav_pt_vnames(lavpartable, "ov.nox",
                                                group = group_values[g])
    } else {
      ov_names <- lav_pt_vnames(lavpartable, "ov",
                                                group = group_values[g])
    }
    if (categorical) {
      ov_names_num <- lav_pt_vnames(lavpartable, "ov.num",
                                                group = group_values[g])
      ov_names_ord <- lav_pt_vnames(lavpartable, "ov.ord",
                                                group = group_values[g])
    } else {
      ov_names_num <- ov_names
    }
    lv_names <- lav_pt_vnames(lavpartable, "lv",
                                                 group = group_values[g])
    lv_names_efa <- lav_pt_vnames(lavpartable, "lv.efa",
                                                 group = group_values[g])
    ov_names_x <- lav_pt_vnames(lavpartable, "ov.x",
                                                 group = group_values[g])
    ov_ind_c <- lav_pt_vnames(lavpartable, "ov.cind",
                                                 group = group_values[g])
    lv_names_c <- lav_pt_vnames(lavpartable, "lv.composite",
                                                 group = group_values[g])

    # just for the nlevels >1 case
    ov_names <- unique(unlist(ov_names))
    ov_names_num <- unique(unlist(ov_names_num))
    lv_names <- unique(unlist(lv_names))
    lv_names_efa <- unique(unlist(lv_names_efa))
    ov_names_x <- unique(unlist(ov_names_x))
    ov_ind_c <- unique(unlist(ov_ind_c))
    lv_names_c <- unique(unlist(lv_names_c))
    lv_names_noc <- lv_names[!lv_names %in% lv_names_c]

    # residual ov variances (including exo/ind, to be overridden)
    ov_var_idx <- which(lavpartable$group == group_values[g] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% ov_names_num &
      lavpartable$lhs == lavpartable$rhs)
    sample_var_idx <- match(lavpartable$lhs[ov_var_idx], ov_names)
    if (model_type == "unrestricted") {
    # this does not work if conditional.x = TRUE...
      if (!is.null(lavh1$implied$cov[[g]])) {
        h1_cov <- lavh1$implied$cov[[g]]
      } else if (!is.null(lavsamplestats@missing.h1[[g]])) {
        h1_cov <- lavsamplestats@missing.h1[[g]]$sigma
      } else {
        h1_cov <- lavsamplestats@cov[[g]]
      }
      start[ov_var_idx] <- diag(h1_cov)[sample_var_idx]
    } else {
      #if (start.initial == "mplus") {
      #  if (conditional.x && nlevels == 1L) {
      #    start[ov.var.idx] <-
      #      (1.0 - 0.50) * lavsamplestats@res.var[[1L]][sample.var.idx]
      #  } else {
      #    start[ov.var.idx] <-
      #      (1.0 - 0.50) * lavsamplestats@var[[1L]][sample.var.idx]
      #  }
      #} else {
        if (conditional_x && nlevels == 1L) {
          start[ov_var_idx] <-
            (1.0 - 0.50) * diag(lavsamplestats@res.cov[[g]])[sample_var_idx]
        } else {
          start[ov_var_idx] <-
            (1.0 - 0.50) * diag(lavsamplestats@cov[[g]])[sample_var_idx]
        }
      #}
      # composite indicators: fill in total variances
      if (composites) {
        start[ov_var_idx] <- diag(lavsamplestats@cov[[g]])[sample_var_idx]
      }
    }

    # 1-fac measurement models: loadings, psi, theta
    if (start_initial %in% c("lavaan") &&
      model_type %in% c("sem", "cfa")) {
      # fabin3 estimator (2sls) of Hagglund (1982) per factor
      for (f in lv_names_noc) {
        # not for efa factors
        if (f %in% lv_names_efa) {
          next
        }
        lambda_idx <- which(lavpartable$lhs == f &
          lavpartable$op == "=~" &
          lavpartable$group == group_values[g])
        # standardized?
        std_lv <- FALSE
        var_f_idx <- which(lavpartable$lhs == f &
          lavpartable$op == "~~" &
          lavpartable$group == group_values[g] &
          lavpartable$rhs == f)
        if (length(var_f_idx) > 0L &&
          all(lavpartable$free[var_f_idx] == 0) &&
          all(lavpartable$ustart[var_f_idx] == 1)) {
          std_lv <- TRUE
        }

        # no second order
        if (any(lavpartable$rhs[lambda_idx] %in% lv_names)) next

        # get observed indicators for this latent variable
        ov_idx <- match(lavpartable$rhs[lambda_idx], ov_names)
        if (length(ov_idx) > 0L && !any(is.na(ov_idx))) {
          if (lavsamplestats@missing.flag && nlevels == 1L) {
            if (!is.null(lavh1$implied$cov[[g]])) {
              h1_cov <- lavh1$implied$cov[[g]]
            } else {
              h1_cov <- lavsamplestats@missing.h1[[g]]$sigma
            }
            cov_1 <- h1_cov[ov_idx, ov_idx, drop = FALSE]
          } else {
            if (conditional_x && nlevels == 1L) {
              cov_1 <- lavsamplestats@res.cov[[g]][ov_idx,
                ov_idx,
                drop = FALSE
              ]
            } else {
              cov_1 <- lavsamplestats@cov[[g]][ov_idx,
                ov_idx,
                drop = FALSE
              ]
            }
          }

          # fabin for 1-factor
          fabin <- lav_cfa_1fac_fabin(cov_1,
            std_lv = std_lv,
            lambda_only = TRUE,
            method = "fabin3"
          )

          # factor loadings
          tmp <- fabin$lambda
          tmp[!is.finite(tmp)] <- 1.0 # just in case (eg 0/0)

          # check for negative triad if nvar=3L (new in 0.6-8)
          if (!is.null(fabin$neg.triad) && fabin$neg.triad) {
            if (std_lv) {
              tmp <- rep(0.7, length(tmp))
            } else {
              tmp <- rep(1.0, length(tmp))
            }
          }

          # check for negative marker
          if (!std_lv && !is.na(lavpartable$ustart[lambda_idx[1]]) &&
              lavpartable$ustart[lambda_idx[1]] < 0) {
            tmp <- -1 * tmp
          }

          start[lambda_idx] <- tmp

          # factor variance
          # if(!std.lv) {
          #    start[var.f.idx] <- fabin$psi
          #    # if residual var, make smaller
          #    y.idx <- which(lavpartable$lhs == f &
          #                   lavpartable$group == group.values[g] &
          #                   lavpartable$op == "~")
          #    if(length(y.idx) > 0L) {
          #        # how much explained variance do we expect?
          #        # we take 0.50
          #        start[var.f.idx] <- 0.5 * start[var.f.idx]
          #    }
          #    # no negative variances (we get these if we have an
          #    # inconsistent triad (eg, covariance signs are +,+,-)
          #    if(start[var.f.idx] < 0) {
          #        start[var.f.idx] <- 0.05
          #    }
          # }

          # NOTE: fabin (sometimes) gives residual variances
          # that are larger than the original variances...

          # residual variances -- order?
          # res.idx <- which(lavpartable$lhs %in% ov.names[ov.idx] &
          #                 lavpartable$op == "~~" &
          #                 lavpartable$group == group.values[g] &
          #                 lavpartable$rhs == lavpartable$lhs)
          # start[res.idx] <- fabin$theta

          # negative variances?
          # neg.idx <- which(start[res.idx] < 0)
          # if(length(neg.idx) > 0L) {
          #    start[res.idx][neg.idx] <- 0.05
          # }
        }
      } # fabin3

      # efa?
      nefa <- lav_pt_nefa(lavpartable)
      if (nefa > 0L) {
        efa_values <- lav_pt_efa_values(lavpartable)

        for (set in seq_len(nefa)) {
          # determine ov idx for this set
          ov_efa <-
            unique(lavpartable$rhs[lavpartable$op == "=~" &
              lavpartable$block == g &
              lavpartable$efa == efa_values[set]])
          lv_efa <-
            unique(lavpartable$lhs[lavpartable$op == "=~" &
              lavpartable$block == g &
              lavpartable$efa == efa_values[set]])
          lambda_idx <- which(lavpartable$lhs %in% lv_efa &
            lavpartable$op == "=~" &
            lavpartable$group == group_values[g])

          theta_idx <- which(lavpartable$lhs %in% ov_efa &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs &
            lavpartable$group == group_values[g])

          # get observed indicators for these EFA lv variables
          ov_idx <- match(
            unique(lavpartable$rhs[lambda_idx]),
            ov_names
          )

          if (length(ov_idx) > 0L && !any(is.na(ov_idx))) {
            if (lavsamplestats@missing.flag && nlevels == 1L) {
              if (!is.null(lavh1$implied$cov[[g]])) {
                h1_cov <- lavh1$implied$cov[[g]]
              } else {
                h1_cov <- lavsamplestats@missing.h1[[g]]$sigma
              }
              cov_1 <- h1_cov[ov_idx, ov_idx, drop = FALSE]
            } else {
              if (conditional_x) {
                cov_1 <- lavsamplestats@res.cov[[g]][ov_idx,
                  ov_idx,
                  drop = FALSE
                ]
              } else {
                cov_1 <- lavsamplestats@cov[[g]][ov_idx,
                  ov_idx,
                  drop = FALSE
                ]
              }
            }

            # EFA solution with zero upper-right corner
            efa_1 <- lav_efa_extraction(
              s = cov_1,
              nfactors = length(lv_efa),
              method = "ML",
              order_lv_by = order_lv_by,
              # order.lv.by = "none",
              # reflect = reflect,
              reflect = FALSE,
              corner = TRUE
            )

            # factor loadings
            tmp <- as.numeric(efa_1$LAMBDA)
            tmp[!is.finite(tmp)] <- 1.0 # just in case (eg 0/0)
            start[lambda_idx] <- tmp

            # residual variances
            tmp <- diag(efa_1$THETA)
            tmp[!is.finite(tmp)] <- 1.0 # just in case
            start[theta_idx] <- tmp
          }
        } # set
      } # efa
    } # factor loadings

    if (model_type == "unrestricted") {
      # fill in 'covariances' from lavsamplestats
      cov_idx <- which(lavpartable$group == group_values[g] &
        lavpartable$op == "~~" &
        lavpartable$lhs != lavpartable$rhs)
      lhs_idx <- match(lavpartable$lhs[cov_idx], ov_names)
      rhs_idx <- match(lavpartable$rhs[cov_idx], ov_names)
      if (!is.null(lavh1$implied$cov[[g]])) {
        h1_cov <- lavh1$implied$cov[[g]]
      } else if (!is.null(lavsamplestats@missing.h1[[g]])) {
        h1_cov <- lavsamplestats@missing.h1[[g]]$sigma
      } else {
        h1_cov <- lavsamplestats@cov[[g]]
      }
      start[cov_idx] <- h1_cov[cbind(lhs_idx, rhs_idx)]
    }

    # composites
    if (composites) {
      std_lv <- FALSE
      var_f_idx <- which(lavpartable$lhs %in% lv_names_c &
        lavpartable$op == "~~" &
        lavpartable$group == group_values[g] &
        lavpartable$rhs %in% lv_names_c)
      if (length(var_f_idx) > 0L &&
        all(lavpartable$free[var_f_idx] == 0) &&
        !all(is.na(lavpartable$ustart[var_f_idx])) &&
        all(lavpartable$ustart[var_f_idx] == 1)) {
        std_lv <- TRUE
      }

      # weights
      cidx <- which(lavpartable$group == group_values[g] &
        lavpartable$op == "<~")
      if (std_lv) {
        start[cidx] <- 0.10
      } else {
        start[cidx] <- 1
      }

      # fill in 'covariances' from lavsamplestats
      cov_idx <- which(lavpartable$group == group_values[g] &
        lavpartable$op == "~~" &
        lavpartable$rhs %in% ov_ind_c &
        lavpartable$lhs != lavpartable$rhs)
      lhs_idx <- match(lavpartable$lhs[cov_idx], ov_names)
      rhs_idx <- match(lavpartable$rhs[cov_idx], ov_names)
      if (!is.null(lavh1$implied$cov[[g]])) {
        h1_cov <- lavh1$implied$cov[[g]]
      } else if (!is.null(lavsamplestats@missing.h1[[g]])) {
        h1_cov <- lavsamplestats@missing.h1[[g]]$sigma
      } else {
        h1_cov <- lavsamplestats@cov[[g]]
      }
      start[cov_idx] <- h1_cov[cbind(lhs_idx, rhs_idx)]
    }

    # variances of ordinal variables - set to 1.0
    if (categorical) {
      ov_var_ord_idx <- which(lavpartable$group == group_values[g] &
        lavpartable$op == "~~" &
        lavpartable$lhs %in% ov_names_ord &
        lavpartable$lhs == lavpartable$rhs)
      start[ov_var_ord_idx] <- 1.0
    }

    # 3g) intercepts/means
    ov_int_idx <- which(lavpartable$group == group_values[g] &
      lavpartable$op == "~1" &
      lavpartable$lhs %in% ov_names)
    sample_int_idx <- match(lavpartable$lhs[ov_int_idx], ov_names)
    if (lavsamplestats@missing.flag && nlevels == 1L) {
      if (!is.null(lavh1$implied$mean[[g]])) {
        h1_mean <- lavh1$implied$mean[[g]]
      } else if (!is.null(lavsamplestats@missing.h1[[g]])) {
        h1_mean <- lavsamplestats@missing.h1[[g]]$mu
      } else {
        h1_mean <- lavsamplestats@mean[[g]]
      }
      start[ov_int_idx] <- h1_mean[sample_int_idx]
    } else {
      if (conditional_x && nlevels == 1L) {
        start[ov_int_idx] <- lavsamplestats@res.int[[g]][sample_int_idx]
      } else {
        start[ov_int_idx] <- lavsamplestats@mean[[g]][sample_int_idx]
      }
    }

    # TODo: if marker.int.zero = TRUE, set lv means to marker means,
    #       and the non-marker means to
    #       lavsamplestats@mean[[g]] - LAMBDA %*% ALPHA
    #       where ALPHA = means of the markers

    # 4g) thresholds
    th_idx <- which(lavpartable$group == group_values[g] &
                      lavpartable$op == "|")
    if (length(th_idx) > 0L) {
      th_names_lavpartable <- paste(lavpartable$lhs[th_idx], "|",
        lavpartable$rhs[th_idx],
        sep = ""
      )
      th_names_sample <-
        lavsamplestats@th.names[[g]][lavsamplestats@th.idx[[g]] > 0L]
      # th.names.sample should identical to
      # lav_pt_vnames(lavpartable, "th", group = group.values[g])
      if (conditional_x && nlevels == 1L) {
        th_values <-
          lavsamplestats@res.th[[g]][lavsamplestats@th.idx[[g]] > 0L]
      } else {
        th_values <-
          lavsamplestats@th[[g]][lavsamplestats@th.idx[[g]] > 0L]
      }
      start[th_idx] <- th_values[match(
        th_names_lavpartable,
        th_names_sample
      )]
    }

    # 5g) exogenous `fixed.x' covariates
    if (length(ov_names_x) > 0) {
      exo_idx <- which(lavpartable$group == group_values[g] &
        lavpartable$op == "~~" &
        lavpartable$lhs %in% ov_names_x &
        lavpartable$rhs %in% ov_names_x)
      if (!conditional_x) {
        row_idx <- match(lavpartable$lhs[exo_idx], ov_names)
        col_idx <- match(lavpartable$rhs[exo_idx], ov_names)
        if (lavsamplestats@missing.flag && nlevels == 1L) {
          if (!is.null(lavh1$implied$cov[[g]])) {
            h1_cov <- lavh1$implied$cov[[g]]
          } else if (!is.null(lavsamplestats@missing.h1[[g]])) {
            h1_cov <- lavsamplestats@missing.h1[[g]]$sigma
          } else {
            h1_cov <- lavsamplestats@cov[[g]]
          }
          start[exo_idx] <- h1_cov[cbind(row_idx, col_idx)]
          # using slightly smaller starting values for free
          # variance/covariances (fixed.x = FALSE);
          # this somehow avoids false convergence in saturated models
          nobs <- lavsamplestats@nobs[[g]]
          this_idx <- which(seq_along(lavpartable$free) %in% exo_idx &
                              lavpartable$free > 0L)
          start[this_idx] <- start[this_idx] * (nobs - 1) / nobs
        } else {
          start[exo_idx] <- lavsamplestats@cov[[g]][cbind(row_idx, col_idx)]
        }
      } else {
        # cov.x
        row_idx <- match(lavpartable$lhs[exo_idx], ov_names_x)
        col_idx <- match(lavpartable$rhs[exo_idx], ov_names_x)
        start[exo_idx] <- lavsamplestats@cov.x[[g]][cbind(
          row_idx,
          col_idx
        )]
        # mean.x
        exo_int_idx <- which(lavpartable$group == group_values[g] &
          lavpartable$op == "~1" &
          lavpartable$lhs %in% ov_names_x)
        int_idx <- match(lavpartable$lhs[exo_int_idx], ov_names_x)
        start[exo_int_idx] <- lavsamplestats@mean.x[[g]][int_idx]
      }
    }

    # 6b. exogenous lv variances if single indicator -- new in 0.5-21
    lv_x <- lav_pt_vnames(lavpartable, "lv.x", group = group_values[g])
    # FIXME: also for multilevel?
    lv_x <- unique(unlist(lv_x))
    if (length(lv_x) > 0L) {
      for (ll in lv_x) {
        ind_idx <- which(lavpartable$op == "=~" &
          lavpartable$lhs == ll &
          lavpartable$group == group_values[g])
        if (length(ind_idx) == 1L) {
          single_ind <- lavpartable$rhs[ind_idx]
          single_fvar_idx <- which(lavpartable$op == "~~" &
            lavpartable$lhs == ll &
            lavpartable$rhs == ll &
            lavpartable$group == group_values[g])
          single_var_idx <- which(lavpartable$op == "~~" &
            lavpartable$lhs == single_ind &
            lavpartable$rhs == single_ind &
            lavpartable$group == group_values[g])
          # user-defined residual variance
          # fixme: we take the first, in case we have multiple matches
          # (eg nlevels)
          single_var <- lavpartable$ustart[single_var_idx[1]]
          if (is.na(single_var)) {
            single_var <- 1
          }
          ov_idx <- match(single_ind, ov_names)
          if (conditional_x && nlevels == 1L) {
            ov_var <- diag(lavsamplestats@res.cov[[g]])[ov_idx]
          } else {
            ov_var <- diag(lavsamplestats@cov[[g]])[ov_idx]
          }
          # take (1 - (rvar/ov.var) * ov.var
          tmp <- (1 - (single_var / ov_var)) * ov_var
          # just in case
          if (is.na(tmp) || tmp < 0.05) {
            tmp <- 0.05
          }
          start[single_fvar_idx] <- tmp
        }
      }
    }

    # 7g) regressions "~" # new in 0.6-10
    if (length(lv_names) == 0L && nlevels == 1L && !conditional_x) {
      # observed only
      reg_idx <- which(lavpartable$group == group_values[g] &
        lavpartable$op == "~")
      if (length(reg_idx) > 0L) {
        eqs_y <- unique(lavpartable$lhs[reg_idx])
        ny <- length(eqs_y)
        for (i in seq_len(ny)) {
          y_name <- eqs_y[i]
          start_idx <- which(lavpartable$group == group_values[g] &
            lavpartable$op == "~" &
            lavpartable$lhs == y_name)
          x_names <- lavpartable$rhs[start_idx]
          cov_1 <- lavsamplestats@cov[[g]]
          y_idx <- match(y_name, ov_names)
          x_idx <- match(x_names, ov_names)
          s_xx <- cov_1[x_idx, x_idx, drop = FALSE]
          s_xy <- cov_1[x_idx, y_idx, drop = FALSE]
          # regression coefficient(s)
          beta_i <- try(solve(s_xx, s_xy), silent = TRUE)
          if (inherits(beta_i, "try-error")) {
            start[start_idx] <- beta_i <- rep(0, length(start_idx))
          } else {
            start[start_idx] <- drop(beta_i)
          }
          # residual variance
          res_idx <- which(lavpartable$group == group_values[g] &
            lavpartable$op == "~~" &
            lavpartable$lhs == y_name &
            lavpartable$rhs == y_name)
          res_val <- cov_1[y_idx, y_idx] - drop(crossprod(beta_i, s_xy))
          if (res_val > 0.001 * cov_1[y_idx, y_idx] &&
            res_val < 0.999 * cov_1[y_idx, y_idx]) {
            start[res_idx] <- res_val
          } else {
            # do nothing (keep what we have)
          }
          # intercept
          int_idx <- which(lavpartable$group == group_values[g] &
            lavpartable$op == "~1" &
            lavpartable$lhs == y_name)
          if (length(int_idx) > 0L) {
            mean_1 <- lavsamplestats@mean[[g]]
            ybar <- mean_1[y_idx]
            xbar <- mean_1[x_idx]
            int_val <- ybar - drop(crossprod(beta_i, xbar))
            if (is.finite(int_val)) {
              start[int_idx] <- int_val
            }
          }
        }
      }
    }

    #  # 8 latent variances (new in 0.6-2)
    #  lv.names.y <- lav_pt_vnames(lavpartable,
    #                               "lv.y", group = group.values[g])
    #  lv.names.x <- lav_pt_vnames(lavpartable,
    #                               "lv.x", group = group.values[g])
    #  # multilevel? take first level only
    #  if(is.list(lv.names.y)) {
    #      lv.names.y <- unlist(lv.names.y) # for now
    #  }
    #  if(is.list(lv.names.x)) {
    #      lv.names.x <- unlist(lv.names.x) # for now
    #  }
    #  lv.names.xy <- unique(c(lv.names.x, lv.names.y))


    #  if(length(lv.names.xy) > 0L) {
    #      free.var.idx <- which(lavpartable$op == "~~" &
    #                            lavpartable$lhs %in% lv.names.xy &
    #                            lavpartable$rhs == lavpartable$lhs &
    #                            lavpartable$group == group.values[g])
    #      if(length(free.var.idx) > 0L) {
    #          this.lv.names <- lavpartable$lhs[free.var.idx]
    #          for(v in seq_len(length(free.var.idx))) {
    #              # single marker item?
    #              ind.idx <- which(lavpartable$op == "=~" &
    #                               lavpartable$lhs %in% this.lv.names[v] &
    #                               #lavpartable$rhs %in% ov.names.num &
    #                               lavpartable$free == 0L &
    #                               lavpartable$group == group.values[g])
    #              if(length(ind.idx) == 0) {
    #                  next
    #              } else if(length(ind.idx) > 1L) {
    #                  # FIXME! perhaps a random effect? do something clever
    #                  next
    #              } else if(length(ind.idx) == 1L) {
    #                  marker.ind <- lavpartable$rhs[ind.idx]
    #                  ov.idx <- match(marker.ind, ov.names)
    #                  if(conditional.x) {
    #                      ov.var <- diag(lavsamplestats@res.cov[[g]])[ov.idx]
    #                  } else {
    #                      ov.var <- diag(lavsamplestats@cov[[g]])[ov.idx]
    #                  }
    #
    #                  # exogenous? assume rel = 0.50
    #                  lambda <- lavpartable$ustart[ind.idx]
    #                  tmp <- (0.50 * ov.var)/lambda^2
    #                  if(this.lv.names[v] %in% lv.names.y) {
    #                      # endogenous, assume R2 = 0.2
    #                      tmp <- 0.8 * tmp
    #                  }
    #                  # within variance?
    #                  if(nlevels > 1L &&
    #                     lavpartable$level[ free.var.idx[v] ] == 1L) {
    #                      tmp <- tmp * 0.75
    #                  }
    #                  # between variance?
    #                  if(nlevels > 1L &&
    #                     lavpartable$level[ free.var.idx[v] ] > 1L) {
    #                      tmp <- tmp * 0.25
    #                  }
    #                  # just in case
    #                  if(is.na(tmp) || tmp < 0.05) {
    #                      tmp <- 0.05
    #                  }
    #                  start[ free.var.idx[v] ] <- tmp
    #              }
    #          } # v
    #      } # free.var.idx
    #  } # lv var

    # nlevels > 1L
    if (nlevels > 1L) {
      level_values <- lav_pt_level_values(lavpartable)
      # Note: ov.names.x contains all levels within a group!
      if (length(ov_names_x) > 0) {
        for (l in 1:nlevels) {
          # block number
          block <- (g - 1L) * nlevels + l

          this_block_x <- lav_pt_vnames(lavpartable, "ov.x",
            block = block
          )
          this_block_ov <- lav_pt_vnames(lavpartable, "ov",
            block = block
          )
          if (length(this_block_x) == 0L) {
            next
          }

          # var/cov
          exo_idx <- which(lavpartable$group == group_values[g] &
            lavpartable$level == level_values[l] &
            lavpartable$op == "~~" &
            lavpartable$lhs %in% this_block_x &
            lavpartable$rhs %in% this_block_x)

          if (is.null(lavh1$implied$cov[[1]])) {
            row_idx <- match(lavpartable$lhs[exo_idx], ov_names)
            col_idx <- match(lavpartable$rhs[exo_idx], ov_names)
            if (l == 1L) {
              cov_1 <- lavsamplestats@YLp[[g]][[2]]$S.PW.start
            } else {
              cov_1 <- lavsamplestats@YLp[[g]][[l]]$Sigma.B
            }
          } else {
            row_idx <- match(lavpartable$lhs[exo_idx], this_block_ov)
            col_idx <- match(lavpartable$rhs[exo_idx], this_block_ov)
            cov_1 <- lavh1$implied$cov[[block]]
          }
          # make sure starting values for variances are positive
          neg_idx <- which(diag(cov_1) < 0.001)
          if (length(neg_idx) > 0L) {
            diag(cov_1)[neg_idx] <- 0.001
          }
          start[exo_idx] <- cov_1[cbind(row_idx, col_idx)]

          # intercepts
          ov_int_idx <- which(lavpartable$group == group_values[g] &
            lavpartable$level == level_values[l] &
            lavpartable$op == "~1" &
            lavpartable$lhs %in% this_block_x)

          if (is.null(lavh1$implied$mean[[1]])) {
            idx <- match(lavpartable$lhs[ov_int_idx], ov_names)
            if (l == 1L) {
              int <- lavsamplestats@YLp[[g]][[2]]$Mu.W
            } else {
              int <- lavsamplestats@YLp[[g]][[l]]$Mu.B.start
            }
          } else {
            idx <- match(lavpartable$lhs[ov_int_idx], this_block_ov)
            int <- lavh1$implied$mean[[block]]
          }
          start[ov_int_idx] <- int[idx]

          # new in 0.6-12
          # very special case: conditional.x with a combination of
          # split-x and regular-x
          # here, we must:
          # 1) replace var/cov of split-x by *residual* varcov
          #    after regressing out regular-x
          # 2) replace means of split-x by intercepts
          # 3) fill split-x ~ regular-x regression coefficients
          if (conditional_x) {
            if (is.null(lavh1$implied$cov[[l]])) {
             lav_msg_stop(gettext(
               "lavh1 information is needed; please rerun with h1 = TRUE"))
            }
            blocks_within_group <- (g - 1L) * nlevels + seq_len(nlevels)
            other_block_names <- lav_pt_vnames(lavpartable, "ov",
                                    block = blocks_within_group[-block])
            ov_names_x_block <- this_block_x
            idx <- which(ov_names_x_block %in% other_block_names)
            if (length(idx) > 0L) {
              ov_names_x_block <- ov_names_x_block[-idx]
            }
            ov_names_x1 <- this_block_x[!this_block_x %in% ov_names_x_block]
            ov_names_x2 <- ov_names_x_block
            nx1 <- length(ov_names_x1) # split x
            nx2 <- length(ov_names_x2) # regular  x
            if (nx1 > 0L && nx2 > 0L) {
              # COV
              c1_idx <- match(ov_names_x1, this_block_ov)
              c2_idx <- match(ov_names_x2, this_block_ov)

              cov_y <- cov_1[c1_idx, c1_idx, drop = FALSE]
              cov_x <- cov_1[c2_idx, c2_idx, drop = FALSE]
              cov_yx <- cov_1[c1_idx, c2_idx, drop = FALSE]
              cov_xy <- cov_1[c2_idx, c1_idx, drop = FALSE]
              cov_xinv_yx <- solve(cov_x, cov_xy)
              res_cov <- cov_y - cov_yx %*% cov_xinv_yx

              res_cov_idx <- which(lavpartable$group == group_values[g] &
                lavpartable$level == level_values[l] &
                lavpartable$op == "~~" &
                lavpartable$lhs %in% ov_names_x1 &
                lavpartable$rhs %in% ov_names_x1)

              row_idx <- match(lavpartable$lhs[res_cov_idx], ov_names_x1)
              col_idx <- match(lavpartable$rhs[res_cov_idx], ov_names_x1)
              start[res_cov_idx] <- res_cov[cbind(row_idx, col_idx)]

              # INT
              int_y <- int[c1_idx]
              int_x <- int[c2_idx]
              res_int <- int_y - t(cov_xinv_yx) %*% int_x

              res_int_idx <- which(lavpartable$group == group_values[g] &
                lavpartable$level == level_values[l] &
                lavpartable$op == "~1" &
                lavpartable$lhs %in% ov_names_x1)
              idx <- match(lavpartable$lhs[res_int_idx], ov_names_x1)
              start[res_int_idx] <- res_int[idx]

              # REG
              reg_idx <- which(lavpartable$group == group_values[g] &
                lavpartable$level == level_values[l] &
                lavpartable$op == "~" &
                lavpartable$lhs %in% ov_names_x1 &
                lavpartable$rhs %in% ov_names_x2)
              row_idx <- match(lavpartable$lhs[reg_idx], ov_names_x1)
              col_idx <- match(lavpartable$rhs[reg_idx], ov_names_x2)
              start[reg_idx] <- t(cov_xinv_yx)[cbind(row_idx, col_idx)]
            } # special case
          } # conditional.x
        } # levels
      } # fixed.x
    } # nlevels > 1L
  } # groups




  # group weights
  group_idx <- which(lavpartable$lhs == "group" &
    lavpartable$op == "%")
  if (length(group_idx) > 0L) {
    ngroups <- length(group_idx)
    # prop <- rep(1/ngroups, ngroups)
    # use last group as reference
    # start[group.idx] <- log(prop/prop[ngroups])

    # poisson version
    start[group_idx] <- log(rep(lavsamplestats@ntotal / ngroups, ngroups))
  }

  # growth models:
  # - compute starting values for mean latent variables
  # - compute starting values for variance latent variables
  if (start_initial %in% c("lavaan") &&
    model_type == "growth") {
    ### DEBUG ONLY
    # lv.var.idx <- which(lavpartable$op == "~~"                &
    #                lavpartable$lhs %in% lv.names &
    #                lavpartable$lhs == lavpartable$rhs)

    ### DEBUG ONLY
    # lv.int.idx <- which(lavpartable$op == "~1"         &
    #                    lavpartable$lhs %in% lv.names)
  }

  # adjust if outside bounds -- new in 0.6-6
  if (!is.null(lavpartable$lower)) {
    bad_idx <- which(start < lavpartable$lower)
    if (length(bad_idx)) {
      start[bad_idx] <- lavpartable$lower[bad_idx]
    }
  }
  if (!is.null(lavpartable$upper)) {
    bad_idx <- which(start > lavpartable$upper)
    if (length(bad_idx)) {
      start[bad_idx] <- lavpartable$upper[bad_idx]
    }
  }


  # override if the model syntax contains explicit starting values (free only)
  # user.idx <- which(!is.na(lavpartable$ustart) &
  #                  lavpartable$user != 7L) # new in 0.6-7, if rotation and
  #                                          # and we change the order of lv's
  user_idx <- which(!is.na(lavpartable$ustart) & lavpartable$free > 0L)
  start[user_idx] <- lavpartable$ustart[user_idx]

  # override if a user list with starting values is provided
  # we only look at the 'est' column for now
  if (!is.null(start_user)) {
    if (is.null(lavpartable$group)) {
      lavpartable$group <- rep(1L, length(lavpartable$lhs))
    }
    if (is.null(start_user$group)) {
      start_user$group <- rep(1L, length(start_user$lhs))
    }

    # FIXME: avoid for loop!!!
    for (i in seq_along(lavpartable$lhs)) {
      # find corresponding parameters
      lhs <- lavpartable$lhs[i]
      op <- lavpartable$op[i]
      rhs <- lavpartable$rhs[i]
      grp <- lavpartable$group[i]

      start_user_idx <- which(start_user$lhs == lhs &
        start_user$op == op &
        start_user$rhs == rhs &
        start_user$group == grp)
      if (length(start_user_idx) == 1L &&
        is.finite(start_user$est[start_user_idx])) {
        start[i] <- start_user$est[start_user_idx]
      }
    }
  }

  # override fixed values with ustart values
  user_idx <- which(!is.na(lavpartable$ustart) & lavpartable$free == 0L)
  start[user_idx] <- lavpartable$ustart[user_idx]


  # new in 0.6-21
  # if any thresholds, make sure they are in increasing order
  # (could be an issue if the first and second were fixed to 0 and 1)
  for (g in 1:ngroups) {
    th_idx <- which(lavpartable$group == group_values[g] &
                    lavpartable$op == "|")
    if (length(th_idx) > 0L) {
      # for every ov.ord, check if t1 < t2 < t3 < ...
      ov_ord <- unique(lavpartable$lhs[th_idx])
      for (oo in ov_ord) {
        this_idx <- which(lavpartable$group == group_values[g] &
                          lavpartable$op == "|" & lavpartable$lhs == oo)
        this_free_idx <- which(lavpartable$group == group_values[g] &
                               lavpartable$op == "|" & lavpartable$lhs == oo &
                               lavpartable$free > 0)
        item_th <- start[this_idx]
        if (length(item_th) > 1L && length(this_free_idx) > 0 &&
                                          !all(diff(item_th) > 0)) {
          start[this_free_idx] <-
            cumsum(abs(item_th))[match(this_free_idx, this_idx)]
        }
      }
    }
  }


  # final check: no NaN or other non-finite values
  bad_idx <- which(!is.finite(start))
  if (length(bad_idx) > 0L) {
    cat("starting values:\n")
    print(start)
     lav_msg_warn(gettext(
       "some starting values are non-finite; replacing them with 0.5;
       please provide better starting values."))
    start[bad_idx] <- 0.5
  }

  if (lav_debug()) {
    cat("lavaan DEBUG: lavaanStart\n")
    print(start)
  }

  start
}

# backwards compatibility
# StartingValues <- lav_start

# sanity check: (user-specified) variances smaller than covariances
# but not for composites, as we have not 'set' their variances yet
lav_start_check_cov <- function(lavpartable = NULL, start = lavpartable$start) {
  nblocks <- lav_pt_nblocks(lavpartable)
  block_values <- lav_pt_block_values(lavpartable)

  for (g in 1:nblocks) {

    lv_names_c <- lav_pt_vnames(lavpartable, "lv.composite", block = g)

    # collect all non-zero covariances
    cov_idx <- which(lavpartable$op == "~~" &
      lavpartable$block == block_values[g] &
      !lavpartable$lhs %in% lv_names_c &
      lavpartable$lhs != lavpartable$rhs &
      !lavpartable$exo &
      start != 0)

    # for each covariance, use corresponding variances to standardize;
    # the end result should not exceed abs(1)
    for (cc in seq_along(cov_idx)) {
      this_cov_idx <- cov_idx[cc]

      # find corresponding variances
      var_lhs <- lavpartable$lhs[this_cov_idx]
      var_rhs <- lavpartable$rhs[this_cov_idx]

      var_lhs_idx <- which(lavpartable$op == "~~" &
        lavpartable$block == block_values[g] &
        lavpartable$lhs == var_lhs &
        lavpartable$lhs == lavpartable$rhs)

      var_rhs_idx <- which(lavpartable$op == "~~" &
        lavpartable$block == block_values[g] &
        lavpartable$lhs == var_rhs &
        lavpartable$lhs == lavpartable$rhs)

      var_lhs_value <- start[var_lhs_idx]
      var_rhs_value <- start[var_rhs_idx]

      block_txt <- ""
      if (nblocks > 1L) {
        block_txt <- paste(" [in block ", g, "]", sep = "")
      }

      # check for zero variances
      if (var_lhs_value == 0 || var_rhs_value == 0) {
        # this can only happen if it is user-specified
        # cov.idx free? set it to zero
        if (start[this_cov_idx] == 0) {
          # nothing to do
        } else if (lavpartable$free[this_cov_idx] > 0L) {
            lav_msg_warn(gettextf(
             "non-zero covariance element set to zero, due to fixed-to-zero
             variances variables involved are: %s", var_lhs), var_rhs,
             block_txt
            )
          start[this_cov_idx] <- 0
        } else {
          lav_msg_stop(gettextf(
            "please provide better fixed values for (co)variances;
            variables involved are: %s ", var_lhs), var_rhs, block_txt
          )
        }
        next
      }

      # which one is the smallest? abs() in case of negative variances
      if (abs(var_lhs_value) < abs(var_rhs_value)) {
        var_min_idx <- var_lhs_idx
        var_max_idx <- var_rhs_idx
      } else {
        var_min_idx <- var_rhs_idx
        var_max_idx <- var_lhs_idx
      }

      # check
      cor_1 <- abs(start[this_cov_idx] / sqrt(var_lhs_value * var_rhs_value))

      # NOTE: we treat this as an unconditional COR!

      if (!is.finite(cor_1)) {
        # force simple values
           lav_msg_warn(gettextf(
            "starting values imply NaN for a correlation value; variables
            involved are: %s", var_lhs), var_rhs, block_txt
          )
        start[var_lhs_idx] <- 1
        start[var_rhs_idx] <- 1
        start[this_cov_idx] <- 0
      } else if (cor_1 > 1) {
        txt <- gettextf(
          "starting values imply a correlation larger than 1; variables
          involved are: %1$s %2$s %3$s", var_lhs, var_rhs, block_txt)

        # three ways to fix it: rescale cov12, var1 or var2

        # we prefer a free parameter, and not user-specified
        if (lavpartable$free[this_cov_idx] > 0L &&
          is.na(lavpartable$ustart[this_cov_idx])) {
             lav_msg_warn(gettext(txt))
          start[this_cov_idx] <- start[this_cov_idx] / (cor_1 * 1.1)
        } else if (lavpartable$free[var_min_idx] > 0L &&
          is.na(lavpartable$ustart[var_min_idx])) {
             lav_msg_warn(gettext(txt))
          start[var_min_idx] <- start[var_min_idx] * (cor_1 * 1.1)^2
        } else if (lavpartable$free[var_max_idx] > 0L &&
          is.na(lavpartable$ustart[var_max_idx])) {
             lav_msg_warn(gettext(txt))
          start[var_max_idx] <- start[var_max_idx] * (cor_1 * 1.1)^2

          # not found? try just a free parameter
        } else if (lavpartable$free[this_cov_idx] > 0L) {
             lav_msg_warn(gettext(txt))
          start[this_cov_idx] <- start[this_cov_idx] / (cor_1 * 1.1)
        } else if (lavpartable$free[var_min_idx] > 0L) {
             lav_msg_warn(gettext(txt))
          start[var_min_idx] <- start[var_min_idx] * (cor_1 * 1.1)^2
        } else if (lavpartable$free[var_max_idx] > 0L) {
             lav_msg_warn(gettext(txt))
          start[var_max_idx] <- start[var_max_idx] * (cor_1 * 1.1)^2

          # nothing? abort or warn (and fail later...): warn
        } else {
             lav_msg_warn(gettext(txt))
          # lav_msg_stop(gettextf(
          # "please provide better fixed values for (co)variances;
          #  variables involved are: %s ", var.lhs), var.rhs, block.txt)
        }
      } # COR > 1
    } # cov.idx
  }
  start
}
