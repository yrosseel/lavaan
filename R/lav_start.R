# lav_start.R: provide starting values for model parameters
#
# YR 30/11/2010: initial version
# YR 08/06/2011: add fabin3 start values for factor loadings
# YR 14 Jan 2014: moved to lav_start.R

# fill in the 'ustart' column in a User data.frame with reasonable
# starting values, using the sample data

lav_start <- function(start.method = "default",
                      lavpartable = NULL,
                      lavsamplestats = NULL,
                      lavh1 = NULL, # fixme: only use lavh1?
                      model.type = "sem",
                      mimic = "lavaan",
                      reflect = FALSE, # rotation only
					  samplestats.flag = TRUE,
                      order.lv.by = "none", # rotation only
                      debug = FALSE) {
  # check arguments
  stopifnot(is.list(lavpartable))

  # categorical?
  categorical <- any(lavpartable$op == "|")

  # correlation structure?
  correlation <- any(lavpartable$op == "~*~")

  # conditional.x?
  conditional.x <- any(lavpartable$exo == 1L &
    lavpartable$op %in% c("~", "<~"))
  # ord.names <- unique(lavpartable$lhs[ lavpartable$op == "|" ])

  # nlevels?
  nlevels <- lav_partable_nlevels(lavpartable)

  # reflect/order.lv.by
  if (is.null(reflect)) {
    reflect <- FALSE
  }
  if (is.null(order.lv.by)) {
    order.lv.by <- "index"
  }

  # check start.method
  if (mimic == "lavaan") {
    start.initial <- "lavaan"
  } else if (mimic == "Mplus") {
    start.initial <- "mplus"
  } else {
    # FIXME: use LISREL/EQS/AMOS/.... schemes
    start.initial <- "lavaan"
  }

  # start.method
  start.user <- NULL
  if (is.character(start.method)) {
    start.method.lc <- tolower(start.method)
	if (start.method.lc != "simple" && !samplestats.flag) {
	    start.method.lc <- start.method <- "simple"
    }
    if (start.method.lc == "default") {
      # nothing to do
    } else if (start.method == "simple") {
      start <- numeric(length(lavpartable$ustart))
      # if(categorical || correlation) {
      start[which(lavpartable$op == "=~")] <- 0.7
      # } else {
      #    start[ which(lavpartable$op == "=~") ] <- 1.0
      # }
      start[which(lavpartable$op == "~*~")] <- 1.0
      ov.names.ord <- vnames(lavpartable, "ov.ord")
      var.idx <- which(lavpartable$op == "~~" &
        lavpartable$lhs == lavpartable$rhs &
        !(lavpartable$lhs %in% ov.names.ord))
      start[var.idx] <- 1.0
      user.idx <- which(!is.na(lavpartable$ustart))
      start[user.idx] <- lavpartable$ustart[user.idx]
      return(start) # assuming fixed.x = FALSE!
    } else if (start.method == "est") {
      return(lavpartable$est)
    } else if (start.method.lc %in% c("simple", "lavaan", "mplus")) {
      start.initial <- start.method.lc
    } else {
     lav_msg_stop(gettext("unknown value for start argument"))
    }
  } else if (is.list(start.method)) {
    start.user <- start.method
  } else if (is.numeric(start.method)) {
    nx.free <- sum(lavpartable$free > 0L)
    if (length(start.method) != nx.free) {
     lav_msg_stop(gettextf(
       "start argument contains %1$s elements; but parameter table
       expects %2$s free parameters.", length(start.method), nx.free))
    }
    lavpartable$ustart[lavpartable$free > 0L] <- start.method
  } else if (inherits(start.method, "lavaan")) {
    start.user <- parTable(start.method)
  }

  # check model list elements, if provided
  if (!is.null(start.user)) {
    if (is.null(start.user$lhs) ||
      is.null(start.user$op) ||
      is.null(start.user$rhs)) {
     lav_msg_stop(gettext(
       "problem with start argument: model list does not contain
       all elements: lhs/op/rhs"))
    }
    if (!is.null(start.user$est)) {
      # excellent, we got an est column; nothing to do
    } else if (!is.null(start.user$start)) {
      # no est column, but we use the start column
      start.user$est <- start.user$start
    } else if (!is.null(start.user$ustart)) {
      # no ideal, but better than nothing
      start.user$est <- start.user$ustart
    } else {
     lav_msg_stop(gettext(
       "problem with start argument: could not find est/start column
       in model list"))
    }
  }


  # global settings
  # 0. everyting is zero
  start <- numeric(length(lavpartable$ustart))

  # 1. =~ factor loadings:
  if (categorical || correlation) {
    # if std.lv=TRUE, 0.8 is too large
    start[which(lavpartable$op == "=~")] <- 0.7
  } else {
    start[which(lavpartable$op == "=~")] <- 1.0
  }

  # 2. (residual) lv variances for latent variables
  lv.names <- vnames(lavpartable, "lv") # all groups
  lv.var.idx <- which(lavpartable$op == "~~" &
    lavpartable$lhs %in% lv.names &
    lavpartable$lhs == lavpartable$rhs)
  start[lv.var.idx] <- 0.05
  # start[lv.var.idx] <- 0.5 # new in 0.6-2? (for optim.parscale = "stand")

  # 3. latent response scales (if any)
  delta.idx <- which(lavpartable$op == "~*~")
  start[delta.idx] <- 1.0


  # group-specific settings
  ngroups <- lav_partable_ngroups(lavpartable)

  # for now, if no group column, add one (again), until we rewrite
  # this function to handle block/group hybrid settings
  if (is.null(lavpartable$group) && ngroups == 1L) {
    lavpartable$group <- rep(1L, length(lavpartable$lhs))
    lavpartable$group[lavpartable$block == 0L] <- 0L
  }


  for (g in 1:ngroups) {
    # group values (not necessarily 1,2,... anymore)
    group.values <- lav_partable_group_values(lavpartable)

    # info from user model for this group
    if (conditional.x) {
      ov.names <- vnames(lavpartable, "ov.nox", group = group.values[g])
    } else {
      ov.names <- vnames(lavpartable, "ov", group = group.values[g])
    }
    if (categorical) {
      ov.names.num <- vnames(lavpartable, "ov.num", group = group.values[g])
      ov.names.ord <- vnames(lavpartable, "ov.ord", group = group.values[g])
    } else {
      ov.names.num <- ov.names
    }
    lv.names <- vnames(lavpartable, "lv", group = group.values[g])
    lv.names.efa <- vnames(lavpartable, "lv.efa", group = group.values[g])
    ov.names.x <- vnames(lavpartable, "ov.x", group = group.values[g])

    # just for the nlevels >1 case
    ov.names <- unique(unlist(ov.names))
    ov.names.num <- unique(unlist(ov.names.num))
    lv.names <- unique(unlist(lv.names))
    lv.names.efa <- unique(unlist(lv.names.efa))
    ov.names.x <- unique(unlist(ov.names.x))


    # residual ov variances (including exo/ind, to be overriden)
    ov.var.idx <- which(lavpartable$group == group.values[g] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% ov.names.num &
      lavpartable$lhs == lavpartable$rhs)
    sample.var.idx <- match(lavpartable$lhs[ov.var.idx], ov.names)
    if (model.type == "unrestricted") {
      if (!is.null(lavsamplestats@missing.h1[[g]])) {
        start[ov.var.idx] <-
          diag(lavsamplestats@missing.h1[[g]]$sigma)[sample.var.idx]
      } else {
        start[ov.var.idx] <-
          diag(lavsamplestats@cov[[g]])[sample.var.idx]
      }
    } else {
      if (start.initial == "mplus") {
        if (conditional.x && nlevels == 1L) {
          start[ov.var.idx] <-
            (1.0 - 0.50) * lavsamplestats@res.var[[1L]][sample.var.idx]
        } else {
          start[ov.var.idx] <-
            (1.0 - 0.50) * lavsamplestats@var[[1L]][sample.var.idx]
        }
      } else {
        if (conditional.x && nlevels == 1L) {
          start[ov.var.idx] <-
            (1.0 - 0.50) * diag(lavsamplestats@res.cov[[g]])[sample.var.idx]
        } else {
          start[ov.var.idx] <-
            (1.0 - 0.50) * diag(lavsamplestats@cov[[g]])[sample.var.idx]
        }
      }
    }

    # 1-fac measurement models: loadings, psi, theta
    if (start.initial %in% c("lavaan", "mplus") &&
      model.type %in% c("sem", "cfa")) {
      # fabin3 estimator (2sls) of Hagglund (1982) per factor
      for (f in lv.names) {
        # not for efa factors
        if (f %in% lv.names.efa) {
          next
        }
        lambda.idx <- which(lavpartable$lhs == f &
          lavpartable$op == "=~" &
          lavpartable$group == group.values[g])
        # standardized?
        std.lv <- FALSE
        var.f.idx <- which(lavpartable$lhs == f &
          lavpartable$op == "~~" &
          lavpartable$group == group.values[g] &
          lavpartable$rhs == f)
        if (length(var.f.idx) > 0L &&
          all(lavpartable$free[var.f.idx] == 0) &&
          all(lavpartable$ustart[var.f.idx] == 1)) {
          std.lv <- TRUE
        }

        # no second order
        if (any(lavpartable$rhs[lambda.idx] %in% lv.names)) next

        # get observed indicators for this latent variable
        ov.idx <- match(lavpartable$rhs[lambda.idx], ov.names)
        if (length(ov.idx) > 0L && !any(is.na(ov.idx))) {
          if (lavsamplestats@missing.flag && nlevels == 1L) {
            COV <- lavsamplestats@missing.h1[[g]]$sigma[ov.idx,
              ov.idx,
              drop = FALSE
            ]
          } else {
            if (conditional.x && nlevels == 1L) {
              COV <- lavsamplestats@res.cov[[g]][ov.idx,
                ov.idx,
                drop = FALSE
              ]
            } else {
              COV <- lavsamplestats@cov[[g]][ov.idx,
                ov.idx,
                drop = FALSE
              ]
            }
          }

          # fabin for 1-factor
          fabin <- lav_cfa_1fac_fabin(COV,
            std.lv = std.lv,
            lambda.only = TRUE,
            method = "fabin3"
          )

          # factor loadings
          tmp <- fabin$lambda
          tmp[!is.finite(tmp)] <- 1.0 # just in case (eg 0/0)

          # check for negative triad if nvar=3L (new in 0.6-8)
          if (!is.null(fabin$neg.triad) && fabin$neg.triad) {
            if (std.lv) {
              tmp <- rep(0.7, length(tmp))
            } else {
              tmp <- rep(1.0, length(tmp))
            }
          }
          start[lambda.idx] <- tmp

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
      nefa <- lav_partable_nefa(lavpartable)
      if (nefa > 0L) {
        efa.values <- lav_partable_efa_values(lavpartable)

        for (set in seq_len(nefa)) {
          # determine ov idx for this set
          ov.efa <-
            unique(lavpartable$rhs[lavpartable$op == "=~" &
              lavpartable$block == g &
              lavpartable$efa == efa.values[set]])
          lv.efa <-
            unique(lavpartable$lhs[lavpartable$op == "=~" &
              lavpartable$block == g &
              lavpartable$efa == efa.values[set]])
          lambda.idx <- which(lavpartable$lhs %in% lv.efa &
            lavpartable$op == "=~" &
            lavpartable$group == group.values[g])

          theta.idx <- which(lavpartable$lhs %in% ov.efa &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs &
            lavpartable$group == group.values[g])

          # get observed indicators for these EFA lv variables
          ov.idx <- match(
            unique(lavpartable$rhs[lambda.idx]),
            ov.names
          )

          if (length(ov.idx) > 0L && !any(is.na(ov.idx))) {
            if (lavsamplestats@missing.flag && nlevels == 1L) {
              COV <- lavsamplestats@missing.h1[[g]]$sigma[ov.idx,
                ov.idx,
                drop = FALSE
              ]
            } else {
              if (conditional.x) {
                COV <- lavsamplestats@res.cov[[g]][ov.idx,
                  ov.idx,
                  drop = FALSE
                ]
              } else {
                COV <- lavsamplestats@cov[[g]][ov.idx,
                  ov.idx,
                  drop = FALSE
                ]
              }
            }

            # EFA solution with zero upper-right corner
            EFA <- lav_efa_extraction(
              S = COV,
              nfactors = length(lv.efa),
              method = "ML",
              order.lv.by = order.lv.by,
              # order.lv.by = "none",
              # reflect = reflect,
              reflect = FALSE,
              corner = TRUE
            )

            # factor loadings
            tmp <- as.numeric(EFA$LAMBDA)
            tmp[!is.finite(tmp)] <- 1.0 # just in case (eg 0/0)
            start[lambda.idx] <- tmp

            # residual variances
            tmp <- diag(EFA$THETA)
            tmp[!is.finite(tmp)] <- 1.0 # just in case
            start[theta.idx] <- tmp
          }
        } # set
      } # efa
    } # factor loadings

    if (model.type == "unrestricted") {
      # fill in 'covariances' from lavsamplestats
      cov.idx <- which(lavpartable$group == group.values[g] &
        lavpartable$op == "~~" &
        lavpartable$lhs != lavpartable$rhs)
      lhs.idx <- match(lavpartable$lhs[cov.idx], ov.names)
      rhs.idx <- match(lavpartable$rhs[cov.idx], ov.names)
      if (!is.null(lavsamplestats@missing.h1[[g]])) {
        start[cov.idx] <- lavsamplestats@missing.h1[[g]]$sigma[
          cbind(lhs.idx, rhs.idx)
        ]
      } else {
        start[cov.idx] <- lavsamplestats@cov[[g]][
          cbind(lhs.idx, rhs.idx)
        ]
      }
    }

    # variances of ordinal variables - set to 1.0
    if (categorical) {
      ov.var.ord.idx <- which(lavpartable$group == group.values[g] &
        lavpartable$op == "~~" &
        lavpartable$lhs %in% ov.names.ord &
        lavpartable$lhs == lavpartable$rhs)
      start[ov.var.ord.idx] <- 1.0
    }

    # 3g) intercepts/means
    ov.int.idx <- which(lavpartable$group == group.values[g] &
      lavpartable$op == "~1" &
      lavpartable$lhs %in% ov.names)
    sample.int.idx <- match(lavpartable$lhs[ov.int.idx], ov.names)
    if (lavsamplestats@missing.flag && nlevels == 1L) {
      start[ov.int.idx] <- lavsamplestats@missing.h1[[g]]$mu[sample.int.idx]
    } else {
      if (conditional.x && nlevels == 1L) {
        start[ov.int.idx] <- lavsamplestats@res.int[[g]][sample.int.idx]
      } else {
        start[ov.int.idx] <- lavsamplestats@mean[[g]][sample.int.idx]
      }
    }

    # TODo: if marker.int.zero = TRUE, set lv means to marker means,
    #       and the non-marker means to
    #       lavsamplestats@mean[[g]] - LAMBDA %*% ALPHA
    #       where ALPHA = means of the markers

    # 4g) thresholds
    th.idx <- which(lavpartable$group == group.values[g] &
                      lavpartable$op == "|")
    if (length(th.idx) > 0L) {
      th.names.lavpartable <- paste(lavpartable$lhs[th.idx], "|",
        lavpartable$rhs[th.idx],
        sep = ""
      )
      th.names.sample <-
        lavsamplestats@th.names[[g]][lavsamplestats@th.idx[[g]] > 0L]
      # th.names.sample should identical to
      # vnames(lavpartable, "th", group = group.values[g])
      if (conditional.x && nlevels == 1L) {
        th.values <-
          lavsamplestats@res.th[[g]][lavsamplestats@th.idx[[g]] > 0L]
      } else {
        th.values <-
          lavsamplestats@th[[g]][lavsamplestats@th.idx[[g]] > 0L]
      }
      start[th.idx] <- th.values[match(
        th.names.lavpartable,
        th.names.sample
      )]
    }

    # 5g) exogenous `fixed.x' covariates
    if (length(ov.names.x) > 0) {
      exo.idx <- which(lavpartable$group == group.values[g] &
        lavpartable$op == "~~" &
        lavpartable$lhs %in% ov.names.x &
        lavpartable$rhs %in% ov.names.x)
      if (!conditional.x) {
        row.idx <- match(lavpartable$lhs[exo.idx], ov.names)
        col.idx <- match(lavpartable$rhs[exo.idx], ov.names)
        if (lavsamplestats@missing.flag && nlevels == 1L) {
          start[exo.idx] <-
            lavsamplestats@missing.h1[[g]]$sigma[cbind(row.idx, col.idx)]
          # using slightly smaller starting values for free
          # variance/covariances (fixed.x = FALSE);
          # this somehow avoids false convergence in saturated models
          nobs <- lavsamplestats@nobs[[g]]
          this.idx <- which(seq_len(length(lavpartable$free)) %in% exo.idx &
                              lavpartable$free > 0L)
          start[this.idx] <- start[this.idx] * (nobs - 1) / nobs
        } else {
          start[exo.idx] <- lavsamplestats@cov[[g]][cbind(row.idx, col.idx)]
        }
      } else {
        # cov.x
        row.idx <- match(lavpartable$lhs[exo.idx], ov.names.x)
        col.idx <- match(lavpartable$rhs[exo.idx], ov.names.x)
        start[exo.idx] <- lavsamplestats@cov.x[[g]][cbind(
          row.idx,
          col.idx
        )]
        # mean.x
        exo.int.idx <- which(lavpartable$group == group.values[g] &
          lavpartable$op == "~1" &
          lavpartable$lhs %in% ov.names.x)
        int.idx <- match(lavpartable$lhs[exo.int.idx], ov.names.x)
        start[exo.int.idx] <- lavsamplestats@mean.x[[g]][int.idx]
      }
    }

    # 6b. exogenous lv variances if single indicator -- new in 0.5-21
    lv.x <- vnames(lavpartable, "lv.x", group = group.values[g])
    # FIXME: also for multilevel?
    lv.x <- unique(unlist(lv.x))
    if (length(lv.x) > 0L) {
      for (ll in lv.x) {
        ind.idx <- which(lavpartable$op == "=~" &
          lavpartable$lhs == ll &
          lavpartable$group == group.values[g])
        if (length(ind.idx) == 1L) {
          single.ind <- lavpartable$rhs[ind.idx]
          single.fvar.idx <- which(lavpartable$op == "~~" &
            lavpartable$lhs == ll &
            lavpartable$rhs == ll &
            lavpartable$group == group.values[g])
          single.var.idx <- which(lavpartable$op == "~~" &
            lavpartable$lhs == single.ind &
            lavpartable$rhs == single.ind &
            lavpartable$group == group.values[g])
          # user-defined residual variance
          # fixme: we take the first, in case we have multiple matches
          # (eg nlevels)
          single.var <- lavpartable$ustart[single.var.idx[1]]
          if (is.na(single.var)) {
            single.var <- 1
          }
          ov.idx <- match(single.ind, ov.names)
          if (conditional.x && nlevels == 1L) {
            ov.var <- diag(lavsamplestats@res.cov[[g]])[ov.idx]
          } else {
            ov.var <- diag(lavsamplestats@cov[[g]])[ov.idx]
          }
          # take (1 - (rvar/ov.var) * ov.var
          tmp <- (1 - (single.var / ov.var)) * ov.var
          # just in case
          if (is.na(tmp) || tmp < 0.05) {
            tmp <- 0.05
          }
          start[single.fvar.idx] <- tmp
        }
      }
    }

    # 7g) regressions "~" # new in 0.6-10
    if (length(lv.names) == 0L && nlevels == 1L && !conditional.x) {
      # observed only
      reg.idx <- which(lavpartable$group == group.values[g] &
        lavpartable$op == "~")
      if (length(reg.idx) > 0L) {
        eqs.y <- unique(lavpartable$lhs[reg.idx])
        ny <- length(eqs.y)
        for (i in seq_len(ny)) {
          y.name <- eqs.y[i]
          start.idx <- which(lavpartable$group == group.values[g] &
            lavpartable$op == "~" &
            lavpartable$lhs == y.name)
          x.names <- lavpartable$rhs[start.idx]
          COV <- lavsamplestats@cov[[g]]
          y.idx <- match(y.name, ov.names)
          x.idx <- match(x.names, ov.names)
          S.xx <- COV[x.idx, x.idx, drop = FALSE]
          S.xy <- COV[x.idx, y.idx, drop = FALSE]
          # regression coefficient(s)
          beta.i <- try(solve(S.xx, S.xy), silent = TRUE)
          if (inherits(beta.i, "try-error")) {
            start[start.idx] <- beta.i <- rep(0, length(start.idx))
          } else {
            start[start.idx] <- drop(beta.i)
          }
          # residual variance
          res.idx <- which(lavpartable$group == group.values[g] &
            lavpartable$op == "~~" &
            lavpartable$lhs == y.name &
            lavpartable$rhs == y.name)
          res.val <- COV[y.idx, y.idx] - drop(crossprod(beta.i, S.xy))
          if (res.val > 0.001 * COV[y.idx, y.idx] &&
            res.val < 0.999 * COV[y.idx, y.idx]) {
            start[res.idx] <- res.val
          } else {
            # do nothing (keep what we have)
          }
          # intercept
          int.idx <- which(lavpartable$group == group.values[g] &
            lavpartable$op == "~1" &
            lavpartable$lhs == y.name)
          if (length(int.idx) > 0L) {
            MEAN <- lavsamplestats@mean[[g]]
            Ybar <- MEAN[y.idx]
            Xbar <- MEAN[x.idx]
            int.val <- Ybar - drop(crossprod(beta.i, Xbar))
            if (is.finite(int.val)) {
              start[int.idx] <- int.val
            }
          }
        }
      }
    }

    #  # 8 latent variances (new in 0.6-2)
    #  lv.names.y <- vnames(lavpartable, "lv.y", group = group.values[g])
    #  lv.names.x <- vnames(lavpartable, "lv.x", group = group.values[g])
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
      level.values <- lav_partable_level_values(lavpartable)
      # Note: ov.names.x contains all levels within a group!
      if (length(ov.names.x) > 0) {
        for (l in 1:nlevels) {
          # block number
          block <- (g - 1L) * nlevels + l

          this.block.x <- lav_partable_vnames(lavpartable, "ov.x",
            block = block
          )
          this.block.ov <- lav_partable_vnames(lavpartable, "ov",
            block = block
          )
          if (length(this.block.x) == 0L) {
            next
          }

          # var/cov
          exo.idx <- which(lavpartable$group == group.values[g] &
            lavpartable$level == level.values[l] &
            lavpartable$op == "~~" &
            lavpartable$lhs %in% this.block.x &
            lavpartable$rhs %in% this.block.x)

          if (is.null(lavh1$implied$cov[[1]])) {
            row.idx <- match(lavpartable$lhs[exo.idx], ov.names)
            col.idx <- match(lavpartable$rhs[exo.idx], ov.names)
            if (l == 1L) {
              COV <- lavsamplestats@YLp[[g]][[2]]$S.PW.start
            } else {
              COV <- lavsamplestats@YLp[[g]][[l]]$Sigma.B
            }
          } else {
            row.idx <- match(lavpartable$lhs[exo.idx], this.block.ov)
            col.idx <- match(lavpartable$rhs[exo.idx], this.block.ov)
            COV <- lavh1$implied$cov[[block]]
          }
          # make sure starting values for variances are positive
          neg.idx <- which(diag(COV) < 0.001)
          if (length(neg.idx) > 0L) {
            diag(COV)[neg.idx] <- 0.001
          }
          start[exo.idx] <- COV[cbind(row.idx, col.idx)]

          # intercepts
          ov.int.idx <- which(lavpartable$group == group.values[g] &
            lavpartable$level == level.values[l] &
            lavpartable$op == "~1" &
            lavpartable$lhs %in% this.block.x)

          if (is.null(lavh1$implied$mean[[1]])) {
            idx <- match(lavpartable$lhs[ov.int.idx], ov.names)
            if (l == 1L) {
              INT <- lavsamplestats@YLp[[g]][[2]]$Mu.W
            } else {
              INT <- lavsamplestats@YLp[[g]][[l]]$Mu.B.start
            }
          } else {
            idx <- match(lavpartable$lhs[ov.int.idx], this.block.ov)
            INT <- lavh1$implied$mean[[block]]
          }
          start[ov.int.idx] <- INT[idx]

          # new in 0.6-12
          # very special case: conditional.x with a combination of
          # splitted-x and regular-x
          # here, we must:
          # 1) replace var/cov of splitted-x by *residual* varcov
          #    after regressing out regular-x
          # 2) replace means of splitted-x by intercepts
          # 3) fill splitted-x ~ regular-x regression coefficients
          if (conditional.x) {
            if (is.null(lavh1$implied$cov[[l]])) {
             lav_msg_stop(gettext(
               "lavh1 information is needed; please rerun with h1 = TRUE"))
            }
            blocks.within.group <- (g - 1L) * nlevels + seq_len(nlevels)
            OTHER.BLOCK.NAMES <- lav_partable_vnames(lavpartable, "ov",
                                    block = blocks.within.group[-block])
            ov.names.x.block <- this.block.x
            idx <- which(ov.names.x.block %in% OTHER.BLOCK.NAMES)
            if (length(idx) > 0L) {
              ov.names.x.block <- ov.names.x.block[-idx]
            }
            ov.names.x1 <- this.block.x[!this.block.x %in% ov.names.x.block]
            ov.names.x2 <- ov.names.x.block
            nx1 <- length(ov.names.x1) # splitted x
            nx2 <- length(ov.names.x2) # regular  x
            if (nx1 > 0L && nx2 > 0L) {
              # COV
              c1.idx <- match(ov.names.x1, this.block.ov)
              c2.idx <- match(ov.names.x2, this.block.ov)

              COV.Y <- COV[c1.idx, c1.idx, drop = FALSE]
              COV.X <- COV[c2.idx, c2.idx, drop = FALSE]
              COV.YX <- COV[c1.idx, c2.idx, drop = FALSE]
              COV.XY <- COV[c2.idx, c1.idx, drop = FALSE]
              COV.XinvYX <- solve(COV.X, COV.XY)
              RES.COV <- COV.Y - COV.YX %*% COV.XinvYX

              res.cov.idx <- which(lavpartable$group == group.values[g] &
                lavpartable$level == level.values[l] &
                lavpartable$op == "~~" &
                lavpartable$lhs %in% ov.names.x1 &
                lavpartable$rhs %in% ov.names.x1)

              row.idx <- match(lavpartable$lhs[res.cov.idx], ov.names.x1)
              col.idx <- match(lavpartable$rhs[res.cov.idx], ov.names.x1)
              start[res.cov.idx] <- RES.COV[cbind(row.idx, col.idx)]

              # INT
              INT.Y <- INT[c1.idx]
              INT.X <- INT[c2.idx]
              RES.INT <- INT.Y - t(COV.XinvYX) %*% INT.X

              res.int.idx <- which(lavpartable$group == group.values[g] &
                lavpartable$level == level.values[l] &
                lavpartable$op == "~1" &
                lavpartable$lhs %in% ov.names.x1)
              idx <- match(lavpartable$lhs[res.int.idx], ov.names.x1)
              start[res.int.idx] <- RES.INT[idx]

              # REG
              reg.idx <- which(lavpartable$group == group.values[g] &
                lavpartable$level == level.values[l] &
                lavpartable$op == "~" &
                lavpartable$lhs %in% ov.names.x1 &
                lavpartable$rhs %in% ov.names.x2)
              row.idx <- match(lavpartable$lhs[reg.idx], ov.names.x1)
              col.idx <- match(lavpartable$rhs[reg.idx], ov.names.x2)
              start[reg.idx] <- t(COV.XinvYX)[cbind(row.idx, col.idx)]
            } # special case
          } # conditional.x
        } # levels
      } # fixed.x
    } # nlevels > 1L
  } # groups




  # group weights
  group.idx <- which(lavpartable$lhs == "group" &
    lavpartable$op == "%")
  if (length(group.idx) > 0L) {
    ngroups <- length(group.idx)
    # prop <- rep(1/ngroups, ngroups)
    # use last group as reference
    # start[group.idx] <- log(prop/prop[ngroups])

    # poisson version
    start[group.idx] <- log(rep(lavsamplestats@ntotal / ngroups, ngroups))
  }

  # growth models:
  # - compute starting values for mean latent variables
  # - compute starting values for variance latent variables
  if (start.initial %in% c("lavaan", "mplus") &&
    model.type == "growth") {
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
    bad.idx <- which(start < lavpartable$lower)
    if (length(bad.idx)) {
      start[bad.idx] <- lavpartable$lower[bad.idx]
    }
  }
  if (!is.null(lavpartable$upper)) {
    bad.idx <- which(start > lavpartable$upper)
    if (length(bad.idx)) {
      start[bad.idx] <- lavpartable$upper[bad.idx]
    }
  }


  # override if the model syntax contains explicit starting values (free only)
  # user.idx <- which(!is.na(lavpartable$ustart) &
  #                  lavpartable$user != 7L) # new in 0.6-7, if rotation and
  #                                          # and we change the order of lv's
  user.idx <- which(!is.na(lavpartable$ustart) & lavpartable$free > 0L)
  start[user.idx] <- lavpartable$ustart[user.idx]

  # override if a user list with starting values is provided
  # we only look at the 'est' column for now
  if (!is.null(start.user)) {
    if (is.null(lavpartable$group)) {
      lavpartable$group <- rep(1L, length(lavpartable$lhs))
    }
    if (is.null(start.user$group)) {
      start.user$group <- rep(1L, length(start.user$lhs))
    }

    # FIXME: avoid for loop!!!
    for (i in seq_along(lavpartable$lhs)) {
      # find corresponding parameters
      lhs <- lavpartable$lhs[i]
      op <- lavpartable$op[i]
      rhs <- lavpartable$rhs[i]
      grp <- lavpartable$group[i]

      start.user.idx <- which(start.user$lhs == lhs &
        start.user$op == op &
        start.user$rhs == rhs &
        start.user$group == grp)
      if (length(start.user.idx) == 1L &&
        is.finite(start.user$est[start.user.idx])) {
        start[i] <- start.user$est[start.user.idx]
      }
    }
  }

  # override fixed values with ustart values
  user.idx <- which(!is.na(lavpartable$ustart) & lavpartable$free == 0L)
  start[user.idx] <- lavpartable$ustart[user.idx]

  # final check: no NaN or other non-finite values
  bad.idx <- which(!is.finite(start))
  if (length(bad.idx) > 0L) {
    cat("starting values:\n")
    print(start)
     lav_msg_warn(gettext(
       "some starting values are non-finite; replacing them with 0.5;
       please provide better starting values."))
    start[bad.idx] <- 0.5
  }

  if (debug) {
    cat("lavaan DEBUG: lavaanStart\n")
    print(start)
  }

  start
}

# backwards compatibility
# StartingValues <- lav_start

# sanity check: (user-specified) variances smaller than covariances
lav_start_check_cov <- function(lavpartable = NULL, start = lavpartable$start,
                                warn = TRUE) {
  nblocks <- lav_partable_nblocks(lavpartable)
  block.values <- lav_partable_block_values(lavpartable)

  for (g in 1:nblocks) {
    # collect all non-zero covariances
    cov.idx <- which(lavpartable$op == "~~" &
      lavpartable$block == block.values[g] &
      lavpartable$lhs != lavpartable$rhs &
      !lavpartable$exo &
      start != 0)

    # for each covariance, use corresponding variances to standardize;
    # the end result should not exceed abs(1)
    for (cc in seq_along(cov.idx)) {
      this.cov.idx <- cov.idx[cc]

      # find corresponding variances
      var.lhs <- lavpartable$lhs[this.cov.idx]
      var.rhs <- lavpartable$rhs[this.cov.idx]

      var.lhs.idx <- which(lavpartable$op == "~~" &
        lavpartable$block == block.values[g] &
        lavpartable$lhs == var.lhs &
        lavpartable$lhs == lavpartable$rhs)

      var.rhs.idx <- which(lavpartable$op == "~~" &
        lavpartable$block == block.values[g] &
        lavpartable$lhs == var.rhs &
        lavpartable$lhs == lavpartable$rhs)

      var.lhs.value <- start[var.lhs.idx]
      var.rhs.value <- start[var.rhs.idx]

      block.txt <- ""
      if (nblocks > 1L) {
        block.txt <- paste(" [in block ", g, "]", sep = "")
      }

      # check for zero variances
      if (var.lhs.value == 0 || var.rhs.value == 0) {
        # this can only happen if it is user-specified
        # cov.idx free? set it to zero
        if (start[this.cov.idx] == 0) {
          # nothing to do
        } else if (lavpartable$free[this.cov.idx] > 0L) {
          if (warn) {
             lav_msg_warn(gettextf(
              "non-zero covariance element set to zero, due to fixed-to-zero
              variances variables involved are: %s", var.lhs), var.rhs,
              block.txt
            )
          }
          start[this.cov.idx] <- 0
        } else {
          lav_msg_stop(gettextf(
            "please provide better fixed values for (co)variances;
            variables involved are: %s ", var.lhs), var.rhs, block.txt
          )
        }
        next
      }

      # which one is the smallest? abs() in case of negative variances
      if (abs(var.lhs.value) < abs(var.rhs.value)) {
        var.min.idx <- var.lhs.idx
        var.max.idx <- var.rhs.idx
      } else {
        var.min.idx <- var.rhs.idx
        var.max.idx <- var.lhs.idx
      }

      # check
      COR <- abs(start[this.cov.idx] / sqrt(var.lhs.value * var.rhs.value))

      # NOTE: we treat this as an unconditional COR!

      if (!is.finite(COR)) {
        # force simple values
        if (warn) {
           lav_msg_warn(gettextf(
            "starting values imply NaN for a correlation value; variables
            involved are: %s", var.lhs), var.rhs, block.txt
          )
        }
        start[var.lhs.idx] <- 1
        start[var.rhs.idx] <- 1
        start[this.cov.idx] <- 0
      } else if (COR > 1) {
        txt <- gettextf(
          "starting values imply a correlation larger than 1; variables
          involved are: %1$s %2$s %3$s", var.lhs, var.rhs, block.txt)

        # three ways to fix it: rescale cov12, var1 or var2

        # we prefer a free parameter, and not user-specified
        if (lavpartable$free[this.cov.idx] > 0L &&
          is.na(lavpartable$ustart[this.cov.idx])) {
          if (warn) {
             lav_msg_warn(gettext(txt))
          }
          start[this.cov.idx] <- start[this.cov.idx] / (COR * 1.1)
        } else if (lavpartable$free[var.min.idx] > 0L &&
          is.na(lavpartable$ustart[var.min.idx])) {
          if (warn) {
             lav_msg_warn(gettext(txt))
          }
          start[var.min.idx] <- start[var.min.idx] * (COR * 1.1)^2
        } else if (lavpartable$free[var.max.idx] > 0L &&
          is.na(lavpartable$ustart[var.max.idx])) {
          if (warn) {
             lav_msg_warn(gettext(txt))
          }
          start[var.max.idx] <- start[var.max.idx] * (COR * 1.1)^2

          # not found? try just a free parameter
        } else if (lavpartable$free[this.cov.idx] > 0L) {
          if (warn) {
             lav_msg_warn(gettext(txt))
          }
          start[this.cov.idx] <- start[this.cov.idx] / (COR * 1.1)
        } else if (lavpartable$free[var.min.idx] > 0L) {
          if (warn) {
             lav_msg_warn(gettext(txt))
          }
          start[var.min.idx] <- start[var.min.idx] * (COR * 1.1)^2
        } else if (lavpartable$free[var.max.idx] > 0L) {
          if (warn) {
             lav_msg_warn(gettext(txt))
          }
          start[var.max.idx] <- start[var.max.idx] * (COR * 1.1)^2

          # nothing? abort or warn (and fail later...): warn
        } else {
          if (warn) {
             lav_msg_warn(gettext(txt))
          }
          # lav_msg_stop(gettextf(
          # "please provide better fixed values for (co)variances;
          #  variables involved are: %s ", var.lhs), var.rhs, block.txt)
        }
      } # COR > 1
    } # cov.idx
  }
  start
}
