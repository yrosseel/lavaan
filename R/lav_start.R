# lav_start.R: provide starting values for model parameters
#
# YR 30/11/2010: initial version
# YR 08/06/2011: add fabin3 start values for factor loadings
# YR 14 Jan 2014: moved to lav_start.R

# fill in the 'ustart' column in a User data.frame with reasonable
# starting values, using the sample data

lav_start <- function(start.method    = "default",
                      lavpartable     = NULL,
                      lavsamplestats  = NULL,
                      model.type      = "sem",
                      mimic           = "lavaan",
                      debug           = FALSE) {

    # check arguments
    stopifnot(is.list(lavpartable))

    # categorical?
    categorical <- any(lavpartable$op == "|")

    # conditional.x?
    conditional.x <- any(lavpartable$exo == 1L & lavpartable$op == "~")
    #ord.names <- unique(lavpartable$lhs[ lavpartable$op == "|" ])

    # nlevels?
    nlevels <- lav_partable_nlevels(lavpartable)

    # shortcut for 'simple'
    # be we should also take care of the 'fixed.x' issue
    if(identical(start.method, "simple")) {
        start <- numeric( length(lavpartable$ustart) )
        start[ which(lavpartable$op == "=~") ] <- 1.0
        start[ which(lavpartable$op == "~*~") ] <- 1.0
        ov.names.ord <- vnames(lavpartable, "ov.ord")
        var.idx <- which(lavpartable$op == "~~" & lavpartable$lhs == lavpartable$rhs &
                         !(lavpartable$lhs %in% ov.names.ord))
        start[var.idx] <- 1.0
        user.idx <- which(!is.na(lavpartable$ustart))
        start[user.idx] <- lavpartable$ustart[user.idx]
        return(start)
    }

    # check start.method
    if(mimic == "lavaan") {
        start.initial <- "lavaan"
    } else if(mimic == "Mplus") {
        start.initial <- "mplus"
    } else {
        # FIXME: use LISREL/EQS/AMOS/.... schems
        start.initial <- "lavaan"
    }
    start.user    <- NULL
    if(is.character(start.method)) {
        start.method. <- tolower(start.method)
        if(start.method. == "default") {
            # nothing to do
        } else if(start.method. %in% c("simple", "lavaan", "mplus")) {
            start.initial <- start.method.
        } else {
            stop("lavaan ERROR: unknown value for start argument")
        }
    } else if(is.list(start.method)) {
        start.user <- start.method
    } else if(inherits(start.method, "lavaan")) {
        start.user <- parTable(start.method)
    }
    # check model list elements, if provided
    if(!is.null(start.user)) {
        if(is.null(start.user$lhs) ||
           is.null(start.user$op)  ||
           is.null(start.user$rhs)) {
            stop("lavaan ERROR: problem with start argument: model list does not contain all elements: lhs/op/rhs")
        }
        if(!is.null(start.user$est)) {
            # excellent, we got an est column; nothing to do
        } else if(!is.null(start.user$start)) {
            # no est column, but we use the start column
            start.user$est <- start.user$start
        } else if(!is.null(start.user$ustart)) {
            # no ideal, but better than nothing
            start.user$est <- start.user$ustart
        } else {
            stop("lavaan ERROR: problem with start argument: could not find est/start column in model list")
        }
    }


    # global settings
    # 0. everyting is zero
    start <- numeric( length(lavpartable$ustart) )

    # 1. =~ factor loadings:
    if(categorical) {
        # if std.lv=TRUE, more likely initial Sigma.hat is positive definite
        # 0.8 is too large
        start[ which(lavpartable$op == "=~") ] <- 0.7
    } else {
        start[ which(lavpartable$op == "=~") ] <- 1.0
    }

    # 2. (residual) lv variances for latent variables
    lv.names    <- vnames(lavpartable, "lv") # all groups
    lv.var.idx <- which(lavpartable$op == "~~"        &
                        lavpartable$lhs %in% lv.names &
                        lavpartable$lhs == lavpartable$rhs)
    start[lv.var.idx] <- 0.05
    #start[lv.var.idx] <- 0.5 # new in 0.6-2? (for optim.parscale = "stand")

    # 3. latent response scales (if any)
    delta.idx <- which(lavpartable$op == "~*~")
    start[delta.idx] <- 1.0


    # group-specific settings
    ngroups <- lav_partable_ngroups(lavpartable)

    # for now, if no group column, add one (again), until we rewrite
    # this function to handle block/group hybrid settings
    if(is.null(lavpartable$group) && ngroups == 1L) {
        lavpartable$group <- rep(1L, length(lavpartable$lhs))
        lavpartable$group[ lavpartable$block == 0L] <- 0L
    }


    for(g in 1:ngroups) {

        # group values (not necessarily 1,2,... anymore)
        group.values <- lav_partable_group_values(lavpartable)

        # info from user model for this group
        if(conditional.x) {
            ov.names     <- vnames(lavpartable, "ov.nox", group = group.values[g])
        } else {
            ov.names     <- vnames(lavpartable, "ov", group = group.values[g])
        }
        if(categorical) {
            ov.names.num <- vnames(lavpartable, "ov.num", group = group.values[g])
            ov.names.ord <- vnames(lavpartable, "ov.ord", group = group.values[g])
        } else {
            ov.names.num <- ov.names
        }
        lv.names    <- vnames(lavpartable, "lv",   group = group.values[g])
        ov.names.x  <- vnames(lavpartable, "ov.x", group = group.values[g])

        # just for the nlevels >1 case
        ov.names <- unique(unlist(ov.names))
        ov.names.num <- unique(unlist(ov.names.num))
        lv.names <- unique(unlist(lv.names))
        ov.names.x <- unique(unlist(ov.names.x))


        # residual ov variances (including exo/ind, to be overriden)
        ov.var.idx <- which(lavpartable$group == group.values[g]             &
                            lavpartable$op    == "~~"          &
                            lavpartable$lhs %in% ov.names.num  &
                            lavpartable$lhs == lavpartable$rhs)
        sample.var.idx <- match(lavpartable$lhs[ov.var.idx], ov.names)
        if(model.type == "unrestricted") {
            start[ov.var.idx] <- diag(lavsamplestats@cov[[g]])[sample.var.idx]
        } else {
            if(start.initial == "mplus") {
                if(conditional.x) {
                    start[ov.var.idx] <-
                       (1.0 - 0.50)*lavsamplestats@res.var[[1L]][sample.var.idx]
                } else {
                    start[ov.var.idx] <-
                       (1.0 - 0.50)*lavsamplestats@var[[1L]][sample.var.idx]
                }
            } else {
                if(conditional.x) {
                    start[ov.var.idx] <-
                  (1.0 - 0.50)*diag(lavsamplestats@res.cov[[g]])[sample.var.idx]
                } else {
                    start[ov.var.idx] <-
                    (1.0 - 0.50)*diag(lavsamplestats@cov[[g]])[sample.var.idx]
                }
            }
        }

        # 1-fac measurement models: loadings, psi, theta
        if(start.initial %in% c("lavaan", "mplus") &&
           model.type %in% c("sem", "cfa") ) {
            # fabin3 estimator (2sls) of Hagglund (1982) per factor
            for(f in lv.names) {
                lambda.idx <- which( lavpartable$lhs == f &
                                     lavpartable$op == "=~" &
                                     lavpartable$group == group.values[g] )
                # standardized?
                std.lv <- FALSE
                var.f.idx <- which(lavpartable$lhs == f &
                                   lavpartable$op == "~~" &
                                   lavpartable$group == group.values[g] &
                                   lavpartable$rhs == f)
                if(length(var.f.idx) > 0L &&
                   lavpartable$free[var.f.idx] == 0 &&
                   lavpartable$ustart[var.f.idx] == 1) {
                    std.lv <- TRUE
                }

                # no second order
                if(any(lavpartable$rhs[lambda.idx] %in% lv.names)) next

                # get observed indicators for this latent variable
                ov.idx <- match(lavpartable$rhs[lambda.idx], ov.names)
                if(length(ov.idx) > 0L && !any(is.na(ov.idx))) {
                    if(lavsamplestats@missing.flag && nlevels == 1L) {
                        COV <- lavsamplestats@missing.h1[[g]]$sigma[ov.idx,
                                                      ov.idx, drop = FALSE]
                    } else {
                        if(conditional.x) {
                            COV <- lavsamplestats@res.cov[[g]][ov.idx,
                                                      ov.idx, drop = FALSE]
                        } else {
                            COV <- lavsamplestats@cov[[g]][ov.idx,
                                                      ov.idx, drop = FALSE]
                        }
                    }

                    # fabin for 1-factor
                    fabin <- lav_cfa_1fac_fabin(COV, std.lv = std.lv,
                                                lambda.only = TRUE,
                                                method = "fabin3")

                    # factor loadings
                    start[lambda.idx] <- fabin$lambda

                    # factor variance
                    #if(!std.lv) {
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
                    #}

                    # NOTE: fabin (sometimes) gives residual variances
                    # that are larger than the original variances...

                    # residual variances -- order?
                    #res.idx <- which(lavpartable$lhs %in% ov.names[ov.idx] &
                    #                 lavpartable$op == "~~" &
                    #                 lavpartable$group == group.values[g] &
                    #                 lavpartable$rhs == lavpartable$lhs)
                    #start[res.idx] <- fabin$theta

                    # negative variances?
                    #neg.idx <- which(start[res.idx] < 0)
                    #if(length(neg.idx) > 0L) {
                    #    start[res.idx][neg.idx] <- 0.05
                    #}
                }
            }
        } # fabin

        if(model.type == "unrestricted") {
           # fill in 'covariances' from lavsamplestats
            cov.idx <- which(lavpartable$group == group.values[g]             &
                             lavpartable$op    == "~~"          &
                             lavpartable$lhs != lavpartable$rhs)
            lhs.idx <- match(lavpartable$lhs[cov.idx], ov.names)
            rhs.idx <- match(lavpartable$rhs[cov.idx], ov.names)
            start[cov.idx] <- lavsamplestats@cov[[g]][ cbind(lhs.idx, rhs.idx) ]
        }

        # variances of ordinal variables - set to 1.0
        if(categorical) {
            ov.var.ord.idx <- which(lavpartable$group == group.values[g]            &
                                    lavpartable$op    == "~~"         &
                                    lavpartable$lhs %in% ov.names.ord &
                                    lavpartable$lhs == lavpartable$rhs)
            start[ov.var.ord.idx] <- 1.0
        }

        # 3g) intercepts/means
        ov.int.idx <- which(lavpartable$group == group.values[g]         &
                            lavpartable$op == "~1"         &
                            lavpartable$lhs %in% ov.names)
        sample.int.idx <- match(lavpartable$lhs[ov.int.idx], ov.names)
        if(lavsamplestats@missing.flag && nlevels == 1L) {
            start[ov.int.idx] <- lavsamplestats@missing.h1[[g]]$mu[sample.int.idx]
        } else {
            if(conditional.x) {
                start[ov.int.idx] <- lavsamplestats@res.int[[g]][sample.int.idx]
            } else {
                start[ov.int.idx] <- lavsamplestats@mean[[g]][sample.int.idx]
            }
        }

        # 4g) thresholds
        th.idx <- which(lavpartable$group == group.values[g] & lavpartable$op == "|")
        if(length(th.idx) > 0L) {
            th.names.lavpartable <- paste(lavpartable$lhs[th.idx], "|",
                                       lavpartable$rhs[th.idx], sep="")
            th.names.sample   <-
                lavsamplestats@th.names[[g]][ lavsamplestats@th.idx[[g]] > 0L ]
            # th.names.sample should identical to
            # vnames(lavpartable, "th", group = group.values[g])
            if(conditional.x) {
                th.values <-
                 lavsamplestats@res.th[[g]][lavsamplestats@th.idx[[g]] > 0L]
            } else {
                th.values <-
                  lavsamplestats@th[[g]][lavsamplestats@th.idx[[g]] > 0L]
            }
            start[th.idx] <- th.values[match(th.names.lavpartable,
                                             th.names.sample)]
        }

        # 5g) exogenous `fixed.x' covariates
        if(length(ov.names.x) > 0) {
            exo.idx <- which(lavpartable$group == group.values[g]          &
                             lavpartable$op == "~~"          &
                             lavpartable$lhs %in% ov.names.x &
                             lavpartable$rhs %in% ov.names.x)
            if(!conditional.x) {
                row.idx <- match(lavpartable$lhs[exo.idx], ov.names)
                col.idx <- match(lavpartable$rhs[exo.idx], ov.names)
                if(lavsamplestats@missing.flag && nlevels == 1L) {
                    start[exo.idx] <-
                    lavsamplestats@missing.h1[[g]]$sigma[cbind(row.idx,col.idx)]
                    # using slightly smaller starting values for free
                    # variance/covariances (fixed.x = FALSE);
                    # this somehow avoids false convergence in saturated models
                    nobs <- lavsamplestats@nobs[[g]]
                    this.idx <- which( seq_len(length(lavpartable$free)) %in% exo.idx & lavpartable$free > 0L )
                    start[this.idx] <- start[this.idx] * (nobs-1)/nobs
                } else {
                    start[exo.idx] <- lavsamplestats@cov[[g]][cbind(row.idx,col.idx)]
                }
            } else {
                # cov.x
                row.idx <- match(lavpartable$lhs[exo.idx], ov.names.x)
                col.idx <- match(lavpartable$rhs[exo.idx], ov.names.x)
                start[exo.idx] <- lavsamplestats@cov.x[[g]][cbind(row.idx,
                                                                  col.idx)]
                # mean.x
                exo.int.idx <- which(lavpartable$group == group.values[g] &
                                     lavpartable$op == "~1"               &
                                     lavpartable$lhs %in% ov.names.x)
                int.idx <- match(lavpartable$lhs[exo.int.idx], ov.names.x)
                start[exo.int.idx] <- lavsamplestats@mean.x[[g]][int.idx]
            }
        }

        # 6b. exogenous lv variances if single indicator -- new in 0.5-21
        lv.x <- vnames(lavpartable, "lv.x", group = group.values[g])
        # FIXME: also for multilevel?
        lv.x <- unique(unlist(lv.x))
        if(length(lv.x) > 0L) {
            for(ll in lv.x) {
                ind.idx <- which(lavpartable$op == "=~" &
                                 lavpartable$lhs == ll,
                                 lavpartable$group == group.values[g])
                if(length(ind.idx) == 1L) {
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
                    if(is.na(single.var)) {
                         single.var <- 1
                    }
                    ov.idx <- match(single.ind, ov.names)
                    if(conditional.x) {
                        ov.var <- diag(lavsamplestats@res.cov[[g]])[ov.idx]
                    } else {
                        ov.var <- diag(lavsamplestats@cov[[g]])[ov.idx]
                    }
                    # take (1 - (rvar/ov.var) * ov.var
                    tmp <- (1 - (single.var/ov.var)) * ov.var
                    # just in case
                    if(is.na(tmp) || tmp < 0.05) {
                        tmp <- 0.05
                    }
                    start[single.fvar.idx] <- tmp
                }
            }
        }

        # 7g) regressions "~"

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

    } # groups


    # nlevels > 1L
    if(nlevels > 1L) {
        for(g in 1:ngroups) {
            group.values <- lav_partable_group_values(lavpartable)
            ov.names.x  <- vnames(lavpartable, "ov.x", group = group.values[g])
            ov.names.x <- unique(unlist(ov.names.x))
            level.values <- lav_partable_level_values(lavpartable)

            if(!conditional.x && length(ov.names.x) > 0) {
                for(l in 1:nlevels) {
                    # var/cov
                    exo.idx <- which(lavpartable$group == group.values[g] &
                                     lavpartable$level == level.values[l] &
                                     lavpartable$op == "~~"               &
                                     lavpartable$lhs %in% ov.names.x      &
                                     lavpartable$rhs %in% ov.names.x)
                    row.idx <- match(lavpartable$lhs[exo.idx], ov.names)
                    col.idx <- match(lavpartable$rhs[exo.idx], ov.names)

                    if(l == 1L) {
                        COV <- lavsamplestats@YLp[[g]][[2]]$S.PW.start
                    } else {
                        COV <- lavsamplestats@YLp[[g]][[l]]$Sigma.B
                    }
                    start[exo.idx] <- COV[ cbind(row.idx, col.idx) ]

                    # intercepts
                    ov.int.idx <- which(lavpartable$group == group.values[g] &
                                        lavpartable$level == level.values[l] &
                                        lavpartable$op == "~1"               &
                                        lavpartable$lhs %in% ov.names.x)
                    idx <- match(lavpartable$lhs[ov.int.idx], ov.names)

                    if(l == 1L) {
                        INT <- lavsamplestats@YLp[[g]][[2]]$Mu.W
                    } else {
                        INT <- lavsamplestats@YLp[[g]][[l]]$Mu.B.start
                    }

                    start[ov.int.idx] <- INT[idx]
                } # levels
            } # fixed.x
        } # groups
    } # nlevels > 1L


    # group weights
    group.idx <- which(lavpartable$lhs == "group" &
                       lavpartable$op  == "%")
    if(length(group.idx) > 0L) {
        ngroups <- length(group.idx)
        #prop <- rep(1/ngroups, ngroups)
        # use last group as reference
        #start[group.idx] <- log(prop/prop[ngroups])

        # poisson version
        start[group.idx] <- log( rep(lavsamplestats@ntotal/ngroups, ngroups) )
    }

    # growth models:
    # - compute starting values for mean latent variables
    # - compute starting values for variance latent variables
    if(start.initial %in% c("lavaan", "mplus") &&
       model.type == "growth") {
        ### DEBUG ONLY
        #lv.var.idx <- which(lavpartable$op == "~~"                &
        #                lavpartable$lhs %in% lv.names &
        #                lavpartable$lhs == lavpartable$rhs)

        ### DEBUG ONLY
        #lv.int.idx <- which(lavpartable$op == "~1"         &
        #                    lavpartable$lhs %in% lv.names)
    }

    # override if a user list with starting values is provided
    # we only look at the 'est' column for now
    if(!is.null(start.user)) {

        if(is.null(lavpartable$group)) {
            lavpartable$group <- rep(1L, length(lavpartable$lhs))
        }
        if(is.null(start.user$group)) {
            start.user$group <- rep(1L, length(start.user$lhs))
        }

        # FIXME: avoid for loop!!!
        for(i in 1:length(lavpartable$lhs)) {
            # find corresponding parameters
            lhs <- lavpartable$lhs[i]
             op <- lavpartable$op[i]
            rhs <- lavpartable$rhs[i]
            grp <- lavpartable$group[i]

            start.user.idx <- which(start.user$lhs == lhs &
                                    start.user$op  ==  op &
                                    start.user$rhs == rhs &
                                    start.user$group == grp)
            if(length(start.user.idx) == 1L &&
               is.finite(start.user$est[start.user.idx])) {
                start[i] <- start.user$est[start.user.idx]
            }
        }
    }

    # override if the model syntax contains explicit starting values
    user.idx <- which(!is.na(lavpartable$ustart))
    start[user.idx] <- lavpartable$ustart[user.idx]

    if(debug) {
        cat("lavaan DEBUG: lavaanStart\n")
        print( start )
    }

    start
}

# backwards compatibility
# StartingValues <- lav_start

# sanity check: (user-specified) variances smaller than covariances
lav_start_check_cov <- function(lavpartable = NULL, start = lavpartable$start) {

    nblocks <- lav_partable_nblocks(lavpartable)

    for(g in 1:nblocks) {

        # block values
        block.values <- lav_partable_block_values(lavpartable)

        # collect all non-zero covariances
        cov.idx <- which(lavpartable$op == "~~" &
                         lavpartable$block == block.values[g] &
                         lavpartable$lhs != lavpartable$rhs &
                         !lavpartable$exo &
                         start != 0)

        # for each covariance, use corresponding variances to standardize;
        # the end result should not exceed abs(1)
        for(cc in seq_along(cov.idx)) {
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
            if(nblocks > 1L) {
                block.txt <- paste(" [in block ", g, "]", sep = "")
            }

            # check for zero variances
            if(var.lhs.value == 0 || var.rhs.value == 0) {
                # this can only happen if it is user-specified
                # cov.idx free? set it to zero
                if(start[this.cov.idx] == 0) {
                    # nothing to do
                } else if(lavpartable$free[this.cov.idx] > 0L) {
                    warning(
  "lavaan WARNING: non-zero covariance element set to zero, due to fixed-to-zero variances\n",
"                  variables involved are: ", var.lhs, " ", var.rhs, block.txt)
                    start[this.cov.idx] <- 0
                } else {
                    stop("lavaan ERROR: please provide better fixed values for (co)variances;\n",
"                variables involved are: ", var.lhs, " ", var.rhs, block.txt)
                }
                next
            }

            # which one is the smallest? abs() in case of negative variances
            if(abs(var.lhs.value) < abs(var.rhs.value)) {
                var.min.idx <- var.lhs.idx
                var.max.idx <- var.rhs.idx
            } else {
                var.min.idx <- var.rhs.idx
                var.max.idx <- var.lhs.idx
            }

            # check
            COR <- start[this.cov.idx] / sqrt(var.lhs.value * var.rhs.value)

            if(!is.finite(COR)) {
                # force simple values
                warning(
  "lavaan WARNING: starting values imply NaN for a correlation value;\n",
"                  variables involved are: ", var.lhs, " ", var.rhs, block.txt)
                start[var.lhs.idx] <- 1
                start[var.rhs.idx] <- 1
                start[this.cov.idx] <- 0
            } else if(abs(COR) > 1) {
                warning(
  "lavaan WARNING: starting values imply a correlation larger than 1;\n",
"                  variables involved are: ", var.lhs, " ", var.rhs, block.txt)

                # three ways to fix it: rescale cov12, var1 or var2

                # we prefer a free parameter, and not user-specified
                if(        lavpartable$free[this.cov.idx] > 0L &&
                   is.na(lavpartable$ustart[this.cov.idx])) {
                    start[this.cov.idx] <- start[this.cov.idx] / (COR * 1.1)
                } else if( lavpartable$free[var.min.idx] > 0L &&
                           is.na(lavpartable$ustart[var.min.idx])) {
                    start[var.min.idx] <- start[var.min.idx] * (COR * 1.1)^2
                } else if(        lavpartable$free[var.max.idx] > 0L &&
                          is.na(lavpartable$ustart[var.max.idx])) {
                    start[var.max.idx] <- start[var.max.idx] * (COR * 1.1)^2

                # not found? try just a free parameter
                } else if      (lavpartable$free[this.cov.idx] > 0L) {
                    start[this.cov.idx] <- start[this.cov.idx] / (COR * 1.1)
                } else if(     lavpartable$free[var.min.idx] > 0L) {
                    start[var.min.idx] <- start[var.min.idx] * (COR * 1.1)^2
                } else if(     lavpartable$free[var.max.idx] > 0L) {
                    start[var.max.idx] <- start[var.max.idx] * (COR * 1.1)^2

                # nothing? abort
                } else {
                    stop("lavaan ERROR: please provide better fixed values for (co)variances;\n",
"                variables involved are: ", var.lhs, " ", var.rhs, block.txt)
                }
            } # COR > 1
        } # cov.idx
    }

    start
}
