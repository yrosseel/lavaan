# compute sample statistics for the unrestricted (h1) model
# and also the logl (if available)
lav_h1_implied_logl <- function(lavdata = NULL,
                                lavsamplestats = NULL,
                                lavpartable = NULL, # multilevel + missing
                                lavoptions = NULL) {
  lavpta <- NULL
  if (!is.null(lavpartable)) {
    lavpta <- lav_partable_attributes(lavpartable)
    lavpartable <- lav_partable_set_cache(lavpartable, lavpta)
  }
  if (lavdata@nlevels == 1L) {
    if (lavsamplestats@missing.flag) {
      if (lavoptions$conditional.x) {
        implied <- list() # not available yet
      } else {
        implied <- list(
          cov = lapply(
            lavsamplestats@missing.h1,
            "[[", "sigma"
          ),
          mean = lapply(
            lavsamplestats@missing.h1,
            "[[", "mu"
          ),
          th = lavsamplestats@th,
          group.w = lavsamplestats@group.w
        )
      }
    } else {
      if (lavoptions$conditional.x) {
        implied <- list(
          res.cov = lavsamplestats@res.cov,
          res.int = lavsamplestats@res.int,
          res.slopes = lavsamplestats@res.slopes,
          cov.x = lavsamplestats@cov.x,
          mean.x = lavsamplestats@mean.x,
          res.th = lavsamplestats@res.th,
          group.w = lavsamplestats@group.w
        )
      } else {
        implied <- list(
          cov = lavsamplestats@cov,
          mean = lavsamplestats@mean,
          th = lavsamplestats@th,
          group.w = lavsamplestats@group.w
        )
      }
    } # complete data

    logl <- lav_h1_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  } else {
    # estimate Mu.B, Mu.W, Sigma.B and Sigma.W for unrestricted model
    ngroups <- lavdata@ngroups
    nlevels <- lavdata@nlevels
    implied <- list(
      cov = vector("list", length = ngroups * nlevels),
      mean = vector("list", length = ngroups * nlevels)
    )
    loglik.group <- numeric(lavdata@ngroups)

    for (g in 1:lavdata@ngroups) {
      if (lavoptions$verbose) {
        cat("\n\nfitting unrestricted (H1) model in group ", g, "\n")
      }
      if (lavsamplestats@missing.flag) {
        # missing data

        # 1. first a few EM iteration faking complete data
        # Y1 <- lavdata@X[[g]]
        # cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
        # Y2.complete <- unname(as.matrix(aggregate(Y1,
        #                by = list(cluster.idx),
        #                FUN = function(x) {
        #                    if( all(is.na(x)) ) { # all elements are NA
        #                        as.numeric(0)     # in this cluster
        #                    } else {
        #                        mean(x, na.rm = TRUE)
        #                    }
        #                })[,-1]))
        # YLp = lavsamplestats@YLp[[g]]
        # YLp[[2]]$Y2 <- Y2.complete
        # OUT <- lav_mvnorm_cluster_em_sat(
        #    YLp      = YLp,
        #    Lp       = lavdata@Lp[[g]],
        #    verbose  = TRUE, # for now
        #    tol      = 1e-04,
        #    min.variance = 1e-05,
        #    max.iter = 5L)
        ## create tmp lav1, only for this group
        # implied$cov[[ (g-1)*nlevels + 1L]] <- OUT$Sigma.W
        # implied$cov[[ (g-1)*nlevels + 2L]] <- OUT$Sigma.B
        # implied$mean[[(g-1)*nlevels + 1L]] <- OUT$Mu.W
        # implied$mean[[(g-1)*nlevels + 2L]] <- OUT$Mu.B
        # loglik.group[g] <- OUT$logl
        # lavh1 <- list(implied = implied, logl = sum(loglik.group))
        # lavpartable <- lav_partable_unrestricted(lavdata = lavdata,
        #     lavsamplestats = lavsamplestats, lavoptions = lavoptions,
        #     lavpta = lavpta, lavh1 = lavh1)
        # lavpartable$lower <- rep(-Inf, length(lavpartable$lhs))
        # var.idx <- which(lavpartable$free > 0L &
        #                 lavpartable$op == "~~" &
        #                 lavpartable$lhs == lavpartable$rhs)
        # lavpartable$lower[var.idx] <- 1e-05

        lavpartable <- lav_partable_unrestricted_chol(
          lavdata = lavdata,
          lavoptions = lavoptions, lavpta = lavpta
        )

        lavoptions2 <- lavoptions
        lavoptions2$se <- "none"
        lavoptions2$test <- "none"
        lavoptions2$do.fit <- TRUE
        # lavoptions2$verbose <- FALSE
        lavoptions2$h1 <- FALSE
        lavoptions2$baseline <- FALSE
        lavoptions2$fixed.x <- FALSE # even if model uses fixed.x=TRUE
        lavoptions2$model.type <- "unrestricted"
        lavoptions2$optim.attempts <- 4L
        lavoptions2$warn <- FALSE
        lavoptions2$check.gradient <- FALSE
        lavoptions2$optim.force.convergence <- TRUE # for now...
        lavoptions2$control <- list(rel.tol = 1e-7)
        # FIT <- lavaan(lavpartable, slotOptions = lavoptions2,
        #          slotSampleStats = lavsamplestats,
        #          slotData = lavdata, sloth1 = lavh1)
        FIT <- lavaan(lavpartable,
          slotOptions = lavoptions2,
          slotSampleStats = lavsamplestats,
          slotData = lavdata
        )
        OUT <- list(
          Sigma.W = FIT@implied$cov[[1]],
          Sigma.B = FIT@implied$cov[[2]],
          Mu.W = FIT@implied$mean[[1]],
          Mu.B = FIT@implied$mean[[2]],
          logl = FIT@loglik$loglik
        )
        # if(lavoptions$fixed.x) {
        #   OUT$logl <- OUT$logl - lavsamplestats@YLp[[g]][[2]]$loglik.x
        # }
      } else {
        # complete data
        OUT <- lav_mvnorm_cluster_em_sat(
          YLp = lavsamplestats@YLp[[g]],
          Lp = lavdata@Lp[[g]],
          verbose = lavoptions$verbose,
          tol = 1e-04, # option?
          min.variance = 1e-05, # option?
          max.iter = 5000L
        ) # option?
      }
      if (lavoptions$verbose) {
        cat("\n")
      }

      # if any near-zero within variance(s), produce warning here
      zero.var <- which(diag(OUT$Sigma.W) <= 1e-05)
      if (length(zero.var)) {
        gtxt <- if (ngroups > 1L) {
          gettextf(" in group %s.", g)
        } else {
          " "
        }
        lav_msg_warn(gettextf(
          "H1 estimation resulted in a within covariance matrix %1$s with
          (near) zero variances for some of the level-1 variables: %2$s",
            gtxt, lav_msg_view(lavdata@ov.names.l[[g]][[1]][zero.var]))
        )
      }

      # store in implied
      implied$cov[[(g - 1) * nlevels + 1L]] <- OUT$Sigma.W
      implied$cov[[(g - 1) * nlevels + 2L]] <- OUT$Sigma.B
      implied$mean[[(g - 1) * nlevels + 1L]] <- OUT$Mu.W
      implied$mean[[(g - 1) * nlevels + 2L]] <- OUT$Mu.B

      # store logl per group
      loglik.group[g] <- OUT$logl
    }

    logl <- list(loglik = sum(loglik.group), loglik.group = loglik.group)
  }

  list(implied = implied, logl = logl)
}
