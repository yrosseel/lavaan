# compute sample statistics for the unrestricted (h1) model
# and also the logl (if available)
lav_h1_implied_logl <- function(lavdata = NULL,
                                lavsamplestats = NULL,
                                lavpartable = NULL, # multilevel + missing
                                lavoptions = NULL) {
  lavpta <- NULL
  if (!is.null(lavpartable)) {
    lavpta <- lav_pt_attributes(lavpartable)
    lavpartable <- lav_pt_set_cache(lavpartable, lavpta)
  }

  # single-level case
  if (lavdata@nlevels == 1L) {

    # missing data
    if (lavsamplestats@missing.flag) {
      if (lavoptions$conditional.x) {
        implied <- list() # not available yet
      } else {
        implied <- list(
          cov = lavsamplestats@cov,
          mean = lavsamplestats@mean,
          th = lavsamplestats@th,
          group.w = lavsamplestats@group.w
        )
        # insert ML estimates
        for (g in 1:lavdata@ngroups) {
          # zero coverage?
          if (any(lav_mat_vech(lavdata@Mp[[g]]$coverage,
                                  diagonal = FALSE) == 0L)) {
            out <- lav_mvn_mi_h1_est_moments_chol(
              lavdata = lavdata, lavsamplestats = lavsamplestats,
              lavoptions = lavoptions, group = g)
            implied$cov[[g]] <- out$Sigma
            implied$mean[[g]] <- out$Mu
          } else {
            # regular EM estimates
            implied$cov[[g]] <- lavsamplestats@missing.h1[[g]]$sigma
            implied$mean[[g]] <- lavsamplestats@missing.h1[[g]]$mu
          }
        }
      }

    # complete data
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
      h1_implied = implied,
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
    loglik_group <- numeric(lavdata@ngroups)

    for (g in 1:lavdata@ngroups) {
      if (lav_verbose()) {
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
        # OUT <- lav_mvn_cl_em_sat(
        #    YLp      = YLp,
        #    Lp       = lavdata@Lp[[g]],
        #    verbose  = TRUE, # for now
        #    tol      = 1e-04,
        #    min.variance = 1e-05,
        #    max.iter = 5L)
        ## create tmp lav1, only for this group
        # implied$cov[[(g-1)*nlevels + 1L]] <- OUT$Sigma.W
        # implied$cov[[(g-1)*nlevels + 2L]] <- OUT$Sigma.B
        # implied$mean[[(g-1)*nlevels + 1L]] <- OUT$Mu.W
        # implied$mean[[(g-1)*nlevels + 2L]] <- OUT$Mu.B
        # loglik.group[g] <- OUT$logl
        # lavh1 <- list(implied = implied, logl = sum(loglik.group))
        # lavpartable <- lav_pt_unrestricted(lavdata = lavdata,
        #     lavsamplestats = lavsamplestats, lavoptions = lavoptions,
        #     lavpta = lavpta, lavh1 = lavh1)
        # lavpartable$lower <- rep(-Inf, length(lavpartable$lhs))
        # var.idx <- which(lavpartable$free > 0L &
        #                 lavpartable$op == "~~" &
        #                 lavpartable$lhs == lavpartable$rhs)
        # lavpartable$lower[var.idx] <- 1e-05

        lavpartable <- lav_pt_unrestricted_chol(
          lavdata = lavdata,
          lavoptions = lavoptions, lavpta = lavpta
        )

        lavoptions2 <- lavoptions
    lavoptions2$estimator <- "ML"
        lavoptions2$se <- "none"
        lavoptions2$test <- "none"
        lavoptions2$do.fit <- TRUE
    lavoptions2$optim.method <- "nlminb"
        lavoptions2$h1 <- FALSE
    lavoptions2$implied <- TRUE
    lavoptions2$loglik <- TRUE
        lavoptions2$baseline <- FALSE
        lavoptions2$fixed.x <- FALSE # even if model uses fixed.x=TRUE
        lavoptions2$model.type <- "unrestricted"
        lavoptions2$optim.attempts <- 4L
        lavoptions2$check.gradient <- FALSE
        lavoptions2$optim.force.converged <- TRUE # for now...
        lavoptions2$control <- list(rel.tol = 1e-7)
        # FIT <- lavaan(lavpartable, slot_options = lavoptions2,
        #          slot_sample_stats = lavsamplestats,
        #          slot_data = lavdata, sloth1 = lavh1)
        fit <- lavaan(lavpartable,
          slot_options = lavoptions2,
          slot_sample_stats = lavsamplestats,
          slot_data = lavdata,
          warn = FALSE
        )
        out_1 <- list(
          Sigma.W = fit@implied$cov[[1]],
          Sigma.B = fit@implied$cov[[2]],
          Mu.W = fit@implied$mean[[1]],
          Mu.B = fit@implied$mean[[2]],
          logl = fit@loglik$loglik
        )
        # if(lavoptions$fixed.x) {
        #   OUT$logl <- OUT$logl - lavsamplestats@YLp[[g]][[2]]$loglik.x
        # }
      } else {
        # complete data
        out_1 <- lav_mvn_cl_em_sat(
          ylp = lavsamplestats@YLp[[g]],
          lp = lavdata@Lp[[g]],
          tol = 1e-04, # option?
          min_variance = 1e-05, # option?
          max_iter = 5000L
        ) # option?
      }
      if (lav_verbose()) {
        cat("\n")
      }

      # if any near-zero within variance(s), produce warning here
      zero_var <- which(diag(out_1$Sigma.W) <= 1e-05)
      if (length(zero_var)) {
        gtxt <- if (ngroups > 1L) {
          gettextf(" in group %s.", g)
        } else {
          " "
        }
        lav_msg_warn(gettextf(
          "H1 estimation resulted in a within covariance matrix %1$s with
          (near) zero variances for some of the level-1 variables: %2$s",
            gtxt, lav_msg_view(lavdata@ov.names.l[[g]][[1]][zero_var]))
        )
      }

      # new in 0.6-18: ensure Mu.W[both.idx] is zero (post-estimation!)
      # (not correctly; fixed in 0.6-19...)
      both_idx   <- lavdata@Lp[[g]]$both.idx[[2]]
      within_idx <- lavdata@Lp[[g]]$within.idx[[2]]
      ov_idx     <- lavdata@Lp[[g]]$ov.idx
      p_tilde <- length(unique(c(ov_idx[[1]], ov_idx[[2]])))
      mu_w_tilde <- mu_b_tilde <- mu_wb_tilde <- numeric(p_tilde)
      mu_w_tilde[ov_idx[[1]]] <- out_1$Mu.W
      mu_b_tilde[ov_idx[[2]]] <- out_1$Mu.B
      mu_wb_tilde[both_idx] <- (mu_b_tilde[both_idx] + mu_w_tilde[both_idx])
      mu_w_tilde[both_idx] <- 0
      mu_b_tilde[both_idx] <- mu_wb_tilde[both_idx]
      out_1$Mu.W <- mu_w_tilde[ov_idx[[1]]]
      out_1$Mu.B <- mu_b_tilde[ov_idx[[2]]]

      # store in implied
      implied$cov[[(g - 1) * nlevels + 1L]] <- out_1$Sigma.W
      implied$cov[[(g - 1) * nlevels + 2L]] <- out_1$Sigma.B
      implied$mean[[(g - 1) * nlevels + 1L]] <- out_1$Mu.W
      implied$mean[[(g - 1) * nlevels + 2L]] <- out_1$Mu.B

      # store logl per group
      loglik_group[g] <- out_1$logl
    }

    logl <- list(loglik = sum(loglik_group), loglik.group = loglik_group)
  }

  list(implied = implied, logl = logl)
}
