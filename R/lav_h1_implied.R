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

    # EM control arguments (fall back to the multilevel defaults if
    # lavoptions is incomplete, e.g. when refitting objects created by
    # older versions)
    em_h1_tol <- lavoptions$em.h1.args$tol
    if (!is.numeric(em_h1_tol)) {
      em_h1_tol <- 1e-04
    }
    em_h1_max_iter <- lavoptions$em.h1.args$max_iter
    if (!is.numeric(em_h1_max_iter)) {
      em_h1_max_iter <- 5000L
    }
    em_h1_min_variance <- lavoptions$em.h1.args$min_variance
    if (!is.numeric(em_h1_min_variance)) {
      em_h1_min_variance <- 1e-05
    }
    non_pd_action <- lavoptions$em.h1.args$non_pd_action
    if (is.null(non_pd_action)) {
      non_pd_action <- "warn" # backwards compatibility
    }
    non_pd_tol <- lavoptions$em.h1.args$non_pd_tol
    if (!is.numeric(non_pd_tol)) {
      non_pd_tol <- 1e-05
    }
    em_h1_accel <- lavoptions$em.h1.args$acceleration
    if (is.null(em_h1_accel)) {
      em_h1_accel <- "none" # backwards compatibility
    }
    em_h1_fused <- lavoptions$em.h1.args$fused
    if (is.null(em_h1_fused)) {
      em_h1_fused <- TRUE
    }
    implied <- list(
      cov = vector("list", length = ngroups * nlevels),
      mean = vector("list", length = ngroups * nlevels)
    )
    loglik_group <- numeric(lavdata@ngroups)

    for (g in 1:lavdata@ngroups) {
      if (lav_verbose()) {
        cat("\n\nfitting unrestricted (H1) model in group ", g, "\n")
      }
      if (lavsamplestats@missing.flag &&
          (is.null(lavoptions$h1.missing.method) || # backwards compatibility
           lavoptions$h1.missing.method == "em")) {
        # missing data: EM algorithm (new in 0.7-1; this mimics Mplus)
        out_1 <- lav_mvn_cl_mi_em_sat(
          y1 = lavdata@X[[g]],
          y2 = lavsamplestats@YLp[[g]][[2]]$Y2,
          lp = lavdata@Lp[[g]],
          mp = lavdata@Mp[[g]],
          loglik_x = lavsamplestats@YLp[[g]][[2]]$loglik.x,
          tol = em_h1_tol,
          max_iter = em_h1_max_iter,
          min_variance = em_h1_min_variance,
          acceleration = em_h1_accel,
          fused = em_h1_fused
        )
      } else if (lavsamplestats@missing.flag) {
        # missing data: h1.missing.method = "fiml"
        # fit the unrestricted model using a cholesky parameterization
        # and nlminb() (this was the default in 0.6-* and 0.7-1)

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
        lavoptions2$fit.by.level <- FALSE
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
          tol = em_h1_tol,
          min_variance = em_h1_min_variance,
          max_iter = em_h1_max_iter,
          acceleration = em_h1_accel
        )
      }
      if (lav_verbose()) {
        cat("\n")
      }

      # check for (near) singular within/between covariance matrices;
      # this typically happens when a (within or between) variance
      # approaches zero, but the full matrix may also be (near) singular
      if (non_pd_action != "none") {
        gtxt <- if (ngroups > 1L) {
          gettextf(" in group %s", g)
        } else {
          ""
        }
        for (l in 1:2) {
          if (l == 1L) {
            sigma_l <- out_1$Sigma.W
            level_txt <- "within"
          } else {
            sigma_l <- out_1$Sigma.B
            level_txt <- "between"
          }
          if (NROW(sigma_l) == 0L) {
            next
          }
          msg <- NULL
          near_zero_idx <- which(diag(sigma_l) <= non_pd_tol)
          if (length(near_zero_idx) > 0L) {
            msg <- gettextf(
              "H1 estimation resulted in a %1$s covariance matrix%2$s with
               (near) zero variances for some of the observed variables:
               %3$s.", level_txt, gtxt,
              lav_msg_view(lavdata@ov.names.l[[g]][[l]][near_zero_idx]))
          } else {
            ev <- eigen(sigma_l, symmetric = TRUE, only.values = TRUE)$values
            if (any(ev < non_pd_tol)) {
              msg <- gettextf(
                "H1 estimation resulted in a (near) singular %1$s covariance
                 matrix%2$s: its smallest eigenvalue is smaller than %3$s.",
                level_txt, gtxt, format(non_pd_tol))
            }
          }
          if (!is.null(msg)) {
            if (non_pd_action == "stop") {
              lav_msg_stop(paste(msg, gettext(
                "Set the non_pd_action element of the em.h1.args= argument
                 to \"warn\" or \"none\" to continue anyway.")))
            } else {
              lav_msg_warn(msg)
            }
          }
        }
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
