# compute the loglikelihood of the data, given the current values of the
# model parameters
lav_model_loglik <- function(lavdata = NULL,
                             lavsamplestats = NULL,
                             lavh1 = NULL,
                             lavimplied = NULL,
                             lavmodel = NULL,
                             lavoptions = NULL) {
  ngroups <- lavdata@ngroups

  logl_group <- rep(as.numeric(NA), ngroups)

  # should compute logl, or return NA?
  logl_ok <- FALSE
  if (lavoptions$estimator %in% c("ML", "MML")) {
    # check if everything is numeric, OR if we have exogenous
    # factor with 2 levels only
    # if(all(lavdata@ov$type == "numeric")) {
    logl_ok <- TRUE
    # } else {
    if (lavoptions$fixed.x == FALSE) {
      exo_idx <- which(lavdata@ov$exo == 1L)
      for (i in exo_idx) {
        if (lavdata@ov$nlev[i] > 1L) {
          logl_ok <- FALSE
        }
      }
    }
    # nlevels + fiml
    # if(lavdata@nlevels > 1L && lavsamplestats@missing.flag) {
    #    logl.ok <- FALSE
    # }
  }

  # lavsamplestats filled in? (not if no data, or samplestats = FALSE)
  if (length(lavsamplestats@ntotal) == 0L ||
      (!is.null(lavoptions$samplestats) && !lavoptions$samplestats)) {
    logl_ok <- FALSE
  }

  # catch all-zero Sigma (new in 0.6-20)
  nblocks <- lavmodel@nblocks
  for (b in seq_len(nblocks)) {
    # except for level-2, where Sigma could be all zero (0.6-22)
    if (lavdata@nlevels > 1L && ((b %% lavdata@nlevels) != 1)) {
      next
    }
    if (lavmodel@conditional.x) {
      if (all(lavimplied$res.cov[[b]] == 0)) {
        logl_ok <- FALSE
      }
    } else {
      if (all(lavimplied$cov[[b]] == 0)) {
        logl_ok <- FALSE
      }
    }
  }


  if (logl_ok) {
    for (g in seq_len(ngroups)) {
      if (lavdata@nlevels > 1L) {
        # here, we assume only 2 levels, at [[1]] and [[2]]
        if (lavmodel@conditional.x) {
          res_sigma_w <- lavimplied$res.cov[[(g - 1) * 2 + 1]]
          res_int_w <- lavimplied$res.int[[(g - 1) * 2 + 1]]
          res_pi_w <- lavimplied$res.slopes[[(g - 1) * 2 + 1]]

          res_sigma_b <- lavimplied$res.cov[[(g - 1) * 2 + 2]]
          res_int_b <- lavimplied$res.int[[(g - 1) * 2 + 2]]
          res_pi_b <- lavimplied$res.slopes[[(g - 1) * 2 + 2]]
        } else {
          sigma_w <- lavimplied$cov[[(g - 1) * 2 + 1]]
          mu_w <- lavimplied$mean[[(g - 1) * 2 + 1]]
          sigma_b <- lavimplied$cov[[(g - 1) * 2 + 2]]
          mu_b <- lavimplied$mean[[(g - 1) * 2 + 2]]
        }

        if (lavsamplestats@missing.flag) {
          if (lavmodel@conditional.x) {
            # TODO
            logl_group[g] <- as.numeric(NA)
          } else {
            logl_group[g] <-
              lav_mvnorm_cluster_missing_loglik_samplestats_2l(
                y1 = lavdata@X[[g]],
                y2 = lavsamplestats@YLp[[g]][[2]]$Y2,
                lp = lavdata@Lp[[g]],
                mp = lavdata@Mp[[g]],
                mu_w = mu_w, sigma_w = sigma_w,
                mu_b = mu_b, sigma_b = sigma_b,
                loglik_x = lavsamplestats@YLp[[g]][[2]]$loglik.x,
                log2pi = TRUE, minus_two = FALSE
              )
          }
        } else {
          # complete case
          if (lavmodel@conditional.x) {
            logl_group[g] <-
              lav_mvreg_cluster_loglik_samplestats_2l(
                ylp          = lavsamplestats@YLp[[g]],
                lp           = lavdata@Lp[[g]],
                res_sigma_w  = res_sigma_w,
                res_int_w    = res_int_w,
                res_pi_w     = res_pi_w,
                res_sigma_b  = res_sigma_b,
                res_int_b    = res_int_b,
                res_pi_b     = res_pi_b,
                sinv_method  = "eigen",
                log2pi       = TRUE,
                minus_two    = FALSE
              )
          } else {
            logl_group[g] <-
              lav_mvnorm_cluster_loglik_samplestats_2l(
                ylp          = lavsamplestats@YLp[[g]],
                lp           = lavdata@Lp[[g]],
                mu_w         = mu_w,
                sigma_w      = sigma_w,
                mu_b         = mu_b,
                sigma_b      = sigma_b,
                sinv_method  = "eigen",
                log2pi       = TRUE,
                minus_two    = FALSE
              )
          }
        } # complete
        # end multilevel
      } else if (lavsamplestats@missing.flag) {
        x_idx <- lavsamplestats@x.idx[[g]]
        x_mean <- x_cov <- NULL
        if (length(x_idx) > 0L) {
          x_mean <- lavh1$implied$mean[[g]][x_idx]
          x_cov <- lavh1$implied$cov[[g]][x_idx, x_idx, drop = FALSE]
        }
        logl_group[g] <- lav_mvnorm_missing_loglik_samplestats(
          yp     = lavsamplestats@missing[[g]],
          mu     = lavimplied$mean[[g]],
          sigma_1  = lavimplied$cov[[g]],
          x_idx  = lavsamplestats@x.idx[[g]],
          x_mean = x_mean, # not needed? should be part of Sigma
          x_cov  = x_cov
        ) # not needed at all!
      } else { # single-level, complete data
        if (lavoptions$conditional.x) {
          # FIXME: use lavh1
          logl_group[g] <- lav_mvreg_loglik_samplestats(
            sample_res_int    = lavsamplestats@res.int[[g]],
            sample_res_slopes = lavsamplestats@res.slopes[[g]],
            sample_res_cov    = lavsamplestats@res.cov[[g]],
            sample_mean_x     = lavsamplestats@mean.x[[g]],
            sample_cov_x      = lavsamplestats@cov.x[[g]],
            sample_nobs       = lavsamplestats@nobs[[g]],
            res_int           = lavimplied$res.int[[g]],
            res_slopes        = lavimplied$res.slopes[[g]],
            res_cov           = lavimplied$res.cov[[g]],
            sinv_method       = "eigen"
          )
        } else {
          if (lavoptions$meanstructure) {
            mu <- lavimplied$mean[[g]]
          } else {
            mu <- lavsamplestats@mean[[g]]
          }
          logl_group[g] <- lav_mvnorm_loglik_samplestats(
            sample_mean = lavsamplestats@mean[[g]],
            sample_cov  = lavsamplestats@cov[[g]],
            sample_nobs = lavsamplestats@nobs[[g]],
            mu          = mu,
            sigma_1     = lavimplied$cov[[g]],
            x_idx       = lavsamplestats@x.idx[[g]],
            x_mean      = lavsamplestats@mean.x[[g]],
            x_cov       = lavsamplestats@cov.x[[g]],
            sinv_method = "eigen",
            sigma_inv   = NULL
          )
        }
      } # complete
    } # g
  } # logl.ok is TRUE

  # logl
  logl <- sum(logl_group)

  # number of parameters, taking into account any equality constraints
  npar <- lavmodel@nx.free
  ceq_simple_only <- lavmodel@ceq.simple.only
  cin_simple_only <- lavmodel@cin.simple.only

  if (ceq_simple_only) {
    # nothing to do
  } else if (!cin_simple_only && nrow(lavmodel@con.jac) > 0L) {
    ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
    if (length(ceq_idx) > 0L) {
      neq <- qr(lavmodel@con.jac[ceq_idx, , drop = FALSE])$rank
      npar <- npar - neq
    }
  }

  # logl
  logl <- sum(logl_group)

  if (logl_ok) {
    # AIC
    aic <- lav_fit_aic(logl = logl, npar = npar)

    # BIC
    bic <- lav_fit_bic(logl = logl, npar = npar, n = lavsamplestats@ntotal)

    # BIC2
    bic2 <- lav_fit_sabic(
      logl = logl, npar = npar,
      n = lavsamplestats@ntotal
    )
  } else {
    aic <- bic <- bic2 <- as.numeric(NA)
  }

  out <- list(
    loglik = logl,
    loglik.group = logl_group,
    npar = npar,
    ntotal = lavsamplestats@ntotal,
    AIC = aic,
    BIC = bic,
    BIC2 = bic2,
    estimator = lavoptions$estimator,
    conditional.x = lavoptions$conditional.x,
    fixed.x = lavoptions$fixed.x
  )

  out
}
