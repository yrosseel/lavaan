# compute the loglikelihood of the data, given the current values of the
# model parameters
lav_model_loglik <- function(lavdata = NULL,
                             lavsamplestats = NULL,
                             lavh1 = NULL,
                             lavimplied = NULL,
                             lavmodel = NULL,
                             lavoptions = NULL) {
  ngroups <- lavdata@ngroups

  logl.group <- rep(as.numeric(NA), ngroups)

  # should compute logl, or return NA?
  logl.ok <- FALSE
  if (lavoptions$estimator %in% c("ML", "MML")) {
    # check if everything is numeric, OR if we have exogenous
    # factor with 2 levels only
    # if(all(lavdata@ov$type == "numeric")) {
    logl.ok <- TRUE
    # } else {
    if (lavoptions$fixed.x == FALSE) {
      exo.idx <- which(lavdata@ov$exo == 1L)
      for (i in exo.idx) {
        if (lavdata@ov$nlev[i] > 1L) {
          logl.ok <- FALSE
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
    logl.ok <- FALSE
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
        logl.ok <- FALSE
      }
    } else {
      if (all(lavimplied$cov[[b]] == 0)) {
        logl.ok <- FALSE
      }
    }
  }


  if (logl.ok) {
    for (g in seq_len(ngroups)) {
      if (lavdata@nlevels > 1L) {
        # here, we assume only 2 levels, at [[1]] and [[2]]
        if (lavmodel@conditional.x) {
          Res.Sigma.W <- lavimplied$res.cov[[(g - 1) * 2 + 1]]
          Res.Int.W <- lavimplied$res.int[[(g - 1) * 2 + 1]]
          Res.Pi.W <- lavimplied$res.slopes[[(g - 1) * 2 + 1]]

          Res.Sigma.B <- lavimplied$res.cov[[(g - 1) * 2 + 2]]
          Res.Int.B <- lavimplied$res.int[[(g - 1) * 2 + 2]]
          Res.Pi.B <- lavimplied$res.slopes[[(g - 1) * 2 + 2]]
        } else {
          Sigma.W <- lavimplied$cov[[(g - 1) * 2 + 1]]
          Mu.W <- lavimplied$mean[[(g - 1) * 2 + 1]]
          Sigma.B <- lavimplied$cov[[(g - 1) * 2 + 2]]
          Mu.B <- lavimplied$mean[[(g - 1) * 2 + 2]]
        }

        if (lavsamplestats@missing.flag) {
          if (lavmodel@conditional.x) {
            # TODO
            logl.group[g] <- as.numeric(NA)
          } else {
            logl.group[g] <-
              lav_mvnorm_cluster_missing_loglik_samplestats_2l(
                y1 = lavdata@X[[g]],
                y2 = lavsamplestats@YLp[[g]][[2]]$Y2,
                lp = lavdata@Lp[[g]],
                mp = lavdata@Mp[[g]],
                mu_w = Mu.W, sigma_w = Sigma.W,
                mu_b = Mu.B, sigma_b = Sigma.B,
                loglik_x = lavsamplestats@YLp[[g]][[2]]$loglik.x,
                log2pi = TRUE, minus_two = FALSE
              )
          }
        } else {
          # complete case
          if (lavmodel@conditional.x) {
            logl.group[g] <-
              lav_mvreg_cluster_loglik_samplestats_2l(
                YLp          = lavsamplestats@YLp[[g]],
                Lp           = lavdata@Lp[[g]],
                Res.Sigma.W  = Res.Sigma.W,
                Res.Int.W    = Res.Int.W,
                Res.Pi.W     = Res.Pi.W,
                Res.Sigma.B  = Res.Sigma.B,
                Res.Int.B    = Res.Int.B,
                Res.Pi.B     = Res.Pi.B,
                Sinv.method  = "eigen",
                log2pi       = TRUE,
                minus.two    = FALSE
              )
          } else {
            logl.group[g] <-
              lav_mvnorm_cluster_loglik_samplestats_2l(
                ylp          = lavsamplestats@YLp[[g]],
                lp           = lavdata@Lp[[g]],
                mu_w         = Mu.W,
                sigma_w      = Sigma.W,
                mu_b         = Mu.B,
                sigma_b      = Sigma.B,
                sinv_method  = "eigen",
                log2pi       = TRUE,
                minus_two    = FALSE
              )
          }
        } # complete
        # end multilevel
      } else if (lavsamplestats@missing.flag) {
        x.idx <- lavsamplestats@x.idx[[g]]
        X.MEAN <- X.COV <- NULL
        if (length(x.idx) > 0L) {
          X.MEAN <- lavh1$implied$mean[[g]][x.idx]
          X.COV <- lavh1$implied$cov[[g]][x.idx, x.idx, drop = FALSE]
        }
        logl.group[g] <- lav_mvnorm_missing_loglik_samplestats(
          yp     = lavsamplestats@missing[[g]],
          mu     = lavimplied$mean[[g]],
          sigma_1  = lavimplied$cov[[g]],
          x_idx  = lavsamplestats@x.idx[[g]],
          x_mean = X.MEAN, # not needed? should be part of Sigma
          x_cov  = X.COV
        ) # not needed at all!
      } else { # single-level, complete data
        if (lavoptions$conditional.x) {
          # FIXME: use lavh1
          logl.group[g] <- lav_mvreg_loglik_samplestats(
            sample.res.int    = lavsamplestats@res.int[[g]],
            sample.res.slopes = lavsamplestats@res.slopes[[g]],
            sample.res.cov    = lavsamplestats@res.cov[[g]],
            sample.mean.x     = lavsamplestats@mean.x[[g]],
            sample.cov.x      = lavsamplestats@cov.x[[g]],
            sample.nobs       = lavsamplestats@nobs[[g]],
            res.int           = lavimplied$res.int[[g]],
            res.slopes        = lavimplied$res.slopes[[g]],
            res.cov           = lavimplied$res.cov[[g]],
            Sinv.method       = "eigen"
          )
        } else {
          if (lavoptions$meanstructure) {
            Mu <- lavimplied$mean[[g]]
          } else {
            Mu <- lavsamplestats@mean[[g]]
          }
          logl.group[g] <- lav_mvnorm_loglik_samplestats(
            sample_mean = lavsamplestats@mean[[g]],
            sample_cov  = lavsamplestats@cov[[g]],
            sample_nobs = lavsamplestats@nobs[[g]],
            mu          = Mu,
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
  logl <- sum(logl.group)

  # number of parameters, taking into account any equality constraints
  npar <- lavmodel@nx.free
  ceq.simple.only <- lavmodel@ceq.simple.only
  cin.simple.only <- lavmodel@cin.simple.only

  if (ceq.simple.only) {
    # nothing to do
  } else if (!cin.simple.only && nrow(lavmodel@con.jac) > 0L) {
    ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
    if (length(ceq.idx) > 0L) {
      neq <- qr(lavmodel@con.jac[ceq.idx, , drop = FALSE])$rank
      npar <- npar - neq
    }
  }

  # logl
  logl <- sum(logl.group)

  if (logl.ok) {
    # AIC
    AIC <- lav_fit_aic(logl = logl, npar = npar)

    # BIC
    BIC <- lav_fit_bic(logl = logl, npar = npar, n = lavsamplestats@ntotal)

    # BIC2
    BIC2 <- lav_fit_sabic(
      logl = logl, npar = npar,
      n = lavsamplestats@ntotal
    )
  } else {
    AIC <- BIC <- BIC2 <- as.numeric(NA)
  }

  out <- list(
    loglik = logl,
    loglik.group = logl.group,
    npar = npar,
    ntotal = lavsamplestats@ntotal,
    AIC = AIC,
    BIC = BIC,
    BIC2 = BIC2,
    estimator = lavoptions$estimator,
    conditional.x = lavoptions$conditional.x,
    fixed.x = lavoptions$fixed.x
  )

  out
}
