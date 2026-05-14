# compute logl for the unrestricted (h1) model -- per group
lav_h1_logl <- function(lavdata = NULL,
                        lavsamplestats = NULL,
                        h1_implied = NULL,
                        lavoptions = NULL) {
  # number of groups
  ngroups <- lavdata@ngroups

  logl_group <- rep(as.numeric(NA), ngroups)

  # should compute logl, or return NA?
  logl_ok <- FALSE
  if (lavoptions$estimator %in% c("ML", "MML")) {
    # check if everything is numeric, OR if we have exogenous
    # factor with 2 levels only
    if (all(lavdata@ov$type == "numeric")) {
      logl_ok <- TRUE
    } else {
      not_idx <- which(lavdata@ov$type != "numeric")
      for (i in not_idx) {
        if (lavdata@ov$type[i] == "factor" &&
          lavdata@ov$exo[i] == 1L &&
          lavdata@ov$nlev[i] == 2L) {
          logl_ok <- TRUE
        } else {
          logl_ok <- FALSE
          break
        }
      }
    }
  }

  # lavsamplestats filled in? (not if no data, or samplestats = FALSE)
  if (length(lavsamplestats@ntotal) == 0L ||
      (!is.null(lavoptions$samplestats) && !lavoptions$samplestats)) {
    logl_ok <- FALSE
  }

  # new in 0.6-9 (so SAM can handle N<P)
  if (!is.null(lavoptions$sample.icov) && !lavoptions$sample.icov) {
    logl_ok <- FALSE
  }

  if (logl_ok) {
    for (g in seq_len(ngroups)) {
      if (lavdata@nlevels > 1L) {
        current_verbose <- lav_verbose()
        if (lav_verbose(FALSE))
            on.exit(lav_verbose(current_verbose), TRUE)
        out_1 <- lav_mvnorm_cluster_em_sat(
          ylp = lavsamplestats@YLp[[g]],
          lp = lavdata@Lp[[g]],
          tol = 1e-04, # option?
          min_variance = 1e-05, # option?
          max_iter = 5000L
        ) # option?
        lav_verbose(current_verbose)
        # store logl per group
        logl_group[g] <- out_1$logl
      } else if (lavsamplestats@missing.flag) {
        logl_group[g] <-
          lav_mvnorm_missing_loglik_samplestats(
            yp     = lavsamplestats@missing[[g]],
            #Mu     = lavsamplestats@missing.h1[[g]]$mu,
            mu     = h1_implied$mean[[g]],
            #Sigma  = lavsamplestats@missing.h1[[g]]$sigma,
            sigma_1  = h1_implied$cov[[g]],
            x_idx  = lavsamplestats@x.idx[[g]],
            x_mean = lavsamplestats@mean.x[[g]],
            x_cov  = lavsamplestats@cov.x[[g]]
          )
      } else { # single-level, complete data
        # all we need is: logdet of covariance matrix, nobs and nvar
        if (lavoptions$conditional.x) {
          logl_group[g] <-
            lav_mvnorm_h1_loglik_samplestats(
              sample_cov_logdet =
                lavsamplestats@res.cov.log.det[[g]],
              sample_nvar =
                NCOL(lavsamplestats@res.cov[[g]]),
              sample_nobs = lavsamplestats@nobs[[g]]
            )
        } else {
          logl_group[g] <-
            lav_mvnorm_h1_loglik_samplestats(
              sample_cov_logdet = lavsamplestats@cov.log.det[[g]],
              sample_nvar       = NCOL(lavsamplestats@cov[[g]]),
              sample_nobs       = lavsamplestats@nobs[[g]],
              x_idx             = lavsamplestats@x.idx[[g]],
              x_cov             = lavsamplestats@cov.x[[g]]
            )
        }
      } # complete
    } # g
  } # logl.ok is TRUE

  out <- list(
    loglik = sum(logl_group),
    loglik.group = logl_group
  )

  out
}
