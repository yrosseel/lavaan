# compute logl for the unrestricted (h1) model -- per group
lav_h1_logl <- function(lavdata = NULL,
                        lavsamplestats = NULL,
                        lavoptions = NULL) {
  # number of groups
  ngroups <- lavdata@ngroups

  logl.group <- rep(as.numeric(NA), ngroups)

  # should compute logl, or return NA?
  logl.ok <- FALSE
  if (lavoptions$estimator %in% c("ML", "MML")) {
    # check if everything is numeric, OR if we have exogenous
    # factor with 2 levels only
    if (all(lavdata@ov$type == "numeric")) {
      logl.ok <- TRUE
    } else {
      not.idx <- which(lavdata@ov$type != "numeric")
      for (i in not.idx) {
        if (lavdata@ov$type[i] == "factor" &&
          lavdata@ov$exo[i] == 1L &&
          lavdata@ov$nlev[i] == 2L) {
          logl.ok <- TRUE
        } else {
          logl.ok <- FALSE
          break
        }
      }
    }
  }

  # lavsamplestats filled in? (not if no data, or samplestats = FALSE)
  if (length(lavsamplestats@ntotal) == 0L ||
      (!is.null(lavoptions$samplestats) && !lavoptions$samplestats)) {
    logl.ok <- FALSE
  }

  # new in 0.6-9 (so SAM can handle N<P)
  if (!is.null(lavoptions$sample.icov) && !lavoptions$sample.icov) {
    logl.ok <- FALSE
  }

  if (logl.ok) {
    for (g in seq_len(ngroups)) {
      if (lavdata@nlevels > 1L) {
        OUT <- lav_mvnorm_cluster_em_sat(
          YLp = lavsamplestats@YLp[[g]],
          Lp = lavdata@Lp[[g]],
          verbose = FALSE,
          tol = 1e-04, # option?
          min.variance = 1e-05, # option?
          max.iter = 5000L
        ) # option?
        # store logl per group
        logl.group[g] <- OUT$logl
      } else if (lavsamplestats@missing.flag) {
        logl.group[g] <-
          lav_mvnorm_missing_loglik_samplestats(
            Yp     = lavsamplestats@missing[[g]],
            Mu     = lavsamplestats@missing.h1[[g]]$mu,
            Sigma  = lavsamplestats@missing.h1[[g]]$sigma,
            x.idx  = lavsamplestats@x.idx[[g]],
            x.mean = lavsamplestats@mean.x[[g]],
            x.cov  = lavsamplestats@cov.x[[g]]
          )
      } else { # single-level, complete data
        # all we need is: logdet of covariance matrix, nobs and nvar
        if (lavoptions$conditional.x) {
          logl.group[g] <-
            lav_mvnorm_h1_loglik_samplestats(
              sample.cov.logdet =
                lavsamplestats@res.cov.log.det[[g]],
              sample.nvar =
                NCOL(lavsamplestats@res.cov[[g]]),
              sample.nobs = lavsamplestats@nobs[[g]]
            )
        } else {
          logl.group[g] <-
            lav_mvnorm_h1_loglik_samplestats(
              sample.cov.logdet = lavsamplestats@cov.log.det[[g]],
              sample.nvar       = NCOL(lavsamplestats@cov[[g]]),
              sample.nobs       = lavsamplestats@nobs[[g]],
              x.idx             = lavsamplestats@x.idx[[g]],
              x.cov             = lavsamplestats@cov.x[[g]]
            )
        }
      } # complete
    } # g
  } # logl.ok is TRUE

  out <- list(
    loglik = sum(logl.group),
    loglik.group = logl.group
  )

  out
}
