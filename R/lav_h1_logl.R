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
        # EM control arguments (fall back to the multilevel defaults if
        # lavoptions is incomplete)
        em_h1_tol <- lavoptions$em.h1.args$tol
        if (!is.numeric(em_h1_tol)) {
          em_h1_tol <- 1e-04
        }
        em_h1_iter_max <- lavoptions$em.h1.args$iter_max
        if (!is.numeric(em_h1_iter_max)) {
          em_h1_iter_max <- 5000L
        }
        em_h1_min_variance <- lavoptions$em.h1.args$min_variance
        if (!is.numeric(em_h1_min_variance)) {
          em_h1_min_variance <- 1e-05
        }
        em_h1_accel <- lavoptions$em.h1.args$acceleration
        if (is.null(em_h1_accel)) {
          em_h1_accel <- "none" # backwards compatibility
        }
        out_1 <- lav_mvn_cl_em_sat(
          ylp = lavsamplestats@YLp[[g]],
          lp = lavdata@Lp[[g]],
          tol = em_h1_tol,
          min_variance = em_h1_min_variance,
          max_iter = em_h1_iter_max,
          acceleration = em_h1_accel
        )
        lav_verbose(current_verbose)
        # store logl per group
        logl_group[g] <- out_1$logl
      } else if (lavsamplestats@missing.flag) {
        logl_group[g] <-
          lav_mvn_mi_loglik_samp(
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
            lav_mvn_h1_loglik_samp(
              sample_cov_logdet =
                lavsamplestats@res.cov.log.det[[g]],
              sample_nvar =
                NCOL(lavsamplestats@res.cov[[g]]),
              sample_nobs = lavsamplestats@nobs[[g]]
            )
        } else {
          logl_group[g] <-
            lav_mvn_h1_loglik_samp(
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
