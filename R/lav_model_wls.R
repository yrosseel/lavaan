# compute WLS.est (as a list per group)
lav_model_wls_est <- function(lavmodel = NULL, glist = NULL,
                              lavimplied = NULL) {
  nblocks <- lavmodel@nblocks
  meanstructure <- lavmodel@meanstructure
  correlation <- lavmodel@correlation
  categorical <- lavmodel@categorical
  group_w_free <- lavmodel@group.w.free
  num_idx <- lavmodel@num.idx

  # model-implied statistics
  if (is.null(lavimplied)) {
    lavimplied <- lav_model_implied(lavmodel, glist = glist)
  }

  wls_est_1 <- vector("list", length = nblocks)
  for (g in 1:nblocks) {
    if (categorical) {
      # order of elements is important here:
      # 1. thresholds + means (interleaved)
      # 2. slopes (if any, columnwise per exo)
      # 3. variances (if any)
      # 4. correlations (no diagonal!)
      if (lavmodel@conditional.x) {
        wls_est <- c(
          lavimplied$res.th[[g]],
          lav_mat_vec(lavimplied$res.slopes[[g]]),
          diag(lavimplied$res.cov[[g]])[num_idx[[g]]],
          lav_mat_vech(lavimplied$res.cov[[g]],
            diagonal = FALSE
          )
        )
      } else {
        wls_est <- c(
          lavimplied$th[[g]],
          diag(lavimplied$cov[[g]])[num_idx[[g]]],
          lav_mat_vech(lavimplied$cov[[g]],
            diagonal = FALSE
          )
        )
      }
    } else {
      # CONTINUOUS
      diag_1 <- TRUE
      if (correlation) {
        diag_1 <- FALSE
      }

      if (lavmodel@conditional.x && lavmodel@nexo[g] > 0L) {
        # order = vec(Beta), where first row are intercepts
        # cbind(res.int, res.slopes) is t(Beta)
        # so we need vecr
        if (meanstructure) {
          wls_est <- c(
            lav_mat_vecr(
              cbind(
                lavimplied$res.int[[g]],
                lavimplied$res.slopes[[g]]
              )
            ),
            lav_mat_vech(lavimplied$res.cov[[g]],
              diagonal = diag_1
            )
          )
        } else {
          wls_est <- c(
            lav_mat_vecr(lavimplied$res.slopes[[g]]),
            lav_mat_vech(lavimplied$res.cov[[g]],
              diagonal = diag_1
            )
          )
        }
      } else {
        if (meanstructure) {
          wls_est <- c(
            lavimplied$mean[[g]],
            lav_mat_vech(lavimplied$cov[[g]],
              diagonal = diag_1
            )
          )
        } else {
          wls_est <- lav_mat_vech(lavimplied$cov[[g]],
            diagonal = diag_1
          )
        }
      } # conditional.x = FALSE
    } # categorical = FALSE

    if (group_w_free) {
      wls_est <- c(lavimplied$group.w[[g]], wls_est)
    }

    wls_est_1[[g]] <- wls_est
  }

  wls_est_1
}

# Note: lav_model_wls_v() is replaced by lav_model_h1_info() in 0.6-1
