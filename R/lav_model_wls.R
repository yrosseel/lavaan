# compute WLS.est (as a list per group)
lav_model_wls_est <- function(lavmodel = NULL, GLIST = NULL,
                              lavimplied = NULL) {
  nblocks <- lavmodel@nblocks
  meanstructure <- lavmodel@meanstructure
  if (.hasSlot(lavmodel, "correlation")) {
    correlation <- lavmodel@correlation
  } else {
    correlation <- FALSE
  }
  categorical <- lavmodel@categorical
  group.w.free <- lavmodel@group.w.free
  num.idx <- lavmodel@num.idx

  # model-implied statistics
  if (is.null(lavimplied)) {
    lavimplied <- lav_model_implied(lavmodel, GLIST = GLIST)
  }

  WLS.est <- vector("list", length = nblocks)
  for (g in 1:nblocks) {
    if (categorical) {
      # order of elements is important here:
      # 1. thresholds + means (interleaved)
      # 2. slopes (if any, columnwise per exo)
      # 3. variances (if any)
      # 4. correlations (no diagonal!)
      if (lavmodel@conditional.x) {
        wls.est <- c(
          lavimplied$res.th[[g]],
          lav_matrix_vec(lavimplied$res.slopes[[g]]),
          diag(lavimplied$res.cov[[g]])[num.idx[[g]]],
          lav_matrix_vech(lavimplied$res.cov[[g]],
            diagonal = FALSE
          )
        )
      } else {
        wls.est <- c(
          lavimplied$th[[g]],
          diag(lavimplied$cov[[g]])[num.idx[[g]]],
          lav_matrix_vech(lavimplied$cov[[g]],
            diagonal = FALSE
          )
        )
      }
    } else {
      # CONTINUOUS
      DIAG <- TRUE
      if (correlation) {
        DIAG <- FALSE
      }

      if (lavmodel@conditional.x && lavmodel@nexo[g] > 0L) {
        # order = vec(Beta), where first row are intercepts
        # cbind(res.int, res.slopes) is t(Beta)
        # so we need vecr
        if (meanstructure) {
          wls.est <- c(
            lav_matrix_vecr(
              cbind(
                lavimplied$res.int[[g]],
                lavimplied$res.slopes[[g]]
              )
            ),
            lav_matrix_vech(lavimplied$res.cov[[g]],
              diagonal = DIAG
            )
          )
        } else {
          wls.est <- c(
            lav_matrix_vecr(lavimplied$res.slopes[[g]]),
            lav_matrix_vech(lavimplied$res.cov[[g]],
              diagonal = DIAG
            )
          )
        }
      } else {
        if (meanstructure) {
          wls.est <- c(
            lavimplied$mean[[g]],
            lav_matrix_vech(lavimplied$cov[[g]],
              diagonal = DIAG
            )
          )
        } else {
          wls.est <- lav_matrix_vech(lavimplied$cov[[g]],
            diagonal = DIAG
          )
        }
      } # conditional.x = FALSE
    } # categorical = FALSE

    if (group.w.free) {
      wls.est <- c(lavimplied$group.w[[g]], wls.est)
    }

    WLS.est[[g]] <- wls.est
  }

  WLS.est
}

# Note: lav_model_wls_v() is replaced by lav_model_h1_information() in 0.6-1
