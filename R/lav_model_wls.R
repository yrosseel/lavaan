# compute WLS.est (as a list per group)
#
# multilevel: the per-block (within/between) vectors are stacked per group:
# c(within, between), matching WLS.obs (lav_samp_wls_2l) and the per-group
# Delta matrix (lav_model_delta)
lav_model_wls_est <- function(lavmodel = NULL, glist = NULL,
                              lavimplied = NULL) {
  # two-level + categorical: the standardized stage-wise layout
  # (see lav_model_wls_2l_cat.R); already stacked per group
  if (lavmodel@multilevel && lavmodel@categorical) {
    return(lav_m07_wls_est(lavmodel = lavmodel, glist = glist))
  }
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
        if (correlation && length(num_idx[[g]]) > 0L) {
          # partial correlation structure: keep the variances of the
          # non-correlation variables (num.idx), drop only the correlation
          # variables' diagonal (mirrors the categorical layout above)
          v_est <- c(
            diag(lavimplied$cov[[g]])[num_idx[[g]]],
            lav_mat_vech(lavimplied$cov[[g]], diagonal = FALSE)
          )
        } else {
          v_est <- lav_mat_vech(lavimplied$cov[[g]], diagonal = diag_1)
        }
        if (meanstructure) {
          wls_est <- c(lavimplied$mean[[g]], v_est)
        } else {
          wls_est <- v_est
        }
      } # conditional.x = FALSE
    } # categorical = FALSE

    if (group_w_free) {
      wls_est <- c(lavimplied$group.w[[g]], wls_est)
    }

    wls_est_1[[g]] <- wls_est
  }

  # multilevel: stack the blocks of each group (within on top)
  if (lavmodel@multilevel) {
    nlevels <- 2L
    ngroups <- nblocks %/% nlevels
    wls_est_2 <- vector("list", length = ngroups)
    for (g in seq_len(ngroups)) {
      wls_est_2[[g]] <- c(
        wls_est_1[[(g - 1) * nlevels + 1L]],
        wls_est_1[[(g - 1) * nlevels + 2L]]
      )
    }
    wls_est_1 <- wls_est_2
  }

  wls_est_1
}

# Note: lav_model_wls_v() is replaced by lav_model_h1_info() in 0.6-1
