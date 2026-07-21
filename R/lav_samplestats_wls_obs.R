lav_samp_wls_obs <- function(mean_g, cov_g, var_g,
                                    th_g, th_idx_g,
                                    res_int_g, res_cov_g, res_var_g, res_th_g,
                                    res_slopes_g,
                                    group_w_g,
                                    categorical = FALSE,
                                    conditional_x = FALSE,
                                    meanstructure = FALSE,
                                    correlation = FALSE,
                                    num_idx = integer(0L),
                                    slopestructure = FALSE,
                                    group_w_free = FALSE) {
  # WLS.obs
  if (categorical) {
    # order of elements is important here:
    # 1. thresholds + (negative) means (interleaved)
    # 2. slopes (if any)
    # 3. variances (if any)
    # 4. covariance matrix (no diagonal!)

    # NOTE: prior to 0.5-17, we had this:
    # TH[ov.types == "numeric"] <- -1*TH[ov.types == "numeric"]
    # which is WRONG if we have more than one threshold per variable
    # (thanks to Sacha Epskamp for spotting this!)
    if (conditional_x) {
      th <- res_th_g
      th[th_idx_g == 0] <- -1 * th[th_idx_g == 0]

      nvar <- length(res_var_g)
      num_idx <- which(!seq_len(nvar) %in% th_idx_g)

      wls_obs <- c(
        th,
        lav_mat_vec(res_slopes_g),
        res_var_g[num_idx],
        lav_mat_vech(res_cov_g, diagonal = FALSE)
      )
    } else {
      th <- th_g
      th[th_idx_g == 0] <- -1 * th[th_idx_g == 0]

      nvar <- length(var_g)
      num_idx <- which(!seq_len(nvar) %in% th_idx_g)

      wls_obs <- c(
        th,
        var_g[num_idx],
        lav_mat_vech(cov_g, diagonal = FALSE)
      )
    }
  } else {
    # CONTINUOUS:
    diag_1 <- TRUE
    if (correlation) {
      diag_1 <- FALSE
    }
    if (conditional_x) {
      if (correlation && length(num_idx) > 0L) {
        # partial correlation structure: the non-standardized y variables
        # keep a residual-variance row (before the off-diagonal elements,
        # mirroring the categorical layout above); num_idx indexes the
        # y (res.cov) block
        v_res <- c(
          diag(res_cov_g)[num_idx],
          lav_mat_vech(res_cov_g, diagonal = FALSE)
        )
      } else {
        v_res <- lav_mat_vech(res_cov_g, diagonal = diag_1)
      }
      if (meanstructure) {
        if (slopestructure) {
          # order = vec(Beta), where first row are intercepts
          # cbind(res.int, res.slopes) is t(Beta)
          # so we need vecr
          wls_obs <- c(
            lav_mat_vecr(cbind(
              res_int_g,
              res_slopes_g
            )),
            v_res
          )
        } else {
          wls_obs <- c(
            res_int_g,
            v_res
          )
        }
      } else {
        if (slopestructure) {
          wls_obs <- c(
            lav_mat_vecr(res_slopes_g),
            v_res
          )
        } else {
          wls_obs <- v_res
        }
      }
    } else {
      if (correlation && length(num_idx) > 0L) {
        # partial correlation structure: keep the variances of the
        # non-correlation variables (num_idx), drop only the correlation
        # variables' diagonal (mirrors the categorical layout above)
        v_obs <- c(var_g[num_idx], lav_mat_vech(cov_g, diagonal = FALSE))
      } else {
        v_obs <- lav_mat_vech(cov_g, diagonal = diag_1)
      }
      if (meanstructure) {
        wls_obs <- c(mean_g, v_obs)
      } else {
        wls_obs <- v_obs
      }
    }
  }

  # group.w.free?
  if (group_w_free) {
    wls_obs <- c(group_w_g, wls_obs)
  }

  wls_obs
}
