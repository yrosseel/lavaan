# create random starting values starting from a parameter table
# - using the lower/upper bounds and runif() for factor loadings
#   and variances
# - using runif(,-0.5,+0.5) for correlations; rescale using variances
# - check if Sigma.hat is PD; if not, try again
#
# YR 26 Feb 2024

lav_partable_random <- function(lavpartable = NULL,
                                # needed if we still need to compute bounds:
                                lavh1 = NULL, lavdata = NULL,
                                lavsamplestats = NULL, lavoptions = NULL) {

  lavpta <- lav_partable_attributes(lavpartable)

  # ALWAYS (recompute) bounds, as user may have provide other
  # bounds (eg "pos.var") (0.6-20)
  lavoptions2 <- lavoptions
  lavoptions2$bounds <- "standard"
  lavoptions2$optim.bounds <-
    list(
      lower = c(
        "ov.var", "lv.var", "loadings",
        "covariances"
      ),
      upper = c(
        "ov.var", "lv.var", "loadings",
        "covariances"
      ),
      lower.factor = c(1.0, 1.0, 1.0, 0.999),
      upper.factor = c(1.0, 1.0, 1.0, 0.999),
      min.reliability.marker = 0.1,
      min.var.lv.endo = 0.005
    )
  lavpartable <- lav_partable_add_bounds(
    partable = lavpartable,
    lavh1 = lavh1, lavdata = lavdata,
    lavsamplestats = lavsamplestats, lavoptions = lavoptions2
  )

  # replace -Inf/Inf by -1/1 * .Machine$double.eps (for runif)
  inf_idx <- which(lavpartable$lower < -1e+16)
  if (length(inf_idx) > 0L) {
    lavpartable$lower[inf_idx] <- -1e+16
  }
  inf_idx <- which(lavpartable$upper > 1e+16)
  if (length(inf_idx) > 0L) {
    lavpartable$upper[inf_idx] <- 1e+16
  }

  # empty lavpartable$start?
  if (is.null(lavpartable$start)) {
    start_1 <- numeric(length(lavpartable$lhs))
    # set loadings to 0.7
    loadings_idx <- which(lavpartable$free > 0L &
      lavpartable$op == "=~")
    if (length(loadings_idx) > 0L) {
      start_1[loadings_idx] <- 0.7
    }
    # set (only) variances to 1
    var_idx <- which(lavpartable$free > 0L &
      lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs)
    if (length(var_idx) > 0L) {
      start_1[var_idx] <- 1
    }

    lavpartable$start <- start_1
  }

  # initial values
  start_1 <- lavpartable$start

  nblocks <- lav_partable_nblocks(lavpartable)
  block_values <- lav_partable_block_values(lavpartable)
  for (b in 1:nblocks) {
    ov_names <- lavpta$vnames$ov[[b]]
    lv_names <- lavpta$vnames$lv[[b]]
    ov_ind_names <- lavpta$vnames$ov.ind[[b]]

    # start with the lv (residual) variances
    lv_var_idx <- which(lavpartable$block == block_values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% lv_names &
      lavpartable$rhs %in% lv_names &
      lavpartable$lhs == lavpartable$rhs)

    if (length(lv_var_idx) > 0L) {
      for (i in lv_var_idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          start_1[i] <- runif(
            n = 1L, min = lavpartable$lower[i],
            max = lavpartable$upper[i]
          )
        }
      }
    }

    # first, we generate lv correlations, and then rescale to covariances
    lv_cov_idx <- which(lavpartable$block == block_values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% lv_names &
      lavpartable$rhs %in% lv_names &
      lavpartable$lhs != lavpartable$rhs)

    if (length(lv_cov_idx) > 0L) {
      for (i in lv_cov_idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          cor_val <- runif(n = 1L, -0.5, +0.5)
          var1_idx <- which(lavpartable$block == block_values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$lhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          var2_idx <- which(lavpartable$block == block_values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          start_1[i] <- cor_val * sqrt(start_1[var1_idx]) *
                                  sqrt(start_1[var2_idx])
        }
      }
    }

    # next, (residual) ov variances
    ov_var_idx <- which(lavpartable$block == block_values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% ov_names &
      lavpartable$rhs %in% ov_names &
      lavpartable$lhs == lavpartable$rhs)

    if (length(ov_var_idx) > 0L) {
      for (i in ov_var_idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          start_1[i] <- runif(
            n = 1L, min = lavpartable$lower[i],
            max = lavpartable$upper[i]
          )
        }
      }
    }

    # (residual) ov covariances (if any)
    ov_cov_idx <- which(lavpartable$block == block_values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% ov_names &
      lavpartable$rhs %in% ov_names &
      lavpartable$lhs != lavpartable$rhs)

    if (length(ov_cov_idx) > 0L) {
      for (i in ov_cov_idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          cor_val <- runif(n = 1L, -0.5, +0.5)
          var1_idx <- which(lavpartable$block == block_values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$lhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          var2_idx <- which(lavpartable$block == block_values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          start_1[i] <- cor_val * sqrt(start_1[var1_idx]) *
                                  sqrt(start_1[var2_idx])
        }
      }
    }

    # finally, the lambda values, keeping in mind that
    # lambda_p^(u) = sqrt( upper(res.var.indicators_p) /
    #                      lower(var.factor) )
    lambda_idx <- which(lavpartable$block == block_values[b] &
      lavpartable$op == "=~" &
      lavpartable$lhs %in% lv_names &
      lavpartable$rhs %in% ov_ind_names)

    if (length(lambda_idx)) {
      for (i in lambda_idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          varov_idx <- which(lavpartable$block == block_values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          varlv_idx <- which(lavpartable$block == block_values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$lhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          lambda_u <- sqrt(start_1[varov_idx] / start_1[varlv_idx])
          start_1[i] <- runif(n = 1, -lambda_u, lambda_u)
        }
      }
    }
  }

  # sanity check; needed?
  current_warn <- lav_warn()
  if (lav_warn(TRUE))
      on.exit(lav_warn(current_warn), TRUE)
  start_1 <- lav_start_check_cov(
    lavpartable = lavpartable, start = start_1
  )

  start_1
}
