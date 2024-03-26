# create random starting values starting from a parameter table
# - using the lower/upper bounds and runif() for factor loadings
#   and variances
# - using runif(,-1,+1) for correlations; rescale using variances
# - check if Sigma.hat is PD; if not, try again
#
# YR 26 Feb 2024

lav_partable_random <- function(lavpartable = NULL,
                                # needed if we still need to compute bounds:
                                lavh1 = NULL, lavdata = NULL,
                                lavsamplestats = NULL, lavoptions = NULL) {
  # check if we have bounds; if not, add them
  lavpta <- lav_partable_attributes(lavpartable)
  if (is.null(lavpartable$lower) ||
    is.null(lavpartable$upper)) {
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
  }

  # replace -Inf/Inf by -1/1 * .Machine$double.eps (for runif)
  inf.idx <- which(lavpartable$lower < -1e+16)
  if (length(inf.idx) > 0L) {
    lavpartable$lower[inf.idx] <- -1e+16
  }
  inf.idx <- which(lavpartable$upper > 1e+16)
  if (length(inf.idx) > 0L) {
    lavpartable$upper[inf.idx] <- 1e+16
  }

  # empty lavpartable$start?
  if (is.null(lavpartable$start)) {
    START <- numeric(length(lavpartable$lhs))
    # set loadings to 0.7
    loadings.idx <- which(lavpartable$free > 0L &
      lavpartable$op == "=~")
    if (length(loadings.idx) > 0L) {
      START[loadings.idx] <- 0.7
    }
    # set (only) variances to 1
    var.idx <- which(lavpartable$free > 0L &
      lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs)
    if (length(var.idx) > 0L) {
      START[var.idx] <- 1
    }

    lavpartable$start <- START
  }

  # initial values
  START <- lavpartable$start

  nblocks <- lav_partable_nblocks(lavpartable)
  block.values <- lav_partable_block_values(lavpartable)
  for (b in 1:nblocks) {
    ov.names <- lavpta$vnames$ov[[b]]
    lv.names <- lavpta$vnames$lv[[b]]
    ov.ind.names <- lavpta$vnames$ov.ind[[b]]

    # start with the lv (residual) variances
    lv.var.idx <- which(lavpartable$block == block.values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% lv.names &
      lavpartable$rhs %in% lv.names &
      lavpartable$lhs == lavpartable$rhs)

    if (length(lv.var.idx) > 0L) {
      for (i in lv.var.idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          START[i] <- runif(
            n = 1L, min = lavpartable$lower[i],
            max = lavpartable$upper[i]
          )
        }
      }
    }

    # first, we generate lv correlations, and then rescale to covariances
    lv.cov.idx <- which(lavpartable$block == block.values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% lv.names &
      lavpartable$rhs %in% lv.names &
      lavpartable$lhs != lavpartable$rhs)

    if (length(lv.cov.idx) > 0L) {
      for (i in lv.cov.idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          cor.val <- runif(n = 1L, -0.5, +0.5)
          var1.idx <- which(lavpartable$block == block.values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$lhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          var2.idx <- which(lavpartable$block == block.values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          START[i] <- cor.val * sqrt(START[var1.idx]) * sqrt(START[var2.idx])
        }
      }
    }

    # next, (residual) ov variances
    ov.var.idx <- which(lavpartable$block == block.values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% ov.names &
      lavpartable$rhs %in% ov.names &
      lavpartable$lhs == lavpartable$rhs)

    if (length(ov.var.idx) > 0L) {
      for (i in ov.var.idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          START[i] <- runif(
            n = 1L, min = lavpartable$lower[i],
            max = lavpartable$upper[i]
          )
        }
      }
    }

    # (residual) ov covariances (if any)
    ov.cov.idx <- which(lavpartable$block == block.values[b] &
      lavpartable$op == "~~" &
      lavpartable$lhs %in% ov.names &
      lavpartable$rhs %in% ov.names &
      lavpartable$lhs != lavpartable$rhs)

    if (length(ov.cov.idx) > 0L) {
      for (i in ov.cov.idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          cor.val <- runif(n = 1L, -0.5, +0.5)
          var1.idx <- which(lavpartable$block == block.values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$lhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          var2.idx <- which(lavpartable$block == block.values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          START[i] <- cor.val * sqrt(START[var1.idx]) * sqrt(START[var2.idx])
        }
      }
    }

    # finally, the lambda values, keeping in mind that
    # lambda_p^(u) = sqrt( upper(res.var.indicators_p) /
    #                      lower(var.factor) )
    lambda.idx <- which(lavpartable$block == block.values[b] &
      lavpartable$op == "=~" &
      lavpartable$lhs %in% lv.names &
      lavpartable$rhs %in% ov.ind.names)

    if (length(lambda.idx)) {
      for (i in lambda.idx) {
        if (lavpartable$free[i] > 0L &&
          (lavpartable$lower[i] < lavpartable$upper[i])) {
          varov.idx <- which(lavpartable$block == block.values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          varlv.idx <- which(lavpartable$block == block.values[b] &
            lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$lhs[i] &
            lavpartable$lhs == lavpartable$rhs)
          lambda.u <- sqrt(START[varov.idx] / START[varlv.idx])
          START[i] <- runif(n = 1, -lambda.u, lambda.u)
        }
      }
    }
  }

  # sanity check; needed?
  START <- lav_start_check_cov(
    lavpartable = lavpartable, start = START,
    warn = TRUE
  )

  START
}
