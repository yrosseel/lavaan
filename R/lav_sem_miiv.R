# place-holder for MIIV estimation

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_miiv_internal <- function(lavmodel = NULL, lavsamplestats = NULL,
                                  lavpartable = NULL,
                                  lavdata = NULL, lavoptions = NULL) {
  # this is the entry-point for MIIV estimation
  lav_msg_note(gettext("** Estimator MIIV is still under development! **\n"))

  lavpta <- lav_partable_attributes(lavpartable)
  nblocks <- lavmodel@nblocks
  # we assume the blocks are independent groups for now
  stopifnot(lavdata@nlevels == 1L)

  # find ALL model-implied instrumental variables (miivs) per equation
  eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpta = lavpta)

  # parameter vector -> all NA (for the free parameters)
  x <- lav_model_get_parameters(lavmodel) * as.numeric(NA)

  # 2SLS per equation
  # for now: - no equality constraints (yet)
  #          - only OLS (not robust, no lasso, ...)
  for (b in seq_len(nblocks)) {
    ov.names <- lavpta$vnames$ov[[b]]
    XY <- lavdata@X[[b]]
    for (j in seq_along(eqs[[b]])) {
      eq <- eqs[[b]][[j]]

      # variable indices
      y.idx <- match(eq$lhs_new, ov.names)
      x.idx <- match(eq$rhs_new, ov.names)
      i.idx <- match(eq$miiv, ov.names)

      # matrices
      imat <- cbind(1, XY[, i.idx, drop = FALSE]) # instruments
      xmat <- XY[, x.idx, drop = FALSE] # x-variables
      yvec <- XY[, y.idx, drop = TRUE] # y-variable (always scalar)
      weights <- lavdata@weights[[b]]

      # sargan vector
      sargan <- rep(as.numeric(NA), 3L)
      names(sargan) <- c("stat", "df", "pvalue")

      # 0. check
      if (length(eq$miiv) < length(eq$rhs)) {
        # what to do? skip, or proceed anyway?
        eqs[[b]][[j]]$df <- as.integer(NA)
        eqs[[b]][[j]]$vcov <- matrix(0, 0L, 0L)
        eqs[[b]][[j]]$sargan <- sargan
        next
      }

      # 1. regress x on instruments
      if (is.null(weights)) {
        fit_x_on_z <- lm.fit(x = imat, y = xmat)
      } else {
        fit_x_on_z <- lm.wfit(x = imat, y = xmat, weights = weights)
      }
      # check for NA's
      if (anyNA(fit_x_on_z$coefficients)) {
        lav_msg_warn(gettextf(
          "regression %s on %s failed (NAs);
                               redundant instruments?",
          paste(ov.names[x.idx], collapse = " + "),
          paste(ov.names[i.idx], collapse = " + ")
        ))
      }
      xhat <- as.matrix(fit_x_on_z$fitted.values)

      # 2. regress y on xhat
      if (is.null(weights)) {
        fit_y_on_xhat <- lm.fit(x = cbind(1, xhat), y = yvec)
      } else {
        fit_y_on_xhat <- lm.wfit(
          x = cbind(1, xhat), y = yvec,
          weights = weights
        )
      }

      # 3. fill estimates in x
      x[lavpartable$free[eq$pt]] <- fit_y_on_xhat$coefficients[-1]
      if (lavmodel@meanstructure) {
        x[lavpartable$free[eq$ptint]] <- fit_y_on_xhat$coefficients[1]
      }
      # 3b. take care of exogenous variables
      # TODO

      # 4. resvar
      notna.idx <- unname(which(!is.na(fit_y_on_xhat$coefficients)))
      ycoef <- fit_y_on_xhat$coefficients[notna.idx]
      res <- yvec - drop(cbind(1, xmat)[, notna.idx, drop = FALSE] %*% ycoef)
      df <- fit_y_on_xhat$df.residual
      if (is.null(weights)) {
        sse <- sum(res * res)
      } else {
        sse <- sum(weights * res * res)
      }
      #resvar <- sse / df # what we should do...
      resvar <- sse / length(res)

      # 5. naive cov (for standard errors) (see summary.lm in base R)
      p1 <- 1L:fit_y_on_xhat$rank
      R <- chol2inv(fit_y_on_xhat$qr$qr[p1, p1, drop = FALSE])
      vcov <- R * resvar # scaled (for now)

      # 6. Sargan test (see summary.ivreg.R 363--371)
      sargan["df"] <- length(i.idx) - length(x.idx)
      if (sargan["df"] > 0L) {
        if (is.null(weights)) {
          fit_yres_on_z <- lm.fit(x = imat, y = res)
          rssr <- sum((res - mean(res))^2)
          sse2 <- sum(fit_yres_on_z$residuals * fit_yres_on_z$residuals)
        } else {
          fit_yres_on_z <- lm.wfit(x = imat, y = res, weights = weights)
          rssr <- sum(weights * (res - weighted.mean(res, weights))^2)
          sse2 <- sum(weights * fit_yres_on_z$residuals *
                        fit_yres_on_z$residuals)
        }
        sargan["stat"] <- length(res) * (1 - sse2 / rssr)
        sargan["pvalue"] <- pchisq(sargan["stat"], sargan["df"],
          lower.tail = FALSE
        )
      }

      # add info to eqs list
      eqs[[b]][[j]]$df <- df
      eqs[[b]][[j]]$vcov <- vcov
      eqs[[b]][[j]]$sargan <- sargan
    } # eqs

    # take care exogenous variables
    if (lavmodel@meanstructure && lavoptions$marker.int.zero) {
      lv.x <- lavpta$vnames$lv.x[[b]]
      if (length(lv.x) > 0L) {
        target.id <- which(lavpartable$op == "~1" &
                             lavpartable$block == b &
                             lavpartable$free != 0 &
                             lavpartable$lhs %in% lv.x)
        ov.idx <- match(lavpta$vnames$lv.marker[[b]][lv.x], ov.names)
        ov.mean <- colMeans(XY[, ov.idx, drop = FALSE])
        x[lavpartable$free[target.id]] <- ov.mean
      }
    }

    # estimate variances/covariances in this block

  } # nblocks

  attr(x, "eqs") <- eqs
  x
}
