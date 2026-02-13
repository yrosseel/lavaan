# IV/MIIV estimation

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_miiv_internal <- function(lavmodel = NULL, lavh1 = NULL,
                                  lavsamplestats = NULL,
                                  lavpartable = NULL,
                                  lavdata = NULL, lavoptions = NULL) {
  # this is the entry-point for MIIV estimation
  lav_msg_note(gettext("** Estimator MIIV is still under development! **\n"))

  # IV options
  iv.method <- toupper(lavoptions$estimator.args$iv.method)
  stopifnot(iv.method %in% "2SLS")
  iv.varcov.method <- toupper(lavoptions$estimator.args$iv.varcov.method)
  stopifnot(iv.varcov.method %in% c("ULS", "GLS", "2RLS", "RLS"))

  # we assume the blocks are independent groups for now
  stopifnot(lavdata@nlevels == 1L)

  # parameter vector -> all NA (for the free parameters)
  x <- lav_model_get_parameters(lavmodel) * as.numeric(NA)

  ###############
  # first stage #
  ###############
  x <- lav_sem_miiv_2sls_rawdata(
    lavmodel = lavmodel, lavpartable = lavpartable,
    lavdata = lavdata, x = x
  )
  eqs <- attr(x, "eqs")
  x <- as.numeric(x) # drop attributes

  ################
  # second stage #
  ################
  x <- lav_sem_miiv_varcov(
    lavmodel = lavmodel, lavpartable = lavpartable,
    lavh1 = lavh1,
    x = x, iv.varcov.method = iv.varcov.method
  )

  attr(x, "eqs") <- eqs
  x
}


lav_sem_miiv_2sls_rawdata <- function(lavmodel = NULL, lavpartable = NULL,
                                      eqs = NULL, lavdata = NULL, x = NULL) {
  # get lavpta
  lavpta <- lav_partable_attributes(lavpartable)

  if (is.null(eqs)) {
    # find ALL model-implied instrumental variables (miivs) per equation
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpta = lavpta)
  }
  if (is.null(x)) {
    x <- lav_model_get_parameters(lavmodel, type = "free")
  }

  # number of blocks
  nblocks <- lavmodel@nblocks

  # parameter vector -> all NA (for the free parameters)
  x <- lav_model_get_parameters(lavmodel) * as.numeric(NA)

  # 2SLS per equation
  # for now: - no equality constraints (yet)
  #          - only OLS (not robust, no lasso, ...)
  for (b in seq_len(nblocks)) {

    # ov.names for this block
    ov.names <- lavpta$vnames$ov[[b]]

    # raw data for this block
    XY <- lavdata@X[[b]]

    # estimation per equation
    for (j in seq_along(eqs[[b]])) {

      # this equation
      eq <- eqs[[b]][[j]]

      # iv_flag?
      iv_flag <- TRUE
      if (!is.null(eq$iv_flag)) {
        iv_flag <- eq$iv_flag
      } else {
        # if rhs_new matches miiv, iv_flag is FALSE
        if (identical(eq$rhs_new, eq$miiv)) {
          iv_flag <- FALSE
        }
      }

      # Y: there is always an y variable
      y.idx <- match(eq$lhs_new, ov.names)
      yvec <- XY[, y.idx, drop = TRUE] # y-variable (always scalar)

      # X: usually, there are x variables (apart from the "1")
      if (identical(eq$rhs_new, "1")) {
        x.idx <- integer(0L)
        xmat <- matrix(0, nrow = nrow(XY), ncol = 0L)
      } else {
        x.idx <- match(eq$rhs_new, ov.names)
        xmat <- XY[, x.idx, drop = FALSE]
      }

      # Z: instruments
      if (iv_flag) {
        i.idx <- match(eq$miiv, ov.names)
        imat <- cbind(1, XY[, i.idx, drop = FALSE]) # instruments
      }

      # weights
      weights <- lavdata@weights[[b]]

      # sargan vector
      sargan <- rep(as.numeric(NA), 3L)
      names(sargan) <- c("stat", "df", "pvalue")

      # 0. check
      if (iv_flag && length(eq$miiv) < length(eq$rhs)) {
        # what to do? skip, or proceed anyway?
        eqs[[b]][[j]]$df <- as.integer(NA)
        eqs[[b]][[j]]$vcov <- matrix(0, 0L, 0L)
        eqs[[b]][[j]]$sargan <- sargan
        next
      }

      # 1. regress x on instruments
      if (iv_flag) {
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
      } else {
        # just plain regression: x == xhat
        xhat <- xmat
      }

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
      # - eq$pt contains partable/user rows
      # - lavpartable$free[eq$pt] should give the free idx
      free.idx <- lavpartable$free[eq$pt]
      if (length(x.idx) > 0L) {
        if (all(free.idx > 0L)) {
          x[free.idx] <- fit_y_on_xhat$coefficients[-1]
        } else {
          # remove non-free elements
          zero.idx <- which(free.idx == 0L)
          x[free.idx[-zero.idx]] <- fit_y_on_xhat$coefficients[-1][-zero.idx]
        }
      }
      if (lavmodel@meanstructure) {
        free.int.idx <- lavpartable$free[eq$ptint]
        if (free.int.idx > 0L) {
          x[free.int.idx] <- fit_y_on_xhat$coefficients[1]
        }
      }

      # 4. resvar
      notna.idx <- unname(which(!is.na(fit_y_on_xhat$coefficients)))
      ycoef <- fit_y_on_xhat$coefficients[notna.idx]
      res <- yvec - drop(cbind(1, xmat)[, notna.idx, drop = FALSE] %*% ycoef)
      df_res <- fit_y_on_xhat$df.residual
      if (is.null(weights)) {
        sse <- sum(res * res)
      } else {
        sse <- sum(weights * res * res)
      }
      # resvar <- sse / df_res # what we should do...
      resvar <- sse / length(res)

      # 5. naive cov (for standard errors) (see summary.lm in base R)
      p1 <- 1L:fit_y_on_xhat$rank
      R <- chol2inv(fit_y_on_xhat$qr$qr[p1, p1, drop = FALSE])
      vcov <- R * resvar # scaled (for now)

      # 6. Sargan test (see summary.ivreg.R 363--371)
      if (iv_flag) {
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
      }

      # add info to eqs list
      eqs[[b]][[j]]$df <- df_res
      eqs[[b]][[j]]$vcov <- vcov
      eqs[[b]][[j]]$sargan <- sargan
    } # eqs

  } # nblocks

  # return partially filled in x
  attr(x, "eqs") <- eqs
  x
}

lav_sem_miiv_varcov <- function(lavmodel = NULL, lavpartable = NULL,
                                lavh1 = NULL,
                                x = NULL, iv.varcov.method = "RLS") {
  # get lavpta
  lavpta <- lav_partable_attributes(lavpartable)

  # we x
  stopifnot(!is.null(x))

  # number of blocks
  nblocks <- lavmodel@nblocks

  # preparations
  lavmodel.tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
  delta_block <- lav_model_delta(
    lavmodel = lavmodel.tmp,
    ceq.simple = lavmodel@ceq.simple.only
  )
  # across all groups! (possible equality constraints)
  undirected.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  free.idx <- unique(lavpartable$free[undirected.idx])

  # create block-diagonal W_2 and s
  w2_block <- vector("list", length = nblocks)
  s_block <- vector("list", length = nblocks)
  delta2_block <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    sample_cov <- lavh1$implied$cov[[b]]
    sample_mean <- lavh1$implied$mean[[b]]

    # delta2
    delta2_block[[b]] <- delta_block[[b]][, free.idx, drop = FALSE]
    if (iv.varcov.method == "GLS") {
      s.inv <- solve(sample_cov)
    } else { # ULS, or starting point for 2RLS/RLS
      s.inv <- diag(1, nrow = nrow(sample_cov))
    }
    w2_22 <- 0.5 * lav_matrix_duplication_pre_post(s.inv %x% s.inv)
    if (lavmodel@meanstructure) {
      w2_11 <- s.inv
      w2_block[[b]] <- lav_matrix_bdiag(w2_11, w2_22)
      s_block[[b]] <- c(sample_mean, lav_matrix_vech(sample_cov))
    } else {
      w2_block[[b]] <- w2_22
      s_block[[b]] <- lav_matrix_vech(sample_cov)
    }
  }
  Delta2 <- do.call("rbind", delta2_block)
  W2 <- lav_matrix_bdiag(w2_block)
  svech <- unlist(s_block)

  # initial estimate for theta.2
  # TODO: if any constraints are needed, we need to augment the information
  theta2 <- drop(solve(
    t(Delta2) %*% W2 %*% Delta2,
    t(Delta2) %*% W2 %*% svech
  ))

  # insert theta2 element in x
  x[free.idx] <- theta2

  # again, but now with Sigma
  if (iv.varcov.method == "2RLS") {
    lavmodel.tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
    sigma <- lav_model_sigma(lavmodel = lavmodel.tmp, extra = FALSE)
    w2_block <- vector("list", length = nblocks)
    for (b in seq_len(nblocks)) {
      # delta2
      s.inv <- solve(sigma[[b]])
      w2_22 <- 0.5 * lav_matrix_duplication_pre_post(s.inv %x% s.inv)
      if (lavmodel@meanstructure) {
        w2_11 <- s.inv
        w2_block[[b]] <- lav_matrix_bdiag(w2_11, w2_22)
      } else {
        w2_block[[b]] <- w2_22
      }
    }
    W2 <- lav_matrix_bdiag(w2_block)
    # final estimate for theta.2
    # TODO: if any constraints are needed, we need to augment the information
    theta2 <- drop(solve(
      t(Delta2) %*% W2 %*% Delta2,
      t(Delta2) %*% W2 %*% svech
    ))
  } else if (iv.varcov.method == "RLS") {
    for (i in 1:200) {
      old_x <- theta2
      lavmodel.tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
      sigma <- lav_model_sigma(lavmodel = lavmodel.tmp, extra = FALSE)
      w2_block <- vector("list", length = nblocks)
      for (b in seq_len(nblocks)) {
        # delta2
        s.inv <- solve(sigma[[b]])
        w2_22 <- 0.5 * lav_matrix_duplication_pre_post(s.inv %x% s.inv)
        if (lavmodel@meanstructure) {
          w2_11 <- s.inv
          w2_block[[b]] <- lav_matrix_bdiag(w2_11, w2_22)
        } else {
          w2_block[[b]] <- w2_22
        }
      }
      W2 <- lav_matrix_bdiag(w2_block)
      # final estimate for theta.2
      # TODO: if any constraints are needed, we need to augment the information
      theta2 <- drop(solve(
        t(Delta2) %*% W2 %*% Delta2,
        t(Delta2) %*% W2 %*% svech
      ))
      x[free.idx] <- theta2
      sse <- sum((old_x - theta2)^2)
      # cat("i = ", i, "sse = ", sse, "\n")
      if (sse < 1e-05) {
        break
      } else if (i == 200L) {
        lav_msg_warn(gettext("RLS for variances/covariances did not converge
                              after 200 iterations."))
      }
    }
  }

  # insert final theta2 element in x
  x[free.idx] <- theta2

  # return x
  x
}
