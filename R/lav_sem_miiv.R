# IV/MIIV estimation

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_miiv_internal <- function(lavmodel = NULL, lavh1 = NULL,
                                  lavsamplestats = NULL,
                                  lavpartable = NULL,
                                  lavdata = NULL, lavoptions = NULL) {
  # IV options
  iv.method <- toupper(lavoptions$estimator.args$iv.method)
  stopifnot(iv.method %in% "2SLS")
  iv.varcov.method <- toupper(lavoptions$estimator.args$iv.varcov.method)
  iv.samplestats <- lavoptions$estimator.args$iv.samplestats
  # just in case
  if (lavmodel@categorical) {
    iv.samplestats <- TRUE
  }
  stopifnot(iv.varcov.method %in% c("ULS", "GLS", "2RLS", "RLS"))

  # get lavpta
  lavpta <- lav_partable_attributes(lavpartable)

  # we assume the blocks are independent groups for now
  stopifnot(lavdata@nlevels == 1L)

  # directed versus undirected (free) parameters
  undirected.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  directed.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~", "~1"))
  free.directed.idx <- unique(lavpartable$free[directed.idx])
  free.undirected.idx <- unique(lavpartable$free[undirected.idx])

  # find ALL model-implied instrumental variables (miivs) per equation
  eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpta = lavpta)

  # parameter vector -> all NA (for the free parameters)
  x <- lav_model_get_parameters(lavmodel)

  ########################################
  # first stage: directed parameter only #
  ########################################
  theta1 <- numeric(0L)
  if (length(free.directed.idx) > 0L) {
    if (iv.samplestats) {
      theta1 <- lav_sem_miiv_2sls_samplestats(
        x = NULL, samplestats = FALSE, eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavdata = lavdata, lavsamplestats = lavsamplestats,
        lavh1 = lavh1, free.directed.idx = free.directed.idx
      )
    } else {
      theta1 <- lav_sem_miiv_2sls(
        eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavdata = lavdata, free.directed.idx = free.directed.idx
      )
    }
    # update equations
    eqs <- attr(theta1, "eqs")
    theta1 <- as.numeric(theta1) # drop attributes
    # store theta1 elements in x
    x[free.directed.idx] <- theta1
  }

  #######################################
  # second stage: undirected parameters #
  #######################################

  # compute theta2 using ULS/GLS/RLS/2RLS
  theta2 <- numeric(0L)
  if (length(free.undirected.idx) > 0L) {
    theta2 <- lav_sem_miiv_varcov(
      x = theta1,
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats,
      lavh1 = lavh1, free.directed.idx = free.directed.idx,
      free.undirected.idx = free.undirected.idx,
      iv.varcov.method = iv.varcov.method
    )
    # store theta2 elements in x
    x[free.undirected.idx] <- theta2
  }

  attr(x, "eqs") <- eqs
  x
}


# stage 1: use 2SLS to find the regression coefficients of all
#          directed effects in the model -- continous/raw-data version
lav_sem_miiv_2sls <- function(eqs = NULL, lavmodel = NULL, lavpartable = NULL,
                              lavdata = NULL, free.directed.idx = NULL) {
  # this function is for continuous/raw-data only
  stopifnot(lavdata@data.type == "full")
  stopifnot(!lavmodel@categorical)

  # get lavpta
  lavpta <- lav_partable_attributes(lavpartable)

  if (is.null(eqs)) {
    # find ALL model-implied instrumental variables (miivs) per equation
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpta = lavpta)
  }

  # number of blocks
  nblocks <- lavmodel@nblocks

  # parameter vector
  x <- lav_model_get_parameters(lavmodel)

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
        eqs[[b]][[j]]$coef <- numeric(0L)
        eqs[[b]][[j]]$nobs <- nrow(xmat)
        eqs[[b]][[j]]$df_res <- as.integer(NA)
        eqs[[b]][[j]]$resvar <- as.numeric(NA)
        eqs[[b]][[j]]$resvar_df_res <- as.numeric(NA)
        eqs[[b]][[j]]$XX <- crossprod(cbind(1, xmat))
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
          free.idx <- free.idx[-zero.idx]
          if (length(free.idx) > 0) {
            x[free.idx] <- fit_y_on_xhat$coefficients[-1][-zero.idx]
          }
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
      resvar_df_res <- sse / df_res # what we should do...
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
      eqs[[b]][[j]]$coef <- unname(fit_y_on_xhat$coefficients)
      eqs[[b]][[j]]$nobs <- nrow(xmat)
      eqs[[b]][[j]]$df_res <- df_res
      eqs[[b]][[j]]$resvar <- resvar
      eqs[[b]][[j]]$resvar_df_res <- resvar_df_res
      eqs[[b]][[j]]$XX <- crossprod(cbind(1, xmat)) # needed? or only vcov?
      eqs[[b]][[j]]$vcov <- vcov
      eqs[[b]][[j]]$sargan <- sargan
    } # eqs
  } # nblocks

  # return theta1
  theta1 <- x[free.directed.idx]

  # add equations as an attribute
  attr(theta1, "eqs") <- eqs

  theta1
}

# stage 1: use 2SLS to find the regression coefficients of all
#          directed effects in the model -- samplestats version
# if samplestats = TRUE, x should contain the sample statistics
# (taken from lavh1$implied)
lav_sem_miiv_2sls_samplestats <- function(x = NULL, samplestats = FALSE,
                                          eqs = NULL, lavmodel = NULL,
                                          lavpartable = NULL, lavdata = NULL,
                                          lavsamplestats = NULL, lavh1 = NULL,
                                          free.directed.idx = NULL) {
  # helper function: first chol, then eigen if chol fails
  solve_spd <- function(A, B, tol = 1e-10) {
    cholA <- try(chol(A), silent = TRUE)
    if (!inherits(cholA, "try-error")) {
      # Cholesky solve
      return(backsolve(cholA, forwardsolve(t(cholA), B)))
    } else {
      # Eigen fallback (pseudo-inverse)
      eig <- eigen(A, symmetric = TRUE)
      keep <- eig$values > tol * max(eig$values)
      Ainv <- eig$vectors[, keep, drop = FALSE] %*%
        diag(1 / eig$values[keep], length(eig$values[keep])) %*%
        t(eig$vectors[, keep, drop = FALSE])
      return(Ainv %*% B)
    }
  }

  # no conditional.x for now!
  stopifnot(!lavmodel@conditional.x)

  # get lavpta
  lavpta <- lav_partable_attributes(lavpartable)

  if (is.null(eqs)) {
    # find ALL model-implied instrumental variables (miivs) per equation
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpta = lavpta)
  }

  if (samplestats) {
    implied <- lav_vec_to_implied(x, lavmodel = lavmodel)
  } else {
    implied <- lavh1$implied
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

    # sample statistics for this block
    sample.cov <- implied$cov[[b]]
    sample.mean <- implied$mean[[b]]

    # nobs for this 'block'
    nobs <- lavsamplestats@nobs[[b]]

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
      y.bar <- sample.mean[y.idx]

      # X: usually, there are x variables (apart from the "1")
      if (identical(eq$rhs_new, "1")) {
        x.idx <- integer(0L)
      } else {
        x.idx <- match(eq$rhs_new, ov.names)
        x.bar <- sample.mean[x.idx]
      }

      # Z: instruments
      S_Xy <- sample.cov[x.idx, y.idx, drop = FALSE]
      S_XX <- sample.cov[x.idx, x.idx, drop = FALSE]
      s_yy <- sample.cov[y.idx, y.idx, drop = FALSE]
      if (iv_flag) {
        i.idx <- match(eq$miiv, ov.names)
        S_XZ <- sample.cov[x.idx, i.idx, drop = FALSE]
        S_ZX <- sample.cov[i.idx, x.idx, drop = FALSE]
        S_ZZ <- sample.cov[i.idx, i.idx, drop = FALSE]
        S_Zy <- sample.cov[i.idx, y.idx, drop = FALSE]
      }

      # sargan vector
      sargan <- rep(as.numeric(NA), 3L)
      names(sargan) <- c("stat", "df", "pvalue")

      # 0. check
      if (iv_flag && length(eq$miiv) < length(eq$rhs)) {
        # what to do? skip, or proceed anyway?
        eqs[[b]][[j]]$coef <- numeric(0L)
        eqs[[b]][[j]]$nobs <- nobs
        eqs[[b]][[j]]$df_res <- as.integer(NA)
        eqs[[b]][[j]]$resvar <- as.numeric(NA)
        eqs[[b]][[j]]$resvar_df_res <- as.numeric(NA)
        eqs[[b]][[j]]$XX <- matrix(0, 0L, 0L)
        eqs[[b]][[j]]$vcov <- matrix(0, 0L, 0L)
        eqs[[b]][[j]]$sargan <- sargan
        next
      }

      fit_y_on_xhat <- list()
      if (iv_flag) {
        # Step 1: compute S_ZZ^{-1} S_ZX and S_ZZ^{-1} S_Zy
        W <- solve_spd(S_ZZ, S_ZX)
        v <- solve_spd(S_ZZ, S_Zy)
        # Step 2: build reduced system
        Amat <- S_XZ %*% W
        bvec <- S_XZ %*% v
        beta_slopes <- drop(solve_spd(Amat, bvec))
        beta0 <- as.vector(y.bar - t(x.bar) %*% beta_slopes)
      } else {
        if (length(x.idx) > 0L) {
          Amat <- S_XX
          bvec <- S_Xy
          beta_slopes <- drop(solve_spd(Amat, bvec))
          beta0 <- as.vector(y.bar - t(x.bar) %*% beta_slopes)
        } else {
          beta_slopes <- numeric(0L)
          beta0 <- y.bar
        }
      }
      fit_y_on_xhat$coefficients <- c(beta0, beta_slopes)

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
          free.idx <- free.idx[-zero.idx]
          if (length(free.idx) > 0) {
            x[free.idx] <- fit_y_on_xhat$coefficients[-1][-zero.idx]
          }
        }
      }
      if (lavmodel@meanstructure) {
        free.int.idx <- lavpartable$free[eq$ptint]
        if (free.int.idx > 0L) {
          x[free.int.idx] <- fit_y_on_xhat$coefficients[1]
        }
      }

      # 4. resvar
      if (length(x.idx) > 0L) {
        tmp <- S_XX %*% beta_slopes
        resvar <- as.numeric(s_yy - 2 * t(beta_slopes) %*% S_Xy +
          t(beta_slopes) %*% tmp)
      } else {
        resvar <- drop(s_yy)
      }
      df_res <- nobs - (length(x.idx) + 1L)
      resvar_df_res <- resvar * nobs / df_res

      # 5. naive cov (for standard errors) (see summary.lm in base R)
      vcov <- NULL
      if (!lavmodel@categorical) {
        if (length(x.idx) > 0L) {
          Ainv <- solve_spd(Amat, diag(x = 1, nrow = nrow(Amat)))
          vcov.slopes <- (resvar / nobs) * Ainv
          vcov.beta0 <- (resvar / nobs) + t(x.bar) %*% vcov.slopes %*% x.bar
          cov.beta0.slopes <- -vcov.slopes %*% x.bar
          vcov <- matrix(0, length(x.idx) + 1L, length(x.idx) + 1L)
          vcov[1, 1] <- vcov.beta0
          vcov[1, -1] <- t(cov.beta0.slopes)
          vcov[-1, 1] <- cov.beta0.slopes
          vcov[-1, -1] <- vcov.slopes
        } else {
          # only intercept
          vcov <- matrix(resvar / nobs, 1L, 1L)
        }
      }

      # 6. Sargan test (see summary.ivreg.R 363--371)
      if (iv_flag) {
        sargan["df"] <- length(i.idx) - length(x.idx)
        if (sargan["df"] > 0L) {
          g <- S_Zy - S_ZX %*% beta_slopes
          w <- solve_spd(S_ZZ, g)
          sargan["stat"] <- as.numeric(nobs * t(g) %*% w / resvar)
          if (!lavmodel@categorical) {
            sargan["pvalue"] <- pchisq(sargan["stat"], sargan["df"],
              lower.tail = FALSE
            )
          }
        }
      }

      # XX (needed?)
      XX <- NULL
      if (!lavmodel@categorical) {
        if (length(x.idx) > 0L) {
          XX <- tcrossprod(c(1, x.bar)) * nobs
          XX[-1, -1] <- (S_XX + tcrossprod(x.bar)) * nobs
        } else {
          XX <- matrix(nobs, 1L, 1L)
        }
      }

      # add info to eqs list
      eqs[[b]][[j]]$coef <- unname(fit_y_on_xhat$coefficients)
      eqs[[b]][[j]]$nobs <- nobs
      eqs[[b]][[j]]$df_res <- df_res
      eqs[[b]][[j]]$resvar <- resvar
      eqs[[b]][[j]]$resvar_df_res <- resvar_df_res
      eqs[[b]][[j]]$XX <- XX
      eqs[[b]][[j]]$vcov <- vcov
      eqs[[b]][[j]]$sargan <- sargan
    } # eqs
  } # nblocks

  # return theta1
  theta1 <- x[free.directed.idx]

  # add equations as an attribute
  if (!samplestats) {
    attr(theta1, "eqs") <- eqs
  }

  theta1
}


# by default: input (x) is theta1 (only)
# BUT if samplestats is TRUE, then x contains the sample statisics
# this allows us to compute the jacobian wrt theta1 elements (keeping
# sample statistics fixed, or the jacobian wrt the sample statistics
# (keeping theta1 fixed)
lav_sem_miiv_varcov <- function(x = NULL, samplestats = FALSE,
                                lavmodel = NULL, lavpartable = NULL,
                                lavsamplestats = NULL,
                                lavh1 = NULL, free.directed.idx = NULL,
                                free.undirected.idx = NULL,
                                iv.varcov.method = "RLS") {
  if (samplestats) {
    implied <- lav_vec_to_implied(x, lavmodel = lavmodel)
    x <- lav_model_get_parameters(lavmodel)
    theta1 <- x[free.directed.idx]
  } else {
    theta1 <- x
    x <- lav_model_get_parameters(lavmodel)
    x[free.directed.idx] <- theta1
    implied <- lavh1$implied
  }

  # number of blocks
  nblocks <- lavmodel@nblocks

  # preparations
  lavmodel.tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
  delta_block <- lav_model_delta(
    lavmodel = lavmodel.tmp,
    ceq.simple = lavmodel@ceq.simple.only
  )

  # create block-diagonal W_2 and s
  w2_block <- vector("list", length = nblocks)
  s_block <- vector("list", length = nblocks)
  delta2_block <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    # delta2
    delta2_block[[b]] <- delta_block[[b]][, free.undirected.idx, drop = FALSE]

    if (lavmodel@categorical) {
      # if (iv.varcov.method == "ULS") {
      w2_block[[b]] <- diag(1, nrow = length(lavsamplestats@WLS.obs[[b]]))
      # } else {
      #  w2_block[[b]] <- diag(lavsamplestats@WLS.VD[[b]])
      # }
      # s_block[[b]] <- lavsamplestats@WLS.obs[[b]]
      s_block[[b]] <- lav_implied_to_vec(implied = implied, lavmodel = lavmodel, drop.list = FALSE)[[b]]
    } else {
      # continuous case
      sample_cov <- implied$cov[[b]]
      sample_mean <- implied$mean[[b]]

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
  } # blocks

  Delta2 <- do.call("rbind", delta2_block)
  W2 <- lav_matrix_bdiag(w2_block)
  svec <- unlist(s_block)

  # initial estimate for theta.2
  # TODO: if any constraints are needed, we need to augment the information
  theta2 <- drop(solve(
    t(Delta2) %*% W2 %*% Delta2,
    t(Delta2) %*% W2 %*% svec
  ))

  # in the categorical case, we are done
  if (lavmodel@categorical) {
    return(theta2)
  }

  # again, but now with Sigma
  if (iv.varcov.method == "2RLS") {
    x[free.undirected.idx] <- theta2
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
      t(Delta2) %*% W2 %*% svec
    ))
  } else if (iv.varcov.method == "RLS") {
    for (i in 1:200) {
      old_x <- theta2
      x[free.undirected.idx] <- theta2
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
        t(Delta2) %*% W2 %*% svec
      ))
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

  theta2
}


# VCOV for free parameters
lav_sem_miiv_vcov <- function(lavmodel = NULL, lavsamplestats = NULL,
                              lavoptions = NULL, lavpartable = NULL,
                              lavimplied = NULL,
                              lavh1 = NULL, eqs = NULL) {
  # iv options
  iv.varcov.se <- lavoptions$estimator.args$iv.varcov.se
  iv.varcov.modelbased <- lavoptions$estimator.args$iv.varcov.modelbased
  iv.varcov.method <- toupper(lavoptions$estimator.args$iv.varcov.method)

  # empty vcov
  vcov <- matrix(0, lavmodel@nx.free, lavmodel@nx.free)

  # nblocks
  nblocks <- lavmodel@nblocks

  # eqs?
  if (is.null(eqs)) {
    lavpta <- lav_partable_attributes(lavpartable)
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpta = lavpta)
  }

  # directed versus undirected (free) parameters
  undirected.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  directed.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~", "~1"))
  free.directed.idx <- unique(lavpartable$free[directed.idx])
  free.undirected.idx <- unique(lavpartable$free[undirected.idx])

  # directed free parameters only
  x <- lav_model_get_parameters(lavmodel, type = "free")
  theta1 <- x[free.directed.idx]

  # do we need gamma_big?
  if (lavmodel@categorical ||
    (iv.varcov.se && length(free.undirected.idx) > 0L)) {
    gamma_g <- vector("list", lavmodel@ngroups)
    for (g in seq_len(lavmodel@ngroups)) {
      if (!is.null(lavsamplestats@NACOV[[g]])) {
        gamma_g[[g]] <- lavsamplestats@NACOV[[g]]
      } else {
        if (iv.varcov.modelbased) {
          mean_g <- lavimplied$mean[[g]]
          cov_g <- lavimplied$cov[[g]]
        } else {
          mean_g <- lavh1$implied$mean[[g]]
          cov_g <- lavh1$implied$cov[[g]]
        }
        # NT version (for now), model-based
        gamma_g[[g]] <- lav_samplestats_Gamma_NT(
          COV = cov_g,
          MEAN = mean_g,
          x.idx = lavsamplestats@x.idx[[g]],
          fixed.x = lavmodel@fixed.x,
          conditional.x = lavmodel@conditional.x,
          meanstructure = lavmodel@meanstructure,
          slopestructure = lavmodel@conditional.x
        )
        # divide by (group) sample size
        gamma_g[[g]] <- gamma_g[[g]] / lavsamplestats@nobs[[g]]
      }
    }
    gamma_big <- lav_matrix_bdiag(gamma_g)
  }

  # stage 1: directed effects
  if (length(free.directed.idx) > 0L) {
    if (lavmodel@categorical) {
      # compute K matrix
      vec <- lav_implied_to_vec(
        implied = lavh1$implied, lavmodel = lavmodel,
        drop.list = TRUE
      )
      jac_k <- numDeriv::jacobian(
        lav_sem_miiv_2sls_samplestats,
        x = vec, samplestats = TRUE, eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavsamplestats = lavsamplestats,
        lavh1 = lavh1, free.directed.idx = free.directed.idx
      )
      vcov[free.directed.idx, free.directed.idx] <-
        (jac_k %*% gamma_big %*% t(jac_k)) / lavsamplestats@ntotal
    } else {
      for (b in seq_len(nblocks)) {
        neqs <- length(eqs[[b]])
        for (j in seq_len(neqs)) {
          eq <- eqs[[b]][[j]]

          eq_vcov <- NULL
          if (!is.null(eq$vcov) && nrow(eq$vcov) > 0L) {
            eq_vcov <- eq$vcov # always includes intercept
          } else if (!is.null(eq$XX) && !is.null(eq$resvar)) {
            # reconstruct using XX and resvar (less stable!)
            eq_vcov <- solve(crossprod(eq$XX)) * eq$resvar # or resvar_df_res
          }
          if (is.null(eq_vcov)) {
            next
          }

          # only keep free parameters
          free.idx <- lavpartable$free[eq$pt]
          if (length(eq$rhs) > 0L) {
            if (all(free.idx > 0L)) {
              vcov[free.idx, free.idx] <- eq_vcov[-1, -1]
            } else {
              # remove non-free elements
              zero.idx <- which(free.idx == 0L)
              free.idx <- free.idx[-zero.idx]
              if (length(free.idx) > 0L) {
                vcov[free.idx, free.idx] <-
                  eq_vcov[-1, -1][-zero.idx, -zero.idx, drop = FALSE]
              }
            }
          }
          if (lavmodel@meanstructure) {
            free.int.idx <- lavpartable$free[eq$ptint]
            if (free.int.idx > 0L) {
              vcov[free.int.idx, free.int.idx] <- eq_vcov[1, 1]
            }
          }
        } # neqs
      } # nblocks
    } # continuous/vcov version
  }

  # stage 2: undirected effects (note: intercepts have no effect!)

  if (iv.varcov.se && length(free.undirected.idx) > 0L) {
    if (lavmodel@categorical) {
      delta_list <- lav_model_delta(lavmodel) # per block/group
      delta <- do.call("rbind", delta_list)
      delta1 <- delta[, free.directed.idx, drop = FALSE]
      delta2 <- delta[, free.undirected.idx, drop = FALSE]
      h2 <- solve(crossprod(delta2), t(delta2))
      idiag <- diag(1, nrow = nrow(delta1))
      if (length(free.directed.idx) > 0L) {
        tmp <- h2 %*% (idiag - delta1 %*% jac_k)
      } else {
        tmp <- h2
      }
      vcov[free.undirected.idx, free.undirected.idx] <-
        (tmp %*% gamma_big %*% t(tmp)) / lavsamplestats@ntotal
    } else {
      # part a: effect of theta1
      jac_a <- try(lav_func_jacobian_complex(
        func = lav_sem_miiv_varcov,
        x = theta1, samplestats = FALSE, lavmodel = lavmodel,
        lavpartable = lavpartable, lavsamplestats = lavsamplestats,
        lavh1 = lavh1, free.directed.idx = free.directed.idx,
        free.undirected.idx = free.undirected.idx,
        iv.varcov.method = iv.varcov.method
      ), silent = TRUE)
      if (inherits(jac_a, "try-error")) {
        jac_a <- numDeriv::jacobian(
          func = lav_sem_miiv_varcov,
          x = theta1, samplestats = FALSE, lavmodel = lavmodel,
          lavpartable = lavpartable, lavsamplestats = lavsamplestats,
          lavh1 = lavh1, free.directed.idx = free.directed.idx,
          free.undirected.idx = free.undirected.idx,
          iv.varcov.method = iv.varcov.method
        )
      }
      vcov_directed <- vcov[free.directed.idx, free.directed.idx, drop = FALSE]
      vcov_a <- jac_a %*% vcov_directed %*% t(jac_a)

      # part b: effect of sample statistics
      # - get jacobian lav_sem_miiv_varcov wrt implied_vec
      # - get Gamma (NT or ADF)
      # - vcov_undirected_b <- jac_b %*% gamma_mat %*% t(jac_b)
      vec <- lav_implied_to_vec(
        implied = lavh1$implied, lavmodel = lavmodel,
        drop.list = TRUE
      )
      jac_b <- try(lav_func_jacobian_complex(
        func = lav_sem_miiv_varcov,
        x = vec, samplestats = TRUE, lavmodel = lavmodel,
        lavpartable = lavpartable, lavsamplestats = lavsamplestats,
        lavh1 = lavh1, free.directed.idx = free.directed.idx,
        free.undirected.idx = free.undirected.idx,
        iv.varcov.method = iv.varcov.method
      ), silent = TRUE)
      if (inherits(jac_b, "try-error")) {
        jac_b <- numDeriv::jacobian(
          func = lav_sem_miiv_varcov,
          x = vec, samplestats = TRUE, lavmodel = lavmodel,
          lavpartable = lavpartable, lavsamplestats = lavsamplestats,
          lavh1 = lavh1, free.directed.idx = free.directed.idx,
          free.undirected.idx = free.undirected.idx,
          iv.varcov.method = iv.varcov.method
        )
      }

      vcov_b <- (jac_b %*% gamma_big %*% t(jac_b)) / lavsamplestats@ntotal
      vcov_ab <- vcov_a + vcov_b
      vcov[free.undirected.idx, free.undirected.idx] <- vcov_ab
    } # continuous
  } # iv.varcov.se = TRUE

  vcov
}
