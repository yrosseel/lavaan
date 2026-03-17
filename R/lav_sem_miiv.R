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
  if (lavdata@data.type == "moment") {
    # force iv.samplestats = TRUE
    iv.samplestats <- TRUE
  }
  iv.vcov.stage1 <- tolower(lavoptions$estimator.args$iv.vcov.stage1)
  iv.vcov.stage2 <- tolower(lavoptions$estimator.args$iv.vcov.stage2)
  iv.sargan <- lavoptions$estimator.args$iv.sargan
  # just in case
  if (lavmodel@categorical) {
    iv.samplestats <- TRUE
  }
  stopifnot(iv.varcov.method %in% c("ULS", "GLS", "2RLS", "RLS", "NONE"))

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
  eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)

  # initial parameter vector
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
        lavh1 = lavh1, free.directed.idx = free.directed.idx,
        iv.vcov.stage1 = iv.vcov.stage1, iv.sargan = iv.sargan,
        iv.vcov.stage2 = iv.vcov.stage2
      )
    } else {
      theta1 <- lav_sem_miiv_2sls(
        eqs = eqs,
        lavmodel = lavmodel, lavpartable = lavpartable,
        lavdata = lavdata, free.directed.idx = free.directed.idx,
        iv.vcov.stage1 = iv.vcov.stage1, iv.sargan = iv.sargan,
        iv.vcov.stage2 = iv.vcov.stage2
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
  if (length(free.undirected.idx) > 0L && iv.varcov.method != "NONE") {
    theta2 <- lav_sem_miiv_varcov(
      x = theta1,
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavsamplestats = lavsamplestats,
      lavh1 = lavh1, free.directed.idx = free.directed.idx,
      free.undirected.idx = free.undirected.idx,
      iv.varcov.method = iv.varcov.method
    )
  } else {
    theta2 <- rep(as.numeric(NA), length(free.undirected.idx))
  }
  # store theta2 elements in x
  x[free.undirected.idx] <- theta2

  attr(x, "eqs") <- eqs
  x
}


# stage 1: use 2SLS to find the regression coefficients of all
#          directed effects in the model -- continous/raw-data version
lav_sem_miiv_2sls <- function(eqs = NULL, lavmodel = NULL, lavpartable = NULL,
                              lavdata = NULL, free.directed.idx = NULL,
                              iv.vcov.stage1 = "lm.vcov.dfres",
                              iv.vcov.stage2 = "delta",
                              iv.sargan = TRUE) {
  # this function is for continuous/raw-data only
  stopifnot(lavdata@data.type == "full")
  stopifnot(!lavmodel@categorical)
  iv.vcov.stage1 <- tolower(iv.vcov.stage1)
  stopifnot(iv.vcov.stage1 %in% c("lm.vcov.dfres", "lm.vcov", "none"))

  # get lavpta
  lavpta <- lav_partable_attributes(lavpartable)

  if (is.null(eqs)) {
    # find ALL model-implied instrumental variables (miivs) per equation
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
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
      if (!is.null(eq$iv_type)) {
        iv_flag <- eq$iv_type != "ols"
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
        i.idx <- match(eq$iv, ov.names)
        imat <- cbind(1, XY[, i.idx, drop = FALSE]) # instruments
      }

      # weights
      weights <- lavdata@weights[[b]]

      # sargan vector
      sargan <- rep(as.numeric(NA), 3L)
      names(sargan) <- c("stat", "df", "pvalue")

      # 0. check
      if (iv_flag && length(eq$iv) < length(eq$rhs)) {
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
      resvar_df_res <- resvar <- df_res <- NULL
      if (iv.vcov.stage1 %in% c("lm.vcov.dfres", "lm.vcov")) {
        df_res <- fit_y_on_xhat$df.residual
        notna.idx <- unname(which(!is.na(fit_y_on_xhat$coefficients)))
        ycoef <- fit_y_on_xhat$coefficients[notna.idx]
        res <- yvec - drop(cbind(1, xmat)[, notna.idx, drop = FALSE] %*% ycoef)
        if (is.null(weights)) {
          sse <- sum(res * res)
        } else {
          sse <- sum(weights * res * res)
        }
        resvar_df_res <- sse / df_res # what we should do...
        resvar <- sse / length(res)
      }

      # 5. naive cov (for standard errors) (see summary.lm in base R)
      vcov <- NULL
      if (iv.vcov.stage1 %in% c("lm.vcov.dfres", "lm.vcov")) {
        p1 <- 1L:fit_y_on_xhat$rank
        R <- chol2inv(fit_y_on_xhat$qr$qr[p1, p1, drop = FALSE])
        if (iv.vcov.stage1 == "lm.vcov") {
          vcov <- R * resvar
        } else {
          vcov <- R * resvar_df_res
        }
      }

      # 6. Sargan test (see summary.ivreg.R 363--371)
      if (iv_flag && iv.sargan && iv.vcov.stage1 != "none") {
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
      # eqs[[b]][[j]]$XX <- crossprod(cbind(1, xmat)) # needed? or only vcov?
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
                                          iv.vcov.stage1 = "lm.vcov.dfres",
                                          iv.vcov.stage2 = "delta",
                                          iv.sargan = TRUE,
                                          free.directed.idx = NULL) {
  # no conditional.x for now!
  stopifnot(!lavmodel@conditional.x)
  iv.vcov.stage1 <- tolower(iv.vcov.stage1)
  stopifnot(iv.vcov.stage1 %in%
    c("lm.vcov.dfres", "lm.vcov", "gamma", "none"))

  # get lavpta
  lavpta <- lav_partable_attributes(lavpartable)

  if (is.null(eqs)) {
    # find ALL model-implied instrumental variables (miivs) per equation
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
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
      if (!is.null(eq$iv_type)) {
        iv_flag <- eq$iv_type != "ols"
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
      nx <- length(x.idx)

      # Z: instruments
      S_Xy <- sample.cov[x.idx, y.idx, drop = FALSE]
      S_XX <- sample.cov[x.idx, x.idx, drop = FALSE]
      s_yy <- sample.cov[y.idx, y.idx, drop = FALSE]
      nz <- 0L
      if (iv_flag) {
        i.idx <- match(eq$iv, ov.names)
        nz <- length(i.idx)
        S_XZ <- sample.cov[x.idx, i.idx, drop = FALSE]
        S_ZX <- sample.cov[i.idx, x.idx, drop = FALSE]
        S_ZZ <- sample.cov[i.idx, i.idx, drop = FALSE]
        S_Zy <- sample.cov[i.idx, y.idx, drop = FALSE]
      }

      # sargan vector
      sargan <- rep(as.numeric(NA), 3L)
      names(sargan) <- c("stat", "df", "pvalue")

      # 0. check
      if (iv_flag && length(eq$iv) < length(eq$rhs)) {
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

      # 1. + 2. OLS + 2SLS
      fit_y_on_xhat <- list()
      if (iv_flag) {
        # Step 1: compute S_ZZ^{-1} S_ZX and S_ZZ^{-1} S_Zy
        W <- lav_matrix_symmetric_solve_spd(S_ZZ, S_ZX)
        v <- lav_matrix_symmetric_solve_spd(S_ZZ, S_Zy)
        # Step 2: build reduced system
        Amat <- S_XZ %*% W
        bvec <- S_XZ %*% v
        beta_slopes <- drop(lav_matrix_symmetric_solve_spd(Amat, bvec))
        beta0 <- as.vector(y.bar - t(x.bar) %*% beta_slopes)
      } else {
        if (nx > 0L) {
          Amat <- S_XX
          bvec <- S_Xy
          beta_slopes <- drop(lav_matrix_symmetric_solve_spd(Amat, bvec))
          beta0 <- as.vector(y.bar - t(x.bar) %*% beta_slopes)
        } else {
          beta_slopes <- numeric(0L)
          beta0 <- y.bar
        }
      }
      fit_y_on_xhat$coefficients <- c(beta0, beta_slopes)

      # 2b. jac_eqs_k
      k_mat <- NULL
      k_mat_int <- NULL
      if (iv.vcov.stage1 == "gamma" || iv.vcov.stage2 == "h2") {
        nvar <- nrow(sample.cov)
        pstar <- nvar * (nvar + 1) / 2
        # Ex[k, j] = 1 iff x.idx[k] == j
        Ex <- matrix(0.0, nx, nvar)
        if (nx > 0L) {
          Ex[cbind(seq_len(nx), x.idx)] <- 1.0
        }
        # ---- Jacobian of beta_slopes w.r.t. vech(S) --------------
        a_vec <- lav_matrix_vech_row_idx(nvar)
        b_vec <- lav_matrix_vech_col_idx(nvar)
        offdiag <- (a_vec > b_vec) # logical: strictly off-diagonal elements

        if (!iv_flag && nx == 0L) {
          J_slopes_S <- matrix(0.0, nrow = 0L, ncol = pstar)
        } else if (!iv_flag) {
          g <- -as.vector(crossprod(Ex, beta_slopes))
          g[y.idx] <- g[y.idx] + 1.0
          # build M (p x m) column-by-column in vectorised form
          M <- sweep(Ex[, b_vec, drop = FALSE], 2L, g[a_vec], `*`) +
            sweep(Ex[, a_vec, drop = FALSE], 2L, g[b_vec] * offdiag, `*`)
          J_slopes_S <- lav_matrix_symmetric_solve_spd(S_XX, M)
        } else {
          r_Z <- drop(S_Zy) - S_ZX %*% beta_slopes
          u <- drop(lav_matrix_symmetric_solve_spd(S_ZZ, r_Z))
          Ez <- matrix(0.0, nz, nvar)
          Ez[cbind(seq_len(nz), i.idx)] <- 1.0
          Ew <- crossprod(W, Ez)
          cu <- as.vector(crossprod(Ez, u))
          cb_x <- as.vector(crossprod(Ex, beta_slopes))
          h <- -cu
          h[y.idx] <- h[y.idx] + 1.0
          # build M (p x m) vectorised
          M <-
            sweep(Ex[, b_vec, drop = FALSE], 2L, cu[a_vec], `*`) +
            sweep(Ew[, a_vec, drop = FALSE], 2L, h[b_vec], `*`) -
            sweep(Ew[, b_vec, drop = FALSE], 2L, cb_x[a_vec], `*`) +
            sweep(Ex[, a_vec, drop = FALSE], 2L, cu[b_vec] * offdiag, `*`) +
            sweep(Ew[, b_vec, drop = FALSE], 2L, h[a_vec] * offdiag, `*`) -
            sweep(Ew[, a_vec, drop = FALSE], 2L, cb_x[b_vec] * offdiag, `*`)
          J_slopes_S <- solve(Amat, M)
        }

        # fill in k_mat_int and k_mat
        k_mat_int <- numeric(0L)
        if (!lavmodel@categorical && lavmodel@meanstructure) {
          k_mat_int <- numeric(nvar + pstar)
          k_mat_int[y.idx] <- 1.0
          if (nx > 0L) {
            k_mat_int[x.idx] <- -beta_slopes
            k_mat_int[nvar + seq_len(pstar)] <- -x.bar %*% J_slopes_S
          }
          k_mat <- matrix(0.0, nx, nvar + pstar)
          if (nx > 0L) {
            k_mat[, nvar + seq_len(pstar)] <- J_slopes_S
          }
        } else if (lavmodel@categorical) {
          # split: var_num, off
          # FIXME: all ordinal for now (ie no variances)
          tmp <- J_slopes_S[, -lav_matrix_diagh_idx(nvar), drop = FALSE]
          nth <- nrow(lavsamplestats@NACOV[[b]]) - ncol(tmp)
          k_mat_int <- numeric(nth + ncol(tmp))
          k_mat <- matrix(0.0, nx, nth + ncol(tmp))
          if (nx > 0L) {
            k_mat[, nth + seq_len(ncol(tmp))] <- tmp
          }
        } else {
          k_mat <- matrix(0.0, nrow = nx, ncol = pstar)
          if (nx > 0L) {
            k_mat <- J_slopes_S
          }
        }
        k_mat
      }

      # 3. fill estimates in x
      # - eq$pt contains partable/user rows
      # - lavpartable$free[eq$pt] should give the free idx
      free.idx <- lavpartable$free[eq$pt]
      if (nx > 0L) {
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
      resvar_df_res <- resvar <- df_res <- NULL
      if (iv.vcov.stage1 %in% c("lm.vcov.dfres", "lm.vcov") || iv.sargan) {
        if (nx > 0L) {
          tmp <- S_XX %*% beta_slopes
          resvar <- as.numeric(s_yy - 2 * t(beta_slopes) %*% S_Xy +
            t(beta_slopes) %*% tmp)
        } else {
          resvar <- drop(s_yy)
        }
        df_res <- nobs - (nx + 1L)
        resvar_df_res <- resvar * nobs / df_res
      }

      # 5. naive cov (for standard errors) (see summary.lm in base R)
      vcov <- NULL
      if (iv.vcov.stage1 %in% c("lm.vcov.dfres", "lm.vcov")) {
        this_resvar <- resvar_df_res
        if (iv.vcov.stage1 == "lm.vcov") {
          this_resvar <- resvar
        }
        if (nx > 0L) {
          Ainv <-
            lav_matrix_symmetric_solve_spd(Amat, diag(x = 1, nrow = nrow(Amat)))
          vcov.slopes <- (this_resvar / nobs) * Ainv
          vcov.beta0 <-
            (this_resvar / nobs) + t(x.bar) %*% vcov.slopes %*% x.bar
          cov.beta0.slopes <- -vcov.slopes %*% x.bar
          vcov <- matrix(0, nx + 1L, nx + 1L)
          vcov[1, 1] <- vcov.beta0
          vcov[1, -1] <- t(cov.beta0.slopes)
          vcov[-1, 1] <- cov.beta0.slopes
          vcov[-1, -1] <- vcov.slopes
        } else {
          # only intercept
          vcov <- matrix(this_resvar / nobs, 1L, 1L)
        }
      }

      # 6. Sargan test (see summary.ivreg.R 363--371)
      if (iv_flag && iv.sargan && iv.vcov.stage1 != "none") {
        sargan["df"] <- nz - nx
        if (sargan["df"] > 0L) {
          g <- S_Zy - S_ZX %*% beta_slopes
          w <- lav_matrix_symmetric_solve_spd(S_ZZ, g)
          sargan["stat"] <- as.numeric(nobs * t(g) %*% w / resvar_df_res)
          if (!lavmodel@categorical) {
            sargan["pvalue"] <- pchisq(sargan["stat"], sargan["df"],
              lower.tail = FALSE
            )
          }
        }
      }

      # XX (needed?)
      # XX <- NULL
      #       if (!lavmodel@categorical) {
      #         if (nx > 0L) {
      #           XX <- tcrossprod(c(1, x.bar)) * nobs
      #           XX[-1, -1] <- (S_XX + tcrossprod(x.bar)) * nobs
      #         } else {
      #           XX <- matrix(nobs, 1L, 1L)
      #         }
      #       }

      # add info to eqs list
      eqs[[b]][[j]]$coef <- unname(fit_y_on_xhat$coefficients)
      eqs[[b]][[j]]$nobs <- nobs
      eqs[[b]][[j]]$df_res <- df_res
      eqs[[b]][[j]]$resvar <- resvar
      eqs[[b]][[j]]$resvar_df_res <- resvar_df_res
      # eqs[[b]][[j]]$XX <- XX
      eqs[[b]][[j]]$vcov <- vcov
      eqs[[b]][[j]]$sargan <- sargan
      eqs[[b]][[j]]$k_mat_int <- k_mat_int
      eqs[[b]][[j]]$k_mat <- k_mat
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


# second stage: compute theta_2 as follows:
# theta_2 = solve(t(Delta2) %*% W2 %*% Delta2) %*% t(Delta2) %*% W2 %*% svec
#
# this is the fast but single-block, no-equality-constraints version
# (for now)
#
# see lav_sem_miiv_varcov_old() for an alternative
lav_sem_miiv_varcov <- function(x = NULL, samplestats = FALSE,
                                lavmodel = NULL, lavpartable = NULL,
                                lavsamplestats = NULL,
                                lavh1 = NULL, free.directed.idx = NULL,
                                free.undirected.idx = NULL,
                                add.h2 = FALSE,
                                iv.varcov.method = "RLS") {

  # x = NULL?
  if (is.null(x)) {
    # calling from vcov?
    x <- lav_model_get_parameters(lavmodel)
    lavmodel.tmp <- lavmodel
    implied <- lavh1$implied
  } else if (samplestats) {
    # only used to compute jac_b numerically
    implied <- lav_vec_to_implied(x, lavmodel = lavmodel)
    x <- lav_model_get_parameters(lavmodel)
    lavmodel.tmp <- lavmodel
  } else {
    # first time, or to compute jac_a numerically
    theta1 <- x
    x <- lav_model_get_parameters(lavmodel)
    x[free.directed.idx] <- theta1
    lavmodel.tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
    implied <- lavh1$implied
  }

  # number of blocks
  nblocks <- b <- 1L # for now

  # delta2
  tmp <- lav_model_delta_lisrel(lavmodel.tmp, block = b)
  delta2 <- tmp[, free.undirected.idx, drop = FALSE]

  # svec
  sample_cov <- implied$cov[[b]]
  sample_mean <- implied$mean[[b]]
  if (lavmodel@categorical) {
    sample_th <- implied$th[[b]]
    svec <- c(sample_th, lav_matrix_vech(sample_cov, diagonal = FALSE))
  } else if (!lavmodel@categorical && lavmodel@meanstructure) {
    svec <- c(sample_mean, lav_matrix_vech(sample_cov))
  } else {
    svec <- lav_matrix_vech(sample_cov)
  }

  # initial estimate for theta.2
  if (iv.varcov.method == "GLS") {
    S <- sample_cov
  } else { # ULS, or starting point for 2RLS/RLS
    S <- diag(1, nrow = nrow(sample_cov))
  }
  theta2 <- lav_utils_wls_linearization(Delta = delta2, S = S,
                                        meanstructure = lavmodel@meanstructure,
                                        categorical = lavmodel@categorical,
                                        svec = svec, return_H = add.h2)

  # if ULS/GLS, we are done

  # 2RLS
  if (iv.varcov.method == "2RLS") {
    x[free.undirected.idx] <- theta2
    lavmodel.tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
    sigma <- lav_model_sigma(lavmodel = lavmodel.tmp, extra = FALSE)[[b]]
    theta2 <-
      lav_utils_wls_linearization(Delta = delta2, S = sigma,
                                  meanstructure = lavmodel@meanstructure,
                                  svec = svec, return_H = add.h2)
  # RLS
  } else if (iv.varcov.method == "RLS") {
    # typically, we need 5-15 iterations, so max = 200 should be enough
    for (i in seq_len(200L)) {
      old_x <- theta2
      x[free.undirected.idx] <- theta2
      lavmodel.tmp <- lav_model_set_parameters(lavmodel = lavmodel, x = x)
      sigma <- lav_model_sigma(lavmodel = lavmodel.tmp, extra = FALSE)[[b]]
      theta2 <-
        lav_utils_wls_linearization(Delta = delta2, S = sigma,
                                    meanstructure = lavmodel@meanstructure,
                                    svec = svec, return_H = add.h2)
      sse <- Re(sum((old_x - theta2)^2))
      if (sse < 1e-12 * Re(1 + sum(theta2^2))) {
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
  iv.vcov.stage1 <- tolower(lavoptions$estimator.args$iv.vcov.stage1)
  iv.vcov.stage2 <- tolower(lavoptions$estimator.args$iv.vcov.stage2)
  iv.samplestats <- lavoptions$estimator.args$iv.samplestats
  iv.vcov.gamma.modelbased <-
    lavoptions$estimator.args$iv.vcov.gamma.modelbased
  iv.varcov.method <- toupper(lavoptions$estimator.args$iv.varcov.method)
  iv.vcov.jack.numerical <- lavoptions$estimator.args$iv.vcov.jack.numerical
  iv.vcov.jaca.numerical <- lavoptions$estimator.args$iv.vcov.jaca.numerical
  iv.vcov.jacb.numerical <- lavoptions$estimator.args$iv.vcov.jacb.numerical

  # empty vcov
  vcov <- matrix(0, lavmodel@nx.free, lavmodel@nx.free)

  # nblocks
  nblocks <- lavmodel@nblocks
  stopifnot(lavmodel@nblocks == 1L) # for now...

  # eqs?
  if (is.null(eqs)) {
    eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpartable = lavpartable)
  }

  # directed versus undirected (free) parameters
  undirected.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "~~")
  directed.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~", "~1")) # WITH intercepts!!
  directed_noint.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op %in% c("=~", "~")) # NO intercepts!!
  if (lavmodel@categorical) {
    th.idx <- which(lavpartable$free > 0L &
    !duplicated(lavpartable$free) & # if ceq.simple
    lavpartable$op == "|")
    if (length(unlist(lavmodel@num.idx)) > 0L) {
      lav_msg_warn(gettext(
        "[IV] no standard errors (yet) when variables are both ordinal
         and continuous"
      ))
      iv.vcov.stage1 <- "none"
      iv.vcov.stage2 <- "none"
    }
  }

  free.directed.idx <- unique(lavpartable$free[directed.idx])
  free.directed_noint.idx <- unique(lavpartable$free[directed_noint.idx])
  free.undirected.idx <- unique(lavpartable$free[undirected.idx])
  if (lavmodel@categorical) {
    free.th.idx <- unique(lavpartable$free[th.idx])
  }

  # directed free parameters only
  x <- lav_model_get_parameters(lavmodel, type = "free")
  theta1 <- x[free.directed.idx]
  theta1_noint <- x[free.directed_noint.idx]
  theta2 <- x[free.undirected.idx] # final estimate

  # switch off iv.vcov.stage2?
  if (length(free.directed.idx) > 0L && iv.vcov.stage1 == "none") {
    iv.vcov.stage2 <- "none"
  }

  # do we need gamma_big?
  if (iv.vcov.stage1 == "gamma" ||
     (iv.vcov.stage2 != "none" && length(free.undirected.idx) > 0L)) {
      gamma_g <- vector("list", lavmodel@ngroups)
    for (g in seq_len(lavmodel@ngroups)) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      if (!is.null(lavsamplestats@NACOV[[g]])) {
        gamma_g[[g]] <- fg * lavsamplestats@NACOV[[g]]
      } else {
        if (iv.vcov.gamma.modelbased) {
          cov_g <- lavimplied$cov[[g]]
        } else {
          cov_g <- lavh1$implied$cov[[g]]
        }
#         # NT version (for now), model-based
#         gamma_g[[g]] <- lav_samplestats_Gamma_NT(
#           COV = cov_g,
#           MEAN = mean_g,
#           x.idx = lavsamplestats@x.idx[[g]],
#           fixed.x = lavmodel@fixed.x,
#           conditional.x = lavmodel@conditional.x,
#           meanstructure = lavmodel@meanstructure,
#           slopestructure = lavmodel@conditional.x
#         )
#         gamma_g[[g]] <- fg * gamma_g[[g]]
      }
    }
    if (!is.null(lavsamplestats@NACOV[[g]])) {
      gamma_big <- lav_matrix_bdiag(gamma_g)
    }
  }

  # compute jac_k if needed
  if (length(free.directed.idx) > 0L && iv.vcov.stage1 != "none"
      && (iv.vcov.stage1 == "gamma" || iv.vcov.stage2 == "h2")) {
    if (iv.vcov.jack.numerical) {
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
    } else {
      # analytic version -- # FIXME: multiple blocks!!!
      stopifnot(lavmodel@nblocks == 1L)
      jac_k <- lav_sem_miiv_utils_jack_eqs(
        eqs = eqs, block = 1L,
        lavmodel = lavmodel, lavpartable = lavpartable,
        free.directed.idx = free.directed.idx
      )
    }
  }

  # stage 1: directed effects
  if (length(free.directed.idx) > 0L && iv.vcov.stage1 != "none") {
    if (iv.vcov.stage1 == "gamma") {
      if (lavmodel@categorical) {
        k_gamma_adf_kt <- jac_k %*% gamma_big %*% t(jac_k)
        vcov[free.directed.idx, free.directed.idx] <-
          k_gamma_adf_kt / lavsamplestats@ntotal
      } else {
        # vcov = K %*% Gamma_NT %*% t(K)
        # k_gammant_kt <- jac_k %*% gamma_big %*% t(jac_k)
        x.idx <- if (lavmodel@fixed.x) lavsamplestats@x.idx[[1]] else integer(0L)
        k_gammant_kt <- lav_matrix_k_gammant_kt(K = jac_k, S = cov_g,
          meanstructure = lavmodel@meanstructure, x.idx = x.idx)
        vcov[free.directed.idx, free.directed.idx] <-
          k_gammant_kt / lavsamplestats@ntotal
      }

    # iv.vcov.stage1 = lm.vcov or lm.vcov.dfres
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
  } # stage 1

  # stage 2: undirected effects (note: intercepts have no effect!)

  if (iv.vcov.stage2 != "none" && length(free.undirected.idx) > 0L &&
    iv.varcov.method != "NONE") {
    if (iv.vcov.stage2 == "h2") {
      # delta_list <- lav_model_delta(lavmodel) # per block/group

      # # block-wise
      # delta1_block <- vector("list", length = nblocks)
      # delta2_block <- vector("list", length = nblocks)
      # w2_block <- vector("list", length = nblocks)
      # for (b in seq_len(nblocks)) {
      #   nvar <- lavmodel@nvar[b]
      #   # delta1 + delta2
      #   delta1_block[[b]] <- delta_list[[b]][, free.directed.idx, drop = FALSE]
      #   delta2_block[[b]] <-
      #     delta_list[[b]][, free.undirected.idx, drop = FALSE]
      #   # w2
      #   if (lavmodel@categorical || iv.varcov.method == "ULS") {
      #     s.inv <- diag(1, nrow = nvar)
      #   } else if (iv.varcov.method == "GLS") {
      #     s.inv <- lavh1$implied$cov[[b]]
      #   } else {
      #     s.inv <- lavimplied$cov[[b]]
      #   }
      #   w2_22 <- 0.5 * lav_matrix_duplication_pre_post(s.inv %x% s.inv)
      #   if (lavmodel@meanstructure) {
      #     w2_block[[b]] <- lav_matrix_bdiag(s.inv, w2_22)
      #   } else {
      #     w2_block[[b]] <- w2_22
      #   }
      # } # blocks
      # delta1 <- do.call("rbind", delta1_block)
      # delta2 <- do.call("rbind", delta2_block)
      # w2 <- lav_matrix_bdiag(w2_block)
      # h2 <- solve(t(delta2) %*% w2 %*% delta2, t(delta2) %*% w2)

      # tmp <- lav_sem_miiv_varcov(
      #   x = NULL, lavmodel = lavmodel,
      #   lavpartable = lavpartable, lavsamplestats = lavsamplestats,
      #   lavh1 = lavh1, free.directed.idx = free.directed.idx,
      #   free.undirected.idx = free.undirected.idx, add.h2 = TRUE,
      #   iv.varcov.method = iv.varcov.method
      # )
      # h2 <- attr(tmp, "H")
      stopifnot(lavmodel@nblocks == 1L)
      delta <- lav_model_delta_lisrel(lavmodel, block = 1L)
      delta2 <- delta[, free.undirected.idx, drop = FALSE]
      svec <- lav_implied_to_vec(
          implied = lavh1$implied, lavmodel = lavmodel,
          drop.list = TRUE
      )
      if (lavmodel@categorical || iv.varcov.method == "ULS") {
        s_mat <- diag(lavmodel@nvar[[1]])
      } else if (iv.varcov.method == "GLS") {
        s_mat <- lavh1$implied$cov[[1]]
      } else {
        s_mat <- lavimplied$cov[[1]]
      }
      tmp <- lav_utils_wls_linearization(Delta = delta2, S = s_mat,
        meanstructure = lavmodel@meanstructure,
        categorical = lavmodel@categorical, svec = svec, return_H = TRUE)
      h2 <- attr(tmp, "H")
      if (length(free.directed.idx) > 0L) {
        # assuming a SINGLE block for now
        delta1 <- delta[, free.directed.idx, drop = FALSE]
        idiag <- diag(1, nrow = nrow(delta1))
        nth <- nrow(delta1) - ncol(h2)
        h2 <- cbind(matrix(0, nrow(h2), ncol = nth), h2)
        tmp <- h2 %*% (idiag - delta1 %*% jac_k)
      } else {
        tmp <- h2
      }
      if (lavmodel@categorical) {
        tmp_gamma_tmpt <- tmp %*% gamma_big %*% t(tmp)
        vcov[free.undirected.idx, free.undirected.idx] <-
          tmp_gamma_tmpt / lavsamplestats@ntotal
        # add thresholds
        nth <- length(free.th.idx)
        if (length(nth) > 0L) {
          delta.th <- delta[, free.th.idx, drop = FALSE]
          # FIXME: remove zero rows from delta
          zero_idx <- which(rowSums(delta.th) == 0)
          if (length(zero_idx) > 0L) {
            delta.th <- delta.th[-zero_idx, , drop = FALSE]
            gamma_th <- gamma_big[-zero_idx, -zero_idx, drop = FALSE]
          }
          th_gamma_tht <- t(delta.th) %*% gamma_th %*% delta.th
          vcov[free.th.idx, free.th.idx] <-
            th_gamma_tht / lavsamplestats@ntotal
        }
      } else {
        x.idx <-
          if (lavmodel@fixed.x) lavsamplestats@x.idx[[1]] else integer(0L)
        tmp_gammant_tmpt <- lav_matrix_k_gammant_kt(K = tmp, S = cov_g,
          meanstructure = lavmodel@meanstructure, x.idx = x.idx)
        #tmp_gammant_tmpt_bis <- tmp %*% gamma_big %*% t(tmp)
        vcov[free.undirected.idx, free.undirected.idx] <-
          tmp_gammant_tmpt / lavsamplestats@ntotal
      }

      # iv.vcov.stage2 == "delta"
    } else {
      # part a: effect of theta1
      if (!iv.vcov.jaca.numerical) {
        # use analytic expressions
        if (iv.varcov.method %in% c("ULS", "GLS")) {
          jac_a <- lav_sem_miiv_utils_jaca_uls_gls(
            lavmodel = lavmodel,
            lavpartable = lavpartable, lavh1 = lavh1,
            free.directed.idx = free.directed_noint.idx,
            free.undirected.idx = free.undirected.idx,
            iv.varcov.method = iv.varcov.method
          )
        } else if (iv.varcov.method == "2RLS") {
          jac_a <- lav_sem_miiv_utils_jaca_2rls(
            lavmodel = lavmodel,
            lavpartable = lavpartable, lavh1 = lavh1,
            free.directed.idx = free.directed_noint.idx,
            free.undirected.idx = free.undirected.idx
          )
        } else if (iv.varcov.method == "RLS") {
          jac_a <- lav_sem_miiv_utils_jaca_rls(
            lavmodel = lavmodel,
            lavpartable = lavpartable, lavh1 = lavh1,
            free.directed.idx = free.directed_noint.idx,
            free.undirected.idx = free.undirected.idx
          )
        }
        # for checking only: numerical approximation
      } else {
        jac_a <- lav_func_jacobian_complex(
          func = lav_sem_miiv_varcov,
          x = theta1_noint, samplestats = FALSE, lavmodel = lavmodel,
          lavpartable = lavpartable, lavsamplestats = lavsamplestats,
          lavh1 = lavh1, free.directed.idx = free.directed_noint.idx,
          free.undirected.idx = free.undirected.idx,
          iv.varcov.method = iv.varcov.method
        )

        # jac_a <- lav_func_jacobian_complex(
        #   func = lav_sem_miiv_varcov_old,
        #   x = theta1_noint, samplestats = FALSE, lavmodel = lavmodel,
        #   lavpartable = lavpartable, lavsamplestats = lavsamplestats,
        #   lavh1 = lavh1, free.directed.idx = free.directed_noint.idx,
        #   free.undirected.idx = free.undirected.idx,
        #   iv.varcov.method = iv.varcov.method
        # )
      }

      # compute vcov_a (due to theta1)
      vcov_directed <-
        vcov[free.directed_noint.idx, free.directed_noint.idx, drop = FALSE]
      vcov_a <- jac_a %*% vcov_directed %*% t(jac_a)

      # part b: effect of sample statistics
      # - get jacobian lav_sem_miiv_varcov wrt samplestats_vec
      # - get Gamma (NT or ADF)
      # - vcov_undirected_b <- jac_b %*% (1/N * gamma_mat) %*% t(jac_b)

      if (!iv.vcov.jaca.numerical) {
        delta_block <- lav_model_delta(
          lavmodel = lavmodel,
          ceq.simple = lavmodel@ceq.simple.only
        )
        jac_b_block <- vector("list", length = nblocks)
        for (b in seq_len(nblocks)) {
          sample_cov <- lavh1$implied$cov[[b]]
          this_delta2 <- delta_block[[b]][, free.undirected.idx, drop = FALSE]
          if (iv.varcov.method == "ULS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_uls(sample_cov, delta2 = this_delta2)
          } else if (iv.varcov.method == "GLS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_gls(sample_cov, delta2 = this_delta2,
                                          theta2 = theta2)
          } else if (iv.varcov.method == "2RLS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_2rls(sample_cov, delta2 = this_delta2,
                                           theta2 = theta2)
          } else if (iv.varcov.method == "RLS") {
            jac_b_block[[b]] <-
              lav_sem_miiv_utils_jacb_rls(sample_cov, delta2 = this_delta2,
                                          theta2 = theta2)
          }
        } # blocks
        jac_b <- lav_matrix_bdiag(jac_b_block)

      # for checking only: numerical approximation (very slow!)
      } else {
        vec <- lav_implied_to_vec(
          implied = lavh1$implied, lavmodel = lavmodel,
          drop.list = TRUE
        )

        jac_b <- lav_func_jacobian_complex(
          func = lav_sem_miiv_varcov,
          x = vec, samplestats = TRUE, lavmodel = lavmodel,
          lavpartable = lavpartable, lavsamplestats = lavsamplestats,
          lavh1 = lavh1, free.directed.idx = free.directed_noint.idx,
          free.undirected.idx = free.undirected.idx,
          iv.varcov.method = iv.varcov.method
        )

        # jac_b <- lav_func_jacobian_complex(
        #   func = lav_sem_miiv_varcov_old,
        #   x = vec, samplestats = TRUE, lavmodel = lavmodel,
        #   lavpartable = lavpartable, lavsamplestats = lavsamplestats,
        #   lavh1 = lavh1, free.directed.idx = free.directed_noint.idx,
        #   free.undirected.idx = free.undirected.idx,
        #   iv.varcov.method = iv.varcov.method
        # )
      }
      x.idx <- if (lavmodel@fixed.x) lavsamplestats@x.idx[[1]] else integer(0L)
      jacb_gammant_jacbt <- lav_matrix_k_gammant_kt(K = jac_b, S = cov_g,
        meanstructure = lavmodel@meanstructure, x.idx = x.idx)
      #jacb_gammant_jacbt_bis <- jac_b %*% gamma_big %*% t(jac_b)
      vcov_b <- jacb_gammant_jacbt / lavsamplestats@ntotal
      vcov_ab <- vcov_a + vcov_b
      vcov[free.undirected.idx, free.undirected.idx] <- vcov_ab
    } # continuous
  } # iv.vcov.stage2 != "none"

  vcov
}
