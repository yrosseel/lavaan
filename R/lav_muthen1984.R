# This function was written in January 2012 -- Yves Rosseel
# First success: Friday 20 Jan 2012: the standard errors for
#                thresholds and polychoric correlations (in an
#                unrestricted/saturated model) are spot on!
# Second success: Saturday 9 June 2012: support for mixed (ordinal + metric)
#                 variables; thanks to the delta method to get the ACOV
#                 right (see H matrix)
# Third success: Monday 2 July 2012: support for fixed.x covariates
#
# Friday 13 July 2012: merge exo + non-exo code
# Monday 16 July 2012: fixed sign numeric in WLS.W; I think we got it right now

# YR 26 Nov 2015: move step1 + step2 to external functions
#
muthen1984 <- function(data_1 = NULL,
                       ov_names = NULL,
                       ov_types = NULL,
                       ov_levels = NULL,
                       ov_names_x = character(0L),
                       exo = NULL,
                       wt = NULL,
                       sampling_weights_type = "design",
                       wls_w = TRUE,
                       zero_add = c(0.5, 0.0),
                       zero_keep_margins = TRUE,
                       zero_cell_warn = FALSE,
                       zero_cell_tables = TRUE,
                       allow_empty_cell = TRUE,
                       group = 1L) { # group only for error messages

  # just in case Data is a vector
  data_1 <- as.matrix(data_1)

  nvar <- NCOL(data_1)
  n <- NROW(data_1)
  num_idx <- which(ov_types == "numeric")
  ord_idx <- which(ov_types == "ordered")
  nexo <- length(ov_names_x)
  if (nexo > 0L) stopifnot(NCOL(exo) == nexo)
  pstar <- nvar * (nvar - 1) / 2

  if (lav_verbose()) {
    cat("\nPreparing for WLS estimation -- STEP 1 + 2\n")
    cat("Number of endogenous variables: ", nvar, "\n")
    cat("Endogenous variable names:\n")
    print(ov_names)
    cat("\n")
    cat("Endogenous ov types:\n")
    print(ov_types)
    cat("\n")
    cat("Endogenous ov levels:\n ")
    print(ov_levels)
    cat("\n")
    cat("Number of exogenous variables: ", nexo, "\n")
    cat("Exogenous variable names:\n")
    print(ov_names_x)
    cat("\n")
  }

  step1 <- lav_samp_step1(
    y = data_1, wt = wt, ov_names = ov_names,
    ov_types = ov_types, ov_levels = ov_levels, ov_names_x = ov_names_x,
    exo = exo, scores_flag = wls_w, allow_empty_cell = allow_empty_cell,
    group = group
  )

  fit <- step1$FIT
  th <- step1$TH
  th_nox <- step1$TH.NOX
  th_idx_1 <- step1$TH.IDX
  th_names <- step1$TH.NAMES
  var_1 <- step1$VAR
  slopes <- step1$SLOPES
  sc_th <- step1$SC.TH
  sc_sl <- step1$SC.SL
  sc_var <- step1$SC.VAR
  th_start_idx <- step1$th.start.idx
  th_end_idx <- step1$th.end.idx

  # rm SC.VAR columns from ordinal variables
  if (wls_w && length(ord_idx) > 0L) {
    sc_var <- sc_var[, -ord_idx, drop = FALSE]
  }


  if (lav_verbose()) {
    cat("STEP 1: univariate statistics\n")
    cat("Threshold + means:\n")
    tthh <- unlist(th)
    names(tthh) <- unlist(th_names)
    print(tthh)
    cat("Slopes (if any):\n")
    colnames(slopes) <- ov_names_x
    rownames(slopes) <- ov_names
    print(slopes)
    cat("Variances:\n")
    names(var_1) <- ov_names
    print(unlist(var_1))
  }

  # stage two -- correlations

  if (lav_verbose()) cat("\n\nSTEP 2: covariances/correlations:\n")
  cor_1 <- lav_samp_step2(
    uni = fit, wt = wt, ov_names = ov_names,
    zero_add = zero_add,
    zero_keep_margins = zero_keep_margins,
    zero_cell_warn = zero_cell_warn,
    zero_cell_tables = zero_cell_tables
  )
  empty_cell_tables <- attr(cor_1, "zero.cell.tables")
  attr(cor_1, "zero.cell.tables") <- NULL

  if (lav_verbose()) {
    colnames(cor_1) <- rownames(cor_1) <- ov_names
    print(cor_1)
  }

  if (!wls_w) { # we do not need the asymptotic variance matrix
    if (any("numeric" %in% ov_types)) {
      cov_1 <- lav_cor2cov(r = cor_1, sds = sqrt(unlist(var_1)))
    } else {
      cov_1 <- cor_1
    }
    out <- list(
      TH = th, SLOPES = slopes, VAR = var_1, COR = cor_1, COV = cov_1,
      SC = NULL, TH.NOX = th_nox, TH.NAMES = th_names, TH.IDX = th_idx_1,
      INNER = NULL, A11 = NULL, A12 = NULL, A21 = NULL, A22 = NULL,
      WLS.W = NULL, H = NULL, zero.cell.tables = matrix("", 0, 2)
    )
    return(out)
  }


  # stage three -- WLS.W
  sc_cor <- matrix(0, n, pstar)
  pstar_1 <- matrix(0, nvar, nvar)
  pstar_1[lav_mat_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar

  a11_size <- NCOL(sc_th) + NCOL(sc_sl) + NCOL(sc_var)








  # A21
  a21 <- matrix(0, pstar, a11_size)
  h22 <- diag(pstar) # for the delta rule
  h21 <- matrix(0, pstar, a11_size)
  # for this one, we need new scores: for each F_ij (cor), the
  # scores with respect to the TH, VAR, ...
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      pstar_idx <- pstar_1[i, j]
      th_idx_i <- th_start_idx[i]:th_end_idx[i]
      th_idx_j <- th_start_idx[j]:th_end_idx[j]
      if (nexo > 0L) {
        sl_idx_i <- NCOL(sc_th) + seq(i, by = nvar, length.out = nexo)
        sl_idx_j <- NCOL(sc_th) + seq(j, by = nvar, length.out = nexo)

        var_idx_i <- NCOL(sc_th) + NCOL(sc_sl) + match(i, num_idx)
        var_idx_j <- NCOL(sc_th) + NCOL(sc_sl) + match(j, num_idx)
      } else {
        var_idx_i <- NCOL(sc_th) + match(i, num_idx)
        var_idx_j <- NCOL(sc_th) + match(j, num_idx)
      }
      if (ov_types[i] == "numeric" && ov_types[j] == "numeric") {
        sc_cor_uni <- lav_bvreg_cor_sc(
          rho = cor_1[i, j],
          fit_y1 = fit[[i]],
          fit_y2 = fit[[j]],
          wt = wt
        )

        # RHO
        if (is.null(wt)) {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho
        } else {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho / wt # unweight
        }

        # TH
        a21[pstar_idx, th_idx_i] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_mu_y1
          )
        a21[pstar_idx, th_idx_j] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_mu_y2
          )
        # SL
        if (nexo > 0L) {
          a21[pstar_idx, sl_idx_i] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y1
            )
          a21[pstar_idx, sl_idx_j] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y2
            )
        }
        # VAR
        a21[pstar_idx, var_idx_i] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_var_y1
          )
        a21[pstar_idx, var_idx_j] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_var_y2
          )
        # H21 only needed for VAR
        h21[pstar_idx, var_idx_i] <-
          (sqrt(var_1[j]) * cor_1[i, j]) / (2 * sqrt(var_1[i]))
        h21[pstar_idx, var_idx_j] <-
          (sqrt(var_1[i]) * cor_1[i, j]) / (2 * sqrt(var_1[j]))
        h22[pstar_idx, pstar_idx] <- sqrt(var_1[i]) * sqrt(var_1[j])
      } else if (ov_types[i] == "numeric" && ov_types[j] == "ordered") {
        sc_cor_uni <- lav_bvmix_cor_sc(
          rho = cor_1[i, j],
          fit_y1 = fit[[i]],
          fit_y2 = fit[[j]],
          wt = wt
        )
        # RHO
        if (is.null(wt)) {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho
        } else {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho / wt # unweight
        }

        # TH
        a21[pstar_idx, th_idx_i] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_mu_y1
          )
        a21[pstar_idx, th_idx_j] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_th_y2
          )
        # SL
        if (nexo > 0L) {
          a21[pstar_idx, sl_idx_i] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y1
            )
          a21[pstar_idx, sl_idx_j] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y2
            )
        }
        # VAR
        a21[pstar_idx, var_idx_i] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_var_y1
          )
        # H21 only need for VAR
        h21[pstar_idx, var_idx_i] <- cor_1[i, j] / (2 * sqrt(var_1[i]))
        h22[pstar_idx, pstar_idx] <- sqrt(var_1[i])
      } else if (ov_types[j] == "numeric" && ov_types[i] == "ordered") {
        sc_cor_uni <- lav_bvmix_cor_sc(
          rho = cor_1[i, j],
          fit_y1 = fit[[j]],
          fit_y2 = fit[[i]],
          wt = wt
        )
        # RHO
        if (is.null(wt)) {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho
        } else {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho / wt # unweight
        }

        # TH
        a21[pstar_idx, th_idx_j] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_mu_y1
          )
        a21[pstar_idx, th_idx_i] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_th_y2
          )
        # SL
        if (nexo > 0L) {
          a21[pstar_idx, sl_idx_j] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y1
            )
          a21[pstar_idx, sl_idx_i] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y2
            )
        }
        # VAR
        a21[pstar_idx, var_idx_j] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_var_y1
          )
        # H21 only for VAR
        h21[pstar_idx, var_idx_j] <- cor_1[i, j] / (2 * sqrt(var_1[j]))
        h22[pstar_idx, pstar_idx] <- sqrt(var_1[j])
      } else if (ov_types[i] == "ordered" && ov_types[j] == "ordered") {
        # polychoric correlation
        sc_cor_uni <- lav_bvord_cor_sc(
          rho = cor_1[i, j],
          fit_y1 = fit[[i]],
          fit_y2 = fit[[j]],
          wt = wt
        )
        # RHO
        if (is.null(wt)) {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho
        } else {
          sc_cor[, pstar_idx] <- sc_cor_uni$dx_rho / wt # unweight
        }

        # TH
        a21[pstar_idx, th_idx_i] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_th_y1
          )
        a21[pstar_idx, th_idx_j] <-
          lav_mat_crossprod(
            sc_cor[, pstar_idx],
            sc_cor_uni$dx_th_y2
          )
        # SL
        if (nexo > 0L) {
          a21[pstar_idx, sl_idx_i] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y1
            )
          a21[pstar_idx, sl_idx_j] <-
            lav_mat_crossprod(
              sc_cor[, pstar_idx],
              sc_cor_uni$dx_sl_y2
            )
        }
        # NO VAR
      }
    }
  }
  if (!is.null(wt)) {
    sc_cor <- sc_cor * wt # reweight
  }





  # stage three

  sc <- cbind(sc_th, sc_sl, sc_var, sc_cor)
  inner <- lav_mat_crossprod(sc)

  # A11
  # new approach (2 June 2012): A11 is just a 'sparse' version of
  # (the left upper block of) INNER
  a11 <- matrix(0, a11_size, a11_size)
  if (!is.null(wt)) {
    inner2 <- lav_mat_crossprod(sc / wt, sc)
  } else {
    inner2 <- inner
  }
  for (i in 1:nvar) {
    th_idx <- th_start_idx[i]:th_end_idx[i]
    sl_idx <- integer(0L)
    var_idx <- integer(0L)
    if (nexo > 0L) {
      sl_idx <- NCOL(sc_th) + seq(i, by = nvar, length.out = nexo)
      # sl.end.idx <- (i*nexo); sl.start.idx <- (i-1L)*nexo + 1L
      # sl.idx <- NCOL(SC.TH) + (sl.start.idx:sl.end.idx)
    }
    if (ov_types[i] == "numeric") {
      var_idx <- NCOL(sc_th) + NCOL(sc_sl) + match(i, num_idx)
    }
    a11_idx <- c(th_idx, sl_idx, var_idx)
    a11[a11_idx, a11_idx] <- inner2[a11_idx, a11_idx]
  }

  ##### DEBUG ######
  #### for numeric VAR only, use hessian to get better residual var value
  ####
  # for(i in 1:nvar) {
  #     if(ov.types[i] == "numeric") {
  #         tmp.npar <- FIT[[i]]$npar
  #         e.var <- FIT[[i]]$theta[ tmp.npar ]
  #         sq.e.var <- sqrt(e.var)
  #         sq.e.var6 <- sq.e.var*sq.e.var*sq.e.var*sq.e.var*sq.e.var*sq.e.var
  #         dx2.var <- N/(2*e.var*e.var) - 1/sq.e.var6 * (e.var * N)
  #
  #         var.idx <- NCOL(SC.TH) + NCOL(SC.SL) + match(i, num.idx)
  #         A11[var.idx, var.idx] <- -1 * dx2.var
  #     }
  # }
  ################
  ################

  # A22 (diagonal)
  a22 <- matrix(0, pstar, pstar)
  for (i in seq_len(pstar)) {
    if (is.null(wt)) {
      a22[i, i] <- sum(sc_cor[, i] * sc_cor[, i], na.rm = TRUE)
    } else {
      a22[i, i] <- sum(sc_cor[, i] * sc_cor[, i] / wt, na.rm = TRUE)
    }
  }

  # A12 (zero)
  a12 <- matrix(0, NROW(a11), NCOL(a22))


  # B <- rbind( cbind(A11,A12),
  #            cbind(A21,A22) )

  # we invert B as a block-triangular matrix (0.5-23)
  #
  # B.inv = A11^{-1}                   0
  #         -A22^{-1} A21 A11^{-1}     A22^{-1}
  #

  # invert A
  a11_inv <- try(solve(a11), silent = TRUE)
  if (inherits(a11_inv, "try-error")) {
    # brute force
    a11_inv <- MASS::ginv(a11)
    lav_msg_warn(gettext("trouble constructing W matrix;
                         used generalized inverse for A11 submatrix"))
  }

  # invert
  da22 <- diag(a22)
  if (any(da22 == 0)) {
    lav_msg_warn(gettext("trouble constructing W matrix;
                         used generalized inverse for A22 submatrix"))
    a22_inv <- MASS::ginv(a22)
  } else {
    a22_inv <- a22
    diag(a22_inv) <- 1 / da22
  }

  # lower-left block
  a21_inv <- -a22_inv %*% a21 %*% a11_inv

  # upper-left block remains zero
  a12_inv <- a12

  # construct B.inv
  b_inv <- rbind(
    cbind(a11_inv, a12_inv),
    cbind(a21_inv, a22_inv)
  )


  #  weight matrix (correlation metric)
  # The bread b_inv is built from sum(wt)-weighted quantities (a11/a22 use
  # 'inner2'). The choice of meat encodes how the sampling weights are
  # interpreted (sampling.weights.type):
  #  - "design"    : meat = 'inner' = crossprod(sc) where sc is wt-scaled,
  #                  i.e. weighted by sum(wt^2). b_inv %*% inner %*% b_inv is
  #                  then the design-based sandwich Gamma (matches Mplus).
  #  - "frequency" : meat = 'inner2' = crossprod(sc/wt, sc), weighted by
  #                  sum(wt); treats the weights as frequencies, so Gamma
  #                  (and the WLS weight matrix and robust SEs) matches a fit
  #                  on the row-replicated data.
  # Without sampling weights inner2 == inner, so the choice is a no-op.
  if (!is.null(wt) && identical(sampling_weights_type, "frequency")) {
    wls_w <- b_inv %*% inner2 %*% t(b_inv)
  } else {
    wls_w <- b_inv %*% inner %*% t(b_inv)
  }

  # COV matrix?
  if (any("numeric" %in% ov_types)) {
    cov_1 <- lav_cor2cov(r = cor_1, sds = sqrt(unlist(var_1)))

    # construct H matrix to apply delta rule (for the transformation
    # of rho_ij to cov_ij)
    h11 <- diag(NROW(a11))
    h12 <- matrix(0, NROW(a11), NCOL(a22))
    # H22 and H21 already filled in
    h <- rbind(
      cbind(h11, h12),
      cbind(h21, h22)
    )

    wls_w <- h %*% wls_w %*% t(h)
  } else {
    cov_1 <- cor_1
    h <- diag(NCOL(wls_w))
  }

  # reverse sign numeric TH (because we provide -mu in WLS.obs)
  # (WOW, it took me a LOOONGGG time to realize this!)
  # YR 16 July 2012

  # NOTE: prior to 0.5-17, we used num.idx (instead of NUM.idx)
  # which is WRONG if we have more than one threshold per variable
  # (thanks to Sacha Epskamp for spotting this!)
  if (length(num_idx) > 0L) {
    num_idx_1 <- which(unlist(th_idx_1) == 0L)
    wls_w[num_idx_1, ] <- -wls_w[num_idx_1, ]
    wls_w[, num_idx_1] <- -wls_w[, num_idx_1]
  }


  out <- list(
    TH = th, SLOPES = slopes, VAR = var_1, COR = cor_1, COV = cov_1,
    SC = sc, TH.NOX = th_nox, TH.NAMES = th_names, TH.IDX = th_idx_1,
    INNER = inner, A11 = a11, A12 = a12, A21 = a21, A22 = a22,
    WLS.W = wls_w, H = h,
    zero.cell.tables = empty_cell_tables
  )
  out
}
