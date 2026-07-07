# the Muthen (1984) three-stage estimator for mixed continuous/ordinal data:
# univariate parameters (step 1), correlations (step 2), and the asymptotic
# covariance matrix WLS.W of all estimates (step 3, Muthen & Satorra 1995)
#
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
# YR/LDW 2026: refactored; step 3 decomposed into named helpers
#              (lav_m84_*); no change in behavior
#
# Step 3 computes the sandwich
#
#     WLS.W = B^{-1} INNER B^{-T}
#
# where the 'bread' B is the (block-triangular) matrix of derivatives of the
# stacked stage-wise estimating equations with respect to the parameters
#
#     B = | A11   0  |     A11: univariate (th/sl/var) information
#         | A21  A22 |     A21: d(cor scores)/d(univariate parameters)
#                          A22: diagonal correlation information
#
# and INNER is the crossproduct of the stacked casewise scores
# (with variants for sampling weights and cluster-robust aggregation).
# The H matrix afterwards converts the correlation metric to the
# covariance metric (delta rule) when numeric variables are present.

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
                       cluster_idx = NULL,
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
  th_keep <- step1$TH.KEEP
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

  # index layout of the stacked parameter vector
  # (th/means | slopes (if any) | variances (numeric only) | correlations)
  layout <- lav_m84_layout(
    nvar = nvar, num_idx = num_idx,
    th_start_idx = th_start_idx, th_end_idx = th_end_idx,
    nth = NCOL(sc_th), nsl = NCOL(sc_sl), nexo = nexo
  )
  a11_size <- layout$a11_size

  # A21 (+ the H21/H22 blocks of the delta-rule matrix, and the
  #      correlation scores)
  tmp <- lav_m84_a21(
    fit = fit, ov_types = ov_types, cor_1 = cor_1, var_1 = var_1,
    wt = wt, n = n, layout = layout, th_keep = th_keep
  )
  sc_cor <- tmp$sc_cor
  a21 <- tmp$a21
  h21 <- tmp$h21
  h22 <- tmp$h22

  if (!is.null(wt)) {
    sc_cor <- sc_cor * wt # reweight
  }

  # stacked scores + meat
  sc <- cbind(sc_th, sc_sl, sc_var, sc_cor)
  inner <- lav_mat_crossprod(sc)

  # A11
  # new approach (2 June 2012): A11 is just a 'sparse' version of
  # (the left upper block of) INNER
  if (!is.null(wt)) {
    inner2 <- lav_mat_crossprod(sc / wt, sc)
  } else {
    inner2 <- inner
  }
  a11 <- lav_m84_a11(inner2 = inner2, layout = layout)

  # A22 (diagonal)
  a22 <- lav_m84_a22(sc_cor = sc_cor, wt = wt)

  # A12 (zero)
  a12 <- matrix(0, NROW(a11), NCOL(a22))

  # invert B as a block-triangular matrix (0.5-23)
  tmp <- lav_m84_b_inv(a11 = a11, a21 = a21, a22 = a22)
  b_inv <- tmp$b_inv
  # (warnings are emitted here, not in the helper, so that the message
  #  mentions muthen1984 as before)
  if (tmp$ginv_a11) {
    lav_msg_warn(gettext("trouble constructing W matrix;
                         used generalized inverse for A11 submatrix"))
  }
  if (tmp$ginv_a22) {
    lav_msg_warn(gettext("trouble constructing W matrix;
                         used generalized inverse for A22 submatrix"))
  }

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

  # cluster-robust meat (issue #254): for cluster-robust Gamma/NACOV we keep
  # the same 'bread' (b_inv), but aggregate the case-wise scores per cluster
  # before forming the outer product -- exactly as in the continuous case
  # (lav_samp_gamma: zc <- rowsum(zc, cluster_idx)). Only affects Gamma/NACOV
  # (hence the robust SEs and scaled test); the WLS weight matrix (and thus the
  # point estimates) is still built from the standard i.i.d. 'wls_w' above.
  wls_w_cl <- NULL
  if (!is.null(cluster_idx)) {
    sc_cl <- rowsum.default(sc, group = cluster_idx, reorder = FALSE)
    inner_cl <- lav_mat_crossprod(sc_cl)
    wls_w_cl <- b_inv %*% inner_cl %*% t(b_inv)
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
  if (!is.null(wls_w_cl)) {
    wls_w_cl <- h %*% wls_w_cl %*% t(h)
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
    if (!is.null(wls_w_cl)) {
      wls_w_cl[num_idx_1, ] <- -wls_w_cl[num_idx_1, ]
      wls_w_cl[, num_idx_1] <- -wls_w_cl[, num_idx_1]
    }
  }


  out <- list(
    TH = th, SLOPES = slopes, VAR = var_1, COR = cor_1, COV = cov_1,
    SC = sc, TH.NOX = th_nox, TH.NAMES = th_names, TH.IDX = th_idx_1,
    INNER = inner, A11 = a11, A12 = a12, A21 = a21, A22 = a22,
    WLS.W = wls_w, WLS.W.cluster = wls_w_cl, H = h,
    zero.cell.tables = empty_cell_tables
  )
  out
}

# index layout of the stacked (univariate | correlation) parameter vector:
# for each variable, the positions of its thresholds/mean, slopes and
# (numeric only) variance within the A11 block; and the pstar mapping of
# the correlations
lav_m84_layout <- function(nvar = NULL, num_idx = NULL,
                           th_start_idx = NULL, th_end_idx = NULL,
                           nth = NULL, nsl = NULL, nexo = NULL) {
  pstar <- nvar * (nvar - 1) / 2
  pstar_1 <- matrix(0, nvar, nvar)
  pstar_1[lav_mat_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar

  th_idx <- vector("list", length = nvar)
  sl_idx <- vector("list", length = nvar)
  var_idx <- vector("list", length = nvar)
  for (i in 1:nvar) {
    th_idx[[i]] <- th_start_idx[i]:th_end_idx[i]
    if (nexo > 0L) {
      sl_idx[[i]] <- nth + seq(i, by = nvar, length.out = nexo)
    } else {
      sl_idx[[i]] <- integer(0L)
    }
    # variance slot (numeric variables only; NA for ordinal)
    var_idx[[i]] <- nth + nsl + match(i, num_idx)
  }

  list(
    nvar = nvar, pstar = pstar, pstar_1 = pstar_1, nexo = nexo,
    th_idx = th_idx, sl_idx = sl_idx, var_idx = var_idx,
    a11_size = nth + nsl + length(num_idx)
  )
}

# the casewise scores of one correlation (pair i > j), together with the
# cross-derivatives with respect to the univariate parameters of both
# variables, mapped to a type-independent layout:
#   dx_rho, dx_th_i/dx_th_j (thresholds or mean), dx_sl_i/dx_sl_j (slopes),
#   dx_var_i/dx_var_j (residual variance; NULL for ordinal variables)
lav_m84_pair_sc <- function(fit = NULL, ov_types = NULL,
                            i = NULL, j = NULL, rho = NULL, wt = NULL) {
  ord_i <- ov_types[i] == "ordered"
  ord_j <- ov_types[j] == "ordered"
  if (!ord_i && !ord_j) {
    # pearson
    sc <- lav_bvreg_cor_sc(
      rho = rho, fit_y1 = fit[[i]], fit_y2 = fit[[j]], wt = wt
    )
    out <- list(
      dx_rho = sc$dx_rho,
      dx_th_i = sc$dx_mu_y1, dx_th_j = sc$dx_mu_y2,
      dx_sl_i = sc$dx_sl_y1, dx_sl_j = sc$dx_sl_y2,
      dx_var_i = sc$dx_var_y1, dx_var_j = sc$dx_var_y2
    )
  } else if (!ord_i && ord_j) {
    # polyserial (i = linear, j = ordinal)
    sc <- lav_bvmix_cor_sc(
      rho = rho, fit_y1 = fit[[i]], fit_y2 = fit[[j]], wt = wt
    )
    out <- list(
      dx_rho = sc$dx_rho,
      dx_th_i = sc$dx_mu_y1, dx_th_j = sc$dx_th_y2,
      dx_sl_i = sc$dx_sl_y1, dx_sl_j = sc$dx_sl_y2,
      dx_var_i = sc$dx_var_y1, dx_var_j = NULL
    )
  } else if (ord_i && !ord_j) {
    # polyserial (j = linear, i = ordinal)
    sc <- lav_bvmix_cor_sc(
      rho = rho, fit_y1 = fit[[j]], fit_y2 = fit[[i]], wt = wt
    )
    out <- list(
      dx_rho = sc$dx_rho,
      dx_th_i = sc$dx_th_y2, dx_th_j = sc$dx_mu_y1,
      dx_sl_i = sc$dx_sl_y2, dx_sl_j = sc$dx_sl_y1,
      dx_var_i = NULL, dx_var_j = sc$dx_var_y1
    )
  } else {
    # polychoric
    sc <- lav_bvord_cor_sc(
      rho = rho, fit_y1 = fit[[i]], fit_y2 = fit[[j]], wt = wt
    )
    out <- list(
      dx_rho = sc$dx_rho,
      dx_th_i = sc$dx_th_y1, dx_th_j = sc$dx_th_y2,
      dx_sl_i = sc$dx_sl_y1, dx_sl_j = sc$dx_sl_y2,
      dx_var_i = NULL, dx_var_j = NULL
    )
  }
  out
}

# the A21 block: for each correlation (row), the crossproduct of its scores
# with the cross-derivatives w.r.t. the univariate parameters (columns);
# also returns the correlation scores sc_cor (unweighted) and the H21/H22
# blocks of the delta-rule matrix (rho -> cov metric, numeric variables)
lav_m84_a21 <- function(fit = NULL, ov_types = NULL,
                        cor_1 = NULL, var_1 = NULL,
                        wt = NULL, n = NULL, layout = NULL,
                        th_keep = NULL) {
  nvar <- layout$nvar
  pstar <- layout$pstar

  sc_cor <- matrix(0, n, pstar)
  a21 <- matrix(0, pstar, layout$a11_size)
  h22 <- diag(pstar) # for the delta rule
  h21 <- matrix(0, pstar, layout$a11_size)

  # for this one, we need new scores: for each F_ij (cor), the
  # scores with respect to the TH, VAR, ...
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      pstar_idx <- layout$pstar_1[i, j]

      pair <- lav_m84_pair_sc(
        fit = fit, ov_types = ov_types, i = i, j = j,
        rho = cor_1[i, j], wt = wt
      )

      # RHO
      if (is.null(wt)) {
        sc_cor[, pstar_idx] <- pair$dx_rho
      } else {
        sc_cor[, pstar_idx] <- pair$dx_rho / wt # unweight
      }

      # TH
      # (th_keep: skip the pseudo-threshold slots of empty categories --
      #  the fit has no scores for them; their a21 entries stay zero)
      a21[pstar_idx, layout$th_idx[[i]][th_keep[[i]]]] <-
        lav_mat_crossprod(sc_cor[, pstar_idx], pair$dx_th_i)
      a21[pstar_idx, layout$th_idx[[j]][th_keep[[j]]]] <-
        lav_mat_crossprod(sc_cor[, pstar_idx], pair$dx_th_j)

      # SL
      if (layout$nexo > 0L) {
        a21[pstar_idx, layout$sl_idx[[i]]] <-
          lav_mat_crossprod(sc_cor[, pstar_idx], pair$dx_sl_i)
        a21[pstar_idx, layout$sl_idx[[j]]] <-
          lav_mat_crossprod(sc_cor[, pstar_idx], pair$dx_sl_j)
      }

      # VAR (numeric variables only) + H21/H22 delta-rule entries
      num_i <- !is.null(pair$dx_var_i)
      num_j <- !is.null(pair$dx_var_j)
      if (num_i) {
        a21[pstar_idx, layout$var_idx[[i]]] <-
          lav_mat_crossprod(sc_cor[, pstar_idx], pair$dx_var_i)
      }
      if (num_j) {
        a21[pstar_idx, layout$var_idx[[j]]] <-
          lav_mat_crossprod(sc_cor[, pstar_idx], pair$dx_var_j)
      }
      if (num_i && num_j) {
        h21[pstar_idx, layout$var_idx[[i]]] <-
          (sqrt(var_1[j]) * cor_1[i, j]) / (2 * sqrt(var_1[i]))
        h21[pstar_idx, layout$var_idx[[j]]] <-
          (sqrt(var_1[i]) * cor_1[i, j]) / (2 * sqrt(var_1[j]))
        h22[pstar_idx, pstar_idx] <- sqrt(var_1[i]) * sqrt(var_1[j])
      } else if (num_i) {
        h21[pstar_idx, layout$var_idx[[i]]] <-
          cor_1[i, j] / (2 * sqrt(var_1[i]))
        h22[pstar_idx, pstar_idx] <- sqrt(var_1[i])
      } else if (num_j) {
        h21[pstar_idx, layout$var_idx[[j]]] <-
          cor_1[i, j] / (2 * sqrt(var_1[j]))
        h22[pstar_idx, pstar_idx] <- sqrt(var_1[j])
      }
    }
  }

  list(sc_cor = sc_cor, a21 = a21, h21 = h21, h22 = h22)
}

# the A11 block: 'sparse' (per-variable block-diagonal) version of the
# left-upper block of INNER (crossprod of the univariate scores)
lav_m84_a11 <- function(inner2 = NULL, layout = NULL) {
  a11 <- matrix(0, layout$a11_size, layout$a11_size)
  for (i in seq_len(layout$nvar)) {
    a11_idx <- c(
      layout$th_idx[[i]], layout$sl_idx[[i]],
      layout$var_idx[[i]][!is.na(layout$var_idx[[i]])]
    )
    a11[a11_idx, a11_idx] <- inner2[a11_idx, a11_idx]
  }
  a11
}

# the A22 block: diagonal of the correlation-score crossproducts
lav_m84_a22 <- function(sc_cor = NULL, wt = NULL) {
  pstar <- NCOL(sc_cor)
  a22 <- matrix(0, pstar, pstar)
  for (i in seq_len(pstar)) {
    if (is.null(wt)) {
      a22[i, i] <- sum(sc_cor[, i] * sc_cor[, i], na.rm = TRUE)
    } else {
      a22[i, i] <- sum(sc_cor[, i] * sc_cor[, i] / wt, na.rm = TRUE)
    }
  }
  a22
}

# invert the block-triangular bread matrix
#
# B <- rbind( cbind(A11,A12),
#            cbind(A21,A22) )
#
# B.inv = A11^{-1}                   0
#         -A22^{-1} A21 A11^{-1}     A22^{-1}
#
lav_m84_b_inv <- function(a11 = NULL, a21 = NULL, a22 = NULL) {
  # invert A11
  ginv_a11 <- FALSE
  a11_inv <- try(solve(a11), silent = TRUE)
  if (inherits(a11_inv, "try-error")) {
    # brute force
    a11_inv <- MASS::ginv(a11)
    ginv_a11 <- TRUE
  }

  # invert A22 (diagonal)
  ginv_a22 <- FALSE
  da22 <- diag(a22)
  if (any(da22 == 0)) {
    ginv_a22 <- TRUE
    a22_inv <- MASS::ginv(a22)
  } else {
    a22_inv <- a22
    diag(a22_inv) <- 1 / da22
  }

  # lower-left block
  a21_inv <- -a22_inv %*% a21 %*% a11_inv

  # upper-right block remains zero
  a12_inv <- matrix(0, NROW(a11), NCOL(a22))

  # construct B.inv
  b_inv <- rbind(
    cbind(a11_inv, a12_inv),
    cbind(a21_inv, a22_inv)
  )

  list(b_inv = b_inv, ginv_a11 = ginv_a11, ginv_a22 = ginv_a22)
}
