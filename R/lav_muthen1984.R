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
muthen1984 <- function(Data = NULL,
                       ov.names = NULL,
                       ov.types = NULL,
                       ov.levels = NULL,
                       ov.names.x = character(0L),
                       eXo = NULL,
                       wt = NULL,
                       verbose = FALSE,
                       WLS.W = TRUE,
                       zero.add = c(0.5, 0.0),
                       zero.keep.margins = TRUE,
                       zero.cell.warn = FALSE,
                       zero.cell.tables = TRUE,
                       group = 1L) { # group only for error messages

  # just in case Data is a vector
  Data <- as.matrix(Data)

  nvar <- NCOL(Data)
  N <- NROW(Data)
  num.idx <- which(ov.types == "numeric")
  ord.idx <- which(ov.types == "ordered")
  nexo <- length(ov.names.x)
  if (nexo > 0L) stopifnot(NCOL(eXo) == nexo)
  pstar <- nvar * (nvar - 1) / 2

  if (verbose) {
    cat("\nPreparing for WLS estimation -- STEP 1 + 2\n")
    cat("Number of endogenous variables: ", nvar, "\n")
    cat("Endogenous variable names:\n")
    print(ov.names)
    cat("\n")
    cat("Endogenous ov types:\n")
    print(ov.types)
    cat("\n")
    cat("Endogenous ov levels:\n ")
    print(ov.levels)
    cat("\n")
    cat("Number of exogenous variables: ", nexo, "\n")
    cat("Exogenous variable names:\n")
    print(ov.names.x)
    cat("\n")
  }

  step1 <- lav_samplestats_step1(
    Y = Data, wt = wt, ov.names = ov.names,
    ov.types = ov.types, ov.levels = ov.levels, ov.names.x = ov.names.x,
    eXo = eXo, scores.flag = WLS.W, group = group
  )

  FIT <- step1$FIT
  TH <- step1$TH
  TH.NOX <- step1$TH.NOX
  TH.IDX <- step1$TH.IDX
  TH.NAMES <- step1$TH.NAMES
  VAR <- step1$VAR
  SLOPES <- step1$SLOPES
  SC.TH <- step1$SC.TH
  SC.SL <- step1$SC.SL
  SC.VAR <- step1$SC.VAR
  th.start.idx <- step1$th.start.idx
  th.end.idx <- step1$th.end.idx

  # rm SC.VAR columns from ordinal variables
  if (WLS.W && length(ord.idx) > 0L) {
    SC.VAR <- SC.VAR[, -ord.idx, drop = FALSE]
  }


  if (verbose) {
    cat("STEP 1: univariate statistics\n")
    cat("Threshold + means:\n")
    TTHH <- unlist(TH)
    names(TTHH) <- unlist(TH.NAMES)
    print(TTHH)
    cat("Slopes (if any):\n")
    colnames(SLOPES) <- ov.names.x
    rownames(SLOPES) <- ov.names
    print(SLOPES)
    cat("Variances:\n")
    names(VAR) <- ov.names
    print(unlist(VAR))
  }

  # stage two -- correlations

  if (verbose) cat("\n\nSTEP 2: covariances/correlations:\n")
  COR <- lav_samplestats_step2(
    UNI = FIT, wt = wt, ov.names = ov.names,
    zero.add = zero.add,
    zero.keep.margins = zero.keep.margins,
    zero.cell.warn = zero.cell.warn,
    zero.cell.tables = zero.cell.tables
  )
  empty.cell.tables <- attr(COR, "zero.cell.tables")
  attr(COR, "zero.cell.tables") <- NULL

  if (verbose) {
    colnames(COR) <- rownames(COR) <- ov.names
    print(COR)
  }

  if (!WLS.W) { # we do not need the asymptotic variance matrix
    if (any("numeric" %in% ov.types)) {
      COV <- cor2cov(R = COR, sds = sqrt(unlist(VAR)))
    } else {
      COV <- COR
    }
    out <- list(
      TH = TH, SLOPES = SLOPES, VAR = VAR, COR = COR, COV = COV,
      SC = NULL, TH.NOX = TH.NOX, TH.NAMES = TH.NAMES, TH.IDX = TH.IDX,
      INNER = NULL, A11 = NULL, A12 = NULL, A21 = NULL, A22 = NULL,
      WLS.W = NULL, H = NULL, zero.cell.tables = matrix("", 0, 2)
    )
    return(out)
  }


  # stage three -- WLS.W
  SC.COR <- matrix(0, N, pstar)
  PSTAR <- matrix(0, nvar, nvar)
  PSTAR[lav_matrix_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar

  A11.size <- NCOL(SC.TH) + NCOL(SC.SL) + NCOL(SC.VAR)








  # A21
  A21 <- matrix(0, pstar, A11.size)
  H22 <- diag(pstar) # for the delta rule
  H21 <- matrix(0, pstar, A11.size)
  # for this one, we need new scores: for each F_ij (cor), the
  # scores with respect to the TH, VAR, ...
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      pstar.idx <- PSTAR[i, j]
      th.idx_i <- th.start.idx[i]:th.end.idx[i]
      th.idx_j <- th.start.idx[j]:th.end.idx[j]
      if (nexo > 0L) {
        sl.idx_i <- NCOL(SC.TH) + seq(i, by = nvar, length.out = nexo)
        sl.idx_j <- NCOL(SC.TH) + seq(j, by = nvar, length.out = nexo)

        var.idx_i <- NCOL(SC.TH) + NCOL(SC.SL) + match(i, num.idx)
        var.idx_j <- NCOL(SC.TH) + NCOL(SC.SL) + match(j, num.idx)
      } else {
        var.idx_i <- NCOL(SC.TH) + match(i, num.idx)
        var.idx_j <- NCOL(SC.TH) + match(j, num.idx)
      }
      if (ov.types[i] == "numeric" && ov.types[j] == "numeric") {
        SC.COR.UNI <- lav_bvreg_cor_scores(
          rho = COR[i, j],
          fit.y1 = FIT[[i]],
          fit.y2 = FIT[[j]],
          wt = wt
        )

        # RHO
        if (is.null(wt)) {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho
        } else {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho / wt # unweight
        }

        # TH
        A21[pstar.idx, th.idx_i] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.mu.y1
          )
        A21[pstar.idx, th.idx_j] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.mu.y2
          )
        # SL
        if (nexo > 0L) {
          A21[pstar.idx, sl.idx_i] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y1
            )
          A21[pstar.idx, sl.idx_j] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y2
            )
        }
        # VAR
        A21[pstar.idx, var.idx_i] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.var.y1
          )
        A21[pstar.idx, var.idx_j] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.var.y2
          )
        # H21 only needed for VAR
        H21[pstar.idx, var.idx_i] <-
          (sqrt(VAR[j]) * COR[i, j]) / (2 * sqrt(VAR[i]))
        H21[pstar.idx, var.idx_j] <-
          (sqrt(VAR[i]) * COR[i, j]) / (2 * sqrt(VAR[j]))
        H22[pstar.idx, pstar.idx] <- sqrt(VAR[i]) * sqrt(VAR[j])
      } else if (ov.types[i] == "numeric" && ov.types[j] == "ordered") {
        SC.COR.UNI <- lav_bvmix_cor_scores(
          rho = COR[i, j],
          fit.y1 = FIT[[i]],
          fit.y2 = FIT[[j]],
          wt = wt
        )
        # RHO
        if (is.null(wt)) {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho
        } else {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho / wt # unweight
        }

        # TH
        A21[pstar.idx, th.idx_i] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.mu.y1
          )
        A21[pstar.idx, th.idx_j] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.th.y2
          )
        # SL
        if (nexo > 0L) {
          A21[pstar.idx, sl.idx_i] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y1
            )
          A21[pstar.idx, sl.idx_j] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y2
            )
        }
        # VAR
        A21[pstar.idx, var.idx_i] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.var.y1
          )
        # H21 only need for VAR
        H21[pstar.idx, var.idx_i] <- COR[i, j] / (2 * sqrt(VAR[i]))
        H22[pstar.idx, pstar.idx] <- sqrt(VAR[i])
      } else if (ov.types[j] == "numeric" && ov.types[i] == "ordered") {
        SC.COR.UNI <- lav_bvmix_cor_scores(
          rho = COR[i, j],
          fit.y1 = FIT[[j]],
          fit.y2 = FIT[[i]],
          wt = wt
        )
        # RHO
        if (is.null(wt)) {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho
        } else {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho / wt # unweight
        }

        # TH
        A21[pstar.idx, th.idx_j] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.mu.y1
          )
        A21[pstar.idx, th.idx_i] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.th.y2
          )
        # SL
        if (nexo > 0L) {
          A21[pstar.idx, sl.idx_j] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y1
            )
          A21[pstar.idx, sl.idx_i] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y2
            )
        }
        # VAR
        A21[pstar.idx, var.idx_j] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.var.y1
          )
        # H21 only for VAR
        H21[pstar.idx, var.idx_j] <- COR[i, j] / (2 * sqrt(VAR[j]))
        H22[pstar.idx, pstar.idx] <- sqrt(VAR[j])
      } else if (ov.types[i] == "ordered" && ov.types[j] == "ordered") {
        # polychoric correlation
        SC.COR.UNI <- lav_bvord_cor_scores(
          rho = COR[i, j],
          fit.y1 = FIT[[i]],
          fit.y2 = FIT[[j]],
          wt = wt
        )
        # RHO
        if (is.null(wt)) {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho
        } else {
          SC.COR[, pstar.idx] <- SC.COR.UNI$dx.rho / wt # unweight
        }

        # TH
        A21[pstar.idx, th.idx_i] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.th.y1
          )
        A21[pstar.idx, th.idx_j] <-
          lav_matrix_crossprod(
            SC.COR[, pstar.idx],
            SC.COR.UNI$dx.th.y2
          )
        # SL
        if (nexo > 0L) {
          A21[pstar.idx, sl.idx_i] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y1
            )
          A21[pstar.idx, sl.idx_j] <-
            lav_matrix_crossprod(
              SC.COR[, pstar.idx],
              SC.COR.UNI$dx.sl.y2
            )
        }
        # NO VAR
      }
    }
  }
  if (!is.null(wt)) {
    SC.COR <- SC.COR * wt # reweight
  }





  # stage three

  SC <- cbind(SC.TH, SC.SL, SC.VAR, SC.COR)
  INNER <- lav_matrix_crossprod(SC)

  # A11
  # new approach (2 June 2012): A11 is just a 'sparse' version of
  # (the left upper block of) INNER
  A11 <- matrix(0, A11.size, A11.size)
  if (!is.null(wt)) {
    INNER2 <- lav_matrix_crossprod(SC / wt, SC)
  } else {
    INNER2 <- INNER
  }
  for (i in 1:nvar) {
    th.idx <- th.start.idx[i]:th.end.idx[i]
    sl.idx <- integer(0L)
    var.idx <- integer(0L)
    if (nexo > 0L) {
      sl.idx <- NCOL(SC.TH) + seq(i, by = nvar, length.out = nexo)
      # sl.end.idx <- (i*nexo); sl.start.idx <- (i-1L)*nexo + 1L
      # sl.idx <- NCOL(SC.TH) + (sl.start.idx:sl.end.idx)
    }
    if (ov.types[i] == "numeric") {
      var.idx <- NCOL(SC.TH) + NCOL(SC.SL) + match(i, num.idx)
    }
    a11.idx <- c(th.idx, sl.idx, var.idx)
    A11[a11.idx, a11.idx] <- INNER2[a11.idx, a11.idx]
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
  A22 <- matrix(0, pstar, pstar)
  for (i in seq_len(pstar)) {
    if (is.null(wt)) {
      A22[i, i] <- sum(SC.COR[, i] * SC.COR[, i], na.rm = TRUE)
    } else {
      A22[i, i] <- sum(SC.COR[, i] * SC.COR[, i] / wt, na.rm = TRUE)
    }
  }

  # A12 (zero)
  A12 <- matrix(0, NROW(A11), NCOL(A22))


  # B <- rbind( cbind(A11,A12),
  #            cbind(A21,A22) )

  # we invert B as a block-triangular matrix (0.5-23)
  #
  # B.inv = A11^{-1}                   0
  #         -A22^{-1} A21 A11^{-1}     A22^{-1}
  #

  # invert A
  A11.inv <- try(solve(A11), silent = TRUE)
  if (inherits(A11.inv, "try-error")) {
    # brute force
    A11.inv <- MASS::ginv(A11)
    lav_msg_warn(gettext("trouble constructing W matrix;
                         used generalized inverse for A11 submatrix"))
  }

  # invert
  da22 <- diag(A22)
  if (any(da22 == 0)) {
    lav_msg_warn(gettext("trouble constructing W matrix;
                         used generalized inverse for A22 submatrix"))
    A22.inv <- MASS::ginv(A22)
  } else {
    A22.inv <- A22
    diag(A22.inv) <- 1 / da22
  }

  # lower-left block
  A21.inv <- -A22.inv %*% A21 %*% A11.inv

  # upper-left block remains zero
  A12.inv <- A12

  # construct B.inv
  B.inv <- rbind(
    cbind(A11.inv, A12.inv),
    cbind(A21.inv, A22.inv)
  )


  #  weight matrix (correlation metric)
  WLS.W <- B.inv %*% INNER %*% t(B.inv)

  # COV matrix?
  if (any("numeric" %in% ov.types)) {
    COV <- cor2cov(R = COR, sds = sqrt(unlist(VAR)))

    # construct H matrix to apply delta rule (for the tranformation
    # of rho_ij to cov_ij)
    H11 <- diag(NROW(A11))
    H12 <- matrix(0, NROW(A11), NCOL(A22))
    # H22 and H21 already filled in
    H <- rbind(
      cbind(H11, H12),
      cbind(H21, H22)
    )

    WLS.W <- H %*% WLS.W %*% t(H)
  } else {
    COV <- COR
    H <- diag(NCOL(WLS.W))
  }

  # reverse sign numeric TH (because we provide -mu in WLS.obs)
  # (WOW, it took me a LOOONGGG time to realize this!)
  # YR 16 July 2012

  # NOTE: prior to 0.5-17, we used num.idx (instead of NUM.idx)
  # which is WRONG if we have more than one threshold per variable
  # (thanks to Sacha Epskamp for spotting this!)
  if (length(num.idx) > 0L) {
    NUM.idx <- which(unlist(TH.IDX) == 0L)
    WLS.W[NUM.idx, ] <- -WLS.W[NUM.idx, ]
    WLS.W[, NUM.idx] <- -WLS.W[, NUM.idx]
  }


  out <- list(
    TH = TH, SLOPES = SLOPES, VAR = VAR, COR = COR, COV = COV,
    SC = SC, TH.NOX = TH.NOX, TH.NAMES = TH.NAMES, TH.IDX = TH.IDX,
    INNER = INNER, A11 = A11, A12 = A12, A21 = A21, A22 = A22,
    WLS.W = WLS.W, H = H,
    zero.cell.tables = empty.cell.tables
  )
  out
}
