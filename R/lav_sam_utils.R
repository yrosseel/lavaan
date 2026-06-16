# utility functions for the sam() function
# YR 4 April 2023

# clamped small-sample correction factor alpha/(n - 1), forced into [0, 1]
lav_sam_alpha_n1 <- function(alpha_correction = 0L, n = 20L) {
  alpha_n1 <- alpha_correction / (n - 1)
  if (!is.finite(alpha_n1)) {
    alpha_n1 <- 0
  }
  min(max(alpha_n1, 0), 1)
}

# row/column index pairs and labels for the elements of (eta* %x% eta*),
# where lv_names already includes the leading "..int.." intercept entry;
# the first nfac labels are the plain (first-order) lv names
lav_sam_lvnames2 <- function(lv_names) {
  nfac <- length(lv_names)
  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)
  names2 <- paste(lv_names[idx1], lv_names[idx2], sep = ":")
  names2[seq_len(nfac)] <- lv_names
  list(idx1 = idx1, idx2 = idx2, names = names2)
}

# construct 'mapping matrix' M using either "ML", "GLS" or "ULS" method
#
# by construction, M %*% LAMBDA = I (the identity matrix)
lav_sam_mapping_mat <- function(mm_lambda = NULL, mm_theta = NULL,
                                   s = NULL, s_inv = NULL,
                                   method = "ML") {

  method <- toupper(method)

  # catch empty columns in LAMBDA (eg higher order!)
  lambda_orig <- mm_lambda
  empty_idx <- which(apply(mm_lambda, 2L, function(x) all(x == 0)))
  if (length(empty_idx) > 0L) {
    mm_lambda <- lambda_orig[, -empty_idx, drop = FALSE]
  }

  # ULS
  # M == solve( t(LAMBDA) %*% LAMBDA ) %*% t(LAMBDA)
  #   == MASS:::ginv(LAMBDA)
  if (method == "ULS") {
    m <- try(tcrossprod(solve(crossprod(mm_lambda)), mm_lambda),
      silent = TRUE
    )
    if (inherits(m, "try-error")) {
      lav_msg_warn(gettext(
        "cannot invert crossprod(LAMBDA); using generalized inverse"))
      m <- MASS::ginv(mm_lambda)
    }


    # GLS
    # M == solve( t(LAMBDA) %*% S.inv %*% LAMBDA ) %*% t(LAMBDA) %*% S.inv
  } else if (method == "GLS") {
    if (is.null(s_inv)) {
      s_inv <- try(solve(s), silent = TRUE)
    }
    if (inherits(s_inv, "try-error")) {
      lav_msg_warn(gettext("S is not invertible; switching to ULS method"))
      m <- lav_sam_mapping_mat(mm_lambda = mm_lambda, method = "ULS")
    } else {
      t_lsinv <- t(mm_lambda) %*% s_inv
      t_lsinv_l <- t_lsinv %*% mm_lambda
      m <- try(solve(t_lsinv_l, t_lsinv), silent = TRUE)
      if (inherits(m, "try-error")) {
        lav_msg_warn(gettext("problem constructing mapping matrix;
                              switching to generalized inverse"))
        m <- MASS::ginv(t_lsinv_l) %*% t_lsinv
      }
    }

    # ML
    # M == solve(t(LAMBDA) %*% THETA.inv %*% LAMBDA) %*% t(LAMBDA) %*% THETA.inv
  } else if (method == "ML") {
    # Problem: if THETA has zero elements on the diagonal, we cannot invert

    # As we do not have access to Sigma(.inv), we cannot
    # use the trick as in lavPredict() (where we replace THETA.inv by
    # Sigma.inv)

    # old method (<0.6-16): remove rows/cols with zero values
    # on the diagonal of THETA, invert the submatrix and
    # set the '0' diagonal elements to one. This resulted in (somewhat)
    # distorted results.

    # new in 0.6-16: we use the Wall & Amemiya (2000) method using
    # the so-called 'T' transformation

    # new in 0.6-18: if we cannot use the marker method (eg growth models),
    # use the 'old' method anyway: remove zero rows/cols and invert submatrix

    zero_theta_idx <- which(abs(diag(mm_theta)) < 1e-4) # be conservative
    if (length(zero_theta_idx) == 0L) {
      # ok, no zero diagonal elements: try to invert THETA
      if (lav_mat_is_diagonal(mm_theta)) {
        theta_inv <- diag(1 / diag(mm_theta), nrow = nrow(mm_theta))
      } else {
        theta_inv <- try(solve(mm_theta), silent = TRUE)
        if (inherits(theta_inv, "try-error")) {
          theta_inv <- NULL
        }
      }
    } else {
      # see if we can use marker method
      marker_idx <- lav_utils_get_marker(mm_lambda = mm_lambda, std_lv = TRUE)
      if (any(is.na(marker_idx))) {
        theta_inv <- try(lav_mat_sym_inverse(mm_theta), silent = TRUE)
        if (inherits(theta_inv, "try-error")) {
          theta_inv <- NULL
        } else {
          diag(theta_inv)[zero_theta_idx] <- 1
        }
      } else {
        # try tmat method later on
        theta_inv <- NULL
      }
    }

    # could we invert THETA?
    if (!is.null(theta_inv)) {
      # ha, all is good; compute M the usual way
      t_lti <- t(mm_lambda) %*% theta_inv
      t_lti_l <- t_lti %*% mm_lambda
      m <- try(solve(t_lti_l, t_lti), silent = TRUE)
      if (inherits(m, "try-error")) {
        lav_msg_warn(gettext(
          "problem constructing ML mapping matrix; switching to ULS"))
        m <- lav_sam_mapping_mat(mm_lambda = mm_lambda, method = "ULS")
      }
    } else {
      # use W&A2000's method using the 'T' transformation
      m <- try(lav_sam_mapping_mat_tmat(
        mm_lambda = mm_lambda,
        mm_theta = mm_theta
      ), silent = TRUE)
      if (inherits(m, "try-error")) {
        lav_msg_warn(gettext(
          "problem constructing ML mapping matrix; switching to ULS"))
        m <- lav_sam_mapping_mat(mm_lambda = mm_lambda, method = "ULS")
      }
    }
  } # ML

  # empty.idx?
  if (length(empty_idx) > 0L) {
    m_full <- m
    m <- matrix(0, nrow = ncol(lambda_orig), ncol = ncol(m_full))
    m[-empty_idx, ] <- m_full
  }

  m
}

# use 'T' transformation to create the 'Bartlett/ML' mapping matrix
# see Wall & Amemiya (2000) eq (7)
# see also Fuller 1987 page 357 (where T is called H), and page 364

# although Fuller/W&A always assumed that THETA is diagonal,
# their method seems to work equally well for non-diagonal THETA
#
# in our implementation:
# - we do NOT reorder the rows of LAMBDA
# - if std.lv = TRUE, we first rescale to 'create' marker indicators
#   and then rescale back at the end
#
lav_sam_mapping_mat_tmat <- function(mm_lambda = NULL,
                                        mm_theta = NULL,
                                        marker_idx = NULL,
                                        std_lv = NULL) {

  mm_lambda <- as.matrix.default(mm_lambda)

  # catch empty columns in LAMBDA (eg higher order!)
  lambda_orig <- mm_lambda
  empty_idx <- which(apply(mm_lambda, 2L, function(x) all(x == 0)))
  if (length(empty_idx) > 0L) {
    mm_lambda <- lambda_orig[, -empty_idx, drop = FALSE]
  }

  # nvar <- nrow(mm_lambda)
  # nfac <- ncol(mm_lambda)

  # do we have marker.idx?
  if (is.null(marker_idx)) {
    # 'marker' indicator has a single non-zero element in a row
    marker_idx <- lav_utils_get_marker(mm_lambda = mm_lambda, std_lv = TRUE)
    if (any(is.na(marker_idx))) {
      lav_msg_stop(gettext("no clear markers in LAMBDA matrix"))
    }
  }

  # std.lv TRUE or FALSE?
  if (is.null(std_lv)) {
    std_lv <- FALSE
    if (any(diag(mm_lambda[marker_idx, , drop = FALSE]) != 1)) {
      std_lv <- TRUE
    }
  }

  # if std.lv = TRUE, rescale
  if (std_lv) {
    marker <- mm_lambda[marker_idx, , drop = FALSE]
    marker_inv <- 1 / diag(marker)
    mm_lambda <- t(t(mm_lambda) * marker_inv)
  }

  # compute 'T' matrix
  tmat <- lav_sam_tmat(
    mm_lambda = mm_lambda, mm_theta = mm_theta,
    marker_idx = marker_idx
  )

  # ML mapping matrix
  m <- tmat[marker_idx, , drop = FALSE]

  if (std_lv) {
    m <- m * marker_inv
  }

  # empty.idx?
  if (length(empty_idx) > 0L) {
    m_full <- m
    m <- matrix(0, nrow = ncol(lambda_orig), ncol = ncol(m_full))
    m[-empty_idx, ] <- m_full
  }

  m
}

# create 'T' matrix (tmat) for T-transformation
#
# Notes: - here we assume that LAMBDA has unity markers (no std.lv = TRUE)
#        - TMAT is NOT symmetric!
#        - Yc %*% t(TMAT) transforms the data in such a way that we get:
#           1)  Bartlett factor scores in the marker columns
#           2) 'V' values in the non-marker columns, where:
#               V = Yc - Yc[,marker.idx] %*% t(LAMBDA)
#
lav_sam_tmat <- function(mm_lambda = NULL,
                         mm_theta = NULL,
                         marker_idx = NULL) {
  mm_lambda <- as.matrix.default(mm_lambda)
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)

  # do we have marker.idx?
  if (is.null(marker_idx)) {
    # 'marker' indicator has a single 1 element in a row
    marker_idx <- lav_utils_get_marker(mm_lambda = mm_lambda, std_lv = FALSE)
    if (any(is.na(marker_idx))) {
      lav_msg_stop(gettext("no clear markers in LAMBDA matrix"))
    }
  }

  # construct 'm_c' matrix
  c2 <- diag(nvar)
  c2[, marker_idx] <- -1 * mm_lambda
  m_c <- c2[-marker_idx, , drop = FALSE]

  # compute Sigma.ve and Sigma.vv
  sigma_ve <- m_c %*% mm_theta
  # Sigma.vv <- m_c %*% THETA %*% t(m_c)
  sigma_vv <- sigma_ve %*% t(m_c)

  # construct 'Gamma' (=gamma_1) (and gamma2) matrix
  # Gamma <- (t(Sigma.ve) %*% solve(Sigma.vv))[marker.idx,, drop = FALSE]
  gamma_1 <- try(t(solve(sigma_vv, sigma_ve)[, marker_idx, drop = FALSE]),
    silent = TRUE
  )
  if (inherits(gamma_1, "try-error")) {
    tmp <- t(sigma_ve) %*% MASS::ginv(sigma_vv)
    gamma_1 <- tmp[marker_idx, , drop = FALSE]
  }
  gamma2 <- matrix(0, nfac, nvar)
  gamma2[, -marker_idx] <- gamma_1
  gamma2[, marker_idx] <- diag(nfac)

  # transformation matrix 'T' (we call it here 'Tmat')
  tmat <- matrix(0, nvar, nvar)
  tmat[-marker_idx, ] <- m_c
  tmat[marker_idx, ] <- -gamma2 %*% c2

  tmat
}


# compute VETA
# - if alpha.correction == 0     -> same as local SAM (or MOC)
# - if alpha.correction == (N-1) -> same as FSR+Bartlett
lav_sam_veta <- function(m = NULL, s = NULL, mm_theta = NULL,
                         alpha_correction = 0L, lambda_correction = TRUE,
                         n = 20L, dummy_lv_idx = integer(0L), extra = FALSE) {

  # catch empty rows in M (higher-order?)
  m_orig <- m
  empty_idx <- which(apply(m_orig, 1L, function(x) all(x == 0)))
  if (length(empty_idx) > 0L) {
    m <- m_orig[-empty_idx, , drop = FALSE]
  }

  # MSM
  msm <- m %*% s %*% t(m)

  # MTM
  mtm <- m %*% mm_theta %*% t(m)

  # note: MTM no longer needs to be regularized (force-pd) here:
  # lav_mat_sym_diff_smallest_root() handles a singular (or even zero)
  # MTM natively, by regressing out its null space (which includes the
  # rows/columns of any dummy lvs, and any empty theta elements)
  if (all(abs(mtm) < sqrt(.Machine$double.eps))) {
    # all zero? certainly no need for alpha or lambda.correction
    alpha_correction <- 0L
    lambda_correction <- FALSE
  }

  # apply small sample correction (if requested)
  alpha <- alpha_correction
  if (alpha_correction > 0) {
    alpha_n1 <- lav_sam_alpha_n1(alpha_correction, n)
    mtm <- (1 - alpha_n1) * mtm
  }

  lambda <- lambda_star <- +Inf
  if (lambda_correction) {
    # use Fuller (1987) approach to ensure VETA is positive
    lambda <- try(lav_mat_sym_diff_smallest_root(msm, mtm),
      silent = TRUE
    )
    if (inherits(lambda, "try-error")) {
      lav_msg_warn(gettext("failed to compute lambda"))
      veta <- msm - mtm # and hope for the best
    } else {
      cutoff <- 1 + 1 / (n - 1)
      if (lambda < cutoff) {
        lambda_star <- lambda - 1 / (n - 1)
        veta <- msm - lambda_star * mtm
      } else {
        veta <- msm - mtm
      }
    }
  } else {
    veta <- msm - mtm
  }

  # empty.idx?
  if (length(empty_idx) > 0L) {
    msm_full <- msm
    mtm_full <- mtm
    veta_full <- veta
    nfac_orig <- nrow(m_orig)

    msm <- mtm <- veta <- matrix(0, nrow = nfac_orig, ncol = nfac_orig)
    msm[-empty_idx, -empty_idx] <- msm_full
    mtm[-empty_idx, -empty_idx] <- mtm_full
    veta[-empty_idx, -empty_idx] <- veta_full
  }

  # extra attributes?
  if (extra) {
    attr(veta, "lambda") <- lambda
    attr(veta, "alpha") <- alpha
    attr(veta, "lambda.star") <- lambda_star
    attr(veta, "MSM") <- msm
    attr(veta, "MTM") <- mtm
  }

  veta
}

# compute EETA = E(Eta) = M %*% [YBAR - NU]
lav_sam_eeta <- function(m = NULL, ybar = NULL, mm_nu = NULL) {
  m %*% (ybar - mm_nu)
}

# compute veta including quadratic/interaction terms
lav_sam_veta2 <- function(fs = NULL, m = NULL,
                          veta = NULL, eeta = NULL, mm_theta = NULL,
                          lv_names = NULL,
                          lv_int_names = NULL,
                          dummy_lv_names = character(0L),
                          alpha_correction = 0L,
                          lambda_correction = TRUE,
                          lambda1 = 1,
                          fs_outlier_idx = integer(0L),
                          return_fs = FALSE,
                          return_cov_iveta2 = TRUE,
                          extra = FALSE) {

  # small utility function: var() divided by N
  varn <- function(x, n) {
    var(x, use = "pairwise.complete.obs") * (n - 1) / n
  }

  if (length(lv_int_names) == 0L) {
    lav_msg_stop(gettext("lv.int.names is empty: no lv quadratic/interaction
                         terms are provided"))
  }

  if (is.null(lv_names)) {
    lv_names <- paste("eta", seq_len(ncol(fs)), sep = "")
  }

  # MTM
  mtm <- m %*% mm_theta %*% t(m)

  # note: MTM no longer needs to be regularized (force-pd) here:
  # lav_mat_sym_diff_smallest_root() handles a singular var.error matrix
  # natively, by regressing out its null space (which includes the
  # rows/columns of any dummy lvs); this also keeps MTM identical to the
  # (raw) version used in lav_sam_gamma_add()

  # augment to include intercept
  fs <- cbind(1, fs)
  n <- nrow(fs)
  mtm <- lav_mat_bdiag(0, mtm)
  veta <- lav_mat_bdiag(0, veta)
  eeta <- c(1, eeta)
  lv_names <- c("..int..", lv_names)
  nfac <- ncol(fs)

  tmp <- lav_sam_lvnames2(lv_names)
  idx1 <- tmp$idx1
  idx2 <- tmp$idx2
  names_1 <- tmp$names

  fs2 <- fs[, idx1] * fs[, idx2]

  k_nfac <- lav_mat_com(nfac, nfac)
  ik <- diag(nfac * nfac) + k_nfac

  eeta <- as.matrix(drop(eeta))
  vetak_mtm <- veta %x% mtm

  # normal version (for now):
  gamma_me22 <- ik %*% (mtm %x% mtm)

  # ingredients (normal ME case)
  var_fs2 <- varn(fs2, n)
  var_etak_me <- (tcrossprod(eeta) %x% mtm + vetak_mtm)
  var_mek_eta <- lav_mat_com_pre_post(var_etak_me)
  var_me2 <- gamma_me22

  cov_etak_me_mek_eta <- lav_mat_com_post(var_etak_me)
  cov_mek_eta_etak_me <- t(cov_etak_me_mek_eta)

  var_error <- (var_etak_me + var_mek_eta + cov_etak_me_mek_eta
    + cov_mek_eta_etak_me + var_me2)

  # select only what we need
  colnames(var_fs2) <- rownames(var_fs2) <- names_1
  colnames(var_error) <- rownames(var_error) <- names_1
  colnames(fs2) <- names_1
  lv_keep <- c(lv_names[-1], lv_int_names)
  var_fs2 <- var_fs2[lv_keep, lv_keep]
  var_error <- var_error[lv_keep, lv_keep]
  fs_mean <- colMeans(fs2, na.rm = TRUE)
  names(fs_mean) <- names_1
  fs2_mean <- fs_mean # all of them
  fs_mean <- fs_mean[lv_keep]

  # apply small sample correction (if requested)
  alpha <- alpha_correction
  alpha_n1 <- 0
  if (alpha_correction > 0) {
    alpha_n1 <- lav_sam_alpha_n1(alpha_correction, n)
    var_error <- (1 - alpha_n1) * var_error
  }

  lambda <- +Inf
  lambda_star <- 1
  if (lambda_correction) {
    # use Fuller (1987) approach to ensure VETA2 is positive
    lambda <- try(lav_mat_sym_diff_smallest_root(
      var_fs2,
      var_error
    ), silent = TRUE)
    if (inherits(lambda, "try-error")) {
      lav_msg_warn(gettext("failed to compute lambda"))
      veta2 <- var_fs2 - var_error # and hope for the best
    } else {
      cutoff <- 1 + 2 / n # be more conservative for VETA2
      if (lambda < cutoff) {
        lambda_star <- max(c(0, lambda - ncol(var_fs2) / (n - 1)))
        veta2 <- var_fs2 - lambda_star * var_error
      } else {
        veta2 <- var_fs2 - var_error
      }
    }
  } else {
    veta2 <- var_fs2 - var_error
  }

  # new in 0.6-20: to compute Gamma.eta
  # the casewise contributions are defined so that their average over the
  # observations equals *exactly* the augmented summary statistics
  # (EETA2, VETA2) that are used in the second step; this requires the
  # effective *first-order* lambda1 inside the decomposition (because
  # EETA2 and var.error are built from the lambda1-corrected first-order
  # VETA), while the *second-order* lambda.star multiplies the whole
  # correction term:
  #   ieeta2_i = fs2[i, ] - lambda1 * vec(MTM)
  #   (the terms involving the mean factor scores cancel out exactly)
  # and
  #   iveta2_i = tcrossprod(fs2[i, ] - fs2.mean) - lambda2.eff * tmp_i
  # where lambda2.eff = lambda.star * (1 - alpha/(n-1)) is the effective
  # multiplier of the *unscaled* var.error (the alpha correction rescales
  # var.error before the lambda step), and
  #       tmp_i = ((F_i - lambda1 * MTM) %x% MTM) +
  #               (MTM %x% (F_i - lambda1 * MTM)) +
  #               lav_mat_com_post((F_i - lambda1 * MTM) %x% MTM) +
  #               lav_mat_com_pre((F_i - lambda1 * MTM) %x% MTM) +
  #               (I + K) %*% (MTM %x% MTM)
  # with F_i = tcrossprod(f_i); both selected in the rows/columns we keep
  # note: cov(iveta2_1) is invariant to the value of lambda1 (the lambda1
  # terms are constant across the observations); but lambda1 does matter
  # for the jacobian of the average casewise contribution computed in
  # lav_sam_gamma_add()
  # rewritten 12 June 2026 to avoid the loop over the observations: every
  # kept element is a linear function of the second-order factor scores
  # fs2[i, ]; elementwise, for kept elements a = (r1, s1) and b = (r2, s2)
  # of (eta* %x% eta*):
  #   tmp_i[a, b] = F_i[r1, r2] MTM[s1, s2] + MTM[r1, r2] F_i[s1, s2] +
  #                 F_i[r1, s2] MTM[s1, r2] + F_i[s1, r2] MTM[r1, s2] +
  #                 (1 - 2 * lambda1) * (MTM[r1, r2] MTM[s1, s2] +
  #                                      MTM[s1, r2] MTM[r1, s2])
  # where F_i[r, s] = fs2[i, (r - 1) * nfac + s]
  if (return_cov_iveta2) {
    lambda2_eff <- lambda_star * (1 - alpha_n1)
    keep_idx <- match(lv_keep, names_1)
    nkeep <- length(keep_idx)
    # row/column pairs of the kept elements of (eta* %x% eta*)
    i1k <- idx1[keep_idx]
    i2k <- idx2[keep_idx]

    # casewise contributions to E(eta* %x% eta*)[keep]
    part1 <- t(t(fs2[, keep_idx, drop = FALSE]) -
               lambda1 * lav_mat_vec(mtm)[keep_idx])

    # casewise contributions to vech(Var(eta* %x% eta*)[keep, keep])
    # vech row (pa) and column (pb) indices over the kept elements
    tmp_mat <- matrix(seq_len(nkeep), nkeep, nkeep)
    pa <- tmp_mat[lower.tri(tmp_mat, diag = TRUE)]
    pb <- t(tmp_mat)[lower.tri(tmp_mat, diag = TRUE)]
    r1 <- i1k[pa]; s1 <- i2k[pa]
    r2 <- i1k[pb]; s2 <- i2k[pb]

    # fs2 columns holding F_i[r1, r2], F_i[s1, s2], F_i[r1, s2], F_i[s1, r2]
    col_11 <- (r1 - 1L) * nfac + r2
    col_22 <- (s1 - 1L) * nfac + s2
    col_12 <- (r1 - 1L) * nfac + s2
    col_21 <- (s1 - 1L) * nfac + r2

    b_11 <- mtm[cbind(r1, r2)]
    b_22 <- mtm[cbind(s1, s2)]
    b_12 <- mtm[cbind(r1, s2)]
    b_21 <- mtm[cbind(s1, r2)]

    fs2kc <- t(t(fs2[, keep_idx, drop = FALSE]) - fs2_mean[keep_idx])
    part2 <- (fs2kc[, pa, drop = FALSE] * fs2kc[, pb, drop = FALSE] -
              lambda2_eff *
              (t(t(fs2[, col_11, drop = FALSE]) * b_22) +
               t(t(fs2[, col_22, drop = FALSE]) * b_11) +
               t(t(fs2[, col_12, drop = FALSE]) * b_21) +
               t(t(fs2[, col_21, drop = FALSE]) * b_12) +
               rep((1 - 2 * lambda1) * (b_11 * b_22 + b_21 * b_12),
                   each = n)))

    iveta2_1 <- cbind(part1, part2)
    cov_iveta2 <- cov(iveta2_1) * (n - 1) / n
  }

  # extra attributes?
  if (extra) {
    attr(veta2, "lambda") <- lambda
    attr(veta2, "alpha") <- alpha
    attr(veta2, "lambda.star") <- lambda_star
    attr(veta2, "MSM") <- var_fs2
    attr(veta2, "MTM") <- var_error
    attr(veta2, "FS.mean") <- fs_mean
  }
  if (return_fs) {
    attr(veta2, "FS") <- fs2[, lv_keep, drop = FALSE]
  }
  if (return_cov_iveta2) {
    attr(veta2, "cov.iveta2") <- cov_iveta2
  }

  veta2
}

lav_sam_eeta2 <- function(eeta = NULL, veta = NULL, lv_names = NULL,
                          lv_int_names = NULL) {
  if (length(lv_int_names) == 0L) {
    lav_msg_stop(gettext("lv.int.names is empty: no lv quadratic/interaction
                         terms are provided"))
  }

  if (is.null(lv_names)) {
    lv_names <- paste("eta", seq_len(ncol(veta)), sep = "")
  }

  nfac <- nrow(veta)
  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)
  names_1 <- c(lv_names, paste(lv_names[idx1], lv_names[idx2], sep = ":"))

  # E(\eta %x% \eta)
  eeta2 <- lav_mat_vec(veta) + eeta %x% eeta

  # add 1st order
  eeta2_aug <- c(eeta, eeta2)

  # select only what we need
  names(eeta2_aug) <- names_1
  lv_keep <- c(lv_names, lv_int_names)
  eeta2_aug <- eeta2_aug[lv_keep]

  eeta2_aug
}

# compute var(fs2) including quadratic/interaction terms
lav_sam_fs2 <- function(fs = NULL, lv_names = NULL, lv_int_names = NULL) {
  varn <- function(x, n) {
    var(x) * (n - 1) / n
  }

  if (length(lv_int_names) == 0L) {
    lav_msg_stop(gettext("lv.int.names is empty: no lv quadratic/interaction
                         terms are provided"))
  }

  if (is.null(lv_names)) {
    lv_names <- paste("eta", seq_len(ncol(fs)), sep = "")
  }

  # augment to include intercept
  fs <- cbind(1, fs)
  n <- nrow(fs)
  lv_names <- c("..int..", lv_names)
  nfac <- ncol(fs)

  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)

  names_1 <- paste(lv_names[idx1], lv_names[idx2], sep = ":")

  fs2 <- fs[, idx1] * fs[, idx2]
  var_fs2 <- varn(fs2, n)

  # select only what we need
  colnames(var_fs2) <- rownames(var_fs2) <- names_1
  lv_main <- paste(lv_names[-1], "..int..", sep = ":")
  lv_keep <- c(lv_main, lv_int_names)
  var_fs2 <- var_fs2[lv_keep, lv_keep]

  var_fs2
}

# create consistent lavaan object, based on (filled in) PT
lav_sam_step3_joint <- function(fit = NULL, pt_1 = NULL, sam_method = "local") {
  lavoptions <- fit@Options

  lavoptions_joint <- lavoptions
  lavoptions_joint$optim.method <- "none"
  lavoptions_joint$optim.parscale <- "none"
  lavoptions_joint$start <- "default"
  lavoptions_joint$optim.force.converged <- TRUE
  lavoptions_joint$check.gradient <- FALSE
  lavoptions_joint$check.start <- FALSE
  lavoptions_joint$check.post <- FALSE
  lavoptions_joint$rotation <- "none"
  lavoptions_joint$se <- "none"
  lavoptions_joint$store.vcov <- FALSE # we do this manually

  if (sam_method %in% c("local", "fsr", "cfsr")) {
    lavoptions_joint$baseline <- FALSE
    lavoptions_joint$sample.icov <- FALSE
    #lavoptions.joint$h1 <- TRUE # we need this if we re-use the sam object
    lavoptions_joint$test <- "none"
    lavoptions_joint$estimator <- "none"
  } else {
    # global: the joint computes the core (base) test(s); "yuan.chan" is not a
    # core lavaan test (it is added afterwards by lav_sam_global_test()), so we
    # strip it here and make sure the unscaled "standard" base is present.
    joint_test <- lavoptions$test[lavoptions$test != "yuan.chan"]
    if (length(joint_test) == 0L || !"standard" %in% joint_test) {
      joint_test <- unique(c("standard", joint_test))
    }
    lavoptions_joint$test <- joint_test
    lavoptions_joint$estimator <- lavoptions$estimator
  }

  # set ustart values
  pt_1$ustart <- pt_1$est # as this is used if optim.method == "none"

  joint <- lavaan::lavaan(pt_1,
    slot_options = lavoptions_joint,
    slot_sample_stats = fit@SampleStats,
    slot_data = fit@Data,
    verbose = FALSE
  )
  joint
}

# (two-step robust) fit measures for the structural part of a local SAM model.
#
# The structural model FIT.PA is fitted with se = "robust.sem" using the
# step-1 induced NACOV (Gamma.eta) of the latent (co)variances/means. Its
# Satorra-Bentler scaled chi-square and the robust fit measures derived from it
# are the correct 'two-step robust' fit measures.
#
# For a single group, and for multiple groups WITHOUT across-group constraints
# in the measurement model (Case A), the per-group Gamma.eta list passed to
# FIT.PA already yields the correct robust test, and we simply reproduce it.
# For multiple groups WITH across-group measurement constraints (Case B), the
# step-1 variability couples the groups, so the relevant Gamma has nonzero
# cross-group blocks (step1$Gamma.eta.full) that a per-group NACOV list cannot
# represent. We then recompute the Satorra-Bentler scaled chi-square ourselves
# -- for the model AND the baseline -- using that full Gamma, and read off the
# robust fit measures.
lav_sam_struc_fit <- function(fit_pa = NULL, step1 = NULL) {
  measures <- c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr")

  fallback <- function() {
    out <- try(fitMeasures(fit_pa, measures), silent = TRUE)
    if (inherits(out, "try-error")) {
      out <- "(unable to obtain fit measures)"
      names(out) <- "warning"
    }
    out
  }

  # do we have a (local) Gamma.eta and a non-saturated structural model?
  have_gamma <- !is.null(step1$Gamma.eta) &&
    length(step1$Gamma.eta) > 0L && !is.null(step1$Gamma.eta[[1]])
  df_struc <- if (length(fit_pa@test) > 0L) fit_pa@test[[1]]$df else 0L
  has_baseline <- !is.null(fit_pa@baseline) &&
    !is.null(fit_pa@baseline$partable)

  if (!have_gamma || df_struc < 1L || !has_baseline) {
    return(fallback())
  }

  out <- try(
    {
      ng <- fit_pa@Data@ngroups
      ntot <- fit_pa@SampleStats@ntotal
      fg <- unlist(fit_pa@SampleStats@nobs) / ntot
      meanstr <- fit_pa@Model@meanstructure

      # full Gamma of the stacked structural statistics (cross-group blocks
      # for Case B; block-diagonal otherwise)
      if (isTRUE(step1$caseB)) {
        gamma_all <- ntot * step1$Gamma.eta.full
      } else {
        gamma_all <- lav_mat_bdiag(lapply(
          seq_len(ng), function(g) step1$Gamma.eta[[g]] / fg[g]
        ))
      }

      # model: Satorra-Bentler scaled test using the full Gamma
      test_model <- lav_test_sb(
        lavobject = fit_pa, test = "satorra.bentler",
        gamma_full = gamma_all
      )

      # baseline: refit the independence model on the latent (co)variances,
      # then apply the same Satorra-Bentler correction with the full Gamma
      fit_base <- lavaan::lavaan(
        model = fit_pa@baseline$partable,
        sample.cov = step1$VETA,
        sample.mean = if (meanstr) step1$EETA else NULL,
        sample.nobs = as.list(unlist(fit_pa@SampleStats@nobs)),
        nacov = step1$Gamma.eta,
        slot_options = fit_pa@Options
      )
      test_base <- lav_test_sb(
        lavobject = fit_base, test = "satorra.bentler",
        gamma_full = gamma_all
      )

      # inject the corrected scaled tests, then read the robust fit measures
      fit_pa@test[["satorra.bentler"]] <- test_model[["satorra.bentler"]]
      base_idx <- which(vapply(fit_pa@baseline$test,
        function(x) x$test, character(1L)) == "satorra.bentler")
      if (length(base_idx) > 0L) {
        fit_pa@baseline$test[[base_idx[1]]] <- test_base[["satorra.bentler"]]
      }

      # use the fitMeasures() labels for the SCALED/ROBUST measures, so the
      # printed header makes clear these are the two-step corrected values
      # (chisq.scaled / pvalue.scaled / cfi.robust / rmsea.robust), exactly as
      # fitMeasures() reports them. (The naive fallback() above keeps the plain
      # chisq/pvalue/cfi/rmsea labels.)
      scaled_labels <- c(
        "chisq.scaled", "df", "pvalue.scaled",
        "cfi.robust", "rmsea.robust", "srmr"
      )
      fm <- fitMeasures(fit_pa, scaled_labels)
      fm <- as.numeric(fm)
      names(fm) <- scaled_labels
      fm
    },
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    return(fallback())
  }
  out
}

# Yuan & Chan (2002) rescaled GLOBAL test statistic for sam.method = "global".
#
# The 'joint' object already carries the base statistic T = n F[S, Sigma(theta,
# gamma.hat)] -- the discrepancy of the FULL (observed-variable) model evaluated
# at the plugged-in SAM estimates -- and the correct reference df = p* - q_total
# (all measurement AND structural parameters counted) in joint@test[[1]]. What is
# missing is the scaling correction tau.hat = tr(Q Gamma)/(p* - q) of eq. (22),
# which accounts for (a) the non-normality of the data (via the ADF Gamma) and
# (b) the fact that the measurement parameters gamma were estimated SEPARATELY in
# step 1 (segregated estimation), not jointly with theta.
#
# Key simplification that lets us reuse lav_test_sb() verbatim. With
#   H.theta = sigma.dot.theta (sigma.dot.theta' W sigma.dot.theta)^-1
#                             sigma.dot.theta' W   (the W-projector onto the
#                                                   structural columns of Delta),
#   D       = I - sigma.dot.gamma P   (P = d gamma.hat / d stats, eq. 12a),
# Yuan & Chan's U (eq. 19) gives  I - sigma.dot U = (I - H.theta) D, hence
#   Q = (I - sigma.dot U)' W (I - sigma.dot U) = D' U.struc D ,
# with U.struc = W - W Delta.theta (Delta.theta' W Delta.theta)^-1 Delta.theta' W
# -- exactly the Satorra-Bentler residual-weight matrix, but built from the
# STRUCTURAL columns of Delta only. By the cyclic trace property
#   tr(Q Gamma) = tr( U.struc . (D Gamma D') ).
# So we hand lav_test_sb():
#   - delta   = Delta.theta              (structural columns of the full Delta)
#   - e_inv   = (sum_g f_g Delta.theta_g' W_g Delta.theta_g)^-1
#   - wls_v   = W                        (the estimator weight of the full model)
#   - gamma_full = D Gamma D'            (the 'measurement-residualized' Gamma)
#   - test_unscaled = the standard (unscaled) full-model test  (T and df)
# and read off the satorra.bentler-scaled statistic, which IS eq. (22).
#
# P (the step-1 measurement influence) is the SAME jacobian already used by the
# twostep.robust SEs (lav_sam_step1_local_jac(p_only = TRUE)). Multigroup support
# (Case A block-diagonal, Case B cross-group) follows from the stacked assembly:
# D and Gamma.tilde are built in the stacked statistics space, so the cross-group
# coupling that across-group measurement constraints induce in P is carried
# through automatically.
#
# Returns a list with two components, to be stored in the joint object:
#   $test          -> joint@test          (standard + the YC-scaled test)
#   $baseline.test -> joint@baseline$test (standard + a scaled baseline test)
# The baseline (independence) model has no measurement model, so the Yuan-Chan
# correction there reduces to the ordinary Satorra-Bentler scaling; we obtain it
# by refitting the baseline with test = "satorra.bentler". Both scaled tests are
# needed so that the robust/scaled fit measures (CFI, RMSEA) can be computed. On
# any failure the original (unscaled) test list is returned unchanged.
# 'test' is the internal SB-family base used to drive lav_test_sb (eq. 22 is the
# satorra.bentler-type rescaling); the resulting entry is renamed to "yuan.chan"
# (and given a Yuan-Chan label) in both the model and the baseline test lists.
lav_sam_global_test <- function(joint = NULL, step1 = NULL, step2 = NULL,
                                fit = NULL, test = "satorra.bentler") {
  test_orig <- joint@test

  # rename a satorra.bentler entry (and its $test field) to "yuan.chan"
  rename_yc <- function(test_list) {
    idx <- which(names(test_list) == "satorra.bentler")
    if (length(idx) == 0L) {
      return(test_list)
    }
    e <- test_list[[idx[1]]]
    e$test  <- "yuan.chan"
    e$label <- "Yuan-Chan (2002) correction"
    names(test_list)[idx[1]] <- "yuan.chan"
    test_list[[idx[1]]] <- e
    test_list
  }
  fallback  <- list(test = test_orig, baseline.test = NULL)

  # locate the standard (unscaled) full-model test = base statistic T
  std_idx <- which(vapply(test_orig, function(x) x$test, character(1L)) ==
                   "standard")
  if (length(std_idx) == 0L) {
    return(fallback) # no usable base statistic (eg test = "none")
  }
  std <- test_orig[[std_idx[1]]]
  if (is.null(std$df) || anyNA(std$df) || std$df < 1L) {
    return(fallback) # saturated full model: no global test
  }

  out <- try(
    {
      ngroups <- joint@Data@ngroups
      ntot    <- joint@SampleStats@ntotal
      fg      <- unlist(joint@SampleStats@nobs) / ntot
      step1_free_idx <- step1$step1.free.idx
      step2_free_idx <- step2$step2.free.idx

      # ingredients of the full (joint) model
      if (!is.null(joint@SampleStats@NACOV[[1]])) {
        gamma <- joint@SampleStats@NACOV
      } else {
        # default: the (unbiased) ADF Gamma. This is unavailable for the
        # fixed.x / conditional.x setting (lav_samp_gamma() with unbiased =
        # TRUE stops), so fall back to the *biased* ADF Gamma there
        # (gamma.unbiased = FALSE) -- exactly the Gamma that simultaneous
        # sem(test = "satorra.bentler") uses for such models. Asymptotically
        # equivalent; the simple (no-x) cases keep the unbiased Gamma so their
        # values are unchanged.
        gamma <- tryCatch(lavTech(joint, "gamma"),
                          error = function(e) NULL)
        if (is.null(gamma)) {
          opts <- joint@Options
          opts$gamma.unbiased <- FALSE
          gamma <- lav_object_gamma(
            lavdata        = joint@Data,
            lavoptions     = opts,
            lavsamplestats = joint@SampleStats,
            lavh1          = joint@h1,
            lavimplied     = joint@implied,
            model_based    = FALSE
          )
        }
      }
      delta <- lavTech(joint, "Delta")
      wls_v <- lavTech(joint, "WLS.V")

      # P = d(gamma.mm)/d(stats), reusing the twostep.robust jacobian; align its
      # rows with step1.free.idx (P rows are in ascending free-parameter order)
      p <- lav_sam_step1_local_jac(step1 = step1, fit = fit, p_only = TRUE)
      p_row_idx <- match(step1_free_idx, attr(p, "free.idx"))
      if (anyNA(p_row_idx)) {
        lav_msg_stop(gettext(
          "internal error: unable to align the step 1 jacobian (P) for the
           global test statistic"))
      }
      p <- p[p_row_idx, , drop = FALSE]

      # per group: structural and measurement columns of Delta, the weight W,
      # and the SB-convention Gamma (Cov of sqrt(N) * per-group statistics)
      delta_theta_list <- vector("list", ngroups)
      delta_gamma_list <- vector("list", ngroups)
      gamma_sb_list    <- vector("list", ngroups)
      a_theta <- matrix(0, length(step2_free_idx), length(step2_free_idx))
      for (g in seq_len(ngroups)) {
        delta_theta_list[[g]] <- delta[[g]][, step2_free_idx, drop = FALSE]
        delta_gamma_list[[g]] <- delta[[g]][, step1_free_idx, drop = FALSE]
        gamma_sb_list[[g]]    <- gamma[[g]] / fg[g]
        # structural information E.theta = sum_g f_g Delta.theta_g' W_g Delta.theta_g
        wd <- lav_sam_wls_delta(wls_v[[g]], delta_theta_list[[g]])
        a_theta <- a_theta + fg[g] * crossprod(delta_theta_list[[g]], wd)
      }
      e_theta_inv <- solve(a_theta)

      # D (stacked) = I - Delta.gamma P ; Gamma.tilde = D Gamma D'
      delta_gamma_all <- do.call(rbind, delta_gamma_list)
      nstat_all <- nrow(delta_gamma_all)
      d_mat <- diag(nstat_all) - delta_gamma_all %*% p
      gamma_sb_all <- lav_mat_bdiag(gamma_sb_list)
      gamma_tilde  <- d_mat %*% gamma_sb_all %*% t(d_mat)

      # reuse lav_test_sb: structural-only Delta/E.inv/W + residualized Gamma
      # (gamma_full forces the assembled 'original' path)
      lav_test_sb(
        lavmodel       = joint@Model,
        lavsamplestats = joint@SampleStats,
        lavoptions     = joint@Options,
        lavimplied     = joint@implied,
        lavdata        = joint@Data,
        test_unscaled  = std,
        e_inv          = e_theta_inv,
        delta          = delta_theta_list,
        wls_v          = wls_v,
        m_gamma        = gamma, # ignored when gamma_full is set; avoids recompute
        gamma_full     = gamma_tilde,
        test           = test
      )
    },
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    # the Yuan-Chan correction reuses the same step-1 jacobian (P) as the robust
    # SEs, so it shares their current limitations: it is not available for models
    # with within-block equality constraints or conditional.x = TRUE (the P-row
    # alignment then fails). Exogenous covariates with the default fixed.x = TRUE
    # ARE supported (the biased ADF Gamma is used). Fall back to the unscaled
    # test, with a note.
    lav_msg_warn(gettext(
      "the Yuan & Chan (2002) scaled test statistic (test = \"yuan.chan\") is
       not available for this model (eg models with within-block equality
       constraints or conditional.x = TRUE); the standard (unscaled) test is
       reported instead."))
    return(fallback) # fall back to the unscaled standard test
  }

  # name the scaled entry "yuan.chan" (the SAM-global analogue of the
  # simultaneous "satorra.bentler"); @test = list(standard, yuan.chan)
  out <- rename_yc(out)

  # baseline (independence) model: no measurement model -> the Yuan-Chan
  # correction reduces to the ordinary Satorra-Bentler scaling. Refit the
  # baseline with the matching scaled test so the robust CFI/RMSEA can be
  # computed (and so fitMeasures()/summary() do not choke on a model that has a
  # scaled test while its baseline does not).
  baseline_test <- NULL
  if (!is.null(joint@baseline$partable)) {
    baseline_test <- try(
      {
        opts <- fit@Options
        opts$se       <- "none"
        opts$test     <- test
        opts$baseline <- FALSE
        opts$estimator <- joint@Model@estimator
        fit_base <- lavaan::lavaan(
          model             = joint@baseline$partable,
          slot_data         = joint@Data,
          slot_sample_stats = joint@SampleStats,
          slot_options      = opts,
          verbose           = FALSE
        )
        rename_yc(fit_base@test)
      },
      silent = TRUE
    )
    if (inherits(baseline_test, "try-error")) {
      baseline_test <- NULL
    }
  }

  list(test = out, baseline.test = baseline_test)
}

lav_sam_table <- function(joint = NULL, step1 = NULL, fit_pa = NULL,
                          cmd = NULL, lavoptions = NULL,
                          mm_args = list(), struc_args = list(),
                          sam_method = "local",
                          local_options = list(), global_options = list()) {
  mm_fit <- step1$MM.FIT

  sam_mm_table <- data.frame(
    Block = seq_along(step1$mm.list),
    Latent = sapply(mm_fit, function(x) {
      paste(unique(unlist(x@pta$vnames$lv)), collapse = ",")
    }),
    Nind = sapply(mm_fit, function(x) {
      length(unique(unlist(x@pta$vnames$ov)))
    }),
    Chisq = sapply(mm_fit, function(x) {
      x@test[[1]]$stat
    }),
    Df = sapply(mm_fit, function(x) {
      x@test[[1]]$df
    })
  )
  class(sam_mm_table) <- c("lavaan.data.frame", "data.frame")

  # extra info for @internal slot
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    # (two-step robust) fit measures for the structural part; for multigroup
    # models with across-group measurement constraints (Case B) these are
    # recomputed with the full cross-group Gamma
    sam_struc_fit <- lav_sam_struc_fit(fit_pa = fit_pa, step1 = step1)
    sam_mm_rel <- step1$REL
  } else {
    sam_struc_fit <- paste0("no local fit measures available for",
                    "structural part if sam.method is global")
    names(sam_struc_fit) <- "warning"
    sam_mm_rel <- numeric(0L)
  }


  sam_1 <- list(
    sam.cmd = cmd,
    sam.method = sam_method,
    sam.local.options = local_options,
    sam.global.options = global_options,
    sam.mm.list = step1$mm.list,
    sam.mm.estimator = mm_fit[[1]]@Model@estimator,
    sam.mm.args = mm_args,
    sam.mm.ov.names = lapply(mm_fit, function(x) {
      x@pta$vnames$ov
    }),
    sam.mm.table = sam_mm_table,
    sam.mm.rel = sam_mm_rel,
    sam.struc.estimator = fit_pa@Model@estimator,
    sam.struc.args = struc_args,
    sam.struc.fit = sam_struc_fit,
    sam.lavoptions = lavoptions
  )
  sam_1
}

lav_sam_get_cov_ybar <- function(fit = NULL, local_options = list(
                                  M.method = "ML",
                                  lambda.correction = TRUE,
                                  alpha.correction = 0L,
                                  twolevel.method = "h1"
                                )) {

  # local.twolevel.method
  local_twolevel_method <- tolower(local_options[["twolevel.method"]])
  if (!local_twolevel_method %in% c("h1", "anova", "mean")) {
    lav_msg_stop(gettext(
      "local option twolevel.method should be one of h1, anova or mean."))
  }

  lavpta <- fit@pta
  nlevels <- lavpta$nlevels
  nblocks <- lavpta$nblocks

  # do we need H1?
  if (nlevels > 1L && local_twolevel_method == "h1") {
    h1 <- lav_h1_implied_logl(
      lavdata = fit@Data,
      lavsamplestats = fit@SampleStats,
      lavoptions = fit@Options
    )
    h1implied <- h1$implied
  } else {
    h1implied <- fit@h1$implied
  }

  # containers
  cov_list  <- vector("list", nblocks)
  ybar_list <- vector("list", nblocks)

  # label
  if (nblocks > 1L) {
    names(cov_list)  <- fit@Data@block.label
    names(ybar_list) <- fit@Data@block.label
  }

  # collect COV/YBAR per block
  for (b in seq_len(nblocks)) {

    # get sample statistics for this block
    if (nlevels > 1L) {
      # block order is group-major: (g1l1, g1l2, g2l1, g2l2, ...)
      this_level <- (b - 1L) %% nlevels + 1L
      this_group <- floor((b - 1L) / nlevels) + 1L
      if (this_level > 2L) {
        lav_msg_stop(gettext("level 3 not supported (yet)."))
      }

      if (local_twolevel_method == "h1") {
        # note: the h1 estimates are already reduced to the
        # level-specific variables -- no need to subset
        cov_1 <- h1implied$cov[[(this_group - 1L) * nlevels + this_level]]
        ybar <- h1implied$mean[[(this_group - 1L) * nlevels + this_level]]
      } else {
        ylp <- fit@SampleStats@YLp[[this_group]][[2]]
        if (this_level == 1L) {
          # "anova" and "mean" use the same within statistics
          cov_1 <- ylp$Sigma.W
          ybar <- ylp$Mu.W
        } else if (local_twolevel_method == "anova") {
          cov_1 <- ylp$Sigma.B
          ybar <- ylp$Mu.B
        } else { # "mean"
          s_pw <- ylp$Sigma.W
          nj <- ylp$s
          y2 <- ylp$Y2
          # grand mean
          mu_y <- ylp$Mu.W + ylp$Mu.B
          y2c <- t(t(y2) - mu_y) # MUST be centered
          yb <- crossprod(y2c) / nrow(y2c)
          cov_1 <- yb - 1 / nj * s_pw
          ybar <- ylp$Mu.B
        }

        # reduce to the level-specific variables
        ov_idx <- fit@Data@Lp[[this_group]]$ov.idx[[this_level]]
        cov_1 <- cov_1[ov_idx, ov_idx, drop = FALSE]
        ybar <- ybar[ov_idx]
      }

      # single level
    } else {
      this_group <- b
      if (fit@Model@conditional.x) {
        ybar <- h1implied$res.int[[b]]
        cov_1  <- h1implied$res.cov[[b]]
      } else {
        ybar <- h1implied$mean[[b]] # EM version if missing="ml"
        cov_1  <- h1implied$cov[[b]]
      }
    } # single level

    cov_list[[b]] <- cov_1
    ybar_list[[b]] <- ybar
  }

  list(COV = cov_list, YBAR = ybar_list)
}

# automatically generate mm.list; we create measurement blocks so that:
# - unlinked factors become a singleton
# - linked factor are joined into a single measurement block
#   based: cross-loadings, or correlated residuals between item errors
lav_sam_get_mmlist <- function(lavobject) {

  lavmodel <- lavobject@Model
  lavpta <- lavobject@pta
  nblocks <- lavpta$nblocks

  lambda_idx <- which(names(lavmodel@GLIST) == "lambda")

  mm_list <- vector("list", length = nblocks)

  glist <- lavTech(lavobject, "partable")
  for (b in seq_len(nblocks)) {
    mm_in_block <- (seq_len(lavmodel@nmat[b]) +
        cumsum(c(0, lavmodel@nmat))[b])
    mlist <- glist[mm_in_block]
    cc <- t(mlist$lambda) %*% mlist$theta %*% mlist$lambda

    # note: CC contains dummy lv's, higher-order, etc... and they
    # must be removed

    # get ALL lv names (including dummy ov.x/ov.y)
    ov_names <- lavmodel@dimNames[[lambda_idx[b]]][[1L]]
    lv_names <- lavmodel@dimNames[[lambda_idx[b]]][[2L]]
    names_1 <- lv_names

    # needs to removed:
    # - all lv.ind variables
    # - all higher-order latent variables
    # - all dummy lv's
    rm_idx <- c(match(lavpta$vnames$lv.ind[[b]], lv_names),
                match(lavpta$vnames$lv.interaction[[b]], lv_names),
                which(lv_names %in% ov_names))
    if (length(rm_idx) > 0L) {
      cc <- cc[-rm_idx, -rm_idx, drop = FALSE]
      names_1 <- lv_names[-rm_idx]
    }

    # cluster membership
    membership <- lav_graph_get_connected_nodes(cc)

    out <- split(names_1, names_1[membership])
    names(out) <- paste("block", seq_along(out), sep = "")

    mm_list[[b]] <- out
  }

  mm_list
}

# labels for the elements of the WLS.obs statistics vector in the
# categorical, conditional.x = FALSE case: 1) thresholds and (negative)
# means, interleaved per variable, 2) variances (numeric variables only),
# 3) correlations (vech, no diagonal); used to map the statistics of a
# measurement block into the statistics vector of the joint model
lav_sam_wls_obs_labels <- function(ov_names, th_idx) {
  nvar <- length(ov_names)

  # variable of each th element: ordered variables have th_idx > 0; the
  # th_idx == 0 elements are the (negative) means of the numeric
  # variables, appearing in variable order
  v_of <- th_idx
  num_var_idx <- which(!seq_len(nvar) %in% th_idx)
  v_of[th_idx == 0L] <- num_var_idx

  th_labels <- character(length(th_idx))
  count <- integer(nvar)
  for (k in seq_along(th_idx)) {
    v <- v_of[k]
    count[v] <- count[v] + 1L
    th_labels[k] <- paste0(ov_names[v], "|t", count[v])
  }

  # variances (numeric variables only)
  var_labels <- character(0L)
  if (length(num_var_idx) > 0L) {
    var_labels <- paste0(ov_names[num_var_idx], "|var")
  }

  # correlations: lav_mat_vech(cov, diagonal = FALSE), column-major;
  # use a canonical (sorted) pair label, as the relative order of two
  # variables may differ between the joint model and a measurement block
  cor_labels <- character(0L)
  if (nvar > 1L) {
    row_idx <- matrix(seq_len(nvar), nvar, nvar)[lower.tri(diag(nvar))]
    col_idx <- matrix(seq_len(nvar), nvar, nvar,
                      byrow = TRUE)[lower.tri(diag(nvar))]
    cor_labels <- vapply(seq_along(row_idx), function(k) {
      paste(sort(c(ov_names[row_idx[k]], ov_names[col_idx[k]])),
            collapse = "~~")
    }, character(1L))
  }

  c(th_labels, var_labels, cor_labels)
}

# merge the per-block mm clusters (as returned by lav_sam_get_mmlist())
# into a single mm.list, transitively joining clusters that share a latent
# variable or an observed indicator; needed for multilevel models, where
# the within and between factors of the same indicators must end up in the
# same (two-level) measurement block
lav_sam_mmlist_merge_blocks <- function(lavobject, mm_list_per_block) {
  pt_1 <- lavobject@ParTable
  nblocks <- lavobject@pta$nblocks
  block_values <- lav_pt_block_values(pt_1)

  # collect all clusters: a set of lv names + their observed indicators
  clusters <- list()
  for (b in seq_len(nblocks)) {
    for (cl in mm_list_per_block[[b]]) {
      lv <- unlist(cl)
      ind_idx <- which(pt_1$op == "=~" & pt_1$block == block_values[b] &
                       pt_1$lhs %in% lv)
      clusters <- c(clusters,
                    list(list(lv = lv, ov = unique(pt_1$rhs[ind_idx]))))
    }
  }

  # merge until no two clusters share an lv or an ov
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    i <- 1L
    while (i < length(clusters)) {
      j <- i + 1L
      while (j <= length(clusters)) {
        if (any(clusters[[i]]$lv %in% clusters[[j]]$lv) ||
            any(clusters[[i]]$ov %in% clusters[[j]]$ov)) {
          clusters[[i]]$lv <- unique(c(clusters[[i]]$lv, clusters[[j]]$lv))
          clusters[[i]]$ov <- unique(c(clusters[[i]]$ov, clusters[[j]]$ov))
          clusters[[j]] <- NULL
          changed <- TRUE
        } else {
          j <- j + 1L
        }
      }
      i <- i + 1L
    }
  }

  lapply(clusters, "[[", "lv")
}

lav_sam_veta_pt <- function(lavobject, block = 1L) {

  lavmodel <- lavobject@Model

  lambda_idx <- which(names(lavmodel@GLIST) == "lambda")

  glist <- lavTech(lavobject, "partable")
  b <- block
  mm_in_block <- (seq_len(lavmodel@nmat[b]) + cumsum(c(0, lavmodel@nmat))[b])
  mlist <- glist[mm_in_block]
  mm_psi    <- (mlist$psi    != 0) + 0L
  nr <- nrow(mlist$psi)
  if (!is.null(mlist$beta)) {
    mm_beta <- (mlist$beta != 0) + 0L
    tmp <- -mm_beta
    tmp[lav_mat_diag_idx(nr)] <- 1
    ib_inv <- solve(tmp)
  } else {
    ib_inv <- diag(nr)
  }
  veta <- ib_inv %*% mm_psi %*% t(ib_inv)
  # get ALL lv names (including dummy ov.x/ov.y)
  lv_names <- lavmodel@dimNames[[lambda_idx[b]]][[2L]]
  colnames(veta) <- rownames(veta) <- lv_names

  veta
}

lav_sam_veta_con <- function(s = NULL, mm_lambda = NULL, mm_theta = NULL,
                             l_veta = NULL, local_m_method = "ML",
                             tol = 1e-07, max_iter = 100L) {

  # first do GLS/ULS
  m <- ncol(mm_lambda)
  if (local_m_method == "ULS") {
    m_w <- diag(m)
  } else {
    m_w <- solve(s)
  }

  wl <- m_w %*% mm_lambda
  t_lw  <- t(wl)
  t_lwl <- t_lw %*% mm_lambda
  tmp1 <- t(l_veta) %*% (t_lwl %x% t_lwl) %*% l_veta
  tmp2 <- t(l_veta) %*% lav_mat_vec(t_lw %*% (s - mm_theta) %*% wl)
  out <- solve(tmp1, tmp2)
  veta_init <- matrix(l_veta %*% out, m, m)

  if (local_m_method != "ML") {
    # we are done!
    return(veta_init)
  }

  veta_new <- veta_init
  for (i in seq_len(max_iter)) {
    sigma_ml <- mm_lambda %*% veta_new %*% t(mm_lambda) + mm_theta
    m_w <- solve(sigma_ml)
    wl <- m_w %*% mm_lambda
    t_lw  <- t(wl)
    t_lwl <- t_lw %*% mm_lambda
    tmp1 <- t(l_veta) %*% (t_lwl %x% t_lwl) %*% l_veta
    tmp2 <- t(l_veta) %*% lav_mat_vec(t_lw %*% (s - mm_theta) %*% wl)
    out <- solve(tmp1, tmp2)
    veta_ml <- matrix(l_veta %*% out, m, m)
    rmsea <- sqrt(sum((veta_new - veta_ml)^2))
    #cat("i = ", i, " rmsea = ", rmsea, "\n")
    veta_new <- veta_ml
    if (rmsea < tol) {
      break
    }
  }

  veta_new
}
