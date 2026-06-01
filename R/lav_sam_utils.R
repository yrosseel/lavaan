# utility functions for the sam() function
# YR 4 April 2023

# construct 'mapping matrix' M using either "ML", "GLS" or "ULS" method
# optionally return MTM (for ML)
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
    # M == solve( t(LAMBDA) %*% LAMBDA ) %*% t(LAMBDA)
    #   == MASS:::ginv(LAMBDA)
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
        if (inherits(mm_theta, "try-error")) {
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

  # empty theta elements?
  empty_theta_idx <- which(diag(mtm) == 0)

  # new in 0.6-16: make sure MTM is pd
  # (otherwise lav_mat_sym_diff_smallest_root will fail)
  theta_rm_idx <- unique(c(dummy_lv_idx, empty_theta_idx))
  if (length(theta_rm_idx) == nrow(mtm)) {
    # all zero?
    # do nothing, but certainly no need for alpha or lambda.correction
    alpha_correction <- 0L
    lambda_correction <- FALSE
  } else if (length(theta_rm_idx) > 0L) {
    mtm_small <- mtm[-theta_rm_idx, -theta_rm_idx, drop = FALSE]
    mtm_small <- zapsmall(lav_mat_sym_force_pd(
      mtm_small,
      tol = 1e-04
    ))
    mtm[-theta_rm_idx, -theta_rm_idx] <- mtm_small
  } else {
    mtm <- zapsmall(lav_mat_sym_force_pd(mtm, tol = 1e-04))
  }

  # apply small sample correction (if requested)
  if (alpha_correction > 0) {
    alpha_n1 <- alpha_correction / (n - 1)
    if (alpha_n1 > 1.0) {
      alpha_n1 <- 1.0
    } else if (alpha_n1 < 0.0) {
      alpha_n1 <- 0.0
    }
    mtm <- (1 - alpha_n1) * mtm
    alpha <- alpha_correction
  } else {
    alpha <- alpha_correction
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
  eeta <- m %*% (ybar - mm_nu)
  eeta
}

# compute veta including quadratic/interaction terms
lav_sam_veta2 <- function(fs = NULL, m = NULL,
                          veta = NULL, eeta = NULL, mm_theta = NULL,
                          lv_names = NULL,
                          lv_int_names = NULL,
                          dummy_lv_names = character(0L),
                          alpha_correction = 0L,
                          lambda_correction = TRUE,
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

  # new in 0.6-16: make sure MTM is pd
  # (otherwise lav_mat_sym_diff_smallest_root will fail)
  dummy_lv_idx <- which(lv_names %in% dummy_lv_names)
  if (length(dummy_lv_idx) > 0L) {
    mtm_nodummy <- mtm[-dummy_lv_idx, -dummy_lv_idx, drop = FALSE]
    mtm_nodummy <- zapsmall(lav_mat_sym_force_pd(
      mtm_nodummy,
      tol = 1e-04
    ))
    mtm[-dummy_lv_idx, -dummy_lv_idx] <- mtm_nodummy
  } else {
    mtm <- zapsmall(lav_mat_sym_force_pd(mtm, tol = 1e-04))
  }

  # augment to include intercept
  fs <- cbind(1, fs)
  n <- nrow(fs)
  mtm <- lav_mat_bdiag(0, mtm)
  veta <- lav_mat_bdiag(0, veta)
  eeta <- c(1, eeta)
  lv_names <- c("..int..", lv_names)
  nfac <- ncol(fs)

  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)

  names_1 <- paste(lv_names[idx1], lv_names[idx2], sep = ":")
  names_1[seq_len(nfac)] <- lv_names

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

  # compute Gamma for FS2[,lv.keep]
  #FS.gamma <- lav_samp_gamma(FS2[,lv.keep, drop  = FALSE],
  #                                  meanstructure = TRUE)

  # apply small sample correction (if requested)
  if (alpha_correction > 0) {
    alpha_n1 <- alpha_correction / (n - 1)
    if (alpha_n1 > 1.0) {
      alpha_n1 <- 1.0
    } else if (alpha_n1 < 0.0) {
      alpha_n1 <- 0.0
    }
    var_error <- (1 - alpha_n1) * var_error
    alpha <- alpha_correction
  } else {
    alpha <- alpha_correction
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
      #cutoff <- 1 + 1 / (N - 1)
    cutoff <- 1 + 2 / n # be more conservative for VETA2
      if (lambda < cutoff) {
        #lambda.star <- lambda - 1 / (N - 1)
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
  if (return_cov_iveta2) {
    iveta2_1 <- matrix(0, n,
                ncol = length(lav_mat_vech(veta2)) + ncol(veta2))
    ef <- colMeans(fs)
    for (i in 1:n) {
      fi <- as.matrix(fs2[i, seq_len(nfac)])
      tmp <- (((tcrossprod(fi) - mtm) %x% mtm) +
           (mtm %x% (tcrossprod(fi) - mtm)) +
           lav_mat_com_post((tcrossprod(fi) - mtm) %x% mtm) +
           lav_mat_com_pre((tcrossprod(fi) - mtm) %x% mtm) +
           (ik %*% (mtm %x% mtm)))
      iveta2 <- tcrossprod(fs2[i, ] - fs2_mean) - lambda_star * tmp
      colnames(iveta2) <- rownames(iveta2) <- names_1

      ieeta2 <- (lav_mat_vec(tcrossprod(fi)) -
                 lav_mat_vec(tcrossprod(ef)) +
                 (ef %x% ef) -
                 lambda_star * lav_mat_vec(mtm))
      names(ieeta2) <- names_1

      iveta2_1[i, ] <- c(ieeta2[lv_keep],
                      lav_mat_vech(iveta2[lv_keep, lv_keep]))
    } # N
    # experimental:
    # remove outliers in FS?
    #if (length(fs.outlier.idx) > 0L) {
    #  IVETA2 <- IVETA2[-fs.outlier.idx, ,drop = FALSE]
    #  N <- N - length(fs.outlier.idx)
    #}
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
    #attr(VETA2, "FS.gamma") <- FS.gamma
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
    lavoptions_joint$test <- lavoptions$test
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
    # Estimator = sapply(MM.FIT, function(x) { x@Model@estimator} ),
    Chisq = sapply(mm_fit, function(x) {
      x@test[[1]]$stat
    }),
    Df = sapply(mm_fit, function(x) {
      x@test[[1]]$df
    })
  )
  # pvalue = sapply(MM.FIT, function(x) {x@test[[1]]$pvalue}) )
  class(sam_mm_table) <- c("lavaan.data.frame", "data.frame")


  # extra info for @internal slot
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    sam_struc_fit <- try(
      fitMeasures(
        fit_pa,
        c(
          "chisq", "df", # "pvalue",
          "cfi", "rmsea", "srmr"
        )
      ),
      silent = TRUE
    )
    if (inherits(sam_struc_fit, "try-error")) {
      sam_struc_fit <- "(unable to obtain fit measures)"
      names(sam_struc_fit) <- "warning"
    }
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

  # local_m_method <- toupper(local_options[["M.method"]])

  lavpta <- fit@pta
  ngroups <- lavpta$ngroups
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
    # if (FIT@Options$conditional.x) {
    #   h1implied <- lav_model_implied_cond2uncond(h1implied)
    # }
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
      if (ngroups > 1L) {
        this_level <- (b - 1L) %% ngroups + 1L
      } else {
        this_level <- b
      }
      this_group <- floor(b / nlevels + 0.5)

      if (this_level == 1L) {
        if (local_twolevel_method == "h1") {
          cov_1 <- h1implied$cov[[1]]
          ybar <- h1implied$mean[[1]]
        } else if (local_twolevel_method == "anova" ||
          local_twolevel_method == "mean") {
          cov_1 <- fit@SampleStats@YLp[[this_group]][[2]]$Sigma.W
          ybar <- fit@SampleStats@YLp[[this_group]][[2]]$Mu.W
        }

        # reduce
        ov_idx <- fit@Data@Lp[[this_group]]$ov.idx[[this_level]]
        cov_1 <- cov_1[ov_idx, ov_idx, drop = FALSE]
        ybar <- ybar[ov_idx]
      } else if (this_level == 2L) {
        if (local_twolevel_method == "h1") {
          cov_1 <- h1implied$cov[[2]]
          ybar <- h1implied$mean[[2]]
        } else if (local_twolevel_method == "anova") {
          cov_1 <- fit@SampleStats@YLp[[this_group]][[2]]$Sigma.B
          ybar <- fit@SampleStats@YLp[[this_group]][[2]]$Mu.B
        } else if (local_twolevel_method == "mean") {
          s_pw <- fit@SampleStats@YLp[[this_group]][[2]]$Sigma.W
          nj <- fit@SampleStats@YLp[[this_group]][[2]]$s
          y2 <- fit@SampleStats@YLp[[this_group]][[2]]$Y2
          # grand mean
          mu_y <- (fit@SampleStats@YLp[[this_group]][[2]]$Mu.W +
                   fit@SampleStats@YLp[[this_group]][[2]]$Mu.B)
          y2c <- t(t(y2) - mu_y) # MUST be centered
          yb <- crossprod(y2c) / nrow(y2c)
          cov_1 <- yb - 1 / nj * s_pw
          ybar <- fit@SampleStats@YLp[[this_group]][[2]]$Mu.B
        }

        # reduce
        ov_idx <- fit@Data@Lp[[this_group]]$ov.idx[[this_level]]
        cov_1 <- cov_1[ov_idx, ov_idx, drop = FALSE]
        ybar <- ybar[ov_idx]
      } else {
        lav_msg_stop(gettext("level 3 not supported (yet)."))
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

  # flags
  # lv_interaction_flag <- FALSE
  # lv_higherorder_flag <- FALSE
  # if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
  #   lv_interaction_flag <- TRUE
  # }
  # if (length(unlist(lavpta$vnames$lv.ind)) > 0L) {
  #   lv_higherorder_flag <- TRUE
  # }

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

lav_sam_veta_pt <- function(lavobject, block = 1L) {

  lavmodel <- lavobject@Model
  lavpta   <- lavobject@pta
  # nblocks  <- lavpta$nblocks

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
