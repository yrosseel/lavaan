# YR 18 Dec 2015
# - functions to (directly) compute the inverse of 'Gamma' (the asymptotic
#   variance matrix of the sample statistics)
# - often used as 'WLS.V' (the weight matrix in WLS estimation)
#   and when computing the expected information matrix

# NOTE:
#  - three types:
#       1) plain         (conditional_x = FALSE, fixed_x = FALSE)
#       2) fixed_x       (conditional_x = FALSE, fixed_x = TRUE)
#       3) conditional_x (conditional_x = TRUE)
#  - if conditional_x = TRUE, we ignore fixed_x (can be TRUE or FALSE)

# GLS (and DLS with dls.a = 1): the normal-theory weight matrix for one
# group -- the inverse of the NT Gamma, in the (partial) correlation metric
# when correlation = TRUE. Shared by lav_samp_from_data() and
# lav_samp_from_moments().
lav_samp_wls_v_nt_g <- function(m_cov = NULL, m_mean = NULL, m_icov = NULL,
                                cor_idx = NULL, correlation = FALSE,
                                x_idx = integer(0L), fixed_x = FALSE,
                                conditional_x = FALSE,
                                meanstructure = FALSE) {
  if (correlation) {
    gamma_nt <- lav_samp_partial_cor_gamma_nt(
      m_cov         = m_cov,
      cor_idx       = cor_idx,
      meanstructure = meanstructure,
      fixed_x       = fixed_x,
      x_idx         = x_idx,
      conditional_x = conditional_x,
      m_mean        = m_mean
    )
    lav_mat_sym_inverse(gamma_nt)
  } else {
    lav_samp_gamma_inverse_nt(
      m_icov         = m_icov,
      m_cov          = m_cov,
      m_mean         = m_mean,
      rescale        = FALSE,
      x_idx          = x_idx,
      fixed_x        = fixed_x,
      conditional_x  = conditional_x,
      meanstructure  = meanstructure,
      slopestructure = conditional_x
    )
  }
}

# GLS: compute WLS.V %*% x -- where WLS.V is the (unconditional)
# normal-theory weight matrix lav_samp_wls_v_nt_g() would return --
# WITHOUT constructing WLS.V. The cov block of WLS.V is
#   0.5 t(D) (S.inv %x% S.inv) D
# so for each column a of x (in vech metric) we have
#   0.5 t(D) vec(S.inv A S.inv)   with A = vech.reverse(a)
# and the (optional) mean block is S.inv itself (rescaled by v11_scale
# for gls.v11.mplus). Cost is O(nvar^3) per column instead of O(nvar^4).
# note: correlation = TRUE and conditional.x are NOT supported here;
# their weight matrices are not this simple Kronecker sandwich
lav_samp_wls_v_nt_prod <- function(m_icov = NULL, x = NULL,
                                   meanstructure = FALSE,
                                   fixed_x = FALSE,
                                   x_idx = integer(0L),
                                   v11_scale = 1.0) {
  x <- as.matrix(x)
  out <- matrix(0, nrow = NROW(x), ncol = NCOL(x))
  nvar <- NROW(m_icov)
  pstar <- (nvar * (nvar + 1L)) %/% 2L

  # mean block: V11 = S.inv (zero rows/cols for the x-variables if fixed.x)
  if (meanstructure) {
    m_v11 <- m_icov * v11_scale
    if (fixed_x && length(x_idx) > 0L) {
      m_v11[x_idx, ] <- 0
      m_v11[, x_idx] <- 0
    }
    out[seq_len(nvar), ] <- m_v11 %*% x[seq_len(nvar), , drop = FALSE]
    cov_idx <- nvar + seq_len(pstar)
  } else {
    cov_idx <- seq_len(pstar)
  }

  # fixed.x: WLS.V has zero rows/cols for the x/x combinations
  # (see lav_samp_gamma_inverse_nt); emulate by zeroing the corresponding
  # elements of x (the zero columns) and of the result (the zero rows)
  zero_idx <- integer(0L)
  if (fixed_x && length(x_idx) > 0L) {
    m_tmp <- matrix(0L, nvar, nvar)
    m_tmp[lav_mat_vech_idx(nvar)] <- seq_len(pstar)
    zero_idx <- lav_mat_vech(m_tmp[x_idx, x_idx, drop = FALSE])
  }

  x_cov <- x[cov_idx, , drop = FALSE]
  if (length(zero_idx) > 0L) {
    x_cov[zero_idx, ] <- 0
  }
  for (j in seq_len(NCOL(x))) {
    m_a <- lav_mat_vech_rev(x_cov[, j])
    out[cov_idx, j] <- 0.5 * lav_mat_dup_pre(
      matrix(m_icov %*% m_a %*% m_icov, ncol = 1L))
  }
  if (length(zero_idx) > 0L) {
    out[cov_idx[zero_idx], ] <- 0
  }

  out
}

# NORMAL-THEORY
lav_samp_gamma_inverse_nt <- function(m_y = NULL,
                                             m_cov = NULL,
                                             m_icov = NULL,
                                             m_mean = NULL,
                                             rescale = TRUE,
                                             x_idx = integer(0L),
                                             fixed_x = FALSE,
                                             conditional_x = FALSE,
                                             meanstructure = FALSE,
                                             slopestructure = FALSE) {
  # check arguments
  if (length(x_idx) == 0L) {
    conditional_x <- FALSE
    fixed_x <- FALSE
  }

  if (is.null(m_icov)) {
    if (is.null(m_cov)) {
      stopifnot(!is.null(m_y))

      # coerce to matrix
      m_y <- unname(as.matrix(m_y))
      n <- nrow(m_y)
      m_cov <- cov(m_y)
      if (rescale) {
        m_cov <- m_cov * (n - 1) / n # ML version
      }
    }

    m_icov <- solve(m_cov)
  }

  # if conditional_x, we may also need m_cov and m_mean
  if (conditional_x && length(x_idx) > 0L &&
    (meanstructure || slopestructure)) {
    if (is.null(m_cov)) {
      stopifnot(!is.null(m_y))

      # coerce to matrix
      m_y <- unname(as.matrix(m_y))
      n <- nrow(m_y)
      m_cov <- cov(m_y)

      if (rescale) {
        m_cov <- m_cov * (n - 1) / n # ML version
      }
    }

    if (is.null(m_mean)) {
      stopifnot(!is.null(m_y))
      m_mean <- unname(colMeans(m_y))
    }
  }

  # rename
  s_inv <- m_icov
  m_s <- m_cov
  m_m <- m_mean

  # unconditional
  if (!conditional_x) {
    # unconditional - stochastic x
    if (!fixed_x) {
      # if (lav_use_lavaanC()) {
      #   gamma_inv <- lavaanC::m_kronecker_dup_pre_post(s_inv,
      #                                                multiplicator = 0.5)
      # } else {
        gamma_inv <- 0.5 * lav_mat_dup_pre_post(s_inv %x% s_inv)
      # }
      if (meanstructure) {
        gamma_inv <- lav_mat_bdiag(s_inv, gamma_inv)
      }

      # unconditional - fixed x
    } else {
      # handle fixed_x = TRUE
      # if (lav_use_lavaanC()) {
      #   gamma_inv <- lavaanC::m_kronecker_dup_pre_post(s_inv,
      #                                                multiplicator = 0.5)
      # } else {
        gamma_inv <- 0.5 * lav_mat_dup_pre_post(s_inv %x% s_inv)
      # }

      # zero rows/cols corresponding with x/x combinations
      nvar <- NROW(m_icov)
      pstar <- nvar * (nvar + 1) / 2
      m_m <- matrix(0, nvar, nvar)
      m_m[lav_mat_vech_idx(nvar)] <- seq_len(pstar)
      zero_idx <- lav_mat_vech(m_m[x_idx, x_idx, drop = FALSE])
      gamma_inv[zero_idx, ] <- 0
      gamma_inv[, zero_idx] <- 0

      if (meanstructure) {
        s_inv_nox <- s_inv
        s_inv_nox[x_idx, ] <- 0
        s_inv_nox[, x_idx] <- 0
        gamma_inv <- lav_mat_bdiag(s_inv_nox, gamma_inv)
      }
    }
  } else {
    # conditional_x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    s11 <- s_inv[-x_idx, -x_idx, drop = FALSE]

    # if (lav_use_lavaanC()) {
    #   gamma_inv <- lavaanC::m_kronecker_dup_pre_post(s11, multiplicator = 0.5)
    # } else {
      gamma_inv <- 0.5 * lav_mat_dup_pre_post(s11 %x% s11)
    # }

    if (meanstructure || slopestructure) {
      m_c <- m_s[x_idx, x_idx, drop = FALSE]
      m_mx <- m_m[x_idx]
      c3 <- lav_samp_gamma_c3(m_mx, m_c)
    }

    # the [int|slopes] block of Gamma is cov_ybarx %x% solve(c3) (see
    # lav_samp_gamma_nt): the statistics are in vecr(cbind(res.int,
    # res.slopes)) order (all coefficients of y1, all coefficients of
    # y2, ...), so the fast-running index is the coefficient index;
    # the inverse is then s11 %x% c3 (idem for the sub-cases)
    if (meanstructure) {
      if (slopestructure) {
        a11 <- s11 %x% c3
      } else {
        c11 <- 1 / solve(c3)[1, 1, drop = FALSE]
        a11 <- s11 %x% c11
      }
    } else {
      if (slopestructure) {
        # (solve(c3)[-1, -1])^{-1} == m_c (Schur complement)
        a11 <- s11 %x% m_c
      } else {
        a11 <- matrix(0, 0, 0)
      }
    }

    if (meanstructure || slopestructure) {
      gamma_inv <- lav_mat_bdiag(a11, gamma_inv)
    }
  }

  gamma_inv
}
