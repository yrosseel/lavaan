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

# NORMAL-THEORY
lav_samplestats_gamma_inverse_nt <- function(m_y = NULL,
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
        gamma_inv <- 0.5 * lav_matrix_duplication_pre_post(s_inv %x% s_inv)
      # }
      if (meanstructure) {
        gamma_inv <- lav_matrix_bdiag(s_inv, gamma_inv)
      }

      # unconditional - fixed x
    } else {
      # handle fixed_x = TRUE
      # if (lav_use_lavaanC()) {
      #   gamma_inv <- lavaanC::m_kronecker_dup_pre_post(s_inv,
      #                                                multiplicator = 0.5)
      # } else {
        gamma_inv <- 0.5 * lav_matrix_duplication_pre_post(s_inv %x% s_inv)
      # }

      # zero rows/cols corresponding with x/x combinations
      nvar <- NROW(m_icov)
      pstar <- nvar * (nvar + 1) / 2
      m_m <- matrix(0, nvar, nvar)
      m_m[lav_matrix_vech_idx(nvar)] <- seq_len(pstar)
      zero_idx <- lav_matrix_vech(m_m[x_idx, x_idx, drop = FALSE])
      gamma_inv[zero_idx, ] <- 0
      gamma_inv[, zero_idx] <- 0

      if (meanstructure) {
        s_inv_nox <- s_inv
        s_inv_nox[x_idx, ] <- 0
        s_inv_nox[, x_idx] <- 0
        gamma_inv <- lav_matrix_bdiag(s_inv_nox, gamma_inv)
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
      gamma_inv <- 0.5 * lav_matrix_duplication_pre_post(s11 %x% s11)
    # }

    if (meanstructure || slopestructure) {
      m_c <- m_s[x_idx, x_idx, drop = FALSE]
      m_mx <- m_m[x_idx]
      c3 <- rbind(
        c(1, m_mx),
        cbind(m_mx, m_c + tcrossprod(m_mx))
      )
    }

    if (meanstructure) {
      if (slopestructure) {
        a11 <- c3 %x% s11
      } else {
        c11 <- 1 / solve(c3)[1, 1, drop = FALSE]
        a11 <- c11 %x% s11
      }
    } else {
      if (slopestructure) {
        a11 <- m_c %x% s11
      } else {
        a11 <- matrix(0, 0, 0)
      }
    }

    if (meanstructure || slopestructure) {
      gamma_inv <- lav_matrix_bdiag(a11, gamma_inv)
    }
  }

  gamma_inv
}
