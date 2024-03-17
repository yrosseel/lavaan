# YR 18 Dec 2015
# - functions to (directly) compute the inverse of 'Gamma' (the asymptotic
#   variance matrix of the sample statistics)
# - often used as 'WLS.V' (the weight matrix in WLS estimation)
#   and when computing the expected information matrix

# NOTE:
#  - three types:
#       1) plain         (conditional.x = FALSE, fixed.x = FALSE)
#       2) fixed.x       (conditional.x = FALSE, fixed.x = TRUE)
#       3) conditional.x (conditional.x = TRUE)
#  - if conditional.x = TRUE, we ignore fixed.x (can be TRUE or FALSE)

# NORMAL-THEORY
lav_samplestats_Gamma_inverse_NT <- function(Y = NULL,
                                             COV = NULL,
                                             ICOV = NULL,
                                             MEAN = NULL,
                                             rescale = TRUE,
                                             x.idx = integer(0L),
                                             fixed.x = FALSE,
                                             conditional.x = FALSE,
                                             meanstructure = FALSE,
                                             slopestructure = FALSE) {
  # check arguments
  if (length(x.idx) == 0L) {
    conditional.x <- FALSE
    fixed.x <- FALSE
  }

  if (is.null(ICOV)) {
    if (is.null(COV)) {
      stopifnot(!is.null(Y))

      # coerce to matrix
      Y <- unname(as.matrix(Y))
      N <- nrow(Y)
      COV <- cov(Y)
      if (rescale) {
        COV <- COV * (N - 1) / N # ML version
      }
    }

    ICOV <- solve(COV)
  }

  # if conditional.x, we may also need COV and MEAN
  if (conditional.x && length(x.idx) > 0L &&
    (meanstructure || slopestructure)) {
    if (is.null(COV)) {
      stopifnot(!is.null(Y))

      # coerce to matrix
      Y <- unname(as.matrix(Y))
      N <- nrow(Y)
      COV <- cov(Y)

      if (rescale) {
        COV <- COV * (N - 1) / N # ML version
      }
    }

    if (is.null(MEAN)) {
      stopifnot(!is.null(Y))
      MEAN <- unname(colMeans(Y))
    }
  }

  # rename
  S.inv <- ICOV
  S <- COV
  M <- MEAN

  # unconditional
  if (!conditional.x) {
    # unconditional - stochastic x
    if (!fixed.x) {
      Gamma.inv <- 0.5 * lav_matrix_duplication_pre_post(S.inv %x% S.inv)
      if (meanstructure) {
        Gamma.inv <- lav_matrix_bdiag(S.inv, Gamma.inv)
      }

      # unconditional - fixed x
    } else {
      # handle fixed.x = TRUE
      Gamma.inv <- 0.5 * lav_matrix_duplication_pre_post(S.inv %x% S.inv)

      # zero rows/cols corresponding with x/x combinations
      nvar <- NROW(ICOV)
      pstar <- nvar * (nvar + 1) / 2
      M <- matrix(0, nvar, nvar)
      M[lav_matrix_vech_idx(nvar)] <- seq_len(pstar)
      zero.idx <- lav_matrix_vech(M[x.idx, x.idx, drop = FALSE])
      Gamma.inv[zero.idx, ] <- 0
      Gamma.inv[, zero.idx] <- 0

      if (meanstructure) {
        S.inv.nox <- S.inv
        S.inv.nox[x.idx, ] <- 0
        S.inv.nox[, x.idx] <- 0
        Gamma.inv <- lav_matrix_bdiag(S.inv.nox, Gamma.inv)
      }
    }
  } else {
    # conditional.x

    # 4 possibilities:
    # - no meanstructure, no slopes
    # -    meanstructure, no slopes
    # - no meanstructure, slopes
    # -    meanstructure, slopes

    S11 <- S.inv[-x.idx, -x.idx, drop = FALSE]

    Gamma.inv <- 0.5 * lav_matrix_duplication_pre_post(S11 %x% S11)

    if (meanstructure || slopestructure) {
      C <- S[x.idx, x.idx, drop = FALSE]
      MY <- M[-x.idx]
      MX <- M[x.idx]
      C3 <- rbind(
        c(1, MX),
        cbind(MX, C + tcrossprod(MX))
      )
    }

    if (meanstructure) {
      if (slopestructure) {
        A11 <- C3 %x% S11
      } else {
        c11 <- 1 / solve(C3)[1, 1, drop = FALSE]
        A11 <- c11 %x% S11
      }
    } else {
      if (slopestructure) {
        A11 <- C %x% S11
      } else {
        A11 <- matrix(0, 0, 0)
      }
    }

    if (meanstructure || slopestructure) {
      Gamma.inv <- lav_matrix_bdiag(A11, Gamma.inv)
    }
  }

  Gamma.inv
}
