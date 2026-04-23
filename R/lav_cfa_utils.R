# utility functions needed for lav_cfa_*

# compute THETA and PSI, given lambda using either ULS or GLS
# this function assumes:
# - THETA is diagonal
# - PSI is unrestricted
# - we assume W = S^{-1}
#
# YR 17 oct 2022: - add lower/upper bounds for theta (only to compute PSI)
#                 - use 'lambda' correction to ensure PSI is positive definite
# YR 02 feb 2023: - add psi.mapping.ML argument
lav_cfa_lambda2thetapsi <- function(lambda = NULL, s = NULL, s_inv = NULL,
                                    gls = FALSE, psi_mapping_ml = FALSE,
                                    nobs = 20L) {
  mm_lambda <- as.matrix(lambda)
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)

  if (gls) {
    # see Browne, 1974 section 4 case II
    if (is.null(s_inv)) {
      m_w <- solve(s)
    } else {
      m_w <- s_inv
    }
    t_lw <- crossprod(mm_lambda, m_w)
    m <- solve(t_lw %*% mm_lambda, t_lw) # GLS mapping
    # d <- W %*% LAMBDA %*% M # symmmetric
    d <- crossprod(m, t_lw)
    # theta <- solve(W*W - d*d, diag(W %*% S %*% W - d %*% S %*% d))
    theta <- try(solve(m_w * m_w - d * d, diag(m_w - d)), # because W == S^{-1}
      silent = TRUE
    )
    if (inherits(theta, "try-error")) {
      # what to do?
      lav_msg_warn(gettext(
        "problem computing THETA values; trying pace algorithm"))
      theta <- lav_efa_pace(S = s, nfactors = nfac, theta.only = TRUE)
    }
  } else {
    # see Hagglund 1982, section 4
    m <- solve(crossprod(mm_lambda), t(mm_lambda)) # ULS mapping function
    d <- mm_lambda %*% m
    theta <- try(solve(diag(nvar) - d * d, diag(s - (d %*% s %*% d))),
      silent = TRUE
    )
    if (inherits(theta, "try-error")) {
      # what to do?
      lav_msg_warn(gettext(
        "problem computing THETA values; trying pace algorithm"))
      theta <- lav_efa_pace(S = s, nfactors = nfac, theta.only = TRUE)
    }
  }
  theta_nobounds <- theta

  # ALWAYS check bounds for theta (only to to compute PSI)!
  theta_bounds <- TRUE
  if (theta_bounds) {
    diag_s <- diag(s)
    # lower bound
    lower_bound <- diag_s * 0 # * 0.01
    too_small_idx <- which(theta < lower_bound)
    if (length(too_small_idx) > 0L) {
      theta[too_small_idx] <- lower_bound[too_small_idx]
    }

    # upper bound
    upper_bound <- diag_s * 1 # * 0.99
    too_large_idx <- which(theta > upper_bound)
    if (length(too_large_idx) > 0L) {
      theta[too_large_idx] <- upper_bound[too_large_idx]
    }
  }

  # psi
  diag_theta <- diag(theta, nvar)
  lambda <- try(lav_matrix_symmetric_diff_smallest_root(s, diag_theta),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    s_min_theta <- s - diag_theta # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda_star <- lambda - 1 / (nobs - 1)
      s_min_theta <- s - lambda_star * diag_theta
    } else {
      s_min_theta <- s - diag_theta
    }
  }

  # just like local SAM
  if (psi_mapping_ml) {
    ti <- 1 / theta
    zero_theta_idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero_theta_idx) > 0L) {
      ti[zero_theta_idx] <- 1
    }
    m <- solve(t(mm_lambda) %*% diag(ti, nvar) %*% mm_lambda) %*%
      t(mm_lambda) %*% diag(ti, nvar)
    mm_psi <- m %*% s_min_theta %*% t(m) # ML
  } else {
    mm_psi <- m %*% s_min_theta %*% t(m) # ULS/GLS
  }

  # we take care of the bounds later!
  list(lambda = mm_lambda, theta = theta_nobounds, psi = mm_psi)
}

# compute PSI, given lambda and theta using either ULS, GLS, ML
# this function assumes:
# - THETA is diagonal
# - PSI is unrestricted
#
# YR 08 Mar 2023: - first version
lav_cfa_lambdatheta2psi <- function(lambda = NULL, theta = NULL, # vector!
                                    s = NULL, s_inv = NULL,
                                    mapping = "ML", nobs = 20L) {
  mm_lambda <- as.matrix(lambda)
  nvar <- nrow(mm_lambda)

  # ALWAYS check bounds for theta to compute PSI
  diag_s <- diag(s)
  # lower bound
  lower_bound <- diag_s * 0 # * 0.01
  too_small_idx <- which(theta < lower_bound)
  if (length(too_small_idx) > 0L) {
    theta[too_small_idx] <- lower_bound[too_small_idx]
  }

  # upper bound
  upper_bound <- diag_s * 1 # * 0.99
  too_large_idx <- which(theta > upper_bound)
  if (length(too_large_idx) > 0L) {
    theta[too_large_idx] <- upper_bound[too_large_idx]
  }

  # psi
  diag_theta <- diag(theta, nvar)
  lambda <- try(lav_matrix_symmetric_diff_smallest_root(s, diag_theta),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    s_min_theta <- s - diag_theta # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda_star <- lambda - 1 / (nobs - 1)
      s_min_theta <- s - lambda_star * diag_theta
    } else {
      s_min_theta <- s - diag_theta
    }
  }

  # mapping matrix
  if (mapping == "ML") {
    ti <- 1 / theta
    zero_theta_idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero_theta_idx) > 0L) {
      ti[zero_theta_idx] <- 1
    }
    m <- solve(t(mm_lambda) %*% diag(ti, nvar) %*% mm_lambda) %*%
         t(mm_lambda) %*% diag(ti, nvar)
  } else if (mapping == "GLS") {
    if (is.null(s_inv)) {
      s_inv <- try(solve(s), silent = TRUE)
    }
    if (inherits(s_inv, "try-error")) {
      m <- tcrossprod(solve(crossprod(mm_lambda)), mm_lambda)
    } else {
      m <- solve(t(mm_lambda) %*% s_inv %*% mm_lambda) %*%
           t(mm_lambda) %*% s_inv
    }
  } else if (mapping == "ULS") {
    m <- tcrossprod(solve(crossprod(mm_lambda)), mm_lambda)
  }

  # compute PSI
  mm_psi <- m %*% s_min_theta %*% t(m)

  mm_psi
}

# compute theta elements for a 1-factor model
lav_cfa_theta_spearman <- function(s, bounds = "wide") {
  p <- ncol(s)
  out <- numeric(p)
  r <- cov2cor(s)
  for (p_idx in seq_len(p)) {
    x <- r[, p_idx][-p_idx]
    aa <- lav_matrix_vech(tcrossprod(x), diagonal = FALSE)
    ss <- lav_matrix_vech(r[-p_idx, -p_idx, drop = FALSE], diagonal = FALSE)
    h2 <- mean(aa / ss) # communaliteit
    if (bounds == "standard") {
      h2[h2 < 0] <- 0
      h2[h2 > 1] <- 1
    } else if (bounds == "wide") {
      h2[h2 < -0.05] <- -0.05 # correponds to lower bound ov.var "wide"
      h2[h2 > +1.20] <- +1.20 # correponds to upper bound ov.var "wide"
    }
    out[p_idx] <- (1 - h2) * s[p_idx, p_idx]
  }
  out
}
