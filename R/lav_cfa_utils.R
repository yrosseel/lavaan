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
lav_cfa_lambda2thetapsi <- function(lambda = NULL, S = NULL, S.inv = NULL,
                                    GLS = FALSE, psi.mapping.ML = FALSE,
                                    nobs = 20L) {
  LAMBDA <- as.matrix(lambda)
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

  if (GLS) {
    # see Browne, 1974 section 4 case II
    if (is.null(S.inv)) {
      W <- solve(S)
    } else {
      W <- S.inv
    }
    tLW <- crossprod(LAMBDA, W)
    M <- solve(tLW %*% LAMBDA, tLW) # GLS mapping
    # D <- W %*% LAMBDA %*% M # symmmetric
    D <- crossprod(M, tLW)
    # theta <- solve(W*W - D*D, diag(W %*% S %*% W - D %*% S %*% D))
    theta <- try(solve(W * W - D * D, diag(W - D)), # because W == S^{-1}
      silent = TRUE
    )
    if (inherits(theta, "try-error")) {
      # what to do?
      lav_msg_warn(gettext(
        "problem computing THETA values; trying pace algorithm"))
      theta <- lav_efa_pace(S = S, nfactors = nfac, theta.only = TRUE)
    }
  } else {
    # see Hagglund 1982, section 4
    M <- solve(crossprod(LAMBDA), t(LAMBDA)) # ULS mapping function
    D <- LAMBDA %*% M
    theta <- try(solve(diag(nvar) - D * D, diag(S - (D %*% S %*% D))),
      silent = TRUE
    )
    if (inherits(theta, "try-error")) {
      # what to do?
      lav_msg_warn(gettext(
        "problem computing THETA values; trying pace algorithm"))
      theta <- lav_efa_pace(S = S, nfactors = nfac, theta.only = TRUE)
    }
  }
  theta.nobounds <- theta

  # ALWAYS check bounds for theta (only to to compute PSI)!
  theta.bounds <- TRUE
  if (theta.bounds) {
    diagS <- diag(S)
    # lower bound
    lower.bound <- diagS * 0 # * 0.01
    too.small.idx <- which(theta < lower.bound)
    if (length(too.small.idx) > 0L) {
      theta[too.small.idx] <- lower.bound[too.small.idx]
    }

    # upper bound
    upper.bound <- diagS * 1 # * 0.99
    too.large.idx <- which(theta > upper.bound)
    if (length(too.large.idx) > 0L) {
      theta[too.large.idx] <- upper.bound[too.large.idx]
    }
  }

  # psi
  diag.theta <- diag(theta, nvar)
  lambda <- try(lav_matrix_symmetric_diff_smallest_root(S, diag.theta),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    SminTheta <- S - diag.theta # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda.star <- lambda - 1 / (nobs - 1)
      SminTheta <- S - lambda.star * diag.theta
    } else {
      SminTheta <- S - diag.theta
    }
  }

  # just like local SAM
  if (psi.mapping.ML) {
    Ti <- 1 / theta
    zero.theta.idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero.theta.idx) > 0L) {
      Ti[zero.theta.idx] <- 1
    }
    M <- solve(t(LAMBDA) %*% diag(Ti, nvar) %*% LAMBDA) %*% t(LAMBDA) %*% diag(Ti, nvar)
    PSI <- M %*% SminTheta %*% t(M) # ML
  } else {
    PSI <- M %*% SminTheta %*% t(M) # ULS/GLS
  }

  # we take care of the bounds later!
  list(lambda = LAMBDA, theta = theta.nobounds, psi = PSI)
}

# compute PSI, given lambda and theta using either ULS, GLS, ML
# this function assumes:
# - THETA is diagonal
# - PSI is unrestricted
#
# YR 08 Mar 2023: - first version
lav_cfa_lambdatheta2psi <- function(lambda = NULL, theta = NULL, # vector!
                                    S = NULL, S.inv = NULL,
                                    mapping = "ML", nobs = 20L) {
  LAMBDA <- as.matrix(lambda)
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

  theta.nobounds <- theta

  # ALWAYS check bounds for theta to compute PSI
  diagS <- diag(S)
  # lower bound
  lower.bound <- diagS * 0 # * 0.01
  too.small.idx <- which(theta < lower.bound)
  if (length(too.small.idx) > 0L) {
    theta[too.small.idx] <- lower.bound[too.small.idx]
  }

  # upper bound
  upper.bound <- diagS * 1 # * 0.99
  too.large.idx <- which(theta > upper.bound)
  if (length(too.large.idx) > 0L) {
    theta[too.large.idx] <- upper.bound[too.large.idx]
  }

  # psi
  diag.theta <- diag(theta, nvar)
  lambda <- try(lav_matrix_symmetric_diff_smallest_root(S, diag.theta),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    SminTheta <- S - diag.theta # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda.star <- lambda - 1 / (nobs - 1)
      SminTheta <- S - lambda.star * diag.theta
    } else {
      SminTheta <- S - diag.theta
    }
  }

  # mapping matrix
  if (mapping == "ML") {
    Ti <- 1 / theta
    zero.theta.idx <- which(abs(theta) < 0.01) # be conservative
    if (length(zero.theta.idx) > 0L) {
      Ti[zero.theta.idx] <- 1
    }
    M <- solve(t(LAMBDA) %*% diag(Ti, nvar) %*% LAMBDA) %*% t(LAMBDA) %*% diag(Ti, nvar)
  } else if (mapping == "GLS") {
    if (is.null(S.inv)) {
      S.inv <- try(solve(S), silent = TRUE)
    }
    if (inherits(S.inv, "try-error")) {
      M <- tcrossprod(solve(crossprod(LAMBDA)), LAMBDA)
    } else {
      M <- solve(t(LAMBDA) %*% S.inv %*% LAMBDA) %*% t(LAMBDA) %*% S.inv
    }
  } else if (mapping == "ULS") {
    M <- tcrossprod(solve(crossprod(LAMBDA)), LAMBDA)
  }

  # compute PSI
  PSI <- M %*% SminTheta %*% t(M)

  PSI
}

# compute theta elements for a 1-factor model
lav_cfa_theta_spearman <- function(S, bounds = "wide") {
  p <- ncol(S)
  out <- numeric(p)
  R <- cov2cor(S)
  for (p.idx in seq_len(p)) {
    var.p <- R[p.idx, p.idx]
    x <- R[, p.idx][-p.idx]
    aa <- lav_matrix_vech(tcrossprod(x), diagonal = FALSE)
    ss <- lav_matrix_vech(R[-p.idx, -p.idx, drop = FALSE], diagonal = FALSE)
    h2 <- mean(aa / ss) # communaliteit
    if (bounds == "standard") {
      h2[h2 < 0] <- 0
      h2[h2 > 1] <- 1
    } else if (bounds == "wide") {
      h2[h2 < -0.05] <- -0.05 # correponds to lower bound ov.var "wide"
      h2[h2 > +1.20] <- +1.20 # correponds to upper bound ov.var "wide"
    }
    out[p.idx] <- (1 - h2) * S[p.idx, p.idx]
  }
  out
}
