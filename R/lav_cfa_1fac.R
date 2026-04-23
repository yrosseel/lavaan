# special functions for the one-factor model
# YR 24 June 2018

# 1-factor model with (only) three indicators:
# no iterations needed; can be solved analytically

# denote s11, s22, s33 the diagonal elements, and
# s21, s31, s32 the off-diagonal elements

# under the 1-factor model; typically, either psi == 1, or l1 == 1
# - s11 == l1^2*psi + theta1
# - s22 == l2^2*psi + theta2
# - s33 == l3^2*psi + theta3
# - s21 == l2*l1*psi
# - s31 == l3*l1*psi
# - s32 == l3*l2*psi
# 6 unknowns, 6 knowns

# note: if the triad of covariances is negative, there is no
# `valid' solution, for example:
#
# > S
#      [,1] [,2] [,3]
# [1,]  1.0  0.6  0.3
# [2,]  0.6  1.0 -0.1
# [3,]  0.3 -0.1  1.0
#
# (note: all eigenvalues are positive)

lav_cfa_1fac_3ind <- function(sample_cov, std_lv = FALSE,
                              warn_neg_triad = TRUE, bounds = TRUE) {
  # check sample cov
  stopifnot(is.matrix(sample_cov))
  n_row <- NROW(sample_cov)
  n_col <- NCOL(sample_cov)
  stopifnot(n_row == n_col, n_row < 4L, n_col < 4L)
  nvar <- n_row

  # we expect a 3x3 sample covariance matrix
  # however, if we get a 2x2 (or 1x1 covariance matrix), do something
  # useful anyways...
  if (nvar == 1L) {
    # lambda = 1, theta = 0, psi = sample.cov[1,1]
    # lambda = 1, theta = 0, psi = 1 (for now, until NlsyLinks is fixed)
    sample_cov <- matrix(1, 3L, 3L) * 1.0
  } else if (nvar == 2L) {
    # hm, we could force both lambda's to be 1, but if the second
    # one is negative, this will surely lead to non-convergence issues
    #
    # just like lavaan < 0.6.2, we will use the regression of y=marker
    # on x=item2
    mean_2var <- mean(diag(sample_cov))
    max_var <- max(diag(sample_cov))
    extra <- c(mean_2var, sample_cov[2, 1])
    sample_cov <- rbind(cbind(sample_cov, extra, deparse.level = 0),
      c(extra, max_var),
      deparse.level = 0
    )
  }

  s11 <- sample_cov[1, 1]
  s22 <- sample_cov[2, 2]
  s33 <- sample_cov[3, 3]
  stopifnot(s11 > 0, s22 > 0, s33 > 0)

  s21 <- sample_cov[2, 1]
  s31 <- sample_cov[3, 1]
  s32 <- sample_cov[3, 2]
  # note: s21*s31*s32 should be positive!
  neg_triad <- FALSE
  if (s21 * s31 * s32 < 0) {
    neg_triad <- TRUE
    if (warn_neg_triad) {
      lav_msg_warn(gettext("product of the three covariances is negative!"))
    }
  }

  # first, we assume l1 = 1
  psi <- (s21 * s31) / s32 # note that we assume |s32|>0
  l1 <- 1
  l2 <- s32 / s31 # l2 <- s21/psi
  l3 <- s32 / s21 # l3 <- s31/psi
  theta1 <- s11 - psi
  theta2 <- s22 - l2 * l2 * psi
  theta3 <- s33 - l3 * l3 * psi

  # sanity check (new in 0.6-11): apply standard bounds
  if (bounds) {
    lower_psi <- s11 - (1 - 0.1) * s11 # we assume REL(y1) >= 0.1
    psi <- min(max(psi, lower_psi), s11)

    l2_bound <- sqrt(s22 / lower_psi)
    l2 <- min(max(-l2_bound, l2), l2_bound)
    l3_bound <- sqrt(s33 / lower_psi)
    l3 <- min(max(-l3_bound, l3), l3_bound)

    theta1 <- min(max(theta1, 0), s11)
    theta2 <- min(max(theta2, 0), s22)
    theta3 <- min(max(theta3, 0), s33)
  }

  lambda <- c(l1, l2, l3)
  theta <- c(theta1, theta2, theta3)


  # std.lv?
  if (std_lv) {
    # we allow for negative psi (if bounds = FALSE)
    lambda <- lambda * sign(psi) * sqrt(abs(psi))
    psi <- 1
  }

  # special cases
  if (nvar == 1L) {
    lambda <- lambda[1]
    theta <- theta[1]
  } else if (nvar == 2L) {
    lambda <- lambda[1:2]
    theta <- theta[1:2]
    psi <- psi / 2 # smaller works better?
  }

  list(lambda = lambda, theta = theta, psi = psi, neg.triad = neg_triad)
}

# FABIN (Hagglund, 1982)
# 1-factor only
lav_cfa_1fac_fabin <- function(s, lambda_only = FALSE, method = "fabin3",
                               std_lv = FALSE, bounds = TRUE) {
  # check arguments
  if (std_lv) {
    lambda_only <- FALSE # we need psi
  }

  nvar <- NCOL(s)

  # catch nvar < 4
  if (nvar < 4L) {
    out <- lav_cfa_1fac_3ind(
      sample_cov = s, std_lv = std_lv,
      warn_neg_triad = FALSE
    )
    return(out)
  }

  # 1. lambda
  lambda <- numeric(nvar)
  lambda[1L] <- 1.0
  for (i in 2:nvar) {
    idx3 <- (1:nvar)[-c(i, 1L)]
    s23 <- s[i, idx3]
    s31 <- s13 <- s[idx3, 1L]
    if (method == "fabin3") {
      s33 <- s[idx3, idx3]
      tmp <- try(solve(s33, s31), silent = TRUE) # GaussJordanPivot is
      # slighty more efficient
      if (inherits(tmp, "try-error")) {
        lambda[i] <- sum(s23 * s31) / sum(s13^2)
      } else {
        lambda[i] <- sum(s23 * tmp) / sum(s13 * tmp)
      }
    } else {
      lambda[i] <- sum(s23 * s31) / sum(s13^2)
    }
  }

  # bounds? (new in 0.6-11)
  if (bounds) {
    s11 <- s[1, 1]
    lower_psi <- s11 - (1 - 0.1) * s11 # we assume REL(y1) >= 0.1
    for (i in 2:nvar) {
      l_bound <- sqrt(s[i, i] / lower_psi)
      lambda[i] <- min(max(-l_bound, lambda[i]), l_bound)
    }
  }

  if (lambda_only) {
    return(list(
      lambda = lambda, psi = as.numeric(NA),
      theta = rep(as.numeric(NA), nvar)
    ))
  }

  # 2. theta

  # GLS version
  # W <- solve(S)
  # LAMBDA <- as.matrix(lambda)
  # A1 <- solve(t(LAMBDA) %*% W %*% LAMBDA) %*% t(LAMBDA) %*% W
  # A2 <- W %*% LAMBDA %*% A1

  # tmp1 <- W*W - A2*A2
  # tmp2 <- diag( W %*% S %*% W - A2 %*% S %*% A2 )
  # theta.diag <- solve(tmp1, tmp2)

  # 'least squares' version, assuming W = I
  d <- tcrossprod(lambda) / sum(lambda^2)
  theta <- solve(diag(nvar) - d * d, diag(s - (d %*% s %*% d)))

  # 3. psi (W=I)
  s1 <- s - diag(theta)
  l2 <- sum(lambda^2)
  psi <- sum(colSums(as.numeric(lambda) * s1) * lambda) / (l2 * l2)

  # std.lv?
  if (std_lv) {
    # we allow for negative psi
    lambda <- lambda * sign(psi) * sqrt(abs(psi))
    psi <- 1
  }

  list(lambda = lambda, theta = theta, psi = psi)
}
