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

lav_cfa_1fac_3ind <- function(sample.cov, std.lv = FALSE,
                              warn.neg.triad = TRUE, bounds = TRUE) {
  # check sample cov
  stopifnot(is.matrix(sample.cov))
  nRow <- NROW(sample.cov)
  nCol <- NCOL(sample.cov)
  stopifnot(nRow == nCol, nRow < 4L, nCol < 4L)
  nvar <- nRow

  # we expect a 3x3 sample covariance matrix
  # however, if we get a 2x2 (or 1x1 covariance matrix), do something
  # useful anyways...
  if (nvar == 1L) {
    # lambda = 1, theta = 0, psi = sample.cov[1,1]
    # lambda = 1, theta = 0, psi = 1 (for now, until NlsyLinks is fixed)
    sample.cov <- matrix(1, 3L, 3L) * 1.0
  } else if (nvar == 2L) {
    # hm, we could force both lambda's to be 1, but if the second
    # one is negative, this will surely lead to non-convergence issues
    #
    # just like lavaan < 0.6.2, we will use the regression of y=marker
    # on x=item2
    mean.2var <- mean(diag(sample.cov))
    max.var <- max(diag(sample.cov))
    extra <- c(mean.2var, sample.cov[2, 1])
    sample.cov <- rbind(cbind(sample.cov, extra, deparse.level = 0),
      c(extra, max.var),
      deparse.level = 0
    )
  }

  s11 <- sample.cov[1, 1]
  s22 <- sample.cov[2, 2]
  s33 <- sample.cov[3, 3]
  stopifnot(s11 > 0, s22 > 0, s33 > 0)

  s21 <- sample.cov[2, 1]
  s31 <- sample.cov[3, 1]
  s32 <- sample.cov[3, 2]
  # note: s21*s31*s32 should be positive!
  neg.triad <- FALSE
  if (s21 * s31 * s32 < 0) {
    neg.triad <- TRUE
    if (warn.neg.triad) {
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
    lower.psi <- s11 - (1 - 0.1) * s11 # we assume REL(y1) >= 0.1
    psi <- min(max(psi, lower.psi), s11)

    l2.bound <- sqrt(s22 / lower.psi)
    l2 <- min(max(-l2.bound, l2), l2.bound)
    l3.bound <- sqrt(s33 / lower.psi)
    l3 <- min(max(-l3.bound, l3), l3.bound)

    theta1 <- min(max(theta1, 0), s11)
    theta2 <- min(max(theta2, 0), s22)
    theta3 <- min(max(theta3, 0), s33)
  }

  lambda <- c(l1, l2, l3)
  theta <- c(theta1, theta2, theta3)


  # std.lv?
  if (std.lv) {
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

  list(lambda = lambda, theta = theta, psi = psi, neg.triad = neg.triad)
}

# FABIN (Hagglund, 1982)
# 1-factor only
lav_cfa_1fac_fabin <- function(S, lambda.only = FALSE, method = "fabin3",
                               std.lv = FALSE, bounds = TRUE) {
  # check arguments
  if (std.lv) {
    lambda.only <- FALSE # we need psi
  }

  nvar <- NCOL(S)

  # catch nvar < 4
  if (nvar < 4L) {
    out <- lav_cfa_1fac_3ind(
      sample.cov = S, std.lv = std.lv,
      warn.neg.triad = FALSE
    )
    return(out)
  }

  # 1. lambda
  lambda <- numeric(nvar)
  lambda[1L] <- 1.0
  for (i in 2:nvar) {
    idx3 <- (1:nvar)[-c(i, 1L)]
    s23 <- S[i, idx3]
    S31 <- S13 <- S[idx3, 1L]
    if (method == "fabin3") {
      S33 <- S[idx3, idx3]
      tmp <- try(solve(S33, S31), silent = TRUE) # GaussJordanPivot is
      # slighty more efficient
      if (inherits(tmp, "try-error")) {
        lambda[i] <- sum(s23 * S31) / sum(S13^2)
      } else {
        lambda[i] <- sum(s23 * tmp) / sum(S13 * tmp)
      }
    } else {
      lambda[i] <- sum(s23 * S31) / sum(S13^2)
    }
  }

  # bounds? (new in 0.6-11)
  if (bounds) {
    s11 <- S[1, 1]
    lower.psi <- s11 - (1 - 0.1) * s11 # we assume REL(y1) >= 0.1
    for (i in 2:nvar) {
      l.bound <- sqrt(S[i, i] / lower.psi)
      lambda[i] <- min(max(-l.bound, lambda[i]), l.bound)
    }
  }

  if (lambda.only) {
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
  D <- tcrossprod(lambda) / sum(lambda^2)
  theta <- solve(diag(nvar) - D * D, diag(S - (D %*% S %*% D)))

  # 3. psi (W=I)
  S1 <- S - diag(theta)
  l2 <- sum(lambda^2)
  psi <- sum(colSums(as.numeric(lambda) * S1) * lambda) / (l2 * l2)

  # std.lv?
  if (std.lv) {
    # we allow for negative psi
    lambda <- lambda * sign(psi) * sqrt(abs(psi))
    psi <- 1
  }

  list(lambda = lambda, theta = theta, psi = psi)
}
