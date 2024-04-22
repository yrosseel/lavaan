# 'robust' mean and (co)variance matrix using Huber weights
#
# see Yuan & Hayashi (2010). Fitting Data to Model: SEM Diagnosis using two
#     scatter plots. Psychological Methods, 15(4), 335-351
#
# this function is based on the 'robmusig' function from K.H. Yuan's website:
# https://www.nd.edu/~kyuan/SEMdiagnosis
# see file CFA.r lines 46--96
#
lav_cov_huber <- function(Y = NULL, prob = 0.95, max.it = 200L, tol = 1e-07) {
  Y <- as.matrix(Y)
  NAMES <- colnames(Y)
  Y <- unname(Y)
  N <- nrow(Y)
  P <- ncol(Y)

  # tuning parameters for Huber's weight
  chip <- qchisq(prob, P)
  ck <- sqrt(chip)
  cbeta <- (P * pchisq(chip, P + 2L) + chip * (1 - prob)) / P

  # initial values
  this.mu <- colMeans(Y, na.rm = TRUE)
  this.sigma <- cov(Y, use = "pairwise.complete.obs")

  for (i in seq_len(max.it)) {
    # store old
    old.mu <- this.mu
    old.sigma <- this.sigma

    # squared Mahalanobis distance
    inv.sigma <- solve(this.sigma)
    Y.c <- t(t(Y) - this.mu)
    mdist2 <- rowSums((Y.c %*% inv.sigma) * Y.c)
    mdist <- sqrt(mdist2)

    # Huber weights
    wt <- ifelse(mdist <= ck, 1, ck / mdist)

    # weighted mean
    this.mu <- apply(Y, 2L, weighted.mean, w = wt, na.rm = TRUE)

    # weighted cov
    Y.c <- t(t(Y) - this.mu)
    this.sigma <- crossprod(Y.c * wt) / (N * cbeta)
    # question: why N, and not sum(wt)?

    # check progress
    diff.mu <- abs(this.mu - old.mu)
    diff.sigma <- abs(this.sigma - old.sigma)
    crit <- max(c(max(diff.mu), max(diff.sigma)))
    if (crit < tol) {
      break
    }
    if (i == max.it) {
      lav_msg_warn(gettext(
        "maximum number of iterations has been reached, without convergence."))
    }
  }

  names(this.mu) <- NAMES
  colnames(this.sigma) <- rownames(this.sigma) <- NAMES

  res <- list(Mu = this.mu, Sigma = this.sigma, niter = i, wt = wt)

  res
}
