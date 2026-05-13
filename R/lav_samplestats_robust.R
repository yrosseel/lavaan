# 'robust' mean and (co)variance matrix using Huber weights
#
# see Yuan & Hayashi (2010). Fitting Data to Model: SEM Diagnosis using two
#     scatter plots. Psychological Methods, 15(4), 335-351
#
# this function is based on the 'robmusig' function from K.H. Yuan's website:
# https://www.nd.edu/~kyuan/SEMdiagnosis
# see file CFA.r lines 46--96
#
lav_cov_huber <- function(y = NULL, prob = 0.95, max_it = 200L, tol = 1e-07) {
  y <- as.matrix(y)
  names_1 <- colnames(y)
  y <- unname(y)
  n <- nrow(y)
  p <- ncol(y)

  # tuning parameters for Huber's weight
  chip <- qchisq(prob, p)
  ck <- sqrt(chip)
  cbeta <- (p * pchisq(chip, p + 2L) + chip * (1 - prob)) / p

  # initial values
  this_mu <- colMeans(y, na.rm = TRUE)
  this_sigma <- cov(y, use = "pairwise.complete.obs")

  for (i in seq_len(max_it)) {
    # store old
    old_mu <- this_mu
    old_sigma <- this_sigma

    # squared Mahalanobis distance
    inv_sigma <- solve(this_sigma)
    y_c <- t(t(y) - this_mu)
    mdist2 <- rowSums((y_c %*% inv_sigma) * y_c)
    mdist <- sqrt(mdist2)

    # Huber weights
    wt <- ifelse(mdist <= ck, 1, ck / mdist)

    # weighted mean
    this_mu <- apply(y, 2L, weighted.mean, w = wt, na.rm = TRUE)

    # weighted cov
    y_c <- t(t(y) - this_mu)
    this_sigma <- crossprod(y_c * wt) / (n * cbeta)
    # question: why N, and not sum(wt)?

    # check progress
    diff_mu <- abs(this_mu - old_mu)
    diff_sigma <- abs(this_sigma - old_sigma)
    crit <- max(c(max(diff_mu), max(diff_sigma)))
    if (crit < tol) {
      break
    }
    if (i == max_it) {
      lav_msg_warn(gettext(
        "maximum number of iterations has been reached, without convergence."))
    }
  }

  names(this_mu) <- names_1
  colnames(this_sigma) <- rownames(this_sigma) <- names_1

  res <- list(Mu = this_mu, Sigma = this_sigma, niter = i, wt = wt)

  res
}
