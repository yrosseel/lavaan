# casewise Mahalanobis-distance diagnostics for models with ordered
# categorical (and possibly mixed continuous/ordinal) indicators, following
# Mansolf & Reise (2017)
#
# the M-distances defined for the continuous case (the leverage distance of
# the factor scores, the residual-based distance, and the distance of the
# fitted values) are generalized to categorical data by taking their
# expectation over the region of the (latent) multivariate normal response
# vector y* that is consistent with the observed response pattern:
#
#     d*_i = E[ d(y*) | y* in R_i ],   with  y* ~ N(mu, Sigma(theta))
#
# where R_i is the rectangle defined by the model-implied thresholds and the
# observed categories of case i; observed continuous coordinates fix the
# corresponding elements of y* (so the expectation is conditional on them),
# and missing coordinates leave the corresponding region unbounded (they are
# integrated out); as a result, when ALL variables are continuous and
# observed, the region degenerates to a single point and the distances
# reduce exactly to their continuous counterparts
#
# the expectation is approximated by Monte Carlo integration, using a
# (vectorized) GHK importance sampler: no burn-in or mixing issues (as with
# a Gibbs sampler), and a deterministic computational cost
#
# YR + Claude, July 2026

# main entry: compute the expected (squared) M-distances for all cases
#
# type = "lv"    -> leverage distance of the (Bartlett/regression) factor
#                   scores (d_f* in Mansolf & Reise, 2017)
# type = "resid" -> residual-based distance (d_r*), based on the casewise
#                   residuals of the latent responses
# type = "yhat"  -> distance of the fitted (model-predicted) values
#
# returns a list with one element per group; each element is a list with
#  - d2:    the expected squared distances E[d^2] (length n)
#  - d:     the expected distances E[d] (length n); note that for the
#           Monte Carlo based values E[d] != sqrt(E[d^2])
#  - resid: (only if type = "resid") the n x p matrix of expected casewise
#           residuals of the latent responses
lav_predict_mdist_cat <- function(lavobject = NULL, # for convenience
                                  # sub objects
                                  lavmodel = NULL, lavdata = NULL,
                                  lavimplied = NULL,
                                  # optional new data
                                  data_obs = NULL, exo = NULL,
                                  # options
                                  type = "lv", method = "ML",
                                  ndraws = 2000L) {
  # full object?
  if (inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
  } else {
    stopifnot(!is.null(lavmodel), !is.null(lavdata))
  }

  if (lavdata@nlevels > 1L) {
    lav_msg_stop(gettext(
      "casewise Mahalanobis distances are not available (yet) for
      multilevel data with categorical variables."))
  }

  type <- tolower(type)
  stopifnot(type %in% c("lv", "resid", "yhat"))

  # Bartlett (ML) or regression (EB/EBM) factor-score mapping?
  bartlett <- tolower(method) %in% c("ml", "bartlett", "bartlet")

  ndraws <- as.integer(ndraws)
  if (length(ndraws) != 1L || is.na(ndraws) || ndraws < 10L) {
    lav_msg_stop(gettext("mdist.draws must be a single integer >= 10."))
  }

  # data
  if (is.null(data_obs)) {
    data_obs <- lavdata@X
  }
  if (is.null(exo)) {
    exo <- lavdata@eXo
  }

  # model-implied moments for the latent responses y*
  # (on the 'marginal' scale, i.e., including the delta scaling, so that
  #  they match the sample polychoric correlations and thresholds)
  if (is.null(lavimplied) || length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel)
  }
  cond_x <- lavmodel@conditional.x
  if (cond_x) {
    sigma_hat <- lavimplied$res.cov
    int_hat <- lavimplied$res.int
    slopes_hat <- lavimplied$res.slopes
    th_hat <- lavimplied$res.th
  } else {
    sigma_hat <- lavimplied$cov
    int_hat <- lavimplied$mean
    th_hat <- lavimplied$th
  }

  mm_lambda <- lav_model_lambda(lavmodel = lavmodel,
    remove_dummy_lv = FALSE, use_wmat = TRUE)
  if (cond_x) {
    # all moments are conditional on the exogenous covariates, so the
    # latent (co)variances must be too: V(eta|x), not V(eta)
    # (this only affects the regression-based factor-score mapping; the
    #  Bartlett mapping does not involve the latent covariance matrix)
    veta <- lav_model_vetax(lavmodel = lavmodel)
  } else {
    veta <- lav_model_veta(lavmodel = lavmodel, remove_dummy_lv = FALSE)
  }
  mm_idx <- lav_model_group_mm_indices(lavmodel@nmat)

  out <- vector("list", length = lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
    p <- ncol(sigma_hat[[g]])
    n <- nrow(data_obs[[g]])

    # delta-scale the factor loadings, so that
    # Sigma = LAMBDA VETA LAMBDA' + DELTA THETA DELTA
    # with Sigma the (delta-scaled) model-implied covariance matrix of y*
    mlist <- lavmodel@GLIST[mm_idx[[g]]]
    lambda_g <- mm_lambda[[g]]
    if (!is.null(mlist$delta)) {
      lambda_g <- mlist$delta[, 1L] * lambda_g
    }

    # the 'regular' latent variables (drop the dummy lv's, if any)
    lv_dummy_idx <- c(
      lavmodel@ov.y.dummy.lv.idx[[g]],
      lavmodel@ov.x.dummy.lv.idx[[g]]
    )
    reg_lv_idx <- seq_len(ncol(lambda_g))
    if (length(lv_dummy_idx) > 0L) {
      reg_lv_idx <- reg_lv_idx[-lv_dummy_idx]
    }

    # model-implied mean of y*: a single vector, or (if conditional.x)
    # an n x p matrix with case-specific conditional means
    if (cond_x && !is.null(exo[[g]]) && ncol(exo[[g]]) > 0L) {
      mu_g <- t(as.numeric(int_hat[[g]]) +
                slopes_hat[[g]] %*% t(exo[[g]]))
    } else if (length(int_hat[[g]]) > 0L) {
      mu_g <- as.numeric(int_hat[[g]])
    } else {
      mu_g <- numeric(p)
    }

    out[[g]] <- lav_predict_mdist_cat_group(
      x = data_obs[[g]], sigma = sigma_hat[[g]], mu = mu_g,
      lambda = lambda_g, veta = veta[[g]],
      th = th_hat[[g]], th_idx = lavmodel@th.idx[[g]],
      num_idx = lavmodel@num.idx[[g]], reg_lv_idx = reg_lv_idx,
      type = type, bartlett = bartlett, ndraws = ndraws
    )
  } # g

  out
}

# single group
lav_predict_mdist_cat_group <- function(x = NULL, sigma = NULL, mu = NULL,
                                        lambda = NULL, veta = NULL,
                                        th = NULL, th_idx = NULL,
                                        num_idx = NULL, reg_lv_idx = NULL,
                                        type = "lv", bartlett = TRUE,
                                        ndraws = 2000L) {
  n <- nrow(x)
  p <- ncol(sigma)
  resid_flag <- (type == "resid")

  # drop factors without any (direct) observed indicator (e.g. higher-order
  # factors): their column in LAMBDA is empty, so they do not span any
  # dimension of y*, and keeping them would miscount the rank of the
  # residual space below
  zero_col <- which(apply(lambda, 2L, function(v) all(v == 0)))
  if (length(zero_col) > 0L) {
    keep <- seq_len(ncol(lambda))[-zero_col]
    lambda <- lambda[, keep, drop = FALSE]
    veta <- veta[keep, keep, drop = FALSE]
    reg_lv_idx <- match(intersect(reg_lv_idx, keep), keep)
  }
  nfac <- ncol(lambda)

  sigma_inv <- lav_predict_solve(sigma)

  # factor-score coefficient matrix (Bartlett or regression), mapping the
  # centered latent responses to (centered) factor scores
  if (bartlett) {
    fsc <- lav_predict_solve(t(lambda) %*% sigma_inv %*% lambda) %*%
      t(lambda) %*% sigma_inv
  } else {
    fsc <- veta %*% t(lambda) %*% sigma_inv
  }

  # for each distance type, reduce the computation to a single quadratic
  # form: d^2(y*_c) = y*_c' Q y*_c, with y*_c the centered latent response
  # vector; this mirrors the (continuous) computations in
  # lav_predict_internal(), but pushed through to the y* metric
  proj <- NULL
  if (type == "lv") {
    fsc_r <- fsc[reg_lv_idx, , drop = FALSE]
    if (nrow(fsc_r) == 0L) {
      # no regular latent variables
      return(list(d2 = rep(as.numeric(NA), n),
                  d = rep(as.numeric(NA), n), resid = NULL))
    }
    fs_cov <- fsc_r %*% sigma %*% t(fsc_r)
    qmat <- t(fsc_r) %*% solve(fs_cov) %*% fsc_r
  } else if (type == "resid") {
    proj <- diag(p) - lambda %*% fsc
    omega_e <- proj %*% sigma %*% t(proj)
    # omega_e is rank-deficient (rank p - nfac); project on the range space
    # (Yuan & Zhong, 2008)
    eig <- eigen(omega_e, symmetric = TRUE)
    a <- eig$vectors[, seq_len(p - nfac), drop = FALSE]
    qmat <- t(proj) %*% a %*% solve(t(a) %*% omega_e %*% a) %*%
      t(a) %*% proj
  } else { # yhat
    pmat <- lambda %*% fsc
    omega_h <- pmat %*% sigma %*% t(pmat)
    eig <- eigen(omega_h, symmetric = TRUE)
    a <- eig$vectors[, seq_len(nfac), drop = FALSE]
    qmat <- t(pmat) %*% a %*% solve(t(a) %*% omega_h %*% a) %*%
      t(a) %*% pmat
  }

  # per-case integration region: for the ordinal variables, the rectangle
  # implied by the thresholds and the observed category; observed continuous
  # variables are 'fixed'; missing values are left unbounded
  lower <- matrix(-Inf, nrow = n, ncol = p)
  upper <- matrix(+Inf, nrow = n, ncol = p)
  fixed <- matrix(FALSE, nrow = n, ncol = p)
  ord_idx <- unique(th_idx[th_idx > 0L])
  for (v in ord_idx) {
    tv <- th[th_idx == v]
    lo_v <- c(-Inf, tv)
    up_v <- c(tv, +Inf)
    k <- x[, v]
    ok <- which(!is.na(k))
    lower[ok, v] <- lo_v[k[ok]]
    upper[ok, v] <- up_v[k[ok]]
  }
  for (v in num_idx) {
    fixed[, v] <- !is.na(x[, v])
  }

  # center everything at the model-implied mean
  if (is.matrix(mu)) {
    mu_case <- mu
  } else {
    mu_case <- matrix(mu, nrow = n, ncol = p, byrow = TRUE)
  }
  lower_c <- lower - mu_case
  upper_c <- upper - mu_case
  xc <- x - mu_case

  # the integral is identical for cases sharing the same response pattern
  # (and, if conditional.x, the same covariate values): compute once per
  # unique combination; key on the exact (round-tripping) hexadecimal
  # representation so that only bitwise-identical rows are merged
  key_mat <- cbind(x, mu_case)
  case_key <- do.call(paste, c(
    lapply(seq_len(ncol(key_mat)), function(j) sprintf("%a", key_mat[, j])),
    sep = "\r"
  ))
  rep_idx <- which(!duplicated(case_key))
  u_of_case <- match(case_key, case_key[rep_idx])

  # the Cholesky factor of (the permuted) Sigma only depends on which
  # coordinates are fixed; cache it per missingness pattern of the
  # continuous variables
  chol_cache <- new.env(parent = emptyenv())

  n_u <- length(rep_idx)
  d2_u <- rep(as.numeric(NA), n_u)
  d_u <- rep(as.numeric(NA), n_u)
  if (resid_flag) {
    resid_u <- matrix(as.numeric(NA), nrow = n_u, ncol = p)
  }

  for (u in seq_len(n_u)) {
    i <- rep_idx[u]

    # empty cases: no information at all
    if (all(is.na(x[i, ]))) {
      next
    }

    fix_i <- which(fixed[i, ])
    samp_i <- seq_len(p)
    if (length(fix_i) > 0L) {
      samp_i <- samp_i[-fix_i]
    }
    nf <- length(fix_i)
    ns <- length(samp_i)
    # variable order: fixed (observed continuous) coordinates first, so
    # that the GHK recursion below samples exactly from the conditional
    # distribution given the observed continuous values (their importance
    # weight contribution is constant across draws, and cancels)
    ord <- c(fix_i, samp_i)

    ckey <- paste0("f", paste(fix_i, collapse = ","))
    lmat <- chol_cache[[ckey]]
    if (is.null(lmat)) {
      lmat <- try(t(chol(sigma[ord, ord, drop = FALSE])), silent = TRUE)
      if (inherits(lmat, "try-error")) {
        lav_msg_stop(gettext(
          "the model-implied covariance matrix of the latent responses is
          not positive definite; Mahalanobis distances cannot be
          computed."))
      }
      chol_cache[[ckey]] <- lmat
    }

    # degenerate region: all coordinates observed and continuous; the
    # distances reduce to their exact continuous counterparts
    if (ns == 0L) {
      zc <- xc[i, ]
      d2_u[u] <- sum((zc %*% qmat) * zc)
      d_u[u] <- sqrt(d2_u[u])
      if (resid_flag) {
        resid_u[u, ] <- drop(proj %*% zc)
      }
      next
    }

    lo_i <- lower_c[i, ord]
    up_i <- upper_c[i, ord]

    # GHK sampler: y*_c = L eta with eta iid standard normal; sample the
    # eta's recursively from their (univariate) truncated conditionals, and
    # accumulate the importance (log-)weights
    emat <- matrix(0, nrow = ndraws, ncol = p)
    if (nf > 0L) {
      eta_f <- forwardsolve(lmat[seq_len(nf), seq_len(nf), drop = FALSE],
                            xc[i, fix_i])
      emat[, seq_len(nf)] <- matrix(eta_f, nrow = ndraws, ncol = nf,
                                    byrow = TRUE)
    }
    logw <- numeric(ndraws)
    for (j in (nf + 1L):p) {
      if (j > 1L) {
        m <- drop(emat[, seq_len(j - 1L), drop = FALSE] %*%
                  lmat[j, seq_len(j - 1L)])
      } else {
        m <- numeric(ndraws)
      }
      ljj <- lmat[j, j]
      lo_j <- (lo_i[j] - m) / ljj
      up_j <- (up_i[j] - m) / ljj
      plo <- pnorm(lo_j)
      pup <- pnorm(up_j)
      logw <- logw + log(pup - plo)
      u01 <- plo + runif(ndraws) * (pup - plo)
      # guard against numerical under/overflow in the extreme tails: keep
      # u01 strictly inside (0, 1), and the draw inside its bounds
      u01 <- pmin(pmax(u01, 1e-12), 1 - 1e-12)
      emat[, j] <- pmin(pmax(qnorm(u01), lo_j), up_j)
    }

    # normalized importance weights (log-scale to avoid underflow when the
    # response pattern has a very small probability)
    mw <- max(logw)
    if (is.finite(mw)) {
      w <- exp(logw - mw)
    } else {
      # the pattern probability underflows completely; fall back to equal
      # weights (the pattern is extremely outlying anyway)
      w <- rep(1, ndraws)
    }
    wsum <- sum(w)

    # draws of the centered y*, back-permuted to the original order
    zmat_ord <- tcrossprod(emat, lmat)
    zmat <- zmat_ord
    zmat[, ord] <- zmat_ord

    dsq <- rowSums((zmat %*% qmat) * zmat)
    d2_u[u] <- sum(w * dsq) / wsum
    d_u[u] <- sum(w * sqrt(dsq)) / wsum
    if (resid_flag) {
      resid_u[u, ] <- colSums(w * tcrossprod(zmat, proj)) / wsum
    }
  } # u

  # map the per-pattern values back to the cases
  ret <- list(
    d2 = d2_u[u_of_case],
    d = d_u[u_of_case],
    resid = NULL
  )
  if (resid_flag) {
    ret$resid <- resid_u[u_of_case, , drop = FALSE]
  }

  ret
}
