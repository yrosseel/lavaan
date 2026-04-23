# James-Stein estimator
#
# Burghgraeve, E., De Neve, J., & Rosseel, Y. (2021). Estimating structural
# equation models using James-Stein type shrinkage estimators. Psychometrika,
# 86(1), 96-130.
#
# YR 08 Feb 2023: - first version in lavaan, cfa only (for now)

lav_cfa_jamesstein <- function(s,
                               y = NULL, # raw data
                               marker_idx = NULL,
                               lambda_nonzero_idx = NULL,
                               theta = NULL, # vector!
                               theta_bounds = TRUE,
                               aggregated = FALSE) { # aggregated?
  # dimensions
  nvar <- ncol(s)
  nfac <- length(marker_idx)
  stopifnot(length(theta) == nvar)
  stopifnot(ncol(y) == nvar)

  # overview of lambda structure
  m_b <- mm_lambda <- b_nomarker <- matrix(0, nvar, nfac)
  lambda_marker_idx <- (seq_len(nfac) - 1L) * nvar + marker_idx
  m_b[lambda_marker_idx] <- mm_lambda[lambda_marker_idx] <- 1L
  m_b[lambda_nonzero_idx] <- b_nomarker[lambda_nonzero_idx] <- 1L

  # Nu
  mm_nu <- numeric(nvar)

  # do we first 'clip' the theta values so they are within standard bounds?
  # (Question: do we need the 0.01 and 0.99 multipliers?)
  diag_s <- diag(s)
  if (theta_bounds) {
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

  # compute conditional expectation conditional on the scaling indicator
  e_js1 <- lav_cfa_jamesstein_ce(
    y = y, marker_idx = marker_idx,
    resvars_markers = theta[marker_idx]
  )

  # compute LAMBDA
  for (f in seq_len(nfac)) {
    nomarker_idx <- which(b_nomarker[, f] == 1)
    y_nomarker_f <- y[, nomarker_idx, drop = FALSE]

    # regress no.marker.idx data on E(\eta|Y)
    fit <- lm(y_nomarker_f ~ e_js1[, f, drop = FALSE])

    # extract 'lambda' values
    mm_lambda[nomarker_idx, f] <- drop(coef(fit)[-1, ])

    # (optional) extract means
    # NU[nomarker.idx] <- drop(coef(fit)[1,])

    if (aggregated) {
      # local copy of 'scaling' LAMBDA
      lambda_scaling <- mm_lambda

      j_1 <- length(nomarker_idx)
      for (j in seq_len(j_1)) {
        # data without this indicator
        j_idx <- nomarker_idx[j]
        no_j_idx <- c(marker_idx[f], nomarker_idx[-j])
        y_agg <- y[, no_j_idx, drop = FALSE]
        y_j <- y[, j_idx, drop = FALSE]

        # retrieve estimated values scaling JS
        lambda_js_scaling <- lambda_scaling[no_j_idx, f, drop = FALSE]

        # optimize the weights
        starting_weights <- rep(1 / j_1, times = j_1)
        w <- optim(
          par = starting_weights,
          fn = lav_cfa_jamesstein_rel,
          data = y_agg,
          resvars = theta[no_j_idx]
        )$par

        # make sure the weights sum up to 1
        w_optim <- w / sum(w)

        # compute aggregated indicator using the optimal weights
        y_agg_1 <- t(t(w_optim) %*% t(y_agg))

        # compute error variance of the aggregated indicator
        var_eps_agg <- drop(t(w_optim) %*%
          diag(theta[no_j_idx], nrow = length(no_j_idx)) %*% w_optim)

        # compute conditional expectation using aggregated indicator
        tmp <- lav_cfa_jamesstein_ce(
          y = y_agg_1, marker_idx = 1L,
          resvars_markers = var_eps_agg
        )
        ce_agg <- tmp / drop(w_optim %*% lambda_js_scaling)

        # compute factor loading
        fit <- lm(y_j ~ ce_agg)
        mm_lambda[j_idx, f] <- drop(coef(fit)[-1])

        # (optional) extract means
        # NU[j.idx] <- drop(coef(fit)[1,])
      } # j
    } # aggregate
  } # f

  list(lambda = mm_lambda, nu = mm_nu)
}

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_cfa_jamesstein_internal <- function(lavobject = NULL, # convenience
                                        # internal slot
                                        lavmodel = NULL,
                                        lavsamplestats = NULL,
                                        lavpartable = NULL,
                                        lavdata = NULL,
                                        lavoptions = NULL,
                                        theta_bounds = TRUE) {
  lavpta <- NULL
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))

    # extract slots
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavpartable <- lav_partable_set_cache(lavobject@ParTable, lavobject@pta)
    lavpta <- lavobject@pta
    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
  }
  if (is.null(lavpta)) {
    lavpta <- lav_partable_attributes(lavpartable)
    lavpartable <- lav_partable_set_cache(lavpartable, lavpta)
  }

  # no structural part!
  if (any(lavpartable$op == "~")) {
    lav_msg_stop(gettext(
      "JS(A) estimator only available for CFA models (for now)"))
  }
  # no BETA matrix! (i.e., no higher-order factors)
  if (!is.null(lavmodel@GLIST$beta)) {
    lav_msg_stop(gettext(
      "JS(A) estimator not available for models that require a BETA matrix"))
  }
  # no std.lv = TRUE for now
  if (lavoptions$std.lv) {
    lav_msg_stop(gettext("S(A) estimator not available if std.lv = TRUE"))
  }

  nblocks <- lav_partable_nblocks(lavpartable)
  stopifnot(nblocks == 1L) # for now
  b <- 1L
  sample_cov <- lavsamplestats@cov[[b]]
  nvar <- nrow(sample_cov)
  lv_names <- lavpta$vnames$lv.regular[[b]]
  nfac <- length(lv_names)
  marker_idx <- lavpta$vidx$lv.marker[[b]]
  lambda_idx <- which(names(lavmodel@GLIST) == "lambda")
  lambda_nonzero_idx <- lavmodel@m.free.idx[[lambda_idx]]
  # only diagonal THETA for now...
  theta_idx <- which(names(lavmodel@GLIST) == "theta") # usually '2'
  m_theta <- lavmodel@m.free.idx[[theta_idx]]
  nondiag_idx <- m_theta[!m_theta %in% lav_matrix_diag_idx(nvar)]
  if (length(nondiag_idx) > 0L) {
    lav_msg_warn(gettext(
      "this implementation of JS/JSA does not handle correlated residuals yet!"
      ))
  }


  # 1. obtain estimate for (diagonal elements of) THETA
  #    for now we use Spearman per factor
  m_b <- matrix(0, nvar, nfac)
  lambda_marker_idx <- (seq_len(nfac) - 1L) * nvar + marker_idx
  m_b[lambda_marker_idx] <- 1L
  m_b[lambda_nonzero_idx] <- 1L
  theta <- numeric(nvar)
  for (f in seq_len(nfac)) {
    ov_idx <- which(m_b[, f] == 1L)
    s_fac <- sample_cov[ov_idx, ov_idx, drop = FALSE]
    theta[ov_idx] <- lav_cfa_theta_spearman(s_fac, bounds = "wide")
  }
  mm_theta <- diag(theta, nrow = nvar)

  # 2. run James-Stein algorithm
  y <- lavdata@X[[1]] # raw data
  aggregated <- FALSE
  if (lavoptions$estimator == "JSA") {
    aggregated <- TRUE
  }
  out <- lav_cfa_jamesstein(
    s = sample_cov, y = y, marker_idx = marker_idx,
    lambda_nonzero_idx = lambda_nonzero_idx,
    theta = theta,
    # experimental
    theta_bounds = theta_bounds,
    #
    aggregated = aggregated
  )
  mm_lambda <- out$lambda

  # 3. PSI
  mm_psi <- lav_cfa_lambdatheta2psi(
    lambda = mm_lambda, theta = theta,
    s = sample_cov, mapping = "ML"
  )

  # store matrices in lavmodel@GLIST
  lavmodel@GLIST$lambda <- mm_lambda
  lavmodel@GLIST$theta <- mm_theta
  lavmodel@GLIST$psi <- mm_psi

  # extract free parameters only
  x <- lav_model_get_parameters(lavmodel)

  # apply bounds (if any)
  if (!is.null(lavpartable$lower)) {
    lower_x <- lavpartable$lower[lavpartable$free > 0]
    too_small_idx <- which(x < lower_x)
    if (length(too_small_idx) > 0L) {
      x[too_small_idx] <- lower_x[too_small_idx]
    }
  }
  if (!is.null(lavpartable$upper)) {
    upper_x <- lavpartable$upper[lavpartable$free > 0]
    too_large_idx <- which(x > upper_x)
    if (length(too_large_idx) > 0L) {
      x[too_large_idx] <- upper_x[too_large_idx]
    }
  }

  x
}


# Conditional expectation (Section 2.1, eq. 10)
lav_cfa_jamesstein_ce <- function(y = NULL,
                                  marker_idx = NULL,
                                  resvars_markers = NULL) {
  y <- as.matrix(y)

  # sample size
  n <- nrow(y)
  n1 <- n - 1
  n3 <- n - 3

  # markers only
  y_marker <- y[, marker_idx, drop = FALSE]

  # means and variances
  mean_1 <- colMeans(y_marker, na.rm = TRUE)
  var_1 <- apply(y_marker, 2, var, na.rm = TRUE)

  # 1 - R per maker
  onemin_r <- n3 * resvars_markers / (n1 * var_1)

  # R per marker
  r <- 1 - onemin_r

  # create E(\eta | Y)
  e_eta_cond_y <- t(t(y_marker) * r + onemin_r * mean_1)

  e_eta_cond_y
}

# Reliability function used to obtain the weights (Section 4, Aggregation)
lav_cfa_jamesstein_rel <- function(w = NULL, data = NULL, resvars = NULL) {
  # construct weight vector
  w <- matrix(w, ncol = 1)

  # construct aggregated indicator: y_agg = t(w) %*% y_i
  y_agg <- t(t(w) %*% t(data))

  # calculate variance of aggregated indicator
  var_y_agg <- var(y_agg)

  # calculate error variance of the aggregated indicator
  var_eps_agg <- t(w) %*% diag(resvars) %*% w

  # reliability function to be maximized
  rel <- (var_y_agg - var_eps_agg) %*% solve(var_y_agg)

  # return value
  -rel
}
