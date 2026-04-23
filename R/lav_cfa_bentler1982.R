# a partial implementation of the Bentler (1982) non-iterative method for CFA
#
# Bentler, P. M. (1982). Confirmatory factor-analysis via noniterative
# estimation - a fast, inexpensive method. Journal of Marketing Research,
# 19(4), 417-424. https://doi.org/10.1177/002224378201900403
#
#
# YR 03 Feb 2023: - first version in lavaan: simple setting only,
#                   no constraints, no 'fixed' (but nonzero) values,
#                   no correlated residuals (ie diagonal-theta only!)
# YR 23 Apr 2023: - quadprog is not needed if we have no (in)equality
#                   constraints

lav_cfa_bentler1982 <- function(s,
                                marker_idx = NULL,
                                lambda_nonzero_idx = NULL,
                                gls = FALSE,
                                bounds = TRUE,
                                min_reliability_marker = 0.1,
                                quadprog = FALSE,
                                nobs = 20L) { # for cutoff
  # dimensions
  nvar <- ncol(s)
  nfac <- length(marker_idx)

  # lambda structure
  m_b <- matrix(0, nvar, nfac)
  lambda_marker_idx <- (seq_len(nfac) - 1L) * nvar + marker_idx
  m_b[lambda_marker_idx] <- 1L
  m_b[lambda_nonzero_idx] <- 1L

  # partition sample covariance matrix: marker vs non-marker
  s_xx <- s[marker_idx, marker_idx, drop = FALSE]
  s_yx <- s[-marker_idx, marker_idx, drop = FALSE]
  s_yy <- s[-marker_idx, -marker_idx, drop = FALSE]
  p <- nvar - nfac
  b_y <- m_b[-marker_idx, , drop = FALSE]

  # check for p = 0?

  # phase 1: initial estimate for Sigma.yx

  # phase 2: using GLS/ULS to obtain PSI and Theta
  if (gls) {
    m_w <- try(solve(s_yy), silent = TRUE)
    if (inherits(m_w, "try-error")) {
      lav_msg_warn(gettext("could not inverte S.yy; switching to ULS"))
      m_w <- diag(p)
    }
    ws_yx <- m_w %*% s_yx
    xy_sws_yx <- crossprod(s_yx, ws_yx)
    g <- ws_yx %*% solve(xy_sws_yx) %*% t(ws_yx)
  } else {
    ip <- diag(p)
    xy_ss_yx <- crossprod(s_yx)
    g <- s_yx %*% solve(xy_ss_yx) %*% t(s_yx)
  }

  # only needed if theta.y is not diagonal:
  # q <- 6 # all free
  # # dimension P: q x p*p where q is the number of free elements theta.y
  # theta.fy <- function(x) {
  #     theta.y <- matrix(0, p, p)
  #     # insert 'free' parameters only
  #     diag(theta.y) <- x
  #     lav_matrix_vec(theta.y)
  # }
  # P <- t(numDeriv::jacobian(func = theta.fy, x = rep(1, q)))

  # tmp1 <- P %*% ((W %x% W) - (G %x% G)) %*% t(P)
  # NOTE:
  # if only the 'diagonal' element of Theta are free (as usual), then we
  # can write tmp1 as
  if (gls) {
    tmp1 <- m_w * m_w - g * g
  } else {
    tmp1 <- ip - g * g
  }

  # only needed if fixed values
  # Theta.F <- matrix(0, p, p) # all free
  # tmp2 <- W %*% (S.yy - Theta.F) %*% W - G %*% (S.yy - Theta.F) %*% G
  if (gls) {
    tmp2 <- m_w %*% s_yy %*% m_w - g %*% s_yy %*% g
  } else {
    tmp2 <- s_yy - g %*% s_yy %*% g
  }

  # Theta.f    <- as.numeric(solve(tmp1) %*% P %*% lav_matrix_vec(tmp2))
  # Note:
  # if only the 'diagonal' element of Theta are free (as usual), then we
  # can write Theta.f as
  theta_f <- solve(tmp1, diag(tmp2))
  theta_f_nobounds <- theta_f # store unbounded Theta.f values

  # ALWAYS apply standard bounds to proceed
  too_small_idx <- which(theta_f < 0)
  if (length(too_small_idx) > 0L) {
    theta_f[too_small_idx] <- 0
  }
  too_large_idx <- which(theta_f > diag(s_yy))
  if (length(too_large_idx) > 0L) {
    theta_f[too_large_idx] <- diag(s_yy)[too_large_idx] * 1
  }

  # create diagonal matrix with Theta.f elements on diagonal
  theta_yhat <- diag(theta_f, p)

  # force (S.yy - Theta.yhat) to be positive definite
  lambda <- try(lav_matrix_symmetric_diff_smallest_root(s_yy, theta_yhat),
    silent = TRUE
  )
  if (inherits(lambda, "try-error")) {
    lav_msg_warn(gettext("failed to compute lambda"))
    s_min_theta <- s_yy - theta_yhat # and hope for the best
  } else {
    cutoff <- 1 + 1 / (nobs - 1)
    if (lambda < cutoff) {
      lambda_star <- lambda - 1 / (nobs - 1)
      s_min_theta <- s_yy - lambda_star * theta_yhat
    } else {
      s_min_theta <- s_yy - theta_yhat
    }
  }

  # estimate Phi
  if (gls) {
    tmp1 <- xy_sws_yx
    tmp2 <- t(ws_yx) %*% s_min_theta %*% ws_yx
  } else {
    tmp1 <- xy_ss_yx
    tmp2 <- t(s_yx) %*% s_min_theta %*% s_yx
  }
  mm_psi <- tmp1 %*% solve(tmp2, tmp1)
  psi_nobounds <- mm_psi

  # ALWAYS apply bounds to proceed
  lower_bounds_psi <- diag(s_xx) - (1 - min_reliability_marker) * diag(s_xx)
  toolow_idx <- which(diag(mm_psi) < lower_bounds_psi)
  if (length(toolow_idx) > 0L) {
    diag(mm_psi)[toolow_idx] <- lower_bounds_psi[toolow_idx]
  }
  too_large_idx <- which(diag(mm_psi) > diag(s_xx))
  if (length(too_large_idx) > 0L) {
    diag(mm_psi)[too_large_idx] <- diag(s_xx)[too_large_idx] * 1
  }

  # in addition, force PSI to be PD
  mm_psi <- lav_matrix_symmetric_force_pd(mm_psi, tol = 1e-04)

  # residual variances markers
  theta_x <- diag(s_xx - mm_psi)

  # create theta vector
  theta_nobounds <- numeric(nvar)
  theta_nobounds[marker_idx] <- theta_x
  theta_nobounds[-marker_idx] <- theta_f_nobounds

  # compute LAMBDA for non-marker items

  if (quadprog) {
    # only really needed if we need to impose (in)equality constraints
    # (TODO)
    dmat <- lav_matrix_bdiag(rep(list(mm_psi), p))
    dvec <- as.vector(t(s_yx))
    eq_idx <- which(t(b_y) != 1) # these must be zero (row-wise!)
    rmat <- diag(nrow(dmat))[eq_idx, , drop = FALSE]
    bvec <- rep(0, length(eq_idx)) # optional, 0=default
    out <- try(quadprog::solve.QP(
      Dmat = dmat, dvec = dvec, Amat = t(rmat),
      meq = length(eq_idx), bvec = bvec
    ), silent = TRUE)
    if (inherits(out, "try-error")) {
      lav_msg_warn(gettext("solve.QP failed to find a solution"))
      lambda_1 <- matrix(0, nvar, nfac)
      lambda_1[marker_idx, ] <- diag(nfac)
      lambda_1[lambda_nonzero_idx] <- as.numeric(NA)
      theta <- numeric(nvar)
      theta[marker_idx] <- theta_x
      theta[-marker_idx] <- theta_f
      return(list(
        lambda = lambda_1, theta = theta_nobounds,
        psi = psi_nobounds
      ))
    } else {
      lambda_y <- matrix(out$solution,
        nrow = p, ncol = nfac,
        byrow = TRUE
      )
      # zap almost zero elements
      lambda_y[abs(lambda_y) < sqrt(.Machine$double.eps)] <- 0
    }
  } else {
    # simple version
    #LAMBDA.y <- t(t(S.yx) / diag(PSI)) * B.y # works only if no crossloadings
    lambda_y <- t(solve(mm_psi, t(s_yx))) * b_y
  }


  # assemble matrices
  mm_lambda <- matrix(0, nvar, nfac)
  mm_lambda[marker_idx, ] <- diag(nfac)
  mm_lambda[-marker_idx, ] <- lambda_y

  list(lambda = mm_lambda, theta = theta_nobounds, psi = psi_nobounds)
}

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_cfa_bentler1982_internal <- function(lavobject = NULL, # convenience
                                         # internal slot
                                         lavmodel = NULL,
                                         lavsamplestats = NULL,
                                         lavpartable = NULL,
                                         lavdata = NULL,
                                         lavoptions = NULL,
                                         gls = TRUE,
                                         min_reliability_marker = 0.1,
                                         quadprog = FALSE,
                                         nobs = 20L) {
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
    lav_msg_stop(gettext("bentler1982 estimator only available for CFA models"))
  }
  # no BETA matrix! (i.e., no higher-order factors)
  if (!is.null(lavmodel@GLIST$beta)) {
    lav_msg_stop(gettext("bentler1982 estimator not available
                 for models that require a BETA matrix"))
  }
  # no std.lv = TRUE for now
  if (lavoptions$std.lv) {
    lav_msg_stop(gettext(
      "bentler1982 estimator not available if std.lv = TRUE"))
  }

  nblocks <- lav_partable_nblocks(lavpartable)
  stopifnot(nblocks == 1L) # for now
  b <- 1L
  sample_cov <- lavsamplestats@cov[[b]]
  nvar <- nrow(sample_cov)
  marker_idx <- lavpta$vidx$lv.marker[[b]]
  lambda_idx <- which(names(lavmodel@GLIST) == "lambda")
  lambda_nonzero_idx <- lavmodel@m.free.idx[[lambda_idx]]
  # only diagonal THETA for now...
  # because if we have correlated residuals, we should remove the
  # corresponding variables as instruments before we estimate lambda...
  # (see MIIV)
  theta_idx <- which(names(lavmodel@GLIST) == "theta") # usually '2'
  m_theta <- lavmodel@m.free.idx[[theta_idx]]
  nondiag_idx <- m_theta[!m_theta %in% lav_matrix_diag_idx(nvar)]
  if (length(nondiag_idx) > 0L) {
    lav_msg_warn(gettext(
      "this implementation of FABIN does not handle correlated residuals yet!"))
  }

  if (!missing(gls)) {
    gls_flag <- gls
  } else {
    gls_flag <- FALSE
    if (!is.null(lavoptions$estimator.args$GLS) &&
      lavoptions$estimator.args$GLS) {
      gls_flag <- TRUE
    }
  }

  if (missing(quadprog) &&
    !is.null(lavoptions$estimator.args$quadprog)) {
    quadprog <- lavoptions$estimator.args$quadprog
  }

  # run bentler1982 non-iterative CFA algorithm
  out <- lav_cfa_bentler1982(
    s = sample_cov, marker_idx = marker_idx,
    lambda_nonzero_idx = lambda_nonzero_idx,
    gls = gls_flag,
    min_reliability_marker = 0.1,
    quadprog = quadprog,
    nobs = lavsamplestats@ntotal
  )
  mm_lambda <- out$lambda
  mm_theta <- diag(out$theta, nvar)
  mm_psi <- out$psi

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
