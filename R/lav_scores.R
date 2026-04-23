# contributed by Ed Merkle (17 Jan 2013)
# WLS version contributed by Franz Classe (March 2024)
# (adapted for inclusion in lavaan by YR)


# YR 12 Feb 2013: small changes to match the results of lav_model_gradient
#                 in the multiple group case
# YR 30 May 2014: handle 1-variable case (fixing apply in lines 56, 62, 108)
# YR 05 Nov 2015: add remove.duplicated = TRUE, to cope with strucchange in
#                 case of simple equality constraints
# YR 19 Nov 2015: if constraints have been used, compute case-wise Lagrange
#                 multipliers, and define the scores as: SC + (t(R) lambda)
# YR 05 Feb 2016: catch conditional.x = TRUE: no support (for now), until
#                 we can use the generic 0.6 infrastructure for scores,
#                 including the missing-values case
# YR 16 Feb 2016: adapt to changed @Mp slot elements; add remove.empty.cases=
#                 argument
# YR 12 Mar 2024: make lintr (more) happy; include WLS code from Franz Classe
#                 move ML-specific code to lav_scores_ml() function
# YR 26 Apr 2025: add lav_scores_gls()

lav_scores <- function(object, scaling = FALSE,                    # nolint start
                                       ignore.constraints = FALSE,
                                       remove.duplicated = TRUE,
                                       remove.empty.cases = TRUE) { # nolint end
  stopifnot(inherits(object, "lavaan"))

  # check object
  object <- lav_object_check_version(object)

  # what if estimator is not ML or WLS?
  # avoid hard error (using stop); throw a warning, and return an empty matrix
  if (!object@Options$estimator %in% c("ML", "WLS", "GLS", "ULS")) {
    lav_msg_warn(gettext("scores only availalbe if estimator is ML"))
    return(matrix(0, 0, 0))
  }

  # check if conditional.x = TRUE
  if (object@Model@conditional.x) {
    lav_msg_stop(gettext("scores not available (yet) if conditional.x = TRUE"))
  }

  # shortcuts
  lavdata <- object@Data
  lavmodel <- object@Model
  lavsamplestats <- object@SampleStats
  lavoptions <- object@Options

  ## number variables/sample size
  # ntab <- unlist(lavsamplestats@nobs)
  ## change in 0.5-17: we keep the 'empty cases'
  ##                   and 'fill' in the scores at their 'case.idx'
  ##                   later, we remove the 'empty rows'
  # ntot <- max( object@Data@case.idx[[ object@Data@ngroups ]] )
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  npar <- lav_object_inspect_npar(object, ceq = FALSE)

  if (object@Options$estimator == "ML") {
    moments <- fitted(object)
    score_matrix <- lav_scores_ml(
      ntab = ntab, ntot = ntot, npar = npar,
      moments = moments, lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavmodel = lavmodel, lavoptions = lavoptions, scaling = scaling
    )
  } else if (object@Options$estimator == "WLS" && lavmodel@categorical) {
    # check if ALL observed variables are ordered
    ov_names <- unlist(lavdata@ov.names)
    ov_idx <- which(lavdata@ov$name %in% ov_names)
    if (!all(lavdata@ov$type[ov_idx] == "ordered")) {
      lav_msg_stop(gettext(
        "WLS scores only available if all observed variables are ordered."))
    }

    # compute WLS scores
    score_matrix <- lav_scores_wls(
      ntab = ntab, ntot = ntot, npar = npar,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavmodel = lavmodel, lavoptions = lavoptions
    )
  } else if (!lavmodel@categorical &&
             object@Options$estimator %in% c("GLS", "ULS", "WLS")) {
    # compute WLS/GLS/ULS `scores'
    score_matrix <- lav_scores_ls(
      ntab = ntab, ntot = ntot, npar = npar,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavmodel = lavmodel, lavoptions = lavoptions
    )
  } else {
    # should not happen
    lav_msg_fixme("this should not happen")
  }

  # handle empty rows
  if (remove.empty.cases) {
    # empty.idx <- which( apply(score_matrix, 1L,
    #                        function(x) sum(is.na(x))) == ncol(score_matrix) )
    empty_idx <- unlist(lapply(lavdata@Mp, "[[", "empty.idx"))
    if (length(empty_idx) > 0L) {
      score_matrix <- score_matrix[-empty_idx, , drop = FALSE]
    }
  }

  # provide column names
  colnames(score_matrix) <- names(lav_object_inspect_coef(object,
    type = "free", add.labels = TRUE
  ))

  # handle general constraints, so that the sum of the columns equals zero
  if (!ignore.constraints &&
    sum(
      lavmodel@ceq.linear.idx, lavmodel@ceq.nonlinear.idx,
      lavmodel@cin.linear.idx, lavmodel@cin.nonlinear.idx
    ) > 0) {
    r_matrix <- object@Model@con.jac[, ]
    pre <- lav_constraints_lambda_pre(object)
    # LAMBDA <- -1 * t(pre %*% t(score_matrix))
    # RLAMBDA <- t(t(r_matrix) %*% t(LAMBDA))
    score_matrix <- score_matrix - t(t(r_matrix) %*% pre %*% t(score_matrix))
  }

  # handle simple equality constraints
  if (remove.duplicated && lavmodel@eq.constraints) {
    simple_flag <- lav_constraints_check_simple(lavmodel)
    if (simple_flag) {
      k_matrix <- lav_constraints_r2k(lavmodel)
      score_matrix <- score_matrix %*% k_matrix
    } else {
      lav_msg_warn(gettext(
        "remove.duplicated is TRUE, but equality constraints do not appear
        to be simple; returning full scores"))
    }
  }

  score_matrix
}
lavScores <- lav_scores       # synonym #nolint
estfun.lavaan <- lav_scores   # synonym

lav_scores_ml <- function(ntab = 0L,
                          ntot = 0L,
                          npar = 0L,
                          moments = NULL,
                          lavdata = NULL,
                          lavsamplestats = NULL,
                          lavmodel = NULL,
                          lavoptions = NULL,
                          scaling = FALSE) {
  score_matrix <- matrix(NA, ntot, npar)

  # Delta matrix
  delta <- lav_model_delta(lavmodel = lavmodel)

  # rename moments
  moments_groups <- moments

  for (g in 1:lavsamplestats@ngroups) {
    if (lavsamplestats@ngroups > 1) {
      moments <- moments_groups[[g]]
    }
    sigma_hat <- moments$cov

    if (lavoptions$likelihood == "wishart") {
      nobs1 <- lavsamplestats@nobs[[g]] / (lavsamplestats@nobs[[g]] - 1)
    } else {
      nobs1 <- 1
    }

    if (!lavsamplestats@missing.flag) { # complete data
      # if(lavmodel@meanstructure) { # mean structure
      nvar <- ncol(lavsamplestats@cov[[g]])
      mu_hat <- moments$mean
      x_1 <- lavdata@X[[g]]
      sigma_inv <- chol2inv(chol(sigma_hat)) # FIXME: check for pd?
      group_w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)

      j <- matrix(1, 1L, ntab[g]) ## FIXME: needed?better maybe rowSums/colSums?
      j2 <- matrix(1, nvar, nvar)
      diag(j2) <- 0.5

      if (lavmodel@meanstructure) {
        ## scores_h1 (H1 = saturated model)
        mean_diff <- t(t(x_1) - mu_hat %*% j)
        dx_mu <- -1 * mean_diff %*% sigma_inv
        dx_sigma <- t(matrix(apply(
          mean_diff, 1L,
          function(x) {
            lav_matrix_vech(-j2 *
              (sigma_inv %*% (tcrossprod(x) * nobs1 - sigma_hat) %*% sigma_inv))
          }
        ), ncol = nrow(mean_diff)))

        scores_h1 <- cbind(dx_mu, dx_sigma)
      } else {
        mean_diff <- t(t(x_1) - lavsamplestats@mean[[g]] %*% j)
        dx_sigma <- t(matrix(apply(
          mean_diff, 1L,
          function(x) {
            lav_matrix_vech(-j2 *
              (sigma_inv %*% (tcrossprod(x) * nobs1 - sigma_hat) %*% sigma_inv))
          }
        ), ncol = nrow(mean_diff)))
        scores_h1 <- dx_sigma
      }
      ## FIXME? Seems like we would need group.w even in the
      ##        complete-data case:
      ## if(scaling){
      ##  scores_h1 <- group.w[g] * scores_h1
      ## }

      # } else {
      #  ## no mean structure
      #  stop("Score calculation with no mean structure is not implemented.")
      # }
    } else { # incomplete data
      nsub <- ntab[g]
      m <- lavsamplestats@missing[[g]]
      mp <- lavdata@Mp[[g]]
      # pat.idx <- match(MP1$id, MP1$order)
      group_w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)

      mu_hat <- moments$mean
      nvar <- ncol(lavsamplestats@cov[[g]])
      score_sigma <- matrix(0, nsub, nvar * (nvar + 1) / 2)
      score_mu <- matrix(0, nsub, nvar)

      for (p in seq_along(m)) {
        ## Data
        # X <- M[[p]][["X"]]
        case_idx <- mp$case.idx[[p]]
        var_idx <- m[[p]][["var.idx"]]
        x_1 <- lavdata@X[[g]][case_idx, var_idx, drop = FALSE]
        nobs <- m[[p]][["freq"]]
        ## Which unique entries of covariance matrix are estimated?
        ## (Used to keep track of scores in score.sigma)
        var_idx_mat <- tcrossprod(var_idx)
        sigma_idx <-
          which(var_idx_mat[lower.tri(var_idx_mat, diag = TRUE)] == 1)

        j <- matrix(1, 1L, nobs) # [var.idx]
        j2 <- matrix(1, nvar, nvar)[var_idx, var_idx, drop = FALSE]
        diag(j2) <- 0.5
        # FIXME: check for pd?
        sigma_inv <- chol2inv(chol(sigma_hat[var_idx, var_idx, drop = FALSE]))
        mu <- mu_hat[var_idx]
        mean_diff <- t(t(x_1) - mu %*% j)

        ## Scores for missing pattern p within group g
        score_mu[case_idx, var_idx] <- -1 * mean_diff %*% sigma_inv
        score_sigma[case_idx, sigma_idx] <- t(matrix(apply(
          mean_diff, 1L,
          function(x) {
            lav_matrix_vech(-j2 *
              (sigma_inv %*% (tcrossprod(x) -
                sigma_hat[var_idx, var_idx, drop = FALSE]) %*% sigma_inv))
          }
        ), ncol = nrow(mean_diff)))
      }

      scores_h1 <- cbind(score_mu, score_sigma)
      if (scaling) {
        scores_h1 <- group_w[g] * scores_h1
      }
    } # missing

    # if(lavmodel@eq.constraints) {
    #    Delta <- Delta %*% lavmodel@eq.constraints.K
    #    #x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
    #    #                lavmodel@eq.constraints.k0
    # }
    wi <- lavdata@case.idx[[g]]
    score_matrix[wi, ] <- -scores_h1 %*% delta[[g]]
    if (scaling) {
      score_matrix[wi, ] <- (-1 / ntot) * score_matrix[wi, ]
    }
  } # g

  score_matrix
}

# this function is based on code originally written by Franz Classe (Munich)
# for categorical data only!
lav_scores_wls <- function(ntab = 0L,
                           ntot = 0L,
                           npar = 0L,
                           lavdata = NULL,
                           lavsamplestats = NULL,
                           lavmodel = NULL,
                           lavoptions = NULL) {
  # internal function
  do_dummy_single_var <- function(x_1, lv, ntot, num) {
    xd <- matrix(NA, nrow = ntot, ncol = lv[num] - 1)
    x <- x_1[, num]
    minx <- min(x)
    categ <- minx - 1
    v <- 1
    while (categ < lv[num] - 1) {
      categ <- categ + 1
      xd[, v] <- ifelse(x > categ, 1, 0)
      v <- v + 1
    }

    xd
  }

  # containere for scores
  score_matrix <- matrix(NA, ntot, npar)

  # Delta matrix
  delta <- lav_model_delta(lavmodel = lavmodel)

  # shortcuts
  lv <- lavdata@ov[["nlev"]]

  for (g in 1:lavsamplestats@ngroups) {
    nvar <- ncol(lavsamplestats@cov[[g]])
    x_1 <- lavdata@X[[g]]

    # convert categorical data to dummy variables
    # FIXME: skip continuous variables
    xd <- do.call(
      cbind,
      lapply(
        1:nvar,
        function(i) do_dummy_single_var(x_1, lv, ntot, i)
      )
    )

    # e1
    musd <- colMeans(xd)
    e1 <- t(t(xd) - musd)

    # e2
    mus <- colMeans(x_1)
    y_minus_mu <- t(apply(x_1, 1L, function(x) x - mus))
    s_vech <- t(apply(y_minus_mu, 1L, function(i) {
      lavaan::lav_matrix_vech(tcrossprod(i), diagonal = FALSE)
    })) # s=c( (y1-mu1)(y2-mu2)....
    sigma <- colMeans(s_vech)
    e2 <- t(apply(s_vech, 1L, function(x) x - sigma))

    # e
    e <- cbind(e1, e2)

    # weight matrix
    m_w <- lavsamplestats@WLS.V[[g]]

    # combine matrices
    wi <- lavdata@case.idx[[g]]
    score_matrix[wi, ] <- t(t(delta[[g]]) %*% m_w %*% t(e))
  } # g

  score_matrix
}


# wls/gls/uls (continuous data only)
lav_scores_ls <- function(ntab = 0L,
                          ntot = 0L,
                          npar = 0L,
                          lavdata = NULL,
                          lavsamplestats = NULL,
                          lavmodel = NULL,
                          lavoptions = NULL) {

  # containere for scores
  score_matrix <- matrix(NA, ntot, npar)

  # estimator
  estimator <- lavoptions$estimator

  # Delta matrix
  delta <- lav_model_delta(lavmodel = lavmodel)

  # implied stats
  implied <- lav_model_implied(lavmodel)

  for (g in 1:lavsamplestats@ngroups) {
    nvar <- ncol(lavsamplestats@cov[[g]])
    nobs <- lavsamplestats@nobs[[g]]
    y <- lavdata@X[[g]]

    # center (not using model-implied!)
    yc <- t(t(y) - colMeans(y, na.rm = TRUE))

    # create Z where the rows_i contain the following elements:
    #  - Y_i (if meanstructure is TRUE)
    #  - vech(Yc_i' %*% Yc_i) where Yc_i are the residuals
    idx1 <- lav_matrix_vech_col_idx(nvar)
    idx2 <- lav_matrix_vech_row_idx(nvar)
    if (lavmodel@meanstructure) {
      z <- cbind(y, yc[, idx1, drop = FALSE] * yc[, idx2, drop = FALSE])
    } else {
      z <- (yc[, idx1, drop = FALSE] * yc[, idx2, drop = FALSE])
    }

    # model-based sample statistics
    if (lavmodel@meanstructure) {
      sigma <- c(as.numeric(implied$mean[[g]]),
                 lav_matrix_vech(implied$cov[[g]]))
    } else {
      sigma <- lav_matrix_vech(implied$cov[[g]])
    }

    # adjust sigma for N-1, so that colMeans(scores) == gradient
    # (not for the means)
    if (lavmodel@meanstructure) {
      z[, -seq_len(nvar)] <-
                z[, -seq_len(nvar), drop = FALSE] * nobs / (nobs - 1)
    } else {
      z <- z * nobs / (nobs - 1)
    }

    # compute Zc
    zc <- t(t(z) - sigma)

    # weight matrix
    if (estimator == "ULS") {
      m_w <- diag(ncol(z))
    } else {
      m_w <- lavsamplestats@WLS.V[[g]]
    }

    # combine matrices
    wi <- lavdata@case.idx[[g]]
    score_matrix[wi, ] <- zc %*% m_w %*% delta[[g]]
  } # g

  score_matrix
}
