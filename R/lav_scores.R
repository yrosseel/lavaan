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
#                 move ML-specific code to lav_sc_ml() function
# YR 26 Apr 2025: add lav_scores_gls()

lav_sc <- function(object, scaling = FALSE,
                           ignore_constraints = FALSE,
                           remove_duplicated = TRUE,
                           remove_empty_cases = TRUE) {
  stopifnot(inherits(object, "lavaan"))

  # check object
  object <- lav_object_check_version(object)

  # what if estimator is not ML or WLS?
  # avoid hard error (using stop); throw a warning, and return an empty matrix
  if (!object@Options$estimator %in% c("ML", "WLS", "GLS", "ULS")) {
    lav_msg_warn(gettext("scores only available if estimator is ML"))
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
  npar <- lav_inspect_npar(object, ceq = FALSE)
  # with simple equality constraints (ceq.simple), the casewise scores are first
  # computed in the larger 'unco' space (one column per non-collapsed parameter,
  # matching Delta); they are collapsed to the free parameters further below
  if (lavmodel@ceq.simple.only) {
    npar <- lavmodel@nx.unco
  }

  if (object@Options$estimator == "ML") {
    moments <- fitted(object)
    score_matrix <- lav_sc_ml(
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
    score_matrix <- lav_sc_wls(
      ntab = ntab, ntot = ntot, npar = npar,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavmodel = lavmodel, lavoptions = lavoptions
    )
  } else if (!lavmodel@categorical &&
             object@Options$estimator %in% c("GLS", "ULS", "WLS")) {
    # compute WLS/GLS/ULS `scores'
    score_matrix <- lav_sc_ls(
      ntab = ntab, ntot = ntot, npar = npar,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavmodel = lavmodel, lavoptions = lavoptions
    )
  } else {
    # should not happen
    lav_msg_fixme("this should not happen")
  }

  # sampling weights? the casewise score contributions become wt_i * s_i,
  # so that the column sums equal the (weighted) gradient (= 0 at the MLE),
  # and crossprod() yields the design-based information that is also used
  # for the (robust) standard errors. This is a no-op without weights.
  if (length(lavdata@sampling.weights) > 0L) {
    wt_full <- numeric(ntot)
    for (g in seq_len(lavdata@ngroups)) {
      wt_full[lavdata@case.idx[[g]]] <- lavdata@weights[[g]]
    }
    score_matrix <- score_matrix * wt_full
  }

  # handle empty rows
  if (remove_empty_cases) {
    # empty.idx <- which( apply(score_matrix, 1L,
    #                        function(x) sum(is.na(x))) == ncol(score_matrix) )
    empty_idx <- unlist(lapply(lavdata@Mp, "[[", "empty.idx"))
    if (length(empty_idx) > 0L) {
      score_matrix <- score_matrix[-empty_idx, , drop = FALSE]
    }
  }

  # provide column names
  colnames(score_matrix) <- names(lav_inspect_coef(object,
    type = "free", add_labels = TRUE
  ))

  # handle general constraints, so that the sum of the columns equals zero
  # note: for ceq.simple.only models there are no explicit ceq/cin indices, but
  # the (synthesized) con.jac encodes the simple equalities in the 'unco' space;
  # projecting here makes the full (remove.duplicated = FALSE) scores respect the
  # constraints, exactly as for explicit equality constraints. This projection
  # does NOT affect the collapsed scores below, since con.jac %*% ceq.simple.K = 0
  if (!ignore_constraints &&
    (sum(
      lavmodel@ceq.linear.idx, lavmodel@ceq.nonlinear.idx,
      lavmodel@cin.linear.idx, lavmodel@cin.nonlinear.idx
    ) > 0 || (lavmodel@ceq.simple.only && nrow(lavmodel@con.jac) > 0L))) {
    r_matrix <- object@Model@con.jac[, , drop = FALSE]
    pre <- lav_con_lambda_pre(object)
    # LAMBDA <- -1 * t(pre %*% t(score_matrix))
    # RLAMBDA <- t(t(r_matrix) %*% t(LAMBDA))
    score_matrix <- score_matrix - t(t(r_matrix) %*% pre %*% t(score_matrix))
  }

  # handle simple equality constraints
  if (remove_duplicated && lavmodel@eq.constraints) {
    simple_flag <- lav_con_check_simple(lavmodel)
    if (simple_flag) {
      k_matrix <- lav_con_r2k(lavmodel)
      score_matrix <- score_matrix %*% k_matrix
    } else {
      lav_msg_warn(gettext(
        "remove.duplicated is TRUE, but equality constraints do not appear
        to be simple; returning full scores"))
    }
  } else if (remove_duplicated && lavmodel@ceq.simple.only) {
    # collapse the 'unco'-space scores to the free parameters: the score for a
    # free parameter is the sum of the scores of the parameters that share it
    # (ceq.simple.K is the nx.unco x nx.free duplication matrix)
    score_matrix <- score_matrix %*% lavmodel@ceq.simple.K
  }

  score_matrix
}

lavScores <- estfun.lavaan <- function(               # nolint
                       object,
                       scaling = FALSE,
                       ignore_constraints = FALSE,
                       remove_duplicated = TRUE,
                       remove_empty_cases = TRUE,
                      ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  lav_sc(object, scaling = scaling,
        ignore_constraints = ignore_constraints,
        remove_duplicated = remove_duplicated,
        remove_empty_cases = remove_empty_cases)
}

# YR 18 Jun 2026: the casewise H1 (saturated-model) scores are now computed by
#                 the (better tested) lav_mvn_sc_mu_sigma() (complete data) and
#                 lav_mvn_mi_sc_mu_sigma() (incomplete data) functions from the
#                 lav_mvnorm*.R files; this function only assembles them into the
#                 model-parameter scores via the Delta matrix. The mvn casewise
#                 score has the opposite sign of the old 'scores_h1' (so the old
#                 '-scores_h1 %*% Delta' becomes 'sc %*% Delta'), and produces
#                 vech(Sigma) with the diagonal already halved -- matching Delta.
lav_sc_ml <- function(ntab = 0L,
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

  # group weights (only used if scaling = TRUE)
  group_w <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal

  for (g in 1:lavsamplestats@ngroups) {
    if (lavsamplestats@ngroups > 1) {
      moments <- moments_groups[[g]]
    }
    sigma_hat <- moments$cov
    nvar <- ncol(lavsamplestats@cov[[g]])

    # H1 (saturated) mean: model-implied if meanstructure, else sample mean
    if (lavmodel@meanstructure) {
      mu_hat <- moments$mean
    } else {
      mu_hat <- lavsamplestats@mean[[g]]
    }

    if (!lavsamplestats@missing.flag) { # complete data
      sc <- lav_mvn_sc_mu_sigma(
        y = lavdata@X[[g]], mu = mu_hat, sigma_1 = sigma_hat
      )

      # wishart likelihood: the second-order (vech(Sigma)) block uses the
      # unbiased cross-product, i.e. tcrossprod(y - mu) is scaled by n/(n-1).
      # This rescales the vech(Sigma) columns of the casewise score, and adds a
      # (constant) shift proportional to vech(Sigma.inv).
      if (lavoptions$likelihood == "wishart") {
        nobs1 <- lavsamplestats@nobs[[g]] / (lavsamplestats@nobs[[g]] - 1)
        j2 <- matrix(1, nvar, nvar)
        diag(j2) <- 0.5
        d_row <- lav_mat_vech(-j2 * chol2inv(chol(sigma_hat)))
        sig_idx <- seq.int(nvar + 1L, ncol(sc))
        sc[, sig_idx] <- nobs1 * sc[, sig_idx, drop = FALSE] -
          matrix((nobs1 - 1) * d_row, nrow(sc), length(d_row), byrow = TRUE)
      }
    } else { # incomplete data
      sc <- lav_mvn_mi_sc_mu_sigma(
        y = lavdata@X[[g]], mp = lavdata@Mp[[g]],
        mu = mu_hat, sigma_1 = sigma_hat
      )
      # lav_mvn_mi_sc_mu_sigma() leaves NA in the score entries that refer to
      # the *unobserved* cells of a case; the (saturated) casewise score for a
      # missing variable is zero, so we zero them here before mapping through
      # Delta (the old implementation initialised these entries to 0).
      sc[is.na(sc)] <- 0
    } # missing

    # no mean structure: drop the (leading) mu columns
    if (!lavmodel@meanstructure) {
      sc <- sc[, -seq_len(nvar), drop = FALSE]
    }

    # scaling? (matches the historical behaviour: group weights are applied to
    # the incomplete-data scores only, then everything is scaled by -1/ntot)
    if (scaling && lavsamplestats@missing.flag) {
      sc <- group_w[g] * sc
    }

    wi <- lavdata@case.idx[[g]]
    score_matrix[wi, ] <- sc %*% delta[[g]]
    if (scaling) {
      score_matrix[wi, ] <- (-1 / ntot) * score_matrix[wi, ]
    }
  } # g

  score_matrix
}

# this function is based on code originally written by Franz Classe (Munich)
# for categorical data only!
lav_sc_wls <- function(ntab = 0L,
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

  # container for scores
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
      lavaan::lav_mat_vech(tcrossprod(i), diagonal = FALSE)
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
lav_sc_ls <- function(ntab = 0L,
                          ntot = 0L,
                          npar = 0L,
                          lavdata = NULL,
                          lavsamplestats = NULL,
                          lavmodel = NULL,
                          lavoptions = NULL) {

  # container for scores
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
    idx1 <- lav_mat_vech_col_idx(nvar)
    idx2 <- lav_mat_vech_row_idx(nvar)
    if (lavmodel@meanstructure) {
      z <- cbind(y, yc[, idx1, drop = FALSE] * yc[, idx2, drop = FALSE])
    } else {
      z <- (yc[, idx1, drop = FALSE] * yc[, idx2, drop = FALSE])
    }

    # model-based sample statistics
    if (lavmodel@meanstructure) {
      sigma <- c(as.numeric(implied$mean[[g]]),
                 lav_mat_vech(implied$cov[[g]]))
    } else {
      sigma <- lav_mat_vech(implied$cov[[g]])
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
