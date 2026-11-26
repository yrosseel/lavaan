# approximate bootstrap draws, without refitting the model
#
# The (ordinary, nonparametric) bootstrap can be written as a reweighting
# scheme: each bootstrap sample corresponds to a vector of (multinomial)
# case weights w, and the estimator is a smooth function theta(w) of these
# weights. Instead of re-optimizing for every w, theta(w) is approximated
# by a Taylor expansion around the original (unit) weights:
#
# - method = "ij": the first-order expansion, or infinitesimal jackknife:
#
#       theta(w) ~ theta + E.inv %*% (1/N) sum_i (w_i - 1) s_i
#
#   where s_i are the casewise scores, and E.inv is the inverse of the
#   observed (unit) information. The covariance matrix of these draws is
#   (up to Monte Carlo error) the Huber-White sandwich estimator.
#
# - method = "hoij2": the second-order expansion (the second-order
#   infinitesimal jackknife). Let w(t) = 1 + t * (w - 1) be the path from
#   the original weights (t = 0) to the bootstrap weights (t = 1), and let
#   G(theta, w) be the gradient of the objective function, evaluated at the
#   reweighted sample statistics. Along the (first-order) path
#
#       r(t) = G(theta + t * theta.dot, w(t))
#
#   we have r(0) = 0 and r'(0) = 0, and (implicit function theorem)
#
#       theta(w) ~ theta + theta.dot - 0.5 * E.inv %*% r''(0)
#
#   where r''(0) is computed by finite differences (four gradient
#   evaluations per bootstrap sample). No casewise Hessians or arrays of
#   third-order derivatives are needed.
#
# The draws are returned in exactly the same format as the draws of
# lav_bootstrap_internal(fun = "coef"), so that everything downstream
# (standard errors, defined parameters, confidence intervals, ...) can be
# used unchanged.
#
# YR/LV 26 Nov 2026: first version (method = "ij")
# YR/LV 26 Nov 2026: add method = "hoij2" (complete data only)

# check if the approximate bootstrap is available for this model; if not,
# stop with an informative message
lav_bootstrap_ij_check <- function(lavmodel = NULL,
                                   lavdata = NULL,
                                   lavoptions = NULL,
                                   method = "ij") {
  not_available <- function(reason) {
    lav_msg_stop(gettextf(
      "The bootstrap method %1$s is not available (yet) %2$s; please use
       bootstrap = list(method = \"refit\") instead.", dQuote(method), reason))
  }
  if (!method %in% c("ij", "hoij2")) {
    lav_msg_stop(gettextf(
      "The bootstrap method %s is not available (yet).", dQuote(method)))
  }
  if (lavmodel@estimator != "ML") {
    not_available(gettext("if the estimator is not ML"))
  }
  if (lavmodel@categorical) {
    not_available(gettext("for categorical data"))
  }
  if (lavdata@nlevels > 1L) {
    not_available(gettext("for multilevel data"))
  }
  if (length(lavdata@cluster) > 0L) {
    not_available(gettext("for clustered data"))
  }
  if (lavmodel@conditional.x) {
    not_available(gettext("if conditional_x = TRUE"))
  }
  if (!lavoptions$missing %in% c("listwise", "ml", "ml.x")) {
    not_available(gettextf("if missing = %s", dQuote(lavoptions$missing)))
  }
  if ("bollen.stine" %in% lavoptions$test) {
    not_available(gettext("in combination with the Bollen-Stine bootstrap"))
  }
  if (length(lavoptions$rotation) > 0L && lavoptions$rotation != "none" &&
      lavmodel@nefa > 0L) {
    not_available(gettext("for models with rotated (EFA/ESEM) blocks"))
  }
  if (!is.null(lavoptions$backend) &&
      !identical(lavoptions$backend, "lavaan")) {
    not_available(gettext("for models estimated by an external backend"))
  }

  # second-order: more restrictions (for now)
  if (method == "hoij2") {
    if (lavoptions$missing != "listwise") {
      not_available(gettextf("if missing = %s", dQuote(lavoptions$missing)))
    }
    if (length(lavdata@sampling.weights) > 0L) {
      not_available(gettext("in combination with sampling weights"))
    }
    if (!is.null(lavoptions$likelihood) &&
        lavoptions$likelihood != "normal") {
      not_available(gettextf("if likelihood = %s",
                             dQuote(lavoptions$likelihood)))
    }
    if (lavmodel@correlation) {
      not_available(gettext("for correlation structures"))
    }
    if (lavmodel@group.w.free) {
      not_available(gettext("if group_w_free = TRUE"))
    }
    if (lavmodel@ceq.simple.only) {
      not_available(gettext("if ceq_simple = TRUE"))
    }
    if (length(lavmodel@ceq.nonlinear.idx) > 0L) {
      not_available(gettext("for models with nonlinear equality constraints"))
    }
    # inequality constraints, or bounds
    n_ceq <- length(lavmodel@ceq.linear.idx)
    if (length(lavmodel@cin.linear.idx) > 0L ||
        length(lavmodel@cin.nonlinear.idx) > 0L ||
        nrow(lavmodel@con.jac) > n_ceq) {
      not_available(
        gettext("for models with inequality constraints or bounds"))
    }
  }

  invisible(TRUE)
}

# the casewise scores (of the loglikelihood), per group; the rows of the
# g-th element correspond to the rows of lavdata@X[[g]]
lav_bootstrap_ij_scores <- function(lavmodel = NULL,
                                    lavsamplestats = NULL,
                                    lavdata = NULL,
                                    lavoptions = NULL,
                                    lavimplied = NULL) {
  ngroups <- lavdata@ngroups
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)

  # with simple equality constraints, the scores (and the information) live
  # in the larger 'unco' space
  npar <- lavmodel@nx.free
  if (lavmodel@ceq.simple.only) {
    npar <- lavmodel@nx.unco
  }

  # (model-implied) moments, in the format expected by lav_sc_ml()
  if (is.null(lavimplied) || length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel)
  }
  moments <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    moments[[g]] <- lapply(lavimplied, "[[", g)
  }
  if (ngroups == 1L) {
    moments <- moments[[1]]
  }

  score_matrix <- lav_sc_ml(
    ntab = ntab, ntot = ntot, npar = npar,
    moments = moments, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavmodel = lavmodel, lavoptions = lavoptions, scaling = FALSE
  )

  out <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    sc_g <- score_matrix[lavdata@case.idx[[g]], , drop = FALSE]
    # sampling weights: the casewise contributions are wt_i * s_i
    if (length(lavdata@sampling.weights) > 0L &&
        !is.null(lavdata@weights[[g]])) {
      sc_g <- sc_g * lavdata@weights[[g]]
    }
    out[[g]] <- sc_g
  }

  out
}

# the bootstrap samples, as frequencies (weights), per group; the r rows
# are the bootstrap samples, the columns are the rows of lavdata@X[[g]]
#
# the samples are drawn in the same way (and in the same order) as in
# lav_bootstrap_internal(): for the same value of iseed, bootstrap sample b
# is identical to the one that is used for (serial) refit b
lav_bootstrap_ij_freq <- function(lavdata = NULL, r = 1000L, iseed = NULL) {
  ngroups <- lavdata@ngroups

  # exclude the empty cases from the resampling pool
  pool <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    pool[[g]] <- seq_len(nrow(lavdata@X[[g]]))
    empty_idx <- lavdata@Mp[[g]]$empty.idx
    if (length(empty_idx) > 0L) {
      pool[[g]] <- pool[[g]][-empty_idx]
    }
  }

  boot_idx <- lav_bootstrap_indices(
    r = r, nobs = lapply(pool, length), iseed = iseed
  )

  freq <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    n_pool <- length(pool[[g]])
    n_g <- nrow(lavdata@X[[g]])
    # frequency table: one row per bootstrap sample
    cell <- (boot_idx[[g]] - 1L) * r + as.vector(row(boot_idx[[g]]))
    freq_pool <- matrix(tabulate(cell, nbins = r * n_pool), r, n_pool)
    freq_g <- matrix(0L, r, n_g)
    freq_g[, pool[[g]]] <- freq_pool
    freq[[g]] <- freq_g
  }

  # cases that can not be selected keep their (unit) weight
  attr(freq, "pool") <- pool
  # the bootstrap samples as row indices (to refit the model, if needed)
  attr(freq, "boot.idx") <- lapply(seq_len(ngroups), function(g) {
    matrix(pool[[g]][boot_idx[[g]]], nrow = r)
  })
  freq
}

# flag the draws with (at least one) negative variance
lav_bootstrap_ij_nonadmissible <- function(coef_boot = NULL,
                                           lavpartable = NULL) {
  if (is.null(lavpartable)) {
    return(integer(0L))
  }
  var_idx <- unique(lavpartable$free[lavpartable$free > 0L &
    lavpartable$op == "~~" & lavpartable$lhs == lavpartable$rhs])
  var_idx <- var_idx[var_idx <= ncol(coef_boot)]
  if (length(var_idx) == 0L) {
    return(integer(0L))
  }
  which(rowSums(coef_boot[, var_idx, drop = FALSE] < 0) > 0L)
}

# main function: return a (r x npar) matrix with approximate bootstrap draws
lav_bootstrap_ij <- function(lavmodel = NULL,
                             lavsamplestats = NULL,
                             lavdata = NULL,
                             lavoptions = NULL,
                             lavpartable = NULL,
                             lavimplied = NULL,
                             lavh1 = NULL,
                             lavcache = NULL,
                             r = 1000L,
                             method = "ij",
                             check_post = TRUE,
                             freq = NULL) {
  lav_bootstrap_ij_check(
    lavmodel = lavmodel, lavdata = lavdata, lavoptions = lavoptions,
    method = method
  )
  r <- as.integer(r)

  # seed: same logic as in lav_bootstrap_internal()
  iseed <- lavoptions[["iseed"]]
  if (is.null(iseed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      runif(1)
    }
    iseed <- runif(1, 0, 999999999)
  }

  # original estimates
  x0 <- lav_model_get_parameters(lavmodel = lavmodel, type = "free")

  # casewise scores (per group)
  sc <- lav_bootstrap_ij_scores(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavoptions = lavoptions, lavimplied = lavimplied
  )

  # inverse of the observed (unit) information; the augmented version takes
  # care of the (equality) constraints, if any
  lavoptions_info <- lavoptions
  lavoptions_info$information[1] <- "observed"
  e_inv <- lav_model_info(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavcache = lavcache, lavimplied = lavimplied, lavh1 = lavh1,
    lavoptions = lavoptions_info, extra = FALSE, augmented = TRUE,
    inverted = TRUE
  )
  if (inherits(e_inv, "try-error") || anyNA(e_inv) ||
      ncol(e_inv) != ncol(sc[[1]])) {
    lav_msg_stop(gettextf(
      "The bootstrap method %s needs the inverse of the observed information
       matrix, but it could not be computed.", dQuote(method)))
  }

  # bootstrap weights; these can also be provided by the caller, as a list
  # (per group) of weight matrices (not necessarily integers; internal use)
  if (is.null(freq)) {
    freq <- lav_bootstrap_ij_freq(lavdata = lavdata, r = r, iseed = iseed)
  } else {
    stopifnot(is.list(freq), length(freq) == lavdata@ngroups)
    r <- nrow(freq[[1]])
    attr(freq, "pool") <- lapply(lavdata@X, function(x) seq_len(nrow(x)))
  }

  # weighted sum of the casewise scores: sum_i (w_i - 1) s_i
  g_delta <- matrix(0, r, ncol(e_inv))
  for (g in seq_len(lavdata@ngroups)) {
    pool <- attr(freq, "pool")[[g]]
    sc_pool <- sc[[g]][pool, , drop = FALSE]
    g_delta <- g_delta + (freq[[g]][, pool, drop = FALSE] %*% sc_pool -
      matrix(colSums(sc_pool), r, ncol(sc_pool), byrow = TRUE))
  }

  # first-order change in the parameters (e_inv is symmetric)
  d_theta <- (g_delta / lavsamplestats@ntotal) %*% e_inv

  # second-order term
  error_idx <- integer(0L)
  guard_dist <- NULL
  if (method == "hoij2") {
    tmp <- lav_bootstrap_hoij2(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavpartable = lavpartable, lavcache = lavcache, freq = freq, x0 = x0,
      d_theta = d_theta, e_inv = e_inv
    )
    d_theta <- tmp$d_theta + tmp$d2_theta
    error_idx <- tmp$error_idx
    guard_dist <- tmp$dist
  }

  # simple equality constraints: from the 'unco' space to the free parameters
  if (lavmodel@ceq.simple.only) {
    k_matrix <- lavmodel@ceq.simple.K
    d_theta <- d_theta %*% k_matrix %*% solve(crossprod(k_matrix))
  }

  coef_boot <- d_theta + matrix(x0, r, length(x0), byrow = TRUE)
  nonadmissible <- lav_bootstrap_ij_nonadmissible(coef_boot,
    lavpartable = lavpartable)

  # guard: refit the model for the bootstrap samples where the approximation
  # is not accurate enough (or failed)
  guard_info <- NULL
  boot_opt <- lavoptions$bootstrap
  if (method == "hoij2" && is.list(boot_opt) && isTRUE(boot_opt$guard) &&
      !is.null(attr(freq, "boot.idx"))) {
    guard_tol <- boot_opt$guard_tol
    if (is.null(guard_tol)) {
      guard_tol <- 0.25
    }
    refit_idx <- which(!is.finite(guard_dist) | guard_dist > guard_tol)
    refit_idx <- sort(union(refit_idx, error_idx))
    if (length(refit_idx) > 0L) {
      coef_refit <- lav_bootstrap_internal(
        object = NULL, lavmodel = lavmodel, lavsamplestats = lavsamplestats,
        lavpartable = lavpartable, lavoptions = lavoptions,
        lavdata = lavdata, r = length(refit_idx), check_post = TRUE,
        type = "ordinary", fun = "coef",
        boot_idx_fixed = lapply(attr(freq, "boot.idx"), function(x) {
          x[refit_idx, , drop = FALSE]
        })
      )
      coef_boot[refit_idx, ] <- coef_refit[, , drop = FALSE]
      error_idx <- refit_idx[attr(coef_refit, "error.idx")]
      nonadmissible <- sort(c(setdiff(nonadmissible, refit_idx),
                              refit_idx[attr(coef_refit, "nonadmissible")]))
    }
    guard_info <- list(tol = guard_tol, refit.idx = refit_idx)
    if (length(refit_idx) > r / 2) {
      lav_msg_note(gettextf(
        "The bootstrap method %1$s was not accurate enough for %2$s of the
         %3$s bootstrap samples (for which the model was refitted instead);
         consider using bootstrap = list(method = \"refit\").",
        dQuote(method), length(refit_idx), r))
    }
  }

  # same attributes as lav_bootstrap_internal()
  attr(coef_boot, "error.idx") <- error_idx
  if (check_post) {
    attr(coef_boot, "nonadmissible") <- nonadmissible
  }
  attr(coef_boot, "guard") <- guard_info
  attr(coef_boot, "seed") <- iseed
  attr(coef_boot, "method") <- method
  attr(coef_boot, "dist") <- guard_dist

  coef_boot
}

# the second-order term of the second-order infinitesimal jackknife
#
# - freq: list (per group) of (r x n_g) weight matrices
# - x0: the original estimates
# - d_theta: the (r x npar) first-order changes (score based)
# - e_inv: the (augmented) inverse of the observed (unit) information
#
# returns a list with
# - d_theta: the first-order changes; these are the score-based ones, unless
#     the model contains fixed exogenous covariates (fixed.x = TRUE): their
#     (co)variances and means are fixed to the sample values, which change
#     under reweighting; the first-order changes are then recomputed (by
#     finite differences) to take this into account
# - d2_theta: the second-order changes
# - error_idx: bootstrap samples for which the gradient could not be
#     computed along the path (the corresponding rows are set to NA)
lav_bootstrap_hoij2 <- function(lavmodel = NULL,
                                lavsamplestats = NULL,
                                lavdata = NULL,
                                lavpartable = NULL,
                                lavcache = NULL,
                                freq = NULL,
                                x0 = NULL,
                                d_theta = NULL,
                                e_inv = NULL,
                                h = 0.05) {
  ngroups <- lavdata@ngroups
  r <- nrow(d_theta)
  npar <- length(x0)

  # per group: the first and second (raw) moments, and their changes (per
  # unit of t) for every bootstrap sample
  mom <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    y <- lavdata@X[[g]]
    n_g <- nrow(y)
    nvar <- ncol(y)
    vech_idx <- which(lower.tri(diag(nvar), diag = TRUE), arr.ind = TRUE)
    yy <- y[, vech_idx[, 1L], drop = FALSE] * y[, vech_idx[, 2L], drop = FALSE]
    dw <- freq[[g]] - 1
    mean_g <- colMeans(y)
    yy_g <- colMeans(yy)
    mom[[g]] <- list(
      nvar = nvar, vech_idx = vech_idx, mean = mean_g, yy = yy_g,
      d_mean = (dw %*% y) / n_g, d_yy = (dw %*% yy) / n_g
    )
    # sanity check: can we reproduce the sample statistics?
    cov_g <- matrix(0, nvar, nvar)
    cov_g[vech_idx] <- yy_g
    cov_g[vech_idx[, 2:1, drop = FALSE]] <- yy_g
    cov_g <- cov_g - tcrossprod(mean_g)
    if (max(abs(cov_g - lavsamplestats@cov[[g]])) >
        1e-8 * max(1, abs(cov_g)) ||
        max(abs(mean_g - lavsamplestats@mean[[g]])) >
        1e-8 * max(1, abs(mean_g))) {
      lav_msg_stop(gettext(
        "The bootstrap method \"hoij2\" could not reproduce the sample
         statistics from the raw data; please use
         bootstrap = list(method = \"refit\") instead."))
    }
  }

  # fixed exogenous covariates: where do their (fixed) sample statistics
  # live in GLIST?
  exo_map <- lav_bootstrap_hoij2_exo_map(
    lavmodel = lavmodel, lavdata = lavdata, lavpartable = lavpartable
  )
  exo_flag <- length(exo_map) > 0L

  # the gradient, at theta, for bootstrap sample b, at position t on the path
  grad_fn <- function(theta, b, t) {
    for (g in seq_len(ngroups)) {
      mm <- mom[[g]]
      mean_t <- mm$mean + t * mm$d_mean[b, ]
      yy_t <- mm$yy + t * mm$d_yy[b, ]
      cov_t <- matrix(0, mm$nvar, mm$nvar)
      cov_t[mm$vech_idx] <- yy_t
      cov_t[mm$vech_idx[, 2:1, drop = FALSE]] <- yy_t
      cov_t <- cov_t - tcrossprod(mean_t)
      lavsamplestats@mean[[g]] <- mean_t
      lavsamplestats@cov[[g]] <- cov_t
    }
    glist <- lav_model_x2glist(lavmodel = lavmodel, x = theta)
    for (el in exo_map) {
      if (el$type == "cov") {
        value <- lavsamplestats@cov[[el$group]][el$idx]
      } else {
        value <- lavsamplestats@mean[[el$group]][el$idx]
      }
      glist[[el$mm]][el$m_idx] <- value
    }
    out <- try(suppressWarnings(lav_model_grad(
      lavmodel = lavmodel, glist = glist, lavsamplestats = lavsamplestats,
      lavdata = lavdata, lavcache = lavcache, type = "free"
    )), silent = TRUE)
    if (inherits(out, "try-error") || length(out) != npar ||
        !all(is.finite(out))) {
      out <- rep(as.numeric(NA), npar)
    }
    as.numeric(out)
  }

  # (should be zero)
  r0 <- grad_fn(x0, 1L, 0)

  # effective number of parameters (rank of e_inv)
  ntotal <- lavsamplestats@ntotal
  npar_eff <- npar
  if (lavmodel@eq.constraints) {
    npar_eff <- ncol(lavmodel@eq.constraints.K)
  }

  d1_theta <- d_theta
  d2_theta <- matrix(as.numeric(NA), r, npar)
  dist <- rep(as.numeric(NA), r)
  for (b in seq_len(r)) {
    # fixed.x: recompute the first-order change, including the effect of
    # the changing sample statistics of the exogenous covariates
    if (exo_flag) {
      g_dot <- (-grad_fn(x0, b, 2 * h) + 8 * grad_fn(x0, b, h) -
        8 * grad_fn(x0, b, -h) + grad_fn(x0, b, -2 * h)) / (12 * h)
      d1_theta[b, ] <- -1 * drop(e_inv %*% g_dot)
    }
    dir <- d1_theta[b, ]
    # second derivative of r(t) at t = 0 (five-point stencil)
    r2 <- (-grad_fn(x0 + 2 * h * dir, b, 2 * h) +
      16 * grad_fn(x0 + h * dir, b, h) - 30 * r0 +
      16 * grad_fn(x0 - h * dir, b, -h) -
      grad_fn(x0 - 2 * h * dir, b, -2 * h)) / (12 * h * h)
    d2_theta[b, ] <- -0.5 * drop(e_inv %*% r2)
    # how far are we from the exact solution for this bootstrap sample?
    # the gradient at the approximate solution (zero if it is exact),
    # expressed as the length of the corresponding Newton step, in the
    # metric of the (total) information: an estimate of the distance to
    # the exact solution in standard error units (root mean square over
    # the parameters)
    g1 <- grad_fn(x0 + d1_theta[b, ] + d2_theta[b, ], b, 1)
    dist[b] <- sqrt(max(0, ntotal * sum(g1 * drop(e_inv %*% g1)) / npar_eff))
  }

  error_idx <- which(!stats::complete.cases(d1_theta + d2_theta))

  list(d_theta = d1_theta, d2_theta = d2_theta, error_idx = error_idx,
       dist = dist)
}

# fixed.x = TRUE: the (co)variances and the means of the exogenous
# covariates are fixed to their sample values, and stored in GLIST; return a
# list with (per model matrix) the elements that need to be updated if the
# sample statistics change
lav_bootstrap_hoij2_exo_map <- function(lavmodel = NULL,
                                        lavdata = NULL,
                                        lavpartable = NULL) {
  out <- list()
  if (!lavmodel@fixed.x || is.null(lavpartable) ||
      is.null(lavpartable$exo)) {
    return(out)
  }
  exo_id <- lavpartable$id[lavpartable$exo == 1L & lavpartable$free == 0L &
    lavpartable$op %in% c("~~", "~1")]
  if (length(exo_id) == 0L) {
    return(out)
  }
  pt_row <- match(exo_id, lavpartable$id)
  nmat <- lavmodel@nmat
  for (mm in seq_along(lavmodel@GLIST)) {
    keep <- which(lavmodel@x.user.idx[[mm]] %in% exo_id)
    if (length(keep) == 0L) {
      next
    }
    g <- which(mm <= cumsum(nmat))[1L]
    ov_names <- lavdata@ov.names[[g]]
    row_idx <- pt_row[match(lavmodel@x.user.idx[[mm]][keep], exo_id)]
    op <- lavpartable$op[row_idx]
    lhs <- match(lavpartable$lhs[row_idx], ov_names)
    rhs <- match(lavpartable$rhs[row_idx], ov_names)
    m_idx <- lavmodel@m.user.idx[[mm]][keep]
    cov_idx <- which(op == "~~")
    if (length(cov_idx) > 0L) {
      out[[length(out) + 1L]] <- list(
        mm = mm, group = g, type = "cov", m_idx = m_idx[cov_idx],
        idx = cbind(lhs[cov_idx], rhs[cov_idx])
      )
    }
    mean_idx <- which(op == "~1")
    if (length(mean_idx) > 0L) {
      out[[length(out) + 1L]] <- list(
        mm = mm, group = g, type = "mean", m_idx = m_idx[mean_idx],
        idx = lhs[mean_idx]
      )
    }
  }

  out
}

# the 'method' element of the bootstrap option
lav_bootstrap_method <- function(lavoptions = NULL) {
  method <- "refit"
  if (is.list(lavoptions$bootstrap) &&
      !is.null(lavoptions$bootstrap$method)) {
    method <- lavoptions$bootstrap$method
  }
  method
}
