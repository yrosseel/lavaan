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
# YR/LV 26 Nov 2026: cluster bootstrap; hoij2: missing = "ml", sampling
#                    weights, bounds and (linear) inequality constraints
# YR/LV 26 Nov 2026: least-squares estimators (DWLS, ULS), continuous or
#                    categorical data (see lav_bootstrap_ij_ls() below)

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
  ls_flag <- lavmodel@estimator %in% c("DWLS", "ULS")
  if (lavmodel@estimator != "ML" && !ls_flag) {
    not_available(gettext("if the estimator is not ML, DWLS or ULS"))
  }
  if (lavmodel@categorical && !ls_flag) {
    not_available(gettext("for categorical data, unless the estimator is
                           DWLS or ULS"))
  }
  # least-squares estimators: more restrictions (for now)
  if (ls_flag) {
    if (lavoptions$missing != "listwise") {
      not_available(gettextf("if missing = %1$s and estimator = %2$s",
                             dQuote(lavoptions$missing),
                             dQuote(lavmodel@estimator)))
    }
    if (length(lavdata@sampling.weights) > 0L) {
      not_available(gettextf(
        "in combination with sampling weights, if estimator = %s",
        dQuote(lavmodel@estimator)))
    }
    if (length(unlist(lavdata@ov.names.x)) > 0L) {
      not_available(gettextf(
        "for models with exogenous covariates, if estimator = %s",
        dQuote(lavmodel@estimator)))
    }
    if (lavmodel@correlation || lavmodel@group.w.free) {
      not_available(gettextf(
        "for correlation structures or free group weights, if
         estimator = %s", dQuote(lavmodel@estimator)))
    }
    if (lavmodel@ceq.simple.only) {
      not_available(gettextf("if ceq_simple = TRUE and estimator = %s",
                             dQuote(lavmodel@estimator)))
    }
  }
  if (lavdata@nlevels > 1L) {
    not_available(gettext("for multilevel data"))
  }
  if (length(lavdata@cluster) > 1L) {
    not_available(gettext("if there is more than one cluster variable"))
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
    if (!ls_flag && !is.null(lavoptions$likelihood) &&
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
    # inequality constraints, or bounds: the active ones are treated as
    # (linear) equality constraints
    if (length(lavmodel@cin.nonlinear.idx) > 0L) {
      not_available(
        gettext("for models with nonlinear inequality constraints"))
    }
    active_idx <- setdiff(attr(lavmodel@con.jac, "cin.idx"),
                          attr(lavmodel@con.jac, "inactive.idx"))
    if (length(active_idx) > 0L) {
      if (length(lavmodel@ceq.nonlinear.idx) > 0L) {
        not_available(gettext(
          "for models with nonlinear equality constraints, if a bound or an
           inequality constraint is active"))
      }
      if (any(active_idx %in% attr(lavmodel@con.jac, "pd.idx"))) {
        not_available(gettext(
          "if the solution is on the boundary of the pd_lv parameter space"))
      }
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

# the (infinitesimal jackknife) influence values of the free parameters, per
# group: the rows of the g-th element correspond to the rows of
# lavdata@X[[g]], the columns to the free parameters; the change in the
# parameters for a set of weights w is sum_i (w_i - 1) * influence_i
lav_bootstrap_ij_influence <- function(lavmodel = NULL,
                                       lavsamplestats = NULL,
                                       lavdata = NULL,
                                       lavoptions = NULL,
                                       lavimplied = NULL,
                                       lavh1 = NULL,
                                       lavcache = NULL) {
  # least-squares estimators: the influence values of case i are the
  # (first-order) changes in the parameters if the weight of case i is
  # increased by one
  if (lavmodel@estimator != "ML") {
    x0 <- lav_model_get_parameters(lavmodel = lavmodel, type = "free")
    n_g <- lapply(lavdata@X, nrow)
    out <- vector("list", length = lavdata@ngroups)
    for (g in seq_len(lavdata@ngroups)) {
      freq <- lapply(n_g, function(n) matrix(1, n_g[[g]], n))
      freq[[g]] <- freq[[g]] + diag(n_g[[g]])
      out[[g]] <- lav_bootstrap_ij_ls(
        lavmodel = lavmodel, lavsamplestats = lavsamplestats,
        lavdata = lavdata, lavoptions = lavoptions, lavcache = lavcache,
        freq = freq, x0 = x0, method = "ij"
      )$d_theta
    }
    return(out)
  }
  sc <- lav_bootstrap_ij_scores(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavoptions = lavoptions, lavimplied = lavimplied
  )
  lavoptions$information[1] <- "observed"
  e_inv <- lav_model_info(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavcache = lavcache, lavimplied = lavimplied, lavh1 = lavh1,
    lavoptions = lavoptions, extra = FALSE, augmented = TRUE,
    inverted = TRUE
  )
  if (inherits(e_inv, "try-error") || anyNA(e_inv) ||
      ncol(e_inv) != ncol(sc[[1]])) {
    lav_msg_stop(gettext(
      "The inverse of the observed information matrix could not be
       computed."))
  }
  if (lavmodel@ceq.simple.only) {
    k_matrix <- lavmodel@ceq.simple.K
    e_inv <- e_inv %*% k_matrix %*% solve(crossprod(k_matrix))
  }
  lapply(sc, function(x) (x %*% e_inv) / lavsamplestats@ntotal)
}

# the bootstrap samples, as frequencies (weights), per group; the r rows
# are the bootstrap samples, the columns are the rows of lavdata@X[[g]]
#
# the samples are drawn in the same way (and in the same order) as in
# lav_bootstrap_internal(): for the same value of iseed, bootstrap sample b
# is identical to the one that is used for (serial) refit b
lav_bootstrap_ij_freq <- function(lavdata = NULL, r = 1000L, iseed = NULL) {
  ngroups <- lavdata@ngroups

  # cluster bootstrap: resample whole clusters (the number of clusters per
  # group is fixed, the number of observations is not)
  if (length(lavdata@cluster) > 0L) {
    cl_idx <- lapply(seq_len(ngroups), function(g) {
      lavdata@Lp[[g]]$cluster.idx[[2L]]
    })
    nclus <- lapply(cl_idx, function(x) length(unique(x)))
    boot_idx <- lav_bootstrap_indices(r = r, nobs = nclus, iseed = iseed)
    freq <- vector("list", length = ngroups)
    for (g in seq_len(ngroups)) {
      cell <- (boot_idx[[g]] - 1L) * r + as.vector(row(boot_idx[[g]]))
      freq_clus <- matrix(tabulate(cell, nbins = r * nclus[[g]]), r,
                          nclus[[g]])
      # the clusters, in the order used by lav_bootstrap_internal()
      cl_id <- as.integer(names(split(seq_along(cl_idx[[g]]), cl_idx[[g]])))
      freq[[g]] <- freq_clus[, match(cl_idx[[g]], cl_id), drop = FALSE]
    }
    attr(freq, "pool") <- lapply(lavdata@X, function(x) seq_len(nrow(x)))
    # the bootstrap samples as cluster indices
    attr(freq, "boot.idx") <- boot_idx
    return(freq)
  }

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
  if (length(lavdata@cluster) > 0L && is.null(freq)) {
    lav_msg_note(gettext(
      "clustered/multilevel data detected: the bootstrap resamples whole
       clusters (not individual observations)."))
  }

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

  # bootstrap weights; these can also be provided by the caller, as a list
  # (per group) of weight matrices (not necessarily integers; internal use)
  if (is.null(freq)) {
    freq <- lav_bootstrap_ij_freq(lavdata = lavdata, r = r, iseed = iseed)
  } else {
    stopifnot(is.list(freq), length(freq) == lavdata@ngroups)
    r <- nrow(freq[[1]])
    attr(freq, "pool") <- lapply(lavdata@X, function(x) seq_len(nrow(x)))
  }

  error_idx <- integer(0L)
  bound_idx <- integer(0L)
  guard_dist <- NULL
  if (lavmodel@estimator != "ML") {
    # least-squares estimators (continuous or categorical data)
    tmp <- lav_bootstrap_ij_ls(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavoptions = lavoptions, lavcache = lavcache, freq = freq, x0 = x0,
      method = method
    )
    d_theta <- tmp$d_theta
    error_idx <- tmp$error_idx
    guard_dist <- tmp$dist
    bound_idx <- tmp$bound_idx
  } else {
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
    if (method == "hoij2") {
      # simple equality constraints: everything in the compact parameter
      # space (theta.unco = K %*% theta, and e_inv = K (K'HK)^-1 K')
      h2_e_inv <- e_inv
      h2_g_delta <- g_delta
      if (lavmodel@ceq.simple.only) {
        k_matrix <- lavmodel@ceq.simple.K
        k_pinv <- solve(crossprod(k_matrix), t(k_matrix))
        h2_e_inv <- k_pinv %*% e_inv %*% t(k_pinv)
        h2_g_delta <- g_delta %*% k_matrix
      }
      tmp <- lav_bootstrap_hoij2(
        lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
        lavpartable = lavpartable, lavcache = lavcache, freq = freq, x0 = x0,
        g_delta = h2_g_delta, e_inv = h2_e_inv
      )
      d_theta <- tmp$d_theta + tmp$d2_theta
      error_idx <- tmp$error_idx
      guard_dist <- tmp$dist
      bound_idx <- tmp$bound_idx
    }

    # simple equality constraints: from the 'unco' space to the free parameters
    if (lavmodel@ceq.simple.only && method != "hoij2") {
      k_matrix <- lavmodel@ceq.simple.K
      d_theta <- d_theta %*% k_matrix %*% solve(crossprod(k_matrix))
    }
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
    refit_idx <- sort(unique(c(refit_idx, error_idx, bound_idx)))
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
    guard_info <- list(tol = guard_tol, refit.idx = refit_idx,
                       bound.idx = bound_idx)
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
# - g_delta: the (r x npar) weighted sums of the casewise scores,
#     sum_i (w_i - 1) s_i
# - e_inv: the (augmented) inverse of the observed (unit) information
#
# nonlinear equality constraints c(theta) = 0: the same reasoning is applied
# to the stationarity conditions of the Lagrangian,
#     F(theta, lambda, w) = [G(theta, w) + t(J(theta)) lambda; c(theta)] = 0
# with J the Jacobian of c(); the path is then a path for (theta, lambda),
# and e_inv is replaced by the inverse of the Jacobian of F, which also
# contains the curvature of the constraints
#
# returns a list with
# - d_theta: the first-order changes; if the model contains fixed
#     exogenous covariates (fixed.x = TRUE), their (co)variances and means
#     are fixed to the sample values, which change under reweighting; the
#     derivative of the gradient along the path is then computed by finite
#     differences (instead of using the scores) to take this into account
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
                                g_delta = NULL,
                                e_inv = NULL,
                                h = 0.05) {
  ngroups <- lavdata@ngroups
  r <- nrow(g_delta)
  npar <- length(x0)
  ntotal <- lavsamplestats@ntotal
  ceq_simple <- lavmodel@ceq.simple.only
  nonlinear <- length(lavmodel@ceq.nonlinear.idx) > 0L

  missing_flag <- lavsamplestats@missing.flag
  not_reproduced <- function() {
    lav_msg_stop(gettext(
      "The bootstrap method \"hoij2\" could not reproduce the sample
       statistics from the raw data; please use
       bootstrap = list(method = \"refit\") instead."))
  }
  differs <- function(x, y) {
    length(x) != length(y) || !all(is.finite(x)) ||
      max(abs(x - y)) > 1e-07 * max(1, abs(y))
  }

  # fixed exogenous covariates: where do their (fixed) sample statistics
  # live in GLIST?
  exo_map <- lav_bootstrap_hoij2_exo_map(
    lavmodel = lavmodel, lavdata = lavdata, lavpartable = lavpartable
  )
  exo_flag <- length(exo_map) > 0L

  # per group: the (weighted) sums that are needed to compute the sample
  # statistics, and their changes (per unit of t) for every bootstrap sample;
  # - complete data: the sums for all the variables
  # - missing data: the sums per missing pattern (and, if needed, the sums
  #   for the fixed exogenous covariates only)
  mom <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    y <- lavdata@X[[g]]
    wt <- lavdata@weights[[g]]
    if (is.null(wt)) {
      wt <- rep(1, nrow(y))
    }
    dw <- freq[[g]] - 1
    rows <- seq_len(nrow(y))
    cols <- seq_len(ncol(y))
    pat <- NULL
    if (missing_flag) {
      mp <- lavdata@Mp[[g]]
      if (length(mp$empty.idx) > 0L) {
        rows <- rows[-mp$empty.idx]
      }
      cols <- sort(unique(unlist(lapply(exo_map, function(el) {
        if (el$group == g) as.vector(el$idx) else integer(0L)
      }))))
      if (length(lavsamplestats@missing[[g]]) != mp$npatterns) {
        not_reproduced()
      }
      pat <- lapply(seq_len(mp$npatterns), function(p) {
        lav_bootstrap_hoij2_sums(y = y, wt = wt, dw = dw,
          rows = mp$case.idx[[p]], cols = which(mp$pat[p, ]))
      })
      # sanity check: can we reproduce the sample statistics?
      for (p in seq_len(mp$npatterns)) {
        st0 <- lav_bootstrap_hoij2_stats(pat[[p]], 1L, 0)
        yp <- lavsamplestats@missing[[g]][[p]]
        sy <- yp$SY
        if (length(sy) == 1L) {
          sy <- matrix(sy, pat[[p]]$nvar, pat[[p]]$nvar)
        }
        if (differs(st0$den, yp$freq) || differs(st0$mean, yp$MY) ||
            differs(st0$cov, sy)) {
          not_reproduced()
        }
      }
    }
    mom[[g]] <- list(
      nvar = ncol(y), cols = cols, pat = pat,
      full = lav_bootstrap_hoij2_sums(y = y, wt = wt, dw = dw, rows = rows,
                                      cols = cols)
    )
    if (!missing_flag) {
      st0 <- lav_bootstrap_hoij2_stats(mom[[g]]$full, 1L, 0)
      if (differs(st0$cov, lavsamplestats@cov[[g]]) ||
          differs(st0$mean, lavsamplestats@mean[[g]]) ||
          differs(st0$den, lavsamplestats@nobs[[g]])) {
        not_reproduced()
      }
    }
  }

  # the gradient, at theta, for bootstrap sample b, at position t on the path
  grad_fn <- function(theta, b, t) {
    cov_t <- mean_t <- vector("list", length = ngroups)
    for (g in seq_len(ngroups)) {
      mm <- mom[[g]]
      if (length(mm$cols) > 0L) {
        st <- lav_bootstrap_hoij2_stats(mm$full, b, t)
        mean_t[[g]] <- rep(as.numeric(NA), mm$nvar)
        mean_t[[g]][mm$cols] <- st$mean
        cov_t[[g]] <- matrix(as.numeric(NA), mm$nvar, mm$nvar)
        cov_t[[g]][mm$cols, mm$cols] <- st$cov
      }
      if (missing_flag) {
        nobs_g <- 0
        for (p in seq_along(mm$pat)) {
          st <- lav_bootstrap_hoij2_stats(mm$pat[[p]], b, t)
          lavsamplestats@missing[[g]][[p]]$SY <- st$cov
          lavsamplestats@missing[[g]][[p]]$MY <- st$mean
          lavsamplestats@missing[[g]][[p]]$freq <- st$den
          nobs_g <- nobs_g + st$den
        }
        lavsamplestats@nobs[[g]] <- nobs_g
      } else {
        lavsamplestats@mean[[g]] <- mean_t[[g]]
        lavsamplestats@cov[[g]] <- cov_t[[g]]
        lavsamplestats@nobs[[g]] <- st$den
      }
    }
    lavsamplestats@ntotal <- sum(unlist(lavsamplestats@nobs))
    glist <- lav_model_x2glist(lavmodel = lavmodel, x = theta)
    for (el in exo_map) {
      if (el$type == "cov") {
        value <- cov_t[[el$group]][el$idx]
      } else {
        value <- mean_t[[el$group]][el$idx]
      }
      glist[[el$mm]][el$m_idx] <- value
    }
    out <- try(suppressWarnings(lav_model_grad(
      lavmodel = lavmodel, glist = glist, lavsamplestats = lavsamplestats,
      lavdata = lavdata, lavcache = lavcache, type = "free",
      group_weight = !missing_flag, ceq_simple = ceq_simple
    )), silent = TRUE)
    if (inherits(out, "try-error") || length(out) != npar ||
        !all(is.finite(out))) {
      out <- rep(as.numeric(NA), npar)
    }
    as.numeric(out)
  }

  # sanity check: can we reproduce the (fixed) values of the exogenous
  # covariates in the model matrices?
  if (exo_flag) {
    for (g in seq_len(ngroups)) {
      st0 <- lav_bootstrap_hoij2_stats(mom[[g]]$full, 1L, 0)
      for (el in exo_map) {
        if (el$group != g) {
          next
        }
        if (el$type == "cov") {
          tmp <- matrix(0, mom[[g]]$nvar, mom[[g]]$nvar)
          tmp[mom[[g]]$cols, mom[[g]]$cols] <- st0$cov
        } else {
          tmp <- numeric(mom[[g]]$nvar)
          tmp[mom[[g]]$cols] <- st0$mean
        }
        if (differs(tmp[el$idx], lavmodel@GLIST[[el$mm]][el$m_idx])) {
          not_reproduced()
        }
      }
    }
  }

  # the derivative of the gradient along the path of the weights; if the
  # model contains fixed exogenous covariates, it is computed by finite
  # differences, to include the effect of the changing sample statistics
  # of these covariates
  g_dot <- NULL
  if (!exo_flag) {
    g_dot <- -1 * g_delta / ntotal
  }

  lav_bootstrap_hoij2_core(
    lavmodel = lavmodel, grad_fn = grad_fn, x0 = x0, r = r, e_inv = e_inv,
    g_dot = g_dot, ntotal = ntotal, h = h
  )
}

# the core of the second-order infinitesimal jackknife
#
# - grad_fn(theta, b, t): the gradient of the objective function at theta,
#     for bootstrap sample b, at position t on the path from the original
#     sample (t = 0) to the bootstrap sample (t = 1)
# - x0: the original estimates
# - e_inv: the (augmented) inverse of the Jacobian of the gradient
# - g_dot: (r x npar) matrix, the derivative of the gradient with respect to
#     t (at t = 0, and theta = x0); if NULL, finite differences are used
# - se: optional, the standard errors of the parameters; if provided, they
#     are used to express the distance to the exact solution in standard
#     error units (if not, the metric of the information matrix is used)
lav_bootstrap_hoij2_core <- function(lavmodel = NULL,
                                     grad_fn = NULL,
                                     x0 = NULL,
                                     r = 0L,
                                     e_inv = NULL,
                                     g_dot = NULL,
                                     ntotal = 1,
                                     se = NULL,
                                     h = 0.05) {
  npar <- length(x0)
  ceq_simple <- lavmodel@ceq.simple.only
  nonlinear <- length(lavmodel@ceq.nonlinear.idx) > 0L
  g_dot_all <- g_dot

  # nonlinear equality constraints: work with z = (theta, lambda)
  ncon <- 0L
  if (nonlinear) {
    con_fn <- lavmodel@ceq.function
    jac_fn <- function(theta) {
      lav_func_jacobian_complex(func = con_fn, x = theta)
    }
    ncon <- length(con_fn(x0))
    jac0 <- jac_fn(x0)
    # Lagrange multipliers at the solution
    lambda0 <- -1 * drop(solve(tcrossprod(jac0),
                               jac0 %*% grad_fn(x0, 1L, 0)))
    f_fn <- function(z, b, t) {
      theta <- z[seq_len(npar)]
      c(grad_fn(theta, b, t) + drop(crossprod(jac_fn(theta),
                                              z[npar + seq_len(ncon)])),
        con_fn(theta))
    }
    z0 <- c(x0, lambda0)
    # Jacobian of F with respect to z (central differences)
    nz <- npar + ncon
    a_mat <- matrix(0, nz, nz)
    for (k in seq_len(nz)) {
      step <- numeric(nz)
      step[k] <- 1e-05 * max(1, abs(z0[k]))
      a_mat[, k] <- (f_fn(z0 + step, 1L, 0) - f_fn(z0 - step, 1L, 0)) /
        (2 * step[k])
    }
    a_inv <- try(solve(a_mat), silent = TRUE)
    if (inherits(a_inv, "try-error") || anyNA(a_inv)) {
      lav_msg_stop(gettext(
        "The bootstrap method \"hoij2\" could not be computed for this
         model (singular system); please use
         bootstrap = list(method = \"refit\") instead."))
    }
    b_mat <- a_mat[seq_len(npar), seq_len(npar), drop = FALSE]
  } else {
    f_fn <- grad_fn
    z0 <- x0
    nz <- npar
    a_inv <- e_inv
  }
  pad <- numeric(ncon)

  # the (linear) equality constraints, and the active bounds and inequality
  # constraints (which are treated as equality constraints)
  # - a draw that violates an inactive bound/inequality constraint, or for
  #   which an active one is no longer active (the sign of its Lagrange
  #   multiplier changes), is flagged (bound_idx)
  con_jac <- lavmodel@con.jac
  cin_idx <- attr(con_jac, "cin.idx")
  inactive_idx <- attr(con_jac, "inactive.idx")
  cin_flag <- !nonlinear && length(cin_idx) > 0L &&
    !is.null(body(lavmodel@cin.function))
  h_mat <- matrix(0, 0L, npar)
  lambda_fn <- NULL
  if (!nonlinear && nrow(con_jac) > 0L) {
    keep_idx <- setdiff(seq_len(nrow(con_jac)), inactive_idx)
    h_mat <- unclass(con_jac)[keep_idx, , drop = FALSE]
    attributes(h_mat) <- list(dim = dim(h_mat))
    if (ceq_simple) {
      h_mat <- h_mat %*% lavmodel@ceq.simple.K
    }
    active_flag <- keep_idx %in% cin_idx
    nonzero_idx <- which(rowSums(abs(h_mat)) > 1e-12)
    h_mat <- h_mat[nonzero_idx, , drop = FALSE]
    active_flag <- active_flag[nonzero_idx]
    if (any(active_flag)) {
      # (least-squares) Lagrange multipliers of the active constraints, as
      # a function of the gradient
      l_mat <- (MASS::ginv(tcrossprod(h_mat)) %*% h_mat)[active_flag, ,
        drop = FALSE]
      lambda_fn <- function(gradient) drop(l_mat %*% gradient)
    }
  }

  # effective number of parameters
  npar_eff <- npar - ncon
  if (nrow(h_mat) > 0L) {
    npar_eff <- npar - qr(h_mat)$rank
  } else if (lavmodel@eq.constraints) {
    npar_eff <- ncol(lavmodel@eq.constraints.K)
  }
  npar_eff <- max(npar_eff, 1L)

  # (should be zero, unless there are active bounds/inequality constraints)
  r0 <- f_fn(z0, 1L, 0)
  lambda0 <- NULL
  if (!is.null(lambda_fn)) {
    lambda0 <- lambda_fn(r0)
  }
  bound_flag <- logical(r)

  d1_theta <- matrix(as.numeric(NA), r, npar)
  d2_theta <- matrix(as.numeric(NA), r, npar)
  dist <- rep(as.numeric(NA), r)
  for (b in seq_len(r)) {
    # derivative of the gradient along the path of the weights
    if (is.null(g_dot_all)) {
      g_dot <- (-grad_fn(x0, b, 2 * h) + 8 * grad_fn(x0, b, h) -
        8 * grad_fn(x0, b, -h) + grad_fn(x0, b, -2 * h)) / (12 * h)
    } else {
      g_dot <- g_dot_all[b, ]
    }
    # first-order change
    z_dot <- -1 * drop(a_inv %*% c(g_dot, pad))
    # second derivative of r(t) at t = 0 (five-point stencil)
    r2 <- (-f_fn(z0 + 2 * h * z_dot, b, 2 * h) +
      16 * f_fn(z0 + h * z_dot, b, h) - 30 * r0 +
      16 * f_fn(z0 - h * z_dot, b, -h) -
      f_fn(z0 - 2 * h * z_dot, b, -2 * h)) / (12 * h * h)
    z_ddot <- -1 * drop(a_inv %*% r2)
    d1_theta[b, ] <- z_dot[seq_len(npar)]
    d2_theta[b, ] <- 0.5 * z_ddot[seq_len(npar)]
    # how far are we from the exact solution for this bootstrap sample?
    # the gradient at the approximate solution (zero if it is exact),
    # expressed as the length of the corresponding Newton step, in the
    # metric of the (total) information: an estimate of the distance to
    # the exact solution in standard error units (root mean square over
    # the parameters)
    f1 <- f_fn(z0 + z_dot + 0.5 * z_ddot, b, 1)
    if (nonlinear) {
      step1 <- drop(a_inv %*% f1)[seq_len(npar)]
      if (!is.null(se)) {
        dist[b] <- sqrt(sum((step1[se > 0] / se[se > 0])^2) / npar_eff)
      } else {
        dist[b] <- sqrt(max(0, ntotal *
          sum(step1 * drop(b_mat %*% step1)) / npar_eff))
      }
    } else if (!is.null(se)) {
      step1 <- drop(e_inv %*% f1)[se > 0]
      dist[b] <- sqrt(sum((step1 / se[se > 0])^2) / npar_eff)
    } else {
      dist[b] <- sqrt(max(0, ntotal * sum(f1 * drop(e_inv %*% f1)) /
        npar_eff))
    }
    # bounds and inequality constraints
    if (cin_flag && all(is.finite(f1))) {
      theta_b <- x0 + d1_theta[b, ] + d2_theta[b, ]
      if (ceq_simple) {
        theta_b <- drop(lavmodel@ceq.simple.K %*% theta_b)
      }
      cin_b <- lavmodel@cin.function(theta_b)
      cin_inactive <- which(cin_idx %in% inactive_idx)
      if (any(cin_b[cin_inactive] < -1e-08)) {
        bound_flag[b] <- TRUE
      }
      if (!is.null(lambda_fn) && any(lambda_fn(f1) * sign(lambda0) <= 0)) {
        bound_flag[b] <- TRUE
      }
    }
  }

  error_idx <- which(!stats::complete.cases(d1_theta + d2_theta))

  list(d_theta = d1_theta, d2_theta = d2_theta, error_idx = error_idx,
       dist = dist, bound_idx = which(bound_flag))
}

# the (weighted) sums of a set of rows/columns of the data that are needed
# to compute the mean vector and the covariance matrix, and the changes in
# these sums (per unit of t) for every bootstrap sample (the rows of dw);
# the case weights at position t on the path are wt * (1 + t * dw)
lav_bootstrap_hoij2_sums <- function(y = NULL, wt = NULL, dw = NULL,
                                     rows = NULL, cols = NULL) {
  nvar <- length(cols)
  y <- y[rows, cols, drop = FALSE]
  if (anyNA(y)) {
    lav_msg_stop(gettext(
      "The bootstrap method \"hoij2\" is not available if the fixed
       exogenous covariates contain missing values; please use
       bootstrap = list(method = \"refit\") instead."))
  }
  vech_idx <- which(lower.tri(diag(nvar), diag = TRUE), arr.ind = TRUE)
  yy <- y[, vech_idx[, 1L], drop = FALSE] * y[, vech_idx[, 2L], drop = FALSE]
  wt <- wt[rows]
  dw <- dw[, rows, drop = FALSE]
  list(
    nvar = nvar, vech_idx = vech_idx,
    den = sum(wt), s1 = colSums(wt * y), s2 = colSums(wt * yy),
    d_den = drop(dw %*% wt), d_s1 = dw %*% (wt * y), d_s2 = dw %*% (wt * yy)
  )
}

# the sum of the weights, the mean vector and the ('ML') covariance matrix
# for bootstrap sample b, at position t on the path
lav_bootstrap_hoij2_stats <- function(sums = NULL, b = 1L, t = 0) {
  nvar <- sums$nvar
  den <- sums$den + t * sums$d_den[b]
  cov_t <- matrix(0, nvar, nvar)
  if (den < 1e-10) {
    # (almost) no observations left
    return(list(den = 0, mean = numeric(nvar), cov = cov_t))
  }
  mean_t <- (sums$s1 + t * sums$d_s1[b, ]) / den
  yy_t <- (sums$s2 + t * sums$d_s2[b, ]) / den
  cov_t[sums$vech_idx] <- yy_t
  cov_t[sums$vech_idx[, 2:1, drop = FALSE]] <- yy_t
  list(den = den, mean = mean_t, cov = cov_t - tcrossprod(mean_t))
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

# ------------------------------------------------------------------------
# least-squares estimators (DWLS, ULS), continuous or categorical data
#
# The estimator is a function of the sample statistics s (means and
# covariances, or thresholds and polychoric correlations) and of the weight
# matrix W: theta = argmin (s - sigma(theta))' W (s - sigma(theta)).
#
# - first stage: the sample statistics of a bootstrap sample are
#   approximated by a first-order expansion,
#       s(w) ~ s + sum_i (w_i - 1) infl_i
#   where infl_i are the casewise influence values of the statistics (their
#   crossproduct is the asymptotic covariance matrix of the statistics)
# - the weight matrix is kept fixed (for ULS, this is exact)
# - second stage, method = "ij": first-order expansion of theta(s), using
#   the observed information (the Jacobian of the gradient); the covariance
#   matrix of the draws is (up to Monte Carlo error) the robust (sandwich)
#   covariance matrix (se = "robust.sem") with information_bread = "observed"
# - second stage, method = "hoij2": second-order expansion of theta(s),
#   along the path s(t) = s + t * (s(w) - s)
#
# Note that "hoij2" is therefore NOT a complete second-order expansion in
# the weights for these estimators: the second-order terms of the first
# stage, and the changes in the weight matrix, are ignored.
# ------------------------------------------------------------------------

# the casewise influence values of the sample statistics (in the order of
# WLS.obs), per group
lav_bootstrap_ij_ls_influence <- function(lavmodel = NULL,
                                          lavsamplestats = NULL,
                                          lavdata = NULL,
                                          lavoptions = NULL) {
  ngroups <- lavdata@ngroups
  out <- vector("list", length = ngroups)
  for (g in seq_len(ngroups)) {
    y <- lavdata@X[[g]]
    n_g <- nrow(y)
    nvar <- ncol(y)
    if (anyNA(y)) {
      lav_msg_stop(gettext(
        "The approximate bootstrap methods need complete data if the
         estimator is DWLS or ULS."))
    }
    if (lavmodel@categorical) {
      ov_idx <- match(lavdata@ov.names[[g]], lavdata@ov$name)
      m84_args <- list(
        data_1 = y, ov_names = lavdata@ov.names[[g]],
        ov_types = lavdata@ov$type[ov_idx],
        ov_levels = lavdata@ov$nlev[ov_idx], wls_w = TRUE, group = g
      )
      if (!is.null(lavoptions$zero.add)) {
        m84_args$zero_add <- lavoptions$zero.add
      }
      if (!is.null(lavoptions$zero.keep.margins)) {
        m84_args$zero_keep_margins <- lavoptions$zero.keep.margins
      }
      if (!is.null(lavoptions$allow.empty.cell)) {
        m84_args$allow_empty_cell <- lavoptions$allow.empty.cell
      }
      cat_1 <- do.call(muthen1984, m84_args)
      b_inv <- lav_m84_b_inv(
        a11 = cat_1$A11, a21 = cat_1$A21, a22 = cat_1$A22
      )$b_inv
      infl <- cat_1$SC %*% t(b_inv) %*% t(cat_1$H)
      # numeric variables: WLS.obs contains -mu
      num_idx <- which(unlist(cat_1$TH.IDX) == 0L)
      if (length(num_idx) > 0L) {
        infl[, num_idx] <- -1 * infl[, num_idx]
      }
      infl[is.na(infl)] <- 0
      # sanity check: the crossproduct must reproduce Gamma
      nacov <- lavsamplestats@NACOV[[g]]
      if (!is.null(nacov) && length(nacov) > 0L &&
          length(lavdata@cluster) == 0L) {
        if (ncol(nacov) != ncol(infl) ||
            max(abs(n_g * crossprod(infl) - nacov)) >
              1e-06 * max(1, abs(nacov))) {
          lav_msg_stop(gettext(
            "The casewise influence values of the sample statistics could
             not be computed for this model; please use
             bootstrap = list(method = \"refit\") instead."))
        }
      }
    } else {
      yc <- t(t(y) - colMeans(y))
      idx1 <- lav_mat_vech_col_idx(nvar)
      idx2 <- lav_mat_vech_row_idx(nvar)
      zc <- yc[, idx1, drop = FALSE] * yc[, idx2, drop = FALSE]
      zc <- t(t(zc) - colMeans(zc))
      # is the divisor of the sample covariance matrix N or N - 1?
      cov_n <- crossprod(yc) / n_g
      ratio <- sum(diag(lavsamplestats@cov[[g]])) / sum(diag(cov_n))
      if (abs(ratio - n_g / (n_g - 1)) < abs(ratio - 1)) {
        zc <- zc * n_g / (n_g - 1)
      }
      if (lavmodel@meanstructure) {
        infl <- cbind(yc, zc) / n_g
      } else {
        infl <- zc / n_g
      }
    }
    if (ncol(infl) != length(lavsamplestats@WLS.obs[[g]])) {
      lav_msg_stop(gettext(
        "The casewise influence values of the sample statistics could not
         be computed for this model; please use
         bootstrap = list(method = \"refit\") instead."))
    }
    out[[g]] <- infl
  }

  out
}

# the changes in the free parameters, for the least-squares estimators
lav_bootstrap_ij_ls <- function(lavmodel = NULL,
                                lavsamplestats = NULL,
                                lavdata = NULL,
                                lavoptions = NULL,
                                lavcache = NULL,
                                freq = NULL,
                                x0 = NULL,
                                method = "ij",
                                h = 0.05) {
  ngroups <- lavdata@ngroups
  r <- nrow(freq[[1]])
  failed <- function() {
    lav_msg_stop(gettextf(
      "The bootstrap method %s could not be computed for this model; please
       use bootstrap = list(method = \"refit\") instead.", dQuote(method)))
  }

  # first stage: the changes in the sample statistics (first-order)
  infl <- lav_bootstrap_ij_ls_influence(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavoptions = lavoptions
  )
  s0 <- lavsamplestats@WLS.obs
  nstat <- lengths(s0)
  stat_idx <- lapply(seq_len(ngroups), function(g) {
    sum(nstat[seq_len(g - 1L)]) + seq_len(nstat[g])
  })
  d_s <- matrix(0, r, sum(nstat))
  for (g in seq_len(ngroups)) {
    d_s[, stat_idx[[g]]] <- (freq[[g]] - 1) %*% infl[[g]]
  }

  # the gradient, as a function of theta and the (stacked) changes in the
  # sample statistics; the weight matrix is fixed
  nvar <- lapply(lavdata@X, ncol)
  grad_fn <- function(theta, delta_s, ceq_simple = FALSE) {
    for (g in seq_len(ngroups)) {
      s_g <- s0[[g]] + delta_s[stat_idx[[g]]]
      lavsamplestats@WLS.obs[[g]] <- s_g
      if (!lavmodel@categorical) {
        # continuous data: the gradient uses the moments themselves
        if (lavmodel@meanstructure) {
          lavsamplestats@mean[[g]] <- s_g[seq_len(nvar[[g]])]
          s_g <- s_g[-seq_len(nvar[[g]])]
        }
        lavsamplestats@cov[[g]] <- lav_mat_vech_rev(s_g)
      }
    }
    out <- try(suppressWarnings(lav_model_grad(
      lavmodel = lavmodel,
      glist = lav_model_x2glist(lavmodel = lavmodel, x = theta),
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavcache = lavcache, type = "free", ceq_simple = ceq_simple
    )), silent = TRUE)
    if (inherits(out, "try-error") || !all(is.finite(out))) {
      out <- rep(as.numeric(NA), length(out))
    }
    as.numeric(out)
  }

  # the gradient is linear in the sample statistics: G = G0 + M %*% delta_s
  zero <- numeric(sum(nstat))
  g0 <- grad_fn(x0, zero)
  m_mat <- matrix(0, length(g0), sum(nstat))
  for (k in seq_len(sum(nstat))) {
    step <- zero
    step[k] <- 1
    m_mat[, k] <- grad_fn(x0, step) - g0
  }
  if (anyNA(m_mat)) {
    failed()
  }
  g_dot <- d_s %*% t(m_mat)

  # observed information (the Jacobian of the gradient with respect to the
  # parameters)
  npar <- length(x0)
  hess <- matrix(0, npar, npar)
  for (k in seq_len(npar)) {
    step <- numeric(npar)
    step[k] <- 1e-05 * max(1, abs(x0[k]))
    hess[, k] <- (grad_fn(x0 + step, zero) - grad_fn(x0 - step, zero)) /
      (2 * step[k])
  }
  if (anyNA(hess)) {
    failed()
  }
  hess <- (hess + t(hess)) / 2
  e_inv <- try(lav_model_info_augment_invert(
    lavmodel = lavmodel, information = hess, inverted = TRUE
  ), silent = TRUE)
  if (inherits(e_inv, "try-error") || anyNA(e_inv)) {
    failed()
  }
  if (method == "ij") {
    return(list(d_theta = -1 * g_dot %*% e_inv, error_idx = integer(0L),
                dist = NULL, bound_idx = integer(0L)))
  }

  # standard errors (of the first-order draws), to scale the distances
  se <- apply(g_dot %*% e_inv, 2L, stats::sd)

  out <- lav_bootstrap_hoij2_core(
    lavmodel = lavmodel,
    grad_fn = function(theta, b, t) grad_fn(theta, t * d_s[b, ]),
    x0 = x0, r = r, e_inv = e_inv, g_dot = g_dot,
    ntotal = lavsamplestats@ntotal, se = se, h = h
  )
  list(d_theta = out$d_theta + out$d2_theta, error_idx = out$error_idx,
       dist = out$dist, bound_idx = out$bound_idx)
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
