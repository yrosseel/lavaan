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
# - method = "hoij2": the second-order expansion (not available yet)
#
# The draws are returned in exactly the same format as the draws of
# lav_bootstrap_internal(fun = "coef"), so that everything downstream
# (standard errors, defined parameters, confidence intervals, ...) can be
# used unchanged.
#
# YR/LV 26 Nov 2026: first version (method = "ij")

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
  if (method != "ij") {
    lav_msg_stop(gettextf(
      "The bootstrap method %s is not available (yet); please use
       bootstrap = list(method = \"ij\") or
       bootstrap = list(method = \"refit\") instead.", dQuote(method)))
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
                             check_post = TRUE) {
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

  # bootstrap weights
  freq <- lav_bootstrap_ij_freq(lavdata = lavdata, r = r, iseed = iseed)

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

  # simple equality constraints: from the 'unco' space to the free parameters
  if (lavmodel@ceq.simple.only) {
    k_matrix <- lavmodel@ceq.simple.K
    d_theta <- d_theta %*% k_matrix %*% solve(crossprod(k_matrix))
  }

  coef_boot <- d_theta + matrix(x0, r, length(x0), byrow = TRUE)

  # same attributes as lav_bootstrap_internal()
  attr(coef_boot, "error.idx") <- integer(0L)
  if (check_post) {
    attr(coef_boot, "nonadmissible") <-
      lav_bootstrap_ij_nonadmissible(coef_boot, lavpartable = lavpartable)
  }
  attr(coef_boot, "seed") <- iseed
  attr(coef_boot, "method") <- method

  coef_boot
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
