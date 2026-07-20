# delta-method standard errors for the noniterative estimators
#
# YR 20 Jul 2026: first version (wired for estimator = "MGM")
#
# Every noniterative estimator is a smooth map of the per-group sample
# moments, x = g(mean_b, vech(S_b)), so the covariance matrix of the
# estimates is the block-diagonal moment sandwich
#
#    vcov = sum_b J_b Gamma_b J_b' / nobs_b
#
# with J_b the Jacobian of the complete estimation map over block b's
# moments, and Gamma_b the asymptotic covariance of those moments. Two
# Gamma flavors are available:
#  - "nt" (default): the normal-theory Gamma;
#  - "adf": the distribution-free (empirical) Gamma. Because the casewise
#    influence of a moment-based estimator is J times the casewise moment
#    contribution, this flavor equals the infinitesimal-jackknife
#    covariance of the estimator.
#
# The Jacobian is computed numerically (numDeriv over the stacked moment
# vector, as in the JS fallback path), so ALL branches of the estimation
# map -- purging, pooled equality constraints, cross-loadings, the std.lv
# rescaling, the restricted PSI fit, the second-stage RLS weights -- are
# propagated automatically. The inner solvers converge to tolerances far
# below the numDeriv step size, so the finite differences stay clean; the
# regime switches (pd correction, clipping, bounds) are measure-zero
# kinks -- when a bound is ACTIVE at the estimate, the corresponding
# standard error degenerates towards zero.
#
# Other noniterative estimators (FABIN, BENTLER1982) can reuse the engine
# by adding a branch to lav_noniter_estimate_from_vec().

# re-run the complete estimation map on a (perturbed) moment vector;
# vec = NULL evaluates the map at the sample moments
lav_noniter_estimate_from_vec <- function(vec = NULL,
                                          lavmodel = NULL,
                                          lavsamplestats = NULL,
                                          lavpartable = NULL,
                                          lavdata = NULL,
                                          lavh1 = NULL,
                                          lavoptions = NULL) {
  if (!is.null(vec)) {
    implied <- lav_vec_to_implied(vec, lavmodel = lavmodel)
    # stage 1 reads lavsamplestats@cov/@mean; the second stage (and the
    # mean solve) read lavh1$implied when available: both must move
    # together
    for (b in seq_len(lavmodel@nblocks)) {
      lavsamplestats@cov[[b]] <- implied$cov[[b]]
      if (lavmodel@meanstructure) {
        lavsamplestats@mean[[b]] <- implied$mean[[b]]
      }
    }
    if (!is.null(lavh1$implied$cov)) {
      lavh1$implied$cov <- implied$cov
      if (lavmodel@meanstructure) {
        lavh1$implied$mean <- implied$mean
      }
    }
  }

  if (lavoptions$estimator == "MGM") {
    x <- lav_cfa_guttman1952_internal(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats,
      lavpartable = lavpartable, lavdata = lavdata,
      lavh1 = lavh1, lavoptions = lavoptions
    )
  } else {
    lav_msg_stop(gettextf(
      "delta-method standard errors are not available for estimator %s.",
      lavoptions$estimator))
  }

  as.numeric(x)
}

lav_noniter_vcov <- function(lavmodel = NULL, lavsamplestats = NULL,
                             lavoptions = NULL, lavpartable = NULL,
                             lavimplied = NULL, lavh1 = NULL,
                             lavdata = NULL) {
  # the engine handles continuous, unconditional moments only
  if (lavmodel@categorical || lavmodel@conditional.x ||
      lavmodel@correlation || lavmodel@group.w.free) {
    lav_msg_stop(gettextf(
      "delta-method standard errors for estimator %s require continuous,
       unconditional sample moments (no categorical, conditional.x,
       correlation or group.w.free structures).", lavoptions$estimator))
  }

  # Gamma flavor; NOTE: read with [[ ]] -- $ would partially match
  gamma_flavor <- lavoptions$estimator.args[["mgm.gamma"]]
  if (is.null(gamma_flavor)) {
    gamma_flavor <- "nt"
  }
  gamma_flavor <- tolower(gamma_flavor)
  if (gamma_flavor == "adf" && lavdata@data.type != "full") {
    lav_msg_warn(gettext(
      "the \"adf\" moment covariance requires complete raw data; using
       the normal-theory moment covariance instead."))
    gamma_flavor <- "nt"
  }

  # the sample moments (per block: mean if meanstructure, then vech(S)),
  # in the same layout that lav_vec_to_implied() unpacks
  implied0 <- lavh1$implied
  if (is.null(implied0$cov)) {
    implied0 <- list(cov = lavsamplestats@cov, mean = lavsamplestats@mean)
  }
  vec0 <- numeric(0L)
  for (b in seq_len(lavmodel@nblocks)) {
    if (lavmodel@meanstructure) {
      vec0 <- c(vec0, implied0$mean[[b]])
    }
    vec0 <- c(vec0, lav_mat_vech(implied0$cov[[b]]))
  }

  # numerical Jacobian of the complete estimation map over the stacked
  # moments; r = 2 Richardson extrapolation is accurate to ~1e-6 relative
  # and halves the number of map evaluations (the map itself is smooth);
  # warnings raised on perturbed evaluations (purge admissibility, RLS
  # iteration counts) merely repeat what the point estimation already
  # reported, so they are suppressed here
  jac <- suppressWarnings(numDeriv::jacobian(
    func = function(v) {
      lav_noniter_estimate_from_vec(
        vec = v, lavmodel = lavmodel, lavsamplestats = lavsamplestats,
        lavpartable = lavpartable, lavdata = lavdata,
        lavh1 = lavh1, lavoptions = lavoptions
      )
    },
    x = vec0, method.args = list(r = 2L)
  ))

  # vcov = sum over blocks of J_b Gamma_b J_b' / nobs_b; for the (default)
  # normal-theory Gamma the sandwich is computed directly from the sample
  # covariance matrix (lav_mat_k_gammant_kt with K = J_b), without ever
  # forming the pstar x pstar Gamma matrix
  npar <- nrow(jac)
  vcov <- matrix(0, nrow = npar, ncol = npar)
  offset <- 0L
  for (b in seq_len(lavmodel@nblocks)) {
    nvar_b <- lavmodel@nvar[b]
    nd_b <- if (lavmodel@meanstructure) {
      nvar_b + nvar_b * (nvar_b + 1L) / 2L
    } else {
      nvar_b * (nvar_b + 1L) / 2L
    }
    jac_b <- jac[, offset + seq_len(nd_b), drop = FALSE]
    offset <- offset + nd_b
    if (gamma_flavor == "adf") {
      gg <- lav_samp_gamma(
        m_y = lavdata@X[[b]],
        meanstructure = lavmodel@meanstructure
      )
      vcov <- vcov + (jac_b %*% gg %*% t(jac_b)) / lavsamplestats@nobs[[b]]
    } else {
      vcov <- vcov + lav_mat_k_gammant_kt(
        m_k = jac_b, s = implied0$cov[[b]],
        meanstructure = lavmodel@meanstructure,
        x_idx = integer(0L)
      ) / lavsamplestats@nobs[[b]]
    }
  }
  vcov <- (vcov + t(vcov)) / 2

  # the rest of lavaan expects the vcov of a ceq.simple.only model in the
  # 'unco' space (one row/column per non-collapsed parameter); the vcov
  # above is built in the compact (nx.free) space, so expand it via the
  # ceq.simple.K mapping (as for the IV and JS estimators)
  if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L) {
    vcov <- lavmodel@ceq.simple.K %*% vcov %*% t(lavmodel@ceq.simple.K)
  }

  vcov
}
