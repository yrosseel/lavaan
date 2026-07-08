# two-level WLS/DWLS/ULS: sample 'statistics' and their asymptotic
# covariance matrix
#
# YR 2026 (two-level DWLS, phase 1: continuous, complete data)
#
# For two-level data the (per-group) unrestricted statistic vector is the
# saturated (h1) model estimate:
#
#     s = c( Mu.W, vech(Sigma.W), Mu.B, vech(Sigma.B) )
#
# in the same order as (1) the h1 parameter vector used by the
# lav_mvn_cl_* information/score functions, and (2) the per-group Delta
# matrix (lav_model_delta() rbinds the within and between blocks, each
# block laid out as means first, then vech of the covariance matrix).
#
# Note that s contains structurally-fixed entries: the within-level means
# of variables that live at both levels are zero by convention (their mean
# goes to Mu.B), and the model fixes the corresponding intercepts to zero
# as well. These entries have zero rows/columns in Gamma and in Delta, and
# get zero weight; their residuals are identically zero.
#
# The asymptotic covariance matrix of s treats the CLUSTERS as the
# independent units (sandwich):
#
#     acov(s) = (1/J) * I1^{-1} J1 I1^{-1}
#
# with I1 the (per-cluster, expected) h1 information and J1 the first-order
# information from the cluster-wise scores. To keep all downstream
# consumers (which scale by the total number of OBSERVATIONS) working
# unchanged, we store
#
#     NACOV = Gamma = nobs * acov(s)
#
# so that acov(s) = Gamma/nobs exactly, whatever the cluster structure.

# Gamma (= NACOV = nobs * acov(s)) for ONE group of a two-level dataset,
# evaluated at the h1 (saturated) estimates; the shared kernel for
# lav_samp_wls_2l() below and the post-fit lav_object_gamma()
lav_samp_gamma_2l_g <- function(lavsamplestats = NULL,
                                lavh1 = NULL,
                                lavdata = NULL,
                                g = 1L) {
  nlevels <- lavdata@nlevels
  lp <- lavdata@Lp[[g]]
  nobs <- lavsamplestats@nobs[[g]]
  nclusters <- lp$nclusters[[2]]

  # no fixed.x exogenous covariates (yet)
  if (length(lavsamplestats@x.idx[[g]]) > 0L) {
    lav_msg_stop(gettext(
      "fixed.x = TRUE is not supported (yet) for the two-level Gamma;
      use fixed.x = FALSE."))
  }

  # h1 (saturated) estimates for this group
  mu_w <- lavh1$implied$mean[[(g - 1) * nlevels + 1L]]
  sigma_w <- lavh1$implied$cov[[(g - 1) * nlevels + 1L]]
  mu_b <- lavh1$implied$mean[[(g - 1) * nlevels + 2L]]
  sigma_b <- lavh1$implied$cov[[(g - 1) * nlevels + 2L]]

  # expected h1 information (per cluster) at the h1 estimates
  i1 <- lav_mvn_cl_info_expected(
    lp = lp,
    mu_w = mu_w, sigma_w = sigma_w,
    mu_b = mu_b, sigma_b = sigma_b,
    x_idx = lavsamplestats@x.idx[[g]]
  )

  # first-order information (per cluster) from the cluster-wise scores
  j1 <- lav_mvn_cl_info_firstorder(
    y1 = lavdata@X[[g]],
    ylp = lavsamplestats@YLp[[g]],
    lp = lp,
    mu_w = mu_w, sigma_w = sigma_w,
    mu_b = mu_b, sigma_b = sigma_b,
    x_idx = lavsamplestats@x.idx[[g]],
    divide_by_two = TRUE
  )

  # invert I1, dropping the structurally-fixed (zero-information)
  # entries (cfr. lav_model_h1_acov)
  diag_i1 <- diag(i1)
  keep <- which(abs(diag_i1) > sqrt(.Machine$double.eps) * max(abs(diag_i1)))
  i1_inv <- matrix(0, nrow(i1), ncol(i1))
  inv_keep <- try(
    lav_mat_sym_inverse(i1[keep, keep, drop = FALSE]),
    silent = TRUE
  )
  if (inherits(inv_keep, "try-error")) {
    lav_msg_stop(gettextf(
      "could not invert the h1 information matrix in group %s", g))
  }
  i1_inv[keep, keep] <- inv_keep

  # sandwich acov of s (clusters are the units), rescaled to the
  # 'Gamma/nobs' convention
  acov <- (i1_inv %*% j1 %*% i1_inv) / nclusters
  gamma_1 <- nobs * acov
  # enforce symmetry (numerical)
  gamma_1 <- (gamma_1 + t(gamma_1)) / 2

  attr(gamma_1, "keep") <- keep
  gamma_1
}

# fill the WLS.obs/WLS.V/WLS.VD/NACOV slots of lavsamplestats for
# two-level data; called after the h1 model has been estimated (step 06)
lav_samp_wls_2l <- function(lavsamplestats = NULL,
                            lavh1 = NULL,
                            lavdata = NULL,
                            lavoptions = NULL) {
  stopifnot(lavdata@nlevels == 2L)
  estimator <- lavoptions$estimator
  ngroups <- lavsamplestats@ngroups
  nlevels <- lavdata@nlevels

  wls_obs <- vector("list", length = ngroups)
  wls_v <- vector("list", length = ngroups)
  wls_vd <- vector("list", length = ngroups)
  nacov <- vector("list", length = ngroups)

  for (g in seq_len(ngroups)) {
    lp <- lavdata@Lp[[g]]
    nclusters <- lp$nclusters[[2]]

    # h1 (saturated) estimates for this group
    mu_w <- lavh1$implied$mean[[(g - 1) * nlevels + 1L]]
    sigma_w <- lavh1$implied$cov[[(g - 1) * nlevels + 1L]]
    mu_b <- lavh1$implied$mean[[(g - 1) * nlevels + 2L]]
    sigma_b <- lavh1$implied$cov[[(g - 1) * nlevels + 2L]]

    # the statistic vector
    wls_obs[[g]] <- c(
      mu_w, lav_mat_vech(sigma_w, diagonal = TRUE),
      mu_b, lav_mat_vech(sigma_b, diagonal = TRUE)
    )

    # Gamma (cluster sandwich at the h1 estimates; shared kernel)
    gamma_1 <- lav_samp_gamma_2l_g(
      lavsamplestats = lavsamplestats,
      lavh1 = lavh1,
      lavdata = lavdata,
      g = g
    )
    keep <- attr(gamma_1, "keep")
    attr(gamma_1, "keep") <- NULL
    nacov[[g]] <- gamma_1

    # weight matrices
    if (estimator == "WLS") {
      if (nclusters < length(keep)) {
        lav_msg_stop(gettextf(
          "estimator WLS for two-level data needs more clusters (%1$s) than
          sample statistics (%2$s) in group %3$s; use estimator WLSMV
          (or ULSMV) instead.", nclusters, length(keep), g))
      }
      w_full <- matrix(0, nrow(gamma_1), ncol(gamma_1))
      w_keep <- try(
        lav_mat_sym_inverse(gamma_1[keep, keep, drop = FALSE]),
        silent = TRUE
      )
      if (inherits(w_keep, "try-error")) {
        lav_msg_stop(gettextf(
          "could not invert Gamma (acov of the two-level sample statistics)
          in group %s; use estimator WLSMV (or ULSMV) instead.", g))
      }
      w_full[keep, keep] <- w_keep
      wls_v[[g]] <- w_full
    } else if (estimator == "DWLS") {
      dacov <- diag(gamma_1)
      idacov <- ifelse(dacov > 0, 1 / dacov, 0)
      wls_v[[g]] <- diag(idacov, nrow = length(idacov))
      wls_vd[[g]] <- idacov
    } else if (estimator == "ULS") {
      wls_v[[g]] <- diag(length(wls_obs[[g]]))
      wls_vd[[g]] <- rep(1, length(wls_obs[[g]]))
    } else {
      lav_msg_stop(gettextf(
        "estimator %s is not supported for two-level WLS estimation.",
        estimator))
    }
  }

  # fill the slots
  lavsamplestats@WLS.obs <- wls_obs
  lavsamplestats@WLS.V <- wls_v
  lavsamplestats@WLS.VD <- wls_vd
  if (!lavsamplestats@NACOV.user) {
    lavsamplestats@NACOV <- nacov
  }

  lavsamplestats
}
