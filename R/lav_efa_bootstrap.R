# Bootstrap inference for rotated EFA/ESEM models.
#
# YR/Claude - 2026
#
# Ordinary (or other) bootstrap inference is not straightforward for models
# with rotated exploratory factors: the optimizer works with the *unrotated*
# solution (with echelon constraints), so the raw bootstrap draws (optim$x) do
# not correspond to the rotated parameters reported to the user. Moreover, even
# the rotated solution is only identified up to a signed permutation of the
# factors, so the rotated factors of each bootstrap sample must be aligned
# (re-ordered and reflected) to the original (full-sample) rotated solution
# before their variability is meaningful.
#
# This file implements that machinery. It is driven from
# lav_lavaan_step16_rotation.R when se = "bootstrap" and the model contains
# (rotated) exploratory factors. See github issue #376.

# extract, per block and per EFA set, the (rotated) loading submatrix involved
# in that set. Returns a list (per block) of lists (per set); entries for
# single-factor sets (or empty sets) are NULL (nothing to align there).
lav_efa_bootstrap_lambda_list <- function(lavmodel = NULL, glist = NULL) {
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }
  out <- vector("list", lavmodel@nblocks)
  for (g in seq_len(lavmodel@nblocks)) {
    mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
    mlist <- glist[mm_in_group]
    set_list <- vector("list", lavmodel@nefa)
    for (set in seq_len(lavmodel@nefa)) {
      ov_idx <- lavmodel@ov.efa.idx[[g]][[set]]
      lv_idx <- lavmodel@lv.efa.idx[[g]][[set]]
      if (length(ov_idx) == 0L || length(lv_idx) < 2L) {
        next
      }
      set_list[[set]] <- mlist$lambda[ov_idx, lv_idx, drop = FALSE]
    }
    out[[g]] <- set_list
  }
  out
}

# given a (rotated) fitted bootstrap model, align its rotated EFA factors to the
# reference solution 'ref_lambda' (a lambda-list as returned by
# lav_efa_bootstrap_lambda_list()), and return the aligned rotated parameter
# vector (free parameters only, in parameter-table order).
#
# Alignment is a signed permutation T per EFA set; it is applied to *all*
# factor-related model matrices (lambda, psi, beta, alpha) exactly as a rotation
# would be (see lav_model_efa_rotate_x()), so that every parameter connected to
# the rotated factors (loadings, factor (co)variances, structural coefficients,
# factor means) is transformed consistently. This makes the approach valid for
# ESEM (factors embedded in a larger structural model), not just plain EFA.
lav_efa_bootstrap_align_coef <- function(lavmodel = NULL, lavpartable = NULL,
                                         ref_lambda = NULL) {
  glist <- lavmodel@GLIST
  for (g in seq_len(lavmodel@nblocks)) {
    mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
    mlist <- glist[mm_in_group]

    # build the block-level signed permutation hg (identity outside the EFA
    # factors, so only the exploratory factors are touched)
    hg <- diag(ncol(mlist$lambda))
    changed <- FALSE
    for (set in seq_len(lavmodel@nefa)) {
      ref <- ref_lambda[[g]][[set]]
      if (is.null(ref)) {
        next
      }
      ov_idx <- lavmodel@ov.efa.idx[[g]][[set]]
      lv_idx <- lavmodel@lv.efa.idx[[g]][[set]]
      lambda_b <- mlist$lambda[ov_idx, lv_idx, drop = FALSE]
      tt <- lav_efa_align_signed_perm(lambda = lambda_b, lambda_ref = ref)
      hg[lv_idx, lv_idx] <- tt
      changed <- TRUE
    }

    if (changed) {
      # hg is an orthonormal signed permutation (solve(t(hg)) == hg), so the
      # general rotation transform reduces to lambda %*% hg etc. Reuse the same
      # helper as lav_model_efa_rotate_x() for a single, consistent transform.
      mlist <- lav_model_efa_apply_hg(mlist, hg)
      glist[mm_in_group] <- mlist
    }
  } # block

  lavmodel@GLIST <- glist
  est_user <- lav_model_get_parameters(
    lavmodel = lavmodel, type = "user", extra = TRUE
  )
  free_idx <- which(lavpartable$free > 0L & !duplicated(lavpartable$free))
  est_user[free_idx]
}

# run the (ordinary) bootstrap for a rotated EFA/ESEM model.
#
# 'lavmodel'/'lavpartable' are the *unrotated* model and parameter table (so
# each bootstrap sample is refit and rotated from scratch), while 'efa_ref'
# carries the reference rotated solution and the rotated free-parameter layout
# (see lav_efa_bootstrap_reference()). Returns the (r x npar.rotated) matrix of
# aligned rotated bootstrap coefficients, with the usual error.idx / seed /
# nonadmissible attributes.
lav_efa_bootstrap_run <- function(lavmodel = NULL, lavpartable = NULL,
                                  lavsamplestats = NULL, lavoptions = NULL,
                                  lavdata = NULL, efa_ref = NULL) {
  # number of bootstrap draws
  r <- lavoptions$bootstrap
  if (is.list(r)) {
    r <- r$R
  }
  if (is.null(r)) {
    r <- 1000L
  }
  r <- as.integer(r)

  coef_boot <- lav_bootstrap_internal(
    object = NULL,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavpartable = lavpartable,
    lavoptions = lavoptions,
    lavdata = lavdata,
    r = r,
    check_post = lavoptions$check.post,
    type = "ordinary",
    fun = "coef",
    efa_reference = efa_ref
  )

  coef_boot
}

# turn the (aligned, rotated) bootstrap coefficient matrix into a vcov + se,
# using the *rotated* model and parameter table. Returns a list with elements
# 'se' (for lavpartable$se), 'vcov' (the rotated bootstrap vcov) and 'coef' (the
# raw bootstrap coefficients, for the @boot slot). Mirrors the scaling used in
# lav_model_nvcov_bootstrap().
lav_efa_bootstrap_vcov <- function(coef_boot = NULL, lavmodel = NULL,
                                   lavpartable = NULL, lavoptions = NULL) {
  coef_orig <- coef_boot

  # drop failed runs
  error_idx <- attr(coef_boot, "error.idx")
  if (length(error_idx) > 0L) {
    coef_boot <- coef_boot[-error_idx, , drop = FALSE] # drops attributes
  }

  nboot <- nrow(coef_boot)
  if (is.null(nboot) || nboot < 2L) {
    lav_msg_stop(gettext(
      "not enough successful bootstrap runs to compute standard errors."))
  }

  # the bootstrap directly estimates Cov(coef), which IS the target vcov of the
  # (rotated) parameters. Use the same (nboot-1)/nboot scaling as the rest of
  # lavaan (see lav_model_nvcov_bootstrap()).
  vcov <- cov(coef_boot) * (nboot - 1) / nboot

  # standard errors for the (rotated) free parameters + defined parameters
  se <- lav_model_vcov_se(
    lavmodel = lavmodel, lavpartable = lavpartable,
    vcov = vcov, boot = coef_orig, lavoptions = lavoptions
  )

  list(se = se, vcov = vcov, coef = coef_orig)
}
