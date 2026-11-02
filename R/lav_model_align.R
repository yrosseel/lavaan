# alignment of a multiple-group confirmatory model (lavOptions alignment =
# "fixed" | "free"): the factor variances and means of the std.lv configural
# model (variances one, means zero in every group) are chosen to minimize
# the alignment loss over the loadings and the intercepts of the indicators
# (Asparouhov & Muthen, 2014); the model is then rescaled (H = diag(scales)
# per group) and translated by the aligned means, so that the model-implied
# moments are unchanged. No rotation is involved; the machinery is the
# scales-only mode of the agreement rotation (lav_mat_rotate_scales_mg()).
#
# YR 02 Nov 2026: first version

# the latent variables and indicators that are aligned, with the checks
lav_model_align_lvs <- function(lavmodel = NULL, lavpartable = NULL,
                                lavoptions = NULL) {
  if (lavmodel@ngroups < 2L) {
    lav_msg_stop(gettext("alignment requires a multiple-group model."))
  }
  if (lavmodel@nefa > 0L) {
    lav_msg_stop(gettext(
      "alignment is not available for models with efa blocks; use
       rotation.args = list(mg_agreement = TRUE, ...) instead."))
  }
  if (lavmodel@nblocks > lavmodel@ngroups) {
    lav_msg_stop(gettext("alignment is not available for multilevel models."))
  }
  if (lavmodel@categorical || lavmodel@conditional.x) {
    lav_msg_stop(gettext(
      "alignment is not available for categorical data or conditional.x =
       TRUE."))
  }
  if (!lavmodel@meanstructure) {
    lav_msg_stop(gettext("alignment requires a mean structure
                          (meanstructure = TRUE)."))
  }
  if (any(c("loadings", "intercepts") %in% lavoptions$group.equal)) {
    lav_msg_stop(gettext(
      "alignment requires a configural model: the loadings and intercepts
       must be free in every group (group.equal)."))
  }
  for (g in seq_len(lavmodel@ngroups)) {
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L ||
        length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L) {
      lav_msg_stop(gettext(
        "alignment is not available when observed variables appear as
         predictors or outcomes in regressions."))
    }
  }
  lv_names <- lav_pt_vnames(lavpartable, type = "lv.regular", block = 1L)
  if (length(lv_names) == 0L) {
    lav_msg_stop(gettext("alignment requires latent variables with
                          indicators."))
  }
  # configural identification: variances one, means zero, in every group
  var_idx <- which(lavpartable$op == "~~" & lavpartable$lhs == lavpartable$rhs &
                   lavpartable$lhs %in% lv_names)
  int_idx <- which(lavpartable$op == "~1" & lavpartable$lhs %in% lv_names)
  # the current values (the est column may not be filled in yet)
  est <- lav_model_get_parameters(lavmodel, type = "user", extra = TRUE)
  if (any(lavpartable$free[var_idx] > 0L) ||
      any(abs(est[var_idx] - 1) > 1e-12) ||
      any(lavpartable$free[int_idx] > 0L) ||
      any(abs(est[int_idx]) > 1e-12)) {
    lav_msg_stop(gettextf(
      "alignment requires a configural model with std.lv = TRUE: the factor
       variances fixed to one and the factor means fixed to zero in every
       group (free variances: %1$s, free means: %2$s, variances: %3$s,
       means: %4$s).",
      sum(lavpartable$free[var_idx] > 0L), sum(lavpartable$free[int_idx] > 0L),
      lav_msg_view(round(est[var_idx], 3), log_sep = "none"),
      lav_msg_view(round(est[int_idx], 3), log_sep = "none")))
  }
  ov_names <- unique(lavpartable$rhs[lavpartable$op == "=~" &
                                     lavpartable$lhs %in% lv_names])
  list(lv_names = lv_names, ov_names = ov_names)
}

# lower-level function, needed for numDeriv: the aligned parameters (type)
# as a function of the configural parameters x. init: list(scales, means)
# to start from (the solution, for the jacobian); extra = TRUE adds the
# attribute "extra" (GLIST, H, lv.order, rotation info)
lav_model_align_x <- function(x, lavmodel = NULL, lavoptions = NULL,
                              lavpartable = NULL, init = NULL,
                              extra = FALSE, type = "free", nobs = NULL) {
  aopts <- lavoptions$alignment.args
  intercepts <- lavoptions$alignment
  ngroups <- lavmodel@ngroups

  lavmodel <- lav_model_set_parameters(lavmodel, x = x)
  glist <- lavmodel@GLIST
  nmat <- lavmodel@nmat
  mm_list <- lapply(seq_len(ngroups), function(g) {
    seq_len(nmat[g]) + cumsum(c(0, nmat))[g]
  })

  sets <- lav_model_align_lvs(lavmodel, lavpartable, lavoptions)
  lv_names <- sets$lv_names
  ov_names <- sets$ov_names
  mlist_1 <- glist[mm_list[[1L]]]
  lv_all <- lav_model_efa_block_lv_names(lavmodel, mlist_1, 1L)
  ov_all <- lav_model_efa_block_ov_names(lavmodel, mlist_1, 1L)
  lv_idx <- match(lv_names, lv_all)
  ov_idx <- match(ov_names, ov_all)
  nlv <- length(lv_all)

  # the configural loadings and intercepts of the indicators
  a_list <- nu_list <- ov_var_list <- mu_list <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    mlist <- glist[mm_list[[g]]]
    a_list[[g]] <- mlist$lambda[ov_idx, lv_idx, drop = FALSE]
    nu_list[[g]] <- as.numeric(mlist$nu[ov_idx, 1L])
    if (isTRUE(aopts$std_ov)) {
      sigma_g <- mlist$lambda %*% mlist$psi %*% t(mlist$lambda) + mlist$theta
      ov_var_list[[g]] <- diag(sigma_g)[ov_idx]
      mu_list[[g]] <- nu_list[[g]]
    }
  }
  # group weights
  w_g <- NULL
  if (aopts$group_weights %in% c("nobs", "sqrt.nobs") ||
      isTRUE(aopts$std_ov)) {
    n_g <- unlist(nobs)
    if (length(n_g) != ngroups || any(!is.finite(n_g)) || any(n_g <= 0)) {
      lav_msg_stop(gettext("alignment needs the sample sizes of all groups."))
    }
  }
  if (aopts$group_weights == "nobs") {
    w_g <- n_g / mean(n_g)
  } else if (aopts$group_weights == "sqrt.nobs") {
    w_g <- sqrt(n_g / mean(n_g))
  }
  # standardization over the entire population
  if (isTRUE(aopts$std_ov)) {
    w_n <- n_g / sum(n_g)
    mu_pooled <- Reduce("+", Map("*", mu_list, w_n))
    var_pooled <- Reduce("+", Map(function(v, mu, w) {
      w * (v + (mu - mu_pooled)^2)
    }, ov_var_list, mu_list, w_n))
    sd_pooled <- sqrt(var_pooled)
    a_list <- lapply(a_list, function(x) x / sd_pooled)
    nu_list <- lapply(nu_list, function(x) x / sd_pooled)
  }

  # the scale split and the means: the scales-only engine of the agreement
  # rotation (no rotation criterion is involved: w = 1); "reference" is the
  # product solution re-expressed with variance one in the first group
  identification <- aopts$identification
  engine_identification <- if (identification == "reference") "product" else identification
  rot_list <- lav_mat_rotate_scales_mg(
    a_list = a_list, phi_list = NULL, orthogonal = FALSE,
    method_fname = lav_mat_rotate_quartimin, method_args = list(),
    mg_agreement_weight = 1, mg_agreement_method = aopts$method,
    mg_agreement_epsilon = aopts$epsilon, group_weights = w_g,
    scales = "free", identification = engine_identification,
    intercepts = intercepts, nu_list = nu_list,
    rstarts = if (is.null(init)) aopts$rstarts else 0L,
    c_init = init$scales, a_init = init$means
  )
  info <- attr(rot_list, "info")
  c_mat <- info$scales
  a_mat <- info$means
  if (identification == "reference") {
    c_first <- c_mat[, 1L]
    c_mat <- c_mat / c_first
    a_mat <- a_mat / c_first
  }

  # rescale (H = diag of the scales) and translate every group
  h <- vector("list", ngroups)
  rotation_info <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    mlist <- glist[mm_list[[g]]]
    hg <- diag(nlv)
    hg[lv_idx, lv_idx] <- diag(c_mat[, g], nrow = length(lv_idx))
    mlist <- lav_model_efa_apply_hg(mlist, hg)
    t_g <- numeric(nlv)
    t_g[lv_idx] <- a_mat[, g]
    if (any(t_g != 0)) {
      mlist <- lav_model_efa_apply_translation(mlist, t_g)
    }
    glist[mm_list[[g]]] <- mlist
    h[[g]] <- hg
    if (extra) {
      scales_g <- c_mat[, g]^2
      names(scales_g) <- lv_names
      means_g <- a_mat[, g]
      names(means_g) <- lv_names
      rotation_info[[g]] <- list(alignment = list(
        alignment = intercepts, lv_names = lv_names, ov_names = ov_names,
        method = aopts$method, epsilon = aopts$epsilon,
        group_weights = aopts$group_weights,
        identification = identification, std_ov = isTRUE(aopts$std_ov),
        converged = info$converged, criterion = info$method.value,
        variances = scales_g, means = means_g
      ))
    }
  }

  x_rot <- lav_model_get_parameters(lavmodel, glist = glist, type = type)
  if (extra) {
    attr(x_rot, "extra") <- list(
      GLIST = glist, H = h,
      lv.order = lapply(seq_len(ngroups), function(g) seq_len(nlv)),
      rotation = rotation_info,
      init = list(scales = info$scales, means = info$means)
    )
  }
  x_rot
}

# align the model: returns GLIST, H, lv.order, the rotation (alignment)
# information and the engine's solution (for warm starts)
lav_model_align <- function(lavmodel = NULL, lavoptions = NULL,
                            lavpartable = NULL, nobs = NULL) {
  x_orig <- lav_model_get_parameters(lavmodel, type = "free", extra = FALSE)
  tmp <- lav_model_align_x(
    x = x_orig, lavmodel = lavmodel, lavoptions = lavoptions,
    lavpartable = lavpartable, extra = TRUE, nobs = nobs
  )
  extra <- attr(tmp, "extra")
  list(GLIST = extra$GLIST, H = extra$H, lv.order = extra$lv.order,
       rotation = extra$rotation, init = extra$init)
}
