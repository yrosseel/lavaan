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
  if (lavmodel@nblocks > lavmodel@ngroups) {
    lav_msg_stop(gettext("alignment is not available for multilevel models."))
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
  # the aligned factors: the regular latent variables that are not part
  # of an efa block (the efa factors are handled by the agreement rotation)
  lv_names <- lav_pt_vnames(lavpartable, type = "lv.regular", block = 1L)
  if (lavmodel@nefa > 0L) {
    lv_efa <- unique(lavpartable$lhs[lavpartable$op == "=~" &
      !is.na(lavpartable$efa) & nzchar(lavpartable$efa)])
    lv_names <- lv_names[!lv_names %in% lv_efa]
  }
  if (length(lv_names) == 0L) {
    lav_msg_stop(gettext("alignment requires (confirmatory) latent
                          variables with indicators."))
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
  # the indicators must not be predictors or outcomes in regressions (their
  # intercepts would live in alpha of a dummy latent variable); other
  # observed variables may (exogenous covariates, outcomes)
  reg_ov <- unique(c(lavpartable$lhs[lavpartable$op == "~"],
                     lavpartable$rhs[lavpartable$op == "~"]))
  if (any(ov_names %in% reg_ov)) {
    lav_msg_stop(gettextf(
      "alignment is not available when an indicator of an aligned factor is
       also a predictor or an outcome in a regression (%s).",
      lav_msg_view(ov_names[ov_names %in% reg_ov], log_sep = "none")))
  }
  list(lv_names = lv_names, ov_names = ov_names)
}

# lower-level function, needed for numDeriv: the aligned parameters (type)
# as a function of the configural parameters x. init: list(scales, means)
# to start from (the solution, for the jacobian); extra = TRUE adds the
# attribute "extra" (GLIST, H, lv.order, rotation info). If the model has
# efa blocks, they are rotated first (lav_model_efa_rotate_x(), with the
# agreement rotation; init_rot warm-starts it for the jacobian) and the
# confirmatory factors are aligned on the rotated model: the two
# transformations act on disjoint factors and compose
lav_model_align_x <- function(x, lavmodel = NULL, lavoptions = NULL,
                              lavpartable = NULL, init = NULL,
                              init_rot = NULL, extra = FALSE,
                              type = "free", nobs = NULL) {
  aopts <- lavoptions$alignment.args
  intercepts <- lavoptions$alignment
  ngroups <- lavmodel@ngroups

  lavmodel <- lav_model_set_parameters(lavmodel, x = x)
  h_efa <- NULL
  rotation_efa <- NULL
  lv_order <- NULL
  if (lavmodel@nefa > 0L && lavoptions$rotation != "none") {
    rot <- lav_model_efa_rotate_x(
      x = x, lavmodel = lavmodel, lavoptions = lavoptions,
      init_rot = init_rot, extra = TRUE, type = "free", nobs = nobs
    )
    rot_extra <- attr(rot, "extra")
    glist <- rot_extra$GLIST
    h_efa <- rot_extra$H
    lv_order <- rot_extra$lv.order
    rotation_efa <- rot_extra$rotation
  } else {
    glist <- lavmodel@GLIST
  }
  nmat <- lavmodel@nmat
  mm_list <- lapply(seq_len(ngroups), function(g) {
    seq_len(nmat[g]) + cumsum(c(0, nmat))[g]
  })

  sets <- lav_model_align_lvs(lavmodel, lavpartable, lavoptions)
  lv_names <- sets$lv_names
  ov_names <- sets$ov_names
  mlist_1 <- glist[mm_list[[1L]]]
  # all latent variables of the block, the dummy ones (covariates,
  # observed outcomes) included: H has that dimension
  lambda_mm <- mm_list[[1L]][which(names(mlist_1) == "lambda")[1L]]
  lv_all <- lavmodel@dimNames[[lambda_mm]][[2L]]
  ov_all <- lav_model_efa_block_ov_names(lavmodel, mlist_1, 1L)
  lv_idx <- match(lv_names, lv_all)
  ov_idx <- match(ov_names, ov_all)
  nlv <- length(lv_all)

  # the configural loadings and intercepts of the indicators (for ordered
  # indicators: the thresholds, see lav_model_efa_intercepts())
  a_list <- nu_list <- ov_var_list <- mu_list <- vector("list", ngroups)
  nu_row_idx <- NULL
  for (g in seq_len(ngroups)) {
    mlist <- glist[mm_list[[g]]]
    a_list[[g]] <- mlist$lambda[ov_idx, lv_idx, drop = FALSE]
    nu_g <- lav_model_efa_intercepts(mlist, ov_idx,
                                     th_idx = lavmodel@th.idx[[g]])
    nu_list[[g]] <- nu_g$values
    nu_row_idx <- nu_g$row_idx # assumed equal across groups!
    if (isTRUE(aopts$std_ov)) {
      # the model-implied variances (of the latent responses, for ordered
      # indicators) and means
      sigma_g <- lav_lisrel_sigma(mlist = mlist, delta = TRUE)
      ov_var_list[[g]] <- diag(sigma_g)[ov_idx]
      mu_list[[g]] <- as.numeric(mlist$nu[ov_idx, 1L])
    }
  }
  if (isTRUE(aopts$std_ov) && lavmodel@conditional.x) {
    lav_msg_stop(gettext(
      "alignment.args = list(std_ov = TRUE) is not available with
       conditional.x = TRUE."))
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
    nu_list <- lapply(nu_list, function(x) x / sd_pooled[nu_row_idx])
  }

  # the scale split and the means: the scales-only engine of the agreement
  # rotation (no rotation criterion is involved: w = 1); "reference" is the
  # product solution re-expressed with variance one in the first group
  identification <- aopts$identification
  engine_identification <- if (identification == "reference") "product" else identification
  # the reference group (variance one under "reference", mean zero under
  # "fixed"): the engine takes the first group as the reference, so the
  # groups are permuted; "auto" (fixed alignment, not inside the jacobian)
  # tries every group and keeps the one with the smallest loss, i.e. the
  # most invariant aligned model
  ref <- aopts$reference_group
  if (identical(ref, "auto")) {
    candidates <- if (intercepts == "fixed" && is.null(init)) {
      seq_len(ngroups)
    } else {
      1L
    }
  } else {
    candidates <- lav_align_reference_group(ref, ngroups)
  }
  run_engine <- function(ref_g) {
    ord <- c(ref_g, seq_len(ngroups)[-ref_g])
    out <- lav_mat_rotate_scales_mg(
      a_list = a_list[ord], phi_list = NULL, orthogonal = FALSE,
      method_fname = lav_mat_rotate_quartimin, method_args = list(),
      mg_agreement_weight = 1, mg_agreement_method = aopts$method,
      mg_agreement_epsilon = aopts$epsilon,
      group_weights = if (is.null(w_g)) NULL else w_g[ord],
      scales = "free", identification = engine_identification,
      intercepts = intercepts, nu_list = nu_list[ord],
      nu_row_idx = nu_row_idx,
      rstarts = if (is.null(init)) aopts$rstarts else 0L,
      c_init = if (is.null(init$scales)) NULL else init$scales[, ord, drop = FALSE],
      a_init = if (is.null(init$means)) NULL else init$means[, ord, drop = FALSE]
    )
    info_g <- attr(out, "info")
    inv_ord <- order(ord)
    info_g$scales <- info_g$scales[, inv_ord, drop = FALSE]
    if (!is.null(info_g$means)) {
      info_g$means <- info_g$means[, inv_ord, drop = FALSE]
    }
    info_g
  }
  info <- NULL
  for (ref_g in candidates) {
    info_g <- run_engine(ref_g)
    if (is.null(info) || isTRUE(info_g$method.value < info$method.value)) {
      info <- info_g
      ref <- ref_g
    }
  }
  c_mat <- info$scales
  a_mat <- info$means
  if (identification == "reference") {
    c_ref <- c_mat[, ref]
    c_mat <- c_mat / c_ref
    a_mat <- a_mat / c_ref
  }

  # rescale (H = diag of the scales) and translate every group; with efa
  # blocks the rescaling composes with their rotation (disjoint factors)
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
      mlist <- lav_model_efa_apply_translation(mlist, t_g,
                                               th_idx = lavmodel@th.idx[[g]])
    }
    glist[mm_list[[g]]] <- mlist
    h[[g]] <- if (is.null(h_efa)) hg else h_efa[[g]] %*% hg
    if (extra && !is.null(rotation_efa)) {
      rotation_info[[g]] <- rotation_efa[[g]]
    }
    if (extra) {
      scales_g <- c_mat[, g]^2
      names(scales_g) <- lv_names
      means_g <- a_mat[, g]
      names(means_g) <- lv_names
      rotation_info[[g]]$alignment <- list(
        alignment = intercepts, lv_names = lv_names, ov_names = ov_names,
        method = aopts$method, epsilon = aopts$epsilon,
        group_weights = aopts$group_weights,
        identification = identification, std_ov = isTRUE(aopts$std_ov),
        reference_group = ref,
        converged = info$converged, criterion = info$method.value,
        variances = scales_g, means = means_g,
        starts = info$starts
      )
    }
  }

  x_rot <- lav_model_get_parameters(lavmodel, glist = glist, type = type)
  if (extra) {
    if (is.null(lv_order)) {
      lv_order <- lapply(seq_len(ngroups), function(g) seq_len(nlv))
    }
    attr(x_rot, "extra") <- list(
      GLIST = glist, H = h, H_efa = h_efa, lv.order = lv_order,
      rotation = rotation_info,
      init = list(scales = info$scales, means = info$means),
      reference_group = ref
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
  list(GLIST = extra$GLIST, H = extra$H, H_efa = extra$H_efa,
       lv.order = extra$lv.order, rotation = extra$rotation,
       init = extra$init, reference_group = extra$reference_group)
}

# the reference group of an alignment: a group number or a group label
# (resolved against group_label); "auto" is passed through if allowed
lav_align_reference_group <- function(ref = 1L, ngroups = 1L,
                                      group_label = NULL,
                                      allow_auto = FALSE) {
  if (is.null(ref)) {
    return(1L)
  }
  if (allow_auto && identical(ref, "auto")) {
    return("auto")
  }
  if (is.character(ref)) {
    if (length(ref) != 1L || is.null(group_label) ||
        !(ref %in% group_label)) {
      lav_msg_stop(gettextf(
        "reference_group %s is not a group label (or \"auto\").",
        dQuote(ref)))
    }
    ref <- match(ref, group_label)
  }
  ref <- as.integer(ref)
  if (length(ref) != 1L || is.na(ref) || ref < 1L || ref > ngroups) {
    lav_msg_stop(gettextf(
      "reference_group must be a group number between 1 and %d, or a group
       label.", ngroups))
  }
  ref
}

# multimodality: the random starts of the alignment of the scales and
# means (starts$value: the loss at the end of every start); a warning when
# the best value was reached by a single start
lav_align_check_starts <- function(starts = NULL, what = "alignment") {
  if (is.null(starts) || length(starts$value) < 2L) {
    return(invisible(NULL))
  }
  v <- starts$value[is.finite(starts$value)]
  if (length(v) < 2L) {
    return(invisible(NULL))
  }
  best <- min(v)
  n_best <- sum(v <= best + 1e-4 * max(1, abs(best)))
  if (n_best < 2L) {
    lav_msg_warn(gettextf(
      "the best solution of the %1$s was found by only one of %2$d starts:
       the loss function has several local optima; consider increasing
       the number of random starts (rstarts) and compare the solutions
       (lavInspect(, \"rotation\")).", what, length(v)))
  }
  invisible(n_best)
}

# identification of a free alignment (the means of all groups are aligned):
# the overall level of the factor means is identified by the loss only,
# and poorly so when there is little non-invariance or few groups. Per
# factor: the standard error of the average of the group means (the level)
# against the largest standard error of a difference between two group
# means (identified by the data); a warning when the level is the less
# precise of the two
lav_align_check_free_means <- function(lavpartable = NULL, vcov = NULL,
                                       lv_names = NULL,
                                       what = "alignment = \"free\"") {
  out <- NULL
  for (lv in lv_names) {
    idx <- which(lavpartable$op == "~1" & lavpartable$lhs == lv &
                 lavpartable$free > 0L)
    if (length(idx) < 2L) {
      next
    }
    v <- vcov[lavpartable$free[idx], lavpartable$free[idx], drop = FALSE]
    if (any(!is.finite(v))) {
      next
    }
    se_level <- sqrt(max(mean(v), 0))
    pairs <- utils::combn(length(idx), 2L)
    se_diff <- sqrt(pmax(v[cbind(pairs[1L, ], pairs[1L, ])] +
                         v[cbind(pairs[2L, ], pairs[2L, ])] -
                         2 * v[cbind(pairs[1L, ], pairs[2L, ])], 0))
    out <- rbind(out, data.frame(lv = lv, se_level = se_level,
                                 se_diff = max(se_diff),
                                 stringsAsFactors = FALSE))
  }
  if (!is.null(out) && any(out$se_level > out$se_diff)) {
    bad <- out$lv[out$se_level > out$se_diff]
    lav_msg_warn(gettextf(
      "the free alignment of the factor means may be poorly identified
       (%1$s): the standard error of the overall level of the means of %2$s
       exceeds that of the differences between the groups; consider the
       fixed alignment (a reference group), more groups, or more
       non-invariance.", what, lav_msg_view(bad, log_sep = "none")))
  }
  invisible(out)
}
