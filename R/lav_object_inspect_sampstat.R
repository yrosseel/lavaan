# lavInspect()/lavTech() accessors for the data, sample statistics and
# h1 (unrestricted model) internals that are otherwise only reachable
# through the lavData, lavSampleStats and h1 slots
#
# - sampling weights (per group)
# - missing data patterns: frequencies and per-pattern statistics
# - response patterns for ordered data
# - sample variances, inverse covariance matrix and its log-determinant
#   (both the unconditional and the residual (conditional_x) versions)
# - thresholds ignoring the covariates (categorical + conditional_x)
# - the diagonal of the WLS weight matrix
# - (log) group weights, per-group loglikelihood of the h0 and h1 models
# - two-level intermediate statistics (pooled within/between, per
#   cluster-size means/covariances, ...)

# helper: turn a per-group list into the user-facing return value
lav_inspect_finalize_list <- function(return_value, object,
                                      drop_list_single_group = FALSE) {
  n_g <- length(return_value)
  if (n_g == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else {
    if (length(object@Data@group.label) > 0L &&
        length(object@Data@group.label) == n_g) {
      names(return_value) <- unlist(object@Data@group.label)
    }
  }
  return_value
}

# helper: same, but for per-block lists
lav_inspect_finalize_block_list <- function(return_value, object,
                                            drop_list_single_group = FALSE) {
  nblocks <- length(return_value)
  if (nblocks == 1L && drop_list_single_group) {
    return_value <- return_value[[1]]
  } else if (nblocks > 1L) {
    names(return_value) <- object@Data@block.label
  }
  return_value
}

# helper: per-group scalar -> (named) numeric vector
lav_inspect_finalize_scalar <- function(return_value, object,
                                        add_labels = FALSE,
                                        add_class = FALSE) {
  return_value <- unlist(return_value)
  if (add_labels && length(object@Data@group.label) > 0L &&
      length(object@Data@group.label) == length(return_value)) {
    names(return_value) <- unlist(object@Data@group.label)
  }
  if (add_class) {
    class(return_value) <- c("lavaan.vector", "numeric")
  }
  return_value
}


#### data ####

# (normalized) sampling weights, per group
lav_inspect_sampling_weights <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavdata <- object@Data
  if (length(lavdata@sampling.weights) == 0L) {
    lav_msg_stop(gettext(
      "no sampling weights were used to fit this model."))
  }

  n_g <- lavdata@ngroups
  return_value <- vector("list", n_g)
  for (g in seq_len(n_g)) {
    return_value[[g]] <- as.numeric(lavdata@weights[[g]])
    if (add_labels) {
      names(return_value[[g]]) <- lavdata@case.idx[[g]]
    }
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.vector", "numeric")
    }
  }

  lav_inspect_finalize_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}

# missing data pattern frequencies, per group; the order of the
# frequencies matches the rows of lavInspect(object, "patterns")
lav_inspect_mi_patterns_freq <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavdata <- object@Data
  n_g <- lavdata@ngroups
  return_value <- vector("list", n_g)

  for (g in seq_len(n_g)) {
    if (!is.null(lavdata@Mp[[g]])) {
      return_value[[g]] <- as.integer(lavdata@Mp[[g]]$freq)
    } else {
      # complete data: a single pattern
      return_value[[g]] <- as.integer(lavdata@nobs[[g]])
    }
    if (add_class) {
      class(return_value[[g]]) <- c("lavaan.vector", "numeric")
    }
  }

  lav_inspect_finalize_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}

# per-pattern sample statistics (missing = "ml", "ml.x", "two.stage" or
# "robust.two.stage"): for each missing data pattern, the mean vector and
# covariance matrix of the observed variables in that pattern, the
# (weighted) number of cases, and the pattern itself
lav_inspect_mi_patterns_stats <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavdata <- object@Data
  lavsamplestats <- object@SampleStats
  n_g <- lavdata@ngroups

  if (length(lavsamplestats@missing) == 0L ||
      all(sapply(lavsamplestats@missing, is.null))) {
    lav_msg_stop(gettext(
      "per-pattern sample statistics are only available if the model was
       fitted with missing = \"ml\", \"ml.x\", \"two.stage\" or
       \"robust.two.stage\"."))
  }

  return_value <- vector("list", n_g)
  for (g in seq_len(n_g)) {
    yp <- lavsamplestats@missing[[g]]
    if (is.null(yp)) {
      next
    }
    ov_names <- object@pta$vnames$ov.model[[g]]
    if (lavdata@nlevels > 1L) {
      # the between-only variables are not part of the level-1 patterns
      between_idx <- lavdata@Lp[[g]]$between.idx[[2]]
      if (length(between_idx) > 0L) {
        ov_names <- ov_names[-between_idx]
      }
    }
    npatterns <- length(yp)
    out <- vector("list", npatterns)
    for (p in seq_len(npatterns)) {
      var_idx <- as.logical(yp[[p]]$var.idx)
      nvar_p <- sum(var_idx)
      cov_p <- yp[[p]]$SY
      # a single case (freq == 1): the covariance matrix is stored as 0
      if (length(cov_p) == 1L && nvar_p > 1L) {
        cov_p <- matrix(0, nvar_p, nvar_p)
      } else {
        cov_p <- matrix(cov_p, nvar_p, nvar_p)
      }
      mean_p <- as.numeric(yp[[p]]$MY)
      if (add_labels) {
        names(var_idx) <- ov_names
        rownames(cov_p) <- colnames(cov_p) <- ov_names[var_idx]
        names(mean_p) <- ov_names[var_idx]
      }
      if (add_class) {
        class(cov_p) <- c("lavaan.matrix.symmetric", "matrix")
        class(mean_p) <- c("lavaan.vector", "numeric")
      }
      out[[p]] <- list(cov = cov_p, mean = mean_p,
                       freq = as.numeric(yp[[p]]$freq), var_idx = var_idx)
    }
    return_value[[g]] <- out
  }

  lav_inspect_finalize_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}

# response patterns for the ordered variables, per group
lav_inspect_resp_patterns <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavdata <- object@Data
  if (length(lavdata@ordered) == 0L) {
    lav_msg_stop(gettext(
      "response patterns are only available for ordered variables."))
  }
  if (lavdata@data.type != "full") {
    lav_msg_stop(gettext(
      "response patterns are only available if raw data was used."))
  }

  n_g <- lavdata@ngroups
  return_value <- vector("list", n_g)
  for (g in seq_len(n_g)) {
    rp <- lavdata@Rp[[g]]
    if (is.null(rp)) {
      next
    }
    pat <- rp$pat
    freq <- as.integer(rownames(pat))
    rownames(pat) <- NULL
    if (add_labels) {
      ov_names <- lavdata@ov.names[[g]]
      ord_idx <- which(ov_names %in% lavdata@ov$name[lavdata@ov$type ==
                                                     "ordered"])
      colnames(pat) <- ov_names[ord_idx]
    } else {
      colnames(pat) <- NULL
    }
    if (add_class) {
      class(pat) <- c("lavaan.matrix", "matrix")
      class(freq) <- c("lavaan.vector", "numeric")
    }
    return_value[[g]] <- list(pat = pat, freq = freq,
      npatterns = as.integer(rp$npatterns),
      total_patterns = as.numeric(rp$total.patterns),
      empty_patterns = as.numeric(rp$empty.patterns))
  }

  lav_inspect_finalize_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}


#### sample statistics ####

# sample variances: the diagonal of the (h1-based, if available) sample
# covariance matrix, as returned by lavInspect(object, "sampstat"); with
# res = TRUE, the residual variances (conditional_x = TRUE)
lav_inspect_sampstat_var <- function(object, res = FALSE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavsamplestats <- object@SampleStats
  lavmodel <- object@Model
  conditional_x <- lavmodel@conditional.x && object@Data@nlevels == 1L

  if (res && !conditional_x) {
    lav_msg_stop(gettext(
      "residual variances are only available if conditional_x = TRUE."))
  }

  sampstat <- lav_inspect_sampstat(object, h1 = TRUE,
    add_labels = add_labels, add_class = FALSE,
    drop_list_single_group = FALSE)
  nblocks <- length(sampstat)

  return_value <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    if (res) {
      return_value[[b]] <- diag(sampstat[[b]]$res.cov)
    } else if (!conditional_x) {
      return_value[[b]] <- diag(sampstat[[b]]$cov)
    } else {
      # conditional_x: the joint (y, x) sample variances
      if (is.null(lavsamplestats@var[[b]])) {
        lav_msg_stop(gettext(
          "the unconditional sample variances are not available for this
           model; use \"res.var\" to obtain the residual variances."))
      }
      return_value[[b]] <- as.numeric(lavsamplestats@var[[b]])
      if (add_labels) {
        names(return_value[[b]]) <- object@pta$vnames$ov[[b]]
      }
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  lav_inspect_finalize_block_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}

# inverse of the sample covariance matrix (icov = TRUE) or the
# log-determinant of the sample covariance matrix (icov = FALSE), per
# group; with res = TRUE, the residual (conditional_x = TRUE) versions;
# if the slot was not filled in, the values are computed from the
# (h1-based) sample covariance matrix
lav_inspect_sampstat_icov <- function(object, res = FALSE, icov = TRUE,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavsamplestats <- object@SampleStats
  lavmodel <- object@Model
  conditional_x <- lavmodel@conditional.x && object@Data@nlevels == 1L

  if (res && !conditional_x) {
    lav_msg_stop(gettext(
      "residual (inverse) covariance matrices are only available if
       conditional_x = TRUE."))
  }
  if (object@Data@nlevels > 1L) {
    lav_msg_stop(gettext(
      "the inverse sample covariance matrix is not available for
       multilevel models."))
  }

  if (res) {
    slot_icov <- lavsamplestats@res.icov
    slot_logdet <- lavsamplestats@res.cov.log.det
  } else {
    slot_icov <- lavsamplestats@icov
    slot_logdet <- lavsamplestats@cov.log.det
  }

  sampstat <- lav_inspect_sampstat(object, h1 = TRUE,
    add_labels = add_labels, add_class = FALSE,
    drop_list_single_group = FALSE)
  nblocks <- length(sampstat)

  return_value <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    if (res) {
      cov_b <- sampstat[[b]]$res.cov
    } else if (!conditional_x) {
      cov_b <- sampstat[[b]]$cov
    } else {
      cov_b <- lavsamplestats@cov[[b]]
      if (is.null(cov_b)) {
        lav_msg_stop(gettext(
          "the unconditional sample covariance matrix is not available for
           this model; use \"res.icov\" or \"res.cov.log.det\" instead."))
      }
      if (add_labels) {
        rownames(cov_b) <- colnames(cov_b) <- object@pta$vnames$ov[[b]]
      }
    }
    if (icov) {
      if (!is.null(slot_icov[[b]])) {
        out <- slot_icov[[b]]
      } else {
        out <- lav_mat_sym_inverse(cov_b, logdet = FALSE)
      }
      if (add_labels) {
        rownames(out) <- colnames(out) <- rownames(cov_b)
      }
      if (add_class) {
        class(out) <- c("lavaan.matrix.symmetric", "matrix")
      }
    } else {
      if (!is.null(slot_logdet[[b]])) {
        out <- as.numeric(slot_logdet[[b]])
      } else {
        out <- as.numeric(determinant(cov_b, logarithm = TRUE)$modulus)
      }
    }
    return_value[[b]] <- out
  }

  if (icov) {
    return_value <- lav_inspect_finalize_block_list(return_value, object,
      drop_list_single_group = drop_list_single_group)
  } else {
    return_value <- lav_inspect_finalize_scalar(return_value, object,
      add_labels = add_labels, add_class = add_class)
  }

  return_value
}

# thresholds ignoring the exogenous covariates (categorical data);
# only differs from "th" if conditional_x = TRUE
lav_inspect_sampstat_th_nox <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavsamplestats <- object@SampleStats
  lavmodel <- object@Model

  if (!lavmodel@categorical) {
    lav_msg_stop(gettext(
      "thresholds ignoring the covariates are only available for
       categorical data."))
  }

  nblocks <- lavmodel@nblocks
  return_value <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    if (!is.null(lavsamplestats@res.th.nox[[b]])) {
      th <- as.numeric(lavsamplestats@res.th.nox[[b]])
    } else if (!is.null(lavsamplestats@th[[b]])) {
      # no covariates (conditional_x = FALSE): the plain thresholds
      th <- as.numeric(lavsamplestats@th[[b]])
    } else {
      lav_msg_stop(gettext(
        "thresholds ignoring the covariates are not available for this
         object."))
    }
    if (length(lavmodel@num.idx[[b]]) > 0L) {
      num_idx <- which(lavmodel@th.idx[[b]] == 0)
      th <- th[-num_idx]
    }
    if (add_labels) {
      names(th) <- object@pta$vnames$th[[b]]
    }
    if (add_class) {
      class(th) <- c("lavaan.vector", "numeric")
    }
    return_value[[b]] <- th
  }

  lav_inspect_finalize_block_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}

# the diagonal of the WLS weight matrix (per block); this is what DWLS
# and ULS use (the full matrix "wls.v" is diagonal in that case)
lav_inspect_wls_vd <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  wls_v <- lav_inspect_wls_v(object,
    add_labels = add_labels, add_class = FALSE,
    drop_list_single_group = FALSE)
  nblocks <- length(wls_v)

  return_value <- vector("list", nblocks)
  for (b in seq_len(nblocks)) {
    return_value[[b]] <- diag(wls_v[[b]])
    if (add_labels && !is.null(rownames(wls_v[[b]]))) {
      names(return_value[[b]]) <- rownames(wls_v[[b]])
    }
    if (add_class) {
      class(return_value[[b]]) <- c("lavaan.vector", "numeric")
    }
  }

  lav_inspect_finalize_block_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}


#### group weights and loglikelihoods ####

# observed group weights: the log of the (effective) number of
# observations in each group, in the same metric as the group.w elements
# of the "implied" and "sampstat" output
lav_inspect_group_w <- function(object,
    add_labels = FALSE, add_class = FALSE) {

  lavsamplestats <- object@SampleStats
  n_g <- lavsamplestats@ngroups
  return_value <- vector("list", n_g)
  for (g in seq_len(n_g)) {
    return_value[[g]] <-
      log(lavsamplestats@group.w[[g]] * lavsamplestats@ntotal)
  }

  lav_inspect_finalize_scalar(return_value, object,
    add_labels = add_labels, add_class = add_class)
}

# per-group loglikelihood of the fitted (h0) model or of the
# unrestricted (h1) model
lav_inspect_logl_group <- function(object, h1 = FALSE,
    add_labels = FALSE, add_class = FALSE) {

  if (h1) {
    if (length(object@h1) == 0L || is.null(object@h1$logl)) {
      lav_msg_stop(gettext("h1 slot is not available; refit with h1 = TRUE"))
    }
    logl <- object@h1$logl
  } else {
    if (length(object@loglik) == 0L) {
      lav_msg_stop(gettext(
        "no loglikelihood information is available for this model."))
    }
    logl <- object@loglik
  }
  return_value <- as.numeric(logl$loglik.group)

  # older objects (or stored h1 slots) may only have the total
  if (length(return_value) == 0L) {
    if (object@Data@ngroups == 1L) {
      return_value <- as.numeric(logl$loglik)
    } else {
      lav_msg_stop(gettext(
        "per-group loglikelihood values are not available for this
         object."))
    }
  }

  lav_inspect_finalize_scalar(as.list(return_value), object,
    add_labels = add_labels, add_class = add_class)
}


#### two-level intermediate statistics ####

# the sufficient statistics of the two-level (level-2) data, per group:
# - y1y1: crossproduct of the level-1 data
# - y2: the cluster means (one row per cluster)
# - sigma_w, mu_w: pooled within-cluster covariance matrix (divided by
#   N - nclusters) and within-level means
# - sigma_b, mu_b: the (moment-based) between-level covariance matrix
#   and means (used as starting values for the h1 model)
# - s_b: the between covariance matrix of the cluster means
# - s: the average cluster size (used in the s_b/sigma_b scaling)
# - mean_d, cov_d: means and covariance matrices of the cluster means,
#   per cluster size (in the order of "cluster.sizes")
# - loglik_x: the loglikelihood contribution of the fixed.x covariates
lav_inspect_sampstat_2l <- function(object,
    add_labels = FALSE, add_class = FALSE, drop_list_single_group = FALSE) {

  lavdata <- object@Data
  lavsamplestats <- object@SampleStats

  if (lavdata@nlevels == 1L) {
    lav_msg_stop(gettext(
      "two-level sample statistics are only available for multilevel
       models."))
  }
  if (length(lavsamplestats@YLp) == 0L ||
      all(sapply(lavsamplestats@YLp, is.null))) {
    lav_msg_stop(gettext(
      "two-level sample statistics are not available for this object."))
  }

  n_g <- lavdata@ngroups
  return_value <- vector("list", n_g)
  for (g in seq_len(n_g)) {
    ylp <- lavsamplestats@YLp[[g]][[2]]
    lp <- lavdata@Lp[[g]]
    if (is.null(ylp)) {
      next
    }

    # the level-1 data columns
    ov_names <- lavdata@ov.names[[g]]
    # the mean_d/cov_d columns are reordered: between-only variables
    # first, then the others (in their original order)
    between_idx <- lp$between.idx[[2]]
    both_idx <- lp$both.idx[[2]]
    within_idx <- lp$within.idx[[2]]
    d_idx <- c(between_idx, sort.int(c(both_idx, within_idx)))
    d_names <- ov_names[d_idx]

    sigma_w <- ylp$Sigma.W
    mu_w <- as.numeric(ylp$Mu.W)
    sigma_b <- ylp$Sigma.B
    mu_b <- as.numeric(ylp$Mu.B)
    s_b <- ylp$S.b
    y1y1 <- ylp$Y1Y1
    y2 <- ylp$Y2
    dimnames(y2) <- NULL
    mean_d <- lapply(ylp$mean.d, as.numeric)
    cov_d <- lapply(ylp$cov.d, function(x) {
      if (is.matrix(x)) {
        dimnames(x) <- NULL
      }
      x
    })

    if (add_labels) {
      rownames(sigma_w) <- colnames(sigma_w) <- ov_names
      rownames(sigma_b) <- colnames(sigma_b) <- ov_names
      rownames(s_b) <- colnames(s_b) <- ov_names
      rownames(y1y1) <- colnames(y1y1) <- ov_names
      colnames(y2) <- ov_names
      rownames(y2) <- lp$cluster.id[[2]]
      names(mu_w) <- names(mu_b) <- ov_names
      for (clz in seq_along(mean_d)) {
        names(mean_d[[clz]]) <- d_names
        if (is.matrix(cov_d[[clz]])) {
          rownames(cov_d[[clz]]) <- colnames(cov_d[[clz]]) <- d_names
        }
      }
      names(mean_d) <- names(cov_d) <- lp$cluster.sizes[[2]]
    }
    if (add_class) {
      class(sigma_w) <- class(sigma_b) <- class(s_b) <- class(y1y1) <-
        c("lavaan.matrix.symmetric", "matrix")
      class(y2) <- c("lavaan.matrix", "matrix")
      class(mu_w) <- class(mu_b) <- c("lavaan.vector", "numeric")
      for (clz in seq_along(mean_d)) {
        class(mean_d[[clz]]) <- c("lavaan.vector", "numeric")
        if (is.matrix(cov_d[[clz]])) {
          class(cov_d[[clz]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
      }
    }

    return_value[[g]] <- list(
      sigma_w = sigma_w, mu_w = mu_w,
      sigma_b = sigma_b, mu_b = mu_b,
      s_b = s_b, s = as.numeric(ylp$s),
      y1y1 = y1y1, y2 = y2,
      mean_d = mean_d, cov_d = cov_d,
      loglik_x = as.numeric(ylp$loglik.x))
  }

  lav_inspect_finalize_list(return_value, object,
    drop_list_single_group = drop_list_single_group)
}
