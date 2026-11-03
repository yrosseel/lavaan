# invariance analysis of the aligned measurement parameters of a
# multiple-group ESEM with an agreement rotation (rotation.args =
# list(mg_agreement = TRUE)), following the post-estimation procedure of the
# alignment method (Asparouhov & Muthen, 2014, SEM journal):
#
# for every (aligned) loading and intercept of the efa blocks,
# 1. test the equality of the parameter for every pair of groups (z test,
#    using the covariance matrix of the aligned estimates), and 'connect'
#    two groups when the p-value exceeds p_connect (0.01);
# 2. the starting invariant set is the largest connected set of groups;
# 3. repeat: compute the average of the parameter over the invariant set,
#    test every group against that average, add the groups with a p-value
#    above p_set (0.001) and remove the others, until the set is stable;
# 4. report, per parameter, the invariant set, the non-invariant groups,
#    the average over the invariant set, the R-square (the proportion of the
#    across-group variance of the unaligned parameter that is explained by
#    the group differences in factor variances and means) and the
#    contribution of the parameter to the agreement (loss) function.
#
# YR 02 Nov 2026: first version

lavTestAlignment <- function(object, p_connect = 0.01, p_set = 0.001, # nolint
                             max_iter = 100L) {
  if (!inherits(object, "lavaan")) {
    lav_msg_stop(gettext("object is not a lavaan object."))
  }
  object <- lav_object_check_version(object)
  lavoptions <- object@Options
  ropts <- lavoptions$rotation.args
  ngroups <- object@Data@ngroups
  if (ngroups < 2L) {
    lav_msg_stop(gettext("lavTestAlignment() requires a multiple-group model."))
  }
  if (object@Data@nlevels > 1L) {
    lav_msg_stop(gettext("lavTestAlignment() is not available for multilevel
                          models."))
  }
  aligned_cfa <- !is.null(lavoptions$alignment) &&
    lavoptions$alignment != "none"
  aligned_efa <- object@Model@nefa > 0L && lavoptions$rotation != "none" &&
    isTRUE(ropts$mg_agreement) && !("loadings" %in% lavoptions$group.equal)
  if (!aligned_cfa && !aligned_efa) {
    lav_msg_stop(gettext("lavTestAlignment() requires efa blocks that were
                          rotated with rotation.args = list(mg_agreement =
                          TRUE), or a confirmatory model with alignment =
                          \"fixed\" or \"free\"."))
  }
  # the loss settings: those of the agreement rotation for the efa
  # parameters, those of the confirmatory alignment for the others
  aopts <- lavoptions$alignment.args
  copts <- NULL
  if (aligned_cfa) {
    copts <- list(
      mg_agreement_method = aopts$method,
      mg_agreement_epsilon = aopts$epsilon,
      mg_agreement_group_weights = aopts$group_weights,
      mg_agreement_intercepts = lavoptions$alignment
    )
  }
  if (!aligned_efa) {
    ropts <- copts
  }
  if (lavoptions$se == "none") {
    lav_msg_stop(gettext("lavTestAlignment() requires standard errors
                          (se is \"none\")."))
  }
  vcov <- lavInspect(object, "vcov")
  vcov <- unclass(vcov)
  attr(vcov, "header") <- NULL
  if (is.null(vcov) || any(!is.finite(diag(vcov)))) {
    lav_msg_stop(gettext("the covariance matrix of the parameter estimates
                          is not available."))
  }
  pt <- object@ParTable
  group_label <- object@Data@group.label
  if (length(group_label) == 0L) {
    group_label <- as.character(seq_len(ngroups))
  }
  lv_efa_only <- unique(pt$lhs[pt$op == "=~" & !is.na(pt$efa) &
                                nzchar(pt$efa)])
  lv_cfa <- character(0L)
  if (aligned_cfa) {
    lv_cfa <- lav_pt_vnames(pt, type = "lv.regular", block = 1L)
    lv_cfa <- lv_cfa[!lv_cfa %in% lv_efa_only]
  }
  if (!aligned_efa) {
    lv_efa_only <- character(0L)
  }
  lv_efa <- c(lv_efa_only, lv_cfa)
  ov_efa_only <- unique(pt$rhs[pt$op == "=~" & pt$lhs %in% lv_efa_only])
  ov_cfa <- unique(pt$rhs[pt$op == "=~" & pt$lhs %in% lv_cfa])
  ov_efa <- unique(c(ov_efa_only, ov_cfa))
  align_int_efa <- aligned_efa && !is.null(ropts$mg_agreement_intercepts) &&
    ropts$mg_agreement_intercepts != "none"

  # the parameters: the loadings of the efa blocks and of the aligned
  # confirmatory factors, and the intercepts of their indicators (those of
  # the efa indicators only if the efa intercepts were aligned)
  par_idx <- which(pt$op == "=~" & pt$lhs %in% lv_efa & pt$group == 1L)
  ov_int <- ov_cfa
  if (align_int_efa) {
    ov_int <- unique(c(ov_int, ov_efa_only))
  }
  if (length(ov_int) > 0L) {
    # the intercepts (free; those of ordered indicators are fixed to zero)
    # and the thresholds of the ordered indicators
    par_idx <- c(par_idx,
                 which(pt$op == "~1" & pt$lhs %in% ov_int & pt$group == 1L &
                       pt$free > 0L),
                 which(pt$op == "|" & pt$lhs %in% ov_int & pt$group == 1L))
  }
  npar <- length(par_idx)
  if (npar == 0L) {
    lav_msg_stop(gettext("no aligned parameters were found."))
  }

  # the loss function (for the contributions): per parameter, that of the
  # agreement rotation (efa parameters) or of the confirmatory alignment
  n_g <- unlist(object@SampleStats@nobs)
  loss_of <- function(opts) {
    w_g <- rep(1, ngroups)
    if (opts$mg_agreement_group_weights %in% c("nobs", "sqrt.nobs")) {
      w_g <- n_g / mean(n_g)
      if (opts$mg_agreement_group_weights == "sqrt.nobs") {
        w_g <- sqrt(w_g)
      }
    }
    f <- function(d) {
      if (opts$mg_agreement_method == "procrustes") {
        d * d
      } else {
        sqrt(sqrt(d * d + opts$mg_agreement_epsilon))
      }
    }
    list(w_g = w_g, f = f,
         relative = opts$mg_agreement_method == "relative.alignment",
         epsilon = opts$mg_agreement_epsilon)
  }
  loss_efa <- loss_of(ropts)
  loss_cfa <- if (!is.null(copts)) loss_of(copts) else loss_efa

  # per-group model matrices (for the R-square)
  est <- lavInspect(object, "est")

  pairs <- utils::combn(ngroups, 2L)
  p_mat <- matrix(NA_real_, npar, ngroups)
  status <- matrix(NA, npar, ngroups)
  average <- contribution <- r2 <- numeric(npar)
  invariant <- noninvariant <- character(npar)
  for (i in seq_len(npar)) {
    lhs <- pt$lhs[par_idx[i]]
    op <- pt$op[par_idx[i]]
    rhs <- pt$rhs[par_idx[i]]
    row_idx <- vapply(seq_len(ngroups), function(g) {
      idx <- which(pt$lhs == lhs & pt$op == op & pt$rhs == rhs &
                   pt$group == g)
      if (length(idx) != 1L) NA_integer_ else idx
    }, integer(1L))
    if (any(is.na(row_idx))) {
      next
    }
    theta <- pt$est[row_idx]
    free_idx <- pt$free[row_idx]
    if (any(free_idx == 0L)) {
      # a fixed parameter in some group: no test
      next
    }
    v_mat <- vcov[free_idx, free_idx, drop = FALSE]

    # 1. pairwise tests + connected sets
    adj <- matrix(FALSE, ngroups, ngroups)
    for (p in seq_len(ncol(pairs))) {
      g <- pairs[1L, p]
      h <- pairs[2L, p]
      v_d <- v_mat[g, g] + v_mat[h, h] - 2 * v_mat[g, h]
      p_val <- if (v_d > 0) {
        2 * stats::pnorm(-abs(theta[g] - theta[h]) / sqrt(v_d))
      } else {
        1
      }
      adj[g, h] <- adj[h, g] <- (p_val > p_connect)
    }
    # 2. the largest connected set (the first one in case of ties)
    visited <- rep(FALSE, ngroups)
    best_set <- integer(0L)
    for (g in seq_len(ngroups)) {
      if (visited[g]) {
        next
      }
      comp <- g
      queue <- g
      visited[g] <- TRUE
      while (length(queue) > 0L) {
        node <- queue[1L]
        queue <- queue[-1L]
        nbrs <- which(adj[node, ] & !visited)
        visited[nbrs] <- TRUE
        comp <- c(comp, nbrs)
        queue <- c(queue, nbrs)
      }
      if (length(comp) > length(best_set)) {
        best_set <- sort(comp)
      }
    }
    # 3. iterate: test every group against the average of the set
    set <- best_set
    seen <- list(set)
    p_g <- rep(NA_real_, ngroups)
    for (iter in seq_len(max_iter)) {
      for (g in seq_len(ngroups)) {
        cvec <- numeric(ngroups)
        cvec[g] <- 1
        cvec[set] <- cvec[set] - 1 / length(set)
        v_d <- as.numeric(crossprod(cvec, v_mat %*% cvec))
        d <- sum(cvec * theta)
        p_g[g] <- if (v_d > 0) 2 * stats::pnorm(-abs(d) / sqrt(v_d)) else 1
      }
      new_set <- which(p_g > p_set)
      if (length(new_set) == 0L) {
        # nothing passes: keep the current set
        new_set <- set
      }
      if (identical(new_set, set)) {
        break
      }
      # a cycle: stop with the current set
      if (any(vapply(seen, identical, logical(1L), new_set))) {
        break
      }
      seen <- c(seen, list(new_set))
      set <- new_set
    }
    p_mat[i, ] <- p_g
    status[i, ] <- seq_len(ngroups) %in% set
    average[i] <- mean(theta[set])
    invariant[i] <- paste(group_label[set], collapse = " ")
    noninvariant[i] <- paste(group_label[-set], collapse = " ")

    # 4. contribution to the loss, and R-square
    is_efa <- (op == "=~" && lhs %in% lv_efa_only) ||
      (op %in% c("~1", "|") && lhs %in% ov_efa_only && !(lhs %in% ov_cfa))
    lo <- if (is_efa) loss_efa else loss_cfa
    theta_c <- theta
    scale_c <- rep(1, ngroups)
    if (lo$relative) {
      # the relative loss: the loadings enter as log(lambda^2 + eps)/2,
      # an intercept difference is divided by the average over the pair of
      # the (damped) norm of the row of loadings of the indicator (over
      # the factors of its efa block, or over all factors for a
      # confirmatory indicator)
      if (op == "=~") {
        theta_c <- 0.5 * log(theta * theta + lo$epsilon)
      } else {
        lv_row <- if (is_efa) {
          unique(pt$lhs[pt$op == "=~" & pt$rhs == lhs &
                        pt$lhs %in% lv_efa_only])
        } else {
          colnames(est[[1L]]$lambda)
        }
        scale_c <- vapply(seq_len(ngroups), function(g) {
          sqrt(sum(est[[g]]$lambda[lhs, lv_row]^2) + lo$epsilon)
        }, numeric(1L))
      }
    }
    contribution[i] <- sum(vapply(seq_len(ncol(pairs)), function(p) {
      g <- pairs[1L, p]
      h <- pairs[2L, p]
      d <- (theta_c[g] - theta_c[h]) / ((scale_c[g] + scale_c[h]) / 2)
      lo$w_g[g] * lo$w_g[h] * lo$f(d)
    }, numeric(1L)))
    theta_bar <- mean(theta)
    if (op == "=~") {
      # unaligned loading: the aligned one times the factor's sd; the
      # residual after absorbing the scale: sd_g (theta_g - mean)
      sd_g <- vapply(seq_len(ngroups), function(g) {
        sqrt(est[[g]]$psi[lhs, lhs])
      }, numeric(1L))
      theta_0 <- theta * sd_g
      resid <- sd_g * (theta - theta_bar)
    } else {
      # unaligned intercept: nu + lambda alpha (threshold: tau - lambda
      # alpha); residual: the aligned value minus the mean
      sign_op <- if (op == "|") -1 else 1
      theta_0 <- vapply(seq_len(ngroups), function(g) {
        theta[g] + sign_op *
          sum(est[[g]]$lambda[lhs, lv_efa] * est[[g]]$alpha[lv_efa, 1L])
      }, numeric(1L))
      resid <- theta - theta_bar
    }
    # (the across-group variance needs at least three groups)
    v_0 <- if (ngroups > 2L) stats::var(theta_0) else NA_real_
    r2[i] <- if (is.finite(v_0) && v_0 > 0) 1 - stats::var(resid) / v_0 else NA_real_
  }

  out <- data.frame(
    lhs = pt$lhs[par_idx], op = pt$op[par_idx], rhs = pt$rhs[par_idx],
    average = average, invariant = invariant, noninvariant = noninvariant,
    r2 = r2, contribution = contribution,
    stringsAsFactors = FALSE
  )
  rownames(p_mat) <- rownames(status) <- trimws(paste(out$lhs, out$op, out$rhs))
  colnames(p_mat) <- colnames(status) <- group_label
  attr(out, "pvalue") <- p_mat
  attr(out, "invariant") <- status
  attr(out, "p_connect") <- p_connect
  attr(out, "p_set") <- p_set
  attr(out, "header") <- sprintf(
    "Approximate measurement invariance of the aligned parameters (pairwise p > %s connects; p > %s keeps a group in the invariant set)",
    format(p_connect), format(p_set))
  class(out) <- c("lavaan.data.frame", "data.frame")
  out
}
