# adapt marker (bad.marker.crit)
#
# By default, lavaan sets the metric of a latent variable by fixing the
# factor loading of the FIRST indicator to 1.0. This works well, except when
# the first indicator happens to be a poor item (i.e., it correlates weakly
# with the other indicators of the same factor). In that case, fixing its
# loading to 1.0 often leads to (serious) convergence problems, and the user
# has no idea why.
#
# When bad.marker.crit > 0 (the default is 0.1, and auto.fix.first = TRUE),
# lavaan inspects -- after the unrestricted (h1) sample statistics are
# available -- whether the first indicator of each latent variable is a poor
# item. If so, a warning is issued and another (better) indicator is used as
# the marker instead. If the first indicator is good enough, nothing happens.
#
# The 'quality' of an indicator is measured by its corrected item-total
# correlation: the correlation between the indicator and the sum score of the
# OTHER indicators of the same factor. This is computed from the (pooled)
# unrestricted (h1) covariance matrix, so the behavior is consistent whether
# or not we have missing data, categorical data, etc. Before the item-total
# correlations are computed, the items are aligned with their dominant common
# direction (the sign pattern of the first eigenvector of their correlation
# matrix): without this, a single reverse-coded item degrades the sum score
# of the 'other' items, and a perfectly good marker item could be flagged as
# poor (github issue #628).
#
# We only switch in the 'clean default' situation, where the user did not
# set up (part of) the scaling or the measurement model deliberately. For
# each latent variable/composite, we require -- in EVERY block -- that
# exactly one indicator has a fixed loading/weight, that it is the FIRST
# indicator, and that its value is 1.0 (i.e., the auto.fix.first default).
# In addition, none of the loading/weight rows may carry user metadata: a
# label (which may tie the parameter to other parameters through equality
# constraints or ':=' definitions), a start value, a prior, or bounds; nor
# may any of them be referenced -- through their plabels -- in user-written
# constraint or definition rows (e.g. via the constraints= argument). See
# again github issue #628: with a labeled marker, switching freed the marker
# loading, but an equality constraint (same label on the -- still fixed --
# markers of other factors) kept it pinned to 1.0, while the new marker was
# also fixed to 1.0, yielding a doubly-scaled factor.
#
# Composites (op == "<~") get a similar treatment, but with a different
# quality measure: the metric of a composite is set by fixing the WEIGHT of
# the first indicator to 1.0. If the ML estimate of that weight is (near)
# zero, the remaining (free) weights must diverge to infinity to reach the
# optimum, and the model never converges. The implied (relative) weights can
# be gauged before estimation: under the composite model, the covariances of
# the indicators of one composite with ANY outside variable are proportional
# to a single vector l = S_jj w (S_jj = the within-block covariance matrix).
# We estimate l as the dominant (rank-1) left singular vector of the
# cross-correlation block between the indicators and all other observed
# variables, solve w = S_jj^{-1} l, and rescale so max |w| = 1. If the first
# indicator's relative weight is below the threshold (and a clearly better
# indicator exists), we switch the marker to the indicator with the largest
# absolute weight.
#
# lav_pt_marker_adapt() returns NULL if no change is needed, or a list
# with:
#   - marker: a named character vector (lv -> new indicator) for the latent
#             variables whose marker should be switched
#   - info:   a data.frame (lv, type, old, new, r.old, r.new) for the
#             warning; type is "factor" or "composite"; for factors r.* are
#             corrected item-total correlations, for composites relative
#             weights
#
# 'threshold' is the value below which the first indicator is considered a
# 'poor' marker (in absolute value of its corrected item-total correlation,
# or of its relative weight for composites); it is set by the
# 'bad.marker.crit' option (default 0.1)
lav_pt_marker_adapt <- function(lavpartable = NULL,
                                lavh1 = NULL,
                                lavoptions = NULL,
                                lavdata = NULL,
                                threshold = 0.1) {

  # we need the unrestricted (h1) covariance matrix/matrices
  implied <- lavh1$implied
  if (is.null(implied)) {
    return(NULL)
  }
  cov_list <- implied[["cov"]]
  if (is.null(cov_list)) {
    # conditional.x: use the residual covariance matrix
    cov_list <- implied[["res.cov"]]
  }
  if (is.null(cov_list) || length(cov_list) == 0L) {
    return(NULL)
  }

  # variable names per block (from the dimnames, with a single-level fallback)
  block_names <- lapply(seq_along(cov_list), function(b) {
    vn <- rownames(cov_list[[b]])
    if (is.null(vn) && lavdata@nlevels == 1L &&
        b <= length(lavdata@ov.names)) {
      vn <- lavdata@ov.names[[b]]
    }
    vn
  })

  # observed indicators (ordered indicators are fine: for those, the h1 'cov'
  # already holds the latent (polychoric) correlations)
  ov_names <- lav_pt_vnames(lavpartable, type = "ov")
  # exclude efa factors (composites are handled separately below:
  # lv.regular only contains =~ definitions)
  lv_regular <- lav_pt_vnames(lavpartable, type = "lv.regular")
  if (!is.null(lavpartable$efa)) {
    lv_efa <- unique(lavpartable$lhs[lavpartable$op == "=~" &
                                     nchar(lavpartable$efa) > 0L])
    lv_regular <- lv_regular[!lv_regular %in% lv_efa]
  }
  # composites
  lv_comp <- unique(lavpartable$lhs[lavpartable$op == "<~"])
  if (length(lv_regular) == 0L && length(lv_comp) == 0L) {
    return(NULL)
  }

  # corrected item-total correlations, given a covariance submatrix C
  #
  # the items are first aligned with their dominant common direction (the
  # sign pattern of the first eigenvector of their correlation matrix), so
  # that a reverse-coded item does not degrade the sum score of the 'other'
  # items; the returned values keep each item's own direction (relative to
  # the majority-positive orientation), so a reversed item still shows a
  # negative value
  cit_from_cov <- function(C) {
    s <- rep(1, nrow(C))
    d <- diag(C)
    if (all(is.finite(C)) && all(d > 0)) {
      isd <- 1 / sqrt(d)
      ev <- try(eigen(C * tcrossprod(isd), symmetric = TRUE),
                silent = TRUE)
      if (!inherits(ev, "try-error")) {
        s <- sign(ev$vectors[, 1L])
        s[s == 0] <- 1
        if (sum(s) < 0) {
          s <- -s
        }
        C <- C * tcrossprod(s)
      }
    }
    R <- rowSums(C)
    S <- sum(C)
    cov_rest <- R - d
    var_rest <- S - 2 * R + d
    out <- s * cov_rest / sqrt(d * var_rest)
    out[!is.finite(out)] <- NA_real_
    # a value outside [-1, 1] is not a correlation; this can only happen if
    # the (pseudo-)correlation matrix is not positive definite, and the
    # quality measure is meaningless for such an item
    out[abs(out) > 1 + sqrt(.Machine$double.eps)] <- NA_real_
    out
  }

  # TRUE if the user attached any metadata to these (loading/weight) rows: a
  # label (which may tie the parameter to other parameters through equality
  # constraints or ':=' definitions), a start value, a prior, or bounds; in
  # that case the user has deliberately set up (part of) the measurement
  # model, and switching the marker could silently change the meaning of the
  # model (github issue #628)
  user_modified_rows <- function(rows) {
    if (!is.null(lavpartable$label)) {
      lab <- lavpartable$label[rows]
      lab <- lab[nchar(lab) > 0L]
      # group.equal= writes the (auto-generated) plabels of the first group
      # into the label column; those are regenerated consistently when the
      # parameter table is rebuilt with a new marker, so only labels the
      # user provided should block the switch
      if (!is.null(lavpartable$plabel)) {
        lab <- lab[!lab %in% lavpartable$plabel]
      }
      if (length(lab) > 0L) {
        return(TRUE)
      }
    }
    # user-provided starting values on free parameters
    if (any(lavpartable$free[rows] != 0L &
            !is.na(lavpartable$ustart[rows]))) {
      return(TRUE)
    }
    if (!is.null(lavpartable$prior) &&
        any(nchar(lavpartable$prior[rows]) > 0L)) {
      return(TRUE)
    }
    if (!is.null(lavpartable$lower) &&
        any(is.finite(lavpartable$lower[rows]))) {
      return(TRUE)
    }
    if (!is.null(lavpartable$upper) &&
        any(is.finite(lavpartable$upper[rows]))) {
      return(TRUE)
    }
    FALSE
  }

  # user-written constraint/definition rows may also reference parameters
  # through their plabels (e.g. via the constraints= argument); rows
  # generated by lavaan itself (group.equal, effect.coding, ...) have
  # user != 1 and are regenerated consistently when the parameter table is
  # rebuilt with a new marker
  con_idx <- which(lavpartable$op %in% c("==", "<", ">", ":=") &
                   lavpartable$user == 1L)
  con_text <- c(lavpartable$lhs[con_idx], lavpartable$rhs[con_idx])
  plabel_referenced <- function(rows) {
    if (length(con_idx) == 0L || is.null(lavpartable$plabel)) {
      return(FALSE)
    }
    plabs <- lavpartable$plabel[rows]
    plabs <- plabs[nchar(plabs) > 0L]
    any(vapply(plabs, function(p) {
      any(grepl(p, con_text, fixed = TRUE))
    }, logical(1L)))
  }

  # TRUE if, in EVERY block, the rows for this lv show the clean default
  # scaling: the same indicators (ind) in the same order, with exactly one
  # fixed loading/weight (free == 0) that belongs to the FIRST indicator and
  # equals 1.0
  clean_default_scaling <- function(all_rows, ind) {
    for (b in unique(lavpartable$block[all_rows])) {
      rows_b <- all_rows[lavpartable$block[all_rows] == b &
                         lavpartable$rhs[all_rows] %in% ind]
      if (!identical(lavpartable$rhs[rows_b], ind)) {
        return(FALSE)
      }
      fixed <- which(lavpartable$free[rows_b] == 0L)
      if (length(fixed) != 1L || fixed != 1L ||
          !isTRUE(lavpartable$ustart[rows_b][1L] == 1)) {
        return(FALSE)
      }
    }
    TRUE
  }

  new_marker <- character(0L)
  info <- list()

  for (lv in lv_regular) {
    all_rows <- which(lavpartable$op == "=~" & lavpartable$lhs == lv)
    # leave the factor alone if the user attached any metadata to its
    # loadings (in any block), or references them in constraints
    if (user_modified_rows(all_rows) || plabel_referenced(all_rows)) {
      next
    }
    # only consider one block for the indicator list (the structure is
    # identical across blocks)
    b1 <- lavpartable$block[all_rows]
    lv_rows <- all_rows[b1 == b1[1L]]
    ind <- lavpartable$rhs[lv_rows]

    # only observed indicators, and at least two of them
    keep <- ind %in% ov_names
    lv_rows <- lv_rows[keep]
    ind <- ind[keep]
    if (length(ind) < 2L) {
      next
    }

    # only adapt the 'clean' default situation (in every block): exactly one
    # indicator has a fixed loading (free == 0), it is the FIRST indicator,
    # and its value is 1.0; otherwise the user set up the scaling and we
    # leave it alone
    if (!clean_default_scaling(all_rows, ind)) {
      next
    }
    cur <- ind[1L]

    # collect corrected item-total correlations across all blocks that
    # contain all the indicators of this factor
    cit_blocks <- list()
    for (b in seq_along(cov_list)) {
      vn <- block_names[[b]]
      if (is.null(vn) || !all(ind %in% vn)) {
        next
      }
      idx <- match(ind, vn) # the h1 cov matrices may have no dimnames
      C <- cov_list[[b]][idx, idx, drop = FALSE]
      cit_blocks[[length(cit_blocks) + 1L]] <- cit_from_cov(C)
    }
    if (length(cit_blocks) == 0L) {
      next
    }
    cit <- colMeans(do.call(rbind, cit_blocks), na.rm = TRUE)
    names(cit) <- ind
    if (all(is.na(cit))) {
      next
    }

    r <- abs(cit)
    r_cur <- r[[1L]]
    best <- which.max(r)
    # switch only if the first indicator is poor AND a clearly better
    # indicator exists
    if (!is.na(r_cur) && r_cur < threshold &&
        ind[best] != cur && r[[best]] >= threshold) {
      new_marker[lv] <- ind[best]
      info[[length(info) + 1L]] <- data.frame(
        lv = lv, type = "factor", old = cur, new = ind[best],
        r.old = round(cit[[1L]], 3), r.new = round(cit[[best]], 3),
        stringsAsFactors = FALSE
      )
    }
  }

  # composites: check the implied relative weight of the first indicator
  for (lv in lv_comp) {
    all_rows <- which(lavpartable$op == "<~" & lavpartable$lhs == lv)
    # leave the composite alone if the user attached any metadata to its
    # weights (in any block), or references them in constraints
    if (user_modified_rows(all_rows) || plabel_referenced(all_rows)) {
      next
    }
    # only consider one block for the indicator list (the structure is
    # identical across blocks)
    b1 <- lavpartable$block[all_rows]
    lv_rows <- all_rows[b1 == b1[1L]]
    ind <- lavpartable$rhs[lv_rows]

    # only observed indicators, and at least two of them
    keep <- ind %in% ov_names
    lv_rows <- lv_rows[keep]
    ind <- ind[keep]
    if (length(ind) < 2L) {
      next
    }

    # only adapt the 'clean' default situation (in every block): exactly one
    # indicator has a fixed weight (free == 0), it is the FIRST indicator,
    # and its value is 1.0; otherwise the user set up the scaling and we
    # leave it alone
    if (!clean_default_scaling(all_rows, ind)) {
      next
    }
    cur <- ind[1L]

    # implied relative weights, collected across all blocks that contain
    # all the indicators of this composite plus at least one other variable
    w_blocks <- list()
    for (b in seq_along(cov_list)) {
      vn <- block_names[[b]]
      if (is.null(vn) || !all(ind %in% vn)) {
        next
      }
      in_idx <- match(ind, vn) # the h1 cov matrices may have no dimnames
      out_idx <- setdiff(seq_along(vn), in_idx)
      if (length(out_idx) == 0L) {
        next
      }
      C <- cov_list[[b]]
      d <- diag(C)
      if (any(!is.finite(d)) || any(d <= 0)) {
        next
      }
      # scale-free: work with correlations
      isd <- 1 / sqrt(d)
      r_jj <- C[in_idx, in_idx, drop = FALSE] *
              tcrossprod(isd[in_idx])
      r_jo <- C[in_idx, out_idx, drop = FALSE] *
              tcrossprod(isd[in_idx], isd[out_idx])
      # dominant direction of the cross-correlations (~ S_jj %*% w)
      sv <- try(svd(r_jo, nu = 1L, nv = 0L), silent = TRUE)
      if (inherits(sv, "try-error") || sv$d[1L] < .Machine$double.eps) {
        next
      }
      w <- try(solve(r_jj, sv$u[, 1L]), silent = TRUE)
      if (inherits(w, "try-error") || any(!is.finite(w))) {
        next
      }
      w_blocks[[length(w_blocks) + 1L]] <- abs(w) / max(abs(w))
    }
    if (length(w_blocks) == 0L) {
      next
    }
    r <- colMeans(do.call(rbind, w_blocks), na.rm = TRUE)
    names(r) <- ind
    if (all(is.na(r))) {
      next
    }

    r_cur <- r[[1L]]
    best <- which.max(r)
    # switch only if the first indicator has a (near) zero implied weight
    # AND a clearly better indicator exists
    if (!is.na(r_cur) && r_cur < threshold &&
        ind[best] != cur && r[[best]] >= threshold) {
      new_marker[lv] <- ind[best]
      info[[length(info) + 1L]] <- data.frame(
        lv = lv, type = "composite", old = cur, new = ind[best],
        r.old = round(r[[1L]], 3), r.new = round(r[[best]], 3),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(new_marker) == 0L) {
    return(NULL)
  }

  list(marker = new_marker, info = do.call(rbind, info))
}
