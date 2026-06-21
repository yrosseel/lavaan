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
# or not we have missing data, categorical data, etc.
#
# lav_pt_marker_adapt() returns NULL if no change is needed, or a list
# with:
#   - marker: a named character vector (lv -> new indicator) for the latent
#             variables whose marker should be switched
#   - info:   a data.frame (lv, old, new, r.old, r.new) for the warning
#
# 'threshold' is the value below which the first indicator is considered a
# 'poor' item (in absolute value of its corrected item-total correlation);
# it is set by the 'bad.marker.crit' option (default 0.1)
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
  # exclude efa and composite factors
  lv_regular <- lav_pt_vnames(lavpartable, type = "lv.regular")
  if (!is.null(lavpartable$efa)) {
    lv_efa <- unique(lavpartable$lhs[lavpartable$op == "=~" &
                                     nchar(lavpartable$efa) > 0L])
    lv_regular <- lv_regular[!lv_regular %in% lv_efa]
  }
  if (length(lv_regular) == 0L) {
    return(NULL)
  }

  # corrected item-total correlations, given a covariance submatrix C
  cit_from_cov <- function(C) {
    d <- diag(C)
    R <- rowSums(C)
    S <- sum(C)
    cov_rest <- R - d
    var_rest <- S - 2 * R + d
    out <- cov_rest / sqrt(d * var_rest)
    out[!is.finite(out)] <- NA_real_
    out
  }

  new_marker <- character(0L)
  info <- list()

  for (lv in lv_regular) {
    lv_rows <- which(lavpartable$op == "=~" & lavpartable$lhs == lv)
    # only consider one block (the structure is identical across blocks)
    b1 <- lavpartable$block[lv_rows]
    lv_rows <- lv_rows[b1 == b1[1L]]
    ind <- lavpartable$rhs[lv_rows]

    # only observed indicators, and at least two of them
    keep <- ind %in% ov_names
    lv_rows <- lv_rows[keep]
    ind <- ind[keep]
    if (length(ind) < 2L) {
      next
    }

    # only adapt the 'clean' default situation: exactly one indicator has a
    # fixed loading (free == 0), it is the FIRST indicator, and its value is
    # 1.0; otherwise the user set up the scaling and we leave it alone
    fixed <- which(lavpartable$free[lv_rows] == 0L)
    if (length(fixed) != 1L || fixed != 1L ||
        !isTRUE(lavpartable$ustart[lv_rows][1L] == 1)) {
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
        lv = lv, old = cur, new = ind[best],
        r.old = round(cit[[1L]], 3), r.new = round(cit[[best]], 3),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(new_marker) == 0L) {
    return(NULL)
  }

  list(marker = new_marker, info = do.call(rbind, info))
}
