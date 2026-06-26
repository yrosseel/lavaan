# utility function related to EFA

# generate 'efa' syntax for a single block of factors
# generate EFA model syntax
#
# 'nfactors' is a per-block vector of (integer) factor counts; the blocks are
# ordered group-major, level-minor, ie block = (g - 1L) * nlevels + level.
# When there is a single block (ngroups == 1 and nlevels == 1), the output is
# just the '=~' lines (no 'group:'/'level:' header), identical to the original
# single-block syntax. With multiple levels and/or groups, the blocks are
# wrapped in 'level:' and/or 'group:' sections.
lav_syntax_efa <- function(ov_names = NULL, nfactors = 1L,
                           nlevels = 1L, ngroups = 1L) {
  # 'twolevel = TRUE' was the old (pre per-block) calling convention
  if (isTRUE(nlevels)) {
    nlevels <- 2L
  }
  nblocks <- ngroups * nlevels

  # recycle a single factor count to all blocks
  if (length(nfactors) == 1L) {
    nfactors <- rep.int(nfactors, nblocks)
  }
  stopifnot(length(nfactors) == nblocks)

  # the '=~' lines for a single block with 'k' factors
  block_lines <- function(k) {
    vapply(seq_len(k), function(f) {
      paste0('efa("efa")*f', f, " =~ ", paste(ov_names, collapse = " + "))
    }, character(1L))
  }

  # single block: no header (backward compatible)
  if (nblocks == 1L) {
    return(block_lines(nfactors[1L]))
  }

  model <- character(0L)
  for (g in seq_len(ngroups)) {
    if (ngroups > 1L) {
      model <- c(model, paste("group:", g))
    }
    for (l in seq_len(nlevels)) {
      block <- (g - 1L) * nlevels + l
      if (nlevels > 1L) {
        model <- c(model, paste("level:", l))
      }
      model <- c(model, block_lines(nfactors[block]))
    }
  }

  model
}

# extract *standardized* loadings from efaList
lav_efa_get_loadings <- function(object, ...) {
  # kill object$loadings if present
  object[["loadings"]] <- NULL

  out <- lapply(object, function(x) {
    std <- lavTech(x, "std",
      add.class = TRUE, add.labels = TRUE,
      list.by.group = FALSE
    )
    lambda_idx <- which(names(std) == "lambda")
    mm_lambda <- std[lambda_idx]
    names(mm_lambda) <- NULL
    # if only single block, drop list
    if (length(mm_lambda) == 1L) {
      mm_lambda <- mm_lambda[[1]]
    } else {
      names(mm_lambda) <- x@Data@block.label
    }
    mm_lambda
  })

  # drop list if only a single model
  if (length(out) == 1L) {
    out <- out[[1]]
  }

  out
}

# determine the column order for a (rotated or unrotated) loading matrix,
# according to 'order_lv_by': "sumofsquares", "index" (Asparouhov & Muthen
# 2009, Appendix D) or "none". Returns an integer permutation vector.
lav_efa_order_idx <- function(mm_lambda = NULL, order_lv_by = "index") {
  if (order_lv_by == "sumofsquares") {
    l2 <- mm_lambda * mm_lambda
    order_idx <- base::order(colSums(l2), decreasing = TRUE)
  } else if (order_lv_by == "index") {
    # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
    max_loading <- apply(abs(mm_lambda), 2, max)
    # 1: per factor, number of the loadings that are at least 0.8 of the
    #    highest loading of the factor
    # 2: mean of the index numbers
    average_index <- sapply(seq_len(ncol(mm_lambda)), function(i) {
      mean(which(abs(mm_lambda[, i]) >= 0.8 * max_loading[i]))
    })
    # order of the factors
    order_idx <- base::order(average_index)
  } else if (order_lv_by == "none") {
    order_idx <- seq_len(ncol(mm_lambda))
  } else {
    lav_msg_stop(gettext("order must be index, sumofsquares or none"))
  }
  order_idx
}

# reflect the columns of a loading matrix so that each column sum is
# non-negative; returns the (possibly sign-flipped) matrix.
lav_efa_reflect_cols <- function(mm_lambda = NULL) {
  neg_idx <- which(colSums(mm_lambda) < 0)
  if (length(neg_idx) > 0L) {
    mm_lambda[, neg_idx] <- -1 * mm_lambda[, neg_idx, drop = FALSE]
  }
  mm_lambda
}

# optionally reflect, then reorder the columns of a loading matrix. This is the
# common post-processing applied by the EFA extraction/PACE routines (the
# rotation engine applies the same ordering, but also needs to permute ROT/PHI,
# so it calls lav_efa_order_idx() directly).
lav_efa_reflect_and_reorder <- function(mm_lambda = NULL, reflect = TRUE,
                                        order_lv_by = "index") {
  if (reflect) {
    mm_lambda <- lav_efa_reflect_cols(mm_lambda)
  }
  order_idx <- lav_efa_order_idx(mm_lambda, order_lv_by)
  mm_lambda[, order_idx, drop = FALSE]
}

# find the best ordering of the columns in lambda_target, to minimize
# the difference with the reference lambda matrix (lambda_ref)
#
# enumerate all m! permutations of seq_len(m), one per row. The row order is the
# one produced by expand.grid(); both callers below rely on it for tie-breaking
# (they keep the first optimum), so do not change it.
lav_efa_permutations <- function(m = 1L) {
  tmp <- unname(as.matrix(expand.grid(rep(list(seq_len(m)), m))))
  # keep only rows where all numbers appear (i.e. actual permutations)
  tmp[apply(tmp, 1L, function(x) length(unique(x)) == m), , drop = FALSE]
}

# return the optimal order of the column for lambda_target
lav_efa_find_best_order <- function(lambda_ref = NULL, lambda_target = NULL,
                                    crit = "rmse") {
  m <- ncol(lambda_ref)
  stopifnot(ncol(lambda_target) == m)

  # all possible permutations
  perm <- lav_efa_permutations(m)

  rmse_perm <- numeric(nrow(perm))
  for (p in seq_len(nrow(perm))) {
    diff <- lambda_ref - lambda_target[, perm[p, ], drop = FALSE]
    diff2 <- diff * diff
    mse <- mean(diff2)
    rmse_perm[p] <- sqrt(mse)
  }
  best_idx <- which.min(rmse_perm)

  # return 'best' permutation
  perm[best_idx, ]
}

# find the signed permutation matrix T (m x m) such that 'lambda %*% T' best
# matches the reference loading matrix 'lambda_ref' (in a least-squares sense).
#
# this resolves the rotational indeterminacy of an EFA solution relative to a
# reference solution: factors may be re-ordered (column permutation) and/or
# reflected (sign change). It is used to align the rotated loadings of each
# bootstrap sample to the original (full-sample) rotated solution, so that the
# bootstrap distribution of the rotated parameters is meaningful.
#
# T is an orthonormal signed permutation matrix (T^-1 == t(T)), with a single
# non-zero (+1 or -1) entry per row/column. Applying it as an extra 'rotation'
# transforms all factor-related model matrices consistently (see
# lav_efa_bootstrap_align_coef()).
lav_efa_align_signed_perm <- function(lambda = NULL, lambda_ref = NULL) {
  m <- ncol(lambda)
  stopifnot(ncol(lambda_ref) == m)

  # single factor: only a possible sign flip
  if (m < 2L) {
    s <- if (sum(lambda * lambda_ref) < 0) -1 else 1
    return(matrix(s, 1L, 1L))
  }

  # all possible column permutations
  perm <- lav_efa_permutations(m)

  best_val <- Inf
  best_t <- diag(m)
  for (p in seq_len(nrow(perm))) {
    this_perm <- perm[p, ]
    lambda_p <- lambda[, this_perm, drop = FALSE]
    # optimal per-column sign: align each column with its reference column
    s <- sign(colSums(lambda_p * lambda_ref))
    s[s == 0] <- 1
    lambda_ps <- sweep(lambda_p, 2L, s, "*")
    val <- sum((lambda_ps - lambda_ref)^2)
    if (val < best_val) {
      best_val <- val
      tt <- matrix(0, m, m)
      for (j in seq_len(m)) {
        tt[this_perm[j], j] <- s[j]
      }
      best_t <- tt
    }
  }

  best_t
}
