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

# find the best ordering of the columns in lambda_target, to minimize
# the difference with the reference lambda matrix (lambda_ref)
#
# return the optimal order of the column for lambda_target
lav_efa_find_best_order <- function(lambda_ref = NULL, lambda_target = NULL,
                                    crit = "rmse") {
  m <- ncol(lambda_ref)
  stopifnot(ncol(lambda_target) == m)

  # all possible permutation
  # FIXEM:we really need a more elegant way to find the permutations...
  tmp <- unname(as.matrix(expand.grid(rep(list(seq_len(m)), m))))
  # select only rows where all numbers appear
  perm <- tmp[apply(tmp, 1L, function(x) {length(unique(x)) == m}), ,
              drop = FALSE]

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
  tmp <- unname(as.matrix(expand.grid(rep(list(seq_len(m)), m))))
  perm <- tmp[apply(tmp, 1L, function(x) length(unique(x)) == m), ,
    drop = FALSE]

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
