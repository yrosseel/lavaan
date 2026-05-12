# utility function related to EFA

# generate 'efa' syntax for a single block of factors
lav_syntax_efa <- function(ov_names = NULL, nfactors = 1L, twolevel = FALSE) {
  if (twolevel) {
    tmp <- lav_syntax_efa(ov_names = ov_names, nfactors = nfactors)
    model <- c("level: 1", tmp, "level: 2", tmp)
  } else {
    model <- character(nfactors)
    for (f in seq_len(nfactors)) {
      txt <- paste('efa("efa")*f', f, " =~ ",
        paste(ov_names, collapse = " + "),
        sep = ""
      )
      model[f] <- txt
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
