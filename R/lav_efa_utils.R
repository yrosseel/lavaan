# utility function related to EFA

# generate 'efa' syntax for a single block of factors
lav_syntax_efa <- function(ov.names = NULL, nfactors = 1L, twolevel = FALSE) {
  if (twolevel) {
    tmp <- lav_syntax_efa(ov.names = ov.names, nfactors = nfactors)
    model <- c("level: 1", tmp, "level: 2", tmp)
  } else {
    model <- character(nfactors)
    for (f in seq_len(nfactors)) {
      txt <- paste('efa("efa")*f', f, " =~ ",
        paste(ov.names, collapse = " + "),
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
    STD <- lavTech(x, "std",
      add.class = TRUE, add.labels = TRUE,
      list.by.group = FALSE
    )
    lambda.idx <- which(names(STD) == "lambda")
    LAMBDA <- STD[lambda.idx]
    names(LAMBDA) <- NULL
    # if only single block, drop list
    if (length(LAMBDA) == 1L) {
      LAMBDA <- LAMBDA[[1]]
    } else {
      names(LAMBDA) <- x@Data@block.label
    }
    LAMBDA
  })

  # drop list if only a single model
  if (length(out) == 1L) {
    out <- out[[1]]
  }

  out
}
