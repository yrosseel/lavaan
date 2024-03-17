# print only (standardized) loadings
print.lavaan.efa <- function(x, nd = 3L, cutoff = 0.3,
                             dot.cutoff = 0.1, alpha.level = 0.01, ...) {
  # unclass
  y <- unclass(x)

  if (!y$header$optim.converged) {
    cat("** WARNING ** Optimizer did not end normally\n")
    cat("** WARNING ** Estimates below are most likely unreliable\n")
  }

  # loadings per block
  for (b in seq_len(y$efa$nblocks)) {
    cat("\n")
    if (length(y$efa$block.label) > 0L) {
      cat(y$efa$block.label[[b]], ":\n\n", sep = "")
    }
    LAMBDA <- unclass(y$efa$lambda[[b]])
    lav_print_loadings(LAMBDA,
      nd = nd, cutoff = cutoff,
      dot.cutoff = dot.cutoff,
      alpha.level = alpha.level,
      x.se = y$efa$lambda.se[[b]]
    )
    cat("\n")
  }

  invisible(LAMBDA)
}

# print efaList
print.efaList <- function(x, nd = 3L, cutoff = 0.3,
                          dot.cutoff = 0.1, alpha.level = 0.01, ...) {
  # unclass
  y <- unclass(x)

  # kill loadings element if present
  y[["loadings"]] <- NULL

  nfits <- length(y)
  RES <- vector("list", nfits)
  for (ff in seq_len(nfits)) {
    res <- lav_object_summary(y[[ff]],
      fit.measures = FALSE,
      estimates = FALSE,
      modindices = FALSE,
      efa = TRUE,
      efa.args = list(
        lambda           = TRUE,
        theta            = FALSE,
        psi              = FALSE,
        eigenvalues      = FALSE,
        sumsq.table      = FALSE,
        lambda.structure = FALSE,
        fs.determinacy   = FALSE,
        se               = FALSE,
        zstat            = FALSE,
        pvalue           = FALSE
      )
    )
    RES[[ff]] <- print.lavaan.efa(res,
      nd = nd, cutoff = cutoff,
      dot.cutoff = dot.cutoff,
      alpha.level = alpha.level, ...
    )
  }

  invisible(RES)
}


# print summary efaList
print.efaList.summary <- function(x, nd = 3L, cutoff = 0.3,
                                  dot.cutoff = 0.1, alpha.level = 0.01,
                                  ...) {
  # unclass
  y <- unclass(x)

  # get nd, if it is stored as an attribute
  ND <- attr(y, "nd")
  if (!is.null(ND) && is.numeric(ND)) {
    nd <- as.integer(ND)
  }
  # get cutoff, if it is stored as an attribute
  CT <- attr(y, "cutoff")
  if (!is.null(CT) && is.numeric(CT)) {
    cutoff <- CT
  }
  # get dot.cutoff, if it is stored as an attribute
  DC <- attr(y, "dot.cutoff")
  if (!is.null(DC) && is.numeric(DC)) {
    dot.cutoff <- DC
  }
  # get alpha.level, if it is stored as an attribute
  AL <- attr(y, "alpha.level")
  if (!is.null(AL) && is.numeric(AL)) {
    alpha.level <- AL
  }

  cat("This is ",
    sprintf("lavaan %s", x$lavaan.version),
    " -- running exploratory factor analysis\n",
    sep = ""
  )

  # everything converged?
  if (!x$converged.flag) {
    cat("lavaan WARNING: not all models did converge!\n")
  }
  cat("\n")


  # estimator
  c1 <- c("Estimator")
  # second column
  tmp.est <- toupper(x$estimator)
  if (tmp.est == "DLS") {
    dls.first.letter <- substr(
      x$estimator.args$dls.GammaNT,
      1L, 1L
    )
    tmp.est <- paste("DLS-", toupper(dls.first.letter), sep = "")
  }
  c2 <- tmp.est

  # additional estimator args
  if (!is.null(x$estimator.args) &&
    length(x$estimator.args) > 0L) {
    if (x$estimator == "DLS") {
      c1 <- c(c1, "Estimator DLS value for a")
      c2 <- c(c2, x$estimator.args$dls.a)
    }
  }

  # rotation method
  c1 <- c(c1, "Rotation method")
  if (x$rotation == "none") {
    MM <- toupper(x$rotation)
  } else if (x$rotation.args$orthogonal) {
    MM <- paste(toupper(x$rotation), " ", "ORTHOGONAL",
      sep = ""
    )
  } else {
    MM <- paste(toupper(x$rotation), " ", "OBLIQUE",
      sep = ""
    )
  }
  c2 <- c(c2, MM)

  if (x$rotation != "none") {
    # method options
    if (x$rotation == "geomin") {
      c1 <- c(c1, "Geomin epsilon")
      c2 <- c(c2, x$rotation.args$geomin.epsilon)
    } else if (x$rotation == "orthomax") {
      c1 <- c(c1, "Orthomax gamma")
      c2 <- c(c2, x$rotation.args$orthomax.gamma)
    } else if (x$rotation == "cf") {
      c1 <- c(c1, "Crawford-Ferguson gamma")
      c2 <- c(c2, x$rotation.args$cf.gamma)
    } else if (x$rotation == "oblimin") {
      c1 <- c(c1, "Oblimin gamma")
      c2 <- c(c2, x$rotation.args$oblimin.gamma)
    } else if (x$rotation == "promax") {
      c1 <- c(c1, "Promax kappa")
      c2 <- c(c2, x$rotation.args$promax.kappa)
    }

    # rotation algorithm
    c1 <- c(c1, "Rotation algorithm (rstarts)")
    tmp <- paste(toupper(x$rotation.args$algorithm),
      " (", x$rotation.args$rstarts, ")",
      sep = ""
    )
    c2 <- c(c2, tmp)

    # Standardized metric (or not)
    c1 <- c(c1, "Standardized metric")
    if (x$rotation.args$std.ov) {
      c2 <- c(c2, "TRUE")
    } else {
      c2 <- c(c2, "FALSE")
    }

    # Row weights
    c1 <- c(c1, "Row weights")
    tmp.txt <- x$rotation.args$row.weights
    c2 <- c(c2, paste(toupper(substring(tmp.txt, 1, 1)),
      substring(tmp.txt, 2),
      sep = ""
    ))
  }

  # format c1/c2
  c1 <- format(c1, width = 33L)
  c2 <- format(c2,
    width = 18L + max(0, (nd - 3L)) * 4L,
    justify = "right"
  )

  # create character matrix
  M <- cbind(c1, c2, deparse.level = 0)
  colnames(M) <- rep("", ncol(M))
  rownames(M) <- rep(" ", nrow(M))

  # print
  write.table(M, row.names = TRUE, col.names = FALSE, quote = FALSE)

  # data
  if (!is.null(x$lavdata)) {
    cat("\n")
    lav_data_print_short(x$lavdata, nd = nd)
  }

  # number of models
  nfits <- length(x$model.list)

  # number of factors
  nfactors <- x$nfactors

  # fit measures
  if (!is.null(x$fit.table)) {
    cat("\n")
    if (nfits > 1L) {
      cat("Overview models:\n")
    } else {
      cat("Fit measures:\n")
    }
    print(x$fit.table, nd = nd, shift = 2L)
  }

  # eigenvalues
  if (!is.null(x$model.list[[1]]$efa$eigvals[[1]])) {
    cat("\n")
    if (x$model.list[[1]]$efa$std.ov) {
      cat("Eigenvalues correlation matrix:\n")
    } else {
      cat("Eigenvalues covariance matrix:\n")
    }
    for (b in seq_len(x$model.list[[1]]$efa$nblocks)) {
      cat("\n")
      if (length(x$model.list[[1]]$efa$block.label) > 0L) {
        cat(x$model.list[[1]]$efa$block.label[[b]], ":\n\n", sep = "")
      }
      print(x$model.list[[1]]$efa$eigvals[[b]], nd = nd, shift = 2L)
    } # blocks
  }

  # print summary for each model
  for (f in seq_len(nfits)) {
    res <- x$model.list[[f]]
    attr(res, "nd") <- nd
    attr(res, "cutoff") <- cutoff
    attr(res, "dot.cutoff") <- dot.cutoff
    attr(res, "alpha.level") <- alpha.level

    if (nfits > 1L) {
      if (f == 1L) {
        cat("\n")
      }
      cat("Number of factors: ", nfactors[f], "\n")
    }
    # print.lavaan.summary() prints the $efa element (only) or res
    print(res)
  }

  invisible(y)
}
