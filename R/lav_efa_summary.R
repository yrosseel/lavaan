# summary information for a single (lavaan) efa model
#
# workflow:
# - summary() first calls summary.efaList()
# - for each model, summary.efaList() calls lav_object_summary() with
#   efa = TRUE and efa.args
# - for each model, lav_object_summary() calls
#   lav_efa_summary(object, efa.args = efa.args) to populate the $efa slot


# efa summary for a single lavaan object
lav_efa_summary <- function(object,
                            efa.args = list(
                              lambda = TRUE,
                              theta = TRUE,
                              psi = TRUE,
                              eigenvalues = TRUE,
                              sumsq.table = TRUE,
                              lambda.structure = FALSE,
                              fs.determinacy = FALSE,
                              se = FALSE,
                              zstat = FALSE,
                              pvalue = FALSE
                            )) {
  stopifnot(inherits(object, "lavaan"))

  nblocks <- object@Model@nblocks
  orthogonal.flag <- object@Options$rotation.args$orthogonal

  # get standardized solution
  LAMBDA <- THETA <- PSI <- NULL
  STD <- lavTech(object, "std",
    add.class = TRUE, add.labels = TRUE,
    list.by.group = FALSE
  )
  lambda.idx <- which(names(STD) == "lambda")
  theta.idx <- which(names(STD) == "theta")
  psi.idx <- which(names(STD) == "psi")

  # LAMBDA
  LAMBDA <- STD[lambda.idx]
  names(LAMBDA) <- NULL

  # THETA
  THETA <- STD[theta.idx]
  # make THETA diagonal
  THETA <- lapply(seq_len(nblocks), function(b) {
    tmp <- diag(THETA[[b]])
    class(tmp) <- c("lavaan.vector", "numeric")
    tmp
  })

  # PSI
  PSI <- STD[psi.idx]
  names(PSI) <- NULL

  # eigenvalues correlation matrix
  std.ov <- object@Options$rotation.args$std.ov
  COV <- object@h1$implied$cov # h1
  if (std.ov) {
    COV <- lapply(COV, cov2cor)
  }
  eigvals <- NULL
  if (efa.args$eigenvalues) {
    eigvals <- lapply(seq_len(nblocks), function(b) {
      tmp <- eigen(COV[[b]], only.values = TRUE)$values
      names(tmp) <- paste("ev", 1:nrow(LAMBDA[[b]]), sep = "")
      class(tmp) <- c("lavaan.vector", "numeric")
      tmp
    })
  }

  fs.determinacy <- NULL
  # Note: these 'determinacy' values are only properly defined for the
  #       'regression' factor scores! (If we would apply the same formulas
  #       for Bartlett factor scores, we would obtain 1's!
  if (efa.args$fs.determinacy) {
    fs.determinacy <- lapply(seq_len(nblocks), function(b) {
      COR <- cov2cor(COV[[b]]) # just in case
      COR.inv <- try(solve(COR), silent = TRUE)
      if (inherits(COR.inv, "try-error")) {
        return(rep(as.numeric(NA), nrow(PSI[[b]])))
      }
      fs <- LAMBDA[[b]] %*% PSI[[b]] # factor structure
      out <- sqrt(diag(t(fs) %*% COR.inv %*% fs))
      class(out) <- c("lavaan.vector", "numeric")
      out
    })
  }

  # sum-of-squares table
  sumsq.table <- NULL
  if (efa.args$sumsq.table) {
    sumsq.table <- lapply(seq_len(nblocks), function(b) {
      nvar <- nrow(LAMBDA[[b]])
      nfactor <- ncol(LAMBDA[[b]])

      # sum of squares:
      # - if orthogonal, this is really the sum of the squared factor
      #   loadings
      # - if oblique, we need to take the correlation into account
      sumsq <- diag(PSI[[b]] %*% crossprod(LAMBDA[[b]]))

      # reorder
      if (nfactor > 1L) {
        # determine order
        order.idx <- sort.int(sumsq, decreasing = TRUE, index.return = TRUE)$ix
        # re-order from large to small
        sumsq <- sumsq[order.idx]
      }

      # Proportion 'explained' (= proportion of total sumsq)
      #   note: sum(sumsq) == sum(communalities)
      propexpl <- sumsq / sum(sumsq)

      # Proportion var (= sumsq/nvar)
      propvar <- sumsq / nrow(LAMBDA[[b]])

      # Cumulative var
      cumvar <- cumsum(propvar)

      # construct table
      tmp <- rbind(sumsq, propexpl, propvar, cumvar)

      # total + colnames
      if (nfactor > 1L) {
        # add total column
        tmp <- cbind(tmp, rowSums(tmp))
        tmp[4, ncol(tmp)] <- tmp[3, ncol(tmp)]
        colnames(tmp) <-
          c(
            colnames(LAMBDA[[b]])[order.idx],
            "total"
          )
      } else {
        colnames(tmp) <- colnames(LAMBDA[[b]])[1]
      }

      # rownames
      if (nfactor == 1L) {
        ssq.label <- "Sum of squared loadings"
      } else if (orthogonal.flag) {
        ssq.label <- "Sum of sq (ortho) loadings"
      } else {
        ssq.label <- "Sum of sq (obliq) loadings"
      }
      rownames(tmp) <- c(
        ssq.label,
        "Proportion of total",
        "Proportion var",
        "Cumulative var"
      )

      # class
      class(tmp) <- c("lavaan.matrix", "matrix")

      tmp
    })
  } # sumsq.table

  # (factor) structure coefficients
  if (efa.args$lambda.structure) {
    lambda.structure <- lapply(seq_len(nblocks), function(b) {
      tmp <- LAMBDA[[b]] %*% PSI[[b]]
      class(tmp) <- c("lavaan.matrix", "matrix")
      tmp
    })
  } else {
    lambda.structure <- NULL
  }

  # standard errors (if any)
  lambda.se <- theta.se <- psi.se <- NULL
  lambda.zstat <- theta.zstat <- psi.zstat <- NULL
  lambda.pval <- theta.pval <- psi.pval <- NULL
  if (object@Options$se != "none") {
    SE <- lavTech(object, "std.se",
      add.class = TRUE, add.labels = TRUE,
      list.by.group = FALSE
    )

    se.flag <- (efa.args$se || efa.args$zstat || efa.args$pvalue)

    # ALWAYS use lambda.se
    if (efa.args$lambda) {
      lambda.se <- SE[lambda.idx]
      names(lambda.se) <- NULL
    }

    # theta.se
    if (se.flag && efa.args$theta) {
      theta.se <- SE[theta.idx]
      # make theta.se diagonal
      theta.se <- lapply(seq_len(nblocks), function(b) {
        tmp <- diag(theta.se[[b]])
        class(tmp) <- c("lavaan.vector", "numeric")
        tmp
      })
    }

    # ALWAYS use psi.se
    if (efa.args$psi) {
      psi.se <- SE[psi.idx]
      names(psi.se) <- NULL
    }

    # compute zstat
    if (efa.args$zstat || efa.args$pvalue) {
      if (efa.args$lambda) {
        lambda.zstat <- lapply(seq_len(nblocks), function(b) {
          tmp.se <- lambda.se[[b]]
          tmp.se[tmp.se < sqrt(.Machine$double.eps)] <-
            as.numeric(NA)
          tmp <- LAMBDA[[b]] / tmp.se
          class(tmp) <- c("lavaan.matrix", "matrix")
          tmp
        })
      }
      if (efa.args$theta) {
        theta.zstat <- lapply(seq_len(nblocks), function(b) {
          tmp.se <- theta.se[[b]]
          tmp.se[tmp.se < sqrt(.Machine$double.eps)] <-
            as.numeric(NA)
          tmp <- THETA[[b]] / tmp.se
          class(tmp) <- c("lavaan.vector", "numeric")
          tmp
        })
      }
      if (efa.args$psi) {
        psi.zstat <- lapply(seq_len(nblocks), function(b) {
          tmp.se <- psi.se[[b]]
          tmp.se[tmp.se < sqrt(.Machine$double.eps)] <-
            as.numeric(NA)
          tmp <- PSI[[b]] / tmp.se
          class(tmp) <- c(
            "lavaan.matrix.symmetric",
            "matrix"
          )
          tmp
        })
      }
    }

    # compute pval
    if (efa.args$pvalue) {
      if (efa.args$lambda) {
        lambda.pval <- lapply(seq_len(nblocks), function(b) {
          tmp <- 2 * (1 - pnorm(abs(lambda.zstat[[b]])))
          class(tmp) <- c("lavaan.matrix", "matrix")
          tmp
        })
      }
      if (efa.args$theta) {
        theta.pval <- lapply(seq_len(nblocks), function(b) {
          tmp <- 2 * (1 - pnorm(abs(theta.zstat[[b]])))
          class(tmp) <- c("lavaan.vector", "numeric")
          tmp
        })
      }
      if (efa.args$psi) {
        psi.pval <- lapply(seq_len(nblocks), function(b) {
          tmp <- 2 * (1 - pnorm(abs(psi.zstat[[b]])))
          class(tmp) <- c(
            "lavaan.matrix.symmetric",
            "matrix"
          )
          tmp
        })
      }
    }
  } # se/zstat/pvalue

  # block.label
  block.label <- object@Data@block.label

  # we remove them here; we may have needed them for other parts
  if (!efa.args$lambda) {
    LAMBDA <- NULL
  }
  if (!efa.args$theta) {
    THETA <- NULL
  }
  if (!efa.args$psi) {
    PSI <- NULL
  }
  if (!efa.args$se) {
    # always keep lambda.se and psi.se (for the signif stars)
    theta.se <- NULL
  }
  if (!efa.args$zstat) {
    lambda.zstat <- theta.zstat <- psi.zstat <- NULL
  }

  res <- list(
    nblocks = nblocks,
    block.label = block.label,
    std.ov = std.ov,
    eigvals = eigvals,
    sumsq.table = sumsq.table,
    orthogonal = object@Options$rotation.args$orthogonal,
    lambda.structure = lambda.structure,
    fs.determinacy = fs.determinacy,
    lambda = LAMBDA,
    theta = THETA,
    psi = PSI,
    lambda.se = lambda.se,
    lambda.zstat = lambda.zstat,
    lambda.pvalue = lambda.pval,
    psi.se = psi.se,
    psi.zstat = psi.zstat,
    psi.pvalue = psi.pval,
    theta.se = theta.se,
    theta.zstat = theta.zstat,
    theta.pvalue = theta.pval
  )

  res
}


# summary efaList
summary.efaList <- function(object, nd = 3L, cutoff = 0.3, dot.cutoff = 0.1,
                            alpha.level = 0.01,
                            lambda = TRUE, theta = TRUE, psi = TRUE,
                            fit.table = TRUE, fs.determinacy = FALSE,
                            eigenvalues = TRUE, sumsq.table = TRUE,
                            lambda.structure = FALSE, se = FALSE,
                            zstat = FALSE, pvalue = FALSE, ...) {
  # kill object$loadings if present
  object[["loadings"]] <- NULL

  # unclass the object
  y <- unclass(object)

  # construct efa.args
  efa.args <- list(
    lambda = lambda, theta = theta, psi = psi,
    eigenvalues = eigenvalues, sumsq.table = sumsq.table,
    lambda.structure = lambda.structure,
    fs.determinacy = fs.determinacy,
    se = se, zstat = zstat, pvalue = pvalue
  )

  # extract useful info from first model
  out <- lav_object_summary(y[[1]],
    header = TRUE, estimates = FALSE,
    efa = FALSE
  )

  # header information
  lavaan.version <- out$header$lavaan.version
  converged.flag <- all(sapply(y, lavInspect, "converged"))

  # estimator
  estimator <- out$optim$estimator
  estimator.args <- out$optim$estimator.args

  # rotation
  rotation <- out$rotation$rotation
  rotation.args <- out$rotation$rotation.args

  # data
  lavdata <- out$data

  # main part: lav_object_summary information per model
  RES <- lapply(y, lav_object_summary,
    header = FALSE,
    fit.measures = FALSE, estimates = TRUE, efa = TRUE,
    efa.args = efa.args
  )

  # number of factors (for ALL blocks)
  nfactors <- sapply(y, function(x) x@pta$nfac[[1]])

  # fit.measures
  Table <- NULL
  if (fit.table) {
    # first, create standard table
    FIT <- fitMeasures(object, fit.measures = "default")
    NAMES <- rownames(FIT)
    idx <- integer(0L)

    # AIC/BIC
    if (all(c("aic", "bic", "bic2") %in% NAMES)) {
      this.idx <- match(c("aic", "bic", "bic2"), NAMES)
      idx <- c(idx, this.idx)
    }

    # chi-square
    if (all(c("chisq.scaled", "df.scaled", "pvalue.scaled") %in% NAMES)) {
      this.idx <- match(
        c("chisq.scaled", "df.scaled", "pvalue.scaled"),
        NAMES
      )
      idx <- c(idx, this.idx)
    } else {
      this.idx <- match(c("chisq", "df", "pvalue"), NAMES)
      idx <- c(idx, this.idx)
    }

    # CFI
    if ("cfi.robust" %in% NAMES && !all(is.na(FIT["cfi.robust", ]))) {
      this.idx <- match("cfi.robust", NAMES)
      idx <- c(idx, this.idx)
    } else if ("cfi.scaled" %in% NAMES) {
      this.idx <- match("cfi.scaled", NAMES)
      idx <- c(idx, this.idx)
    } else if ("cfi" %in% NAMES) {
      this.idx <- match("cfi", NAMES)
      idx <- c(idx, this.idx)
    }

    # RMSEA
    if ("rmsea.robust" %in% NAMES && !all(is.na(FIT["rmsea.robust", ]))) {
      this.idx <- match("rmsea.robust", NAMES)
      idx <- c(idx, this.idx)
    } else if ("rmsea.scaled" %in% NAMES) {
      this.idx <- match("rmsea.scaled", NAMES)
      idx <- c(idx, this.idx)
    } else if ("rmsea" %in% NAMES) {
      this.idx <- match("rmsea", NAMES)
      idx <- c(idx, this.idx)
    }

    # table with fitmeasures
    if (length(idx) > 0L) {
      Table <- t(FIT[idx, , drop = FALSE])
      tmp <- NAMES[idx]
      # strip '.scaled'
      tmp <- gsub(".scaled", "", tmp)
      # replace 'robust' by 'r' (if any)
      tmp <- gsub(".robust", "", tmp)
      # rename "bic2" -> "sabic"
      bic2.idx <- which(tmp == "bic2")
      if (length(bic2.idx) > 0L) {
        tmp[bic2.idx] <- "sabic"
      }
      colnames(Table) <- tmp
    } else {
      Table <- matrix(0, nrow = nfactors, ncol = 0L)
    }
    rownames(Table) <- paste("nfactors = ", nfactors, sep = "")
    class(Table) <- c("lavaan.matrix", "matrix")
  }

  # create return object
  out <- list(
    lavaan.version = lavaan.version,
    converged.flag = converged.flag,
    estimator = estimator,
    estimator.args = estimator.args,
    rotation = rotation,
    rotation.args = rotation.args,
    lavdata = lavdata,
    fit.table = Table,
    nfactors = nfactors,
    model.list = RES
  )

  # add nd, cutoff, dot.cutoff, ... as attributes (for printing)
  attr(out, "nd") <- nd
  attr(out, "cutoff") <- cutoff
  attr(out, "dot.cutoff") <- dot.cutoff
  attr(out, "alpha.level") <- alpha.level

  # create class
  class(out) <- c("efaList.summary", "list")

  out
}
