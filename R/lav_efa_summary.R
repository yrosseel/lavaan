# summary information for a single (lavaan) efa model
#
# workflow:
# - summary() first calls lav_efalist_summary()
# - for each model, lav_efalist_summary() calls lav_object_summary() with
#   efa = TRUE and efa.args
# - for each model, lav_object_summary() calls
#   lav_efa_summary(object, efa.args = efa.args) to populate the $efa slot


# efa summary for a single lavaan object
lav_efa_summary <- function(object,
                            efa_args = list(
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
  orthogonal_flag <- object@Options$rotation.args$orthogonal

  # get standardized solution
  mm_lambda <- mm_theta <- mm_psi <- NULL
  std <- lavTech(object, "std",
    add.class = TRUE, add.labels = TRUE,
    list.by.group = FALSE
  )
  lambda_idx <- which(names(std) == "lambda")
  theta_idx <- which(names(std) == "theta")
  psi_idx <- which(names(std) == "psi")

  # LAMBDA
  mm_lambda <- std[lambda_idx]
  names(mm_lambda) <- NULL

  # THETA
  mm_theta <- std[theta_idx]
  # make THETA diagonal
  mm_theta <- lapply(seq_len(nblocks), function(b) {
    tmp <- diag(mm_theta[[b]])
    class(tmp) <- c("lavaan.vector", "numeric")
    tmp
  })

  # PSI
  mm_psi <- std[psi_idx]
  names(mm_psi) <- NULL

  # eigenvalues correlation matrix
  std_ov <- object@Options$rotation.args$std_ov
  cov_1 <- object@h1$implied$cov # h1
  if (std_ov) {
    cov_1 <- lapply(cov_1, cov2cor)
  }
  eigvals <- NULL
  if (efa_args$eigenvalues) {
    eigvals <- lapply(seq_len(nblocks), function(b) {
      tmp <- eigen(cov_1[[b]], only.values = TRUE)$values
      names(tmp) <- paste("ev", seq_len(nrow(mm_lambda[[b]])), sep = "")
      class(tmp) <- c("lavaan.vector", "numeric")
      tmp
    })
  }

  fs_determinacy <- NULL
  # Note: these 'determinacy' values are only properly defined for the
  #       'regression' factor scores! (If we would apply the same formulas
  #       for Bartlett factor scores, we would obtain 1's!
  if (efa_args$fs.determinacy) {
    fs_determinacy <- lapply(seq_len(nblocks), function(b) {
      cor_1 <- cov2cor(cov_1[[b]]) # just in case
      cor_inv <- try(solve(cor_1), silent = TRUE)
      if (inherits(cor_inv, "try-error")) {
        return(rep(as.numeric(NA), nrow(mm_psi[[b]])))
      }
      fs <- mm_lambda[[b]] %*% mm_psi[[b]] # factor structure
      out <- sqrt(diag(t(fs) %*% cor_inv %*% fs))
      class(out) <- c("lavaan.vector", "numeric")
      out
    })
  }

  # sum-of-squares table
  sumsq_table <- NULL
  if (efa_args$sumsq.table) {
    sumsq_table <- lapply(seq_len(nblocks), function(b) {
      # nvar <- nrow(mm_lambda[[b]])
      nfactor <- ncol(mm_lambda[[b]])

      # sum of squares:
      # - if orthogonal, this is really the sum of the squared factor
      #   loadings
      # - if oblique, we need to take the correlation into account
      sumsq <- diag(mm_psi[[b]] %*% crossprod(mm_lambda[[b]]))

      # reorder
      if (nfactor > 1L) {
        # determine order
        order_idx <- sort.int(sumsq, decreasing = TRUE, index.return = TRUE)$ix
        # re-order from large to small
        sumsq <- sumsq[order_idx]
      }

      # Proportion 'explained' (= proportion of total sumsq)
      #   note: sum(sumsq) == sum(communalities)
      propexpl <- sumsq / sum(sumsq)

      # Proportion var (= sumsq/nvar)
      propvar <- sumsq / nrow(mm_lambda[[b]])

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
            colnames(mm_lambda[[b]])[order_idx],
            "total"
          )
      } else {
        colnames(tmp) <- colnames(mm_lambda[[b]])[1]
      }

      # rownames
      if (nfactor == 1L) {
        ssq_label <- "Sum of squared loadings"
      } else if (orthogonal_flag) {
        ssq_label <- "Sum of sq (ortho) loadings"
      } else {
        ssq_label <- "Sum of sq (obliq) loadings"
      }
      rownames(tmp) <- c(
        ssq_label,
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
  if (efa_args$lambda.structure) {
    lambda_structure <- lapply(seq_len(nblocks), function(b) {
      tmp <- mm_lambda[[b]] %*% mm_psi[[b]]
      class(tmp) <- c("lavaan.matrix", "matrix")
      tmp
    })
  } else {
    lambda_structure <- NULL
  }

  # standard errors (if any)
  lambda_se <- theta_se <- psi_se <- NULL
  lambda_zstat <- theta_zstat <- psi_zstat <- NULL
  lambda_pval <- theta_pval <- psi_pval <- NULL
  if (object@Options$se != "none") {
    se <- lavTech(object, "std.se",
      add.class = TRUE, add.labels = TRUE,
      list.by.group = FALSE
    )

    se_flag <- (efa_args$se || efa_args$zstat || efa_args$pvalue)

    # ALWAYS use lambda.se
    if (efa_args$lambda) {
      lambda_se <- se[lambda_idx]
      names(lambda_se) <- NULL
    }

    # theta.se
    if (se_flag && efa_args$theta) {
      theta_se <- se[theta_idx]
      # make theta.se diagonal
      theta_se <- lapply(seq_len(nblocks), function(b) {
        tmp <- diag(theta_se[[b]])
        class(tmp) <- c("lavaan.vector", "numeric")
        tmp
      })
    }

    # ALWAYS use psi.se
    if (efa_args$psi) {
      psi_se <- se[psi_idx]
      names(psi_se) <- NULL
    }

    # compute zstat
    if (efa_args$zstat || efa_args$pvalue) {
      if (efa_args$lambda) {
        lambda_zstat <- lapply(seq_len(nblocks), function(b) {
          tmp_se <- lambda_se[[b]]
          tmp_se[tmp_se < sqrt(.Machine$double.eps)] <-
            as.numeric(NA)
          tmp <- mm_lambda[[b]] / tmp_se
          class(tmp) <- c("lavaan.matrix", "matrix")
          tmp
        })
      }
      if (efa_args$theta) {
        theta_zstat <- lapply(seq_len(nblocks), function(b) {
          tmp_se <- theta_se[[b]]
          tmp_se[tmp_se < sqrt(.Machine$double.eps)] <-
            as.numeric(NA)
          tmp <- mm_theta[[b]] / tmp_se
          class(tmp) <- c("lavaan.vector", "numeric")
          tmp
        })
      }
      if (efa_args$psi) {
        psi_zstat <- lapply(seq_len(nblocks), function(b) {
          tmp_se <- psi_se[[b]]
          tmp_se[tmp_se < sqrt(.Machine$double.eps)] <-
            as.numeric(NA)
          tmp <- mm_psi[[b]] / tmp_se
          class(tmp) <- c(
            "lavaan.matrix.symmetric",
            "matrix"
          )
          tmp
        })
      }
    }

    # compute pval
    if (efa_args$pvalue) {
      if (efa_args$lambda) {
        lambda_pval <- lapply(seq_len(nblocks), function(b) {
          tmp <- 2 * (1 - pnorm(abs(lambda_zstat[[b]])))
          class(tmp) <- c("lavaan.matrix", "matrix")
          tmp
        })
      }
      if (efa_args$theta) {
        theta_pval <- lapply(seq_len(nblocks), function(b) {
          tmp <- 2 * (1 - pnorm(abs(theta_zstat[[b]])))
          class(tmp) <- c("lavaan.vector", "numeric")
          tmp
        })
      }
      if (efa_args$psi) {
        psi_pval <- lapply(seq_len(nblocks), function(b) {
          tmp <- 2 * (1 - pnorm(abs(psi_zstat[[b]])))
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
  block_label <- object@Data@block.label

  # we remove them here; we may have needed them for other parts
  if (!efa_args$lambda) {
    mm_lambda <- NULL
  }
  if (!efa_args$theta) {
    mm_theta <- NULL
  }
  if (!efa_args$psi) {
    mm_psi <- NULL
  }
  if (!efa_args$se) {
    # always keep lambda.se and psi.se (for the signif stars)
    theta_se <- NULL
  }
  if (!efa_args$zstat) {
    lambda_zstat <- theta_zstat <- psi_zstat <- NULL
  }

  res <- list(
    nblocks = nblocks,
    block.label = block_label,
    std_ov = std_ov,
    eigvals = eigvals,
    sumsq.table = sumsq_table,
    orthogonal = object@Options$rotation.args$orthogonal,
    lambda.structure = lambda_structure,
    fs.determinacy = fs_determinacy,
    lambda = mm_lambda,
    theta = mm_theta,
    psi = mm_psi,
    lambda.se = lambda_se,
    lambda.zstat = lambda_zstat,
    lambda.pvalue = lambda_pval,
    psi.se = psi_se,
    psi.zstat = psi_zstat,
    psi.pvalue = psi_pval,
    theta.se = theta_se,
    theta.zstat = theta_zstat,
    theta.pvalue = theta_pval
  )

  res
}


# summary efaList
lav_efalist_summary <- function(object, nd = 3L, cutoff = 0.3, dot.cutoff = 0.1, # nolint start
                            alpha.level = 0.01,
                            lambda = TRUE, theta = TRUE, psi = TRUE,
                            fit.table = TRUE, fs.determinacy = FALSE,
                            eigenvalues = TRUE, sumsq.table = TRUE,
                            lambda.structure = FALSE, se = FALSE,
                            zstat = FALSE, pvalue = FALSE, ...) {               # nolint end
  # kill object$loadings if present
  object[["loadings"]] <- NULL

  # unclass the object
  y <- unclass(object)

  # construct efa.args
  efa_args <- list(
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
  lavaan_version <- out$header$lavaan.version
  converged_flag <- all(sapply(y, lavInspect, "converged"))

  # estimator
  estimator <- out$optim$estimator
  estimator_args <- out$optim$estimator.args

  # rotation
  rotation <- out$rotation$rotation
  rotation_args <- out$rotation$rotation.args

  # data
  lavdata <- out$data

  # main part: lav_object_summary information per model
  res <- lapply(y, lav_object_summary,
    header = FALSE,
    fit_measures = FALSE, estimates = TRUE, efa = TRUE,
    efa_args = efa_args
  )

  # number of factors (for ALL blocks)
  nfactors <- sapply(y, function(x) x@pta$nfac[[1]])

  # fit.measures
  table_1 <- NULL
  if (fit.table) {
    # first, create standard table
    fit <- fitMeasures(object, fit.measures = "default")
    names_1 <- rownames(fit)
    idx <- integer(0L)

    # AIC/BIC
    if (all(c("aic", "bic", "bic2") %in% names_1)) {
      this_idx <- match(c("aic", "bic", "bic2"), names_1)
      idx <- c(idx, this_idx)
    }

    # chi-square
    if (all(c("chisq.scaled", "df.scaled", "pvalue.scaled") %in% names_1)) {
      this_idx <- match(
        c("chisq.scaled", "df.scaled", "pvalue.scaled"),
        names_1
      )
      idx <- c(idx, this_idx)
    } else {
      this_idx <- match(c("chisq", "df", "pvalue"), names_1)
      idx <- c(idx, this_idx)
    }

    # CFI
    if ("cfi.robust" %in% names_1 && !all(is.na(fit["cfi.robust", ]))) {
      this_idx <- match("cfi.robust", names_1)
      idx <- c(idx, this_idx)
    } else if ("cfi.scaled" %in% names_1) {
      this_idx <- match("cfi.scaled", names_1)
      idx <- c(idx, this_idx)
    } else if ("cfi" %in% names_1) {
      this_idx <- match("cfi", names_1)
      idx <- c(idx, this_idx)
    }

    # RMSEA
    if ("rmsea.robust" %in% names_1 && !all(is.na(fit["rmsea.robust", ]))) {
      this_idx <- match("rmsea.robust", names_1)
      idx <- c(idx, this_idx)
    } else if ("rmsea.scaled" %in% names_1) {
      this_idx <- match("rmsea.scaled", names_1)
      idx <- c(idx, this_idx)
    } else if ("rmsea" %in% names_1) {
      this_idx <- match("rmsea", names_1)
      idx <- c(idx, this_idx)
    }

    # table with fitmeasures
    if (length(idx) > 0L) {
      table_1 <- t(fit[idx, , drop = FALSE])
      tmp <- names_1[idx]
      # strip '.scaled'
      tmp <- gsub(".scaled", "", tmp)
      # replace 'robust' by 'r' (if any)
      tmp <- gsub(".robust", "", tmp)
      # rename "bic2" -> "sabic"
      bic2_idx <- which(tmp == "bic2")
      if (length(bic2_idx) > 0L) {
        tmp[bic2_idx] <- "sabic"
      }
      colnames(table_1) <- tmp
    } else {
      table_1 <- matrix(0, nrow = nfactors, ncol = 0L)
    }
    rownames(table_1) <- paste("nfactors = ", nfactors, sep = "")
    class(table_1) <- c("lavaan.matrix", "matrix")
  }

  # create return object
  out <- list(
    lavaan.version = lavaan_version,
    converged.flag = converged_flag,
    estimator = estimator,
    estimator.args = estimator_args,
    rotation = rotation,
    rotation.args = rotation_args,
    lavdata = lavdata,
    fit.table = table_1,
    nfactors = nfactors,
    model.list = res
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
