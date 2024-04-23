# residual diagnostics

# two types:
# 1) residuals for summary statistics
# 2) case-wise residuals

# this (new) version written around Aug/Sept 2018 for 0.6-3
# - based on obsList (inspect_sampstat) and estList (inspect_implied)
# - pre-scaling for type = "cor.bollen" and type = "cor.bentler"
# - summary statistics: rmr, srmr, crmr, urmr, usrmr, ucrmr; standard errors,
#       confidence intervals (for u(cs)rmr),
#       z-statistics (exact test, close test), p-values
# - type = "normalized" is based on lav_model_h1_acov(), and should now work
#   for all estimators
# - type = "standardized" now uses the correct formula, and should work for
#   for all estimators
# - type = "standardized.mplus" uses the simplified Mplus/LISREL version,
#   often resulting in NAs due to negative var(resid) estimates
#   (this was "standardized" in lavaan < 0.6.3

# WARNING: only partial support for conditional.x = TRUE!!
# - in categorical case: we only compute summary statistics, using cor + th
#   (no var, slopes, ...)
# - twolevel not supported here; see lav_fit_srmr.R, where we convert to
#   the unconditional setting

# - change 0.6-6: we enforce observed.information = "h1" to ensure 'Q' is a
#                 projection matrix (see lav_residuals_acov)

# - change 0.6-13: fixed.x = TRUE is ignored (to conform with 'tradition')

setMethod(
  "residuals", "lavaan",
  function(object, type = "raw", labels = TRUE) {
    # lowercase type
    type <- tolower(type)

    # type = "casewise"
    if (type %in% c("casewise", "case", "obs", "observations", "ov")) {
      return(lav_residuals_casewise(object, labels = labels))
    } else {
      out <- lav_residuals(
        object = object, type = type, h1 = TRUE,
        add.type = TRUE,
        rename.cov.cor = FALSE, # should become FALSE!
        # after packages (eg jmv)
        # have adapted 0.6-3 style
        add.labels = labels, add.class = TRUE,
        drop.list.single.group = TRUE
      )
    }

    out
  }
)

setMethod(
  "resid", "lavaan",
  function(object, type = "raw") {
    residuals(object, type = type, labels = TRUE)
  }
)


# user-visible function
lavResiduals <- function(object, type = "cor.bentler", custom.rmr = NULL,
                         se = FALSE, zstat = TRUE, summary = TRUE,
                         h1.acov = "unstructured",
                         add.type = TRUE, add.labels = TRUE, add.class = TRUE,
                         drop.list.single.group = TRUE,
                         maximum.number = length(res.vech),
                         output = "list") {
  out <- lav_residuals(
    object = object, type = type, h1 = TRUE,
    custom.rmr = custom.rmr, se = se, zstat = zstat,
    summary = summary,
    summary.options = list(
      se = TRUE, zstat = TRUE, pvalue = TRUE,
      unbiased = TRUE, unbiased.se = TRUE, unbiased.ci = TRUE,
      unbiased.ci.level = 0.90, unbiased.zstat = TRUE,
      unbiased.test.val = 0.05, unbiased.pvalue = TRUE
    ),
    h1.acov = h1.acov, add.type = add.type,
    add.labels = add.labels, add.class = add.class,
    drop.list.single.group = drop.list.single.group
  )

  # no pretty printing yet...
  if (output == "table") {
    # new in 0.6-18: handle multiple blocks
    nblocks <- lav_partable_nblocks(object@ParTable)
    out.list <- vector("list", length = nblocks)
    for (block in seq_len(nblocks)) {
      if (nblocks == 1L) {
        res <- out$cov
      } else {
        res <- out[[block]]$cov
      }
      # extract only below-diagonal elements
      res.vech <- lav_matrix_vech(res, diagonal = FALSE)

      # get names
      P <- nrow(res)
      NAMES <- colnames(res)

      nam <- expand.grid(
        NAMES,
        NAMES
      )[lav_matrix_vech_idx(P, diagonal = FALSE), ]
      NAMES.vech <- paste(nam[, 1], "~~", nam[, 2], sep = "")

      # create table
      TAB <- data.frame(
        name = NAMES.vech, res = round(res.vech, 3),
        stringsAsFactors = FALSE
      )

      # sort table
      idx <- sort.int(abs(TAB$res),
        decreasing = TRUE,
        index.return = TRUE
      )$ix
      out.sorted <- TAB[idx, ]

      # show first rows only
      if (maximum.number == 0L || maximum.number > length(res.vech)) {
        maximum.number <- length(res.vech)
      }
      out.list[[block]] <- out.sorted[seq_len(maximum.number), ]
    }
    if (nblocks == 1L) {
      out <- out.list[[1]]
    } else {
      out <- out.list
      names(out) <- object@Data@block.label
    }
  } else {
    # list -> nothing to do
  }

  out
}

# main function
lav_residuals <- function(object, type = "raw", h1 = TRUE, custom.rmr = NULL,
                          se = FALSE, zstat = FALSE, summary = FALSE,
                          summary.options = list(
                            se = TRUE, zstat = TRUE,
                            pvalue = TRUE, unbiased = TRUE, unbiased.se = TRUE,
                            unbiased.ci = TRUE, unbiased.ci.level = 0.90,
                            unbiased.zstat = FALSE, unbiased.test.val = 0.05,
                            unbiased.pvalue = FALSE
                          ),
                          h1.acov = "unstructured", add.type = FALSE,
                          rename.cov.cor = FALSE,
                          add.labels = FALSE, add.class = FALSE,
                          drop.list.single.group = FALSE) {
  # type
  type <- tolower(type)[1]

  # check type
  if (!type %in% c(
    "raw", "cor", "cor.bollen", "cor.bentler", "cor.eqs",
    "rmr", "srmr", "crmr",
    "normalized", "standardized", "standardized.mplus"
  )) {
    lav_msg_stop(gettext("unknown argument for type:"), dQuote(type))
  }

  # if cor, choose 'default'
  if (type == "cor") {
    if (object@Options$mimic == "EQS") {
      type <- "cor.bentler"
    } else {
      type <- "cor.bollen"
    }
  }
  if (type == "cor.eqs") {
    type <- "cor.bentler"
  }
  if (type == "rmr") {
    type <- "raw"
  }
  if (type == "srmr") {
    type <- "cor.bentler"
  }
  if (type == "crmr") {
    type <- "cor.bollen"
  }

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model

  # change options if multilevel (for now)
  if (lavdata@nlevels > 1L) {
    zstat <- se <- FALSE
    summary <- FALSE
  }

  # change options if categorical (for now)
  if (lavmodel@categorical) {
    # only if conditional.x = FALSE AND no continuous endogenous variables
    # -> only the simple setting where we only have thresholds and
    #    correlations

    # As soon as we add continuous variables, we get means/variances too,
    # and we need to decide how WLS.obs/WLS.est/WLS.V will then map to
    # the output of lavInspect(fit, "implied") and
    # lavInspect(fit, "sampstat")

    if (lavmodel@conditional.x || length(unlist(lavmodel@num.idx)) > 0L) {
      zstat <- se <- FALSE
      summary <- FALSE
      summary.options <- list(
        se = FALSE, zstat = FALSE,
        pvalue = FALSE, unbiased = FALSE,
        unbiased.se = FALSE,
        unbiased.ci = FALSE, unbiased.ci.level = 0.90,
        unbiased.zstat = FALSE, unbiased.test.val = 0.05,
        unbiased.pvalue = FALSE
      )
    }
  }

  # change options if conditional.x (for now)
  if (!lavmodel@categorical && lavmodel@conditional.x) {
    zstat <- se <- FALSE
    summary <- FALSE
    summary.options <- list(
      se = FALSE, zstat = FALSE,
      pvalue = FALSE, unbiased = FALSE,
      unbiased.se = FALSE,
      unbiased.ci = FALSE, unbiased.ci.level = 0.90,
      unbiased.zstat = FALSE, unbiased.test.val = 0.05,
      unbiased.pvalue = FALSE
    )
  }

  # observed and fitted sample statistics
  obsList <- lav_object_inspect_sampstat(object,
    h1 = h1,
    add.labels = add.labels, add.class = add.class,
    drop.list.single.group = FALSE
  )
  estList <- lav_object_inspect_implied(object,
    add.labels = add.labels, add.class = add.class,
    drop.list.single.group = FALSE
  )
  # blocks
  nblocks <- length(obsList)

  # pre-scale?
  if (type %in% c("cor.bentler", "cor.bollen")) {
    for (b in seq_len(nblocks)) {
      var.obs <- if (lavmodel@conditional.x) {
        diag(obsList[[b]][["res.cov"]])
      } else {
        diag(obsList[[b]][["cov"]])
      }
      var.est <- if (lavmodel@conditional.x) {
        diag(estList[[b]][["res.cov"]])
      } else {
        diag(estList[[b]][["cov"]])
      }

      # rescale obsList
      obsList[[b]] <-
        lav_residuals_rescale(x = obsList[[b]], diag.cov = var.obs)
      # rescale estList
      if (type == "cor.bentler") { # use obsList
        estList[[b]] <-
          lav_residuals_rescale(x = estList[[b]], diag.cov = var.obs)
      } else if (type == "cor.bollen") { # use estList for COV only
        estList[[b]] <- lav_residuals_rescale(
          x = estList[[b]],
          diag.cov = var.est, diag.cov2 = var.obs
        )
      }
    }
  }

  # compute residuals: (observed - implied)
  resList <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    resList[[b]] <- lapply(seq_len(length(obsList[[b]])),
      FUN = function(el) {
        obsList[[b]][[el]] - estList[[b]][[el]]
      }
    )
    # always name the elements, even if add.labels = FALSE
    NAMES <- names(obsList[[b]])
    names(resList[[b]]) <- NAMES
  }

  # do we need seList?
  if (se || zstat) {
    seList <- lav_residuals_se(object,
      type = type, z.type = "standardized",
      h1.acov = h1.acov,
      add.class = add.class, add.labels = add.labels
    )
  } else if (type %in% c("normalized", "standardized", "standardized.mplus")) {
    seList <- lav_residuals_se(object,
      type = "raw", z.type = type,
      h1.acov = h1.acov,
      add.class = add.class, add.labels = add.labels
    )
  } else {
    seList <- NULL
  }

  # normalize/standardize?
  if (type %in% c("normalized", "standardized", "standardized.mplus")) {
    for (b in seq_len(nblocks)) {
      if (add.labels) {
        NAMES <- names(resList[[b]])
      }
      resList[[b]] <- lapply(seq_len(length(resList[[b]])),
        FUN = function(el) {
          A <- resList[[b]][[el]]
          B <- seList[[b]][[el]]
          near.zero.idx <- which(abs(A) < 1e-05)
          if (length(near.zero.idx) > 0L) {
            B[near.zero.idx] <- 1
          }
          A / B
        }
      )
      if (add.labels) {
        names(resList[[b]]) <- NAMES
      }
    }
  }

  # add se
  resList.orig <- resList
  if (se) {
    for (b in seq_len(nblocks)) {
      NAMES.res <- names(resList[[b]])
      NAMES.se <- paste0(NAMES.res, ".se")
      resList[[b]] <- c(resList[[b]], seList[[b]])
      names(resList[[b]]) <- c(NAMES.res, NAMES.se)
    }
  }


  # add zstat
  if (zstat) {
    for (b in seq_len(nblocks)) {
      NAMES.res <- names(resList[[b]])
      NAMES.z <- paste0(names(resList.orig[[b]]), ".z")
      tmp <- lapply(seq_len(length(resList.orig[[b]])),
        FUN = function(el) {
          A <- resList.orig[[b]][[el]]
          B <- seList[[b]][[el]]
          # NOTE: which threshold should we use?
          # used to be 1e-05
          # changed to 1e-04 in 0.6-4
          near.zero.idx <- which(abs(A) < 1e-04)
          if (length(near.zero.idx) > 0L) {
            # B[near.zero.idx] <- as.numeric(NA)
            B[near.zero.idx] <- 1.0
          }
          A / B
        }
      )
      resList[[b]] <- c(resList[[b]], tmp)
      names(resList[[b]]) <- c(NAMES.res, NAMES.z)
    }
  }

  # add summary statistics (rms, mabs)
  if (summary) {
    args <- c(
      list(
        object = object, type = type, h1.acov = h1.acov,
        add.class = add.class, custom.rmr = custom.rmr
      ),
      summary.options
    )
    sumStat <- do.call("lav_residuals_summary", args)
    for (b in seq_len(nblocks)) {
      NAMES <- names(resList[[b]])
      resList[[b]] <- c(resList[[b]], list(sumStat[[b]][[1]])) # only 1
      NAMES <- c(NAMES, "summary")
      names(resList[[b]]) <- NAMES
    }
  }

  # last: add type
  if (add.type) {
    for (b in seq_len(nblocks)) {
      NAMES <- names(resList[[b]])
      resList[[b]] <- c(type, resList[[b]])
      NAMES <- c("type", NAMES)
      names(resList[[b]]) <- NAMES
    }
  }

  # optional: rename 'cov' to 'cor' (if type = "cor")
  if (rename.cov.cor && type %in% c("cor.bentler", "cor.bollen")) {
    for (b in seq_len(nblocks)) {
      NAMES <- names(resList[[b]])
      NAMES <- gsub("cov", "cor", NAMES)
      names(resList[[b]]) <- NAMES
    }
  }


  # output
  OUT <- resList
  if (nblocks == 1L && drop.list.single.group) {
    OUT <- OUT[[1]]
  } else {
    if (lavdata@nlevels == 1L &&
      length(lavdata@group.label) > 0L) {
      names(OUT) <- unlist(lavdata@group.label)
    } else if (lavdata@nlevels > 1L &&
      length(lavdata@group.label) == 0L) {
      names(OUT) <- lavdata@level.label
    }
  }

  OUT
}

# return ACOV as list per group
lav_residuals_acov <- function(object, type = "raw", z.type = "standardized",
                               h1.acov = "unstructured") {
  # check type
  if (z.type %in% c("normalized", "standardized.mplus") && type != "raw") {
    lav_msg_stop(gettextf(
      "z.type = %1$s can only be used with type = %2$s",
      dQuote(z.type), dQuote("raw")))
  }

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model
  lavsamplestats <- object@SampleStats

  # return list per group
  ACOV.res <- vector("list", length = lavdata@ngroups)

  # compute ACOV for observed h1 sample statistics (ACOV == Gamma/N)
  if (!is.null(lavsamplestats@NACOV[[1]])) {
    NACOV.obs <- lavsamplestats@NACOV # if this changes, tag @TDJorgensen in commit message
    ACOV.obs <- lapply(NACOV.obs, function(x) x / lavsamplestats@ntotal)
  } else {
    ACOV.obs <- lav_model_h1_acov(
      lavobject = object,
      h1.information = h1.acov
    )
  }

  # shortcut for normalized
  if (z.type == "normalized") {
    ACOV.res <- ACOV.obs
    return(ACOV.res)
  } else {
    if (z.type == "standardized") {
      A1 <- lav_model_h1_information(object)
      if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
        # A1 is diagonal matrix
        A1 <- lapply(A1, diag)
      }
      if (type %in% c("cor.bentler", "cor.bollen")) {
        sampstat <- lavTech(object, "sampstat")
      }
    } else if (z.type == "standardized.mplus") {
      VCOV <- lavTech(object, "vcov")
    }
    DELTA <- lavTech(object, "delta")
  }

  # for each group, compute ACOV
  for (g in seq_len(lavdata@ngroups)) {
    # group weight
    gw <- object@SampleStats@nobs[[g]] / object@SampleStats@ntotal  # if this changes, tag @TDJorgensen in commit message

    if (z.type == "standardized.mplus") { # simplified formula
      # also used by LISREL?
      # see https://www.statmodel.com/download/StandardizedResiduals.pdf

      ACOV.est.g <- DELTA[[g]] %*% VCOV %*% t(DELTA[[g]])
      ACOV.res[[g]] <- ACOV.obs[[g]] - ACOV.est.g
    } else if (z.type == "standardized") {
      # see Ogasawara (2001) using Bentler & Dijkstra (1985) eq 1.7.4

      # NVarCov, but always 'not' robust
      #
      # new in 0.6-6: to ensure Q is a projection matrix, we
      #               force observed.information = "h1"
      #               (only needed if information is observed)
      this.options <- object@Options
      this.options$observed.information[1] <- "h1"
      A0.g.inv <- lav_model_information(
        lavmodel = lavmodel,
        lavsamplestats = object@SampleStats,
        lavdata = lavdata,
        lavcache = object@Cache,
        lavimplied = object@implied,
        lavh1 = object@h1,
        lavoptions = this.options,
        extra = FALSE,
        augmented = TRUE,
        inverted = TRUE,
        use.ginv = TRUE
      )

      ACOV.est.g <- gw * (DELTA[[g]] %*% A0.g.inv %*% t(DELTA[[g]]))
      Q <- diag(nrow = nrow(ACOV.est.g)) - ACOV.est.g %*% A1[[g]]
      ACOV.res[[g]] <- Q %*% ACOV.obs[[g]] %*% t(Q)

      # correct ACOV.res for type = "cor.bentler" or type = "cor.bollen"
      if (type == "cor.bentler") {
        if (lavmodel@categorical) {
          if (lavmodel@conditional.x ||
            length(unlist(lavmodel@num.idx)) > 0L) {
            lav_msg_stop(gettext(
              "SE for cor.bentler not available (yet) if categorical = TRUE, and
              conditional.x = TRUE OR some endogenous variables are continuous"))
          } else {
            # nothing to do, as we already are in correlation metric
          }
        } else {
          # Ogasawara (2001), eq (13), or
          # Maydeu-Olivares (2017), eq (16)
          COV <- if (lavmodel@conditional.x) {
            sampstat[[g]][["res.cov"]]
          } else {
            sampstat[[g]][["cov"]]
          }
          SS <- 1 / sqrt(diag(COV))
          tmp <- lav_matrix_vech(tcrossprod(SS))
          G.inv.sqrt <- diag(tmp, nrow = length(tmp))
          if (lavmodel@meanstructure) {
            GG <- lav_matrix_bdiag(
              diag(SS, nrow = length(SS)),
              G.inv.sqrt
            )
          } else {
            GG <- G.inv.sqrt
          }
          ACOV.res[[g]] <- GG %*% ACOV.res[[g]] %*% GG
        } # continuous
      } else if (type == "cor.bollen") {
        if (lavmodel@categorical) {
          if (lavmodel@conditional.x ||
            length(unlist(lavmodel@num.idx)) > 0L) {
            lav_msg_stop(gettext(
              "SE for cor.bentler not available (yet) if categorical = TRUE, and
              conditional.x = TRUE OR some endogenous variables are continuous"))
          } else {
            # nothing to do, as we already are in correlation metric
          }
        } else {
          # here we use the Maydeu-Olivares (2017) approach, see eq 17
          COV <- if (lavmodel@conditional.x) {
            sampstat[[g]][["res.cov"]]
          } else {
            sampstat[[g]][["cov"]]
          }
          F1 <- lav_deriv_cov2corB(COV)
          if (lavmodel@meanstructure) {
            SS <- 1 / sqrt(diag(COV))
            FF <- lav_matrix_bdiag(diag(SS, nrow = length(SS)), F1)
          } else {
            FF <- F1
          }
          ACOV.res[[g]] <- FF %*% ACOV.res[[g]] %*% t(FF)
        } # continuous
      } # cor.bollen
    } # z.type = "standardized"
  } # g

  ACOV.res
}

# return resList with 'se' values for each residual
lav_residuals_se <- function(object, type = "raw", z.type = "standardized",
                             h1.acov = "unstructured",
                             add.class = FALSE, add.labels = FALSE) {
  # slots
  lavdata <- object@Data
  lavmodel <- object@Model
  lavpta <- object@pta

  # return list per group
  seList <- vector("list", length = lavdata@ngroups)

  # get ACOV per group
  ACOV.res <- lav_residuals_acov(
    object = object, type = type,
    z.type = z.type, h1.acov = h1.acov
  )

  # labels
  if (add.labels) {
    ov.names <- object@pta$vnames$ov
    ov.names.res <- object@pta$vnames$ov.nox
    ov.names.x <- object@pta$vnames$ov.x
  }

  # for each group, compute 'se' values, and fill list
  for (g in seq_len(lavdata@ngroups)) {
    nvar <- object@pta$nvar[[g]] # block or group-based?
    diag.ACOV <- diag(ACOV.res[[g]])

    # take care of negative, or non-finite diag.ACOV elements
    diag.ACOV[!is.finite(diag.ACOV)] <- NA
    diag.ACOV[diag.ACOV < 0] <- NA

    # categorical
    if (lavmodel@categorical) {
      if (lavmodel@conditional.x ||
        length(unlist(lavmodel@num.idx)) > 0L) {
        lav_msg_stop(gettext("not ready yet!"))
      }

      # COR
      nth <- length(lavmodel@th.idx[[g]])
      tmp <- sqrt(diag.ACOV[-(1:nth)])
      cov.se <- lav_matrix_vech_reverse(tmp, diagonal = FALSE)

      # MEAN
      mean.se <- rep(as.numeric(NA), nth)

      # TH
      th.se <- sqrt(diag.ACOV[1:nth])

      if (add.class) {
        class(cov.se) <- c("lavaan.matrix.symmetric", "matrix")
        class(mean.se) <- c("lavaan.vector", "numeric")
        class(th.se) <- c("lavaan.vector", "numeric")
      }
      if (add.labels) {
        rownames(cov.se) <- colnames(cov.se) <- ov.names[[g]]
        names(mean.se) <- ov.names[[g]]
        names(th.se) <- lavpta$vnames$th.mean[[g]]
      }
      seList[[g]] <- list(
        cov.se = cov.se, mean.se = mean.se,
        th.se = th.se
      )


      # continuous -- single level
    } else if (lavdata@nlevels == 1L) {
      if (lavmodel@conditional.x) {
        lav_msg_stop(gettext("not ready yet"))
      } else {
        if (lavmodel@meanstructure) {
          tmp <- sqrt(diag.ACOV[-(1:nvar)])
          cov.se <- lav_matrix_vech_reverse(tmp)
          mean.se <- sqrt(diag.ACOV[1:nvar])
          if (add.class) {
            class(cov.se) <- c("lavaan.matrix.symmetric", "matrix")
            class(mean.se) <- c("lavaan.vector", "numeric")
          }
          if (add.labels) {
            rownames(cov.se) <- colnames(cov.se) <- ov.names[[g]]
            names(mean.se) <- ov.names[[g]]
          }
          seList[[g]] <- list(cov.se = cov.se, mean.se = mean.se)
        } else {
          cov.se <- lav_matrix_vech_reverse(sqrt(diag.ACOV))
          if (add.class) {
            class(cov.se) <- c("lavaan.matrix.symmetric", "matrix")
          }
          if (add.labels) {
            rownames(cov.se) <- colnames(cov.se) <- ov.names[[g]]
          }
          seList[[g]] <- list(cov.se = cov.se)
        }
      }

      # continuous -- multilevel
    } else if (lavdata@nlevels > 1L) {
      lav_msg_stop(gettext("not ready yet"))
    }
  } # g

  seList
}

# return summary statistics as list per group
lav_residuals_summary <- function(object, type = c("rmr", "srmr", "crmr"),
                                  h1.acov = "unstructured", custom.rmr = NULL,
                                  se = FALSE, zstat = FALSE, pvalue = FALSE,
                                  unbiased = FALSE, unbiased.se = FALSE,
                                  unbiased.ci = FALSE, unbiased.ci.level = 0.90,
                                  unbiased.zstat = FALSE,
                                  unbiased.test.val = 0.05,
                                  unbiased.pvalue = FALSE,
                                  add.class = FALSE) {
  # arguments
  if (length(custom.rmr)) {
    if (!is.list(custom.rmr)) lav_msg_stop(gettext("custom.rmr must be a list"))
    ## Each custom (S/C)RMR must have a unique name
    customNAMES <- names(custom.rmr)
    if (is.null(customNAMES)) lav_msg_stop(gettext(
      "custom.rmr list must have names"))
    if (length(unique(customNAMES)) < length(custom.rmr)) {
      lav_msg_stop(gettext(
        "custom.rmr must have a unique name for each summary"))
    }
    ## Each list must contain a list consisting of $cov and/or $mean (no $th yet)
    for (i in seq_along(custom.rmr)) {
      if (!is.list(custom.rmr[[i]])) {
        lav_msg_stop(gettext("Each element in custom.rmr must be a list"))
      }
      if (is.null(names(custom.rmr[[i]]))) {
        lav_msg_stop(gettext("The list in custom.rmr must have names"))
      }
      if (!all(names(custom.rmr[[i]]) %in% c("cov", "mean"))) {
        lav_msg_stop(gettext(
          'Elements in custom.rmr must be names "cov" and/or "mean"'))
      }
      ## below, verify dimensions match rmsList.g
    }
    # FIXME: blocks can have unique models, need another layer of lists
    #       between custom summaries and moments
  } else {
    customNAMES <- NULL
  }

  if (pvalue) {
    zstat <- TRUE
  }
  if (zstat) {
    se <- TRUE
  }
  if (unbiased.pvalue) {
    unbiased.zstat <- TRUE
  }
  if (unbiased.zstat) {
    unbiased.se <- TRUE
  }

  if (!all(type %in% c(
    "rmr", "srmr", "crmr",
    "raw", "cor.bentler", "cor.bollen"
  ))) {
    lav_msg_stop(gettext("unknown type:"), dQuote(type))
  }

  # change type name
  idx <- which(type == "raw")
  if (length(idx) > 0L) {
    type[idx] <- "rmr"
  }
  idx <- which(type == "cor.bentler")
  if (length(idx) > 0L) {
    type[idx] <- "srmr"
  }
  idx <- which(type == "cor.bollen")
  if (length(idx) > 0L) {
    type[idx] <- "crmr"
  }

  # slots
  lavdata <- object@Data
  lavmodel <- object@Model

  # fixed.x/conditional.x
  fixed.x <- lavmodel@fixed.x
  conditional.x <- lavmodel@conditional.x


  rmrFlag <- srmrFlag <- crmrFlag <- FALSE
  if ("rmr" %in% type || "raw" %in% type) {
    # FIXME: recursive call to lav_residuals() is summary = TRUE!!
    rmrList <- lav_residuals(object = object, type = "raw")
    if (se || unbiased) {
      rmrList.se <- lav_residuals_acov(
        object = object, type = "raw",
        z.type = "standardized",
        h1.acov = "unstructured"
      )
    }
  }
  if ("srmr" %in% type || "cor.bentler" %in% type || "cor" %in% type) {
    srmrList <- lav_residuals(object = object, type = "cor.bentler")
    if (se || unbiased) {
      srmrList.se <- lav_residuals_acov(
        object = object,
        type = "cor.bentler",
        z.type = "standardized",
        h1.acov = "unstructured"
      )
    }
  }
  if ("crmr" %in% type || "cor.bollen" %in% type) {
    crmrList <- lav_residuals(object = object, type = "cor.bollen")
    if (se || unbiased) {
      crmrList.se <- lav_residuals_acov(
        object = object,
        type = "cor.bollen",
        z.type = "standardized",
        h1.acov = "unstructured"
      )
    }
  }

  # return list per group
  sumStat <- vector("list", length = lavdata@ngroups)

  # for each group, compute ACOV
  for (g in seq_len(lavdata@ngroups)) {
    nvar <- object@pta$nvar[[g]] # block or group-based?

    # categorical single level
    if (lavdata@nlevels == 1L && lavmodel@categorical) {
      if ((se || unbiased) && (conditional.x ||
        length(unlist(lavmodel@num.idx)) > 0L)) {
        lav_msg_stop(gettext("not ready yet"))
      } else {
        # remove fixed.x elements:
        # seems like a good idea, but nobody likes it
        # nvar.x <- pstar.x <- 0L
        # if(lavmodel@fixed.x) {
        #     nvar.x <- lavmodel@nexo[g]
        #     pstar.x <- nvar.x * (nvar.x - 1) / 2 # note '-'
        # }

        OUT <- vector("list", length(type))
        names(OUT) <- type

        for (typ in seq_len(length(type))) {
          if (type[typ] == "rmr") {
            rmsList.g <- rmrList[[g]]
            if (se || unbiased) {
              rmsList.se.g <- rmrList.se[[g]]
            }
          } else if (type[typ] == "srmr") {
            rmsList.g <- srmrList[[g]]
            if (se || unbiased) {
              rmsList.se.g <- srmrList.se[[g]]
            }
          } else if (type[typ] == "crmr") {
            rmsList.g <- crmrList[[g]]
            if (se || unbiased) {
              rmsList.se.g <- crmrList.se[[g]]
            }
          }

          # COR
          nth <- length(lavmodel@th.idx[[g]])
          if (conditional.x) {
            STATS <- lav_matrix_vech(rmsList.g[["res.cov"]],
              diagonal = FALSE
            )
          } else {
            STATS <- lav_matrix_vech(rmsList.g[["cov"]],
              diagonal = FALSE
            )
          }

          # should pstar be p*(p+1)/2 or p*(p-1)/2
          # we use the first for SRMR and the latter for CRMR
          if (type[typ] == "crmr") {
            pstar <- length(STATS)
          } else {
            pstar <- length(STATS) + nvar
          }
          ACOV <- NULL
          if (se || unbiased) {
            ACOV <- rmsList.se.g[-seq_len(nth),
              -seq_len(nth),
              drop = FALSE
            ]
          }
          RMS.COR <- lav_residuals_summary_rms(
            STATS = STATS,
            ACOV = ACOV, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased.se = unbiased.se,
            unbiased.ci = unbiased.ci,
            unbiased.ci.level = unbiased.ci.level,
            unbiased.zstat = unbiased.zstat,
            unbiased.test.val = unbiased.test.val,
            unbiased.pvalue = unbiased.pvalue,
            pstar = pstar, type = type[typ]
          )


          # THRESHOLDS
          if (conditional.x) {
            STATS <- rmsList.g[["res.th"]]
          } else {
            STATS <- rmsList.g[["th"]]
          }
          pstar <- length(STATS)
          ACOV <- NULL
          if (se || unbiased) {
            ACOV <- rmsList.se.g[seq_len(nth),
              seq_len(nth),
              drop = FALSE
            ]
          }
          RMS.TH <- lav_residuals_summary_rms(
            STATS = STATS,
            ACOV = ACOV, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased.se = unbiased.se,
            unbiased.ci = unbiased.ci,
            unbiased.ci.level = unbiased.ci.level,
            unbiased.zstat = unbiased.zstat,
            unbiased.test.val = unbiased.test.val,
            unbiased.pvalue = unbiased.pvalue,
            pstar = pstar, type = type[typ]
          )

          # MEAN
          # STATS <- rmsList.g[["mean"]]
          STATS <- numeric(0L)
          pstar <- length(STATS)
          ACOV <- NULL
          if (se || unbiased) {
            # TODO: extract from rmsList.se.g
          }
          RMS.MEAN <- lav_residuals_summary_rms(
            STATS = STATS,
            ACOV = ACOV, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased.se = unbiased.se,
            unbiased.ci = unbiased.ci,
            unbiased.ci.level = unbiased.ci.level,
            unbiased.zstat = unbiased.zstat,
            unbiased.test.val = unbiased.test.val,
            unbiased.pvalue = unbiased.pvalue,
            pstar = pstar, type = type[typ]
          )

          # VAR (not ready yet)
          # STATS <- diag(rmsList.g[["cov"]])[lavmodel@num.idx[[g]]]
          STATS <- numeric(0L)
          pstar <- length(STATS)
          ACOV <- NULL
          if (se || unbiased) {
            # TODO: extract from rmsList.se.g
          }
          RMS.VAR <- lav_residuals_summary_rms(
            STATS = STATS,
            ACOV = ACOV, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased.se = unbiased.se,
            unbiased.ci = unbiased.ci,
            unbiased.ci.level = unbiased.ci.level,
            unbiased.zstat = unbiased.zstat,
            unbiased.test.val = unbiased.test.val,
            unbiased.pvalue = unbiased.pvalue,
            pstar = pstar, type = type[typ]
          )

          # TOTAL -- FIXME: for conditional.x ....
          if (conditional.x) {
            STATS <- c(
              lav_matrix_vech(rmsList.g[["res.cov"]],
                diagonal = FALSE
              ),
              rmsList.g[["res.th"]]
            )
          } else {
            STATS <- c(
              lav_matrix_vech(rmsList.g[["cov"]],
                diagonal = FALSE
              ),
              rmsList.g[["th"]]
            )
            # rmsList.g[["mean"]],
            # diag(rmsList.g[["cov"]])[lavmodel@num.idx[[g]]])
          }

          # should pstar be p*(p+1)/2 or p*(p-1)/2 for COV/COR?
          # we use the first for SRMR and the latter for CRMR
          if (type[typ] == "crmr") {
            pstar <- length(STATS)
          } else {
            pstar <- length(STATS) + nvar
          }

          # if(lavmodel@fixed.x) {
          #    pstar <- pstar - pstar.x
          # }

          ACOV <- NULL
          if (se || unbiased) {
            ACOV <- rmsList.se.g
          }
          RMS.TOTAL <- lav_residuals_summary_rms(
            STATS = STATS,
            ACOV = ACOV, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased.se = unbiased.se,
            unbiased.ci = unbiased.ci,
            unbiased.ci.level = unbiased.ci.level,
            unbiased.zstat = unbiased.zstat,
            unbiased.test.val = unbiased.test.val,
            unbiased.pvalue = unbiased.pvalue,
            pstar = pstar, type = type[typ]
          )

          TABLE <- as.data.frame(cbind(
            RMS.COR,
            RMS.TH,
            # RMS.MEAN,
            # RMS.VAR,
            RMS.TOTAL
          ))
          # colnames(TABLE) <- c("cor", "thresholds", "mean",
          #                     "var", "total")
          colnames(TABLE) <- c("cor", "thresholds", "total")
          if (add.class) {
            class(TABLE) <- c("lavaan.data.frame", "data.frame")
          }
          OUT[[typ]] <- TABLE
        } # type
      } # not conditional.x or mixed cat/con

      # continuous -- single level
    } else if (lavdata@nlevels == 1L) {
      if ((se || unbiased) && conditional.x) {
        lav_msg_stop(gettext("not ready yet"))
      } else {
        # nvar.x <- pstar.x <- 0L
        # if(lavmodel@fixed.x) {
        #    nvar.x <- lavmodel@nexo[g]
        #    pstar.x <- nvar.x * (nvar.x + 1) / 2
        # }

        OUT <- vector("list", length(type))
        names(OUT) <- type

        for (typ in seq_len(length(type))) {
          if (type[typ] == "rmr") {
            rmsList.g <- rmrList[[g]]
            if (se || unbiased) {
              rmsList.se.g <- rmrList.se[[g]]
            }
          } else if (type[typ] == "srmr") {
            rmsList.g <- srmrList[[g]]
            if (se || unbiased) {
              rmsList.se.g <- srmrList.se[[g]]
            }
          } else if (type[typ] == "crmr") {
            rmsList.g <- crmrList[[g]]
            if (se || unbiased) {
              rmsList.se.g <- crmrList.se[[g]]
            }
          }

          # COV
          if (conditional.x) {
            STATS <- lav_matrix_vech(rmsList.g[["res.cov"]])
          } else {
            STATS <- lav_matrix_vech(rmsList.g[["cov"]])
          }
          # pstar <- ( length(STATS) - pstar.x )
          pstar <- length(STATS)
          if (type[typ] == "crmr") {
            # pstar <- pstar - ( nvar - nvar.x )
            pstar <- pstar - nvar
          }

          ACOV <- NULL
          if (se || unbiased) {
            ACOV <- if (lavmodel@meanstructure) {
              rmsList.se.g[-seq_len(nvar),
                -seq_len(nvar),
                drop = FALSE
              ]
            } else {
              rmsList.se.g
            }
          }
          RMS.COV <- lav_residuals_summary_rms(
            STATS = STATS,
            ACOV = ACOV, se = se, zstat = zstat, pvalue = pvalue,
            unbiased = unbiased, unbiased.se = unbiased.se,
            unbiased.ci = unbiased.ci,
            unbiased.ci.level = unbiased.ci.level,
            unbiased.zstat = unbiased.zstat,
            unbiased.test.val = unbiased.test.val,
            unbiased.pvalue = unbiased.pvalue,
            pstar = pstar, type = type[typ]
          )

          # MEAN
          if (lavmodel@meanstructure) {
            if (conditional.x) {
              STATS <- rmsList.g[["res.int"]]
            } else {
              STATS <- rmsList.g[["mean"]]
            }
            # pstar <- ( length(STATS) - nvar.x )
            pstar <- length(STATS)
            ACOV <- NULL
            if (se || unbiased) {
              ACOV <- rmsList.se.g[seq_len(nvar),
                seq_len(nvar),
                drop = FALSE
              ]
            }
            RMS.MEAN <- lav_residuals_summary_rms(
              STATS = STATS,
              ACOV = ACOV,
              se = se, zstat = zstat, pvalue = pvalue,
              unbiased = unbiased, unbiased.se = unbiased.se,
              unbiased.ci = unbiased.ci,
              unbiased.ci.level = unbiased.ci.level,
              unbiased.zstat = unbiased.zstat,
              unbiased.test.val = unbiased.test.val,
              unbiased.pvalue = unbiased.pvalue,
              pstar = pstar, type = type[typ]
            )
          }

          # TOTAL
          if (lavmodel@meanstructure) {
            if (conditional.x) {
              STATS <- c(
                rmsList.g[["res.int"]],
                lav_matrix_vech(rmsList.g[["res.cov"]])
              )
            } else {
              STATS <- c(
                rmsList.g[["mean"]],
                lav_matrix_vech(rmsList.g[["cov"]])
              )
            }
            # pstar <- ( length(STATS) - ( pstar.x + nvar.x) )
            pstar <- length(STATS)
            if (type[typ] == "crmr") {
              # pstar <- pstar - ( nvar - nvar.x )
              pstar <- pstar - nvar
            }
            ACOV <- NULL
            if (se || unbiased) {
              ACOV <- rmsList.se.g
            }
            RMS.TOTAL <- lav_residuals_summary_rms(
              STATS = STATS,
              ACOV = ACOV,
              se = se, zstat = zstat, pvalue = pvalue,
              unbiased = unbiased, unbiased.se = unbiased.se,
              unbiased.ci = unbiased.ci,
              unbiased.ci.level = unbiased.ci.level,
              unbiased.zstat = unbiased.zstat,
              unbiased.test.val = unbiased.test.val,
              unbiased.pvalue = unbiased.pvalue,
              pstar = pstar, type = type[typ]
            )
          }

          # CUSTOM
          if (length(custom.rmr)) {
            if (lavmodel@fixed.x && !lavmodel@conditional.x) {
              ## save exogenous-variable indices, use to remove or set
              ## FALSE any moments that cannot have nonzero residuals
              x.idx <- which(rownames(rmsList.g$cov) %in% object@Data@ov.names.x[[g]])
            }

            RMS.CUSTOM.LIST <- vector("list", length(customNAMES))

            for (cus in customNAMES) {
              ## in case there is no meanstructure
              STATS <- NULL
              ACOV.idx <- NULL

              # MEANS?
              if (lavmodel@meanstructure) {
                if ("mean" %in% names(custom.rmr[[cus]])) {
                  ## if logical, save numeric indices
                  if (is.logical(custom.rmr[[cus]]$mean)) {
                    ## check length
                    if (length(custom.rmr[[cus]]$mean) != length(rmsList.g[["mean"]])) {
                      lav_msg_stop(gettextf("length(custom.rmr$%s$mean) must
                            match length(lavResiduals(fit)$mean)", cus))
                    }
                    ACOV.idx <- which(custom.rmr[[cus]]$mean)
                    if (lavmodel@fixed.x && !lavmodel@conditional.x) {
                      ACOV.idx[x.idx] <- FALSE
                    }
                  } else if (!is.numeric(custom.rmr[[cus]]$mean)) {
                    lav_msg_stop(gettextf("custom.rmr$%s$mean must contain
                                  logical or numeric indices.", cus))
                  } else {
                    ACOV.idx <- custom.rmr[[cus]]$mean
                    if (lavmodel@fixed.x && !lavmodel@conditional.x) {
                      ACOV.idx <- setdiff(ACOV.idx, x.idx)
                    }
                    ACOV.idx <- ACOV.idx[!is.na(ACOV.idx)] # necessary?
                    if (max(ACOV.idx) > length(rmsList.g[["mean"]])) {
                      lav_msg_stop(gettextf(
                        "custom.rmr$%1$s$mean[%2$s] is an out-of-bounds index",
                        cus, which.max(ACOV.idx))
                      )
                    }
                  }
                  STATS <- rmsList.g[["mean"]][ACOV.idx]
                }
              }
              # (CO)VARIANCES?
              if ("cov" %in% names(custom.rmr[[cus]])) {
                ## if numeric, create a logical matrix to obtain
                ## ACOV.idx and check for x.idx
                if (is.numeric(custom.rmr[[cus]]$cov)) {
                  cusCOV <- rmsList.g[["cov"]] == "start with all FALSE"
                  ## matrix of row/column indices?
                  if (length(dim(custom.rmr[[cus]]$cov))) {
                    if (max(custom.rmr[[cus]]$cov[, 1:2] > nrow(rmsList.g[["cov"]]))) {
                      lav_msg_stop(gettextf(
                        "numeric indices in custom.rmr$%1$s$cov cannot
                        exceed %2$s", cus, nrow(rmsList.g[["cov"]])))
                    }
                    for (RR in 1:nrow(custom.rmr[[cus]]$cov)) {
                      cusCOV[
                        custom.rmr[[cus]]$cov[RR, 1],
                        custom.rmr[[cus]]$cov[RR, 2]
                      ] <- TRUE
                    }
                  } else {
                    ## numeric-vector indices
                    if (max(custom.rmr[[cus]]$cov > length(rmsList.g[["cov"]]))) {
                      lav_msg_stop(gettextf(
                        "numeric indices in custom.rmr$%1$s$cov cannot
                        exceed %2$s", cus, length(rmsList.g[["cov"]])))
                    }
                    cusCOV[custom.rmr[[cus]]$cov] <- TRUE
                  }

                  ## numeric indices no longer needed, use logical
                  custom.rmr[[cus]]$cov <- cusCOV
                } else if (!is.logical(custom.rmr[[cus]]$cov)) {
                  lav_msg_stop(gettextf(
                    "custom.rmr$%s$cov must be a logical square matrix or a
                    numeric matrix of (row/column) indices.", cus))
                }

                ## check dimensions
                if (!all(dim(custom.rmr[[cus]]$cov) == dim(rmsList.g[["cov"]]))) {
                  lav_msg_stop(gettextf(
                    "dim(custom.rmr$%s$cov) must match
                    dim(lavResiduals(fit)$cov)", cus))
                }
                ## users can specify upper.tri or lower.tri indices
                custom.rmr[[cus]]$cov <- custom.rmr[[cus]]$cov | t(custom.rmr[[cus]]$cov)
                ## but ACOV refers to lower.tri indices
                custom.rmr[[cus]]$cov[upper.tri(custom.rmr[[cus]]$cov)] <- FALSE
                ## diagonal relevant?
                if (type[typ] == "crmr") diag(custom.rmr[[cus]]$cov) <- FALSE
                ## extract lower.tri indices
                vech.idx <- which(lav_matrix_vech(custom.rmr[[cus]]$cov))

                ## add residuals to STATS, indices to ACOV.idx
                STATS <- c(STATS, lav_matrix_vech(rmsList.g[["cov"]])[vech.idx])
                ACOV.idx <- c(ACOV.idx, vech.idx)
              }


              ## count residuals in summary (x.idx already removed)
              pstar <- length(STATS)

              ACOV <- NULL
              if (se || unbiased) {
                ACOV <- rmsList.se.g[ACOV.idx, ACOV.idx, drop = FALSE]
              }
              RMS.CUSTOM.LIST[[cus]] <- lav_residuals_summary_rms(
                STATS = STATS,
                ACOV = ACOV,
                se = se, zstat = zstat, pvalue = pvalue,
                unbiased = unbiased, unbiased.se = unbiased.se,
                unbiased.ci = unbiased.ci,
                unbiased.ci.level = unbiased.ci.level,
                unbiased.zstat = unbiased.zstat,
                unbiased.test.val = unbiased.test.val,
                unbiased.pvalue = unbiased.pvalue,
                pstar = pstar, type = type[typ]
              )

              # FIXME: update for categorical
            } # cus
            RMS.CUSTOM <- do.call(rbind, RMS.CUSTOM.LIST)
          } else {
            RMS.CUSTOM <- NULL
          }


          if (lavmodel@meanstructure) {
            TABLE <- as.data.frame(cbind(
              RMS.COV,
              RMS.MEAN,
              RMS.TOTAL,
              RMS.CUSTOM
            ))
            colnames(TABLE) <- c(
              "cov", "mean", "total",
              customNAMES
            )
          } else {
            TABLE <- as.data.frame(cbind(RMS.COV, RMS.CUSTOM))
            colnames(TABLE) <- c("cov", customNAMES)
          }
          if (add.class) {
            class(TABLE) <- c("lavaan.data.frame", "data.frame")
          }
          OUT[[typ]] <- TABLE
        } # type
      } # continuous, single-level, unconditional

      # continuous -- multilevel
    } else if (lavdata@nlevels > 1L) {
      lav_msg_stop(gettext("not ready yet"))
    }

    sumStat[[g]] <- OUT
  } # g

  sumStat
}


lav_residuals_summary_rms <- function(STATS = NULL, ACOV = NULL,
                                      se = FALSE,
                                      level = 0.90,
                                      zstat = FALSE, pvalue = FALSE,
                                      unbiased = FALSE,
                                      unbiased.se = FALSE,
                                      unbiased.ci = FALSE,
                                      unbiased.ci.level = 0.90,
                                      unbiased.zstat = FALSE,
                                      unbiased.test.val = 0.05,
                                      unbiased.pvalue = FALSE,
                                      pstar = 0, type = "rms") {
  OUT <- vector("list", length = 0L)


  # covariance matrix
  if (length(STATS) > 0L) {
    rms <- sqrt(sum(STATS * STATS) / pstar)
  } else {
    rms <- 0
    se <- unbiased <- zstat <- FALSE
  }

  # default is NULL
  rms.se <- rms.z <- rms.pvalue <- NULL
  urms <- urms.se <- urms.z <- urms.pvalue <- NULL
  urms.ci.lower <- urms.ci.upper <- NULL
  if (!unbiased.zstat) {
    unbiased.test.val <- NULL
  }

  if (se || unbiased) {
    TR2 <- sum(diag(ACOV %*% ACOV))
    TR1 <- sum(diag(ACOV))
    if (se) {
      rms.avar <- TR2 / (TR1 * 2 * pstar)
      if (!is.finite(rms.avar) || rms.avar < .Machine$double.eps) {
        rms.se <- as.numeric(NA)
      } else {
        rms.se <- sqrt(rms.avar)
      }
    }
  }

  if (zstat) {
    E.rms <- (sqrt(TR1 / pstar) * (4 * TR1 * TR1 - TR2) / (4 * TR1 * TR1))
    rms.z <- max((rms - E.rms), 0) / rms.se
    if (pvalue) {
      rms.pvalue <- 1 - pnorm(rms.z)
    }
  }

  if (unbiased) {
    T.cov <- as.numeric(crossprod(STATS))
    eVe <- as.numeric(t(STATS) %*% ACOV %*% STATS)
    k.cov <- 1 - (TR2 + 2 * eVe) / (4 * T.cov * T.cov)
    urms <- (1 / k.cov * sqrt(max((T.cov - TR1), 0) / pstar))
    if (unbiased.se) {
      urms.avar <- (1 / (k.cov * k.cov) * (TR2 + 2 * eVe) / (2 * pstar * T.cov))
      if (!is.finite(urms.avar) || urms.avar < .Machine$double.eps) {
        urms.se <- as.numeric(NA)
      } else {
        urms.se <- sqrt(urms.avar)
      }
      if (unbiased.ci) {
        a <- (1 - unbiased.ci.level) / 2
        a <- c(a, 1 - a)
        fac <- stats::qnorm(a)
        urms.ci.lower <- urms + urms.se * fac[1]
        urms.ci.upper <- urms + urms.se * fac[2]
      }
      if (unbiased.zstat) {
        urms.z <- (urms - unbiased.test.val) / urms.se
        if (unbiased.pvalue) {
          urms.pvalue <- 1 - pnorm(urms.z)
        }
      }
    }
  }

  # labels
  if (type == "rmr") {
    OUT <- list(
      rmr = rms, rmr.se = rms.se,
      rmr.exactfit.z = rms.z, rmr.exactfit.pvalue = rms.pvalue,
      urmr = urms, urmr.se = urms.se,
      urmr.ci.lower = urms.ci.lower,
      urmr.ci.upper = urms.ci.upper,
      urmr.closefit.h0.value = unbiased.test.val,
      urmr.closefit.z = urms.z,
      urmr.closefit.pvalue = urms.pvalue
    )
  } else if (type == "srmr") {
    OUT <- list(
      srmr = rms, srmr.se = rms.se,
      srmr.exactfit.z = rms.z, srmr.exactfit.pvalue = rms.pvalue,
      usrmr = urms, usrmr.se = urms.se,
      usrmr.ci.lower = urms.ci.lower,
      usrmr.ci.upper = urms.ci.upper,
      usrmr.closefit.h0.value = unbiased.test.val,
      usrmr.closefit.z = urms.z,
      usrmr.closefit.pvalue = urms.pvalue
    )
  } else if (type == "crmr") {
    OUT <- list(
      crmr = rms, crmr.se = rms.se,
      crmr.exactfit.z = rms.z, crmr.exactfit.pvalue = rms.pvalue,
      ucrmr = urms, ucrmr.se = urms.se,
      ucrmr.ci.lower = urms.ci.lower,
      ucrmr.ci.upper = urms.ci.upper,
      ucrmr.closefit.h0.value = unbiased.test.val,
      ucrmr.closefit.z = urms.z,
      ucrmr.closefit.pvalue = urms.pvalue
    )
  }

  unlist(OUT)
}

# generate summary statistics for the residuals
lav_residuals_summary_old <- function(resList = NULL,
                                      add.class = FALSE, add.labels = FALSE) {
  # per block
  nblocks <- length(resList)

  for (b in seq_len(nblocks)) {
    # create new list, including with summary statistics interleaved
    x <- vector("list", length = 0L)
    nel <- length(resList[[b]])
    NAMES <- names(resList[[b]])

    for (el in seq_len(nel)) {
      EL <- resList[[b]][[el]]
      if (!is.null(NAMES)) {
        NAME <- NAMES[el]
      }

      if (is.character(EL)) {
        new.x <- list(EL)
        if (add.labels) {
          names(new.x) <- "type"
        }
        x <- c(x, new.x)
      } else if (is.matrix(EL) && isSymmetric(EL)) {
        tmp <- na.omit(lav_matrix_vech(EL))
        rms <- sqrt(sum(tmp * tmp) / length(tmp))
        mabs <- mean(abs(tmp))
        tmp2 <- na.omit(lav_matrix_vech(EL, diagonal = FALSE))
        rms.nodiag <- sqrt(sum(tmp2 * tmp2) / length(tmp2))
        mabs.nodiag <- mean(abs(tmp2))
        cov.summary <- c(rms, rms.nodiag, mabs, mabs.nodiag)
        if (add.labels) {
          names(cov.summary) <-
            c("rms", "rms.nodiag", "mabs", "mabs.nodiag")
        }
        if (add.class) {
          class(cov.summary) <- c("lavaan.vector", "numeric")
        }
        new.x <- list(EL, cov.summary)
        if (add.labels && !is.null(NAMES)) {
          names(new.x) <- c(NAME, paste0(NAME, ".summary"))
        }
        x <- c(x, new.x)
      } else {
        tmp <- na.omit(EL)
        rms <- sqrt(sum(tmp * tmp) / length(tmp))
        mabs <- mean(abs(tmp))
        mean.summary <- c(rms, mabs)
        if (add.labels) {
          names(mean.summary) <- c("rms", "mabs")
        }
        if (add.class) {
          class(mean.summary) <- c("lavaan.vector", "numeric")
        }
        new.x <- list(EL, mean.summary)
        if (add.labels && !is.null(NAMES)) {
          names(new.x) <- c(NAME, paste0(NAME, ".summary"))
        }
        x <- c(x, new.x)
      }
    } # nel

    # fill in block including summary statistics
    resList[[b]] <- x
  } # nblocks

  resList
}


# x is a list with sample statistics (eg output of inspect(fit, "sampstat")
# y is another (possibly the same) list with sample statistics
#
# to avoid many 'NAs', we set the scale-factor to 1
# if the to-be-scaled value is < 1e-05 (in absolute value)
lav_residuals_rescale <- function(x, diag.cov = NULL, diag.cov2 = NULL) {
  if (is.null(diag.cov2)) {
    diag.cov2 <- diag.cov
  }

  # make sure we can take the sqrt and invert
  diag.cov[!is.finite(diag.cov)] <- NA
  diag.cov[diag.cov < .Machine$double.eps] <- NA
  scale.cov <- tcrossprod(1 / sqrt(diag.cov))

  # for the mean, we use diag.cov2
  diag.cov2[!is.finite(diag.cov2)] <- NA
  diag.cov2[diag.cov2 < .Machine$double.eps] <- NA
  scale.mean <- 1 / sqrt(diag.cov2)

  # rescale cov
  if (!is.null(x[["cov"]])) {
    # catch (near) zero elements in x$cov
    near.zero.idx <- which(abs(x[["cov"]]) < 1e-05)
    scale.cov[near.zero.idx] <- 1
    x[["cov"]][] <- x[["cov"]] * scale.cov
  }
  if (!is.null(x[["res.cov"]])) {
    # catch (near) zero elements in x$res.cov
    near.zero.idx <- which(abs(x[["res.cov"]]) < 1e-05)
    scale.cov[near.zero.idx] <- 1
    x[["res.cov"]][] <- x[["res.cov"]] * scale.cov
  }

  # rescale int/mean
  if (!is.null(x[["res.int"]])) {
    # catch (near) zero elements in x$res.int
    near.zero.idx <- which(abs(x[["res.int"]]) < 1e-05)
    scale.mean[near.zero.idx] <- 1
    x[["res.int"]] <- x[["res.int"]] * scale.mean
  }

  if (!is.null(x[["mean"]])) {
    # catch (near) zero elements in x$mean
    near.zero.idx <- which(abs(x[["mean"]]) < 1e-05)
    scale.mean[near.zero.idx] <- 1
    x[["mean"]] <- x[["mean"]] * scale.mean
  }

  # FIXME: do something sensible for th, slopes, ...

  x
}
