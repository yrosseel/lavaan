# step 1 in SAM: fitting the measurement blocks

lav_sam_step1 <- function(cmd = "sem", mm.list = NULL, mm.args = list(),
                          FIT = FIT, data = NULL, sam.method = "local") {
  lavoptions <- FIT@Options
  lavpta <- FIT@pta
  PT <- FIT@ParTable
  nblocks <- lavpta$nblocks
  ngroups <- lavpta$ngroups

  if (lavoptions$verbose) {
    cat("Fitting the measurement part:\n")
  }

  # local only -> handle missing data
  if (sam.method %in% c("local", "fsr")) {
    # if missing = "listwise", make data complete, to avoid different
    # datasets per measurement block
    if (lavoptions$missing == "listwise") {
      # FIXME: make this work for multiple groups!!
      OV <- unique(unlist(FIT@pta$vnames$ov))
      # add group/cluster/sample.weights variables (if any)
      OV <- c(
        OV, FIT@Data@group, FIT@Data@cluster,
        FIT@Data@sampling.weights
      )
      data <- na.omit(data[, OV])
    }
  }

  # total number of free parameters
  if (FIT@Model@ceq.simple.only) {
    npar <- FIT@Model@nx.unco
    PT.free <- PT$free
    PT.free[PT.free > 0] <- seq_len(npar)
  } else {
    npar <- FIT@Model@nx.free
    PT.free <- PT$free
  }
  if (npar < 1L) {
    lav_msg_stop(gettext("model does not contain any free parameters"))
  }

  # check for higher-order factors
  # 0.6-11: hard stop for now, as we do not support them (yet)!
  LV.IND.names <- unique(unlist(FIT@pta$vnames$lv.ind))
  if (length(LV.IND.names) > 0L) {
    lav_msg_stop(gettext(
      "model contains indicators that are also latent variables:"),
      lav_msg_view(LV.IND.names, "none"))
    # ind.idx <- match(LV.IND.names, LV.names)
    # LV.names <- LV.names[-ind.idx]
  }

  # do we have at least 1 'regular' (measured) latent variable?
  LV.names <- unique(unlist(FIT@pta$vnames$lv.regular))
  if (length(LV.names) == 0L) {
    lav_msg_stop(gettext("model does not contain any (measured) latent
                         variables; use sem() instead"))
  }


  # how many measurement models?
  if (!is.null(mm.list)) {
    nMMblocks <- length(mm.list)
    # check each measurement block
    for (b in seq_len(nMMblocks)) {
      # check if we can find all lv names in LV.names
      if (!all(unlist(mm.list[[b]]) %in% LV.names)) {
        tmp <- unlist(mm.list[[b]])
        lav_msg_stop(gettext("mm.list contains unknown latent variable(s):"),
          lav_msg_view(tmp[!tmp %in% LV.names], "none"))
      }
      # make list per block
      if (!is.list(mm.list[[b]])) {
        mm.list[[b]] <- rep(list(mm.list[[b]]), nblocks)
      } else {
        if (length(mm.list[[b]]) != nblocks) {
          lav_msg_stop(gettextf(
            "mm.list block %1$s has length %2$s but nblocks = %3$s",
            b, length(mm.list[[b]]), nblocks))
        }
      }
    }
  } else {
    # TODO: here comes the automatic 'detection' of linked
    #       measurement models
    #
    # for now we take a single latent variable per measurement model block
    mm.list <- as.list(LV.names)
    nMMblocks <- length(mm.list)
    for (b in seq_len(nMMblocks)) {
      # make list per block
      mm.list[[b]] <- rep(list(mm.list[[b]]), nblocks)
    }
  }

  # adjust options for measurement models
  lavoptions.mm <- lavoptions
  lavoptions.mm$optim.bounds <- NULL
  if (lavoptions$se == "none") {
    lavoptions.mm$se <- "none"
  } else {
    # categorical?
    if (FIT@Model@categorical) {
      lavoptions.mm$se <- "robust.sem"
    } else if (lavoptions$estimator.orig == "MLM") {
      lavoptions.mm$se <- "robust.sem"
    } else if (lavoptions$estimator.orig == "MLR") {
      lavoptions.mm$se <- "robust.huber.white"
    } else if (lavoptions$estimator.orig == "PML") {
      lavoptions.mm$se <- "robust.huber.white"
    } else {
      lavoptions.mm$se <- "standard" # may be overriden later
    }
  }
  # if(sam.method == "global") {
  #    lavoptions.mm$test <- "none"
  # }
  # we need the tests to create summary info about MM
  lavoptions.mm$debug <- FALSE
  lavoptions.mm$verbose <- FALSE
  lavoptions.mm$check.post <- FALSE # neg lv variances may be overriden
  lavoptions.mm$check.gradient <- FALSE # too sensitive in large model (global)
  lavoptions.mm$baseline <- FALSE
  lavoptions.mm$bounds <- "wide.zerovar"

  # override with user-specified mm.args
  lavoptions.mm <- modifyList(lavoptions.mm, mm.args)

  # create MM slotOptions
  slotOptions.mm <- lav_options_set(lavoptions.mm)

  # we assume the same number/names of lv's per group!!!
  MM.FIT <- vector("list", nMMblocks) # fitted object

  # for joint model later
  if (lavoptions$se != "none") {
    Sigma.11 <- matrix(0, npar, npar)
  }
  step1.free.idx <- integer(0L)


  # NOTE: we should explicitly add zero-constrained LV covariances
  # to PT, and keep them zero in PTM
  if (cmd == "lavaan") {
    add.lv.cov <- FALSE
  } else {
    add.lv.cov <- TRUE
  }

  # fit mm model for each measurement block
  for (mm in seq_len(nMMblocks)) {
    if (lavoptions$verbose) {
      cat(
        "  block ", mm, "[",
        paste(mm.list[[mm]], collapse = " "), "]\n"
      )
    }

    # create parameter table for this measurement block only
    PTM <- lav_partable_subset_measurement_model(
      PT = PT,
      add.lv.cov = add.lv.cov,
      add.idx = TRUE,
      lv.names = mm.list[[mm]]
    )
    mm.idx <- attr(PTM, "idx")
    attr(PTM, "idx") <- NULL
    PTM$est <- NULL
    PTM$se <- NULL

    # check for categorical in PTM in this mm-block
    if (!any(PTM$op == "|")) {
      slotOptions.mm$categorical <- FALSE
      slotOptions.mm$.categorical <- FALSE
    }

    # update slotData for this measurement block
    ov.names.block <- lapply(1:ngroups, function(g) {
      unique(unlist(lav_partable_vnames(PTM, type = "ov", group = g)))
    })
    slotData.block <- lav_data_update_subset(FIT@Data,
      ov.names = ov.names.block
    )
    # handle single block 1-factor CFA with (only) two indicators
    if (length(unlist(ov.names.block)) == 2L && ngroups == 1L) {
      # hard stop for now, unless se = "none"
      if (lavoptions$se != "none") {
        lav_msg_stop(gettextf(
          "measurement block [%1$s", "] (%2$s) contains only two indicators;
          fix both factor loadings to unity, or combine factors into a single
          measurement block.", mm, lav_msg_view(mm.list[[mm]], "none")))
      } else {
        lambda.idx <- which(PTM$op == "=~")
        PTM$free[lambda.idx] <- 0L
        PTM$ustart[lambda.idx] <- 1
        PTM$start[lambda.idx] <- 1
        free.idx <- which(as.logical(PTM$free))
        if (length(free.idx) > 0L) {
          PTM$free[free.idx] <- seq_len(length(free.idx))
        }
        lav_msg_warn(gettextf(
          "measurement block [%1$s] (%2$s) contains only two indicators;
          -> fixing both factor loadings to unity",
          mm, lav_msg_view(mm.list[[mm]], "none")))
      }
    }

    # fit this measurement model only
    # (question: can we re-use even more slots?)
    fit.mm.block <- lavaan(
      model = PTM, slotData = slotData.block,
      slotOptions = slotOptions.mm
    )

    # check convergence
    if (!lavInspect(fit.mm.block, "converged")) {
      # warning for now, but this is not good!
      lav_msg_warn(gettextf(
        "measurement model for %s did not converge!",
        lav_msg_view(mm.list[[mm]], "none")))
    }

    # store fitted measurement model
    MM.FIT[[mm]] <- fit.mm.block

    # fill in point estimates measurement block (including slack values)
    PTM <- MM.FIT[[mm]]@ParTable
    # pt.idx: the row-numbers in PT that correspond to the rows in PTM
    # pt.idx <- lav_partable_map_id_p1_in_p2(p1 = PTM, p2 = PT,
    #             stopifnotfound = TRUE, exclude.nonpar = FALSE)
    # pt.idx == mm.idx
    ptm.idx <- which((PTM$free > 0L | PTM$op %in% c(":=", "<", ">")) &
      PTM$user != 3L)
    PT$est[mm.idx[ptm.idx]] <- PTM$est[ptm.idx]

    # if categorical, add non-free residual variances
    if (fit.mm.block@Model@categorical || fit.mm.block@Model@correlation) {
      extra.idx <- which(PTM$op %in% c("~~", "~*~") &
        PTM$lhs == PTM$rhs &
        PTM$user == 0L &
        PTM$free == 0L &
        PTM$ustart == 1)
      if (length(extra.idx) > 0L) {
        PT$est[mm.idx[extra.idx]] <- PTM$est[extra.idx]
      }
    }
    # if EFA, add user=7 values (but do not add to ptm.idx)
    user7.idx <- which(PTM$user == 7L)
    if (length(user7.idx)) {
      PT$est[mm.idx[user7.idx]] <- PTM$est[user7.idx]
    }

    # fill in standard errors measurement block
    if (lavoptions$se != "none") {
      if (fit.mm.block@Model@ceq.simple.only) {
        PTM.free <- PTM$free
        PTM.free[PTM.free > 0] <- seq_len(fit.mm.block@Model@nx.unco)
      } else {
        PTM.free <- PTM$free
      }

      ptm.se.idx <- which((PTM$free > 0L) & PTM$user != 3L) # no :=, <, >
      # PT$se[ seq_len(length(PT$lhs)) %in% mm.idx & PT$free > 0L ] <-
      #    PTM$se[ PTM$free > 0L & PTM$user != 3L]
      PT$se[mm.idx[ptm.se.idx]] <- PTM$se[ptm.se.idx]

      # compute variance matrix for this measurement block
      sigma.11 <- MM.FIT[[mm]]@vcov$vcov

      # fill in variance matrix
      par.idx <- PT.free[mm.idx[ptm.idx]]
      keep.idx <- PTM.free[ptm.idx]
      # par.idx <- PT.free[ seq_len(length(PT$lhs)) %in% mm.idx &
      #                    PT$free > 0L ]
      # keep.idx <- PTM.free[ PTM$free > 0 & PTM$user != 3L ]
      Sigma.11[par.idx, par.idx] <-
        sigma.11[keep.idx, keep.idx, drop = FALSE]

      # store (ordered) indices in step1.free.idx
      step1.free.idx <- c(step1.free.idx, sort.int(par.idx))
    }
  } # measurement block

  # only keep 'measurement part' parameters in Sigma.11
  if (lavoptions$se != "none") {
    Sigma.11 <- Sigma.11[step1.free.idx, step1.free.idx, drop = FALSE]
  } else {
    Sigma.11 <- NULL
  }

  # create STEP1 list
  STEP1 <- list(
    MM.FIT = MM.FIT, Sigma.11 = Sigma.11,
    step1.free.idx = step1.free.idx, PT.free = PT.free,
    mm.list = mm.list, PT = PT
  )

  STEP1
}

## STEP 1b: compute Var(eta) and E(eta) per block
##          only needed for local/fsr approach!
lav_sam_step1_local <- function(STEP1 = NULL, FIT = NULL,
                                sam.method = "local",
                                local.options = list(
                                  M.method = "ML",
                                  lambda.correction = TRUE,
                                  alpha.correction = 0L,
                                  twolevel.method = "h1"
                                )) {
  # local.M.method
  local.M.method <- toupper(local.options[["M.method"]])
  if (!local.M.method %in% c("GLS", "ML", "ULS")) {
    lav_msg_stop(gettext(
      "local option M.method should be one of GLS, ML or ULS."))
  }

  # local.twolevel.method
  local.twolevel.method <- tolower(local.options[["twolevel.method"]])
  if (!local.twolevel.method %in% c("h1", "anova", "mean")) {
    lav_msg_stop(gettext(
      "local option twolevel.method should be one of h1, anova or mean."))
  }

  lavoptions <- FIT@Options
  lavpta <- FIT@pta

  ngroups <- lavpta$ngroups
  nlevels <- lavpta$nlevels
  nblocks <- lavpta$nblocks
  nMMblocks <- length(STEP1$MM.FIT)
  mm.list <- STEP1$mm.list

  if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
    lv.interaction.flag <- TRUE
  } else {
    lv.interaction.flag <- FALSE
  }

  if (lavoptions$verbose) {
    cat("Constructing the mapping matrix using the ",
      local.M.method, " method ... ",
      sep = ""
    )
  }

  LAMBDA.list <- vector("list", nMMblocks)
  THETA.list <- vector("list", nMMblocks)
  NU.list <- vector("list", nMMblocks)
  DELTA.list <- vector("list", nMMblocks) # correlation/categorical
  LV.idx.list <- vector("list", nMMblocks)
  OV.idx.list <- vector("list", nMMblocks)
  for (mm in seq_len(nMMblocks)) {
    fit.mm.block <- STEP1$MM.FIT[[mm]]

    # LV.idx.list/OV.idx.list: list per block
    LV.idx.list[[mm]] <- vector("list", nblocks)
    OV.idx.list[[mm]] <- vector("list", nblocks)

    # store LAMBDA/THETA
    LAMBDA.list[[mm]] <- computeLAMBDA(fit.mm.block@Model)
    THETA.list[[mm]] <- computeTHETA(fit.mm.block@Model)
    if (fit.mm.block@Model@meanstructure) {
      NU.list[[mm]] <- computeNU(fit.mm.block@Model,
        lavsamplestats = fit.mm.block@SampleStats
      )
    }
    if (fit.mm.block@Model@categorical || fit.mm.block@Model@correlation) {
      delta.idx <- which(names(fit.mm.block@Model@GLIST) == "delta")
      DELTA.list[[mm]] <- fit.mm.block@Model@GLIST[delta.idx]
    }

    for (bb in seq_len(nblocks)) {
      lambda.idx <- which(names(FIT@Model@GLIST) == "lambda")[bb]
      ind.names <- fit.mm.block@pta$vnames$ov.ind[[bb]]
      LV.idx.list[[mm]][[bb]] <- match(
        mm.list[[mm]][[bb]],
        FIT@Model@dimNames[[lambda.idx]][[2]]
      )
      OV.idx.list[[mm]][[bb]] <- match(
        ind.names,
        FIT@Model@dimNames[[lambda.idx]][[1]]
      )
    } # nblocks
  } ## nMMblocks

  # assemble global LAMBDA/THETA (per block)
  LAMBDA <- computeLAMBDA(FIT@Model, handle.dummy.lv = FALSE)
  THETA <- computeTHETA(FIT@Model, fix = FALSE) # keep dummy lv
  if (FIT@Model@meanstructure) {
    NU <- computeNU(FIT@Model, lavsamplestats = FIT@SampleStats)
  }
  if (FIT@Model@categorical || FIT@Model@correlation) {
    delta.idx <- which(names(FIT@Model@GLIST) == "delta")
    DELTA <- FIT@Model@GLIST[delta.idx]
  }
  for (b in seq_len(nblocks)) {
    for (mm in seq_len(nMMblocks)) {
      ov.idx <- OV.idx.list[[mm]][[b]]
      lv.idx <- LV.idx.list[[mm]][[b]]
      LAMBDA[[b]][ov.idx, lv.idx] <- LAMBDA.list[[mm]][[b]]
      THETA[[b]][ov.idx, ov.idx] <- THETA.list[[mm]][[b]]
      # new in 0.6-10: check if any indicators are also involved
      # in the structural part; if so, set THETA row/col to zero
      # and make sure LAMBDA element is correctly set
      # (we also need to adjust M)
      dummy.ov.idx <- FIT@Model@ov.y.dummy.ov.idx[[b]]
      dummy.lv.idx <- FIT@Model@ov.y.dummy.lv.idx[[b]]
      if (length(dummy.ov.idx)) {
        THETA[[b]][dummy.ov.idx, ] <- 0
        THETA[[b]][, dummy.ov.idx] <- 0
        LAMBDA[[b]][dummy.ov.idx, ] <- 0
        LAMBDA[[b]][cbind(dummy.ov.idx, dummy.lv.idx)] <- 1
      }
      if (FIT@Model@meanstructure) {
        NU[[b]][ov.idx, 1] <- NU.list[[mm]][[b]]
        if (length(dummy.ov.idx)) {
          NU[[b]][dummy.ov.idx, 1] <- 0
        }
      }
      if ((FIT@Model@categorical || FIT@Model@correlation) &&
        !is.null(DELTA.list[[mm]][[b]])) { # could be mixed cat/cont
        DELTA[[b]][ov.idx, 1] <- DELTA.list[[mm]][[b]]
      }
    }

    # remove 'lv.interaction' columns from LAMBDA[[b]]
    # if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
    #  LAMBDA[[b]] <- LAMBDA[[b]][, -lavpta$vidx$lv.interaction[[b]]]
    # }

    # check if LAMBDA has full column rank
    this.lambda <- LAMBDA[[b]]
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      this.lambda <- this.lambda[, -lavpta$vidx$lv.interaction[[b]]]
    }
    if (qr(this.lambda)$rank < ncol(this.lambda)) {
      print(this.lambda)
      lav_msg_stop(gettext(
        "LAMBDA has no full column rank. Please use sam.method = global"))
    }
  } # b

  # store LAMBDA/THETA/NU per block
  STEP1$LAMBDA <- LAMBDA
  STEP1$THETA <- THETA
  if (FIT@Model@meanstructure) {
    STEP1$NU <- NU
  }
  if (FIT@Model@categorical || FIT@Model@correlation) {
    STEP1$DELTA <- DELTA
  }

  VETA <- vector("list", nblocks)
  REL <- vector("list", nblocks)
  alpha <- vector("list", nblocks)
  lambda <- vector("list", nblocks)
  if (lavoptions$meanstructure) {
    EETA <- vector("list", nblocks)
  } else {
    EETA <- NULL
  }
  M <- vector("list", nblocks)

  if (lv.interaction.flag) {
    # compute Bartlett factor scores
    FS <- vector("list", nblocks)
    # FS.mm <- lapply(STEP1$MM.FIT, lav_predict_eta_bartlett)
    FS.mm <- lapply(STEP1$MM.FIT, lavPredict,
      method = "Bartlett",
      drop.list.single.group = FALSE
    )
    for (b in seq_len(nblocks)) {
      tmp <- lapply(
        1:length(STEP1$MM.FIT),
        function(x) FS.mm[[x]][[b]]
      )
      LABEL <- unlist(lapply(tmp, colnames))
      FS[[b]] <- do.call("cbind", tmp)
      colnames(FS[[b]]) <- LABEL

      # dummy lv's?
      if (length(dummy.lv.idx) > 0L) {
        FS.obs <- FIT@Data@X[[b]][, dummy.ov.idx, drop = FALSE]
        colnames(FS.obs) <- FIT@Data@ov.names[[b]][dummy.ov.idx]
        FS[[b]] <- cbind(FS[[b]], FS.obs)
      }
    }
  }

  # compute VETA/EETA per block
  if (nlevels > 1L && local.twolevel.method == "h1") {
    H1 <- lav_h1_implied_logl(
      lavdata = FIT@Data,
      lavsamplestats = FIT@SampleStats,
      lavoptions = FIT@Options
    )
  }

  for (b in seq_len(nblocks)) {
    # get sample statistics for this block
    if (nlevels > 1L) {
      if (ngroups > 1L) {
        this.level <- (b - 1L) %% ngroups + 1L
      } else {
        this.level <- b
      }
      this.group <- floor(b / nlevels + 0.5)

      if (this.level == 1L) {
        if (local.twolevel.method == "h1") {
          COV <- H1$implied$cov[[1]]
          YBAR <- H1$implied$mean[[1]]
        } else if (local.twolevel.method == "anova" ||
          local.twolevel.method == "mean") {
          COV <- FIT@SampleStats@YLp[[this.group]][[2]]$Sigma.W
          YBAR <- FIT@SampleStats@YLp[[this.group]][[2]]$Mu.W
        }

        # reduce
        ov.idx <- FIT@Data@Lp[[this.group]]$ov.idx[[this.level]]
        COV <- COV[ov.idx, ov.idx, drop = FALSE]
        YBAR <- YBAR[ov.idx]
      } else if (this.level == 2L) {
        if (local.twolevel.method == "h1") {
          COV <- H1$implied$cov[[2]]
          YBAR <- H1$implied$mean[[2]]
        } else if (local.twolevel.method == "anova") {
          COV <- FIT@SampleStats@YLp[[this.group]][[2]]$Sigma.B
          YBAR <- FIT@SampleStats@YLp[[this.group]][[2]]$Mu.B
        } else if (local.twolevel.method == "mean") {
          S.PW <- FIT@SampleStats@YLp[[this.group]][[2]]$Sigma.W
          NJ <- FIT@SampleStats@YLp[[this.group]][[2]]$s
          Y2 <- FIT@SampleStats@YLp[[this.group]][[2]]$Y2
          # grand mean
          MU.Y <- (FIT@SampleStats@YLp[[this.group]][[2]]$Mu.W + FIT@SampleStats@YLp[[this.group]][[2]]$Mu.B)
          Y2c <- t(t(Y2) - MU.Y) # MUST be centered
          YB <- crossprod(Y2c) / nrow(Y2c)
          COV <- YB - 1 / NJ * S.PW
          YBAR <- FIT@SampleStats@YLp[[this.group]][[2]]$Mu.B
        }

        # reduce
        ov.idx <- FIT@Data@Lp[[this.group]]$ov.idx[[this.level]]
        COV <- COV[ov.idx, ov.idx, drop = FALSE]
        YBAR <- YBAR[ov.idx]
      } else {
        lav_msg_stop(gettext("level 3 not supported (yet)."))
      }

      # single level
    } else {
      this.group <- b
      YBAR <- FIT@h1$implied$mean[[b]] # EM version if missing="ml"
      COV <- FIT@h1$implied$cov[[b]]
      # rescale COV?
      if (FIT@Model@categorical || FIT@Model@correlation) {
        SCALE.vector <- 1 / (drop(DELTA[[b]]))
        COV <- SCALE.vector * COV * rep(SCALE.vector, each = ncol(COV))
        YBAR <- SCALE.vector * YBAR # Checkme!
      }
      # do we need ICOV?
      if (local.M.method == "GLS") {
        if (FIT@Options$sample.cov.rescale) {
          # get unbiased S
          N <- FIT@SampleStats@nobs[[b]]
          COV.unbiased <- COV * N / (N - 1)
          ICOV <- solve(COV.unbiased)
        } else {
          ICOV <- solve(COV)
        }
      }
    }

    # compute mapping matrix 'M'
    this.lambda <- LAMBDA[[b]]
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      this.lambda <- this.lambda[, -lavpta$vidx$lv.interaction[[b]]]
    }
    Mb <- lav_sam_mapping_matrix(
      LAMBDA = this.lambda,
      THETA = THETA[[b]],
      S = COV, S.inv = ICOV,
      method = local.M.method,
      warn = lavoptions$warn
    )
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      tmp <- Mb
      Mb <- matrix(0, nrow = ncol(LAMBDA[[b]]), ncol = nrow(LAMBDA[[b]]))
      Mb[-lavpta$vidx$lv.interaction[[b]], ] <- tmp
    }

    # handle observed-only variables
    dummy.ov.idx <- c(
      FIT@Model@ov.x.dummy.ov.idx[[b]],
      FIT@Model@ov.y.dummy.ov.idx[[b]]
    )
    dummy.lv.idx <- c(
      FIT@Model@ov.x.dummy.lv.idx[[b]],
      FIT@Model@ov.y.dummy.lv.idx[[b]]
    )
    if (length(dummy.ov.idx)) {
      Mb[dummy.lv.idx, ] <- 0
      Mb[cbind(dummy.lv.idx, dummy.ov.idx)] <- 1
    }

    # here, we remove the lv.interaction row(s) from Mb
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      Mb <- Mb[-lavpta$vidx$lv.interaction[[b]], ]
    }

    # compute EETA
    if (lavoptions$meanstructure) {
      EETA[[b]] <- lav_sam_eeta(M = Mb, YBAR = YBAR, NU = NU[[b]])
    }

    # compute VETA
    if (sam.method == "local") {
      tmp <- lav_sam_veta(
        M = Mb, S = COV, THETA = THETA[[b]],
        alpha.correction = local.options[["alpha.correction"]],
        lambda.correction = local.options[["lambda.correction"]],
        N <- FIT@SampleStats@nobs[[this.group]],
        dummy.lv.idx = dummy.lv.idx,
        extra = TRUE
      )
      VETA[[b]] <- tmp[, , drop = FALSE] # drop attributes
      alpha[[b]] <- attr(tmp, "alpha")
      lambda[[b]] <- attr(tmp, "lambda")
    } else {
      # FSR -- no correction
      VETA[[b]] <- Mb %*% COV %*% t(Mb)
    }

    # standardize? not really needed, but we may have 1.0000001
    # as variances, and this may lead to false convergence
    if (FIT@Options$std.lv) {
      VETA[[b]] <- stats::cov2cor(VETA[[b]])
    }

    # lv.names, including dummy-lv covariates
    psi.idx <- which(names(FIT@Model@GLIST) == "psi")[b]
    lv.names <- FIT@Model@dimNames[[psi.idx]][[b]] # including dummy/interact/.
    if (!lv.interaction.flag) {
      dimnames(VETA[[b]]) <- FIT@Model@dimNames[[psi.idx]]
    } else {
      lv.int.names <- lavpta$vnames$lv.interaction[[b]]
      # remove interaction terms
      lv.names1 <- lv.names[!lv.names %in% lv.int.names]
      colnames(VETA[[b]]) <- rownames(VETA[[b]]) <- lv.names1
    }

    # compute model-based RELiability
    MSM <- Mb %*% COV %*% t(Mb)
    # REL[[b]] <- diag(VETA[[b]]] %*% solve(MSM)) # CHECKme!
    REL[[b]] <- diag(VETA[[b]]) / diag(MSM)

    # check for lv.interactions
    if (lv.interaction.flag && length(lv.int.names) > 0L) {
      if (FIT@Model@categorical || FIT@Model@correlation) {
        lav_msg_stop(gettext("SAM + lv interactions do not work (yet) if
                             correlation structures are used."))
      }

      # EETA2
      EETA1 <- EETA[[b]]
      EETA[[b]] <- lav_sam_eeta2(
        EETA = EETA1, VETA = VETA[[b]],
        lv.names = lv.names1,
        lv.int.names = lv.int.names
      )

      # VETA2
      if (sam.method == "local") {
        # reorder FS[[b]] if needed
        FS.b <- FS[[b]][, lv.names1, drop = FALSE]
        tmp <- lav_sam_veta2(
          FS = FS.b, M = Mb,
          VETA = VETA[[b]], EETA = EETA1,
          THETA = THETA[[b]],
          lv.names = lv.names1,
          lv.int.names = lv.int.names,
          dummy.lv.names = lv.names[dummy.lv.idx],
          alpha.correction = local.options[["alpha.correction"]],
          lambda.correction = local.options[["lambda.correction"]],
          extra = TRUE
        )
        VETA[[b]] <- tmp[, , drop = FALSE] # drop attributes
        alpha[[b]] <- attr(tmp, "alpha")
        lambda[[b]] <- attr(tmp, "lambda")
      } else {
        lav_msg_fixme("not ready yet!")
        # FSR -- no correction
        VETA[[b]] <- lav_sam_fs2(
          FS = FS[[b]],
          lv.names = lv.names1, lv.int.names = lv.int.names
        )
      }
    }

    # store Mapping matrix for this block
    M[[b]] <- Mb
  } # blocks

  # label blocks
  if (nblocks > 1L) {
    names(EETA) <- FIT@Data@block.label
    names(VETA) <- FIT@Data@block.label
    names(REL) <- FIT@Data@block.label
  }

  # store EETA/VETA/M/alpha/lambda
  STEP1$VETA <- VETA
  STEP1$EETA <- EETA
  STEP1$REL <- REL
  STEP1$M <- M
  STEP1$lambda <- lambda
  STEP1$alpha <- alpha

  if (lavoptions$verbose) {
    cat("done.\n")
  }

  STEP1
}
