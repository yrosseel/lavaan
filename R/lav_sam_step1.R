# step 1 in SAM: fitting the measurement blocks

lav_sam_step1 <- function(cmd = "sem", mm.list = NULL, mm.args = list(),
                          FIT = FIT, sam.method = "local") {
  lavoptions <- FIT@Options
  lavpta <- FIT@pta
  PT <- FIT@ParTable
  nblocks <- lavpta$nblocks
  ngroups <- lavpta$ngroups

  if (lav_verbose()) {
    cat("Fitting the measurement part:\n")
  }

  # local only -> handle missing data
  # if (sam.method %in% c("local", "fsr")) {
  #   # if missing = "listwise", make data complete, to avoid different
  #   # datasets per measurement block
  #   if (lavoptions$missing == "listwise") {
  #     # FIXME: make this work for multiple groups!!
  #     OV <- unique(unlist(FIT@pta$vnames$ov))
  #     # add group/cluster/sample.weights variables (if any)
  #     OV <- c(
  #       OV, FIT@Data@group, FIT@Data@cluster,
  #       FIT@Data@sampling.weights
  #     )
  #     data <- na.omit(data[, OV])
  #   }
  # }

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

  # do we have at least 1 'regular' (measured) latent variable?
  LV.names <- unique(unlist(FIT@pta$vnames$lv.regular))

  eqs.x <- unlist(FIT@pta$vnames$eqs.x)
  eqs.y <- unlist(FIT@pta$vnames$eqs.y)

  #if (length(eqs.x) == 0L && length(eqs.y) == 0L) {
  #  lav_msg_warn(gettext("the model does seem to contain a structural part
  #                        (i.e., regressions); consider using sem() instead"))
  #} else
  if (length(LV.names) == 0L) {
    lav_msg_warn(gettext("structural model does not contain any (measured)
                          latent variables; consider using sem() instead"))
  }

  # check for higher-order factors
  # 0.6-20: now we do!
  LV.IND.names <- unique(unlist(FIT@pta$vnames$lv.ind))
  lv.higherorder.flag <- FALSE
  if (length(LV.IND.names) > 0L) {
    lv.higherorder.flag <- TRUE
    LV.names <- LV.names[!LV.names %in% LV.IND.names]
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
  if (lavoptions$se %in% c("none", "bootstrap")) {
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
  lavoptions.mm$check.post <- FALSE # neg lv variances may be overriden
  lavoptions.mm$check.gradient <- FALSE # too sensitive in large model (global)
  lavoptions.mm$baseline <- FALSE
  lavoptions.mm$bounds <- "wide.zerovar"

  # ALWAYS conditional.x = FALSE!
  # even if global model uses conditional.x = TRUE
  # this should not affect the measurement models (if the covariates act on
  # the structural part only)
  lavoptions.mm$conditional.x = FALSE

  # override with user-specified mm.args
  lavoptions.mm <- modifyList(lavoptions.mm, mm.args)

  # create MM slotOptions
  slotOptions.mm <- lav_options_set(lavoptions.mm)

  # we assume the same number/names of lv's per group!!!
  MM.FIT <- vector("list", nMMblocks) # fitted object

  # for joint model later
  if (!lavoptions$se %in% c("none", "bootstrap")) {
    Sigma.11 <- matrix(0, npar, npar)
    colnames(Sigma.11) <- rownames(Sigma.11) <-
      lav_partable_labels(FIT@ParTable, type = "free")
  }
  step1.free.idx <- integer(0L)
  block.mm.idx  <- vector("list", length = nMMblocks)
  block.ptm.idx <- vector("list", length = nMMblocks)

  # NOTE: we should explicitly add zero-constrained LV covariances
  # to PT, and keep them zero in PTM
  if (cmd == "lavaan") {
    add.lv.cov <- FALSE
  } else {
    add.lv.cov <- TRUE
  }

  # fit mm model for each measurement block
  for (mm in seq_len(nMMblocks)) {
    if (lav_verbose()) {
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
      lv.names = mm.list[[mm]],
    )
    mm.idx <- attr(PTM, "idx")
    attr(PTM, "idx") <- NULL
    PTM$est <- NULL
    PTM$se <- NULL
    block.mm.idx[[mm]] <- mm.idx

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
    # get rid of ov.names.x
    if (!slotOptions.mm$conditional.x) {
      slotData.block@ov.names.x <-
        lapply(seq_len(nblocks), function(x) character(0L))
      slotData.block@eXo <-
        lapply(seq_len(nblocks), function(x) NULL)
    }

    # if data.type == "moment", (re)create sample.cov and sample.nobs
    if (FIT@Data@data.type == "moment") {
      if (ngroups == 1L) {
        mm.sample.cov <- lavInspect(FIT, "h1")$cov
        mm.sample.mean <- NULL
        if (FIT@Model@meanstructure) {
          mm.sample.mean <- lavInspect(FIT, "h1")$mean
        }
        mm.sample.nobs <- FIT@SampleStats@nobs[[1L]]
      } else {
        cov.list <- lapply(lavTech(FIT, "h1", add.labels = TRUE),
                           "[[", "cov")
        mm.sample.cov <- lapply(seq_len(ngroups),
          function(x) cov.list[[x]][ov.names.block[[x]], ov.names.block[[x]]])
        mm.sample.mean <- NULL
        if (FIT@Model@meanstructure) {
          mean.list <- lapply(lavTech(FIT, "h1", add.labels = TRUE),
                           "[[", "mean")
          mm.sample.mean <- lapply(seq_len(ngroups),
            function(x) mean.list[[x]][ov.names.block[[x]]])
        }
        mm.sample.nobs <- FIT@SampleStats@nobs
      }
    }


    # handle single block 1-factor CFA with (only) two indicators
    if (length(unlist(ov.names.block)) == 2L && ngroups == 1L) {
      lambda.idx <- which(PTM$op == "=~")
	  # check if both factor loadings are fixed
	  # (note: this assumes std.lv = FALSE)
	  if(any(PTM$free[lambda.idx] != 0)) {
        PTM$free[  lambda.idx] <- 0L
        PTM$ustart[lambda.idx] <- 1
        PTM$start[ lambda.idx] <- 1
        free.idx <- which(as.logical(PTM$free))
		# adjust free counter
        if (length(free.idx) > 0L) {
          PTM$free[free.idx] <- seq_len(length(free.idx))
        }
		# warn about it (needed?)
        lav_msg_warn(gettextf(
          "measurement block [%1$s] (%2$s) contains only two indicators;
          -> fixing both factor loadings to unity",
          mm, lav_msg_view(mm.list[[mm]], "none")))
      }
    }

    # fit this measurement model only
    # (question: can we re-use even more slots?)
    if (FIT@Data@data.type == "full") {
      fit.mm.block <- lavaan(
        model = PTM, slotData = slotData.block,
        slotOptions = slotOptions.mm, debug = FALSE, verbose = FALSE
      )
    } else if (FIT@Data@data.type == "moment") {
      slotOptions.mm$sample.cov.rescale <- FALSE
      fit.mm.block <- lavaan(
        model = PTM, slotData = slotData.block,
        sample.cov = mm.sample.cov, sample.mean = mm.sample.mean,
        sample.nobs = mm.sample.nobs,
        slotOptions = slotOptions.mm, debug = FALSE, verbose = FALSE
      )
    }

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
    block.ptm.idx[[mm]] <- ptm.idx
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
    if (!lavoptions$se %in% c("none", "bootstrap")) {
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
      this.mm.idx <- sort.int(par.idx)
      step1.free.idx <- c(step1.free.idx, this.mm.idx) # all combined
    }
  } # measurement block

  # only keep 'measurement part' parameters in Sigma.11
  if (!lavoptions$se %in% c("none", "bootstrap")) {
    Sigma.11 <- Sigma.11[step1.free.idx, step1.free.idx, drop = FALSE]
  } else {
    Sigma.11 <- NULL
  }

  # create STEP1 list
  STEP1 <- list(
    MM.FIT = MM.FIT, Sigma.11 = Sigma.11,
    step1.free.idx = step1.free.idx,
    block.mm.idx = block.mm.idx,
    block.ptm.idx = block.ptm.idx,
    PT.free = PT.free,
    mm.list = mm.list, PT = PT
  )

  STEP1
}

## STEP 1b: compute Var(eta) and E(eta) per block
##          only needed for local/fsr approach!
lav_sam_step1_local <- function(STEP1 = NULL, FIT = NULL, Y = NULL,
                                sam.method = "local",
                                local.options = list(
                                  M.method = "ML",
                                  lambda.correction = TRUE,
                                  alpha.correction = 0L,
                                  twolevel.method = "h1"
                                ),
                                return.cov.iveta2 = TRUE,
                                return.FS = FALSE) {
  # local.M.method
  local.M.method <- toupper(local.options[["M.method"]])
  if (!local.M.method %in% c("GLS", "ML", "ULS")) {
    lav_msg_stop(gettext(
      "local option M.method should be one of ML, GLS or ULS."))
  }

  lavoptions <- FIT@Options
  lavpta <- FIT@pta
  nblocks <- lavpta$nblocks

  # flags
  lv.interaction.flag <- FALSE
  lv.higherorder.flag <- FALSE
  if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
    lv.interaction.flag <- TRUE
  }
  if (length(unlist(lavpta$vnames$lv.ind)) > 0L) {
    lv.higherorder.flag <- TRUE
  }

  if (lav_verbose()) {
    cat("Constructing the mapping matrix using the ",
      local.M.method, " method ... ",
      sep = ""
    )
  }

  # all the measurement parameters are already stored in PT
  PT <- STEP1$PT
  if (FIT@Model@ceq.simple.only) {
    x.free <- PT$est[PT$free > 0 & !duplicated(PT$free)]
  } else {
    x.free <- PT$est[PT$free > 0]
  }
  # check for NA values (eg in BETA); set them to zero
  x.free[!is.finite(x.free)] <- 0

  lavmodel.tmp <- lav_model_set_parameters(FIT@Model, x = x.free)
  LAMBDA <- THETA <- BETA <- PSI <- NU <- DELTA <- NULL

  # create LAMBDA
  lambda.idx <- which(names(FIT@Model@GLIST) == "lambda")
  LAMBDA <- lavmodel.tmp@GLIST[lambda.idx]

  # create THETA
  theta.idx <- which(names(FIT@Model@GLIST) == "theta")
  THETA <- lavmodel.tmp@GLIST[theta.idx]

  # NU
  if (FIT@Model@meanstructure) {
    nu.idx <- which(names(FIT@Model@GLIST) == "nu")
    NU <- lavmodel.tmp@GLIST[nu.idx]
  }

  # DELTA
  if (FIT@Model@categorical || FIT@Model@correlation) {
    delta.idx <- which(names(FIT@Model@GLIST) == "delta")
    DELTA <- lavmodel.tmp@GLIST[delta.idx]
  }

  # BETA/PSI
  if (lv.higherorder.flag) {
    beta.idx <- which(names(FIT@Model@GLIST) == "beta")
    BETA <- lavmodel.tmp@GLIST[beta.idx]
    psi.idx <- which(names(FIT@Model@GLIST) == "psi")
    PSI <- lavmodel.tmp@GLIST[psi.idx]
  }

  # GAMMA (only for names, if conditional.x)
  #if (FIT@Model@conditional.x) {
  #  gamma.idx <- which(names(FIT@Model@GLIST) == "gamma")
  #}

  # handle dummy's + higher-order + rank-deficient
  for (b in seq_len(nblocks)) {
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
      if (length(dummy.ov.idx)) {
        NU[[b]][dummy.ov.idx, 1] <- 0
      }
    }

    # get ALL lv names (including dummy ov.x/ov.y)
    lv.names <- FIT@Model@dimNames[[lambda.idx[b]]][[2L]]

    # if conditional.x, we must add the ov.names.x manually
    # if (FIT@Model@conditional.x) {
    #   exo.names <- FIT@Model@dimNames[[gamma.idx[b]]][[2L]]
    #   lv.names <- c(lv.names, exo.names)
    # }

    # handle higher-order factors here
    if (length(lavpta$vidx$lv.ind[[b]]) > 0L) {
      lv.ind.names <- lavpta$vnames$lv.ind[[b]]
      lv.target <- lv.names[!lv.names %in% lv.ind.names]

      target.idx <- match(lv.target, lv.names)
      other.idx <- seq_len(length(lv.names))[-target.idx]

      this.beta <- BETA[[b]]
      this.beta[is.na(this.beta)] <- 0

      IB <- diag(nrow(this.beta)) - this.beta
      IB.inv <- solve(IB)
      LB.inv <- LAMBDA[[b]] %*% IB.inv

      # replace LAMBDA
      LAMBDA[[b]] <- LB.inv[,target.idx,drop = FALSE]

      PSI.other <- PSI[[b]][other.idx, other.idx, drop = FALSE]
      LB.inv2 <- LB.inv[, other.idx, drop = FALSE]

      # replace THETA
      THETA[[b]] <- LB.inv2 %*% PSI.other %*% t(LB.inv2) + THETA[[b]]
    }

    # check if LAMBDA has full column rank
    this.lambda <- LAMBDA[[b]]
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      if (length(lavpta$vidx$lv.ind[[b]]) > 0L) {
        rm.idx <- c(match(lavpta$vnames$lv.ind[[b]], lv.names),
                    match(lavpta$vnames$lv.interaction[[b]], lv.names))
        this.lambda <- this.lambda[, -rm.idx, drop = FALSE]
      } else {
        rm.idx <- match(lavpta$vnames$lv.interaction[[b]], lv.names)
        this.lambda <- this.lambda[, -rm.idx, drop = FALSE]
      }
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

  VETA     <- vector("list", nblocks)
  MSM..    <- vector("list", nblocks)
  MTM..    <- vector("list", nblocks)
  FS.mean  <- vector("list", nblocks)
  #FS.gamma <- vector("list", nblocks)
  FS       <- vector("list", nblocks)
  COV.IVETA2 <- vector("list", nblocks)
  REL <- vector("list", nblocks)
  alpha <- vector("list", nblocks)
  lambda <- vector("list", nblocks)
  if (lavoptions$meanstructure) {
    EETA <- vector("list", nblocks)
  } else {
    EETA <- NULL
  }
  M <- vector("list", nblocks)
  LV.NAMES <- vector("list", nblocks)

  # if (lv.interaction.flag && is.null(FS)) {
  #   # compute Bartlett factor scores
  #   FS <- vector("list", nblocks)
  #   # FS.mm <- lapply(STEP1$MM.FIT, lav_predict_eta_bartlett)
  #   FS.mm <- lapply(STEP1$MM.FIT, lavPredict,
  #     method = "Bartlett",
  #     drop.list.single.group = FALSE
  #   )
  #   for (b in seq_len(nblocks)) {
  #     tmp <- lapply(
  #       1:length(STEP1$MM.FIT),
  #       function(x) FS.mm[[x]][[b]]
  #     )
  #     LABEL <- unlist(lapply(tmp, colnames))
  #     FS[[b]] <- do.call("cbind", tmp)
  #     colnames(FS[[b]]) <- LABEL
  #     FS[[b]] <- FIT@Data@X[[b]] %*%

  #     # dummy lv's? (both 'x' and 'y'!)
  #     dummy.ov.idx <- c(FIT@Model@ov.y.dummy.ov.idx[[b]],
  #                       FIT@Model@ov.x.dummy.ov.idx[[b]])
  #     dummy.lv.idx <- c(FIT@Model@ov.y.dummy.lv.idx[[b]],
  #                       FIT@Model@ov.x.dummy.lv.idx[[b]])
  #     if (length(dummy.lv.idx) > 0L) {
  #       FS.obs <- FIT@Data@X[[b]][, dummy.ov.idx, drop = FALSE]
  #       colnames(FS.obs) <- FIT@Data@ov.names[[b]][dummy.ov.idx]
  #       FS[[b]] <- cbind(FS[[b]], FS.obs)
  #     }
  #   }
  # }

  # compute VETA/EETA per block
  for (b in seq_len(nblocks)) {

    # which group is this?
    this.group <- floor(b / FIT@Data@nlevels + 0.5)

    # lv.names, including dummy-lv covariates
    psi.idx <- which(names(FIT@Model@GLIST) == "psi")[b]
    lv.names.b <- FIT@Model@dimNames[[psi.idx]][[1L]] # including dummy/inter.
    # if (FIT@Model@conditional.x) {
    #   exo.names <- FIT@Model@dimNames[[gamma.idx[b]]][[2L]]
    #   lv.names.b <- c(lv.names.b, exo.names)
    # }
    rm.idx <- integer(0L)

    # higher-order? remove lower-order factors
    if (lv.higherorder.flag && length(lavpta$vnames$lv.ind[[b]]) > 0L) {
      rm.idx <- c(rm.idx, match(lavpta$vnames$lv.ind[[b]], lv.names.b))
    }

    # interaction terms? remove them for VETA
    if (lv.interaction.flag && length(lavpta$vnames$lv.interaction[[b]]) > 0L) {
      rm.idx <- c(rm.idx, match(lavpta$vnames$lv.interaction[[b]], lv.names.b))
      lv.int.names <- lavpta$vnames$lv.interaction[[b]]
    }

    # final names for EETA/VETA (not including interaction terms!)
    lv.names1 <- lv.names.b
    if (length(rm.idx) > 0L) {
      lv.names1 <- lv.names.b[-rm.idx]
    }
    LV.NAMES[[b]] <- lv.names1

    # get sample statistics for this block
    COV <- STEP1$COV[[b]]
    YBAR <- drop(STEP1$YBAR[[b]])

    # rescale COV?
    if (FIT@Data@nlevels == 1L &&
        (FIT@Model@categorical || FIT@Model@correlation)) {
      SCALE.vector <- 1 / (drop(DELTA[[b]]))
      COV <- SCALE.vector * COV * rep(SCALE.vector, each = ncol(COV))
      YBAR <- SCALE.vector * YBAR # Checkme!
    }

    # do we need ICOV?
    if (local.M.method == "GLS") {
      if (FIT@Options$sample.cov.rescale) {
        # get unbiased S
        N <- FIT@SampleStats@nobs[[this.group]]
        COV.unbiased <- COV * N / (N - 1)
        ICOV <- solve(COV.unbiased)
      } else {
        ICOV <- solve(COV)
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
      method = local.M.method
    )
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      tmp <- Mb
      Mb <- matrix(0, nrow = ncol(LAMBDA[[b]]), ncol = nrow(LAMBDA[[b]]))
      Mb[-lavpta$vidx$lv.interaction[[b]], ] <- tmp
    }

    # handle observed-only variables (needed?)
    dummy.ov.idx <- c(
      FIT@Model@ov.x.dummy.ov.idx[[b]],
      FIT@Model@ov.y.dummy.ov.idx[[b]]
    )
    dummy.lv.idx <- c(
      FIT@Model@ov.x.dummy.lv.idx[[b]],
      FIT@Model@ov.y.dummy.lv.idx[[b]]
    )

    # handle conditional.x
    # if (FIT@Model@conditional.x) {
    #   I0 <- diag(x = 0, nrow = length(exo.names))
    #   I1 <- diag(x = 1, nrow = length(exo.names))
    #   Mb <- lav_matrix_bdiag(Mb, I1)
    #   LAMBDA[[b]] <- lav_matrix_bdiag(LAMBDA[[b]], I1)
    #   THETA[[b]] <- lav_matrix_bdiag(THETA[[b]], I0)
    #   NU[[b]] <- c(drop(NU[[b]]), numeric(length(exo.names)))
    # }

    # fix dummy.lv.idx if we have higher-order factors!
    if (lv.higherorder.flag) {
      dummy.lv.idx <- match(lv.names.b[dummy.lv.idx], lv.names1)
    }

    if (length(dummy.ov.idx)) {
      Mb[dummy.lv.idx, ] <- 0
      Mb[cbind(dummy.lv.idx, dummy.ov.idx)] <- 1
    }

    # here, we remove the lv.interaction row(s) from Mb
    # FIXME: if we have higher order factors!
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      Mb <- Mb[-lavpta$vidx$lv.interaction[[b]], ]
    }

    # compute EETA
    if (lavoptions$meanstructure) {
      EETA[[b]] <- lav_sam_eeta(M = Mb, YBAR = YBAR, NU = NU[[b]])
	  FS.mean[[b]] <- EETA[[b]] # ok if no interaction
    }

    # compute VETA
    if (sam.method == "local") {
      this.group <- floor(b / FIT@Data@nlevels + 0.5)
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
      lambda[[b]] <- attr(tmp, "lambda.star")
	  MSM..[[b]] <- attr(tmp, "MSM")
	  MTM..[[b]] <- attr(tmp, "MTM")
    } else if (sam.method == "cfsr") {
      # first, we need to 'true' VETA (to get Sigma)
      this.group <- floor(b / FIT@Data@nlevels + 0.5)
      tmp <- lav_sam_veta(
        M = Mb, S = COV, THETA = THETA[[b]],
        alpha.correction = 0L,
        lambda.correction = local.options[["lambda.correction"]],
        N <- FIT@SampleStats@nobs[[this.group]],
        dummy.lv.idx = dummy.lv.idx,
        extra = FALSE
      )
      VETA[[b]] <- tmp[, , drop = FALSE]
      # compute 'Sigma'
      Sigma <- this.lambda %*% VETA[[b]] %*% t(this.lambda) + THETA[[b]]
      tmat <- lav_predict_tmat_det_internal(Sigma = Sigma, Veta = VETA[[b]],
                                            Lambda = this.lambda)
      A <- tmat %*% Mb
      VETA[[b]] <- A %*% COV %*% t(A)
    } else {
      # FSR -- no correction
      VETA[[b]] <- Mb %*% COV %*% t(Mb)
    }

    # standardize? not really needed, but we may have 1.0000001
    # as variances, and this may lead to false convergence
    if (FIT@Options$std.lv) {
      # warning: we should only do this for the LVs, not the
      # observed variables
      if (length(dummy.lv.idx) == 0L) {
        VETA[[b]] <- stats::cov2cor(VETA[[b]])
      } else {
        tmp <- VETA[[b]]
        tmp.lv <- stats::cov2cor(VETA[[b]][-dummy.lv.idx,
                                           -dummy.lv.idx, drop = FALSE])
        VETA[[b]][-dummy.lv.idx, -dummy.lv.idx] <- tmp.lv
      }
    }
    colnames(VETA[[b]]) <- rownames(VETA[[b]]) <- lv.names1

    # compute model-based RELiability
    # REL[[b]] <- diag(VETA[[b]]] %*% solve(MSM..[[b]])) # CHECKme! -> done, must be:
    REL[[b]] <- diag(VETA[[b]]) / diag(MSM..[[b]]) #!

    # check for lv.interactions
    if (lv.interaction.flag && length(lv.int.names) > 0L) {
      if (FIT@Model@categorical || FIT@Model@correlation) {
        lav_msg_stop(gettext("SAM + lv interactions do not work (yet) if
                             correlation structures are used."))
      }

      # compute Bartlett factor scores here
      if (is.null(Y)) {
        Yb <- FIT@Data@X[[b]]
      } else {
        Yb <- Y[[b]]
      }
      # center
      Yb.c <- t( t(Yb) - drop(NU[[b]]) )
      FS.b <- Yb.c %*% t(Mb)
      colnames(FS.b) <- lv.names1
      # FIXME: what about observed covariates?

      # EETA2
      EETA1 <- EETA[[b]]
      EETA[[b]] <- lav_sam_eeta2(
        EETA = EETA1, VETA = VETA[[b]],
        lv.names = lv.names1,
        lv.int.names = lv.int.names
      )

      # VETA2
      if (sam.method == "local") {
        tmp <- lav_sam_veta2(
          FS = FS.b, M = Mb,
          VETA = VETA[[b]], EETA = EETA1,
          THETA = THETA[[b]],
          lv.names = lv.names1,
          lv.int.names = lv.int.names,
          dummy.lv.names = lv.names.b[dummy.lv.idx],
          alpha.correction = local.options[["alpha.correction"]],
          lambda.correction = local.options[["lambda.correction"]],
          return.FS = return.FS,
          return.cov.iveta2 = return.cov.iveta2,
          extra = TRUE
        )
        VETA[[b]] <- tmp[, , drop = FALSE] # drop attributes
        alpha[[b]] <- attr(tmp, "alpha")
        lambda[[b]] <- attr(tmp, "lambda.star")
        MSM..[[b]] <- attr(tmp, "MSM")
        MTM..[[b]] <- attr(tmp, "MTM")
		FS.mean[[b]] <- attr(tmp, "FS.mean")
        if (return.FS) {
          FS[[b]] <- attr(tmp, "FS")
        }
        if (return.cov.iveta2) {
          COV.IVETA2[[b]] <- attr(tmp, "cov.iveta2")
        }
        #FS.gamma[[b]] <- attr(tmp, "FS.gamma")
      } else {
        lav_msg_fixme("not ready yet!")
        # FSR -- no correction
        VETA[[b]] <- lav_sam_fs2(
          FS = FS.b,
          lv.names = lv.names1, lv.int.names = lv.int.names
        )
      }
    }

    # store Mapping matrix for this block
    M[[b]] <- Mb
  } # blocks

  # label blocks
  if (nblocks > 1L) {
    names(EETA)       <- FIT@Data@block.label
    names(VETA)       <- FIT@Data@block.label
    names(REL)        <- FIT@Data@block.label
	names(MSM..)      <- FIT@Data@block.label
	names(MTM..)      <- FIT@Data@block.label
	names(FS.mean)    <- FIT@Data@block.label
    names(FS)         <- FIT@Data@block.label
    names(COV.IVETA2) <- FIT@Data@block.label
    #names(FS.gamma) <- FIT@Data@block.label
  }

  # handle conditional.x: add res.slopes, cov.x and mean.x
  if (FIT@Model@conditional.x) {
    res.slopes <- vector("list", length = nblocks)
    for (b in seq_len(nblocks)) {
      res.slopes[[b]] <- M[[b]] %*% FIT@h1$implied$res.slopes[[b]]
    }
    attr(VETA, "res.slopes") <- res.slopes
    attr(VETA, "cov.x") <- FIT@h1$implied$cov.x
    attr(VETA, "mean.x") <- FIT@h1$implied$mean.x
  }

  # store EETA/VETA/M/alpha/lambda
  STEP1$VETA     <- VETA
  STEP1$EETA     <- EETA
  STEP1$REL      <- REL
  STEP1$M        <- M
  STEP1$lambda   <- lambda
  STEP1$alpha    <- alpha
  STEP1$MSM      <- MSM..
  STEP1$MTM      <- MTM..
  STEP1$FS.mean  <- FS.mean
  STEP1$FS       <- FS
  STEP1$COV.IVETA2 <- COV.IVETA2
  #STEP1$FS.gamma <- FS.gamma
  STEP1$LV.NAMES <- LV.NAMES
  # store also sam.method and local.options
  STEP1$sam.method <- sam.method
  STEP1$local.options <- local.options

  if (lav_verbose()) {
    cat("done.\n")
  }

  STEP1
}


lav_sam_step1_local_jac <- function(STEP1 = NULL, FIT = NULL) {

  lavdata <- FIT@Data
  lavsamplestats <- FIT@SampleStats
  lavmodel <- FIT@Model
  lavpta <- FIT@pta
  nblocks <- lavpta$nblocks

  local.options <- STEP1$local.options
  sam.method <- STEP1$sam.method

  ngroups <- lavdata@ngroups
  if (ngroups > 1L) {
    stop("IJ local SEs: not available with multiple groups!\n")
    # if multiple groups:
    # - we have a separate Gamma, h1.expected, delta matrix per group
    # - but we have only 1 observed information matrix, reflecting possible
    #   across-group equality constraints
    # - we may need to use the same procedure as for robust test statistics:
    #   create a (huge) block-diagonal Gamma matrix, and one big 'JAC'
    #   matrix...
  }
  g <- 1L
  if (lavmodel@categorical) {
    stop("IJ local SEs: not available for the categorical setting (yet)!\n")
  }
  nMMblocks <- length(STEP1$MM.FIT)

  # JAC = (JACc %*% JACa) + JACb
  # - rows are the elements of vech(VETA)
  # - cols are the elements of vech(S)

  # JACa: mm.theta x vech(S)
  # JACc: vech(VETA) x mm.theta (keeping S fixed)
  # JACb: vech(VETA) x vech(S)  (keeping mm.theta fixed)

  # JACa: jacobian of theta.mm = f(vech(S))
  JACa <- matrix(0, nrow = length(FIT@ParTable$lhs), # we select later
                    ncol = length(FIT@SampleStats@WLS.obs[[g]]))
  for (mm in seq_len(nMMblocks)) {
    fit.mm.block <- STEP1$MM.FIT[[mm]]
    mm.h1.expected   <- lavTech(fit.mm.block, "h1.information.expected")
    mm.delta         <- lavTech(fit.mm.block, "Delta")
    mm.inv.observed  <- lavTech(fit.mm.block, "inverted.information.observed")

    #h1.info <- matrix(0, nrow(mm.h1.expected[[1]]), ncol(mm.h1.expected[[1]]))
    #for (g in seq_len(ngroups)) {
    #  fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    #  tmp <- fg * (mm.h1.expected[[g]] %*% mm.delta[[g]])
    #  h1.info <- h1.info + tmp
    #}
    mm.jac <- t(mm.h1.expected[[g]] %*% mm.delta[[g]] %*% mm.inv.observed)
    # keep only rows that are also in FIT@ParTable
    mm.keep.idx <- fit.mm.block@ParTable$free[STEP1$block.ptm.idx[[mm]]]
    mm.jac <- mm.jac[mm.keep.idx, , drop = FALSE]

    # select 'S' elements (row index)
    mm.ov.idx <- match(STEP1$MM.FIT[[mm]]@Data@ov.names[[g]],
                       lavdata@ov.names[[g]])
    mm.nvar <- length(lavdata@ov.names[[g]])
    mm.col.idx <- lav_matrix_vech_which_idx(mm.nvar, idx = mm.ov.idx,
      add.idx.at.start = lavmodel@meanstructure)
    mm.row.idx <- STEP1$block.mm.idx[[mm]][STEP1$block.ptm.idx[[mm]]]
    JACa[mm.row.idx, mm.col.idx] <- mm.jac
  }

  # keep only 'LAMBDA/THETA' parameters
  PT <- STEP1$PT
  # only ov.names that are actually used in the measurement models
  ov.names <- unique(unlist(lapply(STEP1$MM.FIT, lavNames, "ov")))
  lambda.idx <- which(PT$op == "=~" & PT$free > 0L & !duplicated(PT$free))
  theta.idx  <- which(PT$op == "~~" & PT$free > 0L & !duplicated(PT$free) &
                      PT$lhs %in% ov.names & PT$rhs %in% ov.names)
  nu.idx <- integer(0L)
  if (lavmodel@meanstructure) {
    nu.idx     <- which(PT$op == "~1" & PT$free > 0L & !duplicated(PT$free) &
                        PT$lhs %in% ov.names)
  }
  delta.idx <- integer(0L)
  if (lavmodel@categorical || lavmodel@correlation) {
    delta.idx <- which(PT$op == "~*~" & PT$free > 0L & !duplicated(PT$free))
  }
  beta.idx <- psi.idx <- integer(0L)
  lv.ind <- unlist(lavpta$vnames$lv.ind)
  if (length(lv.ind) > 0L) {
    beta.idx <- which(PT$op == "=~" & PT$free > 0L & !duplicated(PT$free) &
                      PT$rhs %in% lv.ind)
    psi.idx <- which(PT$op == "~~" & PT$free > 0L & !duplicated(PT$free) &
                     PT$rhs %in% lv.ind & PT$lhs %in% lv.ind)
  }
  # keep only these free parameters (measurement only)
  keep.idx <- sort(c(lambda.idx, theta.idx, nu.idx,
                     delta.idx, beta.idx, psi.idx))
  JACa <- JACa[keep.idx, ,drop = FALSE]

  # JACb: jacobian of the function vech(VETA) = f(vech(S), theta.mm)
  #       (treating theta.mm as fixed)
  if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
    ffb <- function(x) {
      if (lavmodel@meanstructure) {
        nvar <- nrow(FIT@h1$implied$cov[[g]])
        this.ybar <- x[seq_len(nvar)]
        this.cov <- lav_matrix_vech_reverse(x[-seq_len(nvar)])
      } else {
        this.ybar <- FIT@h1$implied$mean[[g]]
        this.cov <- lav_matrix_vech_reverse(x)
      }

      # change COV/YBAR
      step1 <- STEP1
      step1$COV[[1]] <- this.cov
      if (lavmodel@meanstructure) {
        step1$YBAR[[1]] <- this.ybar
      }

      # transform data to comply with the new COV/YBAR
      Y <- FIT@Data@X[[1]]
      Ytrans <- vector("list", nblocks)
      Ytrans[[1]] <- lav_matrix_transform_mean_cov(Y, target.mean = this.ybar,
                                                   target.cov = this.cov)
      colnames(Ytrans[[1]]) <- FIT@pta$vnames$ov[[1]]

      step1 <- lav_sam_step1_local(STEP1 = step1, FIT = FIT, Y = Ytrans,
           sam.method = STEP1$sam.method, local.options = STEP1$local.options)
      if (lavmodel@meanstructure) {
        out <- c(step1$EETA[[1]], lav_matrix_vech(step1$VETA[[1]]))
      } else {
        out <- lav_matrix_vech(step1$VETA[[1]])
      }
      out
    }
    # shut off verbose
    verbose.flag <- lav_verbose()
    lav_verbose(FALSE)
    this.x <- lav_matrix_vech(FIT@h1$implied$cov[[g]])
    if (lavmodel@meanstructure) {
      this.x <- c(FIT@h1$implied$mean[[g]], this.x)
    }
    JACb <- numDeriv::jacobian(func = ffb, x = this.x)
    lav_verbose(verbose.flag)
    # lv.names1 <- STEP1$LV.NAMES[[g]]
    # lv.int.names <- FIT@pta$vnames$lv.interaction[[g]]
    # nfac <- length(lv.names1)

    # idx1 <- rep(seq_len(nfac), each = nfac)
    # idx2 <- rep(seq_len(nfac), times = nfac)

    # NAMES <- paste(lv.names[idx1], lv.names[idx2], sep = ":")



    # JACb <- lav_sam_step1_local_jac_var2(ybar = STEP1$YBAR[[g]],
    #   S = STEP1$COV[[g]], M = STEP1$M[[g]], NU = STEP1$NU[[g]])
    # JACb <- JACb[STEP1$lv.keep2, , drop = FALSE]
    # if (lavmodel@meanstructure) {
    #   tmp <- lav_sam_step1_local_jac_mean2(ybar = STEP1$YBAR[[g]],
    #     M = STEP1$M[[g]], NU = STEP1$NU[[g]])
    #   dd
    #   JACb <- lav_matrix_bdiag(tmp, JACb)
    # }

  } else { # no latent interactions
    Mb <- STEP1$M[[g]]
    MbxMb <- Mb %x% Mb
    row.idx <- lav_matrix_vech_idx(nrow(Mb))
    JACb <- lav_matrix_duplication_post(MbxMb)[row.idx,,drop = FALSE]
    if (lavmodel@meanstructure) {
      JACb <- lav_matrix_bdiag(Mb, JACb)
    }
  }

  # JACc: jacobian of the function vech(VETA) = f(theta.mm, vech(S))
  #       (treating vech(S) as fixed)

  # ffc <- function(x, YBAR = NULL, COV = NULL, b = 1L) {
  #   # x only contains the LAMBDA/THETA/NU elements
  #   PT$est[keep.idx] <- x
  #   # get all free parameters (for lav_model_set_parameters)
  #   x.free <- PT$est[PT$free > 0L & !duplicated(PT$free)]
  #   this.model <- lav_model_set_parameters(lavmodel, x = x.free)
  #   lambda.idx <- which(names(this.model@GLIST) == "lambda")[b]
  #   theta.idx  <- which(names(this.model@GLIST) ==  "theta")[b]
  #   LAMBDA <- this.model@GLIST[[lambda.idx]]
  #   THETA  <- this.model@GLIST[[ theta.idx]]
  #   Mb <- lav_sam_mapping_matrix(LAMBDA = LAMBDA,
  #                                THETA = THETA, S = COV,
  #                                method = local.options$M.method)
  #   # handle observed-only variables
  #   dummy.ov.idx <- c(
  #     FIT@Model@ov.x.dummy.ov.idx[[b]],
  #     FIT@Model@ov.y.dummy.ov.idx[[b]]
  #   )
  #   dummy.lv.idx <- c(
  #     FIT@Model@ov.x.dummy.lv.idx[[b]],
  #     FIT@Model@ov.y.dummy.lv.idx[[b]]
  #   )
  #   if (length(dummy.ov.idx)) {
  #     Mb[dummy.lv.idx, ] <- 0
  #     Mb[cbind(dummy.lv.idx, dummy.ov.idx)] <- 1
  #   }
  #   MSM <- Mb %*% COV %*% t(Mb)
  #   MTM <- Mb %*% THETA %*% t(Mb)
  #   VETA <- MSM - MTM

  #   if (lavmodel@meanstructure) {
  #     nu.idx <- which(names(this.model@GLIST) ==  "nu")[b]
  #     NU <- this.model@GLIST[[nu.idx]]
  #     EETA <- lav_sam_eeta(M = Mb, YBAR = YBAR, NU = NU)
  #     out <- c(EETA, lav_matrix_vech(VETA))
  #   } else {
  #     out <- lav_matrix_vech(VETA)
  #   }

  #   out
  # }

  # current point estimates
  # PT <- STEP1$PT
  # x <- PT$est[keep.idx]
  # if (lavmodel@meanstructure) {
  #   JACc <- numDeriv::jacobian(func = ffc, x = x, YBAR = STEP1$YBAR[[1]],
  #                              COV = STEP1$COV[[1]])
  # } else {
  #   JACc <- numDeriv::jacobian(func = ffc, x = x, COV = STEP1$COV[[1]])
  # }

  # calling lav_sam_step1_local() directly
  ffc <- function(x) {
    step1 <- STEP1

    # fill in new 'x' values in PT
    step1$PT$est[keep.idx] <- x
    step1 <- lav_sam_step1_local(STEP1 = step1, FIT = FIT,
           sam.method = STEP1$sam.method, local.options = STEP1$local.options)
    if (lavmodel@meanstructure) {
      out <- c(step1$EETA[[1]], lav_matrix_vech(step1$VETA[[1]]))
    } else {
      out <- lav_matrix_vech(step1$VETA[[1]])
    }
    out
  }
  # shut off verbose
  verbose.flag <- lav_verbose()
  lav_verbose(FALSE)
  JACc <- numDeriv::jacobian(func = ffc, x = STEP1$PT$est[keep.idx])
  lav_verbose(verbose.flag)

  # assemble JAC
  JAC <- (JACc %*% JACa) + JACb

  # eventually, this will be a list per group/block
  list(JAC)
}

lav_sam_step1_local_jac_var2 <- function(ybar = NULL, S = NULL,
                                         M = NULL, NU = NULL) {

  nfac <- nrow(M)
  K.nfac <- lav_matrix_commutation(nfac, nfac)
  IK <- diag(nfac*nfac) + K.nfac

  MEAN.FS <- M %*% (ybar - NU)
  MSM <- M %*% S %*% t(M)

  part1a <- ( diag(nfac) %x% matrix(lav_matrix_vec(IK %*% (MSM %x% diag(nfac)) %*% K.nfac), nfac^3, nfac) +
        matrix(diag(nfac) %x% (IK %*% (diag(nfac) %x% MSM)), nfac^4, nfac^2) )
  tmp1 <- part1a %*% IK %*% (diag(nfac) %x% MEAN.FS) %*% M

# S part
  A <- tcrossprod(MEAN.FS)
  part1b <- diag(nfac) %x% matrix(lav_matrix_vec(IK %*% (A %x% diag(nfac)) %*% K.nfac), nfac^3, nfac)
  part1c <- matrix(diag(nfac) %x% (IK %*% (diag(nfac) %x% A)), nfac^4, nfac^2)
  tmp2 <- (part1a + part1b + part1c) %*% lav_matrix_duplication_post(M %x% M)

  # together
  JAC.S.analytic <- cbind(tmp1, tmp2)
  JAC.S.analytic
}

lav_sam_step1_local_jac_mean2 <- function(ybar = NULL, M = NULL, NU = NULL) {
  nfac <- nrow(M)
  K.nfac <- lav_matrix_commutation(nfac, nfac)
  IK <- diag(nfac*nfac) + K.nfac

  MEAN.FS <- M %*% (ybar - NU)

  tmp1 <- IK %*% (diag(nfac) %x% MEAN.FS) %*% M
  tmp2 <- lav_matrix_duplication_post(M %x% M)
  JAC.M.analytic <- cbind(tmp1, tmp2)
  JAC.M.analytic
}

lav_sam_gamma_add <- function(STEP1 = NULL, FIT = NULL, group = 1L) {

  lavdata <- FIT@Data
  lavsamplestats <- FIT@SampleStats
  lavmodel <- FIT@Model
  lavpta <- FIT@pta
  nblocks <- lavpta$nblocks

  local.options <- STEP1$local.options
  sam.method <- STEP1$sam.method

  ngroups <- lavdata@ngroups
  if (ngroups > 1L) {
    stop("IJ local SEs: not available with multiple groups!\n")
  }
  g <- group
  Y <- FIT@Data@X[[g]]
  N <- nrow(Y)

  # NAMES + lv.keep
  lv.names <- STEP1$LV.NAMES[[1]]
  lv.names <- c("..int..", lv.names)
  nfac <- length(lv.names)
  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)

  K.nfac <- lav_matrix_commutation(nfac, nfac)
  IK <- diag(nfac * nfac) + K.nfac

  NAMES <- paste(lv.names[idx1], lv.names[idx2], sep = ":")
  NAMES[seq_len(nfac)] <- lv.names
  lv.keep <- colnames(STEP1$VETA[[1]])
  FS.mean <- STEP1$FS.mean[[1]]
  keep.idx <- which(NAMES %in% lv.keep)

  theta.to.eetavetai <- function(x, i = 1L) {
    PT <- STEP1$PT
    PT$est[step1.idx] <- x
    this.lavmodel <- lav_model_set_parameters(lavmodel,
                       x = PT$est[PT$free > 0 & !duplicated(PT$free)])
    this.nu     <- this.lavmodel@GLIST$nu
    rm.idx <- lavpta$vidx$lv.interaction[[1]]
    # no interaction columns!
    this.lambda <- this.lavmodel@GLIST$lambda[, -rm.idx,drop = FALSE]
    this.theta  <- this.lavmodel@GLIST$theta
    this.M <- lav_sam_mapping_matrix(LAMBDA = this.lambda,
                                     THETA = this.theta,
                                     S = STEP1$COV[[1]],
                                     method = STEP1$local.options$M.method)
    MTM <- this.M %*% this.theta %*% t(this.M)
    MTM <- lav_matrix_bdiag(0,MTM)

    fi <- this.M %*% (Y[i,] - this.nu); fi <- rbind(1, fi); fii <- drop(fi)
    fi2 <- (fii[idx1]*fii[idx2])

    tmp <- ( ((tcrossprod(fi) - MTM) %x% MTM) +
             (MTM %x% (tcrossprod(fi) - MTM)) +
             lav_matrix_commutation_post((tcrossprod(fi) - MTM) %x% MTM) +
             lav_matrix_commutation_pre((tcrossprod(fi) - MTM) %x% MTM) +
             (IK %*% (MTM %x% MTM)) )

    # we have no access to FS2.mean, so we need to reduce the matrices
    # right away
    f.star <- tcrossprod(fi2[keep.idx] - FS.mean)
    e.star <- STEP1$lambda[[1]] * tmp[keep.idx, keep.idx, drop = FALSE]
    iveta2 <- f.star - e.star
    #iveta  <- iveta2[seq_len(nfac - 1), seq_len(nfac - 1)]

    ieeta2 <- ( lav_matrix_vec(tcrossprod(fi))[keep.idx] -
                lav_matrix_vec(tcrossprod(FS.mean))[keep.idx] +
                (FS.mean %x% FS.mean)[keep.idx] -
                STEP1$lambda[[1]] * lav_matrix_vec(MTM)[keep.idx] )

    c(ieeta2, lav_matrix_vech(iveta2))
  } # single 'i'

  step1.idx <- which(STEP1$PT$free %in% STEP1$step1.free.idx)
  x.step1 <- STEP1$PT$est[step1.idx]
  try.one <- theta.to.eetavetai(x = x.step1, i = 1)
  CVETA <- matrix(0, nrow = length(try.one), ncol = length(x.step1))
  for(i in 1:N) {
    tmp <- numDeriv::jacobian(func = theta.to.eetavetai, x = x.step1, i = i)
    CVETA <- CVETA + 1/N * tmp
  }

  Gamma.addition <- N * (CVETA %*% STEP1$Sigma.11 %*% t(CVETA))
  Gamma.addition
}
