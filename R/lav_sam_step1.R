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
    mm.list.per.block <- lav_sam_get_mmlist(FIT)
    # we get a list PER BLOCK
    # for now, we only keep the first block
    mm.list <- mm.list.per.block[[1]]
    nMMblocks <- length(mm.list)
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

    # add step1.free.idx
    par.idx <- PT.free[mm.idx[ptm.idx]]
    # store (ordered) indices in step1.free.idx
    this.mm.idx <- sort.int(par.idx)
    step1.free.idx <- c(step1.free.idx, this.mm.idx) # all combined


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
      keep.idx <- PTM.free[ptm.idx]
      # par.idx <- PT.free[ seq_len(length(PT$lhs)) %in% mm.idx &
      #                    PT$free > 0L ]
      # keep.idx <- PTM.free[ PTM$free > 0 & PTM$user != 3L ]
      Sigma.11[par.idx, par.idx] <-
        sigma.11[keep.idx, keep.idx, drop = FALSE]
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

