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
  lsam.analytic.flag <- rep(TRUE, nblocks)
  L.veta <- vector("list", length = nblocks)
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
      if (sam.method == "local" && !lv.interaction.flag) {
        lsam.analytic.flag[b] <- FALSE
        # we will try an iterative solution
      } else {
        # eg cfsr or lv interactions: no idea what to do here (yet)
        print(this.lambda)
        lav_msg_stop(gettext(
          "LAMBDA has no full column rank. Please use sam.method = global"))
      }
    }

    # if lambda has full rank, check if cov.lv is unrestricted
    if (!lsam.analytic.flag[b]) {
      VETA.symbolic <- lav_sam_veta_partable(FIT, block = b)
      # this is the tricky thing: which rows/cols should we remove?
      # none for now
      if (FIT@Options$std.lv) {
        veta.symbolic <- lav_matrix_vech(VETA.symbolic, diagonal = FALSE)
      } else {
        veta.symbolic <- lav_matrix_vech(VETA.symbolic, diagonal = TRUE)
      }
      if (any(veta.symbolic == 0)) {
        lsam.analytic.flag[b] <- FALSE
        nfac <- ncol(VETA.symbolic)
        m.free <- which(VETA.symbolic != 0)
        tmp <- lav_matrix_vech_reverse(lav_matrix_vech_idx(nfac))
        x.free <- tmp[which(VETA.symbolic != 0)]
        unique.idx <- unique(x.free)
        row.idx <- match(x.free, unique.idx)
        L.psi <- matrix(0L, nrow = nfac * nfac, ncol = length(unique.idx))
        IDX <- cbind(m.free, row.idx)
        L.psi[IDX] <- 1L
        L.veta[[b]] <- L.psi
      }
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
    if (lsam.analytic.flag[b]) {
      Mb <- lav_sam_mapping_matrix(
        LAMBDA = this.lambda,
        THETA = THETA[[b]],
        S = COV, S.inv = ICOV,
        method = local.M.method
      )
    } else {
      Mb <- matrix(as.numeric(NA), ncol(this.lambda), nrow(this.lambda))
    }
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
      if (lsam.analytic.flag[b]) {
        EETA[[b]] <- lav_sam_eeta(M = Mb, YBAR = YBAR, NU = NU[[b]])
      } else {
        # EETA is constrained somehow
        lav_msg_stop(gettext("not ready yet"))
        # EETA[[b]] <- lav_sam_eeta_con(YBAR = YBAR, LAMBDA = LAMBDA[[b]],
        #                               THETA = THETA[[b]],
        #                               L.veta = L.veta[[b]])

      }
	  FS.mean[[b]] <- EETA[[b]] # ok if no interaction
    }

    # compute VETA
    if (sam.method == "local") {
      if (lsam.analytic.flag[b]) {
        tmp <- lav_sam_veta(
          M = Mb, S = COV, THETA = THETA[[b]],
          alpha.correction = local.options[["alpha.correction"]],
          lambda.correction = local.options[["lambda.correction"]],
          N <- FIT@SampleStats@nobs[[this.group]],
          dummy.lv.idx = dummy.lv.idx,
          extra = TRUE
        )
        VETA[[b]] <- tmp[, , drop = FALSE] # drop attributes
        alpha[[b]]  <- attr(tmp, "alpha")
        lambda[[b]] <- attr(tmp, "lambda.star")
	    MSM..[[b]]  <- attr(tmp, "MSM")
	    MTM..[[b]]  <- attr(tmp, "MTM")
      } else {
        # VETA is constrained somehow
        VETA[[b]] <- lav_sam_veta_con(S = COV, LAMBDA = LAMBDA[[b]],
                                      THETA = THETA[[b]], L.veta = L.veta[[b]],
                                      local.M.method = local.M.method)
        alpha[[b]]  <- as.numeric(NA)
        lambda[[b]] <- as.numeric(NA)
        MSM..[[b]]  <- matrix(0, 0, 0)
        MTM..[[b]]  <- matrix(0, 0, 0)
      }
    } else if (sam.method == "cfsr") {
      # first, we need to 'true' VETA (to get Sigma)
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
    if (lsam.analytic.flag[b]) {
      REL[[b]] <- diag(VETA[[b]]) / diag(MSM..[[b]]) #!
    } else {
      REL[[b]] <- as.numeric(NA)
    }

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


lav_sam_step1_local_jac <- function(STEP1 = NULL, FIT = NULL, P.only = FALSE,
                                    return.jac = FALSE) {

  lavdata <- FIT@Data
  lavsamplestats <- FIT@SampleStats
  lavmodel <- FIT@Model
  lavpta <- FIT@pta
  nblocks <- lavpta$nblocks

  local.options <- STEP1$local.options
  sam.method <- STEP1$sam.method

  ngroups <- lavdata@ngroups
  if (ngroups > 1L) {
    lav_msg_stop(gettext("IJ local SEs: not available with multiple groups!\n"))
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
    lav_msg_stop(gettext("IJ local SEs: not available for the categorical setting (yet)!\n"))
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
    if (P.only && FIT@Options$information[1] == "expected") {
      # for twostep.robust
      mm.inv.observed <- lavTech(fit.mm.block, "inverted.information.expected")
    } else {
      mm.inv.observed  <- lavTech(fit.mm.block, "inverted.information.observed")
    }

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
  if (P.only) {
    return(JACa)
  }

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

  if (return.jac) {
    attr(JAC, "JACa") <- JACa
    attr(JAC, "JACb") <- JACb
    attr(JAC, "JACc") <- JACc
  }

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

lav_sam_gamma_add_numerical <- function(STEP1 = NULL, FIT = NULL, group = 1L) {

  lavdata <- FIT@Data
  lavsamplestats <- FIT@SampleStats
  lavmodel <- FIT@Model
  lavpta <- FIT@pta
  nblocks <- lavpta$nblocks

  local.options <- STEP1$local.options
  sam.method <- STEP1$sam.method

  ngroups <- lavdata@ngroups
  if (ngroups > 1L) {
    lav_msg_stop(gettext("IJ local SEs: not available with multiple groups!\n"))
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
  keep.idx <- match(lv.keep, NAMES)

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

# semi-analytic version
# YR 4 June 2025: works, but still needs cleanup + avoid redundant calculations
# (eg when multiplied with a matrix with many zero rows/cols)
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
    lav_msg_stop(gettext("IJ local SEs: not available with multiple groups!\n"))
  }
  g <- group
  Y <- FIT@Data@X[[g]]
  N <- nrow(Y)
  P <- ncol(Y)

  # NAMES + lv.keep
  lv.names <- STEP1$LV.NAMES[[1]]
  lv.names <- c("..int..", lv.names)
  nfac <- length(lv.names)
  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)

  K.nfac <- lav_matrix_commutation(nfac, nfac)
  IK <- diag(nfac * nfac) + K.nfac
  D <- lav_matrix_duplication(nfac)

  NAMES <- paste(lv.names[idx1], lv.names[idx2], sep = ":")
  NAMES[seq_len(nfac)] <- lv.names
  lv.keep <- colnames(STEP1$VETA[[1]])
  FS.mean <- STEP1$FS.mean[[1]]
  keep.idx <- match(lv.keep, NAMES)

  # step 1 free parameters
  step1.idx <- which(STEP1$PT$free %in% STEP1$step1.free.idx)
  x.step1 <- STEP1$PT$est[step1.idx]

  this.nu <- STEP1$NU[[1]]
  this.M  <- STEP1$M[[1]]
  this.MTM <- lav_matrix_bdiag(0, STEP1$MTM[[1]][seq_len(nfac - 1), seq_len(nfac - 1)])
  this <- c(this.nu, lav_matrix_vec(this.M), lav_matrix_vech(this.MTM))
  INDEX <- matrix(seq_len(nfac^2 * nfac^2), nfac^2, nfac^2)

  x2this <- function(x) {
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

    out <- c(this.nu, lav_matrix_vec(this.M), lav_matrix_vech(MTM))
    out
  }
  JAC.x2this <- numDeriv::jacobian(func = x2this, x = x.step1)
  # 46 x 24

  # JAC.eetai2.this
  JAC.eeta2.this <- matrix(0, nrow = length(STEP1$EETA[[1]]), ncol = length(this))
  pstar <- nfac * (nfac + 1) / 2
  idx <- P + length(this.M) + seq_len(pstar)
  JAC.eeta2.this[,idx] <- (-diag(nfac^2) %*% D)[keep.idx,]

  # define function for JAC.veta2.fi
  get.JAC.veta2.fi <- function(x, MTM, FS.MEAN = NULL, lambda.star = 1, keep.idx = NULL) {
    fi <- drop(fi)
    P <- length(fi)
    stopifnot(P == ncol(MTM))

    #INDEX <- matrix(seq_len(P^2 * P^2), P^2, P^2)

    # (tcrossprod(fi) - MTM) %x% MTM) -- part A
    block.a <- fi %x% do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% MTM[,i]))
    long.a <- diag(P) %x% lav_matrix_vec(fi %x% MTM)
    big.a <- block.a + long.a
    part.a <- big.a[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # (MTM %x% (tcrossprod(fi) - MTM)) -- part B
    block.b <- lav_matrix_vec(fi %x% MTM) %x% diag(P)
    long.b <- do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% MTM[,i])) %x% fi
    big.b <- block.b + long.b
    part.b <- big.b[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # lav_matrix_commutation_post((tcrossprod(fi) - MTM) %x% MTM) -- part C
    block.c <- do.call("rbind", lapply(seq_len(P), function(i) fi %x% (diag(P) %x% MTM[,i])))
    long.c  <- do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% (fi %x% MTM[,i])))
    big.c <- block.c + long.c
    part.c <- big.c[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # lav_matrix_commutation_pre((tcrossprod(fi) - MTM) %x% MTM) -- part D
    long.d <- diag(P) %x% lav_matrix_vec(MTM %x% fi)
    block.d <- (fi %x% lav_matrix_vec(MTM)) %x% diag(P)
    big.d <- block.d + long.d
    part.d <- big.d[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # t1 (tcrossprod(fi2[keep.idx] - FS.mean) only) -- part E
    idx1 <- rep(seq_len(P), each = P)
    idx2 <- rep(seq_len(P), times = P)
    fi2 <- (fi[idx1]*fi[idx2])
    fik <- fi2[keep.idx]

    fi <- as.matrix(fi)
    tt <- ((fi %x% diag(P))[keep.idx,] + (diag(P) %x% fi)[keep.idx,])
    tmp <- ((fik - FS.mean) %x% tt) + (tt %x% (fik - FS.mean))
    part.e <- tmp[lav_matrix_vech_idx(length(fik)),]

    final <- part.e - lambda.star * (part.a + part.b + part.c + part.d)
    final
  }

  # define function for JAC.veta2.this.i
  get.JAC.veta2.this <- function(x, fi = NULL, lambda.star = 1, keep.idx = NULL) {
    nvar <- ncol(STEP1$M[[1]])
    nfac <- nrow(STEP1$M[[1]])
    this.nu <- x[seq_len(nvar)]; x <- x[-seq_len(nvar)]
    this.M  <- matrix(x[seq_len(nvar*nfac)], nrow = nfac, ncol = nvar)
    x <- x[-seq_len(nvar*nfac)]
    MTM <- lav_matrix_vech_reverse(x)
    fi <- as.matrix(fi)
    P <- nrow(fi)

    # part a: ((tcrossprod(fi) - MTM) %x% MTM)
    long.a <- diag(P) %x% do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% -MTM[,i]))
    block.a <- do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% (tcrossprod(fi) - MTM)[,i])) %x% diag(P)
    big.a.vec <- long.a + block.a
    big.a <- big.a.vec %*% D
    part.a <- big.a[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # part b: (MTM %x% (tcrossprod(fi) - MTM))
    block.b <- do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% -MTM[,i])) %x% diag(P)
    long.b <- diag(P) %x% do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% (tcrossprod(fi) - MTM)[,i]))
    big.b.vec <- block.b + long.b
    big.b <- big.b.vec %*% D
    part.b <- big.b[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # part c: lav_matrix_commutation_post((tcrossprod(fi) - MTM) %x% MTM)
    block.c <- diag(P) %x% do.call("rbind", lapply(seq_len(P), function(i) ((tcrossprod(fi) - MTM)[,i]) %x% diag(P)))
    long.c <- do.call("rbind", lapply(seq_len(P), function(i) diag(P*P) %x% (-MTM[,i]) ))
    big.c.vec <- block.c + long.c
    big.c <- big.c.vec %*% D
    part.c <- big.c[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # part d: lav_matrix_commutation_pre((tcrossprod(fi) - MTM) %x% MTM)
    block.d <- diag(P) %x% do.call("rbind", lapply(seq_len(P), function(i) -MTM[,i] %x% diag(P)))
    long.d <- do.call("rbind", lapply(seq_len(P), function(i) diag(P*P) %x% ((tcrossprod(fi) - MTM)[,i]) ))
    big.d.vec <- block.d + long.d
    big.d <- big.d.vec %*% D
    part.d <- big.d[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    # part e: (IK %*% (MTM %x% MTM))
    K <- lav_matrix_commutation(P, P)
    e1 <- diag(P) %x% do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% MTM[,i]))
    e2 <- do.call("rbind", lapply(seq_len(P), function(i) diag(P*P) %x% MTM[,i] ))
    e3 <- do.call("rbind", lapply(seq_len(P), function(i) diag(P) %x% MTM[,i])) %x% diag(P)
    e4 <- diag(P) %x% do.call("rbind", lapply(seq_len(P), function(i) K %*% (diag(P) %x% MTM[,i])))
    big.e.vec <- e1 + e2 + e3 + e4
    big.e <- big.e.vec %*% D
    part.e <-  big.e[lav_matrix_vech(INDEX[keep.idx, keep.idx]),]

    final <- -1 * lambda.star * (part.a + part.b + part.c + part.d + part.e)
    final
  }

  n_eeta_veta <- length(STEP1$EETA[[1]]) + length(lav_matrix_vech(STEP1$VETA[[1]]))
  CVETA <- matrix(0, nrow = n_eeta_veta, ncol = length(x.step1))
  for(i in 1:N) {
    # factor score
    fi <- rbind(1, this.M %*% (Y[i,] - this.nu))

    #tmp <- numDeriv::jacobian(func = theta.to.eetavetai, x = x.step1, i = i)
    JAC.this2fi.i <- rbind(0, cbind(-1 * this.M,
                           t(as.matrix(Y[i,] - drop(this.nu))) %x% diag(nrow(this.M)),
                           matrix(0, nrow = nrow(this.M), ncol = length(lav_matrix_vech(this.MTM)))))

    JAC.eeta2.fi.i <- ((fi %x% diag(nfac)) + (diag(nfac) %x% fi))[keep.idx,]
    EETA2 <- ((JAC.eeta2.fi.i %*% JAC.this2fi.i) + JAC.eeta2.this) %*% JAC.x2this

    JAC.veta2.fi.i <- get.JAC.veta2.fi(x = fi, MTM = this.MTM, FS.MEAN = FS.mean,
                                       lambda.star = STEP1$lambda[[1]], keep.idx = keep.idx)
    JAC.veta2.this.i <- matrix(0, nrow = nrow(JAC.veta2.fi.i), ncol = ncol(JAC.this2fi.i))
    JAC.veta2.this.i[,idx] <- get.JAC.veta2.this(x = this, fi = fi,
                                                 lambda.star = STEP1$lambda[[1]],
                                                 keep.idx = keep.idx)
    VETA2 <- ((JAC.veta2.fi.i %*% JAC.this2fi.i) + JAC.veta2.this.i) %*% JAC.x2this

    CVETA <- CVETA + 1/N * rbind(EETA2, VETA2)
  }

  Gamma.addition <- N * (CVETA %*% STEP1$Sigma.11 %*% t(CVETA))
  Gamma.addition
}


