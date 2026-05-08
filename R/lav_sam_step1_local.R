## STEP 1b: compute Var(eta) and E(eta) per block
##          only needed for local/fsr approach!
lav_sam_step1_local <- function(step1 = NULL, fit = NULL, y = NULL,
                                sam_method = "local",
                                local_options = list(
                                  M.method = "ML",
                                  lambda.correction = TRUE,
                                  alpha.correction = 0L,
                                  twolevel.method = "h1"
                                ),
                                return_cov_iveta2 = TRUE,
                                return_fs = FALSE) {
  # local.M.method
  local_m_method <- toupper(local_options[["M.method"]])
  if (!local_m_method %in% c("GLS", "ML", "ULS")) {
    lav_msg_stop(gettext(
      "local option M.method should be one of ML, GLS or ULS."))
  }

  lavoptions <- fit@Options
  lavpta <- fit@pta
  nblocks <- lavpta$nblocks

  # flags
  lv_interaction_flag <- FALSE
  lv_higherorder_flag <- FALSE
  if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
    lv_interaction_flag <- TRUE
  }
  if (length(unlist(lavpta$vnames$lv.ind)) > 0L) {
    lv_higherorder_flag <- TRUE
  }

  if (lav_verbose()) {
    cat("Constructing the mapping matrix using the ",
      local_m_method, " method ... ",
      sep = ""
    )
  }

  # all the measurement parameters are already stored in PT
  pt_1 <- step1$PT
  if (fit@Model@ceq.simple.only) {
    x_free <- pt_1$est[pt_1$free > 0 & !duplicated(pt_1$free)]
  } else {
    x_free <- pt_1$est[pt_1$free > 0]
  }
  # check for NA values (eg in BETA); set them to zero
  x_free[!is.finite(x_free)] <- 0

  lavmodel_tmp <- lav_model_set_parameters(fit@Model, x = x_free)
  mm_lambda <- mm_theta <- mm_beta <- mm_psi <- mm_nu <- mm_delta <- NULL

  # create LAMBDA
  lambda_idx <- which(names(fit@Model@GLIST) == "lambda")
  mm_lambda <- lavmodel_tmp@GLIST[lambda_idx]

  # create THETA
  theta_idx <- which(names(fit@Model@GLIST) == "theta")
  mm_theta <- lavmodel_tmp@GLIST[theta_idx]

  # NU
  if (fit@Model@meanstructure) {
    nu_idx <- which(names(fit@Model@GLIST) == "nu")
    mm_nu <- lavmodel_tmp@GLIST[nu_idx]
  }

  # DELTA
  if (fit@Model@categorical || fit@Model@correlation) {
    delta_idx <- which(names(fit@Model@GLIST) == "delta")
    mm_delta <- lavmodel_tmp@GLIST[delta_idx]
  }

  # BETA/PSI
  if (lv_higherorder_flag) {
    beta_idx <- which(names(fit@Model@GLIST) == "beta")
    mm_beta <- lavmodel_tmp@GLIST[beta_idx]
    psi_idx <- which(names(fit@Model@GLIST) == "psi")
    mm_psi <- lavmodel_tmp@GLIST[psi_idx]
  }

  # GAMMA (only for names, if conditional.x)
  #if (FIT@Model@conditional.x) {
  #  gamma.idx <- which(names(FIT@Model@GLIST) == "gamma")
  #}

  # handle dummy's + higher-order + rank-deficient
  lsam_analytic_flag <- rep(TRUE, nblocks)
  l_veta <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    # new in 0.6-10: check if any indicators are also involved
    # in the structural part; if so, set THETA row/col to zero
    # and make sure LAMBDA element is correctly set
    # (we also need to adjust M)
    dummy_ov_idx <- fit@Model@ov.y.dummy.ov.idx[[b]]
    dummy_lv_idx <- fit@Model@ov.y.dummy.lv.idx[[b]]
    if (length(dummy_ov_idx)) {
      mm_theta[[b]][dummy_ov_idx, ] <- 0
      mm_theta[[b]][, dummy_ov_idx] <- 0
      mm_lambda[[b]][dummy_ov_idx, ] <- 0
      mm_lambda[[b]][cbind(dummy_ov_idx, dummy_lv_idx)] <- 1
    }
    if (fit@Model@meanstructure) {
      if (length(dummy_ov_idx)) {
        mm_nu[[b]][dummy_ov_idx, 1] <- 0
      }
    }

    # get ALL lv names (including dummy ov.x/ov.y)
    lv_names <- fit@Model@dimNames[[lambda_idx[b]]][[2L]]

    # if conditional.x, we must add the ov.names.x manually
    # if (FIT@Model@conditional.x) {
    #   exo.names <- FIT@Model@dimNames[[gamma.idx[b]]][[2L]]
    #   lv.names <- c(lv.names, exo.names)
    # }

    # handle higher-order factors here
    if (length(lavpta$vidx$lv.ind[[b]]) > 0L) {
      lv_ind_names <- lavpta$vnames$lv.ind[[b]]
      lv_target <- lv_names[!lv_names %in% lv_ind_names]

      target_idx <- match(lv_target, lv_names)
      other_idx <- seq_along(lv_names)[-target_idx]

      this_beta <- mm_beta[[b]]
      this_beta[is.na(this_beta)] <- 0

      ib <- diag(nrow(this_beta)) - this_beta
      ib_inv <- solve(ib)
      lb_inv <- mm_lambda[[b]] %*% ib_inv

      # replace LAMBDA
      mm_lambda[[b]] <- lb_inv[, target_idx, drop = FALSE]

      psi_other <- mm_psi[[b]][other_idx, other_idx, drop = FALSE]
      lb_inv2 <- lb_inv[, other_idx, drop = FALSE]

      # replace THETA
      mm_theta[[b]] <- lb_inv2 %*% psi_other %*% t(lb_inv2) + mm_theta[[b]]
    }

    # check if LAMBDA has full column rank
    this_lambda <- mm_lambda[[b]]
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      if (length(lavpta$vidx$lv.ind[[b]]) > 0L) {
        rm_idx <- c(match(lavpta$vnames$lv.ind[[b]], lv_names),
                    match(lavpta$vnames$lv.interaction[[b]], lv_names))
        this_lambda <- this_lambda[, -rm_idx, drop = FALSE]
      } else {
        rm_idx <- match(lavpta$vnames$lv.interaction[[b]], lv_names)
        this_lambda <- this_lambda[, -rm_idx, drop = FALSE]
      }
    }
    if (qr(this_lambda)$rank < ncol(this_lambda)) {
      if (sam_method == "local" && !lv_interaction_flag) {
        lsam_analytic_flag[b] <- FALSE
        # we will try an iterative solution
      } else {
        # eg cfsr or lv interactions: no idea what to do here (yet)
        print(this_lambda)
        lav_msg_stop(gettext(
          "LAMBDA has no full column rank. Please use sam.method = global"))
      }
    }

    # if lambda has full rank, check if cov.lv is unrestricted
    if (!lsam_analytic_flag[b]) {
      veta_symbolic_1 <- lav_sam_veta_partable(fit, block = b)
      # this is the tricky thing: which rows/cols should we remove?
      # none for now
      if (fit@Options$std.lv) {
        veta_symbolic <- lav_matrix_vech(veta_symbolic_1, diagonal = FALSE)
      } else {
        veta_symbolic <- lav_matrix_vech(veta_symbolic_1, diagonal = TRUE)
      }
      if (any(veta_symbolic == 0)) {
        lsam_analytic_flag[b] <- FALSE
        nfac <- ncol(veta_symbolic_1)
        m_free <- which(veta_symbolic_1 != 0)
        tmp <- lav_matrix_vech_reverse(lav_matrix_vech_idx(nfac))
        x_free <- tmp[which(veta_symbolic_1 != 0)]
        unique_idx <- unique(x_free)
        row_idx <- match(x_free, unique_idx)
        l_psi <- matrix(0L, nrow = nfac * nfac, ncol = length(unique_idx))
        idx <- cbind(m_free, row_idx)
        l_psi[idx] <- 1L
        l_veta[[b]] <- l_psi
      }
    }
  } # b

  # store LAMBDA/THETA/NU per block
  step1$LAMBDA <- mm_lambda
  step1$THETA <- mm_theta
  if (fit@Model@meanstructure) {
    step1$NU <- mm_nu
  }
  if (fit@Model@categorical || fit@Model@correlation) {
    step1$DELTA <- mm_delta
  }

  veta     <- vector("list", nblocks)
  msm_    <- vector("list", nblocks)
  mtm_    <- vector("list", nblocks)
  fs_mean  <- vector("list", nblocks)
  #FS.gamma <- vector("list", nblocks)
  fs       <- vector("list", nblocks)
  cov_iveta2 <- vector("list", nblocks)
  rel <- vector("list", nblocks)
  alpha <- vector("list", nblocks)
  lambda <- vector("list", nblocks)
  if (lavoptions$meanstructure) {
    eeta <- vector("list", nblocks)
  } else {
    eeta <- NULL
  }
  #fs.outlier.idx <- vector("list", nblocks)
  m <- vector("list", nblocks)
  lv_names_1 <- vector("list", nblocks)

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
    this_group <- floor(b / fit@Data@nlevels + 0.5)

    # lv.names, including dummy-lv covariates
    psi_idx <- which(names(fit@Model@GLIST) == "psi")[b]
    lv_names_b <- fit@Model@dimNames[[psi_idx]][[1L]] # including dummy/inter.
    # if (FIT@Model@conditional.x) {
    #   exo.names <- FIT@Model@dimNames[[gamma.idx[b]]][[2L]]
    #   lv.names.b <- c(lv.names.b, exo.names)
    # }
    rm_idx <- integer(0L)

    # higher-order? remove lower-order factors
    if (lv_higherorder_flag && length(lavpta$vnames$lv.ind[[b]]) > 0L) {
      rm_idx <- c(rm_idx, match(lavpta$vnames$lv.ind[[b]], lv_names_b))
    }

    # interaction terms? remove them for VETA
    if (lv_interaction_flag && length(lavpta$vnames$lv.interaction[[b]]) > 0L) {
      rm_idx <- c(rm_idx, match(lavpta$vnames$lv.interaction[[b]], lv_names_b))
      lv_int_names <- lavpta$vnames$lv.interaction[[b]]
    }

    # final names for EETA/VETA (not including interaction terms!)
    lv_names1 <- lv_names_b
    if (length(rm_idx) > 0L) {
      lv_names1 <- lv_names_b[-rm_idx]
    }
    lv_names_1[[b]] <- lv_names1

    # get sample statistics for this block
    cov_1 <- step1$COV[[b]]
    ybar <- drop(step1$YBAR[[b]])

    # rescale COV?
    if (fit@Data@nlevels == 1L &&
        (fit@Model@categorical || fit@Model@correlation)) {
      scale_vector <- 1 / (drop(mm_delta[[b]]))
      cov_1 <- scale_vector * cov_1 * rep(scale_vector, each = ncol(cov_1))
      ybar <- scale_vector * ybar # Checkme!
    }

    # do we need ICOV?
    if (local_m_method == "GLS") {
      if (fit@Options$sample.cov.rescale) {
        # get unbiased S
        n <- fit@SampleStats@nobs[[this_group]]
        cov_unbiased <- cov_1 * n / (n - 1)
        icov <- solve(cov_unbiased)
      } else {
        icov <- solve(cov_1)
      }
    }

    # compute mapping matrix 'M'
    this_lambda <- mm_lambda[[b]]
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      this_lambda <- this_lambda[, -lavpta$vidx$lv.interaction[[b]]]
    }
    if (lsam_analytic_flag[b]) {
      mb <- lav_sam_mapping_matrix(
        mm_lambda = this_lambda,
        mm_theta = mm_theta[[b]],
        s = cov_1, s_inv = icov,
        method = local_m_method
      )
    } else {
      mb <- matrix(as.numeric(NA), ncol(this_lambda), nrow(this_lambda))
    }
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      tmp <- mb
      mb <- matrix(0, nrow = ncol(mm_lambda[[b]]), ncol = nrow(mm_lambda[[b]]))
      mb[-lavpta$vidx$lv.interaction[[b]], ] <- tmp
    }

    # handle observed-only variables (needed?)
    dummy_ov_idx <- c(
      fit@Model@ov.x.dummy.ov.idx[[b]],
      fit@Model@ov.y.dummy.ov.idx[[b]]
    )
    dummy_lv_idx <- c(
      fit@Model@ov.x.dummy.lv.idx[[b]],
      fit@Model@ov.y.dummy.lv.idx[[b]]
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
    if (lv_higherorder_flag) {
      dummy_lv_idx <- match(lv_names_b[dummy_lv_idx], lv_names1)
    }

    if (length(dummy_ov_idx)) {
      mb[dummy_lv_idx, ] <- 0
      mb[cbind(dummy_lv_idx, dummy_ov_idx)] <- 1
    }

    # here, we remove the lv.interaction row(s) from Mb
    # FIXME: if we have higher order factors!
    if (length(lavpta$vidx$lv.interaction[[b]]) > 0L) {
      mb <- mb[-lavpta$vidx$lv.interaction[[b]], ]
    }

    # compute EETA
    if (lavoptions$meanstructure) {
      if (lsam_analytic_flag[b]) {
        eeta[[b]] <- lav_sam_eeta(m = mb, ybar = ybar, mm_nu = mm_nu[[b]])
      } else {
        # EETA is constrained somehow
        lav_msg_stop(gettext("not ready yet"))
        # EETA[[b]] <- lav_sam_eeta_con(YBAR = YBAR, LAMBDA = LAMBDA[[b]],
        #                               THETA = THETA[[b]],
        #                               L.veta = L.veta[[b]])

      }
    fs_mean[[b]] <- eeta[[b]] # ok if no interaction
    }

    # compute VETA
    if (sam_method == "local") {
      if (lsam_analytic_flag[b]) {
        tmp <- lav_sam_veta(
          m = mb, s = cov_1, mm_theta = mm_theta[[b]],
          alpha_correction = local_options[["alpha.correction"]],
          lambda_correction = local_options[["lambda.correction"]],
          n <- fit@SampleStats@nobs[[this_group]],
          dummy_lv_idx = dummy_lv_idx,
          extra = TRUE
        )
        veta[[b]] <- tmp[, , drop = FALSE] # drop attributes
        alpha[[b]]  <- attr(tmp, "alpha")
        lambda[[b]] <- attr(tmp, "lambda.star")
      msm_[[b]]  <- attr(tmp, "MSM")
      mtm_[[b]]  <- attr(tmp, "MTM")
      } else {
        # VETA is constrained somehow
        veta[[b]] <- lav_sam_veta_con(s = cov_1, mm_lambda = mm_lambda[[b]],
                                      mm_theta = mm_theta[[b]],
                                      l_veta = l_veta[[b]],
                                      local_m_method = local_m_method)
        alpha[[b]]  <- as.numeric(NA)
        lambda[[b]] <- as.numeric(NA)
        msm_[[b]]  <- matrix(0, 0, 0)
        mtm_[[b]]  <- matrix(0, 0, 0)
      }
    } else if (sam_method == "cfsr") {
      # first, we need to 'true' VETA (to get Sigma)
      tmp <- lav_sam_veta(
        m = mb, s = cov_1, mm_theta = mm_theta[[b]],
        alpha_correction = 0L,
        lambda_correction = local_options[["lambda.correction"]],
        n <- fit@SampleStats@nobs[[this_group]],
        dummy_lv_idx = dummy_lv_idx,
        extra = FALSE
      )
      veta[[b]] <- tmp[, , drop = FALSE]
      # compute 'Sigma'
      sigma_1 <- this_lambda %*% veta[[b]] %*% t(this_lambda) + mm_theta[[b]]
      tmat <- lav_predict_tmat_det_internal(Sigma = sigma_1, Veta = veta[[b]],
                                            Lambda = this_lambda)
      a <- tmat %*% mb
      veta[[b]] <- a %*% cov_1 %*% t(a)
    } else {
      # FSR -- no correction
      veta[[b]] <- mb %*% cov_1 %*% t(mb)
    }

    # standardize? not really needed, but we may have 1.0000001
    # as variances, and this may lead to false convergence
    if (fit@Options$std.lv) {
      # warning: we should only do this for the LVs, not the
      # observed variables
      if (length(dummy_lv_idx) == 0L) {
        veta[[b]] <- stats::cov2cor(veta[[b]])
      } else {
        tmp <- veta[[b]]
        tmp_lv <- stats::cov2cor(veta[[b]][-dummy_lv_idx,
                                           -dummy_lv_idx, drop = FALSE])
        veta[[b]][-dummy_lv_idx, -dummy_lv_idx] <- tmp_lv
      }
    }
    colnames(veta[[b]]) <- rownames(veta[[b]]) <- lv_names1

    # compute model-based RELiability
    # REL[[b]] <- diag(VETA[[b]]] %*% solve(MSM..[[b]]))
    #                 CHECKme! -> done, must be:
    if (lsam_analytic_flag[b]) {
      rel[[b]] <- diag(veta[[b]]) / diag(msm_[[b]]) #!
    } else {
      rel[[b]] <- as.numeric(NA)
    }

    # check for lv.interactions
    if (lv_interaction_flag && length(lv_int_names) > 0L) {
      if (fit@Model@categorical || fit@Model@correlation) {
        lav_msg_stop(gettext("SAM + lv interactions do not work (yet) if
                             correlation structures are used."))
      }

      # compute Bartlett factor scores here
      if (is.null(y)) {
        yb <- fit@Data@X[[b]]
      } else {
        yb <- y[[b]]
      }
      # center
      yb_c <- t(t(yb) - drop(mm_nu[[b]]))
      fs_b <- yb_c %*% t(mb)
      colnames(fs_b) <- lv_names1
      # FIXME: what about observed covariates?

      # get (approximate) indices with outliers
      #fs.outlier.idx[[b]] <-
      #        lav_sample_outlier_idx(lav_sample_mdist(FS.b), coef = 1.5)

      # EETA2
      eeta1 <- eeta[[b]]
      eeta[[b]] <- lav_sam_eeta2(
        eeta = eeta1, veta = veta[[b]],
        lv_names = lv_names1,
        lv_int_names = lv_int_names
      )

      # VETA2
      if (sam_method == "local") {
        tmp <- lav_sam_veta2(
          fs = fs_b, m = mb,
          veta = veta[[b]], eeta = eeta1,
          mm_theta = mm_theta[[b]],
          lv_names = lv_names1,
          lv_int_names = lv_int_names,
          dummy_lv_names = lv_names_b[dummy_lv_idx],
          alpha_correction = local_options[["alpha.correction"]],
          lambda_correction = local_options[["lambda.correction"]],
          #fs.outlier.idx = fs.outlier.idx[[b]],
          return_fs = return_fs,
          return_cov_iveta2 = return_cov_iveta2,
          extra = TRUE
        )
        veta[[b]] <- tmp[, , drop = FALSE] # drop attributes
        alpha[[b]] <- attr(tmp, "alpha")
        lambda[[b]] <- attr(tmp, "lambda.star")
        msm_[[b]] <- attr(tmp, "MSM")
        mtm_[[b]] <- attr(tmp, "MTM")
    fs_mean[[b]] <- attr(tmp, "FS.mean")
        if (return_fs) {
          fs[[b]] <- attr(tmp, "FS")
        }
        if (return_cov_iveta2) {
          cov_iveta2[[b]] <- attr(tmp, "cov.iveta2")
        }
        #FS.gamma[[b]] <- attr(tmp, "FS.gamma")
      } else {
        lav_msg_fixme("not ready yet!")
        # FSR -- no correction
        veta[[b]] <- lav_sam_fs2(
          fs = fs_b,
          lv_names = lv_names1, lv_int_names = lv_int_names
        )
      }
    }

    # store Mapping matrix for this block
    m[[b]] <- mb
  } # blocks

  # label blocks
  if (nblocks > 1L) {
    names(eeta)       <- fit@Data@block.label
    names(veta)       <- fit@Data@block.label
    names(rel)        <- fit@Data@block.label
  names(msm_)      <- fit@Data@block.label
  names(mtm_)      <- fit@Data@block.label
  names(fs_mean)    <- fit@Data@block.label
    names(fs)         <- fit@Data@block.label
    names(cov_iveta2) <- fit@Data@block.label
    #names(fs.outlier.idx) <- FIT@Data@block.label
    #names(FS.gamma) <- FIT@Data@block.label
  }

  # handle conditional.x: add res.slopes, cov.x and mean.x
  if (fit@Model@conditional.x) {
    res_slopes <- vector("list", length = nblocks)
    for (b in seq_len(nblocks)) {
      res_slopes[[b]] <- m[[b]] %*% fit@h1$implied$res.slopes[[b]]
    }
    attr(veta, "res.slopes") <- res_slopes
    attr(veta, "cov.x") <- fit@h1$implied$cov.x
    attr(veta, "mean.x") <- fit@h1$implied$mean.x
  }

  # store EETA/VETA/M/alpha/lambda
  step1$VETA     <- veta
  step1$EETA     <- eeta
  step1$REL      <- rel
  step1$M        <- m
  step1$lambda   <- lambda
  step1$alpha    <- alpha
  step1$MSM      <- msm_
  step1$MTM      <- mtm_
  step1$FS.mean  <- fs_mean
  step1$FS       <- fs
  step1$COV.IVETA2 <- cov_iveta2
  #STEP1$fs.outlier.idx <- fs.outlier.idx
  #STEP1$FS.gamma <- FS.gamma
  step1$LV.NAMES <- lv_names_1
  # store also sam.method and local.options
  step1$sam.method <- sam_method
  step1$local.options <- local_options

  if (lav_verbose()) {
    cat("done.\n")
  }

  step1
}


lav_sam_step1_local_jac <- function(step1 = NULL, fit = NULL, p_only = FALSE,
                                    return_jac = FALSE) {

  lavdata <- fit@Data
  # lavsamplestats <- fit@SampleStats
  lavmodel <- fit@Model
  lavpta <- fit@pta
  nblocks <- lavpta$nblocks

  # local_options <- step1$local.options
  # sam_method <- step1$sam.method

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
    lav_msg_stop(gettext(
      "IJ local SEs: not available for the categorical setting (yet)!\n"))
  }
  n_mmblocks <- length(step1$MM.FIT)

  # JAC = (JACc %*% JACa) + JACb
  # - rows are the elements of vech(VETA)
  # - cols are the elements of vech(S)

  # JACa: mm.theta x vech(S)
  # JACc: vech(VETA) x mm.theta (keeping S fixed)
  # JACb: vech(VETA) x vech(S)  (keeping mm.theta fixed)

  # JACa: jacobian of theta.mm = f(vech(S))
  jaca <- matrix(0, nrow = length(fit@ParTable$lhs), # we select later
                    ncol = length(fit@SampleStats@WLS.obs[[g]]))
  for (mm in seq_len(n_mmblocks)) {
    fit_mm_block <- step1$MM.FIT[[mm]]
    mm_h1_expected   <- lavTech(fit_mm_block, "h1.information.expected")
    mm_delta         <- lavTech(fit_mm_block, "Delta")
    if (p_only && fit@Options$information[1] == "expected") {
      # for twostep.robust
      mm_inv_observed <- lavTech(fit_mm_block, "inverted.information.expected")
    } else {
      mm_inv_observed  <- lavTech(fit_mm_block, "inverted.information.observed")
    }

    #h1.info <- matrix(0, nrow(mm.h1.expected[[1]]), ncol(mm.h1.expected[[1]]))
    #for (g in seq_len(ngroups)) {
    #  fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    #  tmp <- fg * (mm.h1.expected[[g]] %*% mm.delta[[g]])
    #  h1.info <- h1.info + tmp
    #}
    mm_jac <- t(mm_h1_expected[[g]] %*% mm_delta[[g]] %*% mm_inv_observed)
    # keep only rows that are also in FIT@ParTable
    mm_keep_idx <- fit_mm_block@ParTable$free[step1$block.ptm.idx[[mm]]]
    mm_jac <- mm_jac[mm_keep_idx, , drop = FALSE]

    # select 'S' elements (row index)
    mm_ov_idx <- match(step1$MM.FIT[[mm]]@Data@ov.names[[g]],
                       lavdata@ov.names[[g]])
    mm_nvar <- length(lavdata@ov.names[[g]])
    mm_col_idx <- lav_matrix_vech_which_idx(mm_nvar, idx = mm_ov_idx,
      add_idx_at_start = lavmodel@meanstructure)
    mm_row_idx <- step1$block.mm.idx[[mm]][step1$block.ptm.idx[[mm]]]
    jaca[mm_row_idx, mm_col_idx] <- mm_jac
  }

  # keep only 'LAMBDA/THETA' parameters
  pt_1 <- step1$PT
  # only ov.names that are actually used in the measurement models
  ov_names <- unique(unlist(lapply(step1$MM.FIT, lav_object_vnames, "ov")))
  lambda_idx <- which(pt_1$op == "=~" & pt_1$free > 0L & !duplicated(pt_1$free))
  theta_idx  <- which(pt_1$op == "~~" & pt_1$free > 0L &
                      !duplicated(pt_1$free) &
                      pt_1$lhs %in% ov_names & pt_1$rhs %in% ov_names)
  nu_idx <- integer(0L)
  if (lavmodel@meanstructure) {
    nu_idx     <- which(pt_1$op == "~1" & pt_1$free > 0L &
                        !duplicated(pt_1$free) &
                        pt_1$lhs %in% ov_names)
  }
  delta_idx <- integer(0L)
  if (lavmodel@categorical || lavmodel@correlation) {
    delta_idx <- which(pt_1$op == "~*~" & pt_1$free > 0L &
                                   !duplicated(pt_1$free))
  }
  beta_idx <- psi_idx <- integer(0L)
  lv_ind <- unlist(lavpta$vnames$lv.ind)
  if (length(lv_ind) > 0L) {
    beta_idx <- which(pt_1$op == "=~" & pt_1$free > 0L &
                      !duplicated(pt_1$free) & pt_1$rhs %in% lv_ind)
    psi_idx <- which(pt_1$op == "~~" & pt_1$free > 0L &
                    !duplicated(pt_1$free) &
                    pt_1$rhs %in% lv_ind & pt_1$lhs %in% lv_ind)
  }
  # keep only these free parameters (measurement only)
  keep_idx <- sort(c(lambda_idx, theta_idx, nu_idx,
                     delta_idx, beta_idx, psi_idx))
  jaca <- jaca[keep_idx, , drop = FALSE]
  if (p_only) {
    return(jaca)
  }

  # JACb: jacobian of the function vech(VETA) = f(vech(S), theta.mm)
  #       (treating theta.mm as fixed)
  if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
    ffb <- function(x) {
      if (lavmodel@meanstructure) {
        nvar <- nrow(fit@h1$implied$cov[[g]])
        this_ybar <- x[seq_len(nvar)]
        this_cov <- lav_matrix_vech_reverse(x[-seq_len(nvar)])
      } else {
        this_ybar <- fit@h1$implied$mean[[g]]
        this_cov <- lav_matrix_vech_reverse(x)
      }

      # change COV/YBAR
      step1_1 <- step1
      step1_1$COV[[1]] <- this_cov
      if (lavmodel@meanstructure) {
        step1_1$YBAR[[1]] <- this_ybar
      }

      # transform data to comply with the new COV/YBAR
      y <- fit@Data@X[[1]]
      ytrans <- vector("list", nblocks)
      ytrans[[1]] <- lav_matrix_transform_mean_cov(y, target_mean = this_ybar,
                                                   target_cov = this_cov)
      colnames(ytrans[[1]]) <- fit@pta$vnames$ov[[1]]

      step1_1 <- lav_sam_step1_local(step1 = step1_1, fit = fit, y = ytrans,
           sam_method = step1$sam.method, local_options = step1$local.options)
      if (lavmodel@meanstructure) {
        out <- c(step1_1$EETA[[1]], lav_matrix_vech(step1_1$VETA[[1]]))
      } else {
        out <- lav_matrix_vech(step1_1$VETA[[1]])
      }
      out
    }
    # shut off verbose
    verbose_flag <- lav_verbose()
    lav_verbose(FALSE)
    this_x <- lav_matrix_vech(fit@h1$implied$cov[[g]])
    if (lavmodel@meanstructure) {
      this_x <- c(fit@h1$implied$mean[[g]], this_x)
    }
    jacb <- numDeriv::jacobian(func = ffb, x = this_x)
    lav_verbose(verbose_flag)
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
    mb <- step1$M[[g]]
    mbx_mb <- mb %x% mb
    row_idx <- lav_matrix_vech_idx(nrow(mb))
    jacb <- lav_matrix_duplication_post(mbx_mb)[row_idx, , drop = FALSE]
    if (lavmodel@meanstructure) {
      jacb <- lav_matrix_bdiag(mb, jacb)
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
    step1_1 <- step1

    # fill in new 'x' values in PT
    step1_1$PT$est[keep_idx] <- x
    step1_1 <- lav_sam_step1_local(step1 = step1_1, fit = fit,
           sam_method = step1$sam.method, local_options = step1$local.options)
    if (lavmodel@meanstructure) {
      out <- c(step1_1$EETA[[1]], lav_matrix_vech(step1_1$VETA[[1]]))
    } else {
      out <- lav_matrix_vech(step1_1$VETA[[1]])
    }
    out
  }
  # shut off verbose
  verbose_flag <- lav_verbose()
  lav_verbose(FALSE)
  jacc <- numDeriv::jacobian(func = ffc, x = step1$PT$est[keep_idx])
  lav_verbose(verbose_flag)

  # assemble JAC
  jac <- (jacc %*% jaca) + jacb

  if (return_jac) {
    attr(jac, "JACa") <- jaca
    attr(jac, "JACb") <- jacb
    attr(jac, "JACc") <- jacc
  }

  # eventually, this will be a list per group/block
  list(jac)
}

# semi-analytic version
# YR 4 June 2025: works, but still needs cleanup + avoid redundant calculations
# (eg when multiplied with a matrix with many zero rows/cols)
lav_sam_gamma_add <- function(step1 = NULL, fit = NULL, group = 1L) {

  lavdata <- fit@Data
  # lavsamplestats <- fit@SampleStats
  lavmodel <- fit@Model
  lavpta <- fit@pta
  # nblocks <- lavpta$nblocks

  # local_options <- step1$local.options
  # sam_method <- step1$sam.method

  ngroups <- lavdata@ngroups
  if (ngroups > 1L) {
    lav_msg_stop(gettext("IJ local SEs: not available with multiple groups!\n"))
  }
  g <- group
  y <- fit@Data@X[[g]]
  n <- nrow(y)
  p <- ncol(y)

  # NAMES + lv.keep
  lv_names <- step1$LV.NAMES[[1]]
  lv_names <- c("..int..", lv_names)
  nfac <- length(lv_names)
  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)

  # k_nfac <- lav_matrix_commutation(nfac, nfac)
  # ik <- diag(nfac * nfac) + k_nfac
  m_d <- lav_matrix_duplication(nfac)

  names_1 <- paste(lv_names[idx1], lv_names[idx2], sep = ":")
  names_1[seq_len(nfac)] <- lv_names
  lv_keep <- colnames(step1$VETA[[1]])
  fs_mean <- step1$FS.mean[[1]]
  keep_idx <- match(lv_keep, names_1)

  # step 1 free parameters
  step1_idx <- which(step1$PT$free %in% step1$step1.free.idx)
  x_step1 <- step1$PT$est[step1_idx]

  this_nu <- step1$NU[[1]]
  this_m  <- step1$M[[1]]
  this_mtm <- lav_matrix_bdiag(0,
                  step1$MTM[[1]][seq_len(nfac - 1), seq_len(nfac - 1)])
  this <- c(this_nu, lav_matrix_vec(this_m), lav_matrix_vech(this_mtm))
  index <- matrix(seq_len(nfac^2 * nfac^2), nfac^2, nfac^2)

  x2this <- function(x) {
    pt_1 <- step1$PT
   pt_1$est[step1_idx] <- x
    this_lavmodel <- lav_model_set_parameters(lavmodel,
                       x = pt_1$est[pt_1$free > 0 & !duplicated(pt_1$free)])
    this_nu     <- this_lavmodel@GLIST$nu
    rm_idx <- lavpta$vidx$lv.interaction[[1]]
    # no interaction columns!
    this_lambda <- this_lavmodel@GLIST$lambda[, -rm_idx, drop = FALSE]
    this_theta  <- this_lavmodel@GLIST$theta
    this_m <- lav_sam_mapping_matrix(mm_lambda = this_lambda,
                                     mm_theta = this_theta,
                                     s = step1$COV[[1]],
                                     method = step1$local.options$M.method)
    mtm <- this_m %*% this_theta %*% t(this_m)
    mtm <- lav_matrix_bdiag(0, mtm)

    out <- c(this_nu, lav_matrix_vec(this_m), lav_matrix_vech(mtm))
    out
  }
  jac_x2this <- numDeriv::jacobian(func = x2this, x = x_step1)
  # 46 x 24

  # JAC.eetai2.this
  jac_eeta2_this <- matrix(0, nrow = length(step1$EETA[[1]]),
                              ncol = length(this))
  pstar <- nfac * (nfac + 1) / 2
  idx <- p + length(this_m) + seq_len(pstar)
  jac_eeta2_this[, idx] <- (-diag(nfac^2) %*% m_d)[keep_idx, ]

  # define function for JAC.veta2.fi
  get_jac_veta2_fi <- function(x, mtm, fs_mean = NULL,
                      lambda_star = 1, keep_idx = NULL) {
    fi <- drop(fi)
    p <- length(fi)
    stopifnot(p == ncol(mtm))

    #INDEX <- matrix(seq_len(P^2 * P^2), P^2, P^2)

    # (tcrossprod(fi) - MTM) %x% MTM) -- part A
    block_a <- fi %x% do.call("rbind",
            lapply(seq_len(p), function(i) diag(p) %x% mtm[, i]))
    long_a <- diag(p) %x% lav_matrix_vec(fi %x% mtm)
    big_a <- block_a + long_a
    part_a <- big_a[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # (MTM %x% (tcrossprod(fi) - MTM)) -- part B
    block_b <- lav_matrix_vec(fi %x% mtm) %x% diag(p)
    long_b <- do.call("rbind",
            lapply(seq_len(p), function(i) diag(p) %x% mtm[, i])) %x% fi
    big_b <- block_b + long_b
    part_b <- big_b[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # lav_matrix_commutation_post((tcrossprod(fi) - MTM) %x% MTM) -- part C
    block_c <- do.call("rbind",
            lapply(seq_len(p), function(i) fi %x% (diag(p) %x% mtm[, i])))
    long_c  <- do.call("rbind",
            lapply(seq_len(p), function(i) diag(p) %x% (fi %x% mtm[, i])))
    big_c <- block_c + long_c
    part_c <- big_c[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # lav_matrix_commutation_pre((tcrossprod(fi) - MTM) %x% MTM) -- part m_d
    long_d <- diag(p) %x% lav_matrix_vec(mtm %x% fi)
    block_d <- (fi %x% lav_matrix_vec(mtm)) %x% diag(p)
    big_d <- block_d + long_d
    part_d <- big_d[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # t1 (tcrossprod(fi2[keep.idx] - FS.mean) only) -- part E
    idx1 <- rep(seq_len(p), each = p)
    idx2 <- rep(seq_len(p), times = p)
    fi2 <- (fi[idx1] * fi[idx2])
    fik <- fi2[keep_idx]

    fi <- as.matrix(fi)
    tt <- ((fi %x% diag(p))[keep_idx, ] + (diag(p) %x% fi)[keep_idx, ])
    tmp <- ((fik - fs_mean) %x% tt) + (tt %x% (fik - fs_mean))
    part_e <- tmp[lav_matrix_vech_idx(length(fik)), ]

    final <- part_e - lambda_star * (part_a + part_b + part_c + part_d)
    final
  }

  # define function for JAC.veta2.this.i
  get_jac_veta2_this <- function(x, fi = NULL,
                        lambda_star = 1, keep_idx = NULL) {
    nvar <- ncol(step1$M[[1]])
    nfac <- nrow(step1$M[[1]])
    # this_nu <- x[seq_len(nvar)]
    x <- x[-seq_len(nvar)]
    # this_m  <- matrix(x[seq_len(nvar * nfac)], nrow = nfac, ncol = nvar)
    x <- x[-seq_len(nvar * nfac)]
    mtm <- lav_matrix_vech_reverse(x)
    fi <- as.matrix(fi)
    p <- nrow(fi)

    # part a: ((tcrossprod(fi) - MTM) %x% MTM)
    long_a <- diag(p) %x% do.call("rbind",
                lapply(seq_len(p), function(i) diag(p) %x% -mtm[, i]))
    block_a <- do.call("rbind",
                lapply(seq_len(p), function(i) {
                  diag(p) %x% (tcrossprod(fi) - mtm)[, i]
                        })) %x% diag(p)
    big_a_vec <- long_a + block_a
    big_a <- big_a_vec %*% m_d
    part_a <- big_a[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # part b: (MTM %x% (tcrossprod(fi) - MTM))
    block_b <- do.call("rbind",
             lapply(seq_len(p), function(i) diag(p) %x% -mtm[, i])) %x% diag(p)
    long_b <- diag(p) %x% do.call("rbind",
             lapply(seq_len(p), function(i) {
               diag(p) %x% (tcrossprod(fi) - mtm)[, i]
                        }))
    big_b_vec <- block_b + long_b
    big_b <- big_b_vec %*% m_d
    part_b <- big_b[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # part c: lav_matrix_commutation_post((tcrossprod(fi) - MTM) %x% MTM)
    block_c <- diag(p) %x% do.call("rbind",
           lapply(seq_len(p),
                       function(i) ((tcrossprod(fi) - mtm)[, i]) %x% diag(p)))
    long_c <- do.call("rbind",
           lapply(seq_len(p), function(i) diag(p * p) %x% (-mtm[, i])))
    big_c_vec <- block_c + long_c
    big_c <- big_c_vec %*% m_d
    part_c <- big_c[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # part d: lav_matrix_commutation_pre((tcrossprod(fi) - MTM) %x% MTM)
    block_d <- diag(p) %x% do.call("rbind",
                  lapply(seq_len(p), function(i) -mtm[, i] %x% diag(p)))
    long_d <- do.call("rbind",
                  lapply(seq_len(p),
                     function(i) diag(p * p) %x% ((tcrossprod(fi) - mtm)[, i])))
    big_d_vec <- block_d + long_d
    big_d <- big_d_vec %*% m_d
    part_d <- big_d[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    # part e: (IK %*% (MTM %x% MTM))
    m_k <- lav_matrix_commutation(p, p)
    e1 <- diag(p) %x% do.call("rbind",
            lapply(seq_len(p), function(i) diag(p) %x% mtm[, i]))
    e2 <- do.call("rbind",
            lapply(seq_len(p), function(i) diag(p * p) %x% mtm[, i]))
    e3 <- do.call("rbind",
            lapply(seq_len(p), function(i) diag(p) %x% mtm[, i])) %x% diag(p)
    e4 <- diag(p) %x% do.call("rbind",
            lapply(seq_len(p), function(i) m_k %*% (diag(p) %x% mtm[, i])))
    big_e_vec <- e1 + e2 + e3 + e4
    big_e <- big_e_vec %*% m_d
    part_e <-  big_e[lav_matrix_vech(index[keep_idx, keep_idx]), ]

    final <- -1 * lambda_star * (part_a + part_b + part_c + part_d + part_e)
    final
                }

  n_eeta_veta <- length(step1$EETA[[1]]) +
                 length(lav_matrix_vech(step1$VETA[[1]]))
  cveta <- matrix(0, nrow = n_eeta_veta, ncol = length(x_step1))
  for (i in 1:n) {
    # experimental: remove fs outliers
    #if (i %in% STEP1$fs.outlier.idx[[1]]) {
    #  next
    #}

    # factor score
    fi <- rbind(1, this_m %*% (y[i, ] - this_nu))

    #tmp <- numDeriv::jacobian(func = theta.to.eetavetai, x = x.step1, i = i)
    jac_this2fi_i <- rbind(0, cbind(-1 * this_m,
                           t(as.matrix(y[i, ] - drop(this_nu))) %x%
                                                         diag(nrow(this_m)),
                           matrix(0, nrow = nrow(this_m),
                                  ncol = length(lav_matrix_vech(this_mtm)))))

    jac_eeta2_fi_i <- ((fi %x% diag(nfac)) + (diag(nfac) %x% fi))[keep_idx, ]
    eeta2 <- ((jac_eeta2_fi_i %*% jac_this2fi_i) +
                                              jac_eeta2_this) %*% jac_x2this

    jac_veta2_fi_i <- get_jac_veta2_fi(x = fi, mtm = this_mtm,
                                       fs_mean = fs_mean,
                                       lambda_star = step1$lambda[[1]],
                                       keep_idx = keep_idx)
    jac_veta2_this_i <- matrix(0, nrow = nrow(jac_veta2_fi_i),
                                                  ncol = ncol(jac_this2fi_i))
    jac_veta2_this_i[, idx] <- get_jac_veta2_this(x = this, fi = fi,
                                                  lambda_star = step1$lambda[[1]],
                                                  keep_idx = keep_idx)
    veta2 <- ((jac_veta2_fi_i %*% jac_this2fi_i) +
                                       jac_veta2_this_i) %*% jac_x2this

    # FIXME: remove outliers somehow?
    cveta <- cveta + 1 / n * rbind(eeta2, veta2)
  }

  gamma_addition <- n * (cveta %*% step1$Sigma.11 %*% t(cveta))
  gamma_addition
}
