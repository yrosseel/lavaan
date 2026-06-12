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
      veta_symbolic_1 <- lav_sam_veta_pt(fit, block = b)
      # this is the tricky thing: which rows/cols should we remove?
      # none for now
      if (fit@Options$std.lv) {
        veta_symbolic <- lav_mat_vech(veta_symbolic_1, diagonal = FALSE)
      } else {
        veta_symbolic <- lav_mat_vech(veta_symbolic_1, diagonal = TRUE)
      }
      if (any(veta_symbolic == 0)) {
        lsam_analytic_flag[b] <- FALSE
        nfac <- ncol(veta_symbolic_1)
        m_free <- which(veta_symbolic_1 != 0)
        tmp <- lav_mat_vech_rev(lav_mat_vech_idx(nfac))
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

  veta       <- vector("list", nblocks)
  msm_list   <- vector("list", nblocks)
  mtm_list   <- vector("list", nblocks)
  fs_mean    <- vector("list", nblocks)
  fs         <- vector("list", nblocks)
  cov_iveta2 <- vector("list", nblocks)
  rel        <- vector("list", nblocks)
  alpha      <- vector("list", nblocks)
  lambda     <- vector("list", nblocks)
  lambda1    <- vector("list", nblocks)
  if (lavoptions$meanstructure) {
    eeta <- vector("list", nblocks)
  } else {
    eeta <- NULL
  }
  m <- vector("list", nblocks)
  lv_names_1 <- vector("list", nblocks)

  # compute VETA/EETA per block
  for (b in seq_len(nblocks)) {

    # which group is this?
    this_group <- floor(b / fit@Data@nlevels + 0.5)

    # lv.names, including dummy-lv covariates
    psi_idx <- which(names(fit@Model@GLIST) == "psi")[b]
    lv_names_b <- fit@Model@dimNames[[psi_idx]][[1L]] # including dummy/inter.
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
      mb <- lav_sam_mapping_mat(
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
          n = fit@SampleStats@nobs[[this_group]],
          dummy_lv_idx = dummy_lv_idx,
          extra = TRUE
        )
        veta[[b]] <- tmp[, , drop = FALSE] # drop attributes
        alpha[[b]]  <- attr(tmp, "alpha")
        lambda[[b]] <- attr(tmp, "lambda.star")
        msm_list[[b]] <- attr(tmp, "MSM")
        mtm_list[[b]] <- attr(tmp, "MTM")
        # effective first-order coefficient on the (pd-forced) MTM matrix,
        # so that VETA = MSM - lambda1 * MTM; this absorbs both the
        # lambda.star correction (1 if the correction was not needed) and
        # the alpha correction (lav_sam_veta() rescales MTM by
        # (1 - alpha/(n-1)) *before* the lambda step); when lv interactions
        # are present, lambda[[b]] will be overwritten by the second-order
        # lambda.star, but the casewise decomposition of EETA2/VETA2 (see
        # lav_sam_veta2() and lav_sam_gamma_add()) needs this first-order
        # value
        nobs_b <- fit@SampleStats@nobs[[this_group]]
        alpha_n1 <- lav_sam_alpha_n1(local_options[["alpha.correction"]],
                                     nobs_b)
        l1 <- if (is.finite(lambda[[b]])) lambda[[b]] else 1
        lambda1[[b]] <- (1 - alpha_n1) * l1
      } else {
        # VETA is constrained somehow
        veta[[b]] <- lav_sam_veta_con(s = cov_1, mm_lambda = mm_lambda[[b]],
                                      mm_theta = mm_theta[[b]],
                                      l_veta = l_veta[[b]],
                                      local_m_method = local_m_method)
        alpha[[b]]   <- as.numeric(NA)
        lambda[[b]]  <- as.numeric(NA)
        lambda1[[b]] <- as.numeric(NA)
        msm_list[[b]] <- matrix(0, 0, 0)
        mtm_list[[b]] <- matrix(0, 0, 0)
      }
    } else if (sam_method == "cfsr") {
      # first, we need the 'true' VETA (to get Sigma)
      tmp <- lav_sam_veta(
        m = mb, s = cov_1, mm_theta = mm_theta[[b]],
        alpha_correction = 0L,
        lambda_correction = local_options[["lambda.correction"]],
        n = fit@SampleStats@nobs[[this_group]],
        dummy_lv_idx = dummy_lv_idx,
        extra = FALSE
      )
      veta[[b]] <- tmp[, , drop = FALSE]
      # compute 'Sigma'
      sigma_1 <- this_lambda %*% veta[[b]] %*% t(this_lambda) + mm_theta[[b]]
      tmat <- lav_predict_tmat_det_internal(sigma_1 = sigma_1, veta = veta[[b]],
                                            lambda = this_lambda)
      a <- tmat %*% mb
      veta[[b]] <- a %*% cov_1 %*% t(a)
    } else {
      # FSR -- no correction
      veta[[b]] <- mb %*% cov_1 %*% t(mb)
    }

    # standardize? not really needed, but we may have 1.0000001
    # as variances, and this may lead to false convergence
    # d_std holds the rescaling factors (1 = no rescaling); when lv
    # interactions are present, we will need them to express the factor
    # scores (and related matrices) in the same standardized metric
    d_std <- rep(1, ncol(veta[[b]]))
    if (fit@Options$std.lv) {
      # we should only standardize the LVs, not the observed variables
      # (dummy lvs); but note that the cross-blocks between the dummy lvs
      # and the latent variables must be rescaled by the factor-side
      # scaling factors as well (new in June 2026; before, these
      # cross-blocks were left unrescaled, and -- when interaction terms
      # were present -- the dummy lvs were accidentally standardized too,
      # because dummy_lv_idx refers to the full set of lv names, while
      # veta[[b]] only contains lv_names1)
      dummy_lv_idx1 <- dummy_lv_idx
      if (!lv_higherorder_flag) { # if higher-order, already remapped
        dummy_lv_idx1 <- match(lv_names_b[dummy_lv_idx], lv_names1)
      }
      if (length(dummy_lv_idx1) == 0L) {
        d_std <- sqrt(diag(veta[[b]]))
        veta[[b]] <- stats::cov2cor(veta[[b]])
      } else {
        d_std[-dummy_lv_idx1] <- sqrt(diag(veta[[b]]))[-dummy_lv_idx1]
        veta[[b]] <- t(veta[[b]] / d_std) / d_std
      }
    }
    colnames(veta[[b]]) <- rownames(veta[[b]]) <- lv_names1

    # compute model-based RELiability
    if (lsam_analytic_flag[b]) {
      rel[[b]] <- diag(veta[[b]]) / diag(msm_list[[b]])
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

      # if std.lv = TRUE, express the factor scores, EETA and the mapping
      # matrix in the same standardized metric as the cov2cor() rescaled
      # VETA above, so that the augmented statistics EETA2/VETA2 (and
      # their casewise decomposition in lav_sam_veta2()) are
      # metric-coherent
      eeta1 <- eeta[[b]]
      mb_v2 <- mb
      if (fit@Options$std.lv) {
        fs_b <- t(t(fs_b) / d_std)
        eeta1 <- eeta1 / d_std
        mb_v2 <- mb / d_std
      }

      # EETA2
      eeta[[b]] <- lav_sam_eeta2(
        eeta = eeta1, veta = veta[[b]],
        lv_names = lv_names1,
        lv_int_names = lv_int_names
      )

      # VETA2
      if (sam_method == "local") {
        tmp <- lav_sam_veta2(
          fs = fs_b, m = mb_v2,
          veta = veta[[b]], eeta = eeta1,
          mm_theta = mm_theta[[b]],
          lv_names = lv_names1,
          lv_int_names = lv_int_names,
          dummy_lv_names = lv_names_b[dummy_lv_idx],
          alpha_correction = local_options[["alpha.correction"]],
          lambda_correction = local_options[["lambda.correction"]],
          lambda1 = lambda1[[b]],
          return_fs = return_fs,
          return_cov_iveta2 = return_cov_iveta2,
          extra = TRUE
        )
        veta[[b]] <- tmp[, , drop = FALSE] # drop attributes
        alpha[[b]] <- attr(tmp, "alpha")
        lambda[[b]] <- attr(tmp, "lambda.star")
        msm_list[[b]] <- attr(tmp, "MSM")
        mtm_list[[b]] <- attr(tmp, "MTM")
        fs_mean[[b]] <- attr(tmp, "FS.mean")
        if (return_fs) {
          fs[[b]] <- attr(tmp, "FS")
        }
        if (return_cov_iveta2) {
          cov_iveta2[[b]] <- attr(tmp, "cov.iveta2")
        }
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
    if (!is.null(eeta)) {
      names(eeta)     <- fit@Data@block.label
    }
    names(veta)       <- fit@Data@block.label
    names(rel)        <- fit@Data@block.label
    names(msm_list)   <- fit@Data@block.label
    names(mtm_list)   <- fit@Data@block.label
    names(fs_mean)    <- fit@Data@block.label
    names(fs)         <- fit@Data@block.label
    names(cov_iveta2) <- fit@Data@block.label
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
  step1$lambda1  <- lambda1
  step1$alpha    <- alpha
  step1$MSM      <- msm_list
  step1$MTM      <- mtm_list
  step1$FS.mean  <- fs_mean
  step1$FS       <- fs
  step1$COV.IVETA2 <- cov_iveta2
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
  if (lavmodel@categorical && !p_only) {
    lav_msg_stop(gettext(
      "IJ local SEs: not available for the categorical setting (yet)!\n"))
  }
  if (lavmodel@categorical && lavmodel@conditional.x) {
    lav_msg_stop(gettext(
      "twostep.robust SEs: not available for the categorical setting in
       combination with conditional.x = TRUE (yet)!\n"))
  }
  n_mmblocks <- length(step1$MM.FIT)

  # categorical: we map the statistics of each measurement block into the
  # WLS statistics vector of the joint model by label
  if (lavmodel@categorical) {
    joint_labels <- lav_sam_wls_obs_labels(
      ov_names = lavdata@ov.names[[g]],
      th_idx = fit@SampleStats@th.idx[[g]]
    )
  }

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

    if (lavmodel@categorical) {
      if (!fit_mm_block@Model@categorical) {
        lav_msg_stop(gettext(
          "twostep.robust SEs: not available (yet) when categorical and
           continuous measurement blocks are mixed"))
      }
      # influence of the (D)WLS estimator: (Delta' W Delta)^-1 Delta' W
      # (we reuse the same ingredients that produced the robust vcov of
      # this measurement block)
      tmp <- lav_model_nvcov_robust_sem(
        lavmodel = fit_mm_block@Model,
        lavsamplestats = fit_mm_block@SampleStats,
        lavcache = fit_mm_block@cache, lavdata = fit_mm_block@Data,
        lavimplied = fit_mm_block@implied, lavh1 = fit_mm_block@h1,
        lavoptions = fit_mm_block@Options, use_ginv = FALSE,
        attr_delta = TRUE, attr_e_inv = TRUE, attr_wls_v = TRUE)
      mm_e_inv <- attr(tmp, "E.inv")
      mm_wls_v <- attr(tmp, "WLS.V")[[g]]
      mm_delta_1 <- attr(tmp, "Delta")[[g]]
      if (is.matrix(mm_wls_v)) {
        mm_jac <- mm_e_inv %*% t(mm_wls_v %*% mm_delta_1)
      } else {
        # DWLS/ULS: WLS.V holds the diagonal of the weight matrix
        mm_jac <- mm_e_inv %*% t(mm_wls_v * mm_delta_1)
      }
    } else {
      mm_h1_expected <- lavTech(fit_mm_block, "h1.information.expected")
      mm_delta       <- lavTech(fit_mm_block, "Delta")
      if (p_only && fit@Options$information[1] == "expected") {
        # for twostep.robust
        mm_inv_observed <-
          lavTech(fit_mm_block, "inverted.information.expected")
      } else {
        mm_inv_observed <-
          lavTech(fit_mm_block, "inverted.information.observed")
      }
      mm_jac <- t(mm_h1_expected[[g]] %*% mm_delta[[g]] %*% mm_inv_observed)
    }

    # keep only rows that are also in FIT@ParTable
    mm_keep_idx <- fit_mm_block@ParTable$free[step1$block.ptm.idx[[mm]]]
    mm_jac <- mm_jac[mm_keep_idx, , drop = FALSE]

    # map the columns of mm_jac (the sample statistics of this
    # measurement block) into the statistics vector of the joint model
    if (lavmodel@categorical) {
      mm_labels <- lav_sam_wls_obs_labels(
        ov_names = fit_mm_block@Data@ov.names[[g]],
        th_idx = fit_mm_block@SampleStats@th.idx[[g]]
      )
      mm_col_idx <- match(mm_labels, joint_labels)
      if (anyNA(mm_col_idx)) {
        lav_msg_stop(gettext(
          "internal error: unable to map the statistics of a measurement
           block into the statistics vector of the joint model"))
      }
    } else {
      mm_ov_idx <- match(step1$MM.FIT[[mm]]@Data@ov.names[[g]],
                         lavdata@ov.names[[g]])
      mm_nvar <- length(lavdata@ov.names[[g]])
      mm_col_idx <- lav_mat_vech_which_idx(mm_nvar, idx = mm_ov_idx,
        add_idx_at_start = lavmodel@meanstructure)
    }
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
  th_idx <- integer(0L)
  if (lavmodel@categorical) {
    th_idx <- which(pt_1$op == "|" & pt_1$free > 0L &
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
                     delta_idx, th_idx, beta_idx, psi_idx))
  jaca <- jaca[keep_idx, , drop = FALSE]
  if (p_only) {
    # the rows of P are in ascending free-parameter order; return the
    # free indices, so the caller can align them with step1.free.idx
    # (which is ordered per measurement block)
    attr(jaca, "free.idx") <- step1$PT.free[keep_idx]
    return(jaca)
  }

  # JACb: jacobian of the function vech(VETA) = f(vech(S), theta.mm)
  #       (treating theta.mm as fixed)
  if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
    ffb <- function(x) {
      if (lavmodel@meanstructure) {
        nvar <- nrow(fit@h1$implied$cov[[g]])
        this_ybar <- x[seq_len(nvar)]
        this_cov <- lav_mat_vech_rev(x[-seq_len(nvar)])
      } else {
        this_ybar <- fit@h1$implied$mean[[g]]
        this_cov <- lav_mat_vech_rev(x)
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
      ytrans[[1]] <- lav_mat_transform_mean_cov(y, target_mean = this_ybar,
                                                   target_cov = this_cov)
      colnames(ytrans[[1]]) <- fit@pta$vnames$ov[[1]]

      step1_1 <- lav_sam_step1_local(step1 = step1_1, fit = fit, y = ytrans,
           sam_method = step1$sam.method, local_options = step1$local.options)
      if (lavmodel@meanstructure) {
        out <- c(step1_1$EETA[[1]], lav_mat_vech(step1_1$VETA[[1]]))
      } else {
        out <- lav_mat_vech(step1_1$VETA[[1]])
      }
      out
    }
    # shut off verbose
    verbose_flag <- lav_verbose()
    lav_verbose(FALSE)
    this_x <- lav_mat_vech(fit@h1$implied$cov[[g]])
    if (lavmodel@meanstructure) {
      this_x <- c(fit@h1$implied$mean[[g]], this_x)
    }
    jacb <- numDeriv::jacobian(func = ffb, x = this_x)
    lav_verbose(verbose_flag)
  } else { # no latent interactions
    mb <- step1$M[[g]]
    mbx_mb <- mb %x% mb
    row_idx <- lav_mat_vech_idx(nrow(mb))
    jacb <- lav_mat_dup_post(mbx_mb)[row_idx, , drop = FALSE]
    if (lavmodel@meanstructure) {
      jacb <- lav_mat_bdiag(mb, jacb)
    }
  }

  # JACc: jacobian of the function vech(VETA) = f(theta.mm, vech(S))
  #       (treating vech(S) as fixed)
  # calling lav_sam_step1_local() directly
  ffc <- function(x) {
    step1_1 <- step1

    # fill in new 'x' values in PT
    step1_1$PT$est[keep_idx] <- x
    step1_1 <- lav_sam_step1_local(step1 = step1_1, fit = fit,
           sam_method = step1$sam.method, local_options = step1$local.options)
    if (lavmodel@meanstructure) {
      out <- c(step1_1$EETA[[1]], lav_mat_vech(step1_1$VETA[[1]]))
    } else {
      out <- lav_mat_vech(step1_1$VETA[[1]])
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
# YR 12 June 2026: rewritten to avoid the loop over the observations:
# because differentiation is linear, the average (over the observations) of
# the casewise jacobians equals the jacobian of the *average* casewise
# contribution, and the latter can be computed with a few (small) matrix
# operations; the centering terms (involving the mean of the second-order
# factor scores) can be ignored, as they cancel out when summing over the
# observations
lav_sam_gamma_add <- function(step1 = NULL, fit = NULL, group = 1L) {

  lavdata <- fit@Data
  lavmodel <- fit@Model
  lavpta <- fit@pta

  ngroups <- lavdata@ngroups
  if (ngroups > 1L) {
    lav_msg_stop(gettext("IJ local SEs: not available with multiple groups!\n"))
  }
  g <- group
  y <- fit@Data@X[[g]]
  n <- nrow(y)

  # NAMES + lv.keep
  lv_names <- c("..int..", step1$LV.NAMES[[1]])
  tmp <- lav_sam_lvnames2(lv_names)
  lv_keep <- colnames(step1$VETA[[1]])
  keep_idx <- match(lv_keep, tmp$names)

  # row/column pairs of the kept elements of (eta* %x% eta*)
  i1k <- tmp$idx1[keep_idx]
  i2k <- tmp$idx2[keep_idx]

  # the second-order lambda.star, and the effective first-order lambda1,
  # are treated as fixed constants (not differentiated)
  lambda_star <- step1$lambda[[1]]
  lambda1 <- step1$lambda1[[1]]
  if (is.null(lambda1) || !is.finite(lambda1)) {
    lambda1 <- 1
  }
  # effective second-order multiplier of the unscaled var.error term,
  # absorbing the alpha correction (see lav_sam_veta2())
  alpha_n1 <- lav_sam_alpha_n1(step1$alpha[[1]], n)
  lambda2_eff <- lambda_star * (1 - alpha_n1)

  # dummy lvs (observed variables in the structural part): mirror the
  # patches of lav_sam_step1_local() inside lbar() below; a dummy lv's
  # factor score is the observed variable itself (no measurement error),
  # and the patched rows/columns are *constants* (zero derivative); for
  # M.method = "ML" the patches are usually redundant (the model matrices
  # already encode lambda = 1/theta = 0, and the mapping matrix reproduces
  # the unit rows), but for M.method = "GLS" the unpatched mapping matrix
  # rows of the dummy lvs are *not* unit vectors
  lambda_gidx <- which(names(lavmodel@GLIST) == "lambda")[1]
  lv_names_full <- lavmodel@dimNames[[lambda_gidx]][[2L]]
  rm_idx <- lavpta$vidx$lv.interaction[[1]]
  lv_names_noint <- lv_names_full
  if (length(rm_idx) > 0L) {
    lv_names_noint <- lv_names_full[-rm_idx]
  }
  dummy_ov_y_idx <- lavmodel@ov.y.dummy.ov.idx[[1]]
  dummy_lv_y_idx <- match(lv_names_full[lavmodel@ov.y.dummy.lv.idx[[1]]],
                          lv_names_noint)
  dummy_ov_idx <- c(lavmodel@ov.x.dummy.ov.idx[[1]], dummy_ov_y_idx)
  dummy_lv_idx <- match(lv_names_full[c(lavmodel@ov.x.dummy.lv.idx[[1]],
                                        lavmodel@ov.y.dummy.lv.idx[[1]])],
                        lv_names_noint)

  # std.lv = TRUE? if so, lav_sam_step1_local() has rescaled VETA
  # (cov2cor) and expressed the factor scores in the same standardized
  # metric; we must mirror this inside lbar() below (the rescaling factors
  # depend on the step 1 parameters, so they are differentiated through)
  std_lv_flag <- fit@Options$std.lv

  # step 1 free parameters
  step1_idx <- which(step1$PT$free %in% step1$step1.free.idx)
  x_step1 <- step1$PT$est[step1_idx]
  pt_1 <- step1$PT

  # the average casewise contribution to the augmented summary statistics
  #   Lbar = [ E(eta* %x% eta*)[keep], vech(Var(eta* %x% eta*)[keep, keep]) ]
  # as a function of the step 1 free parameters (holding the data, and
  # lambda.star/lambda1, fixed); see the casewise contributions (iveta2_1)
  # in lav_sam_veta2(); the averages equal the EETA2/VETA2 statistics that
  # are actually used in step 2:
  # with f_i = (1, M (y_i - nu)), B = bdiag(0, M THETA t(M)),
  # C2 = (1/N) \sum_i f_i t(f_i), E = C2 - lambda1 * B, and writing
  # (r1, s1) and (r2, s2) for the row/column pairs of elements 'a' and 'b'
  # of (eta* %x% eta*), the averages of the casewise contributions are:
  # - first-order part:
  #     Lbar1[a] = (1/N) \sum_i (f_i %x% f_i)[a] - lambda1 * vec(B)[a]
  # - second-order part (before taking vech):
  #     Lbar2[a, b] = Var(fs2)[a, b] - lambda2.eff * tmpbar[a, b]
  #   with lambda2.eff = lambda.star * (1 - alpha/(n-1))
  #   where fs2_i = f_i %x% f_i, and tmpbar is the average of
  #     tmp_i = ((F_i - lambda1 * B) %x% B) + (B %x% (F_i - lambda1 * B)) +
  #             ((F_i - lambda1 * B) %x% B) %*% K +
  #             K %*% ((F_i - lambda1 * B) %x% B) +
  #             (I + K) %*% (B %x% B)              (with F_i = f_i t(f_i))
  #   which depends on the f_i's only through C2; elementwise:
  #     tmpbar[a, b] = E[r1, r2] B[s1, s2] + B[r1, r2] E[s1, s2] +
  #                    E[r1, s2] B[s1, r2] + E[s1, r2] B[r1, s2] +
  #                    B[r1, r2] B[s1, s2] + B[s1, r2] B[r1, s2]
  lbar <- function(x) {
    pt_x <- pt_1
    pt_x$est[step1_idx] <- x
    this_lavmodel <- lav_model_set_parameters(lavmodel,
                       x = pt_x$est[pt_x$free > 0 & !duplicated(pt_x$free)])
    this_nu     <- drop(this_lavmodel@GLIST$nu)
    # no interaction columns!
    this_lambda <- this_lavmodel@GLIST$lambda[, -rm_idx, drop = FALSE]
    this_theta  <- this_lavmodel@GLIST$theta
    # dummy lv patches (cfr. lav_sam_step1_local())
    if (length(dummy_ov_y_idx) > 0L) {
      this_theta[dummy_ov_y_idx, ] <- 0
      this_theta[, dummy_ov_y_idx] <- 0
      this_lambda[dummy_ov_y_idx, ] <- 0
      this_lambda[cbind(dummy_ov_y_idx, dummy_lv_y_idx)] <- 1
      this_nu[dummy_ov_y_idx] <- 0
    }
    this_m <- lav_sam_mapping_mat(mm_lambda = this_lambda,
                                     mm_theta = this_theta,
                                     s = step1$COV[[1]],
                                     method = step1$local.options$M.method)
    if (length(dummy_ov_idx) > 0L) {
      this_m[dummy_lv_idx, ] <- 0
      this_m[cbind(dummy_lv_idx, dummy_ov_idx)] <- 1
    }
    if (std_lv_flag) {
      # mirror the cov2cor() rescaling of the first-order VETA: rescale
      # the rows of the mapping matrix so that the (lambda1-corrected)
      # factor-score covariance matrix has unit diagonal; dummy lvs are
      # never rescaled
      msm_diag <- rowSums((this_m %*% step1$COV[[1]]) * this_m)
      mtm_diag <- rowSums((this_m %*% this_theta) * this_m)
      d_std <- sqrt(msm_diag - lambda1 * mtm_diag)
      d_std[dummy_lv_idx] <- 1
      this_m <- this_m / d_std
    }
    b_aug <- lav_mat_bdiag(0, this_m %*% this_theta %*% t(this_m))

    # (Bartlett) factor scores, augmented with an intercept
    fs <- cbind(1, t(t(y) - this_nu) %*% t(this_m))

    # C2 = (1/N) \sum_i f_i t(f_i)
    c2 <- crossprod(fs) / n
    e_mat <- c2 - lambda1 * b_aug

    # first-order part
    out1 <- e_mat[cbind(i2k, i1k)]

    # second-order part: only the kept columns of fs2 are needed
    fs2k <- fs[, i1k, drop = FALSE] * fs[, i2k, drop = FALSE]
    fs2kc <- t(t(fs2k) - colMeans(fs2k))
    var_fs2k <- crossprod(fs2kc) / n

    tmpbar <- (e_mat[i1k, i1k] * b_aug[i2k, i2k] +
               b_aug[i1k, i1k] * e_mat[i2k, i2k] +
               e_mat[i1k, i2k] * b_aug[i2k, i1k] +
               e_mat[i2k, i1k] * b_aug[i1k, i2k] +
               b_aug[i1k, i1k] * b_aug[i2k, i2k] +
               b_aug[i2k, i1k] * b_aug[i1k, i2k])

    c(out1, lav_mat_vech(var_fs2k - lambda2_eff * tmpbar))
  }

  cveta <- numDeriv::jacobian(func = lbar, x = x_step1)

  gamma_addition <- n * (cveta %*% step1$Sigma.11 %*% t(cveta))
  gamma_addition
}
