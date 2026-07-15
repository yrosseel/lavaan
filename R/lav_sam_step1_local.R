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
  d_std_all <- vector("list", nblocks)

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
    d_std_all[[b]] <- d_std
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
    if (lv_interaction_flag) {
      lav_msg_stop(gettext(
        "sam() does not support latent interactions in combination with
        conditional.x = TRUE (yet)."))
    }
    res_slopes <- vector("list", length = nblocks)
    cov_x  <- fit@h1$implied$cov.x
    mean_x <- fit@h1$implied$mean.x
    for (b in seq_len(nblocks)) {
      x_names <- fit@Data@ov.names.x[[b]]
      # y-level slopes: read from step1$SLOPES (set by lav_sam_get_cov_ybar,
      # and perturbed by the categorical numeric jacobian); fall back to the
      # h1 implied slopes
      slopes_b <- step1$SLOPES[[b]]
      if (is.null(slopes_b)) {
        slopes_b <- fit@h1$implied$res.slopes[[b]]
      }
      # categorical/correlation: rescale to the model's latent-response
      # metric, exactly like COV/YBAR above
      if (fit@Data@nlevels == 1L &&
          (fit@Model@categorical || fit@Model@correlation)) {
        slopes_b <- (1 / drop(mm_delta[[b]])) * slopes_b
      }
      # slopes of the latent variables on the exogenous covariates; if
      # std.lv = TRUE, put them (and EETA, the latent intercepts) on the same
      # standardized metric as the cov2cor() rescaled VETA above, so that the
      # full conditional moment vector (intercepts, slopes, residual
      # covariances) is metric-coherent
      res_slopes[[b]] <- (m[[b]] %*% slopes_b) / d_std_all[[b]]
      if (fit@Options$std.lv && !is.null(eeta)) {
        eeta[[b]] <- eeta[[b]] / d_std_all[[b]]
      }
      # (dim)names, so that the step-2 summary-statistics call can align
      # these statistics with its own internal variable order
      rownames(res_slopes[[b]]) <- colnames(veta[[b]])
      colnames(res_slopes[[b]]) <- x_names
      dimnames(cov_x[[b]]) <- list(x_names, x_names)
      mean_x[[b]] <- drop(mean_x[[b]])
      names(mean_x[[b]]) <- x_names
    }
    attr(veta, "res.slopes") <- res_slopes
    attr(veta, "cov.x") <- cov_x
    attr(veta, "mean.x") <- mean_x
    # needed by the (analytic) conditional.x JACb in the Gamma.eta machinery
    step1$D.STD <- d_std_all
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


# Memory-lean Gamma.eta for group g: JAC %*% Gamma %*% t(JAC) without ever
# forming the p* x p* observed Gamma. Returns NULL when the (continuous, simple)
# setting does not apply, so the caller falls back to the explicit Gamma.
#
# The ADF Gamma of the WLS.obs statistics is the casewise crossprod
#   Gamma = crossprod(z) / n.eff,   z_i = [ y_i (meanstructure),
#                                           vech((y_i - ybar)(y_i - ybar)') ]
# so JAC %*% Gamma %*% t(JAC) = crossprod(z %*% t(JAC)) / n.eff (the U below).
#
# With the unbiased (Browne) Gamma -- the sam() default -- the cov block also
# needs the normal-theory Gamma and a rank-1 term:
#   Gamma_u_cov = c1*Gamma_cov - c2*(Gamma_NT - 2/(n-1)*tcrossprod(vech(S)))
# with c1 = n(n-1)/((n-2)(n-3)), c2 = n/((n-2)(n-3)). The NT sandwich is built
# without the p* x p* Gamma_NT via [JAC Gamma_NT JAC']_rc = 2*tr(G_r S G_c S),
# where G_r is the symmetric matrix with vech G_r = (cov rows of) JAC[r, ]
# (off-diagonals halved). The mean-mean block is unchanged and the mean-cov
# block is scaled by n/(n-2).
lav_sam_gamma_eta_g <- function(fit = NULL, jac_g = NULL, g = 1L) {
  lavoptions <- fit@Options
  lavmodel <- fit@Model
  unbiased <- isTRUE(lavoptions$gamma.unbiased)
  if (lavmodel@categorical || lavmodel@correlation || lavmodel@conditional.x ||
      (isTRUE(lavoptions$fixed.x) && length(fit@Data@ov.names.x[[g]]) > 0L) ||
      fit@Data@nlevels > 1L || length(fit@Data@cluster) > 0L ||
      (unbiased && isTRUE(lavoptions$gamma.n.minus.one))) {
    return(NULL)
  }
  m_y <- unname(as.matrix(fit@Data@X[[g]]))
  if (anyNA(m_y)) {
    return(NULL) # missing data: leave it to the (pairwise) NACOV path
  }
  n <- nrow(m_y)
  p <- ncol(m_y)
  meanstr <- lavmodel@meanstructure
  m_yc <- t(t(m_y) - colMeans(m_y))
  idx1 <- lav_mat_vech_col_idx(p)
  idx2 <- lav_mat_vech_row_idx(p)
  zc_cov <- m_yc[, idx1, drop = FALSE] * m_yc[, idx2, drop = FALSE]
  zc_cov <- t(t(zc_cov) - colMeans(zc_cov))

  if (meanstr) {
    jac_m <- jac_g[, seq_len(p), drop = FALSE]
    jac_c <- jac_g[, -seq_len(p), drop = FALSE]
    u_m <- m_yc %*% t(jac_m)
  } else {
    jac_c <- jac_g
  }
  u_c <- zc_cov %*% t(jac_c)

  # biased ADF Gamma: simple casewise crossprod
  if (!unbiased) {
    n_eff <- if (isTRUE(lavoptions$gamma.n.minus.one)) n - 1L else n
    u <- if (meanstr) u_m + u_c else u_c
    return(crossprod(u) / n_eff)
  }

  # unbiased (Browne): cov block needs the NT Gamma + rank-1 correction
  s_cov <- crossprod(m_yc) / n
  cov_vech <- lav_mat_vech(s_cov)
  mstar <- nrow(jac_g)
  offu <- which(idx1 != idx2)
  # G_mat: rows = vec(G_r), the symmetric matrix with halved off-diagonals
  g_mat <- matrix(0, mstar, p * p)
  lin_11 <- (idx1 - 1L) * p + idx2
  lin_22 <- (idx2 - 1L) * p + idx1
  g_mat[, lin_11] <- jac_c
  g_mat[, lin_22] <- jac_c
  g_mat[, lin_11[offu]] <- jac_c[, offu, drop = FALSE] / 2
  g_mat[, lin_22[offu]] <- jac_c[, offu, drop = FALSE] / 2
  # M_mat: rows = vec(S G_r S)
  m_mat <- t(apply(g_mat, 1L, function(gv) {
    as.vector(s_cov %*% matrix(gv, p, p) %*% s_cov)
  }))
  gamma_nt <- 2 * tcrossprod(g_mat, m_mat) # JAC Gamma_NT JAC'  (mstar x mstar)
  d_vec <- drop(jac_c %*% cov_vech)
  c1 <- n * (n - 1) / ((n - 2) * (n - 3))
  c2 <- n / ((n - 2) * (n - 3))
  ge_cov <- c1 * (crossprod(u_c) / n) -
    c2 * (gamma_nt - 2 / (n - 1) * tcrossprod(d_vec))
  if (!meanstr) {
    return(ge_cov)
  }
  ge_mm <- crossprod(u_m) / n
  ge_mc <- (n / (n - 2)) * (crossprod(u_m, u_c) / n)
  ge_mm + ge_cov + ge_mc + t(ge_mc)
}

# Shared helpers for lav_sam_step1_local_jac() (single group) and
# lav_sam_step1_local_jac_mg() (multigroup) -- they build the same jacobian
# (JAC = JACc %*% JACa + JACb [+ psi channel]); the multigroup version simply
# stacks the per-group statistics and structural moments. Factoring the common
# pieces out keeps the two entry points thin and avoids drift between them.

# Indices (in step1$PT) of the free 'measurement' parameters whose influence
# d(theta.mm)/d(stats) forms the rows of JACa: factor loadings (=~), residual
# (co)variances of observed indicators (~~), intercepts (~1, if meanstructure),
# scaling factors (~*~) and thresholds (|, categorical), plus the higher-order
# loadings (=~ on an lv indicator) and latent (co)variances (~~ between lv
# indicators). One row per free parameter: by default de-duplicated on the
# (constrained) pt_1$free numbering; pass dedup_key = step1$PT.free to de-dup on
# the *unconstrained* numbering instead (needed for the p_only 'P' matrix under
# ceq.simple, where across-group equality-constrained parameters share a
# constrained free number but are distinct rows of P / step1.free.idx).
lav_sam_meas_keep_idx <- function(pt_1, ov_names, lavmodel, lv_ind,
                                  meanstructure, dedup_key = NULL) {
  if (is.null(dedup_key)) {
    dedup_key <- pt_1$free
  }
  dup <- duplicated(dedup_key)
  lambda_idx <- which(pt_1$op == "=~" & pt_1$free > 0L & !dup)
  theta_idx  <- which(pt_1$op == "~~" & pt_1$free > 0L & !dup &
                      pt_1$lhs %in% ov_names & pt_1$rhs %in% ov_names)
  nu_idx <- integer(0L)
  if (meanstructure) {
    nu_idx <- which(pt_1$op == "~1" & pt_1$free > 0L & !dup &
                    pt_1$lhs %in% ov_names)
  }
  delta_idx <- integer(0L)
  if (lavmodel@categorical || lavmodel@correlation) {
    delta_idx <- which(pt_1$op == "~*~" & pt_1$free > 0L & !dup)
  }
  th_idx <- integer(0L)
  if (lavmodel@categorical) {
    th_idx <- which(pt_1$op == "|" & pt_1$free > 0L & !dup)
  }
  beta_idx <- psi_idx <- integer(0L)
  if (length(lv_ind) > 0L) {
    beta_idx <- which(pt_1$op == "=~" & pt_1$free > 0L & !dup &
                      pt_1$rhs %in% lv_ind)
    psi_idx <- which(pt_1$op == "~~" & pt_1$free > 0L & !dup &
                     pt_1$rhs %in% lv_ind & pt_1$lhs %in% lv_ind)
  }
  sort(c(lambda_idx, theta_idx, nu_idx, delta_idx, th_idx, beta_idx, psi_idx))
}

# Stack the structural moments c(EETA_g, vech(VETA_g)) over the requested groups
# (a single g for a per-group JACb, or all groups for JACc/psi). This is the
# function whose jacobian (w.r.t. the sample statistics / measurement params)
# the numeric-difference helpers below compute.
#
# conditional.x: the structural moments follow the conditional WLS.obs layout
# c( vecr(cbind(EETA, RS_eta)), vech(VETA) ), where RS_eta = the res.slopes
# attribute of VETA (slopes of the latent variables on the exogenous
# covariates) and vecr interleaves per variable: (int_1, b_11..b_1q, int_2,
# ...). The variables are in VETA (measurement) order here; the permutation to
# the structural model's internal order happens once, on Gamma.eta (see
# lav_sam_condx_perm_idx()).
lav_sam_veta_vec <- function(step1_obj, groups, meanstructure,
                             conditional_x = FALSE) {
  unlist(lapply(groups, function(g) {
    if (conditional_x) {
      rs_g <- attr(step1_obj$VETA, "res.slopes")[[g]]
      c(lav_mat_vecr(cbind(step1_obj$EETA[[g]], rs_g)),
        lav_mat_vech(step1_obj$VETA[[g]]))
    } else if (meanstructure) {
      c(step1_obj$EETA[[g]], lav_mat_vech(step1_obj$VETA[[g]]))
    } else {
      lav_mat_vech(step1_obj$VETA[[g]])
    }
  }))
}

# Numeric jacobian of the stacked structural moments w.r.t. a set of
# free-parameter rows of step1$PT$est (rows_idx). Used for JACc (rows_idx =
# keep_idx) and the psi channel (rows_idx = the joint factor (co)variance rows).
lav_sam_jac_est_rows <- function(rows_idx, step1, fit, groups, meanstructure) {
  conditional_x <- fit@Model@conditional.x
  ff <- function(x) {
    step1_1 <- step1
    step1_1$PT$est[rows_idx] <- x
    step1_1 <- lav_sam_step1_local(step1 = step1_1, fit = fit,
      sam_method = step1$sam.method, local_options = step1$local.options)
    lav_sam_veta_vec(step1_1, groups, meanstructure, conditional_x)
  }
  verbose_flag <- lav_verbose()
  lav_verbose(FALSE)
  on.exit(lav_verbose(verbose_flag))
  numDeriv::jacobian(func = ff, x = step1$PT$est[rows_idx])
}

# Categorical JACb for one group g: numeric jacobian of vech(VETA_g) w.r.t. the
# off-diagonal (poly)choric correlations and the free continuous variances of
# step1$COV[[g]] (keeping theta.mm fixed), with the columns mapped (by label)
# into group g's joint WLS.obs vector. Thresholds/means stay zero. Under
# conditional.x the y-on-x slopes (step1$SLOPES[[g]], raw probit metric) are
# perturbed as well: the structural moments include RS_eta = M %*% slopes.
lav_sam_jacb_cat_g <- function(g, step1, fit, joint_labels_g, meanstructure) {
  conditional_x <- fit@Model@conditional.x
  cov_g <- step1$COV[[g]]
  low_idx <- lower.tri(cov_g) # strict lower triangle (column-major)
  n_low <- sum(low_idx)
  nvar_g <- ncol(cov_g)
  # continuous variables: their variances (the diagonal) are free statistics
  num_idx <- which(!seq_len(nvar_g) %in% fit@SampleStats@th.idx[[g]])
  x0 <- c(cov_g[low_idx], diag(cov_g)[num_idx])
  n_covvar <- length(x0)
  if (conditional_x) {
    slopes_g <- step1$SLOPES[[g]]
    x0 <- c(x0, as.vector(slopes_g)) # column-major, like lav_mat_vec()
  }
  ffb_cat <- function(xx) {
    step1_1 <- step1
    this_cov <- step1$COV[[g]]
    this_cov[low_idx] <- xx[seq_len(n_low)]
    this_cov[upper.tri(this_cov)] <- t(this_cov)[upper.tri(this_cov)]
    if (length(num_idx) > 0L) {
      dd <- diag(this_cov)
      dd[num_idx] <- xx[n_low + seq_along(num_idx)]
      diag(this_cov) <- dd
    }
    step1_1$COV[[g]] <- this_cov
    if (conditional_x) {
      step1_1$SLOPES[[g]] <- matrix(xx[-seq_len(n_covvar)],
                                    nrow = nrow(slopes_g),
                                    ncol = ncol(slopes_g))
    }
    step1_1 <- lav_sam_step1_local(step1 = step1_1, fit = fit,
      sam_method = step1$sam.method, local_options = step1$local.options)
    lav_sam_veta_vec(step1_1, g, meanstructure, conditional_x)
  }
  verbose_flag <- lav_verbose()
  lav_verbose(FALSE)
  jacb_corr <- numDeriv::jacobian(func = ffb_cat, x = x0)
  lav_verbose(verbose_flag)

  # map the perturbed-statistic columns into group g's joint WLS.obs vector (by
  # label, so any difference in variable ordering is handled)
  ov_g <- colnames(cov_g)
  if (is.null(ov_g)) {
    ov_g <- fit@Data@ov.names[[g]]
  }
  r_idx <- row(cov_g)[low_idx]
  c_idx <- col(cov_g)[low_idx]
  var_labels <- if (length(num_idx) > 0L) {
    paste0(ov_g[num_idx], "|var")
  } else {
    character(0L)
  }
  my_labels <- c(vapply(seq_along(r_idx), function(k) {
    paste(sort(c(ov_g[r_idx[k]], ov_g[c_idx[k]])), collapse = "~~")
  }, character(1L)), var_labels)
  if (conditional_x) {
    x_names_g <- fit@Data@ov.names.x[[g]]
    my_labels <- c(my_labels,
      paste0(rep(ov_g, times = length(x_names_g)), "~",
             rep(x_names_g, each = nvar_g)))
  }
  col_idx <- match(my_labels, joint_labels_g)
  if (anyNA(col_idx)) {
    lav_msg_stop(gettext(
      "internal error: unable to map the polychoric correlations into the
       statistics vector of the joint model (categorical JACb)"))
  }
  jb <- matrix(0, nrow = nrow(jacb_corr),
               ncol = length(fit@SampleStats@WLS.obs[[g]]))
  jb[, col_idx] <- jacb_corr
  jb
}

# Continuous JACb for one group g: vech(VETA_g) = M_g S_g M_g', so
# d vech(VETA_g) / d vech(S_g) is the (row-selected) duplication-post product of
# M_g %x% M_g, with a leading M_g block for the means (if meanstructure).
lav_sam_jacb_cont_g <- function(mb, meanstructure) {
  mbx_mb <- mb %x% mb
  row_idx <- lav_mat_vech_idx(nrow(mb))
  jb <- lav_mat_dup_post(mbx_mb)[row_idx, , drop = FALSE]
  if (meanstructure) {
    jb <- lav_mat_bdiag(mb, jb)
  }
  jb
}

# Continuous JACb for one group g under conditional.x. Rows follow the
# conditional structural-moment layout c(vecr(cbind(EETA, RS_eta)),
# vech(VETA)) of lav_sam_veta_vec(); columns follow the joint conditional
# WLS.obs layout c(vecr(cbind(res.int, res.slopes)), vech(res.cov)).
# With C_eta = cbind(EETA, RS_eta) = M_g %*% C_y - [M_g nu, 0] and
# vecr(A) = vec(A'):
#   d vecr(C_eta) / d vecr(C_y)     = M_g %x% I_(q+1)
#   d vech(VETA)  / d vech(res.cov) = as in the unconditional case
# (theta.mm -- and hence M_g -- is held fixed here; its influence goes
# through the JACc channel). If std.lv = TRUE, the intercept/slope rows are
# on the same standardized metric as VETA (divide the M_g rows by d_std,
# matching the res.slopes attribute); the vech(VETA) block keeps the raw
# M_g %x% M_g map, as in the unconditional case.
lav_sam_jacb_cont_condx_g <- function(mb, nexo, d_std = NULL) {
  mb_m <- mb
  if (!is.null(d_std)) {
    mb_m <- mb / d_std
  }
  jb_mean <- mb_m %x% diag(nexo + 1L)
  row_idx <- lav_mat_vech_idx(nrow(mb))
  jb_cov <- lav_mat_dup_post(mb %x% mb)[row_idx, , drop = FALSE]
  lav_mat_bdiag(jb_mean, jb_cov)
}

# conditional.x, continuous measurement blocks: the blocks themselves are fit
# unconditionally (their unconditional sample moments are exact functions of
# the joint conditional statistics, given the -- fixed -- moments of the
# exogenous covariates):
#   ybar_B = res.int_B + res.slopes_B %*% mean.x
#   S_BB   = res.cov_BB + res.slopes_B %*% cov.x %*% res.slopes_B'
# This helper returns the jacobian T of the block statistics c(ybar_B,
# vech(S_BB)) w.r.t. the joint conditional WLS.obs vector, so that the
# block influence (columns = block statistics) can be mapped into the joint
# conditional statistics space as mm_jac %*% T.
lav_sam_condx_t_g <- function(fit, g, ov_block, meanstructure) {
  ov_nox <- fit@pta$vnames$ov.nox[[g]]
  ov_x <- fit@Data@ov.names.x[[g]]
  ny <- length(ov_nox)
  q <- length(ov_x)
  bs <- fit@h1$implied$res.slopes[[g]] # ny x q
  bsx <- bs %*% fit@h1$implied$cov.x[[g]] # ny x q
  xbar <- drop(fit@h1$implied$mean.x[[g]])

  y_idx <- match(ov_block, ov_nox)
  if (anyNA(y_idx)) {
    lav_msg_stop(gettext(
      "internal error: unable to map the measurement-block variables into
      the joint conditional statistics vector (conditional.x JACa)"))
  }
  nb <- length(y_idx)
  ncol_t <- ny * (1L + q) + ny * (ny + 1L) / 2L
  cov_off <- ny * (1L + q)
  # column index of the intercept / slopes of joint y-variable j
  slope_cols <- function(j) (j - 1L) * (q + 1L) + 1L + seq_len(q)
  # vech position of the joint pair (j, l)
  pos <- matrix(0L, ny, ny)
  pos[lower.tri(pos, diag = TRUE)] <- seq_len(ny * (ny + 1L) / 2L)
  pos[upper.tri(pos)] <- t(pos)[upper.tri(pos)]

  t_mean <- NULL
  if (meanstructure) {
    t_mean <- matrix(0, nb, ncol_t)
    for (i in seq_len(nb)) {
      j <- y_idx[i]
      t_mean[i, (j - 1L) * (q + 1L) + 1L] <- 1
      if (q > 0L) {
        t_mean[i, slope_cols(j)] <- xbar
      }
    }
  }
  r_i <- lav_mat_vech_row_idx(nb)
  c_i <- lav_mat_vech_col_idx(nb)
  t_cov <- matrix(0, nb * (nb + 1L) / 2L, ncol_t)
  for (k in seq_len(nrow(t_cov))) {
    j <- y_idx[r_i[k]]
    l <- y_idx[c_i[k]]
    t_cov[k, cov_off + pos[j, l]] <- 1
    if (q > 0L) {
      t_cov[k, slope_cols(j)] <- t_cov[k, slope_cols(j)] + bsx[l, ]
      t_cov[k, slope_cols(l)] <- t_cov[k, slope_cols(l)] + bsx[j, ]
    }
  }
  rbind(t_mean, t_cov)
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
  # categorical local SEs: supported for the single-group setting. The
  # categorical d(VETA)/d(stats) jacobian needs an extra channel that the
  # continuous case does not: in the delta (correlation) parameterization the
  # residual variances theta are not free but are pinned by the unit-diagonal
  # constraint theta_jj = 1 - (Lambda Psi Lambda')_jj. Hence VETA also depends
  # on the within-block factor (co)variances Psi (through theta, and so through
  # the mapping matrix M). Psi is not a measurement parameter in 'keep_idx', so
  # this 'psi channel' is otherwise dropped; without it the latent (co)variance
  # SEs come out ~0.6x of their true value. It is added to JAC below (the
  # continuous case is unaffected, because there theta is free and
  # d(VETA)/d(Psi) = 0). Block-wise mixed models (some measurement blocks fully
  # categorical, others fully continuous) are supported: the JACa loop picks the
  # influence path per block, and JACb also perturbs the free continuous
  # variances. Higher-order categorical models are supported too (the structural
  # part reduces to the second-order factor (co)variances, as in the continuous
  # case). Under conditional.x the measurement blocks are themselves fitted
  # with conditional.x = TRUE (same latent-response scale as the joint model;
  # see lav_sam_step1()), the joint statistics gain the y-on-x slope entries
  # (label "y~x"), and JACb also perturbs step1$SLOPES. Block-wise mixed +
  # conditional.x is stopped in lav_sam_step1().

  # multiple groups: dispatch to the (single-level, no latent interaction)
  # multigroup implementation, which builds one stacked jacobian over all
  # groups and returns a list of per-group jacobians. Across-group (equality)
  # constraints in the measurement model (Case B) are detected there (via the
  # cross-group blocks of the stacked jacobian) and handled with the full
  # cross-group Gamma in step 2 (SEs and the robust test statistic).
  if (ngroups > 1L) {
    if (lavdata@nlevels > 1L) {
      lav_msg_stop(gettext(
        "local SEs: not available for multigroup multilevel models (yet)!"))
    }
    if (lavmodel@conditional.x) {
      lav_msg_stop(gettext(
        "local SEs: not available with multiple groups in combination with
         conditional.x = TRUE (yet)!"))
    }
    if (length(unlist(lavpta$vnames$lv.interaction)) > 0L) {
      lav_msg_stop(gettext(
        "local SEs: not available for latent interactions with multiple
         groups (yet)!"))
    }
    return(lav_sam_step1_local_jac_mg(step1 = step1, fit = fit,
                                      p_only = p_only,
                                      return_jac = return_jac))
  }

  g <- 1L
  n_mmblocks <- length(step1$MM.FIT)

  # categorical: we map the statistics of each measurement block into the
  # WLS statistics vector of the joint model by label
  if (lavmodel@categorical) {
    joint_labels <- lav_sam_wls_obs_labels(
      ov_names = if (lavmodel@conditional.x) {
        lavpta$vnames$ov.nox[[g]]
      } else {
        lavdata@ov.names[[g]]
      },
      th_idx = fit@SampleStats@th.idx[[g]],
      ov_names_x = if (lavmodel@conditional.x) {
        lavdata@ov.names.x[[g]]
      } else {
        NULL
      }
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
  # 'psi channel' accumulators (categorical only): the influence d(Psi)/d(stats)
  # for the within-block factor (co)variances, and their joint partable rows
  psi_infl_list <- list()
  psi_jrow_list <- integer(0L)
  for (mm in seq_len(n_mmblocks)) {
    fit_mm_block <- step1$MM.FIT[[mm]]

    # influence d(theta.mm)/d(block stats). When the joint model is categorical
    # we use the (D)WLS sandwich ingredients for EVERY block, including a fully
    # continuous one (block-wise mixed): lavTech(., "h1.information.expected")
    # is not available for a continuous measurement block that is embedded in a
    # categorical SAM, whereas lav_model_nvcov_robust_sem() returns a valid
    # Delta/E.inv/WLS.V for it.
    if (lavmodel@categorical) {
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
      # JACa = d(theta.mm)/d(stats) = invI . Delta' . A1 (A1 = saturated/h1
      # weight). Which information to invert depends on the SE method, and the
      # two choices are NOT interchangeable in principle (they coincide only for
      # a just-identified measurement block, where s = sigma):
      #   - se = "twostep.robust" (p_only): the Yuan & Chan (2002, eq. 12a)
      #     ASYMPTOTIC sandwich, whose bread is the population A-matrix -> the
      #     EXPECTED information (Delta'A1Delta)^-1, when information = "expected"
      #     (the lavaan default; honor the option otherwise).
      #   - se = "local"/"ij": the infinitesimal jackknife = an EMPIRICAL,
      #     finite-sample sandwich. Its influence d(theta.mm_hat)/d(stats) is the
      #     exact realized jacobian (implicit function theorem at the solution),
      #     i.e. the OBSERVED information; this also keeps JACa coherent with the
      #     observed/ADF Gamma it is later sandwiched with in Gamma.eta =
      #     JAC . Gamma . JAC'. So local MUST use observed (do not "unify" this
      #     with the twostep.robust choice).
      if (p_only && fit@Options$information[1] == "expected") {
        mm_inv_observed <-
          lavTech(fit_mm_block, "inverted.information.expected")
      } else {
        mm_inv_observed <-
          lavTech(fit_mm_block, "inverted.information.observed")
      }
      mm_jac <- t(mm_h1_expected[[g]] %*% mm_delta[[g]] %*% mm_inv_observed)
    }
    # full (unfiltered) influence: needed for the psi channel below
    mm_jac_full <- mm_jac

    # keep only rows that are also in FIT@ParTable. The rows of mm_jac are in
    # the *unconstrained* parameter order (Delta and the augmented inverted
    # information keep one column/row per unconstrained parameter), so under
    # within-block (ceq.simple) equality constraints the block's constrained
    # 'free' numbering must be remapped to the unconstrained numbering first
    # (the constrained members then each select their own -- identical --
    # influence row); indexing by the constrained numbering would shift every
    # row after the first constrained pair by one. (Same remap as the
    # Sigma.11 assembly in lav_sam_step1().)
    mm_ptm_free <- fit_mm_block@ParTable$free
    if (fit_mm_block@Model@ceq.simple.only) {
      mm_ptm_free[mm_ptm_free > 0L] <- seq_len(fit_mm_block@Model@nx.unco)
    }
    mm_keep_idx <- mm_ptm_free[step1$block.ptm.idx[[mm]]]
    mm_jac <- mm_jac[mm_keep_idx, , drop = FALSE]

    # map the columns of mm_jac (the sample statistics of this
    # measurement block) into the statistics vector of the joint model
    mm_row_idx <- step1$block.mm.idx[[mm]][step1$block.ptm.idx[[mm]]]
    if (lavmodel@categorical) {
      # label-based mapping (handles any variable ordering and a block-wise
      # mix of categorical and continuous measurement blocks)
      mm_col_idx <- lav_sam_mg_cat_col_idx(fit_mm_block, g, joint_labels,
                                           ncol(mm_jac_full))
      # sign flips (continuous conditional.x block: block intercept vs the
      # joint NEGATIVE intercept |t1 entry)
      mm_sign <- attr(mm_col_idx, "sign")
      if (!is.null(mm_sign)) {
        mm_jac <- t(t(mm_jac) * mm_sign)
      }
      jaca[mm_row_idx, mm_col_idx] <- mm_jac
    } else if (lavmodel@conditional.x) {
      # the (unconditional) block statistics are exact functions of the joint
      # conditional statistics: chain through T (see lav_sam_condx_t_g())
      t_block <- lav_sam_condx_t_g(fit, g,
        ov_block = fit_mm_block@Data@ov.names[[g]],
        meanstructure = fit_mm_block@Model@meanstructure)
      jaca[mm_row_idx, ] <- mm_jac %*% t_block
    } else {
      mm_ov_idx <- match(step1$MM.FIT[[mm]]@Data@ov.names[[g]],
                         lavdata@ov.names[[g]])
      mm_nvar <- length(lavdata@ov.names[[g]])
      mm_col_idx <- lav_mat_vech_which_idx(mm_nvar, idx = mm_ov_idx,
        add_idx_at_start = lavmodel@meanstructure)
      jaca[mm_row_idx, mm_col_idx] <- mm_jac
    }

    # psi channel: capture the influence d(Psi)/d(stats) for the within-block
    # factor (co)variances (categorical blocks only), mapped to the joint stats
    if (lavmodel@categorical && fit_mm_block@Model@categorical) {
      blv <- unlist(fit_mm_block@pta$vnames$lv)
      bpt <- fit_mm_block@ParTable
      bpsi <- which(bpt$op == "~~" & bpt$lhs %in% blv & bpt$rhs %in% blv &
                    bpt$free > 0L & !duplicated(bpt$free))
      for (k in bpsi) {
        ja_row <- numeric(ncol(jaca))
        ja_row[mm_col_idx] <- mm_jac_full[mm_ptm_free[k], ]
        # joint partable row for this factor (co)variance (canonical pair)
        jrow <- which(step1$PT$op == "~~" &
          ((step1$PT$lhs == bpt$lhs[k] & step1$PT$rhs == bpt$rhs[k]) |
           (step1$PT$lhs == bpt$rhs[k] & step1$PT$rhs == bpt$lhs[k])))[1]
        psi_infl_list[[length(psi_infl_list) + 1L]] <- ja_row
        psi_jrow_list <- c(psi_jrow_list, jrow)
      }
    }
  }

  # keep only the free 'measurement' parameters
  pt_1 <- step1$PT
  # only ov.names that are actually used in the measurement models
  ov_names <- unique(unlist(lapply(step1$MM.FIT, lav_object_vnames, "ov")))
  lv_ind <- unlist(lavpta$vnames$lv.ind)
  keep_idx <- lav_sam_meas_keep_idx(pt_1, ov_names, lavmodel, lv_ind,
                                    lavmodel@meanstructure)
  if (p_only) {
    # the rows of P are in ascending free-parameter order; return the
    # free indices, so the caller can align them with step1.free.idx
    # (which is ordered per measurement block). Under (within-block)
    # ceq.simple equality constraints, step1.free.idx is in the
    # *unconstrained* numbering (one entry per constrained member), so
    # de-duplicate on step1$PT.free instead of pt_1$free: the constrained
    # members then each keep their own (identical) influence row, exactly
    # as in the multigroup (_mg) path.
    p_keep_idx <- lav_sam_meas_keep_idx(pt_1, ov_names, lavmodel, lv_ind,
                                        lavmodel@meanstructure,
                                        dedup_key = step1$PT.free)
    jaca <- jaca[p_keep_idx, , drop = FALSE]
    attr(jaca, "free.idx") <- step1$PT.free[p_keep_idx]
    return(jaca)
  }
  jaca <- jaca[keep_idx, , drop = FALSE]

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
  } else if (lavmodel@categorical) {
    # categorical (and block-wise mixed): JACb depends on the (poly)choric
    # correlations and the free continuous variances of step1$COV (not the
    # thresholds/means), and through the DELTA rescaling / Croon correction is
    # not a simple M %x% M map; obtain it numerically (see lav_sam_jacb_cat_g()).
    jacb <- lav_sam_jacb_cat_g(g, step1, fit, joint_labels,
                               lavmodel@meanstructure)
  } else if (lavmodel@conditional.x) {
    jacb <- lav_sam_jacb_cont_condx_g(step1$M[[g]],
      nexo = length(fit@Data@ov.names.x[[g]]),
      d_std = step1$D.STD[[g]])
  } else { # no latent interactions
    jacb <- lav_sam_jacb_cont_g(step1$M[[g]], lavmodel@meanstructure)
  }

  # JACc: jacobian of vech(VETA) = f(theta.mm), keeping vech(S) fixed
  jacc <- lav_sam_jac_est_rows(keep_idx, step1, fit, g, lavmodel@meanstructure)

  # assemble JAC
  jac <- (jacc %*% jaca) + jacb

  # psi channel (categorical): add (d VETA / d Psi) %*% (d Psi / d stats), the
  # contribution of the within-block factor (co)variances that enter VETA
  # through the delta constraint theta_jj = 1 - (Lambda Psi Lambda')_jj. Without
  # it the categorical latent (co)variance SEs are ~0.6x too small. For
  # continuous indicators d VETA / d Psi = 0, so this leaves the continuous
  # (and mixed-within-block) results that do not depend on Psi unchanged.
  jacpsi <- NULL
  # drop any factor (co)variance already propagated as a measurement parameter
  # in keep_idx (e.g. the lv.ind (co)variances of a higher-order model), to
  # avoid double counting
  if (length(psi_jrow_list) > 0L) {
    dup_idx <- which(psi_jrow_list %in% keep_idx)
    if (length(dup_idx) > 0L) {
      psi_jrow_list <- psi_jrow_list[-dup_idx]
      psi_infl_list <- psi_infl_list[-dup_idx]
    }
  }
  if (lavmodel@categorical && length(psi_jrow_list) > 0L) {
    jaca_psi <- do.call(rbind, psi_infl_list)
    jacpsi <- lav_sam_jac_est_rows(psi_jrow_list, step1, fit, g,
                                   lavmodel@meanstructure)
    jac <- jac + (jacpsi %*% jaca_psi)
  }

  if (return_jac) {
    attr(jac, "JACa") <- jaca
    attr(jac, "JACb") <- jacb
    attr(jac, "JACc") <- jacc
    attr(jac, "JACpsi") <- jacpsi
  }

  # eventually, this will be a list per group/block
  list(jac)
}

# PML structural-space Gamma.eta (se = "local" for estimator = "PML").
#
# For PML the step-1 estimator is NOT a function of the saturated statistics
# s = (thresholds, polychoric correlations): the pairwise likelihood loads on
# the full bivariate tables. The influence route Gamma.eta = JAC Gamma JAC'
# (which writes theta1.hat as a function of s) therefore does not apply.
# Instead we build the CASEWISE joint influences of the two ingredients of
# the structural statistics c(EETA, vech(VETA)) = f(s, theta1.hat):
#
#   sqrt(N) (stats - target) ~ (1/sqrt(N)) sum_i u_i,
#   u_i = JACb zeta_i + JACc h_i
#
# with
#   zeta_i = casewise influence of the saturated statistics: the muthen1984
#            two-stage estimator's influence function, zeta_i = N B^-1 sc_i
#            (SC/A11/A21/A22 from muthen1984()); their crossprod/N
#            reproduces the categorical Gamma (NACOV) exactly;
#   h_i    = casewise influence of the step-1 (PML) block estimates: per
#            block, I_b^-1 g1_i with g1_i the casewise pairwise scores
#            (lav_sc()) and I_b the unit block information;
#   JACb   = d(stats)/d s'       (theta1 fixed; lav_sam_jacb_cat_g());
#   JACc   = d(stats)/d theta1'  (s fixed; lav_sam_jac_est_rows() over ALL
#            free block parameters, which automatically includes the
#            within-block factor (co)variances -- the 'psi channel' of the
#            (D)WLS categorical path: in the delta parameterization VETA
#            depends on Psi through theta_jj = 1 - (Lambda Psi Lambda')_jj).
#
# Gamma.eta = (1/N) sum_i (u_i - ubar)(u_i - ubar)': the joint casewise
# (infinitesimal-jackknife) covariance, including the cross-covariance
# between s and theta1.hat (same cases). Downstream, step 2 consumes it
# exactly like the (D)WLS Gamma.eta: robust (local) SEs + the corrected
# Satorra-Bentler structural test.
#
# Scope (guarded): single group, single level, complete data, all-ordinal,
# no exogenous covariates, no sampling weights, no equality constraints.
lav_sam_gamma_eta_pml <- function(step1 = NULL, fit = NULL) {
  lavdata <- fit@Data
  lavmodel <- fit@Model
  g <- 1L

  # guards
  if (lavdata@ngroups > 1L) {
    lav_msg_stop(gettext(
      "PML Gamma.eta: not available with multiple groups (yet)!"))
  }
  if (lavdata@nlevels > 1L) {
    lav_msg_stop(gettext(
      "PML Gamma.eta: not available for multilevel models (yet)!"))
  }
  if (lavmodel@conditional.x) {
    lav_msg_stop(gettext(
      "PML Gamma.eta: not available with exogenous covariates
       (conditional.x = TRUE) (yet)!"))
  }
  if (anyNA(lavdata@X[[g]])) {
    lav_msg_stop(gettext(
      "PML Gamma.eta: not available for missing data (yet)!"))
  }
  if (length(lavdata@sampling.weights) > 0L) {
    lav_msg_stop(gettext(
      "PML Gamma.eta: not available with sampling weights (yet)!"))
  }
  if (lavmodel@eq.constraints || lavmodel@ceq.simple.only) {
    lav_msg_stop(gettext(
      "PML Gamma.eta: not available with equality constraints (yet)!"))
  }
  ov_idx <- match(lavdata@ov.names[[g]], lavdata@ov$name)
  ov_types <- lavdata@ov$type[ov_idx]
  if (!all(ov_types == "ordered")) {
    lav_msg_stop(gettext(
      "PML Gamma.eta: all observed variables must be ordered (yet)!"))
  }

  y <- lavdata@X[[g]]
  n <- nrow(y)
  meanstructure <- lavmodel@meanstructure

  # --- 1. zeta_i: casewise influence of the saturated (th, rho) statistics
  cat_1 <- muthen1984(
    data_1 = y,
    ov_names = lavdata@ov.names[[g]],
    ov_types = ov_types,
    ov_levels = lavdata@ov$nlev[ov_idx],
    wls_w = TRUE
  )
  b_inv <- lav_m84_b_inv(
    a11 = cat_1$A11, a21 = cat_1$A21, a22 = cat_1$A22
  )$b_inv
  zeta <- n * (cat_1$SC %*% t(b_inv)) # N x nstats, WLS.obs layout

  # --- 2. JACb: d(stats)/d s' (theta1 fixed)
  joint_labels <- lav_sam_wls_obs_labels(
    ov_names = lavdata@ov.names[[g]],
    th_idx = fit@SampleStats@th.idx[[g]],
    ov_names_x = NULL
  )
  jacb <- lav_sam_jacb_cat_g(g, step1, fit, joint_labels, meanstructure)

  # --- 3. rows: ALL free block parameters, mapped to step1$PT rows
  pt_1 <- step1$PT
  n_mmblocks <- length(step1$MM.FIT)
  rows_list <- cols_list <- vector("list", n_mmblocks)
  for (mm in seq_len(n_mmblocks)) {
    fb <- step1$MM.FIT[[mm]]
    ptm <- fb@ParTable
    ptm_idx <- step1$block.ptm.idx[[mm]]
    ok <- ptm$free[ptm_idx] > 0L
    rows_list[[mm]] <- step1$block.mm.idx[[mm]][ptm_idx][ok]
    cols_list[[mm]] <- ptm$free[ptm_idx][ok]
  }
  rows_idx <- sort(unique(unlist(rows_list)))

  # JACc: d(stats)/d theta1' (s fixed), columns in rows_idx order
  jacc <- lav_sam_jac_est_rows(rows_idx, step1, fit, g, meanstructure)

  # --- 4. h_i: casewise influence of the block (PML) estimates,
  #        columns in rows_idx order
  h <- matrix(0, n, length(rows_idx))
  for (mm in seq_len(n_mmblocks)) {
    fb <- step1$MM.FIT[[mm]]
    # ignore_constraints: the blocks are fitted with bounds =
    # "wide.zerovar", which lavaan represents as (inactive) linear
    # INEQUALITY constraints; lav_sc()'s constraint projection would
    # project the scores onto them and corrupt the casewise influences.
    # Equality constraints are excluded from this route (guard above).
    scb <- lav_sc(fb, remove_empty_cases = FALSE,
                  ignore_constraints = TRUE)
    scb[is.na(scb)] <- 0
    if (nrow(scb) != n) {
      lav_msg_stop(gettext(
        "internal error: case alignment failure between a measurement
         block and the data (PML Gamma.eta)."))
    }
    ib <- lavTech(fb, "information") # unit information
    nb <- fb@SampleStats@ntotal
    hb <- (n / nb) * scb %*% solve(ib)
    h[, match(rows_list[[mm]], rows_idx)] <-
      hb[, cols_list[[mm]], drop = FALSE]
  }

  # --- 5. assemble the joint casewise influences + Gamma.eta
  u <- zeta %*% t(jacb) + h %*% t(jacc)
  u <- t(t(u) - colMeans(u))
  list(crossprod(u) / n)
}

# Two-level (single group) Gamma.eta: the NACOV of the structural statistics
# c(EETA_b, vech(VETA_b)) per level, plus their full (cross-level) covariance.
#
# The 'sample' statistics of the two-level local approach are the h1
# (saturated) estimates s = c(Mu.W, vech(Sigma.W), Mu.B, vech(Sigma.B)), whose
# sampling variability is the cluster-sandwich Gamma of the two-level (D)WLS
# machinery (see lav_samplestats_wls_2l.R):
#   Gamma.s = nobs * (I1^-1 J1 I1^-1) / nclusters
# with I1/J1 the expected/first-order per-cluster h1 informations at the h1
# estimates. The jacobian of the structural statistics w.r.t. s uses the same
# decomposition as the single-level lav_sam_step1_local_jac():
#   JAC_b = (JACc_b %*% JACa) + JACb_b
# where JACa (the influence of s on the step-1 measurement parameters) comes
# from the two-level measurement-block fits -- their Delta rows and h1
# information are laid out exactly like s, so the single-level influence
# formula invI . Delta' . A1 carries over, with the block statistics mapped
# into the joint s per level; JACb_b is the direct map of level b's segment
# of s (VETA_b = M_b S_b M_b'); and JACc_b is the numeric jacobian through
# the measurement estimates, reusing lav_sam_jac_est_rows() (which re-runs
# the two-level lav_sam_step1_local()).
#
# The structural model is refit as a two-'group' model (one group per level,
# nobs = c(N, nclusters)), so the per-level Gamma.eta list must follow
# lavaan's per-group NACOV convention Gamma_g = Cov(sqrt(n_g) stats_g):
#   Gamma.eta[[1]] =         JAC_1 Gamma.s JAC_1'    (n_1 = N; Gamma.s is
#                                                     N-scaled already)
#   Gamma.eta[[2]] = (G/N) * JAC_2 Gamma.s JAC_2'    (n_2 = G = nclusters)
# The two levels' statistics are correlated (they share s and the step 1
# chain), which the per-group list cannot express; we therefore also return
# the plain covariance of the STACKED statistics (Gamma.eta.full), and step 2
# computes the exact cross-level sandwich via the Case B path.
lav_sam_gamma_eta_2l <- function(step1 = NULL, fit = NULL) {
  lavdata <- fit@Data
  lavmodel <- fit@Model
  g <- 1L

  # supported setting: two levels, single group, continuous, complete data,
  # h1 statistics
  if (lavdata@ngroups > 1L) {
    lav_msg_stop(gettext(
      "two-level Gamma.eta: not available for multigroup multilevel
       models (yet)!"))
  }
  if (lavmodel@categorical) {
    lav_msg_stop(gettext(
      "two-level Gamma.eta: not available for categorical data (yet)!"))
  }
  if (length(unlist(fit@pta$vnames$lv.interaction)) > 0L) {
    lav_msg_stop(gettext(
      "two-level Gamma.eta: not available with latent interactions (yet)!"))
  }
  if (!is.null(step1$local.options$twolevel.method) &&
      step1$local.options$twolevel.method != "h1") {
    lav_msg_stop(gettext(
      "two-level Gamma.eta: only available for twolevel.method = \"h1\"."))
  }
  if (anyNA(lavdata@X[[g]])) {
    lav_msg_stop(gettext(
      "two-level Gamma.eta: not available for missing data (yet)!"))
  }

  lp <- lavdata@Lp[[g]]
  meanstructure <- lavmodel@meanstructure # always TRUE for two-level

  # level variable sets, in the order used by the h1 estimates and the
  # Delta rows (the Lp$ov.idx order; note ov.idx[[2]] may be unsorted)
  ov_l <- lapply(1:2, function(l) lavdata@ov.names[[g]][lp$ov.idx[[l]]])
  p_l <- lengths(ov_l)
  seg_len <- p_l + p_l * (p_l + 1L) / 2L
  seg_offset <- c(0L, seg_len[1]) # start of each level's (mean, vech) segment
  n_s <- sum(seg_len)

  # h1 (saturated) estimates: step1$COV/YBAR hold them per level (this is
  # the twolevel.method = "h1" path of lav_sam_get_cov_ybar())
  mu_w <- drop(step1$YBAR[[1]])
  sigma_w <- step1$COV[[1]]
  mu_b <- drop(step1$YBAR[[2]])
  sigma_b <- step1$COV[[2]]
  if (ncol(sigma_w) != p_l[1] || ncol(sigma_b) != p_l[2]) {
    lav_msg_stop(gettext(
      "internal error: two-level Gamma.eta: step1$COV dimensions do not
       match the level variable sets."))
  }

  # Gamma.s: cluster-sandwich NACOV of the h1 estimates (the same
  # computation as lav_samp_wls_2l(), N-scaled convention)
  nobs <- fit@SampleStats@nobs[[g]]
  nclusters <- lp$nclusters[[2]]
  i1 <- lav_mvn_cl_info_expected(
    lp = lp, mu_w = mu_w, sigma_w = sigma_w, mu_b = mu_b, sigma_b = sigma_b,
    x_idx = fit@SampleStats@x.idx[[g]]
  )
  j1 <- lav_mvn_cl_info_firstorder(
    y1 = lavdata@X[[g]], ylp = fit@SampleStats@YLp[[g]], lp = lp,
    mu_w = mu_w, sigma_w = sigma_w, mu_b = mu_b, sigma_b = sigma_b,
    x_idx = fit@SampleStats@x.idx[[g]], divide_by_two = TRUE
  )
  diag_i1 <- diag(i1)
  keep <- which(abs(diag_i1) > sqrt(.Machine$double.eps) * max(abs(diag_i1)))
  i1_inv <- matrix(0, nrow(i1), ncol(i1))
  inv_keep <- try(lav_mat_sym_inverse(i1[keep, keep, drop = FALSE]),
                  silent = TRUE)
  if (inherits(inv_keep, "try-error")) {
    lav_msg_stop(gettext(
      "two-level Gamma.eta: could not invert the h1 information matrix."))
  }
  i1_inv[keep, keep] <- inv_keep
  gamma_s <- nobs * (i1_inv %*% j1 %*% i1_inv) / nclusters
  gamma_s <- (gamma_s + t(gamma_s)) / 2

  # JACa: influence d(theta.mm)/d(s) of the two-level measurement blocks
  jaca <- matrix(0, nrow = length(fit@ParTable$lhs), ncol = n_s)
  for (mm in seq_len(length(step1$MM.FIT))) {
    fit_mm_block <- step1$MM.FIT[[mm]]
    if (fit_mm_block@Data@nlevels != 2L) {
      lav_msg_stop(gettext(
        "two-level Gamma.eta: all measurement blocks must be two-level
         models."))
    }
    # invI . Delta' . A1, exactly as in the single-level case (se = "local"
    # uses the observed information; see lav_sam_step1_local_jac())
    mm_h1_expected <- lavTech(fit_mm_block, "h1.information.expected")
    mm_delta <- lavTech(fit_mm_block, "Delta")
    mm_inv_observed <- lavTech(fit_mm_block, "inverted.information.observed")
    mm_jac <- t(mm_h1_expected[[g]] %*% mm_delta[[g]] %*% mm_inv_observed)

    # keep only rows that are also in FIT@ParTable
    mm_keep_idx <- fit_mm_block@ParTable$free[step1$block.ptm.idx[[mm]]]
    mm_jac <- mm_jac[mm_keep_idx, , drop = FALSE]

    # map the block's statistics (mu_w, vech(Sigma.W), mu_b, vech(Sigma.B)
    # over the block variables) into the joint s, per level
    mm_lp <- fit_mm_block@Data@Lp[[g]]
    mm_col_idx <- unlist(lapply(1:2, function(l) {
      mm_ov_l <- fit_mm_block@Data@ov.names[[g]][mm_lp$ov.idx[[l]]]
      idx_l <- match(mm_ov_l, ov_l[[l]])
      if (anyNA(idx_l)) {
        lav_msg_stop(gettext(
          "internal error: two-level Gamma.eta: unable to map the
           measurement-block variables into the joint statistics vector."))
      }
      seg_offset[l] + lav_mat_vech_which_idx(p_l[l], idx = idx_l,
                                             add_idx_at_start = TRUE)
    }))
    mm_row_idx <- step1$block.mm.idx[[mm]][step1$block.ptm.idx[[mm]]]
    jaca[mm_row_idx, mm_col_idx] <- mm_jac
  }

  # keep only the free 'measurement' parameters
  pt_1 <- step1$PT
  ov_names <- unique(unlist(lapply(step1$MM.FIT, lav_object_vnames, "ov")))
  lv_ind <- unlist(fit@pta$vnames$lv.ind)
  keep_idx <- lav_sam_meas_keep_idx(pt_1, ov_names, lavmodel, lv_ind,
                                    meanstructure)
  jaca <- jaca[keep_idx, , drop = FALSE]

  # per level: JAC_b = (JACc_b %*% JACa) + JACb_b
  jac <- vector("list", 2L)
  for (b in 1:2) {
    # direct part: (EETA_b, vech(VETA_b)) = f(mean_b, vech(S_b)), theta fixed
    jacb_core <- lav_sam_jacb_cont_g(step1$M[[b]], meanstructure)
    jacb <- matrix(0, nrow = nrow(jacb_core), ncol = n_s)
    jacb[, seg_offset[b] + seq_len(seg_len[b])] <- jacb_core
    # indirect part: through the measurement estimates
    jacc <- lav_sam_jac_est_rows(keep_idx, step1, fit, b, meanstructure)
    jac[[b]] <- (jacc %*% jaca) + jacb
  }

  # per-level NACOV (Gamma_g = Cov(sqrt(n_g) stats_g), n = c(N, nclusters))
  # + the full covariance of the stacked statistics (cross-level blocks)
  gamma_eta <- vector("list", 2L)
  gamma_eta[[1]] <- jac[[1]] %*% gamma_s %*% t(jac[[1]])
  gamma_eta[[2]] <- (nclusters / nobs) * (jac[[2]] %*% gamma_s %*% t(jac[[2]]))
  jac_full <- rbind(jac[[1]], jac[[2]])
  gamma_eta_full <- jac_full %*% (gamma_s / nobs) %*% t(jac_full)

  list(Gamma.eta = gamma_eta, Gamma.eta.full = gamma_eta_full, JAC = jac)
}

# multigroup version of lav_sam_step1_local_jac() (single level, no latent
# interactions, continuous data).
#
# Strategy: build ONE 'stacked' jacobian
#   JAC.full = JACc %*% JACa + JACb
# where
#   - the columns run over the stacked sample statistics vech(S_g) (per group),
#   - the rows run over the stacked c(EETA_g, vech(VETA_g)) (per group).
# The sample statistics are independent across groups, so the (full) Gamma of
# the stacked vech(S) is block-diagonal. When the measurement model has NO
# across-group (equality) constraints, the measurement parameters of the
# different groups are estimated independently, hence JAC.full is itself
# block-diagonal across groups, and we can return a simple list of per-group
# jacobians jac[[g]] (the caller then forms Gamma.eta[[g]] = jac[[g]] %*%
# Gamma[[g]] %*% t(jac[[g]]), and feeds it as a per-group NACOV to the
# multigroup structural model).
#
# If the measurement model DOES contain across-group constraints, JAC.full has
# nonzero off-diagonal (cross-group) blocks: vech(VETA_g) then depends on the
# sample statistics of the OTHER groups too. A per-group NACOV list can no
# longer capture this (the structural step would need the full cross-group
# Gamma.eta); we detect this situation and stop with an informative message.
# map the statistics of a (categorical or continuous) measurement block of
# group g into the labels of the joint WLS statistics vector of group g, and
# return the matching column indices (within group g's statistics).
lav_sam_mg_cat_col_idx <- function(fit_mm_block, g, joint_labels, ncol_full) {
  if (fit_mm_block@Model@categorical) {
    mm_labels <- lav_sam_wls_obs_labels(
      ov_names = if (fit_mm_block@Model@conditional.x) {
        fit_mm_block@pta$vnames$ov.nox[[g]]
      } else {
        fit_mm_block@Data@ov.names[[g]]
      },
      th_idx = fit_mm_block@SampleStats@th.idx[[g]],
      ov_names_x = if (fit_mm_block@Model@conditional.x) {
        fit_mm_block@Data@ov.names.x[[g]]
      } else {
        NULL
      })
  } else if (fit_mm_block@Model@conditional.x) {
    # continuous block, fitted with conditional.x = TRUE, within a
    # categorical + conditional.x joint model: the influence columns follow
    # the continuous conditional layout c(vecr(cbind(int, slopes)),
    # vech(res.cov, diag = TRUE)). The joint represents each numeric
    # variable by its NEGATIVE intercept (|t1), its slopes (~x), its
    # residual variance (|var) and residual covariances (~~); the |t1 sign
    # flip is returned in the "sign" attribute.
    block_ov <- fit_mm_block@Data@ov.names[[g]]
    block_x <- fit_mm_block@Data@ov.names.x[[g]]
    nvar_b <- length(block_ov)
    q_b <- length(block_x)
    mean_labels <- unlist(lapply(block_ov, function(y) {
      c(paste0(y, "|t1"), paste0(y, "~", block_x))
    }))
    mean_signs <- rep(c(-1, rep(1, q_b)), times = nvar_b)
    v_lin <- lav_mat_vech_idx(nvar_b)
    zero_b <- matrix(0L, nvar_b, nvar_b)
    v_r <- row(zero_b)[v_lin]
    v_c <- col(zero_b)[v_lin]
    cov_labels <- ifelse(v_r == v_c,
      paste0(block_ov[v_r], "|var"),
      vapply(seq_along(v_r), function(k)
        paste(sort(c(block_ov[v_r[k]], block_ov[v_c[k]])), collapse = "~~"),
        character(1L)))
    mm_labels <- c(mean_labels, cov_labels)
    mm_signs <- c(mean_signs, rep(1, length(cov_labels)))
    col_idx <- match(mm_labels, joint_labels)
    if (anyNA(col_idx)) {
      lav_msg_stop(gettext(
        "internal error: unable to map the statistics of a measurement block
         into the statistics vector of the joint model"))
    }
    attr(col_idx, "sign") <- mm_signs
    return(col_idx)
  } else {
    # continuous block within a categorical joint model: the influence columns
    # are ordered (means, vech(cov)); the joint represents each continuous
    # variable by its mean (|t1), variance (|var) and covariances (~~)
    block_ov <- fit_mm_block@Data@ov.names[[g]]
    nvar_b <- length(block_ov)
    v_lin <- lav_mat_vech_idx(nvar_b)
    zero_b <- matrix(0L, nvar_b, nvar_b)
    v_r <- row(zero_b)[v_lin]
    v_c <- col(zero_b)[v_lin]
    cov_labels <- ifelse(v_r == v_c,
      paste0(block_ov[v_r], "|var"),
      vapply(seq_along(v_r), function(k)
        paste(sort(c(block_ov[v_r[k]], block_ov[v_c[k]])), collapse = "~~"),
        character(1L)))
    if (ncol_full == nvar_b + length(v_lin)) {
      mm_labels <- c(paste0(block_ov, "|t1"), cov_labels) # with means
    } else {
      mm_labels <- cov_labels
    }
  }
  col_idx <- match(mm_labels, joint_labels)
  if (anyNA(col_idx)) {
    lav_msg_stop(gettext(
      "internal error: unable to map the statistics of a measurement block
       into the statistics vector of the joint model"))
  }
  col_idx
}

lav_sam_step1_local_jac_mg <- function(step1 = NULL, fit = NULL,
                                       p_only = FALSE,
                                       return_jac = FALSE) {
  lavdata  <- fit@Data
  lavmodel <- fit@Model
  lavpta   <- fit@pta
  ngroups  <- lavdata@ngroups
  n_mmblocks <- length(step1$MM.FIT)
  meanstructure <- lavmodel@meanstructure

  # per-group size of the (joint) sample-statistics vector vech(S_g)[+ means]
  stat_n <- vapply(seq_len(ngroups),
    function(g) length(fit@SampleStats@WLS.obs[[g]]), integer(1L))
  stat_offset <- c(0L, cumsum(stat_n)) # length ngroups + 1
  stat_tot <- sum(stat_n)

  # NOTE: veta_n / veta_offset (the per-group size of c(EETA_g, vech(VETA_g)))
  # are only needed for the JACb/JACc stacking below; they rely on step1$VETA,
  # which is only computed for the local/fsr/cfsr methods. They are therefore
  # deferred until after the p_only early return (twostep.robust, also used by
  # sam.method = "global", does not need them and does not have step1$VETA).

  # categorical: per-group labels of the joint WLS statistics vector, used to
  # map each measurement block's statistics into the right joint columns
  if (lavmodel@categorical) {
    joint_labels <- lapply(seq_len(ngroups), function(g)
      lav_sam_wls_obs_labels(ov_names = lavdata@ov.names[[g]],
                             th_idx = fit@SampleStats@th.idx[[g]]))
  }

  # JACa: jacobian of theta.mm = f(vech(S_all))
  # rows: all parameters in fit@ParTable (we select 'measurement' rows later)
  # cols: stacked vech(S_g) over the groups
  jaca <- matrix(0, nrow = length(fit@ParTable$lhs), ncol = stat_tot)
  # psi-channel accumulators (categorical only): the influence d(Psi)/d(stats)
  # of the within-block factor (co)variances, and their joint partable rows
  psi_infl_list <- list()
  psi_jrow_list <- integer(0L)
  for (mm in seq_len(n_mmblocks)) {
    fit_mm_block   <- step1$MM.FIT[[mm]]
    if (lavmodel@categorical) {
      # influence of the (D)WLS estimator: (Delta' W Delta)^-1 Delta' W, reusing
      # the robust-vcov ingredients. We use this for EVERY block (incl. a fully
      # continuous one in a block-wise mixed model), because
      # lavTech(., "h1.information.expected") is not available for a continuous
      # block embedded in a categorical SAM.
      tmp <- lav_model_nvcov_robust_sem(
        lavmodel = fit_mm_block@Model,
        lavsamplestats = fit_mm_block@SampleStats,
        lavcache = fit_mm_block@cache, lavdata = fit_mm_block@Data,
        lavimplied = fit_mm_block@implied, lavh1 = fit_mm_block@h1,
        lavoptions = fit_mm_block@Options, use_ginv = FALSE,
        attr_delta = TRUE, attr_e_inv = TRUE, attr_wls_v = TRUE)
      mm_e_inv  <- attr(tmp, "E.inv")
      mm_delta  <- attr(tmp, "Delta")
      mm_wls_v  <- attr(tmp, "WLS.V")
    } else {
      mm_h1_expected <- lavTech(fit_mm_block, "h1.information.expected")
      mm_delta       <- lavTech(fit_mm_block, "Delta")
      # same expected-vs-observed information choice as the single-group path
      # (see lav_sam_step1_local_jac()): EXPECTED for se = "twostep.robust"
      # (p_only; Yuan & Chan asymptotic sandwich, honoring the information
      # option), OBSERVED for se = "local"/"ij" (infinitesimal jackknife =
      # empirical sandwich, coherent with the observed/ADF Gamma).
      if (p_only && fit@Options$information[1] == "expected") {
        mm_inv_observed <-
          lavTech(fit_mm_block, "inverted.information.expected")
      } else {
        mm_inv_observed <-
          lavTech(fit_mm_block, "inverted.information.observed")
      }
    }

    # rows of mm_jac that we keep (the free params present in fit@ParTable),
    # together with the corresponding row numbers in fit@ParTable.
    # NOTE: Delta and (inverted) information have one row/col per *unconstrained*
    # free parameter (nx.unco); with simple equality constraints (ceq.simple,
    # e.g. across-group equal loadings) the ParTable$free numbering (nx.free)
    # is smaller and would mis-index the rows of mm_jac. Remap to the
    # unconstrained position (cf. lav_sam_step1()).
    mm_free <- fit_mm_block@ParTable$free
    if (fit_mm_block@Model@ceq.simple.only) {
      mm_free[mm_free > 0L] <- seq_len(fit_mm_block@Model@nx.unco)
    }
    mm_keep_idx <- mm_free[step1$block.ptm.idx[[mm]]]
    mm_row_idx  <- step1$block.mm.idx[[mm]][step1$block.ptm.idx[[mm]]]

    # psi channel: identify this block's free factor (co)variances (ALL groups)
    # ONCE, with their (unconstrained) block-row index and their joint partable
    # row. We accumulate d(Psi)/d(stats) over the groups below, so that with
    # across-group measurement constraints (Case B) a group's factor variance
    # also picks up its dependence on the OTHER groups' statistics.
    psi_fp <- integer(0L); psi_jr <- integer(0L)
    if (lavmodel@categorical && fit_mm_block@Model@categorical) {
      blv <- unlist(fit_mm_block@pta$vnames$lv)
      bpt <- fit_mm_block@ParTable
      bpsi <- which(bpt$op == "~~" & bpt$lhs %in% blv & bpt$rhs %in% blv &
                    bpt$free > 0L & !duplicated(bpt$free))
      for (k in bpsi) {
        psi_fp <- c(psi_fp, mm_free[k]) # unconstrained block-row index
        jrow <- which(step1$PT$op == "~~" & step1$PT$group == bpt$group[k] &
          ((step1$PT$lhs == bpt$lhs[k] & step1$PT$rhs == bpt$rhs[k]) |
           (step1$PT$lhs == bpt$rhs[k] & step1$PT$rhs == bpt$lhs[k])))[1]
        psi_jr <- c(psi_jr, jrow)
      }
    }
    psi_acc <- matrix(0, nrow = length(psi_fp), ncol = stat_tot)

    # lavaan stores a single 'unit' (per-observation averaged) information
    # matrix for the whole (multigroup) block; the group-g diagonal block of
    # its inverse is therefore (N/n_g) times the inverse one would obtain from
    # group g alone, while Delta[[g]] and h1.information.expected[[g]] are
    # unweighted. To recover the per-group influence d(theta.mm_g)/d(s_g) we
    # undo this normalization with the group proportion w_g = n_g / N.
    mm_ntotal <- fit_mm_block@SampleStats@ntotal

    for (g in seq_len(ngroups)) {
      w_g <- fit_mm_block@SampleStats@nobs[[g]] / mm_ntotal
      # influence of group g's sample statistics on ALL block parameters
      # (in the no-across-group-constraint case, only group g's parameters
      #  have nonzero rows; with across-group constraints, other groups'
      #  parameter rows are nonzero too -> cross-group coupling)
      if (lavmodel@categorical) {
        wv <- mm_wls_v[[g]]
        if (is.matrix(wv)) {
          mm_jac <- w_g * (mm_e_inv %*% t(wv %*% mm_delta[[g]]))
        } else {
          # DWLS/ULS: WLS.V holds the diagonal of the weight matrix
          mm_jac <- w_g * (mm_e_inv %*% t(wv * mm_delta[[g]]))
        }
      } else {
        mm_jac <- w_g *
          t(mm_h1_expected[[g]] %*% mm_delta[[g]] %*% mm_inv_observed)
      }
      # full (unfiltered) influence: needed for the psi channel below
      mm_jac_full <- mm_jac
      mm_jac <- mm_jac[mm_keep_idx, , drop = FALSE]

      # map the columns (block statistics of group g) into the joint
      # statistics vector of group g
      if (lavmodel@categorical) {
        mm_col_idx <- stat_offset[g] +
          lav_sam_mg_cat_col_idx(fit_mm_block, g, joint_labels[[g]],
                                 ncol(mm_jac_full))
      } else {
        mm_ov_idx <- match(fit_mm_block@Data@ov.names[[g]],
                           lavdata@ov.names[[g]])
        mm_nvar   <- length(lavdata@ov.names[[g]])
        mm_col_idx <- stat_offset[g] +
          lav_mat_vech_which_idx(mm_nvar, idx = mm_ov_idx,
            add_idx_at_start = meanstructure)
      }
      jaca[mm_row_idx, mm_col_idx] <- mm_jac

      # psi channel: accumulate d(Psi)/d(group-g stats) for this block's factor
      # (co)variances at group g's columns. Summed over the group loop, this
      # gives the full (cross-group) influence; for Case A only the own-group
      # columns are nonzero.
      if (length(psi_fp) > 0L) {
        psi_acc[, mm_col_idx] <- psi_acc[, mm_col_idx, drop = FALSE] +
          mm_jac_full[psi_fp, , drop = FALSE]
      }
    }

    # store the accumulated psi-channel influence rows for this block
    for (i in seq_along(psi_fp)) {
      psi_infl_list[[length(psi_infl_list) + 1L]] <- psi_acc[i, ]
      psi_jrow_list <- c(psi_jrow_list, psi_jr[i])
    }
  }

  # keep only the free 'measurement' parameters
  pt_1 <- step1$PT
  ov_names <- unique(unlist(lapply(step1$MM.FIT, lav_object_vnames, "ov")))
  lv_ind <- unlist(lavpta$vnames$lv.ind)
  keep_idx <- lav_sam_meas_keep_idx(pt_1, ov_names, lavmodel, lv_ind,
                                    meanstructure)

  if (p_only) {
    # twostep.robust: we only need the 'P' matrix = influence of the
    # measurement (step 1) parameters on the (stacked) sample statistics.
    # The rows are in ascending free-parameter order; the caller aligns them
    # with step1.free.idx and splits the columns per group using stat.offset.
    # (No need for the costly JACb/JACc numeric jacobians below; this branch
    # is therefore also usable for sam.method = "global", which does not
    # populate step1$VETA / step1$M / step1$COV.)
    #
    # keep_idx de-duplicates on pt_1$free (the *constrained* free number). With
    # across-group equality constraints in the measurement model (Case B)
    # several partable rows -- one per group -- share that number but are
    # distinct in step1.free.idx (the *unconstrained* numbering). For the P
    # matrix we need one row per unconstrained step 1 parameter, so de-duplicate
    # on step1$PT.free instead. (When there are no simple equality constraints
    # the two numberings coincide, so p_keep_idx == keep_idx.)
    p_keep_idx <- if (lavmodel@ceq.simple.only) {
      lav_sam_meas_keep_idx(pt_1, ov_names, lavmodel, lv_ind, meanstructure,
                            dedup_key = step1$PT.free)
    } else {
      keep_idx
    }
    pjac <- jaca[p_keep_idx, , drop = FALSE]
    attr(pjac, "free.idx") <- step1$PT.free[p_keep_idx]
    attr(pjac, "stat.offset") <- stat_offset
    return(pjac)
  }

  jaca <- jaca[keep_idx, , drop = FALSE]

  # per-group size of c(EETA_g, vech(VETA_g)) (needed for the JACb/JACc
  # stacking below; relies on step1$VETA, only computed for local/fsr/cfsr)
  veta_n <- vapply(seq_len(ngroups), function(g) {
    p <- ncol(step1$VETA[[g]])
    # NOTE: %/% binds tighter than *, so parenthesize the vech length
    # (p*(p+1)/2); without the parens p*(p+1L)%/%2L is wrong for even p
    (if (meanstructure) p else 0L) + (p * (p + 1L)) %/% 2L
  }, integer(1L))
  veta_offset <- c(0L, cumsum(veta_n))
  veta_tot <- sum(veta_n)

  # JACb: jacobian of vech(VETA) = f(vech(S)), keeping theta.mm fixed
  # (block-diagonal across groups: VETA_g = M_g S_g M_g')
  jacb_list <- vector("list", ngroups)
  if (lavmodel@categorical) {
    for (g in seq_len(ngroups)) {
      jacb_list[[g]] <- lav_sam_jacb_cat_g(g, step1, fit, joint_labels[[g]],
                                           meanstructure)
    }
  } else {
    for (g in seq_len(ngroups)) {
      jacb_list[[g]] <- lav_sam_jacb_cont_g(step1$M[[g]], meanstructure)
    }
  }
  jacb <- lav_mat_bdiag(jacb_list)

  # JACc: jacobian of vech(VETA) = f(theta.mm), keeping vech(S) fixed
  # (over the stacked c(EETA_g, vech(VETA_g)) of all groups)
  jacc <- lav_sam_jac_est_rows(keep_idx, step1, fit, seq_len(ngroups),
                               meanstructure)

  # assemble the stacked JAC
  jac_full <- (jacc %*% jaca) + jacb

  # psi channel (categorical): add (d VETA / d Psi) %*% (d Psi / d stats), the
  # contribution of the within-block factor (co)variances that enter VETA
  # through the delta constraint theta_jj = 1 - (Lambda Psi Lambda')_jj. Without
  # it the categorical latent (co)variance SEs are ~0.6x too small. (For
  # continuous indicators d VETA / d Psi = 0.)
  if (length(psi_jrow_list) > 0L) {
    # drop any factor (co)variance already propagated as a measurement
    # parameter in keep_idx (e.g. lv.ind (co)variances), to avoid double counting
    dup_idx <- which(psi_jrow_list %in% keep_idx)
    if (length(dup_idx) > 0L) {
      psi_jrow_list <- psi_jrow_list[-dup_idx]
      psi_infl_list <- psi_infl_list[-dup_idx]
    }
  }
  if (lavmodel@categorical && length(psi_jrow_list) > 0L) {
    jaca_psi <- do.call(rbind, psi_infl_list)
    jacpsi <- lav_sam_jac_est_rows(psi_jrow_list, step1, fit,
                                   seq_len(ngroups), meanstructure)
    jac_full <- jac_full + (jacpsi %*% jaca_psi)
  }

  # split into per-group jacobians, and detect cross-group (off-diagonal)
  # coupling. Negligible off-diagonal blocks => the measurement model has no
  # across-group constraints (Case A): the per-group jac[[g]] are sufficient,
  # and the caller forms a per-group NACOV list. Non-negligible off-diagonal
  # blocks => across-group measurement constraints (Case B): vech(VETA_g) then
  # depends on the OTHER groups' sample statistics, so the full (stacked)
  # jacobian must be kept and a cross-group sandwich computed in step 2.
  jac <- vector("list", ngroups)
  caseb_flag <- FALSE
  for (g in seq_len(ngroups)) {
    r_idx <- (veta_offset[g] + 1L):veta_offset[g + 1L]
    c_idx <- (stat_offset[g] + 1L):stat_offset[g + 1L]
    jac[[g]] <- jac_full[r_idx, c_idx, drop = FALSE]

    other_c <- if (stat_tot > length(c_idx)) {
      seq_len(stat_tot)[-c_idx]
    } else {
      integer(0L)
    }
    if (length(other_c) > 0L) {
      off_block <- max(abs(jac_full[r_idx, other_c, drop = FALSE]))
      this_block <- max(abs(jac[[g]]), 1)
      # numeric (finite-difference) noise is ~1e-7; genuine across-group
      # coupling is O(1), so a relative tolerance of 1e-4 separates them
      if (off_block > 1e-4 * this_block) {
        caseb_flag <- TRUE
      }
    }
  }

  # always expose the stacked jacobian + block layout, so step 2 can build the
  # cross-group sandwich when needed (Case B)
  attr(jac, "jac.full") <- jac_full
  attr(jac, "caseB") <- caseb_flag
  attr(jac, "stat.offset") <- stat_offset
  attr(jac, "veta.offset") <- veta_offset

  if (return_jac) {
    attr(jac, "JACa") <- jaca
    attr(jac, "JACb") <- jacb
    attr(jac, "JACc") <- jacc
  }

  jac
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
lav_sam_gamma_add <- function(step1 = NULL, fit = NULL, group = 1L,
                              use_analytic = TRUE) {

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
  lbar <- function(x, .return_all = FALSE) {
    pt_x <- pt_1
    pt_x$est[step1_idx] <- x
    this_lavmodel <- lav_model_set_parameters(lavmodel,
                       x = pt_x$est[pt_x$free > 0 & !duplicated(pt_x$free)])
    this_nu     <- drop(this_lavmodel@GLIST$nu)
    # no interaction columns!
    if (length(rm_idx) > 0L) {
      this_lambda <- this_lavmodel@GLIST$lambda[, -rm_idx, drop = FALSE]
    } else {
      this_lambda <- this_lavmodel@GLIST$lambda
    }
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
    m_pre_std <- this_m
    d_std <- NULL
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

    out_vec <- c(out1, lav_mat_vech(var_fs2k - lambda2_eff * tmpbar))
    if (!.return_all) {
      return(out_vec)
    }
    # intermediates needed by the analytic jacobian below
    list(vec = out_vec, nu = this_nu, lambda = this_lambda,
         theta = this_theta, m_pre_std = m_pre_std, m = this_m,
         b = b_aug, e = e_mat, fs = fs, fs2kc = fs2kc, d_std = d_std)
  }

  # jacobian of Lbar w.r.t. the step 1 parameters: cveta[a, k] = d Lbar[a]/d x_k
  # By default we use a closed-form (analytic) jacobian, which avoids the
  # repeated mapping-matrix solves and full factor-score recomputations of the
  # numerical (numDeriv) jacobian. The analytic path applies the chain rule
  # through the model-implied factor scores, with the key building block being
  # the differential of the mapping matrix M = (Lambda' W Lambda)^-1 Lambda' W
  # (W = I, S^-1, Theta^-1 for ULS/GLS/ML):
  #   dM = A dLambda' W Pp  -  M dLambda M  [ -  M dTheta Theta^-1 Pp (ML) ]
  # with A = (Lambda' W Lambda)^-1 and Pp = I - Lambda M. It is only taken for
  # the 'clean' continuous setting (standard mapping matrix, no dummy latent
  # variables, no empty Lambda columns, invertible Theta for ML); otherwise we
  # fall back to numDeriv. See lav_sam_gamma_add_jac().
  cveta <- NULL
  if (use_analytic) {
    cveta <- tryCatch(
      lav_sam_gamma_add_jac(
        base = lbar(x_step1, .return_all = TRUE),
        x_step1 = x_step1, lavmodel = lavmodel, pt_1 = pt_1,
        step1_idx = step1_idx, rm_idx = rm_idx,
        dummy_ov_idx = dummy_ov_idx, s_cov = step1$COV[[1]],
        method = step1$local.options$M.method, std_lv_flag = std_lv_flag,
        lambda1 = lambda1, lambda2_eff = lambda2_eff,
        i1k = i1k, i2k = i2k, y = y, n = n),
      error = function(e) NULL)
  }
  if (is.null(cveta)) {
    cveta <- numDeriv::jacobian(func = lbar, x = x_step1)
  }

  gamma_addition <- n * (cveta %*% step1$Sigma.11 %*% t(cveta))
  gamma_addition
}

# Closed-form jacobian of the augmented summary statistics Lbar w.r.t. the
# step 1 parameters (the analytic counterpart of numDeriv::jacobian(lbar)).
# Returns the d_L x n_par jacobian matrix, or NULL if the 'clean' continuous
# setting does not apply (the caller then falls back to numDeriv).
lav_sam_gamma_add_jac <- function(base = NULL, x_step1 = NULL, lavmodel = NULL,
                                  pt_1 = NULL, step1_idx = NULL, rm_idx = NULL,
                                  dummy_ov_idx = NULL, s_cov = NULL,
                                  method = "ML", std_lv_flag = FALSE,
                                  lambda1 = 1, lambda2_eff = 1,
                                  i1k = NULL, i2k = NULL, y = NULL, n = NULL) {
  Lam0 <- base$lambda; Th0 <- base$theta; nu0 <- base$nu
  M0 <- base$m_pre_std; Mf <- base$m; fs0 <- base$fs
  e0 <- base$e; b0 <- base$b; fs2kc0 <- base$fs2kc; d_std0 <- base$d_std
  S <- s_cov
  p <- nrow(Lam0); mlv <- ncol(Lam0)
  method <- toupper(method)

  # eligibility: no dummy lvs, no empty Lambda columns, standard mapping matrix
  if (length(dummy_ov_idx) > 0L) {
    return(NULL)
  }
  if (any(apply(Lam0, 2L, function(z) all(z == 0)))) {
    return(NULL)
  }

  # weight matrix W and A = (Lambda' W Lambda)^-1
  if (method == "ULS") {
    W <- diag(p)
  } else if (method == "GLS") {
    W <- solve(S)
  } else if (method == "ML") {
    if (any(abs(diag(Th0)) < 1e-4)) {
      return(NULL)
    }
    W <- solve(Th0)
  } else {
    return(NULL)
  }
  LtW <- crossprod(Lam0, W) # mlv x p
  A <- solve(LtW %*% Lam0)
  # is the mapping matrix the standard weighted pseudo-inverse? (if lbar used
  # a fallback -- tmat/ginv/ULS -- the closed form below does not apply)
  if (max(abs(A %*% LtW - M0)) > 1e-7) {
    return(NULL)
  }
  Pp <- diag(p) - Lam0 %*% M0 # residual projector (p x p)
  WPp <- W %*% Pp
  ThiPp <- if (method == "ML") W %*% Pp else NULL # Theta^-1 Pp

  # cheap GLIST extractor (set_parameters only; no mapping/scores). The
  # Lambda/Theta/nu entries are linear in the free parameters, so the forward
  # difference below is exact (up to round-off).
  get_lln <- function(x) {
    pt_x <- pt_1
    pt_x$est[step1_idx] <- x
    ml <- lav_model_set_parameters(lavmodel,
            x = pt_x$est[pt_x$free > 0 & !duplicated(pt_x$free)])
    if (length(rm_idx) > 0L) {
      lam <- ml@GLIST$lambda[, -rm_idx, drop = FALSE]
    } else {
      lam <- ml@GLIST$lambda
    }
    list(nu = drop(ml@GLIST$nu), lambda = lam, theta = ml@GLIST$theta)
  }

  Ymc <- t(t(y) - nu0) # n x p  = Y - 1 nu'
  np <- length(x_step1)
  cv <- matrix(0, length(base$vec), np)
  eps_g <- 1e-6
  for (k in seq_len(np)) {
    xk <- x_step1
    xk[k] <- xk[k] + eps_g
    gk <- get_lln(xk)
    dLam <- (gk$lambda - Lam0) / eps_g
    dTh  <- (gk$theta  - Th0)  / eps_g
    dnu  <- (gk$nu     - nu0)  / eps_g

    # raw mapping-matrix derivative (pre-std rescaling)
    dM0 <- A %*% crossprod(dLam, WPp) - M0 %*% (dLam %*% M0)
    if (!is.null(ThiPp)) {
      dM0 <- dM0 - M0 %*% (dTh %*% ThiPp)
    }

    # std.lv rescaling chain: Mf = M0 / d_std, d_std_j = sqrt(msm_j - l1 mtm_j)
    if (std_lv_flag) {
      d_msm <- 2 * rowSums((dM0 %*% S) * M0)
      d_mtm <- 2 * rowSums((dM0 %*% Th0) * M0) + rowSums((M0 %*% dTh) * M0)
      d_dstd <- (d_msm - lambda1 * d_mtm) / (2 * d_std0)
      dMf <- dM0 / d_std0 - M0 * (d_dstd / d_std0^2)
    } else {
      dMf <- dM0
    }

    # b_aug = bdiag(0, M Theta M'), augmented factor scores fs, C2, E
    dMtMt <- dMf %*% Th0 %*% t(Mf)
    d_b <- lav_mat_bdiag(0, dMtMt + Mf %*% dTh %*% t(Mf) + t(dMtMt))
    d_score <- Ymc %*% t(dMf) - matrix(drop(Mf %*% dnu), n, mlv, byrow = TRUE)
    d_fs <- cbind(0, d_score)
    d_e <- (crossprod(d_fs, fs0) + crossprod(fs0, d_fs)) / n - lambda1 * d_b

    # first-order part
    d_out1 <- d_e[cbind(i2k, i1k)]

    # second-order part
    d_fs2k <- d_fs[, i1k, drop = FALSE] * fs0[, i2k, drop = FALSE] +
              fs0[, i1k, drop = FALSE] * d_fs[, i2k, drop = FALSE]
    d_fs2kc <- t(t(d_fs2k) - colMeans(d_fs2k))
    d_var <- (crossprod(d_fs2kc, fs2kc0) + crossprod(fs2kc0, d_fs2kc)) / n
    d_tmp <- (d_e[i1k, i1k] * b0[i2k, i2k] + e0[i1k, i1k] * d_b[i2k, i2k] +
              d_b[i1k, i1k] * e0[i2k, i2k] + b0[i1k, i1k] * d_e[i2k, i2k] +
              d_e[i1k, i2k] * b0[i2k, i1k] + e0[i1k, i2k] * d_b[i2k, i1k] +
              d_e[i2k, i1k] * b0[i1k, i2k] + e0[i2k, i1k] * d_b[i1k, i2k] +
              d_b[i1k, i1k] * b0[i2k, i2k] + b0[i1k, i1k] * d_b[i2k, i2k] +
              d_b[i2k, i1k] * b0[i1k, i2k] + b0[i2k, i1k] * d_b[i1k, i2k])

    cv[, k] <- c(d_out1, lav_mat_vech(d_var - lambda2_eff * d_tmp))
  }

  cv
}
