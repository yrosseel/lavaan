# SAM step 2: estimate structural part

lav_sam_step2 <- function(step1 = NULL, fit = NULL,
                          sam_method = "local", struc_args = list()) {
  lavoptions <- fit@Options
  lavpta <- fit@pta
  nlevels <- lavpta$nlevels
  pt_1 <- step1$PT

  # Gamma available?
  gamma_flag <- (sam_method %in% c("local", "fsr", "cfsr") &&
                 !is.null(step1$Gamma.eta[[1]]))

  lv_names <- unique(unlist(fit@pta$vnames$lv.regular))

  # adjust options
  lavoptions_pa <- lavoptions
  # "yuan.chan" is a SAM-global test for the JOINT model, computed afterwards in
  # lav_sam_global_test(); the structural fit itself uses the ordinary test
  if (any(lavoptions_pa$test == "yuan.chan")) {
    lavoptions_pa$test <- "standard"
  }
  # the corrected two-step STRUCTURAL test (Satorra-Bentler, using Gamma.eta as
  # the NACOV of vech(VETA)) is the default test for sam.method = local/fsr/cfsr
  # whenever Gamma.eta is available -- INDEPENDENT of the requested SE. (For
  # se = "twostep"/"naive" the FIT.PA SEs below are not the final ones: twostep
  # SEs are recomputed in step 4, naive SEs are FIT.PA's plain vcov.)
  if (gamma_flag) {
    lavoptions_pa$test <- "satorra.bentler"
  }
  if (lavoptions_pa$se == "naive") {
    # naive SEs = FIT.PA's plain (standard) vcov
    lavoptions_pa$se <- "standard"
  } else if (lavoptions_pa$se %in% c("local", "local.nt")) {
    # local SEs ARE FIT.PA's robust.sem vcov (read back in lav_sam_step2_se())
    lavoptions_pa$se <- "robust.sem"
  } else if (gamma_flag) {
    # twostep / twostep.robust: the final SEs are recomputed in step 4
    # (lav_sam_step2_se). Use se = "standard" -- NOT "robust.sem" -- so that
    # FIT.PA's vcov stays the NAIVE (standard) one: the alpha.correction blend
    # in lav_sam_step2_se() reads it as 'vcov_naive'. The corrected structural
    # test is still computed (it is independent of the se).
    lavoptions_pa$se <- "standard"
  } else {
    # twostep or none, without Gamma.eta -> none
    lavoptions_pa$se <- "none"
  }
  if (!lavoptions_pa$conditional.x) {
    lavoptions_pa$fixed.x <- FALSE # until we fix this...
  }
  lavoptions_pa$categorical <- FALSE
  lavoptions_pa$.categorical <- FALSE
  lavoptions_pa$rotation <- "none"
  lavoptions_pa <- modifyList(lavoptions_pa, struc_args)

  if (gamma_flag) {
    lavoptions_pa$check.vcov <- FALSE # always non-pd
                                      # if interactions + fixed.x = FALSE
  }

  # override, no matter what
  lavoptions_pa$do.fit <- TRUE

  if (sam_method %in% c("local", "fsr", "cfsr")) {
    lavoptions_pa$missing <- "listwise"
    lavoptions_pa$sample.cov.rescale <- FALSE
    lavoptions_pa$loglik <- FALSE
    # first.order information (eg estimator = "MLF") needs raw data, which
    # the structural fit does not have (it is fitted from the estimated
    # latent moments VETA/EETA only); use the expected information instead
    if (any(lavoptions_pa$information == "first.order")) {
      lavoptions_pa$information <- rep.int("expected", 2L)
    }
  } else {
    lavoptions_pa$h1 <- FALSE
    lavoptions_pa$loglik <- FALSE
  }

  # construct PTS
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    # extract structural part
    pts <- lav_pt_subset_sm(pt_1,
      add_idx = TRUE,
      add_exo_cov = TRUE,
      fixed_x = lavoptions_pa$fixed.x,
      conditional_x = lavoptions_pa$conditional.x,
      free_fixed_var = TRUE,
      meanstructure = lavoptions_pa$meanstructure
    )

    # any 'extra' parameters: not (free) in PT, but free in PTS (user == 3)
    #  - fixed.x in PT, but fixed.x = FALSE is PTS
    #  - fixed-to-zero intercepts in PT, but free in PTS
    #  - add.exo.cov: absent/fixed-to-zero in PT, but add/free in PTS
    extra_id <- which(pts$user == 3L)

    # remove est/se/start columns
    pts$est <- NULL
    pts$se <- NULL
    pts$start <- NULL

    if (nlevels > 1L) {
      pts$level <- NULL
      pts$group <- NULL
      pts$group <- pts$block
      nobs_1 <- fit@Data@Lp[[1]]$nclusters
    } else {
      nobs_1 <- fit@Data@nobs
    }

    reg_idx <- attr(pts, "idx")
    attr(pts, "idx") <- NULL

    # edge case: conditional.x = TRUE, but the structural part contains no
    # exogenous covariates (eg they only affect indicators directly); the
    # conditional attributes of VETA (res.slopes/cov.x/mean.x) then have no
    # counterpart in the structural model -> drop them and fit the
    # structural part unconditionally
    if (lavoptions_pa$conditional.x &&
        length(unlist(lav_pt_vnames(pts, type = "ov.x"))) == 0L) {
      attr(step1$VETA, "res.slopes") <- NULL
      attr(step1$VETA, "cov.x") <- NULL
      attr(step1$VETA, "mean.x") <- NULL
      lavoptions_pa$conditional.x <- FALSE
      lavoptions_pa$fixed.x <- FALSE
    }
  } else {
    # global SAM

    # the measurement model parameters now become fixed ustart values
    pt_1$ustart[pt_1$free > 0] <- pt_1$est[pt_1$free > 0]

    reg_idx <- lav_pt_subset_sm(
      pt_1 = pt_1,
      idx_only = TRUE
    )

    # remove 'exogenous' factor variances (if any) from reg.idx
    lv_names_x <- lv_names[lv_names %in% unlist(lavpta$vnames$eqs.x) &
      !lv_names %in% unlist(lavpta$vnames$eqs.y)]
    if ((lavoptions_pa$fixed.x || lavoptions_pa$std.lv) &&
        length(lv_names_x) > 0L) {
      var_idx <- which(pt_1$lhs %in% lv_names_x &
        pt_1$op == "~~" &
        pt_1$lhs == pt_1$rhs)
      rm_idx <- which(reg_idx %in% var_idx)
      if (length(rm_idx) > 0L) {
        reg_idx <- reg_idx[-rm_idx]
      }
    }

    # adapt parameter table for structural part
    pts <- pt_1

    # remove constraints we don't need
    con_idx <- which(pts$op %in% c("==", "<", ">", ":="))
    if (length(con_idx) > 0L) {
      needed_idx <- which(con_idx %in% reg_idx)
      if (length(needed_idx) > 0L) {
        con_idx <- con_idx[-needed_idx]
      }
      if (length(con_idx) > 0L) {
        pts <- as.data.frame(pts, stringsAsFactors = FALSE)
        pts <- pts[-con_idx, ]
      }
    }
    pts$est <- NULL
    pts$se <- NULL

    # 'fix' step 1 parameters
    pts$free[!pts$id %in% reg_idx & pts$free > 0L] <- 0L

    # but free up residual variances if fixed (eg std.lv = TRUE) (new in 0.6-20)
    var_idx <- reg_idx[which(pt_1$free[reg_idx] == 0L &
                             pt_1$user[reg_idx] != 1L &
                             pt_1$op[reg_idx] == "~~")] # FIXME: more?
    pts$free[var_idx] <- max(pts$free) + seq_along(var_idx)

    # set 'ustart' values for free FIT.PA parameter to NA
    pts$ustart[pts$free > 0L] <- as.numeric(NA)

    pts <- lav_pt_complete(pts)

    extra_id <- integer(0L)
  } # global

  # fit structural model
  if (lav_verbose()) {
    cat("Fitting the structural part ... \n")
  }
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    if (gamma_flag) {
      nacov <- step1$Gamma.eta
    } else {
      nacov <- NULL
    }
    fit_pa <- tryCatch(
      lavaan::lavaan(pts,
        sample_cov  = step1$VETA,
        sample_mean = step1$EETA,
        sample_nobs = nobs_1,
        nacov       = nacov,
        slot_options = lavoptions_pa,
        verbose     = FALSE
      ),
      error = function(e) e
    )
    if (inherits(fit_pa, "error")) {
      # the corrected two-step STRUCTURAL test (Satorra-Bentler via Gamma.eta)
      # could not be computed for this structural model (eg a bi-factor
      # measurement model, whose VETA jacobian is not conformable here). For
      # se = twostep / twostep.robust / naive the FINAL SEs do not depend on
      # this FIT.PA fit, so degrade gracefully: refit with the standard
      # (uncorrected) structural test. For se = local / local.nt the SEs ARE
      # read from this fit, so we cannot silently degrade -> re-raise.
      if (gamma_flag &&
          lavoptions$se %in% c("twostep", "twostep.robust",
                               "twostep.huber.white", "naive")) {
        lavoptions_pa$test <- "standard"
        fit_pa <- lavaan::lavaan(pts,
          sample_cov  = step1$VETA,
          sample_mean = step1$EETA,
          sample_nobs = nobs_1,
          nacov       = nacov,
          slot_options = lavoptions_pa,
          verbose     = FALSE
        )
        lav_msg_warn(gettext(
          "the two-step corrected structural test could not be computed for
           this model (eg a bi-factor measurement model); the standard
           (uncorrected) structural test is reported instead."))
      } else {
        lav_msg_stop(gettextf(
          "the structural model could not be fitted: %s",
          conditionMessage(fit_pa)))
      }
    }
  } else {
    fit_pa <- lavaan::lavaan(
      model = pts,
      slot_data = fit@Data,
      slot_sample_stats = fit@SampleStats,
      slot_options = lavoptions_pa,
      verbose = FALSE
    )
  }
  if (lav_verbose()) {
    cat("Fitting the structural part ... done.\n")
  }

  # check that the structural part is identified from the latent moments
  # alone: in the SAM approach the structural model is estimated from the
  # (estimated) latent variable moments, so it cannot borrow identification
  # from the measurement part (eg a non-recursive system without
  # instruments may 'fit' in sem() but has more structural parameters than
  # latent moments). Without this check the step-2 vcov machinery fails
  # cryptically further down.
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    pa_df <- fit_pa@test[[1]]$df
    if (!is.null(pa_df) && !is.na(pa_df) && pa_df < 0L) {
      lav_msg_stop(gettextf(
        "the structural part of the model is not identified: it has more
         free parameters than there are (estimated) latent variable moments
         (df = %d). In the SAM approach the structural model must be
         identified from the latent variable moments alone. Consider
         simplifying the structural part, or using sem() instead.", pa_df))
    }
  }

  # which parameters from PTS do we wish to fill in:
  # - all 'free' parameters
  # - :=, <, > (if any)
  # - and NOT element with user=3 (add.exo.cov = TRUE, extra.int.idx)
  pts_idx <- which((pts$free > 0L | (pts$op %in% c(":=", "<", ">"))) &
    !pts$user == 3L)

  # find corresponding rows in PT
  pts2 <- as.data.frame(pts, stringsAsFactors = FALSE)
  pt_idx <- lav_pt_map_id_p1_in_p2(pts2[pts_idx, ], pt_1,
    exclude_nonpar = FALSE
  )
  # fill in
  pt_1$est[pt_idx] <- fit_pa@ParTable$est[pts_idx]

  # create step2.free.idx
  p2_idx <- seq_along(pt_1$lhs) %in% pt_idx & pt_1$free > 0 # no def!
  step2_free_idx <- step1$PT.free[p2_idx]

  # add 'step' column in PT
  pt_1$step <- rep(1L, length(pt_1$lhs))
  pt_1$step[seq_along(pt_1$lhs) %in% reg_idx] <- 2L

  step2 <- list(
    FIT.PA = fit_pa, PT = pt_1, reg.idx = reg_idx,
    step2.free.idx = step2_free_idx, extra.id = extra_id,
    pt.idx = pt_idx, pts.idx = pts_idx
  )

  step2
}
