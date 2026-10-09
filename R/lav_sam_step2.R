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
  # The user may ask for another member of the Satorra-Bentler family instead
  # (test = "mean.var.adjusted", "scaled.shifted", or their ".corrected"
  # versions, see lav_test_hayakawa.R): the same Gamma.eta then feeds that
  # adjustment of the structural test.
  if (gamma_flag) {
    sb_tests <- lavoptions_pa$test[lavoptions_pa$test %in% lav_sam_sb_family]
    if (length(sb_tests) == 0L) {
      sb_tests <- "satorra.bentler"
    }
    # a non-standard base statistic (scaled.test = "browne.residual.nt.model",
    # Hayakawa's RLS version) must stay in front, as lav_options_set() does
    scaled_base <- lavoptions_pa$scaled.test
    if (!is.null(scaled_base) && !scaled_base %in% c("standard", "default")) {
      sb_tests <- unique(c(scaled_base, sb_tests))
    }
    lavoptions_pa$test <- sb_tests
  } else if (sam_method %in% c("local", "fsr", "cfsr") &&
             any(lavoptions_pa$test %in% lav_sam_sb_family)) {
    # no Gamma.eta (eg se = "none"/"bootstrap"): the moments-only structural
    # fit has no NACOV, so none of the scaled tests can be computed
    lav_msg_warn(gettextf(
      "the requested test (%s) needs Gamma.eta, which is not available
       with se = %s; the standard structural test is reported instead.",
      lav_msg_view(lavoptions_pa$test[lavoptions_pa$test %in%
                                      lav_sam_sb_family]),
      dQuote(lavoptions_pa$se, q = FALSE)))
    lavoptions_pa$test <- "standard"
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
    # FIT.PA's vcov stays the NAIVE (standard) one: the alpha_correction blend
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
  if (!is.null(struc_args$test)) {
    lavoptions_pa$test <- lav_test_rename(struc_args$test)
  }

  # the corrected adjusted tests (Hayakawa 2018) need the casewise rows
  # behind Gamma.eta, which the moments-only structural fit does not have:
  # take them out of the structural fit, and add them afterwards (see below)
  corrected_tests <- character(0L)
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    corrected_tests <-
      lavoptions_pa$test[lavoptions_pa$test %in% lav_sam_corrected_tests]
    if (length(corrected_tests) > 0L) {
      # keep the order of the requested tests for the final @test slot
      requested_tests <- lavoptions_pa$test
      lavoptions_pa$test <- setdiff(lavoptions_pa$test, corrected_tests)
      if (length(lavoptions_pa$test) == 0L) {
        lavoptions_pa$test <- "satorra.bentler"
      }
    }
  }
  # information.meat.hc: the leverage adjustment of the local SEs is
  # applied afterwards, in lav_sam_step2_se() (the structural fit has no
  # casewise data of its own); the structural fit keeps the classic meat
  lavoptions_pa$information.meat.hc <- "HC0"

  # new in 0.7-2: the bread of the structural sandwich behind the local
  # standard errors uses the OBSERVED information by default: whenever the
  # structural model constrains the moments of an equation's own predictors
  # (eg two endogenous predictors without a residual covariance), the
  # expected-information bread underestimates the sampling variability under
  # misspecification, while the observed (hessian) bread remains correct;
  # when every equation reproduces its own predictor moments (eg all
  # predictors exogenous), both breads coincide and nothing changes. The
  # classic behavior remains available via
  # struc_args = list(information.bread = "default") (= follow the
  # information option) or "expected". The twostep.robust + conditional.x
  # reroute reads the same FIT.PA sandwich (see tsrobust_condx_flag in
  # lav_sam_step2_se()) and must stay identical to se = "local".
  info_bread_pa <- lavoptions_pa$information.bread
  if (is.null(info_bread_pa)) {
    info_bread_pa <- "default"
  }
  if (info_bread_pa == "default" &&
      is.null(struc_args[["information.bread"]]) &&
      lavoptions_pa$estimator == "ML" &&
      sam_method %in% c("local", "fsr", "cfsr") &&
      (lavoptions$se %in% c("local", "local.nt") ||
       (lavoptions$se == "twostep.robust" && lavoptions_pa$conditional.x &&
        gamma_flag))) {
    lavoptions_pa$information.bread <- "observed"
  }

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
    # reset any stale bounds (a fixed parameter has lower == its fixed
    # value whenever the partable carries bounds columns; see the same
    # fix in lav_pt_subset_sm())
    if (!is.null(pts$lower)) {
      pts$lower[var_idx] <- -Inf
    }
    if (!is.null(pts$upper)) {
      pts$upper[var_idx] <- +Inf
    }

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

  # the corrected adjusted structural tests (Hayakawa 2018): the unbiased
  # estimator of tr(UGamma^2) needs the casewise rows of Gamma.eta (the
  # influence contributions of the cases to the latent moments, built in
  # step 1c), which take the place of the raw data (see lav_test_hayakawa.R).
  # Computed afterwards, so that the structural fit itself stays the regular
  # (moments-only) fit; the resulting entries are added to FIT.PA@test in
  # the order the tests were requested.
  if (length(corrected_tests) > 0L) {
    fit_pa <- lav_sam_step2_corrected_test(
      fit_pa = fit_pa, step1 = step1,
      corrected_tests = corrected_tests, requested_tests = requested_tests
    )
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

# the Satorra-Bentler family of structural tests a local SAM model can
# report (all driven by Gamma.eta as the NACOV of the latent moments)
lav_sam_sb_family <- c(
  "satorra.bentler", "mean.var.adjusted", "scaled.shifted",
  "mean.var.adjusted.corrected", "scaled.shifted.corrected"
)
# the members that need the casewise rows of Gamma.eta (Hayakawa 2018)
lav_sam_corrected_tests <- c(
  "mean.var.adjusted.corrected", "scaled.shifted.corrected"
)

# is a corrected adjusted test requested, either through test = (already
# canonical) or through struc_args = list(test = )?
lav_sam_corrected_test_flag <- function(test = NULL, struc_test = NULL) {
  if (!is.null(struc_test)) {
    struc_test <- lav_test_rename(struc_test)
  }
  any(c(test, struc_test) %in% lav_sam_corrected_tests)
}

# add the corrected adjusted structural tests to FIT.PA (single group only):
# the unbiased tr(UGamma^2) estimator reads the casewise rows of Gamma.eta
# (step1$Gamma.eta.rows). When the rows are not available for this setting,
# the corrected tests are dropped with a warning (the other requested tests
# are reported as usual).
lav_sam_step2_corrected_test <- function(fit_pa = NULL, step1 = NULL,
                                         corrected_tests = character(0L),
                                         requested_tests = character(0L)) {
  rows <- NULL
  if (length(step1$Gamma.eta.rows) > 0L) {
    rows <- step1$Gamma.eta.rows[[1]]
  }
  test_c <- NULL
  if (is.null(rows) || !is.matrix(rows)) {
    lav_msg_warn(gettextf(
      "the corrected adjusted tests (%s) need the casewise contributions
       to Gamma.eta, which are not available for this model (eg se =
       \"local.nt\", or a setting without casewise influence); they are not
       reported.", paste(dQuote(corrected_tests, q = FALSE),
                         collapse = ", ")))
  } else {
    test_c <- tryCatch(
      lav_test_sb(lavobject = fit_pa, test = corrected_tests,
                  gamma_rows = rows),
      error = function(e) e
    )
    if (inherits(test_c, "error")) {
      lav_msg_warn(gettextf(
        "the corrected adjusted tests (%1$s) could not be computed for the
         structural model: %2$s",
        paste(dQuote(corrected_tests, q = FALSE), collapse = ", "),
        conditionMessage(test_c)))
      test_c <- NULL
    }
  }
  if (is.null(test_c)) {
    return(fit_pa)
  }

  # @test: "standard" first, then the requested tests in their order
  keep <- unique(c("standard", requested_tests))
  opts_pa <- fit_pa@Options # the (moments-only) options of the fit
  fit_pa@test <- lav_sam_test_merge(fit_pa@test, test_c, keep)
  fit_pa@Options$test <- names(fit_pa@test)

  # the baseline (independence) model too, so that the scaled fit indices
  # (cfi.scaled, ...) based on the corrected test are available: refit it on
  # the latent moments (as lav_sam_struc_fit_object() does) and apply the
  # same correction
  if (!is.null(fit_pa@baseline$partable) && !is.null(fit_pa@baseline$test)) {
    test_base <- tryCatch({
      meanstr <- fit_pa@Model@meanstructure
      fit_base <- lavaan::lavaan(
        model = fit_pa@baseline$partable,
        sample.cov = step1$VETA,
        sample.mean = if (meanstr) step1$EETA else NULL,
        sample.nobs = as.list(unlist(fit_pa@SampleStats@nobs)),
        nacov = step1$Gamma.eta,
        slot_options = opts_pa
      )
      lav_test_sb(lavobject = fit_base, test = corrected_tests,
                  gamma_rows = rows)
    }, error = function(e) NULL)
    if (!is.null(test_base)) {
      fit_pa@baseline$test <-
        lav_sam_test_merge(fit_pa@baseline$test, test_base, keep)
    }
  }
  fit_pa
}

# merge the entries of test_new (by name) into test_old, and return the
# entries named in keep (in that order), ignoring names that are absent
lav_sam_test_merge <- function(test_old = list(), test_new = list(),
                               keep = character(0L)) {
  for (nm in setdiff(names(test_new), "standard")) {
    test_old[[nm]] <- test_new[[nm]]
  }
  test_old[keep[keep %in% names(test_old)]]
}
