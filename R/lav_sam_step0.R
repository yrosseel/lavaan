# STEP 0: process full model, without fitting
lav_sam_step0 <- function(cmd = "sem", model = NULL, data = NULL,
                          se = "twostep", sam_method = "local",
                          dotdotdot = NULL) {

  # create dotdotdot0 for dummy fit
  dotdotdot0 <- dotdotdot

  # 'parse' the model, so we can inspect a few features; accept the same
  # model= input types as a regular call to lavaan(): model syntax (character),
  # a formula, a fitted lavaan object, a flat model, or a parameter table
  # (issue #514)
  if (is.null(dotdotdot$parser)) {
    useparser <- lav_options_default()$parser
  } else {
    useparser <- dotdotdot$parser
  }
  flat_model <- lav_step01_ovnames_initflat(
    model = model, dotdotdot_parser = useparser
  )

  # composites (<~) have no reflective measurement model, so there is no
  # measurement block that step 1 could estimate in isolation
  if (any(flat_model$op == "<~")) {
    lav_msg_stop(gettext(
      "composites (defined by the <~ operator) are not supported by the
       SAM approach (yet); consider using sem() instead."))
  }

  # remove do.fit option if present
  dotdotdot0$do.fit <- NULL

  if (sam_method %in% c("local", "fsr", "cfsr")) {
    dotdotdot0$sample.icov <- FALSE # if N < nvar
  }
  dotdotdot0$se                   <- "none"
  dotdotdot0$test                 <- "none"
  dotdotdot0$verbose              <- FALSE # no output for this 'dummy' FIT
  dotdotdot0$ceq.simple           <- TRUE # if not the default yet
  dotdotdot0$check.lv.interaction <- FALSE # we allow for it
  # note: setting cat.wls.w = FALSE (no weight matrix if categorical)
  # would break the computation of twostep standard errors

  # single-level clustered data (cluster = , no level: syntax): the classic
  # twostep correction uses the model-based (expected) information for the
  # second-step piece, which is not robust to clustering; switch to the
  # robust (sandwich) variant, which picks up the cluster-robust Gamma
  multilevel_flag <- any(flat_model$op == ":" &
                         tolower(flat_model$lhs) == "level")
  # if model= was a parameter table, the "level:" rows are absent, but the
  # level column is available
  if (!multilevel_flag && !is.null(flat_model$level)) {
    multilevel_flag <-
      length(unique(flat_model$level[flat_model$level > 0L])) > 1L
  }
  if (se == "twostep" && !multilevel_flag && !is.null(dotdotdot$cluster)) {
    se <- "twostep.robust"
  }

  if (se %in% c("local", "ij", "twostep.robust")) {
    dotdotdot0$sample.icov <- TRUE
    # keep the moments of the exogenous covariates in the sample statistics
    # (their sampling variability is part of Gamma) -- but not when
    # conditional.x = TRUE, where the x-moments are fixed by design and
    # fixed.x = FALSE is not supported (see also the retry below, for the
    # settings where conditional.x = TRUE only emerges during the options
    # processing, eg categorical data + exogenous covariates)
    if (!isTRUE(dotdotdot$conditional.x)) {
      dotdotdot0$fixed.x <- FALSE
    }
    # ij/twostep.robust need the stored observed Gamma (NACOV); the
    # "force.model" ov_order sentinel keeps it in model order (it is only
    # honored when NACOV is supplied; see lav_lavaan_step00_init()). se =
    # "local" builds Gamma.eta in (memory-lean) influence form and computes the
    # observed Gamma on demand only when needed (eg categorical), so it needs
    # neither -- and the default ov_order is already "model".
    if (se %in% c("ij", "twostep.robust")) {
      dotdotdot0$NACOV <- TRUE
      dotdotdot0$ov_order <- "force.model" # avoid data ordering...
    }
  }

  # any lv interaction terms?
  if (length(lav_object_vnames(flat_model, "lv.interaction")) > 0L) {
    dotdotdot0$meanstructure <- TRUE
  }

  # initial processing of the model, no fitting
  fit <- tryCatch(
    do.call(cmd,
      args = c(list(
        model  = flat_model,
        data   = data,
        do.fit = FALSE
      ), dotdotdot0)
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    fit0_error <- fit
    fit <- NULL
    if (identical(dotdotdot0$fixed.x, FALSE) ||
        isTRUE(dotdotdot0$gamma.unbiased)) {
      # fixed.x = FALSE (set above for se = local/ij/twostep.robust) and the
      # unbiased Gamma (the sam() default) are both unsupported when
      # conditional.x = TRUE, which may only emerge during the options
      # processing (eg categorical data + exogenous covariates, where
      # conditional.x = TRUE is the default): retry with
      # conditional.x-compatible settings, and keep the result only if
      # conditional.x is indeed TRUE
      dotdotdot0$fixed.x <- NULL
      dotdotdot0$gamma.unbiased <- FALSE
      fit <- tryCatch(
        do.call(cmd,
          args = c(list(
            model  = flat_model,
            data   = data,
            do.fit = FALSE
          ), dotdotdot0)
        ),
        error = function(e) NULL
      )
      if (!is.null(fit) && !fit@Model@conditional.x) {
        fit <- NULL # not the conditional.x situation: keep the original error
      }
    }
    if (is.null(fit)) {
      lav_msg_stop(gettextf(
        "initial processing of the model failed: %s",
        conditionMessage(fit0_error)))
    }
  }

  # restore options

  # do.fit
  fit@Options$do.fit <- TRUE

  # conditional.x: the unbiased Gamma does not exist -- make sure no
  # downstream fit (measurement blocks, structural part, lav_gamma_used())
  # ever requests it (the sam() default is gamma.unbiased = TRUE)
  if (fit@Model@conditional.x && isTRUE(fit@Options$gamma.unbiased)) {
    fit@Options$gamma.unbiased <- FALSE
  }

  # sample.icov
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    fit@Options$sample.icov <- TRUE
  }

  # se
  # PML: the two-step (robust) corrections need ingredients (eg the h1
  # information, the block influence P) that are not available for the
  # pairwise likelihood; only the naive (structural-fit) standard errors
  # can be computed, and only for the local methods
  if (identical(fit@Options$estimator.orig, "PML")) {
    if (sam_method == "global") {
      lav_msg_stop(gettext(
        "sam(sam.method = \"global\") does not support estimator PML (yet);
         use sam.method = \"local\" or sem() instead."))
    }
    # notes:
    # - se = "local" IS available for PML: it uses the casewise
    #   structural-space Gamma.eta (see lav_sam_gamma_eta_pml(); single
    #   group, complete data, all-ordinal -- guarded there);
    # - "twostep.huber.white" is in the fallback list DELIBERATELY, even
    #   though casewise PML scores exist (lav_sc_pml()): for
    #   sam.method = "local" the joint-PML-given-theta1 estimator that the
    #   huber.white sandwich linearizes is NOT the local (VETA-fit) step-2
    #   estimator (their variance estimates differ substantially), so the
    #   sandwich describes the wrong estimator; use se = "local" instead.
    if (se %in% c("twostep", "twostep.robust", "twostep.huber.white",
                  "local.nt")) {
      lav_msg_warn(gettextf(
        "se = \"%s\" is not available (yet) for estimator PML; naive
         standard errors are reported instead
         (tip: se = \"local\" is available for PML).", se))
      se <- "naive"
    }
  }
  if (fit@Model@categorical && se == "twostep") {
    # for categorical data, the classic ('global') two-step correction uses
    # the model-based information matrix, which underestimates the standard
    # errors for the (D)WLS estimator. Use the robust (Yuan & Chan, 2002)
    # correction instead. This needs the 'P' matrix (the influence of the
    # measurement parameters on the sample statistics), which is available for
    # all sam methods (the measurement blocks are always fitted in step 1).
    se <- "twostep.robust"
  }
  if (fit@Model@conditional.x && se == "twostep" &&
      sam_method %in% c("local", "fsr", "cfsr")) {
    # conditional.x + local sam method: the classic two-step correction is
    # built on the joint information matrix, so it linearizes the
    # joint-model-given-theta1 estimator. Under conditional.x that estimator
    # is NOT asymptotically equivalent to the local step-2 estimator (the
    # structural fit is saturated in the latent-on-x slopes and cannot pool
    # information across statistics the way the joint estimator does), and
    # the joint-information formula underestimates the sampling variability
    # (badly so for binary covariates). Use the robust variant instead: under
    # conditional.x it is computed as the structural-space sandwich (see the
    # tsrobust_condx_flag note in lav_sam_step2_se.R), which linearizes the
    # actual two-step estimator (and coincides with se = "local"). For
    # sam.method = "global" step 2 IS the joint estimator, and the classic
    # formula still applies.
    se <- "twostep.robust"
  }
  fit@Options$se <- se

  # test
  if (!is.null(dotdotdot$test)) {
    fit@Options$test <- lav_test_rename(dotdotdot$test)
  } else if (sam_method == "global") {
    # global SAM: the Yuan & Chan (2002) rescaled GLOBAL statistic is the
    # default robust test (the analogue of test = "satorra.bentler" in sem()).
    fit@Options$test <- "yuan.chan"
  } else {
    fit@Options$test <- "standard"
  }

  # adjust parameter table:
  pt_1 <- fit@ParTable

  # check parameter table
  pt_1$est <- pt_1$se <- NULL
  # est equals ustart by default (except exo values)
  pt_1$est <- pt_1$ustart
  if (any(pt_1$exo > 0L)) {
    pt_1$est[pt_1$exo > 0L] <- pt_1$start[pt_1$exo > 0L]
  }
  # thresholds: fill in the sample-based starting values. Thresholds of
  # observed (dummy) variables that are not part of any measurement block
  # (eg an ordered endogenous variable in the structural part) are never
  # (re)estimated in step 1 or step 2 and keep these h1 values; leaving
  # them NA breaks every computation that evaluates the joint model
  # matrices (information, Gamma.eta) before the final assembly. Block
  # thresholds are simply overwritten by the step-1 estimates.
  th_idx <- which(pt_1$op == "|" & is.na(pt_1$est))
  if (length(th_idx) > 0L) {
    pt_1$est[th_idx] <- pt_1$start[th_idx]
  }

  # clear se values (needed here?) only for global approach to compute SE
  pt_1$se <- rep(as.numeric(NA), length(pt_1$lhs))
  pt_1$se[pt_1$free == 0L & !is.na(pt_1$ustart)] <- 0.0

  fit@ParTable <- pt_1


  fit
}
