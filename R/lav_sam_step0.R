# STEP 0: process full model, without fitting
lav_sam_step0 <- function(cmd = "sem", model = NULL, data = NULL,
                          se = "twostep", sam_method = "local",
                          dotdotdot = NULL) {

  # create dotdotdot0 for dummy fit
  dotdotdot0 <- dotdotdot

  # parse model, so we can inspect a few features
  flat_model <- lavParseModelString(model)

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
  if (se == "twostep" && !multilevel_flag && !is.null(dotdotdot$cluster)) {
    se <- "twostep.robust"
  }

  if (se %in% c("local", "ij", "twostep.robust")) {
    dotdotdot0$sample.icov <- TRUE
    dotdotdot0$fixed.x <- FALSE
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
  fit <- do.call(cmd,
    args = c(list(
      model  = flat_model,
      data   = data,
      do.fit = FALSE
    ), dotdotdot0)
  )

  # restore options

  # do.fit
  fit@Options$do.fit <- TRUE

  # sample.icov
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    fit@Options$sample.icov <- TRUE
  }

  # se
  if (fit@Model@categorical && se == "twostep") {
    # for categorical data, the classic ('global') two-step correction uses
    # the model-based information matrix, which underestimates the standard
    # errors for the (D)WLS estimator. Use the robust (Yuan & Chan, 2002)
    # correction instead. This needs the 'P' matrix (the influence of the
    # measurement parameters on the sample statistics), which is available for
    # all sam methods (the measurement blocks are always fitted in step 1).
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

  # clear se values (needed here?) only for global approach to compute SE
  pt_1$se <- rep(as.numeric(NA), length(pt_1$lhs))
  pt_1$se[pt_1$free == 0L & !is.na(pt_1$ustart)] <- 0.0

  fit@ParTable <- pt_1


  fit
}
