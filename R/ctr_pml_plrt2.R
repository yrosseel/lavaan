ctr_pml_plrt2 <- function(lavobject = NULL, lavmodel = NULL, lavdata = NULL,
                          lavsamplestats = NULL, lavpartable = NULL,
                          lavoptions = NULL, x = NULL, VCOV = NULL,
                          lavcache = NULL) {
  lavpta <- NULL
  if (!is.null(lavobject)) {
    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
    lavsamplestats <- lavobject@SampleStats
    lavcache <- lavobject@Cache
    lavpartable <- lav_partable_set_cache(lavobject@ParTable, lavobject@pta)
    lavpta <- lavobject@pta
  }
  if (is.null(lavpta)) {
    lavpta <- lav_partable_attributes(lavpartable)
    lavpartable <- lav_partable_set_cache(lavpartable, lavpta)
  }

  if (is.null(x)) {
    # compute 'fx' = objective function value
    # (NOTE: since 0.5-18, NOT divided by N!!)
    fx <- lav_model_objective(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache
    )
    H0.fx <- as.numeric(fx)
    H0.fx.group <- attr(fx, "fx.group")
  } else {
    H0.fx <- attr(attr(x, "fx"), "fx.pml")
    H0.fx.group <- attr(attr(x, "fx"), "fx.group")
  }

  # fit a saturated model 'fittedSat'
  ModelSat <- lav_partable_unrestricted(
    lavobject = NULL,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavpta = lavpta,
    lavsamplestats = lavsamplestats
  )

  # FIXME: se="none", test="none"??
  Options <- lavoptions
  Options$verbose <- FALSE
  Options$se <- "none"
  Options$test <- "none"
  Options$baseline <- FALSE
  Options$h1 <- FALSE
  fittedSat <- lavaan(ModelSat,
    slotOptions = Options,
    slotSampleStats = lavsamplestats,
    slotData = lavdata, slotCache = lavcache
  )
  fx <- lav_model_objective(
    lavmodel = fittedSat@Model,
    lavsamplestats = fittedSat@SampleStats,
    lavdata = fittedSat@Data,
    lavcache = fittedSat@Cache
  )
  SAT.fx <- as.numeric(fx)
  SAT.fx.group <- attr(fx, "fx.group")

  # we also need a `saturated model', but where the moments are based
  # on the model-implied sample statistics under H0
  ModelSat2 <-
    lav_partable_unrestricted(
      lavobject = NULL,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavpta = lavpta,
      lavsamplestats = NULL,
      sample.cov = computeSigmaHat(lavmodel),
      sample.mean = computeMuHat(lavmodel),
      sample.th = computeTH(lavmodel),
      sample.th.idx = lavsamplestats@th.idx
    )

  Options2 <- Options
  Options2$optim.method <- "none"
  Options2$optim.force.converged <- TRUE
  fittedSat2 <- lavaan(ModelSat2,
    slotOptions = Options2,
    slotSampleStats = lavsamplestats,
    slotData = lavdata, slotCache = lavcache
  )

  # the code below was contributed by Myrsini Katsikatsou (Jan 2015)

  # for now, only a single group is supported:
  # g = 1L


  ########################### The code for PLRT for overall goodness of fit

  ##### Section 1. Compute the asymptotic mean and variance
  #####            of the first quadratic quantity
  # if(is.null(VCOV)) {
  #    VCOV <- lav_model_vcov(lavmodel       = lavmodel,
  #                           lavsamplestats = lavsamplestats,
  #                           lavoptions     = lavoptions,
  #                           lavdata        = lavdata,
  #                           lavpartable    = lavpartable,
  #                           lavcache       = lavcache)
  # }
  # G.inv
  # InvG_attheta0 <- lavsamplestats@ntotal * VCOV[,]
  # Hessian
  # H_attheta0 <- solve(attr(VCOV, "E.inv"))

  # inverted observed information ('H.inv')
  if (is.null(VCOV)) {
    H0.inv <- lav_model_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavcache = lavcache, augmented = TRUE, inverted = TRUE
    )
  } else {
    H0.inv <- attr(VCOV, "E.inv")
  }

  # first order information ('J')
  if (is.null(VCOV)) {
    J0 <- lav_model_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavcache = lavcache
    )[, ]
  } else {
    # we do not get J, but J.group, FIXME?
    J0 <- lav_model_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavcache = lavcache
    )[, ]
  }

  # inverted Godambe information
  G0.inv <- H0.inv %*% J0 %*% H0.inv

  H0tmp_prod1 <- H0.inv %*% J0
  # H0tmp_prod1 <- InvG_attheta0 %*% H_attheta0
  H0tmp_prod2 <- H0tmp_prod1 %*% H0tmp_prod1
  E_tww <- sum(diag(H0tmp_prod1))
  var_tww <- 2 * sum(diag(H0tmp_prod2))

  ##### Section 2: Compute the asymptotic mean and variance
  #####            of the second quadratic quantity.
  tmp.options <- fittedSat2@Options
  tmp.options$se <- "robust.huber.white"
  VCOV.Sat2 <- lav_model_vcov(
    lavmodel = fittedSat2@Model,
    lavsamplestats = fittedSat2@SampleStats,
    lavoptions = tmp.options,
    lavdata = fittedSat2@Data,
    lavpartable = fittedSat2@ParTable,
    lavcache = fittedSat2@Cache
  )
  # G.inv at vartheta_0
  InvG_at_vartheta0 <- lavsamplestats@ntotal * VCOV.Sat2[, ]
  # Hessian at vartheta_0
  H_at_vartheta0 <- solve(attr(VCOV.Sat2, "E.inv")) # should always work
  # H1.inv <- lavTech(fittedSat2, "inverted.information.observed")
  # J1     <- lavTech(fittedSat2, "information.first.order")
  # H1tmp_prod1 <- H1.inv %*% J1
  H1tmp_prod1 <- InvG_at_vartheta0 %*% H_at_vartheta0
  H1tmp_prod2 <- H1tmp_prod1 %*% H1tmp_prod1
  E_tzz <- sum(diag(H1tmp_prod1))
  var_tzz <- 2 * sum(diag(H1tmp_prod2))


  ##### Section 3: Compute the asymptotic covariance
  #####            of the two quadratic quantities

  drhodpsi_MAT <- vector("list", length = lavsamplestats@ngroups)
  group.values <- lav_partable_group_values(fittedSat2@ParTable)
  for (g in 1:lavsamplestats@ngroups) {
    delta.g <- computeDelta(lavmodel)[[g]]
    # order of the rows: first the thresholds, then the correlations
    # we need to map the rows of delta.g to the rows/cols of H_at_vartheta0
    # of H1

    PT <- fittedSat2@ParTable
    PT$label <- lav_partable_labels(PT)
    free.idx <- which(PT$free > 0 & PT$group == group.values[g])
    PARLABEL <- PT$label[free.idx]

    # for now, we can assume that computeDelta will always return
    # the thresholds first, then the correlations
    #
    # later, we should add a (working) add.labels = TRUE option to
    # computeDelta
    th.names <- lavpta$vnames$th[[g]]
    ov.names <- lavpta$vnames$ov[[g]]
    tmp <- utils::combn(ov.names, 2)
    cor.names <- paste(tmp[1, ], "~~", tmp[2, ], sep = "")
    NAMES <- c(th.names, cor.names)
    if (g > 1L) {
      NAMES <- paste(NAMES, ".g", g, sep = "")
    }

    par.idx <- match(PARLABEL, NAMES)
    drhodpsi_MAT[[g]] <- delta.g[par.idx, , drop = FALSE]
  }
  drhodpsi_mat <- do.call(rbind, drhodpsi_MAT)

  # tmp_prod <- ( t(drhodpsi_mat) %*% H_at_vartheta0 %*%
  #                 drhodpsi_mat %*% InvG_attheta0 %*%
  #                 H_attheta0 %*% InvG_attheta0 )
  tmp_prod <- (t(drhodpsi_mat) %*% H_at_vartheta0 %*%
    drhodpsi_mat %*% H0.inv %*% J0 %*% G0.inv)
  cov_tzztww <- 2 * sum(diag(tmp_prod))

  ##### Section 4: compute the adjusted PLRT and its p-value
  PLRTH0Sat <- 2 * (H0.fx - SAT.fx)
  PLRTH0Sat.group <- 2 * (H0.fx.group - SAT.fx.group)
  asym_mean_PLRTH0Sat <- E_tzz - E_tww
  asym_var_PLRTH0Sat <- var_tzz + var_tww - 2 * cov_tzztww
  scaling.factor <- (asym_mean_PLRTH0Sat / (asym_var_PLRTH0Sat / 2))
  FSA_PLRT_SEM <- (asym_mean_PLRTH0Sat / (asym_var_PLRTH0Sat / 2)) * PLRTH0Sat
  adjusted_df <- (asym_mean_PLRTH0Sat * asym_mean_PLRTH0Sat) /
    (asym_var_PLRTH0Sat / 2)
  # In some very few cases (simulations show very few cases
  #                         in small sample sizes)
  # the adjusted_df is a negative number, we should then
  # print a warning like: "The adjusted df is computed to be a negative number
  # and for this the first and second moment adjusted PLRT is not computed." .
  pvalue <- 1 - pchisq(FSA_PLRT_SEM, df = adjusted_df)

  list(
    PLRTH0Sat = PLRTH0Sat, PLRTH0Sat.group = PLRTH0Sat.group,
    stat = FSA_PLRT_SEM, df = adjusted_df, p.value = pvalue,
    scaling.factor = scaling.factor
  )
}
############################################################################
