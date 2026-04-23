lav_pml_plrt2 <- function(lavobject = NULL, lavmodel = NULL, lavdata = NULL,
                          lavsamplestats = NULL, lavpartable = NULL,
                          lavoptions = NULL, x = NULL, vcov_1 = NULL,
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
    h0_fx <- as.numeric(fx)
    h0_fx_group <- attr(fx, "fx.group")
  } else {
    h0_fx <- attr(attr(x, "fx"), "fx.pml")
    h0_fx_group <- attr(attr(x, "fx"), "fx.group")
  }

  # fit a saturated model 'fittedSat'
  model_sat <- lav_partable_unrestricted(
    lavobject = NULL,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavpta = lavpta,
    lavsamplestats = lavsamplestats
  )

  # FIXME: se="none", test="none"??
  options_1 <- lavoptions
  options_1$se <- "none"
  options_1$test <- "none"
  options_1$baseline <- FALSE
  options_1$h1 <- FALSE
  fitted_sat <- lavaan(model_sat,
    slotOptions = options_1, verbose = FALSE,
    slotSampleStats = lavsamplestats,
    slotData = lavdata, slotCache = lavcache
  )
  fx <- lav_model_objective(
    lavmodel = fitted_sat@Model,
    lavsamplestats = fitted_sat@SampleStats,
    lavdata = fitted_sat@Data,
    lavcache = fitted_sat@Cache
  )
  sat_fx <- as.numeric(fx)
  sat_fx_group <- attr(fx, "fx.group")

  # we also need a `saturated model', but where the moments are based
  # on the model-implied sample statistics under H0
  model_sat2 <-
    lav_partable_unrestricted(
      lavobject = NULL,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavpta = lavpta,
      lavsamplestats = NULL,
      sample.cov = lav_model_sigma(lavmodel),
      sample.mean = lav_model_mu(lavmodel),
      sample.th = lav_model_th(lavmodel),
      sample.th.idx = lavsamplestats@th.idx
    )

  options2 <- options_1
  options2$optim.method <- "none"
  options2$optim.force.converged <- TRUE
  fitted_sat2 <- lavaan(model_sat2,
    slotOptions = options2, verbose = FALSE,
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
  if (is.null(vcov_1)) {
    h0_inv <- lav_model_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavoptions = lavoptions,
      lavcache = lavcache, augmented = TRUE, inverted = TRUE
    )
  } else {
    h0_inv <- attr(vcov_1, "E.inv")
  }

  # first order information ('J')
  if (is.null(vcov_1)) {
    j0 <- lav_model_information_firstorder(
      lavmodel = lavmodel, lavoptions = lavoptions,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavcache = lavcache
    )[, ]
  } else {
    # we do not get J, but J.group, FIXME?
    j0 <- lav_model_information_firstorder(
      lavmodel = lavmodel, lavoptions = lavoptions,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavcache = lavcache
    )[, ]
  }

  # inverted Godambe information
  g0_inv <- h0_inv %*% j0 %*% h0_inv

  h0tmp_prod1 <- h0_inv %*% j0
  # H0tmp_prod1 <- InvG_attheta0 %*% H_attheta0
  h0tmp_prod2 <- h0tmp_prod1 %*% h0tmp_prod1
  e_tww <- sum(diag(h0tmp_prod1))
  var_tww <- 2 * sum(diag(h0tmp_prod2))

  ##### Section 2: Compute the asymptotic mean and variance
  #####            of the second quadratic quantity.
  tmp_options <- fitted_sat2@Options
  tmp_options$se <- "robust.huber.white"
  vcov_sat2 <- lav_model_vcov(
    lavmodel = fitted_sat2@Model,
    lavsamplestats = fitted_sat2@SampleStats,
    lavoptions = tmp_options,
    lavdata = fitted_sat2@Data,
    lavpartable = fitted_sat2@ParTable,
    lavcache = fitted_sat2@Cache
  )
  # G.inv at vartheta_0
  inv_g_at_vartheta0 <- lavsamplestats@ntotal * vcov_sat2[, ]
  # Hessian at vartheta_0
  h_at_vartheta0 <- solve(attr(vcov_sat2, "E.inv")) # should always work
  # H1.inv <- lavTech(fittedSat2, "inverted.information.observed")
  # J1     <- lavTech(fittedSat2, "information.first.order")
  # H1tmp_prod1 <- H1.inv %*% J1
  h1tmp_prod1 <- inv_g_at_vartheta0 %*% h_at_vartheta0
  h1tmp_prod2 <- h1tmp_prod1 %*% h1tmp_prod1
  e_tzz <- sum(diag(h1tmp_prod1))
  var_tzz <- 2 * sum(diag(h1tmp_prod2))


  ##### Section 3: Compute the asymptotic covariance
  #####            of the two quadratic quantities

  drhodpsi_mat_1 <- vector("list", length = lavsamplestats@ngroups)
  group_values <- lav_partable_group_values(fitted_sat2@ParTable)
  for (g in 1:lavsamplestats@ngroups) {
    delta_g <- lav_model_delta(lavmodel)[[g]]
    # order of the rows: first the thresholds, then the correlations
    # we need to map the rows of delta.g to the rows/cols of H_at_vartheta0
    # of H1

    pt_1 <- fitted_sat2@ParTable
    pt_1$label <- lav_partable_labels(pt_1)
    free_idx <- which(pt_1$free > 0 & pt_1$group == group_values[g])
    parlabel <- pt_1$label[free_idx]

    # for now, we can assume that lav_model_delta will always return
    # the thresholds first, then the correlations
    #
    # later, we should add a (working) add.labels = TRUE option to
    # lav_model_delta
    th_names <- lavpta$vnames$th[[g]]
    ov_names <- lavpta$vnames$ov[[g]]
    tmp <- utils::combn(ov_names, 2)
    cor_names <- paste(tmp[1, ], "~~", tmp[2, ], sep = "")
    names_1 <- c(th_names, cor_names)
    if (g > 1L) {
      names_1 <- paste(names_1, ".g", g, sep = "")
    }

    par_idx <- match(parlabel, names_1)
    drhodpsi_mat_1[[g]] <- delta_g[par_idx, , drop = FALSE]
  }
  drhodpsi_mat <- do.call(rbind, drhodpsi_mat_1)

  # tmp_prod <- ( t(drhodpsi_mat) %*% H_at_vartheta0 %*%
  #                 drhodpsi_mat %*% InvG_attheta0 %*%
  #                 H_attheta0 %*% InvG_attheta0 )
  tmp_prod <- (t(drhodpsi_mat) %*% h_at_vartheta0 %*%
    drhodpsi_mat %*% h0_inv %*% j0 %*% g0_inv)
  cov_tzztww <- 2 * sum(diag(tmp_prod))

  ##### Section 4: compute the adjusted PLRT and its p-value
  plrth0sat <- 2 * (h0_fx - sat_fx)
  plrth0sat_group <- 2 * (h0_fx_group - sat_fx_group)
  asym_mean_plrth0sat <- e_tzz - e_tww
  asym_var_plrth0sat <- var_tzz + var_tww - 2 * cov_tzztww
  scaling_factor <- (asym_mean_plrth0sat / (asym_var_plrth0sat / 2))
  fsa_plrt_sem <- (asym_mean_plrth0sat / (asym_var_plrth0sat / 2)) * plrth0sat
  adjusted_df <- (asym_mean_plrth0sat * asym_mean_plrth0sat) /
    (asym_var_plrth0sat / 2)
  # In some very few cases (simulations show very few cases
  #                         in small sample sizes)
  # the adjusted_df is a negative number, we should then
  # print a warning like: "The adjusted df is computed to be a negative number
  # and for this the first and second moment adjusted PLRT is not computed." .
  pvalue <- 1 - pchisq(fsa_plrt_sem, df = adjusted_df)

  list(
    PLRTH0Sat = plrth0sat, PLRTH0Sat.group = plrth0sat_group,
    stat = fsa_plrt_sem, df = adjusted_df, p.value = pvalue,
    scaling.factor = scaling_factor
  )
}
############################################################################
