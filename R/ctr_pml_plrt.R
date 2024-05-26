ctr_pml_plrt <- function(lavobject = NULL, lavmodel = NULL, lavdata = NULL,
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
      lavsamplestats = NULL,
      sample.cov = computeSigmaHat(lavmodel),
      sample.mean = computeMuHat(lavmodel),
      sample.th = computeTH(lavmodel),
      sample.th.idx = lavsamplestats@th.idx
    )

  Options2 <- Options
  Options2$optim.method <- "none"
  Options2$optim.force.converged <- TRUE
  Options2$check.start <- FALSE
  Options2$check.gradient <- FALSE
  Options2$check.post <- FALSE
  Options2$check.vcov <- FALSE
  fittedSat2 <- lavaan(ModelSat2,
    slotOptions = Options2,
    slotSampleStats = lavsamplestats,
    slotData = lavdata, slotCache = lavcache
  )

  # the code below was contributed by Myrsini Katsikatsou (Jan 2015)

  # for now, only a single group is supported:
  # g = 1L


  ########################### The code for PLRT for overall goodness of fit

  # First define the number of non-redundant elements of the (fitted)
  # covariance/correlation matrix of the underlying variables.
  # nvar <- lavmodel@nvar[[g]]
  # dSat <- nvar*(nvar-1)/2
  # if(length(lavmodel@num.idx[[g]]) > 0L) {
  #    dSat <- dSat + length(lavmodel@num.idx[[g]])
  # }

  # select `free' parameters (excluding thresholds) from fittedSat2 model
  PT.Sat2 <- fittedSat2@ParTable
  dSat.idx <- PT.Sat2$free[PT.Sat2$free > 0L & PT.Sat2$op != "|"]
                                                        # remove thresholds

  # Secondly, we need to specify the indices of the rows/columns of vcov(),
  # hessian, and variability matrix that refer to all SEM parameters
  # except thresholds.
  PT <- lavpartable
  index.par <- PT$free[PT$free > 0L & PT$op != "|"]

  # Thirdly, specify the sample size.
  # nsize <- lavdata@nobs[[g]]
  nsize <- lavsamplestats@ntotal

  # Now we can proceed to the computation of the quantities needed for PLRT.
  # Briefly, to say that PLRT is equal to the difference of two quadratic forms.
  # To compute the first and second moment adjusted PLRT we should compute
  # the asymptotic mean and variance of each quadratic quantity as well as
  # their asymptotic covariance.

  ##### Section 1. Compute the asymptotic mean and variance
  #####            of the first quadratic quantity
  # Below I assume that lavobject is the output of lavaan function. I guess
  # vcov(lavobject) can be substituted by VCOV object insed lavaan function
  # defined at lines 703 -708. But what is the object inside lavaan function
  # for getHessian(lavobject)?
  if (is.null(VCOV)) {
    lavoptions$se <- "robust.huber.white"
    VCOV <- lav_model_vcov(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions,
      lavdata = lavdata,
      lavpartable = lavpartable,
      lavcache = lavcache
    )
  }
  InvG_to_psipsi_attheta0 <- (lavsamplestats@ntotal * VCOV)[index.par,
                             index.par, drop = FALSE] # G^psipsi(theta0)
  # below the lavaan function getHessian is used
  # Hattheta0 <- (-1) * H0.Hessian
  # Hattheta0 <- H0.Hessian
  # InvHattheta0 <- solve(Hattheta0)
  InvHattheta0 <- attr(VCOV, "E.inv")
  InvH_to_psipsi_attheta0 <- InvHattheta0[index.par, index.par, drop = FALSE]
                                                           # H^psipsi(theta0)
  if (lavmodel@eq.constraints) {
    IN <- InvH_to_psipsi_attheta0
    IN.npar <- ncol(IN)

    # create `bordered' matrix
    if (nrow(lavmodel@con.jac) > 0L) {
      H <- lavmodel@con.jac[, index.par, drop = FALSE]
      inactive.idx <- attr(H, "inactive.idx")
      lambda <- lavmodel@con.lambda # lagrangean coefs
      if (length(inactive.idx) > 0L) {
        H <- H[-inactive.idx, , drop = FALSE]
        lambda <- lambda[-inactive.idx]
      }
      if (nrow(H) > 0L) {
        H0 <- matrix(0, nrow(H), nrow(H))
        H10 <- matrix(0, ncol(IN), nrow(H))
        DL <- 2 * diag(lambda, nrow(H), nrow(H))
        # FIXME: better include inactive + slacks??
        E3 <- rbind(
          cbind(IN, H10, t(H)),
          cbind(t(H10), DL, H0),
          cbind(H, H0, H0)
        )
        Inv_of_InvH_to_psipsi_attheta0 <-
          MASS::ginv(IN)[1:IN.npar, 1:IN.npar, drop = FALSE]
      } else {
        Inv_of_InvH_to_psipsi_attheta0 <- solve(IN)
      }
    }
  } else {
    # YR 26 June 2018: check for empty index.par (eg independence model)
    if (length(index.par) > 0L) {
      Inv_of_InvH_to_psipsi_attheta0 <-
        solve(InvH_to_psipsi_attheta0) # [H^psipsi(theta0)]^(-1)
    } else {
      Inv_of_InvH_to_psipsi_attheta0 <- matrix(0, 0, 0)
    }
  }

  H0tmp_prod1 <- Inv_of_InvH_to_psipsi_attheta0 %*% InvG_to_psipsi_attheta0
  H0tmp_prod2 <- H0tmp_prod1 %*% H0tmp_prod1
  E_tww <- sum(diag(H0tmp_prod1)) # expected mean of first quadratic quantity
  var_tww <- 2 * sum(diag(H0tmp_prod2)) # variance of first quadratic quantity

  ##### Section 2: Compute the asymptotic mean and variance
  #####            of the second quadratic quantity.
  # Now we need to evaluate the fitted (polychoric) correlation/ covariance
  # matrix using the estimates of SEM parameters derived under the fitted model
  # which is the model of the null hypothesis. We also need to compute the
  # vcov matrix of these estimates (estimates of polychoric correlations)
  # as well as the related hessian and variability matrix.
  tmp.options <- fittedSat2@Options
  tmp.options$se <- lavoptions$se
  VCOV.Sat2 <- lav_model_vcov(
    lavmodel = fittedSat2@Model,
    lavsamplestats = fittedSat2@SampleStats,
    lavoptions = tmp.options,
    lavdata = fittedSat2@Data,
    lavpartable = fittedSat2@ParTable,
    lavcache = fittedSat2@Cache,
    use.ginv = TRUE
  )
  InvG_to_sigmasigma_attheta0 <- lavsamplestats@ntotal * VCOV.Sat2[dSat.idx,
                                dSat.idx, drop = FALSE] # G^sigmasigma(theta0)
  # Hattheta0 <- (-1)* getHessian(fittedSat2)
  # Hattheta0 <- getHessian(fittedSat2)
  # InvHattheta0 <- solve(Hattheta0)
  InvHattheta0 <- attr(VCOV.Sat2, "E.inv")
  InvH_to_sigmasigma_attheta0 <- InvHattheta0[dSat.idx, dSat.idx, drop = FALSE]
                                                          # H^sigmasigma(theta0)
  # Inv_of_InvH_to_sigmasigma_attheta0 <- solve(InvH_to_sigmasigma_attheta0)
  #                                                 #[H^sigmasigma(theta0)]^(-1)
  Inv_of_InvH_to_sigmasigma_attheta0 <- MASS::ginv(InvH_to_sigmasigma_attheta0,
    tol = .Machine$double.eps^(3 / 4)
  )
  H1tmp_prod1 <- Inv_of_InvH_to_sigmasigma_attheta0 %*%
                            InvG_to_sigmasigma_attheta0
  H1tmp_prod2 <- H1tmp_prod1 %*% H1tmp_prod1
  E_tzz <- sum(diag(H1tmp_prod1)) # expected mean of the second
                                  # quadratic quantity
  var_tzz <- 2 * sum(diag(H1tmp_prod2)) # variance of the second
                                        # quadratic quantity

  ##### Section 3: Compute the asymptotic covariance of
  #####            the two quadratic quantities

  drhodpsi_MAT <- vector("list", length = lavsamplestats@ngroups)
  group.values <- lav_partable_group_values(fittedSat2@ParTable)
  for (g in 1:lavsamplestats@ngroups) {
    # delta.g <- computeDelta(lavmodel)[[g]] # [[1]] to be substituted by g?
    # The above gives the derivatives of thresholds and polychoric correlations
    # with respect to SEM param (including thresholds) evaluated under H0.
    # From deltamat we need to exclude the rows and columns referring
    # to thresholds.
    # For this:

    # order of the rows: first the thresholds, then the correlations
    # we need to map the rows of delta.g to the rows/cols of H_at_vartheta0
    # of H1

    PT <- fittedSat2@ParTable
    PT$label <- lav_partable_labels(PT)
    free.idx <- which(PT$free > 0 & PT$op != "|" & PT$group == group.values[g])
    PARLABEL <- PT$label[free.idx]

    # for now, we can assume that computeDelta will always return
    # the thresholds first, then the correlations
    #
    # later, we should add a (working) add.labels = TRUE option to
    # computeDelta
    # th.names <- lavobject@pta$vnames$th[[g]]
    # ov.names <- lavobject@pta$vnames$ov[[g]]
    # th.names <- lavNames(lavpartable, "th")
    # ov.names <- lavNames(lavpartable, "ov.nox")
    # ov.names.x <- lavNames(lavpartable, "ov.x")
    # tmp <- utils::combn(ov.names, 2)
    # cor.names <- paste(tmp[1,], "~~", tmp[2,], sep = "")

    # added by YR - 22 Okt 2017 #####################################
    # ov.names.x <- lavNames(lavpartable, "ov.x")
    # if(length(ov.names.x)) {
    #    slope.names <- apply(expand.grid(ov.names, ov.names.x), 1L,
    #                             paste, collapse = "~")
    # } else {
    #    slope.names <- character(0L)
    # }
    #################################################################

    # NAMES <- c(th.names, slope.names, cor.names)

    # added by YR - 26 April 2018, for 0.6-1
    # we now can get 'labelled' delta rownames
    delta.g <- lav_object_inspect_delta_internal(
      lavmodel = lavmodel,
      lavdata = lavdata, lavpartable = lavpartable,
      add.labels = TRUE, add.class = FALSE,
      drop.list.single.group = FALSE
    )[[g]]
    NAMES <- rownames(delta.g)
    if (g > 1L) {
      NAMES <- paste(NAMES, ".g", g, sep = "")
    }

    par.idx <- match(PARLABEL, NAMES)
    if (any(is.na(par.idx))) {
      lav_msg_warn(gettextf(
        "mismatch between DELTA labels and PAR labels!
        PARLABEL: %1$s, DELTA LABELS: %2$s", lav_msg_view(PARLABEL),
         lav_msg_view(NAMES)))
    }

    drhodpsi_MAT[[g]] <- delta.g[par.idx, index.par, drop = FALSE]
  }
  drhodpsi_mat <- do.call(rbind, drhodpsi_MAT)

  tmp_prod <- t(drhodpsi_mat) %*% Inv_of_InvH_to_sigmasigma_attheta0 %*%
    drhodpsi_mat %*% InvG_to_psipsi_attheta0 %*% H0tmp_prod1
  cov_tzztww <- 2 * sum(diag(tmp_prod))

  ##### Section 4: compute the adjusted PLRT and its p-value
  # PLRTH0Sat <- 2*nsize*(lavfit@fx - fittedSat@Fit@fx)
  PLRTH0Sat <- 2 * (H0.fx - SAT.fx)
  PLRTH0Sat.group <- 2 * (H0.fx.group - SAT.fx.group)
  asym_mean_PLRTH0Sat <- E_tzz - E_tww
  # catch zero value for asym_mean_PLRTH0Sat
  if (asym_mean_PLRTH0Sat == 0) {
    asym_var_PLRTH0Sat <- 0
    scaling.factor <- as.numeric(NA)
    FSA_PLRT_SEM <- as.numeric(NA)
    adjusted_df <- as.integer(NA)
    pvalue <- as.numeric(NA)
  } else if (any(is.na(c(var_tzz, var_tww, cov_tzztww)))) {
    asym_var_PLRTH0Sat <- as.numeric(NA)
    scaling.factor <- as.numeric(NA)
    FSA_PLRT_SEM <- as.numeric(NA)
    adjusted_df <- as.integer(NA)
    pvalue <- as.numeric(NA)
  } else {
    asym_var_PLRTH0Sat <- var_tzz + var_tww - 2 * cov_tzztww
    scaling.factor <- (asym_mean_PLRTH0Sat / (asym_var_PLRTH0Sat / 2))
    FSA_PLRT_SEM <- (asym_mean_PLRTH0Sat / (asym_var_PLRTH0Sat / 2)) * PLRTH0Sat
    adjusted_df <- (asym_mean_PLRTH0Sat * asym_mean_PLRTH0Sat) /
      (asym_var_PLRTH0Sat / 2)
    # In some very few cases (simulations show very few cases in small
    # sample sizes) the adjusted_df is a negative number, we should then
    # print a warning like: "The adjusted df is computed to be a negative number
    # and for this the first and second moment adjusted PLRT is not computed."
    if (scaling.factor > 0) {
      pvalue <- 1 - pchisq(FSA_PLRT_SEM, df = adjusted_df)
    } else {
      pvalue <- as.numeric(NA)
    }
  }

  list(
    PLRTH0Sat = PLRTH0Sat, PLRTH0Sat.group = PLRTH0Sat.group,
    stat = FSA_PLRT_SEM, df = adjusted_df, p.value = pvalue,
    scaling.factor = scaling.factor
  )
}
############################################################################


ctr_pml_aic_bic <- function(lavobject) {
  ########################## The code for PL version fo AIC and BIC
  # The following should be done because it is not the pl log-likelihood
  # that is maximized but a fit function that should be minimized. So, we
  # should find the value of log-PL at the estimated parameters through the
  # value of the fitted function.
  # The following may need to be updated if we change the fit function
  # so that it is correct for the case of missing values as well.

  logPL <- lavobject@optim$logl
  nsize <- lavobject@SampleStats@ntotal

  # inverted observed unit information
  H.inv <- lavTech(lavobject, "inverted.information.observed")

  # first order unit information
  J <- lavTech(lavobject, "information.first.order")

  # trace (J %*% H.inv) = sum (J * t(H.inv))
  dimTheta <- sum(J * H.inv)


  # computations of PL versions of AIC and BIC
  PL_AIC <- (-2) * logPL + 2 * dimTheta
  PL_BIC <- (-2) * logPL + dimTheta * log(nsize)

  list(logPL = logPL, PL_AIC = PL_AIC, PL_BIC = PL_BIC)
}
