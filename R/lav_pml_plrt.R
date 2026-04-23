lav_pml_plrt <- function(lavobject = NULL, lavmodel = NULL, lavdata = NULL,
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
      lavsamplestats = NULL,
      sample.cov = lav_model_sigma(lavmodel),
      sample.mean = lav_model_mu(lavmodel),
      sample.th = lav_model_th(lavmodel),
      sample.th.idx = lavsamplestats@th.idx
    )

  options2 <- options_1
  options2$optim.method <- "none"
  options2$optim.force.converged <- TRUE
  options2$check.start <- FALSE
  options2$check.gradient <- FALSE
  options2$check.post <- FALSE
  options2$check.vcov <- FALSE
  fitted_sat2 <- lavaan(model_sat2,
    slotOptions = options2, verbose = FALSE,
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
  pt_sat2 <- fitted_sat2@ParTable
  d_sat_idx <- pt_sat2$free[pt_sat2$free > 0L & pt_sat2$op != "|"]
                                                        # remove thresholds

  # Secondly, we need to specify the indices of the rows/columns of vcov(),
  # hessian, and variability matrix that refer to all SEM parameters
  # except thresholds.
  pt_1 <- lavpartable
  index_par <- pt_1$free[pt_1$free > 0L & pt_1$op != "|"]

  # Thirdly, specify the sample size.
  # nsize <- lavdata@nobs[[g]]
  nsize <- lavsamplestats@ntotal #TODO: check if var nsize needed ?? # nolint

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
  if (is.null(vcov_1)) {
    lavoptions$se <- "robust.huber.white"
    vcov_1 <- lav_model_vcov(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions,
      lavdata = lavdata,
      lavpartable = lavpartable,
      lavcache = lavcache
    )
  }
  inv_g_to_psipsi_attheta0 <- (lavsamplestats@ntotal * vcov_1)[index_par,
                             index_par, drop = FALSE] # G^psipsi(theta0)
  # below the lavaan function getHessian is used
  # Hattheta0 <- (-1) * H0.Hessian
  # Hattheta0 <- H0.Hessian
  # InvHattheta0 <- solve(Hattheta0)
  inv_hattheta0 <- attr(vcov_1, "E.inv")
  inv_h_to_psipsi_attheta0 <- inv_hattheta0[index_par, index_par, drop = FALSE]
                                                           # H^psipsi(theta0)
  if (lavmodel@eq.constraints) {
    in_1 <- inv_h_to_psipsi_attheta0
    in_npar <- ncol(in_1)

    # create `bordered' matrix
    if (nrow(lavmodel@con.jac) > 0L) {
      h <- lavmodel@con.jac[, index_par, drop = FALSE]
      inactive_idx <- attr(h, "inactive.idx")
      lambda <- lavmodel@con.lambda # lagrangean coefs
      if (length(inactive_idx) > 0L) {
        h <- h[-inactive_idx, , drop = FALSE]
        lambda <- lambda[-inactive_idx]
      }
      if (nrow(h) > 0L) {
        h0 <- matrix(0, nrow(h), nrow(h))
        h10 <- matrix(0, ncol(in_1), nrow(h))
        dl <- 2 * diag(lambda, nrow(h), nrow(h))
        # FIXME: better include inactive + slacks??
        e3 <- rbind( #TODO: check if var e3 needed ???  # nolint
          cbind(in_1, h10, t(h)),
          cbind(t(h10), dl, h0),
          cbind(h, h0, h0)
        )
        inv_of_inv_h_to_psipsi_attheta0 <-               # nolint
          MASS::ginv(in_1)[1:in_npar, 1:in_npar, drop = FALSE]
      } else {
        inv_of_inv_h_to_psipsi_attheta0 <- solve(in_1)   # nolint
      }
    }
  } else {
    # YR 26 June 2018: check for empty index.par (eg independence model)
    if (length(index_par) > 0L) {
      inv_of_inv_h_to_psipsi_attheta0 <-                   # nolint
        solve(inv_h_to_psipsi_attheta0) # [H^psipsi(theta0)]^(-1)
    } else {
      inv_of_inv_h_to_psipsi_attheta0 <- matrix(0, 0, 0)   # nolint
    }
  }

  h0tmp_prod1 <- inv_of_inv_h_to_psipsi_attheta0 %*% inv_g_to_psipsi_attheta0
  h0tmp_prod2 <- h0tmp_prod1 %*% h0tmp_prod1
  e_tww <- sum(diag(h0tmp_prod1)) # expected mean of first quadratic quantity
  var_tww <- 2 * sum(diag(h0tmp_prod2)) # variance of first quadratic quantity

  ##### Section 2: Compute the asymptotic mean and variance
  #####            of the second quadratic quantity.
  # Now we need to evaluate the fitted (polychoric) correlation/ covariance
  # matrix using the estimates of SEM parameters derived under the fitted model
  # which is the model of the null hypothesis. We also need to compute the
  # vcov matrix of these estimates (estimates of polychoric correlations)
  # as well as the related hessian and variability matrix.
  tmp_options <- fitted_sat2@Options
  tmp_options$se <- lavoptions$se
  vcov_sat2 <- lav_model_vcov(
    lavmodel = fitted_sat2@Model,
    lavsamplestats = fitted_sat2@SampleStats,
    lavoptions = tmp_options,
    lavdata = fitted_sat2@Data,
    lavpartable = fitted_sat2@ParTable,
    lavcache = fitted_sat2@Cache,
    use.ginv = TRUE
  )
  inv_g_to_sigmasigma_attheta0 <- lavsamplestats@ntotal * vcov_sat2[d_sat_idx,
                                d_sat_idx, drop = FALSE] # G^sigmasigma(theta0)
  # Hattheta0 <- (-1)* getHessian(fittedSat2)
  # Hattheta0 <- getHessian(fittedSat2)
  # InvHattheta0 <- solve(Hattheta0)
  inv_hattheta0 <- attr(vcov_sat2, "E.inv")
  inv_h_to_sigmasigma_attheta0 <-
              inv_hattheta0[d_sat_idx, d_sat_idx, drop = FALSE]
                                                    # H^sigmasigma(theta0)
  # Inv_of_InvH_to_sigmasigma_attheta0 <- solve(InvH_to_sigmasigma_attheta0)
  #                                                 #[H^sigmasigma(theta0)]^(-1)
  inv_of_inv_h_to_sigmasigma_attheta0 <- MASS::ginv(       # nolint
    inv_h_to_sigmasigma_attheta0,
    tol = .Machine$double.eps ^ (3 / 4)
  )
  h1tmp_prod1 <- inv_of_inv_h_to_sigmasigma_attheta0 %*%
                            inv_g_to_sigmasigma_attheta0
  h1tmp_prod2 <- h1tmp_prod1 %*% h1tmp_prod1
  e_tzz <- sum(diag(h1tmp_prod1)) # expected mean of the second
                                  # quadratic quantity
  var_tzz <- 2 * sum(diag(h1tmp_prod2)) # variance of the second
                                        # quadratic quantity

  ##### Section 3: Compute the asymptotic covariance of
  #####            the two quadratic quantities

  drhodpsi_mat_1 <- vector("list", length = lavsamplestats@ngroups)
  group_values <- lav_partable_group_values(fitted_sat2@ParTable)
  for (g in 1:lavsamplestats@ngroups) {
    # delta.g <- lav_model_delta(lavmodel)[[g]] # [[1]] to be substituted by g?
    # The above gives the derivatives of thresholds and polychoric correlations
    # with respect to SEM param (including thresholds) evaluated under H0.
    # From deltamat we need to exclude the rows and columns referring
    # to thresholds.
    # For this:

    # order of the rows: first the thresholds, then the correlations
    # we need to map the rows of delta.g to the rows/cols of H_at_vartheta0
    # of H1

    pt_1 <- fitted_sat2@ParTable
    pt_1$label <- lav_partable_labels(pt_1)
    free_idx <- which(pt_1$free > 0 & pt_1$op != "|" &
                       pt_1$group == group_values[g])
    parlabel <- pt_1$label[free_idx]

    # for now, we can assume that lav_model_delta will always return
    # the thresholds first, then the correlations
    #
    # later, we should add a (working) add.labels = TRUE option to
    # lav_model_delta
    # th.names <- lavobject@pta$vnames$th[[g]]
    # ov.names <- lavobject@pta$vnames$ov[[g]]
    # th.names <- lav_object_vnames(lavpartable, "th")
    # ov.names <- lav_object_vnames(lavpartable, "ov.nox")
    # ov.names.x <- lav_object_vnames(lavpartable, "ov.x")
    # tmp <- utils::combn(ov.names, 2)
    # cor.names <- paste(tmp[1,], "~~", tmp[2,], sep = "")

    # added by YR - 22 Okt 2017 #####################################
    # ov.names.x <- lav_object_vnames(lavpartable, "ov.x")
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
    delta_g <- lav_object_inspect_delta_internal(
      lavmodel = lavmodel,
      lavdata = lavdata, lavpartable = lavpartable,
      add.labels = TRUE, add.class = FALSE,
      drop.list.single.group = FALSE
    )[[g]]
    names_1 <- rownames(delta_g)
    if (g > 1L) {
      names_1 <- paste(names_1, ".g", g, sep = "")
    }

    par_idx <- match(parlabel, names_1)
    if (any(is.na(par_idx))) {
      lav_msg_warn(gettextf(
        "mismatch between DELTA labels and PAR labels!
        PARLABEL: %1$s, DELTA LABELS: %2$s", lav_msg_view(parlabel),
         lav_msg_view(names_1)))
    }

    drhodpsi_mat_1[[g]] <- delta_g[par_idx, index_par, drop = FALSE]
  }
  drhodpsi_mat <- do.call(rbind, drhodpsi_mat_1)

  tmp_prod <- t(drhodpsi_mat) %*% inv_of_inv_h_to_sigmasigma_attheta0 %*%
    drhodpsi_mat %*% inv_g_to_psipsi_attheta0 %*% h0tmp_prod1
  cov_tzztww <- 2 * sum(diag(tmp_prod))

  ##### Section 4: compute the adjusted PLRT and its p-value
  # PLRTH0Sat <- 2*nsize*(lavfit@fx - fittedSat@Fit@fx)
  plrth0sat <- 2 * (h0_fx - sat_fx)
  plrth0sat_group <- 2 * (h0_fx_group - sat_fx_group)
  asym_mean_plrth0sat <- e_tzz - e_tww
  # catch zero value for asym_mean_PLRTH0Sat
  if (asym_mean_plrth0sat == 0) {
    asym_var_plrth0sat <- 0
    scaling_factor <- as.numeric(NA)
    fsa_plrt_sem <- as.numeric(NA)
    adjusted_df <- as.integer(NA)
    pvalue <- as.numeric(NA)
  } else if (any(is.na(c(var_tzz, var_tww, cov_tzztww)))) {
    asym_var_plrth0sat <- as.numeric(NA)
    scaling_factor <- as.numeric(NA)
    fsa_plrt_sem <- as.numeric(NA)
    adjusted_df <- as.integer(NA)
    pvalue <- as.numeric(NA)
  } else {
    asym_var_plrth0sat <- var_tzz + var_tww - 2 * cov_tzztww
    scaling_factor <- (asym_mean_plrth0sat / (asym_var_plrth0sat / 2))
    fsa_plrt_sem <- (asym_mean_plrth0sat / (asym_var_plrth0sat / 2)) * plrth0sat
    adjusted_df <- (asym_mean_plrth0sat * asym_mean_plrth0sat) /
      (asym_var_plrth0sat / 2)
    # In some very few cases (simulations show very few cases in small
    # sample sizes) the adjusted_df is a negative number, we should then
    # print a warning like: "The adjusted df is computed to be a negative number
    # and for this the first and second moment adjusted PLRT is not computed."
    if (scaling_factor > 0) {
      pvalue <- 1 - pchisq(fsa_plrt_sem, df = adjusted_df)
    } else {
      pvalue <- as.numeric(NA)
    }
  }

  list(
    PLRTH0Sat = plrth0sat, PLRTH0Sat.group = plrth0sat_group,
    stat = fsa_plrt_sem, df = adjusted_df, p.value = pvalue,
    scaling.factor = scaling_factor
  )
}
############################################################################


lav_pml_object_aic_bic <- function(lavobject) {
  ########################## The code for PL version fo AIC and BIC
  # The following should be done because it is not the pl log-likelihood
  # that is maximized but a fit function that should be minimized. So, we
  # should find the value of log-PL at the estimated parameters through the
  # value of the fitted function.
  # The following may need to be updated if we change the fit function
  # so that it is correct for the case of missing values as well.

  log_pl <- lavobject@optim$logl
  nsize <- lavobject@SampleStats@ntotal

  # inverted observed unit information
  h_inv <- lavTech(lavobject, "inverted.information.observed")

  # first order unit information
  j <- lavTech(lavobject, "information.first.order")

  # trace (J %*% H.inv) = sum (J * t(H.inv))
  dim_theta <- sum(j * h_inv)


  # computations of PL versions of AIC and BIC
  pl_aic <- (-2) * log_pl + 2 * dim_theta
  pl_bic <- (-2) * log_pl + dim_theta * log(nsize)

  list(logPL = log_pl, PL_AIC = pl_aic, PL_BIC = pl_bic)
}
