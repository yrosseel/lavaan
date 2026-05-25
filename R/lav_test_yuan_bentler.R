# - 0.6-13: fix multiple-group UG^2 bug (reported by Gronneberg, Foldnes and
#           Moss) when Satterthwaite = TRUE, ngroups > 1, and eq constraints.
#
#           Note however that Satterthwaite = FALSE always (for now), so
#           the fix has no (visible) effect

lav_test_yb <- function(lavobject = NULL,
                                  lavsamplestats = NULL,
                                  lavmodel = NULL,
                                  lavimplied = NULL,
                                  lavh1 = NULL,
                                  lavoptions = NULL,
                                  lavdata = NULL,
                                  test_unscaled = NULL,
                                  e_inv = NULL,
                                  b0_group = NULL,
                                  test = "yuan.bentler",
                                  mimic = "lavaan",
                                  # method          = "default",
                                  ug2_old_approach = FALSE,
                                  return_ugamma = FALSE) {
  test_1 <- list()

  if (!is.null(lavobject)) {
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavoptions <- lavobject@Options
    # lavpartable <- lavobject@ParTable
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavdata <- lavobject@Data
    test_1$standard <- lavobject@test[[1]]
  } else {
    test_1$standard <- test_unscaled
  }

  # ug2.old.approach
  if (missing(ug2_old_approach)) {
    if (!is.null(lavoptions$ug2.old.approach)) {
      ug2_old_approach <- lavoptions$ug2.old.approach
    } else {
      ug2_old_approach <- FALSE
    }
  }

  # E.inv ok?
  if (length(lavoptions$information) == 1L &&
    length(lavoptions$h1.information) == 1L &&
    length(lavoptions$observed.information) == 1L) {
    e_inv_recompute <- FALSE
  } else if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    e_inv_recompute <- FALSE
  } else {
    e_inv_recompute <- TRUE
    # change information options
    lavoptions$information[1] <- lavoptions$information[2]
    lavoptions$h1.information[1] <- lavoptions$h1.information[2]
    lavoptions$observed.information[1] <- lavoptions$observed.information[2]
  }
  if (!is.null(e_inv)) {
    e_inv_recompute <- FALSE # user-provided
  }

  # check test
  if (!all(test %in% c(
    "yuan.bentler",
    "yuan.bentler.mplus"
  ))) {
    lav_msg_warn(gettext("test must be one of `yuan.bentler', or
                         `yuan.bentler.mplus'; will use `yuan.bentler' only"))
    test <- "yuan.bentler"
  }

  # information
  # information <- lavoptions$information[1]

  # ndat
  ndat <- numeric(lavsamplestats@ngroups)


  # do we have E.inv?
  if (is.null(e_inv) || e_inv_recompute) {
    e_inv <- try(
      lav_model_info(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavoptions = lavoptions,
        extra = FALSE,
        augmented = TRUE,
        inverted = TRUE
      ),
      silent = TRUE
    )
    if (inherits(e_inv, "try-error")) {
      if (return_ugamma) {
        lav_msg_warn(gettext(
          "could not invert information matrix needed for UGamma"))
        return(NULL)
      } else {
        test_1$standard$stat <- as.numeric(NA)
        test_1$standard$stat.group <- rep(as.numeric(NA), lavdata@ngroups)
        test_1$standard$pvalue <- as.numeric(NA)
        test_1[[test[1]]] <- c(test_1$standard,
          scaling.factor = as.numeric(NA),
          shift.parameter = as.numeric(NA),
          label = character(0)
        )
        lav_msg_warn(gettext("could not invert information matrix needed for
                             robust test statistic"))
        test_1[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when
                       # robust test is detected for some but not all models
        return(test_1)
      }
    }
  }

  # catch df == 0
  if (test_1$standard$df == 0L || test_1$standard$df < 0) {
    test_1[[test[1]]] <- c(test_1$standard,
      scaling.factor = as.numeric(NA),
      label = character(0)
    )
    test_1[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when
                   # robust test is detected for some but not all models
    return(test_1)
  }

  # mean and variance adjusted?
  # satterthwaite <- FALSE # for now
  # if(any(test %in% c("mean.var.adjusted", "scaled.shifted"))) {
  #    Satterthwaite <- TRUE
  # }

  # FIXME: should we not always use 'unstructured' here?
  # if the model is, say, the independence model, the
  # 'structured' information (A1) will be so far away from B1
  # that we will end up with 'NA'
  h1_options <- lavoptions
  if (test == "yuan.bentler.mplus") {
    # always 'unstructured' H1 information
    h1_options$h1.information <- "unstructured"
  }

  # A1 is usually expected or observed
  a1_group <- lav_model_h1_info(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    lavoptions = h1_options
  )
  # B1 is always first.order
  b1_group <- lav_model_h1_info_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    lavoptions = h1_options
  )

  if (test == "yuan.bentler.mplus") {
    if (is.null(b0_group)) {
      b0 <- lav_model_info_firstorder(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        extra = TRUE,
        check_pd = FALSE,
        augmented = FALSE,
        inverted = FALSE
      )
      b0_group <- attr(b0, "B0.group")
    }
    trace_ugamma <-
      lav_test_yb_mplus_trace(
        lavsamplestats = lavsamplestats,
        a1_group = a1_group,
        b1_group = b1_group,
        b0_group = b0_group,
        e_inv = e_inv,
        meanstructure = lavmodel@meanstructure
      )
  } else if (test == "yuan.bentler") {
    # compute Delta
    delta <- lav_model_delta(lavmodel = lavmodel)

    # compute Omega/Gamma
    omega <- lav_model_h1_omega(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions
    )

    # compute trace 'U %*% Gamma' (or 'U %*% Omega')
    trace_ugamma <- lav_test_yb_trace(
      lavsamplestats   = lavsamplestats,
      meanstructure    = lavmodel@meanstructure,
      a1_group         = a1_group,
      b1_group         = b1_group,
      delta            = delta,
      omega            = omega,
      e_inv            = e_inv,
      ug2_old_approach = ug2_old_approach,
      satterthwaite    = FALSE
    ) # for now
  }

  # unscaled test
  df <- test_1$standard$df

  scaling_factor <- trace_ugamma / df
  if (scaling_factor < 0) scaling_factor <- as.numeric(NA)
  chisq_scaled <- test_1$standard$stat / scaling_factor
  pvalue_scaled <- 1 - pchisq(chisq_scaled, df)

  ndat <- sum(attr(trace_ugamma, "h1.ndat"))
  npar <- lavmodel@nx.free

  scaling_factor_h1 <- sum(attr(trace_ugamma, "h1")) / ndat
  scaling_factor_h0 <- sum(attr(trace_ugamma, "h0")) / npar
  trace_ugamma2 <- attr(trace_ugamma, "trace.UGamma2")
  attributes(trace_ugamma) <- NULL

  if ("yuan.bentler" %in% test) {
    test_1$yuan.bentler <-
      list(
        test = test,
        stat = chisq_scaled,
        stat.group = (test_1$standard$stat.group /
          scaling_factor),
        df = df,
        pvalue = pvalue_scaled,
        scaling.factor = scaling_factor,
        scaling.factor.h1 = scaling_factor_h1,
        scaling.factor.h0 = scaling_factor_h0,
        label = "Yuan-Bentler correction",
        trace.UGamma = trace_ugamma,
        trace.UGamma2 = trace_ugamma2,
        scaled.test.stat = test_1$standard$stat,
        scaled.test = test_1$standard$test
      )
  } else if ("yuan.bentler.mplus" %in% test) {
    test_1$yuan.bentler.mplus <-
      list(
        test = test,
        stat = chisq_scaled,
        stat.group = (test_1$standard$stat.group /
          scaling_factor),
        df = df,
        pvalue = pvalue_scaled,
        scaling.factor = scaling_factor,
        scaling.factor.h1 = scaling_factor_h1,
        scaling.factor.h0 = scaling_factor_h0,
        label =
          "Yuan-Bentler correction (Mplus variant)",
        trace.UGamma = trace_ugamma,
        trace.UGamma2 = as.numeric(NA),
        scaled.test.stat = test_1$standard$stat,
        scaled.test = test_1$standard$test
      )
  }

  test_1
}


lav_test_yb_trace <- function(lavsamplestats = lavsamplestats,
                                        meanstructure = TRUE,
                                        a1_group = NULL,
                                        b1_group = NULL,
                                        delta = NULL,
                                        omega = NULL,
                                        e_inv = NULL,
                                        ug2_old_approach = FALSE,
                                        satterthwaite = FALSE) {
  # we always assume a meanstructure (nope, not any longer, since 0.6)
  # meanstructure <- TRUE

  ngroups <- lavsamplestats@ngroups

  trace_h1 <- attr(omega, "trace.h1")
  h1_ndat <- attr(omega, "h1.ndat")

  if (ug2_old_approach || !satterthwaite) {
    trace_ugamma <- numeric(ngroups)
    trace_ugamma2 <- numeric(ngroups)
    trace_h0 <- numeric(ngroups)

    for (g in 1:ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal

      a1 <- a1_group[[g]] * fg
      b1 <- b1_group[[g]] * fg
      mm_delta <- delta[[g]]
      gamma_g <- omega[[g]] / fg

      d_einv_t_d <- mm_delta %*% tcrossprod(e_inv, mm_delta)

      # trace.h1[g] <- sum( B1 * t( A1.inv ) )
      # fg cancels out: trace.h1[g] <- sum( fg*B1 * t( 1/fg*A1.inv ) )
      trace_h0[g] <- sum(b1 * d_einv_t_d)
      # trace.UGamma[g] <- trace.h1[g] - trace.h0[g]
      u <- a1 - a1 %*% d_einv_t_d %*% a1
      trace_ugamma[g] <- sum(u * gamma_g)

      if (satterthwaite) {
        ug <- u %*% gamma_g
        trace_ugamma2[g] <- sum(ug * t(ug))
      }
    } # g
    trace_ugamma <- sum(trace_ugamma)
    attr(trace_ugamma, "h1") <- trace_h1
    attr(trace_ugamma, "h0") <- trace_h0
    attr(trace_ugamma, "h1.ndat") <- h1_ndat
    if (satterthwaite) {
      attr(trace_ugamma, "trace.UGamma2") <- sum(trace_ugamma2)
    }
  } else {
    trace_ugamma <- trace_ugamma2 <- ug <- as.numeric(NA)
    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    # if(Satterthwaite) {

    a1_f <- a1_group
    for (g in 1:ngroups) {
      a1_f[[g]] <- a1_group[[g]] * fg[g]
    }
    a1_all <- lav_mat_bdiag(a1_f)

    b1_f <- b1_group
    for (g in 1:ngroups) {
      b1_f[[g]] <- b1_group[[g]] * fg[g]
    }
    b1_all <- lav_mat_bdiag(b1_f)

    gamma_f <- omega
    for (g in 1:ngroups) {
      gamma_f[[g]] <- 1 / fg[g] * omega[[g]]
    }
    gamma_all <- lav_mat_bdiag(gamma_f)
    delta_all <- do.call("rbind", delta)

    d_einv_t_d <- delta_all %*% tcrossprod(e_inv, delta_all)

    trace_h0 <- sum(b1_all * d_einv_t_d)
    u_all <- a1_all - a1_all %*% d_einv_t_d %*% a1_all
    trace_ugamma <- sum(u_all * gamma_all)

    attr(trace_ugamma, "h1") <- sum(trace_h1)
    attr(trace_ugamma, "h0") <- trace_h0
    attr(trace_ugamma, "h1.ndat") <- sum(h1_ndat)
    if (satterthwaite) {
      ug <- u_all %*% gamma_all
      trace_ugamma2 <- sum(ug * t(ug))
      attr(trace_ugamma, "trace.UGamma2") <- trace_ugamma2
    }

    # } else {
  }

  trace_ugamma
}

lav_test_yb_mplus_trace <- function(lavsamplestats = NULL,
                                              a1_group = NULL,
                                              b1_group = NULL,
                                              b0_group = NULL,
                                              e_inv = NULL,
                                              meanstructure = TRUE) {
  # typical for Mplus:
  # - do NOT use the YB formula, but use an approximation
  #   relying  on A0 ~= Delta' A1 Delta and the same for B0
  #
  # NOTE: if A0 is based on the hessian, then A0 only approximates
  #       Delta' A1 Delta
  #
  # - always use h1.information = "unstructured"!!!

  # ngroups <- lavsamplestats@ngroups

  trace_ugamma <- numeric(lavsamplestats@ngroups)
  trace_h1 <- numeric(lavsamplestats@ngroups)
  trace_h0 <- numeric(lavsamplestats@ngroups)
  h1_ndat <- numeric(lavsamplestats@ngroups)

  for (g in 1:lavsamplestats@ngroups) {
    # group weight
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal

    a1_1 <- a1_group[[g]]
    b1 <- b1_group[[g]]

    # mask independent 'fixed-x' variables
    zero_idx <- which(diag(a1_1) == 0)
    if (length(zero_idx) > 0L) {
      a1_inv_1 <- matrix(0, nrow(a1_1), ncol(a1_1))
      a1 <- a1_1[-zero_idx, -zero_idx]
      a1_inv <- solve(a1)
      a1_inv_1[-zero_idx, -zero_idx] <- a1_inv
    } else {
      a1_inv_1 <- solve(a1_1)
    }
    h1_ndat[g] <- ncol(a1_1) - length(zero_idx)

    # if data is complete, why not just A1 %*% Gamma?
    trace_h1[g] <- sum(b1 * t(a1_inv_1))
    trace_h0[g] <- fg * sum(b0_group[[g]] * t(e_inv))
    trace_ugamma[g] <- (trace_h1[g] - trace_h0[g])
  }

  # we take the sum here
  trace_ugamma <- sum(trace_ugamma)

  attr(trace_ugamma, "h1") <- trace_h1
  attr(trace_ugamma, "h0") <- trace_h0
  attr(trace_ugamma, "h1.ndat") <- h1_ndat

  trace_ugamma
}
