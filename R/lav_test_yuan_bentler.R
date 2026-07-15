# - 0.6-13: fix multiple-group UG^2 bug (reported by Gronneberg, Foldnes and
#           Moss) when Satterthwaite = TRUE, ngroups > 1, and eq constraints.
#
#           Note however that Satterthwaite = FALSE always (for now), so
#           the fix has no (visible) effect
# - 0.7-1:  the traces are computed by the streamed worker
#           lav_test_ug_trace_stream() (lav_test_utils.R): no more
#           pstar x pstar products (A1 %*% ... %*% A1)

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
  lavoptions <- lav_test_change_options(lavoptions)
  e_inv_recompute <- attr(lavoptions, "recompute")
  attr(lavoptions, "recompute") <- NULL
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
        na_out <- lav_test_scaled_na(test_1$standard,
          test1 = test[1],
          ngroups = lavdata@ngroups
        )
        test_1$standard <- na_out$standard
        test_1[[test[1]]] <- na_out$failed
        lav_msg_warn(gettext("could not invert information matrix needed for
                             robust test statistic"))
        return(test_1)
      }
    }
  }

  # catch df == 0
  if (test_1$standard$df == 0L || test_1$standard$df < 0) {
    na_out <- lav_test_scaled_na(test_1$standard,
      test1 = test[1],
      na_standard = FALSE, shift = FALSE
    )
    test_1[[test[1]]] <- na_out$failed
    return(test_1)
  }

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

# trace of U %*% Gamma (U %*% Omega), with U = A1 - A1 Delta E.inv Delta' A1,
# via the streamed worker (same structure as the Satorra-Bentler U, with
# V := A1); trace_h0 = tr(B1 Delta E.inv Delta') is accumulated per group
# as tr(E.inv Delta' B1 Delta)
lav_test_yb_trace <- function(lavsamplestats = lavsamplestats,
                                        meanstructure = TRUE,
                                        a1_group = NULL,
                                        b1_group = NULL,
                                        delta = NULL,
                                        omega = NULL,
                                        e_inv = NULL,
                                        ug2_old_approach = FALSE,
                                        satterthwaite = FALSE) {
  ngroups <- lavsamplestats@ngroups
  fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal

  trace_h1 <- attr(omega, "trace.h1")
  h1_ndat <- attr(omega, "h1.ndat")

  # tr(U Gamma) (and tr((U Gamma)^2) if satterthwaite)
  out <- lav_test_ug_trace_stream(
    m_gamma = omega, delta = delta, wls_v = a1_group, e_inv = e_inv,
    fg = fg, satterthwaite = satterthwaite,
    old_approach = ug2_old_approach
  )

  # trace.h0 per group; fg cancels against the Gamma_g / fg convention:
  # sum( (fg * B1_g) * (Delta E.inv Delta') ) = fg * tr(E.inv Delta' B1 Delta)
  trace_h0 <- numeric(ngroups)
  for (g in seq_len(ngroups)) {
    trace_h0[g] <- fg[g] *
      sum(e_inv * crossprod(delta[[g]], b1_group[[g]] %*% delta[[g]]))
  }

  trace_ugamma <- out$trace.UGamma
  attr(trace_ugamma, "h1") <- trace_h1
  attr(trace_ugamma, "h0") <- trace_h0
  attr(trace_ugamma, "h1.ndat") <- h1_ndat
  if (satterthwaite) {
    attr(trace_ugamma, "trace.UGamma2") <- out$trace.UGamma2
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
