# - 0.6-13: fix multiple-group UG^2 bug (reported by Gronneberg, Foldnes and
#           Moss) when Satterthwaite = TRUE, ngroups > 1, and eq constraints.
#
#           Note however that Satterthwaite = FALSE always (for now), so
#           the fix has no (visible) effect

lav_test_yuan_bentler <- function(lavobject = NULL,
                                  lavsamplestats = NULL,
                                  lavmodel = NULL,
                                  lavimplied = NULL,
                                  lavh1 = NULL,
                                  lavoptions = NULL,
                                  lavdata = NULL,
                                  TEST.unscaled = NULL,
                                  E.inv = NULL,
                                  B0.group = NULL,
                                  test = "yuan.bentler",
                                  mimic = "lavaan",
                                  # method          = "default",
                                  ug2.old.approach = FALSE,
                                  return.ugamma = FALSE) {
  TEST <- list()

  if (!is.null(lavobject)) {
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavoptions <- lavobject@Options
    lavpartable <- lavobject@ParTable
    lavimplied <- lavobject@implied
    lavh1 <- lavobject@h1
    lavdata <- lavobject@Data
    TEST$standard <- lavobject@test[[1]]
  } else {
    TEST$standard <- TEST.unscaled
  }

  # ug2.old.approach
  if (missing(ug2.old.approach)) {
    if (!is.null(lavoptions$ug2.old.approach)) {
      ug2.old.approach <- lavoptions$ug2.old.approach
    } else {
      ug2.old.approach <- FALSE
    }
  }

  # E.inv ok?
  if (length(lavoptions$information) == 1L &&
    length(lavoptions$h1.information) == 1L &&
    length(lavoptions$observed.information) == 1L) {
    E.inv.recompute <- FALSE
  } else if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    E.inv.recompute <- FALSE
  } else {
    E.inv.recompute <- TRUE
    # change information options
    lavoptions$information[1] <- lavoptions$information[2]
    lavoptions$h1.information[1] <- lavoptions$h1.information[2]
    lavoptions$observed.information[1] <- lavoptions$observed.information[2]
  }
  if (!is.null(E.inv)) {
    E.inv.recompute <- FALSE # user-provided
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
  information <- lavoptions$information[1]

  # ndat
  ndat <- numeric(lavsamplestats@ngroups)


  # do we have E.inv?
  if (is.null(E.inv) || E.inv.recompute) {
    E.inv <- try(
      lav_model_information(
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
    if (inherits(E.inv, "try-error")) {
      if (return.ugamma) {
        lav_msg_warn(gettext(
          "could not invert information matrix needed for UGamma"))
        return(NULL)
      } else {
        TEST$standard$stat <- as.numeric(NA)
        TEST$standard$stat.group <- rep(as.numeric(NA), lavdata@ngroups)
        TEST$standard$pvalue <- as.numeric(NA)
        TEST[[test[1]]] <- c(TEST$standard,
          scaling.factor = as.numeric(NA),
          shift.parameter = as.numeric(NA),
          label = character(0)
        )
        lav_msg_warn(gettext("could not invert information [matrix needed for
                             robust test statistic"))
        TEST[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when robust test is detected for some but not all models
        return(TEST)
      }
    }
  }

  # catch df == 0
  if (TEST$standard$df == 0L || TEST$standard$df < 0) {
    TEST[[test[1]]] <- c(TEST$standard,
      scaling.factor = as.numeric(NA),
      label = character(0)
    )
    TEST[[test[1]]]$test <- test[1] # to prevent lavTestLRT error when robust test is detected for some but not all models
    return(TEST)
  }

  # mean and variance adjusted?
  Satterthwaite <- FALSE # for now
  # if(any(test %in% c("mean.var.adjusted", "scaled.shifted"))) {
  #    Satterthwaite <- TRUE
  # }

  # FIXME: should we not always use 'unstructured' here?
  # if the model is, say, the independence model, the
  # 'structured' information (A1) will be so far away from B1
  # that we will end up with 'NA'
  h1.options <- lavoptions
  if (test == "yuan.bentler.mplus") {
    # always 'unstructured' H1 information
    h1.options$h1.information <- "unstructured"
  }

  # A1 is usually expected or observed
  A1.group <- lav_model_h1_information(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    lavoptions = h1.options
  )
  # B1 is always first.order
  B1.group <- lav_model_h1_information_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    lavoptions = h1.options
  )

  if (test == "yuan.bentler.mplus") {
    if (is.null(B0.group)) {
      B0 <- lav_model_information_firstorder(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        extra = TRUE,
        check.pd = FALSE,
        augmented = FALSE,
        inverted = FALSE
      )
      B0.group <- attr(B0, "B0.group")
    }
    trace.UGamma <-
      lav_test_yuan_bentler_mplus_trace(
        lavsamplestats = lavsamplestats,
        A1.group = A1.group,
        B1.group = B1.group,
        B0.group = B0.group,
        E.inv = E.inv,
        meanstructure = lavmodel@meanstructure
      )
  } else if (test == "yuan.bentler") {
    # compute Delta
    Delta <- computeDelta(lavmodel = lavmodel)

    # compute Omega/Gamma
    Omega <- lav_model_h1_omega(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions
    )

    # compute trace 'U %*% Gamma' (or 'U %*% Omega')
    trace.UGamma <- lav_test_yuan_bentler_trace(
      lavsamplestats   = lavsamplestats,
      meanstructure    = lavmodel@meanstructure,
      A1.group         = A1.group,
      B1.group         = B1.group,
      Delta            = Delta,
      Omega            = Omega,
      E.inv            = E.inv,
      ug2.old.approach = ug2.old.approach,
      Satterthwaite    = FALSE
    ) # for now
  }

  # unscaled test
  df <- TEST$standard$df

  scaling.factor <- trace.UGamma / df
  if (scaling.factor < 0) scaling.factor <- as.numeric(NA)
  chisq.scaled <- TEST$standard$stat / scaling.factor
  pvalue.scaled <- 1 - pchisq(chisq.scaled, df)

  ndat <- sum(attr(trace.UGamma, "h1.ndat"))
  npar <- lavmodel@nx.free

  scaling.factor.h1 <- sum(attr(trace.UGamma, "h1")) / ndat
  scaling.factor.h0 <- sum(attr(trace.UGamma, "h0")) / npar
  trace.UGamma2 <- attr(trace.UGamma, "trace.UGamma2")
  attributes(trace.UGamma) <- NULL

  if ("yuan.bentler" %in% test) {
    TEST$yuan.bentler <-
      list(
        test = test,
        stat = chisq.scaled,
        stat.group = (TEST$standard$stat.group /
          scaling.factor),
        df = df,
        pvalue = pvalue.scaled,
        scaling.factor = scaling.factor,
        scaling.factor.h1 = scaling.factor.h1,
        scaling.factor.h0 = scaling.factor.h0,
        label = "Yuan-Bentler correction",
        trace.UGamma = trace.UGamma,
        trace.UGamma2 = trace.UGamma2,
        scaled.test.stat = TEST$standard$stat,
        scaled.test = TEST$standard$test
      )
  } else if ("yuan.bentler.mplus" %in% test) {
    TEST$yuan.bentler.mplus <-
      list(
        test = test,
        stat = chisq.scaled,
        stat.group = (TEST$standard$stat.group /
          scaling.factor),
        df = df,
        pvalue = pvalue.scaled,
        scaling.factor = scaling.factor,
        scaling.factor.h1 = scaling.factor.h1,
        scaling.factor.h0 = scaling.factor.h0,
        label =
          "Yuan-Bentler correction (Mplus variant)",
        trace.UGamma = trace.UGamma,
        trace.UGamma2 = as.numeric(NA),
        scaled.test.stat = TEST$standard$stat,
        scaled.test = TEST$standard$test
      )
  }

  TEST
}


lav_test_yuan_bentler_trace <- function(lavsamplestats = lavsamplestats,
                                        meanstructure = TRUE,
                                        A1.group = NULL,
                                        B1.group = NULL,
                                        Delta = NULL,
                                        Omega = NULL,
                                        E.inv = NULL,
                                        ug2.old.approach = FALSE,
                                        Satterthwaite = FALSE) {
  # we always assume a meanstructure (nope, not any longer, since 0.6)
  # meanstructure <- TRUE

  ngroups <- lavsamplestats@ngroups

  trace.h1 <- attr(Omega, "trace.h1")
  h1.ndat <- attr(Omega, "h1.ndat")

  if (ug2.old.approach || !Satterthwaite) {
    trace.UGamma <- numeric(ngroups)
    trace.UGamma2 <- numeric(ngroups)
    trace.h0 <- numeric(ngroups)

    for (g in 1:ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal

      A1 <- A1.group[[g]] * fg
      B1 <- B1.group[[g]] * fg
      DELTA <- Delta[[g]]
      Gamma.g <- Omega[[g]] / fg

      D.Einv.tD <- DELTA %*% tcrossprod(E.inv, DELTA)

      # trace.h1[g] <- sum( B1 * t( A1.inv ) )
      # fg cancels out: trace.h1[g] <- sum( fg*B1 * t( 1/fg*A1.inv ) )
      trace.h0[g] <- sum(B1 * D.Einv.tD)
      # trace.UGamma[g] <- trace.h1[g] - trace.h0[g]
      U <- A1 - A1 %*% D.Einv.tD %*% A1
      trace.UGamma[g] <- sum(U * Gamma.g)

      if (Satterthwaite) {
        UG <- U %*% Gamma.g
        trace.UGamma2[g] <- sum(UG * t(UG))
      }
    } # g
    trace.UGamma <- sum(trace.UGamma)
    attr(trace.UGamma, "h1") <- trace.h1
    attr(trace.UGamma, "h0") <- trace.h0
    attr(trace.UGamma, "h1.ndat") <- h1.ndat
    if (Satterthwaite) {
      attr(trace.UGamma, "trace.UGamma2") <- sum(trace.UGamma2)
    }
  } else {
    trace.UGamma <- trace.UGamma2 <- UG <- as.numeric(NA)
    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    # if(Satterthwaite) {

    A1.f <- A1.group
    for (g in 1:ngroups) {
      A1.f[[g]] <- A1.group[[g]] * fg[g]
    }
    A1.all <- lav_matrix_bdiag(A1.f)

    B1.f <- B1.group
    for (g in 1:ngroups) {
      B1.f[[g]] <- B1.group[[g]] * fg[g]
    }
    B1.all <- lav_matrix_bdiag(B1.f)

    Gamma.f <- Omega
    for (g in 1:ngroups) {
      Gamma.f[[g]] <- 1 / fg[g] * Omega[[g]]
    }
    Gamma.all <- lav_matrix_bdiag(Gamma.f)
    Delta.all <- do.call("rbind", Delta)

    D.Einv.tD <- Delta.all %*% tcrossprod(E.inv, Delta.all)

    trace.h0 <- sum(B1.all * D.Einv.tD)
    U.all <- A1.all - A1.all %*% D.Einv.tD %*% A1.all
    trace.UGamma <- sum(U.all * Gamma.all)

    attr(trace.UGamma, "h1") <- sum(trace.h1)
    attr(trace.UGamma, "h0") <- trace.h0
    attr(trace.UGamma, "h1.ndat") <- sum(h1.ndat)
    if (Satterthwaite) {
      UG <- U.all %*% Gamma.all
      trace.UGamma2 <- sum(UG * t(UG))
      attr(trace.UGamma, "trace.UGamma2") <- trace.UGamma2
    }

    # } else {
  }

  trace.UGamma
}

lav_test_yuan_bentler_mplus_trace <- function(lavsamplestats = NULL,
                                              A1.group = NULL,
                                              B1.group = NULL,
                                              B0.group = NULL,
                                              E.inv = NULL,
                                              meanstructure = TRUE) {
  # typical for Mplus:
  # - do NOT use the YB formula, but use an approximation
  #   relying  on A0 ~= Delta' A1 Delta and the same for B0
  #
  # NOTE: if A0 is based on the hessian, then A0 only approximates
  #       Delta' A1 Delta
  #
  # - always use h1.information = "unstructured"!!!

  ngroups <- lavsamplestats@ngroups

  trace.UGamma <- numeric(lavsamplestats@ngroups)
  trace.h1 <- numeric(lavsamplestats@ngroups)
  trace.h0 <- numeric(lavsamplestats@ngroups)
  h1.ndat <- numeric(lavsamplestats@ngroups)

  for (g in 1:lavsamplestats@ngroups) {
    # group weight
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal

    A1 <- A1.group[[g]]
    B1 <- B1.group[[g]]

    # mask independent 'fixed-x' variables
    zero.idx <- which(diag(A1) == 0)
    if (length(zero.idx) > 0L) {
      A1.inv <- matrix(0, nrow(A1), ncol(A1))
      a1 <- A1[-zero.idx, -zero.idx]
      a1.inv <- solve(a1)
      A1.inv[-zero.idx, -zero.idx] <- a1.inv
    } else {
      A1.inv <- solve(A1)
    }
    h1.ndat[g] <- ncol(A1) - length(zero.idx)

    # if data is complete, why not just A1 %*% Gamma?
    trace.h1[g] <- sum(B1 * t(A1.inv))
    trace.h0[g] <- fg * sum(B0.group[[g]] * t(E.inv))
    trace.UGamma[g] <- (trace.h1[g] - trace.h0[g])
  }

  # we take the sum here
  trace.UGamma <- sum(trace.UGamma)

  attr(trace.UGamma, "h1") <- trace.h1
  attr(trace.UGamma, "h0") <- trace.h0
  attr(trace.UGamma, "h1.ndat") <- h1.ndat

  trace.UGamma
}
