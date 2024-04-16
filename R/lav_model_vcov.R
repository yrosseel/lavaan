# bootstrap based NVCOV
lav_model_nvcov_bootstrap <- function(lavmodel = NULL,
                                      lavsamplestats = NULL,
                                      lavoptions = NULL,
                                      lavimplied = NULL,
                                      lavh1 = NULL,
                                      lavdata = NULL,
                                      lavcache = NULL,
                                      lavpartable = NULL) {
  # number of bootstrap draws
  if (!is.null(lavoptions$bootstrap)) {
    R <- lavoptions$bootstrap
  } else {
    R <- 1000L
  }

  boot.type <- "ordinary"
  if ("bollen.stine" %in% lavoptions$test) {
    boot.type <- "bollen.stine"
  }

  TEST <- NULL
  COEF <- lav_bootstrap_internal(
    object = NULL,
    lavmodel. = lavmodel,
    lavsamplestats. = lavsamplestats,
    lavpartable. = lavpartable,
    lavoptions. = lavoptions,
    lavdata. = lavdata,
    R = R,
    verbose = lavoptions$verbose,
    check.post = lavoptions$check.post,
    type = boot.type,
    FUN = ifelse(boot.type == "bollen.stine",
      "coeftest", "coef"
    )
  )
  # warn            = -1L)
  COEF.orig <- COEF

  # new in 0.6-12: always warn for failed and nonadmissible
  error.idx <- attr(COEF, "error.idx")
  nfailed <- length(error.idx) # zero if NULL
  if (nfailed > 0L && lavoptions$warn) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs failed or did not converge.", nfailed))
  }

  notok <- length(attr(COEF, "nonadmissible")) # zero if NULL
  if (notok > 0L && lavoptions$warn) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs resulted in nonadmissible solutions.", notok))
  }

  if (length(error.idx) > 0L) {
    # new in 0.6-13: we must still remove them!
    COEF <- COEF[-error.idx, , drop = FALSE]
    # this also drops the attributes
  }

  if (boot.type == "bollen.stine") {
    nc <- ncol(COEF)
    TEST <- COEF[, nc]
    COEF <- COEF[, -nc, drop = FALSE]
  }

  # FIXME: cov rescale? Yes for now
  nboot <- nrow(COEF)
  NVarCov <- lavsamplestats@ntotal * (cov(COEF) * (nboot - 1) / nboot)

  # save COEF and TEST (if any)
  attr(NVarCov, "BOOT.COEF") <- COEF.orig # including attributes
  attr(NVarCov, "BOOT.TEST") <- TEST

  NVarCov
}


# robust `sem' NVCOV (see Browne, 1984,  bentler & dijkstra 1985)
lav_model_nvcov_robust_sem <- function(lavmodel = NULL,
                                       lavsamplestats = NULL,
                                       lavdata = NULL,
                                       lavcache = NULL,
                                       lavimplied = NULL,
                                       lavh1 = NULL,
                                       lavoptions = NULL,
                                       use.ginv = FALSE) {
  # compute inverse of the expected(!) information matrix
  if (lavmodel@estimator == "ML" && lavoptions$mimic == "Mplus") {
    # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
    # - A1 is not based on Sigma.hat and Mu.hat,
    # but on lavsamplestats@cov and lavsamplestats@mean... ('unstructured')
    # - Gamma is not identical to what is used for WLS; closer to EQS
    # - N/N-1 bug in G11 for NVarCov (but not test statistic)
    # - we divide by N-1! (just like EQS)
    E.inv <- lav_model_information_expected_MLM(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      extra = TRUE,
      augmented = TRUE,
      inverted = TRUE,
      use.ginv = use.ginv
    )
  } else {
    E.inv <- lav_model_information(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions,
      extra = TRUE,
      augmented = TRUE,
      inverted = TRUE,
      use.ginv = use.ginv
    )
  }

  # check if E.inv is ok
  if (inherits(E.inv, "try-error")) {
    return(E.inv)
  }

  Delta <- attr(E.inv, "Delta")
  WLS.V <- attr(E.inv, "WLS.V")

  # Gamma
  Gamma <- lavsamplestats@NACOV
  if (lavmodel@estimator == "ML" &&
    lavoptions$mimic == "Mplus" && !lavsamplestats@NACOV.user) {
    # 'fix' G11 part of Gamma (NOTE: this is NOT needed for SB test
    # statistic
    for (g in 1:lavsamplestats@ngroups) {
      gg1 <- (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@nobs[[g]]
      if (lavmodel@conditional.x) {
        nvar <- NCOL(lavsamplestats@res.cov[[g]])
      } else {
        nvar <- NCOL(lavsamplestats@cov[[g]])
      }
      G11 <- Gamma[[g]][1:nvar, 1:nvar, drop = FALSE]
      Gamma[[g]][1:nvar, 1:nvar] <- G11 * gg1
    } # g
  }


  tDVGVD <- matrix(0, ncol = ncol(E.inv), nrow = nrow(E.inv))
  for (g in 1:lavsamplestats@ngroups) {
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    if (lavoptions$mimic == "Mplus") {
      fg1 <- (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@ntotal
    } else {
      # from 0.6 onwards, we use fg1 == fg, to be more consistent with
      # lav_test()
      fg1 <- fg
    }
    # fg twice for WLS.V, 1/fg1 once for GaMMA
    # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
    # t(Delta) * WLS.V %*% Gamma %*% WLS.V %*% Delta
    if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
      # diagonal weight matrix
      WD <- WLS.V[[g]] * Delta[[g]]
    } else {
      # full weight matrix
      WD <- WLS.V[[g]] %*% Delta[[g]]
    }
    tDVGVD <- tDVGVD + fg * fg / fg1 * crossprod(WD, Gamma[[g]] %*% WD)
  } # g
  NVarCov <- (E.inv %*% tDVGVD %*% E.inv)

  # to be reused by lav_test()
  attr(NVarCov, "Delta") <- Delta

  if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    # only when same type of information is used # new in 0.6-6
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "WLS.V") <- WLS.V
  }

  NVarCov
}

lav_model_nvcov_robust_sandwich <- function(lavmodel = NULL,
                                            lavsamplestats = NULL,
                                            lavdata = NULL,
                                            lavoptions = NULL,
                                            lavimplied = NULL,
                                            lavh1 = NULL,
                                            lavcache = NULL,
                                            use.ginv = FALSE) {
  # sandwich estimator: A.inv %*% B %*% t(A.inv)
  # where A.inv == E.inv
  #       B == outer product of case-wise scores

  # inverse observed/expected information matrix
  E.inv <- lav_model_information(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavcache = lavcache,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    lavoptions = lavoptions,
    extra = FALSE,
    augmented = TRUE,
    inverted = TRUE,
    use.ginv = use.ginv
  )

  # check if E.inv is ok
  if (inherits(E.inv, "try-error")) {
    return(E.inv)
  }

  # new in 0.6-6, check for h1.information.meat
  lavoptions2 <- lavoptions
  if (!is.null(lavoptions$information.meat)) {
    lavoptions2$information <- lavoptions$information.meat
  }
  if (!is.null(lavoptions$h1.information.meat)) {
    lavoptions2$h1.information <- lavoptions$h1.information.meat
  }

  # outer product of case-wise scores
  B0 <-
    lav_model_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions2,
      extra = TRUE,
      check.pd = FALSE,
      augmented = FALSE,
      inverted = FALSE,
      use.ginv = use.ginv
    )

  # compute sandwich estimator
  NVarCov <- E.inv %*% B0 %*% E.inv

  attr(NVarCov, "B0.group") <- attr(B0, "B0.group")

  if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    # only when same type of information is used # new in 0.6-6
    attr(NVarCov, "E.inv") <- E.inv
  }

  NVarCov
}

# two stage
# - two.stage: Gamma = I_1^{-1}
# - robust.two.stage: Gamma = incomplete Gamma (I_1^{-1} J_1 I_1^{-1})
# where I_1 and J_1 are based on the (saturated) model h1
# (either unstructured, or structured)
#
# references:
#
# - Savalei \& Bentler (2009) eq (6) for se = "two.stage"
# - Savalei \& Falk (2014) eq  (3)   for se = "robust.two.stage"
# - Yuan \& Bentler (2000)
lav_model_nvcov_two_stage <- function(lavmodel = NULL,
                                      lavsamplestats = NULL,
                                      lavoptions = NULL,
                                      lavimplied = NULL,
                                      lavh1 = NULL,
                                      lavdata = NULL,
                                      use.ginv = FALSE) {
  # expected OR observed, depending on lavoptions$information
  if (is.null(lavoptions) && is.null(lavoptions$information[1])) {
    lavoptions <- list(
      information = "observed",
      observed.information = "h1",
      h1.information = "structured"
    )
  }

  # restrictions:
  # only works if:
  # - information is expected,
  # - or information is observed but with observed.information == "h1"
  if (lavoptions$information[1] == "observed" &&
    lavoptions$observed.information[1] != "h1") {
    lav_msg_stop(
      gettext("two.stage + observed information currently only works
              with observed.information = 'h1'"))
  }
  # no weights (yet)
  if (!is.null(lavdata@weights[[1]])) {
    lav_msg_stop(gettext("two.stage + sampling.weights is not supported yet"))
  }
  # no fixed.x (yet)
  # if(!is.null(lavsamplestats@x.idx) &&
  #   length(lavsamplestats@x.idx[[1]]) > 0L) {
  #    lav_msg_stop(gettext("two.stage + fixed.x = TRUE is not supported yet"))
  # }


  # information matrix
  E.inv <- lav_model_information(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    extra = TRUE,
    augmented = TRUE,
    inverted = TRUE,
    use.ginv = use.ginv
  )
  Delta <- attr(E.inv, "Delta")
  WLS.V <- attr(E.inv, "WLS.V") # this is 'H' or 'A1' in the literature
  attr(E.inv, "Delta") <- NULL
  attr(E.inv, "WLS.V") <- NULL

  # check if E.inv is ok
  if (inherits(E.inv, "try-error")) {
    return(E.inv)
  }

  # check WLS.V = A1
  if (is.null(WLS.V)) {
    lav_msg_stop(gettext("WLS.V/H/A1 is NULL, observed.information = hessian?"))
  }

  # Gamma
  Gamma <- vector("list", length = lavsamplestats@ngroups)

  # handle multiple groups
  tDVGVD <- matrix(0, ncol = ncol(E.inv), nrow = nrow(E.inv))
  for (g in 1:lavsamplestats@ngroups) {
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    # fg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@ntotal
    fg1 <- fg
    # fg twice for WLS.V, 1/fg1 once for GaMMA
    # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
    # t(Delta) * WLS.V %*% Gamma %*% WLS.V %*% Delta
    WD <- WLS.V[[g]] %*% Delta[[g]]

    # to compute (incomplete) GAMMA, should we use
    # structured or unstructured mean/sigma?
    #
    # we use the same setting as to compute 'H' (the h1 information matrix)
    # so that at Omega = H if data is complete
    if (lavoptions$h1.information[1] == "unstructured") {
      MU <- lavsamplestats@missing.h1[[g]]$mu
      SIGMA <- lavsamplestats@missing.h1[[g]]$sigma
    } else {
      MU <- lavimplied$mean[[g]]
      SIGMA <- lavimplied$cov[[g]]
    }

    # compute 'Gamma' (or Omega.beta)
    if (lavoptions$se == "two.stage") {
      # this is Savalei & Bentler (2009)
      if (lavoptions$information[1] == "expected") {
        Info <- lav_mvnorm_missing_information_expected(
          Y = lavdata@X[[g]], Mp = lavdata@Mp[[g]],
          wt = lavdata@weights[[g]],
          Mu = MU, Sigma = SIGMA,
          x.idx = lavsamplestats@x.idx[[g]]
        )
      } else {
        Info <- lav_mvnorm_missing_information_observed_samplestats(
          Yp = lavsamplestats@missing[[g]],
          # wt not needed
          Mu = MU, Sigma = SIGMA,
          x.idx = lavsamplestats@x.idx[[g]]
        )
      }
      Gamma[[g]] <- lav_matrix_symmetric_inverse(Info)
    } else { # we assume "robust.two.stage"
      # NACOV is here incomplete Gamma
      # Savalei & Falk (2014)
      #
      if (length(lavdata@cluster) > 0L) {
        cluster.idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster.idx <- NULL
      }
      Gamma[[g]] <- lav_mvnorm_missing_h1_omega_sw(
        Y =
          lavdata@X[[g]], Mp = lavdata@Mp[[g]],
        Yp = lavsamplestats@missing[[g]],
        wt = lavdata@weights[[g]],
        cluster.idx = cluster.idx,
        Mu = MU, Sigma = SIGMA,
        x.idx = lavsamplestats@x.idx[[g]],
        information = lavoptions$information[1]
      )
    }

    # compute
    tDVGVD <- tDVGVD + fg * fg / fg1 * crossprod(WD, Gamma[[g]] %*% WD)
  } # g

  NVarCov <- (E.inv %*% tDVGVD %*% E.inv)

  # to be reused by lavaanTest
  attr(NVarCov, "Delta") <- Delta
  attr(NVarCov, "Gamma") <- Gamma
  if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    # only when same type of information is used # new in 0.6-6
    attr(NVarCov, "E.inv") <- E.inv
    attr(NVarCov, "WLS.V") <- WLS.V
  }

  NVarCov
}



lav_model_vcov <- function(lavmodel = NULL,
                           lavsamplestats = NULL,
                           lavoptions = NULL,
                           lavdata = NULL,
                           lavpartable = NULL,
                           lavcache = NULL,
                           lavimplied = NULL,
                           lavh1 = NULL,
                           use.ginv = FALSE) {
  likelihood <- lavoptions$likelihood
  information <- lavoptions$information[1] # first one is for vcov
  se <- lavoptions$se
  verbose <- lavoptions$verbose
  mimic <- lavoptions$mimic

  # special cases
  if (se == "none" || se == "external" || se == "twostep") {
    return(matrix(0, 0, 0))
  }

  if (se == "standard") {
    NVarCov <- lav_model_information(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions,
      extra = FALSE,
      augmented = TRUE,
      inverted = TRUE,
      use.ginv = use.ginv
    )
  } else if (se == "first.order") {
    NVarCov <-
      lav_model_information_firstorder(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        extra = TRUE,
        check.pd = FALSE,
        augmented = TRUE,
        inverted = TRUE,
        use.ginv = use.ginv
      )
  } else if (se == "robust.sem" || se == "robust.cluster.sem") {
    NVarCov <-
      lav_model_nvcov_robust_sem(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavcache = lavcache,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        use.ginv = use.ginv
      )
  } else if (se == "robust.huber.white" || se == "robust.cluster") {
    NVarCov <-
      lav_model_nvcov_robust_sandwich(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        use.ginv = use.ginv
      )
  } else if (se %in% c("two.stage", "robust.two.stage")) {
    NVarCov <-
      lav_model_nvcov_two_stage(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavoptions = lavoptions,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        use.ginv = use.ginv
      )
  } else if (se == "bootstrap") {
    NVarCov <- try(
      lav_model_nvcov_bootstrap(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavoptions = lavoptions,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavcache = lavcache,
        lavpartable = lavpartable
      ),
      silent = TRUE
    )
  } else {
    lav_msg_warn(gettextf("unknown se type: %s", se))
  }

  if (!inherits(NVarCov, "try-error")) {
    # denominator!
    if (lavmodel@estimator %in% c("ML", "PML", "FML") &&
      likelihood == "normal") {
      if (lavdata@nlevels == 1L) {
        N <- lavsamplestats@ntotal
        # new in 0.6-9 (to mimic method="lm" in effectLite)
        # special case: univariate regression in each group
        if (lavoptions$mimic == "lm" &&
          .hasSlot(lavmodel, "modprop") &&
          all(lavmodel@modprop$uvreg)) {
          N <- sum(unlist(lavsamplestats@nobs) -
            (unlist(lavmodel@modprop$nexo) + 1L))
          # always adding the intercept (for now)
        }
      } else {
        # total number of clusters (over groups)
        N <- 0
        for (g in 1:lavsamplestats@ngroups) {
          N <- N + lavdata@Lp[[g]]$nclusters[[2]]
        }
      }
    } else {
      N <- lavsamplestats@ntotal - lavsamplestats@ngroups
    }

    VarCov <- 1 / N * NVarCov

    # check if VarCov is pd -- new in 0.6-2
    # mostly important if we have (in)equality constraints (MASS::ginv!)
    if (.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
      # do nothing
    } else if (!is.null(lavoptions$check.vcov) && lavoptions$check.vcov) {
      eigvals <- eigen(VarCov,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      # correct for (in)equality constraints
      neq <- 0L
      niq <- 0L
      if (nrow(lavmodel@con.jac) > 0L) {
        ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
        cin.idx <- attr(lavmodel@con.jac, "cin.idx")
        ina.idx <- attr(lavmodel@con.jac, "inactive.idx")
        if (length(ceq.idx) > 0L) {
          neq <- qr(lavmodel@con.jac[ceq.idx, , drop = FALSE])$rank
        }
        if (length(cin.idx) > 0L) {
          niq <- length(cin.idx) - length(ina.idx) # only active
        }
        # total number of relevant constraints
        neiq <- neq + niq
        if (neiq > 0L) {
          eigvals <- rev(eigvals)[-seq_len(neiq)]
        }
      }
      min.val <- min(eigvals)
      # if(any(eigvals < -1 * sqrt(.Machine$double.eps)) &&
      if (min.val < .Machine$double.eps^(3 / 4) && lavoptions$warn) {
        # VarCov.chol <- suppressWarnings(try(chol(VarCov,
        #                   pivot = TRUE), silent = TRUE))
        # VarCov.rank <- attr(VarCov.chol, "rank")
        # VarCov.pivot <- attr(VarCov.chol, "pivot")
        # VarCov.badidx <- VarCov.pivot[ VarCov.rank + 1L ]
        # pt.idx <- which(lavpartable$free == VarCov.badidx)
        # par.string <- paste(lavpartable$lhs[pt.idx],
        #                    lavpartable$op[ pt.idx],
        #                    lavpartable$rhs[pt.idx])
        # if(lavdata@ngroups > 1L) {
        #    par.string <- paste0(par.string, " in group ",
        #                         lavpartable$group[pt.idx])
        # }
        # if(lavdata@nlevels > 1L) {
        #    par.string <- paste0(par.string, " in level ",
        #                         lavpartable$level[pt.idx])
        # }

        if (min.val > 0) {
          lav_msg_warn(
            gettextf("The variance-covariance matrix of the estimated
                    parameters (vcov) does not appear to be positive
                    definite! The smallest eigenvalue (= %e) is close
                     to zero. This may be a symptom that the model is
                     not identified.", min(min.val)))
        } else {
          lav_msg_warn(
            gettextf("The variance-covariance matrix of the estimated parameters
                     (vcov) does not appear to be positive definite! The smallest
                     eigenvalue (= %e) is smaller than zero. This may be a
                     symptom that the model is not identified.",
                     min(min.val)))
        }
      }
    }
  } else {
    if (lavoptions$warn) {
      lav_msg_warn(
       gettext("Could not compute standard errors! The information matrix
               could not be inverted. This may be a symptom that the model
               is not identified.")
      )
    }
    VarCov <- NULL
  } # could not invert

  VarCov
}

lav_model_vcov_se <- function(lavmodel, lavpartable, VCOV = NULL,
                              BOOT = NULL) {
  # 0. special case
  if (is.null(VCOV)) {
    se <- rep(as.numeric(NA), lavmodel@nx.user)
    se[lavpartable$free == 0L] <- 0.0
    return(se)
  }

  # 1. free parameters only
  x.var <- diag(VCOV)
  # check for negative values (what to do: NA or 0.0?)
  x.var[x.var < 0] <- as.numeric(NA)
  x.se <- sqrt(x.var)
  if (.hasSlot(lavmodel, "ceq.simple.only") && lavmodel@ceq.simple.only) {
    GLIST <- lav_model_x2GLIST(
      lavmodel = lavmodel, x = x.se,
      type = "unco"
    )
  } else {
    GLIST <- lav_model_x2GLIST(
      lavmodel = lavmodel, x = x.se,
      type = "free"
    )
  }

  # se for full parameter table, but with 0.0 entries for def/ceq/cin
  # elements
  se <- lav_model_get_parameters(
    lavmodel = lavmodel, GLIST = GLIST,
    type = "user", extra = FALSE
  )


  # 2. fixed parameters -> se = 0.0
  se[which(lavpartable$free == 0L)] <- 0.0


  # 3. defined parameters:
  def.idx <- which(lavpartable$op == ":=")
  if (length(def.idx) > 0L) {
    if (!is.null(BOOT)) {
      # we must remove the NA rows (and hope we have something left)
      error.idx <- attr(BOOT, "error.idx")
      if (length(error.idx) > 0L) {
        BOOT <- BOOT[-error.idx, , drop = FALSE] # drops attributes
      }

      BOOT.def <- apply(BOOT, 1L, lavmodel@def.function)
      if (length(def.idx) == 1L) {
        BOOT.def <- as.matrix(BOOT.def)
      } else {
        BOOT.def <- t(BOOT.def)
      }
      def.cov <- cov(BOOT.def)
    } else {
      # regular delta method
      x <- lav_model_get_parameters(lavmodel = lavmodel, type = "free")
      JAC <- try(lav_func_jacobian_complex(func = lavmodel@def.function, x = x),
        silent = TRUE
      )
      if (inherits(JAC, "try-error")) { # eg. pnorm()
        JAC <- lav_func_jacobian_simple(func = lavmodel@def.function, x = x)
      }
      if (.hasSlot(lavmodel, "ceq.simple.only") &&
        lavmodel@ceq.simple.only) {
        JAC <- JAC %*% t(lavmodel@ceq.simple.K)
      }
      def.cov <- JAC %*% VCOV %*% t(JAC)
    }
    # check for negative se's
    diag.def.cov <- diag(def.cov)
    diag.def.cov[diag.def.cov < 0] <- as.numeric(NA)
    se[def.idx] <- sqrt(diag.def.cov)
  }

  se
}
