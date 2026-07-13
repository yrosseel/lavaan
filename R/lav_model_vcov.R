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
  if (!is.null(lavoptions$bootstrap) &&
      !is.null(lavoptions$bootstrap$R)) {
    r <- as.integer(lavoptions$bootstrap$R)
    stopifnot(r > 0L)
  } else {
    r <- 1000L
  }

  boot_type <- "ordinary"
  if ("bollen.stine" %in% lavoptions$test) {
    boot_type <- "bollen.stine"
  }

  test <- NULL
  coef_1 <- lav_bootstrap_internal(
    object = NULL,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavpartable = lavpartable,
    lavoptions = lavoptions,
    lavdata = lavdata,
    r = r,
    check_post = lavoptions$check.post,
    type = boot_type,
    fun = ifelse(boot_type == "bollen.stine",
      "coeftest", "coef"
    )
  )
  # warn            = -1L)
  coef_orig <- coef_1

  # new in 0.6-12: always warn for failed and nonadmissible
  error_idx <- attr(coef_1, "error.idx")
  nfailed <- length(error_idx) # zero if NULL
  if (nfailed > 0L) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs failed or did not converge.", nfailed))
  }

  notok <- length(attr(coef_1, "nonadmissible")) # zero if NULL
  if (notok > 0L) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs resulted in nonadmissible solutions.", notok))
  }

  if (length(error_idx) > 0L) {
    # new in 0.6-13: we must still remove them!
    coef_1 <- coef_1[-error_idx, , drop = FALSE]
    # this also drops the attributes
  }

  if (boot_type == "bollen.stine") {
    nc <- ncol(coef_1)
    test <- coef_1[, nc]
    coef_1 <- coef_1[, -nc, drop = FALSE]
    # the last column is the (bollen.stine) test statistic, not a
    # coefficient; remove it from BOOT.COEF too (it is stored separately
    # in BOOT.TEST), but preserve the attributes (error.idx, nonadmissible,
    # seed, ...) that subsetting would otherwise drop
    keep_attr <- attributes(coef_orig)
    keep_attr[c("dim", "dimnames")] <- NULL
    coef_orig <- coef_orig[, -nc, drop = FALSE]
    attributes(coef_orig) <- c(attributes(coef_orig), keep_attr)
  }

  # new in 0.6-20: check for outliers, ie big difference between sd() and mad()
  # see github issue 347
  sd_mad_ratio <- (apply(coef_1, 2,  sd, na.rm = TRUE) /
                   apply(coef_1, 2, mad, na.rm = TRUE))
  crit_ratio <- 5
  if (any(sd_mad_ratio > crit_ratio)) {
    names_1 <- lav_pt_labels(lavpartable, type = "free")
    params_w_outliers <- paste(names_1[sd_mad_ratio > crit_ratio],
                               collapse = " ")
    lav_msg_warn(gettextf(
      "The following bootstrapped free parameters have a high (>5)
      ratio of standard deviation to median absolute deviation: %s.
      P-values and confidence intervals may not match.", params_w_outliers))
  }

  # FIXME: cov rescale? Yes for now
  nboot <- nrow(coef_1)
  # The bootstrap directly estimates Cov(coef), which IS the target VCOV.
  # lav_model_vcov() will afterwards divide this NVCOV by 'n' to obtain the
  # VCOV, so we scale by the same 'n' here, making the two cancel. For
  # multilevel models, lavaan divides by the (total) number of clusters --
  # not by ntotal -- so we must use the same denominator here (otherwise the
  # bootstrap standard errors would be off by a factor sqrt(ntotal/nclusters)).
  # This mirrors the 'n' computation in lav_model_vcov().
  likelihood <- lavoptions$likelihood
  if (lavdata@nlevels > 1L &&
      lavmodel@estimator %in% c("ML", "PML", "FML") &&
      !is.null(likelihood) && likelihood == "normal") {
    n <- 0
    for (g in seq_len(lavsamplestats@ngroups)) {
      n <- n + lavdata@Lp[[g]]$nclusters[[2]]
    }
  } else {
    n <- lavsamplestats@ntotal
  }
  nvar_cov <- n * (cov(coef_1) * (nboot - 1) / nboot)

  # ceq.simple.only models (e.g. estimator = "IV" with equality constraints):
  # the bootstrap draws are in the compact (free) space (one column per
  # non-collapsed free parameter), but the rest of the vcov machinery
  # (lav_model_vcov_se) expects the 'unco' space (one row/col per parameter,
  # with shared rows/cols for parameters constrained equal). Expand via
  # ceq.simple.K so that, e.g., loadings constrained equal get identical
  # standard errors. The BOOT.COEF stored below stays in the compact space
  # (the bootstrap CI code in lavParameterEstimates() expects that).
  if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L &&
      ncol(nvar_cov) == ncol(lavmodel@ceq.simple.K)) {
    tmp_k <- lavmodel@ceq.simple.K
    nvar_cov <- tmp_k %*% nvar_cov %*% t(tmp_k)
  }

  # save COEF and TEST (if any)
  attr(nvar_cov, "BOOT.COEF") <- coef_orig # including attributes
  attr(nvar_cov, "BOOT.TEST") <- test

  nvar_cov
}


# robust `sem' NVCOV (see Browne, 1984,  bentler & dijkstra 1985)
lav_model_nvcov_robust_sem <- function(lavmodel = NULL,
                                       lavsamplestats = NULL,
                                       lavdata = NULL,
                                       lavcache = NULL,
                                       lavimplied = NULL,
                                       lavh1 = NULL,
                                       lavoptions = NULL,
                                       use_ginv = FALSE,
                                       attr_delta = TRUE,
                                       attr_t_dvgvd = FALSE,
                                       attr_e_inv = FALSE,
                                       attr_wls_v = FALSE) {
  # compute inverse of the expected(!) information matrix
  if (lavmodel@estimator == "ML" && lavoptions$information.expected.mplus) {
    # YR - 11 aug 2010 - what Mplus seems to do is (see Muthen apx 4 eq102)
    # - A1 is not based on Sigma.hat and Mu.hat,
    # but on lavsamplestats@cov and lavsamplestats@mean... ('unstructured')
    # - gamma is not identical to what is used for WLS; closer to EQS
    # - N/N-1 bug in G11 for NVarCov (but not test statistic)
    # - we divide by N-1! (just like EQS)
    e_inv <- lav_model_info_expected_mlm(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      extra = TRUE,
      augmented = TRUE,
      inverted = TRUE,
      use_ginv = use_ginv
    )
  } else {
    e_inv <- lav_model_info(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions,
      extra = TRUE,
      augmented = TRUE,
      inverted = TRUE,
      use_ginv = use_ginv
    )
  }

  # check if E.inv is ok
  if (inherits(e_inv, "try-error")) {
    return(e_inv)
  }

  delta <- attr(e_inv, "Delta")
  wls_v <- attr(e_inv, "WLS.V")
  attr(e_inv, "Delta") <- NULL
  attr(e_inv, "WLS.V") <- NULL
  # gamma
  gamma <- lavsamplestats@NACOV
  if (lavmodel@estimator == "ML" &&
    lavoptions$gamma.vcov.mplus && !lavsamplestats@NACOV.user) {
    # 'fix' G11 part of gamma (NOTE: this is NOT needed for SB test
    # statistic
    for (g in 1:lavsamplestats@ngroups) {
      gg1 <- (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@nobs[[g]]
      if (lavmodel@conditional.x) {
        nvar <- NCOL(lavsamplestats@res.cov[[g]])
      } else {
        nvar <- NCOL(lavsamplestats@cov[[g]])
      }
      g11 <- gamma[[g]][1:nvar, 1:nvar, drop = FALSE]
      gamma[[g]][1:nvar, 1:nvar] <- g11 * gg1
    } # g
  }


  t_dvgvd <- matrix(0, ncol = ncol(e_inv), nrow = nrow(e_inv))
  for (g in 1:lavsamplestats@ngroups) {
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    if (lavoptions$gamma.vcov.mplus) {
      fg1 <- (lavsamplestats@nobs[[g]] - 1) / lavsamplestats@ntotal
    } else {
      # from 0.6 onwards, we use fg1 == fg, to be more consistent with
      # lav_test()
      fg1 <- fg
    }
    # fg twice for WLS.V, 1/fg1 once for GaMMA
    # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
    # t(Delta) * WLS.V %*% gamma %*% WLS.V %*% Delta
    if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
      # diagonal weight matrix
      wd <- wls_v[[g]] * delta[[g]]
    } else {
      # full weight matrix
      wd <- wls_v[[g]] %*% delta[[g]]
    }
    t_dvgvd <- t_dvgvd + fg * fg / fg1 * crossprod(wd, gamma[[g]] %*% wd)
  } # g
  nvar_cov <- (e_inv %*% t_dvgvd %*% e_inv)

  # to be reused by lav_test()
  if (attr_delta) {
    attr(nvar_cov, "Delta") <- delta
  }
  # for twostep.robust in sam()
  if (attr_t_dvgvd) {
    attr(nvar_cov, "tDVGVD") <- t_dvgvd
  }

  if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    # only when same type of information is used # new in 0.6-6
    attr(nvar_cov, "E.inv") <- e_inv
    attr(nvar_cov, "WLS.V") <- wls_v
  }

  # user override
  if (attr_e_inv && is.null(attr(nvar_cov, "E.inv"))) {
    attr(nvar_cov, "E.inv") <- e_inv
  }
  if (attr_wls_v && is.null(attr(nvar_cov, "WLS.V"))) {
    attr(nvar_cov, "WLS.V") <- wls_v
  }

  nvar_cov
}

lav_model_nvcov_robust_sandwich <- function(lavmodel = NULL,    # nolint
                                            lavsamplestats = NULL,
                                            lavdata = NULL,
                                            lavoptions = NULL,
                                            lavimplied = NULL,
                                            lavh1 = NULL,
                                            lavcache = NULL,
                                            use_ginv = FALSE) {
  # sandwich estimator: A.inv %*% B %*% t(A.inv)
  # where A.inv == E.inv
  #       B == outer product of case-wise scores

  # inverse observed/expected information matrix
  e_inv <- lav_model_info(
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
    use_ginv = use_ginv
  )

  # check if E.inv is ok
  if (inherits(e_inv, "try-error")) {
    return(e_inv)
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
  b0 <-
    lav_model_info_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavoptions = lavoptions2,
      extra = TRUE,
      check_pd = FALSE,
      augmented = FALSE,
      inverted = FALSE,
      use_ginv = use_ginv
    )

  # compute sandwich estimator
  nvar_cov <- e_inv %*% b0 %*% e_inv

  attr(nvar_cov, "B0.group") <- attr(b0, "B0.group")

  if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    # only when same type of information is used # new in 0.6-6
    attr(nvar_cov, "E.inv") <- e_inv
  }

  nvar_cov
}

# two stage
# - two.stage: gamma = I_1^{-1}
# - robust.two.stage: gamma = incomplete gamma (I_1^{-1} J_1 I_1^{-1})
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
                                      use_ginv = FALSE) {
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
  e_inv <- lav_model_info(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    extra = TRUE,
    augmented = TRUE,
    inverted = TRUE,
    use_ginv = use_ginv
  )
  delta <- attr(e_inv, "Delta")
  wls_v <- attr(e_inv, "WLS.V") # this is 'H' or 'A1' in the literature
  attr(e_inv, "Delta") <- NULL
  attr(e_inv, "WLS.V") <- NULL

  # check if E.inv is ok
  if (inherits(e_inv, "try-error")) {
    return(e_inv)
  }

  # check WLS.V = A1
  if (is.null(wls_v)) {
    lav_msg_stop(gettext("WLS.V/H/A1 is NULL, observed.information = hessian?"))
  }

  # gamma
  gamma <- vector("list", length = lavsamplestats@ngroups)

  # handle multiple groups
  t_dvgvd <- matrix(0, ncol = ncol(e_inv), nrow = nrow(e_inv))
  for (g in 1:lavsamplestats@ngroups) {
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    # fg1 <- (lavsamplestats@nobs[[g]]-1)/lavsamplestats@ntotal
    fg1 <- fg
    # fg twice for WLS.V, 1/fg1 once for GaMMA
    # if fg==fg1, there would be only one fg, as in Satorra 1999 p.8
    # t(Delta) * WLS.V %*% gamma %*% WLS.V %*% Delta
    if (lavmodel@estimator == "DWLS" || lavmodel@estimator == "ULS") {
      # diagonal weight matrix (stored as a vector)
      wd <- wls_v[[g]] * delta[[g]]
    } else {
      # full weight matrix
      wd <- wls_v[[g]] %*% delta[[g]]
    }

    # conditional.x: the stage-1 'saturated' parameters live in the
    # conditional metric (vec(Beta), vech(res.cov)); compute Omega with
    # the lav_mvreg_mi_* kernels (x is complete). Note: with complete x,
    # the joint information is block-diagonal between the conditional
    # parameters and the x moments, so the conditional Omega is the
    # correct stage-1 ACOV.
    if (lavmodel@conditional.x) {
      aux_g2 <- if (length(lavdata@aux) >= g) lavdata@aux[[g]] else NULL
      if (!is.null(aux_g2) && NCOL(aux_g2) > 0L) {
        lav_msg_stop(gettext(
          "auxiliary variables (aux =) are not supported (yet) for
          two-stage standard errors when conditional.x = TRUE."))
      }
      if (lavoptions$h1.information[1] == "unstructured") {
        res_int_g <- lavh1$implied$res.int[[g]]
        res_slopes_g <- lavh1$implied$res.slopes[[g]]
        res_cov_g <- lavh1$implied$res.cov[[g]]
      } else {
        res_int_g <- lavimplied$res.int[[g]]
        res_slopes_g <- lavimplied$res.slopes[[g]]
        res_cov_g <- lavimplied$res.cov[[g]]
      }

      if (lavoptions$se == "two.stage") {
        # Savalei & Bentler (2009), conditional metric
        if (lavoptions$information[1] == "expected") {
          info <- lav_mvreg_mi_info_expected(
            yp = lavsamplestats@missing[[g]],
            res_cov = res_cov_g
          )
        } else {
          info <- lav_mvreg_mi_information_observed_samplestats(
            yp = lavsamplestats@missing[[g]],
            res_int = res_int_g, res_slopes = res_slopes_g,
            res_cov = res_cov_g
          )
        }
        omega_g <- lav_mat_sym_inverse(info)
      } else { # robust.two.stage
        # Savalei & Falk (2014), conditional metric
        if (length(lavdata@cluster) > 0L) {
          cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
        } else {
          cluster_idx <- NULL
        }
        omega_g <- lav_mvreg_mi_h1_omega_sw(
          y = lavdata@X[[g]],
          exo = lavdata@eXo[[g]],
          mp = lavdata@Mp[[g]],
          yp = lavsamplestats@missing[[g]],
          wt = lavdata@weights[[g]],
          cluster_idx = cluster_idx,
          res_int = res_int_g, res_slopes = res_slopes_g,
          res_cov = res_cov_g,
          information = lavoptions$information[1]
        )
      }
      gamma[[g]] <- omega_g

      # compute
      t_dvgvd <- t_dvgvd + fg * fg / fg1 * crossprod(wd, gamma[[g]] %*% wd)
      next
    }

    # to compute (incomplete) GAMMA, should we use
    # structured or unstructured mean/sigma?
    #
    # we use the same setting as to compute 'H' (the h1 information matrix)
    # so that at Omega = H if data is complete
    if (lavoptions$h1.information[1] == "unstructured") {
      #MU <- lavsamplestats@missing.h1[[g]]$mu
      #SIGMA <- lavsamplestats@missing.h1[[g]]$sigma
      mu <- lavh1$implied$mean[[g]]
      sigma_1 <- lavh1$implied$cov[[g]]
    } else {
      mu <- lavimplied$mean[[g]]
      sigma_1 <- lavimplied$cov[[g]]
    }

    # auxiliary variables (aux=): with two-stage estimation, the stage-1
    # saturated moments are estimated over the augmented set [model vars, aux],
    # so the stage-1 ACOV (Omega) must be computed over that augmented set and
    # then reduced to the model-variable moments (Savalei & Bentler, 2009).
    # We only do this when there are no fixed-x covariates (the common case).
    aux_g <- if (length(lavdata@aux) >= g) lavdata@aux[[g]] else NULL
    use_aux <- !is.null(aux_g) && NCOL(aux_g) > 0L &&
      length(lavsamplestats@x.idx[[g]]) == 0L
    if (use_aux) {
      y_g <- cbind(lavdata@X[[g]], aux_g)
      mp_g <- lav_data_mi_patterns(y_g)
      mu_g <- lavsamplestats@missing.h1[[g]]$mu.aug
      sigma_g <- lavsamplestats@missing.h1[[g]]$sigma.aug
      if (is.null(mu_g) || is.null(sigma_g)) {
        em_aug <- lav_mvn_mi_h1_est_moments(y_g,
          mp = mp_g, max_iter = lavoptions$em.h1.args$max_iter,
          tol = lavoptions$em.h1.args$tol,
          non_pd_action = lavoptions$em.h1.args$non_pd_action,
          non_pd_tol = lavoptions$em.h1.args$non_pd_tol,
          acceleration = lavoptions$em.h1.args$acceleration
        )
        mu_g <- em_aug$Mu
        sigma_g <- em_aug$Sigma
      }
      midx <- lav_aux_moment_idx(
        p_model = NCOL(lavdata@X[[g]]), p_aug = NCOL(y_g)
      )
    } else {
      y_g <- lavdata@X[[g]]
      mp_g <- lavdata@Mp[[g]]
      mu_g <- mu
      sigma_g <- sigma_1
      midx <- NULL
    }

    # compute 'gamma' (or Omega.beta)
    if (lavoptions$se == "two.stage") {
      # this is Savalei & Bentler (2009)
      if (lavoptions$information[1] == "expected") {
        info <- lav_mvn_mi_info_expected(
          y = y_g, mp = mp_g,
          wt = lavdata@weights[[g]],
          mu = mu_g, sigma_1 = sigma_g,
          x_idx = lavsamplestats@x.idx[[g]]
        )
      } else {
        info <- lav_mvnorm_missing_information_observed_samplestats(
          yp = if (use_aux) {
            lav_samp_mi_patterns(y = y_g, mp = mp_g)
          } else {
            lavsamplestats@missing[[g]]
          },
          # wt not needed
          mu = mu_g, sigma_1 = sigma_g,
          x_idx = lavsamplestats@x.idx[[g]]
        )
      }
      omega_g <- lav_mat_sym_inverse(info)
    } else { # we assume "robust.two.stage"
      # NACOV is here incomplete gamma
      # Savalei & Falk (2014)
      #
      if (length(lavdata@cluster) > 0L) {
        cluster_idx <- lavdata@Lp[[g]]$cluster.idx[[2]]
      } else {
        cluster_idx <- NULL
      }
      omega_g <- lav_mvn_mi_h1_omega_sw(
        y = y_g,
        mp = mp_g,
        yp = if (use_aux) {
          lav_samp_mi_patterns(y = y_g, mp = mp_g)
        } else {
          lavsamplestats@missing[[g]]
        },
        wt = lavdata@weights[[g]],
        cluster_idx = cluster_idx,
        mu = mu_g, sigma_1 = sigma_g,
        x_idx = lavsamplestats@x.idx[[g]],
        information = lavoptions$information[1]
      )
    }

    # auxiliary variables: reduce the augmented Omega to the model-variable
    # moments (sub-matrix of the inverse information)
    if (use_aux) {
      gamma[[g]] <- omega_g[midx, midx, drop = FALSE]
    } else {
      gamma[[g]] <- omega_g
    }

    # compute
    t_dvgvd <- t_dvgvd + fg * fg / fg1 * crossprod(wd, gamma[[g]] %*% wd)
  } # g

  nvar_cov <- (e_inv %*% t_dvgvd %*% e_inv)

  # to be reused by lavaanTest
  attr(nvar_cov, "Delta") <- delta
  attr(nvar_cov, "Gamma") <- gamma
  if ((lavoptions$information[1] == lavoptions$information[2]) &&
    (lavoptions$h1.information[1] == lavoptions$h1.information[2]) &&
    (lavoptions$information[2] == "expected" ||
      lavoptions$observed.information[1] ==
        lavoptions$observed.information[2])) {
    # only when same type of information is used # new in 0.6-6
    attr(nvar_cov, "E.inv") <- e_inv
    attr(nvar_cov, "WLS.V") <- wls_v
  }

  nvar_cov
}



lav_model_vcov <- function(lavmodel = NULL,
                           lavsamplestats = NULL,
                           lavoptions = NULL,
                           lavdata = NULL,
                           lavpartable = NULL,
                           lavcache = NULL,
                           lavimplied = NULL,
                           lavh1 = NULL,
                           use_ginv = FALSE) {
  likelihood <- lavoptions$likelihood
#  information <- lavoptions$information[1] # first one is for vcov
  se <- lavoptions$se
#  mimic <- lavoptions$mimic

  # special cases
  if (se == "none" || se == "external" || se == "twostep") {
    return(matrix(0, 0, 0))
  }

  if (se == "standard") {
    nvar_cov <- lav_model_info(
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
      use_ginv = use_ginv
    )
  } else if (se == "first.order") {
    nvar_cov <-
      lav_model_info_firstorder(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        extra = TRUE,
        check_pd = FALSE,
        augmented = TRUE,
        inverted = TRUE,
        use_ginv = use_ginv
      )
  } else if (se %in% c("robust.sem", "robust.sem.nt", "robust.cluster.sem")) {
    nvar_cov <-
      lav_model_nvcov_robust_sem(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavcache = lavcache,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        use_ginv = use_ginv
      )
  } else if (se == "robust.huber.white" || se == "robust.cluster") {
    nvar_cov <-
      lav_model_nvcov_robust_sandwich(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavcache = lavcache,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        use_ginv = use_ginv
      )
  } else if (se %in% c("two.stage", "robust.two.stage")) {
    nvar_cov <-
      lav_model_nvcov_two_stage(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavoptions = lavoptions,
        lavdata = lavdata,
        lavimplied = lavimplied,
        lavh1 = lavh1,
        use_ginv = use_ginv
      )
  } else if (se == "bootstrap") {
    nvar_cov <- try(
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

  if (!inherits(nvar_cov, "try-error")) {
    # denominator!
    if (lavmodel@estimator %in% c("ML", "PML", "FML") &&
      likelihood == "normal") {
      if (lavdata@nlevels == 1L) {
        n <- lavsamplestats@ntotal
        # new in 0.6-9 (to mimic method="lm" in effectLite)
        # special case: univariate regression in each group
        if (lavoptions$mimic == "lm" &&
          all(lavmodel@modprop$uvreg)) {
          n <- sum(unlist(lavsamplestats@nobs) -
            (unlist(lavmodel@modprop$nexo) + 1L))
          # always adding the intercept (for now)
        }
      } else {
        # total number of clusters (over groups)
        n <- 0
        for (g in 1:lavsamplestats@ngroups) {
          n <- n + lavdata@Lp[[g]]$nclusters[[2]]
        }
      }
    } else {
      n <- lavsamplestats@ntotal - lavsamplestats@ngroups
    }

    var_cov <- 1 / n * nvar_cov

    # check if VarCov is pd -- new in 0.6-2
    # mostly important if we have (in)equality constraints (MASS::ginv!)
    if (lavmodel@ceq.simple.only) {
      # do nothing
    } else if (!is.null(lavoptions$check.vcov) && lavoptions$check.vcov) {
      eigvals <- eigen(var_cov,
        symmetric = TRUE,
        only.values = TRUE
      )$values
      # correct for (in)equality constraints: the effective number of
      # imposed constraints is the RANK of the stacked equality + active
      # inequality rows (a plain row count over-corrects when rows are
      # redundant, within or across the two sets)
      neiq <- 0L
      if (nrow(lavmodel@con.jac) > 0L) {
        ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
        cin_idx <- attr(lavmodel@con.jac, "cin.idx")
        ina_idx <- attr(lavmodel@con.jac, "inactive.idx")
        act_idx <- c(ceq_idx, setdiff(cin_idx, ina_idx)) # only active cin
        if (length(act_idx) > 0L) {
          neiq <- qr(lavmodel@con.jac[act_idx, , drop = FALSE])$rank
        }
        if (neiq > 0L) {
          eigvals <- rev(eigvals)[-seq_len(neiq)]
        }
      }
      min_val <- min(eigvals)
      # if(any(eigvals < -1 * sqrt(.Machine$double.eps)) &&
      if (min_val < .Machine$double.eps^(3 / 4)) {
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

        if (min_val > 0) {
          lav_msg_warn(
            gettextf("The variance-covariance matrix of the estimated
                    parameters (vcov) does not appear to be positive
                    definite! The smallest eigenvalue (= %e) is close
                     to zero. This may be a symptom that the model is
                     not identified.", min(min_val)))
        } else {
          lav_msg_warn(
            gettextf("The variance-covariance matrix of the estimated parameters
                  (vcov) does not appear to be positive definite! The smallest
                  eigenvalue (= %e) is smaller than zero. This may be a
                  symptom that the model is not identified.",
                  min(min_val)))
        }
      }
    }
  } else {
    lav_msg_warn(
     gettext("Could not compute standard errors! The information matrix
             could not be inverted. This may be a symptom that the model
             is not identified.")
    )
    var_cov <- NULL
  } # could not invert

  var_cov
}

# Generate Monte Carlo draws for the free parameters from a multivariate
# normal distribution with mean = coef_hat (the point estimate) and
# covariance = VCOV. Used for Preacher & Selig (2012)-style Monte Carlo
# confidence intervals for defined parameters.
#
# Returns an R x n_free matrix. The number of draws and (optional) seed
# are taken from lavoptions$monte.carlo.
lav_model_vcov_mc <- function(lavmodel = NULL, vcov = NULL,
                              lavoptions = NULL) {
  if (is.null(vcov) || inherits(vcov, "try-error")) {
    return(NULL)
  }
  r <- 20000L
  seed <- NULL
  if (!is.null(lavoptions) && !is.null(lavoptions$monte.carlo)) {
    if (!is.null(lavoptions$monte.carlo$R)) {
      r <- as.integer(lavoptions$monte.carlo$R)
    }
    if (!is.null(lavoptions$monte.carlo$seed)) {
      seed <- lavoptions$monte.carlo$seed
    }
  }
  stopifnot(r > 0L)

  coef_hat <- lav_model_get_parameters(lavmodel = lavmodel, type = "free")
  # ensure we sample in the dimension of vcov: typically n_free, but
  # with simple equality constraints (ceq.simple.only) vcov is in unco
  # space; in that case lift coef_hat via t(K)
  if (lavmodel@ceq.simple.only && nrow(vcov) == nrow(lavmodel@ceq.simple.K)) {
    coef_hat <- drop(lavmodel@ceq.simple.K %*% coef_hat)
  }

  # save and restore RNG state so mc sampling does not perturb the
  # user's global seed
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    saved_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", saved_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit({
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # symmetrize vcov for safety
  vcov_sym <- 0.5 * (vcov + t(vcov))
  mc_coef <- lav_mvrnorm(n = r, mu = coef_hat, sigma_1 = vcov_sym,
                         check_symmetry = FALSE)
  if (!is.matrix(mc_coef)) {
    mc_coef <- matrix(mc_coef, nrow = r)
  }
  mc_coef
}

# Is the Monte Carlo method requested for defined parameters?
lav_model_vcov_se_mc_active <- function(lavoptions) {
  isTRUE(!is.null(lavoptions) &&
    identical(lavoptions$se.def, "monte.carlo"))
}

# Reduce a vcov matrix from the 'unco' space (one row/column per
# non-collapsed free parameter, nx.unco) to the 'compact' space (one
# row/column per free parameter, nx.free) for ceq.simple.only models. This is
# needed whenever the vcov must be combined with a jacobian/hessian that is
# expressed in the compact (nx.free) space (e.g. defined parameters, Wald
# tests). vcov is returned unchanged if it is not in the unco space.
lav_model_vcov_unco_to_free <- function(lavmodel, vcov) {
  if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L &&
      nrow(vcov) == nrow(lavmodel@ceq.simple.K)) {
    rep_idx <- apply(lavmodel@ceq.simple.K, 2L,
                     function(col) which(col != 0)[1L])
    vcov <- vcov[rep_idx, rep_idx, drop = FALSE]
  }
  vcov
}

lav_model_vcov_se <- function(lavmodel, lavpartable, vcov = NULL,
                              boot = NULL, mc = NULL,
                              lavoptions = NULL, ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  # 0. special case
  if (is.null(vcov)) {
    se <- rep(NA_real_, lavmodel@nx.user)
    se[lavpartable$free == 0L] <- 0.0
    return(se)
  }

  # 1. free parameters only
  x_var <- diag(vcov)
  # check for negative values (what to do: NA or 0.0?)
  x_var[x_var < 0] <- as.numeric(NA)
  x_se <- sqrt(x_var)
  if (lavmodel@ceq.simple.only) {
    glist <- lav_model_x2glist(
      lavmodel = lavmodel, x = x_se,
      type = "unco"
    )
  } else {
    glist <- lav_model_x2glist(
      lavmodel = lavmodel, x = x_se,
      type = "free"
    )
  }

  # se for full parameter table, but with 0.0 entries for def/ceq/cin
  # elements
  se <- lav_model_get_parameters(
    lavmodel = lavmodel, glist = glist,
    type = "user", extra = FALSE
  )


  # 2. fixed parameters -> se = 0.0
  se[which(lavpartable$free == 0L)] <- 0.0


  # 3. defined parameters:
  def_idx <- which(lavpartable$op == ":=")
  if (length(def_idx) > 0L) {
    # if Monte Carlo samples are not supplied but requested, draw them
    if (is.null(mc) && is.null(boot) &&
        lav_model_vcov_se_mc_active(lavoptions)) {
      mc <- lav_model_vcov_mc(lavmodel = lavmodel, vcov = vcov,
                              lavoptions = lavoptions)
    }
    if (!is.null(mc)) {
      mc_def <- apply(mc, 1L, lavmodel@def.function)
      if (length(def_idx) == 1L) {
        mc_def <- as.matrix(mc_def)
      } else {
        mc_def <- t(mc_def)
      }
      # def.function maps invalid evaluations (NaN) to +Inf;
      # treat those draws as missing so that cov() ignores them
      mc_def[!is.finite(mc_def)] <- as.numeric(NA)
      def_cov <- cov(mc_def, use = "pairwise.complete.obs")
      diag_def_cov <- diag(def_cov)
      diag_def_cov[diag_def_cov < 0] <- as.numeric(NA)
      se[def_idx] <- sqrt(diag_def_cov)
      return(se)
    }
    if (!is.null(boot)) {
      # we must remove the NA rows (and hope we have something left)
      error_idx <- attr(boot, "error.idx")
      boot <- boot
      if (length(error_idx) > 0L) {
        boot <- boot[-error_idx, , drop = FALSE] # drops attributes
      }

      boot_def <- apply(boot, 1L, lavmodel@def.function)
      if (length(def_idx) == 1L) {
        boot_def <- as.matrix(boot_def)
      } else {
        boot_def <- t(boot_def)
      }
      # new in 0.6-20: check for outliers, big difference between sd() and mad()
      # see github issue 347
      sd_mad_ratio <- (apply(boot_def, 2,  sd, na.rm = TRUE) /
                       apply(boot_def, 2, mad, na.rm = TRUE))
      crit_ratio <- 5
      if (any(sd_mad_ratio > crit_ratio)) {
        names_1 <- colnames(boot_def)
        def_w_outliers <- paste(names_1[sd_mad_ratio > crit_ratio],
                                collapse = " ")
        lav_msg_warn(gettextf(
          "The following bootstrapped defined parameters have a high (>5)
          ratio of standard deviation to median absolute deviation: %s.
          P-values and confidence intervals may not match.", def_w_outliers))
      }
      nboot <- nrow(boot_def)
      def_cov <- cov(boot_def) * (nboot - 1) / nboot
    } else {
      # regular delta method
      x <- lav_model_get_parameters(lavmodel = lavmodel, type = "free")
      jac <- try(lav_func_jacobian_complex(func = lavmodel@def.function, x = x),
        silent = TRUE
      )
      use_complex <- !inherits(jac, "try-error")
      if (!use_complex) { # eg. pnorm()
        jac <- lav_func_jacobian_simple(func = lavmodel@def.function, x = x)
      }

      # second-order delta method: compute Hessian for each defined param
      hess_list <- NULL
      if (!is.null(lavoptions) &&
        isTRUE(lavoptions$se.delta.second.order)) {
        if (use_complex) {
          hess_list <- try(lav_func_hessian_complex(
            func = lavmodel@def.function, x = x), silent = TRUE)
          if (inherits(hess_list, "try-error")) {
            hess_list <- try(lav_func_hessian_simple(
              func = lavmodel@def.function, x = x), silent = TRUE)
          }
        } else {
          hess_list <- try(lav_func_hessian_simple(
            func = lavmodel@def.function, x = x), silent = TRUE)
        }
        if (inherits(hess_list, "try-error")) {
          lav_msg_warn(gettext(
            "Could not compute the Hessian of the defined parameters;
             falling back to the first-order delta method."))
          hess_list <- NULL
        }
      }

      # jac (and hess_list) are in the compact (nx.free) space; when
      # ceq.simple.only, vcov is in the larger 'unco' space, so reduce vcov to
      # the compact space. NB: previously jac was *expanded* to the unco space
      # via ceq.simple.K, which multiply-counts the shared parameters and
      # inflates the SE of defined parameters that involve them.
      vcov_def <- lav_model_vcov_unco_to_free(lavmodel, vcov)
      def_cov <- jac %*% vcov_def %*% t(jac)

      # add second-order delta method correction:
      # Cov(g_i, g_j) ~ grad_i' V grad_j + 0.5 tr(H_i V H_j V)
      if (!is.null(hess_list)) {
        ndef_1 <- length(hess_list)
        correction <- matrix(0, ndef_1, ndef_1)
        hv_list <- lapply(hess_list, function(h) h %*% vcov_def)
        for (i in seq_len(ndef_1)) {
          for (j in i:ndef_1) {
            # tr(A B) where A = H_i V, B = H_j V
            #   = sum_{a,b} A_{a,b} B_{b,a}
            #   = sum element-wise of A and t(B)
            correction[i, j] <- 0.5 * sum(hv_list[[i]] * t(hv_list[[j]]))
            if (i != j) {
              correction[j, i] <- correction[i, j]
            }
          }
        }
        def_cov <- def_cov + correction
      }
    }
    # check for negative se's
    diag_def_cov <- diag(def_cov)
    diag_def_cov[diag_def_cov < 0] <- as.numeric(NA)
    se[def_idx] <- sqrt(diag_def_cov)
  }

  se
}
