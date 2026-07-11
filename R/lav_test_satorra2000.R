# machinery for the 'robust' versions of restricted tests (score and Wald):
# the scaled, (mean and variance) adjusted, and generalized (asymptotically
# distribution-free) test statistics
#
# reference:
#
# Satorra, A. (2000). Scaled and adjusted restricted tests in multi-sample
# analysis of moment structures. In Heijmans, R.D.H., Pollock, D.S.G. &
# Satorra, A. (Eds.), Innovations in multivariate statistical analysis:
# A festschrift for Heinz Neudecker (pp. 233-247). London: Kluwer.
#
# YR 11 July 2026: initial version

# which sandwich 'flavor' can we use to robustify a restricted test for
# this object?
#
# - if the se's are already sandwich-based, use that flavor
# - otherwise, if a scaled test statistic was requested, the ingredients
#   are available anyway: use Gamma (robust.sem) if we have it, or the
#   casewise scores (robust.huber.white) if we have raw data
# - otherwise, return NULL (only the standard test is available)
lav_test_robust_se_flavor <- function(object) {
  lavoptions <- object@Options

  sandwich_se <- c(
    "robust.sem", "robust.sem.nt", "robust.cluster.sem",
    "robust.huber.white", "robust.cluster"
  )
  if (lavoptions$se %in% sandwich_se) {
    return(lavoptions$se)
  }

  scaled_tests <- c(
    "satorra.bentler", "yuan.bentler", "yuan.bentler.mplus",
    "mean.var.adjusted", "scaled.shifted",
    "mean.var.adjusted.corrected", "scaled.shifted.corrected"
  )
  if (!any(lavoptions$test %in% scaled_tests)) {
    return(NULL)
  }

  gamma <- object@SampleStats@NACOV
  if (length(gamma) > 0L && !is.null(gamma[[1]])) {
    return("robust.sem")
  }
  if (object@Model@estimator %in% c("ML", "PML") &&
    length(object@Data@X) > 0L) {
    return("robust.huber.white")
  }

  NULL
}

# the (unit-scale) 'meat' of the sandwich covariance matrix, in the same
# parameter space (unco when ceq.simple.only) and on the same scale as
# lavTech(object, "information.*"):
#
#   B = sum_g f_g Delta_g' V_g Gamma_g V_g Delta_g   (robust.sem flavors)
#   B = first-order information (casewise scores)    (robust.huber.white)
#
# so that n * Avar(thetahat) = J %*% B %*% J, with J the (constrained)
# inverse of the unit information; returns NULL if the ingredients are
# not available for this fit
lav_test_meat <- function(object, se = NULL) {
  lavoptions <- object@Options
  lavmodel <- object@Model

  if (is.null(se)) {
    se <- lav_test_robust_se_flavor(object)
  }
  if (is.null(se)) {
    return(NULL)
  }

  # the h1 slot may be empty; downstream functions recompute it if NULL
  lavh1 <- object@h1
  if (length(lavh1) == 0L) {
    lavh1 <- NULL
  }

  if (se %in% c("robust.sem", "robust.sem.nt", "robust.cluster.sem")) {
    # Gamma needed
    gamma <- object@SampleStats@NACOV
    if (length(gamma) == 0L || is.null(gamma[[1]])) {
      return(NULL)
    }
    nvcov <- try(
      lav_model_nvcov_robust_sem(
        lavmodel = lavmodel,
        lavsamplestats = object@SampleStats,
        lavdata = object@Data,
        lavcache = object@Cache,
        lavimplied = object@implied,
        lavh1 = lavh1,
        lavoptions = lavoptions,
        attr_delta = FALSE,
        attr_t_dvgvd = TRUE
      ),
      silent = TRUE
    )
    if (inherits(nvcov, "try-error")) {
      return(NULL)
    }
    b_meat <- attr(nvcov, "tDVGVD")
  } else {
    # casewise (or clusterwise) scores
    if (!lavmodel@estimator %in% c("ML", "PML")) {
      return(NULL)
    }
    # allow information.meat/h1.information.meat overrides, as in
    # lav_model_nvcov_robust_sandwich()
    lavoptions2 <- lavoptions
    if (!is.null(lavoptions$information.meat)) {
      lavoptions2$information <- lavoptions$information.meat
    }
    if (!is.null(lavoptions$h1.information.meat)) {
      lavoptions2$h1.information <- lavoptions$h1.information.meat
    }
    b_meat <- try(
      lav_model_info_firstorder(
        lavmodel = lavmodel,
        lavsamplestats = object@SampleStats,
        lavdata = object@Data,
        lavimplied = object@implied,
        lavh1 = lavh1,
        lavcache = object@Cache,
        lavoptions = lavoptions2,
        check_pd = FALSE,
        extra = FALSE,
        augmented = FALSE,
        inverted = FALSE
      ),
      silent = TRUE
    )
    if (inherits(b_meat, "try-error")) {
      return(NULL)
    }
  }

  b_meat
}

# scaled, adjusted, and generalized ('robust') versions of a restricted
# test statistic (score or Wald), following Satorra (2000, section 17.3)
#
# arguments:
#   stat : the standard (normal-theory) test statistic
#   df   : its degrees of freedom (m = number of tested restrictions)
#   m1   : A J A'      (m x m, 'bread-only' middle matrix)
#   m2   : A J B J A'  (m x m, sandwich middle matrix)
#          (equivalently: n * A Avar(thetahat) A')
#   v, n : the generalized statistic is n * t(v) %*% solve(m2) %*% v
#          (score test: v = A J h with n the sample size; Wald test:
#           v = a(thetahat), with m1/m2 based on vcov and n = 1)
#
# the trace quantities follow from eq. (22)-(23):
#   tr(U Gamma) = tr[ (A J B J A') (A J A')^{-1} ] = tr(m2 %*% solve(m1))
lav_test_satorra2000 <- function(stat = NULL, df = NULL,
                                 m1 = NULL, m2 = NULL,
                                 v = NULL, n = 1) {
  out <- list(
    scaling.factor = as.numeric(NA),
    stat.scaled = as.numeric(NA),
    df.scaled = df,
    p.value.scaled = as.numeric(NA),
    stat.adjusted = as.numeric(NA),
    df.adjusted = as.numeric(NA),
    p.value.adjusted = as.numeric(NA),
    stat.robust = as.numeric(NA),
    df.robust = df,
    p.value.robust = as.numeric(NA)
  )

  m1_inv <- try(MASS::ginv(m1), silent = TRUE)
  if (!inherits(m1_inv, "try-error")) {
    ug <- m2 %*% m1_inv
    tr_ug <- sum(diag(ug))
    tr_ug2 <- sum(ug * t(ug))

    # scaled
    if (is.finite(tr_ug) && tr_ug > 0) {
      out$scaling.factor <- tr_ug / df
      out$stat.scaled <- stat / out$scaling.factor
      out$p.value.scaled <- 1 - pchisq(out$stat.scaled, df = df)

      # adjusted (for mean and variance), with fractional df
      if (is.finite(tr_ug2) && tr_ug2 > 0) {
        df_adj <- tr_ug * tr_ug / tr_ug2
        out$stat.adjusted <- stat * df_adj / tr_ug
        out$df.adjusted <- df_adj
        out$p.value.adjusted <- 1 - pchisq(out$stat.adjusted, df = df_adj)
      }
    }
  }

  # generalized ('robust')
  if (!is.null(v)) {
    m2_inv <- try(MASS::ginv(m2), silent = TRUE)
    if (!inherits(m2_inv, "try-error")) {
      out$stat.robust <- as.numeric(n * t(v) %*% m2_inv %*% v)
      out$p.value.robust <- 1 - pchisq(out$stat.robust, df = df)
    }
  }

  out
}
