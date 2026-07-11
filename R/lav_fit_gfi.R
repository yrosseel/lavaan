# functions related to GFI and other 'absolute' fit indices

# lower-level functions:
# - lav_fit_gfi_lisrel    (the 'old' LISREL GFI)
# - lav_fit_agfi_lisrel
# - lav_fit_pgfi
# - lav_fit_gfi           (the 'new' GFI: Maydeu-Olivares et al. 2024)
# - lav_fit_gfi_ci

# higher-level functions:
# - lav_fit_gfi_lavobject

# Y.R. 21 July 2022
# Y.R. 24 June 2026: rename old gfi/agfi to gfi_lisrel/agfi_lisrel and add the
#                    'new' GFI (with confidence interval) of Maydeu-Olivares et
#                    al. (2024)

# ----------------------------------------------------------------------------
# the 'old' (LISREL) GFI and friends
# ----------------------------------------------------------------------------

# original formulas were given in Joreskog and Sorbom (1984) user's guide
# for LISREL VI (one for ML, and another for ULS)

# here we use the more 'general' formulas
# (generalized to allow for meanstructures etc)

# References:

# Mulaik, S. A., James, L. R., Van Alstine, J., Bennett, N., Lind, S., &
# Stilwell, C. D. (1989). Evaluation of goodness-of-fit indices for structural
# equation models. Psychological bulletin, 105(3), 430.

# Tanaka, J. S., & Huba, G. J. (1985). A fit index for covariance structure
# models under arbitrary GLS estimation. British Journal of Mathematical and
# Statistical Psychology, 38,197-201.

lav_fit_gfi_lisrel <- function(wls_obs = NULL, wls_est = NULL, wls_v = NULL,
                               nobs_1 = NULL) {
  # number of groups
  g_1 <- length(wls_obs)

  # compute gfi per group
  gfi_group <- numeric(g_1)
  for (g in 1:g_1) {
    wls_obs_1 <- wls_obs[[g]]
    wls_est_1 <- wls_est[[g]]
    wls_v_1 <- wls_v[[g]]

    if (is.null(wls_v_1)) {
      gfi_group[g] <- as.numeric(NA)
    } else {
      wls_diff <- wls_obs_1 - wls_est_1
      if (is.matrix(wls_v_1)) {
        # full weight matrix
        t1 <- crossprod(wls_diff, wls_v_1) %*% wls_diff
        t2 <- crossprod(wls_obs_1, wls_v_1) %*% wls_obs_1
      } else {
        # diagonal weight matrix
        t1 <- as.numeric(crossprod(wls_diff^2, wls_v_1))
        t2 <- as.numeric(crossprod(wls_obs_1^2, wls_v_1))
      }
      gfi_group[g] <- 1 - t1 / t2
    }
  }

  if (g_1 > 1) {
    ## CHECKME: get the scaling right
    nobs_1 <- unlist(nobs_1)
    gfi <- as.numeric((nobs_1 %*% gfi_group) / sum(nobs_1))
  } else {
    gfi <- gfi_group[1L]
  }

  gfi
}

# 'adjusted' GFI (adjusted for degrees of freedom)
lav_fit_agfi_lisrel <- function(gfi = NULL, nel = NULL, df = NULL) {
  if (!is.finite(gfi) || !is.finite(nel) || !is.finite(df)) {
    agfi <- as.numeric(NA)
  } else if (df > 0) {
    agfi <- 1 - (nel / df) * (1 - gfi)
  } else {
    agfi <- 1
  }

  agfi
}

# PGFI: parsimony goodness-of-fit index

# Mulaik, S. A., James, L. R., Van Alstine, J., Bennett, N., Lind, S., &
# Stilwell, C. D. (1989). Evaluation of goodness-of-fit indices for structural
# equation models. Psychological bulletin, 105(3), 430.

# LISREL formula (Simplis book 2002, p. 126)
lav_fit_pgfi <- function(gfi = NULL, nel = NULL, df = NULL) {
  if (!is.finite(gfi) || !is.finite(nel) || !is.finite(df)) {
    pgfi <- as.numeric(NA)
  } else if (nel == 0) {
    pgfi <- as.numeric(NA)
  } else {
    pgfi <- (df / nel) * gfi
  }

  pgfi
}

# ----------------------------------------------------------------------------
# the 'new' GFI (with confidence interval)
# ----------------------------------------------------------------------------

# Maydeu-Olivares, A., Ximenez, C., & Revuelta, J. (2024). Percentage of
# variance accounted for in structural equation models: The rediscovery of the
# goodness of fit index. Psychological Methods.

# The GFI is interpreted as the percentage of variance (of the sum of squared
# covariances) accounted for by the model; it is the SEM analogue of R^2 in
# regression (without an intercept), using the saturated model as baseline.

# point estimate (Maydeu-Olivares et al., 2024, Eq. 20):
#
#   GFI = min( N*p / (N*p + 2*c_hat*T - 2*c_hat*df), 1 )
#
# where T is the test statistic, p the number of observed variables, N the
# (total) sample size, df the degrees of freedom and c_hat a scaling constant.
#
# Maydeu-Olivares et al. (2024) use the reweighted least squares (RLS)
# statistic T_RLS (== 'browne.residual.nt.model' in lavaan), arguing it offers
# a more reliable assessment of exact fit than the likelihood ratio (LR)
# statistic. We compute the GFI from both the LR statistic (gfi_lrt) and the
# RLS statistic (gfi_rls); the (printed) 'gfi' is the RLS version.
#
# Following the implementation of the (robust) RMSEA in lavaan, we use the
# (unscaled) statistic X2 as 'c_hat*T' (since c_hat*T_M equals the unadjusted
# RLS/LR statistic) and apply c_hat only to the degrees of freedom:
#
#   GFI = min( N*p / (N*p + 2*X2 - 2*c_hat*df), 1 )
#
# - 'standard' (normal-theory, unbiased) version: c_hat = 1
# - 'robust' (to nonnormality) version: c_hat = scaling factor
#
# Note: the scaling constant c (and the shift b for scaled.shifted) is the same
# for T_LR and T_RLS (Maydeu-Olivares et al., 2024, Eq. 19), so the robust RLS
# version reuses the scaling factor(s) obtained for the LR statistic.
#
# Note: this matches the relationship GFI = p / (p + 2*df*RMSEA^2) (their
# Eq. 26), with the corresponding (standard or robust) RMSEA value.
#
# always using N (if a user needs N-1, just replace N by N-1)
lav_fit_gfi <- function(x2 = NULL, df = NULL, n = NULL, p = NULL,
                        c_hat = 1.0) {
  if (length(x2) == 0L || !is.finite(x2) || !is.finite(df) ||
    !is.finite(n) || !is.finite(p)) {
    return(as.numeric(NA))
  }

  np <- n * p
  denom <- np + 2 * x2 - 2 * c_hat * df
  if (!is.finite(denom) || denom <= 0) {
    return(as.numeric(NA))
  }

  gfi <- np / denom

  # see Eq. 20: GFI can not be larger than 1 (eg when the model fits 'too well')
  if (gfi > 1) {
    gfi <- 1
  }

  gfi
}

# confidence interval for the 'new' GFI (Maydeu-Olivares et al., 2024, Eq. 23)
#
# the confidence interval is obtained from the noncentrality parameter of the
# (noncentral) chi-square distribution (exactly as for the RMSEA), and then
# transformed to the GFI metric:
#
#   GFI_bound = min( N*p / (N*p + 2*c_hat*lambda), 1 )
#
# note: the *lower* noncentrality (lambda_l) yields the *upper* GFI bound, and
# the *upper* noncentrality (lambda_u) yields the *lower* GFI bound.
#
# note: for the 'robust' version, x2 should be the scaled test statistic
lav_fit_gfi_ci <- function(x2 = NULL, df = NULL, n = NULL, p = NULL,
                           c_hat = 1, level = 0.90) {
  if (length(x2) == 0L || !is.finite(x2) || !is.finite(df) ||
    !is.finite(n) || !is.finite(p)) {
    return(list(
      gfi.ci.lower = as.numeric(NA),
      gfi.ci.upper = as.numeric(NA)
    ))
  }

  if (!is.finite(level) || level < 0 || level > 1.0) {
    lav_msg_warn(gettextf(
      "invalid level value [%s] set to default 0.90.", level
    ))
    level <- 0.90
  }

  np <- n * p

  upper_perc <- (1 - (1 - level) / 2)
  lower_perc <- (1 - level) / 2

  # internal functions
  lower_lambda <- function(lambda) {
    (pchisq(x2, df = df, ncp = lambda) - upper_perc)
  }
  upper_lambda <- function(lambda) {
    (pchisq(x2, df = df, ncp = lambda) - lower_perc)
  }

  # transform a noncentrality parameter to the GFI metric
  gfi_from_lambda <- function(lambda) {
    if (!is.finite(lambda)) {
      return(as.numeric(NA))
    }
    val <- np / (np + 2 * c_hat * lambda)
    if (val > 1) {
      val <- 1
    }
    val
  }

  # lower noncentrality -> upper GFI bound
  if (df < 1 || lower_lambda(0) < 0.0) {
    gfi_ci_upper <- 1
  } else {
    lambda_l <- try(uniroot(f = lower_lambda, lower = 0, upper = x2)$root,
      silent = TRUE
    )
    if (inherits(lambda_l, "try-error")) {
      lambda_l <- as.numeric(NA)
    }
    gfi_ci_upper <- gfi_from_lambda(lambda_l)
  }

  # upper noncentrality -> lower GFI bound
  n_gfi <- max(n, x2 * 4)
  if (df < 1 || upper_lambda(n_gfi) > 0 || upper_lambda(0) < 0) {
    gfi_ci_lower <- 1
  } else {
    lambda_u <- try(
      uniroot(
        f = upper_lambda, lower = 0,
        upper = n_gfi
      )$root,
      silent = TRUE
    )
    if (inherits(lambda_u, "try-error")) {
      lambda_u <- as.numeric(NA)
    }
    gfi_ci_lower <- gfi_from_lambda(lambda_u)
  }

  list(
    gfi.ci.lower = gfi_ci_lower,
    gfi.ci.upper = gfi_ci_upper
  )
}

# compute a single GFI 'family' (point estimate + confidence interval, and the
# robust counterparts) for a given (unadjusted) test statistic. This is called
# once for the LR statistic and once for the RLS statistic.
lav_fit_gfi_family <- function(x2 = NULL, df = NULL, n = NULL, p = NULL,
                               ci_level = 0.90, robust_flag = FALSE,
                               xx3 = NULL, df3 = NULL, c_hat = NULL,
                               c_hat3 = NULL, xx3_scaled = NULL,
                               n3 = NULL) {
  # n3: the sample size to be used for the *robust* quantities (which may
  # be N rather than N-1 in the categorical case; see lav_fit_rmsea.R)
  if (is.null(n3)) {
    n3 <- n
  }
  out <- list(
    est = as.numeric(NA), ci.lower = as.numeric(NA), ci.upper = as.numeric(NA),
    robust = as.numeric(NA), ci.lower.robust = as.numeric(NA),
    ci.upper.robust = as.numeric(NA)
  )

  # point estimate (normal-theory, unbiased)
  out$est <- lav_fit_gfi(x2 = x2, df = df, n = n, p = p)

  # confidence interval
  ci <- lav_fit_gfi_ci(x2 = x2, df = df, n = n, p = p, level = ci_level)
  out$ci.lower <- ci$gfi.ci.lower
  out$ci.upper <- ci$gfi.ci.upper

  if (robust_flag) {
    out$robust <- lav_fit_gfi(
      x2 = xx3, df = df3, n = n3, p = p, c_hat = c_hat3
    )
    # note: input is the scaled test statistic!
    ci_robust <- lav_fit_gfi_ci(
      x2 = xx3_scaled, df = df3, n = n3, p = p,
      c_hat = c_hat, level = ci_level
    )
    out$ci.lower.robust <- ci_robust$gfi.ci.lower
    out$ci.upper.robust <- ci_robust$gfi.ci.upper
  }

  out
}

# ----------------------------------------------------------------------------
# higher-level function
# ----------------------------------------------------------------------------

lav_fit_gfi_lavobject <- function(lavobject = NULL, fit_measures = "gfi",
                                  standard_test = "standard",
                                  scaled_test = "none",
                                  ci_level = 0.90,
                                  robust = TRUE,
                                  cat_nonpd = "na") {
  # check lavobject
  stopifnot(inherits(lavobject, "lavaan"))

  # check for categorical
  categorical_flag <- lavobject@Model@categorical

  # tests
  test <- lavobject@test
  test_names <- sapply(lavobject@test, "[[", "test")
  if (test_names[1] == "none" || standard_test == "none") {
    return(list())
  }
  test_idx <- which(test_names == standard_test)[1]
  if (length(test_idx) == 0L || is.na(test_idx)) {
    return(list())
  }

  scaled_flag <- FALSE
  if (!scaled_test %in% c("none", "standard", "default")) {
    scaled_idx <- which(test_names == scaled_test)
    if (length(scaled_idx) > 0L) {
      scaled_idx <- scaled_idx[1] # only the first one
      scaled_flag <- TRUE
    }
  }

  # robust?
  robust_flag <- FALSE
  if (robust && scaled_flag &&
    scaled_test %in% c(
      "satorra.bentler", "yuan.bentler.mplus",
      "yuan.bentler", "yuan.chan", "scaled.shifted",
      "scaled.shifted.corrected"
    )) {
    robust_flag <- TRUE
  }

  # FIML? (the FIML-corrected robust indices are only available for
  # single-level models; for multilevel data, we fall back to the regular
  # scaled/robust computation above, which handles missing data too)
  fiml_flag <- FALSE
  if (robust && lavobject@Options$missing %in% c("ml", "ml.x") &&
    lavobject@Data@nlevels == 1L) {
    fiml_flag <- robust_flag <- TRUE
    # check if we can compute corrected values
    if (scaled_flag) {
      version <- "V3"
    } else {
      version <- "V6"
    }
    fiml <- try(
      lav_fit_fiml_corrected(lavobject,
        baseline_model = NULL,
        version = version
      ),
      silent = TRUE
    )
    if (inherits(fiml, "try-error")) {
      lav_msg_warn(gettext("computation of robust GFI failed."))
      fiml <- list(
        XX3 = as.numeric(NA), df3 = as.numeric(NA),
        c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA)
      )
    } else if (anyNA(c(fiml$XX3, fiml$df3, fiml$c.hat3, fiml$XX3.scaled))) {
      lav_msg_warn(gettext(
        "computation of robust GFI resulted in NA values."
      ))
    }
  }

  # supported fit measures in this function
  # 1) the 'old' (LISREL) GFI and friends
  fit_gfi_lisrel <- c("gfi_lisrel", "agfi_lisrel", "pgfi")
  # 2) the 'new' GFI (with confidence interval); 'gfi' is the RLS version
  fit_gfi <- c("gfi", "gfi.ci.lower", "gfi.ci.upper", "gfi.ci.level")
  if (robust_flag) {
    fit_gfi <- c(
      fit_gfi, "gfi.robust",
      "gfi.ci.lower.robust", "gfi.ci.upper.robust"
    )
  }
  # 3) the explicit LR-based (gfi_lrt) and RLS-based (gfi_rls) versions
  fit_gfi_extra <- c(
    "gfi_lrt", "gfi_lrt.ci.lower", "gfi_lrt.ci.upper",
    "gfi_rls", "gfi_rls.ci.lower", "gfi_rls.ci.upper"
  )
  if (robust_flag) {
    fit_gfi_extra <- c(
      fit_gfi_extra,
      "gfi_lrt.robust", "gfi_lrt.ci.lower.robust", "gfi_lrt.ci.upper.robust",
      "gfi_rls.robust", "gfi_rls.ci.lower.robust", "gfi_rls.ci.upper.robust"
    )
  }
  fit_gfi_all <- c(fit_gfi_lisrel, fit_gfi, fit_gfi_extra)

  # which one do we need?
  if (missing(fit_measures)) {
    # default set
    fit_measures <- fit_gfi_all
  } else {
    # remove any not-GFI related index from fit.measures
    rm_idx <- which(!fit_measures %in% fit_gfi_all)
    if (length(rm_idx) > 0L) {
      fit_measures <- fit_measures[-rm_idx]
    }
    if (length(fit_measures) == 0L) {
      return(list())
    }
  }

  # container
  indices <- list()

  # ---------- the 'old' (LISREL) GFI and friends ----------
  if (any(fit_gfi_lisrel %in% fit_measures)) {
    # extract ingredients
    wls_obs <- lav_inspect_wls_obs(lavobject)
    wls_est <- lav_inspect_wls_est(lavobject)
    wls_v <- lav_inspect_wls_v(lavobject)
    nobs_1 <- lavobject@SampleStats@nobs

    # compute (LISREL) GFI
    gfi_lisrel <- lav_fit_gfi_lisrel(
      wls_obs = wls_obs, wls_est = wls_est,
      wls_v = wls_v, nobs_1 = nobs_1
    )

    # total number of modeled sample stats
    nel <- length(unlist(wls_obs))

    # degrees of freedom
    df <- lavobject@test[[1]]$df

    indices["gfi_lisrel"] <- gfi_lisrel
    indices["agfi_lisrel"] <-
      lav_fit_agfi_lisrel(gfi = gfi_lisrel, nel = nel, df = df)
    indices["pgfi"] <- lav_fit_pgfi(gfi = gfi_lisrel, nel = nel, df = df)
  }

  # ---------- the 'new' GFI (with confidence interval) ----------
  if (any(c(fit_gfi, fit_gfi_extra) %in% fit_measures)) {
    # basic ingredients
    x2 <- test[[test_idx]]$stat # LR (standard) chi-square statistic
    df <- test[[test_idx]]$df
    n <- lav_inspect_ntotal(object = lavobject) # N vs N-1
    # number of observed variables (per group)
    p <- lavobject@Model@nvar[[1]]

    # the RLS (reweighted least squares) statistic == 'browne.residual.nt.model'
    # (Maydeu-Olivares et al., 2024, recommend this over the LR statistic). Only
    # needed if an RLS-based measure (gfi, gfi_rls) is requested. There is no
    # separate RLS statistic in the categorical/FIML case (we reuse the LR one).
    rls_names <- c(fit_gfi, fit_gfi_extra[grepl("^gfi_rls", fit_gfi_extra)])
    x2_rls <- x2 # fallback
    if (any(rls_names %in% fit_measures) &&
      !categorical_flag && !fiml_flag) {
      rls <- try(lavTest(lavobject, test = "browne.residual.nt.model"),
        silent = TRUE
      )
      if (!inherits(rls, "try-error") &&
        !is.null(rls[["browne.residual.nt.model"]]) &&
        is.finite(rls[["browne.residual.nt.model"]]$stat)) {
        x2_rls <- rls[["browne.residual.nt.model"]]$stat
      }
    }

    # robust ingredients
    xx3 <- xx3_scaled <- df3 <- c_hat <- c_hat3 <- as.numeric(NA)
    xx3_rls <- xx3_scaled_rls <- as.numeric(NA)
    n3 <- n
    if (robust_flag) {
      if (categorical_flag) {
        out <- try(lav_fit_catml_dwls(lavobject, nonpd = cat_nonpd),
          silent = TRUE
        )
        if (inherits(out, "try-error")) {
          xx3 <- df3 <- c_hat3 <- xx3_scaled <- as.numeric(NA)
          c_hat <- as.numeric(NA)
        } else {
          xx3 <- out$XX3
          df3 <- out$df3
          c_hat3 <- c_hat <- out$c.hat3
          xx3_scaled <- out$XX3.scaled
        }
        # no separate RLS version for categorical data
        xx3_rls <- xx3
        xx3_scaled_rls <- xx3_scaled
        # the robust quantities are ML-based: use N, not N-1
        # (see lav_fit_rmsea.R)
        n3 <- lavobject@SampleStats@ntotal
      } else if (fiml_flag) {
        xx3 <- fiml$XX3
        df3 <- fiml$df3
        c_hat3 <- c_hat <- fiml$c.hat3
        xx3_scaled <- fiml$XX3.scaled
        # no separate RLS version for FIML
        xx3_rls <- xx3
        xx3_scaled_rls <- xx3_scaled
      } else {
        xx3 <- x2
        df3 <- df
        c_hat <- test[[scaled_idx]]$scaling.factor
        shift_par <- 0
        if (scaled_test %in% c(
          "scaled.shifted", "scaled.shifted.corrected"
        )) {
          # compute c.hat from a and b
          a <- test[[scaled_idx]]$scaling.factor
          b <- test[[scaled_idx]]$shift.parameter
          shift_par <- b
          c_hat3 <- a * (df - b) / df
        } else {
          c_hat3 <- c_hat
        }
        xx3_scaled <- test[[scaled_idx]]$stat
        # RLS analogues: the scaling constant c (and the shift b) are identical
        # for T_LR and T_RLS (Maydeu-Olivares et al., 2024, Eq. 19), so we only
        # need the unadjusted RLS statistic and reuse the same scaling factor(s)
        xx3_rls <- x2_rls
        xx3_scaled_rls <- x2_rls / c_hat + shift_par
      }
    }

    # GFI based on the LR statistic
    gfi_lrt <- lav_fit_gfi_family(
      x2 = x2, df = df, n = n, p = p, ci_level = ci_level,
      robust_flag = robust_flag, xx3 = xx3, df3 = df3,
      c_hat = c_hat, c_hat3 = c_hat3, xx3_scaled = xx3_scaled,
      n3 = n3
    )
    # GFI based on the RLS statistic (the recommended version)
    gfi_rls <- lav_fit_gfi_family(
      x2 = x2_rls, df = df, n = n, p = p, ci_level = ci_level,
      robust_flag = robust_flag, xx3 = xx3_rls, df3 = df3,
      c_hat = c_hat, c_hat3 = c_hat3, xx3_scaled = xx3_scaled_rls,
      n3 = n3
    )

    # the default (printed) 'gfi' is the RLS version
    indices["gfi"] <- gfi_rls$est
    indices["gfi.ci.level"] <- ci_level
    indices["gfi.ci.lower"] <- gfi_rls$ci.lower
    indices["gfi.ci.upper"] <- gfi_rls$ci.upper

    # explicit LR-based and RLS-based versions
    indices["gfi_lrt"] <- gfi_lrt$est
    indices["gfi_lrt.ci.lower"] <- gfi_lrt$ci.lower
    indices["gfi_lrt.ci.upper"] <- gfi_lrt$ci.upper
    indices["gfi_rls"] <- gfi_rls$est
    indices["gfi_rls.ci.lower"] <- gfi_rls$ci.lower
    indices["gfi_rls.ci.upper"] <- gfi_rls$ci.upper

    if (robust_flag) {
      indices["gfi.robust"] <- gfi_rls$robust
      indices["gfi.ci.lower.robust"] <- gfi_rls$ci.lower.robust
      indices["gfi.ci.upper.robust"] <- gfi_rls$ci.upper.robust

      indices["gfi_lrt.robust"] <- gfi_lrt$robust
      indices["gfi_lrt.ci.lower.robust"] <- gfi_lrt$ci.lower.robust
      indices["gfi_lrt.ci.upper.robust"] <- gfi_lrt$ci.upper.robust
      indices["gfi_rls.robust"] <- gfi_rls$robust
      indices["gfi_rls.ci.lower.robust"] <- gfi_rls$ci.lower.robust
      indices["gfi_rls.ci.upper.robust"] <- gfi_rls$ci.upper.robust
    }
  }

  # return only those that were requested
  indices[fit_measures]
}
