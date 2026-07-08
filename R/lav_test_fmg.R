# FMG (Foldnes, Moss, Gronneberg 2024) test statistics
# Improved goodness-of-fit tests using eigenvalue methods
#
# References:
# - Foldnes, N., Moss, J., & Gronneberg, S. (2024). Improved goodness of fit
#   procedures for structural equation models. Structural Equation Modeling.
# - Du, H., & Bentler, P. M. (2022). 40-Year Old Unbiased Distribution Free
#   Estimator Reliably Improves SEM Statistics for Nonnormal Data.
# - Wu, H., & Lin, J. (2016). A Scaled F Distribution as an Approximation to
#   the Distribution of Test Statistics in Covariance Structure Analysis.

# =====================================================
# Parser (ported from semTests split_input)
# =====================================================

#' Check if test string is an FMG test
#' @param test Character string specifying the test
#' @return Logical indicating if this is an FMG test
#' @keywords internal
lav_test_fmg_is_fmg <- function(test) {
  test <- tolower(test)
  patterns <- c(
    "^fmg($|_)", "^peba", "^pols", "^pall", "^all$",
    "^all_", "^sb_", "^ss_", "^sf_", "^std_"
  )
  any(vapply(patterns, function(p) grepl(p, test), logical(1)))
}

lav_test_fmg_is_preset <- function(test) {
  grepl("^fmg($|_)", tolower(test))
}

lav_test_fmg_resolve_preset <- function(test, nested = FALSE,
                                        ngroups = 1L) {
  if (!lav_test_fmg_is_preset(test)) {
    return(test)
  }

  test <- tolower(test)
  splitted <- strsplit(test, "_")[[1L]]
  if (splitted[1L] != "fmg" || length(splitted) > 3L) {
    lav_msg_stop(gettextf(
      "FMG preset does not support test= %1$s.",
      dQuote(test)
    ))
  }

  unbiased <- isTRUE(nested)
  chisq <- if (isTRUE(nested) && ngroups > 1L) "ml" else "rls"

  if (length(splitted) > 1L) {
    suffix <- splitted[-1L]
    if (!all(suffix %in% c("ug", "ml", "rls")) ||
        sum(suffix %in% c("ml", "rls")) > 1L ||
        sum(suffix == "ug") > 1L) {
      lav_msg_stop(gettextf(
        "FMG preset does not support test= %1$s.",
        dQuote(test)
      ))
    }
    if ("ug" %in% suffix) {
      unbiased <- TRUE
    }
    if ("ml" %in% suffix) {
      chisq <- "ml"
    } else if ("rls" %in% suffix) {
      chisq <- "rls"
    }
  }

  method <- if (isTRUE(nested)) "pall" else "peba4"
  paste0(
    method,
    if (isTRUE(unbiased)) "_ug" else "",
    "_", chisq
  )
}

#' Parse FMG test string into components
#'
#' Format: {method}{param}_{ug?}_{chisq?}
#' Examples: peba4_ug_rls, eba2, pols2_ml, sb_ug
#'
#' @param string Character string specifying the test
#' @return List with method, param, unbiased, chisq
#' @keywords internal
lav_test_fmg_parse <- function(string) {
  string <- tolower(string)
  splitted <- strsplit(string, "_")[[1]]

  # Defaults
  method <- NULL
  param <- 2L # default for j (eba/peba) or gamma (pols)
  unbiased <- FALSE
  unbiased.explicit <- FALSE
  chisq <- "default"
  chisq.explicit <- FALSE

  type <- splitted[1]

  # Parse unbiased and chisq from suffix parts
  if (length(splitted) == 3L) {
    unbiased <- (splitted[2] == "ug")
    unbiased.explicit <- TRUE
    chisq <- splitted[3]
    chisq.explicit <- TRUE
  } else if (length(splitted) == 2L) {
    if (splitted[2] %in% c("rls", "ml")) {
      chisq <- splitted[2]
      chisq.explicit <- TRUE
    } else if (splitted[2] == "ug") {
      unbiased <- TRUE
      unbiased.explicit <- TRUE
    }
  }

  # Parse method and parameter
  if (startsWith(type, "peba")) {
    method <- "peba"
    param_str <- substring(type, 5)
    if (nchar(param_str) > 0L) param <- as.integer(param_str)
  } else if (startsWith(type, "pols")) {
    method <- "pols"
    param_str <- substring(type, 5)
    if (nchar(param_str) > 0L) param <- as.numeric(param_str)
  } else if (type == "pall") {
    method <- "pall"
  } else if (type == "all") {
    method <- "all"
  } else if (type %in% c("sb", "ss", "sf", "std")) {
    method <- type
  } else {
    lav_msg_stop(gettextf(
      "FMG tests do not support test= %1$s.",
      dQuote(type)
    ))
  }

  list(
    method = method,
    param = param,
    unbiased = unbiased,
    chisq = chisq,
    unbiased.explicit = unbiased.explicit,
    chisq.explicit = chisq.explicit
  )
}

lav_test_fmg_resolve_unbiased <- function(parsed, lavoptions = NULL) {
  if (isTRUE(parsed$unbiased.explicit)) {
    return(parsed$unbiased)
  }
  if (!is.null(lavoptions$gamma.unbiased)) {
    return(isTRUE(lavoptions$gamma.unbiased))
  }
  FALSE
}

lav_test_fmg_resolve_chisq <- function(parsed, lavoptions = NULL) {
  if (isTRUE(parsed$chisq.explicit)) {
    return(parsed$chisq)
  }

  scaled.test <- "standard"
  if (!is.null(lavoptions$scaled.test)) {
    scaled.test <- lavoptions$scaled.test[1L]
  }

  if (scaled.test %in% c("default", "standard", "none")) {
    return("ml")
  } else if (scaled.test == "browne.residual.nt.model") {
    return("rls")
  }

  lav_msg_stop(gettextf(
    "FMG tests require scaled.test= %1$s or %2$s; found %3$s.",
    dQuote("standard"), dQuote("browne.residual.nt.model"),
    dQuote(scaled.test)
  ))
}

lav_test_fmg_check_mi <- function(lavoptions = NULL,
                                       context = "FMG tests") {
  if (is.null(lavoptions) || is.null(lavoptions$missing)) {
    return(invisible(TRUE))
  }

  missing <- lavoptions$missing[1L]
  if (isTRUE(missing == "listwise")) {
    return(invisible(TRUE))
  }

  lav_msg_stop(gettextf(
    "%1$s require missing= %2$s; found %3$s.",
    context, dQuote("listwise"), dQuote(missing)
  ))
}

lav_test_fmg_check_estimator <- function(lavoptions = NULL,
                                         context = "FMG tests") {
  if (is.null(lavoptions)) {
    return(invisible(TRUE))
  }

  estimator <- lavoptions$estimator
  estimator.orig <- lavoptions$estimator.orig
  if (isTRUE(estimator == "ML") &&
      isTRUE(estimator.orig %in% c("ML", "MLM"))) {
    return(invisible(TRUE))
  }

  reported <- if (!is.null(estimator.orig)) estimator.orig else estimator
  lav_msg_stop(gettextf(
    "%1$s require estimator= %2$s or %3$s; found %4$s.",
    context, dQuote("ML"), dQuote("MLM"), dQuote(reported)
  ))
}

lav_test_fmg_ugamma_eigenvalues <- function(ugamma, df,
                                            truncate_negative = TRUE) {
  eigs <- Re(eigen(ugamma, only.values = TRUE)$values)
  lambdas_raw <- sort(eigs, decreasing = TRUE)[seq_len(df)]
  lambdas <- lambdas_raw
  truncated <- rep(FALSE, length(lambdas))

  if (isTRUE(truncate_negative)) {
    truncated <- lambdas < 0
    lambdas[truncated] <- 0
  }

  list(
    values = lambdas,
    raw = lambdas_raw,
    truncated = truncated,
    n_truncated = sum(truncated),
    min_raw = min(lambdas_raw)
  )
}

lav_test_fmg_label <- function(parsed, unbiased = FALSE, chisq = "ml") {
  param <- format(parsed$param, trim = TRUE, scientific = FALSE)
  method_label <- switch(parsed$method,
    "peba" = paste0("pEBA-", param),
    "pols" = paste0("pOLS-", param),
    "pall" = "PALL",
    "all" = "ALL",
    "sb" = "SB",
    "ss" = "scaled-shifted",
    "sf" = "scaled-F",
    "std" = "standard chi-square",
    parsed$method
  )
  chisq_label <- if (chisq == "rls") "RLS" else toupper(chisq)
  gamma_label <- if (isTRUE(unbiased)) "unbiased" else "biased"

  paste0(
    method_label,
    " p-value test (base: ", chisq_label,
    ", gamma: ", gamma_label, ")"
  )
}

lav_test_fmg_chisq_equivalent <- function(pvalue, df) {
  if (is.na(pvalue) || is.na(df) || df <= 0L) {
    return(as.numeric(NA))
  }
  pvalue <- max(0, min(1, pvalue))
  stats::qchisq(pvalue, df = df, lower.tail = FALSE)
}

lav_test_fmg_stat_group_equivalent <- function(stat_group = NULL,
                                               source_stat = NULL,
                                               stat = NULL) {
  if (is.null(stat_group)) {
    return(NULL)
  }
  if (length(stat_group) == 1L) {
    return(stat)
  }
  if (is.na(source_stat) || source_stat == 0 || is.infinite(stat)) {
    return(rep(stat / length(stat_group), length(stat_group)))
  }

  stat_group * stat / source_stat
}

lav_test_fmg_browne_nt_model <- function(lavobject = NULL,
                                         lavdata = NULL,
                                         lavsamplestats = NULL,
                                         lavmodel = NULL,
                                         lavpartable = NULL,
                                         lavoptions = NULL,
                                         lavh1 = NULL,
                                         lavimplied = NULL) {
  out <- try(
    lav_test_browne(
      lavobject = lavobject,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavmodel = lavmodel,
      lavpartable = lavpartable,
      lavoptions = lavoptions,
      lavh1 = lavh1,
      lavimplied = lavimplied,
      adf = FALSE,
      model_based = TRUE
    ),
    silent = TRUE
  )
  if (!inherits(out, "try-error")) {
    return(out)
  }

  if (!is.null(lavobject)) {
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    lavoptions <- lavobject@Options
    lavh1 <- lavobject@h1
    lavimplied <- lavobject@implied
  }

  if (!is.logical(lavoptions$gamma.n.minus.one)) {
    n.minus.one <- !(lavoptions$estimator == "ML" &&
      lavoptions$likelihood == "normal")
  } else {
    n.minus.one <- lavoptions$gamma.n.minus.one
  }

  delta <- lav_model_delta(lavmodel)
  gamma <- lav_object_gamma(
    lavobject = NULL,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavsamplestats = lavsamplestats,
    lavh1 = lavh1,
    lavimplied = lavimplied,
    adf = FALSE,
    model_based = TRUE
  )
  WLS.obs <- lavsamplestats@WLS.obs
  WLS.est <- lav_model_wls_est(lavmodel)
  nobs <- lavsamplestats@nobs
  ntotal <- lavsamplestats@ntotal
  ngroups <- length(WLS.obs)
  stat.group <- numeric(ngroups)

  lineq.flag <- lavmodel@eq.constraints || lavmodel@ceq.simple.only

  if (!lineq.flag) {
    for (g in seq_len(ngroups)) {
      RES <- WLS.obs[[g]] - WLS.est[[g]]
      Delta.c <- lav_mat_ortho_complement(delta[[g]])
      t_dgd <- crossprod(Delta.c, gamma[[g]]) %*% Delta.c
      t_dgd_inv <- lav_mat_sym_inverse(t_dgd)
      Ng <- if (n.minus.one) nobs[[g]] - 1L else nobs[[g]]
      t_res_delta_c <- crossprod(RES, Delta.c)
      stat.group[g] <-
        Ng * drop(t_res_delta_c %*% t_dgd_inv %*% t(t_res_delta_c))
    }
    STAT <- sum(stat.group)
  } else {
    RES.all <- do.call("c", WLS.obs) - do.call("c", WLS.est)
    Delta.all <- do.call("rbind", delta)
    if (lavmodel@eq.constraints) {
      Delta.g <- Delta.all %*% lavmodel@eq.constraints.K
    } else {
      Delta.g <- Delta.all %*% lavmodel@ceq.simple.K
    }
    Gamma.inv.weighted <- vector("list", ngroups)
    for (g in seq_len(ngroups)) {
      Ng <- if (n.minus.one) nobs[[g]] - 1L else nobs[[g]]
      Gamma.inv.temp <- try(solve(gamma[[g]]), silent = TRUE)
      if (inherits(Gamma.inv.temp, "try-error")) {
        Gamma.inv.temp <- MASS::ginv(gamma[[g]])
      }
      Gamma.inv.weighted[[g]] <- Gamma.inv.temp * Ng / ntotal
    }
    GI <- lav_mat_bdiag(Gamma.inv.weighted)
    t_dgid <- t(Delta.g) %*% GI %*% Delta.g
    t_dgid_inv <- MASS::ginv(t_dgid)
    q1 <- drop(t(RES.all) %*% GI %*% RES.all)
    q2 <- drop(t(RES.all) %*%
      GI %*% Delta.g %*% t_dgid_inv %*% t(Delta.g) %*% GI %*%
      RES.all)
    STAT <- ntotal * (q1 - q2)
    stat.group <- STAT * unlist(nobs) / ntotal
  }

  if (!is.null(lavobject)) {
    DF <- lavobject@test[[1]]$df
  } else {
    DF <- lav_pt_df(lavpartable)
    if (!lavmodel@cin.simple.only && nrow(lavmodel@con.jac) > 0L) {
      ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
      if (length(ceq.idx) > 0L) {
        DF <- DF + qr(lavmodel@con.jac[ceq.idx, , drop = FALSE])$rank
      }
    } else if (lavmodel@ceq.simple.only) {
      DF <- lav_pt_ndat(lavpartable) - max(lavpartable$free)
    }
  }

  list(
    test = "browne.residual.nt.model",
    stat = STAT,
    stat.group = stat.group,
    df = DF,
    refdistr = "chisq",
    pvalue = 1 - stats::pchisq(STAT, DF),
    label = "Browne's residual (NT model-based) test"
  )
}

# =====================================================
# Main entry point
# =====================================================

#' Compute FMG test statistics
#'
#' @param lavobject A lavaan object
#' @param lavsamplestats Sample statistics (extracted from lavobject if NULL)
#' @param lavmodel Model object (extracted from lavobject if NULL)
#' @param lavdata Data object (extracted from lavobject if NULL)
#' @param lavoptions Options (extracted from lavobject if NULL)
#' @param test_unscaled Unscaled test from lavobject@test[[1]]
#' @param test Character string specifying the test (e.g., "peba4_ug_rls")
#' @return List with test results including stat, df, pvalue
#' @keywords internal
lav_test_fmg <- function(lavobject = NULL,
                         lavsamplestats = NULL,
                         lavmodel = NULL,
                         lavimplied = NULL,
                         lavdata = NULL,
                         lavoptions = NULL,
                         test_unscaled = NULL,
                         test_chisq = NULL,
                         e_inv = NULL,
                         delta = NULL,
                         wls_v = NULL,
                         gamma = NULL,
                         test = "peba4") {
  # Extract from lavobject if provided
  if (!is.null(lavobject)) {
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavoptions <- lavobject@Options
    lavdata <- lavobject@Data
    if (is.null(test_unscaled)) {
      test_unscaled <- lavobject@test[[1]]
    }
  }

  lav_test_fmg_check_estimator(lavoptions, context = "FMG tests")
  lav_test_fmg_check_mi(lavoptions, context = "FMG tests")

  # Parse test string
  parsed <- lav_test_fmg_parse(test)
  chisq <- lav_test_fmg_resolve_chisq(parsed, lavoptions = lavoptions)
  unbiased <- lav_test_fmg_resolve_unbiased(parsed, lavoptions = lavoptions)
  label <- lav_test_fmg_label(parsed, unbiased = unbiased, chisq = chisq)

  if (is.null(test_chisq)) {
    test_chisq <- test_unscaled
    if (chisq == "rls") {
      if (is.null(lavobject)) {
        return(NULL)
      }
      test_chisq <- lav_test_fmg_browne_nt_model(lavobject = lavobject)
    }
  }
  chisq_stat <- test_chisq$stat

  # Get df
  df <- test_unscaled$df

  # Handle df == 0 case
  if (df == 0L || df < 0L) {
    return(list(
      test = test,
      stat = as.numeric(NA),
      source.stat = chisq_stat,
      source.stat.group = test_chisq$stat.group,
      df = df,
      pvalue = as.numeric(NA),
      refdistr = "chisq",
      method = parsed$method,
      param = parsed$param,
      unbiased = unbiased,
      chisq.type = chisq,
      label = label
    ))
  }

  # Get ugamma matrix (with unbiased option)
  ugamma <- lav_test_fmg_ugamma(
    lavobject = lavobject,
    lavsamplestats = lavsamplestats,
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    lavdata = lavdata,
    lavoptions = lavoptions,
    test_unscaled = test_unscaled,
    e_inv = e_inv,
    delta = delta,
    wls_v = wls_v,
    gamma = gamma,
    unbiased = unbiased
  )

  if (is.null(ugamma)) {
    return(list(
      test = test,
      stat = as.numeric(NA),
      source.stat = chisq_stat,
      source.stat.group = test_chisq$stat.group,
      df = df,
      pvalue = as.numeric(NA),
      refdistr = "chisq",
      method = parsed$method,
      param = parsed$param,
      unbiased = unbiased,
      chisq.type = chisq,
      label = label
    ))
  }

  # Compute eigenvalues (first df eigenvalues)
  eig <- lav_test_fmg_ugamma_eigenvalues(ugamma, df)
  lambdas <- eig$values

  # Compute p-value based on method
  pvalue <- switch(parsed$method,
    "peba" = lav_test_fmg_peba(chisq_stat, lambdas, j = parsed$param),
    "pols" = lav_test_fmg_pols(chisq_stat, lambdas, gamma = parsed$param),
    "pall" = lav_test_fmg_pall(chisq_stat, lambdas),
    "all" = lav_test_fmg_all(chisq_stat, lambdas),
    "sb" = lav_test_fmg_sb(chisq_stat, lambdas),
    "ss" = lav_test_fmg_ss(chisq_stat, lambdas, df),
    "sf" = lav_test_fmg_scaled_f(chisq_stat, lambdas),
    "std" = 1 - stats::pchisq(chisq_stat, df),
    lav_test_fmg_peba(chisq_stat, lambdas, j = 4L) # default
  )

  # Ensure p-value is in valid range
  pvalue <- max(0, min(1, pvalue))
  stat <- lav_test_fmg_chisq_equivalent(pvalue, df)
  stat.group <- lav_test_fmg_stat_group_equivalent(
    stat_group = test_chisq$stat.group,
    source_stat = chisq_stat,
    stat = stat
  )

  # Return TEST structure
  list(
    test = test,
    stat = stat,
    stat.group = stat.group,
    source.stat = chisq_stat,
    source.stat.group = test_chisq$stat.group,
    df = df,
    pvalue = pvalue,
    refdistr = "chisq",
    method = parsed$method,
    param = parsed$param,
    unbiased = unbiased,
    chisq.type = chisq,
    UGamma.eigenvalues = lambdas,
    UGamma.eigenvalues.raw = eig$raw,
    UGamma.eigenvalues.truncated = eig$n_truncated,
    label = label
  )
}

# =====================================================
# ugamma computation
# =====================================================

#' Get ugamma matrix with optional unbiased gamma
#'
#' @param lavobject A lavaan object
#' @param unbiased Logical: use Du & Bentler (2022) unbiased gamma?
#' @return ugamma matrix
#' @keywords internal
lav_test_fmg_ugamma <- function(lavobject = NULL,
                                lavsamplestats = NULL,
                                lavmodel = NULL,
                                lavimplied = NULL,
                                lavdata = NULL,
                                lavoptions = NULL,
                                test_unscaled = NULL,
                                e_inv = NULL,
                                delta = NULL,
                                wls_v = NULL,
                                gamma = NULL,
                                unbiased = FALSE) {
  if (!is.null(lavobject)) {
    lavsamplestats <- lavobject@SampleStats
    lavmodel <- lavobject@Model
    lavimplied <- lavobject@implied
    lavoptions <- lavobject@Options
    lavdata <- lavobject@Data
    if (is.null(test_unscaled)) {
      test_unscaled <- lavobject@test[[1]]
    }
  }

  ngroups <- lavdata@ngroups

  if (unbiased) {
    # Recompute gamma with unbiased = TRUE
    gamma <- vector("list", ngroups)
    for (g in seq_len(ngroups)) {
      gamma[[g]] <- lav_samp_gamma(
        m_y = lavdata@X[[g]],
        meanstructure = lavoptions$meanstructure,
        unbiased = TRUE
      )
    }
  } else {
    if (is.null(gamma)) {
      # stored NACOV, else recompute with the fit-time settings
      if (!is.null(lavobject)) {
        gamma <- lav_gamma_used(lavobject)
      } else {
        gamma <- lav_gamma_used(
          lavdata = lavdata,
          lavoptions = lavoptions,
          lavsamplestats = lavsamplestats,
          lavh1 = NULL,
          lavimplied = lavimplied
        )
      }
    }
    if (!is.list(gamma)) {
      gamma <- list(gamma)
    }
  }

  # Get U matrix and ugamma using existing infrastructure
  out <- lav_test_sb(
    lavobject = lavobject,
    lavsamplestats = lavsamplestats,
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    lavdata = lavdata,
    lavoptions = lavoptions,
    test_unscaled = test_unscaled,
    e_inv = e_inv,
    delta = delta,
    wls_v = wls_v,
    m_gamma = gamma,
    method = "original",
    return_ugamma = TRUE,
    return_u = TRUE
  )

  if (is.null(out)) {
    return(NULL)
  }

  # If unbiased, recompute ugamma with new gamma
  if (unbiased) {
    U <- out$UfromUGamma
    if (is.null(U)) {
      return(NULL)
    }

    # Rescale gamma by group weights
    fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal
    Gamma_scaled <- vector("list", length(gamma))
    for (g in seq_along(gamma)) {
      Gamma_scaled[[g]] <- gamma[[g]] / fg[g]
    }
    Gamma_all <- lav_mat_bdiag(Gamma_scaled)
    ugamma <- U %*% Gamma_all
  } else {
    ugamma <- out$UGamma # field name as returned by lav_test_sb()
  }

  ugamma
}

# =====================================================
# P-value methods
# =====================================================

#' Pure R Imhof p-value for quadratic forms
#'
#' Computes P[Q > q] where Q = sum(lambda_j * chi^2(h_j, delta_j)).
#' This is a small R port of the Imhof integration used by CompQuadForm,
#' using stats::integrate().
#'
#' @param q Observed quadratic-form value
#' @param lambda Eigenvalues/weights
#' @param h Degrees of freedom per component
#' @param delta Noncentrality parameters per component
#' @param epsabs Absolute integration tolerance
#' @param epsrel Relative integration tolerance
#' @param limit Maximum number of subdivisions
#' @return List with Qq and abserr
#' @keywords internal
lav_test_fmg_imhof <- function(q, lambda, h = rep(1, length(lambda)),
                               delta = rep(0, length(lambda)),
                               epsabs = 1e-6, epsrel = 1e-6,
                               limit = 10000L) {
  lambda <- as.numeric(lambda)
  h <- as.numeric(h)
  delta <- as.numeric(delta)

  keep <- abs(lambda) > .Machine$double.eps
  lambda <- lambda[keep]
  h <- h[keep]
  delta <- delta[keep]

  if (length(lambda) == 0L) {
    return(list(Qq = ifelse(q < 0, 1, 0), abserr = 0))
  }

  if (length(lambda) == 1L) {
    if (lambda > 0) {
      Qq <- stats::pchisq(q / lambda,
        df = h, ncp = delta,
        lower.tail = FALSE
      )
    } else if (q < 0) {
      Qq <- stats::pchisq(q / lambda, df = h, ncp = delta)
    } else {
      Qq <- 0
    }
    return(list(Qq = Qq, abserr = 0))
  }

  integrand <- function(u) {
    vapply(u, function(ui) {
      if (ui == 0) {
        return(0.5 * (sum((h + delta) * lambda) - q))
      }

      lu <- lambda * ui
      lu2 <- lu * lu

      theta <- 0.5 * sum(h * atan(lu) + delta * lu / (1 + lu2)) -
        0.5 * q * ui
      log.rho <- sum(0.25 * h * log1p(lu2) +
        0.5 * delta * lu2 / (1 + lu2))

      sin(theta) / (ui * exp(log.rho))
    }, numeric(1L))
  }

  out <- stats::integrate(
    integrand,
    lower = 0,
    upper = Inf,
    rel.tol = epsrel,
    abs.tol = epsabs,
    subdivisions = limit
  )

  list(Qq = 0.5 + out$value / pi, abserr = out$abs.error)
}

#' Penalized EBA p-value (Foldnes et al. 2024)
#'
#' @param chisq Chi-square statistic
#' @param lambdas Eigenvalues of ugamma
#' @param j Number of blocks
#' @return p-value
#' @keywords internal
lav_test_fmg_peba <- function(chisq, lambdas, j = 4L) {
  m <- length(lambdas)
  if (m == 0L) {
    return(as.numeric(NA))
  }

  k <- ceiling(m / j)
  eig <- c(lambdas, rep(NA_real_, k * j - m))
  dim(eig) <- c(k, j)
  eig_means <- colMeans(eig, na.rm = TRUE)
  eig_mean <- mean(lambdas)
  repeated <- rep(eig_means, each = k)[seq_len(m)]
  penalized <- (repeated + eig_mean) / 2

  lav_test_fmg_imhof(chisq, penalized)$Qq
}

#' Penalized OLS p-value (Foldnes et al. 2024)
#'
#' @param chisq Chi-square statistic
#' @param lambdas Eigenvalues of ugamma
#' @param gamma Penalization parameter
#' @return p-value
#' @keywords internal
lav_test_fmg_pols <- function(chisq, lambdas, gamma = 2) {
  m <- length(lambdas)
  if (m == 0L) {
    return(as.numeric(NA))
  }

  if (m == 1L) {
    lambda_hat <- lambdas
  } else {
    x <- seq_along(lambdas)
    beta1_hat <- 1 / gamma * stats::cov(x, lambdas) / stats::var(x)
    beta0_hat <- mean(lambdas) - beta1_hat * mean(x)
    lambda_hat <- pmax(beta0_hat + beta1_hat * x, 0)
  }

  lav_test_fmg_imhof(chisq, lambda_hat)$Qq
}

#' Penalized all eigenvalues (for nested models)
#'
#' @param chisq Chi-square statistic
#' @param lambdas Eigenvalues of ugamma
#' @return p-value
#' @keywords internal
lav_test_fmg_pall <- function(chisq, lambdas) {
  if (length(lambdas) == 0L) {
    return(as.numeric(NA))
  }

  penalized <- lambdas / 2 + mean(lambdas) / 2
  lav_test_fmg_imhof(chisq, penalized)$Qq
}

#' All eigenvalues exact p-value
#'
#' @param chisq Chi-square statistic
#' @param lambdas Eigenvalues of ugamma
#' @return p-value
#' @keywords internal
lav_test_fmg_all <- function(chisq, lambdas) {
  if (length(lambdas) == 0L) {
    return(as.numeric(NA))
  }

  lav_test_fmg_imhof(chisq, lambdas)$Qq
}

#' Satorra-Bentler via eigenvalues
#'
#' @param chisq Chi-square statistic
#' @param lambdas Eigenvalues of ugamma
#' @return p-value
#' @keywords internal
lav_test_fmg_sb <- function(chisq, lambdas) {
  m <- length(lambdas)
  if (m == 0L) {
    return(as.numeric(NA))
  }

  1 - stats::pchisq(chisq * m / sum(lambdas), df = m)
}

#' Scaled and shifted p-value
#'
#' @param chisq Chi-square statistic
#' @param lambdas Eigenvalues of ugamma
#' @param df Degrees of freedom
#' @return p-value
#' @keywords internal
lav_test_fmg_ss <- function(chisq, lambdas, df) {
  if (length(lambdas) == 0L || df == 0L) {
    return(as.numeric(NA))
  }

  tr_ug <- sum(lambdas)
  tr_ug2 <- sum(lambdas^2)

  if (tr_ug2 < .Machine$double.eps) {
    return(as.numeric(NA))
  }

  a <- sqrt(df / tr_ug2)
  b <- df - sqrt(df * tr_ug^2 / tr_ug2)
  t3 <- chisq * a + b

  1 - stats::pchisq(t3, df)
}

#' Scaled F p-value (Wu & Lin 2016)
#'
#' @param chisq Chi-square statistic
#' @param lambdas Eigenvalues of ugamma
#' @return p-value
#' @keywords internal
lav_test_fmg_scaled_f <- function(chisq, lambdas) {
  if (length(lambdas) == 0L) {
    return(as.numeric(NA))
  }

  s1 <- sum(lambdas)
  s2 <- sum(lambdas^2)
  s3 <- sum(lambdas^3)
  denom <- 2 * s1 * s2^2 - s1^2 * s3 + 2 * s2 * s3

  if (denom > 0) {
    d1f3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) / denom
    d2f3 <- (s1^2 * s2 + 2 * s2^2) / (s3 * s1 - s2^2) + 6
    if (d2f3 < 6) d2f3 <- Inf
    cf3 <- s1 * (s1^2 * s2 - 2 * s2^2 + 4 * s1 * s3) /
      (s1^2 * s2 - 4 * s2^2 + 6 * s1 * s3)
  } else {
    d1f3 <- Inf
    d2f3 <- s1^2 / s2 + 4
    cf3 <- s1 * (s1^2 + 2 * s2) / (s1^2 + 4 * s2)
  }

  1 - stats::pf(chisq / cf3, d1f3, d2f3)
}

# =====================================================
# Nested model implementation
# =====================================================

#' FMG test for nested model comparison
#'
#' @param m0 Restricted (null) model
#' @param m1 Unrestricted (alternative) model
#' @param test Test specification string (e.g., "pall_ug_ml")
#' @return List with test results
#' @keywords internal
lav_test_fmg_nested <- function(m0, m1, test = "pall") {
  # Ensure m0 has more df than m1
  if (m0@test[[1]]$df < m1@test[[1]]$df) {
    tmp <- m0
    m0 <- m1
    m1 <- tmp
  }

  lav_test_fmg_check_estimator(m0@Options, context = "FMG nested tests")
  lav_test_fmg_check_estimator(m1@Options, context = "FMG nested tests")
  if (!isTRUE(m0@Options$estimator.orig == m1@Options$estimator.orig)) {
    lav_msg_stop(gettextf(
      "FMG nested tests require both models to use the same estimator; found %1$s and %2$s.",
      dQuote(m0@Options$estimator.orig),
      dQuote(m1@Options$estimator.orig)
    ))
  }
  lav_test_fmg_check_mi(m0@Options, context = "FMG nested tests")
  lav_test_fmg_check_mi(m1@Options, context = "FMG nested tests")

  # Parse test string
  parsed <- lav_test_fmg_parse(test)
  chisq <- lav_test_fmg_resolve_chisq(parsed, lavoptions = m1@Options)
  unbiased <- lav_test_fmg_resolve_unbiased(parsed, lavoptions = m1@Options)
  label <- lav_test_fmg_label(parsed, unbiased = unbiased, chisq = chisq)
  if (!parsed$method %in% c("pall", "all", "peba", "pols")) {
    lav_msg_stop(gettextf(
      "FMG nested tests only support test= %1$s; found %2$s.",
      lav_msg_view(c("pall", "all", "peba", "pols"), "or"),
      dQuote(parsed$method)
    ))
  }

  # Get df difference
  df0 <- m0@test[[1]]$df
  df1 <- m1@test[[1]]$df
  df <- df0 - df1

  if (df == 0L) {
    lav_msg_stop(gettext(
      "Cannot test models with the same degrees of freedom."
    ))
  }

  # Get chi-square difference
  if (chisq == "rls") {
    chisq0 <- lav_test_fmg_browne_nt_model(lavobject = m0)$stat
    chisq1 <- lav_test_fmg_browne_nt_model(lavobject = m1)$stat
  } else {
    chisq0 <- m0@test[[1]]$stat
    chisq1 <- m1@test[[1]]$stat
  }
  chisq_diff <- chisq0 - chisq1

  # Get ugamma for nested comparison (Satorra 2000 projection)
  ugamma <- lav_test_fmg_ugamma_nested(m0, m1, unbiased = unbiased)

  if (is.null(ugamma)) {
    return(list(
      test = test,
      stat = as.numeric(NA),
      source.stat = chisq_diff,
      df = df,
      pvalue = as.numeric(NA),
      label = label
    ))
  }

  # Compute first df eigenvalues
  eig <- lav_test_fmg_ugamma_eigenvalues(ugamma, df)
  lambdas <- eig$values

  # Compute p-value (pall is recommended for nested)
  pvalue <- switch(parsed$method,
    "pall" = lav_test_fmg_pall(chisq_diff, lambdas),
    "all"  = lav_test_fmg_all(chisq_diff, lambdas),
    "peba" = lav_test_fmg_peba(chisq_diff, lambdas, j = parsed$param),
    "pols" = lav_test_fmg_pols(chisq_diff, lambdas, gamma = parsed$param)
  )
  pvalue <- max(0, min(1, pvalue))
  stat <- lav_test_fmg_chisq_equivalent(pvalue, df)

  list(
    test = test,
    stat = stat,
    source.stat = chisq_diff,
    df = df,
    pvalue = pvalue,
    chisq.type = chisq,
    unbiased = unbiased,
    UGamma.eigenvalues = lambdas,
    UGamma.eigenvalues.raw = eig$raw,
    UGamma.eigenvalues.truncated = eig$n_truncated,
    label = label
  )
}

#' Compute ugamma for nested models (Satorra 2000 projection)
#'
#' @param m0 Restricted model
#' @param m1 Unrestricted model
#' @param unbiased Use unbiased gamma?
#' @return ugamma matrix
#' @keywords internal
lav_test_fmg_ugamma_nested <- function(m0, m1, unbiased = FALSE) {
  lavdata <- m1@Data
  ngroups <- lavdata@ngroups

  # Get gamma from m1 (with optional unbiased)
  if (unbiased) {
    gamma <- vector("list", ngroups)
    for (g in seq_len(ngroups)) {
      gamma[[g]] <- lav_samp_gamma(
        m_y = lavdata@X[[g]],
        meanstructure = m1@Options$meanstructure,
        unbiased = TRUE
      )
    }
  } else {
    gamma <- lavTech(m1, "gamma")
    if (is.null(gamma)) {
      gamma <- lavTech(m0, "gamma")
    }
    if (is.null(gamma)) {
      lav_msg_stop(gettextf(
        paste(
          "Could not calculate the gamma matrix.",
          "Use estimator= %1$s or test= %2$s when fitting your lavaan model."
        ),
        dQuote("MLM"), dQuote("satorra.bentler")
      ))
    }
    if (!is.list(gamma)) gamma <- list(gamma)
  }

  wls_v <- lavTech(m1, "WLS.V")
  PI <- lav_model_delta(m1@Model)
  P.inv <- lavTech(m1, "inverted.information")

  if (is.null(P.inv)) {
    return(NULL)
  }

  A <- lav_test_diff_a(m1, m0, method = "delta", reference = "H1")

  # Handle equality constraints
  if (m1@Model@eq.constraints) {
    A <- A %*% t(m1@Model@eq.constraints.K)
  } else if (m1@Model@ceq.simple.only) {
    A <- A %*% t(m1@Model@ceq.simple.K)
  }

  # Safety check: remove zero rows/columns
  APA <- A %*% P.inv %*% t(A)
  col_sums <- colSums(APA)
  row_sums <- rowSums(APA)
  empty.idx <- which(abs(col_sums) < .Machine$double.eps^0.5 &
    abs(row_sums) < .Machine$double.eps^0.5)
  if (length(empty.idx) > 0L) {
    A <- A[-empty.idx, , drop = FALSE]
  }

  if (nrow(A) == 0L) {
    return(NULL)
  }

  # PAAPAAP projection
  PAAPAAP <- P.inv %*% t(A) %*% MASS::ginv(A %*% P.inv %*% t(A)) %*%
    A %*% P.inv

  # Build global matrices
  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal
  Gamma_f <- vector("list", length(gamma))
  for (g in seq_along(gamma)) {
    Gamma_f[[g]] <- gamma[[g]] / fg[g]
  }
  Gamma_all <- lav_mat_bdiag(Gamma_f)

  V_f <- wls_v
  for (g in seq_along(wls_v)) {
    V_f[[g]] <- fg[g] * wls_v[[g]]
  }
  V_all <- lav_mat_bdiag(V_f)

  PI_all <- do.call(rbind, PI)

  # U_global (eq. 22 in Satorra 2000)
  U_all <- V_all %*% PI_all %*% PAAPAAP %*% t(PI_all) %*% V_all
  U_all %*% Gamma_all
}
