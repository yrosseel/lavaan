# utility functions needed to compute various (robust) fit measures:
#
# - lav_fit_catml_dwls (for 'robust' RMSEA/CFI if data is categorical)
# - lav_fit_fiml_corrected (correct RMSEA/CFI if data is incomplete)

# compute scaling-factor (c.hat3) for fit.dwls, using fit.catml ingredients
# see:
#     Savalei, V. (2021) Improving Fit Indices In SEM with categorical data.
#     Multivariate Behavioral Research, 56(3), 390-407.
#

# YR Dec 2022: first version
# YR Jan 2023: catml_dwls should check if the input 'correlation' matrix
#              is positive-definite (or not)
# YR Jul 2026: replace check_pd = TRUE/FALSE by nonpd = "na"/"refit"/"smooth"
#              (what to do when an input correlation matrix is not
#              positive definite):
#                "na"     -> return NA for the robust quantities (default)
#                "refit"  -> smooth the input matrix, and re-estimate the
#                            model parameters using estimator = "catML"
#                            (the former check_pd = FALSE behavior)
#                "smooth" -> smooth the input matrix, but keep the (DWLS)
#                            parameter estimates (as in lavaan <= 0.6-13)
#              if the input matrices are positive definite, all three
#              choices are identical: the ML ingredients are evaluated at
#              the (unchanged) DWLS estimates

lav_fit_catml_dwls <- function(lavobject, nonpd = "na") {
  # empty list
  empty_list <- list(
    XX3 = as.numeric(NA), df3 = as.numeric(NA),
    c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # limitations
  if (!lavobject@Model@categorical ||
    lavobject@Data@nlevels > 1L ||
    lavobject@Options$conditional.x ||
    length(unlist(lavobject@pta$vnames$ov.num)) > 0L) {
    return(empty_list)
  } else {
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
  }

  # check if input matrix (or matrices) are all positive definite
  if (nonpd == "na") {
    for (g in seq_len(lavdata@ngroups)) {
      cor_1 <- lavsamplestats@cov[[g]]
      ev <- eigen(cor_1, symmetric = TRUE, only.values = TRUE)$values
      if (any(ev < .Machine$double.eps^(1 / 2))) {
        # non-pd! return NA
        # should we give a warning here? (not for now)
        # warning("lavaan WARNING: robust RMSEA/CFI could not be computed
        #   because the input correlation matrix is not positive-definite")
        return(empty_list)
      }
    }
  }

  # evaluate the ML ingredients using estimator = "catML": a plain
  # (re)evaluation at the DWLS estimates if the input matrices are
  # positive definite; if not, the input matrices are smoothed first, and
  # the model parameters are re-estimated (nonpd = "refit") or kept at
  # the DWLS values (nonpd = "smooth")
  fit_catml <- try(
    lav_object_catml(lavobject, allow_refit = (nonpd != "smooth")),
    silent = TRUE
  )
  if (inherits(fit_catml, "try-error")) {
    return(empty_list)
  }

  xx3 <- fit_catml@test[[1]]$stat
  df3 <- fit_catml@test[[1]]$df


  # compute 'k'
  v <- lavTech(fit_catml, "wls.v") # NT-ML weight matrix

  w_dwls <- lavTech(lavobject, "wls.v") # DWLS weight matrix
  gamma <- lavTech(lavobject, "gamma") # acov of polychorics
  delta <- lavTech(lavobject, "delta")
  e_inv <- lavTech(lavobject, "inverted.information")

  fg <- unlist(lavsamplestats@nobs) / lavsamplestats@ntotal

  # Fixme: as we only need the trace, perhaps we could do this
  # group-specific? (see lav_test_sb_trace_original)
  v_g <- v
  w_dwls_g <- w_dwls
  gamma_f <- gamma
  delta_g <- delta
  for (g in seq_len(lavdata@ngroups)) {
    ntotal <- nrow(gamma[[g]])
    nvar <- lavobject@Model@nvar[[g]]
    pstar <- nvar * (nvar - 1) / 2
    rm_idx <- seq_len(ntotal - pstar)

    # reduce
    delta_g[[g]] <- delta[[g]][-rm_idx, , drop = FALSE]
    # reduce and weight: Gamma_g / fg paired with fg-weighted V blocks (see
    # the SCALING CONVENTIONS note in lav_samplestats_gamma.R)
    w_dwls_g[[g]] <- fg[g] * w_dwls[[g]][-rm_idx, -rm_idx]
    v_g[[g]] <- fg[g] * v[[g]] # should already have the right dims
    gamma_f[[g]] <- 1 / fg[g] * gamma[[g]][-rm_idx, -rm_idx]
  }
  # create 'big' matrices
  w_dwls_all <- lav_mat_bdiag(w_dwls_g)
  v_all <- lav_mat_bdiag(v_g)
  gamma_all <- lav_mat_bdiag(gamma_f)
  delta_all <- do.call("rbind", delta_g)

  # compute trace
  wi_u_all <- diag(nrow(w_dwls_all)) -
              delta_all %*% e_inv %*% t(delta_all) %*% w_dwls_all
  ks <- sum(diag(t(wi_u_all) %*% v_all %*% wi_u_all %*% gamma_all))

  # convert to lavaan 'scaling.factor'
  c_hat3 <- ks / df3
  xx3_scaled <- xx3 / c_hat3

  # baseline model
  xx3_null <- fit_catml@baseline$test[[1]]$stat
  if (is.null(xx3_null)) {
    xx3_null <- as.numeric(NA)
    df3_null <- as.numeric(NA)
    kbs <- as.numeric(NA)
    c_hat3_null <- as.numeric(NA)
  } else {
    df3_null <- fit_catml@baseline$test[[1]]$df
    kbs <- sum(diag(gamma_all))
    c_hat3_null <- kbs / df3_null
  }

  # return values
  list(
    XX3 = xx3, df3 = df3, c.hat3 = c_hat3, XX3.scaled = xx3_scaled,
    XX3.null = xx3_null, df3.null = df3_null, c.hat3.null = c_hat3_null
  )
}


# compute ingredients to compute FIML-Corrected RMSEA/CFI
# see:
#     Zhang X, Savalei V. (2022). New computations for RMSEA and CFI
#     following FIML and TS estimation with missing data. Psychological Methods.
#
# h1_model: optional user-provided (less restrictive) h1 model fitted to the
# same data. The returned quantities are then those of the h0-vs-h1
# difference test: the complete-data statistics and df are differenced, and
# so are the correction traces k (as for the Satorra-Bentler difference
# test), giving c.hat3 = (k_h0 - k_h1) / (df3_h0 - df3_h1). The same holds
# for the baseline quantities (baseline-vs-h1).

lav_fit_fiml_corrected <- function(lavobject, baseline_model,
                                   version = "V3", h1_model = NULL) {
  version <- toupper(version)
  if (!version %in% c("V3", "V6")) {
    lav_msg_stop(gettext("only FIML-C(V3) and FIML-C(V6) are available."))
  }

  # empty list
  empty_list <- list(
    XX3 = as.numeric(NA), df3 = as.numeric(NA),
    c.hat3 = as.numeric(NA), XX3.scaled = as.numeric(NA),
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # limitations
  if (lavobject@Options$conditional.x ||
    lavobject@Data@nlevels > 1L ||
    is.null(lavobject@h1$implied$cov[[1]])) {
    return(empty_list)
  }

  # h0 model
  h0 <- lav_fit_fiml_k(lavobject, version = version)
  if (is.null(h0)) {
    return(empty_list)
  }

  # user-provided h1 model? (computed from its own saturated-model
  # information, so that its variable order does not need to match the
  # one of lavobject)
  h1 <- NULL
  if (!is.null(h1_model)) {
    stopifnot(inherits(h1_model, "lavaan"))
    if (h1_model@Options$conditional.x ||
      h1_model@Data@nlevels > 1L ||
      is.null(h1_model@h1$implied$cov[[1]])) {
      return(empty_list)
    }
    h1 <- lav_fit_fiml_k(h1_model, version = version)
    if (is.null(h1)) {
      return(empty_list)
    }
  }

  xx3 <- h0$xx3
  df3 <- h0$df3
  k_fimlc <- h0$k
  if (!is.null(h1)) {
    xx3 <- xx3 - h1$xx3
    df3 <- df3 - h1$df3
    k_fimlc <- k_fimlc - h1$k
  }

  # convert to lavaan 'scaling.factor'
  c_hat3 <- k_fimlc / df3
  xx3_scaled <- xx3 / c_hat3

  # collect temp results
  out <- list(
    XX3 = xx3, df3 = df3,
    c.hat3 = c_hat3, XX3.scaled = xx3_scaled,
    XX3.null = as.numeric(NA), df3.null = as.numeric(NA),
    c.hat3.null = as.numeric(NA)
  )

  # baseline model
  if (!is.null(baseline_model)) {
    fit_b <- baseline_model
  } else {
    fit_b <- try(lav_object_independence(lavobject), silent = TRUE)
  }

  if (inherits(fit_b, "try-error")) {
    return(out)
  }

  # the baseline model is fitted to the same data (same variable order)
  # as lavobject: reuse its saturated-model information pieces
  hb <- lav_fit_fiml_k(fit_b, version = version, shared = h0)
  if (is.null(hb)) {
    return(out)
  }

  xx3_null <- hb$xx3
  df3_null <- hb$df3
  kb_fimlc <- hb$k
  if (!is.null(h1)) {
    xx3_null <- xx3_null - h1$xx3
    df3_null <- df3_null - h1$df3
    kb_fimlc <- kb_fimlc - h1$k
  }

  # convert to lavaan 'scaling.factor'
  c_hat3_null <- kb_fimlc / df3_null

  # return values
  list(
    XX3 = xx3, df3 = df3, c.hat3 = c_hat3, XX3.scaled = xx3_scaled,
    XX3.null = xx3_null, df3.null = df3_null, c.hat3.null = c_hat3_null
  )
}

# FIML-C ingredients for a single fitted model: the complete-data ML test
# statistic (xx3) and df (df3) of the model refitted to the EM (saturated)
# sample statistics, and the correction trace k (Zhang & Savalei, 2022;
# version V3 or V6), so that the scaling factor is k / df3
#
# shared: optional output of a previous call for a model fitted to the same
# data in the same variable order; the saturated-model information pieces
# (which do not depend on the fitted model) are then reused
#
# returns NULL if any of the ingredients could not be computed
lav_fit_fiml_k <- function(fit, version = "V3", shared = NULL) {
  if (!is.null(shared)) {
    cov_tilde <- shared$cov_tilde
    mean_tilde <- shared$mean_tilde
    sample_nobs <- shared$sample_nobs
  } else {
    if (is.null(fit@h1$implied$cov[[1]])) {
      return(NULL)
    }
    h1 <- lavTech(fit, "h1", add.labels = TRUE)
    cov_tilde <- lapply(h1, "[[", "cov")
    mean_tilde <- lapply(h1, "[[", "mean")
    sample_nobs <- unlist(fit@SampleStats@nobs)
  }

  # 'refit' using 'tilde' (=EM/saturated) sample statistics
  # re-attach the data-based ov order (ov_order = "data"); parTable() has
  # stripped the "ovda" attribute, so without this fit_tilde would be built
  # in model order while fit (and its delta/information used below) is
  # in data order, yielding an order-dependent scaling factor. Harmless
  # no-op when the data order already equals the model order.
  pt_tilde <- parTable(fit)
  attr(pt_tilde, "ovda") <- fit@Data@ov.names[[1]]
  fit_tilde <- try(lavaan(
    model = pt_tilde,
    sample_cov = cov_tilde,
    sample_mean = mean_tilde,
    sample_nobs = sample_nobs,
    sample.cov.rescale = FALSE,
    information = "observed",
    optim.method = "none",
    se = "none",
    test = "standard",
    baseline = FALSE,
    fit.by.level = FALSE,
    check.post = FALSE
  ), silent = TRUE)
  if (inherits(fit_tilde, "try-error")) {
    return(NULL)
  }

  xx3 <- fit_tilde@test[[1]]$stat
  df3 <- fit_tilde@test[[1]]$df

  # V3/V6: always use h1.information = "unstructured"!!
  fit@Options$h1.information <- c("unstructured", "unstructured")
  fit@Options$observed.information <- c("h1", "h1")

  if (is.null(shared)) {
    fit_tilde@Options$h1.information <- c("unstructured", "unstructured")
    fit_tilde@Options$observed.information <- c("h1", "h1")

    # saturated-model information: missing-data (wm) and complete-data (wc)
    wm <- lav_model_h1_info_observed(fit)
    wc <- wc_g <- lav_model_h1_info_observed(fit_tilde)
    if (version == "V3") {
      jm <- jm_g <- lav_model_h1_info_firstorder(fit)
      gamma_f <- vector("list", length = fit@Data@ngroups)
    }
    wmi <- wmi_g <- try(lapply(wm, lav_mat_sym_inverse), silent = TRUE)
    if (inherits(wmi, "try-error")) {
      return(NULL)
    }

    fg <- unlist(fit@SampleStats@nobs) / fit@SampleStats@ntotal
    # Fixme: as we only need the trace, perhaps we could do this
    # group-specific? (see lav_test_sb_trace_original)
    for (g in seq_len(fit@Data@ngroups)) {
      # group weight
      wc_g[[g]] <- fg[g] * wc[[g]]
      wmi_g[[g]] <- 1 / fg[g] * wmi[[g]]

      # gamma
      if (version == "V3") {
        jm_g[[g]] <- fg[g] * jm[[g]]
        gamma_g <- wmi[[g]] %*% jm[[g]] %*% wmi[[g]]
        gamma_f[[g]] <- 1 / fg[g] * gamma_g
      }
    }
    # create 'big' matrices
    wc_all <- lav_mat_bdiag(wc_g)
    wmi_all <- lav_mat_bdiag(wmi_g)
    jm_all <- NULL
    tr11 <- tr1 <- as.numeric(NA)
    if (version == "V3") {
      gamma_all <- lav_mat_bdiag(gamma_f)
      # VS: Simplification of k.fimlc to minimize matrix multiplication
      #                                                 of big matrices
      jm_all <- lav_mat_bdiag(jm_g)
      # VS: tr11 is also used for baseline
      # VS: tr(AB) = sum(A*t(B)) is more efficient
      tr11 <- sum(wc_all * gamma_all)
    } else {
      # V6
      tr1 <- sum(wc_all * wmi_all)
    }
    shared <- list(
      cov_tilde = cov_tilde, mean_tilde = mean_tilde,
      sample_nobs = sample_nobs, wc_all = wc_all, wmi_all = wmi_all,
      jm_all = jm_all, tr11 = tr11, tr1 = tr1
    )
  }

  # model-specific pieces
  delta <- lavTech(fit, "delta")
  e_inv <- lavTech(fit, "inverted.information")
  delta_all <- do.call("rbind", delta)
  wc_all <- shared$wc_all
  wmi_all <- shared$wmi_all

  # compute trace
  if (version == "V3") {
    jm_all <- shared$jm_all
    tr12 <- sum((t(delta_all) %*% jm_all %*% wmi_all %*%
                 wc_all %*% delta_all) * e_inv)
    tr22 <- sum((t(delta_all) %*% jm_all %*% delta_all %*% e_inv)
             * t(t(delta_all) %*% wc_all %*% delta_all %*% e_inv))
    k <- shared$tr11 - 2 * tr12 + tr22
  } else {
    # V6
    e_comp <- t(delta_all) %*% wc_all %*% delta_all
    k <- shared$tr1 - sum(e_comp * e_inv)
  }

  # (drop the model-specific pieces of a reused 'shared' list first)
  shared[c("xx3", "df3", "k")] <- NULL
  c(shared, list(xx3 = xx3, df3 = df3, k = k))
}
